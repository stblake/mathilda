/* streamplot.c — StreamPlot[{vx,vy}, {x,xmin,xmax}, {y,ymin,ymax}, opts...]
 *
 * Traces streamlines of a 2-D vector field by RK4 integration from a grid
 * of seed points, emitting one Arrow[...] primitive per stream.  The result
 * is returned as a Graphics[...] object (auto-displayed by the REPL). */

#include "streamplot.h"
#include "plot_common.h"
#include "compile/autocompile.h"   /* machine fast path for the field components */
#include "show.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Forward declaration so we can use it inside option parsing          */
/* ------------------------------------------------------------------ */

static bool numericize(Expr* e, double* out) { return numericize_bound(e, out); }

/* ------------------------------------------------------------------ */
/* Vector field evaluation                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr* xvar;  /* borrowed */
    Expr* yvar;  /* borrowed */
    Expr* vx;    /* borrowed */
    Expr* vy;    /* borrowed */
    AutoCompiled* acx;   /* vx compiled as f(x,y); NULL = interpreter */
    AutoCompiled* acy;
} FieldCtx;

/* Compile both field components as f(x, y).  All-or-nothing: one interpreted
 * component costs the evaluator round-trip the fast path exists to avoid.
 * Streamline integration takes several field samples per step, so this is the
 * hottest sampler in the file. */
static void field_compile(FieldCtx* ctx) {
    const Expr* vars[2] = { ctx->xvar, ctx->yvar };
    ctx->acx = autocompile_new(ctx->vx, vars, 2);
    ctx->acy = ctx->acx ? autocompile_new(ctx->vy, vars, 2) : NULL;
    if (!ctx->acy) { autocompiled_free(ctx->acx); ctx->acx = NULL; }
}

static void field_free(FieldCtx* ctx) {
    autocompiled_free(ctx->acx); ctx->acx = NULL;
    autocompiled_free(ctx->acy); ctx->acy = NULL;
}

/* Evaluate the vector field at (x, y), storing results in (*vx_out, *vy_out).
 * Returns false if either component is non-finite or fails to evaluate. */
static bool eval_field(const FieldCtx* ctx, double x, double y,
                       double* vx_out, double* vy_out) {
    if (ctx->acx) {
        double in[2] = { x, y }, dvx, dvy;
        if (autocompiled_eval_real(ctx->acx, in, &dvx) && isfinite(dvx)
            && autocompiled_eval_real(ctx->acy, in, &dvy) && isfinite(dvy)) {
            *vx_out = dvx; *vy_out = dvy;
            return true;
        }
        /* Fall through: a component that declines may still have a real value
         * the interpreter can produce. */
    }
    Expr* xval = expr_new_real(x);
    Expr* yval = expr_new_real(y);
    symtab_add_own_value(ctx->xvar->data.symbol.name, ctx->xvar, xval);
    symtab_add_own_value(ctx->yvar->data.symbol.name, ctx->yvar, yval);

    Expr* rx = evaluate(ctx->vx);
    Expr* ry = evaluate(ctx->vy);

    double dvx, dvy;
    bool ok = expr_to_real_double(rx, &dvx) && isfinite(dvx)
           && expr_to_real_double(ry, &dvy) && isfinite(dvy);
    expr_free(rx);
    expr_free(ry);
    expr_free(xval);
    expr_free(yval);

    if (ok) { *vx_out = dvx; *vy_out = dvy; }
    return ok;
}

/* ------------------------------------------------------------------ */
/* Color function evaluation                                           */
/* ------------------------------------------------------------------ */

/* Evaluate StreamColorFunction at the stream's midpoint.
 * Tries: f[x,y,vx,vy,speed], f[x,y,vx,vy], f[x,y], f[speed],
 * and "Rainbow" (hue = scaled speed). Returns NULL if not set. */
static Expr* eval_stream_color(Expr* cfn,
                               double x, double y,
                               double vx, double vy,
                               double speed,
                               double speed_min, double speed_max) {
    if (!cfn) return NULL;

    if (cfn->type == EXPR_STRING) {
        double t = (speed_max > speed_min)
            ? (speed - speed_min) / (speed_max - speed_min) : 0.5;
        Expr* c = named_color_ramp(cfn->data.string, t);
        if (c) return c;
    }

    /* Try multi-arg forms, falling back to fewer args.
     * f[x,y,vx,vy,speed] → f[x,y,vx,vy] → f[x,y] → f[speed] */
    static const int arities[] = { 5, 4, 2, 1 };
    double vals[5] = { x, y, vx, vy, speed };
    for (size_t ai = 0; ai < 4; ai++) {
        int arity = arities[ai];
        Expr** fargs = malloc(sizeof(Expr*) * (size_t)arity);
        for (int j = 0; j < arity; j++) fargs[j] = expr_new_real(vals[j]);
        Expr* fcall = expr_new_function(expr_copy(cfn), fargs, (size_t)arity);
        free(fargs);
        Expr* result = evaluate(fcall);
        bool is_color = (result->type == EXPR_FUNCTION
                         && result->data.function.head
                         && result->data.function.head->type == EXPR_SYMBOL
                         && (result->data.function.head->data.symbol.name == SYM_RGBColor
                          || result->data.function.head->data.symbol.name == SYM_GrayLevel
                          || result->data.function.head->data.symbol.name == SYM_Hue));
        if (is_color) return result;
        expr_free(result);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Streamline geometry: normalized-field RK4 + evenly-spaced placement */
/* ------------------------------------------------------------------ */

typedef struct { double x, y; } Point2;

/* A growing point buffer for one half of a streamline. */
typedef struct { Point2* p; size_t n, cap; } PBuf;

static void pbuf_push(PBuf* b, double x, double y) {
    if (b->n == b->cap) {
        b->cap = b->cap ? b->cap * 2 : 64;
        b->p = realloc(b->p, sizeof(Point2) * b->cap);
    }
    b->p[b->n].x = x; b->p[b->n].y = y; b->n++;
}

/* Uniform spatial hash grid over the (data-space) domain, cell size = d_sep.
 * It holds every point of every completed streamline, so a candidate point can
 * be proximity-tested against all of them in O(1): any point within d_sep of a
 * query lies in the query's own cell or one of its eight neighbours. This is
 * what turns a hairball of overlapping streamlets into evenly-spaced
 * streamlines (Jobard & Lefebvre, "Creating Evenly-Spaced Streamlines of
 * Arbitrary Density", 1997). */
typedef struct {
    double   ox, oy;   /* domain min corner */
    double   cell;     /* = d_sep */
    int      nx, ny;   /* cell counts */
    Point2** bucket;   /* nx*ny point lists */
    size_t*  cnt;
    size_t*  cap;
} SGrid;

static void sgrid_init(SGrid* g, double xmin, double ymin,
                       double xmax, double ymax, double cell) {
    g->ox = xmin; g->oy = ymin; g->cell = cell;
    g->nx = (int)((xmax - xmin) / cell) + 1;
    g->ny = (int)((ymax - ymin) / cell) + 1;
    if (g->nx < 1) g->nx = 1;
    if (g->ny < 1) g->ny = 1;
    size_t ncell = (size_t)g->nx * (size_t)g->ny;
    g->bucket = calloc(ncell, sizeof(Point2*));
    g->cnt    = calloc(ncell, sizeof(size_t));
    g->cap    = calloc(ncell, sizeof(size_t));
}

static void sgrid_free(SGrid* g) {
    size_t ncell = (size_t)g->nx * (size_t)g->ny;
    for (size_t i = 0; i < ncell; i++) free(g->bucket[i]);
    free(g->bucket); free(g->cnt); free(g->cap);
}

static int sgrid_cellx(const SGrid* g, double x) {
    int c = (int)((x - g->ox) / g->cell);
    if (c < 0) c = 0;
    if (c >= g->nx) c = g->nx - 1;
    return c;
}
static int sgrid_celly(const SGrid* g, double y) {
    int c = (int)((y - g->oy) / g->cell);
    if (c < 0) c = 0;
    if (c >= g->ny) c = g->ny - 1;
    return c;
}

static void sgrid_add(SGrid* g, double x, double y) {
    size_t idx = (size_t)sgrid_celly(g, y) * (size_t)g->nx + (size_t)sgrid_cellx(g, x);
    if (g->cnt[idx] == g->cap[idx]) {
        g->cap[idx] = g->cap[idx] ? g->cap[idx] * 2 : 8;
        g->bucket[idx] = realloc(g->bucket[idx], sizeof(Point2) * g->cap[idx]);
    }
    g->bucket[idx][g->cnt[idx]].x = x;
    g->bucket[idx][g->cnt[idx]].y = y;
    g->cnt[idx]++;
}

/* True if any stored point lies within sqrt(r2) of (x,y). r2 must be <= cell^2
 * so the 3x3 neighbourhood is a sufficient search. */
static bool sgrid_within2(const SGrid* g, double x, double y, double r2) {
    int cx = sgrid_cellx(g, x), cy = sgrid_celly(g, y);
    for (int dy = -1; dy <= 1; dy++) {
        int yy = cy + dy;
        if (yy < 0 || yy >= g->ny) continue;
        for (int dx = -1; dx <= 1; dx++) {
            int xx = cx + dx;
            if (xx < 0 || xx >= g->nx) continue;
            size_t idx = (size_t)yy * (size_t)g->nx + (size_t)xx;
            const Point2* b = g->bucket[idx];
            size_t c = g->cnt[idx];
            for (size_t i = 0; i < c; i++) {
                double ex = b[i].x - x, ey = b[i].y - y;
                if (ex*ex + ey*ey < r2) return true;
            }
        }
    }
    return false;
}

/* Unit field direction at (x,y); false at eval failure or a near-zero field
 * (|v| <= crit -> a critical point, where the direction is undefined). */
static bool unit_field(const FieldCtx* ctx, double x, double y, double crit,
                       double* ux, double* uy) {
    double vx, vy;
    if (!eval_field(ctx, x, y, &vx, &vy)) return false;
    double s = sqrt(vx*vx + vy*vy);
    if (!(s > crit)) return false;
    *ux = vx / s; *uy = vy / s;
    return true;
}

/* One RK4 step of the NORMALIZED field, advancing a fixed arc-length h in
 * direction `sign` (+1 forward, -1 backward). Integrating the unit field (not
 * the raw field) is what makes the point spacing -- and therefore the rendered
 * curvature -- uniform along the whole streamline, independent of local speed;
 * the old raw-field step advanced by h*|v|, so a fast region drew a coarse,
 * jagged arc while a slow one piled up points. Returns false at a critical
 * point (the caller ends the line there). */
static bool rk4_step_unit(const FieldCtx* ctx, double x, double y,
                          double h, double sign, double crit,
                          double* nx, double* ny) {
    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y;
    double s = sign * h;
    if (!unit_field(ctx, x,           y,           crit, &k1x,&k1y)) return false;
    if (!unit_field(ctx, x+0.5*s*k1x, y+0.5*s*k1y, crit, &k2x,&k2y)) return false;
    if (!unit_field(ctx, x+0.5*s*k2x, y+0.5*s*k2y, crit, &k3x,&k3y)) return false;
    if (!unit_field(ctx, x+s*k3x,     y+s*k3y,     crit, &k4x,&k4y)) return false;
    *nx = x + s*(k1x + 2*k2x + 2*k3x + k4x)/6.0;
    *ny = y + s*(k1y + 2*k2y + 2*k3y + k4y)/6.0;
    return true;
}

/* Grow one direction of a streamline from (x0,y0), appending each new point
 * (the seed itself excluded) to `out`. Stops at: the domain boundary, a
 * RegionFunction rejection, a critical point, proximity to an existing
 * streamline (within sqrt(d_test2) of a grid point), a closed-orbit return to
 * the seed (after at least min_loop_arc of travel -- this is what lets a
 * rotational field draw whole closed circles instead of running to the step
 * cap), or max_steps. Sets *closed = true iff it stopped on the closed-orbit
 * return, so the caller can skip the redundant opposite pass. */
static void integrate_dir(const FieldCtx* ctx, const SGrid* grid,
                          double x0, double y0,
                          double xmin, double xmax, double ymin, double ymax,
                          double h, double sign, int max_steps, double crit,
                          double d_test2, double loop2, double min_loop_arc,
                          Expr* region_fn, PBuf* out, bool* closed) {
    double x = x0, y = y0, arc = 0.0;
    if (closed) *closed = false;
    for (int step = 0; step < max_steps; step++) {
        double nx, ny;
        if (!rk4_step_unit(ctx, x, y, h, sign, crit, &nx, &ny)) break;
        if (nx < xmin || nx > xmax || ny < ymin || ny > ymax) break;
        if (region_fn && !eval_region(region_fn, nx, ny)) break;
        if (grid && sgrid_within2(grid, nx, ny, d_test2)) break;
        x = nx; y = ny; arc += h;
        pbuf_push(out, x, y);
        if (arc > min_loop_arc) {
            double ex = x - x0, ey = y - y0;
            if (ex*ex + ey*ey < loop2) {        /* closed orbit */
                if (closed) *closed = true;
                break;
            }
        }
    }
}

/* Build one streamline through a seed. A CLOSED orbit is traced by a single
 * forward pass (the forward loop already returns to the seed, so integrating
 * backward too would just draw the same ring a second time, overlapping it —
 * the ring is closed exactly by repeating the seed at the end). An OPEN line is
 * grown BOTH ways so the seed sits in the middle of a line that extends as far
 * as it can each direction. Caller owns the returned array. */
static Point2* grow_streamline(const FieldCtx* ctx, const SGrid* grid,
                               double sx, double sy,
                               double xmin, double xmax, double ymin, double ymax,
                               double h, int max_steps, double crit,
                               double d_test2, double loop2, double min_loop_arc,
                               Expr* region_fn, size_t* n_out, bool* is_closed) {
    PBuf fwd = {NULL,0,0};
    bool closed = false;
    integrate_dir(ctx, grid, sx, sy, xmin,xmax,ymin,ymax, h, +1.0, max_steps,
                  crit, d_test2, loop2, min_loop_arc, region_fn, &fwd, &closed);
    if (is_closed) *is_closed = closed;

    if (closed) {
        size_t n = 1 + fwd.n + 1;               /* seed, ring, seed (closes it) */
        Point2* s = malloc(sizeof(Point2) * n);
        s[0].x = sx; s[0].y = sy;
        for (size_t i = 0; i < fwd.n; i++) s[1 + i] = fwd.p[i];
        s[n - 1].x = sx; s[n - 1].y = sy;
        free(fwd.p);
        *n_out = n;
        return s;
    }

    PBuf back = {NULL,0,0};
    integrate_dir(ctx, grid, sx, sy, xmin,xmax,ymin,ymax, h, -1.0, max_steps,
                  crit, d_test2, loop2, min_loop_arc, region_fn, &back, NULL);
    size_t n = back.n + 1 + fwd.n;
    Point2* s = malloc(sizeof(Point2) * n);
    size_t k = 0;
    for (size_t i = back.n; i > 0; i--) s[k++] = back.p[i-1];
    s[k].x = sx; s[k].y = sy; k++;
    for (size_t i = 0; i < fwd.n; i++) s[k++] = fwd.p[i];
    free(back.p); free(fwd.p);
    *n_out = n;
    return s;
}

/* Candidate seed, ordered by distance from the domain centre so placement
 * grows outward from the middle. */
typedef struct { double x, y, d2; } SeedCand;

static int seedcand_cmp(const void* a, const void* b) {
    double da = ((const SeedCand*)a)->d2, db = ((const SeedCand*)b)->d2;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* Option parsing                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int    stream_points;     /* controls line density (d_sep = ext/this); default 25 */
    double stream_scale;      /* max arc-length / domain diag; <=0 = run to natural end */
    Expr*  stream_style;      /* borrowed list of directives, or NULL */
    Expr*  color_function;    /* borrowed; NULL = none */
    Expr*  region_function;   /* borrowed; NULL = none */
    Expr*  legends;           /* borrowed; NULL = none */
    ScaleFnType sf_x, sf_y;   /* ScalingFunctions per-axis */
    bool   animate;           /* True = emit AnimatedStreamline instead of Line */
} StreamOpts;

static bool split_stream_options(Expr* res, StreamOpts* so,
                                  Expr*** passthrough_out, size_t* pt_count_out) {
    so->stream_points = 25;   /* line density: d_sep = min(extent)/25, chosen so
                               * a rotational field draws ~12 concentric rings
                               * that nearly touch — matching Mathematica's
                               * default StreamPlot density */
    so->stream_scale  = 0.0;  /* default: run each line to its natural end
                               * (another line / boundary / critical point /
                               * closed orbit), not a fixed short arc */
    so->stream_style  = NULL;
    so->color_function = NULL;
    so->region_function = NULL;
    so->legends = NULL;
    so->sf_x = SF_NONE;
    so->sf_y = SF_NONE;
    so->animate = false;

    size_t argc = res->data.function.arg_count;
    size_t cap = (argc > 3 ? argc - 3 : 0) + 4;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_aspect = false, have_frame = false;

#define SP_FAIL() do { free(passthrough); return false; } while(0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) SP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_StreamPoints || name == SYM_PlotPoints) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_INTEGER && v->data.integer >= 1) {
                so->stream_points = (int)v->data.integer;
                expr_free(v);
            } else {
                /* Automatic / Fine / Coarse keywords */
                expr_free(v);
                /* keep default */
            }
        } else if (name == SYM_StreamScale) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_None) {
                so->stream_scale = -1.0; /* run to natural end */
            } else if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) {
                so->stream_scale = 0.0;  /* run to natural end */
            } else {
                double sv;
                if (numericize(v, &sv) && sv > 0) so->stream_scale = sv;
            }
            expr_free(v);
        } else if (name == SYM_StreamStyle) {
            so->stream_style = rhs; /* borrowed */
        } else if (name == SYM_StreamColorFunction || name == SYM_ColorFunction) {
            so->color_function = rhs; /* borrowed */
        } else if (name == SYM_RegionFunction) {
            so->region_function = rhs; /* borrowed */
        } else if (name == SYM_PlotLegends) {
            so->legends = rhs; /* borrowed */
        } else if (name == SYM_ScalingFunctions) {
            Expr* v = evaluate(expr_copy(rhs));
            parse_scaling_functions(v, &so->sf_x, &so->sf_y);
            expr_free(v);
        } else if (name == SYM_StreamAnimate) {
            Expr* v = evaluate(expr_copy(rhs));
            so->animate = (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True);
            expr_free(v);
        } else {
            /* Pass through to Graphics[...] */
            if (name == SYM_Axes)         have_axes   = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol.name == SYM_False || rhs->data.symbol.name == SYM_None)))
                    have_frame = true;
            }
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* StreamPlot defaults: Frame->True, Axes->False, square aspect. */
    if (!have_axes && !have_frame) {
        Expr* fa[2] = { expr_new_symbol(SYM_Frame), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), fa, 2);
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    } else if (!have_axes) {
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *passthrough_out = passthrough;
    *pt_count_out = n;
    return true;

#undef SP_FAIL
}

/* ------------------------------------------------------------------ */
/* Build Line[...] and Arrow[...] Exprs from a Point2 array            */
/* ------------------------------------------------------------------ */

/* Build Arrow[{{x1,y1},...,{xn,yn}}] from the full point array. */
static Expr* make_arrow(const Point2* pts, size_t n) {
    if (n < 2) return NULL;
    Expr** coord_exprs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* xy[2] = { expr_new_real(pts[i].x), expr_new_real(pts[i].y) };
        coord_exprs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* pt_list = expr_new_function(expr_new_symbol(SYM_List), coord_exprs, n);
    free(coord_exprs);
    Expr* args[1] = { pt_list };
    return expr_new_function(expr_new_symbol(SYM_Arrow), args, 1);
}

/* Build Line[{pts[i0], ..., pts[i1]}] (inclusive) — a streamline segment. */
static Expr* make_line_range(const Point2* pts, size_t i0, size_t i1) {
    if (i1 <= i0) return NULL;
    size_t m = i1 - i0 + 1;
    Expr** coord_exprs = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) {
        Expr* xy[2] = { expr_new_real(pts[i0 + i].x), expr_new_real(pts[i0 + i].y) };
        coord_exprs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* pt_list = expr_new_function(expr_new_symbol(SYM_List), coord_exprs, m);
    free(coord_exprs);
    Expr* args[1] = { pt_list };
    return expr_new_function(expr_new_symbol(SYM_Line), args, 1);
}

/* Build AnimatedStreamline[{{x1,y1},...,{xn,yn}}] — same point list as
 * make_line(), but tagged with a distinct head so the renderer knows to
 * additionally walk it with moving particle dots (StreamAnimate -> True). */
static Expr* make_animated_streamline(const Point2* pts, size_t n) {
    if (n < 2) return NULL;
    Expr** coord_exprs = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* xy[2] = { expr_new_real(pts[i].x), expr_new_real(pts[i].y) };
        coord_exprs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* pt_list = expr_new_function(expr_new_symbol(SYM_List), coord_exprs, n);
    free(coord_exprs);
    Expr* args[1] = { pt_list };
    return expr_new_function(expr_new_symbol(SYM_AnimatedStreamline), args, 1);
}

/* Build a fixed-length Arrow[{p1, p2}] centred at pts[idx], oriented along the
 * local flow direction (pts[idx] - pts[idx-1]).  `len` is in plot coordinates
 * and is the same for every arrow, so a direction chevron in a tight vortex is
 * the same size as one in an open region.  Placing these at regular arc-length
 * spacing along a long streamline (rather than one lone arrow per short arc) is
 * what shows the flow direction without cluttering it. */
static Expr* arrow_at_index(const Point2* pts, size_t n, size_t idx, double len) {
    if (n < 2 || len <= 0 || idx == 0 || idx >= n) return NULL;
    double dx = pts[idx].x - pts[idx - 1].x;
    double dy = pts[idx].y - pts[idx - 1].y;
    double d  = sqrt(dx * dx + dy * dy);
    if (d < 1e-30) return NULL;
    double ux = dx / d, uy = dy / d;

    double cx = pts[idx].x, cy = pts[idx].y;
    double half = len * 0.5;
    Point2 seg[2] = {
        { cx - half * ux, cy - half * uy },
        { cx + half * ux, cy + half * uy },
    };
    return make_arrow(seg, 2);
}

/* Build a small filled arrowHEAD triangle whose tip sits at pts[idx] (the
 * downstream end of a dash), pointing along the local flow.  Emitting the head
 * as a bare Polygon — NOT an Arrow — is the whole point: Arrow[{a,b}] renders a
 * straight shaft PLUS a head, and that shaft, drawn as a chord across a curved
 * streamline, is the "arrow with a straight line" that clutters the plot.  A
 * lone triangle is just the chevron, so the curved dash flows straight into its
 * own arrowhead. `len` is the head length in plot coordinates — the same for
 * every head, independent of local field speed. */
static Expr* chevron_tip(const Point2* pts, size_t idx, double len) {
    if (idx == 0 || len <= 0) return NULL;
    /* Direction = chord from ~`len` of arc upstream to the tip, so the head is
     * aligned with the streamline over its own length rather than with one tiny
     * RK4 step (which would leave the triangle visibly askew on a tight bend). */
    size_t j = idx;
    double acc = 0.0;
    while (j > 0 && acc < len) {
        double sx = pts[j].x - pts[j - 1].x, sy = pts[j].y - pts[j - 1].y;
        acc += sqrt(sx * sx + sy * sy);
        j--;
    }
    double dx = pts[idx].x - pts[j].x;
    double dy = pts[idx].y - pts[j].y;
    double d  = sqrt(dx * dx + dy * dy);
    if (d < 1e-30) return NULL;
    double ux = dx / d, uy = dy / d;   /* flow direction */
    double nx = -uy, ny = ux;          /* unit normal      */
    double hw = len * 0.42;            /* half-width -> slender arrowhead */
    double tipx = pts[idx].x, tipy = pts[idx].y;
    double bx = tipx - len * ux, by = tipy - len * uy;
    Point2 tri[3] = {
        { bx + nx * hw, by + ny * hw },
        { bx - nx * hw, by - ny * hw },
        { tipx,         tipy         },
    };
    Expr** cs = malloc(sizeof(Expr*) * 3);
    for (int i = 0; i < 3; i++) {
        Expr* xy[2] = { expr_new_real(tri[i].x), expr_new_real(tri[i].y) };
        cs[i] = expr_new_function(expr_new_symbol(SYM_List), xy, 2);
    }
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), cs, 3);
    free(cs);
    Expr* args[1] = { lst };
    return expr_new_function(expr_new_symbol(SYM_Polygon), args, 1);
}

/* Append `e` to a growing primitive array (NULL is ignored). */
static void prim_append(Expr*** prims, size_t* n, size_t* cap, Expr* e) {
    if (!e) return;
    if (*n >= *cap) { *cap *= 2; *prims = realloc(*prims, sizeof(Expr*) * (*cap)); }
    (*prims)[(*n)++] = e;
}

/* Emit one streamline in Mathematica's dashed-arrow style: walk the arc length
 * and alternate `dash_len` drawn / `gap_len` skipped; each drawn dash is a
 * curve-following Line (so even a tight inner ring's dash stays curved) capped
 * by a filled arrowhead triangle at its downstream tip.  The gaps are what keep
 * neighbouring lines visually distinct at the high default density. */
static void emit_dashed_stream(Expr*** prims, size_t* nprim, size_t* cap,
                               const Point2* pts, size_t n,
                               double dash_len, double gap_len, double arrow_len) {
    if (n < 2 || !(dash_len > 0) || !(gap_len > 0)) {
        prim_append(prims, nprim, cap, make_line_range(pts, 0, n - 1));
        return;
    }
    double acc = 0.0;
    size_t dash_start = 0;
    bool in_dash = true;
    for (size_t k = 1; k < n; k++) {
        double ex = pts[k].x - pts[k - 1].x, ey = pts[k].y - pts[k - 1].y;
        acc += sqrt(ex * ex + ey * ey);
        if (in_dash && acc >= dash_len) {
            prim_append(prims, nprim, cap, make_line_range(pts, dash_start, k));
            prim_append(prims, nprim, cap, chevron_tip(pts, k, arrow_len));
            in_dash = false;
            acc = 0.0;
        } else if (!in_dash && acc >= gap_len) {
            in_dash = true;
            dash_start = k;
            acc = 0.0;
        }
    }
    /* trailing partial dash */
    if (in_dash && n - 1 > dash_start) {
        prim_append(prims, nprim, cap, make_line_range(pts, dash_start, n - 1));
        prim_append(prims, nprim, cap, chevron_tip(pts, n - 1, arrow_len));
    }
}

/* Default speed-to-color: the system default colormap (Viridis), keyed to the
 * stream's speed normalised across the plot. */
static Expr* default_stream_color(double spd, double spd_min, double spd_max) {
    double t = (spd_max > spd_min + 1e-30)
        ? (spd - spd_min) / (spd_max - spd_min) : 0.5;
    double rv, gv, bv;
    default_ramp_rgb(t, &rv, &gv, &bv);
    Expr* args[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), args, 3);
}

/* ------------------------------------------------------------------ */
/* Main builtin                                                         */
/* ------------------------------------------------------------------ */

Expr* builtin_streamplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;

    /* Arg 0: {vx, vy} — must be a 2-element List (held, not yet evaluated). */
    Expr* field_arg = res->data.function.args[0];
    if (!field_arg || field_arg->type != EXPR_FUNCTION
        || !field_arg->data.function.head
        || field_arg->data.function.head->type != EXPR_SYMBOL
        || field_arg->data.function.head->data.symbol.name != SYM_List
        || field_arg->data.function.arg_count != 2)
        return NULL;

    Expr* vx_body = field_arg->data.function.args[0]; /* borrowed */
    Expr* vy_body = field_arg->data.function.args[1]; /* borrowed */

    /* Arg 1: {x, xmin, xmax} */
    IterSpec xspec;
    if (!iter_spec_parse(res->data.function.args[1], &xspec)) return NULL;
    if (xspec.kind != ITER_KIND_RANGE || !xspec.var) { iter_spec_free(&xspec); return NULL; }
    double xmin, xmax;
    if (!numericize(xspec.imin, &xmin) || !numericize(xspec.imax, &xmax) || !(xmin < xmax)) {
        iter_spec_free(&xspec); return NULL;
    }

    /* Arg 2: {y, ymin, ymax} */
    IterSpec yspec;
    if (!iter_spec_parse(res->data.function.args[2], &yspec)) { iter_spec_free(&xspec); return NULL; }
    if (yspec.kind != ITER_KIND_RANGE || !yspec.var) { iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }
    double ymin, ymax;
    if (!numericize(yspec.imin, &ymin) || !numericize(yspec.imax, &ymax) || !(ymin < ymax)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL;
    }

    /* Options */
    StreamOpts so;
    Expr** passthrough = NULL;
    size_t pt_count = 0;
    if (!split_stream_options(res, &so, &passthrough, &pt_count)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL;
    }

    /* World-space (scaled) axis bounds — integration happens in data space,
     * but seed placement and PlotRange use world space. */
    double u_xmin = scale_apply(so.sf_x, xmin), u_xmax = scale_apply(so.sf_x, xmax);
    double u_ymin = scale_apply(so.sf_y, ymin), u_ymax = scale_apply(so.sf_y, ymax);

    /* Always embed an explicit PlotRange (world coords) so the renderer uses
     * the user's domain.  Skip if the user already passed their own PlotRange. */
    {
        bool have_pr = false;
        for (size_t i = 0; i < pt_count; i++) {
            const Expr* e = passthrough[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol.name == SYM_PlotRange)
                { have_pr = true; break; }
        }
        if (!have_pr) {
            passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
            Expr* xr[2] = { expr_new_real(u_xmin), expr_new_real(u_xmax) };
            Expr* yr[2] = { expr_new_real(u_ymin), expr_new_real(u_ymax) };
            Expr* xrange = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrange = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rlist[2] = { xrange, yrange };
            Expr* pr = expr_new_function(expr_new_symbol(SYM_List), rlist, 2);
            Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
            passthrough[pt_count++] =
                expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    }

    emit_scaling_meta(so.sf_x, so.sf_y, &passthrough, &pt_count);

    /* Shadow both iterator variables so the field sees no stray OwnValues. */
    Rule* old_x = iter_spec_shadow(xspec.var);
    Rule* old_y = iter_spec_shadow(yspec.var);

    FieldCtx ctx = { .xvar = xspec.var, .yvar = yspec.var,
                     .vx = vx_body, .vy = vy_body };
    field_compile(&ctx);

    /* ---- Evenly-spaced streamline placement (Jobard-Lefebvre) ----
     * All geometry (seeds, separation, integration, proximity) is in DATA
     * space; only the finished polylines are mapped to world/scaled space for
     * rendering. `d_sep` is the target spacing between neighbouring lines,
     * driven by StreamPoints. */
    double dx = xmax - xmin, dy = ymax - ymin;
    double ext = (dx < dy) ? dx : dy;
    double ddiag = sqrt(dx * dx + dy * dy);
    int    np = so.stream_points;
    if (np < 1) np = 1;

    double d_sep   = ext / np;                 /* target line spacing */
    if (!(d_sep > 0)) d_sep = ddiag > 0 ? ddiag * 0.05 : 1.0;
    double d_sep2  = d_sep * d_sep;
    double d_test2 = 0.25 * d_sep2;            /* stop ~0.5*d_sep from a line */
    double h       = d_sep / 8.0;              /* fixed arc-length step */
    if (!(h > 0)) h = 1e-4;
    double loop2   = (1.6 * h) * (1.6 * h);    /* closed-orbit return radius^2 */
    double min_loop_arc = 6.0 * h;

    /* World-space spacing for arrow sizing (streams are rendered in world
     * space; for the common identity scaling this equals d_sep). */
    double w_dx = u_xmax - u_xmin, w_dy = u_ymax - u_ymin;
    double w_ext = (w_dx < w_dy) ? w_dx : w_dy;
    double d_sep_world = w_ext / np;
    double arrow_len     = 0.85 * d_sep_world; /* chevron length (animate path) */
    double arrow_spacing = 3.2  * d_sep_world; /* arc between chevrons (animate path) */
    /* Mathematica-style dashed rendering (default, non-animate): each streamline
     * is drawn as dash / gap / dash ... with a small filled arrowhead triangle
     * at every dash tip.  dash_arrow is the arrowHEAD length (a bare chevron,
     * no shaft), kept short so it reads as a slender accent, not a chunky cap. */
    double dash_len   = 1.7   * d_sep_world;
    double gap_len    = 1.0   * d_sep_world;
    double dash_arrow = 0.336 * d_sep_world;  /* 0.8 × the original 0.42 head */

    /* Per-line arc-length cap: StreamScale -> s>0 limits each line to s*diag;
     * otherwise (default) a line runs to its natural end. */
    double max_arc = (so.stream_scale > 0) ? so.stream_scale * ddiag : 4.0 * ddiag;
    int max_steps = (int)(max_arc / h) + 1;
    if (max_steps > 20000) max_steps = 20000;
    if (max_steps < 1) max_steps = 1;

    /* Candidate seeds: a data-space grid ordered from the domain centre
     * outwards.  It is deliberately FINER than d_sep (2 candidates per
     * separation) so the d_sep proximity cull has enough positions to place
     * lines evenly and fill gaps left where a streamline terminated early —
     * candidate density controls coverage, d_sep controls the line spacing.
     * The scan also records the characteristic field speed used to set the
     * critical-point threshold. */
    double cxc = 0.5 * (xmin + xmax), cyc = 0.5 * (ymin + ymax);
    int cn = 2 * np;   /* candidates per axis */
    SeedCand* cand = malloc(sizeof(SeedCand) * (size_t)cn * (size_t)cn);
    size_t ncand = 0;
    double char_speed = 0.0;
    for (int iy = 0; iy < cn; iy++) {
        double sy = ymin + (iy + 0.5) * dy / cn;
        for (int ix = 0; ix < cn; ix++) {
            double sx = xmin + (ix + 0.5) * dx / cn;
            double tvx, tvy;
            if (!eval_field(&ctx, sx, sy, &tvx, &tvy)) continue;
            double sp = sqrt(tvx * tvx + tvy * tvy);
            if (sp > char_speed) char_speed = sp;
            double ex = sx - cxc, ey = sy - cyc;
            cand[ncand].x = sx; cand[ncand].y = sy; cand[ncand].d2 = ex*ex + ey*ey;
            ncand++;
        }
    }
    qsort(cand, ncand, sizeof(SeedCand), seedcand_cmp);
    double crit = (char_speed > 0) ? char_speed * 1e-4 : 1e-12;

    SGrid grid;
    sgrid_init(&grid, xmin, ymin, xmax, ymax, d_sep);

    /* At most one streamline per candidate seed. */
    size_t max_streams = ncand > 0 ? ncand : 1;
    Point2** all_streams = malloc(sizeof(Point2*) * max_streams);
    size_t*  all_lengths = malloc(sizeof(size_t)  * max_streams);
    double*  all_speed   = malloc(sizeof(double)  * max_streams);
    size_t nstreams = 0;

    for (size_t ci = 0; ci < ncand; ci++) {
        double sx = cand[ci].x, sy = cand[ci].y;

        /* Skip critical points and seeds too close to an existing line. */
        double tvx, tvy;
        if (!eval_field(&ctx, sx, sy, &tvx, &tvy)) continue;
        if (sqrt(tvx * tvx + tvy * tvy) <= crit) continue;
        if (so.region_function && !eval_region(so.region_function, sx, sy)) continue;
        if (sgrid_within2(&grid, sx, sy, d_sep2)) continue;

        size_t n_pts;
        bool is_closed = false;
        Point2* pts = grow_streamline(&ctx, &grid, sx, sy,
                                      xmin, xmax, ymin, ymax, h, max_steps, crit,
                                      d_test2, loop2, min_loop_arc,
                                      so.region_function, &n_pts, &is_closed);
        if (n_pts < 3) { free(pts); continue; }   /* discard stubs */

        /* Discard a sub-resolution closed loop: a seed that landed almost on a
         * critical point spins a tiny ring smaller than the line spacing (the
         * "hook" at the centre of a rotational field). Dropping it — without
         * registering it — leaves the next seed out at ~d_sep to draw the true
         * innermost ring, so the centre reads clean. Open lines are never
         * dropped this way (a short arc near the boundary is legitimate). */
        if (is_closed) {
            double xlo = pts[0].x, xhi = pts[0].x, ylo = pts[0].y, yhi = pts[0].y;
            for (size_t k = 1; k < n_pts; k++) {
                if (pts[k].x < xlo) xlo = pts[k].x;
                if (pts[k].x > xhi) xhi = pts[k].x;
                if (pts[k].y < ylo) ylo = pts[k].y;
                if (pts[k].y > yhi) yhi = pts[k].y;
            }
            double bdx = xhi - xlo, bdy = yhi - ylo;
            if (sqrt(bdx*bdx + bdy*bdy) < d_sep) { free(pts); continue; }
        }

        /* Register every point so later lines keep their distance. */
        for (size_t k = 0; k < n_pts; k++) sgrid_add(&grid, pts[k].x, pts[k].y);

        /* Midpoint speed (data space) for colourization. */
        size_t mid = n_pts / 2;
        double svx, svy, spd = 0.0;
        if (eval_field(&ctx, pts[mid].x, pts[mid].y, &svx, &svy))
            spd = sqrt(svx * svx + svy * svy);

        /* Map to world/scaled space for rendering. */
        if (so.sf_x != SF_NONE || so.sf_y != SF_NONE) {
            for (size_t k = 0; k < n_pts; k++) {
                pts[k].x = scale_apply(so.sf_x, pts[k].x);
                pts[k].y = scale_apply(so.sf_y, pts[k].y);
            }
        }

        all_streams[nstreams] = pts;
        all_lengths[nstreams] = n_pts;
        all_speed[nstreams]   = spd;
        nstreams++;
    }

    free(cand);
    sgrid_free(&grid);
    field_free(&ctx);

    /* Restore iterator variable bindings. */
    iter_spec_restore(xspec.var, old_x);
    iter_spec_restore(yspec.var, old_y);

    if (nstreams == 0) {
        free(all_streams); free(all_lengths); free(all_speed);
        iter_spec_free(&xspec); iter_spec_free(&yspec);
        for (size_t i = 0; i < pt_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    /* Speed range for ColorFunction scaling. */
    double spd_min = all_speed[0], spd_max = all_speed[0];
    for (size_t i = 1; i < nstreams; i++) {
        if (all_speed[i] < spd_min) spd_min = all_speed[i];
        if (all_speed[i] > spd_max) spd_max = all_speed[i];
    }

    /* ---- Build primitive list ---- */
    /* Per-stream layout: [color?, Line[all_pts], Arrow, Arrow, ...]
     *
     * color source (in priority order):
     *   1. StreamColorFunction / ColorFunction  — explicit user function
     *   2. speed-based default (blue-violet → orange) — when no StreamStyle
     *   3. StreamStyle directive              — single global style, no per-stream color
     *
     * Each Arrow is a short 2-point chevron placed at regular arc-length
     * spacing along the shaft (see below), a small direction indicator rather
     * than a large cap. */
    bool have_style = (so.stream_style != NULL);
    bool have_cfn   = (so.color_function != NULL);
    bool use_speed_color = !have_style && !have_cfn;

    /* Dashed rendering emits ~2 primitives (Line + chevron) per dash; prim_append
     * grows this as needed, so the initial size is just a reasonable guess. */
    size_t prim_cap = nstreams * 24 + 8;
    Expr** prims = malloc(sizeof(Expr*) * prim_cap);
    size_t nprim = 0;

    /* Global thickness directive (light line weight, closer to Mathematica). */
    {
        Expr* thick_arg[1] = { expr_new_real(0.004375) };  /* 1.25 × the original 0.0035 */
        prims[nprim++] = expr_new_function(expr_new_symbol(SYM_Thickness), thick_arg, 1);
    }

    /* If StreamStyle given, apply it once before all streams. */
    if (have_style) {
        prims[nprim++] = evaluate(expr_copy(so.stream_style));
    }

    for (size_t i = 0; i < nstreams; i++) {
        Point2* pts   = all_streams[i];
        size_t  n_pts = all_lengths[i];
        double  spd   = all_speed[i];

        /* Per-stream color. */
        Expr* col = NULL;
        if (have_cfn) {
            size_t mid = n_pts / 2;
            double svx = 0, svy = 0;
            Rule* sx2 = iter_spec_shadow(xspec.var);
            Rule* sy2 = iter_spec_shadow(yspec.var);
            eval_field(&ctx, pts[mid].x, pts[mid].y, &svx, &svy);
            iter_spec_restore(xspec.var, sx2);
            iter_spec_restore(yspec.var, sy2);
            col = eval_stream_color(so.color_function,
                                    pts[mid].x, pts[mid].y,
                                    svx, svy, spd, spd_min, spd_max);
        } else if (use_speed_color) {
            col = default_stream_color(spd, spd_min, spd_max);
        }

        prim_append(&prims, &nprim, &prim_cap, col);

        if (so.animate) {
            /* Animate: one solid AnimatedStreamline[...] (the renderer walks it
             * with moving particle dots) plus periodic direction chevrons. */
            prim_append(&prims, &nprim, &prim_cap, make_animated_streamline(pts, n_pts));
            if (arrow_len > 0) {
                double acc = arrow_spacing * 0.5;
                size_t placed = 0;
                for (size_t k = 1; k < n_pts; k++) {
                    double ex = pts[k].x - pts[k-1].x, ey = pts[k].y - pts[k-1].y;
                    acc += sqrt(ex * ex + ey * ey);
                    if (acc >= arrow_spacing) {
                        acc = 0.0;
                        prim_append(&prims, &nprim, &prim_cap,
                                    arrow_at_index(pts, n_pts, k, arrow_len));
                        placed++;
                    }
                }
                if (placed == 0)
                    prim_append(&prims, &nprim, &prim_cap,
                                arrow_at_index(pts, n_pts, n_pts / 2 ? n_pts / 2 : 1, arrow_len));
            }
        } else {
            /* Default: Mathematica-style dashed arrow-segments. */
            emit_dashed_stream(&prims, &nprim, &prim_cap, pts, n_pts,
                               dash_len, gap_len, dash_arrow);
        }
        free(pts);
    }
    free(all_streams); free(all_lengths); free(all_speed);

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, nprim);
    free(prims);

    /* Color-scale bar metadata: $StreamColorBar[spd_min, spd_max].
     * Only emitted when PlotLegends is set and speed-based coloring is active.
     * Replaces the discrete swatch legend used by Plot — StreamPlot's legend
     * is a continuous gradient bar, not a labeled-line list. */
    Expr* color_bar_meta = NULL;
    if (so.legends && use_speed_color) {
        Expr* eval_legends = evaluate(expr_copy(so.legends));
        bool want_legend = !(eval_legends->type == EXPR_SYMBOL
                             && (eval_legends->data.symbol.name == SYM_None
                              || eval_legends->data.symbol.name == SYM_False));
        expr_free(eval_legends);
        if (want_legend) {
            Expr* cfn_copy = so.color_function
                             ? expr_copy(so.color_function)
                             : expr_new_symbol(SYM_Automatic);
            Expr* cb_args[3] = { expr_new_real(spd_min), expr_new_real(spd_max), cfn_copy };
            color_bar_meta = expr_new_function(
                expr_new_symbol(SYM_StreamColorBar), cb_args, 3);
        }
    }

    /* Assemble Graphics[primitives, opts..., colorBarMeta?]. */
    size_t gargc = 1 + pt_count + (color_bar_meta ? 1 : 0);
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < pt_count; i++) gargs[1 + i] = passthrough[i];
    if (color_bar_meta) gargs[gargc - 1] = color_bar_meta;
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    iter_spec_free(&xspec);
    iter_spec_free(&yspec);

    return graphics;
}
