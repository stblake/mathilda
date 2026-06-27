/* plot3d.c — Plot3D[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...].
 *
 * Mirrors plot.c's shape as closely as the dimensionality allows: HoldAll
 * (the iterator vars have no value yet), a single split_options3() pass
 * that separates the sampler's own options from a passthrough list copied
 * onto the resulting Graphics3D[...] result -- exactly split_options's role
 * in plot.c, just with the smaller set of options that have a 3D meaning.
 * Anything not recognised here falls into the generic passthrough branch and
 * is inertly ignored by the renderer.
 *
 * The one place 3D genuinely cannot reuse 2D's sampler (sampling.c) is the
 * adaptive refinement itself: that sampler bisects an *ordered* 1D interval,
 * which has no 2D analogue without inventing a per-cell quadtree -- and a
 * quadtree creates T-junction cracks where differently-refined cells meet.
 * Instead MaxRecursion here doubles the *whole* grid's resolution when a
 * cheap flatness spot-check fails, capped at a few levels and a hard
 * point-count ceiling. This stays crack-free (every level is a uniform grid)
 * and gives MaxRecursion real meaning, at the cost of refining more than
 * strictly necessary -- an acceptable trade for "simple and clear" over a
 * true adaptive mesh. */

#include "plot3d.h"
#include "plot_common.h"
#include "sampling.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include <stdlib.h>
#include <math.h>

/* ---- PlotRange: explicit numeric z-band (mirrors plot.c's plotrange_yband) ---- */

static bool list2_nums3(Expr* e, double* a, double* b) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List
        && e->data.function.arg_count == 2
        && numericize_bound(e->data.function.args[0], a)
        && numericize_bound(e->data.function.args[1], b);
}

/* Recognizes a bare {zmin,zmax}, an {xspec,{zmin,zmax}} pair, or a
 * {{xmin,xmax},{ymin,ymax},{zmin,zmax}} triple -- in every case taking the
 * last element if it is itself a 2-number list, else the whole 2-number
 * list. Anything with a non-numeric (Automatic/All) z side leaves the
 * outputs untouched and returns false. */
static bool plotrange_zband(Expr* rhs, double* lo, double* hi) {
    Expr* v = evaluate(expr_copy(rhs));
    double a, b;
    bool ok = false;
    if (v->type == EXPR_FUNCTION && v->data.function.head
        && v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol == SYM_List) {
        size_t n = v->data.function.arg_count;
        if (n >= 2 && list2_nums3(v->data.function.args[n - 1], &a, &b)) ok = true;
        else if (n == 2 && numericize_bound(v->data.function.args[0], &a)
                 && numericize_bound(v->data.function.args[1], &b)) ok = true;
    }
    expr_free(v);
    if (ok && a > b) { double t = a; a = b; b = t; }
    if (ok && a < b) { *lo = a; *hi = b; return true; }
    return false;
}

/* ---- Option parsing ---- */

typedef struct {
    long plot_points;
    int  max_recursion;
    bool mesh;               /* Mesh -> All/True (default): overlay grid wireframe */
    Expr* region_function;   /* borrowed; held; NULL = none */
    Expr* color_function;    /* borrowed; held (function, or string "Rainbow"); NULL = none */
    bool  color_function_scaling; /* default true */
    Expr* exclusion_style;   /* borrowed; held; NULL = none (use default dark-gray line) */
    /* Explicit numeric z-band from PlotRange, fed to the refinement test the
     * same way Plot's yclip feeds sampling.c's flatness test. Degenerate
     * (lo >= hi) means "no explicit PlotRange z" -> use the full data extent. */
    double zclip_lo, zclip_hi;
} Plot3DSampleOpts;

/* Splits res's trailing Rule args (starting at index 3, after f/{x,..}/{y,..})
 * into the sampler options above and a passthrough list of evaluated Rule
 * expressions destined for the Graphics3D[...] result. */
static bool split_options3(Expr* res, Plot3DSampleOpts* sopts,
                            Expr*** passthrough_out, size_t* passthrough_count_out,
                            Expr** single_color_out) {
    sopts->plot_points = 25;
    sopts->max_recursion = 2;
    sopts->mesh = true;
    sopts->region_function = NULL;
    sopts->color_function = NULL;
    sopts->color_function_scaling = true;
    sopts->exclusion_style = NULL;
    sopts->zclip_lo = 0.0;
    sopts->zclip_hi = -1.0;
    *single_color_out = NULL;

    size_t argc = res->data.function.arg_count;
    /* +2 headroom for the Axes/PlotStyle defaults potentially appended
     * below when the caller didn't already supply them. */
    size_t cap = (argc > 3 ? argc - 3 : 0) + 2;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_style = false;

#define FAIL_CLEANUP3() do { free(passthrough); expr_free(*single_color_out); *single_color_out = NULL; return false; } while (0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) FAIL_CLEANUP3();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) FAIL_CLEANUP3();
            sopts->plot_points = v;
        } else if (name == SYM_MaxRecursion) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 0) FAIL_CLEANUP3();
            sopts->max_recursion = (int)v;
        } else if (name == SYM_Mesh) {
            /* Default True (unlike Plot's default None) -- Mathematica's
             * Plot3D shows grid lines on the surface out of the box. */
            Expr* v = evaluate(expr_copy(rhs));
            bool on = (v->type == EXPR_SYMBOL
                       && (v->data.symbol == SYM_All || v->data.symbol == SYM_True
                           || v->data.symbol == SYM_Automatic));
            bool off = (v->type == EXPR_SYMBOL
                        && (v->data.symbol == SYM_None || v->data.symbol == SYM_False));
            expr_free(v);
            if (!on && !off) FAIL_CLEANUP3();
            sopts->mesh = on;
        } else if (name == SYM_RegionFunction) {
            sopts->region_function = rhs;
        } else if (name == SYM_ColorFunction) {
            sopts->color_function = rhs;
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            sopts->color_function_scaling = !(v->type == EXPR_SYMBOL && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_ExclusionStyle) {
            /* Store the raw rhs for use in the sampler; also pass through to
             * the Graphics3D options so the renderer can inspect it if needed. */
            sopts->exclusion_style = rhs;
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else if (name == SYM_PlotStyle) {
            have_style = true;
            if (*single_color_out) expr_free(*single_color_out);
            *single_color_out = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), expr_copy(*single_color_out) };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        } else {
            if (name == SYM_Axes) have_axes = true;
            else if (name == SYM_PlotRange) {
                /* Pins the refinement's z-band; the value still passes
                 * through below for the renderer's box framing. */
                plotrange_zband(rhs, &sopts->zclip_lo, &sopts->zclip_hi);
            }
            /* Plot3D is HoldAll, so this option's value would otherwise
             * reach the Graphics3D[...] result unevaluated -- evaluate it
             * here once, exactly as split_options does for Plot. */
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_style) {
        Expr* rgb_args[3] = { expr_new_real(0.2), expr_new_real(0.4), expr_new_real(0.8) };
        Expr* rgb = expr_new_function(expr_new_symbol(SYM_RGBColor), rgb_args, 3);
        Expr* a[2] = { expr_new_symbol(SYM_PlotStyle), expr_copy(rgb) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        *single_color_out = rgb;
    }

    *passthrough_out = passthrough;
    *passthrough_count_out = n;
    return true;
}
#undef FAIL_CLEANUP3

/* ---- Grid sampling ---- */

typedef struct {
    double x, y, z; /* z = f(x,y) when fn_ok, else 0 */
    bool   valid;   /* fn_ok AND RegionFunction accepted */
    bool   fn_ok;   /* f(x,y) evaluated to a finite real (z is meaningful) */
} GridPt;

typedef struct {
    Expr* varx;            /* iterator symbol, borrowed */
    Expr* vary;             /* iterator symbol, borrowed */
    Expr* body;             /* f, borrowed */
    Expr* region_function;  /* borrowed; NULL = none */
} Plot3DEvalCtx;

/* RegionFunction: f[x,y,z] (3-arg) tried first (Mathematica's Plot3D
 * convention); anything that doesn't resolve to True/False falls back to
 * the shared 2-arg/1-arg idiom (eval_region) Plot already uses, so a
 * RegionFunction written for Plot keeps working under Plot3D too. */
static bool eval_region3(Expr* region_fn, double x, double y, double z) {
    Expr* args3[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    Expr* call3 = expr_new_function(expr_copy(region_fn), args3, 3);
    Expr* r3 = evaluate(call3);
    expr_free(call3);
    bool true3 = (r3->type == EXPR_SYMBOL && r3->data.symbol == SYM_True);
    bool false3 = (r3->type == EXPR_SYMBOL && r3->data.symbol == SYM_False);
    expr_free(r3);
    if (true3) return true;
    if (false3) return false;
    return eval_region(region_fn, x, y);
}

/* Evaluate f(x,y) and store the result in *z_out when finite. Does NOT
 * apply RegionFunction -- use plot3d_eval_fn for the combined check. */
static bool plot3d_eval_z(double x, double y, Plot3DEvalCtx* ctx, double* z_out) {
    Expr* xval = expr_new_real(x);
    Expr* yval = expr_new_real(y);
    symtab_add_own_value(ctx->varx->data.symbol, ctx->varx, xval);
    symtab_add_own_value(ctx->vary->data.symbol, ctx->vary, yval);
    Expr* result = evaluate(ctx->body);
    double z;
    bool ok = expr_to_real_double(result, &z) && isfinite(z);
    expr_free(result);
    expr_free(xval);
    expr_free(yval);
    if (ok) *z_out = z;
    return ok;
}

static bool plot3d_eval_fn(double x, double y, Plot3DEvalCtx* ctx, double* z_out) {
    Expr* xval = expr_new_real(x);
    Expr* yval = expr_new_real(y);
    symtab_add_own_value(ctx->varx->data.symbol, ctx->varx, xval);
    symtab_add_own_value(ctx->vary->data.symbol, ctx->vary, yval);
    Expr* result = evaluate(ctx->body);

    double z;
    bool ok = expr_to_real_double(result, &z) && isfinite(z);
    expr_free(result);

    if (ok && ctx->region_function && !eval_region3(ctx->region_function, x, y, z)) ok = false;

    expr_free(xval);
    expr_free(yval);
    if (ok) *z_out = z;
    return ok;
}

static GridPt* build_grid(Plot3DEvalCtx* ctx, double xmin, double xmax, double ymin, double ymax, long n) {
    GridPt* grid = malloc(sizeof(GridPt) * (size_t)n * (size_t)n);
    for (long i = 0; i < n; i++) {
        double x = xmin + (xmax - xmin) * (double)i / (double)(n - 1);
        if (i == n - 1) x = xmax; /* avoid float drift off the edge */
        for (long j = 0; j < n; j++) {
            double y = ymin + (ymax - ymin) * (double)j / (double)(n - 1);
            if (j == n - 1) y = ymax;
            double z = 0.0;
            bool fn_ok = plot3d_eval_z(x, y, ctx, &z);
            bool valid = fn_ok && (!ctx->region_function ||
                                    eval_region3(ctx->region_function, x, y, z));
            GridPt* p = &grid[i * n + j];
            p->x = x; p->y = y; p->z = fn_ok ? z : 0.0;
            p->valid = valid; p->fn_ok = fn_ok;
        }
    }
    return grid;
}

/* Clamp into the displayed z-band, exactly mirroring sampling.c's
 * clamp_band: a degenerate band (lo >= hi) means "no clip". */
static double clamp_band3(double z, double lo, double hi) {
    if (!(lo < hi)) return z;
    if (z < lo) return lo;
    if (z > hi) return hi;
    return z;
}

/* Maximum deviation of the surface from its bilinear interpolant at each
 * cell center, as a fraction of the displayed z-extent -- the 2D analogue
 * of sampling.c's chord-vs-curve flatness test (see the file header comment
 * for why this drives whole-grid doubling rather than a per-cell quadtree).
 * One extra f[x,y] evaluation per cell. */
#define FLAT_TOL3D 0.0025

static bool surface_is_flat(Plot3DEvalCtx* ctx, const GridPt* grid, long n,
                             double zclip_lo, double zclip_hi, double zspan) {
    for (long i = 0; i + 1 < n; i++) {
        for (long j = 0; j + 1 < n; j++) {
            const GridPt* p00 = &grid[i * n + j],     *p10 = &grid[(i + 1) * n + j];
            const GridPt* p01 = &grid[i * n + j + 1], *p11 = &grid[(i + 1) * n + j + 1];
            if (!p00->valid || !p10->valid || !p01->valid || !p11->valid) continue;

            double cx = (p00->x + p10->x) / 2.0;
            double cy = (p00->y + p01->y) / 2.0;
            double cz;
            if (!plot3d_eval_fn(cx, cy, ctx, &cz)) continue;

            double bilinear = (p00->z + p10->z + p01->z + p11->z) / 4.0;
            double dev = fabs(clamp_band3(cz, zclip_lo, zclip_hi) - clamp_band3(bilinear, zclip_lo, zclip_hi));
            if (dev > FLAT_TOL3D * zspan) return false;
        }
    }
    return true;
}

/* Samples f over [xmin,xmax]x[ymin,ymax] on a uniform grid, doubling the
 * whole grid's resolution (up to max_recursion times, hard-capped at 200
 * points/axis) while surface_is_flat() fails. Returns the final grid at
 * resolution *out_n x *out_n (caller frees with free()). */
static GridPt* sample_surface(Plot3DEvalCtx* ctx, double xmin, double xmax, double ymin, double ymax,
                               long plot_points, int max_recursion,
                               double zclip_lo, double zclip_hi, long* out_n) {
    long n = plot_points;
    GridPt* grid = build_grid(ctx, xmin, xmax, ymin, ymax, n);

    for (int level = 0; level < max_recursion && n < 200; level++) {
        size_t total = (size_t)n * (size_t)n;
        double* zs = malloc(sizeof(double) * total);
        size_t zc = 0;
        for (size_t k = 0; k < total; k++) {
            if (grid[k].valid) zs[zc++] = clamp_band3(grid[k].z, zclip_lo, zclip_hi);
        }
        double zlo, zhi;
        plot_robust_yrange(zs, zc, &zlo, &zhi);
        free(zs);
        double zspan = (zhi > zlo) ? zhi - zlo : 1.0;

        if (surface_is_flat(ctx, grid, n, zclip_lo, zclip_hi, zspan)) break;

        long nn = n * 2;
        if (nn > 200) nn = 200;
        free(grid);
        n = nn;
        grid = build_grid(ctx, xmin, xmax, ymin, ymax, n);
    }

    *out_n = n;
    return grid;
}

/* ---- Primitive assembly ---- */

static Expr* point3(double x, double y, double z) {
    Expr* a[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 3);
}

static Expr* mesh_line(const GridPt* a, const GridPt* b) {
    Expr* pts[2] = { point3(a->x, a->y, a->z), point3(b->x, b->y, b->z) };
    Expr* plist = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
    Expr* largs[1] = { plist };
    return expr_new_function(expr_new_symbol(SYM_Line), largs, 1);
}

/* ---- Sutherland-Hodgman boundary clipping ---- */

typedef struct {
    double x, y, z;
    bool   is_crossing; /* vertex inserted at the RegionFunction boundary */
} ClipPt;

/* Find t in (0,1) where the RegionFunction boundary crosses the segment
 * from V (inside, valid=true) to I (outside, fn_ok=true but valid=false).
 * Uses bisection: tlo stays inside, thi stays outside. */
static double find_crossing_t(Expr* region_fn,
                               const GridPt* V, const GridPt* I, int iters) {
    double tlo = 0.0, thi = 1.0;
    for (int k = 0; k < iters; k++) {
        double tm = (tlo + thi) * 0.5;
        double x  = V->x + tm * (I->x - V->x);
        double y  = V->y + tm * (I->y - V->y);
        double z  = V->z + tm * (I->z - V->z);
        if (eval_region3(region_fn, x, y, z)) tlo = tm; else thi = tm;
    }
    return (tlo + thi) * 0.5;
}

/* Sutherland-Hodgman clip of a quad against the RegionFunction boundary.
 * All four quad vertices must have fn_ok=true.  Vertices q0..q3 are given
 * in CCW order viewed from z+∞: q0=p00, q1=p10, q2=p11, q3=p01.
 * Writes up to 6 ClipPt vertices into `out`; returns the vertex count.
 * Crossing-point vertices have is_crossing=true; original in-region vertices
 * have is_crossing=false. */
static int clip_quad_to_region(Expr* region_fn,
                                const GridPt* q0, const GridPt* q1,
                                const GridPt* q2, const GridPt* q3,
                                ClipPt* out) {
    const GridPt* poly[4] = { q0, q1, q2, q3 };
    int nc = 0;

    for (int i = 0; i < 4; i++) {
        const GridPt* curr = poly[i];
        const GridPt* prev = poly[(i + 3) & 3];

        if (curr->valid) {
            if (!prev->valid) {
                /* Entering valid region: insert crossing on the edge prev→curr.
                 * find_crossing_t(V=curr, I=prev) returns t on curr→prev so
                 * the crossing is lerp(curr, prev, t). */
                double t = find_crossing_t(region_fn, curr, prev, 8);
                out[nc++] = (ClipPt){ curr->x + t*(prev->x - curr->x),
                                      curr->y + t*(prev->y - curr->y),
                                      curr->z + t*(prev->z - curr->z), true };
            }
            out[nc++] = (ClipPt){ curr->x, curr->y, curr->z, false };
        } else {
            /* Region-excluded (fn_ok guaranteed by caller): leaving the region. */
            if (prev->valid) {
                double t = find_crossing_t(region_fn, prev, curr, 8);
                out[nc++] = (ClipPt){ prev->x + t*(curr->x - prev->x),
                                      prev->y + t*(curr->y - prev->y),
                                      prev->z + t*(curr->z - prev->z), true };
            }
        }
    }
    return nc;
}

/* Compute the z-range [*zlo, *zhi] of all valid grid points using the same
 * robust range as the adaptive sampler's flatness test. Returns false when
 * the grid has no valid points (empty zrange, zlo > zhi). */
static bool grid_zrange(const GridPt* grid, long n, double* zlo, double* zhi) {
    size_t total = (size_t)n * (size_t)n;
    double* zs = malloc(sizeof(double) * total);
    size_t zc = 0;
    for (size_t k = 0; k < total; k++) {
        if (grid[k].valid) zs[zc++] = grid[k].z;
    }
    plot_robust_yrange(zs, zc, zlo, zhi);
    free(zs);
    return *zlo < *zhi;
}

/* Return the per-surface color directive for surface `fi` of `nfun`.
 * - single_color NULL -> fall back to palette
 * - single_color is a List -> index into it (cycling if short)
 * - single_color is a direct color -> use it for all surfaces */
static Expr* surface_color(Expr* single_color, size_t fi, size_t nfun) {
    (void)nfun;
    if (!single_color) return palette_color(fi);

    bool is_list = (single_color->type == EXPR_FUNCTION
                    && single_color->data.function.head->type == EXPR_SYMBOL
                    && single_color->data.function.head->data.symbol == SYM_List);
    if (is_list) {
        size_t len = single_color->data.function.arg_count;
        if (len == 0) return palette_color(fi);
        return expr_copy(single_color->data.function.args[fi % len]);
    }

    return expr_copy(single_color);
}

/* Sample every body over the grid and assemble the flat primitive
 * List[...] that becomes Graphics3D[...]'s first argument.
 *
 * Features implemented here:
 *   - Per-surface PlotStyle from single_color (List or direct color)
 *   - ColorFunction using full 3D (x, y, z) coordinates and z-range via
 *     eval_color_function3
 *   - Interior-only mesh lines: boundary perimeter edges are never drawn,
 *     giving a smooth silhouette instead of the jagged sawtooth that appears
 *     when all four cell edges are drawn unconditionally
 *   - ExclusionStyle: boundary edges between valid and excluded grid cells
 *     drawn in the exclusion_style directive when RegionFunction is set */
static Expr* build_surface_primitives(Expr** bodies, size_t nfun, Expr* varx, Expr* vary,
                                       double xmin, double xmax, double ymin, double ymax,
                                       const Plot3DSampleOpts* sopts, Expr* single_color) {
    if (nfun == 0 || !(xmin < xmax) || !(ymin < ymax)) return NULL;

    Rule* old_ownx = iter_spec_shadow(varx);
    Rule* old_owny = iter_spec_shadow(vary);

    Expr** prims = NULL;
    size_t cap = 0, prim_count = 0;
    bool any = false;
    bool multi = (nfun > 1);

#define PUSH(e) do { \
        if (prim_count == cap) { cap = cap ? cap * 2 : 256; prims = realloc(prims, sizeof(Expr*) * cap); } \
        prims[prim_count++] = (e); \
    } while (0)

    for (size_t fi = 0; fi < nfun; fi++) {
        Plot3DEvalCtx ctx = { .varx = varx, .vary = vary, .body = bodies[fi], .region_function = sopts->region_function };
        long n;
        GridPt* grid = sample_surface(&ctx, xmin, xmax, ymin, ymax, sopts->plot_points, sopts->max_recursion,
                                       sopts->zclip_lo, sopts->zclip_hi, &n);

        /* Per-surface color directive from PlotStyle (or palette fallback for
         * multi-surface plots). */
        if (multi) PUSH(surface_color(single_color, fi, nfun));

        /* Compute z-range once per surface for ColorFunction3 scaling. */
        double zlo = 0.0, zhi = 1.0;
        if (sopts->color_function) grid_zrange(grid, n, &zlo, &zhi);

        /* ---- Fill polygons ---- */
        for (long i = 0; i + 1 < n; i++) {
            for (long j = 0; j + 1 < n; j++) {
                GridPt* p00 = &grid[i * n + j],     *p10 = &grid[(i + 1) * n + j];
                GridPt* p01 = &grid[i * n + j + 1], *p11 = &grid[(i + 1) * n + j + 1];

                if (p00->valid && p10->valid && p01->valid && p11->valid) {
                    /* Fast path: fully inside the region. */
                    any = true;
                    if (sopts->color_function) {
                        double cx = (p00->x + p10->x + p01->x + p11->x) / 4.0;
                        double cy = (p00->y + p10->y + p01->y + p11->y) / 4.0;
                        double cz = (p00->z + p10->z + p01->z + p11->z) / 4.0;
                        PUSH(eval_color_function3(sopts->color_function,
                                                  cx, cy, cz,
                                                  xmin, xmax, ymin, ymax, zlo, zhi,
                                                  sopts->color_function_scaling));
                    }
                    Expr* verts[4] = { point3(p00->x, p00->y, p00->z),
                                       point3(p10->x, p10->y, p10->z),
                                       point3(p11->x, p11->y, p11->z),
                                       point3(p01->x, p01->y, p01->z) };
                    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
                    Expr* poly_args[1] = { vlist };
                    PUSH(expr_new_function(expr_new_symbol(SYM_Polygon), poly_args, 1));

                } else if (sopts->region_function
                           && p00->fn_ok && p10->fn_ok && p01->fn_ok && p11->fn_ok) {
                    /* Boundary cell: Sutherland-Hodgman clip against the RegionFunction
                     * boundary to produce a smooth partial polygon instead of a staircase. */
                    ClipPt clipped[8];
                    int nc = clip_quad_to_region(sopts->region_function,
                                                  p00, p10, p11, p01, clipped);
                    if (nc >= 3) {
                        any = true;
                        if (sopts->color_function) {
                            double cx = 0.0, cy = 0.0, cz = 0.0;
                            for (int k = 0; k < nc; k++) {
                                cx += clipped[k].x;
                                cy += clipped[k].y;
                                cz += clipped[k].z;
                            }
                            cx /= nc; cy /= nc; cz /= nc;
                            PUSH(eval_color_function3(sopts->color_function,
                                                      cx, cy, cz,
                                                      xmin, xmax, ymin, ymax, zlo, zhi,
                                                      sopts->color_function_scaling));
                        }
                        {
                            Expr** verts = malloc(sizeof(Expr*) * (size_t)nc);
                            for (int k = 0; k < nc; k++)
                                verts[k] = point3(clipped[k].x, clipped[k].y, clipped[k].z);
                            Expr* vlist = expr_new_function(expr_new_symbol(SYM_List),
                                                             verts, (size_t)nc);
                            free(verts);
                            Expr* poly_args[1] = { vlist };
                            PUSH(expr_new_function(expr_new_symbol(SYM_Polygon), poly_args, 1));
                        }
                    }
                }
                /* fn_ok=false cells: function undefined there — skip entirely. */
            }
        }

        /* ---- Interior mesh lines ----
         *
         * Draw the bottom and left edges of each valid cell, but only when
         * the adjacent cell on that side is also fully valid.  This gives
         * each interior grid line drawn exactly once and suppresses mesh
         * lines at both the rectangular domain perimeter (no cell beyond
         * the outer row/column) and at any RegionFunction exclusion
         * boundary (invalid neighbour → line omitted). */
        if (sopts->mesh) {
            bool mesh_started = false;
            for (long i = 0; i + 1 < n; i++) {
                for (long j = 0; j + 1 < n; j++) {
                    GridPt* p00 = &grid[i * n + j],     *p10 = &grid[(i + 1) * n + j];
                    GridPt* p01 = &grid[i * n + j + 1], *p11 = &grid[(i + 1) * n + j + 1];
                    if (!p00->valid || !p10->valid || !p01->valid || !p11->valid) continue;
                    /* Bottom edge: only when the cell below is also fully valid. */
                    if (j > 0) {
                        GridPt* b0 = &grid[i * n + (j - 1)], *b1 = &grid[(i + 1) * n + (j - 1)];
                        if (b0->valid && b1->valid) {
                            if (!mesh_started) {
                                Expr* _a[1] = { expr_new_real(0.15) };
                                PUSH(expr_new_function(expr_new_symbol(SYM_GrayLevel), _a, 1));
                                mesh_started = true;
                            }
                            PUSH(mesh_line(p00, p10));
                        }
                    }
                    /* Left edge: only when the cell to the left is also fully valid. */
                    if (i > 0) {
                        GridPt* l0 = &grid[(i - 1) * n + j], *l1 = &grid[(i - 1) * n + j + 1];
                        if (l0->valid && l1->valid) {
                            if (!mesh_started) {
                                Expr* _a[1] = { expr_new_real(0.15) };
                                PUSH(expr_new_function(expr_new_symbol(SYM_GrayLevel), _a, 1));
                                mesh_started = true;
                            }
                            PUSH(mesh_line(p01, p00));
                        }
                    }
                }
            }
        }

        /* ---- ExclusionStyle: boundary of the valid region ----
         *
         * Walk every fully-valid cell and draw each of its 4 edges that
         * borders either an invalid/excluded cell or the domain boundary.
         * Using cell-level validity (all 4 corners valid) is critical:
         * invalid grid points carry z=0 (the build_grid fallback), so the
         * old vertex-mismatch approach drew lines from the surface down to
         * the z=0 plane, making them invisible or misplaced.  Cell-level
         * edges always have both endpoints on the actual surface. */
        if (sopts->region_function) {
            bool exc_started = false;

#define ENSURE_EXC_STYLE() do { \
    if (!exc_started) { \
        if (sopts->exclusion_style) { \
            PUSH(expr_copy(sopts->exclusion_style)); \
        } else { \
            Expr* _a[1] = { expr_new_real(0.35) }; \
            PUSH(expr_new_function(expr_new_symbol(SYM_GrayLevel), _a, 1)); \
        } \
        exc_started = true; \
    } \
} while (0)

            /* Pass 1: staircase ExclusionStyle for fully-valid cells adjacent to
             * cells where fn_ok=false (function undefined in excluded region).
             * When fn_ok=true but region-excluded, the boundary is smooth — the
             * clipping pass below draws it precisely; skip it here to avoid jagged
             * staircase lines on top of the smooth outline. */
            for (long i = 0; i + 1 < n; i++) {
                for (long j = 0; j + 1 < n; j++) {
                    GridPt* p00 = &grid[i*n+j],     *p10 = &grid[(i+1)*n+j];
                    GridPt* p01 = &grid[i*n+j+1],   *p11 = &grid[(i+1)*n+j+1];
                    if (!p00->valid || !p10->valid || !p01->valid || !p11->valid) continue;

                    /* Bottom (p00–p10): draw only when the cell below is neither
                     * valid nor fn_ok (the fn_ok case is handled by the smooth pass). */
                    { bool nb = (j > 0)
                             && (grid[i*n+(j-1)].valid   || grid[i*n+(j-1)].fn_ok)
                             && (grid[(i+1)*n+(j-1)].valid || grid[(i+1)*n+(j-1)].fn_ok);
                      if (!nb) { ENSURE_EXC_STYLE(); PUSH(mesh_line(p00, p10)); } }
                    /* Top (p01–p11) */
                    { bool na = (j + 2 < n)
                             && (grid[i*n+(j+2)].valid   || grid[i*n+(j+2)].fn_ok)
                             && (grid[(i+1)*n+(j+2)].valid || grid[(i+1)*n+(j+2)].fn_ok);
                      if (!na) { ENSURE_EXC_STYLE(); PUSH(mesh_line(p01, p11)); } }
                    /* Left (p00–p01) */
                    { bool nl = (i > 0)
                             && (grid[(i-1)*n+j].valid   || grid[(i-1)*n+j].fn_ok)
                             && (grid[(i-1)*n+j+1].valid || grid[(i-1)*n+j+1].fn_ok);
                      if (!nl) { ENSURE_EXC_STYLE(); PUSH(mesh_line(p00, p01)); } }
                    /* Right (p10–p11) */
                    { bool nr = (i + 2 < n)
                             && (grid[(i+2)*n+j].valid   || grid[(i+2)*n+j].fn_ok)
                             && (grid[(i+2)*n+j+1].valid || grid[(i+2)*n+j+1].fn_ok);
                      if (!nr) { ENSURE_EXC_STYLE(); PUSH(mesh_line(p10, p11)); } }
                }
            }

            /* Pass 2: smooth ExclusionStyle from boundary-clipped cells.
             * For each cell where all fn_ok=true but some valid=false, the
             * Sutherland-Hodgman polygon has crossing-point vertices that lie
             * exactly on the RegionFunction boundary.  Edges between consecutive
             * crossing points trace the boundary smoothly (chord approximation
             * within one cell width). */
            for (long i = 0; i + 1 < n; i++) {
                for (long j = 0; j + 1 < n; j++) {
                    GridPt* p00 = &grid[i*n+j],     *p10 = &grid[(i+1)*n+j];
                    GridPt* p01 = &grid[i*n+j+1],   *p11 = &grid[(i+1)*n+j+1];
                    /* Only boundary cells (some invalid, all fn_ok). */
                    if (p00->valid && p10->valid && p01->valid && p11->valid) continue;
                    if (!p00->fn_ok || !p10->fn_ok || !p01->fn_ok || !p11->fn_ok) continue;

                    ClipPt clipped[8];
                    int nc = clip_quad_to_region(sopts->region_function,
                                                  p00, p10, p11, p01, clipped);
                    if (nc < 3) continue;

                    for (int k = 0; k < nc; k++) {
                        int k2 = (k + 1) % nc;
                        if (clipped[k].is_crossing && clipped[k2].is_crossing) {
                            ENSURE_EXC_STYLE();
                            GridPt ga, gb;
                            ga.x = clipped[k].x;  ga.y = clipped[k].y;
                            ga.z = clipped[k].z;  ga.valid = true; ga.fn_ok = true;
                            gb.x = clipped[k2].x; gb.y = clipped[k2].y;
                            gb.z = clipped[k2].z; gb.valid = true; gb.fn_ok = true;
                            PUSH(mesh_line(&ga, &gb));
                        }
                    }
                }
            }
#undef ENSURE_EXC_STYLE
        }

        free(grid);
    }
#undef PUSH

    iter_spec_restore(varx, old_ownx);
    iter_spec_restore(vary, old_owny);

    if (!any) { free(prims); return NULL; }
    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, prim_count);
    free(prims);
    return prim_list;
}

/* ---- Entry point ---- */

Expr* builtin_plot3d(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* xspec = res->data.function.args[1];
    Expr* yspec = res->data.function.args[2];

    IterSpec ix, iy;
    if (!iter_spec_parse(xspec, &ix)) return NULL;
    if (ix.kind != ITER_KIND_RANGE || !ix.var) { iter_spec_free(&ix); return NULL; }
    if (!iter_spec_parse(yspec, &iy)) { iter_spec_free(&ix); return NULL; }
    if (iy.kind != ITER_KIND_RANGE || !iy.var) { iter_spec_free(&ix); iter_spec_free(&iy); return NULL; }

    double xmin, xmax, ymin, ymax;
    if (!numericize_bound(ix.imin, &xmin) || !numericize_bound(ix.imax, &xmax) || !(xmin < xmax)
        || !numericize_bound(iy.imin, &ymin) || !numericize_bound(iy.imax, &ymax) || !(ymin < ymax)) {
        iter_spec_free(&ix);
        iter_spec_free(&iy);
        return NULL;
    }

    Plot3DSampleOpts sopts;
    Expr** passthrough = NULL;
    size_t passthrough_count = 0;
    Expr* single_color = NULL;
    if (!split_options3(res, &sopts, &passthrough, &passthrough_count, &single_color)) {
        iter_spec_free(&ix);
        iter_spec_free(&iy);
        return NULL;
    }

    /* One surface (Plot3D[f, …]) or several (Plot3D[{f1, f2, …}, …]). */
    Expr* single = f;
    Expr** bodies;
    size_t nfun;
    bool is_list = (f->type == EXPR_FUNCTION && f->data.function.head
                    && f->data.function.head->type == EXPR_SYMBOL
                    && f->data.function.head->data.symbol == SYM_List);
    if (is_list) {
        nfun = f->data.function.arg_count;
        bodies = f->data.function.args; /* borrowed */
    } else {
        nfun = 1;
        bodies = &single; /* borrowed */
    }

    if (nfun == 0) { /* Plot3D[{}, …] — nothing to draw */
        iter_spec_free(&ix);
        iter_spec_free(&iy);
        expr_free(single_color);
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    Expr* prim_list = build_surface_primitives(bodies, nfun, ix.var, iy.var, xmin, xmax, ymin, ymax,
                                                &sopts, single_color);
    iter_spec_free(&ix);
    iter_spec_free(&iy);
    expr_free(single_color);

    if (!prim_list) {
        for (size_t i = 0; i < passthrough_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        return NULL;
    }

    size_t gargc = 1 + passthrough_count;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < passthrough_count; i++) gargs[1 + i] = passthrough[i];
    free(passthrough);

    Expr* graphics3d = expr_new_function(expr_new_symbol(SYM_Graphics3D), gargs, gargc);
    free(gargs);

    /* The REPL front end renders any top-level Graphics3D[...] result, just
     * like Graphics[...] -- see repl.c. */
    return graphics3d;
}
