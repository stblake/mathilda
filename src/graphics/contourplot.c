/* contourplot.c — ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
 *
 * Generates iso-contour lines of a 2D function f(x, y) using the marching
 * squares algorithm and returns a Graphics[...] object auto-displayed by the
 * REPL.  ContourPlot is HoldAll: f and the iterator specs are held unevaluated
 * until x and y are bound to numeric values.
 *
 * Algorithm:
 *   1. Evaluate f on a (PlotPoints+1) × (PlotPoints+1) grid.
 *   2. Choose contour levels: explicit list, or N evenly spaced levels.
 *   3. Optionally shade each grid cell by its z value (ContourShading).
 *   4. For each level: run marching squares over the grid, emitting Line[]
 *      segments.  Saddle cells (states 5 and 10) use the bilinear centre
 *      value to pick the correct one of the two possible pairings.
 *   5. Optionally label each level at the midpoint of its first segment.
 *   6. Wrap everything in Graphics[...] with Plot's default options. */

#include "contourplot.h"
#include "plot_common.h"
#include "show.h"
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "print.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Build an RGBColor[] for a normalised z value t ∈ [0,1]. When a custom
 * ColorFunction is given it is tried 1-arg (t) then 2-arg (x,y); falls
 * back to the thermal ramp.  Caller owns the returned Expr. */
static Expr* contour_color(Expr* cfn, double t, double x, double y) {
    if (cfn) {
        if (cfn->type == EXPR_STRING) {
            const char* nm = cfn->data.string;
            if (strcmp(nm, "Rainbow") == 0) {
                Expr* h_arg[1] = { expr_new_real(t * 0.8) };
                return expr_new_function(expr_new_symbol(SYM_Hue), h_arg, 1);
            }
            if (strcmp(nm, "Temperature") == 0 || strcmp(nm, "Thermal") == 0) {
                double rv, gv, bv; thermal_rgb(t, &rv, &gv, &bv);
                Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
                return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
            }
        }
        /* Try f[t] then f[x,y,t] */
        static const int arities[] = { 1, 3 };
        double vals3[3] = { x, y, t };
        for (size_t ai = 0; ai < 2; ai++) {
            int ar = arities[ai];
            Expr** fargs = malloc(sizeof(Expr*) * (size_t)ar);
            if (ar == 1) fargs[0] = expr_new_real(t);
            else { fargs[0] = expr_new_real(x); fargs[1] = expr_new_real(y); fargs[2] = expr_new_real(vals3[2]); }
            Expr* fcall = expr_new_function(expr_copy(cfn), fargs, (size_t)ar);
            free(fargs);
            Expr* result = evaluate(fcall);
            bool is_color = (result->type == EXPR_FUNCTION
                             && result->data.function.head
                             && result->data.function.head->type == EXPR_SYMBOL
                             && (result->data.function.head->data.symbol == SYM_RGBColor
                              || result->data.function.head->data.symbol == SYM_GrayLevel
                              || result->data.function.head->data.symbol == SYM_Hue));
            if (is_color) return result;
            expr_free(result);
        }
    }
    /* Default thermal ramp. */
    double rv, gv, bv; thermal_rgb(t, &rv, &gv, &bv);
    Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

/* ------------------------------------------------------------------ */
/* Option parsing                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    int     plot_points;           /* grid points per axis, default 25 */
    int     n_contours;            /* number of auto-levels, default 10 */
    double* levels;                /* explicit contour values; NULL = auto */
    size_t  n_levels;
    bool    contours_set;          /* true if Contours -> ... was given explicitly */
    bool    show_legend;           /* PlotLegends -> Automatic: color bar / swatches */
    Expr*   contour_style;         /* borrowed; NULL = Automatic (color by level) */
    bool    contour_style_none;    /* ContourStyle -> None: no lines */
    bool    contour_labels;        /* draw Text[level] at midpoint */
    int     shading;               /* 0 = off, 1 = on, -1 = Automatic */
    Expr*   color_function;        /* borrowed; NULL = none */
    bool    color_function_scaling;
    Expr*   region_function;       /* borrowed; NULL = none */
} ContourOpts;

static bool split_contour_options(Expr* res, ContourOpts* co,
                                  Expr*** passthrough_out, size_t* pt_count_out) {
    co->plot_points = 25;
    co->n_contours = 10;
    co->levels = NULL;
    co->n_levels = 0;
    co->contours_set = false;
    co->show_legend = false;
    co->contour_style = NULL;
    co->contour_style_none = false;
    co->contour_labels = false;
    co->shading = -1; /* Automatic */
    co->color_function = NULL;
    co->color_function_scaling = true;
    co->region_function = NULL;

    size_t argc = res->data.function.arg_count;
    /* opts start at index 3 (after f, {x,...}, {y,...}) */
    size_t cap = (argc > 3 ? argc - 3 : 0) + 4;
    Expr** passthrough = malloc(sizeof(Expr*) * cap);
    size_t n = 0;

    bool have_axes = false, have_aspect = false, have_frame = false;

#define CP_FAIL() do { free(passthrough); free(co->levels); co->levels = NULL; return false; } while(0)

    for (size_t i = 3; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) CP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) CP_FAIL();
            co->plot_points = (int)v;
        } else if (name == SYM_Contours) {
            co->contours_set = true;
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_INTEGER && v->data.integer >= 1) {
                co->n_contours = (int)v->data.integer;
                expr_free(v);
            } else if (v->type == EXPR_FUNCTION
                       && v->data.function.head
                       && v->data.function.head->type == EXPR_SYMBOL
                       && v->data.function.head->data.symbol == SYM_List) {
                /* Explicit list of levels. */
                size_t nl = v->data.function.arg_count;
                free(co->levels);
                co->levels = malloc(sizeof(double) * (nl > 0 ? nl : 1));
                co->n_levels = 0;
                for (size_t k = 0; k < nl; k++) {
                    double lv;
                    if (numericize_bound(v->data.function.args[k], &lv))
                        co->levels[co->n_levels++] = lv;
                }
                expr_free(v);
            } else {
                /* Automatic / unknown → keep default. */
                expr_free(v);
            }
        } else if (name == SYM_ContourStyle) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL
                && (v->data.symbol == SYM_None || v->data.symbol == SYM_False)) {
                co->contour_style_none = true;
                expr_free(v);
            } else {
                expr_free(v);
                co->contour_style = rhs; /* borrow the unevaluated value */
            }
        } else if (name == SYM_ContourLabels) {
            Expr* v = evaluate(expr_copy(rhs));
            co->contour_labels = !(v->type == EXPR_SYMBOL
                                   && (v->data.symbol == SYM_False
                                    || v->data.symbol == SYM_None));
            expr_free(v);
        } else if (name == SYM_ContourShading) {
            Expr* v = evaluate(expr_copy(rhs));
            if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic) {
                co->shading = -1;
            } else if (v->type == EXPR_SYMBOL
                       && (v->data.symbol == SYM_False || v->data.symbol == SYM_None)) {
                co->shading = 0;
            } else {
                co->shading = 1;
            }
            expr_free(v);
        } else if (name == SYM_ColorFunction) {
            co->color_function = rhs; /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            co->color_function_scaling = !(v->type == EXPR_SYMBOL
                                           && v->data.symbol == SYM_False);
            expr_free(v);
        } else if (name == SYM_RegionFunction) {
            co->region_function = rhs; /* borrowed */
        } else if (name == SYM_PlotLegends) {
            /* Automatic or True → show colour bar / swatches; None/False → off. */
            Expr* v = evaluate(expr_copy(rhs));
            co->show_legend = !(v->type == EXPR_SYMBOL
                                && (v->data.symbol == SYM_None
                                    || v->data.symbol == SYM_False));
            expr_free(v);
        } else {
            /* Pass through to Graphics[...] */
            if (name == SYM_Axes)        have_axes   = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol == SYM_False || rhs->data.symbol == SYM_None)))
                    have_frame = true;
            }
            Expr* val = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    if (!have_axes && !have_frame) {
        Expr* a[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_True) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        passthrough[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *passthrough_out = passthrough;
    *pt_count_out = n;
    return true;

#undef CP_FAIL
}

/* ------------------------------------------------------------------ */
/* Grid evaluation                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr* xvar;  /* borrowed */
    Expr* yvar;  /* borrowed */
    Expr* body;  /* borrowed */
} GridCtx;

/* Evaluate f at (x, y), returning NaN on failure. */
static double eval_at(const GridCtx* ctx, double x, double y) {
    Expr* xv = expr_new_real(x);
    Expr* yv = expr_new_real(y);
    symtab_add_own_value(ctx->xvar->data.symbol, ctx->xvar, xv);
    symtab_add_own_value(ctx->yvar->data.symbol, ctx->yvar, yv);
    Expr* result = evaluate(ctx->body);
    double v;
    bool ok = expr_to_real_double(result, &v) && isfinite(v);
    expr_free(result);
    expr_free(xv);
    expr_free(yv);
    return ok ? v : NAN;
}

/* ------------------------------------------------------------------ */
/* Marching squares                                                     */
/* ------------------------------------------------------------------ */

/* Edge identifiers for a single cell. */
#define EDGE_B 0   /* bottom: (xi, yj) → (xi+dx, yj) */
#define EDGE_R 1   /* right:  (xi+dx, yj) → (xi+dx, yj+dy) */
#define EDGE_T 2   /* top:    (xi, yj+dy) → (xi+dx, yj+dy) */
#define EDGE_L 3   /* left:   (xi, yj) → (xi, yj+dy) */

typedef struct { double x, y; } Pt2;

/* Linearly interpolate the crossing point on an edge. */
static Pt2 edge_pt(int edge,
                   double xi, double yj, double dx, double dy,
                   double v00, double v10, double v11, double v01, double level) {
    double t;
    Pt2 p;
    switch (edge) {
        case EDGE_B:
            t = (v10 != v00) ? (level - v00) / (v10 - v00) : 0.5;
            t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            p.x = xi + t * dx;  p.y = yj;
            break;
        case EDGE_R:
            t = (v11 != v10) ? (level - v10) / (v11 - v10) : 0.5;
            t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            p.x = xi + dx;  p.y = yj + t * dy;
            break;
        case EDGE_T:
            t = (v11 != v01) ? (level - v01) / (v11 - v01) : 0.5;
            t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            p.x = xi + t * dx;  p.y = yj + dy;
            break;
        case EDGE_L:
            t = (v01 != v00) ? (level - v00) / (v01 - v00) : 0.5;
            t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            p.x = xi;  p.y = yj + t * dy;
            break;
        default: p.x = xi;  p.y = yj;  break;
    }
    return p;
}

/* Append a 2-point Line[{{x1,y1},{x2,y2}}] to prims[*n]. */
static void push_segment(Expr** prims, size_t* n, Pt2 a, Pt2 b) {
    Expr* ca[2] = { expr_new_real(a.x), expr_new_real(a.y) };
    Expr* cb[2] = { expr_new_real(b.x), expr_new_real(b.y) };
    Expr* pts[2] = {
        expr_new_function(expr_new_symbol(SYM_List), ca, 2),
        expr_new_function(expr_new_symbol(SYM_List), cb, 2),
    };
    Expr* plist = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
    Expr* la[1] = { plist };
    prims[(*n)++] = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
}

/* Run marching squares for one cell and emit the resulting segment(s) to
 * prims[].  Returns the count of segments appended (0, 1, or 2). */
static int cell_march(Expr** prims, size_t* n,
                      double xi, double yj, double dx, double dy,
                      double v00, double v10, double v11, double v01, double level) {
    /* Skip any cell with a NaN corner. */
    if (!isfinite(v00) || !isfinite(v10) || !isfinite(v11) || !isfinite(v01)) return 0;

    /* State: bit 0=BL, 1=BR, 2=TR, 3=TL */
    int s = (v00 >= level ? 1 : 0)
          | (v10 >= level ? 2 : 0)
          | (v11 >= level ? 4 : 0)
          | (v01 >= level ? 8 : 0);

    /* Bilinear centre for saddle disambiguation (states 5 and 10). */
    double ctr = (v00 + v10 + v11 + v01) * 0.25;

#define EP(e) edge_pt((e), xi, yj, dx, dy, v00, v10, v11, v01, level)

    switch (s) {
        case 0: case 15: return 0;
        /* Single-segment cases. */
        case  1: push_segment(prims, n, EP(EDGE_B), EP(EDGE_L)); return 1;
        case  2: push_segment(prims, n, EP(EDGE_B), EP(EDGE_R)); return 1;
        case  3: push_segment(prims, n, EP(EDGE_L), EP(EDGE_R)); return 1;
        case  4: push_segment(prims, n, EP(EDGE_R), EP(EDGE_T)); return 1;
        case  6: push_segment(prims, n, EP(EDGE_B), EP(EDGE_T)); return 1;
        case  7: push_segment(prims, n, EP(EDGE_T), EP(EDGE_L)); return 1;
        case  8: push_segment(prims, n, EP(EDGE_T), EP(EDGE_L)); return 1;
        case  9: push_segment(prims, n, EP(EDGE_B), EP(EDGE_T)); return 1;
        case 11: push_segment(prims, n, EP(EDGE_R), EP(EDGE_T)); return 1;
        case 12: push_segment(prims, n, EP(EDGE_L), EP(EDGE_R)); return 1;
        case 13: push_segment(prims, n, EP(EDGE_B), EP(EDGE_R)); return 1;
        case 14: push_segment(prims, n, EP(EDGE_B), EP(EDGE_L)); return 1;
        /* Saddle cases: bilinear centre decides which diagonal connects. */
        case  5: /* BL+TR above */
            if (ctr >= level) {
                push_segment(prims, n, EP(EDGE_B), EP(EDGE_L));
                push_segment(prims, n, EP(EDGE_R), EP(EDGE_T));
            } else {
                push_segment(prims, n, EP(EDGE_B), EP(EDGE_R));
                push_segment(prims, n, EP(EDGE_L), EP(EDGE_T));
            }
            return 2;
        case 10: /* BR+TL above */
            if (ctr >= level) {
                push_segment(prims, n, EP(EDGE_B), EP(EDGE_R));
                push_segment(prims, n, EP(EDGE_L), EP(EDGE_T));
            } else {
                push_segment(prims, n, EP(EDGE_B), EP(EDGE_L));
                push_segment(prims, n, EP(EDGE_R), EP(EDGE_T));
            }
            return 2;
    }
#undef EP
    return 0;
}

/* ------------------------------------------------------------------ */
/* Contour-level label: Text[level, midpoint_of_first_segment]         */
/* ------------------------------------------------------------------ */

static Expr* make_contour_label(double level, Pt2 a, Pt2 b) {
    /* Format the level value as a short string. */
    char buf[48];
    snprintf(buf, sizeof(buf), "%.4g", level);
    Expr* label = expr_new_string(buf);

    double mx = (a.x + b.x) * 0.5;
    double my = (a.y + b.y) * 0.5;
    Expr* pos[2] = { expr_new_real(mx), expr_new_real(my) };
    Expr* pt = expr_new_function(expr_new_symbol(SYM_List), pos, 2);
    Expr* ta[2] = { label, pt };
    return expr_new_function(expr_new_symbol(SYM_Text), ta, 2);
}

/* ------------------------------------------------------------------ */
/* Main builtin                                                         */
/* ------------------------------------------------------------------ */

Expr* builtin_contourplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 3) return NULL;

    Expr* body = res->data.function.args[0]; /* f(x,y) — held */

    /* Detect List body: ContourPlot[{eq1, eq2, ...}, {x,...}, {y,...}]
     * (handled as a separate multi-contour code path below). */
    bool is_list_body = (body->type == EXPR_FUNCTION
        && body->data.function.head
        && body->data.function.head->type == EXPR_SYMBOL
        && body->data.function.head->data.symbol == SYM_List);

    /* Support equation form: ContourPlot[lhs == rhs, ...]
     * Rewrite as body = lhs - rhs and plot the zero level set.
     * effective_body is non-NULL only when we own it and must free it. */
    Expr* effective_body = NULL;
    bool is_equation = false;
    if (!is_list_body
        && body->type == EXPR_FUNCTION
        && body->data.function.head
        && body->data.function.head->type == EXPR_SYMBOL
        && body->data.function.head->data.symbol == SYM_Equal
        && body->data.function.arg_count == 2) {
        is_equation = true;
        Expr* lhs = body->data.function.args[0];
        Expr* rhs = body->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(rhs) };
        Expr* neg_rhs = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        Expr* diff_args[2] = { expr_copy(lhs), neg_rhs };
        effective_body = expr_new_function(expr_new_symbol(SYM_Plus), diff_args, 2);
    }
    Expr* eval_body = effective_body ? effective_body : body;

    /* {x, xmin, xmax} */
    IterSpec xspec;
    if (!iter_spec_parse(res->data.function.args[1], &xspec)) return NULL;
    if (xspec.kind != ITER_KIND_RANGE || !xspec.var) { iter_spec_free(&xspec); return NULL; }
    double xmin, xmax;
    if (!numericize_bound(xspec.imin, &xmin) || !numericize_bound(xspec.imax, &xmax)
        || !(xmin < xmax)) { iter_spec_free(&xspec); return NULL; }

    /* {y, ymin, ymax} */
    IterSpec yspec;
    if (!iter_spec_parse(res->data.function.args[2], &yspec)) {
        iter_spec_free(&xspec); return NULL; }
    if (yspec.kind != ITER_KIND_RANGE || !yspec.var) {
        iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }
    double ymin, ymax;
    if (!numericize_bound(yspec.imin, &ymin) || !numericize_bound(yspec.imax, &ymax)
        || !(ymin < ymax)) { iter_spec_free(&xspec); iter_spec_free(&yspec); return NULL; }

    /* Options */
    ContourOpts co;
    Expr** passthrough = NULL;
    size_t pt_count = 0;
    if (!split_contour_options(res, &co, &passthrough, &pt_count)) {
        iter_spec_free(&xspec); iter_spec_free(&yspec);
        if (effective_body) expr_free(effective_body);
        return NULL; }

    /* Determine shading: Automatic → on if ColorFunction is set */
    bool do_shade = (co.shading == 1) || (co.shading == -1 && co.color_function != NULL);

    int N = co.plot_points;        /* grid is (N+1)×(N+1) points, N×N cells */
    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;

    /* Allocate and evaluate the grid. */
    double* grid = malloc(sizeof(double) * (size_t)(N + 1) * (size_t)(N + 1));
    if (!grid) {
        iter_spec_free(&xspec); iter_spec_free(&yspec);
        free(co.levels);
        for (size_t i = 0; i < pt_count; i++) expr_free(passthrough[i]);
        free(passthrough);
        if (effective_body) expr_free(effective_body);
        return NULL;
    }

    /* ------------------------------------------------------------------ */
    /* Multi-equation List body: ContourPlot[{eq1, eq2, ...}, ...]        */
    /* Each element is plotted as its own zero-contour with a distinct     */
    /* cycling colour.  Returns early — does not fall through to the       */
    /* single-body path below.                                             */
    /* ------------------------------------------------------------------ */
    if (is_list_body) {
        size_t neq = body->data.function.arg_count;

        /* Muted cycling palette (matches Mathematica's default ContourPlot colours). */
        static const double CLR[5][3] = {
            {0.368, 0.507, 0.710},  /* blue   */
            {0.881, 0.407, 0.341},  /* red    */
            {0.560, 0.692, 0.194},  /* green  */
            {0.728, 0.426, 0.811},  /* purple */
            {0.923, 0.694, 0.122},  /* yellow */
        };

        size_t lprim_cap = (neq > 0 ? neq : 1) * (size_t)N * (size_t)N * 3 + 16;
        Expr** lprims = malloc(sizeof(Expr*) * lprim_cap);
        size_t lnprim = 0;

#define LENSURE(extra) do { \
    if (lnprim + (extra) >= lprim_cap) { \
        lprim_cap = lprim_cap * 2 + (extra); \
        lprims = realloc(lprims, sizeof(Expr*) * lprim_cap); \
    } \
} while(0)

        for (size_t eq_i = 0; eq_i < neq; eq_i++) {
            Expr* elem = body->data.function.args[eq_i];
            Expr* sub_body = NULL;

            /* Element is lhs == rhs → evaluate lhs - rhs, seek level 0.
             * Any other expression is used verbatim (level 0 of that function). */
            if (elem->type == EXPR_FUNCTION
                && elem->data.function.head
                && elem->data.function.head->type == EXPR_SYMBOL
                && elem->data.function.head->data.symbol == SYM_Equal
                && elem->data.function.arg_count == 2) {
                Expr* elhs = elem->data.function.args[0];
                Expr* erhs = elem->data.function.args[1];
                Expr* na[2] = { expr_new_integer(-1), expr_copy(erhs) };
                Expr* nr = expr_new_function(expr_new_symbol(SYM_Times), na, 2);
                Expr* da[2] = { expr_copy(elhs), nr };
                sub_body = expr_new_function(expr_new_symbol(SYM_Plus), da, 2);
            } else {
                sub_body = expr_copy(elem);
            }

            /* Evaluate grid for this element. */
            GridCtx lctx = { .xvar = xspec.var, .yvar = yspec.var, .body = sub_body };
            Rule* lox = iter_spec_shadow(xspec.var);
            Rule* loy = iter_spec_shadow(yspec.var);
            for (int iy = 0; iy <= N; iy++) {
                double yj = ymin + iy * dy;
                for (int ix = 0; ix <= N; ix++)
                    grid[iy * (N + 1) + ix] = eval_at(&lctx, xmin + ix * dx, yj);
            }
            iter_spec_restore(xspec.var, lox);
            iter_spec_restore(yspec.var, loy);
            expr_free(sub_body);

            /* Emit colour directive for this equation. */
            const double* clr = CLR[eq_i % 5];
            Expr* ca[3] = { expr_new_real(clr[0]), expr_new_real(clr[1]),
                             expr_new_real(clr[2]) };
            LENSURE(1);
            lprims[lnprim++] = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);

            /* Determine contour levels for this element:
             * user-specified Contours override; otherwise single level at 0. */
            double  zero_lv  = 0.0;
            double* sub_lvs  = (co.levels && co.n_levels > 0) ? co.levels : &zero_lv;
            size_t  sub_nl   = (co.levels && co.n_levels > 0) ? co.n_levels : 1;

            /* Marching squares. */
            for (size_t li = 0; li < sub_nl; li++) {
                double level = sub_lvs[li];
                for (int iy = 0; iy < N; iy++) {
                    double yj = ymin + iy * dy;
                    for (int ix = 0; ix < N; ix++) {
                        double xi  = xmin + ix * dx;
                        double v00 = grid[ iy      * (N + 1) + ix    ];
                        double v10 = grid[ iy      * (N + 1) + ix + 1];
                        double v11 = grid[(iy + 1) * (N + 1) + ix + 1];
                        double v01 = grid[(iy + 1) * (N + 1) + ix    ];
                        LENSURE(4);
                        cell_march(lprims, &lnprim, xi, yj, dx, dy,
                                   v00, v10, v11, v01, level);
                    }
                }
            }
        }

#undef LENSURE

        free(grid);
        free(co.levels);

        /* Embed PlotRange if not already in passthrough. */
        {
            bool have_pr = false;
            for (size_t i = 0; i < pt_count; i++) {
                const Expr* e = passthrough[i];
                if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                    && e->data.function.args[0]->type == EXPR_SYMBOL
                    && e->data.function.args[0]->data.symbol == SYM_PlotRange)
                { have_pr = true; break; }
            }
            if (!have_pr) {
                passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
                Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
                Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
                Expr* xrng = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
                Expr* yrng = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
                Expr* rl[2] = { xrng, yrng };
                Expr* pr   = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
                Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
                passthrough[pt_count++] =
                    expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
            }
        }

        /* PlotLegends -> Automatic: one swatch per equation. */
        if (co.show_legend && neq > 0) {
            Expr** entries = malloc(sizeof(Expr*) * neq);
            for (size_t eq_i = 0; eq_i < neq; eq_i++) {
                const double* clr = CLR[eq_i % 5];
                Expr* ca[3] = { expr_new_real(clr[0]), expr_new_real(clr[1]),
                                 expr_new_real(clr[2]) };
                Expr* color = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);
                char lbuf[16];
                snprintf(lbuf, sizeof(lbuf), "%zu", eq_i + 1);
                Expr* en[2] = { color, expr_new_string(lbuf) };
                entries[eq_i] = expr_new_function(expr_new_symbol(SYM_List), en, 2);
            }
            Expr* legend_meta = expr_new_function(
                expr_new_symbol(SYM_PlotLegendData), entries, neq);
            free(entries);
            passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
            passthrough[pt_count++] = legend_meta;
        }

        Expr* plist = expr_new_function(expr_new_symbol(SYM_List), lprims, lnprim);
        free(lprims);

        size_t gargc = 1 + pt_count;
        Expr** gargs = malloc(sizeof(Expr*) * gargc);
        gargs[0] = plist;
        for (size_t i = 0; i < pt_count; i++) gargs[1 + i] = passthrough[i];
        free(passthrough);

        Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
        free(gargs);

        iter_spec_free(&xspec);
        iter_spec_free(&yspec);
        /* effective_body is always NULL for list bodies. */
        return graphics;
    }

    /* For equation bodies, default to a single contour at level 0.
     * The user can still override with an explicit Contours option. */
    if (is_equation && !co.contours_set) {
        free(co.levels);
        co.levels = malloc(sizeof(double));
        if (co.levels) { co.levels[0] = 0.0; co.n_levels = 1; }
    }

    GridCtx ctx = { .xvar = xspec.var, .yvar = yspec.var, .body = eval_body };
    Rule* old_x = iter_spec_shadow(xspec.var);
    Rule* old_y = iter_spec_shadow(yspec.var);

    double zmin = 1e300, zmax = -1e300;
    for (int iy = 0; iy <= N; iy++) {
        double yj = ymin + iy * dy;
        for (int ix = 0; ix <= N; ix++) {
            double xi = xmin + ix * dx;
            double v = eval_at(&ctx, xi, yj);
            grid[iy * (N + 1) + ix] = v;
            if (isfinite(v)) {
                if (v < zmin) zmin = v;
                if (v > zmax) zmax = v;
            }
        }
    }

    iter_spec_restore(xspec.var, old_x);
    iter_spec_restore(yspec.var, old_y);

    /* Build the contour levels. */
    double* levels;
    size_t n_levels;
    bool levels_owned = false;
    if (co.levels && co.n_levels > 0) {
        levels = co.levels;
        n_levels = co.n_levels;
    } else {
        /* Auto: N evenly spaced levels inside [zmin, zmax]. */
        int nc = co.n_contours;
        levels = malloc(sizeof(double) * (size_t)nc);
        levels_owned = true;
        n_levels = (size_t)nc;
        double zspan = (zmax > zmin) ? (zmax - zmin) : 1.0;
        for (int k = 0; k < nc; k++)
            levels[k] = zmin + (k + 1) * zspan / (nc + 1);
    }

    /* Generous upper bound on primitive count:
     *   shading:       N*N rectangles × 2 (color + Rectangle)
     *   per level:     N*N cells × 2 segments + colour + optional label
     *   scale bar:     not used here
     */
    size_t prim_cap = (size_t)N * N * 2 + n_levels * ((size_t)N * N * 3 + 4) + 8;
    Expr** prims = malloc(sizeof(Expr*) * prim_cap);
    size_t nprim = 0;

#define ENSURE_CAP(extra) do { \
    if (nprim + (extra) >= prim_cap) { \
        prim_cap = prim_cap * 2 + (extra); \
        prims = realloc(prims, sizeof(Expr*) * prim_cap); \
    } \
} while(0)

    /* ---- Shading: colour each grid cell ---- */
    if (do_shade) {
        double zspan = (zmax > zmin) ? (zmax - zmin) : 1.0;
        for (int iy = 0; iy < N; iy++) {
            double yj = ymin + iy * dy;
            for (int ix = 0; ix < N; ix++) {
                double xi = xmin + ix * dx;
                double v00 = grid[iy       * (N + 1) + ix    ];
                double v10 = grid[iy       * (N + 1) + ix + 1];
                double v11 = grid[(iy + 1) * (N + 1) + ix + 1];
                double v01 = grid[(iy + 1) * (N + 1) + ix    ];

                /* Skip cells with any NaN corner. */
                if (!isfinite(v00) || !isfinite(v10) || !isfinite(v11) || !isfinite(v01))
                    continue;

                /* RegionFunction: skip if centre is outside. */
                double cx = xi + dx * 0.5, cy = yj + dy * 0.5;
                if (co.region_function && !eval_region(co.region_function, cx, cy))
                    continue;

                double avg = (v00 + v10 + v11 + v01) * 0.25;
                double t = co.color_function_scaling
                    ? (avg - zmin) / zspan : avg;
                t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);

                ENSURE_CAP(2);
                prims[nprim++] = contour_color(co.color_function, t, cx, cy);

                /* Rectangle[{xi, yj}, {xi+dx, yj+dy}] */
                Expr* p1[2] = { expr_new_real(xi),      expr_new_real(yj) };
                Expr* p2[2] = { expr_new_real(xi + dx), expr_new_real(yj + dy) };
                Expr* ra[2] = {
                    expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                    expr_new_function(expr_new_symbol(SYM_List), p2, 2),
                };
                prims[nprim++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);
            }
        }
        /* Restore full opacity / reset colour for contour lines drawn next. */
        ENSURE_CAP(1);
        Expr* op_arg[1] = { expr_new_integer(1) };
        prims[nprim++] = expr_new_function(expr_new_symbol(SYM_Opacity), op_arg, 1);
    }

    /* ---- Contour lines ---- */
    if (!co.contour_style_none) {
        bool have_custom_style = (co.contour_style != NULL);
        /* If a single style was given, emit it once before all levels. */
        if (have_custom_style && !(co.contour_style->type == EXPR_FUNCTION
              && co.contour_style->data.function.head
              && co.contour_style->data.function.head->type == EXPR_SYMBOL
              && co.contour_style->data.function.head->data.symbol == SYM_List)) {
            ENSURE_CAP(2);
            prims[nprim++] = evaluate(expr_copy(co.contour_style));
        }

        double zspan = (zmax > zmin) ? (zmax - zmin) : 1.0;

        for (size_t li = 0; li < n_levels; li++) {
            double level = levels[li];

            /* Per-level colour: from ContourStyle list or auto thermal. */
            if (have_custom_style
                && co.contour_style->type == EXPR_FUNCTION
                && co.contour_style->data.function.head
                && co.contour_style->data.function.head->type == EXPR_SYMBOL
                && co.contour_style->data.function.head->data.symbol == SYM_List) {
                /* Cycle through the list. */
                size_t list_len = co.contour_style->data.function.arg_count;
                if (list_len > 0) {
                    ENSURE_CAP(1);
                    prims[nprim++] = evaluate(
                        expr_copy(co.contour_style->data.function.args[li % list_len]));
                }
            } else if (!have_custom_style) {
                /* Auto: colour the lines by level height (dark lines over light shading). */
                double t = co.color_function_scaling
                    ? (level - zmin) / zspan : level;
                t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
                ENSURE_CAP(2);
                /* Draw slightly darkened version of the shading colour, or
                 * plain dark gray when not shading. */
                if (do_shade) {
                    double rv, gv, bv;
                    thermal_rgb(t, &rv, &gv, &bv);
                    rv *= 0.6; gv *= 0.6; bv *= 0.6;
                    Expr* ca[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
                    prims[nprim++] = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);
                } else {
                    /* Coloured lines matching the thermal palette — nice for plain contour. */
                    double rv, gv, bv;
                    thermal_rgb(t, &rv, &gv, &bv);
                    Expr* ca[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
                    prims[nprim++] = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);
                }
            }

            /* Marching squares over every cell. */
            bool label_done = false;
            for (int iy = 0; iy < N; iy++) {
                double yj = ymin + iy * dy;
                for (int ix = 0; ix < N; ix++) {
                    double xi = xmin + ix * dx;
                    double v00 = grid[iy       * (N + 1) + ix    ];
                    double v10 = grid[iy       * (N + 1) + ix + 1];
                    double v11 = grid[(iy + 1) * (N + 1) + ix + 1];
                    double v01 = grid[(iy + 1) * (N + 1) + ix    ];

                    ENSURE_CAP(4);
                    size_t before = nprim;
                    cell_march(prims, &nprim, xi, yj, dx, dy, v00, v10, v11, v01, level);

                    /* Capture first segment for the label. */
                    if (co.contour_labels && !label_done && nprim > before) {
                        /* The first emitted segment is a Line[{{x1,y1},{x2,y2}}]; decode it. */
                        Expr* seg = prims[before];
                        if (seg && seg->type == EXPR_FUNCTION
                            && seg->data.function.arg_count == 1
                            && seg->data.function.args[0]->type == EXPR_FUNCTION
                            && seg->data.function.args[0]->data.function.arg_count == 2) {
                            Expr* plist = seg->data.function.args[0];
                            Expr* pa = plist->data.function.args[0];
                            Expr* pb = plist->data.function.args[1];
                            Pt2 a = { 0, 0 }, b = { 0, 0 };
                            if (pa->data.function.arg_count == 2) {
                                expr_to_real_double(pa->data.function.args[0], &a.x);
                                expr_to_real_double(pa->data.function.args[1], &a.y);
                            }
                            if (pb->data.function.arg_count == 2) {
                                expr_to_real_double(pb->data.function.args[0], &b.x);
                                expr_to_real_double(pb->data.function.args[1], &b.y);
                            }
                            ENSURE_CAP(1);
                            prims[nprim++] = make_contour_label(level, a, b);
                            label_done = true;
                        }
                    }
                }
            }
        }
    }

    free(grid);
    if (levels_owned) free(levels);
    free(co.levels);

    /* Always embed an explicit PlotRange so the renderer uses the user's
     * domain, matching the behaviour of StreamPlot and Plot3D. */
    {
        bool have_pr = false;
        for (size_t i = 0; i < pt_count; i++) {
            const Expr* e = passthrough[i];
            if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
                && e->data.function.args[0]->type == EXPR_SYMBOL
                && e->data.function.args[0]->data.symbol == SYM_PlotRange)
                { have_pr = true; break; }
        }
        if (!have_pr) {
            passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
            Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
            Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
            Expr* xrange = expr_new_function(expr_new_symbol(SYM_List), xr, 2);
            Expr* yrange = expr_new_function(expr_new_symbol(SYM_List), yr, 2);
            Expr* rlist[2] = { xrange, yrange };
            Expr* pr = expr_new_function(expr_new_symbol(SYM_List), rlist, 2);
            Expr* ra[2] = { expr_new_symbol(SYM_PlotRange), pr };
            passthrough[pt_count++] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
    }

    /* PlotLegends -> Automatic: embed a $StreamColorBar[zmin, zmax] metadata
     * node so the renderer draws a vertical colour-scale bar. */
    if (co.show_legend) {
        double bar_lo = isfinite(zmin) ? zmin : 0.0;
        double bar_hi = isfinite(zmax) ? zmax : 1.0;
        if (bar_lo == bar_hi) bar_hi = bar_lo + 1.0;
        passthrough = realloc(passthrough, sizeof(Expr*) * (pt_count + 1));
        Expr* cb_args[2] = { expr_new_real(bar_lo), expr_new_real(bar_hi) };
        passthrough[pt_count++] =
            expr_new_function(expr_new_symbol(SYM_StreamColorBar), cb_args, 2);
    }

    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, nprim);
    free(prims);

    size_t gargc = 1 + pt_count;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = prim_list;
    for (size_t i = 0; i < pt_count; i++) gargs[1 + i] = passthrough[i];
    free(passthrough);

    Expr* graphics = expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
    free(gargs);

    iter_spec_free(&xspec);
    iter_spec_free(&yspec);
    if (effective_body) expr_free(effective_body);
    return graphics;
}
