/* complexplot.c — ComplexPlot and ComplexPlot3D.
 *
 * Both functions share the same domain-parsing logic and option set.
 * The only structural difference is what they build from the evaluated
 * grid: ComplexPlot emits Rectangle primitives into Graphics[], and
 * ComplexPlot3D emits Polygon quads (height = |w|, colour = arg(w))
 * into Graphics3D[].
 *
 * Coloring convention: the default maps the phase arg(w) onto the cyclic
 * "Cyclic" ramp (t = (atan2(im, re) + π) / (2π) ∈ [0, 1]) and folds the
 * modulus in as an HSL lightness — zeros fade to black, poles to white — so
 * ComplexPlot[f] and ComplexPlot[f, ColorFunction -> "Cyclic"] are identical.
 * A custom ColorFunction receives the eight Mathematica arguments
 *   Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f]
 * (so #8 is the phase of the value); with ColorFunctionScaling→True (default)
 * each is scaled to [0,1] across the sampled domain.
 *
 * Both are HoldAll: the body and the iterator spec are held unevaluated
 * until z is bound to Complex[x, y] at each grid point.
 *
 * Domain spec:
 *   {z, zmin, zmax}  — z is the complex iterator variable; zmin and zmax
 *   are evaluated and their Re/Im parts define the rectangular plotting
 *   domain: xmin=Re(zmin), xmax=Re(zmax), ymin=Im(zmin), ymax=Im(zmax).
 *   Both endpoints may be real (imaginary part = 0). */

#include "complexplot.h"
#include "plot_common.h"
#include "compile/autocompile.h"   /* machine fast path, complex argument */
#include "iter.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "arithmetic.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>   /* the I in x + y I; included AFTER the headers above */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Option struct shared by both 2D and 3D                               */
/* ------------------------------------------------------------------ */

typedef struct {
    int   plot_points;
    Expr* color_function;         /* borrowed; NULL = thermal default */
    bool  color_function_scaling;
    Expr* region_function;        /* borrowed; NULL = none */
    bool  show_legend;            /* PlotLegends -> Automatic / True */
} CPlotOpts;

/* Parse trailing Rule args starting at index `first_opt_idx`.
 * Fills `o`, builds a passthrough list for Graphics/Graphics3D opts.
 * Returns true on success; on failure frees *pt_out and returns false. */
static bool split_cplot_options(Expr* res, size_t first_opt_idx, int default_plot_points,
                                 CPlotOpts* o,
                                 Expr*** pt_out, size_t* pt_n_out) {
    o->plot_points            = default_plot_points;
    o->color_function         = NULL;
    o->color_function_scaling = true;
    o->region_function        = NULL;
    o->show_legend            = false;

    size_t argc = res->data.function.arg_count;
    size_t extra = (argc > first_opt_idx ? argc - first_opt_idx : 0);
    /* +4 headroom: AspectRatio, Axes, PlotRange, $StreamColorBar defaults */
    Expr** pt = malloc(sizeof(Expr*) * (extra + 5));
    size_t n = 0;

    bool have_axes   = false;
    bool have_aspect = false;
    bool have_frame  = false;

#define CP_FAIL() do { free(pt); return false; } while (0)

    for (size_t i = first_opt_idx; i < argc; i++) {
        Expr* arg = res->data.function.args[i];
        if (!is_rule_arg(arg)) CP_FAIL();
        Expr* lhs = arg->data.function.args[0];
        Expr* rhs = arg->data.function.args[1];
        const char* name = (lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : NULL;

        if (name == SYM_PlotPoints) {
            long v;
            if (!parse_long_value(rhs, &v) || v < 2) CP_FAIL();
            o->plot_points = (int)v;
        } else if (name == SYM_ColorFunction) {
            o->color_function = rhs; /* borrowed */
        } else if (name == SYM_ColorFunctionScaling) {
            Expr* v = evaluate(expr_copy(rhs));
            o->color_function_scaling = !(v->type == EXPR_SYMBOL
                                           && v->data.symbol.name == SYM_False);
            expr_free(v);
        } else if (name == SYM_RegionFunction) {
            o->region_function = rhs; /* borrowed */
        } else if (name == SYM_PlotLegends) {
            Expr* v = evaluate(expr_copy(rhs));
            o->show_legend = !(v->type == EXPR_SYMBOL
                               && (v->data.symbol.name == SYM_None
                                   || v->data.symbol.name == SYM_False));
            expr_free(v);
        } else {
            if      (name == SYM_Axes)        have_axes   = true;
            else if (name == SYM_AspectRatio) have_aspect = true;
            else if (name == SYM_Frame) {
                if (!(rhs->type == EXPR_SYMBOL
                      && (rhs->data.symbol.name == SYM_False
                          || rhs->data.symbol.name == SYM_None)))
                    have_frame = true;
            }
            Expr* val  = evaluate(expr_copy(rhs));
            Expr* a[2] = { expr_copy(lhs), val };
            pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
        }
    }

    /* Inject defaults: raster plots default to Frame->True, Axes->False. */
    if (!have_axes && !have_frame) {
        Expr* fa[2] = { expr_new_symbol(SYM_Frame), expr_new_symbol(SYM_True) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), fa, 2);
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    } else if (!have_axes) {
        Expr* aa[2] = { expr_new_symbol(SYM_Axes), expr_new_symbol(SYM_False) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), aa, 2);
    }
    if (!have_aspect) {
        Expr* a[2] = { expr_new_symbol(SYM_AspectRatio), expr_new_integer(1) };
        pt[n++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
    }

    *pt_out   = pt;
    *pt_n_out = n;
    return true;
#undef CP_FAIL
}

/* ------------------------------------------------------------------ */
/* Domain parsing: {z, zmin, zmax}                                      */
/* ------------------------------------------------------------------ */

/* Extract (re, im) from an already-evaluated expression.
 * Accepts Complex[a,b], a pure real, or an integer; returns false otherwise. */
static bool bound_to_complex(Expr* e, double* re_out, double* im_out) {
    double r;
    if (expr_to_real_double(e, &r) && isfinite(r)) {
        *re_out = r;
        *im_out = 0.0;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double rv, iv;
        if (expr_to_real_double(re, &rv) && isfinite(rv)
            && expr_to_real_double(im, &iv) && isfinite(iv)) {
            *re_out = rv;
            *im_out = iv;
            return true;
        }
    }
    return false;
}

/* Parse iterator arg `{z, zmin, zmax}` and extract:
 *   zvar   — borrowed symbol expression from the iterator
 *   xmin/xmax — Re(zmin), Re(zmax)
 *   ymin/ymax — Im(zmin), Im(zmax)
 * Returns false on any parse / numeric failure. */
static bool parse_complex_iterator(Expr* iter, Expr** zvar_out,
                                    double* xmin, double* xmax,
                                    double* ymin, double* ymax) {
    if (!iter || iter->type != EXPR_FUNCTION
        || !iter->data.function.head
        || iter->data.function.head->type != EXPR_SYMBOL
        || iter->data.function.head->data.symbol.name != SYM_List
        || iter->data.function.arg_count != 3)
        return false;

    Expr* zvar = iter->data.function.args[0];
    if (!zvar || zvar->type != EXPR_SYMBOL) return false;

    Expr* e_zmin = evaluate(expr_copy(iter->data.function.args[1]));
    Expr* e_zmax = evaluate(expr_copy(iter->data.function.args[2]));

    double x0, y0, x1, y1;
    bool ok = bound_to_complex(e_zmin, &x0, &y0)
           && bound_to_complex(e_zmax, &x1, &y1);
    expr_free(e_zmin);
    expr_free(e_zmax);

    if (!ok) return false;
    if (x0 == x1 || y0 == y1) return false; /* degenerate domain */

    *zvar_out = zvar;   /* borrowed from res — do not free */
    *xmin = (x0 < x1 ? x0 : x1);
    *xmax = (x0 < x1 ? x1 : x0);
    *ymin = (y0 < y1 ? y0 : y1);
    *ymax = (y0 < y1 ? y1 : y0);
    return true;
}

/* ------------------------------------------------------------------ */
/* Grid evaluation                                                      */
/* ------------------------------------------------------------------ */

/* Bind z = Complex[x, y], evaluate body, extract (re, im) of result.
 * Returns true and sets *re_out, *im_out on success; returns false on
 * failure (unevaluated, non-numeric, or infinite result). */
static bool cp_eval(const AutoCompiled* ac, Expr* zvar, Expr* body, double x, double y,
                    double* re_out, double* im_out) {
    /* Compiled f(z) with a genuinely COMPLEX argument — the only sampler here
     * whose variable ranges over the plane, so it needs autocompile_new_z rather
     * than a real-input program.  A decline falls through to the interpreter,
     * which also covers every head that has a real kernel but no complex one. */
    if (ac) {
        double _Complex z = x + y * (double _Complex)I, w;
        if (autocompiled_eval_z(ac, &z, &w)) {
            double wr = creal(w), wi = cimag(w);
            /* A pole makes the compiled program yield inf/nan.  Guard here just
             * as the interpreter path below does (ComplexInfinity → invalid):
             * without it the cell is marked valid with |w| = inf and poisons
             * the min/max ranges used for ColorFunctionScaling. */
            if (isfinite(wr) && isfinite(wi)) {
                *re_out = wr; *im_out = wi;
                return true;
            }
            return false;
        }
    }
    Expr* ra[2] = { expr_new_real(x), expr_new_real(y) };
    Expr* zval  = expr_new_function(expr_new_symbol(SYM_Complex), ra, 2);
    symtab_add_own_value(zvar->data.symbol.name, zvar, zval);

    Expr* result = evaluate(expr_copy(body));
    bool ok = false;
    double r = 0.0, i = 0.0;

    double rv;
    Expr *rp, *ip;
    if (expr_to_real_double(result, &rv) && isfinite(rv)) {
        r = rv; i = 0.0; ok = true;
    } else if (is_complex(result, &rp, &ip)) {
        double rv2, iv;
        if (expr_to_real_double(rp, &rv2) && isfinite(rv2)
            && expr_to_real_double(ip, &iv) && isfinite(iv)) {
            r = rv2; i = iv; ok = true;
        }
    }
    expr_free(result);
    expr_free(zval);

    if (ok) { *re_out = r; *im_out = i; }
    return ok;
}

/* ------------------------------------------------------------------ */
/* Coloring                                                             */
/* ------------------------------------------------------------------ */

/* cp_ramp_color — a phase-keyed named ramp × modulus brightness. This is the
 * ONE path used by both the default (ramp "Cyclic") and any explicit string
 * ColorFunction, so `ComplexPlot[f]` and `ComplexPlot[f, ColorFunction ->
 * "Cyclic"]` are byte-for-byte identical (previously the default multiplied in
 * the modulus brightness but a string ramp did not, so the two disagreed).
 *
 * t = (atan2(im, re) + π) / (2π) — the normalised argument, wrapping in [0, 1];
 * a cyclic ramp (e.g. "Cyclic") therefore has no seam at arg = ±π.
 *
 * The modulus is folded in as an HSL-style *lightness*, not a plain brightness
 * multiplier: L = |w|/(1+|w|) ∈ [0, 1), which is exactly 1/2 at |w| = 1. Below
 * 1/2 the saturated hue fades toward BLACK (zeros → black); above 1/2 it fades
 * toward WHITE (poles → white), matching Mathematica's domain colouring. The
 * previous code multiplied the hue by L, so |w| ≈ 1 rendered at only ~50 %
 * brightness and a bounded function (|f| ≈ 1 over most of the plane) came out
 * uniformly dark. Returns NULL if `name` is not a recognised ramp. */
static Expr* cp_ramp_color(const char* name, double re, double im) {
    double arg = atan2(im, re);
    double t   = (arg + M_PI) / (2.0 * M_PI);
    if (t < 0.0) t = 0.0;   /* clamp atan2's (-π, π] edge cases */
    if (t > 1.0) t = 1.0;
    double rv, gv, bv;
    if (!resolve_ramp_to_rgb(name, t, &rv, &gv, &bv)) return NULL;
    double mod = sqrt(re * re + im * im);
    double L   = mod / (1.0 + mod);     /* lightness; 1/2 at |w| = 1 */
    if (L <= 0.5) {
        double s = L / 0.5;             /* 0 at |w|=0 → 1 at |w|=1 (toward black) */
        rv *= s; gv *= s; bv *= s;
    } else {
        double s = (L - 0.5) / 0.5;     /* 0 at |w|=1 → 1 at |w|→∞ (toward white) */
        rv += (1.0 - rv) * s;
        gv += (1.0 - gv) * s;
        bv += (1.0 - bv) * s;
    }
    Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

static bool is_color_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (!(h == SYM_RGBColor || h == SYM_GrayLevel || h == SYM_Hue || h == SYM_CMYKColor))
        return false;
    /* Every component must be a real number to render.  This rejects e.g.
     * Hue[0.5 + #8] left half-symbolic by a ColorFunction that referenced a
     * slot we did not supply — which must fall through to the fallback rather
     * than be handed to the renderer as a "colour". */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        double d;
        if (!expr_to_real_double(e->data.function.args[i], &d) || !isfinite(d))
            return false;
    }
    return true;
}

/* The eight arguments Mathematica supplies to a ComplexPlot ColorFunction, in
 * order: Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f] — where
 * z = x + i·y is the sample point and f = re + i·im is the value there.  So a
 * user's ColorFunction -> (Hue[#8 + 0.5] &) colours by the phase Arg[f]. */
static void cf_eight(double x, double y, double re, double im, double a[8]) {
    a[0] = x;                        /* Re[z]  */
    a[1] = y;                        /* Im[z]  */
    a[2] = sqrt(x * x + y * y);      /* Abs[z] */
    a[3] = atan2(y, x);              /* Arg[z] */
    a[4] = re;                       /* Re[f]  */
    a[5] = im;                       /* Im[f]  */
    a[6] = sqrt(re * re + im * im);  /* Abs[f] */
    a[7] = atan2(im, re);            /* Arg[f] */
}

/* Per-argument [min,max] of the eight cf_eight values across the sampled grid,
 * used to scale them to [0,1] when ColorFunctionScaling→True. */
typedef struct { double lo[8], hi[8]; } CFRange;

/* Resolve color for one grid cell from a custom ColorFunction or the default.
 * (x, y) is the sample point in the z-plane; (re, im) the value f(z) there.
 * A custom function receives the eight cf_eight arguments, each scaled to [0,1]
 * over the grid (via `cfr`) when `scaling` is on — Mathematica's default. */
static Expr* cp_color(Expr* cfn, bool scaling, const CFRange* cfr,
                      double x, double y, double re, double im) {
    /* Default: the "Cyclic" phase ramp — the SAME cp_ramp_color path an explicit
     * ColorFunction -> "Cyclic" takes, so the two render identically. */
    if (!cfn) return cp_ramp_color("Cyclic", re, im);

    if (cfn->type == EXPR_STRING) {
        /* "PhaseRings" needs both re and im — intercept before the 1-D ramp path. */
        if (strcmp(cfn->data.string, "PhaseRings") == 0) {
            double rv, gv, bv;
            phase_rings_rgb(re, im, &rv, &gv, &bv);
            Expr* a[3] = { expr_new_real(rv), expr_new_real(gv), expr_new_real(bv) };
            return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
        }
        /* Any other named ramp: phase-keyed hue × modulus brightness, matching
         * the default. (Domain colouring wants the magnitude to read, so every
         * string ramp gets the brightness — not just the default.) */
        Expr* c = cp_ramp_color(cfn->data.string, re, im);
        if (c) return c;
    }

    /* Custom function: supply the eight Mathematica arguments, scaled to [0,1]
     * across the grid when ColorFunctionScaling is on. */
    double raw[8], sc[8];
    cf_eight(x, y, re, im, raw);
    for (int k = 0; k < 8; k++) {
        sc[k] = raw[k];
        if (scaling && cfr) {
            double span = cfr->hi[k] - cfr->lo[k];
            sc[k] = (span > 0.0) ? (raw[k] - cfr->lo[k]) / span : 0.0;
            if (sc[k] < 0.0) sc[k] = 0.0;
            if (sc[k] > 1.0) sc[k] = 1.0;
        }
    }

    Expr* a8[8];
    for (int k = 0; k < 8; k++) a8[k] = expr_new_real(sc[k]);
    Expr* call8 = expr_new_function(expr_copy(cfn), a8, 8);
    Expr* r8    = evaluate(call8);
    expr_free(call8);
    if (is_color_head(r8)) return r8;
    expr_free(r8);

    /* Fallback for a one-argument colour map (e.g. a ColorData gradient handed
     * in directly): key it to the phase Arg[f], the dominant complex-plot
     * colouring variable. */
    Expr* a1[1] = { expr_new_real(sc[7]) };
    Expr* call1 = expr_new_function(expr_copy(cfn), a1, 1);
    Expr* r1    = evaluate(call1);
    expr_free(call1);
    if (is_color_head(r1)) return r1;
    expr_free(r1);

    /* Fallback: neutral mid-gray */
    Expr* ga[1] = { expr_new_real(0.5) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), ga, 1);
}

/* ------------------------------------------------------------------ */
/* Grid storage                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    double re, im;  /* evaluated complex result */
    bool   valid;   /* evaluation succeeded and RegionFunction passed */
} CGrid;

/* Height-cap for ComplexPlot3D: sort valid |f(z)| values and return the
 * 95th-percentile so that pole spikes are clipped to a sensible ceiling
 * (matching Mathematica's automatic PlotRange clipping behaviour).
 * A cap of exactly 0 is replaced by 1 to avoid a degenerate surface. */
static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

static double compute_height_cap(const CGrid* grid, size_t total) {
    /* Collect valid moduli */
    double* hs = malloc(sizeof(double) * total);
    if (!hs) return 1.0;
    size_t nc = 0;
    for (size_t k = 0; k < total; k++) {
        if (!grid[k].valid) continue;
        double h = sqrt(grid[k].re * grid[k].re + grid[k].im * grid[k].im);
        if (isfinite(h)) hs[nc++] = h;
    }
    if (nc == 0) { free(hs); return 1.0; }
    qsort(hs, nc, sizeof(double), cmp_double);
    /* 95th percentile clips the top 5% of heights (pole spikes) */
    size_t idx = (size_t)(nc * 0.95);
    if (idx >= nc) idx = nc - 1;
    double cap = hs[idx];
    free(hs);
    return (cap > 1e-10) ? cap : 1.0;
}

static CGrid* build_cgrid(Expr* zvar, Expr* body, Expr* region_fn,
                           double xmin, double xmax,
                           double ymin, double ymax, int N) {
    CGrid* grid = malloc(sizeof(CGrid) * (size_t)(N + 1) * (size_t)(N + 1));
    if (!grid) return NULL;

    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;

    const Expr* zv[1] = { zvar };
    AutoCompiled* ac = autocompile_new_z(body, zv, 1);

    /* Sampling a function with poles hits 1/0 at grid points on top of a pole
     * (z = 0 for (z^3-3)/z; z = ±i for 1/(z^2+1)).  Those cells are simply
     * dropped (invalid), so mute the informational Power::infy / Infinity::indet
     * chatter the interpreter would otherwise print for every such point —
     * exactly as Plot and the numeric optimizers do around divergent probes. */
    arith_warnings_mute_push();
    for (int iy = 0; iy <= N; iy++) {
        double y = ymin + iy * dy;
        if (iy == N) y = ymax;
        for (int ix = 0; ix <= N; ix++) {
            double x = xmin + ix * dx;
            if (ix == N) x = xmax;
            CGrid* p = &grid[iy * (N + 1) + ix];
            double re, im;
            bool ok = cp_eval(ac, zvar, body, x, y, &re, &im);
            if (ok && region_fn && !eval_region(region_fn, x, y)) ok = false;
            p->re    = ok ? re : 0.0;
            p->im    = ok ? im : 0.0;
            p->valid = ok;
        }
    }
    arith_warnings_mute_pop();
    autocompiled_free(ac);
    return grid;
}

/* Per-argument [min,max] of the eight ColorFunction arguments across valid grid
 * points, for ColorFunctionScaling→True.  Computed from the sampled grid points
 * (the corners the cell colours are averaged from), which bound the cell-centre
 * values; scaled results are clamped to [0,1] regardless. */
static void grid_cfrange(const CGrid* grid, int N,
                         double xmin, double xmax, double ymin, double ymax,
                         CFRange* r) {
    for (int k = 0; k < 8; k++) { r->lo[k] = 1e300; r->hi[k] = -1e300; }
    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;
    size_t stride = (size_t)(N + 1);
    for (int iy = 0; iy <= N; iy++) {
        double y = ymin + iy * dy;
        for (int ix = 0; ix <= N; ix++) {
            const CGrid* p = &grid[(size_t)iy * stride + (size_t)ix];
            if (!p->valid) continue;
            double a[8];
            cf_eight(xmin + ix * dx, y, p->re, p->im, a);
            for (int k = 0; k < 8; k++) {
                if (a[k] < r->lo[k]) r->lo[k] = a[k];
                if (a[k] > r->hi[k]) r->hi[k] = a[k];
            }
        }
    }
    for (int k = 0; k < 8; k++)
        if (r->lo[k] > r->hi[k]) { r->lo[k] = 0.0; r->hi[k] = 1.0; }
}

/* ------------------------------------------------------------------ */
/* Embed PlotRange matching the domain                                  */
/* ------------------------------------------------------------------ */

static void embed_plot_range(double xmin, double xmax, double ymin, double ymax,
                              Expr*** pt, size_t* pt_n) {
    /* Check if caller already supplied a PlotRange rule */
    for (size_t i = 0; i < *pt_n; i++) {
        const Expr* e = (*pt)[i];
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
            && e->data.function.args[0]->type == EXPR_SYMBOL
            && e->data.function.args[0]->data.symbol.name == SYM_PlotRange)
            return;
    }
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
    Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
    Expr* rl[2] = { expr_new_function(expr_new_symbol(SYM_List), xr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), yr, 2) };
    Expr* pr    = expr_new_function(expr_new_symbol(SYM_List), rl, 2);
    Expr* a[2]  = { expr_new_symbol(SYM_PlotRange), pr };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* 3D variant: x/y domain + explicit z-range so the renderer frames the
 * surface correctly when pole heights are clamped. */
static void embed_plot_range3(double xmin, double xmax,
                               double ymin, double ymax,
                               double zmin, double zmax,
                               Expr*** pt, size_t* pt_n) {
    for (size_t i = 0; i < *pt_n; i++) {
        const Expr* e = (*pt)[i];
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
            && e->data.function.args[0]->type == EXPR_SYMBOL
            && e->data.function.args[0]->data.symbol.name == SYM_PlotRange)
            return;
    }
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* xr[2] = { expr_new_real(xmin), expr_new_real(xmax) };
    Expr* yr[2] = { expr_new_real(ymin), expr_new_real(ymax) };
    Expr* zr[2] = { expr_new_real(zmin), expr_new_real(zmax) };
    Expr* rl[3] = { expr_new_function(expr_new_symbol(SYM_List), xr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), yr, 2),
                    expr_new_function(expr_new_symbol(SYM_List), zr, 2) };
    Expr* pr   = expr_new_function(expr_new_symbol(SYM_List), rl, 3);
    Expr* a[2] = { expr_new_symbol(SYM_PlotRange), pr };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

/* Append $StreamColorBar[-π, π, cfn_or_"Cyclic"] so the renderer draws a
 * vertical phase-angle color scale.  The bar parameter t ∈ [0,1] maps to
 * arg = -π + t·2π, identical to our normalization t = (arg+π)/(2π). The default
 * (cfn == NULL) passes the "Cyclic" named ramp explicitly so the bar's
 * resolve_ramp_to_rgb("Cyclic", t) matches cp_default_color's cyclic_phase_rgb —
 * NOT Automatic, which the renderer would resolve to the Viridis default and
 * drift from the cells. */
static void emit_phase_color_bar(Expr* cfn, Expr*** pt, size_t* pt_n) {
    *pt = realloc(*pt, sizeof(Expr*) * (*pt_n + 1));
    Expr* cfn_copy = cfn ? expr_copy(cfn) : expr_new_string("Cyclic");
    Expr* cb_args[3] = { expr_new_real(-M_PI), expr_new_real(M_PI), cfn_copy };
    (*pt)[(*pt_n)++] = expr_new_function(expr_new_symbol(SYM_StreamColorBar), cb_args, 3);
}

/* ------------------------------------------------------------------ */
/* builtin_complexplot                                                  */
/* ------------------------------------------------------------------ */

Expr* builtin_complexplot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* body = res->data.function.args[0]; /* held */

    /* Parse {z, zmin, zmax} */
    Expr*  zvar;
    double xmin, xmax, ymin, ymax;
    if (!parse_complex_iterator(res->data.function.args[1],
                                 &zvar, &xmin, &xmax, &ymin, &ymax))
        return NULL;

    CPlotOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    /* 2D raster: default 400 grid points/axis. The plane is drawn as a grid of
     * coloured cells, so the cell pitch is the on-screen resolution; 200 was
     * visibly pixelated near zeros/poles where the phase turns fastest. The
     * body is autocompiled per grid (autocompile_new_z), so 400x400 stays well
     * under ~0.5s. */
    if (!split_cplot_options(res, 2, 400, &opts, &pt, &pt_n)) return NULL;

    int N = opts.plot_points;
    size_t stride = (size_t)(N + 1);

    Rule* old_z = iter_spec_shadow(zvar);
    CGrid* grid = build_cgrid(zvar, body, opts.region_function,
                               xmin, xmax, ymin, ymax, N);
    iter_spec_restore(zvar, old_z);

    if (!grid) {
        for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]);
        free(pt);
        return NULL;
    }

    /* Per-argument ranges for ColorFunctionScaling (only a custom ColorFunction
     * consumes them; the default and string ramps key off the raw phase). */
    CFRange cfr = {{0}, {0}};
    if (opts.color_function)
        grid_cfrange(grid, N, xmin, xmax, ymin, ymax, &cfr);

    double dx = (xmax - xmin) / N;
    double dy = (ymax - ymin) / N;

    /* 2 primitives per cell: color directive + Rectangle */
    size_t cap   = (size_t)N * (size_t)N * 2 + 2;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (int iy = 0; iy < N; iy++) {
        double y0 = ymin + iy * dy;
        for (int ix = 0; ix < N; ix++) {
            double x0 = xmin + ix * dx;

            /* Average the four corner values for smooth appearance */
            const CGrid* p00 = &grid[ iy      * (int)stride + ix    ];
            const CGrid* p10 = &grid[ iy      * (int)stride + ix + 1];
            const CGrid* p11 = &grid[(iy + 1) * (int)stride + ix + 1];
            const CGrid* p01 = &grid[(iy + 1) * (int)stride + ix    ];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid)
                continue;

            double re_avg = (p00->re + p10->re + p11->re + p01->re) * 0.25;
            double im_avg = (p00->im + p10->im + p11->im + p01->im) * 0.25;

            /* Sample point at the cell centre — z the ColorFunction sees. */
            double xc = x0 + 0.5 * dx;
            double yc = y0 + 0.5 * dy;
            prims[np++] = cp_color(opts.color_function, opts.color_function_scaling,
                                    &cfr, xc, yc, re_avg, im_avg);

            /* Rectangle in plot coordinates (x = Re axis, y = Im axis).
             * The far corner overlaps one full cell into the +x/+y neighbours
             * (drawn later, so they overdraw the overlap and each cell still
             * shows its own colour), clamped at the plot edges. This closes the
             * sub-pixel seams that otherwise leave the white background showing
             * as thin lines between adjacent Rectangle fills. */
            double x1 = x0 + 2.0 * dx; if (x1 > xmax) x1 = xmax;
            double y1 = y0 + 2.0 * dy; if (y1 > ymax) y1 = ymax;
            Expr* p1[2] = { expr_new_real(x0), expr_new_real(y0) };
            Expr* p2[2] = { expr_new_real(x1), expr_new_real(y1) };
            Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), p1, 2),
                             expr_new_function(expr_new_symbol(SYM_List), p2, 2) };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Rectangle), ra, 2);
        }
    }

    free(grid);

    embed_plot_range(xmin, xmax, ymin, ymax, &pt, &pt_n);
    if (opts.show_legend) emit_phase_color_bar(opts.color_function, &pt, &pt_n);

    Expr* plist  = expr_new_function(expr_new_symbol(SYM_List), prims, np);
    free(prims);

    size_t gargc = 1 + pt_n;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = plist;
    for (size_t i = 0; i < pt_n; i++) gargs[1 + i] = pt[i];
    free(pt);

    return expr_new_function(expr_new_symbol(SYM_Graphics), gargs, gargc);
}

/* ------------------------------------------------------------------ */
/* builtin_complexplot3d                                                */
/* ------------------------------------------------------------------ */

static Expr* point3d(double x, double y, double z) {
    Expr* a[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    return expr_new_function(expr_new_symbol(SYM_List), a, 3);
}

Expr* builtin_complexplot3d(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* body = res->data.function.args[0]; /* held */

    Expr*  zvar;
    double xmin, xmax, ymin, ymax;
    if (!parse_complex_iterator(res->data.function.args[1],
                                 &zvar, &xmin, &xmax, &ymin, &ymax))
        return NULL;

    CPlotOpts opts;
    Expr** pt = NULL; size_t pt_n = 0;
    /* 3D: default passthrough gets Axes → True (reusing split_cplot_options).
     * Keep 200 here — this samples an N×N polygon *mesh*, not a 2D raster, so a
     * 400² grid would be 320k triangles baked into the 3D scene for no
     * comparable gain (the surface is already smoothly shaded). */
    if (!split_cplot_options(res, 2, 200, &opts, &pt, &pt_n)) return NULL;

    int N = opts.plot_points;
    size_t stride = (size_t)(N + 1);

    Rule* old_z = iter_spec_shadow(zvar);
    CGrid* grid = build_cgrid(zvar, body, opts.region_function,
                               xmin, xmax, ymin, ymax, N);
    iter_spec_restore(zvar, old_z);

    if (!grid) {
        for (size_t i = 0; i < pt_n; i++) expr_free(pt[i]);
        free(pt);
        return NULL;
    }

    /* Per-argument ranges for ColorFunctionScaling (custom ColorFunction only). */
    CFRange cfr = {{0}, {0}};
    if (opts.color_function)
        grid_cfrange(grid, N, xmin, xmax, ymin, ymax, &cfr);

    /* Clip heights at the 95th-percentile of |f(z)| so poles appear as
     * flat-topped cylinders rather than infinite spikes (Mathematica style). */
    double hcap = compute_height_cap(grid, stride * stride);

    /* 2 primitives per cell: color + Polygon */
    size_t cap   = (size_t)N * (size_t)N * 2 + 2;
    Expr** prims = malloc(sizeof(Expr*) * cap);
    size_t np    = 0;

    for (int iy = 0; iy < N; iy++) {
        for (int ix = 0; ix < N; ix++) {
            const CGrid* p00 = &grid[ iy      * (int)stride + ix    ];
            const CGrid* p10 = &grid[ iy      * (int)stride + ix + 1];
            const CGrid* p11 = &grid[(iy + 1) * (int)stride + ix + 1];
            const CGrid* p01 = &grid[(iy + 1) * (int)stride + ix    ];

            if (!p00->valid || !p10->valid || !p11->valid || !p01->valid)
                continue;

            /* Heights: |f(z)| clamped to hcap so poles become flat cylinders. */
#define HCLAMP(v) ((v) > hcap ? hcap : (v))
            double h00 = HCLAMP(sqrt(p00->re*p00->re + p00->im*p00->im));
            double h10 = HCLAMP(sqrt(p10->re*p10->re + p10->im*p10->im));
            double h11 = HCLAMP(sqrt(p11->re*p11->re + p11->im*p11->im));
            double h01 = HCLAMP(sqrt(p01->re*p01->re + p01->im*p01->im));
#undef HCLAMP

            /* Color at cell center */
            double re_avg = (p00->re + p10->re + p11->re + p01->re) * 0.25;
            double im_avg = (p00->im + p10->im + p11->im + p01->im) * 0.25;

            /* For 3D, the default color uses arg of the center but without
             * modulus-brightness attenuation (same visual weight as Plot3D). */
            Expr* color;
            if (!opts.color_function) {
                double arg = atan2(im_avg, re_avg);
                double t   = (arg + M_PI) / (2.0 * M_PI);
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
                double r, g, b;
                thermal_rgb(t, &r, &g, &b);
                Expr* ca[3] = { expr_new_real(r), expr_new_real(g), expr_new_real(b) };
                color = expr_new_function(expr_new_symbol(SYM_RGBColor), ca, 3);
            } else {
                double xc = xmin + (ix + 0.5) * (xmax - xmin) / N;
                double yc = ymin + (iy + 0.5) * (ymax - ymin) / N;
                color = cp_color(opts.color_function, opts.color_function_scaling,
                                  &cfr, xc, yc, re_avg, im_avg);
            }
            prims[np++] = color;

            /* x axis = Re(z), y axis = Im(z), z axis = |f(z)| */
            double x0 = xmin + ix       * (xmax - xmin) / N;
            double x1 = xmin + (ix + 1) * (xmax - xmin) / N;
            double y0 = ymin + iy       * (ymax - ymin) / N;
            double y1 = ymin + (iy + 1) * (ymax - ymin) / N;

            Expr* verts[4] = { point3d(x0, y0, h00),
                               point3d(x1, y0, h10),
                               point3d(x1, y1, h11),
                               point3d(x0, y1, h01) };
            Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, 4);
            Expr* pargs[1] = { vlist };
            prims[np++] = expr_new_function(expr_new_symbol(SYM_Polygon), pargs, 1);
        }
    }

    free(grid);

    /* PlotRange for 3D: supply explicit z-extent so the renderer frames
     * the surface at the clamped ceiling rather than the raw pole height. */
    embed_plot_range3(xmin, xmax, ymin, ymax, 0.0, hcap, &pt, &pt_n);
    if (opts.show_legend) emit_phase_color_bar(opts.color_function, &pt, &pt_n);

    Expr* plist  = expr_new_function(expr_new_symbol(SYM_List), prims, np);
    free(prims);

    size_t gargc = 1 + pt_n;
    Expr** gargs = malloc(sizeof(Expr*) * gargc);
    gargs[0] = plist;
    for (size_t i = 0; i < pt_n; i++) gargs[1 + i] = pt[i];
    free(pt);

    return expr_new_function(expr_new_symbol(SYM_Graphics3D), gargs, gargc);
}
