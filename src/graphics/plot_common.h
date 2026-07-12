/* plot_common.h — numeric/option helpers shared by Plot (plot.c) and
 * Plot3D (plot3d.c).
 *
 * Pure Expr/evaluator-level helpers, no Raylib dependency: option-value
 * coercion, the RegionFunction and ColorFunction evaluation idioms, and
 * the multi-curve/multi-surface palette. Keeping these in one place is
 * what lets Plot3D reuse Plot's option semantics verbatim instead of
 * re-implementing them. */
#ifndef MATHILDA_GRAPHICS_PLOT_COMMON_H
#define MATHILDA_GRAPHICS_PLOT_COMMON_H

#include "expr.h"
#include "print.h"
#include <stdbool.h>
#include <stddef.h>

/* Coerce a literal numeric Expr (Integer/Real/BigInt/MPFR/Rational) to a
 * double. Returns false for anything else (including unevaluated symbolic
 * forms -- see numericize_bound for those). */
bool expr_to_real_double(const Expr* e, double* out);

/* Numericize a possibly symbolic-but-numeric bound (2 Pi, E, Sqrt[2], ...)
 * via N[], exactly as a user typing N[expr] would. Returns false if the
 * result isn't a finite real. */
bool numericize_bound(Expr* e, double* out);

/* True if `e` is Rule[_,_] or RuleDelayed[_,_] (a trailing opts... arg). */
bool is_rule_arg(Expr* e);

/* Evaluate `rhs` and require it to be a plain machine integer. */
bool parse_long_value(Expr* rhs, long* out);

/* Distinct, harmonious per-curve/per-surface colours for multi-function
 * plots (Mathematica's ColorData[97] palette), cycled when there are more
 * curves/surfaces than palette entries. Caller owns the returned RGBColor[]
 * Expr. */
Expr* palette_color(size_t i);

/* RegionFunction: f[x,y] (2-arg) or f[x] (1-arg), tried in that order;
 * neither resolving to True/False is treated as "outside the region". */
bool eval_region(Expr* region_fn, double x, double y);

/* Resolves ColorFunction at one sampled point to a concrete color literal
 * Expr (caller owns). "Rainbow" is a built-in Hue sweep over scaled x; a
 * custom function is tried 2-arg (xscaled,yscaled) then 1-arg (xscaled).
 * Falls back to a neutral gray if nothing resolves to a recognized color
 * literal. */
Expr* eval_color_function(Expr* color_fn, double x, double y,
                           double xmin, double xmax, bool scaling);

/* 3D variant: "Rainbow" sweeps over the scaled z (height) instead of x.
 * Custom functions are tried as f[xs,ys,zs] → f[xs,zs] → f[zs] in order.
 * Falls back to neutral gray if nothing resolves to a recognised color. */
Expr* eval_color_function3(Expr* color_fn,
                            double x,    double y,    double z,
                            double xmin, double xmax,
                            double ymin, double ymax,
                            double zmin, double zmax,
                            bool scaling);

/* Named color ramps — all take t ∈ [0,1], write r/g/b ∈ [0,1]. */

/* Thermal: dark blue-purple (t=0) → purple → red → orange → bright yellow (t=1).
 * Matches Mathematica's default StreamPlot speed colormap. */
void thermal_rgb(double t, double* r, double* g, double* b);

/* CoolTones: near-white ice blue (t=0) → sky blue → cornflower → deep navy (t=1). */
void cool_tones_rgb(double t, double* r, double* g, double* b);

/* WarmTones: pale cream (t=0) → amber → orange → deep crimson (t=1). */
void warm_tones_rgb(double t, double* r, double* g, double* b);

/* Resolve a ColorFunction name string + t ∈ [0,1] to a color Expr (caller
 * owns).  Recognised names: "Rainbow", "Temperature"/"Thermal",
 * "CoolTones"/"Cool", "WarmTones"/"Warm",
 * "Greyscale"/"Grayscale"/"Grey"/"Gray".
 * Returns NULL when the name is not recognised. */
Expr* named_color_ramp(const char* name, double t);

/* resolve_ramp_to_rgb — same lookup as named_color_ramp but writes raw RGB
 * doubles [0,1] instead of constructing an Expr.  Returns 1 on success,
 * 0 if the name is not recognised. */
int resolve_ramp_to_rgb(const char* name, double t, double* r, double* g, double* b);

/* ---------------------------------------------------------------------- */
/* Axis scaling (ScalingFunctions option)                                  */
/* ---------------------------------------------------------------------- */

/* Identifies one axis's scaling transform.  SF_NONE is the identity.
 * Values are stored as EXPR_INTEGER inside $ScalingMeta[] metadata so
 * the renderer can read them without any string parsing overhead. */
typedef enum {
    SF_NONE    = 0,
    SF_LOG     = 1,   /* natural log: world = ln(data)  */
    SF_LOG2    = 2,   /* log base 2:  world = log2(data) */
    SF_LOG10   = 3,   /* log base 10: world = log10(data) */
    SF_REVERSE = 4    /* mirror axis: world = -data */
} ScaleFnType;

/* Map a data-space coordinate to world space. */
double scale_apply(ScaleFnType sf, double x);

/* Inverse: world → data space. */
double scale_invert(ScaleFnType sf, double w);

/* Parse a ScalingFunctions spec expression (string or None/Automatic) to
 * a ScaleFnType.  NULL or unrecognised → SF_NONE. */
ScaleFnType parse_scale_fn(Expr* e);

/* Parse ScalingFunctions RHS into (sf_x, sf_y).
 * "Log" → both axes Log; {"Log","Log10"} → per-axis; None/Automatic → SF_NONE. */
void parse_scaling_functions(Expr* rhs, ScaleFnType* sf_x, ScaleFnType* sf_y);

/* Append the $ScalingMeta[sfx, sfy] metadata node to a Graphics options
 * array (pt/pt_n) when at least one axis has a non-identity scale.
 * Reallocates *pt by +1 when needed; caller must pass the current capacity
 * (the array has capacity to hold at least pt_n+1 entries already) or set
 * cap > pt_n to indicate there is room.  Pass cap=0 to always reallocate. */
void emit_scaling_meta(ScaleFnType sf_x, ScaleFnType sf_y,
                       Expr*** pt, size_t* pt_n);

/* ---------------------------------------------------------------------- */
/* Build the $PlotLegendData[{color1,label1}, ...] metadata node that the
 * renderer reads to draw a legend box.  `legends` is the already-evaluated
 * PlotLegends option value; `bodies` / `nfun` supply per-curve body exprs
 * (used as auto-labels for Automatic / "Expressions"); `single_color` is the
 * resolved PlotStyle color for a single-curve plot (NULL → palette_color(0)).
 * Returns NULL if no legend should be drawn (legends is NULL or None). */
Expr* build_legend_meta(Expr* legends, Expr** bodies, size_t nfun, Expr* single_color);

#endif /* MATHILDA_GRAPHICS_PLOT_COMMON_H */
