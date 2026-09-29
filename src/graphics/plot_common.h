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

/* PhaseRings: domain-colouring ramp for ComplexPlot. Takes the raw complex
 * value (re, im) rather than a normalised t because both Arg and |w| are
 * needed: hue = phase, brightness = 0.1 + 0.9·(1+cos(2π·log|w|))/2 (one
 * bright/dark ring per e-fold of |w|, compressing near poles and zeros).
 * The color-bar path (1-D, t only) falls back to a pure hue sweep. */
void phase_rings_rgb(double re, double im, double* r, double* g, double* b);

/* cyclic_phase_rgb — soft cyclic phase colormap (a desaturated HSV-family
 * rainbow), the default for ComplexPlot's phase colouring and the "Cyclic"
 * named ColorFunction ramp. t ∈ [0,1] is the normalised argument
 * t = (arg+π)/(2π): t=0.5 (arg=0) is red, t=0/1 (arg=±π) is cyan. A 32-stop
 * resampling of the reference bar; cyclic (stop[0] ≈ stop[31]) so the ±π wrap
 * has no colour seam. */
void cyclic_phase_rgb(double t, double* r, double* g, double* b);

/* Thermal: dark blue-purple (t=0) → purple → red → orange → bright yellow (t=1).
 * Matches Mathematica's default StreamPlot speed colormap. */
void thermal_rgb(double t, double* r, double* g, double* b);

/* CoolTones: near-white ice blue (t=0) → sky blue → cornflower → deep navy (t=1). */
void cool_tones_rgb(double t, double* r, double* g, double* b);

/* WarmTones: pale cream (t=0) → amber → orange → deep crimson (t=1). */
void warm_tones_rgb(double t, double* r, double* g, double* b);

/* $RaylibVerbose backing flag. When False (the default), the Raylib backend's
 * TraceLog chatter (window/GL init INFO lines) is suppressed; True lets it
 * through. Set/read via the $RaylibVerbose system variable (see eval.c's
 * sysflag table). Lives here because plot_common is compiled regardless of
 * USE_GRAPHICS and carries no Raylib dependency. */
void raylib_verbose_set(bool on);
bool raylib_verbose_enabled(void);

/* Viridis: perceptually-uniform dark-purple (t=0) → teal → green → yellow (t=1).
 * The Berkeley/matplotlib default; a 32-stop resampling of the 256-entry table. */
void viridis_rgb(double t, double* r, double* g, double* b);

/* The system default magnitude colormap (Viridis). Plotters that map a scalar
 * value/speed/height to colour call this when no ColorFunction is given, so the
 * default is set in one place. */
void default_ramp_rgb(double t, double* r, double* g, double* b);

/* Resolve a ColorFunction name string + t ∈ [0,1] to a color Expr (caller
 * owns).  Recognised names: "Viridis", "Magma", "Plasma", "Inferno", "Cividis",
 * "Haze", "Rainbow", "Temperature"/"Thermal", "CoolTones"/"Cool",
 * "WarmTones"/"Warm", "Greyscale"/"Grayscale"/"Grey"/"Gray".
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

/* ---------------------------------------------------------------------- */
/* Style directives shared by every 2D back end (the Raylib renderer, which
 * draws the on-screen window and the PNG/JPEG export, and the headless PDF
 * writer), so the two resolve Thickness/PointSize/Dashing identically.
 *
 * Sizes follow Mathematica: Thickness[r] and PointSize[d] are fractions of
 * the plot width (d is the point's DIAMETER); AbsoluteThickness[p],
 * AbsolutePointSize[p] and AbsoluteDashing[{...}] are printer's points. The
 * named sizes Tiny/Small/Medium/Large (as in Thick = Thickness[Large] and
 * Dashed = Dashing[{Small, Small}]) resolve to fixed point sizes. */

/* Coerce a numeric Expr to a double: the literal forms of
 * expr_to_real_double first (no evaluation), then N[] for an exact
 * symbolic-but-numeric value (Pi/2, Sqrt[2], 1/3 with bigint parts, ...).
 * Returns false for anything that is not a finite real. */
bool gfx_coerce_double(const Expr* e, double* out);

typedef enum { GFX_SIZE_THICKNESS, GFX_SIZE_POINT, GFX_SIZE_DASH } GfxSizeKind;

/* Tiny/Small/Medium/Large -> printer's points for the given kind. */
bool gfx_named_size(const Expr* e, GfxSizeKind kind, double* pts);

/* Line width in points for Thickness[..] / AbsoluteThickness[..] (plot_w is
 * the plot width in points, or pixels for a raster). False if `d` is not
 * one of those two directives or its argument is unreadable. */
bool gfx_thickness_pts(const Expr* d, double plot_w, double* pts);

/* Point RADIUS in points for PointSize[..] / AbsolutePointSize[..]. */
bool gfx_point_radius_pts(const Expr* d, double plot_w, double* pts);

/* Dash pattern in points for Dashing[..] / AbsoluteDashing[..]: writes up to
 * `max` alternating on/off lengths to out[] and sets *n (0 = solid, as in
 * Dashing[{}] or Dashing[None]). Returns false if `d` is not a dashing
 * directive. */
#define GFX_MAX_DASH 8
bool gfx_dash_pts(const Expr* d, double plot_w, double* out, int max, int* n);

/* PlotStyle resolution for curve/dataset i of a plot. `style` is the
 * evaluated PlotStyle value (NULL, None or Automatic mean "unstyled"); a
 * List of styles is cycled, Mathematica's rule, so for a single curve
 * PlotStyle -> {Red, Thick} uses Red alone. `base` (borrowed) is the colour
 * the curve would otherwise get (the palette entry). Returns an owned
 * directive: a plain colour when the style is just a colour, else
 * Directive[base?, s...] (base is omitted when the style names its own
 * colour). *scoped is set when that directive carries non-colour state
 * (dashing, thickness, point size, opacity) which must be confined to the
 * curve with its own List scope. */
Expr* plot_curve_style(const Expr* style, size_t i, const Expr* base, bool* scoped);

/* The colour a curve style draws in (the first colour inside a Directive),
 * else a copy of `base`. Owned. Used for legend swatches. */
Expr* plot_style_color(const Expr* directive, const Expr* base);

/* build_legend_meta with PlotStyle-aware swatch colours: entry i takes the
 * colour of plot_curve_style(style, i, base_i). style may be NULL. */
Expr* build_legend_meta_styled(Expr* legends, Expr** bodies, size_t nfun,
                               Expr* single_color, const Expr* style);

#endif /* MATHILDA_GRAPHICS_PLOT_COMMON_H */
