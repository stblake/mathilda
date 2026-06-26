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
#include <stdbool.h>

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

#endif /* MATHILDA_GRAPHICS_PLOT_COMMON_H */
