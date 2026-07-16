/*
 * integrate_line.h -- Complex line (contour) integration along straight
 * segments and piecewise-linear contours.
 *
 * Given an analytic integrand f, a variable x, and complex points, this module
 * evaluates the contour integral
 *
 *     Integrate[f, {x, z0, z1, ..., zn}]  =  Sum_k Integrate[f, {x, z_k, z_{k+1}}]
 *
 * along the straight segments z_k -> z_{k+1}.  Each segment is parametrised by a
 * real parameter t in [0, 1], gamma(t) = z_k + t (z_{k+1} - z_k), so that
 *
 *   - singularities of f *on the path* become real roots t* in (0, 1) of the
 *     denominator of f(gamma(t)) -- detected exactly as Newton-Leibniz detects
 *     real poles on the real axis; a genuine on-path singularity makes the
 *     contour integral divergent (Integrate::idiv, left unevaluated);
 *
 *   - endpoint limits become real one-sided limits in t (Direction ->
 *     "FromAbove"/"FromBelow"), so the Limit engine can take them even though
 *     the approach is along a complex ray.
 *
 * The antiderivative F = Integrate[f, x] is computed in the friendly variable x
 * (rational/real coefficients).  The segment value is the continuous change of
 * F along the segment.  For rational integrands whose antiderivative carries
 * Log/ArcTan branch cuts, the continuous change is computed branch-correctly by
 * combining logarithms of AFFINE arguments into a single principal Log of a
 * ratio -- exact because a straight segment subtends an angle < pi at any point
 * off it -- so closed-contour residues such as the winding integral of 1/z
 * around the origin come out to 2 Pi I.  Every segment value is numerically
 * cross-checked against a complex quadrature of f(gamma(t)) gamma'(t); a
 * disagreement (an uncorrectable branch crossing) leaves the integral
 * unevaluated rather than return a wrong branch.
 *
 * Reachable three ways:
 *   Integrate[f, {x, a, b}]              (auto-dispatch when a or b is complex)
 *   Integrate[f, {x, z0, ..., zn}]       (piecewise-linear contour)
 *   Integrate`LineIntegral[f, {x, z0, ..., zn}]   (explicit entry point)
 *
 * The on-path singularity detector is exposed for inspection / unit testing:
 *   Integrate`PathSingularPoints[f, {x, z0, ..., zn}]  -> List of the singular
 *                                                         points on the contour.
 */

#ifndef MATHILDA_INTEGRATE_LINE_H
#define MATHILDA_INTEGRATE_LINE_H

#include "expr.h"
#include <stddef.h>
#include <stdbool.h>

/*
 * Contour integral of f in x along the straight segments through pts[0..npts-1]
 * (npts >= 2).  All arguments borrowed (not consumed).  `method` is passed to
 * the inner indefinite Integrate[f, x, Method -> method] call, or NULL for the
 * Automatic cascade.  Returns a freshly-allocated Expr* (the contour value), or
 * NULL to leave the integral unevaluated (unknown antiderivative, on-path
 * singularity / divergence, or an unverifiable branch crossing).
 */
Expr* integrate_line_contour(Expr* f, Expr* x, Expr** pts, size_t npts,
                             const char* method);

/*
 * Contour integral driven from a `{x, z0, z1, ..., zn}` List spec (n >= 1
 * segments).  Borrowed arguments; NULL on malformed spec or unevaluable result.
 */
Expr* integrate_line_from_spec(Expr* f, Expr* spec, const char* method);

/* True iff `e` is a straight-line / contour spec `{x, z0, ..., zn}`: a List
 * with >= 3 elements whose first element is a symbol. */
bool integrate_line_is_contour_spec(const Expr* e);

/* True iff the definite spec `{x, a, b, ...}` should be handled as a complex
 * line integral rather than a real-axis Newton-Leibniz integral: it is a
 * polyline (more than two points) or has at least one non-real endpoint. */
bool integrate_line_spec_is_complex(const Expr* spec);

/* `Integrate`LineIntegral[f, {x, z0, ..., zn}]` builtin. */
Expr* builtin_integrate_line(Expr* res);

/* `Integrate`PathSingularPoints[f, {x, z0, ..., zn}]` builtin. */
Expr* builtin_integrate_path_singular_points(Expr* res);

/* Register the package builtins + attributes + docstrings. */
void integrate_line_init(void);

#endif /* MATHILDA_INTEGRATE_LINE_H */
