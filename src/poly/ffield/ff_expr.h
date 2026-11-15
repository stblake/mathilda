/*
 * ff_expr.h -- Expr <-> elements of K = Q(p_1..p_r)[t]/(mu) and of K[x]
 * (ffield library, src/poly/ffield/; the field itself is ff_fnx.h).
 *
 * Names: x is the polynomial variable, t the generator of the extension (NULL
 * when mu = t, i.e. no extension), p_1..p_r the transcendental parameters.
 * An input expression is a rational function of these symbols with rational
 * coefficients built from Plus, Times, Power (integer exponents), Integer and
 * Rational; anything else (another symbol, a radical, a real) is out of scope
 * and the conversion fails.  A denominator may involve only the parameters
 * (K[x] elements are polynomials in x and t), as for Trager's fieldFn.
 *
 * Outputs are canonical: expanded numerator over the parameters, coefficients
 * of t^j x^d collected, each coefficient a reduced fraction.
 */
#ifndef FF_EXPR_H
#define FF_EXPR_H

#ifdef USE_FLINT
#include "expr.h"
#include "ff_fnx.h"
#include <flint/gr_poly.h>

typedef struct {
    const char* x;            /* interned name of the polynomial variable (may be NULL: scalars only) */
    const char* t;            /* interned name of the extension generator, or NULL */
    const char** p;           /* interned names of the parameters */
    slong np;
} ff_names;

/* e (a polynomial in x and t with coefficients rational in the parameters)
 * -> P in K[x].  Returns 1 on success, 0 out of scope. */
int ff_fnx_poly_from_expr(gr_poly_t P, const Expr* e, const ff_names* nm, gr_ctx_t K);

/* e (free of x) -> c in K.  Returns 1 on success. */
int ff_fnx_elem_from_expr(gr_ptr c, const Expr* e, const ff_names* nm, gr_ctx_t K);

/* P in K[x] -> an expression Sum_d Together[Sum_j c_(d,j) t^j] x^d (unevaluated
 * build, the caller evaluates): the canonical form of Trager's fieldFn
 * (fnCollect), so the package's forms and the kernel's agree identically. */
Expr* ff_fnx_poly_to_expr(const gr_poly_t P, const ff_names* nm, gr_ctx_t K);
Expr* ff_fnx_elem_to_expr(gr_srcptr c, const ff_names* nm, gr_ctx_t K);

#endif /* USE_FLINT */
#endif /* FF_EXPR_H */
