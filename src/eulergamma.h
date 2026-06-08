#ifndef EULERGAMMA_H
#define EULERGAMMA_H

/* Mathilda -- EulerGamma, the Euler-Mascheroni constant gamma.
 *
 *   EulerGamma = lim (H_n - log n)  ~= 0.5772156649015328606...
 *               n->oo
 *
 * EulerGamma is a *constant symbol*, not a builtin function: it carries no
 * C evaluation routine of its own. Its behaviour is supplied by the generic
 * subsystems that already know the symbol:
 *
 *   - N[EulerGamma] / N[EulerGamma, prec] -- machine and arbitrary-precision
 *     numericalisation come from the constant table in src/numeric.c
 *     (the MPFR path uses mpfr_const_euler, so any precision is available).
 *   - D[EulerGamma, x] = 0 -- EulerGamma is in the constant-symbol set used
 *     by the derivative engine (src/calculus/deriv.c).
 *   - NumericQ[EulerGamma] = True -- recognised by is_numeric_quantity
 *     (src/core.c).
 *   - The symbol itself is interned as SYM_EulerGamma (src/sym_names.c).
 *
 * eulergamma_init() exists to make EulerGamma a first-class, REPL-visible
 * symbol: it gives it the Mathematica attributes {Constant, Protected} so
 * that Attributes[EulerGamma] is correct and the constant cannot be
 * accidentally reassigned. The docstring lives in src/info.c. */

void eulergamma_init(void);

#endif /* EULERGAMMA_H */
