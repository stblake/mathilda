/* Mathilda -- machine-precision Mathieu functions (MathieuC / MathieuCPrime).
 *
 * The even Mathieu function with characteristic value a and parameter q,
 *
 *     MathieuC[a, q, z] = sum_k c_k cos((nu + 2k) z),
 *
 * where nu is the characteristic (Floquet) exponent, cos(pi nu) = y_e(pi) for the
 * even solution y_e (y(0) = 1, y'(0) = 0) of y'' + (a - 2 q cos 2z) y = 0, and the
 * Floquet coefficients solve (a - (nu+2k)^2) c_k = q (c_{k-1} + c_{k+1}).
 *
 * Normalization (the analytic convention; see sf_machine_mathieu.c):
 *   generic nu      sum_k c_k^2 = 1   (plain square: analytic in a and q), the
 *                   sign fixed by Re(c_0) > 0, c_0 the cos(nu z) coefficient;
 *   integer nu      (a is a characteristic value) the DLMF normalization
 *                   (int_0^{2 pi} ce^2 = pi): the generic vector times 1/sqrt 2;
 *   q == 0          cos(sqrt(a) z), except a == 0 where ce_0 = 1/sqrt 2.
 *
 * These helpers are the ONE machine implementation behind every surface: the
 * scalar REPL path, the packed / NDArray buffer loop, and the Compile[] VM all
 * call sf_machine_mathieu, so the four agree bit-for-bit. The (a, q) -> (nu, c_k)
 * plan is cached (keyed by the bit patterns), so a buffer of z at fixed (a, q)
 * pays for the coefficients once and then O(K) flops per point.
 */
#ifndef MATHILDA_SF_MACHINE_MATHIEU_H
#define MATHILDA_SF_MACHINE_MATHIEU_H

#include <stdbool.h>

/* MathieuC (deriv == 0) or MathieuCPrime (deriv == 1) at complex (a, q, z),
 * each given as (re, im) machine doubles. Writes the value to (*rr, *ri).
 * Returns false (declines) on: a non-finite input, a coefficient range too large
 * for the machine plan (|a|, |q| beyond ~10^5), a failed exponent refinement, a
 * degenerate (isotropic) normalization, or |Im z| large enough to overflow. The
 * caller then falls back to the arbitrary-precision path.
 *
 * Real inputs: the exact value is then real or purely imaginary (the ODE is real
 * and the normalization constant squares to a real number); the helper returns
 * it with the other part exactly 0.0. */
bool sf_machine_mathieu(int deriv, double ar, double ai, double qr, double qi,
                        double zr, double zi, double* rr, double* ri);

/* The plan data the arbitrary-precision kernel needs as a starting point:
 * the exponent nu (machine estimate), whether it was snapped to an integer
 * (characteristic value), the dominant coefficient index m, and the index range
 * [lo, hi] that held the machine coefficients. Returns false when the machine
 * plan cannot be built (the arb kernel then declines). q must be nonzero. */
bool sf_machine_mathieu_plan_info(double ar, double ai, double qr, double qi,
                                  double* nur, double* nui, int* integer,
                                  int* m, int* lo, int* hi);

#endif /* MATHILDA_SF_MACHINE_MATHIEU_H */
