/* Mathilda -- arbitrary-precision Mathieu functions (arb / acb kernel).
 *
 * The same algorithm as the machine helpers (sf_machine_mathieu.c), carried out
 * in ball arithmetic: the machine plan supplies the exponent nu, the dominant
 * coefficient index and the characteristic-value decision; nu is then refined
 * by precision-doubling Newton on the continued-fraction consistency function,
 * the Floquet coefficients are glued from two continued fractions, normalized
 * (sum c_k^2 = 1, DLMF for integer nu), and summed -- by real Clenshaw when
 * a, q, z and nu are real, by a centred complex Horner otherwise.
 *
 * Used by flint_num_mathieu (flint_num_bridge.c), which runs the precision
 * ladder around it. Only built with FLINT.
 */
#ifndef MATHILDA_MATHIEU_ARB_H
#define MATHILDA_MATHIEU_ARB_H

#ifdef USE_FLINT
#include <flint/acb.h>

/* MathieuC (deriv == 0) or MathieuCPrime (deriv == 1) at (a, q, z), working
 * precision `prec`. Sets `res` to an indeterminate ball when it cannot answer
 * (the machine plan failed, a degenerate normalization, a non-finite input). */
void mathieu_acb(acb_t res, const acb_t a, const acb_t q, const acb_t z,
                 int deriv, slong prec);
#endif

#endif /* MATHILDA_MATHIEU_ARB_H */
