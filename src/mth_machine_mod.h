/*
 * mth_machine_mod.h -- the ONE machine-precision Mod[m, n].
 *
 * Every surface that evaluates Mod on two machine reals -- the scalar builtin
 * (mod_fast, src/core.c), the packed/NDArray kernel (NDKB_Mod,
 * src/ndkernels.c) and the Compile[] VM (which lowers a real Mod to that
 * registered kernel) -- goes through mth_mod_machine, so the answers agree
 * bit-for-bit by construction.
 *
 * Why not m - n Floor[m/n]: the quotient m/n is ROUNDED before the floor, so
 * when the exact quotient lies just below an integer the floor overshoots by
 * one and the remainder collapses. Mod[0.5, 1.*^-10] is the example
 * (numeric-tower review G3): the double 1.*^-10 is slightly above 10^-10, so
 * the exact quotient is 4999999999.9999998..., which rounds to exactly 5*^9,
 * and the formula gave 0. where the true remainder is 9.9999981783901346*^-11.
 * The scalar only looked right because the evaluator's degenerate-result net
 * recomputed its 0.; the packed and compiled surfaces returned 0.
 *
 * fmod is EXACT (C99 F.9.7.1: the remainder of two doubles is always
 * representable), so it is the correctly rounded truncated remainder. Mod is
 * the FLOORED remainder (it carries the divisor's sign, Mod[-7., 3.] = 2.), so
 * a nonzero fmod with the wrong sign moves by one n: that one addition is the
 * only rounding, and the result lands in [0, n] (n > 0), with n itself only
 * when |r| is below half an ulp of n, as the floored definition rounds there.
 * A zero remainder is +0., never -0.
 */
#ifndef MTH_MACHINE_MOD_H
#define MTH_MACHINE_MOD_H

#include <math.h>
#include <stdbool.h>

/* Mod[m, n] for machine reals into *out. False (decline: the caller's exact
 * or symbolic path answers) for n == 0 or a non-finite operand or result. */
static inline bool mth_mod_machine(double m, double n, double* out) {
    if (n == 0.0 || !isfinite(m) || !isfinite(n)) return false;
    double r = fmod(m, n);
    if (r != 0.0 && ((r < 0.0) != (n < 0.0))) r += n;
    if (r == 0.0) r = 0.0;                    /* -0. -> +0. */
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

/* Quotient[m, n] = Floor[m/n] for machine reals, as an integer-valued double:
 * the EXACT floor of the exact quotient wherever |q| < 2^52, so that
 * m == n Quotient[m, n] + Mod[m, n] holds with mth_mod_machine's remainder
 * (numeric-tower review, the G3 note: Floor of the ROUNDED quotient made
 * Quotient[0.5, 1.*^-10] 5000000000 where the floor is 4999999999).
 *
 * floor(m/n) is within one of the answer. fma(-q, n, m) is m - q n computed
 * exactly and rounded once, so its sign is exact: k = Floor[m/n] is the q with
 * m - q n >= 0 > m - (q+1) n for n > 0, and the reverse inequalities for n < 0.
 * Past 2^52, q +- 1 is not representable and the rounded floor stands. False
 * (decline) for n == 0 or a non-finite operand or quotient. */
static inline bool mth_quotient_machine(double m, double n, double* out) {
    if (n == 0.0 || !isfinite(m) || !isfinite(n)) return false;
    double q = floor(m / n);
    if (!isfinite(q)) return false;
    if (fabs(q) < 4503599627370496.0) {                 /* 2^52 */
        double lo = fma(-q, n, m);                      /* m - q n, exact sign */
        if (n > 0.0 ? lo < 0.0 : lo > 0.0) {
            q -= 1.0;
        } else {
            double hi = fma(-(q + 1.0), n, m);          /* m - (q + 1) n */
            if (n > 0.0 ? hi >= 0.0 : hi <= 0.0) q += 1.0;
        }
    }
    if (q == 0.0) q = 0.0;                              /* -0. -> +0. */
    *out = q;
    return true;
}

#endif /* MTH_MACHINE_MOD_H */
