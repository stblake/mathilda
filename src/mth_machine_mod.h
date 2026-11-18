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

#endif /* MTH_MACHINE_MOD_H */
