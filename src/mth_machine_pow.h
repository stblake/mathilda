/*
 * mth_machine_pow.h -- the ONE machine-precision Power.
 *
 * Every surface that raises a machine number to a power goes through these
 * helpers, so the answers agree bit-for-bit by construction on any libm:
 *
 *   - the scalar builtin (builtin_power's machine arm, src/power.c);
 *   - the packed / NDArray kernels (cpow_pair and ndarray_scalar_power's
 *     float64 loop, src/ndarray.c);
 *   - the Compile[] VM (OP_POWI_R / OP_POWI_C / OP_POW_R and their tile forms
 *     OP_VPOWI_R / OP_VPOWI_C / OP_VPOW_R, src/compile/compile_vm.c) and the
 *     optimiser's constant fold of OP_POWI_R (src/compile/optimize.c), which
 *     auto-compiled Table / Sum / NIntegrate bodies reuse;
 *   - the numeric loop fast path's OP_POW (src/numloop.c).
 *
 * Before this header the interpreter took libm pow(x, n), the VM binary
 * powering in plain doubles, and the ND kernel x*x*x for n = 3 and pow
 * otherwise, so x^3 disagreed between the interpreter and Compile[] in the
 * last bit for ~26% of random reals, x^5 for ~51% -- and even x^2, because
 * Apple's pow(x, 2) is not x*x for ~0.1% of inputs.
 *
 * Which answer.  Measured against MPFR on 10^5 random x in [-2, 2]:
 *
 *                     max ulp   correctly rounded
 *   libm pow (macOS)   0.51        99.87 %   (every n, including n = 2)
 *   binary powering    1.25 (n=3) 74 %;  2.9 (n=5) 49 %;  11 (n=16) 13 %
 *   double-double      0.50       100.00 %   (every |n| <= 100 tried)
 *
 * so the real integer power is computed in DOUBLE-DOUBLE binary powering
 * (exact products: Dekker's Veltkamp split, or the hardware fma where C99's
 * FP_FAST_FMA says it is fast -- the same bits either way; the makefile's
 * -ffp-contract=off keeps the error-free transformations exact) and rounded
 * once at the end: correctly rounded except in astronomically rare hard
 * cases, and better than libm pow.  It is also cheap where it matters: for a
 * tile of a Compile[] loop or a float64 buffer (mth_powi_machine_array) the
 * element loop is straight-line and vectorises, ~4 ns per element for x^3
 * against ~1.4 ns for the old x*(x*x) and ~30-50 ns for pow; a lone scalar
 * x^3 or x^4 is ~6 ns without fma, ~3 ns with it.  Scope of the arm:
 *
 *   n = 0, 1, 2, -1      1, x, x*x, 1/x: one correctly rounded operation;
 *   n = -2, 3<=|n|<=32   double-double, when |x| lies in [2^-k, 2^k] with
 *                        k = floor(940/|n|), which keeps every intermediate
 *                        and its low part normal and the Veltkamp split from
 *                        overflowing (sf_dd.h's [2^-960, 2^995] contract);
 *                        a negative n takes the double-double reciprocal;
 *   otherwise            libm pow(x, n), which also carries C99 Annex F's
 *                        special values (0, +-inf, NaN, overflow, underflow).
 *
 * A real exponent with an integer value (x^3.) takes the same path as the
 * Integer 3 (mth_pow_machine): the packed kernel and the VM's real-exponent
 * opcode cannot tell the two apart, so they must answer alike.  x^0.5 of a
 * non-negative x is sqrt(x) (correctly rounded, which pow is not).
 *
 * A complex base with an integer exponent (mth_powi_machine_c) is binary
 * powering in plain doubles -- exact for Gaussian-integer powers (I^2 is
 * -1 + 0 I, not cpow's -1 + 1.2e-16 I) -- and a negative exponent inverts
 * the power at the end by Smith's division, which neither overflows nor
 * underflows |z|^2.  Every surface takes this one routine; before, the
 * interpreter used cpow, the VM inverted the base first through C's complex
 * division, and the ND kernel inverted at the end through conj(z)/|z|^2.
 */
#ifndef MTH_MACHINE_POW_H
#define MTH_MACHINE_POW_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

/* Largest |n| taken by the double-double arm. Beyond it pow costs no more and
 * the safe |x| range has shrunk to [2^-29, 2^29]. */
#define MTH_POWI_DD_MAX 32

/* MTH_POWI_BOUND[m] = 2^floor(940/m): the double-double arm needs
 * 1/bound <= |x| <= bound, so |x|^m stays inside [2^-940, 2^940]. */
static const double MTH_POWI_BOUND[MTH_POWI_DD_MAX + 1] = {
    0.0, 0.0, 0x1p470, 0x1p313, 0x1p235, 0x1p188, 0x1p156, 0x1p134, 0x1p117,
    0x1p104, 0x1p94, 0x1p85, 0x1p78, 0x1p72, 0x1p67, 0x1p62, 0x1p58, 0x1p55,
    0x1p52, 0x1p49, 0x1p47, 0x1p44, 0x1p42, 0x1p40, 0x1p39, 0x1p37, 0x1p36,
    0x1p34, 0x1p33, 0x1p32, 0x1p31, 0x1p30, 0x1p29
};

/* The exact product a*b = h + l. The error of a rounded product is a unique
 * number, so the hardware fused multiply-add (where C99's FP_FAST_FMA says it
 * is fast: arm64, an x86 build with FMA enabled) and Dekker's Veltkamp-split
 * algorithm return the SAME l inside the range the callers gate on -- the
 * choice changes speed, never a bit. */
static inline void mth_pow_two_prod(double a, double b, double* h, double* l) {
#ifdef FP_FAST_FMA
    *h = a * b;
    *l = fma(a, b, -*h);
#else
    const double SPLIT = 134217729.0;                      /* 2^27 + 1 */
    double t = SPLIT * a, ah = t - (t - a), al = a - ah;
    t = SPLIT * b;
    double bh = t - (t - b), bl = b - bh;
    *h = a * b;
    *l = ((ah * bh - *h) + ah * bl + al * bh) + al * bl;
#endif
}

/* The pairs below are NOT renormalised between steps: hi is the rounded
 * product, lo its exact error plus the cross terms, so |lo| stays within about
 * one ulp of hi and the next exact product still captures ~104 bits. Skipping
 * the two-sum renormalisation shortens the dependency chain by a third (x^5
 * measured 27 -> 11 ns per scalar call without fma) and costs nothing in the
 * final rounding: the tests check the result is correctly rounded. */

/* (rh, rl) <- (rh, rl) * (bh, bl). */
static inline void mth_pow_dd_mul(double* rh, double* rl, double bh, double bl) {
    double h, l;
    mth_pow_two_prod(*rh, bh, &h, &l);
    *rl = l + (*rh * bl + *rl * bh);
    *rh = h;
}

/* (bh, bl) <- (bh, bl)^2. */
static inline void mth_pow_dd_sqr(double* bh, double* bl) {
    double h, l;
    mth_pow_two_prod(*bh, *bh, &h, &l);
    *bl = l + 2.0 * *bh * *bl;
    *bh = h;
}

/* Round the double-double r = rh + rl, or its reciprocal when `neg`. */
static inline double mth_pow_dd_finish(double rh, double rl, int neg) {
    if (!neg) return rh + rl;
    double q = 1.0 / rh, h, l;
    mth_pow_two_prod(q, rh, &h, &l);
    double e = ((1.0 - h) - l) - q * rl;                  /* 1 - q r */
    return q + q * e;
}

/* x^3 and x^4 in double-double, straight-line: the exact operations
 * mth_powi_dd's loop performs for m = 3, 4 (its extra terms are products with
 * a zero low part), without the loop -- the shape a scalar Compile[] loop body
 * pays for. */
static inline double mth_pow_dd_cube(double x) {
    double h, l, p, q;
    mth_pow_two_prod(x, x, &h, &l);
    mth_pow_two_prod(x, h, &p, &q);
    return p + (q + x * l);
}
static inline double mth_pow_dd_quart(double x) {
    double h, l, p, q;
    mth_pow_two_prod(x, x, &h, &l);
    mth_pow_two_prod(h, h, &p, &q);
    return p + (q + 2.0 * h * l);
}

/* x^m (or x^-m when neg) in double-double, 2 <= m <= MTH_POWI_DD_MAX and x
 * inside the MTH_POWI_BOUND range. mth_powi_machine_array performs exactly
 * this operation sequence per element. */
static inline double mth_powi_dd(double x, unsigned m, int neg) {
    if (!neg && m == 3) return mth_pow_dd_cube(x);
    if (!neg && m == 4) return mth_pow_dd_quart(x);
    double rh = 1.0, rl = 0.0, bh = x, bl = 0.0;
    int first = 1;
    while (m) {
        if (m & 1u) {
            if (first) { rh = bh; rl = bl; first = 0; }
            else mth_pow_dd_mul(&rh, &rl, bh, bl);
        }
        m >>= 1;
        if (m) mth_pow_dd_sqr(&bh, &bl);
    }
    return mth_pow_dd_finish(rh, rl, neg);
}

/* |n| when 3 <= n <= MTH_POWI_DD_MAX or -MTH_POWI_DD_MAX <= n <= -2, else 0. */
static inline unsigned mth_powi_dd_order(int64_t n) {
    if (n >= 3 && n <= MTH_POWI_DD_MAX) return (unsigned)n;
    if (n <= -2 && n >= -MTH_POWI_DD_MAX) return (unsigned)(-n);
    return 0;
}

/* x^n for a machine real x and an integer n. */
static inline double mth_powi_machine(double x, int64_t n) {
    switch (n) {
        case 0:  return 1.0;
        case 1:  return x;
        case 2:  return x * x;
        case -1: return 1.0 / x;
        default: break;
    }
    unsigned m = mth_powi_dd_order(n);
    if (m) {
        double ax = fabs(x), b = MTH_POWI_BOUND[m];
        if (ax <= b && ax * b >= 1.0) return mth_powi_dd(x, m, n < 0);
    }
    return pow(x, (double)n);
}

/* x^y for machine reals whose result is real (x >= 0, or y an integer). */
static inline double mth_pow_machine(double x, double y) {
    if (y >= -MTH_POWI_DD_MAX && y <= MTH_POWI_DD_MAX && y == (double)(int)y)
        return mth_powi_machine(x, (int64_t)(int)y);
    if (y == 0.5 && x >= 0.0) return sqrt(x);
    return pow(x, y);
}

/* out[k] = mth_powi_machine(in[k], n), bit for bit, for a whole buffer. When
 * the double-double arm applies to every element the loops are straight-line
 * over a block, so they vectorise; in == out is allowed. */
static inline void mth_powi_machine_array(const double* in, double* out, size_t len, int64_t n) {
    size_t k;
    switch (n) {
        case 0:  for (k = 0; k < len; k++) out[k] = 1.0;                    return;
        case 1:  for (k = 0; k < len; k++) out[k] = in[k];                  return;
        case 2:  for (k = 0; k < len; k++) { double t = in[k]; out[k] = t * t; } return;
        case -1: for (k = 0; k < len; k++) out[k] = 1.0 / in[k];            return;
        default: break;
    }
    unsigned m0 = mth_powi_dd_order(n);
    if (!m0) { for (k = 0; k < len; k++) out[k] = pow(in[k], (double)n); return; }
    const double b = MTH_POWI_BOUND[m0];
    const int neg = n < 0;
    enum { BLK = 64 };
    double rh[BLK], rl[BLK], bh[BLK], bl[BLK];
    for (size_t base = 0; base < len; base += BLK) {
        size_t cnt = len - base < BLK ? len - base : BLK;
        const double* src = in + base;
        int inside = 1;
        for (k = 0; k < cnt; k++) {
            double ax = fabs(src[k]);
            inside &= (ax <= b) & (ax * b >= 1.0);
        }
        if (!inside) {                      /* some element needs pow */
            for (k = 0; k < cnt; k++) out[base + k] = mth_powi_machine(src[k], n);
            continue;
        }
        for (k = 0; k < cnt; k++) { bh[k] = src[k]; bl[k] = 0.0; }
        unsigned m = m0;
        int first = 1;
        while (m) {
            if (m & 1u) {
                if (first) { for (k = 0; k < cnt; k++) { rh[k] = bh[k]; rl[k] = bl[k]; } first = 0; }
                else for (k = 0; k < cnt; k++) mth_pow_dd_mul(&rh[k], &rl[k], bh[k], bl[k]);
            }
            m >>= 1;
            if (m) for (k = 0; k < cnt; k++) mth_pow_dd_sqr(&bh[k], &bl[k]);
        }
        for (k = 0; k < cnt; k++) out[base + k] = mth_pow_dd_finish(rh[k], rl[k], neg);
    }
}

/* (ar + ai I)^n for a machine complex base and an integer n: binary powering
 * in plain doubles, a negative n inverting the power by Smith's division. */
static inline void mth_powi_machine_c(double ar, double ai, int64_t n,
                                      double* outr, double* outi) {
    uint64_t m = (n < 0) ? (uint64_t)0 - (uint64_t)n : (uint64_t)n;
    double rr = 1.0, ri = 0.0, br = ar, bi = ai;
    int first = 1;
    while (m) {
        if (m & 1u) {
            if (first) { rr = br; ri = bi; first = 0; }
            else {
                double nr = rr * br - ri * bi;
                double ni = rr * bi + ri * br;
                rr = nr; ri = ni;
            }
        }
        m >>= 1;
        if (m) {
            double sr = br * br - bi * bi;
            double si = 2.0 * br * bi;
            br = sr; bi = si;
        }
    }
    if (n < 0) {                            /* 1 / (rr + ri I), Smith */
        if (fabs(rr) >= fabs(ri)) {
            double t = ri / rr, d = rr + ri * t;
            *outr = 1.0 / d; *outi = -t / d;
        } else {
            double t = rr / ri, d = rr * t + ri;
            *outr = t / d;  *outi = -1.0 / d;
        }
        return;
    }
    *outr = rr; *outi = ri;
}

#endif /* MTH_MACHINE_POW_H */
