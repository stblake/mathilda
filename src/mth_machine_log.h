/*
 * mth_machine_log.h -- the ONE machine-precision natural logarithm.
 *
 * Every surface that evaluates Log at machine precision -- the scalar
 * builtin (builtin_log, src/logexp.c), the packed/NDArray kernels
 * (NDKU_Log / NDKB_Log, src/ndkernels.c) and the Compile[] VM (OP_LOG_C /
 * OP_VLOG_C, src/compile/compile_vm.c; OP_LOG_R / OP_VLOG_R already call
 * libm log) -- goes through mth_log_machine_z, so the answers agree
 * bit-for-bit by construction on any libm.
 *
 * Why the positive-real arm matters: C99 leaves clog's accuracy to the
 * implementation, and creal(clog(x + 0 I)) is NOT required to equal log(x).
 * Apple's libm returns log(x); glibc does not -- for 1 < |x| < 2 it computes
 * log1p((x-1)(x+1)) / 2 (and log1p(y^2)/2 at |x| = 1), which is accurate but
 * rounds differently from log(x) in the last bit. The scalar builtin and the
 * ND kernel used clog while the compiled VM used log, so Log[1.0016]
 * disagreed by one ulp between the interpreter and Compile[] on Linux only.
 * A positive real therefore always takes libm log; everything else (the
 * negative axis, genuinely complex input) takes clog, whose branch conventions
 * (principal branch, Log[-x] = log x + I Pi) every surface already shares.
 */
#ifndef MTH_MACHINE_LOG_H
#define MTH_MACHINE_LOG_H

#include <math.h>
#include <complex.h>

/* Natural log of the machine complex number z on the principal branch.
 * A positive real z (imaginary part +-0) returns log(creal z) + 0 I. */
static inline double complex mth_log_machine_z(double complex z) {
    double re = creal(z), im = cimag(z);
    if (im == 0.0 && re > 0.0) {
        double complex w = log(re);   /* real -> complex: imaginary part +0 */
        return w;
    }
    return clog(z);
}

/* Log2 / Log10 (and any fixed base b with ln b = LNB): a positive real goes
 * through libm's log2 / log10, so exact powers of the base are exact
 * (Log10[1000.] is 3.); anything else is clog(z) / ln b. The ONE form for the
 * scalar builtins (fixed_base_log, src/logexp.c), the packed / NDArray kernels
 * and the Compile[] real arm: the scalar took a negative argument through
 * Log[2, z] = Log[z] / Log[2], which rounded differently in the last bit of the
 * imaginary part (numeric-tower review G4). */
#define MTH_LN10 2.302585092994045684017991454684364208
#define MTH_LN2  0.693147180559945309417232121458176568
static inline double complex mth_logb_fixed_machine_z(double re, double im,
                                                      double (*rfn)(double), double lnb) {
    if (im == 0.0 && re > 0.0) return rfn(re);
    union { double complex z; double p[2]; } u;   /* keep signed zeros */
    u.p[0] = re; u.p[1] = im;
    return clog(u.z) / lnb;
}

/* Split re/im form of mth_log_machine_z, for the ND kernel signatures. */
static inline void mth_log_machine(double ar, double ai, double* rr, double* ri) {
    double complex w = mth_log_machine_z(ar + ai * I);
    *rr = creal(w);
    *ri = cimag(w);
}

#endif /* MTH_MACHINE_LOG_H */
