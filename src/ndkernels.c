/* ndkernels.c — machine-precision element-wise kernels for the elementary
 * functions, mapped over NDArray buffers at C speed by the evaluator's NDArray
 * fast path (see ndarray_map_unary / ndarray_map_binary and the dispatch hook
 * in eval.c). Each kernel mirrors the corresponding scalar builtin's libc
 * computation (csin, 1/ctan, clog, ...) so an NDArray argument evaluates to the
 * same values the List path would produce, without any per-element Expr.
 *
 * Special functions (Gamma, Erf, Bessel, ...) register their own kernels from
 * their own modules; this file owns only the libc-expressible elementary set. */

#include "ndarray.h"
#include "symtab.h"
#include <complex.h>
#include <math.h>
#include <stdbool.h>

/* ---- descriptor generators ---------------------------------------------- */

/* Real-closed unary: real input -> real output. REXPR (in `x`) is the fast real
 * kernel; CEXPR (in `z`) the complex kernel (used for complex input and float32
 * generic path). */
#define UK_CLOSED(NAME, REXPR, CEXPR)                                          \
    static bool ndk_##NAME##_r(double x, double* o) {                          \
        double v = (REXPR); *o = v; return isfinite(v); }                      \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, ndk_##NAME##_r, true, false }

/* Escaping unary: real input may leave the real axis (Log/ArcSin/...). Complex
 * kernel only; the engine promotes the result dtype iff any element is complex. */
#define UK_ESC(NAME, CEXPR)                                                    \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, false }

/* Escaping unary with a branch-cut fix: on the real (1, inf) cut, C99's casin /
 * cacos / catanh land on the opposite side from the scalar builtins (which
 * follow Mathematica). Flip the imaginary sign for real input x > 1 to agree
 * with builtin_arcsin / builtin_arccos / builtin_arctanh (trig.c:1087,
 * 1130; hyperbolic.c:609). */
#define UK_ESC_FLIP(NAME, CEXPR)                                               \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        if (ai == 0.0 && ar > 1.0) *ri = -*ri;                                 \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, false }

/* Projection: always a real result (Abs/Re/Im/Arg), even for complex input.
 * VEXPR (in `ar`, `ai`) is the real-valued result. */
#define UK_PROJ(NAME, VEXPR)                                                   \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        (void)ar; (void)ai;                                                    \
        double v = (VEXPR); *rr = v; *ri = 0.0; return isfinite(v); }          \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, true }

/* ---- trigonometric (real-closed) ---------------------------------------- */
UK_CLOSED(Sin, sin(x),        csin(z));
UK_CLOSED(Cos, cos(x),        ccos(z));
UK_CLOSED(Tan, tan(x),        ctan(z));
UK_CLOSED(Cot, 1.0 / tan(x),  1.0 / ctan(z));
UK_CLOSED(Sec, 1.0 / cos(x),  1.0 / ccos(z));
UK_CLOSED(Csc, 1.0 / sin(x),  1.0 / csin(z));

/* ---- hyperbolic (real-closed) ------------------------------------------- */
UK_CLOSED(Sinh, sinh(x),       csinh(z));
UK_CLOSED(Cosh, cosh(x),       ccosh(z));
UK_CLOSED(Tanh, tanh(x),       ctanh(z));
UK_CLOSED(Coth, 1.0 / tanh(x), 1.0 / ctanh(z));
UK_CLOSED(Sech, 1.0 / cosh(x), 1.0 / ccosh(z));
UK_CLOSED(Csch, 1.0 / sinh(x), 1.0 / csinh(z));

/* ---- exp / log ---------------------------------------------------------- */
UK_CLOSED(Exp, exp(x), cexp(z));
UK_ESC(Log, clog(z));

/* ---- inverse trig ------------------------------------------------------- */
UK_CLOSED(ArcTan, atan(x),       catan(z));
UK_CLOSED(ArcCot, atan(1.0 / x), catan(1.0 / z));
UK_ESC_FLIP(ArcSin, casin(z));
UK_ESC_FLIP(ArcCos, cacos(z));
UK_ESC(ArcSec, cacos(1.0 / z));
UK_ESC(ArcCsc, casin(1.0 / z));

/* ---- inverse hyperbolic ------------------------------------------------- */
UK_CLOSED(ArcSinh, asinh(x),        casinh(z));
UK_CLOSED(ArcCsch, asinh(1.0 / x),  casinh(1.0 / z));
UK_ESC(ArcCosh, cacosh(z));
UK_ESC_FLIP(ArcTanh, catanh(z));
UK_ESC(ArcCoth, catanh(1.0 / z));
UK_ESC(ArcSech, cacosh(1.0 / z));

/* ---- projections & sign ------------------------------------------------- */
UK_PROJ(Abs, hypot(ar, ai));
UK_PROJ(Re,  ar);
UK_PROJ(Im,  ai);
UK_PROJ(Arg, atan2(ai, ar));
UK_CLOSED(Conjugate, x, conj(z));
/* Sign: real -> {-1,0,1}; complex -> z/|z| (0 at the origin). */
static bool ndk_Sign_r(double x, double* o) {
    *o = (x > 0.0) ? 1.0 : (x < 0.0 ? -1.0 : 0.0); return true;
}
static bool ndk_Sign_c(double ar, double ai, double* rr, double* ri) {
    double m = hypot(ar, ai);
    if (m == 0.0) { *rr = 0.0; *ri = 0.0; return true; }
    *rr = ar / m; *ri = ai / m; return isfinite(*rr) && isfinite(*ri);
}
static const NDUnaryKernel NDKU_Sign = { ndk_Sign_c, ndk_Sign_r, true, false };

/* ---- binary (scalar-index) ---------------------------------------------- */

/* Log[b, z] = clog(z)/clog(b). Kernel receives (a0, a1) = (base, arg). May
 * escape to complex, so real_closed = false (the engine narrows real results). */
static bool ndk_LogB_c(double bre, double bim, double zre, double zim,
                       double* rr, double* ri) {
    double complex b = bre + bim * I, z = zre + zim * I;
    double complex lb = clog(b);
    if (lb == 0.0) return false;
    double complex w = clog(z) / lb;
    *rr = creal(w); *ri = cimag(w);
    return isfinite(*rr) && isfinite(*ri);
}
static const NDBinaryKernel NDKB_Log = { ndk_LogB_c, false };

/* ArcTan[x, y] = arg(x + I y). Real for real inputs (real_closed). */
static bool ndk_ArcTan2_c(double xre, double xim, double yre, double yim,
                          double* rr, double* ri) {
    double complex x = xre + xim * I, y = yre + yim * I;
    double complex s = csqrt(x * x + y * y);
    if (s == 0.0) return false;
    double complex w = -I * clog((x + I * y) / s);
    *rr = creal(w); *ri = cimag(w);
    return isfinite(*rr) && isfinite(*ri);
}
static const NDBinaryKernel NDKB_ArcTan = { ndk_ArcTan2_c, true };

/* ---- special functions: real machine kernels ---------------------------- */
/* These cover the common real-array case at C speed via libc; a complex NDArray
 * (no cplx kernel) or an out-of-fast-domain element (order non-integer, pole)
 * declines and the evaluator degrades that call to the exact List path. */

static bool ndk_Gamma_r(double x, double* o) { double v = tgamma(x); *o = v; return isfinite(v); }
static const NDUnaryKernel NDKU_Gamma = { NULL, ndk_Gamma_r, true, false };

/* LogGamma is real only for x > 0; for x <= 0 it is complex (lgamma loses the
 * imaginary part), so decline and let the List path handle it. */
static bool ndk_LogGamma_r(double x, double* o) {
    if (!(x > 0.0)) return false;
    double v = lgamma(x); *o = v; return isfinite(v);
}
static const NDUnaryKernel NDKU_LogGamma = { NULL, ndk_LogGamma_r, true, false };

static bool ndk_Erf_r(double x, double* o)  { double v = erf(x);  *o = v; return isfinite(v); }
static bool ndk_Erfc_r(double x, double* o) { double v = erfc(x); *o = v; return isfinite(v); }
static const NDUnaryKernel NDKU_Erf  = { NULL, ndk_Erf_r,  true, false };
static const NDUnaryKernel NDKU_Erfc = { NULL, ndk_Erfc_r, true, false };

/* Bessel J/Y of integer order over a real array via libc jn/yn. Kernel receives
 * (order, arg); declines non-integer order or any complex operand -> List path. */
static bool ndk_BesselJ_c(double nre, double nim, double zre, double zim,
                          double* rr, double* ri) {
    if (nim != 0.0 || zim != 0.0 || nre != floor(nre)) return false;
    *rr = jn((int)nre, zre); *ri = 0.0; return isfinite(*rr);
}
static bool ndk_BesselY_c(double nre, double nim, double zre, double zim,
                          double* rr, double* ri) {
    if (nim != 0.0 || zim != 0.0 || nre != floor(nre)) return false;
    *rr = yn((int)nre, zre); *ri = 0.0; return isfinite(*rr);
}
static const NDBinaryKernel NDKB_BesselJ = { ndk_BesselJ_c, true };
static const NDBinaryKernel NDKB_BesselY = { ndk_BesselY_c, true };

/* Beta[a, b] = Gamma[a] Gamma[b] / Gamma[a+b] over real arrays via libc. */
static bool ndk_Beta_c(double are, double aim, double bre, double bim,
                       double* rr, double* ri) {
    if (aim != 0.0 || bim != 0.0) return false;
    double v = tgamma(are) * tgamma(bre) / tgamma(are + bre);
    *rr = v; *ri = 0.0; return isfinite(v);
}
static const NDBinaryKernel NDKB_Beta = { ndk_Beta_c, true };

/* Degrade sentinels: no machine kernel available (libc-free algorithms). The
 * NULL cplx makes the map decline for every element, so the evaluator threads
 * the call over the equivalent nested List (correct results, per-element scalar
 * builtin) — NDArray inputs still evaluate, just not through a C buffer loop. */
static const NDUnaryKernel  ND_DEGRADE_U = { NULL, NULL, false, false };
static const NDBinaryKernel ND_DEGRADE_B = { NULL, false };

/* ---- registration ------------------------------------------------------- */

#define REG_U(NAME) symtab_set_ndarray_unary_kernel(#NAME, &NDKU_##NAME)
#define REG_B(NAME, DESC) symtab_set_ndarray_binary_kernel(#NAME, &DESC)
#define REG_DEG_U(NAME) symtab_set_ndarray_unary_kernel((NAME), &ND_DEGRADE_U)
#define REG_DEG_B(NAME) symtab_set_ndarray_binary_kernel((NAME), &ND_DEGRADE_B)

void ndkernels_init(void) {
    REG_U(Sin);  REG_U(Cos);  REG_U(Tan);
    REG_U(Cot);  REG_U(Sec);  REG_U(Csc);
    REG_U(Sinh); REG_U(Cosh); REG_U(Tanh);
    REG_U(Coth); REG_U(Sech); REG_U(Csch);
    REG_U(Exp);  REG_U(Log);
    REG_U(ArcSin); REG_U(ArcCos); REG_U(ArcTan);
    REG_U(ArcCot); REG_U(ArcSec); REG_U(ArcCsc);
    REG_U(ArcSinh); REG_U(ArcCosh); REG_U(ArcTanh);
    REG_U(ArcCoth); REG_U(ArcSech); REG_U(ArcCsch);
    REG_U(Abs); REG_U(Re); REG_U(Im); REG_U(Arg);
    REG_U(Conjugate); REG_U(Sign);

    REG_B(Log, NDKB_Log);      /* two-arg Log[b, x] */
    REG_B(ArcTan, NDKB_ArcTan); /* two-arg ArcTan[x, y] */

    /* Special functions: real C-loop where libc provides it. */
    REG_U(Gamma); REG_U(LogGamma); REG_U(Erf); REG_U(Erfc);
    REG_B(BesselJ, NDKB_BesselJ);
    REG_B(BesselY, NDKB_BesselY);
    REG_B(Beta,    NDKB_Beta);

    /* Special functions with libc-free algorithms: correct results on NDArray
     * via the degrade path (List threading), pending dedicated machine kernels. */
    static const char* deg_u[] = {
        "Erfi", "ExpIntegralEi", "LogIntegral",
        "SinIntegral", "CosIntegral", "SinhIntegral", "CoshIntegral",
        "FresnelS", "FresnelC", "AiryAi", "AiryBi",
        "AiryAiPrime", "AiryBiPrime", "ProductLog", "Zeta", "Sinc",
    };
    for (unsigned i = 0; i < sizeof(deg_u) / sizeof(deg_u[0]); i++)
        REG_DEG_U(deg_u[i]);
    static const char* deg_b[] = {
        "BesselI", "BesselK", "PolyLog", "HurwitzZeta", "LegendreP",
    };
    for (unsigned i = 0; i < sizeof(deg_b) / sizeof(deg_b[0]); i++)
        REG_DEG_B(deg_b[i]);
}
