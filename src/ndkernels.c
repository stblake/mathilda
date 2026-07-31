/* `jn` / `yn` (Bessel functions of integer order) are X/Open extensions to
 * <math.h>, NOT C99. glibc hides them under -std=c99, so the BesselJ / BesselY
 * kernels below fail to compile on Linux with "implicit declaration of function
 * 'jn'" while Darwin — which exposes them regardless of the standard — builds
 * clean, masking the break locally (issue #37). Request the X/Open namespace
 * explicitly, before any include: a no-op on Darwin, the fix on glibc.
 * `make check-c99` enforces this convention tree-wide. */
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

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
#include "gamma.h"
#include "loggamma.h"
/* The complex forms of the exponential-integral family are the interpreter's own
 * ascending series, shared rather than reimplemented — each declines where the
 * series has cancelled away too much of the answer for a double to carry. */
#include "expintegralei.h"
#include "sinintegral.h"
#include "cosintegral.h"
#include "sinhintegral.h"
#include "coshintegral.h"
#include "logintegral.h"
#include "symtab.h"
#include "special_functions/sf_machine.h"
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
        { ndk_##NAME##_c, ndk_##NAME##_r, true, false, NULL, NULL, false }

/* Escaping unary: real input may leave the real axis (Log/ArcSin/...). Complex
 * kernel only; the engine promotes the result dtype iff any element is complex. */
#define UK_ESC(NAME, CEXPR)                                                    \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false }

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
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false }

/* Same fix for the reciprocal inverses (ArcSec/ArcCsc/ArcCoth): here it is
 * 1/x that lands on a cut, so the condition is on the input being real with
 * 0 < |x| < 1 — both cuts, since 1/x is then outside [-1, 1] on either side.
 * Matches builtin_arcsec / builtin_arccsc (trig.c) and builtin_arccoth
 * (hyperbolic.c), which reach the negative side through their Pi-minus / odd
 * folds rather than through this branch. */
#define UK_ESC_RECIP_FLIP(NAME, CEXPR)                                         \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        if (ai == 0.0 && ar != 0.0 && fabs(ar) < 1.0) *ri = -*ri;              \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false }

/* Projection: always a real result (Abs/Re/Im/Arg), even for complex input.
 * VEXPR (in `ar`, `ai`) is the real-valued result. */
#define UK_PROJ(NAME, VEXPR)                                                   \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        (void)ar; (void)ai;                                                    \
        double v = (VEXPR); *rr = v; *ri = 0.0; return isfinite(v); }          \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, true, NULL, NULL, false }

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
UK_ESC_RECIP_FLIP(ArcSec, cacos(1.0 / z));
UK_ESC_RECIP_FLIP(ArcCsc, casin(1.0 / z));

/* ---- inverse hyperbolic ------------------------------------------------- */
UK_CLOSED(ArcSinh, asinh(x),        casinh(z));
UK_CLOSED(ArcCsch, asinh(1.0 / x),  casinh(1.0 / z));
UK_ESC(ArcCosh, cacosh(z));
UK_ESC_FLIP(ArcTanh, catanh(z));
UK_ESC_RECIP_FLIP(ArcCoth, catanh(1.0 / z));
UK_ESC(ArcSech, cacosh(1.0 / z));

/* ---- projections & sign ------------------------------------------------- */
UK_PROJ(Abs, hypot(ar, ai));
/* Abs of an EXACT integer is an exact integer, and hypot() is not that: a
 * projection kernel writes a real dtype, so Abs[{-1, 2}] would be {1., 2.}
 * where the List path gives {1, 2}. The integer arm below makes the buffer
 * answer the same thing, and ndarray_map_unary now falls through to the
 * projection above for every dtype this arm does not cover.
 *
 * |INT64_MIN| is not representable, so that abandons the whole array and the
 * List path answers with the exact bignum -- the same contract as every other
 * integer kernel here. Sequence alignment is what found this: its substitution
 * row is Abs[seq - ch] over an int64 buffer, once per row. */
static bool ndk_Abs_ii(int64_t x, int64_t* o) {
    if (x == INT64_MIN) return false;
    *o = (x < 0) ? -x : x;
    return true;
}
static const NDUnaryKernel NDKU_AbsInt =
    { ndk_Abs_c, NULL, false, true, NULL, ndk_Abs_ii, true };
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
/* Sign's narrowing pair. Real or exact-integer input gives the exact Integer
 * -1/0/1 (Sign[-0.] is 0); COMPLEX input keeps ndk_Sign_c above, which is
 * z/|z| and not an integer at all -- ndarray_map_unary's narrowing branch skips
 * complex dtypes for exactly that reason. Defined here rather than with the
 * other narrowing kernels because this descriptor is declared earlier. */
static bool ndk_Sign_i(double x, int64_t* o) {
    if (!isfinite(x)) return false;
    *o = (x > 0.0) ? 1 : (x < 0.0) ? -1 : 0;        /* -0. gives 0 */
    return true;
}
static bool ndk_Sign_ii(int64_t x, int64_t* o) { *o = (x > 0) - (x < 0); return true; }
static const NDUnaryKernel NDKU_Sign =
    { ndk_Sign_c, ndk_Sign_r, true, false, ndk_Sign_i, ndk_Sign_ii, true };

/* ---- rounding / integer-part (real-closed, real-only) ------------------- */
/* Floor/Ceiling/Round/IntegerPart/FractionalPart over a real array via libc.
 * Real-only (cplx = NULL): a complex NDArray declines and degrades to the List
 * path, where the scalar builtin applies the componentwise complex rule. Round
 * replicates piecewise.c's round_half_even (banker's rounding) exactly rather
 * than trusting the ambient FP rounding mode. */
static bool ndk_Floor_r(double x, double* o)   { double v = floor(x); *o = v; return isfinite(v); }
static bool ndk_Ceiling_r(double x, double* o)  { double v = ceil(x);  *o = v; return isfinite(v); }
static bool ndk_IntegerPart_r(double x, double* o) { double v = trunc(x); *o = v; return isfinite(v); }
static bool ndk_FractionalPart_r(double x, double* o) { double v = x - trunc(x); *o = v; return isfinite(v); }
/* Componentwise, like Floor/Ceiling/Round on a complex: FractionalPart[z] is
 * FractionalPart[Re z] + I FractionalPart[Im z].  Unlike those three it is safe
 * to answer here — they produce Complex[Integer, Integer], which the compile
 * engine's type lattice cannot express, while this one is Complex[Real, Real]. */
static bool ndk_FractionalPart_c(double ar, double ai, double* rr, double* ri) {
    *rr = ar - trunc(ar); *ri = ai - trunc(ai);
    return isfinite(*rr) && isfinite(*ri);
}
static bool ndk_Round_r(double x, double* o) {
    double f = floor(x), r = x - f, v;
    if (r < 0.5)      v = f;
    else if (r > 0.5) v = f + 1.0;
    else              v = (fmod(fabs(f), 2.0) == 0.0) ? f : f + 1.0; /* half -> even */
    *o = v; return isfinite(v);
}
/* ---- narrowing: real in, EXACT INTEGER out ------------------------------ *
 *
 * Floor[{0.25, 0.5, 1.0}] is {0, 0, 1} -- exact Integers, not {0., 0., 1.}. No
 * category above can say that, which is why all of these sat on pack.c's
 * NOT_AWARE list computing nothing, and why UnitStep had no kernel at all and
 * cost ~500 ns/element (Game of Life 658x behind Mathematica, the vectorised
 * Monte Carlo pi 27x).
 *
 * The real forms reuse the very expressions the real-closed kernels above
 * already use, so the two agree on the value by construction and only the output
 * type differs; ndk_narrow_ok is the shared representability guard.
 *
 * Every one of these was checked against the scalar builtin AND Mathematica:
 * Round is banker's (2.5 -> 2, 3.5 -> 4, -2.5 -> -2), IntegerPart truncates
 * toward zero, Sign[-0.] is 0, and UnitStep[0.] and UnitStep[-0.] are both 1
 * (`x >= 0.0` is true for negative zero, which is what makes that fall out). */

/* An exact-integer-valued double that fits int64. Rejects non-finite values and
 * anything outside the range, which abandons the array so the List path answers
 * exactly -- the same abandon-never-wrap contract as ci_*_i64. The bound is
 * written as a power of two because (double)INT64_MAX is not representable and
 * rounds up. */
static bool ndk_narrow_ok(double v, int64_t* out) {
    if (!(v >= -9223372036854775808.0 && v < 9223372036854775808.0)) return false;
    *out = (int64_t)v;
    return true;
}

static bool ndk_Floor_i(double x, int64_t* o)       { return ndk_narrow_ok(floor(x), o); }
static bool ndk_Ceiling_i(double x, int64_t* o)     { return ndk_narrow_ok(ceil(x), o); }
static bool ndk_IntegerPart_i(double x, int64_t* o) { return ndk_narrow_ok(trunc(x), o); }
static bool ndk_Round_i(double x, int64_t* o) {
    double v;
    if (!ndk_Round_r(x, &v)) return false;          /* shares the banker's rule */
    return ndk_narrow_ok(v, o);
}
static bool ndk_UnitStep_i(double x, int64_t* o) {
    if (isnan(x)) return false;
    *o = (x >= 0.0) ? 1 : 0;                        /* -0. >= 0. is true -> 1 */
    return true;
}

/* Exact-integer input. Floor/Ceiling/Round/IntegerPart are the identity on an
 * integer; Sign and UnitStep are not. This arm is what lets the vectorised Game
 * of Life work on the buffer: its neighbour count is an int64 grid, so its
 * UnitStep sees integers, not reals. */
static bool ndk_ident_ii(int64_t x, int64_t* o)    { *o = x; return true; }
static bool ndk_UnitStep_ii(int64_t x, int64_t* o) { *o = (x >= 0) ? 1 : 0; return true; }

/* UnitStep has no real-closed or complex kernel on purpose: a real-typed result
 * would be the wrong head, and there is no complex UnitStep. It is reachable
 * only through the narrowing path. */
static const NDUnaryKernel NDKU_UnitStep =
    { NULL, NULL, false, false, ndk_UnitStep_i, ndk_UnitStep_ii, true };

/* The other five keep their existing real-closed kernels untouched -- the
 * Compile VM and the scalar consumers still use those -- and gain the narrowing
 * pair, which only ndarray_map_unary prefers. */
static const NDUnaryKernel NDKU_Floor          = { NULL, ndk_Floor_r,          true, false, ndk_Floor_i,       ndk_ident_ii, true };
static const NDUnaryKernel NDKU_Ceiling        = { NULL, ndk_Ceiling_r,        true, false, ndk_Ceiling_i,     ndk_ident_ii, true };
static const NDUnaryKernel NDKU_Round          = { NULL, ndk_Round_r,          true, false, ndk_Round_i,       ndk_ident_ii, true };
static const NDUnaryKernel NDKU_IntegerPart    = { NULL, ndk_IntegerPart_r,    true, false, ndk_IntegerPart_i, ndk_ident_ii, true };
static const NDUnaryKernel NDKU_FractionalPart = { ndk_FractionalPart_c, ndk_FractionalPart_r, true, false, NULL, NULL, false };

/* Ramp: the positive part, and the reason it belongs HERE rather than with the
 * narrowing six above is its exactness rule. Floor and UnitStep answer with an
 * exact Integer whatever they are given; Ramp answers at the ARGUMENT's own
 * exactness (Ramp[-1.] is 0., Ramp[-3] is the exact 0), which is precisely what
 * real_closed means. So a Real buffer maps to a Real buffer of the same width
 * and nothing has to be gated.
 *
 * No complex kernel: Ramp[1. + 2. I] is unevaluated in Mathematica, and a NULL
 * cplx makes ndarray_map_unary decline a complex buffer, which degrades to the
 * List path and reproduces that exactly.
 *
 * NaN returns false, abandoning the whole array, because the scalar builtin
 * cannot decide the sign of a NaN either and leaves Ramp[x] unevaluated -- the
 * buffer must not answer where the List would not.
 *
 * No int64 arm: `to_int` is for real-in/integer-out, which is a different
 * function. An int64 buffer therefore materialises, as it does for every head
 * absent from pack.c's INT64_OK list, and the exact Integer answer is unchanged. */
static bool ndk_Ramp_r(double x, double* o) {
    if (isnan(x)) return false;
    *o = (x < 0.0) ? 0.0 : x;
    return true;
}
static const NDUnaryKernel NDKU_Ramp = { NULL, ndk_Ramp_r, true, false, NULL, NULL, false };

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

/* Mod[m, n] = m - n floor(m/n) and Quotient[m, n] = floor(m/n) over real
 * arrays (floored division, result carrying the divisor's sign — matches
 * builtin_mod / builtin_quotient in core.c). Kernel receives (m, n); declines
 * any complex operand or a zero divisor -> List path. real_closed. */
static bool ndk_Mod_c(double mre, double mim, double nre, double nim,
                      double* rr, double* ri) {
    if (mim != 0.0 || nim != 0.0 || nre == 0.0) return false;
    double v = mre - nre * floor(mre / nre);
    *rr = v; *ri = 0.0; return isfinite(v);
}
static bool ndk_Quotient_c(double mre, double mim, double nre, double nim,
                           double* rr, double* ri) {
    if (mim != 0.0 || nim != 0.0 || nre == 0.0) return false;
    double v = floor(mre / nre);
    *rr = v; *ri = 0.0; return isfinite(v);
}
static const NDBinaryKernel NDKB_Mod      = { ndk_Mod_c,      true };
static const NDBinaryKernel NDKB_Quotient = { ndk_Quotient_c, true };

/* ---- special functions: real machine kernels ---------------------------- */
/* These cover the common real-array case at C speed via libc; a complex NDArray
 * (no cplx kernel) or an out-of-fast-domain element (order non-integer, pole)
 * declines and the evaluator degrades that call to the exact List path. */

static bool ndk_Gamma_r(double x, double* o) { double v = tgamma(x); *o = v; return isfinite(v); }
/* Complex Gamma is the interpreter's own Lanczos series, shared rather than
 * reimplemented — see gamma.h.  real_closed stays true, so a real argument
 * still takes tgamma above and only a complex one reaches the Lanczos path. */
static const NDUnaryKernel NDKU_Gamma = { gamma_machine_complex, ndk_Gamma_r, true, false, NULL, NULL, false };

/* Factorial[x] = Gamma[x + 1] over a real array (matches builtin_factorial's
 * real path). A pole (negative integer) yields non-finite -> declines. */
static bool ndk_Factorial_r(double x, double* o) { double v = tgamma(x + 1.0); *o = v; return isfinite(v); }
static const NDUnaryKernel NDKU_Factorial = { NULL, ndk_Factorial_r, true, false, NULL, NULL, false };

/* LogGamma is real only for x > 0; for x <= 0 it is complex (lgamma loses the
 * imaginary part), so decline and let the List path handle it. */
static bool ndk_LogGamma_r(double x, double* o) {
    if (!(x > 0.0)) return false;
    double v = lgamma(x); *o = v; return isfinite(v);
}
/* Complex LogGamma is the CONTINUED log-gamma (imaginary part unbounded in
 * Im z), shared with the interpreter so the branch structure has one home.
 * real_closed stays true: the real kernel above declines for x <= 0, where the
 * true value is complex and a CT_REAL program cannot carry it. */
static const NDUnaryKernel NDKU_LogGamma = { loggamma_machine_complex, ndk_LogGamma_r, true, false, NULL, NULL, false };

static bool ndk_Erf_r(double x, double* o)  { double v = erf(x);  *o = v; return isfinite(v); }
static bool ndk_Erfc_r(double x, double* o) { double v = erfc(x); *o = v; return isfinite(v); }
static const NDUnaryKernel NDKU_Erf  = { NULL, ndk_Erf_r,  true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_Erfc = { NULL, ndk_Erfc_r, true, false, NULL, NULL, false };

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

/* Exponential-integral family and friends.  These modules compute in MPFR and
 * round at the end, so there was no double path to reuse — the kernels in
 * sf_machine.c are written for exactly this.  All are registered
 * real-closed with a `real` function only: the function DECLINES (returns
 * false) for an argument whose true value is complex or infinite, which is the
 * engine's standing contract, and a complex argument has no `cplx` to call so
 * it falls back to the interpreter. */
static bool ndk_ExpIntegralEi_r(double x, double* o) { return sf_machine_ei(x, o); }
static bool ndk_LogIntegral_r(double x, double* o)   { return sf_machine_li(x, o); }
static bool ndk_SinIntegral_r(double x, double* o)   { return sf_machine_si(x, o); }
static bool ndk_CosIntegral_r(double x, double* o)   { return sf_machine_ci(x, o); }
static bool ndk_SinhIntegral_r(double x, double* o)  { return sf_machine_shi(x, o); }
static bool ndk_CoshIntegral_r(double x, double* o)  { return sf_machine_chi(x, o); }
static bool ndk_Sinc_r(double x, double* o)          { return sf_machine_sinc(x, o); }
static bool ndk_InverseErf_r(double x, double* o) {
    if (!(x > -1.0 && x < 1.0)) return false;         /* +-1 are the poles */
    double v = inverf_double(x); *o = v; return isfinite(v);
}
static bool ndk_InverseErfc_r(double x, double* o) {
    if (!(x > 0.0 && x < 2.0)) return false;
    double v = inverfc_double(x); *o = v; return isfinite(v);
}
static const NDUnaryKernel NDKU_ExpIntegralEi = { expintegralei_machine_complex, ndk_ExpIntegralEi_r, true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_LogIntegral   = { logintegral_machine_complex, ndk_LogIntegral_r,   true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_SinIntegral   = { sinintegral_machine_complex, ndk_SinIntegral_r,   true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_CosIntegral   = { cosintegral_machine_complex, ndk_CosIntegral_r,   true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_SinhIntegral  = { sinhintegral_machine_complex, ndk_SinhIntegral_r,  true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_CoshIntegral  = { coshintegral_machine_complex, ndk_CoshIntegral_r,  true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_Sinc          = { NULL, ndk_Sinc_r,          true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_InverseErf    = { NULL, ndk_InverseErf_r,    true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_InverseErfc   = { NULL, ndk_InverseErfc_r,   true, false, NULL, NULL, false };

static bool ndk_Erfi_r(double x, double* o)          { return sf_machine_erfi(x, o); }
static bool ndk_ProductLog_r(double x, double* o)    { return sf_machine_productlog(x, o); }
static bool ndk_FresnelC_r(double x, double* o)      { return sf_machine_fresnel_c(x, o); }
static bool ndk_FresnelS_r(double x, double* o)      { return sf_machine_fresnel_s(x, o); }
static bool ndk_PolyGamma_r(double x, double* o)     { return sf_machine_digamma(x, o); }
static bool ndk_HarmonicNumber_r(double x, double* o){ return sf_machine_harmonic(x, o); }
static bool ndk_Zeta_r(double x, double* o)          { return sf_machine_zeta(x, o); }
static bool ndk_Fibonacci_r(double x, double* o)     { return sf_machine_fibonacci(x, o); }
static bool ndk_LucasL_r(double x, double* o)        { return sf_machine_lucasl(x, o); }
static const NDUnaryKernel NDKU_Erfi           = { NULL, ndk_Erfi_r,           true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_ProductLog     = { NULL, ndk_ProductLog_r,     true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_FresnelC       = { NULL, ndk_FresnelC_r,       true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_FresnelS       = { NULL, ndk_FresnelS_r,       true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_PolyGamma      = { NULL, ndk_PolyGamma_r,      true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_HarmonicNumber = { NULL, ndk_HarmonicNumber_r, true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_Zeta           = { NULL, ndk_Zeta_r,           true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_Fibonacci      = { NULL, ndk_Fibonacci_r,      true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_LucasL         = { NULL, ndk_LucasL_r,         true, false, NULL, NULL, false };
static bool ndk_AiryAi_r(double x, double* o)      { return sf_machine_airy_ai(x, o); }
static bool ndk_AiryBi_r(double x, double* o)      { return sf_machine_airy_bi(x, o); }
static bool ndk_AiryAiPrime_r(double x, double* o) { return sf_machine_airy_ai_prime(x, o); }
static bool ndk_AiryBiPrime_r(double x, double* o) { return sf_machine_airy_bi_prime(x, o); }
static const NDUnaryKernel NDKU_AiryAi      = { NULL, ndk_AiryAi_r,      true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_AiryBi      = { NULL, ndk_AiryBi_r,      true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_AiryAiPrime = { NULL, ndk_AiryAiPrime_r, true, false, NULL, NULL, false };
static const NDUnaryKernel NDKU_AiryBiPrime = { NULL, ndk_AiryBiPrime_r, true, false, NULL, NULL, false };

/* PolyGamma[n, x] and HurwitzZeta[s, a].  PolyGamma matters as a BINARY kernel
 * even when written unary: the evaluator canonicalises PolyGamma[x] to
 * PolyGamma[0, x] before anything downstream sees it. */
static bool ndk_PolyGamma_c(double nre, double nim, double xre, double xim,
                            double* rr, double* ri) {
    if (nim != 0.0 || xim != 0.0) return false;
    double v; if (!sf_machine_polygamma(nre, xre, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_HurwitzZeta_c(double sre, double sim, double are, double aim,
                              double* rr, double* ri) {
    if (sim != 0.0 || aim != 0.0) return false;
    double v; if (!sf_machine_hurwitz_zeta(sre, are, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_Pochhammer_c(double are, double aim, double nre, double nim,
                             double* rr, double* ri) {
    if (aim != 0.0 || nim != 0.0) return false;
    double v; if (!sf_machine_pochhammer(are, nre, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_Binomial_c(double nre, double nim, double kre, double kim,
                           double* rr, double* ri) {
    if (nim != 0.0 || kim != 0.0) return false;
    double v; if (!sf_machine_binomial(nre, kre, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_LegendreP_c(double nre, double nim, double xre, double xim,
                            double* rr, double* ri) {
    if (nim != 0.0 || xim != 0.0) return false;
    double v; if (!sf_machine_legendre_p(nre, xre, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static const NDBinaryKernel NDKB_Pochhammer = { ndk_Pochhammer_c, true };
static const NDBinaryKernel NDKB_Binomial   = { ndk_Binomial_c,   true };
static const NDBinaryKernel NDKB_LegendreP  = { ndk_LegendreP_c,  true };

#define NDK_BIN2(NAME, FN)                                                     \
static bool ndk_##NAME##_c(double are, double aim, double bre, double bim,     \
                           double* rr, double* ri) {                           \
    if (aim != 0.0 || bim != 0.0) return false;                                \
    double v; if (!FN(are, bre, &v)) return false;                             \
    *rr = v; *ri = 0.0; return true;                                           \
}                                                                              \
static const NDBinaryKernel NDKB_##NAME = { ndk_##NAME##_c, true };
NDK_BIN2(BesselI,           sf_machine_bessel_i)
NDK_BIN2(BesselK,           sf_machine_bessel_k)
NDK_BIN2(PolyLog,           sf_machine_polylog)
NDK_BIN2(QPochhammer,       sf_machine_qpochhammer)
NDK_BIN2(Hypergeometric0F1, sf_machine_hyper0f1)
#undef NDK_BIN2

/* N-ary kernels: arity >= 3, which the unary/binary descriptors cannot express.
 * Consumed by the Compile[] VM; the NDArray element-wise path still degrades for
 * these heads, since "one array plus broadcast scalars" is a different question. */
static bool ndk_LerchPhi_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 3 || im[0] || im[1] || im[2]) return false;
    double v; if (!sf_machine_lerchphi(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_Hyper1F1_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 3 || im[0] || im[1] || im[2]) return false;
    double v; if (!sf_machine_hyper1f1(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_Hyper2F1_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 4 || im[0] || im[1] || im[2] || im[3]) return false;
    double v; if (!sf_machine_hyper2f1(re[0], re[1], re[2], re[3], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* PFQ arrives already flattened by the compiler's lowering: re[0] carries p, so
 * q is whatever is left after the p upper parameters and the final z. */
static bool ndk_HypergeometricPFQ_n(const double* re, const double* im, size_t n,
                                    double* rr, double* ri) {
    if (n < 2) return false;
    for (size_t i = 0; i < n; i++) if (im[i] != 0.0) return false;
    size_t p = (size_t)re[0];
    if (n < p + 2) return false;
    size_t q = n - p - 2;
    double v;
    if (!sf_machine_pfq(re + 1, p, re + 1 + p, q, re[n - 1], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static const NDNaryKernel NDKN_HypergeometricPFQ = { ndk_HypergeometricPFQ_n, 0, true };
static const NDNaryKernel NDKN_LerchPhi            = { ndk_LerchPhi_n, 3, true };
static const NDNaryKernel NDKN_Hypergeometric1F1   = { ndk_Hyper1F1_n, 3, true };
static const NDNaryKernel NDKN_Hypergeometric2F1   = { ndk_Hyper2F1_n, 4, true };
#define REG_N(NAME) symtab_set_ndarray_nary_kernel(#NAME, &NDKN_##NAME)
static const NDBinaryKernel NDKB_PolyGamma   = { ndk_PolyGamma_c,   true };
static const NDBinaryKernel NDKB_HurwitzZeta = { ndk_HurwitzZeta_c, true };

/* ---- registration ------------------------------------------------------- */

/* There were DEGRADE SENTINELS here — descriptors with NULL function pointers,
 * registered for heads with no machine numerics so the map declined per element
 * and the evaluator threaded over the equivalent nested List.  Every one of
 * those heads now has a real kernel, so the sentinels, the lists naming them
 * and the loops walking those lists are all gone.  If a head ever again has
 * genuinely unavailable numerics, register `{NULL, NULL, false, false}` for it:
 * an explicit "known, and deliberately has no fast path" reads differently from
 * an absent registry entry, and the coverage audit distinguishes them. */

#define REG_U(NAME) symtab_set_ndarray_unary_kernel(#NAME, &NDKU_##NAME)
#define REG_B(NAME, DESC) symtab_set_ndarray_binary_kernel(#NAME, &DESC)

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
    symtab_set_ndarray_unary_kernel("Abs", &NDKU_AbsInt);
    REG_U(Re); REG_U(Im); REG_U(Arg);
    REG_U(Conjugate); REG_U(Sign);
    REG_U(UnitStep); REG_U(Ramp);
    REG_U(Floor); REG_U(Ceiling); REG_U(Round);
    REG_U(IntegerPart); REG_U(FractionalPart);

    REG_B(Log, NDKB_Log);      /* two-arg Log[b, x] */
    REG_B(ArcTan, NDKB_ArcTan); /* two-arg ArcTan[x, y] */
    REG_B(Mod, NDKB_Mod);       /* Mod[array, n] */
    REG_B(Quotient, NDKB_Quotient); /* Quotient[array, n] */

    /* Special functions: real C-loop where libc provides it. */
    REG_U(Gamma); REG_U(LogGamma); REG_U(Erf); REG_U(Erfc);
    REG_U(Factorial);
    REG_B(BesselJ, NDKB_BesselJ);
    REG_B(BesselY, NDKB_BesselY);
    REG_B(Beta,    NDKB_Beta);

    /* Exponential-integral family (sf_machine.c) + two already-existing
     * double algorithms that were simply never registered. */
    REG_U(ExpIntegralEi); REG_U(LogIntegral);
    REG_U(SinIntegral);   REG_U(CosIntegral);
    REG_U(SinhIntegral);  REG_U(CoshIntegral);
    REG_U(Sinc);          REG_U(InverseErf);   REG_U(InverseErfc);
    REG_U(Erfi);          REG_U(ProductLog);
    REG_U(FresnelC);      REG_U(FresnelS);
    REG_U(PolyGamma);     REG_U(HarmonicNumber);
    REG_U(Zeta);          REG_U(Fibonacci);    REG_U(LucasL);
    REG_B(PolyGamma,   NDKB_PolyGamma);
    REG_B(HurwitzZeta, NDKB_HurwitzZeta);
    REG_B(Pochhammer,  NDKB_Pochhammer);
    REG_B(Binomial,    NDKB_Binomial);
    REG_B(LegendreP,   NDKB_LegendreP);
    REG_U(AiryAi); REG_U(AiryBi); REG_U(AiryAiPrime); REG_U(AiryBiPrime);
    REG_B(BesselI,           NDKB_BesselI);
    REG_B(BesselK,           NDKB_BesselK);
    REG_B(PolyLog,           NDKB_PolyLog);
    REG_B(QPochhammer,       NDKB_QPochhammer);
    REG_B(Hypergeometric0F1, NDKB_Hypergeometric0F1);
    REG_N(LerchPhi); REG_N(Hypergeometric1F1); REG_N(Hypergeometric2F1);
    REG_N(HypergeometricPFQ);

}
