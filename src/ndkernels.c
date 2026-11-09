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
#include "elliptic.h"
#include "symtab.h"
#include "mth_machine_log.h"                /* the one machine Log (scalar/ND/VM) */
#include "piecewise.h"                         /* Unitize machine helpers */
#include "special_functions/sf_machine.h"
#include "special_functions/sf_machine_gamma.h"
#include "special_functions/sf_machine_expint.h"
#include "special_functions/sf_machine_ellk.h"
#include "special_functions/sf_machine_plog.h"
#include "special_functions/sf_machine_zeta.h" /* zeta-family real + complex kernels */
#include "special_functions/sf_machine_hyper.h"  /* complex hypergeometric-family kernels */
#include "special_functions/sf_machine_error.h"   /* complex erf family, Fresnel, Sinc, InverseErf(c) */
#include "special_functions/sf_machine_bessel.h"
#include "special_functions/parabolic_cylinder.h" /* sf_machine_pcfd[_c] */
#include "special_functions/gegenbauer.h"         /* sf_machine_gegenbauer[2]_c */
#include "special_functions/sf_machine_meijerg.h" /* MeijerG flattened n-ary kernel */
#include "special_functions/sf_machine_mathieu.h" /* MathieuC / MathieuCPrime */
#include "special_functions/sf_machine_dirichlet.h" /* DirichletL n-ary kernel */
#include "special_functions/surd.h"          /* Surd binary kernel (sf_surd_r), CubeRoot */
#include "special_functions/riemannr.h"      /* RiemannR unary kernel */
#include "special_functions/besseljzero.h"   /* BesselJZero binary kernel */
#include "special_functions/mathieu_char.h"  /* MathieuCharacteristicA/B kernels */
#include "special_functions/riemann_siegel.h" /* RiemannSiegelZ / Theta unary kernels */
#include "special_functions/sf_machine_nielsen.h" /* Nielsen PolyLog n-ary kernel */
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
        { ndk_##NAME##_c, ndk_##NAME##_r, true, false, NULL, NULL, false, NULL }

/* Escaping unary: real input may leave the real axis (Log/ArcSin/...). Complex
 * kernel only; the engine promotes the result dtype iff any element is complex. */
#define UK_ESC(NAME, CEXPR)                                                    \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        double complex z = ar + ai * I; double complex w = (CEXPR);            \
        *rr = creal(w); *ri = cimag(w);                                        \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false, NULL }

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
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false, NULL }

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
        { ndk_##NAME##_c, NULL, false, false, NULL, NULL, false, NULL }

/* Projection: always a real result (Abs/Re/Im/Arg), even for complex input.
 * VEXPR (in `ar`, `ai`) is the real-valued result. */
#define UK_PROJ(NAME, VEXPR)                                                   \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        (void)ar; (void)ai;                                                    \
        double v = (VEXPR); *rr = v; *ri = 0.0; return isfinite(v); }          \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, NULL, false, true, NULL, NULL, false, NULL }

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
UK_ESC(Log, mth_log_machine_z(z));   /* positive real -> libm log, as the VM */

/* Log10 / Log2: escaping like Log (a negative real lands off the real axis),
 * but a positive real goes through libm's log10 / log2 rather than
 * clog(z)/ln(b), so exact powers of the base are exact -- Log10[1000.] is 3.,
 * where log(1000.)/log(10.) is 2.9999999999999996 -- and every element agrees
 * with builtin_log10 / builtin_log2, which take the same libm path. The real
 * arm is the Compile[] scalar path; it declines off the positive axis so the
 * VM bails to the interpreter instead of inventing a real answer. */
#define ND_LN10 2.302585092994045684017991454684364208
#define ND_LN2  0.693147180559945309417232121458176568
#define UK_FIXED_LOG(NAME, RFN, LNB)                                            \
    static bool ndk_##NAME##_c(double ar, double ai, double* rr, double* ri) { \
        if (ai == 0.0 && ar > 0.0) {                                           \
            *rr = RFN(ar); *ri = 0.0; return isfinite(*rr); }                  \
        double complex w = clog(ar + ai * I) / (LNB);                          \
        *rr = creal(w); *ri = cimag(w);                                        \
        return isfinite(*rr) && isfinite(*ri); }                               \
    static bool ndk_##NAME##_r(double x, double* o) {                          \
        if (!(x > 0.0)) return false;                                          \
        *o = RFN(x); return isfinite(*o); }                                    \
    static const NDUnaryKernel NDKU_##NAME =                                   \
        { ndk_##NAME##_c, ndk_##NAME##_r, false, false, NULL, NULL, false, NULL }
UK_FIXED_LOG(Log10, log10, ND_LN10);
UK_FIXED_LOG(Log2,  log2,  ND_LN2);

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
/* Abs: written by hand (not UK_PROJ) so it leaves no unused projection
 * descriptor -- the registered kernel is NDKU_AbsInt below, which reuses this
 * ndk_Abs_c for the real/complex magnitude and adds the integer arm. */
static bool ndk_Abs_c(double ar, double ai, double* rr, double* ri) {
    (void)ar; (void)ai; double v = hypot(ar, ai); *rr = v; *ri = 0.0; return isfinite(v); }
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
    { ndk_Abs_c, NULL, false, true, NULL, ndk_Abs_ii, true, NULL };
UK_PROJ(Re,  ar);
/* Im: the imaginary part. Written by hand rather than via UK_PROJ because it
 * needs a NARROWING pair and must not leave an unused projection descriptor.
 * Real -- or exact-integer -- input is the exact Integer 0, so the narrowing arm
 * writes an NDT_INT64 zero buffer to match the List path (Im left pack.c's
 * NOT_AWARE list on 2026-08-12, where it had cost ~22x an aware head). It is
 * unconditionally 0: the imaginary part of any real value is 0, finite or not,
 * so unlike Sign there is nothing to reject. COMPLEX input takes ndk_Im_c (the
 * imaginary part, a Real), which ndarray_map_unary's narrowing branch skips for
 * complex dtypes -- exactly as Sign and Abs do. */
static bool ndk_Im_c(double ar, double ai, double* rr, double* ri) {
    (void)ar; *rr = ai; *ri = 0.0; return isfinite(ai); }
static bool ndk_Im_ir(double x, int64_t* o) { (void)x; *o = 0; return true; }
static bool ndk_Im_ii(int64_t x, int64_t* o) { (void)x; *o = 0; return true; }
static const NDUnaryKernel NDKU_ImInt =
    { ndk_Im_c, NULL, false, true, ndk_Im_ir, ndk_Im_ii, true, NULL };
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
    { ndk_Sign_c, ndk_Sign_r, true, false, ndk_Sign_i, ndk_Sign_ii, true, NULL };

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
static bool ndk_UnitBox_i(double x, int64_t* o) {
    if (isnan(x)) return false;
    *o = (x >= -0.5 && x <= 0.5) ? 1 : 0;           /* closed box: UnitBox[1/2] = 1 */
    return true;
}

/* Exact-integer input. Floor/Ceiling/Round/IntegerPart are the identity on an
 * integer; Sign and UnitStep are not. This arm is what lets the vectorised Game
 * of Life work on the buffer: its neighbour count is an int64 grid, so its
 * UnitStep sees integers, not reals. */
static bool ndk_ident_ii(int64_t x, int64_t* o)    { *o = x; return true; }
static bool ndk_UnitStep_ii(int64_t x, int64_t* o) { *o = (x >= 0) ? 1 : 0; return true; }
/* An integer lies in the closed box [-1/2, 1/2] iff it is exactly 0. */
static bool ndk_UnitBox_ii(int64_t x, int64_t* o)  { *o = (x == 0) ? 1 : 0; return true; }

/* UnitStep and UnitBox have no real-closed or complex kernel on purpose: a
 * real-typed result would be the wrong head (both answer with an exact Integer
 * 0/1), and there is no complex form of either. They are reachable only through
 * the narrowing path. */
static const NDUnaryKernel NDKU_UnitStep =
    { NULL, NULL, false, false, ndk_UnitStep_i, ndk_UnitStep_ii, true, NULL };
static const NDUnaryKernel NDKU_UnitBox =
    { NULL, NULL, false, false, ndk_UnitBox_i, ndk_UnitBox_ii, true, NULL };

/* Unitize: narrowing-only for the same reason as UnitStep -- the answer is the
 * exact Integer 0/1 whatever the input precision (Unitize[{1., 0.}] is
 * {1, 0}). The arms ARE the scalar builtin's machine helpers (piecewise.c), so
 * the scalar, List, packed, NDArray and SparseArray surfaces share one
 * implementation. No complex arm: a complex buffer declines and the List path
 * answers (Unitize[1 + I] is 1). NaN declines the whole array. */
static const NDUnaryKernel NDKU_Unitize =
    { NULL, NULL, false, false, unitize_machine_r, unitize_machine_i, true, NULL };
/* Unitize[x, dx]: |x| >= dx -> 1, both exact-integer arms. A threshold that is
 * not a finite positive number makes the helper decline, so the List path
 * issues Unitize::post per element exactly as it would unpacked. */
static const NDBinaryKernel NDKB_Unitize =
    { NULL, false, unitize_machine_tol_r, unitize_machine_tol_i, true, NULL, NULL };

/* The other five keep their existing real-closed kernels untouched -- the
 * Compile VM and the scalar consumers still use those -- and gain the narrowing
 * pair, which only ndarray_map_unary prefers. */
static const NDUnaryKernel NDKU_Floor          = { NULL, ndk_Floor_r,          true, false, ndk_Floor_i,       ndk_ident_ii, true, NULL };
static const NDUnaryKernel NDKU_Ceiling        = { NULL, ndk_Ceiling_r,        true, false, ndk_Ceiling_i,     ndk_ident_ii, true, NULL };
static const NDUnaryKernel NDKU_Round          = { NULL, ndk_Round_r,          true, false, ndk_Round_i,       ndk_ident_ii, true, NULL };
static const NDUnaryKernel NDKU_IntegerPart    = { NULL, ndk_IntegerPart_r,    true, false, ndk_IntegerPart_i, ndk_ident_ii, true, NULL };
static const NDUnaryKernel NDKU_FractionalPart = { ndk_FractionalPart_c, ndk_FractionalPart_r, true, false, NULL, NULL, false, NULL };

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
static const NDUnaryKernel NDKU_Ramp = { NULL, ndk_Ramp_r, true, false, NULL, NULL, false, NULL };

/* ---- binary (scalar-index) ---------------------------------------------- */

/* Log[b, z] = Log(z)/Log(b) through mth_log_machine_z, so a positive real
 * pair is log(z)/log(b) -- what builtin_log and Compile[] compute. Kernel receives (a0, a1) = (base, arg). May
 * escape to complex, so real_closed = false (the engine narrows real results). */
static bool ndk_LogB_c(double bre, double bim, double zre, double zim,
                       double* rr, double* ri) {
    double complex b = bre + bim * I, z = zre + zim * I;
    double complex lb = mth_log_machine_z(b);
    if (lb == 0.0) return false;
    double complex w = mth_log_machine_z(z) / lb;
    *rr = creal(w); *ri = cimag(w);
    return isfinite(*rr) && isfinite(*ri);
}
static const NDBinaryKernel NDKB_Log = { ndk_LogB_c, false, NULL, NULL, false, NULL, NULL };

/* ArcTan[x, y] = arg(x + I y). Real for real inputs (real_closed). */
static bool ndk_ArcTan2_c(double xre, double xim, double yre, double yim,
                          double* rr, double* ri) {
    /* REAL inputs go straight to atan2, which is what the scalar builtin calls
     * (src/trig.c). The complex form below is arg via csqrt and clog, and on
     * real data it lands 1-2 ulp away from atan2 -- 68 of 400 elements differed
     * on the sweep's `arctan2` probe. Under the sweep's 1e-5 tolerance that
     * never showed, but the contract this file states at the top is that a
     * kernel reproduces the scalar builtin's libc computation, and here it
     * simply had not. atan2 also gets the branch at (0, 0), the negative-zero
     * sign convention and the infinities right by definition, where the
     * csqrt/clog route declines on s == 0 and abandoned the whole array for one
     * such pair. */
    if (xim == 0.0 && yim == 0.0) {
        *rr = atan2(yre, xre); *ri = 0.0;
        return isfinite(*rr);
    }
    double complex x = xre + xim * I, y = yre + yim * I;
    double complex s = csqrt(x * x + y * y);
    if (s == 0.0) return false;
    double complex w = -I * clog((x + I * y) / s);
    *rr = creal(w); *ri = cimag(w);
    return isfinite(*rr) && isfinite(*ri);
}
static const NDBinaryKernel NDKB_ArcTan = { ndk_ArcTan2_c, true, NULL, NULL, false, NULL, NULL };

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
/* ---- Mod / Quotient: the EXACT integer arms ------------------------------ *
 *
 * WHY THESE EXIST. `Mod` had a kernel and was packed_aware, but was absent from
 * pack.c's INT64_OK -- correctly, because the double kernel above would write
 * {0., 1.} where the List path gives the exact {0, 1}. The consequence was not
 * that Mod was slow. It was that Mod DESTROYED PACKING, and everything
 * downstream of it paid:
 *
 *     jv = Mod[Range[10^6] 7919, 1000]      came back a plain List of 10^6
 *                                           boxed Integers
 *     Union[kv]        807 ms   ->  4.58 ms once the buffer survives  (176x)
 *     Tally[jv]        68.1 ms  ->  1.89 ms                            (36x)
 *     DeleteDuplicates 67.1 ms  ->  1.80 ms                            (37x)
 *
 * The set operations were never the problem: they are on AWARE and INT64_OK and
 * are fast when handed a buffer. They never were, because the very first
 * operation in the pipeline dropped it. That is the shape of finding this
 * round's audit was written to see -- a head whose own row looks fine and whose
 * cost lands on its consumers.
 *
 * SEMANTICS, matched to builtin_mod / builtin_quotient in core.c rather than to
 * C's operators, which truncate toward zero where these floor:
 *   Quotient[m, n] = Floor[m/n]           (core.c: mpz_fdiv_q)
 *   Mod[m, n]      = m - n Floor[m/n]     so the result carries n's sign
 * Checked against Mathematica: Mod[-7, 3] = 2, Mod[7, -3] = -2,
 * Quotient[-7, 3] = -3, Quotient[7, -3] = -3.
 *
 * Every "return false" abandons the whole array to the List path, which answers
 * with GMP:
 *   n == 0        Mod[m, 0] and Quotient[m, 0] are unevaluated, not a trap
 *   n == -1       m % -1 and m / -1 are UB at m == INT64_MIN. Mod is 0 for all
 *                 m, so that one is answered directly; Quotient[INT64_MIN, -1]
 *                 overflows int64 and declines.
 *   real, unrepresentable  ndk_narrow_ok's range test, as for Floor. */
static bool ndk_Mod_ii(int64_t m, int64_t n, int64_t* o) {
    if (n == 0) return false;
    if (n == -1 || n == 1) { *o = 0; return true; }   /* also dodges the UB */
    int64_t r = m % n;
    if (r != 0 && ((r < 0) != (n < 0))) r += n;       /* truncated -> floored */
    *o = r;
    return true;
}

/* Mod[arr, n] over a whole int64 buffer for a constant n >= 2: the quotient is
 * estimated from a precomputed double reciprocal and corrected exactly in
 * integer arithmetic, so the loop has no hardware division and no per-element
 * call -- it vectorises. For |x| < 2^52 the estimate floor(x * (1/n)) is off by
 * at most one, which the two conditional corrections absorb; a larger element
 * takes the % path. Bit-identical to ndk_Mod_ii (floored remainder in [0, n)).
 * Declines (per-element loop) for n < 2 and for Mod[n, arr]. */
static bool ndk_Mod_whole(const int64_t* in, size_t sz, int64_t n, bool arr_first,
                          int64_t* out) {
    if (!arr_first || n < 2) return false;
    const double inv = 1.0 / (double)n;
    const int64_t lim = (int64_t)1 << 52;
    for (size_t i = 0; i < sz; i++) {
        int64_t x = in[i];
        int64_t r;
        if (x > -lim && x < lim) {
            int64_t q = (int64_t)floor((double)x * inv);
            r = x - q * n;
            r += (r < 0) ? n : 0;
            r -= (r >= n) ? n : 0;
        } else {
            r = x % n;
            if (r < 0) r += n;
        }
        out[i] = r;
    }
    return true;
}

static bool ndk_Quotient_ii(int64_t m, int64_t n, int64_t* o) {
    if (n == 0) return false;
    if (n == -1 && m == INT64_MIN) return false;      /* overflows int64 */
    int64_t q = m / n, r = m % n;
    if (r != 0 && ((r < 0) != (n < 0))) q -= 1;       /* truncated -> floored */
    *o = q;
    return true;
}

/* Quotient NARROWS: real in, exact Integer out, matching builtin_quotient's
 * real fallback, which ends in expr_new_integer((int64_t)floor(...)). Mod has
 * no such arm on purpose -- its result takes the argument's exactness, and
 * `real_closed` already says that. */
static bool ndk_Quotient_ri(double m, double n, int64_t* o) {
    if (n == 0.0) return false;
    return ndk_narrow_ok(floor(m / n), o);
}

static const NDBinaryKernel NDKB_Mod =
    { ndk_Mod_c,      true, NULL,             ndk_Mod_ii,      true, ndk_Mod_whole, NULL };
static const NDBinaryKernel NDKB_Quotient =
    { ndk_Quotient_c, true, ndk_Quotient_ri,  ndk_Quotient_ii, true, NULL, NULL };

/* ---- special functions: real machine kernels ---------------------------- */
/* These cover the common real-array case at C speed via libc; a complex NDArray
 * (no cplx kernel) or an out-of-fast-domain element (order non-integer, pole)
 * declines and the evaluator degrades that call to the exact List path. */

static bool ndk_Gamma_r(double x, double* o) { double v = tgamma(x); *o = v; return isfinite(v); }
/* Complex Gamma is the interpreter's own Lanczos series, shared rather than
 * reimplemented — see gamma.h.  real_closed stays true, so a real argument
 * still takes tgamma above and only a complex one reaches the Lanczos path. */
static const NDUnaryKernel NDKU_Gamma = { gamma_machine_complex, ndk_Gamma_r, true, false, NULL, NULL, false, NULL };

/* Factorial[x] = Gamma[x + 1] over a real array (matches builtin_factorial's
 * real path). A pole (negative integer) yields non-finite -> declines. */
static bool ndk_Factorial_r(double x, double* o) { double v = tgamma(x + 1.0); *o = v; return isfinite(v); }
static const NDUnaryKernel NDKU_Factorial = { NULL, ndk_Factorial_r, true, false, NULL, NULL, false, NULL };

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
static const NDUnaryKernel NDKU_LogGamma = { loggamma_machine_complex, ndk_LogGamma_r, true, false, NULL, NULL, false, NULL };

static bool ndk_Erf_r(double x, double* o)  { return sf_erf_r(x, o); }
static bool ndk_Erfc_r(double x, double* o) { return sf_erfc_r(x, o); }
/* Both arms are the sf_machine_error.c kernels (the scalar builtin and the
 * Compile opcodes OP_ERF_R / OP_ERFC_R / OP_VERF_R / OP_VERFC_R call the same
 * helpers, so every surface agrees bit for bit). */
static const NDUnaryKernel NDKU_Erf  = { sf_erf_c,  ndk_Erf_r,  true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_Erfc = { sf_erfc_c, ndk_Erfc_r, true, false, NULL, NULL, false, NULL };

/* Bessel J/Y/I/K over real AND complex arrays (real order): the shared
 * dispatch in special_functions/sf_machine_bessel.c, which the scalar builtins
 * call too, so the surfaces agree bit for bit. A real element whose value is
 * complex (BesselJ[1/2, -2.]) declines -> List path. */
static bool ndk_BesselJ_c(double nre, double nim, double zre, double zim, double* rr, double* ri)
{ return sf_bessel_nd(SF_BESSEL_J, nre, nim, zre, zim, rr, ri); }
static bool ndk_BesselY_c(double nre, double nim, double zre, double zim, double* rr, double* ri)
{ return sf_bessel_nd(SF_BESSEL_Y, nre, nim, zre, zim, rr, ri); }
static bool ndk_BesselI_c(double nre, double nim, double zre, double zim, double* rr, double* ri)
{ return sf_bessel_nd(SF_BESSEL_I, nre, nim, zre, zim, rr, ri); }
static bool ndk_BesselK_c(double nre, double nim, double zre, double zim, double* rr, double* ri)
{ return sf_bessel_nd(SF_BESSEL_K, nre, nim, zre, zim, rr, ri); }
static const NDBinaryKernel NDKB_BesselJ = { ndk_BesselJ_c, true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_BesselY = { ndk_BesselY_c, true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_BesselI = { ndk_BesselI_c, true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_BesselK = { ndk_BesselK_c, true, NULL, NULL, false, NULL, NULL };

/* Beta[a, b], real and complex: the Gamma family's kernel (sf_machine_gamma.c),
 * shared with the scalar builtin so the surfaces agree bit for bit. */
static bool ndk_Beta_c(double are, double aim, double bre, double bim,
                       double* rr, double* ri) {
    return sfg_beta_c(are, aim, bre, bim, rr, ri);
}
static const NDBinaryKernel NDKB_Beta = { ndk_Beta_c, true, NULL, NULL, false, NULL, NULL };

/* The Legendre elliptic integrals over real arrays, by Carlson's symmetric
 * forms (src/special_functions/elliptic.c). PARAMETER convention m = k^2.
 *
 * Real arguments only: a complex element declines, and so does one outside the
 * real principal domain (m > 1 for the complete forms, or 1 - m Sin[phi]^2 < 0
 * for the incomplete), because there the value is genuinely complex and only the
 * FLINT/Arb path places the branch correctly. Declining ABANDONS the array and
 * the List path answers exactly -- slower, never wrong.
 *
 * EllipticE carries BOTH a unary and a binary kernel, because the head is
 * arity-overloaded: EllipticE[m] is complete and EllipticE[phi, m] incomplete.
 * EllipticPi likewise: a BINARY kernel for the complete form, and an N-ARY one
 * for the incomplete, which the element-wise NDArray layer cannot use (it tops
 * out at arity 2) but Compile[]'s OP_KERNN does -- so registering it is what
 * lowers EllipticPi[n, phi, m] inside Compile[]. Both decline p <= 0, where the
 * value is genuinely complex rather than a real principal value; the builtin's
 * delist fallback keeps the NDArray surface correct there and for arity 3. */
static bool ndk_EllipticK_r(double m, double* o) { return elliptic_machine_k(m, o); }
static bool ndk_EllipticE_r(double m, double* o) { return elliptic_machine_e_complete(m, o); }
/* Complex arms (sf_machine_ellk.c, complex AGM): the same kernels the scalar
 * EllipticK[z] / EllipticE[z] call for a machine-complex z. A REAL buffer keeps
 * real_closed: its real arm declines m > 1 (complex value), which abandons the
 * buffer to the List path, whose scalar call then takes the complex kernel. */
static const NDUnaryKernel NDKU_EllipticK = { sf_machine_ellk_k_c, ndk_EllipticK_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_EllipticE = { sf_machine_ellk_e_c, ndk_EllipticE_r, true, false, NULL, NULL, false, NULL };

static bool ndk_EllipticF_c(double pre, double pim, double mre, double mim,
                            double* rr, double* ri) {
    if (pim != 0.0 || mim != 0.0) return false;
    if (!elliptic_machine_f(pre, mre, rr)) return false;
    *ri = 0.0; return isfinite(*rr);
}
static bool ndk_EllipticEInc_c(double pre, double pim, double mre, double mim,
                               double* rr, double* ri) {
    if (pim != 0.0 || mim != 0.0) return false;
    if (!elliptic_machine_e_inc(pre, mre, rr)) return false;
    *ri = 0.0; return isfinite(*rr);
}
static const NDBinaryKernel NDKB_EllipticF = { ndk_EllipticF_c,    true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_EllipticE = { ndk_EllipticEInc_c, true, NULL, NULL, false, NULL, NULL };

static bool ndk_EllipticPi_c(double nre, double nim, double mre, double mim,
                             double* rr, double* ri) {
    if (nim != 0.0 || mim != 0.0) return false;
    if (!elliptic_machine_pi(nre, mre, rr)) return false;
    *ri = 0.0; return isfinite(*rr);
}
static const NDBinaryKernel NDKB_EllipticPi = { ndk_EllipticPi_c, true, NULL, NULL, false, NULL, NULL };

static bool ndk_EllipticPiInc_n(const double* re, const double* im, size_t n,
                                double* rr, double* ri) {
    if (n != 3) return false;
    for (size_t i = 0; i < 3; i++) if (im[i] != 0.0) return false;
    if (!elliptic_machine_pi_inc(re[0], re[1], re[2], rr)) return false;
    *ri = 0.0; return isfinite(*rr);
}
static const NDNaryKernel NDKN_EllipticPi = { ndk_EllipticPiInc_n, 3, true };

/* Exponential-integral family and friends.  The kernels (sf_machine_expint.c)
 * are the ones the scalar builtins call, so the surfaces agree bit for bit.
 * All are registered real-closed: the `real` function DECLINES (returns false)
 * for an argument whose true value is complex or infinite (Ci/Chi/li of a
 * negative argument, the poles), which abandons the buffer to the List path;
 * a complex buffer goes through the `cplx` kernel. */
static bool ndk_ExpIntegralEi_r(double x, double* o) { return sf_machine_ei(x, o); }
static bool ndk_LogIntegral_r(double x, double* o)   { return sf_machine_li(x, o); }
static bool ndk_SinIntegral_r(double x, double* o)   { return sf_machine_si(x, o); }
static bool ndk_CosIntegral_r(double x, double* o)   { return sf_machine_ci(x, o); }
static bool ndk_SinhIntegral_r(double x, double* o)  { return sf_machine_shi(x, o); }
static bool ndk_CoshIntegral_r(double x, double* o)  { return sf_machine_chi(x, o); }
static bool ndk_Sinc_r(double x, double* o)          { return sf_sinc_r(x, o); }
/* InverseErf/InverseErfc: the piecewise-polynomial kernels (decline at the
 * poles +-1 / 0, 2 and outside the real domain). */
static bool ndk_InverseErf_r(double x, double* o)    { return sf_inverse_erf_r(x, o); }
static bool ndk_InverseErfc_r(double x, double* o)   { return sf_inverse_erfc_r(x, o); }
static const NDUnaryKernel NDKU_ExpIntegralEi = { expintegralei_machine_complex, ndk_ExpIntegralEi_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_LogIntegral   = { logintegral_machine_complex, ndk_LogIntegral_r,   true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_SinIntegral   = { sinintegral_machine_complex, ndk_SinIntegral_r,   true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_CosIntegral   = { cosintegral_machine_complex, ndk_CosIntegral_r,   true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_SinhIntegral  = { sinhintegral_machine_complex, ndk_SinhIntegral_r,  true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_CoshIntegral  = { coshintegral_machine_complex, ndk_CoshIntegral_r,  true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_Sinc          = { sf_sinc_c, ndk_Sinc_r,     true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_InverseErf    = { NULL, ndk_InverseErf_r,    true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_InverseErfc   = { NULL, ndk_InverseErfc_r,   true, false, NULL, NULL, false, NULL };
/* CubeRoot: Surd[x, 3]'s kernel (src/special_functions/surd.c), real-closed;
 * no complex arm (a non-real element is CubeRoot::preal -> the List path). */
/* RiemannR: the double-double Gram series (src/special_functions/riemannr.c).
 * Not real-closed (x < 0 is complex); the real arm declines x < 1/e, where the
 * alternating series cancels, so such a buffer goes to the List path and the
 * scalar builtin's MPFR route -- the same split the scalar path makes. */
static const NDUnaryKernel NDKU_RiemannR      = { riemannr_machine_complex, sf_riemannr_r, false, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_CubeRoot      = { NULL, sf_cuberoot_r,       true, false, NULL, NULL, false, NULL };

static bool ndk_Erfi_r(double x, double* o)          { return sf_erfi_r(x, o); }
static bool ndk_ProductLog_r(double x, double* o)    { return sf_machine_productlog(x, o); }
static bool ndk_FresnelC_r(double x, double* o)      { return sf_fresnel_c_r(x, o); }
static bool ndk_FresnelS_r(double x, double* o)      { return sf_fresnel_s_r(x, o); }
static bool ndk_PolyGamma_r(double x, double* o)     { return sfg_digamma(x, o); }
static bool ndk_PolyGamma_u(double xr, double xi, double* rr, double* ri) {
    return sfg_polygamma_c(0, xr, xi, rr, ri);
}
static bool ndk_HarmonicNumber_r(double x, double* o){ return sf_machine_harmonic(x, o); }
static bool ndk_Zeta_r(double x, double* o)          { return sf_machine_zeta(x, o); }
static bool ndk_Fibonacci_r(double x, double* o)     { return sf_machine_fibonacci(x, o); }
static bool ndk_LucasL_r(double x, double* o)        { return sf_machine_lucasl(x, o); }
static const NDUnaryKernel NDKU_Erfi           = { sf_erfi_c, ndk_Erfi_r,      true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_ProductLog     = { sf_machine_productlog_c, ndk_ProductLog_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_FresnelC       = { sf_fresnel_c_c, ndk_FresnelC_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_FresnelS       = { sf_fresnel_s_c, ndk_FresnelS_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_PolyGamma      = { ndk_PolyGamma_u, ndk_PolyGamma_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_HarmonicNumber = { NULL, ndk_HarmonicNumber_r, true, false, NULL, NULL, false, NULL };
/* Zeta's complex arm: the same sf_zeta_c the scalar builtin calls. */
static bool ndk_Zeta_c(double sr, double si, double* rr, double* ri) { return sf_zeta_c(sr, si, rr, ri); }
static const NDUnaryKernel NDKU_Zeta           = { ndk_Zeta_c, ndk_Zeta_r,     true, false, NULL, NULL, false, NULL };
/* RiemannSiegelZ / RiemannSiegelTheta: the scalar builtins' machine helpers
 * (arb at a 53-bit target), real and complex arms. */
static bool ndk_RSZ_c(double a, double b, double* r, double* i) { return sf_machine_riemann_siegel_z_c(a, b, r, i); }
static bool ndk_RST_c(double a, double b, double* r, double* i) { return sf_machine_riemann_siegel_theta_c(a, b, r, i); }
static const NDUnaryKernel NDKU_RiemannSiegelZ     = { ndk_RSZ_c, sf_machine_riemann_siegel_z,     true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_RiemannSiegelTheta = { ndk_RST_c, sf_machine_riemann_siegel_theta, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_Fibonacci      = { NULL, ndk_Fibonacci_r,      true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_LucasL         = { NULL, ndk_LucasL_r,         true, false, NULL, NULL, false, NULL };
/* Airy: real and complex arms of the shared kernels (sf_machine_bessel.c). */
static bool ndk_AiryAi_c(double a, double b, double* r, double* i)      { return sf_airy_nd(SF_AIRY_AI, a, b, r, i); }
static bool ndk_AiryBi_c(double a, double b, double* r, double* i)      { return sf_airy_nd(SF_AIRY_BI, a, b, r, i); }
static bool ndk_AiryAiPrime_c(double a, double b, double* r, double* i) { return sf_airy_nd(SF_AIRY_AIP, a, b, r, i); }
static bool ndk_AiryBiPrime_c(double a, double b, double* r, double* i) { return sf_airy_nd(SF_AIRY_BIP, a, b, r, i); }
static const NDUnaryKernel NDKU_AiryAi      = { ndk_AiryAi_c,      sf_airy_ai_r,  true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_AiryBi      = { ndk_AiryBi_c,      sf_airy_bi_r,  true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_AiryAiPrime = { ndk_AiryAiPrime_c, sf_airy_aip_r, true, false, NULL, NULL, false, NULL };
static const NDUnaryKernel NDKU_AiryBiPrime = { ndk_AiryBiPrime_c, sf_airy_bip_r, true, false, NULL, NULL, false, NULL };

/* PolyGamma[n, x] and HurwitzZeta[s, a].  PolyGamma matters as a BINARY kernel
 * even when written unary: the evaluator canonicalises PolyGamma[x] to
 * PolyGamma[0, x] before anything downstream sees it. */
/* Integer order n, real or complex x: the Gamma family's kernel. */
static bool ndk_PolyGamma_c(double nre, double nim, double xre, double xim,
                            double* rr, double* ri) {
    if (nim != 0.0 || nre < 0.0 || nre != floor(nre) ||
        nre > (double)SFG_POLYGAMMA_MAX_N) return false;
    return sfg_polygamma_c((int)nre, xre, xim, rr, ri);
}
/* HurwitzZeta: real and complex arms are both sf_hurwitz_c (which hands an
 * all-real call to the real kernel), exactly as the scalar builtin does. */
static bool ndk_HurwitzZeta_c(double sre, double sim, double are, double aim,
                              double* rr, double* ri) {
    return sf_hurwitz_c(sre, sim, are, aim, rr, ri);
}
/* Zeta[s, a]: sf_zeta2_c, the symmetric convention the scalar builtin uses
 * (it differs from HurwitzZeta for a <= 0). */
static bool ndk_Zeta2_c(double sre, double sim, double are, double aim,
                        double* rr, double* ri) {
    return sf_zeta2_c(sre, sim, are, aim, rr, ri);
}
static bool ndk_Pochhammer_c(double are, double aim, double nre, double nim,
                             double* rr, double* ri) {
    return sfg_pochhammer_c(are, aim, nre, nim, rr, ri);
}
static bool ndk_Binomial_c(double nre, double nim, double kre, double kim,
                           double* rr, double* ri) {
    if (nim != 0.0 || kim != 0.0) return false;
    double v; if (!sf_machine_binomial(nre, kre, &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* LegendreP / LegendreQ: the shared real + complex helpers of
 * sf_machine_plog.c (the scalar builtins call the same ones).  Registered
 * real_closed, so a REAL element whose value leaves the real axis -- LegendreQ
 * at a real |x| > 1 -- must DECLINE here rather than have its imaginary part
 * dropped; the List path then answers it (complex) exactly as the scalar does. */
static bool ndk_LegendreP_c(double nre, double nim, double xre, double xim,
                            double* rr, double* ri) {
    if (nim != 0.0) return false;
    if (xim == 0.0 && nre >= 0.0) {          /* the real hot path, one call deep */
        if (!sf_machine_legendre_p(nre, xre, rr)) return false;
        *ri = 0.0; return true;
    }
    return sf_machine_legendre_p_c(nre, xre, xim, rr, ri);
}
static bool ndk_LegendreQ_c(double nre, double nim, double xre, double xim,
                            double* rr, double* ri) {
    if (nim != 0.0) return false;
    if (xim == 0.0) {
        if (!sf_machine_legendre_q(nre, xre, rr)) return false;   /* |x| >= 1 declines */
        *ri = 0.0; return true;
    }
    return sf_machine_legendre_q_c(nre, xre, xim, rr, ri);
}
static const NDBinaryKernel NDKB_Pochhammer = { ndk_Pochhammer_c, true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_Binomial   = { ndk_Binomial_c,   true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_LegendreP  = { ndk_LegendreP_c,  true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_LegendreQ  = { ndk_LegendreQ_c,  true, NULL, NULL, false, NULL, NULL };

/* ExpIntegralE[n, z]: integer order n >= 0 (sf_machine_expint.c, the kernels
 * the scalar builtin calls). Real x > 0 stays real; complex z, and real x < 0
 * (complex-valued), go through the complex kernel. Registered real-closed, so a
 * real buffer holding an x < 0 abandons to the List path, which answers it. */
static bool ndk_ExpIntegralE_c(double nre, double nim, double xre, double xim,
                               double* rr, double* ri) {
    if (nim != 0.0) return false;
    if (xim == 0.0) {
        double v;
        if (sf_machine_expn(nre, xre, &v)) { *rr = v; *ri = 0.0; return true; }
        if (xre >= 0.0) return false;            /* x = 0: exact 1/(n-1) or pole */
    }
    return sf_machine_expn_c(nre, xre, xim, rr, ri);
}
static const NDBinaryKernel NDKB_ExpIntegralE = { ndk_ExpIntegralE_c, true, NULL, NULL, false, NULL, NULL };

/* Anger-Weber and Kelvin: (nu, z) shims over the shared machine helpers of
 * special_functions/sf_machine_angerkelvin.c, which the scalar builtins call
 * too -- so scalar, packed, NDArray and Compile[] agree bit for bit. Each
 * declines where the value is complex (-> List path -> arb). The Kelvin unary
 * kernels are the order-0 one-argument forms KelvinBer[x], ... */
static bool ndk_anger_r(double nu, double z, double* o) { return sf_machine_anger_weber(0, nu, z, o); }
static bool ndk_weber_r(double nu, double z, double* o) { return sf_machine_anger_weber(1, nu, z, o); }
static bool ndk_ber_r(double nu, double x, double* o) { return sf_machine_kelvin(0, nu, x, o); }
static bool ndk_bei_r(double nu, double x, double* o) { return sf_machine_kelvin(1, nu, x, o); }
static bool ndk_ker_r(double nu, double x, double* o) { return sf_machine_kelvin(2, nu, x, o); }
static bool ndk_kei_r(double nu, double x, double* o) { return sf_machine_kelvin(3, nu, x, o); }
#define NDK_KELVIN0(NAME, W)                                                   \
static bool ndk_##NAME##0_r(double x, double* o) { return sf_machine_kelvin(W, 0.0, x, o); } \
static bool ndk_##NAME##0_c(double re, double im, double* ore, double* oim) {  \
    if (im != 0.0) return false;                                               \
    if (!sf_machine_kelvin(W, 0.0, re, ore)) return false;                     \
    *oim = 0.0; return true;                                                   \
}                                                                              \
static const NDUnaryKernel NDKU_##NAME = { ndk_##NAME##0_c, ndk_##NAME##0_r, true, false, NULL, NULL, false, NULL };
NDK_KELVIN0(KelvinBer, 0)
NDK_KELVIN0(KelvinBei, 1)
NDK_KELVIN0(KelvinKer, 2)
NDK_KELVIN0(KelvinKei, 3)
#undef NDK_KELVIN0

#define NDK_BIN2(NAME, FN)                                                     \
static bool ndk_##NAME##_c(double are, double aim, double bre, double bim,     \
                           double* rr, double* ri) {                           \
    if (aim != 0.0 || bim != 0.0) return false;                                \
    double v; if (!FN(are, bre, &v)) return false;                             \
    *rr = v; *ri = 0.0; return true;                                           \
}                                                                              \
static const NDBinaryKernel NDKB_##NAME = { ndk_##NAME##_c, true, NULL, NULL, false, NULL, NULL };
NDK_BIN2(QPochhammer,       sf_machine_qpochhammer)
NDK_BIN2(StruveH,           sf_machine_struve_h)
NDK_BIN2(StruveL,           sf_machine_struve_l)
NDK_BIN2(AngerJ,            ndk_anger_r)
NDK_BIN2(WeberE,            ndk_weber_r)
NDK_BIN2(KelvinBer,         ndk_ber_r)
NDK_BIN2(KelvinBei,         ndk_bei_r)
NDK_BIN2(KelvinKer,         ndk_ker_r)
NDK_BIN2(KelvinKei,         ndk_kei_r)
NDK_BIN2(Surd,              sf_surd_r)          /* Surd[x, n]: real root, integer n */
NDK_BIN2(BesselJZero,       sf_bessel_j_zero_r) /* BesselJZero[n, k]: real n, integer k */
NDK_BIN2(BesselYZero,       sf_bessel_y_zero_r) /* BesselYZero[n, k]: real n >= 0 */
NDK_BIN2(MathieuCharacteristicA, sf_mathieu_char_a_r) /* real r, real q */
NDK_BIN2(MathieuCharacteristicB, sf_mathieu_char_b_r)
#undef NDK_BIN2

/* PolyLog[s, z]: written out because it has a machine COMPLEX arm. A wholly
 * real element uses the real kernel and DECLINES where the value leaves the
 * real axis (z > 1), since the descriptor is real_closed; a complex element
 * uses sf_polylog_c. Same helpers as the scalar builtin. */
static bool ndk_PolyLog_c(double sre, double sim, double zre, double zim,
                          double* rr, double* ri) {
    if (sim == 0.0 && zim == 0.0) {
        double v; if (!sf_machine_polylog(sre, zre, &v)) return false;
        *rr = v; *ri = 0.0; return true;
    }
    return sf_polylog_c(sre, sim, zre, zim, rr, ri);
}
static const NDBinaryKernel NDKB_PolyLog = { ndk_PolyLog_c, true, NULL, NULL, false, NULL, NULL };
/* Hypergeometric0F1[b, z]: written out (not NDK_BIN2) for its machine COMPLEX
 * arm -- a complex z or b goes to sf_machine_hyper0f1_c instead of declining, so
 * a complex buffer computes on the buffer. Same helpers as the scalar builtin. */
static bool ndk_Hypergeometric0F1_c(double are, double aim, double bre, double bim,
                                    double* ore, double* oim) {
    if (aim == 0.0 && bim == 0.0) {
        double v; if (!sf_machine_hyper0f1(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_hyper0f1_c(are, aim, bre, bim, ore, oim);
}
static const NDBinaryKernel NDKB_Hypergeometric0F1 =
    { ndk_Hypergeometric0F1_c, true, NULL, NULL, false, NULL, NULL };

/* ChebyshevT and ChebyshevU are written out rather than generated by NDK_BIN2
 * because each has a machine COMPLEX arm as well (the generated form declines
 * all complex input). The order `a` must still be a real integer; the argument
 * `b` may be complex, so a packed/visible complex buffer computes on the buffer
 * instead of materialising to the scalar path. Same raw-double helpers as the
 * scalar dispatch, so the surfaces agree bit-for-bit. */
static bool ndk_ChebyshevT_c(double are, double aim, double bre, double bim,
                             double* ore, double* oim) {
    if (aim != 0.0) return false;                 /* order must be a real integer */
    if (bim == 0.0) {                              /* real argument */
        double v; if (!sf_machine_chebyshev_t(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_chebyshev_t_c(are, bre, bim, ore, oim);   /* complex argument */
}
static const NDBinaryKernel NDKB_ChebyshevT = { ndk_ChebyshevT_c, true, NULL, NULL, false, NULL, NULL };

static bool ndk_ChebyshevU_c(double are, double aim, double bre, double bim,
                             double* ore, double* oim) {
    if (aim != 0.0) return false;                 /* order must be a real integer */
    if (bim == 0.0) {                              /* real argument */
        double v; if (!sf_machine_chebyshev_u(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_chebyshev_u_c(are, bre, bim, ore, oim);   /* complex argument */
}
static const NDBinaryKernel NDKB_ChebyshevU = { ndk_ChebyshevU_c, true, NULL, NULL, false, NULL, NULL };

/* HermiteH[n, x] (physicists'). Hand-written rather than via NDK_BIN2 because it
 * carries a machine COMPLEX arm: the order `a` must be a real integer, but the
 * argument `b` = x may be complex, so a packed/visible complex buffer computes
 * on the buffer. Same raw-double helpers as the scalar dispatch, so the surfaces
 * agree bit-for-bit. */
static bool ndk_HermiteH_c(double are, double aim, double bre, double bim,
                           double* ore, double* oim) {
    if (aim != 0.0) return false;                 /* order must be a real integer */
    if (bim == 0.0) {                              /* real argument */
        double v; if (!sf_machine_hermite_h(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_hermite_h_c(are, bre, bim, ore, oim);   /* complex argument */
}
static const NDBinaryKernel NDKB_HermiteH = { ndk_HermiteH_c, true, NULL, NULL, false, NULL, NULL };

/* ParabolicCylinderD[nu, z], real or complex order and argument, through
 * sf_machine_pcfd / sf_machine_pcfd_cc -- the scalar builtin's machine helpers,
 * so the surfaces agree bit for bit. Declines only a value outside the normal
 * double range (the buffer is then delisted and the MPFR path answers). */
static bool ndk_ParabolicCylinderD_c(double are, double aim, double bre, double bim,
                                     double* ore, double* oim) {
    if (aim != 0.0) return sf_machine_pcfd_cc(are, aim, bre, bim, ore, oim);
    if (bim == 0.0) {
        double v; if (!sf_machine_pcfd(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_pcfd_c(are, bre, bim, ore, oim);
}
static const NDBinaryKernel NDKB_ParabolicCylinderD = { ndk_ParabolicCylinderD_c, true, NULL, NULL, false, NULL, NULL };

/* GegenbauerC[n, z] = 2 ChebyshevT[n, z]/n: sf_machine_gegenbauer2_c, the
 * scalar builtin's helper. Complex beyond z = -1 for a non-integer n, so not
 * real-closed. */
static bool ndk_GegenbauerC_c(double are, double aim, double bre, double bim,
                              double* ore, double* oim) {
    return sf_machine_gegenbauer2_c(are, aim, bre, bim, ore, oim);
}
static const NDBinaryKernel NDKB_GegenbauerC = { ndk_GegenbauerC_c, false, NULL, NULL, false, NULL, NULL };

/* LaguerreL[n, x] (ordinary Laguerre, a = 0). Hand-written rather than via
 * NDK_BIN2 because it carries a machine COMPLEX arm: the order `a` (here the
 * polynomial degree n) must be a real integer, but the argument `b` = x may be
 * complex, so a packed/visible complex buffer computes on the buffer. Same
 * raw-double helpers as the scalar dispatch, so the surfaces agree bit-for-bit.
 * (The generalized LaguerreL[n, a, x] uses the 3-arg N-ary kernel below.) */
static bool ndk_LaguerreL_c(double are, double aim, double bre, double bim,
                            double* ore, double* oim) {
    if (aim != 0.0) return false;                 /* order must be a real integer */
    if (bim == 0.0) {                              /* real argument */
        double v; if (!sf_machine_laguerre_l(are, bre, &v)) return false;
        *ore = v; *oim = 0.0; return true;
    }
    return sf_machine_laguerre_l_c(are, bre, bim, ore, oim);   /* complex argument */
}
static const NDBinaryKernel NDKB_LaguerreL = { ndk_LaguerreL_c, true, NULL, NULL, false, NULL, NULL };

/* N-ary kernels: arity >= 3, which the unary/binary descriptors cannot express.
 * Consumed by the Compile[] VM; the NDArray element-wise path still degrades for
 * these heads, since "one array plus broadcast scalars" is a different question. */
static bool ndk_LerchPhi_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 3) return false;
    /* sf_lerch_c hands an all-real call to sf_machine_lerchphi, as the scalar
     * builtin does. */
    return sf_lerch_c(re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
/* The hypergeometric n-ary kernels: a real arm (real result, the Compile[] real
 * VM) and, when any argument has an imaginary part, the machine-complex twin
 * (sf_machine_hyper.h). Real arguments never take the complex arm -- a real
 * 2F1 on its cut z > 1 declines here rather than handing a complex value to a
 * real-typed consumer. */
static bool ndk_any_im(const double* im, size_t n) {
    for (size_t i = 0; i < n; i++) if (im[i] != 0.0) return true;
    return false;
}
static bool ndk_Hyper1F1_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 3) return false;
    if (ndk_any_im(im, 3))
        return sf_machine_hyper1f1_c(re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
    double v; if (!sf_machine_hyper1f1(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static bool ndk_Hyper2F1_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 4) return false;
    if (ndk_any_im(im, 4))
        return sf_machine_hyper2f1_c(re[0], im[0], re[1], im[1], re[2], im[2],
                                     re[3], im[3], rr, ri);
    double v; if (!sf_machine_hyper2f1(re[0], re[1], re[2], re[3], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* Generalized LaguerreL[n, a, x]: re[0] = n, re[1] = a, re[2] = x. Real only
 * (the Compile VM consumes it); a complex argument falls to the interpreter. */
static bool ndk_LaguerreL_n(const double* re, const double* im, size_t n,
                            double* rr, double* ri) {
    if (n != 3 || im[0] || im[1] || im[2]) return false;
    double v; if (!sf_machine_laguerre_l_assoc(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* JacobiP[n, a, b, x]: re[0..3] = n, a, b, x. Real n, a, b; x real or complex,
 * through sf_machine_jacobi_p[_c] -- the scalar builtin's machine helper. */
static bool ndk_JacobiP_n(const double* re, const double* im, size_t n,
                          double* rr, double* ri) {
    if (n != 4 || im[0] || im[1] || im[2]) return false;
    if (im[3] != 0.0) return sf_machine_jacobi_p_c(re[0], re[1], re[2], re[3], im[3], rr, ri);
    double v; if (!sf_machine_jacobi_p(re[0], re[1], re[2], re[3], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* WhittakerM[k, m, z]: re[0] = k, re[1] = m, re[2] = z. Real only (the Compile VM
 * consumes it); a complex argument falls to the interpreter. */
static bool ndk_WhittakerM_n(const double* re, const double* im, size_t n,
                             double* rr, double* ri) {
    if (n != 3 || im[0] || im[1]) return false;
    if (im[2] != 0.0) return sf_machine_whittaker_m_c(re[0], re[1], re[2], im[2], rr, ri);
    double v; if (!sf_machine_whittaker_m(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* MathieuC / MathieuCPrime [a, q, z]: complex in every slot. A real call can
 * have a purely imaginary value (an instability band); the kernel returns it
 * whole -- ndarray_map_nary widens the output buffer to complex, and the
 * Compile VM's real-slot KERNN hands such a call to the interpreter (oi != 0),
 * so the imaginary part is never dropped. The plan cache behind
 * sf_machine_mathieu is thread-local, so nd_parallel_for may split the map. */
static bool ndk_mathieu_n(int deriv, const double* re, const double* im, size_t n,
                          double* rr, double* ri) {
    if (n != 3) return false;
    return sf_machine_mathieu(deriv, re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
static bool ndk_MathieuC_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    return ndk_mathieu_n(0, re, im, n, rr, ri);
}
static bool ndk_MathieuCPrime_n(const double* re, const double* im, size_t n,
                                double* rr, double* ri) {
    return ndk_mathieu_n(1, re, im, n, rr, ri);
}
/* MathieuS / MathieuSPrime: the odd twins, the same plan cache. */
static bool ndk_MathieuS_n(const double* re, const double* im, size_t n,
                           double* rr, double* ri) {
    if (n != 3) return false;
    return sf_machine_mathieu_fn(1, 0, re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
static bool ndk_MathieuSPrime_n(const double* re, const double* im, size_t n,
                                double* rr, double* ri) {
    if (n != 3) return false;
    return sf_machine_mathieu_fn(1, 1, re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
/* HypergeometricU[a, b, z]: re[0] = a, re[1] = b, re[2] = z. Real only. */
static bool ndk_HypergeometricU_n(const double* re, const double* im, size_t n,
                                  double* rr, double* ri) {
    if (n != 3 || im[0] || im[1]) return false;
    if (im[2] != 0.0) return sf_machine_hyperu_c(re[0], re[1], re[2], im[2], rr, ri);
    double v; if (!sf_machine_hyperu(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* WhittakerW[k, m, z]: re[0] = k, re[1] = m, re[2] = z. Real only. */
static bool ndk_WhittakerW_n(const double* re, const double* im, size_t n,
                             double* rr, double* ri) {
    if (n != 3 || im[0] || im[1]) return false;
    if (im[2] != 0.0) return sf_machine_whittaker_w_c(re[0], re[1], re[2], im[2], rr, ri);
    double v; if (!sf_machine_whittaker_w(re[0], re[1], re[2], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
/* PFQ arrives already flattened by the compiler's lowering: re[0] carries p, so
 * q is whatever is left after the p upper parameters and the final z. */
static bool ndk_HypergeometricPFQ_n(const double* re, const double* im, size_t n,
                                    double* rr, double* ri) {
    if (n < 2 || im[0] != 0.0) return false;
    size_t p = (size_t)re[0];
    if (n < p + 2) return false;
    size_t q = n - p - 2;
    if (ndk_any_im(im, n))
        return sf_machine_pfq_c(re + 1, im + 1, p, re + 1 + p, im + 1 + p, q,
                                re[n - 1], im[n - 1], rr, ri);
    double v;
    if (!sf_machine_pfq(re + 1, p, re + 1 + p, q, re[n - 1], &v)) return false;
    *rr = v; *ri = 0.0; return true;
}
static const NDNaryKernel NDKN_HypergeometricPFQ = { ndk_HypergeometricPFQ_n, 0, true };
/* MeijerG arrives flattened by the compiler's lowering in the layout
 * sf_mg_eval_flat documents: [shape code, r, a_1..a_p, b_1..b_q, z]. The
 * helper is the one behind the scalar builtin and the buffer path (it falls
 * back to arb internally on a declined point, so it never answers NaN). */
static bool ndk_MeijerG_n(const double* re, const double* im, size_t n,
                          double* rr, double* ri) {
    return sf_mg_eval_flat(re, im, n, rr, ri);
}
static const NDNaryKernel NDKN_MeijerG = { ndk_MeijerG_n, 0, true };
static const NDNaryKernel NDKN_LerchPhi            = { ndk_LerchPhi_n, 3, true };
/* PolyLog[n, p, z], Nielsen's S_{n,p}(z): sf_nielsen_c, the scalar builtin's
 * helper. Complex beyond z = 1 even for real arguments, so not real-closed. */
static bool ndk_PolyLog_n(const double* re, const double* im, size_t n,
                          double* rr, double* ri) {
    if (n != 3) return false;
    return sf_nielsen_c(re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
static const NDNaryKernel NDKN_PolyLog              = { ndk_PolyLog_n, 3, false };
/* DirichletL[k, j, s]: k, j integer labels (an imaginary part, or a
 * non-integer / out-of-range label, declines), s complex. The helper is the one
 * behind the scalar builtin and its buffer loop; a complex character makes the
 * value complex even at real s, so the kernel is not real-closed. */
static bool ndk_DirichletL_n(const double* re, const double* im, size_t n,
                             double* rr, double* ri) {
    if (n != 3 || im[0] != 0.0 || im[1] != 0.0) return false;
    return sf_dirichlet_l_c(re[0], re[1], re[2], im[2], rr, ri);
}
static const NDNaryKernel NDKN_DirichletL          = { ndk_DirichletL_n, 3, false };
static const NDNaryKernel NDKN_Hypergeometric1F1   = { ndk_Hyper1F1_n, 3, true };
static const NDNaryKernel NDKN_Hypergeometric2F1   = { ndk_Hyper2F1_n, 4, true };
static const NDNaryKernel NDKN_LaguerreL           = { ndk_LaguerreL_n, 3, true };
static const NDNaryKernel NDKN_JacobiP             = { ndk_JacobiP_n, 4, true };
static const NDNaryKernel NDKN_WhittakerM          = { ndk_WhittakerM_n, 3, true };
static const NDNaryKernel NDKN_HypergeometricU     = { ndk_HypergeometricU_n, 3, true };
static const NDNaryKernel NDKN_WhittakerW          = { ndk_WhittakerW_n, 3, true };
/* GegenbauerC[n, l, z]: complex in every slot, through sf_machine_gegenbauer_c
 * (the scalar builtin's helper). A real call can be complex (z < -1 with a
 * non-integer n), so not real-closed. */
static bool ndk_GegenbauerC_n(const double* re, const double* im, size_t n,
                              double* rr, double* ri) {
    if (n != 3) return false;
    return sf_machine_gegenbauer_c(re[0], im[0], re[1], im[1], re[2], im[2], rr, ri);
}
static const NDNaryKernel NDKN_GegenbauerC         = { ndk_GegenbauerC_n, 3, false };
static const NDNaryKernel NDKN_MathieuC            = { ndk_MathieuC_n, 3, true };
static const NDNaryKernel NDKN_MathieuCPrime       = { ndk_MathieuCPrime_n, 3, true };
static const NDNaryKernel NDKN_MathieuS            = { ndk_MathieuS_n, 3, true };
static const NDNaryKernel NDKN_MathieuSPrime       = { ndk_MathieuSPrime_n, 3, true };
#define REG_N(NAME) symtab_set_ndarray_nary_kernel(#NAME, &NDKN_##NAME)
static const NDBinaryKernel NDKB_PolyGamma   = { ndk_PolyGamma_c,   true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_HurwitzZeta = { ndk_HurwitzZeta_c, true, NULL, NULL, false, NULL, NULL };
static const NDBinaryKernel NDKB_Zeta        = { ndk_Zeta2_c, true, NULL, NULL, false, NULL, NULL };

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
    REG_U(Exp);  REG_U(Log);  REG_U(Log10);  REG_U(Log2);
    REG_U(ArcSin); REG_U(ArcCos); REG_U(ArcTan);
    REG_U(ArcCot); REG_U(ArcSec); REG_U(ArcCsc);
    REG_U(ArcSinh); REG_U(ArcCosh); REG_U(ArcTanh);
    REG_U(ArcCoth); REG_U(ArcSech); REG_U(ArcCsch);
    symtab_set_ndarray_unary_kernel("Abs", &NDKU_AbsInt);
    REG_U(Re); symtab_set_ndarray_unary_kernel("Im", &NDKU_ImInt); REG_U(Arg);
    REG_U(Conjugate); REG_U(Sign);
    REG_U(UnitStep); REG_U(UnitBox); REG_U(Ramp);
    REG_U(Unitize); REG_B(Unitize, NDKB_Unitize);
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
    REG_U(EllipticK); REG_U(EllipticE);      /* complete: one argument  */
    REG_B(EllipticF, NDKB_EllipticF);        /* incomplete: phi and m   */
    REG_B(EllipticE, NDKB_EllipticE);
    REG_B(EllipticPi, NDKB_EllipticPi);
    REG_N(EllipticPi);
    REG_B(BesselY, NDKB_BesselY);
    REG_B(Beta,    NDKB_Beta);
    REG_B(StruveH, NDKB_StruveH);   /* real z via sf_machine_struve_h; Compile[] too */
    REG_B(StruveL, NDKB_StruveL);   /* real z via sf_machine_struve_l; Compile[] too */
    REG_B(AngerJ,  NDKB_AngerJ);    /* sf_machine_anger_weber */
    REG_B(WeberE,  NDKB_WeberE);
    REG_B(KelvinBer, NDKB_KelvinBer); REG_U(KelvinBer);   /* sf_machine_kelvin; */
    REG_B(KelvinBei, NDKB_KelvinBei); REG_U(KelvinBei);   /* REG_U = order-0   */
    REG_B(KelvinKer, NDKB_KelvinKer); REG_U(KelvinKer);   /* one-arg form      */
    REG_B(KelvinKei, NDKB_KelvinKei); REG_U(KelvinKei);
    REG_B(Surd,    NDKB_Surd);      /* declines even roots of x < 0 -> List path */
    REG_U(CubeRoot);                /* Surd[x, 3]: the same sf_surd_r */
    REG_U(RiemannR);                /* Gram series; declines x < 1/e -> List path */
    REG_B(BesselJZero, NDKB_BesselJZero); /* real n > -1, integer-valued k */
    REG_B(BesselYZero, NDKB_BesselYZero); /* real n >= 0, integer-valued k */
    REG_B(MathieuCharacteristicA, NDKB_MathieuCharacteristicA); /* sf_mathieu_char_r */
    REG_B(MathieuCharacteristicB, NDKB_MathieuCharacteristicB);

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
    REG_U(RiemannSiegelZ); REG_U(RiemannSiegelTheta);
    REG_B(PolyGamma,   NDKB_PolyGamma);
    REG_B(HurwitzZeta, NDKB_HurwitzZeta);
    REG_B(Zeta,        NDKB_Zeta);        /* Zeta[s, a] (array s or a) */
    REG_B(Pochhammer,  NDKB_Pochhammer);
    REG_B(Binomial,    NDKB_Binomial);
    REG_B(LegendreP,   NDKB_LegendreP);
    REG_B(LegendreQ,   NDKB_LegendreQ);
    REG_B(ExpIntegralE, NDKB_ExpIntegralE);
    REG_B(ChebyshevT,  NDKB_ChebyshevT);
    REG_B(ChebyshevU,  NDKB_ChebyshevU);
    REG_B(LaguerreL,   NDKB_LaguerreL);   /* LaguerreL[n, x] (array x) */
    REG_B(HermiteH,    NDKB_HermiteH);    /* HermiteH[n, x] (array x) */
    REG_B(ParabolicCylinderD, NDKB_ParabolicCylinderD); /* D_nu(z), real nu */
    REG_B(GegenbauerC, NDKB_GegenbauerC);   /* GegenbauerC[n, z] = 2 T_n(z)/n */
    REG_N(GegenbauerC);                     /* GegenbauerC[n, l, z] */
    REG_U(AiryAi); REG_U(AiryBi); REG_U(AiryAiPrime); REG_U(AiryBiPrime);
    REG_B(BesselI,           NDKB_BesselI);
    REG_B(BesselK,           NDKB_BesselK);
    REG_B(PolyLog,           NDKB_PolyLog);
    REG_B(QPochhammer,       NDKB_QPochhammer);
    REG_B(Hypergeometric0F1, NDKB_Hypergeometric0F1);
    REG_N(LerchPhi); REG_N(Hypergeometric1F1); REG_N(Hypergeometric2F1);
    REG_N(HypergeometricPFQ);
    REG_N(MeijerG);     /* MeijerG[{{..},{..}}, {{..},{..}}, z(, r)] for Compile[] */
    REG_N(LaguerreL);   /* generalized LaguerreL[n, a, x] for Compile[] */
    REG_N(JacobiP);     /* JacobiP[n, a, b, x]: buffers and Compile[] */
    REG_N(WhittakerM);  /* WhittakerM[k, m, z] for Compile[] */
    REG_N(MathieuC);       /* MathieuC[a, q, z] for Compile[] */
    REG_N(MathieuCPrime);  /* MathieuCPrime[a, q, z] for Compile[] */
    REG_N(MathieuS);       /* MathieuS[a, q, z] for Compile[] */
    REG_N(MathieuSPrime);  /* MathieuSPrime[a, q, z] for Compile[] */
    REG_N(HypergeometricU); /* HypergeometricU[a, b, z] for Compile[] */
    REG_N(WhittakerW);  /* WhittakerW[k, m, z] for Compile[] */
    REG_N(DirichletL);  /* DirichletL[k, j, s] for Compile[] */
    REG_N(PolyLog);     /* Nielsen PolyLog[n, p, z]; PolyLog[s, z] is REG_B */

}
