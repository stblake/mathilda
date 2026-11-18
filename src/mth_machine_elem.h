/*
 * mth_machine_elem.h -- the ONE machine-precision evaluation of each
 * elementary function: Exp, the circular and hyperbolic functions and their
 * inverses (Log has its own header, mth_machine_log.h, on the same rule).
 *
 * Every surface that evaluates these at machine precision -- the scalar
 * builtins (src/trig.c, src/hyperbolic.c), the packed / NDArray kernels
 * (src/ndkernels.c) and the Compile[] VM (src/compile/compile_vm.c) -- goes
 * through these helpers, so the answers agree bit for bit by construction.
 *
 * The rule is the one mth_machine_log.h set for Log: a REAL argument inside the
 * function's real domain takes the real libm function; anything else (a
 * complex argument, or a real one off the domain, where the value is complex)
 * takes the C99 complex function, with the branch-cut side the scalar builtins
 * have always used. C99 leaves the accuracy of both families to the
 * implementation, and they round differently: on Apple's libm
 * creal(catan(0.5)) is 0.46364760900080609 where atan(0.5) is
 * 0.46364760900080615, and creal(casin(0.5)) is 0.52359877559829893 where
 * asin(0.5) is 0.52359877559829882. The scalar builtins took the complex
 * function for real input, the kernels' real arms and the VM the real one, so
 * ArcTan[0.5] differed by one ulp between the interpreter and a packed array
 * (numeric-tower review G4). The forward functions had the same split:
 * creal(csinh(-0.7)) is -0.75858370183953339 where sinh(-0.7) is
 * -0.7585837018395335, and 1/ctanh(2.5) is 1.0135673098126086 where
 * 1/tanh(2.5) is 1.0135673098126083. Neither family is uniformly the closer one; the
 * real one is the faster, and is what a real argument means.
 *
 * Signed zero is IEEE on every surface (the G7 decision): ArcCot[-0.] is
 * atan(1/-0.) = -Pi/2.
 */
#ifndef MTH_MACHINE_ELEM_H
#define MTH_MACHINE_ELEM_H

#include <math.h>
#include <complex.h>

/* The complex number re + im I with BOTH parts exactly as given, signed zeros
 * included. `re + im * I` cannot do that in C99: with I a complex constant the
 * real part is computed as re + 0., and -0. + 0. is +0., so a -0. argument
 * reached every complex path as +0. (ArcCot[-0.] was Pi/2 on the scalar
 * surface, -Pi/2 on the kernels': numeric-tower review G7, IEEE signed zero on
 * every surface). A double complex has the representation of double[2]
 * (C99 6.2.5p13), so the parts are written directly. */
static inline double complex mth_cmplx(double re, double im) {
    union { double complex z; double p[2]; } u;
    u.p[0] = re; u.p[1] = im;
    return u.z;
}
static inline double complex mth_inv_real(double v) { return mth_cmplx(v, 0.0); }

/* ArcSin: real on [-1, 1]. Off it, C99 casin(x + 0 I) for x > 1 lands on the
 * upper side of the cut; Mathilda takes the lower side (the x < -1 cut agrees). */
static inline double complex mth_asin_machine_z(double complex z) {
    double x = creal(z), y = cimag(z);
    if (y == 0.0 && x >= -1.0 && x <= 1.0) return mth_inv_real(asin(x));
    double complex w = casin(z);
    if (y == 0.0 && x > 1.0) w = creal(w) - cimag(w) * I;
    return w;
}

/* ArcCos: real on [-1, 1]; for real x > 1 the +i acosh(x) side. */
static inline double complex mth_acos_machine_z(double complex z) {
    double x = creal(z), y = cimag(z);
    if (y == 0.0 && x >= -1.0 && x <= 1.0) return mth_inv_real(acos(x));
    double complex w = cacos(z);
    if (y == 0.0 && x > 1.0) w = creal(w) - cimag(w) * I;
    return w;
}

/* ArcTan: real on the whole real line. */
static inline double complex mth_atan_machine_z(double complex z) {
    if (cimag(z) == 0.0) return mth_inv_real(atan(creal(z)));
    return catan(z);
}

/* ArcCot[z] = ArcTan[1/z]: real on the real line, +-Pi/2 at +-0. */
static inline double complex mth_acot_machine_z(double complex z) {
    if (cimag(z) == 0.0) return mth_inv_real(atan(1.0 / creal(z)));
    return catan(1.0 / z);
}

/* ArcSinh: real on the real line. */
static inline double complex mth_asinh_machine_z(double complex z) {
    if (cimag(z) == 0.0) return mth_inv_real(asinh(creal(z)));
    return casinh(z);
}

/* ArcCsch[z] = ArcSinh[1/z]: real on the real line away from 0. */
static inline double complex mth_acsch_machine_z(double complex z) {
    if (cimag(z) == 0.0) return mth_inv_real(asinh(1.0 / creal(z)));
    return casinh(1.0 / z);
}

/* ---- forward functions: real on the whole real line ---------------------- *
 * The real arm is the expression the packed kernel's real arm uses (ndkernels.c
 * UK_CLOSED), so a real lane is bit-identical on every surface. */
#define MTH_ELEM_FWD(NAME, REXPR, CEXPR)                                       \
    static inline double complex mth_##NAME##_machine_z(double complex z) {    \
        if (cimag(z) == 0.0) { double x = creal(z); return mth_inv_real(REXPR); } \
        return (CEXPR);                                                        \
    }
MTH_ELEM_FWD(exp,  exp(x),        cexp(z))
MTH_ELEM_FWD(sin,  sin(x),        csin(z))
MTH_ELEM_FWD(cos,  cos(x),        ccos(z))
MTH_ELEM_FWD(tan,  tan(x),        ctan(z))
MTH_ELEM_FWD(cot,  1.0 / tan(x),  1.0 / ctan(z))
MTH_ELEM_FWD(sec,  1.0 / cos(x),  1.0 / ccos(z))
MTH_ELEM_FWD(csc,  1.0 / sin(x),  1.0 / csin(z))
MTH_ELEM_FWD(sinh, sinh(x),       csinh(z))
MTH_ELEM_FWD(cosh, cosh(x),       ccosh(z))
MTH_ELEM_FWD(tanh, tanh(x),       ctanh(z))
MTH_ELEM_FWD(coth, 1.0 / tanh(x), 1.0 / ctanh(z))
MTH_ELEM_FWD(sech, 1.0 / cosh(x), 1.0 / ccosh(z))
MTH_ELEM_FWD(csch, 1.0 / sinh(x), 1.0 / csinh(z))
#undef MTH_ELEM_FWD

#endif /* MTH_MACHINE_ELEM_H */
