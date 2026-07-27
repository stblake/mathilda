/* Machine-precision (double) kernels for the exponential-integral family.
 *
 * The modules that own these functions (expintegralei.c, sinintegral.c, ...)
 * compute everything in MPFR and round at the end, which is right for
 * arbitrary-precision work but leaves no double-only path for the two consumers
 * that need one: the NDArray element-wise fast path, and the Compile[] VM.
 * Both were therefore stuck on the "degrade sentinel" — correct, but an order of
 * magnitude slower, and in the compiler's case enough to make the WHOLE
 * surrounding body fall back to the interpreter.
 *
 * Every entry point has the shape the shared kernel registry expects,
 * `bool f(double in, double* out)`, and returns FALSE for an argument whose true
 * value is complex or infinite — Ci/Chi/Li of a negative argument, Ei at zero.
 * That is the engine's established contract: a real-typed fast path declines
 * exactly where the interpreter would leave the real axis, and the caller falls
 * back rather than inventing a real answer.
 */
#ifndef MATHILDA_EXPINT_MACHINE_H
#define MATHILDA_EXPINT_MACHINE_H

#include <stdbool.h>

bool sf_machine_ei(double x, double* out);    /* ExpIntegralEi */
bool sf_machine_e1(double x, double* out);    /* E_1, x > 0 */
bool sf_machine_si(double x, double* out);    /* SinIntegral */
bool sf_machine_ci(double x, double* out);    /* CosIntegral,  x > 0 */
bool sf_machine_shi(double x, double* out);   /* SinhIntegral */
bool sf_machine_chi(double x, double* out);   /* CoshIntegral, x > 0 */
bool sf_machine_li(double x, double* out);    /* LogIntegral,  x >= 0, x != 1 */
bool sf_machine_sinc(double x, double* out);  /* Sinc */

/* Already-existing double algorithms in inverf.c / inverfc.c, exposed so the
 * same two consumers can use them rather than degrading. */
double inverf_double(double x);
double inverfc_double(double x);

#endif /* MATHILDA_EXPINT_MACHINE_H */
