/* Machine-precision (double) kernels for the special functions.
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
#ifndef MATHILDA_SF_MACHINE_H
#define MATHILDA_SF_MACHINE_H

#include <stdbool.h>
#include <stddef.h>

bool sf_machine_ei(double x, double* out);    /* ExpIntegralEi */
bool sf_machine_e1(double x, double* out);    /* E_1, x > 0 */
bool sf_machine_si(double x, double* out);    /* SinIntegral */
bool sf_machine_ci(double x, double* out);    /* CosIntegral,  x > 0 */
bool sf_machine_shi(double x, double* out);   /* SinhIntegral */
bool sf_machine_chi(double x, double* out);   /* CoshIntegral, x > 0 */
bool sf_machine_li(double x, double* out);    /* LogIntegral,  x >= 0, x != 1 */
bool sf_machine_sinc(double x, double* out);  /* Sinc */
bool sf_machine_erfi(double x, double* out);        /* Erfi */
bool sf_machine_productlog(double x, double* out);  /* ProductLog, x >= -1/e */
bool sf_machine_fresnel_c(double x, double* out);   /* FresnelC */
bool sf_machine_fresnel_s(double x, double* out);   /* FresnelS */
bool sf_machine_digamma(double x, double* out);     /* PolyGamma[x] */
bool sf_machine_harmonic(double x, double* out);    /* HarmonicNumber */
bool sf_machine_zeta(double s, double* out);        /* Zeta, s != 1 */
bool sf_machine_fibonacci(double x, double* out);   /* Fibonacci */
bool sf_machine_lucasl(double x, double* out);      /* LucasL */
bool sf_machine_hurwitz_zeta(double s, double a, double* out);  /* s > 1, a > 0 */
bool sf_machine_polygamma(double n, double x, double* out);     /* PolyGamma[n, x] */
bool sf_machine_pochhammer(double a, double n, double* out);    /* Pochhammer */
bool sf_machine_binomial(double n, double k, double* out);      /* Binomial */
bool sf_machine_legendre_p(double n, double x, double* out);    /* LegendreP, integer n */
bool sf_machine_airy_ai(double x, double* out);        /* AiryAi */
bool sf_machine_airy_bi(double x, double* out);        /* AiryBi */
bool sf_machine_airy_ai_prime(double x, double* out);  /* AiryAiPrime */
bool sf_machine_airy_bi_prime(double x, double* out);  /* AiryBiPrime */
bool sf_machine_bessel_i(double nu, double x, double* out);     /* BesselI */
bool sf_machine_hyper0f1(double a, double z, double* out);      /* Hypergeometric0F1 */
bool sf_machine_qpochhammer(double a, double q, double* out);   /* QPochhammer, |q|<1 */
bool sf_machine_bessel_k(double nu, double x, double* out);     /* BesselK, x > 0 */
bool sf_machine_polylog(double s, double x, double* out);       /* PolyLog, |x| <= 1 */
bool sf_machine_lerchphi(double z, double s, double a, double* out);        /* LerchPhi */
bool sf_machine_hyper1f1(double a, double b, double z, double* out);        /* 1F1 */
bool sf_machine_hyper2f1(double a, double b, double c, double z, double* out); /* 2F1 */
bool sf_machine_pfq(const double* a, size_t p, const double* b, size_t q,
                    double z, double* out);                    /* HypergeometricPFQ */

/* Already-existing double algorithms in inverf.c / inverfc.c, exposed so the
 * same two consumers can use them rather than degrading. */
double inverf_double(double x);
double inverfc_double(double x);

/* Cancellation budget for the ascending series in the exponential-integral
 * family.  See sf_machine.c for why this is measured, not derived from |z|. */
bool sf_series_usable(double peak_term, double result_magnitude);

#endif /* MATHILDA_SF_MACHINE_H */
