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
bool sf_machine_legendre_q(double n, double x, double* out);    /* LegendreQ, integer n, |x|<1 */
bool sf_machine_chebyshev_t(double n, double x, double* out);   /* ChebyshevT, integer n >= 0 */
/* ChebyshevT with a machine COMPLEX argument x = xre + i*xim (integer n >= 0).
 * Same contract/regime as the real helper; result in (*rr, *ri). */
bool sf_machine_chebyshev_t_c(double n, double xre, double xim,
                              double* rr, double* ri);
bool sf_machine_chebyshev_u(double n, double x, double* out);   /* ChebyshevU, integer n >= 0 */
/* ChebyshevU with a machine COMPLEX argument x = xre + i*xim (integer n >= 0).
 * Same contract/regime as the real helper; result in (*rr, *ri). */
bool sf_machine_chebyshev_u_c(double n, double xre, double xim,
                              double* rr, double* ri);
bool sf_machine_hermite_h(double n, double x, double* out);     /* HermiteH, integer n >= 0 */
/* HermiteH with a machine COMPLEX argument x = xre + i*xim (integer n >= 0).
 * Same contract/regime as the real helper; result in (*rr, *ri). */
bool sf_machine_hermite_h_c(double n, double xre, double xim,
                            double* rr, double* ri);
/* LaguerreL: the generalized Laguerre polynomial L_n^(a)(x) for a non-negative
 * integer degree n, by the stable forward recurrence (a = 0 is the ordinary
 * L_n(x)). Declines (returns false) for a non-integer/negative/large n or a
 * non-finite result, so the caller falls back to the MPFR / symbolic path. */
bool sf_machine_laguerre_l_assoc(double n, double a, double x, double* out);
bool sf_machine_laguerre_l(double n, double x, double* out);    /* a = 0 */
/* LaguerreL with a machine COMPLEX argument x = xre + i*xim (order n and
 * parameter a are real); result in (*rr, *ri). */
bool sf_machine_laguerre_l_assoc_c(double n, double a, double xre, double xim,
                                   double* rr, double* ri);
bool sf_machine_laguerre_l_c(double n, double xre, double xim,
                             double* rr, double* ri);
/* JacobiP: the Jacobi polynomial P_n^(a,b)(x) for a non-negative integer degree
 * n <= 1000 and real a, b, by the forward three-term recurrence. Declines for a
 * non-integer/negative/large n, a degenerate recurrence (a + b a negative
 * integer <= -2) or a non-finite result. */
bool sf_machine_jacobi_p(double n, double a, double b, double x, double* out);
/* JacobiP with a machine COMPLEX argument x = xre + i*xim (n, a, b real). */
bool sf_machine_jacobi_p_c(double n, double a, double b, double xre, double xim,
                           double* rr, double* ri);
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
/* WhittakerM M_{k,m}(z) = e^(-z/2) z^(m+1/2) 1F1(m-k+1/2; 2m+1; z), real path.
 * Declines (false) for z <= 0 (the z^(m+1/2) branch cut -> complex), the 2m+1
 * non-positive-integer pole, an un-summable 1F1, or a non-finite result, so the
 * caller falls to the MPFR / symbolic 1F1 template. */
bool sf_machine_whittaker_m(double k, double m, double z, double* out);     /* WhittakerM */
/* HypergeometricU U(a,b,z), the Tricomi confluent hypergeometric function, real
 * machine path for z > 0. Three regimes, each self-checking: a large-z asymptotic
 * series (which also captures the exactly-terminating cases), the integer-b
 * confluent-limit log-series (DLMF 13.2.9, the one regime where the connection
 * formula degenerates), and the non-integer-b Kummer connection to two 1F1's.
 * Declines (false) for z <= 0 (the z^(1-b) term leaves the real axis), a sum lost
 * to cancellation, a Gamma/digamma pole, or a non-finite result, so the caller
 * falls to the FLINT acb_hypgeom_u bridge / symbolic template. */
bool sf_machine_hyperu(double a, double b, double z, double* out);          /* HypergeometricU */
/* WhittakerW W_{k,m}(z) = e^(-z/2) z^(m+1/2) U(m-k+1/2, 2m+1, z), real path.
 * W is even in m, so m is folded to |m| (forcing b = 2|m|+1 >= 1). Same decline
 * contract as sf_machine_hyperu / sf_machine_whittaker_m; shared bit-for-bit with
 * the NDArray kernel, the Compile VM (ndk_WhittakerW_n) and the scalar fast path. */
bool sf_machine_whittaker_w(double k, double m, double z, double* out);     /* WhittakerW */
bool sf_machine_hyper2f1(double a, double b, double c, double z, double* out); /* 2F1 */
bool sf_machine_pfq(const double* a, size_t p, const double* b, size_t q,
                    double z, double* out);                    /* HypergeometricPFQ */
/* StruveH H_nu(z), real machine path. Two regimes, each self-checking: the
 * convergent 1F2-form series (small/moderate |z|, where the double cancellation
 * budget holds) and the Y_nu + K_nu asymptotic (large |z|). Declines (false) in
 * the uniform regime (moderate |z| with large |nu|), for a nu at/past the
 * 1/Gamma(nu+3/2) pole (nu <= -3/2 half/odd), a negative base with non-integer
 * nu+1 (complex), a series lost to cancellation, or a non-finite result, so the
 * caller falls to the FLINT arb bridge. Shared bit-for-bit with the NDArray
 * kernel (ndk_StruveH_c) and the Compile VM. */
bool sf_machine_struve_h(double nu, double z, double* out);    /* StruveH */
/* StruveL L_nu(z), real machine path: the cancellation-free convergent series
 * (DLMF 11.2.2) and, for z > 30, I_{-nu} + the DLMF 11.6.2 asymptotic. Declines
 * where the value is complex (z < 0, non-integer nu), on cancellation
 * (nu < -3/2), Gamma overflow or a non-finite result. Shared bit for bit by the
 * scalar builtin, the NDArray kernel (ndk_StruveL_c) and the Compile VM. */
bool sf_machine_struve_l(double nu, double z, double* out);    /* StruveL */
/* AngerJ J_nu(z) (which == 0) / WeberE E_nu(z) (which == 1), real machine path
 * (sf_machine_angerkelvin.c): the S1/S2 power series (DLMF 11.10.8-9) for small
 * |z|, else J_nu / -Y_nu plus the Gauss-Legendre-integrated Laplace-type
 * integrals of DLMF 11.10.4-11.10.6; negative z through the parity
 * J_nu(-z) = J_{-nu}(z), E_nu(-z) = -E_{-nu}(z). Declines when neither regime
 * meets its error budget (the arb bridge answers). */
bool sf_machine_anger_weber(int which, double nu, double z, double* out);
/* Kelvin functions, real machine path (sf_machine_angerkelvin.c), which =
 * 0 ber, 1 bei, 2 ker, 3 kei. ber+i bei = J_nu(x e^(3 pi i/4)) and
 * ker+i kei = e^(-nu pi i/2) K_nu(x e^(pi i/4)) through the shared complex
 * Bessel kernels. Declines where the value is complex (x < 0 except ber/bei at
 * an integer order) or the kernel declines. */
bool sf_machine_kelvin(int which, double nu, double x, double* out);

/* Already-existing double algorithms in inverf.c / inverfc.c, exposed so the
 * same two consumers can use them rather than degrading. */
double inverf_double(double x);
double inverfc_double(double x);

/* Cancellation budget for the ascending series in the exponential-integral
 * family.  See sf_machine.c for why this is measured, not derived from |z|. */
bool sf_series_usable(double peak_term, double result_magnitude);

#endif /* MATHILDA_SF_MACHINE_H */
