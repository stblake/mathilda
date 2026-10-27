/*
 * flint_num_bridge.h
 * ------------------
 * Boundary between Mathilda's numeric Expr scalars and FLINT's arbitrary-
 * precision complex ball arithmetic (arb/acb), for the transcendental special
 * functions. FLINT's acb_hypgeom / acb_dirichlet kernels are rigorous (correct
 * rounding, full complex plane, standard branch cuts), so they back the FLINT`
 * context wrappers exposed at the REPL and, over time, transparent numeric
 * fast paths for the hand-rolled MPFR kernels.
 *
 * Each entry point takes numeric Expr arguments (Integer / BigInt / Rational /
 * Real / MPFR / Complex of those) and returns a freshly-owned numeric Expr
 * (EXPR_MPFR, or Complex[...] when the imaginary part is non-zero), using the
 * working precision implied by the arguments (machine 53-bit floor). It returns
 * NULL — leaving the caller's FLINT`f[...] unevaluated — when an argument is
 * non-numeric / symbolic, or the result is a pole / non-finite, or FLINT is
 * absent (USE_FLINT undefined). No argument is mutated.
 */
#ifndef FLINT_NUM_BRIDGE_H
#define FLINT_NUM_BRIDGE_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

Expr* flint_num_zeta(const Expr* s);                       /* Zeta[s]            */
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a);/* HurwitzZeta[s, a]  */
Expr* flint_num_polygamma(const Expr* n, const Expr* z);   /* PolyGamma[n, z]    */
Expr* flint_num_polygamma_frac(const Expr* nu, const Expr* z); /* non-integer order */
Expr* flint_num_polylog(const Expr* s, const Expr* z);     /* PolyLog[s, z]      */
Expr* flint_num_airy_ai(const Expr* z);                    /* AiryAi[z]          */
Expr* flint_num_erf(const Expr* z);                        /* Erf[z]  (entire)   */
Expr* flint_num_erfc(const Expr* z);                       /* Erfc[z] (entire)   */
Expr* flint_num_barnesg(const Expr* z);                    /* BarnesG[z] (entire) */
Expr* flint_num_erfi(const Expr* z);                       /* Erfi[z] (entire)   */
Expr* flint_num_airy_ai_prime(const Expr* z);              /* AiryAiPrime[z]     */
Expr* flint_num_airy_bi(const Expr* z);                    /* AiryBi[z]          */
Expr* flint_num_airy_bi_prime(const Expr* z);              /* AiryBiPrime[z]     */
Expr* flint_num_lerchphi(const Expr* z, const Expr* s, const Expr* a); /* LerchPhi[z,s,a] */
Expr* flint_num_stieltjes(const Expr* n, const Expr* a);   /* StieltjesGamma[n,a]*/
Expr* flint_num_expint(const Expr* s, const Expr* z);      /* ExpIntegralE[s, z] */
/* PolyLog[n, p, z], Nielsen's S_{n,p}(z), numerically: the Stirling-number
 * series (positive integer p, small |z|), the generalized-power series
 * (any p, |z| < 1) or tanh-sinh quadrature of
 * 1/(Gamma(n) Gamma(p+1)) Int_0^1 (-Log t)^(n-1) (-Log(1 - z t))^p dt/t (Re n > 0,
 * any z; principal logarithms). The caller decides which calls evaluate. NULL
 * when no method applies or the value is not finite. */
Expr* flint_num_nielsen_polylog(const Expr* n, const Expr* p, const Expr* z);
/* BesselJZero at arbitrary precision: the zero of BesselJ[nu, x] nearest the
 * machine zero `seed`, by Newton in arb, returned at `bits` bits. */
Expr* flint_num_bessel_j_zero(const Expr* nu, double seed, long bits);
/* QPochhammer[a, q] (n == NULL) and QPochhammer[a, q, n] for any n, |q| < 1:
 * (a; q)_n = (a; q)_inf / (a q^n; q)_inf. */
Expr* flint_num_qpochhammer(const Expr* a, const Expr* q, const Expr* n);

/* Legendre functions of a general degree n and order m. `type` follows arb:
 * 0 is the Ferrers form (the reference's types 1 and 2, cut on (-inf, -1] and
 * [1, inf)), 1 the form analytic off (-inf, 1] (the reference's type 3). */
Expr* flint_num_legendre_p(const Expr* n, const Expr* m, const Expr* z, int type);
Expr* flint_num_legendre_q(const Expr* n, const Expr* m, const Expr* z, int type);

/* The numeric tail of a special-function builtin. flint_num_call_inexact is
 * true iff every argument of `call` is a number (Integer, Rational, Real,
 * arbitrary-precision Real, or a Complex of those) and at least one of them
 * is inexact: the case where the arb value is the answer. flint_num_for_call
 * takes an arb result (consumed; NULL passes through) and lowers it to
 * machine Reals when a machine-precision argument governs, so a machine call
 * answers with a machine number. Both always compiled. */
bool  flint_num_call_inexact(const Expr* call);
Expr* flint_num_for_call(Expr* r, const Expr* call);

/* The Legendre elliptic integrals, in the PARAMETER convention m = k^2 that
 * Mathematica uses (not the modulus k) -- which is also Arb's, so no argument
 * translation happens here:
 *
 *   EllipticK[m]           = F(Pi/2 | m)
 *   EllipticE[m]           = E(Pi/2 | m)
 *   EllipticF[phi, m]      = Int_0^phi dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticE[phi, m]      = Int_0^phi Sqrt(1 - m Sin[t]^2) dt
 *   EllipticPi[n, m]       = Pi(n; Pi/2 | m)
 *   EllipticPi[n, phi, m]  = Int_0^phi dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *
 * Arb defines the incomplete forms on -Pi/2 <= Re phi <= Pi/2 and extends them
 * by the quasi-period (F(phi + k Pi | m) = F(phi | m) + 2 k K(m), and likewise
 * for E and Pi), so a phi of any size or complexity is handled here rather than
 * by the caller. EllipticPi with n > 1 crosses the pole at Sin[t]^2 = 1/n and is
 * the Cauchy principal value, which acb_elliptic_pi_inc computes -- a complex
 * value on the principal branch, matching Mathematica. */
Expr* flint_num_elliptic_k(const Expr* m);                              /* EllipticK[m]        */
Expr* flint_num_elliptic_e(const Expr* m);                              /* EllipticE[m]        */
Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m);             /* EllipticF[phi, m]   */
Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m);         /* EllipticE[phi, m]   */
Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m);              /* EllipticPi[n, m]    */
Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi,
                                const Expr* m);                          /* EllipticPi[n,phi,m] */

/* HypergeometricU[a, b, z] -- the Tricomi confluent hypergeometric U, via FLINT's
 * rigorous acb_hypgeom_u (correct rounding, full complex plane, standard branch
 * cut on z from -infinity to 0, integer b handled natively). Used for the MPFR /
 * complex / z<=0 numerics of HypergeometricU and (through the confluent template)
 * WhittakerW. */
Expr* flint_num_hypergeometric_u(const Expr* a, const Expr* b, const Expr* z);

/* StruveH[nu, z] -- the Struve function H_nu(z). FLINT/arb ships no Struve
 * kernel, so this is composed from the rigorous regularized generalized
 * hypergeometric: H_nu(z) = (z/2)^(nu+1) * 1F2~(1; 3/2, nu+3/2; -z^2/4). The
 * regularized form folds the 1/Gamma factors into the sum, so it is entire in
 * nu (no poles) and valid over the whole complex plane. A cancellation
 * allowance (~1.44*|z| bits) is pre-added to the working precision so a large
 * |z| gets its bits up front rather than by blind doubling. */
Expr* flint_num_struveh(const Expr* nu, const Expr* z);     /* StruveH[nu, z]    */

/* Hypergeometric2F1[a, b, c, z] at arbitrary precision (acb_hypgeom_2f1, with
 * its transformations near z = 1). Used by HypergeometricPFQ's MPFR path. */
Expr* flint_num_hyp2f1(const Expr* a, const Expr* b, const Expr* c, const Expr* z);

/* MeijerG[{{a_1..a_n},{a_n+1..a_p}}, {{b_1..b_m},{b_m+1..b_q}}, z, r] via
 * Slater's theorem over arb's regularized pFq (the dual expansion in 1/w for
 * p > q or p == q, |w| > 1; w = z^(1/r)), with the log case (integer-differing
 * b's) handled by perturbation. `a` holds a_1..a_p, `b` b_1..b_q; `r` may be
 * NULL (r = 1). The output precision is the minimum inexact precision of the
 * arguments (machine floor). NULL for non-numeric input, z = 0, a value that
 * does not exist, or p == q >= 3 with |w| = 1 (no continuation implemented). */
Expr* flint_num_meijerg(const Expr* const* a, size_t n, size_t p,
                        const Expr* const* b, size_t m, size_t q,
                        const Expr* z, const Expr* r);

/* The same engine on machine inputs, rounded to double: the fallback behind
 * src/special_functions/sf_machine_meijerg.c when the machine Slater sum
 * declines (log case, cancellation), so every machine surface still answers.
 * ai / bi may be NULL. *is_c is false when the imaginary part is exactly zero.
 * Returns false without FLINT or when the value does not exist. */
bool flint_mg_eval_d(const double* ar, const double* ai, size_t n, size_t p,
                     const double* br, const double* bi, size_t m, size_t q,
                     double r, double zr, double zi,
                     double* rr, double* ri, bool* is_c);

/* Laurent expansion about eps = 0 of
 *
 *   F(eps) = prod_{i < nf} Gamma(x_i + s_i eps)^{e_i},   s_i, e_i in {+1, -1},
 *
 * as F = eps^v (f_0 + f_1 eps + ... + f_{T-1} eps^{T-1} + O(eps^T)). x_i is
 * read EXACTLY from x[4i .. 4i+3] = (re.hi, re.lo, im.hi, im.lo); an x_i with
 * lo parts and imaginary part zero and an integer hi <= 0 is a pole of Gamma
 * (a zero of 1/Gamma) and moves v by -e_i. The coefficients come back as
 * double-doubles, out[4k .. 4k+3] = (re.hi, re.lo, im.hi, im.lo) of f_k,
 * computed in arb to >= 110 certified bits. This is the z-independent start of
 * sf_machine_meijerg.c's log-case residue table. False without FLINT, or when
 * the accuracy is not reached (the caller then has no table and declines). */
bool flint_mg_laurent(size_t nf, const double* x, const int* s, const int* e,
                      size_t T, long* v, double* out);

/* Dirichlet L-function L(s, chi) for a character given as a PRIMITIVE
 * character chi* mod q (conductor q) plus Euler factors: the value is
 * L(s, chi*) * prod_i (1 - chi*(eul_p[i]) eul_p[i]^-s). chi* is described by a
 * callback: fn(ctx, n, &num) sets chi*(n) = exp(2 Pi I num/den) and returns
 * false when gcd(n, q) > 1; chi*(eul_p[i]) = exp(2 Pi I eul_num[i]/den).
 * Precision follows s (floored at 53 bits), with the usual doubling ladder.
 * NULL at a pole or without FLINT. */
typedef bool (*flint_chi_exp_fn)(const void* ctx, uint64_t n, uint64_t* num);
Expr* flint_num_dirichlet_l(const Expr* s, uint64_t q, uint64_t den,
                            flint_chi_exp_fn fn, const void* ctx,
                            int neul, const uint64_t* eul_p, const uint64_t* eul_num);

/* Rigorous real enclosure [*out_lo, *out_hi] of L(s, chi) for s in [lo, hi]
 * (a REAL character, so the value is real). Monotone pieces are certified by
 * the sign of L'(s) on the ball and bounded by their endpoint values; the rest
 * is bisected and bounded by value balls. Outward rounded. False without FLINT
 * or when the range meets a pole. */
bool flint_dirichlet_l_interval(double lo, double hi, uint64_t q, uint64_t den,
                                flint_chi_exp_fn fn, const void* ctx,
                                int neul, const uint64_t* eul_p, const uint64_t* eul_num,
                                double* out_lo, double* out_hi);

/* Certified signs of the real and imaginary parts of an EXACT numeric
 * constant, by rigorous acb ball evaluation (precision ladder 64 -> 256 ->
 * 1024 bits). This is the sign oracle behind Abs / Sign / Re / Im / ReIm /
 * Arg of closed-form constants (Abs[Pi - 4] -> 4 - Pi).
 *
 * Supported language: Integer / BigInt / Rational / Complex of those; the
 * constants Pi, E, EulerGamma, Catalan, GoldenRatio, Degree, Khinchin,
 * Glaisher; Plus, Times, Power, Sqrt, Exp, Log (1 and 2 args), the six trig
 * and six hyperbolic functions, ArcSin/ArcCos/ArcTan/ArcSinh/ArcCosh/ArcTanh,
 * Gamma, Erf, Erfc, Abs, Re, Im, Conjugate. Multi-valued functions are
 * evaluated only where the input ball provably avoids the branch cut (or lies
 * exactly on the negative real axis for Log/Power, where the principal value
 * is unambiguous), so the ball always encloses the value Mathematica means.
 *
 * Returns false when `e` lies outside that language (a free symbol, an
 * inexact number, an unsupported head) or does not evaluate to a finite ball.
 * Otherwise sets *re_sgn and *im_sgn each to -1 / 0 / +1 when certified, or 2
 * when the ball still contains zero at the top precision. 0 is reported only
 * for a ball that is EXACTLY zero (exact cancellation, e.g. the imaginary part
 * of a real expression); an expression that is zero only symbolically, e.g.
 * (1 + Sqrt[2])^2 - 3 - 2 Sqrt[2], reports 2 (undecided) and callers decline.
 * Read-only; false without FLINT. */
bool flint_cert_parts_sign(const Expr* e, int* re_sgn, int* im_sgn);
/* Drop flint_cert_eval_real's variable binding.  It never calls the evaluator,
 * so no TimeConstrained body starts inside it: a timeout that jumps out of it
 * leaves a binding nothing owns, which tc_run_guarded clears. */
void flint_cert_binding_clear(void);

/* Registers the FLINT` context numeric builtins (Zeta / HurwitzZeta /
 * PolyGamma / StieltjesGamma). Called from core_init(). No-op without FLINT. */
void flint_num_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_NUM_BRIDGE_H */
