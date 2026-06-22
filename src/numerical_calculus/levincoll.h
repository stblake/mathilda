/*
 * levincoll.h — Levin collocation rule for highly oscillatory integrals
 *
 * Levin's method computes  I = ∫_a^b f(x) e^{i g(x)} dx  (f slowly varying, the
 * kernel e^{ig} oscillating rapidly) WITHOUT resolving each oscillation.  It
 * seeks an antiderivative-like function p with  d/dx[p e^{ig}] = f e^{ig}, i.e.
 * the collocation ODE
 *
 *     p'(x) + i g'(x) p(x) = f(x).
 *
 * Approximating p in a Chebyshev basis and collocating at Chebyshev–Gauss–
 * Lobatto nodes turns the ODE into a small dense linear system; the integral is
 * then the boundary term  p(b) e^{i g(b)} − p(a) e^{i g(a)}.  Accuracy IMPROVES
 * with the oscillation rate, the opposite of ordinary quadrature.
 *
 * Real trigonometric kernels are handled as the real/imaginary part of the
 * complex-exponential solve:  ∫ f cos g = Re I,  ∫ f sin g = Im I.
 *
 * This module is pure-numeric and sampler-driven: the amplitude f, the phase g
 * and its derivative g' are supplied as callbacks (the caller — nint.c — owns
 * the kernel detection, the symbolic g' = D[g,x], and the variable binding), so
 * the engine is decoupled from the expression layer and unit-testable on its own.
 */
#ifndef MATHILDA_LEVINCOLL_H
#define MATHILDA_LEVINCOLL_H

#include <complex.h>
#include <stdbool.h>

#include "gkadapt.h"   /* GkSampleMachine */

/* Oscillatory kernel form. */
typedef enum {
    LEVIN_KERNEL_EXP = 0,  /* f · e^{i g}                                       */
    LEVIN_KERNEL_COS = 1,  /* f · cos g  (= Re of the EXP solve)                */
    LEVIN_KERNEL_SIN = 2   /* f · sin g  (= Im of the EXP solve)               */
} LevinKernel;

typedef struct {
    bool           have;   /* an estimate was produced at all                  */
    bool           conv;   /* it settled to the requested relative tolerance   */
    double _Complex val;   /* the (kernel-projected) integral estimate         */
    double         err;    /* last successive-order difference magnitude       */
} LevinResult;

/* Levin collocation for  ∫_a^b f · {e^{ig}|cos g|sin g} dx  on a real axis.
 *
 *   amp     — samples f(x)            (slowly varying amplitude)
 *   gprime  — samples g'(x)           (phase derivative; the oscillation rate)
 *   gphase  — samples g(x)            (used only at the two endpoints)
 *
 * Each *_ctx is the opaque context handed to the matching callback.  The order
 * n is doubled (n0, 2·n0, …) until two successive estimates agree to `reltol`
 * or n exceeds `n_max`.  A near-singular collocation matrix (weak oscillation /
 * a stationary phase point) is reported as have=false so the caller can fall
 * back to ordinary quadrature.  Returns have=false when LAPACK is unavailable. */
LevinResult levin_collocation_machine(
    double a, double b,
    GkSampleMachine amp,    void* amp_ctx,
    GkSampleMachine gprime, void* gprime_ctx,
    GkSampleMachine gphase, void* gphase_ctx,
    LevinKernel kernel, double reltol, int n_max);

/* Prepared (factored-once) Levin solver.  The collocation matrix depends only on
 * the nodes and the phase derivative g', NOT on the amplitude f — so when many
 * integrals share the same interval and g' (e.g. the inner axis of a separable
 * multivariate Levin reduction, where g' is constant in the outer variable) the
 * matrix can be factored once and re-used across right-hand sides.  This turns a
 * per-sample O(n^3) factorisation (plus condition estimate) into an O(n^2)
 * back-substitution, which is the difference between a tractable and an
 * intractable multivariate reduction.
 *
 * levin_prepare_machine factors the order-n system; returns NULL on a singular
 * matrix or when LAPACK is unavailable.  levin_prepared_solve evaluates one
 * integral for a given amplitude/phase; the phase g is sampled only at the
 * endpoints, so g MAY depend on parameters even though g' must not. */
typedef struct LevinPrep LevinPrep;

LevinPrep* levin_prepare_machine(double a, double b, int n,
                                 GkSampleMachine gprime, void* gprime_ctx);
bool levin_prepared_solve(const LevinPrep* prep,
                          GkSampleMachine amp,    void* amp_ctx,
                          GkSampleMachine gphase, void* gphase_ctx,
                          LevinKernel kernel, double _Complex* out);
void levin_prepare_free(LevinPrep* prep);

#ifdef USE_MPFR
#include <mpfr.h>

/* Arbitrary-precision sampler: write f(x) into (out_re, out_im) at the
 * precision of those (pre-initialised) outputs; return false if non-finite. */
typedef bool (*LevinSampleMPFR)(void* ctx, const mpfr_t x,
                                mpfr_t out_re, mpfr_t out_im);

/* Arbitrary-precision Levin collocation.  The collocation nodes, Chebyshev
 * basis, complex matrix and its LU solve, and the boundary term are all carried
 * at `bits` of precision.  Writes the (kernel-projected) result into
 * (out_re, out_im) and the convergence verdict into *converged; returns false
 * if no estimate could be formed (detection/sampling failure, singular system,
 * or LAPACK-independent — this path uses an in-house ncpx solve). */
bool levin_collocation_mpfr(double a, double b, long bits,
                            LevinSampleMPFR amp,    void* amp_ctx,
                            LevinSampleMPFR gprime, void* gprime_ctx,
                            LevinSampleMPFR gphase, void* gphase_ctx,
                            LevinKernel kernel, double reltol, int n_max,
                            mpfr_t out_re, mpfr_t out_im, bool* converged);
#endif

#endif /* MATHILDA_LEVINCOLL_H */
