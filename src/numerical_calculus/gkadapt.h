/*
 * gkadapt.h — globally-adaptive Gauss-Kronrod quadrature (machine precision)
 *
 * Evaluates  I = ∫_a^b f(x) dx  for a (possibly complex-valued) integrand on a
 * finite real interval, using the 7-15 Gauss-Kronrod pair (QUADPACK QK15) as
 * the local rule and a global-adaptive bisection strategy: the subinterval
 * carrying the largest local error estimate is repeatedly halved until the
 * accumulated error meets the requested tolerance or a resource limit is hit.
 *
 * The local error estimate is QUADPACK's: |K - G| rescaled by the rule's
 * smoothness gauge resasc, with a roundoff floor from resabs.  For a
 * complex-valued integrand the partial results are complex and the magnitudes
 * (|f|, |K-G|, …) are taken with cabs, so the same heuristic applies to real
 * and complex integrands alike.
 *
 * Optionally (extrapolate = true) Wynn's epsilon-algorithm is run on the
 * sequence of running integral estimates produced as the worst interval is
 * peeled away — the QAGS device that accelerates convergence when the
 * integrand has an integrable endpoint singularity.  The accelerated value is
 * returned when it both improves the error estimate and is consistent with the
 * raw running sum.
 *
 * The integrand is supplied through a sample callback returning machine
 * `double _Complex`, mirroring quadrature.h / dequad.h.  A sample that returns
 * false (non-numeric / non-finite) marks the panel non-evaluable.
 */
#ifndef MATHILDA_GKADAPT_H
#define MATHILDA_GKADAPT_H

#include <complex.h>
#include <stdbool.h>

/* Machine sample: write f(x) into *out; return false if not a finite number. */
typedef bool (*GkSampleMachine)(void* ctx, double x, double _Complex* out);

typedef enum {
    GK_OK = 0,        /* converged to the requested tolerance                 */
    GK_NOCONV,        /* resource limit reached; best estimate still returned  */
    GK_NONNUMERIC     /* integrand not numeric on the interval (no estimate)   */
} GkStatus;

typedef struct {
    double _Complex value;    /* integral estimate                            */
    double          abs_err;  /* estimated absolute error                     */
    GkStatus        status;
    long            n_eval;    /* integrand evaluations performed             */
    int             n_subdiv;  /* bisections performed                        */
} GkResult;

/* Compute ∫_a^b f dx with global-adaptive Gauss-Kronrod.
 *   abstol, reltol  termination: abs_err <= max(abstol, reltol*|I|).
 *   max_subdiv      cap on the number of bisections.
 *   max_eval        cap on integrand evaluations (<=0 => derive from max_subdiv).
 *   extrapolate     enable Wynn epsilon-extrapolation of the running estimate.
 */
GkResult gk_integrate_machine(GkSampleMachine f, void* ctx,
                              double a, double b,
                              double abstol, double reltol,
                              int max_subdiv, long max_eval,
                              bool extrapolate);

#endif /* MATHILDA_GKADAPT_H */
