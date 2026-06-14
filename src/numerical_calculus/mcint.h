/*
 * mcint.h — Monte-Carlo cubature over an axis-aligned box (machine precision)
 *
 * Estimates ∫_box f dV by averaging the integrand over sample points drawn from
 * the box and multiplying by the volume.  Two point sources:
 *
 *   - pseudo-random (plain / "MonteCarlo"): a deterministic SplitMix64 stream
 *     (fixed seed, so results are reproducible), error ~ V·σ/√N;
 *   - quasi-random (Halton low-discrepancy, "QuasiMonteCarlo"): deterministic,
 *     error ~ V·(log N)^d / N — markedly better for smooth integrands.
 *
 * Points are added in batches until the standard-error estimate meets the
 * tolerance or the evaluation budget is exhausted.  This is the engine that
 * NIntegrate uses for high-dimensional integrals and for region (Boole /
 * UnitStep) integrands, where adaptive cubature is impractical.
 *
 * The integrand is supplied through a sample callback taking the d-vector of
 * coordinates; a sample that returns false is treated as a skipped point.
 */
#ifndef MATHILDA_MCINT_H
#define MATHILDA_MCINT_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

typedef bool (*McSampleMachine)(void* ctx, const double* x, size_t d,
                                double _Complex* out);

bool mc_integrate_machine(McSampleMachine f, void* ctx,
                          const double* a, const double* b, size_t d,
                          bool quasi, double abstol, double reltol,
                          long max_points, double _Complex* result,
                          double* abserr);

#endif /* MATHILDA_MCINT_H */
