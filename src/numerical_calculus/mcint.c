/*
 * mcint.c — Monte-Carlo / quasi-Monte-Carlo box cubature.  See mcint.h.
 */

#include "mcint.h"

#include <math.h>
#include <stdint.h>

/* SplitMix64 — a small, fast, deterministic PRNG (fixed seed for reproducible
 * Monte-Carlo results, since the environment provides no entropy source). */
static uint64_t sm_next(uint64_t* s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double sm_unit(uint64_t* s) {
    return (double)(sm_next(s) >> 11) * (1.0 / 9007199254740992.0);  /* [0,1) */
}

/* Radical-inverse (van der Corput) in the given base — the Halton coordinate. */
static double halton(long i, int base) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= base; r += f * (double)(i % base); i /= base; }
    return r;
}

static const int MC_PRIMES[24] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89
};

bool mc_integrate_machine(McSampleMachine f, void* ctx,
                          const double* a, const double* b, size_t d,
                          bool quasi, double abstol, double reltol,
                          long max_points, double _Complex* result,
                          double* abserr) {
    *result = 0.0; *abserr = INFINITY;
    if (d == 0 || d > 24) return false;
    double vol = 1.0;
    for (size_t j = 0; j < d; j++) vol *= (b[j] - a[j]);
    if (max_points <= 0) max_points = 2000000;
    if (reltol <= 0.0) reltol = 1e-2;

    uint64_t rng = 0x123456789ABCDEFULL;
    long halton_i = 1;
    double x[24];

    double _Complex sum = 0.0;
    double sum2 = 0.0;        /* sum of |f|^2, for the variance estimate */
    long n = 0;
    const long batch = 2048;
    double _Complex est = 0.0;
    double err = INFINITY;
    bool converged = false;

    while (n < max_points) {
        for (long k = 0; k < batch; k++) {
            for (size_t j = 0; j < d; j++) {
                double u = quasi ? halton(halton_i, MC_PRIMES[j]) : sm_unit(&rng);
                x[j] = a[j] + u * (b[j] - a[j]);
            }
            halton_i++;
            double _Complex fv;
            if (!f(ctx, x, d, &fv)) continue;
            if (!isfinite(creal(fv)) || !isfinite(cimag(fv))) continue;
            sum += fv;
            sum2 += creal(fv) * creal(fv) + cimag(fv) * cimag(fv);
            n++;
        }
        if (n < 2) continue;
        double _Complex mean = sum / (double)n;
        est = vol * mean;
        /* standard error of the mean, scaled by the volume */
        double var = sum2 / (double)n - (creal(mean) * creal(mean) + cimag(mean) * cimag(mean));
        if (var < 0.0) var = 0.0;
        err = vol * sqrt(var / (double)n);
        if (err <= fmax(abstol, reltol * cabs(est))) { converged = true; break; }
    }

    *result = est;
    *abserr = err;
    return converged;
}
