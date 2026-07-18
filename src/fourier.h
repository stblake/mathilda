#ifndef MATHILDA_FOURIER_H
#define MATHILDA_FOURIER_H

/* ---------------------------------------------------------------------------
 * Fourier / InverseFourier — the discrete Fourier transform.
 *
 * Fourier[list] and InverseFourier[list] compute the (inverse) discrete
 * Fourier transform of a 1-D list or a rectangular nested array. Three numeric
 * regimes are handled:
 *
 *   - machine precision  : O(n log n) via FFTW (when built with USE_FFTW),
 *                          else a naive O(n^2) fallback;
 *   - arbitrary precision: a hand-rolled MPFR-complex FFT — radix-2 for
 *                          power-of-two lengths, Bluestein (chirp-z) otherwise;
 *   - symbolic           : the exact transform built from roots of unity,
 *                          left for the evaluator to simplify.
 *
 * Both accept the FourierParameters -> {a, b} option (default {0, 1}) and the
 * position form Fourier[list, {p1, p2, ...}]. See docs/spec for the conventions.
 * -------------------------------------------------------------------------- */

#include "expr.h"

Expr* builtin_fourier(Expr* res);
Expr* builtin_inverse_fourier(Expr* res);
Expr* builtin_fourier_dct(Expr* res);
Expr* builtin_fourier_dst(Expr* res);

void fourier_init(void);

#endif /* MATHILDA_FOURIER_H */
