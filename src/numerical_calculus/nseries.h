/*
 * nseries.h — NSeries[f, {x, x0, n}, opts]
 *
 * Numerically approximates the Taylor/Laurent series expansion of `f` about
 * x = x0, including the terms (x - x0)^-n through (x - x0)^n, returned as a
 * SeriesData object. f is sampled on a circle in the complex plane centred at
 * x0; a discrete Fourier transform of the samples recovers the Laurent
 * coefficients via Cauchy's integral formula (Lyness & Sande, 1971).
 *
 * Attributes: Protected (NOT Listable — generic threading would split the
 * {x, x0, n} spec).
 */
#ifndef MATHILDA_NSERIES_H
#define MATHILDA_NSERIES_H

#include "expr.h"

Expr* builtin_nseries(Expr* res);
void  nseries_init(void);

#endif /* MATHILDA_NSERIES_H */
