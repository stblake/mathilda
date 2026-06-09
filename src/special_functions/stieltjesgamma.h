#ifndef STIELTJESGAMMA_H
#define STIELTJESGAMMA_H

#include "expr.h"

/* StieltjesGamma[n] -- the n-th Stieltjes constant gamma_n, the coefficients
 * in the Laurent expansion of the Riemann zeta function about s = 1:
 *
 *   zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n / n!) gamma_n (s-1)^n.
 *
 * gamma_0 = EulerGamma. StieltjesGamma is an inert, indexed constant: it
 * carries no closed-form evaluation (it stays symbolic) except for the
 * reduction StieltjesGamma[0] -> EulerGamma. It is produced by the series
 * expansion of Zeta about s = 1.
 *
 * Attributes: Listable, Protected. */
Expr* builtin_stieltjesgamma(Expr* res);

void stieltjesgamma_init(void);

#endif /* STIELTJESGAMMA_H */
