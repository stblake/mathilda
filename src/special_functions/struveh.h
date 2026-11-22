#ifndef STRUVEH_H
#define STRUVEH_H

#include "expr.h"
#include <stdbool.h>

/* StruveH[nu, z] -- the Struve function H_nu(z), the particular solution of
 *   z^2 y'' + z y' + (z^2 - nu^2) y = (4 (z/2)^(nu+1)) / (Sqrt[Pi] Gamma[nu+1/2]).
 * Series (entire in z^2):
 *   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2])) 1F2(1; 3/2, nu+3/2; -z^2/4). */
Expr* builtin_struveh(Expr* res);

/* H_{-n}(z) (modified false) or L_{-n}(z) (modified true) for an exact Integer
 * nu = -n <= -1, reflected to order +n plus an exact Laurent polynomial in z:
 *   H_{-n} = (-1)^n H_n + Sum_{k<n} (-1)^k (z/2)^(2k-n+1)/(Gamma(k+3/2) Gamma(k-n+3/2))
 *   L_{-n} =        L_n + Sum_{k<n}        (z/2)^(2k-n+1)/(Gamma(k+3/2) Gamma(k-n+3/2))
 * Evaluated result, or NULL when nu is not such an integer. Borrows nu, z. */
Expr* struve_negative_integer_order(const Expr* nu, const Expr* z, bool modified);
void  struveh_init(void);

#endif /* STRUVEH_H */
