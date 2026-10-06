#ifndef STRUVEH_H
#define STRUVEH_H

#include "expr.h"

/* StruveH[nu, z] -- the Struve function H_nu(z), the particular solution of
 *   z^2 y'' + z y' + (z^2 - nu^2) y = (4 (z/2)^(nu+1)) / (Sqrt[Pi] Gamma[nu+1/2]).
 * Series (entire in z^2):
 *   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2])) 1F2(1; 3/2, nu+3/2; -z^2/4). */
Expr* builtin_struveh(Expr* res);
void  struveh_init(void);

#endif /* STRUVEH_H */
