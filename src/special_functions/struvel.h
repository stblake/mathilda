#ifndef STRUVEL_H
#define STRUVEL_H

#include "expr.h"

/* StruveL[nu, z] -- the modified Struve function L_nu(z):
 *   L_nu(z) = (z/2)^(nu+1) Sum_k (z/2)^(2k) / (Gamma[k+3/2] Gamma[k+nu+3/2]). */
Expr* builtin_struvel(Expr* res);
void  struvel_init(void);

#endif /* STRUVEL_H */
