#ifndef KELVIN_H
#define KELVIN_H

#include "expr.h"

/* KelvinBer / KelvinBei / KelvinKer / KelvinKei [nu, z] and the order-0
 * one-argument forms:
 *   ber_nu(x) + i bei_nu(x) = J_nu(x e^(3 pi i/4)),
 *   ker_nu(x) + i kei_nu(x) = e^(-nu pi i/2) K_nu(x e^(pi i/4)). */
Expr* builtin_kelvinber(Expr* res);
Expr* builtin_kelvinbei(Expr* res);
Expr* builtin_kelvinker(Expr* res);
Expr* builtin_kelvinkei(Expr* res);
void  kelvin_init(void);

#endif /* KELVIN_H */
