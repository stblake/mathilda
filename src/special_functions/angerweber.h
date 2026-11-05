#ifndef ANGERWEBER_H
#define ANGERWEBER_H

#include "expr.h"

/* AngerJ[nu, z] -- the Anger function J_nu(z) = (1/Pi) Int_0^Pi cos(nu t - z sin t) dt.
 * WeberE[nu, z] -- the Weber function E_nu(z) = (1/Pi) Int_0^Pi sin(nu t - z sin t) dt. */
Expr* builtin_angerj(Expr* res);
Expr* builtin_webere(Expr* res);
void  angerweber_init(void);

#endif /* ANGERWEBER_H */
