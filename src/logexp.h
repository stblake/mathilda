#ifndef LOGEXP_H
#define LOGEXP_H

#include "expr.h"

Expr* builtin_log(Expr* res);
Expr* builtin_exp(Expr* res);
Expr* builtin_log10(Expr* res);
Expr* builtin_log2(Expr* res);

/* exp(z) for an MPFR real or complex z at z's precision (mpfr_exp / the
 * complex exp(a)(cos b + i sin b)); NULL when z carries no MPFR part.  The
 * single kernel behind Exp[z] and E^z, so the two spellings -- and the
 * compiled MPFR exp, which uses the same MPFR operations -- agree bit for
 * bit. */
Expr* logexp_mpfr_exp(const Expr* z);

void logexp_init(void);

#endif // LOGEXP_H
