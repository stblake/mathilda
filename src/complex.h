#ifndef PICO_COMPLEX_H
#define PICO_COMPLEX_H

#pragma GCC system_header
#include_next <complex.h>

#include "expr.h"

Expr* builtin_re(Expr* res);
Expr* builtin_im(Expr* res);
Expr* builtin_reim(Expr* res);
Expr* builtin_abs(Expr* res);
Expr* builtin_sign(Expr* res);
Expr* builtin_conjugate(Expr* res);
Expr* builtin_arg(Expr* res);

void complex_init(void);

#endif // PICO_COMPLEX_H
