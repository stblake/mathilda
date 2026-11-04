/* Mathilda -- MathieuC[a, q, z] / MathieuCPrime[a, q, z] (the even Mathieu
 * function and its z-derivative) and MathieuS / MathieuSPrime (the odd one).
 * See mathieu.c for the evaluation strategy and sf_machine_mathieu.h for the
 * definition and normalization. */
#ifndef MATHILDA_MATHIEU_H
#define MATHILDA_MATHIEU_H

#include "expr.h"

/* Builtin entry points. Take ownership of `res` per the builtin contract. */
Expr* builtin_mathieuc(Expr* res);
Expr* builtin_mathieucprime(Expr* res);
Expr* builtin_mathieus(Expr* res);
Expr* builtin_mathieusprime(Expr* res);

/* Register MathieuC / MathieuCPrime / MathieuS / MathieuSPrime and their
 * attributes. Called from core_init(). */
void mathieu_init(void);

#endif /* MATHILDA_MATHIEU_H */
