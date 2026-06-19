#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H

#include "expr.h"

/* Number-theoretic builtins, split out of arithmetic.c.
 *
 * arithmetic.c retains the core rational/complex constructors and the
 * shared int64 gcd()/lcm() helpers; this module owns the higher-level
 * integer functions: GCD, LCM, ExtendedGCD, PowerMod, Factorial,
 * Factorial2, FactorialPower, Binomial, PrimitiveRoot, PrimitiveRootList,
 * MultiplicativeOrder.
 *
 * Registration, attributes, and docstrings live in numbertheory_init(),
 * called from core_init() alongside the other subsystem initialisers. */
void numbertheory_init(void);

Expr* builtin_gcd(Expr* res);
Expr* builtin_lcm(Expr* res);
Expr* builtin_extendedgcd(Expr* res);
Expr* builtin_powermod(Expr* res);
Expr* builtin_primitiveroot(Expr* res);
Expr* builtin_primitiverootlist(Expr* res);
Expr* builtin_multiplicativeorder(Expr* res);
Expr* builtin_factorial(Expr* res);
Expr* builtin_factorial2(Expr* res);
Expr* builtin_factorialpower(Expr* res);
Expr* builtin_binomial(Expr* res);
Expr* builtin_divisible(Expr* res);

#endif
