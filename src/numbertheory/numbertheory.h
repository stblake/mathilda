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
 * The implementation is split one-file-per-builtin under src/numbertheory/
 * (mirroring src/linalg/); cross-cutting private helpers live in
 * numbertheory_internal.h.  Registration, attributes, and docstrings live
 * in numbertheory_init(), called from core_init() alongside the other
 * subsystem initialisers. */
void numbertheory_init(void);

/* Mod[m, n] / Mod[m, n, d] for Gaussian integers (at least one non-real), with
 * the nearest-quotient convention m - n Round[(m - d)/n]. NULL when
 * an argument is not a Gaussian integer or n is 0. Used by builtin_mod. */
Expr* nt_gaussian_mod(Expr* m, Expr* n, Expr* d);

/* Parse the GaussianIntegers -> True | False option from args[first ..] with
 * head::nonopt / ::optx / ::opttf on error (returns false; leave the call
 * unevaluated). *set reports whether it appeared; the first setting wins. */
#include <stdbool.h>
#include <stddef.h>
bool nt_parse_gaussian_options(Expr* res, size_t first, const char* head,
                               bool* gaussian, bool* set);

Expr* builtin_gcd(Expr* res);
Expr* builtin_lcm(Expr* res);
Expr* builtin_extendedgcd(Expr* res);
Expr* builtin_chineseremainder(Expr* res);
Expr* builtin_powermod(Expr* res);
Expr* builtin_primitiveroot(Expr* res);
Expr* builtin_primitiverootlist(Expr* res);
Expr* builtin_multiplicativeorder(Expr* res);
Expr* builtin_jacobisymbol(Expr* res);
Expr* builtin_factorial(Expr* res);
Expr* builtin_factorial2(Expr* res);
Expr* builtin_factorialpower(Expr* res);
Expr* builtin_binomial(Expr* res);
Expr* builtin_multinomial(Expr* res);
Expr* builtin_divisible(Expr* res);
Expr* builtin_coprimeq(Expr* res);
Expr* builtin_divisors(Expr* res);
Expr* builtin_divisorsigma(Expr* res);
Expr* builtin_moebiusmu(Expr* res);
Expr* builtin_liouvillelambda(Expr* res);
Expr* builtin_mangoldtlambda(Expr* res);
void  mangoldtlambda_init(void);   /* registration, attributes, docstring */
Expr* builtin_prime(Expr* res);
Expr* builtin_primepi(Expr* res);
Expr* builtin_primeomega(Expr* res);
Expr* builtin_primenu(Expr* res);

#endif
