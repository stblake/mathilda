/*
 * ratfun.h -- factored rational functions over Q for Together / Cancel / Apart.
 *
 * One conversion of the expression into FLINT (fmpz_mpoly over the canonical
 * variable order, transcendental kernels as generators), the whole algorithm
 * native, one conversion back.  The answer keeps the denominator FACTORED, as a
 * product of a gcd-free basis of the input's denominator factors, with the
 * numeric content in front -- the Together / Cancel output form:
 *
 *   Together[1/(x-1) + 1/(x+1)]      -> (2 x)/((-1 + x) (1 + x))
 *   Together[1/(2 x + 4) + 1/(6x-3)] -> (1 + 8 x)/(6 (2 + x) (-1 + 2 x))
 *   Cancel[(x^2-1)^2/((x-1)(x+5))]   -> ((-1 + x) (1 + x)^2)/(5 + x)
 *
 * Every entry returns NULL (decline; the caller keeps its other paths) on an
 * input outside Q(vars, kernels): algebraic constants (I, Sqrt[2]), symbolic or
 * fractional powers other than E^u, forward trig, inexact numbers.
 */
#ifndef MATHILDA_RATFUN_H
#define MATHILDA_RATFUN_H

#include <gmp.h>
#include "expr.h"

/* Together[e]: combine every sum of fractions over its least common
 * denominator and cancel. */
Expr* ratfun_together(const Expr* e);

/* Cancel[e]: cancel common factors in each term of a top-level sum (the terms
 * are not combined); nested sums of fractions are combined. */
Expr* ratfun_cancel(const Expr* e);

/* Apart[e, x]: partial fractions in x over Q(the other variables).  `x` is a
 * symbol of the input. */
Expr* ratfun_apart(const Expr* e, const Expr* x);

/* The same over GF(p), for Modulus -> p (p a prime).  Apart: one variable. */
Expr* ratfun_together_mod(const Expr* e, const mpz_t p);
Expr* ratfun_cancel_mod(const Expr* e, const mpz_t p);
Expr* ratfun_apart_mod(const Expr* e, const Expr* x, const mpz_t p);

/* The Modulus option of `res` (a Together/Cancel/Apart call, options after
 * the first argument).  Returns 1 with p set for Modulus -> p (p prime);
 * -1 when the call must stay unevaluated (head::modp for a composite
 * modulus, head::mext with an Extension; the message is emitted); 0
 * otherwise -- no Modulus, a symbolic one, or Modulus -> 0 -- and then, when
 * the call carries Modulus -> 0 or Trig -> False (the defaults), *stripped
 * (caller frees) is the call without them. */
int ratfun_modulus_option(const Expr* res, const char* head, mpz_t p, Expr** stripped);

#endif /* MATHILDA_RATFUN_H */
