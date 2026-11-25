/* linalg_mod.h - exact linear algebra over Z/p for the Modulus -> p option
 * of Det, Inverse, NullSpace and RowReduce (see linalg_mod.c). */
#ifndef LINALG_MOD_H
#define LINALG_MOD_H

#include "expr.h"
#include <stdint.h>

/* Largest modulus handled: products of two residues must fit an int64_t. */
#define LINALG_MOD_MAX INT64_C(2147483647)

/* Classify one option argument.
 *   0   opt is not a Modulus rule (the caller parses it as something else);
 *   1   opt is Modulus -> p with 0 <= p <= LINALG_MOD_MAX (*p set; p == 0 is
 *       the ordinary field, i.e. no modulus);
 *  -1   opt is Modulus -> v with any other v: the call stays unevaluated. */
int linalg_mod_parse_option(const Expr* opt, int64_t* p);

/* Each returns a fresh result, or NULL (call left unevaluated) when an entry
 * is not an integer or a rational with a denominator invertible mod p. p >= 2.
 * `m` is borrowed. */
Expr* linalg_mod_det(Expr* m, int64_t p);
Expr* linalg_mod_inverse(Expr* m, int64_t p);  /* Inverse::sing if singular  */
Expr* linalg_mod_nullspace(Expr* m, int64_t p);
Expr* linalg_mod_rowreduce(Expr* m, int64_t p, const char* head);

#endif /* LINALG_MOD_H */
