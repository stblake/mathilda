#ifndef EXPSPLIT_H
#define EXPSPLIT_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* Powers with a symbolic exponent as polynomial generators (expsplit.c):
 * b^(n + k t) reads as b^n (b^t)^k, as the reference reads it. */
typedef struct {
    Expr** base;   /* b */
    Expr** term;   /* t */
    Expr** sym;    /* the private symbol standing for b^t */
    size_t n, cap;
    unsigned level; /* nesting level of the generator names; 0 = not chosen */
} ExpSplitTab;

/* e with every splittable power rewritten over private generator symbols
 * (evaluated), or NULL when e has none.  Fills t. */
Expr* expsplit_abstract(const Expr* e, ExpSplitTab* t);

/* As expsplit_abstract, but only powers with a numeric base (2^t, (1/3)^t)
 * become generators; x^t and (la^2)^t stay whole. */
Expr* expsplit_abstract_numeric(const Expr* e, ExpSplitTab* t);
/* Every power with a symbolic exponent, a lone (c + d)^n included. */
Expr* expsplit_abstract_all(const Expr* e, ExpSplitTab* t);

/* Map the generators back (K^k -> b^(k t)); consumes e, evaluates. */
Expr* expsplit_restore(Expr* e, const ExpSplitTab* t);

/* Index of generator symbol e in t, or -1; and the generator b^t itself. */
int   expsplit_index(const ExpSplitTab* t, const Expr* e);
Expr* expsplit_generator(const ExpSplitTab* t, size_t i);

void expsplit_free(ExpSplitTab* t);

#endif
