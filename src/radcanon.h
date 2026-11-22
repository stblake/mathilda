/* radcanon.h -- the canonical form of an exact product of rational radicals.
 *
 * One algorithm serves both heads that build such products, Times
 * (2 Sqrt[3/2] -> Sqrt[6]) and Power ((25/12)^(1/3), (2/3)^(5/3)), so a
 * product is a fixed point under each: see radcanon.c for the derivation. */
#ifndef MATHILDA_RADCANON_H
#define MATHILDA_RADCANON_H

#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>
#include "expr.h"

/* v^e with v > 1 and e a nonzero rational. */
typedef struct {
    mpz_t v;
    mpq_t e;
} RcAtom;

typedef struct {
    RcAtom* a;
    size_t n, cap;
} RcAtoms;

void rc_atoms_init(RcAtoms* at);
void rc_atoms_clear(RcAtoms* at);
/* Append v^e; a v of 1 (or e of 0) is dropped. */
void rc_atoms_add(RcAtoms* at, const mpz_t v, const mpq_t e);

/* The canonical factors of
 *     coef * prod(atoms) * (-1)^neg_e
 * coef is a nonzero Gaussian rational re + I im; neg_e (NULL when there is
 * none) is the exponent, 0 < neg_e < 1, of the one negative-base power
 * (-n)^neg_e taking part, whose magnitude n^neg_e is already among the
 * atoms.  Fills *out with a malloc'd array of unevaluated factors (the
 * number first, omitted when it is 1) and returns their count. */
size_t rc_canonical(const RcAtoms* at, const mpq_t coef_re, const mpq_t coef_im,
                    const mpq_t* neg_e, Expr*** out);

#endif
