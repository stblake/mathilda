/*
 * nfunits.{c,h}
 * -------------
 * Fundamental-unit engine (Gate 2 of the Thue solver) for a monogenic
 * number field K = Q(theta) supplied by numberfield.{c,h}.
 *
 * nf_fundamental_units() returns a CERTIFIED system of fundamental units
 * (rank r = r1 + r2 - 1) or NULL (the caller DECLINEs).  Method:
 *
 *   1. small-norm search for units (|N| == 1, exact);
 *   2. greedily pick r multiplicatively-independent ones (log-embedding
 *      independence), preferring the shortest;
 *   3. certify they generate the FULL unit group (index 1) by p-saturation:
 *      for every prime p <= R / R0 (R0 = 0.2, Friedman's unconditional
 *      regulator lower bound), the mod-p matrix of p-th-power characters
 *      at degree-1 primes must reach full column rank r.  Full rank
 *      (ker = {0}) proves p-saturation UNCONDITIONALLY for any prime set;
 *      if it is not reached within budget, or fewer than r independent
 *      units are found, the routine DECLINEs.
 *
 * Rigour: the search and regulator are heuristic (they only *propose*);
 * the saturation rank test and the norm test are exact.  No path returns
 * a non-fundamental system as fundamental -- at worst it declines.
 *
 * Requires FLINT (via numberfield_internal.h); returns NULL otherwise.
 */
#ifndef NFUNITS_H
#define NFUNITS_H

#include "numberfield.h"
#include <gmp.h>
#include <stdbool.h>

typedef struct NFUnits NFUnits;   /* opaque */

/* Certified fundamental units of K, or NULL to DECLINE.  For unit rank 0
 * returns a valid empty system (rank 0). */
NFUnits* nf_fundamental_units(NumberField* K);
void nf_units_free(NFUnits* U);

int    nf_units_rank(const NFUnits* U);
double nf_units_regulator(const NFUnits* U);

/* theta-power NUMERATOR coordinates (length = field degree) of fundamental
 * unit j, 0 <= j < rank: the unit is (1/D) * coords, D = nf_units_denom.
 * D == 1 for a monogenic field.  Borrowed pointer into U; do not free. */
const mpz_t* nf_units_coords(const NFUnits* U, int j);

/* Shared denominator D of the unit coordinates (1 when Z[theta] = O_K). */
void nf_units_denom(const NFUnits* U, mpz_t out);

#endif /* NFUNITS_H */
