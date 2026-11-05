/* diractr.h -- Experimental`FeynCalc`DiracTr, a native Dirac trace kernel.
 *
 * See diractr.c for the definition, the algorithm and the exactness contract.
 */
#ifndef MATHILDA_EXPERIMENTAL_FEYNCALC_DIRACTR_H
#define MATHILDA_EXPERIMENTAL_FEYNCALC_DIRACTR_H

#include "expr.h"

/* Experimental`FeynCalc`DiracTr[a1, ..., an] (builtin; owns nothing of `res`
 * on a NULL return, per the builtin ownership contract). */
Expr* builtin_experimental_diractr(Expr* res);

/* Registers every Experimental`FeynCalc` builtin (called from core_init). */
void experimental_feyncalc_init(void);

#endif
