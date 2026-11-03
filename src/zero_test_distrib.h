/* zero_test_distrib.h -- PossibleZeroQ stage for generalized functions
 * (DiracDelta, Derivative[k][DiracDelta], HeavisideTheta).  See
 * zero_test_distrib.c for the decision procedure. */
#ifndef MATHILDA_ZERO_TEST_DISTRIB_H
#define MATHILDA_ZERO_TEST_DISTRIB_H

#include "expr.h"
#include "zero_test.h"
#include <stdbool.h>

/* True iff `e` mentions DiracDelta, a DiracDelta derivative or HeavisideTheta. */
bool zt_has_distribution(const Expr* e);

/* Decide e == 0 as a distribution in its single variable: every delta-basis
 * coefficient (after sifting) and every regular piece between breakpoints must
 * be zero; each scalar question goes to `decide` (zero_test_decide).  UNKNOWN
 * when `e` is not a distribution of the supported shape.  `e` borrowed. */
ZeroTestResult zt_decide_distribution(const Expr* e,
                                      ZeroTestResult (*decide)(const Expr*));

#endif
