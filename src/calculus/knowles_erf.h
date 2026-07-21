/* knowles_erf.h — Knowles' error-function integration of transcendental
 * Liouvillian functions (KNOWLES_DESIGN.md, K2).  See knowles_erf.c. */
#ifndef MATHILDA_KNOWLES_ERF_H
#define MATHILDA_KNOWLES_ERF_H

#include "expr.h"

/* Integrate f dx over the K0 Liouvillian-primitive tower in terms of elementary
 * functions + Erf/Erfi.  Returns a freshly-owned antiderivative (diff-back
 * verified) or NULL to decline.  Registered as an RT_SPECIAL_FORMS entry. */
Expr* knowles_erf_liouvillian(Expr* f, Expr* x);

#endif /* MATHILDA_KNOWLES_ERF_H */
