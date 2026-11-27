/* intrat_alg.h — native (FLINT) algorithms of Bronstein, Symbolic Integration
 * I, chapter 2, for the rational-function integrator.  See intrat_alg.c.
 *
 * Scope: integrands with exact rational coefficients (A, D in Q[x]).  Every
 * entry point returns NULL to decline -- a coefficient outside Q, a size cap,
 * FLINT compiled out -- and the caller keeps its generic (Expr) path, so
 * declining is always safe.  Arguments are borrowed; results are owned. */
#ifndef MATHILDA_INTRAT_ALG_H
#define MATHILDA_INTRAT_ALG_H

#include <stdbool.h>
#include "expr.h"

/* Hermite reduction variants (Bronstein §2.2, §2.3). */
enum {
    INTRAT_HERMITE_LINEAR = 0,     /* Mack's linear version            */
    INTRAT_HERMITE_QUADRATIC = 1,  /* the original quadratic version   */
    INTRAT_HERMITE_HOROWITZ = 2    /* Horowitz-Ostrogradsky (linear algebra) */
};

/* Log-part algorithms (Bronstein §2.4-2.6). */
enum {
    INTRAT_LOGPART_LRT = 0,         /* Lazard-Rioboo-Trager (subresultants) */
    INTRAT_LOGPART_RT = 1,          /* Rothstein-Trager (resultant, factor, gcd over Q(a)) */
    INTRAT_LOGPART_CZICHOWSKI = 2   /* Czichowski (lex Groebner basis)      */
};

/* HermiteReduce(A/D): {g, h} with A/D = D[g, x] + h, Denominator[h]
 * squarefree, by the given variant.  deg A < deg D is NOT required (the
 * polynomial part is returned in h).  NULL to decline. */
Expr* intrat_alg_hermite(const Expr* A, const Expr* D, const Expr* x, int variant);

/* The log part of A/D (D squarefree, deg A < deg D, gcd(A, D) = 1) as the
 * pair list {{Q_1[t], S_1[t, x]}, ...}: Integral = Sum_i Sum_{Q_i(a)=0} a Log[S_i(a, x)].
 * method is INTRAT_LOGPART_RT or INTRAT_LOGPART_CZICHOWSKI.  NULL to decline. */
Expr* intrat_alg_logpart(const Expr* A, const Expr* D, const Expr* x, const Expr* t, int method);

/* Bronstein's LogToReal (§2.8) for an arbitrary pair list over Q: the real
 * form of Sum_i Sum_{Q_i(a)=0} a Log[S_i(a, x)] -- for each real root a,
 * a Log[S(a, x)]; for each root u + I v with v > 0, u Log[A^2 + B^2] +
 * v LogToAtan(A, B), where S(u + I v, x) = A + I B, by Rioboo's continuous
 * LogToAtan.  The roots are exact (FLINT qqbar).  allow_root_objects: when
 * false, decline as soon as a coefficient is not rational or a quadratic
 * surd (the caller then keeps a RootSum); when true, such a coefficient is
 * written as a Root object.  NULL to decline. */
Expr* intrat_alg_log_to_real(const Expr* pairs, const Expr* x, const Expr* t, bool allow_root_objects);

#endif /* MATHILDA_INTRAT_ALG_H */
