/* cherry_sigma_decomp.h — Cherry 1986 Theorem 4.4 restricted Sigma-decomposition.
 *
 * Decides whether Phi in K(x) has a Sigma-decomposition
 *
 *     Phi = Sum_i b_i * Prod_j f_j^{alpha_ij},     b_i in K,  alpha_i = g(alpha_i1)
 *
 * over the distinct irreducibles Sigma = (f_1, ..., f_m), for the degree-1
 * ALL-EQUAL restriction g(r) = (r, r, ..., r) used by the logarithmic-integral
 * case (Cherry 1986 Thm 5.3/5.4).  With all exponents equal, every term is
 * b_i * P^{r_i} where P = Prod_j f_j, so this decides membership Phi in K[P]
 * (positive powers) via the faithful Thm 4.4 procedure — multiplicity extraction,
 * b_i = (p mod f_1)/(q mod f_1), recursion, the cross-factor consistency check,
 * and the increasing-case degree-overshoot termination that PROVES
 * non-existence (the decision property behind ElementaryIntegralQ -> False for
 * e.g. Integrate[x^2/Log[x^2-1], x]).
 *
 * This is the engine behind the `Integrate`SigmaDecomposition` builtin and the
 * li non-existence decision wired into cherry_li.c.
 */
#ifndef MATHILDA_CHERRY_SIGMA_DECOMP_H
#define MATHILDA_CHERRY_SIGMA_DECOMP_H

#include <stddef.h>
#include "expr.h"

typedef enum {
    SIGMA_EXISTS = 0,     /* a restricted Sigma-decomposition was found       */
    SIGMA_NONEXISTENT,    /* Thm 4.4 termination PROVES none exists (decision) */
    SIGMA_UNKNOWN         /* out of certified scope — decline, do NOT conclude */
} SigmaStatus;

typedef struct {
    SigmaStatus status;
    size_t n;             /* number of terms                                   */
    size_t m;             /* number of factors (length of each exponent vector) */
    Expr** coeffs;        /* b_i, length n (owned)                             */
    long** exps;          /* alpha_i vectors, n x m (owned)                    */
} SigmaDecomp;

/* Decide/compute the all-equal degree-1 Sigma-decomposition of Phi over the
 * given `factors` (distinct irreducibles in x).  Never throws; returns
 * SIGMA_UNKNOWN rather than guessing when the input is out of scope. */
SigmaDecomp cherry_sigma_decompose(Expr* Phi, Expr** factors, size_t m, Expr* x);

void cherry_sigma_free(SigmaDecomp* d);

/* Integrate`SigmaDecomposition[Phi, {f1, ..., fm}, x] — debuggable surface. */
Expr* builtin_sigma_decomposition(Expr* res);

void cherry_sigma_init(void);

#endif /* MATHILDA_CHERRY_SIGMA_DECOMP_H */
