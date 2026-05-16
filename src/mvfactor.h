/*
 * mvfactor.h
 * ----------
 * Multivariate polynomial factorisation: the high-level orchestration
 * layer that ties together the existing univariate Berlekamp-Zassenhaus
 * pipeline (facpoly.c::bz_factor_to_expr) with the new bivariate
 * primitives (bpoly.h, zupoly.h).
 *
 * Currently exposes:
 *   - bpoly_hensel_lift_2: bivariate two-factor Hensel iteration.
 *
 * Future additions (see FACTOR_PLAN.md):
 *   - Multifactor Hensel lift.
 *   - Wang's leading-coefficient correction for non-monic inputs.
 *   - Zassenhaus recombination at the bivariate level.
 *   - n-variate dispatch via specialisation + recursive lifting.
 *   - mvfactor_factor: the public entry point that builtin_factor will
 *     route to once the pipeline is complete.
 */

#ifndef MVFACTOR_H
#define MVFACTOR_H

#include <stdbool.h>
#include "bpoly.h"
#include "zupoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Bivariate two-factor Hensel lift.
 *
 * Given:
 *   - P(x, y) ∈ Z[x, y], monic in x (i.e. lc_x(P) is the constant 1
 *     when viewed as a polynomial in x with y-coefficients).
 *   - u(x), v(x) ∈ Z[x]: monic, coprime over Q[x], with
 *     u(x) * v(x) = P(x, 0).
 *
 * Returns U(x, y), V(x, y) ∈ Z[x, y] satisfying:
 *   - U(x, 0) = u(x), V(x, 0) = v(x).
 *   - U * V = P (over Z[x, y], not just modulo y^B).
 *   - U is monic in x; V is monic in x.
 *
 * The algorithm:
 *   1. Initialise U_0 = u, V_0 = v (each as a BPoly of y-degree 0).
 *   2. For k = 1, 2, ..., deg_y(P):
 *      a. Compute the error E_k = coefficient of y^k in P - U*V.
 *      b. Solve  Δu · v + Δv · u = E_k   (zupoly_diophantine).
 *      c. Update U += y^k · Δu, V += y^k · Δv.
 *   3. Verify U * V == P exactly; return true iff so.
 *
 * Returns false if any Diophantine step fails (typically because the
 * requested factorisation does not exist over Z, or u and v fail to
 * satisfy the monic-coprime preconditions).  On failure *U_out and
 * *V_out are NULL.  On success the caller owns *U_out and *V_out.
 */
bool bpoly_hensel_lift_2(const BPoly* P, const ZUPoly* u, const ZUPoly* v,
                        BPoly** U_out, BPoly** V_out);

/* Bivariate two-factor Hensel lift with PREDICTED leading coefficients
 * (Phase F1 Stage 3 -- Wang's leading-coefficient correction for the
 * polynomial-in-y LC case).
 *
 * Given:
 *   - P(x, y) ∈ Z[x, y] with lc_x(P)(y) = q_u(y) * q_v(y).
 *   - u(x), v(x) ∈ Z[x]: monic, coprime over Q[x], with
 *     u(x) * v(x) = P(x, 0) / (q_u(0) * q_v(0)).
 *   - q_u(y), q_v(y) ∈ Z[y]: predicted leading coefficients with
 *     q_u(0) and q_v(0) being the (constant) values that the lifted
 *     factors' leading-x coefficients take at y = 0.  Their product
 *     must equal lc_x(P)(y) exactly.
 *
 * Returns U(x, y), V(x, y) ∈ Z[x, y] satisfying:
 *   - U(x, 0) = q_u(0) * u(x);  V(x, 0) = q_v(0) * v(x).
 *   - lc_x(U)(y) = q_u(y);  lc_x(V)(y) = q_v(y).
 *   - U * V = P over Z[x, y].
 *
 * The algorithm differs from `bpoly_hensel_lift_2` only at the
 * Diophantine step: each iteration's correction Δu has its leading
 * x-coefficient PINNED to the y^k coefficient of q_u (so the predicted
 * lc_x(U) = q_u is maintained), and only the lower x-degree part is
 * solved via zupoly_diophantine.
 *
 * Initial values absorb q_*(0) into the seed:
 *   U_0(x) = q_u(0) * u(x);  V_0(x) = q_v(0) * v(x).
 * If q_u(0) = q_v(0) = 1 the initial values reduce to u, v as in the
 * ordinary monic lift.
 *
 * Returns false if any Diophantine step fails or the final
 * verification U*V == P does not hold.  On failure *U_out and *V_out
 * are NULL.  On success the caller owns *U_out and *V_out.
 */
bool bpoly_hensel_lift_2_lc(const BPoly* P,
                             const ZUPoly* u, const ZUPoly* v,
                             const ZUPoly* q_u, const ZUPoly* q_v,
                             BPoly** U_out, BPoly** V_out);

/* Multifactor bivariate Hensel lift with subset recombination.
 *
 * Given:
 *   - P(x, y) ∈ Z[x, y], monic in x.
 *   - u[0], u[1], ..., u[r-1]: monic, pairwise coprime polynomials in
 *     Z[x] with prod(u[i]) = P(x, 0).
 *
 * Returns a freshly-allocated array of `*r_out` BPoly* such that:
 *   - prod(U[i]) = P over Z[x, y].
 *   - Each U[i] is monic in x.
 *   - Each U[i](x, 0) is a sub-product of the input u's:
 *       Specifically, the input u's partition into r' = *r_out groups
 *       S_0, S_1, ..., S_{r'-1}, with U[i](x, 0) = prod_{j in S_i} u[j].
 *
 * The output count *r_out may be LESS than r when the bivariate
 * factorisation is coarser than the univariate one.  Example:
 *   P(x,y) = x^4 - y^2 = (x^2 - y)(x^2 + y).
 *   At y=1: P(x,1) = x^4 - 1 = (x-1)(x+1)(x^2+1).  Univariately r=3.
 *   Bivariately r'=2, with (x-1)(x+1) → x^2 - y and (x^2+1) → x^2 + y.
 * In that case Phase 3a's pair-and-recurse fails (peeling off (x-1)
 * vs (x+1)(x^2+1) doesn't lift to a true bivariate factor); Phase 3b
 * tries all 2-subsets and finds the {(x-1), (x+1)} grouping that does.
 *
 * Algorithm:
 *   1. If r == 1, P is the lone factor.
 *   2. For k = 1, 2, ..., r/2:
 *        For each k-subset S of {0, ..., r-1}:
 *          Run two-factor lift with uS = prod(u[i] for i in S) and
 *          uC = prod(u[i] for i NOT in S).
 *          On success: U_S is a true bivariate factor; recurse on
 *          the complement (r - k of the original u's, with V = U_C
 *          as the corresponding bivariate residual).
 *   3. If no subset yields a successful lift, P is irreducible
 *      bivariately under this set of u's; return P as the sole factor.
 *
 * Worst-case complexity O(2^r) in the number of two-factor lifts, but
 * typically dominated by k = 1 (pair-and-recurse) when the univariate
 * and bivariate factorisations agree.  For squarefree images with
 * r ≤ 10 the worst case is tractable.
 *
 * The caller owns the returned array and each BPoly* in it; free
 * them with bpoly_free for each entry, then free the outer array.
 *
 * Returns true on success (always — on irreducible P the result is
 * a single-element array containing a copy of P).  Returns false
 * only on internal allocation failure or invalid input (r <= 0).
 */
bool bpoly_hensel_lift_multi(const BPoly* P, const ZUPoly** us, int r,
                             BPoly*** Us_out, int* r_out);

/* High-level orchestrator: factor a monic-in-x bivariate polynomial
 * P ∈ Z[x, y] into bivariate factors over Z[x, y].
 *
 * Algorithm:
 *   1. Try several integer evaluation points α to find one where:
 *        - lc_x(P) does not vanish at α (preserves x-degree).
 *        - P(x, α) is squarefree.
 *      Skip points that fail either test; bound the number of tries.
 *   2. Factor P(x, α) univariately (caller must supply this; we
 *      receive the factors as ZUPoly inputs).
 *   3. Shift y -> y + α so the lift starts from y = 0, run the
 *      multifactor lift, then shift back.
 *
 * To keep the layering clean -- this module does not know about
 * Mathilda's facpoly.c::bz_factor_to_expr -- the caller passes a
 * factoring callback that produces the univariate factors at a
 * given integer point.
 *
 * Returns true on success.  *factors_out is a freshly-allocated
 * array of `*r_out` BPoly* representing the factorisation; the
 * caller owns the array and each BPoly* in it.
 *
 * Returns false if no good evaluation point was found within the
 * built-in budget, or if the lift failed (non-integer Diophantine
 * solution -- typically because the multivariate factorisation
 * does not exist over Z, e.g. P is irreducible). */
typedef bool (*mvfactor_univariate_cb)(const ZUPoly* image,
                                       ZUPoly*** factors_out,
                                       int* count_out,
                                       void* user_data);

bool mvfactor_try_bivariate_monic(const BPoly* P,
                                  mvfactor_univariate_cb factor_cb,
                                  void* cb_user_data,
                                  BPoly*** factors_out,
                                  int* r_out);

#ifdef __cplusplus
}
#endif

#endif /* MVFACTOR_H */
