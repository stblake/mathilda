/*
 * mvfactor3.h
 * -----------
 * Trivariate (and higher) Hensel factoring on top of the MPoly
 * substrate.  Phase F2 of plans/FACTOR_PLAN.md.
 *
 * Contracts and invariants follow the same pattern as mvfactor.h:
 *   - Functions returning new MPoly* / ZUPoly* always allocate fresh.
 *   - Pointer parameters are borrowed unless otherwise noted.
 *   - Multivariate Hensel preconditions: P monic in the main variable;
 *     factors u, v at the chosen alpha-tuple are pairwise coprime over
 *     the rationals and squarefree as a univariate image.
 *
 * The trivariate lift is built on two new primitives:
 *   - `mpoly_diophantine_2`: bivariate Diophantine solver in
 *     Z[var_main, var_y], implemented as a Hensel iteration in var_y
 *     on top of zupoly_diophantine.
 *   - `mpoly_hensel_lift_3_2`: lifts a bivariate factorisation in
 *     (var_main, var_y) to a trivariate factorisation in (var_main,
 *     var_y, var_z) via a Hensel iteration in var_z.
 *
 * For n > 3 the same pattern composes recursively (lift one new
 * variable at a time), but the MVP only handles n == 3.
 */

#ifndef MVFACTOR3_H
#define MVFACTOR3_H

#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

#include "mpoly.h"
#include "zupoly.h"
#include "bpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------- */
/*  MPoly <-> ZUPoly conversions                                    */
/* ---------------------------------------------------------------- */

/* Project an MPoly that is univariate in var_main (every term has
 * exponent 0 for all other variables) to a ZUPoly.  Returns NULL if
 * the MPoly contains terms with non-zero exponents in vars other than
 * var_main. */
ZUPoly* mpoly_to_zupoly_in(const MPoly* p, int var_main);

/* Inflate a ZUPoly to an n_vars-variable MPoly with all degrees in
 * vars != var_main equal to 0. */
MPoly* zupoly_to_mpoly_in(const ZUPoly* z, int n_vars, int var_main);

/* Project an MPoly to a BPoly{var_main, var_y}.  Every term must have
 * exponent 0 for vars outside {var_main, var_y}; returns NULL if not. */
BPoly* mpoly_to_bpoly_in(const MPoly* p, int var_main, int var_y);

/* Inflate a BPoly to an n_vars-variable MPoly with var_main and var_y
 * occupying the chosen indices and all other vars at exponent 0. */
MPoly* bpoly_to_mpoly_in(const BPoly* b, int n_vars,
                         int var_main, int var_y);

/* ---------------------------------------------------------------- */
/*  Bivariate Diophantine solver                                    */
/* ---------------------------------------------------------------- */

/* Solve  delta_u * V + delta_v * U = E  over Z[var_main, var_y].
 *
 * Preconditions:
 *   - U, V are monic in var_main (their leading var_main-coefficient,
 *     viewed as a polynomial in the other variables, is the constant 1).
 *   - U, V are coprime over Q(var_y)[var_main] (equivalently, U|var_y=alpha_y
 *     and V|var_y=alpha_y are coprime in Q[var_main] for at least one
 *     integer alpha_y).
 *   - E is integer-coefficient.
 *
 * The output satisfies:
 *   - deg_{var_main}(delta_u) < deg_{var_main}(U).
 *   - deg_{var_main}(delta_v) < deg_{var_main}(V).
 *
 * Algorithm: specialise var_y to a chosen integer alpha_y where the
 * specialised U, V remain coprime and degree-preserving in var_main;
 * solve the univariate Diophantine via zupoly_diophantine; lift the
 * solution back to bivariate via Hensel iteration in var_y.
 *
 * Returns false if no good alpha_y is found in the search budget, or
 * if the lift fails to converge to integer coefficients.  On failure
 * *delta_u_out, *delta_v_out are NULL. */
bool mpoly_diophantine_2(const MPoly* U, const MPoly* V, const MPoly* E,
                         int var_main, int var_y,
                         MPoly** delta_u_out, MPoly** delta_v_out);

/* ---------------------------------------------------------------- */
/*  Trivariate two-factor Hensel lift                               */
/* ---------------------------------------------------------------- */

/* Lift a bivariate factorisation U_xy * V_xy = P|var_z=alpha_z  to a
 * trivariate factorisation U * V = P.
 *
 * Preconditions:
 *   - P is monic in var_main (lc_{var_main}(P) is the integer 1 over
 *     Z[var_y, var_z], or equivalently every leading-x term has all
 *     other-variable exponents == 0 and value 1).
 *   - U_xy, V_xy ∈ MPoly{P->n_vars vars}, with var_z-exponent == 0 in
 *     every term, monic in var_main, coprime, and
 *     U_xy * V_xy == P|var_z=alpha_z.
 *
 * Algorithm: shift var_z by alpha_z; iterate var_z-degrees k = 1..B_z
 * computing the residual (P_shifted - U*V)[z^k] and solving the
 * bivariate Diophantine; verify final product.
 *
 * Returns false if the bivariate Diophantine fails or the final
 * verification U*V == P does not hold.  On failure *U_out, *V_out are
 * NULL.  On success the caller owns them. */
bool mpoly_hensel_lift_3_2(const MPoly* P,
                            const MPoly* U_xy, const MPoly* V_xy,
                            int var_main, int var_y, int var_z,
                            int64_t alpha_z,
                            MPoly** U_out, MPoly** V_out);

#ifdef __cplusplus
}
#endif

#endif /* MVFACTOR3_H */
