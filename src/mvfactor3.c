/*
 * mvfactor3.c
 * -----------
 * Trivariate Hensel factoring on top of MPoly.  Phase F2 MVP scope:
 * monic-in-main two-factor lift for n = 3 polynomials.  The bivariate
 * Diophantine inside the lift is itself a Hensel iteration in the
 * secondary variable (var_y) on top of zupoly_diophantine.
 *
 * Layering:
 *   mpoly_hensel_lift_3_2   (this file)
 *     -> mpoly_diophantine_2
 *           -> zupoly_diophantine  (existing, src/zupoly.c)
 *           +  Hensel iteration in var_y
 *
 * See FACTOR_PLAN.md §F2 and SPEC.md for the broader picture.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "mvfactor3.h"
#include "mpoly.h"
#include "zupoly.h"

/* ---------------------------------------------------------------- */
/*  MPoly <-> ZUPoly conversions                                    */
/* ---------------------------------------------------------------- */

ZUPoly* mpoly_to_zupoly_in(const MPoly* p, int var_main) {
    if (var_main < 0 || var_main >= p->n_vars) return NULL;
    /* Validate: every term has exp[v] == 0 for v != var_main. */
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* row = p->exps + i * (size_t)p->n_vars;
        for (int v = 0; v < p->n_vars; v++) {
            if (v != var_main && row[v] != 0) return NULL;
        }
    }
    int d = mpoly_deg_var(p, var_main);
    if (d < 0) return zupoly_zero();
    ZUPoly* z = zupoly_new(d + 1);
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* row = p->exps + i * (size_t)p->n_vars;
        zupoly_setcoef(z, row[var_main], p->coefs[i]);
    }
    zupoly_normalize(z);
    return z;
}

MPoly* zupoly_to_mpoly_in(const ZUPoly* z, int n_vars, int var_main) {
    MPoly* p = mpoly_new(n_vars);
    if (!z || z->deg < 0) return p;
    int row[n_vars > 0 ? n_vars : 1];
    memset(row, 0, sizeof(int) * (size_t)n_vars);
    for (int k = 0; k <= z->deg; k++) {
        const mpz_t* c = zupoly_getcoef(z, k);
        if (!c || mpz_sgn(*c) == 0) continue;
        row[var_main] = k;
        mpoly_push_term(p, row, *c);
    }
    /* push_term appends in increasing var_main order; mpoly storage is
     * lex descending which means var_main first when var_main is at
     * index 0 (highest x_0 first).  To be safe across var_main choices,
     * normalise. */
    mpoly_normalize(p);
    return p;
}

BPoly* mpoly_to_bpoly_in(const MPoly* p, int var_main, int var_y) {
    if (var_main < 0 || var_main >= p->n_vars) return NULL;
    if (var_y    < 0 || var_y    >= p->n_vars) return NULL;
    if (var_main == var_y) return NULL;

    /* Validate: all other vars have exponent 0. */
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* row = p->exps + i * (size_t)p->n_vars;
        for (int v = 0; v < p->n_vars; v++) {
            if (v == var_main || v == var_y) continue;
            if (row[v] != 0) return NULL;
        }
    }

    int dx = mpoly_deg_var(p, var_main);
    if (dx < 0) return bpoly_zero();

    BPoly* b = bpoly_new(dx + 1);

    /* For each var_main-degree i, build the y-coefficient ZUPoly.  We
     * collect them as we walk p (terms come out in lex desc order
     * which ALMOST matches what we want; safest is to allocate
     * per-x-degree slots and fill in any order). */
    ZUPoly** slots = (ZUPoly**)calloc((size_t)(dx + 1), sizeof(ZUPoly*));
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* row = p->exps + t * (size_t)p->n_vars;
        int xi = row[var_main];
        int yi = row[var_y];
        if (!slots[xi]) slots[xi] = zupoly_new(yi + 1);
        zupoly_setcoef(slots[xi], yi, p->coefs[t]);
    }
    for (int i = 0; i <= dx; i++) {
        if (slots[i]) {
            zupoly_normalize(slots[i]);
            if (!zupoly_is_zero(slots[i])) {
                bpoly_set_xcoef(b, i, slots[i]);  /* takes ownership */
            } else {
                zupoly_free(slots[i]);
            }
        }
    }
    free(slots);
    return b;
}

MPoly* bpoly_to_mpoly_in(const BPoly* b, int n_vars,
                         int var_main, int var_y) {
    MPoly* p = mpoly_new(n_vars);
    if (!b || b->deg_x < 0) return p;
    if (var_main < 0 || var_main >= n_vars) return p;
    if (var_y    < 0 || var_y    >= n_vars) return p;

    int row[n_vars > 0 ? n_vars : 1];
    memset(row, 0, sizeof(int) * (size_t)n_vars);
    for (int i = 0; i <= b->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(b, i);
        if (!yi) continue;
        for (int k = 0; k <= yi->deg; k++) {
            const mpz_t* c = zupoly_getcoef(yi, k);
            if (!c || mpz_sgn(*c) == 0) continue;
            row[var_main] = i;
            row[var_y]    = k;
            mpoly_push_term(p, row, *c);
        }
    }
    /* Reset row to all-zero, but it's a stack array -- doesn't matter. */
    mpoly_normalize(p);
    return p;
}

/* ---------------------------------------------------------------- */
/*  Helpers for the bivariate Diophantine                           */
/* ---------------------------------------------------------------- */

/* Extract the y^k coefficient of an MPoly P, viewed as a polynomial
 * in (var_main, var_y).  Returns an MPoly with the same n_vars but
 * with var_y-exponent == 0 in every output term.  Other variables'
 * exponents are preserved -- but for our Diophantine use case the
 * caller has already guaranteed all exponents outside {var_main,
 * var_y} are 0. */
static MPoly* mpoly_yk_coef(const MPoly* P, int var_y, int k) {
    return mpoly_coef_of_var(P, var_y, k);
}

/* Multiply an MPoly by var_y^k: shift each term's var_y-exponent by k. */
static MPoly* mpoly_mul_yk(const MPoly* p, int var_y, int k) {
    if (k == 0 || mpoly_is_zero(p)) return mpoly_copy(p);
    MPoly* q = mpoly_copy(p);
    /* All terms get var_y-exponent shifted by k.  Since the lex order
     * is preserved (we add a constant to a single column), no re-sort
     * is needed.  But we must validate that lex order is actually
     * preserved -- it is, as long as all terms shift identically. */
    for (size_t i = 0; i < q->n_terms; i++) {
        int* row = q->exps + i * (size_t)q->n_vars;
        row[var_y] += k;
    }
    return q;
}

/* In-place: p += q.  Frees q.  Caller still owns p. */
static void mpoly_add_inplace(MPoly** p_ptr, MPoly* q) {
    MPoly* sum = mpoly_add(*p_ptr, q);
    mpoly_free(*p_ptr);
    mpoly_free(q);
    *p_ptr = sum;
}

/* ---------------------------------------------------------------- */
/*  Bivariate Diophantine                                           */
/* ---------------------------------------------------------------- */

/* Internal: try a specific alpha_y for the bivariate Diophantine. */
static bool mpoly_diophantine_2_at_alpha(const MPoly* U, const MPoly* V,
                                         const MPoly* E,
                                         int var_main, int var_y,
                                         int64_t alpha_y,
                                         MPoly** delta_u_out,
                                         MPoly** delta_v_out) {
    *delta_u_out = NULL;
    *delta_v_out = NULL;

    int n = U->n_vars;

    /* Specialise var_y = alpha_y to get univariate U(x), V(x), E(x). */
    MPoly* U_x = mpoly_subst_var_int(U, var_y, alpha_y);
    MPoly* V_x = mpoly_subst_var_int(V, var_y, alpha_y);
    MPoly* E_x = mpoly_subst_var_int(E, var_y, alpha_y);

    /* Convert to ZUPoly.  After substitution, every term's var_y-
     * exponent is 0 by construction.  All other non-main vars must
     * also be 0 -- precondition. */
    ZUPoly* u_zu = mpoly_to_zupoly_in(U_x, var_main);
    ZUPoly* v_zu = mpoly_to_zupoly_in(V_x, var_main);
    ZUPoly* e_zu = mpoly_to_zupoly_in(E_x, var_main);
    mpoly_free(U_x); mpoly_free(V_x); mpoly_free(E_x);

    if (!u_zu || !v_zu || !e_zu) {
        zupoly_free(u_zu); zupoly_free(v_zu); zupoly_free(e_zu);
        return false;
    }

    /* Univariate Diophantine. */
    ZUPoly* du0 = NULL; ZUPoly* dv0 = NULL;
    bool ok = zupoly_diophantine(u_zu, v_zu, e_zu, &du0, &dv0);
    zupoly_free(e_zu);
    if (!ok) {
        zupoly_free(u_zu); zupoly_free(v_zu);
        zupoly_free(du0); zupoly_free(dv0);
        return false;
    }
    zupoly_free(u_zu); zupoly_free(v_zu);

    /* Initial bivariate solution: var_y = alpha_y.  We work in shifted
     * coordinates y' := var_y - alpha_y so the lift starts at y' = 0. */
    MPoly* delta_u = zupoly_to_mpoly_in(du0, n, var_main);
    MPoly* delta_v = zupoly_to_mpoly_in(dv0, n, var_main);
    zupoly_free(du0); zupoly_free(dv0);

    /* Shift U, V, E by var_y -> var_y + alpha_y so the lift sees an
     * "image at var_y = 0" matching the univariate seed. */
    MPoly* U_sh = mpoly_shift_var_int(U, var_y, alpha_y);
    MPoly* V_sh = mpoly_shift_var_int(V, var_y, alpha_y);
    MPoly* E_sh = mpoly_shift_var_int(E, var_y, alpha_y);

    /* Lift in var_y up to deg_{var_y}(E_sh).  Higher degrees of E_sh
     * cannot be cancelled by lower-degree corrections, so this is the
     * right cap. */
    int B_y = mpoly_deg_var(E_sh, var_y);

    /* Maintain the partial product  R = delta_u * V_sh + delta_v * U_sh
     * mod var_y^{k}, expanded incrementally. */
    MPoly* R = mpoly_new(n);
    {
        /* R_0 (at y'=0) = delta_u * V_sh|y=0 + delta_v * U_sh|y=0
         *                = delta_u * V_x + delta_v * U_x = E_x.
         * But we don't actually need R_0 stored as a residual base;
         * the residual at iter k is E_sh[y^k] - cross-terms from prior
         * iterates.  Easier: re-derive at each iteration. */
        mpoly_free(R); R = NULL;
    }

    for (int k = 1; k <= B_y; k++) {
        /* Compute  prod_k = (delta_u * V_sh + delta_v * U_sh)[y^k]. */
        MPoly* duV = mpoly_mul(delta_u, V_sh);
        MPoly* dvU = mpoly_mul(delta_v, U_sh);
        MPoly* prod = mpoly_add(duV, dvU);
        mpoly_free(duV); mpoly_free(dvU);

        MPoly* prod_yk = mpoly_yk_coef(prod, var_y, k);
        mpoly_free(prod);

        MPoly* E_yk = mpoly_yk_coef(E_sh, var_y, k);
        MPoly* res_yk = mpoly_sub(E_yk, prod_yk);
        mpoly_free(E_yk); mpoly_free(prod_yk);

        if (mpoly_is_zero(res_yk)) {
            mpoly_free(res_yk);
            continue;
        }

        /* Solve univariate Diophantine for the residual at y^k. */
        ZUPoly* res_zu = mpoly_to_zupoly_in(res_yk, var_main);
        mpoly_free(res_yk);
        if (!res_zu) {
            mpoly_free(delta_u); mpoly_free(delta_v);
            mpoly_free(U_sh); mpoly_free(V_sh); mpoly_free(E_sh);
            return false;
        }

        ZUPoly* u0_zu = mpoly_to_zupoly_in(
            mpoly_yk_coef(U_sh, var_y, 0), var_main);
        ZUPoly* v0_zu = mpoly_to_zupoly_in(
            mpoly_yk_coef(V_sh, var_y, 0), var_main);

        ZUPoly* du_k = NULL; ZUPoly* dv_k = NULL;
        bool ok_k = u0_zu && v0_zu &&
                    zupoly_diophantine(u0_zu, v0_zu, res_zu, &du_k, &dv_k);
        zupoly_free(u0_zu); zupoly_free(v0_zu); zupoly_free(res_zu);

        if (!ok_k) {
            zupoly_free(du_k); zupoly_free(dv_k);
            mpoly_free(delta_u); mpoly_free(delta_v);
            mpoly_free(U_sh); mpoly_free(V_sh); mpoly_free(E_sh);
            return false;
        }

        /* delta_u += y^k * du_k, delta_v += y^k * dv_k. */
        MPoly* du_k_mp = zupoly_to_mpoly_in(du_k, n, var_main);
        MPoly* dv_k_mp = zupoly_to_mpoly_in(dv_k, n, var_main);
        zupoly_free(du_k); zupoly_free(dv_k);

        MPoly* du_k_yk = mpoly_mul_yk(du_k_mp, var_y, k);
        MPoly* dv_k_yk = mpoly_mul_yk(dv_k_mp, var_y, k);
        mpoly_free(du_k_mp); mpoly_free(dv_k_mp);

        mpoly_add_inplace(&delta_u, du_k_yk);
        mpoly_add_inplace(&delta_v, dv_k_yk);
    }

    /* Final verification: delta_u * V_sh + delta_v * U_sh == E_sh. */
    MPoly* duV = mpoly_mul(delta_u, V_sh);
    MPoly* dvU = mpoly_mul(delta_v, U_sh);
    MPoly* sum = mpoly_add(duV, dvU);
    mpoly_free(duV); mpoly_free(dvU);

    bool match = mpoly_eq(sum, E_sh);
    mpoly_free(sum);
    mpoly_free(U_sh); mpoly_free(V_sh); mpoly_free(E_sh);

    if (!match) {
        mpoly_free(delta_u); mpoly_free(delta_v);
        return false;
    }

    /* Unshift delta_u, delta_v back to original var_y coordinate. */
    MPoly* delta_u_un = mpoly_shift_var_int(delta_u, var_y, -alpha_y);
    MPoly* delta_v_un = mpoly_shift_var_int(delta_v, var_y, -alpha_y);
    mpoly_free(delta_u); mpoly_free(delta_v);

    *delta_u_out = delta_u_un;
    *delta_v_out = delta_v_un;
    return true;
}

bool mpoly_diophantine_2(const MPoly* U, const MPoly* V, const MPoly* E,
                         int var_main, int var_y,
                         MPoly** delta_u_out, MPoly** delta_v_out) {
    /* Try several alpha_y values until one works.  We need:
     *   - U|y=alpha_y has the same var_main-degree as U.
     *   - V|y=alpha_y has the same var_main-degree as V.
     *   - U|y=alpha_y and V|y=alpha_y are coprime over Q[var_main].
     *
     * Heuristic order: 0, 1, -1, 2, -2, ...  Cap at 11 attempts. */
    static const int64_t alpha_order[] = {
        0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5
    };
    enum { ALPHA_COUNT = sizeof(alpha_order) / sizeof(alpha_order[0]) };

    int U_dx = mpoly_deg_var(U, var_main);
    int V_dx = mpoly_deg_var(V, var_main);

    for (size_t i = 0; i < ALPHA_COUNT; i++) {
        int64_t alpha = alpha_order[i];
        MPoly* U_x = mpoly_subst_var_int(U, var_y, alpha);
        MPoly* V_x = mpoly_subst_var_int(V, var_y, alpha);
        bool deg_ok = (mpoly_deg_var(U_x, var_main) == U_dx) &&
                      (mpoly_deg_var(V_x, var_main) == V_dx);
        mpoly_free(U_x); mpoly_free(V_x);
        if (!deg_ok) continue;

        if (mpoly_diophantine_2_at_alpha(U, V, E, var_main, var_y, alpha,
                                         delta_u_out, delta_v_out)) {
            return true;
        }
    }
    return false;
}

/* ---------------------------------------------------------------- */
/*  Trivariate two-factor Hensel lift                               */
/* ---------------------------------------------------------------- */

bool mpoly_hensel_lift_3_2(const MPoly* P,
                            const MPoly* U_xy, const MPoly* V_xy,
                            int var_main, int var_y, int var_z,
                            int64_t alpha_z,
                            MPoly** U_out, MPoly** V_out) {
    *U_out = NULL;
    *V_out = NULL;

    if (!P || mpoly_is_zero(P)) return false;
    if (!U_xy || !V_xy || mpoly_is_zero(U_xy) || mpoly_is_zero(V_xy))
        return false;
    if (var_main < 0 || var_main >= P->n_vars) return false;
    if (var_y < 0 || var_y >= P->n_vars) return false;
    if (var_z < 0 || var_z >= P->n_vars) return false;
    if (var_main == var_y || var_main == var_z || var_y == var_z) return false;

    /* Sanity: U_xy and V_xy must have z-degree 0 (constant in var_z). */
    if (!mpoly_is_constant_in_var(U_xy, var_z)) return false;
    if (!mpoly_is_constant_in_var(V_xy, var_z)) return false;

    int n = P->n_vars;

    /* Shift var_z -> var_z + alpha_z so the lift starts at z = 0. */
    MPoly* P_sh = mpoly_shift_var_int(P, var_z, alpha_z);

    /* Confirm preconditions at z = 0: P_sh|z=0 == U_xy * V_xy. */
    MPoly* P_z0 = mpoly_subst_var_int(P_sh, var_z, 0);
    MPoly* check = mpoly_mul(U_xy, V_xy);
    bool seed_ok = mpoly_eq(P_z0, check);
    mpoly_free(P_z0); mpoly_free(check);
    if (!seed_ok) {
        mpoly_free(P_sh);
        return false;
    }

    /* Initialise U, V as the bivariate seeds (z-degree 0). */
    MPoly* U = mpoly_copy(U_xy);
    MPoly* V = mpoly_copy(V_xy);

    int B_z = mpoly_deg_var(P_sh, var_z);
    if (B_z < 0) B_z = 0;

    for (int k = 1; k <= B_z; k++) {
        /* E_k = (P_sh - U*V)[z^k]. */
        MPoly* prod = mpoly_mul(U, V);
        MPoly* diff = mpoly_sub(P_sh, prod);
        mpoly_free(prod);
        MPoly* E_k = mpoly_yk_coef(diff, var_z, k);
        mpoly_free(diff);

        if (mpoly_is_zero(E_k)) {
            mpoly_free(E_k);
            continue;
        }

        /* E_k is a polynomial in (var_main, var_y) (since we extracted
         * the var_z^k slice; the input P is in n vars, but at this
         * point var_z is gone for E_k).  Solve the bivariate Diophantine. */
        MPoly* delta_u = NULL; MPoly* delta_v = NULL;
        bool ok = mpoly_diophantine_2(U_xy, V_xy, E_k, var_main, var_y,
                                      &delta_u, &delta_v);
        mpoly_free(E_k);
        if (!ok) {
            mpoly_free(U); mpoly_free(V); mpoly_free(P_sh);
            return false;
        }

        /* U += z^k * delta_u, V += z^k * delta_v. */
        MPoly* du_zk = mpoly_mul_yk(delta_u, var_z, k);
        MPoly* dv_zk = mpoly_mul_yk(delta_v, var_z, k);
        mpoly_free(delta_u); mpoly_free(delta_v);

        mpoly_add_inplace(&U, du_zk);
        mpoly_add_inplace(&V, dv_zk);
    }

    /* Final verification. */
    MPoly* product = mpoly_mul(U, V);
    bool match = mpoly_eq(product, P_sh);
    mpoly_free(product);
    mpoly_free(P_sh);

    if (!match) {
        mpoly_free(U); mpoly_free(V);
        return false;
    }

    /* Unshift var_z back. */
    MPoly* U_un = mpoly_shift_var_int(U, var_z, -alpha_z);
    MPoly* V_un = mpoly_shift_var_int(V, var_z, -alpha_z);
    mpoly_free(U); mpoly_free(V);

    *U_out = U_un;
    *V_out = V_un;
    (void)n;  /* currently unused beyond validation */
    return true;
}
