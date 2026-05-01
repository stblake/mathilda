/*
 * mvfactor.c
 * ----------
 * Multivariate factoring orchestration.  Currently implements the
 * bivariate two-factor Hensel lift; additional phases (multifactor
 * lift, recombination, n-variate dispatch) will live here as they
 * land.  See FACTOR_PLAN.md for the broader picture.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <gmp.h>

#include "mvfactor.h"
#include "bpoly.h"
#include "zupoly.h"

/* ====================================================================== */
/*  Bivariate two-factor Hensel lift                                      */
/* ====================================================================== */

/* Promote a univariate ZUPoly u(x) to a BPoly that is u(x) viewed as
 * a bivariate polynomial of y-degree 0 in (x, y).  Each x-coefficient
 * of the result is the constant ZUPoly equal to the corresponding
 * coefficient of u. */
static BPoly* zupoly_to_bpoly_const(const ZUPoly* u) {
    if (u->deg < 0) return bpoly_zero();
    BPoly* b = bpoly_new(u->deg + 1);
    for (int i = 0; i <= u->deg; i++) {
        const mpz_t* c = zupoly_getcoef(u, i);
        if (c && mpz_sgn(*c) != 0) {
            ZUPoly* y_const = zupoly_new(1);
            zupoly_setcoef(y_const, 0, *c);
            bpoly_set_xcoef(b, i, y_const);
        }
    }
    return b;
}

/* Extract the "y^k coefficient" of a BPoly D, viewed as a polynomial
 * in x.  Concretely: returns a fresh ZUPoly whose i-th coefficient is
 * the coefficient of y^k in D.cx[i].  Used during Hensel iteration to
 * isolate the residual term that the next correction must zero out. */
static ZUPoly* bpoly_extract_y_coef(const BPoly* D, int k) {
    if (D->deg_x < 0 || k < 0) return zupoly_zero();
    ZUPoly* r = zupoly_new(D->deg_x + 1);
    for (int i = 0; i <= D->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(D, i);
        if (!yi) continue;
        const mpz_t* cik = zupoly_getcoef(yi, k);
        if (cik && mpz_sgn(*cik) != 0) {
            zupoly_setcoef(r, i, *cik);
        }
    }
    return r;
}

/* Add (in-place) y^k * delta to the BPoly U, where `delta` is a
 * polynomial in x.  Concretely: for each non-zero c_i in delta,
 * augment U.cx[i] by c_i * y^k.  Updates U->deg_x as needed.
 *
 * `delta` is borrowed; the caller retains ownership. */
static void bpoly_add_yk_delta(BPoly* U, int k, const ZUPoly* delta) {
    if (delta->deg < 0) return;
    for (int i = 0; i <= delta->deg; i++) {
        const mpz_t* di = zupoly_getcoef(delta, i);
        if (!di || mpz_sgn(*di) == 0) continue;
        /* Need to update U.cx[i] by adding di * y^k.  If U.cx[i] does
         * not yet exist, create it. */
        const ZUPoly* existing = bpoly_get_xcoef(U, i);
        ZUPoly* updated;
        if (!existing) {
            updated = zupoly_new(k + 1);
            zupoly_setcoef(updated, k, *di);
        } else {
            updated = zupoly_copy(existing);
            /* updated.c[k] += di. */
            const mpz_t* old_k = zupoly_getcoef(updated, k);
            mpz_t sum; mpz_init(sum);
            if (old_k) mpz_add(sum, *old_k, *di);
            else       mpz_set(sum, *di);
            zupoly_setcoef(updated, k, sum);
            mpz_clear(sum);
        }
        /* bpoly_set_xcoef takes ownership of `updated`. */
        bpoly_set_xcoef(U, i, updated);
    }
}

/* Multiply a univariate ZUPoly u by a constant integer k in place,
 * producing a fresh ZUPoly.  k = 1 yields a copy. */
static ZUPoly* zupoly_scale_by_const(const ZUPoly* u, const mpz_t k) {
    if (mpz_cmp_ui(k, 1) == 0) return zupoly_copy(u);
    return zupoly_scale(u, k);
}

/* F3c helper: compute a Mignotte-style upper bound on the integer
 * coefficients of any divisor F of P viewed in Z[x, y].  We use the
 * conservative form
 *
 *     M = 2^{deg_x(P) + deg_y(P)} * ceil(||P||_2)
 *
 * where ||P||_2 is the L2 norm of P's integer coefficient vector.
 * Mignotte's theorem (e.g., Geddes-Czapor-Labahn ch. 4) gives
 * ||F||_inf <= 2^{deg(F)} * ||P||_2 in the univariate case; the same
 * bound applies coefficient-wise in the multivariate case using the
 * total degree of F, which is at most the total degree of P.  Output
 * is initialised inside the function; caller must mpz_clear it. */
static void bpoly_mignotte_bound(const BPoly* P, mpz_t M_out) {
    mpz_t sumsq, sq, l2;
    mpz_init_set_ui(sumsq, 0);
    mpz_init(sq);
    mpz_init(l2);

    for (int i = 0; i <= P->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(P, i);
        if (!yi) continue;
        for (int kk = 0; kk <= yi->deg; kk++) {
            const mpz_t* cik = zupoly_getcoef(yi, kk);
            if (!cik) continue;
            mpz_mul(sq, *cik, *cik);
            mpz_add(sumsq, sumsq, sq);
        }
    }

    mpz_sqrt(l2, sumsq);
    mpz_add_ui(l2, l2, 1);  /* ceil */

    int d_total = P->deg_x + bpoly_deg_y(P);
    if (d_total < 0) d_total = 0;

    mpz_init(M_out);
    mpz_mul_2exp(M_out, l2, (unsigned long)d_total);

    mpz_clear(sumsq);
    mpz_clear(sq);
    mpz_clear(l2);
}

/* F3c helper: return true if any integer coefficient of U exceeds
 * `bound` in absolute value.  Short-circuits on the first violation. */
static bool bpoly_max_abs_coef_exceeds(const BPoly* U, const mpz_t bound) {
    for (int i = 0; i <= U->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(U, i);
        if (!yi) continue;
        for (int kk = 0; kk <= yi->deg; kk++) {
            const mpz_t* cik = zupoly_getcoef(yi, kk);
            if (!cik) continue;
            if (mpz_cmpabs(*cik, bound) > 0) return true;
        }
    }
    return false;
}

/* F3b helper: compute the y^k coefficient of U_{k-1} * V_{k-1} as a
 * polynomial in x.  At the start of iteration k, U has y-degree at
 * most k-1 and V likewise; the y^k coefficient of their product is
 *
 *     sum_{a=1}^{k-1} U[y^a] * V[y^{k-a}]
 *
 * (a = 0 and a = k are skipped because U[y^k] and V[y^k] are still 0
 * at iteration k's start).  Returns a fresh ZUPoly; caller frees. */
static ZUPoly* bpoly_uv_cross_yk(const BPoly* U, const BPoly* V, int k) {
    ZUPoly* acc = zupoly_zero();
    for (int a = 1; a < k; a++) {
        int b = k - a;
        ZUPoly* Ua = bpoly_extract_y_coef(U, a);
        if (zupoly_is_zero(Ua)) { zupoly_free(Ua); continue; }
        ZUPoly* Vb = bpoly_extract_y_coef(V, b);
        if (zupoly_is_zero(Vb)) {
            zupoly_free(Ua); zupoly_free(Vb);
            continue;
        }
        ZUPoly* prod = zupoly_mul(Ua, Vb);
        zupoly_free(Ua);
        zupoly_free(Vb);
        ZUPoly* sum = zupoly_add(acc, prod);
        zupoly_free(acc);
        zupoly_free(prod);
        acc = sum;
    }
    return acc;
}

bool bpoly_hensel_lift_2(const BPoly* P, const ZUPoly* u, const ZUPoly* v,
                        BPoly** U_out, BPoly** V_out) {
    *U_out = NULL;
    *V_out = NULL;

    /* Sanity: P must be non-zero, u and v non-zero, and u*v must
     * equal P(x, 0).  We don't verify u*v == P(x,0) here because the
     * lift loop will fail to converge if the preconditions are wrong;
     * caller is expected to set them up via bz_factor_to_expr. */
    if (P->deg_x < 0) return false;
    if (u->deg < 0 || v->deg < 0) return false;

    /* Initialise U, V as constant-in-y bivariate polynomials. */
    BPoly* U = zupoly_to_bpoly_const(u);
    BPoly* V = zupoly_to_bpoly_const(v);

    /* F3b: maintain UV = U*V mod y^{k+1} incrementally across iterations.
     * Initial state: UV_0 = u * v as a y-constant BPoly (y^0 only). */
    ZUPoly* uv0 = zupoly_mul(u, v);
    BPoly*  UV  = zupoly_to_bpoly_const(uv0);
    zupoly_free(uv0);

    /* F3c: Mignotte coefficient bound.  Any true integer factor of P
     * has all coefficients bounded by this; if U or V develops a
     * coefficient larger than this during lifting, the subset cannot
     * produce a valid integer factor and we abort. */
    mpz_t mignotte_M;
    bpoly_mignotte_bound(P, mignotte_M);

    int B = bpoly_deg_y(P);  /* lift up to y^B inclusive */
    if (B < 0) B = 0;

    for (int k = 1; k <= B; k++) {
        /* Step 1: extend UV from "U_{k-1}*V_{k-1} mod y^k" to
         *         "U_{k-1}*V_{k-1} mod y^{k+1}" by adding the y^k cross
         *         coefficient.  For k = 1 the cross sum is empty. */
        ZUPoly* cross = bpoly_uv_cross_yk(U, V, k);
        if (!zupoly_is_zero(cross)) bpoly_add_yk_delta(UV, k, cross);
        zupoly_free(cross);

        /* Step 2: residual at y^k is (P[y^k] - UV[y^k]). */
        ZUPoly* P_yk  = bpoly_extract_y_coef(P,  k);
        ZUPoly* UV_yk = bpoly_extract_y_coef(UV, k);
        ZUPoly* E_k   = zupoly_sub(P_yk, UV_yk);
        zupoly_free(P_yk);
        zupoly_free(UV_yk);

        if (zupoly_is_zero(E_k)) {
            zupoly_free(E_k);
            continue;
        }

        /* Step 3: solve  delta_u * v + delta_v * u = E_k  (u, v monic
         *         and coprime). */
        ZUPoly *delta_u = NULL, *delta_v = NULL;
        bool ok = zupoly_diophantine(u, v, E_k, &delta_u, &delta_v);
        zupoly_free(E_k);
        if (!ok) {
            bpoly_free(U); bpoly_free(V); bpoly_free(UV);
            mpz_clear(mignotte_M);
            return false;
        }

        /* Step 4: U += y^k * delta_u, V += y^k * delta_v. */
        bpoly_add_yk_delta(U, k, delta_u);
        bpoly_add_yk_delta(V, k, delta_v);

        /* Step 5: UV[y^k] += delta_u * V_0 + U_0 * delta_v
         *         (since V_{k-1}[y^0] = V_0 = v and U_{k-1}[y^0] = u). */
        ZUPoly* du_v = zupoly_mul(delta_u, v);
        ZUPoly* u_dv = zupoly_mul(u, delta_v);
        ZUPoly* corr = zupoly_add(du_v, u_dv);
        zupoly_free(du_v);
        zupoly_free(u_dv);
        if (!zupoly_is_zero(corr)) bpoly_add_yk_delta(UV, k, corr);
        zupoly_free(corr);

        zupoly_free(delta_u);
        zupoly_free(delta_v);

        /* F3c fast-fail: if any U or V coefficient exceeds the
         * Mignotte bound, this subset can't produce an integer factor.
         * Abort the lift early. */
        if (bpoly_max_abs_coef_exceeds(U, mignotte_M) ||
            bpoly_max_abs_coef_exceeds(V, mignotte_M)) {
            bpoly_free(U); bpoly_free(V); bpoly_free(UV);
            mpz_clear(mignotte_M);
            return false;
        }
    }

    mpz_clear(mignotte_M);
    bpoly_free(UV);

    /* Final verification: full U*V must equal P exactly.  (The
     * incremental UV only tracks mod y^{B+1}; high-y-degree
     * cancellation must still be checked.) */
    BPoly* product = bpoly_mul(U, V);
    bool match = bpoly_eq(product, P);
    bpoly_free(product);

    if (!match) {
        bpoly_free(U); bpoly_free(V);
        return false;
    }

    *U_out = U;
    *V_out = V;
    return true;
}

/* ====================================================================== */
/*  Phase F1 Stage 3: predicted-LC two-factor Hensel lift                  */
/*                                                                        */
/*  Wang's leading-coefficient correction for inputs whose lc_x(P)(y) is  */
/*  a non-constant polynomial in y.  Each lifted factor has its leading-  */
/*  x coefficient PINNED to a predetermined polynomial in y (q_u or q_v)  */
/*  whose product equals lc_x(P).                                         */
/*                                                                        */
/*  At each iteration k the correction Δu has the form                    */
/*    Δu = qu_k * x^{deg(u)} + δu_low   (deg(δu_low) < deg(u))           */
/*  -- the leading-x coefficient is forced to be the y^k coefficient of   */
/*  q_u, and only the lower x-degree part is determined by Diophantine.   */
/*  This keeps lc_x(U)(y) growing as exactly q_u(y), and similarly for V. */
/*                                                                        */
/*  The Diophantine receives an "adjusted" residual that subtracts off    */
/*  the contribution of the pinned leading deltas:                        */
/*    E_k_adj = E_k - (qu_k · x^{d_u}) · v - (qv_k · x^{d_v}) · u         */
/*  and is then solved over Z[x] with monic u, v as in bpoly_hensel_lift_2 */
/* ====================================================================== */

bool bpoly_hensel_lift_2_lc(const BPoly* P,
                             const ZUPoly* u, const ZUPoly* v,
                             const ZUPoly* q_u, const ZUPoly* q_v,
                             BPoly** U_out, BPoly** V_out) {
    *U_out = NULL;
    *V_out = NULL;

    if (!P || P->deg_x < 0) return false;
    if (!u || !v || u->deg < 0 || v->deg < 0) return false;
    if (!q_u || !q_v || q_u->deg < 0 || q_v->deg < 0) return false;

    /* The Diophantine step requires monic u and v. */
    if (mpz_cmp_ui(u->c[u->deg], 1) != 0) return false;
    if (mpz_cmp_ui(v->c[v->deg], 1) != 0) return false;

    /* Initial values absorb the constant terms of the predicted LCs:
     *   U_0(x) = q_u(0) * u(x),  V_0(x) = q_v(0) * v(x).
     * This makes lc_x(U)(0) = q_u(0) and lc_x(V)(0) = q_v(0) match the
     * predicted polynomials' constant terms. */
    const mpz_t* qu_0 = zupoly_getcoef(q_u, 0);
    const mpz_t* qv_0 = zupoly_getcoef(q_v, 0);
    if (!qu_0 || !qv_0) return false;
    if (mpz_sgn(*qu_0) == 0 || mpz_sgn(*qv_0) == 0) {
        /* Predicted LC vanishes at y = 0 -- caller must shift first. */
        return false;
    }

    ZUPoly* u_init = zupoly_scale_by_const(u, *qu_0);
    ZUPoly* v_init = zupoly_scale_by_const(v, *qv_0);
    BPoly* U = zupoly_to_bpoly_const(u_init);
    BPoly* V = zupoly_to_bpoly_const(v_init);

    /* F3b: incremental UV = U*V mod y^{k+1}, kept in lockstep with U/V.
     * Initial UV_0 = (qu_0 * u) * (qv_0 * v) lives at y^0. */
    ZUPoly* uv0 = zupoly_mul(u_init, v_init);
    BPoly*  UV  = zupoly_to_bpoly_const(uv0);
    zupoly_free(uv0);

    /* F3c: Mignotte coefficient bound for fast-fail. */
    mpz_t mignotte_M;
    bpoly_mignotte_bound(P, mignotte_M);

    int B = bpoly_deg_y(P);
    if (B < 0) B = 0;

    int d_u = u->deg;
    int d_v = v->deg;

    for (int k = 1; k <= B; k++) {
        /* Extend UV by the y^k coefficient of U_{k-1}*V_{k-1}. */
        ZUPoly* cross = bpoly_uv_cross_yk(U, V, k);
        if (!zupoly_is_zero(cross)) bpoly_add_yk_delta(UV, k, cross);
        zupoly_free(cross);

        /* Residual at y^k. */
        ZUPoly* P_yk  = bpoly_extract_y_coef(P,  k);
        ZUPoly* UV_yk = bpoly_extract_y_coef(UV, k);
        ZUPoly* E_k   = zupoly_sub(P_yk, UV_yk);
        zupoly_free(P_yk);
        zupoly_free(UV_yk);

        /* Pin the leading x-coefficients of Δu and Δv to qu_k and qv_k. */
        const mpz_t* qu_k_p = (k <= q_u->deg) ? zupoly_getcoef(q_u, k) : NULL;
        const mpz_t* qv_k_p = (k <= q_v->deg) ? zupoly_getcoef(q_v, k) : NULL;
        bool qu_k_zero = !qu_k_p || mpz_sgn(*qu_k_p) == 0;
        bool qv_k_zero = !qv_k_p || mpz_sgn(*qv_k_p) == 0;

        ZUPoly* leading_du = zupoly_zero();
        ZUPoly* leading_dv = zupoly_zero();
        if (!qu_k_zero) zupoly_setcoef(leading_du, d_u, *qu_k_p);
        if (!qv_k_zero) zupoly_setcoef(leading_dv, d_v, *qv_k_p);

        /* E_k_adj = E_k - leading_du * v - leading_dv * u. */
        ZUPoly* lu_v = zupoly_mul(leading_du, v);
        ZUPoly* lv_u = zupoly_mul(leading_dv, u);
        ZUPoly* lead_contrib = zupoly_add(lu_v, lv_u);
        zupoly_free(lu_v);
        zupoly_free(lv_u);
        ZUPoly* E_k_adj = zupoly_sub(E_k, lead_contrib);
        zupoly_free(E_k);
        zupoly_free(lead_contrib);

        /* Solve δu * v + δv * u = E_k_adj  with deg(δu) < d_u,
         * deg(δv) < d_v. */
        ZUPoly *delta_u_low = NULL, *delta_v_low = NULL;
        bool ok;
        if (zupoly_is_zero(E_k_adj)) {
            delta_u_low = zupoly_zero();
            delta_v_low = zupoly_zero();
            ok = true;
        } else {
            ok = zupoly_diophantine(u, v, E_k_adj,
                                    &delta_u_low, &delta_v_low);
        }
        zupoly_free(E_k_adj);
        if (!ok) {
            zupoly_free(leading_du);
            zupoly_free(leading_dv);
            zupoly_free(u_init); zupoly_free(v_init);
            bpoly_free(U); bpoly_free(V); bpoly_free(UV);
            mpz_clear(mignotte_M);
            return false;
        }

        /* Δu = leading_du + δu_low (non-overlapping x-degrees: leading
         * lives at x^{d_u}, δu_low has degree < d_u). */
        ZUPoly* delta_u = zupoly_add(leading_du, delta_u_low);
        ZUPoly* delta_v = zupoly_add(leading_dv, delta_v_low);
        zupoly_free(leading_du);
        zupoly_free(leading_dv);
        zupoly_free(delta_u_low);
        zupoly_free(delta_v_low);

        bpoly_add_yk_delta(U, k, delta_u);
        bpoly_add_yk_delta(V, k, delta_v);

        /* F3b: UV[y^k] += delta_u * V_0 + U_0 * delta_v.  V_0 = v_init,
         *      U_0 = u_init (the qu_0/qv_0-scaled univariate factors). */
        ZUPoly* du_v0 = zupoly_mul(delta_u, v_init);
        ZUPoly* u0_dv = zupoly_mul(u_init, delta_v);
        ZUPoly* corr  = zupoly_add(du_v0, u0_dv);
        zupoly_free(du_v0);
        zupoly_free(u0_dv);
        if (!zupoly_is_zero(corr)) bpoly_add_yk_delta(UV, k, corr);
        zupoly_free(corr);

        zupoly_free(delta_u);
        zupoly_free(delta_v);

        /* F3c fast-fail: abort if U or V exceeds the Mignotte bound. */
        if (bpoly_max_abs_coef_exceeds(U, mignotte_M) ||
            bpoly_max_abs_coef_exceeds(V, mignotte_M)) {
            zupoly_free(u_init); zupoly_free(v_init);
            bpoly_free(U); bpoly_free(V); bpoly_free(UV);
            mpz_clear(mignotte_M);
            return false;
        }
    }

    mpz_clear(mignotte_M);
    zupoly_free(u_init);
    zupoly_free(v_init);
    bpoly_free(UV);

    /* Final verification. */
    BPoly* product = bpoly_mul(U, V);
    bool match = bpoly_eq(product, P);
    bpoly_free(product);

    if (!match) {
        bpoly_free(U); bpoly_free(V);
        return false;
    }

    *U_out = U;
    *V_out = V;
    return true;
}

/* ====================================================================== */
/*  Multifactor bivariate Hensel lift with subset recombination           */
/* ====================================================================== */

/* Compute prod(us[i]) over the set bits of `mask` and the cleared bits
 * of `mask` separately.  On return *prod_in points at the product of
 * us[i] for i where bit i of mask is set; *prod_out for cleared bits.
 * The caller must zupoly_free both. */
static void mvfactor_partition_product(const ZUPoly** us, int r,
                                       uint64_t mask,
                                       ZUPoly** prod_in,
                                       ZUPoly** prod_out) {
    ZUPoly *pin = NULL, *pout = NULL;
    for (int i = 0; i < r; i++) {
        bool in_mask = ((mask >> i) & 1u) != 0;
        ZUPoly** target = in_mask ? &pin : &pout;
        if (*target == NULL) {
            *target = zupoly_copy(us[i]);
        } else {
            ZUPoly* nxt = zupoly_mul(*target, us[i]);
            zupoly_free(*target);
            *target = nxt;
        }
    }
    *prod_in  = pin;
    *prod_out = pout;
}

/* Iterate over k-subsets of {0..r-1} via Gosper's hack: starting mask
 * has the k lowest bits set; next_combination advances to the next
 * k-combination in lexicographic order, wrapping past the last when
 * mask >= (1 << r).  Returns 0 when no more k-subsets remain.
 *
 * To avoid trying both halves of every partition (since {S, complement}
 * and {complement, S} produce the same lift via mvfactor_partition_product
 * just with U_in and U_out swapped), we restrict iteration to k-subsets
 * with k <= r/2, and additionally for k == r/2 (only when r is even)
 * we restrict to subsets whose lowest set bit is 0 -- breaking the
 * symmetry. */
static uint64_t gosper_next(uint64_t x) {
    /* Standard Gosper increment for "next number with same popcount". */
    uint64_t c = x & (uint64_t)(-(int64_t)x);
    uint64_t r_ = x + c;
    return (((r_ ^ x) >> 2) / c) | r_;
}

/* Forward declaration. */
static bool lift_multi_internal(const BPoly* P, const ZUPoly** us, int r,
                                BPoly*** Us_out, int* r_out);

/* F3a: sort the borrowed us[] array ascending by univariate x-degree.
 * Stable insertion sort -- r is small (typically <= 10) so the O(r^2)
 * cost is negligible compared to a single Hensel lift.  Putting the
 * lowest-degree factors first means the singleton k = 1 trials at the
 * head of recombination are the most likely to succeed (the true
 * bivariate factor of an input like (1 - x^12)(x - y^13) has the
 * univariate image `x` of degree 1, so subset {x} succeeds on the
 * first attempt instead of the 7th). */
static void mvfactor_sort_us_by_degree(const ZUPoly** us, int r) {
    for (int i = 1; i < r; i++) {
        const ZUPoly* key = us[i];
        int key_deg = key->deg;
        int j = i - 1;
        while (j >= 0 && us[j]->deg > key_deg) {
            us[j + 1] = us[j];
            j--;
        }
        us[j + 1] = key;
    }
}

bool bpoly_hensel_lift_multi(const BPoly* P, const ZUPoly** us, int r,
                             BPoly*** Us_out, int* r_out) {
    *Us_out = NULL;
    *r_out = 0;
    if (r <= 0) return false;

    /* Sort a private copy of the borrowed array.  The complement-side
     * recursion preserves relative order (comp_us iterates cleared bits
     * in increasing index), so a single sort at entry suffices. */
    const ZUPoly** sorted = (const ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)r);
    if (!sorted) return false;
    for (int i = 0; i < r; i++) sorted[i] = us[i];
    mvfactor_sort_us_by_degree(sorted, r);

    bool ok = lift_multi_internal(P, sorted, r, Us_out, r_out);
    free((void*)sorted);
    return ok;
}

static bool lift_multi_internal(const BPoly* P, const ZUPoly** us, int r,
                                BPoly*** Us_out, int* r_out) {
    if (r <= 0) return false;
    if (r == 1) {
        /* Trivial: there's only one factor, and it must equal P. */
        BPoly** result = (BPoly**)malloc(sizeof(BPoly*));
        if (!result) return false;
        result[0] = bpoly_copy(P);
        *Us_out = result;
        *r_out = 1;
        return true;
    }

    /* Try k-subsets in increasing size 1, 2, ..., r/2.
     *
     * For each k-subset S we attempt a two-factor lift_2 with
     *   uS = prod(us[i] : i in S),
     *   uC = prod(us[i] : i not in S).
     * Both products are monic (each us[i] is monic) and their gcd
     * is 1 (the us[i] are pairwise coprime by precondition).
     *
     * On success we found a true bivariate factor U_S and a bivariate
     * residual V_C with U_S * V_C = P.  Recurse on V_C with the
     * complement-side u's.
     *
     * On failure (no subset of size k yields a successful lift), try
     * the next size.  When k > r/2 we've exhausted all distinct
     * partitions -- return P as the sole factor (irreducible). */
    for (int k = 1; k * 2 <= r; k++) {
        /* Starting mask: k lowest bits set. */
        uint64_t mask = ((uint64_t)1 << k) - 1u;
        uint64_t end  = (r >= 64) ? UINT64_MAX : ((uint64_t)1 << r);

        while (mask < end) {
            /* Symmetry-break for k = r/2 with even r: skip subsets
             * whose lowest set bit is NOT 0, since {S, S^c} is the
             * same partition as {S^c, S}. */
            bool skip = false;
            if (k * 2 == r) {
                if ((mask & 1u) == 0) skip = true;
            }

            if (!skip) {
                ZUPoly *uS = NULL, *uC = NULL;
                mvfactor_partition_product(us, r, mask, &uS, &uC);
                /* uS = subset product (mask=1 bits); uC = complement. */
                if (!uS || !uC || uS->deg < 1 || uC->deg < 1) {
                    /* Should not happen with k in [1, r-1]. */
                    zupoly_free(uS); zupoly_free(uC);
                } else {
                                    BPoly *U_part = NULL, *V_part = NULL;
                    bool ok = bpoly_hensel_lift_2(P, uS, uC, &U_part, &V_part);
                                    zupoly_free(uS); zupoly_free(uC);

                    if (ok) {
                        /* Found one true factor.  Recurse on V_part
                         * with the complement-side u's.  Build that
                         * sub-array. */
                        const ZUPoly** comp_us = NULL;
                        int comp_count = r - k;
                        if (comp_count > 0) {
                            comp_us = (const ZUPoly**)malloc(
                                sizeof(ZUPoly*) * (size_t)comp_count);
                            if (!comp_us) {
                                bpoly_free(U_part); bpoly_free(V_part);
                                return false;
                            }
                            int idx = 0;
                            for (int i = 0; i < r; i++) {
                                if (((mask >> i) & 1u) == 0) {
                                    comp_us[idx++] = us[i];
                                }
                            }
                        }

                        BPoly** rest_results = NULL;
                        int rest_count = 0;
                        bool rec_ok = lift_multi_internal(
                            V_part, comp_us, comp_count,
                            &rest_results, &rest_count);
                        free((void*)comp_us);
                        bpoly_free(V_part);

                        if (!rec_ok) {
                            bpoly_free(U_part);
                            return false;
                        }

                        /* Assemble final array [U_part, rest_results...]. */
                        BPoly** result = (BPoly**)malloc(
                            sizeof(BPoly*) * (size_t)(1 + rest_count));
                        if (!result) {
                            bpoly_free(U_part);
                            for (int i = 0; i < rest_count; i++)
                                bpoly_free(rest_results[i]);
                            free(rest_results);
                            return false;
                        }
                        result[0] = U_part;
                        for (int i = 0; i < rest_count; i++) {
                            result[i + 1] = rest_results[i];
                        }
                        free(rest_results);
                        *Us_out = result;
                        *r_out  = 1 + rest_count;
                        return true;
                    }
                }
            }

            /* Advance to next k-subset. */
            uint64_t nxt = gosper_next(mask);
            if (nxt <= mask) break;  /* overflow guard */
            mask = nxt;
        }
    }

    /* No subset yielded a valid lift.  P is irreducible bivariately
     * under this set of u's; return it as the sole factor. */
    BPoly** result = (BPoly**)malloc(sizeof(BPoly*));
    if (!result) return false;
    result[0] = bpoly_copy(P);
    *Us_out = result;
    *r_out = 1;
    return true;
}

/* ====================================================================== */
/*  High-level orchestrator: mvfactor_try_bivariate_monic                 */
/* ====================================================================== */

/* Evaluation points to try, in order of preference.  0 first because
 * it commonly gives a clean image; small magnitudes preferred to keep
 * coefficient growth bounded. */
static const int64_t MV_ALPHA_ORDER[] = {0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
#define MV_ALPHA_COUNT (sizeof(MV_ALPHA_ORDER) / sizeof(MV_ALPHA_ORDER[0]))

bool mvfactor_try_bivariate_monic(const BPoly* P,
                                  mvfactor_univariate_cb factor_cb,
                                  void* cb_user_data,
                                  BPoly*** factors_out,
                                  int* r_out) {
    *factors_out = NULL;
    *r_out = 0;
    if (!P || P->deg_x < 1) return false;

    /* Monic-in-x check: bpoly_lc_x(P) must be the constant 1. */
    const ZUPoly* lcx = bpoly_lc_x(P);
    if (!lcx || lcx->deg != 0 ||
        mpz_cmp_ui(lcx->c[0], 1) != 0) {
        return false;
    }

    int orig_deg_x = P->deg_x;

    for (size_t pi = 0; pi < MV_ALPHA_COUNT; pi++) {
        int64_t alpha = MV_ALPHA_ORDER[pi];
    
        /* Build the univariate image P(x, alpha). */
        ZUPoly* image = bpoly_eval_y_si(P, alpha);
        if (!image || image->deg != orig_deg_x) {
            /* Image lost x-degree -- the chosen alpha killed the
             * leading coefficient.  Skip. */
            zupoly_free(image);
            continue;
        }
    
        /* Squarefree check: gcd(image, image') must be a unit.  The
         * lift's Diophantine step requires pairwise-coprime factors,
         * which holds when the image is squarefree. */
        ZUPoly* image_deriv = zupoly_new(image->deg);
        for (int i = 1; i <= image->deg; i++) {
            const mpz_t* ci = zupoly_getcoef(image, i);
            if (!ci) continue;
            mpz_t scaled; mpz_init(scaled);
            mpz_mul_ui(scaled, *ci, (unsigned long)i);
            zupoly_setcoef(image_deriv, i - 1, scaled);
            mpz_clear(scaled);
        }
            ZUPoly* g = zupoly_gcd(image, image_deriv);
            bool sqf = (g->deg == 0);
        zupoly_free(g); zupoly_free(image_deriv);
        if (!sqf) {
            zupoly_free(image);
            continue;
        }

        /* Factor the image via the caller-supplied callback. */
            ZUPoly** us = NULL;
        int r = 0;
        bool fac_ok = factor_cb(image, &us, &r, cb_user_data);
            zupoly_free(image);
        if (!fac_ok) continue;

        /* Image is irreducible univariately => P is (almost
         * certainly) irreducible bivariately too.  Return P as the
         * sole factor. */
        if (r == 1) {
            for (int j = 0; j < r; j++) zupoly_free(us[j]);
            free(us);
            BPoly** result = (BPoly**)malloc(sizeof(BPoly*));
            result[0] = bpoly_copy(P);
            *factors_out = result;
            *r_out = 1;
            return true;
        }

        /* Shift y -> y + alpha so the lift sees an "image at y = 0"
         * that matches our factoring of P(x, alpha). */
        BPoly* P_shifted = bpoly_shift_y_si(P, alpha);

        BPoly** lifted = NULL;
        int lifted_count = 0;
        bool lift_ok = bpoly_hensel_lift_multi(P_shifted,
                                               (const ZUPoly**)us, r,
                                               &lifted, &lifted_count);
        bpoly_free(P_shifted);
        for (int j = 0; j < r; j++) zupoly_free(us[j]);
        free(us);
        if (!lift_ok) continue;

        /* Undo the shift on each factor: y -> y - alpha.  The lift
         * may have returned fewer factors than `r` when subset
         * recombination merged some univariate factors into a single
         * bivariate one. */
        BPoly** result = (BPoly**)malloc(sizeof(BPoly*) * (size_t)lifted_count);
        for (int j = 0; j < lifted_count; j++) {
            result[j] = bpoly_shift_y_si(lifted[j], -alpha);
            bpoly_free(lifted[j]);
        }
        free(lifted);

        *factors_out = result;
        *r_out = lifted_count;
        return true;
    }

    /* No suitable evaluation point yielded a successful lift. */
    return false;
}
