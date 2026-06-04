/* groebnerwalk.c
 *
 * The Collart–Kalkbrener–Mall Gröbner walk: convert a Gröbner basis from
 * a cheap source order (DegreeReverseLexicographic) to an expensive
 * target order (Lexicographic or a user weight matrix) by walking the
 * Gröbner fan along the weight-vector path
 *
 *     w(t) = (1 - t) * tau + t * sigma,   t in [0, 1]
 *
 * where tau = (1, ..., 1) lies in the (closure of the) grevlex cone and
 * sigma is the first weight row of the target order.  At each cone wall
 * crossed along the path the basis is rewritten through the initial-form
 * ideal: compute in_w(G), take its Gröbner basis H under the target
 * order, then "lift" each h in H back to a full ideal element by dividing
 * it (with quotients) by the initial forms under the *current* order and
 * substituting the full polynomials.
 *
 * Correctness rests on two facts.  (1) If G is a Gröbner basis under an
 * order that refines the weight w, then { in_w(g) } is a Gröbner basis of
 * the initial ideal in_w(I) under that same order, so the lift division
 * always has remainder zero.  (2) The reduced Gröbner basis under a fixed
 * admissible order is unique, so the walk's output is byte-identical to a
 * direct Buchberger run under the target order.  Fact (2) also means any
 * internal safety-guard trip can fall back to gb_buchberger with no loss
 * of correctness -- exactly what this module does on overflow, a failed
 * lift (non-zero remainder from a degenerate path), or a step-count
 * blow-up.
 *
 * All arithmetic is exact (GMP mpq for the path parameter, mpz-guarded
 * int64 for weights); ownership follows the rest of the subsystem
 * (every GBPoly is heap-allocated and released with gb_poly_free).
 */

#include "groebner.h"

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "core.h"        /* tc_check_deadline() */

/* Hard cap on wall crossings; a well-behaved walk crosses far fewer.
 * Hitting it means the path is degenerate -- fall back to Buchberger. */
#define GB_WALK_MAX_STEPS 100000

/* ------------------------------------------------------------------ */
/*  Weight-matrix construction                                         */
/* ------------------------------------------------------------------ */

/* Build the refined order matrix [ w ; target-rows ]: the leading row is
 * the current weight vector, the remaining rows are the target order's
 * matrix (identity for Lexicographic).  Always a valid term order because
 * the target rows alone have full rank.  Caller frees with free_wmat. */
static GBWeightMatrix* build_refine_matrix(const int64_t* w, int n_vars,
                                           GBOrder target_order,
                                           const GBWeightMatrix* target_wmat) {
    int trows = (target_order == GB_ORDER_MATRIX) ? target_wmat->n_rows
                                                  : n_vars; /* lex -> identity */
    int nr = 1 + trows;
    GBWeightMatrix* m = (GBWeightMatrix*)malloc(sizeof(GBWeightMatrix));
    m->n_rows = nr;
    m->n_vars = n_vars;
    m->w = (int64_t*)malloc(sizeof(int64_t) * (size_t)nr * (size_t)n_vars);

    for (int v = 0; v < n_vars; v++) m->w[v] = w[v];   /* row 0 = weight */

    if (target_order == GB_ORDER_MATRIX) {
        for (int r = 0; r < trows; r++)
            for (int v = 0; v < n_vars; v++)
                m->w[(size_t)(1 + r) * n_vars + v]
                    = target_wmat->w[(size_t)r * n_vars + v];
    } else {                                            /* Lexicographic */
        for (int r = 0; r < n_vars; r++)
            for (int v = 0; v < n_vars; v++)
                m->w[(size_t)(1 + r) * n_vars + v] = (r == v) ? 1 : 0;
    }
    return m;
}

static void free_wmat(GBWeightMatrix* m) {
    if (m) { free(m->w); free(m); }
}

/* Set every basis element to the pure target order (re-sorting each).  No
 * weight matrix is allocated -- a Lexicographic basis carries no wmat, a
 * matrix-target basis borrows the caller-owned `target_wmat`. */
static void set_basis_target(GBPoly** G, size_t nG,
                             GBOrder target_order,
                             const GBWeightMatrix* target_wmat) {
    for (size_t i = 0; i < nG; i++) {
        if (target_order == GB_ORDER_MATRIX) {
            gb_poly_set_wmat(G[i], target_wmat);
        } else {
            G[i]->order = target_order;
            G[i]->wmat = NULL;
        }
        gb_poly_normalize(G[i]);
    }
}

/* ------------------------------------------------------------------ */
/*  Path geometry                                                      */
/* ------------------------------------------------------------------ */

/* Integer weight vector at path parameter t = p/q:
 *     w_v = (q - p) * tau_v + p * sigma_v   with tau_v == 1,
 * i.e. (q - p) + p * sigma_v   (scaling a weight row by the positive q
 * does not change the order).  Returns false on int64 overflow, which
 * triggers the Buchberger fallback. */
static bool weight_at_int(const int64_t* sigma, const mpq_t t,
                          int n_vars, int64_t* w_out) {
    mpz_t p, q, qmp, term, acc;
    mpz_init(p); mpz_init(q); mpz_init(qmp); mpz_init(term); mpz_init(acc);
    mpq_get_num(p, t);
    mpq_get_den(q, t);
    mpz_sub(qmp, q, p);
    bool ok = true;
    for (int v = 0; v < n_vars; v++) {
        mpz_set(acc, qmp);                 /* (q - p) * tau_v, tau_v = 1 */
        mpz_set_si(term, (long)sigma[v]);
        mpz_mul(term, term, p);            /* p * sigma_v */
        mpz_add(acc, acc, term);
        if (!mpz_fits_slong_p(acc)) { ok = false; break; }
        w_out[v] = (int64_t)mpz_get_si(acc);
    }
    mpz_clear(p); mpz_clear(q); mpz_clear(qmp); mpz_clear(term); mpz_clear(acc);
    return ok;
}

/* Smallest path parameter strictly greater than `t_cur` (and <= 1) at
 * which the leading term of some basis polynomial changes -- i.e. a
 * non-leading term b ties the current leading term a in w-weight and is
 * about to overtake it.  For a = LT(g) and another term b, with d = a - b:
 *     w(t).d = (1 - t) * (tau.d) + t * (sigma.d) = (1 - t) A + t B,
 * which vanishes at t* = A / (A - B).  The crossing is a genuine wall iff
 * the line is decreasing there (slope B - A < 0, i.e. A - B > 0) and
 * t_cur < t* <= 1.  Returns 1 and sets *t_next on success, 0 if no wall
 * remains, -1 on int64 overflow (fall back). */
static int next_wall(GBPoly* const* G, size_t nG, const int64_t* sigma,
                     int n_vars, const mpq_t t_cur, mpq_t t_next) {
    bool found = false;
    mpq_t cand;
    mpq_init(cand);
    for (size_t gi = 0; gi < nG; gi++) {
        const GBPoly* g = G[gi];
        if (g->n_terms < 2) continue;
        const int* a = g->exps;                 /* leading term (term 0) */
        for (size_t k = 1; k < g->n_terms; k++) {
            const int* b = g->exps + k * (size_t)n_vars;
            int64_t A = 0, B = 0;               /* A = tau.d, B = sigma.d */
            bool overflow = false;
            for (int v = 0; v < n_vars; v++) {
                int64_t d = (int64_t)a[v] - (int64_t)b[v];
                int64_t prod;
                if (__builtin_add_overflow(A, d, &A)
                 || __builtin_mul_overflow(sigma[v], d, &prod)
                 || __builtin_add_overflow(B, prod, &B)) {
                    overflow = true;
                    break;
                }
            }
            if (overflow) { mpq_clear(cand); return -1; }
            int64_t den = A - B;                 /* A - B */
            if (den <= 0) continue;              /* not a forward wall */
            /* t* = A / den, den > 0. */
            mpq_set_si(cand, (long)A, (unsigned long)den);
            mpq_canonicalize(cand);
            if (mpq_cmp(cand, t_cur) > 0 && mpq_cmp_ui(cand, 1, 1) <= 0) {
                if (!found || mpq_cmp(cand, t_next) < 0) {
                    mpq_set(t_next, cand);
                    found = true;
                }
            }
        }
    }
    mpq_clear(cand);
    return found ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/*  Conversion at a weight                                             */
/* ------------------------------------------------------------------ */

/* One walk step at weight `w`.  Replaces *G_io in place with a Gröbner
 * basis under `new_order`/`new_wmat`, computed from the current basis:
 *   1. Iw = { in_w(g) } under the current (previous) order;
 *   2. H  = Gröbner basis of <Iw> under the target order;
 *   3. lift each h by dividing it (with quotients) by Iw under the current
 *      order, then substituting the full g_j for each in_w(g_j);
 *   4. re-sort the lifts to `new_order`.
 * Returns false (caller falls back to Buchberger) if any lift leaves a
 * non-zero remainder -- the signature of a degenerate path where the
 * current order no longer refines w. */
static bool convert_at(GBPoly*** G_io, size_t* nG_io,
                       const int64_t* w, int n_vars,
                       GBOrder new_order, const GBWeightMatrix* new_wmat,
                       GBOrder target_order, const GBWeightMatrix* target_wmat) {
    GBPoly** G = *G_io;
    size_t nG = *nG_io;
    if (nG == 0) return true;

    GBOrder prev_order = G[0]->order;
    int prev_elim = G[0]->elim_pivot;
    const GBWeightMatrix* prev_wmat = G[0]->wmat;

    /* Initial forms: IwPrev keeps the current order (its leading terms
     * form a Gröbner basis of in_w(I) -- used for the lift division);
     * IwTgt is the same set re-sorted to the target order for computing
     * H. */
    GBPoly** IwPrev = (GBPoly**)malloc(sizeof(GBPoly*) * nG);
    GBPoly** IwTgt  = (GBPoly**)malloc(sizeof(GBPoly*) * nG);
    for (size_t j = 0; j < nG; j++) {
        IwPrev[j] = gb_initial_form(G[j], w, n_vars);
        IwTgt[j]  = gb_poly_copy(IwPrev[j]);
        if (target_order == GB_ORDER_MATRIX) {
            gb_poly_set_wmat(IwTgt[j], target_wmat);
        } else {
            IwTgt[j]->order = target_order;
            IwTgt[j]->wmat = NULL;
        }
        gb_poly_normalize(IwTgt[j]);
    }

    size_t nH = 0;
    GBPoly** H = gb_buchberger(IwTgt, nG, &nH);

    GBPoly** G2 = (nH ? (GBPoly**)malloc(sizeof(GBPoly*) * nH) : NULL);
    size_t nG2 = 0;
    bool ok = true;

    for (size_t hi = 0; hi < nH && ok; hi++) {
        GBPoly** quot = NULL;
        GBPoly* rem = gb_divmod(H[hi], IwPrev, nG, &quot);
        if (rem->n_terms != 0) ok = false;       /* degenerate path */
        gb_poly_free(rem);

        if (ok) {
            /* lift = sum_j sum_k quot[j].term_k * g_j  (full polynomials). */
            GBPoly* acc = gb_poly_new(n_vars, prev_order, prev_elim);
            acc->wmat = prev_wmat;
            for (size_t j = 0; j < nG; j++) {
                const GBPoly* qj = quot[j];
                for (size_t k = 0; k < qj->n_terms; k++) {
                    const int* e = qj->exps + k * (size_t)n_vars;
                    GBPoly* tmp = gb_poly_mul_by_monomial(G[j], e, qj->coefs[k]);
                    GBPoly* a2 = gb_poly_add(acc, tmp);
                    gb_poly_free(tmp);
                    gb_poly_free(acc);
                    acc = a2;
                }
            }
            if (new_order == GB_ORDER_MATRIX) {
                gb_poly_set_wmat(acc, new_wmat);
            } else {
                acc->order = new_order;
                acc->wmat = NULL;
            }
            gb_poly_normalize(acc);
            if (acc->n_terms > 0) G2[nG2++] = acc;
            else gb_poly_free(acc);
        }

        if (quot) {
            for (size_t j = 0; j < nG; j++) gb_poly_free(quot[j]);
            free(quot);
        }
    }

    for (size_t j = 0; j < nG; j++) { gb_poly_free(IwPrev[j]); gb_poly_free(IwTgt[j]); }
    free(IwPrev);
    free(IwTgt);
    gb_basis_free(H, nH);

    if (!ok) {
        for (size_t i = 0; i < nG2; i++) gb_poly_free(G2[i]);
        free(G2);
        return false;
    }

    gb_basis_free(G, nG);
    *G_io = G2;
    *nG_io = nG2;
    return true;
}

/* ------------------------------------------------------------------ */
/*  Driver                                                             */
/* ------------------------------------------------------------------ */

GBPoly** gb_groebner_walk(GBPoly* const* F, size_t n,
                          GBOrder target_order,
                          const GBWeightMatrix* target_wmat,
                          size_t* out_n) {
    if (n == 0) { *out_n = 0; return NULL; }
    int n_vars = F[0]->n_vars;

    /* Trivial cases: source cone == target cone, or a single variable
     * (all orders coincide).  Compute directly in the target order. */
    if (target_order == GB_ORDER_GREVLEX || n_vars <= 1) {
        return gb_buchberger(F, n, out_n);
    }

    /* ---- Source basis: DegreeReverseLexicographic (cheap). ---------- */
    GBPoly** Fg = (GBPoly**)malloc(sizeof(GBPoly*) * n);
    for (size_t i = 0; i < n; i++) {
        Fg[i] = gb_poly_copy(F[i]);
        Fg[i]->order = GB_ORDER_GREVLEX;
        Fg[i]->wmat = NULL;
        gb_poly_normalize(Fg[i]);
    }
    size_t nG = 0;
    GBPoly** G = gb_buchberger(Fg, n, &nG);
    for (size_t i = 0; i < n; i++) gb_poly_free(Fg[i]);
    free(Fg);

    if (nG == 0) { *out_n = 0; return G; }

    /* A constant in the source basis means the ideal is <1>; the canonical
     * answer is the same in every order, so defer to a direct run. */
    for (size_t i = 0; i < nG; i++) {
        if (gb_poly_is_constant(G[i])) {
            gb_basis_free(G, nG);
            return gb_buchberger(F, n, out_n);
        }
    }

    /* ---- Path endpoints. -------------------------------------------- */
    int64_t* tau = (int64_t*)malloc(sizeof(int64_t) * (size_t)n_vars);
    int64_t* sigma = (int64_t*)malloc(sizeof(int64_t) * (size_t)n_vars);
    int64_t* wint = (int64_t*)malloc(sizeof(int64_t) * (size_t)n_vars);
    for (int v = 0; v < n_vars; v++) tau[v] = 1;
    if (target_order == GB_ORDER_MATRIX) {
        for (int v = 0; v < n_vars; v++) sigma[v] = target_wmat->w[v]; /* row 0 */
    } else {                                          /* Lexicographic */
        for (int v = 0; v < n_vars; v++) sigma[v] = (v == 0) ? 1 : 0;
    }

    GBWeightMatrix* prev_mat = NULL;       /* matrix the basis currently borrows */
    mpq_t t, t_next;
    mpq_init(t); mpq_init(t_next);
    bool fail = false;

    /* ---- Initial conversion at tau. --------------------------------- */
    /* grevlex refines tau, so in_tau(G) generates in_tau(I); this lands
     * the basis in the order (tau | target). */
    {
        GBWeightMatrix* m = build_refine_matrix(tau, n_vars,
                                                target_order, target_wmat);
        if (!convert_at(&G, &nG, tau, n_vars, GB_ORDER_MATRIX, m,
                        target_order, target_wmat)) {
            free_wmat(m);
            fail = true;
        } else {
            prev_mat = m;
        }
    }

    /* ---- Walk the remaining path. ----------------------------------- */
    int steps = 0;
    while (!fail) {
        if (++steps > GB_WALK_MAX_STEPS) { fail = true; break; }
        tc_check_deadline();

        int wr = next_wall(G, nG, sigma, n_vars, t, t_next);
        if (wr < 0) { fail = true; break; }      /* overflow */
        if (wr == 0) {
            /* No wall ahead: the current cone covers the rest of the path,
             * so the basis is already a Gröbner basis of the target. */
            set_basis_target(G, nG, target_order, target_wmat);
            gb_finalize_basis(G, &nG);
            break;
        }

        bool is_final = (mpq_cmp_ui(t_next, 1, 1) == 0);
        if (!weight_at_int(sigma, t_next, n_vars, wint)) { fail = true; break; }

        GBWeightMatrix* m = NULL;
        GBOrder no = target_order;
        const GBWeightMatrix* nw = target_wmat;
        if (!is_final) {
            m = build_refine_matrix(wint, n_vars, target_order, target_wmat);
            no = GB_ORDER_MATRIX;
            nw = m;
        }
        if (!convert_at(&G, &nG, wint, n_vars, no, nw,
                        target_order, target_wmat)) {
            free_wmat(m);
            fail = true;
            break;
        }
        free_wmat(prev_mat);
        prev_mat = m;
        mpq_set(t, t_next);

        if (is_final) {
            gb_finalize_basis(G, &nG);
            break;
        }
    }

    mpq_clear(t); mpq_clear(t_next);
    free(tau); free(sigma); free(wint);
    free_wmat(prev_mat);

    if (fail) {
        /* Result-identical fallback: a direct Buchberger run in the target
         * order (F already carries that order from the caller). */
        gb_basis_free(G, nG);
        return gb_buchberger(F, n, out_n);
    }

    *out_n = nG;
    return G;
}
