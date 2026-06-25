/*
 * mpoly.c
 * -------
 * Sparse multivariate polynomial over Z with mpz_t coefficients.
 * See mpoly.h for the data layout, invariants, and ownership rules.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "mpoly.h"
#include "expr.h"
#include "sym_names.h"

/* ---------------------------------------------------------------- */
/*  Internal helpers                                                */
/* ---------------------------------------------------------------- */

static inline int* mpoly_exp_at(const MPoly* p, size_t i) {
    return p->exps + i * (size_t)p->n_vars;
}

/* Lex descending compare on exponent tuples.  Returns:
 *   < 0 if a should sort BEFORE b (a is the larger monomial),
 *   > 0 if a should sort AFTER b,
 *     0 if equal. */
static int monomial_cmp_lex(const int* a, const int* b, int n_vars) {
    for (int i = 0; i < n_vars; i++) {
        if (a[i] != b[i]) return b[i] - a[i];
    }
    return 0;
}

/* ---------------------------------------------------------------- */
/*  Construction / destruction                                      */
/* ---------------------------------------------------------------- */

MPoly* mpoly_new(int n_vars) {
    MPoly* p = (MPoly*)malloc(sizeof(MPoly));
    p->n_vars  = n_vars;
    p->exps    = NULL;
    p->coefs   = NULL;
    p->n_terms = 0;
    p->cap     = 0;
    return p;
}

MPoly* mpoly_zero(int n_vars) {
    return mpoly_new(n_vars);
}

void mpoly_reserve(MPoly* p, size_t cap) {
    if (cap <= p->cap) return;
    size_t new_cap = p->cap == 0 ? 8 : p->cap;
    while (new_cap < cap) new_cap *= 2;

    if (p->n_vars > 0) {
        p->exps = (int*)realloc(p->exps,
                                sizeof(int) * new_cap * (size_t)p->n_vars);
    }
    mpz_t* new_coefs = (mpz_t*)realloc(p->coefs, sizeof(mpz_t) * new_cap);
    /* Initialise newly-grown slots so push_term can mpz_set into them. */
    for (size_t i = p->cap; i < new_cap; i++) mpz_init(new_coefs[i]);
    p->coefs = new_coefs;
    p->cap = new_cap;
}

MPoly* mpoly_from_int(int n_vars, int64_t c) {
    MPoly* p = mpoly_new(n_vars);
    if (c == 0) return p;
    mpoly_reserve(p, 1);
    if (n_vars > 0) {
        memset(mpoly_exp_at(p, 0), 0, sizeof(int) * (size_t)n_vars);
    }
    mpz_set_si(p->coefs[0], (long)c);
    p->n_terms = 1;
    return p;
}

MPoly* mpoly_from_mpz(int n_vars, const mpz_t c) {
    MPoly* p = mpoly_new(n_vars);
    if (mpz_sgn(c) == 0) return p;
    mpoly_reserve(p, 1);
    if (n_vars > 0) {
        memset(mpoly_exp_at(p, 0), 0, sizeof(int) * (size_t)n_vars);
    }
    mpz_set(p->coefs[0], c);
    p->n_terms = 1;
    return p;
}

MPoly* mpoly_monomial(int n_vars, int var_idx, int k, int64_t c) {
    MPoly* p = mpoly_new(n_vars);
    if (c == 0 || var_idx < 0 || var_idx >= n_vars) return p;
    mpoly_reserve(p, 1);
    int* row = mpoly_exp_at(p, 0);
    memset(row, 0, sizeof(int) * (size_t)n_vars);
    row[var_idx] = k;
    mpz_set_si(p->coefs[0], (long)c);
    p->n_terms = 1;
    return p;
}

MPoly* mpoly_copy(const MPoly* p) {
    MPoly* q = mpoly_new(p->n_vars);
    if (p->n_terms == 0) return q;
    mpoly_reserve(q, p->n_terms);
    if (p->n_vars > 0) {
        memcpy(q->exps, p->exps,
               sizeof(int) * p->n_terms * (size_t)p->n_vars);
    }
    for (size_t i = 0; i < p->n_terms; i++) mpz_set(q->coefs[i], p->coefs[i]);
    q->n_terms = p->n_terms;
    return q;
}

void mpoly_free(MPoly* p) {
    if (!p) return;
    if (p->coefs) {
        for (size_t i = 0; i < p->cap; i++) mpz_clear(p->coefs[i]);
        free(p->coefs);
    }
    free(p->exps);
    free(p);
}

/* ---------------------------------------------------------------- */
/*  Term manipulation                                               */
/* ---------------------------------------------------------------- */

void mpoly_push_term(MPoly* p, const int* exps, const mpz_t coef) {
    if (mpz_sgn(coef) == 0) return;
    if (p->n_terms == p->cap) mpoly_reserve(p, p->cap == 0 ? 8 : p->cap * 2);
    if (p->n_vars > 0) {
        memcpy(mpoly_exp_at(p, p->n_terms), exps,
               sizeof(int) * (size_t)p->n_vars);
    }
    mpz_set(p->coefs[p->n_terms], coef);
    p->n_terms++;
}

void mpoly_push_term_si(MPoly* p, const int* exps, int64_t coef) {
    if (coef == 0) return;
    if (p->n_terms == p->cap) mpoly_reserve(p, p->cap == 0 ? 8 : p->cap * 2);
    if (p->n_vars > 0) {
        memcpy(mpoly_exp_at(p, p->n_terms), exps,
               sizeof(int) * (size_t)p->n_vars);
    }
    mpz_set_si(p->coefs[p->n_terms], (long)coef);
    p->n_terms++;
}

/* Sort context for normalize (single-threaded). */
static const MPoly* g_sort_ctx = NULL;

static int sort_cmp_idx(const void* a, const void* b) {
    size_t ia = *(const size_t*)a;
    size_t ib = *(const size_t*)b;
    int c = monomial_cmp_lex(
        mpoly_exp_at(g_sort_ctx, ia),
        mpoly_exp_at(g_sort_ctx, ib),
        g_sort_ctx->n_vars);
    if (c != 0) return c;
    /* Stable order on equal monomials: keep insertion order so the
     * merge step combines them deterministically. */
    return (ia < ib) ? -1 : (ia > ib);
}

void mpoly_normalize(MPoly* p) {
    if (p->n_terms < 2) {
        /* Single term: drop if zero. */
        if (p->n_terms == 1 && mpz_sgn(p->coefs[0]) == 0) p->n_terms = 0;
        return;
    }

    /* Build index permutation. */
    size_t* idx = (size_t*)malloc(sizeof(size_t) * p->n_terms);
    for (size_t i = 0; i < p->n_terms; i++) idx[i] = i;

    g_sort_ctx = p;
    qsort(idx, p->n_terms, sizeof(size_t), sort_cmp_idx);
    g_sort_ctx = NULL;

    /* Walk the sorted permutation, materialising into fresh arrays
     * with duplicate-monomial coefficients summed and zero-coef terms
     * dropped. */
    int*   new_exps = (p->n_vars > 0)
                      ? (int*)malloc(sizeof(int) * p->n_terms * (size_t)p->n_vars)
                      : NULL;
    mpz_t* new_coefs = (mpz_t*)malloc(sizeof(mpz_t) * p->n_terms);
    for (size_t i = 0; i < p->n_terms; i++) mpz_init(new_coefs[i]);

    size_t out = 0;
    for (size_t i = 0; i < p->n_terms; i++) {
        size_t src = idx[i];
        if (out > 0 && p->n_vars > 0 &&
            monomial_cmp_lex(new_exps + (out - 1) * (size_t)p->n_vars,
                             mpoly_exp_at(p, src),
                             p->n_vars) == 0) {
            mpz_add(new_coefs[out - 1], new_coefs[out - 1], p->coefs[src]);
            if (mpz_sgn(new_coefs[out - 1]) == 0) {
                /* Cancellation: drop the slot. */
                mpz_set_ui(new_coefs[out - 1], 0);
                out--;
            }
        } else if (out > 0 && p->n_vars == 0) {
            /* All terms in n_vars==0 land share the empty monomial. */
            mpz_add(new_coefs[out - 1], new_coefs[out - 1], p->coefs[src]);
            if (mpz_sgn(new_coefs[out - 1]) == 0) {
                mpz_set_ui(new_coefs[out - 1], 0);
                out--;
            }
        } else {
            if (mpz_sgn(p->coefs[src]) == 0) continue;
            if (p->n_vars > 0) {
                memcpy(new_exps + out * (size_t)p->n_vars,
                       mpoly_exp_at(p, src),
                       sizeof(int) * (size_t)p->n_vars);
            }
            mpz_set(new_coefs[out], p->coefs[src]);
            out++;
        }
    }

    /* Replace storage.  new_coefs was sized to the *original* term count
     * (old_n) with every slot mpz_init'd; only the first `out` are live
     * after merging.  Capture old_n before overwriting p->n_terms so the
     * trailing initialised-but-unused slots can be cleared (otherwise they
     * leak: cap shrinks to `out`, so mpoly_free never reaches them). */
    size_t old_n = p->n_terms;
    for (size_t i = 0; i < p->cap; i++) mpz_clear(p->coefs[i]);
    free(p->coefs);
    free(p->exps);

    p->exps  = new_exps;
    p->coefs = new_coefs;

    /* Clear the unused tail [out, old_n), then trim arrays to `out`. */
    for (size_t i = out; i < old_n; i++) mpz_clear(new_coefs[i]);
    if (out == 0) {
        free(p->coefs);
        free(p->exps);
        p->coefs = NULL;
        p->exps  = NULL;
    } else if (out < old_n) {
        p->coefs = (mpz_t*)realloc(p->coefs, sizeof(mpz_t) * out);
        if (p->n_vars > 0) {
            p->exps = (int*)realloc(p->exps,
                                    sizeof(int) * out * (size_t)p->n_vars);
        }
    }
    p->n_terms = out;
    p->cap     = out;

    free(idx);
}

/* Binary search for the position of `target` in p (sorted lex desc).
 * Returns the index in [0, n_terms] where the target should be
 * inserted (or where it currently lives if found).  Sets *found. */
static size_t mpoly_search(const MPoly* p, const int* target, bool* found) {
    *found = false;
    size_t lo = 0, hi = p->n_terms;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int c = monomial_cmp_lex(mpoly_exp_at(p, mid), target, p->n_vars);
        if (c == 0) { *found = true; return mid; }
        if (c < 0) lo = mid + 1;
        else       hi = mid;
    }
    return lo;
}

void mpoly_set_coef(MPoly* p, const int* exps, const mpz_t coef) {
    bool found;
    size_t pos = mpoly_search(p, exps, &found);
    if (found) {
        if (mpz_sgn(coef) == 0) {
            /* Remove. */
            mpz_clear(p->coefs[pos]);
            mpz_init(p->coefs[pos]);  /* keep slot live */
            /* Shift down. */
            if (p->n_vars > 0) {
                memmove(mpoly_exp_at(p, pos),
                        mpoly_exp_at(p, pos + 1),
                        sizeof(int) * (p->n_terms - pos - 1) *
                            (size_t)p->n_vars);
            }
            for (size_t i = pos; i + 1 < p->n_terms; i++) {
                mpz_swap(p->coefs[i], p->coefs[i + 1]);
            }
            p->n_terms--;
        } else {
            mpz_set(p->coefs[pos], coef);
        }
        return;
    }
    if (mpz_sgn(coef) == 0) return;

    /* Insert at pos. */
    if (p->n_terms == p->cap) {
        mpoly_reserve(p, p->cap == 0 ? 8 : p->cap * 2);
    }
    /* Shift up. */
    if (p->n_vars > 0) {
        memmove(mpoly_exp_at(p, pos + 1),
                mpoly_exp_at(p, pos),
                sizeof(int) * (p->n_terms - pos) * (size_t)p->n_vars);
        memcpy(mpoly_exp_at(p, pos), exps, sizeof(int) * (size_t)p->n_vars);
    }
    /* Shift coef slots up via swap to preserve mpz_init state. */
    for (size_t i = p->n_terms; i > pos; i--) {
        mpz_swap(p->coefs[i], p->coefs[i - 1]);
    }
    mpz_set(p->coefs[pos], coef);
    p->n_terms++;
}

const mpz_t* mpoly_get_coef(const MPoly* p, const int* exps) {
    bool found;
    size_t pos = mpoly_search(p, exps, &found);
    if (!found) return NULL;
    return (const mpz_t*)&p->coefs[pos];
}

/* ---------------------------------------------------------------- */
/*  Predicates and queries                                          */
/* ---------------------------------------------------------------- */

bool mpoly_is_zero(const MPoly* p) {
    return p->n_terms == 0;
}

bool mpoly_eq(const MPoly* a, const MPoly* b) {
    if (a->n_vars != b->n_vars) return false;
    if (a->n_terms != b->n_terms) return false;
    for (size_t i = 0; i < a->n_terms; i++) {
        if (a->n_vars > 0) {
            if (memcmp(mpoly_exp_at(a, i), mpoly_exp_at(b, i),
                       sizeof(int) * (size_t)a->n_vars) != 0)
                return false;
        }
        if (mpz_cmp(a->coefs[i], b->coefs[i]) != 0) return false;
    }
    return true;
}

int mpoly_deg_var(const MPoly* p, int var_idx) {
    if (var_idx < 0 || var_idx >= p->n_vars) return -1;
    if (p->n_terms == 0) return -1;
    int d = 0;
    for (size_t i = 0; i < p->n_terms; i++) {
        int e = mpoly_exp_at((MPoly*)p, i)[var_idx];
        if (e > d) d = e;
    }
    return d;
}

int mpoly_total_deg(const MPoly* p) {
    if (p->n_terms == 0) return -1;
    int td = 0;
    for (size_t i = 0; i < p->n_terms; i++) {
        int sum = 0;
        const int* row = mpoly_exp_at((MPoly*)p, i);
        for (int v = 0; v < p->n_vars; v++) sum += row[v];
        if (sum > td) td = sum;
    }
    return td;
}

bool mpoly_is_constant_in_var(const MPoly* p, int var_idx) {
    if (var_idx < 0 || var_idx >= p->n_vars) return true;
    if (p->n_terms == 0) return true;
    for (size_t i = 0; i < p->n_terms; i++) {
        if (mpoly_exp_at((MPoly*)p, i)[var_idx] != 0) return false;
    }
    return true;
}

/* ---------------------------------------------------------------- */
/*  Arithmetic                                                      */
/* ---------------------------------------------------------------- */

MPoly* mpoly_add(const MPoly* a, const MPoly* b) {
    /* a, b are normalised (sorted desc).  Merge. */
    MPoly* c = mpoly_new(a->n_vars);
    mpoly_reserve(c, a->n_terms + b->n_terms);

    size_t i = 0, j = 0;
    while (i < a->n_terms && j < b->n_terms) {
        int cmp = (a->n_vars > 0)
                  ? monomial_cmp_lex(mpoly_exp_at((MPoly*)a, i),
                                     mpoly_exp_at((MPoly*)b, j),
                                     a->n_vars)
                  : 0;
        if (cmp < 0) {
            mpoly_push_term(c, mpoly_exp_at((MPoly*)a, i), a->coefs[i]);
            i++;
        } else if (cmp > 0) {
            mpoly_push_term(c, mpoly_exp_at((MPoly*)b, j), b->coefs[j]);
            j++;
        } else {
            mpz_t sum; mpz_init(sum);
            mpz_add(sum, a->coefs[i], b->coefs[j]);
            if (mpz_sgn(sum) != 0) {
                mpoly_push_term(c, mpoly_exp_at((MPoly*)a, i), sum);
            }
            mpz_clear(sum);
            i++; j++;
        }
    }
    while (i < a->n_terms) {
        mpoly_push_term(c, mpoly_exp_at((MPoly*)a, i), a->coefs[i]);
        i++;
    }
    while (j < b->n_terms) {
        mpoly_push_term(c, mpoly_exp_at((MPoly*)b, j), b->coefs[j]);
        j++;
    }
    return c;
}

MPoly* mpoly_sub(const MPoly* a, const MPoly* b) {
    MPoly* nb = mpoly_neg(b);
    MPoly* r = mpoly_add(a, nb);
    mpoly_free(nb);
    return r;
}

MPoly* mpoly_neg(const MPoly* a) {
    MPoly* c = mpoly_copy(a);
    for (size_t i = 0; i < c->n_terms; i++) mpz_neg(c->coefs[i], c->coefs[i]);
    return c;
}

MPoly* mpoly_scale(const MPoly* a, const mpz_t k) {
    if (mpz_sgn(k) == 0) return mpoly_zero(a->n_vars);
    if (mpz_cmp_ui(k, 1) == 0) return mpoly_copy(a);
    MPoly* c = mpoly_copy(a);
    for (size_t i = 0; i < c->n_terms; i++) mpz_mul(c->coefs[i], c->coefs[i], k);
    return c;
}

MPoly* mpoly_scale_si(const MPoly* a, int64_t k) {
    if (k == 0) return mpoly_zero(a->n_vars);
    if (k == 1) return mpoly_copy(a);
    MPoly* c = mpoly_copy(a);
    for (size_t i = 0; i < c->n_terms; i++)
        mpz_mul_si(c->coefs[i], c->coefs[i], (long)k);
    return c;
}

MPoly* mpoly_mul(const MPoly* a, const MPoly* b) {
    if (mpoly_is_zero(a) || mpoly_is_zero(b)) return mpoly_zero(a->n_vars);

    MPoly* c = mpoly_new(a->n_vars);
    mpoly_reserve(c, a->n_terms * b->n_terms);

    /* VLA: stack-allocated row buffer for the sum-of-exponents. */
    int n_vars = a->n_vars;
    int row[n_vars > 0 ? n_vars : 1];
    mpz_t prod; mpz_init(prod);

    for (size_t i = 0; i < a->n_terms; i++) {
        for (size_t j = 0; j < b->n_terms; j++) {
            for (int v = 0; v < n_vars; v++) {
                row[v] = mpoly_exp_at((MPoly*)a, i)[v] +
                         mpoly_exp_at((MPoly*)b, j)[v];
            }
            mpz_mul(prod, a->coefs[i], b->coefs[j]);
            mpoly_push_term(c, row, prod);
        }
    }
    mpz_clear(prod);

    mpoly_normalize(c);
    return c;
}

/* ---------------------------------------------------------------- */
/*  Substitution and projection                                     */
/* ---------------------------------------------------------------- */

MPoly* mpoly_subst_var_mpz(const MPoly* p, int var_idx, const mpz_t alpha) {
    if (var_idx < 0 || var_idx >= p->n_vars) return mpoly_copy(p);

    MPoly* q = mpoly_new(p->n_vars);
    mpoly_reserve(q, p->n_terms);

    int n_vars = p->n_vars;
    int row[n_vars > 0 ? n_vars : 1];
    mpz_t scale; mpz_init(scale);
    mpz_t scaled; mpz_init(scaled);

    for (size_t i = 0; i < p->n_terms; i++) {
        const int* src = mpoly_exp_at((MPoly*)p, i);
        int e = src[var_idx];
        memcpy(row, src, sizeof(int) * (size_t)n_vars);
        row[var_idx] = 0;
        mpz_pow_ui(scale, alpha, (unsigned long)e);
        mpz_mul(scaled, p->coefs[i], scale);
        mpoly_push_term(q, row, scaled);
    }

    mpz_clear(scale); mpz_clear(scaled);
    mpoly_normalize(q);
    return q;
}

MPoly* mpoly_subst_var_int(const MPoly* p, int var_idx, int64_t alpha) {
    mpz_t a; mpz_init_set_si(a, (long)alpha);
    MPoly* r = mpoly_subst_var_mpz(p, var_idx, a);
    mpz_clear(a);
    return r;
}

/* Binomial coefficient (n choose k). */
static void mp_binom(mpz_t out, unsigned long n, unsigned long k) {
    mpz_bin_uiui(out, n, k);
}

MPoly* mpoly_shift_var_int(const MPoly* p, int var_idx, int64_t alpha) {
    if (var_idx < 0 || var_idx >= p->n_vars || alpha == 0)
        return mpoly_copy(p);

    MPoly* q = mpoly_new(p->n_vars);
    int n_vars = p->n_vars;
    int row[n_vars > 0 ? n_vars : 1];

    /* For each term c * (prod x_v^{e_v}) with x_{var_idx}^{e_v}, expand
     * (x_{var_idx} + alpha)^{e_v} = sum_{k=0}^{e_v} C(e_v, k) * alpha^{e_v-k}
     *                                  * x_{var_idx}^k. */
    mpz_t bin; mpz_init(bin);
    mpz_t pow_alpha; mpz_init(pow_alpha);
    mpz_t alpha_z; mpz_init_set_si(alpha_z, (long)alpha);
    mpz_t term; mpz_init(term);

    /* Pre-expansion size estimate: each term expands into (e_v + 1)
     * sub-terms; sum_i (e_i + 1) is an upper bound. */
    size_t expand_cap = 0;
    for (size_t i = 0; i < p->n_terms; i++) {
        int e = mpoly_exp_at((MPoly*)p, i)[var_idx];
        expand_cap += (size_t)(e + 1);
    }
    mpoly_reserve(q, expand_cap);

    for (size_t i = 0; i < p->n_terms; i++) {
        const int* src = mpoly_exp_at((MPoly*)p, i);
        int e = src[var_idx];
        memcpy(row, src, sizeof(int) * (size_t)n_vars);
        for (int k = 0; k <= e; k++) {
            row[var_idx] = k;
            mp_binom(bin, (unsigned long)e, (unsigned long)k);
            mpz_pow_ui(pow_alpha, alpha_z, (unsigned long)(e - k));
            mpz_mul(term, p->coefs[i], bin);
            mpz_mul(term, term, pow_alpha);
            mpoly_push_term(q, row, term);
        }
    }

    mpz_clear(bin); mpz_clear(pow_alpha); mpz_clear(alpha_z); mpz_clear(term);
    mpoly_normalize(q);
    return q;
}

MPoly* mpoly_coef_of_var(const MPoly* p, int var_idx, int k) {
    MPoly* q = mpoly_new(p->n_vars);
    if (var_idx < 0 || var_idx >= p->n_vars || k < 0) return q;

    int n_vars = p->n_vars;
    int row[n_vars > 0 ? n_vars : 1];

    for (size_t i = 0; i < p->n_terms; i++) {
        const int* src = mpoly_exp_at((MPoly*)p, i);
        if (src[var_idx] != k) continue;
        memcpy(row, src, sizeof(int) * (size_t)n_vars);
        row[var_idx] = 0;
        mpoly_push_term(q, row, p->coefs[i]);
    }
    /* Already in lex desc order since we walk p in order and skip
     * non-matching terms; but multiple terms may have the same
     * residual monomial only if the original p had duplicates --
     * not possible for a normalised p.  So q is normalised already. */
    return q;
}

MPoly* mpoly_lc_var(const MPoly* p, int var_idx) {
    int d = mpoly_deg_var(p, var_idx);
    if (d < 0) return mpoly_zero(p->n_vars);
    return mpoly_coef_of_var(p, var_idx, d);
}

/* ---------------------------------------------------------------- */
/*  Expr round-trip                                                 */
/* ---------------------------------------------------------------- */

/* Internal: try to interpret a single Expr term as a (coef, exps)
 * pair against the variable list.  Returns false on failure (non-
 * polynomial structure). */
static bool expr_term_to_mpoly(struct Expr* e, struct Expr** vars,
                               int n_vars, mpz_t coef_out, int* exps_out);

/* Attempt to match an Expr against a known variable. */
static int find_var_index(struct Expr* e, struct Expr** vars, int n_vars) {
    for (int i = 0; i < n_vars; i++) {
        if (expr_eq(e, vars[i])) return i;
    }
    return -1;
}

/* Helper: interpret an Expr as an integer; returns false if not a
 * literal integer.  Sets `out` (already initialised) on success. */
static bool expr_int_to_mpz(struct Expr* e, mpz_t out) {
    if (e->type == EXPR_INTEGER) {
        mpz_set_si(out, (long)e->data.integer);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpz_set(out, e->data.bigint);
        return true;
    }
    return false;
}

static bool expr_term_to_mpoly(struct Expr* e, struct Expr** vars,
                               int n_vars, mpz_t coef_out, int* exps_out) {
    /* Initialise outputs. */
    mpz_set_ui(coef_out, 1);
    for (int i = 0; i < n_vars; i++) exps_out[i] = 0;

    /* Base cases. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        return expr_int_to_mpz(e, coef_out);
    }
    if (e->type == EXPR_SYMBOL) {
        int vi = find_var_index(e, vars, n_vars);
        if (vi < 0) return false;
        exps_out[vi] = 1;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return false;

    /* Power[var, k]. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        struct Expr* base = e->data.function.args[0];
        struct Expr* exp  = e->data.function.args[1];
        int vi = find_var_index(base, vars, n_vars);
        if (vi < 0) return false;
        if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return false;
        exps_out[vi] = (int)exp->data.integer;
        return true;
    }

    /* Times[...] -- multiply the contributions of each arg. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times) {
        mpz_t sub_coef; mpz_init(sub_coef);
        int sub_exps[n_vars > 0 ? n_vars : 1];
        for (size_t a = 0; a < e->data.function.arg_count; a++) {
            if (!expr_term_to_mpoly(e->data.function.args[a], vars, n_vars,
                                    sub_coef, sub_exps)) {
                mpz_clear(sub_coef);
                return false;
            }
            mpz_mul(coef_out, coef_out, sub_coef);
            for (int v = 0; v < n_vars; v++) exps_out[v] += sub_exps[v];
        }
        mpz_clear(sub_coef);
        return true;
    }

    /* Some other compound -- not a polynomial term in our vars.
     * (But it might still be a variable Expr we should match before
     * giving up.) */
    int vi = find_var_index(e, vars, n_vars);
    if (vi >= 0) {
        exps_out[vi] = 1;
        return true;
    }
    return false;
}

MPoly* expr_to_mpoly(struct Expr* e, struct Expr** vars, int n_vars) {
    if (!e) return NULL;
    MPoly* p = mpoly_new(n_vars);

    /* Plus[...] -- each arg is a term. */
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Plus) {
        mpz_t coef; mpz_init(coef);
        int exps[n_vars > 0 ? n_vars : 1];
        for (size_t a = 0; a < e->data.function.arg_count; a++) {
            if (!expr_term_to_mpoly(e->data.function.args[a], vars, n_vars,
                                    coef, exps)) {
                mpz_clear(coef);
                mpoly_free(p);
                return NULL;
            }
            mpoly_push_term(p, exps, coef);
        }
        mpz_clear(coef);
        mpoly_normalize(p);
        return p;
    }

    /* Single term. */
    mpz_t coef; mpz_init(coef);
    int exps[n_vars > 0 ? n_vars : 1];
    if (!expr_term_to_mpoly(e, vars, n_vars, coef, exps)) {
        mpz_clear(coef);
        mpoly_free(p);
        return NULL;
    }
    mpoly_push_term(p, exps, coef);
    mpz_clear(coef);
    mpoly_normalize(p);
    return p;
}

struct Expr* mpoly_to_expr(const MPoly* p, struct Expr** vars) {
    if (p->n_terms == 0) return expr_new_integer(0);

    /* Build a Plus of Times-of-(coef, var^exp) terms. */
    struct Expr** terms = (struct Expr**)malloc(
        sizeof(struct Expr*) * p->n_terms);

    for (size_t i = 0; i < p->n_terms; i++) {
        const int* row = mpoly_exp_at((MPoly*)p, i);

        /* Build the variable factors. */
        size_t nfac = 0;
        for (int v = 0; v < p->n_vars; v++) if (row[v] > 0) nfac++;

        /* Coefficient handling: pure 1 omitted unless it's also the
         * only factor (then keep as integer 1). */
        bool coef_is_one = mpz_cmp_ui(p->coefs[i], 1) == 0;
        bool coef_is_neg_one = (mpz_cmp_si(p->coefs[i], -1) == 0);

        size_t total_factors = nfac;
        if (!(coef_is_one || coef_is_neg_one) || nfac == 0) total_factors++;
        if (coef_is_neg_one && nfac > 0) total_factors++;

        struct Expr** factors = (struct Expr**)malloc(
            sizeof(struct Expr*) * (total_factors > 0 ? total_factors : 1));
        size_t fi = 0;

        if (coef_is_neg_one && nfac > 0) {
            factors[fi++] = expr_new_integer(-1);
        } else if (!coef_is_one || nfac == 0) {
            /* Emit the literal coefficient. */
            if (mpz_fits_slong_p(p->coefs[i])) {
                factors[fi++] = expr_new_integer((int64_t)mpz_get_si(p->coefs[i]));
            } else {
                factors[fi++] = expr_new_bigint_from_mpz(p->coefs[i]);
            }
        }

        for (int v = 0; v < p->n_vars; v++) {
            if (row[v] == 0) continue;
            if (row[v] == 1) {
                factors[fi++] = expr_copy(vars[v]);
            } else {
                struct Expr** pa = (struct Expr**)malloc(sizeof(struct Expr*) * 2);
                pa[0] = expr_copy(vars[v]);
                pa[1] = expr_new_integer(row[v]);
                factors[fi++] = expr_new_function(
                    expr_new_symbol(SYM_Power), pa, 2);
            }
        }

        if (fi == 1) {
            terms[i] = factors[0];
            free(factors);
        } else {
            terms[i] = expr_new_function(expr_new_symbol(SYM_Times), factors, fi);
        }
    }

    if (p->n_terms == 1) {
        struct Expr* result = terms[0];
        free(terms);
        return result;
    }
    return expr_new_function(expr_new_symbol(SYM_Plus), terms, p->n_terms);
}
