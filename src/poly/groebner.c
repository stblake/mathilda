/* groebner.c
 *
 * Buchberger Gröbner-basis core.  See groebner.h for the data layout,
 * monomial-order semantics, and ownership rules.
 *
 * The algorithm is the textbook Buchberger ("An Algorithm for Computing
 * a Basis for a Polynomial Ideal", 1965) with Gebauer–Möller pair
 * criteria 1 (coprime leading monomials) and 2 (chain criterion); the
 * resulting basis is then interreduced and made monic, yielding the
 * reduced Gröbner basis -- the unique-up-to-permutation canonical form
 * Mathematica's GroebnerBasis returns.
 *
 * Coefficient arithmetic runs in GMP `mpq_t` throughout (no separate
 * denominator tracking).  This is slower than an `mpz_t` + content-
 * normalisation scheme on dense inputs, but keeps the implementation
 * straightforward and the inputs we accept here are typically small
 * by ideal-membership standards.
 */

#include "groebner.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"        /* tc_check_deadline() */
#include "expr.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ */
/*  Internal helpers                                                   */
/* ------------------------------------------------------------------ */

static inline int* gb_exp_at(const GBPoly* p, size_t i) {
    return p->exps + i * (size_t)p->n_vars;
}

/* Lex descending compare on exponent tuples (largest monomial first). */
static int cmp_lex_desc(const int* a, const int* b, int n_vars) {
    for (int i = 0; i < n_vars; i++) {
        if (a[i] != b[i]) return (a[i] > b[i]) ? -1 : 1;
    }
    return 0;
}

/* Graded reverse lex descending: total degree first, then break ties
 * by lex order *on the reversed exponent tuple* with the *smaller*
 * trailing exponent winning.  Equivalently: if total degrees match,
 * the monomial with the larger trailing exponent comes later. */
static int cmp_grevlex_desc(const int* a, const int* b, int n_vars) {
    int da = 0, db = 0;
    for (int i = 0; i < n_vars; i++) { da += a[i]; db += b[i]; }
    if (da != db) return (da > db) ? -1 : 1;
    for (int i = n_vars - 1; i >= 0; i--) {
        if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
    }
    return 0;
}

/* Elimination order: lex on the first `elim_pivot` variables; ties
 * broken by grevlex on the remaining tail.  This is well-founded and
 * eliminates the leading block whenever an elim-block-free polynomial
 * appears in the ideal. */
static int cmp_elim_desc(const int* a, const int* b, int n_vars, int piv) {
    for (int i = 0; i < piv; i++) {
        if (a[i] != b[i]) return (a[i] > b[i]) ? -1 : 1;
    }
    int da = 0, db = 0;
    for (int i = piv; i < n_vars; i++) { da += a[i]; db += b[i]; }
    if (da != db) return (da > db) ? -1 : 1;
    for (int i = n_vars - 1; i >= piv; i--) {
        if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
    }
    return 0;
}

/* Exact (mpz) fallback for a single weight-matrix row, used only when the
 * int64 accumulation in cmp_matrix_desc overflows. */
static int cmp_matrix_row_mpz(const int* a, const int* b,
                              const int64_t* row, int n_vars) {
    mpz_t wa, wb, t;
    mpz_init_set_ui(wa, 0); mpz_init_set_ui(wb, 0); mpz_init(t);
    for (int v = 0; v < n_vars; v++) {
        mpz_set_si(t, (long)row[v]); mpz_mul_si(t, t, (long)a[v]);
        mpz_add(wa, wa, t);
        mpz_set_si(t, (long)row[v]); mpz_mul_si(t, t, (long)b[v]);
        mpz_add(wb, wb, t);
    }
    int c = mpz_cmp(wa, wb);
    mpz_clear(wa); mpz_clear(wb); mpz_clear(t);
    if (c == 0) return 0;
    return (c > 0) ? -1 : 1;            /* larger weight first (descending) */
}

/* Weight-matrix order: rank monomials by the lexicographic comparison of
 * the integer weight vectors M*a vs M*b (larger weight = larger
 * monomial).  Each row weight accumulates in int64 with overflow
 * detection; on overflow that single row is recomputed exactly in mpz so
 * the comparison never silently wraps. */
static int cmp_matrix_desc(const int* a, const int* b, const GBWeightMatrix* m) {
    int n_vars = m->n_vars;
    for (int r = 0; r < m->n_rows; r++) {
        const int64_t* row = m->w + (size_t)r * (size_t)n_vars;
        int64_t wa = 0, wb = 0;
        bool overflow = false;
        for (int v = 0; v < n_vars; v++) {
            int64_t pa, pb;
            if (__builtin_mul_overflow(row[v], (int64_t)a[v], &pa)
             || __builtin_add_overflow(wa, pa, &wa)
             || __builtin_mul_overflow(row[v], (int64_t)b[v], &pb)
             || __builtin_add_overflow(wb, pb, &wb)) {
                overflow = true;
                break;
            }
        }
        int c;
        if (overflow)        c = cmp_matrix_row_mpz(a, b, row, n_vars);
        else if (wa != wb)   c = (wa > wb) ? -1 : 1;
        else                 c = 0;
        if (c != 0) return c;
    }
    return 0;
}

static int gb_cmp(const int* a, const int* b, const GBPoly* ctx) {
    switch (ctx->order) {
        case GB_ORDER_LEX:     return cmp_lex_desc(a, b, ctx->n_vars);
        case GB_ORDER_GREVLEX: return cmp_grevlex_desc(a, b, ctx->n_vars);
        case GB_ORDER_ELIM:    return cmp_elim_desc(a, b, ctx->n_vars,
                                                     ctx->elim_pivot);
        case GB_ORDER_MATRIX:  return cmp_matrix_desc(a, b, ctx->wmat);
    }
    return 0;
}

/* Validate that an n_rows x n_vars integer matrix defines a global term
 * order.  Two independent, both-required conditions (see groebner.h):
 *   (1) rank(M) == n_vars  -- injective on exponent vectors;
 *   (2) every column's first non-zero entry (top to bottom) is positive
 *       -- each variable exceeds 1, so the order is well-founded.
 * Rank is computed by exact rational Gaussian elimination (mpq); the
 * matrix is tiny so the cost is irrelevant. */
bool gb_wmat_validate(const int64_t* w, int n_rows, int n_vars) {
    if (!w || n_vars <= 0 || n_rows < n_vars) return false;

    /* Condition (2): first non-zero in each column must be > 0. */
    for (int c = 0; c < n_vars; c++) {
        for (int r = 0; r < n_rows; r++) {
            int64_t e = w[(size_t)r * n_vars + c];
            if (e != 0) {
                if (e < 0) return false;
                break;
            }
        }
    }

    /* Condition (1): rank == n_vars via rational row reduction. */
    mpq_t* M = (mpq_t*)malloc(sizeof(mpq_t) * (size_t)n_rows * n_vars);
    for (int i = 0; i < n_rows * n_vars; i++) {
        mpq_init(M[i]);
        mpq_set_si(M[i], (long)w[i], 1);
    }
    int rank = 0;
    mpq_t factor, tmp; mpq_init(factor); mpq_init(tmp);
    for (int col = 0; col < n_vars && rank < n_rows; col++) {
        /* Find a pivot row at or below `rank` with non-zero entry in col. */
        int piv = -1;
        for (int r = rank; r < n_rows; r++) {
            if (mpq_sgn(M[(size_t)r * n_vars + col]) != 0) { piv = r; break; }
        }
        if (piv < 0) continue;
        /* Swap pivot row into position `rank`. */
        if (piv != rank) {
            for (int c = 0; c < n_vars; c++) {
                mpq_swap(M[(size_t)rank * n_vars + c],
                         M[(size_t)piv * n_vars + c]);
            }
        }
        /* Eliminate this column from all other rows. */
        for (int r = 0; r < n_rows; r++) {
            if (r == rank) continue;
            mpq_t* rc = &M[(size_t)r * n_vars + col];
            if (mpq_sgn(*rc) == 0) continue;
            mpq_div(factor, *rc, M[(size_t)rank * n_vars + col]);
            for (int c = 0; c < n_vars; c++) {
                mpq_mul(tmp, factor, M[(size_t)rank * n_vars + c]);
                mpq_sub(M[(size_t)r * n_vars + c],
                        M[(size_t)r * n_vars + c], tmp);
            }
        }
        rank++;
    }
    mpq_clear(factor); mpq_clear(tmp);
    for (int i = 0; i < n_rows * n_vars; i++) mpq_clear(M[i]);
    free(M);
    return rank == n_vars;
}

/* qsort context.  Single-threaded; same convention as MPoly. */
static const GBPoly* g_sort_ctx = NULL;
static int sort_cmp_idx(const void* a, const void* b) {
    size_t ia = *(const size_t*)a;
    size_t ib = *(const size_t*)b;
    int c = gb_cmp(gb_exp_at(g_sort_ctx, ia),
                   gb_exp_at(g_sort_ctx, ib),
                   g_sort_ctx);
    if (c != 0) return c;
    return (ia < ib) ? -1 : (ia > ib);
}

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                         */
/* ------------------------------------------------------------------ */

GBPoly* gb_poly_new(int n_vars, GBOrder order, int elim_pivot) {
    GBPoly* p = (GBPoly*)malloc(sizeof(GBPoly));
    p->n_vars = n_vars;
    p->order = order;
    p->elim_pivot = elim_pivot;
    p->wmat = NULL;
    p->exps = NULL;
    p->coefs = NULL;
    p->n_terms = 0;
    p->cap = 0;
    return p;
}

void gb_poly_reserve(GBPoly* p, size_t cap) {
    if (cap <= p->cap) return;
    size_t new_cap = p->cap == 0 ? 8 : p->cap;
    while (new_cap < cap) new_cap *= 2;
    if (p->n_vars > 0) {
        p->exps = (int*)realloc(p->exps,
                                sizeof(int) * new_cap * (size_t)p->n_vars);
    }
    mpq_t* new_coefs = (mpq_t*)realloc(p->coefs, sizeof(mpq_t) * new_cap);
    for (size_t i = p->cap; i < new_cap; i++) mpq_init(new_coefs[i]);
    p->coefs = new_coefs;
    p->cap = new_cap;
}

void gb_poly_set_wmat(GBPoly* p, const GBWeightMatrix* wmat) {
    p->order = GB_ORDER_MATRIX;
    p->wmat = wmat;
    /* term order changed; caller must gb_poly_normalize() to re-sort */
}

GBPoly* gb_poly_copy(const GBPoly* p) {
    GBPoly* q = gb_poly_new(p->n_vars, p->order, p->elim_pivot);
    q->wmat = p->wmat;
    if (p->n_terms == 0) return q;
    gb_poly_reserve(q, p->n_terms);
    if (p->n_vars > 0) {
        memcpy(q->exps, p->exps,
               sizeof(int) * p->n_terms * (size_t)p->n_vars);
    }
    for (size_t i = 0; i < p->n_terms; i++) mpq_set(q->coefs[i], p->coefs[i]);
    q->n_terms = p->n_terms;
    return q;
}

void gb_poly_free(GBPoly* p) {
    if (!p) return;
    if (p->coefs) {
        for (size_t i = 0; i < p->cap; i++) mpq_clear(p->coefs[i]);
        free(p->coefs);
    }
    free(p->exps);
    free(p);
}

void gb_poly_push_term(GBPoly* p, const int* exps, const mpq_t coef) {
    if (mpz_sgn(mpq_numref(coef)) == 0) return;
    if (p->n_terms == p->cap) gb_poly_reserve(p, p->cap == 0 ? 8 : p->cap * 2);
    if (p->n_vars > 0) {
        memcpy(gb_exp_at(p, p->n_terms), exps,
               sizeof(int) * (size_t)p->n_vars);
    }
    mpq_set(p->coefs[p->n_terms], coef);
    p->n_terms++;
}

void gb_poly_push_term_si(GBPoly* p, const int* exps, int64_t num, int64_t den) {
    if (num == 0) return;
    if (p->n_terms == p->cap) gb_poly_reserve(p, p->cap == 0 ? 8 : p->cap * 2);
    if (p->n_vars > 0) {
        memcpy(gb_exp_at(p, p->n_terms), exps,
               sizeof(int) * (size_t)p->n_vars);
    }
    mpz_set_si(mpq_numref(p->coefs[p->n_terms]), (long)num);
    mpz_set_si(mpq_denref(p->coefs[p->n_terms]), (long)den);
    mpq_canonicalize(p->coefs[p->n_terms]);
    p->n_terms++;
}

void gb_poly_normalize(GBPoly* p) {
    if (p->n_terms < 2) {
        if (p->n_terms == 1 && mpz_sgn(mpq_numref(p->coefs[0])) == 0) {
            p->n_terms = 0;
        }
        return;
    }
    size_t* idx = (size_t*)malloc(sizeof(size_t) * p->n_terms);
    for (size_t i = 0; i < p->n_terms; i++) idx[i] = i;

    g_sort_ctx = p;
    qsort(idx, p->n_terms, sizeof(size_t), sort_cmp_idx);
    g_sort_ctx = NULL;

    int* new_exps = (p->n_vars > 0)
                    ? (int*)malloc(sizeof(int) * p->n_terms * (size_t)p->n_vars)
                    : NULL;
    mpq_t* new_coefs = (mpq_t*)malloc(sizeof(mpq_t) * p->n_terms);
    for (size_t i = 0; i < p->n_terms; i++) mpq_init(new_coefs[i]);

    size_t out = 0;
    for (size_t i = 0; i < p->n_terms; i++) {
        size_t src = idx[i];
        if (out > 0 && p->n_vars > 0 &&
            gb_cmp(new_exps + (out - 1) * (size_t)p->n_vars,
                   gb_exp_at(p, src),
                   p) == 0) {
            mpq_add(new_coefs[out - 1], new_coefs[out - 1], p->coefs[src]);
            if (mpz_sgn(mpq_numref(new_coefs[out - 1])) == 0) {
                /* Collapsed to zero -- back the cursor up. */
                out--;
                mpq_set_ui(new_coefs[out], 0, 1);
            }
        } else {
            if (p->n_vars > 0) {
                memcpy(new_exps + out * (size_t)p->n_vars,
                       gb_exp_at(p, src),
                       sizeof(int) * (size_t)p->n_vars);
            }
            mpq_set(new_coefs[out], p->coefs[src]);
            out++;
        }
    }

    /* Tear down the old arrays and install the new ones. */
    for (size_t i = 0; i < p->cap; i++) mpq_clear(p->coefs[i]);
    free(p->coefs);
    free(p->exps);

    p->exps = new_exps;
    p->coefs = new_coefs;
    p->cap = p->n_terms;
    p->n_terms = out;

    free(idx);
}

/* ------------------------------------------------------------------ */
/*  Predicates and queries                                             */
/* ------------------------------------------------------------------ */

bool gb_poly_is_zero(const GBPoly* p) {
    return p->n_terms == 0;
}

bool gb_poly_is_constant(const GBPoly* p) {
    if (p->n_terms == 0) return true;       /* zero is constant */
    if (p->n_terms > 1) return false;
    const int* row = gb_exp_at(p, 0);
    for (int v = 0; v < p->n_vars; v++) {
        if (row[v] != 0) return false;
    }
    return true;
}

const int* gb_poly_lm(const GBPoly* p) {
    if (p->n_terms == 0) return NULL;
    return gb_exp_at(p, 0);
}

const mpq_t* gb_poly_lc(const GBPoly* p) {
    if (p->n_terms == 0) return NULL;
    return (const mpq_t*)&p->coefs[0];
}

bool gb_poly_free_of_vars(const GBPoly* p, const int* vars, int n) {
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* row = gb_exp_at(p, t);
        for (int k = 0; k < n; k++) {
            if (row[vars[k]] != 0) return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Arithmetic                                                         */
/* ------------------------------------------------------------------ */

GBPoly* gb_poly_neg(const GBPoly* a) {
    GBPoly* r = gb_poly_copy(a);
    for (size_t i = 0; i < r->n_terms; i++) mpq_neg(r->coefs[i], r->coefs[i]);
    return r;
}

GBPoly* gb_poly_add(const GBPoly* a, const GBPoly* b) {
    GBPoly* r = gb_poly_new(a->n_vars, a->order, a->elim_pivot);
    r->wmat = a->wmat;
    gb_poly_reserve(r, a->n_terms + b->n_terms);
    /* Push everything and let normalize do the merge work. */
    for (size_t i = 0; i < a->n_terms; i++) {
        gb_poly_push_term(r, gb_exp_at(a, i), a->coefs[i]);
    }
    for (size_t i = 0; i < b->n_terms; i++) {
        gb_poly_push_term(r, gb_exp_at(b, i), b->coefs[i]);
    }
    gb_poly_normalize(r);
    return r;
}

GBPoly* gb_poly_sub(const GBPoly* a, const GBPoly* b) {
    GBPoly* nb = gb_poly_neg(b);
    GBPoly* r = gb_poly_add(a, nb);
    gb_poly_free(nb);
    return r;
}

GBPoly* gb_poly_scale(const GBPoly* a, const mpq_t c) {
    GBPoly* r = gb_poly_copy(a);
    if (mpz_sgn(mpq_numref(c)) == 0) { r->n_terms = 0; return r; }
    for (size_t i = 0; i < r->n_terms; i++) mpq_mul(r->coefs[i], r->coefs[i], c);
    return r;
}

GBPoly* gb_poly_mul_by_monomial(const GBPoly* a, const int* exps, const mpq_t c) {
    GBPoly* r = gb_poly_new(a->n_vars, a->order, a->elim_pivot);
    r->wmat = a->wmat;
    if (mpz_sgn(mpq_numref(c)) == 0) return r;
    gb_poly_reserve(r, a->n_terms);
    int* tmp = (a->n_vars > 0)
               ? (int*)malloc(sizeof(int) * (size_t)a->n_vars)
               : NULL;
    mpq_t tcoef; mpq_init(tcoef);
    for (size_t i = 0; i < a->n_terms; i++) {
        const int* row = gb_exp_at(a, i);
        for (int v = 0; v < a->n_vars; v++) tmp[v] = row[v] + exps[v];
        mpq_mul(tcoef, a->coefs[i], c);
        gb_poly_push_term(r, tmp, tcoef);
    }
    mpq_clear(tcoef);
    free(tmp);
    /* Multiplying by a single monomial preserves the relative order of
     * existing terms under every supported monomial order, so no need
     * to re-sort.  (Each LM(a_i) gets shifted by the same exponent.) */
    return r;
}

void gb_poly_make_monic(GBPoly* p) {
    if (p->n_terms == 0) return;
    if (mpq_cmp_ui(p->coefs[0], 1, 1) == 0) return;
    mpq_t inv; mpq_init(inv);
    mpq_inv(inv, p->coefs[0]);
    for (size_t i = 0; i < p->n_terms; i++) {
        mpq_mul(p->coefs[i], p->coefs[i], inv);
    }
    mpq_clear(inv);
}

/* Scale `p` so that all coefficients become integers with overall
 * content (GCD of numerators) equal to 1 and the leading coefficient
 * positive.  This is the Mathematica convention for Gröbner basis
 * output. */
static void gb_poly_make_primitive_z(GBPoly* p) {
    if (p->n_terms == 0) return;
    /* Step 1: multiply by LCM of denominators to clear fractions. */
    mpz_t lcm; mpz_init_set_ui(lcm, 1);
    for (size_t i = 0; i < p->n_terms; i++) {
        mpz_lcm(lcm, lcm, mpq_denref(p->coefs[i]));
    }
    if (mpz_cmp_ui(lcm, 1) != 0) {
        for (size_t i = 0; i < p->n_terms; i++) {
            mpz_mul(mpq_numref(p->coefs[i]), mpq_numref(p->coefs[i]), lcm);
            mpq_canonicalize(p->coefs[i]);
        }
    }
    mpz_clear(lcm);
    /* Step 2: divide by GCD of numerators to strip content. */
    mpz_t g; mpz_init_set(g, mpq_numref(p->coefs[0]));
    mpz_abs(g, g);
    for (size_t i = 1; i < p->n_terms; i++) {
        mpz_gcd(g, g, mpq_numref(p->coefs[i]));
        if (mpz_cmp_ui(g, 1) == 0) break;
    }
    if (mpz_cmp_ui(g, 1) != 0) {
        for (size_t i = 0; i < p->n_terms; i++) {
            mpz_divexact(mpq_numref(p->coefs[i]), mpq_numref(p->coefs[i]), g);
        }
    }
    mpz_clear(g);
    /* Step 3: ensure leading coefficient is positive. */
    if (mpz_sgn(mpq_numref(p->coefs[0])) < 0) {
        for (size_t i = 0; i < p->n_terms; i++) {
            mpz_neg(mpq_numref(p->coefs[i]), mpq_numref(p->coefs[i]));
        }
    }
}

/* ------------------------------------------------------------------ */
/*  Buchberger primitives                                              */
/* ------------------------------------------------------------------ */

/* lcm of two exponent vectors -- pointwise max. */
static void lcm_exps(int* out, const int* a, const int* b, int n_vars) {
    for (int i = 0; i < n_vars; i++) out[i] = a[i] > b[i] ? a[i] : b[i];
}

GBPoly* gb_spoly(const GBPoly* f, const GBPoly* g) {
    assert(f->n_terms > 0 && g->n_terms > 0);
    int n = f->n_vars;
    int* lcm = (int*)malloc(sizeof(int) * (size_t)n);
    int* af  = (int*)malloc(sizeof(int) * (size_t)n);
    int* ag  = (int*)malloc(sizeof(int) * (size_t)n);
    const int* lmf = gb_exp_at(f, 0);
    const int* lmg = gb_exp_at(g, 0);
    lcm_exps(lcm, lmf, lmg, n);
    for (int i = 0; i < n; i++) {
        af[i] = lcm[i] - lmf[i];
        ag[i] = lcm[i] - lmg[i];
    }
    mpq_t cf, cg, inv;
    mpq_init(cf); mpq_init(cg); mpq_init(inv);
    /* cf = 1 / lc(f), cg = 1 / lc(g).  S = (lcm/lt(f))*f - (lcm/lt(g))*g
     * with all leading scaling baked in. */
    mpq_inv(cf, f->coefs[0]);
    mpq_inv(cg, g->coefs[0]);
    GBPoly* tf = gb_poly_mul_by_monomial(f, af, cf);
    GBPoly* tg = gb_poly_mul_by_monomial(g, ag, cg);
    GBPoly* s  = gb_poly_sub(tf, tg);
    gb_poly_free(tf);
    gb_poly_free(tg);
    mpq_clear(cf); mpq_clear(cg); mpq_clear(inv);
    free(lcm); free(af); free(ag);
    return s;
}

/* Is the exponent vector `a` divisible by `b` (i.e., a >= b component-
 * wise)?  If yes, write the quotient exponents into `q`. */
static bool divides(const int* b, const int* a, int* q, int n_vars) {
    for (int i = 0; i < n_vars; i++) {
        if (a[i] < b[i]) return false;
        q[i] = a[i] - b[i];
    }
    return true;
}

GBPoly* gb_reduce(const GBPoly* p, GBPoly* const* basis, size_t n) {
    GBPoly* r = gb_poly_copy(p);
    if (r->n_terms == 0) return r;
    int n_vars = r->n_vars;
    int* qexp = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
    mpq_t qc; mpq_init(qc);

    /* Standard multivariate division: repeatedly look for any term of
     * `r` whose monomial is divisible by some basis leading monomial;
     * subtract the appropriate multiple.  Continue until no further
     * reduction is possible. */
    bool reduced;
    do {
        reduced = false;
        for (size_t t = 0; t < r->n_terms; t++) {
            const int* re = gb_exp_at(r, t);
            for (size_t bi = 0; bi < n; bi++) {
                const GBPoly* g = basis[bi];
                if (g->n_terms == 0) continue;
                const int* lm = gb_exp_at(g, 0);
                if (!divides(lm, re, qexp, n_vars)) continue;
                /* qc = r.coef[t] / lc(g) */
                mpq_div(qc, r->coefs[t], g->coefs[0]);
                GBPoly* sub = gb_poly_mul_by_monomial(g, qexp, qc);
                GBPoly* nr  = gb_poly_sub(r, sub);
                gb_poly_free(sub);
                gb_poly_free(r);
                r = nr;
                reduced = true;
                goto restart;
            }
        }
restart: ;
    } while (reduced);

    mpq_clear(qc);
    free(qexp);
    return r;
}

GBPoly* gb_divmod(const GBPoly* p, GBPoly* const* basis, size_t n,
                  GBPoly*** quot_out) {
    GBPoly* r = gb_poly_copy(p);
    int n_vars = r->n_vars;

    /* One quotient accumulator per basis element (terms pushed in
     * arbitrary order, merged by normalize at the end). */
    GBPoly** quot = (GBPoly**)malloc(sizeof(GBPoly*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        quot[i] = gb_poly_new(n_vars, p->order, p->elim_pivot);
        quot[i]->wmat = p->wmat;
    }

    if (r->n_terms == 0) { *quot_out = quot; return r; }

    int* qexp = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
    mpq_t qc; mpq_init(qc);

    /* Same full multivariate division as gb_reduce, but each subtracted
     * multiple  qc * x^qexp * basis[bi]  is also recorded on quot[bi], so
     * that on return  p == sum_i quot[i]*basis[i] + r. */
    bool reduced;
    do {
        reduced = false;
        for (size_t t = 0; t < r->n_terms; t++) {
            const int* re = gb_exp_at(r, t);
            for (size_t bi = 0; bi < n; bi++) {
                const GBPoly* g = basis[bi];
                if (g->n_terms == 0) continue;
                const int* lm = gb_exp_at(g, 0);
                if (!divides(lm, re, qexp, n_vars)) continue;
                mpq_div(qc, r->coefs[t], g->coefs[0]);
                gb_poly_push_term(quot[bi], qexp, qc);
                GBPoly* sub = gb_poly_mul_by_monomial(g, qexp, qc);
                GBPoly* nr  = gb_poly_sub(r, sub);
                gb_poly_free(sub);
                gb_poly_free(r);
                r = nr;
                reduced = true;
                goto restart;
            }
        }
restart: ;
    } while (reduced);

    mpq_clear(qc);
    free(qexp);
    for (size_t i = 0; i < n; i++) gb_poly_normalize(quot[i]);
    *quot_out = quot;
    return r;
}

GBPoly* gb_initial_form(const GBPoly* g, const int64_t* w, int n_vars) {
    assert(g->n_terms > 0);
    GBPoly* r = gb_poly_new(n_vars, g->order, g->elim_pivot);
    r->wmat = g->wmat;

    /* Weight of a term, exact: int64 with mpz fallback on overflow. */
    mpz_t best, cur, prod;
    mpz_init(best); mpz_init(cur); mpz_init(prod);
    bool have_best = false;
    for (size_t t = 0; t < g->n_terms; t++) {
        const int* e = gb_exp_at(g, t);
        mpz_set_ui(cur, 0);
        for (int v = 0; v < n_vars; v++) {
            mpz_set_si(prod, (long)w[v]);
            mpz_mul_si(prod, prod, (long)e[v]);
            mpz_add(cur, cur, prod);
        }
        if (!have_best || mpz_cmp(cur, best) > 0) {
            mpz_set(best, cur);
            have_best = true;
        }
    }
    /* Collect every term whose weight equals the maximum. */
    for (size_t t = 0; t < g->n_terms; t++) {
        const int* e = gb_exp_at(g, t);
        mpz_set_ui(cur, 0);
        for (int v = 0; v < n_vars; v++) {
            mpz_set_si(prod, (long)w[v]);
            mpz_mul_si(prod, prod, (long)e[v]);
            mpz_add(cur, cur, prod);
        }
        if (mpz_cmp(cur, best) == 0) gb_poly_push_term(r, e, g->coefs[t]);
    }
    mpz_clear(best); mpz_clear(cur); mpz_clear(prod);
    gb_poly_normalize(r);
    return r;
}

/* ----- Gebauer–Möller pair management ----- */

/* Each pending S-pair carries pre-computed LCM data so we can:
 *   - rank pairs by "sugar"-free normal strategy (smallest total LCM
 *     degree first; ties broken by lex order of the LCM exponents),
 *   - apply Buchberger criterion 1 (coprime leading monomials) and
 *     Gebauer-Möller criteria M / F / B without re-walking exponents
 *     during every comparison.
 * `lcm` is heap-allocated alongside the pair; `dead` flags eliminated
 * pairs that are compacted out periodically. */
typedef struct {
    size_t i, j;
    int*   lcm;           /* length = n_vars */
    int    total_deg;
    bool   dead;
} GMPair;

/* Coprime LM test: do the leading monomials share no variable?  If so,
 * the S-polynomial reduces to 0 (Buchberger first criterion). */
static bool lm_coprime(const GBPoly* f, const GBPoly* g) {
    const int* a = gb_exp_at(f, 0);
    const int* b = gb_exp_at(g, 0);
    int n = f->n_vars;
    for (int v = 0; v < n; v++) {
        if (a[v] > 0 && b[v] > 0) return false;
    }
    return true;
}

/* lcm of two exponent vectors -- pointwise max.  Returns total degree. */
static int compute_lcm(int* out, const int* a, const int* b, int n_vars) {
    int tot = 0;
    for (int v = 0; v < n_vars; v++) {
        int m = a[v] > b[v] ? a[v] : b[v];
        out[v] = m;
        tot += m;
    }
    return tot;
}

/* Does exponent vector `b` divide `a`?  All-component <= test. */
static bool exp_divides(const int* b, const int* a, int n_vars) {
    for (int v = 0; v < n_vars; v++) if (a[v] < b[v]) return false;
    return true;
}

/* Are two exponent vectors equal? */
static bool exp_equal(const int* a, const int* b, int n_vars) {
    for (int v = 0; v < n_vars; v++) if (a[v] != b[v]) return false;
    return true;
}

static GMPair gm_pair_new(size_t i, size_t j, const GBPoly* gi,
                          const GBPoly* gj, int n_vars) {
    GMPair p;
    p.i = i; p.j = j;
    p.lcm = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
    p.total_deg = compute_lcm(p.lcm, gb_exp_at(gi, 0), gb_exp_at(gj, 0), n_vars);
    p.dead = false;
    return p;
}

static void gm_pair_free(GMPair* p) {
    free(p->lcm); p->lcm = NULL;
}

/* Compact the pair list, discarding dead entries (in place). */
static size_t gm_compact(GMPair* pairs, size_t nP) {
    size_t w = 0;
    for (size_t i = 0; i < nP; i++) {
        if (pairs[i].dead) {
            gm_pair_free(&pairs[i]);
        } else {
            if (w != i) pairs[w] = pairs[i];
            w++;
        }
    }
    return w;
}

/* Normal-strategy pair selection: index of the live pair with smallest
 * total LCM degree (ties broken by lex of LCM exponents, then by pair
 * order to keep results deterministic).  Returns SIZE_MAX if no live
 * pair exists.  O(nP * n_vars) per call -- the cost is dominated by
 * the S-poly reduction so a heap would not pay back its complexity. */
static size_t gm_pick_pair(const GMPair* pairs, size_t nP, int n_vars) {
    size_t best = (size_t)-1;
    for (size_t k = 0; k < nP; k++) {
        if (pairs[k].dead) continue;
        if (best == (size_t)-1) { best = k; continue; }
        if (pairs[k].total_deg != pairs[best].total_deg) {
            if (pairs[k].total_deg < pairs[best].total_deg) best = k;
            continue;
        }
        /* Lex tiebreak: smaller LCM first. */
        for (int v = 0; v < n_vars; v++) {
            if (pairs[k].lcm[v] != pairs[best].lcm[v]) {
                if (pairs[k].lcm[v] < pairs[best].lcm[v]) best = k;
                break;
            }
        }
    }
    return best;
}

/* Gebauer-Möller UPDATE: register a freshly-added basis element h
 * (= G[h_idx]) and refresh the pair list `pairs` / `nP` accordingly.
 * Implements:
 *   - Buchberger criterion 1: candidate pair (h, g) with coprime LMs
 *     is dropped immediately (the S-poly is known to reduce to 0).
 *   - Gebauer-Möller M criterion: among candidates (h, g), drop
 *     (h, g') if some (h, g) has lcm(LM(h), LM(g)) strictly dividing
 *     lcm(LM(h), LM(g')).
 *   - Gebauer-Möller B criterion: drop existing pair (a, b) when
 *     LM(h) divides lcm(LM(a), LM(b)) and the two "lcm-equal" exits
 *     do not apply -- one of the equivalent pairs (a, h) or (b, h)
 *     handles the redundancy.
 * `*capP` is grown in place. */
static void gm_update(GMPair** pairs_p, size_t* nP_p, size_t* capP_p,
                      GBPoly* const* G, size_t h_idx, int n_vars) {
    GMPair* pairs = *pairs_p;
    size_t nP = *nP_p;
    size_t capP = *capP_p;
    const GBPoly* h = G[h_idx];

    /* Step 1: form candidates (h, g) for every g already in G. */
    GMPair* cand = (GMPair*)malloc(sizeof(GMPair) * (h_idx ? h_idx : 1));
    size_t nc = 0;
    for (size_t g = 0; g < h_idx; g++) {
        GMPair p = gm_pair_new(g, h_idx, G[g], h, n_vars);
        /* Criterion 1: skip coprime LMs entirely. */
        if (lm_coprime(G[g], h)) { gm_pair_free(&p); continue; }
        cand[nc++] = p;
    }

    /* Step 2: M criterion -- if some candidate has LCM strictly
     * dividing another's, mark the larger one dead. */
    for (size_t a = 0; a < nc; a++) {
        if (cand[a].dead) continue;
        for (size_t b = 0; b < nc; b++) {
            if (a == b || cand[b].dead) continue;
            if (exp_divides(cand[a].lcm, cand[b].lcm, n_vars)
                && !exp_equal(cand[a].lcm, cand[b].lcm, n_vars)) {
                cand[b].dead = true;
            }
        }
    }
    /* F criterion (de-duplicate LCM-equal pairs): keep the first. */
    for (size_t a = 0; a < nc; a++) {
        if (cand[a].dead) continue;
        for (size_t b = a + 1; b < nc; b++) {
            if (cand[b].dead) continue;
            if (exp_equal(cand[a].lcm, cand[b].lcm, n_vars)) {
                cand[b].dead = true;
            }
        }
    }

    /* Step 3: B criterion -- drop existing pairs (a, b) whose LCM is
     * divisible by LM(h) AND for which (a, h) and (b, h) are not
     * "lcm-equal".  In our notation: drop (a, b) if LM(h) | lcm(a, b)
     * and lcm(a, b) != lcm(a, h) and lcm(a, b) != lcm(b, h). */
    const int* lm_h = gb_exp_at(h, 0);
    for (size_t k = 0; k < nP; k++) {
        if (pairs[k].dead) continue;
        if (!exp_divides(lm_h, pairs[k].lcm, n_vars)) continue;
        /* lcm(a, h) and lcm(b, h) */
        int* lcm_ah = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
        int* lcm_bh = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
        compute_lcm(lcm_ah, gb_exp_at(G[pairs[k].i], 0), lm_h, n_vars);
        compute_lcm(lcm_bh, gb_exp_at(G[pairs[k].j], 0), lm_h, n_vars);
        bool ah_eq = exp_equal(pairs[k].lcm, lcm_ah, n_vars);
        bool bh_eq = exp_equal(pairs[k].lcm, lcm_bh, n_vars);
        if (!ah_eq && !bh_eq) pairs[k].dead = true;
        free(lcm_ah); free(lcm_bh);
    }

    /* Compact + append surviving candidates. */
    nP = gm_compact(pairs, nP);
    size_t need = nP + nc;
    if (need > capP) {
        while (capP < need) capP = capP ? capP * 2 : 16;
        pairs = (GMPair*)realloc(pairs, sizeof(GMPair) * capP);
    }
    for (size_t a = 0; a < nc; a++) {
        if (cand[a].dead) { gm_pair_free(&cand[a]); continue; }
        pairs[nP++] = cand[a];
    }
    free(cand);

    *pairs_p = pairs;
    *nP_p = nP;
    *capP_p = capP;
}

/* Turn a generating set G[0..*nG_io-1] of an ideal (already a Gröbner
 * basis under G[0]'s order, but not yet reduced) into the canonical
 * reduced Gröbner basis: drop elements whose LM is divisible by another's,
 * reduce each survivor by the rest to a fixed point, make each primitive
 * over Z, and sort ascending by leading monomial.  Operates in place; the
 * array is only compacted (never grown) and *nG_io is updated.  Shared by
 * gb_buchberger and the Gröbner walk. */
void gb_finalize_basis(GBPoly** G, size_t* nG_io) {
    size_t nG = *nG_io;
    if (nG == 0) { *nG_io = 0; return; }
    int n_vars = G[0]->n_vars;

    /* Interreduce: discard any g_i whose LM is divisible by LM(g_k)
     * for some k != i, then reduce each survivor by the others. */
    bool* keep = (bool*)malloc(sizeof(bool) * nG);
    for (size_t i = 0; i < nG; i++) keep[i] = true;
    for (size_t i = 0; i < nG; i++) {
        if (!keep[i]) continue;
        const int* lmi = gb_exp_at(G[i], 0);
        for (size_t k = 0; k < nG; k++) {
            if (k == i || !keep[k]) continue;
            const int* lmk = gb_exp_at(G[k], 0);
            bool div_ok = true;
            for (int v = 0; v < n_vars; v++) {
                if (lmi[v] < lmk[v]) { div_ok = false; break; }
            }
            if (div_ok) { keep[i] = false; break; }
        }
    }

    /* Compact in place. */
    size_t out = 0;
    for (size_t i = 0; i < nG; i++) {
        if (keep[i]) {
            G[out++] = G[i];
        } else {
            gb_poly_free(G[i]);
        }
    }
    nG = out;
    free(keep);

    /* Reduce each survivor by the others.  Iterate to a fixed point in
     * case reductions create cross-references. */
    bool changed;
    do {
        changed = false;
        for (size_t i = 0; i < nG; i++) {
            /* Build a temporary basis "all except i". */
            GBPoly** other = (GBPoly**)malloc(sizeof(GBPoly*) * (nG ? nG : 1));
            size_t m = 0;
            for (size_t k = 0; k < nG; k++) if (k != i) other[m++] = G[k];
            GBPoly* nr = gb_reduce(G[i], other, m);
            gb_poly_make_monic(nr);
            free(other);
            /* If the reduction collapsed to zero, drop g_i.  (Should be
             * impossible for a reduced Gröbner basis after the LM-
             * filtering pass, but be defensive.) */
            if (nr->n_terms == 0) {
                gb_poly_free(nr);
                gb_poly_free(G[i]);
                for (size_t k = i; k + 1 < nG; k++) G[k] = G[k + 1];
                nG--;
                i--;
                changed = true;
                continue;
            }
            /* If a strict change happened, install. */
            bool same = (nr->n_terms == G[i]->n_terms);
            if (same) {
                for (size_t t = 0; t < nr->n_terms && same; t++) {
                    if (mpq_cmp(nr->coefs[t], G[i]->coefs[t]) != 0) same = false;
                    if (same && n_vars > 0 && memcmp(gb_exp_at(nr, t),
                                                     gb_exp_at(G[i], t),
                                                     sizeof(int) * (size_t)n_vars)
                                                != 0) same = false;
                }
            }
            if (!same) {
                gb_poly_free(G[i]);
                G[i] = nr;
                changed = true;
            } else {
                gb_poly_free(nr);
            }
        }
    } while (changed);

    /* Convert each basis polynomial to Mathematica's primitive-over-Z
     * form so the output matches its conventions. */
    for (size_t i = 0; i < nG; i++) gb_poly_make_primitive_z(G[i]);

    /* Sort the basis ascending by leading monomial (Mathematica
     * convention: in Lex {x, y, z} the z-only polynomial appears
     * first, the x-leading polynomial last). */
    for (size_t i = 0; i + 1 < nG; i++) {
        size_t pick = i;
        for (size_t j = i + 1; j < nG; j++) {
            if (gb_cmp(gb_exp_at(G[j], 0), gb_exp_at(G[pick], 0), G[0]) > 0) {
                /* gb_cmp returns -1 when a > b (descending sense).
                 * We want ascending: pick the smaller LM. */
                pick = j;
            }
        }
        if (pick != i) { GBPoly* tmp = G[i]; G[i] = G[pick]; G[pick] = tmp; }
    }

    *nG_io = nG;
}

GBPoly** gb_buchberger(GBPoly* const* F, size_t n, size_t* out_n) {
    /* Filter zeros from input, copy the rest into the working basis. */
    GBPoly** G = NULL;
    size_t nG = 0, capG = 0;
    for (size_t i = 0; i < n; i++) {
        if (F[i]->n_terms == 0) continue;
        if (nG == capG) {
            capG = capG ? capG * 2 : 8;
            G = (GBPoly**)realloc(G, sizeof(GBPoly*) * capG);
        }
        G[nG] = gb_poly_copy(F[i]);
        gb_poly_make_monic(G[nG]);
        nG++;
    }

    if (nG == 0) {
        *out_n = 0;
        free(G);
        return NULL;
    }

    int n_vars = G[0]->n_vars;

    /* Seed the pair set by calling gm_update for each element after the
     * first; this naturally applies the M/F/B criteria during the
     * initial pair generation as well as the incremental loop. */
    GMPair* pairs = NULL;
    size_t nP = 0, capP = 0;
    for (size_t h = 1; h < nG; h++) gm_update(&pairs, &nP, &capP, G, h, n_vars);

    /* Main loop: process pairs in normal-strategy order. */
    while (true) {
        size_t pick = gm_pick_pair(pairs, nP, n_vars);
        if (pick == (size_t)-1) break;

        /* Cooperative abort hook -- so TimeConstrained[GroebnerBasis[...],
         * t] siglongjmps out at the next pair instead of running to
         * completion on pathological inputs (the issue-2 hanging case). */
        tc_check_deadline();

        GMPair pr = pairs[pick];
        pairs[pick].dead = true;
        gm_pair_free(&pairs[pick]);

        GBPoly* s = gb_spoly(G[pr.i], G[pr.j]);
        GBPoly* r = gb_reduce(s, G, nG);
        gb_poly_free(s);
        if (r->n_terms == 0) { gb_poly_free(r); continue; }

        gb_poly_make_monic(r);

        if (nG == capG) {
            capG = capG ? capG * 2 : 8;
            G = (GBPoly**)realloc(G, sizeof(GBPoly*) * capG);
        }
        G[nG] = r;
        gm_update(&pairs, &nP, &capP, G, nG, n_vars);
        nG++;
    }

    for (size_t k = 0; k < nP; k++) {
        if (!pairs[k].dead) gm_pair_free(&pairs[k]);
    }
    free(pairs);

    gb_finalize_basis(G, &nG);

    *out_n = nG;
    return G;
}

void gb_basis_free(GBPoly** basis, size_t n) {
    if (!basis) return;
    for (size_t i = 0; i < n; i++) gb_poly_free(basis[i]);
    free(basis);
}

/* ------------------------------------------------------------------ */
/*  RationalFunctions coefficient domain                               */
/* ------------------------------------------------------------------ */
/*
 * A Gröbner basis over the field Q(params)[mainvars] is obtained from the
 * ring Gröbner basis over Q[params][mainvars].  The ring basis is computed
 * with the parameters placed LAST under lexicographic order, so a ring
 * leading term projects onto the field leading term (its main-variable
 * part); the localisation theorem then says the ring basis, autoreduced
 * with the parameter part of every leading coefficient treated as a unit,
 * is a Gröbner basis over the field.  Concretely, the field-redundant
 * generators -- those whose main leading monomial is a multiple of
 * another's -- are removed and surviving generators are tail-reduced.
 *
 * All arithmetic stays integral via PSEUDO-reduction: to cancel f's field
 * leading term using g (with LM_main(g) | LM_main(f)) we form
 *     f := lc_main(g) * f  -  (x^q * lc_main(f)) * g
 * where lc_main(.) is the leading coefficient as a polynomial in the
 * parameters (a unit in the field) and x^q the main-monomial quotient.
 * Each intermediate is made primitive over Z to curb coefficient growth.
 * The first `n_main` exponent slots are the field (main) variables; the
 * trailing slots are the parameters.
 */

/* True iff the main projection of monomial `a` divides that of `b`. */
static bool gb_main_divides(const int* a, const int* b, int n_main) {
    for (int i = 0; i < n_main; i++) if (a[i] > b[i]) return false;
    return true;
}

/* General polynomial product a*b (fresh, normalised). */
static GBPoly* gb_poly_mul(const GBPoly* a, const GBPoly* b) {
    GBPoly* out = gb_poly_new(a->n_vars, a->order, a->elim_pivot);
    out->wmat = a->wmat;
    for (size_t j = 0; j < b->n_terms; j++) {
        GBPoly* t = gb_poly_mul_by_monomial(a, gb_exp_at((GBPoly*)b, j),
                                            b->coefs[j]);
        GBPoly* s = gb_poly_add(out, t);
        gb_poly_free(out);
        gb_poly_free(t);
        out = s;
    }
    return out;
}

/* lc_main(p): the leading coefficient of non-zero `p` as a polynomial in
 * the parameters -- the prefix of terms sharing p's leading main monomial,
 * with the main exponents zeroed.  (Under lex with params last that group
 * is contiguous at the front.) */
static GBPoly* gb_lc_main(const GBPoly* p, int n_main) {
    GBPoly* out = gb_poly_new(p->n_vars, p->order, p->elim_pivot);
    out->wmat = p->wmat;
    const int* lm = gb_exp_at((GBPoly*)p, 0);
    int* row = (int*)malloc(sizeof(int) * (size_t)p->n_vars);
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* e = gb_exp_at((GBPoly*)p, i);
        bool same = true;
        for (int v = 0; v < n_main; v++) if (e[v] != lm[v]) { same = false; break; }
        if (!same) break;
        for (int v = 0; v < p->n_vars; v++) row[v] = (v < n_main) ? 0 : e[v];
        gb_poly_push_term(out, row, p->coefs[i]);
    }
    free(row);
    gb_poly_normalize(out);
    return out;
}

/* The leading main-monomial group of non-zero `p` with its monomials
 * intact (the field leading term, as a sub-polynomial of p). */
static GBPoly* gb_lead_group(const GBPoly* p, int n_main) {
    GBPoly* out = gb_poly_new(p->n_vars, p->order, p->elim_pivot);
    out->wmat = p->wmat;
    const int* lm = gb_exp_at((GBPoly*)p, 0);
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* e = gb_exp_at((GBPoly*)p, i);
        bool same = true;
        for (int v = 0; v < n_main; v++) if (e[v] != lm[v]) { same = false; break; }
        if (!same) break;
        gb_poly_push_term(out, e, p->coefs[i]);
    }
    gb_poly_normalize(out);
    return out;
}

/* Full field reduction of `f` modulo the non-NULL members of
 * basis[0..n-1].  Returns a fresh remainder (possibly zero); coefficients
 * are kept primitive over Z. */
static GBPoly* gb_field_reduce(const GBPoly* f, GBPoly* const* basis,
                               size_t n, int n_main) {
    GBPoly* r = gb_poly_copy(f);
    GBPoly* done = gb_poly_new(f->n_vars, f->order, f->elim_pivot);
    done->wmat = f->wmat;
    int* qexp = (int*)malloc(sizeof(int) * (size_t)f->n_vars);
    mpq_t one; mpq_init(one); mpq_set_ui(one, 1, 1);

    while (!gb_poly_is_zero(r)) {
        const int* lmr = gb_exp_at(r, 0);
        size_t gi = n;
        for (size_t i = 0; i < n; i++) {
            if (!basis[i] || gb_poly_is_zero(basis[i])) continue;
            if (gb_main_divides(gb_exp_at(basis[i], 0), lmr, n_main)) {
                gi = i; break;
            }
        }
        if (gi == n) {
            /* Leading main group is irreducible: retire it to `done`. */
            GBPoly* lead = gb_lead_group(r, n_main);
            GBPoly* rest = gb_poly_sub(r, lead);
            GBPoly* nd   = gb_poly_add(done, lead);
            gb_poly_free(done); done = nd;
            gb_poly_free(lead);
            gb_poly_free(r); r = rest;
            continue;
        }
        /* Pseudo-reduce: r := lc_main(g)*r - (x^q * lc_main(r)) * g. */
        const int* lmg = gb_exp_at(basis[gi], 0);
        for (int v = 0; v < f->n_vars; v++)
            qexp[v] = (v < n_main) ? (lmr[v] - lmg[v]) : 0;
        GBPoly* lcg = gb_lc_main(basis[gi], n_main);
        GBPoly* lcr = gb_lc_main(r, n_main);
        GBPoly* t1  = gb_poly_mul(lcg, r);
        GBPoly* tmp = gb_poly_mul(lcr, basis[gi]);
        GBPoly* t2  = gb_poly_mul_by_monomial(tmp, qexp, one);
        GBPoly* nr  = gb_poly_sub(t1, t2);
        gb_poly_make_primitive_z(nr);
        gb_poly_free(lcg); gb_poly_free(lcr);
        gb_poly_free(t1);  gb_poly_free(tmp); gb_poly_free(t2);
        gb_poly_free(r); r = nr;
    }
    mpq_clear(one);
    free(qexp);
    gb_poly_free(r);
    return done;
}

/* Autoreduce the ring basis G[0..*nG-1] into a Gröbner basis over the
 * field Q(params)[mainvars], in place.  `n_main` is the number of leading
 * (field) variables; the trailing slots are parameters.  NULL entries are
 * compacted away and *nG updated; each survivor is primitive over Z. */
void gb_rational_function_reduce(GBPoly** G, size_t* nG, int n_main) {
    /* Fixed-point autoreduction with a generous safety cap. */
    for (int pass = 0; pass < 1000; pass++) {
        bool changed = false;
        for (size_t i = 0; i < *nG; i++) {
            if (!G[i]) continue;
            /* Reduce G[i] modulo the other (non-NULL) generators. */
            GBPoly** others = (GBPoly**)malloc(sizeof(GBPoly*) * *nG);
            size_t m = 0;
            for (size_t j = 0; j < *nG; j++)
                if (j != i && G[j]) others[m++] = G[j];
            GBPoly* r = gb_field_reduce(G[i], others, m, n_main);
            free(others);

            if (gb_poly_is_zero(r)) {
                gb_poly_free(G[i]); G[i] = NULL;
                gb_poly_free(r);
                changed = true;
            } else {
                gb_poly_make_primitive_z(r);
                GBPoly* d = gb_poly_sub(G[i], r);
                bool same = gb_poly_is_zero(d);
                gb_poly_free(d);
                if (same) {
                    gb_poly_free(r);
                } else {
                    gb_poly_free(G[i]); G[i] = r;
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    /* Compact out the removed generators. */
    size_t k = 0;
    for (size_t i = 0; i < *nG; i++) if (G[i]) G[k++] = G[i];
    *nG = k;
}

/* ------------------------------------------------------------------ */
/*  Expr round-trip                                                    */
/* ------------------------------------------------------------------ */

/* Match a known variable. */
static int find_var_index(struct Expr* e, struct Expr** vars, int n_vars) {
    for (int i = 0; i < n_vars; i++) {
        if (expr_eq(e, vars[i])) return i;
    }
    return -1;
}

/* Interpret e as a rational scalar (Integer, BigInt, Rational[a,b]).
 * Returns true and sets `out` on success.  `out` must already be init'd. */
static bool expr_to_mpq(struct Expr* e, mpq_t out) {
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1);
        mpq_canonicalize(out);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpz_set(mpq_numref(out), e->data.bigint);
        mpz_set_ui(mpq_denref(out), 1);
        return true;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2) {
        struct Expr* a = e->data.function.args[0];
        struct Expr* b = e->data.function.args[1];
        mpq_t qa, qb;
        mpq_init(qa); mpq_init(qb);
        bool ok = expr_to_mpq(a, qa) && expr_to_mpq(b, qb)
                  && mpz_sgn(mpq_numref(qb)) != 0;
        if (ok) mpq_div(out, qa, qb);
        mpq_clear(qa); mpq_clear(qb);
        return ok;
    }
    return false;
}

/* Try to interpret a single Expr as one (coef, exps) term against vars. */
static bool expr_term(struct Expr* e, struct Expr** vars, int n_vars,
                      mpq_t coef_out, int* exps_out) {
    mpq_set_ui(coef_out, 1, 1);
    for (int i = 0; i < n_vars; i++) exps_out[i] = 0;

    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        return expr_to_mpq(e, coef_out);
    }
    if (e->type == EXPR_SYMBOL) {
        int vi = find_var_index(e, vars, n_vars);
        if (vi < 0) return false;
        exps_out[vi] = 1;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return false;

    /* Rational literal at the term level. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational) {
        return expr_to_mpq(e, coef_out);
    }

    /* Power[var, k] with k >= 0. */
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

    /* Times[...] -- multiply contributions. */
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times) {
        mpq_t sub_coef; mpq_init(sub_coef);
        int* sub_exps = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
        bool ok = true;
        for (size_t a = 0; a < e->data.function.arg_count && ok; a++) {
            if (!expr_term(e->data.function.args[a], vars, n_vars,
                           sub_coef, sub_exps)) {
                ok = false;
                break;
            }
            mpq_mul(coef_out, coef_out, sub_coef);
            for (int v = 0; v < n_vars; v++) exps_out[v] += sub_exps[v];
        }
        mpq_clear(sub_coef);
        free(sub_exps);
        return ok;
    }

    /* Last-ditch: maybe it's a known variable Expr that's a function-
     * shaped symbol (rare).  Already covered for plain symbols above. */
    int vi = find_var_index(e, vars, n_vars);
    if (vi >= 0) { exps_out[vi] = 1; return true; }
    return false;
}

GBPoly* gb_from_expr(struct Expr* e, struct Expr** vars, int n_vars,
                     GBOrder order, int elim_pivot,
                     const GBWeightMatrix* wmat) {
    if (!e) return NULL;
    GBPoly* p = gb_poly_new(n_vars, order, elim_pivot);
    if (order == GB_ORDER_MATRIX) p->wmat = wmat;

    if (e->type == EXPR_FUNCTION &&
        e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Plus) {
        mpq_t coef; mpq_init(coef);
        int* exps = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
        for (size_t a = 0; a < e->data.function.arg_count; a++) {
            if (!expr_term(e->data.function.args[a], vars, n_vars, coef, exps)) {
                mpq_clear(coef);
                free(exps);
                gb_poly_free(p);
                return NULL;
            }
            gb_poly_push_term(p, exps, coef);
        }
        mpq_clear(coef);
        free(exps);
        gb_poly_normalize(p);
        return p;
    }

    mpq_t coef; mpq_init(coef);
    int* exps = (int*)malloc(sizeof(int) * (size_t)(n_vars > 0 ? n_vars : 1));
    if (!expr_term(e, vars, n_vars, coef, exps)) {
        mpq_clear(coef);
        free(exps);
        gb_poly_free(p);
        return NULL;
    }
    gb_poly_push_term(p, exps, coef);
    mpq_clear(coef);
    free(exps);
    gb_poly_normalize(p);
    return p;
}

/* Emit a single coefficient as an Expr: Integer / BigInt / Rational. */
static struct Expr* mpq_to_expr(const mpq_t q) {
    /* If denominator is 1, return Integer or BigInt directly. */
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) {
        if (mpz_fits_slong_p(mpq_numref(q))) {
            return expr_new_integer((int64_t)mpz_get_si(mpq_numref(q)));
        }
        return expr_new_bigint_from_mpz(mpq_numref(q));
    }
    /* Rational[num, den]. */
    struct Expr* num = mpz_fits_slong_p(mpq_numref(q))
        ? expr_new_integer((int64_t)mpz_get_si(mpq_numref(q)))
        : expr_new_bigint_from_mpz(mpq_numref(q));
    struct Expr* den = mpz_fits_slong_p(mpq_denref(q))
        ? expr_new_integer((int64_t)mpz_get_si(mpq_denref(q)))
        : expr_new_bigint_from_mpz(mpq_denref(q));
    struct Expr** args = (struct Expr**)malloc(sizeof(struct Expr*) * 2);
    args[0] = num; args[1] = den;
    struct Expr* r = expr_new_function(expr_new_symbol("Rational"), args, 2);
    free(args);
    return r;
}

struct Expr* gb_to_expr(const GBPoly* p, struct Expr** vars) {
    if (p->n_terms == 0) return expr_new_integer(0);

    struct Expr** terms = (struct Expr**)malloc(sizeof(struct Expr*) * p->n_terms);
    for (size_t i = 0; i < p->n_terms; i++) {
        const int* row = gb_exp_at((GBPoly*)p, i);
        size_t nfac = 0;
        for (int v = 0; v < p->n_vars; v++) if (row[v] > 0) nfac++;

        bool coef_is_one = (mpz_cmp_ui(mpq_numref(p->coefs[i]), 1) == 0
                            && mpz_cmp_ui(mpq_denref(p->coefs[i]), 1) == 0);
        bool coef_is_neg_one = (mpz_cmp_si(mpq_numref(p->coefs[i]), -1) == 0
                                && mpz_cmp_ui(mpq_denref(p->coefs[i]), 1) == 0);

        size_t total = nfac;
        if (!(coef_is_one || coef_is_neg_one) || nfac == 0) total++;
        if (coef_is_neg_one && nfac > 0) total++;

        struct Expr** factors = (struct Expr**)malloc(
            sizeof(struct Expr*) * (total > 0 ? total : 1));
        size_t fi = 0;

        if (coef_is_neg_one && nfac > 0) {
            factors[fi++] = expr_new_integer(-1);
        } else if (!coef_is_one || nfac == 0) {
            factors[fi++] = mpq_to_expr(p->coefs[i]);
        }
        for (int v = 0; v < p->n_vars; v++) {
            if (row[v] == 0) continue;
            if (row[v] == 1) {
                factors[fi++] = expr_copy(vars[v]);
            } else {
                struct Expr** pa = (struct Expr**)malloc(sizeof(struct Expr*) * 2);
                pa[0] = expr_copy(vars[v]);
                pa[1] = expr_new_integer(row[v]);
                factors[fi++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
                free(pa);
            }
        }
        if (fi == 1) {
            terms[i] = factors[0];
        } else {
            terms[i] = expr_new_function(expr_new_symbol("Times"), factors, fi);
        }
        free(factors);
    }

    struct Expr* out;
    if (p->n_terms == 1) {
        out = terms[0];
    } else {
        out = expr_new_function(expr_new_symbol("Plus"), terms, p->n_terms);
    }
    free(terms);
    return out;
}

/* Like gb_to_expr but for the InexactNumbers coefficient domain: the
 * polynomial is made monic (leading coefficient 1) and every coefficient
 * -- including a unit leading coefficient -- is emitted as a `bits`-bit
 * MPFR real, matching Mathematica's approximate-arithmetic output
 * (`1.00000000000000000 a y^4`, `0.666... x`).  The caller owns the
 * returned Expr*. */
struct Expr* gb_to_expr_inexact(const GBPoly* p, struct Expr** vars,
                                mpfr_prec_t bits) {
    if (p->n_terms == 0) return expr_new_mpfr_bits(bits);   /* 0.0 */

    GBPoly* m = gb_poly_copy(p);
    gb_poly_make_monic(m);

    struct Expr** terms = (struct Expr**)malloc(sizeof(struct Expr*) * m->n_terms);
    for (size_t i = 0; i < m->n_terms; i++) {
        const int* row = gb_exp_at(m, i);
        size_t nfac = 0;
        for (int v = 0; v < m->n_vars; v++) if (row[v] > 0) nfac++;

        /* Coefficient is always emitted, as a real factor. */
        mpfr_t c; mpfr_init2(c, bits);
        mpfr_set_q(c, m->coefs[i], MPFR_RNDN);
        struct Expr* coef = expr_new_mpfr_move(c);   /* takes ownership of c */

        size_t total = nfac + 1;
        struct Expr** factors = (struct Expr**)malloc(sizeof(struct Expr*) * total);
        size_t fi = 0;
        factors[fi++] = coef;
        for (int v = 0; v < m->n_vars; v++) {
            if (row[v] == 0) continue;
            if (row[v] == 1) {
                factors[fi++] = expr_copy(vars[v]);
            } else {
                struct Expr** pa = (struct Expr**)malloc(sizeof(struct Expr*) * 2);
                pa[0] = expr_copy(vars[v]);
                pa[1] = expr_new_integer(row[v]);
                factors[fi++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
                free(pa);
            }
        }
        if (fi == 1) {
            terms[i] = factors[0];
        } else {
            terms[i] = expr_new_function(expr_new_symbol("Times"), factors, fi);
        }
        free(factors);
    }

    struct Expr* out;
    if (m->n_terms == 1) {
        out = terms[0];
    } else {
        out = expr_new_function(expr_new_symbol("Plus"), terms, m->n_terms);
    }
    free(terms);
    gb_poly_free(m);
    return out;
}
