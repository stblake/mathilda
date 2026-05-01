/*
 * bpoly.c
 * -------
 * Implementation of bivariate polynomials over Z; see bpoly.h for the
 * type contract.
 *
 * Implementation choices:
 *
 * 1. Storage: dense in x (the main variable), sparse-in-NULL for cx[i]
 *    that happens to be the zero polynomial.  Operations always check
 *    `cx[i] != NULL && !zupoly_is_zero(cx[i])` before doing real work.
 *
 * 2. Capacity growth follows the same pattern as ZUPoly: realloc to
 *    double-or-fit, fill new slots with NULL.
 *
 * 3. Bivariate exact division is implemented as long-division viewing
 *    P as Z[y][x] -- at each step the leading x-coefficient of the
 *    running remainder is divided (via zupoly_divexact) by the
 *    leading x-coefficient of the divisor.  This is exactly the same
 *    structure as ZUPoly division but with ZUPoly arithmetic at each
 *    coefficient.
 *
 * 4. Conversion to / from Expr uses the existing poly machinery: we
 *    extract the coefficient of each x^i (which itself is a polynomial
 *    in y), then convert that y-polynomial via expr_to_zupoly.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bpoly.h"
#include "expr.h"
#include "expand.h"
#include "eval.h"
#include "poly.h"

/* ====================================================================== */
/*  Internal helpers                                                      */
/* ====================================================================== */

static void bpoly_reserve(BPoly* p, int new_cap) {
    if (new_cap <= p->cap_x) return;
    int target = (p->cap_x > 0) ? p->cap_x : 1;
    while (target < new_cap) target *= 2;
    p->cx = (ZUPoly**)realloc(p->cx, sizeof(ZUPoly*) * (size_t)target);
    for (int i = p->cap_x; i < target; i++) p->cx[i] = NULL;
    p->cap_x = target;
}

/* True if cx[i] represents zero (NULL or all zero coefficients). */
static inline bool slot_is_zero(const BPoly* p, int i) {
    return !p->cx[i] || zupoly_is_zero(p->cx[i]);
}

/* ====================================================================== */
/*  Construction / destruction                                            */
/* ====================================================================== */

BPoly* bpoly_new(int cap) {
    if (cap < 1) cap = 1;
    BPoly* p = (BPoly*)malloc(sizeof(BPoly));
    p->deg_x = -1;
    p->cap_x = cap;
    p->cx = (ZUPoly**)calloc((size_t)cap, sizeof(ZUPoly*));
    return p;
}

BPoly* bpoly_zero(void) { return bpoly_new(1); }

BPoly* bpoly_copy(const BPoly* p) {
    int cap = (p->deg_x < 0) ? 1 : p->deg_x + 1;
    BPoly* q = bpoly_new(cap);
    for (int i = 0; i <= p->deg_x; i++) {
        if (p->cx[i] && !zupoly_is_zero(p->cx[i])) {
            q->cx[i] = zupoly_copy(p->cx[i]);
        }
    }
    q->deg_x = p->deg_x;
    return q;
}

void bpoly_free(BPoly* p) {
    if (!p) return;
    for (int i = 0; i < p->cap_x; i++) {
        if (p->cx[i]) zupoly_free(p->cx[i]);
    }
    free(p->cx);
    free(p);
}

/* ====================================================================== */
/*  Coefficient access                                                    */
/* ====================================================================== */

void bpoly_set_xcoef(BPoly* p, int i, ZUPoly* c) {
    if (i < 0) { if (c) zupoly_free(c); return; }
    bpoly_reserve(p, i + 1);
    if (p->cx[i]) zupoly_free(p->cx[i]);
    p->cx[i] = (c && !zupoly_is_zero(c)) ? c : NULL;
    if (c && zupoly_is_zero(c)) zupoly_free(c);
    /* Update deg_x. */
    if (p->cx[i]) {
        if (i > p->deg_x) p->deg_x = i;
    } else if (i == p->deg_x) {
        while (p->deg_x >= 0 && slot_is_zero(p, p->deg_x)) p->deg_x--;
    }
}

const ZUPoly* bpoly_get_xcoef(const BPoly* p, int i) {
    if (i < 0 || i > p->deg_x) return NULL;
    return p->cx[i];
}

void bpoly_normalize(BPoly* p) {
    int d = p->cap_x - 1;
    while (d >= 0 && slot_is_zero(p, d)) d--;
    p->deg_x = d;
}

/* ====================================================================== */
/*  Predicates and degree                                                 */
/* ====================================================================== */

bool bpoly_is_zero(const BPoly* p) { return p->deg_x < 0; }

bool bpoly_eq(const BPoly* a, const BPoly* b) {
    if (a->deg_x != b->deg_x) return false;
    for (int i = 0; i <= a->deg_x; i++) {
        bool az = slot_is_zero(a, i), bz = slot_is_zero(b, i);
        if (az != bz) return false;
        if (!az && !zupoly_eq(a->cx[i], b->cx[i])) return false;
    }
    return true;
}

int bpoly_deg_x(const BPoly* p) { return p->deg_x; }

int bpoly_deg_y(const BPoly* p) {
    int d = -1;
    for (int i = 0; i <= p->deg_x; i++) {
        if (p->cx[i] && p->cx[i]->deg > d) d = p->cx[i]->deg;
    }
    return d;
}

const ZUPoly* bpoly_lc_x(const BPoly* p) {
    if (p->deg_x < 0) return NULL;
    return p->cx[p->deg_x];
}

/* ====================================================================== */
/*  Arithmetic                                                            */
/* ====================================================================== */

BPoly* bpoly_add(const BPoly* a, const BPoly* b) {
    int max_d = (a->deg_x > b->deg_x) ? a->deg_x : b->deg_x;
    if (max_d < 0) return bpoly_zero();
    BPoly* r = bpoly_new(max_d + 1);
    for (int i = 0; i <= max_d; i++) {
        const ZUPoly* ai = (i <= a->deg_x) ? a->cx[i] : NULL;
        const ZUPoly* bi = (i <= b->deg_x) ? b->cx[i] : NULL;
        ZUPoly* sum;
        if (!ai && !bi) sum = NULL;
        else if (!ai) sum = zupoly_copy(bi);
        else if (!bi) sum = zupoly_copy(ai);
        else          sum = zupoly_add(ai, bi);
        if (sum && zupoly_is_zero(sum)) {
            zupoly_free(sum); sum = NULL;
        }
        r->cx[i] = sum;
    }
    bpoly_normalize(r);
    return r;
}

BPoly* bpoly_sub(const BPoly* a, const BPoly* b) {
    int max_d = (a->deg_x > b->deg_x) ? a->deg_x : b->deg_x;
    if (max_d < 0) return bpoly_zero();
    BPoly* r = bpoly_new(max_d + 1);
    for (int i = 0; i <= max_d; i++) {
        const ZUPoly* ai = (i <= a->deg_x) ? a->cx[i] : NULL;
        const ZUPoly* bi = (i <= b->deg_x) ? b->cx[i] : NULL;
        ZUPoly* diff;
        if (!ai && !bi)      diff = NULL;
        else if (!ai)        diff = zupoly_neg(bi);
        else if (!bi)        diff = zupoly_copy(ai);
        else                 diff = zupoly_sub(ai, bi);
        if (diff && zupoly_is_zero(diff)) {
            zupoly_free(diff); diff = NULL;
        }
        r->cx[i] = diff;
    }
    bpoly_normalize(r);
    return r;
}

BPoly* bpoly_mul(const BPoly* a, const BPoly* b) {
    if (a->deg_x < 0 || b->deg_x < 0) return bpoly_zero();
    int rd = a->deg_x + b->deg_x;
    BPoly* r = bpoly_new(rd + 1);
    for (int i = 0; i <= a->deg_x; i++) {
        if (slot_is_zero(a, i)) continue;
        for (int j = 0; j <= b->deg_x; j++) {
            if (slot_is_zero(b, j)) continue;
            ZUPoly* term = zupoly_mul(a->cx[i], b->cx[j]);
            if (zupoly_is_zero(term)) { zupoly_free(term); continue; }
            int idx = i + j;
            if (!r->cx[idx]) r->cx[idx] = term;
            else {
                ZUPoly* sum = zupoly_add(r->cx[idx], term);
                zupoly_free(r->cx[idx]); zupoly_free(term);
                r->cx[idx] = zupoly_is_zero(sum) ? (zupoly_free(sum), NULL) : sum;
            }
        }
    }
    bpoly_normalize(r);
    return r;
}

BPoly* bpoly_neg(const BPoly* a) {
    BPoly* r = bpoly_copy(a);
    for (int i = 0; i <= r->deg_x; i++) {
        if (r->cx[i]) {
            ZUPoly* n = zupoly_neg(r->cx[i]);
            zupoly_free(r->cx[i]);
            r->cx[i] = n;
        }
    }
    return r;
}

BPoly* bpoly_mul_truncate_y(const BPoly* a, const BPoly* b,
                            int max_y_plus_one) {
    BPoly* full = bpoly_mul(a, b);
    BPoly* tr   = bpoly_truncate_y(full, max_y_plus_one);
    bpoly_free(full);
    return tr;
}

BPoly* bpoly_mul_zupoly(const BPoly* a, const ZUPoly* s) {
    if (a->deg_x < 0 || zupoly_is_zero(s)) return bpoly_zero();
    BPoly* r = bpoly_new(a->deg_x + 1);
    for (int i = 0; i <= a->deg_x; i++) {
        if (slot_is_zero(a, i)) continue;
        ZUPoly* prod = zupoly_mul(a->cx[i], s);
        if (zupoly_is_zero(prod)) { zupoly_free(prod); continue; }
        r->cx[i] = prod;
    }
    bpoly_normalize(r);
    return r;
}

/* Bivariate exact division viewing P as Z[y][x].
 *
 * Algorithm: textbook polynomial long division.  Each elimination step
 * computes c = leading-x-coeff(remainder) / leading-x-coeff(divisor),
 * a ZUPoly division which must succeed exactly over Z[y].  If any
 * intermediate ZUPoly division is non-exact, the bivariate division
 * fails (returns NULL). */
BPoly* bpoly_divexact(const BPoly* a, const BPoly* b) {
    if (b->deg_x < 0) return NULL;
    if (a->deg_x < 0) return bpoly_zero();
    if (a->deg_x < b->deg_x) return NULL;

    int q_deg = a->deg_x - b->deg_x;
    BPoly* q = bpoly_new(q_deg + 1);
    BPoly* r = bpoly_copy(a);

    for (int i = q_deg; i >= 0; i--) {
        int rd = b->deg_x + i;
        if (rd > r->deg_x || slot_is_zero(r, rd)) continue;

        /* Quotient coefficient at x^i: cx_r[rd] / cx_b[deg_b]. */
        ZUPoly* coef = zupoly_divexact(r->cx[rd], b->cx[b->deg_x]);
        if (!coef) {
            bpoly_free(q); bpoly_free(r);
            return NULL;
        }
        /* Place coef in q at position i. */
        bpoly_reserve(q, i + 1);
        if (q->cx[i]) zupoly_free(q->cx[i]);
        q->cx[i] = coef;
        if (i > q->deg_x) q->deg_x = i;

        /* Subtract coef * b * x^i from r. */
        for (int j = 0; j <= b->deg_x; j++) {
            if (slot_is_zero(b, j)) continue;
            ZUPoly* term = zupoly_mul(coef, b->cx[j]);
            if (zupoly_is_zero(term)) { zupoly_free(term); continue; }
            int idx = i + j;
            ZUPoly* old = r->cx[idx];
            ZUPoly* new_;
            if (!old) {
                new_ = zupoly_neg(term);
            } else {
                new_ = zupoly_sub(old, term);
                zupoly_free(old);
            }
            zupoly_free(term);
            if (zupoly_is_zero(new_)) {
                zupoly_free(new_);
                r->cx[idx] = NULL;
            } else {
                r->cx[idx] = new_;
            }
        }
        bpoly_normalize(r);
    }

    if (!bpoly_is_zero(r)) {
        bpoly_free(q); bpoly_free(r);
        return NULL;
    }
    bpoly_free(r);
    bpoly_normalize(q);
    return q;
}

/* ====================================================================== */
/*  Truncation, evaluation, shift                                         */
/* ====================================================================== */

BPoly* bpoly_truncate_y(const BPoly* p, int k) {
    if (k <= 0) return bpoly_zero();
    BPoly* r = bpoly_new(p->deg_x < 0 ? 1 : p->deg_x + 1);
    for (int i = 0; i <= p->deg_x; i++) {
        if (slot_is_zero(p, i)) continue;
        const ZUPoly* src = p->cx[i];
        /* Truncate src modulo y^k by copying low coefficients only. */
        ZUPoly* tr = zupoly_new(k);
        for (int j = 0; j < k && j <= src->deg; j++) {
            const mpz_t* cj = zupoly_getcoef(src, j);
            if (cj) zupoly_setcoef(tr, j, *cj);
        }
        zupoly_normalize(tr);
        if (zupoly_is_zero(tr)) { zupoly_free(tr); continue; }
        r->cx[i] = tr;
    }
    bpoly_normalize(r);
    return r;
}

ZUPoly* bpoly_eval_y_si(const BPoly* p, int64_t alpha) {
    if (p->deg_x < 0) return zupoly_zero();
    ZUPoly* result = zupoly_new(p->deg_x + 1);
    mpz_t v; mpz_init(v);
    for (int i = 0; i <= p->deg_x; i++) {
        if (slot_is_zero(p, i)) continue;
        zupoly_eval_si(p->cx[i], alpha, v);
        if (mpz_sgn(v) != 0) zupoly_setcoef(result, i, v);
    }
    mpz_clear(v);
    zupoly_normalize(result);
    return result;
}

BPoly* bpoly_shift_y_si(const BPoly* p, int64_t alpha) {
    if (p->deg_x < 0) return bpoly_zero();
    BPoly* r = bpoly_new(p->deg_x + 1);
    for (int i = 0; i <= p->deg_x; i++) {
        if (slot_is_zero(p, i)) continue;
        ZUPoly* shifted = zupoly_shift_si(p->cx[i], alpha);
        if (zupoly_is_zero(shifted)) { zupoly_free(shifted); continue; }
        r->cx[i] = shifted;
    }
    bpoly_normalize(r);
    return r;
}

/* ====================================================================== */
/*  Conversion to / from picocas Expr                                     */
/* ====================================================================== */

BPoly* expr_to_bpoly(const struct Expr* e,
                     const struct Expr* x_var,
                     const struct Expr* y_var) {
    Expr* expanded = expr_expand((Expr*)e);
    int deg_x = get_degree_poly(expanded, (Expr*)x_var);
    if (deg_x < 0) {
        /* The expression is the zero polynomial -- still a valid
         * (trivial) bivariate. */
        expr_free(expanded);
        return bpoly_zero();
    }
    BPoly* p = bpoly_new(deg_x + 1);
    for (int i = 0; i <= deg_x; i++) {
        Expr* ci = get_coeff(expanded, (Expr*)x_var, i);
        /* ci is a polynomial in y (and possibly other variables that
         * we will reject below).  Convert via expr_to_zupoly. */
        ZUPoly* yi = expr_to_zupoly(ci, y_var);
        expr_free(ci);
        if (!yi) {
            bpoly_free(p);
            expr_free(expanded);
            return NULL;
        }
        if (!zupoly_is_zero(yi)) {
            p->cx[i] = yi;
        } else {
            zupoly_free(yi);
        }
    }
    bpoly_normalize(p);
    expr_free(expanded);
    return p;
}

struct Expr* bpoly_to_expr(const BPoly* p,
                           const struct Expr* x_var,
                           const struct Expr* y_var) {
    if (bpoly_is_zero(p)) return expr_new_integer(0);

    /* Build sum of terms: each x^i * (y-poly) contributes a sum over
     * x^i * c_{i,j} * y^j.  We expand by converting each cx[i] to an
     * Expr in y and multiplying by x^i.
     *
     * The multiplication by x^i is done by parsing /constructing an
     * Expr; the evaluator flattens / canonicalises as needed. */
    int n_terms = 0;
    for (int i = 0; i <= p->deg_x; i++) {
        if (!slot_is_zero(p, i)) n_terms++;
    }
    if (n_terms == 0) return expr_new_integer(0);

    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)n_terms);
    int idx = 0;
    for (int i = 0; i <= p->deg_x; i++) {
        if (slot_is_zero(p, i)) continue;
        Expr* y_part = zupoly_to_expr(p->cx[i], y_var);
        Expr* term;
        if (i == 0) {
            term = y_part;
        } else {
            Expr* x_pow;
            if (i == 1) x_pow = expr_copy((Expr*)x_var);
            else {
                Expr* args[2] = { expr_copy((Expr*)x_var), expr_new_integer(i) };
                x_pow = expr_new_function(expr_new_symbol("Power"), args, 2);
                x_pow = eval_and_free(x_pow);
            }
            Expr* args[2] = { y_part, x_pow };
            term = expr_new_function(expr_new_symbol("Times"), args, 2);
            term = eval_and_free(term);
        }
        terms[idx++] = term;
    }
    Expr* result;
    if (n_terms == 1) {
        result = terms[0];
    } else {
        result = expr_new_function(expr_new_symbol("Plus"), terms, n_terms);
        result = eval_and_free(result);
    }
    free(terms);
    return result;
}

/* ====================================================================== */
/*  Debug                                                                 */
/* ====================================================================== */

void bpoly_print(const BPoly* p, const char* x_name, const char* y_name) {
    if (bpoly_is_zero(p)) { fprintf(stderr, "0"); return; }
    bool first = true;
    for (int i = 0; i <= p->deg_x; i++) {
        if (slot_is_zero(p, i)) continue;
        if (!first) fprintf(stderr, " + ");
        first = false;
        if (i == 0) {
            zupoly_print(p->cx[i], y_name);
        } else {
            fprintf(stderr, "(");
            zupoly_print(p->cx[i], y_name);
            fprintf(stderr, ")");
            if (i == 1) fprintf(stderr, "*%s", x_name);
            else        fprintf(stderr, "*%s^%d", x_name, i);
        }
    }
}
