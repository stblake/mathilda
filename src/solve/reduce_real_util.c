/*
 * reduce_real_util.c
 *
 * Shared real-algebraic primitives for Reduce's Reals engines.  See
 * reduce_real_util.h.  Extracted verbatim from the original reduce_univar.c so
 * both the univariate sign diagram and the multivariate CAD share one exact
 * sign/root/sample implementation.
 */
#include "reduce_real_util.h"

#include "eval.h"
#include "sym_names.h"
#include "flint_qqbar.h"

#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ *
 *  Exact real sign oracle                                             *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Sign of a numeric literal; *is_num=false if `e` is not a plain number. */
static int number_sign(const Expr* e, bool* is_num) {
    *is_num = true;
    switch (e->type) {
        case EXPR_INTEGER: return (e->data.integer > 0) - (e->data.integer < 0);
        case EXPR_REAL:    return (e->data.real > 0.0) - (e->data.real < 0.0);
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint);
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_sgn(e->data.mpfr);
#endif
        case EXPR_FUNCTION:
            if (is_head(e, SYM_Rational) && e->data.function.arg_count == 2) {
                const Expr* n = e->data.function.args[0];
                if (n->type == EXPR_INTEGER) return (n->data.integer > 0) - (n->data.integer < 0);
                if (n->type == EXPR_BIGINT)  return mpz_sgn(n->data.bigint);
            }
            *is_num = false; return 0;
        default: *is_num = false; return 0;
    }
}

int rru_sign_of(const Expr* e) {
    bool isn; int s = number_sign(e, &isn);
    if (isn) return s;
    Expr* zero = expr_new_integer(0);
    int r = flint_qqbar_compare(e, zero);   /* sign(e - 0) */
    expr_free(zero);
    return r;
}

int rru_sign_compare(const Expr* a, const Expr* b) {
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
    int s = rru_sign_of(diff);
    expr_free(diff);
    return s;
}

/* ------------------------------------------------------------------ *
 *  Sample-point selection                                             *
 * ------------------------------------------------------------------ */

double rru_approx_double(const Expr* e, bool* ok) {
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy((Expr*)e) }, 1));
    double d = 0.0; *ok = true;
    if (v->type == EXPR_REAL)         d = v->data.real;
    else if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
    else if (v->type == EXPR_BIGINT)  d = mpz_get_d(v->data.bigint);
#ifdef USE_MPFR
    else if (v->type == EXPR_MPFR)    d = mpfr_get_d(v->data.mpfr, MPFR_RNDN);
#endif
    else *ok = false;
    expr_free(v);
    return d;
}

/* num/den as a canonical Integer or Rational. */
static Expr* make_rational(long num, long den) {
    if (den == 1) return expr_new_integer(num);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(num),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_integer(den), expr_new_integer(-1) }, 2) }, 2));
}

Expr* rru_rational_between(const Expr* lo, const Expr* hi) {
    bool ok;
    if (!lo && !hi) return expr_new_integer(0);

    if (!lo) {                              /* (-inf, hi) */
        double b = rru_approx_double(hi, &ok);
        if (!ok) return NULL;
        long s = (long)floor(b) - 1;
        for (int i = 0; i < 64; i++, s--) {
            Expr* c = expr_new_integer(s);
            int cmp = rru_sign_compare(c, hi);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp < 0) return c;
            expr_free(c);
        }
        return NULL;
    }
    if (!hi) {                              /* (lo, +inf) */
        double a = rru_approx_double(lo, &ok);
        if (!ok) return NULL;
        long s = (long)ceil(a) + 1;
        for (int i = 0; i < 64; i++, s++) {
            Expr* c = expr_new_integer(s);
            int cmp = rru_sign_compare(c, lo);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp > 0) return c;
            expr_free(c);
        }
        return NULL;
    }

    double a = rru_approx_double(lo, &ok); if (!ok) return NULL;
    double b = rru_approx_double(hi, &ok); if (!ok) return NULL;
    double mid = 0.5 * (a + b);
    for (long D = 1; D <= (1L << 30); D *= 2) {
        long num = (long)llround(mid * (double)D);
        Expr* c = make_rational(num, D);
        int c1 = rru_sign_compare(c, lo);       /* want > 0 */
        int c2 = rru_sign_compare(c, hi);       /* want < 0 */
        if (c1 == -2 || c2 == -2) { expr_free(c); return NULL; }
        if (c1 > 0 && c2 < 0) return c;
        expr_free(c);
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Polynomial test and point sign                                     *
 * ------------------------------------------------------------------ */

bool rru_is_polynomial(const Expr* poly, const Expr* x) {
    Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)x) }, 2));
    bool r = (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True);
    expr_free(q);
    return r;
}

int rru_poly_sign_at(const Expr* poly, const Expr* x, const Expr* sample) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)sample) }, 2);
    Expr* val = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rule }, 2));
    int s = rru_sign_of(val);
    expr_free(val);
    return s;
}

/* ------------------------------------------------------------------ *
 *  Root collection (with optional provenance)                         *
 * ------------------------------------------------------------------ */

bool rru_collect_roots(const Expr* poly, const Expr* x,
                       Expr*** arr, int* n, int* cap,
                       int** prov, int factor_id) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ eqn, expr_copy((Expr*)x), expr_new_symbol(SYM_Reals) }, 3);
    Expr* sols = eval_and_free(call);

    bool ok = is_head(sols, SYM_List);
    for (size_t i = 0; ok && i < sols->data.function.arg_count; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) { ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) { ok = false; break; }
        Expr* val = rule->data.function.args[1];
        if (is_head(val, SYM_ConditionalExpression)) { ok = false; break; }
        if (*n == *cap) {
            *cap = *cap ? *cap * 2 : 8;
            *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*));
            if (prov) *prov = realloc(*prov, (size_t)*cap * sizeof(int));
        }
        if (prov) (*prov)[*n] = factor_id;
        (*arr)[(*n)++] = expr_copy(val);
    }
    expr_free(sols);
    return ok;
}
