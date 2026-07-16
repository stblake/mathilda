/* ============================================================================
 * series.c - Series and SeriesData
 * ============================================================================
 *
 * This module implements the power-series machinery for Mathilda.
 *
 * SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] is the data head
 * that represents a truncated power series. The i-th coefficient multiplies
 * (x - x0)^((nmin + i)/den) and an O[x - x0]^(nmax/den) term captures the
 * dropped higher-order terms.
 *
 * Series[f, {x, x0, n}]  expands f as a power series in (x - x0) up to
 * order n. Series also accepts the leading-term form
 * Series[f, x -> x0] and the iterated multivariate form
 * Series[f, {x, x0, nx}, {y, y0, ny}, ...]. The algorithm is a recursive
 * "series algebra": primitive subexpressions become SeriesObj's, algebraic
 * heads (Plus, Times, Power) combine them, and elementary heads
 * (Exp, Log, Sin, Cos, Sinh, Cosh, Tan, Tanh) apply their known series
 * kernels. Unknown heads fall back to naive Taylor via D[...]. Expansion
 * about Infinity is handled by substituting x -> 1/u internally and
 * presenting the result with Power[x, -1] as the series variable.
 *
 * Normal[s] drops the O-term from a SeriesData and returns an ordinary sum.
 */

#include "series.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "arithmetic.h"
#include "poly.h"
#include "rationalize.h"
#include "sym_names.h"
#include "ndarray.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

/* ----------------------------------------------------------------------------
 * Tiny Expr helpers
 * -------------------------------------------------------------------------- */

static int64_t abs_i64(int64_t v) { return v < 0 ? -v : v; }

static Expr* mk_symbol(const char* s) { return expr_new_symbol(s); }

/* Build a fresh function expression from owned pieces. The returned node
 * takes ownership of head and of each arg; the caller passes ownership in. */
static Expr* mk_fn(Expr* head, Expr** args, size_t n) {
    return expr_new_function(head, args, n);
}

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return mk_fn(mk_symbol(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return mk_fn(mk_symbol(name), args, 2);
}

static Expr* mk_plus(Expr* a, Expr* b)  { return mk_fn2("Plus",  a, b); }
static Expr* mk_times(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* mk_power(Expr* a, Expr* b) { return mk_fn2("Power", a, b); }

/* Simplify via the evaluator; takes ownership, returns owned.
 * evaluate() does NOT consume its argument, so free `e` after. */
static Expr* simp(Expr* e) { return eval_and_free(e); }

static bool is_lit_zero(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    return false;
}

static bool has_symbol_head(Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Cheap structural test for SeriesData[x, x0, {coefs}, nmin, nmax, den]. Uses
 * the interned SYM_SeriesData pointer so the arithmetic builtins can scan their
 * arguments without a strcmp on the hot path. */
bool is_series_data(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_SeriesData &&
           e->data.function.arg_count == 6;
}

/* Structural containment: true iff target does not appear anywhere in e. */
static bool expr_free_of(Expr* e, Expr* target) {
    if (!e) return true;
    if (expr_eq(e, target)) return false;
    if (e->type == EXPR_FUNCTION) {
        if (!expr_free_of(e->data.function.head, target)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!expr_free_of(e->data.function.args[i], target)) return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------------- *
 * split_two_term: decompose e = a + b*x^(p/q)
 *
 * A Maxima-style probe (cf. split-two-term-poly in series.lisp) that looks
 * at the structural shape of `e` without doing a full series expansion.
 * Returns true when `e` can be written as a sum of a term free of `x` and
 * a single monomial term in `x` with a rational exponent. This is used as
 * a cheap shape check by several fast paths:
 *
 *  - Log[a + b x^(p/q)] can be rewritten as Log[a] + Log[1 + (b/a) x^(p/q)]
 *    and handed to the Log1p kernel without expanding the arg.
 *  - (a + b x^(p/q))^alpha maps onto the monomial binomial fast path.
 *  - Apart preprocessing skips inputs that are already in "a + b x^c" form.
 *
 * On success the caller receives ownership of *a_out and *b_out.
 *
 * Handled forms:
 *   x                     -> (0, 1, 1/1)
 *   free-of-x expression  -> (e, 0, 1/1)
 *   Plus[t_i]             -> sum a-parts; at most one exponent class in b
 *   Times[f_i]            -> at most one non-free factor; others fold into prod
 *   Power[x, p/q]         -> (0, 1, p/q) when the exponent is a rational
 *
 * Returns false (and sets *a_out, *b_out to NULL) for anything else, e.g.
 * (1 + x)(1 + x) or Sin[x] -- those would need full series algebra.
 * ------------------------------------------------------------------------- */
bool series_split_two_term(Expr* e, Expr* x,
                           Expr** a_out, Expr** b_out,
                           int64_t* exp_num, int64_t* exp_den) {
    *a_out = NULL; *b_out = NULL;
    *exp_num = 1; *exp_den = 1;

    if (expr_eq(e, x)) {
        *a_out = expr_new_integer(0);
        *b_out = expr_new_integer(1);
        *exp_num = 1; *exp_den = 1;
        return true;
    }
    /* "Free of x" is tested *after* the e==x check so that nested
     * lookups into the same symbol still bind correctly. */
    if (expr_free_of(e, x)) {
        *a_out = expr_copy(e);
        *b_out = expr_new_integer(0);
        return true;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) {
        return false;
    }
    const char* head = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;

    if (head == SYM_Plus) {
        Expr* a_sum = expr_new_integer(0);
        Expr* b_sum = NULL;
        int64_t cn = 0, cd = 1;
        bool have_b = false;
        for (size_t i = 0; i < n; i++) {
            Expr* aa; Expr* bb; int64_t ccn, ccd;
            if (!series_split_two_term(e->data.function.args[i], x, &aa, &bb, &ccn, &ccd)) {
                expr_free(a_sum);
                if (b_sum) expr_free(b_sum);
                return false;
            }
            a_sum = simp(mk_plus(a_sum, aa));
            if (is_lit_zero(bb)) { expr_free(bb); continue; }
            if (!have_b) {
                b_sum = bb;
                cn = ccn; cd = ccd;
                have_b = true;
            } else {
                /* Require the same rational exponent. Compare p/q == p'/q' as
                 * p*q' == p'*q -- both fit int64 at the call sites we see. */
                if (ccn * cd != cn * ccd) {
                    expr_free(bb); expr_free(a_sum); expr_free(b_sum);
                    return false;
                }
                b_sum = simp(mk_plus(b_sum, bb));
            }
        }
        *a_out = a_sum;
        *b_out = b_sum ? b_sum : expr_new_integer(0);
        *exp_num = have_b ? cn : 1;
        *exp_den = have_b ? cd : 1;
        return true;
    }

    if (head == SYM_Times) {
        Expr* prod = expr_new_integer(1);
        Expr* a_inner = NULL;
        Expr* b_inner = NULL;
        int64_t cn = 0, cd = 1;
        bool have_nf = false;
        for (size_t i = 0; i < n; i++) {
            Expr* aa; Expr* bb; int64_t ccn, ccd;
            if (!series_split_two_term(e->data.function.args[i], x, &aa, &bb, &ccn, &ccd)) {
                expr_free(prod);
                if (a_inner) expr_free(a_inner);
                if (b_inner) expr_free(b_inner);
                return false;
            }
            if (is_lit_zero(bb)) {
                /* Factor has no b part -- multiply into prod. */
                expr_free(bb);
                prod = simp(mk_times(prod, aa));
                continue;
            }
            if (have_nf) {
                /* Two non-free factors would give cross terms we can't represent. */
                expr_free(aa); expr_free(bb); expr_free(prod);
                expr_free(a_inner); expr_free(b_inner);
                return false;
            }
            a_inner = aa; b_inner = bb;
            cn = ccn; cd = ccd;
            have_nf = true;
        }
        if (have_nf) {
            /* (a_inner + b_inner x^c) * prod = prod*a_inner + prod*b_inner * x^c. */
            *a_out = simp(mk_times(expr_copy(prod), a_inner));
            *b_out = simp(mk_times(prod, b_inner));
            *exp_num = cn; *exp_den = cd;
            return true;
        }
        *a_out = prod;
        *b_out = expr_new_integer(0);
        return true;
    }

    if (head == SYM_Power && n == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (expr_eq(base, x) && expr_free_of(exp, x)) {
            int64_t en, ed;
            if (exp->type == EXPR_INTEGER) { en = exp->data.integer; ed = 1; }
            else if (is_rational(exp, &en, &ed)) { /* ok */ }
            else return false;
            if (en == 0) return false;  /* x^0 would have been simplified away */
            *a_out = expr_new_integer(0);
            *b_out = expr_new_integer(1);
            *exp_num = en; *exp_den = ed;
            return true;
        }
        return false;
    }

    return false;
}

/* ReplaceAll helper: substitute `from` with `to` everywhere in `e`.
 * Returns a fresh expression (caller owns); does not free inputs. */
static Expr* replace_all_of(Expr* e, Expr* from, Expr* to) {
    Expr* rule_args[2] = { expr_copy(from), expr_copy(to) };
    Expr* rule = mk_fn(mk_symbol("Rule"), rule_args, 2);
    Expr* ra_args[2] = { expr_copy(e), rule };
    Expr* ra = mk_fn(mk_symbol("ReplaceAll"), ra_args, 2);
    return simp(ra);
}

/* ----------------------------------------------------------------------------
 * SeriesObj -- internal representation during computation
 *
 * A series is sum_{i=0..coef_count-1} coefs[i] * (x - x0)^((nmin+i)/den),
 * with an implied O[x - x0]^(order/den) term. Coefficients are ordinary
 * Expr* nodes (may contain symbolic expressions such as Log[x] or
 * Derivative[k][f][x0]) and are kept in canonical (evaluated) form.
 * -------------------------------------------------------------------------- */

typedef struct {
    Expr*   x;          /* expansion variable, owned */
    Expr*   x0;         /* expansion point, owned */
    Expr**  coefs;      /* owned coefficients, length == coef_count */
    size_t  coef_count;
    int64_t nmin;       /* leading exponent numerator */
    int64_t order;      /* O-term exponent numerator (exclusive) */
    int64_t den;        /* common denominator for exponents (>= 1) */
} SeriesObj;

/* Context threaded through the recursive series_expand. Defined here so the
 * branch-point handlers (which sit above series_expand in the file) can
 * inspect/mutate it directly. See full field documentation at the (former)
 * declaration site just above series_expand. */
typedef struct {
    Expr*   x;
    Expr*   x0;
    int64_t order;
    int64_t target_order;
    bool    allow_branch_wrap;
    Expr*   pending_add_const;
    Expr*   pending_log_coef;
    Expr*   pending_discriminator;
} SeriesCtx;

static SeriesObj* series_expand(Expr* e, SeriesCtx* ctx);
static SeriesObj* so_from_seriesdata(Expr* sd);

static SeriesObj* so_alloc(Expr* x, Expr* x0, int64_t nmin, int64_t order, int64_t den) {
    SeriesObj* s = calloc(1, sizeof(SeriesObj));
    s->x = expr_copy(x);
    s->x0 = expr_copy(x0);
    s->nmin = nmin;
    s->order = order < nmin ? nmin : order;
    s->den = den <= 0 ? 1 : den;
    int64_t count = s->order - s->nmin;
    if (count < 0) count = 0;
    s->coef_count = (size_t)count;
    if (count > 0) {
        s->coefs = calloc((size_t)count, sizeof(Expr*));
        for (int64_t i = 0; i < count; i++) s->coefs[i] = expr_new_integer(0);
    }
    return s;
}

static void so_free(SeriesObj* s) {
    if (!s) return;
    if (s->x) expr_free(s->x);
    if (s->x0) expr_free(s->x0);
    if (s->coefs) {
        for (size_t i = 0; i < s->coef_count; i++) {
            if (s->coefs[i]) expr_free(s->coefs[i]);
        }
        free(s->coefs);
    }
    free(s);
}

/* Replace coefficient at index i, freeing the previous value. */
static void so_set_coef(SeriesObj* s, size_t i, Expr* v) {
    if (i >= s->coef_count) { expr_free(v); return; }
    expr_free(s->coefs[i]);
    s->coefs[i] = v;
}

/* Drop leading zeros (advance nmin) and trailing zeros never happen to
 * affect `order` in the series algebra, so we only trim leading zeros. */
static void so_trim_leading(SeriesObj* s) {
    size_t drop = 0;
    while (drop < s->coef_count && is_lit_zero(s->coefs[drop])) drop++;
    if (drop == 0) return;
    if (drop == s->coef_count) {
        for (size_t i = 0; i < s->coef_count; i++) expr_free(s->coefs[i]);
        free(s->coefs);
        s->coefs = NULL;
        s->coef_count = 0;
        s->nmin = s->order; /* empty series, purely O-term */
        return;
    }
    size_t new_count = s->coef_count - drop;
    Expr** new_coefs = calloc(new_count, sizeof(Expr*));
    for (size_t i = 0; i < drop; i++) expr_free(s->coefs[i]);
    for (size_t i = 0; i < new_count; i++) new_coefs[i] = s->coefs[drop + i];
    free(s->coefs);
    s->coefs = new_coefs;
    s->coef_count = new_count;
    s->nmin += (int64_t)drop;
}

/* Rescale to a new (larger) denominator. The current denominator must
 * divide the new one. Exponents expand by factor k = new_den/den. */
static SeriesObj* so_rescale(SeriesObj* s, int64_t new_den) {
    if (new_den == s->den) {
        /* Return a copy so the caller can always free the result. */
        SeriesObj* c = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) so_set_coef(c, i, expr_copy(s->coefs[i]));
        return c;
    }
    int64_t k = new_den / s->den;
    int64_t new_nmin = s->nmin * k;
    int64_t new_order = s->order * k;
    SeriesObj* r = so_alloc(s->x, s->x0, new_nmin, new_order, new_den);
    /* Only multiples-of-k exponents carry coefficients. */
    for (size_t i = 0; i < s->coef_count; i++) {
        so_set_coef(r, (size_t)(i * k), expr_copy(s->coefs[i]));
    }
    return r;
}

static int64_t gcd_i64(int64_t a, int64_t b) {
    a = abs_i64(a); b = abs_i64(b);
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a ? a : 1;
}

static int64_t lcm_i64(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    return (a / gcd_i64(a, b)) * b;
}

/* Align two series to a common denominator. Returns fresh objects; the
 * caller owns both. */
static void so_align_den(SeriesObj* a, SeriesObj* b, SeriesObj** ra, SeriesObj** rb) {
    int64_t d = lcm_i64(a->den, b->den);
    *ra = so_rescale(a, d);
    *rb = so_rescale(b, d);
}

/* Constant series c + O[(x-x0)^order_over_den]. */
static SeriesObj* so_from_constant(Expr* c, Expr* x, Expr* x0, int64_t order, int64_t den) {
    SeriesObj* s = so_alloc(x, x0, 0, order, den);
    if (s->coef_count > 0) {
        so_set_coef(s, 0, expr_copy(c));
    }
    return s;
}

/* Identity series: x expressed around x0, i.e. x = x0 + 1 * (x - x0). */
static SeriesObj* so_from_variable(Expr* x, Expr* x0, int64_t order, int64_t den) {
    SeriesObj* s = so_alloc(x, x0, 0, order, den);
    if (s->coef_count > 0) so_set_coef(s, 0, expr_copy(x0));
    if (s->coef_count > 1) so_set_coef(s, 1, expr_new_integer(1));
    return s;
}

/* ----------------------------------------------------------------------------
 * Conversion between SeriesObj and SeriesData Expr trees
 * -------------------------------------------------------------------------- */

static Expr* so_to_expr(const SeriesObj* s) {
    Expr** list_args = NULL;
    if (s->coef_count > 0) {
        list_args = calloc(s->coef_count, sizeof(Expr*));
        for (size_t i = 0; i < s->coef_count; i++) list_args[i] = expr_copy(s->coefs[i]);
    }
    Expr* coefs_list = expr_new_function(mk_symbol("List"), list_args, s->coef_count);
    if (list_args) free(list_args);
    Expr* args[6] = {
        expr_copy(s->x),
        expr_copy(s->x0),
        coefs_list,
        expr_new_integer(s->nmin),
        expr_new_integer(s->order),
        expr_new_integer(s->den)
    };
    return mk_fn(mk_symbol("SeriesData"), args, 6);
}

/* ----------------------------------------------------------------------------
 * Basic series arithmetic
 * -------------------------------------------------------------------------- */

static SeriesObj* so_add(SeriesObj* a, SeriesObj* b) {
    SeriesObj *aa, *bb;
    so_align_den(a, b, &aa, &bb);
    int64_t new_nmin  = aa->nmin  < bb->nmin  ? aa->nmin  : bb->nmin;
    int64_t new_order = aa->order < bb->order ? aa->order : bb->order;
    SeriesObj* r = so_alloc(aa->x, aa->x0, new_nmin, new_order, aa->den);
    for (int64_t k = new_nmin; k < new_order; k++) {
        Expr* sum = expr_new_integer(0);
        int64_t ia = k - aa->nmin;
        int64_t ib = k - bb->nmin;
        if (ia >= 0 && (size_t)ia < aa->coef_count) sum = simp(mk_plus(sum, expr_copy(aa->coefs[ia])));
        if (ib >= 0 && (size_t)ib < bb->coef_count) sum = simp(mk_plus(sum, expr_copy(bb->coefs[ib])));
        so_set_coef(r, (size_t)(k - new_nmin), sum);
    }
    so_free(aa); so_free(bb);
    return r;
}

static SeriesObj* so_scalar_mul(Expr* c, SeriesObj* a) {
    SeriesObj* r = so_alloc(a->x, a->x0, a->nmin, a->order, a->den);
    for (size_t i = 0; i < a->coef_count; i++) {
        so_set_coef(r, i, simp(mk_times(expr_copy(c), expr_copy(a->coefs[i]))));
    }
    return r;
}

static SeriesObj* so_neg(SeriesObj* a) {
    Expr* m1 = expr_new_integer(-1);
    SeriesObj* r = so_scalar_mul(m1, a);
    expr_free(m1);
    return r;
}

static SeriesObj* so_sub(SeriesObj* a, SeriesObj* b) {
    SeriesObj* nb = so_neg(b);
    SeriesObj* r = so_add(a, nb);
    so_free(nb);
    return r;
}

static SeriesObj* so_mul(SeriesObj* a, SeriesObj* b) {
    SeriesObj *aa, *bb;
    so_align_den(a, b, &aa, &bb);
    int64_t result_nmin = aa->nmin + bb->nmin;
    /* result_order is min(a.order + b.nmin, b.order + a.nmin): the lowest
     * exponent at which we lose precision from either O-term contribution. */
    int64_t oa = aa->order + bb->nmin;
    int64_t ob = bb->order + aa->nmin;
    int64_t result_order = oa < ob ? oa : ob;
    SeriesObj* r = so_alloc(aa->x, aa->x0, result_nmin, result_order, aa->den);
    for (size_t i = 0; i < aa->coef_count; i++) {
        if (is_lit_zero(aa->coefs[i])) continue;
        for (size_t j = 0; j < bb->coef_count; j++) {
            if (is_lit_zero(bb->coefs[j])) continue;
            int64_t k = aa->nmin + (int64_t)i + bb->nmin + (int64_t)j;
            if (k < result_nmin || k >= result_order) continue;
            size_t idx = (size_t)(k - result_nmin);
            Expr* prod = simp(mk_times(expr_copy(aa->coefs[i]), expr_copy(bb->coefs[j])));
            Expr* sum = simp(mk_plus(expr_copy(r->coefs[idx]), prod));
            so_set_coef(r, idx, sum);
        }
    }
    so_free(aa); so_free(bb);
    return r;
}

/* Multiply a series by (x - x0)^(p/q). In the rescaled (common-den) form
 * this is just a shift in nmin by p*den/q -- which requires q to divide
 * den (or we rescale). Returns a fresh series. */
static SeriesObj* so_shift_by_rational(SeriesObj* s, int64_t p, int64_t q) {
    int64_t g = gcd_i64(abs_i64(p), abs_i64(q));
    if (g) { p /= g; q /= g; }
    if (q < 0) { q = -q; p = -p; }
    int64_t needed_den = lcm_i64(s->den, q);
    SeriesObj* a = so_rescale(s, needed_den);
    int64_t shift = p * (needed_den / q);
    SeriesObj* r = so_alloc(a->x, a->x0, a->nmin + shift, a->order + shift, a->den);
    for (size_t i = 0; i < a->coef_count; i++) {
        so_set_coef(r, i, expr_copy(a->coefs[i]));
    }
    so_free(a);
    return r;
}

/* Make a copy of a with leading zeros trimmed. */
static SeriesObj* so_copy_trimmed(SeriesObj* a) {
    SeriesObj* c = so_alloc(a->x, a->x0, a->nmin, a->order, a->den);
    for (size_t i = 0; i < a->coef_count; i++) so_set_coef(c, i, expr_copy(a->coefs[i]));
    so_trim_leading(c);
    return c;
}

/* Count non-zero coefficients (cheap check for "single-term" exactness). */
static size_t so_nonzero_count(SeriesObj* s) {
    size_t c = 0;
    for (size_t i = 0; i < s->coef_count; i++) if (!is_lit_zero(s->coefs[i])) c++;
    return c;
}

/* Series inverse: 1/a. Requires the leading coefficient to be non-zero
 * after trimming. Handles the exact single-term case without the usual
 * order reduction; otherwise the correction expansion loses 2*nmin of
 * precision against the O-term exponent. */
/* Local copy of a leaf-count walk so so_inv can estimate input size
 * without pulling in limit.c's helper. */
static int64_t so_leaf_count(Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    int64_t c = so_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        c += so_leaf_count(e->data.function.args[i]);
    }
    return c;
}

/* A coefficient is "purely numeric" if it is built from integers, bigints,
 * reals, and the heads Rational/Complex/Times/Plus/Power with numeric
 * operands. Purely-numeric coefficients don't trigger the symbolic-blowup
 * problem the guardrail protects against, so so_inv can compute the full N
 * terms for them regardless of leaf count. */
static bool so_is_purely_numeric(Expr* e) {
    if (!e) return true;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_SYMBOL: case EXPR_STRING: case EXPR_NDARRAY: return false;
        case EXPR_FUNCTION: {
            Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL) return false;
            const char* s = h->data.symbol.name;
            if (s != SYM_Rational && s != SYM_Complex &&
                s != SYM_Times && s != SYM_Plus &&
                s != SYM_Power) return false;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!so_is_purely_numeric(e->data.function.args[i])) return false;
            }
            return true;
        }
    }
    return false;
}

static SeriesObj* so_inv(SeriesObj* a_in) {
    SeriesObj* a = so_copy_trimmed(a_in);
    if (a->coef_count == 0) { so_free(a); return NULL; }

    /* Exact single-term inversion: 1/(c * (x-x0)^nmin) = (1/c) * (x-x0)^{-nmin}
     * and its accuracy is whatever the input's was. */
    if (so_nonzero_count(a) == 1) {
        Expr* inv_c = simp(mk_power(expr_copy(a->coefs[0]), expr_new_integer(-1)));
        SeriesObj* r = so_alloc(a->x, a->x0, -a->nmin, a->order, a->den);
        if (r->coef_count > 0) so_set_coef(r, 0, inv_c);
        else expr_free(inv_c);
        so_free(a);
        return r;
    }

    int64_t ord_rel = a->order - a->nmin;
    size_t N = (size_t)ord_rel;
    if (N == 0) {
        SeriesObj* r = so_alloc(a->x, a->x0, -a->nmin, -a->nmin, a->den);
        so_free(a);
        return r;
    }
    /* Expression-growth guardrail for series inversion with purely symbolic
     * coefficients (e.g. polynomials in Log[2], Log[3] arising from
     * 1/(2^x - 3^x)). Each iteration of the b[k] recurrence multiplies
     * growing expressions, and since Mathilda's evaluator does not fully
     * canonicalise polynomials in independent symbolic constants, the
     * compute time per iteration can grow super-linearly. To keep Series
     * responsive, we cap N at a limit calibrated from the total input
     * coefficient size: small-coefficient inputs (numeric or short
     * symbolic) get the full N, but once the input coefficients are large
     * symbolic expressions we stop early and return a truncated series.
     * The result is still a valid leading-term Laurent expansion; callers
     * that only need the leading behaviour (Limit, 1/f fast paths) are
     * unaffected, and callers that need higher-order expansions will see
     * a smaller-than-requested O-term but not a hang. */
    /* Count leaves in NON-ZERO coefficients only: the guardrail exists to
     * bound symbolic blowup during the b[k] recurrence (each nonzero a[i]
     * feeds into O(N) symbolic multiplies). Zero coefficients contribute
     * nothing to that work but inflate a naive leaf sum for polynomials
     * with mostly-zero coefficients like 1/(1+x^N). Purely-numeric coefs
     * also don't grow under simp() so they don't count against the cap. */
    int64_t total_leaves = 0;
    for (size_t i = 0; i < a->coef_count && i < N; i++) {
        if (is_lit_zero(a->coefs[i])) continue;
        if (so_is_purely_numeric(a->coefs[i])) continue;
        total_leaves += so_leaf_count(a->coefs[i]);
    }
    size_t N_cap;
    if (total_leaves <= 4)        N_cap = 64;   /* numeric / trivial coefs */
    else if (total_leaves <= 20)  N_cap = 16;
    else if (total_leaves <= 60)  N_cap = 6;
    else                           N_cap = 4;
    if (N > N_cap) N = N_cap;
    Expr** b = calloc(N, sizeof(Expr*));
    Expr* inv_a0 = simp(mk_power(expr_copy(a->coefs[0]), expr_new_integer(-1)));
    b[0] = expr_copy(inv_a0);
    /* When every non-zero coefficient of `a` is purely numeric, the b[k]
     * recurrence stays numeric too, and we can accumulate with per-term simp()
     * without fearing O(k^2) symbolic blowup. This both avoids a pathology in
     * the evaluator when it sees a very large flat-Plus tree of numeric
     * products, and lets so_inv compute the full N b[]-values even for high-
     * order Laurent inversions (e.g. 1/Sin[x] at order 9+).
     *
     * For non-numeric inputs we fall back to the batched path, which exists
     * to prevent super-linear evaluation time when coefficients are
     * polynomials in independent symbolic constants like Log[2], Log[3]. */
    bool coefs_all_numeric = true;
    for (size_t i = 0; i < a->coef_count; i++) {
        if (is_lit_zero(a->coefs[i])) continue;
        if (!so_is_purely_numeric(a->coefs[i])) { coefs_all_numeric = false; break; }
    }

    for (size_t k = 1; k < N; k++) {
        if (coefs_all_numeric) {
            Expr* sum = expr_new_integer(0);
            for (size_t i = 1; i <= k; i++) {
                if (i >= a->coef_count) break;
                if (is_lit_zero(a->coefs[i])) continue;
                if (!b[k - i] || is_lit_zero(b[k - i])) continue;
                Expr* term = simp(mk_times(expr_copy(a->coefs[i]), expr_copy(b[k - i])));
                sum = simp(mk_plus(sum, term));
            }
            b[k] = simp(mk_times(expr_new_integer(-1),
                                 mk_times(expr_copy(inv_a0), sum)));
            if (!b[k]) b[k] = expr_new_integer(0);
            continue;
        }
        /* Batched path for symbolic coefficients: build a single flat Plus
         * of all non-zero product contributions and simp() it exactly once.
         * Reasons:
         *   1. Each term simp() call evaluates the partial sum plus a new
         *      product; because Mathilda's evaluator doesn't fully canonicalize
         *      polynomials in symbolic constants (e.g. Log[2], Log[3]), the
         *      inner-loop simp() can blow up to O(k^2) or worse in time.
         *   2. A single end-of-loop simp() lets Plus/Times flattening and
         *      Orderless sorting collapse common sub-expressions in one pass.
         * This matters for shapes like 1/(2^x - 3^x) where each a_i is a
         * polynomial in Log[2], Log[3]; the old loop took seconds per k
         * and never completed past a few iterations. */
        size_t pc = 0;
        Expr** parts = calloc(k, sizeof(Expr*));
        for (size_t i = 1; i <= k; i++) {
            if (i >= a->coef_count) break;
            if (is_lit_zero(a->coefs[i])) continue;
            if (!b[k - i] || is_lit_zero(b[k - i])) continue;
            parts[pc++] = mk_times(expr_copy(a->coefs[i]), expr_copy(b[k - i]));
        }
        Expr* sum;
        if (pc == 0)      sum = expr_new_integer(0);
        else if (pc == 1) sum = parts[0];
        else              sum = expr_new_function(mk_symbol("Plus"), parts, pc);
        free(parts);
        Expr* neg_over_a0 = mk_times(expr_new_integer(-1),
                                     mk_times(expr_copy(inv_a0), sum));
        b[k] = simp(neg_over_a0);
        if (!b[k]) b[k] = expr_new_integer(0);
    }
    expr_free(inv_a0);

    SeriesObj* r = so_alloc(a->x, a->x0, -a->nmin, a->order - 2 * a->nmin, a->den);
    for (size_t i = 0; i < N && i < r->coef_count; i++) so_set_coef(r, i, b[i]);
    free(b);
    so_free(a);
    return r;
}

static SeriesObj* so_div(SeriesObj* a, SeriesObj* b) {
    SeriesObj* binv = so_inv(b);
    if (!binv) return NULL;
    SeriesObj* r = so_mul(a, binv);
    so_free(binv);
    return r;
}

static SeriesObj* so_pow_int(SeriesObj* a, int64_t n) {
    if (n == 0) {
        SeriesObj* one = so_alloc(a->x, a->x0, 0, a->order - a->nmin, a->den);
        if (one->coef_count > 0) so_set_coef(one, 0, expr_new_integer(1));
        return one;
    }
    if (n < 0) {
        SeriesObj* inv = so_inv(a);
        if (!inv) return NULL;
        SeriesObj* p = so_pow_int(inv, -n);
        so_free(inv);
        return p;
    }
    /* Binary exponentiation. */
    SeriesObj* result = NULL;
    SeriesObj* base = so_alloc(a->x, a->x0, a->nmin, a->order, a->den);
    for (size_t i = 0; i < a->coef_count; i++) so_set_coef(base, i, expr_copy(a->coefs[i]));
    int64_t e = n;
    while (e > 0) {
        if (e & 1) {
            if (!result) {
                result = so_alloc(base->x, base->x0, base->nmin, base->order, base->den);
                for (size_t i = 0; i < base->coef_count; i++) so_set_coef(result, i, expr_copy(base->coefs[i]));
            } else {
                SeriesObj* nr = so_mul(result, base);
                so_free(result);
                result = nr;
            }
        }
        e >>= 1;
        if (e > 0) {
            SeriesObj* sq = so_mul(base, base);
            so_free(base);
            base = sq;
        }
    }
    so_free(base);
    return result;
}

/* Monomial fast path for (1 + c*x^e)^alpha.
 *
 * When u has exactly one non-zero coefficient (a monomial), the generic
 * Horner-composition step in so_pow_1plus_alpha spends O(N^2) work on
 * so_mul convolutions that are mostly multiplying by zero. Here we emit
 * coefficients directly: the k-th binomial coefficient sits at exponent
 * k * e, and everything in between is zero.
 *
 *   (1 + c*x^e)^alpha = sum_k C(alpha, k) * c^k * x^(e*k)
 *
 * with C(alpha, k) generated by the recurrence
 *   C(alpha, 0) = 1,  C(alpha, k+1) = C(alpha, k) * (alpha - k) / (k + 1).
 *
 * Returns NULL if u has more than one non-zero coefficient (caller should
 * fall back to the generic path). Preserves u.den and u.order. */
static SeriesObj* so_pow_1plus_alpha_monomial(SeriesObj* u, Expr* alpha) {
    size_t monomial_idx = SIZE_MAX;
    for (size_t i = 0; i < u->coef_count; i++) {
        if (is_lit_zero(u->coefs[i])) continue;
        if (monomial_idx != SIZE_MAX) return NULL;
        monomial_idx = i;
    }

    SeriesObj* r = so_alloc(u->x, u->x0, 0, u->order, u->den);
    if (r->coef_count == 0) return r;

    /* (1 + 0)^alpha = 1 when u is the zero series. */
    if (monomial_idx == SIZE_MAX) {
        so_set_coef(r, 0, expr_new_integer(1));
        return r;
    }

    int64_t e = u->nmin + (int64_t)monomial_idx;   /* monomial exponent numerator */
    if (e <= 0) { so_free(r); return NULL; }       /* u.nmin >= 1 invariant */

    Expr* c = expr_copy(u->coefs[monomial_idx]);
    Expr* b_k = expr_new_integer(1);               /* C(alpha, 0) */
    Expr* c_k = expr_new_integer(1);               /* c^0 */

    for (int64_t k = 0; ; k++) {
        int64_t exp_num = k * e;
        if (exp_num >= u->order) break;
        size_t idx = (size_t)exp_num;
        if (idx >= r->coef_count) break;
        Expr* coef_val = simp(mk_times(expr_copy(b_k), expr_copy(c_k)));
        so_set_coef(r, idx, coef_val);
        /* Advance b_{k+1} = b_k * (alpha - k) / (k + 1). */
        Expr* diff = simp(mk_plus(expr_copy(alpha), expr_new_integer(-k)));
        Expr* numer = simp(mk_times(b_k, diff));
        Expr* rat   = make_rational(1, k + 1);
        b_k = simp(mk_times(numer, rat));
        /* Advance c_{k+1} = c_k * c. */
        Expr* c_next = simp(mk_times(c_k, expr_copy(c)));
        c_k = c_next;
    }
    expr_free(b_k);
    expr_free(c_k);
    expr_free(c);
    return r;
}

/* Compute (1 + u)^alpha, given u with u.nmin >= 1 (u has no constant term),
 * using the recurrence b_{k+1} = (alpha - k)/(k+1) * b_k applied
 * symbolically. alpha can be an arbitrary Expr (rational, symbolic).
 * Returns a series whose exponents are multiples of u.nmin up to
 * u.order (in u.den), i.e. in the original common denominator. */
static SeriesObj* so_pow_1plus_alpha(SeriesObj* u, Expr* alpha) {
    /* Fast path: u is a monomial (0 or 1 non-zero coefficient). */
    {
        SeriesObj* fast = so_pow_1plus_alpha_monomial(u, alpha);
        if (fast) return fast;
    }

    /* Determine how many terms of the 1+u kernel we need: u.nmin is the
     * minimum power contributed, and we need kernel-index k such that
     * k * u.nmin < u.order. */
    int64_t u_nmin = u->nmin;
    int64_t u_order = u->order;
    if (u_nmin <= 0) {
        /* u must have no constant term. Caller's responsibility. */
        return NULL;
    }
    int64_t N = 0;
    while ((N + 1) * u_nmin < u_order) N++;
    /* Build binomial coefficients b_0..b_N symbolically. */
    Expr** bcoef = calloc((size_t)(N + 1), sizeof(Expr*));
    bcoef[0] = expr_new_integer(1);
    for (int64_t k = 0; k < N; k++) {
        /* b_{k+1} = b_k * (alpha - k) / (k+1) */
        Expr* diff = simp(mk_plus(expr_copy(alpha), expr_new_integer(-k)));
        Expr* num = simp(mk_times(expr_copy(bcoef[k]), diff));
        Expr* inv = expr_new_integer(k + 1);
        Expr* inv_rat = simp(mk_power(inv, expr_new_integer(-1)));
        bcoef[k + 1] = simp(mk_times(num, inv_rat));
    }
    /* Compose: sum_k b_k * u^k, via Horner from high-degree to low. */
    SeriesObj* result = so_alloc(u->x, u->x0, 0, u->order, u->den);
    if (result->coef_count > 0) so_set_coef(result, 0, expr_copy(bcoef[N]));
    for (int64_t k = N - 1; k >= 0; k--) {
        SeriesObj* mul = so_mul(result, u);
        so_free(result);
        result = mul;
        /* Add b_k * (x-x0)^0 constant. */
        SeriesObj* c = so_from_constant(bcoef[k], u->x, u->x0, u->order, u->den);
        SeriesObj* sum = so_add(result, c);
        so_free(result); so_free(c);
        result = sum;
    }
    for (int64_t i = 0; i <= N; i++) expr_free(bcoef[i]);
    free(bcoef);
    return result;
}

/* s^alpha for alpha an arbitrary Expr. Factor out the leading monomial:
 *   s = a0 * (x-x0)^(nmin/den) * (1 + u)
 * with u = s/(a0*(x-x0)^(nmin/den)) - 1 having u.nmin >= 1.
 * Then s^alpha = a0^alpha * (x-x0)^(nmin*alpha/den) * (1+u)^alpha. */
static SeriesObj* so_pow_expr(SeriesObj* s_in, Expr* alpha) {
    if (alpha->type == EXPR_INTEGER) {
        return so_pow_int(s_in, alpha->data.integer);
    }

    SeriesObj* s = so_copy_trimmed(s_in);
    if (s->coef_count == 0) { so_free(s); return NULL; }
    Expr* a0 = s->coefs[0];

    /* Build u = s / (a0 * (x-x0)^(nmin/den)) - 1, a series with nmin >= 1. */
    /* Step 1: scalar-divide by a0 (so leading coef becomes 1). */
    Expr* inv_a0 = simp(mk_power(expr_copy(a0), expr_new_integer(-1)));
    SeriesObj* s_scaled = so_scalar_mul(inv_a0, s);
    expr_free(inv_a0);
    /* Step 2: shift so leading exponent is 0 (i.e. subtract nmin from every
     * exponent, then drop the leading 1 and keep u = the rest). */
    int64_t nmin = s->nmin;
    SeriesObj* shifted = so_shift_by_rational(s_scaled, -nmin, s->den);
    so_free(s_scaled);
    /* Step 3: u = shifted - 1 (remove the constant term, should be 1). */
    Expr* one = expr_new_integer(1);
    SeriesObj* cone = so_from_constant(one, shifted->x, shifted->x0, shifted->order, shifted->den);
    SeriesObj* u = so_sub(shifted, cone);
    so_free(cone); so_free(shifted); expr_free(one);
    so_trim_leading(u);
    if (u->coef_count == 0) {
        /* Exact single-term input: s = a0*(x-x0)^(nmin/den), so
         * s^alpha = a0^alpha * (x-x0)^(nmin*alpha/den). If alpha is rational,
         * we can express the shift exactly; otherwise we give up. */
        so_free(u);
        Expr* a0_pow = simp(mk_power(expr_copy(a0), expr_copy(alpha)));
        if (nmin == 0) {
            SeriesObj* r = so_alloc(s->x, s->x0, 0, s->order, s->den);
            if (r->coef_count > 0) so_set_coef(r, 0, a0_pow);
            else expr_free(a0_pow);
            so_free(s);
            return r;
        }
        int64_t p, q;
        if (alpha->type == EXPR_INTEGER) { p = alpha->data.integer; q = 1; }
        else if (is_rational(alpha, &p, &q)) { /* ok */ }
        else { expr_free(a0_pow); so_free(s); return NULL; }
        /* Build SeriesData with the single coefficient a0_pow at exponent
         * nmin*p/(den*q). */
        int64_t out_p = nmin * p;
        int64_t out_q = s->den * q;
        int64_t g = gcd_i64(abs_i64(out_p), abs_i64(out_q));
        if (g) { out_p /= g; out_q /= g; }
        if (out_q < 0) { out_q = -out_q; out_p = -out_p; }
        SeriesObj* r = so_alloc(s->x, s->x0, out_p, s->order, out_q);
        if (r->coef_count > 0) so_set_coef(r, 0, a0_pow);
        else expr_free(a0_pow);
        so_free(s);
        return r;
    }
    /* Step 4: (1+u)^alpha. */
    SeriesObj* oneplus_pow = so_pow_1plus_alpha(u, alpha);
    so_free(u);
    if (!oneplus_pow) { so_free(s); return NULL; }

    /* Step 5: scalar multiply by a0^alpha and shift by nmin*alpha/den. */
    Expr* a0_pow_alpha = simp(mk_power(expr_copy(a0), expr_copy(alpha)));
    SeriesObj* scaled = so_scalar_mul(a0_pow_alpha, oneplus_pow);
    expr_free(a0_pow_alpha);
    so_free(oneplus_pow);

    SeriesObj* final_result = scaled;
    if (nmin != 0) {
        int64_t p, q;
        if (alpha->type == EXPR_INTEGER) { p = alpha->data.integer; q = 1; }
        else if (is_rational(alpha, &p, &q)) { /* ok */ }
        else { so_free(scaled); so_free(s); return NULL; }
        final_result = so_shift_by_rational(scaled, nmin * p, s->den * q);
        so_free(scaled);
    }
    so_free(s);
    return final_result;
}

/* ----------------------------------------------------------------------------
 * Elementary kernels and function application
 * -------------------------------------------------------------------------- */

/* Compose sum_k kernel[k] * u^k, where kernel[] is a length-N array of
 * Expr* scalar coefficients and u is a series with u.nmin >= 1. */
static SeriesObj* so_compose_scalar_kernel(Expr** kernel, size_t N, SeriesObj* u) {
    if (N == 0) {
        return so_from_constant(expr_new_integer(0), u->x, u->x0, u->order, u->den);
    }
    SeriesObj* result = so_from_constant(kernel[N - 1], u->x, u->x0, u->order, u->den);
    /* Iterate k = N-2, N-3, ..., 0 using a size_t-safe reverse loop. */
    for (size_t k = N - 1; k-- > 0; ) {
        SeriesObj* mul = so_mul(result, u);
        so_free(result);
        result = mul;
        SeriesObj* c = so_from_constant(kernel[k], u->x, u->x0, u->order, u->den);
        SeriesObj* sum = so_add(result, c);
        so_free(result); so_free(c);
        result = sum;
    }
    return result;
}

/* Build the Taylor series coefficients of an elementary function f(u) at
 * u = 0, as scalar Expr* values, for u-powers 0..N-1.
 *   Exp:   1, 1, 1/2, 1/6, 1/24, ...          (1/k!)
 *   Log1p: 0, 1, -1/2, 1/3, -1/4, ...         ((-1)^(k-1)/k for k>=1)
 *   Sin:   0, 1, 0, -1/6, 0, 1/120, ...
 *   Cos:   1, 0, -1/2, 0, 1/24, ...
 *   Sinh:  0, 1, 0,  1/6, 0,  1/120, ...
 *   Cosh:  1, 0,  1/2, 0,  1/24, ...
 * Returns a newly allocated array of owned Expr*. */
static Expr** kernel_coefs(const char* name, size_t N) {
    Expr** c = calloc(N, sizeof(Expr*));
    if (strcmp(name, "Exp") == 0) {
        /* c[k] = 1/k! */
        Expr* fact = expr_new_integer(1);
        for (size_t k = 0; k < N; k++) {
            if (k > 0) fact = simp(mk_times(fact, expr_new_integer((int64_t)k)));
            c[k] = simp(mk_power(expr_copy(fact), expr_new_integer(-1)));
        }
        expr_free(fact);
    } else if (strcmp(name, "Log1p") == 0) {
        /* Log[1+u] = sum_{k>=1} (-1)^(k-1) u^k / k */
        c[0] = expr_new_integer(0);
        for (size_t k = 1; k < N; k++) {
            int64_t sign = (k & 1) ? 1 : -1;
            c[k] = make_rational(sign, (int64_t)k);
        }
    } else if (strcmp(name, "Sin") == 0) {
        Expr* fact = expr_new_integer(1);
        for (size_t k = 0; k < N; k++) {
            if (k > 0) fact = simp(mk_times(fact, expr_new_integer((int64_t)k)));
            if (k % 2 == 0) c[k] = expr_new_integer(0);
            else {
                int64_t sign = ((k / 2) & 1) ? -1 : 1;
                Expr* inv = simp(mk_power(expr_copy(fact), expr_new_integer(-1)));
                c[k] = simp(mk_times(expr_new_integer(sign), inv));
            }
        }
        expr_free(fact);
    } else if (strcmp(name, "Cos") == 0) {
        Expr* fact = expr_new_integer(1);
        for (size_t k = 0; k < N; k++) {
            if (k > 0) fact = simp(mk_times(fact, expr_new_integer((int64_t)k)));
            if (k % 2 != 0) c[k] = expr_new_integer(0);
            else {
                int64_t sign = ((k / 2) & 1) ? -1 : 1;
                Expr* inv = simp(mk_power(expr_copy(fact), expr_new_integer(-1)));
                c[k] = simp(mk_times(expr_new_integer(sign), inv));
            }
        }
        expr_free(fact);
    } else if (strcmp(name, "Sinh") == 0) {
        Expr* fact = expr_new_integer(1);
        for (size_t k = 0; k < N; k++) {
            if (k > 0) fact = simp(mk_times(fact, expr_new_integer((int64_t)k)));
            if (k % 2 == 0) c[k] = expr_new_integer(0);
            else c[k] = simp(mk_power(expr_copy(fact), expr_new_integer(-1)));
        }
        expr_free(fact);
    } else if (strcmp(name, "Cosh") == 0) {
        Expr* fact = expr_new_integer(1);
        for (size_t k = 0; k < N; k++) {
            if (k > 0) fact = simp(mk_times(fact, expr_new_integer((int64_t)k)));
            if (k % 2 != 0) c[k] = expr_new_integer(0);
            else c[k] = simp(mk_power(expr_copy(fact), expr_new_integer(-1)));
        }
        expr_free(fact);
    } else if (strcmp(name, "ArcTan") == 0) {
        /* ArcTan[u] = sum_{k odd, k>=1} (-1)^((k-1)/2) u^k / k */
        for (size_t k = 0; k < N; k++) {
            if ((k & 1) == 0) { c[k] = expr_new_integer(0); continue; }
            int64_t sign = (((k - 1) / 2) & 1) ? -1 : 1;
            c[k] = make_rational(sign, (int64_t)k);
        }
    } else if (strcmp(name, "ArcTanh") == 0) {
        /* ArcTanh[u] = sum_{k odd, k>=1} u^k / k */
        for (size_t k = 0; k < N; k++) {
            if ((k & 1) == 0) c[k] = expr_new_integer(0);
            else c[k] = make_rational(1, (int64_t)k);
        }
    } else if (strcmp(name, "ArcSin") == 0) {
        /* ArcSin[u] coefficients via recurrence
         * c_{2m+1} = c_{2m-1} * (2m-1)^2 / (2m * (2m+1)), c_1 = 1. */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(1);
        for (size_t k = 2; k < N; k++) {
            if ((k & 1) == 0) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 1) / 2;
            int64_t num = (2*m - 1) * (2*m - 1);
            int64_t den = 2 * m * (2*m + 1);
            Expr* ratio = make_rational(num, den);
            c[k] = simp(mk_times(expr_copy(c[k - 2]), ratio));
        }
    } else if (strcmp(name, "ArcSinh") == 0) {
        /* Same magnitudes as ArcSin with alternating sign on odd k. */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(1);
        for (size_t k = 2; k < N; k++) {
            if ((k & 1) == 0) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 1) / 2;
            int64_t num = -(2*m - 1) * (2*m - 1);  /* flip sign */
            int64_t den = 2 * m * (2*m + 1);
            Expr* ratio = make_rational(num, den);
            c[k] = simp(mk_times(expr_copy(c[k - 2]), ratio));
        }
    } else if (strcmp(name, "SinIntegral") == 0) {
        /* Si[u] = sum_{m>=0} (-1)^m u^(2m+1) / ((2m+1)(2m+1)!).
         * Odd powers only; c_{2m+1} = c_{2m-1} * -(2m-1)/((2m)(2m+1)^2), c_1 = 1. */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(1);
        for (size_t k = 2; k < N; k++) {
            if ((k & 1) == 0) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 1) / 2;
            int64_t num = -(2*m - 1);
            int64_t den = (2*m) * (2*m + 1) * (2*m + 1);
            Expr* ratio = make_rational(num, den);
            c[k] = simp(mk_times(expr_copy(c[k - 2]), ratio));
        }
    } else if (strcmp(name, "SinhIntegral") == 0) {
        /* Shi[u] = sum_{m>=0} u^(2m+1) / ((2m+1)(2m+1)!).  Same as Si without the
         * alternating sign; c_{2m+1} = c_{2m-1} * (2m-1)/((2m)(2m+1)^2), c_1 = 1. */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(1);
        for (size_t k = 2; k < N; k++) {
            if ((k & 1) == 0) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 1) / 2;
            int64_t num = (2*m - 1);
            int64_t den = (2*m) * (2*m + 1) * (2*m + 1);
            Expr* ratio = make_rational(num, den);
            c[k] = simp(mk_times(expr_copy(c[k - 2]), ratio));
        }
    } else if (strcmp(name, "FresnelC") == 0) {
        /* C[u] = sum_{m>=0} (-1)^m (Pi/2)^(2m) u^(4m+1) / ((2m)!(4m+1)).
         * Powers 4m+1 only; c_1 = 1,
         * c_{4m+1} = c_{4m-3} * Pi^2 * -(4m-3)/(4(2m-1)(2m)(4m+1)). */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(1);
        for (size_t k = 2; k < N; k++) {
            if (k % 4 != 1) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 1) / 4;
            int64_t num = -(4*m - 3);
            int64_t den = 4 * (2*m - 1) * (2*m) * (4*m + 1);
            Expr* ratio = simp(mk_times(mk_power(mk_symbol("Pi"), expr_new_integer(2)),
                                        make_rational(num, den)));
            c[k] = simp(mk_times(expr_copy(c[k - 4]), ratio));
        }
    } else if (strcmp(name, "FresnelS") == 0) {
        /* S[u] = sum_{m>=0} (-1)^m (Pi/2)^(2m+1) u^(4m+3) / ((2m+1)!(4m+3)).
         * Powers 4m+3 only; c_3 = Pi/6,
         * c_{4m+3} = c_{4m-1} * Pi^2 * -(4m-1)/(4(2m)(2m+1)(4m+3)). */
        c[0] = expr_new_integer(0);
        if (N > 1) c[1] = expr_new_integer(0);
        if (N > 2) c[2] = expr_new_integer(0);
        if (N > 3) c[3] = simp(mk_times(make_rational(1, 6), mk_symbol("Pi")));
        for (size_t k = 4; k < N; k++) {
            if (k % 4 != 3) { c[k] = expr_new_integer(0); continue; }
            int64_t m = (int64_t)(k - 3) / 4;
            int64_t num = -(4*m - 1);
            int64_t den = 4 * (2*m) * (2*m + 1) * (4*m + 3);
            Expr* ratio = simp(mk_times(mk_power(mk_symbol("Pi"), expr_new_integer(2)),
                                        make_rational(num, den)));
            c[k] = simp(mk_times(expr_copy(c[k - 4]), ratio));
        }
    } else {
        for (size_t i = 0; i < N; i++) c[i] = expr_new_integer(0);
    }
    return c;
}

static void kernel_coefs_free(Expr** c, size_t N) {
    if (!c) return;
    for (size_t i = 0; i < N; i++) if (c[i]) expr_free(c[i]);
    free(c);
}

/* Split s into constant part and non-constant-part.
 * On return, *c is the coefficient at (x-x0)^0 (or 0 if nmin > 0), and
 * *u is a new series equal to s minus that constant (nmin >= 1 after).
 *
 * Works on a trimmed copy of the input so spurious leading zeros (which
 * routinely arise from e.g. Sin[u]/u, whose so_mul output has coef[0]=0
 * at exp=-1 even though the true series starts at exp=0) don't look like
 * a Laurent input to the caller. */
static void so_split_constant(SeriesObj* s_in, Expr** c_out, SeriesObj** u_out) {
    SeriesObj* s = so_copy_trimmed(s_in);
    Expr* c = expr_new_integer(0);
    SeriesObj* u;
    if (s->nmin > 0 || s->coef_count == 0) {
        /* No constant term; u = s. */
        u = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) so_set_coef(u, i, expr_copy(s->coefs[i]));
    } else if (s->nmin < 0) {
        /* True Laurent: splitting the constant is not meaningful. Return s
         * as u with c = 0 so the caller's `u->nmin < 1` check rejects it. */
        u = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) so_set_coef(u, i, expr_copy(s->coefs[i]));
    } else {
        /* nmin == 0. */
        expr_free(c);
        c = expr_copy(s->coefs[0]);
        u = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) {
            if (i == 0) so_set_coef(u, i, expr_new_integer(0));
            else        so_set_coef(u, i, expr_copy(s->coefs[i]));
        }
        so_trim_leading(u);
    }
    so_free(s);
    *c_out = c;
    *u_out = u;
}

/* Apply Exp to a series. Exp[c + u] = Exp[c] * Exp[u]. */
static SeriesObj* so_apply_exp(SeriesObj* s) {
    Expr* c; SeriesObj* u;
    so_split_constant(s, &c, &u);
    int64_t N = (u->order - u->nmin) / (u->nmin > 0 ? u->nmin : 1) + 1;
    if (N < 1) N = 1;
    if (u->nmin < 1) { expr_free(c); so_free(u); return NULL; }
    Expr** k = kernel_coefs("Exp", (size_t)N);
    SeriesObj* ex_u = so_compose_scalar_kernel(k, (size_t)N, u);
    kernel_coefs_free(k, (size_t)N);
    Expr* ex_c = simp(mk_fn1("Exp", c));
    SeriesObj* r = so_scalar_mul(ex_c, ex_u);
    expr_free(ex_c); so_free(u); so_free(ex_u);
    return r;
}

/* Apply Log to a series. Log[c + u] = Log[c] + Log[1 + u/c] when c != 0. */
static SeriesObj* so_apply_log(SeriesObj* s) {
    Expr* c; SeriesObj* u;
    so_split_constant(s, &c, &u);
    if (is_lit_zero(c)) {
        /* Handled by the recursion above via the Power rewrite; here we
         * simply reject. */
        expr_free(c); so_free(u); return NULL;
    }
    if (u->nmin < 1) { expr_free(c); so_free(u); return NULL; }
    Expr* inv_c = simp(mk_power(expr_copy(c), expr_new_integer(-1)));
    SeriesObj* u_over_c = so_scalar_mul(inv_c, u);
    expr_free(inv_c);
    int64_t N = (u_over_c->order - u_over_c->nmin) / (u_over_c->nmin > 0 ? u_over_c->nmin : 1) + 2;
    if (N < 2) N = 2;
    Expr** k = kernel_coefs("Log1p", (size_t)N);
    SeriesObj* log_term = so_compose_scalar_kernel(k, (size_t)N, u_over_c);
    kernel_coefs_free(k, (size_t)N);
    so_free(u_over_c); so_free(u);
    Expr* log_c = simp(mk_fn1("Log", c));
    SeriesObj* c_series = so_from_constant(log_c, log_term->x, log_term->x0, log_term->order, log_term->den);
    expr_free(log_c);
    SeriesObj* r = so_add(c_series, log_term);
    so_free(c_series); so_free(log_term);
    return r;
}

/* Sin[c + u] = Sin[c] Cos[u] + Cos[c] Sin[u]. */
static SeriesObj* so_apply_sin_or_cos(SeriesObj* s, bool is_sin) {
    Expr* c; SeriesObj* u;
    so_split_constant(s, &c, &u);
    if (u->nmin < 1) { expr_free(c); so_free(u); return NULL; }
    int64_t N = (u->order - u->nmin) / (u->nmin > 0 ? u->nmin : 1) + 2;
    if (N < 2) N = 2;
    Expr** ks = kernel_coefs("Sin", (size_t)N);
    Expr** kc = kernel_coefs("Cos", (size_t)N);
    SeriesObj* su = so_compose_scalar_kernel(ks, (size_t)N, u);
    SeriesObj* cu = so_compose_scalar_kernel(kc, (size_t)N, u);
    kernel_coefs_free(ks, (size_t)N);
    kernel_coefs_free(kc, (size_t)N);
    Expr* sinc = simp(mk_fn1("Sin", expr_copy(c)));
    Expr* cosc = simp(mk_fn1("Cos", c));
    SeriesObj* t1, *t2;
    if (is_sin) { t1 = so_scalar_mul(sinc, cu); t2 = so_scalar_mul(cosc, su); }
    else        { t1 = so_scalar_mul(cosc, cu); t2 = so_scalar_mul(sinc, su); }
    expr_free(sinc); expr_free(cosc); so_free(u); so_free(su); so_free(cu);
    SeriesObj* r;
    if (is_sin) { r = so_add(t1, t2); }
    else        { r = so_sub(t1, t2); }
    so_free(t1); so_free(t2);
    return r;
}

/* Apply any odd-series kernel whose series at 0 is known analytically
 * (ArcTan, ArcTanh, ArcSin, ArcSinh) to a series s with zero constant term.
 * These are all analytic at u=0 with the elementary function's value 0
 * there, and their series in u at c = s(x0) requires computing the value
 * and derivatives of f at c. For c != 0 we fall back to evaluating
 * f at c symbolically and using the chain rule via the D-based approach
 * (which the caller can do). Here we only handle c = 0 directly; for
 * c != 0 we return NULL so the dispatcher can pick another path.
 */
static SeriesObj* so_apply_kernel_at_zero(const char* name, SeriesObj* s) {
    Expr* c; SeriesObj* u;
    so_split_constant(s, &c, &u);
    bool c_is_zero = is_lit_zero(c);
    expr_free(c);
    if (!c_is_zero) { so_free(u); return NULL; }
    if (u->nmin < 1) { so_free(u); return NULL; }
    int64_t N = (u->order - u->nmin) / (u->nmin > 0 ? u->nmin : 1) + 2;
    if (N < 2) N = 2;
    Expr** k = kernel_coefs(name, (size_t)N);
    SeriesObj* r = so_compose_scalar_kernel(k, (size_t)N, u);
    kernel_coefs_free(k, (size_t)N);
    so_free(u);
    return r;
}

/* Detect the constant term of `inner` being a literal integer +1 or -1.
 * Returns +1, -1, or 0 (none). Used to dispatch Puiseux branch-point
 * expansions for ArcSin / ArcCos when the at-zero kernel can't apply. */
static int so_branch_point_sign(SeriesObj* inner) {
    SeriesObj* t = so_copy_trimmed(inner);
    int result = 0;
    if (t->nmin == 0 && t->coef_count >= 1) {
        Expr* c = t->coefs[0];
        if (c->type == EXPR_INTEGER) {
            if (c->data.integer == 1) result = 1;
            else if (c->data.integer == -1) result = -1;
        }
    }
    so_free(t);
    return result;
}

/* Forward declarations for the shared branch-point helpers defined below
 * (after the so_apply_arccos block). so_apply_arc_branch_point — refactored
 * to share the same finalisation path as the new Family A / Family B
 * handlers — calls into them. */
static bool branch_check_linear_inner(SeriesObj* inner, int sign_c,
                                      bool at_imag, Expr** q_out);
static SeriesObj* family_a_alloc(SeriesCtx* ctx, int64_t user_order);
static SeriesObj* family_a_finalize(SeriesObj* puiseux, Expr* add_const,
                                    SeriesCtx* ctx);

/* Puiseux branch-point expansion for ArcSin[c + q*w] or ArcCos[c + q*w],
 * where c = ±1, q is any non-zero Expr free of w, and the inner series is
 * exactly c + q*w (i.e. a constant + a single linear term; higher-order
 * coefficients must be zero).
 *
 * Derivation. Starting from the integral identity
 *   ArcCos[1-s] = Sqrt[2s] * sum_{k>=0} b_k / (2k+1) * s^k,
 *   b_k = (2k)! / (8^k (k!)^2)
 * and ArcSin[x] = Pi/2 - ArcCos[x], plus odd symmetry ArcSin[-y] = -ArcSin[y]
 * and ArcCos[-y] = Pi - ArcCos[y] for the c=-1 case, we get:
 *
 *   at c=+1, with s = -q*w:
 *     ArcCos[1+q*w] = Sqrt[-2q]*Sqrt[w] * sum (-q)^k * b_k/(2k+1) * w^k
 *     ArcSin[1+q*w] = Pi/2 - ArcCos[1+q*w]
 *   at c=-1, with s = q*w:
 *     ArcCos[-1+q*w] = Pi - Sqrt[2q]*Sqrt[w] * sum q^k * b_k/(2k+1) * w^k
 *     ArcSin[-1+q*w] = -Pi/2 + Sqrt[2q]*Sqrt[w] * sum q^k * b_k/(2k+1) * w^k
 *
 * Output uses den = 2 (half-integer exponents in w) with the constant term
 * at exponent 0 and the k-th Puiseux coefficient at exponent (2k+1)/2.
 */
static SeriesObj* so_apply_arc_branch_point(SeriesObj* inner, int sign_c,
                                            bool is_arcsin, int64_t user_order,
                                            SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/false, &q)) return NULL;

    SeriesObj* r = family_a_alloc(ctx, user_order);

    /* s = -sign_c * q, so Sqrt[2s] is real for q in the right domain. */
    Expr* s_expr = (sign_c == 1)
        ? simp(mk_times(expr_new_integer(-1), expr_copy(q)))
        : expr_copy(q);
    Expr* two_s = simp(mk_times(expr_new_integer(2), expr_copy(s_expr)));
    Expr* sqrt_2s = simp(mk_power(two_s, make_rational(1, 2)));
    /* sign_adj: -1 for ArcSin@+1 or ArcCos@-1, +1 otherwise. */
    int sign_adj = is_arcsin ? -sign_c : sign_c;
    Expr* prefactor = (sign_adj == -1)
        ? simp(mk_times(expr_new_integer(-1), sqrt_2s))
        : sqrt_2s;

    /* r_k = b_k / (2k+1) with b_0 = 1, b_{k+1}/b_k = (2k+1)(2k+2) / (8(k+1)^2).
     * Chained with the 1/(2k+1) factor: r_{k+1}/r_k = (2k+1)^2 / (4(k+1)(2k+3)). */
    Expr* s_pow = expr_new_integer(1);
    Expr* r_k   = expr_new_integer(1);
    int64_t max_k = (r->order - 1) / 2;
    for (int64_t k = 0; k <= max_k; k++) {
        size_t idx = (size_t)(2*k + 1);
        if (idx >= r->coef_count) break;
        Expr* inner_prod = simp(mk_times(expr_copy(s_pow), expr_copy(r_k)));
        Expr* coef = simp(mk_times(expr_copy(prefactor), inner_prod));
        so_set_coef(r, idx, coef);
        /* Advance the coefficient recurrence. */
        int64_t num = (2*k + 1) * (2*k + 1);
        int64_t den_val = 4 * (k + 1) * (2*k + 3);
        Expr* ratio = make_rational(num, den_val);
        r_k = simp(mk_times(r_k, ratio));
        /* Advance s_pow = s^k -> s^(k+1). */
        Expr* s_next = simp(mk_times(s_pow, expr_copy(s_expr)));
        s_pow = s_next;
    }
    expr_free(s_pow);
    expr_free(r_k);
    expr_free(s_expr);
    expr_free(q);
    expr_free(prefactor);

    /* Additive constant at u^0: ±Pi/2 (ArcSin) or 0/Pi (ArcCos). */
    Expr* const_term;
    if (is_arcsin) {
        const_term = simp(mk_times(make_rational(sign_c, 2), mk_symbol("Pi")));
    } else {
        const_term = (sign_c == -1) ? mk_symbol("Pi") : expr_new_integer(0);
    }
    return family_a_finalize(r, const_term, ctx);
}

/* ============================================================================
 * Branch-point expansion for inverse trig / hyperbolic functions
 *
 * Two mathematical families:
 *   Family A (square-root): derivative has 1/Sqrt[(x-x0)*linear]
 *      ArcSin / ArcCos at x0 = ±1, ArcSinh at x0 = ±I,  ArcCosh at x0 = ±1.
 *      Output is a Puiseux series with denom 2 (half-integer exponents).
 *   Family B (logarithmic):  derivative has a simple pole at x0
 *      ArcTan / ArcCot at x0 = ±I,  ArcTanh / ArcCoth at x0 = ±1.
 *      Output is `c_log * Log[x-x0] + regular_power_series`.
 *
 * Each handler operates in one of two modes, selected by ctx->allow_branch_wrap:
 *   - Wrap mode (top-level call from do_series_single): stashes the additive
 *     constant, log coefficient, and (-1)^Floor[...] branch discriminator into
 *     ctx->pending_*, returns a SeriesObj with coef[0] = 0. The outermost
 *     epilogue builds the MMA-faithful
 *       Plus[add_const, log_coef*Log[x-x0], Times[disc, SeriesData[...]]]
 *     wrapper around the returned SeriesData.
 *   - Inline mode (composed under another head): returns a SeriesObj with
 *     coef[0] = add_const (+ log_coef*Log[x-x0] for Family B) so callers
 *     like so_apply_sin can compose against the full value. The branch
 *     discriminator is dropped — composed branch cases give a principal-
 *     branch numerical answer without the Floor wrapper.
 * ========================================================================== */

/* Build the imaginary unit times an integer sign: Complex[0, sign]. */
static Expr* make_imag_unit_signed(int sign) {
    Expr* args[2] = { expr_new_integer(0), expr_new_integer(sign) };
    return expr_new_function(mk_symbol("Complex"), args, 2);
}

/* True iff e is structurally equal to sign * I (i.e. Complex[0, sign]
 * after canonicalisation). */
static bool is_imag_unit(Expr* e, int sign) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol.name != SYM_Complex) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* re = e->data.function.args[0];
    Expr* im = e->data.function.args[1];
    if (!is_lit_zero(re)) return false;
    return im->type == EXPR_INTEGER && im->data.integer == sign;
}

/* Detect that the constant term of `inner` is ±I and return that sign.
 * Returns 0 if not a clean ±I literal. Mirror of so_branch_point_sign. */
static int so_branch_point_imag_sign(SeriesObj* inner) {
    SeriesObj* t = so_copy_trimmed(inner);
    int result = 0;
    if (t->nmin == 0 && t->coef_count >= 1) {
        Expr* c = t->coefs[0];
        if (is_imag_unit(c, 1))       result = 1;
        else if (is_imag_unit(c, -1)) result = -1;
    }
    so_free(t);
    return result;
}

/* Build x - x0 as a (simplified) Expr (caller owns). */
static Expr* build_x_minus_x0(Expr* x, Expr* x0) {
    if (is_lit_zero(x0)) return expr_copy(x);
    Expr* neg = simp(mk_times(expr_new_integer(-1), expr_copy(x0)));
    return simp(mk_plus(expr_copy(x), neg));
}

/* MMA-style branch discriminator (-1)^Floor[(Pi/2 - Arg[x - x0])/(2*Pi)].
 * For x on the principal sheet near x0 the Floor evaluates to 0 and the
 * factor is +1; off-sheet evaluation crosses the branch and the factor
 * flips sign. We use a uniform reference angle of Pi/2 for all eight
 * branch points; MMA varies this per head/branch but the principal-branch
 * result is identical (the Floor evaluates to 0 there in every case). */
static Expr* make_branch_discriminator(Expr* x, Expr* x0) {
    Expr* xmx0  = build_x_minus_x0(x, x0);
    Expr* arg_e = mk_fn1("Arg", xmx0);
    Expr* half_pi  = simp(mk_times(make_rational(1, 2), mk_symbol("Pi")));
    Expr* diff     = simp(mk_plus(half_pi,
                                  simp(mk_times(expr_new_integer(-1), arg_e))));
    Expr* two_pi   = simp(mk_times(expr_new_integer(2), mk_symbol("Pi")));
    Expr* inv_2pi  = simp(mk_power(two_pi, expr_new_integer(-1)));
    Expr* ratio    = simp(mk_times(diff, inv_2pi));
    Expr* floor_e  = mk_fn1("Floor", ratio);
    return mk_power(expr_new_integer(-1), floor_e);
}

/* Build `factor * Log[x - x0]` as a (simplified) Expr (caller owns). */
static Expr* build_log_term(Expr* factor, Expr* x, Expr* x0) {
    Expr* base = build_x_minus_x0(x, x0);
    Expr* log_e = simp(mk_fn1("Log", base));
    return simp(mk_times(expr_copy(factor), log_e));
}

/* Stash branch-wrapper metadata onto ctx. Takes ownership of all three
 * Expr arguments (any of which may be NULL). Frees any previous stash. */
static void ctx_set_branch_meta(SeriesCtx* ctx,
                                Expr* add_const,
                                Expr* log_coef,
                                Expr* discriminator) {
    if (ctx->pending_add_const)     expr_free(ctx->pending_add_const);
    if (ctx->pending_log_coef)      expr_free(ctx->pending_log_coef);
    if (ctx->pending_discriminator) expr_free(ctx->pending_discriminator);
    ctx->pending_add_const     = add_const;
    ctx->pending_log_coef      = log_coef;
    ctx->pending_discriminator = discriminator;
}

/* Recursive series_expand call that does NOT inherit the parent's
 * allow_branch_wrap. Used by composition paths (Plus, Times, Power, Log,
 * elementary-head inner argument expansion) so that a nested branch-point
 * case is silently produced as a constant-inside SeriesObj (composition-
 * friendly) rather than corrupting the outer wrap state. */
static SeriesObj* series_expand_nested(Expr* e, SeriesCtx* ctx) {
    bool saved = ctx->allow_branch_wrap;
    ctx->allow_branch_wrap = false;
    SeriesObj* r = series_expand(e, ctx);
    ctx->allow_branch_wrap = saved;
    return r;
}

/* Assemble a Family A handler result from a Puiseux series (built with
 * coef[0] zero and Puiseux terms at coef[1..]) and the additive constant
 * f(x0). In wrap mode, stash the constant + discriminator on ctx and
 * return the SeriesObj as-is. In inline mode, place the constant at
 * coef[0] for composition. The series must already have nmin == 0 and
 * den == 2. Takes ownership of add_const. */
static SeriesObj* family_a_finalize(SeriesObj* puiseux, Expr* add_const,
                                    SeriesCtx* ctx) {
    if (ctx->allow_branch_wrap) {
        Expr* disc = make_branch_discriminator(ctx->x, ctx->x0);
        ctx_set_branch_meta(ctx, add_const, NULL, disc);
        /* coef[0] stays 0; wrapper builds Plus[add_const, disc*SeriesData]. */
        return puiseux;
    }
    /* Inline mode: bake the constant at coef[0] for composition. */
    if (puiseux->coef_count > 0) {
        so_set_coef(puiseux, 0, add_const);
    } else {
        expr_free(add_const);
    }
    return puiseux;
}

/* Assemble a Family B handler result from a regular power series (built
 * with coef[0] zero and regular coefs at coef[1..], den == 1), the
 * additive constant f(x0) (finite part), and the coefficient of the
 * Log[x-x0] singularity. In wrap mode stash all three; in inline mode
 * bake (add_const + log_coef*Log[x-x0]) into coef[0]. Takes ownership of
 * add_const and log_coef. */
static SeriesObj* family_b_finalize(SeriesObj* regular,
                                    Expr* add_const, Expr* log_coef,
                                    SeriesCtx* ctx) {
    if (ctx->allow_branch_wrap) {
        Expr* disc = make_branch_discriminator(ctx->x, ctx->x0);
        ctx_set_branch_meta(ctx, add_const, log_coef, disc);
        return regular;
    }
    /* Inline mode: coef[0] gets add_const + log_coef * Log[x - x0]. */
    Expr* log_term = build_log_term(log_coef, ctx->x, ctx->x0);
    Expr* combined = simp(mk_plus(add_const, log_term));
    expr_free(log_coef);
    if (regular->coef_count > 0) {
        so_set_coef(regular, 0, combined);
    } else {
        expr_free(combined);
    }
    return regular;
}

/* Pre-flight check shared by all new branch-point handlers: `inner` must
 * decompose exactly as `sign_c*branch_val + q*w`, with all higher
 * coefficients zero. Returns true and yields q (caller owns) on success.
 * `branch_val` is either Integer 1 (real branch) or Complex[0, 1] (imag).
 * `inner` is consumed only for inspection; the original is left intact. */
static bool branch_check_linear_inner(SeriesObj* inner, int sign_c,
                                      bool at_imag, Expr** q_out) {
    *q_out = NULL;
    SeriesObj* t = so_copy_trimmed(inner);
    bool ok = (t->nmin == 0 && t->coef_count >= 2);
    if (ok) {
        Expr* c = t->coefs[0];
        if (at_imag) ok = is_imag_unit(c, sign_c);
        else         ok = (c->type == EXPR_INTEGER && c->data.integer == sign_c);
    }
    if (ok) {
        for (size_t i = 2; i < t->coef_count; i++) {
            if (!is_lit_zero(t->coefs[i])) { ok = false; break; }
        }
    }
    if (ok) *q_out = expr_copy(t->coefs[1]);
    so_free(t);
    if (ok && is_lit_zero(*q_out)) { expr_free(*q_out); *q_out = NULL; return false; }
    return ok;
}

/* Allocate a Puiseux SeriesObj for a Family A handler. nmin = 0, den = 2,
 * coef[0] = 0 (caller fills puiseux coefs at coef[1], coef[3], ...). */
static SeriesObj* family_a_alloc(SeriesCtx* ctx, int64_t user_order) {
    int64_t new_order = (user_order + 1) * 2;
    if (new_order < 2) new_order = 2;
    return so_alloc(ctx->x, ctx->x0, 0, new_order, 2);
}

/* Allocate a regular SeriesObj for a Family B handler. nmin = 0, den = 1,
 * coef[0] = 0 (caller fills regular coefs at coef[1], coef[2], ...). */
static SeriesObj* family_b_alloc(SeriesCtx* ctx, int64_t user_order) {
    int64_t new_order = user_order + 1;
    if (new_order < 1) new_order = 1;
    return so_alloc(ctx->x, ctx->x0, 0, new_order, 1);
}

/* ----------------------------------------------------------------------------
 * Family A — Puiseux at branch points
 * ---------------------------------------------------------------------------- */

/* ArcSinh[sign_c*I + q*w] = sign_c*I*Pi/2 + 2 * ArcSinh[Sqrt[q*w/(2*sign_c*I)]].
 * Reuses kernel_coefs("ArcSinh", ...) for the inner ArcSinh expansion.
 *
 * Leading Puiseux coefficient (matches MMA's "1+I" * (-I) = 1-I for sign=+1):
 *   c_0 = 2 / Sqrt[2*sign_c*I] * q^{1/2} = Sqrt[-2*sign_c*I*q]
 */
static SeriesObj* so_apply_arcsinh_branch_point(SeriesObj* inner, int sign_c,
                                                int64_t user_order,
                                                SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/true, &q)) return NULL;

    SeriesObj* r = family_a_alloc(ctx, user_order);
    /* Use ArcSinh's odd kernel coefficients: ArcSinh[w] = sum_{j odd} a_j w^j.
     * With w = Sqrt[q*v/(2*sign_c*I)] where v = w (the series variable),
     * and the outer factor 2, we get terms at v^{(j+1)/2 - 1/2} ... wait.
     * Cleaner: c_k = coefficient of v^((2k+1)/2) in 2*ArcSinh[Sqrt[v/(2c)]]
     *             = 2 * arcsinh[2k+1] * (1/(2c))^{k+1/2}
     *             = 2*arcsinh[2k+1] * (1/(2c))^k * (1/(2c))^{1/2}
     * With an outer monomial q*v -> coefficient at v^((2k+1)/2) scales by
     * q^{k+1/2}: c_k(q) = 2*arcsinh[2k+1] * (q/(2c))^k * (q/(2c))^{1/2}
     * where c = sign_c*I. Set base = q/(2*sign_c*I) = -sign_c*I*q/2.
     */
    int64_t max_k = (r->order - 1) / 2;
    size_t N = (size_t)(2 * max_k + 2);
    if (N < 2) N = 2;
    Expr** kc = kernel_coefs("ArcSinh", N);

    /* base = q / (2 * sign_c * I) = -sign_c * I * q / 2 */
    Expr* iconst = make_imag_unit_signed(-sign_c);
    Expr* iq     = simp(mk_times(iconst, expr_copy(q)));
    Expr* base   = simp(mk_times(make_rational(1, 2), iq));

    /* sqrt_base = Sqrt[base] -- leading factor */
    Expr* sqrt_base = simp(mk_power(expr_copy(base), make_rational(1, 2)));
    Expr* prefactor = simp(mk_times(expr_new_integer(2), sqrt_base));

    Expr* base_pow = expr_new_integer(1);  /* base^0 */
    for (int64_t k = 0; k <= max_k; k++) {
        size_t idx_in_r = (size_t)(2*k + 1);
        if (idx_in_r >= r->coef_count) break;
        /* kc[2k+1] is the ArcSinh coefficient at u^{2k+1}. */
        Expr* ker = (2*(size_t)k + 1 < N) ? kc[2*k + 1] : expr_new_integer(0);
        Expr* coef = simp(mk_times(expr_copy(prefactor),
                          simp(mk_times(expr_copy(ker), expr_copy(base_pow)))));
        so_set_coef(r, idx_in_r, coef);
        /* Advance base_pow *= base for next iteration. */
        Expr* next = simp(mk_times(base_pow, expr_copy(base)));
        base_pow = next;
    }
    expr_free(base_pow);
    expr_free(base);
    expr_free(prefactor);
    expr_free(q);
    kernel_coefs_free(kc, N);

    /* f(x0) = sign_c * I * Pi/2 = sign_c * Complex[0, 1/2] * Pi. */
    Expr* iconst2 = expr_new_function(mk_symbol("Complex"),
                        (Expr*[]){ expr_new_integer(0), make_rational(sign_c, 2) }, 2);
    Expr* add_const = simp(mk_times(iconst2, mk_symbol("Pi")));
    return family_a_finalize(r, add_const, ctx);
}

/* ArcCosh[sign_c + q*w]:
 *   sign_c = +1:  0 + 2*ArcSinh[Sqrt[q*w/2]]
 *   sign_c = -1:  I*Pi + 2*ArcSinh[Sqrt[q*w/(-2)]] - 2*ArcSinh[0] ... no.
 *
 * Use the identities (principal branch):
 *   ArcCosh[1 + u]  = 2*ArcSinh[Sqrt[u/2]]
 *   ArcCosh[-1 + u] = I*Pi - ArcCosh[1 - u] = I*Pi - 2*ArcSinh[Sqrt[-u/2]]
 *
 * Same Sqrt[base] * 2 * ArcSinh kernel as ArcSinh@±I; the only differences
 * are the additive constant and the sign of base. */
static SeriesObj* so_apply_arccosh_branch_point(SeriesObj* inner, int sign_c,
                                                int64_t user_order,
                                                SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/false, &q)) return NULL;

    SeriesObj* r = family_a_alloc(ctx, user_order);
    int64_t max_k = (r->order - 1) / 2;
    size_t N = (size_t)(2 * max_k + 2);
    if (N < 2) N = 2;
    Expr** kc = kernel_coefs("ArcSinh", N);

    /* base = q/2 (sign_c=+1) or -q/2 (sign_c=-1, with overall minus sign on
     * the prefactor: -2*ArcSinh[Sqrt[-u/2]] ... but Sqrt[-u/2] for u>0 real
     * is imaginary; the i absorbs into the leading Sqrt and yields the
     * correct principal-branch result.) */
    Expr* base = (sign_c == 1)
        ? simp(mk_times(make_rational(1, 2), expr_copy(q)))
        : simp(mk_times(make_rational(-1, 2), expr_copy(q)));

    Expr* sqrt_base = simp(mk_power(expr_copy(base), make_rational(1, 2)));
    Expr* sign_factor = (sign_c == 1) ? expr_new_integer(2) : expr_new_integer(-2);
    Expr* prefactor = simp(mk_times(sign_factor, sqrt_base));

    Expr* base_pow = expr_new_integer(1);
    for (int64_t k = 0; k <= max_k; k++) {
        size_t idx_in_r = (size_t)(2*k + 1);
        if (idx_in_r >= r->coef_count) break;
        Expr* ker = (2*(size_t)k + 1 < N) ? kc[2*k + 1] : expr_new_integer(0);
        Expr* coef = simp(mk_times(expr_copy(prefactor),
                          simp(mk_times(expr_copy(ker), expr_copy(base_pow)))));
        so_set_coef(r, idx_in_r, coef);
        Expr* next = simp(mk_times(base_pow, expr_copy(base)));
        base_pow = next;
    }
    expr_free(base_pow);
    expr_free(base);
    expr_free(prefactor);
    expr_free(q);
    kernel_coefs_free(kc, N);

    /* f(x0): 0 at +1, I*Pi at -1. */
    Expr* add_const;
    if (sign_c == 1) {
        add_const = expr_new_integer(0);
    } else {
        Expr* iconst = make_imag_unit_signed(1);
        add_const = simp(mk_times(iconst, mk_symbol("Pi")));
    }
    return family_a_finalize(r, add_const, ctx);
}

/* ----------------------------------------------------------------------------
 * Family B — Logarithmic at branch points
 *
 * Each function is rewritten near its branch point in the form
 *     H[x0 + u] = add_const + log_coef * Log[u] + (factor) * Log[1 ± k*u]
 * where the final Log[1 ± k*u] is expanded by the Log1p kernel. The first
 * two terms become the add_const / log_coef metadata; the last becomes the
 * regular power-series part with coef[0] = 0 by construction.
 * ---------------------------------------------------------------------------- */

/* Compose `factor * Log[1 + alpha*u]` as a regular power series in u
 * (truncated to user_order). u corresponds to (x - x0) with den = 1.
 * alpha can be any Expr (Integer, Complex, etc.) free of u.
 *
 * Returns a SeriesObj with nmin = 0, den = 1, coef[0] = 0, coef[k>=1] =
 * factor * (-1)^(k+1) * alpha^k / k. */
static SeriesObj* family_b_build_log1p_part(SeriesCtx* ctx, int64_t user_order,
                                            Expr* factor, Expr* alpha) {
    SeriesObj* r = family_b_alloc(ctx, user_order);
    Expr* alpha_pow = expr_new_integer(1);  /* alpha^0 */
    for (int64_t k = 1; k < r->order; k++) {
        Expr* next = simp(mk_times(alpha_pow, expr_copy(alpha)));
        alpha_pow = next;
        /* coef = factor * (-1)^(k+1) * alpha^k / k */
        int64_t sign = (k & 1) ? 1 : -1;
        Expr* sk = simp(mk_times(expr_new_integer(sign), expr_copy(alpha_pow)));
        Expr* over_k = simp(mk_times(sk, make_rational(1, k)));
        Expr* coef = simp(mk_times(expr_copy(factor), over_k));
        so_set_coef(r, (size_t)k, coef);
    }
    expr_free(alpha_pow);
    return r;
}

/* ArcTanh[1 + u] near u = 0:
 *   = (1/2) * (Log[2 + u] - Log[-u])
 *   = (Log[2]/2 + I*Pi/2) + (-1/2)*Log[u] + (1/2)*Log[1 + u/2]
 *
 * ArcTanh[-1 + u] near u = 0:
 *   = (1/2) * (Log[u] - Log[2 - u])
 *   = (-Log[2]/2) + (+1/2)*Log[u] + (-1/2)*Log[1 - u/2]
 */
static SeriesObj* so_apply_arctanh_branch_point(SeriesObj* inner, int sign_c,
                                                int64_t user_order,
                                                SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/false, &q)) return NULL;

    Expr* add_const;
    Expr* log_coef;
    Expr* reg_factor;
    Expr* reg_alpha;  /* alpha in Log[1 + alpha * u] */
    if (sign_c == 1) {
        /* add_const = Log[2]/2 + I*Pi/2 */
        Expr* half_log2 = simp(mk_times(make_rational(1, 2),
                                       simp(mk_fn1("Log", expr_new_integer(2)))));
        Expr* ipi2 = simp(mk_times(make_imag_unit_signed(1),
                                  simp(mk_times(make_rational(1, 2), mk_symbol("Pi")))));
        add_const = simp(mk_plus(half_log2, ipi2));
        log_coef  = make_rational(-1, 2);
        reg_factor = make_rational(1, 2);
        reg_alpha  = simp(mk_times(make_rational(1, 2), expr_copy(q)));
    } else {
        /* add_const = -Log[2]/2 */
        add_const = simp(mk_times(make_rational(-1, 2),
                                  simp(mk_fn1("Log", expr_new_integer(2)))));
        log_coef  = make_rational(1, 2);
        reg_factor = make_rational(-1, 2);
        reg_alpha  = simp(mk_times(make_rational(-1, 2), expr_copy(q)));
    }
    SeriesObj* r = family_b_build_log1p_part(ctx, user_order, reg_factor, reg_alpha);
    expr_free(reg_factor); expr_free(reg_alpha); expr_free(q);
    return family_b_finalize(r, add_const, log_coef, ctx);
}

/* ArcCoth[1 + u]  = (1/2) * (Log[2 + u] - Log[u])  (both real for u > 0)
 *                 = Log[2]/2 + (-1/2)*Log[u] + (1/2)*Log[1 + u/2]
 *
 * ArcCoth[-1 + u] = (1/2) * (Log[u] - Log[-2 + u])
 *                 For u > 0 real: -2 + u is negative ⇒ Log[-2+u] = Log[2-u] + I*Pi
 *                 = -Log[2]/2 + I*Pi/2 + (+1/2)*Log[u] + (-1/2)*Log[1 - u/2]
 */
static SeriesObj* so_apply_arccoth_branch_point(SeriesObj* inner, int sign_c,
                                                int64_t user_order,
                                                SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/false, &q)) return NULL;

    Expr* add_const;
    Expr* log_coef;
    Expr* reg_factor;
    Expr* reg_alpha;
    if (sign_c == 1) {
        add_const = simp(mk_times(make_rational(1, 2),
                                  simp(mk_fn1("Log", expr_new_integer(2)))));
        log_coef  = make_rational(-1, 2);
        reg_factor = make_rational(1, 2);
        reg_alpha  = simp(mk_times(make_rational(1, 2), expr_copy(q)));
    } else {
        Expr* neg_half_log2 = simp(mk_times(make_rational(-1, 2),
                                            simp(mk_fn1("Log", expr_new_integer(2)))));
        Expr* ipi2 = simp(mk_times(make_imag_unit_signed(1),
                                  simp(mk_times(make_rational(1, 2), mk_symbol("Pi")))));
        add_const = simp(mk_plus(neg_half_log2, ipi2));
        log_coef  = make_rational(1, 2);
        reg_factor = make_rational(-1, 2);
        reg_alpha  = simp(mk_times(make_rational(-1, 2), expr_copy(q)));
    }
    SeriesObj* r = family_b_build_log1p_part(ctx, user_order, reg_factor, reg_alpha);
    expr_free(reg_factor); expr_free(reg_alpha); expr_free(q);
    return family_b_finalize(r, add_const, log_coef, ctx);
}

/* ArcTan[I + u]  = (Pi/4 + (I/2)*Log[2]) + (-I/2)*Log[u] + (I/2)*Log[1 - I*u/2]
 * ArcTan[-I + u] = (-3*Pi/4 - (I/2)*Log[2]) + (I/2)*Log[u] + (-I/2)*Log[1 + I*u/2]
 */
static SeriesObj* so_apply_arctan_branch_point(SeriesObj* inner, int sign_c,
                                               int64_t user_order,
                                               SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/true, &q)) return NULL;

    Expr* add_const;
    Expr* log_coef;
    Expr* reg_factor;
    Expr* reg_alpha;
    Expr* half_i_log2 = simp(mk_times(make_imag_unit_signed(1),
                              simp(mk_times(make_rational(1, 2),
                                            simp(mk_fn1("Log", expr_new_integer(2)))))));
    if (sign_c == 1) {
        Expr* pi4 = simp(mk_times(make_rational(1, 4), mk_symbol("Pi")));
        add_const = simp(mk_plus(pi4, half_i_log2));
        log_coef  = simp(mk_times(make_rational(-1, 2), make_imag_unit_signed(1)));
        reg_factor = simp(mk_times(make_rational(1, 2), make_imag_unit_signed(1)));
        /* alpha = -I/2 * q */
        Expr* mI = make_imag_unit_signed(-1);
        reg_alpha = simp(mk_times(make_rational(1, 2),
                                  simp(mk_times(mI, expr_copy(q)))));
    } else {
        Expr* neg3pi4 = simp(mk_times(make_rational(-3, 4), mk_symbol("Pi")));
        /* add_const = -3*Pi/4 - (I/2)*Log[2] = neg3pi4 + (-1)*half_i_log2 */
        Expr* neg_half_i_log2 = simp(mk_times(expr_new_integer(-1), half_i_log2));
        add_const = simp(mk_plus(neg3pi4, neg_half_i_log2));
        log_coef  = simp(mk_times(make_rational(1, 2), make_imag_unit_signed(1)));
        reg_factor = simp(mk_times(make_rational(-1, 2), make_imag_unit_signed(1)));
        /* alpha = +I/2 * q */
        Expr* pI = make_imag_unit_signed(1);
        reg_alpha = simp(mk_times(make_rational(1, 2),
                                  simp(mk_times(pI, expr_copy(q)))));
    }
    SeriesObj* r = family_b_build_log1p_part(ctx, user_order, reg_factor, reg_alpha);
    expr_free(reg_factor); expr_free(reg_alpha); expr_free(q);
    return family_b_finalize(r, add_const, log_coef, ctx);
}

/* ArcCot = Pi/2 - ArcTan (Mathilda convention; matches the at-zero rule
 * in so_apply_arccot). At a branch point the constants shift accordingly
 * and the log_coef / regular factors negate. */
static SeriesObj* so_apply_arccot_branch_point(SeriesObj* inner, int sign_c,
                                               int64_t user_order,
                                               SeriesCtx* ctx) {
    Expr* q;
    if (!branch_check_linear_inner(inner, sign_c, /*at_imag=*/true, &q)) return NULL;

    Expr* add_const;
    Expr* log_coef;
    Expr* reg_factor;
    Expr* reg_alpha;
    Expr* half_i_log2 = simp(mk_times(make_imag_unit_signed(1),
                              simp(mk_times(make_rational(1, 2),
                                            simp(mk_fn1("Log", expr_new_integer(2)))))));
    if (sign_c == 1) {
        /* Pi/2 - (Pi/4 + (I/2)*Log[2]) = Pi/4 - (I/2)*Log[2] */
        Expr* pi4 = simp(mk_times(make_rational(1, 4), mk_symbol("Pi")));
        Expr* neg_half_i_log2 = simp(mk_times(expr_new_integer(-1), half_i_log2));
        add_const = simp(mk_plus(pi4, neg_half_i_log2));
        log_coef  = simp(mk_times(make_rational(1, 2), make_imag_unit_signed(1)));
        reg_factor = simp(mk_times(make_rational(-1, 2), make_imag_unit_signed(1)));
        Expr* mI = make_imag_unit_signed(-1);
        reg_alpha = simp(mk_times(make_rational(1, 2),
                                  simp(mk_times(mI, expr_copy(q)))));
    } else {
        /* Pi/2 - (-3*Pi/4 - (I/2)*Log[2]) = 5*Pi/4 + (I/2)*Log[2] */
        Expr* fpi4 = simp(mk_times(make_rational(5, 4), mk_symbol("Pi")));
        add_const = simp(mk_plus(fpi4, half_i_log2));
        log_coef  = simp(mk_times(make_rational(-1, 2), make_imag_unit_signed(1)));
        reg_factor = simp(mk_times(make_rational(1, 2), make_imag_unit_signed(1)));
        Expr* pI = make_imag_unit_signed(1);
        reg_alpha = simp(mk_times(make_rational(1, 2),
                                  simp(mk_times(pI, expr_copy(q)))));
    }
    SeriesObj* r = family_b_build_log1p_part(ctx, user_order, reg_factor, reg_alpha);
    expr_free(reg_factor); expr_free(reg_alpha); expr_free(q);
    return family_b_finalize(r, add_const, log_coef, ctx);
}

/* Apply ArcCos[s] = Pi/2 - ArcSin[s]. Requires ArcSin expansion at s(x0). */
static SeriesObj* so_apply_arccos(SeriesObj* s) {
    SeriesObj* as = so_apply_kernel_at_zero("ArcSin", s);
    if (!as) return NULL;
    /* Negate and add Pi/2 constant. */
    SeriesObj* neg = so_neg(as);
    so_free(as);
    Expr* halfpi = simp(mk_times(make_rational(1, 2), mk_symbol("Pi")));
    SeriesObj* k = so_from_constant(halfpi, neg->x, neg->x0, neg->order, neg->den);
    expr_free(halfpi);
    SeriesObj* r = so_add(k, neg);
    so_free(k); so_free(neg);
    return r;
}

/* Apply ArcCot[s] = Pi/2 - ArcTan[s]. */
static SeriesObj* so_apply_arccot(SeriesObj* s) {
    SeriesObj* at = so_apply_kernel_at_zero("ArcTan", s);
    if (!at) return NULL;
    SeriesObj* neg = so_neg(at);
    so_free(at);
    Expr* halfpi = simp(mk_times(make_rational(1, 2), mk_symbol("Pi")));
    SeriesObj* k = so_from_constant(halfpi, neg->x, neg->x0, neg->order, neg->den);
    expr_free(halfpi);
    SeriesObj* r = so_add(k, neg);
    so_free(k); so_free(neg);
    return r;
}

/* Apply ArcCoth[s] = I*Pi/2 + ArcTanh[s] (principal-branch convention that
 * matches Mathematica's Series[ArcCoth[x], {x, 0, n}] output). */
static SeriesObj* so_apply_arccoth(SeriesObj* s) {
    SeriesObj* at = so_apply_kernel_at_zero("ArcTanh", s);
    if (!at) return NULL;
    /* I*Pi/2 = Complex[0, 1/2] * Pi  -- build as Times[Complex[0, 1/2], Pi]. */
    Expr* half = make_rational(1, 2);
    Expr* iconst = expr_new_function(mk_symbol("Complex"),
                        (Expr*[]){ expr_new_integer(0), half }, 2);
    Expr* ihp = simp(mk_times(iconst, mk_symbol("Pi")));
    SeriesObj* k = so_from_constant(ihp, at->x, at->x0, at->order, at->den);
    expr_free(ihp);
    SeriesObj* r = so_add(k, at);
    so_free(k); so_free(at);
    return r;
}

/* Apply ArcCosh[s] = I*ArcCos[s] (principal-branch identity that holds in a
 * neighbourhood of s(x0) = 0 and gives the series matching Mathematica's
 * Series[ArcCosh[x], {x, 0, n}] output: I*Pi/2 - I*x - I*x^3/6 - ...). */
static SeriesObj* so_apply_arccosh(SeriesObj* s) {
    SeriesObj* acs = so_apply_arccos(s);
    if (!acs) return NULL;
    Expr* iconst = expr_new_function(mk_symbol("Complex"),
                        (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
    SeriesObj* r = so_scalar_mul(iconst, acs);
    expr_free(iconst); so_free(acs);
    return r;
}

/* Build the series representing I*Pi/2 as a SeriesObj-compatible constant. */
static SeriesObj* so_const_half_i_pi(SeriesObj* tmpl) {
    Expr* half = make_rational(1, 2);
    Expr* iconst = expr_new_function(mk_symbol("Complex"),
                        (Expr*[]){ expr_new_integer(0), half }, 2);
    Expr* ihp = simp(mk_times(iconst, mk_symbol("Pi")));
    SeriesObj* k = so_from_constant(ihp, tmpl->x, tmpl->x0, tmpl->order, tmpl->den);
    expr_free(ihp);
    return k;
}

/* Sinh[c + u] = Sinh[c] Cosh[u] + Cosh[c] Sinh[u].
 * Cosh[c + u] = Cosh[c] Cosh[u] + Sinh[c] Sinh[u]. */
static SeriesObj* so_apply_sinh_or_cosh(SeriesObj* s, bool is_sinh) {
    Expr* c; SeriesObj* u;
    so_split_constant(s, &c, &u);
    if (u->nmin < 1) { expr_free(c); so_free(u); return NULL; }
    int64_t N = (u->order - u->nmin) / (u->nmin > 0 ? u->nmin : 1) + 2;
    if (N < 2) N = 2;
    Expr** ks = kernel_coefs("Sinh", (size_t)N);
    Expr** kc = kernel_coefs("Cosh", (size_t)N);
    SeriesObj* su = so_compose_scalar_kernel(ks, (size_t)N, u);
    SeriesObj* cu = so_compose_scalar_kernel(kc, (size_t)N, u);
    kernel_coefs_free(ks, (size_t)N);
    kernel_coefs_free(kc, (size_t)N);
    Expr* sinc = simp(mk_fn1("Sinh", expr_copy(c)));
    Expr* cosc = simp(mk_fn1("Cosh", c));
    SeriesObj* t1, *t2;
    if (is_sinh) { t1 = so_scalar_mul(sinc, cu); t2 = so_scalar_mul(cosc, su); }
    else         { t1 = so_scalar_mul(cosc, cu); t2 = so_scalar_mul(sinc, su); }
    expr_free(sinc); expr_free(cosc); so_free(u); so_free(su); so_free(cu);
    SeriesObj* r = so_add(t1, t2);
    so_free(t1); so_free(t2);
    return r;
}

/* ----------------------------------------------------------------------------
 * series_expand: recursive descent
 * -------------------------------------------------------------------------- */

/* SeriesCtx fields (full documentation; struct defined near SeriesObj):
 *   x, x0        - expansion variable & point (borrowed)
 *   order        - padded internal order; composite series arithmetic
 *                  (so_inv, so_compose_scalar_kernel, ...) needs headroom
 *                  beyond the user-facing order to survive cancellations.
 *   target_order - user-facing order (unpadded). Independent-coef paths
 *                  like series_taylor_via_D use this so they don't compute
 *                  derivatives they'll just truncate away again --
 *                  D[f, x, k=13] for ArcTan[x] at x=2 takes ~12s, versus
 *                  ~0.01s for k=3.
 *
 * MMA-faithful branch-point wrapper machinery:
 *   allow_branch_wrap     - true at the outermost call from do_series_single,
 *                           false in any nested composition (Plus/Times/Power
 *                           args, elementary-head inner). series_expand_nested
 *                           does the save/restore.
 *   pending_add_const     - H[x0]; populated by a branch-point handler in
 *                           wrap mode.
 *   pending_log_coef      - coefficient of Log[x - x0]; NULL for Family A.
 *   pending_discriminator - (-1)^Floor[(Pi/2 - Arg[x-x0])/(2*Pi)]; NULL
 *                           means no wrapper was emitted. The outermost
 *                           do_series_single inspects this to assemble
 *                             Plus[ add_const,
 *                                   Times[log_coef, Log[x-x0]],
 *                                   Times[discriminator, SeriesData[...]] ]
 */

/* Detect Infinity / ComplexInfinity / Indeterminate / DirectedInfinity
 * anywhere inside e. Used to bail out of naive Taylor before it spins. */
static bool has_infinity(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        return (s == SYM_Infinity ||
                s == SYM_ComplexInfinity ||
                s == SYM_Indeterminate);
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_SYMBOL) {
            const char* h = e->data.function.head->data.symbol.name;
            if (h == SYM_DirectedInfinity ||
                h == SYM_Indeterminate) return true;
        }
        if (has_infinity(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (has_infinity(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Naive Taylor via D[]: used for heads we don't recognise and for the
 * canonical Series[f[x], {x, a, n}] output shape. Computes coefficients
 * a_k = (D^k f at x0) / k! for k = 0..order-1.
 *
 * Bails out (returns NULL) if the evaluated f(x0) -- or any of its
 * derivatives -- contains Infinity, ComplexInfinity or Indeterminate, so
 * we do not spin when naive Taylor has hit a singularity of f. Also caps
 * the number of D iterations at MAX_NAIVE_ORDER to avoid the exponential
 * expression-size blow-up that hits non-trivial heads like ArcCos; tests
 * that truly need more terms should go through a direct kernel path. */
#define MAX_NAIVE_ORDER 20
static SeriesObj* series_taylor_via_D(Expr* e, SeriesCtx* ctx) {
    /* Taylor coefficients are computed independently, so there is no benefit
     * to running beyond the user's requested order -- the padding that
     * composite arithmetic (so_inv, so_compose_scalar_kernel) relies on is
     * discarded at the end. Cap at target_order + a couple of slack terms
     * in case a surrounding operation needs a touch more precision. */
    int64_t n_iter = ctx->target_order > 0 ? ctx->target_order + 2 : ctx->order;
    if (n_iter > ctx->order) n_iter = ctx->order;
    if (n_iter > MAX_NAIVE_ORDER) n_iter = MAX_NAIVE_ORDER;
    if (n_iter < 1) n_iter = 1;
    /* Quick singularity check before starting. The raw substitution
     * x -> x0 may leave an *unevaluated* pole (e.g. f[1/x] at x=0 becomes
     * f[1/0] = f[Power[0,-1]]), which has_infinity cannot see and which would
     * later spill a spurious `Power::infy: 1/0` to stderr once a coefficient
     * is simplified. Evaluate the probe (with arithmetic warnings muted, since
     * any 1/0 here is our exploratory substitution, not the user's) so the
     * pole collapses to ComplexInfinity / Indeterminate and is detected. */
    arith_warnings_mute_push();
    Expr* probe = eval_and_free(replace_all_of(e, ctx->x, ctx->x0));
    arith_warnings_mute_pop();
    bool bad = has_infinity(probe);
    expr_free(probe);
    if (bad) return NULL;

    SeriesObj* s = so_alloc(ctx->x, ctx->x0, 0, n_iter, 1);
    Expr* current = expr_copy(e);
    Expr* factorial = expr_new_integer(1);
    for (int64_t k = 0; k < n_iter; k++) {
        if (k > 0) factorial = simp(mk_times(factorial, expr_new_integer(k)));
        /* Evaluate f^(k)(x0) with warnings muted so an unevaluated pole in the
         * substitution collapses to ComplexInfinity / Indeterminate (detected
         * below) instead of leaking a spurious `Power::infy` when the
         * coefficient is later simplified. */
        arith_warnings_mute_push();
        Expr* at_x0 = eval_and_free(replace_all_of(current, ctx->x, ctx->x0));
        arith_warnings_mute_pop();
        if (has_infinity(at_x0)) {
            expr_free(at_x0); expr_free(current); expr_free(factorial);
            so_free(s); return NULL;
        }
        Expr* inv_fact = simp(mk_power(expr_copy(factorial), expr_new_integer(-1)));
        Expr* coef = simp(mk_times(at_x0, inv_fact));
        so_set_coef(s, (size_t)k, coef);
        if (k + 1 < n_iter) {
            Expr* next = simp(mk_fn2("D", current, expr_copy(ctx->x)));
            current = next;
        }
    }
    expr_free(current);
    expr_free(factorial);
    return s;
}
#undef MAX_NAIVE_ORDER

/* True iff `e` is an elementary head we have a direct kernel for. */
static bool is_known_elementary(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.arg_count != 1) return false;
    const char* names[] = {
        "Sin", "Cos", "Tan", "Exp", "Log",
        "Sinh", "Cosh", "Tanh",
        "ArcSin", "ArcCos", "ArcTan", "ArcCot",
        "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth",
        "SinIntegral", "SinhIntegral", "FresnelC", "FresnelS",
        NULL
    };
    for (int i = 0; names[i]; i++) if (has_symbol_head(e, names[i])) return true;
    return false;
}

/* Reciprocal rewrites. Covers:
 *   Forward-reciprocal trig/hyperbolic heads (Sec, Csc, Cot, Sech, Csch, Coth)
 *   which rewrite as 1/f[z] so the kernel path can handle them.
 *   Inverse-reciprocal heads (ArcSec, ArcCsc, ArcSech, ArcCsch) which rewrite
 *   via the identities
 *     ArcSec[z]  = ArcCos[1/z]
 *     ArcCsc[z]  = ArcSin[1/z]
 *     ArcSech[z] = ArcCosh[1/z]
 *     ArcCsch[z] = ArcSinh[1/z]
 *   so that composition with a blowing-up inner series (e.g. z = 1/x)
 *   collapses to a convergent kernel case rather than triggering naive
 *   Taylor's x=x0 probe (which would emit spurious `Power::infy: 1/0`
 *   warnings before the guard bails out).
 *
 * Returns a new expression (caller owns) or NULL if not applicable. */
static Expr* rewrite_reciprocal_head(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return NULL;
    if (e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol.name;
    Expr* arg = expr_copy(e->data.function.args[0]);
    if (h == SYM_Sec) return mk_power(mk_fn1("Cos",  arg), expr_new_integer(-1));
    if (h == SYM_Csc) return mk_power(mk_fn1("Sin",  arg), expr_new_integer(-1));
    if (h == SYM_Cot) return mk_times(mk_fn1("Cos",  expr_copy(arg)),
                                                mk_power(mk_fn1("Sin",  arg), expr_new_integer(-1)));
    if (h == SYM_Sinc) return mk_times(mk_fn1("Sin", expr_copy(arg)),
                                                mk_power(arg, expr_new_integer(-1)));
    if (h == SYM_Sech) return mk_power(mk_fn1("Cosh", arg), expr_new_integer(-1));
    if (h == SYM_Csch) return mk_power(mk_fn1("Sinh", arg), expr_new_integer(-1));
    if (h == SYM_Coth) return mk_times(mk_fn1("Cosh", expr_copy(arg)),
                                                mk_power(mk_fn1("Sinh", arg), expr_new_integer(-1)));
    /* Inverse-reciprocal identities. simp() collapses 1/(1/x) back to x, so
     * e.g. ArcSec[1/x] becomes ArcCos[x] and dispatches to the ArcCos kernel. */
    if (h == SYM_ArcSec) return mk_fn1("ArcCos",  simp(mk_power(arg, expr_new_integer(-1))));
    if (h == SYM_ArcCsc) return mk_fn1("ArcSin",  simp(mk_power(arg, expr_new_integer(-1))));
    if (h == SYM_ArcSech) return mk_fn1("ArcCosh", simp(mk_power(arg, expr_new_integer(-1))));
    if (h == SYM_ArcCsch) return mk_fn1("ArcSinh", simp(mk_power(arg, expr_new_integer(-1))));
    expr_free(arg);
    return NULL;
}

/* Taylor coefficient c_k of Zeta about s = 0 (k = 0..3), matching
 * Series[Zeta[x], {x, 0, 3}]. L = Log[2 Pi]. Returns an owned, evaluated Expr;
 * StieltjesGamma[k], Zeta[3], EulerGamma, Pi and L stay symbolic. */
static Expr* series_zeta0_coef(int64_t k) {
    Expr* L = mk_fn1("Log", mk_times(expr_new_integer(2), mk_symbol("Pi")));
    switch (k) {
    case 0:
        expr_free(L);
        return make_rational(-1, 2);
    case 1:
        /* -1/2 L */
        return simp(mk_times(make_rational(-1, 2), L));
    case 2: {
        /* 1/2 ( 1/2 EulerGamma^2 - 1/24 Pi^2 - 1/2 L^2 + StieltjesGamma[1] ) */
        Expr* t1 = mk_times(make_rational(1, 2),
                            mk_power(mk_symbol("EulerGamma"), expr_new_integer(2)));
        Expr* t2 = mk_times(make_rational(-1, 24),
                            mk_power(mk_symbol("Pi"), expr_new_integer(2)));
        Expr* t3 = mk_times(make_rational(-1, 2),
                            mk_power(expr_copy(L), expr_new_integer(2)));
        Expr* t4 = mk_fn1("StieltjesGamma", expr_new_integer(1));
        expr_free(L);
        Expr* inner = mk_plus(mk_plus(t1, t2), mk_plus(t3, t4));
        return simp(mk_times(make_rational(1, 2), inner));
    }
    case 3: {
        /* 1/48 ( 8 EG^3 + 12 EG^2 L - Pi^2 L - 4 L^3
         *        + 24 EG SG1 + 24 L SG1 + 12 SG2 - 8 Zeta[3] ),  EG = EulerGamma. */
        Expr* a1 = mk_times(expr_new_integer(8),
                            mk_power(mk_symbol("EulerGamma"), expr_new_integer(3)));
        Expr* a2 = mk_times(expr_new_integer(12),
                            mk_times(mk_power(mk_symbol("EulerGamma"), expr_new_integer(2)),
                                     expr_copy(L)));
        Expr* a3 = mk_times(expr_new_integer(-1),
                            mk_times(mk_power(mk_symbol("Pi"), expr_new_integer(2)),
                                     expr_copy(L)));
        Expr* a4 = mk_times(expr_new_integer(-4),
                            mk_power(expr_copy(L), expr_new_integer(3)));
        Expr* a5 = mk_times(expr_new_integer(24),
                            mk_times(mk_symbol("EulerGamma"),
                                     mk_fn1("StieltjesGamma", expr_new_integer(1))));
        Expr* a6 = mk_times(expr_new_integer(24),
                            mk_times(expr_copy(L),
                                     mk_fn1("StieltjesGamma", expr_new_integer(1))));
        Expr* a7 = mk_times(expr_new_integer(12),
                            mk_fn1("StieltjesGamma", expr_new_integer(2)));
        Expr* a8 = mk_times(expr_new_integer(-8),
                            mk_fn1("Zeta", expr_new_integer(3)));
        expr_free(L);
        Expr* sum = mk_plus(mk_plus(mk_plus(a1, a2), mk_plus(a3, a4)),
                            mk_plus(mk_plus(a5, a6), mk_plus(a7, a8)));
        return simp(mk_times(make_rational(1, 48), sum));
    }
    default:
        expr_free(L);
        return expr_new_integer(0);
    }
}

/* Custom series for Zeta[x] about x0 = 1 (the Laurent expansion that defines
 * the Stieltjes constants) and x0 = 0 (low-order Taylor). Returns NULL for any
 * other expansion point, or for x0 = 0 beyond x^3, deferring to the generic
 * differentiation path. */
static SeriesObj* series_zeta_at(SeriesCtx* ctx) {
    Expr* x  = ctx->x;
    Expr* x0 = ctx->x0;
    int64_t order = ctx->order;

    /* x0 = 1:  zeta(s) = 1/(s-1) + Sum_{k>=0} ((-1)^k / k!) gamma_k (s-1)^k,
     *          gamma_0 = EulerGamma, gamma_k = StieltjesGamma[k] (k >= 1). */
    if (x0->type == EXPR_INTEGER && x0->data.integer == 1) {
        SeriesObj* s = so_alloc(x, x0, -1, order, 1);
        for (size_t i = 0; i < s->coef_count; i++) {
            int64_t expo = s->nmin + (int64_t)i;     /* -1, 0, 1, ... */
            Expr* c;
            if (expo == -1) {
                c = expr_new_integer(1);
            } else if (expo == 0) {
                c = mk_symbol("EulerGamma");
            } else {
                /* (-1)^expo / expo! * StieltjesGamma[expo] (Factorial keeps
                 * the coefficient exact for any order via BigInt). */
                Expr* sg   = mk_fn1("StieltjesGamma", expr_new_integer(expo));
                Expr* finv = mk_power(mk_fn1("Factorial", expr_new_integer(expo)),
                                      expr_new_integer(-1));
                Expr* sgn  = expr_new_integer((expo % 2 == 0) ? 1 : -1);
                c = simp(mk_times(sgn, mk_times(finv, sg)));
            }
            so_set_coef(s, i, c);
        }
        return s;
    }

    /* x0 = 0:  low-order Taylor with closed-form coefficients (up to x^3).
     * The cap is on the user-facing target order; padding coefficients beyond
     * x^3 are zero and get truncated away by the outer Series builtin. */
    if (is_lit_zero(x0)) {
        /* target_order is the user n plus 1 (the O-term exponent); n <= 3. */
        if (ctx->target_order > 4) return NULL;      /* beyond x^3: defer */
        SeriesObj* s = so_alloc(x, x0, 0, order, 1);
        for (size_t i = 0; i < s->coef_count; i++) {
            so_set_coef(s, i, series_zeta0_coef((int64_t)i));
        }
        return s;
    }

    return NULL;
}

static SeriesObj* series_expand(Expr* e, SeriesCtx* ctx) {
    /* Early out: free of x => constant series. */
    if (expr_free_of(e, ctx->x)) {
        return so_from_constant(e, ctx->x, ctx->x0, ctx->order, 1);
    }
    /* e == x => identity series. */
    if (expr_eq(e, ctx->x)) {
        return so_from_variable(ctx->x, ctx->x0, ctx->order, 1);
    }

    /* e is itself a SeriesData about the same (x, x0): adopt it directly. This
     * lets expressions that embed a series flow through the expander -- notably
     * the Exp[exp*Log[base]] rewrite series_power uses for series exponents, and
     * Series[] of an expression that already contains a SeriesData. A series
     * about a different variable or point cannot be merged here, so we bail. */
    if (is_series_data(e)) {
        SeriesObj* s = so_from_seriesdata(e);
        if (!s) return NULL;
        if (!expr_eq(s->x, ctx->x) || !expr_eq(s->x0, ctx->x0)) {
            so_free(s);
            return NULL;
        }
        return s;
    }

    if (e->type == EXPR_FUNCTION) {
        const char* head = (e->data.function.head->type == EXPR_SYMBOL)
                               ? e->data.function.head->data.symbol.name : NULL;
        /* ---- Plus ---- */
        if (head && strcmp(head, "Plus") == 0) {
            SeriesObj* acc = NULL;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                SeriesObj* t = series_expand_nested(e->data.function.args[i], ctx);
                if (!t) { if (acc) so_free(acc); return NULL; }
                if (!acc) acc = t;
                else { SeriesObj* sum = so_add(acc, t); so_free(acc); so_free(t); acc = sum; }
            }
            if (!acc) acc = so_from_constant(expr_new_integer(0), ctx->x, ctx->x0, ctx->order, 1);
            return acc;
        }
        /* ---- Times ---- */
        if (head && strcmp(head, "Times") == 0) {
            SeriesObj* acc = NULL;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                SeriesObj* t = series_expand_nested(e->data.function.args[i], ctx);
                if (!t) { if (acc) so_free(acc); return NULL; }
                if (!acc) acc = t;
                else { SeriesObj* mul = so_mul(acc, t); so_free(acc); so_free(t); acc = mul; }
            }
            if (!acc) acc = so_from_constant(expr_new_integer(1), ctx->x, ctx->x0, ctx->order, 1);
            return acc;
        }
        /* ---- Power ---- */
        if (head && strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp  = e->data.function.args[1];
            bool base_has_x = !expr_free_of(base, ctx->x);
            bool exp_has_x  = !expr_free_of(exp, ctx->x);
            if (!base_has_x && !exp_has_x) {
                return so_from_constant(e, ctx->x, ctx->x0, ctx->order, 1);
            }
            if (base_has_x && !exp_has_x) {
                SeriesObj* bs = series_expand_nested(base, ctx);
                if (!bs) return NULL;
                SeriesObj* r = so_pow_expr(bs, exp);
                so_free(bs);
                return r;
            }
            /* Exponent depends on x: rewrite as Exp[exp * Log[base]] and recurse. */
            Expr* log_base = mk_fn1("Log", expr_copy(base));
            Expr* prod = mk_times(expr_copy(exp), log_base);
            Expr* rew = mk_fn1("Exp", prod);
            SeriesObj* r = series_expand_nested(rew, ctx);
            expr_free(rew);
            return r;
        }
        /* ---- Log[x]-at-x=0 symbolic coefficient ---- */
        if (head && strcmp(head, "Log") == 0 && e->data.function.arg_count == 1) {
            Expr* arg = e->data.function.args[0];
            /* Log[c*(x-x0)^r * (1+u)] = r*Log[x-x0] + Log[c] + Log[1+u].
             * For x0 = 0 and arg = x, this is Log[x] -- keep symbolic. */
            if (expr_eq(arg, ctx->x) && is_lit_zero(ctx->x0)) {
                return so_from_constant(e, ctx->x, ctx->x0, ctx->order, 1);
            }
            /* Maxima-style sp2log fast path: if arg = a + b*x^(p/q) with a, b
             * both free of x and both non-zero, rewrite as
             *   Log[a] + Log[1 + (b/a)*x^(p/q)]
             * so the Log1p kernel sees a pure monomial inner and the monomial
             * binomial fast path kicks in. Only applies at x0 = 0 because the
             * rewrite needs x -> x - x0 translation otherwise. */
            if (is_lit_zero(ctx->x0)) {
                Expr* a; Expr* b; int64_t en, ed;
                if (series_split_two_term(arg, ctx->x, &a, &b, &en, &ed)) {
                    /* Guard against self-recursion: when a is already 1, the
                     * rewrite Log[1] + Log[1 + b*x^c] reduces back to the
                     * original Log[1 + b*x^c] and we'd loop. Require a != 1. */
                    bool a_is_one = (a->type == EXPR_INTEGER && a->data.integer == 1);
                    if (!a_is_one && !is_lit_zero(a) && !is_lit_zero(b)) {
                        Expr* inv_a = simp(mk_power(expr_copy(a), expr_new_integer(-1)));
                        Expr* b_over_a = simp(mk_times(b, inv_a));
                        Expr* exp_expr = (ed == 1) ? expr_new_integer(en) : make_rational(en, ed);
                        Expr* x_pow = simp(mk_power(expr_copy(ctx->x), exp_expr));
                        Expr* monom = simp(mk_times(b_over_a, x_pow));
                        Expr* inner = simp(mk_plus(expr_new_integer(1), monom));
                        Expr* log1pu = mk_fn1("Log", inner);
                        Expr* log_a = simp(mk_fn1("Log", a));
                        Expr* rewrite = simp(mk_plus(log_a, log1pu));
                        SeriesObj* r = series_expand_nested(rewrite, ctx);
                        expr_free(rewrite);
                        if (r) return r;
                        /* Fall through to the generic path on failure. */
                    } else {
                        expr_free(a); expr_free(b);
                    }
                }
            }
            /* General: expand arg, split leading monomial if vanishing. */
            SeriesObj* as = series_expand_nested(arg, ctx);
            if (!as) return NULL;
            so_trim_leading(as);
            if (as->coef_count == 0) { so_free(as); return NULL; }
            int64_t nmin = as->nmin;
            if (nmin == 0) {
                SeriesObj* r = so_apply_log(as);
                so_free(as);
                return r;
            }
            /* nmin != 0: arg vanishes (or blows up) at x0. Peel off
             * (x-x0)^(nmin/den) and treat its log as a symbolic constant. */
            Expr* a0 = expr_copy(as->coefs[0]);
            Expr* inv_a0 = simp(mk_power(expr_copy(a0), expr_new_integer(-1)));
            SeriesObj* a_scaled = so_scalar_mul(inv_a0, as);
            expr_free(inv_a0);
            SeriesObj* shifted = so_shift_by_rational(a_scaled, -nmin, as->den);
            so_free(a_scaled);
            SeriesObj* log_rest = so_apply_log(shifted);
            so_free(shifted); so_free(as);
            if (!log_rest) { expr_free(a0); return NULL; }
            /* Add symbolic Log[a0] + (nmin/den) * Log[x - x0]. */
            Expr* log_a0 = simp(mk_fn1("Log", a0));
            Expr* base_sym;
            if (is_lit_zero(ctx->x0)) base_sym = expr_copy(ctx->x);
            else                      base_sym = simp(mk_plus(expr_copy(ctx->x),
                                                              simp(mk_times(expr_new_integer(-1),
                                                                            expr_copy(ctx->x0)))));
            Expr* log_base = simp(mk_fn1("Log", base_sym));
            Expr* n_over_d = make_rational(nmin, log_rest->den);
            Expr* extra = simp(mk_plus(log_a0, simp(mk_times(n_over_d, log_base))));
            SeriesObj* extras = so_from_constant(extra, log_rest->x, log_rest->x0, log_rest->order, log_rest->den);
            expr_free(extra);
            SeriesObj* r = so_add(extras, log_rest);
            so_free(extras); so_free(log_rest);
            return r;
        }
        /* ---- Zeta[x]: custom expansions about x0 = 1 (Laurent, introducing
         * the Stieltjes constants) and x0 = 0 (low-order Taylor). ---- */
        if (head && strcmp(head, "Zeta") == 0 && e->data.function.arg_count == 1 &&
            expr_eq(e->data.function.args[0], ctx->x)) {
            SeriesObj* z = series_zeta_at(ctx);
            if (z) return z;
            /* unsupported point/order: fall through to generic differentiation */
        }
        /* ---- Reciprocal heads (Sec, Csc, Cot, Sech, Csch, Coth) ---- */
        {
            Expr* rewrite = rewrite_reciprocal_head(e);
            if (rewrite) {
                SeriesObj* r = series_expand_nested(rewrite, ctx);
                expr_free(rewrite);
                return r;
            }
        }
        /* ---- Known elementary functions ---- */
        if (is_known_elementary(e)) {
            SeriesObj* inner = series_expand_nested(e->data.function.args[0], ctx);
            if (!inner) return NULL;
            SeriesObj* r = NULL;
            if (strcmp(head, "Exp") == 0)       r = so_apply_exp(inner);
            else if (strcmp(head, "Sin") == 0)  r = so_apply_sin_or_cos(inner, true);
            else if (strcmp(head, "Cos") == 0)  r = so_apply_sin_or_cos(inner, false);
            else if (strcmp(head, "Sinh") == 0) r = so_apply_sinh_or_cosh(inner, true);
            else if (strcmp(head, "Cosh") == 0) r = so_apply_sinh_or_cosh(inner, false);
            else if (strcmp(head, "Tan") == 0) {
                SeriesObj* sn = so_apply_sin_or_cos(inner, true);
                SeriesObj* cs = so_apply_sin_or_cos(inner, false);
                if (sn && cs) r = so_div(sn, cs);
                if (sn) so_free(sn);
                if (cs) so_free(cs);
            } else if (strcmp(head, "Tanh") == 0) {
                SeriesObj* sn = so_apply_sinh_or_cosh(inner, true);
                SeriesObj* cs = so_apply_sinh_or_cosh(inner, false);
                if (sn && cs) r = so_div(sn, cs);
                if (sn) so_free(sn);
                if (cs) so_free(cs);
            } else if (strcmp(head, "Log") == 0) {
                r = so_apply_log(inner);
            } else if (strcmp(head, "ArcTan")  == 0) {
                r = so_apply_kernel_at_zero("ArcTan", inner);
                if (!r) {
                    int sc = so_branch_point_imag_sign(inner);
                    if (sc != 0) r = so_apply_arctan_branch_point(inner, sc, ctx->target_order, ctx);
                }
            }
            else if (strcmp(head, "SinIntegral") == 0) {
                /* Si is entire and odd (Si(0)=0), analytic at u=0 like ArcTan. */
                r = so_apply_kernel_at_zero("SinIntegral", inner);
            }
            else if (strcmp(head, "SinhIntegral") == 0) {
                /* Shi is entire and odd (Shi(0)=0), analytic at u=0 like Si. */
                r = so_apply_kernel_at_zero("SinhIntegral", inner);
            }
            else if (strcmp(head, "FresnelC") == 0) {
                /* FresnelC is entire and odd (C(0)=0), analytic at u=0. */
                r = so_apply_kernel_at_zero("FresnelC", inner);
            }
            else if (strcmp(head, "FresnelS") == 0) {
                /* FresnelS is entire and odd (S(0)=0), analytic at u=0. */
                r = so_apply_kernel_at_zero("FresnelS", inner);
            }
            else if (strcmp(head, "ArcTanh") == 0) {
                r = so_apply_kernel_at_zero("ArcTanh", inner);
                /* Branch points at ±1 first; then the at-infinity rewrite. */
                if (!r) {
                    int sc = so_branch_point_sign(inner);
                    if (sc != 0) r = so_apply_arctanh_branch_point(inner, sc, ctx->target_order, ctx);
                }
                /* If arg blows up at x0, use principal-branch identity
                 * ArcTanh[1/u] = I*Pi/2 + ArcTanh[u]. */
                if (!r && inner->nmin < 0) {
                    SeriesObj* inv = so_inv(inner);
                    if (inv) {
                        SeriesObj* at = so_apply_kernel_at_zero("ArcTanh", inv);
                        so_free(inv);
                        if (at) {
                            SeriesObj* k = so_const_half_i_pi(at);
                            r = so_add(k, at);
                            so_free(k); so_free(at);
                        }
                    }
                }
            }
            else if (strcmp(head, "ArcSin")  == 0) {
                r = so_apply_kernel_at_zero("ArcSin", inner);
                if (!r) {
                    int sc = so_branch_point_sign(inner);
                    if (sc != 0) r = so_apply_arc_branch_point(inner, sc, true, ctx->target_order, ctx);
                }
            }
            else if (strcmp(head, "ArcSinh") == 0) {
                r = so_apply_kernel_at_zero("ArcSinh", inner);
                /* Branch points at ±I (new): square-root Puiseux series. */
                if (!r) {
                    int sc = so_branch_point_imag_sign(inner);
                    if (sc != 0) r = so_apply_arcsinh_branch_point(inner, sc, ctx->target_order, ctx);
                }
                /* Identity at infinity: ArcSinh[1/v] = -Log[v] + Log[1 + Sqrt[1 + v^2]].
                 * Fire when inner actually blows up (nmin < 0); series_expand's
                 * Log-at-x branch absorbs the symbolic -Log[v]. */
                if (!r && inner->nmin < 0) {
                    Expr* arg = e->data.function.args[0];
                    Expr* v = simp(mk_power(expr_copy(arg), expr_new_integer(-1)));
                    Expr* v_sq = simp(mk_power(expr_copy(v), expr_new_integer(2)));
                    Expr* sqrt_term = simp(mk_power(
                        simp(mk_plus(expr_new_integer(1), v_sq)),
                        make_rational(1, 2)));
                    Expr* log2 = simp(mk_fn1("Log",
                                    simp(mk_plus(expr_new_integer(1), sqrt_term))));
                    Expr* neg_log_v = simp(mk_times(expr_new_integer(-1),
                                       simp(mk_fn1("Log", v))));
                    Expr* rewrite = simp(mk_plus(neg_log_v, log2));
                    r = series_expand_nested(rewrite, ctx);
                    expr_free(rewrite);
                }
            }
            else if (strcmp(head, "ArcCos")  == 0) {
                r = so_apply_arccos(inner);
                if (!r) {
                    int sc = so_branch_point_sign(inner);
                    if (sc != 0) r = so_apply_arc_branch_point(inner, sc, false, ctx->target_order, ctx);
                }
            }
            else if (strcmp(head, "ArcCosh") == 0) {
                r = so_apply_arccosh(inner);
                /* Branch points at ±1 (new): square-root Puiseux series. */
                if (!r) {
                    int sc = so_branch_point_sign(inner);
                    if (sc != 0) r = so_apply_arccosh_branch_point(inner, sc, ctx->target_order, ctx);
                }
                /* Identity at infinity: ArcCosh[1/v] = -Log[v] + Log[1 + Sqrt[1 - v^2]]. */
                if (!r && inner->nmin < 0) {
                    Expr* arg = e->data.function.args[0];
                    Expr* v = simp(mk_power(expr_copy(arg), expr_new_integer(-1)));
                    Expr* v_sq = simp(mk_power(expr_copy(v), expr_new_integer(2)));
                    Expr* one_minus_v2 = simp(mk_plus(expr_new_integer(1),
                                          simp(mk_times(expr_new_integer(-1), v_sq))));
                    Expr* sqrt_term = simp(mk_power(one_minus_v2, make_rational(1, 2)));
                    Expr* log2 = simp(mk_fn1("Log",
                                    simp(mk_plus(expr_new_integer(1), sqrt_term))));
                    Expr* neg_log_v = simp(mk_times(expr_new_integer(-1),
                                       simp(mk_fn1("Log", v))));
                    Expr* rewrite = simp(mk_plus(neg_log_v, log2));
                    r = series_expand_nested(rewrite, ctx);
                    expr_free(rewrite);
                }
            }
            else if (strcmp(head, "ArcCot")  == 0) {
                r = so_apply_arccot(inner);
                /* Branch points at ±I: logarithmic series. */
                if (!r) {
                    int sc = so_branch_point_imag_sign(inner);
                    if (sc != 0) r = so_apply_arccot_branch_point(inner, sc, ctx->target_order, ctx);
                }
                /* If arg blows up at x0 (e.g. 1/x), fall back to the
                 * at-infinity branch: ArcCot[u] = ArcTan[1/u]. */
                if (!r && inner->nmin < 0) {
                    SeriesObj* inv = so_inv(inner);
                    if (inv) { r = so_apply_kernel_at_zero("ArcTan", inv); so_free(inv); }
                }
            }
            else if (strcmp(head, "ArcCoth") == 0) {
                r = so_apply_arccoth(inner);
                /* Branch points at ±1: logarithmic series. */
                if (!r) {
                    int sc = so_branch_point_sign(inner);
                    if (sc != 0) r = so_apply_arccoth_branch_point(inner, sc, ctx->target_order, ctx);
                }
                /* If arg blows up at x0 (e.g. 1/x), use ArcCoth[u] = ArcTanh[1/u]. */
                if (!r && inner->nmin < 0) {
                    SeriesObj* inv = so_inv(inner);
                    if (inv) { r = so_apply_kernel_at_zero("ArcTanh", inv); so_free(inv); }
                }
            }
            so_free(inner);
            /* Fall back to naive Taylor-via-D when the kernel path can't
             * handle the expansion point -- e.g. ArcSin[x] at x = 1/2, where
             * the arg's constant term is non-zero so the at-zero kernels
             * return NULL. D-path still fails at true branch points
             * (ArcSin[x] at x = 1) because the derivative blows up there,
             * which Puiseux expansion would handle but is out of scope. */
            if (!r) r = series_taylor_via_D(e, ctx);
            return r;
        }
    }

    /* Fallback: naive Taylor via D. */
    return series_taylor_via_D(e, ctx);
}

/* ----------------------------------------------------------------------------
 * Normal
 * -------------------------------------------------------------------------- */

static Expr* series_build_xmx0(Expr* x, Expr* x0) {
    if ((x0->type == EXPR_INTEGER && x0->data.integer == 0) ||
        (x0->type == EXPR_REAL && x0->data.real == 0.0)) {
        return expr_copy(x);
    }
    Expr* neg = simp(mk_times(expr_new_integer(-1), expr_copy(x0)));
    return simp(mk_plus(expr_copy(x), neg));
}

/* Convert a single validated SeriesData[x, x0, coefs, nmin, nmax, den] node
 * into an ordinary sum, dropping the O-term. Returns NULL if the node is
 * malformed (so the caller can leave it untouched). */
/* ----------------------------------------------------------------------------
 * Calculus on SeriesData (term-by-term differentiation / integration)
 * -------------------------------------------------------------------------- */

/* Read a validated SeriesData[x, x0, {coefs}, nmin, nmax, den] Expr into a
 * fresh SeriesObj. Returns NULL if the shape is unsupported (non-integer
 * nmin/nmax/den, missing List of coefficients, zero denominator); the caller
 * then leaves the enclosing expression unevaluated. */
static SeriesObj* so_from_seriesdata(Expr* sd) {
    if (sd->type != EXPR_FUNCTION ||
        !has_symbol_head(sd, "SeriesData") ||
        sd->data.function.arg_count != 6) return NULL;
    Expr** a = sd->data.function.args;
    Expr* coefs = a[2];
    if (!has_symbol_head(coefs, "List") ||
        a[3]->type != EXPR_INTEGER ||
        a[4]->type != EXPR_INTEGER ||
        a[5]->type != EXPR_INTEGER ||
        a[5]->data.integer <= 0) return NULL;
    int64_t nmin  = a[3]->data.integer;
    int64_t order = a[4]->data.integer;
    int64_t den   = a[5]->data.integer;
    size_t  k     = coefs->data.function.arg_count;
    /* The coefficient list length must agree with order - nmin. */
    if ((int64_t)k != order - nmin) return NULL;
    SeriesObj* s = so_alloc(a[0], a[1], nmin, order, den);
    for (size_t i = 0; i < k && i < s->coef_count; i++) {
        so_set_coef(s, i, expr_copy(coefs->data.function.args[i]));
    }
    return s;
}

/* D[SeriesData[...], var] -- term-by-term differentiation. When var is the
 * series variable, apply the power rule to each term; otherwise (var free of
 * the expansion point) differentiate each coefficient, keeping the powers.
 * Returns NULL when the shape is unsupported, so the caller falls back. */
Expr* series_differentiate(Expr* sd, Expr* var) {
    SeriesObj* s = so_from_seriesdata(sd);
    if (!s) return NULL;

    /* Differentiating with respect to a variable other than the series
     * variable: thread D into the coefficients (powers unchanged). Only valid
     * when the expansion point does not itself depend on var. */
    if (!expr_eq(s->x, var)) {
        if (!expr_free_of(s->x0, var)) { so_free(s); return NULL; }
        SeriesObj* r = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) {
            so_set_coef(r, i, simp(mk_fn2("D", expr_copy(s->coefs[i]), expr_copy(var))));
        }
        so_trim_leading(r);
        Expr* out = so_to_expr(r);
        so_free(s); so_free(r);
        return out;
    }

    /* Power rule: a_i (x-x0)^((nmin+i)/den)
     *   ->  a_i * (nmin+i)/den * (x-x0)^((nmin+i-den)/den). */
    SeriesObj* r = so_alloc(s->x, s->x0, s->nmin - s->den, s->order - s->den, s->den);
    for (size_t i = 0; i < s->coef_count && i < r->coef_count; i++) {
        int64_t num = s->nmin + (int64_t)i;          /* old exponent numerator */
        Expr* factor = make_rational(num, s->den);   /* (nmin+i)/den */
        so_set_coef(r, i, simp(mk_times(expr_copy(s->coefs[i]), factor)));
    }
    so_trim_leading(r);
    Expr* out = so_to_expr(r);
    so_free(s); so_free(r);
    return out;
}

/* Integrate[SeriesData[...], var] -- term-by-term integration. When var is the
 * series variable, integrate each term (constant of integration 0, as in
 * Mathematica); otherwise (var free of the expansion point) integrate each
 * coefficient. Returns NULL when unsupported, including the residue case where
 * a nonzero coefficient sits at exponent -1 (whose integral is a Log term that
 * SeriesData cannot represent). */
Expr* series_integrate(Expr* sd, Expr* var) {
    SeriesObj* s = so_from_seriesdata(sd);
    if (!s) return NULL;

    /* Integration with respect to a variable other than the series variable:
     * thread Integrate into the coefficients (powers unchanged). */
    if (!expr_eq(s->x, var)) {
        if (!expr_free_of(s->x0, var)) { so_free(s); return NULL; }
        SeriesObj* r = so_alloc(s->x, s->x0, s->nmin, s->order, s->den);
        for (size_t i = 0; i < s->coef_count; i++) {
            so_set_coef(r, i, simp(mk_fn2("Integrate", expr_copy(s->coefs[i]), expr_copy(var))));
        }
        Expr* out = so_to_expr(r);
        so_free(s); so_free(r);
        return out;
    }

    /* Residue guard: a nonzero (x-x0)^-1 term would integrate to a Log. */
    for (size_t i = 0; i < s->coef_count; i++) {
        if (s->nmin + (int64_t)i + s->den == 0 && !is_lit_zero(s->coefs[i])) {
            so_free(s);
            return NULL;
        }
    }

    /* Power rule: a_i (x-x0)^((nmin+i)/den)
     *   ->  a_i * den/(nmin+i+den) * (x-x0)^((nmin+i+den)/den). */
    SeriesObj* r = so_alloc(s->x, s->x0, s->nmin + s->den, s->order + s->den, s->den);
    for (size_t i = 0; i < s->coef_count && i < r->coef_count; i++) {
        int64_t denom = s->nmin + (int64_t)i + s->den;  /* new exponent numerator */
        if (denom == 0) {
            /* The (x-x0)^-1 slot: its coefficient is zero (guaranteed by the
             * residue guard above), so the integrated coefficient is zero. */
            so_set_coef(r, i, expr_new_integer(0));
            continue;
        }
        Expr* factor = make_rational(s->den, denom);    /* den/(nmin+i+den) */
        so_set_coef(r, i, simp(mk_times(expr_copy(s->coefs[i]), factor)));
    }
    so_trim_leading(r);
    Expr* out = so_to_expr(r);
    so_free(s); so_free(r);
    return out;
}

/* ----------------------------------------------------------------------------
 * SeriesData arithmetic (Plus / Times / Power dispatch targets)
 *
 * The strategy is uniform: convert every operand to a SeriesObj about the
 * common (x, x0), then drive the existing so_add / so_mul / so_pow_* algebra.
 * A non-series operand is expanded with the full series engine, so constants,
 * the bare variable, polynomials, and transcendental functions all combine the
 * way Mathematica does. Incompatible operands (different x or x0, or something
 * the engine cannot expand) yield NULL so the caller leaves the node alone.
 * -------------------------------------------------------------------------- */

/* Convert a single Plus/Times operand to a SeriesObj about (x, x0), expanding
 * to O-term exponent order_num/den. On an incompatible operand (a SeriesData
 * about a different point, or an expression the engine cannot expand) sets
 * *incompatible and returns NULL. */
/* True iff `e` depends on the series variable `x` through a piecewise-constant
 * / non-analytic head (Floor, Arg, Sign, ...). Such a factor has no power-series
 * expansion about a generic point, and Mathematica keeps it as a symbolic
 * coefficient OUTSIDE the SeriesData rather than folding it in -- e.g. the
 * branch-point wrappers emit Times[(-1)^Floor[(Pi/2 - Arg[x])/(2 Pi)], series].
 * series_combine must detect these and bail to a symbolic Plus/Times: feeding
 * one to series_expand below otherwise loops forever (it never reduces). */
static bool expr_nonanalytic_in(const Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (expr_free_of((Expr*)e, x)) return false;   /* no dependence -> safe */
    Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name) {
        static const char* const bad[] = {
            "Floor", "Ceiling", "Round", "IntegerPart", "FractionalPart",
            "Mod", "Quotient", "Sign", "Abs", "Arg", "UnitStep",
            "KroneckerDelta", "Boole", "Max", "Min", NULL };
        for (size_t i = 0; bad[i]; i++)
            if (strcmp(h->data.symbol.name, bad[i]) == 0) return true;
    }
    if (expr_nonanalytic_in(h, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_nonanalytic_in(e->data.function.args[i], x)) return true;
    return false;
}

/* The underlying expansion variable behind a (possibly compound) formal series
 * variable. Series about a finite point uses an atomic symbol; series at
 * Infinity uses Power[var, -k] (e.g. Power[x,-1] = 1/x). Returns the base symbol
 * (`x`) so callers can tell a genuine constant (free of `x`) from a factor that
 * is structurally free of the formal variable 1/x yet still depends on x -- such
 * as E^x in Times[E^x, SeriesData[1/x, ...]], which must NOT fold into the
 * coefficients. Returns NULL for an unrecognised shape (treated as no base). */
static Expr* series_base_var(Expr* x) {
    if (!x) return NULL;
    if (x->type == EXPR_SYMBOL) return x;
    if (x->type == EXPR_FUNCTION &&
        x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol.name == SYM_Power &&
        x->data.function.arg_count == 2 &&
        x->data.function.args[0]->type == EXPR_SYMBOL)
        return x->data.function.args[0];
    return NULL;
}

static SeriesObj* operand_to_series(Expr* op, Expr* x, Expr* x0,
                                    int64_t order_num, int64_t den,
                                    bool* incompatible) {
    SeriesObj* s = NULL;
    if (is_series_data(op)) {
        s = so_from_seriesdata(op);
        if (!s) { *incompatible = true; return NULL; }
        if (!expr_eq(s->x, x) || !expr_eq(s->x0, x0)) {
            so_free(s);
            *incompatible = true;
            return NULL;
        }
        if (s->den != den) {
            SeriesObj* r = so_rescale(s, lcm_i64(s->den, den));
            so_free(s);
            s = r;
        }
    } else if (expr_free_of(op, x)) {
        /* Free of the formal series variable. For a compound variable (series at
         * Infinity, where x is 1/var), an operand that still depends on the
         * underlying var is not a genuine constant and must not fold into the
         * coefficients -- e.g. Times[E^x, SeriesData[1/x, ...]] stays symbolic. */
        Expr* bv = series_base_var(x);
        if (bv && bv != x && !expr_free_of(op, bv)) {
            *incompatible = true;
            return NULL;
        }
        /* Constant: folds into the a0 term. Fast path for the common
         * `series + scalar` / `scalar * series` shapes. */
        s = so_from_constant(op, x, x0, order_num, den);
    } else if (expr_nonanalytic_in(op, x)) {
        /* Piecewise-constant / non-analytic dependence on x (Floor, Arg, ...):
         * no power-series expansion exists. Bail so the node stays symbolic,
         * matching Mathematica's branch-point wrapper output. */
        *incompatible = true;
        return NULL;
    } else {
        /* Anything else: expand it as a series about (x, x0). Pad the internal
         * order the same way builtin_series does so Laurent/Puiseux operands
         * survive. */
        /* A concrete number -- including Complex[num, num] (every UHP pole) and
         * Rational -- keeps so_inv's convolution coefficients numeric and
         * bounded, so it earns the full Laurent pad. Only a free/symbolic x0
         * risks the O(N^2) symbolic-coefficient blow-up the tight pad guards. */
        bool x0_is_numeric = expr_is_numeric_like(x0);
        int64_t target = order_num > 0 ? order_num : 1;
        int64_t pad = x0_is_numeric ? 12 : 2;
        SeriesCtx ctx = {
            x, x0, target + pad, target,
            /* allow_branch_wrap */ true,
            /* pending_add_const  */ NULL,
            /* pending_log_coef   */ NULL,
            /* pending_discriminator */ NULL,
        };
        s = series_expand(op, &ctx);
        if (!s) { *incompatible = true; return NULL; }
    }

    /* Trim leading zeros so nmin reflects the true leading exponent. so_mul's
     * order tracking (min(a.order+b.nmin, b.order+a.nmin)) depends on accurate
     * nmin: an untrimmed identity series for `x` (nmin 0, coef0 = 0) would
     * otherwise pull the product's O-term one power too low. */
    so_trim_leading(s);
    return s;
}

/* Shared core for Plus and Times: gather the SeriesData operands to find the
 * controlling (x, x0, den) and minimum order, convert every operand, then fold
 * with `combine` (so_add or so_mul). */
static Expr* series_combine(Expr* const* args, size_t n,
                            SeriesObj* (*combine)(SeriesObj*, SeriesObj*)) {
    /* The controlling (x, x0) come from the first SeriesData operand; borrow
     * them straight out of the live argument tree (no copy) so they stay valid
     * for the whole call. The common denominator is the lcm of all SeriesData
     * denominators, and the controlling O-term order is the minimum order over
     * the SeriesData operands, rescaled to that denominator. */
    Expr* x = NULL;
    Expr* x0 = NULL;
    int64_t den = 1;
    int64_t order_num = INT64_MAX;
    for (size_t i = 0; i < n; i++) {
        if (!is_series_data(args[i])) continue;
        Expr** sa = args[i]->data.function.args;
        if (sa[3]->type != EXPR_INTEGER || sa[4]->type != EXPR_INTEGER ||
            sa[5]->type != EXPR_INTEGER || sa[5]->data.integer <= 0) return NULL;
        if (!x) {
            x = sa[0]; x0 = sa[1]; den = sa[5]->data.integer;
        } else {
            if (!expr_eq(sa[0], x) || !expr_eq(sa[1], x0)) return NULL;
            den = lcm_i64(den, sa[5]->data.integer);
        }
    }
    if (!x) return NULL;                      /* no (well-formed) SeriesData */

    /* Second pass: now that den is known, fold in each operand's order. */
    for (size_t i = 0; i < n; i++) {
        if (!is_series_data(args[i])) continue;
        Expr** sa = args[i]->data.function.args;
        int64_t o = sa[4]->data.integer * (den / sa[5]->data.integer);
        if (o < order_num) order_num = o;
    }
    if (order_num == INT64_MAX) order_num = 1;

    bool incompatible = false;
    SeriesObj* acc = NULL;
    for (size_t i = 0; i < n; i++) {
        SeriesObj* t = operand_to_series(args[i], x, x0, order_num, den, &incompatible);
        if (!t) {
            if (acc) so_free(acc);
            return NULL;                     /* incompatible -> stays symbolic */
        }
        if (!acc) {
            acc = t;
        } else {
            SeriesObj* r = combine(acc, t);
            so_free(acc); so_free(t);
            acc = r;
        }
    }
    if (!acc) return NULL;
    so_trim_leading(acc);
    Expr* out = so_to_expr(acc);
    so_free(acc);
    return out;
}

Expr* series_combine_plus(Expr* const* args, size_t n) {
    return series_combine(args, n, so_add);
}

Expr* series_combine_times(Expr* const* args, size_t n) {
    return series_combine(args, n, so_mul);
}

Expr* series_power(Expr* base, Expr* exp) {
    bool base_series = is_series_data(base);
    bool exp_series  = is_series_data(exp);
    if (!base_series && !exp_series) return NULL;

    /* Direct path: SeriesData base raised to a scalar (free of the series
     * variable) exponent. Integer exponents go through so_pow_int (which
     * handles negative powers via so_inv); other scalars through so_pow_expr. */
    if (base_series) {
        SeriesObj* b = so_from_seriesdata(base);
        if (!b) return NULL;
        if (expr_free_of(exp, b->x)) {
            SeriesObj* r = so_pow_expr(b, exp);
            so_free(b);
            if (!r) return NULL;
            so_trim_leading(r);
            Expr* out = so_to_expr(r);
            so_free(r);
            return out;
        }
        so_free(b);
    }

    /* General path: the exponent depends on the series variable, or the base is
     * an ordinary expression raised to a series exponent (e.g. 2^series).
     * Rewrite base^exp = Exp[exp * Log[base]] and re-expand. The controlling
     * (x, x0, order, den) comes from whichever operand is the series. */
    Expr* sd = base_series ? base : exp;
    Expr* x   = sd->data.function.args[0];
    Expr* x0  = sd->data.function.args[1];
    SeriesObj* probe = so_from_seriesdata(sd);
    if (!probe) return NULL;
    int64_t order_num = probe->order;
    int64_t den = probe->den;
    so_free(probe);

    Expr* rewrite = mk_fn1("Exp",
        mk_times(expr_copy(exp), mk_fn1("Log", expr_copy(base))));
    bool incompatible = false;
    SeriesObj* r = operand_to_series(rewrite, x, x0, order_num, den, &incompatible);
    expr_free(rewrite);
    if (!r) return NULL;
    so_trim_leading(r);
    Expr* out = so_to_expr(r);
    so_free(r);
    return out;
}

static Expr* seriesdata_to_normal(Expr* arg) {
    Expr* x      = arg->data.function.args[0];
    Expr* x0     = arg->data.function.args[1];
    Expr* coefs  = arg->data.function.args[2];
    Expr* nmin_e = arg->data.function.args[3];
    Expr* den_e  = arg->data.function.args[5];

    if (!has_symbol_head(coefs, "List") ||
        nmin_e->type != EXPR_INTEGER ||
        den_e->type  != EXPR_INTEGER ||
        den_e->data.integer == 0) return NULL;

    int64_t nmin = nmin_e->data.integer;
    int64_t den  = den_e->data.integer;
    size_t k     = coefs->data.function.arg_count;

    Expr* base = series_build_xmx0(x, x0);
    Expr** terms = calloc(k + 1, sizeof(Expr*));
    size_t tc = 0;
    for (size_t i = 0; i < k; i++) {
        Expr* coef = coefs->data.function.args[i];
        if (is_lit_zero(coef)) continue;
        int64_t num = nmin + (int64_t)i;
        Expr* exp_e = (den == 1) ? expr_new_integer(num) : make_rational(num, den);
        Expr* term;
        if (is_lit_zero(exp_e)) { expr_free(exp_e); term = expr_copy(coef); }
        else {
            Expr* p = simp(mk_power(expr_copy(base), exp_e));
            term = simp(mk_times(expr_copy(coef), p));
        }
        terms[tc++] = term;
    }
    expr_free(base);

    Expr* sum;
    if (tc == 0)      sum = expr_new_integer(0);
    else if (tc == 1) sum = terms[0];
    else {
        sum = expr_new_function(mk_symbol("Plus"), terms, tc);
        terms = NULL;
    }
    if (terms) free(terms);

    return simp(sum);
}

/* Recursively replace every SeriesData subexpression (at any depth) with its
 * O-term-free sum. Series expansions around +-Infinity wrap their SeriesData
 * inside Plus/Times (e.g. the trig-prefactored asymptotic forms of BesselJ,
 * AiryAi, ...), so Normal must see through the wrapping rather than only
 * handling a top-level SeriesData. Returns a freshly allocated tree. */
static Expr* normal_recurse(Expr* e) {
    if (e->type == EXPR_FUNCTION &&
        has_symbol_head(e, "SeriesData") &&
        e->data.function.arg_count == 6) {
        Expr* n = seriesdata_to_normal(e);
        if (n) return n;
        /* Malformed SeriesData: fall through to a structural copy. */
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    Expr* head = normal_recurse(e->data.function.head);
    size_t n = e->data.function.arg_count;
    Expr** args = n ? calloc(n, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) {
        args[i] = normal_recurse(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(head, args, n);
    if (args) free(args);
    return out;
}

Expr* builtin_normal(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Normal[ndarray] converts an NDArray back to its equivalent nested List. */
    if (arg->type == EXPR_NDARRAY) {
        return ndarray_to_nested_list(arg);
    }

    /* Normal[assoc] converts an association to its list of rules. */
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_Association) {
        size_t n = arg->data.function.arg_count;
        Expr** rules = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) rules[i] = expr_copy(arg->data.function.args[i]);
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
        free(rules);
        return list;
    }

    /* Drop the O-term from every SeriesData, including those nested inside the
     * Plus/Times wrappers produced by expansions at +-Infinity. The rebuilt
     * tree is re-evaluated by the evaluator, which recombines the factors. */
    return normal_recurse(arg);
}

/* ----------------------------------------------------------------------------
 * Series builtin
 * -------------------------------------------------------------------------- */

/* Sign of a recognizable real numeric literal: -1, 0, +1; or 2 ("unknown")
 * for anything that is not an obvious real number. Used to interpret the
 * Assumptions option of Series (e.g. the bound in `x < 0`). */
static int lit_real_sign(Expr* e) {
    if (!e) return 2;
    switch (e->type) {
        case EXPR_INTEGER:
            return e->data.integer < 0 ? -1 : (e->data.integer > 0 ? 1 : 0);
        case EXPR_REAL:
            return e->data.real < 0.0 ? -1 : (e->data.real > 0.0 ? 1 : 0);
        case EXPR_BIGINT:
            return mpz_sgn(e->data.bigint);
#ifdef USE_MPFR
        case EXPR_MPFR:
            return mpfr_sgn(e->data.mpfr) < 0 ? -1
                 : (mpfr_sgn(e->data.mpfr) > 0 ? 1 : 0);
#endif
        default: break;
    }
    /* Rational[p, q] has q > 0 by construction, so its sign is sign of p. */
    if (has_symbol_head(e, "Rational") && e->data.function.arg_count == 2)
        return lit_real_sign(e->data.function.args[0]);
    return 2;
}

/* Interpret a Series Assumptions option as the sign of the expansion
 * variable `x`. Returns -1 if the assumption forces x < 0, +1 if it forces
 * x > 0, and 0 if it says nothing definite about x's sign. Recognizes the
 * ordering of x against a real numeric bound (Less/LessEqual/Greater/
 * GreaterEqual, in either argument order) and conjunctions (And) thereof.
 * Series uses this to choose the principal branch of Log[x] in the
 * logarithmic expansions of ExpIntegralEi / LogIntegral at x = 0. */
static int assumption_sign_of(Expr* assm, Expr* x) {
    if (!assm || assm->type != EXPR_FUNCTION) return 0;
    size_t ac = assm->data.function.arg_count;
    /* And[...]: the first component with a definite sign wins. */
    if (has_symbol_head(assm, "And")) {
        for (size_t i = 0; i < ac; i++) {
            int s = assumption_sign_of(assm->data.function.args[i], x);
            if (s) return s;
        }
        return 0;
    }
    bool less    = has_symbol_head(assm, "Less")    || has_symbol_head(assm, "LessEqual");
    bool greater = has_symbol_head(assm, "Greater") || has_symbol_head(assm, "GreaterEqual");
    if ((!less && !greater) || ac != 2) return 0;

    Expr* a = assm->data.function.args[0];
    Expr* b = assm->data.function.args[1];
    Expr* bound; bool x_on_left;
    if (expr_eq(a, x))      { bound = b; x_on_left = true;  }
    else if (expr_eq(b, x)) { bound = a; x_on_left = false; }
    else return 0;

    int bs = lit_real_sign(bound);
    if (bs == 2) return 0;
    /* Normalize "a REL b" with x on one side into "x < bound" vs "x > bound". */
    bool x_less_than_bound = (x_on_left && less) || (!x_on_left && greater);
    if (x_less_than_bound) return (bs <= 0) ? -1 : 0;  /* x < (<=0) => x < 0 */
    else                   return (bs >= 0) ?  1 : 0;  /* x > (>=0) => x > 0 */
}

/* Parse a single spec argument, accepting either `{x, x0, n}` (full form)
 * or `x -> x0` (leading-term form). Returns true on success and populates
 * *x_out / *x0_out / *n_out (borrowed Expr pointers into the spec) and
 * sets *leading_only accordingly. */
static bool parse_series_spec(Expr* spec, Expr** x_out, Expr** x0_out,
                              int64_t* n_out, bool* leading_only) {
    if (has_symbol_head(spec, "List") && spec->data.function.arg_count == 3) {
        Expr* n_e = spec->data.function.args[2];
        if (n_e->type != EXPR_INTEGER) {
            /* Evaluate in case the user passed a computable expression. */
            Expr* ev = eval_and_free(expr_copy(n_e));
            if (ev->type != EXPR_INTEGER) { expr_free(ev); return false; }
            /* Replace into the spec's own slot for the caller to read back;
             * but we just populate outputs and leak the eval? Easier: we
             * return the integer value directly. */
            *x_out         = spec->data.function.args[0];
            *x0_out        = spec->data.function.args[1];
            *n_out         = ev->data.integer;
            *leading_only  = false;
            expr_free(ev);
            return true;
        }
        *x_out         = spec->data.function.args[0];
        *x0_out        = spec->data.function.args[1];
        *n_out         = n_e->data.integer;
        *leading_only  = false;
        return true;
    }
    if (has_symbol_head(spec, "Rule") && spec->data.function.arg_count == 2) {
        *x_out         = spec->data.function.args[0];
        *x0_out        = spec->data.function.args[1];
        *n_out         = 0;
        *leading_only  = true;
        return true;
    }
    return false;
}

/* Detect a single Times-factor of the form Power[x, alpha] where alpha is
 * free of x and is NOT an integer/rational (so the series machinery can't
 * absorb it into the exponent of the monomial). If found, returns a
 * freshly allocated Power[x, alpha] and sets *body_out to the product of
 * the remaining factors (also freshly allocated). Returns NULL otherwise.
 *
 * Mathematica prints Series[x^a Exp[x], {x, 0, 5}] as
 *   x^a (1 + x + x^2/2 + ... + O[x]^6)
 * i.e. the symbolic power is factored out verbatim. We mimic this by
 * expanding `body_out` as an ordinary series and wrapping the SeriesData
 * result in a Times with the pre-factor at the end. */
static Expr* try_factor_power_prefactor(Expr* f, Expr* x, Expr** body_out) {
    if (f->type != EXPR_FUNCTION || !has_symbol_head(f, "Times")) return NULL;
    size_t n = f->data.function.arg_count;
    Expr* prefactor = NULL;
    size_t found_idx = n;
    for (size_t i = 0; i < n; i++) {
        Expr* a = f->data.function.args[i];
        if (!has_symbol_head(a, "Power") || a->data.function.arg_count != 2) continue;
        Expr* base = a->data.function.args[0];
        Expr* exp  = a->data.function.args[1];
        if (!expr_eq(base, x)) continue;
        if (!expr_free_of(exp, x)) continue;
        if (exp->type == EXPR_INTEGER) continue;
        int64_t pp, qq;
        if (is_rational(exp, &pp, &qq)) continue;
        prefactor = expr_copy(a);
        found_idx = i;
        break;
    }
    if (!prefactor) return NULL;
    /* Build the rest of the factors as a new Times (or a single Expr). */
    size_t rc = n > 0 ? n - 1 : 0;
    Expr** rest = rc ? calloc(rc, sizeof(Expr*)) : NULL;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == found_idx) continue;
        rest[k++] = expr_copy(f->data.function.args[i]);
    }
    Expr* body;
    if (rc == 0)      { body = expr_new_integer(1); if (rest) free(rest); }
    else if (rc == 1) { body = rest[0]; free(rest); }
    else {
        /* expr_new_function copies the pointers out of `rest` into its own
         * args buffer; we must still free the caller-owned array. */
        body = expr_new_function(mk_symbol("Times"), rest, rc);
        free(rest);
    }
    *body_out = body;
    return prefactor;
}

/* True iff `e` contains a Power[g(x), k] subterm where g involves x and
 * k is a negative integer. A cheap conservative test used to decide whether
 * it's worth calling Apart on the input. */
static bool has_negative_power_in(Expr* e, Expr* x) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (!expr_free_of(base, x) && exp->type == EXPR_INTEGER &&
            exp->data.integer < 0) {
            return true;
        }
    }
    if (has_negative_power_in(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_negative_power_in(e->data.function.args[i], x)) return true;
    }
    return false;
}

/* True iff every Power[base, negative_integer] subterm of `e` with `base`
 * depending on `x` has `base` as a polynomial in `x`. Non-polynomial
 * denominators (e.g. Exp[x] - 1 - x) would confuse Apart, whose semantics
 * are only well-defined for rational functions. Conservative: if any
 * suspicious subterm fails the polynomial test, return false. */
static bool all_negative_powers_polynomial(Expr* e, Expr* x) {
    if (!e) return true;
    if (e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (!expr_free_of(base, x) && exp->type == EXPR_INTEGER &&
            exp->data.integer < 0) {
            Expr* vars[1] = { x };
            if (!is_polynomial(base, vars, 1)) return false;
        }
    }
    if (!all_negative_powers_polynomial(e->data.function.head, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!all_negative_powers_polynomial(e->data.function.args[i], x)) return false;
    }
    return true;
}

/* Maxima-style rational-function preprocessing. If `f` contains a negative
 * integer power of an x-dependent subexpression (the marker of a rational
 * function in x), call Apart[f, x] to decompose it into partial fractions.
 * Returns a newly owned expression on success (caller must free), or NULL
 * when Apart either couldn't help or returned something equal to the input.
 *
 * The payoff: each partial-fraction term is easier for the series machinery
 * to expand. Terms like c/(a + b*x)^m route through the monomial binomial
 * fast path and produce closed-form coefficients instead of a full Newton
 * inversion over the composite polynomial denominator. */
/* True iff `e` contains a `Power[base, exp]` with base involving x and
 * exp that is neither an integer nor a literal rational. Apart does not
 * handle such subterms (they turn the input into a non-rational
 * function) and yields 0 for some shapes -- catch them here before
 * preprocessing. */
static bool has_non_rational_power_in(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (!expr_free_of(base, x)) {
            if (exp->type != EXPR_INTEGER) {
                int64_t p, q;
                if (!is_rational(exp, &p, &q)) return true;
            }
        }
    }
    if (has_non_rational_power_in(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_non_rational_power_in(e->data.function.args[i], x)) return true;
    }
    return false;
}

static Expr* try_apart_preprocess(Expr* f, Expr* x) {
    if (!has_negative_power_in(f, x)) return NULL;
    /* Apart is only defined for rational functions. If any negative-power
     * subterm has a non-polynomial base (e.g. 1/(Exp[x] - 1 - x), which has
     * a transcendental denominator), bail out rather than risk feeding
     * Apart something it can't cleanly handle. */
    if (!all_negative_powers_polynomial(f, x)) return NULL;
    /* Also refuse when any `Power[base, symbolic_exp]` with x in base has
     * a non-integer, non-rational exponent (e.g. x^a). Apart collapses
     * such inputs to 0 in Mathilda; better to skip preprocessing. */
    if (has_non_rational_power_in(f, x)) return NULL;
    Expr* args[2] = { expr_copy(f), expr_copy(x) };
    Expr* call = expr_new_function(mk_symbol("Apart"), args, 2);
    Expr* result = eval_and_free(call);
    if (!result) return NULL;
    if (expr_eq(result, f)) { expr_free(result); return NULL; }
    return result;
}

/* Scan coefs of s for the first non-zero entry at or beyond exponent
 * `from_exp`. Returns the exponent of that entry, or INT64_MAX if none. */
static int64_t so_first_nonzero_exp(SeriesObj* s, int64_t from_exp) {
    int64_t start_i = from_exp - s->nmin;
    if (start_i < 0) start_i = 0;
    for (size_t i = (size_t)start_i; i < s->coef_count; i++) {
        if (!is_lit_zero(s->coefs[i])) return s->nmin + (int64_t)i;
    }
    return INT64_MAX;
}

/* Asymptotic expansion of ExpIntegralEi[x] at x = Infinity (DLMF 6.12.2):
 *
 *   Ei(x) ~ E^x * Sum_{k>=0} k! / x^(k+1)
 *         = E^x ( 1/x + 1/x^2 + 2/x^3 + 6/x^4 + ... ).
 *
 * The leading E^x is an essential singularity and cannot live inside a
 * power series, so it stays as a symbolic prefactor (matching Mathematica,
 * which returns `E^x (1/x + 1/x^2 + 2/x^3 + O[1/x]^4)`) while the bracket is
 * a clean Laurent series in 1/x. We build the result expression directly:
 *
 *   Times[ Exp[x],
 *          SeriesData[1/x, 0, {0!, 1!, ..., (n-1)!}, 1, n+1, 1] ].
 *
 * Returns NULL (so the caller falls through) unless f is exactly
 * ExpIntegralEi[x] in the expansion variable. */
static Expr* try_series_ei_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "ExpIntegralEi") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* Coefficients k! for k = 0 .. n-1, placed at 1/x exponents 1 .. n. */
    Expr** coefs = calloc((size_t)n, sizeof(Expr*));
    for (int64_t k = 0; k < n; k++)
        coefs[k] = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, (size_t)n);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)       */
    sd[2] = coef_list;                                    /* coefficients      */
    sd[3] = expr_new_integer(1);                          /* nmin              */
    sd[4] = expr_new_integer(n + 1);                      /* nmax (O-term)     */
    sd[5] = expr_new_integer(1);                          /* denominator       */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);

    return mk_times(mk_fn1("Exp", expr_copy(x)), series);
}

/* Asymptotic expansion of SinIntegral[x] at x = Infinity (DLMF 6.12.3):
 *
 *   Si(x) ~ Pi/2 - cos(x) f(x) - sin(x) g(x),
 *     f(x) = Sum_{k>=0} (-1)^k (2k)!   / x^(2k+1),
 *     g(x) = Sum_{k>=0} (-1)^k (2k+1)! / x^(2k+2).
 *
 * cos(x)/sin(x) are essential-singularity prefactors that stay symbolic while
 * each bracket is a Laurent series in 1/x. Folding the leading minus signs in,
 * the result (matching Mathematica's Series[SinIntegral[x],{x,Infinity,k}]) is
 *
 *   Pi/2
 *   + Cos[x] SeriesData[1/x, 0, { (-1)^{k+1}(2k)!   at 1/x^(2k+1) }, 1, .., 1]
 *   + Sin[x] SeriesData[1/x, 0, { (-1)^{k+1}(2k+1)! at 1/x^(2k+2) }, 2, .., 1].
 *
 * Returns NULL unless f is exactly SinIntegral[x] in the expansion variable. */
static Expr* try_series_sinintegral_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "SinIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* f-part: powers 2k+1 (k>=0), populated at 1/x exponents 1, 3, 5, ...
     * coefficient (-1)^{k+1}(2k)!. nmin = 1. */
    size_t ncf = (size_t)n;                       /* 1/x exponents 1 .. n */
    Expr** cf = calloc(ncf, sizeof(Expr*));
    for (size_t i = 0; i < ncf; i++) cf[i] = expr_new_integer(0);
    for (int64_t k = 0; ; k++) {
        int64_t p = 2 * k + 1;                    /* 1/x exponent */
        if (p > n) break;
        Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k)));
        int64_t sign = (k % 2 == 0) ? -1 : 1;     /* (-1)^{k+1} */
        expr_free(cf[(size_t)(p - 1)]);
        cf[(size_t)(p - 1)] = simp(mk_times(expr_new_integer(sign), fact));
    }
    Expr** sdf = calloc(6, sizeof(Expr*));
    sdf[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sdf[1] = expr_new_integer(0);
    sdf[2] = expr_new_function(mk_symbol("List"), cf, ncf);
    sdf[3] = expr_new_integer(1);                 /* nmin = 1 */
    sdf[4] = expr_new_integer(n + 1);             /* O-term  */
    sdf[5] = expr_new_integer(1);
    Expr* seriesF = expr_new_function(mk_symbol("SeriesData"), sdf, 6);
    free(sdf); free(cf);

    Expr* halfpi = simp(mk_times(make_rational(1, 2), mk_symbol("Pi")));
    Expr* result = mk_plus(halfpi, mk_times(mk_fn1("Cos", expr_copy(x)), seriesF));

    /* g-part: powers 2k+2 (k>=0) at 1/x exponents 2, 4, 6, ...; nmin = 2.
     * Only emitted if any g-term fits within order n. */
    if (n >= 2) {
        size_t ncg = (size_t)(n - 1);             /* 1/x exponents 2 .. n, index p-2 */
        Expr** cg = calloc(ncg, sizeof(Expr*));
        for (size_t i = 0; i < ncg; i++) cg[i] = expr_new_integer(0);
        for (int64_t k = 0; ; k++) {
            int64_t p = 2 * k + 2;                /* 1/x exponent */
            if (p > n) break;
            Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k + 1)));
            int64_t sign = (k % 2 == 0) ? -1 : 1; /* (-1)^{k+1} */
            expr_free(cg[(size_t)(p - 2)]);
            cg[(size_t)(p - 2)] = simp(mk_times(expr_new_integer(sign), fact));
        }
        Expr** sdg = calloc(6, sizeof(Expr*));
        sdg[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdg[1] = expr_new_integer(0);
        sdg[2] = expr_new_function(mk_symbol("List"), cg, ncg);
        sdg[3] = expr_new_integer(2);             /* nmin = 2 */
        sdg[4] = expr_new_integer(n + 1);
        sdg[5] = expr_new_integer(1);
        Expr* seriesG = expr_new_function(mk_symbol("SeriesData"), sdg, 6);
        free(sdg); free(cg);
        result = mk_plus(result, mk_times(mk_fn1("Sin", expr_copy(x)), seriesG));
    }
    return result;
}

/* Asymptotic expansion of CosIntegral[x] at x = Infinity (DLMF 6.12.4):
 *
 *   Ci(x) ~ sin(x) f(x) - cos(x) g(x),
 *     f(x) = Sum_{k>=0} (-1)^k (2k)!   / x^(2k+1),
 *     g(x) = Sum_{k>=0} (-1)^k (2k+1)! / x^(2k+2).
 *
 * There is no constant term (unlike Si's Pi/2). Matching Mathematica's
 * Series[CosIntegral[x],{x,Infinity,k}]:
 *
 *   Sin[x] SeriesData[1/x, 0, { (-1)^k (2k)!     at 1/x^(2k+1) }, 1, .., 1]
 *   + Cos[x] SeriesData[1/x, 0, { (-1)^{k+1}(2k+1)! at 1/x^(2k+2) }, 2, .., 1].
 *
 * Returns NULL unless f is exactly CosIntegral[x] in the expansion variable. */
static Expr* try_series_cosintegral_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "CosIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* f-part (attached to Sin): powers 2k+1 at 1/x exponents 1, 3, 5, ...
     * coefficient (-1)^k (2k)!. nmin = 1. */
    size_t ncf = (size_t)n;                       /* 1/x exponents 1 .. n */
    Expr** cf = calloc(ncf, sizeof(Expr*));
    for (size_t i = 0; i < ncf; i++) cf[i] = expr_new_integer(0);
    for (int64_t k = 0; ; k++) {
        int64_t p = 2 * k + 1;                    /* 1/x exponent */
        if (p > n) break;
        Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k)));
        int64_t sign = (k % 2 == 0) ? 1 : -1;     /* (-1)^k */
        expr_free(cf[(size_t)(p - 1)]);
        cf[(size_t)(p - 1)] = simp(mk_times(expr_new_integer(sign), fact));
    }
    Expr** sdf = calloc(6, sizeof(Expr*));
    sdf[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sdf[1] = expr_new_integer(0);
    sdf[2] = expr_new_function(mk_symbol("List"), cf, ncf);
    sdf[3] = expr_new_integer(1);                 /* nmin = 1 */
    sdf[4] = expr_new_integer(n + 1);             /* O-term  */
    sdf[5] = expr_new_integer(1);
    Expr* seriesF = expr_new_function(mk_symbol("SeriesData"), sdf, 6);
    free(sdf); free(cf);

    Expr* result = mk_times(mk_fn1("Sin", expr_copy(x)), seriesF);

    /* g-part (attached to Cos): powers 2k+2 at 1/x exponents 2, 4, 6, ...;
     * coefficient (-1)^{k+1}(2k+1)!. nmin = 2. Only if it fits within order n. */
    if (n >= 2) {
        size_t ncg = (size_t)(n - 1);             /* 1/x exponents 2 .. n, index p-2 */
        Expr** cg = calloc(ncg, sizeof(Expr*));
        for (size_t i = 0; i < ncg; i++) cg[i] = expr_new_integer(0);
        for (int64_t k = 0; ; k++) {
            int64_t p = 2 * k + 2;                /* 1/x exponent */
            if (p > n) break;
            Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k + 1)));
            int64_t sign = (k % 2 == 0) ? -1 : 1; /* (-1)^{k+1} */
            expr_free(cg[(size_t)(p - 2)]);
            cg[(size_t)(p - 2)] = simp(mk_times(expr_new_integer(sign), fact));
        }
        Expr** sdg = calloc(6, sizeof(Expr*));
        sdg[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdg[1] = expr_new_integer(0);
        sdg[2] = expr_new_function(mk_symbol("List"), cg, ncg);
        sdg[3] = expr_new_integer(2);             /* nmin = 2 */
        sdg[4] = expr_new_integer(n + 1);
        sdg[5] = expr_new_integer(1);
        Expr* seriesG = expr_new_function(mk_symbol("SeriesData"), sdg, 6);
        free(sdg); free(cg);
        result = mk_plus(result, mk_times(mk_fn1("Cos", expr_copy(x)), seriesG));
    }
    return result;
}

/* Build the -I Pi/2 Stokes constant shared by the Shi/Chi asymptotic
 * expansions: Times[Complex[0,-1], Rational[1,2], Pi]. */
static Expr* si_hyper_stokes_const(void) {
    return simp(mk_times(make_complex(expr_new_integer(0), expr_new_integer(-1)),
                         mk_times(make_rational(1, 2), mk_symbol("Pi"))));
}

/* Asymptotic expansion of SinhIntegral[x] at x = Infinity:
 *
 *   Shi(x) ~ -I Pi/2 + cosh(x) F(x) + sinh(x) G(x),
 *     F(x) = Sum_{k>=0} (2k)!   / x^(2k+1),
 *     G(x) = Sum_{k>=0} (2k+1)! / x^(2k+2).
 *
 * F, G are Si's f, g without the (-1)^k (all coefficients positive), and the
 * combination pairs cosh with F, sinh with G. Matching Mathematica's
 * Series[SinhIntegral[x],{x,Infinity,k}]:
 *
 *   -I Pi/2
 *   + Cosh[x] SeriesData[1/x, 0, { (2k)!   at 1/x^(2k+1) }, 1, .., 1]
 *   + Sinh[x] SeriesData[1/x, 0, { (2k+1)! at 1/x^(2k+2) }, 2, .., 1].
 *
 * Returns NULL unless f is exactly SinhIntegral[x] in the expansion variable. */
static Expr* try_series_sinhintegral_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "SinhIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* F-part (attached to Cosh): powers 2k+1 at 1/x exponents 1, 3, 5, ...
     * coefficient (2k)!. nmin = 1. */
    size_t ncf = (size_t)n;                       /* 1/x exponents 1 .. n */
    Expr** cf = calloc(ncf, sizeof(Expr*));
    for (size_t i = 0; i < ncf; i++) cf[i] = expr_new_integer(0);
    for (int64_t k = 0; ; k++) {
        int64_t p = 2 * k + 1;                    /* 1/x exponent */
        if (p > n) break;
        Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k)));
        expr_free(cf[(size_t)(p - 1)]);
        cf[(size_t)(p - 1)] = fact;
    }
    Expr** sdf = calloc(6, sizeof(Expr*));
    sdf[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sdf[1] = expr_new_integer(0);
    sdf[2] = expr_new_function(mk_symbol("List"), cf, ncf);
    sdf[3] = expr_new_integer(1);                 /* nmin = 1 */
    sdf[4] = expr_new_integer(n + 1);             /* O-term  */
    sdf[5] = expr_new_integer(1);
    Expr* seriesF = expr_new_function(mk_symbol("SeriesData"), sdf, 6);
    free(sdf); free(cf);

    Expr* result = mk_plus(si_hyper_stokes_const(),
                           mk_times(mk_fn1("Cosh", expr_copy(x)), seriesF));

    /* G-part (attached to Sinh): powers 2k+2 at 1/x exponents 2, 4, 6, ...;
     * coefficient (2k+1)!. nmin = 2. Only emitted if any g-term fits order n. */
    if (n >= 2) {
        size_t ncg = (size_t)(n - 1);             /* 1/x exponents 2 .. n, index p-2 */
        Expr** cg = calloc(ncg, sizeof(Expr*));
        for (size_t i = 0; i < ncg; i++) cg[i] = expr_new_integer(0);
        for (int64_t k = 0; ; k++) {
            int64_t p = 2 * k + 2;                /* 1/x exponent */
            if (p > n) break;
            Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k + 1)));
            expr_free(cg[(size_t)(p - 2)]);
            cg[(size_t)(p - 2)] = fact;
        }
        Expr** sdg = calloc(6, sizeof(Expr*));
        sdg[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdg[1] = expr_new_integer(0);
        sdg[2] = expr_new_function(mk_symbol("List"), cg, ncg);
        sdg[3] = expr_new_integer(2);             /* nmin = 2 */
        sdg[4] = expr_new_integer(n + 1);
        sdg[5] = expr_new_integer(1);
        Expr* seriesG = expr_new_function(mk_symbol("SeriesData"), sdg, 6);
        free(sdg); free(cg);
        result = mk_plus(result, mk_times(mk_fn1("Sinh", expr_copy(x)), seriesG));
    }
    return result;
}

/* Asymptotic expansion of CoshIntegral[x] at x = Infinity:
 *
 *   Chi(x) ~ -I Pi/2 + sinh(x) F(x) + cosh(x) G(x),
 *     F(x) = Sum_{k>=0} (2k)!   / x^(2k+1),
 *     G(x) = Sum_{k>=0} (2k+1)! / x^(2k+2).
 *
 * Same F, G as SinhIntegral; only the combination differs (sinh with F, cosh
 * with G). Matching Mathematica's Series[CoshIntegral[x],{x,Infinity,k}]:
 *
 *   -I Pi/2
 *   + Sinh[x] SeriesData[1/x, 0, { (2k)!   at 1/x^(2k+1) }, 1, .., 1]
 *   + Cosh[x] SeriesData[1/x, 0, { (2k+1)! at 1/x^(2k+2) }, 2, .., 1].
 *
 * Returns NULL unless f is exactly CoshIntegral[x] in the expansion variable. */
static Expr* try_series_coshintegral_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "CoshIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* F-part (attached to Sinh): powers 2k+1 at 1/x exponents 1, 3, 5, ...
     * coefficient (2k)!. nmin = 1. */
    size_t ncf = (size_t)n;                       /* 1/x exponents 1 .. n */
    Expr** cf = calloc(ncf, sizeof(Expr*));
    for (size_t i = 0; i < ncf; i++) cf[i] = expr_new_integer(0);
    for (int64_t k = 0; ; k++) {
        int64_t p = 2 * k + 1;                    /* 1/x exponent */
        if (p > n) break;
        Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k)));
        expr_free(cf[(size_t)(p - 1)]);
        cf[(size_t)(p - 1)] = fact;
    }
    Expr** sdf = calloc(6, sizeof(Expr*));
    sdf[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sdf[1] = expr_new_integer(0);
    sdf[2] = expr_new_function(mk_symbol("List"), cf, ncf);
    sdf[3] = expr_new_integer(1);                 /* nmin = 1 */
    sdf[4] = expr_new_integer(n + 1);             /* O-term  */
    sdf[5] = expr_new_integer(1);
    Expr* seriesF = expr_new_function(mk_symbol("SeriesData"), sdf, 6);
    free(sdf); free(cf);

    Expr* result = mk_plus(si_hyper_stokes_const(),
                           mk_times(mk_fn1("Sinh", expr_copy(x)), seriesF));

    /* G-part (attached to Cosh): powers 2k+2 at 1/x exponents 2, 4, 6, ...;
     * coefficient (2k+1)!. nmin = 2. Only emitted if any g-term fits order n. */
    if (n >= 2) {
        size_t ncg = (size_t)(n - 1);             /* 1/x exponents 2 .. n, index p-2 */
        Expr** cg = calloc(ncg, sizeof(Expr*));
        for (size_t i = 0; i < ncg; i++) cg[i] = expr_new_integer(0);
        for (int64_t k = 0; ; k++) {
            int64_t p = 2 * k + 2;                /* 1/x exponent */
            if (p > n) break;
            Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(2 * k + 1)));
            expr_free(cg[(size_t)(p - 2)]);
            cg[(size_t)(p - 2)] = fact;
        }
        Expr** sdg = calloc(6, sizeof(Expr*));
        sdg[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdg[1] = expr_new_integer(0);
        sdg[2] = expr_new_function(mk_symbol("List"), cg, ncg);
        sdg[3] = expr_new_integer(2);             /* nmin = 2 */
        sdg[4] = expr_new_integer(n + 1);
        sdg[5] = expr_new_integer(1);
        Expr* seriesG = expr_new_function(mk_symbol("SeriesData"), sdg, 6);
        free(sdg); free(cg);
        result = mk_plus(result, mk_times(mk_fn1("Cosh", expr_copy(x)), seriesG));
    }
    return result;
}

/* One Laurent 1/x part of a Fresnel asymptotic expansion. Builds
 *   SeriesData[1/x, 0, {coefs}, nmin, n+1, 1]
 * whose coefficient at 1/x^p (p = 4j + nmin, j = 0, 1, ...) is
 *   base_sign * (-1)^j * D_j / Pi^(2j + (nmin+1)/2),
 * where D_j is (4j-1)!! for the f-series (use_R = 0) or (4j+1)!! for the
 * g-series (use_R = 1). Returns NULL if no term fits within order n. */
static Expr* fresnel_inf_part(Expr* x, int64_t n, int64_t nmin, int use_R,
                              int base_sign) {
    if (n < nmin) return NULL;
    size_t len = (size_t)(n - nmin + 1);
    Expr** c = calloc(len, sizeof(Expr*));
    for (size_t i = 0; i < len; i++) c[i] = expr_new_integer(0);
    int64_t pi_base = (nmin + 1) / 2;          /* 1 for nmin=1, 2 for nmin=3 */
    Expr* dfact = expr_new_integer(1);         /* D_0 = (-1)!! = 1!! = 1 */
    for (int64_t j = 0; ; j++) {
        int64_t p = 4 * j + nmin;
        if (p > n) break;
        int64_t sgn = base_sign * ((j % 2 == 0) ? 1 : -1);
        Expr* pipow = mk_power(mk_symbol("Pi"), expr_new_integer(-(2 * j + pi_base)));
        Expr* coef = simp(mk_times(mk_times(expr_new_integer(sgn), expr_copy(dfact)),
                                   pipow));
        expr_free(c[(size_t)(p - nmin)]);
        c[(size_t)(p - nmin)] = coef;
        /* Advance D_j -> D_{j+1}. */
        int64_t jn = j + 1;
        int64_t mult = use_R ? (4 * jn + 1) * (4 * jn - 1)
                             : (4 * jn - 1) * (4 * jn - 3);
        dfact = simp(mk_times(dfact, expr_new_integer(mult)));
    }
    expr_free(dfact);
    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sd[1] = expr_new_integer(0);
    sd[2] = expr_new_function(mk_symbol("List"), c, len);
    sd[3] = expr_new_integer(nmin);
    sd[4] = expr_new_integer(n + 1);
    sd[5] = expr_new_integer(1);
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd); free(c);
    return series;
}

/* phi = Pi x^2 / 2, the essential-singularity phase of the Fresnel integrals. */
static Expr* fresnel_inf_phase(Expr* x) {
    return simp(mk_times(make_rational(1, 2),
               mk_times(mk_symbol("Pi"), mk_power(expr_copy(x), expr_new_integer(2)))));
}

/* Asymptotic expansion of FresnelC[x] at x = Infinity (DLMF 7.12):
 *   C(x) ~ 1/2 + Sin[Pi x^2/2] f(x) - Cos[Pi x^2/2] g(x),
 *     f(x) = (1/(Pi x))    Sum_j (-1)^j (4j-1)!! / (Pi x^2)^(2j),
 *     g(x) = (1/(Pi^2 x^3)) Sum_j (-1)^j (4j+1)!! / (Pi x^2)^(2j).
 * Sin/Cos of the quadratic phase stay symbolic prefactors on Laurent series in
 * 1/x. Returns NULL unless f is exactly FresnelC[x] in the expansion variable. */
static Expr* try_series_fresnelc_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "FresnelC") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    Expr* result = make_rational(1, 2);
    Expr* fpart = fresnel_inf_part(x, n, 1, 0, +1);   /* +Sin * f */
    if (fpart) result = mk_plus(result, mk_times(mk_fn1("Sin", fresnel_inf_phase(x)), fpart));
    Expr* gpart = fresnel_inf_part(x, n, 3, 1, -1);   /* -Cos * g */
    if (gpart) result = mk_plus(result, mk_times(mk_fn1("Cos", fresnel_inf_phase(x)), gpart));
    return result;
}

/* Asymptotic expansion of FresnelS[x] at x = Infinity (DLMF 7.12):
 *   S(x) ~ 1/2 - Cos[Pi x^2/2] f(x) - Sin[Pi x^2/2] g(x),
 * with the same auxiliary f, g. Returns NULL unless f is exactly FresnelS[x]. */
static Expr* try_series_fresnels_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "FresnelS") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    Expr* result = make_rational(1, 2);
    Expr* fpart = fresnel_inf_part(x, n, 1, 0, -1);   /* -Cos * f */
    if (fpart) result = mk_plus(result, mk_times(mk_fn1("Cos", fresnel_inf_phase(x)), fpart));
    Expr* gpart = fresnel_inf_part(x, n, 3, 1, -1);   /* -Sin * g */
    if (gpart) result = mk_plus(result, mk_times(mk_fn1("Sin", fresnel_inf_phase(x)), gpart));
    return result;
}

/* Shared multiplier for the error-function asymptotic expansions at x = Infinity
 * (DLMF 7.12.1). Each of Erf/Erfc/Erfi is  prefactor(Exp[±x^2]) times a series
 *
 *   Sum_{k>=0} a_k / x^(2k+1),   a_0 = lead_sign / Sqrt[Pi],
 *   a_k = a_{k-1} * ratio_sign * (2k-1)/2.
 *
 * Only the odd 1/x powers (2k+1) are populated; the even slots are 0. The result
 * is a Laurent series in 1/x, matching Mathilda's Infinity convention (expansion
 * variable Power[x,-1], x0 = 0). Order n fills 1/x exponents 1..n with an O-term
 * at n+1, mirroring the ExpIntegralEi hook above. */
static Expr* erf_family_series(Expr* x, int64_t n, int lead_sign, int ratio_sign) {
    if (n < 1) n = 1;
    size_t ncoef = (size_t)n;                 /* 1/x exponents 1 .. n */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    /* a_0 = lead_sign * Pi^(-1/2), living at 1/x exponent 1 (index 0). */
    Expr* ak = eval_and_free(mk_times(expr_new_integer(lead_sign),
                                      mk_power(mk_symbol("Pi"), make_rational(-1, 2))));
    for (int64_t k = 0; ; k++) {
        size_t idx = (size_t)(2 * k + 1 - 1); /* exponent 2k+1, nmin = 1 */
        if (idx >= ncoef) { expr_free(ak); break; }
        expr_free(coefs[idx]);
        coefs[idx] = expr_copy(ak);
        /* a_{k+1} = a_k * ratio_sign * (2(k+1)-1)/2 = a_k * ratio_sign*(2k+1)/2. */
        ak = eval_and_free(mk_times(ak, make_rational(ratio_sign * (2 * k + 1), 2)));
    }

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)       */
    sd[2] = coef_list;                                    /* coefficients      */
    sd[3] = expr_new_integer(1);                          /* nmin              */
    sd[4] = expr_new_integer(n + 1);                      /* nmax (O-term)     */
    sd[5] = expr_new_integer(1);                          /* denominator       */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Build the Exp[±x^2] essential-singularity prefactor. */
static Expr* erf_exp_prefactor(Expr* x, int sign) {
    Expr* xsq = mk_power(expr_copy(x), expr_new_integer(2));
    if (sign < 0) xsq = mk_times(expr_new_integer(-1), xsq);
    return mk_fn1("Exp", xsq);
}

/* Asymptotic expansion of Erf[x] at x = Infinity (DLMF 7.12.1):
 *
 *   Erf(x) ~ 1 - E^(-x^2)/Sqrt[Pi] (1/x - 1/(2x^3) + 3/(4x^5) - ...).
 *
 * The constant 1 is the limit; the E^(-x^2) essential singularity stays a
 * symbolic prefactor multiplying a Laurent series in 1/x, so the result is
 *   Plus[1, Times[Exp[-x^2], SeriesData[1/x, 0, {-1/Sqrt[Pi], 0, ...}, 1, n+1, 1]]].
 * Returns NULL unless f is exactly Erf[x] in the expansion variable. */
static Expr* try_series_erf_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "Erf") || f->data.function.arg_count != 1) return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    Expr* mult = erf_family_series(x, n, -1, -1);   /* -1/Sqrt[Pi], alternating */
    return mk_plus(expr_new_integer(1),
                   mk_times(erf_exp_prefactor(x, -1), mult));
}

/* Erfc(x) = 1 - Erf(x) ~ E^(-x^2)/Sqrt[Pi] (1/x - 1/(2x^3) + ...). No constant
 * term (Erfc -> 0); the multiplier is the negation of Erf's. */
static Expr* try_series_erfc_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "Erfc") || f->data.function.arg_count != 1) return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    Expr* mult = erf_family_series(x, n, +1, -1);   /* +1/Sqrt[Pi], alternating */
    return mk_times(erf_exp_prefactor(x, -1), mult);
}

/* Erfi(x) = -I Erf(I x) ~ E^(x^2)/Sqrt[Pi] (1/x + 1/(2x^3) + 3/(4x^5) + ...).
 * All-positive coefficients, growing prefactor Exp[+x^2] (Erfi -> Infinity). */
static Expr* try_series_erfi_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "Erfi") || f->data.function.arg_count != 1) return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    Expr* mult = erf_family_series(x, n, +1, +1);   /* +1/Sqrt[Pi], all positive */
    return mk_times(erf_exp_prefactor(x, +1), mult);
}

/* Asymptotic expansion of AiryAi[x] at x = Infinity (DLMF 9.7.5):
 *
 *   Ai(x) ~ E^(-zeta) / (2 Sqrt[Pi] x^(1/4)) Sum_{k>=0} (-1)^k u_k / zeta^k,
 *           zeta = (2/3) x^(3/2),  u_k = (6k-5)(6k-3)(6k-1)/(216 k) u_{k-1}.
 *
 * Since zeta^(-k) = (3/2)^k x^(-3k/2), term k contributes
 *   (-1)^k u_k (3/2)^k / (2 Sqrt[Pi]) * x^(-(1+6k)/4),
 * a Puiseux series in (1/x) with denominator 4 (powers 1/4, 7/4, 13/4, ...).
 * The exponential E^(-zeta) is an essential singularity and stays a symbolic
 * prefactor (matching Mathematica), so the result is built directly as
 *
 *   Times[ Exp[-(2/3) x^(3/2)],
 *          SeriesData[1/x, 0, {coefs}, 1, 6(K+1)+1, 4] ],
 *
 * where K is the largest k whose power numerator 1+6k does not exceed 4n.
 * Returns NULL unless f is exactly AiryAi[x] in the expansion variable. */
static Expr* try_series_airyai_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "AiryAi") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* Largest k with x-power numerator (1 + 6k) <= 4n. */
    int64_t K = (4 * n - 1) / 6;
    if (K < 0) K = 0;

    size_t ncoef = (size_t)(6 * K + 1);     /* positions 0..6K, nonzero at 6k */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    Expr* uk = expr_new_integer(1);          /* u_0 = 1 (kept as an evaluated rational) */
    for (int64_t k = 0; k <= K; k++) {
        /* coef_k = (-1)^k u_k (3/2)^k / 2 * Pi^(-1/2). */
        Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? 1 : -1), expr_copy(uk));
        c = mk_times(c, mk_power(make_rational(3, 2), expr_new_integer(k)));
        c = mk_times(c, make_rational(1, 2));
        c = mk_times(c, mk_power(mk_symbol("Pi"), make_rational(-1, 2)));
        expr_free(coefs[(size_t)(6 * k)]);
        coefs[(size_t)(6 * k)] = eval_and_free(c);

        /* u_{k+1} = u_k (6m-5)(6m-3)(6m-1)/((2m-1) 216 m),  m = k+1 (DLMF 9.7.2). */
        if (k < K) {
            int64_t m = k + 1;
            Expr* num = expr_new_integer((6 * m - 5) * (6 * m - 3) * (6 * m - 1));
            Expr* den = mk_power(expr_new_integer((2 * m - 1) * 216 * m), expr_new_integer(-1));
            Expr* nxt = mk_times(uk, mk_times(num, den));
            uk = eval_and_free(nxt);
        } else {
            expr_free(uk);
        }
    }

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)       */
    sd[2] = coef_list;                                    /* coefficients      */
    sd[3] = expr_new_integer(1);                          /* nmin (x^(-1/4))   */
    sd[4] = expr_new_integer(6 * (K + 1) + 1);            /* nmax (O-term)     */
    sd[5] = expr_new_integer(4);                          /* denominator       */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);

    /* Prefactor Exp[-(2/3) x^(3/2)]. */
    Expr* zeta = mk_times(make_rational(-2, 3), mk_power(expr_copy(x), make_rational(3, 2)));
    return mk_times(mk_fn1("Exp", zeta), series);
}

/* Asymptotic expansion of AiryAiPrime[x] at x = Infinity (DLMF 9.7.6):
 *
 *   Ai'(x) ~ -x^(1/4) E^(-zeta) / (2 Sqrt[Pi]) Sum_{k>=0} (-1)^k v_k / zeta^k,
 *            zeta = (2/3) x^(3/2),  v_k = -(6k+1)/(6k-1) u_k,
 *   with u_k the DLMF 9.7.2 coefficients (v_0 = 1).
 *
 * Since zeta^(-k) = (3/2)^k x^(-3k/2), term k contributes
 *   -(-1)^k v_k (3/2)^k / (2 Sqrt[Pi]) * x^((1-6k)/4),
 * a Puiseux series in (1/x) with denominator 4 and powers 1/4, -5/4, -11/4, ...
 * i.e. numerators 1-6k = -1+6m (starting at +1, then descending). Relative to
 * the AiryAi expansion the leading power is x^(+1/4) (nmin = -1, vs +1 for Ai)
 * and the coefficients use v_k. The result is built as
 *
 *   Times[ Exp[-(2/3) x^(3/2)],
 *          SeriesData[1/x, 0, {coefs}, -1, 6(K+1)-1, 4] ].
 *
 * Coefficient at position 6k is -(1/(2 Sqrt[Pi])) (-1)^k v_k (3/2)^k, matching
 * the reference values -1/(2 Sqrt[Pi]), -7/(96 Sqrt[Pi]), 455/(9216 Sqrt[Pi]).
 * Returns NULL unless f is exactly AiryAiPrime[x] in the expansion variable. */
static Expr* try_series_airyaiprime_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "AiryAiPrime") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* Largest k whose x-power numerator (6k-1) <= 4n. */
    int64_t K = (4 * n + 1) / 6;
    if (K < 0) K = 0;

    size_t ncoef = (size_t)(6 * K + 1);     /* positions 0..6K, nonzero at 6k */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    Expr* uk = expr_new_integer(1);          /* u_0 = 1 (kept as an evaluated rational) */
    for (int64_t k = 0; k <= K; k++) {
        /* v_k = -(6k+1)/(6k-1) u_k; coef_k = -(-1)^k v_k (3/2)^k / 2 * Pi^(-1/2). */
        Expr* vk = mk_times(make_rational(-(6 * k + 1), 6 * k - 1), expr_copy(uk));
        Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? -1 : 1), vk);  /* -(-1)^k v_k */
        c = mk_times(c, mk_power(make_rational(3, 2), expr_new_integer(k)));
        c = mk_times(c, make_rational(1, 2));
        c = mk_times(c, mk_power(mk_symbol("Pi"), make_rational(-1, 2)));
        expr_free(coefs[(size_t)(6 * k)]);
        coefs[(size_t)(6 * k)] = eval_and_free(c);

        /* u_{k+1} = u_k (6m-5)(6m-3)(6m-1)/((2m-1) 216 m),  m = k+1 (DLMF 9.7.2). */
        if (k < K) {
            int64_t m = k + 1;
            Expr* num = expr_new_integer((6 * m - 5) * (6 * m - 3) * (6 * m - 1));
            Expr* den = mk_power(expr_new_integer((2 * m - 1) * 216 * m), expr_new_integer(-1));
            Expr* nxt = mk_times(uk, mk_times(num, den));
            uk = eval_and_free(nxt);
        } else {
            expr_free(uk);
        }
    }

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x   */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)         */
    sd[2] = coef_list;                                    /* coefficients        */
    sd[3] = expr_new_integer(-1);                         /* nmin (x^(+1/4))     */
    sd[4] = expr_new_integer(6 * (K + 1) - 1);            /* nmax (O-term)       */
    sd[5] = expr_new_integer(4);                          /* denominator         */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);

    /* Prefactor Exp[-(2/3) x^(3/2)]. */
    Expr* zeta = mk_times(make_rational(-2, 3), mk_power(expr_copy(x), make_rational(3, 2)));
    return mk_times(mk_fn1("Exp", zeta), series);
}

/* Asymptotic expansion of AiryBi[x] at x = Infinity (DLMF 9.7.7):
 *
 *   Bi(x) ~ E^(zeta) / (Sqrt[Pi] x^(1/4)) Sum_{k>=0} u_k / zeta^k,
 *           zeta = (2/3) x^(3/2),  u_k as in DLMF 9.7.2.
 *
 * Same Puiseux structure as the Ai expansion above, but the dominant Bi series
 * carries no (-1)^k sign, a prefactor 1/Sqrt[Pi] (not 1/(2 Sqrt[Pi])), and the
 * exponential E^(+zeta). Term k contributes
 *   u_k (3/2)^k / Sqrt[Pi] * x^(-(1+6k)/4),
 * a Puiseux series in (1/x) with denominator 4. The result is built as
 *
 *   Times[ Exp[(2/3) x^(3/2)],
 *          SeriesData[1/x, 0, {coefs}, 1, 6(K+1)+1, 4] ].
 *
 * Returns NULL unless f is exactly AiryBi[x] in the expansion variable. */
static Expr* try_series_airybi_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "AiryBi") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* Largest k with x-power numerator (1 + 6k) <= 4n. */
    int64_t K = (4 * n - 1) / 6;
    if (K < 0) K = 0;

    size_t ncoef = (size_t)(6 * K + 1);     /* positions 0..6K, nonzero at 6k */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    Expr* uk = expr_new_integer(1);          /* u_0 = 1 (kept as an evaluated rational) */
    for (int64_t k = 0; k <= K; k++) {
        /* coef_k = u_k (3/2)^k * Pi^(-1/2)  (no (-1)^k, full 1/Sqrt[Pi]). */
        Expr* c = expr_copy(uk);
        c = mk_times(c, mk_power(make_rational(3, 2), expr_new_integer(k)));
        c = mk_times(c, mk_power(mk_symbol("Pi"), make_rational(-1, 2)));
        expr_free(coefs[(size_t)(6 * k)]);
        coefs[(size_t)(6 * k)] = eval_and_free(c);

        /* u_{k+1} = u_k (6m-5)(6m-3)(6m-1)/((2m-1) 216 m),  m = k+1 (DLMF 9.7.2). */
        if (k < K) {
            int64_t m = k + 1;
            Expr* num = expr_new_integer((6 * m - 5) * (6 * m - 3) * (6 * m - 1));
            Expr* den = mk_power(expr_new_integer((2 * m - 1) * 216 * m), expr_new_integer(-1));
            Expr* nxt = mk_times(uk, mk_times(num, den));
            uk = eval_and_free(nxt);
        } else {
            expr_free(uk);
        }
    }

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)       */
    sd[2] = coef_list;                                    /* coefficients      */
    sd[3] = expr_new_integer(1);                          /* nmin (x^(-1/4))   */
    sd[4] = expr_new_integer(6 * (K + 1) + 1);            /* nmax (O-term)     */
    sd[5] = expr_new_integer(4);                          /* denominator       */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);

    /* Prefactor Exp[(2/3) x^(3/2)]. */
    Expr* zeta = mk_times(make_rational(2, 3), mk_power(expr_copy(x), make_rational(3, 2)));
    return mk_times(mk_fn1("Exp", zeta), series);
}

/* Asymptotic expansion of AiryBiPrime[x] at x = Infinity (DLMF 9.7.8):
 *
 *   Bi'(x) ~ x^(1/4) E^(zeta) / Sqrt[Pi] Sum_{k>=0} v_k / zeta^k,
 *            zeta = (2/3) x^(3/2),  v_k = -(6k+1)/(6k-1) u_k,  v_0 = 1,
 *   with u_k the DLMF 9.7.2 coefficients.
 *
 * The dominant Bi' series is the AiryBi expansion above but with coefficient
 * v_k (not u_k) and a leading power x^(+1/4) (vs x^(-1/4) for Bi): no (-1)^k,
 * prefactor 1/Sqrt[Pi], exponential E^(+zeta). Since zeta^(-k) = (3/2)^k
 * x^(-3k/2), term k contributes
 *   v_k (3/2)^k / Sqrt[Pi] * x^((1-6k)/4),
 * a Puiseux series in (1/x) with denominator 4 and powers 1/4, -5/4, -11/4, ...
 * i.e. numerators 1-6k = -1+6m (starting at +1, then descending, nmin = -1).
 * The result is built as
 *
 *   Times[ Exp[(2/3) x^(3/2)],
 *          SeriesData[1/x, 0, {coefs}, -1, 6(K+1)-1, 4] ].
 *
 * Coefficient at position 6k is (1/Sqrt[Pi]) v_k (3/2)^k, matching the
 * reference values 1/Sqrt[Pi], -7/(48 Sqrt[Pi]), -455/(4608 Sqrt[Pi]).
 * Returns NULL unless f is exactly AiryBiPrime[x] in the expansion variable. */
static Expr* try_series_airybiprime_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "AiryBiPrime") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* Largest k whose x-power numerator (6k-1) <= 4n. */
    int64_t K = (4 * n + 1) / 6;
    if (K < 0) K = 0;

    size_t ncoef = (size_t)(6 * K + 1);     /* positions 0..6K, nonzero at 6k */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    Expr* uk = expr_new_integer(1);          /* u_0 = 1 (kept as an evaluated rational) */
    for (int64_t k = 0; k <= K; k++) {
        /* v_k = -(6k+1)/(6k-1) u_k; coef_k = v_k (3/2)^k / Sqrt[Pi]. */
        Expr* vk = mk_times(make_rational(-(6 * k + 1), 6 * k - 1), expr_copy(uk));
        Expr* c = mk_times(vk, mk_power(make_rational(3, 2), expr_new_integer(k)));
        c = mk_times(c, mk_power(mk_symbol("Pi"), make_rational(-1, 2)));
        expr_free(coefs[(size_t)(6 * k)]);
        coefs[(size_t)(6 * k)] = eval_and_free(c);

        /* u_{k+1} = u_k (6m-5)(6m-3)(6m-1)/((2m-1) 216 m),  m = k+1 (DLMF 9.7.2). */
        if (k < K) {
            int64_t m = k + 1;
            Expr* num = expr_new_integer((6 * m - 5) * (6 * m - 3) * (6 * m - 1));
            Expr* den = mk_power(expr_new_integer((2 * m - 1) * 216 * m), expr_new_integer(-1));
            Expr* nxt = mk_times(uk, mk_times(num, den));
            uk = eval_and_free(nxt);
        } else {
            expr_free(uk);
        }
    }

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1)); /* expansion var 1/x   */
    sd[1] = expr_new_integer(0);                          /* x0 (in 1/x)         */
    sd[2] = coef_list;                                    /* coefficients        */
    sd[3] = expr_new_integer(-1);                         /* nmin (x^(+1/4))     */
    sd[4] = expr_new_integer(6 * (K + 1) - 1);            /* nmax (O-term)       */
    sd[5] = expr_new_integer(4);                          /* denominator         */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);

    /* Prefactor Exp[(2/3) x^(3/2)]. */
    Expr* zeta = mk_times(make_rational(2, 3), mk_power(expr_copy(x), make_rational(3, 2)));
    return mk_times(mk_fn1("Exp", zeta), series);
}

/* Asymptotic expansion of BesselJ[nu, x] at x = Infinity (DLMF 10.17.3):
 *
 *   J_nu(x) ~ sqrt(2/(pi x)) [ cos(w0) P(1/x) + sin(w0) Q(1/x) ],
 *     w0 = (2 nu + 1) pi/4 - x,
 *     P  = Sum_{k>=0} (-1)^k a_{2k} / x^{2k},
 *     Q  = Sum_{k>=0} (-1)^k a_{2k+1} / x^{2k+1},
 *     a_0 = 1,  a_j = a_{j-1} (mu - (2j-1)^2) / (8 j),  mu = 4 nu^2.
 *
 * Folding sqrt(2/(pi x)) in, the two oscillatory parts become Puiseux series
 * in 1/x with half-integer powers (denominator 2):
 *   P-part powers  1/2, 5/2, 9/2, ...  (nmin = 1, nonzero at index 4k)
 *   Q-part powers  3/2, 7/2, 11/2,...  (nmin = 3, nonzero at index 4k)
 * with coefficient sqrt(2/pi) (-1)^k a_{2k} (resp. a_{2k+1}). The result is
 *
 *   Cos[w0] SeriesData[1/x, 0, {Pcoefs}, 1, ..., 2]
 *   + Sin[w0] SeriesData[1/x, 0, {Qcoefs}, 3, ..., 2].
 *
 * The coefficients are polynomials in nu (via mu), so this works for symbolic
 * order too, matching Mathematica's Series[BesselJ[n,x],{x,Infinity,k}].
 * Returns NULL unless f is exactly BesselJ[nu, x] in the expansion variable. */
static Expr* try_series_besselj_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "BesselJ") || f->data.function.arg_count != 2)
        return NULL;
    Expr* nu = f->data.function.args[0];
    if (!expr_eq(f->data.function.args[1], x)) return NULL;

    /* Largest k with the relevant power <= n: P power (1+4k)/2, Q power (3+4k)/2. */
    int64_t kP_max = -1, kQ_max = -1;
    while (1 + 4 * (kP_max + 1) <= 2 * n) kP_max++;
    while (3 + 4 * (kQ_max + 1) <= 2 * n) kQ_max++;

    int64_t maxidx = 2 * kP_max;               /* highest a-index needed */
    if (kQ_max >= 0 && 2 * kQ_max + 1 > maxidx) maxidx = 2 * kQ_max + 1;

    /* mu = 4 nu^2 (symbolic if nu is). */
    Expr* mu = simp(mk_times(expr_new_integer(4),
                             mk_power(expr_copy(nu), expr_new_integer(2))));
    Expr* sqrt2pi = simp(mk_power(mk_times(expr_new_integer(2),
                             mk_power(mk_symbol("Pi"), expr_new_integer(-1))),
                             make_rational(1, 2)));   /* Sqrt[2/Pi] */

    /* a-sequence a_0..a_maxidx, evaluated. */
    Expr** a = calloc((size_t)maxidx + 1, sizeof(Expr*));
    a[0] = expr_new_integer(1);
    for (int64_t j = 1; j <= maxidx; j++) {
        int64_t odd = 2 * j - 1;
        Expr* num = mk_plus(expr_copy(mu), expr_new_integer(-(odd * odd)));
        Expr* den = mk_power(expr_new_integer(8 * j), expr_new_integer(-1));
        a[j] = simp(mk_times(expr_copy(a[j - 1]), mk_times(num, den)));
    }

    /* P coefficients (positions 4k, k = 0..kP_max). */
    size_t ncoefP = (size_t)(4 * kP_max + 1);
    Expr** coefP = calloc(ncoefP, sizeof(Expr*));
    for (size_t i = 0; i < ncoefP; i++) coefP[i] = expr_new_integer(0);
    for (int64_t k = 0; k <= kP_max; k++) {
        Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? 1 : -1),
                           expr_copy(a[2 * k]));
        c = mk_times(c, expr_copy(sqrt2pi));
        expr_free(coefP[(size_t)(4 * k)]);
        coefP[(size_t)(4 * k)] = simp(c);
    }
    Expr** sdP = calloc(6, sizeof(Expr*));
    sdP[0] = mk_power(expr_copy(x), expr_new_integer(-1));   /* 1/x */
    sdP[1] = expr_new_integer(0);
    sdP[2] = expr_new_function(mk_symbol("List"), coefP, ncoefP);
    sdP[3] = expr_new_integer(1);                            /* nmin = 1 */
    sdP[4] = expr_new_integer(1 + 4 * (kP_max + 1));         /* nmax */
    sdP[5] = expr_new_integer(2);                            /* den  */
    Expr* seriesP = expr_new_function(mk_symbol("SeriesData"), sdP, 6);
    free(sdP); free(coefP);

    /* w0 = (2 nu + 1) pi/4 - x. */
    Expr* twonu1 = mk_plus(mk_times(expr_new_integer(2), expr_copy(nu)),
                           expr_new_integer(1));
    Expr* w0 = simp(mk_plus(mk_times(make_rational(1, 4),
                                     mk_times(twonu1, mk_symbol("Pi"))),
                            mk_times(expr_new_integer(-1), expr_copy(x))));

    Expr* result = mk_times(mk_fn1("Cos", expr_copy(w0)), seriesP);

    /* Q part (Sin), only if it carries any term. */
    if (kQ_max >= 0) {
        size_t ncoefQ = (size_t)(4 * kQ_max + 1);
        Expr** coefQ = calloc(ncoefQ, sizeof(Expr*));
        for (size_t i = 0; i < ncoefQ; i++) coefQ[i] = expr_new_integer(0);
        for (int64_t k = 0; k <= kQ_max; k++) {
            Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? 1 : -1),
                               expr_copy(a[2 * k + 1]));
            c = mk_times(c, expr_copy(sqrt2pi));
            expr_free(coefQ[(size_t)(4 * k)]);
            coefQ[(size_t)(4 * k)] = simp(c);
        }
        Expr** sdQ = calloc(6, sizeof(Expr*));
        sdQ[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdQ[1] = expr_new_integer(0);
        sdQ[2] = expr_new_function(mk_symbol("List"), coefQ, ncoefQ);
        sdQ[3] = expr_new_integer(3);                       /* nmin = 3 */
        sdQ[4] = expr_new_integer(3 + 4 * (kQ_max + 1));    /* nmax */
        sdQ[5] = expr_new_integer(2);
        Expr* seriesQ = expr_new_function(mk_symbol("SeriesData"), sdQ, 6);
        free(sdQ); free(coefQ);
        result = mk_plus(result, mk_times(mk_fn1("Sin", expr_copy(w0)), seriesQ));
    }

    expr_free(w0);
    expr_free(mu);
    expr_free(sqrt2pi);
    for (int64_t j = 0; j <= maxidx; j++) expr_free(a[j]);
    free(a);
    return result;
}

/* Asymptotic expansion of BesselY[nu, x] at x = Infinity (DLMF 10.17.4):
 *
 *   Y_nu(x) ~ sqrt(2/(pi x)) [ sin(w0) P(1/x) + cos(w0) Q(1/x) ],
 *     w0 = (2 nu + 1) pi/4 - x,
 * with the SAME P, Q, a_j as BesselJ (try_series_besselj_at_infinity); only
 * the trig prefactors swap (J uses cos P - ...; Y uses sin P + cos Q). Returns
 * NULL unless f is exactly BesselY[nu, x] in the expansion variable. */
static Expr* try_series_bessely_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "BesselY") || f->data.function.arg_count != 2)
        return NULL;
    Expr* nu = f->data.function.args[0];
    if (!expr_eq(f->data.function.args[1], x)) return NULL;

    /* Largest k with the relevant power <= n: P power (1+4k)/2, Q power (3+4k)/2. */
    int64_t kP_max = -1, kQ_max = -1;
    while (1 + 4 * (kP_max + 1) <= 2 * n) kP_max++;
    while (3 + 4 * (kQ_max + 1) <= 2 * n) kQ_max++;

    int64_t maxidx = 2 * kP_max;               /* highest a-index needed */
    if (kQ_max >= 0 && 2 * kQ_max + 1 > maxidx) maxidx = 2 * kQ_max + 1;

    /* mu = 4 nu^2 (symbolic if nu is). */
    Expr* mu = simp(mk_times(expr_new_integer(4),
                             mk_power(expr_copy(nu), expr_new_integer(2))));
    Expr* sqrt2pi = simp(mk_power(mk_times(expr_new_integer(2),
                             mk_power(mk_symbol("Pi"), expr_new_integer(-1))),
                             make_rational(1, 2)));   /* Sqrt[2/Pi] */

    /* a-sequence a_0..a_maxidx, evaluated. */
    Expr** a = calloc((size_t)maxidx + 1, sizeof(Expr*));
    a[0] = expr_new_integer(1);
    for (int64_t j = 1; j <= maxidx; j++) {
        int64_t odd = 2 * j - 1;
        Expr* num = mk_plus(expr_copy(mu), expr_new_integer(-(odd * odd)));
        Expr* den = mk_power(expr_new_integer(8 * j), expr_new_integer(-1));
        a[j] = simp(mk_times(expr_copy(a[j - 1]), mk_times(num, den)));
    }

    /* P coefficients (positions 4k, k = 0..kP_max). */
    size_t ncoefP = (size_t)(4 * kP_max + 1);
    Expr** coefP = calloc(ncoefP, sizeof(Expr*));
    for (size_t i = 0; i < ncoefP; i++) coefP[i] = expr_new_integer(0);
    for (int64_t k = 0; k <= kP_max; k++) {
        Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? 1 : -1),
                           expr_copy(a[2 * k]));
        c = mk_times(c, expr_copy(sqrt2pi));
        expr_free(coefP[(size_t)(4 * k)]);
        coefP[(size_t)(4 * k)] = simp(c);
    }
    Expr** sdP = calloc(6, sizeof(Expr*));
    sdP[0] = mk_power(expr_copy(x), expr_new_integer(-1));   /* 1/x */
    sdP[1] = expr_new_integer(0);
    sdP[2] = expr_new_function(mk_symbol("List"), coefP, ncoefP);
    sdP[3] = expr_new_integer(1);                            /* nmin = 1 */
    sdP[4] = expr_new_integer(1 + 4 * (kP_max + 1));         /* nmax */
    sdP[5] = expr_new_integer(2);                            /* den  */
    Expr* seriesP = expr_new_function(mk_symbol("SeriesData"), sdP, 6);
    free(sdP); free(coefP);

    /* w0 = (2 nu + 1) pi/4 - x. */
    Expr* twonu1 = mk_plus(mk_times(expr_new_integer(2), expr_copy(nu)),
                           expr_new_integer(1));
    Expr* w0 = simp(mk_plus(mk_times(make_rational(1, 4),
                                     mk_times(twonu1, mk_symbol("Pi"))),
                            mk_times(expr_new_integer(-1), expr_copy(x))));

    /* DLMF 10.17.4: Y_nu ~ sqrt(2/(pi x))[sin(w) P + cos(w) Q], w = x-(2nu+1)pi/4.
     * With w0 = (2nu+1)pi/4 - x = -w: sin(w) = -sin(w0), cos(w) = cos(w0), so
     * the P term carries a leading minus: Y = -sin(w0) P + cos(w0) Q. (For J the
     * P term used cos(w0) = cos(w), so no sign flip was needed there.) */
    Expr* result = mk_times(expr_new_integer(-1),
                            mk_times(mk_fn1("Sin", expr_copy(w0)), seriesP));

    /* Q part (Cos), only if it carries any term. */
    if (kQ_max >= 0) {
        size_t ncoefQ = (size_t)(4 * kQ_max + 1);
        Expr** coefQ = calloc(ncoefQ, sizeof(Expr*));
        for (size_t i = 0; i < ncoefQ; i++) coefQ[i] = expr_new_integer(0);
        for (int64_t k = 0; k <= kQ_max; k++) {
            Expr* c = mk_times(expr_new_integer((k % 2 == 0) ? 1 : -1),
                               expr_copy(a[2 * k + 1]));
            c = mk_times(c, expr_copy(sqrt2pi));
            expr_free(coefQ[(size_t)(4 * k)]);
            coefQ[(size_t)(4 * k)] = simp(c);
        }
        Expr** sdQ = calloc(6, sizeof(Expr*));
        sdQ[0] = mk_power(expr_copy(x), expr_new_integer(-1));
        sdQ[1] = expr_new_integer(0);
        sdQ[2] = expr_new_function(mk_symbol("List"), coefQ, ncoefQ);
        sdQ[3] = expr_new_integer(3);                       /* nmin = 3 */
        sdQ[4] = expr_new_integer(3 + 4 * (kQ_max + 1));    /* nmax */
        sdQ[5] = expr_new_integer(2);
        Expr* seriesQ = expr_new_function(mk_symbol("SeriesData"), sdQ, 6);
        free(sdQ); free(coefQ);
        result = mk_plus(result, mk_times(mk_fn1("Cos", expr_copy(w0)), seriesQ));
    }

    expr_free(w0);
    expr_free(mu);
    expr_free(sqrt2pi);
    for (int64_t j = 0; j <= maxidx; j++) expr_free(a[j]);
    free(a);
    return result;
}

/* Generalized asymptotic series of LogIntegral[x] at x = 0.
 *
 * li(x) = Ei(Log[x]). With L = Log[x] -> -Infinity as x -> 0+, the function
 * lands in Ei's asymptotic regime (DLMF 6.12.2), giving
 *
 *   li(x) ~ E^L * Sum_{k>=0} k!/L^(k+1)
 *         = x * ( 1/L + 1/L^2 + 2/L^3 + 6/L^4 + ... )
 *
 * since E^L = x. The ONLY x-dependence is the prefactor E^L = x, so every
 * term carries exactly x^1: there is no ordinary power series here -- the x^1
 * coefficient is itself a (divergent, asymptotic) series in 1/Log[x]. We keep
 * 2n+2 terms of that 1/L series (k = 0 .. 2n+1), which reproduces
 * Mathematica's Series[LogIntegral[x], {x, 0, n}] output in its
 * Assumptions -> x > 0 form, e.g. for n = 2:
 *
 *   ((120 + 24 Log[x] + 6 Log[x]^2 + 2 Log[x]^3 + Log[x]^4 + Log[x]^5) x)
 *     / Log[x]^6 + O[x]^3.
 *
 * Result shape (only the x^1 slot is non-zero):
 *
 *   SeriesData[x, 0, {Together[Sum k!/Log[x]^(k+1)], 0, ..., 0}, 1, n+1, 1].
 *
 * Mathematica's no-assumptions form wraps this in a Floor[Arg[...]] branch
 * discriminator; we always emit the principal (x > 0) form. Returns NULL
 * unless f is exactly LogIntegral[x] in the expansion variable. */
static Expr* try_series_logintegral_at_zero(Expr* f, Expr* x, int64_t n, int x_sign) {
    if (!has_symbol_head(f, "LogIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    if (n < 1) n = 1;

    /* On the principal (x > 0) branch the log argument is x; when the
     * Assumptions option forces x < 0, li(x) = Ei(Log[x]) picks up the
     * I*Pi from Log[x] = Log[-x] + I*Pi: every Log[x] becomes Log[-x] and
     * an additive I*Pi (the x^0 term) appears in front of the series. */
    bool neg = (x_sign < 0);

    /* x^1 coefficient: Sum_{k=0}^{2n+1} k! / Log[.]^(k+1). */
    int64_t kmax = 2 * n + 1;
    Expr* logx = neg ? mk_fn1("Log", mk_times(expr_new_integer(-1), expr_copy(x)))
                     : mk_fn1("Log", expr_copy(x));
    Expr** terms = calloc((size_t)kmax + 1, sizeof(Expr*));
    for (int64_t k = 0; k <= kmax; k++) {
        Expr* fact = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* logpow = mk_power(expr_copy(logx), expr_new_integer(-(k + 1)));
        terms[k] = mk_times(fact, logpow);
    }
    expr_free(logx);
    Expr* sum = expr_new_function(mk_symbol("Plus"), terms, (size_t)kmax + 1);
    free(terms);
    /* Combine over a common Log[x]^(2n+2) denominator to match MMA's form;
     * Together falls through to the raw sum if it cannot combine. */
    Expr* coef1 = eval_and_free(mk_fn1("Together", sum));

    /* Principal branch: coefficients span exponents 1 .. n (length n, leading
     * x^1). For x < 0 the additive I*Pi is an x^0 term, so the series starts
     * at exponent 0: coefs = {I*Pi, coef1, 0, ..., 0}, length n + 1. */
    int64_t nmin = neg ? 0 : 1;
    size_t ncoef = neg ? (size_t)n + 1 : (size_t)n;
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    size_t slot = 0;
    if (neg) coefs[slot++] = simp(mk_times(mk_symbol("I"), mk_symbol("Pi")));
    coefs[slot++] = coef1;
    for (size_t i = slot; i < ncoef; i++) coefs[i] = expr_new_integer(0);
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);            /* expansion var      */
    sd[1] = expr_new_integer(0);     /* x0                 */
    sd[2] = coef_list;               /* coefficients       */
    sd[3] = expr_new_integer(nmin);  /* nmin               */
    sd[4] = expr_new_integer(n + 1); /* nmax (O-term)      */
    sd[5] = expr_new_integer(1);     /* denominator        */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Logarithmic series of ExpIntegralEi[x] at x = 0 (DLMF 6.6.2):
 *
 *   Ei(x) = EulerGamma + Log[x] + Sum_{k>=1} x^k / (k * k!)
 *         = EulerGamma + Log[x] + x + x^2/4 + x^3/18 + ...
 *
 * This is a convergent series with a single logarithmic term. The branch
 * term Log[x] (and the additive EulerGamma) is baked into the x^0
 * coefficient; the remaining coefficients are the regular Taylor
 * coefficients 1/(k*k!). Naive Taylor-via-D cannot produce it because
 * f(0) = Ei(0) = -Infinity. Result:
 *
 *   SeriesData[x, 0, {EulerGamma + Log[x], 1, 1/4, ...}, 0, n+1, 1].
 *
 * Returns NULL unless f is exactly ExpIntegralEi[x] in the expansion var. */
static Expr* try_series_ei_at_zero(Expr* f, Expr* x, int64_t n, int x_sign) {
    if (!has_symbol_head(f, "ExpIntegralEi") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    if (n < 0) n = 0;

    size_t ncoef = (size_t)n + 1;            /* exponents 0 .. n */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    /* x^0 coefficient: EulerGamma + Log[x] on the principal (x > 0) branch,
     * or EulerGamma + Log[-x] when the Assumptions option forces x < 0. */
    Expr* logarg = (x_sign < 0)
                 ? mk_times(expr_new_integer(-1), expr_copy(x))
                 : expr_copy(x);
    coefs[0] = simp(mk_plus(mk_symbol("EulerGamma"), mk_fn1("Log", logarg)));
    /* x^k coefficient: 1/(k*k!), exact rational via the bigint Factorial. */
    for (int64_t k = 1; k <= n; k++) {
        Expr* kfact = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* denom = eval_and_free(mk_times(expr_new_integer(k), kfact));
        coefs[k] = eval_and_free(mk_power(denom, expr_new_integer(-1)));
    }
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);            /* expansion var */
    sd[1] = expr_new_integer(0);     /* x0            */
    sd[2] = coef_list;               /* coefficients  */
    sd[3] = expr_new_integer(0);     /* nmin          */
    sd[4] = expr_new_integer(n + 1); /* nmax (O-term) */
    sd[5] = expr_new_integer(1);     /* denominator   */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* CosIntegral[x] at x = 0: Ci(0) = -Infinity blocks naive Taylor. Emit the
 * DLMF 6.6.2 logarithmic series directly:
 *
 *   Ci(x) = EulerGamma + Log[x] + Sum_{m>=1} (-1)^m x^(2m) / (2m (2m)!),
 *
 * i.e. the x^0 coefficient carries EulerGamma + Log[x] (or Log[-x] when the
 * Assumptions option forces x < 0) and only even powers appear thereafter. */
static Expr* try_series_cosintegral_at_zero(Expr* f, Expr* x, int64_t n, int x_sign) {
    if (!has_symbol_head(f, "CosIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    if (n < 0) n = 0;

    size_t ncoef = (size_t)n + 1;            /* exponents 0 .. n */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    Expr* logarg = (x_sign < 0)
                 ? mk_times(expr_new_integer(-1), expr_copy(x))
                 : expr_copy(x);
    coefs[0] = simp(mk_plus(mk_symbol("EulerGamma"), mk_fn1("Log", logarg)));
    for (int64_t k = 1; k <= n; k++) {
        if (k % 2 != 0) { coefs[k] = expr_new_integer(0); continue; }
        int64_t m = k / 2;                   /* x^(2m) coefficient (-1)^m/(2m (2m)!) */
        Expr* kfact = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* denom = eval_and_free(mk_times(expr_new_integer(k), kfact));
        Expr* recip = eval_and_free(mk_power(denom, expr_new_integer(-1)));
        coefs[k] = (m % 2 == 0) ? recip
                                : eval_and_free(mk_times(expr_new_integer(-1), recip));
    }
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);            /* expansion var */
    sd[1] = expr_new_integer(0);     /* x0            */
    sd[2] = coef_list;               /* coefficients  */
    sd[3] = expr_new_integer(0);     /* nmin          */
    sd[4] = expr_new_integer(n + 1); /* nmax (O-term) */
    sd[5] = expr_new_integer(1);     /* denominator   */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* CoshIntegral[x] at x = 0: Chi(0) = -Infinity blocks naive Taylor. Emit the
 * logarithmic series directly (the trig series without the alternating sign):
 *
 *   Chi(x) = EulerGamma + Log[x] + Sum_{m>=1} x^(2m) / (2m (2m)!),
 *
 * i.e. the x^0 coefficient carries EulerGamma + Log[x] (or Log[-x] when the
 * Assumptions option forces x < 0) and only even powers appear thereafter, all
 * with positive coefficient. */
static Expr* try_series_coshintegral_at_zero(Expr* f, Expr* x, int64_t n, int x_sign) {
    if (!has_symbol_head(f, "CoshIntegral") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;
    if (n < 0) n = 0;

    size_t ncoef = (size_t)n + 1;            /* exponents 0 .. n */
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    Expr* logarg = (x_sign < 0)
                 ? mk_times(expr_new_integer(-1), expr_copy(x))
                 : expr_copy(x);
    coefs[0] = simp(mk_plus(mk_symbol("EulerGamma"), mk_fn1("Log", logarg)));
    for (int64_t k = 1; k <= n; k++) {
        if (k % 2 != 0) { coefs[k] = expr_new_integer(0); continue; }
        /* x^(2m) coefficient 1/(2m (2m)!) -- positive, no sign alternation. */
        Expr* kfact = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* denom = eval_and_free(mk_times(expr_new_integer(k), kfact));
        coefs[k] = eval_and_free(mk_power(denom, expr_new_integer(-1)));
    }
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);            /* expansion var */
    sd[1] = expr_new_integer(0);     /* x0            */
    sd[2] = coef_list;               /* coefficients  */
    sd[3] = expr_new_integer(0);     /* nmin          */
    sd[4] = expr_new_integer(n + 1); /* nmax (O-term) */
    sd[5] = expr_new_integer(1);     /* denominator   */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Asymptotic series of BesselK[nu, x] at x = Infinity (DLMF 10.40.2):
 *
 *   K_nu(x) ~ sqrt(pi/(2x)) e^{-x} Sum_{k>=0} a_k / x^k,
 *     a_0 = 1,  a_k = a_{k-1} (4 nu^2 - (2k-1)^2) / (8 k).
 *
 * Folding sqrt(pi/(2x)) = sqrt(pi/2) x^{-1/2} in, the sum becomes a Puiseux
 * series in 1/x with half-integer powers (1/2, 3/2, 5/2, ...; denominator 2),
 * coefficient sqrt(pi/2) a_k, and the e^{-x} prefactor is kept symbolic:
 *
 *   E^{-x} SeriesData[1/x, 0, {sqrt(pi/2) a_k at index 2k}, 1, 2n+1, 2].
 *
 * Coefficients are polynomial in nu, so symbolic order works too. Returns NULL
 * unless f is exactly BesselK[nu, x] in the expansion variable. */
static Expr* try_series_besselk_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "BesselK") || f->data.function.arg_count != 2)
        return NULL;
    Expr* nu = f->data.function.args[0];
    if (!expr_eq(f->data.function.args[1], x)) return NULL;

    int64_t kmax = n - 1;                    /* (2k+1)/2 <= n  <=>  k <= n-1/2 */
    if (kmax < 0) kmax = 0;

    /* mu = 4 nu^2 (symbolic if nu is). */
    Expr* mu = simp(mk_times(expr_new_integer(4),
                             mk_power(expr_copy(nu), expr_new_integer(2))));
    Expr* sqrtpi2 = simp(mk_power(mk_times(mk_symbol("Pi"), make_rational(1, 2)),
                                  make_rational(1, 2)));   /* Sqrt[Pi/2] */

    /* a-sequence a_0..a_kmax. */
    Expr** a = calloc((size_t)kmax + 1, sizeof(Expr*));
    a[0] = expr_new_integer(1);
    for (int64_t k = 1; k <= kmax; k++) {
        int64_t odd = 2 * k - 1;
        Expr* num = mk_plus(expr_copy(mu), expr_new_integer(-(odd * odd)));
        Expr* den = mk_power(expr_new_integer(8 * k), expr_new_integer(-1));
        a[k] = simp(mk_times(expr_copy(a[k - 1]), mk_times(num, den)));
    }

    /* Coefficients: index 2k carries sqrt(pi/2) a_k; length 2*(kmax+1). */
    size_t ncoef = (size_t)(2 * kmax + 1);
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);
    for (int64_t k = 0; k <= kmax; k++) {
        expr_free(coefs[(size_t)(2 * k)]);
        coefs[(size_t)(2 * k)] = simp(mk_times(expr_copy(sqrtpi2), expr_copy(a[k])));
    }

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1));   /* 1/x */
    sd[1] = expr_new_integer(0);
    sd[2] = expr_new_function(mk_symbol("List"), coefs, ncoef);
    sd[3] = expr_new_integer(1);                            /* nmin = 1 (x^{-1/2}) */
    sd[4] = expr_new_integer(2 * kmax + 3);                 /* nmax (O-term)       */
    sd[5] = expr_new_integer(2);                            /* den                 */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd); free(coefs);

    Expr* result = mk_times(mk_fn1("Exp", mk_times(expr_new_integer(-1), expr_copy(x))),
                            series);

    expr_free(mu);
    expr_free(sqrtpi2);
    for (int64_t k = 0; k <= kmax; k++) expr_free(a[k]);
    free(a);
    return result;
}

/* Asymptotic series of BesselI[nu, x] at x = Infinity (DLMF 10.40.5):
 *
 *   I_nu(x) ~ e^x/sqrt(2 pi x) Sum_{k} (-1)^k a_k/x^k
 *           + i e^{i nu pi} e^{-x}/sqrt(2 pi x) Sum_{k} a_k/x^k,
 *     a_0 = 1, a_k = a_{k-1}(4 nu^2 - (2k-1)^2)/(8k),
 *
 * using the upper-sign (arg x = 0) branch e^{i(nu+1/2)pi} = i e^{i nu pi}.
 * Both sums fold sqrt(1/(2 pi)) x^{-1/2} into half-integer-power Puiseux series
 * in 1/x (the same shape as BesselK's), with the two exponential prefactors
 * kept symbolic. The result is a Plus of two E^... * SeriesData terms,
 * reproducing Mathematica's Normal form. Coefficients are polynomial in nu, so
 * symbolic order works. Returns NULL unless f is exactly BesselI[nu, x]. */
static Expr* try_series_besseli_at_infinity(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "BesselI") || f->data.function.arg_count != 2)
        return NULL;
    Expr* nu = f->data.function.args[0];
    if (!expr_eq(f->data.function.args[1], x)) return NULL;

    int64_t kmax = n - 1;
    if (kmax < 0) kmax = 0;

    /* mu = 4 nu^2; sqrt(1/(2 pi)) = (2 pi)^{-1/2}. */
    Expr* mu = simp(mk_times(expr_new_integer(4),
                             mk_power(expr_copy(nu), expr_new_integer(2))));
    Expr* sqrt12pi = simp(mk_power(mk_times(expr_new_integer(2), mk_symbol("Pi")),
                                   make_rational(-1, 2)));   /* 1/Sqrt[2 Pi] */

    /* a-sequence a_0..a_kmax. */
    Expr** a = calloc((size_t)kmax + 1, sizeof(Expr*));
    a[0] = expr_new_integer(1);
    for (int64_t k = 1; k <= kmax; k++) {
        int64_t odd = 2 * k - 1;
        Expr* num = mk_plus(expr_copy(mu), expr_new_integer(-(odd * odd)));
        Expr* den = mk_power(expr_new_integer(8 * k), expr_new_integer(-1));
        a[k] = simp(mk_times(expr_copy(a[k - 1]), mk_times(num, den)));
    }

    size_t ncoef = (size_t)(2 * kmax + 1);

    /* Build a SeriesData[1/x, 0, {c_k}, 1, 2kmax+3, 2] from per-k coefficients;
     * `dom` selects (-1)^k a_k sqrt(1/2pi) (dominant) vs a_k sqrt(1/2pi). */
    /* (constructed inline twice below) */
    Expr** coefs_dom = calloc(ncoef, sizeof(Expr*));
    Expr** coefs_sub = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) {
        coefs_dom[i] = expr_new_integer(0);
        coefs_sub[i] = expr_new_integer(0);
    }
    for (int64_t k = 0; k <= kmax; k++) {
        Expr* base = simp(mk_times(expr_copy(sqrt12pi), expr_copy(a[k])));
        expr_free(coefs_sub[(size_t)(2 * k)]);
        coefs_sub[(size_t)(2 * k)] = expr_copy(base);
        expr_free(coefs_dom[(size_t)(2 * k)]);
        coefs_dom[(size_t)(2 * k)] =
            (k & 1) ? simp(mk_times(expr_new_integer(-1), base))
                    : base;
    }

    Expr** sdd = calloc(6, sizeof(Expr*));
    sdd[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sdd[1] = expr_new_integer(0);
    sdd[2] = expr_new_function(mk_symbol("List"), coefs_dom, ncoef);
    sdd[3] = expr_new_integer(1);
    sdd[4] = expr_new_integer(2 * kmax + 3);
    sdd[5] = expr_new_integer(2);
    Expr* series_dom = expr_new_function(mk_symbol("SeriesData"), sdd, 6);
    free(sdd); free(coefs_dom);

    Expr** sds = calloc(6, sizeof(Expr*));
    sds[0] = mk_power(expr_copy(x), expr_new_integer(-1));
    sds[1] = expr_new_integer(0);
    sds[2] = expr_new_function(mk_symbol("List"), coefs_sub, ncoef);
    sds[3] = expr_new_integer(1);
    sds[4] = expr_new_integer(2 * kmax + 3);
    sds[5] = expr_new_integer(2);
    Expr* series_sub = expr_new_function(mk_symbol("SeriesData"), sds, 6);
    free(sds); free(coefs_sub);

    /* Dominant: E^x * series_dom. */
    Expr* dom = mk_times(mk_fn1("Exp", expr_copy(x)), series_dom);

    /* Subdominant: E^{I nu pi - x} * I * series_sub. */
    Expr* exparg = mk_plus(
        mk_times(mk_symbol("I"), mk_times(expr_copy(nu), mk_symbol("Pi"))),
        mk_times(expr_new_integer(-1), expr_copy(x)));
    Expr* sub = mk_times(mk_fn1("Exp", exparg),
                         mk_times(mk_symbol("I"), series_sub));

    Expr* result = mk_plus(dom, sub);

    expr_free(mu);
    expr_free(sqrt12pi);
    for (int64_t k = 0; k <= kmax; k++) expr_free(a[k]);
    free(a);
    return result;
}

/* HarmonicNumber H_m = Sum_{j=1}^{m} 1/j, as an exact rational Expr. */
static Expr* bk_harmonic(long m) {
    if (m <= 0) return expr_new_integer(0);
    Expr** terms = calloc((size_t)m, sizeof(Expr*));
    for (long j = 1; j <= m; j++)
        terms[(size_t)(j - 1)] = mk_power(expr_new_integer(j), expr_new_integer(-1));
    Expr* sum = expr_new_function(mk_symbol("Plus"), terms, (size_t)m);
    free(terms);
    return eval_and_free(sum);
}

/* Logarithmic series of BesselK[p, x] at x = 0 for integer order p (DLMF
 * 10.31.1). K_p is even in p (K_{-p}=K_p), so the caller's |p| is used. The
 * coefficient of x^{p+2k} carries the branch term Log[x] and EulerGamma:
 *
 *   coef(x^{p+2k}) = (-1)^p (-2 EulerGamma + 2 Log[2] - 2 Log[x]
 *                            + H_k + H_{p+k}) / (2^{p+2k+1} k! (p+k)!),
 *
 * plus, for p >= 1, the finite principal part at x^{-p+2k}, k = 0..p-1:
 *
 *   coef(x^{-p+2k}) = (-1)^k (p-k-1)! 2^{p-1-2k} / k!.
 *
 * For p = 0 this reproduces Mathematica's Series[BesselK[0,x],{x,0,n}]
 * exactly: x^0 -> -EulerGamma+Log[2]-Log[x], x^2 -> (1-EulerGamma+Log[2]-
 * Log[x])/4, ... Returns NULL unless f is BesselK[<integer>, x] in the
 * expansion variable. Non-integer order (a pure x^nu Puiseux series, with no
 * logarithm) is left to the elementary half-integer rewrites / generic path. */
static Expr* try_series_besselk_at_zero(Expr* f, Expr* x, int64_t n) {
    if (!has_symbol_head(f, "BesselK") || f->data.function.arg_count != 2)
        return NULL;
    if (!expr_eq(f->data.function.args[1], x)) return NULL;
    Expr* ord = f->data.function.args[0];
    if (ord->type != EXPR_INTEGER) return NULL;       /* integer order only */
    long p = (long)ord->data.integer;
    if (p < 0) p = -p;                                 /* K_{-p} = K_p */
    if (n < 0) n = 0;

    int64_t nmin = (p == 0) ? 0 : -p;
    int64_t nmax = n + 1;
    if (nmax <= nmin) nmax = nmin + 1;
    size_t ncoef = (size_t)(nmax - nmin);
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    long sign_p = (p & 1) ? -1 : 1;

    /* Finite principal part (p >= 1): power -p+2k, k = 0..p-1. */
    for (long k = 0; k <= p - 1; k++) {
        int64_t power = -p + 2 * k;
        if (power < nmin || power >= nmax) continue;
        Expr* fk = eval_and_free(mk_fn1("Factorial", expr_new_integer(p - k - 1)));
        Expr* kf = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* pow2 = eval_and_free(mk_power(expr_new_integer(2),
                                            expr_new_integer(p - 1 - 2 * k)));
        Expr* c = mk_times(expr_new_integer((k & 1) ? -1 : 1),
                           mk_times(fk, mk_times(pow2,
                                    mk_power(kf, expr_new_integer(-1)))));
        size_t idx = (size_t)(power - nmin);
        expr_free(coefs[idx]);
        coefs[idx] = simp(c);
    }

    /* Main (log + regular) part: power p+2k, k = 0,1,... while p+2k < nmax. */
    for (long k = 0; (int64_t)(p + 2 * k) < nmax; k++) {
        int64_t power = p + 2 * k;
        /* numerator: -2 EulerGamma + 2 Log[2] - 2 Log[x] + H_k + H_{p+k}. */
        Expr* num = mk_plus(
            mk_times(expr_new_integer(-2), mk_symbol("EulerGamma")),
            mk_plus(
                mk_times(expr_new_integer(2), mk_fn1("Log", expr_new_integer(2))),
                mk_plus(
                    mk_times(expr_new_integer(-2), mk_fn1("Log", expr_copy(x))),
                    mk_plus(bk_harmonic(k), bk_harmonic(p + k)))));
        /* denominator: 2^{p+2k+1} k! (p+k)!. */
        Expr* pow2 = eval_and_free(mk_power(expr_new_integer(2),
                                            expr_new_integer(power + 1)));
        Expr* kf = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* pkf = eval_and_free(mk_fn1("Factorial", expr_new_integer(p + k)));
        Expr* den = eval_and_free(mk_times(pow2, mk_times(kf, pkf)));
        Expr* c = mk_times(expr_new_integer(sign_p),
                           mk_times(num, mk_power(den, expr_new_integer(-1))));
        size_t idx = (size_t)(power - nmin);
        expr_free(coefs[idx]);
        coefs[idx] = simp(c);
    }

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);
    sd[1] = expr_new_integer(0);
    sd[2] = expr_new_function(mk_symbol("List"), coefs, ncoef);
    sd[3] = expr_new_integer(nmin);
    sd[4] = expr_new_integer(nmax);
    sd[5] = expr_new_integer(1);
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd); free(coefs);
    return series;
}

/* Logarithmic series of BesselY[p, x] at x = 0 for integer order p (DLMF
 * 10.8.1). Y is odd under order reflection (Y_{-p} = (-1)^p Y_p), so the
 * caller's |p| is used and the (-1)^p sign is folded into the coefficients via
 * the reflection in src/internal/bessel.m for negative p; here p >= 0.
 *
 *   coef(x^{p+2k}) = (-1)^k (2 Log[x] - 2 Log[2] + 2 EulerGamma - H_k - H_{p+k})
 *                    / (Pi 2^{p+2k} k! (p+k)!),
 *
 * plus, for p >= 1, the finite principal part at x^{-p+2k}, k = 0..p-1:
 *
 *   coef(x^{-p+2k}) = -(p-k-1)! 2^{p-2k} / (Pi k!).
 *
 * For p = 0 this reproduces Mathematica's Series[BesselY[0,x],{x,0,n}]:
 * x^0 -> 2(EulerGamma-Log[2]+Log[x])/Pi, x^2 -> (1-EulerGamma+Log[2]-Log[x])
 * /(2 Pi), ... Returns NULL unless f is BesselY[<integer>, x] in the expansion
 * variable. Non-integer order (a pure Puiseux series, no logarithm) is left to
 * the elementary half-integer rewrites / generic path. */
static Expr* try_series_bessely_at_zero(Expr* f, Expr* x, int64_t n) {
    if (!has_symbol_head(f, "BesselY") || f->data.function.arg_count != 2)
        return NULL;
    if (!expr_eq(f->data.function.args[1], x)) return NULL;
    Expr* ord = f->data.function.args[0];
    if (ord->type != EXPR_INTEGER) return NULL;       /* integer order only */
    long p = (long)ord->data.integer;
    long psign = 1;                                    /* Y_{-p} = (-1)^p Y_p */
    if (p < 0) { if (p & 1) psign = -1; p = -p; }
    if (n < 0) n = 0;

    int64_t nmin = (p == 0) ? 0 : -p;
    int64_t nmax = n + 1;
    if (nmax <= nmin) nmax = nmin + 1;
    size_t ncoef = (size_t)(nmax - nmin);
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (size_t i = 0; i < ncoef; i++) coefs[i] = expr_new_integer(0);

    /* Finite principal part (p >= 1): power -p+2k, k = 0..p-1.
     * coef = -(p-k-1)! 2^{p-2k} / (Pi k!). */
    for (long k = 0; k <= p - 1; k++) {
        int64_t power = -p + 2 * k;
        if (power < nmin || power >= nmax) continue;
        Expr* fk = eval_and_free(mk_fn1("Factorial", expr_new_integer(p - k - 1)));
        Expr* kf = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* pow2 = eval_and_free(mk_power(expr_new_integer(2),
                                            expr_new_integer(p - 1 - 2 * k)));
        /* 2^{p-2k} = 2 * 2^{p-1-2k}; absorb the leading 2 and the 1/2... keep
         * pow2 = 2^{p-1-2k} then multiply by 2 below to get 2^{p-2k}. */
        Expr* c = mk_times(expr_new_integer(-2),
                    mk_times(fk,
                      mk_times(pow2,
                        mk_times(mk_power(kf, expr_new_integer(-1)),
                                 mk_power(mk_symbol("Pi"), expr_new_integer(-1))))));
        size_t idx = (size_t)(power - nmin);
        expr_free(coefs[idx]);
        coefs[idx] = simp(c);
    }

    /* Main (log + regular) part: power p+2k, k = 0,1,... while p+2k < nmax.
     * coef = (-1)^k psign (2 Log[x] - 2 Log[2] + 2 EulerGamma - H_k - H_{p+k})
     *        / (Pi 2^{p+2k} k! (p+k)!). */
    for (long k = 0; (int64_t)(p + 2 * k) < nmax; k++) {
        int64_t power = p + 2 * k;
        /* numerator: 2 Log[x] - 2 Log[2] + 2 EulerGamma - H_k - H_{p+k}. */
        Expr* num = mk_plus(
            mk_times(expr_new_integer(2), mk_fn1("Log", expr_copy(x))),
            mk_plus(
                mk_times(expr_new_integer(-2), mk_fn1("Log", expr_new_integer(2))),
                mk_plus(
                    mk_times(expr_new_integer(2), mk_symbol("EulerGamma")),
                    mk_plus(mk_times(expr_new_integer(-1), bk_harmonic(k)),
                            mk_times(expr_new_integer(-1), bk_harmonic(p + k))))));
        /* denominator: Pi 2^{p+2k} k! (p+k)!. */
        Expr* pow2 = eval_and_free(mk_power(expr_new_integer(2),
                                            expr_new_integer(power)));
        Expr* kf = eval_and_free(mk_fn1("Factorial", expr_new_integer(k)));
        Expr* pkf = eval_and_free(mk_fn1("Factorial", expr_new_integer(p + k)));
        Expr* den = eval_and_free(mk_times(mk_symbol("Pi"),
                                    mk_times(pow2, mk_times(kf, pkf))));
        Expr* c = mk_times(expr_new_integer((k & 1) ? -psign : psign),
                           mk_times(num, mk_power(den, expr_new_integer(-1))));
        size_t idx = (size_t)(power - nmin);
        expr_free(coefs[idx]);
        coefs[idx] = simp(c);
    }

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);
    sd[1] = expr_new_integer(0);
    sd[2] = expr_new_function(mk_symbol("List"), coefs, ncoef);
    sd[3] = expr_new_integer(nmin);
    sd[4] = expr_new_integer(nmax);
    sd[5] = expr_new_integer(1);
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd); free(coefs);
    return series;
}

/* LerchPhi[x, s, a] at x = 0: the defining series Sum_{k>=0} x^k (k+a)^-s is
 * already a power series, so the coefficient of x^k is the symmetric power
 * ((k+a)^2)^(-s/2) that Mathematica reports.  Emit it directly: naive
 * Taylor-via-D fails because D[LerchPhi[x,s,a],x] =
 * (LerchPhi[x,s-1,a] - a LerchPhi[x,s,a])/x is a 0/0 form at x = 0.  Only fires
 * when the expansion variable is the first argument and s, a are still symbolic
 * (a reducible a -- positive/non-positive integer -- has already collapsed onto
 * PolyLog before this point, and that form expands fine). */
static Expr* try_series_lerchphi_at_zero(Expr* f, Expr* x, int64_t n) {
    if (!f || f->type != EXPR_FUNCTION) return NULL;
    if (!has_symbol_head(f, "LerchPhi")) return NULL;
    if (f->data.function.arg_count != 3) return NULL;
    Expr* z = f->data.function.args[0];
    Expr* s = f->data.function.args[1];
    Expr* a = f->data.function.args[2];
    if (!expr_eq(z, x)) return NULL;
    if (n < 0) return NULL;

    size_t ncoef = (size_t)n + 1;
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (int64_t k = 0; k <= n; k++) {
        /* ((k + a)^2)^(-s/2) */
        Expr* kpa  = mk_plus(expr_new_integer(k), expr_copy(a));
        Expr* sq   = mk_power(kpa, expr_new_integer(2));
        Expr* hns  = mk_times(expr_copy(s),
                              mk_power(expr_new_integer(2), expr_new_integer(-1)));
        Expr* c    = mk_power(sq, mk_times(expr_new_integer(-1), hns));
        coefs[(size_t)k] = simp(c);
    }

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);
    sd[1] = expr_new_integer(0);
    sd[2] = expr_new_function(mk_symbol("List"), coefs, ncoef);
    sd[3] = expr_new_integer(0);
    sd[4] = expr_new_integer(n + 1);
    sd[5] = expr_new_integer(1);
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd); free(coefs);
    return series;
}

/* ----------------------------------------------------------------------------
 * ProductLog (Lambert W) series expansions.
 * ------------------------------------------------------------------------- */

/* Taylor series of ProductLog[x] at x = 0:
 *   W(x) = Sum_{k>=1} (-k)^(k-1)/k! x^k
 *        = x - x^2 + 3/2 x^3 - 8/3 x^4 + 125/24 x^5 - ...
 * The generic Taylor-via-D path cannot produce this (W'(0) is a 0/0 form), so
 * the closed-form coefficients are emitted directly. Principal branch only. */
static Expr* try_series_productlog_at_zero(Expr* f, Expr* x, int64_t n) {
    if (n < 1) n = 1;
    if (!has_symbol_head(f, "ProductLog") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    /* a_k = (-k)^(k-1)/k! at x-power k, k = 1 .. n. */
    Expr** coefs = calloc((size_t)n, sizeof(Expr*));
    for (int64_t k = 1; k <= n; k++) {
        Expr* num = mk_power(expr_new_integer(-k), expr_new_integer(k - 1));
        Expr* den = mk_power(mk_fn1("Factorial", expr_new_integer(k)),
                             expr_new_integer(-1));
        coefs[k - 1] = eval_and_free(mk_times(num, den));
    }
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, (size_t)n);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);                /* expansion variable */
    sd[1] = expr_new_integer(0);         /* x0                 */
    sd[2] = coef_list;                   /* coefficients       */
    sd[3] = expr_new_integer(1);         /* nmin (x^1)         */
    sd[4] = expr_new_integer(n + 1);     /* nmax (O-term)      */
    sd[5] = expr_new_integer(1);         /* denominator        */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Multiply two truncated power series a[0..D], b[0..D] (coefficients are owned
 * Expr*), returning a freshly allocated array c[0..D] with c_p = Sum a_i b_{p-i}
 * simplified. Inputs are left intact (the caller frees them). */
static Expr** pl_poly_mul(Expr** a, Expr** b, int64_t D) {
    Expr** c = calloc((size_t)(D + 1), sizeof(Expr*));
    for (int64_t p = 0; p <= D; p++) {
        Expr* acc = expr_new_integer(0);
        for (int64_t i = 0; i <= p; i++) {
            Expr* term = mk_times(expr_copy(a[i]), expr_copy(b[p - i]));
            acc = eval_and_free(mk_plus(acc, term));
        }
        c[p] = acc;
    }
    return c;
}

/* Compute the rational branch-point coefficients mu_0..mu_K of
 *   W(z) = Sum_{k>=0} mu_k p^k,   p = Sqrt[2(e z + 1)],   z near -1/e,
 * by reverting  p^2/2 = Sum_{j>=2} (j-1)/j! v^j  with v = Sum_{k>=1} mu_k p^k.
 * Returns an array of K+1 owned Expr* (mu_0 = -1, mu_1 = 1, ...). */
static Expr** pl_branch_mu(int64_t K) {
    Expr** mu = calloc((size_t)(K + 1), sizeof(Expr*));
    mu[0] = expr_new_integer(-1);
    if (K >= 1) mu[1] = expr_new_integer(1);

    for (int64_t i = 2; i <= K; i++) {
        int64_t J = i + 1;                       /* solve at p-power J */
        /* v[0..J] from already-known mu_1..mu_{i-1} (mu_i, ... treated as 0). */
        Expr** v = calloc((size_t)(J + 1), sizeof(Expr*));
        for (int64_t t = 0; t <= J; t++)
            v[t] = (t >= 1 && t <= i - 1) ? expr_copy(mu[t]) : expr_new_integer(0);

        /* acc = [p^J] Sum_{j=2}^{J} c_j v^j,  c_j = (j-1)/j!. */
        Expr* acc = expr_new_integer(0);
        Expr** vpow = calloc((size_t)(J + 1), sizeof(Expr*)); /* v^1 */
        for (int64_t t = 0; t <= J; t++) vpow[t] = expr_copy(v[t]);
        for (int64_t j = 2; j <= J; j++) {
            Expr** nxt = pl_poly_mul(vpow, v, J);            /* v^j */
            for (int64_t t = 0; t <= J; t++) expr_free(vpow[t]);
            free(vpow);
            vpow = nxt;
            Expr* cj = eval_and_free(mk_times(expr_new_integer(j - 1),
                          mk_power(mk_fn1("Factorial", expr_new_integer(j)),
                                   expr_new_integer(-1))));
            Expr* contrib = mk_times(cj, expr_copy(vpow[J]));
            acc = eval_and_free(mk_plus(acc, contrib));
        }
        for (int64_t t = 0; t <= J; t++) expr_free(vpow[t]);
        free(vpow);
        for (int64_t t = 0; t <= J; t++) expr_free(v[t]);
        free(v);

        /* mu_i = -[p^J](...). */
        mu[i] = eval_and_free(mk_times(expr_new_integer(-1), acc));
    }
    return mu;
}

/* Puiseux series of ProductLog[x] at the branch point x = -1/E:
 *   W(x) = Sum_{k>=0} mu_k (2E)^(k/2) (x + 1/E)^(k/2)
 *        = -1 + Sqrt[2E] Sqrt[x+1/E] - 2/3 E (x+1/E) + ...
 * Half-integer powers (denominator 2). The O-term convention here keeps powers
 * through (x+1/E)^n (i.e. numerator 2n); it is mathematically a valid degree-n
 * truncation and need not match Mathematica's term count exactly. Principal
 * branch only. */
static Expr* try_series_productlog_at_branchpoint(Expr* f, Expr* x, int64_t n) {
    if (n < 0) n = 0;
    if (!has_symbol_head(f, "ProductLog") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    int64_t K = 2 * n;                          /* highest half-integer numerator */
    if (K < 1) K = 1;
    Expr** mu = pl_branch_mu(K);

    /* coef numerator p -> mu_p (2E)^(p/2) at (x+1/E)^(p/2). */
    size_t ncoef = (size_t)(K + 1);
    Expr** coefs = calloc(ncoef, sizeof(Expr*));
    for (int64_t p = 0; p <= K; p++) {
        Expr* two_e = mk_times(expr_new_integer(2), mk_symbol("E"));
        Expr* fac   = mk_power(two_e, make_rational(p, 2));   /* (2E)^(p/2) */
        coefs[(size_t)p] = eval_and_free(mk_times(expr_copy(mu[p]), fac));
    }
    for (int64_t p = 0; p <= K; p++) expr_free(mu[p]);
    free(mu);

    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, ncoef);
    free(coefs);

    /* x0 = -1/E; the printer forms powers of (x - x0) = x + 1/E. */
    Expr* x0 = mk_times(expr_new_integer(-1),
                        mk_power(mk_symbol("E"), expr_new_integer(-1)));

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = expr_copy(x);                /* expansion variable x       */
    sd[1] = eval_and_free(x0);           /* x0 = -1/E                  */
    sd[2] = coef_list;
    sd[3] = expr_new_integer(0);         /* nmin (numerator)           */
    sd[4] = expr_new_integer(K + 1);     /* nmax (O-term numerator)    */
    sd[5] = expr_new_integer(2);         /* denominator (half powers)  */
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Unsigned Stirling numbers of the first kind c1(m, k) (cycle counts), small m. */
static int64_t pl_stirling1u(int64_t m, int64_t k) {
    if (k < 0 || k > m) return 0;
    if (m == 0) return (k == 0) ? 1 : 0;
    /* c1(m,k) = c1(m-1,k-1) + (m-1) c1(m-1,k). */
    return pl_stirling1u(m - 1, k - 1) + (m - 1) * pl_stirling1u(m - 1, k);
}

/* Asymptotic (nested-logarithm) expansion of ProductLog[x] at x = Infinity:
 *   W(x) ~ L1 - L2 + Sum_{k>=0,m>=1} ((-1)^k/m!) c1(k+m, k+1) L2^m / L1^(k+m),
 *   L1 = Log[x],  L2 = Log[Log[x]].
 * This is the coefficient of x^0; it is not a power series in 1/x, so it is
 * emitted as the single x^0 coefficient of a SeriesData[1/x, ...] with an
 * O[1/x]^1 term (matching Mathematica's order-0 output, which keeps 1/L1 powers
 * through 2). Higher 1/x corrections are not produced. */
static Expr* try_series_productlog_at_infinity(Expr* f, Expr* x, int64_t n) {
    (void)n;
    if (!has_symbol_head(f, "ProductLog") || f->data.function.arg_count != 1)
        return NULL;
    if (!expr_eq(f->data.function.args[0], x)) return NULL;

    Expr* L1 = mk_fn1("Log", expr_copy(x));
    Expr* L2 = mk_fn1("Log", mk_fn1("Log", expr_copy(x)));

    /* block = L1 - L2 + double sum (k+m from 1 to DEPTH). */
    const int64_t DEPTH = 2;            /* max power of 1/L1 (matches WL order 0) */
    Expr* block = mk_plus(expr_copy(L1),
                          mk_times(expr_new_integer(-1), expr_copy(L2)));
    for (int64_t d = 1; d <= DEPTH; d++) {           /* d = k + m = 1/L1 power */
        for (int64_t m = 1; m <= d; m++) {
            int64_t k = d - m;
            int64_t s1 = pl_stirling1u(k + m, k + 1);
            if (s1 == 0) continue;
            /* coeff = (-1)^k c1(k+m,k+1) / m!  (a rational). */
            Expr* rat = eval_and_free(mk_times(
                expr_new_integer(((k % 2) == 0) ? s1 : -s1),
                mk_power(mk_fn1("Factorial", expr_new_integer(m)),
                         expr_new_integer(-1))));
            Expr* term = mk_times(rat,
                mk_times(mk_power(expr_copy(L2), expr_new_integer(m)),
                         mk_power(expr_copy(L1), expr_new_integer(-d))));
            block = mk_plus(block, term);
        }
    }
    expr_free(L1);
    expr_free(L2);
    block = eval_and_free(block);

    Expr** coefs = calloc(1, sizeof(Expr*));
    coefs[0] = block;                                /* coefficient of (1/x)^0 */
    Expr* coef_list = expr_new_function(mk_symbol("List"), coefs, 1);
    free(coefs);

    Expr** sd = calloc(6, sizeof(Expr*));
    sd[0] = mk_power(expr_copy(x), expr_new_integer(-1));  /* expansion var 1/x */
    sd[1] = expr_new_integer(0);
    sd[2] = coef_list;
    sd[3] = expr_new_integer(0);         /* nmin (x^0)   */
    sd[4] = expr_new_integer(1);         /* nmax -> O[1/x]^1 */
    sd[5] = expr_new_integer(1);
    Expr* series = expr_new_function(mk_symbol("SeriesData"), sd, 6);
    free(sd);
    return series;
}

/* Expand f around x=x0 to order n and return SeriesData expression. */
static Expr* do_series_single(Expr* f, Expr* x, Expr* x0, int64_t n, bool leading_only,
                              int x_sign) {
    /* Evaluate f with the series context implicit. Since Series has
     * HoldAll, f has not been evaluated yet; we evaluate now. */
    Expr* f_eval = eval_and_free(expr_copy(f));
    Expr* x0_eval = eval_and_free(expr_copy(x0));

    /* Mathematica convention: Series[c, {x, ...}] where c is free of x
     * returns c verbatim (not SeriesData[...]). Applies to Series[0, ...]
     * and Series[Sin[y], {x, 0, 4}] alike. */
    if (expr_free_of(f_eval, x)) {
        expr_free(x0_eval);
        return f_eval;
    }

    /* Special functions with essential singularities at Infinity (e.g.
     * ExpIntegralEi) have no Laurent series there -- the generic x -> 1/u
     * substitution would hand a pole to naive Taylor. Emit their known
     * asymptotic expansions (with the E^x prefactor kept symbolic) directly. */
    if (x0_eval->type == EXPR_SYMBOL && x0_eval->data.symbol.name == SYM_Infinity) {
        Expr* ei = try_series_ei_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (ei) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return ei;
        }
        Expr* si = try_series_sinintegral_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (si) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return si;
        }
        Expr* ci = try_series_cosintegral_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (ci) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return ci;
        }
        Expr* shi = try_series_sinhintegral_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (shi) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return shi;
        }
        Expr* chi = try_series_coshintegral_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (chi) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return chi;
        }
        Expr* frc = try_series_fresnelc_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (frc) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return frc;
        }
        Expr* frs = try_series_fresnels_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (frs) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return frs;
        }
        Expr* ai = try_series_airyai_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (ai) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return ai;
        }
        Expr* bi = try_series_airybi_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (bi) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return bi;
        }
        Expr* aip = try_series_airyaiprime_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (aip) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return aip;
        }
        Expr* bip = try_series_airybiprime_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (bip) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return bip;
        }
        Expr* bj = try_series_besselj_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (bj) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return bj;
        }
        Expr* bk = try_series_besselk_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (bk) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return bk;
        }
        Expr* besi = try_series_besseli_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (besi) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return besi;
        }
        Expr* by = try_series_bessely_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (by) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return by;
        }
        /* Error functions at Infinity: DLMF 7.12.1 asymptotic expansions,
         * each an Exp[±x^2] essential singularity times a Laurent series. */
        Expr* erf = try_series_erf_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (erf) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return erf;
        }
        Expr* erfc = try_series_erfc_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (erfc) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return erfc;
        }
        Expr* erfi = try_series_erfi_at_infinity(f_eval, x, leading_only ? 1 : n);
        if (erfi) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return erfi;
        }
        /* ProductLog[x] at Infinity: nested-logarithm asymptotic expansion
         * (the x^0 coefficient), not a power series in 1/x. */
        Expr* plw = try_series_productlog_at_infinity(f_eval, x, leading_only ? 0 : n);
        if (plw) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return plw;
        }
    }

    /* LogIntegral[x] = Ei(Log[x]) has no Taylor series at x = 0 -- Log[x]
     * drives it into Ei's asymptotic regime, yielding a generalized series
     * whose x^1 coefficient is a series in 1/Log[x]. Emit it directly so
     * naive Taylor-via-D does not hit (and leak) the 1/0 pole of the second
     * derivative. */
    if (is_lit_zero(x0_eval)) {
        Expr* li = try_series_logintegral_at_zero(f_eval, x, leading_only ? 1 : n, x_sign);
        if (li) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return li;
        }
        /* ExpIntegralEi[x] at x = 0: Ei(0) = -Infinity blocks naive Taylor;
         * emit the DLMF 6.6.2 logarithmic series directly. */
        Expr* ei0 = try_series_ei_at_zero(f_eval, x, leading_only ? 0 : n, x_sign);
        if (ei0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return ei0;
        }
        /* CosIntegral[x] at x = 0: Ci(0) = -Infinity blocks naive Taylor;
         * emit the DLMF 6.6.2 logarithmic series directly. */
        Expr* ci0 = try_series_cosintegral_at_zero(f_eval, x, leading_only ? 0 : n, x_sign);
        if (ci0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return ci0;
        }
        /* CoshIntegral[x] at x = 0: Chi(0) = -Infinity blocks naive Taylor;
         * emit the logarithmic series directly. */
        Expr* chi0 = try_series_coshintegral_at_zero(f_eval, x, leading_only ? 0 : n, x_sign);
        if (chi0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return chi0;
        }
        /* BesselK[p, x] at x = 0 (integer p): K_p has a pole and a logarithmic
         * branch term, so naive Taylor-via-D cannot produce it. Emit the
         * DLMF 10.31.1 logarithmic series directly. */
        Expr* bk0 = try_series_besselk_at_zero(f_eval, x, leading_only ? 0 : n);
        if (bk0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return bk0;
        }
        /* BesselY[p, x] at x = 0 (integer p): logarithmic singularity (DLMF
         * 10.8.1), so naive Taylor-via-D cannot produce it. */
        Expr* by0 = try_series_bessely_at_zero(f_eval, x, leading_only ? 0 : n);
        if (by0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return by0;
        }
        /* LerchPhi[x, s, a] at x = 0 (symbolic s, a): the defining power series,
         * coefficient ((k+a)^2)^(-s/2). Naive Taylor-via-D hits a 0/0 form. */
        Expr* lp0 = try_series_lerchphi_at_zero(f_eval, x, leading_only ? 0 : n);
        if (lp0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return lp0;
        }
        /* ProductLog[x] at x = 0: closed-form Taylor coefficients (-k)^(k-1)/k!
         * (W'(0) is a 0/0 form, so naive Taylor-via-D fails). */
        Expr* plw0 = try_series_productlog_at_zero(f_eval, x, leading_only ? 1 : n);
        if (plw0) {
            expr_free(f_eval);
            expr_free(x0_eval);
            return plw0;
        }
    }

    /* ProductLog[x] at the branch point x = -1/E: Puiseux series in
     * Sqrt[x + 1/E] (W'(-1/e) = Infinity blocks naive Taylor-via-D). */
    {
        Expr* neg_inv_e = eval_and_free(
            mk_times(expr_new_integer(-1),
                     mk_power(mk_symbol("E"), expr_new_integer(-1))));
        bool at_bp = expr_eq(neg_inv_e, x0_eval);
        expr_free(neg_inv_e);
        if (at_bp) {
            Expr* plbp = try_series_productlog_at_branchpoint(
                f_eval, x, leading_only ? 0 : n);
            if (plbp) {
                expr_free(f_eval);
                expr_free(x0_eval);
                return plbp;
            }
        }
    }

    /* Apart preprocessing: if f_eval is (or contains) a rational function
     * in x, decompose it into partial fractions first. Each term in the
     * decomposition is cheaper to expand (c/(a + b*x)^m -> monomial fast
     * path) and avoids the O(N*deg(q)) Newton inversion of the composite
     * denominator. Gated by a cheap check so non-rational inputs pay no
     * overhead. */
    {
        Expr* apart_result = try_apart_preprocess(f_eval, x);
        if (apart_result) {
            expr_free(f_eval);
            f_eval = apart_result;
            /* Re-check the "free of x" early-out after Apart; a rational
             * function whose simplification collapses to a constant would
             * otherwise take the generic path unnecessarily. */
            if (expr_free_of(f_eval, x)) {
                expr_free(x0_eval);
                return f_eval;
            }
        }
    }

    /* Factor out x^alpha for symbolic alpha (e.g. Series[x^a Exp[x], ...]).
     * We detect only the top-level Times case; nested occurrences are rare
     * and harder to match Mathematica's output form on. */
    /* Factor out x^alpha only when the expansion point is 0 or Infinity.
     * At x0 = 0, `x^alpha` literally is the prefactor (Laurent monomial
     * with fractional exponent), and at Infinity we substitute x -> 1/u
     * which turns x^alpha into u^(-alpha), still a clean prefactor. At
     * any other x0, `x^alpha = (x0 + (x-x0))^alpha = x0^alpha (1 + u/x0)^alpha`
     * which the binomial kernel can expand exactly, so pulling `x^alpha`
     * out as a symbolic prefactor would SKIP that expansion and leave
     * trailing factors (e.g. `(x-1)^-2`) computed at the wrong point. */
    Expr* prefactor = NULL;
    Expr* f_used_for_expand = f_eval;
    Expr* f_body_owned = NULL;
    bool at_zero_or_infinity =
        is_lit_zero(x0_eval) ||
        (x0_eval->type == EXPR_SYMBOL &&
         x0_eval->data.symbol.name == SYM_Infinity);
    if (at_zero_or_infinity) {
        Expr* body = NULL;
        Expr* pf = try_factor_power_prefactor(f_eval, x, &body);
        if (pf) {
            prefactor = pf;
            f_body_owned = body;
            f_used_for_expand = body;
        }
    }

    /* Handle expansion at Infinity by substituting x -> 1/u. */
    bool at_infinity = (x0_eval->type == EXPR_SYMBOL &&
                        (x0_eval->data.symbol.name == SYM_Infinity));
    Expr* f_use  = f_used_for_expand;
    Expr* x_use  = x;
    Expr* x0_use = x0_eval;
    Expr* u_sym  = NULL;
    Expr* f_sub_owned = NULL;
    if (at_infinity) {
        u_sym = mk_symbol("$SeriesInvVar$");
        Expr* inv_u = mk_power(expr_copy(u_sym), expr_new_integer(-1));
        f_sub_owned = replace_all_of(f_used_for_expand, x, inv_u);
        expr_free(inv_u);
        f_use = f_sub_owned;
        x_use = u_sym;
        expr_free(x0_eval);
        x0_use = expr_new_integer(0);
    }

    /* target_exp_n is the user-facing O-term exponent, i.e. x^target_exp_n
     * is the smallest power that has been dropped. For Series[f, {x, x0, n}]
     * this is n. For the leading-term form Series[f, x -> x0] we set it to
     * 0 initially, then widen it by scanning the computed series so the
     * O-term lands at the next non-zero exponent (matching Mathematica's
     * extension-on-cancellation behaviour, e.g. Sin[x]-x -> -x^3/6 + O[x]^5). */
    int64_t target_exp_n = leading_only ? 0 : n;
    int64_t order = target_exp_n + 1;
    if (order < 1) order = 1;

    /* Laurent and compound expressions lose O-term accuracy when negative
     * exponents flow through multiplication; we compensate by expanding
     * internally to a padded order, then truncating the result back down
     * to the user-requested O-term. Padding = 12 is chosen to survive
     * deep Laurent expansions like 1/Sin[x]^10 while keeping big-rational
     * coefficient arithmetic (Puiseux, logarithmic) tractable.
     *
     * When x0 is not a literal number, coefficients are symbolic
     * (e.g. Sinh[a], Cosh[a]) and so_inv's O(N^2) symbolic convolution
     * explodes expression size via evaluate(): e.g. Series[Coth[x],
     * {x, a, 1}] with pad=12 would grow unbounded. Cap pad tight in that
     * case -- real poles/roots at a symbolic point are uncommon, and
     * compound expressions still get a couple of correction slots. */
    /* Complex[num, num] (every upper-half-plane pole) and Rational are concrete
     * numbers: so_inv keeps their coefficients numeric and bounded, so they earn
     * the full Laurent pad. Restricting to INTEGER/REAL/BIGINT wrongly demoted a
     * complex pole like I = Complex[0,1] to pad=2, one term short -- the residue
     * at a double/triple complex pole (Fourier/Jordan family) then dropped its
     * product-rule cross term. Only a free/symbolic x0 keeps the tight pad. */
    bool x0_is_numeric = expr_is_numeric_like(x0_use);
    int64_t pad = x0_is_numeric ? 12 : 2;
    int64_t internal_order = order + pad;

    SeriesCtx ctx = {
        x_use, x0_use, internal_order, order,
        /* allow_branch_wrap */ true,
        /* pending_add_const  */ NULL,
        /* pending_log_coef   */ NULL,
        /* pending_discriminator */ NULL,
    };
    SeriesObj* s = series_expand(f_use, &ctx);
    Expr* result = NULL;
    if (s) {
        /* Scan for leading-term extension in the rule form: find the first
         * non-zero coefficient at exponent e1, then the next at e2, and
         * place the O-term at e2 (or e1+1 if no next term within the
         * internal expansion, or at the original target otherwise). */
        if (leading_only) {
            int64_t e1 = so_first_nonzero_exp(s, 0);
            int64_t new_target = 1;
            if (e1 != INT64_MAX) {
                int64_t e2 = so_first_nonzero_exp(s, e1 + 1);
                new_target = (e2 != INT64_MAX) ? e2 : (e1 + 1);
                if (new_target < 1) new_target = 1;
            }
            target_exp_n = new_target;
        }
        /* Truncate to the user-facing O-term exponent.
         * User asked for x^target_exp_n to be the first dropped term.
         * In the series's internal denominator this is target_exp_n*den
         * plus one step: terms at exponents strictly less than
         * target_exp_n + 1/den are kept, and the O-term is placed at
         * the next valid step beyond target_exp_n. */
        int64_t target_num;
        if (leading_only) {
            target_num = target_exp_n * s->den;  /* O-term exactly at target_exp_n. */
        } else {
            target_num = target_exp_n * s->den + 1;
        }
        if (s->order > target_num) {
            SeriesObj* t = so_alloc(s->x, s->x0, s->nmin, target_num, s->den);
            for (size_t i = 0; i < t->coef_count && i < s->coef_count; i++) {
                so_set_coef(t, i, expr_copy(s->coefs[i]));
            }
            so_free(s);
            s = t;
        }
        if (at_infinity) {
            Expr* inv_x = mk_power(expr_copy(x), expr_new_integer(-1));
            expr_free(s->x);
            s->x = inv_x;
            /* Coefficients computed during expansion may still reference the
             * internal inverse variable u = $SeriesInvVar$ -- any part of
             * the input that ended up treated as "constant in u" keeps it
             * in closed form (Log[x] expands to Log[1/u] -> -Log[u], which
             * is constant in u and so lives in the coefficient). Substitute
             * u -> 1/x in every coefficient so $SeriesInvVar$ never leaks
             * back to the user. simp() normalises 1/(1/x) back to x. */
            Expr* u_to_inv_x = mk_power(expr_copy(x), expr_new_integer(-1));
            for (size_t ci = 0; ci < s->coef_count; ci++) {
                Expr* new_c = replace_all_of(s->coefs[ci], u_sym, u_to_inv_x);
                Expr* simp_c = simp(new_c);
                expr_free(s->coefs[ci]);
                s->coefs[ci] = simp_c;
            }
            expr_free(u_to_inv_x);
        }
        result = so_to_expr(s);
        so_free(s);
    }

    /* MMA-faithful branch-point wrapper. A top-level branch-point handler
     * populated ctx.pending_* with (add_const, log_coef, discriminator).
     * Build
     *     Plus[ add_const,
     *           Times[ log_coef, Log[x - x0] ],   (omitted if log_coef NULL)
     *           Times[ discriminator, result ] ]
     * (omitting trivial 0/1 factors). For at_infinity the substituted variable
     * makes the metadata's x/x0 stale; just discard the wrapper there. */
    if (result && !at_infinity && ctx.pending_discriminator) {
        Expr** parts = calloc(3, sizeof(Expr*));
        size_t pc = 0;
        if (ctx.pending_add_const && !is_lit_zero(ctx.pending_add_const)) {
            parts[pc++] = ctx.pending_add_const;
            ctx.pending_add_const = NULL;
        }
        if (ctx.pending_log_coef) {
            parts[pc++] = build_log_term(ctx.pending_log_coef, x, x0_eval);
        }
        Expr* wrapped_series = simp(mk_times(ctx.pending_discriminator, result));
        ctx.pending_discriminator = NULL;
        parts[pc++] = wrapped_series;
        Expr* sum;
        if (pc == 1) { sum = parts[0]; free(parts); }
        else         { sum = expr_new_function(mk_symbol("Plus"), parts, pc); free(parts); }
        result = simp(sum);
    }
    /* Clean up any remaining pending metadata (e.g. if the wrap was skipped). */
    if (ctx.pending_add_const)     { expr_free(ctx.pending_add_const);     ctx.pending_add_const = NULL; }
    if (ctx.pending_log_coef)      { expr_free(ctx.pending_log_coef);      ctx.pending_log_coef = NULL; }
    if (ctx.pending_discriminator) { expr_free(ctx.pending_discriminator); ctx.pending_discriminator = NULL; }

    if (at_infinity) { expr_free(f_sub_owned); expr_free(u_sym); expr_free(x0_use); }
    else             expr_free(x0_eval);
    if (f_body_owned) expr_free(f_body_owned);
    expr_free(f_eval);

    /* Wrap the result in the factored-out prefactor, if any. We keep the
     * factorisation at the *expression* level (Times[x^a, SeriesData[...]])
     * so the SeriesData printer still formats the series and the outer
     * Times adds the symbolic x^a factor. */
    if (prefactor) {
        if (!result) { expr_free(prefactor); return NULL; }
        Expr* wrapped = mk_times(prefactor, result);
        return wrapped;
    }
    return result;
}

Expr* builtin_series(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;

    /* Series uses PolynomialQuotient/Remainder, GCD and Together while
     * extracting the leading-term expansion — coefficients must be
     * rational. Convert inexact inputs, run, and numericalise the
     * resulting series back. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_series);
    }

    Expr* f = res->data.function.args[0];

    /* List threading on the first argument (not ATTR_LISTABLE since the
     * specs must not thread). */
    if (has_symbol_head(f, "List")) {
        size_t n = f->data.function.arg_count;
        Expr** threaded = calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr** new_args = calloc(res->data.function.arg_count, sizeof(Expr*));
            new_args[0] = expr_copy(f->data.function.args[i]);
            for (size_t j = 1; j < res->data.function.arg_count; j++) {
                new_args[j] = expr_copy(res->data.function.args[j]);
            }
            Expr* call = expr_new_function(mk_symbol("Series"), new_args,
                                           res->data.function.arg_count);
            threaded[i] = eval_and_free(call);
        }
        Expr* lst = expr_new_function(mk_symbol("List"), threaded, n);
        free(threaded);
        return lst;
    }

    /* Separate the Assumptions option (if present) from the expansion specs.
     * Series[f, spec..., Assumptions -> assm]: the option may appear in any
     * trailing position. We pull it out so it is not mistaken for an extra
     * expansion spec (a leading-order spec `x -> x0` is also a Rule, so we
     * key on the LHS symbol being exactly Assumptions). `assm` and the
     * entries of `specs` borrow pointers into `res`. */
    size_t raw_count = res->data.function.arg_count - 1;
    Expr** raw = res->data.function.args + 1;
    Expr* assm = NULL;
    Expr** specs = calloc(raw_count, sizeof(Expr*));
    size_t spec_count = 0;
    for (size_t i = 0; i < raw_count; i++) {
        Expr* s = raw[i];
        if (has_symbol_head(s, "Rule") && s->data.function.arg_count == 2 &&
            s->data.function.args[0]->type == EXPR_SYMBOL &&
            strcmp(s->data.function.args[0]->data.symbol.name, "Assumptions") == 0) {
            assm = s->data.function.args[1];
            continue;
        }
        specs[spec_count++] = s;
    }
    if (spec_count == 0) { free(specs); return NULL; }

    /* Multivariate: process right-to-left. Expand f first in the rightmost
     * variable, which produces a series whose coefficients are themselves
     * expressions in the remaining variables. Then expand each coefficient
     * in the next variable. */
    if (spec_count == 1) {
        Expr* x_arg; Expr* x0_arg; int64_t n_val; bool leading;
        if (!parse_series_spec(specs[0], &x_arg, &x0_arg, &n_val, &leading)) {
            free(specs); return NULL;
        }
        int x_sign = assm ? assumption_sign_of(assm, x_arg) : 0;
        Expr* out = do_series_single(f, x_arg, x0_arg, n_val, leading, x_sign);
        free(specs);
        return out;
    }

    /* For multivariate: expand outermost (leftmost) first -- the natural
     * reading matches Mathematica: Series[f, spec_x, spec_y] means
     * expand in x first, getting coefs that are expressions in y, then
     * expand each coef in y. */
    Expr* x_arg; Expr* x0_arg; int64_t n_val; bool leading;
    if (!parse_series_spec(specs[0], &x_arg, &x0_arg, &n_val, &leading)) {
        free(specs); return NULL;
    }
    int x_sign = assm ? assumption_sign_of(assm, x_arg) : 0;
    Expr* outer = do_series_single(f, x_arg, x0_arg, n_val, leading, x_sign);
    if (!outer) { free(specs); return NULL; }

    /* outer is SeriesData[x, x0, {a0, a1, ...}, 0, n+1, 1]. Recurse on
     * each coefficient with the remaining specs (forwarding the Assumptions
     * option so the inner variables see it too). */
    if (!has_symbol_head(outer, "SeriesData") || outer->data.function.arg_count != 6) {
        free(specs);
        return outer;
    }
    Expr* coefs = outer->data.function.args[2];
    for (size_t i = 0; i < coefs->data.function.arg_count; i++) {
        /* Build Series[ai, spec_y, spec_z, ..., Assumptions -> assm]. */
        size_t rest = spec_count - 1;
        size_t extra = assm ? 1 : 0;
        Expr** new_args = calloc(rest + 1 + extra, sizeof(Expr*));
        new_args[0] = expr_copy(coefs->data.function.args[i]);
        for (size_t j = 0; j < rest; j++) new_args[j + 1] = expr_copy(specs[j + 1]);
        if (assm) {
            new_args[rest + 1] = mk_fn2("Rule", mk_symbol("Assumptions"),
                                        expr_copy(assm));
        }
        Expr* call = expr_new_function(mk_symbol("Series"), new_args, rest + 1 + extra);
        free(new_args);
        Expr* expanded = eval_and_free(call);
        expr_free(coefs->data.function.args[i]);
        coefs->data.function.args[i] = expanded;
    }
    free(specs);
    return outer;
}

/* SeriesCoefficient[f, {x, x0, k}] -- the coefficient of (x - x0)^k in the
 * power-series expansion of f about x = x0. Computed by expanding with Series
 * and extracting the k-th coefficient from the resulting SeriesData; this is
 * general (works for any head Series can expand) for a concrete integer index
 * and a finite expansion point. Composite results (a prefactor times a
 * SeriesData, e.g. the asymptotic expansions at Infinity) and non-integer
 * indices are left unevaluated. The symbolic-index general term (a Piecewise)
 * is not produced. */
Expr* builtin_seriescoefficient(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    /* General term SeriesCoefficient[ProductLog[x], {x, 0, n}] with symbolic n:
     *   Piecewise[{{(-n)^(n-1)/n!, n >= 1}}, 0]. */
    {
        Expr* f0  = res->data.function.args[0];
        Expr* sp0 = res->data.function.args[1];
        if (has_symbol_head(f0, "ProductLog") && f0->data.function.arg_count == 1 &&
            sp0->type == EXPR_FUNCTION && has_symbol_head(sp0, "List") &&
            sp0->data.function.arg_count == 3) {
            Expr* xv  = sp0->data.function.args[0];
            Expr* x0v = sp0->data.function.args[1];
            Expr* nv  = sp0->data.function.args[2];
            if (expr_eq(f0->data.function.args[0], xv) &&
                is_lit_zero(x0v) && nv->type == EXPR_SYMBOL) {
                Expr* val = mk_times(
                    mk_power(mk_times(expr_new_integer(-1), expr_copy(nv)),
                             mk_plus(expr_new_integer(-1), expr_copy(nv))),
                    mk_power(mk_fn1("Factorial", expr_copy(nv)),
                             expr_new_integer(-1)));
                Expr* cond = mk_fn2("GreaterEqual", expr_copy(nv), expr_new_integer(1));
                Expr* pair = expr_new_function(mk_symbol("List"),
                                 (Expr*[]){ val, cond }, 2);
                Expr* pairs = expr_new_function(mk_symbol("List"),
                                  (Expr*[]){ pair }, 1);
                Expr* pw = expr_new_function(mk_symbol("Piecewise"),
                               (Expr*[]){ pairs, expr_new_integer(0) }, 2);
                return eval_and_free(pw);
            }
        }
    }

    /* General term SeriesCoefficient[FresnelC/FresnelS[x], {x, 0, n}], symbolic n:
     *   FresnelC: Piecewise[{{ (-1)^((n-1)/4) 2^(1-n) Pi^(n/2) /
     *                          (n Gamma[(1+n)/4] Gamma[(3+n)/4]), Mod[n-1,4]==0 && n>=1 }}, 0]
     *   FresnelS: same value with the (-1) exponent (n-3)/4 and Mod[n-3,4]==0 && n>=3. */
    {
        Expr* f0  = res->data.function.args[0];
        Expr* sp0 = res->data.function.args[1];
        int isC = has_symbol_head(f0, "FresnelC");
        int isS = has_symbol_head(f0, "FresnelS");
        if ((isC || isS) && f0->data.function.arg_count == 1 &&
            sp0->type == EXPR_FUNCTION && has_symbol_head(sp0, "List") &&
            sp0->data.function.arg_count == 3) {
            Expr* xv  = sp0->data.function.args[0];
            Expr* x0v = sp0->data.function.args[1];
            Expr* nv  = sp0->data.function.args[2];
            if (expr_eq(f0->data.function.args[0], xv) &&
                is_lit_zero(x0v) && nv->type == EXPR_SYMBOL) {
                int64_t off = isC ? 1 : 3;                 /* power class 4m+off */
                Expr* nmoff = mk_plus(expr_copy(nv), expr_new_integer(-off));  /* n-off */
                /* (-1)^((n-off)/4) */
                Expr* sign = mk_power(expr_new_integer(-1),
                                 mk_times(make_rational(1, 4), expr_copy(nmoff)));
                /* 2^(1-n) */
                Expr* pow2 = mk_power(expr_new_integer(2),
                                 mk_plus(expr_new_integer(1),
                                         mk_times(expr_new_integer(-1), expr_copy(nv))));
                /* Pi^(n/2) */
                Expr* pipow = mk_power(mk_symbol("Pi"),
                                 mk_times(make_rational(1, 2), expr_copy(nv)));
                /* n Gamma[(1+n)/4] Gamma[(3+n)/4] */
                Expr* gam1 = mk_fn1("Gamma", mk_times(make_rational(1, 4),
                                 mk_plus(expr_new_integer(1), expr_copy(nv))));
                Expr* gam2 = mk_fn1("Gamma", mk_times(make_rational(1, 4),
                                 mk_plus(expr_new_integer(3), expr_copy(nv))));
                Expr* denom = mk_times(expr_copy(nv), mk_times(gam1, gam2));
                Expr* val = mk_times(mk_times(sign, pow2),
                                 mk_times(pipow, mk_power(denom, expr_new_integer(-1))));
                /* Mod[n-off, 4] == 0 && n >= off */
                Expr* cond = mk_fn2("And",
                    mk_fn2("Equal",
                        mk_fn2("Mod", expr_copy(nmoff), expr_new_integer(4)),
                        expr_new_integer(0)),
                    mk_fn2("GreaterEqual", expr_copy(nv), expr_new_integer(off)));
                expr_free(nmoff);
                Expr* pair = expr_new_function(mk_symbol("List"),
                                 (Expr*[]){ val, cond }, 2);
                Expr* pairs = expr_new_function(mk_symbol("List"),
                                  (Expr*[]){ pair }, 1);
                Expr* pw = expr_new_function(mk_symbol("Piecewise"),
                               (Expr*[]){ pairs, expr_new_integer(0) }, 2);
                return eval_and_free(pw);
            }
        }
    }

    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_seriescoefficient);
    }

    Expr* f    = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    Expr* x_arg; Expr* x0_arg; int64_t k; bool leading;
    if (!parse_series_spec(spec, &x_arg, &x0_arg, &k, &leading)) return NULL;
    if (leading) return NULL;  /* {x, x0} 2-arg form has no coefficient index */

    /* Expand to an order that includes (x - x0)^k. */
    int64_t order = k > 0 ? k : 0;
    Expr** sp = calloc(3, sizeof(Expr*));
    sp[0] = expr_copy(x_arg);
    sp[1] = expr_copy(x0_arg);
    sp[2] = expr_new_integer(order);
    Expr* spec2 = expr_new_function(mk_symbol("List"), sp, 3);
    free(sp);
    Expr* call = expr_new_function(mk_symbol("Series"),
                     (Expr*[]){ expr_copy(f), spec2 }, 2);
    Expr* s = eval_and_free(call);
    if (!s) return NULL;

    Expr* result = NULL;
    if (has_symbol_head(s, "SeriesData") && s->data.function.arg_count == 6 &&
        expr_eq(s->data.function.args[0], x_arg)) {
        Expr* coefs = s->data.function.args[2];
        Expr* nminE = s->data.function.args[3];
        Expr* denE  = s->data.function.args[5];
        if (coefs->type == EXPR_FUNCTION &&
            nminE->type == EXPR_INTEGER && denE->type == EXPR_INTEGER) {
            int64_t nmin = nminE->data.integer, den = denE->data.integer;
            long long j = (long long)k * den - nmin;   /* power k = (nmin+j)/den */
            size_t len = coefs->data.function.arg_count;
            if (j >= 0 && (size_t)j < len)
                result = expr_copy(coefs->data.function.args[(size_t)j]);
            else
                result = expr_new_integer(0);
        }
    } else if (expr_free_of(s, x_arg)) {
        /* Series collapsed to something free of x (e.g. a constant). */
        result = (k == 0) ? expr_copy(s) : expr_new_integer(0);
    }
    /* Otherwise (composite prefactor * SeriesData, expansion in 1/x, ...)
     * leave unevaluated. */

    expr_free(s);
    return result;
}

/* ----------------------------------------------------------------------------
 * Initialization
 * -------------------------------------------------------------------------- */

void series_init(void) {
    symtab_get_def("SeriesData")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Series", builtin_series);
    symtab_get_def("Series")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_add_builtin("Normal", builtin_normal);
    symtab_get_def("Normal")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("SeriesCoefficient", builtin_seriescoefficient);
    symtab_get_def("SeriesCoefficient")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
