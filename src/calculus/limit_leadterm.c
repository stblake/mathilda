/*
 * limit_leadterm.c -- limits at a finite point from leading terms.
 *
 * For t -> 0, every expression e built here has a leading term
 *
 *     e = c t^p (1 + o(1)),        c != 0 free of t, p rational,
 *
 * and the limit is read off it: c when p = 0, 0 when p > 0.  The leading
 * term is computed STRUCTURALLY wherever that is exact, and by a truncated
 * Series only where it must be:
 *
 *   constant      c t^0                               (c provably nonzero)
 *   t             1 t^1
 *   a b           lead(a) lead(b)                     (product of nonzero c's)
 *   b^n           lead(b)^n                           (n an integer)
 *   F(G(t))       a b^s t^(s q)  when G = b t^q + ..., q > 0, and
 *                                F(y) = a y^s + ... at y = 0 with s an integer
 *   otherwise     the first provably nonzero term of Series[e, {t, 0, k}],
 *                 k = 4, 8, 16, ... (exact rational arithmetic for the
 *                 elementary functions at 0, so a cancelling sum is resolved
 *                 by expanding just far enough)
 *
 * The composition rule is what makes a deeply cancelling composite cheap: in
 * f(g(t))/g(f(t)) with f ~ a t^15 and g ~ b t^11 the quotient's leading
 * terms sit at t^165, which a direct expansion reaches only by composing
 * series truncated at order 166, while f(g(t)) ~ a (b t^11)^15 needs f to
 * order 15 and g to order 11.  It is sound exactly when its two side
 * conditions hold, and both are checked: the inner G must vanish at the
 * point (q > 0, so F is evaluated near y = 0, where its expansion holds), and
 * F's leading coefficient a must be nonzero (so a b^s t^(s q) is the leading
 * term and not cancelled).  s must be an integer so b^s has no branch.
 *
 * Every coefficient a leading term rests on is tested for zero exactly
 * (literal rationals) or by the zero tester; an undecided coefficient
 * declines the whole computation rather than guessing.
 */
#include "limit_leadterm.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "zero_test.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define LT_SERIES_START   4
#define LT_SERIES_MAX     64
#define LT_MAX_DEPTH      8

/* Leading term c t^p (both owned), or `zero`: e vanishes identically near
 * the point (c, p NULL). */
typedef struct { Expr* c; Expr* p; bool zero; } LeadTerm;

static void lt_clear(LeadTerm* L) {
    if (L->c) expr_free(L->c);
    if (L->p) expr_free(L->p);
    L->c = L->p = NULL;
    L->zero = false;
}

static Expr* fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* ev(Expr* e) { return eval_and_free(e); }

static bool contains(const Expr* e, const Expr* t) {
    if (expr_eq(e, t)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains(e->data.function.head, t)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains(e->data.function.args[i], t)) return true;
    return false;
}

static bool is_rational_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Rational;
}

static bool is_head(const Expr* e, const char* sym) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* 1: provably nonzero; 0: provably zero; -1: undecided. */
static int nonzero(const Expr* c) {
    if (is_rational_number(c))
        return !(c->type == EXPR_INTEGER && c->data.integer == 0);
    ZeroTestResult z = zero_test_decide(c);
    if (z == ZERO_TEST_FALSE) return 1;
    if (z == ZERO_TEST_TRUE) return 0;
    return -1;
}

/* Sign of a rational exponent: -1, 0, +1. */
static int rational_sign(const Expr* p) {
    const Expr* n = is_head(p, SYM_Rational) ? p->data.function.args[0] : p;
    if (n->type == EXPR_INTEGER) return n->data.integer > 0 ? 1 : n->data.integer < 0 ? -1 : 0;
    if (n->type == EXPR_BIGINT) return mpz_sgn(n->data.bigint);
    return 0;
}

static bool lead(Expr* e, Expr* t, LeadTerm* out, int depth);

/* The first provably nonzero term of Series[e, {t, 0, k}], raising k. */
static bool lead_series(Expr* e, Expr* t, LeadTerm* out) {
    for (int64_t k = LT_SERIES_START; k <= LT_SERIES_MAX; k *= 2) {
        Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
                         (Expr*[]){ expr_copy(t), expr_new_integer(0), expr_new_integer(k) }, 3);
        Expr* s = ev(fn2("Series", expr_copy(e), spec));
        if (!s) return false;
        if (!contains(s, t)) {                       /* a constant */
            int nz = nonzero(s);
            if (nz == 1) { out->c = s; out->p = expr_new_integer(0); return true; }
            expr_free(s);
            return false;                            /* zero or undecided */
        }
        if (!is_head(s, SYM_SeriesData) || s->data.function.arg_count != 6) {
            expr_free(s);
            return false;
        }
        Expr** a = s->data.function.args;
        if (!expr_eq(a[0], t) || !(a[1]->type == EXPR_INTEGER && a[1]->data.integer == 0) ||
            !is_head(a[2], SYM_List) || a[3]->type != EXPR_INTEGER ||
            a[5]->type != EXPR_INTEGER || a[5]->data.integer < 1) {
            expr_free(s);
            return false;
        }
        int64_t nmin = a[3]->data.integer, den = a[5]->data.integer;
        Expr* cl = a[2];
        for (size_t i = 0; i < cl->data.function.arg_count; i++) {
            Expr* c = cl->data.function.args[i];
            if (contains(c, t)) { expr_free(s); return false; }   /* Log t, branch phase */
            int nz = nonzero(c);
            if (nz == 0) continue;
            if (nz < 0) { expr_free(s); return false; }
            out->c = expr_copy(c);
            out->p = ev(fn2("Times", expr_new_integer(nmin + (int64_t)i),
                            fn2("Power", expr_new_integer(den), expr_new_integer(-1))));
            expr_free(s);
            return true;
        }
        expr_free(s);                                /* all zero so far: go deeper */
    }
    /* No term through order LT_SERIES_MAX: e is identically zero when the
     * zero tester agrees (Sin[y]^2 + Cos[y]^2 - 1). */
    if (zero_test_decide(e) == ZERO_TEST_TRUE) { out->zero = true; return true; }
    return false;
}

/* Does every occurrence of t in e lie inside a copy of g? */
static bool only_through(const Expr* e, const Expr* g, const Expr* t) {
    if (expr_eq(e, g)) return true;
    if (expr_eq(e, t)) return false;
    if (e->type != EXPR_FUNCTION) return true;
    if (!only_through(e->data.function.head, g, t)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!only_through(e->data.function.args[i], g, t)) return false;
    return true;
}

/* The largest proper subexpression G != t of `root` through which every
 * occurrence of t passes (pre-order: the first one found is maximal).
 * Borrowed from root; NULL when only t itself qualifies. */
static Expr* find_inner(Expr* root, Expr* e, Expr* t) {
    if (e->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* a = e->data.function.args[i];
        if (expr_eq(a, t) || !contains(a, t)) continue;
        if (a->type == EXPR_FUNCTION && only_through(root, a, t)) return a;
        Expr* r = find_inner(root, a, t);
        if (r) return r;
    }
    return NULL;
}

static Expr* replace_all(Expr* e, Expr* from, Expr* to) {
    if (expr_eq(e, from)) return expr_copy(to);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = replace_all(e->data.function.args[i], from, to);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}

/* F(G(t)) with G -> 0: lead(F, y) composed with lead(G, t). */
static bool lead_compose(Expr* e, Expr* t, LeadTerm* out, int depth) {
    Expr* G = find_inner(e, e, t);
    if (!G) return false;
    LeadTerm g = { NULL, NULL, false };
    if (!lead(G, t, &g, depth + 1)) return false;
    if (g.zero || rational_sign(g.p) <= 0) { lt_clear(&g); return false; }   /* G must vanish */
    char name[48];
    snprintf(name, sizeof name, "LeadTerm`$y%d", depth);
    Expr* y = expr_new_symbol(name);
    Expr* F = ev(replace_all(e, G, y));
    LeadTerm f = { NULL, NULL, false };
    bool ok = F && !contains(F, t) && lead(F, y, &f, depth + 1) &&
              (f.zero || f.p->type == EXPR_INTEGER);            /* no branch in b^s */
    if (F) expr_free(F);
    expr_free(y);
    if (!ok) { lt_clear(&f); lt_clear(&g); return false; }
    if (f.zero) { lt_clear(&f); lt_clear(&g); out->zero = true; return true; }
    out->c = ev(fn2("Times", expr_copy(f.c), fn2("Power", expr_copy(g.c), expr_copy(f.p))));
    out->p = ev(fn2("Times", expr_copy(f.p), expr_copy(g.p)));
    lt_clear(&f); lt_clear(&g);
    if (!out->c || !out->p || nonzero(out->c) != 1) { lt_clear(out); return false; }
    return true;
}

static bool lead(Expr* e, Expr* t, LeadTerm* out, int depth) {
    out->c = out->p = NULL;
    out->zero = false;
    if (depth > LT_MAX_DEPTH) return false;
    if (!contains(e, t)) {
        if (nonzero(e) != 1) return false;
        out->c = expr_copy(e); out->p = expr_new_integer(0);
        return true;
    }
    if (expr_eq(e, t)) { out->c = expr_new_integer(1); out->p = expr_new_integer(1); return true; }

    if (is_head(e, SYM_Times)) {
        Expr* c = expr_new_integer(1);
        Expr* p = expr_new_integer(0);
        bool ok = true;
        bool zero = false;
        for (size_t i = 0; ok && i < e->data.function.arg_count; i++) {
            LeadTerm f;
            if (!lead(e->data.function.args[i], t, &f, depth + 1)) { ok = false; break; }
            if (f.zero) { zero = true; continue; }
            c = ev(fn2("Times", c, f.c));
            p = ev(fn2("Plus", p, f.p));
        }
        if (ok && zero) { expr_free(c); expr_free(p); out->zero = true; return true; }
        if (ok) { out->c = c; out->p = p; return true; }
        expr_free(c); expr_free(p);
    } else if (is_head(e, SYM_Power) && e->data.function.arg_count == 2 &&
               e->data.function.args[1]->type == EXPR_INTEGER) {
        LeadTerm b;
        Expr* n = e->data.function.args[1];
        if (lead(e->data.function.args[0], t, &b, depth + 1) && b.zero) {
            if (n->data.integer > 0) { out->zero = true; return true; }
            return false;                                   /* 1/0 */
        }
        if (b.c) {
            out->c = ev(fn2("Power", b.c, expr_copy(n)));
            out->p = ev(fn2("Times", b.p, expr_copy(n)));
            return true;
        }
    } else if (lead_compose(e, t, out, depth)) {
        return true;
    }
    return lead_series(e, t, out);
}

Expr* limit_leading_term_value(Expr* f, Expr* x, Expr* x0) {
    if (!f || !x || !x0 || x->type != EXPR_SYMBOL || contains(x0, x)) return NULL;
    Expr* t = expr_new_symbol("LeadTerm`$t");
    /* Shift to t -> 0: f(x0 + t). */
    Expr* shift = (x0->type == EXPR_INTEGER && x0->data.integer == 0)
                      ? expr_copy(t) : fn2("Plus", expr_copy(x0), expr_copy(t));
    Expr* g = ev(replace_all(f, x, shift));
    expr_free(shift);
    Expr* r = NULL;
    LeadTerm L;
    if (g && lead(g, t, &L, 0)) {
        int s = L.zero ? 1 : rational_sign(L.p);
        if (s == 0 && L.c && !contains(L.c, t)) { r = L.c; L.c = NULL; }
        else if (s > 0) r = expr_new_integer(0);
        lt_clear(&L);
    }
    if (g) expr_free(g);
    expr_free(t);
    return r;
}
