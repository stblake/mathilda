/* expsplit.c -- powers with a symbolic exponent as polynomial generators.
 *
 * The reference reads b^(n + k t) as the monomial b^n (b^t)^k: b^t is a
 * generator, the integer coefficient k its degree, the numeric part n a power
 * of b itself.  So
 *     Variables[x^(1 + e) + x^(2 + e)]   is {x, x^e},
 *     Factor[x^(2 + e) - x^e]            is (-1 + x) x^e (1 + x),
 *     FactorList[(x + y)^(-3 + 3 e)]     is {{1, 1}, {(x + y)^e, 3}, {x + y, -3}}.
 * A rational non-integer coefficient keeps its term whole (x^(e/2) is its own
 * generator, as there).  Mathilda's polynomial code takes every power with a
 * symbolic exponent as one opaque atom, so these builtins abstract first:
 * expsplit_abstract rewrites each such power as b^n K^k, K a private symbol
 * standing for b^t, and the result is mapped back afterwards. */

#include "expsplit.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void expsplit_free(ExpSplitTab* t) {
    for (size_t i = 0; i < t->n; i++) {
        expr_free(t->base[i]); expr_free(t->term[i]); expr_free(t->sym[i]);
    }
    free(t->base); free(t->term); free(t->sym);
    memset(t, 0, sizeof *t);
}

static bool is_head(const Expr* e, const char* h) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* The private symbol standing for base^term, created on first use. */
static Expr* gen_symbol(ExpSplitTab* t, const Expr* base, const Expr* term) {
    for (size_t i = 0; i < t->n; i++)
        if (expr_eq(t->base[i], (Expr*)base) && expr_eq(t->term[i], (Expr*)term))
            return expr_copy(t->sym[i]);
    if (t->n == t->cap) {
        size_t cap = t->cap ? 2 * t->cap : 8;
        t->base = realloc(t->base, cap * sizeof(Expr*));
        t->term = realloc(t->term, cap * sizeof(Expr*));
        t->sym  = realloc(t->sym,  cap * sizeof(Expr*));
        t->cap = cap;
    }
    char name[48];
    snprintf(name, sizeof name, "ExpSplit`g%zu", t->n);
    t->base[t->n] = expr_copy((Expr*)base);
    t->term[t->n] = expr_copy((Expr*)term);
    t->sym[t->n]  = expr_new_symbol(name);
    t->n++;
    return expr_copy(t->sym[t->n - 1]);
}

/* Split a symbolic term c m into its numeric coefficient and monomial. */
static void term_coeff(const Expr* term, Expr** c, Expr** m) {
    if (is_head(term, SYM_Times) && term->data.function.arg_count >= 2 &&
        expr_is_numeric_like(term->data.function.args[0])) {
        *c = expr_copy(term->data.function.args[0]);
        size_t k = term->data.function.arg_count - 1;
        if (k == 1) *m = expr_copy(term->data.function.args[1]);
        else {
            Expr** a = malloc(k * sizeof(Expr*));
            for (size_t i = 0; i < k; i++) a[i] = expr_copy(term->data.function.args[i + 1]);
            *m = expr_new_function(expr_new_symbol(SYM_Times), a, k);
            free(a);
        }
        return;
    }
    *c = expr_new_integer(1);
    *m = expr_copy((Expr*)term);
}

/* b^p rewritten over generators, or NULL when p has nothing to split.  With
 * `all`, a power that is already a generator (x^e) maps to its symbol too, so
 * it meets the split ones: x^(2 + e) - x^e is x^2 K - K. */
static Expr* split_power(const Expr* b, const Expr* p, ExpSplitTab* t, bool all) {
    size_t nt = is_head(p, SYM_Plus) ? p->data.function.arg_count : 1;
    Expr* const* terms = is_head(p, SYM_Plus) ? p->data.function.args : (Expr* const*)&p;
    size_t nsym = 0;
    bool has_num = false, nontrivial = false;
    for (size_t i = 0; i < nt; i++) {
        if (expr_is_numeric_like(terms[i])) { has_num = true; continue; }
        nsym++;
        Expr *c, *m;
        term_coeff(terms[i], &c, &m);
        if (!(c->type == EXPR_INTEGER && c->data.integer == 1) && expr_is_integer_like(c))
            nontrivial = true;                       /* x^(2 e), x^(-e) */
        expr_free(c); expr_free(m);
    }
    if (nsym == 0 || (!all && !has_num && nsym == 1 && !nontrivial)) return NULL;
    Expr** f = malloc((nt + 1) * sizeof(Expr*));
    size_t nf = 0;
    Expr* num = expr_new_integer(0);
    for (size_t i = 0; i < nt; i++) {
        if (expr_is_numeric_like(terms[i])) {
            num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                    (Expr*[]){ num, expr_copy(terms[i]) }, 2));
            continue;
        }
        Expr *c, *m;
        term_coeff(terms[i], &c, &m);
        if (expr_is_integer_like(c)) {
            Expr* g = gen_symbol(t, b, m);
            f[nf++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ g, c }, 2);
            expr_free(m);
        } else {                                     /* x^(e/2): its own generator */
            Expr* whole = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ c, m }, 2);
            f[nf++] = gen_symbol(t, b, whole);
            expr_free(whole);
        }
    }
    f[nf++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy((Expr*)b), num }, 2);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), f, nf);
    free(f);
    return r;
}

/* numeric_only: split only powers whose base is a number (b^t with b in Q),
 * leaving x^t whole for the caller's own power rules. */
static Expr* abstract_rec(const Expr* e, ExpSplitTab* t, bool numeric_only) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = abstract_rec(e->data.function.args[i], t, numeric_only);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (is_head(r, SYM_Power) && r->data.function.arg_count == 2 &&
        (!numeric_only || expr_is_numeric_like(r->data.function.args[0]))) {
        const Expr* p = r->data.function.args[1];
        if (!expr_is_numeric_like(p)) {
            Expr* s = split_power(r->data.function.args[0], p, t, true);
            if (s) { expr_free(r); return s; }
        }
    }
    return r;
}

static bool has_splittable(const Expr* e, bool numeric_only) {
    if (e->type != EXPR_FUNCTION) return false;
    if (is_head(e, SYM_Power) && e->data.function.arg_count == 2 &&
        !expr_is_numeric_like(e->data.function.args[1]) &&
        (!numeric_only || expr_is_numeric_like(e->data.function.args[0]))) {
        ExpSplitTab t = {0};
        Expr* s = split_power(e->data.function.args[0], e->data.function.args[1], &t, false);
        bool yes = s != NULL;
        expr_free(s);
        expsplit_free(&t);
        if (yes) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_splittable(e->data.function.args[i], numeric_only)) return true;
    return false;
}

static Expr* abstract_top(const Expr* e, ExpSplitTab* t, bool numeric_only) {
    static int off = -1;                 /* MATHILDA_NO_EXPSPLIT: A/B kill switch */
    if (off < 0) off = getenv("MATHILDA_NO_EXPSPLIT") != NULL;
    if (off) return NULL;
    if (!has_splittable(e, numeric_only)) return NULL;
    return eval_and_free(abstract_rec(e, t, numeric_only));
}

Expr* expsplit_abstract(const Expr* e, ExpSplitTab* t) {
    return abstract_top(e, t, false);
}

Expr* expsplit_abstract_numeric(const Expr* e, ExpSplitTab* t) {
    return abstract_top(e, t, true);
}

/* Every power with a symbolic exponent becomes generators, a lone (c + d)^n
 * included (expsplit_abstract leaves one that has nothing to split); NULL
 * when e has no such power. */
static bool has_symbolic_power(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (is_head(e, SYM_Power) && e->data.function.arg_count == 2 &&
        !expr_is_numeric_like(e->data.function.args[1]))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_symbolic_power(e->data.function.args[i])) return true;
    return false;
}

Expr* expsplit_abstract_all(const Expr* e, ExpSplitTab* t) {
    if (getenv("MATHILDA_NO_EXPSPLIT") || !has_symbolic_power(e)) return NULL;
    return eval_and_free(abstract_rec(e, t, false));
}

int expsplit_index(const ExpSplitTab* t, const Expr* e) {
    if (e->type != EXPR_SYMBOL) return -1;
    for (size_t i = 0; i < t->n; i++)
        if (t->sym[i]->data.symbol.name == e->data.symbol.name) return (int)i;
    return -1;
}

static Expr* restore_rec(const Expr* e, const ExpSplitTab* t);

Expr* expsplit_generator(const ExpSplitTab* t, size_t i) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ restore_rec(t->base[i], t), restore_rec(t->term[i], t) }, 2));
}

/* A fresh copy of e with K^k -> b^(k t) and K -> b^t.  Built anew rather
 * than rewritten in place: expr_copy shares nodes, and a shared node already
 * marked evaluated would not be evaluated again (x^e x would stay apart). */
static Expr* restore_rec(const Expr* e, const ExpSplitTab* t) {
    /* A generator's base and term can hold earlier generators (E^(E^x) is
     * K1 = E^K0 once the inner E^x is K0), so they are restored too. */
    if (e->type == EXPR_SYMBOL) {
        int i = expsplit_index(t, e);
        if (i < 0) return expr_copy((Expr*)e);
        return expr_new_function(expr_new_symbol(SYM_Power),
                   (Expr*[]){ restore_rec(t->base[i], t), restore_rec(t->term[i], t) }, 2);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (is_head(e, SYM_Power) && e->data.function.arg_count == 2) {
        int i = expsplit_index(t, e->data.function.args[0]);
        if (i >= 0) {
            Expr* k = restore_rec(e->data.function.args[1], t);
            Expr* ex = expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ k, restore_rec(t->term[i], t) }, 2);
            return expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ restore_rec(t->base[i], t), ex }, 2);
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t j = 0; j < n; j++) a[j] = restore_rec(e->data.function.args[j], t);
    Expr* r = expr_new_function(restore_rec(e->data.function.head, t), a, n);
    free(a);
    return r;
}

Expr* expsplit_restore(Expr* e, const ExpSplitTab* t) {
    if (!e) return NULL;
    Expr* r = restore_rec(e, t);
    expr_free(e);
    return eval_and_free(r);
}
