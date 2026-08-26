/*
 * reduce_companions.c
 *
 * Companion builtins for `Reduce` (REDUCE_PLAN.md, Phase 8).  v1: LogicalExpand
 * + a minimal NotElement head.
 *
 * LogicalExpand distributes a logical statement to disjunctive normal form (an
 * Or of Ands of literals), applying idempotence / complementation / absorption
 * contractions, and collapsing to True (tautology) or False (contradiction)
 * when the statement decides.  Every non-connective subexpression is treated as
 * an OPAQUE Boolean atom -- a symbol, a relation `x == a`, a membership
 * `Element[..]` -- with NO domain reasoning, exactly as Mathematica's
 * LogicalExpand does.  Two relational atoms are complementary iff one is the
 * (head-flipped) logical negation of the other (`x==a` / `x!=a`, `x<1` / `x>=1`,
 * `Element` / `NotElement`, `a` / `!a`).
 *
 * The True/False collapse is sound *and* complete without truth-table
 * enumeration: over independent opaque atoms a DNF is unsatisfiable iff every
 * clause holds a complementary pair -- i.e. it distributes to ZERO surviving
 * clauses.  So `phi` empty => False, and `Not[phi]` empty => True (the negation
 * is unsatisfiable, hence `phi` is a tautology).
 */
#include "reduce_companions.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#include "expr.h"
#include "eval.h"
#include "attr.h"
#include "symtab.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ *
 *  Small Expr helpers                                                 *
 * ------------------------------------------------------------------ */

static bool is_head_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == s;
}
static size_t nargs(const Expr* e) { return e->data.function.arg_count; }
static Expr*  argn(const Expr* e, size_t i) { return e->data.function.args[i]; }

/* Build a fresh 2-arg relation `newhead[d.arg0, d.arg1]`. */
static Expr* relhead_swap(const Expr* d, const char* newhead) {
    Expr* a[2] = { expr_copy(argn(d, 0)), expr_copy(argn(d, 1)) };
    return expr_new_function(expr_new_symbol(newhead), a, 2);
}
static Expr* wrap_not(const Expr* d) {
    Expr* a[1] = { expr_copy(d) };
    return expr_new_function(expr_new_symbol(SYM_Not), a, 1);
}

/* The folded logical negation of a *display* atom: fold Not into the relation
 * head where one exists (so output carries no stray Not on a relation), else
 * wrap in Not.  Only binary relations are folded; higher-arity Equal/Unequal
 * etc. are wrapped. */
static Expr* logical_negate(const Expr* d) {
    if (d->type == EXPR_SYMBOL) {
        if (d->data.symbol.name == SYM_True)  return expr_new_symbol(SYM_False);
        if (d->data.symbol.name == SYM_False) return expr_new_symbol(SYM_True);
        return wrap_not(d);
    }
    if (d->type == EXPR_FUNCTION && d->data.function.head
        && d->data.function.head->type == EXPR_SYMBOL) {
        const char* h = d->data.function.head->data.symbol.name;
        if (h == SYM_Not && nargs(d) == 1) return expr_copy(argn(d, 0));
        if (nargs(d) == 2) {
            if (h == SYM_Equal)        return relhead_swap(d, SYM_Unequal);
            if (h == SYM_Unequal)      return relhead_swap(d, SYM_Equal);
            if (h == SYM_Less)         return relhead_swap(d, SYM_GreaterEqual);
            if (h == SYM_GreaterEqual) return relhead_swap(d, SYM_Less);
            if (h == SYM_Greater)      return relhead_swap(d, SYM_LessEqual);
            if (h == SYM_LessEqual)    return relhead_swap(d, SYM_Greater);
            if (h == SYM_Element)      return relhead_swap(d, SYM_NotElement);
            if (h == SYM_NotElement)   return relhead_swap(d, SYM_Element);
        }
    }
    return wrap_not(d);
}

/* ------------------------------------------------------------------ *
 *  DNF data model                                                     *
 * ------------------------------------------------------------------ */

typedef struct { Expr** lit; int n, cap; } Clause; /* AND of display literals   */
typedef struct { Clause* cl; int n, cap; } Dnf;    /* OR of clauses             */
/* False = zero clauses; True = a single clause with zero literals.             */

static void clause_free(Clause* c) {
    for (int i = 0; i < c->n; i++) expr_free(c->lit[i]);
    free(c->lit);
    c->lit = NULL; c->n = c->cap = 0;
}
static void dnf_free(Dnf* d) {
    for (int i = 0; i < d->n; i++) clause_free(&d->cl[i]);
    free(d->cl);
    d->cl = NULL; d->n = d->cap = 0;
}
static void dnf_push(Dnf* d, Clause c) {
    if (d->n == d->cap) { d->cap = d->cap ? d->cap * 2 : 4;
                          d->cl = realloc(d->cl, sizeof(Clause) * d->cap); }
    d->cl[d->n++] = c;
}
static Dnf dnf_false(void) { Dnf d = { NULL, 0, 0 }; return d; }
static Dnf dnf_true(void)  { Dnf d = { NULL, 0, 0 }; Clause e = { NULL, 0, 0 }; dnf_push(&d, e); return d; }
static bool dnf_has_empty(const Dnf* d) { for (int i = 0; i < d->n; i++) if (d->cl[i].n == 0) return true; return false; }

static bool lits_equal(const Expr* p, const Expr* q) { return expr_eq(p, q); }
static bool lits_complementary(const Expr* p, const Expr* q) {
    Expr* nq = logical_negate(q);
    bool r = expr_eq(p, nq);
    expr_free(nq);
    return r;
}

/* Add `lit` (owned) to clause: -1 = complementary pair (lit freed, clause dead),
 * 0 = duplicate (lit freed), 1 = appended. */
static int clause_add(Clause* c, Expr* lit) {
    for (int i = 0; i < c->n; i++) {
        if (lits_complementary(c->lit[i], lit)) { expr_free(lit); return -1; }
        if (lits_equal(c->lit[i], lit))         { expr_free(lit); return 0; }
    }
    if (c->n == c->cap) { c->cap = c->cap ? c->cap * 2 : 4;
                          c->lit = realloc(c->lit, sizeof(Expr*) * c->cap); }
    c->lit[c->n++] = lit;
    return 1;
}
/* Merge two clauses into `out`; false = contradictory (out freed). */
static bool clause_merge(const Clause* a, const Clause* b, Clause* out) {
    out->lit = NULL; out->n = out->cap = 0;
    for (int i = 0; i < a->n; i++)
        if (clause_add(out, expr_copy(a->lit[i])) < 0) { clause_free(out); return false; }
    for (int i = 0; i < b->n; i++)
        if (clause_add(out, expr_copy(b->lit[i])) < 0) { clause_free(out); return false; }
    return true;
}
static bool clause_subset(const Clause* A, const Clause* B) { /* every A-lit in B */
    for (int i = 0; i < A->n; i++) {
        bool found = false;
        for (int j = 0; j < B->n; j++) if (lits_equal(A->lit[i], B->lit[j])) { found = true; break; }
        if (!found) return false;
    }
    return true;
}

/* Dedup + absorption: an empty clause makes the whole DNF True; otherwise drop
 * any clause that is a (proper or duplicate) superset of another. */
static void dnf_absorb(Dnf* d) {
    for (int i = 0; i < d->n; i++) if (d->cl[i].n == 0) {   /* True absorbs all */
        for (int k = 0; k < d->n; k++) clause_free(&d->cl[k]);
        d->n = 0;
        Clause e = { NULL, 0, 0 };
        dnf_push(d, e);
        return;
    }
    bool* rm = calloc(d->n, sizeof(bool));
    for (int i = 0; i < d->n; i++) {
        if (rm[i]) continue;
        for (int j = 0; j < d->n; j++) {
            if (i == j || rm[j]) continue;
            if (clause_subset(&d->cl[j], &d->cl[i]) &&
                (d->cl[j].n < d->cl[i].n || (d->cl[j].n == d->cl[i].n && j < i))) {
                rm[i] = true; break;
            }
        }
    }
    int w = 0;
    for (int i = 0; i < d->n; i++) { if (rm[i]) clause_free(&d->cl[i]); else d->cl[w++] = d->cl[i]; }
    d->n = w;
    free(rm);
}

/* Distributive product; consumes a and b. */
static Dnf dnf_and(Dnf a, Dnf b) {
    Dnf r = dnf_false();
    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < b.n; j++) {
            Clause m;
            if (clause_merge(&a.cl[i], &b.cl[j], &m)) dnf_push(&r, m);
        }
    dnf_free(&a); dnf_free(&b);
    dnf_absorb(&r);
    return r;
}
/* Disjunction; consumes a and b (their clauses are moved into r). */
static Dnf dnf_or(Dnf a, Dnf b) {
    Dnf r = dnf_false();
    for (int i = 0; i < a.n; i++) dnf_push(&r, a.cl[i]);
    for (int j = 0; j < b.n; j++) dnf_push(&r, b.cl[j]);
    free(a.cl); free(b.cl);
    dnf_absorb(&r);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Recursive DNF builder                                              *
 * ------------------------------------------------------------------ */

static Dnf to_dnf(const Expr* e);
static Dnf to_dnf_neg(const Expr* e);   /* DNF of Not[e] */

static Dnf leaf_dnf(const Expr* e, bool neg) {
    Dnf d = dnf_false();
    Clause c = { NULL, 0, 0 };
    clause_add(&c, neg ? logical_negate(e) : expr_copy(e));
    dnf_push(&d, c);
    return d;
}

static bool is_container(const Expr* c) {
    return c && c->type == EXPR_FUNCTION && c->data.function.head
        && c->data.function.head->type == EXPR_SYMBOL
        && (c->data.function.head->data.symbol.name == SYM_Alternatives
         || c->data.function.head->data.symbol.name == SYM_List);
}
static void collect_container(const Expr* c, Expr*** out, int* n, int* cap) {
    if (is_container(c)) { for (size_t i = 0; i < nargs(c); i++) collect_container(argn(c, i), out, n, cap); return; }
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 4; *out = realloc(*out, sizeof(Expr*) * (*cap)); }
    (*out)[(*n)++] = (Expr*)c;
}

/* Element[container, dom] / NotElement[container, dom] with a container arg:
 * positive membership => AND of Element[leaf, dom]; negated => OR of
 * NotElement[leaf, dom]. */
static Dnf element_multi_dnf(const Expr* elem, bool eff_not) {
    const Expr* dom = argn(elem, 1);
    Expr** leaves = NULL; int n = 0, cap = 0;
    collect_container(argn(elem, 0), &leaves, &n, &cap);
    Dnf d = eff_not ? dnf_false() : dnf_true();
    for (int i = 0; i < n; i++) {
        Expr* ea[2] = { expr_copy(leaves[i]), expr_copy(dom) };
        Expr* el = expr_new_function(expr_new_symbol(SYM_Element), ea, 2);
        Dnf leafd = leaf_dnf(el, eff_not);
        expr_free(el);
        d = eff_not ? dnf_or(d, leafd) : dnf_and(d, leafd);
    }
    free(leaves);
    return d;
}

/* Implies[a1, ..., an] = a1 -> (a2 -> ... -> an), right-associative. */
static Dnf implies_dnf(const Expr* e, bool neg) {
    size_t n = nargs(e);
    const Expr* a1 = argn(e, 0);
    Expr* built = NULL;
    const Expr* C;
    if (n == 2) C = argn(e, 1);
    else {
        Expr** ca = malloc(sizeof(Expr*) * (n - 1));
        for (size_t i = 1; i < n; i++) ca[i - 1] = expr_copy(argn(e, i));
        built = expr_new_function(expr_new_symbol(SYM_Implies), ca, n - 1);
        free(ca);
        C = built;
    }
    Dnf d = neg ? dnf_and(to_dnf(a1), to_dnf_neg(C))    /* !(a1 -> C) = a1 && !C */
                : dnf_or(to_dnf_neg(a1), to_dnf(C));     /*  (a1 -> C) = !a1 || C */
    if (built) expr_free(built);
    return d;
}

/* DNF of (neg ? Not[Xor[args[i..n-1]]] : Xor[args[i..n-1]]). */
static Dnf xor_slice(Expr* const* args, int i, int n, bool neg) {
    int len = n - i;
    if (len <= 0) return neg ? dnf_true() : dnf_false();   /* Xor[] = False */
    if (len == 1) return neg ? to_dnf_neg(args[i]) : to_dnf(args[i]);
    const Expr* a = args[i];
    if (!neg) {  /* Xor[a, R] = (a && !R) || (!a && R) */
        Dnf t1 = dnf_and(to_dnf(a),     xor_slice(args, i + 1, n, true));
        Dnf t2 = dnf_and(to_dnf_neg(a), xor_slice(args, i + 1, n, false));
        return dnf_or(t1, t2);
    }             /* !Xor[a, R] = (a && R) || (!a && !R) */
    Dnf t1 = dnf_and(to_dnf(a),     xor_slice(args, i + 1, n, false));
    Dnf t2 = dnf_and(to_dnf_neg(a), xor_slice(args, i + 1, n, true));
    return dnf_or(t1, t2);
}

static Dnf to_dnf(const Expr* e) {
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_True)  return dnf_true();
        if (e->data.symbol.name == SYM_False) return dnf_false();
        return leaf_dnf(e, false);
    }
    if (e->type != EXPR_FUNCTION) return leaf_dnf(e, false);
    if (is_head_sym(e, SYM_And)) { Dnf d = dnf_true();  for (size_t i = 0; i < nargs(e); i++) d = dnf_and(d, to_dnf(argn(e, i))); return d; }
    if (is_head_sym(e, SYM_Or))  { Dnf d = dnf_false(); for (size_t i = 0; i < nargs(e); i++) d = dnf_or(d, to_dnf(argn(e, i)));  return d; }
    if (is_head_sym(e, SYM_Not) && nargs(e) == 1)     return to_dnf_neg(argn(e, 0));
    if (is_head_sym(e, SYM_Implies) && nargs(e) >= 2) return implies_dnf(e, false);
    if (is_head_sym(e, SYM_Xor))                      return xor_slice(e->data.function.args, 0, (int)nargs(e), false);
    if ((is_head_sym(e, SYM_Element) || is_head_sym(e, SYM_NotElement)) && nargs(e) == 2 && is_container(argn(e, 0)))
        return element_multi_dnf(e, is_head_sym(e, SYM_NotElement));
    return leaf_dnf(e, false);
}

static Dnf to_dnf_neg(const Expr* e) {
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_True)  return dnf_false();
        if (e->data.symbol.name == SYM_False) return dnf_true();
        return leaf_dnf(e, true);
    }
    if (e->type != EXPR_FUNCTION) return leaf_dnf(e, true);
    if (is_head_sym(e, SYM_And)) { Dnf d = dnf_false(); for (size_t i = 0; i < nargs(e); i++) d = dnf_or(d, to_dnf_neg(argn(e, i))); return d; }   /* !(&&) = ||!  */
    if (is_head_sym(e, SYM_Or))  { Dnf d = dnf_true();  for (size_t i = 0; i < nargs(e); i++) d = dnf_and(d, to_dnf_neg(argn(e, i))); return d; }  /* !(||) = &&!  */
    if (is_head_sym(e, SYM_Not) && nargs(e) == 1)     return to_dnf(argn(e, 0));                 /* !!a = a */
    if (is_head_sym(e, SYM_Implies) && nargs(e) >= 2) return implies_dnf(e, true);
    if (is_head_sym(e, SYM_Xor))                      return xor_slice(e->data.function.args, 0, (int)nargs(e), true);
    if ((is_head_sym(e, SYM_Element) || is_head_sym(e, SYM_NotElement)) && nargs(e) == 2 && is_container(argn(e, 0)))
        return element_multi_dnf(e, !is_head_sym(e, SYM_NotElement));
    return leaf_dnf(e, true);
}

/* ------------------------------------------------------------------ *
 *  Emission                                                           *
 * ------------------------------------------------------------------ */

static Expr* clause_to_expr(const Clause* c) {
    if (c->n == 1) return expr_copy(c->lit[0]);
    Expr** a = malloc(sizeof(Expr*) * c->n);
    for (int i = 0; i < c->n; i++) a[i] = expr_copy(c->lit[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_And), a, (size_t)c->n);
    free(a);
    return r;
}
static Expr* dnf_to_expr(const Dnf* d) {
    if (d->n == 1) return clause_to_expr(&d->cl[0]);
    Expr** a = malloc(sizeof(Expr*) * d->n);
    for (int i = 0; i < d->n; i++) a[i] = clause_to_expr(&d->cl[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Or), a, (size_t)d->n);
    free(a);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Builtins                                                           *
 * ------------------------------------------------------------------ */

Expr* builtin_logical_expand(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* e = res->data.function.args[0];

    Dnf phi = to_dnf(e);
    Expr* out;
    if (phi.n == 0) {
        out = expr_new_symbol(SYM_False);                 /* contradiction */
    } else if (dnf_has_empty(&phi)) {
        out = expr_new_symbol(SYM_True);                  /* structural True */
    } else {
        Dnf nphi = to_dnf_neg(e);
        bool taut = (nphi.n == 0);                        /* !e unsatisfiable */
        dnf_free(&nphi);
        if (taut) {
            out = expr_new_symbol(SYM_True);
        } else {
            Expr* pre = dnf_to_expr(&phi);
            out = evaluate(pre);                          /* canonicalise And/Or */
            expr_free(pre);
        }
    }
    dnf_free(&phi);
    return out;
}

Expr* builtin_not_element(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* x   = res->data.function.args[0];
    Expr* dom = res->data.function.args[1];

    Expr* ea[2] = { expr_copy(x), expr_copy(dom) };
    Expr* el = expr_new_function(expr_new_symbol(SYM_Element), ea, 2);
    Expr* ev = evaluate(el);
    expr_free(el);

    Expr* out = NULL;
    if (ev && ev->type == EXPR_SYMBOL) {
        if (ev->data.symbol.name == SYM_True)  out = expr_new_symbol(SYM_False);
        else if (ev->data.symbol.name == SYM_False) out = expr_new_symbol(SYM_True);
    }
    if (ev) expr_free(ev);
    return out;   /* NULL keeps NotElement[x, dom] symbolic */
}

void reduce_companions_init(void) {
    symtab_add_builtin("LogicalExpand", builtin_logical_expand);
    symtab_add_builtin("NotElement",    builtin_not_element);

    SymbolDef* d;
    d = symtab_get_def("LogicalExpand"); if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("NotElement");    if (d) d->attributes |= ATTR_PROTECTED;

    symtab_set_docstring("LogicalExpand",
        "LogicalExpand[expr]\n"
        "\tExpands the logical combination expr -- of equations, inequalities\n"
        "\tand Boolean atoms -- into disjunctive normal form (an Or of Ands),\n"
        "\tapplying distributive, De Morgan, idempotence, complementation and\n"
        "\tabsorption laws, and expanding Implies and Xor.  Returns True for a\n"
        "\ttautology and False for a contradiction.  Every non-logical\n"
        "\tsubexpression is treated as an opaque Boolean atom (no domain\n"
        "\treasoning), so e.g. x==a and x!=a are complementary literals.");
    symtab_set_docstring("NotElement",
        "NotElement[x, dom]\n"
        "\tThe statement that x is not an element of the domain dom -- the\n"
        "\tnegation of Element[x, dom].  Decides to True or False when the\n"
        "\tmembership decides, and stays symbolic otherwise.");
}
