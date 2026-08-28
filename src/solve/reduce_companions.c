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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "expr.h"
#include "eval.h"
#include "attr.h"
#include "symtab.h"
#include "sym_names.h"
#include "message.h"           /* mth_msg_suppress_push/pop: quiet internal probes */
#include "reduce_real_util.h"   /* rru_rational_between, rru_sign_of, rru_approx_double */

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
/* expr_copy through a const pointer -- the copy (a refcount bump) never mutates
 * its source, so the cast is safe and keeps the read-only helpers const. */
static Expr* xcopy(const Expr* e) { return expr_copy((Expr*)e); }

/* Build a fresh 2-arg relation `newhead[d.arg0, d.arg1]`. */
static Expr* relhead_swap(const Expr* d, const char* newhead) {
    Expr* a[2] = { expr_copy(argn(d, 0)), expr_copy(argn(d, 1)) };
    return expr_new_function(expr_new_symbol(newhead), a, 2);
}
static Expr* wrap_not(const Expr* d) {
    Expr* a[1] = { xcopy(d) };
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
    bool* rm = calloc((size_t)(d->n > 0 ? d->n : 0), sizeof(bool));
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
    clause_add(&c, neg ? logical_negate(e) : xcopy(e));
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
        Expr* ea[2] = { xcopy(leaves[i]), xcopy(dom) };
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

/* Equivalent[a1, ..., an] -- every argument shares one truth value.  Rewritten to
 * the cyclic conjunction And[Implies[a1,a2], ..., Implies[an,a1]] and recursed. */
static Dnf equivalent_dnf(const Expr* e, bool neg) {
    size_t n = nargs(e);
    if (n <= 1) return neg ? dnf_false() : dnf_true();   /* trivially equivalent */
    Expr** imps = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* ia[2] = { xcopy(argn(e, i)), xcopy(argn(e, (i + 1) % n)) };
        imps[i] = expr_new_function(expr_new_symbol(SYM_Implies), ia, 2);
    }
    Expr* conj = expr_new_function(expr_new_symbol(SYM_And), imps, n);
    free(imps);
    Dnf d = neg ? to_dnf_neg(conj) : to_dnf(conj);
    expr_free(conj);
    return d;
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
    if (is_head_sym(e, SYM_Equivalent))               return equivalent_dnf(e, false);
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
    if (is_head_sym(e, SYM_Equivalent))               return equivalent_dnf(e, true);
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

/* ================================================================== *
 *  FindInstance  (REDUCE_PLAN.md, Phase 8)                            *
 *                                                                     *
 *  Find up to n witness points that make `expr` True over a domain,   *
 *  returned in Solve's form: {{x->v1, ...}, ...}, or {} when the set  *
 *  is provably empty, or unevaluated when a witness can neither be    *
 *  produced nor emptiness proved.                                     *
 *                                                                     *
 *  Strategy (soundness-first, maximal reuse): witnesses are extracted *
 *  from the PUBLIC cylindrical outputs of `Reduce` (the satisfiability *
 *  oracle) and `Solve` (parametric fallback), interval samples come   *
 *  from rru_rational_between, and EVERY candidate is verified against *
 *  the original `expr` (expr /. point === True) -- so a verified point *
 *  is never wrong, and an unprovable one declines.  The Booleans      *
 *  domain reuses the to_dnf engine above for satisfiability.          *
 * ================================================================== */

#include <math.h>   /* ceil, floor for the integer sampler */

static bool fi_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* Evaluate a freshly-built call and free the call tree; returns owned result. */
static Expr* fi_eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ReplaceAll[e, rules] evaluated; e and rules borrowed, result owned. */
static Expr* fi_replace(const Expr* e, const Expr* rules) {
    Expr* a[2] = { xcopy(e), xcopy(rules) };
    return fi_eval_take(expr_new_function(expr_new_symbol(SYM_ReplaceAll), a, 2));
}

/* Rule[copy(var), val] -- takes ownership of val. */
static Expr* fi_rule(const Expr* var, Expr* val) {
    Expr* a[2] = { xcopy(var), val };
    return expr_new_function(expr_new_symbol(SYM_Rule), a, 2);
}

static Expr* fi_empty_list(void) {
    return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
}

/* A "variable" may be a plain symbol (x) or an indexed form (c[1]); both are
 * matched structurally with expr_eq so FindInstance accepts either. */
static bool fi_is_var(const Expr* e, const Expr* v) {
    return e && v && expr_eq((Expr*)e, (Expr*)v);
}
static int fi_var_index(const Expr* e, Expr** V, int nv) {
    if (!e) return -1;
    for (int i = 0; i < nv; i++)
        if (expr_eq((Expr*)e, V[i])) return i;
    return -1;
}
/* Does any listed variable appear anywhere in e?  A listed var is tested as a
 * whole subtree (so `c[1]` matches, while the bare head `c` does not). */
static bool fi_contains_var(const Expr* e, Expr** V, int nv) {
    if (!e) return false;
    if (fi_var_index(e, V, nv) >= 0) return true;
    if (e->type == EXPR_FUNCTION) {
        if (fi_contains_var(e->data.function.head, V, nv)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (fi_contains_var(e->data.function.args[i], V, nv)) return true;
    }
    return false;
}
/* Does the single symbol v appear anywhere in e? */
static bool fi_contains_one(const Expr* e, const Expr* v) {
    Expr* one = (Expr*)v;
    return fi_contains_var(e, &one, 1);
}

/* Build a Reduce/Solve call head[expr, vars, dom?, Modulus->m?] and evaluate. */
static Expr* fi_call(const char* head, const Expr* a0, const Expr* a1,
                     const Expr* dom, const Expr* modulus) {
    int total = 2 + (dom ? 1 : 0) + (modulus ? 1 : 0);
    Expr** a = malloc(sizeof(Expr*) * (size_t)total);
    int k = 0;
    a[k++] = xcopy(a0);
    a[k++] = xcopy(a1);
    if (dom) a[k++] = xcopy(dom);
    if (modulus) {
        Expr* r[2] = { expr_new_symbol(SYM_Modulus), xcopy(modulus) };
        a[k++] = expr_new_function(expr_new_symbol(SYM_Rule), r, 2);
    }
    Expr* call = expr_new_function(expr_new_symbol(head), a, (size_t)total);
    free(a);
    return fi_eval_take(call);
}

/* True if e contains an infinity / indeterminate sentinel anywhere -- the marks
 * of an UNDEFINED evaluation (Log[0] -> -Infinity = Times[-1, Infinity],
 * 1/0 -> ComplexInfinity, 0/0 -> Indeterminate).  A witness point that drives a
 * relation operand to one of these is NOT a true instance: a fold such as
 * -Infinity == -Infinity is undefined, not satisfied. */
static bool fi_has_nonfinite(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        return n == SYM_Infinity || n == SYM_ComplexInfinity
            || n == SYM_Indeterminate || n == SYM_Undefined
            || n == SYM_Overflow;
    }
    if (e->type == EXPR_FUNCTION) {
        if (is_head_sym(e, SYM_DirectedInfinity)) return true;
        if (fi_has_nonfinite(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (fi_has_nonfinite(e->data.function.args[i])) return true;
    }
    return false;
}

/* Substitute one operand at the point and report whether it is DEFINED (no
 * infinity/indeterminate sentinel). */
static bool fi_operand_defined(const Expr* operand, const Expr* pt) {
    Expr* s = fi_replace(operand, pt);
    bool bad = fi_has_nonfinite(s);
    expr_free(s);
    return !bad;
}

/* Three-valued, definedness-aware truth of a (possibly compound) statement at
 * the point: 1 = defined and true, 0 = defined and false, -1 = undefined or
 * undecided.  Mirrors the logical structure of fi_num_true but decides EXACTLY
 * via the evaluator, and rejects any relation whose substituted operand is a
 * non-finite sentinel (so z=0 is not accepted for Log[z^2] == 2 Log[z] + 2 Pi I,
 * where Log[0] = -Infinity makes the equation fold spuriously to True). */
static int fi_defined_truth(const Expr* node, const Expr* pt) {
    if (fi_is_sym(node, SYM_True))  return 1;
    if (fi_is_sym(node, SYM_False)) return 0;
    if (is_head_sym(node, SYM_And)) {
        int r = 1;
        for (size_t i = 0; i < nargs(node); i++) {
            int t = fi_defined_truth(argn(node, i), pt);
            if (t == 0) return 0;
            if (t < 0) r = -1;
        }
        return r;
    }
    if (is_head_sym(node, SYM_Or)) {
        int r = 0;
        for (size_t i = 0; i < nargs(node); i++) {
            int t = fi_defined_truth(argn(node, i), pt);
            if (t == 1) return 1;
            if (t < 0) r = -1;
        }
        return r;
    }
    if (is_head_sym(node, SYM_Not) && nargs(node) == 1) {
        int t = fi_defined_truth(argn(node, 0), pt);
        return t < 0 ? -1 : !t;
    }
    /* Binary relation: both operands must be defined before the relation decides. */
    if (node->type == EXPR_FUNCTION && node->data.function.head->type == EXPR_SYMBOL
        && nargs(node) == 2) {
        const char* h = node->data.function.head->data.symbol.name;
        if (h == SYM_Equal || h == SYM_Unequal || h == SYM_Less || h == SYM_Greater
            || h == SYM_LessEqual || h == SYM_GreaterEqual) {
            if (!fi_operand_defined(argn(node, 0), pt)
                || !fi_operand_defined(argn(node, 1), pt)) return -1;
            Expr* sub = fi_replace(node, pt);
            int r = fi_is_sym(sub, SYM_True) ? 1 : (fi_is_sym(sub, SYM_False) ? 0 : -1);
            expr_free(sub);
            return r;
        }
    }
    /* Chained Inequality[a, op, b, op, c, ...]: every operand must be defined. */
    if (is_head_sym(node, SYM_Inequality) && (nargs(node) & 1)) {
        for (size_t i = 0; i < nargs(node); i += 2)
            if (!fi_operand_defined(argn(node, i), pt)) return -1;
        Expr* sub = fi_replace(node, pt);
        int r = fi_is_sym(sub, SYM_True) ? 1 : (fi_is_sym(sub, SYM_False) ? 0 : -1);
        expr_free(sub);
        return r;
    }
    /* Unknown atom (bare Boolean, Element[..], ...): evaluate as a whole and
     * require a defined True. */
    Expr* sub = fi_replace(node, pt);
    int r = (!fi_has_nonfinite(sub) && fi_is_sym(sub, SYM_True)) ? 1
          : (fi_is_sym(sub, SYM_False) ? 0 : -1);
    expr_free(sub);
    return r;
}

/* The single soundness gate: does expr hold -- and is it DEFINED -- at the
 * candidate point?  Definedness matters because the evaluator folds arithmetic
 * on infinities (Log[0] = -Infinity, -Infinity == -Infinity -> True), which
 * would otherwise admit a point where the statement is undefined. */
static bool fi_verify(const Expr* expr, const Expr* point /* List[Rule..] */) {
    return fi_defined_truth(expr, point) == 1;
}

/* ---- witness accumulator (dedup on FullForm equality) ---------------------- */

typedef struct { Expr** p; int n, cap; } FiWit;
static void fi_wit_add(FiWit* w, Expr* pt) {
    for (int i = 0; i < w->n; i++)
        if (expr_eq(w->p[i], pt)) { expr_free(pt); return; }
    if (w->n == w->cap) { w->cap = w->cap ? w->cap * 2 : 4;
                          w->p = realloc(w->p, sizeof(Expr*) * (size_t)w->cap); }
    w->p[w->n++] = pt;
}
/* Move the first `want` witnesses into a List; free the rest; frees w->p. */
static Expr* fi_wit_take(FiWit* w, long want) {
    int n = (want < w->n) ? (int)want : w->n;
    Expr** a = (n > 0) ? malloc(sizeof(Expr*) * (size_t)n) : NULL;
    for (int i = 0; i < n; i++) a[i] = w->p[i];
    for (int i = n; i < w->n; i++) expr_free(w->p[i]);
    free(w->p);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), a, (size_t)n);
    free(a);
    return r;
}
static void fi_wit_free(FiWit* w) {
    for (int i = 0; i < w->n; i++) expr_free(w->p[i]);
    free(w->p);
}

/* ---- interval sampling ----------------------------------------------------- */

static bool fi_forbidden(const Expr* s, Expr** forb, int nforb) {
    for (int i = 0; i < nforb; i++)
        if (rru_sign_compare(s, forb[i]) == 0) return true;
    return false;
}

/* An integer in [lo,hi] (NULL = unbounded), honouring strictness / exclusions. */
static Expr* fi_sample_int(const Expr* lo, bool lo_strict, const Expr* hi, bool hi_strict,
                           Expr** forb, int nforb) {
    long start; int dir;
    bool okl = false, okh = false; double dl = 0, dh = 0;
    if (lo) dl = rru_approx_double(lo, &okl);
    if (hi) dh = rru_approx_double(hi, &okh);
    if (lo && okl)      { start = (long)floor(dl) - 2; dir = +1; }
    else if (hi && okh) { start = (long)ceil(dh)  + 2; dir = -1; }
    else                { start = 0; dir = +1; }
    for (int step = 0; step < 4000; step++) {
        long cand = start + (long)dir * step;
        Expr* c = expr_new_integer(cand);
        bool good = true;
        if (lo) { int sc = rru_sign_compare(c, lo); good = lo_strict ? (sc > 0) : (sc >= 0); }
        if (good && hi) { int sc = rru_sign_compare(c, hi); good = hi_strict ? (sc < 0) : (sc <= 0); }
        if (good && fi_forbidden(c, forb, nforb)) good = false;
        if (good) return c;
        expr_free(c);
        /* once we have passed the far bound, stop */
        if (dir > 0 && hi && okh && cand > dh + 2) break;
        if (dir < 0 && lo && okl && cand < dl - 2) break;
    }
    return NULL;
}

/* A sample value for a free variable given collected bounds/exclusions. */
static Expr* fi_sample(const Expr* lo, bool lo_strict, const Expr* hi, bool hi_strict,
                       Expr** forb, int nforb, const Expr* dom) {
    if (fi_is_sym(dom, SYM_Integers))
        return fi_sample_int(lo, lo_strict, hi, hi_strict, forb, nforb);
    if (!lo && !hi) {
        for (int k = 0; k <= 20; k++) {
            long cand = (k == 0) ? 0 : ((k & 1) ? (k + 1) / 2 : -(k / 2));
            Expr* c = expr_new_integer(cand);
            if (!fi_forbidden(c, forb, nforb)) return c;
            expr_free(c);
        }
        return NULL;
    }
    Expr* s = rru_rational_between(lo, hi);   /* strict interior, certified */
    if (!s) return NULL;
    if (!fi_forbidden(s, forb, nforb)) return s;
    /* nudge off an excluded point */
    Expr* s2 = rru_rational_between(s, hi);
    if (s2 && !fi_forbidden(s2, forb, nforb)) { expr_free(s); return s2; }
    if (s2) expr_free(s2);
    Expr* s3 = rru_rational_between(lo, s);
    if (s3 && !fi_forbidden(s3, forb, nforb)) { expr_free(s); return s3; }
    if (s3) expr_free(s3);
    expr_free(s);
    return NULL;
}

/* Tighten a lower bound: keep the larger of *lo and cand (borrowed). */
static void fi_tighten_lo(Expr** lo, bool* strict, const Expr* cand, bool cand_strict) {
    if (!*lo) { *lo = (Expr*)cand; *strict = cand_strict; return; }
    int sc = rru_sign_compare(cand, *lo);
    if (sc > 0 || (sc == 0 && cand_strict)) { *lo = (Expr*)cand; *strict = cand_strict; }
}
static void fi_tighten_hi(Expr** hi, bool* strict, const Expr* cand, bool cand_strict) {
    if (!*hi) { *hi = (Expr*)cand; *strict = cand_strict; return; }
    int sc = rru_sign_compare(cand, *hi);
    if (sc < 0 || (sc == 0 && cand_strict)) { *hi = (Expr*)cand; *strict = cand_strict; }
}

/* One relational operand pair `L op R` where exactly one side is the free var v
 * and the other is constant.  Record the implied bound / exclusion. */
static void fi_bound_from(const Expr* L, const char* op, const Expr* R, const Expr* v,
                          Expr** V, int nv,
                          Expr** lo, bool* lostrict, Expr** hi, bool* histrict,
                          Expr*** forb, int* nforb, int* cforb) {
    const Expr* other; bool v_left;
    if (fi_is_var(L, v) && !fi_contains_var(R, V, nv)) { other = R; v_left = true; }
    else if (fi_is_var(R, v) && !fi_contains_var(L, V, nv)) { other = L; v_left = false; }
    else return;
    if (op == SYM_Unequal) {
        if (*nforb == *cforb) { *cforb = *cforb ? *cforb * 2 : 4;
                                *forb = realloc(*forb, sizeof(Expr*) * (size_t)*cforb); }
        (*forb)[(*nforb)++] = (Expr*)other;
        return;
    }
    /* normalise to a statement about v: (v_left) v op other  else  other op v */
    bool upper;   /* other bounds v from above? */
    bool strict;
    if (op == SYM_Less)              { upper =  v_left; strict = true;  }
    else if (op == SYM_LessEqual)    { upper =  v_left; strict = false; }
    else if (op == SYM_Greater)      { upper = !v_left; strict = true;  }
    else if (op == SYM_GreaterEqual) { upper = !v_left; strict = false; }
    else return;
    if (upper) fi_tighten_hi(hi, histrict, other, strict);
    else       fi_tighten_lo(lo, lostrict, other, strict);
}

/* Scan a (substituted, evaluated) node for bounds on the free variable v. */
static void fi_scan_bounds(const Expr* node, const Expr* v, Expr** V, int nv,
                           Expr** lo, bool* lostrict, Expr** hi, bool* histrict,
                           Expr*** forb, int* nforb, int* cforb) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (h->type != EXPR_SYMBOL) return;
    const char* hn = h->data.symbol.name;
    size_t n = node->data.function.arg_count;
    if (hn == SYM_And || hn == SYM_Or) {
        for (size_t i = 0; i < n; i++)
            fi_scan_bounds(node->data.function.args[i], v, V, nv,
                           lo, lostrict, hi, histrict, forb, nforb, cforb);
        return;
    }
    if ((hn == SYM_Less || hn == SYM_LessEqual || hn == SYM_Greater
         || hn == SYM_GreaterEqual || hn == SYM_Unequal) && n == 2) {
        fi_bound_from(node->data.function.args[0], hn, node->data.function.args[1], v,
                      V, nv, lo, lostrict, hi, histrict, forb, nforb, cforb);
        return;
    }
    if (hn == SYM_Inequality && n >= 3 && (n & 1)) {
        for (size_t i = 0; i + 2 < n; i += 2) {
            const Expr* L = node->data.function.args[i];
            const Expr* opE = node->data.function.args[i + 1];
            const Expr* R = node->data.function.args[i + 2];
            if (opE->type == EXPR_SYMBOL)
                fi_bound_from(L, opE->data.symbol.name, R, v, V, nv,
                              lo, lostrict, hi, histrict, forb, nforb, cforb);
        }
    }
}

/* ---- one verified point from a conjunction clause -------------------------- */

/* Extract Solve's first solution value for variable v: sol is
 * List[List[Rule[v, value], ...], ...].  Returns owned value, or NULL. */
static Expr* fi_first_value(const Expr* sol, const Expr* v, Expr** V, int nv) {
    if (!is_head_sym(sol, SYM_List) || nargs(sol) == 0) return NULL;
    const Expr* first = argn(sol, 0);
    if (!is_head_sym(first, SYM_List)) return NULL;
    for (size_t i = 0; i < nargs(first); i++) {
        const Expr* rule = argn(first, i);
        if (is_head_sym(rule, SYM_Rule) && nargs(rule) == 2 && fi_is_var(argn(rule, 0), v)) {
            const Expr* val = argn(rule, 1);
            if (fi_contains_var(val, V, nv)) return NULL;   /* parametric -> skip */
            return xcopy(val);
        }
    }
    return NULL;
}

/* Build List[Rule..] from currently-assigned vars + params, for substitution. */
static Expr* fi_rules_so_far(Expr** V, Expr** val, int nv, Expr** pk, Expr** pv, int np) {
    int cnt = np;
    for (int i = 0; i < nv; i++) if (val[i]) cnt++;
    Expr** a = malloc(sizeof(Expr*) * (size_t)(cnt > 0 ? cnt : 1));
    int k = 0;
    for (int i = 0; i < np; i++) a[k++] = fi_rule(pk[i], expr_copy(pv[i]));
    for (int i = 0; i < nv; i++) if (val[i]) a[k++] = fi_rule(V[i], expr_copy(val[i]));
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), a, (size_t)k);
    free(a);
    return lst;
}

static Expr* fi_clause_point(Expr** atoms, int na, Expr** V, int nv,
                             const Expr* dom, const Expr* modulus) {
    Expr** val = calloc((size_t)nv, sizeof(Expr*));   /* owned concrete values */
    bool* pintgt = calloc((size_t)nv, sizeof(bool));
    Expr** pk = NULL; Expr** pv = NULL; int np = 0, pcap = 0;
    Expr* result = NULL;

    /* parameters: Element[p, _] where p is NOT a listed var -> 0 */
    for (int i = 0; i < na; i++) {
        Expr* a = atoms[i];
        if (is_head_sym(a, SYM_Element) && nargs(a) == 2 && fi_var_index(argn(a, 0), V, nv) < 0
            && !fi_contains_var(argn(a, 0), V, nv)) {
            if (np == pcap) { pcap = pcap ? pcap * 2 : 4;
                             pk = realloc(pk, sizeof(Expr*) * (size_t)pcap);
                             pv = realloc(pv, sizeof(Expr*) * (size_t)pcap); }
            pk[np] = expr_copy(argn(a, 0));
            pv[np] = expr_new_integer(0);
            np++;
        }
    }
    /* pin targets: a listed var isolated on one side of an Equal atom */
    for (int i = 0; i < na; i++) {
        Expr* a = atoms[i];
        if (is_head_sym(a, SYM_Equal) && nargs(a) == 2) {
            int vi = fi_var_index(argn(a, 0), V, nv);
            if (vi < 0) vi = fi_var_index(argn(a, 1), V, nv);
            if (vi >= 0) pintgt[vi] = true;
        }
    }
    /* free vars first (not pin targets), in order: sample from bounds */
    for (int vi = 0; vi < nv; vi++) {
        if (pintgt[vi] || val[vi]) continue;
        /* substitute known values into each atom, scan the (evaluated) result for
         * bounds on this free var; the bound pointers borrow into `subs`, so keep
         * every substituted atom alive until after sampling. */
        Expr* rules = fi_rules_so_far(V, val, nv, pk, pv, np);
        Expr** subs = malloc(sizeof(Expr*) * (size_t)(na > 0 ? na : 1));
        Expr* lo = NULL, *hi = NULL; bool los = false, his = false;
        Expr** forb = NULL; int nforb = 0, cforb = 0;
        for (int i = 0; i < na; i++) {
            subs[i] = fi_replace(atoms[i], rules);
            fi_scan_bounds(subs[i], V[vi], V, nv, &lo, &los, &hi, &his, &forb, &nforb, &cforb);
        }
        Expr* sample = fi_sample(lo, los, hi, his, forb, nforb, dom);
        for (int i = 0; i < na; i++) expr_free(subs[i]);
        free(subs); free(forb); expr_free(rules);
        if (!sample) goto cleanup;   /* cannot sample this free var -> fail clause */
        val[vi] = sample;
    }
    /* resolve pin equations by substitution + Solve (fixpoint) */
    for (int iter = 0; iter <= nv; iter++) {
        bool progressed = false;
        bool all_pins = true;
        Expr* rules = fi_rules_so_far(V, val, nv, pk, pv, np);
        for (int i = 0; i < na; i++) {
            Expr* a = atoms[i];
            if (!(is_head_sym(a, SYM_Equal) && nargs(a) == 2)) continue;
            Expr* a2 = fi_replace(a, rules);
            if (fi_is_sym(a2, SYM_True)) { expr_free(a2); continue; }
            if (fi_is_sym(a2, SYM_False)) { expr_free(a2); expr_free(rules); goto cleanup; }
            int uidx = -1, ucount = 0;
            for (int k = 0; k < nv; k++)
                if (!val[k] && fi_contains_one(a2, V[k])) { ucount++; uidx = k; }
            if (ucount == 1) {
                Expr* sol = fi_call(SYM_Solve, a2, V[uidx], dom, modulus);
                Expr* value = fi_first_value(sol, V[uidx], V, nv);
                expr_free(sol);
                if (value) { val[uidx] = value; progressed = true;
                             /* refresh rules for subsequent atoms this pass */
                             expr_free(rules); rules = fi_rules_so_far(V, val, nv, pk, pv, np); }
            }
            expr_free(a2);
        }
        expr_free(rules);
        for (int k = 0; k < nv; k++) if (pintgt[k] && !val[k]) all_pins = false;
        if (all_pins) break;
        if (!progressed) break;
    }
    /* leftover unassigned vars -> 0 */
    for (int vi = 0; vi < nv; vi++) if (!val[vi]) val[vi] = expr_new_integer(0);

    /* build the point; every value must be free of listed vars */
    {
        Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
        bool ok = true;
        for (int vi = 0; vi < nv; vi++) {
            if (fi_contains_var(val[vi], V, nv)) ok = false;
            outr[vi] = fi_rule(V[vi], expr_copy(val[vi]));
        }
        if (ok) {
            result = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
            free(outr);
        } else {
            for (int vi = 0; vi < nv; vi++) expr_free(outr[vi]);
            free(outr);
        }
    }

cleanup:
    for (int i = 0; i < nv; i++) expr_free(val[i]);
    free(val); free(pintgt);
    for (int i = 0; i < np; i++) { expr_free(pk[i]); expr_free(pv[i]); }
    free(pk); free(pv);
    return result;
}

/* Instantiate a Solve rule-list into a full concrete point, free listed vars -> g. */
static Expr* fi_solve_point(const Expr* rulelist, Expr** V, int nv, long g) {
    if (!is_head_sym(rulelist, SYM_List)) return NULL;
    /* keys present in the rule-list */
    bool* iskey = calloc((size_t)nv, sizeof(bool));
    for (size_t i = 0; i < nargs(rulelist); i++) {
        const Expr* rule = argn(rulelist, i);
        if (is_head_sym(rule, SYM_Rule) && nargs(rule) == 2) {
            int vi = fi_var_index(argn(rule, 0), V, nv);
            if (vi >= 0) iskey[vi] = true;
        }
    }
    /* free-variable rules: every listed var not a key -> g */
    Expr** fr = malloc(sizeof(Expr*) * (size_t)(nv > 0 ? nv : 1));
    int nfr = 0;
    for (int vi = 0; vi < nv; vi++)
        if (!iskey[vi]) fr[nfr++] = fi_rule(V[vi], expr_new_integer(g));
    Expr* frlist = expr_new_function(expr_new_symbol(SYM_List), fr, (size_t)nfr);
    free(fr);

    Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
    bool ok = true;
    for (int vi = 0; vi < nv; vi++) {
        Expr* value = NULL;
        if (iskey[vi]) {
            /* find the rule and evaluate its RHS with free vars fixed */
            for (size_t i = 0; i < nargs(rulelist); i++) {
                const Expr* rule = argn(rulelist, i);
                if (is_head_sym(rule, SYM_Rule) && nargs(rule) == 2
                    && fi_var_index(argn(rule, 0), V, nv) == vi) {
                    value = fi_replace(argn(rule, 1), frlist);
                    break;
                }
            }
        } else {
            value = expr_new_integer(g);
        }
        if (!value || fi_contains_var(value, V, nv)) ok = false;
        outr[vi] = fi_rule(V[vi], value ? value : expr_new_integer(0));
    }
    expr_free(frlist); free(iskey);
    if (!ok) { for (int vi = 0; vi < nv; vi++) expr_free(outr[vi]); free(outr); return NULL; }
    Expr* p = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
    free(outr);
    return p;
}

/* ---- generated-parameter instantiation (parametric Diophantine / periodic) --- *
 *                                                                                *
 *  Reduce/Solve return parametric families carrying a generated parameter C[k]:  *
 *      x -> ConditionalExpression[value(C[1]), Element[C[1], Integers]]          *
 *  A witness needs a CONCRETE parameter.  We (1) try a small integer grid        *
 *  (fundamental solutions usually sit at C=1 -- e.g. Pell x^2-61 y^2==1), then   *
 *  (2) for a single parameter, solve the residual constraint for it over the     *
 *  Reals and pick an integer -- which reaches C[1]=15916 for the periodic        *
 *  Sin[1/x]==0 && 0<x<10^-5.  Every candidate is still verified.                 */

/* small expression builders (each consumes its Expr* arguments) */
static Expr* fi_bin(const char* head, Expr* a, Expr* b) {
    Expr* x[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), x, 2);
}
static Expr* fi_neg(Expr* a)                 { return fi_bin(SYM_Times, expr_new_integer(-1), a); }
static Expr* fi_sub_eval(Expr* a, Expr* b)   { return fi_eval_take(fi_bin(SYM_Plus, a, fi_neg(b))); }
static Expr* fi_add_eval(Expr* a, Expr* b)   { return fi_eval_take(fi_bin(SYM_Plus, a, b)); }
static Expr* fi_mul_eval(Expr* a, Expr* b)   { return fi_eval_take(fi_bin(SYM_Times, a, b)); }

/* ConditionalExpression[v, _] -> v, recursively (owned copy). */
static Expr* fi_strip_ce(const Expr* e) {
    if (!e) return NULL;
    if (is_head_sym(e, SYM_ConditionalExpression) && nargs(e) == 2)
        return fi_strip_ce(argn(e, 0));
    if (e->type != EXPR_FUNCTION) return xcopy(e);
    size_t n = e->data.function.arg_count;
    Expr* h = fi_strip_ce(e->data.function.head);
    Expr** a = n ? malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) a[i] = fi_strip_ce(e->data.function.args[i]);
    Expr* r = expr_new_function(h, a, n);
    free(a);
    return r;
}

static void fi_param_push(Expr*** ps, int* np, int* cp, const Expr* p) {
    for (int i = 0; i < *np; i++) if (expr_eq((*ps)[i], (Expr*)p)) return;
    if (*np == *cp) { *cp = *cp ? *cp * 2 : 4; *ps = realloc(*ps, sizeof(Expr*) * (size_t)*cp); }
    (*ps)[(*np)++] = xcopy(p);
}
/* A plausible generated parameter: a bare symbol (not a listed var, domain or Pi/
 * True/False), or an indexed form C[1] (symbol head, all-integer args). */
static bool fi_is_param_atom(const Expr* e, Expr** V, int nv) {
    if (fi_var_index(e, V, nv) >= 0) return false;
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        return !(n == SYM_Integers || n == SYM_Reals || n == SYM_Rationals
                 || n == SYM_Complexes || n == SYM_Booleans || n == SYM_Pi
                 || n == SYM_True || n == SYM_False);
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.arg_count == 0) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (e->data.function.args[i]->type != EXPR_INTEGER) return false;
        return true;
    }
    return false;
}
/* Collect parameter atoms from a ConditionalExpression condition / Element atom. */
static void fi_params_from_pred(const Expr* pred, Expr** V, int nv,
                                Expr*** ps, int* np, int* cp) {
    if (!pred || pred->type != EXPR_FUNCTION) return;
    const Expr* h = pred->data.function.head;
    if (h->type != EXPR_SYMBOL) return;
    const char* hn = h->data.symbol.name;
    if (hn == SYM_And || hn == SYM_Or || hn == SYM_Not) {
        for (size_t i = 0; i < nargs(pred); i++)
            fi_params_from_pred(argn(pred, i), V, nv, ps, np, cp);
        return;
    }
    if (hn == SYM_Element && nargs(pred) == 2) {
        if (fi_is_param_atom(argn(pred, 0), V, nv)) fi_param_push(ps, np, cp, argn(pred, 0));
        return;
    }
    if (hn == SYM_Less || hn == SYM_LessEqual || hn == SYM_Greater || hn == SYM_GreaterEqual
        || hn == SYM_Equal || hn == SYM_Unequal || hn == SYM_Inequality) {
        for (size_t i = 0; i < nargs(pred); i++)
            if (fi_is_param_atom(argn(pred, i), V, nv)) fi_param_push(ps, np, cp, argn(pred, i));
    }
}
/* Collect generated parameters referenced anywhere in `node` (via its
 * ConditionalExpression conditions and bare Element atoms). */
static void fi_gather_params(const Expr* node, Expr** V, int nv,
                             Expr*** ps, int* np, int* cp) {
    if (!node || node->type != EXPR_FUNCTION) return;
    if (is_head_sym(node, SYM_ConditionalExpression) && nargs(node) == 2)
        fi_params_from_pred(argn(node, 1), V, nv, ps, np, cp);
    if (is_head_sym(node, SYM_Element) && nargs(node) == 2
        && fi_is_param_atom(argn(node, 0), V, nv))
        fi_param_push(ps, np, cp, argn(node, 0));
    for (size_t i = 0; i < node->data.function.arg_count; i++)
        fi_gather_params(node->data.function.args[i], V, nv, ps, np, cp);
}

static const long FI_PGRID[] = { 1, -1, 2, -2, 3, -3, 0 };

/* Substitute params[i] -> vals[i] into rl and evaluate (collapses the
 * ConditionalExpression once its condition holds). */
static Expr* fi_subst_params(const Expr* rl, Expr** params, const long* vals, int np) {
    Expr** rr = malloc(sizeof(Expr*) * (size_t)np);
    for (int i = 0; i < np; i++) rr[i] = fi_rule(params[i], expr_new_integer(vals[i]));
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), rr, (size_t)np);
    free(rr);
    Expr* out = fi_replace(rl, rules);
    expr_free(rules);
    return out;
}
/* Recurse over a small integer grid for each parameter; for every assignment,
 * materialise concrete points from the substituted rule-list and verify. */
static void fi_param_combos(const Expr* rl, Expr** V, int nv, const Expr* origExpr,
                            Expr** params, int np, int gsz, long* chosen, int idx,
                            FiWit* ws, long want) {
    if (ws->n >= want) return;
    if (idx == np) {
        Expr* rl2 = fi_subst_params(rl, params, chosen, np);
        static const long g[] = { 0, 1, -1, 2, -2 };
        for (size_t gi = 0; gi < sizeof g / sizeof g[0] && ws->n < want; gi++) {
            Expr* p = fi_solve_point(rl2, V, nv, g[gi]);
            if (p) { if (fi_verify(origExpr, p)) fi_wit_add(ws, p); else expr_free(p); }
        }
        expr_free(rl2);
        return;
    }
    for (int k = 0; k < gsz && ws->n < want; k++) {
        chosen[idx] = FI_PGRID[k];
        fi_param_combos(rl, V, nv, origExpr, params, np, gsz, chosen, idx + 1, ws, want);
    }
}
/* Ceiling/Floor of a (possibly symbolic, e.g. 50000/Pi) real bound -> the
 * boundary integer, nudged one step in when the bound is a strict integer. */
static bool fi_bound_int(const Expr* b, bool strict, bool lower, long* out) {
    Expr* a[1] = { xcopy(b) };
    Expr* v = fi_eval_take(expr_new_function(
        expr_new_symbol(lower ? SYM_Ceiling : SYM_Floor), a, 1));
    bool ok = false;
    if (v->type == EXPR_INTEGER) {
        long g = v->data.integer;
        if (strict && rru_sign_compare(v, b) == 0) g += lower ? 1 : -1;
        *out = g; ok = true;
    }
    expr_free(v);
    return ok;
}
/* Single-parameter case that the grid missed: solve the residual constraint for
 * the parameter over the Reals, pick an integer, materialise and verify. */
static void fi_solve_param(const Expr* rl, Expr** V, int nv, const Expr* origExpr,
                           const Expr* Rconstr, const Expr* param, FiWit* ws, long want) {
    /* value-map: each key var -> strip_ce(RHS) (parameter retained); others -> 0 */
    Expr** rr = malloc(sizeof(Expr*) * (size_t)nv);
    for (int vi = 0; vi < nv; vi++) {
        Expr* val = NULL;
        for (size_t i = 0; i < nargs(rl); i++) {
            const Expr* rule = argn(rl, i);
            if (is_head_sym(rule, SYM_Rule) && nargs(rule) == 2
                && fi_var_index(argn(rule, 0), V, nv) == vi) { val = fi_strip_ce(argn(rule, 1)); break; }
        }
        rr[vi] = fi_rule(V[vi], val ? val : expr_new_integer(0));
    }
    Expr* vmap = expr_new_function(expr_new_symbol(SYM_List), rr, (size_t)nv);
    free(rr);

    const Expr* base = Rconstr ? Rconstr : origExpr;
    /* rename the (possibly indexed) parameter to a fresh plain symbol Reduce accepts */
    Expr* q = expr_new_symbol("FindInstance$param");
    Expr* qr[1] = { fi_rule(param, xcopy(q)) };
    Expr* qrules = expr_new_function(expr_new_symbol(SYM_List), qr, 1);
    Expr* pred0 = fi_replace(base, vmap);
    Expr* pred = fi_replace(pred0, qrules);
    expr_free(pred0); expr_free(qrules);

    Expr* reals = expr_new_symbol(SYM_Reals);
    Expr* redp = fi_call(SYM_Reduce, pred, q, reals, NULL);
    expr_free(reals); expr_free(pred);

    Expr* lo = NULL, *hi = NULL; bool los = false, his = false;
    Expr** forb = NULL; int nforb = 0, cforb = 0;
    Expr* qv[1] = { q };
    fi_scan_bounds(redp, q, qv, 1, &lo, &los, &hi, &his, &forb, &nforb, &cforb);
    /* Ceiling/Floor of the (possibly symbolic, e.g. 50000/Pi) bound gives the
     * boundary integer; walk inward from there. */
    long base_i = 0; int dir = 0;
    if (lo && fi_bound_int(lo, los, true, &base_i)) dir = +1;
    else if (hi && fi_bound_int(hi, his, false, &base_i)) dir = -1;
    if (dir != 0)
        for (long k = 0; k < want + 4 && ws->n < want; k++) {
            long vals[1] = { base_i + dir * k };
            Expr* one[1] = { (Expr*)param };
            Expr* rl2 = fi_subst_params(rl, one, vals, 1);
            Expr* p = fi_solve_point(rl2, V, nv, 0);
            expr_free(rl2);
            if (p) { if (fi_verify(origExpr, p)) fi_wit_add(ws, p); else expr_free(p); }
        }
    free(forb);
    expr_free(redp); expr_free(vmap); expr_free(q);
}
/* Produce verified points from a Solve rule-list that may carry parameters. */
static void fi_materialize_rulelist(const Expr* rl, Expr** V, int nv, const Expr* origExpr,
                                    const Expr* Rconstr, FiWit* ws, long want) {
    Expr** ps = NULL; int np = 0, cp = 0;
    fi_gather_params(rl, V, nv, &ps, &np, &cp);
    if (np == 0) {
        static const long grid[] = { 0, 1, -1, 2, -2, 3, -3 };
        for (size_t gi = 0; gi < sizeof grid / sizeof grid[0] && ws->n < want; gi++) {
            Expr* p = fi_solve_point(rl, V, nv, grid[gi]);
            if (p) { if (fi_verify(origExpr, p)) fi_wit_add(ws, p); else expr_free(p); }
        }
        free(ps);
        return;
    }
    int gsz = 7;
    while (gsz > 1) { double t = 1; for (int i = 0; i < np; i++) t *= gsz; if (t <= 64) break; gsz--; }
    long* chosen = malloc(sizeof(long) * (size_t)np);
    fi_param_combos(rl, V, nv, origExpr, ps, np, gsz, chosen, 0, ws, want);
    free(chosen);
    if (ws->n < want && np == 1)
        fi_solve_param(rl, V, nv, origExpr, Rconstr, ps[0], ws, want);
    for (int i = 0; i < np; i++) expr_free(ps[i]);
    free(ps);
}

/* Split expr's top-level conjunction into the Equal conjuncts and the rest. */
static Expr* fi_conj(Expr** parts, int n) {
    if (n == 0) { free(parts); return NULL; }
    if (n == 1) { Expr* r = parts[0]; free(parts); return r; }
    Expr* r = expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)n);
    free(parts);
    return r;
}
static void fi_split_eqs(const Expr* expr, Expr** Eout, Expr** Rout) {
    int cap = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    Expr** es = malloc(sizeof(Expr*) * (size_t)cap);
    Expr** rs = malloc(sizeof(Expr*) * (size_t)cap);
    int ne = 0, nr = 0;
    if (is_head_sym(expr, SYM_And)) {
        for (size_t i = 0; i < nargs(expr); i++) {
            const Expr* c = argn(expr, i);
            if (is_head_sym(c, SYM_Equal)) es[ne++] = xcopy(c); else rs[nr++] = xcopy(c);
        }
    } else if (is_head_sym(expr, SYM_Equal)) es[ne++] = xcopy(expr);
    else rs[nr++] = xcopy(expr);
    *Eout = fi_conj(es, ne);
    *Rout = fi_conj(rs, nr);
}

/* ---- bounded integer search (Diophantine witness / finite-domain emptiness) --- *
 *                                                                                *
 *  When Reduce and Solve decline over the Integers, enumerate an integer box.    *
 *  Bounds come from the constraints (fi_scan_bounds).  A finite box is decidable *
 *  (witness, or {} on exhaustion); an open box is best-effort within a budget    *
 *  (witness, or decline -- never {}).  A linear-equality reach-range check gives *
 *  an instant {} for the common knapsack shape.  Every point is verified.        */

typedef struct { long lo, hi; bool has_lo, has_hi; } FiRange;

static bool fi_is_number(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_REAL
                 || e->type == EXPR_BIGINT || is_head_sym(e, SYM_Rational));
}
/* Tightest integer satisfying `strict/non-strict` bound b on the given side. */
static bool fi_tight_int(const Expr* b, bool strict, bool lower, long* out) {
    bool ok = false; double d = rru_approx_double(b, &ok); if (!ok) return false;
    long g = (long)d;
    for (int it = 0; it < 200; it++) {
        Expr* ge = expr_new_integer(g);
        int sc = rru_sign_compare(ge, b); expr_free(ge);   /* sign of (g - b) */
        bool sat = lower ? (strict ? sc > 0 : sc >= 0) : (strict ? sc < 0 : sc <= 0);
        if (lower) {
            if (!sat) { g++; continue; }
            Expr* g1 = expr_new_integer(g - 1); int sc1 = rru_sign_compare(g1, b); expr_free(g1);
            if (strict ? sc1 > 0 : sc1 >= 0) { g--; continue; }
        } else {
            if (!sat) { g--; continue; }
            Expr* g1 = expr_new_integer(g + 1); int sc1 = rru_sign_compare(g1, b); expr_free(g1);
            if (strict ? sc1 < 0 : sc1 <= 0) { g++; continue; }
        }
        *out = g; return true;
    }
    *out = g; return true;
}
/* base^exp evaluated (exp a small nonneg integer). */
static Expr* fi_pow_eval(Expr* base, long exp) {
    Expr* e[2] = { base, expr_new_integer(exp) };
    return fi_eval_take(expr_new_function(expr_new_symbol(SYM_Power), e, 2));
}

/* ---- sum-of-even-powers bound (a positive-definite diagonal equality) ------- *
 *                                                                                *
 *  An equality that rearranges to  c1 v1^(2k1) + ... + cm vm^(2km) == B  with     *
 *  every c_i > 0 and B >= 0 bounds each of its variables: since every term is     *
 *  nonnegative, c_j v_j^(2k_j) <= B, so |v_j| <= (B/c_j)^(1/(2k_j)).  This is a    *
 *  strictly necessary condition, hence a SOUND box for the integer search --      *
 *  the a^2+b^2+c^2+d^2+e^2==5 conjunct pins every variable to [-2,2].  Cross      *
 *  terms (a b, ...) and odd powers make an equality non-conforming; it is then    *
 *  skipped (other conjuncts may still bound).                                     */

/* Largest nonnegative integer m with coeff * m^(2k) <= B (all exact).  Returns
 * false if no finite bound is found within a sane cap. */
static bool fi_sos_solve_m(const Expr* coeff, long two_k, const Expr* B, long* out_m) {
    /* estimate via doubles, then adjust exactly so rounding never under-counts */
    bool okc = false, okb = false;
    double dc = rru_approx_double(coeff, &okc);
    double db = rru_approx_double(B, &okb);
    if (!okc || !okb || dc <= 0.0 || db < 0.0) return false;
    long m = (long)floor(pow(db / dc, 1.0 / (double)two_k)) + 2;
    if (m < 0) m = 0;
    if (m > 200000) return false;                 /* too large -> leave unbounded */
    /* shrink while coeff*m^(2k) > B */
    for (; m > 0; m--) {
        Expr* t = fi_mul_eval(xcopy(coeff), fi_pow_eval(expr_new_integer(m), two_k));
        int sc = rru_sign_compare(t, B); expr_free(t);      /* sign of (t - B) */
        if (sc <= 0) break;
    }
    /* grow while coeff*(m+1)^(2k) <= B (repairs an under-estimate) */
    for (;;) {
        Expr* t = fi_mul_eval(xcopy(coeff), fi_pow_eval(expr_new_integer(m + 1), two_k));
        int sc = rru_sign_compare(t, B); expr_free(t);
        if (sc > 0) break;
        m++;
        if (m > 200000) return false;
    }
    *out_m = m;
    return true;
}

/* If the top-level statement is an And of atoms (or a single atom) with no Or,
 * derive a symmetric bound |v| <= m from any conforming even-power equality.
 * Returns true and sets *out_m when a bound is found. */
static bool fi_sos_upper_bound(const Expr* expr, const Expr* v, Expr** V, int nv, long* out_m) {
    if (is_head_sym(expr, SYM_Or)) return false;      /* per-clause bound only */
    int nc = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    bool found = false; long best = 0;
    for (int c = 0; c < nc; c++) {
        const Expr* conj = is_head_sym(expr, SYM_And) ? argn(expr, c) : expr;
        if (!(is_head_sym(conj, SYM_Equal) && nargs(conj) == 2)) continue;
        /* poly == 0  where poly = L - R */
        Expr* poly = fi_sub_eval(xcopy(argn(conj, 0)), xcopy(argn(conj, 1)));
        int nt = is_head_sym(poly, SYM_Plus) ? (int)nargs(poly) : 1;
        Expr* offset = expr_new_integer(0);           /* sum of constant terms  */
        Expr* coeff_v = NULL; long exp_v = 0;         /* this var's term, if any */
        bool conforming = true;
        for (int i = 0; i < nt && conforming; i++) {
            const Expr* term = is_head_sym(poly, SYM_Plus) ? argn(poly, i) : poly;
            if (fi_is_number(term)) { offset = fi_add_eval(offset, xcopy(term)); continue; }
            /* term must be  coeff * base^(2k)  with base a single listed var */
            Expr* coeff = expr_new_integer(1);
            const Expr* powpart = NULL;
            if (is_head_sym(term, SYM_Power)) {
                powpart = term;
            } else if (is_head_sym(term, SYM_Times)) {
                for (size_t j = 0; j < nargs(term); j++) {
                    const Expr* f = argn(term, j);
                    if (fi_is_number(f)) coeff = fi_mul_eval(coeff, xcopy(f));
                    else if (!powpart && is_head_sym(f, SYM_Power)) powpart = f;
                    else { conforming = false; break; }
                }
            } else { conforming = false; }
            long ex = 0; int vi = -1;
            if (conforming && powpart && nargs(powpart) == 2
                && argn(powpart, 1)->type == EXPR_INTEGER) {
                ex = argn(powpart, 1)->data.integer;
                vi = fi_var_index(argn(powpart, 0), V, nv);
            }
            /* need: a listed var, an even positive exponent, a positive coeff */
            Expr* zero = expr_new_integer(0);
            bool cpos = rru_sign_compare(coeff, zero) > 0; expr_free(zero);
            if (!conforming || !powpart || vi < 0 || ex < 2 || (ex & 1) || !cpos) {
                conforming = false; expr_free(coeff);
            } else {
                if (fi_is_var(argn(powpart, 0), v)) {   /* it is OUR variable */
                    if (coeff_v) expr_free(coeff);       /* v twice -> keep first */
                    else { coeff_v = coeff; exp_v = ex; }
                } else {
                    expr_free(coeff);                    /* another var's term */
                }
            }
        }
        if (conforming && coeff_v) {
            Expr* B = fi_neg(xcopy(offset)); B = fi_eval_take(B);  /* B = -offset */
            Expr* zero = expr_new_integer(0);
            long m;
            if (rru_sign_compare(B, zero) >= 0 && fi_sos_solve_m(coeff_v, exp_v, B, &m)) {
                if (!found || m < best) { best = m; found = true; }
            }
            expr_free(zero); expr_free(B);
        }
        if (coeff_v) expr_free(coeff_v);
        expr_free(offset); expr_free(poly);
    }
    if (found) *out_m = best;
    return found;
}

static void fi_int_range(const Expr* expr, const Expr* v, Expr** V, int nv, FiRange* out) {
    Expr* lo = NULL, *hi = NULL; bool los = false, his = false;
    Expr** forb = NULL; int nf = 0, cf = 0;
    fi_scan_bounds(expr, v, V, nv, &lo, &los, &hi, &his, &forb, &nf, &cf);
    out->has_lo = out->has_hi = false; out->lo = out->hi = 0;
    if (lo) out->has_lo = fi_tight_int(lo, los, true, &out->lo);
    if (hi) out->has_hi = fi_tight_int(hi, his, false, &out->hi);
    free(forb);
    /* Intersect a sum-of-even-powers bound |v| <= m (a sound necessary box). */
    long m;
    if (fi_sos_upper_bound(expr, v, V, nv, &m)) {
        if (!out->has_lo || out->lo < -m) out->lo = -m;
        if (!out->has_hi || out->hi >  m) out->hi =  m;
        out->has_lo = out->has_hi = true;
    }
}
/* Substitute the single variable var -> val into e and evaluate (owned). */
static Expr* fi_subst_one(const Expr* e, const Expr* var, long val) {
    Expr* r1[1] = { fi_rule(var, expr_new_integer(val)) };
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), r1, 1);
    Expr* out = fi_replace(e, rules);
    expr_free(rules);
    return out;
}
/* If a top-level linear equality's reachable range over the finite box excludes
 * its target, the whole system is unsatisfiable -> proven empty. */
static bool fi_linear_reject(const Expr* expr, Expr** V, int nv, FiRange* rg) {
    int nc = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    for (int c = 0; c < nc; c++) {
        const Expr* conj = is_head_sym(expr, SYM_And) ? argn(expr, c) : expr;
        if (!(is_head_sym(conj, SYM_Equal) && nargs(conj) == 2)) continue;
        Expr* poly = fi_sub_eval(xcopy(argn(conj, 0)), xcopy(argn(conj, 1)));
        /* constant term = poly with every var -> 0 */
        Expr** zr = malloc(sizeof(Expr*) * (size_t)nv);
        for (int i = 0; i < nv; i++) zr[i] = fi_rule(V[i], expr_new_integer(0));
        Expr* zrules = expr_new_function(expr_new_symbol(SYM_List), zr, (size_t)nv); free(zr);
        Expr* c0 = fi_replace(poly, zrules); expr_free(zrules);
        bool ok = fi_is_number(c0);
        Expr* mn = ok ? xcopy(c0) : NULL, *mx = ok ? xcopy(c0) : NULL;
        Expr* recon = ok ? xcopy(c0) : NULL;   /* c0 + sum a_i V_i, to confirm linearity */
        for (int i = 0; i < nv && ok; i++) {
            Expr* p1 = fi_subst_one(poly, V[i], 1);
            Expr* p0 = fi_subst_one(poly, V[i], 0);
            Expr* ai = fi_sub_eval(p1, p0);           /* coefficient of V[i] if linear */
            if (!fi_is_number(ai)) { expr_free(ai); ok = false; break; }
            recon = fi_add_eval(recon, fi_mul_eval(xcopy(ai), xcopy(V[i])));
            Expr* zero = expr_new_integer(0);
            int za = rru_sign_compare(ai, zero); expr_free(zero);
            if (za != 0 && !(rg[i].has_lo && rg[i].has_hi)) { expr_free(ai); ok = false; break; }
            Expr* tlo = fi_mul_eval(xcopy(ai), expr_new_integer(rg[i].lo));
            Expr* thi = fi_mul_eval(ai, expr_new_integer(rg[i].hi));
            int sc = rru_sign_compare(tlo, thi);
            Expr* tmin = sc <= 0 ? tlo : thi, *tmax = sc <= 0 ? thi : tlo;
            mn = fi_add_eval(mn, xcopy(tmin));
            mx = fi_add_eval(mx, xcopy(tmax));
            expr_free(tlo); expr_free(thi);
        }
        bool reject = false;
        if (ok) {
            Expr* resid = fi_sub_eval(xcopy(poly), recon); recon = NULL;
            bool linear = fi_is_number(resid) && resid->type == EXPR_INTEGER
                          && resid->data.integer == 0;
            expr_free(resid);
            if (linear) {
                Expr* z1 = expr_new_integer(0), *z2 = expr_new_integer(0);
                if (rru_sign_compare(mn, z1) > 0 || rru_sign_compare(mx, z2) < 0) reject = true;
                expr_free(z1); expr_free(z2);
            }
        }
        if (recon) expr_free(recon);
        if (mn) expr_free(mn);
        if (mx) expr_free(mx);
        expr_free(c0); expr_free(poly);
        if (reject) return true;
    }
    return false;
}
/* Returns 1 = added witness(es), 0 = proven empty, -1 = decline. */
static int fi_integer_search(const Expr* expr, Expr** V, int nv, long want, FiWit* ws) {
    if (nv <= 0 || nv > 24) return -1;
    FiRange* rg = malloc(sizeof(FiRange) * (size_t)nv);
    for (int i = 0; i < nv; i++) fi_int_range(expr, V[i], V, nv, &rg[i]);
    bool all_finite = true;
    for (int i = 0; i < nv; i++) if (!(rg[i].has_lo && rg[i].has_hi)) all_finite = false;
    if (all_finite && fi_linear_reject(expr, V, nv, rg)) { free(rg); return 0; }

    /* Candidate-evaluation budget: bounds the best-effort (open-domain) search and
     * the exhaustive proof of a finite box.  Small solutions surface early; an
     * out-of-reach system (e.g. a^4+b^4+c^4==d^4) declines once it is spent. */
    const long BUDGET = 120000;
    double finite_prod = 1.0; int n_open = 0;
    for (int i = 0; i < nv; i++) {
        if (rg[i].has_lo && rg[i].has_hi) {
            if (rg[i].hi < rg[i].lo) { free(rg); return 0; }   /* empty range */
            finite_prod *= (double)(rg[i].hi - rg[i].lo + 1);
        } else n_open++;
    }
    long span = 1;
    if (n_open > 0)
        while (span < 200000) {
            double t = finite_prod; bool of = false;
            for (int k = 0; k < n_open; k++) { t *= (double)(span + 1); if (t > (double)BUDGET) { of = true; break; } }
            if (of) break;
            span++;
        }
    long* elo = malloc(sizeof(long) * (size_t)nv), *ehi = malloc(sizeof(long) * (size_t)nv);
    for (int i = 0; i < nv; i++) {
        if (rg[i].has_lo && rg[i].has_hi) { elo[i] = rg[i].lo; ehi[i] = rg[i].hi; }
        else if (rg[i].has_lo)            { elo[i] = rg[i].lo; ehi[i] = rg[i].lo + span; }
        else if (rg[i].has_hi)            { ehi[i] = rg[i].hi; elo[i] = rg[i].hi - span; }
        else                              { elo[i] = -span; ehi[i] = span; }
    }
    bool can_prove_empty = all_finite && (finite_prod <= (double)BUDGET);
    long* cur = malloc(sizeof(long) * (size_t)nv);
    for (int i = 0; i < nv; i++) cur[i] = elo[i];
    long evals = 0; bool budget_hit = false, exhausted = false; int start_n = ws->n;
    for (;;) {
        if (evals >= BUDGET) { budget_hit = true; break; }
        evals++;
        Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
        for (int i = 0; i < nv; i++) outr[i] = fi_rule(V[i], expr_new_integer(cur[i]));
        Expr* p = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
        free(outr);
        if (fi_verify(expr, p)) { fi_wit_add(ws, p); if (ws->n >= want) break; }
        else expr_free(p);
        int i = 0;
        for (; i < nv; i++) { cur[i]++; if (cur[i] <= ehi[i]) break; cur[i] = elo[i]; }
        if (i == nv) { exhausted = true; break; }
    }
    free(cur); free(elo); free(ehi); free(rg);
    if (ws->n > start_n) return 1;
    if (exhausted && can_prove_empty && !budget_hit) return 0;
    return -1;
}

/* ---- numerical witness (transcendental / inexact Real systems) -------------- *
 *                                                                                *
 *  Reduce is not a sound decision procedure over transcendental functions or     *
 *  inexact numbers -- e.g. it wrongly returns False for                          *
 *  `0<x<0.001 && Sin[1/x]>0.999`.  For such systems we neither trust its False    *
 *  nor read its cells; we pose the statement as a numerical feasibility problem   *
 *  `NMinimize[{0, expr}, vars]` and VERIFY the returned point.                    */

/* Does e contain an inexact real or a transcendental function head?  Such systems
 * are outside Reduce's sound algebraic reach; a numerical witness is appropriate. */
static bool fi_is_transc_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Sin || n == SYM_Cos || n == SYM_Tan || n == SYM_Cot
            || n == SYM_Sec || n == SYM_Csc || n == SYM_Exp || n == SYM_Log
            || n == SYM_ArcSin || n == SYM_ArcCos || n == SYM_ArcTan || n == SYM_ArcTanh
            || n == SYM_Sinh || n == SYM_Cosh || n == SYM_Tanh
            || n == SYM_Gamma || n == SYM_Erf
            || n == SYM_PolyGamma || n == SYM_LogGamma || n == SYM_Zeta)
            return true;
        /* Exp[x] canonicalises to Power[E, x]; treat a base-E power as
         * exponential (transcendental) so Reduce's False is not trusted. */
        if (n == SYM_Power && e->data.function.arg_count == 2
            && fi_is_sym(e->data.function.args[0], SYM_E))
            return true;
    }
    if (fi_is_transc_inexact(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fi_is_transc_inexact(e->data.function.args[i])) return true;
    return false;
}

/* NMinimize[{0, expr}, vars] -- a feasibility solve; verify the returned point.
 * Returns an owned List[Rule..] witness, or NULL. */
static Expr* fi_numeric_search(const Expr* expr, const Expr* vars) {
    Expr* obj[2] = { expr_new_integer(0), xcopy(expr) };
    Expr* objl = expr_new_function(expr_new_symbol(SYM_List), obj, 2);
    Expr* nm[2] = { objl, xcopy(vars) };
    Expr* res = fi_eval_take(expr_new_function(expr_new_symbol(SYM_NMinimize), nm, 2));
    Expr* witness = NULL;
    if (is_head_sym(res, SYM_List) && nargs(res) == 2) {
        const Expr* rl = argn(res, 1);   /* {var -> val, ...} */
        if (is_head_sym(rl, SYM_List) && fi_verify(expr, rl)) witness = xcopy(rl);
    }
    expr_free(res);
    return witness;
}

/* ---- Groebner emptiness (declined polynomial systems) ----------------------- *
 *                                                                                *
 *  A positive-dimensional polynomial system with a disequation (e.g. the 2x2     *
 *  nilpotent M^2==0 && det!=0, which is empty since nilpotent => det==0) is       *
 *  declined by Solve (`nsdim`).  Rabinowitsch: the equalities together with       *
 *  t*(prod of disequation LHSs) - 1 have Groebner basis {1} exactly when the      *
 *  variety (minus the disequation hypersurfaces) is empty over C -- hence over    *
 *  R and Z.  A one-way, sound emptiness certificate.                             */
static bool fi_groebner_empty(const Expr* expr, const Expr* vars, int nv) {
    if (nv < 1 || nv > 8) return false;
    if (fi_is_transc_inexact(expr)) return false;   /* algebraic systems only */
    int nc = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    Expr** polys = malloc(sizeof(Expr*) * (size_t)(nc + 1));
    Expr** dis   = malloc(sizeof(Expr*) * (size_t)(nc + 1));
    int npoly = 0, ndis = 0;
    bool ok = true;
    for (int c = 0; c < nc && ok; c++) {
        const Expr* a = is_head_sym(expr, SYM_And) ? argn(expr, c) : expr;
        if (is_head_sym(a, SYM_Equal) && nargs(a) == 2)
            polys[npoly++] = fi_sub_eval(xcopy(argn(a, 0)), xcopy(argn(a, 1)));
        else if (is_head_sym(a, SYM_Unequal) && nargs(a) == 2)
            dis[ndis++] = fi_sub_eval(xcopy(argn(a, 0)), xcopy(argn(a, 1)));
        else if (is_head_sym(a, SYM_Less) || is_head_sym(a, SYM_LessEqual)
                 || is_head_sym(a, SYM_Greater) || is_head_sym(a, SYM_GreaterEqual)
                 || is_head_sym(a, SYM_Inequality))
            continue;                                /* inequalities only shrink the set */
        else ok = false;                             /* unrecognised atom -> give up */
    }
    bool empty = false;
    if (ok && npoly >= 1) {
        int ng = npoly + (ndis > 0 ? 1 : 0);
        Expr** gens = malloc(sizeof(Expr*) * (size_t)ng);
        for (int i = 0; i < npoly; i++) gens[i] = polys[i];   /* transfer ownership */
        Expr* tsym = NULL;
        if (ndis > 0) {
            Expr* prod = dis[0];
            for (int i = 1; i < ndis; i++) prod = fi_mul_eval(prod, dis[i]);
            tsym = expr_new_symbol("FindInstance$t");
            Expr* rab = fi_sub_eval(fi_mul_eval(xcopy(tsym), prod), expr_new_integer(1));
            gens[npoly] = rab;
        }
        Expr* genlist = expr_new_function(expr_new_symbol(SYM_List), gens, (size_t)ng);
        free(gens);
        /* variable list: the listed vars, plus t when Rabinowitsch is used */
        int nvl = nv + (tsym ? 1 : 0);
        Expr** vl = malloc(sizeof(Expr*) * (size_t)nvl);
        if (is_head_sym(vars, SYM_List))
            for (int i = 0; i < nv; i++) vl[i] = xcopy(argn(vars, i));
        else vl[0] = xcopy(vars);
        if (tsym) vl[nv] = xcopy(tsym);
        Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vl, (size_t)nvl);
        free(vl);
        if (tsym) expr_free(tsym);
        Expr* gb[2] = { genlist, varlist };
        Expr* res = fi_eval_take(expr_new_function(expr_new_symbol(SYM_GroebnerBasis), gb, 2));
        if (is_head_sym(res, SYM_List) && nargs(res) == 1) {
            const Expr* g0 = argn(res, 0);
            if (g0->type == EXPR_INTEGER && g0->data.integer != 0) empty = true;   /* {c} = {1} */
        }
        expr_free(res);
    } else {
        for (int i = 0; i < npoly; i++) expr_free(polys[i]);
        for (int i = 0; i < ndis; i++) expr_free(dis[i]);
    }
    free(polys); free(dis);
    return empty;
}

/* ---- structured candidate sampling (ℂ / ℝ disequations & inequalities) ------ *
 *                                                                                *
 *  Reduce/Solve decline many statements that a concrete point trivially          *
 *  satisfies -- branch-cut disequations (Sqrt[z^2] != z at z=-1;                  *
 *  Log[x y] != Log[x]+Log[y] at x=y=-1), open regions, and the like.  As a LAST  *
 *  resort over Complexes/Reals we evaluate the statement at a small, ordered      *
 *  grid of interesting values and keep any point that VERIFIES.  Purely additive  *
 *  (it never claims {}), so it cannot corrupt a sound empty/decline; and the      *
 *  fi_verify gate keeps every returned point correct.                            */

/* Build the ordered candidate value list (owned).  -1 and the small magnitudes
 * come first so the common witnesses (e.g. all-vars = -1) are hit immediately. */
static Expr** fi_make_candidates(bool allow_complex, int* n_out) {
    /* integer / rational reals, small first */
    static const long ints[] = { -1, 1, 2, -2, 0, 3, -3, 5, -5 };
    struct { long num, den; } rats[] = { {1,2}, {-1,2}, {3,2} };
    int cap = 32;
    Expr** c = malloc(sizeof(Expr*) * (size_t)cap);
    int n = 0;
    for (size_t i = 0; i < sizeof ints / sizeof ints[0]; i++)
        c[n++] = expr_new_integer(ints[i]);
    for (size_t i = 0; i < sizeof rats / sizeof rats[0]; i++)
        c[n++] = fi_eval_take(fi_bin(SYM_Times, expr_new_integer(rats[i].num),
                                     fi_pow_eval(expr_new_integer(rats[i].den), -1)));
    if (allow_complex) {
        /* I, -I, 2 I, and the four unit diagonals */
        long re[] = { 0,  0,  0,  1,  1, -1, -1 };
        long im[] = { 1, -1,  2,  1, -1,  1, -1 };
        for (size_t i = 0; i < sizeof re / sizeof re[0]; i++) {
            Expr* bI = fi_bin(SYM_Times, expr_new_integer(im[i]), expr_new_symbol(SYM_I));
            c[n++] = fi_eval_take(fi_bin(SYM_Plus, expr_new_integer(re[i]), bI));
        }
    }
    *n_out = n;
    return c;
}

/* Odometer over the candidate grid, capped at BUDGET evaluations.  Adds every
 * point that verifies until `want` are collected. */
static void fi_sample_search(const Expr* expr, Expr** V, int nv, bool allow_complex,
                             long want, FiWit* ws) {
    if (nv < 1 || nv > 6) return;
    int nc = 0;
    Expr** cand = fi_make_candidates(allow_complex, &nc);
    const long BUDGET = 8000;
    int* idx = calloc((size_t)nv, sizeof(int));
    for (long step = 0; step < BUDGET && ws->n < want; step++) {
        Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
        for (int i = 0; i < nv; i++) outr[i] = fi_rule(V[i], expr_copy(cand[idx[i]]));
        Expr* pt = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
        free(outr);
        if (fi_verify(expr, pt)) fi_wit_add(ws, pt); else expr_free(pt);
        int i = 0;                                    /* advance the odometer */
        for (; i < nv; i++) { if (++idx[i] < nc) break; idx[i] = 0; }
        if (i == nv) break;                           /* full grid exhausted */
    }
    free(idx);
    for (int i = 0; i < nc; i++) expr_free(cand[i]);
    free(cand);
}

/* ---- 1-variable real transcendental root search ---------------------------- *
 *                                                                                *
 *  A single real variable and one equation lhs == rhs (plus any bounding          *
 *  inequalities) that Reduce/Solve decline -- Tan[x]==x && x>10^6.  We scan a      *
 *  POLE-FREE refactoring of lhs-rhs for a sign change, refine the bracket by       *
 *  high-precision bisection, and VERIFY against the original statement.           *
 *                                                                                 *
 *  The refactoring multiplies lhs-rhs by the denominators of any tangent/         *
 *  secant/... so the whole oscillating-but-smooth numerator is bracketable:        *
 *  Tan[x]-x has a pole ~10^-6 from its root, but (Tan[x]-x)Cos[x] = Sin[x]-x Cos[x] *
 *  does not, so a coarse grid can straddle it.  Machine trig at 10^6 is            *
 *  meaningless (argument reduction), so the refinement runs at high precision;     *
 *  the coarse scan only needs the sign, which stays reliable while |h| is large.  */

/* Denominator head that clears the pole of a tangent-like head, else NULL. */
static const char* fi_pole_factor_head(const char* h) {
    if (h == SYM_Tan  || h == SYM_Sec)  return SYM_Cos;
    if (h == SYM_Cot  || h == SYM_Csc)  return SYM_Sin;
    if (h == SYM_Tanh || h == SYM_Sech) return SYM_Cosh;
    if (h == SYM_Coth || h == SYM_Csch) return SYM_Sinh;
    return NULL;
}
static bool fi_has_transc_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (n == SYM_Sin || n == SYM_Cos || n == SYM_Tan || n == SYM_Cot
            || n == SYM_Sec || n == SYM_Csc || n == SYM_Exp || n == SYM_Log
            || n == SYM_Sinh || n == SYM_Cosh || n == SYM_Tanh || n == SYM_Coth
            || n == SYM_Sech || n == SYM_Csch || n == SYM_ArcTan || n == SYM_ArcSin
            || n == SYM_ArcCos || n == SYM_ArcTanh)
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fi_has_transc_head(e->data.function.args[i])) return true;
    return false;
}
static void fi_collect_poles(const Expr* node, Expr*** fs, int* nf, int* cf) {
    if (!node || node->type != EXPR_FUNCTION) return;
    const Expr* h = node->data.function.head;
    if (h->type == EXPR_SYMBOL && node->data.function.arg_count == 1) {
        const char* fh = fi_pole_factor_head(h->data.symbol.name);
        if (fh) {
            Expr* a[1] = { xcopy(node->data.function.args[0]) };
            Expr* fac = expr_new_function(expr_new_symbol(fh), a, 1);
            bool dup = false;
            for (int i = 0; i < *nf; i++) if (expr_eq((*fs)[i], fac)) { dup = true; break; }
            if (dup) expr_free(fac);
            else {
                if (*nf == *cf) { *cf = *cf ? *cf * 2 : 4; *fs = realloc(*fs, sizeof(Expr*) * (size_t)*cf); }
                (*fs)[(*nf)++] = fac;
            }
        }
    }
    for (size_t i = 0; i < node->data.function.arg_count; i++)
        fi_collect_poles(node->data.function.args[i], fs, nf, cf);
}
/* h = Expand[(lhs - rhs) * (product of pole denominators)] -- pole-free. */
static Expr* fi_smooth_residual(const Expr* L, const Expr* R) {
    Expr* h = fi_sub_eval(xcopy(L), xcopy(R));
    Expr** fs = NULL; int nf = 0, cf = 0;
    fi_collect_poles(h, &fs, &nf, &cf);
    if (nf == 0) { free(fs); return h; }
    Expr* prod = fs[0];
    for (int i = 1; i < nf; i++) prod = fi_mul_eval(prod, fs[i]);
    free(fs);
    Expr* hp = fi_mul_eval(h, prod);
    Expr* a[1] = { hp };
    return fi_eval_take(expr_new_function(expr_new_symbol(SYM_Expand), a, 1));
}
/* Machine-precision value of hs at v = x (finite iff a real machine number). */
static double fi_eval_real_at(const Expr* hs, const Expr* v, double x, bool* finite) {
    Expr* r1[1] = { fi_rule(v, expr_new_real(x)) };
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), r1, 1);
    Expr* sub = fi_replace(hs, rules); expr_free(rules);
    bool ok = false; double d = rru_approx_double(sub, &ok);
    expr_free(sub);
    *finite = ok && isfinite(d);
    return d;
}
/* Sign of hs at the exact rational point, evaluated to `prec` digits. */
static int fi_hp_sign_at(const Expr* hs, const Expr* v, const Expr* rat, long prec) {
    Expr* r1[1] = { fi_rule(v, xcopy(rat)) };
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), r1, 1);
    Expr* sub = fi_replace(hs, rules); expr_free(rules);
    Expr* na[2] = { sub, expr_new_integer(prec) };
    Expr* nx = fi_eval_take(expr_new_function(expr_new_symbol(SYM_N), na, 2));
    Expr* sa[1] = { nx };
    Expr* sg = fi_eval_take(expr_new_function(expr_new_symbol(SYM_Sign), sa, 1));
    int r = (sg && sg->type == EXPR_INTEGER) ? (int)sg->data.integer : 0;
    if (sg) expr_free(sg);
    return r;
}
static Expr* fi_rationalize(double x) {
    Expr* a[2] = { expr_new_real(x), expr_new_integer(0) };
    return fi_eval_take(expr_new_function(expr_new_symbol("Rationalize"), a, 2));
}
static Expr* fi_rat_mid(const Expr* a, const Expr* b) {
    return fi_mul_eval(fi_add_eval(xcopy(a), xcopy(b)),
                       fi_pow_eval(expr_new_integer(2), -1));
}
/* Refine a machine bracket [xa,xb] to a high-precision root, verify, add. */
static void fi_refine_and_verify(const Expr* hs, const Expr* origExpr, const Expr* v,
                                 double xa, double xb, long prec, FiWit* ws) {
    Expr* a = fi_rationalize(xa);
    Expr* b = fi_rationalize(xb);
    int sa = fi_hp_sign_at(hs, v, a, prec);
    int sb = fi_hp_sign_at(hs, v, b, prec);
    if (sa != 0 && sb != 0 && sa != sb) {
        for (int it = 0; it < 170; it++) {
            Expr* m = fi_rat_mid(a, b);
            int sm = fi_hp_sign_at(hs, v, m, prec);
            if (sm == 0) { expr_free(a); expr_free(b); a = m; b = xcopy(m); break; }
            if (sm == sa) { expr_free(a); a = m; } else { expr_free(b); b = m; }
        }
        Expr* mid = fi_rat_mid(a, b);
        Expr* na[2] = { mid, expr_new_integer(prec) };
        Expr* xstar = fi_eval_take(expr_new_function(expr_new_symbol(SYM_N), na, 2));
        Expr* rr[1] = { fi_rule(v, xstar) };
        Expr* pt = expr_new_function(expr_new_symbol(SYM_List), rr, 1);
        if (fi_verify(origExpr, pt)) fi_wit_add(ws, pt); else expr_free(pt);
    }
    expr_free(a); expr_free(b);
}
static void fi_scan_dir(const Expr* hs, const Expr* origExpr, const Expr* v,
                        double start, int dir, bool has_far, double far,
                        double step, long maxsteps, long prec, long want, FiWit* ws) {
    bool fprev; double hprev = fi_eval_real_at(hs, v, start, &fprev);
    for (long k = 1; k <= maxsteps && ws->n < want; k++) {
        double xcur = start + (double)dir * step * (double)k;
        if (has_far && ((dir > 0 && xcur > far) || (dir < 0 && xcur < far))) break;
        bool fcur; double hcur = fi_eval_real_at(hs, v, xcur, &fcur);
        if (fprev && fcur && hprev != 0.0 && hcur != 0.0
            && ((hprev < 0.0) != (hcur < 0.0)))
            fi_refine_and_verify(hs, origExpr, v, start + (double)dir * step * (double)(k - 1),
                                 xcur, prec, ws);
        hprev = hcur; fprev = fcur;
    }
}
static void fi_real_root_search(const Expr* expr, Expr** V, int nv, const Expr* dom,
                                long want, FiWit* ws) {
    if (nv != 1 || !fi_is_sym(dom, SYM_Reals) || !fi_has_transc_head(expr)) return;
    Expr* E = NULL, *Rp = NULL;
    fi_split_eqs(expr, &E, &Rp);
    if (!E || !(is_head_sym(E, SYM_Equal) && nargs(E) == 2)) {
        if (E) expr_free(E);
        if (Rp) expr_free(Rp);
        return;
    }
    Expr* hs = fi_smooth_residual(argn(E, 0), argn(E, 1));
    Expr* lo = NULL, *hi = NULL; bool los = false, his = false;
    Expr** forb = NULL; int nf = 0, cf = 0;
    fi_scan_bounds(expr, V[0], V, 1, &lo, &los, &hi, &his, &forb, &nf, &cf);
    free(forb);
    bool okl = false, okh = false; double dl = 0, dh = 0;
    if (lo) dl = rru_approx_double(lo, &okl);
    if (hi) dh = rru_approx_double(hi, &okh);
    const double STEP = 0.5; const long MAXS = 4000; const long PREC = 50;
    if (okl && !okh)      fi_scan_dir(hs, expr, V[0], dl, +1, false, 0, STEP, MAXS, PREC, want, ws);
    else if (okh && !okl) fi_scan_dir(hs, expr, V[0], dh, -1, false, 0, STEP, MAXS, PREC, want, ws);
    else if (okl && okh)  fi_scan_dir(hs, expr, V[0], dl, +1, true, dh, STEP, MAXS, PREC, want, ws);
    else {
        fi_scan_dir(hs, expr, V[0], 0.0, +1, false, 0, STEP, MAXS, PREC, want, ws);
        if (ws->n < want) fi_scan_dir(hs, expr, V[0], 0.0, -1, false, 0, STEP, MAXS, PREC, want, ws);
    }
    expr_free(hs); expr_free(E); if (Rp) expr_free(Rp);
}

/* ---- numeric feasibility with tolerance verify (inexact Real systems) ------- *
 *                                                                                *
 *  A system that carries an inexact number (e.g. == 0.1) is not meant to be       *
 *  satisfied EXACTLY -- a machine-precision witness is the expected answer, so     *
 *  the exact `=== True` gate is wrong for it.  We minimise the sum of squared      *
 *  equation residuals with FindMinimum from several seeds (seeds spread the        *
 *  variables so paired unknowns like x and y land on distinct roots), then         *
 *  VERIFY the point NUMERICALLY: equalities to a tolerance, strict inequalities    *
 *  and disequations with a margin.  Reaches the damped-oscillation system          *
 *  Exp[-a x]Cos[b x]==0.1 && Exp[-a y]Cos[b y]==0.1 && x!=y && a>0.               */

#define FI_EQ_TOL   1e-6
#define FI_MARGIN   1e-9

/* Machine value of e at the point (owned point = List[Rule..]); finite iff real. */
static bool fi_point_val(const Expr* e, const Expr* point, double* out) {
    Expr* sub = fi_replace(e, point);
    bool ok = false; double d = rru_approx_double(sub, &ok);
    expr_free(sub);
    *out = d;
    return ok && isfinite(d);
}
/* Numeric three-valued truth of a (possibly compound) statement at a point:
 * 1 = holds within tolerance, 0 = fails, -1 = cannot decide numerically. */
static int fi_num_true(const Expr* node, const Expr* point);
static int fi_num_rel(const char* op, const Expr* L, const Expr* R, const Expr* point) {
    double a, b;
    if (!fi_point_val(L, point, &a) || !fi_point_val(R, point, &b)) return -1;
    double d = a - b;
    if (op == SYM_Equal)        return fabs(d) <= FI_EQ_TOL ? 1 : 0;
    if (op == SYM_Unequal)      return fabs(d) >  FI_EQ_TOL ? 1 : 0;
    if (op == SYM_Less)         return d < -FI_MARGIN ? 1 : (d >  FI_MARGIN ? 0 : -1);
    if (op == SYM_Greater)      return d >  FI_MARGIN ? 1 : (d < -FI_MARGIN ? 0 : -1);
    if (op == SYM_LessEqual)    return d <=  FI_EQ_TOL ? 1 : 0;
    if (op == SYM_GreaterEqual) return d >= -FI_EQ_TOL ? 1 : 0;
    return -1;
}
static int fi_num_true(const Expr* node, const Expr* point) {
    if (fi_is_sym(node, SYM_True))  return 1;
    if (fi_is_sym(node, SYM_False)) return 0;
    if (is_head_sym(node, SYM_And)) {
        int r = 1;
        for (size_t i = 0; i < nargs(node); i++) {
            int t = fi_num_true(argn(node, i), point);
            if (t == 0) return 0;
            if (t < 0) r = -1;
        }
        return r;
    }
    if (is_head_sym(node, SYM_Or)) {
        int r = 0;
        for (size_t i = 0; i < nargs(node); i++) {
            int t = fi_num_true(argn(node, i), point);
            if (t == 1) return 1;
            if (t < 0) r = -1;
        }
        return r;
    }
    if (is_head_sym(node, SYM_Not) && nargs(node) == 1) {
        int t = fi_num_true(argn(node, 0), point); return t < 0 ? -1 : !t;
    }
    if (node->type == EXPR_FUNCTION && node->data.function.head->type == EXPR_SYMBOL
        && nargs(node) == 2) {
        const char* h = node->data.function.head->data.symbol.name;
        if (h == SYM_Equal || h == SYM_Unequal || h == SYM_Less || h == SYM_Greater
            || h == SYM_LessEqual || h == SYM_GreaterEqual)
            return fi_num_rel(h, argn(node, 0), argn(node, 1), point);
    }
    if (is_head_sym(node, SYM_Inequality) && (nargs(node) & 1)) {
        int r = 1;
        for (size_t i = 0; i + 2 < nargs(node); i += 2) {
            const Expr* opE = argn(node, i + 1);
            if (opE->type != EXPR_SYMBOL) return -1;
            int t = fi_num_rel(opE->data.symbol.name, argn(node, i), argn(node, i + 2), point);
            if (t == 0) return 0;
            if (t < 0) r = -1;
        }
        return r;
    }
    return -1;
}

/* Penalty margin: how far strictly inside a strict inequality the optimiser is
 * pushed, so a `<` witness clears the numeric verify tolerances (FI_EQ_TOL,
 * FI_MARGIN) rather than sitting on the boundary. */
#define FI_PEN_MARGIN 1e-4

static void fi_push_term(Expr*** ts, int* nt, int* ct, Expr* term) {
    if (*nt == *ct) { *ct = *ct ? *ct * 2 : 4; *ts = realloc(*ts, sizeof(Expr*) * (size_t)*ct); }
    (*ts)[(*nt)++] = term;
}
/* Max[0, g] squared -- a one-sided hinge penalty, zero when g <= 0.  Takes g. */
static Expr* fi_hinge_sq(Expr* g) {
    Expr* mx[2] = { expr_new_integer(0), g };
    Expr* m = expr_new_function(expr_new_symbol(SYM_Max), mx, 2);
    return fi_pow_eval(m, 2);
}
/* (L - R) + margin, owned and evaluated. */
static Expr* fi_slack(const Expr* L, const Expr* R, double margin) {
    return fi_add_eval(fi_sub_eval(xcopy(L), xcopy(R)), expr_new_real(margin));
}
/* One relation's contribution to the least-infeasibility objective:
 * equality -> squared residual; strict/loose inequality -> one-sided hinge;
 * Unequal -> nothing (open, satisfied almost everywhere). */
static void fi_collect_rel_penalty(const char* op, const Expr* L, const Expr* R,
                                   Expr*** ts, int* nt, int* ct) {
    const double M = FI_PEN_MARGIN;
    if (op == SYM_Equal)
        fi_push_term(ts, nt, ct, fi_pow_eval(fi_sub_eval(xcopy(L), xcopy(R)), 2));
    else if (op == SYM_Less)             /* L < R : drive L-R below -M */
        fi_push_term(ts, nt, ct, fi_hinge_sq(fi_slack(L, R, M)));
    else if (op == SYM_Greater)          /* L > R : drive R-L below -M */
        fi_push_term(ts, nt, ct, fi_hinge_sq(fi_slack(R, L, M)));
    else if (op == SYM_LessEqual)        /* L <= R : boundary allowed */
        fi_push_term(ts, nt, ct, fi_hinge_sq(fi_sub_eval(xcopy(L), xcopy(R))));
    else if (op == SYM_GreaterEqual)     /* L >= R */
        fi_push_term(ts, nt, ct, fi_hinge_sq(fi_sub_eval(xcopy(R), xcopy(L))));
}
/* Least-infeasibility penalty of the whole statement (equalities AND
 * inequalities), collected for the FindMinimum objective.  A point with zero
 * penalty satisfies every equality and inequality -- so, unlike an
 * equality-only residual, the optimiser is kept inside the feasible box (e.g.
 * 0 < t < 0.1) and halts on a feasible open region (e.g. Rastrigin < 0.1)
 * rather than driving to an excluded minimiser. */
static void fi_collect_penalty(const Expr* node, Expr*** ts, int* nt, int* ct) {
    if (is_head_sym(node, SYM_And)) {
        for (size_t i = 0; i < nargs(node); i++) fi_collect_penalty(argn(node, i), ts, nt, ct);
        return;
    }
    if (node->type == EXPR_FUNCTION && node->data.function.head->type == EXPR_SYMBOL
        && nargs(node) == 2) {
        const char* h = node->data.function.head->data.symbol.name;
        if (h == SYM_Equal || h == SYM_Less || h == SYM_Greater
            || h == SYM_LessEqual || h == SYM_GreaterEqual) {
            fi_collect_rel_penalty(h, argn(node, 0), argn(node, 1), ts, nt, ct);
            return;
        }
    }
    if (is_head_sym(node, SYM_Inequality) && (nargs(node) & 1)) {
        for (size_t i = 0; i + 2 < nargs(node); i += 2) {
            const Expr* opE = argn(node, i + 1);
            if (opE->type == EXPR_SYMBOL)
                fi_collect_rel_penalty(opE->data.symbol.name, argn(node, i),
                                       argn(node, i + 2), ts, nt, ct);
        }
        return;
    }
    /* Unequal, Not, bare Boolean atoms: no penalty term. */
}

/* Fold definitional equalities `s == c` (s a symbol that is NOT a listed
 * variable, c a number) into substitutions and drop them, so a statement that
 * pins auxiliary symbols to constants -- d1 == 3.2 && d2 == 2.8 && g == 9.8 &&
 * ... -- reduces to one purely in the solve variables.  Returns an owned,
 * evaluated statement (a plain copy when nothing folds). */
static Expr* fi_fold_aux(const Expr* expr, Expr** V, int nv) {
    int nat = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    Expr** rules = malloc(sizeof(Expr*) * (size_t)nat);
    int nr = 0;
    for (int i = 0; i < nat; i++) {
        const Expr* a = is_head_sym(expr, SYM_And) ? argn(expr, i) : expr;
        if (!is_head_sym(a, SYM_Equal) || nargs(a) != 2) continue;
        const Expr* L = argn(a, 0); const Expr* R = argn(a, 1);
        const Expr* s = NULL, *c = NULL;
        if (L->type == EXPR_SYMBOL && fi_is_number(R)) { s = L; c = R; }
        else if (R->type == EXPR_SYMBOL && fi_is_number(L)) { s = R; c = L; }
        if (s && fi_var_index(s, V, nv) < 0) rules[nr++] = fi_rule(s, xcopy(c));
    }
    if (nr == 0) { free(rules); return xcopy(expr); }
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nr);
    free(rules);
    Expr* out = fi_replace(expr, rl);
    expr_free(rl);
    return out;
}

static Expr* fi_numeric_feasibility(const Expr* expr, Expr** V, int nv) {
    if (nv < 1 || nv > 8) return NULL;
    Expr* work = fi_fold_aux(expr, V, nv);   /* pin auxiliary symbols to constants */
    Expr** ts = NULL; int nt = 0, ct = 0;
    fi_collect_penalty(work, &ts, &nt, &ct);
    if (nt == 0) { free(ts); expr_free(work); return NULL; }
    Expr* obj = (nt == 1) ? ts[0]
              : expr_new_function(expr_new_symbol(SYM_Plus), ts, (size_t)nt);
    free(ts);
    /* Multi-start seeds base + i*spread: spread paired variables onto distinct
     * roots, keep positive-unknown seeds positive, and include small-magnitude
     * patterns so a feasible region hugging the origin (Rastrigin < 0.1) is
     * seeded from inside. */
    static const double base[]   = { 0.5, 1.0, 0.3, 2.0, 1.5, 0.7, 0.01, 0.05, 0.1 };
    static const double spread[] = { 1.0, 1.5, 2.0, 0.9, 0.6, 2.5, 0.003, 0.02, 0.05 };
    const int NP = (int)(sizeof base / sizeof base[0]);
    Expr* witness = NULL;
    for (int p = 0; p < NP && !witness; p++) {
        Expr** specs = malloc(sizeof(Expr*) * (size_t)nv);
        for (int i = 0; i < nv; i++) {
            Expr* s2[2] = { xcopy(V[i]), expr_new_real(base[p] + (double)i * spread[p]) };
            specs[i] = expr_new_function(expr_new_symbol(SYM_List), s2, 2);
        }
        Expr* varspec = expr_new_function(expr_new_symbol(SYM_List), specs, (size_t)nv);
        free(specs);
        Expr* fm2[2] = { xcopy(obj), varspec };
        Expr* fm = fi_eval_take(expr_new_function(expr_new_symbol(SYM_FindMinimum), fm2, 2));
        if (is_head_sym(fm, SYM_List) && nargs(fm) == 2) {
            bool okv = false; double mv = rru_approx_double(argn(fm, 0), &okv);
            const Expr* rl = argn(fm, 1);
            /* The FindMinimum value is only a screen; fi_num_true (verified against
             * the folded statement, at the real tolerances) is the soundness gate. */
            if (okv && mv < 1e-6 && is_head_sym(rl, SYM_List)
                && fi_num_true(work, rl) == 1)
                witness = xcopy(rl);
        }
        expr_free(fm);
    }
    expr_free(obj);
    expr_free(work);
    return witness;
}

/* ---- solve one equation for one variable, sample the rest ------------------ *
 *                                                                               *
 *  Reduce/Solve decline many Real/Complex systems as a WHOLE that become        *
 *  univariate once the other variables are pinned:                              *
 *    c1 e^{-L1 t} + c2 e^{-L2 t} == 0 && c1>0 && c2<0 && L1>L2>0 && t>0          *
 *      -- solve for c1 (or c2), sample the rest;                                 *
 *    (x^2-y^2)/(x^2+y^2) == 1/2 && x^2+y^2 < 10^-10 && x>0 && y>0                *
 *      -- solve for x, sample a tiny y.                                          *
 *  For each equation E and variable xk in it we Solve[E, xk], then odometer a    *
 *  constraint-aware candidate grid over the OTHER variables, compute xk, and     *
 *  VERIFY the full statement.  Budget-capped and verify-gated -- purely          *
 *  additive, so it can only add correct witnesses, never a wrong {}.            */

/* Is atom a univariate constraint in exactly V[i] (mentions V[i] and no other
 * listed variable)?  Such atoms pre-filter the per-variable candidate grid. */
static bool fi_atom_univar_in(const Expr* atom, Expr** V, int nv, int i) {
    if (!fi_contains_one(atom, V[i])) return false;
    for (int j = 0; j < nv; j++)
        if (j != i && fi_contains_one(atom, V[j])) return false;
    return true;
}
/* Candidate cv admissible for variable i: every univariate-in-i constraint atom
 * holds (or stays undecided) at V[i] = cv.  Coupled atoms defer to the verify. */
static bool fi_cand_ok_for_var(const Expr* cv, int i, Expr** V, int nv, const Expr* Rpart) {
    if (!Rpart) return true;
    int nat = is_head_sym(Rpart, SYM_And) ? (int)nargs(Rpart) : 1;
    Expr* onerule[1] = { fi_rule(V[i], xcopy(cv)) };
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), onerule, 1);
    bool ok = true;
    for (int k = 0; k < nat && ok; k++) {
        const Expr* a = is_head_sym(Rpart, SYM_And) ? argn(Rpart, k) : Rpart;
        if (!fi_atom_univar_in(a, V, nv, i)) continue;
        Expr* sub = fi_replace(a, rl);
        if (fi_is_sym(sub, SYM_False)) ok = false;
        expr_free(sub);
    }
    expr_free(rl);
    return ok;
}
/* fi_make_candidates plus small/large magnitudes, so tight coupled bounds
 * (x^2+y^2 < 10^-10) and larger scales are reachable when one var is solved. */
static Expr** fi_make_candidates_ext(bool allow_complex, int* n_out) {
    int nbase = 0;
    Expr** base = fi_make_candidates(allow_complex, &nbase);
    struct { long num, den; } extra[] =
        { {1,10}, {-1,10}, {1,1000}, {1,1000000}, {10,1}, {-10,1}, {7,1} };
    int ne = (int)(sizeof extra / sizeof extra[0]);
    Expr** c = malloc(sizeof(Expr*) * (size_t)(nbase + ne));
    int n = 0;
    for (int i = 0; i < nbase; i++) c[n++] = base[i];
    free(base);
    for (int i = 0; i < ne; i++)
        c[n++] = fi_eval_take(fi_bin(SYM_Times, expr_new_integer(extra[i].num),
                                     fi_pow_eval(expr_new_integer(extra[i].den), -1)));
    *n_out = n;
    return c;
}
/* Extract the value bound to V[k] in a Solve solution rule-list, or NULL. */
static Expr* fi_rl_value(const Expr* RL, const Expr* vk) {
    if (!is_head_sym(RL, SYM_List)) return NULL;
    for (size_t r = 0; r < nargs(RL); r++) {
        const Expr* rr = argn(RL, r);
        if (is_head_sym(rr, SYM_Rule) && nargs(rr) == 2 && fi_is_var(argn(rr, 0), vk))
            return xcopy(argn(rr, 1));
    }
    return NULL;
}
static void fi_solve_one_sample(const Expr* expr, Expr** V, int nv, const Expr* dom,
                                long want, FiWit* ws) {
    if (nv < 1 || nv > 6) return;
    Expr* Epart = NULL, *Rpart = NULL;
    fi_split_eqs(expr, &Epart, &Rpart);
    if (!Epart) { if (Rpart) expr_free(Rpart); return; }

    bool allow_complex = !fi_is_sym(dom, SYM_Reals) && !fi_is_sym(dom, SYM_Rationals)
                       && !fi_is_sym(dom, SYM_Integers);
    int nc = 0;
    Expr** cand = fi_make_candidates_ext(allow_complex, &nc);
    int neq = is_head_sym(Epart, SYM_And) ? (int)nargs(Epart) : 1;
    const long BUDGET = 6000;
    long spent = 0;

    for (int e = 0; e < neq && ws->n < want; e++) {
        const Expr* E = is_head_sym(Epart, SYM_And) ? argn(Epart, e) : Epart;
        if (!is_head_sym(E, SYM_Equal)) continue;
        for (int k = 0; k < nv && ws->n < want && spent < BUDGET; k++) {
            if (!fi_contains_one(E, V[k])) continue;
            Expr* sols = fi_call(SYM_Solve, E, V[k], NULL, NULL);
            if (!is_head_sym(sols, SYM_List)) { expr_free(sols); continue; }
            for (size_t si = 0; si < nargs(sols) && ws->n < want && spent < BUDGET; si++) {
                Expr* g = fi_rl_value(argn(sols, si), V[k]);
                if (!g) continue;
                int nf = nv - 1;
                int* fmap = malloc(sizeof(int) * (size_t)(nf > 0 ? nf : 1));
                { int m = 0; for (int j = 0; j < nv; j++) if (j != k) fmap[m++] = j; }
                /* Compact each free variable's admissible candidates up front, so
                 * a sign-constrained variable shrinks the odometer radix and its
                 * first feasible value sits at index 0 (keeping the winning small
                 * combination shallow within the shared budget). */
                Expr*** fc = malloc(sizeof(Expr**) * (size_t)(nf > 0 ? nf : 1));
                int* fn = malloc(sizeof(int) * (size_t)(nf > 0 ? nf : 1));
                bool empty = false;
                for (int j = 0; j < nf; j++) {
                    fc[j] = malloc(sizeof(Expr*) * (size_t)nc);
                    fn[j] = 0;
                    for (int c = 0; c < nc; c++)
                        if (fi_cand_ok_for_var(cand[c], fmap[j], V, nv, Rpart))
                            fc[j][fn[j]++] = cand[c];
                    if (fn[j] == 0) empty = true;
                }
                if (!empty) {
                    int* idx = calloc((size_t)(nf > 0 ? nf : 1), sizeof(int));
                    bool done = false;
                    while (!done && ws->n < want && spent < BUDGET) {
                        spent++;
                        Expr** frr = malloc(sizeof(Expr*) * (size_t)(nf > 0 ? nf : 1));
                        for (int j = 0; j < nf; j++)
                            frr[j] = fi_rule(V[fmap[j]], xcopy(fc[j][idx[j]]));
                        Expr* frl = expr_new_function(expr_new_symbol(SYM_List), frr, (size_t)nf);
                        free(frr);
                        Expr* xk = fi_replace(g, frl);
                        expr_free(frl);
                        if (!fi_has_nonfinite(xk) && !fi_contains_var(xk, V, nv)) {
                            Expr** allr = malloc(sizeof(Expr*) * (size_t)nv);
                            for (int j = 0; j < nf; j++)
                                allr[j] = fi_rule(V[fmap[j]], xcopy(fc[j][idx[j]]));
                            allr[nf] = fi_rule(V[k], xcopy(xk));
                            Expr* pt = expr_new_function(expr_new_symbol(SYM_List), allr, (size_t)nv);
                            free(allr);
                            if (fi_verify(expr, pt)) fi_wit_add(ws, pt); else expr_free(pt);
                        }
                        expr_free(xk);
                        if (nf == 0) break;
                        int t = 0; for (; t < nf; t++) { if (++idx[t] < fn[t]) break; idx[t] = 0; }
                        if (t == nf) done = true;
                    }
                    free(idx);
                }
                for (int j = 0; j < nf; j++) free(fc[j]);
                free(fc); free(fn); free(fmap);
                expr_free(g);
            }
            expr_free(sols);
        }
    }
    for (int i = 0; i < nc; i++) expr_free(cand[i]);
    free(cand);
    expr_free(Epart);
    if (Rpart) expr_free(Rpart);
}

/* ---- ideal saturation for declined complex polynomial systems -------------- *
 *                                                                               *
 *  Solve declines `nsdim` when the variety carries a positive-dimensional        *
 *  component the != disequations exclude -- e.g. the x=0 / y=0 components of      *
 *    x^4 y^3 - 3 x^2 y + y^4 == 0 && 4 x^3 y^3 - 6 x y == 0 && x != 0 && y != 0.  *
 *  Adjoin a Rabinowitsch slack w with w * prod(disequation LHS) == 1: this        *
 *  saturates the ideal by the disequations, so the remaining (finite) variety is  *
 *  zero-dimensional and Solve returns its roots.  Drop w, verify each root.      */
static void fi_saturate_solve(const Expr* expr, Expr** V, int nv, const Expr* dom,
                              long want, FiWit* ws) {
    if (nv < 1 || nv > 5) return;
    if (fi_is_transc_inexact(expr)) return;   /* algebraic systems only */
    int nc = is_head_sym(expr, SYM_And) ? (int)nargs(expr) : 1;
    Expr** eqs = malloc(sizeof(Expr*) * (size_t)nc);
    Expr** dis = malloc(sizeof(Expr*) * (size_t)nc);
    int neq = 0, ndis = 0;
    for (int c = 0; c < nc; c++) {
        const Expr* a = is_head_sym(expr, SYM_And) ? argn(expr, c) : expr;
        if (is_head_sym(a, SYM_Equal) && nargs(a) == 2)
            eqs[neq++] = fi_sub_eval(xcopy(argn(a, 0)), xcopy(argn(a, 1)));
        else if (is_head_sym(a, SYM_Unequal) && nargs(a) == 2)
            dis[ndis++] = fi_sub_eval(xcopy(argn(a, 0)), xcopy(argn(a, 1)));
        /* inequalities and other atoms only shrink the set: ignored here, and
         * enforced by the final verify against the original statement. */
    }
    if (neq < 1 || ndis < 1) {
        for (int i = 0; i < neq; i++) expr_free(eqs[i]);
        for (int i = 0; i < ndis; i++) expr_free(dis[i]);
        free(eqs); free(dis); return;
    }
    Expr* wsym = expr_new_symbol("FindInstance$w");
    Expr** atoms = malloc(sizeof(Expr*) * (size_t)(neq + 1));
    for (int i = 0; i < neq; i++) {
        Expr* pr[2] = { eqs[i], expr_new_integer(0) };   /* transfers eqs[i] */
        atoms[i] = expr_new_function(expr_new_symbol(SYM_Equal), pr, 2);
    }
    Expr* prod = dis[0];
    for (int i = 1; i < ndis; i++) prod = fi_mul_eval(prod, dis[i]);
    Expr* pr2[2] = { fi_mul_eval(xcopy(wsym), prod), expr_new_integer(1) };
    atoms[neq] = expr_new_function(expr_new_symbol(SYM_Equal), pr2, 2);
    Expr* sys = expr_new_function(expr_new_symbol(SYM_And), atoms, (size_t)(neq + 1));
    free(atoms); free(eqs); free(dis);

    Expr** vl = malloc(sizeof(Expr*) * (size_t)(nv + 1));
    for (int i = 0; i < nv; i++) vl[i] = xcopy(V[i]);
    vl[nv] = xcopy(wsym);
    Expr* varlist = expr_new_function(expr_new_symbol(SYM_List), vl, (size_t)(nv + 1));
    free(vl);

    Expr* sols = fi_call(SYM_Solve, sys, varlist, dom, NULL);
    expr_free(sys); expr_free(varlist);
    if (is_head_sym(sols, SYM_List)) {
        for (size_t si = 0; si < nargs(sols) && ws->n < want; si++) {
            const Expr* RL = argn(sols, si);
            if (!is_head_sym(RL, SYM_List)) continue;
            Expr** pr = malloc(sizeof(Expr*) * (size_t)nv);
            for (int i = 0; i < nv; i++) {
                Expr* val = fi_rl_value(RL, V[i]);
                pr[i] = fi_rule(V[i], val ? val : xcopy(V[i]));
            }
            Expr* pt = expr_new_function(expr_new_symbol(SYM_List), pr, (size_t)nv);
            free(pr);
            if (fi_verify(expr, pt)) fi_wit_add(ws, pt); else expr_free(pt);
        }
    }
    expr_free(sols);
    expr_free(wsym);
}

/* ---- Booleans: SAT via the DNF engine above -------------------------------- */

static Expr* fi_boolean(const Expr* expr, Expr** V, int nv, long want) {
    Dnf phi = to_dnf(expr);
    if (phi.n == 0) { dnf_free(&phi); return fi_empty_list(); }   /* unsatisfiable */
    FiWit ws = { NULL, 0, 0 };
    for (int c = 0; c < phi.n && ws.n < want; c++) {
        Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
        int* tval = malloc(sizeof(int) * (size_t)nv);   /* 0 = False (default), 1 = True */
        for (int i = 0; i < nv; i++) tval[i] = 0;
        bool ok = true;
        for (int l = 0; l < phi.cl[c].n && ok; l++) {
            const Expr* lit = phi.cl[c].lit[l];
            int vi = fi_var_index(lit, V, nv);
            if (vi >= 0) { tval[vi] = 1; continue; }
            if (is_head_sym(lit, SYM_Not) && nargs(lit) == 1) {
                int wi = fi_var_index(argn(lit, 0), V, nv);
                if (wi >= 0) { tval[wi] = 0; continue; }
            }
            ok = false;   /* a literal we cannot read as a Boolean assignment */
        }
        if (!ok) { free(outr); free(tval); continue; }
        for (int i = 0; i < nv; i++)
            outr[i] = fi_rule(V[i], expr_new_symbol(tval[i] ? SYM_True : SYM_False));
        free(tval);
        Expr* point = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
        free(outr);
        if (fi_verify(expr, point)) fi_wit_add(&ws, point);
        else expr_free(point);
    }
    dnf_free(&phi);
    if (ws.n == 0) { fi_wit_free(&ws); return NULL; }   /* had clauses, none usable */
    return fi_wit_take(&ws, want);
}

/* ---- argument parsing helpers (mirror reduce.c) ---------------------------- */

static bool fi_is_domain(const Expr* e) {
    return fi_is_sym(e, SYM_Complexes) || fi_is_sym(e, SYM_Reals)
        || fi_is_sym(e, SYM_Integers)  || fi_is_sym(e, SYM_Rationals)
        || fi_is_sym(e, SYM_Booleans);
}
static bool fi_is_option_name(const char* s) {
    return s == SYM_Modulus || s == SYM_Method
        || s == SYM_WorkingPrecision || s == SYM_RandomSeeding;
}
/* A single variable is a symbol (x) or an indexed form (c[1] = c[_]): a function
 * whose head is a symbol.  Numbers, rules and bare List heads are not variables. */
static bool fi_is_varlike(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return true;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name != SYM_List) return true;
    return false;
}
static bool fi_valid_vars(const Expr* vars) {
    if (!vars) return false;
    if (is_head_sym(vars, SYM_List)) {
        if (nargs(vars) == 0) return false;
        for (size_t i = 0; i < nargs(vars); i++)
            if (!fi_is_varlike(argn(vars, i))) return false;
        return true;
    }
    return fi_is_varlike(vars);
}
static Expr** fi_collect_vars(Expr* vars, int* nv_out) {
    if (!is_head_sym(vars, SYM_List)) {
        Expr** v = malloc(sizeof(Expr*)); v[0] = vars; *nv_out = 1; return v;
    }
    int nv = (int)nargs(vars);
    Expr** v = malloc((size_t)nv * sizeof(Expr*));
    for (int i = 0; i < nv; i++) v[i] = argn(vars, i);
    *nv_out = nv;
    return v;
}

static void fi_warn_optx(const Expr* opt) {
    const Expr* lhs = (opt && opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2)
        ? opt->data.function.args[0] : NULL;
    const char* name = (lhs && lhs->type == EXPR_SYMBOL) ? lhs->data.symbol.name : "?";
    fprintf(stderr, "FindInstance::optx: Unknown option %s in FindInstance.\n", name);
}

/* The search cascade, run under message suppression.  Takes ownership of V. */
static Expr* fi_run_search(Expr* expr, Expr* vars, Expr* dom, long nWanted,
                           Expr* modulus, Expr** V, int nv);

Expr* builtin_find_instance(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* peel trailing options (Modulus honoured; Method/WorkingPrecision/RandomSeeding
     * accepted and ignored -- we are exact & deterministic). */
    Expr* modulus = NULL;
    size_t pos_end = argc;
    while (pos_end > 0) {
        Expr* a = args[pos_end - 1];
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol.name;
            if (fi_is_option_name(name)) {
                if (name == SYM_Modulus) {
                    Expr* mv = a->data.function.args[1];
                    if (!(mv->type == EXPR_INTEGER && mv->data.integer == 0)) modulus = mv;
                }
                pos_end--; continue;
            }
            fi_warn_optx(a);
            return NULL;
        }
        break;
    }
    if (pos_end < 2) return NULL;

    Expr* expr = args[0];
    Expr* vars = args[1];
    Expr* dom = NULL;
    long nWanted = 1;
    for (size_t i = 2; i < pos_end; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_INTEGER) nWanted = a->data.integer;
        else if (fi_is_domain(a)) dom = a;
        else return NULL;   /* unrecognised positional -> leave unevaluated */
    }
    if (!fi_valid_vars(vars)) return NULL;
    if (nWanted <= 0) return fi_empty_list();

    int nv = 0;
    Expr** V = fi_collect_vars(vars, &nv);

    /* Run the cascade quietly: its internal Reduce/Solve/NMinimize/FindMinimum
     * probes are speculative (a sampled division by zero, an unsupported
     * constraint shape) and must not leak diagnostics -- Mathematica evaluates
     * its FindInstance internals under an implicit Quiet in the same way. */
    mth_msg_suppress_push();
    Expr* fi_out = fi_run_search(expr, vars, dom, nWanted, modulus, V, nv);
    mth_msg_suppress_pop();
    return fi_out;
}

static Expr* fi_run_search(Expr* expr, Expr* vars, Expr* dom, long nWanted,
                           Expr* modulus, Expr** V, int nv) {
    if (fi_is_sym(dom, SYM_Booleans)) {
        Expr* out = fi_boolean(expr, V, nv, nWanted);
        free(V);
        return out;
    }

    FiWit ws = { NULL, 0, 0 };

    /* Reduce is NOT a sound decision procedure over transcendental functions or
     * inexact numbers (it wrongly returns False for 0<x<0.001 && Sin[1/x]>0.999),
     * so for such systems its False is not trusted as {} -- a numerical witness
     * search decides instead. */
    bool transc = fi_is_transc_inexact(expr);

    /* Reduce rejects indexed variables (c[1]); Solve and the integer search
     * accept them, so skip the Reduce oracle when any listed var is not a
     * plain symbol (avoids a spurious Reduce::ivar message). */
    bool vars_indexed = false;
    for (int i = 0; i < nv; i++) if (V[i]->type != EXPR_SYMBOL) vars_indexed = true;

    /* Step 1: Reduce is the satisfiability + solution-set oracle. */
    Expr* red = vars_indexed ? NULL : fi_call(SYM_Reduce, expr, vars, dom, modulus);
    bool provably_false = fi_is_sym(red, SYM_False);
    if (!provably_false && red && !is_head_sym(red, SYM_Reduce)) {
        if (fi_is_sym(red, SYM_True)) {
            Expr** outr = malloc(sizeof(Expr*) * (size_t)nv);
            for (int i = 0; i < nv; i++) outr[i] = fi_rule(V[i], expr_new_integer(0));
            Expr* p = expr_new_function(expr_new_symbol(SYM_List), outr, (size_t)nv);
            free(outr);
            if (fi_verify(expr, p)) fi_wit_add(&ws, p); else expr_free(p);
        } else {
            /* walk each top-level Or clause */
            int nclause; Expr** clauses;
            if (is_head_sym(red, SYM_Or)) {
                nclause = (int)nargs(red);
                clauses = malloc(sizeof(Expr*) * (size_t)nclause);
                for (int i = 0; i < nclause; i++) clauses[i] = argn(red, i);
            } else {
                nclause = 1; clauses = malloc(sizeof(Expr*)); clauses[0] = red;
            }
            for (int ci = 0; ci < nclause && ws.n < nWanted; ci++) {
                Expr* clause = clauses[ci];
                int na; Expr** atoms;
                if (is_head_sym(clause, SYM_And)) {
                    na = (int)nargs(clause);
                    atoms = malloc(sizeof(Expr*) * (size_t)na);
                    for (int i = 0; i < na; i++) atoms[i] = argn(clause, i);
                } else {
                    na = 1; atoms = malloc(sizeof(Expr*)); atoms[0] = clause;
                }
                Expr* p = fi_clause_point(atoms, na, V, nv, dom, modulus);
                free(atoms);
                if (p) { if (fi_verify(expr, p)) fi_wit_add(&ws, p); else expr_free(p); }
            }
            free(clauses);
        }
    }
    expr_free(red);

    /* Reduce's False is a proof of emptiness only for exactly-decidable systems;
     * for transcendental/inexact ones it is unsound, so fall through to the
     * numerical search instead of claiming {}. */
    if (provably_false && !transc) { fi_wit_free(&ws); free(V); return fi_empty_list(); }

    /* Non-equation constraints, used to pin generated parameters. */
    Expr* Epart = NULL, *Rpart = NULL;
    fi_split_eqs(expr, &Epart, &Rpart);

    /* Step 2: Solve fallback (covers cases Reduce declines), with generated-
     * parameter instantiation -- reaches parametric Diophantine families such as
     * the Pell equation x^2 - 61 y^2 == 1 (fundamental solution at C[1] == 1). */
    if (ws.n < nWanted) {
        Expr* sols = fi_call(SYM_Solve, expr, vars, dom, modulus);
        if (is_head_sym(sols, SYM_List))
            for (size_t si = 0; si < nargs(sols) && ws.n < nWanted; si++)
                fi_materialize_rulelist(argn(sols, si), V, nv, expr, Rpart, &ws, nWanted);
        expr_free(sols);
    }

    /* Step 3: bounded integer search over the Integers -- a Diophantine witness,
     * or a finite-domain emptiness proof ({} when the box is exhausted). */
    if (ws.n < nWanted && fi_is_sym(dom, SYM_Integers)) {
        int r = fi_integer_search(expr, V, nv, nWanted, &ws);
        if (r == 0 && ws.n == 0) {
            if (Epart) expr_free(Epart);
            if (Rpart) expr_free(Rpart);
            free(V); return fi_empty_list();
        }
    }

    /* Step 4: equations-only retry -- solve the equation part to obtain a
     * parametric family, then pin the generated parameter with the remaining
     * constraints.  Reaches periodic transcendental instances such as
     * Sin[1/x] == 0 && 0 < x < 10^-5 (parameter C[1] == 15916). */
    if (ws.n < nWanted && Epart && Rpart) {
        Expr* sols = fi_call(SYM_Solve, Epart, vars, dom, NULL);
        if (is_head_sym(sols, SYM_List))
            for (size_t si = 0; si < nargs(sols) && ws.n < nWanted; si++)
                fi_materialize_rulelist(argn(sols, si), V, nv, expr, Rpart, &ws, nWanted);
        expr_free(sols);
    }
    if (Epart) expr_free(Epart);
    if (Rpart) expr_free(Rpart);

    /* Step 4b: 1-variable real transcendental root search -- bracket a pole-free
     * refactoring of a declined equation and refine at high precision.  Reaches
     * Tan[x] == x && x > 10^6 (root just below (n+1/2) Pi, n ~ 318310). */
    if (ws.n < nWanted && fi_is_sym(dom, SYM_Reals))
        fi_real_root_search(expr, V, nv, dom, nWanted, &ws);

    /* Step 4c: solve one equation for one variable, sample the rest from a
     * constraint-aware grid -- reaches systems Reduce/Solve decline as a whole
     * but that are univariate once the others are pinned
     * (c1 e^{-L1 t}+c2 e^{-L2 t}==0 && signs; (x^2-y^2)/(x^2+y^2)==1/2 && tiny box). */
    if (ws.n < nWanted && !fi_is_sym(dom, SYM_Integers))
        fi_solve_one_sample(expr, V, nv, dom, nWanted, &ws);

    /* Step 4d: ideal saturation for a declined polynomial system with != atoms
     * -- Rabinowitsch slack turns the saturated (finite) variety zero-dimensional
     * so Solve returns its roots (x^4 y^3-3x^2 y+y^4==0 && 4x^3 y^3-6xy==0 && x!=0 && y!=0). */
    if (ws.n < nWanted && !fi_is_sym(dom, SYM_Integers))
        fi_saturate_solve(expr, V, nv, dom, nWanted, &ws);

    /* Step 6: Groebner emptiness certificate for a declined polynomial system
     * (e.g. the 2x2 nilpotent M^2==0 && det!=0, which is empty). */
    if (ws.n == 0 && fi_groebner_empty(expr, vars, nv)) {
        free(V); return fi_empty_list();
    }

    /* Step 7: structured candidate sampling over Complexes/Reals -- a verify-gated
     * EXACT witness for statements Reduce/Solve decline but a concrete grid point
     * satisfies (branch-cut disequations Sqrt[z^2]!=z, Log[x y]!=Log[x]+Log[y], open
     * regions).  Tried before the numeric feasibility so an exact witness is
     * preferred over an approximate one.  Skipped for Integers (its own search) and
     * Booleans (handled above).  Complex candidates only when the domain admits them. */
    if (ws.n < nWanted && !fi_is_sym(dom, SYM_Integers)) {
        bool allow_complex = !fi_is_sym(dom, SYM_Reals) && !fi_is_sym(dom, SYM_Rationals);
        fi_sample_search(expr, V, nv, allow_complex, nWanted, &ws);
    }

    /* Step 8: numerical feasibility -- the LAST resort, for transcendental / inexact
     * Real systems no exact method or grid point reaches: 0<x<0.001 && Sin[1/x]>0.999;
     * two rational equations pinned to inexact constants (Step 8a folds those
     * definitions); the open Rastrigin < 0.1 region.  NMinimize form first, then the
     * residual/least-infeasibility form for conjunctive constraint shapes. */
    if (ws.n < nWanted && transc && !fi_is_sym(dom, SYM_Integers)
        && !fi_is_sym(dom, SYM_Complexes)) {
        Expr* w = fi_numeric_search(expr, vars);
        if (w) fi_wit_add(&ws, w);
    }
    if (ws.n < nWanted && transc && !fi_is_sym(dom, SYM_Integers)
        && !fi_is_sym(dom, SYM_Complexes)) {
        Expr* w = fi_numeric_feasibility(expr, V, nv);
        if (w) fi_wit_add(&ws, w);
    }

    free(V);
    if (ws.n == 0) { fi_wit_free(&ws); return NULL; }   /* found none, not proven empty */
    return fi_wit_take(&ws, nWanted);
}

/* ================================================================== *
 *  CylindricalDecomposition (REDUCE_PLAN.md, Phase 8)                  *
 *                                                                     *
 *  CylindricalDecomposition[expr, vars] gives a cylindrical algebraic *
 *  decomposition of the REAL solution set of `expr` -- a              *
 *  quantifier-free And/Or formula in which each variable is bounded   *
 *  cylindrically in terms of the earlier ones.  Its only semantic     *
 *  difference from Reduce is that it is Reals-only (Reduce defaults to *
 *  Complexes for equations), so it is a thin front-end: validate,     *
 *  force the Reals domain, and delegate to the Reduce builtin, whose   *
 *  Reals engine (Fourier-Motzkin / CAD / sign diagram) already emits   *
 *  the merged cylindrical formula.  If Reduce declines (leaves itself  *
 *  unevaluated), so does this -- soundness over completeness.          *
 * ================================================================== */

/* A single symbol, or a non-empty List of symbols (mirrors Reduce's
 * reduce_valid_vars; compound / indexed variables are not accepted). */
static bool cad_valid_vars(const Expr* vars) {
    if (!vars) return false;
    if (vars->type == EXPR_SYMBOL) return true;
    if (is_head_sym(vars, SYM_List)) {
        if (nargs(vars) == 0) return false;
        for (size_t i = 0; i < nargs(vars); i++)
            if (argn(vars, i)->type != EXPR_SYMBOL) return false;
        return true;
    }
    return false;
}

Expr* builtin_cylindrical_decomposition(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Peel trailing option Rules (Rule/RuleDelayed with a symbol LHS).  The
     * option NAMES are validated by the delegate Reduce, so we do not check
     * them here -- an unknown one makes Reduce decline and we decline with it. */
    size_t pos_end = argc;
    while (pos_end > 0) {
        Expr* a = args[pos_end - 1];
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol.name == SYM_Rule
                || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            pos_end--; continue;
        }
        break;
    }

    /* Positional forms: [expr, vars] or [expr, vars, Reals].  The domain is
     * always the Reals, so an explicit Reals is redundant-and-accepted; any
     * other third positional (another domain, an operation-direction arg) is
     * not supported and declines soundly (stays unevaluated). */
    if (pos_end < 2 || pos_end > 3) return NULL;
    if (pos_end == 3 && !fi_is_sym(args[2], SYM_Reals)) return NULL;

    Expr* expr = args[0];
    Expr* vars = args[1];
    if (!cad_valid_vars(vars)) return NULL;

    /* Build Reduce[expr, vars, Reals, <trailing option Rules...>] and evaluate.
     * Forwarding the option Rules verbatim reuses all of Reduce's options
     * (Modulus, Cubics, Quartics, WorkingPrecision, ...) with no per-option
     * logic here.  expr_new_function adopts the args array elements. */
    size_t nopt = argc - pos_end;
    size_t total = 3 + nopt;
    Expr** a = malloc(sizeof(Expr*) * total);
    a[0] = xcopy(expr);
    a[1] = xcopy(vars);
    a[2] = expr_new_symbol(SYM_Reals);
    for (size_t i = 0; i < nopt; i++) a[3 + i] = xcopy(args[pos_end + i]);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Reduce), a, total);
    free(a);

    /* Evaluate quietly: a declining Reduce may probe internally, and those
     * diagnostics must not be attributed to CylindricalDecomposition. */
    mth_msg_suppress_push();
    Expr* out = fi_eval_take(call);
    mth_msg_suppress_pop();

    /* Reduce declined (left itself unevaluated) -> so do we: return NULL rather
     * than echo an inner Reduce[...] under the CylindricalDecomposition head. */
    if (is_head_sym(out, SYM_Reduce)) { expr_free(out); return NULL; }
    return out;
}

void reduce_companions_init(void) {
    symtab_add_builtin("LogicalExpand", builtin_logical_expand);
    symtab_add_builtin("NotElement",    builtin_not_element);
    symtab_add_builtin("FindInstance",  builtin_find_instance);
    symtab_add_builtin("CylindricalDecomposition", builtin_cylindrical_decomposition);

    SymbolDef* d;
    d = symtab_get_def("LogicalExpand"); if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("NotElement");    if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("FindInstance");  if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("CylindricalDecomposition"); if (d) d->attributes |= ATTR_PROTECTED;

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
    symtab_set_docstring("FindInstance",
        "FindInstance[expr, vars]\n"
        "\tFinds a single instance of vars satisfying the statement expr -- a\n"
        "\tlogical combination of equations and inequalities -- returned in\n"
        "\tSolve's form {{x -> v, ...}}, or {} if none exists.  The default\n"
        "\tdomain is Complexes, or Reals when expr carries an ordering (as in\n"
        "\tReduce).\n"
        "FindInstance[expr, vars, dom]\n"
        "\tFinds an instance over dom: Complexes, Reals, Integers, Rationals,\n"
        "\tor Booleans (Boolean satisfiability).\n"
        "FindInstance[expr, vars, dom, n]\n"
        "\tFinds up to n instances (fewer if fewer exist).\n"
        "\n"
        "Every instance returned is verified against expr, so it is always a\n"
        "true solution.  Variables may be symbols or indexed forms c[i].\n"
        "FindInstance may find an instance even where Reduce cannot give a\n"
        "complete reduction -- instantiating parametric Diophantine families,\n"
        "searching a bounded integer box over the Integers, and (for\n"
        "transcendental or inexact Real systems) a numerical feasibility search.\n"
        "It returns {} only when the set is provably empty -- including a\n"
        "Groebner certificate for declined polynomial systems -- and stays\n"
        "unevaluated otherwise.  Modulus -> p over Z/pZ.");
    symtab_set_docstring("CylindricalDecomposition",
        "CylindricalDecomposition[expr, vars]\n"
        "\tGives a cylindrical algebraic decomposition of the real solution set\n"
        "\tof expr -- a logical combination of polynomial equations and\n"
        "\tinequalities -- as a quantifier-free And/Or formula in which each\n"
        "\tvariable is bounded cylindrically in terms of the earlier ones, e.g.\n"
        "\tCylindricalDecomposition[x^2 + y^2 <= 1, {x, y}] gives\n"
        "\t-1 <= x <= 1 && -Sqrt[1 - x^2] <= y <= Sqrt[1 - x^2].  The domain is\n"
        "\talways the Reals.  Returns True / False when the statement decides,\n"
        "\tand stays unevaluated when the decomposition cannot be computed\n"
        "\texactly (an undecidable sign, or a positive-dimensional system with\n"
        "\tirrational fibres).  Reduce's options (Modulus, Cubics, Quartics,\n"
        "\tWorkingPrecision, ...) may be given and are forwarded.");
}
