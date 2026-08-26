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

static bool fi_is_var(const Expr* e, const Expr* v) {
    return e && e->type == EXPR_SYMBOL && v && v->type == EXPR_SYMBOL
        && e->data.symbol.name == v->data.symbol.name;
}
static int fi_var_index(const Expr* e, Expr** V, int nv) {
    if (!e || e->type != EXPR_SYMBOL) return -1;
    for (int i = 0; i < nv; i++)
        if (e->data.symbol.name == V[i]->data.symbol.name) return i;
    return -1;
}
/* Does any listed variable appear anywhere in e? */
static bool fi_contains_var(const Expr* e, Expr** V, int nv) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return fi_var_index(e, V, nv) >= 0;
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

/* The single soundness gate: does expr hold at the candidate point? */
static bool fi_verify(const Expr* expr, const Expr* point /* List[Rule..] */) {
    Expr* sub = fi_replace(expr, point);
    bool ok = fi_is_sym(sub, SYM_True);
    expr_free(sub);
    return ok;
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
static bool fi_valid_vars(const Expr* vars) {
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
static Expr** fi_collect_vars(Expr* vars, int* nv_out) {
    if (vars->type == EXPR_SYMBOL) {
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

    if (fi_is_sym(dom, SYM_Booleans)) {
        Expr* out = fi_boolean(expr, V, nv, nWanted);
        free(V);
        return out;
    }

    FiWit ws = { NULL, 0, 0 };

    /* Step 1: Reduce is the satisfiability + solution-set oracle. */
    Expr* red = fi_call(SYM_Reduce, expr, vars, dom, modulus);
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

    if (provably_false) { fi_wit_free(&ws); free(V); return fi_empty_list(); }

    /* Step 2: Solve fallback (covers cases Reduce declines). */
    if (ws.n < nWanted) {
        Expr* sols = fi_call(SYM_Solve, expr, vars, dom, modulus);
        if (is_head_sym(sols, SYM_List)) {
            for (size_t si = 0; si < nargs(sols) && ws.n < nWanted; si++) {
                const Expr* rl = argn(sols, si);
                static const long grid[] = { 0, 1, -1, 2, -2, 3, -3 };
                for (size_t gi = 0; gi < sizeof(grid) / sizeof(grid[0]) && ws.n < nWanted; gi++) {
                    Expr* p = fi_solve_point(rl, V, nv, grid[gi]);
                    if (!p) continue;
                    if (fi_verify(expr, p)) fi_wit_add(&ws, p); else expr_free(p);
                }
            }
        }
        expr_free(sols);
    }

    free(V);
    if (ws.n == 0) { fi_wit_free(&ws); return NULL; }   /* found none, not proven empty */
    return fi_wit_take(&ws, nWanted);
}

void reduce_companions_init(void) {
    symtab_add_builtin("LogicalExpand", builtin_logical_expand);
    symtab_add_builtin("NotElement",    builtin_not_element);
    symtab_add_builtin("FindInstance",  builtin_find_instance);

    SymbolDef* d;
    d = symtab_get_def("LogicalExpand"); if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("NotElement");    if (d) d->attributes |= ATTR_PROTECTED;
    d = symtab_get_def("FindInstance");  if (d) d->attributes |= ATTR_PROTECTED;

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
        "true solution.  FindInstance may find an instance even where Reduce\n"
        "cannot give a complete reduction; it returns {} only when the set is\n"
        "provably empty, and stays unevaluated when it can neither exhibit an\n"
        "instance nor prove emptiness.  Option Modulus -> p solves over Z/pZ.");
}
