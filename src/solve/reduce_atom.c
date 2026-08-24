/*
 * reduce_atom.c
 *
 * Atom canonicalisation and input parsing for `Reduce` (see reduce_form.h).
 *
 *  - reduce_atom_canonicalize:  one relational Expr -> canonical RAtom
 *                               (poly REL 0, with > / >= swapped to < / <=).
 *  - reduce_atom_decide_constant: decide a constant atom via evaluate().
 *  - reduce_form_from_expr:      a logical combination (&&, ||, !, Implies,
 *                               binary Xor, chained Inequality, relational
 *                               leaves) -> DNF RForm, via mutually-recursive
 *                               positive/negative builders (De Morgan).
 *
 * The RForm algebra itself lives in reduce_form.c.
 */
#include "reduce_form.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* True iff `e` mentions the interned symbol `name` anywhere. */
static bool expr_contains_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (expr_contains_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_contains_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* Index of the first ambient variable that appears in `e`, else -1. */
static int first_var_in(const Expr* e, Expr** vars, int nv) {
    for (int i = 0; i < nv; i++) {
        if (vars[i]->type == EXPR_SYMBOL
            && expr_contains_symbol(e, vars[i]->data.symbol.name))
            return i;
    }
    return -1;
}

/* For EQ/NE atoms the overall sign of the polynomial is irrelevant (p == 0 iff
 * -p == 0), so flip a leading negative.  This keeps `-b == 0` from printing
 * differently than `b == 0` and canonicalises the parametric-branch conditions.
 * Consumes `poly`; returns an owned (possibly negated) polynomial. */
static Expr* sign_normalize(Expr* poly) {
    bool neg = false;
    if (poly->type == EXPR_INTEGER) neg = poly->data.integer < 0;
    else if (poly->type == EXPR_REAL) neg = poly->data.real < 0.0;
    else if (poly->type == EXPR_BIGINT) neg = mpz_sgn(poly->data.bigint) < 0;
    else if (is_head(poly, SYM_Rational) && poly->data.function.arg_count == 2) {
        const Expr* num = poly->data.function.args[0];
        if (num->type == EXPR_INTEGER) neg = num->data.integer < 0;
        else if (num->type == EXPR_BIGINT) neg = mpz_sgn(num->data.bigint) < 0;
    } else if (is_head(poly, SYM_Times) && poly->data.function.arg_count >= 1) {
        const Expr* f = poly->data.function.args[0];
        if (f->type == EXPR_INTEGER) neg = f->data.integer < 0;
        else if (f->type == EXPR_REAL) neg = f->data.real < 0.0;
        else if (f->type == EXPR_BIGINT) neg = mpz_sgn(f->data.bigint) < 0;
        else if (is_head(f, SYM_Rational) && f->data.function.arg_count == 2
                 && f->data.function.args[0]->type == EXPR_INTEGER)
            neg = f->data.function.args[0]->data.integer < 0;
    }
    if (!neg) return poly;
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), poly }, 2));
}

/* True iff `e` contains a transcendental head whose `Together` simplification
 * invokes complex-branch identities that are unsound over the reals (e.g.
 * Together[Log[x^2] - 2 Log[-x]] collapses to the constant -2 I Pi, which is
 * WRONG for real x < 0 where the difference is 0).  For such atoms canonical_poly
 * keeps the plain evaluated difference instead of running Together, so the real
 * sign-diagram engine sees the honest expression. */
static bool contains_branch_transcendental(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Log || h == SYM_ArcSin || h == SYM_ArcCos || h == SYM_ArcTan
            || h == SYM_ArcCot || h == SYM_ArcSec || h == SYM_ArcCsc
            || h == SYM_ArcSinh || h == SYM_ArcCosh || h == SYM_ArcTanh
            || h == SYM_ArcCoth || h == SYM_ArcSech || h == SYM_ArcCsch)
            return true;
    }
    if (contains_branch_transcendental(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_branch_transcendental(e->data.function.args[i])) return true;
    return false;
}

/* den := Denominator[Together[lhs - rhs]], evaluated.  Returns the owned
 * denominator polynomial when it is NOT a nonzero numeric constant -- i.e. the
 * relation is a genuine rational function whose poles and denominator sign
 * matter -- and NULL when the denominator is constant (an ordinary polynomial
 * relation).  When non-NULL, clearing the denominator to form the numerator
 * dropped the pole/sign structure, so the multivariate engines must decline and
 * the univariate sign diagram must reason with the returned denominator. */
static Expr* canonical_denom(const Expr* lhs, const Expr* rhs) {
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)lhs), expr_copy((Expr*)rhs) }, 2);
    /* Branch-cut transcendentals: skip Together, treat the denominator as the
     * trivial constant 1 (no pole structure to recover). */
    Expr* d0 = eval_and_free(expr_copy(diff));
    if (contains_branch_transcendental(d0)) { expr_free(d0); expr_free(diff); return NULL; }
    expr_free(d0);
    Expr* tog = expr_new_function(expr_new_symbol(SYM_Together), (Expr*[]){ diff }, 1);
    Expr* den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator),
        (Expr*[]){ tog }, 1));
    bool constant;
    switch (den->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            constant = true; break;
        case EXPR_FUNCTION:
            constant = is_head(den, SYM_Rational);    /* Rational[p,q] is constant */
            break;
        default: constant = false; break;             /* a symbol, etc. */
    }
    if (constant) { expr_free(den); return NULL; }
    return den;
}

/* poly := Numerator[Together[lhs - rhs]], evaluated. */
static Expr* canonical_poly(const Expr* lhs, const Expr* rhs) {
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)lhs), expr_copy((Expr*)rhs) }, 2);
    /* Branch-cut transcendentals (Log, inverse trig/hyperbolic): Together would
     * apply complex-branch identities unsound over the reals, so keep the plain
     * evaluated difference as the atom polynomial. */
    Expr* d0 = eval_and_free(diff);
    if (contains_branch_transcendental(d0)) return d0;
    Expr* tog  = expr_new_function(expr_new_symbol(SYM_Together),
        (Expr*[]){ d0 }, 1);
    Expr* num  = expr_new_function(expr_new_symbol(SYM_Numerator),
        (Expr*[]){ tog }, 1);
    return eval_and_free(num);
}

/* ------------------------------------------------------------------ *
 *  Atom canonicalisation                                              *
 * ------------------------------------------------------------------ */

RAtom reduce_atom_canonicalize(const Expr* rel, Expr** vars, int nv, bool* ok) {
    RAtom a;
    a.poly = NULL; a.denom = NULL; a.elem_dom = NULL; a.display = NULL; a.rel = R_EQ;
    a.nonconst_denom = false;
    a.main_var = -1; a.deg_main = -1; a.is_linear = false;

    if (!rel || rel->type != EXPR_FUNCTION
        || rel->data.function.head->type != EXPR_SYMBOL
        || rel->data.function.arg_count != 2) {
        if (ok) *ok = false;
        return a;
    }
    const char* h = rel->data.function.head->data.symbol.name;
    Expr* l = rel->data.function.args[0];
    Expr* r = rel->data.function.args[1];

    if (h == SYM_Element) {
        a.rel = R_ELEM;
        a.poly = expr_copy(l);
        a.elem_dom = expr_copy(r);
        a.main_var = first_var_in(a.poly, vars, nv);
        return a;
    }

    /* Map head -> relation, swapping operands for the > / >= forms so only
     * LT / LE / EQ / NE ever leave this function. */
    if (h == SYM_Equal)             { a.rel = R_EQ; a.poly = canonical_poly(l, r); }
    else if (h == SYM_Unequal)      { a.rel = R_NE; a.poly = canonical_poly(l, r); }
    else if (h == SYM_Less)         { a.rel = R_LT; a.poly = canonical_poly(l, r); }
    else if (h == SYM_LessEqual)    { a.rel = R_LE; a.poly = canonical_poly(l, r); }
    else if (h == SYM_Greater)      { a.rel = R_LT; a.poly = canonical_poly(r, l); }
    else if (h == SYM_GreaterEqual) { a.rel = R_LE; a.poly = canonical_poly(r, l); }
    else { if (ok) *ok = false; return a; }

    a.denom = canonical_denom(l, r);
    a.nonconst_denom = (a.denom != NULL);
    if (a.rel == R_EQ || a.rel == R_NE) a.poly = sign_normalize(a.poly);

    a.main_var = first_var_in(a.poly, vars, nv);
    return a;
}

RAtom ratom_solved(const Expr* var, const Expr* value, Expr** vars, int nv) {
    /* semantic content: canonical(var - value) rel EQ */
    Expr* eqrel = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)var), expr_copy((Expr*)value) }, 2);
    bool ok = true;
    RAtom a = reduce_atom_canonicalize(eqrel, vars, nv, &ok);
    expr_free(eqrel);
    /* display override: the solved form  var == value  verbatim */
    a.display = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)var), expr_copy((Expr*)value) }, 2);
    return a;
}

/* ------------------------------------------------------------------ *
 *  Constant-atom decision                                             *
 * ------------------------------------------------------------------ */

int reduce_atom_decide_constant(const RAtom* a) {
    /* A relation whose canonicalisation cleared a variable denominator is not
     * equivalent to its polynomial numerator: the denominator's sign and its
     * poles matter.  So a nonzero-constant numerator does NOT decide an ordering
     * (< <=) or an Unequal (p/q != 0 still excludes the poles, e.g. 1/x != 0 is
     * x != 0, not True).  Keep these undecided so the univariate sign diagram
     * resolves them and the multivariate engines can decline.  An Equal
     * (p/q == 0) with a nonzero-constant numerator IS decidably false, and is
     * left to the evaluate() path below. */
    if (a->nonconst_denom && (a->rel == R_LT || a->rel == R_LE || a->rel == R_NE))
        return -1;
    Expr* rel;
    if (a->rel == R_ELEM) {
        rel = expr_new_function(expr_new_symbol(SYM_Element),
            (Expr*[]){ expr_copy(a->poly), expr_copy(a->elem_dom) }, 2);
    } else {
        const char* head;
        switch (a->rel) {
            case R_EQ: head = SYM_Equal;     break;
            case R_NE: head = SYM_Unequal;   break;
            case R_LT: head = SYM_Less;      break;
            case R_LE: head = SYM_LessEqual; break;
            default:   head = SYM_Equal;     break;
        }
        rel = expr_new_function(expr_new_symbol(head),
            (Expr*[]){ expr_copy(a->poly), expr_new_integer(0) }, 2);
    }
    Expr* v = eval_and_free(rel);
    int result = -1;
    if (is_sym(v, SYM_True)) result = 1;
    else if (is_sym(v, SYM_False)) result = 0;
    expr_free(v);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Input parsing: logical combination -> DNF                          *
 * ------------------------------------------------------------------ */

static RForm* build_pos(const Expr* e, Expr** vars, int nv, bool* ok);
static RForm* build_neg(const Expr* e, Expr** vars, int nv, bool* ok);

/* A single relational leaf -> one-atom form (optionally negated). */
static RForm* leaf_form(const Expr* e, bool neg, Expr** vars, int nv, bool* ok) {
    RAtom a = reduce_atom_canonicalize(e, vars, nv, ok);
    if (ok && !*ok) { ratom_free_members(&a); return rform_false(); }
    if (!neg) return rform_single_atom(a);
    RAtom na = ratom_negate(&a, ok);
    ratom_free_members(&a);
    if (ok && !*ok) { ratom_free_members(&na); return rform_false(); }
    return rform_single_atom(na);
}

/* Rewrite a chained Inequality[e0, op0, e1, op1, ...] into an unevaluated
 * And[op0[e0,e1], op1[e1,e2], ...] (freshly owned).  Returns NULL on a
 * malformed chain. */
static Expr* inequality_to_and(const Expr* e) {
    size_t n = e->data.function.arg_count;
    if (n < 3 || (n % 2) == 0) return NULL;
    size_t pairs = (n - 1) / 2;
    Expr** rels = malloc(pairs * sizeof(Expr*));
    for (size_t k = 0; k < pairs; k++) {
        Expr* op    = e->data.function.args[2 * k + 1];
        Expr* left  = e->data.function.args[2 * k];
        Expr* right = e->data.function.args[2 * k + 2];
        rels[k] = expr_new_function(expr_copy(op),
            (Expr*[]){ expr_copy(left), expr_copy(right) }, 2);
    }
    Expr* out;
    if (pairs == 1) { out = rels[0]; }
    else { out = expr_new_function(expr_new_symbol(SYM_And), rels, pairs); }
    free(rels);
    return out;
}

static RForm* build_pos(const Expr* e, Expr** vars, int nv, bool* ok) {
    if (is_sym(e, SYM_True))  return rform_true();
    if (is_sym(e, SYM_False)) return rform_false();

    if (is_head(e, SYM_And)) {
        RForm* acc = rform_true();
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            acc = rform_and(acc, build_pos(e->data.function.args[i], vars, nv, ok));
        return acc;
    }
    if (is_head(e, SYM_Or)) {
        RForm* acc = rform_false();
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            acc = rform_or(acc, build_pos(e->data.function.args[i], vars, nv, ok));
        return acc;
    }
    if (is_head(e, SYM_Not) && e->data.function.arg_count == 1)
        return build_neg(e->data.function.args[0], vars, nv, ok);

    if (is_head(e, SYM_Implies) && e->data.function.arg_count == 2) {
        /* a => b   ==   !a || b */
        RForm* na = build_neg(e->data.function.args[0], vars, nv, ok);
        RForm* b  = build_pos(e->data.function.args[1], vars, nv, ok);
        return rform_or(na, b);
    }
    if (is_head(e, SYM_Xor) && e->data.function.arg_count == 2) {
        /* a XOR b  ==  (a && !b) || (!a && b) */
        Expr* a0 = e->data.function.args[0];
        Expr* b0 = e->data.function.args[1];
        RForm* t1 = rform_and(build_pos(a0, vars, nv, ok), build_neg(b0, vars, nv, ok));
        RForm* t2 = rform_and(build_neg(a0, vars, nv, ok), build_pos(b0, vars, nv, ok));
        return rform_or(t1, t2);
    }
    if (is_head(e, SYM_Inequality)) {
        Expr* conj = inequality_to_and(e);
        if (!conj) { if (ok) *ok = false; return rform_false(); }
        RForm* f = build_pos(conj, vars, nv, ok);
        expr_free(conj);
        return f;
    }

    /* Relational leaf. */
    if (is_head(e, SYM_Equal) || is_head(e, SYM_Unequal)
        || is_head(e, SYM_Less) || is_head(e, SYM_LessEqual)
        || is_head(e, SYM_Greater) || is_head(e, SYM_GreaterEqual)
        || is_head(e, SYM_Element)) {
        return leaf_form(e, false, vars, nv, ok);
    }

    if (ok) *ok = false;
    return rform_false();
}

static RForm* build_neg(const Expr* e, Expr** vars, int nv, bool* ok) {
    if (is_sym(e, SYM_True))  return rform_false();
    if (is_sym(e, SYM_False)) return rform_true();

    if (is_head(e, SYM_And)) {
        /* !(∧ ai) == ∨ !ai */
        RForm* acc = rform_false();
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            acc = rform_or(acc, build_neg(e->data.function.args[i], vars, nv, ok));
        return acc;
    }
    if (is_head(e, SYM_Or)) {
        /* !(∨ ai) == ∧ !ai */
        RForm* acc = rform_true();
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            acc = rform_and(acc, build_neg(e->data.function.args[i], vars, nv, ok));
        return acc;
    }
    if (is_head(e, SYM_Not) && e->data.function.arg_count == 1)
        return build_pos(e->data.function.args[0], vars, nv, ok);

    if (is_head(e, SYM_Implies) && e->data.function.arg_count == 2) {
        /* !(a => b)  ==  a && !b */
        RForm* a = build_pos(e->data.function.args[0], vars, nv, ok);
        RForm* nb = build_neg(e->data.function.args[1], vars, nv, ok);
        return rform_and(a, nb);
    }
    if (is_head(e, SYM_Xor) && e->data.function.arg_count == 2) {
        /* !(a XOR b)  ==  (a && b) || (!a && !b) */
        Expr* a0 = e->data.function.args[0];
        Expr* b0 = e->data.function.args[1];
        RForm* t1 = rform_and(build_pos(a0, vars, nv, ok), build_pos(b0, vars, nv, ok));
        RForm* t2 = rform_and(build_neg(a0, vars, nv, ok), build_neg(b0, vars, nv, ok));
        return rform_or(t1, t2);
    }
    if (is_head(e, SYM_Inequality)) {
        Expr* conj = inequality_to_and(e);
        if (!conj) { if (ok) *ok = false; return rform_false(); }
        RForm* f = build_neg(conj, vars, nv, ok);
        expr_free(conj);
        return f;
    }

    if (is_head(e, SYM_Equal) || is_head(e, SYM_Unequal)
        || is_head(e, SYM_Less) || is_head(e, SYM_LessEqual)
        || is_head(e, SYM_Greater) || is_head(e, SYM_GreaterEqual)
        || is_head(e, SYM_Element)) {
        return leaf_form(e, true, vars, nv, ok);
    }

    if (ok) *ok = false;
    return rform_false();
}

RForm* reduce_form_from_expr(const Expr* e, Expr** vars, int nv, bool* ok) {
    if (ok) *ok = true;
    return build_pos(e, vars, nv, ok);
}
