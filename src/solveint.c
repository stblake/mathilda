/*
 * solveint.c
 *
 * Integer-domain (Diophantine) solving pre-pass for
 * Solve[eqns && constraints, vars, Integers].  See solveint.h for the
 * contract.
 *
 * Phase 1 engine:
 *   Stage A  separate equations from inequality / ordering / disequation
 *            constraints; convert each equation residual to a sparse
 *            integer MPoly (denominators cleared).
 *   Stage B  derive a finite integer box [lo_i, hi_i] per variable by a
 *            fixpoint of: explicit bounds, ordering propagation, and an
 *            interval-positivity rule (a sign-definite term is bounded by
 *            the rest of its (in)equality).  Decline if any variable stays
 *            unbounded (that is a later phase: linear-parametric / Pell /
 *            research-grade forms).
 *   Stage C  recursive elimination over the search variables with an exact
 *            univariate leaf (integer k-th root / quadratic discriminant /
 *            rational-root), every candidate re-verified against the
 *            original conjunction before it is emitted.
 *
 * Only necessary conditions are ever used to tighten a bound, so an
 * exhausted finite search that finds nothing returns the empty List {}
 * (a proof of no integer solutions); an input we cannot bound or evaluate
 * returns NULL (Solve stays unevaluated).
 */

#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (mirror the sibling solve specialists).  *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }
static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}
static Expr* mk_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* ------------------------------------------------------------------ *
 *  Limits.                                                            *
 * ------------------------------------------------------------------ */

#define SI_MAX_VARS   10
/* Refuse (return unevaluated) rather than enumerate a search box with more
 * than this many raw nodes -- prevents a hang without ever truncating a
 * result.  Ordering constraints usually shrink the actual walk well below
 * this.  Sized to admit the largest genuinely-bounded Phase-1 workloads
 * (~10^8 ordered tuples) while declining families that want a closed-form
 * solver instead (a Pell orbit, a wide linear lattice). */
#define SI_MAX_NODES  200000000LL

/* ------------------------------------------------------------------ *
 *  Constraint store.                                                  *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*  poly_src;   /* borrowed: the (in)equation as an Expr (unused after conv) */
    MPoly* Q;          /* residual polynomial; interpret with `kind` */
    int    kind;       /* SI_EQ (Q == 0, both orientations) or SI_LE (Q <= 0) */
} BoundConstraint;

#define SI_EQ 0
#define SI_LE 1

typedef struct {
    Expr** var;                 /* n borrowed symbols */
    int    n;

    int64_t lo[SI_MAX_VARS];
    int64_t hi[SI_MAX_VARS];
    bool    has_lo[SI_MAX_VARS];
    bool    has_hi[SI_MAX_VARS];

    /* orderings: var[a] < var[b] (strict) or <= (non-strict) */
    int   ord_a[SI_MAX_VARS * SI_MAX_VARS];
    int   ord_b[SI_MAX_VARS * SI_MAX_VARS];
    bool  ord_strict[SI_MAX_VARS * SI_MAX_VARS];
    int   n_ord;

    /* abs-orderings: |var[a]| < |var[b]| (strict) or <= (non-strict).  A
     * separate store because it constrains magnitudes, not signed values:
     * bound propagation and candidate filtering treat it via |.|, and it never
     * reduces to a signed ordering (both signs of each variable survive). */
    int   abs_ord_a[SI_MAX_VARS * SI_MAX_VARS];
    int   abs_ord_b[SI_MAX_VARS * SI_MAX_VARS];
    bool  abs_ord_strict[SI_MAX_VARS * SI_MAX_VARS];
    int   n_abs_ord;

    /* disequations: var[a] != var[b] */
    int   neq_a[SI_MAX_VARS * SI_MAX_VARS];
    int   neq_b[SI_MAX_VARS * SI_MAX_VARS];
    int   n_neq;

    /* equation MPolys used both for solving and for bounding */
    MPoly* eq[SI_MAX_VARS * 2];
    int    neq;

    /* every relation normalised for the bounder */
    BoundConstraint bc[SI_MAX_VARS * 4];
    int    nbc;

    /* True iff every constraint is fully represented by the bound / ordering /
     * disequation store, so a candidate can be checked numerically without a
     * symbolic re-evaluation of the original conjunction. */
    bool   all_captured;

    Expr*  original;            /* borrowed: full conjunction, for verification */
} SICtx;

/* Fast per-candidate check (numeric when every constraint is store-captured,
 * else a full symbolic evaluation).  Defined below; forward-declared so the
 * leaf search can use it. */
static bool si_verify(SICtx* c, const int64_t* vals);

/* ------------------------------------------------------------------ *
 *  Small utilities.                                                   *
 * ------------------------------------------------------------------ */

static int find_var_index(const SICtx* c, const char* name) {
    for (int i = 0; i < c->n; i++)
        if (c->var[i]->data.symbol.name == name) return i;
    return -1;
}

static bool is_sym(const Expr* e, const char* head) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == head;
}

static bool is_fun(const Expr* e, const char* head, size_t argc) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == head
        && e->data.function.arg_count == argc;
}

/* Read an integer-valued Expr (Integer or BigInt fitting int64) into *out.
 * Returns false for anything else (Real, Rational, symbol, too-big BigInt). */
static bool expr_as_i64(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (mpz_fits_slong_p(e->data.bigint)) { *out = mpz_get_si(e->data.bigint); return true; }
    }
    return false;
}

/* Evaluate an MPoly at an integer assignment `vals` into `out` (pre-init'd).
 * Pure GMP arithmetic -- the per-candidate hot path uses this instead of a
 * symbolic re-evaluation. */
static void si_eval_mpoly(const MPoly* p, const int64_t* vals, mpz_t out) {
    mpz_set_ui(out, 0);
    mpz_t term; mpz_init(term);
    for (size_t t = 0; t < p->n_terms; t++) {
        const int* ex = p->exps + t * (size_t)p->n_vars;
        mpz_set(term, p->coefs[t]);
        for (int v = 0; v < p->n_vars; v++)
            for (int e = 0; e < ex[v]; e++) mpz_mul_si(term, term, (long)vals[v]);
        mpz_add(out, out, term);
    }
    mpz_clear(term);
}

/* ------------------------------------------------------------------ *
 *  Stage A: residual -> integer MPoly (denominators cleared).         *
 * ------------------------------------------------------------------ */

/* Build Expand[Numerator[Together[a - b]]] and convert to an MPoly over
 * the ctx variable ordering.  Returns NULL if the result is not a
 * polynomial with integer coefficients in those variables. */
static MPoly* relation_to_mpoly(const SICtx* c, Expr* a, Expr* b) {
    Expr* diff = mk_fn2("Plus", expr_copy(a),
                        mk_fn2("Times", mk_int(-1), expr_copy(b)));
    Expr* tog = eval_and_free(internal_together((Expr*[]){ diff }, 1));
    Expr* num = eval_and_free(internal_numerator((Expr*[]){ tog }, 1));
    Expr* exp = eval_and_free(internal_expand((Expr*[]){ num }, 1));
    MPoly* Q = expr_to_mpoly(exp, c->var, c->n);
    expr_free(exp);
    return Q;
}

/* ------------------------------------------------------------------ *
 *  Stage A: parse the conjunction.                                    *
 * ------------------------------------------------------------------ */

static void add_ordering(SICtx* c, int a, int b, bool strict) {
    if (a < 0 || b < 0 || a == b) return;
    if (c->n_ord >= (int)(sizeof(c->ord_a) / sizeof(c->ord_a[0]))) return;
    c->ord_a[c->n_ord] = a;
    c->ord_b[c->n_ord] = b;
    c->ord_strict[c->n_ord] = strict;
    c->n_ord++;
}

/* Record  |var[a]| (< / <=) |var[b]|. */
static void add_abs_ordering(SICtx* c, int a, int b, bool strict) {
    if (a < 0 || b < 0 || a == b) return;
    if (c->n_abs_ord >= (int)(sizeof(c->abs_ord_a) / sizeof(c->abs_ord_a[0]))) return;
    c->abs_ord_a[c->n_abs_ord] = a;
    c->abs_ord_b[c->n_abs_ord] = b;
    c->abs_ord_strict[c->n_abs_ord] = strict;
    c->n_abs_ord++;
}

static void tighten_lo(SICtx* c, int i, int64_t v) {
    if (!c->has_lo[i] || v > c->lo[i]) { c->lo[i] = v; c->has_lo[i] = true; }
}
static void tighten_hi(SICtx* c, int i, int64_t v) {
    if (!c->has_hi[i] || v < c->hi[i]) { c->hi[i] = v; c->has_hi[i] = true; }
}

/* Register a relation `lhs OP rhs` (OP one of <, <=, >, >=) for bounding and
 * for the structured bound store.  `dir` is +1 for a "less" relation
 * (lhs < rhs / lhs <= rhs) and -1 for a "greater" relation. `strict` marks
 * < / >.  Returns true if consumed. */
static void register_inequality(SICtx* c, Expr* lhs, Expr* rhs, int dir, bool strict) {
    /* Normalise to L (<) R with L = the smaller side. */
    Expr* L = (dir > 0) ? lhs : rhs;
    Expr* R = (dir > 0) ? rhs : lhs;
    bool captured = false;

    /* var (rel) const  and  const (rel) var  -> explicit lo/hi. */
    int64_t k;
    if (L->type == EXPR_SYMBOL && expr_as_i64(R, &k)) {
        int i = find_var_index(c, L->data.symbol.name);
        if (i >= 0) { tighten_hi(c, i, strict ? k - 1 : k); captured = true; }  /* var <= R */
    }
    if (R->type == EXPR_SYMBOL && expr_as_i64(L, &k)) {
        int i = find_var_index(c, R->data.symbol.name);
        if (i >= 0) { tighten_lo(c, i, strict ? k + 1 : k); captured = true; }  /* var >= L */
    }
    /* Abs[var] (rel) const -> symmetric box. */
    if (is_fun(L, SYM_Abs, 1) && L->data.function.args[0]->type == EXPR_SYMBOL
        && expr_as_i64(R, &k)) {
        int i = find_var_index(c, L->data.function.args[0]->data.symbol.name);
        if (i >= 0) { int64_t b = strict ? k - 1 : k; tighten_hi(c, i, b); tighten_lo(c, i, -b); captured = true; }
    }
    /* Abs[var] (rel) Abs[var] -> abs-ordering (|a| < |b| or |a| <= |b|).  This
     * is what makes  Abs[x] < Abs[y] < Abs[z] < B  a bounded, ordered search:
     * the magnitude chain propagates a box onto every variable (see
     * propagate_abs_orderings) and filters the result to the ordered subset. */
    if (is_fun(L, SYM_Abs, 1) && is_fun(R, SYM_Abs, 1)
        && L->data.function.args[0]->type == EXPR_SYMBOL
        && R->data.function.args[0]->type == EXPR_SYMBOL) {
        int a = find_var_index(c, L->data.function.args[0]->data.symbol.name);
        int b = find_var_index(c, R->data.function.args[0]->data.symbol.name);
        if (a >= 0 && b >= 0) { add_abs_ordering(c, a, b, strict); captured = true; }
    }
    /* var (rel) var -> ordering. */
    if (L->type == EXPR_SYMBOL && R->type == EXPR_SYMBOL) {
        int a = find_var_index(c, L->data.symbol.name);
        int b = find_var_index(c, R->data.symbol.name);
        add_ordering(c, a, b, strict);
        captured = true;
    }
    /* An inequality on an expression (e.g. x^3 + y^3 < 10^5) still feeds the
     * bounder but is not a store-checkable constraint. */
    if (!captured) c->all_captured = false;
    /* Feed the polynomial form to the bounder as Q <= 0 (Q = L - R). */
    MPoly* Q = relation_to_mpoly(c, L, R);
    if (Q && c->nbc < (int)(sizeof(c->bc) / sizeof(c->bc[0]))) {
        c->bc[c->nbc].Q = Q; c->bc[c->nbc].kind = SI_LE; c->nbc++;
    } else if (Q) {
        mpoly_free(Q);
    }
}

/* Classify one conjunct.  Equations feed both the solver and the bounder;
 * inequalities feed the bounder and the structured store; Unequal / opaque
 * conjuncts are ignored here and caught by the final verification.
 * Returns false only on a hard failure (a non-polynomial equation). */
static bool classify_conjunct(SICtx* c, Expr* e) {
    if (is_sym(e, SYM_True)) return true;

    if (is_fun(e, SYM_Equal, 2)) {
        MPoly* Q = relation_to_mpoly(c, e->data.function.args[0],
                                     e->data.function.args[1]);
        if (!Q) return false;                       /* non-polynomial equation */
        if (mpoly_is_zero(Q)) { mpoly_free(Q); return true; }  /* 0 == 0 */
        if (c->neq >= (int)(sizeof(c->eq) / sizeof(c->eq[0]))) { mpoly_free(Q); return false; }
        c->eq[c->neq++] = mpoly_copy(Q);
        if (c->nbc < (int)(sizeof(c->bc) / sizeof(c->bc[0]))) {
            c->bc[c->nbc].Q = Q; c->bc[c->nbc].kind = SI_EQ; c->nbc++;
        } else { mpoly_free(Q); }
        return true;
    }
    if (is_fun(e, SYM_Less, 2))
        { register_inequality(c, e->data.function.args[0], e->data.function.args[1], +1, true);  return true; }
    if (is_fun(e, SYM_LessEqual, 2))
        { register_inequality(c, e->data.function.args[0], e->data.function.args[1], +1, false); return true; }
    if (is_fun(e, SYM_Greater, 2))
        { register_inequality(c, e->data.function.args[0], e->data.function.args[1], -1, true);  return true; }
    if (is_fun(e, SYM_GreaterEqual, 2))
        { register_inequality(c, e->data.function.args[0], e->data.function.args[1], -1, false); return true; }

    /* Chained Inequality[e0, op1, e1, op2, e2, ...] -> pairwise. */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Inequality
        && e->data.function.arg_count >= 3 && (e->data.function.arg_count % 2) == 1) {
        size_t m = e->data.function.arg_count;
        for (size_t i = 0; i + 2 < m; i += 2) {
            Expr* a = e->data.function.args[i];
            Expr* op = e->data.function.args[i + 1];
            Expr* b = e->data.function.args[i + 2];
            if (is_sym(op, SYM_Less))        register_inequality(c, a, b, +1, true);
            else if (is_sym(op, SYM_LessEqual)) register_inequality(c, a, b, +1, false);
            else if (is_sym(op, SYM_Greater))   register_inequality(c, a, b, -1, true);
            else if (is_sym(op, SYM_GreaterEqual)) register_inequality(c, a, b, -1, false);
        }
        return true;
    }

    /* Unequal[var, var] -> disequation store; other Unequal shapes and any
     * unrecognised conjunct fall back to full symbolic verification. */
    if (is_fun(e, SYM_Unequal, 2)) {
        Expr* a = e->data.function.args[0];
        Expr* b = e->data.function.args[1];
        if (a->type == EXPR_SYMBOL && b->type == EXPR_SYMBOL) {
            int ia = find_var_index(c, a->data.symbol.name);
            int ib = find_var_index(c, b->data.symbol.name);
            if (ia >= 0 && ib >= 0
                && c->n_neq < (int)(sizeof(c->neq_a)/sizeof(c->neq_a[0]))) {
                c->neq_a[c->n_neq] = ia; c->neq_b[c->n_neq] = ib; c->n_neq++;
                return true;
            }
        }
        c->all_captured = false;
        return true;
    }

    c->all_captured = false;     /* opaque conjunct: keep symbolic verification */
    return true;
}

/* ------------------------------------------------------------------ *
 *  Stage B: interval-positivity bounding.                             *
 * ------------------------------------------------------------------ */

/* Lower bound of an MPoly over the current box, EXCLUDING term index `skip`
 * (pass -1 to include all terms).  Assumes every variable that appears has
 * lo >= 0 (checked by the caller).  Writes the bound to `out` (pre-init'd)
 * and returns true when finite; false means a negative-coefficient term met
 * a variable with no finite upper bound. */
static bool poly_lower_bound_nonneg(const MPoly* p, const SICtx* c, int skip, mpz_t out) {
    mpz_set_ui(out, 0);
    mpz_t term, pw;
    mpz_init(term); mpz_init(pw);
    bool finite = true;
    for (size_t t = 0; t < p->n_terms && finite; t++) {
        if ((int)t == skip) continue;
        int sign = mpz_sgn(p->coefs[t]);
        mpz_set(term, p->coefs[t]);                /* term := coef */
        const int* ex = p->exps + t * (size_t)p->n_vars;
        for (int v = 0; v < p->n_vars && finite; v++) {
            int e = ex[v];
            if (e == 0) continue;
            /* Monotone reasoning is valid only for non-negative variables. */
            if (!c->has_lo[v] || c->lo[v] < 0) { finite = false; break; }
            /* coef > 0 -> minimise the monomial -> use lo^e;
             * coef < 0 -> maximise the monomial -> use hi^e (needs finite hi). */
            if (sign > 0) {
                mpz_set_si(pw, c->lo[v]);
            } else {
                if (!c->has_hi[v]) { finite = false; break; }
                mpz_set_si(pw, c->hi[v]);
            }
            for (int k = 0; k < e; k++) mpz_mul(term, term, pw);   /* term *= pw^e */
        }
        if (finite) mpz_add(out, out, term);
    }
    mpz_clear(term); mpz_clear(pw);
    return finite;
}

/* Try to tighten upper bounds using the relation G (interpreted as G <= 0).
 * For every positive-definite term that contains a variable v, that term is
 * <= -(lower bound of the rest), which bounds v.  Returns true if any hi
 * changed. */
static bool try_bound_up(SICtx* c, const MPoly* G) {
    bool changed = false;
    mpz_t rest_lo, urest, denom, q, root;
    mpz_init(rest_lo); mpz_init(urest); mpz_init(denom); mpz_init(q); mpz_init(root);

    for (size_t t = 0; t < G->n_terms; t++) {
        if (mpz_sgn(G->coefs[t]) <= 0) continue;                 /* need coef > 0 */
        const int* ex = G->exps + t * (size_t)G->n_vars;
        /* positive-definite: every variable in the term has lo >= 0. */
        bool posdef = true;
        for (int v = 0; v < G->n_vars; v++)
            if (ex[v] > 0 && !(c->has_lo[v] && c->lo[v] >= 0)) { posdef = false; break; }
        if (!posdef) continue;

        for (int v = 0; v < G->n_vars; v++) {
            int ev = ex[v];
            if (ev <= 0) continue;
            /* other vars in the term must have lo >= 1 (else the term can be 0). */
            bool ok = true;
            mpz_set_ui(denom, 1);
            for (int u = 0; u < G->n_vars && ok; u++) {
                if (u == v || ex[u] == 0) continue;
                if (c->lo[u] < 1) { ok = false; break; }
                for (int k = 0; k < ex[u]; k++) mpz_mul_si(denom, denom, c->lo[u]);
            }
            if (!ok) continue;
            mpz_mul(denom, denom, G->coefs[t]);                  /* coef * prod lo^e */
            if (mpz_sgn(denom) <= 0) continue;

            if (!poly_lower_bound_nonneg(G, c, (int)t, rest_lo)) continue;
            mpz_neg(urest, rest_lo);                             /* term <= urest */
            /* urest == 0 is a sound, tight bound (a positive-definite term
             * pinned to <= 0 forces its variable to 0); only a negative
             * allowance -- an already-infeasible split -- is skipped here. */
            if (mpz_sgn(urest) < 0) continue;

            mpz_fdiv_q(q, urest, denom);                         /* v^ev <= q */
            if (mpz_sgn(q) < 0) continue;
            mpz_root(root, q, (unsigned long)ev);                /* floor(q^(1/ev)) */
            if (mpz_fits_slong_p(root)) {
                int64_t bnd = mpz_get_si(root);
                if (!c->has_hi[v] || bnd < c->hi[v]) { c->hi[v] = bnd; c->has_hi[v] = true; changed = true; }
            }
        }
    }
    mpz_clear(rest_lo); mpz_clear(urest); mpz_clear(denom); mpz_clear(q); mpz_clear(root);
    return changed;
}

/* Upper bound of an MPoly over the current box, EXCLUDING term `skip`.
 * Symmetric to poly_lower_bound_nonneg; assumes every variable that appears
 * has lo >= 0.  Returns true when finite. */
static bool poly_upper_bound_nonneg(const MPoly* p, const SICtx* c, int skip, mpz_t out) {
    mpz_set_ui(out, 0);
    mpz_t term, pw;
    mpz_init(term); mpz_init(pw);
    bool finite = true;
    for (size_t t = 0; t < p->n_terms && finite; t++) {
        if ((int)t == skip) continue;
        int sign = mpz_sgn(p->coefs[t]);
        mpz_set(term, p->coefs[t]);
        const int* ex = p->exps + t * (size_t)p->n_vars;
        for (int v = 0; v < p->n_vars && finite; v++) {
            int e = ex[v];
            if (e == 0) continue;
            if (!c->has_lo[v] || c->lo[v] < 0) { finite = false; break; }
            /* coef > 0 -> maximise -> hi^e (needs finite hi); coef < 0 -> lo^e. */
            if (sign > 0) {
                if (!c->has_hi[v]) { finite = false; break; }
                mpz_set_si(pw, c->hi[v]);
            } else {
                mpz_set_si(pw, c->lo[v]);
            }
            for (int k = 0; k < e; k++) mpz_mul(term, term, pw);
        }
        if (finite) mpz_add(out, out, term);
    }
    mpz_clear(term); mpz_clear(pw);
    return finite;
}

/* Try to tighten LOWER bounds using G (as G <= 0 / G == 0): a positive-
 * definite term T that contains v satisfies T >= -(upper bound of the rest),
 * which -- when that is positive -- forces v up from below.  Returns true if
 * any lo changed. */
static bool try_bound_lo(SICtx* c, const MPoly* G) {
    bool changed = false;
    mpz_t rest_hi, brest, denom, q, root;
    mpz_init(rest_hi); mpz_init(brest); mpz_init(denom); mpz_init(q); mpz_init(root);

    for (size_t t = 0; t < G->n_terms; t++) {
        if (mpz_sgn(G->coefs[t]) <= 0) continue;
        const int* ex = G->exps + t * (size_t)G->n_vars;

        for (int v = 0; v < G->n_vars; v++) {
            int ev = ex[v];
            if (ev <= 0) continue;
            /* Sign feasibility: every OTHER variable in the term must be
             * non-negative (so the monomial's sign is fixed); the target v may
             * lack a prior lower bound only when ev is ODD, in which case a
             * strictly-positive term forces v > 0 (v^ev > 0 <=> v > 0). */
            bool okv = true;
            if (ev % 2 == 0 && !(c->has_lo[v] && c->lo[v] >= 0)) okv = false;
            for (int u = 0; u < G->n_vars && okv; u++) {
                if (u == v || ex[u] == 0) continue;
                if (!(c->has_lo[u] && c->lo[u] >= 0)) { okv = false; break; }
            }
            if (!okv) continue;
            /* other vars in the term contribute at most prod hi^e. */
            bool ok = true;
            mpz_set_ui(denom, 1);
            for (int u = 0; u < G->n_vars && ok; u++) {
                if (u == v || ex[u] == 0) continue;
                if (!c->has_hi[u]) { ok = false; break; }
                for (int k = 0; k < ex[u]; k++) mpz_mul_si(denom, denom, c->hi[u]);
            }
            if (!ok) continue;
            mpz_mul(denom, denom, G->coefs[t]);             /* coef * prod hi^e */
            if (mpz_sgn(denom) <= 0) continue;

            if (!poly_upper_bound_nonneg(G, c, (int)t, rest_hi)) continue;
            mpz_neg(brest, rest_hi);                        /* T >= brest */
            if (mpz_sgn(brest) <= 0) continue;              /* no lower info */

            mpz_cdiv_q(q, brest, denom);                    /* v^ev >= q (ceil) */
            if (mpz_sgn(q) <= 0) continue;
            /* smallest integer v with v^ev >= q is ceil(q^(1/ev)). */
            mpz_root(root, q, (unsigned long)ev);           /* floor root */
            /* bump to ceil if not exact */
            mpz_t chk; mpz_init(chk);
            mpz_pow_ui(chk, root, (unsigned long)ev);
            if (mpz_cmp(chk, q) < 0) mpz_add_ui(root, root, 1);
            mpz_clear(chk);
            if (mpz_fits_slong_p(root)) {
                int64_t bnd = mpz_get_si(root);
                if (!c->has_lo[v] || bnd > c->lo[v]) { c->lo[v] = bnd; c->has_lo[v] = true; changed = true; }
            }
        }
    }
    mpz_clear(rest_hi); mpz_clear(brest); mpz_clear(denom); mpz_clear(q); mpz_clear(root);
    return changed;
}

/* Propagate ordering constraints once; returns true if anything changed. */
static bool propagate_orderings(SICtx* c) {
    bool changed = false;
    for (int k = 0; k < c->n_ord; k++) {
        int a = c->ord_a[k], b = c->ord_b[k];      /* var a < / <= var b */
        int s = c->ord_strict[k] ? 1 : 0;
        if (c->has_hi[b]) {                          /* a <= hi_b - s */
            int64_t nb = c->hi[b] - s;
            if (!c->has_hi[a] || nb < c->hi[a]) { c->hi[a] = nb; c->has_hi[a] = true; changed = true; }
        }
        if (c->has_lo[a]) {                          /* b >= lo_a + s */
            int64_t nb = c->lo[a] + s;
            if (!c->has_lo[b] || nb > c->lo[b]) { c->lo[b] = nb; c->has_lo[b] = true; changed = true; }
        }
    }
    return changed;
}

/* Propagate abs-orderings once; returns true if anything changed.  For each
 * |a| < / <= |b|, if b is two-sided bounded then |b| <= Bb = max(|lo_b|,|hi_b|)
 * and |a| <= Bb - s (s = 1 for strict), giving a a symmetric box [-(Bb-s), Bb-s].
 * Only the smaller side is tightened -- |a| < |b| gives no usable *upper* bound
 * on |b| -- so this is a sound (necessary-condition) narrowing. */
static bool propagate_abs_orderings(SICtx* c) {
    bool changed = false;
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k];
        int s = c->abs_ord_strict[k] ? 1 : 0;
        if (!(c->has_lo[b] && c->has_hi[b])) continue;         /* |b| unbounded */
        int64_t alo = c->lo[b] < 0 ? -c->lo[b] : c->lo[b];
        int64_t ahi = c->hi[b] < 0 ? -c->hi[b] : c->hi[b];
        int64_t Bb  = alo > ahi ? alo : ahi;                   /* |b| <= Bb */
        int64_t Ba  = Bb - s;                                   /* |a| <= Ba */
        if (Ba < 0) Ba = -1;                                    /* empties the box */
        if (!c->has_hi[a] ||  Ba < c->hi[a]) { c->hi[a] =  Ba; c->has_hi[a] = true; changed = true; }
        if (!c->has_lo[a] || -Ba > c->lo[a]) { c->lo[a] = -Ba; c->has_lo[a] = true; changed = true; }
    }
    return changed;
}

/* Run the bound fixpoint (explicit bounds + ordering propagation +
 * interval-positivity).  Bounds are only ever tightened, so the result is a
 * set of necessary conditions. */
static void derive_bounds(SICtx* c) {
    for (int iter = 0; iter < 4 * c->n + 8; iter++) {
        bool changed = false;
        if (propagate_orderings(c)) changed = true;
        if (propagate_abs_orderings(c)) changed = true;
        for (int b = 0; b < c->nbc; b++) {
            if (try_bound_up(c, c->bc[b].Q)) changed = true;
            if (c->bc[b].kind == SI_EQ) {
                MPoly* neg = mpoly_neg(c->bc[b].Q);
                if (try_bound_up(c, neg)) changed = true;
                if (try_bound_lo(c, neg)) changed = true;
                mpoly_free(neg);
                if (try_bound_lo(c, c->bc[b].Q)) changed = true;
            }
        }
        if (!changed) break;
    }
}

/* A variable that appears only with EVEN exponents in every equation is
 * sign-symmetric: the equations constrain its magnitude, not its sign.  When
 * such a variable carries no explicit sign/ordering constraint it is left
 * unbounded by derive_bounds (the monotone bounder needs lo >= 0), so an
 * unconstrained sum of even powers -- the textbook  x^2 + y^2 == N  -- is
 * wrongly declined and Solve falls through to the generic path, which
 * fabricates {} over the Integers.
 *
 * Fix: treat each such variable as non-negative for the (monotone) bounding
 * pass -- v^even is minimised at 0, so a 0 lower bound is a sound lower bound
 * on every monomial it appears in and never loosens another variable's derived
 * bound -- then widen the derived [0, B] to the true symmetric window [-B, B]
 * so the search still covers the negative branch.  Downstream stages only ever
 * assume lo >= 0 inside the bounder, which has already finished here; the leaf
 * solver and enumeration handle negative domains directly. */
static void derive_even_only_bounds(SICtx* c) {
    bool cand[SI_MAX_VARS];
    bool any = false;
    for (int v = 0; v < c->n; v++) {
        cand[v] = false;
        if (c->has_lo[v] || c->has_hi[v]) continue;     /* already bounded */
        bool in_ord = false;                             /* sign matters if ordered */
        for (int k = 0; k < c->n_ord; k++)
            if (c->ord_a[k] == v || c->ord_b[k] == v) { in_ord = true; break; }
        if (in_ord) continue;
        bool appears = false, even_only = true;
        for (int q = 0; q < c->neq && even_only; q++) {
            const MPoly* p = c->eq[q];
            for (size_t t = 0; t < p->n_terms; t++) {
                int e = p->exps[t * (size_t)p->n_vars + v];
                if (e != 0) appears = true;
                if (e & 1) { even_only = false; break; }
            }
        }
        if (appears && even_only) { cand[v] = true; any = true; }
    }
    if (!any) return;

    for (int v = 0; v < c->n; v++)
        if (cand[v]) { c->lo[v] = 0; c->has_lo[v] = true; }   /* shadow non-negative */

    derive_bounds(c);

    for (int v = 0; v < c->n; v++) {
        if (!cand[v]) continue;
        if (c->has_hi[v] && c->hi[v] >= 0) c->lo[v] = -c->hi[v];   /* symmetric window */
        else c->has_lo[v] = false;                                /* not boundable: revert */
    }
}

/* ------------------------------------------------------------------ *
 *  Exact univariate leaf solver.                                      *
 * ------------------------------------------------------------------ */

/* Horner evaluation of a_0..a_d at integer r into `out`. */
static void horner_eval(const mpz_t* a, int d, int64_t r, mpz_t out) {
    mpz_set(out, a[d]);
    for (int k = d - 1; k >= 0; k--) {
        mpz_mul_si(out, out, r);
        mpz_add(out, out, a[k]);
    }
}

/* Collect the integer roots of a_0..a_d that lie in [lo, hi] into `roots`
 * (int64), returning the count.  `roots` must hold at least `d` entries. */
static int univariate_roots(const mpz_t* a, int d, int64_t lo, int64_t hi,
                            int64_t* roots) {
    int nr = 0;
    mpz_t tmp, disc, s, num;
    mpz_init(tmp); mpz_init(disc); mpz_init(s); mpz_init(num);

    if (d == 1) {                                    /* a1 v + a0 = 0 */
        if (mpz_sgn(a[1]) != 0 && mpz_divisible_p(a[0], a[1])) {
            mpz_neg(tmp, a[0]); mpz_divexact(tmp, tmp, a[1]);
            if (mpz_fits_slong_p(tmp)) {
                int64_t r = mpz_get_si(tmp);
                if (r >= lo && r <= hi) roots[nr++] = r;
            }
        }
        goto done;
    }
    if (d == 2) {                                    /* a2 v^2 + a1 v + a0 = 0 */
        mpz_mul(disc, a[1], a[1]);
        mpz_mul(tmp, a[2], a[0]); mpz_mul_ui(tmp, tmp, 4);
        mpz_sub(disc, disc, tmp);                    /* b^2 - 4ac */
        if (mpz_sgn(disc) < 0) goto done;
        if (!mpz_perfect_square_p(disc)) goto done;
        mpz_sqrt(s, disc);
        for (int sign = -1; sign <= 1; sign += 2) {
            mpz_set(num, s); if (sign < 0) mpz_neg(num, num);
            mpz_sub(num, num, a[1]);                 /* -b +/- sqrt */
            mpz_mul_ui(tmp, a[2], 2);
            if (mpz_sgn(tmp) != 0 && mpz_divisible_p(num, tmp)) {
                mpz_divexact(num, num, tmp);
                if (mpz_fits_slong_p(num)) {
                    int64_t r = mpz_get_si(num);
                    if (r >= lo && r <= hi) {
                        bool dup = false;
                        for (int i = 0; i < nr; i++) if (roots[i] == r) dup = true;
                        if (!dup) roots[nr++] = r;
                    }
                }
            }
        }
        goto done;
    }

    /* General degree: pure power a_d v^d + a_0 == 0 fast path. */
    {
        bool pure = true;
        for (int k = 1; k < d; k++) if (mpz_sgn(a[k]) != 0) { pure = false; break; }
        if (pure && mpz_sgn(a[d]) != 0) {
            if (mpz_divisible_p(a[0], a[d])) {
                mpz_neg(tmp, a[0]); mpz_divexact(tmp, tmp, a[d]);   /* tmp = v^d */
                if (d % 2 == 1) {                                   /* one real root */
                    mpz_abs(num, tmp);
                    if (mpz_root(s, num, (unsigned long)d) && mpz_fits_slong_p(s)) {
                        int64_t r = mpz_get_si(s);
                        if (mpz_sgn(tmp) < 0) r = -r;
                        if (r >= lo && r <= hi) roots[nr++] = r;
                    }
                } else if (mpz_sgn(tmp) >= 0) {                     /* +/- root */
                    if (mpz_root(s, tmp, (unsigned long)d) && mpz_fits_slong_p(s)) {
                        int64_t r = mpz_get_si(s);
                        if (r >= lo && r <= hi) roots[nr++] = r;
                        if (r != 0 && -r >= lo && -r <= hi) roots[nr++] = -r;
                    }
                }
            }
            goto done;
        }
    }

    /* General degree fallback: rational-root theorem on the constant term. */
    {
        int m = 0;                                    /* multiplicity of v = 0 */
        while (m <= d && mpz_sgn(a[m]) == 0) m++;
        if (m > 0 && 0 >= lo && 0 <= hi) roots[nr++] = 0;   /* v = 0 is a root */
        if (m <= d) {
            /* nonzero integer roots divide a[m]. */
            mpz_t am; mpz_init(am); mpz_abs(am, a[m]);
            Expr* dl = divisors_ordinary(am);
            if (dl && dl->type == EXPR_FUNCTION) {
                for (size_t i = 0; i < dl->data.function.arg_count; i++) {
                    int64_t dv;
                    if (!expr_as_i64(dl->data.function.args[i], &dv)) continue;
                    for (int sign = 1; sign >= -1; sign -= 2) {
                        int64_t r = sign * dv;
                        if (r < lo || r > hi) continue;
                        horner_eval(a, d, r, tmp);
                        if (mpz_sgn(tmp) == 0) {
                            bool dup = false;
                            for (int j = 0; j < nr; j++) if (roots[j] == r) dup = true;
                            if (!dup) roots[nr++] = r;
                        }
                    }
                }
            }
            if (dl) expr_free(dl);
            mpz_clear(am);
        }
    }

done:
    mpz_clear(tmp); mpz_clear(disc); mpz_clear(s); mpz_clear(num);
    return nr;
}

/* ------------------------------------------------------------------ *
 *  Stage C: recursive elimination search.                             *
 * ------------------------------------------------------------------ */

#define SI_LEAF_MAXDEG 8

typedef struct {
    SICtx*  ctx;
    int     leaf;               /* leaf variable index */
    int     order[SI_MAX_VARS]; /* search-var indices, outer-first */
    int     n_search;
    int64_t val[SI_MAX_VARS];   /* current assignment */
    int64_t* sols;              /* flattened n per solution */
    int      nsol, cap;
    int64_t  visits;            /* leaf nodes visited */
    int64_t  max_visits;        /* runtime backstop */
    bool     overflow;

    /* Per-equation int64 coefficient cache (leaf-independent): the int64 forms
     * of eq[q]->coefs, built once so the hot leaf loop never touches GMP.
     * eqc[q] is NULL / eqc_ok[q] false when a coefficient does not fit int64. */
    int64_t* eqc[SI_MAX_VARS * 2];
    bool     eqc_ok[SI_MAX_VARS * 2];
    bool     eqc_built;

    /* Multi-leaf staged elimination (A2): `n_peel` determined variables, each
     * resolved from a single equation given the enumerated free variables.
     * peel[k] is the k-th determined var; peel_eq[k] its resolving equation. */
    bool     multileaf;
    int      peel[SI_MAX_VARS];
    int      peel_eq[SI_MAX_VARS];
    int      n_peel;
} SearchState;

static void emit_solution(SearchState* st, int64_t leafval) {
    SICtx* c = st->ctx;
    if (st->nsol == st->cap) {
        st->cap = st->cap ? st->cap * 2 : 32;
        st->sols = (int64_t*)realloc(st->sols, sizeof(int64_t) * (size_t)st->cap * (size_t)c->n);
    }
    int64_t* row = st->sols + (size_t)st->nsol * (size_t)c->n;
    for (int i = 0; i < c->n; i++) row[i] = (i == st->leaf) ? leafval : st->val[i];
    st->nsol++;
}

/* Verify a full leaf assignment against the original conjunction. */
static bool verify_candidate(SearchState* st, int64_t leafval) {
    SICtx* c = st->ctx;
    int64_t vals[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) vals[i] = (i == st->leaf) ? leafval : st->val[i];
    return si_verify(c, vals);
}

/* Evaluate equation `eq`'s univariate-in-leaf coefficients a_0..a_d at the
 * current search assignment `st->val` (leaf left symbolic).  Reads the
 * precomputed monomials directly -- no MPoly allocation, the search hot
 * path.  Returns the leaf degree, -1 for identically zero, -2 if the leaf
 * degree exceeds the solver's reach. */
static int eval_leaf_coeffs(SearchState* st, const MPoly* eq, int leaf, mpz_t* a, mpz_t term) {
    SICtx* c = st->ctx;
    int n = c->n;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_set_ui(a[k], 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int el = ex[leaf];
        if (el > SI_LEAF_MAXDEG) return -2;
        mpz_set(term, eq->coefs[t]);
        for (int u = 0; u < n; u++) {
            if (u == leaf) continue;
            int e = ex[u];
            for (int k = 0; k < e; k++) mpz_mul_si(term, term, (long)st->val[u]);
        }
        mpz_add(a[el], a[el], term);
    }
    for (int k = SI_LEAF_MAXDEG; k >= 0; k--) if (mpz_sgn(a[k]) != 0) return k;
    return -1;
}

/* Build the per-equation int64 coefficient cache so the hot leaf loop never
 * calls mpz_get_si / mpz_fits_slong_p per node.  Idempotent. */
static void si_build_leaf_cache(SearchState* st) {
    if (st->eqc_built) return;
    SICtx* c = st->ctx;
    for (int q = 0; q < c->neq; q++) {
        const MPoly* eq = c->eq[q];
        st->eqc[q] = (int64_t*)malloc(sizeof(int64_t) * (eq->n_terms ? eq->n_terms : 1));
        bool ok = true;
        for (size_t t = 0; t < eq->n_terms; t++) {
            if (mpz_fits_slong_p(eq->coefs[t])) st->eqc[q][t] = mpz_get_si(eq->coefs[t]);
            else { ok = false; break; }
        }
        st->eqc_ok[q] = ok;
    }
    st->eqc_built = true;
}

static void si_free_leaf_cache(SearchState* st) {
    if (!st->eqc_built) return;
    for (int q = 0; q < st->ctx->neq; q++) { free(st->eqc[q]); st->eqc[q] = NULL; }
    st->eqc_built = false;
}

/* Exact floor(sqrt(n)) for n >= 0, division-based so it never overflows. */
static int64_t si_isqrt_i64(int64_t n) {
    if (n < 0) return -1;
    if (n == 0) return 0;
    int64_t r = (int64_t)sqrt((double)n);
    if (r < 1) r = 1;
    while (r > n / r) r--;                 /* r*r > n  (safe: r > 0) */
    while ((r + 1) <= n / (r + 1)) r++;    /* (r+1)^2 <= n */
    return r;
}

/* Integer roots of a2 x^2 + a1 x + a0 (degree d in {1,2}) inside [lo,hi], in
 * pure int64.  Sets *overflow and returns 0 when an intermediate would exceed
 * int64 (caller then falls back to the GMP path).  Returns the root count. */
static int si_uniroots_i64(int64_t a2, int64_t a1, int64_t a0, int d,
                           int64_t lo, int64_t hi, int64_t* roots, bool* overflow) {
    int nr = 0;
    *overflow = false;
    if (d == 1) {                                  /* a1 x + a0 = 0 */
        if (a1 != 0 && a0 % a1 == 0) {
            int64_t r = -a0 / a1;
            if (r >= lo && r <= hi) roots[nr++] = r;
        }
        return nr;
    }
    /* d == 2: discriminant a1^2 - 4 a2 a0. */
    int64_t a1sq, fac, disc, den;
    if (ci_mul_i64(a1, a1, &a1sq) || ci_mul_i64(a2, a0, &fac)
        || ci_mul_i64(fac, 4, &fac) || ci_sub_i64(a1sq, fac, &disc)
        || ci_mul_i64(a2, 2, &den)) { *overflow = true; return 0; }
    if (disc < 0) return 0;
    int64_t s = si_isqrt_i64(disc);
    int64_t s2;
    if (ci_mul_i64(s, s, &s2)) { *overflow = true; return 0; }
    if (s2 != disc) return 0;                      /* not a perfect square */
    for (int sign = -1; sign <= 1; sign += 2) {
        int64_t num = -a1 + sign * s;              /* |a1|,|s| already fit; sum fits */
        if (den != 0 && num % den == 0) {
            int64_t r = num / den;
            if (r >= lo && r <= hi) {
                bool dup = false;
                for (int i = 0; i < nr; i++) if (roots[i] == r) dup = true;
                if (!dup) roots[nr++] = r;
            }
        }
    }
    return nr;
}

/* Evaluate equation `eq`'s univariate-in-leaf coefficients at st->val entirely
 * in int64 and return the leaf roots in [lo,hi].  Return codes mirror the GMP
 * path: -1 = eq is identically zero in the leaf (no constraint), -2 = a nonzero
 * constant (infeasible), >=0 = root count.  Sets *need_mpz (caller falls back to
 * the exact GMP path) when the leaf degree exceeds 2 or any product/coefficient
 * would exceed int64.  This is the hot path for large ordered boxes. */
static int si_leaf_roots_i64(SearchState* st, const MPoly* eq, const int64_t* ci,
                             bool ci_ok, int leaf, int64_t lo, int64_t hi,
                             int64_t* roots, bool* need_mpz) {
    SICtx* c = st->ctx;
    int n = c->n;
    int64_t a0 = 0, a1 = 0, a2 = 0;
    *need_mpz = false;
    if (!ci_ok) { *need_mpz = true; return 0; }        /* a coefficient overflows int64 */
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int el = ex[leaf];
        if (el > 2) { *need_mpz = true; return 0; }
        int64_t val = ci[t];                            /* precomputed int64 coefficient */
        for (int u = 0; u < n; u++) {
            if (u == leaf) continue;
            for (int k = 0; k < ex[u]; k++)
                if (ci_mul_i64(val, st->val[u], &val)) { *need_mpz = true; return 0; }
        }
        int64_t* slot = (el == 0) ? &a0 : (el == 1) ? &a1 : &a2;
        if (ci_add_i64(*slot, val, slot)) { *need_mpz = true; return 0; }
    }
    int d = (a2 != 0) ? 2 : (a1 != 0) ? 1 : 0;
    if (d == 0) return (a0 != 0) ? -2 : -1;
    bool ovf;
    int nr = si_uniroots_i64(a2, a1, a0, d, lo, hi, roots, &ovf);
    if (ovf) { *need_mpz = true; return 0; }
    return nr;
}

/* Solve the leaf variable given a full search-var assignment in st->val. */
static void solve_leaf(SearchState* st) {
    SICtx* c = st->ctx;
    int leaf = st->leaf;
    int64_t lo = c->lo[leaf], hi = c->hi[leaf];
    if (++st->visits > st->max_visits) { st->overflow = true; return; }

    /* int64 fast path: no GMP allocation per node (the hot path for large
     * ordered boxes).  Falls back to the exact GMP body below whenever a leaf
     * degree exceeds 2 or an intermediate would exceed int64. */
    {
        int64_t inter[512]; int ninter = 0; bool have = false, feasible = true, need_mpz = false;
        for (int q = 0; q < c->neq && feasible; q++) {
            int64_t roots[8];
            int nr = si_leaf_roots_i64(st, c->eq[q], st->eqc[q], st->eqc_ok[q],
                                       st->leaf, lo, hi, roots, &need_mpz);
            if (need_mpz) break;
            if (nr == -1) continue;                  /* no constraint on leaf */
            if (nr == -2) { feasible = false; break; } /* nonzero constant: infeasible */
            if (!have) {
                for (int i = 0; i < nr && ninter < 512; i++) inter[ninter++] = roots[i];
                have = true;
            } else {
                int64_t keep[512]; int nk = 0;
                for (int i = 0; i < ninter; i++)
                    for (int j = 0; j < nr; j++)
                        if (inter[i] == roots[j]) { keep[nk++] = inter[i]; break; }
                memcpy(inter, keep, sizeof(int64_t) * (size_t)nk);
                ninter = nk;
            }
            if (ninter == 0) feasible = false;
        }
        if (!need_mpz) {
            if (feasible) {
                if (have) {
                    for (int i = 0; i < ninter; i++)
                        if (verify_candidate(st, inter[i])) emit_solution(st, inter[i]);
                } else {
                    if (hi - lo > 10000000LL) st->overflow = true;
                    else for (int64_t r = lo; r <= hi; r++)
                        if (verify_candidate(st, r)) emit_solution(st, r);
                }
            }
            return;
        }
    }

    int64_t inter[512]; int ninter = 0; bool have = false;
    mpz_t a[SI_LEAF_MAXDEG + 1], term;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_init(a[k]);
    mpz_init(term);

    bool feasible = true;
    for (int q = 0; q < c->neq && feasible; q++) {
        int d = eval_leaf_coeffs(st, c->eq[q], st->leaf, a, term);
        if (d == -1) continue;                       /* 0 == 0: no constraint */
        if (d == -2) { st->overflow = true; feasible = false; break; }
        if (d == 0) {                                /* constant == 0 ? */
            if (mpz_sgn(a[0]) != 0) feasible = false;
            continue;
        }
        int64_t roots[SI_LEAF_MAXDEG * 2];
        int nr = univariate_roots(a, d, lo, hi, roots);
        if (!have) {
            for (int i = 0; i < nr && ninter < 512; i++) inter[ninter++] = roots[i];
            have = true;
        } else {
            int64_t keep[512]; int nk = 0;
            for (int i = 0; i < ninter; i++)
                for (int j = 0; j < nr; j++)
                    if (inter[i] == roots[j]) { keep[nk++] = inter[i]; break; }
            memcpy(inter, keep, sizeof(int64_t) * (size_t)nk);
            ninter = nk;
        }
        if (ninter == 0) { feasible = false; }
    }

    if (feasible && !st->overflow) {
        if (have) {
            for (int i = 0; i < ninter; i++)
                if (verify_candidate(st, inter[i])) emit_solution(st, inter[i]);
        } else {
            /* No equation constrains the leaf at this branch: every value in
             * its range satisfying the constraints is a solution.  Only safe
             * for a genuinely small range; otherwise abort. */
            if (hi - lo > 10000000LL) { st->overflow = true; }
            else for (int64_t r = lo; r <= hi; r++)
                if (verify_candidate(st, r)) emit_solution(st, r);
        }
    }
    mpz_clear(term);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_clear(a[k]);
}

/* Effective [lo,hi] for search var at position, tightened by orderings
 * against already-assigned variables. */
static void effective_bounds(SearchState* st, int depth, int vi, int64_t* elo, int64_t* ehi) {
    SICtx* c = st->ctx;
    *elo = c->lo[vi]; *ehi = c->hi[vi];
    (void)depth;
    for (int k = 0; k < c->n_ord; k++) {
        int a = c->ord_a[k], b = c->ord_b[k]; int s = c->ord_strict[k] ? 1 : 0;
        /* var a < / <= var b */
        if (a == vi) {                               /* vi < b: hi <= val_b - s */
            bool assigned = false; int64_t vb = 0;
            for (int p = 0; p < depth; p++) if (st->order[p] == b) { assigned = true; vb = st->val[b]; }
            if (assigned && vb - s < *ehi) *ehi = vb - s;
        }
        if (b == vi) {                               /* a < vi: lo >= val_a + s */
            bool assigned = false; int64_t va = 0;
            for (int p = 0; p < depth; p++) if (st->order[p] == a) { assigned = true; va = st->val[a]; }
            if (assigned && va + s > *elo) *elo = va + s;
        }
    }
    /* Abs-orderings prune only the smaller side: |vi| < |val_b| once the larger
     * variable b is assigned, tightening vi to a symmetric window.  (The larger
     * side is a hole around 0, not a contiguous narrowing, so it is left to the
     * final verify.)  Realised only when b is enumerated before vi -- see the
     * abs-aware search-order builder. */
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k]; int s = c->abs_ord_strict[k] ? 1 : 0;
        if (a != vi) continue;
        bool assigned = false; int64_t vb = 0;
        for (int p = 0; p < depth; p++) if (st->order[p] == b) { assigned = true; vb = st->val[b]; }
        if (!assigned) continue;
        int64_t avb = vb < 0 ? -vb : vb;
        int64_t bnd = avb - s;                       /* |vi| <= bnd */
        if (bnd < 0) bnd = -1;
        if ( bnd < *ehi) *ehi =  bnd;
        if (-bnd > *elo) *elo = -bnd;
    }
}

static void resolve_peeled(SearchState* st, int pi);   /* multi-leaf resolver */

static void search_rec(SearchState* st, int depth) {
    if (st->overflow) return;
    if (depth == st->n_search) {
        if (st->multileaf) resolve_peeled(st, 0); else solve_leaf(st);
        return;
    }

    int vi = st->order[depth];
    int64_t elo, ehi;
    effective_bounds(st, depth, vi, &elo, &ehi);

    for (int64_t r = elo; r <= ehi && !st->overflow; r++) {
        st->val[vi] = r;
        search_rec(st, depth + 1);
    }
}

/* Length of the longest <= / < ordering chain among the search vars `vars`.
 * A chain of length L cuts an otherwise-N^L box to ~N^L / L! ordered tuples,
 * so the search-space guard divides the raw estimate by L! -- otherwise an
 * ordered box like  0 < x <= y <= z < 1000  (raw 10^9, ordered ~1.7*10^8) is
 * wrongly declined.  The visit cap backstops any under-estimate, so this only
 * ever lets a genuinely-ordered box be attempted, never truncates a result. */
static int si_longest_chain(const SICtx* c, const int* vars, int nv) {
    bool in[SI_MAX_VARS];
    int memo[SI_MAX_VARS];
    for (int i = 0; i < SI_MAX_VARS; i++) { in[i] = false; memo[i] = 1; }
    for (int i = 0; i < nv; i++) in[vars[i]] = true;
    /* memo[v] = longest chain starting at v; relax to a fixpoint (DAG, tiny). */
    for (int iter = 0; iter < nv; iter++) {
        for (int i = 0; i < nv; i++) {
            int v = vars[i], m = 1;
            for (int e = 0; e < c->n_ord; e++) {
                if (c->ord_a[e] == v && in[c->ord_b[e]]) {   /* edge v < / <= b */
                    int cand = 1 + memo[c->ord_b[e]];
                    if (cand > m) m = cand;
                }
            }
            for (int e = 0; e < c->n_abs_ord; e++) {         /* |v| < / <= |b| */
                if (c->abs_ord_a[e] == v && in[c->abs_ord_b[e]]) {
                    int cand = 1 + memo[c->abs_ord_b[e]];
                    if (cand > m) m = cand;
                }
            }
            if (m > memo[v]) memo[v] = m;
        }
    }
    int best = 1;
    for (int i = 0; i < nv; i++) if (memo[vars[i]] > best) best = memo[vars[i]];
    return best;
}

/* ------------------------------------------------------------------ *
 *  Result assembly.                                                   *
 * ------------------------------------------------------------------ */

static int sol_cmp(const int64_t* a, const int64_t* b, int n) {
    for (int i = 0; i < n; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* qsort comparator over whole solution rows; row width via a file-scope shim
 * (the solve is single-threaded, so a static is safe and avoids qsort_r). */
static int g_build_row_n = 0;
static int build_row_cmp(const void* a, const void* b) {
    return sol_cmp((const int64_t*)a, (const int64_t*)b, g_build_row_n);
}

static Expr* build_result(SearchState* st) {
    SICtx* c = st->ctx;
    int ns = st->nsol, n = c->n;
    /* Sort the tuple set in place: O(ns log ns), so a large solution family
     * (e.g. a parametric slice enumerated over a wide box) does not blow up. */
    if (ns > 1) {
        g_build_row_n = n;
        qsort(st->sols, (size_t)ns, sizeof(int64_t) * (size_t)n, build_row_cmp);
    }
    Expr** tuples = (Expr**)malloc(sizeof(Expr*) * (size_t)(ns > 0 ? ns : 1));
    int nt = 0;
    for (int i = 0; i < ns; i++) {
        int64_t* row = st->sols + (size_t)i * n;
        if (i > 0 && sol_cmp(st->sols + (size_t)(i - 1) * n, row, n) == 0) continue; /* dedupe */
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int k = 0; k < n; k++)
            rules[k] = mk_rule(expr_copy(c->var[k]), mk_int(row[k]));
        tuples[nt++] = mk_list(rules, (size_t)n);
        free(rules);
    }
    Expr* res = mk_list(tuples, (size_t)nt);
    free(tuples);
    return res;
}

/* ------------------------------------------------------------------ *
 *  Meet-in-the-middle for a separable additive equation.              *
 *                                                                     *
 *  When the single equation is  sum_i g_i(x_i) == TARGET  with each    *
 *  g_i univariate (every term touches at most one variable), the       *
 *  ~N^(n-1) leaf search is replaced by splitting the variables into    *
 *  two groups, tabulating one group's partial sums and binary-         *
 *  searching the other -- ~N^ceil(n/2) work.  Correctness rests on the *
 *  same final verify_candidate against the original conjunction, so    *
 *  cross-group orderings / disequations need no special handling here. *
 * ------------------------------------------------------------------ */

typedef struct { int64_t sum; int64_t vals[SI_MAX_VARS]; } MitmEntry;

#define MITM_HASH_CAP 5000000LL

static int mitm_cmp(const void* pa, const void* pb) {
    int64_t a = ((const MitmEntry*)pa)->sum, b = ((const MitmEntry*)pb)->sum;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Emit / verify a complete assignment (all n variables set in `vals`). */
static bool verify_full(SICtx* c, const int64_t* vals) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)c->n);
    for (int i = 0; i < c->n; i++)
        rules[i] = mk_rule(expr_copy(c->var[i]), mk_int(vals[i]));
    Expr* rl = mk_list(rules, (size_t)c->n);
    free(rules);
    Expr* subbed = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(c->original), rl }, 2));
    bool ok = is_sym(subbed, SYM_True);
    expr_free(subbed);
    return ok;
}

/* Fast candidate check.  When every constraint was captured in the store
 * (bounds / orderings / disequations), a candidate is validated purely
 * numerically: every equation MPoly must vanish and every stored constraint
 * must hold -- no Expr allocation, no symbolic evaluation.  This is what keeps
 * the divisor / reciprocal paths (hundreds of candidates) fast.  Anything with
 * an un-captured constraint falls back to the full symbolic check. */
static bool si_verify(SICtx* c, const int64_t* vals) {
    if (!c->all_captured) return verify_full(c, vals);
    mpz_t r; mpz_init(r);
    for (int q = 0; q < c->neq; q++) {
        si_eval_mpoly(c->eq[q], vals, r);
        if (mpz_sgn(r) != 0) { mpz_clear(r); return false; }
    }
    mpz_clear(r);
    for (int i = 0; i < c->n; i++) {
        if (c->has_lo[i] && vals[i] < c->lo[i]) return false;
        if (c->has_hi[i] && vals[i] > c->hi[i]) return false;
    }
    for (int k = 0; k < c->n_ord; k++) {
        int a = c->ord_a[k], b = c->ord_b[k];
        if (c->ord_strict[k] ? !(vals[a] < vals[b]) : !(vals[a] <= vals[b])) return false;
    }
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k];
        int64_t va = vals[a] < 0 ? -vals[a] : vals[a];
        int64_t vb = vals[b] < 0 ? -vals[b] : vals[b];
        if (c->abs_ord_strict[k] ? !(va < vb) : !(va <= vb)) return false;
    }
    for (int k = 0; k < c->n_neq; k++)
        if (vals[c->neq_a[k]] == vals[c->neq_b[k]]) return false;
    return true;
}

static void emit_full(SearchState* st, const int64_t* vals) {
    SICtx* c = st->ctx;
    if (st->nsol == st->cap) {
        st->cap = st->cap ? st->cap * 2 : 32;
        st->sols = (int64_t*)realloc(st->sols, sizeof(int64_t) * (size_t)st->cap * (size_t)c->n);
    }
    int64_t* row = st->sols + (size_t)st->nsol * (size_t)c->n;
    for (int i = 0; i < c->n; i++) row[i] = vals[i];
    st->nsol++;
}

/* ------------------------------------------------------------------ *
 *  A2: multi-leaf staged elimination.                                 *
 *                                                                     *
 *  The single-leaf search enumerates every variable but one.  A system *
 *  like the Euler brick  x^2+y^2==a^2 && x^2+z^2==b^2 && y^2+z^2==c^2   *
 *  has THREE variables (a,b,c) each determined by a single equation --  *
 *  enumerating even two of them is hopeless.  Staged elimination peels  *
 *  every variable that appears in exactly one equation and is           *
 *  univariate-solvable there, resolving it per free-variable assignment *
 *  (an exact root, branching on multiplicity) instead of enumerating.   *
 *  Only the genuinely-coupled "free" variables are walked.              *
 * ------------------------------------------------------------------ */

/* Resolve peeled variable peel[pi] from its equation given the free-variable
 * assignment (and earlier-resolved peels) in st->val, recursing to the next
 * peel; at the end, verify the full assignment and emit.  Each peel_eq contains
 * only its own peeled variable plus free variables, so any resolution order is
 * valid. */
static void resolve_peeled(SearchState* st, int pi) {
    if (st->overflow) return;
    SICtx* c = st->ctx;
    if (pi == st->n_peel) {
        if (++st->visits > st->max_visits) { st->overflow = true; return; }
        if (si_verify(c, st->val)) emit_full(st, st->val);
        return;
    }
    int v = st->peel[pi], e = st->peel_eq[pi];
    int64_t lo = c->has_lo[v] ? c->lo[v] : -(1LL << 50);
    int64_t hi = c->has_hi[v] ? c->hi[v] :  (1LL << 50);

    /* int64 fast path (peeled var degree <= 2). */
    int64_t roots[8]; bool need_mpz = false;
    int nr = si_leaf_roots_i64(st, c->eq[e], st->eqc[e], st->eqc_ok[e], v, lo, hi, roots, &need_mpz);
    if (!need_mpz) {
        if (nr == -2) return;                       /* infeasible at this branch */
        if (nr == -1) { st->overflow = true; return; } /* v unconstrained here: decline */
        for (int i = 0; i < nr && !st->overflow; i++) { st->val[v] = roots[i]; resolve_peeled(st, pi + 1); }
        return;
    }

    /* exact GMP fallback (peeled var degree up to SI_LEAF_MAXDEG). */
    mpz_t a[SI_LEAF_MAXDEG + 1], term;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_init(a[k]);
    mpz_init(term);
    int d = eval_leaf_coeffs(st, c->eq[e], v, a, term);
    if (d == -2) st->overflow = true;               /* degree beyond the solver's reach */
    else if (d == -1) st->overflow = true;          /* v unconstrained here: decline */
    else if (d == 0) { /* constant: feasible iff 0, but then v is free -> decline */
        if (mpz_sgn(a[0]) != 0) { /* infeasible: no roots */ } else st->overflow = true;
    } else {
        int64_t r2[SI_LEAF_MAXDEG * 2];
        int nr2 = univariate_roots(a, d, lo, hi, r2);
        for (int i = 0; i < nr2 && !st->overflow; i++) { st->val[v] = r2[i]; resolve_peeled(st, pi + 1); }
    }
    mpz_clear(term);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_clear(a[k]);
}

/* Detect a staged-elimination structure and, if the free variables form a
 * tractable ordered box, run the multi-leaf search.  Returns true when it
 * settled the answer (candidates in st->sols); false to fall back. */
static bool si_solve_multileaf(SICtx* c, SearchState* st) {
    int n = c->n;
    if (c->neq < 2) return false;                   /* single-leaf path handles neq == 1 */

    /* appear[v] = number of equations that contain v (deg > 0). */
    int appear[SI_MAX_VARS], only_eq[SI_MAX_VARS];
    for (int v = 0; v < n; v++) { appear[v] = 0; only_eq[v] = -1; }
    for (int q = 0; q < c->neq; q++)
        for (int v = 0; v < n; v++)
            if (mpoly_deg_var(c->eq[q], v) > 0) { appear[v]++; only_eq[v] = q; }

    /* A variable is a peel candidate iff it appears in exactly one equation and
     * is univariate-solvable there.  Each equation is claimed by at most one
     * peeled variable (widest domain / unbounded first); losers become free. */
    bool peeled[SI_MAX_VARS]; int claim[SI_MAX_VARS * 2];
    for (int v = 0; v < n; v++) peeled[v] = false;
    for (int q = 0; q < c->neq; q++) claim[q] = -1;
    for (int v = 0; v < n; v++) {
        if (appear[v] != 1) continue;
        int e = only_eq[v];
        int dg = mpoly_deg_var(c->eq[e], v);
        if (dg < 1 || dg > SI_LEAF_MAXDEG) continue;
        int cur = claim[e];
        if (cur < 0) { claim[e] = v; continue; }
        /* Prefer to peel the wider / unbounded variable (enumerate the narrow). */
        long double wv = (c->has_lo[v] && c->has_hi[v]) ? (long double)(c->hi[v] - c->lo[v]) : 1e30L;
        long double wc = (c->has_lo[cur] && c->has_hi[cur]) ? (long double)(c->hi[cur] - c->lo[cur]) : 1e30L;
        if (wv > wc) claim[e] = v;
    }
    st->n_peel = 0;
    for (int q = 0; q < c->neq; q++)
        if (claim[q] >= 0) { int v = claim[q]; peeled[v] = true; st->peel[st->n_peel] = v; st->peel_eq[st->n_peel] = q; st->n_peel++; }
    if (st->n_peel < 2) return false;               /* single leaf: ordinary path */

    /* Free variables: everything not peeled.  All must be finitely bounded. */
    int free_v[SI_MAX_VARS], nfree = 0;
    for (int v = 0; v < n; v++) if (!peeled[v]) {
        if (!(c->has_lo[v] && c->has_hi[v])) return false;
        free_v[nfree++] = v;
    }
    if (nfree == 0) return false;

    /* Ordered free-box estimate (raw product / longest-chain factorial). */
    long double est = 1.0L;
    for (int i = 0; i < nfree; i++) est *= (long double)(c->hi[free_v[i]] - c->lo[free_v[i]] + 1);
    int L = si_longest_chain(c, free_v, nfree);
    for (int i = 2; i <= L; i++) est /= (long double)i;
    if (est > (long double)SI_MAX_NODES) return false;

    /* Order free vars ascending domain (narrowest innermost). */
    st->n_search = 0;
    for (int i = 0; i < nfree; i++) st->order[st->n_search++] = free_v[i];
    for (int i = 0; i + 1 < st->n_search; i++) {
        int pk = i;
        for (int j = i + 1; j < st->n_search; j++)
            if ((c->hi[st->order[j]] - c->lo[st->order[j]]) < (c->hi[st->order[pk]] - c->lo[st->order[pk]])) pk = j;
        if (pk != i) { int t = st->order[i]; st->order[i] = st->order[pk]; st->order[pk] = t; }
    }
    st->multileaf = true;
    st->max_visits = SI_MAX_NODES;
    si_build_leaf_cache(st);
    for (int i = 0; i < n; i++) st->val[i] = 0;
    search_rec(st, 0);
    si_free_leaf_cache(st);
    return !st->overflow;
}

/* Attempt the meet-in-the-middle path.  Returns true if it ran (results are
 * in st->sols); false to fall back to the general search. */
static bool mitm_solve(SearchState* st) {
    SICtx* c = st->ctx;
    int n = c->n;
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0];

    /* Separability + constant term. */
    mpz_t c0; mpz_init(c0); mpz_set_ui(c0, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nv = 0;
        for (int v = 0; v < n; v++) if (ex[v] > 0) nv++;
        if (nv >= 2) { mpz_clear(c0); return false; }   /* not separable */
        if (nv == 0) mpz_add(c0, c0, eq->coefs[t]);
    }

    /* Fully bounded? */
    for (int i = 0; i < n; i++)
        if (!(c->has_lo[i] && c->has_hi[i]) || c->hi[i] < c->lo[i]) { mpz_clear(c0); return false; }

    /* TARGET = -c0 as int64. */
    mpz_neg(c0, c0);
    if (!mpz_fits_slong_p(c0)) { mpz_clear(c0); return false; }
    int64_t target = mpz_get_si(c0);
    mpz_clear(c0);

    /* Per-variable value tables g_i(v). */
    int64_t* gtab[SI_MAX_VARS];
    for (int i = 0; i < n; i++) gtab[i] = NULL;
    int64_t domain[SI_MAX_VARS];
    long double total_tab = 0.0L;
    bool ok = true;
    mpz_t gv, pw;
    mpz_init(gv); mpz_init(pw);
    for (int i = 0; i < n && ok; i++) {
        domain[i] = c->hi[i] - c->lo[i] + 1;
        total_tab += (long double)domain[i];
        if (total_tab > 100000000.0L) { ok = false; break; }   /* tables too big */
        gtab[i] = (int64_t*)malloc(sizeof(int64_t) * (size_t)domain[i]);
        for (int64_t d = 0; d < domain[i] && ok; d++) {
            int64_t v = c->lo[i] + d;
            mpz_set_ui(gv, 0);
            for (size_t t = 0; t < eq->n_terms; t++) {
                const int* ex = eq->exps + t * (size_t)n;
                if (ex[i] <= 0) continue;
                bool only_i = true;
                for (int u = 0; u < n; u++) if (u != i && ex[u] > 0) { only_i = false; break; }
                if (!only_i) continue;
                mpz_set(pw, eq->coefs[t]);
                for (int k = 0; k < ex[i]; k++) mpz_mul_si(pw, pw, (long)v);
                mpz_add(gv, gv, pw);
            }
            if (!mpz_fits_slong_p(gv)) { ok = false; break; }
            gtab[i][d] = mpz_get_si(gv);
        }
    }
    mpz_clear(gv); mpz_clear(pw);
    if (!ok) { for (int i = 0; i < n; i++) free(gtab[i]); return false; }

    /* Balance variables into two groups by product of domains (greedy). */
    int order[SI_MAX_VARS];
    for (int i = 0; i < n; i++) order[i] = i;
    for (int i = 0; i + 1 < n; i++) {           /* sort desc by domain */
        int pk = i;
        for (int j = i + 1; j < n; j++) if (domain[order[j]] > domain[order[pk]]) pk = j;
        if (pk != i) { int tmp = order[i]; order[i] = order[pk]; order[pk] = tmp; }
    }
    int gA[SI_MAX_VARS], gB[SI_MAX_VARS], nA = 0, nB = 0;
    long double pA = 1.0L, pB = 1.0L;
    for (int k = 0; k < n; k++) {
        int i = order[k];
        if (pA <= pB) { gA[nA++] = i; pA *= (long double)domain[i]; }
        else          { gB[nB++] = i; pB *= (long double)domain[i]; }
    }
    /* Hash the smaller group, iterate the larger. */
    int *hg, *ig, nh, ni;
    if (pA <= pB) { hg = gA; nh = nA; ig = gB; ni = nB; }
    else          { hg = gB; nh = nB; ig = gA; ni = nA; }
    long double p_hash = (pA <= pB) ? pA : pB;
    long double p_iter = (pA <= pB) ? pB : pA;
    if (p_hash > (long double)MITM_HASH_CAP || p_iter > (long double)SI_MAX_NODES) {
        for (int i = 0; i < n; i++) free(gtab[i]);
        return false;
    }

    /* Build the hash-group table via an odometer. */
    size_t hcap = (size_t)(p_hash) + 1, hcnt = 0;
    MitmEntry* H = (MitmEntry*)malloc(sizeof(MitmEntry) * hcap);
    int idx[SI_MAX_VARS];
    for (int j = 0; j < nh; j++) idx[j] = 0;
    for (;;) {
        int64_t sum = 0;
        MitmEntry* e = &H[hcnt];
        for (int j = 0; j < n; j++) e->vals[j] = 0;
        for (int j = 0; j < nh; j++) {
            int vi = hg[j]; int64_t v = c->lo[vi] + idx[j];
            e->vals[vi] = v; sum += gtab[vi][idx[j]];
        }
        e->sum = sum; hcnt++;
        int j = 0;                               /* odometer increment */
        for (; j < nh; j++) { if (++idx[j] < domain[hg[j]]) break; idx[j] = 0; }
        if (j == nh) break;
    }
    qsort(H, hcnt, sizeof(MitmEntry), mitm_cmp);

    /* Iterate the other group; binary-search complements. */
    int64_t full[SI_MAX_VARS];
    for (int j = 0; j < ni; j++) idx[j] = 0;
    if (ni == 0) {
        /* No iterate group: TARGET must be hit by hash entries alone. */
        for (size_t h = 0; h < hcnt; h++)
            if (H[h].sum == target && si_verify(c, H[h].vals)) emit_full(st, H[h].vals);
    } else for (;;) {
        int64_t sum = 0;
        for (int i = 0; i < n; i++) full[i] = 0;
        for (int j = 0; j < ni; j++) {
            int vi = ig[j]; int64_t v = c->lo[vi] + idx[j];
            full[vi] = v; sum += gtab[vi][idx[j]];
        }
        int64_t need = target - sum;
        /* lower_bound on need */
        size_t lo = 0, hi = hcnt;
        while (lo < hi) { size_t mid = (lo + hi) / 2; if (H[mid].sum < need) lo = mid + 1; else hi = mid; }
        for (size_t h = lo; h < hcnt && H[h].sum == need; h++) {
            for (int j = 0; j < nh; j++) full[hg[j]] = H[h].vals[hg[j]];
            if (si_verify(c, full)) emit_full(st, full);
        }
        int j = 0;
        for (; j < ni; j++) { if (++idx[j] < domain[ig[j]]) break; idx[j] = 0; }
        if (j == ni) break;
    }

    free(H);
    for (int i = 0; i < n; i++) free(gtab[i]);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Top-level.                                                         *
 * ------------------------------------------------------------------ */

static void ctx_free(SICtx* c) {
    for (int i = 0; i < c->neq; i++) mpoly_free(c->eq[i]);
    for (int i = 0; i < c->nbc; i++) mpoly_free(c->bc[i].Q);
}

/* Flatten a top-level And / List into a conjunct array (borrowed).  A bare
 * relation is a single conjunct. */
static void flatten_conjuncts(Expr* e, Expr*** out, int* n) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol.name == SYM_And
            || e->data.function.head->data.symbol.name == SYM_List)) {
        *n = (int)e->data.function.arg_count;
        *out = e->data.function.args;
    } else {
        static Expr* one[1];
        one[0] = e;
        *out = one;
        *n = 1;
    }
}

/* --- Solve::svars diagnostic --------------------------------------------- *
 * Collect the distinct symbol atoms of `e` (recursing into function arguments
 * but never the head, so operator heads such as Plus / Greater / Inequality
 * are skipped) into seen[].  Stores interned name pointers; deduplicated. */
static void si_collect_atoms(const Expr* e, const char** seen, int* nseen, int cap) {
    if (!e || *nseen >= cap) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        for (int i = 0; i < *nseen; i++) if (seen[i] == nm) return;
        seen[(*nseen)++] = nm;
        return;
    }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            si_collect_atoms(e->data.function.args[i], seen, nseen, cap);
}

/* Warn (once) when the system mentions a non-constant symbol that is not among
 * the solve variables -- almost always a mistyped variable list, e.g. the free
 * `d` left out of {x, y, z, y}.  Mirrors Mathematica's Solve::svars.  Purely
 * advisory: it does not change the result, and the caller proceeds unchanged.
 * Operator heads and named constants (Pi, E, ...) are Protected and skipped, so
 * only genuine variable-like symbols trip it. */
static void si_warn_free_symbols(const SICtx* c, Expr** conj, int ncj) {
    enum { CAP = 64 };
    const char* seen[CAP]; int nseen = 0;
    /* Scan only inequality / ordering constraints, not equations: a free symbol
     * carrying its OWN bound (`d > 0 && d < 100000`) is a near-certain dropped
     * variable, whereas a bare parameter in an equation (`a x == b`) is a
     * legitimate symbolic solve and must not be flagged. */
    for (int i = 0; i < ncj; i++) {
        Expr* e = conj[i];
        if (is_fun(e, SYM_Less, 2) || is_fun(e, SYM_LessEqual, 2)
         || is_fun(e, SYM_Greater, 2) || is_fun(e, SYM_GreaterEqual, 2)
         || is_fun(e, SYM_Unequal, 2)
         || (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
             && e->data.function.head->data.symbol.name == SYM_Inequality))
            si_collect_atoms(e, seen, &nseen, CAP);
    }
    for (int i = 0; i < nseen; i++) {
        const char* nm = seen[i];
        bool is_var = false;
        for (int v = 0; v < c->n; v++)
            if (c->var[v]->data.symbol.name == nm) { is_var = true; break; }
        if (is_var) continue;
        if (get_attributes(nm) & ATTR_PROTECTED) continue;   /* head / constant */
        /* The evaluator re-enters a NULL-returning builtin once more to confirm
         * the fixed point, which would print the message twice.  Suppress the
         * immediate repeat by remembering the last system we warned about. */
        static uint64_t last_hash = 0; static bool have_last = false;
        uint64_t h = expr_hash(c->original);
        if (have_last && h == last_hash) return;
        last_hash = h; have_last = true;
        fprintf(stderr, "Solve::svars: Equations may not give solutions for "
                        "all \"solve\" variables.\n");
        return;                                              /* emit only once */
    }
}

/* ------------------------------------------------------------------ *
 *  Phase 2 special forms: divisor-factoring / reciprocal solvers.     *
 *                                                                     *
 *  Two shapes that the bounded search cannot reach because positivity *
 *  fails to bound them, but which are finite once the right identity   *
 *  is applied:                                                         *
 *                                                                     *
 *  - A single bilinear equation  a*u*v + b*u + c*v + d == 0  (after    *
 *    eliminating unit-coefficient linear equations) factors as        *
 *    (a*u + c)(a*v + b) = b*c - a*d, so the integer solutions come     *
 *    from the DIVISORS of that constant -- no enumeration of u, v.     *
 *    This is the Pythagorean-perimeter case once z is eliminated.      *
 *                                                                     *
 *  - A sum of unit fractions  sum 1/x_i == R  with an ordering chain   *
 *    x_1 <= ... <= x_k bounds the smallest variable to                *
 *    [ceil(1/R), floor(k/R)] and recurses; the last variable is        *
 *    determined exactly.  This is the Egyptian-fraction case.          *
 * ------------------------------------------------------------------ */

/* Convert a residual Expr to an MPoly, clearing denominators and expanding
 * first (so composed forms like (3000 - x - y)^2 become flat monomials). */
static MPoly* si_resid_to_mpoly(Expr* resid, Expr** vars, int n) {
    Expr* tog = eval_and_free(internal_together((Expr*[]){ expr_copy(resid) }, 1));
    Expr* num = eval_and_free(internal_numerator((Expr*[]){ tog }, 1));
    Expr* exp = eval_and_free(internal_expand((Expr*[]){ num }, 1));
    MPoly* P = expr_to_mpoly(exp, vars, n);
    expr_free(exp);
    return P;
}

/* --- Reciprocal (Egyptian-fraction) recursion. --- */

static void si_recip_rec(SICtx* c, SearchState* st, const int* order, int k,
                         int pos, const mpz_t p, const mpz_t q, int64_t prev,
                         int64_t* full) {
    if (st->overflow) return;
    if (mpz_sgn(p) <= 0) return;                 /* R must stay positive */
    int m = k - pos;                             /* remaining terms */
    int vi = order[pos];

    if (m == 1) {                                /* 1/v == p/q -> v = q/p */
        if (mpz_divisible_p(q, p)) {
            mpz_t vv; mpz_init(vv); mpz_divexact(vv, q, p);
            if (mpz_fits_slong_p(vv)) {
                int64_t v = mpz_get_si(vv);
                if (v >= prev) { full[vi] = v; if (si_verify(c, full)) emit_full(st, full); }
            }
            mpz_clear(vv);
        }
        return;
    }

    /* 1/v < p/q  => v >= floor(q/p)+1;  1/v >= (p/q)/m => v <= floor(m q / p). */
    mpz_t lo_m, hi_m, tmp; mpz_init(lo_m); mpz_init(hi_m); mpz_init(tmp);
    mpz_fdiv_q(lo_m, q, p); mpz_add_ui(lo_m, lo_m, 1);
    mpz_mul_ui(tmp, q, (unsigned long)m); mpz_fdiv_q(hi_m, tmp, p);
    int64_t lo = mpz_fits_slong_p(lo_m) ? mpz_get_si(lo_m) : INT64_MAX;
    int64_t hi = mpz_fits_slong_p(hi_m) ? mpz_get_si(hi_m) : INT64_MAX;
    mpz_clear(lo_m); mpz_clear(hi_m); mpz_clear(tmp);
    if (prev > lo) lo = prev;

    for (int64_t v = lo; v <= hi && !st->overflow; v++) {
        if (++st->visits > st->max_visits) { st->overflow = true; return; }
        mpz_t np, nq, g; mpz_init(np); mpz_init(nq); mpz_init(g);
        mpz_mul_si(np, p, (long)v); mpz_sub(np, np, q);   /* p*v - q */
        mpz_mul_si(nq, q, (long)v);                        /* q*v */
        if (mpz_sgn(np) > 0) {
            mpz_gcd(g, np, nq); mpz_divexact(np, np, g); mpz_divexact(nq, nq, g);
            full[vi] = v;
            si_recip_rec(c, st, order, k, pos + 1, np, nq, v, full);
        }
        mpz_clear(np); mpz_clear(nq); mpz_clear(g);
    }
}

/* Detect the shape  c_full * prod x_i - a * sum_i prod_{j!=i} x_j == 0
 * (equivalently  a * sum 1/x_i == c_full), with all a_i equal.  Fills the
 * active-variable list, c_full and a. */
static bool si_reciprocal_detect(const MPoly* eq, int n, int* active, int* ka_out,
                                 mpz_t c_full, mpz_t a_out) {
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++) if (ex[v] > 1) return false;   /* multilinear only */
    }
    int ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) active[ka++] = v;
    if (ka < 2 || eq->n_terms != (size_t)(ka + 1)) return false;

    int* ex = (int*)calloc((size_t)n, sizeof(int));
    for (int j = 0; j < ka; j++) ex[active[j]] = 1;
    const mpz_t* cf = mpoly_get_coef(eq, ex);
    if (!cf) { free(ex); return false; }
    mpz_set(c_full, *cf);

    bool ok = true, first = true;
    mpz_t a0; mpz_init(a0);
    for (int j = 0; j < ka && ok; j++) {
        ex[active[j]] = 0;
        const mpz_t* cc = mpoly_get_coef(eq, ex);
        ex[active[j]] = 1;
        if (!cc) { ok = false; break; }
        if (first) { mpz_neg(a0, *cc); first = false; }
        else {
            mpz_t ai; mpz_init(ai); mpz_neg(ai, *cc);
            if (mpz_cmp(ai, a0) != 0) ok = false;
            mpz_clear(ai);
        }
    }
    free(ex);
    if (!ok) { mpz_clear(a0); return false; }
    if (mpz_sgn(c_full) < 0) { mpz_neg(c_full, c_full); mpz_neg(a0, a0); }
    if (mpz_sgn(c_full) <= 0 || mpz_sgn(a0) <= 0) { mpz_clear(a0); return false; }
    mpz_set(a_out, a0); mpz_clear(a0);
    *ka_out = ka;
    return true;
}

/* Build a total order (smallest-first) over the active variables from the
 * ordering constraints; returns false if they do not form a chain. */
static bool si_build_total_order(const SICtx* c, const int* active, int ka, int* order) {
    bool used[SI_MAX_VARS]; for (int i = 0; i < SI_MAX_VARS; i++) used[i] = false;
    for (int step = 0; step < ka; step++) {
        int mn = -1, count = 0;
        for (int ai = 0; ai < ka; ai++) {
            int v = active[ai];
            if (used[v]) continue;
            bool has_below = false;
            for (int e = 0; e < c->n_ord; e++) {
                int a = c->ord_a[e], b = c->ord_b[e];   /* a <= b */
                if (b != v || used[a]) continue;
                for (int t = 0; t < ka; t++) if (active[t] == a) { has_below = true; break; }
                if (has_below) break;
            }
            if (!has_below) { mn = v; count++; }
        }
        if (count != 1) return false;
        order[step] = mn; used[mn] = true;
    }
    return true;
}

static bool si_solve_reciprocal(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int active[SI_MAX_VARS], ka = 0;
    mpz_t c_full, a; mpz_init(c_full); mpz_init(a);
    if (!si_reciprocal_detect(c->eq[0], c->n, active, &ka, c_full, a)) {
        mpz_clear(c_full); mpz_clear(a); return false;
    }
    int order[SI_MAX_VARS];
    if (!si_build_total_order(c, active, ka, order)) {
        mpz_clear(c_full); mpz_clear(a); return false;
    }
    int64_t full[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) full[i] = 0;
    /* p/q = R = c_full / a. */
    st->max_visits = SI_MAX_NODES;
    si_recip_rec(c, st, order, ka, 0, c_full, a, 1, full);
    mpz_clear(c_full); mpz_clear(a);
    return true;
}

/* --- Linear elimination + bilinear divisor solver (Pythagorean). --- */

/* Solve the reduced bilinear equation P (only vars u, w free) by factoring
 * M = b*c - a*d over its divisors.  Eliminated variables are reconstructed
 * numerically, in reverse elimination order, from their stored integer-
 * polynomial formulas.  Returns false if P is not a genuine hyperbola. */
static bool si_bilinear_divisor_solve(const MPoly* P, int u, int w, SICtx* c,
                                      SearchState* st, MPoly* const* formula_mp,
                                      const int* elim_order, int n_elim) {
    mpz_t a, b, cc, d, M, tmp; mpz_init(a); mpz_init(b); mpz_init(cc);
    mpz_init(d); mpz_init(M); mpz_init(tmp);
    int* ex = (int*)calloc((size_t)c->n, sizeof(int));
    ex[u] = 1; ex[w] = 1; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(a, *p); }
    ex[u] = 1; ex[w] = 0; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(b, *p); }
    ex[u] = 0; ex[w] = 1; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(cc, *p); }
    ex[u] = 0; ex[w] = 0; { const mpz_t* p = mpoly_get_coef(P, ex); if (p) mpz_set(d, *p); }
    free(ex);
    mpz_mul(M, b, cc); mpz_mul(tmp, a, d); mpz_sub(M, M, tmp);  /* M = bc - ad */
    bool handled = (mpz_sgn(a) != 0 && mpz_sgn(M) != 0);
    if (handled) {
        st->max_visits = SI_MAX_NODES;
        mpz_t absM; mpz_init(absM); mpz_abs(absM, M);
        Expr* dl = divisors_ordinary(absM);
        mpz_clear(absM);
        if (dl && dl->type == EXPR_FUNCTION) {
            mpz_t P1, Q1, uu, ww; mpz_init(P1); mpz_init(Q1); mpz_init(uu); mpz_init(ww);
            for (size_t i = 0; i < dl->data.function.arg_count && !st->overflow; i++) {
                if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                mpz_t dv; mpz_init(dv);
                expr_to_mpz(dl->data.function.args[i], dv);
                for (int sgn = 1; sgn >= -1; sgn -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(P1, dv); if (sgn < 0) mpz_neg(P1, P1);   /* P1 | M */
                    mpz_divexact(Q1, M, P1);                          /* Q1 = M/P1 */
                    mpz_sub(uu, P1, cc); mpz_sub(ww, Q1, b);          /* u=(P1-c)/a, w=(Q1-b)/a */
                    if (!mpz_divisible_p(uu, a) || !mpz_divisible_p(ww, a)) continue;
                    mpz_divexact(uu, uu, a); mpz_divexact(ww, ww, a);
                    if (!mpz_fits_slong_p(uu) || !mpz_fits_slong_p(ww)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int k = 0; k < c->n; k++) vals[k] = 0;
                    vals[u] = mpz_get_si(uu);
                    vals[w] = mpz_get_si(ww);
                    /* Reconstruct eliminated variables in reverse order: a
                     * later-eliminated variable's formula never references an
                     * earlier-eliminated one, so all references are resolved. */
                    bool okrec = true;
                    mpz_t rv; mpz_init(rv);
                    for (int e = n_elim - 1; e >= 0 && okrec; e--) {
                        int v = elim_order[e];
                        si_eval_mpoly(formula_mp[v], vals, rv);
                        if (!mpz_fits_slong_p(rv)) okrec = false;
                        else vals[v] = mpz_get_si(rv);
                    }
                    mpz_clear(rv);
                    if (okrec && si_verify(c, vals)) emit_full(st, vals);
                }
                mpz_clear(dv);
            }
            mpz_clear(P1); mpz_clear(Q1); mpz_clear(uu); mpz_clear(ww);
        }
        if (dl) expr_free(dl);
    }
    mpz_clear(a); mpz_clear(b); mpz_clear(cc); mpz_clear(d); mpz_clear(M); mpz_clear(tmp);
    return handled;
}

/* Try to reduce the equation system to a single bilinear equation in the
 * variables {keepU, keepW} by eliminating every other variable with a
 * unit-coefficient linear equation, then divisor-solve it. */
static bool si_attempt_pair(SICtx* c, SearchState* st, int keepU, int keepW) {
    Expr** conj; int ncj;
    flatten_conjuncts(c->original, &conj, &ncj);
    Expr* eqs[SI_MAX_VARS * 2]; int neq = 0;
    for (int i = 0; i < ncj && neq < (int)(sizeof(eqs)/sizeof(eqs[0])); i++)
        if (is_fun(conj[i], SYM_Equal, 2))
            eqs[neq++] = mk_fn2("Plus", expr_copy(conj[i]->data.function.args[0]),
                mk_fn2("Times", mk_int(-1), expr_copy(conj[i]->data.function.args[1])));

    bool elim_done[SI_MAX_VARS]; for (int i = 0; i < c->n; i++) elim_done[i] = false;
    MPoly* formula_mp[SI_MAX_VARS]; for (int i = 0; i < c->n; i++) formula_mp[i] = NULL;
    int elim_order[SI_MAX_VARS], n_elim = 0;
    bool alive[SI_MAX_VARS * 2]; for (int i = 0; i < neq; i++) alive[i] = true;

    /* Eliminate every non-kept variable that is unit-linear in some equation.
     * The reconstruction formula is kept as an MPoly (evaluated numerically per
     * candidate); only the equation reduction goes through the symbolic
     * substitution, which runs a handful of times, not per candidate. */
    for (bool progress = true; progress; ) {
        progress = false;
        for (int e = 0; e < neq && !progress; e++) {
            if (!alive[e]) continue;
            MPoly* P = si_resid_to_mpoly(eqs[e], c->var, c->n);
            if (!P) continue;
            for (int v = 0; v < c->n && !progress; v++) {
                if (v == keepU || v == keepW || elim_done[v] || mpoly_deg_var(P, v) != 1) continue;
                MPoly* lc = mpoly_coef_of_var(P, v, 1);
                bool unit = (mpoly_total_deg(lc) == 0 && lc->n_terms == 1
                             && (mpz_cmp_si(lc->coefs[0], 1) == 0
                                 || mpz_cmp_si(lc->coefs[0], -1) == 0));
                long coef = (lc->n_terms == 1) ? mpz_get_si(lc->coefs[0]) : 0;
                mpoly_free(lc);
                if (!unit) continue;
                MPoly* rest = mpoly_subst_var_int(P, v, 0);
                MPoly* fpoly = mpoly_scale_si(rest, -coef);   /* v = -(rest)/coef, coef = +/-1 */
                mpoly_free(rest);
                Expr* fexpr = mpoly_to_expr(fpoly, c->var);
                for (int e2 = 0; e2 < neq; e2++) {
                    if (!alive[e2] || e2 == e) continue;
                    Expr* rl = mk_list((Expr*[]){ mk_rule(expr_copy(c->var[v]), expr_copy(fexpr)) }, 1);
                    eqs[e2] = eval_and_free(internal_replace_all((Expr*[]){ eqs[e2], rl }, 2));
                }
                expr_free(fexpr);
                formula_mp[v] = fpoly; elim_order[n_elim++] = v;
                elim_done[v] = true; alive[e] = false; progress = true;
            }
            mpoly_free(P);
        }
    }

    /* Success needs every non-kept variable eliminated and exactly one live
     * equation, bilinear in {keepU, keepW}. */
    bool all_elim = true;
    for (int i = 0; i < c->n; i++)
        if (i != keepU && i != keepW && !elim_done[i]) { all_elim = false; break; }
    int live = -1, nlive = 0;
    for (int e = 0; e < neq; e++) if (alive[e]) { live = e; nlive++; }

    bool handled = false;
    if (all_elim && nlive == 1) {
        MPoly* P = si_resid_to_mpoly(eqs[live], c->var, c->n);
        if (P && mpoly_deg_var(P, keepU) <= 1 && mpoly_deg_var(P, keepW) <= 1
            && mpoly_total_deg(P) <= 2)
            handled = si_bilinear_divisor_solve(P, keepU, keepW, c, st,
                                                formula_mp, elim_order, n_elim);
        if (P) mpoly_free(P);
    }

    for (int i = 0; i < neq; i++) if (eqs[i]) expr_free(eqs[i]);
    for (int i = 0; i < c->n; i++) if (formula_mp[i]) mpoly_free(formula_mp[i]);
    return handled;
}

/* Try every pair of variables to keep; the first that reduces to a genuine
 * bilinear hyperbola wins. */
static bool si_solve_linelim_bilinear(SICtx* c, SearchState* st) {
    for (int u = 0; u < c->n; u++)
        for (int w = u + 1; w < c->n; w++)
            if (si_attempt_pair(c, st, u, w)) return true;
    return false;
}

/* --- Pell equation  x^2 - D y^2 == N  (N = +/-1) via continued fractions. --- */

/* Recognise a single equation whose only terms are xside^2, yside^2 and a
 * constant, with opposite-sign squares and the +square coefficient equal to 1.
 * Fills xside, yside, D (>0, non-square) and N (+/-1). */
static bool si_pell_detect(const MPoly* eq, int n, int* xside, int* yside,
                           mpz_t D, mpz_t N) {
    int active[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) active[ka++] = v;
    if (ka != 2) return false;
    int A = active[0], B = active[1];
    mpz_t cA, cB, e; mpz_init_set_ui(cA, 0); mpz_init_set_ui(cB, 0); mpz_init_set_ui(e, 0);
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++) if (v != A && v != B && ex[v] != 0) ok = false;
        if (!ok) break;
        int dA = ex[A], dB = ex[B];
        if (dA == 2 && dB == 0) mpz_set(cA, eq->coefs[t]);
        else if (dA == 0 && dB == 2) mpz_set(cB, eq->coefs[t]);
        else if (dA == 0 && dB == 0) mpz_set(e, eq->coefs[t]);
        else ok = false;                                     /* linear or cross term */
    }
    if (ok && (mpz_sgn(cA) == 0 || mpz_sgn(cB) == 0 || mpz_sgn(cA) == mpz_sgn(cB)))
        ok = false;
    if (ok) {
        if (mpz_sgn(cA) > 0) {
            if (mpz_cmp_si(cA, 1) != 0) ok = false;
            else { *xside = A; *yside = B; mpz_neg(D, cB); mpz_neg(N, e); }
        } else {
            if (mpz_cmp_si(cB, 1) != 0) ok = false;
            else { *xside = B; *yside = A; mpz_neg(D, cA); mpz_neg(N, e); }
        }
    }
    mpz_clear(cA); mpz_clear(cB); mpz_clear(e);
    if (!ok) return false;
    if (mpz_sgn(D) <= 0 || mpz_perfect_square_p(D)) return false;
    return (mpz_cmp_si(N, 1) == 0 || mpz_cmp_si(N, -1) == 0);
}

/* Continued fraction of sqrt(D): sets the fundamental unit (u,v) of
 * x^2 - D y^2 = 1, and the first convergent (bx,by) reaching value N.
 * Returns true if the N-base was found. */
static bool si_pell_cf(const mpz_t D, const mpz_t N, mpz_t u, mpz_t v,
                       mpz_t bx, mpz_t by) {
    mpz_t a0, m, d, a, hprev, h, kprev, k, val, t2;
    mpz_init(a0); mpz_sqrt(a0, D);
    mpz_init_set_ui(m, 0); mpz_init_set_ui(d, 1); mpz_init_set(a, a0);
    mpz_init_set_ui(hprev, 1); mpz_init_set(h, a0);
    mpz_init_set_ui(kprev, 0); mpz_init_set_ui(k, 1);
    mpz_init(val); mpz_init(t2);
    bool haveU = false, haveB = false;
    for (int iter = 0; iter < 200000 && !haveU; iter++) {
        mpz_mul(val, h, h); mpz_mul(t2, k, k); mpz_mul(t2, t2, D);
        mpz_sub(val, val, t2);                               /* h^2 - D k^2 */
        if (!haveB && mpz_cmp(val, N) == 0) { mpz_set(bx, h); mpz_set(by, k); haveB = true; }
        if (mpz_cmp_si(val, 1) == 0) { mpz_set(u, h); mpz_set(v, k); haveU = true; break; }
        /* next term: m=d*a-m; d=(D-m^2)/d; a=floor((a0+m)/d) */
        mpz_mul(t2, d, a); mpz_sub(m, t2, m);
        mpz_mul(t2, m, m); mpz_sub(t2, D, t2); mpz_divexact(d, t2, d);
        mpz_add(t2, a0, m); mpz_fdiv_q(a, t2, d);
        mpz_mul(t2, a, h); mpz_add(t2, t2, hprev); mpz_set(hprev, h); mpz_set(h, t2);
        mpz_mul(t2, a, k); mpz_add(t2, t2, kprev); mpz_set(kprev, k); mpz_set(k, t2);
    }
    mpz_clear(a0); mpz_clear(m); mpz_clear(d); mpz_clear(a);
    mpz_clear(hprev); mpz_clear(h); mpz_clear(kprev); mpz_clear(k);
    mpz_clear(val); mpz_clear(t2);
    return haveB;
}

/* Emit (+/-X, +/-Y) at (xside, yside), each sign combo verified. */
static void si_pell_emit(SICtx* c, SearchState* st, int xs, int ys,
                         const mpz_t X, const mpz_t Y) {
    if (!mpz_fits_slong_p(X) || !mpz_fits_slong_p(Y)) return;
    int64_t xa = mpz_get_si(X), ya = mpz_get_si(Y);
    for (int sx = 1; sx >= -1; sx -= 2)
        for (int sy = 1; sy >= -1; sy -= 2) {
            int64_t vals[SI_MAX_VARS];
            for (int i = 0; i < c->n; i++) vals[i] = 0;
            vals[xs] = sx * xa; vals[ys] = sy * ya;
            if (si_verify(c, vals)) emit_full(st, vals);
        }
}

static bool si_solve_pell(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_pell_detect(c->eq[0], c->n, &xs, &ys, D, N)) { mpz_clear(D); mpz_clear(N); return false; }
    /* A finite bound is required to terminate the (otherwise infinite) orbit. */
    if (!c->has_hi[xs] && !c->has_hi[ys]) { mpz_clear(D); mpz_clear(N); return false; }

    mpz_t u, v, bx, by; mpz_init(u); mpz_init(v); mpz_init(bx); mpz_init(by);
    bool haveB = si_pell_cf(D, N, u, v, bx, by);
    st->max_visits = SI_MAX_NODES;

    /* Base class of solutions to x^2 - D y^2 = N. */
    mpz_t cx, cy, nx, ny, t1, t2; mpz_init(cx); mpz_init(cy);
    mpz_init(nx); mpz_init(ny); mpz_init(t1); mpz_init(t2);
    bool go = true;
    if (mpz_cmp_si(N, 1) == 0) {
        mpz_set_ui(cx, 1); mpz_set_ui(cy, 0);                /* trivial (1,0) */
        si_pell_emit(c, st, xs, ys, cx, cy);
        mpz_set(cx, u); mpz_set(cy, v);                      /* first nontrivial = U */
    } else if (haveB) {
        mpz_set(cx, bx); mpz_set(cy, by);                    /* N = -1 base */
    } else {
        go = false;                                          /* x^2 - D y^2 = -1 unsolvable */
    }

    for (int guard = 0; go && guard < 100000 && !st->overflow; guard++) {
        bool past = (c->has_hi[xs] && mpz_cmp_si(cx, c->hi[xs]) > 0)
                 || (c->has_hi[ys] && mpz_cmp_si(cy, c->hi[ys]) > 0);
        if (past) break;
        si_pell_emit(c, st, xs, ys, cx, cy);
        /* compose with U: (x,y) -> (x u + D y v, x v + y u) */
        mpz_mul(nx, cx, u); mpz_mul(t1, cy, v); mpz_mul(t1, t1, D); mpz_add(nx, nx, t1);
        mpz_mul(ny, cx, v); mpz_mul(t1, cy, u); mpz_add(ny, ny, t1);
        mpz_set(cx, nx); mpz_set(cy, ny);
    }

    mpz_clear(cx); mpz_clear(cy); mpz_clear(nx); mpz_clear(ny); mpz_clear(t1); mpz_clear(t2);
    mpz_clear(u); mpz_clear(v); mpz_clear(bx); mpz_clear(by);
    mpz_clear(D); mpz_clear(N);
    return true;
}

/* --- Unbounded Pell  x^2 - D y^2 == 1,  x > 0 && y > 0  -> parametric family. ---
 *
 * With no finite bound the positive-orthant solution set is infinite: the k-th
 * solution (k >= 1) is  x_k + y_k sqrt(D) = (x1 + y1 sqrt(D))^k  for the
 * fundamental unit (x1, y1).  Emit the closed form
 *   x -> ((x1+y1 Sqrt[D])^C[1] + (x1-y1 Sqrt[D])^C[1]) / 2,
 *   y -> ((x1+y1 Sqrt[D])^C[1] - (x1-y1 Sqrt[D])^C[1]) / (2 Sqrt[D]),
 * each a ConditionalExpression valid for the integer parameter C[1] >= 1 --
 * mirroring Mathematica's representation.  Returns the owned family, or NULL to
 * decline (bounded, wrong sign pattern, or N = -1, all handled elsewhere). */
static Expr* si_solve_pell_parametric(SICtx* c) {
    if (c->neq != 1) return NULL;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_pell_detect(c->eq[0], c->n, &xs, &ys, D, N)) { mpz_clear(D); mpz_clear(N); return NULL; }
    bool unbounded = !c->has_hi[xs] && !c->has_hi[ys];
    bool posx = c->has_lo[xs] && c->lo[xs] >= 1;
    bool posy = c->has_lo[ys] && c->lo[ys] >= 1;
    if (!(unbounded && posx && posy && mpz_cmp_si(N, 1) == 0)
        || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) {
        mpz_clear(D); mpz_clear(N); return NULL;
    }
    mpz_t u, v, bx, by; mpz_init(u); mpz_init(v); mpz_init(bx); mpz_init(by);
    si_pell_cf(D, N, u, v, bx, by);                  /* fundamental unit (x1,y1) = (u,v) */

    Expr* sqrtD = mk_fn1("Sqrt", mk_mpz(D));
    mpz_t ny; mpz_init(ny); mpz_neg(ny, v);
    /* uu = x1 + y1 Sqrt[D],  vv = x1 - y1 Sqrt[D]. */
    Expr* uu = mk_fn2("Plus", mk_mpz(u), mk_fn2("Times", mk_mpz(v), expr_copy(sqrtD)));
    Expr* vv = mk_fn2("Plus", mk_mpz(u), mk_fn2("Times", mk_mpz(ny), expr_copy(sqrtD)));
    mpz_clear(ny);
    Expr* uk = mk_fn2("Power", uu, mk_fn1("C", mk_int(1)));
    Expr* vk = mk_fn2("Power", vv, mk_fn1("C", mk_int(1)));
    /* xval = (uk + vk)/2 */
    Expr* xval = mk_fn2("Times",
        mk_fn2("Plus", expr_copy(uk), expr_copy(vk)),
        mk_fn2("Power", mk_int(2), mk_int(-1)));
    /* yval = (uk - vk)/(2 Sqrt[D]) */
    Expr* yval = mk_fn2("Times",
        mk_fn2("Plus", uk, mk_fn2("Times", mk_int(-1), vk)),
        mk_fn2("Power", mk_fn2("Times", mk_int(2), expr_copy(sqrtD)), mk_int(-1)));
    expr_free(sqrtD);

    Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(1));
    Expr* cex = mk_fn2("ConditionalExpression", xval, expr_copy(cond));
    Expr* cey = mk_fn2("ConditionalExpression", yval, cond);

    Expr* rules[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) rules[i] = NULL;
    rules[xs] = mk_rule(expr_copy(c->var[xs]), cex);
    rules[ys] = mk_rule(expr_copy(c->var[ys]), cey);
    /* Emit in variable order (xs, ys are the only two active). */
    Expr* rlist[SI_MAX_VARS]; int nr = 0;
    for (int i = 0; i < c->n; i++) if (rules[i]) rlist[nr++] = rules[i];
    Expr* tuple = mk_list(rlist, (size_t)nr);
    Expr* result = eval_and_free(mk_list((Expr*[]){ tuple }, 1));

    mpz_clear(u); mpz_clear(v); mpz_clear(bx); mpz_clear(by);
    mpz_clear(D); mpz_clear(N);
    return result;
}

/* --- Generalised Pell  x^2 - D y^2 == N  (any nonzero N), unbounded positive
 *     orthant -> a parametric family PER solution class. ---
 *
 * Detect x^2 - D y^2 == N with a +1 square coefficient, D > 0 non-square, and
 * ANY nonzero N (si_pell_detect handles only N = +/-1). */
static bool si_genpell_detect(const MPoly* eq, int n, int* xside, int* yside,
                              mpz_t D, mpz_t N) {
    int active[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) active[ka++] = v;
    if (ka != 2) return false;
    int A = active[0], B = active[1];
    mpz_t cA, cB, e; mpz_init_set_ui(cA, 0); mpz_init_set_ui(cB, 0); mpz_init_set_ui(e, 0);
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++) if (v != A && v != B && ex[v] != 0) ok = false;
        if (!ok) break;
        int dA = ex[A], dB = ex[B];
        if (dA == 2 && dB == 0) mpz_set(cA, eq->coefs[t]);
        else if (dA == 0 && dB == 2) mpz_set(cB, eq->coefs[t]);
        else if (dA == 0 && dB == 0) mpz_set(e, eq->coefs[t]);
        else ok = false;                                     /* linear or cross term */
    }
    if (ok && (mpz_sgn(cA) == 0 || mpz_sgn(cB) == 0 || mpz_sgn(cA) == mpz_sgn(cB)))
        ok = false;
    if (ok) {
        if (mpz_sgn(cA) > 0) {
            if (mpz_cmp_si(cA, 1) != 0) ok = false;
            else { *xside = A; *yside = B; mpz_neg(D, cB); mpz_neg(N, e); }
        } else {
            if (mpz_cmp_si(cB, 1) != 0) ok = false;
            else { *xside = B; *yside = A; mpz_neg(D, cA); mpz_neg(N, e); }
        }
    }
    mpz_clear(cA); mpz_clear(cB); mpz_clear(e);
    if (!ok) return false;
    if (mpz_sgn(D) <= 0 || mpz_perfect_square_p(D)) return false;
    return (mpz_sgn(N) != 0);
}

/* Largest y a Nagell fundamental solution can have.  Ceiling over-estimate of
 * u*sqrt(|N|/(2(t +/- 1))) (+ for N>0, - for N<0), a safe bound validated
 * against brute force.  Sets *ok false if it would be impractically large. */
static void si_genpell_ybound(const mpz_t N, const mpz_t t, const mpz_t u,
                              mpz_t out, bool* ok) {
    mpz_t num, den, absN; mpz_init(num); mpz_init(den); mpz_init(absN);
    mpz_abs(absN, N);
    mpz_mul(num, u, u); mpz_mul(num, num, absN);             /* u^2 |N| */
    if (mpz_sgn(N) > 0) { mpz_add_ui(den, t, 1); }           /* 2(t+1) */
    else { mpz_sub_ui(den, t, 1); }                          /* 2(t-1) */
    mpz_mul_ui(den, den, 2);
    if (mpz_sgn(den) <= 0) { mpz_set_ui(out, 0); *ok = true; }
    else {
        mpz_fdiv_q(num, num, den);
        mpz_sqrt(out, num);
        mpz_add_ui(out, out, 2);                             /* safety margin */
    }
    /* Guard against a pathological search (huge fundamental unit and |N|). */
    *ok = (mpz_cmp_ui(out, 50000000UL) <= 0);
    mpz_clear(num); mpz_clear(den); mpz_clear(absN);
}

/* Compose (x, y) with the fundamental unit (t, u): (x,y) <- (x t + D y u,
 * x u + y t). */
static void si_pell_step(mpz_t x, mpz_t y, const mpz_t t, const mpz_t u,
                         const mpz_t D, mpz_t s1, mpz_t s2) {
    mpz_mul(s1, x, t); mpz_mul(s2, D, y); mpz_mul(s2, s2, u); mpz_add(s1, s1, s2); /* x t + D y u */
    mpz_mul(s2, x, u); mpz_mul(x, y, t); mpz_add(s2, s2, x);                       /* x u + y t */
    mpz_set(x, s1); mpz_set(y, s2);
}

/* Fundamental unit (t, u) of x^2 - D y^2 = 1, and the minimal positive-orthant
 * representative (x > 0, y > 0) of every solution class of x^2 - D y^2 == N.
 * The bases are written (deduplicated) into bx[]/by[]; *nb is their count.
 * Returns false to decline (no unit found, or the Nagell box too large). */
static bool si_genpell_bases(const mpz_t D, const mpz_t N, mpz_t t, mpz_t u,
                             mpz_t* bx, mpz_t* by, int* nb, int cap) {
    *nb = 0;
    mpz_t one, dummy1, dummy2;
    mpz_init_set_ui(one, 1); mpz_init(dummy1); mpz_init(dummy2);
    if (!si_pell_cf(D, one, t, u, dummy1, dummy2)) {  /* (t,u): x^2 - D y^2 = 1 */
        mpz_clear(one); mpz_clear(dummy1); mpz_clear(dummy2); return false;
    }
    mpz_clear(one); mpz_clear(dummy1); mpz_clear(dummy2);

    mpz_t Ymax; mpz_init(Ymax); bool okb = true;
    si_genpell_ybound(N, t, u, Ymax, &okb);
    if (!okb || !mpz_fits_slong_p(Ymax)) { mpz_clear(Ymax); return false; }
    int64_t ymax = mpz_get_si(Ymax);
    mpz_clear(Ymax);

    mpz_t val, root, x, y, s1, s2, px, py;
    mpz_inits(val, root, x, y, s1, s2, px, py, NULL);
    for (int64_t yy = 0; yy <= ymax && *nb < cap; yy++) {
        mpz_set_si(y, yy); mpz_mul(val, y, y); mpz_mul(val, val, D);
        mpz_add(val, val, N);                                /* N + D y^2 */
        if (mpz_sgn(val) <= 0) continue;
        if (!mpz_perfect_square_p(val)) continue;
        mpz_sqrt(root, val);                                 /* x0 > 0 */
        /* All four sign combos (+/-x0, +/-yy) are candidate fundamentals of
         * (possibly distinct) classes; advance each into the positive orthant. */
        for (int sx = 1; sx >= -1; sx -= 2)
            for (int sy = 1; sy >= -1; sy -= 2) {
                mpz_set(x, root); if (sx < 0) mpz_neg(x, x);
                mpz_set_si(s1, sy * yy); mpz_set(y, s1);
                int guard = 0;
                while (!(mpz_sgn(x) > 0 && mpz_sgn(y) > 0) && guard < 1000) {
                    si_pell_step(x, y, t, u, D, s1, s2); guard++;
                }
                if (!(mpz_sgn(x) > 0 && mpz_sgn(y) > 0)) continue;
                /* Reduce to the MINIMAL positive-orthant representative: apply
                 * eps^{-1} = (t, -u) while the result stays in the orthant, so
                 * two solutions of the same class (e.g. (1,1) and (5,3)=(1,1)eps
                 * for D=3, N=-2) collapse to one base rather than two overlapping
                 * families.  eps^{-1}: x' = x t - D y u,  y' = -x u + y t. */
                for (int g2 = 0; g2 < 1000; g2++) {
                    mpz_mul(px, x, t); mpz_mul(s2, D, y); mpz_mul(s2, s2, u);
                    mpz_sub(px, px, s2);                       /* x' = x t - D y u */
                    mpz_mul(py, y, t); mpz_mul(s2, x, u);
                    mpz_sub(py, py, s2);                       /* y' = y t - x u */
                    if (mpz_sgn(px) > 0 && mpz_sgn(py) > 0) { mpz_set(x, px); mpz_set(y, py); }
                    else break;
                }
                bool dup = false;
                for (int i = 0; i < *nb; i++)
                    if (mpz_cmp(bx[i], x) == 0 && mpz_cmp(by[i], y) == 0) { dup = true; break; }
                if (dup || *nb >= cap) continue;
                mpz_set(bx[*nb], x); mpz_set(by[*nb], y); (*nb)++;
            }
    }
    mpz_clears(val, root, x, y, s1, s2, px, py, NULL);
    return true;
}

/* Order bases ascending by (x, then y) for a deterministic emission order. */
static void si_genpell_sort_bases(mpz_t* bx, mpz_t* by, int nb) {
    for (int i = 0; i < nb; i++)
        for (int j = i + 1; j < nb; j++) {
            int c = mpz_cmp(bx[i], bx[j]);
            if (c > 0 || (c == 0 && mpz_cmp(by[i], by[j]) > 0)) {
                mpz_swap(bx[i], bx[j]); mpz_swap(by[i], by[j]);
            }
        }
}

/* Build the class family tuple for base (a, b), fundamental unit (t, u):
 *   x -> (P eps^C + Pbar epsbar^C) / 2,
 *   y -> (P eps^C - Pbar epsbar^C) / (2 Sqrt[D]),
 * where P = a + b Sqrt[D], eps = t + u Sqrt[D] (bars negate the Sqrt[D] part),
 * each a ConditionalExpression on C[1] >= 0.  Rules placed at xs/ys. */
static Expr* si_genpell_family(SICtx* c, int xs, int ys, const mpz_t a,
                               const mpz_t b, const mpz_t t, const mpz_t u,
                               const mpz_t D) {
    Expr* sqrtD = mk_fn1("Sqrt", mk_mpz(D));
    mpz_t nb, nu; mpz_init(nb); mpz_init(nu); mpz_neg(nb, b); mpz_neg(nu, u);
    Expr* P    = mk_fn2("Plus", mk_mpz(a), mk_fn2("Times", mk_mpz(b),  expr_copy(sqrtD)));
    Expr* Pbar = mk_fn2("Plus", mk_mpz(a), mk_fn2("Times", mk_mpz(nb), expr_copy(sqrtD)));
    Expr* eps    = mk_fn2("Plus", mk_mpz(t), mk_fn2("Times", mk_mpz(u),  expr_copy(sqrtD)));
    Expr* epsbar = mk_fn2("Plus", mk_mpz(t), mk_fn2("Times", mk_mpz(nu), expr_copy(sqrtD)));
    mpz_clear(nb); mpz_clear(nu);
    Expr* epsC    = mk_fn2("Power", eps,    mk_fn1("C", mk_int(1)));
    Expr* epsbarC = mk_fn2("Power", epsbar, mk_fn1("C", mk_int(1)));
    Expr* Pe    = mk_fn2("Times", P, epsC);
    Expr* Pbare = mk_fn2("Times", Pbar, epsbarC);
    /* xval = (Pe + Pbare)/2 */
    Expr* xval = mk_fn2("Times",
        mk_fn2("Plus", expr_copy(Pe), expr_copy(Pbare)),
        mk_fn2("Power", mk_int(2), mk_int(-1)));
    /* yval = (Pe - Pbare)/(2 Sqrt[D]) */
    Expr* yval = mk_fn2("Times",
        mk_fn2("Plus", Pe, mk_fn2("Times", mk_int(-1), Pbare)),
        mk_fn2("Power", mk_fn2("Times", mk_int(2), expr_copy(sqrtD)), mk_int(-1)));
    expr_free(sqrtD);

    Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(0));
    Expr* cex = mk_fn2("ConditionalExpression", xval, expr_copy(cond));
    Expr* cey = mk_fn2("ConditionalExpression", yval, cond);

    Expr* rules[SI_MAX_VARS];
    for (int i = 0; i < c->n; i++) rules[i] = NULL;
    rules[xs] = mk_rule(expr_copy(c->var[xs]), cex);
    rules[ys] = mk_rule(expr_copy(c->var[ys]), cey);
    Expr* rlist[SI_MAX_VARS]; int nr = 0;
    for (int i = 0; i < c->n; i++) if (rules[i]) rlist[nr++] = rules[i];
    return mk_list(rlist, (size_t)nr);
}

/* Unbounded x^2 - D y^2 == N (N != +1: any N < 0 or N >= 2) with x > 0 && y > 0
 * -> one parametric family per solution class.  Returns the owned family List
 * (empty {} when the equation is provably unsolvable, including negative Pell
 * over a D whose sqrt has even CF period), or NULL to decline (not this shape,
 * bounded, N = +1, or the Nagell search too large -- handled elsewhere). */
static Expr* si_solve_genpell_parametric(SICtx* c) {
    if (c->neq != 1) return NULL;
    int xs, ys; mpz_t D, N; mpz_init(D); mpz_init(N);
    if (!si_genpell_detect(c->eq[0], c->n, &xs, &ys, D, N)
        || mpz_sgn(N) == 0 || mpz_cmp_si(N, 1) == 0) {  /* N = +1 -> pell_parametric */
        mpz_clear(D); mpz_clear(N); return NULL;
    }
    bool unbounded = !c->has_hi[xs] && !c->has_hi[ys];
    bool posx = c->has_lo[xs] && c->lo[xs] >= 1;
    bool posy = c->has_lo[ys] && c->lo[ys] >= 1;
    if (!(unbounded && posx && posy) || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) {
        mpz_clear(D); mpz_clear(N); return NULL;
    }

    enum { GP_CAP = 4096 };
    mpz_t t, u; mpz_init(t); mpz_init(u);
    mpz_t *bx = (mpz_t*)malloc(sizeof(mpz_t) * GP_CAP);
    mpz_t *by = (mpz_t*)malloc(sizeof(mpz_t) * GP_CAP);
    for (int i = 0; i < GP_CAP; i++) { mpz_init(bx[i]); mpz_init(by[i]); }
    int nb = 0;
    bool solved = si_genpell_bases(D, N, t, u, bx, by, &nb, GP_CAP);

    Expr* result = NULL;
    if (solved) {
        si_genpell_sort_bases(bx, by, nb);
        if (nb == 0) {
            result = mk_list(NULL, 0);            /* provably unsolvable -> {} */
        } else {
            Expr** tuples = (Expr**)malloc(sizeof(Expr*) * (size_t)nb);
            for (int i = 0; i < nb; i++)
                tuples[i] = si_genpell_family(c, xs, ys, bx[i], by[i], t, u, D);
            result = eval_and_free(mk_list(tuples, (size_t)nb));
            free(tuples);
        }
    }

    for (int i = 0; i < GP_CAP; i++) { mpz_clear(bx[i]); mpz_clear(by[i]); }
    free(bx); free(by);
    mpz_clear(t); mpz_clear(u); mpz_clear(D); mpz_clear(N);
    return result;                               /* NULL -> decline (unevaluated) */
}

static bool si_linear_detect(const MPoly* eq, int n, mpz_t* a, mpz_t b);

/* Fraction-free (Bareiss) determinant of an m x m integer matrix.  M is
 * destroyed; `out` (pre-init'd) receives the exact determinant. */
static void si_det_bareiss(mpz_t** M, int m, mpz_t out) {
    if (m == 0) { mpz_set_ui(out, 1); return; }
    mpz_t prev, t1, t2; mpz_init_set_ui(prev, 1); mpz_init(t1); mpz_init(t2);
    int sign = 1;
    bool zero = false;
    for (int k = 0; k < m; k++) {
        if (mpz_sgn(M[k][k]) == 0) {                 /* pivot: swap in a nonzero row */
            int r = -1;
            for (int i = k + 1; i < m; i++) if (mpz_sgn(M[i][k]) != 0) { r = i; break; }
            if (r < 0) { zero = true; break; }
            for (int j = 0; j < m; j++) mpz_swap(M[k][j], M[r][j]);
            sign = -sign;
        }
        for (int i = k + 1; i < m; i++)
            for (int j = k + 1; j < m; j++) {
                mpz_mul(t1, M[i][j], M[k][k]);
                mpz_mul(t2, M[i][k], M[k][j]);
                mpz_sub(t1, t1, t2);
                mpz_divexact(M[i][j], t1, prev);     /* Bareiss: exact division */
            }
        mpz_set(prev, M[k][k]);
    }
    if (zero) mpz_set_ui(out, 0);
    else { mpz_set(out, M[m - 1][m - 1]); if (sign < 0) mpz_neg(out, out); }
    mpz_clear(prev); mpz_clear(t1); mpz_clear(t2);
}

/* --- Homogeneous linear SYSTEM with positivity -> parametric ray. ---
 *
 * A system of  m = n-1  homogeneous linear equations in n positive unknowns has
 * (generically) a 1-dimensional integer kernel spanned by a primitive vector v
 * (the generalised cross product: v_j = (-1)^j * det(A with column j removed)).
 * If v (or -v) is entirely positive the positive solutions are exactly
 * { k v : k = C[1] >= 1 }; if v has mixed signs the positive orthant meets the
 * kernel only at 0, so there are no positive solutions.  This settles rational
 * linear systems like  w == 5/6 x + y && x == 9/20 y + z && y == 13/42 z + w
 * (denominators cleared by the MPoly conversion). */
static Expr* si_solve_linear_system_ray(SICtx* c) {
    int n = c->n;
    if (c->neq != n - 1 || n < 2 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < n; i++)                       /* every variable strictly positive */
        if (!(c->has_lo[i] && c->lo[i] >= 1)) return NULL;

    int m = c->neq;
    /* Build the m x n integer coefficient matrix; require every equation to be
     * linear and HOMOGENEOUS (constant term 0). */
    mpz_t* a = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t bb; for (int i = 0; i < n; i++) mpz_init(a[i]);
    mpz_init(bb);
    mpz_t** A = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
    for (int q = 0; q < m; q++) { A[q] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n); for (int i = 0; i < n; i++) mpz_init(A[q][i]); }
    bool ok = true;
    for (int q = 0; q < m && ok; q++) {
        if (!si_linear_detect(c->eq[q], n, a, bb) || mpz_sgn(bb) != 0) ok = false;
        else for (int i = 0; i < n; i++) mpz_set(A[q][i], a[i]);
    }

    Expr* result = NULL;
    if (ok) {
        /* v_j = (-1)^j det(A without column j). */
        mpz_t* v = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init(v[i]);
        mpz_t** Msub = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
        for (int q = 0; q < m; q++) { Msub[q] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m); for (int i = 0; i < m; i++) mpz_init(Msub[q][i]); }
        mpz_t det; mpz_init(det);
        for (int j = 0; j < n; j++) {
            int cc = 0;
            for (int col = 0; col < n; col++) { if (col == j) continue;
                for (int q = 0; q < m; q++) mpz_set(Msub[q][cc], A[q][col]);
                cc++;
            }
            si_det_bareiss(Msub, m, det);
            mpz_set(v[j], det);
            if (j & 1) mpz_neg(v[j], v[j]);
        }
        /* Primitive + sign-normalise (make the first nonzero positive). */
        mpz_t g; mpz_init_set_ui(g, 0);
        for (int i = 0; i < n; i++) mpz_gcd(g, g, v[i]);
        bool nonzero = mpz_sgn(g) != 0;
        if (nonzero) {
            for (int i = 0; i < n; i++) mpz_divexact(v[i], v[i], g);
            int lead = 0; while (lead < n && mpz_sgn(v[lead]) == 0) lead++;
            if (lead < n && mpz_sgn(v[lead]) < 0) for (int i = 0; i < n; i++) mpz_neg(v[i], v[i]);
        }
        bool all_pos = nonzero;
        for (int i = 0; i < n && all_pos; i++) if (mpz_sgn(v[i]) <= 0) all_pos = false;

        if (!nonzero) {
            result = NULL;                            /* degenerate kernel: decline */
        } else if (!all_pos) {
            result = mk_list(NULL, 0);                /* no positive solution: {} */
        } else {
            Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) {
                Expr* term = mpz_cmp_ui(v[i], 1) == 0
                    ? mk_fn1("C", mk_int(1))
                    : mk_fn2("Times", mk_mpz(v[i]), mk_fn1("C", mk_int(1)));
                Expr* cond = mk_fn2("GreaterEqual", mk_fn1("C", mk_int(1)), mk_int(1));
                rules[i] = mk_rule(expr_copy(c->var[i]),
                                   mk_fn2("ConditionalExpression", term, cond));
            }
            Expr* tuple = mk_list(rules, (size_t)n); free(rules);
            result = eval_and_free(mk_list((Expr*[]){ tuple }, 1));
        }
        mpz_clear(det); mpz_clear(g);
        for (int q = 0; q < m; q++) { for (int i = 0; i < m; i++) mpz_clear(Msub[q][i]); free(Msub[q]); }
        free(Msub);
        for (int i = 0; i < n; i++) mpz_clear(v[i]);
        free(v);
    }

    for (int q = 0; q < m; q++) { for (int i = 0; i < n; i++) mpz_clear(A[q][i]); free(A[q]); }
    free(A);
    for (int i = 0; i < n; i++) mpz_clear(a[i]);
    free(a); mpz_clear(bb);
    return result;
}

/* --- Prouhet-Tarry-Escott  ->  {}  via Newton's identities. ---
 *
 * If two disjoint k-element variable groups G1, G2 satisfy the equal power sums
 *   Sum_{G1} v^j == Sum_{G2} v^j   for every j = 1..k,
 * then (Newton's identities) their elementary symmetric functions agree, so the
 * two groups are the SAME multiset.  With a strict ordering inside each group
 * this forces  G1[i] == G2[i]  elementwise, so a disequation like  a != d
 * between the corresponding smallest elements is contradictory: the system has
 * NO solution.  This settles the size-3, degree-1..3 case
 *   a+b+c==d+e+f && a^2+..==d^2+.. && a^3+..==d^3+.. && 0<a<b<c && 0<d<e<f && a!=d
 * as {} despite every variable being unbounded.  Returns the empty List when
 * proved empty, or NULL to decline (no forcing disequation -> parametric). */
static bool si_pset_eq(const bool* x, const bool* y, int n) {
    for (int i = 0; i < n; i++) if (x[i] != y[i]) return false;
    return true;
}
static Expr* si_solve_power_sum_equal(SICtx* c) {
    int k = c->neq, n = c->n;
    if (k < 2 || n != 2 * k || k > SI_MAX_VARS) return NULL;

    bool Pset[SI_MAX_VARS], Mset[SI_MAX_VARS];
    bool have_sets = false;
    bool jseen[SI_MAX_VARS + 1];
    for (int j = 0; j <= k; j++) jseen[j] = false;

    for (int q = 0; q < k; q++) {
        const MPoly* eq = c->eq[q];
        bool p[SI_MAX_VARS], m[SI_MAX_VARS];
        for (int v = 0; v < n; v++) { p[v] = false; m[v] = false; }
        int np = 0, nm = 0, thisj = -1; bool ok = true;
        for (size_t t = 0; t < eq->n_terms && ok; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            int nv = -1, e = 0;
            for (int v = 0; v < n; v++) if (ex[v] > 0) { if (nv >= 0) { ok = false; break; } nv = v; e = ex[v]; }
            if (!ok) break;
            if (nv < 0) { ok = false; break; }             /* constant term: not PTE */
            if (thisj < 0) thisj = e; else if (e != thisj) { ok = false; break; }
            long cf = mpz_cmp_si(eq->coefs[t], 1) == 0 ? 1
                    : mpz_cmp_si(eq->coefs[t], -1) == 0 ? -1 : 0;
            if (cf == 1) { p[nv] = true; np++; }
            else if (cf == -1) { m[nv] = true; nm++; }
            else { ok = false; break; }
        }
        if (!ok || np != k || nm != k || thisj < 1 || thisj > k || jseen[thisj]) return NULL;
        jseen[thisj] = true;
        if (!have_sets) { for (int v = 0; v < n; v++) { Pset[v] = p[v]; Mset[v] = m[v]; } have_sets = true; }
        else if (!((si_pset_eq(p, Pset, n) && si_pset_eq(m, Mset, n))     /* same convention */
                 || (si_pset_eq(p, Mset, n) && si_pset_eq(m, Pset, n))))  /* or the sign flip */
            return NULL;
    }
    for (int j = 1; j <= k; j++) if (!jseen[j]) return NULL;              /* need degrees 1..k */

    int G1[SI_MAX_VARS], G2[SI_MAX_VARS], n1 = 0, n2 = 0;
    for (int v = 0; v < n; v++) { if (Pset[v]) G1[n1++] = v; if (Mset[v]) G2[n2++] = v; }
    if (n1 != k || n2 != k) return NULL;
    int ord1[SI_MAX_VARS], ord2[SI_MAX_VARS];
    if (!si_build_total_order(c, G1, k, ord1) || !si_build_total_order(c, G2, k, ord2)) return NULL;

    /* Newton forces ord1[i] == ord2[i]; a disequation between them is a proof of
     * emptiness. */
    for (int d = 0; d < c->n_neq; d++) {
        int a = c->neq_a[d], b = c->neq_b[d];
        for (int i = 0; i < k; i++)
            if ((a == ord1[i] && b == ord2[i]) || (a == ord2[i] && b == ord1[i]))
                return mk_list(NULL, 0);                                  /* {} */
    }
    return NULL;                                     /* no forcing disequation: decline */
}

static Expr* si_exp_one_tuple(Expr** var, const int64_t* vals, int n);
static bool si_verify_symbolic(Expr* expr, Expr** var, const int64_t* vals, int n);

static long si_lgcd(long a, long b) { a = labs(a); b = labs(b); while (b) { long t = a % b; a = b; b = t; } return a; }

/* Is m (|m|) squarefree? */
static bool si_is_squarefree_l(long m) {
    m = labs(m);
    if (m == 0) return false;
    for (long p = 2; p * p <= m; p++)
        if (m % p == 0) { m /= p; if (m % p == 0) return false; }
    return true;
}

/* Class number of the imaginary quadratic field of fundamental discriminant
 * D = -4c (c >= 1 squarefree, c = 1,2 mod 4), by counting primitive reduced
 * positive-definite binary quadratic forms (A,B,C) with B^2 - 4AC = D:
 * reduced iff |B| <= A <= C and B >= 0 when |B| == A or A == C. */
static long si_class_number_neg4c(long c) {
    long D = -4 * c, h = 0;
    long Amax = (long)floor(sqrt((double)(4 * c) / 3.0)) + 1;
    for (long A = 1; A <= Amax; A++)
        for (long B = -A + 1; B <= A; B++) {
            long num = B * B - D;                     /* = 4AC */
            if (num % (4 * A) != 0) continue;
            long C = num / (4 * A);
            if (C < A) continue;
            if (si_lgcd(si_lgcd(A, B), C) != 1) continue;    /* primitive */
            if ((labs(B) == A || A == C) && B < 0) continue; /* boundary: B >= 0 */
            h++;
        }
    return h;
}

/* The Z[sqrt k] cube descent below is sound AND complete exactly when:
 *   - k < 0 and |k| is squarefree (so D = -4|k| is a fundamental discriminant
 *     and the field is imaginary),
 *   - k = 2,3 (mod 4)  (ring of integers is Z[sqrt k]; unit group {+/-1}, and a
 *     mod-8 argument then forces x ODD and coprime to k, so the two factors of
 *     x^3 = (y - sqrt k)(y + sqrt k) are ALWAYS coprime -- no non-coprime case
 *     analysis is needed),
 *   - 3 does not divide the class number h  (so an ideal whose cube is principal
 *     is itself principal, and every point comes from y + sqrt k = (a+b sqrt k)^3).
 * Outside this set (k = 1 mod 4 -> half-integer ring; 3 | h; k > 0 -> real field
 * with infinite units) the descent is declined, never guessed. */
static bool si_mordell_descent_sound(long k) {
    if (k >= 0 || !si_is_squarefree_l(k)) return false;
    int r = (int)(((k % 4) + 4) % 4);
    if (r != 2 && r != 3) return false;
    return si_class_number_neg4c(-k) % 3 != 0;
}

/* --- Unbounded Mordell  y^2 == x^3 + k  via factorisation in Z[sqrt k]. ---
 *
 * x^3 = y^2 - k = (y - sqrt k)(y + sqrt k).  For the sound set of k (see
 * si_mordell_descent_sound) the factors are coprime, the class number is prime
 * to 3, and the units are {+/-1}, so every point satisfies
 *   y + sqrt k = (a + b sqrt k)^3,   1 = 3 a^2 b + k b^3 = b (3 a^2 + k b^2),
 * forcing b = +/-1 and a fixed, giving the COMPLETE integer-point set
 *   x = a^2 - k b^2,   y = +/- (a^3 + 3 a k b^2).
 * So y^2 = x^3 - 2 -> (3, +/-5), y^2 = x^3 - 1 -> (1, 0), and y^2 = x^3 - 5 ->
 * {} (proved, not merely unfound).  Only engaged when x is unbounded (a bounded
 * box is handled by the ordinary leaf search). */
static Expr* si_solve_mordell(SICtx* c) {
    if (c->neq != 1) return NULL;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return NULL;

    /* Identify Y (only Y^2) and X (only X^3), monic and opposite sign, plus a
     * constant.  Normalise to  Y^2 - X^3 - k == 0  (so k = -const). */
    for (int choice = 0; choice < 2; choice++) {
        int Yv = act[choice], Xv = act[1 - choice];
        mpz_t cY, cX, c0; mpz_init_set_ui(cY, 0); mpz_init_set_ui(cX, 0); mpz_init_set_ui(c0, 0);
        bool shape = true;
        for (size_t t = 0; t < eq->n_terms && shape; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            for (int v = 0; v < n; v++) if (v != Yv && v != Xv && ex[v] != 0) shape = false;
            if (!shape) break;
            int ey = ex[Yv], exx = ex[Xv];
            if (ey == 2 && exx == 0) mpz_add(cY, cY, eq->coefs[t]);
            else if (ey == 0 && exx == 3) mpz_add(cX, cX, eq->coefs[t]);
            else if (ey == 0 && exx == 0) mpz_add(c0, c0, eq->coefs[t]);
            else shape = false;
        }
        bool ok = shape && (mpz_cmpabs_ui(cY, 1) == 0) && (mpz_cmpabs_ui(cX, 1) == 0)
               && (mpz_sgn(cY) == -mpz_sgn(cX));
        long k = 0;
        if (ok) {
            if (mpz_sgn(cY) < 0) { mpz_neg(cY, cY); mpz_neg(cX, cX); mpz_neg(c0, c0); }
            if (mpz_fits_slong_p(c0)) k = -mpz_get_si(c0);      /* y^2 = x^3 + k */
            else ok = false;
        }
        /* Engage on the sound descent set (any imaginary k with k = 2,3 mod 4,
         * |k| squarefree, 3 not dividing the class number), and only when X is
         * unbounded (a bounded box is handled by the ordinary leaf search). */
        bool unbounded = !(c->has_lo[Xv] && c->has_hi[Xv]);
        if (ok && si_mordell_descent_sound(k) && unbounded) {
            Expr** sols = NULL; int nsol = 0, cap = 0;
            for (int b = 1; b >= -1; b -= 2) {
                /* b (3 a^2 + k b^2) = 1  ->  3 a^2 = (1/b) - k b^2 = b - k b^2. */
                long rhs = b - k * b * b;                /* = 3 a^2 */
                if (rhs < 0 || rhs % 3 != 0) continue;
                long a2 = rhs / 3;
                long a0 = (long)llround(sqrt((double)a2));
                for (long a = a0 - 1; a <= a0 + 1; a++) {  /* pin the exact root */
                    if (a < 0 || a * a != a2) continue;
                    for (int as = (a == 0 ? 1 : 1); as >= (a == 0 ? 1 : -1); as -= 2) {
                        long av = as * a;
                        long xv = av * av - k * b * b;             /* norm */
                        long yv = av * av * av + 3 * av * k * b * b;
                        for (int ys = 1; ys >= -1; ys -= 2) {      /* unit +/-1 -> +/- y */
                            int64_t vals[SI_MAX_VARS];
                            for (int i = 0; i < n; i++) vals[i] = 0;
                            vals[Xv] = xv; vals[Yv] = ys * yv;
                            bool dup = false;
                            if (si_verify_symbolic(c->original, c->var, vals, n)) {
                                Expr* cand = si_exp_one_tuple(c->var, vals, n);
                                for (int u = 0; u < nsol && !dup; u++) if (expr_eq(sols[u], cand)) dup = true;
                                if (dup) { expr_free(cand); }
                                else { if (nsol == cap) { cap = cap ? cap * 2 : 4; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); } sols[nsol++] = cand; }
                            }
                        }
                    }
                }
            }
            mpz_clear(cY); mpz_clear(cX); mpz_clear(c0);
            /* sort ascending by (x, y) for a canonical result */
            for (int i = 0; i + 1 < nsol; i++) for (int j = i + 1; j < nsol; j++)
                if (expr_compare(sols[j], sols[i]) < 0) { Expr* tt = sols[i]; sols[i] = sols[j]; sols[j] = tt; }
            Expr* result = mk_list(sols, (size_t)nsol);
            free(sols);
            return result;
        }
        mpz_clear(cY); mpz_clear(cX); mpz_clear(c0);
    }
    return NULL;
}

/* --- Linear Diophantine: parametric (unbounded) family. --- */

/* Detect a single linear equation  sum a_i x_i == b.  Fills a[i] (coef of
 * x_i) and b (= -constant term); false if any term has degree > 1. */
static bool si_linear_detect(const MPoly* eq, int n, mpz_t* a, mpz_t b) {
    for (int i = 0; i < n; i++) mpz_set_ui(a[i], 0);
    mpz_set_ui(b, 0);
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nz = -1, deg = 0;
        for (int v = 0; v < n; v++) { deg += ex[v]; if (ex[v] > 0) { if (nz >= 0) return false; nz = v; } }
        if (deg > 1) return false;
        if (nz < 0) mpz_neg(b, eq->coefs[t]);           /* constant -> b = -const */
        else mpz_set(a[nz], eq->coefs[t]);
    }
    return true;
}

/* Particular solution x0[n] and homogeneous-lattice basis (n-1 vectors,
 * basis[j][i] = component i of vector j) of  sum a_i x_i == b, via the gcd
 * staircase.  Returns false (no solution) when gcd(a) does not divide b.  All
 * mpz arrays are pre-init'd by the caller. */
static bool si_linear_lattice(const mpz_t* a, int n, const mpz_t b,
                              mpz_t* x0, mpz_t** basis) {
    mpz_t g, s, t, ng, q1, q2, tmp;
    mpz_init(g); mpz_init(s); mpz_init(t); mpz_init(ng);
    mpz_init(q1); mpz_init(q2); mpz_init(tmp);
    mpz_t* r = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) mpz_init_set_ui(r[i], 0);

    mpz_set(g, a[0]); mpz_set_ui(r[0], 1);              /* g_0 = a_0, r_0 = e_0 */
    for (int i = 1; i < n; i++) {
        mpz_gcdext(ng, s, t, g, a[i]);                  /* ng = s*g + t*a_i */
        mpz_t* d = basis[i - 1];                        /* d = (a_i/ng) r - (g/ng) e_i */
        if (mpz_sgn(ng) != 0) { mpz_divexact(q1, a[i], ng); mpz_divexact(q2, g, ng); }
        else { mpz_set_ui(q1, 0); mpz_set_ui(q2, 0); }
        for (int k = 0; k < n; k++) mpz_mul(d[k], q1, r[k]);
        mpz_sub(d[i], d[i], q2);
        for (int k = 0; k < n; k++) mpz_mul(r[k], s, r[k]);   /* r = s*r + t*e_i */
        mpz_add(r[i], r[i], t);
        mpz_set(g, ng);
    }

    bool solvable;
    if (mpz_sgn(g) == 0) {                               /* all coefficients zero */
        solvable = (mpz_sgn(b) == 0);
        for (int i = 0; i < n; i++) mpz_set_ui(x0[i], 0);
    } else if (mpz_divisible_p(b, g)) {
        solvable = true;
        mpz_divexact(tmp, b, g);
        for (int i = 0; i < n; i++) mpz_mul(x0[i], tmp, r[i]);   /* x0 = (b/g) r */
    } else solvable = false;

    for (int i = 0; i < n; i++) mpz_clear(r[i]);
    free(r);
    mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(ng);
    mpz_clear(q1); mpz_clear(q2); mpz_clear(tmp);
    return solvable;
}

/* A single linear equation with NO other constraints -> the full parametric
 * family  {{x_i -> x0_i + sum_j basis[j][i] C[j+1]}}, C[k] integer parameters.
 * Returns the owned result Expr, {} for an unsolvable equation, or NULL to
 * decline (not this shape). */
static Expr* si_solve_linear_parametric(SICtx* c) {
    if (c->neq != 1 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < c->n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    int n = c->n;
    mpz_t* a = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t* x0 = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) { mpz_init(a[i]); mpz_init(x0[i]); }
    mpz_t b; mpz_init(b);
    mpz_t** basis = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)(n > 1 ? n - 1 : 1));
    for (int j = 0; j < n - 1; j++) {
        basis[j] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init_set_ui(basis[j][i], 0);
    }

    Expr* result = NULL;
    if (!si_linear_detect(c->eq[0], n, a, b)) {
        result = NULL;
    } else if (!si_linear_lattice(a, n, b, x0, basis)) {
        result = mk_list(NULL, 0);                       /* provably no solution */
    } else {
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
            int nt = 0;
            if (mpz_sgn(x0[i]) != 0) terms[nt++] = mk_mpz(x0[i]);
            for (int j = 0; j < n - 1; j++) {
                if (mpz_sgn(basis[j][i]) == 0) continue;
                Expr* Ck = expr_new_function(mk_sym("C"), (Expr*[]){ mk_int(j + 1) }, 1);
                terms[nt++] = mk_fn2("Times", mk_mpz(basis[j][i]), Ck);
            }
            Expr* val = (nt == 0) ? mk_int(0)
                      : (nt == 1) ? terms[0]
                      : expr_new_function(mk_sym("Plus"), terms, (size_t)nt);
            free(terms);
            val = eval_and_free(val);
            rules[i] = mk_rule(expr_copy(c->var[i]), val);
        }
        Expr* tuple = mk_list(rules, (size_t)n);
        free(rules);
        result = mk_list((Expr*[]){ tuple }, 1);
    }

    for (int j = 0; j < n - 1; j++) { for (int i = 0; i < n; i++) mpz_clear(basis[j][i]); free(basis[j]); }
    free(basis);
    for (int i = 0; i < n; i++) { mpz_clear(a[i]); mpz_clear(x0[i]); }
    free(a); free(x0); mpz_clear(b);
    return result;
}

/* --- General linear Diophantine SYSTEM (m >= 2 equations) via HNF. --- */

/* Solve the unconstrained integer linear system  A x = b  (m = c->neq linear
 * equations in n = c->n unknowns) completely, using the Hermite normal form.
 *
 * With  P A^T = R  (row HNF of A^T, P unimodular n x n, R n x m, rank rr), set
 * U = P^T and H = R^T so that  A U = H  and H is in column-echelon form.  The
 * substitution x = U y turns A x = b into H y = b; forward-substitution over
 * the rr pivot columns yields the pivot coordinates y_0..y_{rr-1} (each step an
 * exact-division / divisibility test -- a failure proves NO integer solution),
 * the free coordinates y_{rr}..y_{n-1} parameterise the kernel, and a check on
 * the non-pivot rows detects inconsistency.  The particular solution is
 * x0 = U y|_{free = 0}; the kernel basis is the free columns of U (= the free
 * rows of P).  Emits the family  {{x_i -> x0_i + sum_k Ker[k][i] C[k+1]}}, the
 * empty List for a provably unsolvable system, or NULL to decline (not the
 * shape: fewer than two equations, any bound / ordering / disequation present,
 * or a non-linear equation).  Single equations keep the dedicated
 * si_solve_linear_parametric path. */
static Expr* si_solve_linear_system_hnf(SICtx* c) {
    int n = c->n, m = c->neq;
    if (m < 2 || c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured) return NULL;
    for (int i = 0; i < n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    /* Build A (m x n) and b (m) from the linear equations; decline on any
     * nonlinear equation (si_linear_detect returns false). */
    mpz_t* A  = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m * (size_t)n);
    mpz_t* bv = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)m);
    for (int i = 0; i < m * n; i++) mpz_init(A[i]);
    for (int i = 0; i < m; i++) mpz_init(bv[i]);
    mpz_t* arow = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) mpz_init(arow[i]);
    mpz_t brow; mpz_init(brow);

    bool shape_ok = true;
    for (int q = 0; q < m && shape_ok; q++) {
        if (!si_linear_detect(c->eq[q], n, arow, brow)) { shape_ok = false; break; }
        for (int i = 0; i < n; i++) mpz_set(A[(size_t)q * n + i], arow[i]);
        mpz_set(bv[q], brow);
    }
    for (int i = 0; i < n; i++) mpz_clear(arow[i]);
    free(arow); mpz_clear(brow);

    if (!shape_ok) {
        for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
        for (int i = 0; i < m; i++) mpz_clear(bv[i]);
        free(A); free(bv);
        return NULL;
    }

    /* AT = A^T  (n x m). */
    mpz_t* AT = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n * (size_t)m);
    for (int i = 0; i < n * m; i++) mpz_init(AT[i]);
    for (int q = 0; q < m; q++)
        for (int i = 0; i < n; i++)
            mpz_set(AT[(size_t)i * m + q], A[(size_t)q * n + i]);

    mpz_t* R = NULL; mpz_t* P = NULL;
    int rr = linalg_hnf(AT, n, m, &R, &P);            /* P (n x n) * AT = R (n x m) */
    for (int i = 0; i < n * m; i++) mpz_clear(AT[i]);
    free(AT);
    if (rr < 0) {
        for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
        for (int i = 0; i < m; i++) mpz_clear(bv[i]);
        free(A); free(bv);
        linalg_hnf_free(R, n * m); linalg_hnf_free(P, n * n);
        return NULL;
    }

    /* Pivot column of each HNF row k (k = 0..rr-1): the H-row where pivot
     * column k of H = R^T leads.  R row k's first nonzero column. */
    int* pivrow = (int*)malloc(sizeof(int) * (size_t)(rr > 0 ? rr : 1));
    for (int k = 0; k < rr; k++) {
        pivrow[k] = -1;
        for (int j = 0; j < m; j++)
            if (mpz_sgn(R[(size_t)k * m + j]) != 0) { pivrow[k] = j; break; }
    }

    /* Forward substitution:  y_k = ( b[pivrow[k]] - sum_{l<k} H[pivrow[k]][l] y_l )
     *                               / H[pivrow[k]][k],   H[i][l] = R[l][i]. */
    mpz_t* y = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)(rr > 0 ? rr : 1));
    for (int k = 0; k < rr; k++) mpz_init(y[k]);
    mpz_t acc, prod; mpz_init(acc); mpz_init(prod);
    bool solvable = true;
    for (int k = 0; k < rr && solvable; k++) {
        int pr = pivrow[k];
        mpz_set(acc, bv[pr]);
        for (int l = 0; l < k; l++) {
            mpz_mul(prod, R[(size_t)l * m + pr], y[l]);   /* H[pr][l] = R[l][pr] */
            mpz_sub(acc, acc, prod);
        }
        const mpz_t* piv = &R[(size_t)k * m + pr];        /* H[pr][k] = R[k][pr] > 0 */
        if (!mpz_divisible_p(acc, *piv)) { solvable = false; break; }
        mpz_divexact(y[k], acc, *piv);
    }

    /* Consistency over every row i: sum_{l<rr} H[i][l] y_l == b[i]. */
    for (int i = 0; i < m && solvable; i++) {
        mpz_set_ui(acc, 0);
        for (int l = 0; l < rr; l++) {
            mpz_mul(prod, R[(size_t)l * m + i], y[l]);
            mpz_add(acc, acc, prod);
        }
        if (mpz_cmp(acc, bv[i]) != 0) solvable = false;
    }

    Expr* result;
    if (!solvable) {
        result = mk_list(NULL, 0);                        /* provably no solution */
    } else {
        /* x0_i = sum_{k<rr} U[i][k] y_k = sum_k P[k][i] y_k.
         * Kernel vector for free coordinate f (rr <= f < n): comp i = P[f][i]. */
        int nfree = n - rr;
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            mpz_t x0i; mpz_init(x0i);
            for (int k = 0; k < rr; k++) {
                mpz_mul(prod, P[(size_t)k * n + i], y[k]);
                mpz_add(x0i, x0i, prod);
            }
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(nfree + 1));
            int nt = 0;
            if (mpz_sgn(x0i) != 0) terms[nt++] = mk_mpz(x0i);
            for (int f = 0; f < nfree; f++) {
                const mpz_t* comp = &P[(size_t)(rr + f) * n + i];   /* P[rr+f][i] */
                if (mpz_sgn(*comp) == 0) continue;
                Expr* Ck = expr_new_function(mk_sym("C"),
                                             (Expr*[]){ mk_int(f + 1) }, 1);
                terms[nt++] = mk_fn2("Times", mk_mpz(*comp), Ck);
            }
            Expr* val = (nt == 0) ? mk_int(0)
                      : (nt == 1) ? terms[0]
                      : expr_new_function(mk_sym("Plus"), terms, (size_t)nt);
            free(terms);
            val = eval_and_free(val);
            rules[i] = mk_rule(expr_copy(c->var[i]), val);
            mpz_clear(x0i);
        }
        Expr* tuple = mk_list(rules, (size_t)n);
        free(rules);
        result = mk_list((Expr*[]){ tuple }, 1);
    }

    for (int k = 0; k < rr; k++) mpz_clear(y[k]);
    free(y); free(pivrow);
    mpz_clear(acc); mpz_clear(prod);
    for (int i = 0; i < m * n; i++) mpz_clear(A[i]);
    for (int i = 0; i < m; i++) mpz_clear(bv[i]);
    free(A); free(bv);
    linalg_hnf_free(R, n * m); linalg_hnf_free(P, n * n);
    return result;
}

/* Invert an m x m matrix G (double) into Ginv via Gauss-Jordan.  Returns false
 * if singular. */
static bool si_mat_inverse(double G[SI_MAX_VARS][SI_MAX_VARS], int m,
                           double Ginv[SI_MAX_VARS][SI_MAX_VARS]) {
    double A[SI_MAX_VARS][2 * SI_MAX_VARS];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) { A[i][j] = G[i][j]; A[i][m + j] = (i == j) ? 1.0 : 0.0; }
    for (int col = 0; col < m; col++) {
        int piv = col; double best = 0;
        for (int r = col; r < m; r++) {
            double av = A[r][col] < 0 ? -A[r][col] : A[r][col];
            if (av > best) { best = av; piv = r; }
        }
        if (best < 1e-9) return false;
        for (int j = 0; j < 2 * m; j++) { double tmp = A[col][j]; A[col][j] = A[piv][j]; A[piv][j] = tmp; }
        double d = A[col][col];
        for (int j = 0; j < 2 * m; j++) A[col][j] /= d;
        for (int r = 0; r < m; r++) if (r != col) {
            double f = A[r][col];
            for (int j = 0; j < 2 * m; j++) A[r][j] -= f * A[col][j];
        }
    }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) Ginv[i][j] = A[i][m + j];
    return true;
}

/* A single linear equation over a finite box.  gcd(a) not dividing b -> {}.
 * Otherwise the solution set is  x0 + Z-span of an (n-1)-vector lattice; the
 * lattice is LLL-reduced (LatticeReduce) and the coefficient box is the exact
 * projection of the value box under the pseudoinverse, so the search box is
 * small when the coefficients are large (few solutions) and provably too large
 * (declined) when the lattice is dense.  Every candidate is checked exactly.
 * Returns true when it settled the answer (emitting nothing means {}). */
static bool si_solve_linear_bounded(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int n = c->n, m = n - 1;
    if (m < 1 || n > SI_MAX_VARS) return false;
    for (int i = 0; i < n; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;

    mpz_t* a  = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    mpz_t* x0 = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
    for (int i = 0; i < n; i++) { mpz_init(a[i]); mpz_init(x0[i]); }
    mpz_t b, g; mpz_init(b); mpz_init_set_ui(g, 0);
    mpz_t** basis = (mpz_t**)malloc(sizeof(mpz_t*) * (size_t)m);
    for (int j = 0; j < m; j++) { basis[j] = (mpz_t*)malloc(sizeof(mpz_t) * (size_t)n);
        for (int i = 0; i < n; i++) mpz_init_set_ui(basis[j][i], 0); }

    bool handled = false;
    if (si_linear_detect(c->eq[0], n, a, b)) {
        for (int i = 0; i < n; i++) mpz_gcd(g, g, a[i]);
        bool solvable = (mpz_sgn(g) == 0) ? (mpz_sgn(b) == 0)
                                          : mpz_divisible_p(b, g);
        if (!solvable) {
            handled = true;                              /* provably {} */
        } else if (si_linear_lattice(a, n, b, x0, basis)) {
            /* LLL-reduce the homogeneous basis via LatticeReduce[...]. */
            Expr** rowsE = (Expr**)malloc(sizeof(Expr*) * (size_t)m);
            for (int j = 0; j < m; j++) {
                Expr** comp = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                for (int i = 0; i < n; i++) comp[i] = mk_mpz(basis[j][i]);
                rowsE[j] = mk_list(comp, (size_t)n); free(comp);
            }
            Expr* red = eval_and_free(expr_new_function(mk_sym("LatticeReduce"),
                                      (Expr*[]){ mk_list(rowsE, (size_t)m) }, 1));
            free(rowsE);

            mpz_t Bz[SI_MAX_VARS][SI_MAX_VARS];
            double Bd[SI_MAX_VARS][SI_MAX_VARS];
            for (int j = 0; j < SI_MAX_VARS; j++) for (int i = 0; i < SI_MAX_VARS; i++) mpz_init(Bz[j][i]);
            int mm = 0;
            if (red && red->type == EXPR_FUNCTION
                && red->data.function.head->type == EXPR_SYMBOL
                && red->data.function.head->data.symbol.name == SYM_List
                && (int)red->data.function.arg_count == m) {
                mm = m;
                for (int j = 0; j < m; j++) {
                    Expr* row = red->data.function.args[j];
                    for (int i = 0; i < n; i++) {
                        mpz_set(Bz[j][i], basis[j][i]);      /* fallback: unreduced */
                        if (row->type == EXPR_FUNCTION && (int)row->data.function.arg_count == n) {
                            int64_t vv;
                            if (expr_as_i64(row->data.function.args[i], &vv)) mpz_set_si(Bz[j][i], vv);
                        }
                        Bd[j][i] = mpz_get_d(Bz[j][i]);
                    }
                }
            }
            if (red) expr_free(red);

            /* P = (B B^T)^{-1} B  (m x n). */
            double Gm[SI_MAX_VARS][SI_MAX_VARS], Ginv[SI_MAX_VARS][SI_MAX_VARS];
            for (int j = 0; j < mm; j++) for (int k = 0; k < mm; k++) {
                double s = 0; for (int i = 0; i < n; i++) s += Bd[j][i] * Bd[k][i]; Gm[j][k] = s;
            }
            if (mm == m && si_mat_inverse(Gm, m, Ginv)) {
                /* Coefficient box = exact projection of the value box. */
                double Lb[SI_MAX_VARS], Ub[SI_MAX_VARS];
                for (int i = 0; i < n; i++) {
                    double x0d = mpz_get_d(x0[i]);
                    Lb[i] = (double)c->lo[i] - x0d;
                    Ub[i] = (double)c->hi[i] - x0d;
                }
                int64_t Clo[SI_MAX_VARS], Chi[SI_MAX_VARS];
                long double space = 1.0L; bool ok = true;
                for (int j = 0; j < m && ok; j++) {
                    double lo = 0, hi = 0;
                    for (int i = 0; i < n; i++) {
                        double P = 0; for (int k = 0; k < m; k++) P += Ginv[j][k] * Bd[k][i];
                        double t1 = P * Lb[i], t2 = P * Ub[i];
                        lo += (t1 < t2) ? t1 : t2;
                        hi += (t1 < t2) ? t2 : t1;
                    }
                    Clo[j] = (int64_t)floor(lo) - 1;
                    Chi[j] = (int64_t)ceil(hi) + 1;
                    space *= (long double)(Chi[j] - Clo[j] + 1);
                    if (space > (long double)SI_MAX_NODES) ok = false;
                }
                if (ok) {
                    st->max_visits = SI_MAX_NODES;
                    int64_t cc[SI_MAX_VARS];
                    for (int j = 0; j < m; j++) cc[j] = Clo[j];
                    mpz_t xi, term; mpz_init(xi); mpz_init(term);
                    for (;;) {
                        int64_t vals[SI_MAX_VARS]; bool inbox = true;
                        for (int i = 0; i < n && inbox; i++) {
                            mpz_set(xi, x0[i]);
                            for (int j = 0; j < m; j++) { mpz_mul_si(term, Bz[j][i], (long)cc[j]); mpz_add(xi, xi, term); }
                            if (!mpz_fits_slong_p(xi)) { inbox = false; break; }
                            int64_t xv = mpz_get_si(xi);
                            if (xv < c->lo[i] || xv > c->hi[i]) inbox = false;
                            vals[i] = xv;
                        }
                        if (inbox && si_verify(c, vals)) emit_full(st, vals);
                        int j = 0;
                        for (; j < m; j++) { if (++cc[j] <= Chi[j]) break; cc[j] = Clo[j]; }
                        if (j == m) break;
                    }
                    mpz_clear(xi); mpz_clear(term);
                    handled = true;
                }
            }
            for (int j = 0; j < SI_MAX_VARS; j++) for (int i = 0; i < SI_MAX_VARS; i++) mpz_clear(Bz[j][i]);
        }
    }

    for (int j = 0; j < m; j++) { for (int i = 0; i < n; i++) mpz_clear(basis[j][i]); free(basis[j]); }
    free(basis);
    for (int i = 0; i < n; i++) { mpz_clear(a[i]); mpz_clear(x0[i]); }
    free(a); free(x0); mpz_clear(b); mpz_clear(g);
    return handled;
}

/* --- Power-sum divisor method for separable additive equations. ---
 *
 * For an equation  coef*(x^e + y^e) + (terms in the other variables) == 0
 * with e ODD, enumerate the "outer" variables over their box; for each, the
 * inner pair satisfies  x^e + y^e == m.  Because e is odd, s = x + y divides m,
 * and the power sum p_e(s, p) (p = x*y) is a degree-e/2 polynomial in p, so for
 * each divisor s of m the integer roots p give (x, y) via t^2 - s t + p = 0.
 * This turns the O(N) inner loop into O(#divisors(m)) -- the difference between
 * O(N^2) brute force and O(N * factoring) for e.g. sums of three cubes. */

/* Solve x^e + y^e == m (e odd) for the inner pair (ip, jp), outer variables
 * already fixed in vals.  Emits verified full assignments. */
static void si_two_power_solve(SICtx* c, SearchState* st, int ip, int jp,
                               int e, const mpz_t m, int64_t* vals) {
    int64_t loi = c->lo[ip], hii = c->hi[ip], loj = c->lo[jp], hij = c->hi[jp];
    /* Tighten the pair windows by any abs-ordering that bounds a pair variable
     * below an already-assigned outer variable (|pair| < |vals[outer]|).  Sound
     * (prune-only; si_verify still backstops), and for the ordered three-cubes
     * query it shrinks the s = x+y range to ~2|outer| as the outer sweeps. */
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k]; int s = c->abs_ord_strict[k] ? 1 : 0;
        if (b == ip || b == jp) continue;            /* larger side must be an outer */
        if (a != ip && a != jp) continue;            /* smaller side must be a pair var */
        int64_t vb = vals[b] < 0 ? -vals[b] : vals[b];
        int64_t bnd = vb - s; if (bnd < 0) bnd = -1; /* |a| <= bnd */
        if (a == ip) { if ( bnd < hii) hii =  bnd; if (-bnd > loi) loi = -bnd; }
        else         { if ( bnd < hij) hij =  bnd; if (-bnd > loj) loj = -bnd; }
    }
    if (mpz_sgn(m) == 0) {                        /* x^e + y^e = 0  ->  y = -x */
        int64_t t0 = (loi > -hij) ? loi : -hij;
        int64_t t1 = (hii < -loj) ? hii : -loj;
        for (int64_t t = t0; t <= t1 && !st->overflow; t++) {
            if (++st->visits > st->max_visits) { st->overflow = true; return; }
            vals[ip] = t; vals[jp] = -t; if (si_verify(c, vals)) emit_full(st, vals);
        }
        return;
    }
    int deg = e / 2;                              /* degree in p = x*y */
    if (deg > SI_LEAF_MAXDEG) return;
    mpz_t absm; mpz_init(absm); mpz_abs(absm, m);
    Expr* dl = divisors_ordinary(absm);
    mpz_clear(absm);
    if (!dl || dl->type != EXPR_FUNCTION) { if (dl) expr_free(dl); return; }

    /* Newton-recurrence scratch for p_k(s, p) as a polynomial in p:
     * p_0 = 2, p_1 = s, p_k = s*p_{k-1} - p*p_{k-2}. */
    mpz_t sz, pm2[SI_LEAF_MAXDEG + 1], pm1[SI_LEAF_MAXDEG + 1], cur[SI_LEAF_MAXDEG + 1];
    mpz_t a[SI_LEAF_MAXDEG + 1], sq, r, num, tmp;
    mpz_init(sz); mpz_init(sq); mpz_init(r); mpz_init(num); mpz_init(tmp);
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) { mpz_init(pm2[k]); mpz_init(pm1[k]); mpz_init(cur[k]); mpz_init(a[k]); }

    /* |s| = |x+y| is at most the box span. */
    int64_t smax = (hii - loj);
    int64_t smin = (loi - hij);
    int64_t sbound = (smax > -smin) ? smax : -smin;

    for (size_t di = 0; di < dl->data.function.arg_count && !st->overflow; di++) {
        int64_t dv;
        if (!expr_as_i64(dl->data.function.args[di], &dv)) continue;
        for (int sgn = 1; sgn >= -1; sgn -= 2) {
            int64_t sval = sgn * dv;
            if (sval < smin || sval > smax) continue;      /* s = x+y out of range */
            if (++st->visits > st->max_visits) { st->overflow = true; break; }
            mpz_set_si(sz, sval);
            /* Build p_e(s, p) via the recurrence (polynomials in p). */
            for (int k = 0; k <= deg; k++) { mpz_set_ui(pm2[k], 0); mpz_set_ui(pm1[k], 0); }
            mpz_set_ui(pm2[0], 2);                          /* p_0 = 2 */
            mpz_set(pm1[0], sz);                            /* p_1 = s */
            int dm2 = 0, dm1 = 0;
            for (int k = 2; k <= e; k++) {
                for (int t = 0; t <= deg; t++) mpz_set_ui(cur[t], 0);
                for (int t = 0; t <= dm1; t++) mpz_mul(cur[t], sz, pm1[t]);        /* s*p_{k-1} */
                for (int t = 0; t <= dm2; t++) mpz_sub(cur[t + 1], cur[t + 1], pm2[t]); /* - p*p_{k-2} */
                int dc = dm1; if (dm2 + 1 > dc) dc = dm2 + 1;
                for (int t = 0; t <= deg; t++) { mpz_set(pm2[t], pm1[t]); mpz_set(pm1[t], cur[t]); }
                dm2 = dm1; dm1 = dc;
            }
            for (int t = 0; t <= deg; t++) mpz_set(a[t], pm1[t]);
            mpz_sub(a[0], a[0], m);                         /* p_e(s,p) - m == 0 */

            int64_t proots[SI_LEAF_MAXDEG * 2];
            int npr = univariate_roots(a, deg, INT64_MIN / 2, INT64_MAX / 2, proots);
            for (int pi = 0; pi < npr; pi++) {
                /* x, y from t^2 - s t + p = 0: disc = s^2 - 4p. */
                mpz_mul(sq, sz, sz);
                mpz_set_si(tmp, proots[pi]); mpz_mul_ui(tmp, tmp, 4); mpz_sub(sq, sq, tmp);
                if (mpz_sgn(sq) < 0 || !mpz_perfect_square_p(sq)) continue;
                mpz_sqrt(r, sq);
                for (int rs = 1; rs >= -1; rs -= 2) {
                    mpz_set(num, sz); if (rs < 0) mpz_sub(num, num, r); else mpz_add(num, num, r);
                    if (!mpz_divisible_ui_p(num, 2)) continue;
                    mpz_divexact_ui(num, num, 2);           /* x = (s +/- r)/2 */
                    if (!mpz_fits_slong_p(num)) continue;
                    int64_t xv = mpz_get_si(num);
                    int64_t yv = sval - xv;
                    if (xv < loi || xv > hii || yv < loj || yv > hij) continue;
                    vals[ip] = xv; vals[jp] = yv;
                    if (si_verify(c, vals)) emit_full(st, vals);
                }
            }
            (void)sbound;
        }
    }
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) { mpz_clear(pm2[k]); mpz_clear(pm1[k]); mpz_clear(cur[k]); mpz_clear(a[k]); }
    mpz_clear(sz); mpz_clear(sq); mpz_clear(r); mpz_clear(num); mpz_clear(tmp);
    expr_free(dl);
}

/* Detect a separable additive equation with two variables sharing an odd
 * exponent, and solve it by the divisor method (outer variables enumerated,
 * inner pair solved per fixed outer assignment).  Returns true if handled. */
static bool si_solve_powersum_divisor(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0];
    int n = c->n;
    for (int i = 0; i < n; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;

    /* Separable: every term touches at most one variable.  Record each
     * variable's single power term (exponent, coefficient), if it is a lone
     * power c*v^e. */
    int vexp[SI_MAX_VARS]; mpz_t vcoef[SI_MAX_VARS];
    int nterm[SI_MAX_VARS];
    for (int i = 0; i < n; i++) { vexp[i] = -1; nterm[i] = 0; mpz_init_set_ui(vcoef[i], 0); }
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nz = -1;
        for (int v = 0; v < n; v++) if (ex[v] > 0) { if (nz >= 0) ok = false; nz = v; }
        if (!ok) break;
        if (nz < 0) continue;                        /* constant term */
        nterm[nz]++; vexp[nz] = ex[nz]; mpz_set(vcoef[nz], eq->coefs[t]);
    }
    int ip = -1, jp = -1;
    if (ok) {
        for (int i = 0; i < n && ip < 0; i++) {
            if (nterm[i] != 1 || vexp[i] < 3 || (vexp[i] % 2) == 0) continue;
            for (int j = i + 1; j < n; j++) {
                if (nterm[j] == 1 && vexp[j] == vexp[i] && mpz_cmp(vcoef[j], vcoef[i]) == 0)
                    { ip = i; jp = j; break; }
            }
        }
    }

    bool handled = false;
    if (ip >= 0) {
        int e = vexp[ip];
        /* Outer variables enumeration size. */
        long double outer = 1.0L;
        for (int i = 0; i < n; i++) if (i != ip && i != jp)
            outer *= (long double)(c->hi[i] - c->lo[i] + 1);
        /* |m| = |x^e + y^e| ~ B^e is the number the method factors once per outer
         * assignment; keep it in the fast-factoring range so the O(N*factoring)
         * budget stays reasonable (this admits cubes over |.|<10^6, and higher
         * powers only over correspondingly smaller boxes). */
        int64_t Bin = 0;
        for (int q = 0; q < 2; q++) {
            int vv = (q == 0) ? ip : jp;
            int64_t b = (c->hi[vv] > -c->lo[vv]) ? c->hi[vv] : -c->lo[vv];
            if (b > Bin) Bin = b;
        }
        long double mmag = 1.0L;
        for (int k = 0; k < e; k++) mmag *= (long double)Bin;
        if (outer <= 3.0e5L && mmag <= 1.0e18L) {    /* factoring budget bounded */
            handled = true;
            st->max_visits = SI_MAX_NODES;
            /* Odometer over the outer variables. */
            int outv[SI_MAX_VARS], nouter = 0;
            for (int i = 0; i < n; i++) if (i != ip && i != jp) outv[nouter++] = i;
            int64_t vals[SI_MAX_VARS];
            for (int i = 0; i < n; i++) vals[i] = 0;
            for (int q = 0; q < nouter; q++) vals[outv[q]] = c->lo[outv[q]];
            mpz_t rest, mm, coef; mpz_init(rest); mpz_init(mm); mpz_init_set(coef, vcoef[ip]);
            for (;;) {
                vals[ip] = 0; vals[jp] = 0;
                si_eval_mpoly(eq, vals, rest);       /* = sum of other terms + const */
                mpz_neg(mm, rest);
                if (mpz_divisible_p(mm, coef)) {
                    mpz_divexact(mm, mm, coef);      /* x^e + y^e == mm */
                    si_two_power_solve(c, st, ip, jp, e, mm, vals);
                }
                int q = 0;
                for (; q < nouter; q++) {
                    if (++vals[outv[q]] <= c->hi[outv[q]]) break;
                    vals[outv[q]] = c->lo[outv[q]];
                }
                if (q == nouter || st->overflow) break;
            }
            mpz_clear(rest); mpz_clear(mm); mpz_clear(coef);
        }
    }
    for (int i = 0; i < n; i++) mpz_clear(vcoef[i]);
    return handled;
}

/* --- Conic  Y^2 == A X^2 + B X + C  with A a perfect square. ---
 *
 * When one variable Y appears only as Y^2 (coefficient +/-1) and the other, X,
 * to degree <= 2 with no cross term, and the X^2 coefficient A is a POSITIVE
 * PERFECT SQUARE p^2, the curve is a "parabola-like" hyperbola with finitely
 * many integer points.  Completing the square,
 *   (2 p Y)^2 - (2 A X + B)^2 = 4 A C - B^2 =: D,
 * so every solution comes from a factorisation  d1 * d2 = D  via
 *   Y = (d1 + d2) / (4 p),  X = (d2 - d1 - 2 B) / (4 A).
 * Enumerating the divisors of |D| is exhaustive, so an empty result is a proof
 * of no solutions.  Euler's  n^2 + n + 41 == y^2  has D = 163 (prime) and the
 * unique positive point n = 40, y = 41.  (A non-square A is a genuine Pell
 * conic -- left to the Pell continued-fraction path.) */
static bool si_solve_conic(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;

    for (int choice = 0; choice < 2; choice++) {
        int Yv = act[choice], Xv = act[1 - choice];
        mpz_t e0, e2, e1, ec; mpz_init_set_ui(e0, 0); mpz_init_set_ui(e2, 0);
        mpz_init_set_ui(e1, 0); mpz_init_set_ui(ec, 0);
        bool shape = true;
        for (size_t t = 0; t < eq->n_terms && shape; t++) {
            const int* ex = eq->exps + t * (size_t)n;
            for (int v = 0; v < n; v++) if (v != Yv && v != Xv && ex[v] != 0) shape = false;
            if (!shape) break;
            int ey = ex[Yv], exx = ex[Xv];
            if (ey == 2 && exx == 0) mpz_add(e0, e0, eq->coefs[t]);
            else if (ey == 0 && exx == 2) mpz_add(e2, e2, eq->coefs[t]);
            else if (ey == 0 && exx == 1) mpz_add(e1, e1, eq->coefs[t]);
            else if (ey == 0 && exx == 0) mpz_add(ec, ec, eq->coefs[t]);
            else shape = false;                       /* Y linear, cross, or higher */
        }
        bool ok = shape && (mpz_cmpabs_ui(e0, 1) == 0);   /* Y^2 coefficient +/-1 */
        if (ok) {
            if (mpz_sgn(e0) < 0) { mpz_neg(e0, e0); mpz_neg(e2, e2); mpz_neg(e1, e1); mpz_neg(ec, ec); }
            /* Y^2 = A X^2 + B X + C. */
            mpz_t A, B, C, D, p, tmp; mpz_init(A); mpz_init(B); mpz_init(C);
            mpz_init(D); mpz_init(p); mpz_init(tmp);
            mpz_neg(A, e2); mpz_neg(B, e1); mpz_neg(C, ec);
            bool go = (mpz_sgn(A) > 0) && mpz_perfect_square_p(A);
            if (go) {
                mpz_sqrt(p, A);
                mpz_mul(D, A, C); mpz_mul_ui(D, D, 4);
                mpz_mul(tmp, B, B); mpz_sub(D, D, tmp);   /* D = 4AC - B^2 */
            }
            if (go && mpz_sgn(D) != 0) {
                st->max_visits = SI_MAX_NODES;
                mpz_t absD; mpz_init(absD); mpz_abs(absD, D);
                Expr* dl = divisors_ordinary(absD);
                mpz_clear(absD);
                mpz_t d1, d2, Y, X, fourp, fourA, twoB, num;
                mpz_init(d1); mpz_init(d2); mpz_init(Y); mpz_init(X);
                mpz_init(fourp); mpz_init(fourA); mpz_init(twoB); mpz_init(num);
                mpz_mul_ui(fourp, p, 4); mpz_mul_ui(fourA, A, 4); mpz_mul_ui(twoB, B, 2);
                if (dl && dl->type == EXPR_FUNCTION) {
                    for (size_t i = 0; i < dl->data.function.arg_count; i++) {
                        if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                        mpz_t g; mpz_init(g); expr_to_mpz(dl->data.function.args[i], g);
                        for (int sg = 1; sg >= -1; sg -= 2) {
                            mpz_set(d1, g); if (sg < 0) mpz_neg(d1, d1);   /* d1 | D */
                            if (!mpz_divisible_p(D, d1)) continue;
                            mpz_divexact(d2, D, d1);
                            mpz_add(num, d1, d2);                          /* Y = (d1+d2)/(4p) */
                            if (!mpz_divisible_p(num, fourp)) continue;
                            mpz_divexact(Y, num, fourp);
                            mpz_sub(num, d2, d1); mpz_sub(num, num, twoB);  /* X = (d2-d1-2B)/(4A) */
                            if (!mpz_divisible_p(num, fourA)) continue;
                            mpz_divexact(X, num, fourA);
                            if (!mpz_fits_slong_p(Y) || !mpz_fits_slong_p(X)) continue;
                            int64_t vals[SI_MAX_VARS];
                            for (int k = 0; k < n; k++) vals[k] = 0;
                            vals[Yv] = mpz_get_si(Y); vals[Xv] = mpz_get_si(X);
                            if (si_verify(c, vals)) emit_full(st, vals);
                        }
                        mpz_clear(g);
                    }
                }
                if (dl) expr_free(dl);
                mpz_clear(d1); mpz_clear(d2); mpz_clear(Y); mpz_clear(X);
                mpz_clear(fourp); mpz_clear(fourA); mpz_clear(twoB); mpz_clear(num);
                mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(D); mpz_clear(p); mpz_clear(tmp);
                mpz_clear(e0); mpz_clear(e2); mpz_clear(e1); mpz_clear(ec);
                return true;                                   /* exhaustive: handled */
            }
            mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(D); mpz_clear(p); mpz_clear(tmp);
        }
        mpz_clear(e0); mpz_clear(e2); mpz_clear(e1); mpz_clear(ec);
    }
    return false;
}

/* --- Factorable binary quadratic  A x^2 + B x y + C y^2 + D x + E y + F == 0
 *     with a CROSS term and perfect-square discriminant delta = B^2 - 4 A C. ---
 *
 * When delta = k^2 > 0 the quadratic part factors into two distinct rational
 * linear forms, so the conic is a hyperbola with rational asymptotes and only
 * finitely many integer points -- Runge's simplest case.  Completing the square
 * (A != 0; swap x,y when A == 0) gives, with U = 2 A x + B y + D,
 *   4 A * eq = U^2 - k^2 y^2 + P y + Q,   P = 4 A E - 2 B D,  Q = 4 A F - D^2,
 * and multiplying by 4 k^2 with V = 2 k^2 y - P collapses it to a difference of
 * squares
 *   (2 k U)^2 - V^2 = W,   W = -(P^2 + 4 k^2 Q),
 * so every solution comes from a factorisation d1 * d2 = W via
 *   U = (d1 + d2) / (4 k),  V = (d2 - d1) / 2,
 *   y = (V + P) / (2 k^2),  x = (U - B y - D) / (2 A).
 * Enumerating the divisors of |W| is exhaustive, so an empty result is a PROOF
 * of no solutions (e.g. (x - y)(x + 2 y) == 15 -> {} by a mod-3 obstruction).
 * Generalises si_solve_conic to the cross-term case and to non-unit square
 * coefficients; the no-cross-term / Pell forms are handled before this.  W == 0
 * (two parallel/coincident lines: a positive-dimensional union) is declined. */
static bool si_solve_factorable_conic(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;

    int Xv = act[0], Yv = act[1];
    mpz_t A, B, C, D, E, F;
    mpz_inits(A, B, C, D, E, F, NULL);
    bool shape = true;
    for (size_t t = 0; t < eq->n_terms && shape; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++)
            if (v != Xv && v != Yv && ex[v] != 0) { shape = false; break; }
        if (!shape) break;
        int ix = ex[Xv], iy = ex[Yv];
        if (ix == 2 && iy == 0) mpz_add(A, A, eq->coefs[t]);
        else if (ix == 1 && iy == 1) mpz_add(B, B, eq->coefs[t]);
        else if (ix == 0 && iy == 2) mpz_add(C, C, eq->coefs[t]);
        else if (ix == 1 && iy == 0) mpz_add(D, D, eq->coefs[t]);
        else if (ix == 0 && iy == 1) mpz_add(E, E, eq->coefs[t]);
        else if (ix == 0 && iy == 0) mpz_add(F, F, eq->coefs[t]);
        else shape = false;                               /* total degree > 2 */
    }

    /* Need A != 0; swap x <-> y when only C is nonzero.  Both zero -> bilinear
     * (handled elsewhere), decline. */
    if (shape && mpz_sgn(A) == 0 && mpz_sgn(C) != 0) {
        mpz_swap(A, C); mpz_swap(D, E); int tmp = Xv; Xv = Yv; Yv = tmp;
    }
    bool go = shape && mpz_sgn(A) != 0;

    mpz_t delta, k, P, Q, W, tmp;
    mpz_inits(delta, k, P, Q, W, tmp, NULL);
    if (go) {
        mpz_mul(delta, B, B); mpz_mul(tmp, A, C); mpz_submul_ui(delta, tmp, 4); /* B^2 - 4AC */
        go = (mpz_sgn(delta) > 0) && mpz_perfect_square_p(delta);
    }
    if (go) {
        mpz_sqrt(k, delta);
        mpz_mul(P, A, E); mpz_mul_ui(P, P, 4);
        mpz_mul(tmp, B, D); mpz_submul_ui(P, tmp, 2);      /* P = 4AE - 2BD */
        mpz_mul(Q, A, F); mpz_mul_ui(Q, Q, 4);
        mpz_mul(tmp, D, D); mpz_sub(Q, Q, tmp);            /* Q = 4AF - D^2 */
        mpz_mul(W, P, P);                                  /* W = -(P^2 + 4 k^2 Q) */
        mpz_mul(tmp, k, k); mpz_mul(tmp, tmp, Q); mpz_addmul_ui(W, tmp, 4);
        mpz_neg(W, W);
        go = (mpz_sgn(W) != 0);
    }

    bool handled = false;
    if (go) {
        handled = true;
        st->max_visits = SI_MAX_NODES;
        mpz_t absW; mpz_init(absW); mpz_abs(absW, W);
        Expr* dl = divisors_ordinary(absW);
        mpz_clear(absW);
        mpz_t d1, d2, U, V, x, y, fourk, twok2, num;
        mpz_inits(d1, d2, U, V, x, y, fourk, twok2, num, NULL);
        mpz_mul_ui(fourk, k, 4);
        mpz_mul(twok2, k, k); mpz_mul_ui(twok2, twok2, 2);
        if (dl && dl->type == EXPR_FUNCTION) {
            for (size_t i = 0; i < dl->data.function.arg_count && !st->overflow; i++) {
                if (!expr_is_integer_like(dl->data.function.args[i])) continue;
                mpz_t g; mpz_init(g); expr_to_mpz(dl->data.function.args[i], g);
                for (int sg = 1; sg >= -1; sg -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(d1, g); if (sg < 0) mpz_neg(d1, d1);   /* d1 | W */
                    if (!mpz_divisible_p(W, d1)) continue;
                    mpz_divexact(d2, W, d1);
                    mpz_add(num, d1, d2);                          /* U = (d1+d2)/(4k) */
                    if (!mpz_divisible_p(num, fourk)) continue;
                    mpz_divexact(U, num, fourk);
                    mpz_sub(V, d2, d1);                            /* V = (d2 - d1)/2 */
                    if (mpz_odd_p(V)) continue;
                    mpz_divexact_ui(V, V, 2);
                    mpz_add(num, V, P);                            /* y = (V+P)/(2k^2) */
                    if (!mpz_divisible_p(num, twok2)) continue;
                    mpz_divexact(y, num, twok2);
                    mpz_mul(num, B, y); mpz_sub(num, U, num); mpz_sub(num, num, D); /* x=(U-By-D)/(2A) */
                    mpz_mul_ui(tmp, A, 2);
                    if (!mpz_divisible_p(num, tmp)) continue;
                    mpz_divexact(x, num, tmp);
                    if (!mpz_fits_slong_p(x) || !mpz_fits_slong_p(y)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int kk = 0; kk < n; kk++) vals[kk] = 0;
                    vals[Xv] = mpz_get_si(x); vals[Yv] = mpz_get_si(y);
                    if (si_verify(c, vals)) emit_full(st, vals);
                }
                mpz_clear(g);
            }
        }
        if (dl) expr_free(dl);
        mpz_clears(d1, d2, U, V, x, y, fourk, twok2, num, NULL);
    }

    mpz_clears(A, B, C, D, E, F, delta, k, P, Q, W, tmp, NULL);
    return handled;
}

/* --- Definite binary quadratic (ellipse)  A x^2 + B x y + C y^2 + D x + E y + F
 *     == 0  with discriminant  delta = B^2 - 4 A C < 0. ---
 *
 * A negative discriminant makes the quadratic part positive- (or negative-)
 * definite, so the conic is a compact ellipse with FINITELY many integer points
 * -- but a rotated one (B != 0) is not bounded by derive_bounds, which only
 * boxes sign-definite / even-only variables.  Treat the equation as a quadratic
 * in x for each fixed y:  A x^2 + (B y + D) x + (C y^2 + E y + F) == 0.  Its
 * x-discriminant  disc_x(y) = delta y^2 + (2 B D - 4 A E) y + (D^2 - 4 A F)  is a
 * DOWNWARD parabola in y (leading coefficient delta < 0), so real x exist only
 * for y in the finite interval between its roots.  Enumerate those y, solve the
 * integer quadratic in x exactly (perfect-square discriminant), and verify.
 * Exhaustive over the (bounded) ellipse, so an empty result is a PROOF. */
static bool si_solve_elliptic_bqf(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    const MPoly* eq = c->eq[0]; int n = c->n;
    int act[SI_MAX_VARS], ka = 0;
    for (int v = 0; v < n; v++) if (mpoly_deg_var(eq, v) >= 1) act[ka++] = v;
    if (ka != 2) return false;
    int Xv = act[0], Yv = act[1];

    mpz_t A, B, C, D, E, F; mpz_inits(A, B, C, D, E, F, NULL);
    bool shape = true;
    for (size_t t = 0; t < eq->n_terms && shape; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        for (int v = 0; v < n; v++)
            if (v != Xv && v != Yv && ex[v] != 0) { shape = false; break; }
        if (!shape) break;
        int ix = ex[Xv], iy = ex[Yv];
        if (ix == 2 && iy == 0) mpz_add(A, A, eq->coefs[t]);
        else if (ix == 1 && iy == 1) mpz_add(B, B, eq->coefs[t]);
        else if (ix == 0 && iy == 2) mpz_add(C, C, eq->coefs[t]);
        else if (ix == 1 && iy == 0) mpz_add(D, D, eq->coefs[t]);
        else if (ix == 0 && iy == 1) mpz_add(E, E, eq->coefs[t]);
        else if (ix == 0 && iy == 0) mpz_add(F, F, eq->coefs[t]);
        else shape = false;                                /* total degree > 2 */
    }

    mpz_t delta, tmp; mpz_inits(delta, tmp, NULL);
    if (shape) {
        mpz_mul(delta, B, B); mpz_mul(tmp, A, C); mpz_submul_ui(delta, tmp, 4);  /* B^2 - 4AC */
    }
    if (!shape || mpz_sgn(delta) >= 0) {                   /* not a definite ellipse */
        mpz_clears(A, B, C, D, E, F, delta, tmp, NULL); return false;
    }
    /* Normalise the leading coefficient positive (delta unchanged by negation). */
    if (mpz_sgn(A) < 0) {
        mpz_neg(A, A); mpz_neg(B, B); mpz_neg(C, C);
        mpz_neg(D, D); mpz_neg(E, E); mpz_neg(F, F);
    }

    /* y-range where disc_x(y) = delta y^2 + beta y + gamma >= 0. */
    mpz_t beta, gamma, inner; mpz_inits(beta, gamma, inner, NULL);
    mpz_mul(beta, B, D); mpz_mul_ui(beta, beta, 2);
    mpz_mul(tmp, A, E); mpz_submul_ui(beta, tmp, 4);       /* beta = 2BD - 4AE */
    mpz_mul(gamma, D, D); mpz_mul(tmp, A, F); mpz_submul_ui(gamma, tmp, 4); /* D^2 - 4AF */
    mpz_mul(inner, beta, beta); mpz_mul(tmp, delta, gamma);
    mpz_submul_ui(inner, tmp, 4);                          /* beta^2 - 4 delta gamma */

    bool handled = true;
    if (mpz_sgn(inner) < 0) {                              /* no real y -> {} proof */
        st->max_visits = SI_MAX_NODES;
    } else {
        double dDelta = mpz_get_d(delta), dBeta = mpz_get_d(beta), dInner = mpz_get_d(inner);
        double sq = sqrt(dInner < 0 ? 0 : dInner);
        double rA = (-dBeta + sq) / (2.0 * dDelta), rB = (-dBeta - sq) / (2.0 * dDelta);
        double dlo = (rA < rB ? rA : rB) - 2.0, dhi = (rA < rB ? rB : rA) + 2.0;
        if (!(dlo > -9e17 && dhi < 9e17) || (dhi - dlo) > (double)SI_MAX_NODES) {
            handled = false;                               /* too wide / overflow -> decline */
        } else {
            st->max_visits = SI_MAX_NODES;
            int64_t ylo = (int64_t)floor(dlo), yhi = (int64_t)ceil(dhi);
            mpz_t a1, a0, disc, root, num, twoA, yz, s2;
            mpz_inits(a1, a0, disc, root, num, twoA, yz, s2, NULL);
            mpz_mul_ui(twoA, A, 2);
            for (int64_t yy = ylo; yy <= yhi && !st->overflow; yy++) {
                mpz_set_si(yz, yy);
                /* a1 = B y + D,  a0 = C y^2 + E y + F. */
                mpz_mul(a1, B, yz); mpz_add(a1, a1, D);
                mpz_mul(a0, C, yz); mpz_add(a0, a0, E); mpz_mul(a0, a0, yz); mpz_add(a0, a0, F);
                mpz_mul(disc, a1, a1); mpz_mul(s2, A, a0); mpz_submul_ui(disc, s2, 4); /* a1^2 - 4A a0 */
                if (mpz_sgn(disc) < 0 || !mpz_perfect_square_p(disc)) continue;
                mpz_sqrt(root, disc);
                for (int sgn = 1; sgn >= -1; sgn -= 2) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    mpz_set(num, root); if (sgn < 0) mpz_neg(num, num);
                    mpz_sub(num, num, a1);                  /* -a1 +/- sqrt(disc) */
                    if (!mpz_divisible_p(num, twoA)) continue;
                    mpz_divexact(num, num, twoA);           /* x */
                    if (!mpz_fits_slong_p(num)) continue;
                    int64_t vals[SI_MAX_VARS];
                    for (int i = 0; i < n; i++) vals[i] = 0;
                    vals[Xv] = mpz_get_si(num); vals[Yv] = yy;
                    if (si_verify(c, vals)) emit_full(st, vals);
                    if (mpz_sgn(root) == 0) break;          /* double root: one x */
                }
            }
            mpz_clears(a1, a0, disc, root, num, twoA, yz, s2, NULL);
        }
    }

    mpz_clears(A, B, C, D, E, F, delta, tmp, beta, gamma, inner, NULL);
    return handled;
}

/* ================================================================== *
 *  Booker cube-root-mod-d engine for  x^3 + y^3 + z^3 == k.
 *
 *  Reference: A. R. Booker, "Cracking the problem with 33" (2019).
 *  For any solution, k - z^3 = x^3 + y^3 = (x+y)(x^2-xy+y^2), so the
 *  divisor d = |x+y| divides k - z^3, i.e.  z^3 ≡ k (mod d).  Instead
 *  of enumerating z and factoring the ~B^3 value k - z^3 (the classical
 *  divisor method, which caps the z-range), we enumerate the SMALL
 *  divisor d and solve the cube-root congruence, recovering z as a union
 *  of arithmetic progressions.  The pair {x,y} then follows in closed
 *  form from d and z, and every hit is verified exactly.
 *
 *  The primitive below computes ALL cube roots of k modulo a machine
 *  integer d (d <= SI_BK_DMAX), the correctness-critical core: a missing
 *  residue would silently drop solutions.  d is bounded so every value
 *  and intermediate fits in int64.
 * ================================================================== */

#define SI_BK_BRUTE_P   4096LL         /* brute-force prime moduli up to this   */
#define SI_BK_DMAX      3000000LL      /* max divisor d (coords to ~1e7, SPF-fed)*/
#define SI_BK_ROOTCAP   4096           /* max cube roots mod d before declining  */
#define SI_BK_SOLCAP    200000         /* decline if the box holds a huge family */
#define SI_BK_MAXNODES  1000000000LL   /* part-2 candidate backstop (decline)    */

/* Modular inverse of a mod m (both positive, gcd(a,m)==1), int64. */
static int64_t si_inv_i64(int64_t a, int64_t m) {
    int64_t t = 0, newt = 1, r = m, newr = a % m;
    if (newr < 0) newr += m;
    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (r != 1) return 0;            /* not invertible */
    if (t < 0) t += m;
    return t;
}

/* (a*b) mod m for 0<=a,b<m<=~2.6e5: product < 7e10 fits int64. */
static int64_t si_mulmod_i64(int64_t a, int64_t b, int64_t m) {
    return (a * b) % m;
}

/* a^e mod m by square-and-multiply, int64 (m small). */
static int64_t si_powmod_i64(int64_t a, int64_t e, int64_t m) {
    int64_t r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) { if (e & 1) r = si_mulmod_i64(r, a, m); a = si_mulmod_i64(a, a, m); e >>= 1; }
    return r;
}

/* Trial-division factorisation of a positive int64 d (d <= SI_BK_DMAX, so the
 * largest prime factor is <= sqrt bound; cheap).  Fills primes/exps, sets *np.
 * Capacity of the arrays must be >= 16 (ample for d <= 2.6e5). */
static void si_factor_i64(int64_t d, int64_t* primes, unsigned* exps, int* np) {
    *np = 0;
    for (int64_t p = 2; p * p <= d; p++) {
        if (d % p) continue;
        unsigned e = 0;
        while (d % p == 0) { d /= p; e++; }
        primes[*np] = p; exps[*np] = e; (*np)++;
    }
    if (d > 1) { primes[*np] = d; exps[*np] = 1; (*np)++; }
}

/* Factor d via a precomputed smallest-prime-factor sieve: O(number of prime
 * factors) instead of O(sqrt d), so the divisor loop can range far wider. */
static void si_factor_spf(int64_t d, const int32_t* spf, int64_t* primes,
                          unsigned* exps, int* np) {
    *np = 0;
    while (d > 1) {
        int64_t p = spf[d]; unsigned e = 0;
        while (d % p == 0) { d /= p; e++; }
        primes[*np] = p; exps[*np] = e; (*np)++;
    }
}

/* Build the smallest-prime-factor table for [0, n]; spf[i] = least prime | i
 * (spf[0]=spf[1]=0).  Caller owns the returned buffer (NULL on OOM). */
static int32_t* si_build_spf(int64_t n) {
    int32_t* spf = (int32_t*)calloc((size_t)n + 1, sizeof(int32_t));
    if (!spf) return NULL;
    for (int64_t i = 2; i <= n; i++) {
        if (spf[i]) continue;                         /* already has a factor */
        for (int64_t j = i; j <= n; j += i)
            if (!spf[j]) spf[j] = (int32_t)i;
    }
    return spf;
}

/* Exact floor(sqrt(n)) for a nonnegative 128-bit n (long double seed + adjust). */
static int64_t si_isqrt_i128(__int128 n) {
    if (n <= 0) return 0;
    int64_t r = (int64_t)sqrtl((long double)n);
    if (r < 1) r = 1;
    while ((__int128)r * r > n) r--;
    while ((__int128)(r + 1) * (r + 1) <= n) r++;
    return r;
}

/* ALL cube roots of a (mod p), p an odd prime, a in [0,p).  Writes the
 * distinct roots into out[0..*n) (out capacity >= 3).  Always succeeds. */
static void si_croots_mod_p(int64_t* out, int* n, int64_t a, int64_t p) {
    *n = 0;
    a %= p; if (a < 0) a += p;
    if (a == 0) { out[(*n)++] = 0; return; }
    if (p == 3) { out[(*n)++] = a; return; }         /* z^3 ≡ z (mod 3) */
    if (p % 3 == 2) {                                /* cubing is a bijection */
        int64_t inv3 = si_inv_i64(3, p - 1);
        out[(*n)++] = si_powmod_i64(a, inv3, p);
        return;
    }
    /* p ≡ 1 (mod 3): 0 or 3 roots. */
    int64_t ecr = (p - 1) / 3;
    if (si_powmod_i64(a, ecr, p) != 1) return;       /* not a cubic residue */
    /* p-1 = 3^s * t, 3 ∤ t. */
    int64_t t = p - 1; int s = 0;
    while (t % 3 == 0) { t /= 3; s++; }
    /* cubic non-residue nn -> 3-Sylow generator b (order 3^s), gamma (order 3). */
    int64_t nn = 2;
    while (si_powmod_i64(nn, ecr, p) == 1) nn++;
    int64_t b = si_powmod_i64(nn, t, p);
    int64_t binv = si_inv_i64(b, p);
    int64_t threepow_sm1 = 1; for (int i = 0; i < s - 1; i++) threepow_sm1 *= 3;
    int64_t gamma = si_powmod_i64(b, threepow_sm1, p);   /* primitive cube root of 1 */
    /* x1 = a^(3^{-1} mod t); then x1^3 = a * beta^w with beta = a^t (in 3-Sylow),
     * w = (3*m - 1)/t.  Correct x1 by b^{-c}, c = dlog_b(beta^w)/3. */
    int64_t m = si_inv_i64(3, t);
    int64_t x1 = si_powmod_i64(a, m, p);
    int64_t beta = si_powmod_i64(a, t, p);
    int64_t w = (3 * m - 1) / t;                     /* exact, w in {0,1,2} */
    int64_t betaw = si_powmod_i64(beta, w, p);
    /* Pohlig-Hellman discrete log E of betaw base b in the 3-Sylow. */
    int64_t E = 0, cur = betaw, pw3 = 1;             /* pw3 = 3^i */
    for (int i = 0; i < s; i++) {
        int64_t expo = 1; for (int j = 0; j < s - 1 - i; j++) expo *= 3;
        int64_t h = si_powmod_i64(cur, expo, p);
        int digit = (h == 1) ? 0 : (h == gamma ? 1 : 2);
        if (digit) {
            E += (int64_t)digit * pw3;
            cur = si_mulmod_i64(cur, si_powmod_i64(binv, (int64_t)digit * pw3, p), p);
        }
        pw3 *= 3;
    }
    if (E % 3 != 0) return;                           /* inconsistent: no root */
    int64_t x0 = si_mulmod_i64(x1, si_powmod_i64(binv, E / 3, p), p);
    /* three roots x0 * {1, gamma, gamma^2}. */
    out[(*n)++] = x0;
    out[(*n)++] = si_mulmod_i64(x0, gamma, p);
    out[(*n)++] = si_mulmod_i64(x0, si_mulmod_i64(gamma, gamma, p), p);
}

/* ALL cube roots of a (mod p^e), a in [0,pe).  A large prime modulus (e==1)
 * takes the analytic path (<=3 roots); every other prime power is brute-forced,
 * which is correct for all p=3 and p|a cases -- but those can yield MANY roots
 * (e.g. z^3 ≡ 0 mod 2^18 has 4096 roots), so the count is capped: returns the
 * number of roots written, or -1 if it would exceed `cap` (caller declines). */
static int si_croots_mod_pe(int64_t* out, int64_t a, int64_t p,
                            unsigned e, int64_t pe, int cap) {
    if (e == 1 && p > SI_BK_BRUTE_P) { int n; si_croots_mod_p(out, &n, a, p); return n; }
    int n = 0;
    a %= pe; if (a < 0) a += pe;
    for (int64_t z = 0; z < pe; z++)
        if (si_mulmod_i64(si_mulmod_i64(z, z, pe), z, pe) == a) {
            if (n >= cap) return -1;
            out[n++] = z;
        }
    return n;
}

/* ALL residues r in [0,d) with r^3 ≡ k (mod d).  Returns the count, or -1 if
 * the root set would exceed `cap` (caller must then decline: never truncate).
 * `out` must have capacity `cap`.  `spf`, if non-NULL, factors d in O(log d)
 * (a smallest-prime-factor table covering d); otherwise trial division is used. */
static int si_all_cube_roots_mod(int64_t* out, int cap, const mpz_t k, int64_t d,
                                 const int32_t* spf) {
    if (d == 1) { if (cap < 1) return -1; out[0] = 0; return 1; }
    if (cap > SI_BK_ROOTCAP) cap = SI_BK_ROOTCAP;
    int64_t kd = mpz_fdiv_ui(k, (unsigned long)d);   /* k mod d, in [0,d) */
    int64_t primes[16]; unsigned exps[16]; int np;
    if (spf) si_factor_spf(d, spf, primes, exps, &np);
    else     si_factor_i64(d, primes, exps, &np);
    /* Accumulate via CRT over prime powers. acc holds residues mod `mod_acc`. */
    static int64_t acc[SI_BK_ROOTCAP], newacc[SI_BK_ROOTCAP], rp[SI_BK_ROOTCAP];
    int nacc = 1; acc[0] = 0; int64_t mod_acc = 1;
    for (int f = 0; f < np; f++) {
        int64_t p = primes[f]; unsigned e = exps[f];
        int64_t pe = 1; for (unsigned j = 0; j < e; j++) pe *= p;
        int nr = si_croots_mod_pe(rp, kd % pe, p, e, pe, SI_BK_ROOTCAP);
        if (nr < 0) return -1;                        /* too many roots mod pe */
        if (nr == 0) return 0;                        /* no root mod pe -> none */
        /* CRT-combine acc (mod mod_acc) with rp (mod pe). */
        if ((long long)nacc * nr > cap) return -1;
        int64_t ninv = si_inv_i64(mod_acc % pe, pe);
        int nnew = 0;
        for (int i = 0; i < nacc; i++)
            for (int j = 0; j < nr; j++) {
                int64_t diff = rp[j] - (acc[i] % pe); diff %= pe; if (diff < 0) diff += pe;
                int64_t t = si_mulmod_i64(diff, ninv, pe);
                newacc[nnew++] = acc[i] + mod_acc * t;   /* < mod_acc*pe <= d */
            }
        for (int i = 0; i < nnew; i++) acc[i] = newacc[i];
        nacc = nnew; mod_acc *= pe;
    }
    if (nacc > cap) return -1;
    for (int i = 0; i < nacc; i++) out[i] = acc[i];
    return nacc;
}

/* Internal test hook: Solve`CubeRootsMod[k, d] -> sorted list of all r in
 * [0,d) with r^3 ≡ k (mod d).  Used by the test suite to cross-check the
 * primitive against brute force. */
static Expr* builtin_cube_roots_mod(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* ke = res->data.function.args[0];
    Expr* de = res->data.function.args[1];
    if (!expr_is_integer_like(ke) || de->type != EXPR_INTEGER || de->data.integer < 1)
        return NULL;
    int64_t d = de->data.integer;
    if (d > SI_BK_DMAX) return NULL;
    mpz_t k; mpz_init(k); expr_to_mpz(ke, k);
    enum { CAP = 4096 };
    int64_t roots[CAP];
    int nr = si_all_cube_roots_mod(roots, CAP, k, d, NULL);
    mpz_clear(k);
    if (nr < 0) return NULL;
    /* sort ascending (simple insertion; nr is small). */
    for (int i = 1; i < nr; i++) {
        int64_t v = roots[i]; int j = i - 1;
        while (j >= 0 && roots[j] > v) { roots[j + 1] = roots[j]; j--; }
        roots[j + 1] = v;
    }
    Expr** elems = (nr > 0) ? (Expr**)malloc((size_t)nr * sizeof(Expr*)) : NULL;
    for (int i = 0; i < nr; i++) elems[i] = expr_new_integer(roots[i]);
    Expr* list = mk_list(elems, (size_t)nr);
    free(elems);
    return list;
}

/* Recover the pair {x,y} with x+y = ±d and x^3+y^3 = m (m != 0, d | m), the
 * closed form  {x,y} = (sgn(m)*d ± sqrt((4|m|/d - d^2)/3)) / 2, and emit any
 * that lands in the box.  vals already carries the modular variable; iv/jv are
 * the pair positions. */
static void si_bk_emit_pair(SICtx* c, SearchState* st, int iv, int jv,
                            int64_t d, __int128 m, int64_t* vals) {
    if (m == 0) return;                              /* x=-y family: part1 covers it */
    __int128 am = (m < 0) ? -m : m;
    if (am % d != 0) return;                          /* not a divisor (guard) */
    __int128 q = am / d;                              /* x^2 - xy + y^2 > 0 */
    __int128 disc = 4 * q - (__int128)d * d;
    if (disc < 0 || disc % 3 != 0) return;
    disc /= 3;
    int64_t R = si_isqrt_i128(disc);
    if ((__int128)R * R != disc) return;              /* not a perfect square */
    int64_t s = (m > 0) ? d : -d;                     /* s = x + y = sgn(m)*d */
    for (int rs = 1; rs >= -1; rs -= 2) {
        int64_t num = s + (int64_t)rs * R;
        if (num & 1) continue;                        /* x must be integral */
        int64_t xv = num / 2, yv = s - xv;
        vals[iv] = xv; vals[jv] = yv;
        if (si_verify(c, vals)) emit_full(st, vals);
        if (R == 0) break;                            /* double root: one pair */
    }
}

/* Detect  x^3 + y^3 + z^3 == k  (unit coefficients, k a nonzero integer, three
 * box-bounded variables) and solve it by the Booker cube-root-mod-d method.
 *
 * For any solution the divisor d = |a+b| of a pair {a,b} divides k - c^3 of the
 * third variable c, i.e. c^3 ≡ k (mod d).  So instead of LINEARLY enumerating
 * the third variable (the classical outer loop, whose range is the bottleneck),
 * we enumerate the SMALL divisor d up to alpha*B (alpha = 2^(1/3)-1), cube-root
 * k mod d to pin c to arithmetic progressions, and recover {a,b} in closed
 * form.  This lifts the reachable coordinate range from the ~3e5 outer cap
 * toward ~1e6.  Completeness (across the three roles):
 *   - part 1: every solution whose smallest |coordinate| <= T = floor(k^{1/3})
 *     is found by enumerating that small coordinate and divisor-solving the pair
 *     (this also covers the k=c^3 => a=-b family);
 *   - part 2: every solution whose coordinates all exceed T satisfies Booker's
 *     bound d < alpha*|c| <= alpha*B for the pairing of its two largest, so it
 *     is found in the role where c is the smallest coordinate.
 * Returns true (a complete solution set emitted) only when the box is provably
 * covered by the budget; otherwise declines so the normal pipeline continues. */
static bool si_solve_three_cubes_booker(SICtx* c, SearchState* st) {
    if (c->n != 3 || c->neq != 1) return false;
    for (int i = 0; i < 3; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;
    const MPoly* eq = c->eq[0];

    /* Shape: each variable appears once as v^3 with coefficient +/-1, plus one
     * constant term -k.  K accumulates -constant = k; sgn[i] records the sign of
     * v_i^3 (a -1 is handled by the substitution u_i = -v_i below, which turns
     * the system into the pure  u_x^3 + u_y^3 + u_z^3 == k  over a mirrored box:
     * e.g.  x^3 + y^3 - z^3 == 227  is the sum of three cubes for 227). */
    mpz_t K; mpz_init(K); mpz_set_ui(K, 0);
    int seen[3] = {0, 0, 0}, sgn[3] = {1, 1, 1}; bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * 3;
        int nz = -1, cnt = 0;
        for (int v = 0; v < 3; v++) if (ex[v] > 0) { nz = v; cnt++; }
        if (cnt == 0) { mpz_sub(K, K, eq->coefs[t]); continue; }   /* constant */
        if (cnt > 1 || ex[nz] != 3 || seen[nz]) { ok = false; break; }
        if (mpz_cmp_ui(eq->coefs[t], 1) == 0) sgn[nz] = 1;
        else if (mpz_cmp_si(eq->coefs[t], -1) == 0) sgn[nz] = -1;
        else { ok = false; break; }                               /* coeff not +/-1 */
        seen[nz] = 1;
    }
    for (int i = 0; i < 3 && ok; i++) if (!seen[i]) ok = false;
    /* |k| must be modest so part 1's pair factoring stays cheap. */
    if (ok && mpz_sizeinbase(K, 2) > 30) ok = false;              /* |k| < ~1e9 */
    if (!ok || mpz_sgn(K) == 0) { mpz_clear(K); return false; }
    bool signed_case = (sgn[0] < 0 || sgn[1] < 0 || sgn[2] < 0);
    /* The u_i = sgn_i v_i substitution rewrites the equation with unit
     * coefficients and mirrors the flipped variable's box; it only stays exact
     * for a pure box (no orderings / disequations, every constraint captured),
     * since a mixed-sign ordering is not an ordering in u-space. */
    if (signed_case && (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured)) {
        mpz_clear(K); return false;
    }
    int64_t Ki = mpz_get_si(K);
    int64_t absk = Ki < 0 ? -Ki : Ki;

    /* Box radius B and divisor bound Dmax = ceil(alpha*B).  The u_i = sgn_i v_i
     * mirror preserves |bounds|, so B is the same in u-space -- compute from c. */
    int64_t B = 0;
    for (int i = 0; i < 3; i++) {
        int64_t a = c->hi[i] > -c->lo[i] ? c->hi[i] : -c->lo[i];
        if (a > B) B = a;
    }
    bool force = getenv("MATHILDA_BK_FORCE") != NULL;
    /* The cube-root-mod-d method's O(alpha*B * roots) work beats the leaf
     * search's O(B^2) and the divisor method's O(B * factoring) across the whole
     * bounded range -- measured ~10-400x faster from B~500 up, and validated
     * to return the identical solution set (Booker vs the leaf/divisor pipeline)
     * over a wide (k, box) sweep incl. the (a,-a,cbrt k) family and signed
     * equations.  So engage it for any non-trivial box; only a truly tiny box
     * (|coord| <= 100), already solved in well under a millisecond by the
     * leaf/mitm paths, is left alone.  The upper bound is Dmax <= SI_BK_DMAX
     * below (past which int64 coverage is not guaranteed). */
    if (!force && 2 * B <= 200) { mpz_clear(K); return false; }
    int64_t Dmax = (int64_t)(0.25992104989 * (double)B) + 2;
    if (Dmax > SI_BK_DMAX) { mpz_clear(K); return false; }

    /* T = floor(|k|^{1/3}). */
    int64_t T = 0; while ((T + 1) * (T + 1) * (T + 1) <= absk) T++;

    /* Working context cc: the pure  u_x^3+u_y^3+u_z^3 == k  problem.  The all-+
     * case uses c directly; a sign flips the corresponding box and swaps the
     * equation MPoly for unit coefficients (freed at the end).  Built after the
     * early declines so those paths need not free it. */
    SICtx cu; SICtx* cc = c; MPoly* newpoly = NULL;
    if (signed_case) {
        cu = *c;
        for (int i = 0; i < 3; i++)
            if (sgn[i] < 0) { cu.lo[i] = -c->hi[i]; cu.hi[i] = -c->lo[i]; }
        MPoly* p = mpoly_monomial(3, 0, 3, 1);
        for (int i = 1; i < 3; i++) {
            MPoly* tm = mpoly_monomial(3, i, 3, 1);
            MPoly* s = mpoly_add(p, tm); mpoly_free(p); mpoly_free(tm); p = s;
        }
        MPoly* kk = mpoly_from_mpz(3, K);
        newpoly = mpoly_sub(p, kk); mpoly_free(p); mpoly_free(kk);
        cu.eq[0] = newpoly; cu.neq = 1;
        cc = &cu;
    }

    /* Smallest-prime-factor sieve over [0, Dmax]: factors every divisor in
     * O(log d), so the enumeration is dominated by the candidate loop, not by
     * factoring.  Sized to Dmax, so small boxes stay cheap. */
    int32_t* spf = si_build_spf(Dmax);
    if (!spf) { if (newpoly) mpoly_free(newpoly); mpz_clear(K); return false; }

    st->max_visits = SI_BK_MAXNODES;
    int64_t vals[SI_MAX_VARS];

    /* ---- part 1: smallest coordinate <= T (per role), divisor-solve pair. ---- */
    for (int r = 0; r < 3 && !st->overflow; r++) {
        int iv = (r + 1) % 3, jv = (r + 2) % 3;
        int64_t lo = cc->lo[r] > -T ? cc->lo[r] : -T;
        int64_t hi = cc->hi[r] <  T ? cc->hi[r] :  T;
        mpz_t m, v3; mpz_init(m); mpz_init(v3);
        for (int64_t v = lo; v <= hi && !st->overflow; v++) {
            mpz_set_si(v3, v); mpz_mul_si(v3, v3, v); mpz_mul_si(v3, v3, v);
            mpz_sub(m, K, v3);                         /* m = k - v^3 */
            for (int q = 0; q < 3; q++) vals[q] = 0;
            vals[r] = v;
            si_two_power_solve(cc, st, iv, jv, 3, m, vals);
        }
        mpz_clears(m, v3, NULL);
    }

    /* A box that holds a huge parametric family (e.g. (a,-a,c) when k=c^3) is
     * better left to a symbolic/parametric treatment than enumerated; decline
     * rather than materialise hundreds of thousands of tuples. */
    bool incomplete = st->nsol > SI_BK_SOLCAP;

    /* ---- part 2: all coordinates > T, Booker divisor enumeration. ---- */
    int64_t roots[SI_BK_ROOTCAP];
    for (int64_t d = 1; d <= Dmax && !st->overflow && !incomplete; d++) {
        if (st->nsol > SI_BK_SOLCAP) { incomplete = true; break; }
        int nr = si_all_cube_roots_mod(roots, SI_BK_ROOTCAP, K, d, spf);
        if (nr < 0) { incomplete = true; break; }     /* too many roots: bail */
        if (nr == 0) continue;
        for (int r = 0; r < 3 && !st->overflow; r++) {
            int iv = (r + 1) % 3, jv = (r + 2) % 3;
            int64_t lo = cc->lo[r], hi = cc->hi[r];
            for (int ridx = 0; ridx < nr; ridx++) {
                int64_t start = lo + ((((roots[ridx] - lo) % d) + d) % d);
                for (int64_t v = start; v <= hi && !st->overflow; v += d) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    if (v >= -T && v <= T) continue;   /* handled by part 1 */
                    __int128 v3 = (__int128)v * v * v; /* |v| up to ~1e7 */
                    __int128 m = (__int128)Ki - v3;
                    for (int q = 0; q < 3; q++) vals[q] = 0;
                    vals[r] = v;
                    si_bk_emit_pair(cc, st, iv, jv, d, m, vals);
                }
            }
        }
    }

    free(spf);
    /* Candidates were enumerated in u-space; map each back to the original
     * variables (v_i = sgn_i u_i) before the result is built. */
    if (!incomplete && !st->overflow && signed_case)
        for (int i = 0; i < st->nsol; i++) {
            int64_t* row = st->sols + (size_t)i * 3;
            for (int q = 0; q < 3; q++) if (sgn[q] < 0) row[q] = -row[q];
        }
    if (newpoly) mpoly_free(newpoly);
    mpz_clear(K);
    if (incomplete || st->overflow) { st->nsol = 0; return false; }
    return true;
}

/* Dispatch the special forms.  Returns true if one handled the input
 * (candidates emitted into st). */
static bool si_try_special_forms(SICtx* c, SearchState* st) {
    if (si_solve_three_cubes_booker(c, st)) return true;
    if (si_solve_pell(c, st)) return true;
    if (si_solve_conic(c, st)) return true;
    if (si_solve_factorable_conic(c, st)) return true;
    if (si_solve_elliptic_bqf(c, st)) return true;
    if (si_solve_reciprocal(c, st)) return true;
    if (si_solve_linelim_bilinear(c, st)) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Exponential Diophantine  (variable exponents).                     *
 *                                                                     *
 *  Equations like  x^a - y^b == 1  cannot be represented as polynomials
 *  (the exponent is a variable), so they are handled here, before the  *
 *  MPoly stage.  A fully bounded box is enumerated exactly; the        *
 *  unbounded Catalan shape  x^a - y^b == +/-1  (bases and exponents    *
 *  >= 2) is settled by Mihailescu's theorem: the only solution is      *
 *  3^2 - 2^3 = 1.                                                      *
 * ------------------------------------------------------------------ */

typedef struct { int coef; int bvar; int64_t bconst; int evar; int64_t econst; } SIExpTerm;

/* Collect additive power terms of `e` (multiplied by outer `sign`) into t[],
 * folding integer constants into K.  Returns false on an unrecognised shape. */
static bool si_exp_collect(Expr* e, int sign, Expr** var, int n,
                           SIExpTerm* t, int* nt, mpz_t K) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Plus) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!si_exp_collect(e->data.function.args[i], sign, var, n, t, nt, K)) return false;
        return true;
    }
    int coef = sign; Expr* powfac = e;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER) {
        coef = sign * (int)e->data.function.args[0]->data.integer;
        powfac = e->data.function.args[1];
    }
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {   /* constant */
        mpz_t v; expr_to_mpz(e, v);
        if (sign < 0) mpz_sub(K, K, v); else mpz_add(K, K, v);
        mpz_clear(v);
        return true;
    }
    /* powfac must be base^exp, a bare variable (exp 1), or Power[base, exp]. */
    Expr *base = NULL, *exq = NULL;
    if (powfac->type == EXPR_SYMBOL) { base = powfac; exq = NULL; }
    else if (powfac->type == EXPR_FUNCTION && powfac->data.function.head->type == EXPR_SYMBOL
             && powfac->data.function.head->data.symbol.name == SYM_Power
             && powfac->data.function.arg_count == 2) {
        base = powfac->data.function.args[0];
        exq  = powfac->data.function.args[1];
    } else return false;

    if (*nt >= 8) return false;
    SIExpTerm* pt = &t[*nt];
    pt->coef = coef; pt->bvar = -1; pt->bconst = 0; pt->evar = -1; pt->econst = 1;
    if (base->type == EXPR_SYMBOL) {
        int bi = -1; for (int i = 0; i < n; i++) if (var[i]->data.symbol.name == base->data.symbol.name) bi = i;
        if (bi < 0) return false;
        pt->bvar = bi;
    } else if (base->type == EXPR_INTEGER) pt->bconst = base->data.integer;
    else return false;
    if (exq == NULL) pt->econst = 1;
    else if (exq->type == EXPR_SYMBOL) {
        int ei = -1; for (int i = 0; i < n; i++) if (var[i]->data.symbol.name == exq->data.symbol.name) ei = i;
        if (ei < 0) return false;
        pt->evar = ei;
    } else if (exq->type == EXPR_INTEGER) pt->econst = exq->data.integer;
    else return false;
    (*nt)++;
    return true;
}

/* Build {{v -> val, ...}} for the single tuple `vals`. */
static Expr* si_exp_one_tuple(Expr** var, const int64_t* vals, int n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) rules[i] = mk_rule(expr_copy(var[i]), mk_int(vals[i]));
    Expr* tup = mk_list(rules, (size_t)n); free(rules);
    return tup;
}

/* Is  target == base^exp  for some exp >= 1?  (base >= 2, target >= 1.) */
static bool si_int_ppow(int64_t base, int64_t target, int64_t* exp) {
    if (base < 2 || target < 2) return false;        /* base^exp >= 2 for exp >= 1 */
    int64_t v = base, e = 1;
    while (v < target) {
        if (v > INT64_MAX / base) return false;
        v *= base; e++;
    }
    if (v == target) { *exp = e; return true; }
    return false;
}

/* Verify a full integer assignment `vals` against the original conjunction
 * `expr` symbolically (used by paths whose equation is not an MPoly). */
static bool si_verify_symbolic(Expr* expr, Expr** var, const int64_t* vals, int n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) rules[i] = mk_rule(expr_copy(var[i]), mk_int(vals[i]));
    Expr* rl = mk_list(rules, (size_t)n); free(rules);
    Expr* out = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(expr), rl }, 2));
    bool ok = is_sym(out, SYM_True);
    expr_free(out);
    return ok;
}

/* Try to solve an exponential Diophantine equation.  Returns the owned result
 * List, or NULL to decline (not this shape, or not decidable here). */
static Expr* si_solve_exponential(Expr* expr, Expr** var, int n) {
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);

    /* Find the (first) equation whose residual carries a variable exponent. */
    Expr* eqn = NULL;
    for (int i = 0; i < ncj; i++)
        if (is_fun(conj[i], SYM_Equal, 2)) { eqn = conj[i]; break; }
    if (!eqn) return NULL;

    SIExpTerm t[8]; int nt = 0;
    mpz_t K; mpz_init_set_ui(K, 0);
    bool okp = si_exp_collect(eqn->data.function.args[0], +1, var, n, t, &nt, K)
            && si_exp_collect(eqn->data.function.args[1], -1, var, n, t, &nt, K);
    mpz_neg(K, K);                                     /* sum of power terms == K */
    bool has_var_exp = false;
    for (int i = 0; i < nt; i++) if (t[i].evar >= 0) has_var_exp = true;
    if (!okp || !has_var_exp) { mpz_clear(K); return NULL; }

    /* Parse the constraints for per-variable bounds (skip the equation). */
    SICtx c; memset(&c, 0, sizeof(c)); c.var = var; c.n = n; c.original = expr; c.all_captured = true;
    for (int i = 0; i < ncj; i++) if (conj[i] != eqn) classify_conjunct(&c, conj[i]);

    Expr* result = NULL;

    /* Catalan / Mihailescu shape: p^m - q^n == +/-1, bases and exponents all
     * variables constrained to >= 2.  Unique solution 3^2 - 2^3 = 1. */
    if (nt == 2 && (mpz_cmp_si(K, 1) == 0 || mpz_cmp_si(K, -1) == 0)) {
        int ip = (t[0].coef > 0) ? 0 : 1, in = 1 - ip;   /* + and - terms */
        bool shape = t[ip].coef == 1 && t[in].coef == -1
            && t[ip].bvar >= 0 && t[ip].evar >= 0 && t[in].bvar >= 0 && t[in].evar >= 0;
        int b9 = -1, e9 = -1, b8 = -1, e8 = -1;          /* the ^ that is 9 (3^2) / 8 (2^3) */
        if (shape) {
            if (mpz_cmp_si(K, 1) == 0) { b9 = t[ip].bvar; e9 = t[ip].evar; b8 = t[in].bvar; e8 = t[in].evar; }
            else                        { b9 = t[in].bvar; e9 = t[in].evar; b8 = t[ip].bvar; e8 = t[ip].evar; }
            /* every base/exponent must admit the required 3,2 / 2,3 and be >= 2. */
            int need_lo[4] = {2, 2, 2, 2};
            int idx[4] = {b9, e9, b8, e8}; int val[4] = {3, 2, 2, 3};
            bool feas = true;
            for (int q = 0; q < 4 && feas; q++) {
                int v = idx[q];
                if (!(c.has_lo[v] && c.lo[v] >= need_lo[q])) feas = false;   /* bases,exps >= 2 */
                if (c.has_lo[v] && val[q] < c.lo[v]) feas = false;
                if (c.has_hi[v] && val[q] > c.hi[v]) feas = false;
            }
            /* the four variables must be distinct (x,y,a,b). */
            if (b9==e9||b9==b8||b9==e8||e9==b8||e9==e8||b8==e8) feas = false;
            if (feas) {
                int64_t vals[SI_MAX_VARS];
                for (int i = 0; i < n; i++) vals[i] = 0;
                vals[b9] = 3; vals[e9] = 2; vals[b8] = 2; vals[e8] = 3;
                Expr* tup = si_exp_one_tuple(var, vals, n);
                result = mk_list((Expr*[]){ tup }, 1);
            } else {
                result = mk_list(NULL, 0);                /* provably no solution */
            }
        }
    }

    /* Fixed-base Pillai:  P^m - Q^n == +/-1  with CONSTANT bases P,Q >= 2 and
     * VARIABLE exponents m,n.  Sound & complete via Mihailescu: solutions with
     * both exponents >= 2 are consecutive perfect powers, i.e. only 3^2 - 2^3 = 1;
     * the exponent-1 cases are  Q^n = P-1  and  P^m = Q+1.  Settles 3^m - 2^n == 1
     * -> {(1,1),(2,3)} even though m,n are unbounded. */
    if (!result && n == 2 && (mpz_cmp_si(K, 1) == 0 || mpz_cmp_si(K, -1) == 0)) {
        int ip = (t[0].coef > 0) ? 0 : 1, in = 1 - ip;
        bool shape = t[ip].coef == 1 && t[in].coef == -1
            && t[ip].bvar < 0 && t[in].bvar < 0          /* constant bases */
            && t[ip].evar >= 0 && t[in].evar >= 0         /* variable exponents */
            && t[ip].bconst >= 2 && t[in].bconst >= 2;
        if (shape) {
            /* Normalise to  P^a - Q^b == 1  (av,bv = exponent variable indices). */
            int64_t P, Q; int av, bv;
            if (mpz_cmp_si(K, 1) == 0) { P = t[ip].bconst; av = t[ip].evar; Q = t[in].bconst; bv = t[in].evar; }
            else                       { P = t[in].bconst; av = t[in].evar; Q = t[ip].bconst; bv = t[ip].evar; }
            if (av != bv) {
                int64_t sa[8], sb[8]; int ns = 0; int64_t e;
                if (P == 3 && Q == 2) { sa[ns] = 2; sb[ns] = 3; ns++; }       /* both exps >= 2 */
                if (si_int_ppow(Q, P - 1, &e)) { sa[ns] = 1; sb[ns] = e; ns++; } /* a == 1 */
                if (si_int_ppow(P, Q + 1, &e)) { sa[ns] = e; sb[ns] = 1; ns++; } /* b == 1 */
                for (int u = 0; u + 1 < ns; u++)         /* sort ascending by (a, b) */
                    for (int w = u + 1; w < ns; w++)
                        if (sa[w] < sa[u] || (sa[w] == sa[u] && sb[w] < sb[u])) {
                            int64_t ta = sa[u], tb = sb[u]; sa[u] = sa[w]; sb[u] = sb[w]; sa[w] = ta; sb[w] = tb;
                        }
                Expr** sols = NULL; int nsol = 0, cap = 0;
                for (int s = 0; s < ns; s++) {
                    int64_t vals[SI_MAX_VARS];
                    for (int i = 0; i < n; i++) vals[i] = 0;
                    vals[av] = sa[s]; vals[bv] = sb[s];
                    bool dup = false;                     /* de-dup (a=b=1 arises twice) */
                    for (int s2 = 0; s2 < s && !dup; s2++) if (sa[s2] == sa[s] && sb[s2] == sb[s]) dup = true;
                    if (dup) continue;
                    if (si_verify_symbolic(expr, var, vals, n)) {
                        if (nsol == cap) { cap = cap ? cap * 2 : 4; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                        sols[nsol++] = si_exp_one_tuple(var, vals, n);
                    }
                }
                result = mk_list(sols, (size_t)nsol);
                free(sols);
            }
        }
    }

    /* Otherwise: fully bounded enumeration.  Every base and exponent variable
     * must have a finite box; enumerate and test exactly. */
    if (!result) {
        int uvar[SI_MAX_VARS], nu = 0;
        bool bounded = true;
        for (int i = 0; i < n; i++) {
            bool used = false;
            for (int k = 0; k < nt; k++) if (t[k].bvar == i || t[k].evar == i) used = true;
            if (!used) continue;
            if (!(c.has_lo[i] && c.has_hi[i])) { bounded = false; break; }
            uvar[nu++] = i;
        }
        long double space = 1.0L;
        for (int q = 0; q < nu; q++) space *= (long double)(c.hi[uvar[q]] - c.lo[uvar[q]] + 1);
        if (bounded && space <= 5.0e7L) {
            Expr** sols = NULL; int nsol = 0, cap = 0;
            int64_t val[SI_MAX_VARS];
            for (int i = 0; i < n; i++) val[i] = 0;
            for (int q = 0; q < nu; q++) val[uvar[q]] = c.lo[uvar[q]];
            mpz_t acc, term, be; mpz_init(acc); mpz_init(term); mpz_init(be);
            for (;;) {
                mpz_set_ui(acc, 0);
                for (int k = 0; k < nt; k++) {
                    int64_t bb = (t[k].bvar >= 0) ? val[t[k].bvar] : t[k].bconst;
                    int64_t ee = (t[k].evar >= 0) ? val[t[k].evar] : t[k].econst;
                    if (ee < 0) { mpz_set_ui(be, 0); }     /* negative exponent -> skip term as 0 (nonintegral) */
                    else { mpz_set_si(be, bb); mpz_pow_ui(be, be, (unsigned long)ee); }
                    mpz_mul_si(term, be, t[k].coef); mpz_add(acc, acc, term);
                }
                if (mpz_cmp(acc, K) == 0) {
                    if (nsol == cap) { cap = cap ? cap * 2 : 8; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                    sols[nsol++] = si_exp_one_tuple(var, val, n);
                }
                int q = 0;
                for (; q < nu; q++) { if (++val[uvar[q]] <= c.hi[uvar[q]]) break; val[uvar[q]] = c.lo[uvar[q]]; }
                if (q == nu) break;
            }
            mpz_clear(acc); mpz_clear(term); mpz_clear(be);
            result = mk_list(sols, (size_t)nsol);
            free(sols);
        }
    }

    for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q);
    for (int i = 0; i < c.neq; i++) mpoly_free(c.eq[i]);
    mpz_clear(K);
    return result;
}

/* ------------------------------------------------------------------ *
 *  A4: non-polynomial bounded power-leaf.                             *
 *                                                                     *
 *  An equation like  n! + 1 == m^2  cannot become an MPoly (Factorial *
 *  is not a monomial), so Stage A declines it.  When one side is a     *
 *  pure power  m^e  of a leaf variable m that appears nowhere else,    *
 *  and every OTHER variable is bounded to a small box, we enumerate    *
 *  those variables, evaluate the other side numerically through the    *
 *  interpreter (which knows Factorial, Binomial, ...), and solve m by  *
 *  an exact integer e-th root.  m itself may be unbounded (it is       *
 *  determined, not enumerated), so its value may be a bignum.          *
 * ------------------------------------------------------------------ */

static bool si_expr_mentions(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type != EXPR_FUNCTION) return false;
    if (si_expr_mentions(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (si_expr_mentions(e->data.function.args[i], name)) return true;
    return false;
}

/* Is `side` a bare  Power[m, e]  with m a solve variable (index -> *leaf),
 * e a positive integer literal (-> *e), and m absent from `other`? */
static bool si_power_leaf_side(Expr* side, Expr* other, Expr** var, int n,
                               int* leaf, int* e) {
    if (!is_fun(side, SYM_Power, 2)) return false;
    Expr* base = side->data.function.args[0];
    Expr* exq  = side->data.function.args[1];
    if (base->type != EXPR_SYMBOL || exq->type != EXPR_INTEGER || exq->data.integer < 1)
        return false;
    int mi = -1;
    for (int i = 0; i < n; i++)
        if (var[i]->data.symbol.name == base->data.symbol.name) mi = i;
    if (mi < 0 || si_expr_mentions(other, base->data.symbol.name)) return false;
    *leaf = mi; *e = (int)exq->data.integer;
    return true;
}

static Expr* si_solve_bounded_powerleaf(Expr* expr, Expr** var, int n) {
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
    Expr* eqn = NULL;
    for (int i = 0; i < ncj; i++)
        if (is_fun(conj[i], SYM_Equal, 2)) { eqn = conj[i]; break; }
    if (!eqn) return NULL;

    SICtx c; memset(&c, 0, sizeof(c)); c.var = var; c.n = n; c.original = expr; c.all_captured = true;

    /* Gate: engage only for a NON-polynomial equation (the polynomial leaf
     * search is better and already handles the polynomial power-leaf case). */
    MPoly* testQ = relation_to_mpoly(&c, eqn->data.function.args[0], eqn->data.function.args[1]);
    if (testQ) { mpoly_free(testQ); for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Identify the power-leaf side and the evaluable target side. */
    Expr *A = eqn->data.function.args[0], *B = eqn->data.function.args[1];
    int leaf = -1, e = 0; Expr* target = NULL;
    if (si_power_leaf_side(A, B, var, n, &leaf, &e)) target = B;
    else if (si_power_leaf_side(B, A, var, n, &leaf, &e)) target = A;
    if (leaf < 0) { for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Parse constraints for the OTHER variables' bounds (skip the equation). */
    for (int i = 0; i < ncj; i++) if (conj[i] != eqn) classify_conjunct(&c, conj[i]);

    /* Every non-leaf variable must be finitely bounded and enumerable. */
    int ov[SI_MAX_VARS], nov = 0; long double box = 1.0L;
    bool ok = true;
    for (int i = 0; i < n; i++) {
        if (i == leaf) continue;
        if (!(c.has_lo[i] && c.has_hi[i]) || c.hi[i] < c.lo[i]) { ok = false; break; }
        ov[nov++] = i;
        box *= (long double)(c.hi[i] - c.lo[i] + 1);
    }
    if (!ok || box > 5.0e6L) { for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q); return NULL; }

    /* Enumerate the other variables; solve m by an exact e-th root. */
    Expr** sols = NULL; int nsol = 0, cap = 0;
    int64_t val[SI_MAX_VARS];
    for (int i = 0; i < n; i++) val[i] = 0;
    for (int q = 0; q < nov; q++) val[ov[q]] = c.lo[ov[q]];
    mpz_t T, root, chk; mpz_init(root); mpz_init(chk);  /* T is init-set by expr_to_mpz */
    for (;;) {
        /* target evaluated at the current other-variable assignment. */
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)(nov ? nov : 1));
        for (int q = 0; q < nov; q++) rules[q] = mk_rule(expr_copy(var[ov[q]]), mk_int(val[ov[q]]));
        Expr* rl = mk_list(rules, (size_t)nov); free(rules);
        Expr* tv = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(target), rl }, 2));
        bool got = expr_is_integer_like(tv);
        if (got) expr_to_mpz(tv, T);    /* initialises T */
        expr_free(tv);

        if (got) {
            /* m^e == T.  Even e needs T >= 0 and yields +/- root; odd e keeps sign. */
            bool have_m = false; int msign = 1;
            if (e % 2 == 0) {
                if (mpz_sgn(T) >= 0 && mpz_root(root, T, (unsigned long)e)) have_m = true;
            } else {
                mpz_abs(chk, T);
                if (mpz_root(root, chk, (unsigned long)e)) { have_m = true; msign = mpz_sgn(T) < 0 ? -1 : 1; }
            }
            if (have_m) {
                for (int sg = 1; sg >= -1; sg -= 2) {
                    if (e % 2 == 1 && sg != msign) continue;      /* odd: single sign */
                    if (e % 2 == 0 && sg < 0 && mpz_sgn(root) == 0) continue; /* avoid -0 dup */
                    mpz_t mval; mpz_init(mval); mpz_set(mval, root); if (sg < 0) mpz_neg(mval, mval);
                    /* Verify the full conjunction symbolically (handles Factorial etc.). */
                    Expr** vr = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                    for (int i = 0; i < n; i++)
                        vr[i] = (i == leaf) ? mk_rule(expr_copy(var[i]), mk_mpz(mval))
                                            : mk_rule(expr_copy(var[i]), mk_int(val[i]));
                    Expr* vrl = mk_list(vr, (size_t)n); free(vr);
                    Expr* chk2 = eval_and_free(internal_replace_all((Expr*[]){ expr_copy(expr), vrl }, 2));
                    bool good = is_sym(chk2, SYM_True);
                    expr_free(chk2);
                    if (good) {
                        Expr** tr = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
                        for (int i = 0; i < n; i++)
                            tr[i] = (i == leaf) ? mk_rule(expr_copy(var[i]), mk_mpz(mval))
                                                : mk_rule(expr_copy(var[i]), mk_int(val[i]));
                        if (nsol == cap) { cap = cap ? cap * 2 : 8; sols = realloc(sols, sizeof(Expr*) * (size_t)cap); }
                        sols[nsol++] = mk_list(tr, (size_t)n); free(tr);
                    }
                    mpz_clear(mval);
                }
            }
        }
        if (got) mpz_clear(T);          /* T was init-set this iteration */
        int q = 0;
        for (; q < nov; q++) { if (++val[ov[q]] <= c.hi[ov[q]]) break; val[ov[q]] = c.lo[ov[q]]; }
        if (q == nov) break;
    }
    mpz_clear(root); mpz_clear(chk);
    for (int i = 0; i < c.nbc; i++) mpoly_free(c.bc[i].Q);
    for (int i = 0; i < c.neq; i++) mpoly_free(c.eq[i]);

    Expr* result = mk_list(sols, (size_t)nsol);
    free(sols);
    return result;
}

/* Fermat's Last Theorem short-circuit.  a*x^n + a*y^n - a*z^n == 0 (equal
 * coefficient magnitudes, exponent n >= 3, no constant term) with x, y, z all
 * strictly positive has NO solutions (Wiles 1995) -- return {} immediately
 * rather than search the box, and independently of any upper bound so the
 * unbounded "indefinite" form is decided too.  Returns {} on a match, else NULL
 * (fall through to the ordinary machinery). */
static Expr* si_solve_fermat(SICtx* c) {
    if (c->neq != 1 || c->n != 3) return NULL;
    const MPoly* eq = c->eq[0];
    int vexp[3] = {0, 0, 0}, vsgn[3] = {0, 0, 0}, seen[3] = {0, 0, 0};
    mpz_t mag; mpz_init(mag); bool have_mag = false, ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * 3;
        int nz = -1, cnt = 0;
        for (int v = 0; v < 3; v++) if (ex[v] > 0) { nz = v; cnt++; }
        if (cnt == 0) { ok = false; break; }          /* constant term: k != 0 */
        if (cnt > 1 || seen[nz]) { ok = false; break; }
        int s = mpz_sgn(eq->coefs[t]);
        if (s == 0) { ok = false; break; }
        mpz_t a; mpz_init(a); mpz_abs(a, eq->coefs[t]);
        if (!have_mag) { mpz_set(mag, a); have_mag = true; }
        else if (mpz_cmp(mag, a) != 0) ok = false;
        mpz_clear(a);
        seen[nz] = 1; vexp[nz] = ex[nz]; vsgn[nz] = s;
    }
    mpz_clear(mag);
    if (!ok || !seen[0] || !seen[1] || !seen[2]) return NULL;
    int n = vexp[0];
    if (n < 3 || vexp[1] != n || vexp[2] != n) return NULL;   /* FLT is n >= 3 */
    int ssum = vsgn[0] + vsgn[1] + vsgn[2];
    if (ssum != 1 && ssum != -1) return NULL;         /* must be a^n + b^n == c^n */
    for (int i = 0; i < 3; i++)
        if (!(c->has_lo[i] && c->lo[i] >= 1)) return NULL;    /* strictly positive */
    return mk_list(NULL, 0);                          /* provably no solutions */
}

/* ------------------------------------------------------------------ *
 *  Thue equations  F(x, y) == m  (F irreducible homogeneous, deg >=3). *
 *                                                                     *
 *  Detect a single homogeneous binary form of degree n >= 3 equal to a *
 *  constant, extract the form coefficients a_j (of x^(n-j) y^j) and m, *
 *  and hand them to the Tzanakis-de Weger engine (src/solvethue.c).    *
 *  Returns the complete finite solution set (sorted), the empty list   *
 *  (a proof), or NULL to decline -- the engine itself declines unless  *
 *  it can certify completeness.                                        */
static Expr* si_solve_thue(SICtx* c) {
    if (c->neq != 1 || c->n != 2) return NULL;
    /* pure equation only: extra bound / ordering / disequation constraints
     * would need a bignum-aware post-filter (a later extension). */
    if (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0) return NULL;
    for (int i = 0; i < c->n; i++) if (c->has_lo[i] || c->has_hi[i]) return NULL;

    const MPoly* eq = c->eq[0];
    /* active variables */
    int act[SI_MAX_VARS], na = 0;
    for (int v = 0; v < c->n; v++) if (mpoly_deg_var(eq, v) >= 1) act[na++] = v;
    if (na != 2) return NULL;
    int xv = act[0], yv = act[1];

    /* determine the form degree (max total degree of a non-constant term) */
    int deg = 0;
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)eq->n_vars;
        for (int v = 0; v < c->n; v++) if (v != xv && v != yv && ex[v] != 0) return NULL;
        int d = ex[xv] + ex[yv];
        if (d > deg) deg = d;
    }
    if (deg < 3) return NULL;

    /* form[j] = coeff of x^(n-j) y^j (j = 0..deg); constant term = -m */
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(deg + 1));
    for (int j = 0; j <= deg; j++) mpz_init(form[j]);
    mpz_t m; mpz_init(m);
    bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * (size_t)eq->n_vars;
        int dx = ex[xv], dy = ex[yv], d = dx + dy;
        if (d == 0) { mpz_neg(m, eq->coefs[t]); }         /* residual is F - m */
        else if (d == deg) { mpz_set(form[dy], eq->coefs[t]); }
        else ok = false;                                   /* inhomogeneous -> not Thue */
    }
    if (ok && (mpz_sgn(form[0]) == 0 || mpz_sgn(form[deg]) == 0)) ok = false;
    if (ok && mpz_sgn(m) == 0) ok = false;                 /* F == 0 is not a Thue eq */

    Expr* result = NULL;
    if (ok) {
        ThueSol* sols = NULL;
        int cnt = thue_solve_binary_form((const mpz_t*)form, deg, m, &sols);
        if (cnt >= 0) {
            /* sort by (x, y) ascending for a canonical result */
            for (int i = 1; i < cnt; i++) {
                ThueSol key = sols[i]; int j = i - 1;
                while (j >= 0 && (mpz_cmp(sols[j].x, key.x) > 0 ||
                       (mpz_cmp(sols[j].x, key.x) == 0 && mpz_cmp(sols[j].y, key.y) > 0)))
                    { sols[j + 1] = sols[j]; j--; }
                sols[j + 1] = key;
            }
            Expr** rows = (cnt > 0) ? malloc(sizeof(Expr*) * (size_t)cnt) : NULL;
            for (int i = 0; i < cnt; i++) {
                Expr* rx = mk_rule(expr_copy(c->var[xv]), mk_mpz(sols[i].x));
                Expr* ry = mk_rule(expr_copy(c->var[yv]), mk_mpz(sols[i].y));
                rows[i] = mk_list((Expr*[]){ rx, ry }, 2);
            }
            result = mk_list(rows, (size_t)cnt);
            free(rows);
            thue_sols_free(sols, cnt);
        }
    }

    for (int j = 0; j <= deg; j++) mpz_clear(form[j]);
    free(form); mpz_clear(m);
    return result;
}

Expr* solveint_solve_integer(Expr* expr, Expr* vars, Expr* dom) {
    if (!expr || !vars || !dom) return NULL;
    if (!(dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Integers)) return NULL;

    /* Parse the variable list into a symbol array. */
    SICtx c;
    memset(&c, 0, sizeof(c));
    Expr* var_storage[SI_MAX_VARS];
    if (vars->type == EXPR_SYMBOL) {
        var_storage[0] = vars; c.n = 1;
    } else if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List) {
        if (vars->data.function.arg_count < 1
            || vars->data.function.arg_count > SI_MAX_VARS) return NULL;
        c.n = (int)vars->data.function.arg_count;
        for (int i = 0; i < c.n; i++) {
            if (vars->data.function.args[i]->type != EXPR_SYMBOL) return NULL;
            var_storage[i] = vars->data.function.args[i];
        }
    } else return NULL;
    c.var = var_storage;
    c.original = expr;
    c.all_captured = true;      /* cleared by any constraint the store can't hold */

    /* This pre-pass engages when there is at least one inequality / ordering /
     * disequation constraint, OR the equation is multivariable (so the
     * parametric linear path can produce a solution family).  A bare
     * single-variable equation with no constraints is left to the ordinary
     * polynomial dispatch (which, with the Integers reality + integer filter,
     * already handles x^2 == 4). */
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
    si_warn_free_symbols(&c, conj, ncj);   /* Solve::svars for stray symbols */
    bool has_constraint = false, has_equation = false;
    for (int i = 0; i < ncj; i++) {
        Expr* e = conj[i];
        if (is_fun(e, SYM_Equal, 2)) has_equation = true;
        else if (is_fun(e, SYM_Less, 2) || is_fun(e, SYM_LessEqual, 2)
              || is_fun(e, SYM_Greater, 2) || is_fun(e, SYM_GreaterEqual, 2)
              || is_fun(e, SYM_Unequal, 2)
              || (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
                  && e->data.function.head->data.symbol.name == SYM_Inequality))
            has_constraint = true;
    }
    if (!has_equation || (!has_constraint && c.n < 2)) return NULL;

    /* Exponential Diophantine (variable exponents) is handled before the MPoly
     * stage, which cannot represent x^a. */
    {
        Expr* ex = si_solve_exponential(expr, c.var, c.n);
        if (ex) return ex;
    }

    /* Non-polynomial bounded power-leaf (n! + 1 == m^2, ...): also before the
     * MPoly stage, which cannot represent Factorial / Binomial / ... */
    {
        Expr* pl = si_solve_bounded_powerleaf(expr, c.var, c.n);
        if (pl) return pl;
    }

    /* Stage A. */
    for (int i = 0; i < ncj; i++) {
        if (!classify_conjunct(&c, conj[i])) { ctx_free(&c); return NULL; }
    }
    if (c.neq == 0) { ctx_free(&c); return NULL; }

    /* Stage B. */
    derive_bounds(&c);
    derive_even_only_bounds(&c);   /* sign-symmetric even-power vars -> [-B, B] */

    /* Fermat's Last Theorem: x^n + y^n == z^n (n >= 3, all positive) has no
     * solutions -- decide it instantly, before any (possibly unbounded) search. */
    { Expr* flt = si_solve_fermat(&c); if (flt) { ctx_free(&c); return flt; } }

    /* Per-variable degree (max over equations) and whether it is solvable as
     * an exact leaf. */
    int maxdeg[SI_MAX_VARS];
    for (int i = 0; i < c.n; i++) {
        maxdeg[i] = 0;
        for (int q = 0; q < c.neq; q++) {
            int dg = mpoly_deg_var(c.eq[q], i);
            if (dg > maxdeg[i]) maxdeg[i] = dg;
        }
    }

    /* An unbounded variable is admissible ONLY as the leaf (it is solved
     * exactly, never enumerated).  Count them: two or more and the box is not
     * finite -> decline (later phases handle parametric / Pell). */
    int n_unbounded = 0, unbounded_var = -1;
    for (int i = 0; i < c.n; i++)
        if (!(c.has_lo[i] && c.has_hi[i])) { n_unbounded++; unbounded_var = i; }

    SearchState st; memset(&st, 0, sizeof(st));
    st.ctx = &c;

    /* Unconstrained single linear equation -> parametric family (symbolic,
     * so it bypasses the numeric candidate machinery). */
    {
        Expr* lin = si_solve_linear_parametric(&c);
        if (lin) { ctx_free(&c); return lin; }
    }

    /* Unbounded positive Pell -> parametric fundamental-unit family (symbolic). */
    {
        Expr* pell = si_solve_pell_parametric(&c);
        if (pell) { ctx_free(&c); return pell; }
    }

    /* Unbounded positive generalised Pell  x^2 - D y^2 == N  (|N| >= 2) -> one
     * parametric family per solution class (Nagell fundamentals + orbit). */
    {
        Expr* gp = si_solve_genpell_parametric(&c);
        if (gp) { ctx_free(&c); return gp; }
    }

    /* Homogeneous linear system with positivity -> parametric ray (symbolic). */
    {
        Expr* ray = si_solve_linear_system_ray(&c);
        if (ray) { ctx_free(&c); return ray; }
    }

    /* General unconstrained linear system (m >= 2 equations, unbounded) -> the
     * complete integer solution via HNF: particular solution + kernel lattice
     * as C[k], or {} when provably unsolvable.  (Replaces the silent wrong `{}`
     * the Complexes-oriented linear-system dispatch produced for underdetermined
     * integer systems.) */
    {
        Expr* hs = si_solve_linear_system_hnf(&c);
        if (hs) { ctx_free(&c); return hs; }
    }

    /* Prouhet-Tarry-Escott equal-power-sum system with a forcing disequation
     * -> provably {} via Newton's identities (even when unbounded). */
    {
        Expr* pte = si_solve_power_sum_equal(&c);
        if (pte) { ctx_free(&c); return pte; }
    }

    /* Unbounded Mordell y^2 == x^3 + k (k = -1, -2): the complete integer-point
     * set via factorisation in Z[sqrt k]. */
    {
        Expr* mor = si_solve_mordell(&c);
        if (mor) { ctx_free(&c); return mor; }
    }

    /* Thue equation F(x,y) == m (irreducible homogeneous binary form,
     * deg >= 3): the complete finite solution set via Tzanakis-de Weger,
     * or a decline (the engine never returns an unproven set). */
    {
        Expr* th = si_solve_thue(&c);
        if (th) { ctx_free(&c); return th; }
    }

    /* Special forms first: divisor-factoring bilinear and unit-fraction
     * recursion are exact and O(#divisors) / O(bounded), so they beat the
     * enumerative fallback whenever they match -- including fully bounded
     * systems whose box would otherwise force a large leaf search (e.g. the
     * Pythagorean-perimeter case, bounded by its linear equation). */
    if (si_try_special_forms(&c, &st)) {
        if (st.overflow) { free(st.sols); ctx_free(&c); return NULL; }
        Expr* result = build_result(&st);
        free(st.sols); ctx_free(&c);
        return result;
    }
    /* Multi-leaf staged elimination: a system whose "determined" variables each
     * fall out of a single equation (e.g. the Euler brick's a,b,c) -- these may
     * be individually unbounded, so this runs BEFORE the unbounded-box decline.
     * Only the coupled free variables are enumerated. */
    if (si_solve_multileaf(&c, &st)) {
        Expr* result = build_result(&st);
        free(st.sols); ctx_free(&c);
        return result;
    }
    st.multileaf = false;   /* reset in case the attempt set up partial state */
    st.nsol = 0;

    /* Unbounded box and no special form fit -> later phases (Pell, lattice). */
    if (n_unbounded >= 2) { free(st.sols); ctx_free(&c); return NULL; }

    /* Meet-in-the-middle fast path: a single separable additive equation is
     * solved in ~N^ceil(n/2) instead of the ~N^(n-1) leaf search. */
    if (n_unbounded == 0 && mitm_solve(&st)) {
        Expr* result = build_result(&st);
        free(st.sols);
        ctx_free(&c);
        return result;
    }

    int leaf;
    if (n_unbounded == 1) {
        /* Must be the leaf; it has to appear in an equation at a solvable
         * degree, else we cannot pin it. */
        leaf = unbounded_var;
        if (maxdeg[leaf] <= 0 || maxdeg[leaf] > SI_LEAF_MAXDEG) { ctx_free(&c); return NULL; }
        /* Give the leaf a wide finite window for exact-root filtering; the
         * final verification enforces the true (possibly one-sided) bounds. */
        const int64_t SI_WIDE = 1LL << 50;
        if (!c.has_lo[leaf]) { c.lo[leaf] = -SI_WIDE; c.has_lo[leaf] = true; }
        if (!c.has_hi[leaf]) { c.hi[leaf] =  SI_WIDE; c.has_hi[leaf] = true; }
    } else {
        /* Fully bounded: leaf = widest-domain variable that appears in some
         * equation at a solvable degree (so the widest range is never
         * enumerated); ties broken toward the lower leaf degree. */
        leaf = -1; int64_t best = -1; int best_deg = 1 << 30;
        for (int i = 0; i < c.n; i++) {
            if (maxdeg[i] <= 0 || maxdeg[i] > SI_LEAF_MAXDEG) continue;
            int64_t w = c.hi[i] - c.lo[i];
            if (w > best || (w == best && maxdeg[i] < best_deg)) {
                best = w; best_deg = maxdeg[i]; leaf = i;
            }
        }
        if (leaf < 0) leaf = 0;                       /* degenerate: no eqn var */
    }

    for (int i = 0; i < c.n; i++)
        if (c.lo[i] > c.hi[i]) { ctx_free(&c); return mk_list(NULL, 0); }  /* empty box */

    st.leaf = leaf;
    int64_t domain[SI_MAX_VARS];
    for (int i = 0; i < c.n; i++) domain[i] = c.hi[i] - c.lo[i] + 1;
    st.n_search = 0;
    for (int i = 0; i < c.n; i++) if (i != leaf) st.order[st.n_search++] = i;
    for (int i = 0; i + 1 < st.n_search; i++) {          /* sort ascending domain */
        int pick = i;
        for (int j = i + 1; j < st.n_search; j++)
            if (domain[st.order[j]] < domain[st.order[pick]]) pick = j;
        if (pick != i) { int t = st.order[i]; st.order[i] = st.order[pick]; st.order[pick] = t; }
    }
    /* Refine so every abs-ordering |a| < |b| enumerates the larger b before the
     * smaller a -- only then does effective_bounds' abs pruning fire and the box
     * collapse to ~N^L/L! nodes (matching si_longest_chain).  Stable topological
     * pass over the abs-ord DAG, seeded by the ascending-domain order above so
     * ties and unconstrained variables keep their heuristic placement. */
    if (c.n_abs_ord > 0 && st.n_search > 1) {
        int seed[SI_MAX_VARS];
        for (int i = 0; i < st.n_search; i++) seed[i] = st.order[i];
        bool placed[SI_MAX_VARS];
        for (int i = 0; i < c.n; i++) placed[i] = false;
        int out = 0;
        for (int pass = 0; pass < st.n_search && out < st.n_search; pass++) {
            for (int i = 0; i < st.n_search; i++) {
                int v = seed[i];
                if (placed[v]) continue;
                bool ready = true;                       /* all larger-|.| vars placed? */
                for (int e = 0; e < c.n_abs_ord; e++)
                    if (c.abs_ord_a[e] == v) {
                        int b = c.abs_ord_b[e];
                        if (b != leaf && !placed[b]) { ready = false; break; }
                    }
                if (ready) { st.order[out++] = v; placed[v] = true; }
            }
        }
        for (int i = 0; i < st.n_search && out < st.n_search; i++)   /* cycle leftover */
            if (!placed[seed[i]]) { st.order[out++] = seed[i]; placed[seed[i]] = true; }
    }

    /* Search-space guard: decline rather than enumerate an intractable box
     * (e.g. a Pell family or a wide linear lattice -- those are later,
     * closed-form phases).  `raw_est` is the box product over the search vars;
     * `est` divides it by the factorial of the longest ordering chain, since
     * x <= y <= z  walks ~N^3/6 nodes, not N^3 (the visit cap backstops any
     * under-estimate). */
    long double raw_est = 1.0L;
    for (int i = 0; i < st.n_search; i++) raw_est *= (long double)domain[st.order[i]];
    long double est = raw_est;
    {
        int L = si_longest_chain(&c, st.order, st.n_search);
        long double fact = 1.0L;
        for (int i = 2; i <= L; i++) fact *= (long double)i;
        est /= fact;
    }
    /* Two number-theoretic methods turn the O(N^k) search into far less: a
     * bounded linear equation via its LLL-reduced lattice, and a separable
     * odd-power sum via the divisor method (s = x+y divides m).  Gate them on
     * the RAW box size, not the ordering-pruned `est`: an ordering's L! division
     * must not drop a large separable box just under the leaf budget and into
     * the O(est) leaf enumeration -- that is exactly what made the ordered
     * three-cubes query hang.  They decline (return false) when inapplicable,
     * falling through to the leaf search only if the pruned box fits. */
    if (raw_est > (long double)SI_MAX_NODES) {
        bool done = si_solve_linear_bounded(&c, &st)
                 || si_solve_powersum_divisor(&c, &st);
        if (done && !st.overflow) {
            Expr* result = build_result(&st);
            free(st.sols); ctx_free(&c);
            return result;
        }
        if (st.overflow) { free(st.sols); ctx_free(&c); return NULL; }
        st.nsol = 0; st.multileaf = false;          /* reset any partial state */
        if (est > (long double)SI_MAX_NODES) {       /* pruned box still too big */
            free(st.sols); ctx_free(&c); return NULL;
        }
        /* else: the ordering-pruned box fits the leaf budget -- fall through. */
    }

    /* Stage C. */
    st.max_visits = SI_MAX_NODES;      /* runtime backstop for non-ordered boxes */
    si_build_leaf_cache(&st);          /* int64 coefficient cache for the hot loop */
    search_rec(&st, 0);
    si_free_leaf_cache(&st);

    if (st.overflow) {                 /* hit the solver's degree/size limit */
        free(st.sols);
        ctx_free(&c);
        return NULL;                   /* Solve stays unevaluated -- never wrong */
    }

    Expr* result = build_result(&st);
    free(st.sols);
    ctx_free(&c);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry: Solve`SolveIntegers[eqns, vars]          *
 * ------------------------------------------------------------------ */

static Expr* builtin_solve_integers(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* eqns = res->data.function.args[0];
    Expr* vars = res->data.function.args[1];
    Expr* integers = mk_sym(SYM_Integers);
    Expr* out = solveint_solve_integer(eqns, vars, integers);
    expr_free(integers);
    return out;
}

void solveint_init(void) {
    symtab_add_builtin("Solve`SolveIntegers", builtin_solve_integers);
    symtab_add_builtin("Solve`CubeRootsMod", builtin_cube_roots_mod);
    symtab_set_docstring("Solve`CubeRootsMod",
        "Solve`CubeRootsMod[k, d]\n"
        "\tInternal: the sorted list of all r in [0, d) with r^3 == k (mod d),\n"
        "\tvia factorisation + per-prime-power roots + CRT.  Backs the Booker\n"
        "\tcube-root-mod-d path for x^3 + y^3 + z^3 == k; exposed for testing.");
    symtab_set_docstring("Solve`SolveIntegers",
        "Solve`SolveIntegers[eqns, vars]\n"
        "\tInternal: solves a system of polynomial equations with\n"
        "\tinequality / ordering constraints over the integers by bound\n"
        "\tpropagation and exhaustive elimination.  Returns\n"
        "\t{{v -> n, ...}, ...} ascending, {} when there are provably no\n"
        "\tsolutions, or is left unevaluated when the variables cannot be\n"
        "\tbounded to a finite box.");
}
