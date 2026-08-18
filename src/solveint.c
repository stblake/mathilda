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
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (mirror the sibling solve specialists).  *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }
static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
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
            if (mpz_sgn(urest) <= 0) continue;

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

/* Run the bound fixpoint (explicit bounds + ordering propagation +
 * interval-positivity).  Bounds are only ever tightened, so the result is a
 * set of necessary conditions. */
static void derive_bounds(SICtx* c) {
    for (int iter = 0; iter < 4 * c->n + 8; iter++) {
        bool changed = false;
        if (propagate_orderings(c)) changed = true;
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
static int eval_leaf_coeffs(SearchState* st, const MPoly* eq, mpz_t* a, mpz_t term) {
    SICtx* c = st->ctx;
    int leaf = st->leaf, n = c->n;
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

/* Solve the leaf variable given a full search-var assignment in st->val. */
static void solve_leaf(SearchState* st) {
    SICtx* c = st->ctx;
    int leaf = st->leaf;
    int64_t lo = c->lo[leaf], hi = c->hi[leaf];
    if (++st->visits > st->max_visits) { st->overflow = true; return; }

    int64_t inter[512]; int ninter = 0; bool have = false;
    mpz_t a[SI_LEAF_MAXDEG + 1], term;
    for (int k = 0; k <= SI_LEAF_MAXDEG; k++) mpz_init(a[k]);
    mpz_init(term);

    bool feasible = true;
    for (int q = 0; q < c->neq && feasible; q++) {
        int d = eval_leaf_coeffs(st, c->eq[q], a, term);
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
}

static void search_rec(SearchState* st, int depth) {
    if (st->overflow) return;
    if (depth == st->n_search) { solve_leaf(st); return; }

    int vi = st->order[depth];
    int64_t elo, ehi;
    effective_bounds(st, depth, vi, &elo, &ehi);

    for (int64_t r = elo; r <= ehi && !st->overflow; r++) {
        st->val[vi] = r;
        search_rec(st, depth + 1);
    }
}

/* ------------------------------------------------------------------ *
 *  Result assembly.                                                   *
 * ------------------------------------------------------------------ */

static int sol_cmp(const int64_t* a, const int64_t* b, int n) {
    for (int i = 0; i < n; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

static Expr* build_result(SearchState* st) {
    SICtx* c = st->ctx;
    int ns = st->nsol, n = c->n;
    /* selection sort the (small) tuple set */
    int* order = (int*)malloc(sizeof(int) * (size_t)(ns > 0 ? ns : 1));
    for (int i = 0; i < ns; i++) order[i] = i;
    for (int i = 0; i + 1 < ns; i++) {
        int pick = i;
        for (int j = i + 1; j < ns; j++)
            if (sol_cmp(st->sols + (size_t)order[j] * n, st->sols + (size_t)order[pick] * n, n) < 0)
                pick = j;
        if (pick != i) { int t = order[i]; order[i] = order[pick]; order[pick] = t; }
    }
    Expr** tuples = (Expr**)malloc(sizeof(Expr*) * (size_t)(ns > 0 ? ns : 1));
    int nt = 0;
    for (int i = 0; i < ns; i++) {
        int64_t* row = st->sols + (size_t)order[i] * n;
        if (i > 0 && sol_cmp(st->sols + (size_t)order[i - 1] * n, row, n) == 0) continue; /* dedupe */
        Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int k = 0; k < n; k++)
            rules[k] = mk_rule(expr_copy(c->var[k]), mk_int(row[k]));
        tuples[nt++] = mk_list(rules, (size_t)n);
        free(rules);
    }
    Expr* res = mk_list(tuples, (size_t)nt);
    free(tuples);
    free(order);
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

/* Dispatch the special forms.  Returns true if one handled the input
 * (candidates emitted into st). */
static bool si_try_special_forms(SICtx* c, SearchState* st) {
    if (si_solve_pell(c, st)) return true;
    if (si_solve_reciprocal(c, st)) return true;
    if (si_solve_linelim_bilinear(c, st)) return true;
    return false;
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

    /* This pre-pass only engages when there is at least one inequality /
     * ordering / disequation constraint; a bare polynomial equation with no
     * constraints is left to the ordinary polynomial dispatch (which, with
     * the Integers reality + integer filter, already handles x^2 == 4). */
    Expr** conj; int ncj;
    flatten_conjuncts(expr, &conj, &ncj);
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
    if (!has_equation || !has_constraint) return NULL;

    /* Stage A. */
    for (int i = 0; i < ncj; i++) {
        if (!classify_conjunct(&c, conj[i])) { ctx_free(&c); return NULL; }
    }
    if (c.neq == 0) { ctx_free(&c); return NULL; }

    /* Stage B. */
    derive_bounds(&c);

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

    /* Search-space guard (raw product over search vars): decline rather than
     * enumerate an intractable box (e.g. a Pell family or a wide linear
     * lattice -- those are later, closed-form phases). */
    long double est = 1.0L;
    for (int i = 0; i < st.n_search; i++) est *= (long double)domain[st.order[i]];
    if (est > (long double)SI_MAX_NODES) { ctx_free(&c); return NULL; }

    /* Stage C. */
    st.max_visits = SI_MAX_NODES;      /* runtime backstop for non-ordered boxes */
    search_rec(&st, 0);

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
    symtab_set_docstring("Solve`SolveIntegers",
        "Solve`SolveIntegers[eqns, vars]\n"
        "\tInternal: solves a system of polynomial equations with\n"
        "\tinequality / ordering constraints over the integers by bound\n"
        "\tpropagation and exhaustive elimination.  Returns\n"
        "\t{{v -> n, ...}, ...} ascending, {} when there are provably no\n"
        "\tsolutions, or is left unevaluated when the variables cannot be\n"
        "\tbounded to a finite box.");
}
