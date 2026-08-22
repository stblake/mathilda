/*
 * solve_common.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
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
#include "solveint_internal.h"


/* Evaluate an MPoly at an integer assignment `vals` into `out` (pre-init'd).
 * Pure GMP arithmetic -- the per-candidate hot path uses this instead of a
 * symbolic re-evaluation. */
void si_eval_mpoly(const MPoly* p, const int64_t* vals, mpz_t out) {
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
MPoly* relation_to_mpoly(const SICtx* c, Expr* a, Expr* b) {
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
bool classify_conjunct(SICtx* c, Expr* e) {
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
void derive_bounds(SICtx* c) {
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
void derive_even_only_bounds(SICtx* c) {
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
int univariate_roots(const mpz_t* a, int d, int64_t lo, int64_t hi,
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


/* Exact floor(sqrt(n)) for n >= 0, division-based so it never overflows. */
int64_t si_isqrt_i64(int64_t n) {
    if (n < 0) return -1;
    if (n == 0) return 0;
    int64_t r = (int64_t)sqrt((double)n);
    if (r < 1) r = 1;
    while (r > n / r) r--;                 /* r*r > n  (safe: r > 0) */
    while ((r + 1) <= n / (r + 1)) r++;    /* (r+1)^2 <= n */
    return r;
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

Expr* build_result(SearchState* st) {
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
bool si_verify(SICtx* c, const int64_t* vals) {
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

void emit_full(SearchState* st, const int64_t* vals) {
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
 *  Top-level.                                                         *
 * ------------------------------------------------------------------ */
void ctx_free(SICtx* c) {
    for (int i = 0; i < c->neq; i++) mpoly_free(c->eq[i]);
    for (int i = 0; i < c->nbc; i++) mpoly_free(c->bc[i].Q);
}


/* Flatten a top-level And / List into a conjunct array (borrowed).  A bare
 * relation is a single conjunct. */
void flatten_conjuncts(Expr* e, Expr*** out, int* n) {
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
void si_warn_free_symbols(const SICtx* c, Expr** conj, int ncj) {
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
MPoly* si_resid_to_mpoly(Expr* resid, Expr** vars, int n) {
    Expr* tog = eval_and_free(internal_together((Expr*[]){ expr_copy(resid) }, 1));
    Expr* num = eval_and_free(internal_numerator((Expr*[]){ tog }, 1));
    Expr* exp = eval_and_free(internal_expand((Expr*[]){ num }, 1));
    MPoly* P = expr_to_mpoly(exp, vars, n);
    expr_free(exp);
    return P;
}


/* Build a total order (smallest-first) over the active variables from the
 * ordering constraints; returns false if they do not form a chain. */
bool si_build_total_order(const SICtx* c, const int* active, int ka, int* order) {
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


/* Exact floor(sqrt(n)) for a nonnegative 128-bit n (long double seed + adjust). */
int64_t si_isqrt_i128(__int128 n) {
    if (n <= 0) return 0;
    int64_t r = (int64_t)sqrtl((long double)n);
    if (r < 1) r = 1;
    while ((__int128)r * r > n) r--;
    while ((__int128)(r + 1) * (r + 1) <= n) r++;
    return r;
}
