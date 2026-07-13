/*
 * solverad.c
 *
 * Radical-equation specialist.  Backs `Solve` for single-variable
 * equations containing fractional-power subterms (Sqrt, x^(p/q), nested
 * radicals).  Algorithm summary: solverad.h.  Memory contract: every
 * helper returns a freshly-owned Expr* (or Expr**); inputs are borrowed
 * and deep-copied wherever they appear in the output.  On NULL return
 * from the top-level dispatcher, the caller leaves its enclosing Solve
 * unevaluated.
 */

#include "solverad.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly.h"
#include "solvepoly.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand.  The mk_* helpers take ownership of *
 *  their pointer arguments; they do not deep-copy.  Callers pass      *
 *  freshly built or expr_copy()'d Expr*'s.                            *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}

static Expr* mk_pow(Expr* base, Expr* exp) { return mk_fn2("Power", base, exp); }
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }

static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}

/* ------------------------------------------------------------------ *
 *  Numeric predicates.                                                *
 * ------------------------------------------------------------------ */

/* True iff `e` is structurally a *concrete* zero.  Mirrors the helper
 * in solvepoly.c.  Used by the verification pass to decide whether the
 * back-substitution of a candidate proves it spurious. */
static bool is_definite_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        return is_definite_zero(e->data.function.args[0]);
    }
    return false;
}

/* If `e` is a *concrete* numeric leaf (Integer, Real, BigInt, Rational
 * with integer parts, MPFR), write |e| (in double) to *mag and return
 * true.  Returns false for symbolic / Complex / unhandled forms. */
static bool numeric_magnitude(const Expr* e, double* mag) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *mag = fabs((double)e->data.integer); return true; }
    if (e->type == EXPR_REAL)    { *mag = fabs(e->data.real); return true; }
    if (e->type == EXPR_BIGINT)  {
        double v = mpz_get_d(e->data.bigint);
        *mag = fabs(v);
        return true;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *mag = fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
        return true;
    }
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        double n, d;
        if (numeric_magnitude(e->data.function.args[0], &n)
            && numeric_magnitude(e->data.function.args[1], &d)
            && d != 0.0) {
            *mag = n / d;
            return true;
        }
    }
    /* Complex[re, im] — magnitude = sqrt(re^2 + im^2). */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2) {
        double re, im;
        if (numeric_magnitude(e->data.function.args[0], &re)
            && numeric_magnitude(e->data.function.args[1], &im)) {
            *mag = sqrt(re*re + im*im);
            return true;
        }
    }
    return false;
}

/* True iff the head of `e` is the symbol with the given interned name.
 * Callers should pass one of the cached SYM_* pointers from sym_names.h
 * -- this routine uses pointer equality against the global intern table
 * to keep the check O(1). */
static bool head_eq(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* ------------------------------------------------------------------ *
 *  Radical-atom detection.                                            *
 * ------------------------------------------------------------------ */

/* Return true iff `e` is a Power[base, p/q] node with q > 1 and `base`
 * depending on `var`.  On success populate *base_out (borrowed, NOT
 * deep-copied) and *p_out / *q_out. */
static bool is_radical_atom(const Expr* e, Expr* var,
                            Expr** base_out, int64_t* p_out, int64_t* q_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol.name != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* exp = e->data.function.args[1];
    int64_t p, q;
    if (!is_rational(exp, &p, &q)) return false;
    /* Canonical Rational has q > 0; we only treat the magnitude > 1
     * as "radical".  Integer exponents (q == 1) are ordinary
     * polynomial powers and must not be treated as radicals. */
    int64_t q_abs = q < 0 ? -q : q;
    if (q_abs <= 1) return false;
    Expr* b = e->data.function.args[0];
    if (!contains_any_symbol_from(b, var)) return false;
    *base_out = b;
    *p_out = p;
    *q_out = q;
    return true;
}

/* Walk `e` depth-first; on the first radical-atom hit return its base
 * (DEEP-COPIED into *base_out so the walker can free its inputs without
 * affecting the caller). */
static bool find_first_radical(Expr* e, Expr* var, Expr** base_out) {
    if (!e) return false;
    Expr* b; int64_t p, q;
    if (is_radical_atom(e, var, &b, &p, &q)) {
        *base_out = expr_copy(b);
        return true;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (find_first_radical(e->data.function.head, var, base_out)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (find_first_radical(e->data.function.args[i], var, base_out)) return true;
    }
    return false;
}

/* Walk `e`; for every Power[base, k/q] (structural equality on base),
 * fold q into the running LCM.  Returns the result. */
static int64_t collect_lcm_for_base(Expr* e, Expr* base) {
    if (!e) return 1;
    int64_t L = 1;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && expr_eq(e->data.function.args[0], base)) {
        int64_t p, q;
        if (is_rational(e->data.function.args[1], &p, &q)) {
            int64_t qa = q < 0 ? -q : q;
            if (qa > 1) L = lcm(L, qa);
        }
    }
    if (e->type == EXPR_FUNCTION) {
        L = lcm(L, collect_lcm_for_base(e->data.function.head, base));
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            L = lcm(L, collect_lcm_for_base(e->data.function.args[i], base));
        }
    }
    return L;
}

/* Replace every Power[base, p/q] with Power[u, p*L/q] in `e`.  Returns
 * a freshly-owned new tree.  `u_name` is the interned generator-symbol
 * name. */
static Expr* subst_radical_to_u(Expr* e, Expr* base, int64_t L, const char* u_name) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && expr_eq(e->data.function.args[0], base)) {
        int64_t p, q;
        if (is_rational(e->data.function.args[1], &p, &q)) {
            int64_t qa = q < 0 ? -q : q;
            if (qa > 1 && (L % qa) == 0) {
                /* k = p * L / q (preserving sign of p/q). */
                int64_t sign = (q < 0) ? -1 : 1;
                int64_t k = sign * p * (L / qa);
                if (k == 0) return mk_int(1);            /* u^0 == 1 */
                if (k == 1) return mk_sym(u_name);       /* u^1 == u */
                return mk_pow(mk_sym(u_name), mk_int(k));
            }
        }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = subst_radical_to_u(e->data.function.args[i], base, L, u_name);
    }
    Expr* new_head = subst_radical_to_u(e->data.function.head, base, L, u_name);
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    /* Evaluate so that newly-introduced u^k terms collapse cleanly
     * (Power[u, 2] * Power[u, 3] -> Power[u, 5], etc.) and so that
     * Plus / Times re-canonicalise. */
    return eval_and_free(out);
}

/* True iff `e` contains the bare symbol with the given interned name. */
static bool walk_contains_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type != EXPR_FUNCTION) return false;
    if (walk_contains_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (walk_contains_name(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Verification.                                                      *
 * ------------------------------------------------------------------ */

/* Verification result. */
typedef enum {
    VERIFY_ACCEPT,    /* candidate provably satisfies the equation */
    VERIFY_REJECT,    /* candidate provably fails */
    VERIFY_UNKNOWN    /* indeterminate (symbolic parameters) */
} VerifyResult;

/* Decide whether `cand` is a valid root of the *original* equation.
 *
 * Verification order is N[]-first, Simplify-fallback.  Simplify on an
 * algebraic-coefficient substitution (e.g. cand involving Sqrt[2] or
 * the imaginary part of a non-real quadratic root) can run for seconds
 * per candidate, while N[] on the same residual evaluates in
 * microseconds.  We only pay the Simplify cost when the numerical
 * pass cannot decide -- typically because the candidate carries free
 * parameters that survive substitution. */
static VerifyResult verify_candidate(Expr* e_orig, Expr* var, Expr* cand) {
    /* Root[] objects came out of the polynomial specialist with their
     * own verification semantics.  Accept them without further checks
     * -- they describe the unique algebraic root of an irreducible
     * polynomial factor and are not amenable to back-substitution. */
    if (head_eq(cand, SYM_Root)) return VERIFY_ACCEPT;

    /* Substitute var -> cand without Simplify; the bare evaluator
     * collapses anything obvious already, and N[] handles structural
     * cancellations (Sqrt[(Sqrt[3]-1)^2] etc.) numerically. */
    Expr* rule = mk_rule(expr_copy(var), expr_copy(cand));
    Expr* sub  = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(e_orig), rule }, 2));

    if (is_definite_zero(sub)) {
        expr_free(sub);
        return VERIFY_ACCEPT;
    }

    /* Numerical check.  N[sub] yields a concrete numeric leaf when
     * the candidate is closed-form (no free parameters).  Compare
     * magnitude against a generous absolute threshold -- the residual
     * is dimensionless in the equation's units so 1e-9 catches the
     * roundoff floor of double precision while comfortably rejecting
     * spurious roots, which typically miss by O(1). */
    Expr* nval = eval_and_free(mk_fn1("N", expr_copy(sub)));
    double mag;
    if (numeric_magnitude(nval, &mag)) {
        expr_free(nval);
        expr_free(sub);
        return (mag < 1.0e-9) ? VERIFY_ACCEPT : VERIFY_REJECT;
    }
    expr_free(nval);

    /* Slow path: numerical evaluation could not decide (free parameters,
     * Indeterminate from a removable singularity, etc.).  A symbolic
     * Simplify pass catches structural cancellations like
     * Sqrt[(Sqrt[3]-1)^2] -> Sqrt[3]-1 that the bare evaluator leaves
     * alone, then a final N[] retry handles parameter-free residuals
     * that only Simplify could unwrap. */
    Expr* simp = eval_and_free(mk_fn1("Simplify", sub));
    if (is_definite_zero(simp)) {
        expr_free(simp);
        return VERIFY_ACCEPT;
    }
    Expr* nval2 = eval_and_free(mk_fn1("N", expr_copy(simp)));
    bool concrete = numeric_magnitude(nval2, &mag);
    expr_free(nval2);
    expr_free(simp);

    if (concrete) {
        return (mag < 1.0e-9) ? VERIFY_ACCEPT : VERIFY_REJECT;
    }

    /* Residual still symbolic -- candidate carries free parameters
     * neither N[] nor Simplify could discharge.  Keep, with nongen
     * flagging. */
    return VERIFY_UNKNOWN;
}

/* ------------------------------------------------------------------ *
 *  Generator-variable management.                                     *
 * ------------------------------------------------------------------ */

typedef struct {
    char* name;       /* interned symbol pointer (intern_symbol) */
    Expr* side_eq;    /* u^L - base (post-substitution; polynomial in
                       * x and previously-introduced u's) */
} UVar;

typedef struct {
    UVar* items;
    size_t count;
    size_t cap;
} UVarList;

static void uvars_init(UVarList* L) { L->items = NULL; L->count = 0; L->cap = 0; }
static void uvars_push(UVarList* L, char* name, Expr* side_eq) {
    if (L->count == L->cap) {
        size_t nc = L->cap ? L->cap * 2 : 4;
        L->items = (UVar*)realloc(L->items, sizeof(UVar) * nc);
        L->cap = nc;
    }
    L->items[L->count].name = name;
    L->items[L->count].side_eq = side_eq;
    L->count++;
}
static void uvars_free(UVarList* L) {
    for (size_t i = 0; i < L->count; i++) {
        expr_free(L->items[i].side_eq);
        /* `name` points into the global intern table -- do NOT free. */
    }
    free(L->items);
    L->items = NULL;
    L->count = L->cap = 0;
}

/* Produce a fresh, never-before-used generator-symbol name.  Returns
 * the interned pointer (suitable for direct == comparison with symbol
 * names elsewhere in the system).  The chosen names follow the same
 * `$radu<n>$` template used by poly_make_fresh_gen, but namespaced
 * separately so we never clash with the polynomial GCD's helper
 * generators. */
static const char* fresh_u_name(Expr* witness1, Expr* witness2, int* counter) {
    char buf[32];
    for (;;) {
        snprintf(buf, sizeof(buf), "$radu%d$", (*counter)++);
        /* Intern first so we can do pointer-equality membership checks. */
        const char* interned = intern_symbol(buf);
        if (!walk_contains_name(witness1, interned)
            && !walk_contains_name(witness2, interned)) {
            return interned;
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Main entry.                                                        *
 * ------------------------------------------------------------------ */

/* Hard ceiling on the number of radical atoms we are willing to
 * eliminate.  Each pass adds one Resultant call, which is polynomial in
 * the inputs' Sylvester dimension -- big problems explode rapidly. */
#define SOLVERAD_MAX_GENS 12
/* Cap on the lcm of exponent denominators (= the radical substitution exponent
 * u^L).  The cleared polynomial degree scales with L, so an absurd denominator
 * such as x^(p/67890) would otherwise build a degree-67890 polynomial and hang.
 * Real radical equations use small denominators (2, 3, 4, 6, 12, ...). */
#define SOLVERAD_MAX_LCM 120

Expr* solverad_solve_radicals_equality(Expr* equation, Expr* var, Expr* dom) {
    if (!equation || !var) return NULL;
    if (var->type != EXPR_SYMBOL) return NULL;
    if (equation->type != EXPR_FUNCTION
        || equation->data.function.head->type != EXPR_SYMBOL
        || equation->data.function.head->data.symbol.name != SYM_Equal
        || equation->data.function.arg_count != 2) {
        return NULL;
    }

    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];

    /* e_orig = lhs - rhs.  Retained for verification at the tail. */
    Expr* e_orig = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs),
        mk_neg(expr_copy(rhs))));

    /* Together canonicalises rational radical expressions like
     *   (x + Sqrt[x]) / Sqrt[x] + Sqrt[x] / (x + Sqrt[x]) - 4
     * into a single fraction whose numerator we can hand to the
     * substitution pipeline.  The denominator is implicitly cleared by
     * cross-multiplying -- exactly the same trick solvepoly uses, with
     * the same caveat: extraneous roots at the denominator's zeros
     * will be filtered out at verification time. */
    Expr* e_together = eval_and_free(internal_together(
        (Expr*[]){ expr_copy(e_orig) }, 1));
    Expr* current = eval_and_free(internal_numerator(
        (Expr*[]){ expr_copy(e_together) }, 1));
    expr_free(e_together);

    /* Bail out cleanly if there are no radicals at all -- in that case
     * the polynomial specialist is the correct dispatch and we should
     * not steal the equation from it. */
    {
        Expr* probe = NULL;
        if (!find_first_radical(current, var, &probe)) {
            expr_free(current);
            expr_free(e_orig);
            return NULL;
        }
        expr_free(probe);
    }

    /* ------ Substitution loop. ------ */
    UVarList uvars; uvars_init(&uvars);
    int u_counter = 0;
    bool overflow = false;

    for (;;) {
        Expr* base = NULL;
        bool found_in_main = false;
        size_t found_in_side = (size_t)-1;

        /* Scan main first, then each side eq in introduction order. */
        if (find_first_radical(current, var, &base)) {
            found_in_main = true;
        } else {
            for (size_t i = 0; i < uvars.count; i++) {
                if (find_first_radical(uvars.items[i].side_eq, var, &base)) {
                    found_in_side = i;
                    break;
                }
            }
        }
        if (!base) break;
        (void)found_in_main; (void)found_in_side;

        if (uvars.count >= SOLVERAD_MAX_GENS) {
            expr_free(base);
            overflow = true;
            break;
        }

        /* L = lcm of all exponent denominators across every equation
         * that contains a Power[base, _] subterm. */
        int64_t L = collect_lcm_for_base(current, base);
        for (size_t i = 0; i < uvars.count; i++) {
            L = lcm(L, collect_lcm_for_base(uvars.items[i].side_eq, base));
        }
        if (L <= 1) {
            /* Shouldn't happen -- is_radical_atom rejects q<=1.  Defensive. */
            expr_free(base);
            overflow = true;
            break;
        }
        /* Guard against an exponent-denominator explosion: the substitution
         * u^L makes the cleared polynomial degree scale with L, so a base like
         * x^(123451/67890) (L = 67890) would build an astronomically high-degree
         * polynomial and hang the resultant/root step.  Real radical equations
         * use small denominators; bail (leave Solve unevaluated) past the cap. */
        if (L > SOLVERAD_MAX_LCM) {
            expr_free(base);
            overflow = true;
            break;
        }

        const char* u_name = fresh_u_name(current, e_orig, &u_counter);

        /* Substitute Power[base, p/q] -> u^(p*L/q) everywhere. */
        Expr* new_main = subst_radical_to_u(current, base, L, u_name);
        expr_free(current);
        current = new_main;
        for (size_t i = 0; i < uvars.count; i++) {
            Expr* ns = subst_radical_to_u(uvars.items[i].side_eq, base, L, u_name);
            expr_free(uvars.items[i].side_eq);
            uvars.items[i].side_eq = ns;
        }

        /* New side eq: u^L - base.  Note: `base` here was captured BEFORE
         * the substitution pass above, so it may still hold radicals
         * (the inner ones for nested cases).  Apply the same
         * substitution we just applied to `current` so future loop
         * iterations can find the now-fresh inner atoms. */
        Expr* base_subst = subst_radical_to_u(base, base, L, u_name);
        /* In normal use base_subst != base only if base structurally
         * matched itself as a power, which is impossible (base is the
         * raw base of a Power node, not a Power node).  Keep the
         * defensive subst to handle pathological inputs. */

        Expr* new_side = eval_and_free(mk_fn2("Plus",
            mk_pow(mk_sym(u_name), mk_int(L)),
            mk_neg(base_subst)));
        char* name_cstr = (char*)u_name;  /* interned */
        uvars_push(&uvars, name_cstr, new_side);
        expr_free(base);
    }

    if (overflow) {
        uvars_free(&uvars);
        expr_free(current);
        expr_free(e_orig);
        return NULL;
    }

    /* ------ Resultant chain. ------ *
     *
     * `current` is now a polynomial in `var` and every u_i; each u_i
     * has a side equation that is polynomial in `var` and a *subset* of
     * the u_j's introduced earlier in the loop.  Eliminate them in
     * introduction order: side_eq[i] is linear-or-higher in u_i and
     * polynomial in the rest, so Resultant_{u_i}(current, side_eq[i],
     * u_i) drops u_i while introducing whatever u_j's appear in
     * side_eq[i] (already in `current` anyway).  After processing every
     * u_i in this way the result is in `var` alone. */
    Expr* poly = eval_and_free(internal_expand((Expr*[]){ current }, 1));
    current = NULL;

    for (size_t i = 0; i < uvars.count; i++) {
        Expr* uvar = mk_sym(uvars.items[i].name);
        Expr* side = eval_and_free(internal_expand(
            (Expr*[]){ expr_copy(uvars.items[i].side_eq) }, 1));
        /* internal_resultant takes ownership of its arg array entries. */
        Expr* r = internal_resultant(
            (Expr*[]){ poly, side, uvar }, 3);
        if (!r) {
            uvars_free(&uvars);
            expr_free(e_orig);
            return NULL;
        }
        poly = eval_and_free(internal_expand((Expr*[]){ r }, 1));
        /* Sanity: if Resultant degenerated to a non-zero constant, the
         * system is inconsistent (no solutions).  Return the empty
         * list rather than handing a constant equation to solvepoly. */
        if (poly && poly->type == EXPR_FUNCTION
            && contains_any_symbol_from(poly, var) == false) {
            if (is_definite_zero(poly)) {
                /* 0 == 0: every x is a solution -- but we never see
                 * this for a finite radical system; pass through to
                 * solvepoly which knows how to render it as {{}} too. */
            } else {
                /* Concrete non-zero: empty solution set. */
                expr_free(poly);
                uvars_free(&uvars);
                expr_free(e_orig);
                return mk_list(NULL, 0);
            }
        } else if (poly && (poly->type == EXPR_INTEGER || poly->type == EXPR_REAL
                            || poly->type == EXPR_BIGINT)) {
            if (is_definite_zero(poly)) {
                /* fall through to solvepoly */
            } else {
                expr_free(poly);
                uvars_free(&uvars);
                expr_free(e_orig);
                return mk_list(NULL, 0);
            }
        }
    }
    uvars_free(&uvars);

    /* ------ Polynomial solve. ------ */
    Expr* eqn = mk_fn2("Equal", poly, mk_int(0));
    SolvePolyOpts opts = { false, false };
    Expr* candidates = solvepoly_solve_polynomial_equality(eqn, var, dom, &opts);
    expr_free(eqn);

    if (!candidates) {
        expr_free(e_orig);
        return NULL;
    }
    if (candidates->type != EXPR_FUNCTION
        || candidates->data.function.head->type != EXPR_SYMBOL
        || candidates->data.function.head->data.symbol.name != SYM_List) {
        /* Should never happen, but be defensive. */
        expr_free(e_orig);
        return candidates;
    }

    /* ------ Verification + dedupe. ------ *
     *
     * The resultant chain can introduce extraneous roots in two ways:
     *
     *   (1) cross-multiplication by the Together denominator
     *       (denominator zeros become spurious solutions);
     *   (2) the substitution u^L = base hands every L-th root of `base`
     *       to the polynomial solver, while only the principal branch
     *       satisfies the original Sqrt / x^(1/L) on the real line.
     *
     * Per-candidate verification rejects (1) and (2) numerically;
     * candidates that depend on free parameters we cannot decide
     * symbolically are kept and trigger Solve::nongen, matching
     * Mathematica's convention.  We also deduplicate the surviving
     * candidates (structural equality) because the resultant can
     * produce repeated factors -- e.g. an irrelevant x^2 factor in the
     * rational test case.  Multiplicity > 1 in a *valid* root would
     * already have been preserved by solvepoly's Factor / multiplicity
     * pipeline. */
    size_t nc = candidates->data.function.arg_count;
    Expr** good = nc ? (Expr**)malloc(sizeof(Expr*) * nc) : NULL;
    size_t kept = 0;
    bool any_nongen = false;

    for (size_t i = 0; i < nc; i++) {
        Expr* sol = candidates->data.function.args[i];
        /* sol is expected to be List[Rule[var, val]]. */
        if (sol->type != EXPR_FUNCTION
            || !head_eq(sol, SYM_List)
            || sol->data.function.arg_count != 1
            || !head_eq(sol->data.function.args[0], SYM_Rule)) {
            /* Unusual shape -- pass through unverified, conservatively. */
            good[kept++] = expr_copy(sol);
            any_nongen = true;
            continue;
        }
        Expr* rule = sol->data.function.args[0];
        Expr* val  = rule->data.function.args[1];

        VerifyResult v = verify_candidate(e_orig, var, val);
        if (v == VERIFY_REJECT) continue;

        /* Dedupe against already-kept candidates. */
        bool dup = false;
        for (size_t j = 0; j < kept; j++) {
            if (expr_eq(good[j], sol)) { dup = true; break; }
        }
        if (dup) continue;

        if (v == VERIFY_UNKNOWN) any_nongen = true;
        good[kept++] = expr_copy(sol);
    }

    if (any_nongen) {
        fprintf(stderr,
            "Solve::nongen: There may be values of the parameters for which "
            "some or all solutions are not valid.\n");
    }

    expr_free(candidates);
    expr_free(e_orig);

    return mk_list(good, kept);
}

/* ------------------------------------------------------------------ *
 *  Builtin entry & init.                                              *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_radicals_equality(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* dom      = (argc >= 3) ? res->data.function.args[2] : NULL;
    return solverad_solve_radicals_equality(equation, var, dom);
}

void solverad_init(void) {
    symtab_add_builtin("Solve`SolveRadicalsEquality",
                       builtin_solve_radicals_equality);
    SymbolDef* def = symtab_get_def("Solve`SolveRadicalsEquality");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveRadicalsEquality",
        "Solve`SolveRadicalsEquality[lhs == rhs, var]\n"
        "Solve`SolveRadicalsEquality[lhs == rhs, var, dom]\n"
        "\tThe radical-equation specialist used by Solve.  Handles\n"
        "\tequations containing Sqrt and fractional powers x^(p/q),\n"
        "\tincluding multi-base and nested cases such as\n"
        "\tSqrt[x + 1] + Sqrt[x - 1] == 3 and\n"
        "\tSqrt[x + Sqrt[x + Sqrt[x]]] == 2.\n"
        "\n"
        "\tIntroduces fresh generators u_i for each distinct radical\n"
        "\tbase, builds the elimination system, chains Resultant calls\n"
        "\tto reduce to a polynomial in var, and verifies each candidate\n"
        "\tby back-substitution into the original equation.  Candidates\n"
        "\tthat depend on free parameters and cannot be verified\n"
        "\tsymbolically trigger Solve::nongen.\n"
        "\n"
        "\tReturns NULL (i.e. leaves the call unevaluated) when the\n"
        "\tequation has no radical atom in var; the polynomial\n"
        "\tspecialist is the correct dispatch in that case.");
}
