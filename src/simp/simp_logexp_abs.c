#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* Log/Power rewriter (positive-real cascade, v1)                          */
/* ----------------------------------------------------------------------- */

/*
 * The strict-positive cascade implements the Log/Power identities that
 * are sound under positivity / reality assumptions on the operands.
 * Identities cover (1) log of products and quotients, (2) log of a power
 * of a positive base, (3) power of a product, and (4) tower-of-powers
 * collapse for a positive base.
 *
 * The general-real and general-complex branches of the user's cascade
 * (with Boole / Floor / Ceiling phase corrections) are deliberately not
 * implemented in v1; see Mathilda_spec.md for v2 scope.
 *
 * Implementation: a bottom-up structural walker that consults the
 * AssumeCtx for positivity/reality of operands. Each top-level rewrite
 * emits a freshly evaluated tree, so e.g. nested Log[Times[x, 1/y]] ->
 * Log[x] + Log[1/y] -> Log[x] - Log[y] (via the Power[..., -1] case)
 * stabilises after a small fixed number of passes.
 */

static Expr* logexp_top_rewrite(const Expr* e, const AssumeCtx* ctx);

/* Returns NULL if the recursive walk produced no change. Otherwise returns
 * a newly owned, evaluated tree. */
static Expr* logexp_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) {
        return logexp_top_rewrite(e, ctx);
    }

    /* First rewrite children. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = NULL;
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        Expr* r = logexp_walk(e->data.function.args[i], ctx);
        if (r) {
            if (!new_args) {
                new_args = (Expr**)calloc(n, sizeof(Expr*));
                for (size_t j = 0; j < i; j++) new_args[j] = expr_copy(e->data.function.args[j]);
            }
            new_args[i] = r;
            any = true;
        } else if (new_args) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* current_owned = NULL;
    const Expr* target;
    if (any) {
        Expr* head_copy = expr_copy(e->data.function.head);
        Expr* rebuilt = expr_new_function(head_copy, new_args, n);
        free(new_args);
        current_owned = evaluate(rebuilt);
        expr_free(rebuilt);
        target = current_owned;
    } else {
        target = e;
    }

    Expr* top = logexp_top_rewrite(target, ctx);
    if (top) {
        if (current_owned) expr_free(current_owned);
        return top;
    }
    return current_owned;  /* may be NULL if no change anywhere */
}

static Expr* build_unary(const char* head, Expr* owned_arg) {
    Expr* a[1] = { owned_arg };
    return expr_new_function(expr_new_symbol(head), a, 1);
}

static Expr* build_binary(const char* head, Expr* a0, Expr* a1) {
    Expr* a[2] = { a0, a1 };
    return expr_new_function(expr_new_symbol(head), a, 2);
}

static Expr* logexp_top_rewrite(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol;
    Expr** a = e->data.function.args;
    size_t n = e->data.function.arg_count;

    /* Log[Times[u1,...,un]] -> Sum Log[ui]  when every ui is positive.
     * Log[Power[x, p]]      -> p Log[x]      when x positive and p real. */
    if (h == SYM_Log && n == 1) {
        Expr* inner = a[0];
        if (inner->type == EXPR_FUNCTION &&
            inner->data.function.head &&
            inner->data.function.head->type == EXPR_SYMBOL) {
            const char* ih = inner->data.function.head->data.symbol;
            size_t in = inner->data.function.arg_count;
            Expr** ia = inner->data.function.args;

            if (ih == SYM_Times && in > 0) {
                bool all_pos = true;
                for (size_t i = 0; i < in; i++) {
                    if (!prov_pos(ctx, ia[i])) { all_pos = false; break; }
                }
                if (all_pos) {
                    Expr** logs = (Expr**)calloc(in, sizeof(Expr*));
                    for (size_t i = 0; i < in; i++) {
                        logs[i] = build_unary("Log", expr_copy(ia[i]));
                    }
                    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), logs, in);
                    free(logs);
                    Expr* canon = evaluate(sum);
                    expr_free(sum);
                    return canon;
                }
            }
            if (ih == SYM_Power && in == 2) {
                Expr* base = ia[0];
                Expr* p    = ia[1];
                if (prov_pos(ctx, base) && prov_re(ctx, p)) {
                    Expr* logx = build_unary("Log", expr_copy(base));
                    Expr* mul  = build_binary("Times", expr_copy(p), logx);
                    Expr* canon = evaluate(mul);
                    expr_free(mul);
                    return canon;
                }
            }
        }
    }

    /* Power[-1, k] reductions: even exponent -> 1, odd exponent -> -1.
     * Catches `(-1)^(2 n)` with `n` provably integer, etc. The Times
     * propagation in prov_even (a product of integers with at least one
     * even factor is even) handles the canonical user input. */
    if (h == SYM_Power && n == 2) {
        Expr* base = a[0];
        Expr* exp_  = a[1];
        if (base->type == EXPR_INTEGER && base->data.integer == -1) {
            if (prov_even(ctx, exp_)) return expr_new_integer(1);
        }
    }

    /* Exp distribute: Power[E, Plus[t1,...,tn]] -> Product Power[E, ti].
     * Always sound (E^(x+y) = E^x · E^y for any x, y). When the exponent
     * is a product like `Times[c, Plus[...]]`, expand it first so the
     * Plus surfaces. After distribution, individual Power[E, c·Log u]
     * subterms collapse to `u^c` via Mathilda's existing
     * Power[E, c·Log[u]] -> u^c rule, completing identities like
     *
     *   Exp[3 (Log[a] + Log[b])]  ->  a^3 · b^3
     *   Exp[y (Log[a] + Log[b])]  ->  a^y · b^y    (a > 0, b > 0)
     */
    if (h == SYM_Power && n == 2) {
        Expr* base = a[0];
        Expr* exp_  = a[1];

        if (base->type == EXPR_SYMBOL && base->data.symbol == SYM_E) {
            /* Try to expand the exponent so Plus distributes through any
             * outer Times. Cheap, idempotent, and only used to detect
             * Plus structure. */
            Expr* expanded_exp = expr_expand(exp_);
            if (expanded_exp &&
                expanded_exp->type == EXPR_FUNCTION &&
                expanded_exp->data.function.head &&
                expanded_exp->data.function.head->type == EXPR_SYMBOL &&
                expanded_exp->data.function.head->data.symbol == SYM_Plus &&
                expanded_exp->data.function.arg_count > 1) {
                size_t en = expanded_exp->data.function.arg_count;
                Expr** factors = (Expr**)calloc(en, sizeof(Expr*));
                for (size_t i = 0; i < en; i++) {
                    factors[i] = build_binary("Power", expr_new_symbol(SYM_E),
                                              expr_copy(expanded_exp->data.function.args[i]));
                }
                Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), factors, en);
                free(factors);
                Expr* canon = evaluate(prod);
                expr_free(prod);
                expr_free(expanded_exp);
                return canon;
            }
            if (expanded_exp) expr_free(expanded_exp);
        }

        if (base->type == EXPR_FUNCTION &&
            base->data.function.head &&
            base->data.function.head->type == EXPR_SYMBOL) {
            const char* bh = base->data.function.head->data.symbol;
            size_t bn = base->data.function.arg_count;
            Expr** ba = base->data.function.args;

            if (bh == SYM_Times && bn > 0) {
                bool all_pos = true;
                for (size_t i = 0; i < bn; i++) {
                    if (!prov_pos(ctx, ba[i])) { all_pos = false; break; }
                }
                if (all_pos) {
                    Expr** powers = (Expr**)calloc(bn, sizeof(Expr*));
                    for (size_t i = 0; i < bn; i++) {
                        powers[i] = build_binary("Power", expr_copy(ba[i]), expr_copy(exp_));
                    }
                    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), powers, bn);
                    free(powers);
                    Expr* canon = evaluate(prod);
                    expr_free(prod);
                    return canon;
                }
            }
            if (bh == SYM_Power && bn == 2) {
                Expr* xx = ba[0];
                Expr* pp = ba[1];
                if (prov_pos(ctx, xx) && prov_re(ctx, pp)) {
                    Expr* prod = build_binary("Times", expr_copy(pp), expr_copy(exp_));
                    Expr* prod_canon = evaluate(prod);
                    expr_free(prod);
                    Expr* pow_ = build_binary("Power", expr_copy(xx), prod_canon);
                    Expr* canon = evaluate(pow_);
                    expr_free(pow_);
                    return canon;
                }
            }
        }
    }

    return NULL;
}

/* Apply the rewriter to a fixed point. Returns NULL if unchanged.
 * Bounded iteration count protects against pathological alternations
 * with the evaluator's canonicalisation. */
Expr* apply_logexp_rules(const Expr* input, const AssumeCtx* ctx) {
    /* NULL ctx is treated as an empty context. The positivity-aware
     * Log[a*b] -> Log[a]+Log[b] etc. rewrites simply won't fire
     * (prov_pos returns false), but the unconditional
     * Power[E, Plus[...]] distribute rule still does its job. */
    Expr* current = expr_copy((Expr*)input);
    bool changed = false;
    for (int iter = 0; iter < 8; iter++) {
        Expr* r = logexp_walk(current, ctx);
        if (!r) break;
        if (expr_eq(r, current)) { expr_free(r); break; }
        expr_free(current);
        current = r;
        changed = true;
    }
    if (!changed) {
        expr_free(current);
        return NULL;
    }
    if (expr_eq(current, input)) {
        expr_free(current);
        return NULL;
    }
    return current;
}

/* ----------------------------------------------------------------------- */
/* Abs simplification: structural rewrites over Abs[...] subexpressions   */
/* ----------------------------------------------------------------------- */

/* Cheap pre-check: skip the walker when the input is Abs-free. */
bool contains_abs(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Abs) return true;
    if (contains_abs(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_abs(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff `h` is a 1-arg head whose presence makes a transform that
 * targets trig or hyperbolic functions potentially fire. Covers the six
 * canonical pairs and their inverses. */
static bool head_is_trig_or_hyperbolic(const char* h) {
    static const char* const TRIG_HEADS[] = {
        "Sin","Cos","Tan","Cot","Sec","Csc",
        "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
        "Sinh","Cosh","Tanh","Coth","Sech","Csch",
        "ArcSinh","ArcCosh","ArcTanh","ArcCoth","ArcSech","ArcCsch",
        NULL
    };
    for (size_t i = 0; TRIG_HEADS[i]; i++) {
        if (strcmp(h, TRIG_HEADS[i]) == 0) return true;
    }
    return false;
}

bool contains_trig_or_hyperbolic(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        head_is_trig_or_hyperbolic(e->data.function.head->data.symbol)) return true;
    if (contains_trig_or_hyperbolic(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_trig_or_hyperbolic(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff `e` contains a Power[E, _] subexpression -- i.e. an Exp atom in
 * exponential form (E^x, E^(-x), E^(I Pi), ...). Used to gate the
 * ExpToTrig seed: pure-Exp inputs miss every trig-gated transform
 * (TrigRoundtrip, PythagReduce, TrigFactor, ...) because their gates
 * check for Cos/Sin/Cosh/Sinh heads. ExpToTrig converts the Exp form
 * to Cosh/Sinh, opening those transforms to subsequent rounds. */
bool contains_exp_form(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base && base->type == EXPR_SYMBOL && base->data.symbol == SYM_E) {
            return true;
        }
    }
    if (contains_exp_form(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_exp_form(e->data.function.args[i])) return true;
    }
    return false;
}

bool contains_log(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Log) return true;
    if (contains_log(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_log(e->data.function.args[i])) return true;
    }
    return false;
}

bool contains_power(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power) return true;
    if (contains_power(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_power(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff a Plus, Times, or Power head appears anywhere in `e`. Used by
 * the TrigReduce gate to short-circuit on a bare single trig call (no
 * product or power means no product-to-sum work). Power is included
 * because Sin[x]^2 is the canonical TrigReduce input. */
bool contains_plus_or_times(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Plus ||
            h == SYM_Times ||
            h == SYM_Power) return true;
    }
    if (contains_plus_or_times(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_plus_or_times(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff any non-numeric-constant symbol leaf appears anywhere in `e`.
 * Pi, E, EulerGamma, Degree, Catalan, Glaisher, Khinchin do not count --
 * they are positive numeric constants. Used to short-circuit transforms
 * that have nothing to do on a purely numeric input (Factor, Apart, ...). */
bool contains_variable(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        return !is_real_constant_symbol(e->data.symbol);
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_variable(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_variable(e->data.function.args[i])) return true;
    }
    return false;
}

/* Number of distinct non-constant symbol leaves in `e`, capped at `cap`.
 * Returns as soon as the count reaches `cap`, so callers that only need
 * "0 / 1 / >=2" can pass cap=2 and early-out. Constant symbols (Pi, E,
 * ...) are excluded, matching contains_variable. */
static size_t expr_variables_count_capped_walk(const Expr* e,
                                               char** seen, size_t* nseen,
                                               size_t cap) {
    if (!e || *nseen >= cap) return *nseen;
    if (e->type == EXPR_SYMBOL) {
        if (is_real_constant_symbol(e->data.symbol)) return *nseen;
        for (size_t i = 0; i < *nseen; i++) {
            if (strcmp(seen[i], e->data.symbol) == 0) return *nseen;
        }
        seen[*nseen] = e->data.symbol;
        (*nseen)++;
        return *nseen;
    }
    if (e->type != EXPR_FUNCTION) return *nseen;
    expr_variables_count_capped_walk(e->data.function.head, seen, nseen, cap);
    for (size_t i = 0; i < e->data.function.arg_count && *nseen < cap; i++) {
        expr_variables_count_capped_walk(e->data.function.args[i],
                                         seen, nseen, cap);
    }
    return *nseen;
}

size_t expr_variables_count_capped(const Expr* e, size_t cap) {
    if (cap == 0) return 0;
    char* seen[8];  /* cap is at most 2 in our call sites; 8 is a safe ceiling */
    size_t nseen = 0;
    if (cap > 8) cap = 8;
    expr_variables_count_capped_walk(e, seen, &nseen, cap);
    return nseen;
}

/* True iff the assumption ctx has at least one usable fact. NULL ctx, an
 * empty fact list, or an inconsistent ctx all return false -- no
 * assumption-driven rewrite can do anything in those cases. */
bool ctx_has_facts(const AssumeCtx* ctx) {
    return ctx != NULL && ctx->count > 0 && !ctx->inconsistent;
}

/* Try to simplify a single Abs[arg] node. `arg` is the inner expression
 * (i.e. the argument to Abs). Returns a new Expr* on success, NULL if no
 * rule fires. */
static Expr* try_simp_abs(const Expr* arg, const AssumeCtx* ctx) {
    /* Universal: idempotency Abs[Abs[x]] -> Abs[x]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Abs &&
        arg->data.function.arg_count == 1) {
        return expr_copy((Expr*)arg);
    }

    /* Universal: conjugate symmetry Abs[Conjugate[x]] -> Abs[x]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Conjugate &&
        arg->data.function.arg_count == 1) {
        Expr* a[1] = { expr_copy(arg->data.function.args[0]) };
        return expr_new_function(expr_new_symbol(SYM_Abs), a, 1);
    }

    /* Universal: Abs[E^z] -> E^Re[z]. The magnitude of any complex
     * exponential is e^(real part of the exponent). */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        arg->data.function.args[0]->type == EXPR_SYMBOL &&
        arg->data.function.args[0]->data.symbol == SYM_E) {
        Expr* re_in[1] = { expr_copy(arg->data.function.args[1]) };
        Expr* re_call = expr_new_function(expr_new_symbol(SYM_Re), re_in, 1);
        Expr* pa[2] = { expr_new_symbol(SYM_E), re_call };
        return expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    }

    /* Universal: split products. Abs[Times[a, b, ...]] -> Abs[a] Abs[b] ...
     * Captures both Abs[c x] (numeric coefficient extraction) and the
     * Abs[x/y] case since x/y is Times[x, Power[y, -1]]. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Times &&
        arg->data.function.arg_count >= 2) {
        size_t n = arg->data.function.arg_count;
        Expr** new_args = (Expr**)calloc(n, sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            Expr* a[1] = { expr_copy(arg->data.function.args[i]) };
            new_args[i] = expr_new_function(expr_new_symbol(SYM_Abs), a, 1);
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_Times), new_args, n);
        free(new_args);
        return result;
    }

    /* Universal: integer-power split. Abs[x^n] -> Abs[x]^n for integer n.
     * For complex x and integer n the identity |x^n| = |x|^n is exact;
     * for non-integer n it can fail (branch-cut), so the unconditional
     * rule applies only to integer exponents. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        (arg->data.function.args[1]->type == EXPR_INTEGER ||
         arg->data.function.args[1]->type == EXPR_BIGINT)) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* a[1] = { expr_copy(base) };
        Expr* abs_call = expr_new_function(expr_new_symbol(SYM_Abs), a, 1);
        Expr* pa[2] = { abs_call, expr_copy(exp) };
        return expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    }

    /* The remaining rules need an assumption context. */
    if (!ctx) return NULL;

    /* Cascading: Abs[x] -> x  if x >= 0 (provably nonnegative). */
    if (assume_known_nonneg(ctx, arg)) {
        return expr_copy((Expr*)arg);
    }

    /* Cascading: Abs[x] -> -x  if x <= 0 (provably nonpositive). */
    if (assume_known_nonpos(ctx, arg)) {
        Expr* na[2] = { expr_new_integer(-1), expr_copy((Expr*)arg) };
        return expr_new_function(expr_new_symbol(SYM_Times), na, 2);
    }

    /* Cascading: Abs[x^y] -> Abs[x]^y if y is real. The integer-power
     * rule above handles n in Z; this generalises to any real y under
     * an Element[y, Reals] assumption. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        assume_known_real(ctx, arg->data.function.args[1])) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* a[1] = { expr_copy(base) };
        Expr* abs_call = expr_new_function(expr_new_symbol(SYM_Abs), a, 1);
        Expr* pa[2] = { abs_call, expr_copy(exp) };
        return expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    }

    /* Cascading: Abs[x^y] -> x^Re[y] if x > 0 (strictly positive).
     * Proof: for x > 0, x^y = x^(Re[y] + I Im[y]) = x^Re[y] * Exp[I Im[y]
     * Log[x]] and the second factor has unit modulus. */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_Power &&
        arg->data.function.arg_count == 2 &&
        assume_known_positive(ctx, arg->data.function.args[0])) {
        Expr* base = arg->data.function.args[0];
        Expr* exp  = arg->data.function.args[1];
        Expr* re_in[1] = { expr_copy(exp) };
        Expr* re_call = expr_new_function(expr_new_symbol(SYM_Re), re_in, 1);
        Expr* pa[2] = { expr_copy(base), re_call };
        return expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    }

    /* The Abs[Sin[x]] -> Sign[Sin[x]] Sin[x] rule from the user-provided
     * cascade is omitted: the rewrite expands leaf count (3 -> 6) and only
     * pays off when a downstream Sign-folding pass narrows Sign[Sin[x]] on
     * a known interval, which Mathilda does not currently perform. Adding
     * it without that infrastructure produces a strictly larger expression
     * with no observable benefit. */
    return NULL;
}

/* Bottom-up walker that rewrites Abs[...] subexpressions. Returns a new
 * Expr* if any rewrite fired anywhere in the tree, NULL otherwise. */
static Expr* abs_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;

    size_t argc = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(argc ? argc : 1, sizeof(Expr*));
    bool any = false;
    for (size_t i = 0; i < argc; i++) {
        Expr* sub = abs_walk(e->data.function.args[i], ctx);
        if (sub) {
            new_args[i] = sub;
            any = true;
        } else {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }

    Expr* this_form;
    if (any) {
        this_form = expr_new_function(expr_copy(e->data.function.head),
                                       new_args, argc);
    } else {
        for (size_t i = 0; i < argc; i++) expr_free(new_args[i]);
        this_form = NULL;
    }
    free(new_args);

    /* Rule fires only on Abs[_]. */
    bool is_abs = e->data.function.head &&
                  e->data.function.head->type == EXPR_SYMBOL &&
                  e->data.function.head->data.symbol == SYM_Abs &&
                  e->data.function.arg_count == 1;
    if (is_abs) {
        const Expr* inner = this_form ? this_form->data.function.args[0]
                                       : e->data.function.args[0];
        Expr* simp = try_simp_abs(inner, ctx);
        if (simp) {
            if (this_form) expr_free(this_form);
            return simp;
        }
    }

    return this_form;
}

/* Returns a rewritten copy of `input` if any Abs simplification fired,
 * NULL otherwise. ctx may be NULL (universal rules still fire). */
Expr* apply_abs_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!contains_abs(input)) return NULL;
    return abs_walk(input, ctx);
}

/* ----------------------------------------------------------------------- */
/* Sqrt[e^2] simplification                                                */
/*                                                                          */
/* Companion walker for Sqrt[_^2] (i.e. Power[Power[base, 2], Rational[1,2]])
 * subexpressions. The per-symbol rule cluster in apply_assumption_rules
 * only handles bare symbols; this walker covers general subexpressions and
 * relies on the prov_pos / prov_neg numeric-sign fallback to decide the
 * sign of `base` for closed-form numeric expressions like
 * `1 + Cos[2] + Cos[3] + Sqrt[1 + Sqrt[3]]`.                              */
/* ----------------------------------------------------------------------- */

/* True iff e is the literal `Rational[1, 2]`. */
static bool is_rational_half(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1
        && e->data.function.args[1]->type == EXPR_INTEGER
        && e->data.function.args[1]->data.integer == 2;
}

/* True iff e has the shape Power[Power[_, 2], Rational[1, 2]] = Sqrt[X^2]. */
static bool is_sqrt_of_square(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol != SYM_Power ||
        e->data.function.arg_count != 2) return false;
    Expr* outer_base = e->data.function.args[0];
    Expr* outer_exp = e->data.function.args[1];
    if (!is_rational_half(outer_exp)) return false;
    if (outer_base->type != EXPR_FUNCTION ||
        !outer_base->data.function.head ||
        outer_base->data.function.head->type != EXPR_SYMBOL ||
        outer_base->data.function.head->data.symbol != SYM_Power ||
        outer_base->data.function.arg_count != 2) return false;
    Expr* inner_exp = outer_base->data.function.args[1];
    return inner_exp->type == EXPR_INTEGER && inner_exp->data.integer == 2;
}

bool contains_sqrt_of_square(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_sqrt_of_square(e)) return true;
    if (contains_sqrt_of_square(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_sqrt_of_square(e->data.function.args[i])) return true;
    }
    return false;
}

/* Rewrite Power[Power[base, 2], Rational[1, 2]] (= Sqrt[base^2]) when the
 * sign of `base` is determinable. Returns a fresh Expr on success, NULL
 * if nothing fires (so the outer caller can keep the unsimplified form). */
static Expr* try_simp_sqrt_of_square(const Expr* sqrt_node, const AssumeCtx* ctx) {
    if (!is_sqrt_of_square(sqrt_node)) return NULL;
    const Expr* base = sqrt_node->data.function.args[0]->data.function.args[0];
    if (prov_nn(ctx, base)) {
        /* base >= 0 → Sqrt[base^2] = base. */
        return expr_copy((Expr*)base);
    }
    if (prov_np(ctx, base)) {
        /* base <= 0 → Sqrt[base^2] = -base. */
        Expr* na[2] = { expr_new_integer(-1), expr_copy((Expr*)base) };
        return expr_new_function(expr_new_symbol(SYM_Times), na, 2);
    }
    if (prov_re(ctx, base)) {
        /* base real but sign undetermined → Sqrt[base^2] = Abs[base].
         * Downstream Abs simplification may reduce further (or stop here
         * — Abs[real] is at least as canonical as Sqrt[real^2]). */
        Expr* a[1] = { expr_copy((Expr*)base) };
        return expr_new_function(expr_new_symbol(SYM_Abs), a, 1);
    }
    return NULL;
}

static Expr* sqrt_of_square_walk(const Expr* e, const AssumeCtx* ctx) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    size_t argc = e->data.function.arg_count;
    Expr** new_args = (Expr**)calloc(argc ? argc : 1, sizeof(Expr*));
    bool any = false;
    for (size_t i = 0; i < argc; i++) {
        Expr* sub = sqrt_of_square_walk(e->data.function.args[i], ctx);
        if (sub) {
            new_args[i] = sub;
            any = true;
        } else {
            new_args[i] = expr_copy(e->data.function.args[i]);
        }
    }
    Expr* this_form;
    if (any) {
        this_form = expr_new_function(expr_copy(e->data.function.head),
                                       new_args, argc);
    } else {
        for (size_t i = 0; i < argc; i++) expr_free(new_args[i]);
        this_form = NULL;
    }
    free(new_args);
    const Expr* candidate = this_form ? this_form : e;
    if (is_sqrt_of_square(candidate)) {
        Expr* simp = try_simp_sqrt_of_square(candidate, ctx);
        if (simp) {
            if (this_form) expr_free(this_form);
            return simp;
        }
    }
    return this_form;
}

/* Returns a rewritten copy of `input` if any Sqrt[_^2] rewrite fired,
 * NULL otherwise. */
Expr* apply_sqrt_of_square_rules(const Expr* input, const AssumeCtx* ctx) {
    if (!contains_sqrt_of_square(input)) return NULL;
    return sqrt_of_square_walk(input, ctx);
}

