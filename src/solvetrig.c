/*
 * solvetrig.c
 *
 * Trig canonicalisation pre-pass: see solvetrig.h for the algorithm
 * narrative.  Reuses the solveinv specialist's peel_exp machinery to
 * unwind each u-root into the periodic var-family.
 */

#include "solvetrig.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly.h"
#include "solveinv.h"
#include "solvepoly.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (mirrors solveinv.c style).             *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_pow(Expr* base, Expr* exp) {
    return mk_fn2("Power", base, exp);
}
static Expr* mk_neg(Expr* e) {
    return mk_fn2("Times", mk_int(-1), e);
}

static bool head_is_sym(const Expr* e, const char* interned_head) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == interned_head;
}

/* ------------------------------------------------------------------ *
 *  Trig-head detection: structural walk.                             *
 * ------------------------------------------------------------------ */

static bool is_trig_head_name(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Tan
        || h == SYM_Cot || h == SYM_Sec || h == SYM_Csc
        || h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh
        || h == SYM_Coth || h == SYM_Sech || h == SYM_Csch;
}

static bool var_in(const Expr* e, const Expr* var) {
    if (!e || !var) return false;
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (var_in(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (var_in(e->data.function.args[i], var)) return true;
    }
    return false;
}

bool solvetrig_has_trig(const Expr* expr, const Expr* var) {
    if (!expr || !var) return false;
    if (expr->type == EXPR_FUNCTION) {
        if (expr->data.function.head->type == EXPR_SYMBOL) {
            const char* h = expr->data.function.head->data.symbol.name;
            if (is_trig_head_name(h)
                && expr->data.function.arg_count == 1
                && var_in(expr->data.function.args[0], var)) {
                return true;
            }
        }
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            if (solvetrig_has_trig(expr->data.function.args[i], var))
                return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Exponent extraction.                                              *
 *                                                                     *
 *  After TrigToExp + Together + Numerator, the residual is a Plus    *
 *  of Times-terms in which every var-bearing factor is               *
 *  Power[E, Times[Complex[0, n], var]] for some integer n (possibly  *
 *  with extra free-of-var Times factors).  This helper extracts the  *
 *  integer n from such a Power[E, e] node.  Returns false if e is   *
 *  not of the expected `n I var` shape.                              *
 * ------------------------------------------------------------------ */

static bool extract_exp_n(const Expr* e_expr, const Expr* var,
                          int64_t* n_out) {
    /* Accept either:
     *   Times[Complex[0, n], var]
     *   Times[Complex[0, n], var, free, free, ...]    -- no
     *   Times[I, var]   (less canonical; Mathilda usually canonicalises
     *                    I to Complex[0, 1], but be robust)
     *   Complex[0, n] * var: matches above
     *   The bare e == var (n=1, with implicit I missing) -- never, the
     *   exponent always carries the I factor after TrigToExp. */
    if (!e_expr || !var) return false;
    if (!head_is_sym(e_expr, SYM_Times)) return false;
    if (e_expr->data.function.arg_count < 2) return false;

    /* Walk the Times args.  We expect exactly one var occurrence and
     * exactly one Complex[0, n] / Times[I, ...] coefficient; everything
     * else must be free-of-var literals (but in canonical form there
     * shouldn't be any extras). */
    bool seen_var = false;
    int64_t n_acc = 0;
    bool have_n = false;
    for (size_t i = 0; i < e_expr->data.function.arg_count; i++) {
        const Expr* a = e_expr->data.function.args[i];
        if (a->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
            && a->data.symbol.name == var->data.symbol.name) {
            if (seen_var) return false;
            seen_var = true;
            continue;
        }
        /* Complex[0, n] with n an integer. */
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol.name == SYM_Complex
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_INTEGER
            && a->data.function.args[0]->data.integer == 0
            && a->data.function.args[1]->type == EXPR_INTEGER) {
            if (have_n) return false;
            n_acc = a->data.function.args[1]->data.integer;
            have_n = true;
            continue;
        }
        /* I as a bare symbol -- shouldn't appear post-canonicalisation
         * but handle defensively. */
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_I) {
            if (have_n) return false;
            n_acc = 1;
            have_n = true;
            continue;
        }
        return false;
    }
    if (!seen_var || !have_n) return false;
    *n_out = n_acc;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Substitute Power[E, m*I*var] -> Power[var, m]                     *
 *  (using `var` as the temporary u-symbol).  Updates *min_m / *max_m *
 *  with the range of exponents seen.  Returns NULL if the input      *
 *  contains a var occurrence outside an admissible Power[E, ...].    *
 * ------------------------------------------------------------------ */

typedef struct {
    int64_t min_m;
    int64_t max_m;
    bool    any_seen;
    bool    bad;       /* var occurred outside an admissible Power[E, _] */
} SubstState;

static Expr* substitute_walk(const Expr* e, const Expr* var,
                             SubstState* st) {
    if (!e) return NULL;
    if (st->bad) return NULL;

    /* Atomic var -- bad: every var occurrence must be inside an
     * admissible Power[E, m*I*var]. */
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) {
        st->bad = true;
        return NULL;
    }

    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Recognise Power[E, exp]. */
    if (head_is_sym(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_ = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E
            && var_in(exp_, var)) {
            int64_t m;
            if (!extract_exp_n(exp_, var, &m)) {
                st->bad = true;
                return NULL;
            }
            if (!st->any_seen) {
                st->min_m = m;
                st->max_m = m;
                st->any_seen = true;
            } else {
                if (m < st->min_m) st->min_m = m;
                if (m > st->max_m) st->max_m = m;
            }
            return mk_pow(expr_copy((Expr*)var), mk_int(m));
        }
        /* Power with var-bearing base/exponent that isn't E^(m I var):
         * fall through to the structural recursion. */
    }

    /* Generic recursion over args. */
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = substitute_walk(e->data.function.args[i], var, st);
        if (st->bad) {
            for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
            free(new_args);
            return NULL;
        }
    }
    Expr* head = expr_copy(e->data.function.head);
    Expr* out = expr_new_function(head, new_args, n);
    free(new_args);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Unwind a u-root into the var-branch via the Exp[I var] == u peel. *
 * ------------------------------------------------------------------ */

static Expr* unwind_u_root(Expr* u_val, Expr* var, Expr* dom,
                           const SolveInvOpts* opts) {
    Expr* exp_eq = mk_fn2("Equal",
        mk_pow(mk_sym("E"),
               mk_fn2("Times", mk_sym("I"), expr_copy(var))),
        expr_copy(u_val));
    exp_eq = eval_and_free(exp_eq);
    Expr* sols = solveinv_solve_inverse_equality(exp_eq, var, dom, opts);
    expr_free(exp_eq);
    return sols;
}

static Expr* concat_solutions(Expr* a, Expr* b) {
    if (!a) return b;
    if (!b) return a;
    if (!head_is_sym(a, SYM_List)) { expr_free(b); return a; }
    if (!head_is_sym(b, SYM_List)) { expr_free(b); return a; }
    size_t na = a->data.function.arg_count;
    size_t nb = b->data.function.arg_count;
    size_t n  = na + nb;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < na; i++)
        args[i] = expr_copy(a->data.function.args[i]);
    for (size_t i = 0; i < nb; i++)
        args[na + i] = expr_copy(b->data.function.args[i]);
    Expr* out = expr_new_function(mk_sym("List"), args, n);
    free(args);
    expr_free(a); expr_free(b);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                     *
 * ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ *
 *  Output cleanup: TrigToExp leaves each unwound branch as a complex  *
 *  logarithm (e.g. -I Log[(-1/2 - I/2) Sqrt[2]]); Simplify collapses  *
 *  those to the canonical Pi-multiple form (Pi/4 + Pi C[k]).          *
 * ------------------------------------------------------------------ */

static size_t sl_leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    size_t n = sl_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n += sl_leaf_count(e->data.function.args[i]);
    return n;
}

/* Count occurrences of the interned head `head` anywhere in `e`. */
static size_t sl_count_head(const Expr* e, const char* head) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    size_t n = (e->data.function.head->type == EXPR_SYMBOL
                && e->data.function.head->data.symbol.name == head) ? 1 : 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n += sl_count_head(e->data.function.args[i], head);
    return n;
}

/* Simplify `value` ONLY to strip the complex logarithms TrigToExp
 * leaves behind (e.g. -I Log[(-1/2 - I/2) Sqrt[2]] -> Pi/4), returning
 * the simplified form only when it removes a Log head without growing.
 * A value with no Log is already in principal+period form and is left
 * untouched, so readable sums are not churned into equivalent factored
 * forms.  Takes ownership of `value`; returns a value the caller owns. */
static Expr* simplify_smaller(Expr* value) {
    size_t logs = sl_count_head(value, SYM_Log);
    if (logs == 0) return value;
    Expr* simp = eval_and_free(mk_fn1("Simplify", expr_copy(value)));
    if (!simp) return value;
    if (sl_count_head(simp, SYM_Log) < logs
        && sl_leaf_count(simp) <= sl_leaf_count(value)) {
        expr_free(value);
        return simp;
    }
    expr_free(simp);
    return value;
}

/* In-place: Simplify the RHS value of every Rule in a solution
 * List[Rule[var, rhs]], descending into ConditionalExpression to
 * simplify its value while preserving the Element[C[k],Integers]
 * (or branch-cut) condition. */
static void simplify_solution_inplace(Expr* sol) {
    if (!head_is_sym(sol, SYM_List)) return;
    for (size_t i = 0; i < sol->data.function.arg_count; i++) {
        Expr* rule = sol->data.function.args[i];
        if (!head_is_sym(rule, SYM_Rule)
            || rule->data.function.arg_count != 2) continue;
        Expr* rhs = rule->data.function.args[1];
        if (head_is_sym(rhs, SYM_ConditionalExpression)
            && rhs->data.function.arg_count == 2) {
            Expr* val = rhs->data.function.args[0];
            rhs->data.function.args[0] = simplify_smaller(val);
        } else {
            rule->data.function.args[1] = simplify_smaller(rhs);
        }
    }
}

Expr* solvetrig_solve_trig_equality(Expr* equation, Expr* var,
                                    Expr* dom,
                                    const SolveInvOpts* opts) {
    if (!equation || !var) return NULL;
    if (!head_is_sym(equation, SYM_Equal)
        || equation->data.function.arg_count != 2) return NULL;

    if (!solvetrig_has_trig(equation, var)) return NULL;

    /* residual = lhs - rhs. */
    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];
    Expr* residual = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs), mk_neg(expr_copy(rhs))));

    /* 1. TrigToExp: every trig head over var becomes Exp[m*I*var]. */
    Expr* exp_form = eval_and_free(mk_fn1("TrigToExp", residual));

    /* 2. Together: collapse Tan/Cot/Sec/Csc-introduced denominators
     *    into a single Numerator/Denominator pair. */
    Expr* together = eval_and_free(
        internal_together((Expr*[]){ exp_form }, 1));

    /* 3. Take numerator. */
    Expr* numer = eval_and_free(
        internal_numerator((Expr*[]){ expr_copy(together) }, 1));
    expr_free(together);

    /* 4. Substitute every Power[E, m*I*var] with Power[var, m]. */
    SubstState st = { 0, 0, false, false };
    Expr* in_u = substitute_walk(numer, var, &st);
    expr_free(numer);
    if (!in_u || st.bad || !st.any_seen) {
        if (in_u) expr_free(in_u);
        return NULL;
    }

    /* 5. Multiply by var^(-min_m) to clear negative powers. */
    if (st.min_m < 0) {
        Expr* shift = mk_pow(expr_copy(var), mk_int(-st.min_m));
        in_u = eval_and_free(mk_fn2("Times", in_u, shift));
    }
    /* Expand so the polynomial detector sees an unfactored Plus. */
    in_u = eval_and_free(
        internal_expand((Expr*[]){ in_u }, 1));

    /* If the substitution collapsed to a constant (no var left), the
     * equation reduces to "constant == 0".  Solve trivially: nonzero
     * constant -> empty list, zero constant -> tautology. */
    if (!var_in(in_u, var)) {
        bool zero = (in_u->type == EXPR_INTEGER
                     && in_u->data.integer == 0);
        expr_free(in_u);
        if (zero) {
            Expr* empty = expr_new_function(mk_sym("List"), NULL, 0);
            return expr_new_function(mk_sym("List"),
                                     (Expr*[]){ empty }, 1);
        }
        return expr_new_function(mk_sym("List"), NULL, 0);
    }

    /* 6. Solve poly_in_u == 0 for var (treated as `u`). */
    Expr* u_eq = eval_and_free(
        mk_fn2("Equal", in_u, mk_int(0)));
    SolvePolyOpts polyopts = { false, false };
    Expr* u_solutions = solvepoly_solve_polynomial_equality(
        u_eq, var, NULL /* Complexes */, &polyopts);
    expr_free(u_eq);
    if (!u_solutions) return NULL;
    if (!head_is_sym(u_solutions, SYM_List)) {
        expr_free(u_solutions);
        return NULL;
    }

    /* 7. Unwind each u-root via Exp[I var] == u_i. */
    Expr* aggregate = expr_new_function(mk_sym("List"), NULL, 0);
    for (size_t i = 0; i < u_solutions->data.function.arg_count; i++) {
        Expr* sol = u_solutions->data.function.args[i];
        if (!head_is_sym(sol, SYM_List)
            || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (!head_is_sym(rule, SYM_Rule)
            || rule->data.function.arg_count != 2) continue;
        Expr* u_val = rule->data.function.args[1];
        Expr* branch = unwind_u_root(u_val, var, dom, opts);
        if (!branch) continue;
        aggregate = concat_solutions(aggregate, branch);
    }
    expr_free(u_solutions);

    /* 8. Canonicalise each branch: TrigToExp unwinding leaves complex
     *    logarithms that Simplify collapses to Pi-multiple form. */
    if (head_is_sym(aggregate, SYM_List)) {
        for (size_t i = 0; i < aggregate->data.function.arg_count; i++)
            simplify_solution_inplace(aggregate->data.function.args[i]);
    }
    return aggregate;
}

/* ================================================================== *
 *  Polynomial in a single transcendental kernel.                     *
 *                                                                    *
 *  Generalises the trig-to-exp path to any equation that is a        *
 *  polynomial in one invertible kernel g(x): substitute u = g(x),    *
 *  solve poly(u) == 0 via the polynomial specialist, then unwind     *
 *  each root u0 via g(x) == u0 through the inverse-function          *
 *  specialist.  Two kernel shapes are handled:                       *
 *                                                                    *
 *    exponential   E^(c x)  ->  u^c     (E^(2x)-3E^x+2 -> u^2-3u+2)  *
 *    generic head  H[x]^k   ->  u^k     (Log[x]^2-3Log[x]+2)         *
 *                                                                    *
 *  Trig kernels are already consumed upstream by the trig-to-exp     *
 *  pass; radicals by the radicals specialist.                        *
 * ------------------------------------------------------------------ */

/* Fresh temporary symbol standing for the kernel u = g(x). */
#define KERNEL_USYM "Solve`kernelU"

/* Extract the (real) integer exponent multiplier n from a Power[E, e]
 * exponent `e` of the shape `var` (n=1) or `Times[int..., var]`.  Any
 * non-integer, non-var factor makes this fail. */
static bool extract_real_exp_n(const Expr* e, const Expr* var, int64_t* n) {
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) {
        *n = 1;
        return true;
    }
    if (!head_is_sym(e, SYM_Times)) return false;
    bool seen_var = false, have = false;
    int64_t acc = 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* a = e->data.function.args[i];
        if (a->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
            && a->data.symbol.name == var->data.symbol.name) {
            if (seen_var) return false;
            seen_var = true;
        } else if (a->type == EXPR_INTEGER) {
            acc *= a->data.integer;
            have = true;
        } else {
            return false;   /* non-integer, non-var factor */
        }
    }
    if (!seen_var) return false;
    *n = have ? acc : 1;
    return true;
}

typedef struct {
    int64_t min_m, max_m;
    bool any_seen, bad;
} SubstReal;

/* Replace every Power[E, c*var] by Power[usym, c]; any bare var (or var
 * outside an admissible exponential) sets st->bad. */
static Expr* subst_exp_real_walk(const Expr* e, const Expr* var,
                                 const Expr* usym, SubstReal* st) {
    if (!e || st->bad) return NULL;
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) {
        st->bad = true;
        return NULL;
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (head_is_sym(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_ = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E
            && var_in(exp_, var)) {
            int64_t m;
            if (!extract_real_exp_n(exp_, var, &m)) { st->bad = true; return NULL; }
            if (!st->any_seen) { st->min_m = st->max_m = m; st->any_seen = true; }
            else { if (m < st->min_m) st->min_m = m; if (m > st->max_m) st->max_m = m; }
            return mk_pow(expr_copy((Expr*)usym), mk_int(m));
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = subst_exp_real_walk(e->data.function.args[i], var, usym, st);
        if (st->bad) { for (size_t j = 0; j < i; j++) expr_free(args[j]); free(args); return NULL; }
    }
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return out;
}

/* Solve poly(usym) == 0 and unwind each root u0 through `kernel == u0`
 * (kernel is Power[E,var] for the exponential path, or H[x] for the
 * generic path).  Returns the concatenated solution List, or NULL if any
 * root fails to unwind (so the caller can fall through). */
static Expr* solve_u_and_unwind(Expr* poly_in_u, Expr* usym, Expr* kernel,
                                Expr* var, Expr* dom, const SolveInvOpts* opts) {
    Expr* u_eq = eval_and_free(mk_fn2("Equal", expr_copy(poly_in_u), mk_int(0)));
    SolvePolyOpts po = { false, false };
    Expr* u_sols = solvepoly_solve_polynomial_equality(u_eq, usym, NULL, &po);
    expr_free(u_eq);
    if (!u_sols || !head_is_sym(u_sols, SYM_List)
        || u_sols->data.function.arg_count == 0) {
        if (u_sols) expr_free(u_sols);
        return NULL;
    }
    Expr* agg = expr_new_function(mk_sym("List"), NULL, 0);
    bool ok = true;
    for (size_t i = 0; i < u_sols->data.function.arg_count; i++) {
        Expr* sol = u_sols->data.function.args[i];
        if (!head_is_sym(sol, SYM_List) || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (!head_is_sym(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
        Expr* u_val = rule->data.function.args[1];
        Expr* eq = eval_and_free(mk_fn2("Equal", expr_copy(kernel), expr_copy(u_val)));
        Expr* branch = solveinv_solve_inverse_equality(eq, var, dom, opts);
        expr_free(eq);
        if (!branch) { ok = false; break; }
        agg = concat_solutions(agg, branch);
    }
    expr_free(u_sols);
    if (!ok) { expr_free(agg); return NULL; }
    return agg;
}

/* Exponential kernel: E^(c x).  Returns NULL if the residual is not a
 * polynomial in E^x. */
static Expr* try_exp_kernel(Expr* residual, Expr* var, Expr* usym,
                            Expr* dom, const SolveInvOpts* opts) {
    SubstReal st = { 0, 0, false, false };
    Expr* in_u = subst_exp_real_walk(residual, var, usym, &st);
    if (!in_u || st.bad || !st.any_seen) { if (in_u) expr_free(in_u); return NULL; }
    if (st.min_m < 0) {
        in_u = eval_and_free(mk_fn2("Times", in_u,
            mk_pow(expr_copy(usym), mk_int(-st.min_m))));
    }
    in_u = eval_and_free(internal_expand((Expr*[]){ in_u }, 1));
    if (var_in(in_u, var) || !var_in(in_u, usym)) { expr_free(in_u); return NULL; }
    Expr* kernel = mk_pow(mk_sym("E"), expr_copy(var));   /* E^var */
    Expr* out = solve_u_and_unwind(in_u, usym, kernel, var, dom, opts);
    expr_free(kernel);
    expr_free(in_u);
    return out;
}

/* Is `h` an invertible single-argument head whose inverse the inverse
 * specialist can peel? */
static bool is_kernel_head(const char* h) {
    return h == SYM_Log
        || h == SYM_Sin  || h == SYM_Cos  || h == SYM_Tan
        || h == SYM_Cot  || h == SYM_Sec  || h == SYM_Csc
        || h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh
        || h == SYM_Coth || h == SYM_Sech || h == SYM_Csch
        || h == SYM_ArcSin  || h == SYM_ArcCos  || h == SYM_ArcTan
        || h == SYM_ArcCot  || h == SYM_ArcSec  || h == SYM_ArcCsc
        || h == SYM_ArcSinh || h == SYM_ArcCosh || h == SYM_ArcTanh
        || h == SYM_ArcCoth || h == SYM_ArcSech || h == SYM_ArcCsch;
}

/* Collect the distinct invertible-head subexpressions of `e` that
 * contain `var`, into `ks` (deduplicated, up to `max`). */
static void collect_kernels(const Expr* e, const Expr* var,
                            Expr** ks, size_t* n, size_t max) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && is_kernel_head(e->data.function.head->data.symbol.name)
        && var_in(e, var)) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++)
            if (expr_eq(ks[i], (Expr*)e)) { dup = true; break; }
        if (!dup && *n < max) ks[(*n)++] = expr_copy((Expr*)e);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_kernels(e->data.function.args[i], var, ks, n, max);
}

/* Generic-head kernel: substitute one invertible-head subexpression
 * g(x) -> u and require the residual to become a polynomial in u free of
 * var.  Returns NULL if no single kernel works. */
static Expr* try_generic_kernel(Expr* residual, Expr* var, Expr* usym,
                                Expr* dom, const SolveInvOpts* opts) {
    Expr* ks[8];
    size_t nk = 0;
    collect_kernels(residual, var, ks, &nk, 8);
    Expr* result = NULL;
    for (size_t i = 0; i < nk && !result; i++) {
        Expr* rule = mk_fn2("Rule", expr_copy(ks[i]), expr_copy(usym));
        Expr* subbed = eval_and_free(internal_replace_all(
            (Expr*[]){ expr_copy(residual), rule }, 2));
        if (var_in(subbed, var) || !var_in(subbed, usym)
            || !is_polynomial(subbed, &usym, 1)) {
            expr_free(subbed);
            continue;
        }
        result = solve_u_and_unwind(subbed, usym, ks[i], var, dom, opts);
        expr_free(subbed);
    }
    for (size_t i = 0; i < nk; i++) expr_free(ks[i]);
    return result;
}

Expr* solvetrig_solve_poly_in_kernel(Expr* equation, Expr* var,
                                     Expr* dom, const SolveInvOpts* opts) {
    if (!equation || !var) return NULL;
    if (!head_is_sym(equation, SYM_Equal)
        || equation->data.function.arg_count != 2) return NULL;

    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];
    Expr* residual = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs), mk_neg(expr_copy(rhs))));
    Expr* usym = mk_sym(KERNEL_USYM);

    Expr* out = try_exp_kernel(residual, var, usym, dom, opts);
    if (!out) out = try_generic_kernel(residual, var, usym, dom, opts);

    /* Canonicalise unwound branches (strip complex logs). */
    if (out && head_is_sym(out, SYM_List)) {
        for (size_t i = 0; i < out->data.function.arg_count; i++)
            simplify_solution_inplace(out->data.function.args[i]);
    }

    expr_free(usym);
    expr_free(residual);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry: Solve`SolveTrigEquality                  *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_trig_equality(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* dom      = (argc >= 3) ? res->data.function.args[2] : NULL;
    SolveInvOpts opts = { true, intern_symbol("C") };
    return solvetrig_solve_trig_equality(equation, var, dom, &opts);
}

void solvetrig_init(void) {
    symtab_add_builtin("Solve`SolveTrigEquality",
                       builtin_solve_trig_equality);
    SymbolDef* def = symtab_get_def("Solve`SolveTrigEquality");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveTrigEquality",
        "Solve`SolveTrigEquality[lhs == rhs, var]\n"
        "Solve`SolveTrigEquality[lhs == rhs, var, dom]\n"
        "\tThe trig canonicalisation pre-pass used by Solve.  Converts\n"
        "\tevery trig head over the variable to Power[E, m I var] via\n"
        "\tTrigToExp, substitutes u = Exp[I var] to obtain a polynomial\n"
        "\tin u, solves it, then unwinds each u-root through the\n"
        "\tinverse-function specialist to produce the periodic\n"
        "\tvar-branch family.  Mirrors Maxima's trig-cannon /\n"
        "\ttrig-subst-p path but routes through Exp instead of sin/cos.");
}
