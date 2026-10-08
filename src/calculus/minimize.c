/* minimize.c — Minimize / Maximize: exact (symbolic) global optimization.
 *
 * The exact counterpart to NMinimize/NMaximize (nm_driver.c). Where those run a
 * stochastic numeric global search, Minimize finds a GLOBAL optimum exactly and
 * returns {f_opt, {x -> x_opt, ...}} with exact (Integer / Rational / radical /
 * Root[]) values, delegating to NMinimize only when the input is inexact.
 *
 * Soundness invariant (every path honours it): a result is emitted only when its
 * value, attainment and feasibility are each backed by an EXACT oracle —
 *   - Solve[..., Reals] / the zero-dimensional engine for the candidate set,
 *   - the FLINT qqbar sign oracle (rru_sign_of / rru_sign_compare) for ordering,
 *   - Reduce[..., Reals] for the global lower-bound certificate and emptiness.
 * Any "don't know" from those oracles (rru_sign == -2, an unclean Solve/Reduce
 * result) is a hard stop: the path is abandoned and the expression is left
 * unevaluated (NULL) rather than guessed. The ONLY licence to return an inexact
 * answer is inexact INPUT, which routes to NMinimize from the top.
 *
 * Scope this release (see docs/spec/builtins/numerical_calculus.md):
 *   (a) univariate polynomial objective (unbounded/attainment via the degree
 *       parity + leading-coefficient tail theorem);
 *   (b) multivariate polynomial objective with isolated critical points;
 *   (c) polynomial objective under polynomial constraints over the Reals
 *       (KKT/active-set enumeration + a Reduce lower-bound certificate; small
 *       linear programs fall out as the degenerate case).
 * Deferred (declines, never guesses): transcendental closed forms, parametric
 * Piecewise answers, positive-dimensional minimizer sets, general unbounded/
 * not-attained via quantifier elimination, exact Integers/ILP optimisation.
 *
 * Minimize/Maximize are Protected but NOT HoldAll (matching Mathematica): the
 * variables are unbound symbols that evaluate to themselves, so the objective
 * and constraints arrive in symbolic form and the engine proceeds by
 * value-preserving substitution (D, Solve, Reduce). */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "message.h"            /* mth_message + mth_msg_suppress_push/pop */
#include "reduce_real_util.h"   /* rru_is_polynomial / rru_sign_of / rru_sign_compare */

/* ------------------------------------------------------------------ *
 *  Small expression builders (the SYM_* names are interned)           *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}

static bool mz_is_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}
static bool mz_is_true(const Expr* e) { return mz_is_sym(e, SYM_True); }

static const char* mz_noun(const char* head) {
    return (strcmp(head, "Maximize") == 0) ? "maximum" : "minimum";
}

/* An option argument: Rule[...] or RuleDelayed[...] (Method -> ..., etc.). */
static bool mz_is_option_arg(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}

/* A relational / logical tree — distinguishes {f, cons} from a vector
 * objective (same head set as NMinimize's nm_is_constraint_tree). */
static bool mz_is_constraint_tree(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_And || h == SYM_Or || h == SYM_Not || h == SYM_Xor ||
           h == SYM_Implies || h == SYM_Equal || h == SYM_Unequal ||
           h == SYM_Less || h == SYM_LessEqual || h == SYM_Greater ||
           h == SYM_GreaterEqual || h == SYM_Inequality || h == SYM_Element;
}

/* A fresh, collision-free symbol (for Lagrange multipliers in the KKT system). */
static Expr* mz_fresh_symbol(void) {
    static uint64_t ctr = 0;
    char buf[64];
    for (;;) {
        snprintf(buf, sizeof(buf), "Minimize$%llu", (unsigned long long)ctr++);
        if (!symtab_lookup(buf)) break;
    }
    return expr_new_symbol(buf);
}

/* ±Infinity as canonical Expr values. */
static Expr* mz_pos_inf(void) { return expr_new_symbol(SYM_Infinity); }
static Expr* mz_neg_inf(void) {
    return eval_and_free(fn2(SYM_Times, mk_int(-1), expr_new_symbol(SYM_Infinity)));
}

/* a - b, evaluated. */
static Expr* mz_sub(const Expr* a, const Expr* b) {
    return eval_and_free(fn2(SYM_Plus, expr_copy((Expr*)a),
                             fn2(SYM_Times, mk_int(-1), expr_copy((Expr*)b))));
}

/* D[e, x], evaluated. */
static Expr* mz_deriv(const Expr* e, const Expr* x) {
    return eval_and_free(fn2(SYM_D, expr_copy((Expr*)e), expr_copy((Expr*)x)));
}

/* List[v1, ..., vn] of copies. */
static Expr* mz_varlist(Expr* const* vars, size_t n) {
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(vars[i]);
    Expr* L = expr_new_function(expr_new_symbol(SYM_List), a, n);
    free(a);
    return L;
}

/* True iff `e` mentions any of the given symbols. */
static bool mz_contains_var(const Expr* e, Expr* const* vars, size_t n) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; i++)
            if (e->data.symbol.name == vars[i]->data.symbol.name) return true;
        return false;
    }
    if (e->type == EXPR_FUNCTION) {
        if (mz_contains_var(e->data.function.head, vars, n)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (mz_contains_var(e->data.function.args[i], vars, n)) return true;
    }
    return false;
}

/* True iff `e` carries an inexact (machine or MPFR) real anywhere. */
static bool mz_has_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION) {
        if (mz_has_inexact(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (mz_has_inexact(e->data.function.args[i])) return true;
    }
    return false;
}

/* PolynomialQ[e, {vars}] -> True. */
static bool mz_poly_in_vars(const Expr* e, Expr* const* vars, size_t n) {
    if (n == 1) return rru_is_polynomial(e, vars[0]);
    Expr* q = eval_and_free(fn2(SYM_PolynomialQ, expr_copy((Expr*)e),
                               mz_varlist(vars, n)));
    bool r = mz_is_true(q);
    expr_free(q);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Substitution, result builders                                      *
 * ------------------------------------------------------------------ */

/* ReplaceAll[e, {var_i -> val_i}], evaluated. */
static Expr* mz_subst_eval(const Expr* e, Expr* const* vars,
                           Expr* const* vals, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        rules[i] = fn2(SYM_Rule, expr_copy(vars[i]), expr_copy(vals[i]));
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    mth_msg_suppress_push();                /* internal probe: quiet numeric-domain noise */
    Expr* r = eval_and_free(fn2(SYM_ReplaceAll, expr_copy((Expr*)e), rl));
    mth_msg_suppress_pop();
    return r;
}

/* {value, {var_i -> val_i}} — consumes `value`, copies vars/vals. */
static Expr* mz_result(Expr* value, Expr* const* vars, Expr* const* vals, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        rules[i] = fn2(SYM_Rule, expr_copy(vars[i]), expr_copy(vals[i]));
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    return fn2(SYM_List, value, rl);
}

/* {value, {var_i -> Indeterminate}} — consumes `value`. */
static Expr* mz_result_indet(Expr* value, Expr* const* vars, size_t n) {
    Expr** vals = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) vals[i] = expr_new_symbol(SYM_Indeterminate);
    Expr* r = mz_result(value, vars, vals, n);
    for (size_t i = 0; i < n; i++) expr_free(vals[i]);
    free(vals);
    return r;
}

/* Negate the value field of a {value, rules} result — consumes and rebuilds. */
static Expr* mz_negate_value(Expr* result) {
    if (!result || result->type != EXPR_FUNCTION ||
        !mz_is_sym(result->data.function.head, SYM_List) ||
        result->data.function.arg_count != 2)
        return result;
    Expr* v = expr_copy(result->data.function.args[0]);
    Expr* rules = expr_copy(result->data.function.args[1]);
    expr_free(result);
    Expr* nv = eval_and_free(fn2(SYM_Times, mk_int(-1), v));
    return fn2(SYM_List, nv, rules);
}

/* ------------------------------------------------------------------ *
 *  Candidate-point extraction and selection                          *
 * ------------------------------------------------------------------ */

/* Free a points table: npts rows, each an Expr*[n]. */
static void mz_free_points(Expr*** pts, size_t npts, size_t n) {
    if (!pts) return;
    for (size_t i = 0; i < npts; i++) {
        if (!pts[i]) continue;
        for (size_t j = 0; j < n; j++) expr_free(pts[i][j]);
        free(pts[i]);
    }
    free(pts);
}

/* Parse a Solve result L = List[ List[Rule[v,val] ...] ... ] into a points
 * table, extracting the value of each of `xvars` (order preserved) and
 * requiring every extracted value to be free of all `forbid` symbols (so a
 * parametric / underdetermined branch is rejected). Returns 1 on success
 * (fills pts_out and npts_out; the count may be 0 for a clean empty set), or 0
 * to decline (an unclean shape: unevaluated Solve, ConditionalExpression, a
 * missing or parametric coordinate). */
static int mz_parse_points(const Expr* L, Expr* const* xvars, size_t nx,
                           Expr* const* forbid, size_t nforbid,
                           Expr**** pts_out, size_t* npts_out) {
    if (!L || L->type != EXPR_FUNCTION || !mz_is_sym(L->data.function.head, SYM_List))
        return 0;
    size_t m = L->data.function.arg_count;
    Expr*** pts = m ? (Expr***)malloc(sizeof(Expr**) * m) : NULL;
    size_t cnt = 0;
    for (size_t i = 0; i < m; i++) {
        const Expr* sol = L->data.function.args[i];
        if (!sol || sol->type != EXPR_FUNCTION ||
            !mz_is_sym(sol->data.function.head, SYM_List))
            goto bad;
        Expr** row = (Expr**)malloc(sizeof(Expr*) * (nx ? nx : 1));
        for (size_t j = 0; j < nx; j++) row[j] = NULL;
        for (size_t j = 0; j < nx; j++) {
            const Expr* val = NULL;
            for (size_t r = 0; r < sol->data.function.arg_count; r++) {
                const Expr* rule = sol->data.function.args[r];
                if (rule->type == EXPR_FUNCTION &&
                    mz_is_sym(rule->data.function.head, SYM_Rule) &&
                    rule->data.function.arg_count == 2 &&
                    rule->data.function.args[0]->type == EXPR_SYMBOL &&
                    rule->data.function.args[0]->data.symbol.name ==
                        xvars[j]->data.symbol.name) {
                    val = rule->data.function.args[1];
                    break;
                }
            }
            if (!val || mz_contains_var(val, forbid, nforbid)) {
                for (size_t k = 0; k < nx; k++) if (row[k]) expr_free(row[k]);
                free(row);
                goto bad;
            }
            row[j] = expr_copy((Expr*)val);
        }
        pts[cnt++] = row;
    }
    *pts_out = pts;
    *npts_out = cnt;
    return 1;
bad:
    mz_free_points(pts, cnt, nx);
    return 0;
}

/* Solve[system, {vars}, Reals], evaluated — consumes `system`. Internal probe:
 * an unsolvable / parametric subproblem is expected (and skipped), so Solve's
 * own diagnostics (Solve::nsdim, ...) are suppressed. */
static Expr* mz_solve_real(Expr* system, Expr* const* vars, size_t n) {
    mth_msg_suppress_push();
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ system, mz_varlist(vars, n), expr_new_symbol(SYM_Reals) }, 3));
    mth_msg_suppress_pop();
    return r;
}

/* Evaluate f at each of `nrows` points (each Expr*[n]) and select the least by
 * exact comparison. On success sets *best_val (owned) and *best_idx, returns 1;
 * returns 0 to decline (an undecidable comparison). Requires nrows >= 1. */
static int mz_pick_min(const Expr* f, Expr* const* vars, size_t n,
                       Expr*** rows, size_t nrows,
                       Expr** best_val_out, size_t* best_idx_out) {
    Expr* bestv = mz_subst_eval(f, vars, rows[0], n);
    size_t bi = 0;
    for (size_t i = 1; i < nrows; i++) {
        Expr* v = mz_subst_eval(f, vars, rows[i], n);
        int c = rru_sign_compare(v, bestv);   /* sign(v - bestv) */
        if (c == -2) { expr_free(v); expr_free(bestv); return 0; }
        if (c < 0) { expr_free(bestv); bestv = v; bi = i; }
        else expr_free(v);
    }
    *best_val_out = bestv;
    *best_idx_out = bi;
    return 1;
}

/* ------------------------------------------------------------------ *
 *  Reduce-based oracles                                               *
 * ------------------------------------------------------------------ */

/* 1 iff  A ⇒ P  is provable over the Reals (Reduce[A && !P] is False), else 0. */
static int mz_entails(const Expr* A, const Expr* P, Expr* const* vars, size_t n) {
    Expr* notP = eval_and_free(fn1(SYM_Not, expr_copy((Expr*)P)));
    Expr* stmt = (!A || mz_is_true(A)) ? notP
               : fn2(SYM_And, expr_copy((Expr*)A), notP);
    mth_msg_suppress_push();
    Expr* red = eval_and_free(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ stmt, mz_varlist(vars, n), expr_new_symbol(SYM_Reals) }, 3));
    mth_msg_suppress_pop();
    int r = mz_is_sym(red, SYM_False) ? 1 : 0;
    expr_free(red);
    return r;
}

/* True iff Reduce proves the constraint region empty over the Reals. */
static bool mz_region_empty(const Expr* cons, Expr* const* vars, size_t n) {
    mth_msg_suppress_push();
    Expr* red = eval_and_free(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ expr_copy((Expr*)cons), mz_varlist(vars, n),
                   expr_new_symbol(SYM_Reals) }, 3));
    mth_msg_suppress_pop();
    bool r = mz_is_sym(red, SYM_False);
    expr_free(red);
    return r;
}

/* Feasibility of a point in `cons`: 1 True, 0 False, -1 undecided. */
static int mz_feasible_at(const Expr* cons, Expr* const* vars,
                          Expr* const* vals, size_t n) {
    Expr* r = mz_subst_eval(cons, vars, vals, n);
    int res = mz_is_true(r) ? 1 : (mz_is_sym(r, SYM_False) ? 0 : -1);
    expr_free(r);
    return res;
}

/* Topological closure of a constraint tree: strict <, > become <=, >=. */
static Expr* mz_closure(const Expr* c) {
    if (c->type != EXPR_FUNCTION) return expr_copy((Expr*)c);
    const char* h = (c->data.function.head->type == EXPR_SYMBOL)
        ? c->data.function.head->data.symbol.name : NULL;
    const char* nh = h;
    if (h == SYM_Less) nh = SYM_LessEqual;
    else if (h == SYM_Greater) nh = SYM_GreaterEqual;
    size_t n = c->data.function.arg_count;
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = mz_closure(c->data.function.args[i]);
    Expr* head = nh ? expr_new_symbol(nh) : expr_copy(c->data.function.head);
    Expr* r = expr_new_function(head, a, n);
    free(a);
    return r;
}

/* ------------------------------------------------------------------ *
 *  (a) Univariate unconstrained polynomial                            *
 * ------------------------------------------------------------------ */

static Expr* mz_univar_poly(const Expr* f, Expr* x, const char* head) {
    Expr* vars1[1] = { x };
    if (!rru_is_polynomial(f, x)) return NULL;

    Expr* cl = eval_and_free(fn2(SYM_CoefficientList, expr_copy((Expr*)f),
                                 expr_copy(x)));
    if (!cl || cl->type != EXPR_FUNCTION ||
        !mz_is_sym(cl->data.function.head, SYM_List) ||
        cl->data.function.arg_count == 0) { expr_free(cl); return NULL; }
    size_t d1 = cl->data.function.arg_count;    /* degree + 1 */
    int d = (int)d1 - 1;

    if (d == 0) {                               /* constant objective */
        Expr* val = expr_copy(cl->data.function.args[0]);
        Expr* zero = mk_int(0);
        Expr* vals[1] = { zero };
        Expr* r = mz_result(val, vars1, vals, 1);
        expr_free(zero);
        expr_free(cl);
        return r;
    }

    int s = rru_sign_of(cl->data.function.args[d1 - 1]);   /* leading coeff sign */
    expr_free(cl);
    if (s == -2) return NULL;                   /* parametric leading coeff: decline */

    /* Tail theorem: odd degree, or even degree with a negative leading
     * coefficient, is unbounded below (above, for Maximize of -f). */
    if (d % 2 == 1 || s < 0) {
        mth_message(head, "natt",
                "The %s is not attained at any point satisfying the given "
                "constraints.", mz_noun(head));
        return mz_result_indet(mz_neg_inf(), vars1, 1);
    }

    /* Even degree, positive leading coefficient: coercive, so the global
     * minimum is attained at a real stationary point. */
    Expr* fp = mz_deriv(f, x);
    Expr* eqn = fn2(SYM_Equal, fp, mk_int(0));
    Expr* L = mz_solve_real(eqn, vars1, 1);
    Expr*** pts = NULL; size_t npts = 0;
    int ok = mz_parse_points(L, vars1, 1, vars1, 1, &pts, &npts);
    expr_free(L);
    if (!ok) return NULL;
    if (npts == 0) { mz_free_points(pts, npts, 1); return NULL; }

    Expr* bestv = NULL; size_t bi = 0;
    if (!mz_pick_min(f, vars1, 1, pts, npts, &bestv, &bi)) {
        mz_free_points(pts, npts, 1);
        return NULL;
    }
    Expr* pt[1] = { pts[bi][0] };
    Expr* r = mz_result(bestv, vars1, pt, 1);   /* bestv consumed; pt copied */
    mz_free_points(pts, npts, 1);
    return r;
}

/* ------------------------------------------------------------------ *
 *  (b)/(c) Multivariate / constrained polynomial over the Reals       *
 * ------------------------------------------------------------------ */

#define MZ_MAX_INEQ 6               /* 2^MZ_MAX_INEQ active-set subproblems */

typedef struct { Expr* z; bool strict; } MzIneq;   /* feasible <=> z <= 0 (< if strict) */

/* Append a borrowed conjunct pointer. */
static void mz_push_conj(const Expr*** arr, size_t* n, size_t* cap, const Expr* c) {
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 8;
        *arr = (const Expr**)realloc(*arr, sizeof(const Expr*) * *cap); }
    (*arr)[(*n)++] = c;
}
static void mz_collect_and(const Expr* c, const Expr*** arr, size_t* n, size_t* cap) {
    if (c->type == EXPR_FUNCTION && mz_is_sym(c->data.function.head, SYM_And))
        for (size_t i = 0; i < c->data.function.arg_count; i++)
            mz_collect_and(c->data.function.args[i], arr, n, cap);
    else
        mz_push_conj(arr, n, cap, c);
}

/* Owned deep copy of a constraint tree with chained Inequality[e0, op1, e1,
 * op2, e2, ...] rewritten to the And of its pairwise comparisons op_i[e_i,
 * e_{i+1}], so the per-conjunct classifier sees only two-argument relations. */
static Expr* mz_normalize_cons(const Expr* c) {
    if (c->type != EXPR_FUNCTION) return expr_copy((Expr*)c);
    const Expr* h = c->data.function.head;
    size_t m = c->data.function.arg_count;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Inequality &&
        m >= 3 && (m % 2) == 1) {
        size_t k = (m - 1) / 2;
        Expr** parts = (Expr**)malloc(sizeof(Expr*) * k);
        for (size_t i = 0; i < k; i++) {
            const Expr* op = c->data.function.args[2 * i + 1];
            Expr* lhs = mz_normalize_cons(c->data.function.args[2 * i]);
            Expr* rhs = mz_normalize_cons(c->data.function.args[2 * i + 2]);
            parts[i] = (op->type == EXPR_SYMBOL)
                ? fn2(op->data.symbol.name, lhs, rhs)
                : expr_new_function(expr_copy((Expr*)op), (Expr*[]){ lhs, rhs }, 2);
        }
        Expr* r = (k == 1) ? parts[0]
                           : expr_new_function(expr_new_symbol(SYM_And), parts, k);
        free(parts);
        return r;
    }
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) a[i] = mz_normalize_cons(c->data.function.args[i]);
    Expr* r = expr_new_function(expr_copy((Expr*)h), a, m);
    free(a);
    return r;
}

static Expr* mz_exact_poly(const Expr* f, const Expr* cons,
                           Expr* const* vars, size_t n, const char* head) {
    Expr* result = NULL;
    const Expr** conj = NULL; size_t nconj = 0, ccap = 0;
    Expr** EQ = NULL; size_t neq = 0;             /* equality zero-exprs (z == 0) */
    MzIneq* IN = NULL; size_t nin = 0;            /* inequalities z (<=) 0 */
    Expr** df = NULL;                             /* ∂f/∂x_j */
    Expr*** dzEQ = NULL; Expr*** dzIN = NULL;     /* constraint jacobian rows */
    Expr*** cand = NULL; size_t ncand = 0, candcap = 0;
    Expr*** feas = NULL; size_t nfeas = 0;
    Expr* CL = NULL;
    Expr* consN = NULL;                           /* cons with Inequality chains split */
    bool constrained = (cons != NULL) && !mz_is_true(cons);

    if (!mz_poly_in_vars(f, vars, n)) return NULL;

    /* Parse the constraint conjunction into equalities and inequalities. */
    if (constrained) {
        consN = mz_normalize_cons(cons);
        mz_collect_and(consN, &conj, &nconj, &ccap);
        EQ = (Expr**)malloc(sizeof(Expr*) * (nconj ? nconj : 1));
        IN = (MzIneq*)malloc(sizeof(MzIneq) * (nconj ? nconj : 1));
        for (size_t i = 0; i < nconj; i++) {
            const Expr* c = conj[i];
            if (c->type != EXPR_FUNCTION || c->data.function.arg_count != 2 ||
                c->data.function.head->type != EXPR_SYMBOL) goto cleanup;
            const char* h = c->data.function.head->data.symbol.name;
            const Expr* a = c->data.function.args[0];
            const Expr* b = c->data.function.args[1];
            Expr* z;
            if (h == SYM_Equal)        { z = mz_sub(a, b); EQ[neq++] = z; }
            else if (h == SYM_LessEqual || h == SYM_Less) {
                z = mz_sub(a, b); IN[nin].z = z; IN[nin].strict = (h == SYM_Less); nin++;
            } else if (h == SYM_GreaterEqual || h == SYM_Greater) {
                z = mz_sub(b, a); IN[nin].z = z; IN[nin].strict = (h == SYM_Greater); nin++;
            } else goto cleanup;      /* Inequality chain / Or / Element: decline */
        }
        if (nin > MZ_MAX_INEQ) goto cleanup;
        for (size_t e = 0; e < neq; e++) if (!mz_poly_in_vars(EQ[e], vars, n)) goto cleanup;
        for (size_t k = 0; k < nin; k++) if (!mz_poly_in_vars(IN[k].z, vars, n)) goto cleanup;
    }

    /* Objective and constraint gradients. */
    df = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t j = 0; j < n; j++) df[j] = mz_deriv(f, vars[j]);
    if (neq) { dzEQ = (Expr***)malloc(sizeof(Expr**) * neq);
        for (size_t e = 0; e < neq; e++) { dzEQ[e] = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < n; j++) dzEQ[e][j] = mz_deriv(EQ[e], vars[j]); } }
    if (nin) { dzIN = (Expr***)malloc(sizeof(Expr**) * nin);
        for (size_t k = 0; k < nin; k++) { dzIN[k] = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < n; j++) dzIN[k][j] = mz_deriv(IN[k].z, vars[j]); } }

    /* Enumerate active sets: all equalities are always active; each subset of
     * the inequalities is additionally activated (set to its boundary z == 0).
     * The KKT stationarity + boundary system for each active set is solved over
     * the Reals; a parametric / unsolvable subset is SKIPPED (not fatal) — the
     * Reduce lower-bound certificate below is the ultimate soundness gate. */
    size_t nsub = (size_t)1 << nin;
    for (size_t mask = 0; mask < nsub; mask++) {
        size_t nact = neq;
        for (size_t k = 0; k < nin; k++) if (mask & ((size_t)1 << k)) nact++;

        /* Active constraint z-exprs and jacobian rows: EQ first, then active IN. */
        Expr** actz = (Expr**)malloc(sizeof(Expr*) * (nact ? nact : 1));
        Expr*** actdz = (Expr***)malloc(sizeof(Expr**) * (nact ? nact : 1));
        size_t ai = 0;
        for (size_t e = 0; e < neq; e++) { actz[ai] = EQ[e]; actdz[ai] = dzEQ[e]; ai++; }
        for (size_t k = 0; k < nin; k++)
            if (mask & ((size_t)1 << k)) { actz[ai] = IN[k].z; actdz[ai] = dzIN[k]; ai++; }

        /* Variables: the x_j plus one fresh Lagrange multiplier per active c. */
        size_t nall = n + nact;
        Expr** allv = (Expr**)malloc(sizeof(Expr*) * nall);
        Expr** lam = (Expr**)malloc(sizeof(Expr*) * (nact ? nact : 1));
        for (size_t j = 0; j < n; j++) allv[j] = expr_copy(vars[j]);
        for (size_t c = 0; c < nact; c++) { lam[c] = mz_fresh_symbol();
            allv[n + c] = expr_copy(lam[c]); }

        /* System: ∂f/∂x_j - Σ_c λ_c ∂z_c/∂x_j == 0 (n eqns), then z_c == 0. */
        Expr** sys = (Expr**)malloc(sizeof(Expr*) * (n + nact));
        size_t si = 0;
        for (size_t j = 0; j < n; j++) {
            Expr* sum;
            if (nact == 0) sum = mk_int(0);
            else {
                Expr** terms = (Expr**)malloc(sizeof(Expr*) * nact);
                for (size_t c = 0; c < nact; c++)
                    terms[c] = fn2(SYM_Times, expr_copy(lam[c]),
                                   expr_copy(actdz[c][j]));
                sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, nact);
                free(terms);
            }
            Expr* lhs = fn2(SYM_Plus, expr_copy(df[j]),
                            fn2(SYM_Times, mk_int(-1), sum));
            sys[si++] = fn2(SYM_Equal, lhs, mk_int(0));
        }
        for (size_t c = 0; c < nact; c++)
            sys[si++] = fn2(SYM_Equal, expr_copy(actz[c]), mk_int(0));
        Expr* system = expr_new_function(expr_new_symbol(SYM_List), sys, si);
        free(sys);

        Expr* L = mz_solve_real(system, allv, nall);
        Expr*** pp = NULL; size_t npp = 0;
        int ok = mz_parse_points(L, vars, n, allv, nall, &pp, &npp);
        expr_free(L);
        for (size_t c = 0; c < nact; c++) expr_free(lam[c]);
        for (size_t j = 0; j < nall; j++) expr_free(allv[j]);
        free(lam); free(allv); free(actz); free(actdz);

        if (!ok) continue;                      /* parametric/unclean subset: skip */
        for (size_t i = 0; i < npp; i++) {
            if (ncand == candcap) { candcap = candcap ? candcap * 2 : 8;
                cand = (Expr***)realloc(cand, sizeof(Expr**) * candcap); }
            cand[ncand++] = pp[i];
        }
        free(pp);                               /* rows moved into cand */
    }

    /* Feasibility filter against the closure of the region. */
    CL = constrained ? mz_closure(consN) : expr_new_symbol(SYM_True);
    feas = (Expr***)malloc(sizeof(Expr**) * (ncand ? ncand : 1));
    for (size_t i = 0; i < ncand; i++)
        if (mz_feasible_at(CL, vars, cand[i], n) == 1) feas[nfeas++] = cand[i];

    if (nfeas == 0) {
        if (constrained && mz_region_empty(consN, vars, n)) {
            mth_message(head, "infeas",
                    "The constraints are infeasible; the feasible region is empty.");
            result = mz_result_indet(mz_pos_inf(), vars, n);
        }
        goto cleanup;                           /* else: unbounded/undecided → decline */
    }

    {
        Expr* bestv = NULL; size_t bi = 0;
        if (!mz_pick_min(f, vars, n, feas, nfeas, &bestv, &bi)) goto cleanup;

        /* Global lower-bound certificate: closure(cons) ⇒ f >= bestv. */
        Expr* ge = fn2(SYM_GreaterEqual, expr_copy((Expr*)f), expr_copy(bestv));
        int proven = mz_entails(CL, ge, vars, n);
        expr_free(ge);
        if (!proven) { expr_free(bestv); goto cleanup; }

        /* Attainment: the minimising point must satisfy the ORIGINAL (strict)
         * constraints. A point feasible only in the closure (excluded strict
         * boundary) is an un-attained infimum — declined this release. */
        int attained = constrained ? mz_feasible_at(consN, vars, feas[bi], n) : 1;
        if (attained != 1) { expr_free(bestv); goto cleanup; }

        Expr** pt = feas[bi];                   /* borrowed row of n values */
        result = mz_result(bestv, vars, pt, n); /* bestv consumed; pt copied */
    }

cleanup:
    free(conj);
    for (size_t e = 0; e < neq; e++) expr_free(EQ[e]);
    free(EQ);
    for (size_t k = 0; k < nin; k++) expr_free(IN[k].z);
    free(IN);
    if (df) {
        for (size_t j = 0; j < n; j++) expr_free(df[j]);
        free(df);
    }
    if (dzEQ) {
        for (size_t e = 0; e < neq; e++) {
            for (size_t j = 0; j < n; j++) expr_free(dzEQ[e][j]);
            free(dzEQ[e]);
        }
        free(dzEQ);
    }
    if (dzIN) {
        for (size_t k = 0; k < nin; k++) {
            for (size_t j = 0; j < n; j++) expr_free(dzIN[k][j]);
            free(dzIN[k]);
        }
        free(dzIN);
    }
    mz_free_points(cand, ncand, n);
    free(feas);                                 /* rows owned by cand, freed above */
    expr_free(CL);
    expr_free(consN);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Numeric fallback: delegate inexact input to NMinimize/NMaximize    *
 * ------------------------------------------------------------------ */

static Expr* mz_numeric_fallback(Expr* arg0, Expr* arg1, const char* domname,
                                 Expr** opts, size_t nopts, bool is_max) {
    size_t k = 0;
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (2 + nopts));
    a[k++] = expr_copy(arg0);
    if (domname == SYM_Integers)              /* NMinimize wants Element[vars, dom] */
        a[k++] = fn2(SYM_Element, expr_copy(arg1), expr_new_symbol(SYM_Integers));
    else
        a[k++] = expr_copy(arg1);
    for (size_t i = 0; i < nopts; i++) a[k++] = expr_copy(opts[i]);
    const char* h = is_max ? SYM_NMaximize : SYM_NMinimize;
    Expr* call = expr_new_function(expr_new_symbol(h), a, k);
    free(a);
    return eval_and_free(call);
}

/* ------------------------------------------------------------------ *
 *  Dispatcher                                                         *
 * ------------------------------------------------------------------ */

static Expr* mz_run(Expr* res, const char* head, bool is_max) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** A = res->data.function.args;

    size_t pos_end = argc;
    while (pos_end > 0 && mz_is_option_arg(A[pos_end - 1])) pos_end--;
    if (pos_end < 2 || pos_end > 3) {
        mth_message(head, "argt",
                "needs 2 or 3 positional arguments; got %zu", pos_end);
        return NULL;
    }

    Expr* arg0 = A[0];
    Expr* arg1 = A[1];
    Expr* dom  = (pos_end == 3) ? A[2] : NULL;
    const char* domname = SYM_Reals;
    if (dom) {
        if (dom->type == EXPR_SYMBOL &&
            (dom->data.symbol.name == SYM_Reals ||
             dom->data.symbol.name == SYM_Integers ||
             dom->data.symbol.name == SYM_Complexes))
            domname = dom->data.symbol.name;
        else return NULL;                       /* unknown domain: leave unevaluated */
    }

    /* Split the objective from any constraints: a first argument {f, cons, ...}
     * whose second element is a relational/logical tree is a constrained
     * problem; a bare expression is unconstrained. */
    Expr* f = arg0;
    Expr* cons = NULL; bool cons_owned = false;
    if (arg0->type == EXPR_FUNCTION && mz_is_sym(arg0->data.function.head, SYM_List)) {
        size_t m = arg0->data.function.arg_count;
        if (m >= 2 && mz_is_constraint_tree(arg0->data.function.args[1])) {
            f = arg0->data.function.args[0];
            if (m == 2) cons = arg0->data.function.args[1];
            else {
                Expr** cc = (Expr**)malloc(sizeof(Expr*) * (m - 1));
                for (size_t i = 0; i < m - 1; i++)
                    cc[i] = expr_copy(arg0->data.function.args[1 + i]);
                cons = expr_new_function(expr_new_symbol(SYM_And), cc, m - 1);
                free(cc);
                cons_owned = true;
            }
        } else {
            return NULL;                        /* vector objective: not supported */
        }
    }

    Expr** vars = NULL; size_t n = 0;
    if (!(arg1->type == EXPR_SYMBOL ||
          (arg1->type == EXPR_FUNCTION && mz_is_sym(arg1->data.function.head, SYM_List)))) {
        if (cons_owned) expr_free(cons);
        return NULL;
    }
    if (arg1->type == EXPR_SYMBOL) {
        vars = (Expr**)malloc(sizeof(Expr*));
        vars[0] = expr_copy(arg1);
        n = 1;
    } else {
        size_t m = arg1->data.function.arg_count;
        if (m == 0) { if (cons_owned) expr_free(cons); return NULL; }
        for (size_t i = 0; i < m; i++)
            if (arg1->data.function.args[i]->type != EXPR_SYMBOL) {
                if (cons_owned) expr_free(cons);
                return NULL;                    /* indexed/vector vars: not supported */
            }
        vars = (Expr**)malloc(sizeof(Expr*) * m);
        for (size_t i = 0; i < m; i++) vars[i] = expr_copy(arg1->data.function.args[i]);
        n = m;
    }

    /* Maximize minimises -f, then negates the reported value. */
    Expr* fobj = is_max ? eval_and_free(fn2(SYM_Times, mk_int(-1), expr_copy(f)))
                        : expr_copy(f);

    Expr* result = NULL;

    if (mz_has_inexact(fobj) || (cons && mz_has_inexact(cons))) {
        /* Inexact input: hand the ORIGINAL problem to NMinimize/NMaximize. */
        result = mz_numeric_fallback(arg0, arg1, domname,
                                     &A[pos_end], argc - pos_end, is_max);
    } else if (domname == SYM_Reals) {
        bool unconstrained = (cons == NULL);
        if (unconstrained && n == 1)
            result = mz_univar_poly(fobj, vars[0], head);
        else
            result = mz_exact_poly(fobj, cons, vars, n, head);
        if (result && is_max) result = mz_negate_value(result);
    }
    /* Integers / Complexes with exact input: deferred — leave unevaluated. */

    expr_free(fobj);
    for (size_t i = 0; i < n; i++) expr_free(vars[i]);
    free(vars);
    if (cons_owned) expr_free(cons);
    return result;
}

Expr* builtin_minimize(Expr* res) { return mz_run(res, "Minimize", false); }
Expr* builtin_maximize(Expr* res) { return mz_run(res, "Maximize", true); }

/* ------------------------------------------------------------------ *
 *  Registration (called from core_init, among the calculus modules)   *
 * ------------------------------------------------------------------ */

void minimize_init(void) {
    /* Protected, NOT HoldAll (matching Mathematica's Attributes[Minimize] ==
     * {Protected, ReadProtected}): the variables are unbound symbols that
     * evaluate to themselves, so the objective and constraints arrive in
     * symbolic form and the exact engine proceeds by value-preserving
     * substitution (D, Solve, Reduce), needing no Block binding. */
    symtab_add_builtin("Minimize", builtin_minimize);
    symtab_get_def("Minimize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Maximize", builtin_maximize);
    symtab_get_def("Maximize")->attributes |= ATTR_PROTECTED;
}
