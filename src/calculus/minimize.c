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
#include <time.h>               /* clock() — C99, for the per-call deadline */

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

/* Read an integer-valued option (name -> k) out of the trailing option args.
 * Returns `dflt` when absent or malformed, and 0 ("no limit") for name ->
 * Infinity. A later occurrence wins (matching OptionValue semantics). */
static int mz_option_int(Expr* const* opts, size_t nopts,
                         const char* name, int dflt) {
    int found = dflt;
    for (size_t i = 0; i < nopts; i++) {
        const Expr* o = opts[i];
        if (!mz_is_option_arg(o) || o->data.function.arg_count != 2) continue;
        const Expr* lhs = o->data.function.args[0];
        const Expr* rhs = o->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL || strcmp(lhs->data.symbol.name, name) != 0)
            continue;
        if (rhs->type == EXPR_INTEGER)      found = (int)rhs->data.integer;
        else if (rhs->type == EXPR_REAL)    found = (int)rhs->data.real;
        else if (mz_is_sym(rhs, SYM_Infinity)) found = 0;   /* unlimited */
    }
    return found;
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

/* Bounded evaluation: run `e` under TimeConstrained[e, budget, $Aborted] so a
 * blow-up in an internal Solve/Reduce probe declines gracefully instead of
 * hanging. `budget` <= 0 means "no limit" (TimeConstraint -> Infinity). `e` is
 * consumed. On timeout the result is the $Aborted symbol, which every caller
 * treats as "unclean / not proven" — never as a positive answer — so the
 * soundness invariant is preserved (a timeout only ever *removes* a candidate or
 * a certificate, never fabricates one). Nesting clamps correctly: an outer
 * TimeConstrained[Minimize[...], T] still bounds these inner probes. */
static Expr* mz_beval(Expr* e, int budget) {
    if (budget <= 0) return eval_and_free(e);
    Expr* g = expr_new_function(expr_new_symbol(SYM_TimeConstrained),
                  (Expr*[]){ e, mk_int(budget),
                             expr_new_symbol(SYM_DollarAborted) }, 3);
    return eval_and_free(g);
}

/* Per-CALL deadline: a problem with many active-set probes must not cost
 * budget × probes. Given the call start `t0` and the total `budget` (seconds;
 * <= 0 = unlimited), return the seconds left to hand the NEXT probe: 0 means
 * unlimited, a positive value is the remaining budget, and -1 means the deadline
 * is spent (the caller declines). So the whole call is bounded by `budget`. */
static int mz_rem_budget(clock_t t0, int budget) {
    if (budget <= 0) return 0;                         /* unlimited */
    double el = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    int rem = budget - (int)el;
    return rem > 0 ? rem : -1;
}

/* Solve[system, {vars}, Reals], evaluated under the time budget — consumes
 * `system`. Internal probe: an unsolvable / parametric subproblem is expected
 * (and skipped), so Solve's own diagnostics (Solve::nsdim, ...) are suppressed;
 * a budget abort returns $Aborted, which mz_parse_points rejects as unclean. */
static Expr* mz_solve_real(Expr* system, Expr* const* vars, size_t n, int budget) {
    mth_msg_suppress_push();
    Expr* r = mz_beval(expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ system, mz_varlist(vars, n), expr_new_symbol(SYM_Reals) }, 3),
        budget);
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

/* 1 iff  A ⇒ P  is provable over the Reals (Reduce[A && !P] is False), else 0.
 * A budget abort returns $Aborted (not False) → 0 = "not proven" → decline. */
static int mz_entails(const Expr* A, const Expr* P, Expr* const* vars, size_t n,
                      int budget) {
    Expr* notP = eval_and_free(fn1(SYM_Not, expr_copy((Expr*)P)));
    Expr* stmt = (!A || mz_is_true(A)) ? notP
               : fn2(SYM_And, expr_copy((Expr*)A), notP);
    mth_msg_suppress_push();
    Expr* red = mz_beval(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ stmt, mz_varlist(vars, n), expr_new_symbol(SYM_Reals) }, 3),
        budget);
    mth_msg_suppress_pop();
    int r = mz_is_sym(red, SYM_False) ? 1 : 0;
    expr_free(red);
    return r;
}

/* True iff Reduce proves the constraint region empty over the Reals. A budget
 * abort returns $Aborted (not False) → false = "not provably empty" → decline. */
static bool mz_region_empty(const Expr* cons, Expr* const* vars, size_t n,
                            int budget) {
    mth_msg_suppress_push();
    Expr* red = mz_beval(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ expr_copy((Expr*)cons), mz_varlist(vars, n),
                   expr_new_symbol(SYM_Reals) }, 3), budget);
    mth_msg_suppress_pop();
    bool r = mz_is_sym(red, SYM_False);
    expr_free(red);
    return r;
}

/* True iff Reduce proves the equality variety {g == 0} lies inside some ball
 * Σ x_i^2 <= R^2 (rational R) — i.e. Reduce[g==0 && Σx_i^2 > R^2] === False for
 * some R on a geometric ladder. Every coefficient is rational, so this stays in
 * the CAD's supported regime (it never hits the algebraic-coefficient wall that
 * stops the lower-bound certificate). Bounded + closed (an equality variety is
 * closed) ⇒ compact. Each probe uses a small sub-budget so an unbounded variety
 * (no R ever certifies) cannot run away. */
static bool mz_bounded_region(const Expr* g, Expr* const* vars, size_t n,
                              int budget) {
    int sub = (budget > 0 && budget < 4) ? budget : 4;   /* per-probe cap */
    Expr** sq = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t j = 0; j < n; j++)
        sq[j] = fn2(SYM_Power, expr_copy(vars[j]), mk_int(2));
    Expr* sumsq = expr_new_function(expr_new_symbol(SYM_Plus), sq, n);
    free(sq);
    Expr* geq = fn2(SYM_Equal, expr_copy((Expr*)g), mk_int(0));
    bool bounded = false;
    for (int e = 2; e <= 20 && !bounded; e += 2) {        /* R^2 = 2^2 .. 2^20 */
        Expr* R2 = eval_and_free(fn2(SYM_Power, mk_int(2), mk_int(e)));
        Expr* gt = fn2(SYM_Greater, expr_copy(sumsq), R2);
        Expr* stmt = fn2(SYM_And, expr_copy(geq), gt);
        mth_msg_suppress_push();
        Expr* red = mz_beval(expr_new_function(expr_new_symbol(SYM_Reduce),
            (Expr*[]){ stmt, mz_varlist(vars, n),
                       expr_new_symbol(SYM_Reals) }, 3), sub);
        mth_msg_suppress_pop();
        if (mz_is_sym(red, SYM_False)) bounded = true;
        expr_free(red);
    }
    expr_free(sumsq);
    expr_free(geq);
    return bounded;
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
 *  Radical clearing for the lower-bound certificate                   *
 * ------------------------------------------------------------------ */

/* The CAD behind Reduce declines on irrational-algebraic input coefficients, so
 * the lower-bound certificate `cons => f >= bestv` cannot be decided when bestv
 * is irrational (e.g. 14 - 2 Sqrt[13]). We rewrite bestv into a RATIONAL-
 * coefficient system by replacing each radical with a fresh variable pinned by a
 * polynomial defining relation, then reduce over the enlarged variable set — a
 * statement logically identical to the original, which Reduce CAN decide. */
typedef struct {
    Expr** node; Expr** val;  size_t nrad, radcap;   /* dedup: radical node -> u^p */
    Expr** aux;  size_t naux,  auxcap;    /* fresh u symbols (owned copies) */
    Expr** defs; size_t ndefs, defcap;    /* defining relations (owned) */
} MzRadCtx;

static void mz_rad_push_aux(MzRadCtx* c, Expr* u) {
    if (c->naux == c->auxcap) { c->auxcap = c->auxcap ? c->auxcap * 2 : 4;
        c->aux = (Expr**)realloc(c->aux, sizeof(Expr*) * c->auxcap); }
    c->aux[c->naux++] = u;
}
static void mz_rad_push_def(MzRadCtx* c, Expr* d) {
    if (c->ndefs == c->defcap) { c->defcap = c->defcap ? c->defcap * 2 : 4;
        c->defs = (Expr**)realloc(c->defs, sizeof(Expr*) * c->defcap); }
    c->defs[c->ndefs++] = d;
}
/* Record that radical `node` is replaced by `val`; both copied in. */
static void mz_rad_push_pair(MzRadCtx* c, const Expr* node, const Expr* val) {
    if (c->nrad == c->radcap) { c->radcap = c->radcap ? c->radcap * 2 : 4;
        c->node = (Expr**)realloc(c->node, sizeof(Expr*) * c->radcap);
        c->val  = (Expr**)realloc(c->val,  sizeof(Expr*) * c->radcap); }
    c->node[c->nrad] = expr_copy((Expr*)node);
    c->val [c->nrad] = expr_copy((Expr*)val);
    c->nrad++;
}
static void mz_rad_free(MzRadCtx* c) {
    for (size_t i = 0; i < c->naux;  i++) expr_free(c->aux[i]);
    for (size_t i = 0; i < c->ndefs; i++) expr_free(c->defs[i]);
    for (size_t i = 0; i < c->nrad;  i++) { expr_free(c->node[i]); expr_free(c->val[i]); }
    free(c->aux); free(c->defs); free(c->node); free(c->val);
}

/* Replace each radical Power[c, p/q] (c a provably-POSITIVE constant, q >= 2) in
 * `e` with a fresh u: record u^q == c and u >= 0 (the principal real branch; for
 * c > 0 these are satisfiable and pin u = c^(1/q) exactly), and let the radical's
 * value c^(p/q) become u^p. Nested radicals in the base are cleared first, so
 * every emitted relation is polynomial in the fresh variables. A radical whose
 * base is not a provably-positive constant is left untouched — the certificate
 * then simply fails to prove and the caller declines, so this never fabricates a
 * proof. Returns a radical-cleared copy of `e`. */
static Expr* mz_clear_radicals(const Expr* e, MzRadCtx* c) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (mz_is_sym(e->data.function.head, SYM_Power) &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* ex   = e->data.function.args[1];
        if (ex->type == EXPR_FUNCTION &&
            mz_is_sym(ex->data.function.head, SYM_Rational) &&
            ex->data.function.arg_count == 2 &&
            ex->data.function.args[0]->type == EXPR_INTEGER &&
            ex->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t p = ex->data.function.args[0]->data.integer;
            int64_t q = ex->data.function.args[1]->data.integer;
            if (q >= 2 && rru_sign_of(base) == 1) {       /* positive base only */
                /* Dedup: a radical already cleared reuses its variable, so a
                 * value with the same radical twice (common for an optimum like
                 * 14 - 2 Sqrt[13]) does NOT add a second CAD variable — the extra
                 * dimension makes the certificate blow up. */
                for (size_t i = 0; i < c->nrad; i++)
                    if (expr_eq(e, c->node[i])) return expr_copy(c->val[i]);
                Expr* baseC = mz_clear_radicals(base, c);        /* nested first */
                Expr* u = mz_fresh_symbol();
                mz_rad_push_def(c, fn2(SYM_Equal,
                    fn2(SYM_Power, expr_copy(u), mk_int(q)), baseC));
                mz_rad_push_def(c, fn2(SYM_GreaterEqual, expr_copy(u), mk_int(0)));
                mz_rad_push_aux(c, expr_copy(u));
                Expr* val = (p == 1) ? expr_copy(u)
                                     : fn2(SYM_Power, expr_copy(u), mk_int(p));
                expr_free(u);
                mz_rad_push_pair(c, e, val);
                return val;
            }
        }
    }
    size_t m = e->data.function.arg_count;
    Expr* headC = mz_clear_radicals(e->data.function.head, c);
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) a[i] = mz_clear_radicals(e->data.function.args[i], c);
    Expr* r = expr_new_function(headC, a, m);
    free(a);
    return r;
}

/* Prove the global lower bound  closure(cons) => f >= bestv  over the Reals.
 * When bestv is rational (or an atom we cannot clear, e.g. a Root[]) the direct
 * certificate is exact and is the only attempt. When bestv carries radicals the
 * direct certificate is skipped — the CAD would burn the whole budget declining
 * on the algebraic coefficient — and we reduce the logically-identical
 * radical-cleared system over vars ∪ aux instead. Returns 1 iff proven. */
static int mz_prove_lower_bound(const Expr* CL, const Expr* f, const Expr* bestv,
                                Expr* const* vars, size_t n, int budget) {
    MzRadCtx c = {0};
    Expr* bestvC = mz_clear_radicals(bestv, &c);
    int proven;
    if (c.naux == 0) {
        Expr* ge = fn2(SYM_GreaterEqual, expr_copy((Expr*)f), expr_copy((Expr*)bestv));
        proven = mz_entails(CL, ge, vars, n, budget);
        expr_free(ge);
    } else {
        size_t nall = n + c.naux;
        Expr** vall = (Expr**)malloc(sizeof(Expr*) * nall);
        for (size_t j = 0; j < n; j++)       vall[j]     = expr_copy(vars[j]);
        for (size_t k = 0; k < c.naux; k++)  vall[n + k] = expr_copy(c.aux[k]);
        /* A' = closure(cons) && defining relations */
        size_t na = 1 + c.ndefs;
        Expr** parts = (Expr**)malloc(sizeof(Expr*) * na);
        parts[0] = expr_copy((Expr*)CL);
        for (size_t k = 0; k < c.ndefs; k++) parts[1 + k] = expr_copy(c.defs[k]);
        Expr* Aext = expr_new_function(expr_new_symbol(SYM_And), parts, na);
        free(parts);
        Expr* geC = fn2(SYM_GreaterEqual, expr_copy((Expr*)f), expr_copy(bestvC));
        proven = mz_entails(Aext, geC, vall, nall, budget);
        expr_free(geC); expr_free(Aext);
        for (size_t j = 0; j < nall; j++) expr_free(vall[j]);
        free(vall);
    }
    expr_free(bestvC);
    mz_rad_free(&c);
    return proven;
}

/* ------------------------------------------------------------------ *
 *  Radical / fractional-power OBJECTIVES and CONSTRAINTS (Tier 2a)     *
 * ------------------------------------------------------------------ */

/* Rewrite each radical Power[g, p/q] (q >= 2, g any expression) in an objective
 * or constraint as u^p, adjoining the defining relations u^q == g_cleared and
 * u >= 0 to the problem. Because u >= 0 and u^q == g together force g >= 0, the
 * enlarged feasible region is exactly where the real principal branch is defined
 * (Mathilda's Power[g, p/q] is complex for g < 0 at a non-integer exponent), so
 * the value and the real domain are both preserved. Nested radicals in the base
 * are cleared first and identical radicals are deduplicated (a shared variable),
 * so the result is polynomial in the original vars ∪ the fresh u's. Unlike
 * mz_clear_radicals (the certificate, which needs a provably-positive CONSTANT
 * base for soundness) this accepts variable bases, because here the relations are
 * ADDED to the problem rather than used to refute a bound. */
static Expr* mz_rad_walk(const Expr* e, MzRadCtx* c) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (mz_is_sym(e->data.function.head, SYM_Power) &&
        e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        if (ex->type == EXPR_FUNCTION &&
            mz_is_sym(ex->data.function.head, SYM_Rational) &&
            ex->data.function.arg_count == 2 &&
            ex->data.function.args[0]->type == EXPR_INTEGER &&
            ex->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t p = ex->data.function.args[0]->data.integer;
            int64_t q = ex->data.function.args[1]->data.integer;
            if (q >= 2) {
                for (size_t i = 0; i < c->nrad; i++)
                    if (expr_eq(e, c->node[i])) return expr_copy(c->val[i]);
                Expr* gC = mz_rad_walk(e->data.function.args[0], c);  /* nested first */
                Expr* u = mz_fresh_symbol();
                mz_rad_push_def(c, fn2(SYM_Equal,
                    fn2(SYM_Power, expr_copy(u), mk_int(q)), gC));      /* u^q == g */
                mz_rad_push_def(c, fn2(SYM_GreaterEqual, expr_copy(u), mk_int(0)));
                mz_rad_push_aux(c, expr_copy(u));
                Expr* val = (p == 1) ? expr_copy(u)
                                     : fn2(SYM_Power, expr_copy(u), mk_int(p));
                expr_free(u);
                mz_rad_push_pair(c, e, val);
                return val;
            }
        }
    }
    size_t m = e->data.function.arg_count;
    Expr* headC = mz_rad_walk(e->data.function.head, c);
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) a[i] = mz_rad_walk(e->data.function.args[i], c);
    Expr* r = expr_new_function(headC, a, m);
    free(a);
    return r;
}

/* Drop the auxiliary-variable rules (u -> ...) from a {value, {rules}} result,
 * leaving only the original problem variables. Consumes and rebuilds `result`. */
static Expr* mz_strip_aux_rules(Expr* result, Expr* const* aux, size_t naux) {
    if (!result || result->type != EXPR_FUNCTION ||
        !mz_is_sym(result->data.function.head, SYM_List) ||
        result->data.function.arg_count != 2) return result;
    const Expr* value = result->data.function.args[0];
    const Expr* rules = result->data.function.args[1];
    if (!rules || rules->type != EXPR_FUNCTION ||
        !mz_is_sym(rules->data.function.head, SYM_List)) return result;
    size_t m = rules->data.function.arg_count;
    Expr** keep = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
    size_t nk = 0;
    for (size_t i = 0; i < m; i++) {
        const Expr* rule = rules->data.function.args[i];
        bool is_aux = false;
        if (rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2 &&
            rule->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* nm = rule->data.function.args[0]->data.symbol.name;
            for (size_t k = 0; k < naux; k++)
                if (nm == aux[k]->data.symbol.name) { is_aux = true; break; }
        }
        if (!is_aux) keep[nk++] = expr_copy((Expr*)rule);
    }
    Expr* newrules = expr_new_function(expr_new_symbol(SYM_List), keep, nk);
    free(keep);
    Expr* newval = expr_copy((Expr*)value);
    expr_free(result);
    return fn2(SYM_List, newval, newrules);
}

/* ------------------------------------------------------------------ *
 *  (a) Univariate unconstrained polynomial                            *
 * ------------------------------------------------------------------ */

static Expr* mz_univar_poly(const Expr* f, Expr* x, const char* head, int budget) {
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
    Expr* L = mz_solve_real(eqn, vars1, 1, budget);
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
 *  Univariate Abs / piecewise-polynomial objective (Tier 2b)          *
 * ------------------------------------------------------------------ */

/* Collect the distinct Abs[g] arguments g (g polynomial in x) appearing in e. */
static void mz_collect_abs(const Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (mz_is_sym(e->data.function.head, SYM_Abs) &&
        e->data.function.arg_count == 1 &&
        rru_is_polynomial(e->data.function.args[0], x)) {
        const Expr* g = e->data.function.args[0];
        for (size_t i = 0; i < *n; i++) if (expr_eq(g, (*arr)[i])) return;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
            *arr = (Expr**)realloc(*arr, sizeof(Expr*) * *cap); }
        (*arr)[(*n)++] = expr_copy((Expr*)g);
        return;                              /* g is polynomial: no nested Abs */
    }
    mz_collect_abs(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        mz_collect_abs(e->data.function.args[i], x, arr, n, cap);
}

/* Replace each Abs[ag[i]] by sg[i]*ag[i] (sg[i] in {+1,-1}); other Abs untouched. */
static Expr* mz_resolve_abs(const Expr* e, Expr* const* ag, const int* sg, size_t k) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (mz_is_sym(e->data.function.head, SYM_Abs) && e->data.function.arg_count == 1)
        for (size_t i = 0; i < k; i++)
            if (expr_eq(e->data.function.args[0], ag[i])) {
                Expr* g = expr_copy(ag[i]);
                return (sg[i] < 0) ? eval_and_free(fn2(SYM_Times, mk_int(-1), g)) : g;
            }
    size_t m = e->data.function.arg_count;
    Expr* h = mz_resolve_abs(e->data.function.head, ag, sg, k);
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) a[i] = mz_resolve_abs(e->data.function.args[i], ag, sg, k);
    Expr* r = expr_new_function(h, a, m);
    free(a);
    return r;
}

/* Degree and leading-coeff sign of polynomial p in x; 1 on success, 0 if the
 * coefficient list is unavailable or the leading sign is undecidable. */
static int mz_poly_tail(const Expr* p, Expr* x, int* deg, int* lsign) {
    Expr* cl = eval_and_free(fn2(SYM_CoefficientList, expr_copy((Expr*)p), expr_copy(x)));
    if (!cl || cl->type != EXPR_FUNCTION || !mz_is_sym(cl->data.function.head, SYM_List) ||
        cl->data.function.arg_count == 0) { expr_free(cl); return 0; }
    size_t d1 = cl->data.function.arg_count;
    int s = rru_sign_of(cl->data.function.args[d1 - 1]);
    expr_free(cl);
    if (s == -2) return 0;
    *deg = (int)d1 - 1; *lsign = s;
    return 1;
}

/* Sign of polynomial g at x = tp (exact): +1/-1/0, or -2 undecidable. tp consumed. */
static int mz_sign_at(const Expr* g, Expr* x, Expr* tp) {
    Expr* xv[1] = { x }; Expr* vv[1] = { tp };
    Expr* val = mz_subst_eval(g, xv, vv, 1);
    int s = rru_sign_of(val);
    expr_free(val); expr_free(tp);
    return s;
}

static void mz_push_cand1(Expr**** cand, size_t* nc, size_t* cap, Expr* val) {
    if (*nc == *cap) { *cap = *cap ? *cap * 2 : 8;
        *cand = (Expr***)realloc(*cand, sizeof(Expr**) * *cap); }
    (*cand)[*nc] = (Expr**)malloc(sizeof(Expr*));
    (*cand)[*nc][0] = val;
    (*nc)++;
}

/* Univariate objective that is piecewise-polynomial through Abs[poly] terms.
 * The objective is smooth on each sign cell of the Abs arguments; its global
 * minimum is attained at a breakpoint (a real root of some Abs argument) or at a
 * stationary point of a piece. We therefore evaluate the ORIGINAL objective at
 * the union {breakpoints} ∪ {real roots of each piece's derivative} and take the
 * exact least — sound because every candidate is a genuine real point and the
 * true minimiser is always in the set. Unboundedness is read off the two end
 * pieces' leading terms. Declines (NULL) if a root solve or an exact sign/compare
 * is not clean, or if an Abs wraps a non-polynomial argument. */
static Expr* mz_univar_piecewise(const Expr* f, Expr* x, const char* head, int budget) {
    Expr* vars1[1] = { x };
    Expr** ag = NULL; size_t nag = 0, agcap = 0;
    mz_collect_abs(f, x, &ag, &nag, &agcap);
    if (nag == 0) { free(ag); return NULL; }        /* no Abs: ordinary path */

    Expr* result = NULL;
    int* signs = (int*)malloc(sizeof(int) * nag);
    Expr** bp = NULL; size_t nbp = 0;
    Expr*** cand = NULL; size_t ncand = 0, candcap = 0;
    bool unbounded = false;
    clock_t t0 = clock();                           /* per-call deadline anchor */

    /* Resolving every Abs (+1) must leave a polynomial — else an Abs wraps a
     * non-polynomial argument and the piecewise reduction does not apply. */
    for (size_t i = 0; i < nag; i++) signs[i] = 1;
    { Expr* chk = eval_and_free(mz_resolve_abs(f, ag, signs, nag));
      bool poly = rru_is_polynomial(chk, x); expr_free(chk);
      if (!poly) goto cleanup; }

    /* Breakpoints: the sorted, deduplicated real roots of the Abs arguments. */
    { Expr** raw = NULL; size_t nraw = 0, rawcap = 0; bool bad = false;
      for (size_t i = 0; i < nag && !bad; i++) {
        Expr* eqn = fn2(SYM_Equal, expr_copy(ag[i]), mk_int(0));
        int rem = mz_rem_budget(t0, budget);
        Expr* L = (rem < 0) ? (expr_free(eqn), NULL) : mz_solve_real(eqn, vars1, 1, rem);
        Expr*** pts = NULL; size_t npts = 0;
        if (!L || !mz_parse_points(L, vars1, 1, vars1, 1, &pts, &npts)) bad = true;
        expr_free(L);
        if (!bad) for (size_t p = 0; p < npts; p++) {
            if (nraw == rawcap) { rawcap = rawcap ? rawcap * 2 : 8;
                raw = (Expr**)realloc(raw, sizeof(Expr*) * rawcap); }
            raw[nraw++] = expr_copy(pts[p][0]);
        }
        mz_free_points(pts, npts, 1);
      }
      for (size_t i = 1; i < nraw && !bad; i++) { Expr* key = raw[i]; size_t j = i;
        while (j > 0) { int c = rru_sign_compare(raw[j-1], key);
            if (c == -2) { bad = true; break; }
            if (c > 0) { raw[j] = raw[j-1]; j--; } else break; }
        raw[j] = key; }
      if (!bad) { bp = (Expr**)malloc(sizeof(Expr*) * (nraw ? nraw : 1));
        for (size_t i = 0; i < nraw; i++) {
            if (nbp > 0 && rru_sign_compare(bp[nbp-1], raw[i]) == 0) { expr_free(raw[i]); continue; }
            bp[nbp++] = raw[i];
        } }
      else for (size_t q = 0; q < nraw; q++) expr_free(raw[q]);
      free(raw);
      if (bad) goto cleanup;
    }

    /* No real breakpoints: the Abs signs are globally constant, so a single
     * polynomial piece covers the whole line — hand it to the plain engine. */
    if (nbp == 0) {
        for (size_t i = 0; i < nag; i++)
            if ((signs[i] = mz_sign_at(ag[i], x, mk_int(0))) == -2 || signs[i] == 0) goto cleanup;
        Expr* piece = eval_and_free(mz_resolve_abs(f, ag, signs, nag));
        result = mz_univar_poly(piece, x, head, budget);
        expr_free(piece);
        goto cleanup;
    }

    for (size_t i = 0; i < nbp; i++) mz_push_cand1(&cand, &ncand, &candcap, expr_copy(bp[i]));

    /* Each interval j = 0..nbp: resolve its piece and add the piece-derivative's
     * real roots; read unboundedness off the two unbounded end pieces. */
    for (size_t j = 0; j <= nbp; j++) {
        Expr* tp;
        if (j == 0) {                                   /* bp[0] - 1 */
            Expr* one = mk_int(1); tp = mz_sub(bp[0], one); expr_free(one);
        } else if (j == nbp) {                          /* bp[nbp-1] + 1 */
            tp = eval_and_free(fn2(SYM_Plus, expr_copy(bp[nbp-1]), mk_int(1)));
        } else {
            tp = rru_rational_between(bp[j-1], bp[j]);
        }
        if (!tp) goto cleanup;
        bool sbad = false;
        for (size_t i = 0; i < nag; i++) {
            signs[i] = mz_sign_at(ag[i], x, expr_copy(tp));
            if (signs[i] == -2 || signs[i] == 0) sbad = true;
        }
        expr_free(tp);
        if (sbad) goto cleanup;

        Expr* piece = eval_and_free(mz_resolve_abs(f, ag, signs, nag));
        if (!mz_contains_var(piece, vars1, 1)) { expr_free(piece); continue; }
        int deg = 0, ls = 0;                       /* constant piece: no crit pts,
                                                    * bounded; endpoint is a cand */
        if (!mz_poly_tail(piece, x, &deg, &ls)) { expr_free(piece); goto cleanup; }

        /* unbounded-below detection on the two semi-infinite end pieces */
        if (deg >= 1) {
            if (j == 0   && ((deg % 2 == 1 && ls > 0) || (deg % 2 == 0 && ls < 0))) unbounded = true;
            if (j == nbp && ls < 0) unbounded = true;
        }

        if (deg >= 1) {           /* non-constant piece: add stationary points */
            Expr* fp = mz_deriv(piece, x);
            Expr* eqn = fn2(SYM_Equal, fp, mk_int(0));
            int rem = mz_rem_budget(t0, budget);
            Expr* L = (rem < 0) ? (expr_free(eqn), NULL) : mz_solve_real(eqn, vars1, 1, rem);
            Expr*** pts = NULL; size_t npts = 0;
            int ok = L ? mz_parse_points(L, vars1, 1, vars1, 1, &pts, &npts) : 0;
            expr_free(L);
            if (!ok) { expr_free(piece); goto cleanup; }
            for (size_t p = 0; p < npts; p++)
                mz_push_cand1(&cand, &ncand, &candcap, expr_copy(pts[p][0]));
            mz_free_points(pts, npts, 1);
        }
        expr_free(piece);
    }

    if (unbounded) {
        mth_message(head, "natt",
                "The %s is not attained at any point satisfying the given "
                "constraints.", mz_noun(head));
        result = mz_result_indet(mz_neg_inf(), vars1, 1);
        goto cleanup;
    }

    { Expr* bestv = NULL; size_t bi = 0;
      if (!mz_pick_min(f, vars1, 1, cand, ncand, &bestv, &bi)) goto cleanup;
      Expr* pt[1] = { cand[bi][0] };
      result = mz_result(bestv, vars1, pt, 1);      /* bestv consumed; pt copied */
    }

cleanup:
    for (size_t i = 0; i < nag; i++) expr_free(ag[i]);
    free(ag);
    free(signs);
    for (size_t i = 0; i < nbp; i++) expr_free(bp[i]);
    free(bp);
    mz_free_points(cand, ncand, 1);
    return result;
}

/* ------------------------------------------------------------------ *
 *  (b)/(c) Multivariate / constrained polynomial over the Reals       *
 * ------------------------------------------------------------------ */

#define MZ_MAX_INEQ 6               /* 2^MZ_MAX_INEQ active-set subproblems */
#define MZ_DEFAULT_TIMECONSTRAINT 30  /* seconds per internal Solve/Reduce probe */

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
                           Expr* const* vars, size_t n, const char* head,
                           int budget) {
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
    bool kkt_all_clean = true;                    /* every active-set solve was clean */
    bool shortcut_proven = false;                 /* compact-region EVT certificate */
    clock_t t0 = clock();                         /* per-call deadline anchor */

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

        int rem = mz_rem_budget(t0, budget);
        if (rem < 0) {                          /* per-call deadline spent */
            for (size_t c = 0; c < nact; c++) expr_free(lam[c]);
            for (size_t j = 0; j < nall; j++) expr_free(allv[j]);
            free(lam); free(allv); free(actz); free(actdz);
            expr_free(system);
            goto cleanup;                       /* incomplete enumeration: decline */
        }
        Expr* L = mz_solve_real(system, allv, nall, rem);
        Expr*** pp = NULL; size_t npp = 0;
        int ok = mz_parse_points(L, vars, n, allv, nall, &pp, &npp);
        expr_free(L);
        for (size_t c = 0; c < nact; c++) expr_free(lam[c]);
        for (size_t j = 0; j < nall; j++) expr_free(allv[j]);
        free(lam); free(allv); free(actz); free(actdz);

        if (!ok) { kkt_all_clean = false; continue; }  /* parametric/unclean: skip */
        for (size_t i = 0; i < npp; i++) {
            if (ncand == candcap) { candcap = candcap ? candcap * 2 : 8;
                cand = (Expr***)realloc(cand, sizeof(Expr**) * candcap); }
            cand[ncand++] = pp[i];
        }
        free(pp);                               /* rows moved into cand */
    }

    /* Compact-region extreme-value shortcut (pure single-equality case).
     * For a smooth equality variety V = {g == 0}, every local extremum of f on
     * V is either a REGULAR point (the KKT/Lagrange system, already solved
     * above) or a SINGULAR point of V ({g == 0, ∇g == 0}) — Fritz–John. If V is
     * compact (closed, since it is an equality set, + bounded, certified by a
     * rational ball probe) then f attains its global minimum on V, and it is the
     * least value over the enumerated regular ∪ singular candidates. In that
     * case the Reduce lower-bound certificate below is REDUNDANT and is skipped
     * (which also dodges the algebraic-coefficient wall that makes it decline on
     * irrational optima). Soundness rests on the solvenlsys contract (a complete
     * finite solution set or a decline — never a false empty set): the shortcut
     * fires only when BOTH the regular and the singular solves are clean and the
     * region is provably bounded; any gap falls through to the certificate. */
    if (nin == 0 && neq == 1 && kkt_all_clean) {
        size_t nsys = 1 + n;
        Expr** sys = (Expr**)malloc(sizeof(Expr*) * nsys);
        sys[0] = fn2(SYM_Equal, expr_copy(EQ[0]), mk_int(0));
        for (size_t j = 0; j < n; j++)
            sys[1 + j] = fn2(SYM_Equal, expr_copy(dzEQ[0][j]), mk_int(0));
        Expr* singsys = expr_new_function(expr_new_symbol(SYM_List), sys, nsys);
        free(sys);
        int rem = mz_rem_budget(t0, budget);
        Expr* Ls = (rem < 0) ? (expr_free(singsys), NULL)
                             : mz_solve_real(singsys, vars, n, rem);
        Expr*** sp = NULL; size_t nsp = 0;
        int sing_clean = Ls ? mz_parse_points(Ls, vars, n, vars, n, &sp, &nsp) : 0;
        expr_free(Ls);
        if (sing_clean) {
            for (size_t i = 0; i < nsp; i++) {
                if (ncand == candcap) { candcap = candcap ? candcap * 2 : 8;
                    cand = (Expr***)realloc(cand, sizeof(Expr**) * candcap); }
                cand[ncand++] = sp[i];
            }
            free(sp);                           /* rows moved into cand */
            int rem2 = mz_rem_budget(t0, budget);
            if (rem2 >= 0 && mz_bounded_region(EQ[0], vars, n, rem2))
                shortcut_proven = true;
        }
    }

    /* Feasibility filter against the closure of the region. */
    CL = constrained ? mz_closure(consN) : expr_new_symbol(SYM_True);
    feas = (Expr***)malloc(sizeof(Expr**) * (ncand ? ncand : 1));
    for (size_t i = 0; i < ncand; i++)
        if (mz_feasible_at(CL, vars, cand[i], n) == 1) feas[nfeas++] = cand[i];

    if (nfeas == 0) {
        int reme = mz_rem_budget(t0, budget);
        if (constrained && reme >= 0 && mz_region_empty(consN, vars, n, reme)) {
            mth_message(head, "infeas",
                    "The constraints are infeasible; the feasible region is empty.");
            result = mz_result_indet(mz_pos_inf(), vars, n);
        }
        goto cleanup;                           /* else: unbounded/undecided → decline */
    }

    {
        Expr* bestv = NULL; size_t bi = 0;
        if (!mz_pick_min(f, vars, n, feas, nfeas, &bestv, &bi)) goto cleanup;

        /* Global lower-bound certificate: closure(cons) ⇒ f >= bestv. Skipped
         * when the compact-region shortcut already proved global optimality by
         * the extreme-value theorem; otherwise proved by Reduce, with radicals in
         * bestv cleared to a rational-coefficient system so an irrational optimum
         * is decidable. */
        int remc = mz_rem_budget(t0, budget);
        int proven = shortcut_proven ? 1
                   : (remc < 0 ? 0 : mz_prove_lower_bound(CL, f, bestv, vars, n, remc));
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
 *  Exact Integers-domain optimisation over a finite integer set       *
 * ------------------------------------------------------------------ */

/* Minimise f over the integer points enumerated by Solve[cons, vars, Integers]
 * (equality / Diophantine regions). Returns the exact least over a clean, finite
 * set of all-integer tuples, else NULL (parametric / unbounded / undecidable). */
static Expr* mz_integer_solve(const Expr* f, const Expr* cons, Expr* const* vars,
                              size_t n, int budget) {
    mth_msg_suppress_push();
    Expr* L = mz_beval(expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ expr_copy((Expr*)cons), mz_varlist(vars, n),
                   expr_new_symbol(SYM_Integers) }, 3), budget);
    mth_msg_suppress_pop();
    Expr*** pts = NULL; size_t npts = 0;
    int ok = mz_parse_points(L, vars, n, vars, n, &pts, &npts);
    expr_free(L);
    if (!ok) return NULL;
    if (npts == 0) { mz_free_points(pts, npts, n); return NULL; }
    for (size_t i = 0; i < npts; i++)               /* reject a parametric family */
        for (size_t j = 0; j < n; j++)
            if (!expr_is_integer_like(pts[i][j])) {
                mz_free_points(pts, npts, n); return NULL;
            }
    Expr* bestv = NULL; size_t bi = 0;
    if (!mz_pick_min(f, vars, n, pts, npts, &bestv, &bi)) {
        mz_free_points(pts, npts, n); return NULL;
    }
    Expr* r = mz_result(bestv, vars, pts[bi], n);   /* bestv consumed; pt copied */
    mz_free_points(pts, npts, n);
    return r;
}

/* Continuous bound of a single variable v over cons as an integer: Ceiling of
 * the relaxed min (is_max=false) or Floor of the relaxed max (is_max=true),
 * reusing the Reals engine. 1 on success; 0 if unbounded, undecided, or the
 * bound does not fit a machine integer. */
static int mz_cont_int_bound(const Expr* v, const Expr* cons, Expr* const* vars,
                             size_t n, bool is_max, int budget, int64_t* out) {
    Expr* obj = is_max ? eval_and_free(fn2(SYM_Times, mk_int(-1), expr_copy((Expr*)v)))
                       : expr_copy((Expr*)v);
    mth_msg_suppress_push();
    Expr* res = mz_exact_poly(obj, cons, vars, n, "Minimize", budget);
    mth_msg_suppress_pop();
    expr_free(obj);
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) {
        expr_free(res); return 0; }
    Expr* val = expr_copy(res->data.function.args[0]);          /* min of obj */
    expr_free(res);
    if (is_max) val = eval_and_free(fn2(SYM_Times, mk_int(-1), val));  /* max(v) */
    Expr* ib = eval_and_free(fn1(is_max ? SYM_Floor : SYM_Ceiling, val));
    int okr = (ib && ib->type == EXPR_INTEGER);
    if (okr) *out = ib->data.integer;
    expr_free(ib);
    return okr;
}

/* Minimise f over the integer points of a box whose per-variable bounds come
 * from the continuous relaxation, keeping only those satisfying cons. The box is
 * enumerated exactly (sound). Declines if a variable cannot be bounded, the box
 * exceeds a size cap, or an exact feasibility / comparison is undecided. */
static Expr* mz_integer_box(const Expr* f, const Expr* cons, Expr* const* vars,
                            size_t n, int budget) {
    int64_t* lo = (int64_t*)malloc(sizeof(int64_t) * n);
    int64_t* hi = (int64_t*)malloc(sizeof(int64_t) * n);
    int64_t* cur = (int64_t*)malloc(sizeof(int64_t) * n);
    Expr* bestv = NULL; Expr** bestpt = NULL;
    clock_t t0 = clock();                            /* deadline across all bounds */
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        int r1 = mz_rem_budget(t0, budget);          /* remaining before the min */
        if (r1 < 0 || !mz_cont_int_bound(vars[i], cons, vars, n, false, r1, &lo[i]))
            { ok = false; break; }
        int r2 = mz_rem_budget(t0, budget);          /* refreshed before the max */
        if (r2 < 0 || !mz_cont_int_bound(vars[i], cons, vars, n, true, r2, &hi[i]) ||
            hi[i] < lo[i]) { ok = false; break; }
    }
    if (ok) { double sz = 1.0;
        for (size_t i = 0; i < n; i++) {
            sz *= (double)(hi[i] - lo[i] + 1);
            if (sz > 2.0e6) { ok = false; break; } } }
    if (!ok) { free(lo); free(hi); free(cur); return NULL; }

    for (size_t i = 0; i < n; i++) cur[i] = lo[i];
    bool aborted = false, done = false;
    while (!done && !aborted) {
        Expr** pt = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) pt[i] = mk_int(cur[i]);
        bool kept = false;
        int fe = mz_feasible_at(cons, vars, pt, n);
        if (fe == -1) aborted = true;                /* undecided: stay sound */
        else if (fe == 1) {
            Expr* val = mz_subst_eval(f, vars, pt, n);
            if (!bestv) { bestv = val; bestpt = pt; kept = true; }
            else { int c = rru_sign_compare(val, bestv);
                if (c == -2) { expr_free(val); aborted = true; }
                else if (c < 0) { expr_free(bestv); bestv = val;
                    for (size_t i = 0; i < n; i++) expr_free(bestpt[i]);
                    free(bestpt); bestpt = pt; kept = true; }
                else expr_free(val); }
        }
        if (!kept) { for (size_t i = 0; i < n; i++) expr_free(pt[i]); free(pt); }
        size_t k = 0;                                /* odometer increment */
        for (; k < n; k++) { if (++cur[k] <= hi[k]) break; cur[k] = lo[k]; }
        if (k == n) done = true;
    }
    free(lo); free(hi); free(cur);
    Expr* r = NULL;
    if (!aborted && bestv) { r = mz_result(bestv, vars, bestpt, n); bestv = NULL; }
    if (bestv) expr_free(bestv);
    if (bestpt) { for (size_t i = 0; i < n; i++) expr_free(bestpt[i]); free(bestpt); }
    return r;
}

/* Exact Integers-domain optimisation. Tries the Diophantine enumeration first
 * (equality regions), then a box enumeration bounded by the continuous
 * relaxation (inequality regions). An absent constraint (no bound), or an
 * unbounded / undecidable region — e.g. x^3+y^3+z^3==33, whose solutions are
 * enormous — declines (NULL), never guesses. All variables are taken integer;
 * mixed integer/continuous problems are out of scope. */
static Expr* mz_integer_opt(const Expr* f, const Expr* cons, Expr* const* vars,
                            size_t n, const char* head, int budget) {
    (void)head;
    if (!cons || mz_is_true(cons)) return NULL;     /* unconstrained: unbounded */
    if (!mz_poly_in_vars(f, vars, n)) return NULL;
    clock_t t0 = clock();
    Expr* r = mz_integer_solve(f, cons, vars, n, budget);
    if (r) return r;
    int rem = mz_rem_budget(t0, budget);            /* box shares the call budget */
    return (rem < 0) ? NULL : mz_integer_box(f, cons, vars, n, rem);
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

    /* TimeConstraint -> t (seconds) bounds every internal Solve/Reduce probe so
     * a CAD/Gröbner blow-up declines gracefully instead of hanging; -> Infinity
     * removes the bound. Default MZ_DEFAULT_TIMECONSTRAINT. */
    int budget = mz_option_int(&A[pos_end], argc - pos_end,
                               "TimeConstraint", MZ_DEFAULT_TIMECONSTRAINT);

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
        /* Tier 2a: if the objective or constraints carry radicals / fractional
         * powers, polynomialize by adjoining fresh real variables, solve the
         * enlarged polynomial problem, then strip the auxiliary variables from
         * the reported point. */
        MzRadCtx rc = {0};
        Expr* fpoly = mz_rad_walk(fobj, &rc);
        Expr* cpoly = cons ? mz_rad_walk(cons, &rc) : NULL;
        if (rc.naux > 0) {
            size_t na = (cpoly ? 1u : 0u) + rc.ndefs;
            Expr** parts = (Expr**)malloc(sizeof(Expr*) * (na ? na : 1));
            size_t pi = 0;
            if (cpoly) parts[pi++] = expr_copy(cpoly);
            for (size_t k = 0; k < rc.ndefs; k++) parts[pi++] = expr_copy(rc.defs[k]);
            Expr* cons2 = expr_new_function(expr_new_symbol(SYM_And), parts, pi);
            free(parts);
            size_t n2 = n + rc.naux;
            Expr** vars2 = (Expr**)malloc(sizeof(Expr*) * n2);
            for (size_t j = 0; j < n; j++)         vars2[j]     = expr_copy(vars[j]);
            for (size_t k = 0; k < rc.naux; k++)   vars2[n + k] = expr_copy(rc.aux[k]);
            result = mz_exact_poly(fpoly, cons2, vars2, n2, head, budget);
            if (result) result = mz_strip_aux_rules(result, rc.aux, rc.naux);
            if (result && is_max) result = mz_negate_value(result);
            expr_free(cons2);
            for (size_t j = 0; j < n2; j++) expr_free(vars2[j]);
            free(vars2);
        } else {
            bool unconstrained = (cons == NULL);
            if (unconstrained && n == 1) {
                result = mz_univar_piecewise(fobj, vars[0], head, budget);
                if (!result) result = mz_univar_poly(fobj, vars[0], head, budget);
            } else
                result = mz_exact_poly(fobj, cons, vars, n, head, budget);
            if (result && is_max) result = mz_negate_value(result);
        }
        expr_free(fpoly);
        if (cpoly) expr_free(cpoly);
        mz_rad_free(&rc);
    } else if (domname == SYM_Integers) {
        result = mz_integer_opt(fobj, cons, vars, n, head, budget);
        if (result && is_max) result = mz_negate_value(result);
    }
    /* Complexes with exact input: deferred — leave unevaluated. */

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
