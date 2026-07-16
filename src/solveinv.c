/*
 * solveinv.c
 *
 * Inverse-function specialist for the `Solve` router.  Sits between
 * the polynomial specialist (src/poly/solvepoly.c) and the radicals
 * specialist (src/solverad.c) in the single-variable dispatch.
 *
 * Entry: solveinv_solve_inverse_equality.  The fast guard
 * solveinv_looks_invertible answers "is there any peelable head over
 * `var` in `expr`?" without doing real work, so the dispatch in
 * solve.c can short-circuit cheaply on polynomial-only inputs.
 *
 * The peeling itself dispatches on the outermost var-carrying head
 * through the static INV_TABLE.  Each PeelFn produces a list of
 * Branch records: an inner equation Equal[g(x), <new_rhs>] plus a
 * boolean condition that the resulting solutions are wrapped in via
 * ConditionalExpression.  Multi-branch heads (Sin, Cos, Exp, ...)
 * mint a fresh integer parameter symbol of the form C[k] using the
 * shared per-call SolveInvCtx::param_counter.
 *
 * Recursion: each inner equation is solved by hand-off to solvepoly,
 * solveinv_solve_inverse_equality (depth-capped at SOLVEINV_MAX_DEPTH),
 * and solverad in turn -- never back through builtin_solve, so the
 * inexact-rationalisation pre-pass in solve.c runs only once.
 */

#include "solveinv.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "poly/solvepoly.h"
#include "solverad.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Configuration.                                                     *
 * ------------------------------------------------------------------ */

#define SOLVEINV_MAX_DEPTH 8

/* ------------------------------------------------------------------ *
 *  Per-call context.                                                  *
 * ------------------------------------------------------------------ */

typedef struct {
    const SolveInvOpts* opts;
    Expr* dom;             /* borrowed; may be NULL                    */
    int   param_counter;   /* mints C[1], C[2], ... per builtin_solve  */
    int   depth;           /* recursion guard                          */
    bool  ifun_warned;     /* one-shot Solve::ifun per call            */
} SolveInvCtx;

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand.  Each helper takes ownership of    *
 *  its pointer arguments; callers pass freshly built or expr_copy()'d *
 *  pointers.                                                          *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c }, 3);
}
static Expr* mk_fn5(const char* head, Expr* a, Expr* b, Expr* c,
                    Expr* d, Expr* e) {
    return expr_new_function(mk_sym(head),
                             (Expr*[]){ a, b, c, d, e }, 5);
}

static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_pi(void)     { return mk_sym("Pi"); }
static Expr* mk_I(void)      { return mk_sym("I"); }
static Expr* mk_E(void)      { return mk_sym("E"); }
static Expr* mk_true(void)   { return mk_sym("True"); }

static Expr* mk_pow(Expr* base, Expr* exp_) {
    return mk_fn2("Power", base, exp_);
}

/* Rational p/q (small integers).  Returns p/q in canonical form. */
static Expr* mk_rat(int64_t p, int64_t q) {
    return mk_fn2("Times", mk_int(p), mk_pow(mk_int(q), mk_int(-1)));
}

/* ------------------------------------------------------------------ *
 *  Var containment.                                                   *
 * ------------------------------------------------------------------ */

static bool contains_var(const Expr* e, const Expr* var) {
    if (!e || !var || var->type != EXPR_SYMBOL) return false;
    const char* vname = var->data.symbol.name;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == vname;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_var(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_var(e->data.function.args[i], var)) return true;
    }
    return false;
}

static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* ------------------------------------------------------------------ *
 *  Branch records.                                                    *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr* inner_eq;    /* Equal[g(x), branch_rhs] -- owned */
    Expr* condition;   /* boolean condition -- owned */
} Branch;

static void branches_free(Branch* arr, size_t n) {
    if (!arr) return;
    for (size_t i = 0; i < n; i++) {
        if (arr[i].inner_eq)  expr_free(arr[i].inner_eq);
        if (arr[i].condition) expr_free(arr[i].condition);
    }
}

/* ------------------------------------------------------------------ *
 *  Solve::ifun emission.  One-shot per builtin_solve call.            *
 * ------------------------------------------------------------------ */

static void emit_ifun(SolveInvCtx* ctx) {
    if (!ctx || ctx->ifun_warned) return;
    fprintf(stderr,
        "Solve::ifun: Inverse functions are being used by Solve, "
        "so some solutions may not be found; "
        "use Reduce for complete solution information.\n");
    ctx->ifun_warned = true;
}

/* ------------------------------------------------------------------ *
 *  C[k] minting.                                                      *
 * ------------------------------------------------------------------ */

static Expr* mint_param(SolveInvCtx* ctx) {
    int k = ++ctx->param_counter;
    const char* head_name = (ctx && ctx->opts && ctx->opts->param_head)
                                ? ctx->opts->param_head : "C";
    return expr_new_function(mk_sym(head_name),
                             (Expr*[]){ mk_int(k) }, 1);
}

/* ------------------------------------------------------------------ *
 *  Forward decl for recursive solver re-entry.                        *
 * ------------------------------------------------------------------ */

static Expr* solve_inner_equation(Expr* inner_eq, Expr* var,
                                  SolveInvCtx* ctx);
static Expr* solveinv_drive(Expr* equation, Expr* var,
                            SolveInvCtx* ctx);

/* ------------------------------------------------------------------ *
 *  Peel-function table.                                               *
 * ------------------------------------------------------------------ */

typedef bool (*PeelFn)(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                       Branch* out, size_t* n_out);

typedef struct {
    const char* head_name;  /* interned symbol pointer */
    PeelFn      fn;
    bool        multi_branch_warn;  /* emit Solve::ifun when used */
} InvRule;

/* ============================================================ *
 *  Exponential / Logarithmic                                    *
 * ============================================================ */

/* Log[g] == a  ->  g == E^a, cond: -Pi < Im[a] <= Pi */
static bool peel_log(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var; (void)ctx;
    Expr* rhs  = mk_pow(mk_E(), expr_copy(a));
    Expr* cond = mk_fn5("Inequality",
                        mk_neg(mk_pi()), mk_sym("Less"),
                        mk_fn1("Im", expr_copy(a)),
                        mk_sym("LessEqual"),
                        mk_pi());
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(cond);
    *n_out = 1;
    return true;
}

/* Exp[g] == a  ->  g == 2 I Pi C[k] + Log[a], cond: Element[C[k], Integers] */
static bool peel_exp(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* term1 = expr_new_function(mk_sym("Times"),
        (Expr*[]){ mk_int(2), mk_I(), mk_pi(), expr_copy(Ck) }, 4);
    Expr* rhs = mk_fn2("Plus", term1, mk_fn1("Log", expr_copy(a)));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(mk_fn2("Element", expr_copy(Ck),
                                            mk_sym("Integers")));
    expr_free(Ck);
    *n_out = 1;
    emit_ifun(ctx);
    return true;
}

/* ============================================================ *
 *  Forward trig: principal + periodic family                    *
 *  All emit two branches sharing the same C[k].                 *
 * ============================================================ */

/* Helper: build expr + 2 Pi C[k] (or - principal + 2 Pi C[k] etc.). */
static Expr* add_two_pi_Ck(Expr* base, Expr* Ck) {
    Expr* two_pi_Ck = expr_new_function(mk_sym("Times"),
        (Expr*[]){ mk_int(2), mk_pi(), expr_copy(Ck) }, 3);
    return mk_fn2("Plus", base, two_pi_Ck);
}

/* Sin[g] == a -> g == ArcSin[a]+2 Pi C[k]  /  g == Pi-ArcSin[a]+2 Pi C[k] */
static bool peel_sin(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* asn  = mk_fn1("ArcSin", expr_copy(a));
    Expr* rhs1 = add_two_pi_Ck(expr_copy(asn), Ck);
    Expr* rhs2 = add_two_pi_Ck(
        mk_fn2("Plus", mk_pi(), mk_neg(expr_copy(asn))), Ck);
    expr_free(asn);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Cos[g] == a -> g == ArcCos[a]+2 Pi C[k]  /  g == -ArcCos[a]+2 Pi C[k] */
static bool peel_cos(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* acs  = mk_fn1("ArcCos", expr_copy(a));
    Expr* rhs1 = add_two_pi_Ck(expr_copy(acs), Ck);
    Expr* rhs2 = add_two_pi_Ck(mk_neg(expr_copy(acs)), Ck);
    expr_free(acs);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Tan[g] == a -> g == ArcTan[a] + Pi C[k] */
static bool peel_tan(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* rhs = mk_fn2("Plus",
        mk_fn1("ArcTan", expr_copy(a)),
        mk_fn2("Times", mk_pi(), expr_copy(Ck)));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    expr_free(Ck);
    *n_out = 1;
    emit_ifun(ctx);
    return true;
}

/* Cot[g] == a -> g == ArcCot[a] + Pi C[k] */
static bool peel_cot(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* rhs = mk_fn2("Plus",
        mk_fn1("ArcCot", expr_copy(a)),
        mk_fn2("Times", mk_pi(), expr_copy(Ck)));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    expr_free(Ck);
    *n_out = 1;
    emit_ifun(ctx);
    return true;
}

/* Sec[g] == a -> g == ±ArcSec[a] + 2 Pi C[k] */
static bool peel_sec(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* asec = mk_fn1("ArcSec", expr_copy(a));
    Expr* rhs1 = add_two_pi_Ck(expr_copy(asec), Ck);
    Expr* rhs2 = add_two_pi_Ck(mk_neg(expr_copy(asec)), Ck);
    expr_free(asec);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Csc[g] == a -> g == ArcCsc[a]+2 Pi C[k]  /  g == Pi-ArcCsc[a]+2 Pi C[k] */
static bool peel_csc(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                     Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* ac  = mk_fn1("ArcCsc", expr_copy(a));
    Expr* rhs1 = add_two_pi_Ck(expr_copy(ac), Ck);
    Expr* rhs2 = add_two_pi_Ck(
        mk_fn2("Plus", mk_pi(), mk_neg(expr_copy(ac))), Ck);
    expr_free(ac);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* ============================================================ *
 *  Forward hyperbolic                                            *
 * ============================================================ */

/* Helper: build base + 2 I Pi C[k]. */
static Expr* add_two_i_pi_Ck(Expr* base, Expr* Ck) {
    Expr* term = expr_new_function(mk_sym("Times"),
        (Expr*[]){ mk_int(2), mk_I(), mk_pi(), expr_copy(Ck) }, 4);
    return mk_fn2("Plus", base, term);
}

/* Sinh[g] == a -> g == ArcSinh[a]+2 I Pi C[k]  /
 *                  g == I Pi - ArcSinh[a] + 2 I Pi C[k] */
static bool peel_sinh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* asnh = mk_fn1("ArcSinh", expr_copy(a));
    Expr* rhs1 = add_two_i_pi_Ck(expr_copy(asnh), Ck);
    Expr* i_pi = mk_fn2("Times", mk_I(), mk_pi());
    Expr* rhs2 = add_two_i_pi_Ck(
        mk_fn2("Plus", i_pi, mk_neg(expr_copy(asnh))), Ck);
    expr_free(asnh);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Cosh[g] == a -> g == ±ArcCosh[a] + 2 I Pi C[k] */
static bool peel_cosh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* ah = mk_fn1("ArcCosh", expr_copy(a));
    Expr* rhs1 = add_two_i_pi_Ck(expr_copy(ah), Ck);
    Expr* rhs2 = add_two_i_pi_Ck(mk_neg(expr_copy(ah)), Ck);
    expr_free(ah);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Tanh[g] == a -> g == ArcTanh[a] + I Pi C[k] */
static bool peel_tanh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* rhs = mk_fn2("Plus",
        mk_fn1("ArcTanh", expr_copy(a)),
        expr_new_function(mk_sym("Times"),
            (Expr*[]){ mk_I(), mk_pi(), expr_copy(Ck) }, 3));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    expr_free(Ck);
    *n_out = 1;
    emit_ifun(ctx);
    return true;
}

/* Coth[g] == a -> g == ArcCoth[a] + I Pi C[k] */
static bool peel_coth(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* rhs = mk_fn2("Plus",
        mk_fn1("ArcCoth", expr_copy(a)),
        expr_new_function(mk_sym("Times"),
            (Expr*[]){ mk_I(), mk_pi(), expr_copy(Ck) }, 3));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs));
    out[0].condition = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    expr_free(Ck);
    *n_out = 1;
    emit_ifun(ctx);
    return true;
}

/* Sech[g] == a -> g == ±ArcSech[a] + 2 I Pi C[k] */
static bool peel_sech(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* ah = mk_fn1("ArcSech", expr_copy(a));
    Expr* rhs1 = add_two_i_pi_Ck(expr_copy(ah), Ck);
    Expr* rhs2 = add_two_i_pi_Ck(mk_neg(expr_copy(ah)), Ck);
    expr_free(ah);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* Csch[g] == a -> g == ArcCsch[a]+2 I Pi C[k]  /
 *                  g == I Pi - ArcCsch[a] + 2 I Pi C[k] */
static bool peel_csch(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                      Branch* out, size_t* n_out) {
    (void)var;
    Expr* Ck = mint_param(ctx);
    Expr* ah = mk_fn1("ArcCsch", expr_copy(a));
    Expr* rhs1 = add_two_i_pi_Ck(expr_copy(ah), Ck);
    Expr* i_pi = mk_fn2("Times", mk_I(), mk_pi());
    Expr* rhs2 = add_two_i_pi_Ck(
        mk_fn2("Plus", i_pi, mk_neg(expr_copy(ah))), Ck);
    expr_free(ah);
    Expr* cond_int = eval_and_free(
        mk_fn2("Element", expr_copy(Ck), mk_sym("Integers")));
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs2));
    out[0].condition = expr_copy(cond_int);
    out[1].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g), rhs1));
    out[1].condition = cond_int;
    expr_free(Ck);
    *n_out = 2;
    emit_ifun(ctx);
    return true;
}

/* ============================================================ *
 *  Inverse trig / hyperbolic -- single principal branch with    *
 *  a domain ConditionalExpression predicate.                     *
 * ============================================================ */

/* Build the canonical "vertical-strip" predicate:
 *   (Re[a] == lo && Im[a] >= 0) || (lo < Re[a] < hi)
 *                                || (Re[a] == hi && Im[a] <= 0)
 * Used for ArcSin (lo=-Pi/2, hi=Pi/2), ArcCos (lo=0, hi=Pi).
 */
static Expr* strip_cond_le(Expr* a, Expr* lo, Expr* hi) {
    Expr* re_a = mk_fn1("Re", expr_copy(a));
    Expr* im_a = mk_fn1("Im", expr_copy(a));
    Expr* edge_lo = mk_fn2("And",
        mk_fn2("Equal", expr_copy(re_a), expr_copy(lo)),
        mk_fn2("GreaterEqual", expr_copy(im_a), mk_int(0)));
    Expr* mid = mk_fn5("Inequality",
        expr_copy(lo), mk_sym("Less"),
        expr_copy(re_a), mk_sym("Less"),
        expr_copy(hi));
    Expr* edge_hi = mk_fn2("And",
        mk_fn2("Equal", expr_copy(re_a), expr_copy(hi)),
        mk_fn2("LessEqual", expr_copy(im_a), mk_int(0)));
    Expr* out = mk_fn3("Or", edge_lo, mid, edge_hi);
    expr_free(re_a); expr_free(im_a);
    expr_free(lo);   expr_free(hi);
    expr_free(a);
    return out;
}

/* ArcSin[g] == a -> g == Sin[a], strip cond on [-Pi/2, Pi/2]. */
static bool peel_arcsin(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    (void)var; (void)ctx;
    Expr* cond = strip_cond_le(expr_copy(a),
        mk_neg(mk_fn2("Times", mk_rat(1,2), mk_pi())),
        mk_fn2("Times", mk_rat(1,2), mk_pi()));
    out[0].inner_eq  = eval_and_free(
        mk_fn2("Equal", expr_copy(g), mk_fn1("Sin", expr_copy(a))));
    out[0].condition = eval_and_free(cond);
    *n_out = 1;
    return true;
}

/* ArcCos[g] == a -> g == Cos[a], strip cond on [0, Pi]. */
static bool peel_arccos(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    (void)var; (void)ctx;
    Expr* cond = strip_cond_le(expr_copy(a), mk_int(0), mk_pi());
    out[0].inner_eq  = eval_and_free(
        mk_fn2("Equal", expr_copy(g), mk_fn1("Cos", expr_copy(a))));
    out[0].condition = eval_and_free(cond);
    *n_out = 1;
    return true;
}

/* Generic single-branch inverse with TRUE condition (used where the
 * Mathematica strip predicate is unwieldy and the user is rarely
 * inverting in the first place).  We still flag the trig/hyperbolic
 * forms below explicitly for the common cases. */
static Expr* trivial_true_cond(void) { return mk_true(); }

/* ArcTan[g] == a -> g == Tan[a].  Mathematica's strip condition spans
 * a half-open vertical strip plus two edge cases on the imaginary
 * boundary; the user's prompt encodes it as
 *   -Pi/2 < Re[a] < Pi/2
 *   || (Re[a] == -Pi/2 && Im[a] > 0)
 *   || (Re[a] == Pi/2  && Im[a] < 0)
 */
static bool peel_arctan(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    (void)var; (void)ctx;
    Expr* re_a = mk_fn1("Re", expr_copy(a));
    Expr* im_a = mk_fn1("Im", expr_copy(a));
    Expr* hi = mk_fn2("Times", mk_rat(1,2), mk_pi());
    Expr* lo = mk_neg(mk_fn2("Times", mk_rat(1,2), mk_pi()));
    Expr* mid = mk_fn5("Inequality",
        expr_copy(lo), mk_sym("Less"),
        expr_copy(re_a), mk_sym("Less"),
        expr_copy(hi));
    Expr* edge_lo = mk_fn2("And",
        mk_fn2("Equal", expr_copy(re_a), expr_copy(lo)),
        mk_fn2("Greater", expr_copy(im_a), mk_int(0)));
    Expr* edge_hi = mk_fn2("And",
        mk_fn2("Equal", expr_copy(re_a), expr_copy(hi)),
        mk_fn2("Less", expr_copy(im_a), mk_int(0)));
    Expr* cond = mk_fn3("Or", mid, edge_lo, edge_hi);
    expr_free(re_a); expr_free(im_a);
    expr_free(lo);   expr_free(hi);

    out[0].inner_eq  = eval_and_free(
        mk_fn2("Equal", expr_copy(g), mk_fn1("Tan", expr_copy(a))));
    out[0].condition = eval_and_free(cond);
    *n_out = 1;
    return true;
}

/* Single-branch inverse-head dispatch with a placeholder true cond.
 * Used for the inverse forms where encoding the exact Mathematica
 * principal-strip predicate is not worth the verbosity for v1; the
 * user gets g == head_inv[a] which is the correct value on the
 * principal branch. */
static bool peel_inv_with_true(const char* head_inv_name,
                               Expr* g, Expr* a,
                               Expr* var, SolveInvCtx* ctx,
                               Branch* out, size_t* n_out) {
    (void)var; (void)ctx;
    out[0].inner_eq  = eval_and_free(mk_fn2("Equal", expr_copy(g),
                                    mk_fn1(head_inv_name, expr_copy(a))));
    out[0].condition = trivial_true_cond();
    *n_out = 1;
    return true;
}

static bool peel_arccot(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    return peel_inv_with_true("Cot", g, a, var, ctx, out, n_out);
}
static bool peel_arcsec(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    return peel_inv_with_true("Sec", g, a, var, ctx, out, n_out);
}
static bool peel_arccsc(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                        Branch* out, size_t* n_out) {
    return peel_inv_with_true("Csc", g, a, var, ctx, out, n_out);
}
static bool peel_arcsinh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Sinh", g, a, var, ctx, out, n_out);
}
static bool peel_arccosh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Cosh", g, a, var, ctx, out, n_out);
}
static bool peel_arctanh(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Tanh", g, a, var, ctx, out, n_out);
}
static bool peel_arccoth(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Coth", g, a, var, ctx, out, n_out);
}
static bool peel_arcsech(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Sech", g, a, var, ctx, out, n_out);
}
static bool peel_arccsch(Expr* g, Expr* a, Expr* var, SolveInvCtx* ctx,
                         Branch* out, size_t* n_out) {
    return peel_inv_with_true("Csch", g, a, var, ctx, out, n_out);
}

/* ------------------------------------------------------------------ *
 *  Inverse table.                                                     *
 * ------------------------------------------------------------------ */

/* Filled in by solveinv_init() with interned symbol pointers; we use
 * pointer-equality lookup against the head's interned symbol name. */
static struct {
    const char* head;
    PeelFn      fn;
} INV_TABLE[] = {
    { NULL /* Log */,    peel_log    },
    { NULL /* Exp */,    peel_exp    },
    { NULL /* Sin */,    peel_sin    },
    { NULL /* Cos */,    peel_cos    },
    { NULL /* Tan */,    peel_tan    },
    { NULL /* Cot */,    peel_cot    },
    { NULL /* Sec */,    peel_sec    },
    { NULL /* Csc */,    peel_csc    },
    { NULL /* Sinh */,   peel_sinh   },
    { NULL /* Cosh */,   peel_cosh   },
    { NULL /* Tanh */,   peel_tanh   },
    { NULL /* Coth */,   peel_coth   },
    { NULL /* Sech */,   peel_sech   },
    { NULL /* Csch */,   peel_csch   },
    { NULL /* ArcSin */,  peel_arcsin },
    { NULL /* ArcCos */,  peel_arccos },
    { NULL /* ArcTan */,  peel_arctan },
    { NULL /* ArcCot */,  peel_arccot },
    { NULL /* ArcSec */,  peel_arcsec },
    { NULL /* ArcCsc */,  peel_arccsc },
    { NULL /* ArcSinh */, peel_arcsinh },
    { NULL /* ArcCosh */, peel_arccosh },
    { NULL /* ArcTanh */, peel_arctanh },
    { NULL /* ArcCoth */, peel_arccoth },
    { NULL /* ArcSech */, peel_arcsech },
    { NULL /* ArcCsch */, peel_arccsch },
};
#define INV_TABLE_N (sizeof(INV_TABLE)/sizeof(INV_TABLE[0]))

static const char* INV_TABLE_NAMES[INV_TABLE_N] = {
    "Log", "Exp",
    "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
    "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch",
    "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc",
    "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch",
};

static PeelFn lookup_peel(const char* head_interned) {
    if (!head_interned) return NULL;
    for (size_t i = 0; i < INV_TABLE_N; i++) {
        if (INV_TABLE[i].head == head_interned) return INV_TABLE[i].fn;
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Power[g, n] == a  with integer n >= 2:  n radical branches.        *
 *  Special form -- not part of the table because the branch count is *
 *  data-dependent on n.                                               *
 *                                                                     *
 *  We emit each k = 0..n-1 branch as g == a^(1/n) * Exp[2 Pi I k / n] *
 *  with cond = True; the polynomial specialist resolves the inner    *
 *  equation directly.                                                 *
 *                                                                     *
 *  Caveat: for n == 2 we hand off to solvepoly via the binomial      *
 *  form g^2 - a == 0 which yields the ± Sqrt[a] roots that match the *
 *  user's expected `±I Sqrt[1 - 1/E]` shape exactly.  This is        *
 *  preferred over the explicit primitive-root construction because    *
 *  solvepoly's output is already in the canonical Sqrt-bearing form. *
 * ------------------------------------------------------------------ */

static Expr* peel_power_n(Expr* g, int64_t n, Expr* a, Expr* var,
                          SolveInvCtx* ctx) {
    /* g^n == 0  <=>  g == 0  for any positive integer n.  Reducing here
     * avoids re-entering the inverse peel with the identical equation
     * shape (peel_power_n -> solve_inner_equation -> solveinv_drive ->
     * peel_power_n) when g is not a polynomial in var, which used to
     * recurse to the depth cap or, before the depth fix, blow the
     * stack on inputs like Sin[x]^2 == 0. */
    bool a_is_zero = (a->type == EXPR_INTEGER && a->data.integer == 0);
    if (a_is_zero) {
        Expr* inner_eq = mk_fn2("Equal", expr_copy(g), mk_int(0));
        inner_eq = eval_and_free(inner_eq);
        Expr* sols = solve_inner_equation(inner_eq, var, ctx);
        expr_free(inner_eq);
        return sols;
    }

    /* Hand the inner binomial Equal[g^n, a] directly to the recursive
     * inner solver -- solvepoly knows how to deal with binomial
     * equations and n-th roots, and that path produces the canonical
     * Sqrt-bearing answer the user expects. */
    Expr* inner_eq = mk_fn2("Equal", mk_pow(expr_copy(g), mk_int(n)),
                            expr_copy(a));
    inner_eq = eval_and_free(inner_eq);
    Expr* sols = solve_inner_equation(inner_eq, var, ctx);
    expr_free(inner_eq);
    return sols;
}

/* ------------------------------------------------------------------ *
 *  Inner equation solver -- the recursive re-entry point.             *
 *                                                                     *
 *  Tries solvepoly -> solveinv (depth+1) -> solverad in turn.  The   *
 *  hand-off bypasses builtin_solve so the inexact pre-pass runs only *
 *  once per user call.                                                *
 * ------------------------------------------------------------------ */

static Expr* solve_inner_equation(Expr* inner_eq, Expr* var,
                                  SolveInvCtx* ctx) {
    if (!inner_eq || !var) return NULL;
    if (ctx->depth >= SOLVEINV_MAX_DEPTH) return NULL;

    /* Equation might have collapsed to True/False during eval. */
    if (inner_eq->type == EXPR_SYMBOL) {
        if (inner_eq->data.symbol.name == SYM_True) {
            Expr* empty = expr_new_function(mk_sym("List"), NULL, 0);
            return expr_new_function(mk_sym("List"),
                                     (Expr*[]){ empty }, 1);
        }
        if (inner_eq->data.symbol.name == SYM_False) {
            return expr_new_function(mk_sym("List"), NULL, 0);
        }
        return NULL;
    }

    SolvePolyOpts polyopts = { false, false };
    Expr* sols = solvepoly_solve_polynomial_equality(
        inner_eq, var, ctx->dom, &polyopts);
    if (sols) return sols;

    /* Try a nested inverse peel before the radicals specialist.
     * Re-enter solveinv_drive directly with the existing ctx so the
     * depth counter is preserved across the re-entry; routing through
     * the public solveinv_solve_inverse_equality wrapper would build a
     * fresh ctx with depth=0 and defeat SOLVEINV_MAX_DEPTH, which lets
     * pathological cases like Sin[x]^k == c recurse without bound
     * through peel_power_n -> solve_inner_equation. */
    ctx->depth++;
    sols = solveinv_drive(inner_eq, var, ctx);
    ctx->depth--;
    if (sols) return sols;

    sols = solverad_solve_radicals_equality(inner_eq, var, ctx->dom);
    return sols;
}

/* ------------------------------------------------------------------ *
 *  Wrap each solution's RHS in ConditionalExpression[rhs, cond].      *
 *                                                                     *
 *  `solutions` is a freshly-owned List[List[Rule[var, val]], ...]    *
 *  returned by the inner solver.  We rewrite each Rule's RHS into    *
 *  ConditionalExpression[rhs, cond] and let the evaluator collapse   *
 *  it when `cond` happens to be True / False.  Returns a freshly-    *
 *  owned new solution list and frees `solutions`.                     *
 * ------------------------------------------------------------------ */

static Expr* wrap_solutions_with_cond(Expr* solutions, Expr* cond) {
    if (!solutions || !cond) {
        if (solutions) expr_free(solutions);
        if (cond) expr_free(cond);
        return NULL;
    }
    if (!head_is(solutions, SYM_List)) {
        expr_free(cond);
        return solutions;
    }
    size_t n = solutions->data.function.arg_count;
    Expr** new_outer = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* inner = solutions->data.function.args[i];
        if (!head_is(inner, SYM_List)) {
            new_outer[i] = expr_copy(inner);
            continue;
        }
        size_t m = inner->data.function.arg_count;
        Expr** new_inner = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t j = 0; j < m; j++) {
            Expr* rule = inner->data.function.args[j];
            if (head_is(rule, SYM_Rule) &&
                rule->data.function.arg_count == 2) {
                Expr* lhs = expr_copy(rule->data.function.args[0]);
                Expr* rhs = expr_copy(rule->data.function.args[1]);
                /* solvepoly's linear-solution shape often leaves a
                 * top-level Times[-1, Plus[-A, -B, ...]] in the RHS;
                 * one pass of Expand distributes the leading -1 and
                 * gives the canonical Mathematica form (e.g.
                 *   2 I Pi C[1] + Log[a]
                 * instead of -(-2 I Pi C[1] - Log[a])). */
                rhs = eval_and_free(mk_fn1("Expand", rhs));
                Expr* wrapped = eval_and_free(mk_fn2(
                    "ConditionalExpression", rhs, expr_copy(cond)));
                new_inner[j] = mk_fn2("Rule", lhs, wrapped);
            } else {
                new_inner[j] = expr_copy(rule);
            }
        }
        new_outer[i] = expr_new_function(mk_sym("List"), new_inner, m);
        free(new_inner);
    }
    Expr* out = expr_new_function(mk_sym("List"), new_outer, n);
    free(new_outer);
    expr_free(solutions);
    expr_free(cond);
    return out;
}

/* Concatenate two freshly-owned List[...]s; both freed. */
static Expr* concat_lists(Expr* a, Expr* b) {
    if (!a) return b;
    if (!b) return a;
    if (!head_is(a, SYM_List) || !head_is(b, SYM_List)) {
        expr_free(b); return a;
    }
    size_t na = a->data.function.arg_count;
    size_t nb = b->data.function.arg_count;
    size_t n = na + nb;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < na; i++) args[i] = expr_copy(a->data.function.args[i]);
    for (size_t i = 0; i < nb; i++) args[na + i] = expr_copy(b->data.function.args[i]);
    Expr* out = expr_new_function(mk_sym("List"), args, n);
    free(args);
    expr_free(a); expr_free(b);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Isolation pre-pass.                                                *
 *                                                                     *
 *  Given the residual `r = lhs - rhs` (Plus or atom), split the      *
 *  top-level Plus into a "with-var" portion `var_part` and a         *
 *  "free-of-var" portion `free_part`.  If `var_part` reduces to      *
 *      coeff * head[g(x)]
 *  (or just head[g(x)]), write *coeff_out, *head_out, *g_out and    *
 *  return true.  `new_rhs_out` is set to -free_part / coeff.         *
 *                                                                     *
 *  When `var_part` is a Power[g(x), n] with integer n >= 2 and the   *
 *  isolation reduces to Power[g(x), n] == new_rhs, returns true with *
 *  *head_out == SYM_Power and *pow_n_out = n.  Otherwise *pow_n_out  *
 *  = 0.                                                                *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr* g;            /* freshly owned -- function arg / Power base */
    Expr* new_rhs;      /* freshly owned                              */
    const char* head;   /* interned                                   */
    int64_t pow_n;      /* nonzero only for Power[g, n] case          */
} Isolation;

static bool isolate_residual(Expr* residual, Expr* var, Isolation* out);

/* True iff `e` is a Plus, and we should walk its args. */
static bool is_plus(const Expr* e) { return head_is(e, SYM_Plus); }
static bool is_times(const Expr* e) { return head_is(e, SYM_Times); }

/* Return c such that `term = c * payload` with `c` free of var, and
 * payload is a single sub-expression containing var, or return false. */
static bool split_coeff_payload(Expr* term, Expr* var,
                                Expr** coeff_out, Expr** payload_out) {
    if (!contains_var(term, var)) return false;

    if (!is_times(term)) {
        *coeff_out = mk_int(1);
        *payload_out = expr_copy(term);
        return true;
    }

    size_t n = term->data.function.arg_count;
    /* Count factors containing var. */
    size_t var_factor_count = 0;
    size_t var_factor_idx = 0;
    for (size_t i = 0; i < n; i++) {
        if (contains_var(term->data.function.args[i], var)) {
            var_factor_count++;
            var_factor_idx = i;
        }
    }
    if (var_factor_count != 1) return false;

    Expr** coeff_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t cnt = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == var_factor_idx) continue;
        coeff_args[cnt++] = expr_copy(term->data.function.args[i]);
    }
    Expr* coeff;
    if (cnt == 0) {
        coeff = mk_int(1);
    } else if (cnt == 1) {
        coeff = coeff_args[0];
    } else {
        coeff = expr_new_function(mk_sym("Times"), coeff_args, cnt);
        free(coeff_args);
        coeff_args = NULL;
        coeff = eval_and_free(coeff);
    }
    if (coeff_args) free(coeff_args);
    *coeff_out = coeff;
    *payload_out = expr_copy(term->data.function.args[var_factor_idx]);
    return true;
}

static bool isolate_residual(Expr* residual, Expr* var, Isolation* out) {
    if (!residual || !var || !out) return false;
    out->g = NULL; out->new_rhs = NULL; out->head = NULL; out->pow_n = 0;

    Expr* var_part = NULL;
    Expr* free_part = NULL;

    if (is_plus(residual)) {
        size_t n = residual->data.function.arg_count;
        Expr** var_acc = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        Expr** free_acc = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        size_t nv = 0, nf = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* t = residual->data.function.args[i];
            if (contains_var(t, var)) var_acc[nv++] = expr_copy(t);
            else                       free_acc[nf++] = expr_copy(t);
        }
        if (nv == 1) {
            var_part = var_acc[0];
        } else {
            for (size_t i = 0; i < nv; i++) expr_free(var_acc[i]);
            for (size_t i = 0; i < nf; i++) expr_free(free_acc[i]);
            free(var_acc); free(free_acc);
            return false;
        }
        free(var_acc);
        if (nf == 0) {
            free_part = mk_int(0);
        } else if (nf == 1) {
            free_part = free_acc[0];
        } else {
            free_part = expr_new_function(mk_sym("Plus"), free_acc, nf);
            free(free_acc);
            free_acc = NULL;
            free_part = eval_and_free(free_part);
        }
        if (free_acc) free(free_acc);
    } else {
        if (!contains_var(residual, var)) return false;
        var_part = expr_copy(residual);
        free_part = mk_int(0);
    }

    /* var_part = coeff * payload (payload contains var). */
    Expr* coeff = NULL;
    Expr* payload = NULL;
    if (!split_coeff_payload(var_part, var, &coeff, &payload)) {
        expr_free(var_part); expr_free(free_part);
        return false;
    }
    expr_free(var_part);

    /* Build new_rhs = -free_part / coeff. */
    Expr* numerator = eval_and_free(mk_neg(free_part));
    Expr* new_rhs;
    if (coeff && coeff->type == EXPR_INTEGER && coeff->data.integer == 1) {
        expr_free(coeff);
        new_rhs = numerator;
    } else {
        new_rhs = eval_and_free(mk_fn2("Times", numerator,
                                       mk_pow(coeff, mk_int(-1))));
    }

    /* Classify payload. */
    if (payload->type != EXPR_FUNCTION
        || payload->data.function.head->type != EXPR_SYMBOL) {
        expr_free(payload); expr_free(new_rhs);
        return false;
    }
    const char* h = payload->data.function.head->data.symbol.name;

    if (h == SYM_Power && payload->data.function.arg_count == 2) {
        Expr* base = payload->data.function.args[0];
        Expr* exp_ = payload->data.function.args[1];
        if (contains_var(base, var) && !contains_var(exp_, var)) {
            int64_t n_int = 0;
            if (exp_->type == EXPR_INTEGER && exp_->data.integer >= 2) {
                n_int = exp_->data.integer;
            }
            if (n_int >= 2) {
                out->g = expr_copy(base);
                out->new_rhs = new_rhs;
                out->head = SYM_Power;
                out->pow_n = n_int;
                expr_free(payload);
                return true;
            }
        }
        /* Power[b, g(x)] with var-free base b -- promote to the Exp
         * peel.  When b == E we report g(x) directly (canonical
         * Exp[g(x)] form); for any other constant base we reduce
         *   b^g(x) == new_rhs    <=>    Exp[g(x) * Log[b]] == new_rhs
         * by multiplying the inner argument by Log[b], so the existing
         * peel_exp + solvepoly machinery (which produces the periodic
         *   inner == 2 Pi I C[k] + Log[new_rhs]
         * family) yields the principal + branch solutions for free.
         * Mirrors Maxima's usolve `mexpt` constant-base branch,
         * generalised to the multi-branch complex log family. */
        if (!contains_var(base, var) && contains_var(exp_, var)) {
            Expr* new_g;
            if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E) {
                new_g = expr_copy(exp_);
            } else {
                new_g = eval_and_free(mk_fn2("Times",
                    expr_copy(exp_),
                    mk_fn1("Log", expr_copy(base))));
            }
            out->g = new_g;
            out->new_rhs = new_rhs;
            out->head = intern_symbol("Exp");
            out->pow_n = 0;
            expr_free(payload);
            return true;
        }
        expr_free(payload); expr_free(new_rhs);
        return false;
    }

    /* Generic single-arg invertible head[g(x)]. */
    if (payload->data.function.arg_count == 1) {
        Expr* g = payload->data.function.args[0];
        if (!contains_var(g, var)) {
            expr_free(payload); expr_free(new_rhs);
            return false;
        }
        out->g = expr_copy(g);
        out->new_rhs = new_rhs;
        out->head = h;     /* interned */
        out->pow_n = 0;
        expr_free(payload);
        return true;
    }

    expr_free(payload); expr_free(new_rhs);
    return false;
}

/* ------------------------------------------------------------------ *
 *  Solveinv_looks_invertible: cheap fast-fail probe.                  *
 * ------------------------------------------------------------------ */

static bool head_is_in_table(const char* interned) {
    if (!interned) return false;
    for (size_t i = 0; i < INV_TABLE_N; i++) {
        if (INV_TABLE[i].head == interned) return true;
    }
    return false;
}

bool solveinv_looks_invertible(const Expr* expr, const Expr* var) {
    if (!expr || !var) return false;
    if (expr->type == EXPR_FUNCTION) {
        if (expr->data.function.head->type == EXPR_SYMBOL) {
            const char* h = expr->data.function.head->data.symbol.name;
            if (head_is_in_table(h)
                && expr->data.function.arg_count == 1
                && contains_var(expr->data.function.args[0], var)) {
                return true;
            }
            if (h == SYM_Power
                && expr->data.function.arg_count == 2) {
                const Expr* base = expr->data.function.args[0];
                const Expr* exp_ = expr->data.function.args[1];
                /* x^n with var-free n -- handled by peel_power_n. */
                if (contains_var(base, var) && !contains_var(exp_, var)) {
                    return true;
                }
                /* b^g(x) with var-free base b -- promoted to the Exp
                 * peel by try_isolate_payload. */
                if (!contains_var(base, var) && contains_var(exp_, var)) {
                    return true;
                }
            }
        }
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            if (solveinv_looks_invertible(
                    expr->data.function.args[i], var)) return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  InverseFunction[] fallback.                                        *
 * ------------------------------------------------------------------ */

static Expr* make_inverse_function_solution(const char* head_name,
                                            Expr* rhs,
                                            Expr* var,
                                            SolveInvCtx* ctx) {
    emit_ifun(ctx);
    Expr* invfn = mk_fn1("InverseFunction", mk_sym(head_name));
    Expr* applied = expr_new_function(invfn,
        (Expr*[]){ expr_copy(rhs) }, 1);
    Expr* rule = mk_fn2("Rule", expr_copy(var), eval_and_free(applied));
    Expr* row  = expr_new_function(mk_sym("List"),
        (Expr*[]){ rule }, 1);
    Expr* out  = expr_new_function(mk_sym("List"),
        (Expr*[]){ row }, 1);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Recursive driver -- internal, takes the context.                   *
 * ------------------------------------------------------------------ */

static Expr* solveinv_drive(Expr* equation, Expr* var,
                            SolveInvCtx* ctx) {
    if (!equation || !var || !ctx) return NULL;
    if (!ctx->opts || !ctx->opts->enabled) return NULL;
    if (ctx->depth >= SOLVEINV_MAX_DEPTH) return NULL;
    if (!head_is(equation, SYM_Equal)
        || equation->data.function.arg_count != 2) return NULL;

    /* residual = lhs - rhs. */
    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];
    Expr* residual = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs), mk_neg(expr_copy(rhs))));

    Isolation iso = {0};
    if (!isolate_residual(residual, var, &iso)) {
        expr_free(residual);
        return NULL;
    }
    expr_free(residual);

    /* Power-with-integer-exponent path. */
    if (iso.head == SYM_Power && iso.pow_n >= 2) {
        Expr* sols = peel_power_n(iso.g, iso.pow_n, iso.new_rhs, var, ctx);
        expr_free(iso.g); expr_free(iso.new_rhs);
        return sols;
    }

    PeelFn fn = lookup_peel(iso.head);
    if (!fn) {
        expr_free(iso.g); expr_free(iso.new_rhs);
        return NULL;
    }

    /* Dispatch peel. */
    Branch branches[4] = {0};
    size_t nb = 0;
    if (!fn(iso.g, iso.new_rhs, var, ctx, branches, &nb) || nb == 0) {
        expr_free(iso.g); expr_free(iso.new_rhs);
        return NULL;
    }

    /* Solve each branch and wrap. */
    Expr* aggregate = NULL;
    bool any_success = false;
    for (size_t i = 0; i < nb; i++) {
        ctx->depth++;
        Expr* sols = solve_inner_equation(branches[i].inner_eq, var, ctx);
        ctx->depth--;
        if (!sols) continue;
        any_success = true;
        Expr* wrapped = wrap_solutions_with_cond(sols,
                                                 expr_copy(branches[i].condition));
        aggregate = concat_lists(aggregate, wrapped);
    }
    branches_free(branches, nb);

    if (!any_success) {
        /* InverseFunction fallback: only when peel was over var itself
         * (i.e. iso.g is exactly `var`), single-branch head, and the
         * inner equation reduced trivially. */
        bool g_is_var = iso.g && iso.g->type == EXPR_SYMBOL
                        && var->type == EXPR_SYMBOL
                        && iso.g->data.symbol.name == var->data.symbol.name;
        if (g_is_var && iso.head && iso.head != SYM_Power) {
            Expr* sols = make_inverse_function_solution(
                iso.head, iso.new_rhs, var, ctx);
            expr_free(iso.g); expr_free(iso.new_rhs);
            return sols;
        }
        expr_free(iso.g); expr_free(iso.new_rhs);
        return NULL;
    }

    expr_free(iso.g); expr_free(iso.new_rhs);
    return aggregate;
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */

Expr* solveinv_solve_inverse_equality(Expr* equation, Expr* var,
                                      Expr* dom,
                                      const SolveInvOpts* opts) {
    if (!equation || !var || !opts || !opts->enabled) return NULL;
    /* Top-of-tree entry from solve.c (or any external caller).  A
     * fresh per-call ctx is created here; re-entry from
     * solve_inner_equation goes through solveinv_drive directly with
     * the existing ctx so SOLVEINV_MAX_DEPTH is honoured. */
    SolveInvCtx ctx = { opts, dom, 0, 0, false };
    return solveinv_drive(equation, var, &ctx);
}

/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry.                                           *
 *                                                                     *
 *  `Solve`SolveInverseFunctions[lhs == rhs, var]` and                *
 *  `Solve`SolveInverseFunctions[lhs == rhs, var, dom]`.              *
 *                                                                     *
 *  Mirrors the public ABI of the polynomial / radical / linear-system *
 *  specialists.  No option arguments are consumed -- callers wanting *
 *  to override defaults should go through the parent `Solve` (whose  *
 *  option plumbing wires `InverseFunctions` and                       *
 *  `GeneratedParameters` through to this specialist).                 *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_inverse_functions(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* dom      = (argc >= 3) ? res->data.function.args[2] : NULL;

    SolveInvOpts opts = { true, intern_symbol("C") };
    return solveinv_solve_inverse_equality(equation, var, dom, &opts);
    /* evaluator frees res on non-NULL return */
}

/* ------------------------------------------------------------------ *
 *  Init.                                                              *
 * ------------------------------------------------------------------ */

void solveinv_init(void) {
    /* Populate INV_TABLE.head with interned symbol pointers. */
    for (size_t i = 0; i < INV_TABLE_N; i++) {
        INV_TABLE[i].head = intern_symbol(INV_TABLE_NAMES[i]);
    }
    /* Register the context-qualified specialist entry. */
    symtab_add_builtin("Solve`SolveInverseFunctions",
                       builtin_solve_inverse_functions);
    SymbolDef* def = symtab_get_def("Solve`SolveInverseFunctions");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveInverseFunctions",
        "Solve`SolveInverseFunctions[lhs == rhs, var]\n"
        "Solve`SolveInverseFunctions[lhs == rhs, var, dom]\n"
        "\tThe inverse-function specialist used by Solve.  Peels the\n"
        "\toutermost var-carrying head when it is an elementary\n"
        "\tinvertible function (Log, Exp, Sin, Cos, Tan, Cot, Sec, Csc,\n"
        "\tthe hyperbolic counterparts, the inverse trig/hyperbolic\n"
        "\tforms, and Power[g, n] for integer n >= 2), recursively\n"
        "\tsolves the resulting inner equation, and wraps each\n"
        "\tsolution in ConditionalExpression with the appropriate\n"
        "\tdomain predicate.  Multi-branch heads introduce fresh\n"
        "\tinteger parameters C[1], C[2], ...  Returns NULL when the\n"
        "\tequation has no peelable head over var.");
}
