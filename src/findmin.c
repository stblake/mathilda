/*
 * findmin.c
 *
 * FindMinimum / FindMaximum — Mathematica-compatible local numerical
 * optimization. Both have HoldAll | Protected attributes and use a
 * Block-style snapshot/restore of the search variables' OwnValues so
 * that user-level definitions of those names are not perturbed during
 * iteration.
 *
 * Supported forms
 * ---------------
 *   FindMinimum[f,           {x, x0}]                  1D, Brent default
 *   FindMinimum[f,           {x, x0, x1}]              1D, two-start bracket
 *   FindMinimum[f,           {x, xstart, xmin, xmax}]  1D, bracket
 *   FindMinimum[f,           {{x, x0}, {y, y0}, ...}]  n-D, QuasiNewton default
 *   FindMinimum[f,           {x, y, ...}]              n-D, auto start = 0
 *   FindMinimum[{f, cons},   vars]                     constrained
 *
 * Options (Rule[...] in trailing position, any order):
 *   Method            -> Automatic | "Brent" | "Newton" | "QuasiNewton"
 *                                  | "ConjugateGradient"
 *   WorkingPrecision  -> MachinePrecision | digits   (MPFR for Brent + BFGS)
 *   MaxIterations     -> positive integer (default 500)
 *   AccuracyGoal      -> Automatic | Infinity | digits
 *   PrecisionGoal     -> Automatic | Infinity | digits
 *   Gradient          -> Automatic | { dfdx1, dfdx2, ... }
 *   StepMonitor       -> :> body
 *   EvaluationMonitor -> :> body
 *
 * Constraints (inside the {f, cons} form): boolean tree of comparisons.
 *   Box  ( a <= x <= b , x >= a , x <= b , etc. on a bare variable )
 *     → enforced by projection after each iterate.
 *   General ( g(x) <= 0 , h(x) == 0 , etc. )
 *     → quadratic-penalty wrapper around the inner solver; outer μ schedule.
 *   Or[...] / Element / Integers → emit FindMinimum::nimpl and return NULL.
 *
 * Output: { f_min, { x -> x_min, y -> y_min, ... } }.
 * FindMaximum returns { f_max, ... } via a thin wrapper that minimises −f
 * and negates the first component of the result.
 *
 * Returns NULL (unevaluated) on any failure — variable bindings are always
 * restored to their pre-call state, even on the error path.
 */

#include "findmin.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Types                                                              *
 * ------------------------------------------------------------------ */

typedef enum {
    FM_METHOD_AUTOMATIC = 0,
    FM_METHOD_BRENT,
    FM_METHOD_QUASINEWTON,   /* BFGS */
    FM_METHOD_CONJGRAD,      /* Polak-Ribière+ */
    FM_METHOD_NEWTON
} FmMethod;

typedef enum {
    FM_PREC_MACHINE = 0
#ifdef USE_MPFR
    , FM_PREC_MPFR
#endif
} FmPrecMode;

typedef enum {
    FM_SPEC_BAD = 0,
    FM_SPEC_VAR_ONLY,   /* {x}                       -> x0 = 0           */
    FM_SPEC_SINGLE,     /* {x, x0}                                       */
    FM_SPEC_TWO_START,  /* {x, x0, x1}               -> derivative free  */
    FM_SPEC_BRACKET     /* {x, xstart, xmin, xmax}                       */
} FmSpecKind;

typedef struct {
    FmMethod   method;
    FmPrecMode prec_mode;
    long       wp_bits;          /* MPFR bits when prec_mode == MPFR    */
    int64_t    max_iter;         /* default 500                          */
    double     acc_goal_digits;  /* filled at WP/2 if Automatic          */
    double     prec_goal_digits;
    Expr*      gradient;         /* borrowed; user-supplied list or NULL */
    Expr*      step_monitor;     /* borrowed; or NULL                    */
    Expr*      eval_monitor;     /* borrowed; or NULL                    */
} FmOpts;

/* Per-variable snapshot for Block-style binding. */
typedef struct {
    const char* name;            /* interned                            */
    Rule*       saved_own;
    uint32_t    saved_attrs;
    bool        valid;
} FmVarBind;

/* Per-variable box (used both for 4-elt specs and parsed inequalities).
 * `has_lo`/`has_hi` flag presence. */
typedef struct {
    bool   has_lo, has_hi;
    double lo, hi;
} FmBox;

/* General (non-box) inequality g(x) <= 0 or equality h(x) == 0. The
 * objective during the outer μ loop is f(x) + μ * Σ max(0,g_i)^2
 *                                              + μ * Σ h_j^2.
 * For the inner solver we evaluate each constraint expression directly. */
typedef struct {
    Expr*  expr;     /* owned: rewritten so feasible = (expr <= 0) or (expr == 0) */
    bool   equality; /* true → equality constraint                                */
} FmGenCon;

/* ------------------------------------------------------------------ *
 *  Diagnostic helper                                                  *
 * ------------------------------------------------------------------ */

static void fm_warn(const char* fn, const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "%s::%s: ", fn, tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* The driver function name is captured by the outer call (FindMinimum vs
 * FindMaximum) and threaded through so all diagnostics carry the right
 * tag. */
static const char* g_fm_name = "FindMinimum";

/* ------------------------------------------------------------------ *
 *  Numeric extraction / construction                                  *
 * ------------------------------------------------------------------ */

static bool fm_expr_to_double_real(Expr* e, double* out) {
    if (!e) return false;
    int64_t rn, rd;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;    return true;
        case EXPR_REAL:    *out = e->data.real;               return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);  return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION: {
            if (is_rational(e, &rn, &rd)) { *out = (double)rn / (double)rd; return true; }
            /* Tolerate Complex[a, eps] residues (numericalize occasionally
             * emits them when subtraction cancellation generates a tiny
             * imaginary part). */
            Expr* re_e;
            Expr* im_e;
            if (is_complex(e, &re_e, &im_e)) {
                double re_d, im_d;
                if (fm_expr_to_double_real(re_e, &re_d)
                    && fm_expr_to_double_real(im_e, &im_d)) {
                    double mag_re = fabs(re_d);
                    if (fabs(im_d) <= 1e-12 * (1.0 + mag_re)) {
                        *out = re_d;
                        return true;
                    }
                }
            }
            return false;
        }
        default: return false;
    }
}

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */

static bool fm_is_known_option_name(const char* s) {
    return s == SYM_Method
        || s == SYM_WorkingPrecision
        || s == SYM_MaxIterations
        || s == SYM_AccuracyGoal
        || s == SYM_PrecisionGoal
        || s == SYM_Gradient
        || s == SYM_StepMonitor
        || s == SYM_EvaluationMonitor;
}

static bool fm_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fm_is_known_option_name(lhs->data.symbol);
}

static bool fm_parse_working_precision(Expr* val,
                                       FmPrecMode* mode, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mode = FM_PREC_MACHINE; *bits = 0; return true;
    }
    double digits = 0.0;
    int64_t rn, rd;
    if (val->type == EXPR_INTEGER)         digits = (double)val->data.integer;
    else if (val->type == EXPR_REAL)       digits = val->data.real;
    else if (is_rational(val, &rn, &rd))   digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        *mode = FM_PREC_MACHINE; *bits = 0;
    } else {
        *mode = FM_PREC_MPFR;
        *bits = numeric_digits_to_bits(digits);
    }
    return true;
#else
    *mode = FM_PREC_MACHINE; *bits = 0; return true;
#endif
}

static bool fm_parse_goal(Expr* val, double* digits_out) {
    if (val->type == EXPR_SYMBOL) {
        if (val->data.symbol == SYM_Automatic) { *digits_out = -1.0; return true; }
        if (val->data.symbol == SYM_Infinity)  { *digits_out = INFINITY; return true; }
        return false;
    }
    return fm_expr_to_double_real(val, digits_out);
}

static bool fm_apply_option(Expr* rule, FmOpts* opts) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
            opts->method = FM_METHOD_AUTOMATIC; return true;
        }
        if (rhs->type == EXPR_STRING) {
            const char* s = rhs->data.string;
            if (strcmp(s, "Brent") == 0)             { opts->method = FM_METHOD_BRENT;        return true; }
            if (strcmp(s, "QuasiNewton") == 0)       { opts->method = FM_METHOD_QUASINEWTON;  return true; }
            if (strcmp(s, "ConjugateGradient") == 0) { opts->method = FM_METHOD_CONJGRAD;     return true; }
            if (strcmp(s, "Newton") == 0)            { opts->method = FM_METHOD_NEWTON;       return true; }
            if (strcmp(s, "LevenbergMarquardt") == 0
             || strcmp(s, "PrincipalAxis") == 0
             || strcmp(s, "InteriorPoint") == 0
             || strcmp(s, "LinearProgramming") == 0) {
                fm_warn(g_fm_name, "nimpl", "Method \"%s\" is not yet implemented", s);
                return false;
            }
        }
        fm_warn(g_fm_name, "badmeth", "unknown Method value");
        return false;
    }
    if (name == SYM_WorkingPrecision) {
        if (!fm_parse_working_precision(rhs, &opts->prec_mode, &opts->wp_bits)) {
            fm_warn(g_fm_name, "badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_MaxIterations) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer > 0) {
            opts->max_iter = rhs->data.integer;
            return true;
        }
        fm_warn(g_fm_name, "badopt", "MaxIterations must be a positive integer");
        return false;
    }
    if (name == SYM_AccuracyGoal)  return fm_parse_goal(rhs, &opts->acc_goal_digits);
    if (name == SYM_PrecisionGoal) return fm_parse_goal(rhs, &opts->prec_goal_digits);
    if (name == SYM_Gradient) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
            opts->gradient = NULL; return true;
        }
        opts->gradient = rhs; /* borrowed; verified by caller */
        return true;
    }
    if (name == SYM_StepMonitor)       { opts->step_monitor = rhs; return true; }
    if (name == SYM_EvaluationMonitor) { opts->eval_monitor = rhs; return true; }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Variable binding (Block semantics) — copied from findroot          *
 * ------------------------------------------------------------------ */

static void fm_bind_snapshot(FmVarBind* b, const char* name) {
    b->name = name;
    SymbolDef* def = symtab_get_def(name);
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

static void fm_bind_set(FmVarBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
    Expr* sym = expr_new_symbol(b->name);
    symtab_add_own_value(b->name, sym, value);
    expr_free(sym);
}

/* Free any temp OwnValue we installed but DO NOT yet restore the saved
 * chain — leaves the symbol unbound (free) so that subsequent expression
 * construction can copy it without triggering OwnValue replacement. */
static void fm_bind_clear_temp(FmVarBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = NULL;
    eval_clock_bump();
}

static void fm_bind_restore(FmVarBind* b) {
    if (!b->valid) return;
    SymbolDef* def = symtab_get_def(b->name);
    Rule* curr = def->own_values;
    while (curr) {
        Rule* next = curr->next;
        expr_free(curr->pattern);
        expr_free(curr->replacement);
        free(curr);
        curr = next;
    }
    def->own_values = b->saved_own;
    def->attributes = b->saved_attrs;
    b->valid = false;
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Evaluation helpers                                                 *
 * ------------------------------------------------------------------ */

static void fm_fire_monitor(Expr* monitor) {
    if (!monitor) return;
    Expr* tmp = eval_and_free(expr_copy(monitor));
    expr_free(tmp);
}

static NumericSpec fm_numeric_spec(const FmOpts* opts) {
    NumericSpec s;
#ifdef USE_MPFR
    if (opts->prec_mode == FM_PREC_MPFR) {
        s.mode = NUMERIC_MODE_MPFR;
        s.bits = opts->wp_bits;
        return s;
    }
#else
    (void)opts;
#endif
    return numeric_machine_spec();
}

/* Evaluate `f` with the variable bindings installed at `values[i]`, fire
 * EvaluationMonitor, and numericalize the result at the requested
 * precision so that Power[E, 1.0]-style residues collapse to numbers.
 * Returns NULL on failure; caller owns the return. */
static Expr* fm_eval_with_bindings(Expr* f, FmVarBind* binds,
                                   Expr* const* values, size_t n,
                                   Expr* eval_monitor,
                                   NumericSpec spec) {
    for (size_t i = 0; i < n; i++) fm_bind_set(&binds[i], values[i]);
    eval_clock_bump();
    fm_fire_monitor(eval_monitor);
    Expr* raw = eval_and_free(expr_copy(f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, spec);
    expr_free(raw);
    return num;
}

/* Evaluate the bound objective and return a double; NULL on failure. */
static bool fm_eval_scalar(Expr* f, FmVarBind* binds,
                           const double* x, size_t n,
                           const FmOpts* opts, double* out) {
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x[i]);
    Expr* res = fm_eval_with_bindings(f, binds, xv, n,
                                      opts->eval_monitor,
                                      fm_numeric_spec(opts));
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    if (!res) return false;
    bool ok = fm_expr_to_double_real(res, out);
    expr_free(res);
    return ok;
}

/* Build (and evaluate) D[f, x]. */
static Expr* fm_compute_partial(Expr* f, Expr* var) {
    Expr* args[2] = { expr_copy(f), expr_copy(var) };
    Expr* call = expr_new_function(expr_new_symbol("D"), args, 2);
    return eval_and_free(call);
}

/* Compute the symbolic gradient — list of D[f, x_i]. Returns NULL if any
 * partial fails to evaluate to something usable; caller takes ownership. */
static Expr** fm_compute_gradient(Expr* f, Expr** vars, size_t n) {
    Expr** g = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) g[i] = NULL;
    for (size_t i = 0; i < n; i++) {
        g[i] = fm_compute_partial(f, vars[i]);
        if (!g[i]) {
            for (size_t j = 0; j <= i; j++) expr_free(g[j]);
            free(g);
            return NULL;
        }
    }
    return g;
}

/* Compute the symbolic Hessian — n×n array of D[D[f, x_i], x_j]. */
static Expr*** fm_compute_hessian(Expr* f, Expr** vars, size_t n) {
    Expr*** H = (Expr***)malloc(sizeof(Expr**) * n);
    for (size_t i = 0; i < n; i++) H[i] = (Expr**)calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr* dfi = fm_compute_partial(f, vars[i]);
        if (!dfi) goto fail;
        for (size_t j = 0; j < n; j++) {
            H[i][j] = fm_compute_partial(dfi, vars[j]);
            if (!H[i][j]) { expr_free(dfi); goto fail; }
        }
        expr_free(dfi);
    }
    return H;
fail:
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) expr_free(H[i][j]);
        free(H[i]);
    }
    free(H);
    return NULL;
}

/* Numeric gradient via central differences when the symbolic one fails
 * or produces non-numeric residues. */
static bool fm_grad_finite_diff(Expr* f, FmVarBind* binds,
                                const double* x, size_t n,
                                const FmOpts* opts, double* g_out) {
    const double h = 1e-7;
    double* xp = (double*)malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++) xp[i] = x[i];
    for (size_t i = 0; i < n; i++) {
        double xi = x[i];
        double s = (fabs(xi) > 1.0 ? fabs(xi) : 1.0) * h;
        xp[i] = xi + s; double f1;
        if (!fm_eval_scalar(f, binds, xp, n, opts, &f1)) { free(xp); return false; }
        xp[i] = xi - s; double f0;
        if (!fm_eval_scalar(f, binds, xp, n, opts, &f0)) { free(xp); return false; }
        g_out[i] = (f1 - f0) / (2.0 * s);
        xp[i] = xi;
    }
    free(xp);
    return true;
}

/* Evaluate the symbolic gradient g_expr[i] at the current point. Returns
 * false if any component is non-numeric (caller may retry via FD). */
static bool fm_eval_gradient(Expr** g_exprs, FmVarBind* binds,
                             const double* x, size_t n,
                             const FmOpts* opts, double* g_out) {
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x[i]);
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        Expr* gi = fm_eval_with_bindings(g_exprs[i], binds, xv, n,
                                         opts->eval_monitor,
                                         fm_numeric_spec(opts));
        if (!gi || !fm_expr_to_double_real(gi, &g_out[i])) { ok = false; expr_free(gi); break; }
        expr_free(gi);
    }
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    return ok;
}

/* Evaluate symbolic Hessian H_exprs[i][j] at x. */
static bool fm_eval_hessian(Expr*** H_exprs, FmVarBind* binds,
                            const double* x, size_t n,
                            const FmOpts* opts, double* H_out /* n*n */) {
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x[i]);
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        for (size_t j = 0; j < n && ok; j++) {
            Expr* hij = fm_eval_with_bindings(H_exprs[i][j], binds, xv, n,
                                              opts->eval_monitor,
                                              fm_numeric_spec(opts));
            if (!hij || !fm_expr_to_double_real(hij, &H_out[i*n + j])) {
                ok = false; expr_free(hij); break;
            }
            expr_free(hij);
        }
    }
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    return ok;
}

/* Build the final result list  { fmin, { x->v1, y->v2, ... } } . The vars
 * are borrowed; vals_doubles is consumed. */
static Expr* fm_build_result(double fmin, Expr** vars, const double* vals,
                             size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r_args[2] = { expr_copy(vars[i]), expr_new_real(vals[i]) };
        rules[i] = expr_new_function(expr_new_symbol("Rule"), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol("List"), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_real(fmin), rule_list };
    return expr_new_function(expr_new_symbol("List"), top_args, 2);
}

/* ------------------------------------------------------------------ *
 *  Variable spec parsing                                              *
 * ------------------------------------------------------------------ */

static FmSpecKind fm_parse_var_spec(Expr* spec, Expr** var_out,
                                    Expr** x0_out, Expr** x1_out,
                                    Expr** xmin_out, Expr** xmax_out) {
    *var_out = NULL;
    *x0_out = *x1_out = *xmin_out = *xmax_out = NULL;
    if (!spec) return FM_SPEC_BAD;
    if (spec->type == EXPR_SYMBOL) {
        /* Bare variable, e.g. FindMinimum[f, {x, y, ...}] entry. */
        *var_out = spec;
        *x0_out = expr_new_real(0.0);
        return FM_SPEC_VAR_ONLY;
    }
    if (spec->type != EXPR_FUNCTION) return FM_SPEC_BAD;
    if (spec->data.function.head->type != EXPR_SYMBOL) return FM_SPEC_BAD;
    if (spec->data.function.head->data.symbol != SYM_List) return FM_SPEC_BAD;

    size_t n = spec->data.function.arg_count;
    if (n < 1 || n > 4) return FM_SPEC_BAD;

    Expr* var = spec->data.function.args[0];
    if (var->type != EXPR_SYMBOL) return FM_SPEC_BAD;
    *var_out = var;

    if (n == 1) {
        *x0_out = expr_new_real(0.0);
        return FM_SPEC_VAR_ONLY;
    }
    Expr* x0_raw = spec->data.function.args[1];
    if (x0_raw->type == EXPR_FUNCTION
        && x0_raw->data.function.head->type == EXPR_SYMBOL
        && x0_raw->data.function.head->data.symbol == SYM_List) {
        fm_warn(g_fm_name, "vecvar", "vector-valued variables are not yet supported");
        return FM_SPEC_BAD;
    }
    *x0_out = eval_and_free(expr_copy(x0_raw));
    if (n == 2) return FM_SPEC_SINGLE;
    if (n == 3) {
        *x1_out = eval_and_free(expr_copy(spec->data.function.args[2]));
        return FM_SPEC_TWO_START;
    }
    *xmin_out = eval_and_free(expr_copy(spec->data.function.args[2]));
    *xmax_out = eval_and_free(expr_copy(spec->data.function.args[3]));
    return FM_SPEC_BRACKET;
}

/* ------------------------------------------------------------------ *
 *  Constraint parsing                                                 *
 * ------------------------------------------------------------------ */

/* Try to interpret `cmp` as `var op const` or `const op var`, where var is
 * one of the iteration variables. On success update box[vidx]. Returns
 * true if recognised AS a box constraint. */
static bool fm_try_box_from_compare(Expr* cmp, Expr** vars, size_t nvars,
                                    FmBox* boxes) {
    if (cmp->type != EXPR_FUNCTION || cmp->data.function.arg_count != 2) return false;
    Expr* head = cmp->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;
    const char* op = head->data.symbol;
    if (op != SYM_Less && op != SYM_LessEqual
     && op != SYM_Greater && op != SYM_GreaterEqual) return false;
    Expr* a = cmp->data.function.args[0];
    Expr* b = cmp->data.function.args[1];

    /* Identify which side is a variable. */
    int64_t var_idx = -1;
    Expr* c_side = NULL;
    bool var_left = false;
    if (a->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < nvars; i++) {
            if (vars[i]->data.symbol == a->data.symbol) {
                var_idx = (int64_t)i; c_side = b; var_left = true; break;
            }
        }
    }
    if (var_idx < 0 && b->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < nvars; i++) {
            if (vars[i]->data.symbol == b->data.symbol) {
                var_idx = (int64_t)i; c_side = a; var_left = false; break;
            }
        }
    }
    if (var_idx < 0) return false;
    double c;
    if (!fm_expr_to_double_real(c_side, &c)) return false;

    /* Re-interpret op so that we describe the variable's allowed range.
     * If a bound was already recorded, tighten (intersection) rather than
     * replace. */
    bool is_le_lt;
    if (var_left) {
        is_le_lt = (op == SYM_LessEqual || op == SYM_Less);
        if (is_le_lt) {
            if (!boxes[var_idx].has_hi || c < boxes[var_idx].hi) boxes[var_idx].hi = c;
            boxes[var_idx].has_hi = true;
        } else {
            if (!boxes[var_idx].has_lo || c > boxes[var_idx].lo) boxes[var_idx].lo = c;
            boxes[var_idx].has_lo = true;
        }
    } else {
        /* c op var */
        is_le_lt = (op == SYM_LessEqual || op == SYM_Less);
        if (is_le_lt) {
            /* c <= var → var >= c → lower bound */
            if (!boxes[var_idx].has_lo || c > boxes[var_idx].lo) boxes[var_idx].lo = c;
            boxes[var_idx].has_lo = true;
        } else {
            /* c >= var → var <= c → upper bound */
            if (!boxes[var_idx].has_hi || c < boxes[var_idx].hi) boxes[var_idx].hi = c;
            boxes[var_idx].has_hi = true;
        }
    }
    return true;
}

/* Build an Expr g_expr such that the constraint cmp is equivalent to
 *   g_expr <= 0     (inequality)
 *   g_expr == 0     (equality)
 * for use inside the penalty wrapper. Caller takes ownership of *expr_out.
 * Returns false if the constraint shape is not understood. */
static bool fm_constraint_to_g(Expr* cmp, Expr** expr_out, bool* equality_out) {
    if (cmp->type != EXPR_FUNCTION || cmp->data.function.arg_count != 2) return false;
    Expr* head = cmp->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;
    const char* op = head->data.symbol;
    Expr* lhs = cmp->data.function.args[0];
    Expr* rhs = cmp->data.function.args[1];

    /* lhs - rhs */
    Expr* diff_args[2] = { expr_copy(lhs), expr_copy(rhs) };
    Expr* lhs_minus_rhs = expr_new_function(expr_new_symbol("Subtract"), diff_args, 2);

    if (op == SYM_LessEqual || op == SYM_Less) {
        /* lhs op rhs  →  lhs - rhs <= 0 */
        *expr_out = lhs_minus_rhs;
        *equality_out = false;
        return true;
    }
    if (op == SYM_GreaterEqual || op == SYM_Greater) {
        /* lhs op rhs  →  rhs - lhs <= 0  →  -(lhs - rhs) <= 0 */
        Expr* neg_args[2] = { expr_new_integer(-1), lhs_minus_rhs };
        *expr_out = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
        *equality_out = false;
        return true;
    }
    if (op == SYM_Equal) {
        *expr_out = lhs_minus_rhs;
        *equality_out = true;
        return true;
    }
    expr_free(lhs_minus_rhs);
    return false;
}

/* Walk an And[...] tree, an Inequality[...] node, or a single binary
 * comparison and accumulate boxes / general constraints. Returns false if
 * any branch is unsupported. */
static bool fm_collect_constraints(Expr* cons, Expr** vars, size_t nvars,
                                   FmBox* boxes,
                                   FmGenCon** gens_inout, size_t* ngens_inout,
                                   size_t* gens_cap_inout) {
    if (!cons) return true;
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol == SYM_And) {
        for (size_t i = 0; i < cons->data.function.arg_count; i++) {
            if (!fm_collect_constraints(cons->data.function.args[i], vars, nvars,
                                        boxes, gens_inout, ngens_inout,
                                        gens_cap_inout)) return false;
        }
        return true;
    }
    /* Inequality[v0, op0, v1, op1, ...] — the canonical chained-comparison
     * form produced by the parser. Treat each adjacent (v_i, op_i, v_{i+1})
     * triple as a separate binary comparison and reuse the existing
     * box-or-general classifier on each. */
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol == SYM_Inequality
        && cons->data.function.arg_count >= 3
        && (cons->data.function.arg_count & 1u) == 1) {
        size_t npairs = (cons->data.function.arg_count - 1) / 2;
        for (size_t k = 0; k < npairs; k++) {
            Expr* a  = cons->data.function.args[2*k];
            Expr* op = cons->data.function.args[2*k + 1];
            Expr* b  = cons->data.function.args[2*k + 2];
            if (op->type != EXPR_SYMBOL) {
                fm_warn(g_fm_name, "nimpl", "unsupported constraint shape");
                return false;
            }
            Expr* pair_args[2] = { expr_copy(a), expr_copy(b) };
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol),
                                           pair_args, 2);
            bool ok = fm_collect_constraints(pair, vars, nvars, boxes,
                                             gens_inout, ngens_inout,
                                             gens_cap_inout);
            expr_free(pair);
            if (!ok) return false;
        }
        return true;
    }
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol == SYM_Or) {
        fm_warn(g_fm_name, "nimpl", "disjunctive (Or) constraints are not yet supported");
        return false;
    }
    /* Reject Element[...] (e.g. x ∈ Integers) outright. */
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol == SYM_Element) {
        fm_warn(g_fm_name, "nimpl", "Element / domain constraints are not yet supported");
        return false;
    }
    /* Single comparison. Try as box first. */
    if (fm_try_box_from_compare(cons, vars, nvars, boxes)) return true;
    /* Otherwise it's a general inequality or equality. */
    Expr* g_expr = NULL;
    bool eq = false;
    if (!fm_constraint_to_g(cons, &g_expr, &eq)) {
        fm_warn(g_fm_name, "nimpl", "unsupported constraint shape");
        return false;
    }
    if (*ngens_inout == *gens_cap_inout) {
        size_t nc = *gens_cap_inout ? (*gens_cap_inout) * 2 : 4;
        *gens_inout = (FmGenCon*)realloc(*gens_inout, sizeof(FmGenCon) * nc);
        *gens_cap_inout = nc;
    }
    (*gens_inout)[*ngens_inout].expr = g_expr;
    (*gens_inout)[*ngens_inout].equality = eq;
    (*ngens_inout)++;
    return true;
}

/* Project x in-place to the box. */
static void fm_project_box(double* x, size_t n, const FmBox* boxes) {
    for (size_t i = 0; i < n; i++) {
        if (boxes[i].has_lo && x[i] < boxes[i].lo) x[i] = boxes[i].lo;
        if (boxes[i].has_hi && x[i] > boxes[i].hi) x[i] = boxes[i].hi;
    }
}

/* Evaluate Σ max(0, g_i(x))^2 + Σ h_j(x)^2 over the general constraint
 * set. Returns false if any constraint cannot be evaluated. */
static bool fm_eval_penalty(const FmGenCon* gens, size_t ngens,
                            FmVarBind* binds, const double* x, size_t n,
                            const FmOpts* opts, double* pen_out) {
    if (ngens == 0) { *pen_out = 0.0; return true; }
    double total = 0.0;
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x[i]);
    for (size_t k = 0; k < ngens; k++) {
        Expr* v = fm_eval_with_bindings(gens[k].expr, binds, xv, n,
                                        opts->eval_monitor,
                                        fm_numeric_spec(opts));
        double d;
        if (!v || !fm_expr_to_double_real(v, &d)) {
            expr_free(v);
            for (size_t i = 0; i < n; i++) expr_free(xv[i]);
            free(xv);
            return false;
        }
        expr_free(v);
        if (gens[k].equality)        total += d * d;
        else if (d > 0.0)            total += d * d;
    }
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    *pen_out = total;
    return true;
}

/* Augmented objective at x: f(x) + μ * penalty(x). */
static bool fm_eval_augmented(Expr* f, FmVarBind* binds,
                              const double* x, size_t n,
                              const FmGenCon* gens, size_t ngens,
                              double mu, const FmOpts* opts, double* out) {
    double fv;
    if (!fm_eval_scalar(f, binds, x, n, opts, &fv)) return false;
    double pen;
    if (!fm_eval_penalty(gens, ngens, binds, x, n, opts, &pen)) return false;
    *out = fv + mu * pen;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Line search (backtracking with Armijo)                             *
 * ------------------------------------------------------------------ */

/* Given x, descent direction d, gradient g, and current f0, find an
 * alpha satisfying the Armijo condition. The augmented objective handles
 * f + μ * penalty; pass mu == 0.0 and gens == NULL for plain f.
 * Result: *alpha_out, *f_out (function value at x + alpha*d).
 * Optionally projects after the step (box constraints). */
static bool fm_line_search(Expr* f, FmVarBind* binds, size_t n,
                           const double* x, const double* d,
                           double f0, double g_dot_d,
                           const FmGenCon* gens, size_t ngens, double mu,
                           const FmBox* boxes,
                           const FmOpts* opts,
                           double* alpha_out, double* f_out, double* x_out) {
    const double c1 = 1e-4;
    /* Cap initial step so ||alpha*d|| <= 1: protects against huge initial
     * gradients in the penalty path. */
    double dnorm = 0.0;
    for (size_t i = 0; i < n; i++) dnorm += d[i] * d[i];
    dnorm = sqrt(dnorm);
    double alpha = (dnorm > 1.0) ? 1.0 / dnorm : 1.0;
    for (int k = 0; k < 30; k++) {
        for (size_t i = 0; i < n; i++) x_out[i] = x[i] + alpha * d[i];
        if (boxes) fm_project_box(x_out, n, boxes);
        double fnew;
        bool ok = (mu > 0.0 && gens)
            ? fm_eval_augmented(f, binds, x_out, n, gens, ngens, mu, opts, &fnew)
            : fm_eval_scalar(f, binds, x_out, n, opts, &fnew);
        if (!ok) { alpha *= 0.5; continue; }
        /* When projection is in effect the Armijo bound uses g . (x_proj - x). */
        if (boxes) {
            double gd = 0.0;
            for (size_t i = 0; i < n; i++) gd += d[i] * (x_out[i] - x[i]);
            /* gd is along the projected step; reuse Armijo with g_dot_d scaled */
            (void)gd; /* Simple acceptance: just need f decrease */
            if (fnew <= f0 - 1e-12 * fabs(f0) - 1e-300) {
                *alpha_out = alpha; *f_out = fnew; return true;
            }
        } else if (fnew <= f0 + c1 * alpha * g_dot_d) {
            *alpha_out = alpha; *f_out = fnew; return true;
        }
        alpha *= 0.5;
        if (alpha < 1e-20) break;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  1D bracketing (mnbrak-style)                                        *
 * ------------------------------------------------------------------ */

/* Given a single start x0, find a, b, c with f(a) > f(b) and f(c) > f(b),
 * bracketing a local minimum. Honors optional box constraints. */
static bool fm_bracket(Expr* f, FmVarBind* binds, const FmOpts* opts,
                       double x0, const FmBox* box1,
                       double* a_out, double* b_out, double* c_out) {
    double a = x0, b, fa, fb;
    if (!fm_eval_scalar(f, binds, &a, 1, opts, &fa)) return false;
    double h = (fabs(a) > 1.0 ? fabs(a) : 1.0) * 1e-2;
    b = a + h;
    if (box1 && box1->has_hi && b > box1->hi) b = (a + box1->hi) * 0.5;
    if (box1 && box1->has_lo && b < box1->lo) b = (a + box1->lo) * 0.5;
    if (!fm_eval_scalar(f, binds, &b, 1, opts, &fb)) return false;
    if (fb > fa) {
        /* Step the other way. */
        double t = a; a = b; b = t;
        t = fa; fa = fb; fb = t;
    }
    double c = b + 1.618 * (b - a);
    if (box1 && box1->has_hi && c > box1->hi) c = box1->hi;
    if (box1 && box1->has_lo && c < box1->lo) c = box1->lo;
    double fc;
    if (!fm_eval_scalar(f, binds, &c, 1, opts, &fc)) return false;
    for (int k = 0; k < 100 && fc <= fb; k++) {
        a = b; fa = fb;
        b = c; fb = fc;
        c = b + 1.618 * (b - a);
        if (box1 && box1->has_hi && c >= box1->hi) {
            c = box1->hi;
            if (!fm_eval_scalar(f, binds, &c, 1, opts, &fc)) return false;
            break;
        }
        if (box1 && box1->has_lo && c <= box1->lo) {
            c = box1->lo;
            if (!fm_eval_scalar(f, binds, &c, 1, opts, &fc)) return false;
            break;
        }
        if (!fm_eval_scalar(f, binds, &c, 1, opts, &fc)) return false;
    }
    if (a > c) {
        double t = a; a = c; c = t;
    }
    *a_out = a; *b_out = b; *c_out = c;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Brent's minimisation (1D, machine precision)                        *
 * ------------------------------------------------------------------ */

#define FM_CGOLD 0.3819660112501051
#define FM_ZEPS  1.0e-12

static bool fm_brent_min(Expr* f, FmVarBind* bind, const FmOpts* opts,
                         double a, double b, double c,
                         const FmBox* box1,
                         double* x_out, double* fx_out) {
    if (a > c) { double t = a; a = c; c = t; }
    double tol = pow(10.0, -opts->prec_goal_digits);
    double tol_acc = pow(10.0, -opts->acc_goal_digits);
    double e = 0.0, d = 0.0;
    double x, w, v;
    x = w = v = b;
    /* Ensure x is inside [a, c]. */
    if (x < a || x > c) x = w = v = 0.5 * (a + c);
    double fx;
    if (!fm_eval_scalar(f, bind, &x, 1, opts, &fx)) return false;
    double fw = fx, fv = fx;
    for (int64_t k = 0; k < opts->max_iter; k++) {
        double xm = 0.5 * (a + c);
        double tol1 = tol * fabs(x) + FM_ZEPS;
        double tol2 = 2.0 * tol1;
        if (fabs(x - xm) <= tol2 - 0.5 * (c - a)
            || fabs(fx) < tol_acc * (1.0 + fabs(fx))) {
            *x_out = x; *fx_out = fx; return true;
        }
        double u;
        if (fabs(e) > tol1) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            double etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp)
                || p <= q * (a - x) || p >= q * (c - x)) {
                e = (x >= xm) ? (a - x) : (c - x);
                d = FM_CGOLD * e;
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || c - u < tol2)
                    d = (xm - x >= 0.0) ? tol1 : -tol1;
            }
        } else {
            e = (x >= xm) ? (a - x) : (c - x);
            d = FM_CGOLD * e;
        }
        u = (fabs(d) >= tol1) ? (x + d) : (x + ((d >= 0.0) ? tol1 : -tol1));
        if (box1) {
            if (box1->has_lo && u < box1->lo) u = box1->lo;
            if (box1->has_hi && u > box1->hi) u = box1->hi;
        }
        double fu;
        if (!fm_eval_scalar(f, bind, &u, 1, opts, &fu)) return false;
        fm_fire_monitor(opts->step_monitor);
        if (fu <= fx) {
            if (u >= x) a = x; else c = x;
            v = w; w = x; x = u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) a = u; else c = u;
            if (fu <= fw || w == x) { v = w; w = u; fv = fw; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    fm_warn(g_fm_name, "cvmit", "Brent failed to converge within %lld iterations",
            (long long)opts->max_iter);
    *x_out = x; *fx_out = fx;
    return true;
}

/* ------------------------------------------------------------------ *
 *  In-place modified Cholesky for Newton                              *
 * ------------------------------------------------------------------ */

/* Factor H + τI in-place: result stored in lower-triangle. Returns false
 * if SPD factorisation fails. */
static bool fm_chol_factor(double* H, size_t n, double tau) {
    /* Make a working copy because we may need to retry with larger τ. */
    for (size_t i = 0; i < n; i++) H[i*n + i] += tau;
    for (size_t j = 0; j < n; j++) {
        double sum = H[j*n + j];
        for (size_t k = 0; k < j; k++) sum -= H[j*n + k] * H[j*n + k];
        if (sum <= 1e-14) return false;
        H[j*n + j] = sqrt(sum);
        for (size_t i = j + 1; i < n; i++) {
            double s = H[i*n + j];
            for (size_t k = 0; k < j; k++) s -= H[i*n + k] * H[j*n + k];
            H[i*n + j] = s / H[j*n + j];
        }
    }
    return true;
}

static void fm_chol_solve(const double* L, size_t n,
                          const double* b, double* x) {
    /* Forward: L y = b. */
    double* y = (double*)malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++) {
        double s = b[i];
        for (size_t k = 0; k < i; k++) s -= L[i*n + k] * y[k];
        y[i] = s / L[i*n + i];
    }
    /* Backward: L^T x = y. */
    for (size_t i_p1 = n; i_p1 > 0; i_p1--) {
        size_t i = i_p1 - 1;
        double s = y[i];
        for (size_t k = i + 1; k < n; k++) s -= L[k*n + i] * x[k];
        x[i] = s / L[i*n + i];
    }
    free(y);
}

/* ------------------------------------------------------------------ *
 *  BFGS (machine precision)                                            *
 * ------------------------------------------------------------------ */

static bool fm_run_bfgs(Expr* f, Expr** vars, size_t n,
                        FmVarBind* binds, Expr** g_exprs,
                        double* x, /* in/out */
                        const FmGenCon* gens, size_t ngens, double mu,
                        const FmBox* boxes,
                        const FmOpts* opts,
                        double* fx_out) {
    (void)vars;
    /* Inverse Hessian approximation H, stored row-major. Initialise to I. */
    double* H = (double*)calloc(n * n, sizeof(double));
    double* g = (double*)malloc(sizeof(double) * n);
    double* g_new = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* s = (double*)malloc(sizeof(double) * n);
    double* y = (double*)malloc(sizeof(double) * n);
    double* Hy = (double*)malloc(sizeof(double) * n);
    bool ok = false;

    for (size_t i = 0; i < n; i++) H[i*n + i] = 1.0;
    if (boxes) fm_project_box(x, n, boxes);

    double fx;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    if (augmented) {
        if (!fm_eval_augmented(f, binds, x, n, gens, ngens, mu, opts, &fx)) goto cleanup;
    } else {
        if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) goto cleanup;
    }

    bool got_grad = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
    if (!got_grad) {
        if (!fm_grad_finite_diff(f, binds, x, n, opts, g)) {
            fm_warn(g_fm_name, "nlnum", "gradient evaluation failed at start point");
            goto cleanup;
        }
    }

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* Gradient norm convergence. */
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { ok = true; break; }

        /* d = -H g. */
        for (size_t i = 0; i < n; i++) {
            double s_ = 0.0;
            for (size_t j = 0; j < n; j++) s_ += H[i*n + j] * g[j];
            d[i] = -s_;
        }
        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) {
            /* Not a descent direction — reset H to I and use steepest. */
            for (size_t i = 0; i < n*n; i++) H[i] = 0.0;
            for (size_t i = 0; i < n; i++) { H[i*n + i] = 1.0; d[i] = -g[i]; }
            g_dot_d = 0.0; for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        }

        double alpha, fx_new;
        bool ls_ok;
        if (augmented) {
            ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                   gens, ngens, mu, boxes, opts,
                                   &alpha, &fx_new, x_new);
        } else {
            ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                   NULL, 0, 0.0, boxes, opts,
                                   &alpha, &fx_new, x_new);
        }
        if (!ls_ok) {
            fm_warn(g_fm_name, "lstol", "line search failed at iter %lld",
                    (long long)k);
            /* Return current best (still ok=false, but caller may use x). */
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        /* Step magnitude check (PrecisionGoal). */
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
        }

        /* Compute new gradient. */
        bool ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
        if (!ng_ok) {
            if (!fm_grad_finite_diff(f, binds, x_new, n, opts, g_new)) {
                fm_warn(g_fm_name, "nlnum", "gradient evaluation failed in iteration");
                /* Take the step and stop. */
                for (size_t i = 0; i < n; i++) x[i] = x_new[i];
                fx = fx_new;
                break;
            }
        }

        /* BFGS update: s = x_new - x; y = g_new - g; ρ = 1 / (y . s). */
        for (size_t i = 0; i < n; i++) { s[i] = x_new[i] - x[i]; y[i] = g_new[i] - g[i]; }
        double sy = 0.0;
        for (size_t i = 0; i < n; i++) sy += s[i] * y[i];
        if (sy > 1e-12) {
            double rho = 1.0 / sy;
            /* Hy = H y. */
            for (size_t i = 0; i < n; i++) {
                double t = 0.0;
                for (size_t j = 0; j < n; j++) t += H[i*n + j] * y[j];
                Hy[i] = t;
            }
            double yHy = 0.0;
            for (size_t i = 0; i < n; i++) yHy += y[i] * Hy[i];
            /* H ← H + ((sy + yHy) ρ²) s s^T − ρ (Hy s^T + s (Hy)^T). */
            double coef = (sy + yHy) * rho * rho;
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    H[i*n + j] += coef * s[i] * s[j]
                                - rho * (Hy[i] * s[j] + s[i] * Hy[j]);
                }
            }
        }

        for (size_t i = 0; i < n; i++) { x[i] = x_new[i]; g[i] = g_new[i]; }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }
    }
    if (!ok) {
        /* Either max iters or line search exhausted — still report best. */
    }
    *fx_out = fx;
    /* Always return true so the driver gets the best iterate; warnings
     * already emitted above when convergence failed. */
    ok = true;
cleanup:
    free(H); free(g); free(g_new); free(d); free(x_new);
    free(s); free(y); free(Hy);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Conjugate gradient (Polak-Ribière+ with restart)                    *
 * ------------------------------------------------------------------ */

static bool fm_run_cg(Expr* f, Expr** vars, size_t n,
                      FmVarBind* binds, Expr** g_exprs,
                      double* x,
                      const FmGenCon* gens, size_t ngens, double mu,
                      const FmBox* boxes,
                      const FmOpts* opts,
                      double* fx_out) {
    (void)vars;
    double* g = (double*)malloc(sizeof(double) * n);
    double* g_new = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    bool ok = false;

    if (boxes) fm_project_box(x, n, boxes);
    double fx;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    if (augmented) {
        if (!fm_eval_augmented(f, binds, x, n, gens, ngens, mu, opts, &fx)) goto cleanup;
    } else {
        if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) goto cleanup;
    }
    bool got_grad = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
    if (!got_grad) {
        if (!fm_grad_finite_diff(f, binds, x, n, opts, g)) {
            fm_warn(g_fm_name, "nlnum", "gradient failed at start point");
            goto cleanup;
        }
    }
    for (size_t i = 0; i < n; i++) d[i] = -g[i];

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { ok = true; break; }

        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) {
            /* Restart with steepest descent. */
            for (size_t i = 0; i < n; i++) d[i] = -g[i];
            g_dot_d = 0.0; for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        }
        double alpha, fx_new;
        bool ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                    augmented ? gens : NULL,
                                    augmented ? ngens : 0,
                                    augmented ? mu : 0.0,
                                    boxes, opts, &alpha, &fx_new, x_new);
        if (!ls_ok) { fm_warn(g_fm_name, "lstol", "line search failed"); break; }
        fm_fire_monitor(opts->step_monitor);

        bool ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
        if (!ng_ok) {
            if (!fm_grad_finite_diff(f, binds, x_new, n, opts, g_new)) {
                for (size_t i = 0; i < n; i++) x[i] = x_new[i]; fx = fx_new; break;
            }
        }
        /* Polak-Ribière+. */
        double num = 0.0, den = 0.0;
        for (size_t i = 0; i < n; i++) {
            num += g_new[i] * (g_new[i] - g[i]);
            den += g[i] * g[i];
        }
        double beta = (den > 0.0) ? num / den : 0.0;
        if (beta < 0.0) beta = 0.0;
        if ((k + 1) % n == 0) beta = 0.0;
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
        }
        for (size_t i = 0; i < n; i++) {
            d[i] = -g_new[i] + beta * d[i];
            x[i] = x_new[i];
            g[i] = g_new[i];
        }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }
    }
    *fx_out = fx;
    ok = true;
cleanup:
    free(g); free(g_new); free(d); free(x_new);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Newton (machine precision, modified Cholesky)                       *
 * ------------------------------------------------------------------ */

static bool fm_run_newton(Expr* f, Expr** vars, size_t n,
                          FmVarBind* binds, Expr** g_exprs, Expr*** H_exprs,
                          double* x,
                          const FmGenCon* gens, size_t ngens, double mu,
                          const FmBox* boxes,
                          const FmOpts* opts,
                          double* fx_out) {
    (void)vars;
    double* g = (double*)malloc(sizeof(double) * n);
    double* d = (double*)malloc(sizeof(double) * n);
    double* x_new = (double*)malloc(sizeof(double) * n);
    double* H = (double*)malloc(sizeof(double) * n * n);
    double* Hcopy = (double*)malloc(sizeof(double) * n * n);
    double* neg_g = (double*)malloc(sizeof(double) * n);
    bool ok = false;
    bool augmented = (mu > 0.0 && gens && ngens > 0);
    if (boxes) fm_project_box(x, n, boxes);

    double fx;
    if (augmented) {
        if (!fm_eval_augmented(f, binds, x, n, gens, ngens, mu, opts, &fx)) goto cleanup;
    } else {
        if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) goto cleanup;
    }

    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        bool gok = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
        if (!gok) {
            if (!fm_grad_finite_diff(f, binds, x, n, opts, g)) {
                fm_warn(g_fm_name, "nlnum", "gradient failed during Newton");
                goto cleanup;
            }
        }
        double gnorm = 0.0;
        for (size_t i = 0; i < n; i++) gnorm += g[i] * g[i];
        gnorm = sqrt(gnorm);
        if (gnorm < tol_acc) { ok = true; break; }

        bool Hok = H_exprs && fm_eval_hessian(H_exprs, binds, x, n, opts, H);
        if (!Hok) {
            /* Fall back to BFGS-style steepest. */
            for (size_t i = 0; i < n; i++) d[i] = -g[i];
        } else {
            /* Try Cholesky with increasing τ. */
            double tau = 0.0;
            bool factored = false;
            for (int t = 0; t < 30 && !factored; t++) {
                for (size_t i = 0; i < n*n; i++) Hcopy[i] = H[i];
                factored = fm_chol_factor(Hcopy, n, tau);
                if (!factored) tau = (tau == 0.0) ? 1e-3 : tau * 2.0;
            }
            if (!factored) {
                fm_warn(g_fm_name, "dsing", "Hessian not positive definite");
                for (size_t i = 0; i < n; i++) d[i] = -g[i];
            } else {
                for (size_t i = 0; i < n; i++) neg_g[i] = -g[i];
                fm_chol_solve(Hcopy, n, neg_g, d);
            }
        }
        double g_dot_d = 0.0;
        for (size_t i = 0; i < n; i++) g_dot_d += g[i] * d[i];
        if (g_dot_d >= 0.0) { for (size_t i = 0; i < n; i++) d[i] = -g[i];
                              g_dot_d = 0.0;
                              for (size_t i = 0; i < n; i++) g_dot_d += g[i]*d[i]; }
        double alpha, fx_new;
        bool ls_ok = fm_line_search(f, binds, n, x, d, fx, g_dot_d,
                                    augmented ? gens : NULL,
                                    augmented ? ngens : 0,
                                    augmented ? mu : 0.0,
                                    boxes, opts, &alpha, &fx_new, x_new);
        if (!ls_ok) { fm_warn(g_fm_name, "lstol", "Newton line search failed"); break; }
        fm_fire_monitor(opts->step_monitor);
        double max_step = 0.0, max_x = 0.0;
        for (size_t i = 0; i < n; i++) {
            double ds = fabs(x_new[i] - x[i]);
            if (ds > max_step) max_step = ds;
            if (fabs(x_new[i]) > max_x) max_x = fabs(x_new[i]);
            x[i] = x_new[i];
        }
        fx = fx_new;
        if (max_step < tol_prec * (max_x + 1e-300)) { ok = true; break; }
    }
    *fx_out = fx;
    ok = true;
cleanup:
    free(g); free(d); free(x_new); free(H); free(Hcopy); free(neg_g);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Penalty outer loop                                                  *
 * ------------------------------------------------------------------ */

static bool fm_run_penalty(Expr* f, Expr** vars, size_t n,
                           FmVarBind* binds, FmMethod method,
                           Expr** g_exprs, Expr*** H_exprs,
                           double* x, /* in/out */
                           const FmGenCon* gens, size_t ngens,
                           const FmBox* boxes,
                           const FmOpts* opts,
                           double* fx_out) {
    /* Outer μ schedule: 1 → 10 → ... up to 10^8 or until feasible.
     * Starting small keeps the augmented-objective gradient manageable
     * for infeasible starts. */
    double mu = 1.0;
    double fx = 0.0;
    bool feas = false;
    for (int round = 0; round < 9; round++) {
        bool ok;
        switch (method) {
            case FM_METHOD_QUASINEWTON:
                ok = fm_run_bfgs(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_CONJGRAD:
                ok = fm_run_cg(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            case FM_METHOD_NEWTON:
                ok = fm_run_newton(f, vars, n, binds, g_exprs, H_exprs, x, gens, ngens, mu, boxes, opts, &fx);
                break;
            default:
                ok = fm_run_bfgs(f, vars, n, binds, g_exprs, x, gens, ngens, mu, boxes, opts, &fx);
        }
        if (!ok) return false;
        double pen;
        if (!fm_eval_penalty(gens, ngens, binds, x, n, opts, &pen)) return false;
        if (pen < 1e-12) { feas = true; break; }
        mu *= 10.0;
    }
    if (!feas) {
        fm_warn(g_fm_name, "infeas", "could not satisfy constraints to tolerance");
    }
    /* Report unaugmented objective value at the final point. */
    if (!fm_eval_scalar(f, binds, x, n, opts, &fx)) return false;
    *fx_out = fx;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                              *
 * ------------------------------------------------------------------ */

static Expr* findmin_driver(Expr* res, const char* fn_name) {
    g_fm_name = fn_name;
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn(fn_name, "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && fm_is_option_arg(res->data.function.args[pos_end - 1])) pos_end--;
    for (size_t i = pos_end; i < argc; i++) {
        if (!fm_is_option_arg(res->data.function.args[i])) {
            fm_warn(fn_name, "badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) {
        fm_warn(fn_name, "argt", "needs exactly 2 positional arguments (got %zu)", pos_end);
        return NULL;
    }
    FmOpts opts;
    opts.method = FM_METHOD_AUTOMATIC;
    opts.prec_mode = FM_PREC_MACHINE;
    opts.wp_bits = 0;
    opts.max_iter = 500;
    opts.acc_goal_digits = -1.0;
    opts.prec_goal_digits = -1.0;
    opts.gradient = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;
    for (size_t i = pos_end; i < argc; i++) {
        if (!fm_apply_option(res->data.function.args[i], &opts)) return NULL;
    }
    double wp_digits = (opts.prec_mode == FM_PREC_MACHINE)
        ? NUMERIC_MACHINE_PRECISION_DIGITS
#ifdef USE_MPFR
        : numeric_bits_to_digits(opts.wp_bits);
#else
        : NUMERIC_MACHINE_PRECISION_DIGITS;
#endif
    if (opts.acc_goal_digits  < 0.0) opts.acc_goal_digits  = wp_digits / 2.0;
    if (opts.prec_goal_digits < 0.0) opts.prec_goal_digits = wp_digits / 2.0;

    /* Detect {f, cons} form. */
    Expr* f_arg = res->data.function.args[0];
    Expr* var_arg = res->data.function.args[1];
    Expr* f_raw = f_arg;
    Expr* cons = NULL;
    if (f_arg->type == EXPR_FUNCTION
        && f_arg->data.function.head->type == EXPR_SYMBOL
        && f_arg->data.function.head->data.symbol == SYM_List
        && f_arg->data.function.arg_count == 2) {
        f_raw = f_arg->data.function.args[0];
        cons = f_arg->data.function.args[1];
    }

    /* Parse variables. var_arg may be:
     *   {x}  /  {x, x0}  /  {x, x0, x1}  /  {x, xstart, xmin, xmax}      (scalar)
     *   {{x, ...}, {y, ...}, ...}                                        (vector)
     *   {x, y, z}  (each scalar element is a bare symbol, treated as {x_i, 0}) */
    bool is_system = false;
    if (var_arg->type == EXPR_FUNCTION
        && var_arg->data.function.head->type == EXPR_SYMBOL
        && var_arg->data.function.head->data.symbol == SYM_List
        && var_arg->data.function.arg_count > 0) {
        size_t na = var_arg->data.function.arg_count;
        bool any_inner = false, all_inner_or_sym = true;
        for (size_t i = 0; i < na; i++) {
            Expr* e = var_arg->data.function.args[i];
            bool is_inner = (e->type == EXPR_FUNCTION
                && e->data.function.head->type == EXPR_SYMBOL
                && e->data.function.head->data.symbol == SYM_List);
            bool is_sym = (e->type == EXPR_SYMBOL);
            if (is_inner) any_inner = true;
            if (!is_inner && !is_sym) all_inner_or_sym = false;
        }
        /* {{x,x0},{y,y0}} → system; {x, y, z} (all bare symbols) → system n-D */
        if (any_inner && all_inner_or_sym) is_system = true;
        else if (na >= 2 && var_arg->data.function.args[0]->type == EXPR_SYMBOL) {
            /* Could be either {x, x0} (scalar) or {x, y} (multi-symbol). Check 2nd arg. */
            Expr* a1 = var_arg->data.function.args[1];
            if (a1->type == EXPR_SYMBOL && na > 1) {
                /* {x, y, ...}: all bare symbols → multi-var auto-start. */
                bool all_sym = true;
                for (size_t i = 0; i < na; i++) {
                    if (var_arg->data.function.args[i]->type != EXPR_SYMBOL) {
                        all_sym = false; break;
                    }
                }
                if (all_sym) is_system = true;
            }
        }
    }

    Expr** vars = NULL;
    double* x_vec = NULL;
    FmBox* boxes = NULL;
    FmVarBind* binds = NULL;
    FmGenCon* gens = NULL;
    size_t ngens = 0, gcap = 0;
    Expr** g_exprs = NULL;
    Expr*** H_exprs = NULL;
    Expr* result_out = NULL;
    size_t n = 0;

    if (is_system) {
        n = var_arg->data.function.arg_count;
        vars = (Expr**)calloc(n, sizeof(Expr*));
        x_vec = (double*)calloc(n, sizeof(double));
        boxes = (FmBox*)calloc(n, sizeof(FmBox));
        for (size_t i = 0; i < n; i++) {
            Expr* sub = var_arg->data.function.args[i];
            Expr *u, *x0e = NULL, *x1e = NULL, *xmin = NULL, *xmax = NULL;
            FmSpecKind k = fm_parse_var_spec(sub, &u, &x0e, &x1e, &xmin, &xmax);
            if (k != FM_SPEC_VAR_ONLY && k != FM_SPEC_SINGLE && k != FM_SPEC_BRACKET) {
                fm_warn(fn_name, "ivar", "variable spec %zu malformed", i);
                expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
                goto cleanup;
            }
            vars[i] = u;
            if (!x0e || !fm_expr_to_double_real(x0e, &x_vec[i])) x_vec[i] = 0.0;
            if (k == FM_SPEC_BRACKET) {
                double lo, hi;
                if (fm_expr_to_double_real(xmin, &lo) && fm_expr_to_double_real(xmax, &hi)) {
                    boxes[i].has_lo = true; boxes[i].lo = lo;
                    boxes[i].has_hi = true; boxes[i].hi = hi;
                }
            }
            expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
        }
    } else {
        n = 1;
        vars = (Expr**)calloc(1, sizeof(Expr*));
        x_vec = (double*)calloc(1, sizeof(double));
        boxes = (FmBox*)calloc(1, sizeof(FmBox));
        Expr *u, *x0e = NULL, *x1e = NULL, *xmin = NULL, *xmax = NULL;
        FmSpecKind k = fm_parse_var_spec(var_arg, &u, &x0e, &x1e, &xmin, &xmax);
        if (k == FM_SPEC_BAD) {
            fm_warn(fn_name, "ivar", "variable spec must be {x}, {x, x0}, {x, x0, x1}, or {x, xstart, xmin, xmax}");
            expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
            goto cleanup;
        }
        vars[0] = u;
        if (!x0e || !fm_expr_to_double_real(x0e, &x_vec[0])) x_vec[0] = 0.0;
        if (k == FM_SPEC_BRACKET) {
            double lo, hi;
            if (fm_expr_to_double_real(xmin, &lo) && fm_expr_to_double_real(xmax, &hi)) {
                boxes[0].has_lo = true; boxes[0].lo = lo;
                boxes[0].has_hi = true; boxes[0].hi = hi;
            }
        }
        /* Smuggle TWO_START into method selection via custom path: we
         * encode this by using boxes (lo=x0, hi=x1) if Brent and the
         * caller gave {var, x0, x1}. We'll handle it during method
         * dispatch below by detecting TWO_START separately. */
        if (k == FM_SPEC_TWO_START) {
            double a, b;
            if (fm_expr_to_double_real(x0e, &a) && fm_expr_to_double_real(x1e, &b)) {
                if (a > b) { double t = a; a = b; b = t; }
                boxes[0].has_lo = true; boxes[0].lo = a;
                boxes[0].has_hi = true; boxes[0].hi = b;
                /* For TWO_START with Automatic method we want Brent. */
                if (opts.method == FM_METHOD_AUTOMATIC) opts.method = FM_METHOD_BRENT;
            }
        }
        expr_free(x0e); expr_free(x1e); expr_free(xmin); expr_free(xmax);
    }

    /* Now bind variables. */
    binds = (FmVarBind*)calloc(n, sizeof(FmVarBind));
    for (size_t i = 0; i < n; i++) fm_bind_snapshot(&binds[i], vars[i]->data.symbol);

    /* Constraints. */
    if (cons) {
        if (!fm_collect_constraints(cons, vars, n, boxes, &gens, &ngens, &gcap))
            goto cleanup;
    }

    /* Method selection. */
    FmMethod method = opts.method;
    if (method == FM_METHOD_AUTOMATIC) {
        method = (n == 1) ? FM_METHOD_BRENT : FM_METHOD_QUASINEWTON;
    }

    /* Compute symbolic gradient/Hessian when needed. */
    bool needs_grad = (method == FM_METHOD_QUASINEWTON
                    || method == FM_METHOD_CONJGRAD
                    || method == FM_METHOD_NEWTON);
    bool needs_hess = (method == FM_METHOD_NEWTON);
    if (needs_grad) {
        if (opts.gradient
            && opts.gradient->type == EXPR_FUNCTION
            && opts.gradient->data.function.head->type == EXPR_SYMBOL
            && opts.gradient->data.function.head->data.symbol == SYM_List
            && opts.gradient->data.function.arg_count == n) {
            g_exprs = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) g_exprs[i] = expr_copy(opts.gradient->data.function.args[i]);
        } else {
            g_exprs = fm_compute_gradient(f_raw, vars, n);
            if (!g_exprs) {
                /* OK — will fall back to finite differences inside the solver. */
            }
        }
    }
    if (needs_hess) {
        H_exprs = fm_compute_hessian(f_raw, vars, n);
        /* OK if NULL — Newton will fall back to BFGS-style steepest. */
    }

    /* Dispatch. */
    double fx_min = 0.0;
    bool ok = true;
    bool has_general_cons = (ngens > 0);
    if (method == FM_METHOD_BRENT) {
        if (n != 1) {
            fm_warn(fn_name, "badmeth", "Method \"Brent\" requires a single variable");
            goto cleanup;
        }
        if (has_general_cons) {
            fm_warn(fn_name, "nimpl", "general constraints with Brent are not supported");
            goto cleanup;
        }
        double a, b, c;
        if (boxes[0].has_lo && boxes[0].has_hi) {
            /* Use the supplied bounds; choose interior start. */
            a = boxes[0].lo; c = boxes[0].hi; b = (a + c) * 0.5;
            /* If user gave a starting point inside, use it. */
            if (x_vec[0] > a && x_vec[0] < c) b = x_vec[0];
        } else {
            if (!fm_bracket(f_raw, binds, &opts, x_vec[0], &boxes[0], &a, &b, &c)) {
                fm_warn(fn_name, "nlnum", "bracket-finding failed");
                ok = false; goto run_done;
            }
        }
        double xm, fm;
        ok = fm_brent_min(f_raw, binds, &opts, a, b, c, &boxes[0], &xm, &fm);
        if (ok) { x_vec[0] = xm; fx_min = fm; }
    } else if (method == FM_METHOD_QUASINEWTON) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_QUASINEWTON,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_bfgs(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_CONJGRAD) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_CONJGRAD,
                                g_exprs, NULL, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_cg(f_raw, vars, n, binds, g_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else if (method == FM_METHOD_NEWTON) {
        if (has_general_cons) {
            ok = fm_run_penalty(f_raw, vars, n, binds, FM_METHOD_NEWTON,
                                g_exprs, H_exprs, x_vec, gens, ngens, boxes, &opts, &fx_min);
        } else {
            ok = fm_run_newton(f_raw, vars, n, binds, g_exprs, H_exprs, x_vec, NULL, 0, 0.0, boxes, &opts, &fx_min);
        }
    } else {
        fm_warn(fn_name, "nimpl", "method not implemented");
        ok = false;
    }
run_done:
    /* Clear temp bindings first so the variable symbol stays free during
     * Rule construction (otherwise `Rule[x, v]` would re-evaluate x to its
     * pre-call value once we restore). */
    if (binds) {
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    }
    if (ok) result_out = fm_build_result(fx_min, vars, x_vec, n);

cleanup:
    if (binds) {
        for (size_t i = 0; i < n; i++) fm_bind_restore(&binds[i]);
        free(binds);
    }
    if (g_exprs) { for (size_t i = 0; i < n; i++) expr_free(g_exprs[i]); free(g_exprs); }
    if (H_exprs) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) expr_free(H_exprs[i][j]);
            free(H_exprs[i]);
        }
        free(H_exprs);
    }
    if (gens) {
        for (size_t k = 0; k < ngens; k++) expr_free(gens[k].expr);
        free(gens);
    }
    free(vars);
    free(x_vec);
    free(boxes);
    return result_out;
}

Expr* builtin_findminimum(Expr* res) {
    return findmin_driver(res, "FindMinimum");
}

/* FindMaximum: build FindMinimum[-f, vars, opts...] internally. */
Expr* builtin_findmaximum(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn("FindMaximum", "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Negate the objective (or the f inside {f, cons}). */
    Expr* f_orig = res->data.function.args[0];
    Expr* neg_f;
    Expr* new_first;
    if (f_orig->type == EXPR_FUNCTION
        && f_orig->data.function.head->type == EXPR_SYMBOL
        && f_orig->data.function.head->data.symbol == SYM_List
        && f_orig->data.function.arg_count == 2) {
        /* Wrap inner f only. */
        Expr* inner_f = f_orig->data.function.args[0];
        Expr* cons = f_orig->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(inner_f) };
        neg_f = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
        Expr* list_args[2] = { neg_f, expr_copy(cons) };
        new_first = expr_new_function(expr_new_symbol("List"), list_args, 2);
    } else {
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(f_orig) };
        neg_f = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
        new_first = neg_f;
    }
    /* Construct synthetic FindMinimum[new_first, vars, opts...]. */
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * argc);
    new_args[0] = new_first;
    for (size_t i = 1; i < argc; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    Expr* synthetic = expr_new_function(expr_new_symbol("FindMinimum"), new_args, argc);
    free(new_args);
    /* Drive findmin directly so the diagnostic tag is FindMaximum. */
    Expr* min_result = findmin_driver(synthetic, "FindMaximum");
    expr_free(synthetic);
    if (!min_result) return NULL;
    /* min_result is {fmin, {rules}}; negate fmin. */
    if (min_result->type == EXPR_FUNCTION
        && min_result->data.function.arg_count == 2) {
        Expr* fmin_e = min_result->data.function.args[0];
        double fmin;
        if (fm_expr_to_double_real(fmin_e, &fmin)) {
            expr_free(fmin_e);
            min_result->data.function.args[0] = expr_new_real(-fmin);
        }
    }
    return min_result;
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void findmin_init(void) {
    symtab_add_builtin("FindMinimum", builtin_findminimum);
    symtab_get_def("FindMinimum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_add_builtin("FindMaximum", builtin_findmaximum);
    symtab_get_def("FindMaximum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
