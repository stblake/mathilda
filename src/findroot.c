/*
 * findroot.c
 *
 * FindRoot[f, {x, x0}] — Mathematica-compatible iterative numerical root
 * finding. Has `HoldAll, Protected` attributes; uses Block-style local
 * binding of the search variables so the user's global symbol table is
 * not perturbed during iteration.
 *
 * Supported forms
 * ---------------
 *   FindRoot[f,        {x, x0}]               Newton from a single start
 *   FindRoot[lhs==rhs, {x, x0}]               equation form (lhs - rhs)
 *   FindRoot[f,        {x, x0, x1}]           secant from two starts
 *   FindRoot[f,        {x, xstart, xmin, xmax}]  Brent bracket
 *   FindRoot[{f1,...}, {{x, x0}, {y, y0}, ...}] system Newton
 *
 * Options (Rule[…] in the trailing positions, in any order):
 *   Method            -> Automatic | "Newton" | "Secant" | "Brent"
 *   WorkingPrecision  -> MachinePrecision | digit count
 *   MaxIterations     -> integer (default 100)
 *   AccuracyGoal      -> Automatic | Infinity | digit count
 *   PrecisionGoal     -> Automatic | Infinity | digit count
 *   DampingFactor     -> positive number (default 1)
 *   Jacobian          -> user-supplied derivative / Jacobian matrix
 *   StepMonitor       -> :> expression  (held; evaluated each step)
 *   EvaluationMonitor -> :> expression  (held; evaluated each f-eval)
 *
 * Output: { var -> value } or { var -> v1, ... } for a system. Returns
 * NULL (unevaluated) when arguments cannot be interpreted as numeric.
 *
 * Memory: receives `res` owned by caller (the evaluator). Returns a
 * freshly allocated Expr* on success, or NULL on failure. Never frees
 * `res` itself. All temporary OwnValues installed for search variables
 * are removed before returning, even on the error path.
 */

#include "findroot.h"

#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include "arithmetic.h"   /* is_complex, make_complex, is_rational */
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"
#include "symtab.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ *
 *  Options + types                                                    *
 * ------------------------------------------------------------------ */

typedef enum {
    FR_METHOD_AUTOMATIC = 0,
    FR_METHOD_NEWTON,
    FR_METHOD_SECANT,
    FR_METHOD_BRENT
} FrMethod;

typedef enum {
    FR_PREC_MACHINE = 0
#ifdef USE_MPFR
    , FR_PREC_MPFR
#endif
} FrPrecMode;

typedef enum {
    FR_SPEC_BAD = 0,
    FR_SPEC_SINGLE,     /* {x, x0}            */
    FR_SPEC_TWO_START,  /* {x, x0, x1}        */
    FR_SPEC_BRACKET     /* {x, x0, xmin, xmax} */
} FrSpecKind;

typedef struct {
    FrMethod   method;
    FrPrecMode prec_mode;
    long       wp_bits;          /* MPFR bits, when prec_mode==MPFR     */
    int64_t    max_iter;         /* default 100                          */
    double     acc_goal_digits;  /* WorkingPrecision/2 by default        */
    double     prec_goal_digits; /* WorkingPrecision/2 by default        */
    double     damping;          /* default 1.0                          */
    Expr*      jacobian;         /* borrowed; user-supplied or NULL      */
    Expr*      step_monitor;     /* borrowed; or NULL                    */
    Expr*      eval_monitor;     /* borrowed; or NULL                    */
} FrOpts;

/* Per-variable temporary OwnValue snapshot, restored on exit. */
typedef struct {
    const char* name;            /* interned                            */
    Rule*       saved_own;       /* prior chain (borrowed)              */
    uint32_t    saved_attrs;
    bool        valid;
} FrVarBind;

/* ------------------------------------------------------------------ *
 *  Diagnostic helper                                                  *
 * ------------------------------------------------------------------ */

static void fr_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    fprintf(stderr, "FindRoot::%s: ", tag);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------ *
 *  Option parsing                                                     *
 * ------------------------------------------------------------------ */

static bool fr_is_known_option_name(const char* s) {
    return s == SYM_Method
        || s == SYM_WorkingPrecision
        || s == SYM_MaxIterations
        || s == SYM_AccuracyGoal
        || s == SYM_PrecisionGoal
        || s == SYM_DampingFactor
        || s == SYM_Jacobian
        || s == SYM_StepMonitor
        || s == SYM_EvaluationMonitor;
}

static bool fr_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fr_is_known_option_name(lhs->data.symbol);
}

/* Extract a double from an Integer / Real / Rational / BigInt /
 * (optionally) MPFR leaf. Returns false when the expression cannot be
 * cast to a real number. */
static bool fr_expr_to_double_real(Expr* e, double* out) {
    if (!e) return false;
    int64_t rn, rd;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;  return true;
        case EXPR_REAL:    *out = e->data.real;             return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational((Expr*)e, &rn, &rd)) {
                *out = (double)rn / (double)rd;
                return true;
            }
            return false;
        default: return false;
    }
}

/* Promote any numeric leaf to a (re, im) double-complex value. Returns
 * false when the expression is not numeric / not Complex[numeric, numeric]. */
static bool fr_expr_to_complex(Expr* e, double complex* out) {
    if (!e) return false;
    double rv;
    if (fr_expr_to_double_real(e, &rv)) { *out = rv + 0.0 * I; return true; }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        double r, i;
        if (fr_expr_to_double_real(re, &r) && fr_expr_to_double_real(im, &i)) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

/* True iff e represents a non-trivially complex value (im != 0). */
static bool fr_is_complex_value(Expr* e) {
    Expr *re, *im;
    if (!is_complex((Expr*)e, &re, &im)) return false;
    double i;
    if (!fr_expr_to_double_real(im, &i)) return true; /* symbolic imag treated as complex */
    return i != 0.0;
}

/* Build a numeric Expr from a (re, im) double pair. Pure-real
 * (im == 0) collapses to an EXPR_REAL leaf. */
static Expr* fr_expr_from_complex_d(double complex c) {
    double r = creal(c);
    double i = cimag(c);
    if (i == 0.0) return expr_new_real(r);
    return make_complex(expr_new_real(r), expr_new_real(i));
}

#ifdef USE_MPFR
/* Build an MPFR or Complex[MPFR, MPFR] expression. Pure-real (im == 0)
 * collapses to a plain MPFR leaf. */
static Expr* fr_expr_from_complex_mpfr(const mpfr_t re, const mpfr_t im) {
    if (mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
}
#endif

/* Decode the working-precision option value into (mode, bits).
 * Returns false on malformed values. */
static bool fr_parse_working_precision(Expr* val,
                                       FrPrecMode* mode, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_MachinePrecision) {
        *mode = FR_PREC_MACHINE;
        *bits = 0;
        return true;
    }
    double digits = 0.0;
    int64_t rn, rd;
    if (val->type == EXPR_INTEGER)      digits = (double)val->data.integer;
    else if (val->type == EXPR_REAL)    digits = val->data.real;
    else if (is_rational((Expr*)val, &rn, &rd)) digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;
#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        *mode = FR_PREC_MACHINE;
        *bits = 0;
    } else {
        *mode = FR_PREC_MPFR;
        *bits = numeric_digits_to_bits(digits);
    }
    return true;
#else
    *mode = FR_PREC_MACHINE;
    *bits = 0;
    return true;
#endif
}

/* Decode an AccuracyGoal / PrecisionGoal value to a digit count.
 * Special values: Automatic → -1 (caller fills in WP/2), Infinity →
 * +∞ (caller can treat as relaxed criterion). */
static bool fr_parse_goal(Expr* val, double* digits_out) {
    if (val->type == EXPR_SYMBOL) {
        if (val->data.symbol == SYM_Automatic) { *digits_out = -1.0; return true; }
        if (val->data.symbol == SYM_Infinity)  { *digits_out = INFINITY; return true; }
        return false;
    }
    return fr_expr_to_double_real(val, digits_out);
}

/* Apply a recognised Rule[opt, val] to `opts`. Returns false if the
 * value is structurally invalid for the option's role. */
static bool fr_apply_option(Expr* rule, FrOpts* opts) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
            opts->method = FR_METHOD_AUTOMATIC; return true;
        }
        if (rhs->type == EXPR_STRING) {
            const char* s = rhs->data.string;
            if (strcmp(s, "Newton") == 0) { opts->method = FR_METHOD_NEWTON; return true; }
            if (strcmp(s, "Secant") == 0) { opts->method = FR_METHOD_SECANT; return true; }
            if (strcmp(s, "Brent")  == 0) { opts->method = FR_METHOD_BRENT;  return true; }
        }
        fr_warn("badmeth", "unknown Method value");
        return false;
    }
    if (name == SYM_WorkingPrecision) {
        if (!fr_parse_working_precision(rhs, &opts->prec_mode, &opts->wp_bits)) {
            fr_warn("badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_MaxIterations) {
        if (rhs->type == EXPR_INTEGER && rhs->data.integer > 0) {
            opts->max_iter = rhs->data.integer;
            return true;
        }
        fr_warn("badopt", "MaxIterations must be a positive integer");
        return false;
    }
    if (name == SYM_AccuracyGoal)   return fr_parse_goal(rhs, &opts->acc_goal_digits);
    if (name == SYM_PrecisionGoal)  return fr_parse_goal(rhs, &opts->prec_goal_digits);
    if (name == SYM_DampingFactor) {
        double d;
        if (!fr_expr_to_double_real(rhs, &d) || d <= 0.0) {
            fr_warn("badopt", "DampingFactor must be a positive number");
            return false;
        }
        opts->damping = d;
        return true;
    }
    if (name == SYM_Jacobian)         { opts->jacobian = (Expr*)rhs;     return true; }
    if (name == SYM_StepMonitor)      { opts->step_monitor = (Expr*)rhs; return true; }
    if (name == SYM_EvaluationMonitor){ opts->eval_monitor = (Expr*)rhs; return true; }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Variable spec parsing                                              *
 * ------------------------------------------------------------------ */

/* Parse `{var, x0}` / `{var, x0, x1}` / `{var, x0, xmin, xmax}`.
 *
 *   *var_out  ← the bare variable symbol (borrowed pointer into spec).
 *   x0/x1/xmin/xmax are caller-owned freshly allocated Exprs (numeric
 *   evaluation of the user's start values). Caller frees what's set.
 *
 * Returns FR_SPEC_BAD on any malformed input. Unused outputs are NULL. */
static FrSpecKind fr_parse_var_spec(Expr* spec, Expr** var_out,
                                    Expr** x0_out, Expr** x1_out,
                                    Expr** xmin_out, Expr** xmax_out) {
    *var_out = NULL;
    *x0_out = *x1_out = *xmin_out = *xmax_out = NULL;
    if (!spec || spec->type != EXPR_FUNCTION) return FR_SPEC_BAD;
    if (spec->data.function.head->type != EXPR_SYMBOL) return FR_SPEC_BAD;
    if (spec->data.function.head->data.symbol != SYM_List) return FR_SPEC_BAD;

    size_t n = spec->data.function.arg_count;
    if (n < 2 || n > 4) return FR_SPEC_BAD;

    Expr* var = spec->data.function.args[0];
    if (var->type != EXPR_SYMBOL) return FR_SPEC_BAD;
    *var_out = var;

    /* Vector-valued start (deferred): {x, {1, 2, 3}} */
    Expr* x0_raw = spec->data.function.args[1];
    if (x0_raw->type == EXPR_FUNCTION
        && x0_raw->data.function.head->type == EXPR_SYMBOL
        && x0_raw->data.function.head->data.symbol == SYM_List) {
        fr_warn("vecvar", "vector-valued variables are not yet supported");
        return FR_SPEC_BAD;
    }

    *x0_out = eval_and_free(expr_copy(x0_raw));
    if (n == 2) return FR_SPEC_SINGLE;
    if (n == 3) {
        *x1_out = eval_and_free(expr_copy(spec->data.function.args[2]));
        return FR_SPEC_TWO_START;
    }
    /* n == 4: {x, xstart, xmin, xmax}. xstart goes in x0_out. */
    *xmin_out = eval_and_free(expr_copy(spec->data.function.args[2]));
    *xmax_out = eval_and_free(expr_copy(spec->data.function.args[3]));
    return FR_SPEC_BRACKET;
}

/* ------------------------------------------------------------------ *
 *  Variable binding (Block semantics)                                 *
 * ------------------------------------------------------------------ */

/* Snapshot the current OwnValues of `name` and clear them so that
 * a temporary value can be installed by fr_bind_set. */
static void fr_bind_snapshot(FrVarBind* b, const char* name) {
    b->name = name;
    SymbolDef* def = symtab_get_def(name);
    b->saved_own = def->own_values;
    b->saved_attrs = def->attributes;
    def->own_values = NULL;
    b->valid = true;
}

/* Replace the symbol's current single OwnValue (if any) with one that
 * binds `name -> value`. Deep-copies `value` (symtab does this); the
 * caller still owns the input. */
static void fr_bind_set(FrVarBind* b, Expr* value) {
    SymbolDef* def = symtab_get_def(b->name);
    /* Free whatever temp OwnValue we may have installed previously. */
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
    symtab_add_own_value(b->name, sym, (Expr*)value);
    expr_free(sym);
}

/* Free any temp OwnValue we installed and restore the saved chain. */
static void fr_bind_restore(FrVarBind* b) {
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
    /* Bind chain changed → eval cache must be invalidated. */
    eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Function evaluation helpers                                        *
 * ------------------------------------------------------------------ */

/* Normalise `lhs == rhs` into `Subtract[lhs, rhs]`. Pass-through for
 * non-Equal forms. Caller owns the returned tree. */
static Expr* fr_normalize_eqn(Expr* e) {
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Equal
        && e->data.function.arg_count == 2) {
        Expr* args[2] = {
            expr_copy(e->data.function.args[0]),
            expr_copy(e->data.function.args[1])
        };
        return expr_new_function(expr_new_symbol(SYM_Subtract), args, 2);
    }
    return expr_copy(e);
}

/* Run the user's EvaluationMonitor if one was supplied. The monitor is
 * an RHS expression (e.g. `s++`) that, by Mathematica convention, is
 * passed as `EvaluationMonitor :> body`. Since it arrived under HoldAll
 * inside a held RuleDelayed, we simply evaluate `body`. */
static void fr_fire_monitor(Expr* monitor) {
    if (!monitor) return;
    Expr* tmp = eval_and_free(expr_copy(monitor));
    expr_free(tmp);
}

/* Build a NumericSpec from FrOpts (used to drive numericalize after
 * each evaluation so constants like E and Pi reduce to numbers even
 * when they sit under a Power head). */
static NumericSpec fr_numeric_spec(const FrOpts* opts) {
    NumericSpec s;
#ifdef USE_MPFR
    if (opts->prec_mode == FR_PREC_MPFR) {
        s.mode = NUMERIC_MODE_MPFR;
        s.bits = opts->wp_bits;
        return s;
    }
#else
    (void)opts;
#endif
    return numeric_machine_spec();
}

/* Evaluate `f` after temporarily setting the variable bindings to the
 * given values. The result is numericalized at the requested precision
 * so that exact subexpressions like `Power[E, 1.0]` collapse to a
 * concrete numeric leaf. Returns a freshly-allocated Expr* (caller
 * frees) or NULL. */
static Expr* fr_eval_with_bindings(Expr* f, FrVarBind* binds,
                                   Expr* const* values, size_t n,
                                   Expr* eval_monitor,
                                   NumericSpec spec) {
    for (size_t i = 0; i < n; i++) {
        fr_bind_set(&binds[i], values[i]);
    }
    eval_clock_bump();
    fr_fire_monitor(eval_monitor);
    Expr* raw = eval_and_free(expr_copy(f));
    if (!raw) return NULL;
    Expr* num = numericalize(raw, spec);
    expr_free(raw);
    return num;
}

/* Construct (and evaluate) the symbolic derivative D[f, var]. Caller
 * owns the return value. Returns NULL on failure. */
static Expr* fr_compute_derivative(Expr* f, Expr* var) {
    Expr* args[2] = { expr_copy(f), expr_copy(var) };
    Expr* call = expr_new_function(expr_new_symbol(SYM_D), args, 2);
    Expr* d = eval_and_free(call);
    return d;
}

/* Build the final result {var1 -> val1, var2 -> val2, ...}.
 * `vars[i]` are borrowed; `vals[i]` are owned by this function on
 * entry (consumed into the result tree, NOT freed by the caller). */
static Expr* fr_build_rule_list(Expr* const* vars, Expr** vals, size_t n) {
    Expr** rules = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r_args[2] = { expr_copy(vars[i]), vals[i] };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
        vals[i] = NULL; /* ownership transferred */
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    return list;
}

/* ------------------------------------------------------------------ *
 *  Machine-precision scalar Newton (real)                             *
 * ------------------------------------------------------------------ */

/* Returns NULL on failure (NaN, non-numeric eval, MaxIterations with
 * no fallback). On success, returns a fresh EXPR_REAL. */
static Expr* fr_run_newton_real(Expr* f, Expr* df,
                                FrVarBind* bind, double x0,
                                const FrOpts* opts) {
    double x = x0;
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_real(x);
        Expr* arr[1] = { xv };
        Expr* fv_expr = fr_eval_with_bindings(f, bind, arr, 1,
                                              opts->eval_monitor,
                                              fr_numeric_spec(opts));
        expr_free(xv);
        if (!fv_expr) { fr_warn("nlnum", "could not evaluate f"); return NULL; }
        double fv;
        bool ok_f = fr_expr_to_double_real(fv_expr, &fv);
        expr_free(fv_expr);
        if (!ok_f) { fr_warn("nlnum", "non-real f during iteration"); return NULL; }
        /* Pre-step convergence check: if we are already at the root, stop.
         * Also avoids the f'(x*) = 0 case for repeated roots. */
        if (fabs(fv) < tol_acc) return expr_new_real(x);

        xv = expr_new_real(x);
        Expr* arr2[1] = { xv };
        Expr* dv_expr = fr_eval_with_bindings(df, bind, arr2, 1,
                                              opts->eval_monitor,
                                              fr_numeric_spec(opts));
        expr_free(xv);
        if (!dv_expr) { fr_warn("nlnum", "could not evaluate derivative"); return NULL; }
        double dv;
        bool ok_d = fr_expr_to_double_real(dv_expr, &dv);
        expr_free(dv_expr);
        if (!ok_d) { fr_warn("nlnum", "non-real derivative"); return NULL; }
        if (dv == 0.0) { fr_warn("dsing", "derivative vanished"); return NULL; }
        double step = opts->damping * fv / dv;
        if (!isfinite(step)) { fr_warn("noconv", "divergence (non-finite step)"); return NULL; }
        x -= step;
        fr_fire_monitor(opts->step_monitor);

        if (fabs(step) < (fabs(x) + 1e-300) * tol_prec)
            return expr_new_real(x);
    }
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    return expr_new_real(x);
}

/* ------------------------------------------------------------------ *
 *  Machine-precision scalar Newton (complex)                          *
 * ------------------------------------------------------------------ */

static Expr* fr_run_newton_complex(Expr* f, Expr* df,
                                   FrVarBind* bind, double complex z0,
                                   const FrOpts* opts) {
    double complex z = z0;
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* zv = fr_expr_from_complex_d(z);
        Expr* arr[1] = { zv };
        Expr* fv_expr = fr_eval_with_bindings(f, bind, arr, 1,
                                              opts->eval_monitor,
                                              fr_numeric_spec(opts));
        expr_free(zv);
        if (!fv_expr) { fr_warn("nlnum", "could not evaluate f"); return NULL; }
        double complex fv;
        bool ok_f = fr_expr_to_complex(fv_expr, &fv);
        expr_free(fv_expr);
        if (!ok_f) { fr_warn("nlnum", "non-numeric value"); return NULL; }
        if (cabs(fv) < tol_acc) return fr_expr_from_complex_d(z);

        zv = fr_expr_from_complex_d(z);
        Expr* arr2[1] = { zv };
        Expr* dv_expr = fr_eval_with_bindings(df, bind, arr2, 1,
                                              opts->eval_monitor,
                                              fr_numeric_spec(opts));
        expr_free(zv);
        if (!dv_expr) { fr_warn("nlnum", "could not evaluate derivative"); return NULL; }
        double complex dv;
        bool ok_d = fr_expr_to_complex(dv_expr, &dv);
        expr_free(dv_expr);
        if (!ok_d) { fr_warn("nlnum", "non-numeric derivative"); return NULL; }
        if (dv == 0.0) { fr_warn("dsing", "derivative vanished"); return NULL; }
        double complex step = opts->damping * fv / dv;
        if (!isfinite(creal(step)) || !isfinite(cimag(step))) {
            fr_warn("noconv", "divergence (non-finite step)");
            return NULL;
        }
        z -= step;
        double abs_step = cabs(step);
        fr_fire_monitor(opts->step_monitor);

        if (abs_step < (cabs(z) + 1e-300) * tol_prec)
            return fr_expr_from_complex_d(z);
    }
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    return fr_expr_from_complex_d(z);
}

/* ------------------------------------------------------------------ *
 *  Secant (machine-precision real)                                    *
 * ------------------------------------------------------------------ */

static Expr* fr_run_secant_real(Expr* f, FrVarBind* bind,
                                double x0, double x1,
                                const FrOpts* opts) {
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    double xa = x0, xb = x1;
    Expr* xv0 = expr_new_real(xa);
    Expr* a0[1] = { xv0 };
    Expr* fa_e = fr_eval_with_bindings(f, bind, a0, 1, opts->eval_monitor, fr_numeric_spec(opts));
    expr_free(xv0);
    double fa = 0.0;
    if (!fa_e || !fr_expr_to_double_real(fa_e, &fa)) {
        expr_free(fa_e); fr_warn("nlnum", "non-real f at start point"); return NULL;
    }
    expr_free(fa_e);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_real(xb);
        Expr* a1[1] = { xv };
        Expr* fb_e = fr_eval_with_bindings(f, bind, a1, 1, opts->eval_monitor, fr_numeric_spec(opts));
        expr_free(xv);
        double fb;
        if (!fb_e || !fr_expr_to_double_real(fb_e, &fb)) {
            expr_free(fb_e); fr_warn("nlnum", "non-real f during iteration"); return NULL;
        }
        expr_free(fb_e);
        double denom = fb - fa;
        if (denom == 0.0) { fr_warn("noconv", "secant denominator vanished"); return NULL; }
        double step = opts->damping * fb * (xb - xa) / denom;
        if (!isfinite(step)) { fr_warn("noconv", "divergence (non-finite step)"); return NULL; }
        xa = xb;
        fa = fb;
        xb -= step;
        fr_fire_monitor(opts->step_monitor);

        bool acc_ok  = (fabs(fb) < tol_acc) && (fabs(step) < tol_acc);
        bool prec_ok = (fabs(step) < (fabs(xb) + 1e-300) * tol_prec);
        if (acc_ok || prec_ok) return expr_new_real(xb);
    }
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    return expr_new_real(xb);
}

/* ------------------------------------------------------------------ *
 *  Brent's method (machine-precision real)                            *
 * ------------------------------------------------------------------ */

static Expr* fr_run_brent_real(Expr* f, FrVarBind* bind,
                               double xa, double xb, const FrOpts* opts) {
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);

    Expr* xv = expr_new_real(xa);
    Expr* a0[1] = { xv };
    Expr* fa_e = fr_eval_with_bindings(f, bind, a0, 1, opts->eval_monitor, fr_numeric_spec(opts));
    expr_free(xv);
    xv = expr_new_real(xb);
    Expr* a1[1] = { xv };
    Expr* fb_e = fr_eval_with_bindings(f, bind, a1, 1, opts->eval_monitor, fr_numeric_spec(opts));
    expr_free(xv);
    double fa, fb;
    if (!fa_e || !fb_e
        || !fr_expr_to_double_real(fa_e, &fa)
        || !fr_expr_to_double_real(fb_e, &fb)) {
        expr_free(fa_e); expr_free(fb_e);
        fr_warn("nlnum", "non-real f at bracket endpoints");
        return NULL;
    }
    expr_free(fa_e); expr_free(fb_e);

    if (fa * fb > 0.0) {
        fr_warn("brnoth", "function values at bracket endpoints have the same sign");
        return NULL;
    }
    if (fabs(fa) < fabs(fb)) { double t=xa; xa=xb; xb=t; t=fa; fa=fb; fb=t; }
    double xc = xa, fc = fa;
    double s = xb, fs = 0.0;
    double d = 0.0;
    bool used_bisection = true;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        if (fa != fc && fb != fc) {
            /* Inverse quadratic interpolation. */
            s = xa*fb*fc / ((fa-fb)*(fa-fc))
              + xb*fa*fc / ((fb-fa)*(fb-fc))
              + xc*fa*fb / ((fc-fa)*(fc-fb));
        } else {
            /* Secant. */
            s = xb - fb * (xb - xa) / (fb - fa);
        }
        double mid = (3*xa + xb) / 4.0;
        bool cond1 = (s < fmin(mid, xb)) || (s > fmax(mid, xb));
        bool cond2 = used_bisection && fabs(s - xb) >= fabs(xb - xc) / 2.0;
        bool cond3 = !used_bisection && fabs(s - xb) >= fabs(xc - d) / 2.0;
        bool cond4 = used_bisection && fabs(xb - xc) < tol_acc;
        bool cond5 = !used_bisection && fabs(xc - d) < tol_acc;
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = (xa + xb) / 2.0;
            used_bisection = true;
        } else {
            used_bisection = false;
        }
        Expr* sv = expr_new_real(s);
        Expr* aS[1] = { sv };
        Expr* fs_e = fr_eval_with_bindings(f, bind, aS, 1, opts->eval_monitor, fr_numeric_spec(opts));
        expr_free(sv);
        if (!fs_e || !fr_expr_to_double_real(fs_e, &fs)) {
            expr_free(fs_e);
            fr_warn("nlnum", "non-real f during iteration"); return NULL;
        }
        expr_free(fs_e);
        d  = xc;
        xc = xb; fc = fb;
        if (fa * fs < 0.0) { xb = s; fb = fs; }
        else               { xa = s; fa = fs; }
        if (fabs(fa) < fabs(fb)) { double t=xa; xa=xb; xb=t; t=fa; fa=fb; fb=t; }
        fr_fire_monitor(opts->step_monitor);
        if (fb == 0.0 || fabs(xb - xa) < tol_acc) return expr_new_real(xb);
    }
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    return expr_new_real(xb);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Newton at arbitrary precision (real)                               *
 * ------------------------------------------------------------------ */

static Expr* fr_run_newton_mpfr_real(Expr* f, Expr* df,
                                     FrVarBind* bind, Expr* x0_expr,
                                     const FrOpts* opts) {
    long bits = opts->wp_bits;
    mpfr_t x, fv, dv, step, tmp, tol_acc, tol_prec;
    mpfr_init2(x, bits);
    mpfr_init2(fv, bits);
    mpfr_init2(dv, bits);
    mpfr_init2(step, bits);
    mpfr_init2(tmp, bits);
    mpfr_init2(tol_acc, bits);
    mpfr_init2(tol_prec, bits);
    if (isinf(opts->acc_goal_digits)) {
        mpfr_set_zero(tol_acc, +1);
    } else {
        mpfr_set_ui(tol_acc, 10, MPFR_RNDN);
        mpfr_pow_si(tol_acc, tol_acc, -(long)opts->acc_goal_digits, MPFR_RNDN);
    }
    if (isinf(opts->prec_goal_digits)) {
        mpfr_set_zero(tol_prec, +1);
    } else {
        mpfr_set_ui(tol_prec, 10, MPFR_RNDN);
        mpfr_pow_si(tol_prec, tol_prec, -(long)opts->prec_goal_digits, MPFR_RNDN);
    }

    bool inexact_dummy;
    mpfr_t im_dummy;
    mpfr_init2(im_dummy, bits);
    if (!get_approx_mpfr(x0_expr, x, im_dummy, &inexact_dummy)) {
        mpfr_clears(x, fv, dv, step, tmp, tol_acc, tol_prec, im_dummy, (mpfr_ptr)0);
        fr_warn("nlnum", "could not extract MPFR start value"); return NULL;
    }
    mpfr_clear(im_dummy);

    Expr* result = NULL;
    bool converged = false;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_mpfr_copy(x);
        Expr* a0[1] = { xv };
        Expr* fv_e = fr_eval_with_bindings(f, bind, a0, 1,
                                           opts->eval_monitor,
                                           fr_numeric_spec(opts));
        expr_free(xv);
        if (!fv_e) {
            fr_warn("nlnum", "could not evaluate f at MPFR precision");
            goto cleanup;
        }
        mpfr_t im_f;
        mpfr_init2(im_f, bits);
        bool ok_f = get_approx_mpfr(fv_e, fv, im_f, &inexact_dummy);
        bool fim_zero = mpfr_zero_p(im_f);
        mpfr_clear(im_f);
        expr_free(fv_e);
        if (!ok_f || !fim_zero) {
            fr_warn("nlnum", "non-real f during iteration"); goto cleanup;
        }
        /* Pre-step convergence check. */
        mpfr_abs(tmp, fv, MPFR_RNDN);
        if (mpfr_cmp(tmp, tol_acc) < 0) {
            result = expr_new_mpfr_copy(x);
            converged = true;
            break;
        }

        xv = expr_new_mpfr_copy(x);
        Expr* a1[1] = { xv };
        Expr* dv_e = fr_eval_with_bindings(df, bind, a1, 1,
                                           opts->eval_monitor,
                                           fr_numeric_spec(opts));
        expr_free(xv);
        if (!dv_e) {
            fr_warn("nlnum", "could not evaluate derivative at MPFR precision");
            goto cleanup;
        }
        mpfr_t im_d;
        mpfr_init2(im_d, bits);
        bool ok_d = get_approx_mpfr(dv_e, dv, im_d, &inexact_dummy);
        bool dim_zero = mpfr_zero_p(im_d);
        mpfr_clear(im_d);
        expr_free(dv_e);
        if (!ok_d || !dim_zero) {
            fr_warn("nlnum", "non-real derivative during iteration"); goto cleanup;
        }
        if (mpfr_zero_p(dv)) { fr_warn("dsing", "derivative vanished"); goto cleanup; }
        mpfr_div(step, fv, dv, MPFR_RNDN);
        mpfr_mul_d(step, step, opts->damping, MPFR_RNDN);
        if (mpfr_nan_p(step) || mpfr_inf_p(step)) {
            fr_warn("noconv", "divergence (non-finite step)"); goto cleanup;
        }
        mpfr_sub(x, x, step, MPFR_RNDN);
        mpfr_abs(tmp, step, MPFR_RNDN);
        mpfr_t scale; mpfr_init2(scale, bits);
        mpfr_abs(scale, x, MPFR_RNDN);
        mpfr_mul(scale, scale, tol_prec, MPFR_RNDN);
        bool step_small_prec = mpfr_cmp(tmp, scale) < 0;
        mpfr_clear(scale);
        fr_fire_monitor(opts->step_monitor);
        if (step_small_prec) {
            result = expr_new_mpfr_copy(x);
            converged = true;
            break;
        }
    }
    if (!converged) {
        fr_warn("cvmit", "failed to converge within %lld iterations",
                (long long)opts->max_iter);
        result = expr_new_mpfr_copy(x);
    }
cleanup:
    mpfr_clears(x, fv, dv, step, tmp, tol_acc, tol_prec, (mpfr_ptr)0);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Newton at arbitrary precision (complex)                            *
 * ------------------------------------------------------------------ */

/* Compute z = z - damping * fv/dv, where (zr,zi), (fr_,fi), (dr,di)
 * are MPFR real/imag pairs. step_out_(r,i) receive the step magnitude
 * components. Returns false on non-finite or zero-denom. */
static bool fr_mpc_newton_step(mpfr_t zr, mpfr_t zi,
                               const mpfr_t fr_, const mpfr_t fi,
                               const mpfr_t dr, const mpfr_t di,
                               double damping,
                               mpfr_t step_r, mpfr_t step_i, long bits) {
    mpfr_t denom, num_r, num_i, t1, t2;
    mpfr_init2(denom, bits);
    mpfr_init2(num_r, bits);
    mpfr_init2(num_i, bits);
    mpfr_init2(t1, bits);
    mpfr_init2(t2, bits);
    /* denom = dr^2 + di^2 */
    mpfr_mul(t1, dr, dr, MPFR_RNDN);
    mpfr_mul(t2, di, di, MPFR_RNDN);
    mpfr_add(denom, t1, t2, MPFR_RNDN);
    if (mpfr_zero_p(denom)) {
        mpfr_clears(denom, num_r, num_i, t1, t2, (mpfr_ptr)0); return false;
    }
    /* num = f * conj(d) = (fr*dr + fi*di) + i*(fi*dr - fr*di) */
    mpfr_mul(t1, fr_, dr, MPFR_RNDN);
    mpfr_mul(t2, fi, di, MPFR_RNDN);
    mpfr_add(num_r, t1, t2, MPFR_RNDN);
    mpfr_mul(t1, fi, dr, MPFR_RNDN);
    mpfr_mul(t2, fr_, di, MPFR_RNDN);
    mpfr_sub(num_i, t1, t2, MPFR_RNDN);
    /* step = damping * num / denom */
    mpfr_div(step_r, num_r, denom, MPFR_RNDN);
    mpfr_div(step_i, num_i, denom, MPFR_RNDN);
    mpfr_mul_d(step_r, step_r, damping, MPFR_RNDN);
    mpfr_mul_d(step_i, step_i, damping, MPFR_RNDN);
    if (mpfr_nan_p(step_r) || mpfr_inf_p(step_r)
        || mpfr_nan_p(step_i) || mpfr_inf_p(step_i)) {
        mpfr_clears(denom, num_r, num_i, t1, t2, (mpfr_ptr)0); return false;
    }
    mpfr_sub(zr, zr, step_r, MPFR_RNDN);
    mpfr_sub(zi, zi, step_i, MPFR_RNDN);
    mpfr_clears(denom, num_r, num_i, t1, t2, (mpfr_ptr)0);
    return true;
}

static Expr* fr_run_newton_mpfr_complex(Expr* f, Expr* df,
                                        FrVarBind* bind, Expr* z0_expr,
                                        const FrOpts* opts) {
    long bits = opts->wp_bits;
    mpfr_t zr, zi, fr_, fi, dr, di, step_r, step_i, tol_acc, tol_prec, tmp1, tmp2;
    mpfr_init2(zr, bits);     mpfr_init2(zi, bits);
    mpfr_init2(fr_, bits);    mpfr_init2(fi, bits);
    mpfr_init2(dr, bits);     mpfr_init2(di, bits);
    mpfr_init2(step_r, bits); mpfr_init2(step_i, bits);
    mpfr_init2(tol_acc, bits); mpfr_init2(tol_prec, bits);
    mpfr_init2(tmp1, bits);   mpfr_init2(tmp2, bits);
    if (isinf(opts->acc_goal_digits)) {
        mpfr_set_zero(tol_acc, +1);
    } else {
        mpfr_set_ui(tol_acc, 10, MPFR_RNDN);
        mpfr_pow_si(tol_acc, tol_acc, -(long)opts->acc_goal_digits, MPFR_RNDN);
    }
    if (isinf(opts->prec_goal_digits)) {
        mpfr_set_zero(tol_prec, +1);
    } else {
        mpfr_set_ui(tol_prec, 10, MPFR_RNDN);
        mpfr_pow_si(tol_prec, tol_prec, -(long)opts->prec_goal_digits, MPFR_RNDN);
    }

    bool inexact_dummy;
    if (!get_approx_mpfr(z0_expr, zr, zi, &inexact_dummy)) {
        mpfr_clears(zr, zi, fr_, fi, dr, di, step_r, step_i,
                    tol_acc, tol_prec, tmp1, tmp2, (mpfr_ptr)0);
        fr_warn("nlnum", "could not extract MPFR/complex start"); return NULL;
    }
    Expr* result = NULL;
    bool converged = false;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* zv = fr_expr_from_complex_mpfr(zr, zi);
        Expr* a0[1] = { zv };
        Expr* fv_e = fr_eval_with_bindings(f, bind, a0, 1,
                                           opts->eval_monitor,
                                           fr_numeric_spec(opts));
        expr_free(zv);
        if (!fv_e) {
            fr_warn("nlnum", "could not evaluate f at MPFR precision");
            goto cleanup;
        }
        bool ok_f = get_approx_mpfr(fv_e, fr_, fi, &inexact_dummy);
        expr_free(fv_e);
        if (!ok_f) {
            fr_warn("nlnum", "non-numeric f during iteration"); goto cleanup;
        }
        /* Pre-step convergence on |f|. */
        mpfr_mul(tmp1, fr_, fr_, MPFR_RNDN);
        mpfr_mul(tmp2, fi, fi, MPFR_RNDN);
        mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
        mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);
        if (mpfr_cmp(tmp1, tol_acc) < 0) {
            result = fr_expr_from_complex_mpfr(zr, zi);
            converged = true;
            break;
        }

        zv = fr_expr_from_complex_mpfr(zr, zi);
        Expr* a1[1] = { zv };
        Expr* dv_e = fr_eval_with_bindings(df, bind, a1, 1,
                                           opts->eval_monitor,
                                           fr_numeric_spec(opts));
        expr_free(zv);
        if (!dv_e) {
            fr_warn("nlnum", "could not evaluate derivative at MPFR precision");
            goto cleanup;
        }
        bool ok_d = get_approx_mpfr(dv_e, dr, di, &inexact_dummy);
        expr_free(dv_e);
        if (!ok_d) {
            fr_warn("nlnum", "non-numeric derivative"); goto cleanup;
        }
        if (!fr_mpc_newton_step(zr, zi, fr_, fi, dr, di,
                                opts->damping, step_r, step_i, bits)) {
            fr_warn("dsing", "derivative vanished or divergence"); goto cleanup;
        }
        /* |step| = sqrt(step_r^2 + step_i^2). */
        mpfr_mul(tmp1, step_r, step_r, MPFR_RNDN);
        mpfr_mul(tmp2, step_i, step_i, MPFR_RNDN);
        mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
        mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);          /* tmp1 = |step| */
        mpfr_t mag, scale;
        mpfr_init2(mag, bits); mpfr_init2(scale, bits);
        mpfr_mul(mag, zr, zr, MPFR_RNDN);
        mpfr_mul(scale, zi, zi, MPFR_RNDN);
        mpfr_add(mag, mag, scale, MPFR_RNDN);
        mpfr_sqrt(mag, mag, MPFR_RNDN);
        mpfr_mul(scale, mag, tol_prec, MPFR_RNDN);
        bool step_small_prec = mpfr_cmp(tmp1, scale) < 0;
        mpfr_clear(mag); mpfr_clear(scale);
        fr_fire_monitor(opts->step_monitor);
        if (step_small_prec) {
            result = fr_expr_from_complex_mpfr(zr, zi);
            converged = true;
            break;
        }
    }
    if (!converged) {
        fr_warn("cvmit", "failed to converge within %lld iterations",
                (long long)opts->max_iter);
        result = fr_expr_from_complex_mpfr(zr, zi);
    }
cleanup:
    mpfr_clears(zr, zi, fr_, fi, dr, di, step_r, step_i,
                tol_acc, tol_prec, tmp1, tmp2, (mpfr_ptr)0);
    return result;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  System of equations Newton (machine precision)                     *
 * ------------------------------------------------------------------ */

/* Build a List[Subtract[lhs, rhs] | f_i] from args[0]. Caller owns. */
static Expr* fr_normalize_eqn_list(Expr* eqlist) {
    if (!eqlist || eqlist->type != EXPR_FUNCTION) return NULL;
    if (eqlist->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (eqlist->data.function.head->data.symbol != SYM_List) return NULL;
    size_t n = eqlist->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        items[i] = fr_normalize_eqn(eqlist->data.function.args[i]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return out;
}

static Expr* fr_run_newton_system_real(Expr* flist_normalized,
                                       Expr** vars, size_t n,
                                       FrVarBind* binds,
                                       double* x_vec,
                                       const FrOpts* opts) {
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);

    /* Pre-compute Jacobian symbolically: matrix[i][j] = D[f_i, var_j]. */
    Expr*** jac = malloc(sizeof(Expr**) * n);
    for (size_t i = 0; i < n; i++) {
        jac[i] = malloc(sizeof(Expr*) * n);
        Expr* fi = flist_normalized->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            if (opts->jacobian) {
                /* User-supplied Jacobian must be an n×n nested list. */
                Expr* row = opts->jacobian->data.function.args[i];
                jac[i][j] = expr_copy(row->data.function.args[j]);
            } else {
                jac[i][j] = fr_compute_derivative(fi, vars[j]);
            }
        }
    }

    Expr* result = NULL;
    bool converged = false;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr** xv = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x_vec[i]);

        /* fvec[i] = numeric f_i at current point */
        double* fvec = malloc(sizeof(double) * n);
        double max_f_pre = 0.0;
        for (size_t i = 0; i < n; i++) {
            Expr* fi = flist_normalized->data.function.args[i];
            Expr* fv_e = fr_eval_with_bindings(fi, binds, xv, n, opts->eval_monitor, fr_numeric_spec(opts));
            if (!fv_e || !fr_expr_to_double_real(fv_e, &fvec[i])) {
                expr_free(fv_e);
                fr_warn("nlnum", "non-real f_%zu during iteration", i);
                free(fvec);
                for (size_t q = 0; q < n; q++) expr_free(xv[q]);
                free(xv);
                goto cleanup;
            }
            expr_free(fv_e);
            if (fabs(fvec[i]) > max_f_pre) max_f_pre = fabs(fvec[i]);
        }
        /* Pre-step convergence check on residual norm.  Avoids touching the
         * Jacobian when we're already at the root (which would otherwise
         * trigger a "singular Jacobian" error for systems with f'(x*) = 0). */
        if (max_f_pre < tol_acc) {
            free(fvec);
            for (size_t q = 0; q < n; q++) expr_free(xv[q]);
            free(xv);
            converged = true;
            break;
        }
        /* Jmat[i][j] = D[f_i, var_j] evaluated at current point. */
        double* J = malloc(sizeof(double) * n * n);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Expr* dv_e = fr_eval_with_bindings(jac[i][j], binds, xv, n,
                                                   opts->eval_monitor,
                                                   fr_numeric_spec(opts));
                if (!dv_e || !fr_expr_to_double_real(dv_e, &J[i*n + j])) {
                    expr_free(dv_e);
                    fr_warn("nlnum", "non-real Jacobian[%zu,%zu]", i, j);
                    free(J); free(fvec);
                    for (size_t q = 0; q < n; q++) expr_free(xv[q]);
                    free(xv);
                    goto cleanup;
                }
                expr_free(dv_e);
            }
        }
        for (size_t q = 0; q < n; q++) expr_free(xv[q]);
        free(xv);

        /* Solve J · dx = fvec via Gaussian elimination with partial pivoting.
         * Then x := x - damping * dx. */
        double* A = malloc(sizeof(double) * n * (n + 1));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) A[i*(n+1) + j] = J[i*n + j];
            A[i*(n+1) + n] = fvec[i];
        }
        free(J); free(fvec);
        for (size_t col = 0; col < n; col++) {
            size_t piv = col;
            double maxv = fabs(A[col*(n+1) + col]);
            for (size_t r = col+1; r < n; r++) {
                double v = fabs(A[r*(n+1) + col]);
                if (v > maxv) { maxv = v; piv = r; }
            }
            if (maxv == 0.0) {
                free(A);
                fr_warn("dsing", "Jacobian is singular at iteration %lld", (long long)k);
                goto cleanup;
            }
            if (piv != col) {
                for (size_t j = col; j <= n; j++) {
                    double t = A[col*(n+1)+j];
                    A[col*(n+1)+j] = A[piv*(n+1)+j];
                    A[piv*(n+1)+j] = t;
                }
            }
            for (size_t r = col+1; r < n; r++) {
                double m = A[r*(n+1)+col] / A[col*(n+1)+col];
                for (size_t j = col; j <= n; j++)
                    A[r*(n+1)+j] -= m * A[col*(n+1)+j];
            }
        }
        double* dx = malloc(sizeof(double) * n);
        for (size_t i = n; i-- > 0; ) {
            double s = A[i*(n+1)+n];
            for (size_t j = i+1; j < n; j++) s -= A[i*(n+1)+j] * dx[j];
            dx[i] = s / A[i*(n+1)+i];
        }
        free(A);

        double max_step = 0.0, max_x = 0.0, max_f = 0.0;
        for (size_t i = 0; i < n; i++) {
            double step = opts->damping * dx[i];
            if (!isfinite(step)) {
                free(dx);
                fr_warn("noconv", "divergence during system Newton");
                goto cleanup;
            }
            x_vec[i] -= step;
            if (fabs(step) > max_step) max_step = fabs(step);
            if (fabs(x_vec[i]) > max_x) max_x = fabs(x_vec[i]);
        }
        /* Compute residual norm at new x */
        Expr** xv2 = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) xv2[i] = expr_new_real(x_vec[i]);
        for (size_t i = 0; i < n; i++) {
            Expr* fi = flist_normalized->data.function.args[i];
            Expr* fv_e = fr_eval_with_bindings(fi, binds, xv2, n, opts->eval_monitor, fr_numeric_spec(opts));
            double v;
            if (fv_e && fr_expr_to_double_real(fv_e, &v) && fabs(v) > max_f) {
                max_f = fabs(v);
            }
            expr_free(fv_e);
        }
        for (size_t q = 0; q < n; q++) expr_free(xv2[q]);
        free(xv2);
        free(dx);

        fr_fire_monitor(opts->step_monitor);
        bool acc_ok = (max_f < tol_acc) && (max_step < tol_acc);
        bool prec_ok = (max_step < (max_x + 1e-300) * tol_prec);
        if (acc_ok || prec_ok) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        fr_warn("cvmit", "failed to converge within %lld iterations",
                (long long)opts->max_iter);
    }
    /* Build the result list of Rules (own the new value Exprs). */
    Expr** vals = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) vals[i] = expr_new_real(x_vec[i]);
    result = fr_build_rule_list(vars, vals, n);
    free(vals);

cleanup:
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) expr_free(jac[i][j]);
        free(jac[i]);
    }
    free(jac);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

Expr* builtin_findroot(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fr_warn("argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 0 && fr_is_option_arg(res->data.function.args[pos_end - 1])) {
        pos_end--;
    }
    /* Reject any trailing Rule[sym, _] whose name is not recognised. */
    for (size_t i = pos_end; i < argc; i++) {
        if (!fr_is_option_arg(res->data.function.args[i])) {
            fr_warn("badopt", "unrecognised option in trailing position");
            return NULL;
        }
    }
    if (pos_end != 2) {
        fr_warn("argt", "needs exactly 2 positional arguments (got %zu)", pos_end);
        return NULL;
    }

    /* Default options. */
    FrOpts opts;
    opts.method = FR_METHOD_AUTOMATIC;
    opts.prec_mode = FR_PREC_MACHINE;
    opts.wp_bits = 0;
    opts.max_iter = 100;
    opts.acc_goal_digits = -1.0;
    opts.prec_goal_digits = -1.0;
    opts.damping = 1.0;
    opts.jacobian = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;

    for (size_t i = pos_end; i < argc; i++) {
        if (!fr_apply_option(res->data.function.args[i], &opts)) return NULL;
    }

    /* Fill goal defaults: half of working precision in decimal digits. */
    double wp_digits;
    if (opts.prec_mode == FR_PREC_MACHINE) wp_digits = NUMERIC_MACHINE_PRECISION_DIGITS;
    else                                    wp_digits = numeric_bits_to_digits(opts.wp_bits);
    if (opts.acc_goal_digits  < 0.0) opts.acc_goal_digits  = wp_digits / 2.0;
    if (opts.prec_goal_digits < 0.0) opts.prec_goal_digits = wp_digits / 2.0;

    Expr* f_raw    = res->data.function.args[0];
    Expr* var_spec = res->data.function.args[1];

    /* Detect system vs scalar by inspecting var_spec: a List whose
     * elements are themselves Lists ⇒ system. */
    bool is_system = false;
    if (var_spec->type == EXPR_FUNCTION
        && var_spec->data.function.head->type == EXPR_SYMBOL
        && var_spec->data.function.head->data.symbol == SYM_List
        && var_spec->data.function.arg_count > 0) {
        Expr* first = var_spec->data.function.args[0];
        if (first->type == EXPR_FUNCTION
            && first->data.function.head->type == EXPR_SYMBOL
            && first->data.function.head->data.symbol == SYM_List) {
            is_system = true;
        }
    }

    /* ============================================================== *
     *  System path                                                    *
     * ============================================================== */
    if (is_system) {
        size_t nvars = var_spec->data.function.arg_count;
        Expr** vars = malloc(sizeof(Expr*) * nvars);
        Expr** x0s  = malloc(sizeof(Expr*) * nvars);
        for (size_t i = 0; i < nvars; i++) { vars[i] = NULL; x0s[i] = NULL; }
        FrVarBind* binds = calloc(nvars, sizeof(FrVarBind));
        Expr* fnorm = NULL;
        Expr* out = NULL;
        double* xvec = NULL;

        for (size_t i = 0; i < nvars; i++) {
            Expr* sub = var_spec->data.function.args[i];
            Expr *u, *x1, *xmin, *xmax;
            FrSpecKind k = fr_parse_var_spec(sub, &u, &x0s[i], &x1, &xmin, &xmax);
            expr_free(x1); expr_free(xmin); expr_free(xmax);
            if (k != FR_SPEC_SINGLE) {
                fr_warn("ivar", "system variable spec %zu malformed (must be {var, x0})", i);
                goto sys_cleanup;
            }
            vars[i] = u;
        }
        /* Snapshot+clear OwnValues for each variable BEFORE installing
         * starting bindings, so cross-variable references in start
         * expressions were already resolved during fr_parse_var_spec. */
        for (size_t i = 0; i < nvars; i++) {
            fr_bind_snapshot(&binds[i], vars[i]->data.symbol);
        }

        fnorm = fr_normalize_eqn_list(f_raw);
        if (!fnorm
            || fnorm->data.function.arg_count != nvars) {
            fr_warn("systno", "system size %zu does not match %zu variables",
                    fnorm ? fnorm->data.function.arg_count : 0, nvars);
            goto sys_cleanup;
        }

        /* Check for complex starts (system complex path not in v1). */
        for (size_t i = 0; i < nvars; i++) {
            if (fr_is_complex_value(x0s[i])) {
                fr_warn("syscplx", "complex starts in systems are not yet supported");
                goto sys_cleanup;
            }
        }

        xvec = malloc(sizeof(double) * nvars);
        for (size_t i = 0; i < nvars; i++) {
            if (!fr_expr_to_double_real(x0s[i], &xvec[i])) {
                fr_warn("nlnum", "start value %zu is not numeric", i);
                goto sys_cleanup;
            }
        }
        out = fr_run_newton_system_real(fnorm, vars, nvars, binds, xvec, &opts);

    sys_cleanup:
        for (size_t i = 0; i < nvars; i++) fr_bind_restore(&binds[i]);
        for (size_t i = 0; i < nvars; i++) expr_free(x0s[i]);
        free(x0s);
        free(vars);
        free(binds);
        free(xvec);
        expr_free(fnorm);
        return out;
    }

    /* ============================================================== *
     *  Scalar path                                                    *
     * ============================================================== */
    Expr* var = NULL;
    Expr *x0_e = NULL, *x1_e = NULL, *xmin_e = NULL, *xmax_e = NULL;
    FrSpecKind kind = fr_parse_var_spec(var_spec, &var,
                                        &x0_e, &x1_e, &xmin_e, &xmax_e);
    if (kind == FR_SPEC_BAD) {
        fr_warn("ivar", "variable spec must be {var, x0}, {var, x0, x1}, or {var, xstart, xmin, xmax}");
        expr_free(x0_e); expr_free(x1_e); expr_free(xmin_e); expr_free(xmax_e);
        return NULL;
    }

    FrVarBind bind = {0};
    fr_bind_snapshot(&bind, var->data.symbol);

    Expr* fnorm = fr_normalize_eqn(f_raw);
    Expr* result = NULL;

    /* Choose method. */
    FrMethod method = opts.method;
    if (method == FR_METHOD_AUTOMATIC) {
        if      (kind == FR_SPEC_BRACKET)   method = FR_METHOD_BRENT;
        else if (kind == FR_SPEC_TWO_START) method = FR_METHOD_SECANT;
        else                                method = FR_METHOD_NEWTON;
    }

    /* Detect complex search. */
    bool want_complex = fr_is_complex_value(x0_e);

    if (method == FR_METHOD_BRENT) {
        /* Brent accepts a {var, x0, x1} pair as the bracket, or a
         * full {var, xstart, xmin, xmax} spec. */
        double a, b;
        if (kind == FR_SPEC_BRACKET) {
            if (!fr_expr_to_double_real(xmin_e, &a)
                || !fr_expr_to_double_real(xmax_e, &b)) {
                fr_warn("nlnum", "Brent endpoints must be real numbers");
                goto scalar_cleanup;
            }
        } else if (kind == FR_SPEC_TWO_START) {
            if (!fr_expr_to_double_real(x0_e, &a)
                || !fr_expr_to_double_real(x1_e, &b)) {
                fr_warn("nlnum", "Brent endpoints must be real numbers");
                goto scalar_cleanup;
            }
        } else {
            fr_warn("brnoth", "Brent requires {var, a, b} or {var, xstart, xmin, xmax}");
            goto scalar_cleanup;
        }
        result = fr_run_brent_real(fnorm, &bind, a, b, &opts);
    } else if (method == FR_METHOD_SECANT) {
        if (kind != FR_SPEC_TWO_START) {
            fr_warn("argt", "Secant requires a {var, x0, x1} spec");
            goto scalar_cleanup;
        }
        double x0, x1;
        if (!fr_expr_to_double_real(x0_e, &x0)
            || !fr_expr_to_double_real(x1_e, &x1)) {
            fr_warn("nlnum", "Secant starting points must be real numbers");
            goto scalar_cleanup;
        }
        result = fr_run_secant_real(fnorm, &bind, x0, x1, &opts);
    } else {
        /* Newton scalar (kind ignored beyond x0). */
        Expr* deriv = opts.jacobian ? expr_copy(opts.jacobian)
                                    : fr_compute_derivative(fnorm, var);
        if (!deriv) {
            fr_warn("nlnum", "could not compute symbolic derivative");
            goto scalar_cleanup;
        }
#ifdef USE_MPFR
        if (opts.prec_mode == FR_PREC_MPFR) {
            if (want_complex) {
                result = fr_run_newton_mpfr_complex(fnorm, deriv, &bind, x0_e, &opts);
            } else {
                result = fr_run_newton_mpfr_real(fnorm, deriv, &bind, x0_e, &opts);
            }
        } else
#endif
        {
            if (want_complex) {
                double complex z0;
                if (!fr_expr_to_complex(x0_e, &z0)) {
                    fr_warn("nlnum", "start point not numeric");
                    expr_free(deriv); goto scalar_cleanup;
                }
                result = fr_run_newton_complex(fnorm, deriv, &bind, z0, &opts);
            } else {
                double x0d;
                if (!fr_expr_to_double_real(x0_e, &x0d)) {
                    fr_warn("nlnum", "start point not numeric");
                    expr_free(deriv); goto scalar_cleanup;
                }
                result = fr_run_newton_real(fnorm, deriv, &bind, x0d, &opts);
            }
        }
        expr_free(deriv);
    }

scalar_cleanup:
    fr_bind_restore(&bind);
    expr_free(fnorm);
    expr_free(x0_e); expr_free(x1_e); expr_free(xmin_e); expr_free(xmax_e);

    if (!result) return NULL;
    /* Wrap in {var -> result}. */
    Expr* vars1[1] = { var };
    Expr* vals1[1] = { result };
    return fr_build_rule_list(vars1, vals1, 1);
}

/* ------------------------------------------------------------------ *
 *  Registration                                                       *
 * ------------------------------------------------------------------ */

void findroot_init(void) {
    symtab_add_builtin("FindRoot", builtin_findroot);
    symtab_get_def("FindRoot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
}
