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
#include "../modular.h"

#include <float.h>
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
#include "nc_accuracy.h"  /* shared AccuracyGoal/PrecisionGoal handling */
#include "compile/autocompile.h"
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
    /* Auto-compile fast path for the scalar, machine-precision, real case: the
     * compiled main function and the exact Expr pointer it was built from, so
     * fr_eval_with_bindings only takes the fast path for that function (not for
     * a symbolic derivative or a finite-difference of some other body). */
    Expr*         main_f;        /* borrowed; the function ac_f compiles */
    AutoCompiled* ac_f;          /* owned; NULL unless scalar/machine/real */
    /* Compiled MPFR fast path for a high-WorkingPrecision scalar solve: the
     * compiled main function (and, through the same funnel, its finite-difference
     * derivative) is evaluated at the working precision instead of interpreted per
     * Newton step.  NULL unless scalar/MPFR/real.  Same faithful-degrade contract. */
    AutoCompiled* ac_f_prec;
    /* Set by a machine-precision real solver that met a function or derivative
     * value outside the double range (an MPFR / BigInt value whose magnitude
     * overflows to inf); the entry point then reruns the solve in MPFR.  NULL
     * when no escalation is possible (complex search, already MPFR). */
    bool*         overflow;
    /* True for the MPFR rerun of a machine solve that overflowed. */
    bool          escalated;
} FrOpts;

/* Working precision of that rerun: a little over twice machine precision. */
#define FR_ESCALATION_BITS 106

/* Per-variable temporary OwnValue snapshot, restored on exit. */
typedef struct {
    const char* name;            /* interned                            */
    int         token;     /* binding frame (modular.h mth_bind_push) */
    bool        valid;
} FrVarBind;

/* ------------------------------------------------------------------ *
 *  Diagnostic helper                                                  *
 * ------------------------------------------------------------------ */

static void fr_warn(const char* tag, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    mth_message_v(0, "FindRoot", tag, fmt, ap);
    va_end(ap);
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
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fr_is_known_option_name(lhs->data.symbol.name);
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
    double digits;
    bool machine;
    if (!nc_parse_working_precision(val, &digits, &machine)) return false;
#ifdef USE_MPFR
    if (machine || digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        *mode = FR_PREC_MACHINE;
        *bits = 0;
    } else {
        *mode = FR_PREC_MPFR;
        *bits = numeric_digits_to_bits(digits);
    }
#else
    (void)digits; (void)machine;
    *mode = FR_PREC_MACHINE;
    *bits = 0;
#endif
    return true;
}

/* Decode an AccuracyGoal / PrecisionGoal value to a digit count.  Delegates to
 * the shared parser so the whole numerical-function family agrees on the goal
 * grammar: Automatic → NC_GOAL_AUTO (-1, caller fills in WP/2), Infinity →
 * NC_GOAL_OFF (+∞, criterion disabled), MachinePrecision →
 * NUMERIC_MACHINE_PRECISION_DIGITS, positive number → itself. */
static bool fr_parse_goal(Expr* val, double* digits_out) {
    return nc_parse_goal(val, digits_out);
}

/* Apply a recognised Rule[opt, val] to `opts`. Returns false if the
 * value is structurally invalid for the option's role. */
static bool fr_apply_option_value(const char* name, Expr* rhs, FrOpts* opts);

static bool fr_apply_option(Expr* rule, FrOpts* opts) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    /* FindRoot is HoldAll, so option values arrive unevaluated.  The scalar
     * settings are evaluated first (MaxIterations -> 2 n, WorkingPrecision ->
     * $MachinePrecision + 10); Jacobian and the monitors stay held -- they are
     * expressions in the search variables, evaluated per point. */
    if (name == SYM_Method || name == SYM_WorkingPrecision
        || name == SYM_MaxIterations || name == SYM_AccuracyGoal
        || name == SYM_PrecisionGoal || name == SYM_DampingFactor) {
        Expr* v = eval_and_free(expr_copy(rhs));
        bool ok = v && fr_apply_option_value(name, v, opts);
        expr_free(v);
        return ok;
    }
    return fr_apply_option_value(name, rhs, opts);
}

static bool fr_apply_option_value(const char* name, Expr* rhs, FrOpts* opts) {

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
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
        if (nc_parse_positive_count(rhs, &opts->max_iter)) return true;
        fr_warn("badopt", "MaxIterations must be a positive integer");
        return false;
    }
    if (name == SYM_AccuracyGoal || name == SYM_PrecisionGoal) {
        double* slot = (name == SYM_AccuracyGoal) ? &opts->acc_goal_digits
                                                  : &opts->prec_goal_digits;
        if (fr_parse_goal(rhs, slot)) return true;
        fr_warn("badopt", "invalid %s value", name);
        return false;
    }
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
    if (spec->data.function.head->data.symbol.name != SYM_List) return FR_SPEC_BAD;

    size_t n = spec->data.function.arg_count;
    if (n < 2 || n > 4) return FR_SPEC_BAD;

    Expr* var = spec->data.function.args[0];
    if (var->type != EXPR_SYMBOL) return FR_SPEC_BAD;
    *var_out = var;

    /* Vector-valued start (deferred): {x, {1, 2, 3}} */
    Expr* x0_raw = spec->data.function.args[1];
    if (x0_raw->type == EXPR_FUNCTION
        && x0_raw->data.function.head->type == EXPR_SYMBOL
        && x0_raw->data.function.head->data.symbol.name == SYM_List) {
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
    /* Unwind-safe: a TimeConstrained timeout drains the frame. */
    b->token = mth_bind_push_name(name, MTH_BIND_ATTRS);
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
    mth_bind_pop(b->token);   /* frees the temporary value */
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
        && e->data.function.head->data.symbol.name == SYM_Equal
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

/* Run the StepMonitor with the search variable bound to the NEW iterate (the
 * last binding is otherwise a finite-difference probe point or the previous
 * iterate).  `value` is borrowed. */
static void fr_fire_step_monitor(const FrOpts* opts, FrVarBind* bind, Expr* value) {
    if (!opts->step_monitor) return;
    fr_bind_set(bind, value);
    eval_clock_bump();
    fr_fire_monitor(opts->step_monitor);
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

/* Working precision of the current run in decimal digits, used to resolve the
 * shared AccuracyGoal/PrecisionGoal tolerance (nc_combined_tol). */
static double fr_wp_digits(const FrOpts* opts) {
#ifdef USE_MPFR
    if (opts->prec_mode == FR_PREC_MPFR) return numeric_bits_to_digits(opts->wp_bits);
#endif
    (void)opts;
    return NUMERIC_MACHINE_PRECISION_DIGITS;
}

/* On a MaxIterations shortfall, emit the shared "goal not met" warning for the
 * best iterate `x_mag` with error estimate `resid` (the residual |f| / step /
 * bracket width the solver last measured). Returns the best iterate unchanged;
 * this only decides whether the diagnostic fires. */
static void fr_warn_goal(const FrOpts* opts, double resid, double x_mag) {
    double tol = nc_combined_tol(opts->acc_goal_digits, opts->prec_goal_digits,
                                 x_mag, fr_wp_digits(opts));
    if (resid > tol) nc_warn_goal("FindRoot", resid, tol);
}

/* Evaluate `f` after temporarily setting the variable bindings to the
 * given values. The result is numericalized at the requested precision
 * so that exact subexpressions like `Power[E, 1.0]` collapse to a
 * concrete numeric leaf. Returns a freshly-allocated Expr* (caller
 * frees) or NULL. */
#define FR_MAXV 32   /* compiled fast path caps the variable count here */

/* Evaluate a compiled program `ac` at the real-valued `values` (machine spec,
 * matching arity).  Returns a freshly boxed numeric Expr, or NULL to signal
 * "use the interpreter" — a non-numeric argument, a non-finite result, or a spec
 * mismatch.  Shared by the scalar and system code. */
static Expr* fr_try_compiled(AutoCompiled* ac, Expr* const* values, size_t n,
                             NumericSpec spec) {
    if (!ac || spec.mode != NUMERIC_MODE_MACHINE
        || n == 0 || n > FR_MAXV || n != autocompiled_num_vars(ac)) return NULL;
    double xs[FR_MAXV];
    for (size_t i = 0; i < n; i++)
        if (!fr_expr_to_double_real(values[i], &xs[i])) return NULL;
    double _Complex z;
    if (!autocompiled_eval_complex(ac, xs, &z)
        || !isfinite(creal(z)) || !isfinite(cimag(z))) return NULL;
    return (cimag(z) == 0.0)
        ? expr_new_real(creal(z))
        : make_complex(expr_new_real(creal(z)), expr_new_real(cimag(z)));
}

/* Machine-precision compiled evaluation straight to a real double.  The system
 * Newton solver's hot path is n function values + n^2 Jacobian entries per
 * iteration; routing each through fr_try_compiled boxes the result to an Expr
 * that the caller immediately unboxes, and copies the point into a stack buffer
 * capped at FR_MAXV.  Here `x` is the current point as doubles already, so we
 * call the compiled body directly -- no Expr round-trip, and no FR_MAXV cap, so
 * a system with more than FR_MAXV variables (e.g. Broyden N=40) uses the
 * compiled path instead of falling entirely to the interpreter.  Returns false
 * (caller falls back to the interpreter) when ac is absent, the arity
 * mismatches, or the value is not a finite real. */
static bool fr_eval_compiled_double(const AutoCompiled* ac, const double* x,
                                    size_t n, double* out) {
    if (!ac || n == 0 || n != autocompiled_num_vars(ac)) return false;
    double _Complex z;
    if (!autocompiled_eval_complex(ac, x, &z)) return false;
    if (!isfinite(creal(z)) || !isfinite(cimag(z)) || cimag(z) != 0.0) return false;
    *out = creal(z);
    return true;
}

static Expr* fr_eval_with_bindings(Expr* f, FrVarBind* binds,
                                   Expr* const* values, size_t n,
                                   const FrOpts* opts) {
    NumericSpec spec = fr_numeric_spec(opts);
    /* Compiled fast path — the main function only (pointer identity).  A
     * non-finite / non-real compiled result falls through to the interpreter
     * (which supplies the complex/singular value the iteration then handles,
     * exactly as before).  Not taken under an EvaluationMonitor: the monitor
     * must run once per evaluation with the variables bound, which only the
     * interpreted path does. */
    if (f == opts->main_f && !opts->eval_monitor) {
        Expr* c = fr_try_compiled(opts->ac_f, values, n, spec);
        if (c) return c;
#ifdef USE_MPFR
        /* Compiled MPFR fast path: `values` are numeric Exprs (EXPR_MPFR); the
         * compiled body returns f at the working precision, or NULL to interpret. */
        if (opts->ac_f_prec && spec.mode == NUMERIC_MODE_MPFR
            && n == autocompiled_num_vars(opts->ac_f_prec)) {
            Expr* cm = autocompiled_eval_mpfr(opts->ac_f_prec, (const Expr* const*)values);
            if (cm) return cm;
        }
#endif
    }
    for (size_t i = 0; i < n; i++) {
        fr_bind_set(&binds[i], values[i]);
    }
    eval_clock_bump();
    fr_fire_monitor(opts->eval_monitor);
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

/* True if `e` still contains an inert Derivative[...][...] application -- i.e. D
 * could not produce a closed-form derivative (e.g. Zeta', PolyGamma order
 * derivatives, Derivative[1][f]).  Such a derivative never numericalizes, so
 * Newton must fall back to a finite-difference derivative. */
static bool fr_has_inert_derivative(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* head = e->data.function.head;
    if (head && head->type == EXPR_FUNCTION &&
        head->data.function.head &&
        head->data.function.head->type == EXPR_SYMBOL &&
        head->data.function.head->data.symbol.name == SYM_Derivative)
        return true;
    if (fr_has_inert_derivative(head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fr_has_inert_derivative(e->data.function.args[i])) return true;
    return false;
}

/* The machine finite-difference increment at a point of magnitude `mag`, given
 * the size of Newton's previous step (INFINITY before the first).  Far from
 * the root it is (|x| + 1) 1e-7, robust against a function that varies on a
 * larger scale; once the steps are shorter than that it follows them down to
 * (|x| + 1) sqrt(eps), so the difference quotient tracks the derivative at the
 * root and the iteration converges superlinearly instead of linearly (with a
 * fixed 2e-7 increment and a root at 1 - 3e-7, ChebyshevT[2001, x] took ~13
 * steps to reach its residual goal). */
static double fr_fd_increment(double mag, double prev_step) {
    double h = (mag + 1.0) * 1.0e-7;
    if (prev_step < h) {
        double floor_h = (mag + 1.0) * 1.4901161193847656e-8;   /* sqrt(eps) */
        h = prev_step > floor_h ? prev_step : floor_h;
    }
    return h;
}

/* Central finite-difference derivative f'(x) at a real machine point, used when
 * no usable symbolic derivative is available, with increment h.  Returns false
 * if f cannot be evaluated to a real number at the two probe points. */
static bool fr_fd_deriv_real(Expr* f, FrVarBind* bind, double x, double h,
                             const FrOpts* opts, double* out) {
    double fph, fmh;
    Expr* arr[1];
    Expr* e;
    arr[0] = expr_new_real(x + h);
    e = fr_eval_with_bindings(f, bind, arr, 1, opts);
    expr_free(arr[0]);
    if (!e) return false;
    bool ok = fr_expr_to_double_real(e, &fph); expr_free(e);
    if (!ok) return false;
    arr[0] = expr_new_real(x - h);
    e = fr_eval_with_bindings(f, bind, arr, 1, opts);
    expr_free(arr[0]);
    if (!e) return false;
    ok = fr_expr_to_double_real(e, &fmh); expr_free(e);
    if (!ok) return false;
    *out = (fph - fmh) / (2.0 * h);
    return true;
}

/* Central finite-difference derivative f'(z) at a complex machine point. */
static bool fr_fd_deriv_complex(Expr* f, FrVarBind* bind, double complex z, double h,
                                const FrOpts* opts, double complex* out) {
    double complex fph, fmh;
    Expr* arr[1];
    Expr* e;
    arr[0] = fr_expr_from_complex_d(z + h);
    e = fr_eval_with_bindings(f, bind, arr, 1, opts);
    expr_free(arr[0]);
    if (!e) return false;
    bool ok = fr_expr_to_complex(e, &fph); expr_free(e);
    if (!ok) return false;
    arr[0] = fr_expr_from_complex_d(z - h);
    e = fr_eval_with_bindings(f, bind, arr, 1, opts);
    expr_free(arr[0]);
    if (!e) return false;
    ok = fr_expr_to_complex(e, &fmh); expr_free(e);
    if (!ok) return false;
    *out = (fph - fmh) / (2.0 * h);
    return true;
}

#ifdef USE_MPFR
/* Central finite-difference derivative f'(x) at an MPFR real point (out must be
 * pre-init2'd).  Returns false on a non-real / non-numeric probe value. */
static bool fr_fd_deriv_mpfr_real(Expr* f, FrVarBind* bind, const mpfr_t x,
                                  const FrOpts* opts, long bits, mpfr_t out) {
    mpfr_t h, xp, fph, fmh, im;
    mpfr_inits2(bits, h, xp, fph, fmh, im, (mpfr_ptr)0);
    mpfr_abs(h, x, MPFR_RNDN);
    mpfr_add_ui(h, h, 1, MPFR_RNDN);
    /* (|x|+1) * 2^(-p/2), p the PRECISION of the working width: the guard
     * bits are no digits of f to difference against (numeric.h), and a step
     * sized on them is lost in an operand of lower precision (Zeta[s] == 1.05
     * at WorkingPrecision -> 30 differenced to an exact zero) */
    mpfr_mul_2si(h, h, -(numeric_width_precision_bits(bits) / 2), MPFR_RNDN);
    bool dummy, ok = true;
    Expr* a[1];
    Expr* e;
    mpfr_add(xp, x, h, MPFR_RNDN);
    a[0] = expr_new_mpfr_copy(xp);
    e = fr_eval_with_bindings(f, bind, a, 1, opts);
    expr_free(a[0]);
    ok = e && get_approx_mpfr(e, fph, im, &dummy) && mpfr_zero_p(im);
    if (e) expr_free(e);
    if (ok) {
        mpfr_sub(xp, x, h, MPFR_RNDN);
        a[0] = expr_new_mpfr_copy(xp);
        e = fr_eval_with_bindings(f, bind, a, 1, opts);
        expr_free(a[0]);
        ok = e && get_approx_mpfr(e, fmh, im, &dummy) && mpfr_zero_p(im);
        if (e) expr_free(e);
    }
    if (ok) {
        mpfr_sub(out, fph, fmh, MPFR_RNDN);
        mpfr_mul_2ui(h, h, 1, MPFR_RNDN);                /* 2h */
        mpfr_div(out, out, h, MPFR_RNDN);
    }
    mpfr_clears(h, xp, fph, fmh, im, (mpfr_ptr)0);
    return ok;
}

/* Central finite-difference complex derivative at an MPFR point (dr, di pre-
 * init2'd).  Probes along the real axis (f is analytic). */
static bool fr_fd_deriv_mpfr_complex(Expr* f, FrVarBind* bind,
                                     const mpfr_t zr, const mpfr_t zi,
                                     const FrOpts* opts, long bits,
                                     mpfr_t dr, mpfr_t di) {
    mpfr_t h, t, fpr, fpi, fmr, fmi;
    mpfr_inits2(bits, h, t, fpr, fpi, fmr, fmi, (mpfr_ptr)0);
    mpfr_hypot(h, zr, zi, MPFR_RNDN);
    mpfr_add_ui(h, h, 1, MPFR_RNDN);
    mpfr_mul_2si(h, h, -(numeric_width_precision_bits(bits) / 2), MPFR_RNDN);  /* as above */
    bool dummy, ok = true;
    Expr* a[1];
    Expr* e;
    mpfr_add(t, zr, h, MPFR_RNDN);
    a[0] = fr_expr_from_complex_mpfr(t, zi);
    e = fr_eval_with_bindings(f, bind, a, 1, opts);
    expr_free(a[0]);
    ok = e && get_approx_mpfr(e, fpr, fpi, &dummy);
    if (e) expr_free(e);
    if (ok) {
        mpfr_sub(t, zr, h, MPFR_RNDN);
        a[0] = fr_expr_from_complex_mpfr(t, zi);
        e = fr_eval_with_bindings(f, bind, a, 1, opts);
        expr_free(a[0]);
        ok = e && get_approx_mpfr(e, fmr, fmi, &dummy);
        if (e) expr_free(e);
    }
    if (ok) {
        mpfr_mul_2ui(h, h, 1, MPFR_RNDN);                /* 2h */
        mpfr_sub(dr, fpr, fmr, MPFR_RNDN); mpfr_div(dr, dr, h, MPFR_RNDN);
        mpfr_sub(di, fpi, fmi, MPFR_RNDN); mpfr_div(di, di, h, MPFR_RNDN);
    }
    mpfr_clears(h, t, fpr, fpi, fmr, fmi, (mpfr_ptr)0);
    return ok;
}
#endif /* USE_MPFR */

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

/* A machine solver's value check: true when `v` (a function or derivative
 * value extracted from a numeric result) is not a finite double.  Mathilda's
 * machine reals are always finite, so a non-finite `v` is an arbitrary-
 * precision value -- |L_1000(4002.)| ~ 1e867 -- whose magnitude overflowed the
 * conversion.  The solve cannot continue in doubles; record it so the entry
 * point reruns at arbitrary precision, or report divergence when it cannot. */
static bool fr_overflowed(double v, const FrOpts* opts) {
    if (isfinite(v)) return false;
    if (opts->overflow) *opts->overflow = true;
    else fr_warn("noconv", "divergence (value outside the machine range)");
    return true;
}

/* ------------------------------------------------------------------ *
 *  Machine-precision scalar Newton (real)                             *
 * ------------------------------------------------------------------ */

/* Newton's stopping rule (machine and MPFR, real and complex).
 *
 * A point x_k is accepted when BOTH tests pass there: the residual
 * |f(x_k)| < 10^-AccuracyGoal, and the step that arrived at it,
 * |x_k - x_{k-1}| < 10^-AccuracyGoal + |x_k| 10^-PrecisionGoal (or f(x_k) is
 * exactly 0).  Either test alone stops too early:
 *   - the residual alone accepts the iterate BEFORE the last useful step:
 *     Sin[x] from 3. reaches f = 2.9e-10 at 3.1415926533 (8 digits of Pi)
 *     while the next step would give all 16;
 *   - the step alone only says the iteration has slowed down, and how far x
 *     still is from the root depends on the rate.  A start within the step
 *     tolerance of the root passes it after one step, and near the root at
 *     1 - 3.08e-7 of ChebyshevT[2001, x] Newton's error constant f''/2f' ~
 *     1/(1 - x) ~ 3e6 leaves ~1e-12 after it: 3 correct digits of 1 - x of
 *     the 6 a double resolves.
 * So a passed residual test with a long arriving step takes one more step,
 * and a passed step test with a failing residual POLISHES: it keeps stepping
 * until both pass, or the step reaches the working precision's resolution
 * (|dx| <= FR_NOISE_ULPS ulps of x: further steps cannot move x), or the steps
 * stop shrinking (|dx| > |previous dx|/2: x is in the noise of f's evaluation,
 * as it is for |L_1000| ~ 1e867, whose residual never reaches any absolute
 * goal), or FR_POLISH_MAX steps have been taken.  An ordinary quadratically
 * convergent solve passes both tests one step after the step test first
 * passes, so the rule costs one extra evaluation of f. */
#define FR_NOISE_ULPS 16
#define FR_POLISH_MAX 8

/* True when a polishing step should end the iteration (see above): `step` and
 * `prev` are |dx| now and before, `noise` the resolution |x| * ulp. */
static bool fr_polish_done(double step, double prev, double noise, int polished) {
    return step <= noise || step > 0.5 * prev || polished >= FR_POLISH_MAX;
}

/* Returns NULL on failure (NaN, non-numeric eval, MaxIterations with
 * no fallback). On success, returns a fresh EXPR_REAL. */
static Expr* fr_run_newton_real(Expr* f, Expr* df,
                                FrVarBind* bind, double x0,
                                const FrOpts* opts) {
    double x = x0;
    double resid = INFINITY;    /* last |f(x)|, the goal-shortfall error estimate */
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    int polished = -1;          /* < 0: the step test has not passed yet */
    double prev_step = INFINITY;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_real(x);
        Expr* arr[1] = { xv };
        Expr* fv_expr = fr_eval_with_bindings(f, bind, arr, 1, opts);
        expr_free(xv);
        if (!fv_expr) { fr_warn("nlnum", "could not evaluate f"); return NULL; }
        double fv;
        bool ok_f = fr_expr_to_double_real(fv_expr, &fv);
        expr_free(fv_expr);
        if (!ok_f) { fr_warn("nlnum", "non-real f during iteration"); return NULL; }
        if (fr_overflowed(fv, opts)) return NULL;
        resid = fabs(fv);
        /* Both tests at x (see above); prev_step is the step that arrived. */
        bool resid_ok = fabs(fv) < tol_acc;
        if (fv == 0.0 || (resid_ok && prev_step < tol_acc + fabs(x) * tol_prec))
            return expr_new_real(x);

        double dv;
        if (df) {
            xv = expr_new_real(x);
            Expr* arr2[1] = { xv };
            Expr* dv_expr = fr_eval_with_bindings(df, bind, arr2, 1, opts);
            expr_free(xv);
            if (!dv_expr) { fr_warn("nlnum", "could not evaluate derivative"); return NULL; }
            bool ok_d = fr_expr_to_double_real(dv_expr, &dv);
            expr_free(dv_expr);
            if (!ok_d) { fr_warn("nlnum", "non-real derivative"); return NULL; }
        } else if (!fr_fd_deriv_real(f, bind, x, fr_fd_increment(fabs(x), prev_step),
                                     opts, &dv)) {
            fr_warn("nlnum", "could not evaluate finite-difference derivative");
            return NULL;
        }
        if (fr_overflowed(dv, opts)) return NULL;
        if (dv == 0.0) {
            if (polished >= 0 || resid_ok) return expr_new_real(x);   /* keep x */
            fr_warn("dsing", "derivative vanished"); return NULL;
        }
        double step = opts->damping * fv / dv;
        if (!isfinite(step)) {
            if (polished >= 0 || resid_ok) return expr_new_real(x);
            fr_warn("noconv", "divergence (non-finite step)"); return NULL;
        }
        x -= step;
        if (opts->step_monitor) {
            Expr* xv_new = expr_new_real(x);
            fr_fire_step_monitor(opts, bind, xv_new);
            expr_free(xv_new);
        }

        double astep = fabs(step);
        if (polished >= 0) {
            polished++;
            if (fr_polish_done(astep, prev_step, FR_NOISE_ULPS * DBL_EPSILON * fabs(x),
                               polished))
                return expr_new_real(x);
        } else if (astep < tol_acc + fabs(x) * tol_prec) {
            polished = 0;                         /* step test passed: polish */
        }
        prev_step = astep;
    }
    if (polished >= 0) return expr_new_real(x);   /* the step test had passed */
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    fr_warn_goal(opts, resid, fabs(x));
    return expr_new_real(x);
}

/* ------------------------------------------------------------------ *
 *  Machine-precision scalar Newton (complex)                          *
 * ------------------------------------------------------------------ */

static Expr* fr_run_newton_complex(Expr* f, Expr* df,
                                   FrVarBind* bind, double complex z0,
                                   const FrOpts* opts) {
    double complex z = z0;
    double resid = INFINITY;    /* last |f(z)|, the goal-shortfall error estimate */
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    int polished = -1;          /* see fr_polish_done */
    double prev_step = INFINITY;

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* zv = fr_expr_from_complex_d(z);
        Expr* arr[1] = { zv };
        Expr* fv_expr = fr_eval_with_bindings(f, bind, arr, 1, opts);
        expr_free(zv);
        if (!fv_expr) { fr_warn("nlnum", "could not evaluate f"); return NULL; }
        double complex fv;
        bool ok_f = fr_expr_to_complex(fv_expr, &fv);
        expr_free(fv_expr);
        if (!ok_f) { fr_warn("nlnum", "non-numeric value"); return NULL; }
        resid = cabs(fv);
        bool resid_ok = cabs(fv) < tol_acc;       /* both tests, as in the real case */
        if (fv == 0.0 || (resid_ok && prev_step < tol_acc + cabs(z) * tol_prec))
            return fr_expr_from_complex_d(z);

        double complex dv;
        if (df) {
            zv = fr_expr_from_complex_d(z);
            Expr* arr2[1] = { zv };
            Expr* dv_expr = fr_eval_with_bindings(df, bind, arr2, 1, opts);
            expr_free(zv);
            if (!dv_expr) { fr_warn("nlnum", "could not evaluate derivative"); return NULL; }
            bool ok_d = fr_expr_to_complex(dv_expr, &dv);
            expr_free(dv_expr);
            if (!ok_d) { fr_warn("nlnum", "non-numeric derivative"); return NULL; }
        } else if (!fr_fd_deriv_complex(f, bind, z, fr_fd_increment(cabs(z), prev_step),
                                        opts, &dv)) {
            fr_warn("nlnum", "could not evaluate finite-difference derivative");
            return NULL;
        }
        if (dv == 0.0) {
            if (polished >= 0 || resid_ok) return fr_expr_from_complex_d(z);
            fr_warn("dsing", "derivative vanished"); return NULL;
        }
        double complex step = opts->damping * fv / dv;
        if (!isfinite(creal(step)) || !isfinite(cimag(step))) {
            if (polished >= 0 || resid_ok) return fr_expr_from_complex_d(z);
            fr_warn("noconv", "divergence (non-finite step)");
            return NULL;
        }
        z -= step;
        double abs_step = cabs(step);
        if (opts->step_monitor) {
            Expr* zv_new = fr_expr_from_complex_d(z);
            fr_fire_step_monitor(opts, bind, zv_new);
            expr_free(zv_new);
        }

        if (polished >= 0) {
            polished++;
            if (fr_polish_done(abs_step, prev_step, FR_NOISE_ULPS * DBL_EPSILON * cabs(z),
                               polished))
                return fr_expr_from_complex_d(z);
        } else if (abs_step < tol_acc + cabs(z) * tol_prec) {
            polished = 0;
        }
        prev_step = abs_step;
    }
    if (polished >= 0) return fr_expr_from_complex_d(z);
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    fr_warn_goal(opts, resid, cabs(z));
    return fr_expr_from_complex_d(z);
}

/* ------------------------------------------------------------------ *
 *  Brent's method (machine-precision real)                            *
 * ------------------------------------------------------------------ */

/* Evaluate f at the machine point x into *out.  False (after a message, or
 * after recording an overflow for escalation) when the value is not a finite
 * real double. */
static bool fr_eval_real_d(Expr* f, FrVarBind* bind, double x,
                           const FrOpts* opts, double* out) {
    Expr* a[1] = { expr_new_real(x) };
    Expr* e = fr_eval_with_bindings(f, bind, a, 1, opts);
    expr_free(a[0]);
    bool ok = e && fr_expr_to_double_real(e, out);
    expr_free(e);
    if (!ok) { fr_warn("nlnum", "non-real f during iteration"); return false; }
    return !fr_overflowed(*out, opts);
}

/* Brent's method in doubles from a bracket [a, b] with fa, fb already
 * evaluated and of opposite signs: the machine twin of fr_brent_mpfr_core,
 * with the same steps and the same stopping rules (see there), the working
 * precision's resolution being DBL_EPSILON.  `iters` caps the evaluations. */
static Expr* fr_brent_real_core(Expr* f, FrVarBind* bind,
                                double a, double fa, double b, double fb,
                                int64_t iters, const FrOpts* opts) {
    double tol_acc  = isinf(opts->acc_goal_digits)  ? 0.0 : pow(10.0, -opts->acc_goal_digits);
    double tol_prec = isinf(opts->prec_goal_digits) ? 0.0 : pow(10.0, -opts->prec_goal_digits);
    double c = b, fc = fb, d = b - a, e = d;
    for (int64_t k = 0; ; k++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa; d = b - a; e = d;
        }
        if (fabs(fc) < fabs(fb)) {
            a = b; fa = fb; b = c; fb = fc; c = a; fc = fa;
        }
        double goal = tol_acc + fabs(b) * tol_prec;
        double tol1 = 2.0 * DBL_EPSILON * fabs(b);
        double xm = 0.5 * (c - b);
        if (fb == 0.0 || fabs(xm) <= tol1 + 0.5 * goal) return expr_new_real(b);
        if (k >= iters) {
            fr_warn("cvmit", "failed to converge within %lld iterations",
                    (long long)opts->max_iter);
            fr_warn_goal(opts, 2.0 * fabs(xm), fabs(b));
            return expr_new_real(b);
        }
        bool interpolate = false;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            double p, q, r, sr = fb / fa;
            if (a == c) {                                   /* secant            */
                p = 2.0 * xm * sr;
                q = 1.0 - sr;
            } else {                                        /* inverse quadratic */
                q = fa / fc;
                r = fb / fc;
                p = sr * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (sr - 1.0);
            }
            if (p > 0.0) q = -q;
            p = fabs(p);
            double min1 = 3.0 * xm * q - fabs(tol1 * q), min2 = fabs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2) && q != 0.0) {
                e = d;
                d = p / q;
                interpolate = true;
            }
        }
        if (!interpolate) { d = xm; e = d; }                /* bisection         */
        bool final_step = interpolate && fabs(d) <= goal && fabs(fb) <= tol_acc;
        a = b; fa = fb;
        if (fabs(d) > tol1) b += d;
        else                b += (xm > 0.0 ? tol1 : -tol1);
        if (!fr_eval_real_d(f, bind, b, opts, &fb)) return NULL;
        if (final_step) return expr_new_real(b);
        fr_fire_monitor(opts->step_monitor);
    }
}

static Expr* fr_run_brent_real(Expr* f, FrVarBind* bind,
                               double xa, double xb, const FrOpts* opts) {
    double fa, fb;
    if (!fr_eval_real_d(f, bind, xa, opts, &fa)
        || !fr_eval_real_d(f, bind, xb, opts, &fb)) return NULL;
    if (fa * fb > 0.0) {
        fr_warn("brnoth", "function values at bracket endpoints have the same sign");
        return NULL;
    }
    return fr_brent_real_core(f, bind, xa, fa, xb, fb, opts->max_iter, opts);
}

/* ------------------------------------------------------------------ *
 *  Secant (machine-precision real)                                    *
 * ------------------------------------------------------------------ */

/* The secant iteration, safeguarded by bracketing exactly as the MPFR one
 * (fr_run_secant_mpfr_real): a sign change between consecutive iterates hands
 * the bracket to Brent's method. */
static Expr* fr_run_secant_real(Expr* f, FrVarBind* bind,
                                double x0, double x1,
                                const FrOpts* opts) {
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    double xa = x0, xb = x1;
    double resid = INFINITY;    /* last |f(xb)|, the goal-shortfall error estimate */
    Expr* xv0 = expr_new_real(xa);
    Expr* a0[1] = { xv0 };
    Expr* fa_e = fr_eval_with_bindings(f, bind, a0, 1, opts);
    expr_free(xv0);
    double fa = 0.0;
    if (!fa_e || !fr_expr_to_double_real(fa_e, &fa)) {
        expr_free(fa_e); fr_warn("nlnum", "non-real f at start point"); return NULL;
    }
    expr_free(fa_e);
    if (fr_overflowed(fa, opts)) return NULL;
    if (fa == 0.0) return expr_new_real(xa);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_real(xb);
        Expr* a1[1] = { xv };
        Expr* fb_e = fr_eval_with_bindings(f, bind, a1, 1, opts);
        expr_free(xv);
        double fb;
        if (!fb_e || !fr_expr_to_double_real(fb_e, &fb)) {
            expr_free(fb_e); fr_warn("nlnum", "non-real f during iteration"); return NULL;
        }
        expr_free(fb_e);
        if (fr_overflowed(fb, opts)) return NULL;
        resid = fabs(fb);
        if (fb == 0.0) return expr_new_real(xb);
        /* Bracketed (see fr_run_secant_mpfr_real): continue by Brent. */
        if ((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0))
            return fr_brent_real_core(f, bind, xa, fa, xb, fb, opts->max_iter - k, opts);
        double denom = fb - fa;
        if (denom == 0.0) {
            if (xa == xb) return expr_new_real(xb);
            fr_warn("noconv", "secant denominator vanished");
            return NULL;
        }
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
    fr_warn_goal(opts, resid, fabs(xb));
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
    double resid = INFINITY;    /* last |f(x)|, the goal-shortfall error estimate */
    int polished = -1;          /* see fr_polish_done */
    mpfr_t prev_step;
    mpfr_init2(prev_step, bits);
    mpfr_set_inf(prev_step, +1);

    /* An escalated solve (the machine run overflowed) is by construction one
     * where f takes extreme values, and there an automatic symbolic derivative
     * is often numerically useless: D[LaguerreL[1000, x], x] is the expanded
     * degree-999 power-basis polynomial, whose terms near x = 4000 reach 1e1120
     * while their sum is 1e867 -- 250 digits of cancellation.  Check it once
     * against a central difference at the start point and use the difference
     * quotient throughout when the two disagree. */
    if (opts->escalated && df && !opts->jacobian) {
        mpfr_t fd;
        mpfr_init2(fd, bits);
        bool fd_ok = fr_fd_deriv_mpfr_real(f, bind, x, opts, bits, fd);
        bool sym_ok = false;
        if (fd_ok) {
            Expr* xv = expr_new_mpfr_copy(x);
            Expr* a1[1] = { xv };
            Expr* dv_e = fr_eval_with_bindings(df, bind, a1, 1, opts);
            expr_free(xv);
            mpfr_t im_d;
            mpfr_init2(im_d, bits);
            if (dv_e && get_approx_mpfr(dv_e, dv, im_d, &inexact_dummy)
                && mpfr_zero_p(im_d) && mpfr_number_p(dv)) {
                mpfr_sub(tmp, dv, fd, MPFR_RNDN);       /* |sym - fd| <= 1e-6 |fd| */
                mpfr_abs(tmp, tmp, MPFR_RNDN);
                mpfr_abs(step, fd, MPFR_RNDN);
                mpfr_mul_d(step, step, 1e-6, MPFR_RNDN);
                sym_ok = mpfr_cmp(tmp, step) <= 0;
            }
            mpfr_clear(im_d);
            expr_free(dv_e);
        }
        mpfr_clear(fd);
        if (fd_ok && !sym_ok) df = NULL;
    }

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* xv = expr_new_mpfr_copy(x);
        Expr* a0[1] = { xv };
        Expr* fv_e = fr_eval_with_bindings(f, bind, a0, 1, opts);
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
        /* Both tests at x (see fr_polish_done); prev_step arrived here. */
        mpfr_abs(tmp, fv, MPFR_RNDN);
        resid = mpfr_get_d(tmp, MPFR_RNDN);
        bool resid_ok = mpfr_cmp(tmp, tol_acc) < 0;
        if (resid_ok || mpfr_zero_p(fv)) {
            mpfr_abs(step, x, MPFR_RNDN);
            mpfr_mul(step, step, tol_prec, MPFR_RNDN);
            mpfr_add(step, step, tol_acc, MPFR_RNDN);       /* step tolerance */
            if (mpfr_zero_p(fv) || mpfr_cmp(prev_step, step) < 0) {
                result = expr_new_mpfr_copy(x);
                converged = true;
                break;
            }
        }

        if (df) {
            xv = expr_new_mpfr_copy(x);
            Expr* a1[1] = { xv };
            Expr* dv_e = fr_eval_with_bindings(df, bind, a1, 1, opts);
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
        } else if (!fr_fd_deriv_mpfr_real(f, bind, x, opts, bits, dv)) {
            fr_warn("nlnum", "could not evaluate finite-difference derivative");
            goto cleanup;
        }
        if (mpfr_zero_p(dv)) {
            if (polished >= 0 || resid_ok) {
                result = expr_new_mpfr_copy(x); converged = true; break;
            }
            fr_warn("dsing", "derivative vanished"); goto cleanup;
        }
        mpfr_div(step, fv, dv, MPFR_RNDN);
        mpfr_mul_d(step, step, opts->damping, MPFR_RNDN);
        if (mpfr_nan_p(step) || mpfr_inf_p(step)) {
            if (polished >= 0 || resid_ok) {
                result = expr_new_mpfr_copy(x); converged = true; break;
            }
            fr_warn("noconv", "divergence (non-finite step)"); goto cleanup;
        }
        mpfr_sub(x, x, step, MPFR_RNDN);
        mpfr_abs(tmp, step, MPFR_RNDN);
        mpfr_t scale; mpfr_init2(scale, bits);
        mpfr_abs(scale, x, MPFR_RNDN);
        bool stop = false;
        if (polished >= 0) {
            /* fr_polish_done in MPFR: at the resolution, stalled, or capped. */
            /* 16 ulps of x -- of the double it is returned as, for an escalated
             * machine solve: bits beyond 53 are rounded away. */
            polished++;
            mpfr_mul_2si(scale, scale, 4 - (opts->escalated ? 53 : bits), MPFR_RNDN);
            mpfr_mul_2si(prev_step, prev_step, -1, MPFR_RNDN);    /* prev / 2     */
            stop = mpfr_cmp(tmp, scale) <= 0 || mpfr_cmp(tmp, prev_step) > 0
                   || polished >= FR_POLISH_MAX;
        } else {
            mpfr_mul(scale, scale, tol_prec, MPFR_RNDN);
            mpfr_add(scale, scale, tol_acc, MPFR_RNDN);
            if (mpfr_cmp(tmp, scale) < 0) polished = 0;          /* step test: polish */
        }
        mpfr_set(prev_step, tmp, MPFR_RNDN);
        mpfr_clear(scale);
        if (opts->step_monitor) {
            Expr* xv_new = expr_new_mpfr_copy(x);
            fr_fire_step_monitor(opts, bind, xv_new);
            expr_free(xv_new);
        }
        if (stop) {
            result = expr_new_mpfr_copy(x);
            converged = true;
            break;
        }
    }
    if (!converged && polished >= 0) {             /* the step test had passed */
        result = expr_new_mpfr_copy(x);
        converged = true;
    }
    if (!converged) {
        fr_warn("cvmit", "failed to converge within %lld iterations",
                (long long)opts->max_iter);
        result = expr_new_mpfr_copy(x);
        mpfr_abs(tmp, x, MPFR_RNDN);
        fr_warn_goal(opts, resid, mpfr_get_d(tmp, MPFR_RNDN));
    }
cleanup:
    mpfr_clears(x, fv, dv, step, tmp, tol_acc, tol_prec, prev_step, (mpfr_ptr)0);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Secant at arbitrary precision (real)                               *
 * ------------------------------------------------------------------ */

/* Set the absolute / relative tolerances 10^-AccuracyGoal, 10^-PrecisionGoal
 * (an Infinity goal disables its criterion: tolerance 0). */
static void fr_mpfr_goal_tols(const FrOpts* opts, mpfr_t tol_acc, mpfr_t tol_prec) {
    if (isinf(opts->acc_goal_digits)) mpfr_set_zero(tol_acc, +1);
    else {
        mpfr_set_ui(tol_acc, 10, MPFR_RNDN);
        mpfr_pow_si(tol_acc, tol_acc, -(long)opts->acc_goal_digits, MPFR_RNDN);
    }
    if (isinf(opts->prec_goal_digits)) mpfr_set_zero(tol_prec, +1);
    else {
        mpfr_set_ui(tol_prec, 10, MPFR_RNDN);
        mpfr_pow_si(tol_prec, tol_prec, -(long)opts->prec_goal_digits, MPFR_RNDN);
    }
}

/* Evaluate f at the real MPFR point x into `out` (pre-init2'd).  False when
 * the value is not a real number. */
static bool fr_eval_mpfr_real(Expr* f, FrVarBind* bind, const mpfr_t x,
                              const FrOpts* opts, long bits, mpfr_t out) {
    Expr* a[1] = { expr_new_mpfr_copy(x) };
    Expr* e = fr_eval_with_bindings(f, bind, a, 1, opts);
    expr_free(a[0]);
    if (!e) return false;
    mpfr_t im;
    mpfr_init2(im, bits);
    bool dummy;
    bool ok = get_approx_mpfr(e, out, im, &dummy) && mpfr_zero_p(im)
              && mpfr_number_p(out);
    mpfr_clear(im);
    expr_free(e);
    return ok;
}

/* Brent's method at the working precision, from a bracket [a, b] whose
 * function values fa, fb (already evaluated) have opposite signs.  The
 * classical algorithm (inverse quadratic interpolation / secant, each step
 * accepted only while it stays inside the bracket and shrinks it fast enough,
 * else bisection), so it can never leave the bracket and always converges to a
 * root inside it.  With tol = 10^-AccuracyGoal + |b| 10^-PrecisionGoal it stops
 *   - when the residual |f(b)| is under 10^-AccuracyGoal and the next
 *     interpolation step is shorter than tol: that final step is taken and the
 *     point it reaches returned.  The iteration is then in its superlinear
 *     phase, so that point is far more accurate than tol (a bracket-width test
 *     alone stops at ~tol: half the working precision under the default
 *     goals);
 *   - when the bracket half-width falls under 2^(1-bits) |b| + tol/2, at the
 *     goal or at the working precision, whichever comes first -- an
 *     absolute-only test never terminates for an AccuracyGoal beyond what the
 *     working precision can resolve at |b|.
 * `iters` caps the evaluations.  Returns the best iterate, or NULL (with a
 * message) on failure. */
static Expr* fr_brent_mpfr_core(Expr* f, FrVarBind* bind,
                                const mpfr_t a0, const mpfr_t fa0,
                                const mpfr_t b0, const mpfr_t fb0,
                                int64_t iters, const FrOpts* opts) {
    long bits = opts->wp_bits;
    mpfr_t a, b, c, fa, fb, fc, d, e, p, q, r, s, tol1, xm, t1, t2,
           tol_acc, tol_prec, eps, goal;
    mpfr_inits2(bits, a, b, c, fa, fb, fc, d, e, p, q, r, s, tol1, xm, t1, t2,
                tol_acc, tol_prec, eps, goal, (mpfr_ptr)0);
    fr_mpfr_goal_tols(opts, tol_acc, tol_prec);
    mpfr_set_ui_2exp(eps, 1, 1 - bits, MPFR_RNDN);
    mpfr_set(a, a0, MPFR_RNDN);  mpfr_set(fa, fa0, MPFR_RNDN);
    mpfr_set(b, b0, MPFR_RNDN);  mpfr_set(fb, fb0, MPFR_RNDN);
    mpfr_set(c, b, MPFR_RNDN);   mpfr_set(fc, fb, MPFR_RNDN);
    mpfr_sub(d, b, a, MPFR_RNDN); mpfr_set(e, d, MPFR_RNDN);
    Expr* result = NULL;

    for (int64_t k = 0; ; k++) {
        if (mpfr_sgn(fb) * mpfr_sgn(fc) > 0) {          /* re-bracket with a  */
            mpfr_set(c, a, MPFR_RNDN); mpfr_set(fc, fa, MPFR_RNDN);
            mpfr_sub(d, b, a, MPFR_RNDN); mpfr_set(e, d, MPFR_RNDN);
        }
        if (mpfr_cmpabs(fc, fb) < 0) {                  /* b = best iterate   */
            mpfr_set(a, b, MPFR_RNDN); mpfr_set(fa, fb, MPFR_RNDN);
            mpfr_swap(b, c);           mpfr_swap(fb, fc);
            mpfr_set(c, a, MPFR_RNDN); mpfr_set(fc, fa, MPFR_RNDN);
        }
        /* goal = tol_acc + |b| tol_prec.  tol1 = eps |b| is the resolution of
         * the working precision: the shortest step taken and the scale of the
         * interpolation acceptance tests.  It deliberately excludes the goal --
         * a goal-sized minimum step would push the final, superlinearly
         * accurate interpolation step off its target by ~goal. */
        mpfr_abs(t1, b, MPFR_RNDN);
        mpfr_mul(goal, t1, tol_prec, MPFR_RNDN);
        mpfr_add(goal, goal, tol_acc, MPFR_RNDN);
        mpfr_mul(tol1, t1, eps, MPFR_RNDN);
        mpfr_sub(xm, c, b, MPFR_RNDN);
        mpfr_div_2ui(xm, xm, 1, MPFR_RNDN);
        mpfr_div_2ui(t2, goal, 1, MPFR_RNDN);
        mpfr_add(t2, t2, tol1, MPFR_RNDN);               /* tol1 + goal / 2 */
        if (mpfr_zero_p(fb) || mpfr_cmpabs(xm, t2) <= 0) {
            result = expr_new_mpfr_copy(b);
            break;
        }
        if (k >= iters) {
            fr_warn("cvmit", "failed to converge within %lld iterations",
                    (long long)opts->max_iter);
            mpfr_abs(t1, b, MPFR_RNDN);
            mpfr_abs(t2, xm, MPFR_RNDN);
            fr_warn_goal(opts, 2.0 * mpfr_get_d(t2, MPFR_RNDN), mpfr_get_d(t1, MPFR_RNDN));
            result = expr_new_mpfr_copy(b);
            break;
        }
        bool interpolate = false;
        if (mpfr_cmpabs(e, tol1) >= 0 && mpfr_cmpabs(fa, fb) > 0) {
            mpfr_div(s, fb, fa, MPFR_RNDN);
            if (mpfr_equal_p(a, c)) {                   /* secant             */
                mpfr_mul(p, xm, s, MPFR_RNDN); mpfr_mul_2ui(p, p, 1, MPFR_RNDN);
                mpfr_ui_sub(q, 1, s, MPFR_RNDN);
            } else {                                    /* inverse quadratic  */
                mpfr_div(q, fa, fc, MPFR_RNDN);
                mpfr_div(r, fb, fc, MPFR_RNDN);
                /* p = s (2 xm q (q - r) - (b - a)(r - 1)) */
                mpfr_sub(t1, q, r, MPFR_RNDN);
                mpfr_mul(t1, t1, q, MPFR_RNDN);
                mpfr_mul(t1, t1, xm, MPFR_RNDN);
                mpfr_mul_2ui(t1, t1, 1, MPFR_RNDN);
                mpfr_sub(t2, b, a, MPFR_RNDN);
                mpfr_sub_ui(p, r, 1, MPFR_RNDN);
                mpfr_mul(t2, t2, p, MPFR_RNDN);
                mpfr_sub(t1, t1, t2, MPFR_RNDN);
                mpfr_mul(p, s, t1, MPFR_RNDN);
                /* q = (q - 1)(r - 1)(s - 1) */
                mpfr_sub_ui(q, q, 1, MPFR_RNDN);
                mpfr_sub_ui(t1, r, 1, MPFR_RNDN);
                mpfr_mul(q, q, t1, MPFR_RNDN);
                mpfr_sub_ui(t1, s, 1, MPFR_RNDN);
                mpfr_mul(q, q, t1, MPFR_RNDN);
            }
            if (mpfr_sgn(p) > 0) mpfr_neg(q, q, MPFR_RNDN);
            mpfr_abs(p, p, MPFR_RNDN);
            /* accept if 2p < min(3 xm q - |tol1 q|, |e q|) */
            mpfr_mul(t1, xm, q, MPFR_RNDN); mpfr_mul_ui(t1, t1, 3, MPFR_RNDN);
            mpfr_mul(t2, tol1, q, MPFR_RNDN); mpfr_abs(t2, t2, MPFR_RNDN);
            mpfr_sub(t1, t1, t2, MPFR_RNDN);
            mpfr_mul(t2, e, q, MPFR_RNDN); mpfr_abs(t2, t2, MPFR_RNDN);
            if (mpfr_cmp(t2, t1) < 0) mpfr_set(t1, t2, MPFR_RNDN);
            mpfr_mul_2ui(t2, p, 1, MPFR_RNDN);
            if (mpfr_cmp(t2, t1) < 0 && !mpfr_zero_p(q)) {
                mpfr_set(e, d, MPFR_RNDN);
                mpfr_div(d, p, q, MPFR_RNDN);
                interpolate = true;
            }
        }
        if (!interpolate) {                              /* bisection          */
            mpfr_set(d, xm, MPFR_RNDN);
            mpfr_set(e, d, MPFR_RNDN);
        }
        bool final_step = interpolate && mpfr_cmpabs(d, goal) <= 0
                          && mpfr_cmpabs(fb, tol_acc) <= 0;
        mpfr_set(a, b, MPFR_RNDN); mpfr_set(fa, fb, MPFR_RNDN);
        if (mpfr_cmpabs(d, tol1) > 0) mpfr_add(b, b, d, MPFR_RNDN);
        else if (mpfr_sgn(xm) > 0)    mpfr_add(b, b, tol1, MPFR_RNDN);
        else                          mpfr_sub(b, b, tol1, MPFR_RNDN);
        if (!fr_eval_mpfr_real(f, bind, b, opts, bits, fb)) {
            fr_warn("nlnum", "non-real f during iteration");
            break;
        }
        if (final_step) { result = expr_new_mpfr_copy(b); break; }
        fr_fire_monitor(opts->step_monitor);
    }
    mpfr_clears(a, b, c, fa, fb, fc, d, e, p, q, r, s, tol1, xm, t1, t2,
                tol_acc, tol_prec, eps, goal, (mpfr_ptr)0);
    return result;
}

/* Method -> "Brent" at the working precision: evaluate the two bracket
 * endpoints and run fr_brent_mpfr_core. */
static Expr* fr_run_brent_mpfr_real(Expr* f, FrVarBind* bind,
                                    Expr* a_expr, Expr* b_expr,
                                    const FrOpts* opts) {
    long bits = opts->wp_bits;
    mpfr_t a, b, fa, fb, im;
    mpfr_inits2(bits, a, b, fa, fb, im, (mpfr_ptr)0);
    bool dummy;
    Expr* result = NULL;
    if (!get_approx_mpfr(a_expr, a, im, &dummy) || !mpfr_zero_p(im)
        || !get_approx_mpfr(b_expr, b, im, &dummy) || !mpfr_zero_p(im)) {
        fr_warn("nlnum", "Brent endpoints must be real numbers");
    } else if (!fr_eval_mpfr_real(f, bind, a, opts, bits, fa)
               || !fr_eval_mpfr_real(f, bind, b, opts, bits, fb)) {
        fr_warn("nlnum", "non-real f at bracket endpoints");
    } else if (mpfr_sgn(fa) * mpfr_sgn(fb) > 0) {
        fr_warn("brnoth", "function values at bracket endpoints have the same sign");
    } else {
        result = fr_brent_mpfr_core(f, bind, a, fa, b, fb, opts->max_iter, opts);
    }
    mpfr_clears(a, b, fa, fb, im, (mpfr_ptr)0);
    return result;
}

/* The secant iteration x_{k+1} = x_k - f_k (x_k - x_{k-1}) / (f_k - f_{k-1}),
 * safeguarded by bracketing: as soon as two consecutive iterates have function
 * values of opposite sign (the two starting values included), a root is known
 * to lie between them and the search continues by Brent's method on that
 * bracket, which cannot leave it.  An unguarded secant step from a bracket can
 * extrapolate across many sign changes and land on a far-away root -- or on
 * none.  Carried entirely in MPFR at the working precision.  Every quantity -- the
 * function values, their difference, the step and the convergence tests --
 * stays an MPFR number: a function whose magnitude leaves the double range
 * (L_1000 is ~1e867 near its largest zero) is handled as any other, where a
 * double intermediate overflows to inf/inf.  Same stopping rule as the
 * machine secant: the residual and the step both under 10^-AccuracyGoal, or
 * the step under |x| 10^-PrecisionGoal. */
static Expr* fr_run_secant_mpfr_real(Expr* f, FrVarBind* bind,
                                     Expr* x0_expr, Expr* x1_expr,
                                     const FrOpts* opts) {
    long bits = opts->wp_bits;
    mpfr_t xa, xb, fa, fb, denom, step, tmp, scale, im, tol_acc, tol_prec;
    mpfr_inits2(bits, xa, xb, fa, fb, denom, step, tmp, scale, im,
                tol_acc, tol_prec, (mpfr_ptr)0);
    fr_mpfr_goal_tols(opts, tol_acc, tol_prec);
    Expr* result = NULL;
    bool dummy;
    double resid = INFINITY;    /* last |f(xb)|, the goal-shortfall error estimate */

    if (!get_approx_mpfr(x0_expr, xa, im, &dummy) || !mpfr_zero_p(im)
        || !get_approx_mpfr(x1_expr, xb, im, &dummy) || !mpfr_zero_p(im)) {
        fr_warn("nlnum", "Secant starting points must be real numbers");
        goto cleanup;
    }
    if (!fr_eval_mpfr_real(f, bind, xa, opts, bits, fa)) {
        fr_warn("nlnum", "non-real f at start point");
        goto cleanup;
    }
    if (mpfr_zero_p(fa)) { result = expr_new_mpfr_copy(xa); goto cleanup; }

    for (int64_t k = 0; k < opts->max_iter; k++) {
        if (!fr_eval_mpfr_real(f, bind, xb, opts, bits, fb)) {
            fr_warn("nlnum", "non-real f during iteration");
            goto cleanup;
        }
        mpfr_abs(tmp, fb, MPFR_RNDN);
        resid = mpfr_get_d(tmp, MPFR_RNDN);
        if (mpfr_zero_p(fb)) { result = expr_new_mpfr_copy(xb); goto cleanup; }
        if (mpfr_sgn(fa) * mpfr_sgn(fb) < 0) {          /* bracketed: Brent */
            result = fr_brent_mpfr_core(f, bind, xa, fa, xb, fb,
                                        opts->max_iter - k, opts);
            goto cleanup;
        }
        mpfr_sub(denom, fb, fa, MPFR_RNDN);
        if (mpfr_zero_p(denom)) {
            /* f_k == f_{k-1}: either the two iterates coincide (converged to
             * the working precision) or the secant is horizontal. */
            if (mpfr_equal_p(xa, xb)) { result = expr_new_mpfr_copy(xb); goto cleanup; }
            fr_warn("noconv", "secant denominator vanished");
            goto cleanup;
        }
        mpfr_sub(step, xb, xa, MPFR_RNDN);
        mpfr_mul(step, step, fb, MPFR_RNDN);
        mpfr_div(step, step, denom, MPFR_RNDN);
        mpfr_mul_d(step, step, opts->damping, MPFR_RNDN);
        if (!mpfr_number_p(step)) {
            fr_warn("noconv", "divergence (non-finite step)");
            goto cleanup;
        }
        mpfr_swap(xa, xb);              /* xa <- x_k                       */
        mpfr_swap(fa, fb);              /* fa <- f_k                       */
        mpfr_sub(xb, xa, step, MPFR_RNDN);
        fr_fire_monitor(opts->step_monitor);

        mpfr_abs(step, step, MPFR_RNDN);
        bool acc_ok = mpfr_cmp(tmp, tol_acc) < 0 && mpfr_cmp(step, tol_acc) < 0;
        mpfr_abs(scale, xb, MPFR_RNDN);
        mpfr_mul(scale, scale, tol_prec, MPFR_RNDN);
        bool prec_ok = mpfr_cmp(step, scale) < 0;
        if (acc_ok || prec_ok) { result = expr_new_mpfr_copy(xb); goto cleanup; }
    }
    fr_warn("cvmit", "failed to converge within %lld iterations",
            (long long)opts->max_iter);
    result = expr_new_mpfr_copy(xb);
    mpfr_abs(tmp, xb, MPFR_RNDN);
    fr_warn_goal(opts, resid, mpfr_get_d(tmp, MPFR_RNDN));
cleanup:
    mpfr_clears(xa, xb, fa, fb, denom, step, tmp, scale, im,
                tol_acc, tol_prec, (mpfr_ptr)0);
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
    double resid = INFINITY;    /* last |f(z)|, the goal-shortfall error estimate */
    int polished = -1;          /* see fr_polish_done */
    mpfr_t prev_step;
    mpfr_init2(prev_step, bits);
    mpfr_set_inf(prev_step, +1);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr* zv = fr_expr_from_complex_mpfr(zr, zi);
        Expr* a0[1] = { zv };
        Expr* fv_e = fr_eval_with_bindings(f, bind, a0, 1, opts);
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
        resid = mpfr_get_d(tmp1, MPFR_RNDN);
        bool resid_ok = mpfr_cmp(tmp1, tol_acc) < 0;
        bool f_zero = mpfr_zero_p(fr_) && mpfr_zero_p(fi);
        if (resid_ok || f_zero) {                 /* both tests, as in the real case */
            mpfr_hypot(tmp2, zr, zi, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, tol_prec, MPFR_RNDN);
            mpfr_add(tmp2, tmp2, tol_acc, MPFR_RNDN);
            if (f_zero || mpfr_cmp(prev_step, tmp2) < 0) {
                result = fr_expr_from_complex_mpfr(zr, zi);
                converged = true;
                break;
            }
        }

        if (df) {
            zv = fr_expr_from_complex_mpfr(zr, zi);
            Expr* a1[1] = { zv };
            Expr* dv_e = fr_eval_with_bindings(df, bind, a1, 1, opts);
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
        } else if (!fr_fd_deriv_mpfr_complex(f, bind, zr, zi, opts, bits, dr, di)) {
            fr_warn("nlnum", "could not evaluate finite-difference derivative");
            goto cleanup;
        }
        if (!fr_mpc_newton_step(zr, zi, fr_, fi, dr, di,
                                opts->damping, step_r, step_i, bits)) {
            if (polished >= 0 || resid_ok) {
                result = fr_expr_from_complex_mpfr(zr, zi); converged = true; break;
            }
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
        bool stop = false;
        if (polished >= 0) {
            polished++;
            mpfr_mul_2si(scale, mag, 4 - bits, MPFR_RNDN);       /* 16 ulps of |z| */
            mpfr_mul_2si(prev_step, prev_step, -1, MPFR_RNDN);
            stop = mpfr_cmp(tmp1, scale) <= 0 || mpfr_cmp(tmp1, prev_step) > 0
                   || polished >= FR_POLISH_MAX;
        } else {
            mpfr_mul(scale, mag, tol_prec, MPFR_RNDN);
            mpfr_add(scale, scale, tol_acc, MPFR_RNDN);
            if (mpfr_cmp(tmp1, scale) < 0) polished = 0;
        }
        mpfr_set(prev_step, tmp1, MPFR_RNDN);
        mpfr_clear(mag); mpfr_clear(scale);
        if (opts->step_monitor) {
            Expr* zv_new = fr_expr_from_complex_mpfr(zr, zi);
            fr_fire_step_monitor(opts, bind, zv_new);
            expr_free(zv_new);
        }
        if (stop) {
            result = fr_expr_from_complex_mpfr(zr, zi);
            converged = true;
            break;
        }
    }
    if (!converged && polished >= 0) {
        result = fr_expr_from_complex_mpfr(zr, zi);
        converged = true;
    }
    if (!converged) {
        fr_warn("cvmit", "failed to converge within %lld iterations",
                (long long)opts->max_iter);
        result = fr_expr_from_complex_mpfr(zr, zi);
        mpfr_hypot(tmp1, zr, zi, MPFR_RNDN);
        fr_warn_goal(opts, resid, mpfr_get_d(tmp1, MPFR_RNDN));
    }
cleanup:
    mpfr_clears(zr, zi, fr_, fi, dr, di, step_r, step_i,
                tol_acc, tol_prec, tmp1, tmp2, prev_step, (mpfr_ptr)0);
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
    if (eqlist->data.function.head->data.symbol.name != SYM_List) return NULL;
    size_t n = eqlist->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        items[i] = fr_normalize_eqn(eqlist->data.function.args[i]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, n);
    free(items);
    return out;
}

/* True if `e` is a List of exactly `n` elements. */
static bool fr_is_list_n(const Expr* e, size_t n) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List
        && e->data.function.arg_count == n;
}

/* A user Jacobian for n equations in n unknowns is an n x n matrix (a List of
 * n Lists of n entries).  The entries are held expressions in the variables,
 * evaluated per point, so they need only be numeric once the variables are
 * bound (g[x] with a ?NumberQ guard).  Reports FindRoot::jacsh otherwise. */
static bool fr_jacobian_shape_ok(const Expr* J, size_t n) {
    bool ok = fr_is_list_n(J, n);
    for (size_t i = 0; ok && i < n; i++)
        ok = fr_is_list_n(J->data.function.args[i], n);
    if (!ok)
        fr_warn("jacsh", "the Jacobian must be a %zu x %zu matrix for %zu "
                "equation%s", n, n, n, n == 1 ? "" : "s");
    return ok;
}

/* The derivative expression for one equation: the 1 x 1 matrix form {{d}} is
 * the portable spelling (it is what a system of one equation takes), the bare
 * expression d the scalar one.  Returns a fresh copy, or NULL after a jacsh
 * message for any other List shape. */
static Expr* fr_scalar_jacobian(const Expr* J) {
    if (fr_is_list_n(J, 1) && fr_is_list_n(J->data.function.args[0], 1))
        return expr_copy(J->data.function.args[0]->data.function.args[0]);
    if (J->type == EXPR_FUNCTION && J->data.function.head->type == EXPR_SYMBOL
        && J->data.function.head->data.symbol.name == SYM_List) {
        fr_jacobian_shape_ok(J, 1);
        return NULL;
    }
    return expr_copy((Expr*)J);
}

static Expr* fr_run_newton_system_real(Expr* flist_normalized,
                                       Expr** vars, size_t n,
                                       FrVarBind* binds,
                                       double* x_vec,
                                       const FrOpts* opts) {
    double tol_acc  = pow(10.0, -opts->acc_goal_digits);
    double tol_prec = pow(10.0, -opts->prec_goal_digits);
    if (opts->jacobian && !fr_jacobian_shape_ok(opts->jacobian, n)) return NULL;

    /* Pre-compute Jacobian symbolically: matrix[i][j] = D[f_i, var_j].
     * jz[i*n+j] flags the STRUCTURALLY-ZERO entries: for a sparse system (a
     * tridiagonal Broyden field has 3 non-zeros per row, not n) that is the
     * overwhelming majority, and skipping their compile + per-iteration eval is
     * the difference between O(n^2) and O(nnz) work. */
    Expr*** jac = malloc(sizeof(Expr**) * n);
    bool* jz = malloc(sizeof(bool) * n * n);
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
            Expr* je = jac[i][j];
            jz[i * n + j] = je &&
                ((je->type == EXPR_INTEGER && je->data.integer == 0) ||
                 (je->type == EXPR_REAL && je->data.real == 0.0));
        }
    }

    /* Compile each component f_i and Jacobian entry J_ij as a function of all n
     * variables (machine fast path); a NULL entry falls back to the interpreter.
     * Structurally-zero entries are never compiled. Systems are
     * machine-precision, so the spec is fixed here. */
    AutoCompiled**  ac_f   = malloc(sizeof(*ac_f) * n);
    AutoCompiled*** ac_jac = malloc(sizeof(*ac_jac) * n);
    for (size_t i = 0; i < n; i++) {
        ac_f[i] = autocompile_new(flist_normalized->data.function.args[i],
                                  (const Expr* const*)vars, n);
        ac_jac[i] = malloc(sizeof(**ac_jac) * n);
        for (size_t j = 0; j < n; j++)
            ac_jac[i][j] = jz[i * n + j]
                ? NULL
                : autocompile_new(jac[i][j], (const Expr* const*)vars, n);
    }

    Expr* result = NULL;
    bool converged = false;
    double resid = INFINITY;    /* last residual norm max_i |f_i|, goal estimate */
    double xmag  = 0.0;         /* last max_i |x_i|, for the relative tolerance   */

    for (int64_t k = 0; k < opts->max_iter; k++) {
        Expr** xv = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) xv[i] = expr_new_real(x_vec[i]);

        /* fvec[i] = numeric f_i at current point */
        double* fvec = malloc(sizeof(double) * n);
        double max_f_pre = 0.0;
        for (size_t i = 0; i < n; i++) {
            Expr* fi = flist_normalized->data.function.args[i];
            if (!fr_eval_compiled_double(ac_f[i], x_vec, n, &fvec[i])) {
                Expr* fv_e = fr_eval_with_bindings(fi, binds, xv, n, opts);
                if (!fv_e || !fr_expr_to_double_real(fv_e, &fvec[i])) {
                    expr_free(fv_e);
                    fr_warn("nlnum", "non-real f_%zu during iteration", i);
                    free(fvec);
                    for (size_t q = 0; q < n; q++) expr_free(xv[q]);
                    free(xv);
                    goto cleanup;
                }
                expr_free(fv_e);
            }
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
                if (jz[i*n + j]) { J[i*n + j] = 0.0; continue; }
                if (!fr_eval_compiled_double(ac_jac[i][j], x_vec, n, &J[i*n + j])) {
                    Expr* dv_e = fr_eval_with_bindings(jac[i][j], binds, xv, n, opts);
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
        /* Compute residual norm at new x.  Compiled components evaluate
         * straight to a double; xv2 is built only if some component needs
         * the interpreter (never, for a fully-compiled system). */
        Expr** xv2 = NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* fi = flist_normalized->data.function.args[i];
            double v;
            if (fr_eval_compiled_double(ac_f[i], x_vec, n, &v)) {
                if (fabs(v) > max_f) max_f = fabs(v);
                continue;
            }
            if (!xv2) {
                xv2 = malloc(sizeof(Expr*) * n);
                for (size_t q = 0; q < n; q++) xv2[q] = expr_new_real(x_vec[q]);
            }
            Expr* fv_e = fr_eval_with_bindings(fi, binds, xv2, n, opts);
            if (fv_e && fr_expr_to_double_real(fv_e, &v) && fabs(v) > max_f) {
                max_f = fabs(v);
            }
            expr_free(fv_e);
        }
        if (xv2) {
            for (size_t q = 0; q < n; q++) expr_free(xv2[q]);
            free(xv2);
        }
        free(dx);

        fr_fire_monitor(opts->step_monitor);
        resid = max_f;
        xmag  = max_x;
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
        fr_warn_goal(opts, resid, xmag);
    }
    /* Build the result list of Rules (own the new value Exprs). */
    Expr** vals = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) vals[i] = expr_new_real(x_vec[i]);
    result = fr_build_rule_list(vars, vals, n);
    free(vals);

cleanup:
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            expr_free(jac[i][j]);
            autocompiled_free(ac_jac[i][j]);
        }
        free(jac[i]);
        free(ac_jac[i]);
        autocompiled_free(ac_f[i]);
    }
    free(jac);
    free(ac_jac);
    free(ac_f);
    free(jz);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Entry point                                                        *
 * ------------------------------------------------------------------ */

/* One scalar root-finding problem, as the entry point decoded it.  All
 * pointers borrowed. */
typedef struct {
    FrMethod   method;
    FrSpecKind kind;
    Expr*      f;                /* normalised function                  */
    Expr *x0, *x1, *xmin, *xmax; /* evaluated spec values (some NULL)    */
    Expr*      deriv;            /* Newton derivative, or NULL for FD    */
    bool       want_complex;
} FrScalarProblem;

/* Run the chosen scalar method at opts->prec_mode.  Returns the root value
 * (caller owns) or NULL after a message -- or, for a machine solve that set
 * *opts->overflow, NULL silently so the caller can rerun in MPFR. */
static Expr* fr_solve_scalar(const FrScalarProblem* P, FrVarBind* bind,
                             const FrOpts* opts) {
    if (P->method == FR_METHOD_BRENT) {
        /* Brent accepts a {var, x0, x1} pair as the bracket, or a
         * full {var, xstart, xmin, xmax} spec. */
        Expr *ea, *eb;
        if (P->kind == FR_SPEC_BRACKET)        { ea = P->xmin; eb = P->xmax; }
        else if (P->kind == FR_SPEC_TWO_START) { ea = P->x0;   eb = P->x1;   }
        else {
            fr_warn("brnoth", "Brent requires {var, a, b} or {var, xstart, xmin, xmax}");
            return NULL;
        }
#ifdef USE_MPFR
        if (opts->prec_mode == FR_PREC_MPFR) {
            return fr_run_brent_mpfr_real(P->f, bind, ea, eb, opts);
        }
#endif
        double a, b;
        if (!fr_expr_to_double_real(ea, &a) || !fr_expr_to_double_real(eb, &b)) {
            fr_warn("nlnum", "Brent endpoints must be real numbers");
            return NULL;
        }
        return fr_run_brent_real(P->f, bind, a, b, opts);
    } else if (P->method == FR_METHOD_SECANT) {
        if (P->kind != FR_SPEC_TWO_START) {
            fr_warn("argt", "Secant requires a {var, x0, x1} spec");
            return NULL;
        }
#ifdef USE_MPFR
        if (opts->prec_mode == FR_PREC_MPFR) {
            return fr_run_secant_mpfr_real(P->f, bind, P->x0, P->x1, opts);
        }
#endif
        double x0, x1;
        if (!fr_expr_to_double_real(P->x0, &x0)
            || !fr_expr_to_double_real(P->x1, &x1)) {
            fr_warn("nlnum", "Secant starting points must be real numbers");
            return NULL;
        }
        return fr_run_secant_real(P->f, bind, x0, x1, opts);
    } else {
#ifdef USE_MPFR
        if (opts->prec_mode == FR_PREC_MPFR) {
            if (P->want_complex)
                return fr_run_newton_mpfr_complex(P->f, P->deriv, bind, P->x0, opts);
            return fr_run_newton_mpfr_real(P->f, P->deriv, bind, P->x0, opts);
        }
#endif
        if (P->want_complex) {
            double complex z0;
            if (!fr_expr_to_complex(P->x0, &z0)) {
                fr_warn("nlnum", "start point not numeric");
                return NULL;
            }
            return fr_run_newton_complex(P->f, P->deriv, bind, z0, opts);
        }
        double x0d;
        if (!fr_expr_to_double_real(P->x0, &x0d)) {
            fr_warn("nlnum", "start point not numeric");
            return NULL;
        }
        return fr_run_newton_real(P->f, P->deriv, bind, x0d, opts);
    }
}

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
    /* Uniform numerical-function contract: AccuracyGoal defaults to
     * MachinePrecision, PrecisionGoal to Automatic (both == WorkingPrecision/2,
     * filled in below). Both must be negative sentinels so the resolution at
     * the bottom of this block scales them with WorkingPrecision -- a positive
     * literal here pinned the accuracy goal at ~16 digits regardless of
     * WorkingPrecision, so FindRoot[..., WorkingPrecision -> 100] stopped its
     * MPFR Newton at |f| < 10^-16 and returned ~19 correct digits while
     * labelling the result Precision 100. */
    opts.acc_goal_digits = NC_GOAL_MACHINE;
    opts.prec_goal_digits = NC_GOAL_AUTO;
    opts.damping = 1.0;
    opts.jacobian = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;
    opts.main_f = NULL;
    opts.ac_f = NULL; opts.ac_f_prec = NULL;
    opts.overflow = NULL;
    opts.escalated = false;

    for (size_t i = pos_end; i < argc; i++) {
        if (!fr_apply_option(res->data.function.args[i], &opts)) return NULL;
    }

    /* Resolve any Automatic (NC_GOAL_AUTO) goal to half the working precision in
     * decimal digits; MachinePrecision / explicit / Infinity goals pass through. */
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
        && var_spec->data.function.head->data.symbol.name == SYM_List
        && var_spec->data.function.arg_count > 0) {
        Expr* first = var_spec->data.function.args[0];
        if (first->type == EXPR_FUNCTION
            && first->data.function.head->type == EXPR_SYMBOL
            && first->data.function.head->data.symbol.name == SYM_List) {
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
            fr_bind_snapshot(&binds[i], vars[i]->data.symbol.name);
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
    fr_bind_snapshot(&bind, var->data.symbol.name);

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

    /* Auto-compile the function for the scalar, machine-precision, real path.
     * The compiled program has real inputs, so complex/MPFR searches never use
     * it (fr_eval_with_bindings guards on spec.mode and a real argument). */
    if (!want_complex) {
        const Expr* vs[1] = { var };
        if (opts.prec_mode == FR_PREC_MACHINE) {
            opts.ac_f = autocompile_new(fnorm, vs, 1);
            opts.main_f = fnorm;
        }
#ifdef USE_MPFR
        else if (opts.prec_mode == FR_PREC_MPFR && opts.wp_bits > 0) {
            opts.ac_f_prec = autocompile_new_prec(fnorm, vs, 1, (long)opts.wp_bits);
            opts.main_f = fnorm;
        }
#endif
    }

    /* Newton's derivative, computed once (an escalated rerun reuses it).  When
     * the symbolic derivative is unavailable -- D returns an inert
     * Derivative[..] (e.g. Zeta', PolyGamma order derivatives) -- Newton falls
     * back to a finite-difference derivative (deriv == NULL signals that). */
    Expr* deriv = NULL;
    if (method == FR_METHOD_NEWTON) {
        if (opts.jacobian) {
            deriv = fr_scalar_jacobian(opts.jacobian);
            if (!deriv) goto scalar_cleanup;
        } else {
            deriv = fr_compute_derivative(fnorm, var);
            if (deriv && fr_has_inert_derivative(deriv)) {
                expr_free(deriv);
                deriv = NULL;
            }
        }
    }

    FrScalarProblem prob = { method, kind, fnorm, x0_e, x1_e, xmin_e, xmax_e,
                             deriv, want_complex };
    bool overflow = false;
    if (!want_complex && opts.prec_mode == FR_PREC_MACHINE) opts.overflow = &overflow;
    result = fr_solve_scalar(&prob, &bind, &opts);
#ifdef USE_MPFR
    if (!result && overflow) {
        /* A machine solve met a value outside the double range (|L_1000| ~
         * 1e867 near its largest zero): rerun the same method in MPFR, whose
         * exponent range is unbounded, at a little over twice machine
         * precision, and hand back a machine number.  The goals stay the
         * machine ones resolved above. */
        opts.overflow = NULL;
        opts.prec_mode = FR_PREC_MPFR;
        opts.wp_bits = FR_ESCALATION_BITS;
        opts.escalated = true;
        Expr* r = fr_solve_scalar(&prob, &bind, &opts);
        if (r && r->type == EXPR_MPFR && mpfr_number_p(r->data.mpfr)
            && isfinite(mpfr_get_d(r->data.mpfr, MPFR_RNDN))) {
            result = expr_new_real(mpfr_get_d(r->data.mpfr, MPFR_RNDN));
            expr_free(r);
        } else {
            result = r;          /* a root beyond the double range stays MPFR */
        }
    }
#endif
    expr_free(deriv);

scalar_cleanup:
    fr_bind_restore(&bind);
    autocompiled_free(opts.ac_f);
    autocompiled_free(opts.ac_f_prec);
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
