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
 *   FindMinimum[f,           {x, y, ...}]              n-D, auto start = 1
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
#include "compile/compile.h"
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
    FM_SPEC_VAR_ONLY,   /* {x}                       -> x0 = 1 (MMA-default) */
    FM_SPEC_SINGLE,     /* {x, x0}                                           */
    FM_SPEC_TWO_START,  /* {x, x0, x1}               -> derivative free      */
    FM_SPEC_BRACKET     /* {x, xstart, xmin, xmax}                           */
} FmSpecKind;

typedef struct {
    FmMethod   method;
    FmPrecMode prec_mode;
    long       wp_bits;          /* MPFR bits when prec_mode == MPFR    */
    int64_t    max_iter;         /* default 500                          */
    bool       max_iter_set;     /* true if MaxIterations given explicitly */
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
 * For the inner solver we evaluate each constraint expression directly,
 * and (when present) its symbolic gradient — needed so the augmented
 * objective is differentiated consistently with its value. */
typedef struct {
    Expr*  expr;       /* owned: feasible ≡ (expr <= 0) or (expr == 0)        */
    Expr** grad_exprs; /* owned: ∇expr w.r.t. vars (length n), or NULL → FD   */
    bool   equality;   /* true → equality constraint                          */
} FmGenCon;

/* A disjunctive (Or) constraint: feasible ≡ at least ONE branch is feasible.
 * Stored as its boolean-of-comparisons subtree over the effective variables.
 * The penalty contribution is the MINIMUM branch penalty (fm_bool_penalty), so
 * a point satisfying any one branch scores zero and the (total ≤ NM_FEAS_EPS)
 * feasibility test still means "feasible". Because the min() is non-smooth it is
 * consumed only by the derivative-free global search (NMinimize); the smooth
 * local polish leaves disjunction feasibility to the Deb accept/reject gate, and
 * FindMinimum's gradient penalty method rejects Or outright. Evaluated by the
 * interpreter rather than a compiled program — disjunctive problems are rare. */
typedef struct {
    Expr*  expr;       /* owned: the Or[...] constraint subtree                */
} FmDisjunction;

/* ------------------------------------------------------------------ *
 *  Diagnostic helper                                                  *
 * ------------------------------------------------------------------ */

/* When true, fm_warn is a no-op. NMinimize sets it around its internal
 * local-solver calls so the penalty/line-search chatter that is expected
 * during global search does not reach the user (NMinimize reports feasibility
 * itself via its {Infinity, ...} result). */
static bool g_fm_quiet = false;

static void fm_warn(const char* fn, const char* tag, const char* fmt, ...) {
    if (g_fm_quiet) return;
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
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fm_is_known_option_name(lhs->data.symbol.name);
}

static bool fm_parse_working_precision(Expr* val,
                                       FmPrecMode* mode, long* bits) {
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
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
        if (val->data.symbol.name == SYM_Automatic) { *digits_out = -1.0; return true; }
        if (val->data.symbol.name == SYM_Infinity)  { *digits_out = INFINITY; return true; }
        return false;
    }
    return fm_expr_to_double_real(val, digits_out);
}

static bool fm_apply_option(Expr* rule, FmOpts* opts) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;

    if (name == SYM_Method) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
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
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
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
    /* Mute expected numeric-domain diagnostics (Power::infy from a 1/0 in a
     * gradient term on a non-differentiable ridge, Infinity::indet, ...): the
     * search evaluates the user function at many trial points and treats any
     * non-finite result as a bad point, so these messages are pure noise.
     * Matches Mathematica, which quiets NMinimize's internal evaluation. */
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(f));
    Expr* num = raw ? numericalize(raw, spec) : NULL;
    arith_warnings_mute_pop();
    expr_free(raw);
    return num;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  MPFR scalar-evaluation core                                        *
 * ------------------------------------------------------------------ *
 * The MPFR optimizer paths reuse `fm_eval_with_bindings` — which already
 * accepts any Expr* as a substitution — by handing it `expr_new_mpfr_copy(x)`
 * instead of `expr_new_real(x)`. The numericalize call inside that helper
 * then keeps the entire arithmetic chain at the requested precision. */

/* Extract a real-valued result from an evaluated MPFR expression into
 * `out`. Tolerates a tiny imaginary residue (within ~4 digits of the
 * working precision) to mirror the double-path's behaviour around
 * subtraction cancellation in Complex[] evaluations. */
static bool fm_mpfr_extract_real(Expr* res, long bits, mpfr_t out) {
    if (!res) return false;
    mpfr_t im;
    mpfr_init2(im, bits);
    bool inexact = false;
    bool ok = get_approx_mpfr(res, out, im, &inexact);
    if (ok && !mpfr_zero_p(im)) {
        mpfr_t abs_im, abs_re, thresh;
        mpfr_init2(abs_im, bits);
        mpfr_init2(abs_re, bits);
        mpfr_init2(thresh, bits);
        mpfr_abs(abs_im, im, MPFR_RNDN);
        mpfr_abs(abs_re, out, MPFR_RNDN);
        mpfr_add_ui(abs_re, abs_re, 1, MPFR_RNDN);
        long sub = (long)numeric_bits_to_digits(bits) - 4;
        if (sub < 1) sub = 1;
        mpfr_set_ui(thresh, 10, MPFR_RNDN);
        mpfr_pow_si(thresh, thresh, -sub, MPFR_RNDN);
        mpfr_mul(thresh, thresh, abs_re, MPFR_RNDN);
        if (mpfr_cmp(abs_im, thresh) > 0) ok = false;
        mpfr_clears(abs_im, abs_re, thresh, (mpfr_ptr)0);
    }
    mpfr_clear(im);
    return ok;
}

/* Evaluate the bound objective at MPFR precision using a caller-built
 * array of MPFR-leaf Expr substitutions. Routing through `Expr*`
 * avoids the `mpfr_t*`/parameter-decay hazards: in C, `mpfr_t` is a
 * 1-element array, so taking `&local` of a parameter typed `mpfr_t a`
 * yields a pointer-to-pointer, not an "array of mpfr_t" — passing such
 * a pointer to a function expecting `mpfr_t const*` segfaults the
 * moment the callee dereferences it. */
static bool fm_eval_scalar_mpfr_exprs(Expr* f, FmVarBind* binds,
                                      Expr* const* xv, size_t n,
                                      const FmOpts* opts, mpfr_t out) {
    long bits = opts->wp_bits;
    Expr* res = fm_eval_with_bindings(f, binds, xv, n,
                                      opts->eval_monitor,
                                      fm_numeric_spec(opts));
    if (!res) return false;
    bool ok = fm_mpfr_extract_real(res, bits, out);
    expr_free(res);
    return ok;
}

/* 1D convenience: build one MPFR leaf from `x`, evaluate, extract. */
static bool fm_eval_scalar_mpfr_1d(Expr* f, FmVarBind* binds,
                                   const mpfr_t x, const FmOpts* opts,
                                   mpfr_t out) {
    Expr* xv = expr_new_mpfr_copy(x);
    Expr* arr[1] = { xv };
    bool ok = fm_eval_scalar_mpfr_exprs(f, binds, arr, 1, opts, out);
    expr_free(xv);
    return ok;
}

/* n-D convenience: build a fresh MPFR-leaf array from the iterate, then
 * evaluate/extract. Accepts the iterate as `mpfr_t* x_vec` (an
 * `__mpfr_struct (*)[1]`) and indexes it with `x_vec[i]`, which decays
 * to `__mpfr_struct *` (= `mpfr_t` argument) cleanly. */
static bool fm_eval_scalar_mpfr(Expr* f, FmVarBind* binds,
                                mpfr_t* x_vec, size_t n,
                                const FmOpts* opts, mpfr_t out) {
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_mpfr_copy(x_vec[i]);
    bool ok = fm_eval_scalar_mpfr_exprs(f, binds, xv, n, opts, out);
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    return ok;
}

/* Set `out = 10^-digits` at the current MPFR precision. Treats Infinity
 * (digits == +inf) and large finite values uniformly: anything past the
 * representable exponent becomes +0. */
static void fm_tol_from_digits(mpfr_t out, double digits) {
    if (isinf(digits) || digits > 1e9) { mpfr_set_zero(out, +1); return; }
    mpfr_set_ui(out, 10, MPFR_RNDN);
    mpfr_pow_si(out, out, -(long)digits, MPFR_RNDN);
}
#endif /* USE_MPFR */

/* Machine-precision fast path for the local solvers. When a solve registers its
 * objective here (nm_minimize_driver), any fm_eval_scalar on that *exact*
 * objective expression — pointer identity, matching arity — is served by the
 * compiled program the global search already uses (NmDriver.f_prog) instead of
 * the interpreter. This is what makes RandomSearch's per-start local polish run
 * compiled rather than re-binding 20 OwnValues + deep-copying + evaluating +
 * numericalizing per point. Constraint and gradient sub-expressions are distinct
 * Expr* (and may have a different arity), so they correctly fall through to the
 * interpreter; a non-finite compiled result also falls through, so the point is
 * scored exactly as before. Registered/deregistered around each solve in
 * nm_minimize_driver (a plain reset — see the note there on why save/restore is
 * unnecessary). */
static Expr*            g_fm_obj_expr  = NULL;
static CompiledProgram* g_fm_obj_prog  = NULL;
static size_t           g_fm_obj_nargs = 0;

/* Companion registry for the exact symbolic gradient (FindMinimum). When a solve
 * registers its gradient-component array here, `fm_eval_gradient` evaluates each
 * component through its compiled program instead of the interpreter — the SAME
 * symbolic ∂f/∂x_i, just lowered, so the gradient stays exact (no finite
 * differences) while running on the register machine. Keyed by the g_exprs array
 * pointer + arity, so constraint gradients (a different array) stay on the
 * interpreter; a NULL or non-finite component falls back per-component. */
static Expr**            g_fm_grad_exprs = NULL;
static CompiledProgram** g_fm_grad_progs = NULL;   /* len g_fm_grad_n, entries may be NULL */
static size_t            g_fm_grad_n     = 0;

/* Evaluate the bound objective and return a double; NULL on failure. */
static bool fm_eval_scalar(Expr* f, FmVarBind* binds,
                           const double* x, size_t n,
                           const FmOpts* opts, double* out) {
    if (g_fm_obj_prog && f == g_fm_obj_expr && n == g_fm_obj_nargs
        && compiled_eval_real(g_fm_obj_prog, x, out) && isfinite(*out))
        return true;
    Expr** xv = (Expr**)calloc(n ? n : 1, sizeof(Expr*));
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
    Expr* call = expr_new_function(expr_new_symbol(SYM_D), args, 2);
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
    /* Compiled fast path (see g_fm_grad_*): each component whose program is
     * registered is evaluated on the register machine. The interpreter value
     * bindings `xv` are built lazily — only if some component has no program or
     * returns a non-finite result — so an all-compiled gradient allocates
     * nothing. */
    bool reg = (g_exprs == g_fm_grad_exprs && n == g_fm_grad_n && g_fm_grad_progs);
    Expr** xv = NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        if (reg && g_fm_grad_progs[i]
            && compiled_eval_real(g_fm_grad_progs[i], x, &g_out[i])
            && isfinite(g_out[i]))
            continue;
        if (!xv) {
            xv = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < n; j++) xv[j] = expr_new_real(x[j]);
        }
        Expr* gi = fm_eval_with_bindings(g_exprs[i], binds, xv, n,
                                         opts->eval_monitor,
                                         fm_numeric_spec(opts));
        if (!gi || !fm_expr_to_double_real(gi, &g_out[i])) { ok = false; expr_free(gi); break; }
        expr_free(gi);
    }
    if (xv) { for (size_t i = 0; i < n; i++) expr_free(xv[i]); free(xv); }
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
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_real(fmin), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}

#ifdef USE_MPFR
/* MPFR analogue: the result components are stored as EXPR_MPFR leaves at
 * the working precision rather than EXPR_REAL. */
static Expr* fm_build_result_mpfr(const mpfr_t fmin, Expr** vars,
                                  mpfr_t const* vals, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r_args[2] = { expr_copy(vars[i]), expr_new_mpfr_copy(vals[i]) };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_mpfr_copy(fmin), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}

/* Allocate `count` MPFR scalars at precision `bits`. */
static mpfr_t* fm_mpfr_array(size_t count, long bits) {
    mpfr_t* arr = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(arr[i], bits);
    return arr;
}

static void fm_mpfr_array_free(mpfr_t* arr, size_t count) {
    if (!arr) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(arr[i]);
    free(arr);
}

/* MPFR central-difference gradient. h_rel ~ 10^-(digits/2). */
static bool fm_grad_finite_diff_mpfr(Expr* f, FmVarBind* binds,
                                     mpfr_t const* x, size_t n,
                                     const FmOpts* opts, mpfr_t* g_out) {
    long bits = opts->wp_bits;
    mpfr_t* xp = fm_mpfr_array(n, bits);
    for (size_t i = 0; i < n; i++) mpfr_set(xp[i], x[i], MPFR_RNDN);
    mpfr_t step, scale, f0, f1, denom, h_rel;
    mpfr_init2(step, bits); mpfr_init2(scale, bits);
    mpfr_init2(f0, bits); mpfr_init2(f1, bits);
    mpfr_init2(denom, bits); mpfr_init2(h_rel, bits);
    /* h_rel = 10^-(digits/2), capped at 1e-7 to keep the difference
     * informative even at moderate precision. */
    double hd = numeric_bits_to_digits(bits) / 2.0;
    if (hd < 7.0) hd = 7.0;
    mpfr_set_ui(h_rel, 10, MPFR_RNDN);
    mpfr_pow_si(h_rel, h_rel, -(long)hd, MPFR_RNDN);
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        mpfr_abs(scale, x[i], MPFR_RNDN);
        mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
        if (mpfr_cmp(scale, one) < 0) mpfr_set(scale, one, MPFR_RNDN);
        mpfr_clear(one);
        mpfr_mul(step, scale, h_rel, MPFR_RNDN);
        mpfr_add(xp[i], x[i], step, MPFR_RNDN);
        if (!fm_eval_scalar_mpfr(f, binds, xp, n, opts, f1)) { ok = false; break; }
        mpfr_sub(xp[i], x[i], step, MPFR_RNDN);
        if (!fm_eval_scalar_mpfr(f, binds, xp, n, opts, f0)) { ok = false; break; }
        mpfr_sub(g_out[i], f1, f0, MPFR_RNDN);
        mpfr_mul_ui(denom, step, 2, MPFR_RNDN);
        mpfr_div(g_out[i], g_out[i], denom, MPFR_RNDN);
        mpfr_set(xp[i], x[i], MPFR_RNDN);
    }
    fm_mpfr_array_free(xp, n);
    mpfr_clears(step, scale, f0, f1, denom, h_rel, (mpfr_ptr)0);
    return ok;
}

/* MPFR symbolic gradient evaluator. */
static bool fm_eval_gradient_mpfr(Expr** g_exprs, FmVarBind* binds,
                                  mpfr_t const* x, size_t n,
                                  const FmOpts* opts, mpfr_t* g_out) {
    long bits = opts->wp_bits;
    Expr** xv = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) xv[i] = expr_new_mpfr_copy(x[i]);
    bool ok = true;
    mpfr_t im;
    mpfr_init2(im, bits);
    for (size_t i = 0; i < n; i++) {
        Expr* gi = fm_eval_with_bindings(g_exprs[i], binds, xv, n,
                                         opts->eval_monitor,
                                         fm_numeric_spec(opts));
        if (!gi) { ok = false; break; }
        bool inexact = false;
        bool got = get_approx_mpfr(gi, g_out[i], im, &inexact);
        if (!got || !mpfr_zero_p(im)) {
            /* Tolerate tiny imaginary residue (matches scalar path). */
            if (got) {
                mpfr_t abs_im, abs_re, thresh;
                mpfr_init2(abs_im, bits);
                mpfr_init2(abs_re, bits);
                mpfr_init2(thresh, bits);
                mpfr_abs(abs_im, im, MPFR_RNDN);
                mpfr_abs(abs_re, g_out[i], MPFR_RNDN);
                mpfr_add_ui(abs_re, abs_re, 1, MPFR_RNDN);
                long sub = (long)numeric_bits_to_digits(bits) - 4;
                if (sub < 1) sub = 1;
                mpfr_set_ui(thresh, 10, MPFR_RNDN);
                mpfr_pow_si(thresh, thresh, -sub, MPFR_RNDN);
                mpfr_mul(thresh, thresh, abs_re, MPFR_RNDN);
                if (mpfr_cmp(abs_im, thresh) > 0) got = false;
                mpfr_clears(abs_im, abs_re, thresh, (mpfr_ptr)0);
            }
            if (!got) { ok = false; expr_free(gi); break; }
        }
        expr_free(gi);
    }
    mpfr_clear(im);
    for (size_t i = 0; i < n; i++) expr_free(xv[i]);
    free(xv);
    return ok;
}

/* MPFR project x into box. */
static void fm_project_box_mpfr(mpfr_t* x, size_t n, const FmBox* boxes) {
    long bits = mpfr_get_prec(x[0]);
    mpfr_t bnd; mpfr_init2(bnd, bits);
    for (size_t i = 0; i < n; i++) {
        if (boxes[i].has_lo) {
            mpfr_set_d(bnd, boxes[i].lo, MPFR_RNDN);
            if (mpfr_cmp(x[i], bnd) < 0) mpfr_set(x[i], bnd, MPFR_RNDN);
        }
        if (boxes[i].has_hi) {
            mpfr_set_d(bnd, boxes[i].hi, MPFR_RNDN);
            if (mpfr_cmp(x[i], bnd) > 0) mpfr_set(x[i], bnd, MPFR_RNDN);
        }
    }
    mpfr_clear(bnd);
}

/* MPFR Armijo line search. mu == 0 → plain f; otherwise augmented (not
 * supported yet — penalty/MPFR is deferred). */
static bool fm_line_search_mpfr(Expr* f, FmVarBind* binds, size_t n,
                                mpfr_t const* x, mpfr_t const* d,
                                const mpfr_t f0, const mpfr_t g_dot_d,
                                const FmBox* boxes,
                                const FmOpts* opts,
                                mpfr_t alpha_out, mpfr_t f_out, mpfr_t* x_out) {
    long bits = opts->wp_bits;
    mpfr_t dnorm, alpha, fnew, accept, c1, alpha_d, candidate;
    mpfr_init2(dnorm, bits); mpfr_init2(alpha, bits);
    mpfr_init2(fnew, bits); mpfr_init2(accept, bits);
    mpfr_init2(c1, bits); mpfr_init2(alpha_d, bits);
    mpfr_init2(candidate, bits);
    mpfr_set_d(c1, 1e-4, MPFR_RNDN);
    mpfr_set_zero(dnorm, +1);
    for (size_t i = 0; i < n; i++) {
        mpfr_mul(candidate, d[i], d[i], MPFR_RNDN);
        mpfr_add(dnorm, dnorm, candidate, MPFR_RNDN);
    }
    mpfr_sqrt(dnorm, dnorm, MPFR_RNDN);
    mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
    if (mpfr_cmp(dnorm, one) > 0) mpfr_ui_div(alpha, 1, dnorm, MPFR_RNDN);
    else                            mpfr_set(alpha, one, MPFR_RNDN);
    mpfr_clear(one);
    bool found = false;
    for (int k = 0; k < 60; k++) {
        for (size_t i = 0; i < n; i++) {
            mpfr_mul(candidate, alpha, d[i], MPFR_RNDN);
            mpfr_add(x_out[i], x[i], candidate, MPFR_RNDN);
        }
        if (boxes) fm_project_box_mpfr(x_out, n, boxes);
        if (!fm_eval_scalar_mpfr(f, binds, x_out, n, opts, fnew)) {
            mpfr_div_ui(alpha, alpha, 2, MPFR_RNDN);
            continue;
        }
        if (boxes) {
            /* Projected-step acceptance: just need f decrease. */
            mpfr_abs(accept, f0, MPFR_RNDN);
            mpfr_t rhs; mpfr_init2(rhs, bits);
            mpfr_set_d(rhs, 1e-12, MPFR_RNDN);
            mpfr_mul(accept, accept, rhs, MPFR_RNDN);
            mpfr_sub(rhs, f0, accept, MPFR_RNDN);
            mpfr_clear(accept); mpfr_init2(accept, bits);
            mpfr_set(accept, rhs, MPFR_RNDN);
            mpfr_clear(rhs);
            if (mpfr_cmp(fnew, accept) <= 0) {
                mpfr_set(alpha_out, alpha, MPFR_RNDN);
                mpfr_set(f_out, fnew, MPFR_RNDN);
                found = true; break;
            }
        } else {
            /* Standard Armijo: f(x + α d) ≤ f0 + c1·α·(g·d). */
            mpfr_mul(alpha_d, c1, alpha, MPFR_RNDN);
            mpfr_mul(alpha_d, alpha_d, g_dot_d, MPFR_RNDN);
            mpfr_add(accept, f0, alpha_d, MPFR_RNDN);
            if (mpfr_cmp(fnew, accept) <= 0) {
                mpfr_set(alpha_out, alpha, MPFR_RNDN);
                mpfr_set(f_out, fnew, MPFR_RNDN);
                found = true; break;
            }
        }
        mpfr_div_ui(alpha, alpha, 2, MPFR_RNDN);
        /* Stop if alpha < 10^-(digits) — anything finer is below the
         * representable resolution at the working precision. */
        mpfr_t floor_alpha; mpfr_init2(floor_alpha, bits);
        long edig = (long)numeric_bits_to_digits(bits);
        mpfr_set_ui(floor_alpha, 10, MPFR_RNDN);
        mpfr_pow_si(floor_alpha, floor_alpha, -edig - 5, MPFR_RNDN);
        bool tiny = (mpfr_cmpabs(alpha, floor_alpha) < 0);
        mpfr_clear(floor_alpha);
        if (tiny) break;
    }
    mpfr_clears(dnorm, alpha, fnew, accept, c1, alpha_d, candidate, (mpfr_ptr)0);
    return found;
}

/* BFGS at MPFR precision. Constraints / penalty paths are NOT supported
 * yet at MPFR (a follow-up will lift the existing penalty machinery
 * to mpfr_t when there's user demand); callers must route through the
 * machine-precision path when general constraints are present. */
static bool fm_run_bfgs_mpfr(Expr* f, Expr** vars, size_t n,
                             FmVarBind* binds, Expr** g_exprs,
                             mpfr_t* x, /* in/out */
                             const FmBox* boxes,
                             const FmOpts* opts,
                             mpfr_t fx_out) {
    (void)vars;
    long bits = opts->wp_bits;
    mpfr_t* H = fm_mpfr_array(n * n, bits);
    mpfr_t* g = fm_mpfr_array(n, bits);
    mpfr_t* g_new = fm_mpfr_array(n, bits);
    mpfr_t* d = fm_mpfr_array(n, bits);
    mpfr_t* x_new = fm_mpfr_array(n, bits);
    mpfr_t* s_v = fm_mpfr_array(n, bits);
    mpfr_t* y_v = fm_mpfr_array(n, bits);
    mpfr_t* Hy = fm_mpfr_array(n, bits);
    mpfr_t fx, fx_new, alpha, g_dot_d, gnorm, tol_acc, tol_prec;
    mpfr_t tmp, sy, rho, yHy, coef, max_step, max_x;
    mpfr_init2(fx, bits); mpfr_init2(fx_new, bits);
    mpfr_init2(alpha, bits); mpfr_init2(g_dot_d, bits);
    mpfr_init2(gnorm, bits); mpfr_init2(tol_acc, bits); mpfr_init2(tol_prec, bits);
    mpfr_init2(tmp, bits); mpfr_init2(sy, bits); mpfr_init2(rho, bits);
    mpfr_init2(yHy, bits); mpfr_init2(coef, bits);
    mpfr_init2(max_step, bits); mpfr_init2(max_x, bits);

    fm_tol_from_digits(tol_acc, opts->acc_goal_digits);
    fm_tol_from_digits(tol_prec, opts->prec_goal_digits);

    /* H ← I. */
    for (size_t i = 0; i < n * n; i++) mpfr_set_zero(H[i], +1);
    for (size_t i = 0; i < n; i++) mpfr_set_ui(H[i * n + i], 1, MPFR_RNDN);
    if (boxes) fm_project_box_mpfr(x, n, boxes);

    bool ok = false;
    if (!fm_eval_scalar_mpfr(f, binds, x, n, opts, fx)) goto cleanup;

    bool got_grad = g_exprs
        && fm_eval_gradient_mpfr(g_exprs, binds, (mpfr_t const*)x, n, opts, g);
    if (!got_grad) got_grad = fm_grad_finite_diff_mpfr(f, binds, (mpfr_t const*)x, n, opts, g);
    if (!got_grad) {
        fm_warn(g_fm_name, "nlnum", "MPFR gradient evaluation failed at start point");
        goto cleanup;
    }

    for (int64_t k = 0; k < opts->max_iter; k++) {
        /* ‖g‖₂ < tol_acc. */
        mpfr_set_zero(gnorm, +1);
        for (size_t i = 0; i < n; i++) {
            mpfr_mul(tmp, g[i], g[i], MPFR_RNDN);
            mpfr_add(gnorm, gnorm, tmp, MPFR_RNDN);
        }
        mpfr_sqrt(gnorm, gnorm, MPFR_RNDN);
        if (mpfr_cmp(gnorm, tol_acc) < 0) { ok = true; break; }

        /* d = -H g. */
        for (size_t i = 0; i < n; i++) {
            mpfr_set_zero(d[i], +1);
            for (size_t j = 0; j < n; j++) {
                mpfr_mul(tmp, H[i * n + j], g[j], MPFR_RNDN);
                mpfr_add(d[i], d[i], tmp, MPFR_RNDN);
            }
            mpfr_neg(d[i], d[i], MPFR_RNDN);
        }
        mpfr_set_zero(g_dot_d, +1);
        for (size_t i = 0; i < n; i++) {
            mpfr_mul(tmp, g[i], d[i], MPFR_RNDN);
            mpfr_add(g_dot_d, g_dot_d, tmp, MPFR_RNDN);
        }
        if (mpfr_sgn(g_dot_d) >= 0) {
            /* Reset H to I, fall back to steepest descent. */
            for (size_t i = 0; i < n * n; i++) mpfr_set_zero(H[i], +1);
            for (size_t i = 0; i < n; i++) mpfr_set_ui(H[i * n + i], 1, MPFR_RNDN);
            for (size_t i = 0; i < n; i++) mpfr_neg(d[i], g[i], MPFR_RNDN);
            mpfr_set_zero(g_dot_d, +1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(tmp, g[i], d[i], MPFR_RNDN);
                mpfr_add(g_dot_d, g_dot_d, tmp, MPFR_RNDN);
            }
        }

        bool ls_ok = fm_line_search_mpfr(f, binds, n,
                                         (mpfr_t const*)x, (mpfr_t const*)d,
                                         fx, g_dot_d, boxes, opts,
                                         alpha, fx_new, x_new);
        if (!ls_ok) {
            fm_warn(g_fm_name, "lstol",
                    "line search (MPFR) failed at iter %lld", (long long)k);
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        /* max step + max |x|. */
        mpfr_set_zero(max_step, +1);
        mpfr_set_zero(max_x, +1);
        for (size_t i = 0; i < n; i++) {
            mpfr_sub(tmp, x_new[i], x[i], MPFR_RNDN); mpfr_abs(tmp, tmp, MPFR_RNDN);
            if (mpfr_cmp(tmp, max_step) > 0) mpfr_set(max_step, tmp, MPFR_RNDN);
            mpfr_abs(tmp, x_new[i], MPFR_RNDN);
            if (mpfr_cmp(tmp, max_x) > 0) mpfr_set(max_x, tmp, MPFR_RNDN);
        }

        bool got_ng = g_exprs
            && fm_eval_gradient_mpfr(g_exprs, binds, (mpfr_t const*)x_new, n, opts, g_new);
        if (!got_ng) got_ng = fm_grad_finite_diff_mpfr(f, binds, (mpfr_t const*)x_new, n, opts, g_new);
        if (!got_ng) {
            fm_warn(g_fm_name, "nlnum", "MPFR gradient failed in iteration");
            for (size_t i = 0; i < n; i++) mpfr_set(x[i], x_new[i], MPFR_RNDN);
            mpfr_set(fx, fx_new, MPFR_RNDN);
            break;
        }

        /* s = x_new - x; y = g_new - g; sy = s·y. */
        for (size_t i = 0; i < n; i++) {
            mpfr_sub(s_v[i], x_new[i], x[i], MPFR_RNDN);
            mpfr_sub(y_v[i], g_new[i], g[i], MPFR_RNDN);
        }
        mpfr_set_zero(sy, +1);
        for (size_t i = 0; i < n; i++) {
            mpfr_mul(tmp, s_v[i], y_v[i], MPFR_RNDN);
            mpfr_add(sy, sy, tmp, MPFR_RNDN);
        }
        mpfr_t sy_thresh; mpfr_init2(sy_thresh, bits);
        mpfr_set_d(sy_thresh, 1e-12, MPFR_RNDN);
        if (mpfr_cmp(sy, sy_thresh) > 0) {
            mpfr_ui_div(rho, 1, sy, MPFR_RNDN);
            /* Hy = H y. */
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(Hy[i], +1);
                for (size_t j = 0; j < n; j++) {
                    mpfr_mul(tmp, H[i * n + j], y_v[j], MPFR_RNDN);
                    mpfr_add(Hy[i], Hy[i], tmp, MPFR_RNDN);
                }
            }
            mpfr_set_zero(yHy, +1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(tmp, y_v[i], Hy[i], MPFR_RNDN);
                mpfr_add(yHy, yHy, tmp, MPFR_RNDN);
            }
            /* coef = (sy + yHy) * rho^2. */
            mpfr_add(coef, sy, yHy, MPFR_RNDN);
            mpfr_mul(coef, coef, rho, MPFR_RNDN);
            mpfr_mul(coef, coef, rho, MPFR_RNDN);
            /* H ← H + coef * s s^T - rho * (Hy s^T + s Hy^T). */
            mpfr_t a_, b_, c_;
            mpfr_init2(a_, bits); mpfr_init2(b_, bits); mpfr_init2(c_, bits);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    mpfr_mul(a_, coef, s_v[i], MPFR_RNDN);
                    mpfr_mul(a_, a_, s_v[j], MPFR_RNDN);
                    mpfr_mul(b_, Hy[i], s_v[j], MPFR_RNDN);
                    mpfr_mul(c_, s_v[i], Hy[j], MPFR_RNDN);
                    mpfr_add(b_, b_, c_, MPFR_RNDN);
                    mpfr_mul(b_, b_, rho, MPFR_RNDN);
                    mpfr_add(H[i * n + j], H[i * n + j], a_, MPFR_RNDN);
                    mpfr_sub(H[i * n + j], H[i * n + j], b_, MPFR_RNDN);
                }
            }
            mpfr_clears(a_, b_, c_, (mpfr_ptr)0);
        }
        mpfr_clear(sy_thresh);

        for (size_t i = 0; i < n; i++) {
            mpfr_set(x[i], x_new[i], MPFR_RNDN);
            mpfr_set(g[i], g_new[i], MPFR_RNDN);
        }
        mpfr_set(fx, fx_new, MPFR_RNDN);

        /* PrecisionGoal: |step| < tol_prec * |x|. */
        mpfr_t scale; mpfr_init2(scale, bits);
        mpfr_mul(scale, tol_prec, max_x, MPFR_RNDN);
        bool small = (mpfr_cmp(max_step, scale) < 0);
        mpfr_clear(scale);
        if (small) { ok = true; break; }
    }
    mpfr_set(fx_out, fx, MPFR_RNDN);
    ok = true;
cleanup:
    fm_mpfr_array_free(H, n * n);
    fm_mpfr_array_free(g, n);
    fm_mpfr_array_free(g_new, n);
    fm_mpfr_array_free(d, n);
    fm_mpfr_array_free(x_new, n);
    fm_mpfr_array_free(s_v, n);
    fm_mpfr_array_free(y_v, n);
    fm_mpfr_array_free(Hy, n);
    mpfr_clears(fx, fx_new, alpha, g_dot_d, gnorm, tol_acc, tol_prec,
                tmp, sy, rho, yHy, coef, max_step, max_x, (mpfr_ptr)0);
    return ok;
}
#endif /* USE_MPFR */

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
        /* Bare variable, e.g. FindMinimum[f, {x, y, ...}] entry.
         * Default x0 = 1.0 to match Mathematica and to avoid the common
         * pitfall of starting at the saddle/critical point of oscillatory
         * functions (Sin[x] Sin[2y], Cos[...]+Sin[...] etc. all have a
         * vanishing gradient at the origin, which trivially "converges"
         * the inner solver). */
        *var_out = spec;
        *x0_out = expr_new_real(1.0);
        return FM_SPEC_VAR_ONLY;
    }
    if (spec->type != EXPR_FUNCTION) return FM_SPEC_BAD;
    if (spec->data.function.head->type != EXPR_SYMBOL) return FM_SPEC_BAD;
    if (spec->data.function.head->data.symbol.name != SYM_List) return FM_SPEC_BAD;

    size_t n = spec->data.function.arg_count;
    if (n < 1 || n > 4) return FM_SPEC_BAD;

    Expr* var = spec->data.function.args[0];
    if (var->type != EXPR_SYMBOL) return FM_SPEC_BAD;
    *var_out = var;

    if (n == 1) {
        /* {x} with no initial value: same default as bare-symbol form. */
        *x0_out = expr_new_real(1.0);
        return FM_SPEC_VAR_ONLY;
    }
    Expr* x0_raw = spec->data.function.args[1];
    if (x0_raw->type == EXPR_FUNCTION
        && x0_raw->data.function.head->type == EXPR_SYMBOL
        && x0_raw->data.function.head->data.symbol.name == SYM_List) {
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
    const char* op = head->data.symbol.name;
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
            if (vars[i]->data.symbol.name == a->data.symbol.name) {
                var_idx = (int64_t)i; c_side = b; var_left = true; break;
            }
        }
    }
    if (var_idx < 0 && b->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < nvars; i++) {
            if (vars[i]->data.symbol.name == b->data.symbol.name) {
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
    const char* op = head->data.symbol.name;
    Expr* lhs = cmp->data.function.args[0];
    Expr* rhs = cmp->data.function.args[1];

    /* lhs - rhs */
    Expr* diff_args[2] = { expr_copy(lhs), expr_copy(rhs) };
    Expr* lhs_minus_rhs = expr_new_function(expr_new_symbol(SYM_Subtract), diff_args, 2);

    if (op == SYM_LessEqual || op == SYM_Less) {
        /* lhs op rhs  →  lhs - rhs <= 0 */
        *expr_out = lhs_minus_rhs;
        *equality_out = false;
        return true;
    }
    if (op == SYM_GreaterEqual || op == SYM_Greater) {
        /* lhs op rhs  →  rhs - lhs <= 0  →  -(lhs - rhs) <= 0 */
        Expr* neg_args[2] = { expr_new_integer(-1), lhs_minus_rhs };
        *expr_out = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
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

/* True if `c` is a boolean-of-comparisons tree the disjunction penalty
 * evaluator (fm_bool_penalty) understands: And / Or of {binary comparison,
 * Inequality chain}. Validated when a disjunction is collected so an
 * unsupported shape is reported as nimpl at parse time rather than silently
 * scoring every trial point infeasible during the search. */
static bool fm_bool_supported(Expr* c) {
    if (!c || c->type != EXPR_FUNCTION
        || c->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = c->data.function.head->data.symbol.name;
    size_t ac = c->data.function.arg_count;
    if (h == SYM_And || h == SYM_Or) {
        if (ac == 0) return false;
        for (size_t i = 0; i < ac; i++)
            if (!fm_bool_supported(c->data.function.args[i])) return false;
        return true;
    }
    if (h == SYM_Inequality) {
        if (ac < 3 || (ac & 1u) != 1) return false;
        for (size_t k = 0; 2 * k + 2 < ac; k++)
            if (c->data.function.args[2 * k + 1]->type != EXPR_SYMBOL) return false;
        return true;
    }
    /* single binary comparison: reuse fm_constraint_to_g's acceptance test */
    Expr* g = NULL; bool eq = false;
    if (!fm_constraint_to_g(c, &g, &eq)) return false;
    expr_free(g);
    return true;
}

/* Walk an And[...] tree, an Inequality[...] node, or a single binary
 * comparison and accumulate boxes / general constraints. A top-level (or
 * And-nested) Or[...] is collected into the disjunction sink when one is
 * provided; if `disj_inout` is NULL (FindMinimum), Or is rejected as nimpl.
 * Returns false if any branch is unsupported. */
static bool fm_collect_constraints(Expr* cons, Expr** vars, size_t nvars,
                                   FmBox* boxes,
                                   FmGenCon** gens_inout, size_t* ngens_inout,
                                   size_t* gens_cap_inout,
                                   FmDisjunction** disj_inout,
                                   size_t* ndisj_inout, size_t* disj_cap_inout) {
    if (!cons) return true;
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol.name == SYM_And) {
        for (size_t i = 0; i < cons->data.function.arg_count; i++) {
            if (!fm_collect_constraints(cons->data.function.args[i], vars, nvars,
                                        boxes, gens_inout, ngens_inout,
                                        gens_cap_inout, disj_inout,
                                        ndisj_inout, disj_cap_inout)) return false;
        }
        return true;
    }
    /* Inequality[v0, op0, v1, op1, ...] — the canonical chained-comparison
     * form produced by the parser. Treat each adjacent (v_i, op_i, v_{i+1})
     * triple as a separate binary comparison and reuse the existing
     * box-or-general classifier on each. */
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol.name == SYM_Inequality
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
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol.name),
                                           pair_args, 2);
            bool ok = fm_collect_constraints(pair, vars, nvars, boxes,
                                             gens_inout, ngens_inout,
                                             gens_cap_inout, disj_inout,
                                             ndisj_inout, disj_cap_inout);
            expr_free(pair);
            if (!ok) return false;
        }
        return true;
    }
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol.name == SYM_Or) {
        /* FindMinimum's gradient penalty method cannot use a non-smooth min
         * penalty, so it passes no sink and Or stays unsupported there. */
        if (!disj_inout) {
            fm_warn(g_fm_name, "nimpl",
                    "disjunctive (Or) constraints are not yet supported");
            return false;
        }
        if (!fm_bool_supported(cons)) {
            fm_warn(g_fm_name, "nimpl", "unsupported disjunctive constraint shape");
            return false;
        }
        if (*ndisj_inout == *disj_cap_inout) {
            size_t nc = *disj_cap_inout ? (*disj_cap_inout) * 2 : 2;
            *disj_inout = (FmDisjunction*)realloc(*disj_inout,
                                                  sizeof(FmDisjunction) * nc);
            *disj_cap_inout = nc;
        }
        (*disj_inout)[*ndisj_inout].expr = expr_copy(cons);
        (*ndisj_inout)++;
        return true;
    }
    /* Reject Element[...] (e.g. x ∈ Integers) outright. */
    if (cons->type == EXPR_FUNCTION
        && cons->data.function.head->type == EXPR_SYMBOL
        && cons->data.function.head->data.symbol.name == SYM_Element) {
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
    (*gens_inout)[*ngens_inout].grad_exprs = NULL;
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

/* Gradient of the augmented objective:
 *   ∇[f + μ · (Σ_k max(0,g_k)^2 + Σ_j h_j^2)]
 * = ∇f + 2μ · Σ_k [equality ? g_k : (g_k > 0 ? g_k : 0)] · ∇g_k.
 *
 * Each constraint contributes only when "active": equalities always; pure
 * inequalities only when violated. Symbolic gradients are used when
 * available (set up by the driver); otherwise per-constraint central
 * differences fill in. The base ∇f gradient mirrors the existing solver
 * code: symbolic g_exprs if non-NULL, else FD on f. */
static bool fm_eval_aug_gradient(Expr* f, Expr** g_exprs,
                                 const FmGenCon* gens, size_t ngens,
                                 double mu,
                                 FmVarBind* binds, const double* x, size_t n,
                                 const FmOpts* opts, double* g_out) {
    bool got = false;
    if (g_exprs) got = fm_eval_gradient(g_exprs, binds, x, n, opts, g_out);
    if (!got) {
        if (!fm_grad_finite_diff(f, binds, x, n, opts, g_out)) return false;
    }
    if (mu <= 0.0 || !gens || ngens == 0) return true;

    double* gk_grad = (double*)malloc(sizeof(double) * n);
    if (!gk_grad) return false;
    for (size_t k = 0; k < ngens; k++) {
        double gk;
        if (!fm_eval_scalar(gens[k].expr, binds, x, n, opts, &gk)) {
            free(gk_grad);
            return false;
        }
        /* Active-set: inequalities contribute only when violated. */
        if (!gens[k].equality && gk <= 0.0) continue;
        bool grad_ok = false;
        if (gens[k].grad_exprs) {
            grad_ok = fm_eval_gradient(gens[k].grad_exprs, binds, x, n, opts, gk_grad);
        }
        if (!grad_ok) {
            grad_ok = fm_grad_finite_diff(gens[k].expr, binds, x, n, opts, gk_grad);
        }
        if (!grad_ok) {
            free(gk_grad);
            return false;
        }
        double scale = 2.0 * mu * gk;
        for (size_t i = 0; i < n; i++) g_out[i] += scale * gk_grad[i];
    }
    free(gk_grad);
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

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Brent's minimisation (1D, MPFR)                                    *
 * ------------------------------------------------------------------ *
 * One-for-one transliteration of fm_brent_min — same parabolic-fit
 * acceptance test, same golden-section fallback, same convergence
 * predicates — with every double replaced by mpfr_t and every constant
 * built fresh from base-10 powers at the requested precision. Box
 * constraints clamp candidate iterates via box->lo / box->hi (those are
 * machine-precision doubles, which is fine: the user's input box
 * already lives at that resolution). */
static bool fm_bracket_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
                            const mpfr_t x0, const FmBox* box1,
                            mpfr_t a, mpfr_t b, mpfr_t c,
                            mpfr_t fa, mpfr_t fb, mpfr_t fc) {
    long bits = opts->wp_bits;
    mpfr_set(a, x0, MPFR_RNDN);
    if (!fm_eval_scalar_mpfr_1d(f, bind, a, opts, fa)) return false;
    /* h = max(|a|, 1) * 1e-2 */
    mpfr_t h;
    mpfr_init2(h, bits);
    mpfr_abs(h, a, MPFR_RNDN);
    {
        mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
        if (mpfr_cmp(h, one) < 0) mpfr_set(h, one, MPFR_RNDN);
        mpfr_clear(one);
    }
    mpfr_mul_d(h, h, 1e-2, MPFR_RNDN);
    mpfr_add(b, a, h, MPFR_RNDN);
    mpfr_clear(h);
    if (box1 && box1->has_hi) {
        mpfr_t bhi; mpfr_init2(bhi, bits); mpfr_set_d(bhi, box1->hi, MPFR_RNDN);
        if (mpfr_cmp(b, bhi) > 0) {
            mpfr_add(b, a, bhi, MPFR_RNDN);
            mpfr_div_ui(b, b, 2, MPFR_RNDN);
        }
        mpfr_clear(bhi);
    }
    if (box1 && box1->has_lo) {
        mpfr_t blo; mpfr_init2(blo, bits); mpfr_set_d(blo, box1->lo, MPFR_RNDN);
        if (mpfr_cmp(b, blo) < 0) {
            mpfr_add(b, a, blo, MPFR_RNDN);
            mpfr_div_ui(b, b, 2, MPFR_RNDN);
        }
        mpfr_clear(blo);
    }
    if (!fm_eval_scalar_mpfr_1d(f, bind, b, opts, fb)) return false;
    if (mpfr_cmp(fb, fa) > 0) { mpfr_swap(a, b); mpfr_swap(fa, fb); }
    /* c = b + 1.618 * (b - a) */
    mpfr_t diff; mpfr_init2(diff, bits);
    mpfr_sub(diff, b, a, MPFR_RNDN);
    mpfr_mul_d(diff, diff, 1.618, MPFR_RNDN);
    mpfr_add(c, b, diff, MPFR_RNDN);
    mpfr_clear(diff);
    if (box1 && box1->has_hi) {
        mpfr_t bhi; mpfr_init2(bhi, bits); mpfr_set_d(bhi, box1->hi, MPFR_RNDN);
        if (mpfr_cmp(c, bhi) > 0) mpfr_set(c, bhi, MPFR_RNDN);
        mpfr_clear(bhi);
    }
    if (box1 && box1->has_lo) {
        mpfr_t blo; mpfr_init2(blo, bits); mpfr_set_d(blo, box1->lo, MPFR_RNDN);
        if (mpfr_cmp(c, blo) < 0) mpfr_set(c, blo, MPFR_RNDN);
        mpfr_clear(blo);
    }
    if (!fm_eval_scalar_mpfr_1d(f, bind, c, opts, fc)) return false;
    for (int k = 0; k < 100 && mpfr_cmp(fc, fb) <= 0; k++) {
        mpfr_swap(a, b); mpfr_swap(fa, fb);
        mpfr_swap(b, c); mpfr_swap(fb, fc);
        mpfr_init2(diff, bits);
        mpfr_sub(diff, b, a, MPFR_RNDN);
        mpfr_mul_d(diff, diff, 1.618, MPFR_RNDN);
        mpfr_add(c, b, diff, MPFR_RNDN);
        mpfr_clear(diff);
        bool hit_bound = false;
        if (box1 && box1->has_hi) {
            mpfr_t bhi; mpfr_init2(bhi, bits); mpfr_set_d(bhi, box1->hi, MPFR_RNDN);
            if (mpfr_cmp(c, bhi) >= 0) { mpfr_set(c, bhi, MPFR_RNDN); hit_bound = true; }
            mpfr_clear(bhi);
        }
        if (box1 && box1->has_lo) {
            mpfr_t blo; mpfr_init2(blo, bits); mpfr_set_d(blo, box1->lo, MPFR_RNDN);
            if (mpfr_cmp(c, blo) <= 0) { mpfr_set(c, blo, MPFR_RNDN); hit_bound = true; }
            mpfr_clear(blo);
        }
        if (!fm_eval_scalar_mpfr_1d(f, bind, c, opts, fc)) return false;
        if (hit_bound) break;
    }
    if (mpfr_cmp(a, c) > 0) { mpfr_swap(a, c); mpfr_swap(fa, fc); }
    return true;
}

static bool fm_brent_min_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
                              const mpfr_t a_in, const mpfr_t b_in, const mpfr_t c_in,
                              const FmBox* box1,
                              mpfr_t x_out, mpfr_t fx_out) {
    long bits = opts->wp_bits;
    mpfr_t a, c; mpfr_init2(a, bits); mpfr_init2(c, bits);
    if (mpfr_cmp(a_in, c_in) <= 0) { mpfr_set(a, a_in, MPFR_RNDN); mpfr_set(c, c_in, MPFR_RNDN); }
    else                            { mpfr_set(a, c_in, MPFR_RNDN); mpfr_set(c, a_in, MPFR_RNDN); }
    mpfr_t tol, tol_acc, zeps;
    mpfr_init2(tol, bits); mpfr_init2(tol_acc, bits); mpfr_init2(zeps, bits);
    fm_tol_from_digits(tol, opts->prec_goal_digits);
    fm_tol_from_digits(tol_acc, opts->acc_goal_digits);
    /* zeps tracks the MPFR working precision rather than a fixed 1e-12,
     * so x's last few representable bits are still allowed to settle. */
    long zdig = (long)numeric_bits_to_digits(bits) - 1;
    if (zdig < 1) zdig = 1;
    mpfr_set_ui(zeps, 10, MPFR_RNDN);
    mpfr_pow_si(zeps, zeps, -zdig, MPFR_RNDN);

    mpfr_t e_step, d, x, w, v, fx, fw, fv;
    mpfr_t xm, tol1, tol2, u, fu, p, q, r, etemp;
    mpfr_init2(e_step, bits); mpfr_init2(d, bits);
    mpfr_init2(x, bits); mpfr_init2(w, bits); mpfr_init2(v, bits);
    mpfr_init2(fx, bits); mpfr_init2(fw, bits); mpfr_init2(fv, bits);
    mpfr_init2(xm, bits); mpfr_init2(tol1, bits); mpfr_init2(tol2, bits);
    mpfr_init2(u, bits); mpfr_init2(fu, bits);
    mpfr_init2(p, bits); mpfr_init2(q, bits); mpfr_init2(r, bits);
    mpfr_init2(etemp, bits);

    mpfr_set_zero(e_step, +1);
    mpfr_set_zero(d, +1);
    mpfr_set(x, b_in, MPFR_RNDN);
    mpfr_set(w, x, MPFR_RNDN);
    mpfr_set(v, x, MPFR_RNDN);
    if (mpfr_cmp(x, a) < 0 || mpfr_cmp(x, c) > 0) {
        mpfr_add(x, a, c, MPFR_RNDN); mpfr_div_ui(x, x, 2, MPFR_RNDN);
        mpfr_set(w, x, MPFR_RNDN); mpfr_set(v, x, MPFR_RNDN);
    }

    bool ok = false;
    bool converged = false;
    if (!fm_eval_scalar_mpfr_1d(f, bind, x, opts, fx)) goto cleanup;
    mpfr_set(fw, fx, MPFR_RNDN);
    mpfr_set(fv, fx, MPFR_RNDN);

    mpfr_t diff_xm, half_ca, crit, abs_fx, thresh;
    mpfr_init2(diff_xm, bits); mpfr_init2(half_ca, bits); mpfr_init2(crit, bits);
    mpfr_init2(abs_fx, bits); mpfr_init2(thresh, bits);

    for (int64_t k = 0; k < opts->max_iter; k++) {
        mpfr_add(xm, a, c, MPFR_RNDN); mpfr_div_ui(xm, xm, 2, MPFR_RNDN);
        mpfr_abs(tol1, x, MPFR_RNDN); mpfr_mul(tol1, tol1, tol, MPFR_RNDN);
        mpfr_add(tol1, tol1, zeps, MPFR_RNDN);
        mpfr_mul_ui(tol2, tol1, 2, MPFR_RNDN);
        /* convergence: |x - xm| <= tol2 - (c - a) / 2 */
        mpfr_sub(diff_xm, x, xm, MPFR_RNDN); mpfr_abs(diff_xm, diff_xm, MPFR_RNDN);
        mpfr_sub(half_ca, c, a, MPFR_RNDN);  mpfr_div_ui(half_ca, half_ca, 2, MPFR_RNDN);
        mpfr_sub(crit, tol2, half_ca, MPFR_RNDN);
        bool cvg1 = (mpfr_cmp(diff_xm, crit) <= 0);
        mpfr_abs(abs_fx, fx, MPFR_RNDN);
        mpfr_add_ui(thresh, abs_fx, 1, MPFR_RNDN);
        mpfr_mul(thresh, thresh, tol_acc, MPFR_RNDN);
        bool cvg2 = (mpfr_cmp(abs_fx, thresh) < 0);
        if (cvg1 || cvg2) { converged = true; break; }

        bool used_parabolic = false;
        if (mpfr_cmpabs(e_step, tol1) > 0) {
            mpfr_t xw, xv_d, fxfv, fxfw, t1, t2;
            mpfr_init2(xw, bits); mpfr_init2(xv_d, bits);
            mpfr_init2(fxfv, bits); mpfr_init2(fxfw, bits);
            mpfr_init2(t1, bits); mpfr_init2(t2, bits);
            mpfr_sub(xw, x, w, MPFR_RNDN);
            mpfr_sub(xv_d, x, v, MPFR_RNDN);
            mpfr_sub(fxfv, fx, fv, MPFR_RNDN);
            mpfr_sub(fxfw, fx, fw, MPFR_RNDN);
            mpfr_mul(r, xw, fxfv, MPFR_RNDN);
            mpfr_mul(q, xv_d, fxfw, MPFR_RNDN);
            mpfr_mul(t1, xv_d, q, MPFR_RNDN);
            mpfr_mul(t2, xw, r, MPFR_RNDN);
            mpfr_sub(p, t1, t2, MPFR_RNDN);
            mpfr_sub(q, q, r, MPFR_RNDN);
            mpfr_mul_ui(q, q, 2, MPFR_RNDN);
            if (mpfr_sgn(q) > 0) mpfr_neg(p, p, MPFR_RNDN);
            mpfr_abs(q, q, MPFR_RNDN);
            mpfr_set(etemp, e_step, MPFR_RNDN);
            mpfr_set(e_step, d, MPFR_RNDN);
            mpfr_t lower, upper, abs_p, halfq_e;
            mpfr_init2(lower, bits); mpfr_init2(upper, bits);
            mpfr_init2(abs_p, bits); mpfr_init2(halfq_e, bits);
            mpfr_sub(lower, a, x, MPFR_RNDN); mpfr_mul(lower, lower, q, MPFR_RNDN);
            mpfr_sub(upper, c, x, MPFR_RNDN); mpfr_mul(upper, upper, q, MPFR_RNDN);
            mpfr_abs(abs_p, p, MPFR_RNDN);
            mpfr_mul(halfq_e, q, etemp, MPFR_RNDN); mpfr_abs(halfq_e, halfq_e, MPFR_RNDN);
            mpfr_div_ui(halfq_e, halfq_e, 2, MPFR_RNDN);
            bool reject = (mpfr_cmp(abs_p, halfq_e) >= 0
                        || mpfr_cmp(p, lower) <= 0
                        || mpfr_cmp(p, upper) >= 0);
            if (!reject) {
                mpfr_div(d, p, q, MPFR_RNDN);
                mpfr_add(u, x, d, MPFR_RNDN);
                mpfr_t ua, cu;
                mpfr_init2(ua, bits); mpfr_init2(cu, bits);
                mpfr_sub(ua, u, a, MPFR_RNDN);
                mpfr_sub(cu, c, u, MPFR_RNDN);
                if (mpfr_cmp(ua, tol2) < 0 || mpfr_cmp(cu, tol2) < 0) {
                    int s = (mpfr_cmp(xm, x) >= 0) ? 1 : -1;
                    if (s > 0) mpfr_set(d, tol1, MPFR_RNDN);
                    else       mpfr_neg(d, tol1, MPFR_RNDN);
                }
                mpfr_clears(ua, cu, (mpfr_ptr)0);
                used_parabolic = true;
            }
            mpfr_clears(xw, xv_d, fxfv, fxfw, t1, t2,
                        lower, upper, abs_p, halfq_e, (mpfr_ptr)0);
        }
        if (!used_parabolic) {
            if (mpfr_cmp(x, xm) >= 0) mpfr_sub(e_step, a, x, MPFR_RNDN);
            else                       mpfr_sub(e_step, c, x, MPFR_RNDN);
            mpfr_mul_d(d, e_step, FM_CGOLD, MPFR_RNDN);
        }
        if (mpfr_cmpabs(d, tol1) >= 0) { mpfr_add(u, x, d, MPFR_RNDN); }
        else if (mpfr_sgn(d) >= 0)     { mpfr_add(u, x, tol1, MPFR_RNDN); }
        else                            { mpfr_sub(u, x, tol1, MPFR_RNDN); }
        if (box1) {
            if (box1->has_lo) {
                mpfr_t blo; mpfr_init2(blo, bits); mpfr_set_d(blo, box1->lo, MPFR_RNDN);
                if (mpfr_cmp(u, blo) < 0) mpfr_set(u, blo, MPFR_RNDN);
                mpfr_clear(blo);
            }
            if (box1->has_hi) {
                mpfr_t bhi; mpfr_init2(bhi, bits); mpfr_set_d(bhi, box1->hi, MPFR_RNDN);
                if (mpfr_cmp(u, bhi) > 0) mpfr_set(u, bhi, MPFR_RNDN);
                mpfr_clear(bhi);
            }
        }
        if (!fm_eval_scalar_mpfr_1d(f, bind, u, opts, fu)) goto cleanup_inner;
        fm_fire_monitor(opts->step_monitor);
        if (mpfr_cmp(fu, fx) <= 0) {
            if (mpfr_cmp(u, x) >= 0) mpfr_set(a, x, MPFR_RNDN);
            else                      mpfr_set(c, x, MPFR_RNDN);
            mpfr_set(v, w, MPFR_RNDN); mpfr_set(w, x, MPFR_RNDN); mpfr_set(x, u, MPFR_RNDN);
            mpfr_set(fv, fw, MPFR_RNDN); mpfr_set(fw, fx, MPFR_RNDN); mpfr_set(fx, fu, MPFR_RNDN);
        } else {
            if (mpfr_cmp(u, x) < 0) mpfr_set(a, u, MPFR_RNDN);
            else                     mpfr_set(c, u, MPFR_RNDN);
            if (mpfr_cmp(fu, fw) <= 0 || mpfr_equal_p(w, x)) {
                mpfr_set(v, w, MPFR_RNDN); mpfr_set(w, u, MPFR_RNDN);
                mpfr_set(fv, fw, MPFR_RNDN); mpfr_set(fw, fu, MPFR_RNDN);
            } else if (mpfr_cmp(fu, fv) <= 0 || mpfr_equal_p(v, x) || mpfr_equal_p(v, w)) {
                mpfr_set(v, u, MPFR_RNDN); mpfr_set(fv, fu, MPFR_RNDN);
            }
        }
    }
    if (!converged) {
        fm_warn(g_fm_name, "cvmit",
                "Brent (MPFR) failed to converge within %lld iterations",
                (long long)opts->max_iter);
    }
    mpfr_set(x_out, x, MPFR_RNDN);
    mpfr_set(fx_out, fx, MPFR_RNDN);
    ok = true;
cleanup_inner:
    mpfr_clears(diff_xm, half_ca, crit, abs_fx, thresh, (mpfr_ptr)0);
cleanup:
    mpfr_clears(a, c, tol, tol_acc, zeps,
                e_step, d, x, w, v, fx, fw, fv,
                xm, tol1, tol2, u, fu,
                p, q, r, etemp, (mpfr_ptr)0);
    return ok;
}
#endif /* USE_MPFR */

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

    bool got_grad;
    if (augmented) {
        got_grad = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                        binds, x, n, opts, g);
    } else {
        got_grad = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
        if (!got_grad) got_grad = fm_grad_finite_diff(f, binds, x, n, opts, g);
    }
    if (!got_grad) {
        fm_warn(g_fm_name, "nlnum", "gradient evaluation failed at start point");
        goto cleanup;
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
            /* Line-search exhaustion is expected at high μ in the penalty
             * schedule (steep walls, large directional curvature). The
             * outer fm_run_penalty loop's feasibility check is the
             * authoritative signal in that case, so stay silent here and
             * let it speak instead. */
            if (!augmented) {
                fm_warn(g_fm_name, "lstol", "line search failed at iter %lld",
                        (long long)k);
            }
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
        bool ng_ok;
        if (augmented) {
            ng_ok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                         binds, x_new, n, opts, g_new);
        } else {
            ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
            if (!ng_ok) ng_ok = fm_grad_finite_diff(f, binds, x_new, n, opts, g_new);
        }
        if (!ng_ok) {
            fm_warn(g_fm_name, "nlnum", "gradient evaluation failed in iteration");
            /* Take the step and stop. */
            for (size_t i = 0; i < n; i++) x[i] = x_new[i];
            fx = fx_new;
            break;
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
    bool got_grad;
    if (augmented) {
        got_grad = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                        binds, x, n, opts, g);
    } else {
        got_grad = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
        if (!got_grad) got_grad = fm_grad_finite_diff(f, binds, x, n, opts, g);
    }
    if (!got_grad) {
        fm_warn(g_fm_name, "nlnum", "gradient failed at start point");
        goto cleanup;
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
        if (!ls_ok) {
            if (!augmented) fm_warn(g_fm_name, "lstol", "line search failed");
            break;
        }
        fm_fire_monitor(opts->step_monitor);

        bool ng_ok;
        if (augmented) {
            ng_ok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                         binds, x_new, n, opts, g_new);
        } else {
            ng_ok = g_exprs && fm_eval_gradient(g_exprs, binds, x_new, n, opts, g_new);
            if (!ng_ok) ng_ok = fm_grad_finite_diff(f, binds, x_new, n, opts, g_new);
        }
        if (!ng_ok) {
            for (size_t i = 0; i < n; i++) x[i] = x_new[i];
            fx = fx_new;
            break;
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
        bool gok;
        if (augmented) {
            gok = fm_eval_aug_gradient(f, g_exprs, gens, ngens, mu,
                                       binds, x, n, opts, g);
        } else {
            gok = g_exprs && fm_eval_gradient(g_exprs, binds, x, n, opts, g);
            if (!gok) gok = fm_grad_finite_diff(f, binds, x, n, opts, g);
        }
        if (!gok) {
            fm_warn(g_fm_name, "nlnum", "gradient failed during Newton");
            goto cleanup;
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
        if (!ls_ok) {
            if (!augmented) fm_warn(g_fm_name, "lstol", "Newton line search failed");
            break;
        }
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
    /* Outer μ schedule: 1 → 10 → ... up to 10^8.
     *
     * The classical penalty-method convergence theorem requires μ → ∞:
     * for problems with an active constraint at the optimum, the unaugmented
     * gradient is non-zero on the constraint surface, so the augmented
     * minimizer sits a 1/μ-sized step inside the infeasible region. Stopping
     * as soon as the iterate becomes feasible (the previous behaviour) can
     * therefore terminate well before the constrained optimum — e.g. with
     *
     *   FindMinimum[{x + y, 3 x + 2 y >= 7 && x >= 0 && y >= 0}, {x, y}]
     *
     * the BFGS inner solver lands at (x ≈ 2.45, 0) once it becomes feasible,
     * while the true minimum is (7/3, 0). Continuing to ramp μ pulls the
     * iterate down the active boundary to the corner.
     *
     * Termination: stop once both (a) the iterate is feasible and (b)
     * inter-round movement is small (i.e. the augmented minimum has
     * stabilised), or the schedule is exhausted. */
    double mu = 1.0;
    double fx = 0.0;
    bool feas = false;
    double* x_prev = (double*)malloc(sizeof(double) * n);
    if (!x_prev) return false;
    const double rel_tol = pow(10.0, -opts->prec_goal_digits);
    for (int round = 0; round < 9; round++) {
        for (size_t i = 0; i < n; i++) x_prev[i] = x[i];
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
        if (!ok) { free(x_prev); return false; }
        double pen;
        if (!fm_eval_penalty(gens, ngens, binds, x, n, opts, &pen)) {
            free(x_prev); return false;
        }
        if (pen < 1e-12) feas = true;
        if (feas) {
            /* Stop only if the iterate has stabilised between rounds. */
            double max_step = 0.0, max_x = 0.0;
            for (size_t i = 0; i < n; i++) {
                double ds = fabs(x[i] - x_prev[i]);
                if (ds > max_step) max_step = ds;
                double xa = fabs(x[i]);
                if (xa > max_x) max_x = xa;
            }
            if (max_step < rel_tol * (max_x + 1.0)) break;
        }
        mu *= 10.0;
    }
    free(x_prev);
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
    opts.max_iter_set = false;
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
        && f_arg->data.function.head->data.symbol.name == SYM_List
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
        && var_arg->data.function.head->data.symbol.name == SYM_List
        && var_arg->data.function.arg_count > 0) {
        size_t na = var_arg->data.function.arg_count;
        bool any_inner = false, all_inner_or_sym = true;
        for (size_t i = 0; i < na; i++) {
            Expr* e = var_arg->data.function.args[i];
            bool is_inner = (e->type == EXPR_FUNCTION
                && e->data.function.head->type == EXPR_SYMBOL
                && e->data.function.head->data.symbol.name == SYM_List);
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
    CompiledProgram*  f_prog = NULL;      /* machine-precision compiled objective */
    CompiledProgram** grad_progs = NULL;  /* per-component compiled exact gradient */
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
    for (size_t i = 0; i < n; i++) fm_bind_snapshot(&binds[i], vars[i]->data.symbol.name);

    /* Constraints. FindMinimum passes no disjunction sink: its smooth gradient
     * penalty method cannot use the non-smooth min penalty an Or requires. */
    if (cons) {
        if (!fm_collect_constraints(cons, vars, n, boxes, &gens, &ngens, &gcap,
                                    NULL, NULL, NULL))
            goto cleanup;
        /* Best-effort symbolic gradient of each constraint expression. The
         * penalty solver needs ∇(f + μ·Σ penalty) — using a stale ∇f alone
         * gives the inner BFGS/CG an inconsistent value/gradient pair and
         * the penalty term loses all influence on the descent direction. */
        for (size_t k = 0; k < ngens; k++) {
            gens[k].grad_exprs = fm_compute_gradient(gens[k].expr, vars, n);
            /* NULL is fine — fm_eval_aug_gradient will FD that constraint. */
        }
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
            && opts.gradient->data.function.head->data.symbol.name == SYM_List
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

    /* Machine-precision auto-compilation of the objective. The local solvers
     * evaluate f at every trial point (line search, bracketing, function
     * values); lowering it to bytecode once over the variables and running the
     * register machine per point is far cheaper than the interpreter
     * (expr_copy + n OwnValue installs + evaluate + numericalize). Registered in
     * the g_fm_obj_* slots that fm_eval_scalar consults; deregistered and freed
     * at cleanup. Compiled *after* fm_bind_snapshot (which cleared the vars'
     * OwnValues) so they lower as the argument symbols, not folded constants.
     * The symbolic gradient is left exact — only the value path is compiled — so
     * FindMinimum's precision is unchanged. MPFR keeps the exact interpreter
     * path (its solvers never call fm_eval_scalar). A body Compile can't lower
     * stays NULL and the interpreter is used, and every per-point call falls
     * back on a non-finite compiled result, so this is a pure speedup.
     *
     * The exact symbolic gradient `g_exprs` is compiled the same way (each
     * component is a function of all the variables). This is what actually
     * accelerates the QuasiNewton/CG/Newton loop, whose cost is dominated by the
     * per-iteration gradient — the same ∂f/∂x_i, lowered, so the gradient stays
     * exact (no finite differences) and FindMinimum's precision is unchanged.
     * A component Compile can't lower stays NULL and falls back per-component.
     *
     * Skipped when an "EvaluationMonitor" is set: the monitor fires inside
     * fm_eval_with_bindings (per interpreter evaluation of f), which the compiled
     * path bypasses, so compiling would silently stop the monitor from firing.
     * Monitoring is a debugging aid, not a performance path, so falling back to
     * the interpreter there is the right trade. */
    if (opts.prec_mode == FM_PREC_MACHINE && n > 0 && !opts.eval_monitor) {
        const char** cnames = (const char**)malloc(sizeof(char*) * n);
        CompileType* ctypes = (CompileType*)malloc(sizeof(CompileType) * n);
        for (size_t i = 0; i < n; i++) { cnames[i] = vars[i]->data.symbol.name; ctypes[i] = CT_REAL; }
        f_prog = compile_expr_ex(f_raw, cnames, ctypes, n, COMPILE_FOLD_GLOBALS);
        if (f_prog && compiled_result_type(f_prog) != CT_REAL) { compiled_free(f_prog); f_prog = NULL; }
        if (g_exprs) {
            grad_progs = (CompiledProgram**)calloc(n, sizeof(CompiledProgram*));
            for (size_t i = 0; i < n; i++) {
                CompiledProgram* p = compile_expr_ex(g_exprs[i], cnames, ctypes, n,
                                                     COMPILE_FOLD_GLOBALS);
                if (p && compiled_result_type(p) != CT_REAL) { compiled_free(p); p = NULL; }
                grad_progs[i] = p;
            }
        }
        free(cnames); free(ctypes);
    }
    g_fm_obj_expr = f_raw; g_fm_obj_prog = f_prog; g_fm_obj_nargs = n;
    g_fm_grad_exprs = g_exprs; g_fm_grad_progs = grad_progs; g_fm_grad_n = n;

    /* Dispatch. */
    double fx_min = 0.0;
    bool ok = true;
    bool has_general_cons = (ngens > 0);
#ifdef USE_MPFR
    bool mpfr_result = false;
    mpfr_t* x_vec_mpfr = NULL;
    mpfr_t fx_min_mpfr;
    bool use_mpfr = (opts.prec_mode == FM_PREC_MPFR);
    /* Penalty path is not lifted to MPFR yet — fall back to machine
     * precision in that case rather than silently dropping the constraint. */
    if (use_mpfr && has_general_cons) {
        fm_warn(fn_name, "nimpl",
                "general (non-box) constraints at WorkingPrecision > MachinePrecision "
                "are not yet supported; falling back to machine precision");
        use_mpfr = false;
    }
    if (use_mpfr) {
        mpfr_init2(fx_min_mpfr, opts.wp_bits);
        x_vec_mpfr = fm_mpfr_array(n, opts.wp_bits);
        for (size_t i = 0; i < n; i++) mpfr_set_d(x_vec_mpfr[i], x_vec[i], MPFR_RNDN);
        if (method == FM_METHOD_BRENT) {
            if (n != 1) {
                fm_warn(fn_name, "badmeth", "Method \"Brent\" requires a single variable");
                ok = false; goto run_done;
            }
            mpfr_t a_m, b_m, c_m, fa_m, fb_m, fc_m;
            mpfr_init2(a_m, opts.wp_bits); mpfr_init2(b_m, opts.wp_bits);
            mpfr_init2(c_m, opts.wp_bits);
            mpfr_init2(fa_m, opts.wp_bits); mpfr_init2(fb_m, opts.wp_bits);
            mpfr_init2(fc_m, opts.wp_bits);
            bool bracketed = false;
            if (boxes[0].has_lo && boxes[0].has_hi) {
                mpfr_set_d(a_m, boxes[0].lo, MPFR_RNDN);
                mpfr_set_d(c_m, boxes[0].hi, MPFR_RNDN);
                mpfr_add(b_m, a_m, c_m, MPFR_RNDN); mpfr_div_ui(b_m, b_m, 2, MPFR_RNDN);
                /* If user gave a start inside, use it. */
                if (x_vec[0] > boxes[0].lo && x_vec[0] < boxes[0].hi)
                    mpfr_set(b_m, x_vec_mpfr[0], MPFR_RNDN);
                bracketed = true;
            } else {
                bracketed = fm_bracket_mpfr(f_raw, binds, &opts, x_vec_mpfr[0],
                                            &boxes[0], a_m, b_m, c_m, fa_m, fb_m, fc_m);
                if (!bracketed) fm_warn(fn_name, "nlnum", "MPFR bracket-finding failed");
            }
            if (bracketed) {
                mpfr_t xm_m, fmin_m;
                mpfr_init2(xm_m, opts.wp_bits); mpfr_init2(fmin_m, opts.wp_bits);
                ok = fm_brent_min_mpfr(f_raw, binds, &opts, a_m, b_m, c_m,
                                       &boxes[0], xm_m, fmin_m);
                if (ok) {
                    mpfr_set(x_vec_mpfr[0], xm_m, MPFR_RNDN);
                    mpfr_set(fx_min_mpfr, fmin_m, MPFR_RNDN);
                }
                mpfr_clears(xm_m, fmin_m, (mpfr_ptr)0);
            } else {
                ok = false;
            }
            mpfr_clears(a_m, b_m, c_m, fa_m, fb_m, fc_m, (mpfr_ptr)0);
        } else {
            /* n-D path: BFGS handles QuasiNewton; Newton/CG fall back to
             * BFGS at MPFR with a one-shot diagnostic. */
            if (method == FM_METHOD_NEWTON || method == FM_METHOD_CONJGRAD) {
                fm_warn(fn_name, "nimpl",
                        "Method \"%s\" at WorkingPrecision > MachinePrecision is not yet "
                        "supported; falling back to QuasiNewton",
                        method == FM_METHOD_NEWTON ? "Newton" : "ConjugateGradient");
            }
            ok = fm_run_bfgs_mpfr(f_raw, vars, n, binds, g_exprs,
                                  x_vec_mpfr, boxes, &opts, fx_min_mpfr);
        }
        mpfr_result = ok;
    } else {
#endif
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
#ifdef USE_MPFR
    }
#endif
run_done:
    /* Clear temp bindings first so the variable symbol stays free during
     * Rule construction (otherwise `Rule[x, v]` would re-evaluate x to its
     * pre-call value once we restore). */
    if (binds) {
        for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);
    }
    if (ok) {
#ifdef USE_MPFR
        if (mpfr_result) {
            result_out = fm_build_result_mpfr(fx_min_mpfr, vars,
                                              (mpfr_t const*)x_vec_mpfr, n);
        } else {
            result_out = fm_build_result(fx_min, vars, x_vec, n);
        }
#else
        result_out = fm_build_result(fx_min, vars, x_vec, n);
#endif
    }
#ifdef USE_MPFR
    if (x_vec_mpfr) fm_mpfr_array_free(x_vec_mpfr, n);
    if (use_mpfr)   mpfr_clear(fx_min_mpfr);
#endif

cleanup:
    g_fm_obj_expr = NULL;      /* deregister objective + gradient before freeing */
    g_fm_obj_prog = NULL;
    g_fm_obj_nargs = 0;
    g_fm_grad_exprs = NULL;
    g_fm_grad_progs = NULL;
    g_fm_grad_n = 0;
    if (f_prog) compiled_free(f_prog);
    if (grad_progs) {
        for (size_t i = 0; i < n; i++) if (grad_progs[i]) compiled_free(grad_progs[i]);
        free(grad_progs);
    }
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
        for (size_t k = 0; k < ngens; k++) {
            expr_free(gens[k].expr);
            if (gens[k].grad_exprs) {
                for (size_t i = 0; i < n; i++) expr_free(gens[k].grad_exprs[i]);
                free(gens[k].grad_exprs);
            }
        }
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
        && f_orig->data.function.head->data.symbol.name == SYM_List
        && f_orig->data.function.arg_count == 2) {
        /* Wrap inner f only. */
        Expr* inner_f = f_orig->data.function.args[0];
        Expr* cons = f_orig->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(inner_f) };
        neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        Expr* list_args[2] = { neg_f, expr_copy(cons) };
        new_first = expr_new_function(expr_new_symbol(SYM_List), list_args, 2);
    } else {
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(f_orig) };
        neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        new_first = neg_f;
    }
    /* Construct synthetic FindMinimum[new_first, vars, opts...]. */
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * argc);
    new_args[0] = new_first;
    for (size_t i = 1; i < argc; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    Expr* synthetic = expr_new_function(expr_new_symbol(SYM_FindMinimum), new_args, argc);
    free(new_args);
    /* Drive findmin directly so the diagnostic tag is FindMaximum. */
    Expr* min_result = findmin_driver(synthetic, "FindMaximum");
    expr_free(synthetic);
    if (!min_result) return NULL;
    /* min_result is {fmin, {rules}}; negate fmin while preserving its
     * numeric type so a WorkingPrecision -> N run keeps the N-digit MPFR
     * head instead of collapsing back to machine precision. */
    if (min_result->type == EXPR_FUNCTION
        && min_result->data.function.arg_count == 2) {
        Expr* fmin_e = min_result->data.function.args[0];
#ifdef USE_MPFR
        if (fmin_e && fmin_e->type == EXPR_MPFR) {
            long bits = mpfr_get_prec(fmin_e->data.mpfr);
            mpfr_t neg; mpfr_init2(neg, bits);
            mpfr_neg(neg, fmin_e->data.mpfr, MPFR_RNDN);
            expr_free(fmin_e);
            min_result->data.function.args[0] = expr_new_mpfr_copy(neg);
            mpfr_clear(neg);
        } else
#endif
        {
            double fmin;
            if (fm_expr_to_double_real(fmin_e, &fmin)) {
                expr_free(fmin_e);
                min_result->data.function.args[0] = expr_new_real(-fmin);
            }
        }
    }
    return min_result;
}

/* ================================================================== *
 *  NMinimize / NMaximize — numerical GLOBAL optimization             *
 * ================================================================== *
 *
 *  Layered directly on the FindMinimum machinery above. NMinimize reuses
 *  the same Block-style variable binding (fm_bind_*), objective and
 *  constraint evaluation (fm_eval_scalar / fm_eval_penalty), constraint
 *  parsing (fm_collect_constraints → boxes + general FmGenCon[] +
 *  disjunctions), local polishers (fm_run_bfgs / fm_run_penalty, plus
 *  fm_run_bfgs_mpfr for WorkingPrecision), and the result builders. What is new
 *  here is the stochastic *global* search — DifferentialEvolution (Automatic
 *  default), NelderMead, RandomSearch, SimulatedAnnealing — mixed-integer
 *  handling via Element[x, Integers], and the infeasible {Infinity, ...} return.
 *
 *  Constraints are handled with Deb's feasibility rules during the global
 *  search (no penalty-weight tuning): a feasible point always beats an
 *  infeasible one; among feasible points the smaller objective wins; among
 *  infeasible points the smaller total violation wins. The global best is
 *  then polished with the exact local solver and re-checked for feasibility.
 *
 *  Disjunctive (Or) constraints — feasible iff at least one branch holds — are
 *  supported here (but not in FindMinimum's smooth gradient path): each adds its
 *  minimum-branch penalty (fm_bool_penalty), 0 when any branch is satisfied, so
 *  the same Deb gate selects it with no extra tuning. The derivative-free global
 *  search consumes the non-smooth min directly; the local polish optimises only
 *  the conjunctive part and the post-polish feasibility gate keeps the result
 *  inside the disjunctive-feasible region.
 */

#define NM_DEFAULT_SPAN   10.0     /* half-width of the default search box   */
#define NM_BOUND_SPAN     20.0     /* span added when only one bound is known*/
#define NM_FEAS_EPS       1.0e-8   /* penalty ≤ this ⇒ feasible (selection)  */
#define NM_FEAS_FINAL     1.0e-6   /* final feasible-vs-Infinity threshold   */
#define NM_PENALTY_MU     1.0e6    /* fixed penalty weight for NelderMead/SA */
#define NM_SA_TOTAL_CAP   120000   /* SA aggregate iteration cap across chains
                                    * when "SearchPoints" -> K > 1 restarts     */
#define NM_SA_BURN_IN     30       /* trial steps probed per chain to measure
                                    * the objective scale for the temperature   */
#define NM_DEFAULT_SEED   20260814ULL
#define NM_MAX_REGION_EXPAND 4     /* infeasible-region rescue: grow +-SPAN by
                                    * 10^k for k = 1..this (up to +-1e5)       */

enum { NM_AUTO = 0, NM_DE, NM_NELDERMEAD, NM_RANDOMSEARCH, NM_SA };

typedef struct {
    int      method;         /* NM_AUTO / NM_DE / ...                        */
    int      search_points;  /* 0 ⇒ auto                                     */
    double   F;              /* DE scaling factor;   <0 ⇒ auto               */
    double   CR;             /* DE crossover prob.;  <0 ⇒ auto               */
    double   reflect_ratio;  /* NelderMead reflection coeff;  <0 ⇒ default 1  */
    double   expand_ratio;   /* NelderMead expansion coeff;   <0 ⇒ default 2  */
    double   contract_ratio; /* NelderMead contraction coeff; <0 ⇒ default .5 */
    double   shrink_ratio;   /* NelderMead shrink coeff;      <0 ⇒ default .5 */
    double   tolerance;      /* simplex convergence tolerance; <0 ⇒ default    */
    int      post_process;   /* -1 auto (on) / 1 on / 0 off (skip polish)    */
    Expr*    init_points;    /* "InitialPoints" -> {{...},...}, borrowed / NULL */
    Expr*    penalty_fn;     /* "PenaltyFunction" -> f, borrowed / NULL ⇒ auto */
    double   perturb_scale;  /* SA "PerturbationScale"; <0 ⇒ default 1.0       */
    Expr*    boltzmann_fn;   /* SA "BoltzmannExponent" -> f, borrowed / NULL ⇒ auto */
    int      level_iterations;/* SA "LevelIterations"; <=0 ⇒ auto (50)         */
    uint64_t seed;
} NmConfig;

/* SplitMix64 — a small, fast, fully deterministic PRNG (mirrors mcint). A
 * fixed default seed makes the stochastic search reproducible so the unit
 * tests are stable; "RandomSeed" overrides it. */
typedef struct { uint64_t s; } NmRng;

static void nm_rng_seed(NmRng* r, uint64_t seed) { r->s = seed; }

static uint64_t nm_rng_next(NmRng* r) {
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double nm_rng_unif(NmRng* r) {            /* [0, 1)                    */
    return (double)(nm_rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
static double nm_rng_range(NmRng* r, double lo, double hi) {
    return lo + (hi - lo) * nm_rng_unif(r);
}
static double nm_rng_normal(NmRng* r) {          /* standard normal (Box-Muller) */
    double u1 = nm_rng_unif(r), u2 = nm_rng_unif(r);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(6.28318530717958647692 * u2);
}

static bool nm_is_head(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

/* A variable "atom" is a bare symbol (x) or an indexed form (x[i], x[i,j], ...):
 * anything a Table[x[i], {i, ...}] spec can produce. A {x, lo, hi} spec (head
 * List) is NOT an atom — it is a bounded-variable spec handled separately. */
static bool nm_is_var_atom(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return true;
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name != SYM_List
        && e->data.function.head->data.symbol.name != SYM_Element;
}

/* Is `cons` already a boolean/relational constraint tree, or a wrapper (Table,
 * List, ...) that must be evaluated first to expand into one? */
static bool nm_is_constraint_tree(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    return h == SYM_And || h == SYM_Or || h == SYM_Not || h == SYM_Xor
        || h == SYM_Implies || h == SYM_Equal || h == SYM_Unequal
        || h == SYM_Less || h == SYM_LessEqual || h == SYM_Greater
        || h == SYM_GreaterEqual || h == SYM_Inequality || h == SYM_Element;
}

/* Structural substitution: a fresh copy of `e` with every subtree structurally
 * equal to from[i] replaced by a copy of to[i]. Used to rewrite indexed
 * optimisation variables (x[1], x[2], ...) to fresh scalar symbols so the whole
 * symbol-keyed solver machinery applies unchanged. */
static Expr* nm_subst(Expr* e, Expr* const* from, Expr* const* to, size_t n) {
    if (!e) return NULL;
    for (size_t i = 0; i < n; i++)
        if (from[i] && expr_eq(e, from[i])) return expr_copy(to[i]);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t ac = e->data.function.arg_count;
    Expr* head = nm_subst(e->data.function.head, from, to, n);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t j = 0; j < ac; j++)
        args[j] = nm_subst(e->data.function.args[j], from, to, n);
    Expr* r = expr_new_function(head, args, ac);
    free(args);
    return r;
}

/* Generate a fresh, unused scalar symbol to stand in for an indexed variable.
 * The name is interned and stable; the caller removes it (symtab_remove_symbol)
 * once the optimisation and result construction are done. */
static Expr* nm_fresh_symbol(void) {
    static uint64_t ctr = 0;
    char buf[64];
    for (;;) {
        snprintf(buf, sizeof(buf), "NMinimize$%llu", (unsigned long long)ctr++);
        if (!symtab_lookup(buf)) break;
    }
    return expr_new_symbol(buf);
}

/* Block-localise a set of head symbols: snapshot and clear their Own/DownValues
 * so an objective/constraint/variable expression that mentions x[i] evaluates
 * with x free — expanding Table/Sum symbolically without capturing any stray
 * user definition. Restored by nm_heads_restore. */
typedef struct { const char* name; Rule* own; Rule* down; bool valid; } NmHeadSave;

static void nm_free_rule_chain(Rule* r) {
    while (r) {
        Rule* nx = r->next;
        expr_free(r->pattern);
        expr_free(r->replacement);
        free(r);
        r = nx;
    }
}

static void nm_heads_localize(NmHeadSave* sv, const char** names, size_t m) {
    for (size_t i = 0; i < m; i++) {
        SymbolDef* def = symtab_get_def(names[i]);
        sv[i].name  = names[i];
        sv[i].own   = def->own_values;
        sv[i].down  = def->down_values;
        sv[i].valid = true;
        def->own_values  = NULL;
        def->down_values = NULL;
    }
    if (m) eval_clock_bump();
}

static void nm_heads_restore(NmHeadSave* sv, size_t m) {
    for (size_t i = 0; i < m; i++) {
        if (!sv[i].valid) continue;
        SymbolDef* def = symtab_get_def(sv[i].name);
        nm_free_rule_chain(def->own_values);
        nm_free_rule_chain(def->down_values);
        def->own_values  = sv[i].own;
        def->down_values = sv[i].down;
        sv[i].valid = false;
    }
    if (m) eval_clock_bump();
}

/* ------------------------------------------------------------------ *
 *  Driver context bundle + point evaluation / comparison             *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr*        f_raw;          /* objective (borrowed)                     */
    Expr**       vars;           /* variable symbols (borrowed)             */
    size_t       n;
    FmVarBind*   binds;
    Expr**       g_exprs;        /* symbolic ∇f, or NULL → finite diff       */
    FmGenCon*    gens;
    size_t       ngens;
    FmBox*       boxes;
    const FmOpts* opts;
    const bool*  is_int;
    bool         any_int;
    const double* reg_lo;
    const double* reg_hi;
    /* Machine-precision fast path: a bytecode-compiled objective and per-general-
     * constraint programs over the effective variables, or NULL to use the
     * interpreter. Each point evaluation prefers the compiled program and falls
     * back to the interpreter when it is absent or reports a domain/non-finite
     * result — so the compiled path is a pure speedup, never a correctness risk. */
    CompiledProgram*  f_prog;    /* objective, or NULL                          */
    CompiledProgram** g_progs;   /* per-constraint (len ngens), entries may NULL */
    /* "PenaltyFunction" -> f. Borrowed / NULL ⇒ Automatic (built-in squared
     * penalty). When set, the global-search violation of each active constraint
     * is scored as f[violation] rather than violation^2; see nm_eval_pen. */
    Expr*             penalty_fn;
    /* Disjunctive (Or) constraints. Each contributes min-over-branch penalty via
     * the interpreter (fm_bool_penalty); len ndisj, borrowed from the setup. */
    FmDisjunction*    disj;
    size_t            ndisj;
} NmDriver;

/* Objective at the (already integer-rounded) point xr: the compiled program if
 * present and it returns a finite real, else the interpreter. */
static bool nm_eval_obj(NmDriver* D, const double* xr, double* out) {
    if (D->f_prog && compiled_eval_real(D->f_prog, xr, out) && isfinite(*out))
        return true;
    return fm_eval_scalar(D->f_raw, D->binds, xr, D->n, D->opts, out) && isfinite(*out);
}

/* Apply a user "PenaltyFunction" f to one nonnegative constraint-violation
 * magnitude m, returning f[m] as a double. Evaluator numeric diagnostics are
 * muted, as everywhere else in the trial-point loop. A non-numeric, non-finite,
 * or negative result is rejected so the caller can fall back to the built-in
 * m^2 for that term (a custom penalty must be a usable nonnegative score). */
static bool nm_apply_penalty_fn(Expr* pf, double m, double* out) {
    Expr* arg  = expr_new_real(m);
    Expr* call = expr_new_function(expr_copy(pf), &arg, 1);
    arith_warnings_mute_push();
    Expr* v = eval_and_free(call);
    arith_warnings_mute_pop();
    bool ok = v && fm_expr_to_double_real(v, out) && isfinite(*out) && *out >= 0.0;
    expr_free(v);
    return ok;
}

/* Penalty of a boolean-of-comparisons constraint tree at x, used for
 * disjunctive (Or) constraints:
 *   And[c...]        → Σ penalty(c)            (all must hold)
 *   Or[c...]         → min penalty(c)          (any one holding scores 0)
 *   Inequality[...]  → Σ over adjacent pairs   (a chained conjunction)
 *   binary compare   → squared violation: max(0,g)^2 for g<=0, h^2 for h==0,
 *                      or penalty_fn[violation] when a custom "PenaltyFunction"
 *                      is supplied (matching nm_eval_pen's per-term rule).
 * A satisfied leaf contributes 0, so the whole expression is 0 iff feasible and
 * the (total ≤ NM_FEAS_EPS) feasibility test carries over unchanged. The tree
 * shape is pre-validated by fm_bool_supported at collection time; this returns
 * false only if a leaf cannot be evaluated to a finite real at this point. */
static bool fm_bool_penalty(Expr* c, FmVarBind* binds, const double* x, size_t n,
                            const FmOpts* opts, Expr* penalty_fn, double* out) {
    const char* h = c->data.function.head->data.symbol.name;
    size_t ac = c->data.function.arg_count;
    if (h == SYM_And) {
        double total = 0.0;
        for (size_t i = 0; i < ac; i++) {
            double t;
            if (!fm_bool_penalty(c->data.function.args[i], binds, x, n, opts,
                                 penalty_fn, &t)) return false;
            total += t;
        }
        *out = total;
        return true;
    }
    if (h == SYM_Or) {
        double best = -1.0;
        for (size_t i = 0; i < ac; i++) {
            double t;
            if (!fm_bool_penalty(c->data.function.args[i], binds, x, n, opts,
                                 penalty_fn, &t)) return false;
            if (best < 0.0 || t < best) best = t;
        }
        *out = (best < 0.0) ? 0.0 : best;
        return true;
    }
    if (h == SYM_Inequality) {
        double total = 0.0;
        size_t npairs = (ac - 1) / 2;
        for (size_t k = 0; k < npairs; k++) {
            Expr* a  = c->data.function.args[2 * k];
            Expr* op = c->data.function.args[2 * k + 1];
            Expr* b  = c->data.function.args[2 * k + 2];
            Expr* pair_args[2] = { expr_copy(a), expr_copy(b) };
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol.name),
                                           pair_args, 2);
            double t;
            bool ok = fm_bool_penalty(pair, binds, x, n, opts, penalty_fn, &t);
            expr_free(pair);
            if (!ok) return false;
            total += t;
        }
        *out = total;
        return true;
    }
    /* single binary comparison → squared / custom violation */
    Expr* g = NULL; bool eq = false;
    if (!fm_constraint_to_g(c, &g, &eq)) return false;
    double d;
    bool ok = fm_eval_scalar(g, binds, x, n, opts, &d);
    expr_free(g);
    if (!ok || !isfinite(d)) return false;
    double m = eq ? fabs(d) : (d > 0.0 ? d : 0.0);
    if (m == 0.0) { *out = 0.0; return true; }
    double term;
    if (!penalty_fn || !nm_apply_penalty_fn(penalty_fn, m, &term)) term = m * m;
    *out = term;
    return true;
}

/* Σ pen(g_i) over the general constraints at xr, each constraint via its
 * compiled program if present, else the interpreter. The per-constraint term is
 * the built-in squared violation — max(0, g)^2 for an inequality g <= 0, h^2 for
 * an equality h == 0 — unless a "PenaltyFunction" f was supplied, in which case a
 * *violated* constraint contributes f[violation] instead (Automatic ≡ #^2 &, so
 * this generalises the default exactly). A satisfied inequality always
 * contributes 0, keeping the feasibility test (total ≤ NM_FEAS_EPS) intact.
 * Returns false if any constraint cannot be evaluated at all. */
static bool nm_eval_pen(NmDriver* D, const double* xr, double* out) {
    double total = 0.0;
    /* General (conjunctive) constraints: the compiled/penalty-fn path when either
     * is present, else the shared squared-penalty evaluator. */
    if (D->ngens > 0) {
        if (!D->g_progs && !D->penalty_fn) {
            double base;
            if (!fm_eval_penalty(D->gens, D->ngens, D->binds, xr, D->n, D->opts, &base))
                return false;
            total += base;
        } else {
            for (size_t k = 0; k < D->ngens; k++) {
                double d;
                bool got = D->g_progs && D->g_progs[k]
                        && compiled_eval_real(D->g_progs[k], xr, &d) && isfinite(d);
                if (!got && (!fm_eval_scalar(D->gens[k].expr, D->binds, xr, D->n,
                                             D->opts, &d) || !isfinite(d)))
                    return false;
                double m = D->gens[k].equality ? fabs(d) : (d > 0.0 ? d : 0.0);
                if (m == 0.0) continue;
                double term;
                if (!D->penalty_fn || !nm_apply_penalty_fn(D->penalty_fn, m, &term))
                    term = m * m;
                total += term;
            }
        }
    }
    /* Disjunctive (Or) constraints: each adds its minimum-branch penalty, which
     * is 0 exactly when at least one branch is satisfied. */
    for (size_t d = 0; d < D->ndisj; d++) {
        double dp;
        if (!fm_bool_penalty(D->disj[d].expr, D->binds, xr, D->n, D->opts,
                             D->penalty_fn, &dp))
            return false;
        total += dp;
    }
    *out = total;
    return true;
}

/* Objective value and total constraint violation at x. Integer coordinates
 * are rounded before evaluation. A non-evaluable objective or constraint is
 * treated as maximally bad so the search steers away from it. */
static void nm_eval(NmDriver* D, const double* x, double* f_out, double* pen_out) {
    size_t n = D->n;
    double* xr = (double*)malloc(sizeof(double) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) xr[i] = D->is_int[i] ? round(x[i]) : x[i];
    double fx;
    if (!nm_eval_obj(D, xr, &fx)) fx = 1e300;
    double pen = 0.0;
    if (!nm_eval_pen(D, xr, &pen)) pen = 1e300;
    free(xr);
    *f_out = fx;
    *pen_out = pen;
}

/* Deb's feasibility rules: is (fa, pa) a better candidate than (fb, pb)? */
static bool nm_better(double fa, double pa, double fb, double pb) {
    bool fa_feas = (pa <= NM_FEAS_EPS);
    bool fb_feas = (pb <= NM_FEAS_EPS);
    if (fa_feas && fb_feas) return fa < fb;
    if (fa_feas != fb_feas) return fa_feas;
    return pa < pb;
}

/* Clamp x into the search region and snap integer coordinates. */
static void nm_project(NmDriver* D, double* x) {
    for (size_t j = 0; j < D->n; j++) {
        if (x[j] < D->reg_lo[j]) x[j] = D->reg_lo[j];
        if (x[j] > D->reg_hi[j]) x[j] = D->reg_hi[j];
        if (D->is_int[j]) x[j] = round(x[j]);
    }
}

/* Penalized scalar objective used by NelderMead / SimulatedAnnealing. */
static double nm_phi(NmDriver* D, const double* x) {
    double f, p;
    nm_eval(D, x, &f, &p);
    return f + NM_PENALTY_MU * p;
}

/* ------------------------------------------------------------------ *
 *  Local polish of a candidate (shared by RandomSearch + the driver) *
 * ------------------------------------------------------------------ */

/* Mixed/integer coordinate descent: step each coordinate (±1, ±2 on integer
 * dims; scaled steps on continuous dims), accept any Deb-improvement. */
static void nm_int_descent(NmDriver* D, double* x, double* f_io, double* pen_io) {
    size_t n = D->n;
    double* t = (double*)malloc(sizeof(double) * n);
    bool improved = true;
    int iter = 0;
    while (improved && iter++ < 300) {
        improved = false;
        for (size_t j = 0; j < n; j++) {
            double steps[4];
            if (D->is_int[j]) {
                steps[0] = 1; steps[1] = -1; steps[2] = 2; steps[3] = -2;
            } else {
                double h = (D->reg_hi[j] - D->reg_lo[j]) * 0.05 + 1e-3;
                steps[0] = h; steps[1] = -h; steps[2] = 4 * h; steps[3] = -4 * h;
            }
            for (int s = 0; s < 4; s++) {
                for (size_t k = 0; k < n; k++) t[k] = x[k];
                t[j] = x[j] + steps[s];
                nm_project(D, t);
                double f, p;
                nm_eval(D, t, &f, &p);
                if (nm_better(f, p, *f_io, *pen_io)) {
                    for (size_t k = 0; k < n; k++) x[k] = t[k];
                    *f_io = f; *pen_io = p;
                    improved = true;
                }
            }
        }
    }
    free(t);
}

/* Deep-free a general-constraint array built by nm_polish_gens. */
static void fm_free_gens(FmGenCon* g, size_t ng, size_t n) {
    if (!g) return;
    for (size_t k = 0; k < ng; k++) {
        expr_free(g[k].expr);
        if (g[k].grad_exprs) {
            for (size_t i = 0; i < n; i++) expr_free(g[k].grad_exprs[i]);
            free(g[k].grad_exprs);
        }
    }
    free(g);
}

/* Flatten one disjunction branch (And / Inequality / single comparison — never
 * Or, since a top-level Or's args are the individual disjuncts) into owned
 * general constraints appended to *g, each with its symbolic gradient. Returns
 * false if any leaf is not a convertible comparison. */
static bool nm_collect_branch_gens(Expr* c, Expr** vars, size_t n,
                                   FmGenCon** g, size_t* ng, size_t* cap) {
    const char* h = c->data.function.head->data.symbol.name;
    size_t ac = c->data.function.arg_count;
    if (h == SYM_And) {
        for (size_t i = 0; i < ac; i++)
            if (!nm_collect_branch_gens(c->data.function.args[i], vars, n, g, ng, cap))
                return false;
        return true;
    }
    if (h == SYM_Inequality) {
        size_t npairs = (ac - 1) / 2;
        for (size_t k = 0; k < npairs; k++) {
            Expr* a  = c->data.function.args[2 * k];
            Expr* op = c->data.function.args[2 * k + 1];
            Expr* b  = c->data.function.args[2 * k + 2];
            Expr* pa[2] = { expr_copy(a), expr_copy(b) };
            Expr* pair = expr_new_function(expr_new_symbol(op->data.symbol.name), pa, 2);
            bool ok = nm_collect_branch_gens(pair, vars, n, g, ng, cap);
            expr_free(pair);
            if (!ok) return false;
        }
        return true;
    }
    Expr* ge = NULL; bool eq = false;
    if (!fm_constraint_to_g(c, &ge, &eq)) return false;
    if (*ng == *cap) {
        size_t nc = *cap ? (*cap) * 2 : 4;
        *g = (FmGenCon*)realloc(*g, sizeof(FmGenCon) * nc);
        *cap = nc;
    }
    (*g)[*ng].expr = ge;
    /* Finite-difference gradient (NULL ⇒ fm_run_penalty FDs on demand): the
     * optimisation variables carry transient value-bindings during the search,
     * so symbolic D[...] taken here differentiates a constant and yields 0. */
    (*g)[*ng].grad_exprs = NULL;
    (*g)[*ng].equality = eq;
    (*ng)++;
    return true;
}

/* Effective smooth-constraint set for a local polish at x when disjunctions are
 * present: deep copies of the conjunctive base D->gens, plus — for each
 * disjunction — the constraints of its currently-active (minimum-penalty) branch
 * at x. Folding in the active branch turns the non-smooth Or into a smooth local
 * problem the BFGS penalty solver can descend, so the polish refines *within* the
 * feasible region x already occupies instead of ignoring the Or and drifting into
 * the infeasible gap between branches (which left RandomSearch, whose only
 * descent is this polish, stranded). Returns NULL (and does not touch *nout) when
 * D->ndisj == 0, signalling the caller to use D->gens directly; otherwise returns
 * a newly-allocated array of length *nout, freed with fm_free_gens. */
static FmGenCon* nm_polish_gens(NmDriver* D, const double* x, size_t* nout) {
    if (D->ndisj == 0) return NULL;
    FmGenCon* g = NULL; size_t ng = 0, cap = 0;
    for (size_t k = 0; k < D->ngens; k++) {
        if (ng == cap) {
            size_t nc = cap ? cap * 2 : 4;
            g = (FmGenCon*)realloc(g, sizeof(FmGenCon) * nc);
            cap = nc;
        }
        g[ng].expr = expr_copy(D->gens[k].expr);
        g[ng].grad_exprs = NULL;   /* FD during search; see nm_collect_branch_gens */
        g[ng].equality = D->gens[k].equality;
        ng++;
    }
    for (size_t d = 0; d < D->ndisj; d++) {
        Expr* orx = D->disj[d].expr;          /* Or[branch, ...] */
        size_t ac = orx->data.function.arg_count;
        size_t best = 0; double bestp = -1.0;
        for (size_t i = 0; i < ac; i++) {
            double p;
            if (!fm_bool_penalty(orx->data.function.args[i], D->binds, x, D->n,
                                 D->opts, D->penalty_fn, &p))
                p = 1e300;
            if (bestp < 0.0 || p < bestp) { bestp = p; best = i; }
        }
        /* On an unconvertible branch, drop its augmentation (the global search's
         * exact penalty still gates feasibility); rewind any partial append. */
        size_t ng_save = ng;
        if (!nm_collect_branch_gens(orx->data.function.args[best], D->vars, D->n,
                                    &g, &ng, &cap)) {
            for (size_t k = ng_save; k < ng; k++) {
                expr_free(g[k].expr);
                if (g[k].grad_exprs) {
                    for (size_t i = 0; i < D->n; i++) expr_free(g[k].grad_exprs[i]);
                    free(g[k].grad_exprs);
                }
            }
            ng = ng_save;
        }
    }
    *nout = ng;
    return g;
}

/* Run the exact continuous local solver over all variables from x, confined
 * only by the real box constraints (D->boxes) — never the DE sampling region,
 * so a coordinate with no explicit bound is free to leave the default span
 * (the freedom the pure-continuous polish already has). When pin_int is true
 * each integer coordinate is pinned to its rounded value with a degenerate
 * [v, v] box (fm_line_search projects every trial into the box, so it stays
 * fixed); when false the integer coordinates are relaxed and solved as
 * continuous. x is overwritten with the refined point and (f, penalty) at the
 * integer-rounded point is returned; the caller decides whether to keep it. */
static void nm_continuous_solve(NmDriver* D, double* x, bool pin_int,
                                double* f_out, double* pen_out) {
    size_t n = D->n;
    FmBox* tb = (FmBox*)malloc(sizeof(FmBox) * n);
    if (!tb) { nm_eval(D, x, f_out, pen_out); return; }
    for (size_t j = 0; j < n; j++) {
        if (pin_int && D->is_int[j]) {
            double v = round(x[j]);
            x[j] = v;
            tb[j].has_lo = tb[j].has_hi = true;
            tb[j].lo = tb[j].hi = v;
        } else {
            tb[j] = D->boxes[j];             /* real constraint bounds, else free    */
        }
    }
    double fx = 0.0;
    bool saved_quiet = g_fm_quiet;
    g_fm_quiet = true;                        /* silence internal solver chatter      */
    size_t png = D->ngens;
    FmGenCon* pg = nm_polish_gens(D, x, &png);           /* NULL ⇒ use D->gens */
    FmGenCon* use = pg ? pg : D->gens;
    /* With a compiled objective, take the gradient by finite differences off it
     * (fm_grad_finite_diff → fm_eval_scalar → compiled) rather than evaluating
     * the symbolic gradient through the interpreter — 2n compiled evals beats
     * n symbolic evals each re-binding n OwnValues. */
    Expr** ge = D->f_prog ? NULL : D->g_exprs;
    if (png > 0)
        (void)fm_run_penalty(D->f_raw, D->vars, n, D->binds,
                             FM_METHOD_QUASINEWTON, ge, NULL, x,
                             use, png, tb, D->opts, &fx);
    else
        (void)fm_run_bfgs(D->f_raw, D->vars, n, D->binds, ge, x,
                          NULL, 0, 0.0, tb, D->opts, &fx);
    if (pg) fm_free_gens(pg, png, n);
    g_fm_quiet = saved_quiet;
    free(tb);
    nm_eval(D, x, f_out, pen_out);
}

static void nm_local_polish(NmDriver* D, double* x, double* f_out, double* pen_out) {
    if (D->any_int) {
        nm_project(D, x);
        nm_eval(D, x, f_out, pen_out);
        nm_int_descent(D, x, f_out, pen_out);
        /* A mixed-integer problem whose feasible region lies outside the DE
         * sampling region (e.g. the pressure-vessel MINLP, feasible near
         * x3 ~ 52 with the default +-10 span) is invisible to the region-bound
         * integer descent above. Recover it with the continuous-relaxation +
         * rounding heuristic: solve the continuous relaxation (integers relaxed,
         * every coordinate free of the sampling region) to locate the basin,
         * round the integer coordinates, then refine the continuous coordinates
         * with the integers pinned. Adopt the result only when it is a
         * Deb-improvement, so this can never worsen the region-descent answer.
         * Skipped when every variable is integer (no continuous coordinate to
         * free), where the region descent is already the whole story. */
        bool has_cont = false;
        for (size_t j = 0; j < D->n; j++) if (!D->is_int[j]) { has_cont = true; break; }
        if (has_cont) {
            double* xr = (double*)malloc(sizeof(double) * D->n);
            if (xr) {
                for (size_t j = 0; j < D->n; j++) xr[j] = x[j];
                double fr, pr;
                nm_continuous_solve(D, xr, false, &fr, &pr);     /* relaxation      */
                for (size_t j = 0; j < D->n; j++)
                    if (D->is_int[j]) xr[j] = round(xr[j]);
                nm_continuous_solve(D, xr, true, &fr, &pr);      /* pin + refine    */
                if (nm_better(fr, pr, *f_out, *pen_out)) {
                    for (size_t j = 0; j < D->n; j++) x[j] = xr[j];
                    *f_out = fr; *pen_out = pr;
                }
                free(xr);
            }
            nm_int_descent(D, x, f_out, pen_out);   /* re-settle integers in-region */
        }
    } else {
        double fx = 0.0;
        bool saved_quiet = g_fm_quiet;
        g_fm_quiet = true;                 /* silence internal solver chatter */
        size_t png = D->ngens;
        FmGenCon* pg = nm_polish_gens(D, x, &png);       /* NULL ⇒ use D->gens */
        FmGenCon* use = pg ? pg : D->gens;
        /* Compiled objective ⇒ finite-difference gradient off it (see the
         * matching note in nm_continuous_solve). */
        Expr** ge = D->f_prog ? NULL : D->g_exprs;
        if (png > 0)
            (void)fm_run_penalty(D->f_raw, D->vars, D->n, D->binds,
                                 FM_METHOD_QUASINEWTON, ge, NULL, x,
                                 use, png, D->boxes, D->opts, &fx);
        else
            (void)fm_run_bfgs(D->f_raw, D->vars, D->n, D->binds, ge, x,
                              NULL, 0, 0.0, D->boxes, D->opts, &fx);
        if (pg) fm_free_gens(pg, png, D->n);
        g_fm_quiet = saved_quiet;
        nm_eval(D, x, f_out, pen_out);
    }
}

/* ------------------------------------------------------------------ *
 *  Global engines                                                     *
 * ------------------------------------------------------------------ */

/* Seed the leading members of a DE population from a user "InitialPoints" list
 * and return how many were seeded (the rest are filled at random by the caller).
 * Each seed must be a length-n list of reals (evaluated + numericalized, so Pi
 * etc. are fine); it is projected into the search region and its integer
 * coordinates are rounded. Unlike nm_simplex_from_points this is per-member
 * rather than all-or-nothing: a malformed point ends the seeding and the rest of
 * the population is random, so one bad point is not fatal. Extra points beyond
 * the population size are ignored. Returns 0 (→ fully random) when the list is
 * absent or empty. */
static size_t nm_population_from_points(NmDriver* D, const Expr* pts,
                                        size_t n, size_t NP, double* pop) {
    if (!pts || !nm_is_head(pts, SYM_List) || pts->data.function.arg_count == 0)
        return 0;
    size_t npts = pts->data.function.arg_count;
    size_t use  = npts < NP ? npts : NP;
    size_t seeded = 0;
    for (size_t i = 0; i < use; i++) {
        Expr* p = pts->data.function.args[i];
        if (!nm_is_head(p, SYM_List) || p->data.function.arg_count != n) break;
        double* row = &pop[i * n];
        bool ok = true;
        for (size_t j = 0; j < n; j++) {
            Expr* c  = eval_and_free(expr_copy(p->data.function.args[j]));
            Expr* cn = c ? numericalize(c, numeric_machine_spec()) : NULL;
            double v = 0.0;
            ok = cn && fm_expr_to_double_real(cn, &v);
            expr_free(c); expr_free(cn);
            if (!ok) break;
            row[j] = v;
        }
        if (!ok) break;
        nm_project(D, row);
        for (size_t j = 0; j < n; j++)
            if (D->is_int[j]) row[j] = round(row[j]);
        seeded++;
    }
    return seeded;
}

/* DifferentialEvolution, DE/rand/1/bin, with Deb-rule selection. */
static void nm_de(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* An explicit "SearchPoints" is honored verbatim (only floored at the DE
     * minimum of 4 members its DE/rand/1 mutation needs). The automatic
     * population is Storn & Price's 10n — the textbook default — clamped to
     * [15, 200]: the floor keeps a tiny problem's population workable, and the
     * ceiling 200 = 10·20 bounds the per-generation cost on high-dimensional
     * problems. This sizing is the same whether the method is an explicit
     * "DifferentialEvolution" or Method -> Automatic (there is no longer a lower
     * historical clamp for the explicit form). NelderMead restarts and
     * RandomSearch starts honor SearchPoints verbatim. */
    size_t NP;
    if (nc->search_points > 0) {
        NP = (size_t)nc->search_points;
        if (NP < 4) NP = 4;
    } else {
        NP = 10 * n;
        if (NP < 15)  NP = 15;
        if (NP > 200) NP = 200;
    }
    double F  = nc->F  > 0.0 ? nc->F  : 0.6;
    double CR = nc->CR >= 0.0 ? nc->CR : 0.9;
    int64_t maxgen = D->opts->max_iter > 0 ? D->opts->max_iter : 100;
    /* Method -> Automatic with no explicit MaxIterations scales the generation
     * budget with dimension (150n). A deceptive multimodal landscape such as
     * Michalewicz-10 (Sin[...]^20 ridges, near-flat elsewhere, so the post-polish
     * cannot rescue a weak basin) needs far more than the flat 100 generations to
     * find a good basin. This is free on easy problems: the convergence
     * early-break below stops as soon as the population collapses. An explicit
     * "DifferentialEvolution" (or a user MaxIterations) keeps the flat budget. */
    if (nc->method == NM_AUTO && !D->opts->max_iter_set)
        maxgen = 150 * (int64_t)n;

    double* pop   = (double*)malloc(sizeof(double) * NP * n);
    double* fpop  = (double*)malloc(sizeof(double) * NP);
    double* ppop  = (double*)malloc(sizeof(double) * NP);
    double* trial = (double*)malloc(sizeof(double) * n);

    /* "InitialPoints" seeds the leading population members; the remainder are
     * random. When no list is supplied seeded == 0 and every member is random —
     * the RNG stream is then identical to before, so seeded runs reproduce. */
    size_t seeded = nm_population_from_points(D, nc->init_points, n, NP, pop);
    for (size_t p = 0; p < NP; p++) {
        if (p >= seeded) {
            for (size_t j = 0; j < n; j++) {
                double v = nm_rng_range(rng, rlo[j], rhi[j]);
                if (D->is_int[j]) v = round(v);
                pop[p * n + j] = v;
            }
        }
        nm_eval(D, &pop[p * n], &fpop[p], &ppop[p]);
    }
    size_t bi = 0;
    for (size_t p = 1; p < NP; p++)
        if (nm_better(fpop[p], ppop[p], fpop[bi], ppop[bi])) bi = p;
    for (size_t j = 0; j < n; j++) xbest[j] = pop[bi * n + j];
    *fbest = fpop[bi];
    *penbest = ppop[bi];

    for (int64_t g = 0; g < maxgen; g++) {
        for (size_t p = 0; p < NP; p++) {
            size_t r1, r2, r3;
            if (NP < 4) break;
            do { r1 = nm_rng_next(rng) % NP; } while (r1 == p);
            do { r2 = nm_rng_next(rng) % NP; } while (r2 == p || r2 == r1);
            do { r3 = nm_rng_next(rng) % NP; } while (r3 == p || r3 == r1 || r3 == r2);
            size_t jr = nm_rng_next(rng) % n;
            for (size_t j = 0; j < n; j++) {
                if (nm_rng_unif(rng) < CR || j == jr) {
                    double base = pop[r1 * n + j];
                    double v = base + F * (pop[r2 * n + j] - pop[r3 * n + j]);
                    /* Bounce-back on a bound violation instead of clamping to
                     * the bound. Clamping strands the search: once several
                     * members share the exact boundary value for coordinate j,
                     * their pairwise mutation differentials for j collapse to
                     * zero, no mutant can lift the coordinate off the wall, and
                     * it freezes there for the rest of the run — the Schwefel
                     * failure, where half the coordinates pinned to ±500 while
                     * the rest found the true basin. Bounce-back (Price, Storn &
                     * Lampinen 2005) reflects the overshoot to a random point
                     * between the base vector and the violated bound: feasible,
                     * but almost never exactly on the boundary, so the
                     * population keeps its diversity. */
                    if (v < rlo[j])      v = base + nm_rng_unif(rng) * (rlo[j] - base);
                    else if (v > rhi[j]) v = base + nm_rng_unif(rng) * (rhi[j] - base);
                    if (D->is_int[j]) v = round(v);
                    trial[j] = v;
                } else {
                    trial[j] = pop[p * n + j];
                }
            }
            double ft, pt;
            nm_eval(D, trial, &ft, &pt);
            if (nm_better(ft, pt, fpop[p], ppop[p])) {
                for (size_t j = 0; j < n; j++) pop[p * n + j] = trial[j];
                fpop[p] = ft; ppop[p] = pt;
                if (nm_better(ft, pt, *fbest, *penbest)) {
                    for (size_t j = 0; j < n; j++) xbest[j] = trial[j];
                    *fbest = ft; *penbest = pt;
                }
            }
        }
        /* Converged if the best is feasible and the feasible sub-population's
         * objective spread has collapsed to the requested tolerance. */
        if (*penbest <= NM_FEAS_EPS) {
            double fmin = 1e300, fmax = -1e300;
            size_t cnt = 0;
            for (size_t p = 0; p < NP; p++) {
                if (ppop[p] <= NM_FEAS_EPS) {
                    if (fpop[p] < fmin) fmin = fpop[p];
                    if (fpop[p] > fmax) fmax = fpop[p];
                    cnt++;
                }
            }
            /* Convergence tolerance on the feasible sub-population's objective
             * spread. An explicit "Tolerance" sub-option sets it directly;
             * otherwise it is derived from AccuracyGoal / PrecisionGoal. */
            double tol = nc->tolerance > 0.0
                       ? nc->tolerance
                       : pow(10.0, -D->opts->acc_goal_digits)
                         + (1.0 + fabs(*fbest)) * pow(10.0, -D->opts->prec_goal_digits);
            if (cnt >= NP / 2 && (fmax - fmin) <= tol) break;
        }
    }

    /* Multi-start local polish over the final population. Polishing only the
     * single global best strands DE in whichever basin that one member happened
     * to occupy, so a larger population or a shorter run — both of which leave
     * the population less converged — can report a WORSE optimum (Griewank-10:
     * SearchPoints -> 100 gave 0.197 where the default 40 gave 0.044). This is
     * the same non-monotonicity the SimulatedAnnealing per-chain polish removed.
     * Polish the best Min[2 n, 50] *distinct* members — skipping any that sits in
     * an already-polished basin — and keep the deepest local minimum, ranking
     * basins by their minima rather than by the raw fitness of the member that
     * found them. Gated by "PostProcess" -> False (then the raw global best is
     * carried through as before) and confined to continuous problems, where the
     * BFGS polish is cheap; a mixed-integer run keeps its single driver polish so
     * its integer-descent cost is unchanged. */
    if (nc->post_process != 0 && !D->any_int && NP >= 1) {
        int64_t polish_cap = 2 * (int64_t)n < 50 ? 2 * (int64_t)n : 50;
        if (polish_cap > (int64_t)NP) polish_cap = (int64_t)NP;
        if (polish_cap < 1) polish_cap = 1;

        /* Two members belong to the same basin for dedup purposes when they agree
         * to within 1e-3 of the search span on every coordinate (floored, so a
         * degenerate zero-width span still separates distinct points). */
        double* btol = (double*)malloc(sizeof(double) * n);
        for (size_t j = 0; j < n; j++) {
            double t = 1e-3 * (rhi[j] - rlo[j]);
            btol[j] = t > 1e-6 ? t : 1e-6;
        }
        unsigned char* used = (unsigned char*)calloc(NP, 1);
        double* seeds = (double*)malloc(sizeof(double) * (size_t)polish_cap * n);
        double* xp = (double*)malloc(sizeof(double) * n);
        size_t nseeds = 0;
        int64_t done = 0;

        while (used && seeds && xp && btol && done < polish_cap) {
            /* best not-yet-considered member by Deb's rules */
            size_t best = NP;
            for (size_t p = 0; p < NP; p++) {
                if (used[p]) continue;
                if (best == NP || nm_better(fpop[p], ppop[p], fpop[best], ppop[best]))
                    best = p;
            }
            if (best == NP) break;
            used[best] = 1;

            /* Skip a member already represented by a polished basin seed. */
            bool dup = false;
            for (size_t s = 0; s < nseeds && !dup; s++) {
                bool same = true;
                for (size_t j = 0; j < n; j++)
                    if (fabs(pop[best * n + j] - seeds[s * n + j]) > btol[j]) {
                        same = false; break;
                    }
                if (same) dup = true;
            }
            if (dup) continue;   /* does not count against polish_cap */
            for (size_t j = 0; j < n; j++) seeds[nseeds * n + j] = pop[best * n + j];
            nseeds++;

            for (size_t j = 0; j < n; j++) xp[j] = pop[best * n + j];
            double fp = fpop[best], pp = ppop[best];
            nm_local_polish(D, xp, &fp, &pp);
            if (nm_better(fpop[best], ppop[best], fp, pp)) {  /* overshoot guard */
                for (size_t j = 0; j < n; j++) xp[j] = pop[best * n + j];
                fp = fpop[best]; pp = ppop[best];
            }
            if (nm_better(fp, pp, *fbest, *penbest)) {
                for (size_t j = 0; j < n; j++) xbest[j] = xp[j];
                *fbest = fp; *penbest = pp;
            }
            done++;
        }
        free(btol); free(used); free(seeds); free(xp);
    }

    free(pop); free(fpop); free(ppop); free(trial);
}

/* Build the (n+1)-vertex initial simplex from a user "InitialPoints" list. Each
 * element must be a length-n list of reals (evaluated + numericalized, so Pi
 * etc. are fine). If fewer than n+1 points are given, the remaining vertices are
 * seeded by perturbing the first point along successive axes; extra points are
 * ignored. Returns false (→ fall back to a random simplex) if the list is
 * malformed or any used point has the wrong dimension or a non-numeric entry. */
static bool nm_simplex_from_points(NmDriver* D, const Expr* pts, size_t n, double* V) {
    if (!pts || !nm_is_head(pts, SYM_List) || pts->data.function.arg_count == 0)
        return false;
    size_t npts = pts->data.function.arg_count;
    size_t use  = npts < (n + 1) ? npts : (n + 1);
    for (size_t i = 0; i < use; i++) {
        Expr* p = pts->data.function.args[i];
        if (!nm_is_head(p, SYM_List) || p->data.function.arg_count != n) return false;
        for (size_t j = 0; j < n; j++) {
            Expr* c  = eval_and_free(expr_copy(p->data.function.args[j]));
            Expr* cn = c ? numericalize(c, numeric_machine_spec()) : NULL;
            double v;
            bool ok = cn && fm_expr_to_double_real(cn, &v);
            expr_free(c); expr_free(cn);
            if (!ok) return false;
            V[i * n + j] = v;
        }
    }
    for (size_t i = use; i <= n; i++) {
        for (size_t j = 0; j < n; j++) V[i * n + j] = V[j];
        size_t d = i - 1;
        double step = (D->reg_hi[d] - D->reg_lo[d]) * 0.1;
        if (step == 0.0) step = 1.0;
        V[i * n + d] += step;
    }
    for (size_t i = 0; i <= n; i++) nm_project(D, &V[i * n]);
    return true;
}

/* NelderMead downhill simplex on the penalized objective, with restarts. Each
 * restart's converged best vertex is polished into its basin minimum and the
 * restarts are ranked by those minima (default Min[2 n, 20] restarts), so a
 * multimodal surface improves with more restarts instead of being decided by the
 * raw simplex vertices. Gated by "PostProcess" -> False. */
static void nm_neldermead(NmDriver* D, const NmConfig* nc, NmRng* rng,
                          double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* Explicit "SearchPoints" is honored verbatim (it was silently clamped at 20
     * before, contradicting the documented "honored verbatim"). The automatic
     * default runs Min[2 n, 20] random-simplex restarts rather than a flat 4, so
     * a multimodal surface is probed from several basins — the same "sample many
     * independent starts, keep the deepest" principle the other engines use. */
    int restarts;
    if (nc->search_points > 0) restarts = nc->search_points;
    else restarts = n > 1 ? (2 * (int)n < 20 ? 2 * (int)n : 20) : 2;
    if (restarts < 1) restarts = 1;
    /* Simplex coefficients: reflection (ReflectRatio, default 1), expansion
     * (ExpandRatio, default 2), contraction toward the centroid (ContractRatio,
     * default 0.5), shrink toward the best vertex (ShrinkRatio, default 0.5).
     * Tolerance is the objective-spread convergence threshold. */
    double rr  = nc->reflect_ratio  > 0.0 ? nc->reflect_ratio  : 1.0;
    double er  = nc->expand_ratio   > 0.0 ? nc->expand_ratio   : 2.0;
    double cr  = nc->contract_ratio > 0.0 ? nc->contract_ratio : 0.5;
    double sr  = nc->shrink_ratio   > 0.0 ? nc->shrink_ratio   : 0.5;
    double tol = nc->tolerance      > 0.0 ? nc->tolerance      : 1e-12;
    int64_t maxit = D->opts->max_iter > 0 ? D->opts->max_iter * (int64_t)(5 * n)
                                          : 200 * (int64_t)n;
    if (maxit < 100) maxit = 100;

    /* Domain-convergence scale: the simplex must shrink to a small fraction of
     * the search region — not just reach a flat objective — before we declare
     * convergence. Without this, a broad plateau (the flat outer region of the
     * Easom function, where f ≈ 0 everywhere away from a narrow spike) trips the
     * objective-spread test on the first iteration and the simplex never moves. */
    double rscale = 1.0;
    for (size_t j = 0; j < n; j++) { double e = rhi[j] - rlo[j]; if (e > rscale) rscale = e; }
    double xdtol = 1e-6 * rscale;

    double* V  = (double*)malloc(sizeof(double) * (n + 1) * n);
    double* fv = (double*)malloc(sizeof(double) * (n + 1));
    double* xc = (double*)malloc(sizeof(double) * n);
    double* xr = (double*)malloc(sizeof(double) * n);
    double* xe = (double*)malloc(sizeof(double) * n);
    bool have = false;

    for (int rs = 0; rs < restarts; rs++) {
        /* Restart 0 uses the user's "InitialPoints" simplex when supplied and
         * valid; every other restart (and the fallback) is a random simplex. */
        bool seeded = (rs == 0) && nm_simplex_from_points(D, nc->init_points, n, V);
        if (!seeded) {
            for (size_t j = 0; j < n; j++) V[j] = nm_rng_range(rng, rlo[j], rhi[j]);
            nm_project(D, &V[0]);
            for (size_t i = 1; i <= n; i++) {
                for (size_t j = 0; j < n; j++) V[i * n + j] = V[j];
                size_t d = i - 1;
                double step = (rhi[d] - rlo[d]) * 0.1;
                if (step == 0.0) step = 1.0;
                V[i * n + d] += step;
                nm_project(D, &V[i * n]);
            }
        }
        for (size_t i = 0; i <= n; i++) fv[i] = nm_phi(D, &V[i * n]);

        for (int64_t it = 0; it < maxit; it++) {
            size_t lo = 0, hi = 0, nh = 0;
            for (size_t i = 1; i <= n; i++) {
                if (fv[i] < fv[lo]) lo = i;
                if (fv[i] > fv[hi]) hi = i;
            }
            nh = (hi == 0) ? 1 : 0;
            for (size_t i = 0; i <= n; i++)
                if (i != hi && fv[i] > fv[nh]) nh = i;
            double xspread = 0.0;
            for (size_t i = 0; i <= n; i++)
                for (size_t j = 0; j < n; j++) {
                    double dd = fabs(V[i * n + j] - V[lo * n + j]);
                    if (dd > xspread) xspread = dd;
                }
            if (fabs(fv[hi] - fv[lo]) <= tol * (1.0 + fabs(fv[lo])) && xspread <= xdtol)
                break;

            for (size_t j = 0; j < n; j++) xc[j] = 0.0;
            for (size_t i = 0; i <= n; i++)
                if (i != hi)
                    for (size_t j = 0; j < n; j++) xc[j] += V[i * n + j];
            for (size_t j = 0; j < n; j++) xc[j] /= (double)n;

            for (size_t j = 0; j < n; j++) xr[j] = xc[j] + rr * (xc[j] - V[hi * n + j]);
            nm_project(D, xr);
            double frr = nm_phi(D, xr);
            if (frr < fv[lo]) {
                for (size_t j = 0; j < n; j++)
                    xe[j] = xc[j] + er * (xc[j] - V[hi * n + j]);
                nm_project(D, xe);
                double fe = nm_phi(D, xe);
                if (fe < frr) { for (size_t j = 0; j < n; j++) V[hi * n + j] = xe[j]; fv[hi] = fe; }
                else          { for (size_t j = 0; j < n; j++) V[hi * n + j] = xr[j]; fv[hi] = frr; }
            } else if (frr < fv[nh]) {
                for (size_t j = 0; j < n; j++) V[hi * n + j] = xr[j];
                fv[hi] = frr;
            } else {
                for (size_t j = 0; j < n; j++)
                    xe[j] = xc[j] + cr * (V[hi * n + j] - xc[j]);
                nm_project(D, xe);
                double fc = nm_phi(D, xe);
                if (fc < fv[hi]) { for (size_t j = 0; j < n; j++) V[hi * n + j] = xe[j]; fv[hi] = fc; }
                else {
                    for (size_t i = 0; i <= n; i++) {
                        if (i == lo) continue;
                        for (size_t j = 0; j < n; j++)
                            V[i * n + j] = V[lo * n + j] + sr * (V[i * n + j] - V[lo * n + j]);
                        nm_project(D, &V[i * n]);
                        fv[i] = nm_phi(D, &V[i * n]);
                    }
                }
            }
        }
        size_t lo = 0;
        for (size_t i = 1; i <= n; i++) if (fv[i] < fv[lo]) lo = i;
        double f, p;
        for (size_t j = 0; j < n; j++) xr[j] = V[lo * n + j];  /* xr: polish buffer */
        nm_eval(D, xr, &f, &p);
        /* Polish this restart's best vertex into its basin minimum before ranking
         * restarts, rather than ranking them by their raw simplex vertices and
         * polishing only the single winner afterward. The lowest converged vertex
         * across restarts need not sit in the deepest basin, so ranking by local
         * minima is what lets more restarts help rather than hurt. Skipped under
         * "PostProcess" -> False (then the raw vertex is ranked, as before); a
         * BFGS overshoot falls back to the raw vertex. */
        if (nc->post_process != 0) {
            double fr = f, pr = p;
            nm_local_polish(D, xr, &fr, &pr);
            if (nm_better(f, p, fr, pr)) {
                for (size_t j = 0; j < n; j++) xr[j] = V[lo * n + j];
            } else { f = fr; p = pr; }
        }
        if (!have || nm_better(f, p, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = xr[j];
            *fbest = f; *penbest = p; have = true;
        }
    }
    free(V); free(fv); free(xc); free(xr); free(xe);
}

/* RandomSearch: multiple random starts, each refined by the local solver, best
 * local minimum kept. This is already the "polish each start, rank by basin
 * depth" pattern; the two other multi-start engines were brought in line with it.
 * Pure multi-start local search has no global move, so on a search box far wider
 * than the optimum's basin (e.g. Griewank over [-600, 600], where the central
 * bowl is ~1e-11 of the volume in 10-D) no attainable number of random starts
 * reaches it — DifferentialEvolution / SimulatedAnnealing are the engines for
 * that shape. */
static void nm_randomsearch(NmDriver* D, const NmConfig* nc, NmRng* rng,
                            double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    /* Explicit "SearchPoints" is honored verbatim — it was silently capped at 40,
     * so SearchPoints -> 1000 was a no-op that returned the 40-start result. The
     * automatic default keeps a bound (runtime is one local solve per start). */
    int K;
    if (nc->search_points > 0) {
        K = nc->search_points;          /* honored verbatim (was capped at 40) */
        if (K < 1) K = 1;
    } else {
        K = n > 1 ? (int)(8 * n) : 12;  /* automatic default, bounded for runtime */
        if (K < 4)  K = 4;
        if (K > 40) K = 40;
    }
    double* x = (double*)malloc(sizeof(double) * n);
    bool have = false;
    for (int k = 0; k < K; k++) {
        for (size_t j = 0; j < n; j++) x[j] = nm_rng_range(rng, rlo[j], rhi[j]);
        nm_project(D, x);
        double f, p;
        nm_local_polish(D, x, &f, &p);
        if (!have || nm_better(f, p, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = x[j];
            *fbest = f; *penbest = p; have = true;
        }
    }
    free(x);
}

/* SimulatedAnnealing "BoltzmannExponent" -> f. Evaluate f[i, df, f0] to the
 * real exponent whose Exp is the Metropolis acceptance probability for an
 * uphill move: i is the (1-based) iteration, df ≥ 0 the objective increase, f0
 * the current objective. Evaluator numeric diagnostics are muted, as everywhere
 * in the trial-point loop. Returns false — so the caller falls back to the
 * built-in geometric-cooling exponent — if f does not yield a finite real. */
static bool nm_boltzmann_exponent(Expr* bf, int64_t i, double df, double f0,
                                  double* out) {
    Expr* a[3];
    a[0] = expr_new_integer(i);
    a[1] = expr_new_real(df);
    a[2] = expr_new_real(f0);
    Expr* call = expr_new_function(expr_copy(bf), a, 3);
    arith_warnings_mute_push();
    Expr* v = eval_and_free(call);
    arith_warnings_mute_pop();
    bool ok = v && fm_expr_to_double_real(v, out) && isfinite(*out);
    expr_free(v);
    return ok;
}

/* SimulatedAnnealing with geometric cooling; tracks the best point seen.
 *
 * Honors three "SimulatedAnnealing" sub-options:
 *   "SearchPoints" -> K       run K independent annealing chains from random
 *                             starts and keep the best local minimum (default
 *                             Automatic = Min[2 n, 50], following Mathematica);
 *   "PerturbationScale" -> s  multiply the trial-step size by s (default 1.0);
 *   "BoltzmannExponent" -> f  use Exp[f[i, df, f0]] as the acceptance
 *                             probability for an uphill move (default -df/T).
 *
 * A rugged, many-basin landscape (Griewank, Rastrigin, ...) is not solved by a
 * single annealing walk: which basin one walk lands in is close to luck. So the
 * default runs Min[2 n, 50] independent chains and — crucially — polishes the
 * best raw point of *each* chain into its basin minimum before ranking them,
 * exactly as RandomSearch does per restart. Ranking chains by their local
 * minima rather than by their random-walk lows is what makes the reported
 * optimum improve, not degrade, as SearchPoints grows: the lowest point a walk
 * happens to visit is often in a shallower basin than a slightly higher point
 * that sits above a deeper one. The RNG is still drawn in the original order
 * within each chain, so a seeded single-chain run (SearchPoints -> 1) anneals
 * identically to before and reaches the same polished point the driver's
 * post-process already produced. */
static void nm_sa(NmDriver* D, const NmConfig* nc, NmRng* rng,
                  double* xbest, double* fbest, double* penbest) {
    size_t n = D->n;
    const double* rlo = D->reg_lo;
    const double* rhi = D->reg_hi;
    double pscale = nc->perturb_scale > 0.0 ? nc->perturb_scale : 1.0;
    Expr*  bf     = nc->boltzmann_fn;

    /* "SearchPoints" -> K independent annealing chains. Default Automatic =
     * Min[Max[2 n, 12], 50]. The 2 n scaling follows Mathematica, but a rugged,
     * many-basin landscape in low dimension (Eggholder, Schwefel, Griewank, ...)
     * has far more basins than 2 n = 4 starts can cover, so which basin the run
     * reports is left to luck. A floor of 12 independent starts is what turns
     * these from "sometimes finds the global" into "reliably finds it" at a cost
     * of a few hundredths of a second — SA is opt-in (the automatic method is
     * DifferentialEvolution), so only a caller who asked for it pays. */
    int64_t K;
    if (nc->search_points > 0) {
        K = (int64_t)nc->search_points;
    } else {
        K = 2 * (int64_t)n;
        if (K < 12) K = 12;
        if (K > 50) K = 50;
    }
    if (K < 1) K = 1;

    /* Per-chain iteration budget = "LevelIterations" trial moves at each of the
     * MaxIterations temperature levels, so per_chain = MaxIterations *
     * LevelIterations — Mathematica's semantics, where LevelIterations is the
     * dwell at each level. The old hard-coded 50 was exactly that implicit
     * default. An explicit "LevelIterations" is honored verbatim (like
     * "SearchPoints"): the caller asked for that budget, so the automatic
     * runtime caps below are skipped. The automatic budget keeps runtime bounded
     * — a single chain keeps the original schedule, and many search points share
     * a bounded aggregate, with a floor so each chain still anneals. */
    int64_t level_it  = nc->level_iterations > 0 ? (int64_t)nc->level_iterations : 50;
    int64_t base_iter = D->opts->max_iter > 0 ? D->opts->max_iter : 100;
    int64_t per_chain = base_iter * level_it;
    if (nc->level_iterations <= 0) {
        if (per_chain > 20000) per_chain = 20000;
        if (K > 1 && K * per_chain > NM_SA_TOTAL_CAP) {
            per_chain = NM_SA_TOTAL_CAP / K;
            if (per_chain < 300) per_chain = 300;
        }
    }

    double* x  = (double*)malloc(sizeof(double) * n);
    double* xn = (double*)malloc(sizeof(double) * n);
    double* xc = (double*)malloc(sizeof(double) * n);  /* best raw point in chain */
    bool have = false;

    for (int64_t chain = 0; chain < K; chain++) {
        for (size_t j = 0; j < n; j++) x[j] = nm_rng_range(rng, rlo[j], rhi[j]);
        nm_project(D, x);
        double fx, px;
        nm_eval(D, x, &fx, &px);
        double phi = fx + NM_PENALTY_MU * px;
        /* Track this chain's best raw point separately, so it can be polished on
         * its own and compared against the other chains' basin minima. */
        for (size_t j = 0; j < n; j++) xc[j] = x[j];
        double fc = fx, pc = px;

        /* Adaptive temperature scale. The acceptance exponent -d/T only anneals
         * when T sits on the scale of the objective differences d. With a fixed
         * T in [1e-4, 1] and an objective ranging over hundreds (Eggholder,
         * Schwefel, ...), every uphill move is rejected and the "anneal"
         * degenerates into greedy descent from the start point — so only many
         * random restarts, never the walk itself, ever crossed a ridge into a
         * deeper basin. Probe a handful of full-temperature trial steps to
         * measure the typical |Δφ| and use it as the scale, so a move that
         * worsens φ by one characteristic step is accepted with probability
         * e^{-1/T}. The probe draws from its own RNG so the annealing walk's
         * stream is untouched — a seeded chain takes exactly the trajectory it
         * did before the scale was introduced. xn is scratch, overwritten by the
         * main loop's first proposal. */
        double scale = 0.0;
        {
            NmRng prng;
            nm_rng_seed(&prng, nc->seed ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(chain + 1)));
            int nb = 0;
            for (int b = 0; b < NM_SA_BURN_IN; b++) {
                for (size_t j = 0; j < n; j++) {
                    double span = rhi[j] - rlo[j];
                    xn[j] = x[j] + pscale * (span * 0.1 * 1.1 * nm_rng_normal(&prng));
                }
                nm_project(D, xn);
                double fb, pb;
                nm_eval(D, xn, &fb, &pb);
                double db = fabs((fb + NM_PENALTY_MU * pb) - phi);
                if (isfinite(db)) { scale += db; nb++; }
            }
            scale = nb > 0 ? scale / nb : 0.0;
            if (!(scale > 0.0) || !isfinite(scale))
                scale = fabs(phi) > 1.0 ? fabs(phi) : 1.0;
        }

        double T = 1.0;
        for (int64_t it = 0; it < per_chain; it++) {
            for (size_t j = 0; j < n; j++) {
                double span = rhi[j] - rlo[j];
                xn[j] = x[j]
                      + pscale * (span * 0.1 * (0.1 + T) * nm_rng_normal(rng));
            }
            nm_project(D, xn);
            double fn2, pn2;
            nm_eval(D, xn, &fn2, &pn2);
            double phin = fn2 + NM_PENALTY_MU * pn2;
            double d = phin - phi;
            bool accept;
            if (d < 0.0) {
                accept = true;
            } else {
                double expo;
                if (!bf || !nm_boltzmann_exponent(bf, it + 1, d, phi, &expo))
                    expo = -d / (T * scale + 1e-12);
                accept = nm_rng_unif(rng) < exp(expo);   /* NaN prob ⇒ reject */
            }
            if (accept) {
                for (size_t j = 0; j < n; j++) x[j] = xn[j];
                phi = phin; fx = fn2; px = pn2;
                if (nm_better(fx, px, fc, pc)) {
                    for (size_t j = 0; j < n; j++) xc[j] = x[j];
                    fc = fx; pc = px;
                }
            }
            T *= 0.995;
            if (T < 1e-4) T = 1e-4;
        }

        /* Polish this chain's best raw point into its basin minimum, then keep
         * the global best of the polished candidates. A BFGS step can overshoot
         * a bound-projected point, so fall back to the pre-polish value if the
         * polish came out worse by Deb's rules (never worsens the chain).
         *
         * Skipped when the caller disabled polishing with "PostProcess" -> False:
         * then the per-chain raw best is carried straight through and the global
         * best over chains is the global raw best — bit-for-bit what a single
         * global-raw-best pass over the same RNG sequence produced before. */
        double fp = fc, pp = pc;
        for (size_t j = 0; j < n; j++) x[j] = xc[j];   /* x reused as polish buffer */
        if (nc->post_process != 0) {
            nm_local_polish(D, x, &fp, &pp);
            if (nm_better(fc, pc, fp, pp)) {
                for (size_t j = 0; j < n; j++) x[j] = xc[j];
                fp = fc; pp = pc;
            }
        }
        if (!have || nm_better(fp, pp, *fbest, *penbest)) {
            for (size_t j = 0; j < n; j++) xbest[j] = x[j];
            *fbest = fp; *penbest = pp; have = true;
        }
    }
    free(x); free(xn); free(xc);
}

/* ------------------------------------------------------------------ *
 *  Option / method / variable / constraint parsing                    *
 * ------------------------------------------------------------------ */

static bool nm_method_from_string(const char* s, int* out) {
    if (strcmp(s, "DifferentialEvolution") == 0) { *out = NM_DE;           return true; }
    if (strcmp(s, "NelderMead") == 0)            { *out = NM_NELDERMEAD;   return true; }
    if (strcmp(s, "RandomSearch") == 0)          { *out = NM_RANDOMSEARCH; return true; }
    if (strcmp(s, "SimulatedAnnealing") == 0)    { *out = NM_SA;           return true; }
    return false;
}

/* The Method sub-option LHS may be a string ("SearchPoints") or a symbol. */
static const char* nm_option_name(Expr* e) {
    if (e->type == EXPR_STRING) return e->data.string;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name;
    return NULL;
}

static bool nm_parse_method(Expr* rhs, NmConfig* nc, const char* fn) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
        nc->method = NM_AUTO;
        return true;
    }
    if (rhs->type == EXPR_STRING) {
        int m;
        if (nm_method_from_string(rhs->data.string, &m)) { nc->method = m; return true; }
        fm_warn(fn, "nimpl", "Method \"%s\" is not supported", rhs->data.string);
        return false;
    }
    if (nm_is_head(rhs, SYM_List) && rhs->data.function.arg_count >= 1) {
        Expr* h = rhs->data.function.args[0];
        int m;
        if (h->type != EXPR_STRING || !nm_method_from_string(h->data.string, &m)) {
            fm_warn(fn, "badmeth", "Method list must begin with a method-name string");
            return false;
        }
        nc->method = m;
        for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
            Expr* r = rhs->data.function.args[i];
            if (!nm_is_head(r, SYM_Rule) && !nm_is_head(r, SYM_RuleDelayed)) continue;
            if (r->data.function.arg_count != 2) continue;
            const char* on = nm_option_name(r->data.function.args[0]);
            Expr* ov = r->data.function.args[1];
            if (!on) continue;
            if (strcmp(on, "SearchPoints") == 0) {
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0)
                    nc->search_points = (int)ov->data.integer;
            } else if (strcmp(on, "ScalingFactor") == 0) {
                /* DE differential weight F, mutation scale in DE/rand/1. A real
                 * in (0, 2]; an out-of-range or non-real value warns and keeps
                 * the default 0.6. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0 && dv <= 2.0)
                    nc->F = dv;
                else
                    fm_warn(fn, "sopt",
                            "ScalingFactor must be a real in (0, 2]; using the default");
            } else if (strcmp(on, "CrossProbability") == 0) {
                /* DE crossover probability CR, the per-coordinate chance of
                 * taking the mutant. A real in [0, 1]; an out-of-range or
                 * non-real value warns and keeps the default 0.9. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv >= 0.0 && dv <= 1.0)
                    nc->CR = dv;
                else
                    fm_warn(fn, "sopt",
                            "CrossProbability must be a real in [0, 1]; using the default");
            } else if (strcmp(on, "RandomSeed") == 0) {
                if (ov->type == EXPR_INTEGER && ov->data.integer >= 0)
                    nc->seed = (uint64_t)ov->data.integer;
            } else if (strcmp(on, "ReflectRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->reflect_ratio = dv;
            } else if (strcmp(on, "ExpandRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->expand_ratio = dv;
            } else if (strcmp(on, "ContractRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->contract_ratio = dv;
            } else if (strcmp(on, "ShrinkRatio") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->shrink_ratio = dv;
            } else if (strcmp(on, "Tolerance") == 0) {
                double dv; if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->tolerance = dv;
            } else if (strcmp(on, "PostProcess") == 0) {
                /* Every Mathematica setting is accepted. True / Automatic, a
                 * named local method as a string ("InteriorPoint",
                 * "FindMinimum", "KKT", ...), or {"Method", subopts...} enable
                 * the final exact polish; False / None disable it. Mathilda has a
                 * single FindMinimum-style local polish (BFGS for continuous/box,
                 * quadratic penalty for general constraints) that already selects
                 * the right inner solver for the problem, so a named method turns
                 * post-processing on rather than picking a distinct algorithm. */
                if (ov->type == EXPR_STRING) {
                    nc->post_process = 1;
                } else if (ov->type == EXPR_SYMBOL) {
                    const char* s = ov->data.symbol.name;
                    if (s == SYM_True || s == SYM_Automatic)      nc->post_process = 1;
                    else if (s == SYM_False || s == SYM_None)     nc->post_process = 0;
                    else fm_warn(fn, "pmeth",
                                 "PostProcess -> %s not recognised; using Automatic", s);
                } else if (nm_is_head(ov, SYM_List)
                           && ov->data.function.arg_count >= 1
                           && ov->data.function.args[0]->type == EXPR_STRING) {
                    nc->post_process = 1;           /* {"InteriorPoint", opts...} */
                } else {
                    fm_warn(fn, "pmeth", "invalid PostProcess value; using Automatic");
                }
            } else if (strcmp(on, "InitialPoints") == 0) {
                /* A list of starting points {{x1,...}, {x2,...}, ...}; borrowed
                 * and validated/consumed by the engine, where the dimension n is
                 * known. Anything else is ignored (falls back to random starts). */
                if (nm_is_head(ov, SYM_List) && ov->data.function.arg_count > 0)
                    nc->init_points = ov;
            } else if (strcmp(on, "PerturbationScale") == 0) {
                /* SimulatedAnnealing: multiplies the size of the random step
                 * used to generate a new trial point (default 1.0). A larger
                 * scale explores more widely; must be a positive finite real. */
                double dv;
                if (fm_expr_to_double_real(ov, &dv) && isfinite(dv) && dv > 0.0)
                    nc->perturb_scale = dv;
                else
                    fm_warn(fn, "sopt",
                            "PerturbationScale must be a positive real; using 1.0");
            } else if (strcmp(on, "LevelIterations") == 0) {
                /* SimulatedAnnealing: the number of trial moves spent at each
                 * temperature level. The per-chain iteration budget is
                 * MaxIterations * LevelIterations (default 50, the value the
                 * previous fixed multiplier hard-coded), so this scales how long
                 * each chain anneals. An explicit value is honored verbatim, like
                 * "SearchPoints"; Automatic / None restore the default. */
                if (ov->type == EXPR_INTEGER && ov->data.integer > 0
                    && ov->data.integer <= 2000000000LL)
                    nc->level_iterations = (int)ov->data.integer;
                else if (ov->type == EXPR_SYMBOL
                         && (ov->data.symbol.name == SYM_Automatic
                             || ov->data.symbol.name == SYM_None))
                    nc->level_iterations = 0;
                else
                    fm_warn(fn, "sopt",
                            "LevelIterations must be a positive integer; using Automatic");
            } else if (strcmp(on, "BoltzmannExponent") == 0) {
                /* SimulatedAnnealing: the exponent of the Metropolis acceptance
                 * probability for an uphill move. A function f is called as
                 * f[i, df, f0] (iteration i≥1, objective difference df≥0, current
                 * value f0) and the point is accepted with probability Exp[f[...]];
                 * Automatic / None keep the built-in geometric-cooling exponent
                 * -df/T. Borrowed from the held method list. */
                if (ov->type == EXPR_SYMBOL
                    && (ov->data.symbol.name == SYM_Automatic
                        || ov->data.symbol.name == SYM_None)) {
                    nc->boltzmann_fn = NULL;
                } else if (ov->type == EXPR_FUNCTION || ov->type == EXPR_SYMBOL) {
                    nc->boltzmann_fn = ov;
                } else {
                    fm_warn(fn, "bexp",
                            "invalid BoltzmannExponent value; using Automatic");
                }
            } else if (strcmp(on, "PenaltyFunction") == 0) {
                /* A function applied to each constraint's violation to score
                 * infeasible points during the global search; Automatic ≡ #^2 &.
                 * Automatic / None keep the built-in squared penalty; a pure
                 * function or a function symbol (#^2 &, (10 #) &, Sqrt, ...) is
                 * stored (borrowed from the held method list) and applied in
                 * nm_eval_pen. It affects only the global-search feasibility
                 * scoring — the final local polish keeps the differentiable
                 * squared penalty its analytic gradient assumes. */
                if (ov->type == EXPR_SYMBOL
                    && (ov->data.symbol.name == SYM_Automatic
                        || ov->data.symbol.name == SYM_None)) {
                    nc->penalty_fn = NULL;
                } else if (ov->type == EXPR_FUNCTION || ov->type == EXPR_SYMBOL) {
                    nc->penalty_fn = ov;
                } else {
                    fm_warn(fn, "penf",
                            "invalid PenaltyFunction value; using Automatic");
                }
            }
        }
        return true;
    }
    fm_warn(fn, "badmeth", "invalid Method value");
    return false;
}

static bool nm_apply_option(Expr* rule, FmOpts* opts, NmConfig* nc, const char* fn) {
    Expr* lhs = rule->data.function.args[0];
    Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol.name;
    if (name == SYM_Method) return nm_parse_method(rhs, nc, fn);
    if (name == SYM_WorkingPrecision) {
        if (!fm_parse_working_precision(rhs, &opts->prec_mode, &opts->wp_bits)) {
            fm_warn(fn, "badopt", "invalid WorkingPrecision value");
            return false;
        }
        return true;
    }
    if (name == SYM_MaxIterations) {
        if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
            return true;                       /* keep the NMinimize default */
        if (rhs->type == EXPR_INTEGER && rhs->data.integer > 0) {
            opts->max_iter = rhs->data.integer;
            opts->max_iter_set = true;
            return true;
        }
        fm_warn(fn, "badopt", "MaxIterations must be a positive integer or Automatic");
        return false;
    }
    if (name == SYM_AccuracyGoal)  return fm_parse_goal(rhs, &opts->acc_goal_digits);
    if (name == SYM_PrecisionGoal) return fm_parse_goal(rhs, &opts->prec_goal_digits);
    if (name == SYM_EvaluationMonitor) { opts->eval_monitor = rhs; return true; }
    if (name == SYM_StepMonitor)       { opts->step_monitor = rhs; return true; }
    if (name == SYM_Gradient)          return true;   /* accepted, unused */
    fm_warn(fn, "badopt", "unrecognised option");
    return false;
}

/* Variable set: symbols (borrowed), integer-domain mask, and per-variable
 * search-interval hints parsed from {x, lo, hi} / {x, x0, lo, hi} specs. */
typedef struct {
    size_t  n;
    Expr**  vars;
    bool*   is_int;
    double* rlo;
    double* rhi;
    bool*   has_rlo;
    bool*   has_rhi;
    bool    any_int;
} NmVarSet;

static void nm_varset_free(NmVarSet* vs) {
    free(vs->vars); free(vs->is_int);
    free(vs->rlo); free(vs->rhi); free(vs->has_rlo); free(vs->has_rhi);
    vs->vars = NULL;
}

/* Parse one variable spec element (bare symbol, {x,...} list, or
 * Element[x, Integers|Reals]). var_out is borrowed from the input tree. */
static bool nm_one_var(Expr* sub, Expr** var_out, bool* is_int_out,
                       bool* has_lo, double* lo, bool* has_hi, double* hi,
                       const char* fn) {
    *is_int_out = false; *has_lo = false; *has_hi = false;
    if (nm_is_head(sub, SYM_Element) && sub->data.function.arg_count == 2) {
        Expr* v = sub->data.function.args[0];
        Expr* dom = sub->data.function.args[1];
        if (v->type != EXPR_SYMBOL) {
            fm_warn(fn, "ivar", "Element variable must be a symbol");
            return false;
        }
        if (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Integers) {
            *var_out = v; *is_int_out = true; return true;
        }
        if (dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals) {
            *var_out = v; return true;
        }
        fm_warn(fn, "nimpl", "only the Integers and Reals domains are supported");
        return false;
    }
    /* Indexed variable atom (x[i], x[i,j], ...) from an expanded Table/Array
     * spec: an unbounded variable, no starting interval. Rewritten to a fresh
     * scalar symbol by the driver before the solver machinery runs. */
    if (sub->type == EXPR_FUNCTION
        && sub->data.function.head->type == EXPR_SYMBOL
        && sub->data.function.head->data.symbol.name != SYM_List) {
        *var_out = sub;
        return true;
    }
    Expr *u, *x0 = NULL, *x1 = NULL, *xmn = NULL, *xmx = NULL;
    FmSpecKind k = fm_parse_var_spec(sub, &u, &x0, &x1, &xmn, &xmx);
    bool ok = (k != FM_SPEC_BAD);
    if (ok) {
        *var_out = u;
        double a, b;
        if (k == FM_SPEC_TWO_START && x0 && x1
            && fm_expr_to_double_real(x0, &a) && fm_expr_to_double_real(x1, &b)) {
            if (a > b) { double t = a; a = b; b = t; }
            *has_lo = true; *lo = a; *has_hi = true; *hi = b;
        } else if (k == FM_SPEC_BRACKET && xmn && xmx
            && fm_expr_to_double_real(xmn, &a) && fm_expr_to_double_real(xmx, &b)) {
            if (a > b) { double t = a; a = b; b = t; }
            *has_lo = true; *lo = a; *has_hi = true; *hi = b;
        }
    } else {
        fm_warn(fn, "ivar", "variable specification malformed");
    }
    expr_free(x0); expr_free(x1); expr_free(xmn); expr_free(xmx);
    return ok;
}

static bool nm_parse_vars(Expr* var_arg, NmVarSet* vs, const char* fn) {
    bool system = false;
    size_t na = 0;
    if (nm_is_head(var_arg, SYM_List) && var_arg->data.function.arg_count > 0) {
        na = var_arg->data.function.arg_count;
        bool any_sub = false, all_subsym = true, all_atom = true;
        for (size_t i = 0; i < na; i++) {
            Expr* e = var_arg->data.function.args[i];
            bool inner = nm_is_head(e, SYM_List);
            bool elem  = nm_is_head(e, SYM_Element);
            bool atom  = nm_is_var_atom(e);         /* symbol or x[i] */
            if (inner || elem) any_sub = true;
            if (!(inner || elem || atom)) all_subsym = false;
            if (!atom) all_atom = false;
        }
        if (any_sub && all_subsym) system = true;
        else if (all_atom)         system = true;   /* {x, y, ...} / {x[1], x[2], ...} */
    }

    size_t n = system ? na : 1;
    vs->n = n;
    vs->vars    = (Expr**)calloc(n, sizeof(Expr*));
    vs->is_int  = (bool*)calloc(n, sizeof(bool));
    vs->rlo     = (double*)calloc(n, sizeof(double));
    vs->rhi     = (double*)calloc(n, sizeof(double));
    vs->has_rlo = (bool*)calloc(n, sizeof(bool));
    vs->has_rhi = (bool*)calloc(n, sizeof(bool));
    vs->any_int = false;

    bool ok = true;
    if (system) {
        for (size_t i = 0; i < n && ok; i++) {
            ok = nm_one_var(var_arg->data.function.args[i], &vs->vars[i],
                            &vs->is_int[i], &vs->has_rlo[i], &vs->rlo[i],
                            &vs->has_rhi[i], &vs->rhi[i], fn);
        }
    } else {
        ok = nm_one_var(var_arg, &vs->vars[0], &vs->is_int[0],
                        &vs->has_rlo[0], &vs->rlo[0],
                        &vs->has_rhi[0], &vs->rhi[0], fn);
    }
    if (ok) {
        for (size_t i = 0; i < n; i++)
            if (!nm_is_var_atom(vs->vars[i])) { ok = false; break; }
    }
    if (ok) for (size_t i = 0; i < n; i++) if (vs->is_int[i]) vs->any_int = true;
    if (!ok) nm_varset_free(vs);
    return ok;
}

/* Pull Element[x, Integers|Reals] domain declarations out of the constraint
 * tree (marking is_int for Integers), returning the remaining constraint
 * expression (owned) or NULL if none remain. The declared operand may be a
 * single variable (Element[x, Integers]) or a set of them written as x|y|...
 * (Alternatives) or {x, y, ...} (List): the domain applies to every member,
 * matching Mathematica. Unsupported Element domains — and declarations naming
 * a non-optimization symbol — are left in place so fm_collect_constraints
 * rejects them with its own message. */
static Expr* nm_filter_int(Expr* cons, Expr** vars, size_t n, bool* is_int) {
    if (!cons) return NULL;
    if (nm_is_head(cons, SYM_And)) {
        size_t cnt = cons->data.function.arg_count;
        Expr** kids = (Expr**)malloc(sizeof(Expr*) * (cnt ? cnt : 1));
        size_t m = 0;
        for (size_t i = 0; i < cnt; i++) {
            Expr* c = nm_filter_int(cons->data.function.args[i], vars, n, is_int);
            if (c) kids[m++] = c;
        }
        Expr* r;
        if (m == 0)      { r = NULL; }
        else if (m == 1) { r = kids[0]; }
        else             { r = expr_new_function(expr_new_symbol(SYM_And), kids, m); }
        free(kids);
        return r;
    }
    if (nm_is_head(cons, SYM_Element) && cons->data.function.arg_count == 2) {
        Expr* v   = cons->data.function.args[0];
        Expr* dom = cons->data.function.args[1];
        /* Only the Integers and Reals domains are absorbed here; any other
         * domain (or a non-symbol domain) falls through to be left in place. */
        if (dom->type == EXPR_SYMBOL
            && (dom->data.symbol.name == SYM_Integers
                || dom->data.symbol.name == SYM_Reals)) {
            bool integers = (dom->data.symbol.name == SYM_Integers);
            /* Operand list: the bare symbol, or the members of an
             * Alternatives / List container. */
            Expr** ops; size_t nops;
            if (v->type == EXPR_SYMBOL) { ops = &v; nops = 1; }
            else if (nm_is_head(v, SYM_Alternatives) || nm_is_head(v, SYM_List)) {
                ops = v->data.function.args; nops = v->data.function.arg_count;
            } else { ops = NULL; nops = 0; }
            /* Absorb the declaration only when every operand is one of the
             * optimization variables; otherwise leave the whole node in place
             * (a domain assertion on some other symbol is not enforceable). */
            bool all_vars = (nops > 0);
            for (size_t k = 0; k < nops && all_vars; k++) {
                if (ops[k]->type != EXPR_SYMBOL) { all_vars = false; break; }
                bool found = false;
                for (size_t i = 0; i < n; i++)
                    if (vars[i]->data.symbol.name == ops[k]->data.symbol.name) {
                        found = true; break;
                    }
                all_vars = found;
            }
            if (all_vars) {
                if (integers)
                    for (size_t k = 0; k < nops; k++)
                        for (size_t i = 0; i < n; i++)
                            if (vars[i]->data.symbol.name == ops[k]->data.symbol.name)
                                is_int[i] = true;
                return NULL;   /* domain declaration absorbed */
            }
        }
        return expr_copy(cons);
    }
    return expr_copy(cons);
}

/* ------------------------------------------------------------------ *
 *  Result construction                                                *
 * ------------------------------------------------------------------ */

static Expr* nm_build_result(double fmin, Expr** vars, const double* vals,
                             const bool* is_int, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* val = is_int[i] ? expr_new_integer((int64_t)llround(vals[i]))
                              : expr_new_real(vals[i]);
        Expr* r_args[2] = { expr_copy(vars[i]), val };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_real(fmin), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}

static Expr* nm_build_infeasible(Expr** vars, size_t n) {
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* r_args[2] = { expr_copy(vars[i]), expr_new_symbol(SYM_Indeterminate) };
        rules[i] = expr_new_function(expr_new_symbol(SYM_Rule), r_args, 2);
    }
    Expr* rule_list = expr_new_function(expr_new_symbol(SYM_List), rules, n);
    free(rules);
    Expr* top_args[2] = { expr_new_symbol(SYM_Infinity), rule_list };
    return expr_new_function(expr_new_symbol(SYM_List), top_args, 2);
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

static Expr* nm_minimize_driver(Expr* res, const char* fn_name) {
    g_fm_name = fn_name;
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn(fn_name, "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    /* Peel trailing options (same recogniser FindMinimum uses). */
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
    opts.max_iter = 100;                       /* NMinimize default */
    opts.max_iter_set = false;
    opts.acc_goal_digits = -1.0;
    opts.prec_goal_digits = -1.0;
    opts.gradient = NULL;
    opts.step_monitor = NULL;
    opts.eval_monitor = NULL;

    NmConfig nc;
    nc.method = NM_AUTO;
    nc.search_points = 0;
    nc.F = -1.0;
    nc.CR = -1.0;
    nc.reflect_ratio = -1.0;
    nc.expand_ratio = -1.0;
    nc.contract_ratio = -1.0;
    nc.shrink_ratio = -1.0;
    nc.tolerance = -1.0;
    nc.post_process = -1;
    nc.init_points = NULL;
    nc.penalty_fn = NULL;
    nc.perturb_scale = -1.0;
    nc.boltzmann_fn = NULL;
    nc.level_iterations = 0;
    nc.seed = NM_DEFAULT_SEED;

    for (size_t i = pos_end; i < argc; i++) {
        if (!nm_apply_option(res->data.function.args[i], &opts, &nc, fn_name))
            return NULL;
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

    /* Parse the variables first — they do not depend on the objective — so a
     * problem passed as a bare symbol (prob = {f, cons}; NMinimize[prob, vars])
     * can be resolved below with those variables protected. */
    Expr* var_arg = res->data.function.args[1];

    /* Parse the variable specification. NMinimize is not HoldAll, so a generator
     * such as Table[x[i], {i, 1, 10}] or Array[x, 10] has already expanded to the
     * variable list before we run, and a symbol bound to a list has already
     * resolved. The eval_and_free calls below are therefore idempotent
     * normalizers (they no-op on an already-evaluated {...}/Element, and still
     * cover a spec reaching the driver unevaluated). A {...} list or Element[...]
     * is used directly. A bare symbol is ambiguous: it may be a single
     * optimization variable (NMinimize[f, x], x unbound) or a symbol that still
     * resolves to a variable list — evaluate it, and if it yields a List/Element
     * use that, otherwise keep the symbol itself as the single variable. */
    NmVarSet vs;
    Expr* var_list_eval = NULL;
    {
        Expr* var_spec = var_arg;
        if (var_arg->type == EXPR_SYMBOL) {
            Expr* ev = eval_and_free(expr_copy(var_arg));
            if (ev && (nm_is_head(ev, SYM_List) || nm_is_head(ev, SYM_Element))) {
                var_list_eval = ev;
                var_spec = ev;
            } else {
                expr_free(ev);          /* unbound / scalar: use the symbol */
            }
        } else if (!nm_is_head(var_arg, SYM_List) && !nm_is_head(var_arg, SYM_Element)) {
            var_list_eval = eval_and_free(expr_copy(var_arg));
            var_spec = var_list_eval;
        }
        if (!var_spec || !nm_parse_vars(var_spec, &vs, fn_name)) {
            expr_free(var_list_eval);
            return NULL;
        }
    }
    size_t n = vs.n;

    /* Effective (scalar-symbol) variables for the solver machinery, and the
     * original variable expressions for the result rules. Indexed variables
     * (x[1], x[2], ...) are rewritten to fresh scalar symbols so the entire
     * symbol-keyed solver applies unchanged; plain symbols pass through. */
    bool indexed = false;
    for (size_t i = 0; i < n; i++)
        if (vs.vars[i]->type != EXPR_SYMBOL) indexed = true;

    Expr**       eff_vars  = (Expr**)calloc(n, sizeof(Expr*));
    Expr**       orig_vars = (Expr**)calloc(n, sizeof(Expr*));
    const char** synth     = (const char**)calloc(n, sizeof(char*));
    const char** heads      = (const char**)malloc(sizeof(char*) * (n ? n : 1));
    size_t       nheads     = 0;
    for (size_t i = 0; i < n; i++) {
        orig_vars[i] = expr_copy(vs.vars[i]);
        if (vs.vars[i]->type == EXPR_SYMBOL) {
            eff_vars[i] = expr_copy(vs.vars[i]);
        } else {
            eff_vars[i] = nm_fresh_symbol();
            synth[i]    = eff_vars[i]->data.symbol.name;
        }
        const char* hn = (vs.vars[i]->type == EXPR_SYMBOL)
            ? vs.vars[i]->data.symbol.name
            : vs.vars[i]->data.function.head->data.symbol.name;
        bool seen = false;
        for (size_t j = 0; j < nheads; j++) if (heads[j] == hn) { seen = true; break; }
        if (!seen) heads[nheads++] = hn;
    }

    /* Resolve and split {f, cons}. NMinimize is not HoldAll, so the problem
     * argument normally arrives already in structural form (an inline {f, cons}
     * list or a scalar objective) and is used directly. The bare-symbol branch
     * below is a defensive resolver: if the argument still reaches the driver as
     * a symbol, it is evaluated once with the variable heads localized — so its
     * structure is exposed without capturing any global variable values — and the
     * {f, cons} list is then split. */
    Expr* f_arg  = res->data.function.args[0];
    Expr* f_eval = NULL;               /* owned resolved objective, or NULL     */
    if (f_arg->type == EXPR_SYMBOL) {
        NmHeadSave* hs = (NmHeadSave*)calloc(nheads ? nheads : 1, sizeof(NmHeadSave));
        nm_heads_localize(hs, heads, nheads);
        f_eval = eval_and_free(expr_copy(f_arg));
        nm_heads_restore(hs, nheads);
        free(hs);
        if (f_eval) f_arg = f_eval;
    }
    Expr* f_raw = f_arg;
    Expr* cons = NULL;         /* borrowed, or points at cons_built           */
    Expr* cons_built = NULL;   /* owned And[...] for a multi-constraint list  */
    if (nm_is_head(f_arg, SYM_List) && f_arg->data.function.arg_count >= 2) {
        /* {f, cons} or {f, c1, c2, ...}: the objective is the first element
         * and every remaining element is a constraint, implicitly And-ed. */
        f_raw = f_arg->data.function.args[0];
        size_t ncons = f_arg->data.function.arg_count - 1;
        if (ncons == 1) {
            cons = f_arg->data.function.args[1];
        } else {
            Expr** cc = (Expr**)malloc(sizeof(Expr*) * ncons);
            for (size_t i = 0; i < ncons; i++)
                cc[i] = expr_copy(f_arg->data.function.args[1 + i]);
            cons_built = expr_new_function(expr_new_symbol(SYM_And), cc, ncons);
            free(cc);
            cons = cons_built;
        }
    }

    /* Expand a held Table/Sum constraint and/or rewrite indexed variables. The
     * objective stays held for the plain path (evaluated per point); it is
     * pre-expanded only when we must rewrite indexed vars inside it. */
    Expr* f_eff = f_raw;   bool f_owned = false;
    Expr* cons_eff = cons; bool cons_owned = false;
    bool infeasible_pre = false;
    bool expand_cons = (cons != NULL) && !nm_is_constraint_tree(cons);
    if (indexed || expand_cons) {
        NmHeadSave* hs = (NmHeadSave*)calloc(nheads ? nheads : 1, sizeof(NmHeadSave));
        nm_heads_localize(hs, heads, nheads);

        if (cons) {
            Expr* ce = eval_and_free(expr_copy(cons));   /* expand, vars free */
            if (ce && ce->type == EXPR_SYMBOL && ce->data.symbol.name == SYM_True) {
                expr_free(ce); ce = NULL;
            } else if (ce && ce->type == EXPR_SYMBOL && ce->data.symbol.name == SYM_False) {
                expr_free(ce); ce = NULL; infeasible_pre = true;
            } else if (nm_is_head(ce, SYM_List)) {
                /* A list of constraints is an implicit And; drop trivially-true
                 * entries, and a False entry makes the whole system infeasible. */
                size_t m = ce->data.function.arg_count;
                Expr** cc = (Expr**)malloc(sizeof(Expr*) * (m ? m : 1));
                size_t kept = 0;
                for (size_t i = 0; i < m; i++) {
                    Expr* el = ce->data.function.args[i];
                    if (el->type == EXPR_SYMBOL && el->data.symbol.name == SYM_True) continue;
                    if (el->type == EXPR_SYMBOL && el->data.symbol.name == SYM_False) {
                        infeasible_pre = true; continue;
                    }
                    cc[kept++] = expr_copy(el);
                }
                expr_free(ce);
                if (kept == 0)      { ce = NULL; }
                else if (kept == 1) { ce = cc[0]; }
                else                { ce = expr_new_function(expr_new_symbol(SYM_And), cc, kept); }
                free(cc);
            }
            if (ce && indexed) {
                Expr* cs = nm_subst(ce, vs.vars, eff_vars, n);
                expr_free(ce); ce = cs;
            }
            cons_eff = ce; cons_owned = true;
        }

        if (indexed) {
            Expr* fe = eval_and_free(expr_copy(f_raw));
            f_eff = nm_subst(fe, vs.vars, eff_vars, n);
            expr_free(fe);
            f_owned = true;
        }

        nm_heads_restore(hs, nheads);
        free(hs);
    }

    /* Bind variables (Block semantics). */
    FmVarBind* binds = (FmVarBind*)calloc(n, sizeof(FmVarBind));
    for (size_t i = 0; i < n; i++)
        fm_bind_snapshot(&binds[i], eff_vars[i]->data.symbol.name);

    FmBox* boxes = (FmBox*)calloc(n, sizeof(FmBox));
    FmGenCon* gens = NULL;
    size_t ngens = 0, gcap = 0;
    FmDisjunction* disj = NULL;
    size_t ndisj = 0, dcap = 0;
    Expr** g_exprs = NULL;
    Expr* cons2 = NULL;
    double* reg_lo = NULL;
    double* reg_hi = NULL;
    bool*   used_default = NULL;   /* dim used the default +-SPAN (fully unbounded) */
    double* xbest = NULL;
    Expr* result_out = NULL;
    CompiledProgram*  f_prog = NULL;   /* compiled objective (machine prec)     */
    CompiledProgram** g_progs = NULL;  /* compiled general constraints          */

    /* Extract integer/real domain declarations, then collect the remaining
     * constraints into boxes + general FmGenCon[] + disjunctions. */
    if (cons_eff) {
        cons2 = nm_filter_int(cons_eff, eff_vars, n, vs.is_int);
        for (size_t i = 0; i < n; i++) if (vs.is_int[i]) vs.any_int = true;
        if (cons2 && !fm_collect_constraints(cons2, eff_vars, n, boxes,
                                             &gens, &ngens, &gcap,
                                             &disj, &ndisj, &dcap))
            goto cleanup;
        for (size_t k = 0; k < ngens; k++)
            gens[k].grad_exprs = fm_compute_gradient(gens[k].expr, eff_vars, n);
    }

    /* Resolve the per-dimension search box: box constraints tighten it,
     * else a starting-interval hint, else a default span. Contradictory box
     * bounds (lo > hi, e.g. x > 2 && x < 1) mean an empty feasible set. */
    bool infeasible_box = false;
    reg_lo = (double*)malloc(sizeof(double) * n);
    reg_hi = (double*)malloc(sizeof(double) * n);
    used_default = (bool*)calloc(n ? n : 1, sizeof(bool));
    bool any_default = false;
    for (size_t i = 0; i < n; i++) {
        bool klo = boxes[i].has_lo || vs.has_rlo[i];
        bool khi = boxes[i].has_hi || vs.has_rhi[i];
        double lo = boxes[i].has_lo ? boxes[i].lo : vs.rlo[i];
        double hi = boxes[i].has_hi ? boxes[i].hi : vs.rhi[i];
        if (boxes[i].has_lo && boxes[i].has_hi && boxes[i].lo > boxes[i].hi)
            infeasible_box = true;
        if (klo && khi)      { reg_lo[i] = lo; reg_hi[i] = hi; }
        else if (klo)        { reg_lo[i] = lo; reg_hi[i] = lo + NM_BOUND_SPAN; }
        else if (khi)        { reg_hi[i] = hi; reg_lo[i] = hi - NM_BOUND_SPAN; }
        else                 { reg_lo[i] = -NM_DEFAULT_SPAN; reg_hi[i] = NM_DEFAULT_SPAN;
                               used_default[i] = true; any_default = true; }
        if (reg_hi[i] <= reg_lo[i]) {
            double m = 0.5 * (reg_lo[i] + reg_hi[i]);
            reg_lo[i] = m - 0.5; reg_hi[i] = m + 0.5;
        }
        if (vs.is_int[i]) { reg_lo[i] = floor(reg_lo[i]); reg_hi[i] = ceil(reg_hi[i]); }
    }

    /* Objective gradient (for the continuous local polish; NULL → FD). */
    if (!vs.any_int) g_exprs = fm_compute_gradient(f_eff, eff_vars, n);

    /* Machine-precision auto-compilation: the global search evaluates the
     * objective (and each general constraint) at hundreds–thousands of trial
     * points, so lowering them to bytecode over the effective variables once and
     * running the register machine per point is far cheaper than the interpreter
     * (expr_copy + evaluate + numericalize each call). The variables are already
     * unbound (Block snapshot cleared their OwnValues), so they compile as the
     * argument symbols; COMPILE_FOLD_GLOBALS folds any other machine-valued
     * symbol — safe because these programs live only for this call. A body with
     * a construct Compile can't lower stays NULL and uses the interpreter, and
     * every per-point call falls back to the interpreter on a domain/non-finite
     * result, so this is a pure speedup with no change in answer. MPFR
     * (WorkingPrecision > MachinePrecision) keeps the exact interpreter path. */
    if (opts.prec_mode == FM_PREC_MACHINE) {
        const char** cnames = (const char**)malloc(sizeof(char*) * (n ? n : 1));
        CompileType* ctypes = (CompileType*)malloc(sizeof(CompileType) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            cnames[i] = eff_vars[i]->data.symbol.name;
            ctypes[i] = CT_REAL;
        }
        f_prog = compile_expr_ex(f_eff, cnames, ctypes, n, COMPILE_FOLD_GLOBALS);
        if (f_prog && compiled_result_type(f_prog) != CT_REAL) {
            compiled_free(f_prog); f_prog = NULL;
        }
        if (ngens > 0) {
            g_progs = (CompiledProgram**)calloc(ngens, sizeof(CompiledProgram*));
            for (size_t k = 0; k < ngens; k++) {
                CompiledProgram* p = compile_expr_ex(gens[k].expr, cnames, ctypes,
                                                     n, COMPILE_FOLD_GLOBALS);
                if (p && compiled_result_type(p) != CT_REAL) { compiled_free(p); p = NULL; }
                g_progs[k] = p;
            }
        }
        free(cnames); free(ctypes);
    }

    NmDriver D;
    D.f_raw = f_eff; D.vars = eff_vars; D.n = n; D.binds = binds;
    D.g_exprs = g_exprs; D.gens = gens; D.ngens = ngens; D.boxes = boxes;
    D.opts = &opts; D.is_int = vs.is_int; D.any_int = vs.any_int;
    D.reg_lo = reg_lo; D.reg_hi = reg_hi;
    D.f_prog = f_prog; D.g_progs = g_progs;
    D.penalty_fn = nc.penalty_fn;
    D.disj = disj; D.ndisj = ndisj;

    /* Serve the local polish's objective evaluations from the compiled program
     * instead of the interpreter. RandomSearch runs one local solve per
     * SearchPoint, so this is the difference between 1000 compiled polishes and
     * 1000 interpreter polishes (see g_fm_obj_* and fm_eval_scalar). Cleared to
     * the inactive sentinel below, before f_prog is freed. A plain reset is
     * enough rather than save/restore: the registration is live only while
     * f_prog != NULL, and an objective compilable enough to have a program
     * cannot contain a nested NMinimize call (NMinimize is not a compilable
     * head), so no active registration can ever be clobbered by re-entry. */
    g_fm_obj_expr  = f_eff;
    g_fm_obj_prog  = f_prog;
    g_fm_obj_nargs = n;

    xbest = (double*)malloc(sizeof(double) * n);
    double fbest = 1e300, penbest = 1e300;
    int method = (nc.method == NM_AUTO) ? NM_DE : nc.method;
    bool do_post = (nc.post_process != 0);

    /* Adaptive search-region expansion. If the default +-SPAN sampling region
     * contains no feasible point, the fully-unbounded coordinates are grown by
     * successive powers of ten and the search is retried — so a feasible region
     * whose location is implied by nonlinear constraints rather than stated as
     * variable bounds (e.g. the pressure-vessel MINLP, feasible near x3 ~ 52) is
     * still found instead of reporting {Infinity, ...}. Only fully-unbounded
     * coordinates grow; a coordinate carrying a box bound or a starting-interval
     * hint keeps its resolved region. Attempt 0 is the base region with the base
     * seed, so a problem already feasible there is solved identically to before;
     * expansion triggers only to rescue infeasibility, and stops as soon as a
     * feasible point is found (the smallest region that yields feasibility,
     * which keeps the search from drifting into far, non-physical basins). */
    int max_attempt = (any_default && !infeasible_box && !infeasible_pre)
                    ? NM_MAX_REGION_EXPAND : 0;
    double* xattempt = (double*)malloc(sizeof(double) * n);
    for (int attempt = 0; attempt <= max_attempt; attempt++) {
        if (attempt > 0) {
            double span = NM_DEFAULT_SPAN * pow(10.0, (double)attempt);
            for (size_t i = 0; i < n; i++) if (used_default[i]) {
                reg_lo[i] = vs.is_int[i] ? floor(-span) : -span;
                reg_hi[i] = vs.is_int[i] ? ceil(span)   :  span;
            }
        }
        double fa = 1e300, pa = 1e300;
        NmRng rng;
        nm_rng_seed(&rng, nc.seed + (uint64_t)attempt * 0x100000001B3ULL);
        switch (method) {
            case NM_NELDERMEAD:   nm_neldermead(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_RANDOMSEARCH: nm_randomsearch(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_SA:           nm_sa(&D, &nc, &rng, xattempt, &fa, &pa); break;
            case NM_DE:
            default:              nm_de(&D, &nc, &rng, xattempt, &fa, &pa); break;
        }
        /* Polish this attempt's best with the exact local solver, unless the
         * caller disabled it with "PostProcess" -> False. Guard against a
         * penalty/BFGS step that overshoots: if the polished point is worse by
         * Deb's rules than the pre-polish best, keep the latter. */
        if (do_post) {
            double* xsave = (double*)malloc(sizeof(double) * n);
            for (size_t i = 0; i < n; i++) xsave[i] = xattempt[i];
            double fsave = fa, psave = pa;
            nm_local_polish(&D, xattempt, &fa, &pa);
            if (nm_better(fsave, psave, fa, pa)) {
                for (size_t i = 0; i < n; i++) xattempt[i] = xsave[i];
                fa = fsave; pa = psave;
            }
            free(xsave);
        }
        if (attempt == 0 || nm_better(fa, pa, fbest, penbest)) {
            for (size_t i = 0; i < n; i++) xbest[i] = xattempt[i];
            fbest = fa; penbest = pa;
        }
        if (penbest <= NM_FEAS_FINAL) break;   /* feasible: stop expanding */
    }
    free(xattempt);
    bool feasible = !infeasible_box && !infeasible_pre && (penbest <= NM_FEAS_FINAL);

    /* Optional MPFR refinement for WorkingPrecision > MachinePrecision on
     * continuous, general-constraint-free problems (reuses fm_run_bfgs_mpfr). */
#ifdef USE_MPFR
    bool mpfr_built = false;
    bool want_mpfr = (opts.prec_mode == FM_PREC_MPFR);
    bool mpfr_eligible = want_mpfr && do_post && !vs.any_int && ngens == 0;
    mpfr_t* xm = NULL;
    mpfr_t fmv;
    if (want_mpfr && !mpfr_eligible && do_post)
        fm_warn(fn_name, "nimpl",
                "WorkingPrecision > MachinePrecision with general constraints or "
                "integer domains is not supported; using machine precision");
    if (feasible && mpfr_eligible) {
        Expr** gm = fm_compute_gradient(f_eff, eff_vars, n);
        xm = fm_mpfr_array(n, opts.wp_bits);
        for (size_t i = 0; i < n; i++) mpfr_set_d(xm[i], xbest[i], MPFR_RNDN);
        mpfr_init2(fmv, opts.wp_bits);
        bool saved_quiet = g_fm_quiet;
        g_fm_quiet = true;
        bool mok = fm_run_bfgs_mpfr(f_eff, eff_vars, n, binds, gm, xm, boxes, &opts, fmv);
        g_fm_quiet = saved_quiet;
        if (mok)
            mpfr_built = true;
        else { fm_mpfr_array_free(xm, n); xm = NULL; mpfr_clear(fmv); }
        if (gm) { for (size_t i = 0; i < n; i++) expr_free(gm[i]); free(gm); }
    }
#endif

    /* Free the temporary bindings so the variable symbols are unbound while
     * we build the result rules (Rule[x, v] must not re-evaluate x). */
    for (size_t i = 0; i < n; i++) fm_bind_clear_temp(&binds[i]);

#ifdef USE_MPFR
    if (mpfr_built) {
        result_out = fm_build_result_mpfr(fmv, orig_vars, (const mpfr_t*)xm, n);
        fm_mpfr_array_free(xm, n);
        mpfr_clear(fmv);
    } else
#endif
    if (feasible) result_out = nm_build_result(fbest, orig_vars, xbest, vs.is_int, n);
    else          result_out = nm_build_infeasible(orig_vars, n);

cleanup:
    for (size_t i = 0; i < n; i++) fm_bind_restore(&binds[i]);
    free(binds);
    if (g_exprs) { for (size_t i = 0; i < n; i++) expr_free(g_exprs[i]); free(g_exprs); }
    if (gens) {
        for (size_t k = 0; k < ngens; k++) {
            expr_free(gens[k].expr);
            if (gens[k].grad_exprs) {
                for (size_t i = 0; i < n; i++) expr_free(gens[k].grad_exprs[i]);
                free(gens[k].grad_exprs);
            }
        }
        free(gens);
    }
    if (disj) {
        for (size_t k = 0; k < ndisj; k++) expr_free(disj[k].expr);
        free(disj);
    }
    expr_free(cons2);
    if (cons_owned) expr_free(cons_eff);
    if (f_owned)    expr_free(f_eff);
    expr_free(cons_built);
    expr_free(f_eval);         /* resolved bare-symbol objective; frees f_raw/cons borrows */
    g_fm_obj_expr  = NULL;     /* deregister the objective before f_prog is freed */
    g_fm_obj_prog  = NULL;
    g_fm_obj_nargs = 0;
    if (f_prog) compiled_free(f_prog);
    if (g_progs) {
        for (size_t k = 0; k < ngens; k++) if (g_progs[k]) compiled_free(g_progs[k]);
        free(g_progs);
    }
    free(boxes);
    free(reg_lo);
    free(reg_hi);
    free(used_default);
    free(xbest);
    if (eff_vars)  { for (size_t i = 0; i < n; i++) expr_free(eff_vars[i]);  free(eff_vars); }
    if (orig_vars) { for (size_t i = 0; i < n; i++) expr_free(orig_vars[i]); free(orig_vars); }
    if (synth) {
        for (size_t i = 0; i < n; i++) if (synth[i]) symtab_remove_symbol(synth[i]);
        free(synth);
    }
    free(heads);
    expr_free(var_list_eval);
    nm_varset_free(&vs);
    return result_out;
}

Expr* builtin_nminimize(Expr* res) {
    return nm_minimize_driver(res, "NMinimize");
}

/* NMaximize: minimise −f and negate the reported optimum. Mirrors the
 * FindMaximum → FindMinimum wrapper above. */
Expr* builtin_nmaximize(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        fm_warn("NMaximize", "argt", "needs at least 2 arguments; got %zu", argc);
        return NULL;
    }
    Expr* f_orig = res->data.function.args[0];
    Expr* new_first;
    if (nm_is_head(f_orig, SYM_List) && f_orig->data.function.arg_count == 2) {
        Expr* inner_f = f_orig->data.function.args[0];
        Expr* cons = f_orig->data.function.args[1];
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(inner_f) };
        Expr* neg_f = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
        Expr* list_args[2] = { neg_f, expr_copy(cons) };
        new_first = expr_new_function(expr_new_symbol(SYM_List), list_args, 2);
    } else {
        Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(f_orig) };
        new_first = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    }
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * argc);
    new_args[0] = new_first;
    for (size_t i = 1; i < argc; i++) new_args[i] = expr_copy(res->data.function.args[i]);
    Expr* synthetic = expr_new_function(expr_new_symbol(SYM_NMinimize), new_args, argc);
    free(new_args);
    Expr* min_result = nm_minimize_driver(synthetic, "NMaximize");
    expr_free(synthetic);
    if (!min_result) return NULL;
    /* Negate the reported optimum value while preserving its numeric type. */
    if (min_result->type == EXPR_FUNCTION && min_result->data.function.arg_count == 2) {
        Expr* fmin_e = min_result->data.function.args[0];
#ifdef USE_MPFR
        if (fmin_e && fmin_e->type == EXPR_MPFR) {
            long bits = mpfr_get_prec(fmin_e->data.mpfr);
            mpfr_t neg; mpfr_init2(neg, bits);
            mpfr_neg(neg, fmin_e->data.mpfr, MPFR_RNDN);
            expr_free(fmin_e);
            min_result->data.function.args[0] = expr_new_mpfr_copy(neg);
            mpfr_clear(neg);
        } else
#endif
        {
            double fmin;
            if (fm_expr_to_double_real(fmin_e, &fmin)) {
                expr_free(fmin_e);
                min_result->data.function.args[0] = expr_new_real(-fmin);
            }
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
    /* NMinimize/NMaximize are Protected but NOT HoldAll (matching Mathematica,
     * whose Attributes[NMinimize] is {Protected}). Their variables are ordinary
     * unbound symbols that evaluate to themselves, and the objective is
     * re-evaluated per trial point under a Block-style binding of those symbols,
     * so holding the arguments is unnecessary — and holding them is what made a
     * Method sub-option value such as "RandomSeed" -> s (s a Do/Table iterator or
     * any expression) arrive unevaluated and get silently dropped. FindMinimum
     * stays HoldAll: its {x, x0} specs pair a variable with an initial value that
     * must not evaluate. */
    symtab_add_builtin("NMinimize", builtin_nminimize);
    symtab_get_def("NMinimize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("NMaximize", builtin_nmaximize);
    symtab_get_def("NMaximize")->attributes |= ATTR_PROTECTED;
}
