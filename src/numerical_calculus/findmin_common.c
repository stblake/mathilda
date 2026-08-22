/* findmin_common.c — shared local core: eval/grad/hessian, options, constraints, line search, Brent, Cholesky.
 * Split from the original findmin.c; shared declarations in
 * findmin_internal.h. Do not add cross-file helpers here without a
 * prototype in that header. */
#include "findmin_internal.h"


/* ------------------------------------------------------------------ *
 *  Diagnostic helper                                                  *
 * ------------------------------------------------------------------ */

/* When true, fm_warn is a no-op. NMinimize sets it around its internal
 * local-solver calls so the penalty/line-search chatter that is expected
 * during global search does not reach the user (NMinimize reports feasibility
 * itself via its {Infinity, ...} result). */
bool g_fm_quiet = false;

/* Augmented-Lagrangian multipliers for the constrained local solver. When
 * g_fm_al_lambda is non-NULL and g_fm_al_gens matches the general-constraint
 * array currently being minimized, fm_eval_augmented / fm_eval_aug_gradient use
 * the PHR (Powell–Hestenes–Rockafellar) augmented-Lagrangian terms
 *   equality  h_k:  λ_k·h_k + μ·h_k²          (always active)
 *   ineq g_k≤0:      λ_k·g_k + μ·g_k²  if s>0, else −λ_k²/(4μ),  s = λ_k+2μ·g_k
 * instead of the pure quadratic penalty μ·(Σmax(0,g)²+Σh²). This lets a moderate
 * μ reach feasibility (multipliers absorb the active-constraint gradient) rather
 * than needing μ→∞, which is ill-conditioned and strands nonlinear equalities
 * (e.g. a bilinear pooling constraint). CRUCIAL invariant: with every λ_k == 0
 * the AL terms are algebraically identical to the quadratic penalty, so a solve
 * that never updates λ (g_fm_al_lambda == NULL, the historical path) is
 * bit-for-bit unchanged. The gens-pointer identity guard keeps a nested solver
 * call on a different constraint set from mis-indexing this array. */
const double* g_fm_al_lambda = NULL;
const FmGenCon* g_fm_al_gens = NULL;

void fm_warn(const char* fn, const char* tag, const char* fmt, ...) {
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
const char* g_fm_name = "FindMinimum";

/* ------------------------------------------------------------------ *
 *  Numeric extraction / construction                                  *
 * ------------------------------------------------------------------ */

bool fm_expr_to_double_real(Expr* e, double* out) {
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

bool fm_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fm_is_known_option_name(lhs->data.symbol.name);
}

bool fm_parse_working_precision(Expr* val,
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

bool fm_parse_goal(Expr* val, double* digits_out) {
    if (val->type == EXPR_SYMBOL) {
        if (val->data.symbol.name == SYM_Automatic) { *digits_out = -1.0; return true; }
        if (val->data.symbol.name == SYM_Infinity)  { *digits_out = INFINITY; return true; }
        return false;
    }
    return fm_expr_to_double_real(val, digits_out);
}

bool fm_apply_option(Expr* rule, FmOpts* opts) {
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
            if (strcmp(s, "LBFGSB") == 0
             || strcmp(s, "LBFGS") == 0
             || strcmp(s, "LimitedMemoryBFGS") == 0) { opts->method = FM_METHOD_LBFGSB;       return true; }
            if (strcmp(s, "Powell") == 0
             || strcmp(s, "PrincipalAxis") == 0)      { opts->method = FM_METHOD_POWELL;       return true; }
            if (strcmp(s, "NelderMead") == 0)         { opts->method = FM_METHOD_NELDERMEAD;   return true; }
            if (strcmp(s, "TNC") == 0
             || strcmp(s, "TruncatedNewton") == 0)    { opts->method = FM_METHOD_TNC;          return true; }
            if (strcmp(s, "SLSQP") == 0
             || strcmp(s, "SequentialQuadraticProgramming") == 0) { opts->method = FM_METHOD_SLSQP; return true; }
            if (strcmp(s, "COBYLA") == 0)             { opts->method = FM_METHOD_COBYLA;      return true; }
            if (strcmp(s, "COBYQA") == 0)             { opts->method = FM_METHOD_COBYQA;      return true; }
            if (strcmp(s, "NewtonCG") == 0
             || strcmp(s, "Newton-CG") == 0)          { opts->method = FM_METHOD_NEWTONCG;    return true; }
            if (strcmp(s, "Dogleg") == 0
             || strcmp(s, "dogleg") == 0)             { opts->method = FM_METHOD_DOGLEG;      return true; }
            if (strcmp(s, "TrustNCG") == 0
             || strcmp(s, "trust-ncg") == 0
             || strcmp(s, "TrustRegionNewtonCG") == 0){ opts->method = FM_METHOD_TRUSTNCG;    return true; }
            if (strcmp(s, "TrustExact") == 0
             || strcmp(s, "trust-exact") == 0)        { opts->method = FM_METHOD_TRUSTEXACT;  return true; }
            if (strcmp(s, "TrustKrylov") == 0
             || strcmp(s, "trust-krylov") == 0)       { opts->method = FM_METHOD_TRUSTKRYLOV; return true; }
            if (strcmp(s, "LevenbergMarquardt") == 0
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

void fm_bind_snapshot(FmVarBind* b, const char* name) {
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
void fm_bind_clear_temp(FmVarBind* b) {
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

void fm_bind_restore(FmVarBind* b) {
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

void fm_fire_monitor(Expr* monitor) {
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
Expr* fm_eval_with_bindings(Expr* f, FmVarBind* binds,
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
Expr*            g_fm_obj_expr  = NULL;
CompiledProgram* g_fm_obj_prog  = NULL;
size_t           g_fm_obj_nargs = 0;

/* Companion registry for the exact symbolic gradient (FindMinimum). When a solve
 * registers its gradient-component array here, `fm_eval_gradient` evaluates each
 * component through its compiled program instead of the interpreter — the SAME
 * symbolic ∂f/∂x_i, just lowered, so the gradient stays exact (no finite
 * differences) while running on the register machine. Keyed by the g_exprs array
 * pointer + arity, so constraint gradients (a different array) stay on the
 * interpreter; a NULL or non-finite component falls back per-component. */
Expr**            g_fm_grad_exprs = NULL;
CompiledProgram** g_fm_grad_progs = NULL;   /* len g_fm_grad_n, entries may be NULL */
size_t            g_fm_grad_n     = 0;

/* Evaluate the bound objective and return a double; NULL on failure. */
bool fm_eval_scalar(Expr* f, FmVarBind* binds,
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
Expr** fm_compute_gradient(Expr* f, Expr** vars, size_t n) {
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
Expr*** fm_compute_hessian(Expr* f, Expr** vars, size_t n) {
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
bool fm_grad_finite_diff(Expr* f, FmVarBind* binds,
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
bool fm_eval_gradient(Expr** g_exprs, FmVarBind* binds,
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
bool fm_eval_hessian(Expr*** H_exprs, FmVarBind* binds,
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
Expr* fm_build_result(double fmin, Expr** vars, const double* vals,
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
Expr* fm_build_result_mpfr(const mpfr_t fmin, Expr** vars,
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
mpfr_t* fm_mpfr_array(size_t count, long bits) {
    mpfr_t* arr = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(arr[i], bits);
    return arr;
}

void fm_mpfr_array_free(mpfr_t* arr, size_t count) {
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
bool fm_run_bfgs_mpfr(Expr* f, Expr** vars, size_t n,
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

FmSpecKind fm_parse_var_spec(Expr* spec, Expr** var_out,
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
bool fm_constraint_to_g(Expr* cmp, Expr** expr_out, bool* equality_out) {
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
bool fm_bool_supported(Expr* c) {
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
bool fm_collect_constraints(Expr* cons, Expr** vars, size_t nvars,
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
void fm_project_box(double* x, size_t n, const FmBox* boxes) {
    for (size_t i = 0; i < n; i++) {
        if (boxes[i].has_lo && x[i] < boxes[i].lo) x[i] = boxes[i].lo;
        if (boxes[i].has_hi && x[i] > boxes[i].hi) x[i] = boxes[i].hi;
    }
}

/* Evaluate Σ max(0, g_i(x))^2 + Σ h_j(x)^2 over the general constraint
 * set. Returns false if any constraint cannot be evaluated. */
bool fm_eval_penalty(const FmGenCon* gens, size_t ngens,
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
bool fm_eval_augmented(Expr* f, FmVarBind* binds,
                              const double* x, size_t n,
                              const FmGenCon* gens, size_t ngens,
                              double mu, const FmOpts* opts, double* out) {
    double fv;
    if (!fm_eval_scalar(f, binds, x, n, opts, &fv)) return false;
    /* Augmented-Lagrangian branch: active only when multipliers are installed
     * for exactly this constraint set (see g_fm_al_lambda). Otherwise fall
     * through to the historical quadratic-penalty value, byte-for-byte. */
    if (g_fm_al_lambda && gens == g_fm_al_gens && ngens > 0 && mu > 0.0) {
        double aug = 0.0;
        for (size_t k = 0; k < ngens; k++) {
            double c;
            if (!fm_eval_scalar(gens[k].expr, binds, x, n, opts, &c)) return false;
            double lam = g_fm_al_lambda[k];
            if (gens[k].equality) {
                aug += lam * c + mu * c * c;
            } else {
                double s = lam + 2.0 * mu * c;
                if (s > 0.0) aug += lam * c + mu * c * c;
                else         aug += -lam * lam / (4.0 * mu);
            }
        }
        *out = fv + aug;
        return true;
    }
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
bool fm_eval_aug_gradient(Expr* f, Expr** g_exprs,
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

    bool use_al = (g_fm_al_lambda && gens == g_fm_al_gens);
    double* gk_grad = (double*)malloc(sizeof(double) * n);
    if (!gk_grad) return false;
    for (size_t k = 0; k < ngens; k++) {
        double gk;
        if (!fm_eval_scalar(gens[k].expr, binds, x, n, opts, &gk)) {
            free(gk_grad);
            return false;
        }
        /* Augmented-Lagrangian active value a_k, whose ∇ contribution is a_k·∇c_k.
         * With λ_k == 0 this is exactly the quadratic-penalty coefficient
         * 2μ·gk (equality) / 2μ·gk when gk>0 (inequality), so the non-AL path is
         * unchanged. Inequalities contribute only when s_k = λ_k + 2μ·gk > 0. */
        double lam = use_al ? g_fm_al_lambda[k] : 0.0;
        double a_k;
        if (gens[k].equality) {
            a_k = lam + 2.0 * mu * gk;
        } else {
            double s = lam + 2.0 * mu * gk;
            if (s <= 0.0) continue;
            a_k = s;
        }
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
        for (size_t i = 0; i < n; i++) g_out[i] += a_k * gk_grad[i];
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
bool fm_line_search(Expr* f, FmVarBind* binds, size_t n,
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
bool fm_bracket(Expr* f, FmVarBind* binds, const FmOpts* opts,
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
bool fm_bracket_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
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

bool fm_brent_min_mpfr(Expr* f, FmVarBind* bind, const FmOpts* opts,
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

bool fm_brent_min(Expr* f, FmVarBind* bind, const FmOpts* opts,
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
bool fm_chol_factor(double* H, size_t n, double tau) {
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

void fm_chol_solve(const double* L, size_t n,
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

double fm_dot(const double* a, const double* b, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* out = B·v for a dense symmetric B (row-major). */
void fm_matvec(const double* B, size_t n, const double* v, double* out) {
    for (size_t i = 0; i < n; i++) {
        double s = 0.0;
        const double* Bi = B + i * n;
        for (size_t j = 0; j < n; j++) s += Bi[j] * v[j];
        out[i] = s;
    }
}

/* out = B·v, dense if available else a finite-difference Hessian-vector
 * product on the compiled gradient. Returns false only on a non-finite Hv. */
bool fm_quad_matvec(const FmQuad* q, const double* v, double* out) {
    if (q->B) { fm_matvec(q->B, q->n, v, out); return true; }
    return fm_tnc_hessvec(q->c, q->xbase, q->gbase, v, q->active0,
                          q->xpert, q->gpert, out);
}

/* ------------------------------------------------------------------ *
 *  Nelder-Mead downhill simplex (derivative-free)                      *
 * ------------------------------------------------------------------ *
 * The Nelder & Mead (1965) simplex: n+1 vertices reflected, expanded,
 * contracted and shrunk toward the best, using function values only.
 * Robust on non-smooth / noisy objectives where the gradient methods (and
 * even Powell's line search) can struggle, at the cost of slower
 * convergence on smooth ones. This is the LOCAL simplex (a single start,
 * matching scipy's _minimize_neldermead) -- distinct from the restarted,
 * box-sampling nm_neldermead the global NMinimize driver uses. Standard
 * non-adaptive coefficients rho/chi/psi/sigma = 1/2/0.5/0.5 (scipy's
 * default), so the two agree on the minimiser. Machine precision only
 * (rides the compiled objective through fm_eval_scalar); box bounds by
 * projecting every candidate vertex; n==1 delegates to the exact Brent
 * path. Exposed as Method -> "NelderMead". */

/* Insertion-sort the vertex index permutation by ascending fsim. np = n+1
 * is small (Nelder-Mead is only practical for modest n), so an O(np^2)
 * sort of a size_t permutation — not the vertex rows — is the cheap choice. */
void fm_nm_sort_idx(size_t* idx, const double* fsim, size_t np) {
    for (size_t i = 1; i < np; i++) {
        size_t key = idx[i];
        double fk = fsim[key];
        size_t j = i;
        while (j > 0 && fsim[idx[j-1]] > fk) { idx[j] = idx[j-1]; j--; }
        idx[j] = key;
    }
}
