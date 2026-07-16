/* fit.c — linear least-squares regression (Fit) and DesignMatrix.
 *
 * This module implements Mathematica's `Fit` builtin: it fits a linear
 * combination  a1 f1 + ... + an fn  of basis functions to data, plus the
 * companion `DesignMatrix` (the matrix of basis functions evaluated at the
 * data coordinates).
 *
 * Call forms
 * ----------
 *   Fit[data, {f1,...,fn}, vars]
 *       Fits a1 f1 + ... + an fn to `data`.  `vars` is a single symbol `x`
 *       or a list {x, y, ...}.  Returns the symbolic fit expression.
 *   Fit[{m, v}]
 *       Given a design matrix `m` and response vector `v`, returns the
 *       coefficient vector a minimising ||m.a - v||.
 *   DesignMatrix[data, {f1,...,fn}, vars]
 *       Returns the design matrix m_ij = f_i(coords_j).
 *
 * Data shapes (3-argument form)
 * -----------------------------
 *   {v1,...,vn}            equivalent to {{1,v1},...,{n,vn}}.
 *   {{x1,v1},...}          univariate: coordinate x_i, response v_i.
 *   {{x1,...,xk,v1},...}   multivariate: leading k coordinates, last value.
 *
 * Options
 * -------
 *   WorkingPrecision -> Automatic  (default; exact input -> machine reals)
 *                     -> n         (n-digit MPFR arithmetic)
 *                     -> Infinity  (exact rational arithmetic)
 *   FitRegularization -> {"Tikhonov"|"L2"|"RidgeRegression", lambda}
 *                          minimise ||m.a-v||^2 + lambda ||a||^2 (ridge).
 *                     -> {"LASSO"|"L1", lambda}
 *                          minimise ||m.a-v||^2 + lambda ||a||_1.
 *   NormFunction -> Function[Norm[#,p]]
 *                     minimise normf[m.a - v] instead of the 2-norm.
 *
 * Solvers (reuse-first)
 * ---------------------
 *   * Plain L2 and ridge route through the existing LeastSquares builtin,
 *     which already supports exact (rational), machine (Real) and MPFR
 *     arithmetic.  Ridge is reduced to ordinary least squares on the
 *     augmented system [m; sqrt(lambda) I] / [v; 0].
 *   * LASSO uses cyclic coordinate descent with soft-thresholding
 *     (machine precision).
 *   * NormFunction -> Norm[#,1] (least absolute deviations) uses iteratively
 *     reweighted least squares (IRLS, machine precision).
 *   * Any other norm (or a norm combined with regularisation) falls back to
 *     the FindMinimum builtin, warm-started from the L2 solution.
 *
 * Memory ownership follows the standard builtin contract: the evaluator
 * owns `res` and frees it; this file never frees `res` or its argument
 * subtrees.  Every intermediate built here is freed on every return path.
 */

#include "fit.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "match.h"
#include "numeric.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================== *
 *  Diagnostics                                                        *
 * ================================================================== */

static void fit_warn(const char* head, const char* tag, const char* fmt, ...) {
    fprintf(stderr, "%s::%s: ", head, tag);
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* ================================================================== *
 *  Small structural helpers                                           *
 * ================================================================== */

static bool fit_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

static size_t fit_len(const Expr* e) {
    return (e && e->type == EXPR_FUNCTION) ? e->data.function.arg_count : 0;
}

static Expr* fit_elem(const Expr* e, size_t i) {
    return e->data.function.args[i];
}

static bool fit_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Infinity as a WorkingPrecision value: the bare symbol Infinity (before
 * evaluation) or DirectedInfinity[1] (after). */
static bool fit_is_infinity(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity &&
        e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        return a->type == EXPR_INTEGER && a->data.integer == 1;
    }
    return false;
}

/* Build List[args...] taking ownership of the element pointers. */
static Expr* fit_list(Expr** elems, size_t n) {
    Expr* l = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    return l;
}

/* Numeric atom -> double.  Handles Integer, Real, BigInt, MPFR and
 * Rational[n,d].  Returns false for anything non-numeric. */
static bool fit_to_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (fit_head_is(e, "Rational") && e->data.function.arg_count == 2) {
                double n, d;
                if (fit_to_double(e->data.function.args[0], &n) &&
                    fit_to_double(e->data.function.args[1], &d) && d != 0.0) {
                    *out = n / d;
                    return true;
                }
            }
            return false;
        default: return false;
    }
}

/* ================================================================== *
 *  Options                                                            *
 * ================================================================== */

typedef enum { FIT_WP_AUTO, FIT_WP_MACHINE, FIT_WP_MPFR, FIT_WP_INFINITY } FitWPMode;
typedef enum { FIT_REG_NONE, FIT_REG_TIKHONOV, FIT_REG_LASSO } FitRegKind;

typedef struct {
    FitWPMode  wp_mode;
    long       wp_bits;     /* when wp_mode == FIT_WP_MPFR */
    FitRegKind reg_kind;
    Expr*      reg_lambda;  /* borrowed; numeric scalar */
    Expr*      norm_fun;    /* borrowed; pure Function or NULL (=> 2-norm) */
} FitOpts;

static bool fit_parse_working_precision(Expr* v, FitWPMode* mode, long* bits) {
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) {
        *mode = FIT_WP_AUTO; *bits = 0; return true;
    }
    if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_MachinePrecision) {
        *mode = FIT_WP_MACHINE; *bits = 0; return true;
    }
    if (fit_is_infinity(v)) { *mode = FIT_WP_INFINITY; *bits = 0; return true; }

    double digits = 0.0;
    if (v->type == EXPR_INTEGER)   digits = (double)v->data.integer;
    else if (v->type == EXPR_REAL) digits = v->data.real;
    else if (!fit_to_double(v, &digits)) return false;
    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        *mode = FIT_WP_MACHINE; *bits = 0;
    } else {
        *mode = FIT_WP_MPFR; *bits = numeric_digits_to_bits(digits);
    }
#else
    *mode = FIT_WP_MACHINE; *bits = 0;
#endif
    return true;
}

static bool fit_parse_regularization(Expr* v, FitRegKind* kind, Expr** lambda) {
    /* Expect {"name", lambda}. */
    if (!fit_is_list(v) || fit_len(v) != 2) return false;
    Expr* name = fit_elem(v, 0);
    Expr* lam  = fit_elem(v, 1);
    if (name->type != EXPR_STRING) return false;
    const char* s = name->data.string;
    if (strcmp(s, "Tikhonov") == 0 || strcmp(s, "L2") == 0 ||
        strcmp(s, "RidgeRegression") == 0) {
        *kind = FIT_REG_TIKHONOV;
    } else if (strcmp(s, "LASSO") == 0 || strcmp(s, "L1") == 0) {
        *kind = FIT_REG_LASSO;
    } else {
        return false;
    }
    *lambda = lam;  /* borrowed */
    return true;
}

static bool fit_is_option_name(const char* s) {
    return s == SYM_WorkingPrecision ||
           s == SYM_FitRegularization ||
           s == SYM_NormFunction;
}

static bool fit_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    return lhs->type == EXPR_SYMBOL && fit_is_option_name(lhs->data.symbol.name);
}

/* Parse trailing options.  *npos receives the number of leading positional
 * arguments.  Returns false (and emits a diagnostic) on an invalid value. */
static bool fit_parse_options(Expr** args, size_t argc, size_t* npos, FitOpts* opts) {
    opts->wp_mode = FIT_WP_AUTO;
    opts->wp_bits = 0;
    opts->reg_kind = FIT_REG_NONE;
    opts->reg_lambda = NULL;
    opts->norm_fun = NULL;

    size_t n = argc;
    while (n > 0 && fit_is_option_arg(args[n - 1])) n--;
    *npos = n;

    for (size_t i = n; i < argc; i++) {
        Expr* rule = args[i];
        const char* name = rule->data.function.args[0]->data.symbol.name;
        Expr* val = rule->data.function.args[1];
        if (name == SYM_WorkingPrecision) {
            if (!fit_parse_working_precision(val, &opts->wp_mode, &opts->wp_bits)) {
                fit_warn("Fit", "wprec", "invalid WorkingPrecision value");
                return false;
            }
        } else if (name == SYM_FitRegularization) {
            if (!fit_parse_regularization(val, &opts->reg_kind, &opts->reg_lambda)) {
                fit_warn("Fit", "reg",
                         "FitRegularization must be {\"name\", lambda}");
                return false;
            }
        } else if (name == SYM_NormFunction) {
            opts->norm_fun = val;  /* borrowed */
        }
    }
    return true;
}

/* ================================================================== *
 *  Data / variable normalisation                                      *
 * ================================================================== */

static void fit_free_array(Expr** a, size_t n) {
    if (!a) return;
    for (size_t i = 0; i < n; i++) expr_free(a[i]);
    free(a);
}

/* Normalise `data` into a flat coordinate array (npts*ncoord, row-major)
 * and a response array (npts).  Each entry is an owned copy.  Returns
 * false (with a diagnostic) on malformed data. */
static bool fit_normalize_data(const Expr* data, const char* who,
                               size_t* npts_out, size_t* ncoord_out,
                               Expr*** coords_out, Expr*** resp_out) {
    if (!fit_is_list(data)) {
        fit_warn(who, "notdata", "first argument must be a list of data");
        return false;
    }
    size_t npts = fit_len(data);
    if (npts == 0) {
        fit_warn(who, "empty", "data must be non-empty");
        return false;
    }

    Expr* first = fit_elem(data, 0);
    Expr** coords = NULL;
    Expr** resp = NULL;
    size_t ncoord;

    if (fit_is_list(first)) {
        size_t rowlen = fit_len(first);
        if (rowlen < 2) {
            fit_warn(who, "shape", "each data row needs at least 2 entries");
            return false;
        }
        ncoord = rowlen - 1;
        coords = calloc(npts * ncoord, sizeof(Expr*));
        resp = calloc(npts, sizeof(Expr*));
        for (size_t r = 0; r < npts; r++) {
            Expr* row = fit_elem(data, r);
            if (!fit_is_list(row) || fit_len(row) != rowlen) {
                fit_warn(who, "shape", "data rows must all have the same length");
                fit_free_array(coords, npts * ncoord);
                fit_free_array(resp, npts);
                return false;
            }
            for (size_t k = 0; k < ncoord; k++)
                coords[r * ncoord + k] = expr_copy(fit_elem(row, k));
            resp[r] = expr_copy(fit_elem(row, rowlen - 1));
        }
    } else {
        /* Flat list of values: coordinates are 1..npts. */
        ncoord = 1;
        coords = calloc(npts, sizeof(Expr*));
        resp = calloc(npts, sizeof(Expr*));
        for (size_t r = 0; r < npts; r++) {
            Expr* el = fit_elem(data, r);
            if (fit_is_list(el)) {
                fit_warn(who, "shape", "data mixes scalars and lists");
                fit_free_array(coords, npts);
                fit_free_array(resp, npts);
                return false;
            }
            coords[r] = expr_new_integer((int64_t)(r + 1));
            resp[r] = expr_copy(el);
        }
    }

    *npts_out = npts;
    *ncoord_out = ncoord;
    *coords_out = coords;
    *resp_out = resp;
    return true;
}

/* Normalise `vars` (a symbol or a list of symbols) into an owned array.
 * Validates the count against ncoord. */
static bool fit_normalize_vars(Expr* vars, const char* who, size_t ncoord,
                               Expr*** vars_out, size_t* nvar_out) {
    Expr** vs = NULL;
    size_t nvar;
    if (fit_is_list(vars)) {
        nvar = fit_len(vars);
        vs = calloc(nvar ? nvar : 1, sizeof(Expr*));
        for (size_t k = 0; k < nvar; k++) {
            Expr* s = fit_elem(vars, k);
            if (s->type != EXPR_SYMBOL) {
                fit_warn(who, "var", "variables must be symbols");
                fit_free_array(vs, k);
                return false;
            }
            vs[k] = expr_copy(s);
        }
    } else if (vars->type == EXPR_SYMBOL) {
        nvar = 1;
        vs = calloc(1, sizeof(Expr*));
        vs[0] = expr_copy(vars);
    } else {
        fit_warn(who, "var", "variable specification must be a symbol or list");
        return false;
    }
    if (nvar != ncoord) {
        fit_warn(who, "vararity",
                 "number of variables (%zu) does not match data coordinates (%zu)",
                 nvar, ncoord);
        fit_free_array(vs, nvar);
        return false;
    }
    *vars_out = vs;
    *nvar_out = nvar;
    return true;
}

/* ================================================================== *
 *  Design matrix / response                                           *
 * ================================================================== */

/* Build the npts x nfun design matrix m_rj = funs[j](coords[r]).  Returns an
 * owned List of Lists, or NULL if a basis function fails to evaluate to a
 * usable value (it is left symbolic, which is allowed). */
static Expr* fit_build_design_matrix(Expr** funs, size_t nfun,
                                     Expr** coords, size_t npts, size_t ncoord,
                                     Expr** vars) {
    Expr** rows = calloc(npts, sizeof(Expr*));
    for (size_t r = 0; r < npts; r++) {
        MatchEnv* env = env_new();
        for (size_t k = 0; k < ncoord; k++)
            env_set(env, vars[k]->data.symbol.name, coords[r * ncoord + k]);
        Expr** cells = calloc(nfun, sizeof(Expr*));
        for (size_t j = 0; j < nfun; j++) {
            Expr* sub = replace_bindings(funs[j], env);
            cells[j] = eval_and_free(sub);
        }
        rows[r] = fit_list(cells, nfun);
        free(cells);
        env_free(env);
    }
    Expr* m = fit_list(rows, npts);
    free(rows);
    return m;
}

static Expr* fit_build_response(Expr** resp, size_t npts) {
    Expr** cells = calloc(npts, sizeof(Expr*));
    for (size_t r = 0; r < npts; r++) cells[r] = expr_copy(resp[r]);
    Expr* v = fit_list(cells, npts);
    free(cells);
    return v;
}

/* ================================================================== *
 *  Precision realisation                                              *
 * ================================================================== */

#ifdef USE_MPFR
/* Largest MPFR precision (in bits) appearing in `e`, or 0 if none. */
static long fit_max_mpfr_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
    if (e->type == EXPR_FUNCTION) {
        long m = fit_max_mpfr_bits(e->data.function.head);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            long b = fit_max_mpfr_bits(e->data.function.args[i]);
            if (b > m) m = b;
        }
        return m;
    }
    return 0;
}
#endif

static Expr* fit_realize_precision(Expr* e, const FitOpts* opts) {
    switch (opts->wp_mode) {
        case FIT_WP_INFINITY:
            return expr_copy(e);
#ifdef USE_MPFR
        case FIT_WP_MPFR: {
            NumericSpec spec;
            spec.mode = NUMERIC_MODE_MPFR;
            spec.bits = opts->wp_bits;
            return numericalize(e, spec);
        }
        case FIT_WP_AUTO: {
            /* Exact input -> machine, but already-approximate MPFR input keeps
             * its own precision (WorkingPrecision -> Automatic semantics). */
            long bits = fit_max_mpfr_bits(e);
            if (bits > 0) {
                NumericSpec spec;
                spec.mode = NUMERIC_MODE_MPFR;
                spec.bits = bits;
                return numericalize(e, spec);
            }
            return numericalize(e, numeric_machine_spec());
        }
#endif
        case FIT_WP_MACHINE:
        default:
            return numericalize(e, numeric_machine_spec());
    }
}

/* ================================================================== *
 *  Solvers                                                            *
 * ================================================================== */

/* Call LeastSquares[m, v]; returns owned coefficient vector or NULL if the
 * call comes back unevaluated. */
static Expr* fit_leastsquares(Expr* m, Expr* v) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_LeastSquares),
                                   (Expr*[]){ expr_copy(m), expr_copy(v) }, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    if (!r) return NULL;
    if (fit_head_is(r, "LeastSquares")) { expr_free(r); return NULL; }
    return r;
}

/* Plain L2: realise precision then LeastSquares. */
static Expr* fit_solve_l2(Expr* M, Expr* V, const FitOpts* opts) {
    Expr* mR = fit_realize_precision(M, opts);
    Expr* vR = fit_realize_precision(V, opts);
    Expr* coeffs = fit_leastsquares(mR, vR);
    expr_free(mR);
    expr_free(vR);
    return coeffs;
}

/* Ridge / Tikhonov: augment with sqrt(lambda) I and zeros, then L2. */
static Expr* fit_solve_tikhonov(Expr* M, Expr* V, const FitOpts* opts) {
    size_t npts = fit_len(M);
    if (npts == 0) return NULL;
    size_t n = fit_len(fit_elem(M, 0));  /* columns = basis functions */

    Expr* sqrtlam = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Sqrt), (Expr*[]){ expr_copy(opts->reg_lambda) }, 1));

    size_t total = npts + n;
    Expr** rows = calloc(total, sizeof(Expr*));
    for (size_t r = 0; r < npts; r++) rows[r] = expr_copy(fit_elem(M, r));
    for (size_t i = 0; i < n; i++) {
        Expr** cells = calloc(n, sizeof(Expr*));
        for (size_t c = 0; c < n; c++)
            cells[c] = (c == i) ? expr_copy(sqrtlam) : expr_new_integer(0);
        rows[npts + i] = fit_list(cells, n);
        free(cells);
    }
    Expr* mp = fit_list(rows, total);
    free(rows);
    expr_free(sqrtlam);

    Expr** vcells = calloc(total, sizeof(Expr*));
    for (size_t r = 0; r < npts; r++) vcells[r] = expr_copy(fit_elem(V, r));
    for (size_t i = 0; i < n; i++) vcells[npts + i] = expr_new_integer(0);
    Expr* vp = fit_list(vcells, total);
    free(vcells);

    Expr* mpR = fit_realize_precision(mp, opts);
    Expr* vpR = fit_realize_precision(vp, opts);
    expr_free(mp);
    expr_free(vp);
    Expr* coeffs = fit_leastsquares(mpR, vpR);
    expr_free(mpR);
    expr_free(vpR);
    return coeffs;
}

/* ---- double-precision extraction + dense solve (for LASSO/IRLS) ---- */

static double* fit_matrix_doubles(const Expr* m, size_t* rows, size_t* cols) {
    if (!fit_is_list(m) || fit_len(m) == 0) return NULL;
    size_t R = fit_len(m);
    size_t C = fit_len(fit_elem(m, 0));
    if (C == 0) return NULL;
    double* A = malloc(sizeof(double) * R * C);
    for (size_t i = 0; i < R; i++) {
        Expr* row = fit_elem(m, i);
        if (!fit_is_list(row) || fit_len(row) != C) { free(A); return NULL; }
        for (size_t j = 0; j < C; j++) {
            if (!fit_to_double(fit_elem(row, j), &A[i * C + j])) { free(A); return NULL; }
        }
    }
    *rows = R; *cols = C;
    return A;
}

static double* fit_vector_doubles(const Expr* v, size_t* n) {
    if (!fit_is_list(v)) return NULL;
    size_t N = fit_len(v);
    double* b = malloc(sizeof(double) * (N ? N : 1));
    for (size_t i = 0; i < N; i++) {
        if (!fit_to_double(fit_elem(v, i), &b[i])) { free(b); return NULL; }
    }
    *n = N;
    return b;
}

/* Solve the n x n system A x = b by Gaussian elimination with partial
 * pivoting.  A (row-major) and b are overwritten.  Returns false if the
 * matrix is singular. */
static bool fit_solve_dense(size_t n, double* A, double* b, double* x) {
    for (size_t col = 0; col < n; col++) {
        size_t piv = col;
        double best = fabs(A[col * n + col]);
        for (size_t r = col + 1; r < n; r++) {
            double v = fabs(A[r * n + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-300) return false;
        if (piv != col) {
            for (size_t c = 0; c < n; c++) {
                double t = A[col * n + c]; A[col * n + c] = A[piv * n + c]; A[piv * n + c] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        double diag = A[col * n + col];
        for (size_t r = col + 1; r < n; r++) {
            double f = A[r * n + col] / diag;
            if (f == 0.0) continue;
            for (size_t c = col; c < n; c++) A[r * n + c] -= f * A[col * n + c];
            b[r] -= f * b[col];
        }
    }
    for (size_t ii = n; ii-- > 0; ) {
        double s = b[ii];
        for (size_t c = ii + 1; c < n; c++) s -= A[ii * n + c] * x[c];
        x[ii] = s / A[ii * n + ii];
    }
    return true;
}

/* Build a {r1, r2, ...} List of machine reals from a double array. */
static Expr* fit_doubles_to_vector(const double* x, size_t n) {
    Expr** cells = calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) cells[i] = expr_new_real(x[i]);
    Expr* v = fit_list(cells, n);
    free(cells);
    return v;
}

/* LASSO via cyclic coordinate descent (machine precision).  Minimises
 * ||A a - b||^2 + lambda ||a||_1.  The subgradient of the smooth part is
 * 2 A^T (A a - b), so the soft-threshold level is lambda/2. */
static Expr* fit_solve_lasso(Expr* M, Expr* V, const FitOpts* opts) {
    if (opts->wp_mode == FIT_WP_INFINITY || opts->wp_mode == FIT_WP_MPFR) {
        fit_warn("Fit", "lassonum",
                 "LASSO requires machine WorkingPrecision");
        return NULL;
    }
    double lambda;
    if (!fit_to_double(opts->reg_lambda, &lambda) || lambda < 0.0) {
        fit_warn("Fit", "reglam", "regularization parameter must be non-negative");
        return NULL;
    }

    Expr* mR = numericalize(M, numeric_machine_spec());
    Expr* vR = numericalize(V, numeric_machine_spec());
    size_t npts = 0, n = 0, nb = 0;
    double* A = fit_matrix_doubles(mR, &npts, &n);
    double* b = fit_vector_doubles(vR, &nb);
    expr_free(mR);
    expr_free(vR);
    if (!A || !b || nb != npts) {
        free(A); free(b);
        fit_warn("Fit", "numdata", "data must be numeric for this method");
        return NULL;
    }

    double* a = calloc(n, sizeof(double));
    double* col_sq = calloc(n, sizeof(double));   /* ||A_j||^2 */
    double* r = malloc(sizeof(double) * npts);    /* residual A a - b */
    for (size_t j = 0; j < n; j++) {
        double s = 0.0;
        for (size_t i = 0; i < npts; i++) s += A[i * n + j] * A[i * n + j];
        col_sq[j] = s;
    }
    for (size_t i = 0; i < npts; i++) r[i] = -b[i];  /* a = 0 */

    const int MAXIT = 5000;
    for (int it = 0; it < MAXIT; it++) {
        double max_delta = 0.0, max_abs = 0.0;
        for (size_t j = 0; j < n; j++) {
            if (col_sq[j] < 1e-300) continue;
            /* rho = A_j^T (b - A a_{-j}) = a_j col_sq - A_j^T r */
            double dot = 0.0;
            for (size_t i = 0; i < npts; i++) dot += A[i * n + j] * r[i];
            double rho = a[j] * col_sq[j] - dot;
            double thr = lambda / 2.0;
            double a_new;
            if (rho > thr)       a_new = (rho - thr) / col_sq[j];
            else if (rho < -thr) a_new = (rho + thr) / col_sq[j];
            else                 a_new = 0.0;
            double delta = a_new - a[j];
            if (delta != 0.0) {
                for (size_t i = 0; i < npts; i++) r[i] += A[i * n + j] * delta;
                a[j] = a_new;
            }
            if (fabs(delta) > max_delta) max_delta = fabs(delta);
            if (fabs(a[j]) > max_abs) max_abs = fabs(a[j]);
        }
        if (max_delta < 1e-12 * (1.0 + max_abs)) break;
    }

    Expr* result = fit_doubles_to_vector(a, n);
    free(A); free(b); free(a); free(col_sq); free(r);
    return result;
}

/* Least absolute deviations via iteratively reweighted least squares
 * (machine precision).  Minimises sum_i |A a - b|_i. */
static Expr* fit_solve_irls_l1(Expr* M, Expr* V, const FitOpts* opts) {
    if (opts->wp_mode == FIT_WP_INFINITY || opts->wp_mode == FIT_WP_MPFR) {
        fit_warn("Fit", "l1num",
                 "the 1-norm fit requires machine WorkingPrecision");
        return NULL;
    }
    Expr* mR = numericalize(M, numeric_machine_spec());
    Expr* vR = numericalize(V, numeric_machine_spec());
    size_t npts = 0, n = 0, nb = 0;
    double* A = fit_matrix_doubles(mR, &npts, &n);
    double* b = fit_vector_doubles(vR, &nb);
    expr_free(mR);
    expr_free(vR);
    if (!A || !b || nb != npts || npts < n) {
        free(A); free(b);
        fit_warn("Fit", "numdata", "data must be numeric for this method");
        return NULL;
    }

    double* a = calloc(n, sizeof(double));
    double* w = malloc(sizeof(double) * npts);
    double* ATA = malloc(sizeof(double) * n * n);
    double* ATb = malloc(sizeof(double) * n);
    double* anew = malloc(sizeof(double) * n);
    for (size_t i = 0; i < npts; i++) w[i] = 1.0;  /* first pass: ordinary LS */

    const int MAXIT = 200;
    bool ok = true;
    for (int it = 0; it < MAXIT; it++) {
        /* Weighted normal equations: (A^T W A) a = A^T W b. */
        for (size_t p = 0; p < n; p++) {
            for (size_t q = 0; q < n; q++) {
                double s = 0.0;
                for (size_t i = 0; i < npts; i++) s += w[i] * A[i * n + p] * A[i * n + q];
                ATA[p * n + q] = s;
            }
            double s = 0.0;
            for (size_t i = 0; i < npts; i++) s += w[i] * A[i * n + p] * b[i];
            ATb[p] = s;
        }
        if (!fit_solve_dense(n, ATA, ATb, anew)) { ok = false; break; }

        double max_delta = 0.0, max_abs = 0.0;
        for (size_t j = 0; j < n; j++) {
            double d = fabs(anew[j] - a[j]);
            if (d > max_delta) max_delta = d;
            if (fabs(anew[j]) > max_abs) max_abs = fabs(anew[j]);
            a[j] = anew[j];
        }
        /* Update weights from the new residuals: w_i = 1/max(|r_i|, eps). */
        for (size_t i = 0; i < npts; i++) {
            double ri = -b[i];
            for (size_t j = 0; j < n; j++) ri += A[i * n + j] * a[j];
            double ar = fabs(ri);
            w[i] = 1.0 / (ar > 1e-8 ? ar : 1e-8);
        }
        if (it > 0 && max_delta < 1e-10 * (1.0 + max_abs)) break;
    }

    Expr* result = ok ? fit_doubles_to_vector(a, n) : NULL;
    if (!ok) fit_warn("Fit", "l1sing", "1-norm fit failed (singular system)");
    free(A); free(b); free(a); free(w); free(ATA); free(ATb); free(anew);
    return result;
}

/* General NormFunction / combined objective via FindMinimum, warm-started
 * from the L2 solution.  Builds fresh coefficient symbols c$1..c$n. */
static Expr* fit_solve_findmin(Expr* M, Expr* V, const FitOpts* opts) {
    size_t npts = fit_len(M);
    if (npts == 0) return NULL;
    size_t n = fit_len(fit_elem(M, 0));

    /* Warm start: machine L2 solution. */
    Expr* mMach = numericalize(M, numeric_machine_spec());
    Expr* vMach = numericalize(V, numeric_machine_spec());
    Expr* l2 = fit_leastsquares(mMach, vMach);
    expr_free(mMach);
    expr_free(vMach);
    double* x0 = calloc(n, sizeof(double));
    if (l2 && fit_is_list(l2) && fit_len(l2) == n) {
        for (size_t j = 0; j < n; j++) fit_to_double(fit_elem(l2, j), &x0[j]);
    }
    if (l2) expr_free(l2);

    /* Coefficient symbol vector a = {c$1, ..., c$n}. */
    Expr** csyms = calloc(n, sizeof(Expr*));
    Expr** avec = calloc(n, sizeof(Expr*));
    for (size_t j = 0; j < n; j++) {
        char buf[32];
        snprintf(buf, sizeof(buf), "c$%zu", j + 1);
        csyms[j] = expr_new_symbol(buf);
        avec[j] = expr_new_symbol(buf);
    }
    Expr* aexpr = fit_list(avec, n);
    free(avec);

    /* Residual = M . a - V  (realised to the requested precision). */
    Expr* mR = fit_realize_precision(M, opts);
    Expr* vR = fit_realize_precision(V, opts);
    Expr* dotma = expr_new_function(expr_new_symbol(SYM_Dot),
                                    (Expr*[]){ mR, aexpr }, 2);
    Expr* negv = expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), vR }, 2);
    Expr* residual = expr_new_function(expr_new_symbol(SYM_Plus),
                                       (Expr*[]){ dotma, negv }, 2);

    /* Objective = normf[residual] + reg(a). */
    Expr* normf = opts->norm_fun
        ? expr_new_function(expr_copy(opts->norm_fun),
                            (Expr*[]){ residual }, 1)
        : expr_new_function(expr_new_symbol(SYM_Norm),
                            (Expr*[]){ residual }, 1);
    Expr* objective = normf;
    if (opts->reg_kind == FIT_REG_TIKHONOV) {
        /* + lambda * (a.a) */
        Expr* aa = expr_new_function(expr_new_symbol(SYM_Dot),
            (Expr*[]){ expr_copy(aexpr), expr_copy(aexpr) }, 2);
        Expr* term = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(opts->reg_lambda), aa }, 2);
        objective = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ normf, term }, 2);
    } else if (opts->reg_kind == FIT_REG_LASSO) {
        Expr* anorm = expr_new_function(expr_new_symbol(SYM_Norm),
            (Expr*[]){ expr_copy(aexpr), expr_new_integer(1) }, 2);
        Expr* term = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(opts->reg_lambda), anorm }, 2);
        objective = expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ normf, term }, 2);
    }

    /* Variable spec {{c$1, x0_1}, ...}. */
    Expr** boxes = calloc(n, sizeof(Expr*));
    for (size_t j = 0; j < n; j++) {
        Expr** pair = calloc(2, sizeof(Expr*));
        pair[0] = expr_copy(csyms[j]);
        pair[1] = expr_new_real(x0[j]);
        boxes[j] = fit_list(pair, 2);
        free(pair);
    }
    Expr* varspec = fit_list(boxes, n);
    free(boxes);

    Expr* call = expr_new_function(expr_new_symbol(SYM_FindMinimum),
                                   (Expr*[]){ objective, varspec }, 2);
    Expr* fm = evaluate(call);
    expr_free(call);

    /* Result is {fmin, {c$j -> val, ...}}; extract values in symbol order. */
    Expr* coeffs = NULL;
    if (fm && fit_is_list(fm) && fit_len(fm) == 2 && fit_is_list(fit_elem(fm, 1))) {
        Expr* rules = fit_elem(fm, 1);
        Expr** cells = calloc(n, sizeof(Expr*));
        bool good = (fit_len(rules) == n);
        for (size_t j = 0; j < n && good; j++) {
            /* Find the rule whose LHS matches c$j. */
            Expr* val = NULL;
            for (size_t r = 0; r < fit_len(rules); r++) {
                Expr* rule = fit_elem(rules, r);
                if (rule->type == EXPR_FUNCTION && rule->data.function.arg_count == 2) {
                    Expr* lhs = rule->data.function.args[0];
                    if (lhs->type == EXPR_SYMBOL &&
                        lhs->data.symbol.name == csyms[j]->data.symbol.name) {
                        val = rule->data.function.args[1];
                        break;
                    }
                }
            }
            if (val) cells[j] = expr_copy(val);
            else good = false;
        }
        if (good) { coeffs = fit_list(cells, n); }
        else { for (size_t j = 0; j < n; j++) expr_free(cells[j]); free(cells); }
    }
    if (!coeffs) fit_warn("Fit", "nmin", "objective minimisation failed");
    if (fm) expr_free(fm);
    expr_free(aexpr);
    fit_free_array(csyms, n);
    free(x0);
    return coeffs;
}

/* ---- norm classification ---- */

typedef enum { FIT_NORM_L2, FIT_NORM_L1, FIT_NORM_OTHER } FitNormKind;

static FitNormKind fit_classify_norm(const Expr* fn) {
    if (!fn || fn->type != EXPR_FUNCTION ||
        fn->data.function.head->type != EXPR_SYMBOL ||
        fn->data.function.head->data.symbol.name != SYM_Function ||
        fn->data.function.arg_count < 1) {
        return FIT_NORM_OTHER;
    }
    Expr* body = fn->data.function.args[fn->data.function.arg_count - 1];
    if (!fit_head_is(body, "Norm")) return FIT_NORM_OTHER;
    size_t na = body->data.function.arg_count;
    if (na == 1) return FIT_NORM_L2;               /* Norm[#] is the 2-norm */
    if (na == 2) {
        Expr* p = body->data.function.args[1];
        if (p->type == EXPR_INTEGER) {
            if (p->data.integer == 2) return FIT_NORM_L2;
            if (p->data.integer == 1) return FIT_NORM_L1;
        }
        double pv;
        if (fit_to_double(p, &pv)) {
            if (pv == 2.0) return FIT_NORM_L2;
            if (pv == 1.0) return FIT_NORM_L1;
        }
        return FIT_NORM_OTHER;
    }
    return FIT_NORM_OTHER;
}

/* Top-level solver dispatch.  Returns an owned coefficient vector or NULL. */
static Expr* fit_dispatch_solve(Expr* M, Expr* V, const FitOpts* opts) {
    FitNormKind nk = opts->norm_fun ? fit_classify_norm(opts->norm_fun) : FIT_NORM_L2;

    if (opts->reg_kind == FIT_REG_LASSO) {
        if (opts->norm_fun && nk != FIT_NORM_L2)
            return fit_solve_findmin(M, V, opts);
        return fit_solve_lasso(M, V, opts);
    }
    if (nk == FIT_NORM_L1 && opts->reg_kind == FIT_REG_NONE)
        return fit_solve_irls_l1(M, V, opts);
    if (nk == FIT_NORM_OTHER || (nk == FIT_NORM_L1 && opts->reg_kind != FIT_REG_NONE))
        return fit_solve_findmin(M, V, opts);
    /* nk == L2 (or no norm): plain L2 or ridge. */
    if (opts->reg_kind == FIT_REG_TIKHONOV)
        return fit_solve_tikhonov(M, V, opts);
    return fit_solve_l2(M, V, opts);
}

/* ================================================================== *
 *  Output assembly                                                    *
 * ================================================================== */

static Expr* fit_assemble_output(const Expr* coeffs, Expr** funs, size_t nfun) {
    if (!fit_is_list(coeffs) || fit_len(coeffs) != nfun) return NULL;
    Expr** terms = calloc(nfun, sizeof(Expr*));
    for (size_t j = 0; j < nfun; j++) {
        terms[j] = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_copy(fit_elem(coeffs, j)), expr_copy(funs[j]) }, 2);
    }
    Expr* plus = expr_new_function(expr_new_symbol(SYM_Plus), terms, nfun);
    free(terms);
    return eval_and_free(plus);
}

/* ================================================================== *
 *  Builtins                                                           *
 * ================================================================== */

/* Fit[{m, v}] form: m a design matrix, v a response vector. */
static Expr* fit_from_design(Expr* arg0, const FitOpts* opts) {
    if (!fit_is_list(arg0) || fit_len(arg0) != 2) {
        fit_warn("Fit", "mv", "Fit[{m, v}] expects a matrix and a vector");
        return NULL;
    }
    Expr* m = fit_elem(arg0, 0);
    Expr* v = fit_elem(arg0, 1);
    if (!fit_is_list(m) || fit_len(m) == 0 || !fit_is_list(fit_elem(m, 0)) ||
        !fit_is_list(v)) {
        fit_warn("Fit", "mv", "Fit[{m, v}] expects a matrix and a vector");
        return NULL;
    }
    return fit_dispatch_solve(m, v, opts);
}

Expr* builtin_fit(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;
    Expr** args = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    FitOpts opts;
    size_t npos;
    if (!fit_parse_options(args, argc, &npos, &opts)) return NULL;

    if (npos == 1) {
        return fit_from_design(args[0], &opts);
    }
    if (npos != 3) {
        return NULL;  /* leave unevaluated */
    }

    Expr* data = args[0];
    Expr* funs_list = args[1];
    Expr* vars = args[2];
    if (!fit_is_list(funs_list) || fit_len(funs_list) == 0) {
        fit_warn("Fit", "funs", "second argument must be a non-empty list of functions");
        return NULL;
    }

    size_t npts = 0, ncoord = 0, nvar = 0;
    Expr** coords = NULL;
    Expr** resp = NULL;
    Expr** varsyms = NULL;
    Expr* M = NULL;
    Expr* V = NULL;
    Expr* coeffs = NULL;
    Expr* out = NULL;

    if (!fit_normalize_data(data, "Fit", &npts, &ncoord, &coords, &resp))
        return NULL;
    if (!fit_normalize_vars(vars, "Fit", ncoord, &varsyms, &nvar))
        goto cleanup;

    size_t nfun = fit_len(funs_list);
    M = fit_build_design_matrix(funs_list->data.function.args, nfun,
                                coords, npts, ncoord, varsyms);
    V = fit_build_response(resp, npts);

    coeffs = fit_dispatch_solve(M, V, &opts);
    if (!coeffs) goto cleanup;

    out = fit_assemble_output(coeffs, funs_list->data.function.args, nfun);

cleanup:
    fit_free_array(coords, npts * ncoord);
    fit_free_array(resp, npts);
    fit_free_array(varsyms, nvar);
    if (M) expr_free(M);
    if (V) expr_free(V);
    if (coeffs) expr_free(coeffs);
    return out;
}

Expr* builtin_designmatrix(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    Expr** args = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    FitOpts opts;
    size_t npos;
    if (!fit_parse_options(args, argc, &npos, &opts)) return NULL;
    if (npos != 3) return NULL;

    Expr* data = args[0];
    Expr* funs_list = args[1];
    Expr* vars = args[2];
    if (!fit_is_list(funs_list) || fit_len(funs_list) == 0) {
        fit_warn("DesignMatrix", "funs",
                 "second argument must be a non-empty list of functions");
        return NULL;
    }

    size_t npts = 0, ncoord = 0, nvar = 0;
    Expr** coords = NULL;
    Expr** resp = NULL;
    Expr** varsyms = NULL;
    Expr* M = NULL;
    Expr* out = NULL;

    if (!fit_normalize_data(data, "DesignMatrix", &npts, &ncoord, &coords, &resp))
        return NULL;
    if (!fit_normalize_vars(vars, "DesignMatrix", ncoord, &varsyms, &nvar))
        goto cleanup;

    M = fit_build_design_matrix(funs_list->data.function.args, fit_len(funs_list),
                                coords, npts, ncoord, varsyms);
    /* Apply WorkingPrecision if explicitly requested; otherwise keep exact. */
    if (opts.wp_mode == FIT_WP_MACHINE || opts.wp_mode == FIT_WP_MPFR) {
        out = fit_realize_precision(M, &opts);
    } else {
        out = M;
        M = NULL;  /* transfer ownership */
    }

cleanup:
    fit_free_array(coords, npts * ncoord);
    fit_free_array(resp, npts);
    fit_free_array(varsyms, nvar);
    if (M) expr_free(M);
    return out;
}

/* ================================================================== *
 *  Initialisation                                                     *
 * ================================================================== */

void fit_init(void) {
    symtab_add_builtin("Fit", builtin_fit);
    symtab_get_def("Fit")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("DesignMatrix", builtin_designmatrix);
    symtab_get_def("DesignMatrix")->attributes |= ATTR_PROTECTED;
}
