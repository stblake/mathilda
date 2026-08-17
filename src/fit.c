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
#include "ndarray.h"

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
 * false (with a diagnostic) on malformed data.  Operates on a List; the
 * NDArray surface is handled by the fit_normalize_data wrapper below. */
static bool fit_normalize_data_impl(const Expr* data, const char* who,
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

/* Normalise `data`, accepting the List surface directly and a packed
 * EXPR_NDARRAY by delisting it once.  (The numeric fast path reads NDArray
 * buffers without this delist; this is the symbolic-path / DesignMatrix
 * safety net that also fixes NDArray input, which previously errored.) */
static bool fit_normalize_data(const Expr* data, const char* who,
                               size_t* npts_out, size_t* ncoord_out,
                               Expr*** coords_out, Expr*** resp_out) {
    if (is_ndarray(data)) {
        Expr* lst = ndarray_to_nested_list(data);
        bool r = fit_normalize_data_impl(lst, who, npts_out, ncoord_out,
                                         coords_out, resp_out);
        expr_free(lst);
        return r;
    }
    return fit_normalize_data_impl(data, who, npts_out, ncoord_out,
                                   coords_out, resp_out);
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

/* LASSO coordinate-descent core, operating directly on raw double buffers
 * (row-major A is npts x n, b is length npts).  Minimises
 * ||A a - b||^2 + lambda ||a||_1.  The subgradient of the smooth part is
 * 2 A^T (A a - b), so the soft-threshold level is lambda/2.  Returns a freshly
 * malloc'd coefficient array of length n (never NULL); the buffers A and b are
 * read only, not freed.  Shared by the boxed wrapper and the packed fast path
 * so neither rebuilds the matrix the other already has. */
static double* fit_lasso_dense(const double* A, const double* b,
                               size_t npts, size_t n, double lambda) {
    double* a = calloc(n, sizeof(double));
    double* col_sq = calloc(n, sizeof(double));   /* ||A_j||^2 */
    double* r = malloc(sizeof(double) * (npts ? npts : 1));  /* residual A a - b */
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
    free(col_sq); free(r);
    return a;
}

/* LASSO via cyclic coordinate descent (machine precision).  Boxed-input
 * wrapper: extracts the double buffers then defers to fit_lasso_dense. */
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

    double* a = fit_lasso_dense(A, b, npts, n, lambda);
    Expr* result = fit_doubles_to_vector(a, n);
    free(A); free(b); free(a);
    return result;
}

/* Least-absolute-deviations IRLS core on raw double buffers (row-major A is
 * npts x n, b is length npts).  Minimises sum_i |A a - b|_i.  Returns a freshly
 * malloc'd coefficient array of length n, or NULL on a singular weighted system
 * (*ok set accordingly).  A and b are read only, not freed.  Shared by the
 * boxed wrapper and the packed fast path. */
static double* fit_irls_l1_dense(const double* A, const double* b,
                                 size_t npts, size_t n, bool* ok) {
    double* a = calloc(n, sizeof(double));
    double* w = malloc(sizeof(double) * (npts ? npts : 1));
    double* ATA = malloc(sizeof(double) * n * n);
    double* ATb = malloc(sizeof(double) * n);
    double* anew = malloc(sizeof(double) * n);
    for (size_t i = 0; i < npts; i++) w[i] = 1.0;  /* first pass: ordinary LS */

    const int MAXIT = 200;
    *ok = true;
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
        if (!fit_solve_dense(n, ATA, ATb, anew)) { *ok = false; break; }

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
    free(w); free(ATA); free(ATb); free(anew);
    if (!*ok) { free(a); return NULL; }
    return a;
}

/* Least absolute deviations via IRLS (machine precision).  Boxed-input
 * wrapper: extracts the double buffers then defers to fit_irls_l1_dense. */
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

    bool ok = true;
    double* a = fit_irls_l1_dense(A, b, npts, n, &ok);
    Expr* result = ok ? fit_doubles_to_vector(a, n) : NULL;
    if (!ok) fit_warn("Fit", "l1sing", "1-norm fit failed (singular system)");
    free(A); free(b); free(a);
    return result;
}

/* Design-matrix dimensions, accepting both a boxed List-of-Lists and a packed
 * rank-2 NDArray (the shape the machine fast path hands the solvers). */
static bool fit_matrix_dims(const Expr* M, size_t* rows, size_t* cols) {
    if (is_ndarray(M) && M->data.ndarray.rank == 2) {
        *rows = (size_t)M->data.ndarray.dims[0];
        *cols = (size_t)M->data.ndarray.dims[1];
        return true;
    }
    if (fit_is_list(M) && fit_len(M) > 0 && fit_is_list(fit_elem(M, 0))) {
        *rows = fit_len(M);
        *cols = fit_len(fit_elem(M, 0));
        return true;
    }
    return false;
}

/* General NormFunction / combined objective via FindMinimum, warm-started
 * from the L2 solution.  Builds fresh coefficient symbols c$1..c$n.  M may be a
 * boxed List-of-Lists or a packed rank-2 NDArray. */
static Expr* fit_solve_findmin(Expr* M, Expr* V, const FitOpts* opts) {
    size_t npts = 0, n = 0;
    if (!fit_matrix_dims(M, &npts, &n) || npts == 0) return NULL;

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
 *  Numeric machine-precision fast path                                *
 * ================================================================== *
 *
 * When the data is fully numeric and WorkingPrecision resolves to machine
 * precision, the design matrix is built column-wise — each basis function is
 * evaluated ONCE with the variables bound to the packed coordinate *vectors*,
 * reusing the packed elementwise kernels (Power/Sin/... ) exactly the way
 * numpy's column_stack does — and assembled into a genuinely packed float64
 * NDArray.  That trips the LAPACK fast path inside LeastSquares (and hands the
 * L1/LASSO/FindMinimum solvers their double buffer directly), turning the old
 * npts*nfun per-cell evaluate() + non-packed symbolic PseudoInverse solve into
 * nfun vectorized evaluations + one LAPACK gesdd.  The boxed per-cell path
 * (below) remains the fallback for exact/MPFR/symbolic fits and for any basis
 * column that does not reduce to numbers.
 */

/* Read numeric data into column-major coordinate columns (ncoord*npts) and a
 * response buffer (npts).  Accepts a plain List (flat values or {..,v} rows) or
 * a rank-1/rank-2 real EXPR_NDARRAY.  Returns false (freeing nothing it did not
 * allocate) if the data is not all-numeric or not a supported shape, so the
 * caller falls back to the symbolic path. */
static bool fit_read_numeric(const Expr* data, size_t* npts_out, size_t* ncoord_out,
                             double** coords_out, double** resp_out) {
    size_t npts, ncoord;
    double* coords;
    double* resp;

    if (is_ndarray(data)) {
        const NDArrayData* nd = &data->data.ndarray;
        if (ndt_is_complex(nd->dtype)) return false;
        if (nd->rank == 1) {
            npts = (size_t)nd->dims[0];
            if (npts == 0) return false;
            ncoord = 1;
            coords = malloc(sizeof(double) * npts);
            resp   = malloc(sizeof(double) * npts);
            for (size_t r = 0; r < npts; r++) {
                double re, im; ndt_get(nd->data, r, nd->dtype, &re, &im);
                coords[r] = (double)(r + 1);
                resp[r] = re;
            }
        } else if (nd->rank == 2) {
            npts = (size_t)nd->dims[0];
            size_t rowlen = (size_t)nd->dims[1];
            if (npts == 0 || rowlen < 2) return false;
            ncoord = rowlen - 1;
            coords = malloc(sizeof(double) * npts * ncoord);
            resp   = malloc(sizeof(double) * npts);
            for (size_t r = 0; r < npts; r++) {
                for (size_t k = 0; k < ncoord; k++) {
                    double re, im; ndt_get(nd->data, r * rowlen + k, nd->dtype, &re, &im);
                    coords[k * npts + r] = re;
                }
                double re, im; ndt_get(nd->data, r * rowlen + (rowlen - 1), nd->dtype, &re, &im);
                resp[r] = re;
            }
        } else {
            return false;
        }
        *npts_out = npts; *ncoord_out = ncoord;
        *coords_out = coords; *resp_out = resp;
        return true;
    }

    if (!fit_is_list(data)) return false;
    npts = fit_len(data);
    if (npts == 0) return false;
    Expr* first = fit_elem(data, 0);

    if (fit_is_list(first)) {
        size_t rowlen = fit_len(first);
        if (rowlen < 2) return false;
        ncoord = rowlen - 1;
        coords = malloc(sizeof(double) * npts * ncoord);
        resp   = malloc(sizeof(double) * npts);
        for (size_t r = 0; r < npts; r++) {
            Expr* row = fit_elem(data, r);
            if (!fit_is_list(row) || fit_len(row) != rowlen) { free(coords); free(resp); return false; }
            for (size_t k = 0; k < ncoord; k++)
                if (!fit_to_double(fit_elem(row, k), &coords[k * npts + r])) {
                    free(coords); free(resp); return false;
                }
            if (!fit_to_double(fit_elem(row, rowlen - 1), &resp[r])) {
                free(coords); free(resp); return false;
            }
        }
    } else {
        ncoord = 1;
        coords = malloc(sizeof(double) * npts);
        resp   = malloc(sizeof(double) * npts);
        for (size_t r = 0; r < npts; r++) {
            Expr* el = fit_elem(data, r);
            if (fit_is_list(el) || !fit_to_double(el, &resp[r])) { free(coords); free(resp); return false; }
            coords[r] = (double)(r + 1);
        }
    }
    *npts_out = npts; *ncoord_out = ncoord;
    *coords_out = coords; *resp_out = resp;
    return true;
}

/* Wrap a contiguous double column as a length-npts float64 NDArray vector
 * (copies the data; the NDArray takes ownership of the copy). */
static Expr* fit_col_vector(const double* col, size_t npts) {
    double* buf = malloc(sizeof(double) * (npts ? npts : 1));
    memcpy(buf, col, sizeof(double) * npts);
    int64_t dims[1] = { (int64_t)npts };
    return expr_new_ndarray_raw(1, dims, buf, NDT_FLOAT64);
}

/* Wrap a row-major double buffer as a packed float64 NDArray (matrix/vector). */
static Expr* fit_pack_matrix(const double* A, size_t rows, size_t cols) {
    size_t nel = rows * cols;
    double* buf = malloc(sizeof(double) * (nel ? nel : 1));
    memcpy(buf, A, sizeof(double) * nel);
    int64_t dims[2] = { (int64_t)rows, (int64_t)cols };
    return expr_new_ndarray_raw(2, dims, buf, NDT_FLOAT64);
}
static Expr* fit_pack_vector(const double* b, size_t n) {
    double* buf = malloc(sizeof(double) * (n ? n : 1));
    memcpy(buf, b, sizeof(double) * n);
    int64_t dims[1] = { (int64_t)n };
    return expr_new_ndarray_raw(1, dims, buf, NDT_FLOAT64);
}

/* Read a coefficient vector (a List of numbers or a real rank-1 NDArray) of
 * length n into out[].  Returns false on a type/shape mismatch. */
static bool fit_coeffs_to_doubles(const Expr* c, size_t n, double* out) {
    if (is_ndarray(c) && c->data.ndarray.rank == 1 &&
        (size_t)c->data.ndarray.dims[0] == n && !ndt_is_complex(c->data.ndarray.dtype)) {
        const void* cb = c->data.ndarray.data; NDType dt = c->data.ndarray.dtype;
        for (size_t j = 0; j < n; j++) { double re, im; ndt_get(cb, j, dt, &re, &im); out[j] = re; }
        return true;
    }
    if (fit_is_list(c) && fit_len(c) == n) {
        for (size_t j = 0; j < n; j++)
            if (!fit_to_double(fit_elem(c, j), &out[j])) return false;
        return true;
    }
    return false;
}

/* Build the npts x nfun design matrix as a row-major double buffer by
 * evaluating each basis function ONCE with the variables bound to the packed
 * coordinate vectors.  Returns a malloc'd buffer, or NULL if a column does not
 * reduce to npts numbers or a broadcastable scalar (caller falls back to the
 * symbolic builder).  varsyms are borrowed. */
static double* fit_build_design_fast(Expr** funs, size_t nfun,
                                     const double* coords, size_t npts, size_t ncoord,
                                     Expr** varsyms) {
    Expr** cvec = calloc(ncoord ? ncoord : 1, sizeof(Expr*));
    MatchEnv* env = env_new();
    for (size_t k = 0; k < ncoord; k++) {
        cvec[k] = fit_col_vector(coords + k * npts, npts);
        env_set(env, varsyms[k]->data.symbol.name, cvec[k]);  /* env copies */
    }

    size_t nel = npts * nfun;
    double* A = malloc(sizeof(double) * (nel ? nel : 1));
    bool ok = true;
    for (size_t j = 0; j < nfun && ok; j++) {
        Expr* col = eval_and_free(replace_bindings(funs[j], env));
        if (is_ndarray(col) && col->data.ndarray.rank == 1 &&
            (size_t)col->data.ndarray.dims[0] == npts &&
            !ndt_is_complex(col->data.ndarray.dtype)) {
            const void* cb = col->data.ndarray.data;
            NDType dt = col->data.ndarray.dtype;
            for (size_t r = 0; r < npts; r++) {
                double re, im; ndt_get(cb, r, dt, &re, &im);
                A[r * nfun + j] = re;
            }
        } else if (fit_is_list(col) && fit_len(col) == npts) {
            for (size_t r = 0; r < npts; r++)
                if (!fit_to_double(fit_elem(col, r), &A[r * nfun + j])) { ok = false; break; }
        } else {
            double s;
            if (fit_to_double(col, &s))                  /* scalar -> broadcast column */
                for (size_t r = 0; r < npts; r++) A[r * nfun + j] = s;
            else
                ok = false;                              /* symbolic / wrong shape */
        }
        expr_free(col);
    }

    for (size_t k = 0; k < ncoord; k++) expr_free(cvec[k]);
    free(cvec);
    env_free(env);
    if (!ok) { free(A); return NULL; }
    return A;
}

/* Column equilibration (numpy.polyfit style): scale each column of the
 * npts x nfun buffer to unit L2 norm in place; return a malloc'd scale[nfun]
 * of the original norms so coefficients can be unscaled.
 *
 * A column that is *numerically zero* -- norm negligible relative to `ref` (the
 * response norm), i.e. a basis function that vanishes over the data, such as
 * Sin[2 x] at multiples of pi -- is zeroed out and flagged with scale 0.0 (a
 * sentinel).  Scaling it to unit norm would amplify rounding noise into an
 * apparently independent direction, and unscaling would then multiply the fitted
 * coefficient by 1/norm ~ 1e16; zeroing it makes the SVD drop it cleanly and its
 * coefficient is forced to exactly 0 on unscale (matching the reference). */
static double* fit_column_scale(double* A, size_t npts, size_t nfun, double ref) {
    double* scale = malloc(sizeof(double) * (nfun ? nfun : 1));
    double zero_tol = 1.0e-13 * (ref > 0.0 ? ref : 1.0);
    for (size_t j = 0; j < nfun; j++) {
        double s = 0.0;
        for (size_t r = 0; r < npts; r++) { double v = A[r * nfun + j]; s += v * v; }
        s = sqrt(s);
        if (!(s > zero_tol)) {                       /* numerically-zero column */
            scale[j] = 0.0;                          /* sentinel -> force coeff 0 */
            for (size_t r = 0; r < npts; r++) A[r * nfun + j] = 0.0;
            continue;
        }
        scale[j] = s;
        double inv = 1.0 / s;
        for (size_t r = 0; r < npts; r++) A[r * nfun + j] *= inv;
    }
    return scale;
}

/* Reciprocal-condition proxy of the column-scaled design matrix, over the
 * GENUINE columns only (scale[j] != 0; the zeroed ones are already known
 * degenerate).  Forms the small Gram G = A^T A of those columns and LU-factors
 * it with partial pivoting, returning min|pivot| / max|pivot|.  ~1 for a
 * well-conditioned basis, -> 0 as columns become linearly dependent.  Since
 * cond(G) ~ cond(A)^2, a tiny ratio flags a fit that even scaling cannot
 * rescue (e.g. a very high-degree monomial basis). */
static double fit_gram_rcond(const double* A, size_t npts, size_t nfun,
                             const double* scale) {
    /* Compact the genuine column indices. */
    size_t* idx = malloc(sizeof(size_t) * (nfun ? nfun : 1));
    size_t m = 0;
    for (size_t j = 0; j < nfun; j++) if (scale[j] != 0.0) idx[m++] = j;
    if (m == 0) { free(idx); return 1.0; }

    double* G = malloc(sizeof(double) * m * m);
    for (size_t p = 0; p < m; p++)
        for (size_t q = p; q < m; q++) {
            double s = 0.0;
            for (size_t r = 0; r < npts; r++) s += A[r * nfun + idx[p]] * A[r * nfun + idx[q]];
            G[p * m + q] = s; G[q * m + p] = s;
        }
    free(idx);
    double pmin = 0.0, pmax = 0.0; bool first = true;
    for (size_t col = 0; col < m; col++) {
        size_t piv = col; double best = fabs(G[col * m + col]);
        for (size_t r = col + 1; r < m; r++) {
            double v = fabs(G[r * m + col]); if (v > best) { best = v; piv = r; }
        }
        if (piv != col)
            for (size_t c = 0; c < m; c++) {
                double t = G[col * m + c]; G[col * m + c] = G[piv * m + c]; G[piv * m + c] = t;
            }
        double diag = G[col * m + col];
        double ad = fabs(diag);
        if (first) { pmin = pmax = ad; first = false; }
        else { if (ad < pmin) pmin = ad; if (ad > pmax) pmax = ad; }
        if (ad < 1e-300) { free(G); return 0.0; }
        for (size_t r = col + 1; r < m; r++) {
            double f = G[r * m + col] / diag;
            if (f == 0.0) continue;
            for (size_t c = col; c < m; c++) G[r * m + c] -= f * G[col * m + c];
        }
    }
    free(G);
    return pmax > 0.0 ? pmin / pmax : 0.0;
}

/* Plain L2 via the packed LAPACK LeastSquares path.
 *
 * Column scaling (numpy.polyfit style) is applied only when the system is
 * OVERdetermined (npts >= nfun): there the least-squares solution is unique, so
 * scaling is numerically transparent -- it merely lets the SVD cut-off see all
 * basis directions on a comparable footing, which is what keeps a high-degree
 * monomial fit accurate.  For an UNDERdetermined system (npts < nfun) the
 * minimiser is not unique and scaling would change which minimum-norm
 * representative LAPACK returns, so it is skipped and the unscaled min-norm
 * solution (matching the reference) is used. */
static bool fit_fast_l2(double* A, const double* b, size_t npts, size_t nfun,
                        double* cf, const char* who) {
    double* scale = NULL;
    if (npts >= nfun) {
        double bnorm = 0.0;
        for (size_t r = 0; r < npts; r++) bnorm += b[r] * b[r];
        bnorm = sqrt(bnorm);
        scale = fit_column_scale(A, npts, nfun, bnorm);
        if (fit_gram_rcond(A, npts, nfun, scale) < 1e-14)
            fit_warn(who, "cond",
                     "the fit may be poorly conditioned (basis functions nearly "
                     "dependent over the data); coefficients may be unreliable");
    }
    Expr* M = fit_pack_matrix(A, npts, nfun);
    Expr* V = fit_pack_vector(b, npts);
    Expr* r = fit_leastsquares(M, V);
    expr_free(M); expr_free(V);
    bool ok = r && fit_coeffs_to_doubles(r, nfun, cf);
    if (r) expr_free(r);
    if (ok && scale)                          /* unscale; scale 0 -> zeroed column */
        for (size_t j = 0; j < nfun; j++) cf[j] = scale[j] != 0.0 ? cf[j] / scale[j] : 0.0;
    free(scale);
    return ok;
}

/* Ridge/Tikhonov: solve [A; sqrt(lam) I] a = [b; 0] via the packed LeastSquares
 * path.  No column scaling (it would distort the ||a||^2 penalty). */
static bool fit_fast_ridge(const double* A, const double* b, size_t npts, size_t nfun,
                           double lambda, double* cf) {
    /* Ridge normal equations: (A^T A + lambda I) c = A^T b -- exactly
     * numpy/sklearn ridge.  lambda I regularises the Gram matrix, so
     * squaring the condition number (the usual objection to normal
     * equations) is not a concern here; and it avoids forming and SVD-ing
     * the (npts+nfun) x nfun augmented system, which for a tall design
     * matrix cost ~8x scipy.  Same construction fit_irls_l1_dense uses. */
    double lam = (lambda < 0.0) ? 0.0 : lambda;
    double* G   = malloc(sizeof(double) * (nfun ? nfun * nfun : 1));
    double* rhs = malloc(sizeof(double) * (nfun ? nfun : 1));
    if (!G || !rhs) { free(G); free(rhs); return false; }
    for (size_t p = 0; p < nfun; p++) {
        for (size_t q = 0; q < nfun; q++) {
            double s = 0.0;
            for (size_t i = 0; i < npts; i++) s += A[i * nfun + p] * A[i * nfun + q];
            G[p * nfun + q] = s + (p == q ? lam : 0.0);
        }
        double s = 0.0;
        for (size_t i = 0; i < npts; i++) s += A[i * nfun + p] * b[i];
        rhs[p] = s;
    }
    bool ok = fit_solve_dense(nfun, G, rhs, cf);
    free(G); free(rhs);
    return ok;
}

/* Arbitrary norm / combined objective: fast-built matrix, FindMinimum solve. */
static bool fit_fast_findmin(const double* A, const double* b, size_t npts, size_t nfun,
                             const FitOpts* opts, double* cf) {
    Expr* M = fit_pack_matrix(A, npts, nfun);
    Expr* V = fit_pack_vector(b, npts);
    Expr* c = fit_solve_findmin(M, V, opts);
    expr_free(M); expr_free(V);
    bool ok = c && fit_coeffs_to_doubles(c, nfun, cf);
    if (c) expr_free(c);
    return ok;
}

/* Clean coefficients that are pure numerical noise relative to the solution
 * scale (|cf| <= 1e-12 * max|cf|) to exactly 0.  A machine least-squares solve
 * leaves ~1e-16 dust in a coefficient that is mathematically zero, and its exact
 * value depends on the BLAS/LAPACK build; snapping keeps the output clean and
 * platform-stable (matching Mathematica's 0.).  The threshold is far below any
 * meaningful coefficient, so genuine terms are untouched. */
static void fit_snap_coeffs(double* cf, size_t nfun) {
    double mx = 0.0;
    for (size_t j = 0; j < nfun; j++) { double a = fabs(cf[j]); if (a > mx) mx = a; }
    double tol = 1.0e-12 * mx;
    for (size_t j = 0; j < nfun; j++) if (fabs(cf[j]) <= tol) cf[j] = 0.0;
}

/* Solve for coefficients given a row-major npts x nfun design buffer A and a
 * length-npts response b, dispatching by norm/regularization exactly as
 * fit_dispatch_solve does but entirely on double buffers.  A may be modified in
 * place (plain-L2 column scaling).  Returns a List of reals, or NULL on solver
 * failure.  Shared by the 3-arg fast path and the Fit[{m,v}] numeric path. */
static Expr* fit_solve_buffers(double* A, const double* b, size_t npts, size_t nfun,
                               const FitOpts* opts, const char* who) {
    double* cf = malloc(sizeof(double) * (nfun ? nfun : 1));
    bool ok = false;
    FitNormKind nk = opts->norm_fun ? fit_classify_norm(opts->norm_fun) : FIT_NORM_L2;

    if (opts->reg_kind == FIT_REG_LASSO) {
        if (opts->norm_fun && nk != FIT_NORM_L2) {
            ok = fit_fast_findmin(A, b, npts, nfun, opts, cf);
        } else {
            double lambda;
            if (fit_to_double(opts->reg_lambda, &lambda) && lambda >= 0.0) {
                double* a = fit_lasso_dense(A, b, npts, nfun, lambda);
                memcpy(cf, a, sizeof(double) * nfun); free(a); ok = true;
            } else {
                fit_warn(who, "reglam", "regularization parameter must be non-negative");
            }
        }
    } else if (nk == FIT_NORM_L1 && opts->reg_kind == FIT_REG_NONE) {
        bool sol_ok = true;
        double* a = fit_irls_l1_dense(A, b, npts, nfun, &sol_ok);
        if (sol_ok) { memcpy(cf, a, sizeof(double) * nfun); free(a); ok = true; }
        else fit_warn(who, "l1sing", "1-norm fit failed (singular system)");
    } else if (nk == FIT_NORM_OTHER || (nk == FIT_NORM_L1 && opts->reg_kind != FIT_REG_NONE)) {
        ok = fit_fast_findmin(A, b, npts, nfun, opts, cf);
    } else if (opts->reg_kind == FIT_REG_TIKHONOV) {
        double lambda;
        if (fit_to_double(opts->reg_lambda, &lambda) && lambda >= 0.0)
            ok = fit_fast_ridge(A, b, npts, nfun, lambda, cf);
        else
            fit_warn(who, "reglam", "regularization parameter must be non-negative");
    } else {
        ok = fit_fast_l2(A, b, npts, nfun, cf, who);
    }

    if (!ok) { free(cf); return NULL; }
    fit_snap_coeffs(cf, nfun);
    Expr* coeffs = fit_doubles_to_vector(cf, nfun);
    free(cf);
    return coeffs;
}

/* Read a design matrix (List-of-Lists or real rank-2 NDArray) into a fresh
 * row-major double buffer; NULL if not numeric.  Caller frees. */
static double* fit_read_matrix_buffer(const Expr* m, size_t* rows, size_t* cols) {
    if (is_ndarray(m) && m->data.ndarray.rank == 2 && !ndt_is_complex(m->data.ndarray.dtype)) {
        size_t R = (size_t)m->data.ndarray.dims[0], C = (size_t)m->data.ndarray.dims[1];
        if (R == 0 || C == 0) return NULL;
        double* A = malloc(sizeof(double) * R * C);
        for (size_t i = 0; i < R * C; i++) {
            double re, im; ndt_get(m->data.ndarray.data, i, m->data.ndarray.dtype, &re, &im);
            A[i] = re;
        }
        *rows = R; *cols = C;
        return A;
    }
    return fit_matrix_doubles(m, rows, cols);
}

/* Read a response vector (List or real rank-1 NDArray) into a fresh double
 * buffer; NULL if not numeric.  Caller frees. */
static double* fit_read_vector_buffer(const Expr* v, size_t* n) {
    if (is_ndarray(v) && v->data.ndarray.rank == 1 && !ndt_is_complex(v->data.ndarray.dtype)) {
        size_t N = (size_t)v->data.ndarray.dims[0];
        double* b = malloc(sizeof(double) * (N ? N : 1));
        for (size_t i = 0; i < N; i++) {
            double re, im; ndt_get(v->data.ndarray.data, i, v->data.ndarray.dtype, &re, &im);
            b[i] = re;
        }
        *n = N;
        return b;
    }
    return fit_vector_doubles(v, n);
}

/* Attempt the machine-precision packed fast path for the 3-arg Fit form.  On
 * success sets *handled=true and returns the coefficient vector (a List of
 * reals).  If the fast path does not apply — exact/MPFR precision, non-numeric
 * data, a variable/coordinate mismatch, or a basis column that stays symbolic —
 * sets *handled=false and returns NULL so the caller runs the symbolic path
 * (which emits any diagnostics).  On a genuine numeric solver failure sets
 * *handled=true and returns NULL. */
static Expr* fit_fast_3arg(const char* who, const Expr* data, Expr* vars,
                           Expr** funs, size_t nfun, const FitOpts* opts,
                           bool* handled) {
    *handled = false;
    if (opts->wp_mode == FIT_WP_INFINITY || opts->wp_mode == FIT_WP_MPFR)
        return NULL;
#ifdef USE_MPFR
    /* WorkingPrecision -> Automatic preserves the precision of MPFR input, so
     * MPFR data must take the symbolic/MPFR path, not the machine fast path.
     * (An explicit WorkingPrecision -> MachinePrecision still fast-paths.) */
    if (opts->wp_mode == FIT_WP_AUTO && fit_max_mpfr_bits(data) > 0)
        return NULL;
#endif

    size_t npts = 0, ncoord = 0;
    double* coords = NULL; double* resp = NULL;
    if (!fit_read_numeric(data, &npts, &ncoord, &coords, &resp))
        return NULL;

    /* Extract variable symbols quietly; the symbolic path emits diagnostics
     * on a mismatch.  varsyms are borrowed pointers into `vars`. */
    size_t nvar = fit_is_list(vars) ? fit_len(vars) : 1;
    if (nvar != ncoord) { free(coords); free(resp); return NULL; }
    Expr** varsyms = calloc(ncoord ? ncoord : 1, sizeof(Expr*));
    bool vok = true;
    if (fit_is_list(vars)) {
        for (size_t k = 0; k < ncoord; k++) {
            Expr* s = fit_elem(vars, k);
            if (s->type != EXPR_SYMBOL) { vok = false; break; }
            varsyms[k] = s;
        }
    } else if (vars->type == EXPR_SYMBOL) {
        varsyms[0] = vars;
    } else {
        vok = false;
    }
    if (!vok) { free(varsyms); free(coords); free(resp); return NULL; }

    double* A = fit_build_design_fast(funs, nfun, coords, npts, ncoord, varsyms);
    free(varsyms);
    free(coords);
    if (!A) { free(resp); return NULL; }

    /* From here the fast path owns the outcome. */
    *handled = true;
    Expr* coeffs = fit_solve_buffers(A, resp, npts, nfun, opts, who);
    free(A);
    free(resp);
    return coeffs;
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

    /* Numeric machine fast path: read m/v (List-of-Lists or NDArray) into double
     * buffers and solve on the buffer, tripping the LAPACK/dense solvers. */
    if (opts->wp_mode != FIT_WP_INFINITY && opts->wp_mode != FIT_WP_MPFR) {
        size_t rows = 0, cols = 0, vn = 0;
        double* A = fit_read_matrix_buffer(m, &rows, &cols);
        double* b = A ? fit_read_vector_buffer(v, &vn) : NULL;
        if (A && b && rows > 0 && cols > 0 && vn == rows) {
            Expr* coeffs = fit_solve_buffers(A, b, rows, cols, opts, "Fit");
            free(A); free(b);
            return coeffs;   /* NULL only on a genuine solver failure */
        }
        free(A); free(b);
    }

    /* Fallback: exact/MPFR, or a symbolic design matrix. */
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

    size_t nfun = fit_len(funs_list);

    /* Machine-precision packed fast path (numeric data, machine WorkingPrecision).
     * Builds a packed float64 design matrix column-wise and routes the solve
     * through LAPACK/dense buffers; falls through to the symbolic path below on
     * exact/MPFR/symbolic input. */
    {
        bool handled = false;
        Expr* fastc = fit_fast_3arg("Fit", data, vars,
                                    funs_list->data.function.args, nfun, &opts, &handled);
        if (handled) {
            Expr* out = fastc
                ? fit_assemble_output(fastc, funs_list->data.function.args, nfun)
                : NULL;
            if (fastc) expr_free(fastc);
            return out;
        }
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
