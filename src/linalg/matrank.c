/* matrank.c
 *
 * MatrixRank[m]                       -- rank of m (exact path through
 *                                       RowReduce when m has no Real /
 *                                       MPFR / Complex-of-Real entries).
 * MatrixRank[m, Method  -> "<name>"]  -- explicit RREF method dispatch
 *                                       (exact path only).
 * MatrixRank[m, Tolerance -> t]       -- treat |entry| <= t as zero
 *                                       during pivot selection.
 * MatrixRank[m, Method->..., Tolerance->...]
 *                                     -- both options simultaneously.
 *
 * Algorithm dispatch:
 *
 *   - Numerical path (matrix has any Real / MPFR leaf, or every leaf
 *     converts to a complex double, OR the user supplies Tolerance):
 *     run a partial-pivot Gaussian forward-elimination over
 *     `double _Complex` (here modelled as a {re, im} struct so we stay
 *     portable to compilers without _Complex).  A column is skipped
 *     when its largest sub-pivot |entry| is <= tolerance; the rank is
 *     the number of accepted pivots.
 *
 *   - Exact path (every leaf is exact AND no Tolerance supplied):
 *     route through `RowReduce[m, Method->...]` and count the number
 *     of RREF rows whose entries are not all structurally zero.  This
 *     gives the rank with no precision concerns.
 *
 * Default tolerance (when Tolerance -> Automatic, the default):
 *   - If the matrix contains any Real / MPFR leaf:
 *       tol = max(rows, cols) * DBL_EPSILON * max(|entries|)
 *     (the standard "rank-by-SVD" surrogate; we substitute max(|entries|)
 *     for the largest singular value, which is exact in the row /
 *     column-scaled cases and at most a small constant factor off in
 *     general).
 *   - Otherwise tol = 0, so the numerical path agrees with the exact
 *     path on integer / rational matrices.
 *
 * Memory ownership: standard builtin contract.  This file does NOT
 * free `res` on success or failure -- the evaluator owns it (see
 * MEMORY.md / SPEC.md §4.1).
 */

#include "matrank.h"
#include "linalg.h"
#include "linsolve.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "sym_names.h"
#include "flint_mat_bridge.h"

#include <float.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ *
 *  Lightweight complex-double type and arithmetic.                    *
 * ------------------------------------------------------------------ */
typedef struct { double re, im; } cplx_t;

static cplx_t c_make(double re, double im) { cplx_t z = { re, im }; return z; }
static cplx_t c_sub(cplx_t a, cplx_t b)    { return c_make(a.re - b.re, a.im - b.im); }
static cplx_t c_mul(cplx_t a, cplx_t b)    {
    return c_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static cplx_t c_div(cplx_t a, cplx_t b) {
    double d = b.re * b.re + b.im * b.im;
    return c_make((a.re * b.re + a.im * b.im) / d,
                  (a.im * b.re - a.re * b.im) / d);
}
static double c_abs(cplx_t a) { return hypot(a.re, a.im); }

/* ------------------------------------------------------------------ *
 *  Option parsing.  We accept any mixture of                          *
 *      Method    -> "<name>" | Automatic                              *
 *      Tolerance -> Automatic | <non-negative number>                 *
 *  in any order, but at most one of each.  Unrecognised options leave *
 *  the call unevaluated (return false here).                          *
 * ------------------------------------------------------------------ */
typedef struct {
    MatsolMethod method;
    bool         tol_is_automatic;
    double       tolerance;
} MatrankOpts;

static bool parse_tolerance_value(Expr* rhs, bool* is_auto_out, double* tol_out) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
        *is_auto_out = true;
        *tol_out = 0.0;
        return true;
    }
    if (rhs->type == EXPR_INTEGER) {
        if (rhs->data.integer < 0) return false;
        *is_auto_out = false;
        *tol_out = (double)rhs->data.integer;
        return true;
    }
    if (rhs->type == EXPR_REAL) {
        if (rhs->data.real < 0.0) return false;
        *is_auto_out = false;
        *tol_out = rhs->data.real;
        return true;
    }
    if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol == SYM_Rational
        && rhs->data.function.arg_count == 2) {
        Expr* n = rhs->data.function.args[0];
        Expr* d = rhs->data.function.args[1];
        if (n->type != EXPR_INTEGER || d->type != EXPR_INTEGER) return false;
        if (n->data.integer < 0 || d->data.integer <= 0) return false;
        *is_auto_out = false;
        *tol_out = (double)n->data.integer / (double)d->data.integer;
        return true;
    }
    /* Power[10, k] / Power[2, k] for small k -- common in the spec
     * examples (Tolerance -> 10^-8).  Fall back to evaluating N[rhs]
     * to a real. */
    Expr* n = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy(rhs) }, 1));
    bool ok = false;
    double v = 0.0;
    if (n->type == EXPR_REAL)            { v = n->data.real;                ok = (v >= 0.0); }
    else if (n->type == EXPR_INTEGER)    { v = (double)n->data.integer;     ok = (v >= 0.0); }
    else if (n->type == EXPR_BIGINT)     { v = mpz_get_d(n->data.bigint);   ok = (v >= 0.0); }
#ifdef USE_MPFR
    else if (n->type == EXPR_MPFR)       { v = mpfr_get_d(n->data.mpfr, MPFR_RNDN); ok = (v >= 0.0); }
#endif
    expr_free(n);
    if (!ok) return false;
    *is_auto_out = false;
    *tol_out = v;
    return true;
}

static bool parse_options(Expr* res, MatrankOpts* opts) {
    opts->method            = MATSOL_AUTOMATIC;
    opts->tol_is_automatic  = true;
    opts->tolerance         = 0.0;

    size_t argc = res->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type != EXPR_FUNCTION
            || opt->data.function.head->type != EXPR_SYMBOL
            || opt->data.function.arg_count != 2) return false;
        const char* hd = opt->data.function.head->data.symbol;
        if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL) return false;

        if (lhs->data.symbol == SYM_Method) {
            MatsolMethod m = matsol_parse_method_option(opt);
            if (m == MATSOL_INVALID) return false;
            opts->method = m;
        } else if (lhs->data.symbol == SYM_Tolerance) {
            bool is_auto;
            double tol;
            if (!parse_tolerance_value(rhs, &is_auto, &tol)) return false;
            opts->tol_is_automatic = is_auto;
            opts->tolerance        = tol;
        } else {
            return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Inexact-leaf detection and entry-to-cplx conversion.               *
 * ------------------------------------------------------------------ */
static bool has_inexact_leaf(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    if (has_inexact_leaf(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_inexact_leaf(e->data.function.args[i])) return true;
    }
    return false;
}

/* Convert an Expr leaf to a cplx_t.  Returns false on any symbolic
 * (non-numeric) input.  Accepts:
 *   - EXPR_INTEGER, EXPR_BIGINT, EXPR_REAL, EXPR_MPFR
 *   - Rational[p, q] with integer or bigint p, q
 *   - Complex[re, im] where re, im recursively convert
 *   - the symbol I -> (0, 1)
 *   - Times[a, b, ...] / Plus[a, b, ...] of the above (evaluated to a
 *     concrete number by Mathilda before MatrixRank is called)
 */
static bool entry_to_cplx(Expr* e, cplx_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = c_make((double)e->data.integer, 0); return true; }
    if (e->type == EXPR_REAL)    { *out = c_make(e->data.real, 0);            return true; }
    if (e->type == EXPR_BIGINT)  { *out = c_make(mpz_get_d(e->data.bigint), 0); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = c_make(mpfr_get_d(e->data.mpfr, MPFR_RNDN), 0); return true; }
#endif
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol, "I") == 0) {
        *out = c_make(0.0, 1.0);
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* hd = e->data.function.head->data.symbol;
        if (hd == SYM_Rational && e->data.function.arg_count == 2) {
            cplx_t p, q;
            if (!entry_to_cplx(e->data.function.args[0], &p)) return false;
            if (!entry_to_cplx(e->data.function.args[1], &q)) return false;
            if (q.im != 0.0 || q.re == 0.0) return false;
            *out = c_make(p.re / q.re, p.im / q.re);
            return true;
        }
        if (hd == SYM_Complex && e->data.function.arg_count == 2) {
            cplx_t a, b;
            if (!entry_to_cplx(e->data.function.args[0], &a)) return false;
            if (!entry_to_cplx(e->data.function.args[1], &b)) return false;
            /* Complex[a, b] = a + b * I.  In a well-formed Complex
             * literal a and b are both real, but be defensive. */
            *out = c_make(a.re - b.im, a.im + b.re);
            return true;
        }
    }
    /* Fall back to N[expr] for symbolic-but-evaluates-to-number cases
     * such as Pi, E, Sqrt[2], etc.  N returns a Real or Complex literal
     * which we then re-feed through this function. */
    Expr* n = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy(e) }, 1));
    bool ok = false;
    /* Prevent infinite recursion: only recurse when N produced
     * something simpler than `e`. */
    if (n->type == EXPR_REAL
        || n->type == EXPR_INTEGER
        || n->type == EXPR_BIGINT
#ifdef USE_MPFR
        || n->type == EXPR_MPFR
#endif
        || (n->type == EXPR_FUNCTION
            && n->data.function.head->type == EXPR_SYMBOL
            && n->data.function.head->data.symbol == SYM_Complex)) {
        ok = entry_to_cplx(n, out);
    }
    expr_free(n);
    return ok;
}

/* Try to convert every entry of a flattened matrix to cplx_t.  On
 * success returns a freshly-malloced array of length rows*cols (caller
 * frees).  On any non-numeric entry returns NULL. */
static cplx_t* matrix_to_cplx(Expr** flat, int rows, int cols) {
    cplx_t* M = malloc(sizeof(cplx_t) * (size_t)rows * (size_t)cols);
    if (!M) return NULL;
    for (int i = 0; i < rows * cols; i++) {
        if (!entry_to_cplx(flat[i], &M[i])) {
            free(M);
            return NULL;
        }
    }
    return M;
}

/* ------------------------------------------------------------------ *
 *  Numerical rank via partial-pivot Gaussian forward elimination.     *
 * ------------------------------------------------------------------ */
static int gauss_rank_cplx(cplx_t* M, int rows, int cols, double tol) {
    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int pivot_row = -1;
        double max_abs = tol;
        for (int i = r; i < rows; i++) {
            double v = c_abs(M[i * cols + c]);
            if (v > max_abs) {
                max_abs = v;
                pivot_row = i;
            }
        }
        if (pivot_row < 0) continue;

        if (pivot_row != r) {
            for (int j = c; j < cols; j++) {
                cplx_t tmp = M[r * cols + j];
                M[r * cols + j]         = M[pivot_row * cols + j];
                M[pivot_row * cols + j] = tmp;
            }
        }
        cplx_t pivot = M[r * cols + c];
        for (int i = r + 1; i < rows; i++) {
            cplx_t factor = c_div(M[i * cols + c], pivot);
            M[i * cols + c] = c_make(0.0, 0.0);
            for (int j = c + 1; j < cols; j++) {
                M[i * cols + j] = c_sub(M[i * cols + j],
                                        c_mul(factor, M[r * cols + j]));
            }
        }
        r++;
    }
    return r;
}

/* ------------------------------------------------------------------ *
 *  Exact path: RowReduce -> count non-zero rows.                      *
 * ------------------------------------------------------------------ */
static Expr* call_rowreduce(Expr* m, MatsolMethod method) {
    if (method == MATSOL_AUTOMATIC) {
        Expr** args = malloc(sizeof(Expr*) * 1);
        args[0] = expr_copy(m);
        Expr* call = expr_new_function(expr_new_symbol(SYM_RowReduce), args, 1);
        free(args);
        return eval_and_free(call);
    }
    const char* name = NULL;
    switch (method) {
        case MATSOL_DIVFREE:  name = "DivisionFreeRowReduction"; break;
        case MATSOL_ONESTEP:  name = "OneStepRowReduction";      break;
        case MATSOL_COFACTOR: name = "CofactorExpansion";        break;
        default:              name = "Automatic";                break;
    }
    Expr** rule_args = malloc(sizeof(Expr*) * 2);
    rule_args[0] = expr_new_symbol(SYM_Method);
    rule_args[1] = expr_new_string(name);
    Expr* opt = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
    free(rule_args);

    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(m);
    args[1] = opt;
    Expr* call = expr_new_function(expr_new_symbol(SYM_RowReduce), args, 2);
    free(args);
    return eval_and_free(call);
}

static int count_nonzero_rows(Expr* rref, int rows, int cols) {
    int r = 0;
    for (int i = 0; i < rows; i++) {
        Expr* row = rref->data.function.args[i];
        bool nonzero = false;
        for (int j = 0; j < cols; j++) {
            if (!is_zero_poly(row->data.function.args[j])) {
                nonzero = true;
                break;
            }
        }
        if (nonzero) r++;
    }
    return r;
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */
Expr* builtin_matrixrank(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* m = res->data.function.args[0];

    MatrankOpts opts;
    if (!parse_options(res, &opts)) {
        /* Either an unknown Method value or an unknown option key.
         * Mirror NullSpace's behaviour: emit a one-shot diagnostic for
         * Method, and leave the call unevaluated either way. */
        static uint64_t last_warned = 0;
        matsol_warn_once(&last_warned, res,
            "MatrixRank::opt: Option is not one of Method -> "
            "{\"Automatic\", \"DivisionFreeRowReduction\", "
            "\"OneStepRowReduction\", \"CofactorExpansion\"} or "
            "Tolerance -> Automatic | non-negative number.\n");
        return NULL;
    }

    int64_t dims[64];
    int rank2 = get_tensor_dims(m, dims);
    if (rank2 != 2 || dims[0] == 0 || dims[1] == 0) {
        char* s = expr_to_string(m);
        fprintf(stderr,
                "MatrixRank::matrix: Argument %s at position 1 is not a "
                "non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }

    int rows = (int)dims[0];
    int cols = (int)dims[1];

    /* Decide which path to take.  Numerical when:
     *   - User explicitly supplied a finite Tolerance (>= 0), OR
     *   - The matrix has any inexact leaf and the user did not
     *     opt out via Tolerance -> 0.                       */
    bool inexact = has_inexact_leaf(m);
    bool numerical_path = !opts.tol_is_automatic || inexact;

    Expr** flat = malloc(sizeof(Expr*) * (size_t)rows * (size_t)cols);
    {
        size_t idx = 0;
        flatten_tensor(m, flat, &idx);
    }

    int rank_value = -1;

    if (numerical_path) {
        cplx_t* M = matrix_to_cplx(flat, rows, cols);
        if (M) {
            double tol = opts.tolerance;
            if (opts.tol_is_automatic) {
                /* Default for inexact matrices: max(m, n) * eps *
                 * max(|entries|).  For an all-zero matrix this is 0,
                 * matching the exact path. */
                double max_abs = 0.0;
                for (int i = 0; i < rows * cols; i++) {
                    double v = c_abs(M[i]);
                    if (v > max_abs) max_abs = v;
                }
                int dim = rows > cols ? rows : cols;
                tol = (double)dim * DBL_EPSILON * max_abs;
            }
            rank_value = gauss_rank_cplx(M, rows, cols, tol);
            free(M);
        }
        /* If conversion failed, the matrix has symbolic entries.
         * Fall through to the exact path; the user-supplied tolerance
         * (if any) is irrelevant for symbolic entries. */
    }

    if (rank_value < 0) {
        /* Exact path.  An all-rational matrix gets its rank directly from
         * FLINT (fmpq_mat_rref) without materialising the full RREF Expr;
         * rank is basis-independent so the value matches the classical count.
         * Non-rational (symbolic) matrices return -1 and use RowReduce. */
        int frank = flint_mat_rank(flat, rows, cols);
        if (frank >= 0) rank_value = frank;
    }

    if (rank_value < 0) {
        Expr* rref = call_rowreduce(m, opts.method);
        if (rref) {
            int64_t rref_dims[64];
            int rref_rank = get_tensor_dims(rref, rref_dims);
            if (rref_rank == 2 && rref_dims[0] == rows && rref_dims[1] == cols) {
                rank_value = count_nonzero_rows(rref, rows, cols);
            }
            expr_free(rref);
        }
    }

    for (int i = 0; i < rows * cols; i++) expr_free(flat[i]);
    free(flat);

    if (rank_value < 0) return NULL;
    return expr_new_integer((int64_t)rank_value);
}

void matrank_init(void) {
    symtab_add_builtin("MatrixRank", builtin_matrixrank);
    symtab_get_def("MatrixRank")->attributes |= ATTR_PROTECTED;
}
