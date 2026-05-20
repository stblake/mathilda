#include "eigen.h"
#include "eigen_internal.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "poly.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "common.h"
#include "numeric.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ============================================================ *
 *  Eigenvalues and Eigenvectors                                  *
 *                                                                 *
 *  Strategy: compute the characteristic polynomial p(λ) =         *
 *  Det[m - λ*I] (or Det[m - λ*a] for generalised eigenvalues),    *
 *  then route to the public `Solve` builtin so its rationalise -> *
 *  solve -> numericalize pipeline handles approximate matrices    *
 *  automatically.                                                 *
 *                                                                 *
 *  The internal-only lambda symbol is context-qualified           *
 *  ("Eigenvalues`Lambda") so it does not collide with user names. *
 * ============================================================ */

/* Canonical interned name of the lambda symbol used as the char-poly
 * variable.  Pointer is stable across all calls (interner is idempotent). */
const char* eigen_lambda_name(void) {
    static const char* s = NULL;
    if (!s) s = intern_symbol("Eigenvalues`Lambda");
    return s;
}

/* True iff `m` is a non-empty n×n list-of-lists matrix.  Writes n to
 * *n_out on success. */
static bool eigen_is_square_matrix(Expr* m, int64_t* n_out) {
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) return false;
    *n_out = dims[0];
    return true;
}

/* Classify `arg`:
 *   {m, a} with both m and a square n×n  -> generalised eigenvalue case
 *   single square n×n                    -> ordinary eigenvalue case
 * On success, *m_out, *a_out (NULL for ordinary) and *n_out are set.
 * Returns false when `arg` is not a recognisable matrix shape. */
bool eigen_extract_matrix_pair(Expr* arg, Expr** m_out, Expr** a_out,
                                       int64_t* n_out) {
    if (arg->type == EXPR_FUNCTION
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol == SYM_List
        && arg->data.function.arg_count == 2) {
        Expr* m = arg->data.function.args[0];
        Expr* a = arg->data.function.args[1];
        int64_t nm, na;
        if (eigen_is_square_matrix(m, &nm)
            && eigen_is_square_matrix(a, &na) && nm == na) {
            *m_out = m; *a_out = a; *n_out = nm;
            return true;
        }
    }
    int64_t n;
    if (eigen_is_square_matrix(arg, &n)) {
        *m_out = arg; *a_out = NULL; *n_out = n;
        return true;
    }
    return false;
}

/* Build the n×n matrix m - λ*a (or m - λ*I when a_or_null == NULL) where
 * each entry is an evaluated polynomial-in-lambda expression.  Caller
 * owns the returned List-of-Lists. */
Expr* eigen_build_lambda_matrix(Expr* m, Expr* a_or_null,
                                        const char* lambda_name, int64_t n) {
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int64_t i = 0; i < n; i++) {
        Expr* m_row = m->data.function.args[i];
        Expr* a_row = a_or_null ? a_or_null->data.function.args[i] : NULL;
        Expr** elems = malloc(sizeof(Expr*) * n);
        for (int64_t j = 0; j < n; j++) {
            Expr* mij = m_row->data.function.args[j];
            Expr* sub;
            if (a_or_null) {
                Expr* aij = a_row->data.function.args[j];
                Expr* neg_lam_a = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1),
                               expr_new_symbol(lambda_name),
                               expr_copy(aij) }, 3));
                sub = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij), neg_lam_a }, 2));
            } else if (i == j) {
                Expr* neg_lam = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1),
                               expr_new_symbol(lambda_name) }, 2));
                sub = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij), neg_lam }, 2));
            } else {
                sub = expr_copy(mij);
            }
            elems[j] = sub;
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), elems, n);
        free(elems);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return result;
}

/* Compute Det of the matrix `matrix` (n×n) using existing Laplace expansion.
 * Returns a freshly allocated Expr*. */
Expr* eigen_compute_det(Expr* matrix, int n) {
    Expr** flat = malloc(sizeof(Expr*) * n * n);
    size_t idx = 0;
    flatten_tensor(matrix, flat, &idx);
    int* cols = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) cols[i] = i;
    Expr* det_val = laplace_det(flat, n, n, 0, cols);
    free(cols);
    for (size_t i = 0; i < idx; i++) expr_free(flat[i]);
    free(flat);
    return det_val;
}

/* ---- Faddeev-Leverrier char-poly fast path for the ordinary case ----
 *
 * The Laplace-expansion det used by eigen_compute_det is O(n!) -- usable
 * up to about n = 8 but quickly intolerable beyond.  When the user just
 * wants the eigenvalues of an n×n matrix (no generalised `a`), we can run
 * Faddeev-Leverrier instead: it computes the coefficients of
 * det(λI - A) directly in O(n^4) matrix operations on the constant
 * matrix A, never touching a polynomial-in-λ entry.  The resulting
 * polynomial has the same roots as det(A - λI) (they differ only by the
 * (-1)^n sign, which Solve ignores).
 */

/* Trace of an n×n matrix expressed as List of Lists. */
static Expr* eigen_mat_trace(Expr* M, int n) {
    Expr** terms = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        terms[i] = expr_copy(
            M->data.function.args[i]->data.function.args[i]);
    }
    Expr* res = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), terms, n));
    free(terms);
    return res;
}

/* Return M - s*I (entrywise) as a freshly allocated matrix. */
static Expr* eigen_mat_minus_scalar_id(Expr* M, Expr* s, int n) {
    Expr** rows = malloc(sizeof(Expr*) * n);
    for (int i = 0; i < n; i++) {
        Expr** row = malloc(sizeof(Expr*) * n);
        for (int j = 0; j < n; j++) {
            Expr* mij = M->data.function.args[i]->data.function.args[j];
            if (i == j) {
                row[j] = eval_and_free(expr_new_function(
                    expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(mij),
                               eval_and_free(expr_new_function(
                                   expr_new_symbol("Times"),
                                   (Expr*[]){ expr_new_integer(-1),
                                              expr_copy(s) }, 2)) }, 2));
            } else {
                row[j] = expr_copy(mij);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row, n);
        free(row);
    }
    Expr* result = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return result;
}

/* Matrix multiply via dot2 + evaluate (so the entries are normalised). */
static Expr* eigen_mat_mul(Expr* A, Expr* B) {
    bool err = false;
    Expr* prod = dot2(A, B, &err);
    if (!prod) return NULL;
    return eval_and_free(prod);
}

/* Faddeev-Leverrier-Souriau characteristic polynomial.
 *
 *     M_1 = A,    p_1 = Tr(M_1)
 *     for k = 2, ..., n:
 *         M_k = A . (M_{k-1} − p_{k-1} I)
 *         p_k = Tr(M_k) / k
 *
 * Yields p(λ) = det(λI − A) = λ^n − p_1 λ^{n-1} − p_2 λ^{n-2} − … − p_n.
 * Every coefficient other than the leading one is −p_k (no alternating sign).
 *
 * Returns the polynomial in the lambda variable.  Caller owns the result. */
Expr* eigen_char_poly_faddeev(Expr* A, const char* lambda_name, int n) {
    Expr* M = expr_copy(A);                              /* M_1 = A */
    Expr* p_prev = eigen_mat_trace(M, n);                /* p_1 */

    Expr** coeffs = malloc(sizeof(Expr*) * (n + 1));
    coeffs[n] = expr_new_integer(1);
    coeffs[n - 1] = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), expr_copy(p_prev) }, 2));

    for (int k = 2; k <= n; k++) {
        Expr* shifted = eigen_mat_minus_scalar_id(M, p_prev, n);
        expr_free(M);
        M = eigen_mat_mul(A, shifted);
        expr_free(shifted);
        if (!M) {
            expr_free(p_prev);
            expr_free(coeffs[n]);
            expr_free(coeffs[n - 1]);
            free(coeffs);
            return NULL;
        }
        Expr* tr = eigen_mat_trace(M, n);
        Expr* p_k = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ tr,
                       eval_and_free(expr_new_function(
                           expr_new_symbol("Power"),
                           (Expr*[]){ expr_new_integer(k),
                                      expr_new_integer(-1) }, 2)) }, 2));
        expr_free(p_prev);
        p_prev = p_k;

        /* coefficient of λ^{n-k} in det(λI − A) is −p_k. */
        coeffs[n - k] = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(p_k) }, 2));
    }

    expr_free(p_prev);
    expr_free(M);

    /* Build polynomial in lambda_name. */
    Expr** terms = malloc(sizeof(Expr*) * (n + 1));
    size_t tcount = 0;
    for (int k = 0; k <= n; k++) {
        if (k == 0) {
            terms[tcount++] = coeffs[0];
            continue;
        }
        Expr* lam_pow;
        if (k == 1) {
            lam_pow = expr_new_symbol(lambda_name);
        } else {
            lam_pow = eval_and_free(expr_new_function(
                expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol(lambda_name),
                           expr_new_integer(k) }, 2));
        }
        terms[tcount++] = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ coeffs[k], lam_pow }, 2));
    }
    free(coeffs);
    Expr* poly = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), terms, tcount));
    free(terms);
    return poly;
}

/* Solve poly == 0 for the lambda variable via the public Solve builtin so
 * the inexact-preprocessing pipeline (rationalise -> solve -> numericalize)
 * runs automatically.  Cubics/Quartics options are forwarded so the user
 * can request held Root[] objects via Cubics -> False / Quartics -> False.
 * Returns the solution List, or NULL on failure. */
Expr* eigen_solve_poly(Expr* poly, const char* lambda_name,
                              bool cubics_radical, bool quartics_radical) {
    Expr* eq = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ expr_copy(poly), expr_new_integer(0) }, 2);
    Expr* lam = expr_new_symbol(lambda_name);
    Expr* opt_cubics = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol("Cubics"),
                   expr_new_symbol(cubics_radical ? "True" : "False") }, 2);
    Expr* opt_quartics = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol("Quartics"),
                   expr_new_symbol(quartics_radical ? "True" : "False") }, 2);
    Expr* solve_call = expr_new_function(expr_new_symbol("Solve"),
        (Expr*[]){ eq, lam, opt_cubics, opt_quartics }, 4);
    return eval_and_free(solve_call);
}

/* Extract eigenvalues from Solve's output:
 *   {{λ -> v1}, {λ -> v2}, ...}  ->  freshly-owned array [v1, v2, ...]
 * Empty solution `{{}}` (tautology -- e.g. the input was the zero matrix
 * and m == a generalised case) yields a NULL value placeholder so callers
 * can pad with Indeterminate. */
Expr** eigen_extract_values(Expr* solutions, size_t* count_out) {
    *count_out = 0;
    if (!solutions || solutions->type != EXPR_FUNCTION
        || solutions->data.function.head->type != EXPR_SYMBOL
        || solutions->data.function.head->data.symbol != SYM_List) {
        return NULL;
    }
    size_t n = solutions->data.function.arg_count;
    if (n == 0) return NULL;
    Expr** out = malloc(sizeof(Expr*) * n);
    size_t out_count = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* sol = solutions->data.function.args[i];
        if (sol->type != EXPR_FUNCTION
            || sol->data.function.head->type != EXPR_SYMBOL
            || sol->data.function.head->data.symbol != SYM_List
            || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (rule->type != EXPR_FUNCTION
            || rule->data.function.head->type != EXPR_SYMBOL
            || rule->data.function.head->data.symbol != SYM_Rule
            || rule->data.function.arg_count != 2) continue;
        out[out_count++] = expr_copy(rule->data.function.args[1]);
    }
    *count_out = out_count;
    return out;
}

/* Chop very small inexact imaginary / real parts.  Threshold is relative
 * to the magnitude of the value: anything below 1e-10 * |val| (or below
 * 1e-12 absolute when |val| is near zero) is dropped.  Used to clean up
 * numerical noise introduced by Cardano-style closed-form root formulas
 * applied to real polynomials. */
static double eigen_part_to_double(Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_BIGINT) return mpz_get_d(e->data.bigint);
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        double p = eigen_part_to_double(e->data.function.args[0]);
        double q = eigen_part_to_double(e->data.function.args[1]);
        return q == 0 ? 0 : p / q;
    }
    return NAN;
}

Expr* eigen_chop(Expr* val) {
    if (!val) return NULL;
    /* Complex[re, im] - drop im (or re) if small. */
    if (val->type == EXPR_FUNCTION
        && val->data.function.head->type == EXPR_SYMBOL
        && val->data.function.head->data.symbol == SYM_Complex
        && val->data.function.arg_count == 2) {
        Expr* re = val->data.function.args[0];
        Expr* im = val->data.function.args[1];
        double rd = eigen_part_to_double(re);
        double id = eigen_part_to_double(im);
        if (!isnan(rd) && !isnan(id)) {
            double mag = fabs(rd) + fabs(id);
            double thresh = 1e-10 * mag + 1e-12;
            bool drop_im = fabs(id) < thresh;
            bool drop_re = fabs(rd) < thresh;
            if (drop_im && drop_re) return expr_new_real(0.0);
            if (drop_im) return expr_copy(re);
            if (drop_re) {
                /* Pure imaginary: i * im */
                return eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(im), expr_new_symbol("I") }, 2));
            }
        }
    }
    return expr_copy(val);
}

/* Reduce val to a concrete-real double via N[Abs[val]].  Returns false
 * when the reduction does not collapse to a single Real / Integer /
 * Rational / MPFR / BigInt -- i.e. the value carries symbolic terms. */
static bool eigen_abs_to_double(Expr* val, double* out) {
    Expr* abs_e = eval_and_free(expr_new_function(
        expr_new_symbol("Abs"), (Expr*[]){ expr_copy(val) }, 1));
    Expr* n_abs = eval_and_free(expr_new_function(
        expr_new_symbol("N"), (Expr*[]){ abs_e }, 1));
    bool ok = true;
    double d = 0;
    if (n_abs->type == EXPR_REAL) d = n_abs->data.real;
    else if (n_abs->type == EXPR_INTEGER) d = (double)n_abs->data.integer;
    else if (n_abs->type == EXPR_BIGINT) d = mpz_get_d(n_abs->data.bigint);
    else if (n_abs->type == EXPR_MPFR)
        d = mpfr_get_d(n_abs->data.mpfr, MPFR_RNDN);
    else ok = false;
    expr_free(n_abs);
    *out = d;
    return ok;
}

typedef struct {
    Expr*  val;       /* borrowed during sort */
    double abs_d;
    size_t orig_idx;  /* tiebreaker to keep stable order */
} EigenSortKey;

static int eigen_sort_cmp_desc(const void* a, const void* b) {
    const EigenSortKey* ka = (const EigenSortKey*)a;
    const EigenSortKey* kb = (const EigenSortKey*)b;
    if (ka->abs_d > kb->abs_d) return -1;
    if (ka->abs_d < kb->abs_d) return 1;
    if (ka->orig_idx < kb->orig_idx) return -1;
    if (ka->orig_idx > kb->orig_idx) return 1;
    return 0;
}

/* Sort vals[] by descending |λ| if every entry reduces to a concrete real
 * via N[Abs[...]] (and is not Infinity/Indeterminate).  Otherwise leave
 * vals[] in Solve's natural order. */
void eigen_sort_by_abs_desc(Expr** vals, size_t n) {
    if (n <= 1) return;
    EigenSortKey* keys = malloc(sizeof(EigenSortKey) * n);
    bool all_numeric = true;
    for (size_t i = 0; i < n; i++) {
        keys[i].val = vals[i];
        keys[i].orig_idx = i;
        /* Infinity should sort first regardless. */
        if (vals[i]->type == EXPR_SYMBOL
            && vals[i]->data.symbol == SYM_Infinity) {
            keys[i].abs_d = INFINITY;
            continue;
        }
        if (vals[i]->type == EXPR_SYMBOL
            && vals[i]->data.symbol == SYM_Indeterminate) {
            keys[i].abs_d = -INFINITY; /* sort last */
            continue;
        }
        if (!eigen_abs_to_double(vals[i], &keys[i].abs_d)) {
            all_numeric = false;
            break;
        }
    }
    if (all_numeric) {
        qsort(keys, n, sizeof(EigenSortKey), eigen_sort_cmp_desc);
        for (size_t i = 0; i < n; i++) vals[i] = keys[i].val;
    }
    free(keys);
}

/* Compute null-space basis of n×n matrix M.  Uses RowReduce then constructs
 * a basis vector for each free column.  Returns an array of `*count_out`
 * Expr* vectors (each a List of n elements). */
Expr** eigen_null_space(Expr* M, int n, size_t* count_out) {
    *count_out = 0;
    Expr* rr_call = expr_new_function(expr_new_symbol("RowReduce"),
        (Expr*[]){ expr_copy(M) }, 1);
    Expr* R = eval_and_free(rr_call);
    if (!R || R->type != EXPR_FUNCTION
        || R->data.function.head->type != EXPR_SYMBOL
        || R->data.function.head->data.symbol != SYM_List
        || (int)R->data.function.arg_count != n) {
        if (R) expr_free(R);
        return NULL;
    }

    bool* pivot_col = calloc(n, sizeof(bool));
    int* row_pivot = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        row_pivot[i] = -1;
        Expr* row = R->data.function.args[i];
        for (int j = 0; j < n; j++) {
            if (!is_zero_poly(row->data.function.args[j])) {
                row_pivot[i] = j;
                pivot_col[j] = true;
                break;
            }
        }
    }

    Expr** basis = malloc(sizeof(Expr*) * n);
    size_t bc = 0;
    for (int free_col = 0; free_col < n; free_col++) {
        if (pivot_col[free_col]) continue;
        Expr** vec = malloc(sizeof(Expr*) * n);
        for (int k = 0; k < n; k++) vec[k] = expr_new_integer(0);
        expr_free(vec[free_col]);
        vec[free_col] = expr_new_integer(1);
        for (int i = 0; i < n; i++) {
            int pc = row_pivot[i];
            if (pc < 0) continue;
            Expr* r_val =
                R->data.function.args[i]->data.function.args[free_col];
            Expr* neg = eval_and_free(expr_new_function(
                expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(r_val) }, 2));
            expr_free(vec[pc]);
            vec[pc] = neg;
        }
        basis[bc++] = expr_new_function(expr_new_symbol("List"), vec, n);
        free(vec);
    }

    free(pivot_col);
    free(row_pivot);
    expr_free(R);
    *count_out = bc;
    return basis;
}

/* Detect inexact (Real / MPFR) leaves anywhere in a matrix. */
bool eigen_matrix_is_inexact(Expr* m) {
    if (!m) return false;
    if (m->type == EXPR_REAL) return true;
    if (m->type == EXPR_MPFR) return true;
    if (m->type != EXPR_FUNCTION) return false;
    if (eigen_matrix_is_inexact(m->data.function.head)) return true;
    for (size_t i = 0; i < m->data.function.arg_count; i++) {
        if (eigen_matrix_is_inexact(m->data.function.args[i])) return true;
    }
    return false;
}

/* Common option/positional argument parsing for Eigenvalues / Eigenvectors. */

/* True iff `e` is the symbol True. */
static bool eigen_is_true(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* Compare an EXPR_STRING value to a literal; returns true when they match. */
static bool eigen_string_eq(Expr* e, const char* lit) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, lit) == 0;
}

/* Public: parse the right-hand side of a Method -> <value> rule.
 * See mateigen.h for the accepted forms. */
MateigenMethod mateigen_parse_method_value(Expr* v) {
    if (!v) return MATEIGEN_AUTOMATIC;
    /* Bare Automatic symbol. */
    if (v->type == EXPR_SYMBOL && v->data.symbol == SYM_Automatic)
        return MATEIGEN_AUTOMATIC;
    /* Bare method-name symbol (Mathematica usually wraps in a string, but
     * tolerate the symbol form too). */
    if (v->type == EXPR_SYMBOL) {
        if (v->data.symbol == SYM_Direct)  return MATEIGEN_DIRECT;
        if (v->data.symbol == SYM_Arnoldi) return MATEIGEN_ARNOLDI;
        if (v->data.symbol == SYM_Banded)  return MATEIGEN_BANDED;
        if (v->data.symbol == SYM_FEAST)   return MATEIGEN_FEAST;
        return MATEIGEN_METHOD_UNKNOWN;
    }
    /* String form: "Direct", "Arnoldi", "Banded", "FEAST". */
    if (v->type == EXPR_STRING) {
        if (eigen_string_eq(v, "Direct"))  return MATEIGEN_DIRECT;
        if (eigen_string_eq(v, "Arnoldi")) return MATEIGEN_ARNOLDI;
        if (eigen_string_eq(v, "Banded"))  return MATEIGEN_BANDED;
        if (eigen_string_eq(v, "FEAST"))   return MATEIGEN_FEAST;
        return MATEIGEN_METHOD_UNKNOWN;
    }
    /* List form whose head element is the method name string, e.g.
     *   Method -> {"Arnoldi", "MaxIterations" -> 100}
     * Only the head element classifies the method here; sub-options are
     * interpreted later (Phase 3+). */
    if (v->type == EXPR_FUNCTION
        && v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol == SYM_List
        && v->data.function.arg_count >= 1) {
        return mateigen_parse_method_value(v->data.function.args[0]);
    }
    return MATEIGEN_METHOD_UNKNOWN;
}

/* Parse Eigenvalues/Eigenvectors arguments.  Returns false on shape error. */
bool eigen_parse_args(Expr* res, EigenOpts* opts) {
    if (!res || res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return false;

    opts->arg0 = res->data.function.args[0];
    opts->k_spec = NULL;
    /* Default Cubics/Quartics: True so radicals are emitted by default --
     * essential for the numerical path where Root[] objects cannot be
     * numericalised.  Spec lists False as the default Solve option;
     * Eigenvalues overrides to keep the closed-form pipeline functional. */
    opts->cubics_radical = true;
    opts->quartics_radical = true;
    opts->method_given = false;
    opts->method = MATEIGEN_AUTOMATIC;
    opts->method_value = NULL;

    /* Peel trailing options. */
    size_t pos_end = argc;
    while (pos_end > 1) {
        Expr* a = res->data.function.args[pos_end - 1];
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol == SYM_Rule
                || a->data.function.head->data.symbol == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol;
            Expr* rhs = a->data.function.args[1];
            if (name == SYM_Cubics) {
                opts->cubics_radical = eigen_is_true(rhs);
                pos_end--; continue;
            }
            if (name == SYM_Quartics) {
                opts->quartics_radical = eigen_is_true(rhs);
                pos_end--; continue;
            }
            if (name == SYM_Method) {
                opts->method_given = true;
                opts->method_value = rhs;
                opts->method = mateigen_parse_method_value(rhs);
                pos_end--; continue;
            }
        }
        break;
    }
    if (pos_end > 2) return false;
    if (pos_end == 2) opts->k_spec = res->data.function.args[1];
    return true;
}

/* ============================================================ *
 *  Phase 2: numerical "Direct" method                            *
 *                                                                 *
 *  Hand-rolled hot-path kernels that operate on flat dense        *
 *  numerical arrays so the inner loops never touch Expr nodes.    *
 *                                                                 *
 *  Each kernel has exactly one entry point and is annotated with  *
 *  a LAPACK-HOOK comment for future drop-in replacement under a   *
 *  USE_LAPACK build flag (see lovely-roaming-diffie.md).          *
 *                                                                 *
 *  Step 2.1 implementation (this file):                            *
 *    - Real symmetric matrices at machine precision (double).      *
 *    - Tridiagonalisation via Householder reflectors.              *
 *    - Symmetric tridiagonal QR with Wilkinson shift.              *
 *    - Both eigenvalues and orthonormal eigenvectors.              *
 *                                                                 *
 *  Subsequent commits add:                                         *
 *    - 2.2 real non-symmetric (Hessenberg + Francis QR)            *
 *    - 2.3 complex matrices                                        *
 *    - 2.4 MPFR (arbitrary-precision) kernels                      *
 * ============================================================ */

/* Machine-precision dense matrix workspace.  Row-major.  Always
 * heap-allocated; matD_free must be called regardless of how the
 * matrix was populated. */

void matD_free(MatD* M) {
    if (!M) return;
    free(M->re); free(M->im);
    M->re = M->im = NULL;
    M->n = 0; M->is_complex = 0;
}

/* Coerce an Expr leaf to a double.  Handles Integer, Real, BigInt,
 * MPFR, Rational[p, q], and Complex[re, im] (returning real part; the
 * imaginary part is returned via *out_im when non-NULL).  Returns
 * false if the leaf isn't a recognisable numeric value. */
static bool eigen_leaf_to_double(Expr* e, double* out_re, double* out_im) {
    if (out_im) *out_im = 0.0;
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out_re = e->data.real;                return true; }
    if (e->type == EXPR_INTEGER) { *out_re = (double)e->data.integer;     return true; }
    if (e->type == EXPR_BIGINT)  { *out_re = mpz_get_d(e->data.bigint);   return true; }
    if (e->type == EXPR_MPFR)    { *out_re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational && e->data.function.arg_count == 2) {
            double p, q;
            if (eigen_leaf_to_double(e->data.function.args[0], &p, NULL)
                && eigen_leaf_to_double(e->data.function.args[1], &q, NULL)
                && q != 0.0) {
                *out_re = p / q;
                return true;
            }
        }
        if (h == SYM_Complex && e->data.function.arg_count == 2) {
            double r, i;
            if (eigen_leaf_to_double(e->data.function.args[0], &r, NULL)
                && eigen_leaf_to_double(e->data.function.args[1], &i, NULL)) {
                *out_re = r;
                if (out_im) *out_im = i;
                return true;
            }
        }
    }
    return false;
}

/* True iff `e` is a numeric leaf that the dispatcher can convert to a
 * double (possibly with a non-zero imaginary part). */
static bool eigen_leaf_is_complex(Expr* e) {
    double r, i;
    return eigen_leaf_to_double(e, &r, &i) && i != 0.0;
}

/* Walk an n x n list-of-lists matrix expression and copy its numeric
 * entries into a freshly-allocated MatD.  Returns false (with *out
 * zeroed) if any entry isn't a numeric leaf -- the symbolic path
 * should handle that case. */
bool matD_load(Expr* m, size_t n, MatD* out) {
    out->n = n;
    out->is_complex = 0;
    out->re = (double*)malloc(sizeof(double) * n * n);
    out->im = NULL;
    /* First pass: detect any complex entries so we can allocate im[]. */
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            if (eigen_leaf_is_complex(row->data.function.args[j])) {
                out->is_complex = 1;
                break;
            }
        }
        if (out->is_complex) break;
    }
    if (out->is_complex) {
        out->im = (double*)calloc(n * n, sizeof(double));
    }
    /* Second pass: copy values. */
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            double r, im_v = 0.0;
            if (!eigen_leaf_to_double(row->data.function.args[j],
                                      &r, &im_v)) {
                matD_free(out);
                return false;
            }
            out->re[i * n + j] = r;
            if (out->is_complex) out->im[i * n + j] = im_v;
        }
    }
    return true;
}

/* Infinity-norm of an n x n real matrix. */
double matD_norm_inf_real(const double* A, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < n; j++) row_sum += fabs(A[i * n + j]);
        if (row_sum > m) m = row_sum;
    }
    return m;
}

/* True iff a real n x n matrix is symmetric to within `tol` (absolute
 * threshold; the caller multiplies by ||A||_inf and a small fudge). */
bool matD_is_real_symmetric(const MatD* A, double tol) {
    if (A->is_complex) return false;
    size_t n = A->n;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(A->re[i * n + j] - A->re[j * n + i]) > tol) return false;
        }
    }
    return true;
}

/* Infinity-norm of an n x n complex matrix: max over rows of the row
 * sum of the per-entry complex magnitudes |a_re + i a_im|. */
double matD_norm_inf_complex(const MatD* A) {
    size_t n = A->n;
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (size_t j = 0; j < n; j++) {
            row_sum += hypot(A->re[i * n + j],
                              A->is_complex ? A->im[i * n + j] : 0.0);
        }
        if (row_sum > m) m = row_sum;
    }
    return m;
}

/* True iff a complex n x n matrix is Hermitian to within `tol` (absolute
 * threshold on each entry; the caller multiplies by ||A||_inf and a
 * small fudge).  Hermitian = self-conjugate: A_ij = conj(A_ji), which
 * implies A_ii is purely real.  A purely real symmetric matrix passes
 * this test too. */
bool matD_is_hermitian(const MatD* A, double tol) {
    size_t n = A->n;
    if (!A->is_complex) return matD_is_real_symmetric(A, tol);
    /* Diagonal must be real. */
    for (size_t i = 0; i < n; i++) {
        if (fabs(A->im[i * n + i]) > tol) return false;
    }
    /* A_ij = conj(A_ji) for i != j. */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(A->re[i * n + j] - A->re[j * n + i]) > tol) return false;
            if (fabs(A->im[i * n + j] + A->im[j * n + i]) > tol) return false;
        }
    }
    return true;
}

/* Build a List of k_vectors lists of length n_components, each entry
 * Real or Complex[re,im].  Mirrors direct_build_complex_eigenvector_list
 * but allows rectangular shape (k_vectors != n_components) for the
 * Arnoldi case where mu < n eigenvectors have been recovered. */
Expr* mateigen_build_complex_eigvec_list_rect(const double* V_re,
                                                       const double* V_im,
                                                       size_t k_vectors,
                                                       size_t n_components) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * k_vectors);
    for (size_t k = 0; k < k_vectors; k++) {
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n_components);
        for (size_t i = 0; i < n_components; i++) {
            double r = V_re[k * n_components + i];
            double m = V_im[k * n_components + i];
            if (m == 0.0) {
                comps[i] = expr_new_real(r);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_real(r);
                args[1] = expr_new_real(m);
                comps[i] = expr_new_function(expr_new_symbol("Complex"), args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n_components);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, k_vectors);
    free(rows);
    return out;
}

/* Detect half-bandwidth of a real n*n matrix: the largest k such that
 * |A[i, j]| > tol for some |i - j| = k.  Returns n - 1 (dense) when
 * the corner entries exceed `tol`. */
size_t matD_bandwidth_real(const double* A, size_t n, double tol) {
    if (n <= 1) return 0;
    size_t b = 0;
    /* Walk diagonals from the outermost inward; first one with a
     * significant entry is the half-bandwidth. */
    for (size_t k = n - 1; k > 0; k--) {
        bool found = false;
        for (size_t i = 0; i + k < n; i++) {
            double v1 = A[i * n + (i + k)];
            double v2 = A[(i + k) * n + i];
            if (fabs(v1) > tol || fabs(v2) > tol) { found = true; break; }
        }
        if (found) { b = k; break; }
    }
    return b;
}

/* Half-bandwidth for a complex n*n matrix: any nonzero (re or im)
 * entry at distance k from the diagonal counts. */
size_t matD_bandwidth_complex(const MatD* A, double tol) {
    size_t n = A->n;
    if (n <= 1) return 0;
    size_t b = 0;
    for (size_t k = n - 1; k > 0; k--) {
        bool found = false;
        for (size_t i = 0; i + k < n; i++) {
            double r1 = A->re[i * n + (i + k)];
            double i1 = A->im[i * n + (i + k)];
            double r2 = A->re[(i + k) * n + i];
            double i2 = A->im[(i + k) * n + i];
            if (hypot(r1, i1) > tol || hypot(r2, i2) > tol) {
                found = true; break;
            }
        }
        if (found) { b = k; break; }
    }
    return b;
}
#ifdef USE_MPFR


/* ============================================================ *
 *  Phase 2 step 2d: MPFR (arbitrary-precision) Direct kernels.   *
 *                                                                 *
 *  Mirrors the machine-precision kernels above, but every double  *
 *  arithmetic op is replaced by an mpfr_* call at the matrix's    *
 *  combined precision.  Algorithm and control flow are identical; *
 *  only the working type and rounding discipline change.          *
 *                                                                 *
 *  Scratch discipline:                                            *
 *    - Each kernel pre-initialises its workspace once on entry    *
 *      and clears it once on exit; no mpfr_init / mpfr_clear in   *
 *      inner loops.                                               *
 *    - All rounding is MPFR_RNDN (round to nearest, ties to even).*
 *                                                                 *
 *  No LAPACK-HOOK annotations on this side -- LAPACK has no MPFR  *
 *  variant (MPLAPACK / MPACK is a separate effort, out of scope). *
 * ============================================================ */

/* Arbitrary-precision dense matrix workspace.  Row-major.  Each
 * mpfr_t in `re` / `im` is pre-initialised to `bits` on alloc and
 * cleared by matM_free.  Always heap-allocated. */

/* Free a heap-allocated mpfr_t[] of length `count`, clearing each cell. */
void mpfr_array_free(mpfr_t* a, size_t count) {
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* Allocate and pre-initialise an mpfr_t[] of length `count` to `bits`. */
mpfr_t* mpfr_array_alloc(size_t count, mpfr_prec_t bits) {
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

void matM_free(MatM* M) {
    if (!M) return;
    mpfr_array_free(M->re, M->n * M->n);
    mpfr_array_free(M->im, M->n * M->n);
    M->re = M->im = NULL;
    M->n = 0; M->is_complex = 0; M->bits = 0;
}

/* Load an Expr matrix into a freshly-allocated MatM at the given bit
 * precision.  Returns false (leaving *out zeroed) if any entry isn't a
 * numeric leaf -- the symbolic path should handle that case.  Sets
 * is_complex iff at least one entry has a non-zero imaginary part. */
bool matM_load(Expr* m, size_t n, mpfr_prec_t bits, MatM* out) {
    out->n = n;
    out->is_complex = 0;
    out->bits = bits;
    out->re = mpfr_array_alloc(n * n, bits);
    out->im = NULL;

    /* Probe for any imaginary content so we can allocate im[]. */
    for (size_t i = 0; i < n && !out->is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        for (size_t j = 0; j < n; j++) {
            Expr* cell = row->data.function.args[j];
            if (eigen_leaf_is_complex(cell)) { out->is_complex = 1; break; }
        }
    }
    if (out->is_complex) out->im = mpfr_array_alloc(n * n, bits);

    /* Fill cells via get_approx_mpfr (handles Integer / Rational /
     * Real / MPFR / Complex[...] uniformly at the target precision). */
    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (size_t j = 0; j < n; j++) {
            Expr* cell = row->data.function.args[j];
            bool is_inexact = false;
            if (!get_approx_mpfr(cell, out->re[i * n + j], tmp_im,
                                  &is_inexact)) {
                mpfr_clear(tmp_im);
                matM_free(out);
                return false;
            }
            if (out->is_complex) mpfr_set(out->im[i * n + j], tmp_im, MPFR_RNDN);
        }
    }
    mpfr_clear(tmp_im);
    return true;
}

/* Infinity-norm of a real n x n MPFR matrix.  Writes result to `out`. */
void matM_norm_inf_real(const mpfr_t* A, size_t n,
                                mpfr_prec_t bits, mpfr_t out) {
    mpfr_t row_sum, abs_a;
    mpfr_init2(row_sum, bits);
    mpfr_init2(abs_a,   bits);
    mpfr_set_zero(out, 1);
    for (size_t i = 0; i < n; i++) {
        mpfr_set_zero(row_sum, 1);
        for (size_t j = 0; j < n; j++) {
            mpfr_abs(abs_a, A[i * n + j], MPFR_RNDN);
            mpfr_add(row_sum, row_sum, abs_a, MPFR_RNDN);
        }
        if (mpfr_cmp(row_sum, out) > 0) mpfr_set(out, row_sum, MPFR_RNDN);
    }
    mpfr_clear(row_sum);
    mpfr_clear(abs_a);
}

/* True iff `A` is a real (no imag part) matrix that is symmetric to
 * within absolute tolerance `tol`. */
bool matM_is_real_symmetric(const MatM* A, const mpfr_t tol) {
    if (A->is_complex) return false;
    mpfr_t d;
    mpfr_init2(d, A->bits);
    for (size_t i = 0; i < A->n; i++) {
        for (size_t j = i + 1; j < A->n; j++) {
            mpfr_sub(d, A->re[i * A->n + j], A->re[j * A->n + i], MPFR_RNDN);
            mpfr_abs(d, d, MPFR_RNDN);
            if (mpfr_cmp(d, tol) > 0) { mpfr_clear(d); return false; }
        }
    }
    mpfr_clear(d);
    return true;
}

/* Infinity-norm of an n x n complex MPFR matrix.  Writes |A|_inf to `out`. */
void matM_norm_inf_complex(const MatM* A, mpfr_t out) {
    mpfr_t row_sum, mag;
    mpfr_init2(row_sum, A->bits);
    mpfr_init2(mag,     A->bits);
    mpfr_set_zero(out, 1);
    for (size_t i = 0; i < A->n; i++) {
        mpfr_set_zero(row_sum, 1);
        for (size_t j = 0; j < A->n; j++) {
            mpfr_hypot(mag, A->re[i * A->n + j],
                            A->im[i * A->n + j], MPFR_RNDN);
            mpfr_add(row_sum, row_sum, mag, MPFR_RNDN);
        }
        if (mpfr_cmp(row_sum, out) > 0) mpfr_set(out, row_sum, MPFR_RNDN);
    }
    mpfr_clear(row_sum);
    mpfr_clear(mag);
}

/* True iff `A` is Hermitian (A[i,j] == conj(A[j,i])) to within absolute
 * tolerance `tol`. */
bool matM_is_hermitian(const MatM* A, const mpfr_t tol) {
    if (!A->is_complex) return matM_is_real_symmetric(A, tol);
    mpfr_t d, e;
    mpfr_init2(d, A->bits);
    mpfr_init2(e, A->bits);
    for (size_t i = 0; i < A->n; i++) {
        /* Diagonal must be real (imag part within tol of zero). */
        mpfr_abs(d, A->im[i * A->n + i], MPFR_RNDN);
        if (mpfr_cmp(d, tol) > 0) { mpfr_clear(d); mpfr_clear(e); return false; }
        for (size_t j = i + 1; j < A->n; j++) {
            mpfr_sub(d, A->re[i * A->n + j], A->re[j * A->n + i], MPFR_RNDN);
            mpfr_abs(d, d, MPFR_RNDN);
            if (mpfr_cmp(d, tol) > 0) { mpfr_clear(d); mpfr_clear(e); return false; }
            mpfr_add(e, A->im[i * A->n + j], A->im[j * A->n + i], MPFR_RNDN);
            mpfr_abs(e, e, MPFR_RNDN);
            if (mpfr_cmp(e, tol) > 0) { mpfr_clear(d); mpfr_clear(e); return false; }
        }
    }
    mpfr_clear(d); mpfr_clear(e);
    return true;
}

/* ===================================================================
 * Phase 3f: "Arnoldi" Eigenvalues / Eigenvectors at MPFR precision.
 *
 * MPFR translation of the machine-precision Arnoldi kernels above.
 * Every double op becomes an mpfr_* call rounded MPFR_RNDN at the
 * matrix's combined precision.  The structure is identical:
 *   - Real general Arnoldi: paired arrays carry only real parts.
 *     H_m diagonalisation reuses direct_qr_real_general_M from 2d-B
 *     (and schur_compute_eigvecs_M for the eigenvector lift).
 *   - Complex general Arnoldi: paired re/im mpfr_t for V_m and H_m.
 *     The H_m diagonalisation builds a 2mu x 2mu real block embedding
 *     and routes through the same direct_qr_real_general_M / schur
 *     pipeline, then grouped complex Gram-Schmidt at MPFR precision
 *     deduplicates the doubled spectrum back to mu A-eigenpairs.
 *
 * Scratch discipline: all mpfr_t are pre-initialised on entry and
 * cleared on exit -- never inside the inner loops.
 * =================================================================== */

/* Rectangular variant of direct_build_complex_eigenvector_list_M:
 * emits a List of k_vectors lists of length n_components.  Used by
 * Arnoldi where mu < n eigenvectors have been recovered. */
Expr* mateigen_build_complex_eigvec_list_rect_M(const mpfr_t* V_re,
                                                         const mpfr_t* V_im,
                                                         size_t k_vectors,
                                                         size_t n_components) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * k_vectors);
    for (size_t k = 0; k < k_vectors; k++) {
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n_components);
        for (size_t i = 0; i < n_components; i++) {
            if (mpfr_zero_p(V_im[k * n_components + i])) {
                comps[i] = expr_new_mpfr_copy(V_re[k * n_components + i]);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_mpfr_copy(V_re[k * n_components + i]);
                args[1] = expr_new_mpfr_copy(V_im[k * n_components + i]);
                comps[i] = expr_new_function(expr_new_symbol("Complex"),
                                              args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n_components);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, k_vectors);
    free(rows);
    return out;
}

/* ============================================================ *
 *  Phase 4 step 4d: MPFR Banded kernels.                          *
 *                                                                 *
 *  Same algorithm as the machine kernels (Givens band reduction   *
 *  with bulge chasing) but every double op is an mpfr_*_RNDN.     *
 *  Scratch mpfr_t cells are pre-initialised and reused -- no      *
 *  init/clear inside the hot Givens loop.                         *
 *                                                                 *
 *  Bandwidth detection uses absolute tolerance = 2^{-bits+4}      *
 *  * n * ||A||_inf, matching direct_dispatch_mpfr's Hermitian-    *
 *  detection tolerance.                                           *
 * ============================================================ */

/* Half-bandwidth of a real MPFR matrix to within absolute tolerance
 * `tol`.  Returns 0 when the matrix is diagonal, n-1 when fully dense. */
size_t matM_bandwidth_real(const mpfr_t* A, size_t n,
                                    const mpfr_t tol) {
    if (n <= 1) return 0;
    size_t b = 0;
    mpfr_t a;
    mpfr_init2(a, mpfr_get_prec(tol));
    for (size_t k = n - 1; k > 0; k--) {
        bool found = false;
        for (size_t i = 0; i + k < n; i++) {
            mpfr_abs(a, A[i * n + (i + k)], MPFR_RNDN);
            if (mpfr_cmp(a, tol) > 0) { found = true; break; }
            mpfr_abs(a, A[(i + k) * n + i], MPFR_RNDN);
            if (mpfr_cmp(a, tol) > 0) { found = true; break; }
        }
        if (found) { b = k; break; }
    }
    mpfr_clear(a);
    return b;
}

/* Half-bandwidth of a complex MPFR matrix (any nonzero re/im at
 * distance > tol from the diagonal counts). */
size_t matM_bandwidth_complex(const MatM* A, const mpfr_t tol) {
    size_t n = A->n;
    if (n <= 1) return 0;
    size_t b = 0;
    mpfr_t mag;
    mpfr_init2(mag, A->bits);
    for (size_t k = n - 1; k > 0; k--) {
        bool found = false;
        for (size_t i = 0; i + k < n; i++) {
            mpfr_hypot(mag, A->re[i * n + (i + k)],
                            A->im[i * n + (i + k)], MPFR_RNDN);
            if (mpfr_cmp(mag, tol) > 0) { found = true; break; }
            mpfr_hypot(mag, A->re[(i + k) * n + i],
                            A->im[(i + k) * n + i], MPFR_RNDN);
            if (mpfr_cmp(mag, tol) > 0) { found = true; break; }
        }
        if (found) { b = k; break; }
    }
    mpfr_clear(mag);
    return b;
}
#endif /* USE_MPFR */
