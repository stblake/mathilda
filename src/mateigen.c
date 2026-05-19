#include "mateigen.h"
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
static const char* eigen_lambda_name(void) {
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
static bool eigen_extract_matrix_pair(Expr* arg, Expr** m_out, Expr** a_out,
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
static Expr* eigen_build_lambda_matrix(Expr* m, Expr* a_or_null,
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
static Expr* eigen_compute_det(Expr* matrix, int n) {
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
static Expr* eigen_char_poly_faddeev(Expr* A, const char* lambda_name, int n) {
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
static Expr* eigen_solve_poly(Expr* poly, const char* lambda_name,
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
static Expr** eigen_extract_values(Expr* solutions, size_t* count_out) {
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
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
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

static Expr* eigen_chop(Expr* val) {
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
#ifdef USE_MPFR
    else if (n_abs->type == EXPR_MPFR)
        d = mpfr_get_d(n_abs->data.mpfr, MPFR_RNDN);
#endif
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
static void eigen_sort_by_abs_desc(Expr** vals, size_t n) {
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
static Expr** eigen_null_space(Expr* M, int n, size_t* count_out) {
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
static bool eigen_matrix_is_inexact(Expr* m) {
    if (!m) return false;
    if (m->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (m->type == EXPR_MPFR) return true;
#endif
    if (m->type != EXPR_FUNCTION) return false;
    if (eigen_matrix_is_inexact(m->data.function.head)) return true;
    for (size_t i = 0; i < m->data.function.arg_count; i++) {
        if (eigen_matrix_is_inexact(m->data.function.args[i])) return true;
    }
    return false;
}

/* Common option/positional argument parsing for Eigenvalues / Eigenvectors. */
typedef struct {
    Expr* arg0;             /* m or {m, a}                          */
    Expr* k_spec;           /* Integer k, or UpTo[k], or NULL       */
    bool  cubics_radical;
    bool  quartics_radical;
    bool  method_given;     /* user supplied Method -> ...           */
    MateigenMethod method;  /* parsed Method, or AUTOMATIC if unset  */
    Expr* method_value;     /* original Method RHS (for sub-options) */
} EigenOpts;

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
static bool eigen_parse_args(Expr* res, EigenOpts* opts) {
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
typedef struct {
    size_t  n;
    int     is_complex;   /* 0: real only, 1: complex (im[] non-NULL) */
    double* re;           /* length n*n */
    double* im;           /* length n*n, NULL when !is_complex */
} MatD;

static void matD_free(MatD* M) {
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
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out_re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
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
static bool matD_load(Expr* m, size_t n, MatD* out) {
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
static double matD_norm_inf_real(const double* A, size_t n) {
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
static bool matD_is_real_symmetric(const MatD* A, double tol) {
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
static double matD_norm_inf_complex(const MatD* A) {
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
static bool matD_is_hermitian(const MatD* A, double tol) {
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

/* LAPACK-HOOK: replace with dsytrd when USE_LAPACK is set.
 *
 * Householder tridiagonalisation of a real symmetric n x n matrix `A`
 * (row-major, modified in place).  On return:
 *   diag[i]     = T_ii      for 0 <= i < n
 *   subdiag[i]  = T_{i+1,i} for 0 <= i < n-1
 *   Q (n*n)     = orthogonal accumulated reflectors, when `want_Q`.
 *
 * Algorithm: classic Householder reflectors applied symmetrically
 * via the rank-2 update A <- A - u q^T - q u^T where p = A u, K = u^T p / 2,
 * q = p - K u.  O(2 n^3 / 3) flops.  See Golub & Van Loan, Alg 8.3.1.
 *
 * Scratch buffers `u`, `p`, `q` (size n each) are caller-provided so
 * the inner loops never touch malloc/free.
 */
static void direct_tridiag_real_sym(double* A, size_t n,
                                     double* diag, double* subdiag,
                                     double* Q, bool want_Q,
                                     double* u, double* p, double* q) {
    if (want_Q) {
        /* Q starts at identity; reflectors are applied from the right
         * so the columns of Q become the orthogonal eigenvectors of A
         * once symmetric QR finishes. */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) Q[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* Compute Householder vector u for A[k+1:n, k]. */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double v = A[i * n + k];
            sigma += v * v;
        }
        double xk1 = A[(k + 1) * n + k];
        if (sigma == 0.0) {
            subdiag[k] = xk1;
            continue;
        }
        /* sigma already contains xk1*xk1 (the loop runs i = k+1..n-1, and
         * xk1 = A[k+1, k] is the first such entry).  The cancellation-safe
         * form: alpha = -sign(xk1) * ||x||. */
        double norm_x = sqrt(sigma);
        double alpha = (xk1 >= 0.0) ? -norm_x : norm_x;

        /* u[k+1] = xk1 - alpha; u[i>k+1] = A[i, k];  then normalise. */
        u[k + 1] = xk1 - alpha;
        for (size_t i = k + 2; i < n; i++) u[i] = A[i * n + k];
        double unorm2 = u[k + 1] * u[k + 1];
        for (size_t i = k + 2; i < n; i++) unorm2 += u[i] * u[i];
        if (unorm2 == 0.0) {
            subdiag[k] = xk1;
            continue;
        }
        double unorm = sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) u[i] /= unorm;

        /* p = A_22 u   (only on the trailing sub-block) */
        for (size_t i = k + 1; i < n; i++) {
            double s = 0.0;
            for (size_t j = k + 1; j < n; j++) s += A[i * n + j] * u[j];
            p[i] = 2.0 * s;
        }
        /* K = u^T A u = (u^T p) / 2 (since p = 2 A u). */
        double K = 0.0;
        for (size_t i = k + 1; i < n; i++) K += u[i] * p[i];
        K *= 0.5;
        /* q = p - 2 K u.  The rank-2 update A <- A - u q^T - q u^T then
         * equals H A H for H = I - 2 u u^T (Golub & Van Loan Alg 8.3.1). */
        for (size_t i = k + 1; i < n; i++) q[i] = p[i] - 2.0 * K * u[i];

        /* A_22 -= u q^T + q u^T */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                A[i * n + j] -= u[i] * q[j] + q[i] * u[j];
            }
        }

        /* Set the new subdiagonal element and clear the eliminated
         * column / row entries explicitly (drift control). */
        subdiag[k] = alpha;
        A[(k + 1) * n + k] = alpha;
        A[k * n + (k + 1)] = alpha;
        for (size_t i = k + 2; i < n; i++) {
            A[i * n + k] = 0.0;
            A[k * n + i] = 0.0;
        }

        /* Q <- Q * H_k   (right multiplication; column i replaced by
         * Q_i - 2 (Q row . u) u_i).  Apply only to Q[*, k+1..n-1]. */
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = k + 1; j < n; j++) s += Q[i * n + j] * u[j];
                s *= 2.0;
                for (size_t j = k + 1; j < n; j++) Q[i * n + j] -= s * u[j];
            }
        }
    }

    /* Extract diagonal. */
    for (size_t i = 0; i < n; i++) diag[i] = A[i * n + i];
    /* The n-2 -> n-1 subdiagonal is already in A. */
    if (n >= 2) subdiag[n - 2] = A[(n - 1) * n + (n - 2)];
}

/* Implicit-shift symmetric tridiagonal QR with Wilkinson shift.
 *
 * LAPACK-HOOK: replace with dsteqr (or dstedc) when USE_LAPACK is set.
 *
 *   diag[0..n-1]    : in-place diagonal of the tridiagonal matrix
 *   sub [0..n-2]    : in-place sub/super-diagonal
 *   Q   (n*n)       : in-place orthogonal eigenvector accumulator,
 *                     when `want_Q`.  Caller initialises Q (typically
 *                     to the orthogonal matrix from the tridiag step).
 *
 * Iterates over the active sub-block, deflating when |sub[i]| falls
 * below |diag[i]|+|diag[i+1]| * relative_tol.  Returns 0 on success,
 * -1 if the maximum number of sweeps (30*n) is exceeded.  Stagnation
 * is exceptionally rare for symmetric tridiagonal inputs but we cap
 * to avoid theoretical hangs.
 */
static int direct_symtridiag_qr(double* diag, double* sub, size_t n,
                                 double* Q, bool want_Q) {
    if (n == 0) return 0;
    const double rel_tol = 1e-14;   /* much tighter than chop threshold */
    const size_t max_sweeps = 30 * n;
    size_t sweeps = 0;

    size_t end = n;  /* active sub-block is [0..end-1]. */
    while (end > 1) {
        /* Find the largest m such that sub[m..end-2] are all "significant". */
        size_t m = end - 1;
        while (m > 0) {
            double tol = rel_tol * (fabs(diag[m - 1]) + fabs(diag[m]));
            if (fabs(sub[m - 1]) <= tol) { sub[m - 1] = 0.0; break; }
            m--;
        }
        if (m == end - 1) { end--; continue; }  /* deflated bottom */

        if (++sweeps > max_sweeps) return -1;

        /* Wilkinson shift on the trailing 2x2 block. */
        double d = (diag[end - 2] - diag[end - 1]) * 0.5;
        double e = sub[end - 2];
        double t = (d == 0.0) ? fabs(e)
                              : fabs(d) + sqrt(d * d + e * e);
        double sign_d = (d >= 0.0) ? 1.0 : -1.0;
        double mu = diag[end - 1] - sign_d * (e * e) / t;

        /* Implicit QR sweep on [m..end-1] using Givens rotations. */
        double x = diag[m] - mu;
        double z = sub[m];
        for (size_t k = m; k < end - 1; k++) {
            double c, s;
            double r = hypot(x, z);
            if (r == 0.0) { c = 1.0; s = 0.0; }
            else { c = x / r; s = z / r; }

            if (k > m) sub[k - 1] = r;

            double d_k    = diag[k];
            double d_k1   = diag[k + 1];
            double e_k    = sub[k];

            /* Two-sided Givens rotation Q^T A Q with Q = [c -s; s c]
             * (the rotation whose transpose annihilates z in [x; z]):
             *   d_k'   = c^2 d_k + 2 c s e_k + s^2 d_k1
             *   d_k1'  = s^2 d_k - 2 c s e_k + c^2 d_k1
             *   e_k'   = c s (d_k1 - d_k) + (c^2 - s^2) e_k
             */
            diag[k]     = c * c * d_k + 2.0 * c * s * e_k + s * s * d_k1;
            diag[k + 1] = s * s * d_k - 2.0 * c * s * e_k + c * c * d_k1;
            sub[k]      = c * s * (d_k1 - d_k) + (c * c - s * s) * e_k;

            /* Chase the bulge: rotation also affects sub[k+1] for k < end-2. */
            if (k + 1 < end - 1) {
                double t_next = sub[k + 1];
                x = sub[k];
                z = s * t_next;
                sub[k + 1] = c * t_next;
            }

            /* Update Q: post-multiply by the Givens rotation in columns
             * k and k+1.   Q_col_k  <- c Q_col_k + s Q_col_{k+1}
             *               Q_col_k1 <- -s Q_col_k + c Q_col_{k+1}    */
            if (want_Q) {
                for (size_t i = 0; i < n; i++) {
                    double qk  = Q[i * n + k];
                    double qk1 = Q[i * n + (k + 1)];
                    Q[i * n + k]       =  c * qk + s * qk1;
                    Q[i * n + (k + 1)] = -s * qk + c * qk1;
                }
            }
        }
    }
    return 0;
}

/* Sort eigenvalue indices into descending |lambda|, stable on ties.
 * Writes the permutation into `perm[0..n-1]`. */
static void direct_sort_perm_desc_abs(const double* vals, size_t n,
                                       size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    /* Insertion sort: n is small for the matrix sizes we care about
     * and the comparison is cheap. */
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        double ac = fabs(vals[cur]);
        size_t j = i;
        while (j > 0) {
            double ap = fabs(vals[perm[j - 1]]);
            if (ap > ac || (ap == ac && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
    }
}

/* Build the final List[Real,...] of eigenvalues from a sorted permutation. */
static Expr* direct_build_real_eigenvalue_list(const double* vals, size_t n,
                                                const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        items[i] = expr_new_real(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* Build the final List[List[Real,...], ...] of eigenvectors from a sorted
 * permutation.  Q is the n x n column-major-of-eigenvectors matrix in
 * row-major storage: Q[i, p] is the i-th component of the p-th
 * eigenvector. */
static Expr* direct_build_real_eigenvector_list(const double* Q, size_t n,
                                                 const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            comps[i] = expr_new_real(Q[i * n + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Apply k-spec selection to an already-sorted (descending |lambda|)
 * result.  Returns a fresh List with the trimmed entries.  Mirrors
 * eigen_apply_k_spec but for Expr trees we already own. */
static Expr* direct_apply_k_spec_list(Expr* full_list, Expr* k_spec) {
    if (!k_spec) return full_list;
    size_t count = full_list->data.function.arg_count;
    size_t result_count = count;
    bool from_end = false;
    if (k_spec->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.integer;
        if (k >= 0) {
            result_count = ((size_t)k < count) ? (size_t)k : count;
        } else {
            int64_t abs_k = -k;
            result_count = ((size_t)abs_k < count) ? (size_t)abs_k : count;
            from_end = true;
        }
    } else if (k_spec->type == EXPR_FUNCTION
        && k_spec->data.function.head->type == EXPR_SYMBOL
        && k_spec->data.function.head->data.symbol == SYM_UpTo
        && k_spec->data.function.arg_count == 1
        && k_spec->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.function.args[0]->data.integer;
        result_count = ((size_t)k < count) ? (size_t)k : count;
    }

    Expr** items = result_count
        ? (Expr**)malloc(sizeof(Expr*) * result_count) : NULL;
    if (from_end) {
        size_t start = count - result_count;
        for (size_t i = 0; i < start; i++)
            expr_free(full_list->data.function.args[i]);
        for (size_t i = 0; i < result_count; i++)
            items[i] = full_list->data.function.args[start + i];
    } else {
        for (size_t i = 0; i < result_count; i++)
            items[i] = full_list->data.function.args[i];
        for (size_t i = result_count; i < count; i++)
            expr_free(full_list->data.function.args[i]);
    }
    free(full_list->data.function.args);
    full_list->data.function.args = items;
    full_list->data.function.arg_count = result_count;
    return full_list;
}

/* ---------- Real non-symmetric Direct: Hessenberg + Francis QR ------ *
 *                                                                       *
 *  For a general real n x n matrix, eigenvalues are extracted by:       *
 *    1. Householder reduction to upper Hessenberg form.                 *
 *    2. Implicit double-shift QR (Francis) sweeps over the Hessenberg,  *
 *       deflating off-diagonal entries until the matrix is quasi-       *
 *       triangular (Schur form: 1x1 diagonal blocks for real            *
 *       eigenvalues, 2x2 blocks for complex conjugate pairs).           *
 *    3. Read eigenvalues off the Schur form.                            *
 *                                                                       *
 *  Eigenvectors require Q accumulation through both reduction and QR    *
 *  sweeps and arrive in step 2b (this commit covers eigenvalues only).  *
 *                                                                       *
 *  LAPACK-HOOK: dgehrd for step 1, dhseqr for step 2, dtrevc3 +         *
 *  dorghr for the eigenvector back-transform in step 2b.                *
 * --------------------------------------------------------------------- */

/* Householder reduction of an n x n matrix `A` (row-major, modified in
 * place) to upper Hessenberg form.  When `Q` is non-NULL the caller
 * passes a pre-initialised n x n identity matrix that this routine
 * post-multiplies by each Householder reflector, producing the
 * orthogonal back-transformation Q with Q^T A_in Q = A_out + (Schur).
 *
 * Scratch buffer `u` (size n) is caller-provided.
 *
 * LAPACK-HOOK: dgehrd (+ dorghr to materialise Q from the reflectors). */
static void direct_hessenberg_real(double* A, size_t n, double* u, double* Q) {
    for (size_t k = 0; k + 2 < n; k++) {
        /* Householder vector for column k below the sub-diagonal. */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double v = A[i * n + k];
            sigma += v * v;
        }
        if (sigma == 0.0) continue;
        double xk1 = A[(k + 1) * n + k];
        double norm_x = sqrt(sigma);
        double alpha = (xk1 >= 0.0) ? -norm_x : norm_x;
        u[k + 1] = xk1 - alpha;
        for (size_t i = k + 2; i < n; i++) u[i] = A[i * n + k];
        double unorm2 = u[k + 1] * u[k + 1];
        for (size_t i = k + 2; i < n; i++) unorm2 += u[i] * u[i];
        if (unorm2 == 0.0) continue;
        double unorm = sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) u[i] /= unorm;

        /* Left multiply: A[k+1..n-1, :] <- H * A[k+1..n-1, :] */
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t i = k + 1; i < n; i++) s += u[i] * A[i * n + j];
            s *= 2.0;
            for (size_t i = k + 1; i < n; i++) A[i * n + j] -= s * u[i];
        }
        /* Right multiply: A[:, k+1..n-1] <- A[:, k+1..n-1] * H */
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t j = k + 1; j < n; j++) s += A[i * n + j] * u[j];
            s *= 2.0;
            for (size_t j = k + 1; j < n; j++) A[i * n + j] -= s * u[j];
        }

        /* Q <- Q * H_k (right-multiply by reflector on cols k+1..n-1). */
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = k + 1; j < n; j++) s += Q[i * n + j] * u[j];
                s *= 2.0;
                for (size_t j = k + 1; j < n; j++) Q[i * n + j] -= s * u[j];
            }
        }

        /* Tidy up: explicit subdiag and zero below it. */
        A[(k + 1) * n + k] = alpha;
        for (size_t i = k + 2; i < n; i++) A[i * n + k] = 0.0;
    }
}

/* One Francis double-shift QR step applied to the active sub-block
 * H[q..p-1, q..p-1] of an upper Hessenberg matrix.  The shift pair is
 * derived from the eigenvalues of the trailing 2x2 block H[p-2:p-1,
 * p-2:p-1] (real or complex; real arithmetic throughout since the
 * shift pair is a conjugate pair when complex).
 *
 * Bulge chasing uses Householder reflectors of size 3 (size 2 at the
 * very end of the chase).  Sub-block size p - q must be >= 3.
 *
 * When `Q` is non-NULL each bulge-chase reflector also right-multiplies
 * Q so the caller can extract eigenvectors from the final Schur form.
 *
 * LAPACK-HOOK: this is the inner loop body of dhseqr / dlahqr. */
static void direct_francis_step(double* H, size_t n, size_t q, size_t p,
                                 double* Q) {
    double h11 = H[(p - 2) * n + (p - 2)];
    double h12 = H[(p - 2) * n + (p - 1)];
    double h21 = H[(p - 1) * n + (p - 2)];
    double h22 = H[(p - 1) * n + (p - 1)];
    double s = h11 + h22;                       /* trace of trailing 2x2 */
    double t = h11 * h22 - h12 * h21;           /* det of trailing 2x2   */

    /* First three entries of M's first column for the active block. */
    double g11 = H[q * n + q];
    double g12 = H[q * n + (q + 1)];
    double g21 = H[(q + 1) * n + q];
    double g22 = H[(q + 1) * n + (q + 1)];
    double g32 = H[(q + 2) * n + (q + 1)];
    double x = g11 * g11 + g12 * g21 - s * g11 + t;
    double y = g21 * (g11 + g22 - s);
    double z = g21 * g32;

    for (size_t k = q; k + 2 < p; k++) {
        /* 3-element Householder for (x, y, z). */
        double sig = x * x + y * y + z * z;
        if (sig == 0.0) {
            /* Skip degenerate step; advance and continue. */
            x = H[(k + 1) * n + k];
            y = H[(k + 2) * n + k];
            z = (k + 3 < p) ? H[(k + 3) * n + k] : 0.0;
            continue;
        }
        double norm_v = sqrt(sig);
        double a = (x >= 0.0) ? -norm_v : norm_v;
        double u1 = x - a;
        double u2 = y;
        double u3 = z;
        double un2 = u1 * u1 + u2 * u2 + u3 * u3;
        double un = sqrt(un2);
        u1 /= un; u2 /= un; u3 /= un;

        /* Apply P from the left to rows (k, k+1, k+2). */
        size_t col_start = (k > q) ? k - 1 : q;
        for (size_t j = col_start; j < n; j++) {
            double s0 = (u1 * H[k * n + j]
                       + u2 * H[(k + 1) * n + j]
                       + u3 * H[(k + 2) * n + j]) * 2.0;
            H[k * n + j]       -= s0 * u1;
            H[(k + 1) * n + j] -= s0 * u2;
            H[(k + 2) * n + j] -= s0 * u3;
        }
        /* Apply P from the right to cols (k, k+1, k+2). */
        size_t row_end = (k + 3 < p) ? (k + 3) : (p - 1);
        for (size_t i = 0; i <= row_end; i++) {
            double s0 = (H[i * n + k]       * u1
                       + H[i * n + (k + 1)] * u2
                       + H[i * n + (k + 2)] * u3) * 2.0;
            H[i * n + k]       -= s0 * u1;
            H[i * n + (k + 1)] -= s0 * u2;
            H[i * n + (k + 2)] -= s0 * u3;
        }
        /* Q <- Q * P (right-multiply by reflector on cols k, k+1, k+2). */
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s0 = (Q[i * n + k]       * u1
                           + Q[i * n + (k + 1)] * u2
                           + Q[i * n + (k + 2)] * u3) * 2.0;
                Q[i * n + k]       -= s0 * u1;
                Q[i * n + (k + 1)] -= s0 * u2;
                Q[i * n + (k + 2)] -= s0 * u3;
            }
        }

        x = H[(k + 1) * n + k];
        y = H[(k + 2) * n + k];
        z = (k + 3 < p) ? H[(k + 3) * n + k] : 0.0;
    }

    /* Final 2-element Householder on (x, y) at rows (p-2, p-1). */
    {
        size_t k = p - 2;
        double norm_v = sqrt(x * x + y * y);
        if (norm_v == 0.0) return;
        double a = (x >= 0.0) ? -norm_v : norm_v;
        double u1 = x - a;
        double u2 = y;
        double un2 = u1 * u1 + u2 * u2;
        double un = sqrt(un2);
        u1 /= un; u2 /= un;
        size_t col_start = k - 1;
        for (size_t j = col_start; j < n; j++) {
            double s0 = (u1 * H[k * n + j]
                       + u2 * H[(k + 1) * n + j]) * 2.0;
            H[k * n + j]       -= s0 * u1;
            H[(k + 1) * n + j] -= s0 * u2;
        }
        for (size_t i = 0; i < p; i++) {
            double s0 = (H[i * n + k] * u1
                       + H[i * n + (k + 1)] * u2) * 2.0;
            H[i * n + k]       -= s0 * u1;
            H[i * n + (k + 1)] -= s0 * u2;
        }
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s0 = (Q[i * n + k]       * u1
                           + Q[i * n + (k + 1)] * u2) * 2.0;
                Q[i * n + k]       -= s0 * u1;
                Q[i * n + (k + 1)] -= s0 * u2;
            }
        }
    }
}

/* When a trailing 2x2 block has REAL eigenvalues we need to triangularise
 * it (drive the sub-diagonal to zero) before deflating, otherwise the
 * Schur form is non-triangular and back-substitution returns the wrong
 * eigenvector.  Apply a Givens rotation whose first column is the
 * eigenvector for one of the two real eigenvalues; after the similarity
 * the (i+1, i) entry is zero and the two real eigenvalues sit on the
 * diagonal.  Caller deflates them as two 1x1 blocks on the next QR
 * iteration.
 *
 * The block is at rows / cols (p-2, p-1).  Q is right-multiplied when
 * non-NULL so the caller can recover eigenvectors via back-transform. */
static void direct_split_2x2_real(double* H, size_t n, size_t p, double* Q) {
    size_t i = p - 2;
    double a = H[i * n + i];
    double b = H[i * n + (i + 1)];
    double c = H[(i + 1) * n + i];
    double d = H[(i + 1) * n + (i + 1)];
    double tr = a + d;
    double det = a * d - b * c;
    double disc = tr * tr - 4.0 * det;
    if (disc < 0.0) return;             /* complex; caller handles it */
    double sq = sqrt(disc);
    double lam = (tr + sq) * 0.5;       /* larger eigenvalue */
    double v0 = lam - d;
    double v1 = c;
    double r = sqrt(v0 * v0 + v1 * v1);
    if (r == 0.0) return;
    double cs = v0 / r;
    double sn = v1 / r;
    /* Left-multiply rows (i, i+1) of H by G^T = [[cs, sn], [-sn, cs]]. */
    for (size_t j = 0; j < n; j++) {
        double r0 = H[i * n + j];
        double r1 = H[(i + 1) * n + j];
        H[i * n + j]       =  cs * r0 + sn * r1;
        H[(i + 1) * n + j] = -sn * r0 + cs * r1;
    }
    /* Right-multiply cols (i, i+1) of H by G = [[cs, -sn], [sn, cs]]. */
    for (size_t k = 0; k < n; k++) {
        double c0 = H[k * n + i];
        double c1 = H[k * n + (i + 1)];
        H[k * n + i]       =  cs * c0 + sn * c1;
        H[k * n + (i + 1)] = -sn * c0 + cs * c1;
    }
    if (Q) {
        for (size_t k = 0; k < n; k++) {
            double c0 = Q[k * n + i];
            double c1 = Q[k * n + (i + 1)];
            Q[k * n + i]       =  cs * c0 + sn * c1;
            Q[k * n + (i + 1)] = -sn * c0 + cs * c1;
        }
    }
    H[(i + 1) * n + i] = 0.0;
}

/* Drive Francis QR sweeps on an upper Hessenberg matrix `H` until it
 * reaches quasi-triangular (Schur) form.  Eigenvalues are written into
 * eval_re / eval_im in Schur position order (NOT in sort order -- the
 * caller is responsible for sorting via direct_sort_perm_desc_abs_complex).
 *
 * When `Q` is non-NULL each bulge-chase reflector right-multiplies Q so
 * the caller can extract eigenvectors from the final Schur form via
 * back-substitution + back-transformation.
 *
 * Returns 0 on success, -1 if the maximum number of total Francis
 * sweeps (30 * n) is exceeded without convergence.
 *
 * LAPACK-HOOK: dhseqr. */
static int direct_qr_real_general(double* H, size_t n,
                                    double* eval_re, double* eval_im,
                                    double* Q) {
    const double eps = 1e-14;
    const size_t max_iter = 30 * n;
    size_t iter = 0;
    size_t p = n;

    while (p > 0) {
        if (p == 1) {
            eval_re[0] = H[0];
            eval_im[0] = 0.0;
            break;
        }
        if (iter++ > max_iter) return -1;

        /* Find largest q in [1..p-1] such that H[q, q-1] is negligible
         * (deflation point).  q = 0 means the whole [0..p-1] block is
         * unreduced. */
        size_t q = p - 1;
        while (q > 0) {
            double tol = eps * (fabs(H[(q - 1) * n + (q - 1)])
                              + fabs(H[q * n + q]));
            if (fabs(H[q * n + (q - 1)]) <= tol) {
                H[q * n + (q - 1)] = 0.0;
                break;
            }
            q -= 1;
        }

        if (q == p - 1) {
            /* Trailing 1x1 block deflated -- real eigenvalue. */
            eval_re[p - 1] = H[(p - 1) * n + (p - 1)];
            eval_im[p - 1] = 0.0;
            p -= 1;
            iter = 0;
            continue;
        }
        if (q == p - 2) {
            /* Trailing 2x2 block. */
            double a = H[(p - 2) * n + (p - 2)];
            double b = H[(p - 2) * n + (p - 1)];
            double c = H[(p - 1) * n + (p - 2)];
            double d = H[(p - 1) * n + (p - 1)];
            double tr = a + d;
            double det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc < 0.0) {
                /* Complex conjugate pair: deflate the entire 2x2 block. */
                double sq = sqrt(-disc);
                eval_re[p - 2] = tr * 0.5;
                eval_im[p - 2] =  sq * 0.5;
                eval_re[p - 1] = tr * 0.5;
                eval_im[p - 1] = -sq * 0.5;
                p -= 2;
                iter = 0;
                continue;
            }
            /* Real eigenvalues.  Store the analytic values (more
             * accurate than re-reading after the Givens roundoff),
             * then triangularise the block via a Givens similarity
             * so the Schur form is properly upper triangular and
             * back-substitution can recover eigenvectors. */
            double sq2 = sqrt(disc);
            eval_re[p - 2] = (tr + sq2) * 0.5;
            eval_im[p - 2] = 0.0;
            eval_re[p - 1] = (tr - sq2) * 0.5;
            eval_im[p - 1] = 0.0;
            direct_split_2x2_real(H, n, p, Q);
            /* After the split H[p-2, p-2] = larger eigenvalue and
             * H[p-1, p-1] = smaller, with H[p-1, p-2] = 0.  Overwrite
             * the diagonal entries with the analytic eigenvalues so
             * back-substitution uses values matching eval_re. */
            H[(p - 2) * n + (p - 2)] = eval_re[p - 2];
            H[(p - 1) * n + (p - 1)] = eval_re[p - 1];
            p -= 2;
            iter = 0;
            continue;
        }

        direct_francis_step(H, n, q, p, Q);
    }
    return 0;
}

/* Sort permutation by descending |lambda| (stable) where each lambda
 * is given by its real and imaginary parts. */
static void direct_sort_perm_desc_abs_complex(const double* re,
                                                const double* im,
                                                size_t n, size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        double ac = hypot(re[cur], im[cur]);
        size_t j = i;
        while (j > 0) {
            size_t prev = perm[j - 1];
            double ap = hypot(re[prev], im[prev]);
            if (ap > ac || (ap == ac && prev < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
    }
}

/* Build a List of eigenvalues from real / imaginary parts in the order
 * given by `perm`.  Real eigenvalues become EXPR_REAL; complex pairs
 * become Complex[re, im]. */
static Expr* direct_build_complex_eigenvalue_list(const double* re,
                                                    const double* im,
                                                    size_t n,
                                                    const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t idx = perm[i];
        if (im[idx] == 0.0) {
            items[i] = expr_new_real(re[idx]);
        } else {
            Expr** comp_args = (Expr**)malloc(sizeof(Expr*) * 2);
            comp_args[0] = expr_new_real(re[idx]);
            comp_args[1] = expr_new_real(im[idx]);
            items[i] = expr_new_function(expr_new_symbol("Complex"),
                                          comp_args, 2);
            free(comp_args);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* ---- Eigenvectors from Schur form ----------------------------------- *
 *                                                                       *
 *  After Hessenberg + Francis QR we have H = Q^T A Q in quasi-          *
 *  triangular Schur form (1x1 blocks for real eigenvalues, 2x2 blocks   *
 *  for complex conjugate pairs).  Eigenvectors of A are obtained by:    *
 *    1. Solving (H - lambda I) v_schur = 0 by back-substitution in the  *
 *       Schur basis (real and complex variants below).                  *
 *    2. Back-transforming: v_A = Q v_schur.                             *
 *    3. Normalising to unit 2-norm.                                     *
 *                                                                       *
 *  LAPACK-HOOK: dtrevc3 handles steps 1-2 in one call when USE_LAPACK   *
 *  is wired (dorghr materialises Q from the Hessenberg reflectors).     *
 * --------------------------------------------------------------------- */

/* Back-substitute for the eigenvector of a Schur quasi-triangular matrix
 * H corresponding to a real eigenvalue at Schur position k.  Writes the
 * Schur-basis eigenvector into v (length n; v[i] = 0 for i > k). */
static void schur_eigvec_real(const double* H, size_t n, size_t k,
                                double lambda, double* v) {
    for (size_t i = 0; i < n; i++) v[i] = 0.0;
    v[k] = 1.0;
    if (k == 0) return;

    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && (H[i * n + (i - 1)] != 0.0);
        if (is_2x2) {
            /* Solve 2x2 system coupling v[i-1] and v[i] (the rows of the
             * 2x2 block).  Right-hand side = -sum of already-known
             * components weighted by the off-diagonal entries. */
            double rhs1 = 0.0, rhs2 = 0.0;
            for (size_t j = i + 1; j <= k; j++) {
                rhs1 += H[(i - 1) * n + j] * v[j];
                rhs2 += H[i * n + j] * v[j];
            }
            rhs1 = -rhs1; rhs2 = -rhs2;
            double a = H[(i - 1) * n + (i - 1)] - lambda;
            double b = H[(i - 1) * n + i];
            double c = H[i * n + (i - 1)];
            double d = H[i * n + i] - lambda;
            double det = a * d - b * c;
            if (fabs(det) < 1e-300) {
                v[i - 1] = 0.0; v[i] = 0.0;
            } else {
                v[i - 1] = (d * rhs1 - b * rhs2) / det;
                v[i]     = (a * rhs2 - c * rhs1) / det;
            }
            i--;                            /* skip the row we just used */
        } else {
            double rhs = 0.0;
            for (size_t j = i + 1; j <= k; j++) rhs += H[i * n + j] * v[j];
            double diag = H[i * n + i] - lambda;
            if (fabs(diag) < 1e-300) {
                v[i] = 0.0;                 /* defective: leave zero */
            } else {
                v[i] = -rhs / diag;
            }
        }
    }
}

/* Back-substitute for the complex eigenvector of a Schur quasi-triangular
 * H corresponding to the complex eigenvalue lambda = a + i b (b > 0) at
 * Schur positions (k, k+1).  Writes real / imaginary parts of the Schur-
 * basis eigenvector into v_re / v_im (length n; both 0 above k+1). */
static void schur_eigvec_complex(const double* H, size_t n, size_t k,
                                   double a, double b,
                                   double* v_re, double* v_im) {
    for (size_t i = 0; i < n; i++) { v_re[i] = 0.0; v_im[i] = 0.0; }
    /* Initial values from the 2x2 block at (k, k+1):
     *   2x2 = [[alpha, beta], [gamma, delta]], eigenvalues lambda, conj(lambda)
     *   eigenvector for lambda: x_1 = 1, x_0 = (lambda - delta) / gamma.
     *   lambda - delta = (a - delta) + i b. */
    double delta = H[(k + 1) * n + (k + 1)];
    double gamma = H[(k + 1) * n + k];
    if (gamma == 0.0) {
        /* Shouldn't happen for an unreduced 2x2 block; bail. */
        v_re[k] = 1.0;
        v_re[k + 1] = 0.0;
        return;
    }
    v_re[k]     = (a - delta) / gamma;
    v_im[k]     = b / gamma;
    v_re[k + 1] = 1.0;
    v_im[k + 1] = 0.0;

    if (k == 0) return;
    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && (H[i * n + (i - 1)] != 0.0);
        if (is_2x2) {
            /* Two complex equations -> 2x2 complex system in
             *   (v[i-1], v[i]).  Solve via complex 2x2 inverse. */
            double rhs_u_top = 0.0, rhs_v_top = 0.0;
            double rhs_u_bot = 0.0, rhs_v_bot = 0.0;
            for (size_t j = i + 1; j <= k + 1; j++) {
                rhs_u_top += H[(i - 1) * n + j] * v_re[j];
                rhs_v_top += H[(i - 1) * n + j] * v_im[j];
                rhs_u_bot += H[i * n + j] * v_re[j];
                rhs_v_bot += H[i * n + j] * v_im[j];
            }
            /* Coeff matrix (complex):
             *   [(a11 - lambda)  a12 ]
             *   [   a21         (a22 - lambda)]
             * = [[A_re + i A_im, B_re + i B_im],
             *    [C_re + i C_im, D_re + i D_im]] */
            double A_re = H[(i - 1) * n + (i - 1)] - a, A_im = -b;
            double B_re = H[(i - 1) * n + i],            B_im = 0.0;
            double C_re = H[i * n + (i - 1)],            C_im = 0.0;
            double D_re = H[i * n + i] - a,              D_im = -b;
            /* det = A D - B C  (complex) */
            double det_re = (A_re * D_re - A_im * D_im)
                          - (B_re * C_re - B_im * C_im);
            double det_im = (A_re * D_im + A_im * D_re)
                          - (B_re * C_im + B_im * C_re);
            double det_mag2 = det_re * det_re + det_im * det_im;
            if (det_mag2 < 1e-300) {
                v_re[i - 1] = 0.0; v_im[i - 1] = 0.0;
                v_re[i] = 0.0; v_im[i] = 0.0;
            } else {
                /* RHS (negated). */
                double r1_re = -rhs_u_top, r1_im = -rhs_v_top;
                double r2_re = -rhs_u_bot, r2_im = -rhs_v_bot;
                /* v[i-1] = (D r1 - B r2) / det */
                double num1_re = (D_re * r1_re - D_im * r1_im)
                               - (B_re * r2_re - B_im * r2_im);
                double num1_im = (D_re * r1_im + D_im * r1_re)
                               - (B_re * r2_im + B_im * r2_re);
                /* v[i]   = (A r2 - C r1) / det */
                double num2_re = (A_re * r2_re - A_im * r2_im)
                               - (C_re * r1_re - C_im * r1_im);
                double num2_im = (A_re * r2_im + A_im * r2_re)
                               - (C_re * r1_im + C_im * r1_re);
                /* (p + iq) / (det_re + i det_im) = (p det_re + q det_im
                 *   + i (q det_re - p det_im)) / det_mag2 */
                v_re[i - 1] = (num1_re * det_re + num1_im * det_im) / det_mag2;
                v_im[i - 1] = (num1_im * det_re - num1_re * det_im) / det_mag2;
                v_re[i]     = (num2_re * det_re + num2_im * det_im) / det_mag2;
                v_im[i]     = (num2_im * det_re - num2_re * det_im) / det_mag2;
            }
            i--;
        } else {
            /* Single complex equation: (H[i,i] - lambda) v[i] = -rhs. */
            double rhs_u = 0.0, rhs_v = 0.0;
            for (size_t j = i + 1; j <= k + 1; j++) {
                rhs_u += H[i * n + j] * v_re[j];
                rhs_v += H[i * n + j] * v_im[j];
            }
            double diag_re = H[i * n + i] - a;
            double diag_im = -b;
            double denom = diag_re * diag_re + diag_im * diag_im;
            if (denom < 1e-300) {
                v_re[i] = 0.0; v_im[i] = 0.0;
            } else {
                double num_re = -rhs_u, num_im = -rhs_v;
                v_re[i] = (num_re * diag_re + num_im * diag_im) / denom;
                v_im[i] = (num_im * diag_re - num_re * diag_im) / denom;
            }
        }
    }
}

/* Compute eigenvectors of A from the Schur form H and back-transformation
 * Q, in sorted (descending |lambda|) order.  V_re / V_im are n x n
 * row-major arrays where row k = sorted-k-th eigenvector. */
static void schur_compute_eigvecs(const double* H, const double* Q,
                                    size_t n,
                                    const double* eval_re, const double* eval_im,
                                    const size_t* perm,
                                    double* V_re, double* V_im) {
    double* v_schur_re = (double*)malloc(sizeof(double) * n);
    double* v_schur_im = (double*)malloc(sizeof(double) * n);
    double* w_re       = (double*)malloc(sizeof(double) * n);
    double* w_im       = (double*)malloc(sizeof(double) * n);

    size_t* inv_perm = (size_t*)malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; i++) inv_perm[perm[i]] = i;

    size_t k = 0;
    while (k < n) {
        if (eval_im[k] == 0.0) {
            /* Real eigenvalue at Schur position k. */
            schur_eigvec_real(H, n, k, eval_re[k], v_schur_re);
            /* w = Q . v_schur. */
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = 0; j <= k; j++) s += Q[i * n + j] * v_schur_re[j];
                w_re[i] = s;
            }
            double norm2 = 0.0;
            for (size_t i = 0; i < n; i++) norm2 += w_re[i] * w_re[i];
            double inv = (norm2 > 0.0) ? 1.0 / sqrt(norm2) : 1.0;
            size_t sp = inv_perm[k];
            for (size_t i = 0; i < n; i++) {
                V_re[sp * n + i] = w_re[i] * inv;
                V_im[sp * n + i] = 0.0;
            }
            k++;
        } else {
            /* Complex pair at Schur positions k, k+1.  eval_im[k] is the
             * "+imag" root (the QR loop writes them in this order). */
            double a = eval_re[k];
            double b = fabs(eval_im[k]);
            schur_eigvec_complex(H, n, k, a, b, v_schur_re, v_schur_im);
            for (size_t i = 0; i < n; i++) {
                double s_re = 0.0, s_im = 0.0;
                for (size_t j = 0; j <= k + 1; j++) {
                    s_re += Q[i * n + j] * v_schur_re[j];
                    s_im += Q[i * n + j] * v_schur_im[j];
                }
                w_re[i] = s_re;
                w_im[i] = s_im;
            }
            double norm2 = 0.0;
            for (size_t i = 0; i < n; i++) {
                norm2 += w_re[i] * w_re[i] + w_im[i] * w_im[i];
            }
            double inv = (norm2 > 0.0) ? 1.0 / sqrt(norm2) : 1.0;
            size_t sp1 = inv_perm[k];       /* +imag eigenvalue */
            size_t sp2 = inv_perm[k + 1];   /* -imag eigenvalue */
            for (size_t i = 0; i < n; i++) {
                double r = w_re[i] * inv;
                double m = w_im[i] * inv;
                V_re[sp1 * n + i] = r;  V_im[sp1 * n + i] =  m;
                V_re[sp2 * n + i] = r;  V_im[sp2 * n + i] = -m;
            }
            k += 2;
        }
    }

    free(v_schur_re); free(v_schur_im);
    free(w_re); free(w_im);
    free(inv_perm);
}

/* Emit a List of List of (Real or Complex[re, im]) for the eigenvector
 * matrix V (rows = eigenvectors).  Entries with V_im[i,j] == 0 become
 * Real; others become Complex[re, im]. */
static Expr* direct_build_complex_eigenvector_list(const double* V_re,
                                                     const double* V_im,
                                                     size_t n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            double r = V_re[k * n + i];
            double m = V_im[k * n + i];
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
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level "Direct" kernel for real non-symmetric machine-precision
 * input.  When WANT_VALUES, returns a List of eigenvalues sorted by
 * descending |lambda| (Real or Complex[re, im]).  When WANT_VECTORS,
 * returns a List of List of eigenvector components in the same sorted
 * order.  Returns NULL on convergence failure -- caller falls back to
 * the symbolic path. */
static Expr* direct_real_general_machine(const MatD* A, MateigenWant want,
                                          Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    double* H = (double*)malloc(sizeof(double) * n * n);
    memcpy(H, A->re, sizeof(double) * n * n);
    double* Q = NULL;
    if (want_Q) {
        Q = (double*)malloc(sizeof(double) * n * n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) Q[i * n + j] = (i == j) ? 1.0 : 0.0;
    }
    double* u_buf = (double*)malloc(sizeof(double) * n);
    direct_hessenberg_real(H, n, u_buf, Q);
    free(u_buf);

    double* eval_re = (double*)calloc(n, sizeof(double));
    double* eval_im = (double*)calloc(n, sizeof(double));
    int qr_status = direct_qr_real_general(H, n, eval_re, eval_im, Q);

    if (qr_status != 0) {
        free(H); free(eval_re); free(eval_im);
        if (Q) free(Q);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex(eval_re, eval_im, n, perm);

    Expr* out;
    if (want_Q) {
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        schur_compute_eigvecs(H, Q, n, eval_re, eval_im, perm, V_re, V_im);
        out = direct_build_complex_eigenvector_list(V_re, V_im, n);
        free(V_re); free(V_im);
    } else {
        out = direct_build_complex_eigenvalue_list(eval_re, eval_im, n, perm);
    }

    free(H); free(eval_re); free(eval_im); free(perm);
    if (Q) free(Q);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Top-level "Direct" kernel for real symmetric machine-precision input.
 * Returns a freshly-allocated List[Real,...] (eigenvalues) or
 * List[List[Real,...],...] (eigenvectors), or NULL when the input is
 * outside this kernel's supported domain so the dispatcher can fall
 * back to a different kernel / the symbolic path.  Caller owns the
 * result. */
static Expr* direct_real_sym_machine(const MatD* A, MateigenWant want,
                                       Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    /* Working copy of A (tridiag step modifies in place). */
    double* W = (double*)malloc(sizeof(double) * n * n);
    memcpy(W, A->re, sizeof(double) * n * n);

    double* diag    = (double*)malloc(sizeof(double) * n);
    double* sub     = (double*)calloc(n, sizeof(double));  /* length n-1, +1 slack */
    double* u       = (double*)malloc(sizeof(double) * n);
    double* p_buf   = (double*)malloc(sizeof(double) * n);
    double* q_buf   = (double*)malloc(sizeof(double) * n);
    bool want_Q     = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q       = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    direct_tridiag_real_sym(W, n, diag, sub, Q, want_Q, u, p_buf, q_buf);

    int qr_status = direct_symtridiag_qr(diag, sub, n, Q, want_Q);

    free(u); free(p_buf); free(q_buf); free(W);

    if (qr_status != 0) {
        free(diag); free(sub);
        if (Q) free(Q);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs(diag, n, perm);

    Expr* out = (want_Q)
        ? direct_build_real_eigenvector_list(Q, n, perm)
        : direct_build_real_eigenvalue_list(diag, n, perm);

    free(diag); free(sub); free(perm);
    if (Q) free(Q);

    return direct_apply_k_spec_list(out, k_spec);
}

/* ---------- Complex Hermitian Direct: Householder + phase + sym QR -- *
 *                                                                       *
 *  For a complex Hermitian n x n matrix A all eigenvalues are real and  *
 *  the eigenvectors form a complex unitary basis.  We reduce as:        *
 *                                                                       *
 *    1. Complex Householder reflectors zero the sub-column below the    *
 *       sub-diagonal at each step.  The result is a Hermitian tri-      *
 *       diagonal T with real diagonal but generally COMPLEX             *
 *       sub-diagonal.  Q (n x n complex) accumulates the reflectors.    *
 *    2. Diagonal phase correction D = diag(d_0, ..., d_{n-1}), |d_k|=1, *
 *       chosen so D^H T D has real-positive sub-diagonal.  T becomes a  *
 *       real symmetric tridiagonal.  Q is updated to Q D.               *
 *    3. The existing real symmetric tridiag QR (direct_symtridiag_qr)   *
 *       finds eigenvalues and the real orthogonal accumulator Z so that *
 *       Z^T T_real Z = Lambda.                                          *
 *    4. Final eigenvectors V = Q Z (complex Q times real Z).            *
 *                                                                       *
 *  This avoids needing a separate complex tridiagonal QR while costing  *
 *  only an O(n^2) diagonal-phase application step.  See Wilkinson 1965, *
 *  "The Algebraic Eigenvalue Problem", section 5.45.                    *
 *                                                                       *
 *  LAPACK-HOOK: this whole block maps to zhetrd (step 1) + the implicit *
 *  phase reduction (handled internally by LAPACK's zhetrd via the tau   *
 *  scalar) + dstedc/dsteqr (step 3), or simply zheevd as a one-call     *
 *  wrapper when USE_LAPACK is wired.                                    *
 * --------------------------------------------------------------------- */

/* Hermitian Householder tridiagonalisation.
 *
 * Input:  A (row-major n*n; A_re holds real parts, A_im holds imag parts).
 * Output: diag[i]     = real diagonal of T  (0 <= i < n)
 *         sub_re[k], sub_im[k] = complex sub-diagonal T_{k+1,k}, 0 <= k < n-1
 *         Q (complex, n*n via Q_re/Q_im) when want_Q: unitary accumulator
 *
 * Scratch buffers u, v, q (each length n complex via paired arrays) are
 * caller-provided so the inner loops never touch malloc/free. */
static void direct_tridiag_complex_hermitian(double* A_re, double* A_im,
                                              size_t n,
                                              double* diag,
                                              double* sub_re, double* sub_im,
                                              double* Q_re, double* Q_im,
                                              bool want_Q,
                                              double* u_re, double* u_im,
                                              double* v_re, double* v_im,
                                              double* q_re, double* q_im) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Q_re[i * n + j] = (i == j) ? 1.0 : 0.0;
                Q_im[i * n + j] = 0.0;
            }
        }
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* Norm of column-k tail x = A[k+1:n, k] (complex). */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double r = A_re[i * n + k];
            double m = A_im[i * n + k];
            sigma += r * r + m * m;
        }
        double xk1_re = A_re[(k + 1) * n + k];
        double xk1_im = A_im[(k + 1) * n + k];
        double xk1_abs = hypot(xk1_re, xk1_im);

        if (sigma == 0.0) {
            sub_re[k] = xk1_re;
            sub_im[k] = xk1_im;
            continue;
        }
        double norm_x = sqrt(sigma);

        /* alpha = -phase(x_{k+1}) * norm_x.  Phase = x/|x| (or 1 if x=0).
         * Choosing alpha "opposite" the leading entry maximises |u_0|^2 so
         * the reflector is numerically well-conditioned. */
        double alpha_re, alpha_im;
        if (xk1_abs == 0.0) {
            alpha_re = -norm_x;
            alpha_im = 0.0;
        } else {
            alpha_re = -xk1_re / xk1_abs * norm_x;
            alpha_im = -xk1_im / xk1_abs * norm_x;
        }

        /* u = x with u_{k+1} -= alpha.  Householder reflector
         * H = I - 2 u u^H / ||u||^2; we normalise u so that tau = 2. */
        u_re[k + 1] = xk1_re - alpha_re;
        u_im[k + 1] = xk1_im - alpha_im;
        for (size_t i = k + 2; i < n; i++) {
            u_re[i] = A_re[i * n + k];
            u_im[i] = A_im[i * n + k];
        }
        double unorm2 = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            unorm2 += u_re[i] * u_re[i] + u_im[i] * u_im[i];
        }
        if (unorm2 == 0.0) {
            sub_re[k] = xk1_re;
            sub_im[k] = xk1_im;
            continue;
        }
        double inv_unorm = 1.0 / sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) {
            u_re[i] *= inv_unorm;
            u_im[i] *= inv_unorm;
        }

        /* v = A u   on the trailing sub-block (A is Hermitian; we read
         * the full sub-block to be agnostic to whether the upper triangle
         * has been kept in sync). */
        for (size_t i = k + 1; i < n; i++) {
            double s_re = 0.0, s_im = 0.0;
            for (size_t j = k + 1; j < n; j++) {
                double ar = A_re[i * n + j];
                double ai = A_im[i * n + j];
                double ur = u_re[j];
                double ui = u_im[j];
                /* s += A_ij * u_j */
                s_re += ar * ur - ai * ui;
                s_im += ar * ui + ai * ur;
            }
            v_re[i] = s_re;
            v_im[i] = s_im;
        }

        /* alpha_v = u^H v = sum conj(u_i) v_i.  Real for Hermitian A. */
        double alpha_v = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            alpha_v += u_re[i] * v_re[i] + u_im[i] * v_im[i];
        }

        /* q = 2 v - 2 alpha_v u   so that  A <- A - u q^H - q u^H = H A H. */
        for (size_t i = k + 1; i < n; i++) {
            q_re[i] = 2.0 * v_re[i] - 2.0 * alpha_v * u_re[i];
            q_im[i] = 2.0 * v_im[i] - 2.0 * alpha_v * u_im[i];
        }

        /* Rank-2 update: A_ij -= u_i conj(q_j) + q_i conj(u_j). */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                double ur = u_re[i], ui = u_im[i];
                double qr = q_re[i], qi = q_im[i];
                double uj_r = u_re[j], uj_i = u_im[j];
                double qj_r = q_re[j], qj_i = q_im[j];
                /* u_i * conj(q_j) = (ur + i ui)(qj_r - i qj_i) */
                double t1_re = ur * qj_r + ui * qj_i;
                double t1_im = ui * qj_r - ur * qj_i;
                /* q_i * conj(u_j) = (qr + i qi)(uj_r - i uj_i) */
                double t2_re = qr * uj_r + qi * uj_i;
                double t2_im = qi * uj_r - qr * uj_i;
                A_re[i * n + j] -= t1_re + t2_re;
                A_im[i * n + j] -= t1_im + t2_im;
            }
        }

        /* Force the new sub-column / sub-row to their analytic values
         * to suppress drift.  Sub-column has only one non-zero entry
         * after H is applied: A[k+1, k] = alpha (complex).  All A[i, k]
         * for i > k+1 are exactly zero; conjugate row equally. */
        sub_re[k] = alpha_re;
        sub_im[k] = alpha_im;
        A_re[(k + 1) * n + k] = alpha_re;
        A_im[(k + 1) * n + k] = alpha_im;
        A_re[k * n + (k + 1)] = alpha_re;   /* T_{k,k+1} = conj(T_{k+1,k}) */
        A_im[k * n + (k + 1)] = -alpha_im;
        for (size_t i = k + 2; i < n; i++) {
            A_re[i * n + k] = 0.0;
            A_im[i * n + k] = 0.0;
            A_re[k * n + i] = 0.0;
            A_im[k * n + i] = 0.0;
        }

        /* Q <- Q H (right-multiply).  Column j -> col j - 2 (Q row * u) u_j.
         * In complex: Q_ij -= 2 (Q row . u)_i * conj(u_j) -- wait, this
         * isn't right.  H = I - 2 u u^H, so Q H = Q - 2 (Q u) u^H, which
         * in scalar form is Q_ij <- Q_ij - 2 (Q u)_i conj(u_j). */
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                /* (Q u)_i = sum_j Q_ij u_j, restricted to j >= k+1. */
                double s_re = 0.0, s_im = 0.0;
                for (size_t j = k + 1; j < n; j++) {
                    double qr2 = Q_re[i * n + j];
                    double qi2 = Q_im[i * n + j];
                    double ur2 = u_re[j];
                    double ui2 = u_im[j];
                    s_re += qr2 * ur2 - qi2 * ui2;
                    s_im += qr2 * ui2 + qi2 * ur2;
                }
                s_re *= 2.0;
                s_im *= 2.0;
                /* Q_ij -= s * conj(u_j) = s * (u_re[j] - i u_im[j]) */
                for (size_t j = k + 1; j < n; j++) {
                    double ur2 = u_re[j];
                    double ui2 = u_im[j];
                    /* (s_re + i s_im)(ur2 - i ui2) */
                    double pr = s_re * ur2 + s_im * ui2;
                    double pi = s_im * ur2 - s_re * ui2;
                    Q_re[i * n + j] -= pr;
                    Q_im[i * n + j] -= pi;
                }
            }
        }
    }

    /* Extract diagonal (real for Hermitian, but defensively take the
     * real part). */
    for (size_t i = 0; i < n; i++) diag[i] = A_re[i * n + i];
    /* The (n-2 -> n-1) sub-diagonal entry hasn't been overwritten yet. */
    if (n >= 2) {
        sub_re[n - 2] = A_re[(n - 1) * n + (n - 2)];
        sub_im[n - 2] = A_im[(n - 1) * n + (n - 2)];
    }
}

/* Diagonal phase correction.  Multiplies the complex Hermitian tri-
 * diagonal T by a diagonal unitary D = diag(d_0, ..., d_{n-1}), |d_k|=1,
 * chosen so D^H T D has real-positive sub-diagonal.  Updates Q (the
 * Householder accumulator) by post-multiplication: Q <- Q D. */
static void direct_phase_correct_tridiag(double* sub_re, double* sub_im,
                                           size_t n,
                                           double* Q_re, double* Q_im,
                                           bool want_Q) {
    /* d_0 = 1; d_{k+1} = d_k * sub[k] / |sub[k]|. */
    double d_re = 1.0, d_im = 0.0;
    /* k = 0 -> column 1 of Q gets multiplied by d_1, column 0 stays. */
    for (size_t k = 0; k + 1 < n; k++) {
        double sr = sub_re[k], si = sub_im[k];
        double mag = hypot(sr, si);
        double phase_re, phase_im;
        if (mag == 0.0) { phase_re = 1.0; phase_im = 0.0; }
        else            { phase_re = sr / mag; phase_im = si / mag; }
        /* d_{k+1} = d_k * phase */
        double new_d_re = d_re * phase_re - d_im * phase_im;
        double new_d_im = d_re * phase_im + d_im * phase_re;
        d_re = new_d_re;
        d_im = new_d_im;
        sub_re[k] = mag;
        sub_im[k] = 0.0;
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                double qr = Q_re[i * n + (k + 1)];
                double qi = Q_im[i * n + (k + 1)];
                Q_re[i * n + (k + 1)] = qr * d_re - qi * d_im;
                Q_im[i * n + (k + 1)] = qr * d_im + qi * d_re;
            }
        }
    }
}

/* Compose complex Q (n*n) with real Z (n*n) into complex V (n*n) via
 * V = Q Z.  V_re/V_im are caller-allocated. */
static void direct_compose_complex_Q_real_Z(const double* Q_re,
                                              const double* Q_im,
                                              const double* Z, size_t n,
                                              double* V_re, double* V_im) {
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            double sr = 0.0, si = 0.0;
            for (size_t j = 0; j < n; j++) {
                double z = Z[j * n + k];
                sr += Q_re[i * n + j] * z;
                si += Q_im[i * n + j] * z;
            }
            V_re[i * n + k] = sr;
            V_im[i * n + k] = si;
        }
    }
}

/* Build a List of List of (Real or Complex[re, im]) from a complex
 * eigenvector matrix V (n*n; V[i,k] is the i-th component of the k-th
 * eigenvector) in the order given by `perm` (perm[r] -> column of V). */
static Expr* direct_build_complex_hermitian_eigvec_list(const double* V_re,
                                                          const double* V_im,
                                                          size_t n,
                                                          const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t r = 0; r < n; r++) {
        size_t col = perm[r];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            double rv = V_re[i * n + col];
            double iv = V_im[i * n + col];
            if (iv == 0.0) {
                comps[i] = expr_new_real(rv);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_real(rv);
                args[1] = expr_new_real(iv);
                comps[i] = expr_new_function(expr_new_symbol("Complex"), args, 2);
                free(args);
            }
        }
        rows[r] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level "Direct" kernel for complex Hermitian machine-precision
 * input.  Eigenvalues are real (sorted by descending |lambda|);
 * eigenvectors are complex unitary.  Returns NULL on convergence
 * failure so the caller falls back to the symbolic path. */
static Expr* direct_complex_hermitian_machine(const MatD* A, MateigenWant want,
                                                Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    /* Working complex copy of A (tridiag step modifies in place). */
    double* W_re = (double*)malloc(sizeof(double) * n * n);
    double* W_im = (double*)malloc(sizeof(double) * n * n);
    memcpy(W_re, A->re, sizeof(double) * n * n);
    memcpy(W_im, A->im, sizeof(double) * n * n);

    double* diag   = (double*)malloc(sizeof(double) * n);
    double* sub_re = (double*)calloc(n, sizeof(double));   /* len n-1 + slack */
    double* sub_im = (double*)calloc(n, sizeof(double));

    double* u_re = (double*)malloc(sizeof(double) * n);
    double* u_im = (double*)malloc(sizeof(double) * n);
    double* v_re = (double*)malloc(sizeof(double) * n);
    double* v_im = (double*)malloc(sizeof(double) * n);
    double* q_re = (double*)malloc(sizeof(double) * n);
    double* q_im = (double*)malloc(sizeof(double) * n);

    bool want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q_re  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;
    double* Q_im  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    direct_tridiag_complex_hermitian(W_re, W_im, n,
                                      diag, sub_re, sub_im,
                                      Q_re, Q_im, want_Q,
                                      u_re, u_im, v_re, v_im, q_re, q_im);

    free(W_re); free(W_im);
    free(u_re); free(u_im); free(v_re); free(v_im); free(q_re); free(q_im);

    /* Phase-correct so sub-diagonal becomes real positive. */
    direct_phase_correct_tridiag(sub_re, sub_im, n, Q_re, Q_im, want_Q);

    /* Real symmetric tridiagonal QR.  We accumulate the real orthogonal
     * Z separately when eigenvectors are wanted -- it composes with the
     * complex Q from the Hermitian Householder step. */
    double* Z = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) Z[i * n + j] = (i == j) ? 1.0 : 0.0;
    }
    int qr_status = direct_symtridiag_qr(diag, sub_re, n, Z, want_Q);

    free(sub_re); free(sub_im);

    if (qr_status != 0) {
        free(diag);
        if (Q_re) { free(Q_re); free(Q_im); }
        if (Z) free(Z);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs(diag, n, perm);

    Expr* out;
    if (want_Q) {
        /* V = Q Z (complex). */
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        direct_compose_complex_Q_real_Z(Q_re, Q_im, Z, n, V_re, V_im);
        out = direct_build_complex_hermitian_eigvec_list(V_re, V_im, n, perm);
        free(V_re); free(V_im);
    } else {
        out = direct_build_real_eigenvalue_list(diag, n, perm);
    }

    free(diag); free(perm);
    if (Q_re) { free(Q_re); free(Q_im); }
    if (Z) free(Z);

    return direct_apply_k_spec_list(out, k_spec);
}

/* ---------- Complex general Direct: real block-embedding ------------- *
 *                                                                       *
 *  A complex n x n matrix A = R + i S is reduced to a real 2n x 2n      *
 *  block matrix                                                         *
 *                                                                       *
 *      M = [[ R, -S ],                                                  *
 *           [ S,  R ]]                                                  *
 *                                                                       *
 *  whose spectrum is spec(A) ∪ conj(spec(A)) as a multiset, with each   *
 *  pair (λ, conj λ) sharing a 2-dim real invariant subspace of M.       *
 *                                                                       *
 *  We dispatch M through the existing real Hessenberg + Francis QR      *
 *  pipeline (direct_hessenberg_real + direct_qr_real_general +          *
 *  schur_compute_eigvecs) and then unpair:                              *
 *                                                                       *
 *    - For each complex M-eigenvalue with positive imaginary part:      *
 *      find its conjugate partner among the M-eigenvalues, mark both    *
 *      as used, and emit the positive-imag one as A's eigenvalue.       *
 *    - For each real M-eigenvalue: pair it with the closest unused      *
 *      real M-eigenvalue, mark both as used, emit one copy.             *
 *                                                                       *
 *  The complex A-eigenvector is recovered from M's complex eigenvector  *
 *  w = w_re + i w_im (length 2n) via                                    *
 *                                                                       *
 *      x = (a - d) + i (b + c)                                          *
 *                                                                       *
 *  where a, b, c, d are the n-vector top/bottom splits of w_re / w_im   *
 *  respectively.  For real M-eigenvalues w_im = 0 and the formula       *
 *  collapses to x = top(w_re) + i * bot(w_re).  See test                *
 *  test_direct_general_complex_machine_2x2_diagonal_imag for a worked   *
 *  example.                                                             *
 *                                                                       *
 *  Cost: O((2n)^3) = 8 O(n^3) -- nominally 8x more flops than a native  *
 *  complex Hessenberg + complex QR at the same n.  In return the entire *
 *  numerical machinery is reused (real LAPACK-mappable kernels) and     *
 *  there are no separate complex tridiagonal / Givens implementations   *
 *  to maintain.  When USE_LAPACK is wired, the obvious replacement is   *
 *  a native zgehrd + zhseqr + ztrevc3 path that won't pay the 8x.       *
 *                                                                       *
 *  LAPACK-HOOK: zgehrd + zhseqr + ztrevc3 (or zgeev as a one-call       *
 *  wrapper) can replace this whole block when USE_LAPACK is set.        *
 * --------------------------------------------------------------------- */

static Expr* direct_complex_general_machine(const MatD* A, MateigenWant want,
                                              Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;
    size_t N = 2 * n;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    /* Build real 2n x 2n block embedding. */
    double* H = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double r = A->re[i * n + j];
            double s = A->im[i * n + j];
            H[i * N + j]                   =  r;
            H[i * N + (j + n)]             = -s;
            H[(i + n) * N + j]             =  s;
            H[(i + n) * N + (j + n)]       =  r;
        }
    }

    /* Hessenberg reduction.  We ALWAYS accumulate Q -- M's eigenvectors
     * are needed for the values-only path too, to disambiguate which of
     * the conjugate pair (mu, conj mu) is in spec(A) vs spec(conj A). */
    double* Q = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) Q[i * N + j] = (i == j) ? 1.0 : 0.0;
    double* u_buf = (double*)malloc(sizeof(double) * N);
    direct_hessenberg_real(H, N, u_buf, Q);
    free(u_buf);

    /* Francis QR -> Schur form + complex-conjugate eigenvalue pairs. */
    double* M_eval_re = (double*)calloc(N, sizeof(double));
    double* M_eval_im = (double*)calloc(N, sizeof(double));
    int qr_status = direct_qr_real_general(H, N, M_eval_re, M_eval_im, Q);
    if (qr_status != 0) {
        free(H); free(Q); free(M_eval_re); free(M_eval_im);
        return NULL;
    }

    /* Eigenvectors of M (in the original basis), one row per eigenvalue. */
    double* M_evec_re = (double*)malloc(sizeof(double) * N * N);
    double* M_evec_im = (double*)malloc(sizeof(double) * N * N);
    size_t* identity_perm = (size_t*)malloc(sizeof(size_t) * N);
    for (size_t i = 0; i < N; i++) identity_perm[i] = i;
    schur_compute_eigvecs(H, Q, N, M_eval_re, M_eval_im, identity_perm,
                            M_evec_re, M_evec_im);
    free(identity_perm);
    free(H); free(Q);

    /* Recover A's eigenvalues from M's spectrum:
     *
     * For each M-eigenvector w of eigenvalue mu, the vector
     *
     *     x = (a - d) + i (b + c)
     *
     * with w_re = [a; b] and w_im = [c; d] satisfies A x = mu x.  When mu
     * is in spec(A) the candidate x is non-zero; when mu is only in
     * spec(conj A) (i.e., conj mu in spec(A) but mu itself is not)
     * the candidate x is zero.  This disambiguates the conjugate-pair
     * doubling of spec(M) = spec(A) U conj(spec(A)).
     *
     * Grouped Gram-Schmidt handles algebraic multiplicity:
     *   1. Walk M-eigenvalues; group adjacent (after stable Schur order)
     *      values that fall within `group_tol`.
     *   2. Within each group, project each candidate x against the
     *      already-emitted vectors of the group and emit it iff the
     *      remaining norm stays above `extract_threshold`.
     *
     * This produces m_A(mu) ortho-normal eigenvectors per distinct mu --
     * the rank of the +J subspace in M's mu-eigenspace -- which sums to
     * exactly n. */
    double spec_norm = 0.0;
    for (size_t i = 0; i < N; i++) {
        double a = hypot(M_eval_re[i], M_eval_im[i]);
        if (a > spec_norm) spec_norm = a;
    }
    double group_tol = 1e-8 * (spec_norm == 0.0 ? 1.0 : spec_norm) * (double)N;
    /* Norm threshold: a "valid" candidate (in +J subspace) has |x| in
     * the range [eps, sqrt(2)]; a "wrong" candidate (in -J subspace) has
     * |x| ~ machine_eps.  1e-8 is safely between these. */
    double extract_threshold = sqrt((double)n) * 1e-9;

    int* used = (int*)calloc(N, sizeof(int));
    double* A_eval_re = (double*)malloc(sizeof(double) * n);
    double* A_eval_im = (double*)malloc(sizeof(double) * n);
    double* A_evec_re = (double*)calloc(n * n, sizeof(double));
    double* A_evec_im = (double*)calloc(n * n, sizeof(double));
    double* cand_re = (double*)malloc(sizeof(double) * n);
    double* cand_im = (double*)malloc(sizeof(double) * n);
    size_t out = 0;

    for (size_t i = 0; i < N && out < n; i++) {
        if (used[i]) continue;
        used[i] = 1;
        size_t group_start = out;

        /* Process group member j; iteration includes i itself. */
        for (size_t j = i; j < N && out < n; j++) {
            if (j != i) {
                if (used[j]) continue;
                double dr = M_eval_re[j] - M_eval_re[i];
                double di = M_eval_im[j] - M_eval_im[i];
                if (hypot(dr, di) > group_tol) continue;
                used[j] = 1;
            }
            /* Candidate x = (a - d) + i (b + c). */
            for (size_t l = 0; l < n; l++) {
                double a = M_evec_re[j * N + l];
                double b = M_evec_re[j * N + (l + n)];
                double c = M_evec_im[j * N + l];
                double d_im = M_evec_im[j * N + (l + n)];
                cand_re[l] = a - d_im;
                cand_im[l] = b + c;
            }
            /* Orthogonalise against already-emitted vectors in this
             * group (complex modified Gram-Schmidt).  Twice for numerical
             * stability per the "twice is enough" rule. */
            for (int pass = 0; pass < 2; pass++) {
                for (size_t f = group_start; f < out; f++) {
                    double pr = 0.0, pi = 0.0;
                    for (size_t l = 0; l < n; l++) {
                        double vr = A_evec_re[f * n + l];
                        double vi = A_evec_im[f * n + l];
                        /* conj(V_f) . cand = (vr - i vi)(cand_re + i cand_im) */
                        pr += vr * cand_re[l] + vi * cand_im[l];
                        pi += vr * cand_im[l] - vi * cand_re[l];
                    }
                    for (size_t l = 0; l < n; l++) {
                        double vr = A_evec_re[f * n + l];
                        double vi = A_evec_im[f * n + l];
                        /* (pr + i pi)(vr + i vi) */
                        double pvr = pr * vr - pi * vi;
                        double pvi = pr * vi + pi * vr;
                        cand_re[l] -= pvr;
                        cand_im[l] -= pvi;
                    }
                }
            }
            double norm2 = 0.0;
            for (size_t l = 0; l < n; l++) {
                norm2 += cand_re[l] * cand_re[l] + cand_im[l] * cand_im[l];
            }
            if (norm2 < extract_threshold * extract_threshold) continue;
            double inv = 1.0 / sqrt(norm2);
            for (size_t l = 0; l < n; l++) {
                A_evec_re[out * n + l] = cand_re[l] * inv;
                A_evec_im[out * n + l] = cand_im[l] * inv;
            }
            A_eval_re[out] = M_eval_re[i];
            A_eval_im[out] = M_eval_im[i];
            out++;
        }
    }
    free(used); free(M_eval_re); free(M_eval_im);
    free(M_evec_re); free(M_evec_im);
    free(cand_re); free(cand_im);

    if (out != n) {
        /* Extraction under-produced -- bail out so the symbolic path
         * can take over.  Most commonly hit when extract_threshold is
         * too tight or the matrix is wildly ill-conditioned. */
        free(A_eval_re); free(A_eval_im);
        free(A_evec_re); free(A_evec_im);
        return NULL;
    }

    /* Sort by descending |lambda| and emit. */
    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex(A_eval_re, A_eval_im, n, perm);

    Expr* result;
    if (want_Q) {
        /* Permute eigenvector rows into sort order. */
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        for (size_t k = 0; k < n; k++) {
            size_t src = perm[k];
            for (size_t l = 0; l < n; l++) {
                V_re[k * n + l] = A_evec_re[src * n + l];
                V_im[k * n + l] = A_evec_im[src * n + l];
            }
        }
        result = direct_build_complex_eigenvector_list(V_re, V_im, n);
        free(V_re); free(V_im);
    } else {
        result = direct_build_complex_eigenvalue_list(A_eval_re, A_eval_im,
                                                        n, perm);
    }

    free(A_eval_re); free(A_eval_im); free(perm);
    free(A_evec_re); free(A_evec_im);

    return direct_apply_k_spec_list(result, k_spec);
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
typedef struct {
    size_t      n;
    int         is_complex;     /* 0: real only, 1: complex (im[] non-NULL) */
    mpfr_prec_t bits;
    mpfr_t*     re;             /* length n*n, mpfr_init2'd to bits */
    mpfr_t*     im;             /* length n*n, mpfr_init2'd, NULL when !is_complex */
} MatM;

/* Free a heap-allocated mpfr_t[] of length `count`, clearing each cell. */
static void mpfr_array_free(mpfr_t* a, size_t count) {
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* Allocate and pre-initialise an mpfr_t[] of length `count` to `bits`. */
static mpfr_t* mpfr_array_alloc(size_t count, mpfr_prec_t bits) {
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void matM_free(MatM* M) {
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
static bool matM_load(Expr* m, size_t n, mpfr_prec_t bits, MatM* out) {
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
static void matM_norm_inf_real(const mpfr_t* A, size_t n,
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
static bool matM_is_real_symmetric(const MatM* A, const mpfr_t tol) {
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

/* MPFR variant of `direct_tridiag_real_sym` -- see that routine for
 * the algorithm.  Workspace mpfr_t* arrays are caller-supplied and
 * pre-initialised; `tmp` holds 8 scratch cells used in the inner
 * loops. */
static void direct_tridiag_real_sym_M(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                       mpfr_t* diag, mpfr_t* sub,
                                       mpfr_t* Q, bool want_Q,
                                       mpfr_t* u, mpfr_t* p, mpfr_t* q,
                                       mpfr_t* tmp /* >= 8 cells */) {
    (void)bits;
    mpfr_t* sigma   = &tmp[0];
    mpfr_t* xk1     = &tmp[1];
    mpfr_t* norm_x  = &tmp[2];
    mpfr_t* alpha   = &tmp[3];
    mpfr_t* unorm2  = &tmp[4];
    mpfr_t* unorm   = &tmp[5];
    mpfr_t* K       = &tmp[6];
    mpfr_t* s_acc   = &tmp[7];

    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_si(Q[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* sigma = sum_{i=k+1..n-1} A[i,k]^2 */
        mpfr_set_zero(*sigma, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, A[i * n + k], A[i * n + k], MPFR_RNDN);
            mpfr_add(*sigma, *sigma, *s_acc, MPFR_RNDN);
        }
        mpfr_set(*xk1, A[(k + 1) * n + k], MPFR_RNDN);
        if (mpfr_zero_p(*sigma)) {
            mpfr_set(sub[k], *xk1, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(*norm_x, *sigma, MPFR_RNDN);
        /* alpha = -sign(xk1) * ||x||  (cancellation-safe) */
        if (mpfr_sgn(*xk1) >= 0) mpfr_neg(*alpha, *norm_x, MPFR_RNDN);
        else                     mpfr_set(*alpha, *norm_x, MPFR_RNDN);

        /* u[k+1] = xk1 - alpha;  u[i>k+1] = A[i,k] */
        mpfr_sub(u[k + 1], *xk1, *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) mpfr_set(u[i], A[i * n + k], MPFR_RNDN);

        /* unorm2 = sum u^2 */
        mpfr_mul(*unorm2, u[k + 1], u[k + 1], MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_mul(*s_acc, u[i], u[i], MPFR_RNDN);
            mpfr_add(*unorm2, *unorm2, *s_acc, MPFR_RNDN);
        }
        if (mpfr_zero_p(*unorm2)) {
            mpfr_set(sub[k], *xk1, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(*unorm, *unorm2, MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++)
            mpfr_div(u[i], u[i], *unorm, MPFR_RNDN);

        /* p[i] = 2 * sum_j A[i,j] u[j]   on the trailing block */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_set_zero(*s_acc, 1);
            for (size_t j = k + 1; j < n; j++) {
                mpfr_t prod;
                mpfr_init2(prod, mpfr_get_prec(*s_acc));
                mpfr_mul(prod, A[i * n + j], u[j], MPFR_RNDN);
                mpfr_add(*s_acc, *s_acc, prod, MPFR_RNDN);
                mpfr_clear(prod);
            }
            mpfr_mul_2si(p[i], *s_acc, 1, MPFR_RNDN);   /* p[i] = 2 * s_acc */
        }
        /* K = (u^T p) / 2 */
        mpfr_set_zero(*K, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, u[i], p[i], MPFR_RNDN);
            mpfr_add(*K, *K, *s_acc, MPFR_RNDN);
        }
        mpfr_div_2si(*K, *K, 1, MPFR_RNDN);
        /* q[i] = p[i] - 2 K u[i] */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, *K, u[i], MPFR_RNDN);
            mpfr_mul_2si(*s_acc, *s_acc, 1, MPFR_RNDN);
            mpfr_sub(q[i], p[i], *s_acc, MPFR_RNDN);
        }
        /* A_22 -= u q^T + q u^T */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                mpfr_mul(*s_acc, u[i], q[j], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *s_acc, MPFR_RNDN);
                mpfr_mul(*s_acc, q[i], u[j], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *s_acc, MPFR_RNDN);
            }
        }
        /* Set subdiag + clear eliminated entries (drift control). */
        mpfr_set(sub[k], *alpha, MPFR_RNDN);
        mpfr_set(A[(k + 1) * n + k], *alpha, MPFR_RNDN);
        mpfr_set(A[k * n + (k + 1)], *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_set_zero(A[i * n + k], 1);
            mpfr_set_zero(A[k * n + i], 1);
        }

        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(*s_acc, 1);
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_t prod;
                    mpfr_init2(prod, mpfr_get_prec(*s_acc));
                    mpfr_mul(prod, Q[i * n + j], u[j], MPFR_RNDN);
                    mpfr_add(*s_acc, *s_acc, prod, MPFR_RNDN);
                    mpfr_clear(prod);
                }
                mpfr_mul_2si(*s_acc, *s_acc, 1, MPFR_RNDN);  /* 2 * (Q_row . u) */
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_t prod;
                    mpfr_init2(prod, mpfr_get_prec(*s_acc));
                    mpfr_mul(prod, *s_acc, u[j], MPFR_RNDN);
                    mpfr_sub(Q[i * n + j], Q[i * n + j], prod, MPFR_RNDN);
                    mpfr_clear(prod);
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) mpfr_set(diag[i], A[i * n + i], MPFR_RNDN);
    if (n >= 2) mpfr_set(sub[n - 2], A[(n - 1) * n + (n - 2)], MPFR_RNDN);
}

/* MPFR variant of `direct_symtridiag_qr`.  `tmp` provides 12+ scratch
 * cells.  rel_tol scales as 2^{-bits+3}: a few ULPs at the working
 * precision, mirroring the 1e-14 ~= 2^-46 (for bits=53) used in the
 * machine version. */
static int direct_symtridiag_qr_M(mpfr_t* diag, mpfr_t* sub, size_t n,
                                    mpfr_prec_t bits,
                                    mpfr_t* Q, bool want_Q,
                                    mpfr_t* tmp /* >= 12 cells */) {
    if (n == 0) return 0;
    const size_t max_sweeps = 30 * n;
    size_t sweeps = 0;

    mpfr_t* tol     = &tmp[0];
    mpfr_t* d       = &tmp[1];
    mpfr_t* e       = &tmp[2];
    mpfr_t* t       = &tmp[3];
    mpfr_t* mu      = &tmp[4];
    mpfr_t* x       = &tmp[5];
    mpfr_t* z       = &tmp[6];
    mpfr_t* c       = &tmp[7];
    mpfr_t* s       = &tmp[8];
    mpfr_t* r       = &tmp[9];
    mpfr_t* scratch1= &tmp[10];
    mpfr_t* scratch2= &tmp[11];

    /* Relative tolerance: ~ 8 * 2^-bits == a few ULPs. */
    mpfr_t rel_tol;
    mpfr_init2(rel_tol, bits);
    mpfr_set_ui(rel_tol, 1, MPFR_RNDN);
    mpfr_div_2si(rel_tol, rel_tol, (long)bits - 3, MPFR_RNDN);

    size_t end = n;
    while (end > 1) {
        /* Find largest m s.t. sub[m..end-2] are all significant. */
        size_t m = end - 1;
        while (m > 0) {
            mpfr_abs(*scratch1, diag[m - 1], MPFR_RNDN);
            mpfr_abs(*scratch2, diag[m], MPFR_RNDN);
            mpfr_add(*scratch1, *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*tol, rel_tol, *scratch1, MPFR_RNDN);
            mpfr_abs(*scratch1, sub[m - 1], MPFR_RNDN);
            if (mpfr_cmp(*scratch1, *tol) <= 0) {
                mpfr_set_zero(sub[m - 1], 1);
                break;
            }
            m--;
        }
        if (m == end - 1) { end--; continue; }

        if (++sweeps > max_sweeps) {
            mpfr_clear(rel_tol);
            return -1;
        }

        /* Wilkinson shift on trailing 2x2 block. */
        mpfr_sub(*d, diag[end - 2], diag[end - 1], MPFR_RNDN);
        mpfr_div_2si(*d, *d, 1, MPFR_RNDN);                   /* d = (d_{e-2} - d_{e-1}) / 2 */
        mpfr_set(*e, sub[end - 2], MPFR_RNDN);
        if (mpfr_zero_p(*d)) {
            mpfr_abs(*t, *e, MPFR_RNDN);
        } else {
            mpfr_hypot(*scratch1, *d, *e, MPFR_RNDN);
            mpfr_abs(*scratch2, *d, MPFR_RNDN);
            mpfr_add(*t, *scratch2, *scratch1, MPFR_RNDN);
        }
        /* mu = d_{e-1} - sign(d) * e^2 / t */
        mpfr_mul(*scratch1, *e, *e, MPFR_RNDN);
        if (mpfr_zero_p(*t)) {
            mpfr_set_zero(*scratch1, 1);
        } else {
            mpfr_div(*scratch1, *scratch1, *t, MPFR_RNDN);
        }
        if (mpfr_sgn(*d) < 0) mpfr_neg(*scratch1, *scratch1, MPFR_RNDN);
        mpfr_sub(*mu, diag[end - 1], *scratch1, MPFR_RNDN);

        /* Implicit QR sweep using Givens rotations on [m..end-1]. */
        mpfr_sub(*x, diag[m], *mu, MPFR_RNDN);
        mpfr_set(*z, sub[m], MPFR_RNDN);
        for (size_t k = m; k < end - 1; k++) {
            mpfr_hypot(*r, *x, *z, MPFR_RNDN);
            if (mpfr_zero_p(*r)) {
                mpfr_set_ui(*c, 1, MPFR_RNDN);
                mpfr_set_zero(*s, 1);
            } else {
                mpfr_div(*c, *x, *r, MPFR_RNDN);
                mpfr_div(*s, *z, *r, MPFR_RNDN);
            }

            if (k > m) mpfr_set(sub[k - 1], *r, MPFR_RNDN);

            /* Snapshot d_k, d_k+1, e_k. */
            mpfr_t d_k, d_k1, e_k;
            mpfr_init2(d_k,  bits); mpfr_set(d_k,  diag[k],     MPFR_RNDN);
            mpfr_init2(d_k1, bits); mpfr_set(d_k1, diag[k + 1], MPFR_RNDN);
            mpfr_init2(e_k,  bits); mpfr_set(e_k,  sub[k],      MPFR_RNDN);

            /* diag[k]   = c^2 d_k + 2 c s e_k + s^2 d_k1 */
            mpfr_mul(*scratch1, *c, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_mul_2si(*scratch2, *scratch2, 1, MPFR_RNDN);
            mpfr_add(diag[k], *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*scratch1, *s, *s, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k1, MPFR_RNDN);
            mpfr_add(diag[k], diag[k], *scratch1, MPFR_RNDN);

            /* diag[k+1] = s^2 d_k - 2 c s e_k + c^2 d_k1 */
            mpfr_mul(*scratch1, *s, *s, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_mul_2si(*scratch2, *scratch2, 1, MPFR_RNDN);
            mpfr_sub(diag[k + 1], *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*scratch1, *c, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k1, MPFR_RNDN);
            mpfr_add(diag[k + 1], diag[k + 1], *scratch1, MPFR_RNDN);

            /* sub[k] = c s (d_k1 - d_k) + (c^2 - s^2) e_k */
            mpfr_sub(*scratch1, d_k1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *c, MPFR_RNDN);
            mpfr_mul(d_k, *s, *s, MPFR_RNDN);                 /* reusing d_k as throwaway */
            mpfr_sub(*scratch2, *scratch2, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_add(sub[k], *scratch1, *scratch2, MPFR_RNDN);

            mpfr_clear(d_k); mpfr_clear(d_k1); mpfr_clear(e_k);

            /* Bulge chase: next x, z. */
            if (k + 1 < end - 1) {
                mpfr_t t_next;
                mpfr_init2(t_next, bits);
                mpfr_set(t_next, sub[k + 1], MPFR_RNDN);
                mpfr_set(*x, sub[k], MPFR_RNDN);
                mpfr_mul(*z, *s, t_next, MPFR_RNDN);
                mpfr_mul(sub[k + 1], *c, t_next, MPFR_RNDN);
                mpfr_clear(t_next);
            }

            /* Q <- Q * Givens (post-multiply cols k, k+1). */
            if (want_Q) {
                for (size_t i = 0; i < n; i++) {
                    mpfr_t qk, qk1;
                    mpfr_init2(qk,  bits); mpfr_set(qk,  Q[i * n + k],     MPFR_RNDN);
                    mpfr_init2(qk1, bits); mpfr_set(qk1, Q[i * n + (k + 1)], MPFR_RNDN);
                    mpfr_mul(*scratch1, *c, qk,  MPFR_RNDN);
                    mpfr_mul(*scratch2, *s, qk1, MPFR_RNDN);
                    mpfr_add(Q[i * n + k], *scratch1, *scratch2, MPFR_RNDN);
                    mpfr_mul(*scratch1, *s, qk,  MPFR_RNDN);
                    mpfr_mul(*scratch2, *c, qk1, MPFR_RNDN);
                    mpfr_sub(Q[i * n + (k + 1)], *scratch2, *scratch1, MPFR_RNDN);
                    mpfr_clear(qk); mpfr_clear(qk1);
                }
            }
        }
    }
    mpfr_clear(rel_tol);
    return 0;
}

/* Sort permutation by descending |vals[i]|, stable. */
static void direct_sort_perm_desc_abs_M(const mpfr_t* vals, size_t n,
                                         size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        mpfr_t ac;
        mpfr_init2(ac, mpfr_get_prec(vals[0]));
        mpfr_abs(ac, vals[cur], MPFR_RNDN);
        size_t j = i;
        while (j > 0) {
            mpfr_t ap;
            mpfr_init2(ap, mpfr_get_prec(vals[0]));
            mpfr_abs(ap, vals[perm[j - 1]], MPFR_RNDN);
            int cmp = mpfr_cmp(ap, ac);
            mpfr_clear(ap);
            if (cmp > 0 || (cmp == 0 && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
        mpfr_clear(ac);
    }
}

/* Build a List[MPFR, ...] of eigenvalues from a sorted permutation. */
static Expr* direct_build_real_eigenvalue_list_M(const mpfr_t* vals, size_t n,
                                                  const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        items[i] = expr_new_mpfr_copy(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* Build the eigenvector list (rows of MPFR n-vectors). */
static Expr* direct_build_real_eigenvector_list_M(const mpfr_t* Q, size_t n,
                                                   const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            comps[i] = expr_new_mpfr_copy(Q[i * n + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Real symmetric MPFR entry point. */
static Expr* direct_real_sym_mpfr(const MatM* A, MateigenWant want,
                                    Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    if (n == 0) {
        Expr* empty = expr_new_function(expr_new_symbol("List"), NULL, 0);
        return direct_apply_k_spec_list(empty, k_spec);
    }

    /* Workspace: copy of A (gets destroyed), tridiag arrays, optional Q,
     * Householder vectors u/p/q, and scratch tmp[12]. */
    mpfr_t* Awork = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) mpfr_set(Awork[i], A->re[i], MPFR_RNDN);

    mpfr_t* diag = mpfr_array_alloc(n, bits);
    mpfr_t* sub  = mpfr_array_alloc(n > 0 ? n - 1 + 1 : 0, bits);
    /* (allocate n cells so we can address sub[n-2] cleanly even for n==1) */
    mpfr_t* Q    = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* u    = mpfr_array_alloc(n, bits);
    mpfr_t* p    = mpfr_array_alloc(n, bits);
    mpfr_t* q    = mpfr_array_alloc(n, bits);
    mpfr_t* tmp  = mpfr_array_alloc(12, bits);

    /* Special case n == 1: trivial.  Skip tridiag. */
    if (n == 1) {
        mpfr_set(diag[0], Awork[0], MPFR_RNDN);
        if (want_Q) mpfr_set_ui(Q[0], 1, MPFR_RNDN);
    } else {
        direct_tridiag_real_sym_M(Awork, n, bits, diag, sub, Q, want_Q,
                                   u, p, q, tmp);
        direct_symtridiag_qr_M(diag, sub, n, bits, Q, want_Q, tmp);
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_M(diag, n, perm);

    Expr* result;
    if (want_Q) {
        result = direct_build_real_eigenvector_list_M(Q, n, perm);
    } else {
        result = direct_build_real_eigenvalue_list_M(diag, n, perm);
    }
    free(perm);

    mpfr_array_free(Awork, n * n);
    mpfr_array_free(diag, n);
    mpfr_array_free(sub, n > 0 ? n - 1 + 1 : 0);
    if (Q) mpfr_array_free(Q, n * n);
    mpfr_array_free(u, n);
    mpfr_array_free(p, n);
    mpfr_array_free(q, n);
    mpfr_array_free(tmp, 12);

    return direct_apply_k_spec_list(result, k_spec);
}

/* Dispatcher: route an MPFR-precision matrix through the appropriate
 * Direct kernel.  Returns NULL when the matrix shape isn't yet wired
 * (in 2d-A: only real symmetric); the outer dispatcher falls through
 * to the machine kernel (downgraded precision) and ultimately the
 * symbolic path. */
static Expr* direct_dispatch_mpfr(Expr* m, Expr* a, int64_t n,
                                    mpfr_prec_t bits,
                                    MateigenWant want, Expr* k_spec) {
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    MatM A;
    if (!matM_load(m, (size_t)n, bits, &A)) return NULL;

    Expr* out = NULL;
    if (!A.is_complex) {
        /* Symmetry tolerance: n * 2^{-bits+4} * ||A||_inf. */
        mpfr_t norm, tol, factor;
        mpfr_init2(norm,   bits);
        mpfr_init2(tol,    bits);
        mpfr_init2(factor, bits);
        matM_norm_inf_real(A.re, A.n, bits, norm);
        if (mpfr_zero_p(norm)) mpfr_set_ui(norm, 1, MPFR_RNDN);
        mpfr_set_ui(factor, 1, MPFR_RNDN);
        mpfr_div_2si(factor, factor, (long)bits - 4, MPFR_RNDN);
        mpfr_mul_ui(factor, factor, (unsigned long)A.n, MPFR_RNDN);
        mpfr_mul(tol, norm, factor, MPFR_RNDN);
        bool sym = matM_is_real_symmetric(&A, tol);
        mpfr_clear(norm); mpfr_clear(tol); mpfr_clear(factor);

        if (sym) out = direct_real_sym_mpfr(&A, want, k_spec);
    }
    /* Complex / real-general MPFR paths arrive in 2d-{B,C,D}. */

    matM_free(&A);
    return out;
}

#endif /* USE_MPFR */

/* Dispatcher entry point: route a numeric matrix through the
 * appropriate "Direct" kernel.  Returns NULL when the matrix shape
 * isn't yet supported by a numerical kernel so the caller can fall
 * back to the symbolic path.  This NULL return is also used for
 * Eigenvalues / Eigenvectors combined with a generalised pencil
 * ({m, a}) -- generalised numeric eigenvalues are not part of the
 * current numerical scope.
 *
 * Implemented kernels:
 *   - Real symmetric (machine precision):           values + vectors.
 *   - Real non-symmetric (machine precision):       values + vectors.
 *   - Complex Hermitian (machine precision):        values + vectors.
 *   - Complex non-Hermitian (machine precision):    values + vectors.
 *   - Real symmetric MPFR (step 2d-A):              values + vectors.
 *
 * Falls back to symbolic for:
 *   - MPFR real non-symmetric (step 2d-B).
 *   - MPFR complex (step 2d-{C,D}).
 */
static Expr* direct_dispatch_machine(Expr* m, Expr* a, int64_t n,
                                       MateigenWant want, Expr* k_spec) {
    if (a != NULL) return NULL;          /* generalised: symbolic only */
    if (n <= 0)    return NULL;

    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return NULL;

    Expr* out = NULL;
    if (A.is_complex) {
        double norm = matD_norm_inf_complex(&A);
        double herm_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_hermitian(&A, herm_tol)) {
            out = direct_complex_hermitian_machine(&A, want, k_spec);
        } else {
            out = direct_complex_general_machine(&A, want, k_spec);
        }
    } else {
        double norm = matD_norm_inf_real(A.re, A.n);
        double sym_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_real_symmetric(&A, sym_tol)) {
            out = direct_real_sym_machine(&A, want, k_spec);
        } else {
            out = direct_real_general_machine(&A, want, k_spec);
        }
    }
    matD_free(&A);
    return out;
}

/* Top-level "Direct" dispatcher.  Picks the MPFR kernel when any
 * input leaf carries MPFR precision (per common_scan_inexact); falls
 * back to the machine-precision kernel otherwise, or when the MPFR
 * kernel for the matrix shape isn't yet wired. */
static Expr* direct_dispatch(Expr* m, Expr* a, int64_t n,
                               MateigenWant want, Expr* k_spec) {
#ifdef USE_MPFR
    CommonInexactInfo info = common_scan_inexact(m);
    if (a) {
        CommonInexactInfo ia = common_scan_inexact(a);
        if (ia.has_inexact && (!info.has_inexact || ia.min_bits < info.min_bits))
            info = ia;
    }
    if (info.has_inexact && info.min_bits > 53) {
        Expr* out = direct_dispatch_mpfr(m, a, n, (mpfr_prec_t)info.min_bits,
                                          want, k_spec);
        if (out) return out;
        /* MPFR kernel not yet wired for this matrix shape -- fall
         * through to the machine kernel.  The machine kernel will
         * coerce the MPFR cells to doubles via eigen_leaf_to_double,
         * which is the closest behaviour-preserving fallback. */
    }
#endif
    return direct_dispatch_machine(m, a, n, want, k_spec);
}

/* Emit a once-per-method stderr warning that the requested Method is
 * not yet implemented and the symbolic path will be used instead.  The
 * full numerical kernels arrive in Phases 2-5.  Returning the message
 * suppresses repeats within a single process to avoid log spam in
 * tight loops. */
static void eigen_warn_unimplemented_method(MateigenMethod m) {
    static bool warned[MATEIGEN_METHOD_UNKNOWN + 1] = { false };
    if ((int)m < 0 || (int)m > (int)MATEIGEN_METHOD_UNKNOWN) return;
    if (warned[m]) return;
    warned[m] = true;
    const char* name = "?";
    switch (m) {
        case MATEIGEN_DIRECT:          name = "Direct";  break;
        case MATEIGEN_ARNOLDI:         name = "Arnoldi"; break;
        case MATEIGEN_BANDED:          name = "Banded";  break;
        case MATEIGEN_FEAST:           name = "FEAST";   break;
        case MATEIGEN_METHOD_UNKNOWN:  name = "<unknown>"; break;
        default:                       return;
    }
    fprintf(stderr,
        "Eigenvalues::method: Method -> \"%s\" is not yet implemented; "
        "using the symbolic characteristic-polynomial pipeline.\n", name);
}

/* Apply k-spec (Integer k, -k, or UpTo[k]) to vals[0..count].  Returns a
 * freshly allocated trimmed array; the caller frees the originals it owns
 * for the values not selected.  *out_count holds the new size. */
static Expr** eigen_apply_k_spec(Expr** vals, size_t count, Expr* k_spec,
                                  size_t* out_count) {
    size_t result_count = count;
    bool from_end = false;
    if (k_spec) {
        if (k_spec->type == EXPR_INTEGER) {
            int64_t k = k_spec->data.integer;
            if (k >= 0) {
                result_count = ((size_t)k < count) ? (size_t)k : count;
            } else {
                int64_t abs_k = -k;
                result_count = ((size_t)abs_k < count) ? (size_t)abs_k : count;
                from_end = true;
            }
        } else if (k_spec->type == EXPR_FUNCTION
            && k_spec->data.function.head->type == EXPR_SYMBOL
            && k_spec->data.function.head->data.symbol == SYM_UpTo
            && k_spec->data.function.arg_count == 1
            && k_spec->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t k = k_spec->data.function.args[0]->data.integer;
            result_count = ((size_t)k < count) ? (size_t)k : count;
        }
    }
    Expr** result = result_count ? malloc(sizeof(Expr*) * result_count) : NULL;
    if (from_end) {
        size_t start = count - result_count;
        for (size_t i = 0; i < start; i++) expr_free(vals[i]);
        for (size_t i = 0; i < result_count; i++) result[i] = vals[start + i];
    } else {
        for (size_t i = 0; i < result_count; i++) result[i] = vals[i];
        for (size_t i = result_count; i < count; i++) expr_free(vals[i]);
    }
    *out_count = result_count;
    return result;
}

/* Compute the eigenvalue list (padded to n entries with Infinity for the
 * generalised degree-drop case).  Returns a freshly allocated array; the
 * caller owns the entries and the array.  *out_count is set to n on success
 * or to 0 if the calculation cannot be completed (caller should leave the
 * outer Eigenvalues call unevaluated). */
static Expr** eigen_compute_eigenvalues_full(Expr* m, Expr* a,
                                              int64_t n,
                                              bool cubics_radical,
                                              bool quartics_radical,
                                              size_t* out_count) {
    *out_count = 0;
    const char* lam = eigen_lambda_name();
    Expr* poly;
    if (a == NULL) {
        /* Ordinary case: Faddeev-Leverrier in O(n^4) matrix multiplications.
         * Far cheaper than Laplace expansion of the polynomial-entry
         * matrix det(m − λI) (which is O(n!)) once n grows past ~8. */
        poly = eigen_char_poly_faddeev(m, lam, (int)n);
        if (!poly) {
            return NULL;
        }
    } else {
        /* Generalised case: still use Laplace expansion of det(m − λa).
         * Acceptable in practice: generalised eigenproblems in the test
         * corpus are small (≤ 3×3). */
        Expr* M = eigen_build_lambda_matrix(m, a, lam, n);
        poly = eigen_compute_det(M, (int)n);
        expr_free(M);
    }

    Expr* sols = eigen_solve_poly(poly, lam,
                                  cubics_radical, quartics_radical);
    expr_free(poly);
    if (!sols) return NULL;

    size_t val_count = 0;
    Expr** vals = eigen_extract_values(sols, &val_count);
    expr_free(sols);

    /* For generalised eigenvalues, pad short result with Infinity. */
    if ((size_t)n > val_count) {
        Expr** padded = realloc(vals, sizeof(Expr*) * n);
        if (padded) vals = padded;
        else if (!vals) {
            vals = malloc(sizeof(Expr*) * n);
        }
        for (size_t i = val_count; i < (size_t)n; i++) {
            vals[i] = expr_new_symbol("Infinity");
        }
        val_count = (size_t)n;
    }
    *out_count = val_count;
    return vals;
}

Expr* builtin_eigenvalues(Expr* res) {
    EigenOpts opts;
    if (!eigen_parse_args(res, &opts)) return NULL;

    Expr *m, *a; int64_t n;
    if (!eigen_extract_matrix_pair(opts.arg0, &m, &a, &n)) return NULL;

    bool inexact = eigen_matrix_is_inexact(m)
                || (a && eigen_matrix_is_inexact(a));

    /* Numerical Direct dispatch: Automatic and "Direct" route through
     * the hand-rolled hot-path kernels in this file when the input is
     * an inexact ordinary eigenproblem with a supported shape.  The
     * dispatcher returns NULL for shapes that aren't yet wired (step
     * 2.1 covers real symmetric only), in which case we fall through
     * to the symbolic characteristic-polynomial pipeline below. */
    if (inexact && (opts.method == MATEIGEN_AUTOMATIC
                    || opts.method == MATEIGEN_DIRECT)) {
        Expr* out = direct_dispatch(m, a, n,
                                     MATEIGEN_WANT_VALUES,
                                     opts.k_spec);
        if (out) return out;
    }

    /* Method warnings for the not-yet-implemented kernels.  These
     * always fall back to the symbolic path so the call still
     * produces a result. */
    if (inexact && opts.method_given
        && opts.method != MATEIGEN_AUTOMATIC
        && opts.method != MATEIGEN_DIRECT) {
        eigen_warn_unimplemented_method(opts.method);
    }

    size_t val_count = 0;
    Expr** vals = eigen_compute_eigenvalues_full(m, a, n,
        opts.cubics_radical, opts.quartics_radical, &val_count);
    if (!vals || val_count == 0) {
        free(vals);
        return NULL;
    }

    /* For inexact input, chop numerical-noise imaginary parts so real
     * eigenvalues of real matrices appear as real numbers. */
    if (inexact) {
        for (size_t i = 0; i < val_count; i++) {
            Expr* c = eigen_chop(vals[i]);
            expr_free(vals[i]);
            vals[i] = c;
        }
    }

    eigen_sort_by_abs_desc(vals, val_count);

    size_t result_count = val_count;
    Expr** result_vals = eigen_apply_k_spec(vals, val_count, opts.k_spec,
                                             &result_count);
    free(vals);

    Expr* out = expr_new_function(expr_new_symbol("List"),
                                  result_vals, result_count);
    free(result_vals);
    return out;
}

/* ---------------- Eigenvectors ---------------- */

/* Substitute the lambda symbol in `M` with `val`, evaluating each entry.
 * Returns a freshly allocated matrix. */
static Expr* eigen_subst_lambda(Expr* M, const char* lambda_name, Expr* val) {
    /* Use ReplaceAll: M /. lambda_name -> val */
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol(lambda_name), expr_copy(val) }, 2);
    Expr* replaced = eval_and_free(expr_new_function(
        expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(M), rule }, 2));
    return replaced;
}

/* Build the m - λ*a matrix, then substitute λ = `val`.  Returns a numerical
 * (or symbolic-residual) n×n matrix.  Used by the eigenvector routine. */
static Expr* eigen_residual_matrix(Expr* m, Expr* a, Expr* val, int64_t n) {
    /* For ordinary case: just m - val*I, computed directly.
     * For generalised case: m - val*a. */
    const char* lam = eigen_lambda_name();
    Expr* M_lambda = eigen_build_lambda_matrix(m, a, lam, n);
    Expr* M_sub = eigen_subst_lambda(M_lambda, lam, val);
    expr_free(M_lambda);
    /* Try to canonicalise: expand each entry so RowReduce sees them in a
     * normalised form.  Many tests rely on integer-like cancellation here
     * (e.g. (7/2 - 4)*v1 + (1/2)*v3 -> -1/2 v1 + 1/2 v3). */
    return M_sub;
}

/* Build a length-n vector of zeros. */
static Expr* eigen_zero_vector(int64_t n) {
    Expr** v = malloc(sizeof(Expr*) * n);
    for (int64_t i = 0; i < n; i++) v[i] = expr_new_integer(0);
    Expr* out = expr_new_function(expr_new_symbol("List"), v, (size_t)n);
    free(v);
    return out;
}

/* Normalise vector `v` (a List of n elements) by dividing by its Norm.
 * Used by the numerical eigenvector path to emit unit vectors. */
static Expr* eigen_normalize_vector(Expr* v) {
    Expr* norm = eval_and_free(expr_new_function(
        expr_new_symbol("Norm"), (Expr*[]){ expr_copy(v) }, 1));
    /* If Norm is zero or symbolic, skip normalisation. */
    if (is_zero_poly(norm)) { expr_free(norm); return expr_copy(v); }
    Expr* inv = eval_and_free(expr_new_function(
        expr_new_symbol("Power"),
        (Expr*[]){ norm, expr_new_integer(-1) }, 2));
    Expr* scaled = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){ inv, expr_copy(v) }, 2));
    return scaled;
}

Expr* builtin_eigenvectors(Expr* res) {
    EigenOpts opts;
    if (!eigen_parse_args(res, &opts)) return NULL;

    Expr *m, *a; int64_t n;
    if (!eigen_extract_matrix_pair(opts.arg0, &m, &a, &n)) return NULL;

    bool inexact = eigen_matrix_is_inexact(m)
                || (a && eigen_matrix_is_inexact(a));

    /* Numerical Direct dispatch (mirrors builtin_eigenvalues). */
    if (inexact && (opts.method == MATEIGEN_AUTOMATIC
                    || opts.method == MATEIGEN_DIRECT)) {
        Expr* out = direct_dispatch(m, a, n,
                                     MATEIGEN_WANT_VECTORS,
                                     opts.k_spec);
        if (out) return out;
    }

    if (inexact && opts.method_given
        && opts.method != MATEIGEN_AUTOMATIC
        && opts.method != MATEIGEN_DIRECT) {
        eigen_warn_unimplemented_method(opts.method);
    }

    /* Inexact input: rationalise once up front so we can perform the
     * eigenvalue / null-space arithmetic in exact form -- numerical
     * RowReduce on a (m − λI) substituted with inexact λ would zero out
     * the rank-defect we need to discover the eigenvector.  Numericalize
     * + normalise the result at the end. */
    Expr* m_orig = m;
    Expr* a_orig = a;
    long prec_bits = 53;
    Expr* m_rat = NULL;
    Expr* a_rat = NULL;
    if (inexact) {
        CommonInexactInfo info = common_scan_inexact(m);
        if (a) {
            CommonInexactInfo info_a = common_scan_inexact(a);
            if (info_a.has_inexact
                && (!info.has_inexact || info_a.min_bits < info.min_bits)) {
                info = info_a;
            }
        }
        prec_bits = info.min_bits ? info.min_bits : 53;
        m_rat = common_rationalize_input(m, prec_bits);
        if (a) a_rat = common_rationalize_input(a, prec_bits);
        m = m_rat;
        if (a) a = a_rat;
    }

    /* Compute eigenvalues in the same arithmetic form as the matrix. */
    size_t val_count = 0;
    Expr** vals = eigen_compute_eigenvalues_full(m, a, n,
        opts.cubics_radical, opts.quartics_radical, &val_count);
    if (!vals || val_count == 0) {
        free(vals);
        if (m_rat) expr_free(m_rat);
        if (a_rat) expr_free(a_rat);
        return NULL;
    }
    eigen_sort_by_abs_desc(vals, val_count);

    /* Collect eigenvectors -- traverse eigenvalues in order, collapsing
     * runs of equal values into a single null-space computation that
     * yields up to `multiplicity` vectors. */
    Expr** vectors = malloc(sizeof(Expr*) * n);
    size_t vc = 0;

    size_t i = 0;
    while (i < val_count && vc < (size_t)n) {
        /* Determine run of equal eigenvalues. */
        size_t j = i + 1;
        while (j < val_count && expr_eq(vals[j], vals[i])) j++;
        int64_t mult = (int64_t)(j - i);

        Expr* val = vals[i];

        /* Special handling: Infinity eigenvalues correspond to the null
         * space of `a` (the generalised pencil's "infinite" branch). */
        bool is_inf = (val->type == EXPR_SYMBOL
                       && val->data.symbol == SYM_Infinity);

        Expr* residual = NULL;
        if (is_inf && a) {
            residual = expr_copy(a);
        } else {
            residual = eigen_residual_matrix(m, a, val, n);
        }

        size_t basis_count = 0;
        Expr** basis = eigen_null_space(residual, (int)n, &basis_count);
        expr_free(residual);

        /* Take up to `mult` vectors. */
        size_t take = (basis_count < (size_t)mult) ? basis_count : (size_t)mult;
        for (size_t k = 0; k < take && vc < (size_t)n; k++) {
            Expr* v = basis[k];
            if (inexact) v = eigen_normalize_vector(v);
            else v = expr_copy(v);
            vectors[vc++] = v;
        }
        /* If null-space gives fewer vectors than multiplicity, the matrix
         * is defective for this eigenvalue: pad the shortfall in-line with
         * zero vectors so the i-th eigenvector still lines up positionally
         * with the i-th eigenvalue. */
        for (size_t k = take; k < (size_t)mult && vc < (size_t)n; k++) {
            vectors[vc++] = eigen_zero_vector(n);
        }
        for (size_t k = 0; k < basis_count; k++) expr_free(basis[k]);
        free(basis);

        i = j;
    }

    /* Pad with zero vectors. */
    while (vc < (size_t)n) {
        vectors[vc++] = eigen_zero_vector(n);
    }

    /* Free unused eigenvalues. */
    for (size_t k = 0; k < val_count; k++) expr_free(vals[k]);
    free(vals);

    size_t result_count = vc;
    Expr** result_vecs = eigen_apply_k_spec(vectors, vc, opts.k_spec,
                                             &result_count);
    free(vectors);

    Expr* out = expr_new_function(expr_new_symbol("List"),
                                  result_vecs, result_count);
    free(result_vecs);

    /* Inexact path: numericalize the eigenvectors back to the original
     * precision.  We did all of the null-space and normalisation work in
     * exact (rational) form so the rank-deficient direction was correctly
     * captured.  Free the rationalised matrices we created. */
    if (inexact) {
        Expr* numeric = common_numericalize_result(out, prec_bits);
        expr_free(out);
        out = numeric;
        if (m_rat) expr_free(m_rat);
        if (a_rat) expr_free(a_rat);
    }
    (void)m_orig; (void)a_orig;
    return out;
}

void mateigen_init(void) {
    symtab_add_builtin("Eigenvalues", builtin_eigenvalues);
    symtab_get_def("Eigenvalues")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Eigenvectors", builtin_eigenvectors);
    symtab_get_def("Eigenvectors")->attributes |= ATTR_PROTECTED;
}
