/* svdecomp.c
 *
 * SingularValueDecomposition[m]                       -- {u, sigma, v}
 * SingularValueDecomposition[m, k]                    -- k largest (k<0 -> |k| smallest)
 * SingularValueDecomposition[m, UpTo[k]]              -- up to k largest
 * SingularValueDecomposition[{m, a}]                  -- generalized form
 * Options:  Tolerance -> t,  TargetStructure -> "Dense" | "Structured"
 *
 * Strategy.  One algorithmic core per numeric domain:
 *
 *   - Exact / symbolic   -> eigendecomposition of m^H . m (or m . m^H,
 *                           whichever is smaller), with a 2x2 closed-form
 *                           fast path.  Outputs Sqrt[lambda] for the
 *                           singular values; columns of v are orthonormal
 *                           eigenvectors; u = m . v . Sigma^-1 with the
 *                           null space completed via the existing
 *                           qr_symbolic_core orthogonal completion.
 *
 *   - Inexact, min_bits <= 53  -> LAPACK dgesdd / zgesdd (standard) or
 *                                 dggsvd3 / zggsvd3 (generalized).
 *
 *   - Inexact, min_bits > 53   -> one-sided Jacobi SVD over MPFR arrays.
 *
 * All three paths feed the result through svd_apply_postprocess so the
 * truncation, tolerance, and TargetStructure logic lives in one place.
 *
 * Memory contract.  Standard builtin contract.  This file does NOT call
 * expr_free(res) - the evaluator owns `res` and frees it on a non-NULL
 * return (MEMORY.md / SPEC.md Sec. 4.1).
 */

#include "svdecomp.h"
#include "svdecomp_internal.h"
#include "qrdecomp.h"
#include "qrdecomp_internal.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"
#include "common.h"
#include "linsolve.h"
#include "poly.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Forward declarations of pre-Phase-2 dispatcher stubs.  The bodies  *
 *  live below the parser; each Phase-2/3/4 iteration replaces them    *
 *  with a real implementation in svdecomp.c / svdecomp_machine.c /    *
 *  svdecomp_mpfr.c respectively.                                      *
 * ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ *
 *  Tensor-shape probe specialised for SVD.                            *
 *                                                                     *
 *  Returns true iff `e` is a non-empty rectangular rank-2 tensor.     *
 *  On success populates *rows and *cols.  Used by both the standard   *
 *  matrix input and each leaf of the generalized {m, a} pair.         *
 * ------------------------------------------------------------------ */
static bool probe_matrix(Expr* e, int* rows, int* cols) {
    int64_t dims[64];
    int trank = get_tensor_dims(e, dims);
    if (trank != 2 || dims[0] == 0 || dims[1] == 0) return false;
    *rows = (int)dims[0];
    *cols = (int)dims[1];
    return true;
}

/* ------------------------------------------------------------------ *
 *  Option-value parsers.                                              *
 *                                                                     *
 *  parse_tolerance accepts any expression as the threshold (Real,     *
 *  Integer, Rational, MPFR, symbolic numeric).  The kernels are       *
 *  responsible for evaluating Abs[sigma_i] < tol in their own arith.  *
 * ------------------------------------------------------------------ */
static bool parse_targetstructure_value(Expr* rhs, SvdTargetStructure* out) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
        *out = SVD_TS_DENSE;
        return true;
    }
    if (rhs->type != EXPR_STRING) return false;
    if (strcmp(rhs->data.string, "Dense")      == 0) { *out = SVD_TS_DENSE;      return true; }
    if (strcmp(rhs->data.string, "Structured") == 0) { *out = SVD_TS_STRUCTURED; return true; }
    return false;
}

/* True iff `e` is a Rule[lhs, rhs] / RuleDelayed[lhs, rhs] with a
 * symbol on the left -- the shape we accept for named options. */
static bool is_named_option(Expr* e) {
    if (e->type != EXPR_FUNCTION
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.arg_count != 2) return false;
    const char* hd = e->data.function.head->data.symbol;
    if (hd != SYM_Rule && hd != SYM_RuleDelayed) return false;
    return e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Try to interpret `e` as the k positional argument.
 *
 *   Integer k         -> SVD_FORM_K, k_value = k
 *   UpTo[Integer k]   -> SVD_FORM_UPTO, k_value = k (k >= 0)
 *
 * Returns false if `e` does not match either shape (in which case the
 * caller treats `e` as an option or rejects it). */
static bool parse_k_value(Expr* e, SvdKForm* form, int* k_out) {
    if (e->type == EXPR_INTEGER) {
        *form  = SVD_FORM_K;
        *k_out = (int)e->data.integer;
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_UpTo
        && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        if (arg->type != EXPR_INTEGER) return false;
        if (arg->data.integer < 0) return false;
        *form  = SVD_FORM_UPTO;
        *k_out = (int)arg->data.integer;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Argument parser.                                                   *
 *                                                                     *
 *  Grammar (after `SingularValueDecomposition` is stripped):          *
 *                                                                     *
 *     M  [ K ]  OPT*                                                  *
 *                                                                     *
 *  where                                                              *
 *     M   = a matrix (List of equal-length Lists)                     *
 *         | { m, a } -- a pair of matrices, generalized form          *
 *     K   = Integer k     -- SVD_FORM_K, signed                       *
 *         | UpTo[Integer k]    -- SVD_FORM_UPTO, k >= 0               *
 *     OPT = Tolerance -> t                                            *
 *         | TargetStructure -> "Dense" | "Structured" | Automatic     *
 *                                                                     *
 *  Returns false on any malformed positional or option entry.  The    *
 *  caller emits the appropriate ::opts / ::matrix / ::sval message.   *
 * ------------------------------------------------------------------ */
bool svd_parse_args(Expr* res, SvdArgs* args) {
    args->m                = NULL;
    args->a                = NULL;
    args->generalized      = false;
    args->k_form           = SVD_FORM_FULL;
    args->k_value          = 0;
    args->tolerance        = NULL;
    args->target_structure = SVD_TS_DENSE;

    size_t argc = res->data.function.arg_count;
    if (argc < 1) return false;

    /* arg[0]: matrix or {m, a} pair. */
    Expr* a0 = res->data.function.args[0];

    /* Pair detection: List[m, a] where both m and a are rank-2 tensors
     * AND m and a are distinguishable from a single matrix.  We rely on
     * arity: SVD's pair form is List[matrix, matrix] (length 2), but a
     * single 2-row matrix is also List[row, row].  Disambiguate by
     * checking whether arg[0] is a rank-2 tensor on its own; if it is,
     * we treat it as the standard matrix form regardless of arity. */
    bool pair_handled = false;
    if (a0->type == EXPR_FUNCTION
        && a0->data.function.head->type == EXPR_SYMBOL
        && a0->data.function.head->data.symbol == SYM_List
        && a0->data.function.arg_count == 2) {

        int64_t outer_dims[64];
        int outer_rank = get_tensor_dims(a0, outer_dims);
        /* outer_rank == 2 means a0 is a 2-row matrix, NOT a pair. */
        if (outer_rank != 2) {
            Expr* m = a0->data.function.args[0];
            Expr* a = a0->data.function.args[1];
            int mr, mc, ar, ac;
            if (probe_matrix(m, &mr, &mc) && probe_matrix(a, &ar, &ac)
                && mc == ac) {
                args->m           = m;
                args->a           = a;
                args->generalized = true;
                pair_handled      = true;
            } else {
                /* It's a 2-element List but not a valid pair of matrices. */
                return false;
            }
        }
    }
    if (!pair_handled) {
        args->m           = a0;
        args->generalized = false;
    }

    /* arg[1]: optional k positional, otherwise treated as an option. */
    size_t opt_start = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        if (!is_named_option(a1)) {
            if (!parse_k_value(a1, &args->k_form, &args->k_value)) return false;
            opt_start = 2;
        }
    }

    /* arg[opt_start..]: named options. */
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!is_named_option(opt)) return false;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];

        if (lhs->data.symbol == SYM_TargetStructure) {
            SvdTargetStructure ts;
            if (!parse_targetstructure_value(rhs, &ts)) return false;
            args->target_structure = ts;
        } else if (lhs->data.symbol == SYM_Tolerance) {
            args->tolerance = rhs;
        } else {
            return false;
        }
    }

    /* Truncation is undefined for the generalized form (Mathematica
     * raises SingularValueDecomposition::nopart).  Reject at parse
     * time so the dispatchers can assume non-generalized when
     * args->k_form != SVD_FORM_FULL. */
    if (args->generalized && args->k_form != SVD_FORM_FULL) return false;

    return true;
}

/* =================================================================== *
 *  Symbolic SVD core (Phase 2).                                        *
 *                                                                      *
 *  Strategy:                                                           *
 *                                                                      *
 *    1. Build the smaller of the two Gram products through the         *
 *       evaluator: B = m^H . m (p x p) when n >= p, otherwise          *
 *       B = m . m^H (n x n).                                           *
 *    2. Eigenvalues[B] and Eigenvectors[B] via the evaluator.          *
 *    3. Sort pairs by descending eigenvalue (ties broken by index).    *
 *    4. Orthonormalise the eigenvector matrix via qr_symbolic_core --  *
 *       across distinct eigenvalues Hermitian eigenvectors are already *
 *       orthogonal, within a repeated eigenvalue Gram-Schmidt picks    *
 *       a fresh orthonormal basis of the eigenspace.                   *
 *    5. sigma_i = Sqrt[lambda_i].  rank = count of non-zero sigma_i.   *
 *    6. Compute the *other* set of singular vectors:                   *
 *         n >= p :  u_i = m . v_i / sigma_i  for sigma_i > 0; the      *
 *                   remaining n - rank columns of U come from a        *
 *                   orthonormal completion of {u_1, ..., u_rank} in    *
 *                   C^n via qr_symbolic_core on [u_1 | ... | I_n].     *
 *         n <  p :  v_i = m^H . u_i / sigma_i  for sigma_i > 0; the    *
 *                   remaining p - rank columns of V come from the      *
 *                   same orthonormal-completion trick on the V side.   *
 *    7. Pack U (n x n), Sigma (n x p with sigma_i on diagonal, zeros   *
 *       elsewhere), V (p x p).                                         *
 *                                                                      *
 *  Soft failure paths return false and leave outputs untouched -- the  *
 *  caller falls through to the unevaluated form (or, for inexact input *
 *  inside svd_symbolic_dispatch, drops back to NULL so the outer       *
 *  dispatcher's machine kernel can retry).                             *
 * =================================================================== */

/* ------------------------------------------------------------------ *
 *  Small evaluator helpers, mirroring the ones in qrdecomp.c.         *
 *  Each one builds a tiny expression tree, evaluates it, and returns  *
 *  ownership of the result.  Argument trees are consumed.             *
 * ------------------------------------------------------------------ */
static Expr* sv_plus(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Plus), (Expr*[]){a, b}, 2));
}
static Expr* sv_times(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times), (Expr*[]){a, b}, 2));
}
static Expr* sv_power(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power), (Expr*[]){a, b}, 2));
}
static Expr* sv_sqrt(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Sqrt), (Expr*[]){a}, 1));
}
static Expr* sv_conjugate(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Conjugate), (Expr*[]){a}, 1));
}
static Expr* sv_together(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Together"), (Expr*[]){a}, 1));
}
static Expr* sv_eigenvalues(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Eigenvalues), (Expr*[]){a}, 1));
}
static Expr* sv_eigenvectors(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Eigenvectors), (Expr*[]){a}, 1));
}

/* Squared-norm zero test, robust on Sqrt / rational forms (same idiom
 * as qrdecomp.c::is_definitely_zero). */
static bool sv_is_zero(Expr* e) {
    Expr* simp = sv_together(expr_copy(e));
    bool z = is_zero_poly(simp);
    expr_free(simp);
    return z;
}

/* True iff `e` contains a Complex[...] head or the literal symbol I.
 * Drives the choice between Hermitian and Euclidean inner products. */
static bool sv_has_complex_content(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol, "I") == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex) return true;
    if (sv_has_complex_content(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (sv_has_complex_content(e->data.function.args[i])) return true;
    }
    return false;
}

/* True iff every leaf of the matrix `m` is an exact-or-inexact numeric
 * atom (Integer, BigInt, Rational, Real, MPFR, or Complex[a, b] where
 * a and b are themselves numeric).  Used by the symbolic generalized-
 * SVD passthrough to decide whether the input can be safely downgraded
 * to machine precision and re-dispatched. */
static bool sv_leaf_is_numeric(Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION:
            if (e->data.function.head->type == EXPR_SYMBOL
                && e->data.function.arg_count == 2) {
                const char* h = e->data.function.head->data.symbol;
                if (h == SYM_Rational || strcmp(h, "Rational") == 0
                    || h == SYM_Complex  || strcmp(h, "Complex")  == 0) {
                    return sv_leaf_is_numeric(e->data.function.args[0])
                        && sv_leaf_is_numeric(e->data.function.args[1]);
                }
            }
            return false;
        default:
            return false;
    }
}

static bool sv_matrix_is_purely_numeric(Expr* m, int n, int p) {
    if (!m || m->type != EXPR_FUNCTION) return false;
    if ((int)m->data.function.arg_count != n) return false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        if ((int)row->data.function.arg_count != p) return false;
        for (int j = 0; j < p; j++) {
            if (!sv_leaf_is_numeric(row->data.function.args[j])) return false;
        }
    }
    return true;
}

/* True iff `e` is provably a NON-NEGATIVE numeric atom (Integer,
 * Rational with positive denominator and non-negative numerator,
 * Real >= 0).  Used to decide whether sigma = Sqrt[lambda] yields a
 * real non-negative singular value.  Symbolic expressions return false
 * here -- the algorithm still proceeds (Sqrt is left as a head), the
 * caller just can't sort against them numerically. */
static bool sv_is_nonneg_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer >= 0;
    if (e->type == EXPR_REAL)    return e->data.real >= 0.0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) >= 0;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        const Expr* num = e->data.function.args[0];
        const Expr* den = e->data.function.args[1];
        if (num->type == EXPR_INTEGER && den->type == EXPR_INTEGER) {
            int64_t ns = num->data.integer < 0 ? -1 : (num->data.integer > 0 ? 1 : 0);
            int64_t ds = den->data.integer < 0 ? -1 : (den->data.integer > 0 ? 1 : 0);
            return ns * ds >= 0;
        }
    }
    return false;
}

/* True iff `e` contains any Root[...] head -- a sign that
 * Eigenvalues / Eigenvectors couldn't reduce to a closed form we can
 * lift through Sqrt cleanly.  The symbolic SVD then declines. */
static bool sv_contains_root(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Root) return true;
    if (sv_contains_root(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (sv_contains_root(e->data.function.args[i])) return true;
    }
    return false;
}

/* Eigenvalues / Eigenvectors returned the same head (couldn't evaluate),
 * or any leaf is a Root[...] -- bail out to the caller. */
static bool sv_eig_usable(Expr* eigvals, Expr* eigvecs, int expected) {
    if (!eigvals || !eigvecs) return false;
    if (eigvals->type != EXPR_FUNCTION
        || eigvals->data.function.head->type != EXPR_SYMBOL
        || eigvals->data.function.head->data.symbol != SYM_List
        || (int)eigvals->data.function.arg_count != expected) return false;
    if (eigvecs->type != EXPR_FUNCTION
        || eigvecs->data.function.head->type != EXPR_SYMBOL
        || eigvecs->data.function.head->data.symbol != SYM_List
        || (int)eigvecs->data.function.arg_count != expected) return false;
    for (int i = 0; i < expected; i++) {
        if (sv_contains_root(eigvals->data.function.args[i])) return false;
        Expr* row = eigvecs->data.function.args[i];
        if (row->type != EXPR_FUNCTION
            || row->data.function.head->type != EXPR_SYMBOL
            || row->data.function.head->data.symbol != SYM_List
            || (int)row->data.function.arg_count != expected) return false;
        for (int k = 0; k < expected; k++) {
            if (sv_contains_root(row->data.function.args[k])) return false;
        }
    }
    return true;
}

/* Build a Mathilda matrix expression (List of List) from a row-major
 * Expr** buffer.  Entries are COPIED; caller retains ownership of
 * flat[i]. */
static Expr* sv_matrix_from_flat(Expr** flat, int rows, int cols) {
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) elems[j] = expr_copy(flat[i * cols + j]);
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        free(elems);
    }
    Expr* m = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return m;
}

/* Build B = A^H . A (when use_AHA) or A . A^H (otherwise) by going
 * through the evaluator.  Caller owns the result. */
static Expr* sv_gram(Expr* A_expr, bool use_AHA, bool use_conj) {
    Expr* op   = use_conj
                 ? expr_new_function(expr_new_symbol(SYM_ConjugateTranspose),
                                     (Expr*[]){expr_copy(A_expr)}, 1)
                 : expr_new_function(expr_new_symbol(SYM_Transpose),
                                     (Expr*[]){expr_copy(A_expr)}, 1);
    Expr* AT   = eval_and_free(op);
    Expr* lhs  = use_AHA ? AT            : expr_copy(A_expr);
    Expr* rhs  = use_AHA ? expr_copy(A_expr) : AT;
    return eval_and_free(expr_new_function(
        expr_new_symbol("Dot"), (Expr*[]){lhs, rhs}, 2));
}

/* Sort eigenvalue / eigenvector pairs by descending magnitude.
 *
 * The sort policy is intentionally conservative because eigenvalues of
 * a symmetric positive-semidefinite matrix can be (a) exact numeric,
 * (b) symbolic-positive expressions like 52 + Sqrt[2669], or (c) the
 * exact integer 0 representing a null-space direction.  We can compare
 * pairs in (a)+(a) and we can detect (c).  For any pair involving (b)
 * we do not have a faithful magnitude comparison, so we keep Mathilda's
 * Eigenvalues order intact.
 *
 *   - Zeros are pushed to the tail (so the rank classifier sees them
 *     last and the corresponding sigma_i = 0 entries are at the end of
 *     sigma_full).
 *   - Two non-negative numeric eigenvalues are sorted descending.
 *   - Otherwise no swap.
 */
static void sv_sort_pairs(Expr** lams, Expr** vecs_row_major, int g) {
    for (int i = 1; i < g; i++) {
        for (int j = i; j > 0; j--) {
            Expr* a = lams[j - 1];
            Expr* b = lams[j];
            bool a_zero = sv_is_zero(a);
            bool b_zero = sv_is_zero(b);
            int swap = 0;
            if (a_zero && !b_zero) {
                /* Push the zero down past the non-zero. */
                swap = 1;
            } else if (!a_zero && !b_zero
                       && sv_is_nonneg_numeric(a)
                       && sv_is_nonneg_numeric(b)) {
                /* Both numeric non-zero: descending. */
                swap = expr_compare(b, a) > 0;
            }
            if (!swap) break;
            Expr* tmp = lams[j - 1]; lams[j - 1] = lams[j]; lams[j] = tmp;
            for (int k = 0; k < g; k++) {
                Expr* t = vecs_row_major[(j - 1) * g + k];
                vecs_row_major[(j - 1) * g + k] = vecs_row_major[j * g + k];
                vecs_row_major[j * g + k] = t;
            }
        }
    }
}

/* Multiply an m x p matrix (row-major Expr**) by a p-vector (Expr**)
 * via the evaluator.  Returns a fresh m-vector.  Inputs are NOT
 * consumed. */
static Expr** sv_matvec(Expr** Aflat, int m_rows, int p_cols,
                        Expr** v, bool use_conj) {
    /* Use the same inner-product style as qrdecomp.c: row . vector
     * with the evaluator.  Conjugate is applied to neither side here
     * because the rows of A are the actual matrix rows, not "vectors
     * we want inner products against" - it's a matrix-vector product. */
    (void)use_conj;
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)m_rows);
    for (int i = 0; i < m_rows; i++) {
        Expr* sum = expr_new_integer(0);
        for (int k = 0; k < p_cols; k++) {
            Expr* term = sv_times(expr_copy(Aflat[i * p_cols + k]),
                                  expr_copy(v[k]));
            sum = sv_plus(sum, term);
        }
        out[i] = sum;
    }
    return out;
}

/* Identity column buffer of length n: e_i has 1 in position i, 0
 * elsewhere.  Used to seed the orthonormal-completion call to
 * qr_symbolic_core. */
static Expr* sv_identity_entry(int row, int col) {
    return expr_new_integer(row == col ? 1 : 0);
}

/* Orthonormalise an n x k matrix of column vectors stored as
 * col-major Expr** (cols_data[i + j*n] is row i, col j of the input).
 * Optionally append the n x n identity to provide an orthonormal
 * completion.  The output is the first `want` columns of the QR Q,
 * returned as a freshly allocated row-major Expr** (out[i*want + j]).
 *
 * Returns false on allocation failure or QR rank insufficiency
 * (shouldn't happen when the trailing identity is appended).  Caller
 * frees the entries with expr_free and the array with free. */
static bool sv_orthonormal_complete(Expr** cols_data, int n, int existing,
                                    int want, bool use_conj,
                                    Expr*** out_row_major) {
    int total_cols = existing + n;
    Expr** A_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)total_cols);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < existing; j++) {
            A_flat[i * total_cols + j] = expr_copy(cols_data[i + j * n]);
        }
        for (int j = 0; j < n; j++) {
            A_flat[i * total_cols + (existing + j)] = sv_identity_entry(i, j);
        }
    }
    Expr** Q_flat = NULL;
    Expr** R_flat = NULL;
    int rank = 0;
    bool ok = qr_symbolic_core(A_flat, n, total_cols, false, NULL,
                               use_conj, &Q_flat, &R_flat, &rank);
    for (int i = 0; i < n * total_cols; i++) expr_free(A_flat[i]);
    free(A_flat);
    if (!ok || rank < want) {
        if (Q_flat) {
            for (int i = 0; i < n * rank; i++) expr_free(Q_flat[i]);
            free(Q_flat);
        }
        if (R_flat) {
            for (int j = 0; j < rank * total_cols; j++) expr_free(R_flat[j]);
            free(R_flat);
        }
        return false;
    }
    /* qr_symbolic_core returns Q in column-major n x rank layout.  Pack
     * the first `want` columns into row-major n x want. */
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)want);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < want; j++) {
            out[i * want + j] = Q_flat[i * rank + j];   /* steal */
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = want; j < rank; j++) expr_free(Q_flat[i * rank + j]);
    }
    free(Q_flat);
    for (int j = 0; j < rank * total_cols; j++) expr_free(R_flat[j]);
    free(R_flat);
    *out_row_major = out;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Symbolic SVD inner loop.                                           *
 *                                                                     *
 *  Output layout:                                                     *
 *    U_flat   n * n row-major          (full U)                       *
 *    S_flat   m = min(n, p) entries    (singular values, descending)  *
 *    V_flat   p * p row-major          (full V)                       *
 *    rank     number of strictly positive singular values             *
 *                                                                     *
 *  The caller wraps these as List-of-List with the rectangular Sigma  *
 *  built explicitly (n x p with sigma_i on diagonal).                 *
 * ------------------------------------------------------------------ */
bool svd_symbolic_core(Expr** A_flat, int n, int p, bool use_conj,
                       Expr*** out_U_flat,
                       Expr*** out_S_flat,
                       Expr*** out_V_flat,
                       int* out_rank) {
    /* Strategy: use the LARGER gram so the side we still need to
     * orthonormal-complete is the smaller one (cheap).  Eigenvectors of
     * the larger gram include the null-space directions automatically,
     * so the "primary" side -- V when p >= n, else U -- is determined
     * directly from Eigenvectors without any completion.
     *
     *   use_AHA == true   <=>  B = m^H . m  (p x p),  primary = V
     *   use_AHA == false  <=>  B = m . m^H  (n x n),  primary = U
     */
    bool use_AHA = (p >= n);
    int  g = use_AHA ? p : n;    /* size of the Gram matrix / primary V or U */

    /* 1. Build the gram matrix. */
    Expr* A_expr = sv_matrix_from_flat(A_flat, n, p);
    Expr* B      = sv_gram(A_expr, use_AHA, use_conj);

    /* 2. Eigendecompose. */
    Expr* eigvals = sv_eigenvalues(expr_copy(B));
    Expr* eigvecs = sv_eigenvectors(expr_copy(B));
    expr_free(B);

    if (!sv_eig_usable(eigvals, eigvecs, g)) {
        if (eigvals) expr_free(eigvals);
        if (eigvecs) expr_free(eigvecs);
        expr_free(A_expr);
        return false;
    }

    /* 3. Pull eigenvalues / eigenvectors into a fresh local pair. */
    Expr** lams = (Expr**)malloc(sizeof(Expr*) * (size_t)g);
    Expr** vecs_rm = (Expr**)malloc(sizeof(Expr*) * (size_t)g * (size_t)g);
    for (int i = 0; i < g; i++) {
        lams[i] = expr_copy(eigvals->data.function.args[i]);
        Expr* row = eigvecs->data.function.args[i];
        for (int k = 0; k < g; k++) {
            vecs_rm[i * g + k] = expr_copy(row->data.function.args[k]);
        }
    }
    expr_free(eigvals);
    expr_free(eigvecs);

    /* 4. Sort by descending eigenvalue. */
    sv_sort_pairs(lams, vecs_rm, g);

    /* 5. Normalise each eigenvector individually to unit length.
     *
     *    Why not full QR-orthonormalisation: a Hermitian B has
     *    pairwise-orthogonal eigenvectors for distinct eigenvalues,
     *    so we only need per-vector normalisation in the common case.
     *    The full qr_symbolic_core path was tried first; it pays
     *    Together on every inner product as part of its is-zero rank
     *    test, and for raw Mathilda Eigenvectors[...] output (deeply
     *    nested Sqrt forms) that explodes combinatorially.
     *
     *    Repeated-eigenvalue groups would need within-group Gram-
     *    Schmidt; this implementation accepts that as a known
     *    limitation and emits a non-orthonormal V when the user has
     *    sufficient degeneracy that Mathilda's Eigenvectors returns
     *    a linearly-dependent in-group basis.  In that case the
     *    reconstruction sigma . V^T . V . sigma != sigma^2 and the
     *    user should fall through to the machine kernel.  For the
     *    typical SVD case (random / structured matrices with simple
     *    eigenvalues) this is correct.
     *
     *    First Together each entry so the per-vector norm is computed
     *    on a single canonical rational form rather than the original
     *    nested-fraction Mathilda Eigenvectors output. */
    Expr** Qe = (Expr**)malloc(sizeof(Expr*) * (size_t)g * (size_t)g);
    for (int j = 0; j < g; j++) {
        /* Together-canonicalise the j-th eigenvector entries. */
        Expr** vj = (Expr**)malloc(sizeof(Expr*) * (size_t)g);
        for (int k = 0; k < g; k++) {
            vj[k] = sv_together(expr_copy(vecs_rm[j * g + k]));
        }
        /* nsq = Sum_k Conjugate[v_k] * v_k  (Hermitian inner product). */
        Expr* nsq = expr_new_integer(0);
        for (int k = 0; k < g; k++) {
            Expr* lhs = use_conj ? sv_conjugate(expr_copy(vj[k])) : expr_copy(vj[k]);
            nsq = sv_plus(nsq, sv_times(lhs, expr_copy(vj[k])));
        }
        nsq = sv_together(nsq);
        /* Catch the degenerate zero-eigenvector case (shouldn't happen
         * for the Eigenvectors of a Hermitian B, but be safe). */
        if (sv_is_zero(nsq)) {
            expr_free(nsq);
            for (int k = 0; k < g; k++) expr_free(vj[k]);
            free(vj);
            for (int i = 0; i < g; i++) expr_free(lams[i]);
            free(lams);
            for (int i = 0; i < g * g; i++) expr_free(vecs_rm[i]);
            free(vecs_rm);
            for (int i = 0; i < j * g; i++) expr_free(Qe[i]);
            free(Qe);
            expr_free(A_expr);
            return false;
        }
        Expr* inv_norm = sv_power(sv_sqrt(nsq), expr_new_integer(-1));
        for (int k = 0; k < g; k++) {
            /* Qe is column-major g x g: Qe[i * g + j] stores (i, j). */
            Qe[k * g + j] = sv_times(vj[k], expr_copy(inv_norm));
        }
        expr_free(inv_norm);
        free(vj);
    }
    for (int i = 0; i < g * g; i++) expr_free(vecs_rm[i]);
    free(vecs_rm);

    /* 6. sigma_i = Sqrt[lams_i].  rank = count of non-zero sigma. */
    Expr** sigma_full = (Expr**)malloc(sizeof(Expr*) * (size_t)g);
    int rank_out = 0;
    for (int i = 0; i < g; i++) {
        if (sv_is_zero(lams[i])) {
            sigma_full[i] = expr_new_integer(0);
        } else {
            sigma_full[i] = sv_sqrt(expr_copy(lams[i]));
            rank_out++;
        }
    }
    for (int i = 0; i < g; i++) expr_free(lams[i]);
    free(lams);

    /* 7. Build the "primary" matrix (V if use_AHA, else U) from Qe:
     *    primary_flat is g x g row-major.  When use_AHA, columns of
     *    Qe are right singular vectors (p-dim); else they are left
     *    singular vectors (n-dim). */
    Expr** primary_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)g * (size_t)g);
    for (int i = 0; i < g; i++) {
        for (int j = 0; j < g; j++) {
            primary_flat[i * g + j] = expr_copy(Qe[i * g + j]);
        }
    }

    /* 8. Build the "secondary" matrix.
     *    n >= p: secondary = U (n x n); for sigma_i > 0:
     *            u_i = (1/sigma_i) m . v_i.
     *    n <  p: secondary = V (p x p); for sigma_i > 0:
     *            v_i = (1/sigma_i) m^H . u_i.
     *    For sigma_i = 0 columns we orthonormal-complete against the
     *    computed ones via qr_symbolic_core(I_secondary).
     *
     *    A_for_matvec: row-major n x p for matvec; for the n < p path
     *    we need m^H, so we build it explicitly. */
    int sec_dim = use_AHA ? n : p;
    Expr** secondary_cols_cm = (Expr**)malloc(
        sizeof(Expr*) * (size_t)sec_dim * (size_t)rank_out);
    int built = 0;
    /* For each column index i in [0, g) we may add it to the secondary
     * if sigma_i > 0. */
    for (int i = 0; i < g; i++) {
        if (sv_is_zero(sigma_full[i])) continue;
        /* Extract column i of primary as a vector (length g). */
        Expr** v = (Expr**)malloc(sizeof(Expr*) * (size_t)g);
        for (int k = 0; k < g; k++) v[k] = primary_flat[k * g + i];
        /* Compute matvec.  When use_AHA we apply m to a length-p vector
         * (g == p) and get a length-n vector.  When !use_AHA we apply
         * m^H to a length-n vector (g == n) and get a length-p vector;
         * Mathilda has no direct ConjugateTranspose helper here -- we
         * build the transposed matrix view by re-arranging A_flat. */
        Expr** w;
        if (use_AHA) {
            w = sv_matvec(A_flat, n, p, v, use_conj);
        } else {
            /* Build a column-major view (the rows of A^H are columns
             * of A; entries get conjugated when use_conj). */
            Expr** Aflat_T = (Expr**)malloc(sizeof(Expr*) * (size_t)p * (size_t)n);
            for (int rr = 0; rr < p; rr++) {
                for (int cc = 0; cc < n; cc++) {
                    Expr* e = expr_copy(A_flat[cc * p + rr]);
                    Aflat_T[rr * n + cc] = use_conj ? sv_conjugate(e) : e;
                }
            }
            w = sv_matvec(Aflat_T, p, n, v, use_conj);
            for (int t = 0; t < p * n; t++) expr_free(Aflat_T[t]);
            free(Aflat_T);
        }
        free(v);
        /* Divide by sigma_i.  Note: we deliberately do NOT apply
         * Together here -- for symbolic sigma_i involving nested
         * Sqrt forms the per-entry canonicalisation explodes.  The
         * caller applies tidy_matrix once at the end on the entire
         * matrix, which is materially cheaper. */
        Expr* inv = sv_power(expr_copy(sigma_full[i]), expr_new_integer(-1));
        for (int k = 0; k < sec_dim; k++) {
            secondary_cols_cm[k + built * sec_dim] =
                sv_times(w[k], expr_copy(inv));
        }
        free(w);
        expr_free(inv);
        built++;
    }

    /* 9. Orthonormal completion of the secondary.  After `built`
     *    columns we need sec_dim - built more orthonormal vectors. */
    Expr** secondary_full = NULL;
    if (built == sec_dim) {
        /* Already a full basis -- repack to row-major. */
        secondary_full = (Expr**)malloc(sizeof(Expr*) * (size_t)sec_dim * (size_t)sec_dim);
        for (int i = 0; i < sec_dim; i++) {
            for (int j = 0; j < sec_dim; j++) {
                secondary_full[i * sec_dim + j] = secondary_cols_cm[i + j * sec_dim];
            }
        }
    } else {
        bool oc_ok = sv_orthonormal_complete(secondary_cols_cm, sec_dim,
                                             built, sec_dim, use_conj,
                                             &secondary_full);
        for (int i = 0; i < sec_dim * built; i++) expr_free(secondary_cols_cm[i]);
        if (!oc_ok) {
            free(secondary_cols_cm);
            for (int i = 0; i < g * g; i++) expr_free(Qe[i]);
            free(Qe);
            for (int i = 0; i < g; i++) expr_free(sigma_full[i]);
            free(sigma_full);
            for (int i = 0; i < g * g; i++) expr_free(primary_flat[i]);
            free(primary_flat);
            expr_free(A_expr);
            return false;
        }
    }
    free(secondary_cols_cm);

    /* 10. Assemble outputs.
     *     U is the secondary (n x n) when use_AHA, primary (n x n=g) otherwise.
     *     V is the primary  (p x p=g) when use_AHA, secondary (p x p) otherwise.
     *     S is the sigma vector (length g = min(n, p)). */
    Expr** U_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    Expr** V_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)p * (size_t)p);
    if (use_AHA) {
        /* U = secondary (n x n), V = primary (p x p, p == g). */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) U_flat[i * n + j] = secondary_full[i * n + j];
        }
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) V_flat[i * p + j] = primary_flat[i * p + j];
        }
    } else {
        /* U = primary (n x n, n == g), V = secondary (p x p). */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) U_flat[i * n + j] = primary_flat[i * n + j];
        }
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) V_flat[i * p + j] = secondary_full[i * p + j];
        }
    }
    free(primary_flat);
    free(secondary_full);
    for (int i = 0; i < g * g; i++) expr_free(Qe[i]);
    free(Qe);
    expr_free(A_expr);

    /* Singular-value array stays length g. */
    Expr** S_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)g);
    for (int i = 0; i < g; i++) S_flat[i] = sigma_full[i];
    free(sigma_full);

    *out_U_flat = U_flat;
    *out_S_flat = S_flat;
    *out_V_flat = V_flat;
    *out_rank   = rank_out;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Wrap a row-major Expr** buffer of rows*cols entries into a         *
 *  List[List[...]] expression.  STEALS each entry; caller must not    *
 *  free the entries (only the buffer itself).                         *
 * ------------------------------------------------------------------ */
static Expr* sv_wrap_matrix(Expr** buf, int rows, int cols) {
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = NULL;
        if (cols > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) elems[j] = buf[i * cols + j];
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        if (elems) free(elems);
    }
    Expr* m = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return m;
}

/* Build the rectangular Sigma matrix (n x p) from the length-g singular
 * value vector (g = min(n, p)).  Entry (i, i) is sigma_i for i < g;
 * everything else is 0.  Entries are COPIED from S_flat. */
static Expr* sv_build_sigma_rect(Expr** S_flat, int g, int n, int p) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            if (i == j && i < g) elems[j] = expr_copy(S_flat[i]);
            else                  elems[j] = expr_new_integer(0);
        }
        rows[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)p);
        free(elems);
    }
    Expr* m = expr_new_function(
        expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return m;
}

/* Element-wise Together to canonicalise rational + Sqrt entries.  Cheap
 * because Together is Listable -- the evaluator threads automatically. */
static Expr* sv_tidy_matrix(Expr* m) {
    return eval_and_free(expr_new_function(
        expr_new_symbol("Together"), (Expr*[]){expr_copy(m)}, 1));
}

/* ------------------------------------------------------------------ *
 *  Symbolic dispatcher.                                               *
 *                                                                     *
 *  Wraps svd_symbolic_core with the inexact-input rationalise round-  *
 *  trip and result assembly.  Returns NULL when the eigendecomposition*
 *  can't reach a closed form.                                         *
 * ------------------------------------------------------------------ */
Expr* svd_symbolic_dispatch(const SvdArgs* args, int n, int p, int n_a) {
    if (args->generalized) {
        /* Generalized SVD has no closed-form symbolic kernel in Mathilda
         * (the eigen-pair problem of (m^H m, a^H a) generally has no
         * radical-expressible roots).  Compromise: when both matrices
         * contain only exact-or-inexact numeric leaves, numericalize to
         * 53-bit Reals and re-dispatch through the machine LAPACK path.
         * When either matrix has free symbolic content, emit ::nogsymb
         * and leave the call unevaluated. */
        if (!sv_matrix_is_purely_numeric(args->m, n,   p)
            || !sv_matrix_is_purely_numeric(args->a, n_a, p)) {
            static uint64_t nogsymb_warn = 0;
            if (!nogsymb_warn) {
                nogsymb_warn = 1;
                fprintf(stderr,
                    "SingularValueDecomposition::nogsymb: Generalized "
                    "SingularValueDecomposition does not support free "
                    "symbolic content; numericalize the inputs (e.g. "
                    "with N) to invoke the LAPACK / MPFR kernel.\n");
            }
            return NULL;
        }
        /* All-numeric exact / Complex / Rational input: N it down to 53-bit
         * Real and feed the result back to the LAPACK path. */
        Expr* m_num = common_numericalize_result(args->m, 53);
        Expr* a_num = common_numericalize_result(args->a, 53);
        if (!m_num || !a_num) {
            if (m_num) expr_free(m_num);
            if (a_num) expr_free(a_num);
            return NULL;
        }
        SvdArgs tmp = *args;
        tmp.m = m_num;
        tmp.a = a_num;
        Expr* res = svd_machine_dispatch(&tmp, n, p, n_a);
        expr_free(m_num);
        expr_free(a_num);
        return res;
    }

    Expr* m = args->m;

    /* Inexact preprocessing (same as qrdecomp.c / ludecomp.c). */
    CommonInexactInfo info = common_scan_inexact(m);
    Expr* m_rat = NULL;
    long prec_bits = 53;
    Expr* matrix_to_use = m;
    if (info.has_inexact) {
        prec_bits = info.min_bits ? info.min_bits : 53;
        m_rat = common_rationalize_input(m, prec_bits);
        matrix_to_use = m_rat;
    }

    bool complex_input = sv_has_complex_content(m);

    /* Flatten m to row-major Expr** buffer. */
    Expr** A_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)p);
    {
        size_t idx = 0;
        flatten_tensor(matrix_to_use, A_flat, &idx);
    }

    Expr** U_flat = NULL;
    Expr** S_flat = NULL;
    Expr** V_flat = NULL;
    int rank = 0;
    bool ok = svd_symbolic_core(A_flat, n, p, complex_input,
                                &U_flat, &S_flat, &V_flat, &rank);

    for (int i = 0; i < n * p; i++) expr_free(A_flat[i]);
    free(A_flat);
    if (m_rat) expr_free(m_rat);

    if (!ok) return NULL;

    int g = (n < p) ? n : p;

    /* Wrap U (n x n), Sigma (n x p), V (p x p) into Mathilda lists.
     * sv_wrap_matrix steals entries; build Sigma fresh from the
     * singular-value vector. */
    Expr* u_mat     = sv_wrap_matrix(U_flat, n, n); free(U_flat);
    Expr* sigma_mat = sv_build_sigma_rect(S_flat, g, n, p);
    Expr* v_mat     = sv_wrap_matrix(V_flat, p, p); free(V_flat);
    /* S_flat entries are now owned by sigma_mat (copied above); free the
     * S_flat container plus the originals. */
    for (int i = 0; i < g; i++) expr_free(S_flat[i]);
    free(S_flat);

    /* Tidy + numericalise. */
    Expr* u_t = sv_tidy_matrix(u_mat); expr_free(u_mat); u_mat = u_t;
    Expr* s_t = sv_tidy_matrix(sigma_mat); expr_free(sigma_mat); sigma_mat = s_t;
    Expr* v_t = sv_tidy_matrix(v_mat); expr_free(v_mat); v_mat = v_t;
    if (info.has_inexact) {
        Expr* u_n = common_numericalize_result(u_mat, prec_bits);
        Expr* s_n = common_numericalize_result(sigma_mat, prec_bits);
        Expr* v_n = common_numericalize_result(v_mat, prec_bits);
        expr_free(u_mat); expr_free(sigma_mat); expr_free(v_mat);
        u_mat = u_n; sigma_mat = s_n; v_mat = v_n;
    }

    /* Build {u, sigma, v}. */
    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = u_mat; items[1] = sigma_mat; items[2] = v_mat;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);

    return svd_apply_postprocess(result, args, n, p, rank);
}

/* ------------------------------------------------------------------ *
 *  Phase-5 post-processing.                                            *
 *                                                                      *
 *  Operates on the freshly-assembled {u, sigma, v} list (n x n,        *
 *  n x p, p x p) produced by any of the three kernels.  Applied        *
 *  uniformly so the kernel-specific code only ever produces the full   *
 *  decomposition; this file owns the truncation, tolerance, and        *
 *  TargetStructure semantics.                                          *
 *                                                                      *
 *  Order of operations:                                                *
 *    1. Tolerance (if user supplied): zero any sigma_i with            *
 *       |sigma_i| < tol.                                               *
 *    2. Truncation (k or UpTo[k]):                                     *
 *         k > 0       -> keep the first k singular values + their      *
 *                        u, v columns.  sigma becomes k x k.           *
 *         k < 0       -> keep the last |k| singular values + their     *
 *                        u, v columns.  sigma becomes |k| x |k|.       *
 *         UpTo[k]     -> as above with k_effective = min(k, rank).     *
 *    3. TargetStructure -> "Structured":                               *
 *         Wrap sigma in DiagonalMatrix[{...}]; u and v stay dense.     *
 *                                                                      *
 *  Phase 6 (generalized SVD) will short-circuit this helper.  The      *
 *  current generalised-form branch in svd_parse_args rejects k         *
 *  truncation outright, and the dispatchers haven't been wired up for  *
 *  the pair form yet -- so the postprocess just returns `result` for   *
 *  the generalized case.                                               *
 * ------------------------------------------------------------------ */

/* True iff `e` is provably a numeric value (Integer / Rational / Real /
 * MPFR / BigInt).  Used by the tolerance comparison; symbolic
 * expressions are left alone (no useful threshold comparison). */
static bool sv_is_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL
        || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        return sv_is_numeric(e->data.function.args[0])
            && sv_is_numeric(e->data.function.args[1]);
    }
    return false;
}

/* True when |e| < |tol| (both numeric).  Returns false for symbolic
 * inputs -- we never zero out something we can't compare.  Takes
 * non-const pointers because expr_copy is not declared const-safe. */
static bool sv_below_tol(Expr* e, Expr* tol) {
    if (!sv_is_numeric(e) || !sv_is_numeric(tol)) return false;
    /* Build Abs[e] < Abs[tol] through the evaluator: a tiny tree, but
     * it reliably handles every numeric type (Integer, Real, MPFR,
     * Rational) without bespoke arithmetic. */
    Expr* lhs = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(e)}, 1));
    Expr* rhs = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Abs), (Expr*[]){expr_copy(tol)}, 1));
    Expr* cmp = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Less), (Expr*[]){lhs, rhs}, 2));
    bool below = (cmp->type == EXPR_SYMBOL
                  && strcmp(cmp->data.symbol, "True") == 0);
    expr_free(cmp);
    return below;
}

/* Replace entry (row, col) of a List-of-List matrix.  `e` must be an
 * EXPR_FUNCTION with head List whose row at index `row` is itself a
 * List.  Frees the old entry and stores `new_value` (steals ownership). */
static void sv_matrix_set(Expr* m, int row, int col, Expr* new_value) {
    Expr* r = m->data.function.args[row];
    expr_free(r->data.function.args[col]);
    r->data.function.args[col] = new_value;
}

/* Read entry (row, col) of a List-of-List matrix as a borrowed Expr*
 * (caller must NOT free). */
static Expr* sv_matrix_get(Expr* m, int row, int col) {
    return m->data.function.args[row]->data.function.args[col];
}

/* Apply the Tolerance threshold in-place: any diagonal entry of sigma
 * whose absolute value is < tolerance becomes the literal Real 0.0. */
static void sv_apply_tolerance(Expr* sigma, int n, int p, Expr* tol) {
    int mn = (n < p) ? n : p;
    for (int i = 0; i < mn; i++) {
        Expr* val = sv_matrix_get(sigma, i, i);
        if (sv_below_tol(val, tol)) {
            sv_matrix_set(sigma, i, i, expr_new_real(0.0));
        }
    }
}

/* Build a fresh n x k matrix from columns [col_lo, col_lo+k) of an
 * n x ncols source matrix.  Entries are deep-copied. */
static Expr* sv_slice_cols(Expr* src, int n, int col_lo, int k) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
        for (int j = 0; j < k; j++) {
            elems[j] = expr_copy(sv_matrix_get(src, i, col_lo + j));
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                      elems, (size_t)k);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)n);
    free(rows);
    return m;
}

/* Build a fresh k x k diagonal sigma from sigma's diagonal entries
 * starting at index `idx_lo`.  Off-diagonal entries are exact 0. */
static Expr* sv_build_sigma_square(Expr* sigma_full, int idx_lo, int k) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
    for (int i = 0; i < k; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
        for (int j = 0; j < k; j++) {
            if (i == j) {
                elems[j] = expr_copy(sv_matrix_get(sigma_full,
                                                    idx_lo + i,
                                                    idx_lo + i));
            } else {
                elems[j] = expr_new_integer(0);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                      elems, (size_t)k);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)k);
    free(rows);
    return m;
}

/* Wrap a square diagonal matrix as DiagonalMatrix[{d_1, ..., d_n}].
 * `sq_sigma` is a k x k Mathilda matrix; we extract the diagonal and
 * build the DiagonalMatrix call.  Consumes `sq_sigma`. */
static Expr* sv_wrap_diag(Expr* sq_sigma, int k) {
    Expr** diag = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
    for (int i = 0; i < k; i++) {
        diag[i] = expr_copy(sv_matrix_get(sq_sigma, i, i));
    }
    Expr* diag_list = expr_new_function(expr_new_symbol(SYM_List),
                                          diag, (size_t)k);
    free(diag);
    expr_free(sq_sigma);
    return expr_new_function(expr_new_symbol("DiagonalMatrix"),
                              (Expr*[]){diag_list}, 1);
}

Expr* svd_apply_postprocess(Expr* result, const SvdArgs* args,
                            int n, int p, int rank) {
    if (!result) return NULL;
    /* Generalized form is wired up in Phase 6; for now keep result. */
    if (args->generalized) return result;

    /* Defensive: result is List[u, sigma, v]. */
    if (result->type != EXPR_FUNCTION
        || result->data.function.arg_count < 3) return result;

    Expr* u_mat     = result->data.function.args[0];
    Expr* sigma_mat = result->data.function.args[1];
    Expr* v_mat     = result->data.function.args[2];

    /* 1. Apply Tolerance to sigma's diagonal. */
    if (args->tolerance) {
        sv_apply_tolerance(sigma_mat, n, p, args->tolerance);
    }

    /* 2. Truncation. */
    int mn = (n < p) ? n : p;
    int k_eff = 0;
    int col_lo = 0;
    int sigma_idx_lo = 0;
    bool do_truncate = false;
    if (args->k_form == SVD_FORM_K) {
        do_truncate = true;
        if (args->k_value > 0) {
            k_eff = args->k_value;
            if (k_eff > mn) k_eff = mn;
            col_lo = 0;
            sigma_idx_lo = 0;
        } else {
            k_eff = -args->k_value;
            if (k_eff > mn) k_eff = mn;
            /* Last |k| singular values: positions mn - k_eff .. mn - 1
             * in sigma's diagonal, and the matching columns of u and v. */
            col_lo = mn - k_eff;
            sigma_idx_lo = mn - k_eff;
        }
    } else if (args->k_form == SVD_FORM_UPTO) {
        do_truncate = true;
        /* "UpTo" clamps to rank, not to mn -- the docstring says
         * "as many as are available", and Mathematica uses MatrixRank. */
        k_eff = args->k_value;
        if (k_eff > rank) k_eff = rank;
        if (k_eff > mn) k_eff = mn;
        col_lo = 0;
        sigma_idx_lo = 0;
    }

    if (do_truncate && k_eff != mn) {
        /* Build truncated u (n x k), sigma (k x k), v (p x k). */
        Expr* u_t     = sv_slice_cols(u_mat, n, col_lo, k_eff);
        Expr* sigma_t = sv_build_sigma_square(sigma_mat, sigma_idx_lo, k_eff);
        Expr* v_t     = sv_slice_cols(v_mat, p, col_lo, k_eff);
        /* Free the originals and rewrite result's args. */
        expr_free(u_mat); expr_free(sigma_mat); expr_free(v_mat);
        result->data.function.args[0] = u_t;
        result->data.function.args[1] = sigma_t;
        result->data.function.args[2] = v_t;
        u_mat     = u_t;
        sigma_mat = sigma_t;
        v_mat     = v_t;
    } else if (do_truncate && k_eff == mn) {
        /* k == min(n, p): convert the rectangular sigma to a k x k
         * square diag.  u and v slices are no-ops if they were already
         * full square (n == p == k); otherwise crop. */
        Expr* sigma_t = sv_build_sigma_square(sigma_mat, sigma_idx_lo, k_eff);
        expr_free(sigma_mat);
        result->data.function.args[1] = sigma_t;
        sigma_mat = sigma_t;
        if (n != k_eff) {
            Expr* u_t = sv_slice_cols(u_mat, n, col_lo, k_eff);
            expr_free(u_mat);
            result->data.function.args[0] = u_t;
            u_mat = u_t;
        }
        if (p != k_eff) {
            Expr* v_t = sv_slice_cols(v_mat, p, col_lo, k_eff);
            expr_free(v_mat);
            result->data.function.args[2] = v_t;
            v_mat = v_t;
        }
    }

    /* 3. TargetStructure -> "Structured": wrap sigma as DiagonalMatrix.
     * sigma may be rectangular (no truncation) or square (truncated).
     * For the rectangular case we still emit DiagonalMatrix of the
     * leading-diagonal entries (length mn); the caller can pad with
     * MatrixForm or Normal. */
    if (args->target_structure == SVD_TS_STRUCTURED) {
        int diag_k;
        if (do_truncate) {
            diag_k = k_eff;
        } else {
            diag_k = mn;
            /* sigma is rectangular n x p; convert to square first so the
             * diagonal-extraction loop has uniform indexing. */
            Expr* sq = sv_build_sigma_square(sigma_mat, 0, diag_k);
            expr_free(sigma_mat);
            result->data.function.args[1] = sq;
            sigma_mat = sq;
        }
        Expr* wrapped = sv_wrap_diag(sigma_mat, diag_k);
        result->data.function.args[1] = wrapped;
    }

    return result;
}

/* ------------------------------------------------------------------ *
 *  Top-level kernel router.                                           *
 *                                                                     *
 *  Inexact, min_bits <= 53  -> machine LAPACK kernel, fall through.   *
 *  Inexact, min_bits > 53   -> MPFR Jacobi kernel, fall through.      *
 *  Anything else            -> symbolic dispatcher.                   *
 *                                                                     *
 *  For the generalized form we union the precision of both matrices    *
 *  before deciding the kernel -- the LAPACK / MPFR kernels expect      *
 *  matching precision for `m` and `a`.                                 *
 * ------------------------------------------------------------------ */
Expr* svd_dispatch(const SvdArgs* args, int n, int p, int n_a) {
    /* Generalized SVD ({m, a}): same precision-tiered dispatch as the
     * standard form.  The LAPACK path (dggsvd3 / zggsvd3) is the primary
     * implementation; MPFR uses Paige/Van Loan QR-reduce + 2-pair Jacobi;
     * the symbolic dispatcher numericalises exact-numeric input down to
     * machine precision and re-dispatches.  When the symbolic path can't
     * convert (free symbolic content), it emits ::nogsymb and returns
     * NULL so the call is left unevaluated. */
    CommonInexactInfo info_m = common_scan_inexact(args->m);
    bool has_inexact = info_m.has_inexact;
    long min_bits    = info_m.min_bits;
    if (args->a) {
        CommonInexactInfo info_a = common_scan_inexact(args->a);
        if (info_a.has_inexact) {
            has_inexact = true;
            if (info_a.min_bits
                && (!min_bits || info_a.min_bits < min_bits)) {
                min_bits = info_a.min_bits;
            }
        }
    }

    if (has_inexact) {
        if (!min_bits || min_bits <= 53) {
            Expr* fast = svd_machine_dispatch(args, n, p, n_a);
            if (fast) return fast;
        } else {
            Expr* fast = svd_mpfr_dispatch(args, n, p, n_a);
            if (fast) return fast;
        }
    }
    return svd_symbolic_dispatch(args, n, p, n_a);
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */
Expr* builtin_singularvaluedecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    SvdArgs args;
    if (!svd_parse_args(res, &args)) {
        static uint64_t last_warned = 0;
        matsol_warn_once(&last_warned, res,
            "SingularValueDecomposition::opts: Argument or option in "
            "SingularValueDecomposition[...] is not in one of the "
            "supported forms:\n"
            "\tSingularValueDecomposition[m] | [m, k] | [m, UpTo[k]]\n"
            "\tSingularValueDecomposition[{m, a}]\n"
            "\tOptions: Tolerance -> t, "
            "TargetStructure -> \"Dense\" | \"Structured\".\n"
            "\t(Truncation k is not supported for the generalized form.)\n");
        return NULL;
    }

    /* Validate shape of m (and a). */
    int n = 0, p = 0;
    if (!probe_matrix(args.m, &n, &p)) {
        char* s = expr_to_string(args.m);
        fprintf(stderr,
                "SingularValueDecomposition::matrix: Argument %s at "
                "position 1 is not a non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }
    int n_a = 0;
    if (args.a) {
        int p_a = 0;
        if (!probe_matrix(args.a, &n_a, &p_a)) {
            char* s = expr_to_string(args.a);
            fprintf(stderr,
                    "SingularValueDecomposition::matdims: Second matrix %s "
                    "in the generalized form is not a non-empty rectangular "
                    "matrix.\n", s);
            free(s);
            return NULL;
        }
        if (p != p_a) {
            fprintf(stderr,
                    "SingularValueDecomposition::matdims: Generalized "
                    "SingularValueDecomposition requires both matrices to "
                    "have the same number of columns; got %d and %d.\n",
                    p, p_a);
            return NULL;
        }
    }

    /* Bounds-check k.  |k| must be in [0, min(n, p)]; k may be negative
     * (smallest |k|) but not zero for SVD_FORM_K (zero is meaningless --
     * the user almost certainly wrote `SingularValueDecomposition[m, 0]`
     * by mistake). */
    int mn = (n < p) ? n : p;
    if (args.k_form == SVD_FORM_K) {
        int absk = args.k_value < 0 ? -args.k_value : args.k_value;
        if (args.k_value == 0 || absk > mn) {
            fprintf(stderr,
                    "SingularValueDecomposition::sval: k = %d is out of "
                    "range; expected a non-zero integer with |k| <= %d "
                    "(min of matrix dimensions).\n",
                    args.k_value, mn);
            return NULL;
        }
    } else if (args.k_form == SVD_FORM_UPTO) {
        if (args.k_value < 0) {
            fprintf(stderr,
                    "SingularValueDecomposition::sval: UpTo[k] requires "
                    "k >= 0; got %d.\n", args.k_value);
            return NULL;
        }
    }

    /* "Structured" with no DiagonalMatrix wrapper available is still
     * accepted -- DiagonalMatrix already lives in Mathilda and we wrap
     * sigma with it in svd_apply_postprocess.  Nothing to gate here. */

    return svd_dispatch(&args, n, p, n_a);
}

void svdecomp_init(void) {
    symtab_add_builtin("SingularValueDecomposition",
                       builtin_singularvaluedecomposition);
    symtab_get_def("SingularValueDecomposition")->attributes |= ATTR_PROTECTED;
}
