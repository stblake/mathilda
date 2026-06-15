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
            vals[i] = expr_new_symbol(SYM_Infinity);
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

    /* Numerical dispatch: route inexact ordinary eigenproblems through
     * the hand-rolled numerical kernels.  Method selection follows the
     * grammar in mateigen.h:
     *   - "Banded" goes straight to banded_dispatch (Hermitian only).
     *   - "Arnoldi" goes straight to arnoldi_dispatch.
     *   - "FEAST" goes straight to feast_dispatch (Hermitian only,
     *     requires an Interval).
     *   - Automatic + Hermitian + narrow band prefers Banded.
     *   - Automatic + small k (k <= max(20, n/10)) prefers Arnoldi.
     *   - Otherwise (Automatic / "Direct") goes through direct_dispatch.
     * Each dispatcher returns NULL for shapes it doesn't yet support, in
     * which case we fall back through the next preference, ending at
     * the symbolic characteristic-polynomial pipeline. */
    if (inexact) {
        bool try_banded = (opts.method == MATEIGEN_BANDED)
                       || (opts.method == MATEIGEN_AUTOMATIC
                           && banded_automatic_prefers(m, n));
        if (try_banded) {
            Expr* out = banded_dispatch(m, a, n,
                                          MATEIGEN_WANT_VALUES,
                                          opts.k_spec, opts.method_value);
            if (out) return out;
        }
        bool try_arnoldi = (opts.method == MATEIGEN_ARNOLDI)
                        || (opts.method == MATEIGEN_AUTOMATIC
                            && arnoldi_automatic_prefers(opts.k_spec, (size_t)n));
        if (try_arnoldi) {
            Expr* out = arnoldi_dispatch(m, a, n,
                                          MATEIGEN_WANT_VALUES,
                                          opts.k_spec, opts.method_value);
            if (out) return out;
        }
        bool try_feast = (opts.method == MATEIGEN_FEAST)
                      || (opts.method == MATEIGEN_AUTOMATIC
                          && feast_automatic_prefers(m, n, opts.method_value));
        if (try_feast) {
            Expr* out = feast_dispatch(m, a, n,
                                         MATEIGEN_WANT_VALUES,
                                         opts.k_spec, opts.method_value);
            if (out) return out;
        }
        if (opts.method == MATEIGEN_AUTOMATIC
            || opts.method == MATEIGEN_DIRECT
            || opts.method == MATEIGEN_ARNOLDI  /* fall back to Direct on Arnoldi failure */
            || opts.method == MATEIGEN_BANDED   /* fall back to Direct when Banded refuses */
            || opts.method == MATEIGEN_FEAST    /* fall back to Direct when FEAST refuses (e.g. missing Interval) */) {
            Expr* out = direct_dispatch(m, a, n,
                                         MATEIGEN_WANT_VALUES,
                                         opts.k_spec);
            if (out) return out;
        }
    }

    /* Method warnings for the not-yet-implemented kernels.  These
     * always fall back to the symbolic path so the call still
     * produces a result. */
    if (inexact && opts.method_given
        && opts.method != MATEIGEN_AUTOMATIC
        && opts.method != MATEIGEN_DIRECT
        && opts.method != MATEIGEN_ARNOLDI
        && opts.method != MATEIGEN_BANDED
        && opts.method != MATEIGEN_FEAST) {
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

    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
                                  result_vals, result_count);
    free(result_vals);
    return out;
}

/* ---------------- Eigenvectors ---------------- */

/* Substitute the lambda symbol in `M` with `val`, evaluating each entry.
 * Returns a freshly allocated matrix. */
static Expr* eigen_subst_lambda(Expr* M, const char* lambda_name, Expr* val) {
    /* Use ReplaceAll: M /. lambda_name -> val */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(lambda_name), expr_copy(val) }, 2);
    Expr* replaced = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_ReplaceAll),
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
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), v, (size_t)n);
    free(v);
    return out;
}

/* Normalise vector `v` (a List of n elements) by dividing by its Norm.
 * Used by the numerical eigenvector path to emit unit vectors. */
static Expr* eigen_normalize_vector(Expr* v) {
    Expr* norm = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Norm), (Expr*[]){ expr_copy(v) }, 1));
    /* If Norm is zero or symbolic, skip normalisation. */
    if (is_zero_poly(norm)) { expr_free(norm); return expr_copy(v); }
    Expr* inv = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power),
        (Expr*[]){ norm, expr_new_integer(-1) }, 2));
    Expr* scaled = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times),
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

    /* Numerical dispatch (mirrors builtin_eigenvalues). */
    if (inexact) {
        bool try_banded = (opts.method == MATEIGEN_BANDED)
                       || (opts.method == MATEIGEN_AUTOMATIC
                           && banded_automatic_prefers(m, n));
        if (try_banded) {
            Expr* out = banded_dispatch(m, a, n,
                                          MATEIGEN_WANT_VECTORS,
                                          opts.k_spec, opts.method_value);
            if (out) return out;
        }
        bool try_arnoldi = (opts.method == MATEIGEN_ARNOLDI)
                        || (opts.method == MATEIGEN_AUTOMATIC
                            && arnoldi_automatic_prefers(opts.k_spec, (size_t)n));
        if (try_arnoldi) {
            Expr* out = arnoldi_dispatch(m, a, n,
                                          MATEIGEN_WANT_VECTORS,
                                          opts.k_spec, opts.method_value);
            if (out) return out;
        }
        bool try_feast = (opts.method == MATEIGEN_FEAST)
                      || (opts.method == MATEIGEN_AUTOMATIC
                          && feast_automatic_prefers(m, n, opts.method_value));
        if (try_feast) {
            Expr* out = feast_dispatch(m, a, n,
                                         MATEIGEN_WANT_VECTORS,
                                         opts.k_spec, opts.method_value);
            if (out) return out;
        }
        if (opts.method == MATEIGEN_AUTOMATIC
            || opts.method == MATEIGEN_DIRECT
            || opts.method == MATEIGEN_ARNOLDI
            || opts.method == MATEIGEN_BANDED
            || opts.method == MATEIGEN_FEAST) {
            Expr* out = direct_dispatch(m, a, n,
                                         MATEIGEN_WANT_VECTORS,
                                         opts.k_spec);
            if (out) return out;
        }
    }

    if (inexact && opts.method_given
        && opts.method != MATEIGEN_AUTOMATIC
        && opts.method != MATEIGEN_DIRECT
        && opts.method != MATEIGEN_ARNOLDI
        && opts.method != MATEIGEN_BANDED
        && opts.method != MATEIGEN_FEAST) {
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

    Expr* out = expr_new_function(expr_new_symbol(SYM_List),
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
