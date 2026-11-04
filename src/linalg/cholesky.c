/*
 * src/linalg/cholesky.c -- CholeskyDecomposition[m].
 *
 * Gives the upper-triangular c with ConjugateTranspose[c] . c == m for a
 * Hermitian (or real symmetric) positive-definite matrix m.
 *
 * Dispatch by the precision class of m (la_prec_class):
 *
 *   machine  (NDArray, packed, or boxed <= 53-bit)  -> chol_machine:
 *            LAPACK dpotrf / zpotrf on the upper triangle (an in-house
 *            column Cholesky when LAPACK is absent).  The result is a packed
 *            float64 / complex64 matrix with exact-zero bits below the
 *            diagonal.  This one helper serves the boxed List, the packed
 *            List, the visible NDArray and Compile[] (via ndla_cholesky), so
 *            the four surfaces agree bit for bit.
 *   MPFR     (> 53 bits)                            -> chol_generic at the
 *            input precision;
 *   exact / symbolic                                -> chol_generic over
 *            Exprs, with Sqrt and Conjugate, matching the reference's form:
 *            c[i][i] = Sqrt[m[i][i] - Sum c[k][i] Conjugate[c[k][i]]] and
 *            c[i][j] = (m[i][j] - Sum Conjugate[c[k][i]] c[k][j]) / c[i][i].
 *
 * Messages (reference text):
 *   ::matsq  a non-square or empty matrix;
 *   ::herm   not Hermitian or real symmetric;
 *   ::npdef  an exact matrix that is not positive definite;
 *   ::posdef an inexact matrix that is not (sufficiently) positive definite.
 * A symbolic matrix with a non-positive numeric pivot is left unevaluated
 * without a message, as the reference does.  Option: TargetStructure
 * (Automatic | "Dense" | "Structured" -- all give a dense List).
 *
 * Memory contract: standard builtin ownership (SPEC.md section 4).
 */

#include "la_dense.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "lapack.h"
#include "linalg.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "numeric.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define CHOL_HEAD "CholeskyDecomposition"

/* Outcomes of a kernel, so the caller can emit the matching message. */
typedef enum {
    CHOL_OK = 0,
    CHOL_NOTNUM,     /* not a machine-numeric matrix: try the generic path   */
    CHOL_NOTHERM,    /* not Hermitian / real symmetric                       */
    CHOL_NOTPD,      /* not positive definite (numeric)                      */
    CHOL_DECLINE     /* symbolic, cannot decide: silently unevaluated        */
} CholStatus;

/* ------------------------------------------------------------------------ *
 *  Machine kernel                                                          *
 * ------------------------------------------------------------------------ */

/* Relative Hermitian test on a column-major n x n buffer: every
 * |a_ij - conj(a_ji)| <= 1e-12 * max|a|.  The tolerance is the one the
 * machine Eigenvalues path uses to classify a matrix as symmetric. */
static bool chol_is_hermitian(const double* A, int n, bool cplx) {
    double maxa = 0.0;
    size_t nn = (size_t)n * (size_t)n;
    for (size_t k = 0; k < nn; k++) {
        double v = cplx ? hypot(A[2 * k], A[2 * k + 1]) : fabs(A[k]);
        if (v > maxa) maxa = v;
    }
    double tol = 1e-12 * (maxa > 0.0 ? maxa : 1.0);
    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) {
            size_t ij = (size_t)i + (size_t)j * (size_t)n;
            size_t ji = (size_t)j + (size_t)i * (size_t)n;
            if (!cplx) {
                if (fabs(A[ij] - A[ji]) > tol) return false;
            } else {
                double dr = A[2 * ij] - A[2 * ji];
                double di = A[2 * ij + 1] + A[2 * ji + 1];
                if (hypot(dr, di) > tol) return false;
            }
        }
    }
    return true;
}

/* In-house upper Cholesky (column-major, upper triangle), used only when the
 * LAPACK stubs report "unavailable".  Returns 0 on success or k > 0 when the
 * k-th leading minor is not positive definite (LAPACK's convention). */
static int chol_fallback(double* A, int n, bool cplx) {
    for (int j = 0; j < n; j++) {
        /* diagonal: a_jj - sum_k |u_kj|^2 */
        double d = cplx ? A[2 * ((size_t)j + (size_t)j * n)] : A[(size_t)j + (size_t)j * n];
        for (int k = 0; k < j; k++) {
            size_t kj = (size_t)k + (size_t)j * n;
            d -= cplx ? (A[2 * kj] * A[2 * kj] + A[2 * kj + 1] * A[2 * kj + 1])
                      : A[kj] * A[kj];
        }
        if (!(d > 0.0)) return j + 1;
        double ujj = sqrt(d);
        size_t jj = (size_t)j + (size_t)j * n;
        if (cplx) { A[2 * jj] = ujj; A[2 * jj + 1] = 0.0; } else A[jj] = ujj;
        for (int i = j + 1; i < n; i++) {
            /* u_ji = (a_ji - sum_k conj(u_kj) u_ki) / u_jj */
            size_t ji = (size_t)j + (size_t)i * n;
            if (!cplx) {
                double s = A[ji];
                for (int k = 0; k < j; k++)
                    s -= A[(size_t)k + (size_t)j * n] * A[(size_t)k + (size_t)i * n];
                A[ji] = s / ujj;
            } else {
                double sr = A[2 * ji], si = A[2 * ji + 1];
                for (int k = 0; k < j; k++) {
                    size_t kj = (size_t)k + (size_t)j * n, ki = (size_t)k + (size_t)i * n;
                    double ar = A[2 * kj], ai = -A[2 * kj + 1];      /* conj(u_kj) */
                    double br = A[2 * ki], bi = A[2 * ki + 1];
                    sr -= ar * br - ai * bi;
                    si -= ar * bi + ai * br;
                }
                A[2 * ji] = sr / ujj; A[2 * ji + 1] = si / ujj;
            }
        }
    }
    return 0;
}

/* The machine Cholesky shared by every machine surface.  `m` is an NDArray
 * or a boxed matrix of order n.  On CHOL_OK *out is a packed upper-triangular
 * factor carrying m's presentation. */
static CholStatus chol_machine(const Expr* m, int n, Expr** out) {
    *out = NULL;
    int rows = 0, cols = 0;
    double* A = NULL;
    bool cplx = false;
    if (!na_load_matrix((Expr*)m, false, true, &rows, &cols, &A)) {
        free(A); A = NULL;
        if (!na_load_matrix((Expr*)m, true, true, &rows, &cols, &A)) {
            free(A);
            return CHOL_NOTNUM;
        }
        cplx = true;
    }
    if (rows != n || cols != n) { free(A); return CHOL_NOTNUM; }
    if (!chol_is_hermitian(A, n, cplx)) { free(A); return CHOL_NOTHERM; }

    int info = cplx ? mat_lapack_zpotrf('U', n, A, n) : mat_lapack_dpotrf('U', n, A, n);
    if (info < 0) info = chol_fallback(A, n, cplx);     /* LAPACK unavailable */
    if (info != 0) { free(A); return CHOL_NOTPD; }

    /* Zero the strict lower triangle (dpotrf leaves it untouched). */
    for (int j = 0; j < n; j++)
        for (int i = j + 1; i < n; i++) {
            size_t ij = (size_t)i + (size_t)j * (size_t)n;
            if (cplx) { A[2 * ij] = 0.0; A[2 * ij + 1] = 0.0; } else A[ij] = 0.0;
        }
    *out = na_build_matrix_as(A, n, n, cplx, true, na_result_presentation(m));
    free(A);
    return *out ? CHOL_OK : CHOL_NOTNUM;
}

/* ------------------------------------------------------------------------ *
 *  Generic kernel: exact, symbolic and arbitrary precision                 *
 * ------------------------------------------------------------------------ */

/* Exact-or-MPFR Hermitian test.  Exact/symbolic: m == ConjugateTranspose[m],
 * or m == Transpose[m] with no explicit complex entry (the reference accepts
 * a symbolic symmetric matrix, but not a complex symmetric one).  MPFR: the
 * machine relative test on the doubles. */
static bool chol_generic_hermitian(const Expr* m, int n, LaPrecClass pc) {
    if (pc == LA_MPFR) {
        bool cplx = la_has_complex(m);
        double* A = (double*)calloc((size_t)(cplx ? 2 : 1) * (size_t)n * (size_t)n, sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double re = 0.0, im = 0.0;
                if (!na_read_scalar(la_at(m, i, j), &re, &im)) { free(A); return false; }
                size_t k = (size_t)i + (size_t)j * (size_t)n;
                if (cplx) { A[2 * k] = re; A[2 * k + 1] = im; } else A[k] = re;
            }
        bool h = chol_is_hermitian(A, n, cplx);
        free(A);
        return h;
    }
    Expr* ct = la_ev1(SYM_ConjugateTranspose, expr_copy((Expr*)m));
    bool h = expr_eq(ct, (Expr*)m);
    expr_free(ct);
    if (h) return true;
    if (la_has_complex(m)) return false;
    Expr* tr = la_ev1(SYM_Transpose, expr_copy((Expr*)m));
    h = expr_eq(tr, (Expr*)m);
    expr_free(tr);
    return h;
}

/* Sign test on a pivot: true when it is a number that is not positive (an
 * exact or numerical zero, a negative number, or a non-real number). */
static bool chol_pivot_bad(const Expr* s) {
    if (s->type == EXPR_INTEGER) return s->data.integer <= 0;
    if (!expr_is_numeric_quantity(s)) return false;  /* symbolic: undecidable */
    Expr* ns = la_ev1(SYM_N, expr_copy((Expr*)s));
    double re = 0.0, im = 0.0;
    bool ok = na_read_scalar(ns, &re, &im);
    expr_free(ns);
    if (!ok) return false;
    return !(re > 0.0) || fabs(im) > 1e-12 * fabs(re);
}

/* Generic up-looking Cholesky over Exprs.  `m` is a boxed n x n matrix. */
static CholStatus chol_generic(const Expr* m, int n, LaPrecClass pc, Expr** out) {
    *out = NULL;
    if (!chol_generic_hermitian(m, n, pc)) return CHOL_NOTHERM;

    size_t nn = (size_t)n * (size_t)n;
    Expr** c = (Expr**)calloc(nn, sizeof(Expr*));
    CholStatus st = CHOL_OK;
    for (int i = 0; i < n && st == CHOL_OK; i++) {
        /* pivot s = m_ii - Sum_k c_ki Conjugate[c_ki] */
        Expr* s = expr_copy(la_at(m, i, i));
        for (int k = 0; k < i; k++) {
            Expr* cki = c[(size_t)k * n + i];
            Expr* t = la_ev2(SYM_Times, expr_copy(cki),
                             la_call1(SYM_Conjugate, expr_copy(cki)));
            s = la_ev2(SYM_Plus, s, la_ev2(SYM_Times, expr_new_integer(-1), t));
        }
        s = eval_and_free(s);
        if (chol_pivot_bad(s)) {
            expr_free(s);
            st = (pc == LA_EXACT && !la_all_numeric(m)) ? CHOL_DECLINE : CHOL_NOTPD;
            break;
        }
        Expr* cii = la_ev1(SYM_Sqrt, s);
        c[(size_t)i * n + i] = cii;
        for (int j = i + 1; j < n; j++) {
            Expr* num = expr_copy(la_at(m, i, j));
            for (int k = 0; k < i; k++) {
                Expr* t = la_ev2(SYM_Times,
                                 la_call1(SYM_Conjugate, expr_copy(c[(size_t)k * n + i])),
                                 expr_copy(c[(size_t)k * n + j]));
                num = la_ev2(SYM_Plus, num, la_ev2(SYM_Times, expr_new_integer(-1), t));
            }
            c[(size_t)i * n + j] = la_ev2(SYM_Times, num,
                                          la_call2(SYM_Power, expr_copy(cii),
                                                   expr_new_integer(-1)));
        }
        for (int j = 0; j < i; j++) c[(size_t)i * n + j] = expr_new_integer(0);
    }
    if (st != CHOL_OK) {
        for (size_t k = 0; k < nn; k++) if (c[k]) expr_free(c[k]);
        free(c);
        return st;
    }
    *out = la_matrix_from(c, n, n);
    free(c);
    return CHOL_OK;
}

/* ------------------------------------------------------------------------ *
 *  Messages and entry points                                               *
 * ------------------------------------------------------------------------ */

static void chol_report(CholStatus st, const Expr* m) {
    if (st == CHOL_OK || st == CHOL_DECLINE || st == CHOL_NOTNUM) return;
    char* s = la_str(m);
    if (st == CHOL_NOTHERM)
        mth_message(CHOL_HEAD, "herm", "The matrix %s is not Hermitian or real and symmetric.", s);
    else if (la_prec_class(m, NULL) == LA_EXACT)
        mth_message(CHOL_HEAD, "npdef", "The matrix %s is not positive definite.", s);
    else
        mth_message(CHOL_HEAD, "posdef", "The matrix %s is not sufficiently positive "
                    "definite to complete the Cholesky decomposition to reasonable "
                    "accuracy.", s);
    free(s);
}

/* Accept the TargetStructure option (any of its values gives a dense List);
 * report anything else as the reference does. */
static bool chol_check_options(const Expr* res) {
    size_t argc = res->data.function.arg_count;
    for (size_t i = 1; i < argc; i++) {
        const Expr* o = res->data.function.args[i];
        bool rule = o->type == EXPR_FUNCTION && o->data.function.arg_count == 2
                 && (head_is(o, SYM_Rule) || head_is(o, SYM_RuleDelayed));
        if (!rule) {
            builtin_arg_error(CHOL_HEAD, argc, 1, 1);
            return false;
        }
        const Expr* key = o->data.function.args[0];
        if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_TargetStructure)
            continue;
        char* k = la_str(key);
        char* call = la_str(res);
        mth_message(CHOL_HEAD, "optx", "Unknown option %s in %s.", k, call);
        free(k); free(call);
        return false;
    }
    return true;
}

/* Order of a square non-empty matrix argument, or -1 (after the ::matsq
 * message for a matrix-shaped argument). */
static int chol_order(const Expr* m) {
    int r = 0, c = 0;
    if (la_matrix_dims(m, &r, &c) && r > 0 && r == c) return r;
    if (la_looks_like_matrix(m)) {
        char* s = la_str(m);
        mth_message(CHOL_HEAD, "matsq", "Argument %s at position 1 is not a non-empty square matrix.", s);
        free(s);
    }
    return -1;
}

/* Compile[]'s direct kernel: the machine Cholesky of an NDArray, or NULL --
 * silently, so a compiled degrade does not print the message the interpreter
 * fallback is about to print. */
Expr* la_cholesky_machine(const Expr* m) {
    int r = 0, c = 0;
    if (!la_matrix_dims(m, &r, &c) || r < 1 || r != c) return NULL;
    Expr* out = NULL;
    if (chol_machine(m, r, &out) != CHOL_OK) { if (out) expr_free(out); return NULL; }
    return out;
}

/* The NDArray entry point: CholeskyDecomposition[nd] on a machine
 * buffer.  Never delists and re-evaluates: everything it cannot answer is
 * reported (matsq/herm/posdef) and left unevaluated, which is exactly what
 * the boxed path does for the same matrix. */
Expr* ndla_choleskydecomposition(Expr* res) {
    if (res->data.function.arg_count < 1 || !chol_check_options(res)) return NULL;
    const Expr* m = res->data.function.args[0];
    if (!is_ndarray(m)) return NULL;
    int n = chol_order(m);
    if (n < 1) return NULL;
    Expr* out = NULL;
    CholStatus st = chol_machine(m, n, &out);
    if (st == CHOL_NOTNUM) return linalg_delist_and_reeval(res);
    chol_report(st, m);
    return out;
}

Expr* builtin_choleskydecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return builtin_arg_error(CHOL_HEAD, argc, 1, 1);
    if (linalg_call_has_ndarray(res)) return ndla_choleskydecomposition(res);
    if (!chol_check_options(res)) return NULL;

    const Expr* m = res->data.function.args[0];
    int n = chol_order(m);
    if (n < 1) return NULL;

    long bits = 53;
    LaPrecClass pc = la_prec_class(m, &bits);
    Expr* out = NULL;
    CholStatus st = CHOL_NOTNUM;
    if (pc == LA_MACHINE) st = chol_machine(m, n, &out);
    if (st == CHOL_NOTNUM) {
        /* MPFR, exact, symbolic -- or a machine matrix with a symbolic leaf. */
        Expr* boxed = la_boxed_copy(m);
        st = chol_generic(boxed, n, pc, &out);
        expr_free(boxed);
    }
    chol_report(st, m);
    return out;
}

void cholesky_init(void) {
    symtab_add_builtin(CHOL_HEAD, builtin_choleskydecomposition);
    symtab_get_def(CHOL_HEAD)->attributes |= ATTR_PROTECTED;
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ la_call2(SYM_Rule, expr_new_symbol(SYM_TargetStructure),
                            expr_new_symbol(SYM_Automatic)) }, 1);
    symtab_set_options(CHOL_HEAD, opts);
}
