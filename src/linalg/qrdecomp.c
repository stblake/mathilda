/* qrdecomp.c
 *
 * QRDecomposition[m]                   -- {q, r} thin-QR factorisation.
 * QRDecomposition[m, Pivoting -> True] -- {q, r, p}, m . p == q^H . r.
 *
 * Strategy.  One algorithmic core - Modified Gram-Schmidt on the
 * columns of m, driven through the Mathilda evaluator - serves
 * every input family:
 *
 *   - exact integer / rational / complex / free-symbolic matrices
 *     run the pipeline as-is.  The output is exact (Sqrt[...] in
 *     the norms, Rational / symbolic entries elsewhere).
 *
 *   - inexact matrices (Real or MPFR leaves) follow the
 *     rationalise -> exact pipeline -> numericalise round-trip used
 *     by PseudoInverse / Eigenvalues / Solve.  The output precision
 *     matches the minimum precision present in the input, mirroring
 *     the inexact-in / inexact-out contract advertised across the
 *     rest of the system.
 *
 * Conventions.  We work internally with a standard "thin" QR
 *
 *     A = Q . R          Q n x r orthonormal-columns,  R r x p upper
 *                        trapezoidal,  r = MatrixRank[A]
 *
 * and at the end return  q = ConjugateTranspose[Q],  r = R  so the
 * Mathematica identity  m == ConjugateTranspose[q] . r  holds.
 * Because ConjugateTranspose is involutive this matches the spec
 * convention exactly: Length[q] == Length[r] == r (the rank), the
 * rows of q are orthonormal in the complex inner product, and r has
 * zeros below the leading-diagonal echelon.
 *
 * Modified Gram-Schmidt loop, column k = 0 .. p-1:
 *
 *     v = A[:, k]
 *     for each existing orthonormal column Q[:, j]:
 *         coeff = <Q[:, j], v> = Sum_i Conjugate[Q[i, j]] * v[i]
 *         R[j, k] = coeff
 *         v -= coeff * Q[:, j]
 *     norm_sq = <v, v>
 *     if norm_sq == 0: column is dependent, skip (no new orthonormal row)
 *     norm = Sqrt[norm_sq]
 *     R[rk, k] = norm
 *     Q[:, rk] = v / norm
 *     rk += 1
 *
 * After the loop q is built as  q[j, i] = Conjugate[Q[i, j]] - this
 * collapses to Q^T for real matrices and gives the proper
 * conjugate-transpose for complex matrices.
 *
 * Pivoting (when Pivoting -> True).  At the start of each step we
 * pick, among the remaining columns of A, the one whose residual
 * orthogonal projection (after subtracting components along the
 * already-built Q columns) has the largest squared norm.  This is
 * exactly Householder column pivoting expressed in MGS form and
 * makes the diagonal of R appear in order of decreasing magnitude,
 * matching the Mathematica example.  The permutation array is then
 * inflated into a p x p permutation matrix p with
 *
 *     p[perm[j], j] = 1
 *
 * so that A . p picks the columns in the chosen order, satisfying
 * A . p == ConjugateTranspose[q] . r.
 *
 * Memory contract.  Standard builtin contract.  This file does NOT
 * call expr_free(res) - the evaluator owns `res` and frees it on a
 * non-NULL return (MEMORY.md / SPEC.md §4.1).  Every intermediate
 * allocation is tracked: the Q and R working buffers are freed after
 * the q / r List wrappers have stolen / copied their entries.
 */

#include "qrdecomp.h"
#include "qrdecomp_internal.h"
#include "linalg.h"
#include "linsolve.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"
#include "common.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Small evaluator helpers.  Each one builds a tiny expression tree   *
 *  for a primitive operation, evaluates it, and returns ownership of  *
 *  the result to the caller.  The argument trees are consumed.        *
 * ------------------------------------------------------------------ */
static Expr* eval_plus(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Plus), (Expr*[]){a, b}, 2));
}
static Expr* eval_times(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times), (Expr*[]){a, b}, 2));
}
static Expr* eval_power(Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power), (Expr*[]){a, b}, 2));
}
static Expr* eval_conjugate(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Conjugate), (Expr*[]){a}, 1));
}
static Expr* eval_sqrt(Expr* a) {
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Sqrt), (Expr*[]){a}, 1));
}
static Expr* eval_simplify(Expr* a) {
    /* Together is enough for canonicalising zero-detection on
     * rational + Sqrt expressions, and is materially cheaper than
     * Simplify.  is_zero_poly works correctly on the result. */
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Together), (Expr*[]){a}, 1));
}

/* ------------------------------------------------------------------ *
 *  Vector arithmetic.  Each helper allocates a fresh length-n array   *
 *  and returns it; the caller owns every entry.                       *
 * ------------------------------------------------------------------ */

/* dst -= coeff * src.  Mutates dst in place, consumes the old dst[i]. */
static void vec_axpy_negate(Expr** dst, Expr* coeff, Expr** src, int n) {
    for (int i = 0; i < n; i++) {
        /* prod = coeff * src[i]  (both copied; coeff is reused next iter) */
        Expr* prod = eval_times(expr_copy(coeff), expr_copy(src[i]));
        Expr* neg  = eval_times(expr_new_integer(-1), prod);
        Expr* new_v = eval_plus(dst[i], neg);   /* dst[i] consumed */
        dst[i] = new_v;
    }
}

/* <a, b>.  Two paths:
 *
 *   - Complex inner product (use_conj=true):
 *         Sum_i Conjugate[a[i]] * b[i]
 *     The standard Hermitian inner product on C^n.
 *
 *   - Real inner product  (use_conj=false):
 *         Sum_i a[i] * b[i]
 *     Used for matrices with no Complex head and no I symbol
 *     anywhere - including symbolic free variables, which Mathematica
 *     treats as potentially complex but Mathilda elects to handle
 *     under a real-valued assumption.  This avoids spurious
 *     Conjugate[a] residues in the printed q and r entries.
 *
 * Returns a freshly allocated Expr*. */
static Expr* inner_product(Expr** a, Expr** b, int n, bool use_conj) {
    Expr* sum = expr_new_integer(0);
    for (int i = 0; i < n; i++) {
        Expr* lhs  = use_conj ? eval_conjugate(expr_copy(a[i]))
                              : expr_copy(a[i]);
        Expr* term = eval_times(lhs, expr_copy(b[i]));
        sum        = eval_plus(sum, term);
    }
    return sum;
}

/* <v, v> for the convergence / rank test.  Same as inner_product
 * with a == b. */
static Expr* norm_squared(Expr** v, int n, bool use_conj) {
    return inner_product(v, v, n, use_conj);
}

/* dst[i] = src[i] / scalar.  Allocates a fresh array; consumes no
 * input ownership.  The src entries are NOT consumed - the caller
 * still owns them. */
static Expr** vec_div_scalar(Expr** src, Expr* scalar, int n) {
    Expr* inv = eval_power(expr_copy(scalar), expr_new_integer(-1));
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        out[i] = eval_times(expr_copy(src[i]), expr_copy(inv));
    }
    expr_free(inv);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Squared-norm zero test.                                            *
 *                                                                     *
 *  is_zero_poly is robust on integer / rational / polynomial          *
 *  expressions but is conservatively "false" on Sqrt-bearing forms.   *
 *  We therefore canonicalise via Together first and then test.        *
 *  Returns true only when the expression is provably zero;            *
 *  ambiguous symbolic forms keep building the QR (which is correct -  *
 *  the generic-rank assumption matches Mathematica's behaviour).      *
 * ------------------------------------------------------------------ */
static bool is_definitely_zero(Expr* e) {
    Expr* simp = eval_simplify(expr_copy(e));
    bool z = is_zero_poly(simp);
    expr_free(simp);
    return z;
}

/* ------------------------------------------------------------------ *
 *  Option parsing.                                                    *
 *                                                                     *
 *  The QrOpts struct and its enum live in qrdecomp_internal.h so the  *
 *  Phase 3 / 4 kernel dispatchers can share the parsed form without   *
 *  re-parsing the option list.  The exported entry point is           *
 *  qr_parse_options.                                                  *
 * ------------------------------------------------------------------ */
static bool parse_bool_value(Expr* rhs, bool* out) {
    if (rhs->type != EXPR_SYMBOL) return false;
    const char* s = rhs->data.symbol;
    if (strcmp(s, "True")  == 0) { *out = true;  return true; }
    if (strcmp(s, "False") == 0) { *out = false; return true; }
    return false;
}

static bool parse_targetstructure_value(Expr* rhs, QrTargetStructure* out) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
        *out = QR_TS_DENSE;
        return true;
    }
    if (rhs->type != EXPR_STRING) return false;
    if (strcmp(rhs->data.string, "Dense")      == 0) { *out = QR_TS_DENSE;      return true; }
    if (strcmp(rhs->data.string, "Structured") == 0) { *out = QR_TS_STRUCTURED; return true; }
    return false;
}

bool qr_parse_options(Expr* res, QrOpts* opts) {
    opts->pivoting         = false;
    opts->target_structure = QR_TS_DENSE;

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

        if (lhs->data.symbol == SYM_Pivoting) {
            bool v;
            if (!parse_bool_value(rhs, &v)) return false;
            opts->pivoting = v;
        } else if (lhs->data.symbol == SYM_TargetStructure) {
            QrTargetStructure ts;
            if (!parse_targetstructure_value(rhs, &ts)) return false;
            opts->target_structure = ts;
        } else {
            return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Modified Gram-Schmidt core (Phase 2 rename of mgs_core).           *
 *                                                                     *
 *  Header decl is in qrdecomp_internal.h so Phase 3 / 4 kernel        *
 *  dispatchers can fall back into the symbolic path when LAPACK / MPFR *
 *  is unavailable, without re-implementing the inner loop.            *
 * ------------------------------------------------------------------ */
bool qr_symbolic_core(Expr** A_flat, int n, int p,
                      bool with_pivoting, int* perm,
                      bool use_conj,
                      Expr*** out_Q_flat, Expr*** out_R_flat, int* out_rank) {
    /* Q has at most min(n, p) columns -- the maximum possible rank. */
    int rank_cap = (n < p) ? n : p;

    /* Working column slots for Q: n rows x rank_cap columns. */
    Expr** Q = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)rank_cap);
    /* Working R: rank_cap rows x p cols, initialised to zero. */
    Expr** R = (Expr**)malloc(sizeof(Expr*) * (size_t)rank_cap * (size_t)p);
    for (int t = 0; t < rank_cap * p; t++) R[t] = expr_new_integer(0);

    /* `pending[k]` is the current residual of original column k after
     * projecting away the orthonormal vectors already built.  We keep
     * it across iterations so the pivot-selection step doesn't have
     * to re-orthogonalise every remaining column from scratch.  Each
     * entry is a length-n vector of freshly-allocated Expr*. */
    Expr*** pending = (Expr***)malloc(sizeof(Expr**) * (size_t)p);
    for (int k = 0; k < p; k++) {
        pending[k] = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) pending[k][i] = expr_copy(A_flat[i * p + k]);
    }
    bool* done = (bool*)calloc((size_t)p, sizeof(bool));

    int rk = 0;
    for (int step = 0; step < p; step++) {
        /* Pick the next column.  Without pivoting we just take the
         * leftmost not-yet-handled column; with pivoting we take the
         * one whose residual squared norm is the largest. */
        int chosen = -1;
        Expr* chosen_norm_sq = NULL;
        if (with_pivoting) {
            for (int k = 0; k < p; k++) {
                if (done[k]) continue;
                Expr* nsq = norm_squared(pending[k], n, use_conj);
                if (chosen < 0) {
                    chosen = k;
                    chosen_norm_sq = nsq;
                    continue;
                }
                /* Compare nsq against chosen_norm_sq via expr_compare;
                 * a positive result means nsq is "greater" in canonical
                 * order.  Numeric atoms compare correctly under
                 * expr_compare; for symbolic norms this still produces
                 * a deterministic (if not necessarily magnitude-faithful)
                 * pivot order, matching the non-pivoting fallback. */
                if (expr_compare(nsq, chosen_norm_sq) > 0) {
                    expr_free(chosen_norm_sq);
                    chosen = k;
                    chosen_norm_sq = nsq;
                } else {
                    expr_free(nsq);
                }
            }
        } else {
            for (int k = 0; k < p; k++) {
                if (done[k]) continue;
                chosen = k;
                break;
            }
        }
        if (chosen < 0) break;
        if (perm) perm[step] = chosen;
        done[chosen] = true;

        /* R coefficient indices.  When pivoting is off, `step ==
         * chosen` and R is indexed by original column.  When pivoting
         * is on, R is indexed by *pivot step* (the order in which
         * columns were chosen) - this is the column ordering of the
         * permuted matrix m.p, and the spec identity
         *     m . p  ==  ConjugateTranspose[q] . r
         * requires r's columns to follow m.p's columns, not m's. */
        int rcol = with_pivoting ? step : chosen;

        /* R[j, rcol] = <Q[:, j], A[:, chosen]> (exact-arithmetic
         * identity: m's chosen column equals sum_j R[j, rcol] Q[:, j]).
         * We compute against the *original* A[:, chosen] not the
         * residual so the equality holds before any rounding. */
        for (int j = 0; j < rk; j++) {
            Expr** qcol = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) qcol[i] = Q[i * rank_cap + j];  /* borrowed */
            Expr** ak_orig = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) ak_orig[i] = A_flat[i * p + chosen];  /* borrowed */
            Expr* coeff = inner_product(qcol, ak_orig, n, use_conj);
            free(qcol);
            free(ak_orig);
            expr_free(R[j * p + rcol]);
            R[j * p + rcol] = coeff;
        }

        /* Compute the post-projection norm of pending[chosen] (this is
         * the residual we use to decide rank deficiency). */
        Expr* nsq;
        if (chosen_norm_sq) { nsq = chosen_norm_sq; chosen_norm_sq = NULL; }
        else                  nsq = norm_squared(pending[chosen], n, use_conj);

        if (is_definitely_zero(nsq)) {
            /* Column is in the span of the existing Q's.  R's chosen
             * column already has its coefficients; no new row added. */
            expr_free(nsq);
            for (int i = 0; i < n; i++) expr_free(pending[chosen][i]);
            free(pending[chosen]);
            pending[chosen] = NULL;
            continue;
        }

        Expr* norm = eval_sqrt(nsq);          /* nsq consumed */

        /* R[rk, rcol] = norm.  Same indexing rule as above: rcol is
         * the pivot step when pivoting is on, otherwise the original
         * column index. */
        expr_free(R[rk * p + rcol]);
        R[rk * p + rcol] = expr_copy(norm);

        /* Q[:, rk] = pending[chosen] / norm. */
        Expr** new_col = vec_div_scalar(pending[chosen], norm, n);
        for (int i = 0; i < n; i++) Q[i * rank_cap + rk] = new_col[i];
        free(new_col);

        /* Update every still-active pending column by subtracting its
         * component along the new Q[:, rk]. */
        for (int k = 0; k < p; k++) {
            if (done[k]) continue;
            Expr** qcol = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (int i = 0; i < n; i++) qcol[i] = Q[i * rank_cap + rk];   /* borrowed */
            Expr* coeff = inner_product(qcol, pending[k], n, use_conj);
            vec_axpy_negate(pending[k], coeff, qcol, n);
            free(qcol);
            expr_free(coeff);
        }

        /* pending[chosen] is no longer needed - we have its Q column. */
        for (int i = 0; i < n; i++) expr_free(pending[chosen][i]);
        free(pending[chosen]);
        pending[chosen] = NULL;

        expr_free(norm);
        rk++;
    }

    /* Tidy pending. */
    for (int k = 0; k < p; k++) {
        if (pending[k]) {
            for (int i = 0; i < n; i++) expr_free(pending[k][i]);
            free(pending[k]);
        }
    }
    free(pending);
    free(done);

    /* Compact Q to its first rk columns (still column-major n x rank_cap;
     * we just expose rk to the caller).  Allocate a tight n*rk buffer
     * so the caller doesn't have to know rank_cap. */
    Expr** Q_tight = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)rk);
    for (int j = 0; j < rk; j++) {
        for (int i = 0; i < n; i++) {
            Q_tight[i * rk + j] = Q[i * rank_cap + j];   /* steal */
        }
    }
    /* Free the leftover unused Q slots (j >= rk). */
    for (int j = rk; j < rank_cap; j++) {
        for (int i = 0; i < n; i++) expr_free(Q[i * rank_cap + j]);
    }
    free(Q);

    /* Tight R: rk x p row-major. */
    Expr** R_tight = (Expr**)malloc(sizeof(Expr*) * (size_t)rk * (size_t)p);
    for (int j = 0; j < rk; j++) {
        for (int k = 0; k < p; k++) {
            R_tight[j * p + k] = R[j * p + k];           /* steal */
        }
    }
    for (int j = rk; j < rank_cap; j++) {
        for (int k = 0; k < p; k++) expr_free(R[j * p + k]);
    }
    free(R);

    *out_Q_flat = Q_tight;
    *out_R_flat = R_tight;
    *out_rank   = rk;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Wrap a flat row-major buffer of length rows*cols into a            *
 *  List[List[...], ...] expression.  Steals each entry: the caller   *
 *  must not free buf's entries (only the buf array itself).           *
 * ------------------------------------------------------------------ */
static Expr* wrap_matrix(Expr** buf, int rows, int cols) {
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = NULL;
        if (cols > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) elems[j] = buf[i * cols + j];   /* steal */
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        if (elems) free(elems);
    }
    Expr* out = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return out;
}

/* Build the q matrix (rank x n) from a column-major Q buffer
 * (n x rank).  Two paths:
 *
 *   - Complex input:  q[j, i] = Conjugate[Q[i, j]]  -- the algebraic
 *     definition that makes ConjugateTranspose[q] . r == m hold.
 *
 *   - Real input (no Complex head, no I symbol anywhere in the
 *     original matrix): q[j, i] = Q[i, j].  ConjugateTranspose
 *     collapses to plain Transpose for real entries, so this satisfies
 *     the identity without ever introducing Conjugate[...] heads -
 *     a critical optimisation because Mathilda's Conjugate does not
 *     simplify Conjugate[Conjugate[x]] for symbolic real x, which
 *     would otherwise pollute the printed form. */
static Expr* build_q_from_Q(Expr** Q_flat, int n, int rank,
                             bool complex_input) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)rank);
    for (int j = 0; j < rank; j++) {
        Expr** elems = NULL;
        if (n > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            Expr* entry = expr_copy(Q_flat[i * rank + j]);
            elems[i] = complex_input ? eval_conjugate(entry) : entry;
        }
        rows[j] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)n);
        if (elems) free(elems);
    }
    Expr* q = expr_new_function(
        expr_new_symbol(SYM_List), rows, (size_t)rank);
    free(rows);
    return q;
}

/* Walk a matrix expression for any `Complex[...]` head or the literal
 * symbol `I`.  Used to decide whether QRDecomposition should apply
 * Conjugate when assembling q from the orthonormal column buffer Q. */
static bool has_complex_content(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol, "I") == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex) return true;
    if (has_complex_content(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_complex_content(e->data.function.args[i])) return true;
    }
    return false;
}

/* Build a p x p permutation matrix from a length-p permutation array
 * `perm`, where column k of m.p is column perm[k] of m.  This is
 * equivalent to  P[i, j] = 1 iff i == perm[j], 0 otherwise. */
static Expr* build_perm_matrix(const int* perm, int p) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
    for (int i = 0; i < p; i++) {
        Expr** elems = NULL;
        if (p > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            elems[j] = expr_new_integer(perm[j] == i ? 1 : 0);
        }
        rows[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)p);
        if (elems) free(elems);
    }
    Expr* P = expr_new_function(
        expr_new_symbol(SYM_List), rows, (size_t)p);
    free(rows);
    return P;
}

/* ------------------------------------------------------------------ *
 *  Apply Expand element-wise to canonicalise small symbolic           *
 *  cancellations in q and r.  Cheap and keeps the printed form tidy.  *
 * ------------------------------------------------------------------ */
static Expr* tidy_matrix(Expr* m) {
    /* Expand is Listable -> threads automatically. */
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Expand), (Expr*[]){expr_copy(m)}, 1));
}

/* ------------------------------------------------------------------ *
 *  Symbolic kernel dispatcher.                                        *
 *                                                                     *
 *  Wraps qr_symbolic_core with                                        *
 *    1. inexact-input rationalisation (Real / MPFR leaves promoted to *
 *       exact rationals at min(input bits)),                          *
 *    2. flatten -> work -> wrap conversion between List-of-List and   *
 *       the row-major Expr** buffer the inner loop expects,           *
 *    3. tidy_matrix (element-wise Expand) on q and r,                 *
 *    4. numericalisation back to input precision when rationalisation *
 *       ran, and                                                      *
 *    5. assembly of {q, r} or {q, r, p}.                              *
 *                                                                     *
 *  Phase 3 / 4 sibling dispatchers (qr_machine_dispatch,              *
 *  qr_mpfr_dispatch) will share this entry point as their fallback    *
 *  when LAPACK / MPFR is unavailable.                                 *
 * ------------------------------------------------------------------ */
Expr* qr_symbolic_dispatch(Expr* m, int n, int p, const QrOpts* opts) {
    /* Inexact-preprocessing pipeline (same shape as PseudoInverse).
     * Real / MPFR leaves are rationalised at the minimum input
     * precision; the exact pipeline runs; we numericalise the result
     * back to that precision so the output is Real / MPFR at the
     * input precision.  Pure-exact inputs (no Real / MPFR leaves) run
     * unmodified. */
    CommonInexactInfo info = common_scan_inexact(m);
    Expr* m_rat = NULL;
    long prec_bits = 53;
    Expr* matrix_to_use = m;
    if (info.has_inexact) {
        prec_bits = info.min_bits ? info.min_bits : 53;
        m_rat = common_rationalize_input(m, prec_bits);
        matrix_to_use = m_rat;
    }

    /* Flatten m (row-major) into Expr** for the MGS worker. */
    Expr** A_flat = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)p);
    {
        size_t idx = 0;
        flatten_tensor(matrix_to_use, A_flat, &idx);
    }

    int* perm = NULL;
    if (opts->pivoting) {
        perm = (int*)malloc(sizeof(int) * (size_t)p);
        for (int k = 0; k < p; k++) perm[k] = -1;
    }

    /* Complex content drives the Hermitian-vs-real arithmetic choice
     * and the q-vs-conj(q) assembly choice.  Matrices we rationalised
     * already had any I leaves promoted away (rationalisation only
     * touches Real / MPFR leaves), so we probe the *original* m. */
    bool complex_input = has_complex_content(m);

    Expr** Q_flat = NULL;
    Expr** R_flat = NULL;
    int    rank   = 0;
    bool ok = qr_symbolic_core(A_flat, n, p, opts->pivoting, perm,
                               complex_input,
                               &Q_flat, &R_flat, &rank);

    for (int t = 0; t < n * p; t++) expr_free(A_flat[t]);
    free(A_flat);
    if (m_rat) expr_free(m_rat);

    if (!ok) {
        if (perm) free(perm);
        return NULL;
    }

    /* Wrap q (rank x n) and r (rank x p) into List-of-List expressions.
     * For the zero-rank case (all-zero matrix) we still produce well-
     * formed empty matrices: {} for q and {} for r. */
    Expr* q;
    Expr* r;
    if (rank == 0) {
        /* Free the (empty) buffers and emit `{}` / `{}`. */
        free(Q_flat);
        free(R_flat);
        q = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        r = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else {
        q = build_q_from_Q(Q_flat, n, rank, complex_input);
        for (int t = 0; t < n * rank; t++) expr_free(Q_flat[t]);
        free(Q_flat);

        r = wrap_matrix(R_flat, rank, p);          /* steals entries */
        free(R_flat);
    }

    /* Tidy symbolic entries and numericalise if we rationalised. */
    Expr* q_t = tidy_matrix(q); expr_free(q); q = q_t;
    Expr* r_t = tidy_matrix(r); expr_free(r); r = r_t;
    if (info.has_inexact) {
        Expr* q_n = common_numericalize_result(q, prec_bits);
        Expr* r_n = common_numericalize_result(r, prec_bits);
        expr_free(q); expr_free(r);
        q = q_n; r = r_n;
    }

    /* Build the result list. */
    Expr* result;
    if (opts->pivoting) {
        Expr* P = build_perm_matrix(perm, p);
        free(perm);
        /* Permutation matrix entries are exact 0 / 1; no
         * numericalisation needed even when info.has_inexact. */
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
        items[0] = q; items[1] = r; items[2] = P;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
        free(items);
    } else {
        if (perm) free(perm);
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 2);
        items[0] = q; items[1] = r;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 2);
        free(items);
    }
    return result;
}

/* ------------------------------------------------------------------ *
 *  Top-level kernel router.                                           *
 *                                                                     *
 *  Phase 3 enabled the LAPACK fast path: any inexact input whose      *
 *  minimum leaf precision is <= 53 bits (i.e. dominated by IEEE        *
 *  doubles) is sent to qr_machine_dispatch first.                     *
 *                                                                     *
 *  Phase 4 adds the MPFR Householder fast path: any inexact input     *
 *  whose minimum leaf precision is > 53 bits is sent to               *
 *  qr_mpfr_dispatch first.                                            *
 *                                                                     *
 *  On any soft failure (a stray non-numeric leaf, LAPACK info > 0,    *
 *  USE_LAPACK=0, USE_MPFR=0, rank-deficient + no-pivot) the chosen    *
 *  kernel returns NULL and we fall through to the symbolic            *
 *  dispatcher -- which still understands the inexact input via its    *
 *  rationalise / numericalise round-trip.  The public contract --     *
 *  "QRDecomposition works on any matrix" -- is preserved regardless   *
 *  of how the host BLAS/LAPACK/MPFR is configured.                    */
Expr* qr_dispatch(Expr* m, int n, int p, const QrOpts* opts) {
    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact) {
        if (info.min_bits <= 53) {
            Expr* fast = qr_machine_dispatch(m, n, p, opts);
            if (fast) return fast;
        } else {
            Expr* fast = qr_mpfr_dispatch(m, n, p, opts);
            if (fast) return fast;
        }
    }
    return qr_symbolic_dispatch(m, n, p, opts);
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */
Expr* builtin_qrdecomposition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    QrOpts opts;
    if (!qr_parse_options(res, &opts)) {
        static uint64_t last_warned = 0;
        matsol_warn_once(&last_warned, res,
            "QRDecomposition::opts: Options are not one of "
            "Pivoting -> True|False, "
            "TargetStructure -> \"Dense\"|\"Structured\".\n");
        return NULL;
    }
    /* "Structured" output requires builtin heads (OrthogonalMatrix,
     * UnitaryMatrix, UpperTriangularMatrix, PermutationMatrix) that
     * Mathilda doesn't yet provide.  Leave the call unevaluated so
     * the user sees the unstructured QR via TargetStructure -> "Dense"
     * or the default. */
    if (opts.target_structure == QR_TS_STRUCTURED) return NULL;

    Expr* m = res->data.function.args[0];

    int64_t dims[64];
    int trank = get_tensor_dims(m, dims);
    if (trank != 2 || dims[0] == 0 || dims[1] == 0) {
        char* s = expr_to_string(m);
        fprintf(stderr,
                "QRDecomposition::matrix: Argument %s at position 1 is "
                "not a non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }
    int n = (int)dims[0];
    int p = (int)dims[1];

    return qr_dispatch(m, n, p, &opts);
}

void qrdecomp_init(void) {
    symtab_add_builtin("QRDecomposition", builtin_qrdecomposition);
    symtab_get_def("QRDecomposition")->attributes |= ATTR_PROTECTED;
}
