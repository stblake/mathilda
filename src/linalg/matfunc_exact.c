/*
 * matfunc_exact.c -- the exact / symbolic engine for the matrix functions.
 *
 * For a square matrix m = S J S^-1 in Jordan form, a primary matrix function is
 *
 *     f(m) = sum over the distinct eigenvalues lambda, k = 0 .. (largest
 *            block of lambda) - 1, of  f^(k)(lambda) / k!  C_{lambda,k},
 *
 *     C_{lambda,k} = sum over the Jordan blocks B of lambda, r = 0 .. |B|-1-k,
 *                    of  S[:, B_r] (S^-1)[B_{r+k}, :]
 *
 * (the Hermite form of the Sylvester / Lagrange interpolation formula).  Each
 * coefficient entry is simplified on its own -- it is an algebraic function of
 * the matrix entries, independent of f -- so the answer reads as a short sum
 * f(lambda_1) c_1 + f(lambda_2) c_2 + ..., the shape the reference prints:
 * MatrixExp[{{1, 1}, {0, 1}}] is {{E, E}, {0, E}} and MatrixFunction[f,
 * {{1, 1}, {0, 1}}] is {{f[1], f'[1]}, {0, f[1]}}.
 *
 * Derivatives come from Derivative[k][f], so an undefined f keeps its
 * f'[lambda] and a non-analytic one (Abs) keeps Abs'[lambda], as in the
 * reference.  A value that is infinite or indeterminate (Log at a zero
 * eigenvalue, Sqrt' at a zero eigenvalue in a Jordan block) means the matrix
 * function does not exist; the caller's message is emitted and NULL returned.
 *
 * A real exact matrix whose eigenvalues are complex gets a REAL MatrixExp:
 * ComplexExpand turns e^(a + i b) pairs into e^a Cos[b], e^a Sin[b], which is
 * how the reference writes MatrixExp[{{0, 1}, {-1, 0}}] ({{Cos[1], Sin[1]},
 * {-Sin[1], Cos[1]}}).
 *
 * Memory: every intermediate Expr is owned and freed here; inputs are borrowed.
 */

#include "matfunc_internal.h"
#include "linalg.h"
#include "eval.h"
#include "expr.h"
#include "print.h"
#include "sym_names.h"
#include "message.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ---- tiny builders ------------------------------------------------------- */

static Expr* mk1(const char* h, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(h), args, 1);
}
static Expr* mk2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), args, 2);
}
static Expr* mk3(const char* h, Expr* a, Expr* b, Expr* c) {
    Expr* args[3] = { a, b, c };
    return expr_new_function(expr_new_symbol(h), args, 3);
}

static bool mf_head_is(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* An n x n List-of-Lists? */
static bool is_square(const Expr* e, int n) {
    if (!mf_head_is(e, SYM_List) || (int)e->data.function.arg_count != n) return false;
    for (int i = 0; i < n; i++) {
        const Expr* r = e->data.function.args[i];
        if (!mf_head_is(r, SYM_List) || (int)r->data.function.arg_count != n) return false;
    }
    return true;
}

static const Expr* at(const Expr* m, int i, int j) {
    return m->data.function.args[i]->data.function.args[j];
}

static bool is_zero_int(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Does `e` contain an infinity or Indeterminate (the "f is not defined here"
 * signal)? */
static bool has_singular(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_ComplexInfinity || e->data.symbol.name == SYM_Indeterminate
            || e->data.symbol.name == SYM_Infinity;
    if (e->type != EXPR_FUNCTION) return false;
    if (mf_head_is(e, SYM_DirectedInfinity)) return true;
    if (has_singular(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_singular(e->data.function.args[i])) return true;
    return false;
}

/* Does `e` contain an explicit Complex number? */
static bool has_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (mf_head_is(e, SYM_Complex)) return true;
    if (has_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_complex(e->data.function.args[i])) return true;
    return false;
}

/* The k-th derivative of f as a function: Function[x, D[f[x], {x, k}]] with
 * a private x (k = 0: f itself).  D, not Derivative[k][f], because the
 * evaluator leaves Derivative[1][Exp] and Derivative[1][Log] unevaluated; D on
 * an undefined f still gives Derivative[k][f][x], so f'[lambda] survives as in
 * the reference. */
Expr* mf_derivative_function(const Expr* f, int k) {
    if (k == 0) return expr_copy((Expr*)f);
    Expr* x = expr_new_symbol("MatrixFunction`Private`x");
    Expr* a[1] = { expr_copy(x) };
    Expr* fx = expr_new_function(expr_copy((Expr*)f), a, 1);
    Expr* spec = mk2(SYM_List, expr_copy(x), expr_new_integer(k));
    Expr* d = eval_and_free(mk2(SYM_D, fx, spec));
    return mk2(SYM_Function, x, d);
}

/* f^(k)(lambda) / k! */
static Expr* fderiv(const Expr* f, const Expr* lambda, int k) {
    Expr* fk = mf_derivative_function(f, k);
    Expr* a[1] = { expr_copy((Expr*)lambda) };
    Expr* call = expr_new_function(fk, a, 1);
    if (k >= 2) {
        int64_t fact = 1;
        for (int j = 2; j <= k; j++) fact *= j;
        call = mk2(SYM_Times, mk2(SYM_Power, expr_new_integer(fact), expr_new_integer(-1)), call);
    }
    return eval_and_free(call);
}

/* Simplify one coefficient entry: Together, then Simplify. */
static Expr* simplify_coef(Expr* c) {
    Expr* t = eval_and_free(mk1(SYM_Together, c));
    return eval_and_free(mk1(SYM_Simplify, t));
}

typedef struct { int start, size, group; } Block;

Expr* mf_exact_apply(const Expr* m, int n, const Expr* f, MfxKind kind, const char* head) {
    Expr* jd = eval_and_free(mk1(SYM_JordanDecomposition, expr_copy((Expr*)m)));
    if (!mf_head_is(jd, SYM_List) || jd->data.function.arg_count != 2
        || !is_square(jd->data.function.args[0], n) || !is_square(jd->data.function.args[1], n)) {
        expr_free(jd);
        return NULL;
    }
    const Expr* S = jd->data.function.args[0];
    const Expr* J = jd->data.function.args[1];
    Expr* Sinv = eval_and_free(mk1(SYM_Inverse, expr_copy((Expr*)S)));
    if (!is_square(Sinv, n)) { expr_free(Sinv); expr_free(jd); return NULL; }

    /* Jordan blocks, grouped by eigenvalue. */
    Block* blocks = (Block*)malloc(sizeof(Block) * (size_t)n);
    const Expr** lam = (const Expr**)malloc(sizeof(Expr*) * (size_t)n);
    int* gmax = (int*)calloc((size_t)n, sizeof(int));
    int nb = 0, ng = 0;
    for (int i = 0; i < n; ) {
        int size = 1;
        while (i + size < n && !is_zero_int(at(J, i + size - 1, i + size))
               && expr_eq(at(J, i + size, i + size), at(J, i, i)))
            size++;
        int g = -1;
        for (int q = 0; q < ng; q++) if (expr_eq(lam[q], at(J, i, i))) { g = q; break; }
        if (g < 0) { g = ng; lam[ng++] = at(J, i, i); }
        if (size > gmax[g]) gmax[g] = size;
        blocks[nb].start = i; blocks[nb].size = size; blocks[nb].group = g; nb++;
        i += size;
    }

    bool real_input = !has_complex(m);
    bool any_complex_eig = false;
    for (int g = 0; g < ng; g++) if (has_complex(lam[g])) any_complex_eig = true;

    /* terms[i][j]: the Plus of coefficient * value pieces. */
    size_t nn = (size_t)n * (size_t)n;
    Expr*** terms = (Expr***)calloc(nn, sizeof(Expr**));
    int* nterm = (int*)calloc(nn, sizeof(int));
    int cap = 0;
    for (int g = 0; g < ng; g++) cap += gmax[g];
    for (size_t q = 0; q < nn; q++) terms[q] = (Expr**)malloc(sizeof(Expr*) * (size_t)(cap > 0 ? cap : 1));

    bool ok = true;
    for (int g = 0; g < ng && ok; g++) {
        for (int k = 0; k < gmax[g] && ok; k++) {
            /* A singular value (Log at 0, Sqrt' at 0) is reported by this
             * engine's own message, so the scalar evaluator's Power::infy /
             * Infinity::indet are muted while the value is formed. */
            mth_msg_suppress_push();
            Expr* val = fderiv(f, lam[g], k);
            mth_msg_suppress_pop();
            if (has_singular(val)) {
                bool zero = lam[g]->type == EXPR_INTEGER && lam[g]->data.integer == 0;
                char* fs = expr_to_string((Expr*)f);
                char* ls = expr_to_string((Expr*)lam[g]);
                if (kind == MFX_LOG) {
                    if (gmax[g] > 1)
                        mth_message(head, "fnanc", "The function %s is not analytic at %s.", fs, ls);
                    else
                        mth_message(head, "fnand", "The function %s is not analytic or defined at %s.", fs, ls);
                } else if (kind == MFX_FUNCTION) {
                    mth_message(head, "nosol", "Cannot compute the matrix function for the function %s.", fs);
                }
                (void)zero;
                free(fs); free(ls);
                expr_free(val);
                ok = false;
                break;
            }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    /* C_{g,k}[i][j] */
                    int cnt = 0;
                    Expr** parts = NULL;
                    for (int b = 0; b < nb; b++) {
                        if (blocks[b].group != g || blocks[b].size <= k) continue;
                        for (int r = 0; r + k < blocks[b].size; r++) {
                            parts = (Expr**)realloc(parts, sizeof(Expr*) * (size_t)(cnt + 1));
                            parts[cnt++] = mk2(SYM_Times,
                                expr_copy((Expr*)at(S, i, blocks[b].start + r)),
                                expr_copy((Expr*)at(Sinv, blocks[b].start + r + k, j)));
                        }
                    }
                    if (cnt == 0) { free(parts); continue; }
                    Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), parts, (size_t)cnt));
                    free(parts);
                    if (is_zero_int(c)) { expr_free(c); continue; }
                    c = simplify_coef(c);
                    if (is_zero_int(c)) { expr_free(c); continue; }
                    size_t q = (size_t)i * (size_t)n + (size_t)j;
                    terms[q][nterm[q]++] = mk2(SYM_Times, c, expr_copy(val));
                }
            expr_free(val);
        }
    }

    Expr* result = NULL;
    if (ok) {
        Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                size_t q = (size_t)i * (size_t)n + (size_t)j;
                Expr* e = nterm[q] == 0 ? expr_new_integer(0)
                        : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms[q], (size_t)nterm[q]));
                nterm[q] = 0;
                if (kind == MFX_EXP && real_input && any_complex_eig) {
                    e = eval_and_free(mk1(SYM_ComplexExpand, e));
                    e = eval_and_free(mk1(SYM_Simplify, e));
                }
                cells[j] = e;
            }
            rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        }
        result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
        free(rows); free(cells);
    }
    for (size_t q = 0; q < nn; q++) {
        for (int t = 0; t < nterm[q]; t++) expr_free(terms[q][t]);
        free(terms[q]);
    }
    free(terms); free(nterm);
    free(blocks); free(lam); free(gmax);
    expr_free(Sinv);
    expr_free(jd);
    return result;
}

Expr* mf_eigen_apply(const Expr* m, const Expr* f) {
    Expr* vals = eval_and_free(mk1(SYM_Eigenvalues, expr_copy((Expr*)m)));
    Expr* vecs = eval_and_free(mk1(SYM_Eigenvectors, expr_copy((Expr*)m)));
    int n = mf_head_is(vals, SYM_List) ? (int)vals->data.function.arg_count : -1;
    if (n < 1 || !is_square(vecs, n)) { expr_free(vals); expr_free(vecs); return NULL; }
    Expr* V = eval_and_free(mk1(SYM_Transpose, vecs));
    /* A defective matrix has a dependent eigenvector set: no diagonalisation. */
    Expr* rk = eval_and_free(mk1(SYM_MatrixRank, expr_copy(V)));
    bool full = rk->type == EXPR_INTEGER && rk->data.integer == n;
    expr_free(rk);
    if (!full) { expr_free(V); expr_free(vals); return NULL; }
    Expr* fv = eval_and_free(mk2(SYM_Map, expr_copy((Expr*)f), vals));
    Expr* D = eval_and_free(mk1(SYM_DiagonalMatrix, fv));
    Expr* Vi = eval_and_free(mk1(SYM_Inverse, expr_copy(V)));
    Expr* r = eval_and_free(mk3(SYM_Dot, V, D, Vi));
    if (!is_square(r, n)) { expr_free(r); return NULL; }
    return r;
}
