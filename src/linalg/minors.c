/*
 * src/linalg/minors.c -- Minors and Adjugate.
 *
 *   Minors[m]          the (n-1) x (n-1) minors, n = min(rows, cols)
 *   Minors[m, k]       the k x k minors: a C(rows, k) x C(cols, k) matrix whose
 *                      (r, c) entry is Det of the rows in the r-th and the
 *                      columns in the c-th k-subset, subsets in LEXICOGRAPHIC
 *                      order -- the reference's convention (probed:
 *                      Minors[{{a,b,c},{d,e,f},{g,h,i}}][[1]] is
 *                      {a e - b d, a f - c d, b f - c e}).  k = 0 gives {{1}},
 *                      k < 0 or k > n gives {}.
 *   Minors[m, k, f]    f applied to each k x k submatrix instead of Det.
 *   Adjugate[m]        the classical adjoint: Adjugate[m] . m == Det[m] I,
 *                      defined for a singular m too.
 *
 * Dispatch:
 *   machine (List, packed, NDArray)  -> raw-double kernels (minors_machine:
 *       one partial-pivot LU determinant per submatrix; adjugate_machine: the
 *       SVD form adj(U S V^H) = det(U) det(V^H) V adj(S) U^H, stable for a
 *       singular matrix where Det[m] Inverse[m] is not).  These are the ONLY
 *       machine helpers, shared by List / packed / NDArray / Compile[].
 *   exact rational, nonsingular, n >= 4  -> Det[m] Inverse[m] (FLINT fmpq).
 *   everything else (exact, symbolic, MPFR) -> cofactors through Det, at the
 *       input precision.
 */

#include "la_dense.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "lapack.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------ *
 *  k-subsets in lexicographic order                                        *
 * ------------------------------------------------------------------------ */

/* C(n, k), saturating at INT64_MAX/2 (a size that large is refused anyway). */
static int64_t mn_binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    int64_t r = 1;
    for (int i = 1; i <= k; i++) {
        if (r > INT64_MAX / 2 / (n - k + i)) return INT64_MAX / 2;
        r = r * (n - k + i) / i;
    }
    return r;
}

/* All k-subsets of {0..n-1}, lexicographic, flattened (count * k ints). */
static int* mn_subsets(int n, int k, int64_t* count) {
    *count = mn_binom(n, k);
    int* out = (int*)malloc(sizeof(int) * (size_t)(*count * (k > 0 ? k : 1)));
    if (!out) return NULL;
    if (k == 0) return out;
    int* cur = (int*)malloc(sizeof(int) * (size_t)k);
    for (int i = 0; i < k; i++) cur[i] = i;
    for (int64_t s = 0; s < *count; s++) {
        memcpy(out + s * k, cur, sizeof(int) * (size_t)k);
        int i = k - 1;
        while (i >= 0 && cur[i] == n - k + i) i--;
        if (i < 0) break;
        cur[i]++;
        for (int j = i + 1; j < k; j++) cur[j] = cur[j - 1] + 1;
    }
    free(cur);
    return out;
}

/* ------------------------------------------------------------------------ *
 *  Machine kernels                                                         *
 * ------------------------------------------------------------------------ */

/* Determinant of the k x k row-major work matrix W (destroyed), partial
 * pivoting.  Real: W has k*k doubles; complex: 2*k*k interleaved. */
static void mn_det(double* W, int k, bool cplx, double* dre, double* dim) {
    double pr = 1.0, pi = 0.0;
    for (int c = 0; c < k; c++) {
        int piv = c;
        double best = -1.0;
        for (int r = c; r < k; r++) {
            size_t o = (size_t)r * k + c;
            double v = cplx ? hypot(W[2 * o], W[2 * o + 1]) : fabs(W[o]);
            if (v > best) { best = v; piv = r; }
        }
        if (best == 0.0) { *dre = 0.0; *dim = 0.0; return; }
        if (piv != c) {
            for (int j = 0; j < k; j++) {
                size_t a = (size_t)c * k + j, b = (size_t)piv * k + j;
                if (cplx) {
                    double t0 = W[2 * a], t1 = W[2 * a + 1];
                    W[2 * a] = W[2 * b]; W[2 * a + 1] = W[2 * b + 1];
                    W[2 * b] = t0; W[2 * b + 1] = t1;
                } else { double t = W[a]; W[a] = W[b]; W[b] = t; }
            }
            pr = -pr; pi = -pi;
        }
        size_t cc = (size_t)c * k + c;
        double ar = cplx ? W[2 * cc] : W[cc], ai = cplx ? W[2 * cc + 1] : 0.0;
        double nr = pr * ar - pi * ai, ni = pr * ai + pi * ar;
        pr = nr; pi = ni;
        double den = ar * ar + ai * ai;
        for (int r = c + 1; r < k; r++) {
            size_t rc = (size_t)r * k + c;
            double fr, fi;
            if (cplx) {
                double br = W[2 * rc], bi = W[2 * rc + 1];
                fr = (br * ar + bi * ai) / den; fi = (bi * ar - br * ai) / den;
            } else { fr = W[rc] / ar; fi = 0.0; }
            for (int j = c + 1; j < k; j++) {
                size_t rj = (size_t)r * k + j, cj = (size_t)c * k + j;
                if (cplx) {
                    double ur = W[2 * cj], ui = W[2 * cj + 1];
                    W[2 * rj]     -= fr * ur - fi * ui;
                    W[2 * rj + 1] -= fr * ui + fi * ur;
                } else W[rj] -= fr * W[cj];
            }
        }
    }
    *dre = pr; *dim = pi;
}

/* Load a machine matrix row-major (real first, else complex). */
static double* mn_load(const Expr* m, int* r, int* c, bool* cplx) {
    double* A = NULL;
    *cplx = false;
    if (na_load_matrix((Expr*)m, false, false, r, c, &A)) return A;
    free(A); A = NULL;
    if (na_load_matrix((Expr*)m, true, false, r, c, &A)) { *cplx = true; return A; }
    free(A);
    return NULL;
}

/* The k x k minors of a machine matrix as a packed matrix, or NULL. */
static Expr* minors_machine(const Expr* m, int k) {
    int r = 0, c = 0;
    bool cplx;
    double* A = mn_load(m, &r, &c, &cplx);
    if (!A) return NULL;
    int64_t nr, nc;
    int* rs = mn_subsets(r, k, &nr);
    int* cs = mn_subsets(c, k, &nc);
    size_t comps = cplx ? 2 : 1;
    double* out = (double*)malloc(sizeof(double) * comps * (size_t)(nr * nc));
    double* W = (double*)malloc(sizeof(double) * comps * (size_t)(k * k > 0 ? k * k : 1));
    double* L = (double*)malloc(sizeof(double) * comps * (size_t)(k * k > 0 ? k * k : 1));
    int* ipiv = (int*)malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
    for (int64_t a = 0; a < nr; a++)
        for (int64_t b = 0; b < nc; b++) {
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++) {
                    size_t src = (size_t)rs[a * k + i] * (size_t)c + (size_t)cs[b * k + j];
                    size_t dst = (size_t)i * k + j;
                    if (cplx) { W[2 * dst] = A[2 * src]; W[2 * dst + 1] = A[2 * src + 1]; }
                    else W[dst] = A[src];
                }
            /* the reference takes each minor with Det, which warns Det::luc
             * for a badly conditioned (nonsingular) submatrix: the same
             * infinity-norm rcond test, on a column-major copy */
            if (k >= 2 && mathilda_lapack_probe()) {
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < k; j++) {
                        size_t src = (size_t)i * k + j, dst = (size_t)j * k + i;
                        if (cplx) { L[2 * dst] = W[2 * src]; L[2 * dst + 1] = W[2 * src + 1]; }
                        else L[dst] = W[src];
                    }
                double anorm = nd_mat_norm('I', L, k, cplx);
                int info = cplx ? mat_lapack_zgetrf(k, k, L, k, ipiv)
                                : mat_lapack_dgetrf(k, k, L, k, ipiv);
                if (info == 0 && nd_lu_rcond('I', L, k, cplx, anorm) < ND_LUC_RCOND) {
                    Expr* sub = na_build_matrix(W, k, k, cplx, false);
                    Expr* box = la_boxed_copy(sub);
                    nd_warn_luc("Det", box);
                    expr_free(box); expr_free(sub);
                }
            }
            double dr, di;
            mn_det(W, k, cplx, &dr, &di);
            size_t o = (size_t)(a * nc + b);
            if (cplx) { out[2 * o] = dr; out[2 * o + 1] = di; } else out[o] = dr;
        }
    Expr* res = na_build_matrix_as(out, (int)nr, (int)nc, cplx, false, na_result_presentation(m));
    free(out); free(W); free(L); free(ipiv); free(rs); free(cs); free(A);
    return res;
}

/* Complex determinant of an n x n column-major LAPACK buffer (destroyed). */
static void mn_det_cm(double* M, int n, bool cplx, double* dre, double* dim) {
    /* transpose-invariant: det of the column-major buffer read row-major */
    mn_det(M, n, cplx, dre, dim);
}

/* Machine adjugate through the SVD (see the file comment), or NULL. */
static Expr* adjugate_machine(const Expr* m) {
    int r = 0, c = 0;
    bool cplx;
    double* Arm = mn_load(m, &r, &c, &cplx);
    if (!Arm || r != c || r < 1) { free(Arm); return NULL; }
    int n = r;
    size_t comps = cplx ? 2 : 1, nn = (size_t)n * (size_t)n;
    double* out = (double*)calloc(comps * nn, sizeof(double));
    if (n == 1) {
        out[0] = 1.0;
    } else {
        /* column-major copy for LAPACK */
        double* A = (double*)malloc(sizeof(double) * comps * nn);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                size_t s = (size_t)i * n + j, d = (size_t)i + (size_t)j * n;
                if (cplx) { A[2 * d] = Arm[2 * s]; A[2 * d + 1] = Arm[2 * s + 1]; }
                else A[d] = Arm[s];
            }
        double* S  = (double*)malloc(sizeof(double) * (size_t)n);
        double* U  = (double*)malloc(sizeof(double) * comps * nn);
        double* VT = (double*)malloc(sizeof(double) * comps * nn);
        int info = cplx ? mat_lapack_zgesdd('A', n, n, A, n, S, U, n, VT, n)
                        : mat_lapack_dgesdd('A', n, n, A, n, S, U, n, VT, n);
        free(A);
        if (info != 0) { free(S); free(U); free(VT); free(out); free(Arm); return NULL; }
        /* adj(S)_k = prod_{j != k} s_j */
        double* D = (double*)malloc(sizeof(double) * (size_t)n);
        for (int k = 0; k < n; k++) {
            double p = 1.0;
            for (int j = 0; j < n; j++) if (j != k) p *= S[j];
            D[k] = p;
        }
        /* coef = det(U) det(VT) */
        double* T = (double*)malloc(sizeof(double) * comps * nn);
        double ur, ui, vr, vi;
        memcpy(T, U, sizeof(double) * comps * nn);  mn_det_cm(T, n, cplx, &ur, &ui);
        memcpy(T, VT, sizeof(double) * comps * nn); mn_det_cm(T, n, cplx, &vr, &vi);
        free(T);
        double cr = ur * vr - ui * vi, ci = ur * vi + ui * vr;
        /* adj[i][j] = coef * sum_k conj(VT[k][i]) D_k conj(U[j][k]) */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double sr = 0.0, si = 0.0;
                for (int k = 0; k < n; k++) {
                    size_t vki = (size_t)k + (size_t)i * n;   /* VT[k][i] col-major */
                    size_t ujk = (size_t)j + (size_t)k * n;   /* U[j][k]           */
                    if (cplx) {
                        double ar = VT[2 * vki], ai = -VT[2 * vki + 1];
                        double br = U[2 * ujk],  bi = -U[2 * ujk + 1];
                        sr += D[k] * (ar * br - ai * bi);
                        si += D[k] * (ar * bi + ai * br);
                    } else sr += VT[vki] * D[k] * U[ujk];
                }
                size_t o = (size_t)i * n + j;
                if (cplx) { out[2 * o] = cr * sr - ci * si; out[2 * o + 1] = cr * si + ci * sr; }
                else out[o] = cr * sr;
            }
        free(D); free(S); free(U); free(VT);
    }
    Expr* res = na_build_matrix_as(out, n, n, cplx, false, na_result_presentation(m));
    free(out); free(Arm);
    return res;
}

Expr* la_minors_machine(const Expr* m) {
    int r = 0, c = 0;
    if (!is_ndarray(m) || !la_matrix_dims(m, &r, &c) || r < 1 || c < 1) return NULL;
    int n = r < c ? r : c;
    if (n - 1 < 1) return NULL;   /* {{1}} is exact: let the interpreter answer */
    return minors_machine(m, n - 1);
}

Expr* la_adjugate_machine(const Expr* m) {
    int r = 0, c = 0;
    if (!is_ndarray(m) || !la_matrix_dims(m, &r, &c) || r < 1 || r != c) return NULL;
    return adjugate_machine(m);
}

/* ------------------------------------------------------------------------ *
 *  Generic kernels (exact, symbolic, MPFR, user function)                  *
 * ------------------------------------------------------------------------ */

/* The submatrix of boxed `m` on the given rows/cols (k each). */
static Expr* mn_sub(const Expr* m, const int* rows, const int* cols, int k) {
    Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)(k * k > 0 ? k * k : 1));
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            cells[(size_t)i * k + j] = expr_copy(la_at(m, rows[i], cols[j]));
    Expr* s = la_matrix_from(cells, k, k);
    free(cells);
    return s;
}

static Expr* minors_generic(const Expr* m, int r, int c, int k, const Expr* f) {
    int64_t nr, nc;
    int* rs = mn_subsets(r, k, &nr);
    int* cs = mn_subsets(c, k, &nc);
    Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)(nr * nc > 0 ? nr * nc : 1));
    for (int64_t a = 0; a < nr; a++)
        for (int64_t b = 0; b < nc; b++) {
            Expr* sub = mn_sub(m, rs + a * k, cs + b * k, k);
            Expr* call = f ? expr_new_function(expr_copy((Expr*)f), (Expr*[]){ sub }, 1)
                           : la_call1(SYM_Det, sub);
            cells[a * nc + b] = eval_and_free(call);
        }
    Expr* out = la_matrix_from(cells, (int)nr, (int)nc);
    free(cells); free(rs); free(cs);
    return out;
}

/* True for an exact matrix of Integers / Rationals (FLINT's domain). */
static bool mn_all_rational(const Expr* m, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            const Expr* e = la_at(m, i, j);
            if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) continue;
            if (head_is(e, SYM_Rational)) continue;
            return false;
        }
    return true;
}

static Expr* adjugate_generic(const Expr* m, int n) {
    if (n == 1) return la_matrix_from((Expr*[]){ expr_new_integer(1) }, 1, 1);
    if (n >= 4 && mn_all_rational(m, n)) {
        Expr* d = la_ev1(SYM_Det, expr_copy((Expr*)m));
        if (!(d->type == EXPR_INTEGER && d->data.integer == 0)) {
            Expr* inv = la_ev1(SYM_Inverse, expr_copy((Expr*)m));
            if (head_is(inv, SYM_List)) return la_ev2(SYM_Times, d, inv);
            expr_free(inv);
        }
        expr_free(d);
    }
    /* cofactors: adj[i][j] = (-1)^(i+j) Det[m without row j and column i] */
    int k = n - 1;
    int* rows = (int*)malloc(sizeof(int) * (size_t)k);
    int* cols = (int*)malloc(sizeof(int) * (size_t)k);
    Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n * (size_t)n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            for (int t = 0, w = 0; t < n; t++) if (t != j) rows[w++] = t;
            for (int t = 0, w = 0; t < n; t++) if (t != i) cols[w++] = t;
            Expr* d = la_ev1(SYM_Det, mn_sub(m, rows, cols, k));
            if ((i + j) & 1) d = la_ev2(SYM_Times, expr_new_integer(-1), d);
            cells[(size_t)i * n + j] = d;
        }
    free(rows); free(cols);
    Expr* out = la_matrix_from(cells, n, n);
    free(cells);
    return out;
}

/* ------------------------------------------------------------------------ *
 *  Entry points                                                            *
 * ------------------------------------------------------------------------ */

static void mn_matrix_msg(const char* head, const Expr* m, bool square) {
    char* s = la_str(m);
    if (square)
        mth_message(head, "matsq", "Argument %s at position 1 is not a non-empty square matrix.", s);
    else
        mth_message(head, "matrix", "Argument %s at position 1 is not a non-empty rectangular matrix.", s);
    free(s);
}

static Expr* minors_run(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return builtin_arg_error("Minors", argc, 1, 3);
    const Expr* m = res->data.function.args[0];
    int r = 0, c = 0;
    if (!la_matrix_dims(m, &r, &c) || r < 1 || c < 1) {
        mn_matrix_msg("Minors", m, false);
        return NULL;
    }
    int n = r < c ? r : c;
    int64_t k = n - 1;
    if (argc >= 2) {
        const Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) return NULL;     /* silent, as the reference */
        k = ke->data.integer;
    }
    const Expr* f = argc == 3 ? res->data.function.args[2] : NULL;
    if (k < 0 || k > n) return la_empty_list();
    if (k == 0) {
        Expr* one = f ? eval_and_free(expr_new_function(expr_copy((Expr*)f),
                            (Expr*[]){ la_empty_list() }, 1))
                      : expr_new_integer(1);
        return la_matrix_from(&one, 1, 1);
    }
    if (!f && la_prec_class(m, NULL) == LA_MACHINE) {
        Expr* out = minors_machine(m, (int)k);
        if (out) return out;
    }
    Expr* boxed = la_boxed_copy(m);
    Expr* out = minors_generic(boxed, r, c, (int)k, f);
    expr_free(boxed);
    return out;
}

Expr* ndla_minors(Expr* res) { return minors_run(res); }

Expr* builtin_minors(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_minors(res);
    return minors_run(res);
}

static Expr* adjugate_run(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return builtin_arg_error("Adjugate", argc, 1, 1);
    if (argc > 1) {
        char* s = la_str(res);
        mth_message("Adjugate", "intnm", "Non-negative machine-sized integer expected at "
                    "position 2 in %s.", s);
        free(s);
        return NULL;
    }
    const Expr* m = res->data.function.args[0];
    int r = 0, c = 0;
    if (!la_matrix_dims(m, &r, &c) || r < 1 || r != c) {
        mn_matrix_msg("Adjugate", m, true);
        return NULL;
    }
    if (la_prec_class(m, NULL) == LA_MACHINE) {
        Expr* out = adjugate_machine(m);
        if (out) return out;
    }
    Expr* boxed = la_boxed_copy(m);
    Expr* out = adjugate_generic(boxed, r);
    expr_free(boxed);
    return out;
}

Expr* ndla_adjugate(Expr* res) { return adjugate_run(res); }

Expr* builtin_adjugate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (linalg_call_has_ndarray(res)) return ndla_adjugate(res);
    return adjugate_run(res);
}

void minors_init(void) {
    symtab_add_builtin("Minors", builtin_minors);
    symtab_get_def("Minors")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Adjugate", builtin_adjugate);
    symtab_get_def("Adjugate")->attributes |= ATTR_PROTECTED;
}
