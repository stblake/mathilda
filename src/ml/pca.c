/* pca.c -- column statistics, Standardize and PrincipalComponents.
 *
 * See pca.h for why the kernels are buffer-level. This file also carries the two
 * builtins, since they are thin: read a matrix, call a kernel, build a matrix.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numarray.h"      /* na_load_matrix / na_build_matrix -- the machine bridge */
#include "lapack.h"        /* mat_lapack_dsyev; a stub returning nonzero without LAPACK */
#include "pca.h"

void ml_column_mean(const double* x, size_t n, size_t dim, double* mean) {
    for (size_t j = 0; j < dim; j++) mean[j] = 0.0;
    if (n == 0) return;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < dim; j++) mean[j] += x[i * dim + j];
    for (size_t j = 0; j < dim; j++) mean[j] /= (double)n;
}

void ml_column_sd(const double* x, size_t n, size_t dim, const double* mean,
                  double* sd) {
    for (size_t j = 0; j < dim; j++) sd[j] = 0.0;
    if (n < 2) return;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < dim; j++) {
            double d = x[i * dim + j] - mean[j];
            sd[j] += d * d;
        }
    for (size_t j = 0; j < dim; j++) sd[j] = sqrt(sd[j] / (double)(n - 1));
}

void ml_standardize(double* x, size_t n, size_t dim, bool rescale) {
    double* mean = malloc(sizeof(double) * dim);
    double* sd   = rescale ? malloc(sizeof(double) * dim) : NULL;
    if (!mean || (rescale && !sd)) { free(mean); free(sd); return; }
    ml_column_mean(x, n, dim, mean);
    if (rescale) ml_column_sd(x, n, dim, mean, sd);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < dim; j++) {
            double v = x[i * dim + j] - mean[j];
            /* A constant column stays at zero rather than becoming NaN -- zero
             * variance carries no information, so "no deviation from the mean" is
             * the honest value, and a NaN would poison every downstream row
             * reduction. */
            if (rescale) v = (sd[j] > 0.0) ? v / sd[j] : 0.0;
            x[i * dim + j] = v;
        }
    free(mean); free(sd);
}

/* Cyclic Jacobi for a real symmetric matrix. Used when LAPACK is not linked.
 *
 * Written here rather than pulled from linalg because the in-tree symmetric solvers
 * are static to their translation units (the same wall that put the GMM fit in
 * src/ml). Jacobi is chosen over anything cleverer for the reason the Cholesky in
 * gmm.c is hand-written: `dim` is a feature count, this runs once per call, and the
 * alternative is a conditional dependency. */
static bool ml_jacobi(double* a, size_t dim, double* eval, double* q) {
    for (size_t i = 0; i < dim; i++)
        for (size_t j = 0; j < dim; j++) q[i * dim + j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0.0;
        for (size_t i = 0; i < dim; i++)
            for (size_t j = i + 1; j < dim; j++) off += a[i * dim + j] * a[i * dim + j];
        if (off <= 1e-30) break;
        for (size_t p = 0; p < dim; p++)
            for (size_t r = p + 1; r < dim; r++) {
                double apr = a[p * dim + r];
                if (fabs(apr) < 1e-300) continue;
                double app = a[p * dim + p], arr = a[r * dim + r];
                double theta = 0.5 * (arr - app) / apr;
                double t = (theta >= 0.0 ? 1.0 : -1.0)
                         / (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (size_t m = 0; m < dim; m++) {
                    double amp = a[m * dim + p], amr = a[m * dim + r];
                    a[m * dim + p] = c * amp - s * amr;
                    a[m * dim + r] = s * amp + c * amr;
                }
                for (size_t m = 0; m < dim; m++) {
                    double apm = a[p * dim + m], arm = a[r * dim + m];
                    a[p * dim + m] = c * apm - s * arm;
                    a[r * dim + m] = s * apm + c * arm;
                }
                for (size_t m = 0; m < dim; m++) {
                    double qmp = q[m * dim + p], qmr = q[m * dim + r];
                    q[m * dim + p] = c * qmp - s * qmr;
                    q[m * dim + r] = s * qmp + c * qmr;
                }
            }
    }
    for (size_t i = 0; i < dim; i++) eval[i] = a[i * dim + i];
    return true;
}

bool ml_sym_eigen_desc(const double* a, size_t dim, double* eval, double* evec) {
    if (!a || dim == 0 || !eval || !evec) return false;
    double* w = malloc(sizeof(double) * dim);
    double* m = malloc(sizeof(double) * dim * dim);
    if (!w || !m) { free(w); free(m); return false; }
    memcpy(m, a, sizeof(double) * dim * dim);

    /* dsyev takes a column-major matrix, but the input is symmetric, so its
     * row-major and column-major layouts are the same bytes -- no transpose needed.
     * Eigenvectors come back as COLUMNS of m, ascending by eigenvalue. */
    bool ok = false;
    if (mat_lapack_dsyev((int)dim, m, (int)dim, w) == 0) {
        for (size_t j = 0; j < dim; j++)                    /* column j -> row j */
            for (size_t r = 0; r < dim; r++)
                evec[j * dim + r] = m[r + j * dim];
        for (size_t j = 0; j < dim; j++) eval[j] = w[j];
        ok = true;
    } else {
        double* q = malloc(sizeof(double) * dim * dim);
        if (q) {
            memcpy(m, a, sizeof(double) * dim * dim);
            ok = ml_jacobi(m, dim, w, q);
            if (ok)
                for (size_t j = 0; j < dim; j++)
                    for (size_t r = 0; r < dim; r++)
                        evec[j * dim + r] = q[r * dim + j];
            for (size_t j = 0; j < dim && ok; j++) eval[j] = w[j];
        }
        free(q);
    }
    free(m);
    if (!ok) { free(w); return false; }

    /* Reverse to DESCENDING. dsyev gives ascending and Jacobi gives no order at all,
     * so sort rather than assume: selection sort, since dim is small and this must be
     * a total order regardless of which backend ran. */
    for (size_t i = 0; i < dim; i++) {
        size_t best = i;
        for (size_t j = i + 1; j < dim; j++) if (eval[j] > eval[best]) best = j;
        if (best != i) {
            double tv = eval[i]; eval[i] = eval[best]; eval[best] = tv;
            for (size_t r = 0; r < dim; r++) {
                double t = evec[i * dim + r];
                evec[i * dim + r] = evec[best * dim + r];
                evec[best * dim + r] = t;
            }
        }
    }

    /* Canonicalise signs. An eigenvector is defined only up to sign, and LAPACK and
     * Jacobi disagree on which one they hand back -- so without this the same input
     * gives sign-flipped components depending on how the binary was linked, and no
     * output could be pinned in a test. Flip each row so its largest-magnitude
     * component is positive; ties go to the earliest index, which keeps it total. */
    for (size_t j = 0; j < dim; j++) {
        size_t big = 0; double bv = -1.0;
        for (size_t r = 0; r < dim; r++) {
            double av = fabs(evec[j * dim + r]);
            if (av > bv) { bv = av; big = r; }
        }
        if (evec[j * dim + big] < 0.0)
            for (size_t r = 0; r < dim; r++) evec[j * dim + r] = -evec[j * dim + r];
    }

    free(w);
    return true;
}

bool ml_pca(const double* x, size_t n, size_t dim, bool correlation,
            double* out, double* eval, double* evec) {
    if (!x || n == 0 || dim == 0) return false;
    double* xc = malloc(sizeof(double) * n * dim);
    double* cv = malloc(sizeof(double) * dim * dim);
    double* ev = malloc(sizeof(double) * dim);
    double* vc = malloc(sizeof(double) * dim * dim);
    if (!xc || !cv || !ev || !vc) { free(xc); free(cv); free(ev); free(vc); return false; }
    memcpy(xc, x, sizeof(double) * n * dim);
    ml_standardize(xc, n, dim, correlation);

    /* Covariance of the centred data, or correlation when the columns were also
     * scaled to unit variance -- the two differ only in that scaling, which is why
     * one code path serves both. */
    double den = (n > 1) ? (double)(n - 1) : 1.0;
    for (size_t a = 0; a < dim; a++)
        for (size_t b = 0; b <= a; b++) {
            double s = 0.0;
            for (size_t i = 0; i < n; i++) s += xc[i * dim + a] * xc[i * dim + b];
            cv[a * dim + b] = cv[b * dim + a] = s / den;
        }

    bool ok = ml_sym_eigen_desc(cv, dim, ev, vc);
    if (ok) {
        if (eval) memcpy(eval, ev, sizeof(double) * dim);
        if (evec) memcpy(evec, vc, sizeof(double) * dim * dim);
        if (out)
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < dim; j++) {
                    double s = 0.0;
                    for (size_t r = 0; r < dim; r++)
                        s += xc[i * dim + r] * vc[j * dim + r];
                    out[i * dim + j] = s;
                }
    }
    free(xc); free(cv); free(ev); free(vc);
    return ok;
}

/* ------------------------------------------------------------------------- */
/* Dimensionality reduction                                                   */
/* ------------------------------------------------------------------------- */

/* Classical MDS builds an n x n matrix, so it needs a ceiling the other two do not.
 * Same order as the Spectral method's, and for the same reason: the memory is
 * quadratic in the sample size, and refusing is better than appearing to hang. */
#define ML_MDS_MAX_N 2000

bool ml_reduce(const double* x, size_t n, size_t dim, size_t target,
               MlReduceMethod method, double* out) {
    if (!x || !out || n == 0 || dim == 0 || target == 0) return false;

    if (method == ML_REDUCE_MDS) {
        if (n > ML_MDS_MAX_N) return false;
        /* Classical (Torgerson) scaling. B = -0.5 * J D^2 J with J = I - 11'/n,
         * which is the double-centring that turns squared distances back into inner
         * products; its eigenvectors scaled by sqrt(eigenvalue) are the coordinates.
         *
         * Worth knowing, and asserted in the tests: on EUCLIDEAN distances this is
         * mathematically the same embedding as PCA. That is not a redundancy -- it is
         * the check that both are right -- and it is also why MDS earns its keep only
         * when the distances come from somewhere else. */
        double* b = malloc(sizeof(double) * n * n);
        double* ev = malloc(sizeof(double) * n);
        double* vc = malloc(sizeof(double) * n * n);
        double* rm = malloc(sizeof(double) * n);
        if (!b || !ev || !vc || !rm) { free(b); free(ev); free(vc); free(rm); return false; }

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j <= i; j++) {
                double sq = 0.0;
                for (size_t a = 0; a < dim; a++) {
                    double d = x[i * dim + a] - x[j * dim + a];
                    sq += d * d;
                }
                b[i * n + j] = b[j * n + i] = sq;
            }
        double gm = 0.0;
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t j = 0; j < n; j++) s += b[i * n + j];
            rm[i] = s / (double)n;
            gm += s;
        }
        gm /= (double)(n * n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                b[i * n + j] = -0.5 * (b[i * n + j] - rm[i] - rm[j] + gm);

        bool ok = ml_sym_eigen_desc(b, n, ev, vc);
        if (ok)
            for (size_t i = 0; i < n; i++)
                for (size_t t = 0; t < target; t++) {
                    /* A non-positive eigenvalue means that axis carries no real
                     * structure, so it is zero rather than the square root of a
                     * negative number. */
                    double lam = (t < n && ev[t] > 0.0) ? sqrt(ev[t]) : 0.0;
                    out[i * target + t] = (t < n) ? vc[t * n + i] * lam : 0.0;
                }
        free(b); free(ev); free(vc); free(rm);
        return ok;
    }

    /* PCA and LSA differ only in whether the columns are centred first, so one path
     * serves both -- and that single `if` is the entire mathematical difference
     * between "principal components" and "truncated SVD". */
    size_t keep = (target < dim) ? target : dim;
    double* xc = malloc(sizeof(double) * n * dim);
    double* cv = malloc(sizeof(double) * dim * dim);
    double* ev = malloc(sizeof(double) * dim);
    double* vc = malloc(sizeof(double) * dim * dim);
    if (!xc || !cv || !ev || !vc) { free(xc); free(cv); free(ev); free(vc); return false; }
    memcpy(xc, x, sizeof(double) * n * dim);
    if (method == ML_REDUCE_PCA) ml_standardize(xc, n, dim, false);

    for (size_t a = 0; a < dim; a++)
        for (size_t b2 = 0; b2 <= a; b2++) {
            double s = 0.0;
            for (size_t i = 0; i < n; i++) s += xc[i * dim + a] * xc[i * dim + b2];
            cv[a * dim + b2] = cv[b2 * dim + a] = s;
        }

    bool ok = ml_sym_eigen_desc(cv, dim, ev, vc);
    if (ok)
        for (size_t i = 0; i < n; i++)
            for (size_t t = 0; t < target; t++) {
                if (t >= keep) { out[i * target + t] = 0.0; continue; }
                double s = 0.0;
                for (size_t r = 0; r < dim; r++) s += xc[i * dim + r] * vc[t * dim + r];
                out[i * target + t] = s;
            }
    free(xc); free(cv); free(ev); free(vc);
    return ok;
}

/* ------------------------------------------------------------------------- */
/* Builtins                                                                   */
/* ------------------------------------------------------------------------- */

/* Build a plain List (of Lists) of machine reals.
 *
 * NOT na_build_matrix, deliberately, and this is a correctness point rather than a
 * style one. na_build_matrix hands back a VISIBLE NDArray, whose head is NDArray, so
 * `Standardize[data] === {{...}}` compares False against the literal a user would
 * write -- while Inverse, Dot and LinearSolve all return head List and compare True.
 * Building a List keeps this surface consistent with the rest of the system and lets
 * the evaluator's own packing gate decide whether the result is held as a buffer,
 * which is exactly what those three rely on. */
static Expr* ml_list_of_reals(const double* v, size_t n) {
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) a[i] = expr_new_real(v[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), a, n);
    free(a);
    return out;
}

static Expr* ml_list_matrix(const double* x, size_t n, size_t dim) {
    Expr** rows = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!rows) return NULL;
    for (size_t i = 0; i < n; i++) rows[i] = ml_list_of_reals(x + i * dim, dim);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    return out;
}

/* Read arg 0 as an n x dim machine matrix, or as a single vector treated as one
 * column of n observations -- Standardize[{1, 2, 3}] is a list of scalars, which is
 * n observations of ONE variable, not one observation of three. */
static bool ml_read_data(Expr* e, size_t* n, size_t* dim, double** buf,
                         bool* was_vector) {
    int r = 0, c = 0; double* b = NULL;
    if (na_load_matrix(e, 0, 0, &r, &c, &b) && r > 0 && c > 0) {
        *n = (size_t)r; *dim = (size_t)c; *buf = b; *was_vector = false;
        return true;
    }
    int len = 0; double* v = NULL;
    if (na_load_vector(e, 0, &len, &v) && len > 0) {
        *n = (size_t)len; *dim = 1; *buf = v; *was_vector = true;
        return true;
    }
    return false;
}

static Expr* builtin_standardize(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return NULL;
    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    ml_standardize(x, n, dim, true);
    Expr* out = vec ? ml_list_of_reals(x, n) : ml_list_matrix(x, n, dim);
    free(x);
    return out;
}

static Expr* builtin_principal_components(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    /* Method -> "Covariance" (default) or "Correlation". Anything else declines
     * rather than silently choosing, so a typo is visible. */
    bool corr = false;
    if (argc == 2) {
        Expr* o = res->data.function.args[1];
        if (!o || o->type != EXPR_FUNCTION || o->data.function.arg_count != 2)
            return NULL;
        Expr* h = o->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return NULL;
        const char* hn = h->data.symbol.name;
        if (hn != SYM_Rule && hn != SYM_RuleDelayed) return NULL;
        Expr* lhs = o->data.function.args[0];
        Expr* rhs = o->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method)
            return NULL;
        if (!rhs || rhs->type != EXPR_STRING) return NULL;
        const char* m = rhs->data.string;
        if (strcmp(m, "Correlation") == 0)      corr = true;
        else if (strcmp(m, "Covariance") == 0)  corr = false;
        else return NULL;
    }

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    if (vec) { free(x); return NULL; }   /* one variable has no components to rotate */

    double* out = malloc(sizeof(double) * n * dim);
    if (!out) { free(x); return NULL; }
    Expr* r = NULL;
    if (ml_pca(x, n, dim, corr, out, NULL, NULL))
        r = ml_list_matrix(out, n, dim);
    free(out); free(x);
    return r;
}

/* DimensionReduce[data, k] and DimensionReduce[data, k, Method -> m].
 *
 * Wolfram's DimensionReduce can also be called without a target dimension (choosing
 * one itself) and can return a DimensionReducerFunction applicable to LATER data.
 * Neither is implemented here, and the second is the substantive gap: a reusable
 * reducer is a trained model, and the model representation is being designed with the
 * Predict family rather than invented twice. The data-in/data-out form is what this
 * provides, and an omitted dimension declines rather than guessing. */
static Expr* builtin_dimension_reduce(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    MlReduceMethod method = ML_REDUCE_PCA;
    if (argc == 3) {
        Expr* o = res->data.function.args[2];
        if (!o || o->type != EXPR_FUNCTION || o->data.function.arg_count != 2)
            return NULL;
        Expr* h = o->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return NULL;
        const char* hn = h->data.symbol.name;
        if (hn != SYM_Rule && hn != SYM_RuleDelayed) return NULL;
        Expr* lhs = o->data.function.args[0];
        Expr* rhs = o->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method)
            return NULL;
        if (!rhs || rhs->type != EXPR_STRING) return NULL;
        const char* m = rhs->data.string;
        if      (strcmp(m, "PrincipalComponentsAnalysis") == 0) method = ML_REDUCE_PCA;
        else if (strcmp(m, "LatentSemanticAnalysis") == 0)      method = ML_REDUCE_LSA;
        else if (strcmp(m, "MultidimensionalScaling") == 0)     method = ML_REDUCE_MDS;
        else return NULL;   /* an unknown method declines rather than defaulting */
    }

    Expr* kexpr = res->data.function.args[1];
    if (!kexpr || kexpr->type != EXPR_INTEGER) return NULL;
    if (kexpr->data.integer <= 0) return NULL;
    size_t target = (size_t)kexpr->data.integer;

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    if (vec) { free(x); return NULL; }

    /* Asking for more dimensions than the data has is a question with no answer, so
     * it declines instead of padding with zeros -- padding would look like a
     * successful reduction to a caller checking only the shape. */
    size_t cap = (method == ML_REDUCE_MDS) ? (n > 0 ? n - 1 : 0) : dim;
    if (target > cap) { free(x); return NULL; }

    double* out = malloc(sizeof(double) * n * target);
    if (!out) { free(x); return NULL; }
    Expr* r = NULL;
    if (ml_reduce(x, n, dim, target, method, out)) {
        Expr** rows = malloc(sizeof(Expr*) * (n ? n : 1));
        if (rows) {
            for (size_t i = 0; i < n; i++) rows[i] = ml_list_of_reals(out + i * target, target);
            r = expr_new_function(expr_new_symbol(SYM_List), rows, n);
            free(rows);
        }
    }
    free(out); free(x);
    return r;
}

void ml_init(void) {
    symtab_add_builtin("DimensionReduce", builtin_dimension_reduce);
    symtab_get_def("DimensionReduce")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DimensionReduce",
        "DimensionReduce[data, k] reduces each row of data to k dimensions. "
        "Method -> \"PrincipalComponentsAnalysis\" (default) centres the columns and "
        "projects onto the leading eigenvectors of the covariance; "
        "\"LatentSemanticAnalysis\" skips the centring, giving a truncated SVD, which "
        "is what a sparse non-negative term-document matrix wants; "
        "\"MultidimensionalScaling\" double-centres the squared distance matrix "
        "(classical Torgerson scaling) and is capped at 2000 rows, its matrix being "
        "n x n. Asking for more dimensions than the data supports returns unevaluated "
        "rather than padding with zeros.");

    symtab_add_builtin("Standardize", builtin_standardize);
    symtab_get_def("Standardize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Standardize",
        "Standardize[data] shifts each column of data to zero mean and rescales it "
        "to unit sample standard deviation (divisor n-1, matching "
        "StandardDeviation). A flat list is treated as n observations of one "
        "variable. A constant column becomes exactly 0 rather than Indeterminate.");

    symtab_add_builtin("PrincipalComponents", builtin_principal_components);
    symtab_get_def("PrincipalComponents")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PrincipalComponents",
        "PrincipalComponents[matrix] gives the rows of matrix in "
        "principal-component coordinates, components ordered by decreasing "
        "variance. Rows are observations and columns are variables. "
        "Method -> \"Correlation\" standardises each variable to unit variance "
        "first, which is what you want when the columns have different units; the "
        "default \"Covariance\" does not.");
}
