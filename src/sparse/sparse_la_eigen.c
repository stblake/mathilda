/* Mathilda -- sparse linear algebra: a few eigenpairs and singular values.
 *
 * Native handlers (symtab_set_sparse_native) for Eigenvalues, Eigenvectors
 * and Eigensystem with a count -- Eigenvalues[s, k], [s, -k], [s, UpTo[k]],
 * the generalized [{a, b}, k] -- and for SingularValueList[s, k] and
 * SingularValueDecomposition[s, k] on a machine SparseArray (real or complex
 * values, zero background). The same driver serves the dense
 * Method -> "Arnoldi" of Eigenvalues / Eigenvectors (spla_eigen_dense_arnoldi,
 * called from src/linalg/eigen_arnoldi.c), which before used a single,
 * un-restarted Arnoldi pass with residuals of 1e-2.
 *
 * Reference semantics (probed; docs/spec/builtins/linear-algebra.md):
 *
 *   - Which algorithm. A sparse machine matrix is converted to dense (and the
 *     whole spectrum taken) when n <= 14, when more than a fifth of it is
 *     stored (5 nnz > n^2), or when the count is more than a fifth of the
 *     dimension (5 |k| > n); the last also says ::arh once n > 100 (for the
 *     singular-value heads: min(m, n) and max(m, n) > 100). Otherwise, and
 *     always with Method -> "Arnoldi", a restarted Krylov method answers
 *     (ARPACK, or the native Krylov-Schur of src/linalg/krylov_schur.c).
 *     The generalized problem {a, b} is dense under the densify cap.
 *   - k > 0: the k largest in magnitude; -k: the k smallest, by
 *     shift-and-invert about 0 ((A - sigma B)^-1 applied through an
 *     SpFactor: CHOLMOD, UMFPACK or dense LAPACK LU). A singular matrix is
 *     shifted by a tiny amount instead (sigma = 1e-7 phi ||A||; the
 *     eigenvalues are still sigma + 1/theta exactly), as the reference also
 *     answers there.
 *   - Method -> {"Arnoldi", sub-options}: "Shift" -> s (the eigenvalues
 *     nearest s, for k and -k alike; ::ssing, unevaluated, when A - s B is
 *     singular), "Criteria" -> "Magnitude" | "RealPart" | "ImaginaryPart" |
 *     "BothEnds" (applied to the shift-inverted values 1/(lambda - s) when a
 *     shift or -k is in force; "BothEnds" only for a real symmetric matrix,
 *     else ::arbsym; ::eslct for another name), "BasisSize" (an integer > 3
 *     and <= n, else ::arbss), "MaxIterations" (::maxit2 when reached: the
 *     current values are returned), "Tolerance", "StartingVector" (::svecnm);
 *     ::moptx for any other name, ::wrgopt for a value of the wrong kind.
 *   - Ordering: by decreasing |lambda| (ties: the more negative real part,
 *     then the positive imaginary part, first -- the dense order), a
 *     conjugate pair split by k as the dense answer splits it. Eigenvectors
 *     have unit 2-norm and are made real and positive in their largest
 *     component (the Krylov sign / phase is otherwise arbitrary; the
 *     reference's own sparse and dense signs differ).
 *   - ::take when |k| exceeds the dimension (unevaluated).
 *   - SingularValueList[s, k]: the k largest (-k: smallest) singular values,
 *     decreasing, from the Hermitian eigenproblem of A^H A (or A A^H, the
 *     smaller), (A^H A)^-1 through a factor of A for -k on a square matrix;
 *     each sigma is refined as ||A v||. SingularValueDecomposition[s] and
 *     [s, k] answer {List, SparseArray, List}: the singular-value matrix is
 *     a SparseArray with background 0. (machine input).
 *
 * Accuracy: every pair the Krylov layer returns has a confirmed residual
 * ||OP x - theta x|| <= ~1e3 eps ||OP|| (see krylov_schur.c); for A itself
 * that is ||A x - lambda x|| <= ~1e-12 ||A||.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "krylov_schur.h"
#include "numarray.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "common.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Densify cap (elements), the gate's. */
#define EIG_DENSE_MAX ((int64_t)1 << 27)

typedef enum { EK_VALUES = 1, EK_VECTORS = 2, EK_SYSTEM = 3, EK_SVL = 4, EK_SVD = 5 } EigKind;

static const char* kind_head(EigKind k) {
    switch (k) {
        case EK_VALUES: return "Eigenvalues";
        case EK_VECTORS: return "Eigenvectors";
        case EK_SYSTEM: return "Eigensystem";
        case EK_SVL: return "SingularValueList";
        default: return "SingularValueDecomposition";
    }
}

/* ---- small helpers -------------------------------------------------------- */

static bool is_str(const Expr* e, const char* s) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, s) == 0;
}
static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_head(const Expr* e, const char* interned) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == interned;
}
static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           (is_head(e, SYM_Rule) || is_head(e, SYM_RuleDelayed));
}
static const Expr* rule_rhs(const Expr* r) { return r->data.function.args[1]; }
static bool key_named(const Expr* k, const char* name) { return is_sym(k, name) || is_str(k, name); }
static bool rule_named(const Expr* r, const char* name) { return key_named(r->data.function.args[0], name); }

static char* out_str(const Expr* e) {
    char* s = expr_to_string_outputform((Expr*)e);
    return s;
}

/* head::tag with one expression spliced into fmt's single %s */
static void msg1(const char* head, const char* tag, const char* fmt, const Expr* e) {
    char* s = out_str(e);
    mth_message(head, tag, fmt, s ? s : "");
    free(s);
}

/* ---- the machine class of an array ---------------------------------------- */

/* A SparseArray whose values are machine numbers (or exact ones mixed with
 * at least one machine number: the contagion makes the matrix a machine
 * one), over a zero background. */
static bool sparse_is_machine(const SparseArray* sa) {
    if (!sa || sa->rank != 2 || !sp_is_zero(sa->bg)) return false;
    if (sa->kind == SPV_F64 || sa->kind == SPV_C128) return true;
    if (sa->kind != SPV_EXPR) return false;
    Expr* const* v = sa->vals;
    bool inexact = false;
    for (int64_t k = 0; k < sa->pat->nnz; k++) {
        CommonInexactInfo ci = common_scan_inexact(v[k]);
        if (ci.has_inexact) {
            if (ci.min_bits > 53) return false;
            inexact = true;
        }
    }
    return inexact;
}

static bool dense_is_machine(const Expr* e) {
    if (e->type == EXPR_NDARRAY) return true;
    CommonInexactInfo ci = common_scan_inexact(e);
    return ci.has_inexact && ci.min_bits <= 53;
}

/* ---- SpMat utilities ---------------------------------------------------- */

/* Allocate an owned SpMat (rp, ci, x in one block, released by
 * spla_view_free through owned_ci -- the spla_from_dense convention). */
static bool spmat_alloc(SpMat* A, int64_t m, int64_t n, int64_t nnz, bool cplx) {
    memset(A, 0, sizeof *A);
    int w = cplx ? 2 : 1;
    size_t cap = (size_t)(nnz > 0 ? nnz : 1);
    size_t bytes = sizeof(int64_t) * ((size_t)m + 1) + sizeof(int64_t) * cap + sizeof(double) * (size_t)w * cap;
    int64_t* blk = malloc(bytes);
    if (!blk) return false;
    A->m = m;
    A->n = n;
    A->nnz = nnz;
    A->rp = blk;
    A->ci = blk + m + 1;
    A->x = (double*)(blk + m + 1 + cap);
    A->cplx = cplx;
    A->owned_ci = blk;
    return true;
}
#define RP_W(A) ((int64_t*)(A)->rp)
#define CI_W(A) ((int64_t*)(A)->ci)
#define X_W(A) ((double*)(A)->x)

/* Promote a real view to complex (owned interleaved copy of the values). */
static bool view_to_complex(SpMat* A) {
    if (A->cplx) return true;
    double* x = malloc(sizeof(double) * 2 * (size_t)(A->nnz > 0 ? A->nnz : 1));
    if (!x) return false;
    for (int64_t k = 0; k < A->nnz; k++) { x[2 * k] = A->x[k]; x[2 * k + 1] = 0.0; }
    if (A->owned_x) free(A->owned_x);
    A->owned_x = x;
    A->x = x;
    A->cplx = true;
    return true;
}

/* Value of entry (i, j), 0 when not stored (binary search; columns sorted). */
static void spmat_get(const SpMat* A, int64_t i, int64_t j, double* re, double* im) {
    int64_t lo = A->rp[i], hi = A->rp[i + 1];
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (A->ci[mid] < j) lo = mid + 1;
        else hi = mid;
    }
    if (lo < A->rp[i + 1] && A->ci[lo] == j) {
        *re = A->cplx ? A->x[2 * lo] : A->x[lo];
        *im = A->cplx ? A->x[2 * lo + 1] : 0.0;
    } else {
        *re = *im = 0.0;
    }
}

/* Symmetric (real) / Hermitian (complex), to the bit. */
static bool spmat_hermitian(const SpMat* A) {
    if (A->m != A->n) return false;
    for (int64_t i = 0; i < A->m; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t j = A->ci[k];
            double ar = A->cplx ? A->x[2 * k] : A->x[k];
            double ai = A->cplx ? A->x[2 * k + 1] : 0.0;
            if (j == i) { if (ai != 0.0) return false; continue; }
            double br, bi;
            spmat_get(A, j, i, &br, &bi);
            if (br != ar || bi != -ai) return false;
        }
    return true;
}

/* max_i sum_j |a_ij| */
static double spmat_norm_inf(const SpMat* A) {
    double mx = 0.0;
    for (int64_t i = 0; i < A->m; i++) {
        double s = 0.0;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++)
            s += A->cplx ? hypot(A->x[2 * k], A->x[2 * k + 1]) : fabs(A->x[k]);
        if (s > mx) mx = s;
    }
    return mx;
}

/* out = A - sigma B (B NULL: the identity). Complex when A or B is, or
 * sigma has an imaginary part. Rows merged in column order. */
static bool spmat_shifted(const SpMat* A, double sr, double si, const SpMat* B, SpMat* out) {
    bool c = A->cplx || (B && B->cplx) || si != 0.0;
    int64_t n = A->n, cap = A->nnz + (B ? B->nnz : A->m);
    if (!spmat_alloc(out, A->m, n, cap, c)) return false;
    int64_t* rp = RP_W(out);
    int64_t* ci = CI_W(out);
    double* x = X_W(out);
    int64_t q = 0;
    rp[0] = 0;
    for (int64_t i = 0; i < A->m; i++) {
        int64_t ka = A->rp[i], ea = A->rp[i + 1];
        int64_t kb = B ? B->rp[i] : 0, eb = B ? B->rp[i + 1] : 0;
        bool diag_done = B != NULL;
        while (ka < ea || kb < eb || !diag_done) {
            int64_t ja = ka < ea ? A->ci[ka] : INT64_MAX;
            int64_t jb = B ? (kb < eb ? B->ci[kb] : INT64_MAX) : (diag_done ? INT64_MAX : i);
            int64_t j = ja < jb ? ja : jb;
            double vr = 0.0, vi = 0.0;
            if (ja == j) {
                vr = A->cplx ? A->x[2 * ka] : A->x[ka];
                vi = A->cplx ? A->x[2 * ka + 1] : 0.0;
                ka++;
            }
            if (jb == j) {
                double br = 1.0, bi = 0.0;
                if (B) {
                    br = B->cplx ? B->x[2 * kb] : B->x[kb];
                    bi = B->cplx ? B->x[2 * kb + 1] : 0.0;
                    kb++;
                } else {
                    diag_done = true;
                }
                vr -= sr * br - si * bi;
                vi -= sr * bi + si * br;
            }
            ci[q] = j;
            if (c) { x[2 * q] = vr; x[2 * q + 1] = vi; }
            else x[q] = vr;
            q++;
        }
        rp[i + 1] = q;
    }
    out->nnz = q;
    return true;
}

/* out = A^T (conj: A^H), CSR; counting sort by column. */
static bool spmat_transpose(const SpMat* A, bool conj, SpMat* out) {
    if (!spmat_alloc(out, A->n, A->m, A->nnz, A->cplx)) return false;
    int64_t* rp = RP_W(out);
    int64_t* ci = CI_W(out);
    double* x = X_W(out);
    memset(rp, 0, sizeof(int64_t) * (size_t)(A->n + 1));
    for (int64_t k = 0; k < A->nnz; k++) rp[A->ci[k] + 1]++;
    for (int64_t j = 0; j < A->n; j++) rp[j + 1] += rp[j];
    int64_t* pos = malloc(sizeof(int64_t) * (size_t)(A->n > 0 ? A->n : 1));
    if (!pos) { spla_view_free(out); return false; }
    memcpy(pos, rp, sizeof(int64_t) * (size_t)A->n);
    for (int64_t i = 0; i < A->m; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t j = A->ci[k], d = pos[j]++;
            ci[d] = i;
            if (A->cplx) { x[2 * d] = A->x[2 * k]; x[2 * d + 1] = conj ? -A->x[2 * k + 1] : A->x[2 * k + 1]; }
            else x[d] = A->x[k];
        }
    free(pos);
    return true;
}

/* ---- the operator ------------------------------------------------------- */

typedef enum {
    OPM_A,        /* A x                      */
    OPM_FINV,     /* F^-1 x                   */
    OPM_FINV_B,   /* F^-1 B x   (pencil, shift-invert) */
    OPM_BINV_A,   /* F_B^-1 A x (pencil, regular)      */
    OPM_NORMAL,   /* A2 (A1 x): A^H A or A A^H through the stored adjoint */
    OPM_INV_AHA   /* (A^H A)^-1 x = F^-1 F^-H x        */
} OpMode;

typedef struct {
    OpMode mode;
    const SpMat* A1;
    const SpMat* A2;
    SpFactor* F;
    double* t1;
    double* t2;
    int64_t tn;      /* scratch length in doubles */
    bool c;
} OpCtx;

static bool op_apply(void* vc, const double* x, double* y) {
    OpCtx* o = vc;
    switch (o->mode) {
        case OPM_A:
            spla_spmv(o->A1, x, y, false);
            return true;
        case OPM_FINV:
            return spf_solve(o->F, x, y, 1, false);
        case OPM_FINV_B:
        case OPM_BINV_A:
            spla_spmv(o->A1, x, o->t1, false);
            return spf_solve(o->F, o->t1, y, 1, false);
        case OPM_NORMAL:
            spla_spmv(o->A1, x, o->t1, false);
            spla_spmv(o->A2, o->t1, y, false);
            return true;
        case OPM_INV_AHA: {
            /* z = A^-H x: A^T conj(z) = conj(x) */
            int64_t n = o->F ? spf_n(o->F) : 0;
            int w = o->c ? 2 : 1;
            memcpy(o->t1, x, sizeof(double) * (size_t)w * (size_t)n);
            if (o->c) for (int64_t i = 0; i < n; i++) o->t1[2 * i + 1] = -o->t1[2 * i + 1];
            if (!spf_solve(o->F, o->t1, o->t2, 1, true)) return false;
            if (o->c) for (int64_t i = 0; i < n; i++) o->t2[2 * i + 1] = -o->t2[2 * i + 1];
            return spf_solve(o->F, o->t2, y, 1, false);
        }
    }
    return false;
}

/* ---- options -------------------------------------------------------------- */

typedef struct {
    bool has_k;
    int64_t k;           /* |k| */
    bool from_end;       /* -k */
    bool upto;
    bool arnoldi;        /* Method -> "Arnoldi" given */
    bool shift;
    double sr, si;
    KsWhich crit;
    const Expr* basis;   /* BasisSize value (validated against n later) */
    int64_t maxit;
    double tol;
    const Expr* start;
    double sv_tol;       /* SingularValueList Tolerance (< 0: Automatic) */
} EigOpts;

enum { P_OK = 0, P_DECLINE, P_BAD };

static bool parse_k(const Expr* e, EigOpts* o) {
    if (e->type == EXPR_INTEGER) {
        int64_t k = e->data.integer;
        if (k == 0) return false;
        o->has_k = true;
        o->from_end = k < 0;
        o->k = k < 0 ? -k : k;
        return true;
    }
    if (is_head(e, SYM_UpTo) && e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER && e->data.function.args[0]->data.integer > 0) {
        o->has_k = true;
        o->upto = true;
        o->k = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

static int bad_value(const char* head, const Expr* v, const char* opt) {
    char* s = out_str(v);
    mth_message(head, "wrgopt", "The value %s of the option %s is incorrect.", s ? s : "", opt);
    free(s);
    return P_BAD;
}

/* Method -> value (value borrowed). */
static int parse_method(const Expr* v, const char* head, EigOpts* o) {
    if (is_sym(v, "Automatic")) return P_OK;
    const Expr* name = v;
    size_t nsub = 0;
    Expr* const* sub = NULL;
    if (is_head(v, SYM_List) && v->data.function.arg_count >= 1) {
        name = v->data.function.args[0];
        sub = v->data.function.args + 1;
        nsub = v->data.function.arg_count - 1;
    }
    if (!is_str(name, "Arnoldi")) return P_DECLINE;   /* Direct, Banded, FEAST: the dense paths */
    o->arnoldi = true;
    for (size_t i = 0; i < nsub; i++) {
        const Expr* r = sub[i];
        if (!is_rule(r)) return P_DECLINE;
        const Expr* key = r->data.function.args[0];
        const Expr* val = rule_rhs(r);
        bool autov = is_sym(val, "Automatic");
        if (key_named(key, "Shift")) {
            if (autov) continue;
            double re, im;
            if (!na_read_scalar(val, &re, &im) || !isfinite(re) || !isfinite(im)) return bad_value(head, val, "Shift");
            o->shift = true;
            o->sr = re;
            o->si = im;
        } else if (key_named(key, "Tolerance")) {
            if (autov) continue;
            double re, im;
            if (!na_read_scalar(val, &re, &im) || im != 0.0 || !(re >= 0.0)) return bad_value(head, val, "Tolerance");
            o->tol = re;
        } else if (key_named(key, "MaxIterations")) {
            if (autov) continue;
            if (val->type != EXPR_INTEGER || val->data.integer < 1) return bad_value(head, val, "MaxIterations");
            o->maxit = val->data.integer;
        } else if (key_named(key, "BasisSize")) {
            if (!autov) o->basis = val;
        } else if (key_named(key, "Criteria")) {
            if (autov || is_str(val, "Magnitude")) o->crit = KS_LM;
            else if (is_str(val, "RealPart")) o->crit = KS_LR;
            else if (is_str(val, "ImaginaryPart")) o->crit = KS_LI;
            else if (is_str(val, "BothEnds")) o->crit = KS_BE;
            else {
                mth_message(head, "eslct",
                            "A value of the option Criteria must be one of \"Magnitude\", \"RealPart\", "
                            "\"ImaginaryPart\", or \"BothEnds\", with \"BothEnds\" appropriate for real "
                            "symmetric eigenvalue problems only.");
                return P_BAD;
            }
        } else if (key_named(key, "StartingVector")) {
            if (!autov) o->start = val;
        } else {
            char* ks = out_str(key);
            char* vs = out_str(v);
            mth_message(head, "moptx",
                        "Method option %s in %s is not one of {Shift, Tolerance, BasisSize, MaxIterations, "
                        "Criteria, StartingVector}.", ks ? ks : "", vs ? vs : "");
            free(ks);
            free(vs);
            return P_BAD;
        }
    }
    return P_OK;
}

/* Options after the positional arguments. */
static int parse_opts(Expr* const* opts, size_t nopts, const char* head, bool svd, EigOpts* o) {
    for (size_t i = 0; i < nopts; i++) {
        const Expr* r = opts[i];
        if (!is_rule(r)) return P_DECLINE;
        const Expr* v = rule_rhs(r);
        if (rule_named(r, "Method")) {
            int st = parse_method(v, head, o);
            if (st != P_OK) return st;
        } else if (rule_named(r, "Cubics") || rule_named(r, "Quartics")) {
            continue;              /* radicals: no bearing on machine numbers */
        } else if (svd && rule_named(r, "Tolerance")) {
            double re, im;
            if (is_sym(v, "Automatic")) continue;
            if (!na_read_scalar(v, &re, &im) || im != 0.0 || !(re >= 0.0)) return P_DECLINE;
            o->sv_tol = re;
        } else if (svd && rule_named(r, "TargetStructure")) {
            if (!is_sym(v, "Automatic")) return P_DECLINE;
        } else {
            return P_DECLINE;
        }
    }
    return P_OK;
}

static void opts_init(EigOpts* o) {
    memset(o, 0, sizeof *o);
    o->crit = KS_LM;
    o->sv_tol = -1.0;
}

/* Checks that need the dimension: BasisSize and StartingVector. */
static int check_dim_opts(const char* head, const EigOpts* o, int64_t n, int* basis) {
    *basis = 0;
    if (o->basis) {
        const Expr* b = o->basis;
        if (b->type != EXPR_INTEGER || b->data.integer <= 3 || b->data.integer > n) {
            char* s = out_str(b);
            mth_message(head, "arbss",
                        "The value of the option BasisSize -> %s should be an integer greater than 3, and "
                        "cannot exceed the dimension of the matrix, %lld.", s ? s : "", (long long)n);
            free(s);
            return P_BAD;
        }
        *basis = (int)b->data.integer;
    }
    if (o->start) {
        int len = 0;
        double* buf = NULL;
        bool ok = na_load_vector(o->start, true, &len, &buf) && len == n;
        free(buf);
        if (!ok) {
            msg1(head, "svecnm",
                 "A value of the option StartingVector -> %s must be a vector of machine numbers, of the "
                 "same length as the row dimension of the matrix.", o->start);
            return P_BAD;
        }
    }
    return P_OK;
}

static void msg_take(EigKind kind, const EigOpts* o, int64_t total) {
    const char* what = kind == EK_VALUES ? "eigenvalues"
                     : kind == EK_VECTORS ? "eigenvectors"
                     : kind == EK_SYSTEM ? "eigenvectors and eigenvalues"
                     : kind == EK_SVL ? "singular values" : "singular values and singular vectors";
    if (o->from_end)
        mth_message(kind_head(kind), "take", "Cannot take %s -%lld through -1 out of the total of %lld %s.",
                    what, (long long)o->k, (long long)total, what);
    else
        mth_message(kind_head(kind), "take", "Cannot take %s 1 through %lld out of the total of %lld %s.",
                    what, (long long)o->k, (long long)total, what);
}

/* ::arh (machine) / ::arhm (exact): the sparse array is about to be
 * densified for a count the Krylov methods would not serve well. */
static void msg_arh(EigKind kind, int64_t k, int64_t n, bool machine) {
    bool sv = kind == EK_SVL || kind == EK_SVD;
    const char* what = sv ? "singular values and/or singular vectors" : "eigenvalues and/or eigenvectors";
    if (machine)
        mth_message(kind_head(kind), "arh",
                    "Because finding %lld out of the %lld %s is likely to be faster with dense matrix methods, "
                    "the sparse input matrix will be converted. If fewer %s would be sufficient, consider "
                    "restricting this number using the second argument to %s.",
                    (long long)k, (long long)n, what, what, kind_head(kind));
    else
        mth_message(kind_head(kind), "arhm",
                    "Because finding %lld out of the %lld %s is likely to be faster with dense matrix methods, "
                    "the sparse input matrix will be converted. If fewer %s at machine precision would be "
                    "sufficient, consider using N on the matrix and restricting this number using the second "
                    "argument to %s.",
                    (long long)k, (long long)n, what, what, kind_head(kind));
}

/* Exact entries throughout (no inexact number at all). */
static bool sparse_is_exact(const SparseArray* sa) {
    if (!sa || sa->rank != 2) return false;
    if (common_scan_inexact(sa->bg).has_inexact) return false;
    if (sa->kind == SPV_INT64) return true;
    if (sa->kind != SPV_EXPR) return false;
    Expr* const* v = sa->vals;
    for (int64_t k = 0; k < sa->pat->nnz; k++)
        if (common_scan_inexact(v[k]).has_inexact) return false;
    return true;
}

/* The reference's densify notice for a request it answers densely: n > 100
 * (rows or columns, for the singular-value heads) and at most a fifth of the
 * matrix stored -- a fuller array is dense in all but name, and converting
 * it is not news. `k` is the count (the dimension for the whole spectrum). */
static void maybe_arh(EigKind kind, const SparseArray* sa, int64_t k, int64_t kdim, bool machine) {
    int64_t m = sa->dims[0], n = sa->dims[1];
    int64_t big = m > n ? m : n;
    if (big <= 100 || 5 * sa->pat->nnz > m * n) return;
    msg_arh(kind, k, kdim, machine);
}

static void msg_maxit(const char* head, int64_t lim) {
    mth_message(head, "maxit2",
                "Warning: maximum number of iterations, %lld, has been reached by the Arnoldi algorithm "
                "without convergence to the specified tolerance, but the current best computed value has "
                "been returned. You can use method options with Method -> {Arnoldi, opts} to increase the "
                "size of basis vectors, the maximum number of iterations, reduce the tolerance, or use an "
                "estimate as a shift, any of which may help.", (long long)lim);
}

/* ---- results ---------------------------------------------------------------- */

/* Unit 2-norm, then real and positive in the first component whose modulus
 * is (to 1e-8 relative) the largest -- one fixed representative of the
 * eigenvector's sign / phase class. */
static void normalise_vector(double* x, int64_t n, bool c) {
    double mx = 0.0, nr = 0.0;
    for (int64_t i = 0; i < n; i++) {
        double a = c ? hypot(x[2 * i], x[2 * i + 1]) : fabs(x[i]);
        nr += a * a;
        if (a > mx) mx = a;
    }
    if (!(mx > 0.0)) return;
    nr = sqrt(nr);
    int64_t p = 0;
    for (int64_t i = 0; i < n; i++) {
        double a = c ? hypot(x[2 * i], x[2 * i + 1]) : fabs(x[i]);
        if (a >= mx * (1.0 - 1e-8)) { p = i; break; }
    }
    if (!c) {
        double s = (x[p] < 0.0 ? -1.0 : 1.0) / nr;
        for (int64_t i = 0; i < n; i++) x[i] *= s;
        return;
    }
    double a = hypot(x[2 * p], x[2 * p + 1]);
    double ur = x[2 * p] / a / nr, ui = -x[2 * p + 1] / a / nr;   /* conj(x_p)/|x_p| / ||x|| */
    for (int64_t i = 0; i < n; i++) {
        double xr = x[2 * i], xi = x[2 * i + 1];
        x[2 * i] = xr * ur - xi * ui;
        x[2 * i + 1] = xr * ui + xi * ur;
    }
    x[2 * p + 1] = 0.0;
}

/* perm[0..k-1]: decreasing |lambda|; ties by increasing real part, then
 * decreasing imaginary part (the dense LAPACK order) */
static void sort_desc_abs(const double* re, const double* im, int k, int* perm) {
    for (int i = 0; i < k; i++) perm[i] = i;
    for (int i = 1; i < k; i++) {        /* pre-order: re ascending, im descending */
        int cur = perm[i], j = i;
        while (j > 0) {
            int pv = perm[j - 1];
            bool before = re[cur] < re[pv] || (re[cur] == re[pv] && im[cur] > im[pv]);
            if (!before) break;
            perm[j] = pv;
            j--;
        }
        perm[j] = cur;
    }
    for (int i = 1; i < k; i++) {        /* stable: |lambda| descending */
        int cur = perm[i], j = i;
        double ac = hypot(re[cur], im[cur]);
        while (j > 0 && hypot(re[perm[j - 1]], im[perm[j - 1]]) < ac) { perm[j] = perm[j - 1]; j--; }
        perm[j] = cur;
    }
}

static Expr* build_values(const double* re, const double* im, const int* perm, int k) {
    Expr** it = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
    if (!it) return NULL;
    for (int i = 0; i < k; i++) it[i] = na_scalar(re[perm[i]], im[perm[i]]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)k);
    free(it);
    return out;
}

/* k vectors of n entries (vl doubles apart: n real, 2n complex). A packed
 * matrix when every vector is real, or when the problem is complex; else
 * (a real nonsymmetric matrix with complex eigenvectors) a List of rows,
 * real ones packed and complex ones as Complex entries, as the dense path
 * gives. */
static Expr* build_vectors(const double* X, bool xc, int64_t n, const int* perm, int k, bool problem_cplx) {
    bool all_real = true;
    if (xc)
        for (int q = 0; q < k && all_real; q++) {
            const double* x = X + (size_t)2 * (size_t)n * (size_t)perm[q];
            for (int64_t i = 0; i < n; i++) if (x[2 * i + 1] != 0.0) { all_real = false; break; }
        }
    if (all_real || problem_cplx) {
        bool cplx = !all_real;
        int w = cplx ? 2 : 1;
        double* buf = malloc(sizeof(double) * (size_t)w * (size_t)n * (size_t)(k > 0 ? k : 1));
        if (!buf) return NULL;
        for (int q = 0; q < k; q++) {
            const double* x = X + (size_t)(xc ? 2 : 1) * (size_t)n * (size_t)perm[q];
            double* d = buf + (size_t)w * (size_t)n * (size_t)q;
            if (cplx) memcpy(d, x, sizeof(double) * 2 * (size_t)n);
            else for (int64_t i = 0; i < n; i++) d[i] = xc ? x[2 * i] : x[i];
        }
        Expr* out = na_build_matrix_as(buf, k, (int)n, cplx, false, NDA_HEAD_LIST);
        free(buf);
        return out;
    }
    Expr** rows = malloc(sizeof(Expr*) * (size_t)k);
    if (!rows) return NULL;
    for (int q = 0; q < k; q++) {
        const double* x = X + (size_t)2 * (size_t)n * (size_t)perm[q];
        bool real = true;
        for (int64_t i = 0; i < n; i++) if (x[2 * i + 1] != 0.0) { real = false; break; }
        if (real) {
            double* d = malloc(sizeof(double) * (size_t)n);
            if (!d) { for (int t = 0; t < q; t++) expr_free(rows[t]); free(rows); return NULL; }
            for (int64_t i = 0; i < n; i++) d[i] = x[2 * i];
            rows[q] = na_build_vector_as(d, (int)n, false, NDA_HEAD_LIST);
            free(d);
        } else {
            rows[q] = na_build_vector_as(x, (int)n, true, NDA_HEAD_LIST);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)k);
    free(rows);
    return out;
}

/* ---- the eigenvalue driver ---------------------------------------------- */

typedef struct {
    Expr* result;       /* the answer, or NULL */
    bool settled;       /* a message was given: leave the call unevaluated */
} EigOut;

/* Factor A - sigma B; on an exactly singular matrix with an implicit
 * sigma (the -k case) retry with a tiny shift. */
static SpFactor* factor_shifted(const SpMat* A, const SpMat* B, double* sr, double si, bool implicit,
                                bool dense_src, SpfStatus* st) {
    for (int attempt = 0; attempt < 3; attempt++) {
        SpMat M;
        if (!spmat_shifted(A, *sr, si, B, &M)) { *st = SPF_NOMEM; return NULL; }
        SpfWant want = (dense_src || A->n <= 128) ? SPF_WANT_DENSE : SPF_WANT_AUTO;
        SpFactor* f = spf_factor(&M, want, st);
        if (!f && *st == SPF_TOOBIG && want == SPF_WANT_DENSE) f = spf_factor(&M, SPF_WANT_AUTO, st);
        spla_view_free(&M);
        if (f) {
            double rc = spf_rcond(f);
            if (!(implicit && rc >= 0.0 && rc < 1e-14)) return f;
            spf_free(f);
            *st = SPF_SINGULAR;
        }
        if (*st != SPF_SINGULAR || !implicit) return NULL;
        /* singular at the implicit shift 0: move off it by a tiny,
         * irrational fraction of ||A||; the eigenvalues nearest 0 are
         * still sigma + 1/theta exactly */
        double na = spmat_norm_inf(A);
        if (!(na > 0.0)) na = 1.0;
        *sr = (attempt == 0 ? 1e-7 : 1e-5) * 0.6180339887498949 * na;
    }
    return NULL;
}

static EigOut eigen_krylov(EigKind kind, const SpMat* A0, const SpMat* B0, const EigOpts* o, int basis,
                           bool dense_src) {
    EigOut out = { NULL, false };
    const char* head = kind_head(kind);
    int64_t n = A0->n;
    int nev_user = (int)(o->k < n ? o->k : n);
    bool herm = spmat_hermitian(A0);
    bool herm_b = B0 ? spmat_hermitian(B0) : true;
    bool pencil_hd = B0 && herm && herm_b && spla_hermitian_posdiag(B0);
    bool use_sigma = o->shift || o->from_end;
    double sr = o->shift ? o->sr : 0.0, si = o->shift ? o->si : 0.0;

    if (o->crit == KS_BE && !(herm && !A0->cplx && !B0)) {
        mth_message(head, "arbsym", "You can only choose the Arnoldi method option Criteria -> BothEnds for "
                    "real symmetric matrices.");
        out.settled = true;
        return out;
    }
    bool c = A0->cplx || (B0 && B0->cplx) || si != 0.0;
    /* working copies (promoted to complex when the problem is) */
    SpMat A = *A0, B;
    A.owned_ci = NULL;
    A.owned_x = NULL;
    memset(&B, 0, sizeof B);
    if (B0) { B = *B0; B.owned_ci = NULL; B.owned_x = NULL; }
    if (c) {
        if (!view_to_complex(&A) || (B0 && !view_to_complex(&B))) {
            free(A.owned_x); free(B.owned_x);
            return out;
        }
    }
    int w = c ? 2 : 1;
    OpCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.c = c;
    ctx.tn = (int64_t)w * n;
    ctx.t1 = malloc(sizeof(double) * (size_t)ctx.tn);
    double* v0 = NULL;
    KsResult res;
    memset(&res, 0, sizeof res);
    double* lre = NULL;
    double* lim = NULL;
    int* perm = NULL;
    if (!ctx.t1) goto done;

    bool op_herm;
    if (use_sigma) {
        SpfStatus st;
        ctx.F = factor_shifted(&A, B0 ? &B : NULL, &sr, si, !o->shift, dense_src, &st);
        if (!ctx.F) {
            if (st == SPF_SINGULAR && o->shift) {
                Expr* sv = na_scalar(o->sr, o->si);
                msg1(head, "ssing",
                     "Failed to factor the matrix involving the shift. This may be caused by the shift, %s, "
                     "corresponding exactly to an eigenvalue. Changing the shift value slightly might help.",
                     sv);
                expr_free(sv);
                out.settled = true;
            }
            goto done;
        }
        if (B0) { ctx.mode = OPM_FINV_B; ctx.A1 = &B; }
        else ctx.mode = OPM_FINV;
        op_herm = !B0 && herm && si == 0.0;
    } else if (B0) {
        SpfStatus st;
        ctx.F = spf_factor(&B, (dense_src || n <= 128) ? SPF_WANT_DENSE : SPF_WANT_AUTO, &st);
        if (!ctx.F) goto done;                 /* singular b: the dense QZ answers */
        ctx.mode = OPM_BINV_A;
        ctx.A1 = &A;
        op_herm = false;
    } else {
        ctx.mode = OPM_A;
        ctx.A1 = &A;
        op_herm = herm;
    }

    KsProblem p;
    memset(&p, 0, sizeof p);
    p.n = n;
    p.cplx = c;
    p.hermitian = op_herm;
    p.op = op_apply;
    p.ctx = &ctx;
    p.nev = nev_user;
    if (!op_herm && !c && p.nev < n) p.nev++;      /* complete a conjugate pair at the edge */
    p.which = o->crit;
    if (op_herm && o->crit == KS_LI) p.which = KS_LM;
    if (op_herm && o->crit == KS_LR) p.which = KS_LR;
    p.ncv = basis;
    if (p.ncv > 0 && p.ncv <= p.nev) p.ncv = (int)(p.nev + 1 <= n ? p.nev + 1 : n);
    p.maxit = o->maxit;
    p.tol = o->tol;
    /* vectors always: Eigenvalues, Eigenvectors and Eigensystem then make
     * the identical Krylov run, so their values agree to the bit */
    p.want_vectors = true;
    if (o->start) {
        int len = 0;
        if (!na_load_vector(o->start, c, &len, &v0) || len != n) {
            free(v0);
            v0 = NULL;
            if (!c) goto done;
        }
        p.v0 = v0;
    }
    KsStatus ks = ks_solve(&p, &res);
    if (ks == KS_FAIL || ks == KS_NOMEM) goto done;
    if (ks == KS_MAXIT) msg_maxit(head, o->maxit > 0 ? o->maxit : 1000);

    /* back to eigenvalues of the matrix: lambda = sigma + 1/theta */
    int got = res.nconv < nev_user ? res.nconv : nev_user;
    lre = malloc(sizeof(double) * (size_t)(got > 0 ? got : 1));
    lim = malloc(sizeof(double) * (size_t)(got > 0 ? got : 1));
    perm = malloc(sizeof(int) * (size_t)(got > 0 ? got : 1));
    if (!lre || !lim || !perm) goto done;
    bool real_vals = (op_herm && !use_sigma) || (herm && !B0 && si == 0.0) || pencil_hd;
    for (int i = 0; i < got; i++) {
        double tr = res.re[i], ti = res.im[i];
        if (use_sigma) {
            double d = tr * tr + ti * ti;
            lre[i] = sr + (d > 0.0 ? tr / d : INFINITY);
            lim[i] = si + (d > 0.0 ? -ti / d : 0.0);
        } else {
            lre[i] = tr;
            lim[i] = ti;
        }
        if (real_vals) lim[i] = 0.0;
    }
    sort_desc_abs(lre, lim, got, perm);
    if (kind != EK_VALUES && res.vec) {
        for (int i = 0; i < got; i++)
            normalise_vector(res.vec + (size_t)(res.vec_cplx ? 2 : 1) * (size_t)n * (size_t)i, n, res.vec_cplx);
        if (pencil_hd && !c && res.vec_cplx) {
            /* a real definite pencil: the eigenvectors are real */
            for (int i = 0; i < got; i++) {
                double* x = res.vec + (size_t)2 * (size_t)n * (size_t)i;
                for (int64_t j = 0; j < n; j++) x[2 * j + 1] = 0.0;
            }
        }
    }
    Expr* vals = (kind != EK_VECTORS) ? build_values(lre, lim, perm, got) : NULL;
    Expr* vecs = (kind != EK_VALUES) ? build_vectors(res.vec, res.vec_cplx, n, perm, got, c) : NULL;
    if (kind == EK_VALUES) out.result = vals;
    else if (kind == EK_VECTORS) out.result = vecs;
    else if (vals && vecs) {
        Expr* pr[2] = { vals, vecs };
        out.result = expr_new_function(expr_new_symbol(SYM_List), pr, 2);
        vals = vecs = NULL;
    }
    if (kind == EK_SYSTEM) { expr_free(vals); expr_free(vecs); }

done:
    ks_result_free(&res);
    if (ctx.F) spf_free(ctx.F);
    free(ctx.t1);
    free(v0);
    free(lre);
    free(lim);
    free(perm);
    free(A.owned_x);
    free(B.owned_x);
    return out;
}

/* ---- singular values -------------------------------------------------------- */

/* Dot[ConjugateTranspose[s], s] (or s . ConjugateTranspose[s]) as a machine
 * view: the normal matrix, formed sparse by the native SpGEMM. */
static bool normal_matrix(const Expr* s, bool left, Expr** keep, SpMat* N) {
    Expr* ct_args[1] = { expr_copy((Expr*)s) };
    Expr* ct = expr_new_function(expr_new_symbol("ConjugateTranspose"), ct_args, 1);
    Expr* dargs[2];
    if (left) { dargs[0] = ct; dargs[1] = expr_copy((Expr*)s); }
    else { dargs[0] = expr_copy((Expr*)s); dargs[1] = ct; }
    Expr* call = expr_new_function(expr_new_symbol("Dot"), dargs, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    if (!r || r->type != EXPR_SPARSE || !spla_view(sp_of(r), true, N)) { expr_free(r); return false; }
    *keep = r;
    return true;
}

static Expr* sparse_diag_matrix(const double* sv, int k, int64_t rows, int64_t cols) {
    Expr** rules = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
    if (!rules) return NULL;
    int nr = 0;
    for (int i = 0; i < k; i++) {
        if (sv[i] == 0.0) continue;
        Expr* pos[2] = { expr_new_integer(i + 1), expr_new_integer(i + 1) };
        Expr* ra[2] = { expr_new_function(expr_new_symbol(SYM_List), pos, 2), expr_new_real(sv[i]) };
        rules[nr++] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    }
    Expr* dims[2] = { expr_new_integer(rows), expr_new_integer(cols) };
    Expr* args[3] = { expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nr),
                      expr_new_function(expr_new_symbol(SYM_List), dims, 2), expr_new_real(0.0) };
    free(rules);
    Expr* call = expr_new_function(expr_new_symbol("SparseArray"), args, 3);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static EigOut svd_krylov(EigKind kind, const Expr* se, const SpMat* A0, const EigOpts* o, int basis) {
    EigOut out = { NULL, false };
    const char* head = kind_head(kind);
    int64_t m = A0->m, n = A0->n;
    bool left = m >= n;                 /* N = A^H A (n x n): vectors v; else A A^H: vectors u */
    int64_t dn = left ? n : m;
    int kk = (int)(o->k < dn ? o->k : dn);
    bool c = A0->cplx;
    int w = c ? 2 : 1;
    SpMat At;                           /* A^H */
    memset(&At, 0, sizeof At);
    Expr* nkeep = NULL;
    SpMat N;
    memset(&N, 0, sizeof N);
    OpCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.c = c;
    KsResult res;
    memset(&res, 0, sizeof res);
    double* sv = NULL;
    double* U = NULL;
    double* Vv = NULL;
    int* perm = NULL;
    double sigma = 0.0;
    bool use_sigma = false;
    if (!spmat_transpose(A0, true, &At)) return out;
    int64_t big = m > n ? m : n;
    ctx.t1 = malloc(sizeof(double) * (size_t)w * (size_t)big);
    ctx.t2 = malloc(sizeof(double) * (size_t)w * (size_t)big);
    if (!ctx.t1 || !ctx.t2) goto done;

    if (!o->from_end) {
        ctx.mode = OPM_NORMAL;
        ctx.A1 = left ? A0 : &At;
        ctx.A2 = left ? &At : A0;
    } else {
        SpfStatus st = SPF_FAILED;
        if (m == n) {
            ctx.F = spf_factor(A0, n <= 128 ? SPF_WANT_DENSE : SPF_WANT_LU, &st);
            if (ctx.F && spf_rcond(ctx.F) >= 0.0 && spf_rcond(ctx.F) < 1e-14) { spf_free(ctx.F); ctx.F = NULL; }
            if (ctx.F) ctx.mode = OPM_INV_AHA;
        }
        if (!ctx.F) {
            /* the normal matrix itself, shifted off 0 when singular */
            if (!normal_matrix(se, left, &nkeep, &N)) goto done;
            double sr = 0.0;
            ctx.F = factor_shifted(&N, NULL, &sr, 0.0, true, false, &st);
            if (!ctx.F) goto done;
            sigma = sr;
            use_sigma = true;
            ctx.mode = OPM_FINV;
        }
    }
    KsProblem p;
    memset(&p, 0, sizeof p);
    p.n = dn;
    p.cplx = c;
    p.hermitian = true;
    p.op = op_apply;
    p.ctx = &ctx;
    p.nev = kk;
    p.which = o->from_end ? KS_LM : KS_LR;
    p.ncv = basis;
    if (p.ncv > 0 && p.ncv <= p.nev) p.ncv = (int)(p.nev + 1 <= dn ? p.nev + 1 : dn);
    p.maxit = o->maxit;
    p.tol = o->tol;
    p.want_vectors = true;
    KsStatus ks = ks_solve(&p, &res);
    if (ks == KS_FAIL || ks == KS_NOMEM) goto done;
    if (ks == KS_MAXIT) msg_maxit(head, o->maxit > 0 ? o->maxit : 1000);
    int got = res.nconv;
    sv = malloc(sizeof(double) * (size_t)(got > 0 ? got : 1));
    perm = malloc(sizeof(int) * (size_t)(got > 0 ? got : 1));
    U = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)(got > 0 ? got : 1));
    Vv = malloc(sizeof(double) * (size_t)w * (size_t)n * (size_t)(got > 0 ? got : 1));
    if (!sv || !perm || !U || !Vv) goto done;
    int64_t xl = (int64_t)(res.vec_cplx ? 2 : 1) * dn;
    for (int i = 0; i < got; i++) {
        const double* x = res.vec + (size_t)xl * (size_t)i;
        /* the eigenvector in the problem's arithmetic */
        double* xv = left ? Vv + (size_t)w * (size_t)n * (size_t)i : U + (size_t)w * (size_t)m * (size_t)i;
        for (int64_t j = 0; j < dn; j++) {
            if (c) { xv[2 * j] = x[2 * j]; xv[2 * j + 1] = x[2 * j + 1]; }
            else xv[j] = res.vec_cplx ? x[2 * j] : x[j];
        }
        /* the partner vector: u = A v / sigma, or v = A^H u / sigma */
        double* yv = left ? U + (size_t)w * (size_t)m * (size_t)i : Vv + (size_t)w * (size_t)n * (size_t)i;
        spla_spmv(left ? A0 : &At, xv, yv, false);
        int64_t yl = (int64_t)w * (left ? m : n);
        double s2 = 0.0;
        for (int64_t j = 0; j < yl; j++) s2 += yv[j] * yv[j];
        double s = sqrt(s2);
        if (o->from_end && !(s > 0.0)) {
            /* fall back to the eigenvalue */
            double th = res.re[i];
            double lam = use_sigma ? sigma + 1.0 / th : 1.0 / th;
            s = lam > 0.0 ? sqrt(lam) : 0.0;
        }
        sv[i] = s;
        if (s > 0.0) for (int64_t j = 0; j < yl; j++) yv[j] /= s;
        else memset(yv, 0, sizeof(double) * (size_t)yl);
    }
    /* decreasing */
    for (int i = 0; i < got; i++) perm[i] = i;
    for (int i = 1; i < got; i++) {
        int cur = perm[i], j = i;
        while (j > 0 && sv[perm[j - 1]] < sv[cur]) { perm[j] = perm[j - 1]; j--; }
        perm[j] = cur;
    }
    /* drop numerically zero singular values (and below Tolerance) */
    double smax = 0.0;
    for (int i = 0; i < got; i++) smax = fmax(smax, sv[i]);
    double floor_ = o->sv_tol >= 0.0 ? o->sv_tol * smax : (double)big * DBL_EPSILON * smax;
    /* a value under the tolerance is reported as 0. -- the count reaches
     * it, as in the dense SingularValueList */
    int keep = got;
    for (int i = 0; i < got; i++)
        if (o->sv_tol >= 0.0 ? sv[i] < floor_ : sv[i] <= floor_) sv[i] = 0.0;
    if (kind == EK_SVL) {
        Expr** it = malloc(sizeof(Expr*) * (size_t)(keep > 0 ? keep : 1));
        if (!it) goto done;
        for (int i = 0; i < keep; i++) it[i] = expr_new_real(sv[perm[i]]);
        out.result = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)keep);
        free(it);
        goto done;
    }
    /* {U (m x k), SparseArray (k x k), V (n x k)} */
    {
        double* Ub = malloc(sizeof(double) * (size_t)w * (size_t)m * (size_t)(keep > 0 ? keep : 1));
        double* Vb = malloc(sizeof(double) * (size_t)w * (size_t)n * (size_t)(keep > 0 ? keep : 1));
        double* sk = malloc(sizeof(double) * (size_t)(keep > 0 ? keep : 1));
        if (!Ub || !Vb || !sk) { free(Ub); free(Vb); free(sk); goto done; }
        for (int q = 0; q < keep; q++) {
            int i = perm[q];
            sk[q] = sv[i];
            /* fix the phase on v (the eigenvector), carry it to u */
            double* vv = Vv + (size_t)w * (size_t)n * (size_t)i;
            double* uu = U + (size_t)w * (size_t)m * (size_t)i;
            double* fixv = left ? vv : uu;
            int64_t fl = left ? n : m;
            /* phase of the largest component */
            double mx = 0.0;
            int64_t pidx = 0;
            for (int64_t j = 0; j < fl; j++) {
                double a = c ? hypot(fixv[2 * j], fixv[2 * j + 1]) : fabs(fixv[j]);
                if (a > mx * (1.0 + 1e-8)) { mx = a; pidx = j; }
            }
            double ur = 1.0, ui = 0.0;
            if (mx > 0.0) {
                if (c) { double a = hypot(fixv[2 * pidx], fixv[2 * pidx + 1]); ur = fixv[2 * pidx] / a; ui = -fixv[2 * pidx + 1] / a; }
                else ur = fixv[pidx] < 0.0 ? -1.0 : 1.0;
            }
            for (int pass = 0; pass < 2; pass++) {
                double* z = pass == 0 ? vv : uu;
                int64_t zl = pass == 0 ? n : m;
                double* d = pass == 0 ? Vb : Ub;
                for (int64_t j = 0; j < zl; j++) {
                    /* column q of an (zl x keep) column-major matrix */
                    if (c) {
                        double zr = z[2 * j], zi = z[2 * j + 1];
                        d[2 * ((size_t)j + (size_t)q * (size_t)zl)] = zr * ur - zi * ui;
                        d[2 * ((size_t)j + (size_t)q * (size_t)zl) + 1] = zr * ui + zi * ur;
                    } else {
                        d[(size_t)j + (size_t)q * (size_t)zl] = z[j] * ur;
                    }
                }
            }
        }
        Expr* Ue = na_build_matrix_as(Ub, (int)m, keep, c, true, NDA_HEAD_LIST);
        Expr* Ve = na_build_matrix_as(Vb, (int)n, keep, c, true, NDA_HEAD_LIST);
        Expr* We = sparse_diag_matrix(sk, keep, keep, keep);
        free(Ub); free(Vb); free(sk);
        if (!Ue || !Ve || !We) { expr_free(Ue); expr_free(Ve); expr_free(We); goto done; }
        Expr* tr[3] = { Ue, We, Ve };
        out.result = expr_new_function(expr_new_symbol(SYM_List), tr, 3);
    }

done:
    ks_result_free(&res);
    if (ctx.F) spf_free(ctx.F);
    free(ctx.t1);
    free(ctx.t2);
    spla_view_free(&At);
    if (nkeep) { spla_view_free(&N); expr_free(nkeep); }
    free(sv); free(U); free(Vv); free(perm);
    return out;
}

/* SingularValueDecomposition of a machine SparseArray on the dense path,
 * with the singular-value matrix handed back as a SparseArray (background
 * 0.), as the reference answers {List, SparseArray, List}. */
static Expr* svd_dense_sparse_middle(const Expr* call) {
    size_t na = call->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * na);
    if (!args) return NULL;
    Expr* nargs[1] = { expr_copy(call->data.function.args[0]) };
    Expr* ncall = expr_new_function(expr_new_symbol("Normal"), nargs, 1);
    args[0] = evaluate(ncall);
    expr_free(ncall);
    for (size_t i = 1; i < na; i++) args[i] = expr_copy(call->data.function.args[i]);
    Expr* svd = expr_new_function(expr_new_symbol("SingularValueDecomposition"), args, na);
    free(args);
    Expr* r = evaluate(svd);
    expr_free(svd);
    if (r && r->type == EXPR_NDARRAY) {
        /* a square answer {u, w, v} packs into one 3 x n x n buffer:
         * split it into its three (packed) factors, r[[1]], r[[2]], r[[3]] */
        Expr* parts[3];
        for (int i = 0; i < 3; i++) {
            Expr* pa[2] = { expr_copy(r), expr_new_integer(i + 1) };
            Expr* pc = expr_new_function(expr_new_symbol("Part"), pa, 2);
            parts[i] = evaluate(pc);
            expr_free(pc);
        }
        expr_free(r);
        r = expr_new_function(expr_new_symbol(SYM_List), parts, 3);
    }
    if (!r || !is_head(r, SYM_List) || r->data.function.arg_count != 3) { expr_free(r); return NULL; }
    Expr* mid = r->data.function.args[1];
    int rows = 0, cols = 0;
    double* buf = NULL;
    if (!na_load_matrix(mid, false, false, &rows, &cols, &buf)) { expr_free(r); return NULL; }
    int k = rows < cols ? rows : cols;
    double* d = malloc(sizeof(double) * (size_t)(k > 0 ? k : 1));
    if (!d) { free(buf); expr_free(r); return NULL; }
    for (int i = 0; i < k; i++) d[i] = buf[(size_t)i * (size_t)cols + (size_t)i];
    free(buf);
    Expr* sm = sparse_diag_matrix(d, k, rows, cols);
    free(d);
    if (!sm) { expr_free(r); return NULL; }
    expr_free(r->data.function.args[1]);
    r->data.function.args[1] = sm;
    return r;
}

/* The exact zero matrix with a count: every eigenvalue 0, and the
 * reference's eigenvectors -- the unit vectors e_n, e_(n-1), ..., e_1, of
 * which k takes the first and -k the last. Answered here because the dense
 * exact path builds a characteristic polynomial of degree n first. */
static Expr* exact_zero_eigen(EigKind kind, const EigOpts* o, int64_t n) {
    int64_t k = o->k < n ? o->k : n;
    Expr* vals = NULL;
    Expr* vecs = NULL;
    if (kind != EK_VECTORS) {
        Expr** it = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
        if (!it) return NULL;
        for (int64_t i = 0; i < k; i++) it[i] = expr_new_integer(0);
        vals = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)k);
        free(it);
    }
    if (kind != EK_VALUES) {
        Expr** rows = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
        Expr** comp = malloc(sizeof(Expr*) * (size_t)n);
        if (!rows || !comp) { free(rows); free(comp); expr_free(vals); return NULL; }
        for (int64_t q = 0; q < k; q++) {
            /* position in the full list {e_n, ..., e_1} */
            int64_t pos = o->from_end ? n - k + q : q;
            int64_t one = n - 1 - pos;                 /* 0-based index of the 1 */
            for (int64_t i = 0; i < n; i++) comp[i] = expr_new_integer(i == one ? 1 : 0);
            rows[q] = expr_new_function(expr_new_symbol(SYM_List), comp, (size_t)n);
        }
        vecs = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)k);
        free(rows);
        free(comp);
    }
    if (kind == EK_VALUES) return vals;
    if (kind == EK_VECTORS) return vecs;
    Expr* pr[2] = { vals, vecs };
    return expr_new_function(expr_new_symbol(SYM_List), pr, 2);
}

/* ---- the native handlers ------------------------------------------------- */

/* Eigenvalues / Eigenvectors / Eigensystem [s or {a, b}, k, opts]. */
static Expr* native_eigen(const Expr* call, EigKind kind) {
    size_t na = call->data.function.arg_count;
    Expr* const* a = call->data.function.args;
    if (na < 1) return NULL;
    const char* head = kind_head(kind);
    EigOpts o;
    opts_init(&o);
    if (na < 2 || is_rule(a[1])) {
        /* the whole spectrum: dense, with the reference's notice */
        const SparseArray* sa = sp_of(a[0]);
        if (sa && sa->rank == 2 && sa->dims[0] == sa->dims[1]) {
            if (sparse_is_machine(sa)) maybe_arh(kind, sa, sa->dims[0], sa->dims[0], true);
            else if (sparse_is_exact(sa)) maybe_arh(kind, sa, sa->dims[0], sa->dims[0], false);
        }
        return NULL;
    }
    if (!parse_k(a[1], &o)) return NULL;
    /* the matrix, or the pencil {a, b} */
    const Expr* ae = a[0];
    const Expr* be = NULL;
    if (is_head(ae, SYM_List) && ae->data.function.arg_count == 2 &&
        (ae->data.function.args[0]->type == EXPR_SPARSE || ae->data.function.args[1]->type == EXPR_SPARSE)) {
        be = ae->data.function.args[1];
        ae = ae->data.function.args[0];
    } else if (ae->type != EXPR_SPARSE) {
        return NULL;
    }
    bool machine_a = ae->type == EXPR_SPARSE ? sparse_is_machine(sp_of(ae)) : dense_is_machine(ae);
    bool machine_b = !be || (be->type == EXPR_SPARSE ? sparse_is_machine(sp_of(be)) : dense_is_machine(be));
    if (!machine_a && !(be && machine_b)) {
        const SparseArray* sa = sp_of(ae);
        if (!be && sa && sparse_is_exact(sa) && sa->dims[0] == sa->dims[1]) {
            int64_t n = sa->dims[0];
            if (o.from_end && sa->pat->nnz == 0 && o.k <= n)
                /* the reference's shift-invert about 0 meets the zero matrix
                 * (and it still answers, densely) */
                mth_message(head, "ssing",
                            "Failed to factor the matrix involving the shift. This may be caused by the "
                            "shift, 0, corresponding exactly to an eigenvalue. Changing the shift value "
                            "slightly might help.");
            else if (5 * o.k > n && o.k <= n)
                maybe_arh(kind, sa, o.k, n, false);
            if (sa->pat->nnz == 0 && sp_is_zero(sa->bg) && na == 2 && (o.k <= n || o.upto))
                return exact_zero_eigen(kind, &o, n);
        }
        return NULL;
    }
    int ps = parse_opts(a + 2, na - 2, head, false, &o);
    if (ps == P_DECLINE) return NULL;
    if (ps == P_BAD) return expr_copy((Expr*)call);

    SpMat A, B;
    memset(&B, 0, sizeof B);
    bool okA = ae->type == EXPR_SPARSE ? spla_view(sp_of(ae), true, &A) : spla_from_dense(ae, &A);
    if (!okA) return NULL;
    if (be) {
        bool okB = be->type == EXPR_SPARSE ? spla_view(sp_of(be), true, &B) : spla_from_dense(be, &B);
        if (!okB) { spla_view_free(&A); return NULL; }
        if (B.m != A.m || B.n != A.n) { spla_view_free(&A); spla_view_free(&B); return NULL; }
    }
    Expr* result = NULL;
    int64_t n = A.n;
    if (A.m != n) goto decline;
    if (o.k > n && !o.upto) {
        msg_take(kind, &o, n);
        result = expr_copy((Expr*)call);
        goto finish;
    }
    if (o.upto && o.k > n) o.k = n;
    int basis = 0;
    if (check_dim_opts(head, &o, n, &basis) == P_BAD) { result = expr_copy((Expr*)call); goto finish; }
    if (!o.arnoldi) {
        /* the reference's choice of the dense method */
        if (be) { if (n * n <= EIG_DENSE_MAX) goto decline; }
        else if (n <= 14 || 5 * A.nnz > n * n || 5 * o.k > n) {
            if (5 * o.k > n && ae->type == EXPR_SPARSE) maybe_arh(kind, sp_of(ae), o.k, n, true);
            goto decline;
        }
    }
    {
        EigOut r = eigen_krylov(kind, &A, be ? &B : NULL, &o, basis, false);
        if (r.result) result = r.result;
        else if (r.settled) result = expr_copy((Expr*)call);
    }
    goto finish;
decline:
    result = NULL;
finish:
    spla_view_free(&A);
    if (be) spla_view_free(&B);
    return result;
}

static Expr* native_eigenvalues(const Expr* call) { return native_eigen(call, EK_VALUES); }
static Expr* native_eigenvectors(const Expr* call) { return native_eigen(call, EK_VECTORS); }
static Expr* native_eigensystem(const Expr* call) { return native_eigen(call, EK_SYSTEM); }

/* SingularValueList / SingularValueDecomposition [s, k, opts], and
 * SingularValueDecomposition[s] (dense, with a sparse middle factor). */
static Expr* native_singular(const Expr* call, EigKind kind) {
    size_t na = call->data.function.arg_count;
    Expr* const* a = call->data.function.args;
    if (na < 1 || a[0]->type != EXPR_SPARSE) return NULL;
    if (!sparse_is_machine(sp_of(a[0]))) {
        const SparseArray* sa = sp_of(a[0]);
        if (sa->rank == 2 && sparse_is_exact(sa) && (na < 2 || is_rule(a[1]))) {
            int64_t kmin = sa->dims[0] < sa->dims[1] ? sa->dims[0] : sa->dims[1];
            maybe_arh(kind, sa, kmin, kmin, false);
        }
        /* the exact zero matrix with a count: k exact zeros (reference) */
        EigOpts oz;
        opts_init(&oz);
        if (kind == EK_SVL && na == 2 && sa->rank == 2 && sparse_is_exact(sa) && sa->pat->nnz == 0 &&
            sp_is_zero(sa->bg) && parse_k(a[1], &oz)) {
            int64_t kmin = sa->dims[0] < sa->dims[1] ? sa->dims[0] : sa->dims[1];
            if (oz.k > kmin && !oz.upto) return NULL;          /* ::take: the dense path */
            int64_t k = oz.k < kmin ? oz.k : kmin;
            Expr** it = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
            if (!it) return NULL;
            for (int64_t i = 0; i < k; i++) it[i] = expr_new_integer(0);
            Expr* r = expr_new_function(expr_new_symbol(SYM_List), it, (size_t)k);
            free(it);
            return r;
        }
        return NULL;
    }
    const char* head = kind_head(kind);
    EigOpts o;
    opts_init(&o);
    size_t first_opt = 1;
    if (na >= 2 && !is_rule(a[1])) {
        if (!parse_k(a[1], &o)) return NULL;
        first_opt = 2;
    }
    int ps = parse_opts(a + first_opt, na - first_opt, head, true, &o);
    if (ps == P_DECLINE) return NULL;
    if (ps == P_BAD) return expr_copy((Expr*)call);
    if (!o.has_k) {
        const SparseArray* sa = sp_of(a[0]);
        int64_t kmin = sa->dims[0] < sa->dims[1] ? sa->dims[0] : sa->dims[1];
        maybe_arh(kind, sa, kmin, kmin, true);
        return kind == EK_SVD ? svd_dense_sparse_middle(call) : NULL;
    }
    SpMat A;
    if (!spla_view(sp_of(a[0]), true, &A)) return NULL;
    Expr* result = NULL;
    int64_t kmin = A.m < A.n ? A.m : A.n;
    if (o.k > kmin && !o.upto) {
        msg_take(kind, &o, kmin);
        spla_view_free(&A);
        return expr_copy((Expr*)call);
    }
    if (o.upto && o.k > kmin) o.k = kmin;
    int basis = 0;
    if (check_dim_opts(head, &o, A.m >= A.n ? A.n : A.m, &basis) == P_BAD) {
        spla_view_free(&A);
        return expr_copy((Expr*)call);
    }
    bool dense = !o.arnoldi && (kmin <= 14 || 5 * A.nnz > A.m * A.n || 5 * o.k > kmin);
    if (dense) {
        if (!o.arnoldi && 5 * o.k > kmin) maybe_arh(kind, sp_of(a[0]), o.k, kmin, true);
        spla_view_free(&A);
        return kind == EK_SVD ? svd_dense_sparse_middle(call) : NULL;
    }
    EigOut r = svd_krylov(kind, a[0], &A, &o, basis);
    spla_view_free(&A);
    if (r.result) result = r.result;
    else if (r.settled) result = expr_copy((Expr*)call);
    else if (kind == EK_SVD) result = svd_dense_sparse_middle(call);
    return result;
}

static Expr* native_svl(const Expr* call) { return native_singular(call, EK_SVL); }
static Expr* native_svd(const Expr* call) { return native_singular(call, EK_SVD); }

/* ---- the dense Method -> "Arnoldi" ---------------------------------------- */

Expr* spla_eigen_dense_arnoldi(Expr* m, Expr* a, int64_t n, int want, Expr* k_spec, Expr* method_value) {
    if (!k_spec || n < 1) return NULL;                 /* the whole spectrum: Direct */
    EigKind kind = want == 3 ? EK_SYSTEM : want == 2 ? EK_VECTORS : EK_VALUES;
    EigOpts o;
    opts_init(&o);
    if (!parse_k(k_spec, &o)) return NULL;
    if (o.k > n) o.k = n;
    if (o.k >= n - 1) return NULL;                    /* (nearly) everything: Direct */
    if (method_value && parse_method(method_value, kind_head(kind), &o) != P_OK) return NULL;
    o.arnoldi = true;
    SpMat A, B;
    memset(&B, 0, sizeof B);
    if (!spla_from_dense(m, &A)) return NULL;
    if (A.m != A.n) { spla_view_free(&A); return NULL; }
    if (a) {
        if (!spla_from_dense(a, &B) || B.m != A.m || B.n != A.n) {
            spla_view_free(&A);
            spla_view_free(&B);
            return NULL;
        }
    }
    int basis = 0;
    Expr* result = NULL;
    if (check_dim_opts(kind_head(kind), &o, n, &basis) == P_OK) {
        EigOut r = eigen_krylov(kind, &A, a ? &B : NULL, &o, basis, true);
        result = r.result;
    }
    spla_view_free(&A);
    if (a) spla_view_free(&B);
    return result;
}

void sparse_la_eigen_init(void) {
    symtab_set_sparse_native("Eigenvalues", (const void*)native_eigenvalues);
    symtab_set_sparse_native("Eigenvectors", (const void*)native_eigenvectors);
    symtab_set_sparse_native("Eigensystem", (const void*)native_eigensystem);
    symtab_set_sparse_native("SingularValueList", (const void*)native_svl);
    symtab_set_sparse_native("SingularValueDecomposition", (const void*)native_svd);
    /* the pencil {a, b}: arrays one List level down reach the gate too */
    symtab_set_sparse_mode("Eigenvalues", SPG_DENSE | SPG_DEEP);
    symtab_set_sparse_mode("Eigenvectors", SPG_DENSE | SPG_DEEP);
    symtab_set_sparse_mode("Eigensystem", SPG_DENSE | SPG_DEEP);
}
