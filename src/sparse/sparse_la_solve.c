/* Mathilda -- sparse linear algebra: native LinearSolve, Det and Inverse.
 *
 * These handlers run before the transparency gate (symtab_set_sparse_native)
 * and take the MACHINE case: a square SparseArray of reals or complexes with a
 * zero background. Everything else -- exact and MPFR entries, a nonzero
 * background, options they do not know -- is declined (NULL) and the gate
 * answers on the dense form, which is what the reference does too: exact and
 * arbitrary-precision sparse systems behave like dense ones.
 *
 * Reference semantics for a machine SparseArray that differ from the dense
 * matrix (probed; see docs/spec/builtins/linear-algebra.md):
 *   - a singular matrix: LinearSolve::sing, and the call stays unevaluated
 *     even when b is consistent (dense: a particular solution);
 *   - a rectangular matrix: LinearSolve::sqmat, unevaluated;
 *   - Method -> "Multifrontal" (LU), "Cholesky" (::npdef when not positive
 *     definite), "Krylov" with Method -> "ConjugateGradient" | "BiCGSTAB" |
 *     "GMRES" and "Preconditioner" -> "ILU0" | "ILUT" | "ILUTP" | None | f,
 *     "MaxIterations", Tolerance; ::kryme / ::krypre for unknown names and
 *     ::krymit when the iteration limit is reached (the last iterate is
 *     returned).
 * The answer is a dense (packed) vector or matrix, as in the reference.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "numarray.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "pack.h"
#include "linsolvefn.h"
#include "common.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Up to this order a machine sparse matrix is factored densely with LAPACK:
 * the answer then agrees bit-for-bit with the dense matrix's, and a dense LU
 * of a small matrix is faster than any sparse one. */
#define SPLA_DENSE_N 128

/* ---- small helpers -------------------------------------------------------- */

static bool is_str(const Expr* e, const char* s) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, s) == 0;
}
static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}
static const Expr* rule_lhs(const Expr* r) { return r->data.function.args[0]; }
static const Expr* rule_rhs(const Expr* r) { return r->data.function.args[1]; }
static bool rule_named(const Expr* r, const char* name) {
    const Expr* l = rule_lhs(r);
    return is_sym(l, name) || is_str(l, name);
}

static void msg_with_expr(const char* head, const char* tag, const char* fmt, const Expr* e) {
    char* s = expr_to_string((Expr*)e);
    mth_message(head, tag, fmt, s ? s : "");
    free(s);
}

/* Promote a real view to complex (owned interleaved copy). */
static bool view_to_complex(SpMat* A) {
    if (A->cplx) return true;
    double* x = malloc(sizeof(double) * 2 * (size_t)(A->nnz > 0 ? A->nnz : 1));
    if (!x) return false;
    for (int64_t k = 0; k < A->nnz; k++) { x[2 * k] = A->x[k]; x[2 * k + 1] = 0.0; }
    free(A->owned_x);
    A->owned_x = x;
    A->x = x;
    A->cplx = true;
    return true;
}

/* The numeric class of an operand: what decides which LinearSolve semantics
 * apply. A machine number anywhere (or a machine buffer) makes the system a
 * machine one -- the reference's contagion: an integer SparseArray with a
 * machine right-hand side is solved, and fails, like a machine matrix. */
typedef enum { NC_EXACT = 0, NC_MACHINE, NC_BIGNUM } NumClass;

static NumClass num_class(const Expr* e) {
    if (e->type == EXPR_NDARRAY) return NC_MACHINE;
    if (e->type == EXPR_SPARSE) {
        const SparseArray* sa = sp_of(e);
        if (sa->kind == SPV_F64 || sa->kind == SPV_C128) return NC_MACHINE;
        CommonInexactInfo bi = common_scan_inexact(sa->bg);
        long bits = bi.has_inexact ? bi.min_bits : 0;
        if (sa->kind == SPV_EXPR) {
            Expr* const* v = sa->vals;
            for (int64_t k = 0; k < sa->pat->nnz; k++) {
                CommonInexactInfo vi = common_scan_inexact(v[k]);
                if (vi.has_inexact && (bits == 0 || vi.min_bits < bits)) bits = vi.min_bits;
            }
        }
        return bits == 0 ? NC_EXACT : bits <= 53 ? NC_MACHINE : NC_BIGNUM;
    }
    CommonInexactInfo inf = common_scan_inexact(e);
    if (!inf.has_inexact) return NC_EXACT;
    return inf.min_bits <= 53 ? NC_MACHINE : NC_BIGNUM;
}

/* A machine view of the coefficient operand: a SparseArray (zero background;
 * machine values, or integers when `promote` -- the right-hand side is
 * machine) or a dense machine matrix. */
static bool coeff_view(const Expr* A, bool promote, SpMat* V) {
    if (A->type == EXPR_SPARSE) return spla_view(sp_of(A), promote, V);
    return spla_from_dense(A, V);
}

/* LinearSolve[a, b] (no options), evaluated; a and b borrowed. */
static Expr* plain_solve(const Expr* a, const Expr* b) {
    Expr* args[2] = { expr_copy((Expr*)a), expr_copy((Expr*)b) };
    Expr* call = expr_new_function(expr_new_symbol("LinearSolve"), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* N[e] (machine precision), evaluated. */
static Expr* machine_n(const Expr* e) {
    Expr* args[1] = { expr_copy((Expr*)e) };
    Expr* call = expr_new_function(expr_new_symbol("N"), args, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ---- right-hand sides ---------------------------------------------------- */

typedef struct {
    double* buf;     /* column blocks: nrhs columns of n entries (w doubles each) */
    int64_t n, nrhs;
    bool matrix;     /* b was a matrix (answer a matrix) */
    bool cplx;
} Rhs;

static void rhs_free(Rhs* r) { free(r->buf); memset(r, 0, sizeof *r); }

/* Load b (dense or sparse vector / matrix of machine-convertible numbers). */
static bool rhs_load(const Expr* b, bool want_cplx, Rhs* r) {
    memset(r, 0, sizeof *r);
    Expr* dense = NULL;
    const Expr* src = b;
    if (b->type == EXPR_SPARSE) {
        const SparseArray* sb = sp_of(b);
        if (sb->rank < 1 || sb->rank > 2) return false;
        dense = sp_normal(sb);
        if (!dense) return false;
        src = dense;
    }
    int n = 0, rows = 0, cols = 0;
    bool ok = false;
    for (int pass = want_cplx ? 1 : 0; pass < 2 && !ok; pass++) {
        bool cplx = pass == 1;
        if (na_load_vector(src, cplx, &n, &r->buf)) {
            r->n = n; r->nrhs = 1; r->matrix = false; r->cplx = cplx;
            ok = true;
        } else if (na_load_matrix(src, cplx, true, &rows, &cols, &r->buf)) {
            r->n = rows; r->nrhs = cols; r->matrix = true; r->cplx = cplx;
            ok = true;
        }
    }
    expr_free(dense);
    return ok;
}

static Expr* rhs_result(const double* X, const Rhs* r, bool cplx, NDPresentation pres) {
    if (!r->matrix) return na_build_vector_as(X, (int)r->n, cplx, pres);
    return na_build_matrix_as(X, (int)r->n, (int)r->nrhs, cplx, true, pres);
}

/* ---- options ------------------------------------------------------------- */

typedef enum { MS_AUTO, MS_MULTIFRONTAL, MS_CHOLESKY, MS_KRYLOV, MS_OTHER } MethodSel;

typedef struct {
    MethodSel sel;
    KryOpts kry;
    const Expr* start;        /* "StartingVector" */
    const Expr* user_prec;    /* "Preconditioner" -> f */
    bool bad;                 /* a message was given; leave the call */
} SolveOpts;

/* Parse "Krylov" sub-options (the rules after "Krylov" in the Method list). */
static bool parse_krylov(const Expr* const* sub, size_t nsub, SolveOpts* o) {
    for (size_t i = 0; i < nsub; i++) {
        const Expr* r = sub[i];
        if (!is_rule(r)) return false;
        const Expr* v = rule_rhs(r);
        if (rule_named(r, "Method")) {
            if (is_str(v, "ConjugateGradient")) o->kry.method = KRY_CG;
            else if (is_str(v, "BiCGSTAB")) o->kry.method = KRY_BICGSTAB;
            else if (is_str(v, "GMRES")) o->kry.method = KRY_GMRES;
            else if (is_sym(v, "Automatic")) o->kry.method = KRY_AUTO;
            else {
                if (v->type == EXPR_STRING) mth_message("LinearSolve", "kryme", "%s is not a supported Krylov method.", v->data.string);
                else msg_with_expr("LinearSolve", "kryme", "%s is not a supported Krylov method.", v);
                o->bad = true;
                return true;
            }
        } else if (rule_named(r, "Preconditioner")) {
            const Expr* name = v;
            const Expr* const* popts = NULL;
            size_t npopts = 0;
            if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
                v->data.function.head->data.symbol.name == SYM_List && v->data.function.arg_count >= 1) {
                name = v->data.function.args[0];
                popts = (const Expr* const*)v->data.function.args + 1;
                npopts = v->data.function.arg_count - 1;
            }
            if (is_str(name, "ILU0")) o->kry.prec = KPC_ILU0;
            else if (is_str(name, "ILUT")) o->kry.prec = KPC_ILUT;
            else if (is_str(name, "ILUTP")) o->kry.prec = KPC_ILUTP;
            else if (is_sym(name, "None")) o->kry.prec = KPC_NONE;
            else if (is_sym(name, "Automatic")) o->kry.prec = KPC_AUTO;
            else if (is_str(name, "Jacobi") || name->type == EXPR_STRING) {
                mth_message("LinearSolve", "krypre",
                            "Value of option Preconditioner -> %s is not ILU0, ILUT, ILUTP, or a user-defined function.",
                            name->data.string);
                o->bad = true;
                return true;
            } else {
                o->kry.prec = KPC_USER;
                o->user_prec = name;
            }
            for (size_t k = 0; k < npopts; k++) {
                const Expr* pr = popts[k];
                if (!is_rule(pr)) continue;
                const Expr* pv = rule_rhs(pr);
                double d;
                if ((rule_named(pr, "Tolerance") || rule_named(pr, "DropTolerance")) && na_read_scalar(pv, &d, NULL))
                    o->kry.ilut_droptol = d;
                else if ((rule_named(pr, "FillIn") || rule_named(pr, "Fill")) && pv->type == EXPR_INTEGER)
                    o->kry.ilut_fill = (int)pv->data.integer;
            }
        } else if (rule_named(r, "MaxIterations")) {
            if (v->type == EXPR_INTEGER && v->data.integer > 0) o->kry.max_iter = v->data.integer;
            else if (!is_sym(v, "Automatic")) return false;
        } else if (rule_named(r, "Tolerance")) {
            double d, im;
            if (na_read_scalar(v, &d, &im) && d > 0.0) o->kry.tol = d;
            else if (!is_sym(v, "Automatic")) return false;
        } else if (rule_named(r, "StartingVector")) {
            o->start = v;
        } else if (rule_named(r, "ResidualNormFunction")) {
            /* the 2-norm is the only one the iterations test against */
            if (!is_sym(v, "Automatic") && !(v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Norm))
                return false;
        } else if (rule_named(r, "BasisSize") || rule_named(r, "Restart")) {
            if (v->type == EXPR_INTEGER && v->data.integer > 0) o->kry.restart = (int)v->data.integer;
        } else {
            return false;
        }
    }
    return true;
}

/* False: an option this module does not implement (decline). */
static bool parse_opts(const Expr* const* opts, size_t nopts, SolveOpts* o) {
    memset(o, 0, sizeof *o);
    kry_defaults(&o->kry);
    o->sel = MS_AUTO;
    for (size_t i = 0; i < nopts; i++) {
        const Expr* r = opts[i];
        if (!is_rule(r)) return false;
        if (!rule_named(r, "Method")) {
            /* Modulus -> 0 and ZeroTest are no-ops for machine numbers */
            if (rule_named(r, "Modulus") && rule_rhs(r)->type == EXPR_INTEGER && rule_rhs(r)->data.integer == 0)
                continue;
            if (rule_named(r, "ZeroTest")) continue;
            return false;
        }
        const Expr* v = rule_rhs(r);
        const Expr* name = v;
        const Expr* const* sub = NULL;
        size_t nsub = 0;
        if (v->type == EXPR_FUNCTION && v->data.function.head->type == EXPR_SYMBOL &&
            v->data.function.head->data.symbol.name == SYM_List && v->data.function.arg_count >= 1) {
            name = v->data.function.args[0];
            sub = (const Expr* const*)v->data.function.args + 1;
            nsub = v->data.function.arg_count - 1;
        }
        if (is_sym(name, "Automatic")) o->sel = MS_AUTO;
        else if (is_str(name, "Multifrontal")) o->sel = MS_MULTIFRONTAL;
        else if (is_str(name, "Cholesky")) o->sel = MS_CHOLESKY;
        else if (is_str(name, "Krylov")) {
            o->sel = MS_KRYLOV;
            if (!parse_krylov(sub, nsub, o)) return false;
        } else {
            o->sel = MS_OTHER;        /* Banded, row-reduction methods, ... */
            return false;
        }
    }
    return true;
}

/* ---- a user preconditioner ---------------------------------------------- */

typedef struct { const Expr* f; int64_t n; bool cplx; } UserPrecCtx;

static bool user_prec_cb(void* c, const double* r, double* z) {
    const UserPrecCtx* u = c;
    Expr* v = na_build_vector_as(r, (int)u->n, u->cplx, NDA_HEAD_LIST);
    if (!v) return false;
    Expr* args[1] = { v };
    Expr* call = expr_new_function(expr_copy((Expr*)u->f), args, 1);
    Expr* res = evaluate(call);
    expr_free(call);
    int n = 0;
    double* buf = NULL;
    bool ok = na_load_vector(res, u->cplx, &n, &buf) && n == u->n;
    if (ok) memcpy(z, buf, sizeof(double) * (u->cplx ? 2 : 1) * (size_t)n);
    free(buf);
    expr_free(res);
    return ok;
}

/* ---- LinearSolve --------------------------------------------------------- */

Expr* spla_linearsolve_call(const Expr* call, const Expr* Ae, const Expr* Be,
                            const Expr* const* opts, size_t nopts) {
    SolveOpts o;
    if (!parse_opts(opts, nopts, &o)) return NULL;
    if (o.bad) return expr_copy((Expr*)call);

    NumClass ca = num_class(Ae), cb = num_class(Be);
    bool machine = ca == NC_MACHINE || cb == NC_MACHINE;
    if (!machine) {
        /* exact or arbitrary-precision system (reference behaviour):
         * "Multifrontal" refuses (::exopt1); "Cholesky" numericalises;
         * "Krylov" warns (::kryinfp, exact input) and the direct solve
         * answers; Automatic is the dense path's business */
        switch (o.sel) {
            case MS_MULTIFRONTAL:
                mth_message("LinearSolve", "exopt1",
                            "The option setting Method -> Multifrontal cannot be used with "
                            "arbitrary-precision or exact arguments.");
                return expr_copy((Expr*)call);
            case MS_CHOLESKY: {
                if (ca == NC_BIGNUM || cb == NC_BIGNUM) return plain_solve(Ae, Be);
                Expr* an = machine_n(Ae);
                Expr* bn = machine_n(Be);
                Expr* r = (an && bn && num_class(an) == NC_MACHINE)
                        ? spla_linearsolve_call(call, an, bn, opts, nopts) : NULL;
                if (r && expr_eq(r, call)) { expr_free(an); expr_free(bn); return r; }
                if (!r) r = plain_solve(an, bn);
                expr_free(an);
                expr_free(bn);
                return r;
            }
            case MS_KRYLOV:
                if (ca == NC_EXACT && cb == NC_EXACT)
                    mth_message("LinearSolve", "kryinfp",
                                "Input to Krylov algorithm is of infinite precision. Iterations may "
                                "become very slow due to expression size growth.");
                return plain_solve(Ae, Be);
            default:
                return NULL;
        }
    }

    SpMat A;
    if (!coeff_view(Ae, true, &A)) return NULL;
    if (A.m != A.n) {
        if (Ae->type != EXPR_SPARSE) { spla_view_free(&A); return NULL; }   /* dense: particular solution */
        spla_view_free(&A);
        msg_with_expr("LinearSolve", "sqmat",
                      "The matrix %s is not square. A square matrix is needed to compute a factorization.", Ae);
        return expr_copy((Expr*)call);
    }
    Rhs b;
    if (!rhs_load(Be, A.cplx, &b) || b.n != A.n) {
        if (b.buf) rhs_free(&b);
        spla_view_free(&A);
        return NULL;                         /* the dense path reports the shape */
    }
    bool cplx = A.cplx || b.cplx;
    if (cplx && !A.cplx && !view_to_complex(&A)) { rhs_free(&b); spla_view_free(&A); return NULL; }
    if (cplx && !b.cplx) {
        /* reload b as complex */
        rhs_free(&b);
        if (!rhs_load(Be, true, &b)) { spla_view_free(&A); return NULL; }
    }
    int w = cplx ? 2 : 1;
    size_t len = (size_t)w * (size_t)A.n;
    double* X = malloc(sizeof(double) * len * (size_t)(b.nrhs > 0 ? b.nrhs : 1));
    if (!X) { rhs_free(&b); spla_view_free(&A); return NULL; }
    Expr* result = NULL;

    if (o.sel == MS_KRYLOV) {
        UserPrecCtx uc = { o.user_prec, A.n, cplx };
        if (o.kry.prec == KPC_USER) { o.kry.user = user_prec_cb; o.kry.user_ctx = &uc; }
        double* x0 = NULL;
        if (o.start) {
            int n0 = 0;
            if (!na_load_vector(o.start, cplx, &n0, &x0) || n0 != A.n) {
                free(x0);
                x0 = NULL;
            }
        }
        bool any_maxit = false, fail = false;
        int64_t used = 0;
        for (int64_t r = 0; r < b.nrhs && !fail; r++) {
            double* xr = X + (size_t)r * len;
            if (x0) memcpy(xr, x0, sizeof(double) * len);
            else memset(xr, 0, sizeof(double) * len);
            int64_t it = 0;
            KryStatus st = kry_solve(&A, b.buf + (size_t)r * len, xr, &o.kry, &it);
            if (it > used) used = it;
            if (st == KRY_MAXIT) any_maxit = true;
            else if (st != KRY_CONVERGED) fail = true;
        }
        free(x0);
        if (fail) {
            free(X); rhs_free(&b); spla_view_free(&A);
            return NULL;     /* breakdown or preconditioner failure: the direct solve answers */
        }
        if (any_maxit) {
            int64_t lim = o.kry.max_iter > 0 ? o.kry.max_iter : used;
            mth_message("LinearSolve", "krymit",
                        "The maximum number of iterations, %lld, has been reached by the Krylov algorithm "
                        "without convergence to a solution within the required tolerance. The currently "
                        "computed solution has been returned.", (long long)lim);
        }
        result = rhs_result(X, &b, cplx, na_result_presentation(Be));
    } else {
        SpfWant want = o.sel == MS_CHOLESKY ? SPF_WANT_CHOL
                     : o.sel == MS_MULTIFRONTAL ? SPF_WANT_LU
                     : (A.n <= SPLA_DENSE_N ? SPF_WANT_DENSE : SPF_WANT_AUTO);
        SpfStatus st;
        SpFactor* f = spf_factor(&A, want, &st);
        if (!f && st == SPF_TOOBIG && o.sel == MS_AUTO) {
            /* no sparse direct solver in this build and the matrix is too
             * large to densify: the Krylov default */
            free(X); rhs_free(&b); spla_view_free(&A);
            Expr* ka[2] = { expr_new_symbol(SYM_Method), expr_new_string("Krylov") };
            Expr* kopt = expr_new_function(expr_new_symbol(SYM_Rule), ka, 2);
            const Expr* kops[1] = { kopt };
            Expr* r = spla_linearsolve_call(call, Ae, Be, kops, 1);
            expr_free(kopt);
            return r;
        }
        if (!f) {
            free(X); rhs_free(&b); spla_view_free(&A);
            if (st == SPF_SINGULAR) {
                if (Ae->type != EXPR_SPARSE) return NULL;     /* dense: its own semantics */
                msg_with_expr("LinearSolve", "sing", "Matrix %s is singular.", Ae);
                return expr_copy((Expr*)call);
            }
            if (st == SPF_NOTPD) {
                msg_with_expr("LinearSolve", "npdef", "The matrix %s is not positive definite.", Ae);
                return expr_copy((Expr*)call);
            }
            return NULL;
        }
        bool ok = spf_solve(f, b.buf, X, b.nrhs, false);
        double rc = spf_rcond(f);
        spf_free(f);
        if (!ok) { free(X); rhs_free(&b); spla_view_free(&A); return NULL; }
        if (rc >= 0.0 && rc < 1e-16)
            msg_with_expr("LinearSolve", "luc",
                          "Result for LinearSolve of badly conditioned matrix %s may contain significant numerical errors.",
                          Ae);
        result = rhs_result(X, &b, cplx, na_result_presentation(Be));
    }
    free(X);
    rhs_free(&b);
    spla_view_free(&A);
    return result;
}

static Expr* native_linearsolve(const Expr* call) {
    size_t n = call->data.function.arg_count;
    Expr* const* a = call->data.function.args;
    if (n < 1 || a[0]->type != EXPR_SPARSE) return NULL;
    if (n == 1 || (n >= 2 && is_rule(a[1]))) {
        /* LinearSolve[s, opts...]: a LinearSolveFunction */
        return linsolvefn_make(call);
    }
    for (size_t i = 2; i < n; i++)
        if (!is_rule(a[i])) return NULL;
    return spla_linearsolve_call(call, a[0], a[1], (const Expr* const*)a + 2, n - 2);
}

/* ---- Det ----------------------------------------------------------------- */

static Expr* native_det(const Expr* call) {
    if (call->data.function.arg_count != 1) return NULL;
    const SparseArray* sa = sp_of(call->data.function.args[0]);
    SpMat A;
    if (!sa || !spla_view(sa, false, &A)) return NULL;
    if (A.m != A.n || A.n <= SPLA_DENSE_N) { spla_view_free(&A); return NULL; }
    SpfStatus st;
    SpFactor* f = spf_factor(&A, SPF_WANT_LU, &st);
    bool cplx = A.cplx;
    spla_view_free(&A);
    if (!f) {
        if (st == SPF_SINGULAR) return cplx ? na_scalar(0.0, 0.0) : expr_new_real(0.0);
        return NULL;
    }
    double re, im;
    bool ok = spf_det(f, &re, &im);
    spf_free(f);
    if (!ok) return NULL;           /* overflow: the dense path answers in arbitrary precision */
    return cplx ? na_scalar(re, im) : expr_new_real(re);
}

/* ---- Inverse ------------------------------------------------------------- */

static Expr* native_inverse(const Expr* call) {
    if (call->data.function.arg_count != 1) return NULL;
    const Expr* Ae = call->data.function.args[0];
    const SparseArray* sa = sp_of(Ae);
    SpMat A;
    if (!sa || !spla_view(sa, false, &A)) return NULL;
    if (A.m != A.n || A.n <= SPLA_DENSE_N || A.n > ((int64_t)1 << 27) / A.n) {
        spla_view_free(&A);
        return NULL;
    }
    SpfStatus st;
    SpFactor* f = spf_factor(&A, SPF_WANT_AUTO, &st);
    int64_t n = A.n;
    bool cplx = A.cplx;
    spla_view_free(&A);
    if (!f) {
        if (st == SPF_SINGULAR) {
            msg_with_expr("Inverse", "sing", "Matrix %s is singular.", Ae);
            return expr_copy((Expr*)call);
        }
        return NULL;
    }
    int w = cplx ? 2 : 1;
    double* E = calloc((size_t)w * (size_t)n * (size_t)n, sizeof(double));
    if (!E) { spf_free(f); return NULL; }
    for (int64_t j = 0; j < n; j++) E[(size_t)w * ((size_t)j * (size_t)n + (size_t)j)] = 1.0;
    bool ok = spf_solve(f, E, E, n, false);
    spf_free(f);
    Expr* r = ok ? na_build_matrix_as(E, (int)n, (int)n, cplx, true, NDA_HEAD_LIST) : NULL;
    free(E);
    return r;
}

void sparse_la_init(void) {
    symtab_set_sparse_native("LinearSolve", (const void*)native_linearsolve);
    symtab_set_sparse_native("Det", (const void*)native_det);
    symtab_set_sparse_native("Inverse", (const void*)native_inverse);
}
