/*
 * src/linalg/orthogonalize.c -- Orthogonalize.
 *
 *   Orthogonalize[{v1, v2, ...}]       Gram-Schmidt orthonormalisation under
 *                                       the inner product u . Conjugate[v]
 *   Orthogonalize[{e1, e2, ...}, f]     under the inner product f[u, v] (the
 *                                       elements need not be vectors)
 *   Method -> Automatic | "GramSchmidt" | "ModifiedGramSchmidt" |
 *             "Householder" | "Reorthogonalization"
 *   Tolerance -> Automatic | t           a residual whose norm is below
 *                                       t * Norm[v_k] becomes a zero vector
 *
 * Reference semantics (wolframscript): a vector linearly dependent on the
 * earlier ones becomes a zero vector (machine: 0.); there are as many outputs as
 * inputs.  Symbolic and exact input give exact answers in the classical form
 *     w_k = v_k - Sum_j (<v_k, w_j> / <w_j, w_j>) w_j,   e_k = w_k / Norm[w_k]
 * (so Orthogonalize[{{a,b},{c,d}}] carries Abs[a]^2 + Abs[b]^2 and
 * c Conjugate[a] + d Conjugate[b]).  The Method only shapes a machine
 * computation: exact, symbolic and arbitrary-precision input always take the
 * classical form above, as the reference does (Householder on an exact matrix
 * gives the Gram-Schmidt answer).
 *
 * Machine numbers (default inner product): orth_machine over raw doubles --
 * classical GS, modified GS, MGS with one full reorthogonalisation pass
 * (Automatic and "Reorthogonalization"), or Householder via LAPACK geqrf/orgqr.
 * The Automatic kernel is the ONE helper behind List, packed, NDArray and
 * Compile[], so they agree bit for bit.
 */

#include "la_dense.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "lapack.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "numeric.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"
#include "parse.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define ORTH_HEAD "Orthogonalize"

typedef enum {
    ORTH_AUTO = 0, ORTH_GS, ORTH_MGS, ORTH_HOUSEHOLDER, ORTH_REORTH
} OrthMethod;

typedef enum { OM_OK = 0, OM_DECLINE, OM_FAIL } OrthStatus;

/* ------------------------------------------------------------------------ *
 *  Machine kernel                                                          *
 * ------------------------------------------------------------------------ */

/* <u, v> = Sum u_i conj(v_i) over n entries (interleaved when cplx). */
static void om_dot(const double* u, const double* v, int n, bool cplx,
                   double* re, double* im) {
    double sr = 0.0, si = 0.0;
    if (!cplx) { for (int i = 0; i < n; i++) sr += u[i] * v[i]; }
    else {
        for (int i = 0; i < n; i++) {
            double ar = u[2 * i], ai = u[2 * i + 1], br = v[2 * i], bi = -v[2 * i + 1];
            sr += ar * br - ai * bi;
            si += ar * bi + ai * br;
        }
    }
    *re = sr; *im = si;
}

static double om_norm(const double* u, int n, bool cplx) {
    double re, im;
    om_dot(u, u, n, cplx, &re, &im);
    return sqrt(re > 0.0 ? re : 0.0);
}

/* w -= c * e */
static void om_axpy(double* w, const double* e, int n, bool cplx, double cr, double ci) {
    if (!cplx) { for (int i = 0; i < n; i++) w[i] -= cr * e[i]; return; }
    for (int i = 0; i < n; i++) {
        double er = e[2 * i], ei = e[2 * i + 1];
        w[2 * i]     -= cr * er - ci * ei;
        w[2 * i + 1] -= cr * ei + ci * er;
    }
}

/* Gram-Schmidt family over the k x n row-major buffer V (rows = vectors),
 * writing orthonormal rows (or zero rows) into E. */
static void om_gram_schmidt(const double* V, double* E, int k, int n, bool cplx,
                            OrthMethod meth, double tol) {
    size_t comps = cplx ? 2 : 1, rl = comps * (size_t)n;
    double* w = (double*)malloc(sizeof(double) * rl);
    bool* live = (bool*)calloc((size_t)(k > 0 ? k : 1), sizeof(bool));
    int passes = (meth == ORTH_AUTO || meth == ORTH_REORTH) ? 2 : 1;
    for (int r = 0; r < k; r++) {
        const double* v = V + (size_t)r * rl;
        memcpy(w, v, sizeof(double) * rl);
        double nv = om_norm(v, n, cplx);
        for (int pass = 0; pass < passes; pass++) {
            if (meth == ORTH_GS) {
                /* classical: all coefficients from the original vector */
                double* c = (double*)malloc(sizeof(double) * 2 * (size_t)(r > 0 ? r : 1));
                for (int j = 0; j < r; j++)
                    if (live[j]) om_dot(v, E + (size_t)j * rl, n, cplx, &c[2 * j], &c[2 * j + 1]);
                for (int j = 0; j < r; j++)
                    if (live[j]) om_axpy(w, E + (size_t)j * rl, n, cplx, c[2 * j], c[2 * j + 1]);
                free(c);
            } else {
                for (int j = 0; j < r; j++) {
                    if (!live[j]) continue;
                    double cr, ci;
                    om_dot(w, E + (size_t)j * rl, n, cplx, &cr, &ci);
                    om_axpy(w, E + (size_t)j * rl, n, cplx, cr, ci);
                }
            }
        }
        double nw = om_norm(w, n, cplx);
        double* e = E + (size_t)r * rl;
        if (nv == 0.0 || nw <= tol * nv) {
            memset(e, 0, sizeof(double) * rl);
        } else {
            for (size_t i = 0; i < rl; i++) e[i] = w[i] / nw;
            live[r] = true;
        }
    }
    free(w); free(live);
}

/* Householder: Q of the QR factorisation of the n x k matrix whose columns are
 * the vectors.  False when the vectors are dependent (or k > n). */
static bool om_householder(const double* V, double* E, int k, int n, bool cplx,
                           double tol) {
    if (k > n) return false;
    size_t comps = cplx ? 2 : 1;
    /* The row-major k x n buffer IS the column-major n x k matrix. */
    double* A = (double*)malloc(sizeof(double) * comps * (size_t)n * (size_t)k);
    memcpy(A, V, sizeof(double) * comps * (size_t)n * (size_t)k);
    double* tau = (double*)malloc(sizeof(double) * comps * (size_t)k);
    int info = cplx ? mat_lapack_zgeqrf(n, k, A, n, tau) : mat_lapack_dgeqrf(n, k, A, n, tau);
    bool ok = info == 0;
    for (int i = 0; ok && i < k; i++) {
        size_t d = (size_t)i + (size_t)i * (size_t)n;
        double rii = cplx ? hypot(A[2 * d], A[2 * d + 1]) : fabs(A[d]);
        double nv = om_norm(V + comps * (size_t)i * (size_t)n, n, cplx);
        if (!(rii > tol * nv)) ok = false;
    }
    if (ok) {
        info = cplx ? mat_lapack_zungqr(n, k, k, A, n, tau) : mat_lapack_dorgqr(n, k, k, A, n, tau);
        ok = info == 0;
    }
    if (ok) memcpy(E, A, sizeof(double) * comps * (size_t)n * (size_t)k);
    free(A); free(tau);
    return ok;
}

/* The machine kernel: `m` is an NDArray or a boxed k x n machine matrix.
 * tol < 0 selects the default. */
static OrthStatus orth_machine(const Expr* m, OrthMethod meth, double tol, Expr** out) {
    *out = NULL;
    int k = 0, n = 0;
    double* V = NULL;
    bool cplx = false;
    if (!na_load_matrix((Expr*)m, false, false, &k, &n, &V)) {
        free(V); V = NULL;
        if (!na_load_matrix((Expr*)m, true, false, &k, &n, &V)) { free(V); return OM_DECLINE; }
        cplx = true;
    }
    if (k < 1 || n < 1) { free(V); return OM_DECLINE; }
    if (tol < 0.0) tol = 64.0 * DBL_EPSILON * (double)(k > n ? k : n);
    size_t comps = cplx ? 2 : 1;
    double* E = (double*)calloc(comps * (size_t)k * (size_t)n, sizeof(double));
    if (meth == ORTH_HOUSEHOLDER) {
        if (!om_householder(V, E, k, n, cplx, tol)) {
            free(V); free(E);
            mth_message(ORTH_HEAD, "ornhsl", "The Householder method called by the option "
                        "Method -> Householder cannot orthogonalize this set of vectors. "
                        "Try another orthogonalization method.");
            return OM_FAIL;
        }
    } else {
        om_gram_schmidt(V, E, k, n, cplx, meth, tol);
    }
    *out = na_build_matrix_as(E, k, n, cplx, false, na_result_presentation(m));
    free(V); free(E);
    return *out ? OM_OK : OM_DECLINE;
}

Expr* la_orthogonalize_machine(const Expr* m) {
    if (!is_ndarray(m) || m->data.ndarray.rank != 2) return NULL;
    Expr* out = NULL;
    if (orth_machine(m, ORTH_AUTO, -1.0, &out) != OM_OK) { if (out) expr_free(out); return NULL; }
    return out;
}

/* ------------------------------------------------------------------------ *
 *  Generic kernel: exact, symbolic, MPFR, or a user inner product          *
 * ------------------------------------------------------------------------ */

/* ip(u, v): u . Conjugate[v] by default, else f[u, v]; evaluated. */
static Expr* og_ip(const Expr* f, const Expr* u, const Expr* v) {
    if (f) return eval_and_free(expr_new_function(expr_copy((Expr*)f),
                                (Expr*[]){ expr_copy((Expr*)u), expr_copy((Expr*)v) }, 2));
    return la_ev2(SYM_Dot, expr_copy((Expr*)u), la_call1(SYM_Conjugate, expr_copy((Expr*)v)));
}

/* Norm squared under the inner product (default: Norm[w]^2, which keeps the
 * reference's Abs[a]^2 + Abs[b]^2 form). */
static Expr* og_nsq(const Expr* f, const Expr* w) {
    if (f) return og_ip(f, w, w);
    return la_ev2(SYM_Power, la_ev1(SYM_Norm, expr_copy((Expr*)w)), expr_new_integer(2));
}

/* Numeric value of a scalar, if it has one. */
static bool og_num(const Expr* e, double* out) {
    if (!la_all_numeric(e)) return false;
    Expr* nv = la_ev1(SYM_N, expr_copy((Expr*)e));
    double re = 0.0, im = 0.0;
    bool ok = na_read_scalar(nv, &re, &im);
    expr_free(nv);
    if (ok) *out = re;
    return ok;
}

/* True when the residual w is zero: an exact zero norm, a numeric norm under
 * `tol` relative to the vector's (tol: the Tolerance, or a precision floor
 * for an inexact input), or PossibleZeroQ for an exact numeric one. */
static bool og_is_zero(const Expr* nsq_w, const Expr* nsq_v, double tol, bool inexact) {
    if (nsq_w->type == EXPR_INTEGER && nsq_w->data.integer == 0) return true;
    double a, b;
    if (og_num(nsq_w, &a) && og_num(nsq_v, &b)) {
        if (b <= 0.0) return true;
        if (tol >= 0.0) return sqrt(fabs(a)) < tol * sqrt(b);
        if (inexact) return sqrt(fabs(a)) <= 1e-10 * sqrt(b) && a < 1e-24 * b;
        if (a > 1e-12 * b) return false;
        Expr* z = la_ev1(SYM_PossibleZeroQ, expr_copy((Expr*)nsq_w));
        bool zero = z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_True;
        expr_free(z);
        return zero;
    }
    return false;
}

/* True when `e` is an Integer, a Rational, or a Complex of those. */
static bool og_gauss_rational(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || head_is(e, SYM_Rational)) return true;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2)
        return og_gauss_rational(e->data.function.args[0])
            && og_gauss_rational(e->data.function.args[1]);
    return false;
}

/* An exact residual with (Gaussian) rational entries, rescaled to coprime
 * Gaussian-integer entries before it is normalised -- the direction is all
 * that matters, and dividing a primitive integer vector by the Sqrt of its
 * norm gives the reference's form ((2 + 3 I)/Sqrt[26], not
 * (1 + 3 I/2)/Sqrt[13/2]).  Anything else is returned unchanged (adopted). */
static Expr* og_primitive(Expr* w) {
    if (!head_is(w, SYM_List) || w->data.function.arg_count == 0) return w;
    for (size_t i = 0; i < w->data.function.arg_count; i++)
        if (!og_gauss_rational(w->data.function.args[i])) return w;
    static Expr* prim = NULL;
    if (!prim)
        prim = parse_expression(
            "With[{d = LCM @@ Denominator[Flatten[{Re[#], Im[#]}]]}, "
            "With[{v = d #}, v / GCD @@ Flatten[{Re[v], Im[v]}]]] &");
    if (!prim) return w;
    Expr* r = eval_and_free(expr_new_function(expr_copy(prim), (Expr*[]){ w }, 1));
    return r;
}

/* A zero "vector" shaped like v: 0 * v, evaluated (a list of zeros, or 0). */
static Expr* og_zero_like(const Expr* v) {
    return la_ev2(SYM_Times, expr_new_integer(0), expr_copy((Expr*)v));
}

static Expr* orth_generic(const Expr* vs, const Expr* f, double tol, bool inexact) {
    size_t k = vs->data.function.arg_count;
    Expr** w = (Expr**)calloc(k ? k : 1, sizeof(Expr*));     /* residuals    */
    Expr** wn = (Expr**)calloc(k ? k : 1, sizeof(Expr*));    /* <w_j, w_j>   */
    Expr** e = (Expr**)calloc(k ? k : 1, sizeof(Expr*));     /* outputs      */
    for (size_t r = 0; r < k; r++) {
        const Expr* v = vs->data.function.args[r];
        Expr* cur = expr_copy((Expr*)v);
        for (size_t j = 0; j < r; j++) {
            if (!w[j]) continue;
            Expr* c = la_ev2(SYM_Times, og_ip(f, v, w[j]),
                             la_call2(SYM_Power, expr_copy(wn[j]), expr_new_integer(-1)));
            cur = la_ev2(SYM_Plus, cur,
                         la_ev2(SYM_Times, la_ev2(SYM_Times, expr_new_integer(-1), c),
                                expr_copy(w[j])));
        }
        Expr* nsq = og_nsq(f, cur);
        Expr* nsq_v = og_nsq(f, v);
        if (og_is_zero(nsq, nsq_v, tol, inexact)) {
            e[r] = og_zero_like(v);
            expr_free(cur); expr_free(nsq);
        } else {
            if (!f && !inexact) {
                Expr* pw = og_primitive(expr_copy(cur));
                Expr* pn = og_nsq(NULL, pw);
                e[r] = la_ev2(SYM_Times, pw,
                              la_call2(SYM_Power, la_ev1(SYM_Sqrt, pn), expr_new_integer(-1)));
            } else {
                e[r] = la_ev2(SYM_Times, expr_copy(cur),
                              la_call2(SYM_Power, la_ev1(SYM_Sqrt, expr_copy(nsq)),
                                       expr_new_integer(-1)));
            }
            w[r] = cur; wn[r] = nsq;
        }
        expr_free(nsq_v);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), e, k);
    for (size_t r = 0; r < k; r++) { if (w[r]) expr_free(w[r]); if (wn[r]) expr_free(wn[r]); }
    free(w); free(wn); free(e);
    return out;
}

/* ------------------------------------------------------------------------ *
 *  Entry points                                                            *
 * ------------------------------------------------------------------------ */

static bool orth_is_rule(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed));
}

static Expr* orth_run(Expr* res) {
    size_t argc = res->data.function.arg_count;
    size_t npos = argc;
    while (npos > 1 && orth_is_rule(res->data.function.args[npos - 1])) npos--;
    if (npos < 1 || npos > 2) return builtin_arg_error(ORTH_HEAD, npos, 1, 2);

    OrthMethod meth = ORTH_AUTO;
    double tol = -1.0;
    for (size_t i = npos; i < argc; i++) {
        const Expr* key = res->data.function.args[i]->data.function.args[0];
        const Expr* rhs = res->data.function.args[i]->data.function.args[1];
        const char* kn = key->type == EXPR_SYMBOL ? key->data.symbol.name : NULL;
        if (kn == SYM_Method) {
            if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) meth = ORTH_AUTO;
            else if (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "GramSchmidt") == 0) meth = ORTH_GS;
            else if (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "ModifiedGramSchmidt") == 0) meth = ORTH_MGS;
            else if (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "Householder") == 0) meth = ORTH_HOUSEHOLDER;
            else if (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "Reorthogonalization") == 0) meth = ORTH_REORTH;
            else {
                char* s = rhs->type == EXPR_STRING ? NULL : la_str(rhs);
                mth_message(ORTH_HEAD, "orbop", "%s is not a valid value for the Method option.",
                            s ? s : rhs->data.string);
                free(s);
                return NULL;
            }
        } else if (kn == SYM_Tolerance) {
            if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) tol = -1.0;
            else {
                double t;
                if (!og_num(rhs, &t) || t < 0.0) {
                    char* s = la_str(rhs);
                    mth_message(ORTH_HEAD, "tol", "Value of option Tolerance -> %s is not a "
                                "non-negative real number or Automatic.", s);
                    free(s);
                    return NULL;
                }
                tol = t;
            }
        } else {
            char* k = la_str(key);
            char* call = la_str(res);
            mth_message(ORTH_HEAD, "optx", "Unknown option %s in %s.", k, call);
            free(k); free(call);
            return NULL;
        }
    }

    const Expr* vs = res->data.function.args[0];
    const Expr* f = npos == 2 ? res->data.function.args[1] : NULL;
    /* Dot is the default inner product's real form; keep the machine path. */
    if (f && f->type == EXPR_SYMBOL && f->data.symbol.name == SYM_Dot && meth != ORTH_HOUSEHOLDER
        && !la_has_complex(vs))
        f = NULL;
    if (f && meth == ORTH_HOUSEHOLDER) {
        mth_message(ORTH_HEAD, "ornhip", "The option Method -> Householder can only be used "
                    "with the default (Automatic) inner product.");
        return NULL;
    }
    if (!is_ndarray(vs) && !head_is(vs, SYM_List)) return NULL;
    if (head_is(vs, SYM_List) && vs->data.function.arg_count == 0) return la_empty_list();

    int k = 0, n = 0;
    bool is_matrix = la_matrix_dims(vs, &k, &n);
    if (!f) {
        if (!is_matrix) {
            char* s = la_str(vs);
            if (la_looks_like_matrix(vs))
                mth_message(ORTH_HEAD, "ornt", "The first argument %s is incorrect. It should "
                            "be a vector, matrix, or tensor.", s);
            else
                mth_message(ORTH_HEAD, "orwdt", "At least one element of the first argument "
                            "cannot be used with the Dot inner product.");
            free(s);
            return NULL;
        }
        long bits = 53;
        LaPrecClass pc = la_prec_class(vs, &bits);
        if (pc == LA_MACHINE) {
            Expr* out = NULL;
            OrthStatus st = orth_machine(vs, meth, tol, &out);
            if (st == OM_OK) return out;
            if (st == OM_FAIL) return NULL;
        }
        /* arbitrary precision: a residual below ~2^(8 - bits) of its vector
         * is rounding, not direction */
        if (pc == LA_MPFR && tol < 0.0) tol = ldexp(1.0, -(int)(bits - 8));
        Expr* boxed = la_boxed_copy(vs);
        Expr* out = orth_generic(boxed, NULL, tol, pc != LA_EXACT);
        expr_free(boxed);
        return out;
    }

    /* A user inner product: it must give a number or a symbol, and a
     * non-negative real for a numeric element paired with itself. */
    Expr* boxed = la_boxed_copy(vs);
    Expr* probe = og_ip(f, boxed->data.function.args[0], boxed->data.function.args[0]);
    double pv = 0.0;
    bool bad = head_is(probe, SYM_List)
            || (la_all_numeric(boxed->data.function.args[0]) && la_all_numeric(probe)
                && (!og_num(probe, &pv) || pv < 0.0));
    if (!bad && la_all_numeric(probe) && la_has_complex(probe)) bad = true;
    expr_free(probe);
    if (bad) {
        char* s = la_str(f);
        mth_message(ORTH_HEAD, "ornfa", "The second argument %s is not an inner product "
                    "function, which always should return a number or symbol. If its two "
                    "numeric arguments are the same, it should return a non-negative real "
                    "number.", s);
        free(s);
        expr_free(boxed);
        return NULL;
    }
    Expr* out = orth_generic(boxed, f, tol, la_prec_class(vs, NULL) != LA_EXACT);
    expr_free(boxed);
    return out;
}

Expr* ndla_orthogonalize(Expr* res) { return orth_run(res); }

Expr* builtin_orthogonalize(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1) return builtin_arg_error(ORTH_HEAD, 0, 1, 2);
    if (linalg_call_has_ndarray(res)) return ndla_orthogonalize(res);
    return orth_run(res);
}

void orthogonalize_init(void) {
    symtab_add_builtin(ORTH_HEAD, builtin_orthogonalize);
    symtab_get_def(ORTH_HEAD)->attributes |= ATTR_PROTECTED;
    Expr* r[2] = {
        la_call2(SYM_Rule, expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic)),
        la_call2(SYM_Rule, expr_new_symbol(SYM_Tolerance), expr_new_symbol(SYM_Automatic)),
    };
    symtab_set_options(ORTH_HEAD, expr_new_function(expr_new_symbol(SYM_List), r, 2));
}
