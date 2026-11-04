/*
 * src/linalg/singularvaluelist.c -- SingularValueList.
 *
 *   SingularValueList[m]          the nonzero singular values, decreasing
 *   SingularValueList[m, k]       the first k (k < 0: the last |k|) of ALL
 *                                 min(rows, cols) values, zeros included;
 *                                 also {k}, {i, j}, {i, j, s}, UpTo[k]
 *   SingularValueList[{m, a}]     generalized singular values (machine only)
 *   Options: Tolerance, Method (Automatic | "Arnoldi"), Cubics, Quartics.
 *
 * Reference semantics (wolframscript): a value s is dropped when
 * s < tol * s_max (tol: the Tolerance value, or max(rows, cols) * eps for a
 * machine matrix and 0 for an exact one); an all-zero matrix keeps its zeros
 * ({0, 0}); with an explicit k every value is a candidate and asking for more
 * than exist gives ::take.
 *
 * Dispatch by precision class:
 *   machine -> svl_machine_values: LAPACK gesdd with jobz = 'N' (values only,
 *              no U or V).  The one helper behind the boxed List, packed List,
 *              visible NDArray and Compile[], so they agree bit for bit.
 *   MPFR    -> svd_mpfr_values: one-sided Jacobi with V not accumulated.
 *   exact / symbolic -> Sqrt of the Eigenvalues of the smaller Gram matrix
 *              (ConjugateTranspose[m] . m or m . ConjugateTranspose[m]).
 */

#include "la_dense.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "lapack.h"
#include "svdecomp_internal.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#define SVL_HEAD "SingularValueList"

/* Parsed call. */
typedef struct {
    const Expr* m;          /* matrix (or first of the pair)             */
    const Expr* a;          /* second of the generalized pair, or NULL   */
    bool   has_k;           /* a positional sequence spec was given      */
    bool   upto;            /* ... and it was UpTo[k]                    */
    long long k_lo, k_hi, k_step;  /* 1-based (negative: from the end)   */
    long long upto_k;
    bool   tol_auto;
    double tol;             /* relative tolerance when !tol_auto         */
    bool   cubics, quartics;
} SvlArgs;

/* ------------------------------------------------------------------------ *
 *  Argument parsing                                                        *
 * ------------------------------------------------------------------------ */

static bool svl_is_rule(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed));
}

static void svl_seqs_msg(const Expr* res) {
    char* s = la_str(res);
    mth_message(SVL_HEAD, "seqs", "Sequence specification (+n, -n, {+n}, {-n}, {m, n}, "
                "or {m, n, s}) expected at position 2 in %s.", s);
    free(s);
}

static bool svl_int(const Expr* e, long long* v) {
    if (e->type != EXPR_INTEGER) return false;
    *v = (long long)e->data.integer;
    return true;
}

/* k, {k}, {i, j}, {i, j, s}, UpTo[k] */
static bool svl_parse_seq(const Expr* spec, SvlArgs* a) {
    long long k;
    if (svl_int(spec, &k)) {
        if (k >= 0) { a->k_lo = 1; a->k_hi = k; } else { a->k_lo = k; a->k_hi = -1; }
        a->k_step = 1;
        return true;
    }
    if (head_is(spec, SYM_UpTo)) {
        bool up;
        if (!la_parse_take_spec(spec, &k, &up)) return false;
        a->upto = true; a->upto_k = k;
        return true;
    }
    if (head_is(spec, SYM_List)) {
        size_t c = spec->data.function.arg_count;
        long long v[3];
        if (c < 1 || c > 3) return false;
        for (size_t i = 0; i < c; i++)
            if (!svl_int(spec->data.function.args[i], &v[i])) return false;
        if (c == 1) { a->k_lo = v[0]; a->k_hi = v[0]; a->k_step = 1; }
        else { a->k_lo = v[0]; a->k_hi = v[1]; a->k_step = c == 3 ? v[2] : 1; }
        if (a->k_step == 0) return false;
        return true;
    }
    return false;
}

static bool svl_bool_opt(const Expr* rhs, bool* out) {
    if (rhs->type != EXPR_SYMBOL) return false;
    if (rhs->data.symbol.name == SYM_True)  { *out = true;  return true; }
    if (rhs->data.symbol.name == SYM_False) { *out = false; return true; }
    return false;
}

static bool svl_parse(const Expr* res, SvlArgs* a) {
    memset(a, 0, sizeof *a);
    a->tol_auto = true;
    size_t argc = res->data.function.arg_count;
    size_t npos = argc;
    while (npos > 1 && svl_is_rule(res->data.function.args[npos - 1])) npos--;
    if (npos < 1 || npos > 2) { builtin_arg_error(SVL_HEAD, npos, 1, 2); return false; }

    for (size_t i = npos; i < argc; i++) {
        const Expr* o = res->data.function.args[i];
        const Expr* key = o->data.function.args[0];
        const Expr* rhs = o->data.function.args[1];
        const char* kn = key->type == EXPR_SYMBOL ? key->data.symbol.name : NULL;
        if (kn == SYM_Tolerance) {
            if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
                a->tol_auto = true;
            } else {
                double re = 0.0, im = 0.0;
                Expr* nv = la_ev1(SYM_N, expr_copy((Expr*)rhs));
                bool ok = na_read_scalar(nv, &re, &im) && im == 0.0 && re >= 0.0;
                expr_free(nv);
                if (!ok) {
                    char* s = la_str(rhs);
                    mth_message(SVL_HEAD, "tol", "Value of option Tolerance -> %s is not a "
                                "non-negative real number or Automatic.", s);
                    free(s);
                    return false;
                }
                a->tol_auto = false; a->tol = re;
            }
        } else if (kn == SYM_Method) {
            bool ok = (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
                   || (rhs->type == EXPR_STRING && strcmp(rhs->data.string, "Arnoldi") == 0)
                   || (head_is(rhs, SYM_List) && rhs->data.function.arg_count >= 1
                       && rhs->data.function.args[0]->type == EXPR_STRING
                       && strcmp(rhs->data.function.args[0]->data.string, "Arnoldi") == 0);
            if (!ok) {
                char* s = rhs->type == EXPR_STRING ? NULL : la_str(rhs);
                mth_message(SVL_HEAD, "emeth", "The method specified by Method -> %s should "
                            "be either Automatic or Arnoldi.",
                            s ? s : rhs->data.string);
                free(s);
                return false;
            }
        } else if (kn == SYM_Cubics) {
            if (!svl_bool_opt(rhs, &a->cubics)) return false;
        } else if (kn == SYM_Quartics) {
            if (!svl_bool_opt(rhs, &a->quartics)) return false;
        } else {
            char* k = la_str(key);
            char* call = la_str(res);
            mth_message(SVL_HEAD, "optx", "Unknown option %s in %s.", k, call);
            free(k); free(call);
            return false;
        }
    }

    const Expr* arg0 = res->data.function.args[0];
    /* {m, a}: a length-2 List whose elements are both matrices */
    int r1, c1, r2, c2;
    if (head_is(arg0, SYM_List) && arg0->data.function.arg_count == 2
        && la_matrix_dims(arg0->data.function.args[0], &r1, &c1)
        && la_matrix_dims(arg0->data.function.args[1], &r2, &c2)) {
        a->m = arg0->data.function.args[0];
        a->a = arg0->data.function.args[1];
    } else {
        a->m = arg0;
    }

    if (npos == 2) {
        a->has_k = true;
        if (!svl_parse_seq(res->data.function.args[1], a)) { svl_seqs_msg(res); return false; }
        if (a->upto) a->has_k = false;
    }
    return true;
}

/* ------------------------------------------------------------------------ *
 *  Selection: tolerance / sequence spec over the descending full list      *
 * ------------------------------------------------------------------------ */

/* Indices (0-based, into the full descending list of `count`) selected by
 * the call.  `nz` is how many leading values survive the tolerance.  Returns
 * the selection length, or -1 after a ::take message. */
static long long svl_select(const SvlArgs* a, int count, int nz, int* idx) {
    if (!a->has_k) {
        /* UpTo[k] draws from ALL the values, zeros included, as an explicit k
         * does (reference: SingularValueList[{{1, 2}, {2, 4}}, UpTo[2]] is
         * {5, 0}); with no spec only the nonzero values. */
        long long take = a->upto ? (a->upto_k < count ? a->upto_k : count) : nz;
        for (long long i = 0; i < take; i++) idx[i] = (int)i;
        return take;
    }
    long long lo = a->k_lo, hi = a->k_hi;
    if (a->k_lo >= 1 && a->k_hi == 0) return 0;               /* k = 0 */
    long long plo = lo < 0 ? count + lo + 1 : lo;
    long long phi = hi < 0 ? count + hi + 1 : hi;
    if (plo < 1 || phi < 1 || plo > count || phi > count) {
        mth_message(SVL_HEAD, "take", "Cannot take singular values %lld through %lld out "
                    "of the total of %d singular values.", lo, hi, count);
        return -1;
    }
    long long n = 0;
    if (a->k_step > 0) for (long long i = plo; i <= phi; i += a->k_step) idx[n++] = (int)(i - 1);
    else               for (long long i = plo; i >= phi; i += a->k_step) idx[n++] = (int)(i - 1);
    return n;
}

/* Number of leading values of the descending S[0..count) that survive the
 * tolerance; an all-zero list keeps every value, as the reference does. */
static int svl_nonzero(const double* S, int count, double tol_rel) {
    if (count == 0 || S[0] == 0.0) return count;
    double thr = tol_rel * S[0];
    int nz = 0;
    while (nz < count && !(S[nz] < thr) && S[nz] != 0.0) nz++;
    return nz;
}

/* ------------------------------------------------------------------------ *
 *  Machine kernel                                                          *
 * ------------------------------------------------------------------------ */

/* Values only: gesdd jobz = 'N'.  Returns a malloc'd descending array of
 * min(rows, cols) values, or NULL when `m` is not machine-numeric or LAPACK
 * is unavailable. */
static double* svl_machine_values(const Expr* m, int* rows_out, int* cols_out) {
    int r = 0, c = 0;
    double* A = NULL;
    bool cplx = false;
    if (!na_load_matrix((Expr*)m, false, true, &r, &c, &A)) {
        free(A); A = NULL;
        if (!na_load_matrix((Expr*)m, true, true, &r, &c, &A)) { free(A); return NULL; }
        cplx = true;
    }
    if (r < 1 || c < 1) { free(A); return NULL; }
    int mn = r < c ? r : c;
    double* S = (double*)malloc(sizeof(double) * (size_t)mn);
    double dummy[2] = { 0.0, 0.0 };
    int info = cplx ? mat_lapack_zgesdd('N', r, c, A, r, S, dummy, 1, dummy, 1)
                    : mat_lapack_dgesdd('N', r, c, A, r, S, dummy, 1, dummy, 1);
    free(A);
    if (info != 0) { free(S); return NULL; }
    *rows_out = r; *cols_out = c;
    return S;
}

/* The machine answer for parsed args over a machine matrix: a packed real
 * vector, NULL to fall through (not machine-numeric), or *failed after a
 * message. */
static Expr* svl_machine(const SvlArgs* a, const Expr* m, bool* failed) {
    int r = 0, c = 0;
    double* S = svl_machine_values(m, &r, &c);
    if (!S) return NULL;
    int mn = r < c ? r : c;
    double tol = a->tol_auto ? (double)(r > c ? r : c) * DBL_EPSILON : a->tol;
    int nz = svl_nonzero(S, mn, tol);
    for (int i = nz; i < mn; i++) S[i] = 0.0;   /* under the tolerance: 0., as the reference */
    int* idx = (int*)malloc(sizeof(int) * (size_t)(mn > 0 ? mn : 1) * 2 + sizeof(int));
    long long n = svl_select(a, mn, nz, idx);
    Expr* out = NULL;
    if (n < 0) { *failed = true; }
    else {
        double* v = (double*)malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
        for (long long i = 0; i < n; i++) v[i] = S[idx[i]];
        out = n > 0 ? na_build_vector_as(v, (int)n, false, na_result_presentation(m))
                    : la_empty_list();
        free(v);
    }
    free(idx); free(S);
    return out;
}

Expr* la_singularvaluelist_machine(const Expr* m) {
    if (!is_ndarray(m) || m->data.ndarray.rank != 2) return NULL;
    SvlArgs a; memset(&a, 0, sizeof a); a.tol_auto = true;
    bool failed = false;
    Expr* out = svl_machine(&a, m, &failed);
    if (out && !is_ndarray(out)) { expr_free(out); return NULL; }  /* {} : not an array */
    return out;
}

/* Generalized singular values of {m, a} (machine only): finite alpha/beta
 * from ggsvd3, descending. */
static Expr* svl_generalized(const SvlArgs* args, bool* failed) {
    int r1, c1, r2, c2;
    double *A = NULL, *B = NULL;
    bool cplx = false;
    bool machine = la_prec_class(args->m, NULL) == LA_MACHINE
                && la_prec_class(args->a, NULL) == LA_MACHINE;
    if (machine) {
        if (!(na_load_matrix((Expr*)args->m, false, true, &r1, &c1, &A)
              && na_load_matrix((Expr*)args->a, false, true, &r2, &c2, &B))) {
            free(A); free(B); A = B = NULL;
            cplx = true;
            if (!(na_load_matrix((Expr*)args->m, true, true, &r1, &c1, &A)
                  && na_load_matrix((Expr*)args->a, true, true, &r2, &c2, &B))) {
                free(A); free(B); A = B = NULL;
            }
        }
    }
    if (!A || !B || c1 != c2) {
        free(A); free(B);
        mth_message(SVL_HEAD, "gfargs", "Generalized SingularValueList arguments accept only "
                    "matrices with machine real and complex elements.");
        *failed = true;
        return NULL;
    }
    int n = c1, k = 0, l = 0;
    double* alpha = (double*)calloc((size_t)n, sizeof(double));
    double* beta  = (double*)calloc((size_t)n, sizeof(double));
    double dummy[2] = { 0.0, 0.0 };
    int info = cplx
        ? mat_lapack_zggsvd3('N', 'N', 'N', r1, n, r2, &k, &l, A, r1, B, r2,
                             alpha, beta, dummy, 1, dummy, 1, dummy, 1)
        : mat_lapack_dggsvd3('N', 'N', 'N', r1, n, r2, &k, &l, A, r1, B, r2,
                             alpha, beta, dummy, 1, dummy, 1, dummy, 1);
    free(A); free(B);
    if (info != 0) { free(alpha); free(beta); return NULL; }
    double* S = (double*)malloc(sizeof(double) * (size_t)(l > 0 ? l : 1));
    int cnt = 0;
    for (int i = k; i < k + l && i < n; i++)
        if (beta[i] > 0.0) S[cnt++] = alpha[i] / beta[i];
    free(alpha); free(beta);
    for (int i = 1; i < cnt; i++)
        for (int j = i; j > 0 && S[j - 1] < S[j]; j--) { double t = S[j]; S[j] = S[j - 1]; S[j - 1] = t; }
    double tol = args->tol_auto ? (double)(r1 > n ? r1 : n) * DBL_EPSILON : args->tol;
    int nz = svl_nonzero(S, cnt, tol);
    for (int i = nz; i < cnt; i++) S[i] = 0.0;
    int* idx = (int*)malloc(sizeof(int) * (size_t)(cnt + 1) * 2);
    long long sel = svl_select(args, cnt, nz, idx);
    Expr* out = NULL;
    if (sel < 0) *failed = true;
    else {
        double* v = (double*)malloc(sizeof(double) * (size_t)(sel + 1));
        for (long long i = 0; i < sel; i++) v[i] = S[idx[i]];
        out = sel > 0 ? na_build_vector_as(v, (int)sel, false, na_result_presentation(args->m))
                      : la_empty_list();
        free(v);
    }
    free(idx); free(S);
    return out;
}

/* ------------------------------------------------------------------------ *
 *  MPFR kernel                                                             *
 * ------------------------------------------------------------------------ */

#ifdef USE_MPFR
static Expr* svl_mpfr(const SvlArgs* a, const Expr* m, int r, int c, long bits,
                      bool* failed) {
    Expr* boxed = la_boxed_copy(m);
    if (c > r) {                          /* Jacobi wants a tall matrix */
        boxed = la_ev1(SYM_ConjugateTranspose, boxed);
        int t = r; r = c; c = t;
    }
    mpfr_t* S = svd_mpfr_values(boxed, r, c, (mpfr_prec_t)bits);
    expr_free(boxed);
    if (!S) return NULL;
    int mn = c;
    double* Sd = (double*)malloc(sizeof(double) * (size_t)mn);
    for (int i = 0; i < mn; i++) Sd[i] = mpfr_get_d(S[i], MPFR_RNDN);
    /* default: 8 max(r, c) 2^-bits relative -- the Jacobi floor at `bits` */
    double tol = a->tol_auto ? 8.0 * (double)r * ldexp(1.0, -(int)bits) : a->tol;
    int nz = svl_nonzero(Sd, mn, tol);
    for (int i = nz; i < mn; i++) mpfr_set_zero(S[i], 1);
    int* idx = (int*)malloc(sizeof(int) * (size_t)(mn + 1) * 2);
    long long n = svl_select(a, mn, nz, idx);
    Expr* out = NULL;
    if (n < 0) *failed = true;
    else {
        Expr** v = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
        for (long long i = 0; i < n; i++) v[i] = expr_new_mpfr_copy(S[idx[i]]);
        out = expr_new_function(expr_new_symbol(SYM_List), v, (size_t)n);
        free(v);
    }
    free(idx); free(Sd);
    svd_mpfr_values_free(S, mn);
    return out;
}
#endif

/* ------------------------------------------------------------------------ *
 *  Exact / symbolic kernel                                                 *
 * ------------------------------------------------------------------------ */

static Expr* svl_exact(const SvlArgs* a, const Expr* m, int r, int c, bool* failed) {
    Expr* boxed = la_boxed_copy(m);
    Expr* ct = la_ev1(SYM_ConjugateTranspose, expr_copy(boxed));
    Expr* g = (r >= c) ? la_ev2(SYM_Dot, ct, boxed) : la_ev2(SYM_Dot, boxed, ct);
    Expr* ev = eval_and_free(expr_new_function(expr_new_symbol(SYM_Eigenvalues),
        (Expr*[]){ g,
                   la_call2(SYM_Rule, expr_new_symbol(SYM_Cubics),
                            expr_new_symbol(a->cubics ? SYM_True : SYM_False)),
                   la_call2(SYM_Rule, expr_new_symbol(SYM_Quartics),
                            expr_new_symbol(a->quartics ? SYM_True : SYM_False)) }, 3));
    int mn = r < c ? r : c;
    if (!head_is(ev, SYM_List) || (int)ev->data.function.arg_count != mn) {
        expr_free(ev);
        return NULL;
    }
    Expr** vals = (Expr**)malloc(sizeof(Expr*) * (size_t)mn);
    double* num = (double*)calloc((size_t)(mn > 0 ? mn : 1), sizeof(double));
    bool all_numeric = true;
    for (int i = 0; i < mn; i++) {
        /* Together cancels the content a radical eigenvalue often carries,
         * (30 + 2 Sqrt[221])/2 -> 15 + Sqrt[221], the reference's form. */
        Expr* lam = expr_copy(ev->data.function.args[i]);
        if (la_all_numeric(lam)) lam = la_ev1(SYM_Together, lam);
        vals[i] = la_ev1(SYM_Sqrt, lam);
        num[i] = 0.0;
        double re = 0.0, im = 0.0;
        if (la_all_numeric(vals[i])) {
            Expr* nv = la_ev1(SYM_N, expr_copy(vals[i]));
            if (na_read_scalar(nv, &re, &im)) num[i] = re; else all_numeric = false;
            expr_free(nv);
        } else {
            all_numeric = false;
        }
    }
    expr_free(ev);
    int nz = mn;
    if (all_numeric) {
        /* descending by value (stable insertion sort) */
        for (int i = 1; i < mn; i++)
            for (int j = i; j > 0 && num[j - 1] < num[j]; j--) {
                double t = num[j]; num[j] = num[j - 1]; num[j - 1] = t;
                Expr* e = vals[j]; vals[j] = vals[j - 1]; vals[j - 1] = e;
            }
        /* exact zeros are dropped (an all-zero list kept); a Tolerance is
         * applied to the numerical values */
        if (num[0] != 0.0 || a->tol > 0.0) {
            nz = 0;
            double thr = a->tol_auto ? 0.0 : a->tol * num[0];
            bool all_zero = true;
            for (int i = 0; i < mn; i++) if (!(vals[i]->type == EXPR_INTEGER && vals[i]->data.integer == 0)) all_zero = false;
            if (all_zero) nz = mn;
            else while (nz < mn && !(vals[nz]->type == EXPR_INTEGER && vals[nz]->data.integer == 0)
                        && !(num[nz] < thr)) nz++;
        }
    } else {
        /* symbolic: drop only literal zeros */
        int w = 0;
        for (int i = 0; i < mn; i++) {
            bool zero = vals[i]->type == EXPR_INTEGER && vals[i]->data.integer == 0;
            if (zero && !a->has_k) { expr_free(vals[i]); continue; }
            vals[w++] = vals[i];
        }
        if (w == 0) {      /* an all-zero symbolic matrix: keep its zeros */
            for (int i = 0; i < mn; i++) vals[i] = expr_new_integer(0);
            w = mn;
        }
        mn = w; nz = w;
    }
    int* idx = (int*)malloc(sizeof(int) * (size_t)(mn + 1) * 2);
    long long n = svl_select(a, mn, nz, idx);
    Expr* out = NULL;
    if (n < 0) *failed = true;
    else {
        Expr** v = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
        for (long long i = 0; i < n; i++) v[i] = expr_copy(vals[idx[i]]);
        out = expr_new_function(expr_new_symbol(SYM_List), v, (size_t)n);
        free(v);
    }
    for (int i = 0; i < mn; i++) expr_free(vals[i]);
    free(vals); free(num); free(idx);
    return out;
}

/* ------------------------------------------------------------------------ *
 *  Entry points                                                            *
 * ------------------------------------------------------------------------ */

static Expr* svl_run(Expr* res) {
    SvlArgs a;
    if (!svl_parse(res, &a)) return NULL;
    bool failed = false;
    if (a.a) return svl_generalized(&a, &failed);

    int r = 0, c = 0;
    if (!la_matrix_dims(a.m, &r, &c) || r < 1 || c < 1) {
        char* s = la_str(a.m);
        mth_message(SVL_HEAD, "matrix", "Argument %s at position 1 is not a non-empty "
                    "rectangular matrix.", s);
        free(s);
        return NULL;
    }
    long bits = 53;
    LaPrecClass pc = la_prec_class(a.m, &bits);
    Expr* out = NULL;
    if (pc == LA_MACHINE) {
        out = svl_machine(&a, a.m, &failed);
        if (out || failed) return out;
        /* a symbolic leaf in an inexact matrix, or no LAPACK: numericalise
         * through the exact path at machine precision */
        Expr* nm = la_ev1(SYM_N, la_boxed_copy(a.m));
        out = svl_exact(&a, nm, r, c, &failed);
        expr_free(nm);
        return out;
    }
#ifdef USE_MPFR
    if (pc == LA_MPFR) {
        out = svl_mpfr(&a, a.m, r, c, bits, &failed);
        if (out || failed) return out;
    }
#endif
    return svl_exact(&a, a.m, r, c, &failed);
}

Expr* ndla_singularvaluelist(Expr* res) {
    return svl_run(res);
}

Expr* builtin_singularvaluelist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1)
        return builtin_arg_error(SVL_HEAD, 0, 1, 2);
    if (linalg_call_has_ndarray(res)) return ndla_singularvaluelist(res);
    return svl_run(res);
}

void singularvaluelist_init(void) {
    symtab_add_builtin(SVL_HEAD, builtin_singularvaluelist);
    symtab_get_def(SVL_HEAD)->attributes |= ATTR_PROTECTED;
    Expr* r[4] = {
        la_call2(SYM_Rule, expr_new_symbol(SYM_Cubics), expr_new_symbol(SYM_False)),
        la_call2(SYM_Rule, expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic)),
        la_call2(SYM_Rule, expr_new_symbol(SYM_Quartics), expr_new_symbol(SYM_False)),
        la_call2(SYM_Rule, expr_new_symbol(SYM_Tolerance), expr_new_symbol(SYM_Automatic)),
    };
    symtab_set_options(SVL_HEAD, expr_new_function(expr_new_symbol(SYM_List), r, 4));
}
