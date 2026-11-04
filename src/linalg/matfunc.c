/*
 * matfunc.c -- dispatchers for MatrixExp, MatrixLog, MatrixFunction and the
 * non-integer exponents of MatrixPower.
 *
 * Each head validates its arguments (the reference's messages: ::matsq for a
 * List that is not a non-empty square matrix, ::vector / ::lslc for
 * MatrixExp's second argument, ::argr / ::argt for a wrong argument count;
 * anything that is not a List at all stays unevaluated silently), classifies
 * the precision of its input with common_scan_inexact, and routes:
 *
 *   machine (<= 53 bits, a List of machine numbers or a float64/complex64
 *   buffer)      -> matfunc_machine.c, through na_load_matrix (which reads a
 *                   packed or visible NDArray directly) and back through
 *                   na_build_matrix_as with the input's presentation, so a
 *                   List, a packed List and an NDArray give the same bits;
 *   MPFR (> 53)  -> matfunc_mpfr.c at the input's precision (the answer keeps
 *                   it); a case the MPFR kernels decline (a function with no
 *                   MPFR kernel, a singular matrix) runs the exact engine on the
 *                   rationalised input and numericalises the result back;
 *   exact        -> matfunc_exact.c (Jordan form).
 *
 * The ND entry points ndla_matrixexp / ndla_matrixlog / ndla_matrixfunction
 * (declared in ndlinalg.h) are these same routines: the builtins dispatch to
 * them first when an argument is an NDArray, and Compile[] delegates to them
 * (src/compile/compile_ndtables.c), so the compiled answer is the interpreted
 * one by construction.
 */

#include "matfunc.h"
#include "matfunc_internal.h"
#include "linalg.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "ndarray.h"
#include "common.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"
#include "lapack.h"
#include "numeric.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ===================================================================== *
 *  Small helpers                                                        *
 * ===================================================================== */

static Expr* mk1(const char* h, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(h), args, 1);
}
static Expr* mk2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), args, 2);
}

static bool mf_head_is(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

static bool is_option(const Expr* e) {
    return (mf_head_is(e, SYM_Rule) || mf_head_is(e, SYM_RuleDelayed)) && e->data.function.arg_count == 2;
}

static bool is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

static void msg_matsq(const char* head, const Expr* m, int pos) {
    char* s = expr_to_string((Expr*)m);
    mth_message(head, "matsq", "Argument %s at position %d is not a non-empty square matrix.", s ? s : "?", pos);
    free(s);
}

static void msg_zegev(const char* head, const Expr* m) {
    char* s = expr_to_string((Expr*)m);
    mth_message(head, "zegev", "The matrix %s has a zero eigenvalue.", s ? s : "?");
    free(s);
}

/* Precision classes. */
enum { MC_EXACT = 0, MC_MACHINE = 1, MC_MPFR = 2 };

/* A validated square-matrix operand.  `m` is the matrix to compute on: the
 * argument itself (borrowed), or an owned replacement (an integer buffer
 * delisted for the exact engine, or a Parameter-scaled matrix). */
typedef struct {
    const Expr* m;
    Expr* owned;
    int n;
    int cls;
    long bits;
} MfOperand;

static void operand_free(MfOperand* op) { expr_free(op->owned); op->owned = NULL; }

static void operand_classify(MfOperand* op) {
    const Expr* m = op->m;
    if (is_ndarray(m)) {
        NDType dt = m->data.ndarray.dtype;
        bool fl = dt == NDT_FLOAT64 || dt == NDT_FLOAT32 || ndt_is_complex(dt);
        op->cls = fl ? MC_MACHINE : MC_EXACT;
        op->bits = 53;
        return;
    }
    CommonInexactInfo info = common_scan_inexact(m);
    if (!info.has_inexact) { op->cls = MC_EXACT; op->bits = 0; }
    else if (info.min_bits <= 53) { op->cls = MC_MACHINE; op->bits = 53; }
    else { op->cls = MC_MPFR; op->bits = info.min_bits; }
}

/* 1: a square matrix (op filled); 0: not an array at all (silent); -1: an
 * array of the wrong shape (::matsq emitted). */
static int operand_get(const Expr* arg, const char* head, int pos, MfOperand* op) {
    memset(op, 0, sizeof *op);
    int64_t dims[64];
    int rank;
    if (is_ndarray(arg)) {
        rank = arg->data.ndarray.rank;
        for (int i = 0; i < rank && i < 64; i++) dims[i] = arg->data.ndarray.dims[i];
    } else {
        if (!mf_head_is(arg, SYM_List)) return 0;
        rank = get_tensor_dims((Expr*)arg, dims);
    }
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        msg_matsq(head, arg, pos);
        return -1;
    }
    op->n = (int)dims[0];
    op->m = arg;
    if (is_ndarray(arg)) {
        NDType dt = arg->data.ndarray.dtype;
        if (!(dt == NDT_FLOAT64 || dt == NDT_FLOAT32 || ndt_is_complex(dt))) {
            op->owned = ndarray_to_nested_list(arg);       /* an integer buffer: exact */
            op->m = op->owned;
        }
    }
    operand_classify(op);
    return 1;
}

/* Replace the operand by an evaluated expression (owned). */
static void operand_replace(MfOperand* op, Expr* m) {
    expr_free(op->owned);
    op->owned = m;
    op->m = m;
    operand_classify(op);
}

/* Load a machine operand column-major; *cx tells whether it is complex.  A
 * matrix with symbolic-but-numeric entries (Pi, Sqrt[2]) is numericalised
 * first.  NULL when an entry is not numeric. */
static double* load_machine(const Expr* m, int n, bool* cx) {
    double* buf = NULL;
    int r = 0, c = 0;
    if (is_ndarray(m)) {
        *cx = ndt_is_complex(m->data.ndarray.dtype);
        if (na_load_matrix((Expr*)m, *cx, true, &r, &c, &buf) && r == n && c == n) return buf;
        free(buf);
        return NULL;
    }
    for (int pass = 0; pass < 2; pass++) {
        Expr* src = pass == 0 ? (Expr*)m : eval_and_free(mk1(SYM_N, expr_copy((Expr*)m)));
        bool ok = false;
        for (int want = 0; want < 2 && !ok; want++) {
            buf = NULL;
            if (na_load_matrix(src, want == 1, true, &r, &c, &buf) && r == n && c == n) {
                ok = true;
                *cx = want == 1;
            } else { free(buf); buf = NULL; }
        }
        if (pass == 1) expr_free(src);
        if (ok) return buf;
    }
    return NULL;
}

static Expr* build_machine(const double* F, int n, bool cx, const Expr* like) {
    return na_build_matrix_as(F, n, n, cx, true, na_result_presentation(like));
}

/* Drop the imaginary parts of a complex column-major buffer into a real one. */
static double* realify(const double* F, int n) {
    size_t nn = (size_t)n * (size_t)n;
    double* R = (double*)malloc(sizeof(double) * (nn ? nn : 1));
    if (!R) return NULL;
    for (size_t k = 0; k < nn; k++) R[k] = F[2 * k];
    return R;
}

/* The MPFR route through the exact engine: rationalise at `bits`, apply, and
 * numericalise the answer back to `bits`. */
static Expr* via_exact_numeric(const Expr* m, int n, long bits, const Expr* f,
                               MfxKind kind, const char* head) {
    Expr* r = common_rationalize_input(m, bits);
    Expr* ex = mf_exact_apply(r, n, f, kind, head);
    expr_free(r);
    if (!ex) return NULL;
    Expr* out = common_numericalize_result(ex, bits);
    expr_free(ex);
    return eval_and_free(out);
}

/* Parse the trailing options of a matrix-function call: Method -> anything is
 * accepted (one algorithm per precision class here), Parameter -> t (MatrixExp
 * only) scales the matrix.  false for an unknown option. */
static bool parse_options(const Expr* res, size_t first, bool allow_param, const Expr** param) {
    *param = NULL;
    for (size_t i = first; i < res->data.function.arg_count; i++) {
        const Expr* o = res->data.function.args[i];
        if (!is_option(o)) return false;
        const Expr* l = o->data.function.args[0];
        if (is_symbol(l, SYM_Method)) continue;
        if (allow_param && is_symbol(l, SYM_Parameter)) { *param = o->data.function.args[1]; continue; }
        return false;
    }
    return true;
}

static size_t count_positional(const Expr* res) {
    size_t k = 0;
    while (k < res->data.function.arg_count && !is_option(res->data.function.args[k])) k++;
    return k;
}

/* ===================================================================== *
 *  MatrixExp                                                            *
 * ===================================================================== */

/* Machine exp(A), with an MPFR-at-53-bits retry when the machine result
 * overflows (exp of a large-norm matrix). */
static Expr* expm_machine(const MfOperand* op) {
    bool cx = false;
    double* A = load_machine(op->m, op->n, &cx);
    if (!A) return NULL;
    size_t sz = (cx ? 2 : 1) * (size_t)op->n * (size_t)op->n;
    double* E = (double*)malloc(sizeof(double) * sz);
    Expr* r = NULL;
    if (E && mf_expm(A, op->n, cx, E) == MF_OK) r = build_machine(E, op->n, cx, op->m);
    free(A); free(E);
#ifdef USE_MPFR
    if (!r && !is_ndarray(op->m)) {
        MfMat W;
        if (mf_mat_load(op->m, op->n, 53, &W)) {
            if (mf_mpfr_expm(&W) == MF_OK) r = mf_mat_to_expr(&W, 53, !cx);
            mf_mat_free(&W);
        }
    }
#endif
    return r;
}

static Expr* expm_any(const MfOperand* op) {
    if (op->cls == MC_MACHINE) return expm_machine(op);
    if (op->cls == MC_MPFR) {
#ifdef USE_MPFR
        MfMat W;
        Expr* r = NULL;
        if (mf_mat_load(op->m, op->n, (mpfr_prec_t)op->bits, &W)) {
            bool real_in = true;
            size_t nn = (size_t)op->n * (size_t)op->n;
            for (size_t k = 0; k < nn; k++) if (!mpfr_zero_p(W.im[k])) { real_in = false; break; }
            if (mf_mpfr_expm(&W) == MF_OK) r = mf_mat_to_expr(&W, (mpfr_prec_t)op->bits, real_in);
            mf_mat_free(&W);
        }
        if (r) return r;
#endif
        Expr* f = expr_new_symbol(SYM_Exp);
        Expr* r2 = via_exact_numeric(op->m, op->n, op->bits, f, MFX_EXP, "MatrixExp");
        expr_free(f);
        return r2;
    }
    Expr* f = expr_new_symbol(SYM_Exp);
    Expr* r = mf_exact_apply(op->m, op->n, f, MFX_EXP, "MatrixExp");
    expr_free(f);
    return r;
}

static Expr* matrixexp_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t npos = count_positional(res);
    if (npos < 1 || npos > 2) {
        mth_message("MatrixExp", "argt", "MatrixExp called with %zu arguments; 1 or 2 arguments are expected.", npos);
        return NULL;
    }
    const Expr* param = NULL;
    if (!parse_options(res, npos, true, &param)) return NULL;

    MfOperand op;
    int g = operand_get(res->data.function.args[0], "MatrixExp", 1, &op);
    if (g <= 0) return NULL;
    if (param && !(param->type == EXPR_INTEGER && param->data.integer == 1))
        operand_replace(&op, eval_and_free(mk2(SYM_Times, expr_copy((Expr*)param), expr_copy((Expr*)op.m))));

    Expr* result = NULL;
    if (npos == 1) {
        result = expm_any(&op);
        operand_free(&op);
        return result;
    }

    /* MatrixExp[m, v]: v a non-empty vector of length n. */
    const Expr* v = res->data.function.args[1];
    int64_t vd[64];
    int vr = is_ndarray(v) ? v->data.ndarray.rank : (mf_head_is(v, SYM_List) ? get_tensor_dims((Expr*)v, vd) : -1);
    if (is_ndarray(v) && vr >= 1) vd[0] = v->data.ndarray.dims[0];
    if (vr != 1 || vd[0] == 0) {
        char* s = expr_to_string((Expr*)v);
        mth_message("MatrixExp", "vector", "Argument %s at position 2 is not a non-empty vector.", s ? s : "?");
        free(s);
        operand_free(&op);
        return NULL;
    }
    if (vd[0] != op.n) {
        mth_message("MatrixExp", "lslc", "Coefficient matrix and target vector(s) or matrix do not have the same dimensions.");
        operand_free(&op);
        return NULL;
    }

    /* The pair is numeric if either operand is. */
    int cls = op.cls;
    long bits = op.bits;
    if (!is_ndarray(v)) {
        CommonInexactInfo vi = common_scan_inexact(v);
        if (vi.has_inexact) {
            long vb = vi.min_bits <= 53 ? 53 : vi.min_bits;
            if (cls == MC_EXACT || vb < bits) { bits = vb; cls = vb <= 53 ? MC_MACHINE : MC_MPFR; }
        }
    } else if (cls == MC_EXACT) { cls = MC_MACHINE; bits = 53; }

    if (cls == MC_MACHINE) {
        bool cx = false;
        double* A = load_machine(op.m, op.n, &cx);
        double* b = NULL;
        int nv = 0;
        bool bcx = false;
        if (A) {
            if (na_load_vector((Expr*)v, false, &nv, &b)) bcx = false;
            else { free(b); b = NULL; if (na_load_vector((Expr*)v, true, &nv, &b)) bcx = true; else { free(b); b = NULL; } }
        }
        if (A && b && nv == op.n) {
            bool z = cx || bcx;
            double* out = (double*)malloc(sizeof(double) * (z ? 2 : 1) * (size_t)op.n);
            if (out && mf_expmv(A, op.n, cx, b, bcx, out) == MF_OK)
                result = na_build_vector_as(out, op.n, z, na_result_presentation(is_ndarray(op.m) ? op.m : v));
            free(out);
        }
        free(A); free(b);
        if (result) { operand_free(&op); return result; }
    }
    /* MPFR / exact / machine fallback: exp(m) . v */
    if (cls != op.cls) {
        operand_replace(&op, eval_and_free(common_numericalize_result(op.m, bits)));
    }
    Expr* E = expm_any(&op);
    if (E) result = eval_and_free(mk2(SYM_Dot, E, expr_copy((Expr*)v)));
    operand_free(&op);
    return result;
}

/* ===================================================================== *
 *  Schur-based machine functions (log, square root, power, general f)   *
 * ===================================================================== */

/* f through the machine Schur kernel.  *status receives the kernel code. */
static Expr* schur_machine(const MfOperand* op, int kind, mf_scalar_fn fn, void* ctx, int* status) {
    bool cx = false;
    double* A = load_machine(op->m, op->n, &cx);
    *status = MF_FAIL;
    if (!A) return NULL;
    size_t nn = (size_t)op->n * (size_t)op->n;
    double* F = (double*)malloc(sizeof(double) * 2 * nn);
    Expr* r = NULL;
    bool is_real = false;
    if (F) *status = mf_schur_fun(A, op->n, cx, kind, fn, ctx, F, &is_real);
    if (F && *status == MF_OK) {
        if (is_real) {
            double* R = realify(F, op->n);
            if (R) r = build_machine(R, op->n, false, op->m);
            free(R);
        } else {
            r = build_machine(F, op->n, true, op->m);
        }
    }
    free(A); free(F);
    return r;
}

/* The generic scalar function: f^(k)(z) by evaluating Derivative[k][f][z]. */
typedef struct {
    const Expr* f;
    Expr* deriv[64];
    int bad_k;          /* the derivative order that failed to be numeric */
} SymFnCtx;

static int fn_symbolic(void* vctx, double zr, double zi, int k, double* fr, double* fi) {
    SymFnCtx* c = (SymFnCtx*)vctx;
    if (k >= 64) return MF_FAIL;
    if (!c->deriv[k]) {
        c->deriv[k] = mf_derivative_function(c->f, k);
    }
    Expr* z = zi == 0.0 ? expr_new_real(zr)
                        : mk2(SYM_Complex, expr_new_real(zr), expr_new_real(zi));
    Expr* a[1] = { z };
    Expr* val = eval_and_free(expr_new_function(expr_copy(c->deriv[k]), a, 1));
    double re = 0.0, im = 0.0;
    bool ok = na_read_scalar(val, &re, &im);
    expr_free(val);
    if (!ok) { c->bad_k = k; return MF_NONNUM; }
    if (!isfinite(re) || !isfinite(im)) return k == 0 ? MF_ZEROEIG : MF_NOSOL;
    *fr = re; *fi = im;
    return MF_OK;
}

static void symfn_free(SymFnCtx* c) {
    for (int k = 0; k < 64; k++) expr_free(c->deriv[k]);
}

/* ===================================================================== *
 *  MatrixLog                                                            *
 * ===================================================================== */

static Expr* logm_any(const MfOperand* op, const char* head) {
    if (op->cls == MC_MACHINE) {
        int st;
        Expr* r = schur_machine(op, MF_KIND_LOG, NULL, NULL, &st);
        if (st == MF_ZEROEIG) msg_zegev(head, op->m);
        return r;
    }
    if (op->cls == MC_MPFR) {
#ifdef USE_MPFR
        MfMat W;
        if (mf_mat_load(op->m, op->n, (mpfr_prec_t)op->bits, &W)) {
            int st = mf_mpfr_logm(&W);
            Expr* r = st == MF_OK ? mf_mat_to_expr(&W, (mpfr_prec_t)op->bits, false) : NULL;
            mf_mat_free(&W);
            if (r) return r;
            if (st == MF_ZEROEIG) { msg_zegev(head, op->m); return NULL; }
        }
#endif
        Expr* f = expr_new_symbol(SYM_Log);
        Expr* r = via_exact_numeric(op->m, op->n, op->bits, f, MFX_LOG, head);
        expr_free(f);
        return r;
    }
    Expr* f = expr_new_symbol(SYM_Log);
    Expr* r = mf_exact_apply(op->m, op->n, f, MFX_LOG, head);
    expr_free(f);
    return r;
}

static Expr* matrixlog_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t npos = count_positional(res);
    if (npos != 1) {
        mth_message("MatrixLog", "argx", "MatrixLog called with %zu arguments; 1 argument is expected.", npos);
        return NULL;
    }
    const Expr* param = NULL;
    if (!parse_options(res, 1, false, &param)) return NULL;
    MfOperand op;
    if (operand_get(res->data.function.args[0], "MatrixLog", 1, &op) <= 0) return NULL;
    Expr* r = logm_any(&op, "MatrixLog");
    operand_free(&op);
    return r;
}

/* ===================================================================== *
 *  Principal power (MatrixPower[m, s] for non-integer s)                *
 * ===================================================================== */

/* Read an inexact-or-exact numeric exponent as a double pair. */
static bool read_exponent(const Expr* s, double* sr, double* si) {
    return na_read_scalar(s, sr, si) && isfinite(*sr) && isfinite(*si);
}

static Expr* pure_power(const Expr* s) {
    Expr* slot = mk1(SYM_Slot, expr_new_integer(1));
    return mk1(SYM_Function, mk2(SYM_Power, slot, expr_copy((Expr*)s)));
}

static void msg_power_nosol(double sr, double si) {
    if (sr == 0.5 && si == 0.0)
        mth_message("MatrixPower", "nosol", "Cannot compute the matrix square root.");
    else
        mth_message("MatrixPower", "nosol", "Cannot compute the matrix power.");
}

static Expr* powm_machine(const MfOperand* op, double sr, double si) {
    double s[2] = { sr, si };
    int st;
    int kind = (sr == 0.5 && si == 0.0) ? MF_KIND_SQRT : MF_KIND_FUNC;
    Expr* r = schur_machine(op, kind, mf_fn_pow, s, &st);
    if (st == MF_ZEROEIG) msg_zegev("MatrixPower", op->m);
    else if (st == MF_NOSOL) msg_power_nosol(sr, si);
    return r;
}

Expr* matfunc_matrixpower(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    MfOperand op;
    if (operand_get(res->data.function.args[0], "MatrixPower", 1, &op) <= 0) return NULL;
    const Expr* s = res->data.function.args[1];
    const Expr* v = argc == 3 ? res->data.function.args[2] : NULL;

    CommonInexactInfo si_info = common_scan_inexact(s);
    double sr = 0.0, sim = 0.0;
    bool s_num = read_exponent(s, &sr, &sim);
    Expr* result = NULL;

    if (op.cls == MC_EXACT && s_num && si_info.has_inexact) {
        /* exact matrix, inexact exponent: compute at the exponent's precision */
        long b = si_info.min_bits <= 53 ? 53 : si_info.min_bits;
        operand_replace(&op, eval_and_free(common_numericalize_result(op.m, b)));
    }

    if (op.cls == MC_MACHINE) {
        if (s_num) result = powm_machine(&op, sr, sim);
        else {
            Expr* f = pure_power(s);
            Expr* list = is_ndarray(op.m) ? ndarray_to_nested_list(op.m) : expr_copy((Expr*)op.m);
            result = mf_eigen_apply(list, f);
            expr_free(list);
            expr_free(f);
        }
    } else if (op.cls == MC_MPFR) {
        if (s_num && si_info.has_inexact && si_info.min_bits <= 53) {
            operand_replace(&op, eval_and_free(common_numericalize_result(op.m, 53)));
            result = powm_machine(&op, sr, sim);
        } else if (s_num) {
#ifdef USE_MPFR
            MfMat W;
            int st = MF_FAIL;
            if (mf_mat_load(op.m, op.n, (mpfr_prec_t)op.bits, &W)) {
                mpfr_t pr, pi;
                mpfr_inits2((mpfr_prec_t)op.bits, pr, pi, (mpfr_ptr)0);
                bool inexact = false;
                if (get_approx_mpfr(s, pr, pi, &inexact)) {
                    st = mf_mpfr_powm(&W, pr, pi);
                    if (st == MF_OK) result = mf_mat_to_expr(&W, (mpfr_prec_t)op.bits, false);
                }
                mpfr_clears(pr, pi, (mpfr_ptr)0);
                mf_mat_free(&W);
            }
            if (!result && st == MF_ZEROEIG && sr <= 0.0) msg_zegev("MatrixPower", op.m);
            else if (!result) {
                Expr* f = pure_power(s);
                result = via_exact_numeric(op.m, op.n, op.bits, f, MFX_POWER, "MatrixPower");
                expr_free(f);
            }
#endif
        }
    } else {
        Expr* f = pure_power(s);
        result = mf_exact_apply(op.m, op.n, f, MFX_POWER, "MatrixPower");
        expr_free(f);
    }
    operand_free(&op);

    if (result && v) {
        int64_t vd[64];
        int vr = is_ndarray(v) ? v->data.ndarray.rank : get_tensor_dims((Expr*)v, vd);
        if (is_ndarray(v) && vr >= 1) vd[0] = v->data.ndarray.dims[0];
        if (vr != 1 || vd[0] != (int64_t)op.n) {
            char* vs = expr_to_string((Expr*)v);
            mth_message("MatrixPower", "vecsh", "Vector %s has incompatible length for matrix of size %d.", vs ? vs : "?", op.n);
            free(vs);
            expr_free(result);
            return NULL;
        }
        result = eval_and_free(mk2(SYM_Dot, result, expr_copy((Expr*)v)));
    }
    return result;
}

/* ===================================================================== *
 *  MatrixFunction                                                       *
 * ===================================================================== */

static Expr* matrixfunction_impl(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t npos = count_positional(res);
    if (npos != 2) {
        mth_message("MatrixFunction", "argr", "MatrixFunction called with %zu argument%s; 2 arguments are expected.",
                    npos, npos == 1 ? "" : "s");
        return NULL;
    }
    const Expr* param = NULL;
    if (!parse_options(res, 2, false, &param)) return NULL;
    const Expr* f = res->data.function.args[0];
    MfOperand op;
    if (operand_get(res->data.function.args[1], "MatrixFunction", 2, &op) <= 0) return NULL;

    Expr* r = NULL;
    if (is_symbol(f, SYM_Exp)) {
        r = expm_any(&op);
    } else if (is_symbol(f, SYM_Log) && op.cls != MC_EXACT) {
        r = logm_any(&op, "MatrixFunction");
    } else if (is_symbol(f, SYM_Sqrt) && op.cls == MC_MACHINE) {
        r = powm_machine(&op, 0.5, 0.0);
    } else if (op.cls == MC_MACHINE) {
        SymFnCtx ctx;
        memset(&ctx, 0, sizeof ctx);
        ctx.f = f;
        ctx.bad_k = -1;
        int st;
        r = schur_machine(&op, MF_KIND_FUNC, fn_symbolic, &ctx, &st);
        symfn_free(&ctx);
        if (!r && st == MF_NONNUM && ctx.bad_k == 0) {
            Expr* list = is_ndarray(op.m) ? ndarray_to_nested_list(op.m) : expr_copy((Expr*)op.m);
            r = mf_eigen_apply(list, f);
            expr_free(list);
        } else if (!r && st == MF_NONNUM) {
            char* fs = expr_to_string((Expr*)f);
            mth_message("MatrixFunction", "drvnnum",
                        "Unable to compute the matrix function because the subsequent derivative of the function %s at a numeric value is not numeric.",
                        fs ? fs : "?");
            free(fs);
        } else if (!r && (st == MF_ZEROEIG || st == MF_NOSOL)) {
            char* fs = expr_to_string((Expr*)f);
            mth_message("MatrixFunction", "nosol", "Cannot compute the matrix function for the function %s.", fs ? fs : "?");
            free(fs);
        }
    } else if (op.cls == MC_MPFR) {
#ifdef USE_MPFR
        if (is_symbol(f, SYM_Sqrt)) {
            MfMat W;
            if (mf_mat_load(op.m, op.n, (mpfr_prec_t)op.bits, &W)) {
                if (mf_mpfr_sqrtm(&W) == MF_OK) r = mf_mat_to_expr(&W, (mpfr_prec_t)op.bits, false);
                mf_mat_free(&W);
            }
        }
#endif
        if (!r) r = via_exact_numeric(op.m, op.n, op.bits, f, MFX_FUNCTION, "MatrixFunction");
    } else {
        r = mf_exact_apply(op.m, op.n, f, MFX_FUNCTION, "MatrixFunction");
    }
    operand_free(&op);
    return r;
}

/* ===================================================================== *
 *  ND entry points and builtins                                         *
 * ===================================================================== */

Expr* ndla_matrixexp(Expr* res)      { return matrixexp_impl(res); }
Expr* ndla_matrixlog(Expr* res)      { return matrixlog_impl(res); }
Expr* ndla_matrixfunction(Expr* res) { return matrixfunction_impl(res); }

Expr* builtin_matrixexp(Expr* res) {
    if (linalg_call_has_ndarray(res)) return ndla_matrixexp(res);
    return matrixexp_impl(res);
}

Expr* builtin_matrixlog(Expr* res) {
    if (linalg_call_has_ndarray(res)) return ndla_matrixlog(res);
    return matrixlog_impl(res);
}

Expr* builtin_matrixfunction(Expr* res) {
    if (linalg_call_has_ndarray(res)) return ndla_matrixfunction(res);
    return matrixfunction_impl(res);
}

void matfunc_init(void) {
    symtab_add_builtin("MatrixExp", builtin_matrixexp);
    symtab_get_def("MatrixExp")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatrixLog", builtin_matrixlog);
    symtab_get_def("MatrixLog")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("MatrixFunction", builtin_matrixfunction);
    symtab_get_def("MatrixFunction")->attributes |= ATTR_PROTECTED;
}
