/* Projection[u, v]     the projection of the vector u onto v,
 *                      ((u . Conjugate[v]) / (v . Conjugate[v])) v
 * Projection[u, v, f]  with the inner product f: f[u, v] / f[v, v] v
 *
 * Reference semantics (probed): a zero v gives the zero vector; vectors of
 * different length, or a vector beside a non-vector, give Projection::prnv and
 * the call stays unevaluated; two symbols evaluate to the formal
 * (y Conjugate[y] . x) / Conjugate[y] . y; a SparseArray v gives a SparseArray
 * (the answer is a multiple of v), a SparseArray u beside a dense v a dense
 * vector. Exact, symbolic and MPFR entries keep their exactness / precision
 * through the generic formula.
 *
 * Machine vectors (Lists, packed lists, NDArrays) take ONE kernel,
 * la_projection_machine, from every surface -- the boxed List path, the packed
 * and NDArray path (ndla_projection) and Compile (ND_FN2S row) -- so their
 * answers agree bit for bit.
 */
#include "linalg.h"
#include "ndlinalg.h"
#include "numarray.h"
#include "ndarray.h"
#include "common.h"
#include "sparse/sparse.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "message.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* out = (<u, v> / <v, v>) v with <x, y> = sum x_i conj(y_i); w = 1 (real) or 2
 * (interleaved complex). A zero v gives the zero vector. */
void la_projection_machine(const double* u, const double* v, int n, bool cplx, double* out) {
    if (!cplx) {
        double uv = 0.0, vv = 0.0;
        for (int i = 0; i < n; i++) { uv += u[i] * v[i]; vv += v[i] * v[i]; }
        double c = vv == 0.0 ? 0.0 : uv / vv;
        for (int i = 0; i < n; i++) out[i] = c * v[i];
        return;
    }
    double ur = 0.0, ui = 0.0, vv = 0.0;
    for (int i = 0; i < n; i++) {
        double ar = u[2 * i], ai = u[2 * i + 1], br = v[2 * i], bi = v[2 * i + 1];
        ur += ar * br + ai * bi;          /* u conj(v) */
        ui += ai * br - ar * bi;
        vv += br * br + bi * bi;
    }
    double cr = vv == 0.0 ? 0.0 : ur / vv, ci = vv == 0.0 ? 0.0 : ui / vv;
    for (int i = 0; i < n; i++) {
        double br = v[2 * i], bi = v[2 * i + 1];
        out[2 * i] = cr * br - ci * bi;
        out[2 * i + 1] = cr * bi + ci * br;
    }
}

/* The machine path for two machine vectors (any surface); NULL to decline. */
static Expr* projection_machine(const Expr* u, const Expr* v) {
    if (u->type == EXPR_SPARSE || v->type == EXPR_SPARSE) return NULL;
    if (!is_ndarray(u) && !is_ndarray(v)) {
        CommonInexactInfo iu = common_scan_inexact(u), iv = common_scan_inexact(v);
        bool machine = (iu.has_inexact && iu.min_bits <= 53) || (iv.has_inexact && iv.min_bits <= 53);
        bool bignum = (iu.has_inexact && iu.min_bits > 53) || (iv.has_inexact && iv.min_bits > 53);
        if (!machine || bignum) return NULL;
    }
    int nu = 0, nv = 0;
    double *bu = NULL, *bv = NULL;
    bool cplx = false;
    if (!na_load_vector(u, false, &nu, &bu) || !na_load_vector(v, false, &nv, &bv)) {
        free(bu); free(bv);
        bu = bv = NULL;
        if (!na_load_vector(u, true, &nu, &bu) || !na_load_vector(v, true, &nv, &bv)) {
            free(bu); free(bv);
            return NULL;
        }
        cplx = true;
    }
    if (nu != nv) { free(bu); free(bv); return NULL; }
    double* out = malloc(sizeof(double) * (cplx ? 2 : 1) * (size_t)(nu > 0 ? nu : 1));
    if (!out) { free(bu); free(bv); return NULL; }
    la_projection_machine(bu, bv, nu, cplx, out);
    const Expr* src = is_ndarray(v) ? v : u;
    Expr* r = na_build_vector_as(out, nu, cplx, na_result_presentation(src));
    free(bu); free(bv); free(out);
    return r;
}

Expr* ndla_projection(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return projection_machine(res->data.function.args[0], res->data.function.args[1]);
}

/* A vector's length: List, packed / NDArray rank 1, SparseArray rank 1; -1
 * for anything else. */
static int64_t vec_len(const Expr* e) {
    if (e->type == EXPR_SPARSE) return sp_of(e)->rank == 1 ? sp_of(e)->dims[0] : -1;
    if (e->type == EXPR_NDARRAY) return e->data.ndarray.rank == 1 ? e->data.ndarray.dims[0] : -1;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            const Expr* a = e->data.function.args[i];
            if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
                a->data.function.head->data.symbol.name == SYM_List) return -1;
            if (a->type == EXPR_NDARRAY || a->type == EXPR_SPARSE) return -1;
        }
        return (int64_t)e->data.function.arg_count;
    }
    return -1;
}

static Expr* call2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), args, 2);
}
static Expr* call1(const char* h, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(h), args, 1);
}
static Expr* eval_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}
static bool is_exact_zero(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer == 0) ||
           (e->type == EXPR_REAL && e->data.real == 0.0);
}

Expr* builtin_projection(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na < 2 || na > 3) return NULL;
    Expr* u = res->data.function.args[0];
    Expr* v = res->data.function.args[1];
    if (linalg_call_has_ndarray(res) && na == 2) return ndla_projection(res);

    int64_t lu = vec_len(u), lv = vec_len(v);
    bool atoms = lu < 0 && lv < 0 && u->type != EXPR_FUNCTION && v->type != EXPR_FUNCTION &&
                 u->type != EXPR_NDARRAY && v->type != EXPR_NDARRAY &&
                 u->type != EXPR_SPARSE && v->type != EXPR_SPARSE;
    if (!atoms && (lu < 0 || lv < 0 || lu != lv)) {
        mth_message("Projection", "prnv",
                    "The first or second argument or both are not vectors, or they are not vectors of equal length.");
        return NULL;
    }
    if (na == 2) {
        Expr* m = projection_machine(u, v);
        if (m) return m;
    }
    /* the generic formula: f[u, v] / f[v, v] v, with the default
     * f = #1 . Conjugate[#2] & (written Conjugate[v] . u for symbols, the
     * reference's form) */
    Expr *num, *den;
    if (na == 3) {
        Expr* f = res->data.function.args[2];
        num = eval_owned(expr_new_function(expr_copy(f), (Expr*[]){ expr_copy(u), expr_copy(v) }, 2));
        den = eval_owned(expr_new_function(expr_copy(f), (Expr*[]){ expr_copy(v), expr_copy(v) }, 2));
    } else if (atoms) {
        num = eval_owned(call2("Dot", call1("Conjugate", expr_copy(v)), expr_copy(u)));
        den = eval_owned(call2("Dot", call1("Conjugate", expr_copy(v)), expr_copy(v)));
    } else {
        num = eval_owned(call2("Dot", expr_copy(u), call1("Conjugate", expr_copy(v))));
        den = eval_owned(call2("Dot", expr_copy(v), call1("Conjugate", expr_copy(v))));
    }
    if (!num || !den) { expr_free(num); expr_free(den); return NULL; }
    Expr* coef;
    if (is_exact_zero(den)) {
        expr_free(num);
        coef = den;                                   /* a zero v: 0 v */
    } else {
        coef = call2("Times", num, call2("Power", den, expr_new_integer(-1)));
    }
    return eval_owned(call2("Times", atoms ? expr_copy(v) : coef, atoms ? coef : expr_copy(v)));
}
