/* Mathilda -- LinearSolveFunction (see linsolvefn.h for the contract).
 *
 * The factor cache is a small fixed array searched linearly by key (it holds a
 * handful of entries), evicting the least recently used. A machine matrix is
 * factored once at construction (SpFactor: UMFPACK / CHOLMOD / LAPACK, see
 * src/sparse/sparse_la_factor.c); each application loads the right-hand side,
 * runs the triangular solves, and builds a packed answer.
 */
#include "linsolvefn.h"
#include "sparse/sparse_la.h"
#include "sparse/sparse_internal.h"
#include "numarray.h"
#include "common.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include "message.h"
#include "print.h"
#include "attr.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- the factor cache ---------------------------------------------------- */

#define LSF_CACHE 8

typedef struct {
    int64_t key;
    SpFactor* f;
    uint64_t used;
} CacheEnt;

static CacheEnt g_cache[LSF_CACHE];
static uint64_t g_tick = 0;
static int64_t g_next_key = 1;
static bool g_atexit = false;

void linsolvefn_cache_clear(void) {
    for (int i = 0; i < LSF_CACHE; i++) {
        spf_free(g_cache[i].f);
        g_cache[i].f = NULL;
        g_cache[i].key = 0;
    }
}

static SpFactor* cache_get(int64_t key) {
    for (int i = 0; i < LSF_CACHE; i++)
        if (g_cache[i].f && g_cache[i].key == key) {
            g_cache[i].used = ++g_tick;
            return g_cache[i].f;
        }
    return NULL;
}

static void cache_put(int64_t key, SpFactor* f) {
    if (!g_atexit) { atexit(linsolvefn_cache_clear); g_atexit = true; }
    int victim = 0;
    for (int i = 0; i < LSF_CACHE; i++) {
        if (!g_cache[i].f) { victim = i; break; }
        if (g_cache[i].used < g_cache[victim].used) victim = i;
    }
    spf_free(g_cache[victim].f);
    g_cache[victim].key = key;
    g_cache[victim].f = f;
    g_cache[victim].used = ++g_tick;
}

/* ---- helpers ------------------------------------------------------------- */

static bool is_str(const Expr* e, const char* s) {
    return e && e->type == EXPR_STRING && strcmp(e->data.string, s) == 0;
}
static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}
static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed);
}
static void msg_with_expr(const char* head, const char* tag, const char* fmt, const Expr* e) {
    char* s = expr_to_string((Expr*)e);
    mth_message(head, tag, fmt, s ? s : "");
    free(s);
}

/* The matrix's shape, or false when m is not a matrix. */
static bool mat_dims(const Expr* m, int64_t* r, int64_t* c) {
    if (m->type == EXPR_SPARSE) {
        const SparseArray* sa = sp_of(m);
        if (sa->rank != 2) return false;
        *r = sa->dims[0];
        *c = sa->dims[1];
        return *r > 0 && *c > 0;
    }
    if (m->type == EXPR_NDARRAY) {
        if (m->data.ndarray.rank != 2) return false;
        *r = m->data.ndarray.dims[0];
        *c = m->data.ndarray.dims[1];
        return *r > 0 && *c > 0;
    }
    if (!is_list(m) || m->data.function.arg_count == 0) return false;
    size_t rows = m->data.function.arg_count;
    const Expr* r0 = m->data.function.args[0];
    if (r0->type == EXPR_NDARRAY) {
        if (r0->data.ndarray.rank != 1) return false;
        *c = r0->data.ndarray.dims[0];
    } else {
        if (!is_list(r0) || r0->data.function.arg_count == 0) return false;
        *c = (int64_t)r0->data.function.arg_count;
    }
    for (size_t i = 1; i < rows; i++) {
        const Expr* ri = m->data.function.args[i];
        int64_t ci = ri->type == EXPR_NDARRAY && ri->data.ndarray.rank == 1 ? ri->data.ndarray.dims[0]
                   : is_list(ri) ? (int64_t)ri->data.function.arg_count : -1;
        if (ci != *c) return false;
    }
    *r = (int64_t)rows;
    return true;
}

/* A machine view of m when its entries are machine numbers. */
static bool machine_view(const Expr* m, SpMat* A) {
    if (m->type == EXPR_SPARSE) return spla_view(sp_of(m), false, A);
    if (m->type != EXPR_NDARRAY) {
        CommonInexactInfo inf = common_scan_inexact(m);
        if (!inf.has_inexact || inf.min_bits > 53) return false;
    }
    return spla_from_dense(m, A);
}

/* Only Method -> Automatic / "Multifrontal" / "Cholesky" factor at
 * construction; any other option set solves afresh on application. */
static bool opts_factorable(Expr* const* opts, size_t n, SpfWant* want) {
    *want = SPF_WANT_AUTO;
    for (size_t i = 0; i < n; i++) {
        const Expr* r = opts[i];
        if (!is_rule(r)) return false;
        const Expr* l = r->data.function.args[0];
        const Expr* v = r->data.function.args[1];
        if (l->type != EXPR_SYMBOL || l->data.symbol.name != SYM_Method) return false;
        if (v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "Automatic") == 0) continue;
        if (is_str(v, "Multifrontal")) { *want = SPF_WANT_LU; continue; }
        if (is_str(v, "Cholesky")) { *want = SPF_WANT_CHOL; continue; }
        return false;
    }
    return true;
}

static SpFactor* factor_of(const Expr* m, Expr* const* opts, size_t nopts, SpfStatus* st, bool* tried) {
    *tried = false;
    *st = SPF_FAILED;
    SpfWant want;
    if (!opts_factorable(opts, nopts, &want)) return NULL;
    SpMat A;
    if (!machine_view(m, &A)) return NULL;
    *tried = true;
    if (want == SPF_WANT_AUTO && A.n <= 128) want = SPF_WANT_DENSE;
    SpFactor* f = spf_factor(&A, want, st);
    spla_view_free(&A);
    return f;
}

/* ---- construction ------------------------------------------------------- */

/* MatrixRank[m] < n for a matrix with no machine factor (exact, MPFR or
 * symbolic entries): the singularity test behind LinearSolve::sing1. */
static bool exact_rank_deficient(const Expr* m, int64_t n) {
    Expr* a[1] = { expr_copy((Expr*)m) };
    Expr* call = expr_new_function(expr_new_symbol("MatrixRank"), a, 1);
    Expr* r = evaluate(call);
    expr_free(call);
    bool def = r && r->type == EXPR_INTEGER && r->data.integer < n;
    expr_free(r);
    return def;
}

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

/* Evaluate and release a freshly built expression (evaluate borrows). */
static Expr* evaluate_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static Expr* build_object(int64_t r, int64_t c, int64_t key, const Expr* m, Expr* const* opts, size_t nopts) {
    Expr* dims[2] = { mk_int(r), mk_int(c) };
    Expr* dl = expr_new_function(expr_new_symbol(SYM_List), dims, 2);
    Expr** oa = malloc(sizeof(Expr*) * (nopts ? nopts : 1));
    if (!oa) { expr_free(dl); return NULL; }
    for (size_t i = 0; i < nopts; i++) oa[i] = expr_copy(opts[i]);
    Expr* ol = expr_new_function(expr_new_symbol(SYM_List), oa, nopts);
    free(oa);
    Expr* data[3] = { mk_int(key), expr_copy((Expr*)m), ol };
    Expr* dat = expr_new_function(expr_new_symbol(SYM_List), data, 3);
    Expr* args[2] = { dl, dat };
    return expr_new_function(expr_new_symbol(SYM_LinearSolveFunction), args, 2);
}

Expr* linsolvefn_make(const Expr* call) {
    size_t n = call->data.function.arg_count;
    if (n < 1) return NULL;
    const Expr* m = call->data.function.args[0];
    Expr* const* opts = call->data.function.args + 1;
    size_t nopts = n - 1;
    for (size_t i = 0; i < nopts; i++)
        if (!is_rule(opts[i])) return NULL;
    int64_t r, c;
    if (!mat_dims(m, &r, &c)) return NULL;
    if (r != c) {
        msg_with_expr("LinearSolve", "sqmat1",
                      "The matrix %s is not square so a factorization will not be saved.", m);
        return build_object(r, c, 0, m, opts, nopts);
    }
    SpfStatus st;
    bool tried;
    SpFactor* f = factor_of(m, opts, nopts, &st, &tried);
    if (!f) {
        if (tried && st == SPF_NOTPD) {
            msg_with_expr("LinearSolve", "npdef", "The matrix %s is not positive definite.", m);
            return expr_copy((Expr*)call);
        }
        bool singular = tried && st == SPF_SINGULAR;
        if (!tried && exact_rank_deficient(m, r))
            singular = true;          /* exact / MPFR / symbolic: the reference says so too */
        if (singular)
            msg_with_expr("LinearSolve", "sing1",
                          "The matrix %s is singular so a factorization will not be saved.", m);
        return build_object(r, c, 0, m, opts, nopts);
    }
    {
        double rc = spf_rcond(f);
        if (rc >= 0.0 && rc < 1e-16)
            msg_with_expr("LinearSolve", "luc",
                          "Result for LinearSolve of badly conditioned matrix %s may contain significant numerical errors.",
                          m);
    }
    int64_t key = g_next_key++;
    cache_put(key, f);
    return build_object(r, c, key, m, opts, nopts);
}

/* ---- application -------------------------------------------------------- */

bool linsolvefn_is(const Expr* h) {
    if (!h || h->type != EXPR_FUNCTION || h->data.function.arg_count != 2) return false;
    if (h->data.function.head->type != EXPR_SYMBOL ||
        h->data.function.head->data.symbol.name != SYM_LinearSolveFunction) return false;
    const Expr* d = h->data.function.args[0];
    const Expr* dat = h->data.function.args[1];
    return is_list(d) && d->data.function.arg_count == 2 && is_list(dat) &&
           dat->data.function.arg_count == 3 && dat->data.function.args[0]->type == EXPR_INTEGER &&
           is_list(dat->data.function.args[2]);
}

/* The system a mode asks for, as a matrix expression built from m. */
static Expr* mode_matrix(const Expr* m, char mode) {
    const char* h = mode == 'T' ? "Transpose" : mode == 'C' ? "ConjugateTranspose" : "Conjugate";
    Expr* a[1] = { expr_copy((Expr*)m) };
    return evaluate_owned(expr_new_function(expr_new_symbol(h), a, 1));
}

/* The general path: LinearSolve[m', b, opts], evaluated. */
static Expr* solve_afresh(const Expr* m, char mode, const Expr* b, const Expr* ol) {
    size_t no = ol->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (2 + no));
    if (!a) return NULL;
    a[0] = mode == 'N' ? expr_copy((Expr*)m) : mode_matrix(m, mode);
    a[1] = expr_copy((Expr*)b);
    for (size_t i = 0; i < no; i++) a[2 + i] = expr_copy(ol->data.function.args[i]);
    Expr* call = expr_new_function(expr_new_symbol("LinearSolve"), a, 2 + no);
    free(a);
    Expr* r = evaluate_owned(call);
    if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
        strcmp(r->data.function.head->data.symbol.name, "LinearSolve") == 0) {
        expr_free(r);
        return NULL;
    }
    return r;
}

static Expr* property(const Expr* head, const char* p, SpFactor* f) {
    const Expr* d = head->data.function.args[0];
    if (strcmp(p, "MatrixDimensions") == 0) return expr_copy((Expr*)d);
    if (strcmp(p, "Properties") == 0) {
        const char* names[] = { "ConditionNumber", "getL", "getPermutations", "getU",
                                "MatrixDimensions", "Properties" };
        Expr* a[6];
        for (int i = 0; i < 6; i++) a[i] = expr_new_string(names[i]);
        return expr_new_function(expr_new_symbol(SYM_List), a, 6);
    }
    if (strcmp(p, "ConditionNumber") == 0) {
        double rc = f ? spf_rcond(f) : -1.0;
        if (rc > 0.0) return expr_new_real(1.0 / rc);
        if (rc == 0.0) return expr_new_symbol("Infinity");
    }
    if (strcmp(p, "ConditionNumber") == 0 || strcmp(p, "getL") == 0 || strcmp(p, "getU") == 0 ||
        strcmp(p, "getPermutations") == 0) {
        Expr* a[1] = { expr_new_string("NotAvailable") };
        return expr_new_function(expr_new_symbol(SYM_Missing), a, 1);
    }
    return NULL;
}

Expr* linsolvefn_apply(const Expr* head, Expr* const* args, size_t nargs) {
    if (!linsolvefn_is(head) || nargs < 1 || nargs > 2) return NULL;
    const Expr* dat = head->data.function.args[1];
    int64_t key = dat->data.function.args[0]->data.integer;
    const Expr* m = dat->data.function.args[1];
    const Expr* ol = dat->data.function.args[2];
    int64_t rows = head->data.function.args[0]->data.function.args[0]->type == EXPR_INTEGER
                 ? head->data.function.args[0]->data.function.args[0]->data.integer : -1;

    SpFactor* f = key ? cache_get(key) : NULL;
    if (key && !f) {
        /* evicted: factor again from the matrix the object carries */
        SpfStatus st;
        bool tried;
        f = factor_of(m, ol->data.function.args, ol->data.function.arg_count, &st, &tried);
        if (f) cache_put(key, f);
    }

    const Expr* b = args[0];
    if (nargs == 1 && b->type == EXPR_STRING) return property(head, b->data.string, f);
    char mode = 'N';
    if (nargs == 2) {
        const Expr* s = args[1];
        if (is_str(s, "N")) mode = 'N';
        else if (is_str(s, "T")) mode = 'T';
        else if (is_str(s, "C")) mode = 'C';
        else if (is_str(s, "J")) mode = 'J';
        else {
            mth_message("LinearSolve", "notrhsm",
                        "Incorrect value of the second argument. Only the strings \"N\", \"T\", \"C\", or \"J\" are valid.");
            return NULL;
        }
    }

    /* the right-hand side's leading dimension must match */
    int64_t blen = -1;
    if (b->type == EXPR_SPARSE) blen = sp_of(b)->rank >= 1 ? sp_of(b)->dims[0] : -1;
    else if (b->type == EXPR_NDARRAY) blen = b->data.ndarray.dims[0];
    else if (is_list(b)) blen = (int64_t)b->data.function.arg_count;
    int64_t want = mode == 'N' || mode == 'J' ? rows
                 : head->data.function.args[0]->data.function.args[1]->data.integer;
    if (blen < 0) return NULL;
    if (blen != want) {
        mth_message("LinearSolveFunction", "lslc",
                    "Coefficient matrix and target vector(s) or matrix do not have the same dimensions.");
        return NULL;
    }

    if (!f) return solve_afresh(m, mode, b, ol);

    /* a machine factor: load b (complex when either side is) */
    bool cplx = spf_cplx(f);
    int64_t n = spf_n(f);
    NDPresentation pres = na_result_presentation(b);
    Expr* dense = NULL;
    const Expr* src = b;
    if (b->type == EXPR_SPARSE) { dense = sp_normal(sp_of(b)); src = dense; }
    if (!src) return NULL;
    int vn = 0, br = 0, bc = 0;
    double* buf = NULL;
    bool is_mat = false;
    bool bcplx = cplx;
    if (!na_load_vector(src, bcplx, &vn, &buf)) {
        if (!bcplx && na_load_vector(src, true, &vn, &buf)) bcplx = true;
        else if (na_load_matrix(src, bcplx, true, &br, &bc, &buf)) is_mat = true;
        else if (!bcplx && na_load_matrix(src, true, true, &br, &bc, &buf)) { is_mat = true; bcplx = true; }
    }
    expr_free(dense);
    if (!buf) return solve_afresh(m, mode, b, ol);
    if (bcplx && !cplx) {
        /* a complex right-hand side on a real factor: solve the real and
         * imaginary parts separately */
        int64_t nr = is_mat ? bc : 1;
        double* re = malloc(sizeof(double) * (size_t)n * (size_t)nr);
        double* im = malloc(sizeof(double) * (size_t)n * (size_t)nr);
        if (!re || !im) { free(re); free(im); free(buf); return NULL; }
        for (int64_t k = 0; k < n * nr; k++) { re[k] = buf[2 * k]; im[k] = buf[2 * k + 1]; }
        bool tr = mode == 'T' || mode == 'C';
        bool ok = spf_solve(f, re, re, nr, tr) && spf_solve(f, im, im, nr, tr);
        for (int64_t k = 0; k < n * nr; k++) { buf[2 * k] = re[k]; buf[2 * k + 1] = im[k]; }
        free(re); free(im);
        Expr* out = ok ? (is_mat ? na_build_matrix_as(buf, (int)n, (int)nr, true, true, pres)
                                 : na_build_vector_as(buf, (int)n, true, pres)) : NULL;
        free(buf);
        return out;
    }
    int64_t nr = is_mat ? bc : 1;
    /* conjugate modes: conj(A) x = b  <=>  A conj(x) = conj(b) */
    bool conj = cplx && (mode == 'C' || mode == 'J');
    if (conj) for (int64_t k = 0; k < n * nr; k++) buf[2 * k + 1] = -buf[2 * k + 1];
    bool ok = spf_solve(f, buf, buf, nr, mode == 'T' || mode == 'C');
    if (conj) for (int64_t k = 0; k < n * nr; k++) buf[2 * k + 1] = -buf[2 * k + 1];
    Expr* out = ok ? (is_mat ? na_build_matrix_as(buf, (int)n, (int)nr, cplx, true, pres)
                             : na_build_vector_as(buf, (int)n, cplx, pres)) : NULL;
    free(buf);
    return out;
}

void linsolvefn_init(void) {
    SymbolDef* d = symtab_get_def("LinearSolveFunction");
    d->attributes |= ATTR_PROTECTED;
}
