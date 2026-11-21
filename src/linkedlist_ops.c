/*
 * linkedlist_ops.c -- native numeric algorithms over LinkedList storage.
 *
 * Every routine here walks the chunks directly through the ll_span_at cursor:
 * no List is built, no Expr is created per element, nothing is copied to a
 * contiguous buffer. A full scan costs O(n) slot reads plus O(log n) cursor
 * calls (chunk capacities double), so it runs at the speed of the packed
 * kernels it mirrors.
 *
 * Results agree with the same operation on Normal[ll] to the bit. A machine
 * float sum reproduces ndreduce.c's summation tree exactly -- the same
 * nd_parallel_reduce split into thread ranges, the same pairwise halving down
 * to ND_PAIRWISE_BLOCK, and the same nd_sum_block_f64 leaf -- only over
 * logical indices instead of buffer offsets. A leaf block that straddles two
 * chunks (or lies in the downward-growing front run) is gathered into a
 * 128-element stack buffer first; every other leaf sums the chunk in place.
 */
#include "linkedlist.h"
#include "ndarray.h"           /* NDUnaryKernel, ndk_unary_lane_mode, nd_lane_degenerate1 */
#include "symtab.h"            /* symtab_lookup */
#include "ndreduce.h"          /* ND_PAIRWISE_BLOCK, nd_sum_block_f64, ndred_result_degenerate */
#include "ndarray_internal.h"  /* nd_parallel_reduce */
#include "checked_int.h"       /* ci_add_i64 */
#include "sym_names.h"
#include "eval.h"              /* evaluate, eval_is_inflight_unwind */
#include "attr.h"              /* ATTR_LISTABLE */
#include "arithmetic.h"        /* power_scalar_components */
#include "common.h"            /* head_is */
#include "numloop.h"           /* NumFn, NumPred */
#include "pack.h"              /* pack_enabled, pack_min_elements */

#include <float.h>
#include <gmp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------
 *  Machine float sum (LL_REAL, LL_CPLX)
 * ---------------------------------------------------------------------- */

/* Elements [base, base + count) of a REAL (w = 1) or CPLX (w = 2) list as a
 * contiguous run of doubles in logical order: the slot itself when one chunk
 * holds them in ascending order, else gathered into `tmp`. */
static const double* ll_run_f64(const LinkedList* ll, int64_t base, int64_t count,
                                size_t w, double* tmp) {
    const void* p;
    int dir;
    int64_t k = ll_span_at(ll, base, &p, &dir);
    if (dir > 0 && k >= count) return (const double*)p;
    int64_t got = 0;
    while (got < count) {
        k = ll_span_at(ll, base + got, &p, &dir);
        if (k > count - got) k = count - got;
        const double* d = (const double*)p;
        if (dir > 0) {
            memcpy(tmp + (size_t)got * w, d, (size_t)k * w * sizeof(double));
        } else {
            for (int64_t t = 0; t < k; t++)
                memcpy(tmp + (size_t)(got + t) * w, d - (ptrdiff_t)((size_t)t * w),
                       w * sizeof(double));
        }
        got += k;
    }
    return tmp;
}

/* ndreduce.c's nd_sum_strided (stride 1), over logical indices. */
static void ll_sum_range(const LinkedList* ll, bool cplx, int64_t base, int64_t count,
                         double* re, double* im) {
    if ((size_t)count <= ND_PAIRWISE_BLOCK) {
        double tmp[2 * ND_PAIRWISE_BLOCK];
        size_t w = cplx ? 2 : 1;
        const double* d = ll_run_f64(ll, base, count, w, tmp);
        if (!cplx) {
            *re = nd_sum_block_f64(d, (size_t)count);
            *im = 0.0;
        } else {
            double sr = 0.0, si = 0.0;
            for (int64_t i = 0; i < count; i++) { sr += d[2 * i]; si += d[2 * i + 1]; }
            *re = sr; *im = si;
        }
        return;
    }
    int64_t half = count / 2;
    double r1, i1, r2, i2;
    ll_sum_range(ll, cplx, base, half, &r1, &i1);
    ll_sum_range(ll, cplx, base + half, count - half, &r2, &i2);
    *re = r1 + r2; *im = i1 + i2;
}

typedef struct { const LinkedList* ll; bool cplx; } LLSumCtx;

static void ll_sum_reduce(void* c, size_t lo, size_t hi, double* slot) {
    const LLSumCtx* x = (const LLSumCtx*)c;
    ll_sum_range(x->ll, x->cplx, (int64_t)lo, (int64_t)(hi - lo), &slot[0], &slot[1]);
}

/* ------------------------------------------------------------------------
 *  Exact int64 sum (LL_I64), promoting to GMP on overflow
 * ---------------------------------------------------------------------- */

static Expr* ll_sum_i64(const LinkedList* ll) {
    int64_t n = ll_length(ll), i = 0, acc = 0;
    while (i < n) {
        const void* p;
        int dir;
        int64_t k = ll_span_at(ll, i, &p, &dir);
        const int64_t* v = (const int64_t*)p;
        for (int64_t t = 0; t < k; t++) {
            int64_t x = v[dir > 0 ? t : -t], nx;
            if (!ci_add_i64(acc, x, &nx)) { acc = nx; continue; }
            {
                /* Overflow: carry on exactly in GMP from this element. */
                mpz_t z, y;
                mth_mpz_init_set_i64(z, acc);
                mpz_init(y);
                for (int64_t j = i + t; j < n; ) {
                    const void* q;
                    int dj;
                    int64_t kj = ll_span_at(ll, j, &q, &dj);
                    const int64_t* w = (const int64_t*)q;
                    for (int64_t s = 0; s < kj; s++) {
                        mth_mpz_set_i64(y, w[dj > 0 ? s : -s]);
                        mpz_add(z, z, y);
                    }
                    j += kj;
                }
                Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(z));
                mpz_clear(z);
                mpz_clear(y);
                return r;
            }
        }
        i += k;
    }
    return expr_new_integer(acc);
}

Expr* ll_sum_machine(const LinkedList* ll) {
    switch (ll_kind(ll)) {
        case LL_I64:
            return ll_sum_i64(ll);
        case LL_REAL:
        case LL_CPLX: {
            LLSumCtx c = { ll, ll_kind(ll) == LL_CPLX };
            double slots[NDARRAY_MAX_THREADS * 2];
            int k = nd_parallel_reduce((size_t)ll_length(ll), ll_sum_reduce, &c, 2, slots);
            double re = 0.0, im = 0.0;
            for (int t = 0; t < k; t++) { re += slots[2 * t]; im += slots[2 * t + 1]; }
            Expr* r;
            if (c.cplx) {
                Expr* a[2] = { expr_new_real(re), expr_new_real(im) };
                r = expr_new_function(expr_new_symbol(SYM_Complex), a, 2);
            } else {
                r = expr_new_real(re);
            }
            /* Overflow / subnormal: decline, as the packed Total does, so the
             * caller's Plus keeps the value in range at machine precision. */
            if (ndred_result_degenerate(r)) { expr_free(r); return NULL; }
            return r;
        }
        default:
            return NULL;
    }
}

/* ------------------------------------------------------------------------
 *  Listable threading straight from the chunks into the chunks
 *
 *  Sin[ll], ll + 1., 2. ll, ll1 ll2 over machine reals. The streamed path
 *  (ll_thread_listable) gathers the list into a packed block, runs the packed
 *  kernel and copies the answer back: two extra passes over memory, which for a
 *  cheap kernel cost more than the kernel. Here the kernel reads the source
 *  slots and writes the destination slots directly, split across threads by
 *  nd_parallel_for exactly as the packed map is.
 *
 *  Bit-identical to the packed List by construction: the unary arm is
 *  ndarray.c's ndu_hot_chunk loop (the same k->real per element and the same
 *  lane check), and the binary arm is one IEEE add or multiply per element,
 *  which is all the packed Plus/Times do on float64 buffers. Any element that
 *  path would not answer at machine precision -- a kernel decline, a result
 *  outside the normal range (the numeric tower then wants MPFR) -- declines
 *  the whole call, and the streamed path answers exactly as the List would.
 * ---------------------------------------------------------------------- */


typedef struct {
    const LinkedList* a;     /* source (machine reals)                     */
    const LinkedList* b;     /* second operand list, or NULL               */
    LinkedList*       out;   /* reserved, same length                      */
    const NDUnaryKernel* k;  /* unary arm, or NULL for the binary arm      */
    int    mode;             /* the kernel's lane mode                     */
    bool   times;            /* binary: multiply, else add                 */
    double s;                /* binary: the scalar operand (when !b)       */
} LLMapCtx;

static inline bool ll_real_degenerate(double y) {
    double m = fabs(y);
    return !(m >= DBL_MIN && m <= DBL_MAX) && y != 0.0;
}

static bool ll_map_chunk(void* c, size_t lo, size_t hi) {
    const LLMapCtx* x = (const LLMapCtx*)c;
    for (int64_t i = (int64_t)lo; i < (int64_t)hi; ) {
        const void *pa, *pb = NULL, *po;
        int da, db = 1, dout;
        int64_t k = ll_span_at(x->a, i, &pa, &da);
        int64_t ko = ll_span_at(x->out, i, &po, &dout);      /* back run: dout = +1 */
        if (ko < k) k = ko;
        if (x->b) {
            int64_t kb = ll_span_at(x->b, i, &pb, &db);
            if (kb < k) k = kb;
        }
        if (k > (int64_t)hi - i) k = (int64_t)hi - i;
        const double* sa = (const double*)pa;
        const double* sb = (const double*)pb;
        double* o = (double*)po;
        if (x->k) {
            for (int64_t t = 0; t < k; t++) {
                double in = sa[da > 0 ? t : -t];
                if (!x->k->real(in, &o[t])) return false;
                if (x->mode != ND_LANE_NONE) {
                    double y = o[t], m = fabs(y);
                    if (!(m >= DBL_MIN && m <= DBL_MAX) &&
                        nd_lane_degenerate1(x->mode, in, 0.0, y, 0.0)) return false;
                }
            }
        } else {
            for (int64_t t = 0; t < k; t++) {
                double u = sa[da > 0 ? t : -t];
                double v = sb ? sb[db > 0 ? t : -t] : x->s;
                double y = x->times ? u * v : u + v;
                if (ll_real_degenerate(y)) return false;
                o[t] = y;
            }
        }
        i += k;
    }
    return true;
}

static const LinkedList* ll_real_operand(const Expr* e) {
    const LinkedList* l = expr_linkedlist(e);
    return (l && ll_kind(l) == LL_REAL) ? l : NULL;
}

/* A machine scalar for the binary arm: a machine Real, or an Integer the
 * conversion carries exactly. Exact 0 declines: Plus[x, 0] is x itself (it
 * keeps -0.) and Times[0, x] the exact 0, which the general path supplies. */
static bool ll_real_scalar(const Expr* e, double* v) {
    if (e->type == EXPR_REAL) { *v = e->data.real; return true; }
    if (e->type == EXPR_INTEGER && e->data.integer != 0 &&
        e->data.integer >= -(INT64_C(1) << 53) && e->data.integer <= (INT64_C(1) << 53)) {
        *v = (double)e->data.integer;
        return true;
    }
    return false;
}

Expr* ll_listable_direct(Expr* res) {
    const Expr* h = res->data.function.head;
    if (h->type != EXPR_SYMBOL) return NULL;
    const char* name = h->data.symbol.name;
    size_t argc = res->data.function.arg_count;
    LLMapCtx c;
    memset(&c, 0, sizeof c);
    if (argc == 1) {
        const SymbolDef* d = symtab_lookup(name);
        c.a = ll_real_operand(res->data.function.args[0]);
        if (!c.a || !d || d->down_values || !d->ndarray_unary_kernel) return NULL;
        c.k = (const NDUnaryKernel*)d->ndarray_unary_kernel;
        if (!c.k->real || c.k->to_int || c.k->to_real) return NULL;
        c.mode = ndk_unary_lane_mode(c.k);
    } else if (argc == 2 && (name == SYM_Plus || name == SYM_Times)) {
        const Expr* x = res->data.function.args[0];
        const Expr* y = res->data.function.args[1];
        c.times = name == SYM_Times;
        if ((c.a = ll_real_operand(x)) != NULL) {
            if (!(c.b = ll_real_operand(y)) && !ll_real_scalar(y, &c.s)) return NULL;
        } else if ((c.a = ll_real_operand(y)) != NULL) {
            if (!ll_real_scalar(x, &c.s)) return NULL;
        } else {
            return NULL;
        }
        if (c.b && ll_length(c.b) != ll_length(c.a)) return NULL;   /* tdlen: general path */
    } else {
        return NULL;
    }
    int64_t n = ll_length(c.a);
    LinkedList* out = ll_new();
    if (n > 0) {
        ll_reserve_typed(out, LL_REAL, n);
        c.out = out;
        if (!nd_parallel_for((size_t)n, ll_map_chunk, &c)) { ll_release(out); return NULL; }
    }
    return expr_new_linkedlist(out);
}

/* ------------------------------------------------------------------------
 *  Listable threading through the machine engines
 *
 *  Everything ll_listable_direct does not cover -- every registered unary,
 *  binary and n-ary kernel, Plus and Times of any machine kind (exact int64
 *  included), Power -- runs the very engine a packed List reaches, through its
 *  container-free core (ndarray.h, NDOperand). A list held in one ascending
 *  chunk span is handed over in place; any other is gathered into one raw
 *  scratch buffer of its slots (no List, packed List or NDArray is built), and
 *  the answer buffer becomes the new list's chunks. Same engine, same
 *  category decisions, same lane checks: bit-identical to the packed List, and
 *  declining exactly where it declines.
 * ---------------------------------------------------------------------- */


static NDType ll_dtype(LLKind k) {
    return k == LL_I64 ? NDT_INT64 : k == LL_REAL ? NDT_FLOAT64 : NDT_COMPLEX64;
}

static bool ll_kind_of_dtype(NDType dt, LLKind* k) {
    if (dt == NDT_INT64)     { *k = LL_I64;  return true; }
    if (dt == NDT_FLOAT64)   { *k = LL_REAL; return true; }
    if (dt == NDT_COMPLEX64) { *k = LL_CPLX; return true; }
    return false;
}

/* One engine operand. A typed LinkedList: in place when one ascending span
 * holds it, else gathered into *scratch (the caller frees). A rank-1 array, or
 * a scalar, as nd_operand_of gives it. False for a boxed list, a plain List or
 * a higher-rank array: the element path answers those. */
static bool ll_nd_operand(const Expr* a, NDOperand* o, void** scratch, bool* is_buf) {
    *scratch = NULL;
    *is_buf = false;
    const LinkedList* l = expr_linkedlist(a);
    if (l) {
        LLKind k = ll_kind(l);
        if (k != LL_I64 && k != LL_REAL && k != LL_CPLX) return false;
        int64_t n = ll_length(l);
        const void* p;
        int dir;
        o->dt = ll_dtype(k);
        o->sz = (size_t)n;
        o->scalar = NULL;
        if (ll_span_at(l, 0, &p, &dir) >= n && dir > 0) {
            o->data = p;
        } else {
            *scratch = malloc((size_t)n * ll_slot_size(l));
            if (!*scratch) return false;
            ll_gather(l, 0, n, *scratch);
            o->data = *scratch;
        }
        *is_buf = true;
        return true;
    }
    if (is_ndarray(a)) {
        if (a->data.ndarray.rank != 1) return false;
        *o = nd_operand_of(a);
        *is_buf = true;
        return true;
    }
    if (a->type == EXPR_FUNCTION && head_is(a, SYM_List)) return false;
    *o = nd_operand_of(a);
    return true;
}

/* The engine the evaluator would pick for `res` were its LinkedList operands
 * packed Lists: a registered kernel (eval.c, step 4b), else the Plus / Times /
 * Power builtins' buffer arithmetic. NULL declines. */
static Expr* ll_listable_engine(Expr* res, int64_t len) {
    const Expr* h = res->data.function.head;
    if (h->type != EXPR_SYMBOL) return NULL;
    const char* name = h->data.symbol.name;
    const SymbolDef* d = symtab_lookup(name);
    if (!d || d->down_values) return NULL;          /* a user rule fires first */
    size_t argc = res->data.function.arg_count;
    if (argc == 0 || argc > NDN_MAX_ARGS) return NULL;
    bool plus_times = (name == SYM_Plus || name == SYM_Times);
    bool power = (name == SYM_Power && argc == 2);
    if (!plus_times && !power &&
        !(argc == 1 && d->ndarray_unary_kernel) &&
        !(argc == 2 && d->ndarray_binary_kernel) &&
        !(argc >= 3 && d->ndarray_nary_kernel && (d->attributes & ATTR_LISTABLE)))
        return NULL;

    NDOperand ops[NDN_MAX_ARGS];
    void* scratch[NDN_MAX_ARGS];
    bool buf[NDN_MAX_ARGS];
    size_t got = 0;
    bool ok = true;
    for (; got < argc && ok; got++)
        ok = ll_nd_operand(res->data.function.args[got], &ops[got], &scratch[got], &buf[got]);
    void* out = NULL;
    NDType dto = NDT_FLOAT64;
    if (ok) {
        if (plus_times) {
            bool degenerate = false;
            ok = nd_op_elementwise(ops, argc, name == SYM_Plus, &degenerate, &out, &dto);
        } else if (power) {
            double r, i;
            if (buf[0] && buf[1])
                ok = nd_op_elementwise_power(&ops[0], &ops[1], &out, &dto);
            else if (buf[0])
                ok = power_scalar_components(ops[1].scalar, &r, &i) &&
                     nd_op_scalar_power(&ops[0], r, i, &out, &dto);
            else
                ok = power_scalar_components(ops[0].scalar, &r, &i) &&
                     nd_op_base_scalar_power(r, i, &ops[1], &out, &dto);
            (void)ndarray_take_degenerate();        /* the element path answers */
        } else if (argc == 1) {
            ok = nd_op_map_unary(&ops[0], (const NDUnaryKernel*)d->ndarray_unary_kernel, &out, &dto);
        } else if (argc == 2) {
            ok = nd_op_map_binary(&ops[0], &ops[1], (const NDBinaryKernel*)d->ndarray_binary_kernel,
                                  &out, &dto);
        } else {
            ok = nd_op_map_nary(ops, argc, (const NDNaryKernel*)d->ndarray_nary_kernel, &out, &dto);
        }
    }
    for (size_t i = 0; i < got; i++) free(scratch[i]);
    LLKind k;
    if (!ok || !ll_kind_of_dtype(dto, &k)) { free(out); return NULL; }
    LinkedList* r = ll_from_buffer(k, out, len);
    free(out);
    return expr_new_linkedlist(r);
}

static Expr* ll_thread_arg_any(const Expr* a, int64_t i);

int64_t ll_seq_length(const Expr* a) {
    if (a->type == EXPR_LINKEDLIST) return ll_length(a->data.llist);
    if (a->type == EXPR_FUNCTION && head_is(a, SYM_List)) return (int64_t)a->data.function.arg_count;
    if (is_ndarray(a) && a->data.ndarray.rank >= 1) return a->data.ndarray.dims[0];
    return -1;
}

Expr* ll_seq_element(const Expr* a, int64_t i) {
    if (a->type == EXPR_LINKEDLIST) return ll_get(a->data.llist, i);
    return ll_thread_arg_any(a, i);
}

/* Element i of a threaded operand, as List threading hands it to the head:
 * slot i of a LinkedList (read from the snapshot), part i of a List, element
 * (or row) i of an array; a scalar operand as it stands. New reference. */
static Expr* ll_thread_arg(const Expr* a, const LinkedList* snap, int64_t i) {
    if (snap) return ll_get(snap, i);
    return ll_thread_arg_any(a, i);
}

static Expr* ll_thread_arg_any(const Expr* a, int64_t i) {
    if (a->type == EXPR_FUNCTION && head_is(a, SYM_List)) return expr_copy(a->data.function.args[i]);
    if (is_ndarray(a) && a->data.ndarray.rank >= 1) {
        if (a->data.ndarray.rank == 1) return ndarray_element_to_expr(a, (size_t)i);
        Expr* ix = expr_new_integer(i + 1);
        bool degrade = false;
        Expr* row = ndarray_part(a, &ix, 1, &degrade);
        expr_free(ix);
        if (row) return row;
        Expr* nested = ndarray_to_nested_list(a);
        Expr* e = expr_copy(nested->data.function.args[i]);
        expr_free(nested);
        return e;
    }
    return expr_copy((Expr*)a);
}

/* The element path: head[a1_i, a2_i, ...] evaluated for each i, exactly the
 * calls List threading makes, read from the slots of a snapshot of each
 * LinkedList operand (the head may run user code that changes the list). */
static Expr* ll_listable_elements(Expr* res, int64_t len) {
    size_t argc = res->data.function.arg_count;
    LinkedList** snap = calloc(argc ? argc : 1, sizeof(LinkedList*));
    Expr** a = malloc(sizeof(Expr*) * (argc ? argc : 1));
    if (!snap || !a) { free(snap); free(a); return NULL; }
    for (size_t j = 0; j < argc; j++) {
        const LinkedList* l = expr_linkedlist(res->data.function.args[j]);
        if (l) snap[j] = ll_copy(l);
    }
    LinkedList* out = ll_new();
    Expr* unwind = NULL;
    for (int64_t i = 0; i < len && !unwind; i++) {
        for (size_t j = 0; j < argc; j++)
            a[j] = ll_thread_arg(res->data.function.args[j], snap[j], i);
        Expr* call = expr_new_function(expr_copy(res->data.function.head), a, argc);
        Expr* v = evaluate(call);
        expr_free(call);
        if (eval_is_inflight_unwind(v)) { unwind = v; break; }
        ll_push_back(out, v);
    }
    for (size_t j = 0; j < argc; j++) if (snap[j]) ll_release(snap[j]);
    free(snap);
    free(a);
    if (unwind) { ll_release(out); return unwind; }
    return expr_new_linkedlist(out);
}

Expr* ll_listable_threaded(Expr* res, int64_t len) {
    if (len == 0) return expr_new_linkedlist(ll_new());
    Expr* r = ll_listable_direct(res);
    if (r) return r;
    r = ll_listable_engine(res, len);
    if (r) return r;
    return ll_listable_elements(res, len);
}

/* ------------------------------------------------------------------------
 *  Map with a compiled numeric function (numloop.h, NumFn)
 * ---------------------------------------------------------------------- */

Expr* ll_map_compiled(const Expr* f, const LinkedList* src) {
    LLKind k = ll_kind(src);
    int64_t n = ll_length(src);
    if ((k != LL_REAL && k != LL_I64) || n == 0) return NULL;
    NumFn* fn = numloop_fn_compile(f);
    if (!fn) return NULL;
    /* numloop_map's rule: every result must be inexact -- the element is a
     * Real, or the body carries a Real literal. */
    if (k == LL_I64 && !numloop_fn_body_inexact(fn)) { numloop_fn_free(fn); return NULL; }
    LinkedList* out = ll_new();
    ll_reserve_typed(out, LL_REAL, n);
    bool ok = true;
    for (int64_t i = 0; ok && i < n; ) {
        /* the source span and the (reserved, ascending) destination span */
        const void *p, *q;
        int dir, dq;
        int64_t m = ll_span_at(src, i, &p, &dir);
        int64_t mq = ll_span_at(out, i, &q, &dq);
        if (mq < m) m = mq;
        double* o = (double*)q;
        for (int64_t t = 0; t < m; t++) {
            double x = k == LL_REAL ? LL_SLOT(double, p, dir, t)
                                    : (double)LL_SLOT(int64_t, p, dir, t);
            if (!numloop_fn_run(fn, x, &o[t])) { ok = false; break; }
        }
        i += m;
    }
    numloop_fn_free(fn);
    if (!ok) { ll_release(out); return NULL; }
    return expr_new_linkedlist(out);
}

/* ------------------------------------------------------------------------
 *  Select with a compiled predicate (numloop.h, NumPred) -- machine reals only,
 *  the contract numloop_pred_run carries.
 * ---------------------------------------------------------------------- */
Expr* ll_select_compiled(const Expr* crit, const LinkedList* src) {
    int64_t n = ll_length(src);
    if (ll_kind(src) != LL_REAL || n == 0) return NULL;
    NumPred* pr = numloop_pred_compile(crit);
    if (!pr) return NULL;
    LinkedList* out = ll_new();
    LLCursor c;
    const void* p;
    int dir;
    int64_t base, m;
    bool ok = true;
    ll_cursor(&c, src, 0, n);
    while (ok && (m = ll_cursor_next(&c, &p, &dir, &base)) > 0) {
        for (int64_t t = 0; t < m; t++) {
            double x = LL_SLOT(double, p, dir, t);
            bool keep;
            if (!numloop_pred_run(pr, x, &keep)) { ok = false; break; }
            if (keep) ll_push_back_real(out, x);
        }
    }
    numloop_pred_free(pr);
    if (!ok) { ll_release(out); return NULL; }
    return expr_new_linkedlist(out);
}

/* ------------------------------------------------------------------------
 *  Rebuilding walkers (linkedlist.h, ll_rebuild)
 * ---------------------------------------------------------------------- */
Expr* ll_rebuild(const LinkedList* src, Expr* (*fn)(Expr* elem, int64_t i, void* ctx),
                 void* ctx, Expr** unwind) {
    *unwind = NULL;
    LinkedList* snap = ll_copy(src);
    int64_t n = ll_length(snap);
    LinkedList* out = NULL;
    for (int64_t i = 0; i < n; i++) {
        bool own;
        Expr* el = ll_elem_view(snap, i, &own);
        Expr* r = fn(el, i, ctx);
        ll_elem_done(el, own);
        if (r && eval_is_inflight_unwind(r)) { *unwind = r; break; }
        if (r && !out) {
            /* the first change: the unchanged prefix, slot for slot */
            out = ll_new();
            ll_append_range(out, snap, 0, i);
        }
        if (out) {
            if (r) ll_push_back(out, r);
            else ll_push_back_from(out, snap, i);
        }
    }
    ll_release(snap);
    if (*unwind) { if (out) ll_release(out); return NULL; }
    return out ? expr_new_linkedlist(out) : NULL;
}

/* ------------------------------------------------------------------------
 *  Element vectors (linkedlist.h, LLElems)
 * ---------------------------------------------------------------------- */
bool ll_is_seq(const Expr* e) {
    return e->type == EXPR_LINKEDLIST || (e->type == EXPR_FUNCTION && head_is(e, SYM_List));
}

bool ll_elems_open(const Expr* seq, LLElems* out) {
    if (seq->type == EXPR_FUNCTION && head_is(seq, SYM_List)) {
        out->v = seq->data.function.args;
        out->n = seq->data.function.arg_count;
        out->owned = false;
        return true;
    }
    const LinkedList* ll = expr_linkedlist(seq);
    if (!ll) return false;
    int64_t n = ll_length(ll);
    Expr** v = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    if (!v) abort();                                  /* out of memory */
    if (ll_kind(ll) == LL_BOXED) {
        ll_gather(ll, 0, n, v);                       /* borrowed slots ... */
        for (int64_t i = 0; i < n; i++) v[i] = expr_copy(v[i]);   /* ... held */
    } else {
        for (int64_t i = 0; i < n; i++) v[i] = ll_get(ll, i);
    }
    out->v = v;
    out->n = (size_t)n;
    out->owned = true;
    return true;
}

void ll_elems_close(LLElems* e) {
    if (e->owned) {
        for (size_t i = 0; i < e->n; i++) expr_free(e->v[i]);
        free(e->v);
    }
    e->v = NULL;
    e->n = 0;
    e->owned = false;
}

Expr* ll_seq_like(const Expr* proto, Expr** out, size_t n) {
    if (proto->type != EXPR_LINKEDLIST)
        return expr_new_function(expr_new_symbol(SYM_List), out, n);
    LinkedList* r = ll_new();
    for (size_t i = 0; i < n; i++) ll_push_back(r, out[i]);
    return expr_new_linkedlist(r);
}

/* ------------------------------------------------------------------------
 *  Pairwise reductions with any leaf (linkedlist.h, ll_pairwise_reduce)
 * ---------------------------------------------------------------------- */
typedef struct { const LinkedList* ll; size_t w; LLRunLeaf leaf; void* ctx; int ncomp; } LLPairCtx;

/* ndreduce.c's strided recursion (stride 1) over logical indices. */
static void ll_pairwise_range(const LLPairCtx* x, int64_t base, int64_t count, double* out) {
    if ((size_t)count <= ND_PAIRWISE_BLOCK) {
        double tmp[2 * ND_PAIRWISE_BLOCK];
        const double* d = ll_run_f64(x->ll, base, count, x->w, tmp);
        x->leaf(d, (size_t)count, x->ctx, out);
        return;
    }
    int64_t half = count / 2;
    double a[2] = { 0.0, 0.0 }, b[2] = { 0.0, 0.0 };
    ll_pairwise_range(x, base, half, a);
    ll_pairwise_range(x, base + half, count - half, b);
    for (int c = 0; c < x->ncomp; c++) out[c] = a[c] + b[c];
}

static void ll_pairwise_chunk(void* c, size_t lo, size_t hi, double* slot) {
    const LLPairCtx* x = (const LLPairCtx*)c;
    ll_pairwise_range(x, (int64_t)lo, (int64_t)(hi - lo), slot);
}

void ll_pairwise_reduce(const LinkedList* ll, LLRunLeaf leaf, void* ctx, int ncomp,
                        double* out) {
    LLPairCtx x = { ll, ll_kind(ll) == LL_CPLX ? 2u : 1u, leaf, ctx, ncomp };
    double slots[NDARRAY_MAX_THREADS * 2];
    int k = nd_parallel_reduce((size_t)ll_length(ll), ll_pairwise_chunk, &x, ncomp, slots);
    for (int c = 0; c < ncomp; c++) out[c] = 0.0;
    for (int t = 0; t < k; t++)
        for (int c = 0; c < ncomp; c++) out[c] += slots[ncomp * t + c];
}

/* ------------------------------------------------------------------------
 *  Packed-or-plain counterparts (linkedlist.h)
 * ---------------------------------------------------------------------- */
bool ll_normal_packs(const LinkedList* ll) {
    LLKind k = ll_kind(ll);
    if (k != LL_I64 && k != LL_REAL && k != LL_CPLX) return false;
    return pack_enabled() && (size_t)ll_length(ll) >= pack_min_elements();
}

LinkedList* ll_new_boxed(Expr** v, size_t n) {
    LinkedList* ll = ll_new();
    if (n == 0) return ll;
    /* The first nonconforming element converts a list to boxed storage for
     * good, and a list that is not emptied keeps its kind: admit a symbol,
     * then the elements, then drop the symbol. */
    ll_push_back(ll, expr_new_symbol(SYM_Null));
    for (size_t i = 0; i < n; i++) ll_push_back(ll, expr_copy(v[i]));
    expr_free(ll_pop_front(ll));
    return ll;
}

/* ------------------------------------------------------------------------
 *  Element views and per-element threading (linkedlist.h)
 * ---------------------------------------------------------------------- */


bool ll_elems_of(const Expr* e, LLElemsC* out) {
    out->v = NULL; out->n = 0; out->owned = false; out->ll = false;
    if (!e) return false;
    if (e->type == EXPR_LINKEDLIST) {
        const LinkedList* ll = e->data.llist;
        int64_t n = ll_length(ll);
        out->ll = true;
        out->owned = true;
        out->n = n;
        out->v = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
        if (ll_kind(ll) == LL_BOXED) {
            ll_gather(ll, 0, n, out->v);              /* borrowed slots ... */
            for (int64_t i = 0; i < n; i++) out->v[i] = expr_copy(out->v[i]); /* ... owned */
        } else {
            for (int64_t i = 0; i < n; i++) out->v[i] = ll_get(ll, i);
        }
        return true;
    }
    if (head_is(e, SYM_List)) {
        out->v = e->data.function.args;
        out->n = (int64_t)e->data.function.arg_count;
        return true;
    }
    return false;
}

void ll_elems_release(LLElemsC* x) {
    if (x->owned) {
        for (int64_t i = 0; i < x->n; i++) expr_free(x->v[i]);
        free(x->v);
    }
    x->v = NULL; x->n = 0; x->owned = false;
}

Expr* ll_seq_view(const Expr* e, int64_t i, bool* owned) {
    *owned = false;
    if (e->type == EXPR_LINKEDLIST) return ll_elem_view(e->data.llist, i, owned);
    if (head_is(e, SYM_List)) return e->data.function.args[i];
    return NULL;
}

typedef struct { const Expr* res; size_t k; } LLMapArgCtx;

static Expr* ll_map_arg_elem(Expr* elem, int64_t i, void* ctx) {
    (void)i;
    const LLMapArgCtx* c = (const LLMapArgCtx*)ctx;
    size_t argc = c->res->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * argc);
    for (size_t j = 0; j < argc; j++)
        a[j] = expr_copy(j == c->k ? elem : c->res->data.function.args[j]);
    Expr* call = expr_new_function(expr_copy(c->res->data.function.head), a, argc);
    free(a);
    Expr* v = evaluate(call);
    expr_free(call);
    if (eval_is_inflight_unwind(v)) return v;
    if (expr_eq(v, elem)) { expr_free(v); return NULL; }   /* unchanged: keep the slot */
    return v;
}

Expr* ll_map_arg(const Expr* res, size_t k) {
    if (res->type != EXPR_FUNCTION || k >= res->data.function.arg_count) return NULL;
    const Expr* a = res->data.function.args[k];
    if (a->type != EXPR_LINKEDLIST) return NULL;
    LLMapArgCtx c = { res, k };
    Expr* unwind = NULL;
    Expr* r = ll_rebuild(a->data.llist, ll_map_arg_elem, &c, &unwind);
    if (unwind) return unwind;
    return r ? r : expr_copy((Expr*)a);
}

Expr* ll_from_list_result(Expr* list) {
    if (!list) return NULL;
    if (head_is(list, SYM_List)) {
        LinkedList* ll = ll_new();
        for (size_t i = 0; i < list->data.function.arg_count; i++)
            ll_push_back(ll, expr_copy(list->data.function.args[i]));
        expr_free(list);
        return expr_new_linkedlist(ll);
    }
    if (is_ndarray(list) && list->data.ndarray.rank == 1) {
        LinkedList* ll = ll_new();
        if (ll_append_elements(ll, list)) { expr_free(list); return expr_new_linkedlist(ll); }
        ll_release(ll);
    }
    return list;
}

bool ll_expr_has_linkedlist(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_LINKEDLIST) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (ll_expr_has_linkedlist(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ll_expr_has_linkedlist(e->data.function.args[i])) return true;
    return false;
}

typedef struct { Expr* (*fn)(const Expr*, void*); void* ctx; bool failed; } LLWalkCtx;

static Expr* ll_walk_elem(Expr* elem, int64_t i, void* ctx) {
    (void)i;
    LLWalkCtx* c = (LLWalkCtx*)ctx;
    if (c->failed) return NULL;
    Expr* r = c->fn(elem, c->ctx);
    if (!r) { c->failed = true; return NULL; }
    if (r == elem || (r->type == elem->type && r->type != EXPR_FUNCTION && expr_eq(r, elem))) {
        expr_free(r);
        return NULL;                                   /* unchanged: keep the slot */
    }
    return r;
}

Expr* ll_walk(const Expr* e, Expr* (*fn)(const Expr* el, void* ctx), void* ctx) {
    if (!e || e->type != EXPR_LINKEDLIST) return NULL;
    LLWalkCtx c = { fn, ctx, false };
    Expr* unwind = NULL;
    Expr* r = ll_rebuild(e->data.llist, ll_walk_elem, &c, &unwind);
    if (unwind) { if (r) expr_free(r); return unwind; }
    if (c.failed) { if (r) expr_free(r); return NULL; }
    return r ? r : expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------------
 *  Rectangular nests of Lists and LinkedLists (linkedlist.h)
 * ---------------------------------------------------------------------- */

static int64_t seq_len_any(const Expr* e) {
    if (e->type == EXPR_LINKEDLIST) return ll_length(e->data.llist);
    if (head_is(e, SYM_List)) return (int64_t)e->data.function.arg_count;
    return -1;
}

int ll_tensor_dims(const Expr* e, int64_t* dims, int maxrank) {
    if (maxrank <= 0) return 0;
    if (is_ndarray(e)) {
        int r = e->data.ndarray.rank;
        if (r > maxrank) r = maxrank;
        for (int k = 0; k < r; k++) dims[k] = e->data.ndarray.dims[k];
        return r;
    }
    int64_t n = seq_len_any(e);
    if (n < 0) return 0;
    dims[0] = n;
    if (n == 0) return 1;
    if (e->type == EXPR_LINKEDLIST && ll_kind(e->data.llist) != LL_BOXED) return 1;  /* numbers */
    int64_t sub[64], cur[64];
    int lim = maxrank - 1 > 63 ? 63 : maxrank - 1;
    bool own;
    Expr* el = ll_seq_view(e, 0, &own);
    int depth = ll_tensor_dims(el, sub, lim);
    ll_elem_done(el, own);
    if (depth < 0) return -1;
    for (int64_t i = 1; i < n; i++) {
        el = ll_seq_view(e, i, &own);
        int d = ll_tensor_dims(el, cur, lim);
        ll_elem_done(el, own);
        if (d != depth) return -1;
        for (int k = 0; k < depth; k++) if (cur[k] != sub[k]) return -1;
    }
    for (int k = 0; k < depth; k++) dims[k + 1] = sub[k];
    return depth + 1;
}

void ll_tensor_leaves(const Expr* e, int rank, Expr** out, size_t* idx) {
    if (rank == 0) { out[(*idx)++] = expr_copy((Expr*)e); return; }
    if (is_ndarray(e)) {
        /* a packed block below a LinkedList level: its leaves boxed */
        Expr* l = ndarray_to_nested_list(e);
        ll_tensor_leaves(l, rank, out, idx);
        expr_free(l);
        return;
    }
    int64_t n = seq_len_any(e);
    if (e->type == EXPR_LINKEDLIST && rank == 1) {
        for (int64_t i = 0; i < n; i++) out[(*idx)++] = ll_get(e->data.llist, i);
        return;
    }
    for (int64_t i = 0; i < n; i++) {
        bool own;
        Expr* el = ll_seq_view(e, i, &own);
        ll_tensor_leaves(el, rank - 1, out, idx);
        ll_elem_done(el, own);
    }
}

void ll_tensor_levels(const Expr* e, int rank, bool* lvl) {
    Expr* cur = expr_copy((Expr*)e);
    for (int k = 0; k < rank; k++) {
        lvl[k] = cur->type == EXPR_LINKEDLIST;
        int64_t n = seq_len_any(cur);
        if (n <= 0) { for (int j = k + 1; j < rank; j++) lvl[j] = false; break; }
        Expr* next = cur->type == EXPR_LINKEDLIST ? ll_get(cur->data.llist, 0)
                                                  : expr_copy(cur->data.function.args[0]);
        expr_free(cur);
        cur = next;
    }
    expr_free(cur);
}

static Expr* ll_tensor_build_rec(const int64_t* dims, int rank, Expr** leaves, size_t* idx,
                                 const bool* lvl) {
    if (rank == 0) return leaves[(*idx)++];
    int64_t n = dims[0];
    if (lvl && lvl[0]) {
        LinkedList* ll = ll_new();
        for (int64_t i = 0; i < n; i++)
            ll_push_back(ll, ll_tensor_build_rec(dims + 1, rank - 1, leaves, idx, lvl + 1));
        return expr_new_linkedlist(ll);
    }
    Expr** a = malloc(sizeof(Expr*) * (size_t)(n ? n : 1));
    for (int64_t i = 0; i < n; i++)
        a[i] = ll_tensor_build_rec(dims + 1, rank - 1, leaves, idx, lvl ? lvl + 1 : NULL);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), a, (size_t)n);
    free(a);
    return r;
}

Expr* ll_tensor_build(const int64_t* dims, int rank, Expr** leaves, const bool* lvl) {
    size_t idx = 0;
    return ll_tensor_build_rec(dims, rank, leaves, &idx, lvl);
}

bool ll_has_ll_level1(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_LINKEDLIST) return true;
    if (!head_is(e, SYM_List)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (e->data.function.args[i]->type == EXPR_LINKEDLIST) return true;
    return false;
}

bool ll_normal_is_packed(const LinkedList* ll) {
    LLKind k = ll_kind(ll);
    if (k == LL_EMPTY || k == LL_BOXED || !pack_enabled()) return false;
    return (size_t)ll_length(ll) >= pack_min_elements();
}
