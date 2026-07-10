/* NDArray — see ndarray.h for the design rationale. */

#include "ndarray.h"
#include "sym_names.h"
#include "attr.h"
#include "symtab.h"
#include "common.h"
#include "arithmetic.h"   /* arith_warnings_muted() */
#include "eval.h"         /* eval_clock_get() — dedup warnings within one eval */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

bool is_ndarray(const Expr* e) {
    return e && e->type == EXPR_NDARRAY;
}

size_t ndarray_size(const Expr* a) {
    size_t n = 1;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        n *= (size_t)a->data.ndarray.dims[i];
    }
    return n;
}

/* ------------------------------- dtype helpers --------------------------- */

bool ndt_is_complex(NDType dt) {
    return dt == NDT_COMPLEX64 || dt == NDT_COMPLEX32;
}

int ndt_components(NDType dt) {
    return ndt_is_complex(dt) ? 2 : 1;
}

size_t ndt_comp_size(NDType dt) {
    return (dt == NDT_FLOAT64 || dt == NDT_COMPLEX64) ? sizeof(double)
                                                      : sizeof(float);
}

size_t ndt_elem_size(NDType dt) {
    return ndt_comp_size(dt) * (size_t)ndt_components(dt);
}

void ndt_get(const void* buf, size_t k, NDType dt, double* re, double* im) {
    switch (dt) {
        case NDT_FLOAT64:
            *re = ((const double*)buf)[k]; *im = 0.0; break;
        case NDT_FLOAT32:
            *re = (double)((const float*)buf)[k]; *im = 0.0; break;
        case NDT_COMPLEX64:
            *re = ((const double*)buf)[2 * k];
            *im = ((const double*)buf)[2 * k + 1]; break;
        case NDT_COMPLEX32:
            *re = (double)((const float*)buf)[2 * k];
            *im = (double)((const float*)buf)[2 * k + 1]; break;
    }
}

void ndt_set(void* buf, size_t k, NDType dt, double re, double im) {
    switch (dt) {
        case NDT_FLOAT64:
            ((double*)buf)[k] = re; break;
        case NDT_FLOAT32:
            ((float*)buf)[k] = (float)re; break;
        case NDT_COMPLEX64:
            ((double*)buf)[2 * k] = re;
            ((double*)buf)[2 * k + 1] = im; break;
        case NDT_COMPLEX32:
            ((float*)buf)[2 * k] = (float)re;
            ((float*)buf)[2 * k + 1] = (float)im; break;
    }
}

bool ndt_from_string(const char* s, NDType* out) {
    if (!s) return false;
    if (strcmp(s, "float64") == 0) { *out = NDT_FLOAT64; return true; }
    if (strcmp(s, "float32") == 0) { *out = NDT_FLOAT32; return true; }
    if (strcmp(s, "complex64") == 0) { *out = NDT_COMPLEX64; return true; }
    if (strcmp(s, "complex32") == 0) { *out = NDT_COMPLEX32; return true; }
    return false;
}

const char* ndt_to_string(NDType dt) {
    switch (dt) {
        case NDT_FLOAT64:   return "float64";
        case NDT_FLOAT32:   return "float32";
        case NDT_COMPLEX64: return "complex64";
        case NDT_COMPLEX32: return "complex32";
    }
    return "float64";
}

NDType ndt_promote(NDType a, NDType b) {
    bool cplx = ndt_is_complex(a) || ndt_is_complex(b);
    bool wide = ndt_comp_size(a) == sizeof(double) ||
                ndt_comp_size(b) == sizeof(double);
    if (cplx) return wide ? NDT_COMPLEX64 : NDT_COMPLEX32;
    return wide ? NDT_FLOAT64 : NDT_FLOAT32;
}

/* Move a dtype onto the complex axis, preserving its component width: float32 ->
 * complex32, float64 -> complex64 (complex dtypes unchanged). Used when a
 * complex scalar or a real->complex escape must not widen the array's precision
 * (numpy value-based casting: a float32 array meeting a complex value stays at
 * float32 components, i.e. our complex32). */
NDType ndt_as_complex(NDType dt) {
    return (ndt_comp_size(dt) == sizeof(double)) ? NDT_COMPLEX64 : NDT_COMPLEX32;
}

Expr* ndarray_element_to_expr(const Expr* a, size_t k) {
    NDType dt = a->data.ndarray.dtype;
    double re, im;
    ndt_get(a->data.ndarray.data, k, dt, &re, &im);
    if (ndt_is_complex(dt)) {
        Expr* args[2] = { expr_new_real(re), expr_new_real(im) };
        return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    }
    return expr_new_real(re);
}

/* Integer or Real scalar → double. Used for the components of a Complex[] leaf
 * as well as bare real leaves. */
static bool leaf_real(const Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL) { *out = e->data.real; return true; }
    return false;
}

/* Machine-precision-packable leaf for dtype `dt`, decoded into a (re, im) pair.
 * Real dtypes accept Integer/Real only. Complex dtypes additionally accept a
 * Complex[re, im] literal (and a bare real, im=0). BigInt/Rational/MPFR/symbolic
 * all fail packing so the caller falls back to a nested List. */
static bool leaf_to_component(const Expr* e, NDType dt, double* re, double* im) {
    *im = 0.0;
    if (leaf_real(e, re)) return true;
    if (ndt_is_complex(dt) && head_is(e, SYM_Complex) &&
        e->data.function.arg_count == 2) {
        return leaf_real(e->data.function.args[0], re) &&
               leaf_real(e->data.function.args[1], im);
    }
    return false;
}

/* Recursive rank/shape probe, mirroring linalg/util.c's get_tensor_dims but
 * over `int64_t dims[NDARRAY_MAX_RANK]` and rejecting leaves not packable at
 * dtype `dt` up front. Returns 0 for a non-List leaf, -1 for jagged/non-numeric,
 * else the rank. */
static int probe_dims(const Expr* e, int64_t* dims, NDType dt) {
    if (!head_is(e, SYM_List)) return 0;
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return -1; /* empty: nothing to pack, degrade to List */

    int sub_rank = probe_dims(e->data.function.args[0], dims + 1, dt);
    if (sub_rank == 0) {
        double r, im;
        if (!leaf_to_component(e->data.function.args[0], dt, &r, &im)) return -1;
    } else if (sub_rank < 0) {
        return -1;
    }
    for (int64_t i = 1; i < len; i++) {
        int64_t cur[NDARRAY_MAX_RANK];
        int cur_rank = probe_dims(e->data.function.args[i], cur, dt);
        if (cur_rank != sub_rank) return -1;
        if (cur_rank == 0) {
            double r, im;
            if (!leaf_to_component(e->data.function.args[i], dt, &r, &im)) return -1;
        }
        for (int j = 0; j < sub_rank; j++) {
            if (cur[j] != dims[j + 1]) return -1;
        }
    }
    return sub_rank + 1;
}

static void flatten_into(const Expr* e, void* flat, size_t* idx, NDType dt) {
    if (head_is(e, SYM_List)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_into(e->data.function.args[i], flat, idx, dt);
        }
    } else {
        double re = 0.0, im = 0.0;
        leaf_to_component(e, dt, &re, &im); /* already validated by probe_dims */
        ndt_set(flat, (*idx)++, dt, re, im);
    }
}

Expr* ndarray_from_nested_list(const Expr* list, NDType dtype) {
    if (!head_is(list, SYM_List)) return NULL;
    int64_t dims[NDARRAY_MAX_RANK];
    int rank = probe_dims(list, dims, dtype);
    if (rank <= 0 || rank >= NDARRAY_MAX_RANK) return NULL;

    size_t n = 1;
    for (int i = 0; i < rank; i++) n *= (size_t)dims[i];
    void* data = malloc(ndt_elem_size(dtype) * n);
    if (!data) return NULL;
    size_t idx = 0;
    flatten_into(list, data, &idx, dtype);

    return expr_new_ndarray(rank, dims, data, dtype); /* takes ownership of data */
}

/* Recursive rebuild of one axis-level of the nested List tree, mirroring
 * linalg/dot.c's build_tensor but reading from a flat double* rather than
 * an Expr** buffer. `idx` is the running cursor into `data`. */
static Expr* rebuild_level(const int64_t* dims, int rank, int level,
                            const void* data, NDType dt, size_t* idx) {
    if (level == rank) {
        size_t k = (*idx)++;
        double re, im;
        ndt_get(data, k, dt, &re, &im);
        if (ndt_is_complex(dt)) {
            Expr* cargs[2] = { expr_new_real(re), expr_new_real(im) };
            return expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
        }
        return expr_new_real(re);
    }
    int64_t len = dims[level];
    Expr** args = malloc(sizeof(Expr*) * (size_t)len);
    for (int64_t i = 0; i < len; i++) {
        args[i] = rebuild_level(dims, rank, level + 1, data, dt, idx);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)len);
    free(args);
    return result;
}

Expr* ndarray_to_nested_list(const Expr* a) {
    size_t idx = 0;
    return rebuild_level(a->data.ndarray.dims, a->data.ndarray.rank, 0,
                          a->data.ndarray.data, a->data.ndarray.dtype, &idx);
}

Expr* ndarray_dot2(const Expr* a, const Expr* b, bool* shape_error) {
    int rankA = a->data.ndarray.rank;
    int rankB = b->data.ndarray.rank;
    if (rankA < 1 || rankA > 2 || rankB < 1 || rankB > 2) return NULL;

    const int64_t* dimsA = a->data.ndarray.dims;
    const int64_t* dimsB = b->data.ndarray.dims;
    int64_t K = dimsA[rankA - 1];
    if (K != dimsB[0]) { *shape_error = true; return NULL; }

    int64_t Ra = (rankA == 2) ? dimsA[0] : 1;
    int64_t Sb = (rankB == 2) ? dimsB[1] : 1;

    NDType dta = a->data.ndarray.dtype, dtb = b->data.ndarray.dtype;
    NDType dtc = ndt_promote(dta, dtb);
    const void* pa = a->data.ndarray.data;
    const void* pb = b->data.ndarray.data;
    size_t count = (size_t)(Ra * Sb);
    void* out = malloc(ndt_elem_size(dtc) * count);

    if (dta == NDT_FLOAT64 && dtb == NDT_FLOAT64) {
        /* Preserved hot path: real double triple loop, no per-element widening. */
        double* od = (double*)out;
        for (int64_t r = 0; r < Ra; r++)
            for (int64_t s = 0; s < Sb; s++) {
                double sum = 0.0;
                for (int64_t k = 0; k < K; k++)
                    sum += ((const double*)pa)[r * K + k] *
                           ((const double*)pb)[k * Sb + s];
                od[r * Sb + s] = sum;
            }
    } else {
        for (int64_t r = 0; r < Ra; r++)
            for (int64_t s = 0; s < Sb; s++) {
                double sre = 0.0, sim = 0.0;
                for (int64_t k = 0; k < K; k++) {
                    double ar, ai, br, bi;
                    ndt_get(pa, (size_t)(r * K + k), dta, &ar, &ai);
                    ndt_get(pb, (size_t)(k * Sb + s), dtb, &br, &bi);
                    sre += ar * br - ai * bi;
                    sim += ar * bi + ai * br;
                }
                ndt_set(out, (size_t)(r * Sb + s), dtc, sre, sim);
            }
    }

    int rankC = rankA + rankB - 2;
    if (rankC == 0) {
        double re, im;
        ndt_get(out, 0, dtc, &re, &im);
        free(out);
        if (ndt_is_complex(dtc)) {
            Expr* cargs[2] = { expr_new_real(re), expr_new_real(im) };
            return expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
        }
        return expr_new_real(re);
    }
    int64_t dimsC[2];
    if (rankC == 1) dimsC[0] = (rankA == 2) ? Ra : Sb;
    else { dimsC[0] = Ra; dimsC[1] = Sb; }
    return expr_new_ndarray(rankC, dimsC, out, dtc); /* takes ownership of out */
}

static bool same_shape(const Expr* a, const Expr* b) {
    if (a->data.ndarray.rank != b->data.ndarray.rank) return false;
    for (int i = 0; i < a->data.ndarray.rank; i++) {
        if (a->data.ndarray.dims[i] != b->data.ndarray.dims[i]) return false;
    }
    return true;
}

/* Decode a numeric scalar operand (Integer / Real / Complex[re,im]) into a
 * (re, im) pair. Returns false for a symbolic or non-machine scalar, so a Plus /
 * Times containing such a term falls through to the generic symbolic path. */
static bool scalar_value(const Expr* e, double* re, double* im, bool* is_cplx) {
    *im = 0.0; *is_cplx = false;
    if (leaf_real(e, re)) return true;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2 &&
        leaf_real(e->data.function.args[0], re) &&
        leaf_real(e->data.function.args[1], im)) {
        *is_cplx = true;
        return true;
    }
    return false;
}

/* Elementwise Plus (is_plus) / Times over a mix of NDArray and numeric-scalar
 * operands. numpy-style scalar broadcasting: a scalar operand (Integer / Real /
 * Complex) applies to every element; at least one operand must be an NDArray and
 * all NDArray operands must share a shape. A symbolic operand, or no NDArray at
 * all, returns NULL so the caller uses the generic symbolic path. Scalars are
 * "weak": they set the complex axis (a Complex scalar promotes the result to
 * complex) but do not widen the float precision, which is taken from the NDArray
 * operands (float32 array + 1 stays float32; float32 array + I -> complex32). */
Expr* ndarray_elementwise(Expr** args, size_t n, bool is_plus) {
    if (n == 0) return NULL;

    const Expr* first_nd = NULL;
    NDType dtc = NDT_FLOAT64;
    bool have_dt = false, scalar_complex = false, all_nd_f64 = true;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type == EXPR_NDARRAY) {
            if (!first_nd) first_nd = args[i];
            else if (!same_shape(first_nd, args[i])) return NULL;
            NDType d = args[i]->data.ndarray.dtype;
            if (d != NDT_FLOAT64) all_nd_f64 = false;
            dtc = have_dt ? ndt_promote(dtc, d) : d;
            have_dt = true;
        } else {
            double re, im; bool isc;
            if (!scalar_value(args[i], &re, &im, &isc)) return NULL;
            if (isc) scalar_complex = true;
        }
    }
    if (!first_nd) return NULL;  /* no array operand — not our job */
    if (scalar_complex && !ndt_is_complex(dtc)) dtc = ndt_as_complex(dtc);

    size_t sz = ndarray_size(first_nd);
    void* out = malloc(ndt_elem_size(dtc) * sz);

    if (all_nd_f64 && dtc == NDT_FLOAT64) {
        /* Preserved hot path: real double buffers, scalars folded as doubles. */
        double* od = (double*)out;
        for (size_t k = 0; k < sz; k++) {
            double acc = is_plus ? 0.0 : 1.0;
            for (size_t i = 0; i < n; i++) {
                double v;
                if (args[i]->type == EXPR_NDARRAY)
                    v = ((const double*)args[i]->data.ndarray.data)[k];
                else { double im; bool isc; scalar_value(args[i], &v, &im, &isc); }
                acc = is_plus ? acc + v : acc * v;
            }
            od[k] = acc;
        }
    } else {
        for (size_t k = 0; k < sz; k++) {
            double are = is_plus ? 0.0 : 1.0, aim = 0.0;
            for (size_t i = 0; i < n; i++) {
                double vre, vim;
                if (args[i]->type == EXPR_NDARRAY)
                    ndt_get(args[i]->data.ndarray.data, k,
                            args[i]->data.ndarray.dtype, &vre, &vim);
                else { bool isc; scalar_value(args[i], &vre, &vim, &isc); }
                if (is_plus) { are += vre; aim += vim; }
                else {
                    double nre = are * vre - aim * vim;
                    double nim = are * vim + aim * vre;
                    are = nre; aim = nim;
                }
            }
            ndt_set(out, k, dtc, are, aim);
        }
    }
    /* takes ownership of out */
    return expr_new_ndarray(first_nd->data.ndarray.rank,
                            first_nd->data.ndarray.dims, out, dtc);
}

/* --------------------------- elementwise Power --------------------------- */

/* Complex power (br+bi i) = (ar+ai i)^(er+ei i) via polar form, avoiding
 * <complex.h>. Real, non-negative real base with integer/real exponent stays
 * on the real axis; a negative real base with a non-integer exponent yields a
 * genuinely complex value (caller promotes the result dtype accordingly). */
static void cpow_pair(double ar, double ai, double er, double ei,
                      double* outr, double* outi) {
    if (ai == 0.0 && ei == 0.0 && (ar >= 0.0 || er == floor(er))) {
        /* Stays real: pow() is exact-signed and avoids polar round-trip. */
        *outr = pow(ar, er); *outi = 0.0; return;
    }
    /* Complex base with an integer exponent: binary exponentiation gives an
     * exact result (e.g. I^2 == -1 + 0 I) without polar-form round-off. */
    if (ei == 0.0 && er == floor(er) && fabs(er) <= 1024.0) {
        long p = (long)er;
        unsigned long m = (p < 0) ? (unsigned long)(-p) : (unsigned long)p;
        double rr = 1.0, ri = 0.0, br = ar, bi = ai;
        while (m) {
            if (m & 1ul) {
                double nr = rr * br - ri * bi;
                double ni = rr * bi + ri * br;
                rr = nr; ri = ni;
            }
            double sr = br * br - bi * bi;
            double si = 2.0 * br * bi;
            br = sr; bi = si;
            m >>= 1;
        }
        if (p < 0) {
            double d = rr * rr + ri * ri;
            *outr = rr / d; *outi = -ri / d;
        } else { *outr = rr; *outi = ri; }
        return;
    }
    if (ar == 0.0 && ai == 0.0) { *outr = 0.0; *outi = 0.0; return; }
    double logr = 0.5 * log(ar * ar + ai * ai);   /* log|a|            */
    double th = atan2(ai, ar);                     /* arg(a)            */
    /* log(a) = logr + i th ; (er+ei i)*log(a) = p + i q */
    double p = er * logr - ei * th;
    double q = er * th + ei * logr;
    double mag = exp(p);
    *outr = mag * cos(q);
    *outi = mag * sin(q);
}

/* True when raising a real-dtype array (base) to `(er,ei)` can leave the real
 * axis: a complex exponent, or a non-integer exponent applied to a possibly
 * negative base. Used to decide whether the result must promote to complex. */
static bool real_power_escapes(const Expr* base, double er, double ei) {
    if (ei != 0.0) return true;
    if (er == floor(er)) return false;             /* integer power stays real */
    size_t sz = ndarray_size(base);
    NDType dt = base->data.ndarray.dtype;
    for (size_t k = 0; k < sz; k++) {
        double re, im;
        ndt_get(base->data.ndarray.data, k, dt, &re, &im);
        if (re < 0.0) return true;
    }
    return false;
}

Expr* ndarray_elementwise_power(const Expr* a, const Expr* b) {
    if (!is_ndarray(a) || !is_ndarray(b) || !same_shape(a, b)) return NULL;
    NDType dta = a->data.ndarray.dtype, dtb = b->data.ndarray.dtype;
    size_t sz = ndarray_size(a);

    /* Result promotes to complex if either operand is complex, or a real base
     * escapes the real axis at any element. */
    NDType dtc = ndt_promote(dta, dtb);
    if (!ndt_is_complex(dtc)) {
        for (size_t k = 0; k < sz; k++) {
            double ar, ai, er, ei;
            ndt_get(a->data.ndarray.data, k, dta, &ar, &ai);
            ndt_get(b->data.ndarray.data, k, dtb, &er, &ei);
            if (er != floor(er) && ar < 0.0) { dtc = ndt_as_complex(dtc); break; }
        }
    }

    void* out = malloc(ndt_elem_size(dtc) * sz);
    for (size_t k = 0; k < sz; k++) {
        double ar, ai, er, ei, rr, ri;
        ndt_get(a->data.ndarray.data, k, dta, &ar, &ai);
        ndt_get(b->data.ndarray.data, k, dtb, &er, &ei);
        cpow_pair(ar, ai, er, ei, &rr, &ri);
        ndt_set(out, k, dtc, rr, ri);
    }
    return expr_new_ndarray(a->data.ndarray.rank, a->data.ndarray.dims, out, dtc);
}

Expr* ndarray_base_scalar_power(double br, double bi, const Expr* exp_arr) {
    if (!is_ndarray(exp_arr)) return NULL;
    NDType dte = exp_arr->data.ndarray.dtype;
    size_t sz = ndarray_size(exp_arr);

    /* Result is complex if the base is complex, the exponent array is complex,
     * or a negative real base is raised to a non-integer exponent anywhere. */
    NDType dtc = dte;
    if (bi != 0.0) dtc = ndt_as_complex(dtc);
    else if (!ndt_is_complex(dte) && br < 0.0) {
        for (size_t k = 0; k < sz; k++) {
            double er, ei;
            ndt_get(exp_arr->data.ndarray.data, k, dte, &er, &ei);
            if (er != floor(er)) { dtc = ndt_as_complex(dtc); break; }
        }
    }

    void* out = malloc(ndt_elem_size(dtc) * sz);
    for (size_t k = 0; k < sz; k++) {
        double er, ei, rr, ri;
        ndt_get(exp_arr->data.ndarray.data, k, dte, &er, &ei);
        cpow_pair(br, bi, er, ei, &rr, &ri);
        ndt_set(out, k, dtc, rr, ri);
    }
    return expr_new_ndarray(exp_arr->data.ndarray.rank,
                            exp_arr->data.ndarray.dims, out, dtc);
}

Expr* ndarray_scalar_power(const Expr* a, double er, double ei) {
    if (!is_ndarray(a)) return NULL;
    NDType dta = a->data.ndarray.dtype;
    size_t sz = ndarray_size(a);

    NDType dtc = dta;
    if (ei != 0.0) dtc = ndt_as_complex(dtc);
    else if (!ndt_is_complex(dta) && real_power_escapes(a, er, ei))
        dtc = ndt_as_complex(dtc);

    void* out = malloc(ndt_elem_size(dtc) * sz);
    for (size_t k = 0; k < sz; k++) {
        double ar, ai, rr, ri;
        ndt_get(a->data.ndarray.data, k, dta, &ar, &ai);
        cpow_pair(ar, ai, er, ei, &rr, &ri);
        ndt_set(out, k, dtc, rr, ri);
    }
    return expr_new_ndarray(a->data.ndarray.rank, a->data.ndarray.dims, out, dtc);
}

/* Format an NDArray's shape as "{d0, d1, ...}" into buf. */
static void format_shape(const Expr* a, char* buf, size_t buflen) {
    size_t off = 0;
    int written = snprintf(buf + off, buflen - off, "{");
    if (written > 0) off += (size_t)written;
    for (int i = 0; i < a->data.ndarray.rank && off < buflen; i++) {
        written = snprintf(buf + off, buflen - off, "%s%lld",
                           i ? ", " : "", (long long)a->data.ndarray.dims[i]);
        if (written > 0) off += (size_t)written;
    }
    if (off < buflen) snprintf(buf + off, buflen - off, "}");
}

/* Emit `msg` to stderr at most once per evaluation. The fixed-point evaluator
 * visits the same node more than once per evaluate(), so dedup on the eval
 * clock (constant within one evaluate(), bumped between REPL inputs) plus a
 * hash of the message: an identical warning prints once per evaluation but
 * still fires on a later input. Respects the arithmetic-warning mute. */
static void ndarray_warn_once(const char* msg) {
    if (arith_warnings_muted()) return;
    static uint64_t last_clock = 0;
    static uint64_t last_key = 0;
    uint64_t key = 1469598103934665603ULL; /* FNV-1a over the text */
    for (const char* c = msg; *c; c++)
        key = (key ^ (unsigned char)*c) * 1099511628211ULL;
    uint64_t clock = eval_clock_get();
    if (clock == last_clock && key == last_key) return;
    fputs(msg, stderr);
    last_clock = clock;
    last_key = key;
}

bool ndarray_warn_shape_mismatch(Expr** args, size_t n, const char* verb) {
    if (n < 2) return false;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type != EXPR_NDARRAY) return false; /* not a pure-NDArray op */
    }
    for (size_t i = 1; i < n; i++) {
        if (!same_shape(args[0], args[i])) {
            char s0[256], si[256], msg[600];
            format_shape(args[0], s0, sizeof s0);
            format_shape(args[i], si, sizeof si);
            snprintf(msg, sizeof msg,
                "NDArray::shape: NDArray objects of shapes %s and %s cannot "
                "be %s elementwise.\n", s0, si, verb);
            ndarray_warn_once(msg);
            return true;
        }
    }
    return false;
}

bool ndarray_warn_symbolic(Expr** args, size_t n, const char* verb) {
    bool has_nd = false, has_sym = false;
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type == EXPR_NDARRAY) has_nd = true;
        else if (!expr_is_numeric_like(args[i])) has_sym = true;
    }
    if (!has_nd || !has_sym) return false;
    char msg[300];
    snprintf(msg, sizeof msg,
        "NDArray::sym: NDArray objects are purely numeric and cannot be %s "
        "with a symbolic expression; the expression is left unevaluated.\n",
        verb);
    ndarray_warn_once(msg);
    return true;
}

/* Structural shape probe independent of leaf numeric-ness: returns the depth
 * of `e` treating any non-List as a depth-0 leaf, writing the shape into dims,
 * or -1 if `e` is internally ragged (siblings of differing shape, or a mix of
 * List and non-List siblings). A bare non-List returns 0. Unlike probe_dims,
 * this does NOT reject symbolic leaves — so NDArray[{1, x}] (rectangular but
 * symbolic) is not ragged, while NDArray[{1, {2}}] is. */
static int shape_probe(const Expr* e, int64_t* dims, int max) {
    if (!head_is(e, SYM_List)) return 0;   /* leaf (any non-List) */
    if (max <= 0) return 0;                /* too deep to record; stop */
    int64_t len = (int64_t)e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return 1;                /* empty list: shape {0}, not ragged */

    int64_t sub[NDARRAY_MAX_RANK];
    int sub_rank = shape_probe(e->data.function.args[0], sub, max - 1);
    if (sub_rank < 0) return -1;           /* first child already ragged */
    for (int64_t i = 1; i < len; i++) {
        int64_t cur[NDARRAY_MAX_RANK];
        int cur_rank = shape_probe(e->data.function.args[i], cur, max - 1);
        if (cur_rank != sub_rank) return -1;               /* mixed depth */
        for (int j = 0; j < sub_rank; j++)
            if (cur[j] != sub[j]) return -1;               /* unequal dims */
    }
    for (int j = 0; j < sub_rank; j++) dims[j + 1] = sub[j];
    return sub_rank + 1;
}

/* True iff `list` is a List that is internally ragged (non-rectangular). */
static bool ndarray_is_ragged(const Expr* list) {
    if (!head_is(list, SYM_List)) return false;
    int64_t dims[NDARRAY_MAX_RANK];
    return shape_probe(list, dims, NDARRAY_MAX_RANK) < 0;
}

/* Strip a trailing `DataType -> "..."` option from res's args (rightmost wins),
 * modelled on options.c's extract_extension_option_full. Sets *dtype (default
 * float64), decrements *argc to the count of non-option args, and returns true
 * unless an explicit DataType option carried an unknown string (in which case
 * *bad is set true and the caller leaves the expression unevaluated). */
static bool extract_datatype_option(const Expr* res, size_t* argc,
                                    NDType* dtype, bool* bad) {
    *dtype = NDT_FLOAT64;
    *bad = false;
    size_t n = res->data.function.arg_count;
    bool seen = false;
    while (n > 0) {
        const Expr* opt = res->data.function.args[n - 1];
        if (opt && opt->type == EXPR_FUNCTION && opt->data.function.head
            && opt->data.function.head->type == EXPR_SYMBOL
            && (opt->data.function.head->data.symbol == SYM_Rule
                || opt->data.function.head->data.symbol == SYM_RuleDelayed)
            && opt->data.function.arg_count == 2) {
            const Expr* lhs = opt->data.function.args[0];
            const Expr* rhs = opt->data.function.args[1];
            if (lhs && lhs->type == EXPR_SYMBOL
                && lhs->data.symbol == SYM_DataType) {
                if (!seen) {
                    if (rhs && rhs->type == EXPR_STRING
                        && ndt_from_string(rhs->data.string, dtype)) {
                        /* ok */
                    } else {
                        *bad = true;
                    }
                    seen = true;
                }
                n--;
                continue;
            }
        }
        break;
    }
    *argc = n;
    return true;
}

Expr* builtin_ndarray(Expr* res) {
    size_t argc;
    NDType dtype;
    bool bad;
    extract_datatype_option(res, &argc, &dtype, &bad);
    if (bad) {
        ndarray_warn_once("NDArray::dtype: Unknown DataType; expected one of "
                          "\"float64\", \"float32\", \"complex64\", "
                          "\"complex32\".\n");
        return NULL;
    }
    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* packed = ndarray_from_nested_list(arg, dtype);
    if (packed) return packed;
    /* A ragged (non-rectangular) list can never form an array, so reject it
     * with a warning rather than leaving a valid-looking NDArray[...] around.
     * A non-ragged but non-machine-precision list (e.g. a symbolic entry)
     * stays unevaluated silently — it may become packable after evaluation. */
    if (ndarray_is_ragged(arg)) {
        ndarray_warn_once("NDArray::ragged: The list is not rectangular "
                          "(ragged) and cannot form an NDArray.\n");
    }
    return NULL; /* leave unevaluated */
}

Expr* builtin_datatype(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg)) return NULL;  /* stays symbolic on non-arrays */
    return expr_new_string(ndt_to_string(arg->data.ndarray.dtype));
}

Expr* builtin_ndarrayq(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return expr_new_symbol(is_ndarray(res->data.function.args[0]) ? SYM_True : SYM_False);
}

void ndarray_init(void) {
    symtab_add_builtin("NDArrayQ", builtin_ndarrayq);
    symtab_get_def("NDArrayQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDArrayQ",
        "NDArrayQ[expr]\n"
        "\tGives True if expr is an NDArray value, else False.");

    symtab_add_builtin("NDArray", builtin_ndarray);
    symtab_get_def("NDArray")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NDArray",
        "NDArray[nested_list]\n"
        "\tPacks a rectangular, machine-precision (Integer/Real) nested\n"
        "\tlist into a dense N-dimensional array (numpy ndarray style).\n"
        "\tVisibly distinct from List: Head, ListQ, and printing never\n"
        "\ttreat an NDArray as a List. Dimensions gives its shape,\n"
        "\tArrayDepth its rank, Length its leading-axis length. Builtins\n"
        "\tthat recognize NDArray (Dot, Plus, Times) use a fast C-level\n"
        "\tpath; results that would need a non-machine-precision entry\n"
        "\tauto-degrade to an ordinary nested List.\n"
        "NDArray[nested_list, DataType -> \"float32\"]\n"
        "\tPacks at the given element type: \"float64\" (default), \"float32\",\n"
        "\t\"complex64\", or \"complex32\". DataType[a] gives an array's type.\n"
        "\tA ragged (non-rectangular) list is rejected with an\n"
        "\tNDArray::ragged warning; an empty or non-machine-precision\n"
        "\tlist stays unevaluated.");
    /* Default option surfaced by Options[NDArray]. */
    {
        Expr* rule_args[2] = { expr_new_symbol(SYM_DataType),
                               expr_new_string("float64") };
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule), rule_args, 2);
        Expr* opts_args[1] = { rule };
        Expr* opts = expr_new_function(expr_new_symbol(SYM_List), opts_args, 1);
        symtab_set_options("NDArray", opts);  /* takes ownership */
    }

    symtab_add_builtin("DataType", builtin_datatype);
    symtab_get_def("DataType")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (centralized convention). */
}
