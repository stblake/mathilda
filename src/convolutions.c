/* Mathilda — ListConvolve / ListCorrelate.
 *
 * See convolutions.h for the high-level description.
 *
 * With kernel K_r and list a_s:
 *   ListConvolve  computes  Sum_r K_r a_{s-r}
 *   ListCorrelate computes  Sum_r K_r a_{s+r}
 * over the alignment window fixed by the overhang parameters {kL, kR}.
 *
 * Alignment (derived to match the Wolfram Language). Per axis, normalise kL,kR
 * to positive kernel indices KL,KR in 1..m (a negative k maps to m+1+k), with a
 * single integer k meaning {k,k}. For output element t (1-based, 1..L) and
 * kernel element r (1-based, 1..m) the list index touched is
 *
 *   correlate: j = (t - KL) + r,   L = n + KL - KR
 *   convolve:  j = (t + KL) - r,   L = n - KL + KR
 *
 * result[t] = h( g(K_r, listval(j)) for r=1..m ), default g=Times, h=Plus, with
 * the h-arguments ordered by ascending j. listval(j) resolves out-of-range j via
 * the padding: cyclic list (default), a constant, a cyclic pad list, or empty
 * (the missing list factor is dropped, giving a single-argument g term).
 *
 * A fully general direct engine handles every case (symbolic / exact / numeric,
 * every padding and overhang, generalized g/h, and n dimensions). For large
 * numeric inputs with the default Times/Plus a separable FFT fast path is used
 * instead — FFTW for machine precision and the MPFR FFT for arbitrary precision,
 * in both 1-D and n-D (see fourier.h for the shared primitives). The fast path
 * materialises the padded list over the exact index window it needs, which
 * reduces every padding mode to a plain linear convolution computed by a
 * zero-padded FFT product; it then slices out the L outputs. */

#include "convolutions.h"
#include "sym_names.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"
#include "arithmetic.h"
#include "numeric.h"
#include "fourier.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifdef USE_MPFR
#include "numeric_complex.h"
#include <mpfr.h>
#endif

#define CONV_MAX_RANK 8
/* Below this much direct work (prod(L) * prod(m)) small numeric inputs go
 * through the direct engine — the FFT setup/roundoff is not worth it. */
#define CONV_FFT_MIN_WORK 4096

/* ==================================================================== *
 *  Small helpers
 * ==================================================================== */

static int64_t floormod(int64_t a, int64_t m) {
    int64_t r = a % m;
    if (r < 0) r += m;
    return r;
}

/* Build Head[args...] and free the args array (elements are adopted). If
 * `custom` is non-NULL it is copied as the head, else `def_sym` is used. */
static Expr* mk_head_fn(const char* def_sym, const Expr* custom, Expr** args, size_t n) {
    Expr* head = custom ? expr_copy((Expr*)custom) : expr_new_symbol(def_sym);
    Expr* r = expr_new_function(head, args, n);
    free(args);
    return r;
}

static void row_major_strides(const int64_t* dims, int rank, int64_t* strides) {
    strides[rank - 1] = 1;
    for (int i = rank - 2; i >= 0; i--) strides[i] = strides[i + 1] * dims[i + 1];
}

/* Row-major rebuild of a nested list from row-major leaves (ownership in). */
static Expr* build_nested(Expr** leaves, const int64_t* dims, int rank,
                          int level, size_t* idx);

/* Extract one numeric leaf as a machine complex; false for anything non-numeric. */
static bool leaf_to_cdouble(const Expr* e, double* re, double* im) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *re = (double)e->data.integer; *im = 0.0; return true;
        case EXPR_REAL:    *re = e->data.real;            *im = 0.0; return true;
        case EXPR_BIGINT:  *re = mpz_get_d(e->data.bigint); *im = 0.0; return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *re = mpfr_get_d(e->data.mpfr, MPFR_RNDN); *im = 0.0; return true;
#endif
        default: break;
    }
    int64_t rn, rd;
    if (is_rational(e, &rn, &rd)) { *re = (double)rn / (double)rd; *im = 0.0; return true; }
    Expr *cr, *ci;
    if (is_complex((Expr*)e, &cr, &ci)) {
        double r1, i1, r2, i2;
        if (!leaf_to_cdouble(cr, &r1, &i1)) return false;
        if (!leaf_to_cdouble(ci, &r2, &i2)) return false;
        *re = r1 - i2;
        *im = i1 + r2;
        return true;
    }
    return false;
}

/* ==================================================================== *
 *  Rectangular nested-list shape
 * ==================================================================== */

static bool nest_verify(const Expr* e, const int64_t* dims, int level, int rank) {
    if (level == rank) return !head_is(e, SYM_List);
    if (!head_is(e, SYM_List)) return false;
    if ((int64_t)e->data.function.arg_count != dims[level]) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!nest_verify(e->data.function.args[i], dims, level + 1, rank)) return false;
    return true;
}

/* Fill dims for a rectangular nested List of scalar leaves. false if ragged,
 * not a List, or deeper than CONV_MAX_RANK. */
static bool nest_dims(const Expr* e, int64_t* dims, int* rank_out) {
    int rank = 0;
    const Expr* cur = e;
    while (head_is(cur, SYM_List)) {
        if (rank >= CONV_MAX_RANK) return false;
        dims[rank++] = (int64_t)cur->data.function.arg_count;
        if (cur->data.function.arg_count == 0) { *rank_out = rank; return true; }
        cur = cur->data.function.args[0];
    }
    if (rank == 0) return false;
    *rank_out = rank;
    return nest_verify(e, dims, 0, rank);
}

static void nest_flatten(const Expr* e, int level, int rank, Expr** out, size_t* idx) {
    if (level == rank) { out[(*idx)++] = (Expr*)e; return; }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nest_flatten(e->data.function.args[i], level + 1, rank, out, idx);
}

/* ==================================================================== *
 *  Convolution specification
 * ==================================================================== */

typedef enum { PAD_LIST, PAD_CONST, PAD_EMPTY } PadKind;

typedef struct {
    ConvMode mode;
    int      rank;
    int64_t  kdims[CONV_MAX_RANK];   /* kernel dims */
    int64_t  ldims[CONV_MAX_RANK];   /* list dims */
    int64_t  Ldims[CONV_MAX_RANK];   /* output dims */
    int64_t  KL[CONV_MAX_RANK];      /* normalised positive kernel indices */
    int64_t  KR[CONV_MAX_RANK];
    int64_t  kstride[CONV_MAX_RANK];
    int64_t  lstride[CONV_MAX_RANK];
    Expr**   ker_leaves;             /* borrowed, row-major */
    Expr**   list_leaves;            /* borrowed, row-major */

    PadKind  pad_kind;
    bool     cyclic_default;         /* padding absent => list treated cyclic */
    Expr*    pad_scalar;             /* PAD_CONST (borrowed) */
    int64_t  pdims[CONV_MAX_RANK];   /* PAD_LIST */
    int64_t  pstride[CONV_MAX_RANK];
    Expr**   pad_leaves;             /* PAD_LIST (borrowed) */

    const Expr* g_head;              /* NULL => Times */
    const Expr* h_head;              /* NULL => Plus */
} ConvSpec;

/* Per-axis list index touched by output t and kernel r (1-based). */
static int64_t axis_j(const ConvSpec* s, int ax, int64_t t, int64_t r) {
    return (s->mode == CONV_MODE_CORRELATE) ? (t - s->KL[ax]) + r
                                            : (t + s->KL[ax]) - r;
}

/* Resolve the list value at multi-index jvec. Returns a borrowed Expr* leaf, or
 * NULL with *dropped = true when the position is out of range under empty
 * padding (the list factor is dropped from the term). */
static Expr* listval_leaf(const ConvSpec* s, const int64_t* jvec, bool* dropped) {
    *dropped = false;
    bool inrange = true;
    for (int ax = 0; ax < s->rank; ax++)
        if (jvec[ax] < 1 || jvec[ax] > s->ldims[ax]) { inrange = false; break; }
    if (inrange) {
        int64_t flat = 0;
        for (int ax = 0; ax < s->rank; ax++) flat += (jvec[ax] - 1) * s->lstride[ax];
        return s->list_leaves[flat];
    }
    if (s->pad_kind == PAD_EMPTY) { *dropped = true; return NULL; }
    if (s->pad_kind == PAD_CONST) return s->pad_scalar;
    /* PAD_LIST */
    int64_t flat = 0;
    for (int ax = 0; ax < s->rank; ax++) {
        int64_t P = s->pdims[ax], j = jvec[ax], n = s->ldims[ax], idx;
        if (s->cyclic_default || s->rank > 1) {
            idx = floormod(j - 1, P);
        } else { /* explicit 1-D pad list: independent left/right tilings */
            idx = (j > n) ? floormod(j - n - 1, P) : floormod(j - 1, P);
        }
        flat += idx * s->pstride[ax];
    }
    return s->pad_leaves[flat];
}

/* ==================================================================== *
 *  Direct engine (all cases)
 * ==================================================================== */

/* Build one g-term for a fully chosen kernel multi-index rvec at output tvec. */
static Expr* make_term(const ConvSpec* s, const int64_t* tvec, const int64_t* rvec) {
    int64_t jvec[CONV_MAX_RANK];
    int64_t kflat = 0;
    for (int ax = 0; ax < s->rank; ax++) {
        jvec[ax] = axis_j(s, ax, tvec[ax], rvec[ax]);
        kflat += (rvec[ax] - 1) * s->kstride[ax];
    }
    Expr* kleaf = s->ker_leaves[kflat];
    bool dropped;
    Expr* val = listval_leaf(s, jvec, &dropped);
    if (dropped) {
        Expr** ga = malloc(sizeof(Expr*));
        ga[0] = expr_copy(kleaf);
        return mk_head_fn(SYM_Times, s->g_head, ga, 1);
    }
    Expr** ga = malloc(2 * sizeof(Expr*));
    ga[0] = expr_copy(kleaf);
    ga[1] = expr_copy(val);
    return mk_head_fn(SYM_Times, s->g_head, ga, 2);
}

/* Combine terms with h, nested by kernel axis, ordered by ascending list index
 * j on each axis (ascending r for correlate, descending r for convolve). */
static Expr* combine_axis(const ConvSpec* s, const int64_t* tvec, int kax, int64_t* rvec) {
    if (kax == s->rank) return make_term(s, tvec, rvec);
    int64_t m = s->kdims[kax];
    Expr** terms = malloc((size_t)m * sizeof(Expr*));
    for (int64_t idx = 0; idx < m; idx++) {
        int64_t r = (s->mode == CONV_MODE_CORRELATE) ? (idx + 1) : (m - idx);
        rvec[kax] = r;
        terms[idx] = combine_axis(s, tvec, kax + 1, rvec);
    }
    return mk_head_fn(SYM_Plus, s->h_head, terms, (size_t)m);
}

/* Recurse over the output multi-index, building the nested result list. */
static Expr* build_output(const ConvSpec* s, int oax, int64_t* tvec) {
    if (oax == s->rank) {
        int64_t rvec[CONV_MAX_RANK];
        return eval_and_free(combine_axis(s, tvec, 0, rvec));
    }
    int64_t L = s->Ldims[oax];
    Expr** els = malloc((size_t)L * sizeof(Expr*));
    for (int64_t t = 1; t <= L; t++) {
        tvec[oax] = t;
        els[t - 1] = build_output(s, oax + 1, tvec);
    }
    return mk_head_fn(SYM_List, NULL, els, (size_t)L);
}

static Expr* conv_direct(const ConvSpec* s) {
    int64_t tvec[CONV_MAX_RANK];
    return build_output(s, 0, tvec);
}

/* ==================================================================== *
 *  Direct numeric engine (default Times/Plus, inexact numeric data)
 *
 *  A tight O(L*m) complex sum with no Expr allocation or evaluator round
 *  trip — the fast path for numeric inputs below the FFT threshold. Returns
 *  NULL (fall back to the general engine) if any leaf is not directly
 *  representable as a machine complex (e.g. a symbolic constant like Pi).
 * ==================================================================== */

static Expr* num_leaf(double re, double im, bool any_complex) {
    if (!any_complex) return expr_new_real(re);
    if (fabs(im) <= 16.0 * DBL_EPSILON * (fabs(re) + fabs(im))) return expr_new_real(re);
    return make_complex(expr_new_real(re), expr_new_real(im));
}

static Expr* conv_direct_machine(const ConvSpec* s, bool any_complex) {
    int64_t Ltot = 1, Lstr[CONV_MAX_RANK];
    for (int ax = 0; ax < s->rank; ax++) Ltot *= s->Ldims[ax];
    row_major_strides(s->Ldims, s->rank, Lstr);
    int64_t mtot = 1;
    for (int ax = 0; ax < s->rank; ax++) mtot *= s->kdims[ax];

    Expr** leaves = malloc((size_t)Ltot * sizeof(Expr*));
    for (int64_t o = 0; o < Ltot; o++) {
        int64_t tvec[CONV_MAX_RANK];
        for (int ax = 0; ax < s->rank; ax++) tvec[ax] = (o / Lstr[ax]) % s->Ldims[ax] + 1;
        double ar = 0.0, ai = 0.0;
        for (int64_t kf = 0; kf < mtot; kf++) {
            int64_t jvec[CONV_MAX_RANK];
            for (int ax = 0; ax < s->rank; ax++) {
                int64_t r = (kf / s->kstride[ax]) % s->kdims[ax] + 1;
                jvec[ax] = axis_j(s, ax, tvec[ax], r);
            }
            double kr, ki;
            if (!leaf_to_cdouble(s->ker_leaves[kf], &kr, &ki)) goto fail;
            bool dropped;
            Expr* v = listval_leaf(s, jvec, &dropped);
            if (dropped) { ar += kr; ai += ki; continue; }   /* Times[ker] = ker */
            double vr, vi;
            if (!leaf_to_cdouble(v, &vr, &vi)) goto fail;
            ar += kr * vr - ki * vi;
            ai += kr * vi + ki * vr;
        }
        leaves[o] = num_leaf(ar, ai, any_complex);
        continue;
    fail:
        for (int64_t q = 0; q < o; q++) expr_free(leaves[q]);
        free(leaves);
        return NULL;
    }
    size_t idx = 0;
    Expr* result = build_nested(leaves, s->Ldims, s->rank, 0, &idx);
    free(leaves);
    return result;
}

#ifdef USE_MPFR
/* Arbitrary-precision analogue of conv_direct_machine (ncpx accumulation). */
static Expr* conv_direct_mpfr(const ConvSpec* s, bool any_complex, long target_bits) {
    int64_t Ltot = 1, Lstr[CONV_MAX_RANK];
    for (int ax = 0; ax < s->rank; ax++) Ltot *= s->Ldims[ax];
    row_major_strides(s->Ldims, s->rank, Lstr);
    int64_t mtot = 1;
    for (int ax = 0; ax < s->rank; ax++) mtot *= s->kdims[ax];

    long guard = 16;
    { int64_t t = mtot; while (t > 1) { guard += 1; t >>= 1; } }
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    ncpx acc, kv, vv, prod;
    ncpx_init(&acc, wp); ncpx_init(&kv, wp); ncpx_init(&vv, wp); ncpx_init(&prod, wp);
    mpfr_t re, im; mpfr_init2(re, wp); mpfr_init2(im, wp);

    Expr** leaves = malloc((size_t)Ltot * sizeof(Expr*));
    int ok = 1;
    int64_t o = 0;
    for (; o < Ltot && ok; o++) {
        int64_t tvec[CONV_MAX_RANK];
        for (int ax = 0; ax < s->rank; ax++) tvec[ax] = (o / Lstr[ax]) % s->Ldims[ax] + 1;
        ncpx_set_ui(&acc, 0);
        for (int64_t kf = 0; kf < mtot; kf++) {
            int64_t jvec[CONV_MAX_RANK];
            for (int ax = 0; ax < s->rank; ax++) {
                int64_t r = (kf / s->kstride[ax]) % s->kdims[ax] + 1;
                jvec[ax] = axis_j(s, ax, tvec[ax], r);
            }
            bool inexact;
            if (!get_approx_mpfr(s->ker_leaves[kf], re, im, &inexact)) { ok = 0; break; }
            mpfr_set(kv.re, re, MPFR_RNDN); mpfr_set(kv.im, im, MPFR_RNDN);
            bool dropped;
            Expr* v = listval_leaf(s, jvec, &dropped);
            if (dropped) { ncpx_add(&acc, &acc, &kv); continue; }
            if (!get_approx_mpfr(v, re, im, &inexact)) { ok = 0; break; }
            mpfr_set(vv.re, re, MPFR_RNDN); mpfr_set(vv.im, im, MPFR_RNDN);
            ncpx_mul(&prod, &kv, &vv, wp);
            ncpx_add(&acc, &acc, &prod);
        }
        if (!ok) break;
        mpfr_t rre, rim; mpfr_init2(rre, target_bits); mpfr_init2(rim, target_bits);
        mpfr_set(rre, acc.re, MPFR_RNDN);
        if (!any_complex) {
            mpfr_set_zero(rim, 1);
        } else {
            mpfr_set(rim, acc.im, MPFR_RNDN);
            if (!mpfr_zero_p(rim)) {
                mpfr_t mag, thr; mpfr_init2(mag, target_bits); mpfr_init2(thr, target_bits);
                mpfr_hypot(mag, rre, rim, MPFR_RNDN);
                mpfr_mul_2si(thr, mag, -(long)target_bits, MPFR_RNDN);
                if (mpfr_cmpabs(rim, thr) <= 0) mpfr_set_zero(rim, 1);
                mpfr_clear(mag); mpfr_clear(thr);
            }
        }
        leaves[o] = numeric_mpfr_make_complex(rre, rim);
        mpfr_clear(rre); mpfr_clear(rim);
    }

    ncpx_clear(&acc); ncpx_clear(&kv); ncpx_clear(&vv); ncpx_clear(&prod);
    mpfr_clear(re); mpfr_clear(im);

    if (!ok) {
        for (int64_t q = 0; q < o; q++) expr_free(leaves[q]);
        free(leaves);
        return NULL;
    }
    size_t idx = 0;
    Expr* result = build_nested(leaves, s->Ldims, s->rank, 0, &idx);
    free(leaves);
    return result;
}
#endif /* USE_MPFR */

/* ==================================================================== *
 *  FFT fast path (default Times/Plus, inexact numeric data)
 * ==================================================================== */

/* Per-axis materialisation window: A[p] = listval(jmin + p), p = 0..La-1. */
static void fft_window(const ConvSpec* s, int64_t* jmin, int64_t* La,
                       int64_t* Pdims) {
    for (int ax = 0; ax < s->rank; ax++) {
        int64_t m = s->kdims[ax], L = s->Ldims[ax], KL = s->KL[ax];
        jmin[ax] = (s->mode == CONV_MODE_CORRELATE) ? (2 - KL) : (1 + KL - m);
        La[ax]   = L + m - 1;
        Pdims[ax] = La[ax] + m - 1;      /* linear-convolution length */
    }
}

/* Machine-precision FFT convolution. Returns a nested List, or NULL to fall back. */
static Expr* conv_fft_machine(const ConvSpec* s, bool any_complex) {
    int64_t jmin[CONV_MAX_RANK], La[CONV_MAX_RANK], Pdims[CONV_MAX_RANK];
    fft_window(s, jmin, La, Pdims);
    int64_t P = 1;
    for (int ax = 0; ax < s->rank; ax++) P *= Pdims[ax];
    if (P <= 0) return NULL;

    int64_t Pstr[CONV_MAX_RANK];
    row_major_strides(Pdims, s->rank, Pstr);

    double* Kbuf = calloc((size_t)P * 2, sizeof(double));
    double* Abuf = calloc((size_t)P * 2, sizeof(double));
    if (!Kbuf || !Abuf) { free(Kbuf); free(Abuf); return NULL; }

    /* Place the kernel (reversed on every axis for correlate). */
    int64_t mtot = 1;
    for (int ax = 0; ax < s->rank; ax++) mtot *= s->kdims[ax];
    for (int64_t o = 0; o < mtot; o++) {
        int64_t pos = 0, rem = o;
        for (int ax = 0; ax < s->rank; ax++) {
            int64_t iax = (rem / s->kstride[ax]) % s->kdims[ax];
            int64_t p = (s->mode == CONV_MODE_CORRELATE) ? (s->kdims[ax] - 1 - iax) : iax;
            pos += p * Pstr[ax];
        }
        double re, im;
        if (!leaf_to_cdouble(s->ker_leaves[o], &re, &im)) { free(Kbuf); free(Abuf); return NULL; }
        Kbuf[2 * pos] = re; Kbuf[2 * pos + 1] = im;
    }

    /* Materialise the padded list window into Abuf. */
    int64_t Atot = 1, Astr[CONV_MAX_RANK];
    for (int ax = 0; ax < s->rank; ax++) Atot *= La[ax];
    row_major_strides(La, s->rank, Astr);
    for (int64_t o = 0; o < Atot; o++) {
        int64_t jvec[CONV_MAX_RANK], pos = 0;
        for (int ax = 0; ax < s->rank; ax++) {
            int64_t p = (o / Astr[ax]) % La[ax];
            jvec[ax] = jmin[ax] + p;
            pos += p * Pstr[ax];
        }
        bool dropped;
        Expr* v = listval_leaf(s, jvec, &dropped);
        double re, im;
        if (!v || !leaf_to_cdouble(v, &re, &im)) { free(Kbuf); free(Abuf); return NULL; }
        Abuf[2 * pos] = re; Abuf[2 * pos + 1] = im;
    }

    /* Circular convolution via the FFT product; zero padding => linear conv. */
    fourier_fft_machine(Kbuf, Pdims, s->rank, +1);
    fourier_fft_machine(Abuf, Pdims, s->rank, +1);
    for (int64_t o = 0; o < P; o++) {
        double kr = Kbuf[2 * o], ki = Kbuf[2 * o + 1];
        double ar = Abuf[2 * o], ai = Abuf[2 * o + 1];
        Kbuf[2 * o]     = kr * ar - ki * ai;
        Kbuf[2 * o + 1] = kr * ai + ki * ar;
    }
    fourier_fft_machine(Kbuf, Pdims, s->rank, -1);
    double invP = 1.0 / (double)P;

    /* Extract the L outputs at index (o_ax + m_ax - 1) per axis. */
    int64_t Ltot = 1, Lstr[CONV_MAX_RANK];
    for (int ax = 0; ax < s->rank; ax++) Ltot *= s->Ldims[ax];
    row_major_strides(s->Ldims, s->rank, Lstr);
    Expr** leaves = malloc((size_t)Ltot * sizeof(Expr*));
    for (int64_t o = 0; o < Ltot; o++) {
        int64_t pos = 0;
        for (int ax = 0; ax < s->rank; ax++) {
            int64_t oc = (o / Lstr[ax]) % s->Ldims[ax];
            pos += (oc + s->kdims[ax] - 1) * Pstr[ax];
        }
        double re = Kbuf[2 * pos] * invP, im = Kbuf[2 * pos + 1] * invP;
        if (!any_complex) {
            leaves[o] = expr_new_real(re);
        } else if (fabs(im) <= 16.0 * DBL_EPSILON * (fabs(re) + fabs(im))) {
            leaves[o] = expr_new_real(re);
        } else {
            leaves[o] = make_complex(expr_new_real(re), expr_new_real(im));
        }
    }
    free(Kbuf); free(Abuf);

    size_t idx = 0;
    Expr* result = build_nested(leaves, s->Ldims, s->rank, 0, &idx);
    free(leaves);
    return result;
}

static Expr* build_nested(Expr** leaves, const int64_t* dims, int rank,
                          int level, size_t* idx) {
    if (level == rank) return leaves[(*idx)++];
    int64_t n = dims[level];
    Expr** a = malloc((size_t)n * sizeof(Expr*));
    for (int64_t i = 0; i < n; i++)
        a[i] = build_nested(leaves, dims, rank, level + 1, idx);
    return mk_head_fn(SYM_List, NULL, a, (size_t)n);
}

#ifdef USE_MPFR
/* Per-axis strided in-place MPFR FFT over a row-major block. */
static void mpfr_fft_block(ncpx* buf, const int64_t* dims, int rank, int sign,
                           mpfr_prec_t wp) {
    int64_t N = 1;
    for (int i = 0; i < rank; i++) N *= dims[i];
    int64_t strides[CONV_MAX_RANK];
    row_major_strides(dims, rank, strides);
    for (int ax = 0; ax < rank; ax++) {
        int64_t n = dims[ax], sstep = strides[ax];
        if (n <= 1) continue;
        ncpx* line = malloc((size_t)n * sizeof(ncpx));
        for (int64_t k = 0; k < n; k++) ncpx_init(&line[k], wp);
        for (int64_t o = 0; o < N; o++) {
            if ((o / sstep) % n != 0) continue;
            for (int64_t r = 0; r < n; r++) ncpx_set(&line[r], &buf[o + r * sstep]);
            fourier_fft_mpfr(line, (size_t)n, sign, wp);
            for (int64_t r = 0; r < n; r++) ncpx_set(&buf[o + r * sstep], &line[r]);
        }
        for (int64_t k = 0; k < n; k++) ncpx_clear(&line[k]);
        free(line);
    }
}

static Expr* conv_fft_mpfr(const ConvSpec* s, bool any_complex, long target_bits) {
    int64_t jmin[CONV_MAX_RANK], La[CONV_MAX_RANK], Pdims[CONV_MAX_RANK];
    fft_window(s, jmin, La, Pdims);
    int64_t P = 1;
    for (int ax = 0; ax < s->rank; ax++) P *= Pdims[ax];
    if (P <= 0) return NULL;

    long guard = 32;
    { int64_t t = P; while (t > 1) { guard += 2; t >>= 1; } }
    mpfr_prec_t wp = (mpfr_prec_t)(target_bits + guard);

    int64_t Pstr[CONV_MAX_RANK];
    row_major_strides(Pdims, s->rank, Pstr);

    ncpx* Kbuf = malloc((size_t)P * sizeof(ncpx));
    ncpx* Abuf = malloc((size_t)P * sizeof(ncpx));
    for (int64_t o = 0; o < P; o++) {
        ncpx_init(&Kbuf[o], wp); ncpx_set_ui(&Kbuf[o], 0);
        ncpx_init(&Abuf[o], wp); ncpx_set_ui(&Abuf[o], 0);
    }

    int ok = 1;
    mpfr_t re, im; mpfr_init2(re, wp); mpfr_init2(im, wp);

    int64_t mtot = 1;
    for (int ax = 0; ax < s->rank; ax++) mtot *= s->kdims[ax];
    for (int64_t o = 0; o < mtot && ok; o++) {
        int64_t pos = 0;
        for (int ax = 0; ax < s->rank; ax++) {
            int64_t iax = (o / s->kstride[ax]) % s->kdims[ax];
            int64_t p = (s->mode == CONV_MODE_CORRELATE) ? (s->kdims[ax] - 1 - iax) : iax;
            pos += p * Pstr[ax];
        }
        bool inexact;
        if (!get_approx_mpfr(s->ker_leaves[o], re, im, &inexact)) { ok = 0; break; }
        mpfr_set(Kbuf[pos].re, re, MPFR_RNDN);
        mpfr_set(Kbuf[pos].im, im, MPFR_RNDN);
    }

    int64_t Atot = 1, Astr[CONV_MAX_RANK];
    for (int ax = 0; ax < s->rank; ax++) Atot *= La[ax];
    row_major_strides(La, s->rank, Astr);
    for (int64_t o = 0; o < Atot && ok; o++) {
        int64_t jvec[CONV_MAX_RANK], pos = 0;
        for (int ax = 0; ax < s->rank; ax++) {
            int64_t p = (o / Astr[ax]) % La[ax];
            jvec[ax] = jmin[ax] + p;
            pos += p * Pstr[ax];
        }
        bool dropped, inexact;
        Expr* v = listval_leaf(s, jvec, &dropped);
        if (!v || !get_approx_mpfr(v, re, im, &inexact)) { ok = 0; break; }
        mpfr_set(Abuf[pos].re, re, MPFR_RNDN);
        mpfr_set(Abuf[pos].im, im, MPFR_RNDN);
    }
    mpfr_clear(re); mpfr_clear(im);

    Expr* result = NULL;
    if (ok) {
        mpfr_fft_block(Kbuf, Pdims, s->rank, +1, wp);
        mpfr_fft_block(Abuf, Pdims, s->rank, +1, wp);
        ncpx tmp; ncpx_init(&tmp, wp);
        for (int64_t o = 0; o < P; o++) {
            ncpx_mul(&tmp, &Kbuf[o], &Abuf[o], wp);
            ncpx_set(&Kbuf[o], &tmp);
        }
        ncpx_clear(&tmp);
        mpfr_fft_block(Kbuf, Pdims, s->rank, -1, wp);

        mpfr_t invP; mpfr_init2(invP, wp);
        mpfr_set_si(invP, (long)P, MPFR_RNDN);
        mpfr_ui_div(invP, 1, invP, MPFR_RNDN);

        int64_t Ltot = 1, Lstr[CONV_MAX_RANK];
        for (int ax = 0; ax < s->rank; ax++) Ltot *= s->Ldims[ax];
        row_major_strides(s->Ldims, s->rank, Lstr);
        Expr** leaves = malloc((size_t)Ltot * sizeof(Expr*));
        for (int64_t o = 0; o < Ltot; o++) {
            int64_t pos = 0;
            for (int ax = 0; ax < s->rank; ax++) {
                int64_t oc = (o / Lstr[ax]) % s->Ldims[ax];
                pos += (oc + s->kdims[ax] - 1) * Pstr[ax];
            }
            mpfr_t rre, rim; mpfr_init2(rre, target_bits); mpfr_init2(rim, target_bits);
            mpfr_mul(rre, Kbuf[pos].re, invP, MPFR_RNDN);
            mpfr_mul(rim, Kbuf[pos].im, invP, MPFR_RNDN);
            if (!any_complex) {
                mpfr_set_zero(rim, 1);
            } else if (!mpfr_zero_p(rim)) {
                mpfr_t mag, thr; mpfr_init2(mag, target_bits); mpfr_init2(thr, target_bits);
                mpfr_hypot(mag, rre, rim, MPFR_RNDN);
                mpfr_mul_2si(thr, mag, -(long)target_bits, MPFR_RNDN);
                if (mpfr_cmpabs(rim, thr) <= 0) mpfr_set_zero(rim, 1);
                mpfr_clear(mag); mpfr_clear(thr);
            }
            leaves[o] = numeric_mpfr_make_complex(rre, rim);
            mpfr_clear(rre); mpfr_clear(rim);
        }
        mpfr_clear(invP);
        size_t idx = 0;
        result = build_nested(leaves, s->Ldims, s->rank, 0, &idx);
        free(leaves);
    }

    for (int64_t o = 0; o < P; o++) { ncpx_clear(&Kbuf[o]); ncpx_clear(&Abuf[o]); }
    free(Kbuf); free(Abuf);
    return result;
}
#endif /* USE_MPFR */

/* ==================================================================== *
 *  Classification + dispatch
 * ==================================================================== */

typedef enum { CLS_SYMBOLIC, CLS_EXACT, CLS_MACHINE, CLS_ARB } NumClass;

/* An inexact numeric leaf in the Mathematica sense (Real / MPFR, incl. inside
 * Complex[...]). Build-safe: MPFR is only referenced under USE_MPFR. */
static bool leaf_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *cr, *ci;
    if (is_complex((Expr*)e, &cr, &ci))
        return leaf_is_inexact(cr) || leaf_is_inexact(ci);
    return false;
}

static NumClass classify_leaves(Expr** a, int64_t na, Expr** b, int64_t nb,
                                bool* any_complex, long* arb_bits) {
    bool all_num = true, any_inexact = false;
    *any_complex = false;
    (void)arb_bits;
    Expr** arrs[2] = { a, b };
    int64_t counts[2] = { na, nb };
#ifdef USE_MPFR
    bool any_mpfr = false;
    long minb = 0;
#endif
    for (int k = 0; k < 2; k++) {
        for (int64_t i = 0; i < counts[k]; i++) {
            Expr* e = arrs[k][i];
            if (!expr_is_numeric_like(e)) all_num = false;
            Expr *cr, *ci;
            if (is_complex(e, &cr, &ci)) *any_complex = true;
            if (leaf_is_inexact(e)) any_inexact = true;
#ifdef USE_MPFR
            if (numeric_expr_is_mpfr(e)) {
                any_mpfr = true;
                long ib = numeric_min_inexact_bits(e);
                if (ib > 0 && (minb == 0 || ib < minb)) minb = ib;
            }
#endif
        }
    }
    if (!all_num) return CLS_SYMBOLIC;
#ifdef USE_MPFR
    if (any_mpfr) { if (minb < 53) minb = 53; *arb_bits = minb; return CLS_ARB; }
#endif
    if (any_inexact) return CLS_MACHINE;
    return CLS_EXACT;
}

/* ==================================================================== *
 *  Argument parsing
 * ==================================================================== */

static bool get_int(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        *out = mpz_get_si(e->data.bigint); return true;
    }
    return false;
}

/* Normalise kL,kR to a positive kernel index in 1..m. Returns 0 on invalid. */
static int64_t norm_k(int64_t k, int64_t m) {
    if (k > 0) return k;
    if (k < 0) return m + 1 + k;
    return 0;
}

/* Parse the klist argument into per-axis kL/kR (raw, pre-normalisation). */
static bool parse_klist(const Expr* kl, ConvMode mode, int rank,
                        int64_t* kLraw, int64_t* kRraw) {
    if (kl == NULL) {
        int64_t dL = (mode == CONV_MODE_CONVOLVE) ? -1 : 1;
        int64_t dR = (mode == CONV_MODE_CONVOLVE) ?  1 : -1;
        for (int ax = 0; ax < rank; ax++) { kLraw[ax] = dL; kRraw[ax] = dR; }
        return true;
    }
    int64_t k;
    if (get_int(kl, &k)) {
        for (int ax = 0; ax < rank; ax++) { kLraw[ax] = k; kRraw[ax] = k; }
        return true;
    }
    if (!head_is(kl, SYM_List) || kl->data.function.arg_count != 2) return false;
    Expr* e0 = kl->data.function.args[0];
    Expr* e1 = kl->data.function.args[1];
    int64_t a, b;
    if (get_int(e0, &a) && get_int(e1, &b)) {           /* {kL, kR} broadcast */
        for (int ax = 0; ax < rank; ax++) { kLraw[ax] = a; kRraw[ax] = b; }
        return true;
    }
    /* {{kL1,..},{kR1,..}} per axis */
    if (!head_is(e0, SYM_List) || !head_is(e1, SYM_List)) return false;
    if ((int)e0->data.function.arg_count != rank ||
        (int)e1->data.function.arg_count != rank) return false;
    for (int ax = 0; ax < rank; ax++) {
        if (!get_int(e0->data.function.args[ax], &kLraw[ax])) return false;
        if (!get_int(e1->data.function.args[ax], &kRraw[ax])) return false;
    }
    return true;
}

/* ==================================================================== *
 *  Engine entry point
 * ==================================================================== */

Expr* conv_engine(Expr* res, ConvMode mode) {
    const char* name = (mode == CONV_MODE_CONVOLVE) ? "ListConvolve" : "ListCorrelate";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 7) return builtin_arg_error(name, argc, 2, 7);

    Expr* ker  = res->data.function.args[0];
    Expr* list = res->data.function.args[1];
    Expr* klist_arg = (argc >= 3) ? res->data.function.args[2] : NULL;
    Expr* pad_arg   = (argc >= 4) ? res->data.function.args[3] : NULL;
    Expr* g_arg     = (argc >= 5) ? res->data.function.args[4] : NULL;
    Expr* h_arg     = (argc >= 6) ? res->data.function.args[5] : NULL;
    Expr* lev_arg   = (argc >= 7) ? res->data.function.args[6] : NULL;

    ConvSpec s;
    memset(&s, 0, sizeof(s));
    s.mode = mode;
    s.g_head = g_arg;
    s.h_head = h_arg;

    /* Shapes. */
    int krank = 0, lrank = 0;
    if (!nest_dims(ker, s.kdims, &krank)) return NULL;
    if (!nest_dims(list, s.ldims, &lrank)) return NULL;
    if (krank != lrank) return NULL;          /* differing ranks: unsupported */
    s.rank = krank;
    for (int ax = 0; ax < s.rank; ax++)
        if (s.kdims[ax] <= 0 || s.ldims[ax] <= 0) return NULL;

    /* Optional level argument: only lev == rank(ker) is supported. */
    if (lev_arg) {
        int64_t lev;
        if (!get_int(lev_arg, &lev) || lev != s.rank) return NULL;
    }

    /* Overhang parameters -> normalised positive kernel indices + output dims. */
    int64_t kLraw[CONV_MAX_RANK], kRraw[CONV_MAX_RANK];
    if (!parse_klist(klist_arg, mode, s.rank, kLraw, kRraw)) return NULL;
    for (int ax = 0; ax < s.rank; ax++) {
        s.KL[ax] = norm_k(kLraw[ax], s.kdims[ax]);
        s.KR[ax] = norm_k(kRraw[ax], s.kdims[ax]);
        if (s.KL[ax] == 0 || s.KR[ax] == 0) return NULL;
        s.Ldims[ax] = (mode == CONV_MODE_CORRELATE)
                        ? s.ldims[ax] + s.KL[ax] - s.KR[ax]
                        : s.ldims[ax] - s.KL[ax] + s.KR[ax];
        if (s.Ldims[ax] < 1) return NULL;     /* empty / invalid window */
    }
    row_major_strides(s.kdims, s.rank, s.kstride);
    row_major_strides(s.ldims, s.rank, s.lstride);

    /* Kernel + list leaves. */
    int64_t ktot = 1, ltot = 1;
    for (int ax = 0; ax < s.rank; ax++) { ktot *= s.kdims[ax]; ltot *= s.ldims[ax]; }
    s.ker_leaves  = malloc((size_t)ktot * sizeof(Expr*));
    s.list_leaves = malloc((size_t)ltot * sizeof(Expr*));
    size_t kidx = 0, lidx = 0;
    nest_flatten(ker, 0, s.rank, s.ker_leaves, &kidx);
    nest_flatten(list, 0, s.rank, s.list_leaves, &lidx);

    /* Padding. */
    Expr** pad_leaves_owned = NULL;
    if (pad_arg == NULL) {
        s.pad_kind = PAD_LIST;
        s.cyclic_default = true;
        memcpy(s.pdims, s.ldims, sizeof(s.pdims));
        memcpy(s.pstride, s.lstride, sizeof(s.pstride));
        s.pad_leaves = s.list_leaves;
    } else if (head_is(pad_arg, SYM_List)) {
        if (pad_arg->data.function.arg_count == 0) {
            s.pad_kind = PAD_EMPTY;
        } else {
            int prank = 0;
            if (!nest_dims(pad_arg, s.pdims, &prank) || prank != s.rank) {
                free(s.ker_leaves); free(s.list_leaves); return NULL;
            }
            s.pad_kind = PAD_LIST;
            int64_t ptot = 1;
            for (int ax = 0; ax < s.rank; ax++) {
                if (s.pdims[ax] <= 0) { free(s.ker_leaves); free(s.list_leaves); return NULL; }
                ptot *= s.pdims[ax];
            }
            row_major_strides(s.pdims, s.rank, s.pstride);
            pad_leaves_owned = malloc((size_t)ptot * sizeof(Expr*));
            size_t pidx = 0;
            nest_flatten(pad_arg, 0, s.rank, pad_leaves_owned, &pidx);
            s.pad_leaves = pad_leaves_owned;
        }
    } else {
        s.pad_kind = PAD_CONST;
        s.pad_scalar = pad_arg;
    }

    /* Decide numeric regime. FFT only for default Times/Plus, non-empty pad. */
    Expr* result = NULL;
    bool default_gh = (g_arg == NULL && h_arg == NULL);
    int64_t work = 1;
    for (int ax = 0; ax < s.rank; ax++) work *= s.Ldims[ax] * s.kdims[ax];

    if (default_gh) {
        bool any_complex = false;
        long arb_bits = 0;
        NumClass cls = classify_leaves(s.ker_leaves, ktot, s.list_leaves, ltot,
                                       &any_complex, &arb_bits);
        bool big = (s.pad_kind != PAD_EMPTY && work >= CONV_FFT_MIN_WORK);
        if (cls == CLS_MACHINE) {
            /* FFT above the threshold, tight O(L*m) numeric loop below it;
             * either may return NULL (unrepresentable leaf) -> fall through. */
            if (big) result = conv_fft_machine(&s, any_complex);
            if (!result) result = conv_direct_machine(&s, any_complex);
        } else if (cls == CLS_ARB) {
#ifdef USE_MPFR
            if (big) result = conv_fft_mpfr(&s, any_complex, arb_bits);
            if (!result) result = conv_direct_mpfr(&s, any_complex, arb_bits);
#endif
        }
    }

    /* Symbolic, exact, custom g/h, or a numeric leaf the fast paths could not
     * represent: the general Expr-building engine (exact stays exact). */
    if (!result) result = conv_direct(&s);

    free(s.ker_leaves);
    free(s.list_leaves);
    free(pad_leaves_owned);
    return result;
}

Expr* builtin_list_convolve(Expr* res) { return conv_engine(res, CONV_MODE_CONVOLVE); }

void convolutions_init(void) {
    symtab_add_builtin("ListConvolve", builtin_list_convolve);
    symtab_get_def("ListConvolve")->attributes |= ATTR_PROTECTED;
}
