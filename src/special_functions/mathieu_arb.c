/* Mathilda -- arbitrary-precision Mathieu functions (arb / acb kernel).
 *
 * See mathieu_arb.h and sf_machine_mathieu.c (the machine twin, which documents
 * the algorithm and the normalization). Outline at working precision wp:
 *
 *   0. The odd functions (MathieuS, odd == 1) share every step: the same nu
 *      and c_k, a sine sum in place of the cosine sum (C + i S is the Floquet
 *      solution sum c_k e^{i (nu+2k) z}).
 *   1. q == 0 exactly: cos(sqrt(a) z) / -sqrt(a) sin(sqrt(a) z) (sin / sqrt(a)
 *      cos for S); a == 0 -> 1/sqrt2 (0 for S).
 *   2. The machine plan (sf_machine_mathieu_plan_info) gives nu0, the dominant
 *      index m, and whether a is a characteristic value (nu snapped to an
 *      integer n; then nu = n exactly and the DLMF factor 1/sqrt 2 applies,
 *      and the function of the other parity -- S at a = a_n, C at a = b_n -- is
 *      exactly 0).
 *   3. Otherwise nu is refined by Newton on
 *          g(nu) = (a - (nu+2m)^2) - q (c_{m+1}/c_m + c_{m-1}/c_m),
 *      the ratios from continued fractions. Precision doubles per level from 64
 *      bits with ONE Newton step per level (quadratic convergence: a nu good to
 *      p/2 bits becomes good to p), so the whole refinement costs about two
 *      continued-fraction passes at wp. Arithmetic is on ball midpoints during
 *      the iteration; the converged nu gets an explicit error radius.
 *   4. Coefficients: the same continued fractions (balls), glued at m;
 *      normalized by sum c_k^2 = 1 with Re(c_0) > 0 (the sign decided on
 *      midpoints, exactly as the machine plan decides it).
 *   5. Sum: real a, q, z, nu -> real Clenshaw on cos 2kz / sin 2kz; otherwise the
 *      centred complex Horner of the machine helper.
 *
 * PRECISION PROFILE (what makes 10^4 digits cheap). |c_k| falls off super-
 * exponentially -- about (q/4)^k / (k!)^2 -- so a tail coefficient only
 * contributes its leading wp - b_k bits, b_k = log2(max term / term k). Every
 * operation at index k (the continued-fraction step, the product, the square in
 * the norm, the Clenshaw / Horner step) runs at p_k = wp + 32 - b_k bits (floor
 * 64). The continued fraction is self-correcting toward m (an error in r_{k+1}
 * enters r_k damped by |q r_{k+1}/d_k| ~ |c_{k+1}/c_{k-1}|), so the reduced tail
 * precision never leaks into the dominant coefficients. The profile b_k comes
 * from the asymptotic ratio |q| / |(nu+2k)^2 - a| capped at 1 (an upper bound on
 * |c_k|, so p_k is never too small), weighted by e^{2|k||Im z|} for complex z.
 */
/* pthread_sigmask / sigfillset (the threaded halves below) are POSIX: glibc
 * hides them under -std=c99 unless this precedes every #include. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE 1
#endif

#ifdef USE_FLINT

#include "mathieu_arb.h"
#include "sf_machine_mathieu.h"
#include "core.h"            /* tc_async_defer_push / pop, tc_check_deadline */
#include "term/term_progress.h"   /* worker count for the REPL spinner */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>   /* Bessel J, Y: the far-from-real-axis expansion */

#ifdef MATHILDA_THREADS
#include <pthread.h>
#include <signal.h>
#endif

/* ------------------------------------------------------------------ */
/* Two-way parallelism                                                 */
/* ------------------------------------------------------------------ */

/* The kernel has three pairs of independent halves -- the two Miller sweeps
 * (from +infinity and from -infinity toward m), the two Clenshaw sums (cosine
 * and sine), the coefficient corrections on either side of m -- each a long run
 * of big-number multiplications. Above MATHIEU_PAR_BITS they run on two worker
 * threads under the discipline of src/numbertheory/factorinteger/fi_par.c:
 * workers start with every signal blocked, the region is bracketed by
 * tc_async_defer_push/pop, and the calling thread only waits in pthread_join, so
 * a TimeConstrained timeout is taken at its next safe point after both workers
 * finish (the guarded GMP/FLINT allocators only jump on the main thread). A
 * worker that cannot be created runs inline after the joins. */
#define MATHIEU_PAR_BITS 4096

#ifdef MATHILDA_THREADS
typedef struct { void (*fn)(void*); void* arg; } MathieuJob;
static void* mathieu_job_run(void* p) {
    MathieuJob* j = (MathieuJob*)p;
    j->fn(j->arg);
    flint_cleanup();                       /* this thread's FLINT caches */
    return NULL;
}
#endif

static void mathieu_par2(void (*fn)(void*), void* a0, void* a1, int par) {
#ifdef MATHILDA_THREADS
    if (par) {
        pthread_t th[2];
        MathieuJob job[2];
        int started[2];
        sigset_t all, old;
        sigfillset(&all);
        memset(&old, 0, sizeof old);       /* written by the first call */
        pthread_sigmask(SIG_SETMASK, &all, &old);
        tc_async_defer_push();
        job[0].fn = fn; job[0].arg = a0;
        job[1].fn = fn; job[1].arg = a1;
        mth_progress_workers_add(2);
        for (int i = 0; i < 2; i++)
            started[i] = pthread_create(&th[i], NULL, mathieu_job_run, &job[i]) == 0;
        pthread_sigmask(SIG_SETMASK, &old, NULL);
        for (int i = 0; i < 2; i++) if (started[i]) pthread_join(th[i], NULL);
        mth_progress_workers_add(-2);
        for (int i = 0; i < 2; i++) if (!started[i]) fn(job[i].arg);
        tc_async_defer_pop();
        tc_check_deadline();               /* a timeout due now unwinds here */
        return;
    }
#else
    (void)par;
#endif
    fn(a0);
    fn(a1);
}

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* ------------------------------------------------------------------ */
/* The window and its precision profile                                */
/* ------------------------------------------------------------------ */

typedef struct {
    slong lo, hi;       /* coefficient indices lo .. hi (lo < 0 < hi) */
    slong* p;           /* p[k - lo]: working precision at index k */
} Profile;

/* One side of the profile: log upper bounds t_j (natural log, weighted by
 * e^{2|k|y}) for k = m + sg*j, j = 1.., stopping `bits + 40` bits below the
 * running maximum. Returns the count; *out is malloc'd (caller frees). */
static slong side_profile(double nur, double nui, double ar, double ai,
                          double qr, double qi, slong m, int sg, double y,
                          slong bits, double* tmax, double** out) {
    double aq = hypot(qr, qi);
    double L = 0.0, best = *tmax;
    double drop = ((double)bits + 40.0) * M_LN2;
    slong cap = 64, j;
    double* t = (double*)malloc(sizeof(double) * (size_t)cap);
    if (!t) { *out = NULL; return 0; }
    for (j = 1; j < 4000000; j++) {
        double k = (double)(m + sg * j);
        double tr = nur + 2.0 * k, ti = nui;
        double dr = tr * tr - ti * ti - ar, di = 2.0 * tr * ti - ai;  /* (nu+2k)^2 - a */
        double den = hypot(dr, di);
        double ratio = den > aq ? aq / den : 1.0;
        L += log(ratio);
        double tj = L + 2.0 * fabs(k) * y;
        if (j >= cap) {
            cap *= 2;
            double* nt = (double*)realloc(t, sizeof(double) * (size_t)cap);
            if (!nt) { free(t); *out = NULL; return 0; }
            t = nt;
        }
        t[j] = tj;
        if (tj > best) best = tj;
        if (tj < best - drop && ratio < 0.5) break;
    }
    /* two guard terms beyond the cut */
    slong n = j + 2;
    if (n >= cap) {
        double* nt = (double*)realloc(t, sizeof(double) * (size_t)(n + 1));
        if (!nt) { free(t); *out = NULL; return 0; }
        t = nt;
    }
    for (slong i = j + 1; i <= n; i++) t[i] = t[j] - 1.0;
    *tmax = best;
    *out = t;
    return n;
}

/* Build the window [lo, hi] around m (always containing k = -1 .. 1, the sums
 * are centred on 0) and the per-index precisions for `bits` bits. */
static int profile_build(Profile* P, double nur, double nui, double ar, double ai,
                         double qr, double qi, slong m, double y, slong bits) {
    double tmax = 2.0 * fabs((double)m) * y;
    double *tp = NULL, *tm = NULL;
    slong np = side_profile(nur, nui, ar, ai, qr, qi, m, +1, y, bits, &tmax, &tp);
    slong nm = side_profile(nur, nui, ar, ai, qr, qi, m, -1, y, bits, &tmax, &tm);
    if (!tp || !tm) { free(tp); free(tm); return 0; }
    P->lo = m - nm; P->hi = m + np;
    if (P->lo > -1) P->lo = -1;
    if (P->hi < 1) P->hi = 1;
    slong n = P->hi - P->lo + 1;
    P->p = (slong*)malloc(sizeof(slong) * (size_t)n);
    if (!P->p) { free(tp); free(tm); return 0; }
    for (slong k = P->lo; k <= P->hi; k++) {
        double tk;
        slong j = k - m;
        if (j == 0) tk = 2.0 * fabs((double)m) * y;
        else if (j > 0) tk = j <= np ? tp[j] : tp[np];
        else tk = -j <= nm ? tm[-j] : tm[nm];
        double b = (tmax - tk) / M_LN2;                  /* bits below the max */
        slong pk = (slong)((double)bits + 32.0 - (b > 0.0 ? b : 0.0));
        P->p[k - P->lo] = pk < 64 ? 64 : pk;
    }
    free(tp); free(tm);
    return 1;
}

static void profile_clear(Profile* P) { free(P->p); P->p = NULL; }

#define PREC_AT(P, k) ((P)->p[(k) - (P)->lo])

/* ------------------------------------------------------------------ */
/* Miller sweeps                                                       */
/* ------------------------------------------------------------------ */

/* The recurrence  q c_{k-1} = d_k c_k - q c_{k+1},  d_k = a - (nu+2k)^2, run
 * backward from both ends of the window toward m (Miller's algorithm for the
 * minimal solution): no division per index, unlike the continued fraction of
 * ratios. The rounding at a tail index is harmless -- a perturbation of the
 * pair (e_k, e_{k+1}) splits into the minimal solution (a rescale) and the other
 * one, which DECAYS in the backward direction -- so each step runs at its
 * profile precision.
 *
 *   d_k = (a - nu^2) - 4k nu - 4k^2      (small-integer multiplies only),
 *   d'_k = -2 nu - 4k                    (d/dnu),
 *
 * and, when `de` is given, the nu-derivative is carried in forward mode at
 * precision <= pder (half the value precision suffices: see refine / correct).
 *
 * The sweeps run on MIDPOINTS (radii dropped each step). A ball would carry the
 * far end's 64-bit start as a 2^-64 relative radius all the way to m, although
 * that error is exactly a rescale of the minimal solution (cancelled by the
 * normalization) plus a component that decays backward. The caller attaches the
 * true bound, 2^-(p_k - 10) relative per coefficient, after normalizing. */
typedef struct {
    acb_t A0;          /* a - nu^2 */
    acb_t nu4;         /* 4 nu */
    acb_t mnu2;        /* -2 nu */
    acb_t qinv;        /* 1/q when q is not a short exact number */
    int qshort;        /* q exact with short mantissas: divide by q directly */
} MillerCtx;

static int acb_is_short_exact(const acb_t x) {
    return acb_is_exact(x) &&
           arf_bits(arb_midref(acb_realref(x))) <= 64 &&
           arf_bits(arb_midref(acb_imagref(x))) <= 64;
}

static void mctx_init(MillerCtx* M, const acb_t a, const acb_t q, const acb_t nu, slong prec) {
    acb_init(M->A0); acb_init(M->nu4); acb_init(M->mnu2); acb_init(M->qinv);
    acb_sqr(M->A0, nu, prec);
    acb_sub(M->A0, a, M->A0, prec);
    acb_mul_2exp_si(M->nu4, nu, 2);
    acb_mul_2exp_si(M->mnu2, nu, 1);
    acb_neg(M->mnu2, M->mnu2);
    M->qshort = acb_is_short_exact(q);
    if (!M->qshort) acb_inv(M->qinv, q, prec);
}
static void mctx_clear(MillerCtx* M) {
    acb_clear(M->A0); acb_clear(M->nu4); acb_clear(M->mnu2); acb_clear(M->qinv);
}
static inline void mdivq(acb_t x, const acb_t q, const MillerCtx* M, slong p) {
    if (M->qshort) acb_div(x, x, q, p);
    else acb_mul(x, x, M->qinv, p);
}
static inline void m_dk(acb_t d, const MillerCtx* M, slong k, slong p) {
    acb_mul_si(d, M->nu4, k, p);
    acb_sub(d, M->A0, d, p);
    acb_sub_si(d, d, 4 * k * k, p);
}

/* One side of the sweep pair (see miller_sweeps). */
typedef struct {
    acb_ptr e, de;               /* the shared vectors (each side writes its half) */
    acb_ptr em_minus, dem_minus; /* the - side's own e_m (and derivative) */
    const MillerCtx* M;
    acb_srcptr q;
    slong m;
    const Profile* P;
    slong pmax, pder;
    int side;                    /* 0: from hi down to m; 1: from lo up to m - 1 */
} SweepJob;

static void miller_side(void* arg) {
    const SweepJob* J = (const SweepJob*)arg;
    acb_ptr e = J->e, de = J->de;
    const MillerCtx* M = J->M;
    acb_srcptr q = J->q;
    slong m = J->m, pmax = J->pmax, pder = J->pder;
    const Profile* P = J->P;
    slong lo = P->lo, hi = P->hi;
    int side = J->side;
    acb_t d, dd, t, prev, dprev;
    acb_init(d); acb_init(dd); acb_init(t); acb_init(prev); acb_init(dprev);
    {
        int sg = side == 0 ? -1 : +1;              /* direction of travel */
        slong kfar = side == 0 ? hi : lo;
        /* e_far = 1, e_beyond = 0 */
        acb_one(e + (kfar - lo));
        if (de) acb_zero(de + (kfar - lo));
        acb_zero(prev); acb_zero(dprev);           /* e_{k - sg} (the one behind) */
        for (slong k = kfar; k != m; k += sg) {
            slong kn = k + sg;                     /* index being produced */
            slong pk = PREC_AT(P, kn) < pmax ? PREC_AT(P, kn) : pmax;
            acb_srcptr ek = e + (k - lo);
            m_dk(d, M, k, pk);
            acb_mul(t, d, ek, pk);
            if (!acb_is_zero(prev)) {
                acb_mul(dd, q, prev, pk);
                acb_sub(t, t, dd, pk);
            }
            mdivq(t, q, M, pk);
            acb_get_mid(t, t);
            if (de) {
                /* the derivative needs pder bits where the value needs pmax:
                 * keep the same offset down the tail (floor 64) */
                slong pd = pk - (pmax - pder);
                if (pd < 64) pd = 64;
                acb_srcptr dek = de + (k - lo);
                acb_t u; acb_init(u);
                acb_sub_si(dd, M->mnu2, 4 * k, pd);            /* d'_k */
                acb_mul(u, dd, ek, pd);
                acb_addmul(u, d, dek, pd);
                acb_submul(u, q, dprev, pd);
                mdivq(u, q, M, pd);
                acb_get_mid(u, u);
                acb_set(dprev, dek);
                if (kn == m && side == 1) acb_swap(J->dem_minus, u);
                else acb_swap(de + (kn - lo), u);
                acb_clear(u);
            }
            acb_set(prev, ek);
            if (kn == m && side == 1) acb_swap(J->em_minus, t);
            else acb_swap(e + (kn - lo), t);
        }
    }
    acb_clear(d); acb_clear(dd); acb_clear(t); acb_clear(prev); acb_clear(dprev);
}

/* Fill e[k - lo] (and de) for k in [lo, hi]: the + side sweep gives k >= m (its
 * e_m in e[m - lo]), the - side gives k < m and its own e_m in *em_minus. The
 * two sides carry independent scales, write disjoint entries, and run in
 * parallel at high precision. */
static void miller_sweeps(acb_ptr e, acb_ptr de, acb_t em_minus, acb_t dem_minus,
                          const MillerCtx* M, const acb_t q, slong m,
                          const Profile* P, slong pmax, slong pder) {
    SweepJob J[2];
    for (int s = 0; s < 2; s++) {
        J[s].e = e; J[s].de = de; J[s].em_minus = em_minus; J[s].dem_minus = dem_minus;
        J[s].M = M; J[s].q = q; J[s].m = m; J[s].P = P;
        J[s].pmax = pmax; J[s].pder = pder; J[s].side = s;
    }
    mathieu_par2(miller_side, &J[0], &J[1], pmax >= MATHIEU_PAR_BITS);
}

/* g(nu) = d_m - q (e+_{m+1}/e+_m + e-_{m-1}/e-_m) and (when de given) g'(nu). */
static void miller_g(acb_t g, acb_t dg, acb_srcptr e, acb_srcptr de,
                     const acb_t em_minus, const acb_t dem_minus,
                     const MillerCtx* M, const acb_t q, const acb_t nu, slong m,
                     slong lo, slong p, slong pder) {
    acb_t r1, r2, t, u;
    acb_init(r1); acb_init(r2); acb_init(t); acb_init(u);
    acb_srcptr ep = e + (m - lo), ep1 = e + (m + 1 - lo), em1 = e + (m - 1 - lo);
    acb_div(r1, ep1, ep, p);
    acb_div(r2, em1, em_minus, p);
    acb_add(t, r1, r2, p);
    acb_mul(t, t, q, p);
    m_dk(g, M, m, p);
    acb_sub(g, g, t, p);
    if (dg) {
        /* r' = (de_{m+-1} e_m - e_{m+-1} de_m) / e_m^2, at pder */
        acb_mul(t, de + (m + 1 - lo), ep, pder);
        acb_submul(t, ep1, de + (m - lo), pder);
        acb_div(t, t, ep, pder); acb_div(t, t, ep, pder);
        acb_mul(u, de + (m - 1 - lo), em_minus, pder);
        acb_submul(u, em1, dem_minus, pder);
        acb_div(u, u, em_minus, pder); acb_div(u, u, em_minus, pder);
        acb_add(t, t, u, pder);
        acb_mul(t, t, q, pder);
        acb_add_si(dg, nu, 2 * m, pder);
        acb_mul_2exp_si(dg, dg, 1);
        acb_neg(dg, dg);                                  /* d'_m */
        acb_sub(dg, dg, t, pder);
    }
    acb_clear(r1); acb_clear(r2); acb_clear(t); acb_clear(u);
}

/* ------------------------------------------------------------------ */
/* Summation                                                            */
/* ------------------------------------------------------------------ */

/* Real Clenshaw: A = sum_{k=0}^{N} al_k cos(k th), B = sum_{k=1}^{N} be_k sin(k th),
 * with c2 = 2 cos th; step k runs at pk[k] bits (the profile, folded onto |k|).
 *
 * The recurrence is run on MIDPOINTS and a rigorous bound is attached at the end.
 * Ball arithmetic cannot be used step by step here: the radius obeys
 * r_k = 2|cos th| r_{k+1} + r_{k+2}, which grows like ((|x| + sqrt(x^2+1)))^N --
 * about 2^950 at N = 1634 -- although the true error stays tiny (the Clenshaw
 * error for a cosine sum is bounded by sum_k |e_k| |U_{k-1}| <= sum (k+1) |e_k|).
 * The bound: (N+1)^2 (max rounding + max input radius) + rad(c2) N^3 max|coeff|. */
typedef struct {
    arb_ptr out;                 /* A (pass 0) or B (pass 1) */
    arb_srcptr co;               /* al or be */
    arb_srcptr al;               /* al_0, for the cosine pass */
    slong N;
    const slong* pk;
    arb_srcptr c2, cs, sn;
    slong prec;
    int pass;
} ClenJob;

static void clenshaw_pass(void* arg) {
    const ClenJob* J = (const ClenJob*)arg;
    arb_srcptr co = J->co;
    slong N = J->N;
    const slong* pk = J->pk;
    int pass = J->pass;
    arb_t b1, b2, t, x2;
    arb_init(b1); arb_init(b2); arb_init(t); arb_init(x2);
    arb_get_mid_arb(x2, J->c2);
    mag_t amax, rmax, rnd, err, tmp;
    mag_init(amax); mag_init(rmax); mag_init(rnd); mag_init(err); mag_init(tmp);
    arb_zero(b1); arb_zero(b2);
    for (slong k = N; k >= (pass == 0 ? 0 : 1); k--) {
        arb_get_mag(tmp, co + k);
        mag_max(amax, amax, tmp);
        mag_max(rmax, rmax, arb_radref(co + k));
        if (k == 0) break;
        arb_mul(t, x2, b1, pk[k]);
        arb_sub(t, t, b2, pk[k]);
        arb_add_arf(t, t, arb_midref(co + k), pk[k]);
        /* the rounding of this step, ~2^-pk |b_k|, folded into rnd */
        mag_max(rnd, rnd, arb_radref(t));
        mag_zero(arb_radref(t));
        arb_swap(b2, b1); arb_swap(b1, t);
    }
    if (pass == 0) {                     /* f = al_0 + b1 cos th - b2 */
        arb_mul(J->out, b1, J->cs, J->prec);
        arb_sub(J->out, J->out, b2, J->prec);
        arb_add(J->out, J->out, J->al + 0, J->prec);
    } else {                             /* g = b1 sin th */
        arb_mul(J->out, b1, J->sn, J->prec);
    }
    /* rigorous a-priori error: (N+1)^2 (rnd + rmax) + rad(c2) N^3 amax */
    mag_add(err, rnd, rmax);
    mag_mul_ui(err, err, (ulong)((N + 1) * (N + 1)));
    mag_mul_ui(tmp, arb_radref(J->c2), (ulong)(N + 1));
    mag_mul_ui(tmp, tmp, (ulong)((N + 1) * (N + 1)));
    mag_mul(tmp, tmp, amax);
    mag_add(err, err, tmp);
    arb_add_error_mag(J->out, err);
    mag_clear(amax); mag_clear(rmax); mag_clear(rnd); mag_clear(err); mag_clear(tmp);
    arb_clear(b1); arb_clear(b2); arb_clear(t); arb_clear(x2);
}

/* The two passes (cosine sum A, sine sum B) are independent: in parallel at
 * high precision. */
static void clenshaw_cs(arb_t A, arb_t B, arb_srcptr al, arb_srcptr be, slong N,
                        const slong* pk, const arb_t c2, const arb_t cs,
                        const arb_t sn, slong prec) {
    ClenJob J[2];
    for (int pass = 0; pass < 2; pass++) {
        J[pass].out = pass == 0 ? A : B;
        J[pass].co = pass == 0 ? al : be;
        J[pass].al = al;
        J[pass].N = N; J[pass].pk = pk;
        J[pass].c2 = c2; J[pass].cs = cs; J[pass].sn = sn;
        J[pass].prec = prec; J[pass].pass = pass;
    }
    mathieu_par2(clenshaw_pass, &J[0], &J[1], prec >= MATHIEU_PAR_BITS);
}

static void q0_closed_form(acb_t res, const acb_t a, const acb_t z, int odd, int deriv,
                           slong prec) {
    if (acb_is_zero(a)) {
        if (deriv || odd) acb_zero(res);
        else { acb_set_si(res, 2); acb_rsqrt(res, res, prec); }
        return;
    }
    acb_t s, u, sn, cs;
    acb_init(s); acb_init(u); acb_init(sn); acb_init(cs);
    acb_sqrt(s, a, prec);
    acb_mul(u, s, z, prec);
    acb_sin_cos(sn, cs, u, prec);
    if (odd) {
        if (deriv) acb_mul(res, s, cs, prec);
        else acb_set(res, sn);
        /* real a, z: sin / sqrt(a) cos are real for a > 0, i sinh / i cosh for a < 0 */
        if (acb_is_real(a) && acb_is_real(z)) {
            if (arb_is_negative(acb_realref(a))) arb_zero(acb_realref(res));
            else if (arb_is_positive(acb_realref(a))) arb_zero(acb_imagref(res));
        }
    } else {
        if (deriv) { acb_mul(res, s, sn, prec); acb_neg(res, res); }
        else acb_set(res, cs);
        if (acb_is_real(a) && acb_is_real(z)) arb_zero(acb_imagref(res));   /* cos / cosh */
    }
    acb_clear(s); acb_clear(u); acb_clear(sn); acb_clear(cs);
}

/* ------------------------------------------------------------------ */
/* The kernel                                                           */
/* ------------------------------------------------------------------ */

/* The final pass's per-coefficient work on c[i0 .. i1): the first-order Newton
 * correction c_k += delta c_k' with its second-order remainder (when delta is
 * given), then the glue onto the + side's scale (when ratio is given). */
typedef struct {
    acb_ptr c;
    acb_srcptr dc;
    const Profile* P;
    slong m, wp;
    acb_srcptr delta;
    const mag_struct* d2;        /* |delta|^2 */
    acb_srcptr ratio;
    slong i0, i1;
} CorrJob;

static void corr_job_set(CorrJob* J, acb_ptr c, acb_srcptr dc, const Profile* P, slong m,
                         slong wp, acb_srcptr delta, const mag_struct* d2,
                         acb_srcptr ratio, slong i0, slong i1) {
    J->c = c; J->dc = dc; J->P = P; J->m = m; J->wp = wp;
    J->delta = delta; J->d2 = d2; J->ratio = ratio; J->i0 = i0; J->i1 = i1;
}

static void corr_range(void* arg) {
    const CorrJob* J = (const CorrJob*)arg;
    const Profile* P = J->P;
    slong lo = P->lo, wp = J->wp;
    acb_t dt; acb_init(dt);
    mag_t mk; mag_init(mk);
    for (slong i = J->i0; i < J->i1; i++) {
        if (J->delta) {
            /* delta ~ 2^-(wp/2): the product delta e_k' needs only ~p_k - wp/2
             * bits to be accurate to 2^-p_k relative to e_k */
            slong pc = P->p[i] - wp / 2 + 48;
            if (pc < 64) pc = 64;
            if (pc > P->p[i]) pc = P->p[i];
            acb_mul(dt, J->delta, J->dc + i, pc);
            acb_add(J->c + i, J->c + i, dt, P->p[i]);
            /* second-order remainder of the first-order correction:
             * ~ |delta|^2 |e_k''| <~ |delta|^2 (8 (|k - m| + 2))^2 |e_k| */
            slong dk = lo + i - J->m;
            if (dk < 0) dk = -dk;
            acb_get_mag(mk, J->c + i);
            mag_mul(mk, mk, J->d2);
            mag_mul_ui(mk, mk, (ulong)(64 * (dk + 2) * (dk + 2)));
            acb_add_error_mag(J->c + i, mk);
        }
        if (J->ratio) acb_mul(J->c + i, J->c + i, J->ratio, P->p[i]);
    }
    mag_clear(mk);
    acb_clear(dt);
}

/* Newton for nu up to precision `ptop` (about half the final precision; the
 * final pass in mathieu_acb applies the last correction to the coefficients).
 * Precision doubles per level, ending at ptop; g' is carried at half the level's
 * precision -- an approximate derivative with relative error eps keeps Newton
 * at e_{n+1} ~ eps e_n + C e_n^2, still 2^-pw from 2^-(pw/2). Returns 0 on
 * failure. Arithmetic is on midpoints. */
static int refine_nu(acb_t nu, const acb_t a, const acb_t q, slong m,
                     double nur, double nui, double ar, double ai, double qr, double qi,
                     slong ptop, slong loss, int real_nu) {
    acb_t am, qm, g, dg, step, emm, demm;
    acb_init(am); acb_init(qm); acb_init(g); acb_init(dg); acb_init(step);
    acb_init(emm); acb_init(demm);
    acb_get_mid(am, a); acb_get_mid(qm, q);
    mag_t last;
    mag_init(last);
    int ok = 0;
    /* The level schedule, built top-down (ptop, ~ptop/2, ..., <= 96) so that it
     * ends exactly at ptop: doubling up from 64 would end with a level just
     * above a power of two, a near-duplicate of the full-cost level before it.
     * The floor is raised by twice the normalization loss: g' carries the
     * factor sum c_k^2 / c_m^2 (2^-loss of its scale), so at a level of pw bits
     * Newton's steps sit on a noise floor of ~2^-(pw - loss) and could never
     * pass the 2^-(pw/2) test below unless pw > 2 loss. */
    slong lev[64];
    int nlev = 0;
    slong pfloor = 96 + 2 * loss;
    for (slong pw = ptop; ; pw = pw / 2 + 8) {
        lev[nlev++] = pw;
        if (pw <= pfloor || pw / 2 + 8 < pfloor || nlev == 64) break;
    }
    for (int li = nlev - 1; li >= 0; li--) {
        slong pw = lev[li];
        Profile P;
        if (!profile_build(&P, nur, nui, ar, ai, qr, qi, m, 0.0, pw)) break;
        slong n = P.hi - P.lo + 1;
        acb_ptr e = _acb_vec_init(n), de = _acb_vec_init(n);
        int lvl_ok = 0;
        for (int it = 0; it < 40; it++) {
            MillerCtx M;
            mctx_init(&M, am, qm, nu, pw + 16);
            /* g' itself cancels by 2^-loss: carry it with that much more */
            slong pd = pw / 2 + 32 + loss;
            if (pd > pw + 16) pd = pw + 16;
            miller_sweeps(e, de, emm, demm, &M, qm, m, &P, pw + 16, pd);
            miller_g(g, dg, e, de, emm, demm, &M, qm, nu, m, P.lo, pw + 16, pd);
            mctx_clear(&M);
            if (acb_contains_zero(dg)) break;
            acb_div(step, g, dg, pw + 16);
            acb_get_mid(step, step);
            acb_sub(nu, nu, step, pw + 16);
            acb_get_mid(nu, nu);
            if (real_nu) arb_zero(acb_imagref(nu));
            acb_get_mag(last, step);
            if (mag_cmp_2exp_si(last, -(pw / 2 - 8)) <= 0) { lvl_ok = 1; break; }
        }
        _acb_vec_clear(e, n); _acb_vec_clear(de, n);
        profile_clear(&P);
        if (!lvl_ok) break;
        if (li == 0) ok = 1;
    }
    mag_clear(last);
    acb_clear(am); acb_clear(qm); acb_clear(g); acb_clear(dg); acb_clear(step);
    acb_clear(emm); acb_clear(demm);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Far from the real axis: the Bessel-function expansion               */
/* ------------------------------------------------------------------ */
/* The Fourier (Floquet) sum needs every k with |c_k| e^{2|k||Im z|} above the
 * target, about sqrt|q| e^{|Im z|} terms: 7e12 of them at Im z = 30, which is
 * what made MathieuCPrime[1.5, 2.0, 0.7 + 30 I] run forever (audit
 * MathieuCPrime#4). Far from the axis the same coefficients feed a series of
 * cylinder functions that converges FASTER the farther out z is:
 *
 *     Phi_C(w) = sum_k (-1)^k c_k C_{nu+2k}(2 h cos w),   h^2 = q,
 *
 * solves the Mathieu equation y'' + (a - 2q cos 2w) y = 0 for any cylinder
 * function C (DLMF 28.23, in the modified variable cosh z = cos w), and
 * converges for |cos w| > 1 -- geometrically, with ratio 1/|cos w|^2, on the
 * side where the order grows against the argument. Phi_J and Phi_Y are two
 * solutions, so the requested function (any of ce, se and their derivatives)
 * is alpha Phi_J + beta Phi_Y, with alpha, beta fixed by its value and
 * derivative at w0 = Re z + i y0, close enough to the axis for the Fourier sum
 * to be cheap (y0 <= 3) and far enough for the series to converge (y0 >= 1.5,
 * |cos w0| >= sinh 1.5 > 2). Both points lie on one vertical segment, along
 * which 2 h cos w keeps clear of the cut of C (h or -h is chosen so): the
 * identity is an analytic continuation, not an asymptotic one, so the result
 * is exact up to the ball arithmetic that bounds it. */

#define MATHIEU_FAR_TERMS 1024.0   /* Fourier terms above which the far path runs */

/* Beyond this |Im z| the far path declines: the phase of C_mu(X) alone needs
 * |Im z|/log 2 working bits (|X| ~ e^{|Im z|}), and every value is far outside
 * the MPFR exponent range long before (e^{-|Im z|/2} on the imaginary axis,
 * e^{|h| e^{|Im z|}} off it). 2048 costs ~3 s at machine precision. */
#define MATHIEU_FAR_MAX_IM 2048.0

/* Is z far enough out for the Bessel expansion to be the cheaper one? Sets the
 * connection height *y0 (> 0). -1: too far for either path (decline). */
static int far_plan(double qr, double qi, double zi, double* y0) {
    double sq = sqrt(hypot(qr, qi)), ay = fabs(zi);
    if (!(sq > 0.0) || isnan(ay)) return 0;
    if (ay > MATHIEU_FAR_MAX_IM) return -1;
    double c = log(64.0 / sq);
    if (c > 3.0) c = 3.0;
    if (c < 1.5) c = 1.5;
    if (ay <= c + 1.0) return 0;
    if (sq * exp(ay) <= MATHIEU_FAR_TERMS) return 0;
    *y0 = c;
    return 1;
}

/* The far window: the real-axis window around m, widened on both sides until
 * the geometric tail (ratio 1/|cos w0|^2 <= 1/sinh^2 y0) is 2^-(bits+40) below
 * its head, at the uniform precision bits + 32 -- the tail coefficients are
 * multiplied by cylinder functions that grow as fast as they shrink, so they
 * are needed to full RELATIVE accuracy (the Miller start-up error at the
 * window ends decays super-exponentially inward; the last few terms are below
 * the target whatever their error). */
static int profile_far(Profile* P, double nur, double nui, double ar, double ai,
                       double qr, double qi, slong m, slong bits, double y0) {
    if (!profile_build(P, nur, nui, ar, ai, qr, qi, m, 0.0, bits)) return 0;
    double lc = log(sinh(y0));                         /* log |cos w0| lower bound */
    slong E = (slong)(((double)bits + 40.0) * M_LN2 / (2.0 * lc)) + 8;
    slong lo = P->lo < m - E ? P->lo : m - E;
    slong hi = P->hi > m + E ? P->hi : m + E;
    slong n = hi - lo + 1;
    slong* p = (slong*)malloc(sizeof(slong) * (size_t)n);
    if (!p) { profile_clear(P); return 0; }
    for (slong i = 0; i < n; i++) p[i] = bits + 32;
    free(P->p);
    P->p = p; P->lo = lo; P->hi = hi;
    return 1;
}

/* Does the arc of arg(2 h cos(x + i t)), t between ta and tb (same sign), stay
 * off the negative real axis? h enters through its argument `ah`. The point
 * cos x cosh t - i sin x sinh t moves monotonically along one branch of a
 * hyperbola, inside one closed quadrant, so the arc is the interval between
 * the endpoint arguments (no wrap: they share the quadrant). */
static int far_arc_clear(double ah, double x, double ta, double tb) {
    const double PI = 3.14159265358979323846;
    /* divided through by cosh t > 0: no overflow at a large t */
    double fa = atan2(-sin(x) * tanh(ta), cos(x)) + ah;
    double fb = atan2(-sin(x) * tanh(tb), cos(x)) + ah;
    double lo = fa < fb ? fa : fb, hi = fa < fb ? fb : fa;
    /* any odd multiple of pi in [lo - margin, hi + margin]? */
    double margin = 1e-6;
    double j = ceil((lo - margin - PI) / (2.0 * PI));
    return !(PI + 2.0 * PI * j <= hi + margin);
}

/* C_{nu+j}(X) for j = jlo .. jhi (jlo < 0 < jhi) into v[j - jlo], C = J
 * (bessel_y = 0) or Y, from two direct evaluations and the three-term
 * recurrence C_{m-1} + C_{m+1} = (2m/X) C_m, each run in its stable
 * direction: J downward from the top (it is the minimal solution there, and
 * grows as the order falls toward and through -|X| for a non-integer nu), Y
 * outward from j = 0 (it grows away from the oscillatory band either way).
 * An integer nu uses C_{-m} = (-1)^m C_m instead of the negative orders, where
 * J would be recurred in its unstable direction. Four Bessel calls per point
 * where direct evaluation took four per coefficient. */
static void far_orders(acb_ptr v, int bessel_y, const acb_t nu, int integer,
                       const acb_t X, slong jlo, slong jhi, slong prec) {
    slong N = jhi - jlo + 1;
    acb_t m, t, Xi;
    acb_init(m); acb_init(t); acb_init(Xi);
    acb_inv(Xi, X, prec);
    acb_mul_2exp_si(Xi, Xi, 1);                          /* 2/X */
    if (integer) {
        /* orders n0 + j; |order| <= top */
        slong n0 = arf_get_si(arb_midref(acb_realref(nu)), ARF_RND_NEAR);
        slong olo = n0 + jlo, ohi = n0 + jhi;
        slong top = ohi > -olo ? ohi : -olo;
        if (top < 1) top = 1;
        acb_ptr u = _acb_vec_init(top + 1);              /* C_0 .. C_top */
        if (bessel_y) {
            acb_set_si(m, 0); acb_hypgeom_bessel_y(u + 0, m, X, prec);
            acb_set_si(m, 1); acb_hypgeom_bessel_y(u + 1, m, X, prec);
            for (slong o = 1; o < top; o++) {             /* upward */
                acb_mul_si(t, Xi, o, prec);
                acb_mul(u + o + 1, t, u + o, prec);
                acb_sub(u + o + 1, u + o + 1, u + o - 1, prec);
            }
        } else {
            acb_set_si(m, top);     acb_hypgeom_bessel_j(u + top, m, X, prec);
            acb_set_si(m, top - 1); acb_hypgeom_bessel_j(u + top - 1, m, X, prec);
            for (slong o = top - 1; o > 0; o--) {         /* downward */
                acb_mul_si(t, Xi, o, prec);
                acb_mul(u + o - 1, t, u + o, prec);
                acb_sub(u + o - 1, u + o - 1, u + o + 1, prec);
            }
        }
        for (slong j = jlo; j <= jhi; j++) {
            slong o = n0 + j;
            if (o >= 0) acb_set(v + (j - jlo), u + o);
            else if ((-o) & 1) acb_neg(v + (j - jlo), u - o);
            else acb_set(v + (j - jlo), u - o);
        }
        _acb_vec_clear(u, top + 1);
    } else if (bessel_y) {
        slong i0 = -jlo;                                  /* index of j = 0 */
        acb_hypgeom_bessel_y(v + i0, nu, X, prec);
        acb_add_si(m, nu, 1, prec);
        acb_hypgeom_bessel_y(v + i0 + 1, m, X, prec);
        for (slong i = i0 + 1; i < N - 1; i++) {          /* upward */
            acb_add_si(m, nu, jlo + i, prec);
            acb_mul(t, Xi, m, prec);
            acb_mul(v + i + 1, t, v + i, prec);
            acb_sub(v + i + 1, v + i + 1, v + i - 1, prec);
        }
        for (slong i = i0; i > 0; i--) {                  /* downward */
            acb_add_si(m, nu, jlo + i, prec);
            acb_mul(t, Xi, m, prec);
            acb_mul(v + i - 1, t, v + i, prec);
            acb_sub(v + i - 1, v + i - 1, v + i + 1, prec);
        }
    } else {
        acb_add_si(m, nu, jhi, prec);
        acb_hypgeom_bessel_j(v + N - 1, m, X, prec);
        acb_add_si(m, nu, jhi - 1, prec);
        acb_hypgeom_bessel_j(v + N - 2, m, X, prec);
        for (slong i = N - 2; i > 0; i--) {               /* downward */
            acb_add_si(m, nu, jlo + i, prec);
            acb_mul(t, Xi, m, prec);
            acb_mul(v + i - 1, t, v + i, prec);
            acb_sub(v + i - 1, v + i - 1, v + i + 1, prec);
        }
    }
    acb_clear(m); acb_clear(t); acb_clear(Xi);
}

/* Phi_J, Phi_J', Phi_Y, Phi_Y' at w (derivatives in w).
 *
 * `tanform` (an odd function at a characteristic value, integer nu = n): the
 * coefficients are antisymmetric, c_{-n-k} = -c_k, so the plain series above
 * vanishes identically term by term. The solutions are then the second
 * expansion (A&S 20.6, Se_n in the modified variable, tanh z = i tan w)
 *     Psi_C(w) = tan w sum_k (-1)^k (n + 2k) c_k C_{n+2k}(2 h cos w),
 *     Psi_C'   = sec^2 w T + tan w T',  T the sum. */
static void far_phi(acb_t PJ, acb_t dPJ, acb_t PY, acb_t dPY, acb_srcptr c,
                    slong lo, slong n, const acb_t nu, int integer, int tanform,
                    const acb_t h, const acb_t w, slong prec) {
    /* orders nu + j, j = 2 lo - 1 .. 2 hi: C_mu and C_{mu-1} for mu = nu + 2k */
    slong jlo = 2 * lo - 1, jhi = 2 * (lo + n - 1), N = jhi - jlo + 1;
    acb_ptr J = _acb_vec_init(N), Y = _acb_vec_init(N);
    acb_t X, dX, mu, t, r, dj, dy, sw, cw;
    acb_init(X); acb_init(dX); acb_init(mu); acb_init(t); acb_init(r);
    acb_init(dj); acb_init(dy); acb_init(sw); acb_init(cw);
    acb_sin_cos(sw, cw, w, prec);
    acb_mul(X, cw, h, prec); acb_mul_2exp_si(X, X, 1);        /* 2 h cos w  */
    acb_mul(dX, sw, h, prec); acb_mul_2exp_si(dX, dX, 1);
    acb_neg(dX, dX);                                          /* -2 h sin w */
    far_orders(J, 0, nu, integer, X, jlo, jhi, prec);
    far_orders(Y, 1, nu, integer, X, jlo, jhi, prec);
    acb_zero(PJ); acb_zero(dPJ); acb_zero(PY); acb_zero(dPY);
    for (slong i = 0; i < n; i++) {
        slong k = lo + i, ix = 2 * k - jlo;                   /* index of nu + 2k */
        acb_add_si(mu, nu, 2 * k, prec);
        acb_div(r, mu, X, prec);
        acb_set(t, c + i);
        if (tanform) acb_mul(t, t, mu, prec);                 /* (n + 2k) c_k */
        if (k & 1) acb_neg(t, t);                             /* (-1)^k */
        acb_addmul(PJ, t, J + ix, prec);
        acb_addmul(PY, t, Y + ix, prec);
        /* C'_mu = C_{mu-1} - (mu/X) C_mu */
        acb_set(dj, J + ix - 1); acb_submul(dj, r, J + ix, prec);
        acb_set(dy, Y + ix - 1); acb_submul(dy, r, Y + ix, prec);
        acb_addmul(dPJ, t, dj, prec);
        acb_addmul(dPY, t, dy, prec);
    }
    acb_mul(dPJ, dPJ, dX, prec);
    acb_mul(dPY, dPY, dX, prec);
    if (tanform) {
        acb_t tw, s2;
        acb_init(tw); acb_init(s2);
        acb_div(tw, sw, cw, prec);                            /* tan w   */
        acb_sqr(s2, cw, prec); acb_inv(s2, s2, prec);         /* sec^2 w */
        acb_mul(dPJ, dPJ, tw, prec); acb_addmul(dPJ, PJ, s2, prec);
        acb_mul(dPY, dPY, tw, prec); acb_addmul(dPY, PY, s2, prec);
        acb_mul(PJ, PJ, tw, prec);
        acb_mul(PY, PY, tw, prec);
        acb_clear(tw); acb_clear(s2);
    }
    _acb_vec_clear(J, N); _acb_vec_clear(Y, N);
    acb_clear(X); acb_clear(dX); acb_clear(mu); acb_clear(t); acb_clear(r);
    acb_clear(dj); acb_clear(dy); acb_clear(sw); acb_clear(cw);
}

/* The extra working bits the far path needs, added to wp BEFORE the
 * coefficients are computed (their accuracy bounds everything after):
 * the phase of C_mu(X), |X| ~ |h| e^{|Im z|}, and the connection's
 * cancellation -- off the imaginary axis J and Y both grow like e^{|Im X|}
 * with J + iY or J - iY recessive, so at w0 the Wronskian and the numerators
 * of alpha, beta cancel by ~e^{2|Im X0|}, a loss fixed by y0 and q alone,
 * whatever z (52 bits for q = 2, Re z = 0.7) -- plus the ~16 bits the
 * geometric Y tail cancels in the sums at w0. */
static slong far_guard(double qr, double qi, double zr, double zi, double y0) {
    double r = hypot(qr, qi);
    double hr = sqrt(0.5 * (r + qr)), hi = sqrt(0.5 * (r - qr));
    double x = zi < 0.0 ? -zr : zr;
    double imx0 = 2.0 * fabs(hi * cos(x) * cosh(y0)) + 2.0 * fabs(hr * sin(x) * sinh(y0));
    return 32 + (slong)((fabs(zi) + log(1.0 + 2.0 * sqrt(r))) / M_LN2)
              + (slong)(2.0 * imx0 / M_LN2);
}

static int far_eval(acb_t res, acb_srcptr c, slong lo, slong n, const acb_t nu,
                    int integer, const acb_t a, const acb_t q, const acb_t z, int odd, int deriv,
                    slong wp, double y0) {
    const double PI = 3.14159265358979323846;
    /* Im z < 0: ce and se' are even, se and ce' odd */
    double zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    double x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    int flip = zi < 0.0;
    int tanform = integer && odd;     /* antisymmetric c: see far_phi */
    acb_t zz, w0, h, V0, D0, PJ0, dPJ0, PY0, dPY0, PJ, dPJ, PY, dPY, D, al, be;
    acb_init(zz); acb_init(w0); acb_init(h); acb_init(V0); acb_init(D0);
    acb_init(PJ0); acb_init(dPJ0); acb_init(PY0); acb_init(dPY0);
    acb_init(PJ); acb_init(dPJ); acb_init(PY); acb_init(dPY);
    acb_init(D); acb_init(al); acb_init(be);
    int ok = 0;
    if (flip) { acb_neg(zz, z); x = -x; zi = -zi; } else acb_set(zz, z);

    acb_sqrt(h, q, wp);
    double hr = arf_get_d(arb_midref(acb_realref(h)), ARF_RND_NEAR);
    double hi = arf_get_d(arb_midref(acb_imagref(h)), ARF_RND_NEAR);
    double ah = atan2(hi, hr);
    if (!far_arc_clear(ah, x, y0, zi)) {
        acb_neg(h, h);                                /* rotate X by pi */
        ah += ah > 0.0 ? -PI : PI;
        if (!far_arc_clear(ah, x, y0, zi)) goto out;
    }

    {
        /* wp already carries far_guard's bits (the caller added them before
         * computing c, which bounds every accuracy here) */
        slong pf = wp;
        arb_set(acb_realref(w0), acb_realref(zz));
        arb_set_d(acb_imagref(w0), y0);
        mathieu_acb(V0, a, q, w0, odd, 0, pf);
        mathieu_acb(D0, a, q, w0, odd, 1, pf);
        if (!acb_is_finite(V0) || !acb_is_finite(D0)) goto out;

        far_phi(PJ0, dPJ0, PY0, dPY0, c, lo, n, nu, integer, tanform, h, w0, pf);
        far_phi(PJ, dPJ, PY, dPY, c, lo, n, nu, integer, tanform, h, zz, pf);

        /* alpha Phi_J + beta Phi_Y matches (V0, D0) at w0 */
        acb_mul(D, PJ0, dPY0, pf);
        acb_submul(D, PY0, dPJ0, pf);
        if (acb_contains_zero(D)) goto out;
        acb_mul(al, V0, dPY0, pf);
        acb_submul(al, PY0, D0, pf);
        acb_div(al, al, D, pf);
        acb_mul(be, PJ0, D0, pf);
        acb_submul(be, V0, dPJ0, pf);
        acb_div(be, be, D, pf);

        if (deriv) { acb_mul(res, al, dPJ, pf); acb_addmul(res, be, dPY, pf); }
        else       { acb_mul(res, al, PJ, pf);  acb_addmul(res, be, PY, pf); }
        if (flip && (odd ^ deriv)) acb_neg(res, res);
        ok = 1;
    }
out:
    acb_clear(zz); acb_clear(w0); acb_clear(h); acb_clear(V0); acb_clear(D0);
    acb_clear(PJ0); acb_clear(dPJ0); acb_clear(PY0); acb_clear(dPY0);
    acb_clear(PJ); acb_clear(dPJ); acb_clear(PY); acb_clear(dPY);
    acb_clear(D); acb_clear(al); acb_clear(be);
    return ok;
}

void mathieu_acb(acb_t res, const acb_t a, const acb_t q, const acb_t z,
                 int odd, int deriv, slong prec) {
    if (!acb_is_finite(a) || !acb_is_finite(q) || !acb_is_finite(z)) {
        acb_indeterminate(res); return;
    }
    odd = odd ? 1 : 0;
    deriv = deriv ? 1 : 0;
    if (acb_is_zero(q)) { q0_closed_form(res, a, z, odd, deriv, prec); return; }

    double ar = arf_get_d(arb_midref(acb_realref(a)), ARF_RND_NEAR);
    double ai = arf_get_d(arb_midref(acb_imagref(a)), ARF_RND_NEAR);
    double qr = arf_get_d(arb_midref(acb_realref(q)), ARF_RND_NEAR);
    double qi = arf_get_d(arb_midref(acb_imagref(q)), ARF_RND_NEAR);
    double zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    double zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    double nur, nui; int integer, parity, mi, lo_m, hi_m, rindex;
    if (!sf_machine_mathieu_plan_info(ar, ai, qr, qi, &nur, &nui, &integer,
                                      &parity, &mi, &lo_m, &hi_m, &rindex)) {
        acb_indeterminate(res); return;
    }
    (void)lo_m; (void)hi_m;
    /* the other parity at a characteristic value: identically 0 */
    if (integer && ((odd && parity > 0) || (!odd && parity < 0))) {
        acb_zero(res); return;
    }
    slong m = mi;
    int realpq = acb_is_real(a) && acb_is_real(q);
    int real_nu = integer || (realpq && nui == 0.0);

    /* Guard bits for the phase nu z and the coefficient weights (nu+2k). */
    double zabs = hypot(zr, zi), nuabs = hypot(nur, nui);
    slong wp = prec + 24 + (slong)(log2(1.0 + zabs * (1.0 + nuabs)) + 0.5);
    /* the odd functions (se, ce') vanish like z at 0, so a small |z| costs
     * log2(1/|z|) bits of cancellation in the Floquet sum */
    if ((odd ^ deriv) && zabs > 0.0 && zabs < 1.0) wp += (slong)(-log2(zabs) + 1.0);
    /* The normalization sum c_k^2 can cancel (a complex nu at a large q:
     * 2^-34 of sum |c_k|^2 at MathieuC[1, 1000, 0.1`20], audit
     * MathieuCPrime#2). It costs those bits twice: nu is pinned by a g whose
     * derivative carries the same factor, and the coefficients nu's error
     * leaves behind are amplified by it again in S. Before this the ball
     * kernel ran out of bits, S contained 0, and the call declined. */
    slong loss = (slong)ceil(sf_machine_mathieu_plan_loss(ar, ai, qr, qi));
    wp += 2 * loss;
    /* far from the real axis: the Bessel expansion (far_eval), its guard bits */
    double y0 = 0.0;
    int far = far_plan(qr, qi, zi, &y0);
    if (far < 0) { acb_indeterminate(res); return; }
    if (far) wp += far_guard(qr, qi, zr, zi, y0);

    acb_t nu;
    acb_init(nu);
    if (integer) {
        acb_set_si(nu, (slong)floor(nur + 0.5));
    } else {
        acb_set_d_d(nu, nur, real_nu ? 0.0 : nui);
        if (!refine_nu(nu, a, q, m, nur, nui, ar, ai, qr, qi, wp / 2 + 48 + loss, loss,
                       real_nu)) {
            acb_indeterminate(res); acb_clear(nu); return;
        }
    }

    /* the window at full precision, widened for a complex z -- or, far from
     * the real axis, the real-axis window widened for the Bessel expansion */
    Profile P;
    if (!(far ? profile_far(&P, nur, nui, ar, ai, qr, qi, m, wp, y0)
              : profile_build(&P, nur, nui, ar, ai, qr, qi, m, fabs(zi), wp))) {
        acb_indeterminate(res); acb_clear(nu); return;
    }
    slong lo = P.lo, hi = P.hi, n = hi - lo + 1;
    acb_ptr c = _acb_vec_init(n);
    acb_ptr dc = integer ? NULL : _acb_vec_init(n);
    acb_t t, S, s, emm, demm, g, dg, delta, ratio;
    acb_init(t); acb_init(S); acb_init(s); acb_init(ratio);
    acb_init(emm); acb_init(demm); acb_init(g); acb_init(dg); acb_init(delta);

    /* Final pass at wp: Miller sweeps (values at the profile precision, the
     * nu-derivative at wp/2), the last Newton correction delta = -g/g' (nu is
     * already good to ~wp/2 bits, so delta ~ 2^-(wp/2) and its own error is
     * O(delta^2) ~ 2^-wp), and the coefficients corrected to first order,
     * e_k(nu + delta) = e_k + delta e_k' + O(delta^2).
     *
     * Glue: the - side (k < m) is brought onto the + side's scale in the same
     * pass (e_k *= e+_m/e-_m), so c is ONE vector with c_m = e+_m. Its
     * normalization -- sum c_k^2 = 1, Re(c_0) > 0, DLMF 1/sqrt2 for integer nu
     * -- is a single scalar `s` applied to the final sum, not to every
     * coefficient. */
    {
        /* On the MIDPOINTS of a and q, as in refine_nu: their radii are
         * attached at the end (rin below). Fed into g as balls, a's radius
         * entered g only through d_m, not through the continued-fraction
         * ratios whose a-dependence cancels it, so it reached delta divided
         * by g' -- 2^-loss small at a large q -- and swamped every
         * coefficient: MathieuC[777/10, 10000, 3/10] at 30 digits carried a
         * 2^-102 relative radius at 438 bits and declined. */
        acb_t am, qm;
        acb_init(am); acb_init(qm);
        acb_get_mid(am, a); acb_get_mid(qm, q);
        MillerCtx M;
        mctx_init(&M, am, qm, nu, wp);
        slong pd = wp / 2 + 64 + loss;
        if (pd > wp) pd = wp;
        miller_sweeps(c, dc, emm, demm, &M, qm, m, &P, wp, pd);
        if (!integer) {
            miller_g(g, dg, c, dc, emm, demm, &M, qm, nu, m, lo, wp, pd);
            if (acb_contains_zero(dg)) { mctx_clear(&M); acb_clear(am); acb_clear(qm); acb_indeterminate(res); goto done; }
            acb_div(delta, g, dg, wp);
            acb_neg(delta, delta);
            mag_t d2;
            mag_init(d2);
            acb_get_mag(d2, delta);
            mag_mul(d2, d2, d2);                   /* |delta|^2 */
            /* c_m and the - side's e_m first: the glue ratio needs both */
            CorrJob one;
            corr_job_set(&one, c, dc, &P, m, wp, delta, d2, NULL, m - lo, m - lo + 1);
            corr_range(&one);
            acb_addmul(emm, delta, demm, wp);
            acb_div(ratio, c + (m - lo), emm, wp);
            /* the two sides of m: correction (and on the - side the glue) */
            CorrJob J[2];
            corr_job_set(&J[0], c, dc, &P, m, wp, delta, d2, ratio, 0, m - lo);
            corr_job_set(&J[1], c, dc, &P, m, wp, delta, d2, NULL, m - lo + 1, n);
            mathieu_par2(corr_range, &J[0], &J[1], wp >= MATHIEU_PAR_BITS);
            acb_add(nu, nu, delta, wp);
            /* nu's own error: the correction is exact to O(delta^2) */
            mag_mul_2exp_si(d2, d2, 8);
            if (real_nu) arb_add_error_mag(acb_realref(nu), d2);
            else acb_add_error_mag(nu, d2);
            mag_clear(d2);
        } else {
            /* integer nu: no correction, only the glue of the - side */
            acb_div(ratio, c + (m - lo), emm, wp);
            CorrJob J;
            corr_job_set(&J, c, NULL, &P, m, wp, NULL, NULL, ratio, 0, m - lo);
            corr_range(&J);
        }
        mctx_clear(&M);
        acb_clear(am); acb_clear(qm);
    }
    /* S = sum c_k^2: arb's fused dot product (one rounding, the operands' own
     * mantissa lengths), much cheaper than n squarings at the profile. */
    acb_dot(S, NULL, 0, c, 1, c, 1, n, wp);
    if (acb_contains_zero(S)) { acb_indeterminate(res); goto done; }
    acb_sqrt(s, S, wp);
    acb_div(t, c + (0 - lo), s, wp);
    {
        double c0r = arf_get_d(arb_midref(acb_realref(t)), ARF_RND_NEAR);
        double c0i = arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_NEAR);
        /* the tie first: c_0 / s purely imaginary (exactly, for real (a, q)
         * with Re nu = 0) leaves a rounding-level real part whose sign
         * follows the noise in S -- MathieuS[-30, 500, 1/10] flipped between
         * 114 and 208 bits (see sf_machine_mathieu.c) */
        int flip = fabs(c0r) <= 1e-12 * fabs(c0i) ? c0i < 0.0 : c0r < 0.0;
        if (flip) acb_neg(s, s);
    }
    if (integer && rindex >= 0 && realpq) {
        /* the DLMF sign, by the machine plan's rule (audit MathieuC#2/#3,
         * MathieuS#3): sum_k w(nu + 2k) c_k / s > 0, read on midpoints */
        acb_t T, wk; acb_init(T); acb_init(wk);
        slong nu0 = (slong)floor(nur + 0.5);
        for (slong i = 0; i < n; i++) {
            double w = sf_machine_mathieu_sign_weight(rindex, parity, qr < 0.0,
                                                      (long)(nu0 + 2 * (lo + i)));
            acb_mul_si(wk, c + i, (slong)w, 64);
            acb_add(T, T, wk, 64);
        }
        acb_div(T, T, s, 64);
        if (arf_sgn(arb_midref(acb_realref(T))) < 0) acb_neg(s, s);
        acb_clear(T); acb_clear(wk);
    }
    acb_inv(s, s, wp);
    if (integer) {
        acb_t h; acb_init(h);
        acb_set_si(h, 2); acb_rsqrt(h, h, wp);
        acb_mul(s, s, h, wp);
        acb_clear(h);
    }
    /* The midpoint sweeps' error bound (see miller_sweeps): 2^-(p_k - 10)
     * relative, plus the input radii of a and q propagated with a generous
     * condition factor (they are 0 for exact / exactly-converted inputs). A
     * relative bound, so it is attached before the common scale. */
    {
        mag_t rin, mk, e;
        mag_init(rin); mag_init(mk); mag_init(e);
        mag_add(rin, arb_radref(acb_realref(a)), arb_radref(acb_imagref(a)));
        mag_add(rin, rin, arb_radref(acb_realref(q)));
        mag_add(rin, rin, arb_radref(acb_imagref(q)));
        mag_mul_ui(rin, rin, (ulong)(1024 * (n + 1)));
        for (slong i = 0; i < n; i++) {
            acb_get_mag(mk, c + i);
            mag_one(e); mag_mul_2exp_si(e, e, -(P.p[i] - 10));
            mag_add(e, e, rin);
            mag_mul(e, e, mk);
            acb_add_error_mag(c + i, e);
        }
        mag_clear(rin); mag_clear(mk); mag_clear(e);
    }

    int imag_kind = 0;
    if (realpq && !real_nu) {
        /* Real (a, q): ce = N y_e(z) with y_e real and N^2 real, so the value
         * (and ce' = N y_e') is real or purely imaginary; read which from
         * N = ce(0) = s sum c_k, on midpoints. */
        /* in ball arithmetic: the unscaled Miller values can exceed the
         * double range, so scale before reading the parts */
        /* N_C = ce(0) = s sum c_k; N_S = se'(0) = s sum (nu + 2k) c_k */
        acb_t N, wk; acb_init(N); acb_init(wk);
        for (slong i = 0; i < n; i++) {
            if (odd) {
                acb_add_si(wk, nu, 2 * (lo + i), 64);
                acb_mul(wk, wk, c + i, 64);
                acb_add(N, N, wk, 64);
            } else {
                acb_add(N, N, c + i, 64);
            }
        }
        acb_clear(wk);
        acb_mul(N, N, s, 64);
        double vr = arf_get_d(arb_midref(acb_realref(N)), ARF_RND_NEAR);
        double vi = arf_get_d(arb_midref(acb_imagref(N)), ARF_RND_NEAR);
        acb_clear(N);
        imag_kind = fabs(vi) > fabs(vr);
    }
    if (far) {
        if (!far_eval(res, c, lo, n, nu, integer, a, q, z, odd, deriv, wp, y0)) {
            acb_indeterminate(res); goto done;
        }
        /* Real (a, q) on the imaginary axis: the parity rule of the Fourier
         * branch below (the far path never runs on the real axis). */
        if (realpq && arb_is_zero(acb_realref(z))) {
            int nreal = real_nu || !imag_kind;
            int vreal = !(odd ^ deriv) ? nreal : !nreal;
            if (vreal) arb_zero(acb_imagref(res));
            else arb_zero(acb_realref(res));
        }
        acb_set_round(res, res, prec);
        goto done;
    }
    if (deriv) {                                    /* weights (nu + 2k) */
        for (slong i = 0; i < n; i++) {
            acb_add_si(t, nu, 2 * (lo + i), P.p[i]);
            acb_mul(c + i, c + i, t, P.p[i]);
        }
    }

    if (acb_is_real(z) && real_nu && realpq) {
        /* Real Clenshaw: ce = cos(nu z) A - sin(nu z) B (ce': -(sin A + cos B))
         * with A = sum_{k>=0} al_k cos 2kz, al_0 = c_0, al_k = c_k + c_-k,
         *      B = sum_{k>=1} be_k sin 2kz, be_k = c_k - c_-k. */
        slong N = hi > -lo ? hi : -lo;
        arb_ptr al = _arb_vec_init(N + 1), be = _arb_vec_init(N + 1);
        slong* pk = (slong*)malloc(sizeof(slong) * (size_t)(N + 1));
        for (slong k = 0; k <= N; k++) {
            arb_srcptr cp = (k <= hi) ? acb_realref(c + (k - lo)) : NULL;
            arb_srcptr cm = (k >= 1 && -k >= lo) ? acb_realref(c + (-k - lo)) : NULL;
            slong pp = cp ? P.p[k - lo] : 64, pm = cm ? P.p[-k - lo] : 64;
            pk[k] = pp > pm ? pp : pm;
            if (k == 0) { arb_set(al + 0, cp); continue; }
            if (cp) { arb_set(al + k, cp); arb_set(be + k, cp); }
            if (cm) { arb_add(al + k, al + k, cm, pk[k]); arb_sub(be + k, be + k, cm, pk[k]); }
        }
        arb_t th, c2, cs, sn, A, B, ph, sp, cp2;
        arb_init(th); arb_init(c2); arb_init(cs); arb_init(sn);
        arb_init(A); arb_init(B); arb_init(ph); arb_init(sp); arb_init(cp2);
        arb_mul_2exp_si(th, acb_realref(z), 1);
        arb_sin_cos(sn, cs, th, wp);
        arb_mul_2exp_si(c2, cs, 1);
        clenshaw_cs(A, B, al, be, N, pk, c2, cs, sn, wp);
        arb_mul(ph, acb_realref(nu), acb_realref(z), wp);
        arb_sin_cos(sp, cp2, ph, wp);
        if (odd) {
            /* se = sin(nu z) A + cos(nu z) B; se' (weighted) = cos A - sin B */
            if (!deriv) {
                arb_mul(acb_realref(res), sp, A, wp);
                arb_addmul(acb_realref(res), cp2, B, wp);
            } else {
                arb_mul(acb_realref(res), cp2, A, wp);
                arb_submul(acb_realref(res), sp, B, wp);
            }
        } else if (!deriv) {
            arb_mul(acb_realref(res), cp2, A, wp);
            arb_submul(acb_realref(res), sp, B, wp);
        } else {
            arb_mul(acb_realref(res), sp, A, wp);
            arb_addmul(acb_realref(res), cp2, B, wp);
            arb_neg(acb_realref(res), acb_realref(res));
        }
        arb_mul(acb_realref(res), acb_realref(res), acb_realref(s), wp);
        arb_zero(acb_imagref(res));
        free(pk);
        _arb_vec_clear(al, N + 1); _arb_vec_clear(be, N + 1);
        arb_clear(th); arb_clear(c2); arb_clear(cs); arb_clear(sn);
        arb_clear(A); arb_clear(B); arb_clear(ph); arb_clear(sp); arb_clear(cp2);
    } else {
        /* Centred complex Horner (see eval_plan in sf_machine_mathieu.c). */
        acb_t w, wi, Pw, Pwi, Mw, Mwi, e1, e2, A, B, inz;
        acb_init(w); acb_init(wi); acb_init(Pw); acb_init(Pwi); acb_init(Mw);
        acb_init(Mwi); acb_init(e1); acb_init(e2); acb_init(A); acb_init(B);
        acb_init(inz);
        acb_mul_2exp_si(t, z, 1);
        acb_mul_onei(t, t);                          /* 2 i z */
        acb_exp(w, t, wp);
        acb_neg(t, t);
        acb_exp(wi, t, wp);
        slong i0 = -lo;
        acb_set(Pw, c + (n - 1)); acb_set(Pwi, c + (n - 1));
        for (slong j = n - 2; j >= i0; j--) {
            slong pj = P.p[j];
            acb_mul(Pw, Pw, w, pj);   acb_add(Pw, Pw, c + j, pj);
            acb_mul(Pwi, Pwi, wi, pj); acb_add(Pwi, Pwi, c + j, pj);
        }
        acb_set(Mw, c + 0); acb_set(Mwi, c + 0);
        for (slong j = 1; j < i0; j++) {
            slong pj = P.p[j];
            acb_mul(Mw, Mw, w, pj);   acb_add(Mw, Mw, c + j, pj);
            acb_mul(Mwi, Mwi, wi, pj); acb_add(Mwi, Mwi, c + j, pj);
        }
        acb_mul(Mw, Mw, w, wp); acb_mul(Mwi, Mwi, wi, wp);
        acb_mul(inz, nu, z, wp);
        acb_mul_onei(inz, inz);
        acb_exp(e1, inz, wp);
        acb_neg(inz, inz);
        acb_exp(e2, inz, wp);
        acb_add(A, Pw, Mwi, wp); acb_mul(A, A, e1, wp);
        acb_add(B, Pwi, Mw, wp); acb_mul(B, B, e2, wp);
        /* A = F(z), B = F(-z): ce = (A+B)/2, ce' = (i/2)(A-B),
         * se = -(i/2)(A-B), se' = (A+B)/2 (the derivatives with weighted c) */
        int sum_form = odd ? deriv : !deriv;
        if (sum_form) {
            acb_add(res, A, B, wp);
            acb_mul_2exp_si(res, res, -1);
        } else {
            acb_sub(res, A, B, wp);
            if (odd) acb_div_onei(res, res); else acb_mul_onei(res, res);
            acb_mul_2exp_si(res, res, -1);
        }
        acb_mul(res, res, s, wp);
        /* Real (a, q): the value is N y(z), y = y_e, y_e', y_o, y_o' real on the
         * real axis (the equation is real there). On the imaginary axis the
         * equation is real too, so the EVEN ones (y_e, y_o') are real there and
         * the ODD ones (y_e', y_o) are i times a real. */
        if (realpq && (acb_is_real(z) || arb_is_zero(acb_realref(z)))) {
            int nreal = real_nu || !imag_kind;           /* N real, else imaginary */
            int vreal = acb_is_real(z) || !(odd ^ deriv) ? nreal : !nreal;
            if (vreal) arb_zero(acb_imagref(res));
            else arb_zero(acb_realref(res));
        }
        acb_clear(w); acb_clear(wi); acb_clear(Pw); acb_clear(Pwi); acb_clear(Mw);
        acb_clear(Mwi); acb_clear(e1); acb_clear(e2); acb_clear(A); acb_clear(B);
        acb_clear(inz);
    }
    acb_set_round(res, res, prec);

done:
    acb_clear(t); acb_clear(S); acb_clear(s); acb_clear(ratio);
    acb_clear(emm); acb_clear(demm); acb_clear(g); acb_clear(dg); acb_clear(delta);
    _acb_vec_clear(c, n);
    if (dc) _acb_vec_clear(dc, n);
    profile_clear(&P);
    acb_clear(nu);
}

#endif /* USE_FLINT */
