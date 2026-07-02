/*
 * oscde.c — Ooura–Mori double-exponential quadrature for Fourier integrals.
 * See oscde.h for the method and references.
 *
 * Quadrature (a = 0, ω > 0):
 *
 *     I = ∫_0^∞ F(x) dx ≈ Σ_k w_k F(x_k),
 *     x_k = (M/ω) φ(t_k),   w_k = (π/ω) φ'(t_k),   M = π/h,
 *     t_k = k h                (sine:   nodes land on the zeros kπ/ω),
 *     t_k = (k + ½) h          (cosine: zeros at (k+½)π/ω),
 *
 * derived from I = ∫ (M/ω) φ'(t) F(x(t)) dt with the trapezoidal step h and the
 * identity Mh = π (so the h cancels in w_k).  φ, φ', α, β are as in oscde.h.
 */

#include "oscde.h"

#ifdef USE_MPFR

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Scratch temporaries reused across every node of a level (no per-node MPFR
 * alloc/free). */
typedef struct {
    mpfr_t alpha, beta;                /* transform parameters                */
    mpfr_t Mabsc;                      /* M = π/h  (abscissa scale)           */
    mpfr_t Wc;                         /* π        (weight scale)             */
    mpfr_t t, em, ep, psi, epsi, D, Dp, psip, phi, phip;
    mpfr_t x, w, fr, fi;
    mpfr_t sr, si;                     /* running complex sum                 */
    mpfr_t tmp1, tmp2;
} OscDeScratch;

static void oscde_scratch_init(OscDeScratch* s, mpfr_prec_t prec) {
    mpfr_inits2(prec, s->alpha, s->beta, s->Mabsc, s->Wc, s->t, s->em, s->ep,
                s->psi, s->epsi, s->D, s->Dp, s->psip, s->phi, s->phip, s->x,
                s->w, s->fr, s->fi, s->sr, s->si, s->tmp1, s->tmp2,
                (mpfr_ptr)0);
    mpfr_const_pi(s->Wc, MPFR_RNDN);
    mpfr_set_d(s->beta, 0.25, MPFR_RNDN);
}

static void oscde_scratch_clear(OscDeScratch* s) {
    mpfr_clears(s->alpha, s->beta, s->Mabsc, s->Wc, s->t, s->em, s->ep, s->psi,
                s->epsi, s->D, s->Dp, s->psip, s->phi, s->phip, s->x, s->w,
                s->fr, s->fi, s->sr, s->si, s->tmp1, s->tmp2, (mpfr_ptr)0);
}

/* Set (phi, phip) = (φ(t), φ'(t)) for the current s->t.  Returns false when the
 * node is beyond the double-exponential horizon (φ' has underflowed). */
static bool oscde_phi(OscDeScratch* s) {
    if (mpfr_zero_p(s->t)) {
        /* Limits at t = 0 (where 1 - e^{ψ(0)} = 0):  with sig = 2 + α + β,
         *   φ(0)  = 1/sig,
         *   φ'(0) = ((α - β) + sig²) / (2 sig²). */
        mpfr_add(s->tmp1, s->alpha, s->beta, MPFR_RNDN);
        mpfr_add_ui(s->tmp1, s->tmp1, 2, MPFR_RNDN);          /* sig */
        mpfr_ui_div(s->phi, 1, s->tmp1, MPFR_RNDN);           /* φ(0) */
        mpfr_sub(s->tmp2, s->alpha, s->beta, MPFR_RNDN);      /* α-β */
        mpfr_mul(s->Dp, s->tmp1, s->tmp1, MPFR_RNDN);         /* sig² */
        mpfr_add(s->tmp2, s->tmp2, s->Dp, MPFR_RNDN);         /* (α-β)+sig² */
        mpfr_mul_ui(s->Dp, s->Dp, 2, MPFR_RNDN);              /* 2 sig² */
        mpfr_div(s->phip, s->tmp2, s->Dp, MPFR_RNDN);         /* φ'(0) */
        return true;
    }
    mpfr_exp(s->ep, s->t, MPFR_RNDN);                         /* e^t  */
    mpfr_ui_div(s->em, 1, s->ep, MPFR_RNDN);                  /* e^-t */
    /* ψ = -2t - α(1 - e^-t) - β(e^t - 1). */
    mpfr_mul_si(s->psi, s->t, -2, MPFR_RNDN);
    mpfr_ui_sub(s->tmp1, 1, s->em, MPFR_RNDN);
    mpfr_mul(s->tmp1, s->tmp1, s->alpha, MPFR_RNDN);
    mpfr_sub(s->psi, s->psi, s->tmp1, MPFR_RNDN);
    mpfr_sub_ui(s->tmp2, s->ep, 1, MPFR_RNDN);
    mpfr_mul(s->tmp2, s->tmp2, s->beta, MPFR_RNDN);
    mpfr_sub(s->psi, s->psi, s->tmp2, MPFR_RNDN);
    if (mpfr_cmp_d(s->psi, 1.0e6) > 0) return false;          /* DE horizon */
    mpfr_exp(s->epsi, s->psi, MPFR_RNDN);                     /* e^ψ */
    mpfr_ui_sub(s->D, 1, s->epsi, MPFR_RNDN);                /* D = 1 - e^ψ */
    if (mpfr_zero_p(s->D)) return false;
    mpfr_div(s->phi, s->t, s->D, MPFR_RNDN);                 /* φ = t/D */
    /* ψ' = -2 - α e^-t - β e^t. */
    mpfr_mul(s->psip, s->alpha, s->em, MPFR_RNDN);
    mpfr_mul(s->tmp1, s->beta, s->ep, MPFR_RNDN);
    mpfr_add(s->psip, s->psip, s->tmp1, MPFR_RNDN);
    mpfr_add_ui(s->psip, s->psip, 2, MPFR_RNDN);
    mpfr_neg(s->psip, s->psip, MPFR_RNDN);                   /* ψ' */
    /* D' = -e^ψ ψ';  φ' = (D - t D') / D². */
    mpfr_mul(s->Dp, s->epsi, s->psip, MPFR_RNDN);
    mpfr_neg(s->Dp, s->Dp, MPFR_RNDN);
    mpfr_mul(s->tmp1, s->t, s->Dp, MPFR_RNDN);
    mpfr_sub(s->tmp1, s->D, s->tmp1, MPFR_RNDN);
    mpfr_mul(s->tmp2, s->D, s->D, MPFR_RNDN);
    mpfr_div(s->phip, s->tmp1, s->tmp2, MPFR_RNDN);          /* φ' */
    return true;
}

/* Accumulate the node at s->t into (sr, si); returns the term magnitude as a
 * double (0 when the node is beyond the horizon or non-numeric). */
static double oscde_node(OscDeScratch* s, DeQuadSampleMPFR f, void* ctx,
                         double omega) {
    if (!oscde_phi(s)) return 0.0;
    mpfr_mul(s->x, s->Mabsc, s->phi, MPFR_RNDN);   /* x = (M/ω) φ */
    mpfr_div_d(s->x, s->x, omega, MPFR_RNDN);
    if (!f(ctx, s->x, s->fr, s->fi)) return 0.0;
    if (!mpfr_number_p(s->fr) || !mpfr_number_p(s->fi)) return 0.0;
    mpfr_mul(s->w, s->Wc, s->phip, MPFR_RNDN);     /* w = (π/ω) φ' */
    mpfr_div_d(s->w, s->w, omega, MPFR_RNDN);
    mpfr_mul(s->fr, s->fr, s->w, MPFR_RNDN);
    mpfr_mul(s->fi, s->fi, s->w, MPFR_RNDN);
    mpfr_add(s->sr, s->sr, s->fr, MPFR_RNDN);
    mpfr_add(s->si, s->si, s->fi, MPFR_RNDN);
    return hypot(mpfr_get_d(s->fr, MPFR_RNDN), mpfr_get_d(s->fi, MPFR_RNDN));
}

/* One level: (offset) trapezoidal sum with step h.  Writes (out_re,out_im). */
static void oscde_level(OscDeScratch* s, DeQuadSampleMPFR f, void* ctx,
                        double omega, OscDeKind kind, double h, double reltol,
                        mpfr_t out_re, mpfr_t out_im) {
    /* α (double is fine: it only tunes the convergence rate — the continuous
     * substitution ∫(M/ω)φ'F dt = ∫F dx holds for any monotone φ_α, so α's
     * precision never limits the h→0 limit).  The abscissa scale M = π/h,
     * however, multiplies φ directly: a double-rounded M makes x_k = kπ(1+ε)
     * and floors sin(x_k) at ~1e-16 in the tail (capping accuracy AND stopping
     * the tail from decaying), so it must be formed at full precision. */
    double M = M_PI / h;
    double alpha = 0.25 / sqrt(1.0 + M * log(1.0 + M) / (2.0 * M_PI));
    mpfr_set_d(s->alpha, alpha, MPFR_RNDN);
    mpfr_div_d(s->Mabsc, s->Wc, h, MPFR_RNDN);   /* Mabsc = π/h (Wc = π) */
    mpfr_set_zero(s->sr, 1);
    mpfr_set_zero(s->si, 1);

    double off = (kind == OSCDE_COS) ? 0.5 : 0.0;
    /* Truncate a term once it falls a little below reltol·(peak term); the
     * dropped double-exponential tail is then well under the target error. */
    double drop = reltol * 1.0e-2;
    double maxterm = 0.0;
    const int QUIT = 4;
    const int HARD = 500000;

    /* Rightward march (k = 0, 1, 2, …). */
    int streak = 0;
    for (int k = 0; k < HARD; k++) {
        mpfr_set_d(s->t, (k + off) * h, MPFR_RNDN);
        double tm = oscde_node(s, f, ctx, omega);
        if (tm > maxterm) maxterm = tm;
        if (tm <= drop * maxterm) { if (++streak >= QUIT) break; }
        else streak = 0;
    }
    /* Leftward march (k = -1, -2, …). */
    streak = 0;
    for (int k = -1; k > -HARD; k--) {
        mpfr_set_d(s->t, (k + off) * h, MPFR_RNDN);
        double tm = oscde_node(s, f, ctx, omega);
        if (tm > maxterm) maxterm = tm;
        if (tm <= drop * maxterm) { if (++streak >= QUIT) break; }
        else streak = 0;
    }

    mpfr_set(out_re, s->sr, MPFR_RNDN);
    mpfr_set(out_im, s->si, MPFR_RNDN);
}

bool oscde_fourier_mpfr(DeQuadSampleMPFR f, void* ctx, double omega,
                        OscDeKind kind, long bits, double reltol,
                        int max_levels, mpfr_t out_re, mpfr_t out_im,
                        double* abserr) {
    if (!(omega > 0.0)) return false;
    if (reltol <= 0.0) reltol = 1e-10;
    if (max_levels < 3) max_levels = 3;

    OscDeScratch s;
    oscde_scratch_init(&s, (mpfr_prec_t)bits);

    mpfr_t cur_re, cur_im, prev_re, prev_im;
    mpfr_inits2((mpfr_prec_t)bits, cur_re, cur_im, prev_re, prev_im,
                (mpfr_ptr)0);

    double h = 1.0;                       /* halved each level */
    double err = INFINITY;
    bool conv = false, have_prev = false;

    mpfr_set_zero(out_re, 1);
    mpfr_set_zero(out_im, 1);

    for (int level = 0; level < max_levels; level++, h *= 0.5) {
        oscde_level(&s, f, ctx, omega, kind, h, reltol, cur_re, cur_im);
        mpfr_set(out_re, cur_re, MPFR_RNDN);
        mpfr_set(out_im, cur_im, MPFR_RNDN);
        if (have_prev) {
            /* L1 level-to-level change and |S|₁, both formed in MPFR so the
             * comparison stays meaningful at the high-precision tail where a
             * double subtraction would underflow. */
            mpfr_sub(s.tmp1, cur_re, prev_re, MPFR_RNDN);
            mpfr_sub(s.w,    cur_im, prev_im, MPFR_RNDN);
            mpfr_abs(s.tmp1, s.tmp1, MPFR_RNDN);
            mpfr_abs(s.w,    s.w,    MPFR_RNDN);
            mpfr_add(s.tmp1, s.tmp1, s.w, MPFR_RNDN);         /* |ΔS|₁ */
            err = mpfr_get_d(s.tmp1, MPFR_RNDN);
            mpfr_abs(s.tmp2, cur_re, MPFR_RNDN);
            mpfr_abs(s.w,    cur_im, MPFR_RNDN);
            mpfr_add(s.tmp2, s.tmp2, s.w, MPFR_RNDN);         /* |S|₁ */
            mpfr_mul_d(s.tmp2, s.tmp2, reltol, MPFR_RNDN);
            if (level >= 2 && mpfr_cmp(s.tmp1, s.tmp2) <= 0) conv = true;
        }
        mpfr_set(prev_re, cur_re, MPFR_RNDN);
        mpfr_set(prev_im, cur_im, MPFR_RNDN);
        have_prev = true;
        if (conv) break;
    }

    *abserr = err;
    mpfr_clears(cur_re, cur_im, prev_re, prev_im, (mpfr_ptr)0);
    oscde_scratch_clear(&s);
    return conv;
}

#endif /* USE_MPFR */
