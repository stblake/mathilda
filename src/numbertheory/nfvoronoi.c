/*
 * nfvoronoi.c
 * -----------
 * Voronoi's algorithm for the fundamental unit of a RANK-1 complex cubic
 * field K = Q(rho), signature (r1, r2) = (1, 1).  This is the fallback the
 * coefficient-box search in nfunits.c cannot cover: a fundamental unit's
 * coordinates grow with the regulator, so a bounded box misses the units of
 * large-regulator fields (Q(cbrt15) reg 9.7 has a coordinate 30; Q(cbrt41)
 * reg 56 has 12-digit coordinates).  Voronoi walks the CHAIN OF MINIMA of
 * O_K, whose length is polynomial in the regulator, and returns the first
 * unit it reaches.
 *
 * Contract (identical to the rest of Gate 2): this routine only PROPOSES the
 * fundamental unit's integer coordinates.  Correctness is the caller's exact
 * p-saturation + |N| == 1 certifier (nfunits.c); a wrong proposal is caught
 * there and DECLINEs -- never a wrong answer.  Everything geometric is done
 * in arb/acb interval arithmetic with adaptive precision; every strict
 * comparison that gates the walk is interval-safe (widen precision on a tie).
 *
 * Method (Voronoi; Fung, "Computational Problems in Complex Cubic Fields";
 * Buchmann-Pohst-Schmettow, Math. Comp. 53 (1989)).  Embed O_K in R^3 by the
 * Minkowski map psi(a) = (sigma1(a), Re sigma2(a), Im sigma2(a)); the norm is
 * N(a) = sigma1(a)*(Re^2 + Im^2) = sigma1(a)*|sigma2(a)|^2.  A nonzero minimum
 * theta of O_K has no nonzero b in O_K with |sigma1(b)| < |sigma1(theta)| AND
 * |sigma2(b)| < |sigma2(theta)|.  The directed chain theta_0 = 1, theta_1, ...
 * takes at each step the SECOND-KIND adjacent minimum:
 *      theta' = argmin { sigma1(phi) : phi in O_K,
 *                        sigma1(phi) > sigma1(theta), |sigma2(phi)| < |sigma2(theta)| }.
 * sigma1 strictly increases along the chain; the first theta_k (k >= 1) with
 * |N(theta_k)| == 1 is the fundamental unit eps (> 1 in the real embedding).
 *
 * Well-conditioned neighbour search: the adjacent minimum theta' = psi*theta
 * with psi a SHORT vector of the fractional ideal theta^{-1} O_K.  We LLL-
 * reduce that ideal's Minkowski lattice each step (so psi is a small integer
 * combination of the reduced basis), recovering the reduced basis' O_K
 * preimages g_m via an appended identity block: rows
 *      [ round(2^P * mink(rho^j / theta)) (3) | e_j (3) ],
 * LLL over the (huge, metric-dominated) integer lattice records the unimodular
 * transform T in the identity columns, and g_m has power-basis coordinates
 * equal to row m of T.  Then theta' = sum_m e_m g_m has exact integer
 * coordinates -- no field division, no ideal-HNF bookkeeping.  theta_k is
 * tracked exactly in mpz coordinates (FLINT/GMP absorb the <= 12-digit values).
 *
 * Requires FLINT.  Without it nf_voronoi_unit_cubic11 returns false.
 */
#include <stdio.h>       /* before gmp.h so gmp_fprintf is declared */
#include "numberfield_internal.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef USE_FLINT

bool nf_voronoi_unit_cubic11(struct NumberField* K, mpz_t* out) {
    (void)K; (void)out; return false;
}

#else

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/fmpz.h>
#include "linalg.h"     /* lll_reduce_q */

static int vdbg(void) { static int v = -1; if (v < 0) v = getenv("THUE_DEBUG") ? 1 : 0; return v; }

/* sigma1(coords) as a real arb (the real embedding; imag part ~ 0). */
static void emb_real(const struct NumberField* K, const mpz_t* c, arb_t out) {
    acb_t v; acb_init(v);
    nf_embed_int(K, c, 0, v);
    arb_set(out, acb_realref(v));
    acb_clear(v);
}

/* |sigma2(coords)| as a real arb (modulus of the complex embedding). */
static void emb_absc(const struct NumberField* K, const mpz_t* c, arb_t out) {
    acb_t v; acb_init(v);
    nf_embed_int(K, c, 1, v);
    acb_abs(out, v, K->prec + 16);
    acb_clear(v);
}

/* Bit-magnitude of a positive arb (approx log2), clamped to >= 1. */
static long arb_log2_mag(const arb_t x) {
    mag_t mg; mag_init(mg);
    arb_get_mag(mg, x);
    double d = mag_get_d(mg);
    mag_clear(mg);
    if (!(d > 1.0) || d != d) return 1;   /* <=1 or NaN */
    if (d > 1e300) return 1000;           /* far beyond any cubic we handle */
    long b = (long)ceil(log2(d));
    return b < 1 ? 1 : b;
}

/*
 * One chain step: replace th[0..2] (an O_K minimum) by the second-kind
 * adjacent minimum.  *wp is the working arb precision (grown as sigma1
 * increases).  Returns true on success, false to DECLINE (no valid neighbour
 * found within the enumeration budget, or precision blew past the cap).
 */
static bool voronoi_step(struct NumberField* K, mpz_t* th, slong* wp) {
    const int deg = 3;

    /* --- pick a precision that keeps the 2^P-scaled metric faithful --- */
    /* A first cheap embedding at the current precision to size sigma1. */
    if (!nf_ensure_prec(K, *wp)) return false;
    arb_t curS1, curS2; arb_init(curS1); arb_init(curS2);
    emb_real(K, (const mpz_t*)th, curS1);
    emb_absc(K, (const mpz_t*)th, curS2);
    long s1bits = arb_log2_mag(curS1);

    /* P: scale so the smallest scaled entry (~2^-s1bits) survives rounding;
     * WP: enough headroom above the largest (~2^(P+deg*mag(rho))). */
    long P = s1bits + 80;
    if (P < 128) P = 128;
    slong needWP = (slong)(P + 128);
    if (needWP > (slong)200000) { arb_clear(curS1); arb_clear(curS2); return false; }
    if (needWP > *wp) {
        *wp = needWP;
        if (!nf_ensure_prec(K, *wp)) { arb_clear(curS1); arb_clear(curS2); return false; }
        emb_real(K, (const mpz_t*)th, curS1);   /* recompute at the new prec */
        emb_absc(K, (const mpz_t*)th, curS2);
    }
    slong prec = K->prec + 16;

    /* --- basis embeddings sigma_i(rho^j), j = 0..2 --- */
    /* rho^j is the coordinate vector e_j; embed it directly. */
    acb_t bemb1[3], bemb2[3];               /* real & complex embedding of rho^j */
    for (int j = 0; j < deg; j++) {
        acb_init(bemb1[j]); acb_init(bemb2[j]);
        mpz_t ej[3]; for (int i = 0; i < deg; i++) mpz_init_set_ui(ej[i], i == j ? 1u : 0u);
        nf_embed_int(K, (const mpz_t*)ej, 0, bemb1[j]);
        nf_embed_int(K, (const mpz_t*)ej, 1, bemb2[j]);
        for (int i = 0; i < deg; i++) mpz_clear(ej[i]);
    }

    /* sigma1(theta), sigma2(theta) as acb for the per-embedding scaling. */
    acb_t S1c, S2c; acb_init(S1c); acb_init(S2c);
    nf_embed_int(K, (const mpz_t*)th, 0, S1c);
    nf_embed_int(K, (const mpz_t*)th, 1, S2c);

    /* --- persistent scratch across the U-growth loop --- */
    mpq_t* rows = malloc(sizeof(mpq_t) * (size_t)(deg * 6));
    for (int i = 0; i < deg * 6; i++) mpq_init(rows[i]);
    mpz_t T[3][3];
    for (int mm = 0; mm < deg; mm++) for (int c = 0; c < deg; c++) mpz_init(T[mm][c]);
    mpz_t cand[3], best[3];
    for (int i = 0; i < deg; i++) { mpz_init(cand[i]); mpz_init(best[i]); }
    arb_t re, im, scaled, candS1, candS2, bestS1;
    arb_init(re); arb_init(im); arb_init(scaled);
    arb_init(candS1); arb_init(candS2); arb_init(bestS1);
    acb_t z; acb_init(z);
    fmpz_t zi; fmpz_init(zi);
    mpz_t tmp; mpz_init(tmp);
    bool have_best = false;
    int found_u = -1, dep = 0;

    /*
     * The second-kind neighbour minimises sigma1 SUBJECT TO |sigma2| < 1 (in the
     * theta^{-1}O_K frame) -- an anisotropic box, thin in sigma2, extended in
     * sigma1.  It is NOT a short Minkowski vector, so we rescale the real column
     * by 1/U (U = 2^u) to make the box isotropic, then LLL + a small enumeration
     * find it.  Grow U until the box (sigma1 < U) first captures a valid
     * neighbour; the minimal-sigma1 lattice point in that box is the adjacent
     * minimum (min over (1,U] == global min over (1,inf) once the box is
     * non-empty).  U stays bounded by the field's Voronoi constant, so this is a
     * handful of iterations per step.
     */
    const int UMAX_U = 60;
    for (int u = 1; u <= UMAX_U && !have_best; u++) {
        bool ok = true;
        for (int j = 0; j < deg; j++) {
            /* col 0: Re(sigma1(rho^j)/sigma1(theta)) scaled by 2^(P-u) (== /U) */
            acb_div(z, bemb1[j], S1c, prec);
            arb_set(scaled, acb_realref(z));
            arb_mul_2exp_si(scaled, scaled, P - u);
            arf_get_fmpz(zi, arb_midref(scaled), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * 6 + 0], tmp);

            /* cols 1,2: Re/Im of sigma2(rho^j)/sigma2(theta) scaled by 2^P */
            acb_div(z, bemb2[j], S2c, prec);
            arb_set(re, acb_realref(z)); arb_set(im, acb_imagref(z));
            arb_mul_2exp_si(re, re, P);
            arf_get_fmpz(zi, arb_midref(re), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * 6 + 1], tmp);
            arb_mul_2exp_si(im, im, P);
            arf_get_fmpz(zi, arb_midref(im), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * 6 + 2], tmp);

            /* cols 3..5: identity, passively records the unimodular transform */
            for (int c = 0; c < deg; c++) mpq_set_si(rows[j * 6 + 3 + c], c == j ? 1 : 0, 1);
        }

        /* LLL: metric block (~2^P) dominates; identity block -> transform T. */
        dep = lll_reduce_q(rows, deg, 6, NULL);
        if (dep != 0) continue;
        for (int m = 0; m < deg && ok; m++)
            for (int c = 0; c < deg; c++) {
                mpq_ptr e = rows[m * 6 + 3 + c];
                if (mpz_cmp_ui(mpq_denref(e), 1) != 0) { ok = false; break; }
                mpz_set(T[m][c], mpq_numref(e));   /* g_m power-basis coords */
            }
        if (!ok) continue;

        /* enumerate psi = sum_m e_m g_m over |e_m| <= E; keep the minimal-sigma1
         * point with sigma1 > sigma1(theta) and |sigma2| < |sigma2(theta)|. */
        for (int E = 1; E <= 4; E++) {
            int W = 2 * E + 1;
            long total = 1; for (int i = 0; i < deg; i++) total *= W;
            for (long idx = 0; idx < total; idx++) {
                long t = idx; int e[3]; bool nz = false;
                for (int i = 0; i < deg; i++) { e[i] = (int)(t % W) - E; t /= W; if (e[i]) nz = true; }
                if (!nz) continue;
                if (E > 1) { bool hi = false; for (int i = 0; i < deg; i++) if (abs(e[i]) == E) hi = true; if (!hi) continue; }

                for (int c = 0; c < deg; c++) {
                    mpz_set_ui(cand[c], 0);
                    for (int m = 0; m < deg; m++) {
                        if (e[m] > 0)      mpz_addmul_ui(cand[c], T[m][c], (unsigned long)e[m]);
                        else if (e[m] < 0) mpz_submul_ui(cand[c], T[m][c], (unsigned long)(-e[m]));
                    }
                }

                emb_real(K, (const mpz_t*)cand, candS1);
                if (!arb_lt(curS1, candS1)) continue;             /* need sigma1 > cur */
                emb_absc(K, (const mpz_t*)cand, candS2);
                if (!arb_lt(candS2, curS2)) continue;             /* need |sigma2| < cur */
                if (!have_best || arb_lt(candS1, bestS1)) {
                    for (int c = 0; c < deg; c++) mpz_set(best[c], cand[c]);
                    arb_set(bestS1, candS1);
                    have_best = true;
                }
            }
        }
        if (have_best) found_u = u;      /* first U whose box is non-empty wins */
    }

    if (vdbg()) {
        fprintf(stderr, "  [voro-step] P=%ld wp=%ld found_u=%d have_best=%d",
                P, (long)*wp, found_u, (int)have_best);
        if (have_best) {
            double d1 = arf_get_d(arb_midref(bestS1), ARF_RND_NEAR);
            gmp_fprintf(stderr, " -> th'=[%Zd,%Zd,%Zd] s1=%.6g", best[0], best[1], best[2], d1);
        }
        fprintf(stderr, "\n");
    }

    if (have_best) for (int c = 0; c < deg; c++) mpz_set(th[c], best[c]);

    /* --- cleanup --- */
    for (int i = 0; i < deg * 6; i++) mpq_clear(rows[i]);
    free(rows);
    for (int j = 0; j < deg; j++) { acb_clear(bemb1[j]); acb_clear(bemb2[j]); }
    acb_clear(S1c); acb_clear(S2c); acb_clear(z);
    arb_clear(re); arb_clear(im); arb_clear(scaled);
    arb_clear(curS1); arb_clear(curS2);
    arb_clear(candS1); arb_clear(candS2); arb_clear(bestS1);
    fmpz_clear(zi); mpz_clear(tmp);
    for (int i = 0; i < deg; i++) { mpz_clear(cand[i]); mpz_clear(best[i]); }
    for (int mm = 0; mm < deg; mm++) for (int c = 0; c < deg; c++) mpz_clear(T[mm][c]);

    return have_best;
}

bool nf_voronoi_unit_cubic11(struct NumberField* K, mpz_t* out) {
    if (K->deg != 3 || K->r1 != 1 || K->r2 != 1) return false;
    const int deg = 3;

    mpz_t th[3];
    mpz_init_set_ui(th[0], 1); mpz_init_set_ui(th[1], 0); mpz_init_set_ui(th[2], 0);

    slong wp = 320;
    bool found = false;
    const int MAXSTEP = 50000;
    mpz_t nrm; mpz_init(nrm);

    for (int step = 0; step < MAXSTEP && !found; step++) {
        if (!voronoi_step(K, th, &wp)) { if (vdbg()) fprintf(stderr, "  [voro] step %d FAILED\n", step); break; }
        bool isint = nf_norm_int(K, (const mpz_t*)th, nrm);
        if (vdbg() && (step < 3 || (step % 50) == 0))
            gmp_fprintf(stderr, "  [voro] step %d N=%Zd\n", step, nrm);
        if (isint && mpz_cmpabs_ui(nrm, 1) == 0) {
            for (int i = 0; i < deg; i++) mpz_set(out[i], th[i]);
            found = true;
            if (vdbg()) gmp_fprintf(stderr, "  [voro] UNIT at step %d: [%Zd,%Zd,%Zd]\n", step, th[0], th[1], th[2]);
        }
    }

    mpz_clear(nrm);
    for (int i = 0; i < deg; i++) mpz_clear(th[i]);
    return found;
}

#endif /* USE_FLINT */
