/*
 * nfvoronoi2.c
 * ------------
 * Rank-2 unit finding by directed Voronoi minima walks, for signature (2,1)
 * quartics K = Q(theta) with theta^2 in a real quadratic subfield -- the case
 * the coefficient-box search in nfunits.c and the rank-1 chain in nfvoronoi.c
 * cannot cover.  The archetype is Q(10^{1/4}) (x^4-10y^4 == +-1), unit rank 2,
 * regulator 25.3: its subfield unit eps = theta^2-3 is small (found by the box)
 * but the SECOND fundamental unit is intrinsically large in the Minkowski
 * embedding (|sigma0| ~ 1036), unreachable by any box or short-vector search.
 *
 * Contract (identical to nfvoronoi.c / the rest of Gate 2): this routine only
 * PROPOSES candidate units' integer coordinates.  Correctness is the caller's
 * exact p-saturation + |N|==1 certifier (nfunits.c): a wrong or non-fundamental
 * proposal is caught there and DECLINEs -- never a wrong answer.  So the walk
 * may be heuristic; only the certifier must be exact (it is).
 *
 * Geometry.  The three archimedean places are sigma0, sigma1 (real, theta ->
 * +-10^{1/4}) and sigma2 (complex).  N(a) = sigma0(a)*sigma1(a)*|sigma2(a)|^2.
 * In the log plane (log|sigma0|, log|sigma1|) the unit lattice is 2-D: a unit
 * of the quadratic subfield F=Q(theta^2) has sigma0 = sigma1 (both restrict to
 * the same real place of F), so eps lies along the (1,1) direction; the
 * "relative" unit eta (N_{K/F}(eta)=+-1, so |sigma2(eta)|=1 and
 * sigma0(eta)sigma1(eta)=+-1) lies along (1,-1).  {eps, eta} is therefore a
 * full-rank system.  We reach each by a DIRECTED minima chain:
 *   chain (1,1): grow sigma0 AND sigma1 together, shrink |sigma2|   -> eps
 *   chain (1,-1): grow sigma0, shrink |sigma1|                      -> eta
 * Each step takes the adjacent minimum of O_K in that direction, found by
 * LLL-reducing the theta^{-1}O_K Minkowski lattice (per-column rescaled so the
 * target box is isotropic) and enumerating short combinations -- the exact
 * rank-1 technique of nfvoronoi.c, with the gate generalised to the direction.
 * theta_k is tracked exactly in mpz (power-basis / O_K-basis numerator coords).
 *
 * STATUS: working for the archetype, one chain of two.  The (1,-1) "relative"
 * chain CONVERGES (`chain_first_unit(K, -1, ...)` reaches a unit whose log lies
 * off the subfield's (1,1) direction -- e.g. Q(10^{1/4}): [-1597,898,-505,284]
 * at step 3).  The caller (nfunits.c) pairs it with the subfield unit already
 * in its box-candidate list and p-saturation certifies the pair: for Q(10^{1/4})
 * the regulator comes out 25.2535 (== PARI bnfinit reg) and `x^4-10y^4==+-1`
 * now solves.  The (1,1) subfield chain is NOT walked (the box supplies that
 * unit); its directed gate does not converge, and is left for the full 2-D
 * Voronoi complex.  A field where the (1,-1) chain does not converge blows its
 * coords past the guard and DECLINEs fast (contract-safe -- p-saturation, not
 * the walk, is the proof, so a bad proposal is never a wrong answer).  Extending
 * coverage beyond this archetype = the remaining M4 work
 * (docs/design/thue_completion_plan.md).
 *
 * Requires FLINT.  Without it nf_voronoi_units_sig21 returns 0.
 */
#include <stdio.h>       /* before gmp.h so gmp_fprintf is declared */
#include "numberfield_internal.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef USE_FLINT

int nf_voronoi_units_sig21(struct NumberField* K, mpz_t** out) {
    (void)K; (void)out; return 0;
}

#else

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/fmpz.h>
#include "linalg.h"     /* lll_reduce_q */

static int v2dbg(void) { static int v = -1; if (v < 0) v = getenv("THUE_DEBUG") ? 1 : 0; return v; }

/* |sigma_k(coords)| as a real arb (k in [0, r1+r2)); modulus for the complex
 * place, absolute value for a real one. */
static void emb_abs(const struct NumberField* K, const mpz_t* c, int k, arb_t out) {
    acb_t v; acb_init(v);
    nf_embed_int(K, c, k, v);
    acb_abs(out, v, K->prec + 16);
    acb_clear(v);
}

/* Bit-magnitude of a positive arb (~log2), clamped to >= 1. */
static long arb_log2_mag(const arb_t x) {
    mag_t mg; mag_init(mg);
    arb_get_mag(mg, x);
    double d = mag_get_d(mg);
    mag_clear(mg);
    if (!(d > 1.0) || d != d) return 1;
    if (d > 1e300) return 1000;
    long b = (long)ceil(log2(d));
    return b < 1 ? 1 : b;
}

/*
 * One directed chain step for signature (2,1), deg 4.  Replace th[0..3] (an
 * O_K minimum, power-basis numerator coords) by the adjacent minimum in the
 * requested direction:
 *   dir = +1 : grow sigma0, grow sigma1  (the (1,1) / subfield direction)
 *   dir = -1 : grow sigma0, shrink sigma1 (the (1,-1) / relative direction)
 * In both cases sigma0 strictly increases; the OTHER real place and the complex
 * place are held small (rescaled to make the target box isotropic).  *wp is the
 * working precision, grown as sigma0 increases.  Returns true on success, false
 * to DECLINE (no valid neighbour within budget / precision capped).
 */
static bool step_sig21(struct NumberField* K, mpz_t* th, int dir, slong* wp) {
    const int deg = 4;
    const int m   = 3;                 /* r1 + r2 embeddings: 0,1 real; 2 complex */

    if (!nf_ensure_prec(K, *wp)) return false;
    arb_t curA[3];                     /* |sigma_k(th)|, k=0,1,2 */
    for (int k = 0; k < m; k++) { arb_init(curA[k]); emb_abs(K, (const mpz_t*)th, k, curA[k]); }
    long s0bits = arb_log2_mag(curA[0]);

    long P = s0bits + 80;
    if (P < 128) P = 128;
    slong needWP = (slong)(P + 160);
    if (needWP > (slong)200000) { for (int k = 0; k < m; k++) arb_clear(curA[k]); return false; }
    if (needWP > *wp) {
        *wp = needWP;
        if (!nf_ensure_prec(K, *wp)) { for (int k = 0; k < m; k++) arb_clear(curA[k]); return false; }
        for (int k = 0; k < m; k++) emb_abs(K, (const mpz_t*)th, k, curA[k]);
    }
    slong prec = K->prec + 16;

    /* basis embeddings sigma_k(theta^j), and sigma_k(theta) for rescaling. */
    acb_t bemb[3][4], Sc[3];
    for (int k = 0; k < m; k++) {
        acb_init(Sc[k]); nf_embed_int(K, (const mpz_t*)th, k, Sc[k]);
        for (int j = 0; j < deg; j++) {
            acb_init(bemb[k][j]);
            mpz_t ej[4]; for (int i = 0; i < deg; i++) mpz_init_set_ui(ej[i], i == j ? 1u : 0u);
            nf_embed_int(K, (const mpz_t*)ej, k, bemb[k][j]);
            for (int i = 0; i < deg; i++) mpz_clear(ej[i]);
        }
    }

    /* LLL rows: 4 basis vectors x (m real Minkowski coords, split complex into
     * Re/Im => 4 metric columns) + identity(4).  Column scaling makes the box
     * "sigma0 < U, others small" isotropic: sigma0 column x 2^(P-u) (== /U);
     * the sigma1 column x 2^P when we SHRINK sigma1 (dir=-1) or x 2^(P-u) when
     * we GROW it with sigma0 (dir=+1); the complex Re/Im x 2^P (always shrink). */
    mpq_t* rows = malloc(sizeof(mpq_t) * (size_t)(deg * (4 + deg)));
    for (int i = 0; i < deg * (4 + deg); i++) mpq_init(rows[i]);
    mpz_t T[4][4]; for (int a = 0; a < deg; a++) for (int b = 0; b < deg; b++) mpz_init(T[a][b]);
    mpz_t cand[4], best[4]; for (int i = 0; i < deg; i++) { mpz_init(cand[i]); mpz_init(best[i]); }
    arb_t re, im, scaled, candA[3], bestS0;
    arb_init(re); arb_init(im); arb_init(scaled); arb_init(bestS0);
    for (int k = 0; k < m; k++) arb_init(candA[k]);
    acb_t z; acb_init(z); fmpz_t zi; fmpz_init(zi); mpz_t tmp; mpz_init(tmp);
    bool have_best = false; int found_u = -1;
    const int NCOL = 4 + deg;

    const int UMAX_U = 64;
    for (int u = 1; u <= UMAX_U && !have_best; u++) {
        for (int j = 0; j < deg; j++) {
            /* metric col 0: Re(sigma0(theta^j)/sigma0(theta)) scaled /U (2^(P-u)) */
            acb_div(z, bemb[0][j], Sc[0], prec);
            arb_set(scaled, acb_realref(z)); arb_mul_2exp_si(scaled, scaled, P - u);
            arf_get_fmpz(zi, arb_midref(scaled), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * NCOL + 0], tmp);

            /* metric col 1: sigma1(theta^j)/sigma1(theta); grow with sigma0
             * (dir>0 -> also /U) or shrink (dir<0 -> full 2^P). */
            long e1 = (dir > 0) ? (P - u) : P;
            acb_div(z, bemb[1][j], Sc[1], prec);
            arb_set(scaled, acb_realref(z)); arb_mul_2exp_si(scaled, scaled, e1);
            arf_get_fmpz(zi, arb_midref(scaled), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * NCOL + 1], tmp);

            /* metric cols 2,3: Re/Im of sigma2(theta^j)/sigma2(theta), scaled 2^P */
            acb_div(z, bemb[2][j], Sc[2], prec);
            arb_set(re, acb_realref(z)); arb_set(im, acb_imagref(z));
            arb_mul_2exp_si(re, re, P);
            arf_get_fmpz(zi, arb_midref(re), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * NCOL + 2], tmp);
            arb_mul_2exp_si(im, im, P);
            arf_get_fmpz(zi, arb_midref(im), ARF_RND_NEAR);
            fmpz_get_mpz(tmp, zi); mpq_set_z(rows[j * NCOL + 3], tmp);

            for (int c = 0; c < deg; c++) mpq_set_si(rows[j * NCOL + 4 + c], c == j ? 1 : 0, 1);
        }

        if (lll_reduce_q(rows, deg, NCOL, NULL) != 0) continue;
        bool ok = true;
        for (int a = 0; a < deg && ok; a++)
            for (int c = 0; c < deg; c++) {
                mpq_ptr e = rows[a * NCOL + 4 + c];
                if (mpz_cmp_ui(mpq_denref(e), 1) != 0) { ok = false; break; }
                mpz_set(T[a][c], mpq_numref(e));
            }
        if (!ok) continue;

        for (int E = 1; E <= 3; E++) {
            int W = 2 * E + 1; long total = 1; for (int i = 0; i < deg; i++) total *= W;
            for (long idx = 0; idx < total; idx++) {
                long t = idx; int e[4]; bool nz = false, hi = false;
                for (int i = 0; i < deg; i++) { e[i] = (int)(t % W) - E; t /= W; if (e[i]) nz = true; if (abs(e[i]) == E) hi = true; }
                if (!nz || (E > 1 && !hi)) continue;
                for (int c = 0; c < deg; c++) {
                    mpz_set_ui(cand[c], 0);
                    for (int a = 0; a < deg; a++) {
                        if (e[a] > 0)      mpz_addmul_ui(cand[c], T[a][c], (unsigned long)e[a]);
                        else if (e[a] < 0) mpz_submul_ui(cand[c], T[a][c], (unsigned long)(-e[a]));
                    }
                }
                emb_abs(K, (const mpz_t*)cand, 0, candA[0]);
                if (!arb_lt(curA[0], candA[0])) continue;             /* sigma0 must grow */
                emb_abs(K, (const mpz_t*)cand, 1, candA[1]);
                if (dir > 0) { if (!arb_lt(curA[1], candA[1])) continue; }   /* grow sigma1 */
                else         { if (!arb_lt(candA[1], curA[1])) continue; }   /* shrink sigma1 */
                emb_abs(K, (const mpz_t*)cand, 2, candA[2]);
                if (!arb_lt(candA[2], curA[2])) continue;             /* |sigma2| must shrink */
                if (!have_best || arb_lt(candA[0], bestS0)) {
                    for (int c = 0; c < deg; c++) mpz_set(best[c], cand[c]);
                    arb_set(bestS0, candA[0]); have_best = true;
                }
            }
        }
        if (have_best) found_u = u;
    }

    if (v2dbg()) {
        fprintf(stderr, "    [voro2-step dir=%+d] P=%ld found_u=%d best=%d", dir, P, found_u, (int)have_best);
        if (have_best) gmp_fprintf(stderr, " -> [%Zd,%Zd,%Zd,%Zd]", best[0], best[1], best[2], best[3]);
        fprintf(stderr, "\n");
    }
    if (have_best) for (int c = 0; c < deg; c++) mpz_set(th[c], best[c]);

    for (int i = 0; i < deg * NCOL; i++) mpq_clear(rows[i]);
    free(rows);
    for (int k = 0; k < m; k++) { acb_clear(Sc[k]); for (int j = 0; j < deg; j++) acb_clear(bemb[k][j]); arb_clear(curA[k]); arb_clear(candA[k]); }
    acb_clear(z); arb_clear(re); arb_clear(im); arb_clear(scaled); arb_clear(bestS0);
    fmpz_clear(zi); mpz_clear(tmp);
    for (int i = 0; i < deg; i++) { mpz_clear(cand[i]); mpz_clear(best[i]); }
    for (int a = 0; a < deg; a++) for (int b = 0; b < deg; b++) mpz_clear(T[a][b]);
    return have_best;
}

/* Largest coordinate bit-size the chain tolerates before giving up.  A CORRECT
 * walk reaches the fundamental unit at small coordinates (Q(10^{1/4})'s second
 * unit has coords ~259, ~9 bits); the current first-cut gates do NOT converge
 * (the minima structure is genuinely 2-D -- see the file header), so the coords
 * blow up.  This cap makes a non-converging walk DECLINE fast (contract-safe:
 * declining is never wrong) without blocking a corrected walk later. */
#define VORO2_COORD_BITS 48

static bool coords_too_big(const mpz_t* th, int deg) {
    for (int i = 0; i < deg; i++)
        if (mpz_sizeinbase(th[i], 2) > VORO2_COORD_BITS) return true;
    return false;
}

/* Walk one directed chain from theta=1, returning the FIRST unit (|N|==1,
 * != +-1) reached, in out[0..deg-1] (pre-init'd mpz).  Returns true on find. */
static bool chain_first_unit(struct NumberField* K, int dir, mpz_t* out) {
    const int deg = 4;
    mpz_t th[4]; for (int i = 0; i < deg; i++) mpz_init_set_ui(th[i], i == 0 ? 1u : 0u);
    slong wp = 400;
    mpz_t nrm; mpz_init(nrm);
    bool found = false;
    const int MAXSTEP = 40;
    for (int step = 0; step < MAXSTEP && !found; step++) {
        if (!step_sig21(K, th, dir, &wp)) break;
        if (coords_too_big((const mpz_t*)th, deg)) {          /* non-convergent -> bail fast */
            if (v2dbg()) fprintf(stderr, "  [voro2 dir=%+d] coords exceeded %d bits at step %d -> decline\n", dir, VORO2_COORD_BITS, step);
            break;
        }
        bool isint = nf_norm_int(K, (const mpz_t*)th, nrm);
        if (v2dbg() && (step < 4 || (step % 25) == 0))
            gmp_fprintf(stderr, "  [voro2 dir=%+d] step %d N=%Zd\n", dir, step, nrm);
        /* a non-trivial O_K element whose norm has |N| == ok_den^deg is a unit */
        if (isint) {
            /* is it +-1 (torsion)?  skip; else it's a genuine unit proposal */
            bool nonzero_hi = false;
            for (int i = 1; i < deg; i++) if (mpz_sgn(th[i]) != 0) nonzero_hi = true;
            if (mpz_cmpabs_ui(nrm, 1) == 0 && nonzero_hi) {
                for (int i = 0; i < deg; i++) mpz_set(out[i], th[i]);
                found = true;
                if (v2dbg()) gmp_fprintf(stderr, "  [voro2 dir=%+d] UNIT at step %d: [%Zd,%Zd,%Zd,%Zd]\n",
                                         dir, step, th[0], th[1], th[2], th[3]);
            }
        }
    }
    mpz_clear(nrm);
    for (int i = 0; i < deg; i++) mpz_clear(th[i]);
    return found;
}

/*
 * Propose up to 2 independent units for a signature-(2,1) rank-2 field by two
 * directed minima chains (the (1,1) subfield direction and the (1,-1) relative
 * direction).  Fills out[0], out[1] (each `deg` pre-init'd mpz numerator coords
 * in the theta-power basis) and returns the number proposed (0, 1, or 2).  The
 * caller certifies (p-saturation); a non-fundamental or dependent proposal is
 * rejected there and DECLINEs.
 */
int nf_voronoi_units_sig21(struct NumberField* K, mpz_t** out) {
    if (K->deg != 4 || K->r1 != 2 || K->r2 != 1) return 0;
    const int deg = 4;
    int n = 0;
    /* NOTE: monogenic-only for now (O_K = Z[theta]); the chain tracks theta-power
     * coords directly.  Non-monogenic sig-(2,1) fields decline here (safe). */
    if (nf_ok_index(K) != 1) { if (v2dbg()) fprintf(stderr, "  [voro2] non-monogenic sig(2,1): decline\n"); return 0; }

    mpz_t u2[4];
    for (int i = 0; i < deg; i++) mpz_init(u2[i]);

    /* The (1,-1) "relative" chain converges (it reaches a unit whose log lies
     * off the subfield's (1,1) log-direction).  The subfield unit itself is
     * already in the caller's box-candidate list, so we do NOT run the (1,1)
     * chain here (it does not converge yet -- see the header); the caller
     * combines this unit with the box's subfield unit and certifies the pair. */
    if (chain_first_unit(K, -1, u2)) { for (int i = 0; i < deg; i++) mpz_set(out[0][i], u2[i]); n = 1; }

    for (int i = 0; i < deg; i++) mpz_clear(u2[i]);
    if (v2dbg()) fprintf(stderr, "  [voro2] proposed %d unit(s)\n", n);
    return n;
}

#endif /* USE_FLINT */
