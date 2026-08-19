/*
 * nfunits.c -- see nfunits.h.
 *
 * Fundamental units of a monogenic number field, certified by
 * p-saturation.  All completeness-affecting decisions are exact:
 *   - unit test |N(b)| == 1 is exact (nf_norm_int);
 *   - p-saturation is a rank test of an exact mod-p character matrix,
 *     and rank == r proves saturation unconditionally.
 * The search box and the (double) regulator only PROPOSE / bound; they
 * never decide a unit is fundamental on their own.
 */
#include "nfunits.h"
#include "numberfield_internal.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

struct NFUnits {
    int      rank;
    int      deg;
    mpz_t**  coords;     /* rank vectors of `deg` mpz coordinates */
    double   regulator;
};

int    nf_units_rank(const NFUnits* U)      { return U->rank; }
double nf_units_regulator(const NFUnits* U) { return U->regulator; }
const mpz_t* nf_units_coords(const NFUnits* U, int j) {
    return (const mpz_t*)U->coords[j];
}

void nf_units_free(NFUnits* U) {
    if (!U) return;
    for (int j = 0; j < U->rank; j++) {
        for (int i = 0; i < U->deg; i++) mpz_clear(U->coords[j][i]);
        free(U->coords[j]);
    }
    free(U->coords);
    free(U);
}

#ifndef USE_FLINT
NFUnits* nf_fundamental_units(NumberField* K) { (void)K; return NULL; }
#else

#include <flint/arb.h>
#include <flint/arf.h>

/* ---- small integer helpers (all q < ~10^6, products fit in 64 bits) ---- */

static int64_t mulmod_i64(int64_t a, int64_t b, int64_t m) {
    return (int64_t)(((__int128)a * (__int128)b) % m);
}
static int64_t powmod_i64(int64_t a, int64_t e, int64_t m) {
    int64_t r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) { if (e & 1) r = mulmod_i64(r, a, m); a = mulmod_i64(a, a, m); e >>= 1; }
    return r;
}
static bool is_prime_i64(int64_t n) {
    if (n < 2) return false;
    for (int64_t d = 2; d * d <= n; d++) if (n % d == 0) return false;
    return true;
}

/* Evaluate the integer polynomial with coords c[0..deg-1] at `a` mod q.
 * Coordinates are mpz (a fundamental unit's grow with the regulator and can far
 * exceed int64 -- Q(cbrt41)'s is 24 digits), reduced mod q exactly. */
static int64_t eval_coords_mod(const mpz_t* c, int deg, int64_t a, int64_t q) {
    int64_t v = 0, ap = 1 % q;
    for (int j = 0; j < deg; j++) {
        int64_t cj = (int64_t)mpz_fdiv_ui(c[j], (unsigned long)q);   /* c[j] mod q, >= 0 */
        v = (v + mulmod_i64(cj, ap, q)) % q;
        ap = mulmod_i64(ap, a, q);
    }
    return v;
}

/* ---- collected units ---- */

typedef struct {
    int64_t* c;             /* deg coords (small, box-bounded) */
    double*  logv;          /* r1+r2 log-embedding coords */
    double   nrm2;          /* squared l2 norm of logv */
} Cand;

/* Log-embedding of integer coords into a freshly-malloc'd double[r1+r2]. */
static double* log_embed(NumberField* K, const int64_t* c) {
    int deg = nf_degree(K), m = nf_r1(K) + nf_r2(K);
    mpz_t* mc = malloc(sizeof(mpz_t) * (size_t)deg);
    for (int i = 0; i < deg; i++) mpz_init_set_si(mc[i], (long)c[i]);
    arb_ptr lv = _arb_vec_init(m);
    nf_log_embedding_int(K, (const mpz_t*)mc, lv);
    double* out = malloc(sizeof(double) * (size_t)m);
    for (int i = 0; i < m; i++) out[i] = arf_get_d(arb_midref(lv + i), ARF_RND_NEAR);
    _arb_vec_clear(lv, m);
    for (int i = 0; i < deg; i++) mpz_clear(mc[i]);
    free(mc);
    return out;
}

/* Exact unit test |N(coords)| == 1. */
static bool is_unit(NumberField* K, const int64_t* c) {
    int deg = nf_degree(K);
    mpz_t* mc = malloc(sizeof(mpz_t) * (size_t)deg);
    for (int i = 0; i < deg; i++) mpz_init_set_si(mc[i], (long)c[i]);
    mpz_t nrm; mpz_init(nrm);
    bool ok = nf_norm_int(K, (const mpz_t*)mc, nrm);
    bool unit = ok && (mpz_cmpabs_ui(nrm, 1) == 0);
    mpz_clear(nrm);
    for (int i = 0; i < deg; i++) mpz_clear(mc[i]);
    free(mc);
    return unit;
}

/* First nonzero coordinate positive (canonical rep, drops the u/-u pair). */
static bool canonical_sign(const int64_t* c, int deg) {
    for (int i = 0; i < deg; i++) { if (c[i] > 0) return true; if (c[i] < 0) return false; }
    return false;   /* all zero */
}

/* ---- double linear algebra (independence + regulator) ---- */

/* Try to append row v (length r) to an r-independent set held in `basis`
 * (nb rows, each length r, row-echelon). Returns true and stores if it
 * raises the rank. */
static bool add_if_independent(double* basis, int* nb, const double* v, int r) {
    double* w = malloc(sizeof(double) * (size_t)r);
    memcpy(w, v, sizeof(double) * (size_t)r);
    for (int i = 0; i < *nb; i++) {
        const double* bi = basis + (size_t)i * r;
        /* pivot column of bi */
        int pc = -1;
        for (int k = 0; k < r; k++) if (fabs(bi[k]) > 1e-9) { pc = k; break; }
        if (pc < 0) continue;
        double f = w[pc] / bi[pc];
        for (int k = 0; k < r; k++) w[k] -= f * bi[k];
    }
    double mx = 0; for (int k = 0; k < r; k++) mx = fmax(mx, fabs(w[k]));
    if (mx < 1e-7) { free(w); return false; }   /* dependent */
    memcpy(basis + (size_t)(*nb) * r, w, sizeof(double) * (size_t)r);
    (*nb)++;
    free(w);
    return true;
}

/* |det| of the r x r matrix whose rows are the first r coords of the
 * selected units' log-vectors -- the regulator. */
static double regulator_det(double* const* logs, int r) {
    double* M = malloc(sizeof(double) * (size_t)r * (size_t)r);
    for (int i = 0; i < r; i++)
        for (int k = 0; k < r; k++) M[(size_t)i * r + k] = logs[i][k];
    double det = 1.0;
    for (int col = 0; col < r; col++) {
        int piv = -1; double best = 1e-9;
        for (int row = col; row < r; row++) {
            double a = fabs(M[(size_t)row * r + col]);
            if (a > best) { best = a; piv = row; }
        }
        if (piv < 0) { free(M); return 0.0; }
        if (piv != col) { for (int k = 0; k < r; k++) { double t = M[(size_t)piv*r+k]; M[(size_t)piv*r+k]=M[(size_t)col*r+k]; M[(size_t)col*r+k]=t; } det = -det; }
        det *= M[(size_t)col * r + col];
        for (int row = col + 1; row < r; row++) {
            double f = M[(size_t)row*r+col] / M[(size_t)col*r+col];
            for (int k = col; k < r; k++) M[(size_t)row*r+k] -= f * M[(size_t)col*r+k];
        }
    }
    free(M);
    return fabs(det);
}

/* ---- mod-p rank of a character matrix (rows appended incrementally) ---- */

static int rank_mod_p(const int* rows, int nrows, int ncols, int p) {
    int* M = malloc(sizeof(int) * (size_t)nrows * (size_t)ncols);
    for (int i = 0; i < nrows * ncols; i++) M[i] = ((rows[i] % p) + p) % p;
    int rank = 0;
    for (int col = 0; col < ncols && rank < nrows; col++) {
        int piv = -1;
        for (int row = rank; row < nrows; row++) if (M[(size_t)row*ncols+col] % p != 0) { piv = row; break; }
        if (piv < 0) continue;
        for (int k = 0; k < ncols; k++) { int t = M[(size_t)piv*ncols+k]; M[(size_t)piv*ncols+k]=M[(size_t)rank*ncols+k]; M[(size_t)rank*ncols+k]=t; }
        /* normalise pivot to 1 */
        int inv = 1, a = M[(size_t)rank*ncols+col] % p;
        for (int e = 1; e < p; e++) if ((a * e) % p == 1) { inv = e; break; }
        for (int k = 0; k < ncols; k++) M[(size_t)rank*ncols+k] = (M[(size_t)rank*ncols+k] * inv) % p;
        for (int row = 0; row < nrows; row++) {
            if (row == rank) continue;
            int f = M[(size_t)row*ncols+col] % p;
            if (f) for (int k = 0; k < ncols; k++)
                M[(size_t)row*ncols+k] = ((M[(size_t)row*ncols+k] - f * M[(size_t)rank*ncols+k]) % p + p) % p;
        }
        rank++;
    }
    free(M);
    return rank;
}

/* Certify p-saturation: append p-th-power character rows at degree-1 primes
 * q ≡ 1 (mod p) until the r columns reach full rank r (SATURATED, return
 * true) or the prime budget is exhausted (return false -> DECLINE). */
static bool p_saturate(NumberField* K, mpz_t** ucoords, int r, int p, const mpz_t disc) {
    enum { MAXQ = 200 };
    int* rows = malloc(sizeof(int) * (size_t)MAXQ * (size_t)r);
    int nrows = 0;
    int deg = nf_degree(K);
    bool saturated = false;
    /* A FUNDAMENTAL (saturated) set reaches full rank r within ~r + a few
     * degree-1 primes; a NON-fundamental set plateaus below r forever, so scan
     * every prime up to MAXQ is wasted.  Early-out once the rank has stalled
     * below r for `stall_max` consecutive primes -> return "not saturated"
     * (safe: rank==r is the only path to a positive verdict, and this only
     * shortens a set we could not certify anyway). */
    int last_rank = 0, stall = 0;
    const int stall_max = 12;

    for (int64_t q = 2; q < 200000 && nrows < MAXQ && !saturated; q++) {
        if ((q - 1) % p != 0) continue;
        if (!is_prime_i64(q)) continue;
        if (mpz_fdiv_ui(disc, (unsigned long)q) == 0) continue;   /* q | disc */

        /* find a root a of f mod q (f = the field's defining polynomial) */
        int64_t a = -1;
        {
            int64_t fq[64];
            for (int j = 0; j <= deg; j++) fq[j] = (int64_t)mpz_fdiv_ui(K->coeffs[j], (unsigned long)q);
            for (int64_t x = 0; x < q; x++) {
                int64_t v = 0, ap = 1 % q;
                for (int j = 0; j <= deg; j++) { v = (v + mulmod_i64(fq[j], ap, q)) % q; ap = mulmod_i64(ap, x, q); }
                if (v == 0) { a = x; break; }
            }
        }
        if (a < 0) continue;   /* f has no root mod q -> no degree-1 prime here */

        /* generator h of the order-p subgroup of F_q^* */
        int64_t ex = (q - 1) / p, h = 1;
        for (int64_t b = 2; b < q; b++) { int64_t t = powmod_i64(b, ex, q); if (t != 1) { h = t; break; } }
        if (h == 1) continue;   /* should not happen for q ≡ 1 mod p */

        /* character row */
        int* row = rows + (size_t)nrows * r;
        for (int i = 0; i < r; i++) {
            int64_t val = eval_coords_mod(ucoords[i], deg, a, q);
            int64_t t = powmod_i64(val, ex, q);
            int e = 0; int64_t hp = 1;
            for (int k = 0; k < p; k++) { if (hp == t) { e = k; break; } hp = mulmod_i64(hp, h, q); }
            row[i] = e;
        }
        nrows++;
        int rk = rank_mod_p(rows, nrows, r, p);
        if (rk == r) { saturated = true; break; }
        if (rk > last_rank) { last_rank = rk; stall = 0; }
        else if (++stall >= stall_max) break;           /* rank stalled below r */
    }
    free(rows);
    return saturated;
}

/* Certify a proposed independent unit system `sel` (r vectors) as FUNDAMENTAL
 * by p-saturation at every prime p <= R/R0 (Friedman R0 = 0.2).  Full mod-p
 * rank r proves saturation unconditionally; anything short DECLINEs.  Shared
 * by the box search and the Voronoi fallback. */
static bool cert_saturate(NumberField* K, mpz_t** sel, int r, double R) {
    if (R <= 1e-9) return false;
    mpz_t disc; mpz_init(disc); nf_disc(K, disc);
    int Pmax = (int)floor(R / 0.2) + 1; if (Pmax < 2) Pmax = 1;
    bool cert = true;
    for (int p = 2; p <= Pmax && cert; p++) {
        if (!is_prime_i64(p)) continue;
        if (!p_saturate(K, sel, r, p, disc)) cert = false;
    }
    mpz_clear(disc);
    return cert;
}

NFUnits* nf_fundamental_units(NumberField* K) {
    int deg = nf_degree(K);
    int r   = nf_unit_rank(K);
    int m   = nf_r1(K) + nf_r2(K);

    NFUnits* U = malloc(sizeof(NFUnits));
    U->deg = deg; U->rank = r; U->regulator = 0.0;
    U->coords = (r > 0) ? malloc(sizeof(mpz_t*) * (size_t)r) : NULL;

    if (r == 0) return U;   /* imaginary quadratic etc.: no units to certify */

    if (deg > 32) { free(U->coords); free(U); return NULL; }

    bool udbg = getenv("THUE_DEBUG") != NULL;
    clock_t ut = clock();

    /* --- collect units in a growing coefficient box ---
     * NOTE (benchmark 88): a fundamental unit's integer coordinates grow with
     * the regulator, so this box search only finds units of SMALL-regulator
     * fields.  E.g. Q(cbrt2) unit = 1-theta (coords <=1, found), Q(cbrt6) unit
     * has a coord 6, Q(cbrt15) reaches 30, Q(cbrt41) has 11-digit coords and is
     * unreachable by ANY box.  Fields whose units exceed the box DECLINE at
     * Gate 2 (safe, never wrong).  The proper fix is Voronoi's algorithm /
     * reduction-based unit finding; the degree-adaptive cap below is a cheap
     * catch for the smallest-regulator gaps only.  The search STOPS the moment
     * it holds r independent units (they are the shortest, found at small B),
     * so an easy field costs a handful of norm tests, not the whole box -- a
     * blanket `ncand >= r+8` used to scan all 25^3 = 15625 points of a cubic
     * even when the fundamental unit sat at B=1. */
    const int BOX_MAX = (deg == 3) ? 12 : (deg == 4) ? 6 : 4;
    enum { CAND_CAP = 40000 };
    Cand* cand = malloc(sizeof(Cand) * CAND_CAP);
    int ncand = 0;
    int64_t* c = malloc(sizeof(int64_t) * (size_t)deg);

    /* Grow the box one shell at a time; after each shell, greedily pick the
     * SHORTEST r independent units and try to CERTIFY them by p-saturation.
     * Stop the moment certification succeeds -- a rank-1 cubic whose unit sits
     * at B=1 costs one shell, while a higher-rank field grows only as far as
     * needed for a FUNDAMENTAL (index-1) set: a non-fundamental set found early
     * fails saturation, so we keep growing rather than decline prematurely. */
    int64_t** sel    = malloc(sizeof(int64_t*) * (size_t)r);
    double**  sellog = malloc(sizeof(double*)  * (size_t)r);
    /* mpz mirror of the selected coords: the certifier and U's storage read
     * this, so a Voronoi unit whose coordinates exceed int64 is handled the
     * same as a small box unit. */
    mpz_t**   selz   = malloc(sizeof(mpz_t*) * (size_t)r);
    for (int j = 0; j < r; j++) {
        selz[j] = malloc(sizeof(mpz_t) * (size_t)deg);
        for (int i = 0; i < deg; i++) mpz_init(selz[j][i]);
    }
    bool ok = false;
    double R = 0.0;

    for (int B = 1; B <= BOX_MAX && !ok; B++) {
        int W = 2 * B + 1;
        /* iterate the full box [-B,B]^deg via a mixed-radix counter */
        long total = 1; for (int i = 0; i < deg; i++) total *= W;
        for (long idx = 0; idx < total && ncand < CAND_CAP; idx++) {
            long t = idx; bool has_hi = false;
            for (int i = 0; i < deg; i++) { int d = (int)(t % W) - B; c[i] = d; if (abs(d) == B) has_hi = true; t /= W; }
            if (B > 1 && !has_hi) continue;             /* only new shell for B>1 */
            if (!canonical_sign(c, deg)) continue;      /* skips 0 and -u dups */
            if (!is_unit(K, c)) continue;
            /* dedupe exact coords */
            bool dup = false;
            for (int k = 0; k < ncand && !dup; k++) {
                dup = true; for (int i = 0; i < deg; i++) if (cand[k].c[i] != c[i]) { dup = false; break; }
            }
            if (dup) continue;
            Cand* cd = &cand[ncand++];
            cd->c = malloc(sizeof(int64_t) * (size_t)deg);
            memcpy(cd->c, c, sizeof(int64_t) * (size_t)deg);
            cd->logv = log_embed(K, c);
            cd->nrm2 = 0; for (int i = 0; i < m; i++) cd->nrm2 += cd->logv[i] * cd->logv[i];
        }
        if (ncand < r) continue;                        /* not enough units yet */

        /* sort candidates shortest-first, greedily pick r independent */
        for (int i = 1; i < ncand; i++) { Cand key = cand[i]; int j = i - 1; while (j >= 0 && cand[j].nrm2 > key.nrm2) { cand[j+1] = cand[j]; j--; } cand[j+1] = key; }
        double* basis = malloc(sizeof(double) * (size_t)r * (size_t)r);
        int nb = 0, nsel = 0;
        for (int k = 0; k < ncand && nsel < r; k++)
            if (add_if_independent(basis, &nb, cand[k].logv, r)) { sel[nsel] = cand[k].c; sellog[nsel] = cand[k].logv; nsel++; }
        free(basis);
        if (nsel < r) continue;                         /* independent set not full -> grow */

        /* certify by p-saturation (Friedman R0 = 0.2 bounds the prime list) */
        R = regulator_det(sellog, r);
        for (int j = 0; j < r; j++) for (int i = 0; i < deg; i++) mpz_set_si(selz[j][i], (long)sel[j][i]);
        if (cert_saturate(K, selz, r, R)) ok = true;    /* certified fundamental -> stop */
    }
    free(c);
    if (udbg) fprintf(stderr, "[units-prof] search+certify (%d cand, ok=%d): %.1fms\n", ncand, ok, (clock()-ut)*1000.0/CLOCKS_PER_SEC);

    /* Fallback for large-regulator RANK-1 complex cubics whose fundamental
     * unit exceeds the coefficient box: Voronoi's chain of minima proposes the
     * unit, and the SAME p-saturation test certifies it (contract unchanged).
     * The proposed unit can be huge (Q(cbrt41): 24-digit coords), so it lives in
     * mpz throughout -- no int64 truncation. */
    if (!ok && deg == 3 && nf_r1(K) == 1 && nf_r2(K) == 1 && r == 1) {
        clock_t vt = clock();
        bool proposed = nf_voronoi_unit_cubic11(K, selz[0]);
        double voro_ms = (clock() - vt) * 1000.0 / CLOCKS_PER_SEC;
        clock_t ct = clock();
        if (proposed) {
            /* rank-1 regulator = |log|sigma1(eps)|| from the exact log-embedding */
            arb_ptr lv = _arb_vec_init(m);
            nf_log_embedding_int(K, (const mpz_t*)selz[0], lv);
            R = fabs(arf_get_d(arb_midref(lv + 0), ARF_RND_NEAR));
            _arb_vec_clear(lv, m);
            if (cert_saturate(K, selz, 1, R)) ok = true;
        }
        if (udbg) fprintf(stderr, "[units-prof] voronoi fallback ok=%d R=%.3f  voronoi=%.2fms certify=%.2fms\n",
                          ok, R, voro_ms, (clock() - ct) * 1000.0 / CLOCKS_PER_SEC);
    }

    if (ok) {
        U->regulator = R;
        for (int j = 0; j < r; j++) {
            U->coords[j] = malloc(sizeof(mpz_t) * (size_t)deg);
            for (int i = 0; i < deg; i++) mpz_init_set(U->coords[j][i], selz[j][i]);
        }
    }

    /* cleanup candidates */
    for (int k = 0; k < ncand; k++) { free(cand[k].c); free(cand[k].logv); }
    free(cand); free(sel); free(sellog);
    for (int j = 0; j < r; j++) { for (int i = 0; i < deg; i++) mpz_clear(selz[j][i]); free(selz[j]); }
    free(selz);

    if (!ok) { free(U->coords); free(U); return NULL; }
    return U;
}

#endif /* USE_FLINT */
