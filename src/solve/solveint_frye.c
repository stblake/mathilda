/*
 * solveint_frye.c
 *
 * Part of the Solve[..., Integers] engine; split out of solveint.c.
 * See solveint_internal.h for the shared SICtx/SearchState substrate.
 */
#include "solveint.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"
#include "checked_int.h"
#include "poly/mpoly.h"
#include "numbertheory/numbertheory_internal.h"
#include "linalg/hnf.h"
#include "solvethue.h"
#include "solveint_internal.h"


/* ================================================================== *
 *  Ordering-aware, 128-bit meet-in-the-middle for a separable additive
 *  equation over a fully-bounded box.
 *
 *  Generalises the plain mitm_solve in two ways that unlock power-sum
 *  searches the int64 / product-estimate rejects:
 *    - it splits the variables along their ordering chain and enumerates
 *      only INTRA-GROUP-ORDERED tuples (combinations, not the Cartesian
 *      product), so an ordered box  0<x<y<z<w<r<1000  costs
 *      ~C(1000,2)+C(1000,3) instead of 1000^5 -- this is exactly what
 *      the Lander-Parkin quintic  x^5+y^5+z^5+w^5 == r^5  needs;
 *    - partial sums are __int128, so k-th powers past 2^63 do not force a
 *      decline.
 *  Every hit is re-verified against the original conjunction by si_verify
 *  before it is emitted, so cross-group orderings / disequations need no
 *  special handling here (a group enumerates a SUPERSET of its projection
 *  of the true solution set, never a subset, so nothing is missed).
 *
 *  Scope / contract.  Engages only when it adds capability the plain
 *  mitm_solve lacks (values exceed int64, or an ordering chain makes an
 *  otherwise-too-big iterate side fit); otherwise it declines so the
 *  existing path is unchanged.  It only takes boxes it can EXHAUST within
 *  budget, so it returns the complete set (or a proven {}); a box too big
 *  for the budget declines (false), never an unproven answer.
 * ================================================================== */

#define SI_PSM_HASH_CAP   5000000LL      /* max stored hash-group tuples   */
#define SI_PSM_MAXNODES   1200000000LL   /* max iterate-group tuples        */
#define SI_PSM_TABSUM     50000000LL     /* max total value-table entries   */

typedef struct { __int128 sum; int64_t vals[SI_MAX_VARS]; } MitmEntry128;
static int mitm128_cmp(const void* pa, const void* pb) {
    __int128 a = ((const MitmEntry128*)pa)->sum, b = ((const MitmEntry128*)pb)->sum;
    return a < b ? -1 : (a > b ? 1 : 0);
}


/* mpz -> __int128 (LP64: mpz_get_ui yields a 64-bit limb).  Declines any value
 * wider than `bits` so a group's partial sums cannot overflow int128. */
static bool si_mpz_to_i128(const mpz_t z, unsigned bits, __int128* out) {
    if (mpz_sizeinbase(z, 2) > bits) return false;
    mpz_t t, lo, hi; mpz_init_set(t, z);
    bool neg = mpz_sgn(t) < 0; if (neg) mpz_neg(t, t);
    mpz_init(lo); mpz_init(hi);
    mpz_tdiv_r_2exp(lo, t, 64);
    mpz_tdiv_q_2exp(hi, t, 64);
    unsigned long lo64 = mpz_get_ui(lo);
    unsigned long hi64 = mpz_get_ui(hi);
    unsigned __int128 u = ((unsigned __int128)hi64 << 64) | (unsigned __int128)lo64;
    __int128 r = (__int128)u;
    if (neg) r = -r;
    *out = r;
    mpz_clears(t, lo, hi, NULL);
    return true;
}


/* [elo,ehi] for var vi tightened by orderings whose OTHER endpoint is already
 * assigned (intra-group during a group's own enumeration; cross-group endpoints
 * are unset, so they do not fire and are left to si_verify). */
static void mitm128_bounds(const SICtx* c, int vi, const int64_t* vals,
                           const bool* set, int64_t* elo, int64_t* ehi) {
    *elo = c->lo[vi]; *ehi = c->hi[vi];
    for (int k = 0; k < c->n_ord; k++) {
        int a = c->ord_a[k], b = c->ord_b[k], s = c->ord_strict[k] ? 1 : 0;
        if (a == vi && set[b] && vals[b] - s < *ehi) *ehi = vals[b] - s;
        if (b == vi && set[a] && vals[a] + s > *elo) *elo = vals[a] + s;
    }
    for (int k = 0; k < c->n_abs_ord; k++) {
        int a = c->abs_ord_a[k], b = c->abs_ord_b[k], s = c->abs_ord_strict[k] ? 1 : 0;
        if (a == vi && set[b]) {
            int64_t avb = vals[b] < 0 ? -vals[b] : vals[b];
            int64_t bnd = avb - s; if (bnd < 0) bnd = -1;
            if ( bnd < *ehi) *ehi =  bnd;
            if (-bnd > *elo) *elo = -bnd;
        }
    }
}


typedef struct {
    SICtx*   c;
    __int128* gtab[SI_MAX_VARS];      /* value table per var (index v-lo[i]) */
    int      gv[SI_MAX_VARS];         /* current group's vars, topo order    */
    int      ng;
    int64_t  vals[SI_MAX_VARS];       /* working assignment                  */
    bool     set[SI_MAX_VARS];
    MitmEntry128* H; size_t hcnt, hcap; bool hoverflow;   /* collect mode */
    __int128* keys;                                       /* sorted sums, cache-hot */
    __int128 target;                                      /* probe mode   */
    /* Modular residue sieve: skip the (cache-missing) binary search whenever the
     * complement hash-sum residue is unreachable.  M<=1 disables it. */
    int64_t  M, tmod;
    const unsigned char* HR;                              /* achievable hash residues */
    const int* gmod[SI_MAX_VARS];                         /* g_i(v) mod M per var    */
    SearchState* st;
    int64_t  nodes, maxnodes; bool overflow;
} PsmRun;

/* Image ratio |{v^k mod q}| / q, minimised over the exponents present: a proxy
 * for how strongly q prunes (smaller = stronger).  q <= 127. */
static double psm_strength(int q, const int* exps, int nexp) {
    double best = 1.0;
    for (int e = 0; e < nexp; e++) {
        int k = exps[e];
        bool seen[128]; for (int s = 0; s < 128; s++) seen[s] = false;
        int img = 0;
        for (int v = 0; v < q; v++) {
            int p = 1 % q; for (int t = 0; t < k; t++) p = (p * v) % q;
            if (!seen[p]) { seen[p] = true; img++; }
        }
        double r = (double)img / (double)q;
        if (r < best) best = r;
    }
    return best;
}


/* Choose a sieve modulus for a separable power-sum: greedily multiply in the
 * strongest-pruning coprime prime powers (smallest image ratio first) whose
 * product keeps the residue bitset cache-resident.  Returns 1 when no small
 * modulus meaningfully restricts the powers present. */
static int64_t psm_pick_modulus(const int* exps, int nexp) {
    static const int cand[] = {5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,
                               49,61,64,71,81,101,121,125,127};
    int nc = (int)(sizeof(cand) / sizeof(cand[0]));
    double str[64]; int idx[64];
    for (int i = 0; i < nc; i++) { str[i] = psm_strength(cand[i], exps, nexp); idx[i] = i; }
    for (int i = 0; i < nc; i++)                 /* selection sort by strength asc */
        for (int j = i + 1; j < nc; j++)
            if (str[idx[j]] < str[idx[i]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    const int64_t CAP = 1 << 22;                 /* 4M-entry byte bitset = 4 MB */
    int64_t M = 1;
    for (int r = 0; r < nc; r++) {
        int q = cand[idx[r]];
        if (str[idx[r]] > 0.55) break;           /* remaining candidates too weak */
        int64_t a = M, b = q; while (b) { int64_t t = a % b; a = b; b = t; }
        if (a != 1) continue;                    /* keep the prime powers coprime */
        if (M > CAP / q) continue;
        M *= q;
    }
    return M;
}


/* Collect ordered tuples of the hash group into R->H (partial sum accumulated).
 * The table grows on demand (the ordered-count estimate is only approximate),
 * declining once it would exceed the hard hash cap. */
static void psm_collect(PsmRun* R, int gi, __int128 partial) {
    if (R->hoverflow) return;
    if (gi == R->ng) {
        if (R->hcnt >= R->hcap) {
            size_t ncap = R->hcap ? R->hcap * 2 : 1024;
            if (ncap > (size_t)SI_PSM_HASH_CAP) {
                if (R->hcap >= (size_t)SI_PSM_HASH_CAP) { R->hoverflow = true; return; }
                ncap = (size_t)SI_PSM_HASH_CAP;
            }
            MitmEntry128* nh = (MitmEntry128*)realloc(R->H, sizeof(MitmEntry128) * ncap);
            if (!nh) { R->hoverflow = true; return; }
            R->H = nh; R->hcap = ncap;
        }
        MitmEntry128* e = &R->H[R->hcnt++];
        e->sum = partial;
        for (int i = 0; i < R->c->n; i++) e->vals[i] = R->set[i] ? R->vals[i] : 0;
        return;
    }
    int vi = R->gv[gi];
    int64_t elo, ehi; mitm128_bounds(R->c, vi, R->vals, R->set, &elo, &ehi);
    const __int128* tab = R->gtab[vi]; int64_t base = R->c->lo[vi];
    for (int64_t v = elo; v <= ehi && !R->hoverflow; v++) {
        R->vals[vi] = v; R->set[vi] = true;
        psm_collect(R, gi + 1, partial + tab[v - base]);
    }
    R->set[vi] = false;
}


/* Enumerate ordered tuples of the iterate group; binary-search the hash table
 * for the complementary partial sum; verify + emit each full assignment.
 * `pmod` carries the iterate partial sum mod R->M for the residue sieve. */
static void psm_probe(PsmRun* R, int gi, __int128 partial, int64_t pmod) {
    if (R->overflow) return;
    if (gi == R->ng) {
        if (++R->nodes > R->maxnodes) { R->overflow = true; return; }
        if (R->M > 1) {
            int64_t needm = (R->tmod - pmod) % R->M; if (needm < 0) needm += R->M;
            if (!R->HR[needm]) return;              /* complement residue unreachable */
        }
        __int128 need = R->target - partial;
        const __int128* K = R->keys;
        size_t lo = 0, hi = R->hcnt;
        while (lo < hi) { size_t mid = lo + (hi - lo) / 2; if (K[mid] < need) lo = mid + 1; else hi = mid; }
        for (size_t j = lo; j < R->hcnt && K[j] == need; j++) {
            int64_t full[SI_MAX_VARS];
            for (int i = 0; i < R->c->n; i++)
                full[i] = R->set[i] ? R->vals[i] : R->H[j].vals[i];
            if (si_verify(R->c, full)) emit_full(R->st, full);
        }
        return;
    }
    int vi = R->gv[gi];
    int64_t elo, ehi; mitm128_bounds(R->c, vi, R->vals, R->set, &elo, &ehi);
    const __int128* tab = R->gtab[vi]; int64_t base = R->c->lo[vi];
    const int* mtab = R->M > 1 ? R->gmod[vi] : NULL;

    /* Innermost iterate variable: a tight inline loop (no recursion, no
     * mitm128_bounds per node) over the ~10^8 hot leaves.  The residue reduction
     * is hoisted -- one division for the whole loop, then a subtract per node. */
    if (gi == R->ng - 1) {
        R->set[vi] = true;
        const __int128* K = R->keys; size_t hc = R->hcnt;
        const int64_t M = R->M, tm = R->tmod; const unsigned char* HR = R->HR;
        int64_t base_needm = 0;
        if (M > 1) { base_needm = (tm - pmod % M) % M; if (base_needm < 0) base_needm += M; }
        int64_t nodes = R->nodes, maxn = R->maxnodes;
        for (int64_t v = elo; v <= ehi; v++) {
            if (++nodes > maxn) { R->overflow = true; break; }
            if (M > 1) {
                int64_t needm = base_needm - mtab[v - base]; if (needm < 0) needm += M;
                if (!HR[needm]) continue;
            }
            __int128 need = R->target - (partial + tab[v - base]);
            size_t lo = 0, hi = hc;
            while (lo < hi) { size_t mid = lo + (hi - lo) / 2; if (K[mid] < need) lo = mid + 1; else hi = mid; }
            if (lo < hc && K[lo] == need) {
                R->vals[vi] = v;
                for (size_t j = lo; j < hc && K[j] == need; j++) {
                    int64_t full[SI_MAX_VARS];
                    for (int i = 0; i < R->c->n; i++)
                        full[i] = R->set[i] ? R->vals[i] : R->H[j].vals[i];
                    if (si_verify(R->c, full)) emit_full(R->st, full);
                }
            }
        }
        R->nodes = nodes;
        R->set[vi] = false;
        return;
    }

    for (int64_t v = elo; v <= ehi && !R->overflow; v++) {
        R->vals[vi] = v; R->set[vi] = true;
        int64_t nm = mtab ? pmod + mtab[v - base] : 0;
        psm_probe(R, gi + 1, partial + tab[v - base], nm);
    }
    R->set[vi] = false;
}


/* Ordered-count estimate for a set of vars: product of domains divided by the
 * longest intra-set ordering chain's factorial (exact for a full sub-chain of a
 * total order, the case that matters; a harmless over-estimate otherwise -- the
 * node budget backstops any under-estimate, so nothing is ever truncated). */
static long double psm_ordered_est(const SICtx* c, const int* vars, int nv,
                                   const int64_t* domain) {
    long double p = 1.0L;
    for (int i = 0; i < nv; i++) p *= (long double)domain[vars[i]];
    int L = si_longest_chain(c, vars, nv);
    for (int k = 2; k <= L; k++) p /= (long double)k;
    return p;
}


/* ================================================================== *
 *  Frye's biquadrate search:  x^4 + y^4 + z^4 == w^4.
 *
 *  Reference: R. E. Frye, "Finding 95800^4 + 217519^4 + 414560^4 =
 *  422481^4 on the Connection Machine" (1988) -- the minimal counter-
 *  example to Euler's conjecture, w = 422481.  A box up to 10^6 cannot be
 *  exhaustively verified interactively, so this is a WITNESS search: it
 *  finds and returns the (minimal) solution, ascending in w, and declines
 *  (never a spurious {}) if the node budget is spent with nothing found.
 *
 *  Number theory (Frye Sec.3).  For a primitive solution the fourth powers
 *  mod 5 (each 0 or 1) force exactly one summand C to be != 0 mod 5, the
 *  other two (A,B) to be multiples of 5, and the target D=w to be != 0
 *  mod 5.  Then A^4 + B^4 = D^4 - C^4 is a multiple of 5^4 = 625, so
 *  625 | (w^4 - C^4).  Writing A = 5a, B = 5b gives N = (w^4 - C^4)/625 =
 *  a^4 + b^4, decomposed by scanning a over [~0.841 N^{1/4}, N^{1/4}] and
 *  testing whether N - a^4 is a fourth power.  Extra moduli m coprime to 5
 *  prune the (w,C) pairs before the expensive decompose: w^4 - C^4 mod m
 *  must be a sum of two fourth powers mod m (a necessary condition, so no
 *  real primitive solution is ever dropped).  __int128 throughout
 *  (w^4 <= 10^24 < 1.7e38).
 * ================================================================== */

/* floor(T^{1/4}); sets *root and returns true iff T is a perfect fourth power. */
static bool si_is_4th_power_i128(__int128 T, int64_t* root) {
    if (T < 0) return false;
    if (T == 0) { *root = 0; return true; }
    int64_t s = si_isqrt_i128(T);            /* floor(sqrt T) */
    int64_t r = si_isqrt_i64(s);             /* floor(sqrt(floor(sqrt T))) = floor(T^1/4) */
    for (int64_t rr = (r > 1 ? r - 1 : 0); rr <= r + 1; rr++) {
        __int128 rr2 = (__int128)rr * rr;
        if (rr2 * rr2 == T) { *root = rr; return true; }
    }
    return false;
}


#define SI_FRYE_MODS 6
bool si_solve_biquadrate_frye(SICtx* c, SearchState* st) {
    bool DBG = getenv("FRYE_DEBUG") != NULL;
    if (c->n != 4 || c->neq != 1) return false;
    for (int i = 0; i < 4; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;
    const MPoly* eq = c->eq[0];

    /* Shape: each variable once as v^4 with coefficient +/-1, no constant. */
    int sgn[4] = {0,0,0,0}; bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * 4;
        int nz = -1, cnt = 0;
        for (int v = 0; v < 4; v++) if (ex[v] > 0) { nz = v; cnt++; }
        if (cnt == 0) { ok = false; break; }                     /* no constant term */
        if (cnt != 1 || ex[nz] != 4 || sgn[nz] != 0) { ok = false; break; }
        if (mpz_cmp_ui(eq->coefs[t], 1) == 0) sgn[nz] = 1;
        else if (mpz_cmp_si(eq->coefs[t], -1) == 0) sgn[nz] = -1;
        else { ok = false; break; }
    }
    for (int i = 0; i < 4 && ok; i++) if (sgn[i] == 0) ok = false;
    if (!ok) return false;
    /* target = the odd-one-out sign; the other three are the summands. */
    int nplus = 0; for (int i = 0; i < 4; i++) if (sgn[i] > 0) nplus++;
    int wv = -1;
    if (nplus == 1)      { for (int i = 0; i < 4; i++) if (sgn[i] > 0) wv = i; }
    else if (nplus == 3) { for (int i = 0; i < 4; i++) if (sgn[i] < 0) wv = i; }
    else return false;
    int sm[3], ns = 0; for (int i = 0; i < 4; i++) if (i != wv) sm[ns++] = i;

    /* Positive box, and w the largest.  Engage only on a box too big for the
     * exhaustive separable MITM (small boxes are handled completely there); a
     * tiny box declines so this witness search never hides a real {} proof. */
    int64_t whi = c->hi[wv];
    if (whi < 20000) return false;                    /* leave to the exhaustive path */
    for (int i = 0; i < 4; i++) if (c->lo[i] < 0) return false;   /* positive summands */
    int64_t wlo = c->lo[wv] > 1 ? c->lo[wv] : 1;

    /* Residue tables. */
    int q4_625[625]; for (int r = 0; r < 625; r++) { long long v = ((long long)r*r % 625); v = v*v % 625; q4_625[r] = (int)v; }
    /* Fast 4th-power reject: the low 16 bits of any fourth power lie in a sparse
     * set, so `Tb & 0xFFFF` filters ~98% of non-fourth-powers before the exact
     * (costly) integer 4th-root test. */
    static unsigned char is4_16bit[65536]; static bool is4_built = false;
    if (!is4_built) { for (int b = 0; b < 65536; b++) is4_16bit[(int)(((unsigned)b*b*b*b) & 0xFFFF)] = 1; is4_built = true; }
    /* Valid C residues per target residue: units c (c%5!=0) with c^4 == w^4 mod 625.
     * Secondary moduli: primes with 4 | (p-1) (fourth powers strongly clustered)
     * plus 16 and 9. */
    static const int mods[SI_FRYE_MODS] = {16, 9, 13, 17, 29, 41};
    int q4m[SI_FRYE_MODS][64];
    unsigned char s2[SI_FRYE_MODS][64];               /* sum of two 4th powers mod m */
    for (int mi = 0; mi < SI_FRYE_MODS; mi++) {
        int m = mods[mi];
        for (int r = 0; r < m; r++) { long long v = ((long long)r*r % m); v = v*v % m; q4m[mi][r] = (int)v; }
        for (int r = 0; r < m; r++) s2[mi][r] = 0;
        for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) s2[mi][(q4m[mi][a] + q4m[mi][b]) % m] = 1;
    }

    int topo_rank[4];                                 /* for assigning sorted summands */
    { for (int i = 0; i < 4; i++) { int rk = 0; for (int k = 0; k < c->n_ord; k++) if (c->ord_b[k] == i) rk++; topo_rank[i] = rk; } }

    int64_t maxnodes = 300000000000LL;                /* ~ reaches w=422481 */
    { const char* e = getenv("MATHILDA_FRYE_MAXNODES"); if (e) maxnodes = strtoll(e, NULL, 10); }
    int64_t nodes = 0;
    bool found = false;
    int64_t sol[4];

    for (int64_t w = wlo; w <= whi && !found; w++) {
        int wq = (int)(w % 5); if (wq == 0) continue;              /* w != 0 mod 5 */
        int w625 = (int)(w % 625); int w4r = q4_625[w625];
        int w4m[SI_FRYE_MODS]; for (int mi = 0; mi < SI_FRYE_MODS; mi++) w4m[mi] = q4m[mi][(int)(w % mods[mi])];
        __int128 w4 = (__int128)w*w; w4 = w4*w4;
        /* Enumerate C != 0 mod 5 with 625 | w^4 - C^4, C in [1, w-1]. */
        for (int cr = 1; cr < 625; cr++) {
            if (cr % 5 == 0) continue;
            if (q4_625[cr] != w4r) continue;                      /* C^4 == w^4 mod 625 */
            for (int64_t C = cr; C < w; C += 625) {
                /* secondary modular sieve: w^4 - C^4 mod m must be a sum of two 4th powers */
                bool pass = true;
                for (int mi = 0; mi < SI_FRYE_MODS; mi++) {
                    int d = w4m[mi] - q4m[mi][(int)(C % mods[mi])]; if (d < 0) d += mods[mi];
                    if (!s2[mi][d]) { pass = false; break; }
                }
                if (!pass) continue;
                __int128 C4 = (__int128)C*C; C4 = C4*C4;
                __int128 diff = w4 - C4;                          /* = 625 N */
                if (diff <= 0) continue;
                if (diff % 625 != 0) continue;                    /* guard */
                __int128 N = diff / 625;                          /* = a^4 + b^4 */
                /* Prime-factor constraint (Frye Sec.3.3): in N = a^4 + b^4 every odd
                 * prime P != 1 (mod 8) divides to an exponent = 0 (mod 4) (it must
                 * divide both a and b); a small-prime trial division rejecting any
                 * P whose exponent is not a multiple of 4 eliminates most pairs
                 * before the decompose, and never drops a real primitive solution. */
                {
                    static const int badp[] = {3,7,11,13,19,23,29,31,37,43,47,53,59,61,67,71,79,83,
                        101,103,107,109,127,131,139,149,151,157,163,167,173,179,181,191,197,199,211,223,227,229,239,251};
                    __int128 t = N; bool reject = false;
                    for (size_t pi = 0; pi < sizeof(badp)/sizeof(badp[0]); pi++) {
                        int P = badp[pi];
                        if (t % P != 0) continue;
                        int e = 0; while (t % P == 0) { t /= P; e++; }
                        if (e % 4 != 0) { reject = true; break; }
                    }
                    if (reject) continue;
                }
                /* decompose: scan a in [ceil(0.8409 N^1/4), N^1/4]; b^4 = N - a^4 */
                int64_t amax; { int64_t s = si_isqrt_i128(N); amax = si_isqrt_i64(s); }
                int64_t amin = (int64_t)(0.84089641525L * (long double)amax);
                if (amin < 1) amin = 1;
                for (int64_t a = amax; a >= amin; a--) {
                    if (++nodes > maxnodes) { w = whi + 1; found = false; goto frye_done; }
                    __int128 a4 = (__int128)a*a; a4 = a4*a4;
                    __int128 Tb = N - a4;
                    if (Tb < 0) continue;
                    if (!is4_16bit[(unsigned)Tb & 0xFFFF]) continue;   /* fast reject */
                    int64_t b;
                    if (si_is_4th_power_i128(Tb, &b) && b <= a) {
                        int64_t A = 5*a, B = 5*b, Cv = C;         /* summand values */
                        int64_t s3[3] = {A, B, Cv};
                        /* sort ascending */
                        for (int i2 = 0; i2 < 3; i2++) for (int j2 = i2+1; j2 < 3; j2++) if (s3[j2] < s3[i2]) { int64_t t = s3[i2]; s3[i2] = s3[j2]; s3[j2] = t; }
                        /* assign to summand positions by topo rank (smallest value -> smallest-rank var) */
                        int pos[3] = {sm[0], sm[1], sm[2]};
                        for (int i2 = 0; i2 < 3; i2++) for (int j2 = i2+1; j2 < 3; j2++) if (topo_rank[pos[j2]] < topo_rank[pos[i2]]) { int t = pos[i2]; pos[i2] = pos[j2]; pos[j2] = t; }
                        int64_t vals[SI_MAX_VARS]; for (int i2 = 0; i2 < 4; i2++) vals[i2] = 0;
                        vals[wv] = w; for (int i2 = 0; i2 < 3; i2++) vals[pos[i2]] = s3[i2];
                        if (si_verify(c, vals)) {
                            for (int i2 = 0; i2 < 4; i2++) sol[i2] = vals[i2];
                            found = true; goto frye_done;
                        }
                    }
                }
            }
        }
    }
frye_done:
    if (DBG) fprintf(stderr, "[frye] nodes=%lld found=%d\n", (long long)nodes, (int)found);
    if (!found) return false;                         /* budget spent / none: decline */
    emit_full(st, sol);
    return true;
}

bool si_solve_separable_mitm(SICtx* c, SearchState* st) {
    if (c->neq != 1) return false;
    int n = c->n;
    if (n < 2 || n > SI_MAX_VARS) return false;
    const MPoly* eq = c->eq[0];

    /* Separability + constant term, full boundedness, and the set of exponents
     * present (for the modular sieve). */
    mpz_t c0; mpz_init(c0); mpz_set_ui(c0, 0);
    int exps_present[SI_MAX_VARS * 4], nexp = 0;
    for (size_t t = 0; t < eq->n_terms; t++) {
        const int* ex = eq->exps + t * (size_t)n;
        int nv = 0, kexp = 0;
        for (int v = 0; v < n; v++) if (ex[v] > 0) { nv++; kexp = ex[v]; }
        if (nv >= 2) { mpz_clear(c0); return false; }        /* not separable */
        if (nv == 0) { mpz_add(c0, c0, eq->coefs[t]); continue; }
        bool have = false;
        for (int e = 0; e < nexp; e++) if (exps_present[e] == kexp) { have = true; break; }
        if (!have && nexp < (int)(sizeof(exps_present)/sizeof(exps_present[0]))) exps_present[nexp++] = kexp;
    }
    int64_t domain[SI_MAX_VARS];
    long double tabsum = 0.0L;
    for (int i = 0; i < n; i++) {
        if (!(c->has_lo[i] && c->has_hi[i]) || c->hi[i] < c->lo[i]) { mpz_clear(c0); return false; }
        domain[i] = c->hi[i] - c->lo[i] + 1;
        tabsum += (long double)domain[i];
    }
    if (tabsum > (long double)SI_PSM_TABSUM) { mpz_clear(c0); return false; }

    /* Topological order of the variables consistent with the < / <= chain, so a
     * contiguous prefix/suffix split keeps ordering-adjacent variables together
     * (the only split that turns a chain into cheap combinations). */
    int topo[SI_MAX_VARS], nt = 0; bool placed[SI_MAX_VARS];
    for (int i = 0; i < n; i++) placed[i] = false;
    for (int iter = 0; iter < n; iter++) {
        int pick = -1;
        for (int v = 0; v < n; v++) {
            if (placed[v]) continue;
            bool ready = true;                    /* all a with a<v already placed */
            for (int k = 0; k < c->n_ord; k++)
                if (c->ord_b[k] == v && !placed[c->ord_a[k]]) { ready = false; break; }
            if (ready) { pick = v; break; }
        }
        if (pick < 0) { for (int v = 0; v < n; v++) if (!placed[v]) { pick = v; break; } } /* cycle guard */
        placed[pick] = true; topo[nt++] = pick;
    }

    /* Choose the contiguous split minimising the (larger) iterate est, subject
     * to: smaller group <= hash cap, larger group <= node budget.  Also compute
     * the best UNORDERED iterate est to decide whether we add value over the
     * plain mitm_solve. */
    int best_p = -1; long double best_iter = 0.0L;
    long double best_unord_iter = 1e300L;
    for (int p = 1; p < n; p++) {
        int A[SI_MAX_VARS], B[SI_MAX_VARS], na = 0, nb = 0;
        for (int i = 0; i < p; i++) A[na++] = topo[i];
        for (int i = p; i < n; i++) B[nb++] = topo[i];
        long double eA = psm_ordered_est(c, A, na, domain);
        long double eB = psm_ordered_est(c, B, nb, domain);
        long double small = eA < eB ? eA : eB, large = eA < eB ? eB : eA;
        /* unordered products (no chain discount) for the value-add test. */
        long double uA = 1.0L, uB = 1.0L;
        for (int i = 0; i < na; i++) uA *= (long double)domain[A[i]];
        for (int i = 0; i < nb; i++) uB *= (long double)domain[B[i]];
        long double ularge = uA < uB ? uB : uA;
        if (ularge < best_unord_iter) best_unord_iter = ularge;
        if (small > (long double)SI_PSM_HASH_CAP) continue;
        if (large > (long double)SI_PSM_MAXNODES) continue;
        if (best_p < 0 || large < best_iter) { best_p = p; best_iter = large; }
    }
    if (best_p < 0) { mpz_clear(c0); return false; }         /* no feasible split */

    /* Build value tables (exact via mpz, then narrowed to int128).  Track whether
     * any value would overflow int64 -- part of the "adds value" test. */
    __int128* gtab[SI_MAX_VARS];
    for (int i = 0; i < n; i++) gtab[i] = NULL;
    bool ok = true, exceeds_i64 = false;
    mpz_t gv, pw; mpz_init(gv); mpz_init(pw);
    for (int i = 0; i < n && ok; i++) {
        gtab[i] = (__int128*)malloc(sizeof(__int128) * (size_t)domain[i]);
        if (!gtab[i]) { ok = false; break; }
        for (int64_t d = 0; d < domain[i] && ok; d++) {
            int64_t v = c->lo[i] + d;
            mpz_set_ui(gv, 0);
            for (size_t t = 0; t < eq->n_terms; t++) {
                const int* ex = eq->exps + t * (size_t)n;
                if (ex[i] <= 0) continue;
                bool only_i = true;
                for (int u = 0; u < n; u++) if (u != i && ex[u] > 0) { only_i = false; break; }
                if (!only_i) continue;
                mpz_set(pw, eq->coefs[t]);
                for (int k = 0; k < ex[i]; k++) mpz_mul_si(pw, pw, (long)v);
                mpz_add(gv, gv, pw);
            }
            if (!mpz_fits_slong_p(gv)) exceeds_i64 = true;
            __int128 g; if (!si_mpz_to_i128(gv, 115, &g)) { ok = false; break; }
            gtab[i][d] = g;
        }
    }
    __int128 target = 0;
    if (ok) { mpz_neg(c0, c0); if (!si_mpz_to_i128(c0, 120, &target)) ok = false; }
    mpz_clear(gv); mpz_clear(pw); mpz_clear(c0);
    if (!ok) { for (int i = 0; i < n; i++) free(gtab[i]); return false; }

    /* Adds-value gate: engage only if the plain int64 mitm_solve would decline
     * this box -- i.e. values overflow int64, or its unordered iterate side
     * exceeds SI_MAX_NODES while ours (ordered) fits.  Otherwise defer. */
    if (!exceeds_i64 && best_unord_iter <= (long double)SI_MAX_NODES) {
        for (int i = 0; i < n; i++) free(gtab[i]);
        return false;
    }

    /* Split into hash (smaller) and iterate (larger) contiguous groups. */
    int A[SI_MAX_VARS], B[SI_MAX_VARS], na = 0, nb = 0;
    for (int i = 0; i < best_p; i++) A[na++] = topo[i];
    for (int i = best_p; i < n; i++) B[nb++] = topo[i];
    long double eA = psm_ordered_est(c, A, na, domain), eB = psm_ordered_est(c, B, nb, domain);
    int *hg, nh, *ig, ni;
    if (eA <= eB) { hg = A; nh = na; ig = B; ni = nb; } else { hg = B; nh = nb; ig = A; ni = na; }

    size_t hcap = (size_t)((eA <= eB ? eA : eB)) + 64;
    MitmEntry128* H = (MitmEntry128*)malloc(sizeof(MitmEntry128) * hcap);
    if (!H) { for (int i = 0; i < n; i++) free(gtab[i]); return false; }

    PsmRun R;
    R.c = c; R.st = st; R.H = H; R.hcnt = 0; R.hcap = hcap; R.hoverflow = false;
    R.keys = NULL; R.M = 1; R.tmod = 0; R.HR = NULL;
    R.target = target; R.nodes = 0; R.maxnodes = SI_PSM_MAXNODES; R.overflow = false;
    for (int i = 0; i < n; i++) { R.gtab[i] = gtab[i]; R.gmod[i] = NULL; R.vals[i] = 0; R.set[i] = false; }

    /* Collect the hash group. */
    for (int i = 0; i < nh; i++) R.gv[i] = hg[i];
    R.ng = nh;
    psm_collect(&R, 0, 0);
    if (R.hoverflow) {                                       /* hash group exceeds cap */
        free(R.H); for (int i = 0; i < n; i++) free(gtab[i]); st->nsol = 0; return false;
    }
    qsort(R.H, R.hcnt, sizeof(MitmEntry128), mitm128_cmp);

    /* Compact sorted key array (just the sums) so the ~10^8 probe binary
     * searches stride 16 B, not the 96 B full entries -- keeps the working set
     * cache-resident. */
    R.keys = (__int128*)malloc(sizeof(__int128) * (R.hcnt ? R.hcnt : 1));
    if (!R.keys) { free(R.H); for (int i = 0; i < n; i++) free(gtab[i]); return false; }
    for (size_t j = 0; j < R.hcnt; j++) R.keys[j] = R.H[j].sum;

    /* Modular residue sieve: the achievable residues of the hash-group sum mod M.
     * Built exactly from the collected sums, so it is a sound necessary condition
     * (never prunes a real solution); it pays off only when the power residues
     * cluster, so it is disabled if fewer than ~40% of classes are pruned. */
    int* gmod[SI_MAX_VARS]; for (int i = 0; i < n; i++) gmod[i] = NULL;
    unsigned char* HR = NULL;
    int64_t M = psm_pick_modulus(exps_present, nexp);
    if (M > 1) {
        HR = (unsigned char*)calloc((size_t)M, 1);
        if (HR) {
            for (size_t j = 0; j < R.hcnt; j++) {
                __int128 s = R.H[j].sum % M; int64_t m = (int64_t)s; if (m < 0) m += M;
                HR[m] = 1;
            }
            int64_t reach = 0; for (int64_t m = 0; m < M; m++) reach += HR[m];
            if ((long double)reach > 0.60L * (long double)M) {   /* sieve won't help */
                free(HR); HR = NULL; M = 1;
            }
        } else { M = 1; }
    }
    if (M > 1) {
        __int128 tm = target % M; R.tmod = (int64_t)tm; if (R.tmod < 0) R.tmod += M;
        bool ok2 = true;
        for (int gi = 0; gi < ni && ok2; gi++) {
            int vi = ig[gi];
            gmod[vi] = (int*)malloc(sizeof(int) * (size_t)domain[vi]);
            if (!gmod[vi]) { ok2 = false; break; }
            for (int64_t d = 0; d < domain[vi]; d++) {
                __int128 g = R.gtab[vi][d] % M; int64_t m = (int64_t)g; if (m < 0) m += M;
                gmod[vi][d] = (int)m;
            }
            R.gmod[vi] = gmod[vi];
        }
        if (!ok2) { for (int i = 0; i < n; i++) { free(gmod[i]); } free(HR); HR = NULL; M = 1; }
    }
    R.M = M; R.HR = HR;

    /* Probe with the iterate group. */
    for (int i = 0; i < ni; i++) R.gv[i] = ig[i];
    R.ng = ni;
    for (int i = 0; i < n; i++) R.set[i] = false;
    psm_probe(&R, 0, 0, 0);

    free(R.keys);
    free(R.H);
    free(HR);
    for (int i = 0; i < n; i++) { free(gtab[i]); free(gmod[i]); }
    if (R.overflow) { st->nsol = 0; return false; }          /* did not exhaust */
    return true;                                             /* complete set (or {}) */
}
