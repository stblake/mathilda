/*
 * solveint_cubes.c
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
 *  Booker cube-root-mod-d engine for  x^3 + y^3 + z^3 == k.
 *
 *  Reference: A. R. Booker, "Cracking the problem with 33" (2019).
 *  For any solution, k - z^3 = x^3 + y^3 = (x+y)(x^2-xy+y^2), so the
 *  divisor d = |x+y| divides k - z^3, i.e.  z^3 ≡ k (mod d).  Instead
 *  of enumerating z and factoring the ~B^3 value k - z^3 (the classical
 *  divisor method, which caps the z-range), we enumerate the SMALL
 *  divisor d and solve the cube-root congruence, recovering z as a union
 *  of arithmetic progressions.  The pair {x,y} then follows in closed
 *  form from d and z, and every hit is verified exactly.
 *
 *  The primitive below computes ALL cube roots of k modulo a machine
 *  integer d (d <= SI_BK_DMAX), the correctness-critical core: a missing
 *  residue would silently drop solutions.  d is bounded so every value
 *  and intermediate fits in int64.
 * ================================================================== */

#define SI_BK_BRUTE_P   4096LL         /* brute-force prime moduli up to this   */
#define SI_BK_DMAX      3000000LL      /* max divisor d (coords to ~1e7, SPF-fed)*/
#define SI_BK_ROOTCAP   4096           /* max cube roots mod d before declining  */
#define SI_BK_SOLCAP    200000         /* decline if the box holds a huge family */
#define SI_BK_MAXNODES  1000000000LL   /* part-2 candidate backstop (decline)    */

/* Modular inverse of a mod m (both positive, gcd(a,m)==1), int64. */
static int64_t si_inv_i64(int64_t a, int64_t m) {
    int64_t t = 0, newt = 1, r = m, newr = a % m;
    if (newr < 0) newr += m;
    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (r != 1) return 0;            /* not invertible */
    if (t < 0) t += m;
    return t;
}


/* (a*b) mod m for 0<=a,b<m<=~2.6e5: product < 7e10 fits int64. */
static int64_t si_mulmod_i64(int64_t a, int64_t b, int64_t m) {
    return (a * b) % m;
}


/* a^e mod m by square-and-multiply, int64 (m small). */
static int64_t si_powmod_i64(int64_t a, int64_t e, int64_t m) {
    int64_t r = 1 % m; a %= m; if (a < 0) a += m;
    while (e > 0) { if (e & 1) r = si_mulmod_i64(r, a, m); a = si_mulmod_i64(a, a, m); e >>= 1; }
    return r;
}


/* Trial-division factorisation of a positive int64 d (d <= SI_BK_DMAX, so the
 * largest prime factor is <= sqrt bound; cheap).  Fills primes/exps, sets *np.
 * Capacity of the arrays must be >= 16 (ample for d <= 2.6e5). */
static void si_factor_i64(int64_t d, int64_t* primes, unsigned* exps, int* np) {
    *np = 0;
    for (int64_t p = 2; p * p <= d; p++) {
        if (d % p) continue;
        unsigned e = 0;
        while (d % p == 0) { d /= p; e++; }
        primes[*np] = p; exps[*np] = e; (*np)++;
    }
    if (d > 1) { primes[*np] = d; exps[*np] = 1; (*np)++; }
}


/* Factor d via a precomputed smallest-prime-factor sieve: O(number of prime
 * factors) instead of O(sqrt d), so the divisor loop can range far wider. */
static void si_factor_spf(int64_t d, const int32_t* spf, int64_t* primes,
                          unsigned* exps, int* np) {
    *np = 0;
    while (d > 1) {
        int64_t p = spf[d]; unsigned e = 0;
        while (d % p == 0) { d /= p; e++; }
        primes[*np] = p; exps[*np] = e; (*np)++;
    }
}


/* Build the smallest-prime-factor table for [0, n]; spf[i] = least prime | i
 * (spf[0]=spf[1]=0).  Caller owns the returned buffer (NULL on OOM). */
static int32_t* si_build_spf(int64_t n) {
    int32_t* spf = (int32_t*)calloc((size_t)n + 1, sizeof(int32_t));
    if (!spf) return NULL;
    for (int64_t i = 2; i <= n; i++) {
        if (spf[i]) continue;                         /* already has a factor */
        for (int64_t j = i; j <= n; j += i)
            if (!spf[j]) spf[j] = (int32_t)i;
    }
    return spf;
}


/* ALL cube roots of a (mod p), p an odd prime, a in [0,p).  Writes the
 * distinct roots into out[0..*n) (out capacity >= 3).  Always succeeds. */
static void si_croots_mod_p(int64_t* out, int* n, int64_t a, int64_t p) {
    *n = 0;
    a %= p; if (a < 0) a += p;
    if (a == 0) { out[(*n)++] = 0; return; }
    if (p == 3) { out[(*n)++] = a; return; }         /* z^3 ≡ z (mod 3) */
    if (p % 3 == 2) {                                /* cubing is a bijection */
        int64_t inv3 = si_inv_i64(3, p - 1);
        out[(*n)++] = si_powmod_i64(a, inv3, p);
        return;
    }
    /* p ≡ 1 (mod 3): 0 or 3 roots. */
    int64_t ecr = (p - 1) / 3;
    if (si_powmod_i64(a, ecr, p) != 1) return;       /* not a cubic residue */
    /* p-1 = 3^s * t, 3 ∤ t. */
    int64_t t = p - 1; int s = 0;
    while (t % 3 == 0) { t /= 3; s++; }
    /* cubic non-residue nn -> 3-Sylow generator b (order 3^s), gamma (order 3). */
    int64_t nn = 2;
    while (si_powmod_i64(nn, ecr, p) == 1) nn++;
    int64_t b = si_powmod_i64(nn, t, p);
    int64_t binv = si_inv_i64(b, p);
    int64_t threepow_sm1 = 1; for (int i = 0; i < s - 1; i++) threepow_sm1 *= 3;
    int64_t gamma = si_powmod_i64(b, threepow_sm1, p);   /* primitive cube root of 1 */
    /* x1 = a^(3^{-1} mod t); then x1^3 = a * beta^w with beta = a^t (in 3-Sylow),
     * w = (3*m - 1)/t.  Correct x1 by b^{-c}, c = dlog_b(beta^w)/3. */
    int64_t m = si_inv_i64(3, t);
    int64_t x1 = si_powmod_i64(a, m, p);
    int64_t beta = si_powmod_i64(a, t, p);
    int64_t w = (3 * m - 1) / t;                     /* exact, w in {0,1,2} */
    int64_t betaw = si_powmod_i64(beta, w, p);
    /* Pohlig-Hellman discrete log E of betaw base b in the 3-Sylow. */
    int64_t E = 0, cur = betaw, pw3 = 1;             /* pw3 = 3^i */
    for (int i = 0; i < s; i++) {
        int64_t expo = 1; for (int j = 0; j < s - 1 - i; j++) expo *= 3;
        int64_t h = si_powmod_i64(cur, expo, p);
        int digit = (h == 1) ? 0 : (h == gamma ? 1 : 2);
        if (digit) {
            E += (int64_t)digit * pw3;
            cur = si_mulmod_i64(cur, si_powmod_i64(binv, (int64_t)digit * pw3, p), p);
        }
        pw3 *= 3;
    }
    if (E % 3 != 0) return;                           /* inconsistent: no root */
    int64_t x0 = si_mulmod_i64(x1, si_powmod_i64(binv, E / 3, p), p);
    /* three roots x0 * {1, gamma, gamma^2}. */
    out[(*n)++] = x0;
    out[(*n)++] = si_mulmod_i64(x0, gamma, p);
    out[(*n)++] = si_mulmod_i64(x0, si_mulmod_i64(gamma, gamma, p), p);
}


/* ALL cube roots of a (mod p^e), a in [0,pe).  A large prime modulus (e==1)
 * takes the analytic path (<=3 roots); every other prime power is brute-forced,
 * which is correct for all p=3 and p|a cases -- but those can yield MANY roots
 * (e.g. z^3 ≡ 0 mod 2^18 has 4096 roots), so the count is capped: returns the
 * number of roots written, or -1 if it would exceed `cap` (caller declines). */
static int si_croots_mod_pe(int64_t* out, int64_t a, int64_t p,
                            unsigned e, int64_t pe, int cap) {
    if (e == 1 && p > SI_BK_BRUTE_P) { int n; si_croots_mod_p(out, &n, a, p); return n; }
    int n = 0;
    a %= pe; if (a < 0) a += pe;
    for (int64_t z = 0; z < pe; z++)
        if (si_mulmod_i64(si_mulmod_i64(z, z, pe), z, pe) == a) {
            if (n >= cap) return -1;
            out[n++] = z;
        }
    return n;
}


/* ALL residues r in [0,d) with r^3 ≡ k (mod d).  Returns the count, or -1 if
 * the root set would exceed `cap` (caller must then decline: never truncate).
 * `out` must have capacity `cap`.  `spf`, if non-NULL, factors d in O(log d)
 * (a smallest-prime-factor table covering d); otherwise trial division is used. */
static int si_all_cube_roots_mod(int64_t* out, int cap, const mpz_t k, int64_t d,
                                 const int32_t* spf) {
    if (d == 1) { if (cap < 1) return -1; out[0] = 0; return 1; }
    if (cap > SI_BK_ROOTCAP) cap = SI_BK_ROOTCAP;
    int64_t kd = mpz_fdiv_ui(k, (unsigned long)d);   /* k mod d, in [0,d) */
    int64_t primes[16]; unsigned exps[16]; int np;
    if (spf) si_factor_spf(d, spf, primes, exps, &np);
    else     si_factor_i64(d, primes, exps, &np);
    /* Accumulate via CRT over prime powers. acc holds residues mod `mod_acc`. */
    static int64_t acc[SI_BK_ROOTCAP], newacc[SI_BK_ROOTCAP], rp[SI_BK_ROOTCAP];
    int nacc = 1; acc[0] = 0; int64_t mod_acc = 1;
    for (int f = 0; f < np; f++) {
        int64_t p = primes[f]; unsigned e = exps[f];
        int64_t pe = 1; for (unsigned j = 0; j < e; j++) pe *= p;
        int nr = si_croots_mod_pe(rp, kd % pe, p, e, pe, SI_BK_ROOTCAP);
        if (nr < 0) return -1;                        /* too many roots mod pe */
        if (nr == 0) return 0;                        /* no root mod pe -> none */
        /* CRT-combine acc (mod mod_acc) with rp (mod pe). */
        if ((long long)nacc * nr > cap) return -1;
        int64_t ninv = si_inv_i64(mod_acc % pe, pe);
        int nnew = 0;
        for (int i = 0; i < nacc; i++)
            for (int j = 0; j < nr; j++) {
                int64_t diff = rp[j] - (acc[i] % pe); diff %= pe; if (diff < 0) diff += pe;
                int64_t t = si_mulmod_i64(diff, ninv, pe);
                newacc[nnew++] = acc[i] + mod_acc * t;   /* < mod_acc*pe <= d */
            }
        for (int i = 0; i < nnew; i++) acc[i] = newacc[i];
        nacc = nnew; mod_acc *= pe;
    }
    if (nacc > cap) return -1;
    for (int i = 0; i < nacc; i++) out[i] = acc[i];
    return nacc;
}


/* Internal test hook: Solve`CubeRootsMod[k, d] -> sorted list of all r in
 * [0,d) with r^3 ≡ k (mod d).  Used by the test suite to cross-check the
 * primitive against brute force. */
Expr* builtin_cube_roots_mod(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* ke = res->data.function.args[0];
    Expr* de = res->data.function.args[1];
    if (!expr_is_integer_like(ke) || de->type != EXPR_INTEGER || de->data.integer < 1)
        return NULL;
    int64_t d = de->data.integer;
    if (d > SI_BK_DMAX) return NULL;
    mpz_t k; mpz_init(k); expr_to_mpz(ke, k);
    enum { CAP = 4096 };
    int64_t roots[CAP];
    int nr = si_all_cube_roots_mod(roots, CAP, k, d, NULL);
    mpz_clear(k);
    if (nr < 0) return NULL;
    /* sort ascending (simple insertion; nr is small). */
    for (int i = 1; i < nr; i++) {
        int64_t v = roots[i]; int j = i - 1;
        while (j >= 0 && roots[j] > v) { roots[j + 1] = roots[j]; j--; }
        roots[j + 1] = v;
    }
    Expr** elems = (nr > 0) ? (Expr**)malloc((size_t)nr * sizeof(Expr*)) : NULL;
    for (int i = 0; i < nr; i++) elems[i] = expr_new_integer(roots[i]);
    Expr* list = mk_list(elems, (size_t)nr);
    free(elems);
    return list;
}


/* Recover the pair {x,y} with x+y = ±d and x^3+y^3 = m (m != 0, d | m), the
 * closed form  {x,y} = (sgn(m)*d ± sqrt((4|m|/d - d^2)/3)) / 2, and emit any
 * that lands in the box.  vals already carries the modular variable; iv/jv are
 * the pair positions. */
static void si_bk_emit_pair(SICtx* c, SearchState* st, int iv, int jv,
                            int64_t d, __int128 m, int64_t* vals) {
    if (m == 0) return;                              /* x=-y family: part1 covers it */
    __int128 am = (m < 0) ? -m : m;
    if (am % d != 0) return;                          /* not a divisor (guard) */
    __int128 q = am / d;                              /* x^2 - xy + y^2 > 0 */
    __int128 disc = 4 * q - (__int128)d * d;
    if (disc < 0 || disc % 3 != 0) return;
    disc /= 3;
    int64_t R = si_isqrt_i128(disc);
    if ((__int128)R * R != disc) return;              /* not a perfect square */
    int64_t s = (m > 0) ? d : -d;                     /* s = x + y = sgn(m)*d */
    for (int rs = 1; rs >= -1; rs -= 2) {
        int64_t num = s + (int64_t)rs * R;
        if (num & 1) continue;                        /* x must be integral */
        int64_t xv = num / 2, yv = s - xv;
        vals[iv] = xv; vals[jv] = yv;
        if (si_verify(c, vals)) emit_full(st, vals);
        if (R == 0) break;                            /* double root: one pair */
    }
}


/* Detect  x^3 + y^3 + z^3 == k  (unit coefficients, k a nonzero integer, three
 * box-bounded variables) and solve it by the Booker cube-root-mod-d method.
 *
 * For any solution the divisor d = |a+b| of a pair {a,b} divides k - c^3 of the
 * third variable c, i.e. c^3 ≡ k (mod d).  So instead of LINEARLY enumerating
 * the third variable (the classical outer loop, whose range is the bottleneck),
 * we enumerate the SMALL divisor d up to alpha*B (alpha = 2^(1/3)-1), cube-root
 * k mod d to pin c to arithmetic progressions, and recover {a,b} in closed
 * form.  This lifts the reachable coordinate range from the ~3e5 outer cap
 * toward ~1e6.  Completeness (across the three roles):
 *   - part 1: every solution whose smallest |coordinate| <= T = floor(k^{1/3})
 *     is found by enumerating that small coordinate and divisor-solving the pair
 *     (this also covers the k=c^3 => a=-b family);
 *   - part 2: every solution whose coordinates all exceed T satisfies Booker's
 *     bound d < alpha*|c| <= alpha*B for the pairing of its two largest, so it
 *     is found in the role where c is the smallest coordinate.
 * Returns true (a complete solution set emitted) only when the box is provably
 * covered by the budget; otherwise declines so the normal pipeline continues. */
bool si_solve_three_cubes_booker(SICtx* c, SearchState* st) {
    if (c->n != 3 || c->neq != 1) return false;
    for (int i = 0; i < 3; i++) if (!(c->has_lo[i] && c->has_hi[i])) return false;
    const MPoly* eq = c->eq[0];

    /* Shape: each variable appears once as v^3 with coefficient +/-1, plus one
     * constant term -k.  K accumulates -constant = k; sgn[i] records the sign of
     * v_i^3 (a -1 is handled by the substitution u_i = -v_i below, which turns
     * the system into the pure  u_x^3 + u_y^3 + u_z^3 == k  over a mirrored box:
     * e.g.  x^3 + y^3 - z^3 == 227  is the sum of three cubes for 227). */
    mpz_t K; mpz_init(K); mpz_set_ui(K, 0);
    int seen[3] = {0, 0, 0}, sgn[3] = {1, 1, 1}; bool ok = true;
    for (size_t t = 0; t < eq->n_terms && ok; t++) {
        const int* ex = eq->exps + t * 3;
        int nz = -1, cnt = 0;
        for (int v = 0; v < 3; v++) if (ex[v] > 0) { nz = v; cnt++; }
        if (cnt == 0) { mpz_sub(K, K, eq->coefs[t]); continue; }   /* constant */
        if (cnt > 1 || ex[nz] != 3 || seen[nz]) { ok = false; break; }
        if (mpz_cmp_ui(eq->coefs[t], 1) == 0) sgn[nz] = 1;
        else if (mpz_cmp_si(eq->coefs[t], -1) == 0) sgn[nz] = -1;
        else { ok = false; break; }                               /* coeff not +/-1 */
        seen[nz] = 1;
    }
    for (int i = 0; i < 3 && ok; i++) if (!seen[i]) ok = false;
    /* |k| must be modest so part 1's pair factoring stays cheap. */
    if (ok && mpz_sizeinbase(K, 2) > 30) ok = false;              /* |k| < ~1e9 */
    if (!ok || mpz_sgn(K) == 0) { mpz_clear(K); return false; }
    bool signed_case = (sgn[0] < 0 || sgn[1] < 0 || sgn[2] < 0);
    /* The u_i = sgn_i v_i substitution rewrites the equation with unit
     * coefficients and mirrors the flipped variable's box; it only stays exact
     * for a pure box (no orderings / disequations, every constraint captured),
     * since a mixed-sign ordering is not an ordering in u-space. */
    if (signed_case && (c->n_ord != 0 || c->n_neq != 0 || c->n_abs_ord != 0 || !c->all_captured)) {
        mpz_clear(K); return false;
    }
    int64_t Ki = mpz_get_si(K);
    int64_t absk = Ki < 0 ? -Ki : Ki;

    /* Box radius B and divisor bound Dmax = ceil(alpha*B).  The u_i = sgn_i v_i
     * mirror preserves |bounds|, so B is the same in u-space -- compute from c. */
    int64_t B = 0;
    for (int i = 0; i < 3; i++) {
        int64_t a = c->hi[i] > -c->lo[i] ? c->hi[i] : -c->lo[i];
        if (a > B) B = a;
    }
    bool force = getenv("MATHILDA_BK_FORCE") != NULL;
    /* The cube-root-mod-d method's O(alpha*B * roots) work beats the leaf
     * search's O(B^2) and the divisor method's O(B * factoring) across the whole
     * bounded range -- measured ~10-400x faster from B~500 up, and validated
     * to return the identical solution set (Booker vs the leaf/divisor pipeline)
     * over a wide (k, box) sweep incl. the (a,-a,cbrt k) family and signed
     * equations.  So engage it for any non-trivial box; only a truly tiny box
     * (|coord| <= 100), already solved in well under a millisecond by the
     * leaf/mitm paths, is left alone.  The upper bound is Dmax <= SI_BK_DMAX
     * below (past which int64 coverage is not guaranteed). */
    if (!force && 2 * B <= 200) { mpz_clear(K); return false; }
    int64_t Dmax = (int64_t)(0.25992104989 * (double)B) + 2;
    if (Dmax > SI_BK_DMAX) { mpz_clear(K); return false; }

    /* T = floor(|k|^{1/3}). */
    int64_t T = 0; while ((T + 1) * (T + 1) * (T + 1) <= absk) T++;

    /* Working context cc: the pure  u_x^3+u_y^3+u_z^3 == k  problem.  The all-+
     * case uses c directly; a sign flips the corresponding box and swaps the
     * equation MPoly for unit coefficients (freed at the end).  Built after the
     * early declines so those paths need not free it. */
    SICtx cu; SICtx* cc = c; MPoly* newpoly = NULL;
    if (signed_case) {
        cu = *c;
        for (int i = 0; i < 3; i++)
            if (sgn[i] < 0) { cu.lo[i] = -c->hi[i]; cu.hi[i] = -c->lo[i]; }
        MPoly* p = mpoly_monomial(3, 0, 3, 1);
        for (int i = 1; i < 3; i++) {
            MPoly* tm = mpoly_monomial(3, i, 3, 1);
            MPoly* s = mpoly_add(p, tm); mpoly_free(p); mpoly_free(tm); p = s;
        }
        MPoly* kk = mpoly_from_mpz(3, K);
        newpoly = mpoly_sub(p, kk); mpoly_free(p); mpoly_free(kk);
        cu.eq[0] = newpoly; cu.neq = 1;
        cc = &cu;
    }

    /* Smallest-prime-factor sieve over [0, Dmax]: factors every divisor in
     * O(log d), so the enumeration is dominated by the candidate loop, not by
     * factoring.  Sized to Dmax, so small boxes stay cheap. */
    int32_t* spf = si_build_spf(Dmax);
    if (!spf) { if (newpoly) mpoly_free(newpoly); mpz_clear(K); return false; }

    st->max_visits = SI_BK_MAXNODES;
    int64_t vals[SI_MAX_VARS];

    /* ---- part 1: smallest coordinate <= T (per role), divisor-solve pair. ---- */
    for (int r = 0; r < 3 && !st->overflow; r++) {
        int iv = (r + 1) % 3, jv = (r + 2) % 3;
        int64_t lo = cc->lo[r] > -T ? cc->lo[r] : -T;
        int64_t hi = cc->hi[r] <  T ? cc->hi[r] :  T;
        mpz_t m, v3; mpz_init(m); mpz_init(v3);
        for (int64_t v = lo; v <= hi && !st->overflow; v++) {
            mpz_set_si(v3, v); mpz_mul_si(v3, v3, v); mpz_mul_si(v3, v3, v);
            mpz_sub(m, K, v3);                         /* m = k - v^3 */
            for (int q = 0; q < 3; q++) vals[q] = 0;
            vals[r] = v;
            si_two_power_solve(cc, st, iv, jv, 3, m, vals);
        }
        mpz_clears(m, v3, NULL);
    }

    /* A box that holds a huge parametric family (e.g. (a,-a,c) when k=c^3) is
     * better left to a symbolic/parametric treatment than enumerated; decline
     * rather than materialise hundreds of thousands of tuples. */
    bool incomplete = st->nsol > SI_BK_SOLCAP;

    /* ---- part 2: all coordinates > T, Booker divisor enumeration. ---- */
    int64_t roots[SI_BK_ROOTCAP];
    for (int64_t d = 1; d <= Dmax && !st->overflow && !incomplete; d++) {
        if (st->nsol > SI_BK_SOLCAP) { incomplete = true; break; }
        int nr = si_all_cube_roots_mod(roots, SI_BK_ROOTCAP, K, d, spf);
        if (nr < 0) { incomplete = true; break; }     /* too many roots: bail */
        if (nr == 0) continue;
        for (int r = 0; r < 3 && !st->overflow; r++) {
            int iv = (r + 1) % 3, jv = (r + 2) % 3;
            int64_t lo = cc->lo[r], hi = cc->hi[r];
            for (int ridx = 0; ridx < nr; ridx++) {
                int64_t start = lo + ((((roots[ridx] - lo) % d) + d) % d);
                for (int64_t v = start; v <= hi && !st->overflow; v += d) {
                    if (++st->visits > st->max_visits) { st->overflow = true; break; }
                    if (v >= -T && v <= T) continue;   /* handled by part 1 */
                    __int128 v3 = (__int128)v * v * v; /* |v| up to ~1e7 */
                    __int128 m = (__int128)Ki - v3;
                    for (int q = 0; q < 3; q++) vals[q] = 0;
                    vals[r] = v;
                    si_bk_emit_pair(cc, st, iv, jv, d, m, vals);
                }
            }
        }
    }

    free(spf);
    /* Candidates were enumerated in u-space; map each back to the original
     * variables (v_i = sgn_i u_i) before the result is built. */
    if (!incomplete && !st->overflow && signed_case)
        for (int i = 0; i < st->nsol; i++) {
            int64_t* row = st->sols + (size_t)i * 3;
            for (int q = 0; q < 3; q++) if (sgn[q] < 0) row[q] = -row[q];
        }
    if (newpoly) mpoly_free(newpoly);
    mpz_clear(K);
    if (incomplete || st->overflow) { st->nsol = 0; return false; }
    return true;
}
