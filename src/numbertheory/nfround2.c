/*
 * nfround2.c
 * ----------
 * Round 2 (Pohst-Zassenhaus) maximal-order computation for the number-field
 * layer: given a monic irreducible f of degree n, enlarge the equation order
 * Z[theta] to the ring of integers O_K, returned as an integral basis
 * (1/D)*W (W an n x n integer HNF matrix of theta-power numerators, D a shared
 * denominator) together with the index [O_K : Z[theta]] and d_K.
 *
 * This is Gate 1's constructive replacement: where nf_dedekind_p_maximal only
 * DECIDES p-maximality, this ENLARGES the order at each non-p-maximal prime.
 * The monogenic case (Z[theta] = O_K) returns W = I, D = 1 unchanged.
 *
 * Method (Cohen CCANT Alg. 6.1.3, Thm 6.1.3).  For each p with p^2 | disc(f):
 *   1. structure constants m_{ijk} of O (exact integers, w-basis), mod p.
 *   2. p-radical I_p = ker(Frobenius^k), k minimal with p^k >= n; a Z-basis of
 *      I_p by HNF of [ p*I ; lifted kernel ].
 *   3. ring of multipliers O' = (I_p : I_p): kernel over F_p of the map
 *      a -> (multiply-by-a on I_p/pI_p), assembled and HNF'd; O' = (1/p)*H'.
 *   repeat until the index stops dropping (p-maximal), then move to the next p.
 * All the delicate change-of-basis steps are exact over Q (H_I is singular mod
 * p); only the algebra maps are reduced mod p.
 *
 * Everything here is exact (fmpz/fmpq/nmod matrices); a proposal that cannot be
 * certified simply is not produced and the caller keeps declining.  Requires
 * FLINT.
 */
#include "numberfield_internal.h"

#include <stdlib.h>

#ifdef USE_FLINT
#include <stdio.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

static int r2dbg(void) { static int v = -1; if (v < 0) v = getenv("THUE_DEBUG") ? 1 : 0; return v; }

/* Extract the n non-zero rows of an m x n HNF matrix H (rank n) into the n x n
 * matrix out.  FLINT's HNF may place the zero rows at the top; we copy the rows
 * that are not entirely zero, in order.  Returns true iff exactly n found. */
static bool hnf_extract_square(fmpz_mat_t out, const fmpz_mat_t H, int n) {
    int m = fmpz_mat_nrows(H), got = 0;
    for (int i = 0; i < m && got < n; i++) {
        bool zero = true;
        for (int j = 0; j < n; j++) if (!fmpz_is_zero(fmpz_mat_entry(H, i, j))) { zero = false; break; }
        if (zero) continue;
        for (int j = 0; j < n; j++) fmpz_set(fmpz_mat_entry(out, got, j), fmpz_mat_entry(H, i, j));
        got++;
    }
    return got == n;
}

/* Exact structure constants of O = (1/D)W in the omega-basis:
 * Mc[i] is n x n with (Mc[i])[k][j] = m_{ijk}, where omega_i*omega_j =
 * sum_k m_{ijk} omega_k.  All integers.  f = defining poly (fmpz_poly).
 * Returns false on an unexpected non-integral solve (should not happen). */
static bool struct_consts(fmpz_mat_t* Mc, const fmpz_mat_t W, const fmpz_t D,
                          int n, const fmpz_poly_t f) {
    /* W as fmpq for the exact solves: solve  m * W = P/D  <=>  W^T m^T = (P/D)^T. */
    fmpq_mat_t Wt_q, rhs, sol;
    fmpq_mat_init(Wt_q, n, n); fmpq_mat_init(rhs, n, 1); fmpq_mat_init(sol, n, 1);
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) {   /* Wt_q = W^T */
        fmpz_set(fmpq_numref(fmpq_mat_entry(Wt_q, a, b)), fmpz_mat_entry(W, b, a));
        fmpz_one(fmpq_denref(fmpq_mat_entry(Wt_q, a, b)));
    }

    fmpz_poly_t vi, vj, prod; fmpz_poly_init(vi); fmpz_poly_init(vj); fmpz_poly_init(prod);
    bool ok = true;
    for (int i = 0; i < n && ok; i++) {
        for (int j = 0; j < n && ok; j++) {
            /* P = (row i of W as poly) * (row j) mod f */
            fmpz_poly_zero(vi); fmpz_poly_zero(vj);
            for (int a = 0; a < n; a++) { fmpz_poly_set_coeff_fmpz(vi, a, fmpz_mat_entry(W, i, a));
                                          fmpz_poly_set_coeff_fmpz(vj, a, fmpz_mat_entry(W, j, a)); }
            fmpz_poly_mul(prod, vi, vj);
            fmpz_poly_rem(prod, prod, f);                 /* mod f (f monic) */
            /* rhs = (P/D)^T, length n */
            for (int a = 0; a < n; a++) {
                fmpz_t Pa; fmpz_init(Pa); fmpz_poly_get_coeff_fmpz(Pa, prod, a);
                fmpz_set(fmpq_numref(fmpq_mat_entry(rhs, a, 0)), Pa);
                fmpz_set(fmpq_denref(fmpq_mat_entry(rhs, a, 0)), D);
                fmpq_canonicalise(fmpq_mat_entry(rhs, a, 0));
                fmpz_clear(Pa);
            }
            if (!fmpq_mat_solve(sol, Wt_q, rhs)) { ok = false; break; }
            /* sol = m^T (must be integral); (Mc[i])[k][j] = m_k */
            for (int k = 0; k < n; k++) {
                if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(sol, k, 0)))) { ok = false; break; }
                fmpz_set(fmpz_mat_entry(Mc[i], k, j), fmpq_numref(fmpq_mat_entry(sol, k, 0)));
            }
        }
    }
    fmpz_poly_clear(vi); fmpz_poly_clear(vj); fmpz_poly_clear(prod);
    fmpq_mat_clear(Wt_q); fmpq_mat_clear(rhs); fmpq_mat_clear(sol);
    return ok;
}

/* Product u*v in O/pO (omega-coords over F_p): (u*v)_k = sum_i u_i (L[i] v)_k. */
static void alg_mul(nmod_mat_t out /*n x 1*/, const nmod_mat_t* L, int n,
                    const nmod_mat_t u, const nmod_mat_t v, nmod_t mod) {
    nmod_mat_t tmp; nmod_mat_init(tmp, n, 1, mod.n);
    nmod_mat_zero(out);
    for (int i = 0; i < n; i++) {
        ulong ui = nmod_mat_entry(u, i, 0);
        if (ui == 0) continue;
        nmod_mat_mul(tmp, L[i], v);                       /* L[i] * v */
        for (int k = 0; k < n; k++) {
            ulong t = nmod_mul(ui, nmod_mat_entry(tmp, k, 0), mod);
            nmod_mat_entry(out, k, 0) = nmod_add(nmod_mat_entry(out, k, 0), t, mod);
        }
    }
    nmod_mat_clear(tmp);
}

/* One p-maximization round.  On success sets (Wout,Dout) to O' (>= O), *index to
 * [O':O] (a power of p; 1 == O is p-maximal), and returns true.  false on an
 * internal failure (caller declines). */
static bool p_maximize_once(const fmpz_mat_t W, const fmpz_t D, int n,
                            const fmpz_poly_t f, ulong p,
                            fmpz_mat_t Wout, fmpz_t Dout, fmpz_t index) {
    bool ok = true;
    nmod_t mod; nmod_init(&mod, p);

    /* --- 1. structure constants (exact) + mod p --- */
    fmpz_mat_t* Mc = malloc(sizeof(fmpz_mat_t) * (size_t)n);
    for (int i = 0; i < n; i++) fmpz_mat_init(Mc[i], n, n);
    if (!struct_consts(Mc, W, D, n, f)) { ok = false; goto cleanup_mc; }

    nmod_mat_t* L = malloc(sizeof(nmod_mat_t) * (size_t)n);
    for (int i = 0; i < n; i++) {
        nmod_mat_init(L[i], n, n, p);
        for (int k = 0; k < n; k++) for (int j = 0; j < n; j++)
            nmod_mat_entry(L[i], k, j) = fmpz_fdiv_ui(fmpz_mat_entry(Mc[i], k, j), p);
    }

    /* --- 2. Frobenius matrix F (col j = omega_j^p), A = F^k, radical = ker A --- */
    int k = 1; { ulong q = p; while (q < (ulong)n) { q *= p; k++; } }   /* p^k >= n */
    nmod_mat_t F; nmod_mat_init(F, n, n, p);
    {
        /* col j of F = omega_j^p in O/pO: raise e_j to the p-th power by repeated
         * algebra-multiplication (p is small, p^2 | disc). */
        nmod_mat_t base, res, sq; nmod_mat_init(base, n, 1, p); nmod_mat_init(res, n, 1, p); nmod_mat_init(sq, n, 1, p);
        for (int j = 0; j < n; j++) {
            nmod_mat_zero(base); nmod_mat_entry(base, j, 0) = 1 % p;
            nmod_mat_set(res, base);                        /* res = base^1 */
            for (ulong e = 1; e < p; e++) { alg_mul(sq, (const nmod_mat_t*)L, n, res, base, mod); nmod_mat_set(res, sq); }
            for (int r = 0; r < n; r++) nmod_mat_entry(F, r, j) = nmod_mat_entry(res, r, 0);
        }
        nmod_mat_clear(base); nmod_mat_clear(res); nmod_mat_clear(sq);
    }
    nmod_mat_t A; nmod_mat_init(A, n, n, p);
    if (k == 1) nmod_mat_set(A, F); else nmod_mat_pow(A, F, (ulong)k);

    nmod_mat_t kerI; nmod_mat_init(kerI, n, n, p);
    slong d = nmod_mat_nullspace(kerI, A);               /* kerI cols 0..d-1 = I_p/pO basis */

    /* Z-basis of I_p: HNF of [ p*I_n ; lifted kernel cols ] */
    fmpz_mat_t GI; fmpz_mat_init(GI, n + (int)d, n);
    for (int i = 0; i < n; i++) fmpz_set_ui(fmpz_mat_entry(GI, i, i), p);
    for (int t = 0; t < (int)d; t++)
        for (int r = 0; r < n; r++) fmpz_set_ui(fmpz_mat_entry(GI, n + t, r), nmod_mat_entry(kerI, r, t));
    fmpz_mat_t HIh; fmpz_mat_init(HIh, n + (int)d, n); fmpz_mat_hnf(HIh, GI);
    fmpz_mat_t HI; fmpz_mat_init(HI, n, n);
    if (!hnf_extract_square(HI, HIh, n)) { ok = false; goto cleanup_2; }

    /* --- 3. ring of multipliers: M_i = mult-by-omega_i on I_p/pI_p (beta-basis) --- */
    /* beta_l = row l of HI (omega-coords).  M_i col l = (omega_i * beta_l) in
     * beta-coords mod p; omega_i*beta_l = Mc[i] * beta_l (exact), then solve
     * c * HI = that  (HI singular mod p -> exact solve, reduce after). */
    fmpq_mat_t HIt_q; fmpq_mat_init(HIt_q, n, n);        /* HI^T over Q */
    for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) {
        fmpz_set(fmpq_numref(fmpq_mat_entry(HIt_q, a, b)), fmpz_mat_entry(HI, b, a));
        fmpz_one(fmpq_denref(fmpq_mat_entry(HIt_q, a, b)));
    }
    nmod_mat_t* M = malloc(sizeof(nmod_mat_t) * (size_t)n);
    for (int i = 0; i < n; i++) nmod_mat_init(M[i], n, n, p);
    {
        fmpz_mat_t betal, prod_i; fmpz_mat_init(betal, n, 1); fmpz_mat_init(prod_i, n, 1);
        fmpq_mat_t rhs, csol; fmpq_mat_init(rhs, n, 1); fmpq_mat_init(csol, n, 1);
        for (int i = 0; i < n && ok; i++) {
            for (int l = 0; l < n && ok; l++) {
                for (int r = 0; r < n; r++) fmpz_set(fmpz_mat_entry(betal, r, 0), fmpz_mat_entry(HI, l, r));
                fmpz_mat_mul(prod_i, Mc[i], betal);       /* omega_i * beta_l, omega-coords (exact) */
                for (int r = 0; r < n; r++) { fmpz_set(fmpq_numref(fmpq_mat_entry(rhs, r, 0)), fmpz_mat_entry(prod_i, r, 0));
                                              fmpz_one(fmpq_denref(fmpq_mat_entry(rhs, r, 0))); }
                if (!fmpq_mat_solve(csol, HIt_q, rhs)) { ok = false; break; }
                for (int r = 0; r < n; r++) {
                    if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(csol, r, 0)))) { ok = false; break; }
                    nmod_mat_entry(M[i], r, l) = fmpz_fdiv_ui(fmpq_numref(fmpq_mat_entry(csol, r, 0)), p);
                }
            }
        }
        fmpz_mat_clear(betal); fmpz_mat_clear(prod_i);
        fmpq_mat_clear(rhs); fmpq_mat_clear(csol);
    }

    /* Phi : (n*n) x n over F_p, Phi[(r,c), i] = (M_i)[r][c].  ker = L/pO. */
    nmod_mat_t Phi, kerPhi; nmod_mat_init(Phi, n * n, n, p);
    for (int i = 0; i < n; i++)
        for (int r = 0; r < n; r++) for (int c = 0; c < n; c++)
            nmod_mat_entry(Phi, r * n + c, i) = nmod_mat_entry(M[i], r, c);
    nmod_mat_init(kerPhi, n, n, p);
    slong e = ok ? nmod_mat_nullspace(kerPhi, Phi) : 0;

    /* pO' = HNF([ p*I_n ; lifted kerPhi cols ]); O' = (1/p) H'. */
    fmpz_mat_t Gp; fmpz_mat_init(Gp, n + (int)e, n);
    for (int i = 0; i < n; i++) fmpz_set_ui(fmpz_mat_entry(Gp, i, i), p);
    for (int t = 0; t < (int)e; t++)
        for (int r = 0; r < n; r++) fmpz_set_ui(fmpz_mat_entry(Gp, n + t, r), nmod_mat_entry(kerPhi, r, t));
    fmpz_mat_t Hph; fmpz_mat_init(Hph, n + (int)e, n); fmpz_mat_hnf(Hph, Gp);
    fmpz_mat_t Hp; fmpz_mat_init(Hp, n, n);
    if (ok && !hnf_extract_square(Hp, Hph, n)) ok = false;

    if (ok) {
        /* index = p^n / det(Hp) */
        fmpz_t detHp, pn; fmpz_init(detHp); fmpz_init(pn);
        fmpz_mat_det(detHp, Hp);
        fmpz_abs(detHp, detHp);
        fmpz_ui_pow_ui(pn, p, (ulong)n);
        fmpz_divexact(index, pn, detHp);
        fmpz_clear(detHp); fmpz_clear(pn);

        /* O' in theta-coords: (1/(p*D)) * (Hp * W). */
        fmpz_mat_t WrawM; fmpz_mat_init(WrawM, n, n); fmpz_mat_mul(WrawM, Hp, W);
        fmpz_t Draw; fmpz_init(Draw); fmpz_mul_ui(Draw, D, p);
        /* reduce common gcd */
        fmpz_t g; fmpz_init_set(g, Draw);
        for (int a = 0; a < n && !fmpz_is_one(g); a++) for (int b = 0; b < n; b++) fmpz_gcd(g, g, fmpz_mat_entry(WrawM, a, b));
        if (fmpz_sgn(g) != 0 && !fmpz_is_one(g)) {
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) fmpz_divexact(fmpz_mat_entry(WrawM, a, b), fmpz_mat_entry(WrawM, a, b), g);
            fmpz_divexact(Draw, Draw, g);
        }
        /* re-canonicalize W to HNF */
        fmpz_mat_t Wh; fmpz_mat_init(Wh, n, n); fmpz_mat_hnf(Wh, WrawM);
        if (!hnf_extract_square(Wout, Wh, n)) ok = false;
        fmpz_set(Dout, Draw);
        fmpz_mat_clear(Wh); fmpz_clear(g); fmpz_clear(Draw); fmpz_mat_clear(WrawM);
    }

    if (r2dbg() && ok) { flint_printf("[round2] p=%wu k=%d dimI=%wd dimL=%wd index=", p, k, d, e); fmpz_print(index); flint_printf("\n"); }

    /* cleanup step 3 */
    nmod_mat_clear(Phi); nmod_mat_clear(kerPhi);
    fmpz_mat_clear(Gp); fmpz_mat_clear(Hph); fmpz_mat_clear(Hp);
    for (int i = 0; i < n; i++) nmod_mat_clear(M[i]);
    free(M);
    fmpq_mat_clear(HIt_q);
cleanup_2:
    fmpz_mat_clear(HI); fmpz_mat_clear(HIh); fmpz_mat_clear(GI);
    nmod_mat_clear(kerI); nmod_mat_clear(A); nmod_mat_clear(F);
    for (int i = 0; i < n; i++) nmod_mat_clear(L[i]);
    free(L);
cleanup_mc:
    for (int i = 0; i < n; i++) fmpz_mat_clear(Mc[i]);
    free(Mc);
    return ok;
}

/*
 * Public entry.  coeffs[0..n] the monic defining poly; primes[0..nprimes-1] the
 * primes with p^2 | disc(f) (single-word).  On success fills W_out (malloc'd
 * mpz_t[n*n], row-major, theta-power numerators), *D_out, d_K_out, *index_out,
 * and returns true.  false to DECLINE (some prime too large, or an internal
 * certification failure).  The monogenic case yields W = I, D = 1, index = 1.
 */
bool nf_round2_maximal_order(const mpz_t* coeffs, int n,
                             const mpz_t* primes, const int64_t* pexp, int nprimes,
                             mpz_t* W_out, mpz_t D_out, mpz_t dK_out, int64_t* index_out) {
    fmpz_poly_t f; fmpz_poly_init(f);
    { fmpz_t c; fmpz_init(c); for (int i = 0; i <= n; i++) { fmpz_set_mpz(c, coeffs[i]); fmpz_poly_set_coeff_fmpz(f, i, c); } fmpz_clear(c); }

    fmpz_mat_t W; fmpz_mat_init(W, n, n);
    for (int i = 0; i < n; i++) fmpz_one(fmpz_mat_entry(W, i, i));
    fmpz_t D; fmpz_init_set_ui(D, 1);

    bool ok = true;
    for (int pi = 0; pi < nprimes && ok; pi++) {
        if (pexp[pi] < 2) continue;                       /* p^2 must divide disc */
        if (!mpz_fits_ulong_p(primes[pi])) { ok = false; break; }   /* prime too large */
        ulong p = mpz_get_ui(primes[pi]);
        for (;;) {
            fmpz_mat_t Wn; fmpz_mat_init(Wn, n, n); fmpz_t Dn, idx; fmpz_init(Dn); fmpz_init(idx);
            bool r = p_maximize_once(W, D, n, f, p, Wn, Dn, idx);
            if (!r) { fmpz_mat_clear(Wn); fmpz_clear(Dn); fmpz_clear(idx); ok = false; break; }
            bool maximal = fmpz_is_one(idx);
            if (!maximal) { fmpz_mat_set(W, Wn); fmpz_set(D, Dn); }
            fmpz_mat_clear(Wn); fmpz_clear(Dn); fmpz_clear(idx);
            if (maximal) break;
        }
    }

    if (ok) {
        /* index = D^n / det(W); d_K = disc(f) / index^2. */
        fmpz_t detW, Dn_, idxT; fmpz_init(detW); fmpz_init(Dn_); fmpz_init(idxT);
        fmpz_mat_det(detW, W); fmpz_abs(detW, detW);
        fmpz_pow_ui(Dn_, D, (ulong)n);
        fmpz_divexact(idxT, Dn_, detW);
        *index_out = fmpz_get_si(idxT);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) fmpz_get_mpz(W_out[i*n + j], fmpz_mat_entry(W, i, j));
        fmpz_get_mpz(D_out, D);
        /* d_K = disc(f)/index^2 */
        fmpz_t discf; fmpz_init(discf); fmpz_poly_discriminant(discf, f);
        fmpz_t idx2; fmpz_init(idx2); fmpz_mul(idx2, idxT, idxT);
        fmpz_divexact(discf, discf, idx2);
        fmpz_get_mpz(dK_out, discf);
        fmpz_clear(discf); fmpz_clear(idx2); fmpz_clear(detW); fmpz_clear(Dn_); fmpz_clear(idxT);
    }

    fmpz_mat_clear(W); fmpz_clear(D); fmpz_poly_clear(f);
    return ok;
}

#else  /* !USE_FLINT */
bool nf_round2_maximal_order(const mpz_t* coeffs, int n,
                             const mpz_t* primes, const int64_t* pexp, int nprimes,
                             mpz_t* W_out, mpz_t D_out, mpz_t dK_out, int64_t* index_out) {
    (void)coeffs; (void)n; (void)primes; (void)pexp; (void)nprimes;
    (void)W_out; (void)D_out; (void)dK_out; (void)index_out; return false;
}
#endif
