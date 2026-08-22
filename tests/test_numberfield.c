/*
 * test_numberfield.c
 * ------------------
 * Unit tests for the number-field layer's Gate 1 (maximal-order
 * certification via Dedekind's criterion) and the exact Sturm signature.
 * Directly exercises the C API in src/numbertheory/numberfield.{c,h};
 * no evaluator is involved.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "numberfield.h"
#include "nfunits.h"
#ifdef USE_FLINT
#include "numberfield_internal.h"   /* nf_norm_int, for the exact |N|=1 recheck */
#endif

/* Build a coefficient array {c0, c1, ..., c_deg} from int64 values.
 * Returns a malloc'd mpz_t[deg+1]; free with free_coeffs(). */
static mpz_t* make_coeffs(const long* vals, int deg) {
    mpz_t* c = malloc(sizeof(mpz_t) * (size_t)(deg + 1));
    for (int i = 0; i <= deg; i++) mpz_init_set_si(c[i], vals[i]);
    return c;
}
static void free_coeffs(mpz_t* c, int deg) {
    for (int i = 0; i <= deg; i++) mpz_clear(c[i]);
    free(c);
}

/* Signature: number of real roots via the exact Sturm sequence. */
static void test_signature(void) {
    int r1;
    { const long v[] = {-2, 0, 0, 1};      mpz_t* c = make_coeffs(v, 3);  /* t^3 - 2   */
      assert(nf_signature(c, 3, &r1) && r1 == 1); free_coeffs(c, 3); }
    { const long v[] = {1, -3, 0, 1};      mpz_t* c = make_coeffs(v, 3);  /* t^3-3t+1  */
      assert(nf_signature(c, 3, &r1) && r1 == 3); free_coeffs(c, 3); }
    { const long v[] = {-2, 0, 0, 0, 1};   mpz_t* c = make_coeffs(v, 4);  /* t^4 - 2   */
      assert(nf_signature(c, 4, &r1) && r1 == 2); free_coeffs(c, 4); }
    { const long v[] = {1, 0, 1};          mpz_t* c = make_coeffs(v, 2);  /* t^2 + 1   */
      assert(nf_signature(c, 2, &r1) && r1 == 0); free_coeffs(c, 2); }
    { const long v[] = {-2, 0, 1};         mpz_t* c = make_coeffs(v, 2);  /* t^2 - 2   */
      assert(nf_signature(c, 2, &r1) && r1 == 2); free_coeffs(c, 2); }
    { const long v[] = {-1, 0, 0, 0, 0, 1};mpz_t* c = make_coeffs(v, 5);  /* t^5 - 1   */
      assert(nf_signature(c, 5, &r1) && r1 == 1); free_coeffs(c, 5); }
    printf("  signature: OK\n");
}

/* Dedekind primitive on the classic non-monogenic cubic. */
static void test_dedekind(void) {
    if (!nf_available()) { printf("  dedekind: SKIP (no FLINT)\n"); return; }
    bool ok;
    /* t^3 - 2: Z[cbrt 2] is maximal, so p-maximal at every p (incl. 2,3). */
    { const long v[] = {-2, 0, 0, 1}; mpz_t* c = make_coeffs(v, 3);
      mpz_t p; mpz_init_set_ui(p, 2);
      assert(nf_dedekind_p_maximal(c, 3, p, &ok) && ok);
      mpz_set_ui(p, 3);
      assert(nf_dedekind_p_maximal(c, 3, p, &ok) && ok);
      mpz_clear(p); free_coeffs(c, 3); }
    /* t^3 - t^2 - 2t - 8 (disc = -2012 = -4*503): index 2 at p=2 => NOT 2-maximal. */
    { const long v[] = {-8, -2, -1, 1}; mpz_t* c = make_coeffs(v, 3);
      mpz_t p; mpz_init_set_ui(p, 2);
      assert(!nf_dedekind_p_maximal(c, 3, p, &ok) && ok);   /* proven non-maximal */
      mpz_clear(p); free_coeffs(c, 3); }
    printf("  dedekind: OK\n");
}

/* Full field creation: monogenic accepts (with correct invariants),
 * non-monogenic and reducible decline. */
static void test_create(void) {
    if (!nf_available()) { printf("  create: SKIP (no FLINT)\n"); return; }
    mpz_t d; mpz_init(d);

    { const long v[] = {-2, 0, 0, 1}; mpz_t* c = make_coeffs(v, 3);   /* Q(cbrt 2) */
      NumberField* K = nf_field_create(c, 3);
      assert(K != NULL);
      assert(nf_degree(K) == 3 && nf_r1(K) == 1 && nf_r2(K) == 1 && nf_unit_rank(K) == 1);
      nf_disc(K, d); assert(mpz_cmp_si(d, -108) == 0);
      nf_field_free(K); free_coeffs(c, 3); }

    { const long v[] = {1, -3, 0, 1}; mpz_t* c = make_coeffs(v, 3);   /* cyclic cubic, cond 9 */
      NumberField* K = nf_field_create(c, 3);
      assert(K != NULL);
      assert(nf_degree(K) == 3 && nf_r1(K) == 3 && nf_r2(K) == 0 && nf_unit_rank(K) == 2);
      nf_disc(K, d); assert(mpz_cmp_si(d, 81) == 0);
      nf_field_free(K); free_coeffs(c, 3); }

    { const long v[] = {-2, 0, 0, 0, 1}; mpz_t* c = make_coeffs(v, 4); /* Q(2^{1/4}) */
      NumberField* K = nf_field_create(c, 4);
      assert(K != NULL);
      assert(nf_degree(K) == 4 && nf_r1(K) == 2 && nf_r2(K) == 1 && nf_unit_rank(K) == 2);
      nf_field_free(K); free_coeffs(c, 4); }

    /* Non-monogenic: Dedekind non-maximal at 2, now enlarged by Round 2 (M3):
     * creates with d_K = -503 (== disc(f)/index^2, index 2). */
    { const long v[] = {-8, -2, -1, 1}; mpz_t* c = make_coeffs(v, 3);
      NumberField* K = nf_field_create(c, 3);
      assert(K != NULL);
      nf_disc(K, d); assert(mpz_cmp_si(d, -503) == 0);
      nf_field_free(K); free_coeffs(c, 3); }

    /* Reducible: t^3 - t = t(t-1)(t+1) => DECLINE. */
    { const long v[] = {0, -1, 0, 1}; mpz_t* c = make_coeffs(v, 3);
      NumberField* K = nf_field_create(c, 3);
      assert(K == NULL);
      free_coeffs(c, 3); }

    /* Non-monic input => DECLINE. */
    { const long v[] = {-2, 0, 0, 2}; mpz_t* c = make_coeffs(v, 3);
      NumberField* K = nf_field_create(c, 3);
      assert(K == NULL);
      free_coeffs(c, 3); }

    mpz_clear(d);
    printf("  create: OK\n");
}

/* Certified fundamental units + regulator (Gate 2). */
static void check_units(const long* v, int deg, int exp_rank,
                        double rlo, double rhi, const char* name) {
    mpz_t* c = make_coeffs(v, deg);
    NumberField* K = nf_field_create(c, deg);
    assert(K != NULL);
    NFUnits* U = nf_fundamental_units(K);
    assert(U != NULL);                                   /* must certify */
    int rank = nf_units_rank(U);
    double R = nf_units_regulator(U);
    printf("  %s: rank=%d regulator=%.6f\n", name, rank, R);
    assert(rank == exp_rank);
    assert(R > rlo && R < rhi);
    /* every returned unit must genuinely have norm +/-1 (exact recheck) */
    for (int j = 0; j < rank; j++) {
        const mpz_t* uc = nf_units_coords(U, j);
        int nz = 0; for (int i = 0; i < deg; i++) if (mpz_sgn(uc[i]) != 0) nz = 1;
        assert(nz);
#ifdef USE_FLINT
        mpz_t nrm; mpz_init(nrm);
        assert(nf_norm_int(K, uc, nrm));                 /* norm is a rational integer */
        assert(mpz_cmpabs_ui(nrm, 1) == 0);              /* |N| == 1: genuinely a unit */
        mpz_clear(nrm);
#endif
    }
    nf_units_free(U);
    nf_field_free(K);
    free_coeffs(c, deg);
}

static void test_units(void) {
    if (!nf_available()) { printf("  units: SKIP (no FLINT)\n"); return; }
    /* --- small-regulator fields: found by the coefficient-box search --- */
    /* Q(cbrt 2): rank 1, regulator ~ 1.3474. */
    { const long v[] = {-2, 0, 0, 1};    check_units(v, 3, 1, 1.30, 1.40, "Q(2^1/3)"); }
    /* cyclic cubic t^3-3t+1, disc 81: rank 2, regulator ~ 0.8494. */
    { const long v[] = {1, -3, 0, 1};    check_units(v, 3, 2, 0.80, 0.90, "cyclic cubic"); }
    /* Q(2^{1/4}): rank 2. */
    { const long v[] = {-2, 0, 0, 0, 1}; check_units(v, 4, 2, 0.5, 20.0, "Q(2^1/4)"); }

    /* --- large-regulator complex cubics: the box FAILS, Voronoi's chain of
     * minima proposes the unit, p-saturation certifies it.  Tight regulator
     * windows (PARI/GP bnfinit .reg) catch any non-fundamental power (which
     * would give k*R) -- the exact test that Voronoi found the RIGHT unit. */
    /* Q(cbrt 15): fund. unit has a coord 30 (> box 12); R = 9.6929517. */
    { const long v[] = {-15, 0, 0, 1};   check_units(v, 3, 1, 9.6929, 9.6930, "Q(15^1/3)"); }
    /* Q(cbrt 42): coord 42; R = 11.0589054. */
    { const long v[] = {-42, 0, 0, 1};   check_units(v, 3, 1, 11.0588, 11.0590, "Q(42^1/3)"); }
    /* Q(cbrt 97): ~10-digit coords; R = 49.4921861. */
    { const long v[] = {-97, 0, 0, 1};   check_units(v, 3, 1, 49.491, 49.493, "Q(97^1/3)"); }
    /* Q(cbrt 41): 24-digit coords, R = 56.2893702 -- the stress case. */
    { const long v[] = {-41, 0, 0, 1};   check_units(v, 3, 1, 56.288, 56.290, "Q(41^1/3)"); }
    printf("  units: OK\n");
}

#ifdef USE_FLINT
/* Round 2 maximal order: check index and d_K against PARI/known values.
 * primes/exps give the disc factorization (only exp>=2 primes are processed). */
static void check_round2(const long* v, int deg, const long* pr, const long* ex, int np,
                         long exp_index, long exp_dK, const char* name) {
    mpz_t* c = make_coeffs(v, deg);
    mpz_t* primes = malloc(sizeof(mpz_t) * (size_t)np);
    int64_t* pexp = malloc(sizeof(int64_t) * (size_t)np);
    for (int i = 0; i < np; i++) { mpz_init_set_si(primes[i], pr[i]); pexp[i] = ex[i]; }
    mpz_t* W = malloc(sizeof(mpz_t) * (size_t)(deg * deg));
    for (int i = 0; i < deg * deg; i++) mpz_init(W[i]);
    mpz_t D, dK; mpz_init(D); mpz_init(dK);
    int64_t index = -1;
    bool ok = nf_round2_maximal_order((const mpz_t*)c, deg, (const mpz_t*)primes, pexp, np, W, D, dK, &index);
    printf("  %s: ok=%d index=%lld d_K=%s\n", name, (int)ok, (long long)index, mpz_get_str(NULL, 10, dK));
    assert(ok);
    assert(index == exp_index);
    assert(mpz_cmp_si(dK, exp_dK) == 0);
    for (int i = 0; i < deg * deg; i++) mpz_clear(W[i]);
    free(W);
    for (int i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes); free(pexp);
    mpz_clear(D); mpz_clear(dK);
    free_coeffs(c, deg);
}

static void test_round2(void) {
    if (!nf_available()) { printf("  round2: SKIP (no FLINT)\n"); return; }
    /* Dedekind cubic t^3-t^2-2t-8: disc -2012 = -2^2*503, index 2, d_K -503. */
    { const long v[] = {-8,-2,-1,1}; const long pr[] = {2,503}; const long ex[] = {2,1};
      check_round2(v, 3, pr, ex, 2, 2, -503, "t^3-t^2-2t-8"); }
    /* Q(cbrt20): disc -10800 = -2^4*3^3*5^2, index 2, d_K -2700. */
    { const long v[] = {-20,0,0,1}; const long pr[] = {2,3,5}; const long ex[] = {4,3,2};
      check_round2(v, 3, pr, ex, 3, 2, -2700, "Q(20^1/3)"); }
    /* Q(cbrt10): disc -2700 = -2^2*3^3*5^2, index 3, d_K -300. */
    { const long v[] = {-10,0,0,1}; const long pr[] = {2,3,5}; const long ex[] = {2,3,2};
      check_round2(v, 3, pr, ex, 3, 3, -300, "Q(10^1/3)"); }
    /* Q(cbrt17): disc(x^3-17) = -27*289 = -3^3*17^2, index 3, d_K -867. */
    { const long v[] = {-17,0,0,1}; const long pr[] = {3,17}; const long ex[] = {3,2};
      check_round2(v, 3, pr, ex, 2, 3, -867, "Q(17^1/3)"); }
    /* Monogenic sanity: Q(cbrt2), disc -108, index 1, d_K -108. */
    { const long v[] = {-2,0,0,1}; const long pr[] = {2,3}; const long ex[] = {2,3};
      check_round2(v, 3, pr, ex, 2, 1, -108, "Q(2^1/3) monogenic"); }
    /* Quartics -- multi-round maximization at one prime. */
    /* x^4-5: disc -32000 = -2^8*5^3, index 4, d_K -2000. */
    { const long v[] = {-5,0,0,0,1}; const long pr[] = {2,5}; const long ex[] = {8,3};
      check_round2(v, 4, pr, ex, 2, 4, -2000, "Q(5^1/4)"); }
    /* x^4-12: disc -442368 = -2^14*3^3, index 16 (four p=2 enlargements), d_K -1728. */
    { const long v[] = {-12,0,0,0,1}; const long pr[] = {2,3}; const long ex[] = {14,3};
      check_round2(v, 4, pr, ex, 2, 16, -1728, "Q(12^1/4)"); }
    printf("  round2: OK\n");
}
#endif

int main(void) {
    printf("Running test: numberfield (Gate 1) + nfunits (Gate 2)\n");
    test_signature();
    test_dedekind();
    test_create();
#ifdef USE_FLINT
    test_round2();
#endif
    test_units();
    printf("numberfield: ALL PASS\n");
    return 0;
}
