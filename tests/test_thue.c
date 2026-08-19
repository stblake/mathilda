/*
 * test_thue.c
 * -----------
 * Validates the Thue-equation reconstruction engine (src/solvethue.c):
 * building the number field of the binary form, enumerating unit
 * exponents, reconstructing (x, y) from the resulting unit's
 * coordinates, and verifying F(x,y) == m exactly.
 *
 * Uses thue_solve_binary_form_bounded (explicit exponent bound) since the
 * automatic Baker/de-Weger bound is not yet wired; the target solution
 * sets are the known-complete literature values, so an adequate bound
 * must reproduce them exactly.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "solvethue.h"
#include "numberfield.h"     /* nf_available */

/* form[] : a0..an (coeff of x^(n-j) y^j).  Check the solution set equals
 * the expected pairs (order-independent). */
static void check(const long* form_l, int n, long m_l, int bound,
                  const long expected[][2], int nexp, const char* name) {
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) mpz_init_set_si(form[i], form_l[i]);
    mpz_t m; mpz_init_set_si(m, m_l);

    ThueSol* sols = NULL;
    int cnt = thue_solve_binary_form_bounded((const mpz_t*)form, n, m, bound, &sols);
    printf("  %-22s: %d solutions\n", name, cnt);
    assert(cnt == nexp);

    /* every expected pair is present */
    for (int e = 0; e < nexp; e++) {
        int found = 0;
        for (int i = 0; i < cnt; i++)
            if (mpz_cmp_si(sols[i].x, expected[e][0]) == 0 &&
                mpz_cmp_si(sols[i].y, expected[e][1]) == 0) { found = 1; break; }
        if (!found) fprintf(stderr, "  MISSING (%ld,%ld) in %s\n", expected[e][0], expected[e][1], name);
        assert(found);
    }

    thue_sols_free(sols, cnt);
    for (int i = 0; i <= n; i++) mpz_clear(form[i]);
    free(form); mpz_clear(m);
}

/* An input that should DECLINE (out of the monic/|m|=1 scope): cnt == -1. */
static void check_decline(const long* form_l, int n, long m_l, int bound, const char* name) {
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) mpz_init_set_si(form[i], form_l[i]);
    mpz_t m; mpz_init_set_si(m, m_l);
    ThueSol* sols = NULL;
    int cnt = thue_solve_binary_form_bounded((const mpz_t*)form, n, m, bound, &sols);
    printf("  %-22s: decline=%d\n", name, cnt);
    assert(cnt == -1);
    (void)sols;
    for (int i = 0; i <= n; i++) mpz_clear(form[i]);
    free(form); mpz_clear(m);
}

/* Rigorous path: the automatic Baker/de-Weger bound (no explicit bound).
 * expected_count >= 0 asserts the complete set size; expect_decline asserts
 * a safe decline (-1). */
static void check_rigorous(const long* form_l, int n, long m_l,
                           int expected_count, int expect_decline, const char* name) {
    mpz_t* form = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (int i = 0; i <= n; i++) mpz_init_set_si(form[i], form_l[i]);
    mpz_t m; mpz_init_set_si(m, m_l);
    ThueSol* sols = NULL;
    int cnt = thue_solve_binary_form((const mpz_t*)form, n, m, &sols);
    printf("  %-28s: cnt=%d\n", name, cnt);
    if (expect_decline) assert(cnt == -1);
    else assert(cnt == expected_count);
    if (cnt > 0) thue_sols_free(sols, cnt);
    for (int i = 0; i <= n; i++) mpz_clear(form[i]);
    free(form); mpz_clear(m);
}

static void test_rigorous(void) {
    if (!nf_available()) { printf("  rigorous: SKIP (no FLINT)\n"); return; }
    /* The three target equations, solved via the automatic bound. */
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 1,  2, 0, "x^3-2y^3==1 (auto)"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, -1, 2, 0, "x^3-2y^3==-1 (auto)"); }
    { const long f[] = {1,0,0,-7};   check_rigorous(f, 3, 1,  2, 0, "x^3-7y^3==1 (auto)"); }
    { const long f[] = {1,0,-3,1};   check_rigorous(f, 3, 1,  6, 0, "x^3-3xy^2+y^3==1 (auto)"); }
    { const long f[] = {1,1,-2,-1};  check_rigorous(f, 3, 1,  9, 0, "Thomas cubic (auto)"); }
    { const long f[] = {1,0,0,0,-2}; check_rigorous(f, 4, -1, 4, 0, "x^4-2y^4==-1 (auto)"); }
    /* Large-regulator complex cubics: the coefficient box cannot reach the
     * fundamental unit, so these DECLINED until Voronoi's chain of minima (M1).
     * Complete set is {(1,0)} (== PARI/GP thue).  Q(cbrt41) has a 24-digit unit
     * (R=56.3), Q(cbrt97) ~10 digits (R=49.5): the stress cases. */
    { const long f[] = {1,0,0,-15};  check_rigorous(f, 3, 1,  1, 0, "x^3-15y^3==1 (Voronoi)"); }
    { const long f[] = {1,0,0,-42};  check_rigorous(f, 3, 1,  1, 0, "x^3-42y^3==1 (Voronoi)"); }
    { const long f[] = {1,0,0,-97};  check_rigorous(f, 3, 1,  1, 0, "x^3-97y^3==1 (Voronoi)"); }
    { const long f[] = {1,0,0,-41};  check_rigorous(f, 3, 1,  1, 0, "x^3-41y^3==1 (Voronoi R=56)"); }
    { const long f[] = {1,0,0,-41};  check_rigorous(f, 3, -1, 1, 0, "x^3-41y^3==-1 (Voronoi)"); }
    /* General |m| != 1 (M2): beta = mu * unit, mu over bounded-norm reps.
     * Complete sets cross-checked against PARI/GP thue().  Rank-1 complex cubic. */
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 2,   1, 0, "x^3-2y^3==2  (M2, {(0,-1)})"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 3,   2, 0, "x^3-2y^3==3  (M2, 2 sols)"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 10,  2, 0, "x^3-2y^3==10 (M2, 2 sols)"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 4,   0, 0, "x^3-2y^3==4  (M2, proven {})"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 5,   0, 0, "x^3-2y^3==5  (M2, proven {})"); }
    { const long f[] = {1,0,0,-2};   check_rigorous(f, 3, 73,  0, 0, "x^3-2y^3==73 (M2, no norm-73, {})"); }
    { const long f[] = {1,0,0,-7};   check_rigorous(f, 3, 6,   1, 0, "x^3-7y^3==6  (M2, {(-1,-1)})"); }
    /* Non-monogenic fields now solve via Round 2 (M3), cross-checked vs PARI:
     * Q(cbrt17) index 3 -> {(1,0),(18,7)}, Q(cbrt10) index 3 -> {(1,0)},
     * Q(cbrt20) index 2 -> {(1,0),(-19,-7)}, the Dedekind cubic index 2 -> {(1,0)}. */
    { const long f[] = {1,0,0,-17};  check_rigorous(f, 3, 1,  2, 0, "x^3-17y^3==1 (M3 non-monogenic)"); }
    { const long f[] = {1,0,0,-10};  check_rigorous(f, 3, 1,  1, 0, "x^3-10y^3==1 (M3 non-monogenic)"); }
    { const long f[] = {1,0,0,-20};  check_rigorous(f, 3, 1,  2, 0, "x^3-20y^3==1 (M3 non-monogenic)"); }
    { const long f[] = {1,-1,-2,-8}; check_rigorous(f, 3, 1,  1, 0, "Dedekind cubic (M3 non-monogenic)"); }
    /* Still out of scope -> DECLINE: rank-2 totally-real |m| != 1 (M2b). */
    { const long f[] = {1,-3,0,1};   check_rigorous(f, 3, 2,  0, 1, "x^3-3xy^2+y^3==2 (rank-2 |m|!=1: decline)"); }
    printf("  rigorous: OK\n");
}

int main(void) {
    printf("Running test: thue (reconstruction engine)\n");
    if (!nf_available()) { printf("thue: SKIP (no FLINT)\n"); return 0; }

    /* x^3 - 2 y^3 == 1  ->  {(1,0),(-1,-1)} */
    { const long f[] = {1,0,0,-2}; const long e[][2] = {{1,0},{-1,-1}};
      check(f, 3, 1, 20, e, 2, "x^3-2y^3==1"); }
    /* x^3 - 2 y^3 == -1  ->  {(-1,0),(1,1)} */
    { const long f[] = {1,0,0,-2}; const long e[][2] = {{-1,0},{1,1}};
      check(f, 3, -1, 20, e, 2, "x^3-2y^3==-1"); }
    /* x^3 - 7 y^3 == 1  ->  {(1,0),(2,1)} */
    { const long f[] = {1,0,0,-7}; const long e[][2] = {{1,0},{2,1}};
      check(f, 3, 1, 20, e, 2, "x^3-7y^3==1"); }
    /* x^3 - 3 x y^2 + y^3 == 1 (cyclic cubic, conductor 9) -> 6 solutions */
    { const long f[] = {1,0,-3,1};
      const long e[][2] = {{1,0},{0,1},{-1,-1},{2,-1},{1,3},{-3,-2}};
      check(f, 3, 1, 20, e, 6, "x^3-3xy^2+y^3==1"); }
    /* x^4 - 2 y^4 == -1  ->  {(1,1),(1,-1),(-1,1),(-1,-1)} */
    { const long f[] = {1,0,0,0,-2};
      const long e[][2] = {{1,1},{1,-1},{-1,1},{-1,-1}};
      check(f, 4, -1, 20, e, 4, "x^4-2y^4==-1"); }

    /* |m| != 1 now solves via mu-enumeration (M2): x^3-2y^3==4 -> {} (proven). */
    { const long f[] = {1,0,0,-2}; const long e[][2] = {{0,0}};
      check(f, 3, 4, 20, e, 0, "x^3-2y^3==4 (M2, {})"); }
    /* Still out of scope -> DECLINE: |a0| != 1. */
    { const long f[] = {2,0,0,-3}; check_decline(f, 3, 1, 20, "2x^3-3y^3==1 (|a0|!=1)"); }

    test_rigorous();

    printf("thue: ALL PASS\n");
    return 0;
}
