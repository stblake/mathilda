/* Tests for src/gmp_compat.h -- specifically the mathilda_mpz_prevprime
 * fallback used when the system GMP predates 6.3.0 (which first shipped
 * mpz_prevprime).  This is exactly the path a user on Debian 12 / Ubuntu 22.04
 * (GMP 6.2.x) hits, and it is dead code on every machine that has 6.3.0+.
 *
 * To exercise the fallback regardless of the build machine's GMP, this TU
 * forces it on with MATHILDA_FORCE_MPZ_PREVPRIME_FALLBACK before including the
 * header.  The raw native mpz_prevprime (a plain gmp.h macro on 6.3.0+) stays
 * reachable under its own name, so where it exists we cross-check the fallback
 * against it directly. */
#define MATHILDA_FORCE_MPZ_PREVPRIME_FALLBACK 1
#include "../src/gmp_compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GMP >= 6.3.0 -> native mpz_prevprime available to cross-check against.
 * Resolve to a plain 0/1 object macro here (evaluating defined() at directive
 * time), so later `#if GMP_AT_LEAST_630` never expands `defined` from a macro. */
#if defined(__GNU_MP_VERSION) &&                                              \
    (__GNU_MP_VERSION * 10000 + __GNU_MP_VERSION_MINOR * 100 +               \
     __GNU_MP_VERSION_PATCHLEVEL) >= 60300
#  define GMP_AT_LEAST_630 1
#else
#  define GMP_AT_LEAST_630 0
#endif

static int failures = 0;

static void fail(const char* msg) {
    fprintf(stderr, "FAIL: %s\n", msg);
    failures++;
}

/* 1. Known-answer table -- independent of GMP version. Triples are
 * {op, expected greatest prime < op, expected return code}. */
static void test_known_answers(void) {
    struct { long op; const char* prev; int rc; } cases[] = {
        {0, "0", 0}, {1, "1", 0}, {2, "2", 0},   /* no prime precedes; rop=op */
        {3, "2", 2}, {4, "3", 2}, {5, "3", 2}, {6, "5", 2},
        {7, "5", 2}, {8, "7", 2}, {9, "7", 2},
        {100, "97", 2}, {101, "97", 2}, {1000, "997", 2}, {7920, "7919", 2},
    };
    mpz_t op, r;
    mpz_inits(op, r, NULL);
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        mpz_set_si(op, cases[i].op);
        int rc = mathilda_mpz_prevprime(r, op);
        char* got = mpz_get_str(NULL, 10, r);
        if (rc != cases[i].rc || strcmp(got, cases[i].prev) != 0) {
            fprintf(stderr, "  op=%ld: got rc=%d r=%s, want rc=%d r=%s\n",
                    cases[i].op, rc, got, cases[i].rc, cases[i].prev);
            fail("known-answer mismatch");
        }
        free(got);
    }
    mpz_clears(op, r, NULL);
}

/* 2. Self-consistency: for every op in [3, N], the result must be prime, be
 * strictly less than op, and have no prime strictly between it and op. */
static void test_self_consistent(void) {
    mpz_t op, r, t;
    mpz_inits(op, r, t, NULL);
    for (long v = 3; v <= 3000; v++) {
        mpz_set_si(op, v);
        int rc = mathilda_mpz_prevprime(r, op);
        if (rc == 0) { fail("unexpected rc=0 for op>=3"); break; }
        if (mpz_cmp(r, op) >= 0) { fail("result not < op"); break; }
        if (!mpz_probab_prime_p(r, 25)) { fail("result not prime"); break; }
        for (mpz_add_ui(t, r, 1); mpz_cmp(t, op) < 0; mpz_add_ui(t, t, 1)) {
            if (mpz_probab_prime_p(t, 25)) { fail("skipped a prime"); break; }
        }
    }
    mpz_clears(op, r, t, NULL);
}

/* 3. Aliasing: callers pass the same mpz_t for rop and op (in-place). */
static void test_in_place(void) {
    mpz_t x;
    mpz_init_set_ui(x, 100);
    mathilda_mpz_prevprime(x, x);
    if (mpz_cmp_ui(x, 97) != 0) fail("in-place prevprime(100) != 97");
    mpz_clear(x);
}

/* 4. A value well beyond int64, and a direct native cross-check where the
 * native routine exists (i.e. on the CI machines, which run 6.3.0+). */
static void test_large_and_native(void) {
    mpz_t op, r;
    mpz_inits(op, r, NULL);
    mpz_set_str(op, "1000000000000000000000000000000", 10);   /* 10^30 */
    mathilda_mpz_prevprime(r, op);
    char* got = mpz_get_str(NULL, 10, r);
    if (strcmp(got, "999999999999999999999999999989") != 0)
        fail("prevprime(10^30) wrong");
    free(got);

#if GMP_AT_LEAST_630
    /* Fallback must agree with native across a wide range and a large value. */
    mpz_t rn;
    mpz_init(rn);
    for (long v = 0; v <= 4000; v++) {
        mpz_set_si(op, v);
        int rc_fb = mathilda_mpz_prevprime(r, op);
        int rc_nat = mpz_prevprime(rn, op);   /* native gmp.h macro */
        if ((rc_fb != 0) != (rc_nat != 0) ||
            (rc_nat != 0 && mpz_cmp(r, rn) != 0)) {
            fprintf(stderr, "  op=%ld fb(rc=%d) native(rc=%d) differ\n",
                    v, rc_fb, rc_nat);
            fail("fallback disagrees with native mpz_prevprime");
            break;
        }
    }
    mpz_set_str(op, "340282366920938463463374607431768211456", 10);  /* 2^128 */
    mathilda_mpz_prevprime(r, op);
    mpz_prevprime(rn, op);
    if (mpz_cmp(r, rn) != 0) fail("fallback != native at 2^128");
    mpz_clear(rn);
#endif
    mpz_clears(op, r, NULL);
}

int main(void) {
    printf("Running test: test_known_answers\n");   test_known_answers();
    printf("Running test: test_self_consistent\n"); test_self_consistent();
    printf("Running test: test_in_place\n");        test_in_place();
    printf("Running test: test_large_and_native\n");test_large_and_native();
#if GMP_AT_LEAST_630
    printf("(cross-checked fallback against native mpz_prevprime)\n");
#else
    printf("(native mpz_prevprime absent -- this GMP is the real fallback case)\n");
#endif
    if (failures) { printf("%d gmp_compat failure(s)\n", failures); return 1; }
    printf("All gmp_compat tests passed!\n");
    return 0;
}
