/* gmp_compat.h -- portability shims for GMP functions that are absent from
 * older system GMP releases.
 *
 * mpz_prevprime was only added in GMP 6.3.0 (2023-07-30).  Current long-term
 * Linux distributions still ship GMP 6.2.x (e.g. Debian 12, Ubuntu 22.04 LTS),
 * where the symbol is undefined and the tree fails to link.  Callers use the
 * uniform wrapper mathilda_mpz_prevprime(), which forwards to the native
 * routine on GMP >= 6.3.0 and runs an identical fallback on anything older, so
 * no call site needs a version guard of its own.
 *
 * This header is include-order independent: it pulls in <gmp.h> itself. */
#ifndef MATHILDA_GMP_COMPAT_H
#define MATHILDA_GMP_COMPAT_H

#include <gmp.h>

/* Use the native mpz_prevprime when the system GMP is new enough, unless a
 * build explicitly forces the fallback (the unit test does this to exercise the
 * fallback on machines that have 6.3.0+; collapse the version triple into one
 * comparable integer: 6.3.0 -> 60300, 6.2.1 -> 60201). */
#if !defined(MATHILDA_FORCE_MPZ_PREVPRIME_FALLBACK) &&                        \
    defined(__GNU_MP_VERSION) &&                                             \
    (__GNU_MP_VERSION * 10000 + __GNU_MP_VERSION_MINOR * 100 +               \
     __GNU_MP_VERSION_PATCHLEVEL) >= 60300
#  define MATHILDA_HAVE_NATIVE_MPZ_PREVPRIME 1
#else
#  define MATHILDA_HAVE_NATIVE_MPZ_PREVPRIME 0
#endif

/* Set rop to the greatest prime strictly less than op.  Return value matches
 * GMP 6.3.0's native mpz_prevprime exactly:
 *   2  if rop is certainly prime,
 *   1  if rop is probably prime,
 *   0  if op <= 2 -- no prime precedes it, and rop is left equal to op.
 * Safe when rop and op alias (the callers pass the same mpz_t for both): op is
 * copied before rop is written. */
static inline int mathilda_mpz_prevprime(mpz_t rop, const mpz_t op) {
#if MATHILDA_HAVE_NATIVE_MPZ_PREVPRIME
    return mpz_prevprime(rop, op);
#else
    if (mpz_cmp_ui(op, 2) <= 0) {
        mpz_set(rop, op);
        return 0;
    }
    if (mpz_cmp_ui(op, 3) == 0) {   /* the prime before 3 is 2 */
        mpz_set_ui(rop, 2);
        return 2;
    }
    mpz_t cand;
    mpz_init_set(cand, op);
    mpz_sub_ui(cand, cand, 1);
    if (mpz_even_p(cand)) mpz_sub_ui(cand, cand, 1);   /* land on an odd */
    while (mpz_cmp_ui(cand, 2) > 0) {
        int prime = mpz_probab_prime_p(cand, 25);
        if (prime) {
            mpz_set(rop, cand);
            mpz_clear(cand);
            return prime == 2 ? 2 : 1;
        }
        mpz_sub_ui(cand, cand, 2);
    }
    /* Unreachable for op >= 4 (3 is odd and prime), but a safe backstop. */
    mpz_clear(cand);
    mpz_set_ui(rop, 2);
    return 2;
#endif
}

#endif  /* MATHILDA_GMP_COMPAT_H */
