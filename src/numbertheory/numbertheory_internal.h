#ifndef NUMBERTHEORY_INTERNAL_H
#define NUMBERTHEORY_INTERNAL_H

/* Cross-cutting helpers shared between the per-builtin translation units of
 * src/numbertheory/.  Everything declared here is internal to the
 * number-theory subsystem and must NOT leak into any other compilation unit
 * (mirrors the eigen_internal.h / ludecomp_internal.h pattern in
 * src/linalg/).  Each helper is defined in exactly one .c file as noted. */

#include "numbertheory.h"   /* Expr */
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <gmp.h>

/* Maximum distinct prime factors collected by pr_collect_distinct_primes;
 * the caller-supplied mpz_t arrays must have this capacity.  Shared by
 * primitiveroot.c and multiplicativeorder.c. */
#define PR_MAX_DISTINCT_PRIMES 256

/* --- nt_util.c : generic GCD/LCM/rational helpers ---------------------- */
Expr*   single_arg_abs_or_copy(Expr* arg);
void    rational_like_to_mpz_pair(const Expr* e, mpz_t num, mpz_t den);
bool    rational_like_needs_bigint(const Expr* e);
Expr*   mpz_pair_to_rational_expr(const mpz_t num_in, const mpz_t den_in);
int64_t lcm_checked(int64_t a, int64_t b, bool *overflow);

/* --- primitiveroot.c (shared with multiplicativeorder.c) -------------- */
/* Gather the distinct prime factors of m_in into out[0..*count); out must
 * have PR_MAX_DISTINCT_PRIMES capacity and *count must be 0 on entry. */
bool    pr_collect_distinct_primes(const mpz_t m_in, mpz_t* out, size_t* count);

/* --- nt_gaussian.c : Gaussian-integer & divisor enumeration ----------- */
bool    is_gaussian_integer(Expr* e);
bool    df_factor_mpz(const mpz_t m, mpz_t** out_primes,
                      unsigned long** out_exps, size_t* out_n);
bool    df_to_gaussian(Expr* e, mpz_t a, mpz_t b);
bool    df_gaussian_prime_factor(const mpz_t a_in, const mpz_t b_in,
                                 mpz_t** out_gu, mpz_t** out_gv,
                                 unsigned long** out_ge, size_t* out_gn);
void    df_normalize_quadrant(mpz_t re, mpz_t im);
Expr*   divisors_ordinary(const mpz_t nabs);
Expr*   divisors_gaussian(const mpz_t a_in, const mpz_t b_in);

#endif
