/* moebiusmu.c -- MoebiusMu[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "internal.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

/* ---- MoebiusMu ---------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic for MoebiusMu and return NULL so the
 * call is left unevaluated.  WL expects exactly one argument. */
static Expr* moebiusmu_emit_argx(size_t argc) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "MoebiusMu::argx: MoebiusMu called with %zu argument%s; "
                "1 argument is expected.\n",
                argc, argc == 1 ? "" : "s");
    }
    return NULL;
}

/* Reduce a prime-power factorisation to the Mobius value: 0 when any exponent is
 * >= 2 (a squared prime factor), otherwise (-1)^(number of distinct primes).
 * The empty product (n a unit) gives 1. */
static Expr* moebiusmu_from_exps(const unsigned long* exps, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (exps[i] >= 2) return expr_new_integer(0);
    return expr_new_integer((n % 2 == 0) ? 1 : -1);
}

/* MoebiusMu[n] gives the Mobius function mu(n).  For n = u * Prod p_i^e_i with u
 * a unit and p_i primes, mu(n) is 0 if any e_i >= 2, otherwise (-1)^m where m is
 * the number of distinct primes; mu(1) = 1.  The sign of n is ignored, matching
 * mu(-n) = mu(n).  A non-real Gaussian-integer argument Complex[a, b] is handled
 * over Z[i] by factoring into Gaussian primes (the unit factor does not count).
 * MoebiusMu is Listable, so list arguments are threaded by the evaluator first.
 * Non-integer or zero n is left unevaluated; a wrong argument count emits
 * MoebiusMu::argx. */
Expr* builtin_moebiusmu(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return moebiusmu_emit_argx(argc);
    Expr* n = res->data.function.args[0];

    /* Gaussian-integer path, auto-detected from non-real Complex input. */
    if (is_gaussian_integer(n)) {
        mpz_t a, b;
        if (!df_to_gaussian(n, a, b)) return NULL;       /* a, b init'd on success */
        mpz_t* gu;
        mpz_t* gv;
        unsigned long* ge;
        size_t gn;
        bool fok = df_gaussian_prime_factor(a, b, &gu, &gv, &ge, &gn);
        mpz_clears(a, b, NULL);
        if (!fok) return NULL;                           /* z == 0 / factoring fail */
        Expr* out = moebiusmu_from_exps(ge, gn);
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return out;
    }

    /* Ordinary integer path (covers EXPR_INTEGER and EXPR_BIGINT). */
    if (!expr_is_integer_like(n)) return NULL;           /* symbolic n stays put */
    mpz_t v;
    expr_to_mpz(n, v);                                   /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }  /* MoebiusMu[0] undefined */
    mpz_abs(v, v);                                       /* sign of n is ignored */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(v, &primes, &exps, &np);
    mpz_clear(v);
    if (!fok) return NULL;

    Expr* out = moebiusmu_from_exps(exps, np);
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}
