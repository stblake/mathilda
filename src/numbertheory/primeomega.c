/* primeomega.c -- PrimeOmega[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout.
 *
 * PrimeOmega[n] = Omega(n), the number of prime factors of n counted with
 * multiplicity (the sum of the exponents in the prime factorization).  This is
 * the quantity LiouvilleLambda computes internally before taking (-1)^Omega, so
 * the two share the same factoring machinery and argument handling; PrimeOmega
 * simply returns Omega itself. */

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

/* ---- PrimeOmega --------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic for PrimeOmega and return NULL so
 * the call is left unevaluated. */
static Expr* primeomega_emit_argt(size_t argc) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "PrimeOmega::argt: PrimeOmega called with %zu "
                "argument%s; 1 or 2 arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
    }
    return NULL;
}

/* Omega(n) is the total number of prime factors counted with multiplicity,
 * i.e. the sum of the exponents. */
static Expr* primeomega_from_exps(const unsigned long* exps, size_t n) {
    unsigned long omega = 0;
    for (size_t i = 0; i < n; i++) omega += exps[i];
    return expr_new_integer((long)omega);
}

/* PrimeOmega[n] gives the number of prime factors of n counted with
 * multiplicity, Omega(n).  With GaussianIntegers -> True (or a non-real
 * Gaussian-integer n) the count runs over the Gaussian prime factors of n.
 * PrimeOmega is Listable, so list arguments are threaded by the evaluator
 * first.  PrimeOmega[1] (and PrimeOmega[-1]) is 0; a non-integer or zero n is
 * left unevaluated; a wrong argument count emits PrimeOmega::argt. */
Expr* builtin_primeomega(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc < 1 || argc > 2) return primeomega_emit_argt(argc);

    /* Separate the positional argument (n) from a GaussianIntegers option. */
    bool gaussian = false, gaussian_set = false;
    Expr* posargs[2] = { NULL, NULL };
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol == SYM_True || val->data.symbol == SYM_False)) {
                gaussian = (val->data.symbol == SYM_True);
                gaussian_set = true;
                continue;
            }
            return NULL;  /* malformed / unknown option: leave unevaluated */
        }
        if (npos < 2) posargs[npos] = a;
        npos++;
    }

    if (npos != 1) return NULL;   /* not exactly one positional: leave unevaluated */
    Expr* n = posargs[0];

    /* Explicit option wins; otherwise auto-enable Z[i] for non-real input. */
    bool use_gaussian = gaussian_set ? gaussian : is_gaussian_integer(n);

    if (use_gaussian) {
        mpz_t a, b;
        if (!df_to_gaussian(n, a, b)) return NULL;        /* a, b init'd on success */
        mpz_t* gu;
        mpz_t* gv;
        unsigned long* ge;
        size_t gn;
        bool fok = df_gaussian_prime_factor(a, b, &gu, &gv, &ge, &gn);
        mpz_clears(a, b, NULL);
        if (!fok) return NULL;                            /* n == 0 / factoring fail */
        Expr* out = primeomega_from_exps(ge, gn);
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return out;
    }

    /* Ordinary integer path (covers EXPR_INTEGER and EXPR_BIGINT). */
    if (!expr_is_integer_like(n)) return NULL;            /* symbolic n stays put */
    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* PrimeOmega[0] undefined */
    mpz_abs(v, v);                                        /* Omega(-n) = Omega(n) */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(v, &primes, &exps, &np);
    mpz_clear(v);
    if (!fok) return NULL;

    Expr* out = primeomega_from_exps(exps, np);
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}
