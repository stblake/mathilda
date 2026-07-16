/* primenu.c -- PrimeNu[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout.
 *
 * PrimeNu[n] = nu(n), the number of DISTINCT prime factors of n.  It is the
 * additive companion to PrimeOmega (which counts prime factors with
 * multiplicity): for n = u p_1^k_1 ... p_m^k_m with u a unit and p_i distinct
 * primes, PrimeNu[n] returns m.  nu and Omega coincide exactly when n is
 * square-free.  PrimeNu shares all factoring machinery and argument handling
 * with PrimeOmega/LiouvilleLambda; it simply returns the count of factors
 * rather than the sum of the exponents. */

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

/* ---- PrimeNu ------------------------------------------------------------ */

/* Emit the wrong-argument-count diagnostic for PrimeNu and return NULL so the
 * call is left unevaluated. */
static Expr* primenu_emit_argt(size_t argc) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "PrimeNu::argt: PrimeNu called with %zu "
                "argument%s; 1 or 2 arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
    }
    return NULL;
}

/* nu(n) is the number of DISTINCT prime factors -- the count of factors,
 * independent of their exponents. */
static Expr* primenu_from_count(size_t n) {
    return expr_new_integer((long)n);
}

/* PrimeNu[n] gives the number of distinct prime factors of n, nu(n).  With
 * GaussianIntegers -> True (or a non-real Gaussian-integer n) the count runs
 * over the distinct Gaussian prime factors of n.  PrimeNu is Listable, so list
 * arguments are threaded by the evaluator first.  PrimeNu[1] (and PrimeNu[-1])
 * is 0; a non-integer or zero n is left unevaluated; a wrong argument count
 * emits PrimeNu::argt. */
Expr* builtin_primenu(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc < 1 || argc > 2) return primenu_emit_argt(argc);

    /* Separate the positional argument (n) from a GaussianIntegers option. */
    bool gaussian = false, gaussian_set = false;
    Expr* posargs[2] = { NULL, NULL };
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol.name == SYM_True || val->data.symbol.name == SYM_False)) {
                gaussian = (val->data.symbol.name == SYM_True);
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
        Expr* out = primenu_from_count(gn);
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return out;
    }

    /* Ordinary integer path (covers EXPR_INTEGER and EXPR_BIGINT). */
    if (!expr_is_integer_like(n)) return NULL;            /* symbolic n stays put */
    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* PrimeNu[0] undefined */
    mpz_abs(v, v);                                        /* nu(-n) = nu(n) */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(v, &primes, &exps, &np);
    mpz_clear(v);
    if (!fok) return NULL;

    Expr* out = primenu_from_count(np);
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}
