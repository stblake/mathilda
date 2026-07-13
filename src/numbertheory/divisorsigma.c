/* divisorsigma.c -- DivisorSigma[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* ---- DivisorSigma ------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic for DivisorSigma and return NULL so
 * the call is left unevaluated. */
static Expr* divisorsigma_emit_argrx(size_t npos) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "DivisorSigma::argrx: DivisorSigma called with %zu argument%s; "
                "2 arguments are expected.\n",
                npos, npos == 1 ? "" : "s");
    }
    return NULL;
}

/* Prime atom for the multiplicative sigma_k formula: a rational integer (when
 * v == 0) or the Gaussian integer Complex[u, v]. */
static Expr* ds_prime_atom(const mpz_t u, const mpz_t v) {
    if (mpz_sgn(v) == 0)
        return expr_bigint_normalize(expr_new_bigint_from_mpz(u));
    Expr* cargs[2] = {
        expr_bigint_normalize(expr_new_bigint_from_mpz(u)),
        expr_bigint_normalize(expr_new_bigint_from_mpz(v))
    };
    return expr_new_function(expr_new_symbol(SYM_Complex), cargs, 2);
}

/* Build the Expr for one multiplicative factor of sigma_k(n):
 *   (q^((e+1) k) - 1) / (q^k - 1).
 * Ownership of q is taken (consumed); k is copied. */
static Expr* ds_build_factor(Expr* q, unsigned long e, const Expr* k) {
    /* numerator: q^((e+1) k) - 1 */
    Expr* nexp_args[2] = { expr_new_integer((int64_t)(e + 1)), expr_copy((Expr*)k) };
    Expr* nexp = expr_new_function(expr_new_symbol(SYM_Times), nexp_args, 2);
    Expr* npow_args[2] = { expr_copy(q), nexp };
    Expr* npow = expr_new_function(expr_new_symbol(SYM_Power), npow_args, 2);
    Expr* num_args[2] = { npow, expr_new_integer(-1) };
    Expr* num = expr_new_function(expr_new_symbol(SYM_Plus), num_args, 2);

    /* denominator: q^k - 1 */
    Expr* dpow_args[2] = { q, expr_copy((Expr*)k) };   /* consumes q */
    Expr* dpow = expr_new_function(expr_new_symbol(SYM_Power), dpow_args, 2);
    Expr* den_args[2] = { dpow, expr_new_integer(-1) };
    Expr* den = expr_new_function(expr_new_symbol(SYM_Plus), den_args, 2);

    /* (q^k - 1)^(-1) */
    Expr* inv_args[2] = { den, expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);

    Expr* fac_args[2] = { num, inv };
    return expr_new_function(expr_new_symbol(SYM_Times), fac_args, 2);
}

/* Assemble and evaluate the multiplicative product Prod_i factor(atoms[i]).
 * atoms[i] are consumed; an empty product (n == 0) is the integer 1. */
static Expr* ds_assemble(Expr** atoms, const unsigned long* exps, size_t n,
                         const Expr* k) {
    if (n == 0) return expr_new_integer(1);
    Expr** facs = (Expr**)malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        facs[i] = ds_build_factor(atoms[i], exps[i], k);   /* consumes atoms[i] */
    Expr* prod = (n == 1) ? facs[0]
                          : expr_new_function(expr_new_symbol(SYM_Times), facs, n);
    free(facs);
    return eval_and_free(prod);
}

/* The number of divisors, Prod_i (exps[i] + 1), as an Integer/BigInt.  This is
 * sigma_0(n) — the k == 0 case where the multiplicative formula degenerates. */
static Expr* ds_divisor_count(const unsigned long* exps, size_t n) {
    mpz_t c, t;
    mpz_init_set_ui(c, 1);
    mpz_init(t);
    for (size_t i = 0; i < n; i++) {
        mpz_set_ui(t, (unsigned long)exps[i] + 1);
        mpz_mul(c, c, t);
    }
    mpz_clear(t);
    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(c));
    mpz_clear(c);
    return out;
}

/* DivisorSigma[k, n] gives the divisor function sigma_k(n), the sum of the k-th
 * powers of the divisors of n.  It is computed from the multiplicative formula
 *   sigma_k(n) = Prod_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)   for n = Prod p_i^e_i,
 * built as an expression and evaluated, so a single path serves integer k
 * (exact integers / rationals), and rational or symbolic k (radical / symbolic
 * forms).  k == 0 returns the divisor count.  With GaussianIntegers -> True (or
 * a non-real Gaussian-integer n) the product runs over the first-quadrant
 * associates of the Gaussian prime factors.  DivisorSigma is Listable, so list
 * arguments are threaded by the evaluator first.  Non-integer or zero n is left
 * unevaluated; a wrong argument count emits DivisorSigma::argrx. */
Expr* builtin_divisorsigma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Separate the positional arguments (k, n) from a GaussianIntegers option. */
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

    if (npos != 2) return divisorsigma_emit_argrx(npos);
    Expr* k = posargs[0];
    Expr* n = posargs[1];
    bool k_is_zero = (k->type == EXPR_INTEGER && k->data.integer == 0);

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

        Expr* out;
        if (k_is_zero) {
            out = ds_divisor_count(ge, gn);
        } else {
            Expr** atoms = (gn > 0) ? (Expr**)malloc(gn * sizeof(Expr*)) : NULL;
            for (size_t i = 0; i < gn; i++) {
                df_normalize_quadrant(gu[i], gv[i]);      /* first-quadrant associate */
                atoms[i] = ds_prime_atom(gu[i], gv[i]);
            }
            out = ds_assemble(atoms, ge, gn, k);
            free(atoms);
        }
        for (size_t i = 0; i < gn; i++) mpz_clears(gu[i], gv[i], NULL);
        free(gu); free(gv); free(ge);
        return out;
    }

    /* Ordinary integer path. */
    if (!expr_is_integer_like(n)) return NULL;            /* symbolic n stays put */
    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* DivisorSigma[k, 0] */
    mpz_abs(v, v);                                        /* sign of n is ignored */

    mpz_t* primes;
    unsigned long* exps;
    size_t np;
    bool fok = df_factor_mpz(v, &primes, &exps, &np);
    mpz_clear(v);
    if (!fok) return NULL;

    Expr* out;
    if (k_is_zero) {
        out = ds_divisor_count(exps, np);
    } else {
        Expr** atoms = (np > 0) ? (Expr**)malloc(np * sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < np; i++)
            atoms[i] = expr_bigint_normalize(expr_new_bigint_from_mpz(primes[i]));
        out = ds_assemble(atoms, exps, np, k);
        free(atoms);
    }
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes);
    free(exps);
    return out;
}
