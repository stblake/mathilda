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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "ndinteger.h" /* ndint_divisor_sigma_i64: the shared machine helper */
#include "checked_int.h"
#include "numeric.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>
#include <math.h>
#include <limits.h>

/* ---- DivisorSigma ------------------------------------------------------- */

/* Emit the wrong-argument-count diagnostic and return NULL so the call is left
 * unevaluated. The tag is `argr` for exactly one argument and `argrx`
 * otherwise. Gated on the raw arithmetic probe mute, so an internal probe that
 * pokes at a malformed call is invisible to Check[] as well as silent. */
static Expr* divisorsigma_emit_argcount(size_t npos) {
    if (npos == 1)
        mth_message_gated(g_arith_warnings_muted, "DivisorSigma", "argr",
            "DivisorSigma called with 1 argument; 2 arguments are expected.");
    else
        mth_message_gated(g_arith_warnings_muted, "DivisorSigma", "argrx",
            "DivisorSigma called with %zu arguments; 2 arguments are expected.", npos);
    return NULL;
}

static bool ds_is_rule(const Expr* a) {
    return a->type == EXPR_FUNCTION
        && a->data.function.head->type == EXPR_SYMBOL
        && (a->data.function.head->data.symbol.name == SYM_Rule
            || a->data.function.head->data.symbol.name == SYM_RuleDelayed)
        && a->data.function.arg_count == 2;
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

/* ---- exact and machine fast paths ---------------------------------------- *
 * The multiplicative formula built as an expression (ds_assemble) is the one
 * path that serves EVERY k -- rational, symbolic, complex, arbitrary precision
 * -- and it is also the slow one: a Times/Power/Plus tree per call, then a full
 * evaluation of it (2x the target on Do[DivisorSigma[1, i], {i, 10^5}]). The
 * two kinds of k that dominate real workloads, an exact machine Integer and a
 * machine Real, are answered here directly instead. Each path returns NULL to
 * decline, and the caller then takes the expression path exactly as before. */

/* |n| at or below which a machine n is factored by the int64 trial division
 * the NDArray kernel uses (ndint_divisor_sigma_i64) rather than by
 * FactorInteger: at most 512 odd trial divisors, under a microsecond, with no
 * GMP allocation. Past it FactorInteger's ladder is the faster factoriser. */
#define DS_MACHINE_TRIAL_MAX ((int64_t)1 << 20)

/* sigma_|k|(n) as an exact mpz from n's factorisation, by the same Horner sum
 * as ndint_sigma_from_factors: prod_i (1 + q_i (1 + q_i (...))), q_i = p_i^|k|. */
static void ds_sigma_mpz(mpz_t out, unsigned long ak, mpz_t* primes,
                         const unsigned long* exps, size_t np) {
    mpz_t q, term;
    mpz_inits(q, term, NULL);
    mpz_set_ui(out, 1);
    for (size_t i = 0; i < np; i++) {
        mpz_pow_ui(q, primes[i], ak);
        mpz_set_ui(term, 1);
        for (unsigned long j = 0; j < exps[i]; j++) {
            mpz_mul(term, term, q);
            mpz_add_ui(term, term, 1);
        }
        mpz_mul(out, out, term);
    }
    mpz_clears(q, term, NULL);
}

/* DivisorSigma[k, n] for a machine Integer k and an Integer n != 0 (sign of n
 * ignored). k >= 0 is an Integer; k < 0 uses sigma_-k(n) = sigma_k(n) / n^k and
 * is a Rational in lowest terms. Declines (NULL) only when the factorisation
 * itself fails, or |k| is past what mpz_pow_ui's exponent can carry. */
static Expr* ds_integer_k(int64_t k, const Expr* n) {
    /* The common case: a small n and a k whose answer fits int64 -- the very
     * helper behind the packed kernel and Compile, so the surfaces agree. */
    if (k == INT64_MIN) return NULL;
    if (n->type == EXPR_INTEGER) {
        int64_t an;
        if (!ci_abs_i64(n->data.integer, &an) && an <= DS_MACHINE_TRIAL_MAX) {
            int64_t r, d;
            if (k >= 0) {
                if (ndint_divisor_sigma_i64(k, an, &r)) return expr_new_integer(r);
            } else if (ndint_divisor_sigma_i64(-k, an, &r) && !ci_powi_i64(an, -k, &d)) {
                return make_rational(r, d);         /* sigma_-k(n) = sigma_k(n)/n^k */
            }
        }
    }
    uint64_t ak = (uint64_t)(k < 0 ? -k : k);
    if (ak > (uint64_t)ULONG_MAX) return NULL;           /* LLP64: 32-bit ulong */

    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
    mpz_abs(v, v);
    mpz_t* primes; unsigned long* exps; size_t np;
    if (!df_factor_mpz(v, &primes, &exps, &np)) { mpz_clear(v); return NULL; }

    Expr* out;
    if (k == 0) {
        out = ds_divisor_count(exps, np);
    } else {
        mpz_t num;
        mpz_init(num);
        ds_sigma_mpz(num, (unsigned long)ak, primes, exps, np);
        if (k > 0) {
            out = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
        } else {
            mpz_t den;
            mpz_init(den);
            mpz_pow_ui(den, v, (unsigned long)ak);
            out = make_rational_mpz(num, den);
            mpz_clear(den);
        }
        mpz_clear(num);
    }
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes); free(exps);
    mpz_clear(v);
    return out;
}

/* DivisorSigma[k, n] for a machine Real k and an Integer n != 0, in double.
 *
 * Per prime the geometric sum 1 + q + ... + q^e (q = p^k) is summed by Horner
 * rather than as (q^(e+1) - 1)/(q - 1). Measured against a 40-digit reference
 * over n <= 3000 and k in {0.5, 1.5, 2., -1.5, 3.7}, Horner's worst case is
 * 5.8 ulp, where the closed form (and the expression path it replaces) reached
 * 31 ulp at k = 3.7 -- and the closed form is 0/0 at k = 0., where Horner gives
 * the divisor count. n = 1 is the empty product, the exact Integer 1
 * (DivisorSigma[2., {1, 2, 3}] is {1, 5., 10.}). Declines a
 * non-finite result, so an overflow reaches the expression path, which carries
 * it into an extended-exponent MPFR real. */
static Expr* ds_real_k(double k, const Expr* n) {
    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
    mpz_abs(v, v);
    mpz_t* primes; unsigned long* exps; size_t np;
    bool fok = df_factor_mpz(v, &primes, &exps, &np);
    mpz_clear(v);
    if (!fok) return NULL;
    Expr* out = NULL;
    if (np == 0) {
        out = expr_new_integer(1);
    } else {
        double total = 1.0;
        for (size_t i = 0; i < np; i++) {
            double q = pow(numeric_mpz_to_double(primes[i]), k), term = 1.0;
            for (unsigned long j = 0; j < exps[i]; j++) term = term * q + 1.0;
            total *= term;
        }
        if (isfinite(total)) out = expr_new_real(total);
    }
    for (size_t i = 0; i < np; i++) mpz_clear(primes[i]);
    free(primes); free(exps);
    return out;
}

/* DivisorSigma[k, n] gives the divisor function sigma_k(n), the sum of the k-th
 * powers of the divisors of n.
 *
 * Three paths, chosen by k (n an Integer, not Gaussian):
 *   machine Integer k   exact, directly from the factorisation (ds_integer_k):
 *                       an Integer for k >= 0, a Rational for k < 0.
 *   machine Real k      in double by Horner (ds_real_k).
 *   anything else       the multiplicative formula
 *                         sigma_k(n) = Prod_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)
 *                       built as an expression and evaluated (ds_assemble), so
 *                       rational / symbolic / complex / arbitrary-precision k
 *                       give radical and symbolic forms.
 * With GaussianIntegers -> True (or a non-real Gaussian-integer n) the product
 * runs over the first-quadrant associates of the Gaussian prime factors, on the
 * expression path for every k. DivisorSigma is Listable, so list arguments are
 * threaded by the evaluator first (an int64 buffer never reaches here: its
 * NDArray kernel in ndinteger.c answers, or declines back to the List).
 * Non-integer or zero n is left unevaluated; argument-count and option errors
 * draw argr / argrx / nonopt / optx / opttf. */
Expr* builtin_divisorsigma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Positional arguments are the leading non-rules; anything after the
     * second is an option (or a nonopt error). */
    size_t npos = 0;
    while (npos < argc && npos < 2 && !ds_is_rule(args[npos])) npos++;
    if (npos != 2) return divisorsigma_emit_argcount(npos);

    bool gaussian = false, gaussian_set = false;
    if (!nt_parse_gaussian_options(res, 2, "DivisorSigma", &gaussian, &gaussian_set))
        return NULL;

    Expr* k = args[0];
    Expr* n = args[1];
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
    if (n->type == EXPR_INTEGER ? n->data.integer == 0
                                : mpz_sgn(n->data.bigint) == 0)
        return NULL;                                      /* DivisorSigma[k, 0] */

    if (k->type == EXPR_INTEGER) {
        Expr* fast = ds_integer_k(k->data.integer, n);
        if (fast) return fast;
    } else if (k->type == EXPR_REAL) {
        Expr* fast = ds_real_k(k->data.real, n);
        if (fast) return fast;
    }

    mpz_t v;
    expr_to_mpz(n, v);                                    /* inits v */
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
