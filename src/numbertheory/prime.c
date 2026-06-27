/* prime.c -- Prime[n] (the nth prime) and PrimePi[x] (the prime-counting
 * function).
 *
 * The prime-counting algorithms live in primecount.c; this file holds the two
 * builtins.  PrimePi[x] accepts a Method option selecting the algorithm:
 *   Automatic (default), "Sieve", "Legendre", "Meissel", "Lehmer", "LMO",
 *   "DelegliseRivat", "LucyHedgehog".
 *
 * Prime[n] is the functional inverse of PrimePi: small n are read straight from
 * the sieve table; large n are found by seeding Cipolla's asymptotic estimate,
 * refining it with a Newton step driven by the exact counter, then walking with
 * GMP's nextprime/prevprime to land exactly on p_n. */

#include "numbertheory.h"
#include "primecount.h"
#include "arithmetic.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "common.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>

/* During the Newton search for p_n, stop iterating once pi(x) is within this
 * many primes of n, then finish with a nextprime/prevprime walk. */
#define PRIME_WALK_LIMIT 1000000LL
#define PRIME_NEWTON_MAX 12

/* ------------------------------------------------------------------ */
/* PrimePi[x, Method -> ...]                                           */
/* ------------------------------------------------------------------ */

/* Map a Method option value (a string, or the symbol Automatic) to a method.
 * Returns 1 on success, 0 if the value is not a recognised method. */
static int prime_parse_method(const Expr* val, PrimeCountMethod *out) {
    if (val->type == EXPR_SYMBOL && val->data.symbol == SYM_Automatic) {
        *out = PC_AUTOMATIC; return 1;
    }
    if (val->type != EXPR_STRING) return 0;
    const char *s = val->data.string;
    if      (!strcmp(s, "Automatic"))      *out = PC_AUTOMATIC;
    else if (!strcmp(s, "Sieve"))          *out = PC_SIEVE;
    else if (!strcmp(s, "Legendre"))       *out = PC_LEGENDRE;
    else if (!strcmp(s, "Meissel"))        *out = PC_MEISSEL;
    else if (!strcmp(s, "Lehmer"))         *out = PC_LEHMER;
    else if (!strcmp(s, "LMO"))            *out = PC_LMO;
    else if (!strcmp(s, "DelegliseRivat")) *out = PC_DR;
    else if (!strcmp(s, "LucyHedgehog"))   *out = PC_LUCY;
    else return 0;
    return 1;
}

Expr* builtin_primepi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    Expr* x_expr = res->data.function.args[0];

    /* Parse trailing Method -> _ options; anything else leaves it unevaluated. */
    PrimeCountMethod method = PC_AUTOMATIC;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if ((head_is(opt, SYM_Rule) || head_is(opt, SYM_RuleDelayed))
            && opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2) {
            Expr* lhs = opt->data.function.args[0];
            Expr* rhs = opt->data.function.args[1];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == SYM_Method) {
                if (!prime_parse_method(rhs, &method)) {
                    if (!arith_warnings_muted()) {
                        char* s = expr_to_string(rhs);
                        fprintf(stderr,
                                "PrimePi::method: %s is not a recognised setting "
                                "for Method.\n", s ? s : "?");
                        free(s);
                    }
                    return NULL;
                }
                continue;
            }
        }
        return NULL;   /* unrecognised extra argument */
    }

    /* Resolve x to an int64 floor. */
    int64_t x_val;
    if (x_expr->type == EXPR_INTEGER) {
        x_val = x_expr->data.integer;
    } else if (x_expr->type == EXPR_BIGINT) {
        mpz_t z;
        mpz_init(z);
        expr_to_mpz(x_expr, z);
        if (mpz_sgn(z) <= 0) { mpz_clear(z); return expr_new_integer(0); }
        if (!mpz_fits_slong_p(z)) { mpz_clear(z); return NULL; }
        x_val = (int64_t)mpz_get_si(z);
        mpz_clear(z);
    } else if (x_expr->type == EXPR_REAL) {
        x_val = (int64_t)floor(x_expr->data.real);
    } else {
        int64_t n, d;
        if (is_rational(x_expr, &n, &d)) {
            x_val = n / d;
            if (n < 0 && n % d != 0) x_val--;
        } else {
            return NULL;
        }
    }

    if (x_val < 2) return expr_new_integer(0);
    int64_t count = prime_count(x_val, method);
    if (count < 0) return NULL;   /* beyond the method's supported range */
    return expr_new_integer(count);
}

/* ------------------------------------------------------------------ */
/* Prime[n]                                                            */
/* ------------------------------------------------------------------ */

/* Portable mpz <- int64 (avoids assuming a 64-bit long). */
static void mpz_set_i64(mpz_t r, int64_t v) {
    int neg = (v < 0);
    uint64_t u = neg ? (uint64_t)(-(v + 1)) + 1u : (uint64_t)v;
    mpz_set_ui(r, (unsigned long)(u >> 32));
    mpz_mul_2exp(r, r, 32);
    mpz_add_ui(r, r, (unsigned long)(u & 0xFFFFFFFFULL));
    if (neg) mpz_neg(r, r);
}

/* Cipolla's asymptotic approximation to p_n (accurate for large n). */
static double prime_estimate(double n) {
    double L  = log(n);
    double LL = log(L);
    return n * (L + LL - 1.0
                + (LL - 2.0) / L
                - (LL * LL - 6.0 * LL + 11.0) / (2.0 * L * L));
}

/* Given an anchor x with pi(x) == c, walk to the exact nth prime via GMP. */
static Expr* prime_walk(int64_t x, int64_t c, int64_t n) {
    mpz_t cur;
    mpz_init(cur);
    mpz_set_i64(cur, x);

    if (n > c) {
        for (int64_t k = n - c; k > 0; k--) mpz_nextprime(cur, cur);
    } else {
        mpz_add_ui(cur, cur, 1);
        mpz_prevprime(cur, cur);
        for (int64_t k = c - n; k > 0; k--) mpz_prevprime(cur, cur);
    }

    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(cur));
    mpz_clear(cur);
    return out;
}

/* Emit `Prime::intpp: Positive integer argument expected in <call>.` */
static Expr* prime_emit_intpp(Expr* res) {
    if (!arith_warnings_muted()) {
        char* call = expr_to_string(res);
        fprintf(stderr,
                "Prime::intpp: Positive integer argument expected in %s.\n",
                call ? call : "?");
        free(call);
    }
    return NULL;
}

Expr* builtin_prime(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        if (!arith_warnings_muted()) {
            fprintf(stderr,
                    "Prime::argx: Prime called with %zu argument%s; "
                    "1 argument is expected.\n",
                    argc, argc == 1 ? "" : "s");
        }
        return NULL;
    }

    Expr* arg = res->data.function.args[0];

    /* Resolve the index n as a positive int64, or diagnose / defer. */
    int64_t n;
    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer <= 0) return prime_emit_intpp(res);
        n = arg->data.integer;
    } else if (arg->type == EXPR_BIGINT) {
        mpz_t z;
        mpz_init(z);
        expr_to_mpz(arg, z);
        if (mpz_sgn(z) <= 0) { mpz_clear(z); return prime_emit_intpp(res); }
        int fits = mpz_fits_slong_p(z);
        n = fits ? (int64_t)mpz_get_si(z) : 0;
        mpz_clear(z);
        if (!fits) return NULL;   /* valid but far beyond computable range */
    } else if (expr_is_numeric_like(arg)) {
        return prime_emit_intpp(res);
    } else {
        return NULL;              /* symbolic: leave unevaluated */
    }

    /* Small n: read the answer straight out of the prime table. */
    int64_t table_size = primecount_small_table_size();
    if (n <= table_size)
        return expr_new_integer((int64_t)primecount_small_prime(n));

    /* Large n: estimate, Newton-refine against the exact counter, then walk. */
    double est = prime_estimate((double)n);
    if (!(est <= (double)PI_COUNT_MAX)) return NULL;   /* out of range / NaN */

    int64_t x = (int64_t)(est + 0.5);
    if (x < 2) x = 2;
    if (x > PI_COUNT_MAX) x = PI_COUNT_MAX;

    int64_t c = prime_count(x, PC_AUTOMATIC);
    if (c < 0) return NULL;

    for (int iter = 0;
         iter < PRIME_NEWTON_MAX && llabs(n - c) > PRIME_WALK_LIMIT;
         iter++) {
        double lnx  = log((double)x);
        double step = (double)(n - c) * lnx;
        int64_t nx  = x + (int64_t)(step >= 0 ? step + 0.5 : step - 0.5);
        if (nx < 2) nx = 2;
        if (nx > PI_COUNT_MAX) nx = PI_COUNT_MAX;
        if (nx == x) break;
        x = nx;
        c = prime_count(x, PC_AUTOMATIC);
        if (c < 0) return NULL;
    }

    return prime_walk(x, c, n);
}
