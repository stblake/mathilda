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
#include "ndarray.h"      /* is_ndarray, ndarray_delist_and_reeval */
#include "ndinteger.h"   /* NDArray buffer fast path */
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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "numeric.h"

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
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic) {
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

/* PrimePi[x]: the number of primes <= x, for any real x (a non-integer is
 * floored; an exact numeric such as Pi or E^10 is floored exactly; Infinity
 * is Infinity and -Infinity, like any x < 2, gives 0).
 * PrimePi[] draws PrimePi::argx, a non-rule extra argument ::nonopt and an
 * unknown option ::optx, all leaving the call unevaluated; an unrecognised
 * Method value leaves it unevaluated silently; ProgressReporting is accepted
 * and ignored. Symbolic and complex x stay unevaluated. */
Expr* builtin_primepi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) {
        mth_message_gated(g_arith_warnings_muted, "PrimePi", "argx",
                          "PrimePi called with 0 arguments; 1 argument is expected.");
        return NULL;
    }
    Expr* x_expr = res->data.function.args[0];

    PrimeCountMethod method = PC_AUTOMATIC;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (head_is(opt, SYM_Rule) || head_is(opt, SYM_RuleDelayed))
                       && opt->data.function.arg_count == 2;
        if (!is_rule) {
            char* bad = expr_to_string(opt);
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "PrimePi", "nonopt",
                "Options expected (instead of %s) beyond position 1 in %s. "
                "An option must be a rule or a list of rules.", bad ? bad : "?", call ? call : "?");
            free(bad); free(call);
            return NULL;
        }
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        const char* name = lhs->type == EXPR_SYMBOL ? lhs->data.symbol.name : NULL;
        if (name == SYM_Method) {
            if (!prime_parse_method(rhs, &method)) return NULL;
        } else if (!(name && name == SYM_ProgressReporting)) {
            char* kname = expr_to_string(lhs);
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "PrimePi", "optx",
                "Unknown option %s in %s.", kname ? kname : "?", call ? call : "?");
            free(kname); free(call);
            return NULL;
        }
    }

    /* Infinity and -Infinity. */
    if (x_expr->type == EXPR_SYMBOL && x_expr->data.symbol.name == SYM_Infinity)
        return expr_new_symbol(SYM_Infinity);
    if (head_is(x_expr, SYM_DirectedInfinity) && x_expr->data.function.arg_count == 1) {
        int64_t dir = 0;
        Expr* d = x_expr->data.function.args[0];
        if (d->type == EXPR_INTEGER) dir = d->data.integer;
        if (dir == 1) return expr_new_symbol(SYM_Infinity);
        if (dir == -1) return expr_new_integer(0);
        return NULL;
    }
    if (head_is(x_expr, SYM_Times) && x_expr->data.function.arg_count == 2 &&
        x_expr->data.function.args[0]->type == EXPR_INTEGER &&
        x_expr->data.function.args[0]->data.integer == -1 &&
        x_expr->data.function.args[1]->type == EXPR_SYMBOL &&
        x_expr->data.function.args[1]->data.symbol.name == SYM_Infinity)
        return expr_new_integer(0);

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
        } else if (expr_is_numeric_quantity(x_expr)) {
            /* An exact real numeric (Pi, E^10): its exact Floor. A complex
             * value leaves Floor unevaluated or non-integer, so the call
             * stays as it is. */
            Expr* fl_args[1] = { expr_copy(x_expr) };
            Expr* fl = eval_and_free(expr_new_function(expr_new_symbol(SYM_Floor), fl_args, 1));
            bool ok = fl->type == EXPR_INTEGER;
            x_val = ok ? fl->data.integer : 0;
            expr_free(fl);
            if (!ok) return NULL;
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

/* Cipolla's asymptotic approximation to p_n (accurate for large n). */
/* li(x) by Ramanujan's series, fast-converging for x > 1. */
static double prime_li(double x) {
    const double gamma = 0.57721566490153286061;
    double L = log(x), sum = 0.0, term = 1.0, inner = 0.0;
    for (int k = 1; k < 200; k++) {
        term *= L / (double)k;                              /* L^k / k! */
        if ((k - 1) % 2 == 0) inner += 1.0 / (double)k;     /* sum_{j <= (k-1)/2} 1/(2j+1) */
        double t = ((k % 2) ? 1.0 : -1.0) * term / ldexp(1.0, k - 1) * inner;
        sum += t;
        if (fabs(t) < 1e-17 * fabs(sum)) break;
    }
    return gamma + log(L) + sqrt(x) * sum;
}

/* Riemann's R(x) to three terms: li(x) - li(x^1/2)/2 - li(x^1/3)/3. Its error
 * against pi(x) is O(sqrt(x)/ln x), far below li's own bias, so the n-th prime
 * it predicts is a few sqrt(p_n) numbers from the true one. */
static double prime_R(double x) {
    double r = prime_li(x);
    if (x > 4.0) r -= prime_li(sqrt(x)) / 2.0;
    if (x > 8.0) r -= prime_li(cbrt(x)) / 3.0;
    return r;
}

/* p_n estimated by Newton on R(x) = n, from the classical asymptotic start. */
static double prime_estimate(double n) {
    double L  = log(n);
    double LL = log(L);
    double x = n * (L + LL - 1.0 + (LL - 2.0) / L);
    for (int i = 0; i < 8; i++) {
        double nx = x - (prime_R(x) - n) * log(x);
        if (fabs(nx - x) < 0.5) { x = nx; break; }
        x = nx;
    }
    return x;
}

/* Emit `Prime::intpp: Positive integer argument expected in <call>.` */
static Expr* prime_emit_intpp(Expr* res) {
    if (!g_arith_warnings_muted) {
        char* call = expr_to_string(res);
        mth_message("Prime", "intpp", "Positive integer argument expected in %s.", call ? call : "?");
        free(call);
    }
    return NULL;
}

Expr* builtin_prime(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* An array of indices is one SIEVE, not one Meissel count per element.
     *
     * The degrade is not optional. Marking Prime packed-aware stops the gate
     * materialising its argument, so the buffer now reaches here on BOTH
     * surfaces -- and every branch below tests for EXPR_INTEGER / EXPR_BIGINT
     * and would answer an NDArray with the unevaluated call. Whatever the sieve
     * declines has to go back through the List path explicitly. */
    if (res->data.function.arg_count == 1 &&
        is_ndarray(res->data.function.args[0])) {
        Expr* nd = ndint_prime(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        if (!g_arith_warnings_muted) {
            mth_message("Prime", "argx", "Prime called with %zu argument%s; 1 argument is expected.", argc, argc == 1 ? "" : "s");
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

    /* Large n: estimate p_n from Riemann's R, count exactly there (a second
     * count only if the estimate was unexpectedly far off), then sieve. */
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

    /* Sieve the gap (O(sqrt x) numbers from the R estimate) to the n-th prime. */
    int64_t p = primecount_nth_from(x, c, n);
    return p > 0 ? expr_new_integer(p) : NULL;
}
