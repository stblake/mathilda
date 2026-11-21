/*
 * facint.c — PrimeQ, PrimePowerQ, NextPrime, EulerPhi, and the
 * facint_factor_complete API.  Integer factorisation itself (FactorInteger and every
 * FactorInteger`Method) lives in src/numbertheory/factorinteger/.
 */
#include "linkedlist.h"
#include "facint.h"
#include "numbertheory.h"   /* nt_parse_gaussian_options */
#include "fi_internal.h"
#include "arithmetic.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "gmp_compat.h"   /* mpz_prevprime fallback for GMP < 6.3.0 */
#include "common.h"       /* builtin_arg_error */
#include "message.h"
#include "print.h"
#include "ndarray.h"
#include "ndinteger.h"
#include "numeric.h"     /* expr_is_numeric_quantity */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpz.h>
#endif
/* See facint.h. The 64-bit arm skips GMP entirely: no limb allocation, one
 * Montgomery setup and at most seven witnesses, versus mpz_probab_prime_p's
 * trial division + BPSW on an mpz. Both are exact below 2^64 (BPSW has no
 * counterexample there; the witness sets are proven), so PrimeQ's answers are
 * unchanged -- only its cost on machine integers drops. */
int mth_is_prime_u64(uint64_t n) {
    return fi_u64_is_prime(n);
}

int mth_is_prime_mpz(const mpz_t n) {
    if (mpz_sgn(n) <= 0) return 0;
    if (mpz_sizeinbase(n, 2) <= 64) {
        /* mpz_export, not mpz_get_ui: unsigned long is 32 bits on LLP64. */
        uint64_t v = 0;
        mpz_export(&v, NULL, -1, sizeof v, 0, 0, n);
        return fi_u64_is_prime(v);
    }
    return mpz_probab_prime_p(n, 25) != 0;
}

/* Gaussian-integer primality test for a complex number a + b*I with
 * integer parts. Sets *out_is_gprime to 1/0. Inputs `a` and `b` may carry
 * any sign; they are taken in absolute value internally. */
static void gaussian_prime_test(const mpz_t a_in, const mpz_t b_in, int* out_is_gprime) {
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);
    mpz_abs(a, a_in);
    mpz_abs(b, b_in);
    int is_gprime = 0;
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
        is_gprime = 0;  /* 0 + 0i is not prime */
    } else if (mpz_sgn(b) == 0) {
        /* Pure real: Gaussian prime iff |a| is prime and |a| ≡ 3 (mod 4). */
        if (mpz_probab_prime_p(a, 25) && mpz_fdiv_ui(a, 4) == 3) {
            is_gprime = 1;
        }
    } else if (mpz_sgn(a) == 0) {
        /* Pure imaginary: Gaussian prime iff |b| is prime and |b| ≡ 3 (mod 4). */
        if (mpz_probab_prime_p(b, 25) && mpz_fdiv_ui(b, 4) == 3) {
            is_gprime = 1;
        }
    } else {
        /* Both nonzero: Gaussian prime iff a^2 + b^2 is prime. */
        mpz_t norm;
        mpz_init(norm);
        mpz_mul(norm, a, a);
        mpz_addmul(norm, b, b);
        is_gprime = mpz_probab_prime_p(norm, 25);
        mpz_clear(norm);
    }
    mpz_clear(a);
    mpz_clear(b);
    *out_is_gprime = is_gprime;
}

/* Parse PrimeQ option arguments after the first. Returns true on success,
 * filling *out_gaussian. Returns false if any option is unrecognized or
 * malformed (caller should return False from the predicate in that case
 * to honour the *Q always-boolean contract). */
/* PrimeQ[n] / PrimeQ[n, GaussianIntegers -> True]. True or False for every
 * well-formed call -- anything that is not a prime integer (or, in Z[i], a
 * Gaussian prime) is False. A MALFORMED call is not a question about n:
 * PrimeQ[] draws PrimeQ::argx and a bad option ::nonopt / ::optx
 * / ::opttf, and those calls stay unevaluated. */
Expr* builtin_primeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return expr_new_symbol(SYM_False);
    if (res->data.function.arg_count < 1) {
        mth_message_gated(g_arith_warnings_muted, "PrimeQ", "argx",
                          "PrimeQ called with 0 arguments; 1 argument is expected.");
        return NULL;
    }
    Expr* arg = res->data.function.args[0];

    bool gaussian = false, gaussian_set = false;
    if (!nt_parse_gaussian_options(res, 1, "PrimeQ", &gaussian, &gaussian_set))
        return NULL;

    /* A packed / visible int64 buffer: one pass, a packed bool List. The
     * Gaussian question, and every shape the kernel declines, degrades. */
    if (is_ndarray(arg)) {
        Expr* nd = (!gaussian && res->data.function.arg_count == 1) ? ndint_primeq(res) : NULL;
        return nd ? nd : ndarray_delist_and_reeval(res);
    }

    if (arg->type == EXPR_INTEGER || arg->type == EXPR_BIGINT) {
        mpz_t n;
        expr_to_mpz(arg, n);
        int is_prime = 0;
        if (gaussian) {
            /* A rational integer n is a Gaussian prime iff |n| is a prime
             * in Z AND |n| ≡ 3 (mod 4). (n=2 = -i(1+i)^2, and primes ≡ 1
             * mod 4 split as (a+bi)(a-bi).) */
            mpz_t zero;
            mpz_init_set_ui(zero, 0);
            gaussian_prime_test(n, zero, &is_prime);
            mpz_clear(zero);
        } else {
            mpz_abs(n, n);
            if (mpz_sgn(n) > 0) {
                is_prime = mth_is_prime_mpz(n);
            }
        }
        mpz_clear(n);
        return expr_new_symbol(is_prime ? "True" : "False");
    }

    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
        const char* head = arg->data.function.head->data.symbol.name;
        /* Gaussian prime test for Complex[a, b] with integer parts. Works
         * in both option-on and option-off modes since membership in Z[i]
         * is the relevant universe in either case. */
        if (head == SYM_Complex && arg->data.function.arg_count == 2) {
            Expr* re_expr = arg->data.function.args[0];
            Expr* im_expr = arg->data.function.args[1];
            bool re_int = (re_expr->type == EXPR_INTEGER || re_expr->type == EXPR_BIGINT);
            bool im_int = (im_expr->type == EXPR_INTEGER || im_expr->type == EXPR_BIGINT);
            if (!re_int || !im_int) {
                return expr_new_symbol(SYM_False);
            }
            mpz_t a, b;
            expr_to_mpz(re_expr, a);
            expr_to_mpz(im_expr, b);
            int is_gprime = 0;
            gaussian_prime_test(a, b, &is_gprime);
            mpz_clear(a);
            mpz_clear(b);
            return expr_new_symbol(is_gprime ? "True" : "False");
        }
    }

    /* Everything else (symbols, reals, rationals, strings, symbolic
     * functions, etc.) is not a (Gaussian) prime. *Q predicates must
     * always return True or False — never remain unevaluated. */
    return expr_new_symbol(SYM_False);
}

/* --------------------------------------------------------------------- *
 *  PrimePowerQ                                                          *
 *                                                                       *
 *  n is a prime power iff n = p^k for a prime p and k >= 1.  Every test *
 *  below reduces n to its PRIMITIVE ROOT (the b with n = b^e, b not a   *
 *  perfect power) and then makes exactly one primality test, on b --    *
 *  the smallest candidate there is.  Before that, trial division by the *
 *  small primes settles almost every input outright: a small prime p | n *
 *  answers "is n a power of p?" with one mpz_remove, and once no prime  *
 *  below B divides n the base is > B, which bounds the exponent search  *
 *  by log_B(n).                                                         *
 * --------------------------------------------------------------------- */

/* The odd primes below 64: trial divisors for the machine path.  Once none
 * divides n, any prime base is >= 67, so n < 2^64 forces k <= 10 and the
 * only prime exponents left to try are 2, 3, 5, 7. */
static const uint8_t ppq_odd_primes_u64[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61
};

/* True iff r^k > n, without overflow. */
static bool u64_pow_exceeds(uint64_t r, unsigned k, uint64_t n) {
    if (r <= 1) return r > n;
    uint64_t acc = 1;
    for (unsigned i = 0; i < k; i++) {
        if (acc > n / r) return true;
        acc *= r;
    }
    return acc > n;
}

/* floor(n^(1/k)) for k >= 2: a double estimate, then exact correction. */
static uint64_t u64_iroot(uint64_t n, unsigned k) {
    uint64_t r = (uint64_t)pow((double)n, 1.0 / (double)k);
    while (r > 0 && u64_pow_exceeds(r, k, n)) r--;
    while (!u64_pow_exceeds(r + 1, k, n)) r++;
    return r;
}

/* See facint.h. */
bool mth_prime_power_base_u64(uint64_t n, uint64_t* p_out) {
    for (;;) {
        if (n < 2) return false;
        if ((n & 1) == 0) {
            /* Power of two iff a single bit is set. */
            if (n & (n - 1)) return false;
            if (p_out) *p_out = 2;
            return true;
        }
        for (size_t i = 0; i < sizeof ppq_odd_primes_u64; i++) {
            uint64_t p = ppq_odd_primes_u64[i];
            if (n % p == 0) {
                do n /= p; while (n % p == 0);
                if (n != 1) return false;
                if (p_out) *p_out = p;
                return true;
            }
        }
        if (fi_u64_is_prime(n)) {
            if (p_out) *p_out = n;
            return true;
        }
        /* Composite with every prime factor >= 67: a prime power only as a
         * perfect power with exponent <= 10.  Strip one prime exponent and
         * go round again (the root may itself be a power, 3^(2*5) say). */
        static const unsigned ks[] = { 2, 3, 5, 7 };
        bool reduced = false;
        for (size_t i = 0; i < sizeof ks / sizeof ks[0]; i++) {
            uint64_t r = u64_iroot(n, ks[i]);
            if (!u64_pow_exceeds(r, ks[i], n) && u64_pow_exceeds(r, ks[i], n - 1)) {
                /* r^k <= n and r^k > n-1, i.e. r^k == n. */
                n = r;
                reduced = true;
                break;
            }
        }
        if (!reduced) return false;
    }
}

/* Primes below 2^10, the bignum path's trial divisors (filled by
 * facint_init).  Past them every prime base is >= 1031 > 2^10, so a
 * bits-bit n = b^k has k <= bits/10. */
#define PPQ_TRIAL_BOUND 1024
static uint16_t ppq_trial_primes[PPQ_TRIAL_BOUND / 4];
static size_t   ppq_n_trial_primes = 0;

static void ppq_trial_primes_init(void) {
    if (ppq_n_trial_primes) return;
    unsigned char composite[PPQ_TRIAL_BOUND];
    memset(composite, 0, sizeof composite);
    for (unsigned i = 2; i < PPQ_TRIAL_BOUND; i++) {
        if (composite[i]) continue;
        ppq_trial_primes[ppq_n_trial_primes++] = (uint16_t)i;
        for (unsigned j = i * i; j < PPQ_TRIAL_BOUND; j += i) composite[j] = 1;
    }
}

/* Replace b (>= 2, no prime factor below PPQ_TRIAL_BOUND) by its primitive
 * root.  Neither FLINT's fmpz_is_perfect_power nor a single mpz_root
 * promises the LARGEST exponent, so both arms iterate until b is no
 * longer a perfect power. */
static void ppq_primitive_root(mpz_t b) {
#ifdef USE_FLINT
    fmpz_t f, r;
    fmpz_init(f);
    fmpz_init(r);
    fmpz_set_mpz(f, b);
    while (fmpz_is_perfect_power(r, f) && fmpz_cmp_ui(r, 1) > 0
           && fmpz_cmp(r, f) < 0) {
        fmpz_swap(f, r);
    }
    fmpz_get_mpz(b, f);
    fmpz_clear(f);
    fmpz_clear(r);
#else
    mpz_t r;
    mpz_init(r);
    for (;;) {
        size_t bits = mpz_sizeinbase(b, 2);
        /* mpz_perfect_power_p's scratch is partly alloca'd, proportional to
         * the operand (see find_min_perfect_base in power.c): use the fast
         * gate only below a size the C stack certainly survives. */
        if (bits <= 200000 && !mpz_perfect_power_p(b)) break;
        bool reduced = false;
        for (unsigned long k = 2; k <= bits / 10; k++) {
            if (!fi_u64_is_prime(k)) continue;
            if (mpz_root(r, b, k)) {
                mpz_swap(b, r);
                reduced = true;
                break;
            }
        }
        if (!reduced) break;
    }
    mpz_clear(r);
#endif
}

/* See facint.h. */
bool mth_prime_power_base_mpz(const mpz_t n, mpz_t p_out) {
    if (mpz_sgn(n) <= 0) return false;
    if (mpz_sizeinbase(n, 2) <= 64) {
        /* mpz_export/import, not mpz_get_ui/set_ui: unsigned long is 32 bits
         * on LLP64. */
        uint64_t v = 0, p = 0;
        mpz_export(&v, NULL, -1, sizeof v, 0, 0, n);
        if (!mth_prime_power_base_u64(v, &p)) return false;
        if (p_out) mpz_import(p_out, 1, -1, sizeof p, 0, 0, &p);
        return true;
    }
    /* Even: a power of two iff exactly one bit is set. */
    if (mpz_even_p(n)) {
        if (mpz_popcount(n) != 1) return false;
        if (p_out) mpz_set_ui(p_out, 2);
        return true;
    }
    ppq_trial_primes_init();
    mpz_t b;
    mpz_init(b);
    for (size_t i = 1; i < ppq_n_trial_primes; i++) {      /* [0] is 2 */
        unsigned long p = ppq_trial_primes[i];
        if (mpz_divisible_ui_p(n, p)) {
            /* mpz_remove strips p^k by repeated squaring, not k divisions. */
            mpz_t pz;
            mpz_init_set_ui(pz, p);
            mpz_remove(b, n, pz);
            bool ok = mpz_cmp_ui(b, 1) == 0;
            if (ok && p_out) mpz_set(p_out, pz);
            mpz_clear(pz);
            mpz_clear(b);
            return ok;
        }
    }
    mpz_set(b, n);
    ppq_primitive_root(b);
    bool ok = mth_is_prime_mpz(b) != 0;
    if (ok && p_out) mpz_set(p_out, b);
    mpz_clear(b);
    return ok;
}

/* Is a + b I (exact integer parts, any signs) a Gaussian prime power, i.e.
 * a unit times pi^k for one Gaussian prime pi and k >= 1?  The primes of
 * Z[i] are 1+I (over 2), each rational prime q = 3 (mod 4) (inert, norm
 * q^2), and the two conjugate factors pi, conj(pi) of each p = 1 (mod 4).
 *
 *   on an axis (a b = 0): z = unit * m, m = |a| + |b|.  m is a Gaussian
 *     prime power iff m = 2^k (= unit * (1+I)^(2k)) or m = q^j with
 *     q = 3 (mod 4); m = p^j with p = 1 (mod 4) is pi^j conj(pi)^j.
 *   off the axes: N = a^2 + b^2 must be p^k.  p = 2 always works (1+I is
 *     the only prime over 2); p = 3 (mod 4) cannot occur off the axes;
 *     p = 1 (mod 4) works iff z is not divisible by p = pi conj(pi),
 *     i.e. iff p does not divide both parts. */
static bool gaussian_prime_power(const mpz_t a, const mpz_t b) {
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) return false;
    mpz_t m, p;
    mpz_init(m);
    mpz_init(p);
    bool ok = false;
    if (mpz_sgn(a) == 0 || mpz_sgn(b) == 0) {
        mpz_abs(m, mpz_sgn(a) ? a : b);
        if (mth_prime_power_base_mpz(m, p))
            ok = mpz_cmp_ui(p, 2) == 0 || mpz_fdiv_ui(p, 4) == 3;
    } else {
        mpz_mul(m, a, a);
        mpz_addmul(m, b, b);
        if (mth_prime_power_base_mpz(m, p)) {
            if (mpz_cmp_ui(p, 2) == 0) ok = true;
            else if (mpz_fdiv_ui(p, 4) == 1)
                ok = !(mpz_divisible_p(a, p) && mpz_divisible_p(b, p));
        }
    }
    mpz_clear(m);
    mpz_clear(p);
    return ok;
}

/* An inexact number: machine real, MPFR real, or Complex with an inexact
 * part -- the arguments that draw PrimePowerQ::exact. */
static bool ppq_is_inexact_number(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2)
        return ppq_is_inexact_number(e->data.function.args[0])
            || ppq_is_inexact_number(e->data.function.args[1]);
    return false;
}

/* Read an exact Integer or Rational[p, q] into an initialised mpq. */
static bool ppq_to_mpq(const Expr* e, mpq_t out) {
    if (expr_is_integer_like(e)) {
        mpz_t z;
        expr_to_mpz(e, z);
        mpq_set_z(out, z);
        mpz_clear(z);
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpz_t num, den;
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        if (mpz_sgn(den) == 0) { mpz_clear(num); mpz_clear(den); return false; }
        mpq_set_num(out, num);
        mpq_set_den(out, den);
        mpq_canonicalize(out);
        mpz_clear(num);
        mpz_clear(den);
        return true;
    }
    return false;
}

/* PPQ_GI_DELAYED: a well-formed GaussianIntegers :> v.  Mathematica 14
 * answers False for every argument then (PrimePowerQ[4, GaussianIntegers :>
 * True] is False, though PrimeQ honours :>); matched deliberately, since the
 * reference behaviour is the specification. */
typedef enum { PPQ_GI_AUTOMATIC, PPQ_GI_TRUE, PPQ_GI_FALSE, PPQ_GI_DELAYED } PpqGaussian;

/* Is the rational number x (+ y I when y != 0, or when `gaussian`) a prime
 * power?  A rational or Gaussian rational z is one iff z or 1/z is an
 * (ordinary or Gaussian) integer prime power -- 1/25 and I/2 = 1/(-2 I)
 * are, 2/3 and 3/(2+I) are not. */
static bool ppq_exact(mpq_t x, mpq_t y, bool gaussian) {
    if (mpz_cmp_ui(mpq_denref(x), 1) != 0 || mpz_cmp_ui(mpq_denref(y), 1) != 0) {
        /* Not integral: replace z by 1/z = (x - y I) / (x^2 + y^2). */
        mpq_t nrm, t;
        mpq_init(nrm);
        mpq_init(t);
        mpq_mul(nrm, x, x);
        mpq_mul(t, y, y);
        mpq_add(nrm, nrm, t);
        mpq_div(x, x, nrm);
        mpq_div(y, y, nrm);
        mpq_neg(y, y);
        mpq_clear(nrm);
        mpq_clear(t);
        if (mpz_cmp_ui(mpq_denref(x), 1) != 0 || mpz_cmp_ui(mpq_denref(y), 1) != 0)
            return false;
    }
    if (gaussian || mpq_sgn(y) != 0)
        return gaussian_prime_power(mpq_numref(x), mpq_numref(y));
    mpz_t n;
    mpz_init(n);
    mpz_abs(n, mpq_numref(x));
    bool ok = mth_prime_power_base_mpz(n, NULL);
    mpz_clear(n);
    return ok;
}

/* Parse PrimePowerQ's optional second argument, GaussianIntegers -> v.
 * Emits the Wolfram message and returns false (leave the call unevaluated)
 * for a non-rule (nonopt), an unknown option name (optx) or a value other
 * than True / False / Automatic (opttfa). */
static bool ppq_parse_option(Expr* res, PpqGaussian* out) {
    Expr* opt = res->data.function.args[1];
    bool is_rule = opt->type == EXPR_FUNCTION
        && opt->data.function.head->type == EXPR_SYMBOL
        && (opt->data.function.head->data.symbol.name == SYM_Rule
            || opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
        && opt->data.function.arg_count == 2;
    if (!is_rule) {
        char* bad = expr_to_string(opt);
        char* call = expr_to_string(res);
        mth_message("PrimePowerQ", "nonopt",
                    "Options expected (instead of %s) beyond position 1 in %s. "
                    "An option must be a rule or a list of rules.",
                    bad ? bad : "?", call ? call : "?");
        free(bad);
        free(call);
        return false;
    }
    Expr* key = opt->data.function.args[0];
    Expr* val = opt->data.function.args[1];
    if (!(key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_GaussianIntegers)) {
        char* bad = expr_to_string(opt);
        mth_message("PrimePowerQ", "optx", "Unknown option %s in PrimePowerQ.",
                    bad ? bad : "?");
        free(bad);
        return false;
    }
    const char* v = val->type == EXPR_SYMBOL ? val->data.symbol.name : NULL;
    bool valid = v == SYM_True || v == SYM_False || v == SYM_Automatic;
    if (valid && opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
        *out = PPQ_GI_DELAYED;
    else if (v == SYM_True)      *out = PPQ_GI_TRUE;
    else if (v == SYM_False)     *out = PPQ_GI_FALSE;
    else if (v == SYM_Automatic) *out = PPQ_GI_AUTOMATIC;
    else {
        char* bad = expr_to_string(opt);
        mth_message("PrimePowerQ", "opttfa",
                    "Value of option %s should be True, False, or Automatic.",
                    bad ? bad : "?");
        free(bad);
        return false;
    }
    return true;
}

/* PrimePowerQ[expr] / PrimePowerQ[expr, GaussianIntegers -> v].
 *
 * True iff expr is p^k (k >= 1) for a prime p, up to sign: integers over Z
 * (Gaussian integers over Z[i]), and rationals / Gaussian rationals whose
 * reciprocal is one.  Any other exact expression is False; an inexact
 * number draws PrimePowerQ::exact and stays unevaluated, as in Mathematica. */
Expr* builtin_primepowerq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("PrimePowerQ", argc, 1, 2);

    PpqGaussian mode = PPQ_GI_AUTOMATIC;
    if (argc == 2 && !ppq_parse_option(res, &mode)) return NULL;

    Expr* arg = res->data.function.args[0];

    /* A packed / visible int64 buffer answers in one pass straight off the
     * buffer; every shape the kernel declines degrades element by element. */
    if (is_ndarray(arg)) {
        Expr* nd = ndint_prime_power_predicate(res);
        if (nd) return nd;
        return ndarray_delist_and_reeval(res);
    }

    if (ppq_is_inexact_number(arg)) {
        char* s = expr_to_string(arg);
        mth_message("PrimePowerQ", "exact",
                    "Argument %s in PrimePowerQ is not an exact number.", s ? s : "?");
        free(s);
        return NULL;
    }
    if (mode == PPQ_GI_DELAYED) return expr_new_symbol(SYM_False);

    /* Machine integers skip GMP entirely (the common case, e.g. Select over
     * a Range). The magnitude is taken in uint64 so INT64_MIN = -2^63 works. */
    if (arg->type == EXPR_INTEGER && mode != PPQ_GI_TRUE) {
        int64_t v = arg->data.integer;
        uint64_t m = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
        return expr_new_symbol(mth_prime_power_base_u64(m, NULL) ? SYM_True : SYM_False);
    }

    bool is_complex = arg->type == EXPR_FUNCTION
        && arg->data.function.head->type == EXPR_SYMBOL
        && arg->data.function.head->data.symbol.name == SYM_Complex
        && arg->data.function.arg_count == 2;
    if (is_complex && mode == PPQ_GI_FALSE) {
        mth_message("PrimePowerQ", "par",
                    "Inappropriate parameter: GaussianIntegers -> False.");
        return NULL;
    }

    mpq_t x, y;
    mpq_init(x);
    mpq_init(y);
    bool exact = is_complex
        ? (ppq_to_mpq(arg->data.function.args[0], x)
           && ppq_to_mpq(arg->data.function.args[1], y))
        : ppq_to_mpq(arg, x);
    bool ok = exact && ppq_exact(x, y, mode == PPQ_GI_TRUE);
    mpq_clear(x);
    mpq_clear(y);
    return expr_new_symbol(ok ? SYM_True : SYM_False);
}

/* The signed primes, in signed order: ..., -7, -5, -3, -2, 2, 3,
 * 5, 7, ... (|p| prime). np_up sets cur to the least such prime > cur,
 * np_down to the greatest < cur. */
static void np_up(mpz_t cur) {
    if (mpz_cmp_si(cur, 2) >= 0) { mpz_nextprime(cur, cur); return; }
    if (mpz_cmp_si(cur, -2) >= 0) { mpz_set_ui(cur, 2); return; }   /* -2 .. 1 -> 2 */
    mpz_neg(cur, cur);                                             /* cur < -2 */
    mathilda_mpz_prevprime(cur, cur);                              /* largest prime < |cur| */
    mpz_neg(cur, cur);
}

static void np_down(mpz_t cur) {
    if (mpz_cmp_si(cur, 2) > 0) { mathilda_mpz_prevprime(cur, cur); return; }
    if (mpz_cmp_si(cur, -2) > 0) { mpz_set_si(cur, -2); return; }  /* -1 .. 2 -> -2 */
    mpz_neg(cur, cur);                                             /* cur <= -2 */
    mpz_nextprime(cur, cur);
    mpz_neg(cur, cur);
}

/* NextPrime[x] / NextPrime[x, k]: the k-th prime after x (k > 0) or before it
 * (k < 0), counting the NEGATIVE primes too:
 * NextPrime[-10] is -7, NextPrime[3, -2] is -2. NextPrime[x, 0] is x when x is
 * a (signed) prime and unevaluated otherwise. x may be any real numeric
 * quantity -- an exact rational of any size, a real, Pi, Sqrt[2], ... (a
 * non-integer starts from Floor[x] going up and Ceiling[x] going down); k must
 * be an integer. Not Listable, but a List in either
 * position is threaded explicitly. A wrong argument count draws
 * NextPrime::argt. */
Expr* builtin_nextprime(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) {
        mth_message_gated(g_arith_warnings_muted, "NextPrime", "argt",
            "NextPrime called with %zu arguments; 1 or 2 arguments are expected.", argc);
        return NULL;
    }
    Expr** args = res->data.function.args;

    /* A LinkedList in either position threads as a List does, every List /
     * LinkedList argument in step; the answer is a LinkedList. */
    for (size_t i = 0; i < argc; i++)
        if (args[i]->type == EXPR_LINKEDLIST) {
            Expr* t = ll_thread_listable(res);
            if (t && t->type != EXPR_LINKEDLIST && expr_eq(t, res)) { expr_free(t); return NULL; }
            return t;
        }

    /* Thread over a List in either position. (Not on pack.c's AWARE list, so a
     * packed buffer has already been unpacked by the gate.) */
    for (size_t i = 0; i < argc; i++) {
        if (args[i]->type == EXPR_FUNCTION && args[i]->data.function.head->type == EXPR_SYMBOL &&
            args[i]->data.function.head->data.symbol.name == SYM_List) {
            size_t len = args[i]->data.function.arg_count;
            for (size_t j = 0; j < argc; j++) {
                Expr* o = args[j];
                if (j != i && o->type == EXPR_FUNCTION && o->data.function.head->type == EXPR_SYMBOL &&
                    o->data.function.head->data.symbol.name == SYM_List &&
                    o->data.function.arg_count != len) return NULL;
            }
            Expr** items = malloc((len ? len : 1) * sizeof(Expr*));
            for (size_t e = 0; e < len; e++) {
                Expr* call_args[2];
                for (size_t j = 0; j < argc; j++) {
                    Expr* o = args[j];
                    bool lst = o->type == EXPR_FUNCTION && o->data.function.head->type == EXPR_SYMBOL &&
                               o->data.function.head->data.symbol.name == SYM_List;
                    call_args[j] = expr_copy(lst ? o->data.function.args[e] : o);
                }
                items[e] = expr_new_function(expr_new_symbol(SYM_NextPrime), call_args, argc);
            }
            Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_List), items, len));
            free(items);
            return out;
        }
    }

    Expr* x_expr = args[0];
    int64_t k = 1;
    if (argc == 2) {
        if (args[1]->type != EXPR_INTEGER) return NULL;    /* k must be an integer */
        k = args[1]->data.integer;
    }

    /* x as an integer start: Floor going up, Ceiling going down. */
    mpz_t cur;
    bool x_int = false;
    if (x_expr->type == EXPR_INTEGER || x_expr->type == EXPR_BIGINT) {
        expr_to_mpz(x_expr, cur);
        x_int = true;
    } else if (x_expr->type == EXPR_REAL && isfinite(x_expr->data.real)) {
        mpz_init_set_d(cur, k >= 0 ? floor(x_expr->data.real) : ceil(x_expr->data.real));
    } else {
        mpq_t q;
        mpq_init(q);
        if (ppq_to_mpq(x_expr, q)) {
            /* Exact rational of any size: GMP floor/ceiling division. */
            mpz_init(cur);
            if (k >= 0) mpz_fdiv_q(cur, mpq_numref(q), mpq_denref(q));
            else        mpz_cdiv_q(cur, mpq_numref(q), mpq_denref(q));
            mpq_clear(q);
        } else {
            mpq_clear(q);
            /* Any other real numeric quantity (Pi, Sqrt[2], an MPFR real, ...):
             * let Floor/Ceiling decide it exactly. A symbolic or complex x
             * does not round to an integer, so NextPrime stays unevaluated. */
            if (!expr_is_numeric_quantity(x_expr)) return NULL;
            Expr* fa[1] = { expr_copy(x_expr) };
            Expr* r = eval_and_free(expr_new_function(
                expr_new_symbol(k >= 0 ? SYM_Floor : SYM_Ceiling), fa, 1));
            if (!expr_is_integer_like(r)) { expr_free(r); return NULL; }
            expr_to_mpz(r, cur);
            expr_free(r);
        }
    }

    if (k == 0) {
        /* Only a (signed) prime is its own 0-th next prime. */
        if (!x_int) { mpz_clear(cur); return NULL; }
        mpz_t a;
        mpz_init(a);
        mpz_abs(a, cur);
        bool prime = mpz_cmp_ui(a, 2) >= 0 && mpz_probab_prime_p(a, 25) > 0;
        mpz_clear(a);
        mpz_clear(cur);
        return prime ? expr_copy(x_expr) : NULL;
    }

    if (k > 0) for (int64_t i = 0; i < k; i++) np_up(cur);
    else       for (int64_t i = 0; i > k; i--) np_down(cur);

    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(cur));
    mpz_clear(cur);
    return result;
}

/* --------------------------------------------------------------------- *
 *  Public complete-factorisation helper (number-field maximal-order gate *
 *  src/numbertheory/numberfield.c, Dirichlet characters).                *
 *                                                                        *
 *  Fully factor |n| with the Automatic method.  The distinct primes are  *
 *  written in ascending order to primes[] (each mpz_init'd on success -- *
 *  caller must mpz_clear) and their exponents to exps[], with *count set *
 *  to the number of distinct primes.                                     *
 *                                                                        *
 *  Returns true iff the factorisation is COMPLETE (no composite survived *
 *  the search) and the distinct-prime count fits in `cap`.  Returns      *
 *  false -- writing nothing, *count = 0 -- otherwise; the caller must    *
 *  then DECLINE rather than trust a possibly-composite "prime".          *
 *  n == 0 returns false; |n| == 1 returns true with *count = 0.          */
bool facint_factor_complete(const mpz_t n, mpz_t* primes, int64_t* exps,
                            int cap, int* count) {
    *count = 0;
    if (mpz_sgn(n) == 0) return false;
    if (mpz_cmpabs_ui(n, 1) == 0) return true;
    FiFactors f;
    fi_factors_init(&f);
    bool complete = fi_factor_complete_mpz(n, &f) && f.n <= (size_t)cap;
    if (complete) {
        for (size_t i = 0; i < f.n; i++) {
            mpz_init_set(primes[i], f.v[i].p);
            exps[i] = f.v[i].e;
        }
        *count = (int)f.n;
    }
    fi_factors_clear(&f);
    return complete;
}

/* EulerPhi[n]: the number of 1 <= k <= n coprime to n. EulerPhi[0] is 0 and
 * the sign of n is ignored. A small machine n goes through
 * ndint_euler_phi_i64 -- the helper behind the NDArray kernel and Compile's
 * OP_PHI_I -- so the three agree and no GMP is touched; otherwise FactorInteger
 * and GMP. A non-integer NUMBER (EulerPhi[1/2], EulerPhi[2.], EulerPhi[3 + I])
 * draws EulerPhi::int and a wrong argument count EulerPhi::argx, both leaving
 * the call unevaluated; a symbolic argument stays unevaluated silently. */
#define PHI_MACHINE_TRIAL_MAX ((int64_t)1 << 20)

Expr* builtin_eulerphi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        mth_message_gated(g_arith_warnings_muted, "EulerPhi", "argx",
            "EulerPhi called with %zu arguments; 1 argument is expected.", argc);
        return NULL;
    }
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_INTEGER && arg->data.integer >= -PHI_MACHINE_TRIAL_MAX &&
        arg->data.integer <= PHI_MACHINE_TRIAL_MAX) {
        int64_t r;
        if (ndint_euler_phi_i64(arg->data.integer, &r)) return expr_new_integer(r);
    }
    if (!expr_is_integer_like(arg)) {
        Expr *re, *im;
        bool number = arg->type == EXPR_REAL || is_rational_like(arg) ||
                      is_complex(arg, &re, &im)
#ifdef USE_MPFR
                      || arg->type == EXPR_MPFR
#endif
                      ;
        if (number) {
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "EulerPhi", "int",
                "Integer expected at position 1 in %s.", call ? call : "?");
            free(call);
        }
        return NULL;
    }

    mpz_t n;
    expr_to_mpz(arg, n);
    if (mpz_sgn(n) == 0) { mpz_clear(n); return expr_new_integer(0); }
    mpz_abs(n, n);                                          /* phi(-n) = phi(n) */
    if (mpz_cmp_ui(n, 1) == 0) { mpz_clear(n); return expr_new_integer(1); }

    FiFactors f;
    fi_factors_init(&f);
    if (!fi_factor_complete_mpz(n, &f)) {
        /* A composite survived the search: phi would be wrong, decline. */
        fi_factors_clear(&f);
        mpz_clear(n);
        return NULL;
    }

    /* phi(n) = prod p^(e-1) (p - 1), applied as (phi / p) * (p - 1). */
    mpz_t phi, pm1;
    mpz_init_set(phi, n);
    mpz_init(pm1);
    for (size_t i = 0; i < f.n; i++) {
        mpz_divexact(phi, phi, f.v[i].p);
        mpz_sub_ui(pm1, f.v[i].p, 1);
        mpz_mul(phi, phi, pm1);
    }
    fi_factors_clear(&f);
    mpz_clears(pm1, n, NULL);
    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(phi));
    mpz_clear(phi);
    return out;
}

void fi_init(void);

void facint_init(void) {
    symtab_add_builtin("PrimeQ", builtin_primeq);
    symtab_add_builtin("PrimePowerQ", builtin_primepowerq);
    symtab_add_builtin("NextPrime", builtin_nextprime);
    symtab_add_builtin("EulerPhi", builtin_eulerphi);

    symtab_get_def("PrimeQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("PrimePowerQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    ppq_trial_primes_init();
    /* PrimePi is registered in numbertheory_init() (src/numbertheory/prime.c). */
    symtab_get_def("EulerPhi")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    /* NextPrime is NOT Listable (it threads lists explicitly). */
    symtab_get_def("NextPrime")->attributes |= ATTR_PROTECTED;

    /* FactorInteger and FactorInteger`<Method> (src/numbertheory/factorinteger/). */
    fi_init();
}
