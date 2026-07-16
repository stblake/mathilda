#include "facint.h"
#include "arithmetic.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define METHOD_AUTOMATIC 0
#define METHOD_SQUFOF 10
#define METHOD_RBD 9
#define METHOD_DIXON 8
#define METHOD_CFRAC 5
#define METHOD_POLLARD_P1 6
#define METHOD_WILLIAMS_P1 7
#define METHOD_FERMAT 4
#define METHOD_TRIAL 1
#define METHOD_POLLARD_RHO 2
#define METHOD_ECM 3


#define MAX_PRIME_LIMIT 1000000

/* Small-prime table shared by FactorInteger / trial division below.  The
 * prime-counting function PrimePi (and Prime[n]) now live in
 * src/numbertheory/prime.c with their own self-contained infrastructure. */
static uint32_t *pi_primes = NULL;
static uint32_t total_pi_primes = 0;

static void init_primepi() {
    if (pi_primes) return;

    uint8_t *is_prime = calloc(MAX_PRIME_LIMIT + 1, sizeof(uint8_t));
    for (uint32_t i = 2; i <= MAX_PRIME_LIMIT; i++) is_prime[i] = 1;
    for (uint32_t p = 2; p * p <= MAX_PRIME_LIMIT; p++) {
        if (is_prime[p]) {
            for (uint32_t i = p * p; i <= MAX_PRIME_LIMIT; i += p)
                is_prime[i] = 0;
        }
    }
    for (uint32_t i = 2; i <= MAX_PRIME_LIMIT; i++) {
        if (is_prime[i]) total_pi_primes++;
    }
    pi_primes = malloc(total_pi_primes * sizeof(uint32_t));
    uint32_t idx = 0;
    for (uint32_t i = 2; i <= MAX_PRIME_LIMIT; i++) {
        if (is_prime[i]) pi_primes[idx++] = i;
    }
    free(is_prime);
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
static bool primeq_parse_options(const Expr* res, size_t first_opt,
                                 bool* out_gaussian) {
    *out_gaussian = false;
    size_t argc = res->data.function.arg_count;
    for (size_t i = first_opt; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type != EXPR_FUNCTION ||
            a->data.function.head->type != EXPR_SYMBOL ||
            a->data.function.head->data.symbol.name != SYM_Rule ||
            a->data.function.arg_count != 2) {
            return false;
        }
        Expr* key = a->data.function.args[0];
        Expr* val = a->data.function.args[1];
        if (key->type != EXPR_SYMBOL) return false;
        if (key->data.symbol.name == SYM_GaussianIntegers) {
            if (val->type != EXPR_SYMBOL) return false;
            if (strcmp(val->data.symbol.name, "True") == 0) *out_gaussian = true;
            else if (strcmp(val->data.symbol.name, "False") == 0) *out_gaussian = false;
            else return false;
        } else {
            return false;  /* unknown option */
        }
    }
    return true;
}

Expr* builtin_primeq(Expr* res) {
    /* *Q predicates must always evaluate to True or False — never remain
     * unevaluated. Zero-arg or non-function res (the latter unreachable
     * via the evaluator, kept defensive) is therefore not prime. */
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) {
        return expr_new_symbol(SYM_False);
    }
    Expr* arg = res->data.function.args[0];

    /* Parse optional GaussianIntegers->True/False. A malformed option list
     * makes the call evaluate to False (per *Q always-boolean contract). */
    bool gaussian = false;
    if (res->data.function.arg_count > 1) {
        if (!primeq_parse_options(res, 1, &gaussian)) {
            return expr_new_symbol(SYM_False);
        }
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
                is_prime = mpz_probab_prime_p(n, 25);
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

Expr* builtin_nextprime(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;

    Expr* x_expr = res->data.function.args[0];
    int64_t k = 1;
    if (res->data.function.arg_count == 2) {
        if (res->data.function.args[1]->type != EXPR_INTEGER) return NULL;
        k = res->data.function.args[1]->data.integer;
    }

    if (k == 0) return expr_copy(x_expr);

    mpz_t current;
    if (x_expr->type == EXPR_INTEGER || x_expr->type == EXPR_BIGINT) {
        expr_to_mpz(x_expr, current);
    } else if (x_expr->type == EXPR_REAL) {
        mpz_init_set_d(current, floor(x_expr->data.real));
    } else {
        int64_t n, d;
        if (is_rational(x_expr, &n, &d)) {
            int64_t start_val = n / d;
            if (n < 0 && n % d != 0) start_val--; 
            mpz_init_set_si(current, start_val);
        } else {
            return NULL;
        }
    }

    if (k > 0) {
        for (int64_t i = 0; i < k; i++) {
            mpz_nextprime(current, current);
        }
    } else {
        int64_t count = -k;
        for (int64_t i = 0; i < count; i++) {
            if (mpz_cmp_ui(current, 2) <= 0) {
                mpz_clear(current);
                return NULL;
            }
            int found = mpz_prevprime(current, current);
            if (!found) {
                mpz_clear(current);
                return NULL;
            }
        }
    }

    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(current));
    mpz_clear(current);
    return result;
}

#include <gmp.h>
#ifndef NO_ECM
#include "ecm.h"
#endif

#include <gmp.h>
#ifndef NO_ECM
#include "ecm.h"
#endif

typedef struct {
    mpz_t p;
    int64_t count;
} FactorMpz;

static void add_factor_mpz(FactorMpz* factors, int* num_factors, mpz_t p, int64_t count) {
    for (int i = 0; i < *num_factors; i++) {
        if (mpz_cmp(factors[i].p, p) == 0) {
            factors[i].count += count;
            return;
        }
    }
    mpz_init_set(factors[*num_factors].p, p);
    factors[*num_factors].count = count;
    (*num_factors)++;
}

/* Warn that a known-composite number is being returned unfactored.  The
 * Automatic method reaches this only for hard composites (e.g. a semiprime of
 * two large distinct primes) that survive trial division, Pollard rho and ECM
 * within the search bounds: perfect powers are peeled off earlier, so `n` here
 * is genuinely composite but with no factor found.  It is added with exponent
 * 1, which looks like a prime — the warning makes that non-obvious compromise
 * explicit rather than silent. */
static void factint_warn_incomplete(const mpz_t n) {
    if (arith_warnings_muted()) return;
    char* s = mpz_get_str(NULL, 10, n);
    fprintf(stderr,
            "FactorInteger::nofac: %s is composite but no factor was found "
            "within the search bounds; it is returned unfactored with "
            "exponent 1.\n",
            s ? s : "The argument");
    free(s);
}


static void pollard_rho_brent_mpz(mpz_t factor, mpz_t n, int method) {
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }
    if (mpz_probab_prime_p(n, 25)) {
        mpz_set(factor, n);
        return;
    }

    mpz_t x, y, g, r, q, ys, c, diff;
    mpz_inits(x, y, g, r, q, ys, c, diff, NULL);

    unsigned long y_start = 2;
    unsigned long c_val = 1;
    unsigned long m = 128;

    while (y_start < (unsigned long)((method == METHOD_POLLARD_RHO) ? 10 : 100)) {
        mpz_set_ui(x, y_start);
        mpz_set_ui(y, y_start);
        mpz_set_ui(c, c_val);
        mpz_set_ui(g, 1);
        mpz_set_ui(r, 1);
        mpz_set_ui(q, 1);

        int max_iters = (method == METHOD_POLLARD_RHO) ? 28 : 14; 
        while (mpz_cmp_ui(g, 1) == 0 && max_iters-- > 0) {
            mpz_set(x, y);
            unsigned long r_val = mpz_get_ui(r);
            for (unsigned long i = 0; i < r_val; i++) {
                mpz_mul(y, y, y);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
            }

            unsigned long k = 0;
            while (k < r_val && mpz_cmp_ui(g, 1) == 0) {
                mpz_set(ys, y);
                unsigned long limit = (m < r_val - k) ? m : r_val - k;
                for (unsigned long i = 0; i < limit; i++) {
                    mpz_mul(y, y, y);
                    mpz_add(y, y, c);
                    mpz_mod(y, y, n);

                    mpz_sub(diff, x, y);
                    mpz_abs(diff, diff);
                    mpz_mul(q, q, diff);
                    mpz_mod(q, q, n);
                }
                mpz_gcd(g, q, n);
                k += limit;
            }
            mpz_mul_ui(r, r, 2);
        }

        if (mpz_cmp(g, n) == 0) {
            do {
                mpz_mul(ys, ys, ys);
                mpz_add(ys, ys, c);
                mpz_mod(ys, ys, n);

                mpz_sub(diff, x, ys);
                mpz_abs(diff, diff);
                mpz_gcd(g, diff, n);
            } while (mpz_cmp_ui(g, 1) == 0);
        }

        if (mpz_cmp(g, n) != 0 && mpz_cmp_ui(g, 1) > 0) {
            mpz_set(factor, g);
            mpz_clears(x, y, g, r, q, ys, c, diff, NULL);
            return;
        }

        y_start++;
        c_val++;
    }

    mpz_set_ui(factor, 0);
    mpz_clears(x, y, g, r, q, ys, c, diff, NULL);
}


static void fermat_factor_mpz(mpz_t factor, mpz_t n) {
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }
    
    mpz_t a, b2, tmp;
    mpz_inits(a, b2, tmp, NULL);
    
    mpz_sqrt(a, n);
    mpz_mul(tmp, a, a);
    if (mpz_cmp(tmp, n) < 0) {
        mpz_add_ui(a, a, 1);
    }
    
    while (1) {
        mpz_mul(tmp, a, a);
        mpz_sub(b2, tmp, n);
        if (mpz_perfect_square_p(b2)) {
            mpz_sqrt(tmp, b2);
            mpz_sub(factor, a, tmp);
            if (mpz_cmp_ui(factor, 1) == 0) {
                mpz_add(factor, a, tmp);
            }
            break;
        }
        mpz_add_ui(a, a, 1);
    }
    
    mpz_clears(a, b2, tmp, NULL);
}


static void cfrac_factor_mpz(mpz_t factor, mpz_t n) {
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }
    if (mpz_perfect_square_p(n)) {
        mpz_sqrt(factor, n);
        return;
    }

    mpz_t P, Q, a, a0, A, A_prev, S, g, tmp, next_A;
    mpz_inits(P, Q, a, a0, A, A_prev, S, g, tmp, next_A, NULL);

    mpz_sqrt(a0, n);

    mpz_set_ui(P, 0);
    mpz_set_ui(Q, 1);
    mpz_set(a, a0);

    mpz_set_ui(A_prev, 1);
    mpz_mod(A, a0, n);

    unsigned long k = 1;
    mpz_set_ui(factor, 0);

    while (true) {
        mpz_mul(tmp, a, Q);
        mpz_sub(tmp, tmp, P);
        mpz_set(P, tmp);
        
        mpz_mul(tmp, P, P);
        mpz_sub(tmp, n, tmp);
        mpz_divexact(Q, tmp, Q);
        
        mpz_add(tmp, a0, P);
        mpz_fdiv_q(a, tmp, Q);
        
        mpz_mul(tmp, a, A);
        mpz_add(tmp, tmp, A_prev);
        mpz_mod(next_A, tmp, n);
        
        mpz_set(A_prev, A);
        mpz_set(A, next_A);
        
        if (k % 2 == 0) {
            if (mpz_perfect_square_p(Q)) {
                mpz_sqrt(S, Q);
                
                mpz_sub(tmp, A_prev, S);
                mpz_gcd(g, tmp, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                    mpz_set(factor, g);
                    break;
                }
                
                mpz_add(tmp, A_prev, S);
                mpz_gcd(g, tmp, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                    mpz_set(factor, g);
                    break;
                }
            }
        }
        
        k++;
        
        if (mpz_cmp_ui(Q, 1) == 0 && k > 1) {
            break;
        }
    }

    mpz_clears(P, Q, a, a0, A, A_prev, S, g, tmp, next_A, NULL);
}


#define DIXON_K 400

typedef struct {
    mpz_t a;
    int exp_parity[DIXON_K];
    int exp_full[DIXON_K];
} DixonRelation;

static void dixon_factor_mpz(mpz_t factor, mpz_t n) {
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }
    if (mpz_perfect_square_p(n)) {
        mpz_sqrt(factor, n);
        return;
    }

    init_primepi();
    int K = DIXON_K;
    if ((int)total_pi_primes < K) K = total_pi_primes;

    DixonRelation* rels = calloc((K + 16), sizeof(DixonRelation));
    for (int i = 0; i < K + 16; i++) {
        mpz_init(rels[i].a);
    }

    mpz_t z, a, tmp, sqrt_n;
    mpz_inits(z, a, tmp, sqrt_n, NULL);

    int found_rels = 0;
    mpz_sqrt(sqrt_n, n);
    mpz_add_ui(a, sqrt_n, 1);

    int attempts = 0;
    while (found_rels < K + 8 && attempts < 500000) {
        mpz_mul(z, a, a);
        mpz_mod(z, z, n);
        mpz_set(tmp, z);

        int exp_parity[DIXON_K] = {0};
        int exp_full[DIXON_K] = {0};

        for (int i = 0; i < K; i++) {
            uint32_t p = pi_primes[i];
            while (mpz_divisible_ui_p(tmp, p)) {
                mpz_divexact_ui(tmp, tmp, p);
                exp_parity[i] ^= 1;
                exp_full[i]++;
            }
            if (mpz_cmp_ui(tmp, 1) == 0) break;
        }

        if (mpz_cmp_ui(tmp, 1) == 0) {
            mpz_set(rels[found_rels].a, a);
            memcpy(rels[found_rels].exp_parity, exp_parity, sizeof(int) * K);
            memcpy(rels[found_rels].exp_full, exp_full, sizeof(int) * K);
            found_rels++;
        }

        mpz_add_ui(a, a, 1);
        attempts++;
    }

    mpz_clears(z, a, tmp, sqrt_n, NULL);

    if (found_rels < K + 1) {
        mpz_set_ui(factor, 0);
        for (int i = 0; i < K + 16; i++) mpz_clear(rels[i].a);
        free(rels);
        return;
    }

    uint64_t* mat = calloc(found_rels, sizeof(uint64_t) * ((K + 63) / 64));
    uint64_t* history = calloc(found_rels, sizeof(uint64_t) * ((found_rels + 63) / 64));

    int words_k = (K + 63) / 64;
    int words_r = (found_rels + 63) / 64;

    for (int i = 0; i < found_rels; i++) {
        for (int j = 0; j < K; j++) {
            if (rels[i].exp_parity[j]) {
                mat[i * words_k + (j / 64)] |= (1ULL << (j % 64));
            }
        }
        history[i * words_r + (i / 64)] |= (1ULL << (i % 64));
    }

    int pivot_row = 0;
    for (int c = 0; c < K && pivot_row < found_rels; c++) {
        int r = pivot_row;
        while (r < found_rels && (mat[r * words_k + (c / 64)] & (1ULL << (c % 64))) == 0) {
            r++;
        }
        if (r == found_rels) continue;

        if (r != pivot_row) {
            for (int w = 0; w < words_k; w++) {
                uint64_t t = mat[pivot_row * words_k + w];
                mat[pivot_row * words_k + w] = mat[r * words_k + w];
                mat[r * words_k + w] = t;
            }
            for (int w = 0; w < words_r; w++) {
                uint64_t t = history[pivot_row * words_r + w];
                history[pivot_row * words_r + w] = history[r * words_r + w];
                history[r * words_r + w] = t;
            }
        }

        for (int i = 0; i < found_rels; i++) {
            if (i != pivot_row && (mat[i * words_k + (c / 64)] & (1ULL << (c % 64)))) {
                for (int w = 0; w < words_k; w++) {
                    mat[i * words_k + w] ^= mat[pivot_row * words_k + w];
                }
                for (int w = 0; w < words_r; w++) {
                    history[i * words_r + w] ^= history[pivot_row * words_r + w];
                }
            }
        }
        pivot_row++;
    }

    mpz_t X, Y, gcd_res, p_pow;
    mpz_inits(X, Y, gcd_res, p_pow, NULL);

    mpz_set_ui(factor, 0);

    for (int i = pivot_row; i < found_rels; i++) {
        mpz_set_ui(X, 1);
        int total_exp[DIXON_K] = {0};

        for (int j = 0; j < found_rels; j++) {
            if (history[i * words_r + (j / 64)] & (1ULL << (j % 64))) {
                mpz_mul(X, X, rels[j].a);
                mpz_mod(X, X, n);
                for (int k = 0; k < K; k++) {
                    total_exp[k] += rels[j].exp_full[k];
                }
            }
        }

        mpz_set_ui(Y, 1);
        for (int k = 0; k < K; k++) {
            if (total_exp[k] > 0) {
                mpz_ui_pow_ui(p_pow, pi_primes[k], total_exp[k] / 2);
                mpz_mul(Y, Y, p_pow);
                mpz_mod(Y, Y, n);
            }
        }

        mpz_sub(tmp, X, Y);
        mpz_abs(tmp, tmp);
        mpz_gcd(gcd_res, tmp, n);

        if (mpz_cmp_ui(gcd_res, 1) > 0 && mpz_cmp(gcd_res, n) < 0) {
            mpz_set(factor, gcd_res);
            break;
        }
    }

    mpz_clears(X, Y, gcd_res, p_pow, NULL);
    free(mat);
    free(history);
    for (int i = 0; i < K + 16; i++) mpz_clear(rels[i].a);
    free(rels);
}


/* Try a single (a, b) pair: assume p ~ c*(a/b)^k, descend Q <- floor(b*Q/a)
 * and probe gcd(Q +/- j, n) for small j. Returns 1 if a non-trivial factor
 * is set in `factor`, 0 otherwise. Requires a > b > 0. */
static int rbd_try_pair(mpz_t factor, const mpz_t n,
                        const mpz_t a, const mpz_t b) {
    if (mpz_sgn(b) <= 0 || mpz_cmp(a, b) <= 0) return 0;

    mpz_t diff, search_limit_z, tmp;
    mpz_inits(diff, search_limit_z, tmp, NULL);
    mpz_sub(diff, a, b);

    /* search_limit ~ (2a - b - 1)/(a - b) + 2: small for a/b far from 1,
     * larger when the descent step is gentle. */
    mpz_add(tmp, a, diff);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_fdiv_q(search_limit_z, tmp, diff);
    mpz_add_ui(search_limit_z, search_limit_z, 2);

    if (!mpz_fits_ulong_p(search_limit_z)) {
        mpz_clears(diff, search_limit_z, tmp, NULL);
        return 0;
    }
    unsigned long search_limit = mpz_get_ui(search_limit_z);

    mpz_t Q, g;
    mpz_init_set(Q, n);
    mpz_init(g);

    int found = 0;
    while (1) {
        mpz_mul(tmp, b, Q);
        mpz_fdiv_q(Q, tmp, a);

        if (mpz_cmp_ui(Q, 2) < 0) break;

        for (unsigned long j = 0; j <= search_limit; j++) {
            mpz_sub_ui(tmp, Q, j);
            mpz_gcd(g, tmp, n);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                mpz_set(factor, g);
                found = 1;
                goto done;
            }

            if (j > 0) {
                mpz_add_ui(tmp, Q, j);
                mpz_gcd(g, tmp, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                    goto done;
                }
            }
        }
    }

done:
    mpz_clears(diff, search_limit_z, tmp, Q, g, NULL);
    return found;
}

static void rbd_factor_mpz(mpz_t factor, mpz_t n, Expr* method_opt) {
    mpz_set_ui(factor, 0);
    if (mpz_cmp_ui(n, 4) < 0) return;
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }
    if (mpz_perfect_square_p(n)) {
        mpz_sqrt(factor, n);
        return;
    }

    mpz_t a, b;
    mpz_inits(a, b, NULL);
    mpz_set_ui(b, 1);

    bool have_base = false;
    if (method_opt) {
        if (method_opt->type == EXPR_FUNCTION
            && method_opt->data.function.head->type == EXPR_SYMBOL
            && method_opt->data.function.head->data.symbol.name == SYM_Rational
            && method_opt->data.function.arg_count == 2) {
            expr_to_mpz(method_opt->data.function.args[0], a);
            expr_to_mpz(method_opt->data.function.args[1], b);
            have_base = true;
        } else if (method_opt->type == EXPR_INTEGER
                   || method_opt->type == EXPR_BIGINT) {
            expr_to_mpz(method_opt, a);
            have_base = true;
        }
    }

    if (have_base) {
        rbd_try_pair(factor, n, a, b);
        mpz_clears(a, b, NULL);
        return;
    }

    /* No explicit base: walk j = 3, 4, ... and try every coprime partition
     * a + b = j with a > b >= 1.  Mirrors the reference Python
     * factor_rational_base() search, but bounded by a pair budget so a hard
     * composite cannot hang the REPL (the reference runs until interrupted).
     * On exhaustion `factor` stays 0 and the caller treats n as unfactored. */
    const unsigned long RBD_MAX_PAIRS = 200000UL;
    unsigned long pairs_tried = 0;
    for (unsigned long j = 3; pairs_tried < RBD_MAX_PAIRS; j++) {
        for (unsigned long bb = 1; 2 * bb < j; bb++) {
            unsigned long aa = j - bb;
            /* coprimality */
            unsigned long u = aa, v = bb;
            while (v) { unsigned long t = u % v; u = v; v = t; }
            if (u != 1) continue;

            mpz_set_ui(a, aa);
            mpz_set_ui(b, bb);
            if (rbd_try_pair(factor, n, a, b)) {
                mpz_clears(a, b, NULL);
                return;
            }
            if (++pairs_tried >= RBD_MAX_PAIRS) break;
        }
    }

    mpz_clears(a, b, NULL);
}

static int squfof_single(mpz_t f, const mpz_t n, unsigned long k, unsigned long max_iters) {
    mpz_t N, S, P, Q, Q_next, b, P_next, Q_next_next, r;
    mpz_t b0, P2, Q2, Q2_next, b2, P2_next, Q2_next_next;
    mpz_t temp;
    int success = 0;
    
    mpz_inits(N, S, P, Q, Q_next, b, P_next, Q_next_next, r, NULL);
    mpz_inits(b0, P2, Q2, Q2_next, b2, P2_next, Q2_next_next, temp, NULL);
    
    mpz_mul_ui(N, n, k);
    mpz_sqrt(S, N);
    
    mpz_mul(temp, S, S);
    if (mpz_cmp(temp, N) == 0) {
        mpz_gcd(f, n, S);
        if (mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, n) != 0) {
            success = 1;
        }
        goto cleanup;
    }
    
    mpz_set(P, S);
    mpz_set_ui(Q, 1);
    
    mpz_mul(temp, P, P);
    mpz_sub(Q_next, N, temp); // Q_next = N - P^2
    
    unsigned long i = 0;
    while (i < max_iters) {
        if (mpz_cmp_ui(Q_next, 0) == 0) goto cleanup;
        
        mpz_add(temp, S, P);
        mpz_fdiv_q(b, temp, Q_next); // b = (S + P) / Q_next
        
        mpz_mul(temp, b, Q_next);
        mpz_sub(P_next, temp, P); // P_next = b * Q_next - P
        
        mpz_sub(temp, P, P_next);
        mpz_mul(temp, b, temp);
        mpz_add(Q_next_next, Q, temp); // Q_next_next = Q + b * (P - P_next)
        
        mpz_set(Q, Q_next);
        mpz_set(Q_next, Q_next_next);
        mpz_set(P, P_next);
        i++;
        
        if (i % 2 == 0) {
            if (mpz_perfect_square_p(Q)) {
                mpz_sqrt(r, Q);
                
                // Try Phase 2
                mpz_sub(temp, S, P);
                mpz_fdiv_q(b0, temp, r); // b0 = (S - P) / r
                
                mpz_mul(temp, b0, r);
                mpz_add(P2, temp, P); // P2 = b0 * r + P
                
                mpz_set(Q2, r); // Q2 = r
                
                mpz_mul(temp, P2, P2);
                mpz_sub(temp, N, temp);
                if (mpz_cmp_ui(Q2, 0) == 0) continue;
                mpz_fdiv_q(Q2_next, temp, Q2); // Q2_next = (N - P2^2) / Q2
                
                unsigned long j = 0;
                while (j < max_iters) {
                    if (mpz_cmp_ui(Q2_next, 0) == 0) { break; }
                    
                    mpz_add(temp, S, P2);
                    mpz_fdiv_q(b2, temp, Q2_next);
                    
                    mpz_mul(temp, b2, Q2_next);
                    mpz_sub(P2_next, temp, P2);
                    
                    if (mpz_cmp(P2, P2_next) == 0) {
                        mpz_gcd(temp, n, P2);
                        if (mpz_cmp_ui(temp, 1) != 0 && mpz_cmp(temp, n) != 0) {
                            mpz_set(f, temp);
                            success = 1;
                            goto cleanup;
                        }
                        break;
                    }
                    
                    mpz_sub(temp, P2, P2_next);
                    mpz_mul(temp, b2, temp);
                    mpz_add(Q2_next_next, Q2, temp);
                    
                    mpz_set(Q2, Q2_next);
                    mpz_set(Q2_next, Q2_next_next);
                    mpz_set(P2, P2_next);
                    j++;
                }
            }
        }
    }
    
cleanup:
    mpz_clears(N, S, P, Q, Q_next, b, P_next, Q_next_next, r, NULL);
    mpz_clears(b0, P2, Q2, Q2_next, b2, P2_next, Q2_next_next, temp, NULL);
    return success;
}

void squfof_factor_mpz(mpz_t f, const mpz_t n) {
    mpz_set_ui(f, 0);
    unsigned long multipliers[] = {1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 1155};
    int num_mults = sizeof(multipliers)/sizeof(multipliers[0]);
    // The iterations are set to 2e6 which is roughly what would factor numbers up to 10^24 cleanly.
    for (int i=0; i<num_mults; i++) {
        if (squfof_single(f, n, multipliers[i], 200000)) {
            return;
        }
    }
}

static void factorize_mpz(mpz_t n, FactorMpz* factors, int* num_factors, int* k_limit, int method, Expr* method_opt) {
    if (mpz_cmp_ui(n, 1) <= 0) return;
    if (*k_limit > 0 && *num_factors >= *k_limit) return;

    if (mpz_probab_prime_p(n, 25)) {
        add_factor_mpz(factors, num_factors, n, 1);
        return;
    }

    /* Perfect-power reduction.  A prime power p^k of a large prime defeats the
     * general search methods below: Pollard rho collides on the repeated factor
     * (gcd steps yield only 1 or n), and ECM has no small cofactor to grab and
     * cannot crack the large p directly.  Without this branch the fallbacks
     * below would mislabel a composite such as (2^61-1)^2 as prime^1.  Peel the
     * exponent off explicitly: reduce n = base^mult to a primitive base (not
     * itself a perfect power), factor that once, and scale every exponent by
     * mult. */
    if (mpz_perfect_power_p(n)) {
        int64_t mult = 1;
        mpz_t base, root;
        mpz_init_set(base, n);
        mpz_init(root);
        /* Extract exact e-th roots to a fixpoint.  Restarting the e-scan after
         * every hit accumulates the full exponent (e.g. p^6 -> p via 2 then 3)
         * and leaves base primitive. */
        for (int reduced = 1; reduced; ) {
            reduced = 0;
            unsigned long emax = (unsigned long)mpz_sizeinbase(base, 2);
            for (unsigned long e = 2; e <= emax; e++) {
                if (mpz_root(root, base, e)) {   /* exact -> base is an e-th power */
                    mpz_set(base, root);
                    mult *= (int64_t)e;
                    reduced = 1;
                    break;
                }
            }
        }
        mpz_clear(root);

        /* Factor the primitive base into a private list (full factorization,
         * independent of the caller's k_limit), then merge scaled counts.  Using
         * a private list keeps add_factor_mpz's merge-by-prime correct even when
         * the base shares primes with factors already collected. */
        FactorMpz* bfac = malloc(sizeof(FactorMpz) * 1024);
        int bn = 0;
        int bl = -1;
        if (bfac) {
            factorize_mpz(base, bfac, &bn, &bl, method, method_opt);
            for (int i = 0; i < bn; i++) {
                /* Honour the caller's factor-count cap (FactorInteger[n, k]). */
                if (!(*k_limit > 0 && *num_factors >= *k_limit))
                    add_factor_mpz(factors, num_factors, bfac[i].p, bfac[i].count * mult);
                mpz_clear(bfac[i].p);
            }
            free(bfac);
        } else {
            /* Allocation failure: fall back to the (imperfect) whole-n factor. */
            add_factor_mpz(factors, num_factors, n, 1);
        }
        mpz_clear(base);
        return;
    }


    if (method == METHOD_SQUFOF) {
        mpz_t f;
        mpz_init(f);
        squfof_factor_mpz(f, n);
        
        if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);
            
            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
            mpz_clear(f);
            return;
        }
        
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }

    if (method == METHOD_FERMAT) {
        mpz_t f;
        mpz_init(f);
        fermat_factor_mpz(f, n);
        
        if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);
            
            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
            mpz_clear(f);
            return;
        }
        
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }


    if (method == METHOD_CFRAC) {
        mpz_t f;
        mpz_init(f);
        cfrac_factor_mpz(f, n);
        
        if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);
            
            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
            mpz_clear(f);
            return;
        }
        
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }


    if (method == METHOD_DIXON) {
        mpz_t f;
        mpz_init(f);
        dixon_factor_mpz(f, n);
        
        if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);
            
            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
            mpz_clear(f);
            return;
        }
        
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }


    if (method == METHOD_RBD) {
        mpz_t f;
        mpz_init(f);
        rbd_factor_mpz(f, n, method_opt);
        
        if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);
            
            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
            mpz_clear(f);
            return;
        }
        
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }

    // Trial division for small primes
    if (method == METHOD_AUTOMATIC || method == METHOD_TRIAL) {

    init_primepi();
    uint32_t n_small = (method == METHOD_TRIAL) ? 1000 : 25;
    if (n_small > total_pi_primes) n_small = total_pi_primes;

    for (uint32_t i = 0; i < n_small; i++) {
        uint32_t p = pi_primes[i];
        if (mpz_divisible_ui_p(n, p)) {
            int64_t count = 0;
            while (mpz_divisible_ui_p(n, p)) {
                count++;
                mpz_divexact_ui(n, n, p);
            }
            mpz_t sp;
            mpz_init_set_ui(sp, p);
            add_factor_mpz(factors, num_factors, sp, count);
            mpz_clear(sp);
            
            if (mpz_cmp_ui(n, 1) == 0) return;
            if (mpz_probab_prime_p(n, 25)) {
                add_factor_mpz(factors, num_factors, n, 1);
                return;
            }
            if (*k_limit > 0 && *num_factors >= *k_limit) return;
        }
    }


    }
    if (method == METHOD_TRIAL) {
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        return;
    }
    mpz_t f;
    mpz_init(f);

    if (method == METHOD_AUTOMATIC || method == METHOD_POLLARD_RHO) {
        pollard_rho_brent_mpz(f, n, method);

    if (mpz_cmp_ui(f, 0) > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
        mpz_t n_f;
        mpz_init(n_f);
        mpz_divexact(n_f, n, f);
        
        factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
        factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
        mpz_clear(n_f);
        mpz_clear(f);
        return;
    }

    }
    if (method == METHOD_POLLARD_RHO) {
        if (mpz_cmp_ui(n, 1) > 0) add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }
    
#ifndef NO_ECM
    if (method == METHOD_AUTOMATIC || method == METHOD_ECM || method == METHOD_POLLARD_P1 || method == METHOD_WILLIAMS_P1) {
        ecm_params params;
        ecm_init(params);
        if (method == METHOD_POLLARD_P1) params->method = ECM_PM1;
        else if (method == METHOD_WILLIAMS_P1) params->method = ECM_PP1;
        
        params->B1done = 0.0; /* We can do a quick check with small B1 */

        int found = 0;
        // Try with increasing bounds
        double B1_bounds[] = {2000.0, 11000.0, 50000.0, 250000.0, 1000000.0, 3000000.0, 11000000.0};
        int num_bounds = sizeof(B1_bounds)/sizeof(B1_bounds[0]);
        
        for (int i = 0; i < num_bounds && !found; i++) {
            for (int tries = 0; tries < 10 && !found; tries++) {
                int ret = ecm_factor(f, n, B1_bounds[i], params);
                if (ret > 0) {
                    found = 1;
                }
                // PM1 is deterministic for a given B1, so no need to try 10 times with the same bound.
                if (method == METHOD_POLLARD_P1) break;
            }
        }
        
        ecm_clear(params);

        if (found) {
            if (mpz_cmp_ui(f, 1) == 0 || mpz_cmp(f, n) == 0) {
                if (method == METHOD_AUTOMATIC) factint_warn_incomplete(n);
                add_factor_mpz(factors, num_factors, n, 1);
                mpz_clear(f);
                return;
            }
            mpz_t n_f;
            mpz_init(n_f);
            mpz_divexact(n_f, n, f);

            factorize_mpz(f, factors, num_factors, k_limit, method, method_opt);
            factorize_mpz(n_f, factors, num_factors, k_limit, method, method_opt);
            mpz_clear(n_f);
        } else {
            // If it fails to find a factor within reasonable bounds, just add n as a factor
            if (method == METHOD_AUTOMATIC) factint_warn_incomplete(n);
            add_factor_mpz(factors, num_factors, n, 1);
        }
    }
    mpz_clear(f);
#else
    if (method == METHOD_ECM || method == METHOD_POLLARD_P1 || method == METHOD_WILLIAMS_P1) {
        add_factor_mpz(factors, num_factors, n, 1);
        mpz_clear(f);
        return;
    }
    if (method == METHOD_AUTOMATIC) {
        factint_warn_incomplete(n);
        add_factor_mpz(factors, num_factors, n, 1);
    }
    mpz_clear(f);
#endif
}

static int compare_factors_mpz(const void* a, const void* b) {
    const FactorMpz* fa = (const FactorMpz*)a;
    const FactorMpz* fb = (const FactorMpz*)b;
    return mpz_cmp(fa->p, fb->p);
}

Expr* builtin_factorinteger(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count < 1 || res->data.function.arg_count > 3)) return NULL;
    
    Expr* n_expr = res->data.function.args[0];
    mpz_t num_orig;
    mpz_init(num_orig);
    if (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_BIGINT) expr_to_mpz(n_expr, num_orig);
    mpz_t num, den;
    bool is_rat = false;

    if (n_expr->type == EXPR_INTEGER || n_expr->type == EXPR_BIGINT) {
        expr_to_mpz(n_expr, num);
        mpz_init_set_ui(den, 1);
        is_rat = true;
    } else {
        int64_t n_num, n_den;
        if (is_rational(n_expr, &n_num, &n_den)) {
            mpz_init_set_si(num, n_num);
            mpz_init_set_si(den, n_den);
            is_rat = true;
        }
    }

    if (!is_rat) return NULL;

    int method = METHOD_AUTOMATIC;
    Expr* method_opt = NULL;
    int k_limit = -1;

    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        Expr* arg = res->data.function.args[i];
        if (arg->type == EXPR_INTEGER) {
            k_limit = (int)arg->data.integer;
        } else if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_Automatic) {
            // keep default
        } else if ((arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && (arg->data.function.head->data.symbol.name == SYM_Rule || arg->data.function.head->data.symbol.name == SYM_RuleDelayed)) && arg->data.function.arg_count == 2) {
            Expr* lhs = arg->data.function.args[0];
            Expr* rhs = arg->data.function.args[1];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Method) {
                if (rhs->type == EXPR_STRING) {
                    if (strcmp(rhs->data.string, "TrialDivision") == 0) method = METHOD_TRIAL;
                    else if (strcmp(rhs->data.string, "PollardRho") == 0) method = METHOD_POLLARD_RHO;
                    else if (strcmp(rhs->data.string, "ECM") == 0) method = METHOD_ECM;
                    else if (strcmp(rhs->data.string, "Fermat") == 0) method = METHOD_FERMAT;
                    else if (strcmp(rhs->data.string, "CFRAC") == 0) method = METHOD_CFRAC;
                    else if (strcmp(rhs->data.string, "PollardP-1") == 0) method = METHOD_POLLARD_P1;
                    else if (strcmp(rhs->data.string, "WilliamsP+1") == 0) method = METHOD_WILLIAMS_P1;
                    else if (strcmp(rhs->data.string, "Dixon") == 0) method = METHOD_DIXON;
                    else if (strcmp(rhs->data.string, "ShanksSquareForms") == 0) method = METHOD_SQUFOF;
                    else if (strcmp(rhs->data.string, "BlakeRationalBaseDescent") == 0) method = METHOD_RBD;
                    else if (strcmp(rhs->data.string, "Automatic") == 0) method = METHOD_AUTOMATIC;
                } else if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
                    method = METHOD_AUTOMATIC;
                } else if (rhs->type == EXPR_FUNCTION && rhs->data.function.head->data.symbol.name == SYM_List) {
                    if (rhs->data.function.arg_count == 2) {
                        Expr* m_name = rhs->data.function.args[0];
                        Expr* m_opt = rhs->data.function.args[1];
                        if (m_name->type == EXPR_STRING && strcmp(m_name->data.string, "BlakeRationalBaseDescent") == 0) {
                            method = METHOD_RBD;
                            if (m_opt->type == EXPR_FUNCTION && m_opt->data.function.head->type == EXPR_SYMBOL && (m_opt->data.function.head->data.symbol.name == SYM_Rule || m_opt->data.function.head->data.symbol.name == SYM_RuleDelayed) && m_opt->data.function.arg_count == 2) {
                                Expr* opt_lhs = m_opt->data.function.args[0];
                                Expr* opt_rhs = m_opt->data.function.args[1];
                                if (opt_lhs->type == EXPR_STRING && strcmp(opt_lhs->data.string, "Base") == 0) {
                                    method_opt = expr_copy(opt_rhs);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            mpz_clear(num); mpz_clear(den); mpz_clear(num_orig);
            return NULL;
        }
    }


    FactorMpz factors[1024];
    int num_factors = 0;

    if (mpz_cmp_ui(num, 0) < 0) {
        mpz_t minus_one;
        mpz_init_set_si(minus_one, -1);
        add_factor_mpz(factors, &num_factors, minus_one, 1);
        mpz_clear(minus_one);
        mpz_neg(num, num);
    } else if (mpz_cmp_ui(num, 0) == 0) {
        mpz_t zero;
        mpz_init_set_ui(zero, 0);
        add_factor_mpz(factors, &num_factors, zero, 1);
        mpz_clear(zero);
        mpz_set_ui(num, 1);
    }

    // Since ECM extracts small factors first, Automatic doesn't need to be treated differently.
    factorize_mpz(num, factors, &num_factors, &k_limit, method, method_opt);
    if (method_opt) expr_free(method_opt);
    
    if (mpz_cmp_ui(den, 1) > 0) {
        FactorMpz d_factors[1024];
        int d_num_factors = 0;
        int d_limit = -1;
        factorize_mpz(den, d_factors, &d_num_factors, &d_limit, method, method_opt);
        for (int i = 0; i < d_num_factors; i++) {
            add_factor_mpz(factors, &num_factors, d_factors[i].p, -d_factors[i].count);
            mpz_clear(d_factors[i].p);
        }
    }


    qsort(factors, num_factors, sizeof(FactorMpz), compare_factors_mpz);

    Expr** list_args = malloc(sizeof(Expr*) * num_factors);
    for (int i = 0; i < num_factors; i++) {
        Expr** pair = malloc(sizeof(Expr*) * 2);
        pair[0] = expr_bigint_normalize(expr_new_bigint_from_mpz(factors[i].p));
        pair[1] = expr_new_integer(factors[i].count);
        list_args[i] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        free(pair);
        mpz_clear(factors[i].p);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), list_args, num_factors);
    free(list_args);
    mpz_clear(num); mpz_clear(den); mpz_clear(num_orig);
    return result;
}

Expr* builtin_eulerphi(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (!expr_is_integer_like(arg)) return NULL;

    mpz_t n;
    mpz_init(n);
    expr_to_mpz(arg, n);

    if (mpz_sgn(n) == 0) { mpz_clear(n); return expr_new_integer(0); }
    mpz_abs(n, n);                                          /* phi(-n) = phi(n) */
    if (mpz_cmp_ui(n, 1) == 0) { mpz_clear(n); return expr_new_integer(1); }

    /* Factor a working copy: factorize_mpz mutates its first argument. */
    mpz_t work;
    mpz_init_set(work, n);

    FactorMpz factors[1024];
    int num_factors = 0;
    int k_limit = -1;
    factorize_mpz(work, factors, &num_factors, &k_limit, METHOD_AUTOMATIC, NULL);
    mpz_clear(work);

    /* phi(n) = n * prod (1 - 1/p_i) = (n / prod p_i) * prod (p_i - 1).
     * Apply per-prime as (phi / p) * (p - 1) to keep intermediates exact
     * and avoid materializing prod p_i. */
    mpz_t phi, pm1;
    mpz_init_set(phi, n);
    mpz_init(pm1);
    for (int i = 0; i < num_factors; i++) {
        mpz_divexact(phi, phi, factors[i].p);
        mpz_sub_ui(pm1, factors[i].p, 1);
        mpz_mul(phi, phi, pm1);
        mpz_clear(factors[i].p);
    }
    mpz_clear(pm1);
    mpz_clear(n);

    Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(phi));
    mpz_clear(phi);
    return out;
}

void facint_init(void) {
    symtab_add_builtin("PrimeQ", builtin_primeq);
    symtab_add_builtin("FactorInteger", builtin_factorinteger);
    symtab_add_builtin("NextPrime", builtin_nextprime);
    symtab_add_builtin("EulerPhi", builtin_eulerphi);

    symtab_get_def("PrimeQ")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    /* PrimePi is registered in numbertheory_init() (src/numbertheory/prime.c). */
    symtab_get_def("FactorInteger")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("EulerPhi")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_get_def("NextPrime")->attributes |= (ATTR_PROTECTED | ATTR_READPROTECTED | ATTR_LISTABLE);
}
