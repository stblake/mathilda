#include "eval.h"
#include "power.h"
#include "arithmetic.h"
#include "times.h"
#include "numeric.h"
#include "common.h"
#include "sym_names.h"
#include "trig_canon.h"
#include "internal.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <gmp.h>

/* Strict-positive constant predicate. Returns true for positive integer/
 * rational/real, Pi, E, and Power[B, integer] where B is one of those
 * (covers Sqrt[2], Pi^2, etc. as factors). Used by the rational-exponent
 * Times-base distribution to gate when the base may be cleanly split
 * into "imaginary part" and "positive part" pieces. */
static bool is_known_positive_pwr(Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) > 0;
    if (e->type == EXPR_REAL)    return e->data.real > 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return n > 0;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol;
        if (s == SYM_E)  return true;
        if (s == SYM_Pi) return true;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        return is_known_positive_pwr(e->data.function.args[0]);
    }
    return false;
}

static int64_t ipow(int64_t base, int64_t exp, bool* overflow) {
    if (exp < 0) return 0;
    if (exp == 0) return 1;
    if (base == 0) return 0;
    if (base == 1) return 1;
    if (base == -1) return (exp % 2 == 0) ? 1 : -1;

    __int128_t res = 1;
    __int128_t b = base;
    while (exp > 0) {
        if (exp % 2 == 1) {
            res *= b;
            if (res > INT64_MAX || res < INT64_MIN) {
                *overflow = true;
                return 0;
            }
        }
        exp /= 2;
        if (exp > 0) {
            if (b > INT64_MAX || b < INT64_MIN) {
                *overflow = true;
                return 0;
            }
            b *= b;
        }
    }
    return (int64_t)res;
}

static void factor_out_kth_power(int64_t n, int64_t k, int64_t* out_m, int64_t* out_r) {
    int64_t m = 1;
    int64_t r = n;
    int64_t d = 2;
    while (d * d <= r) {
        bool ov = false;
        int64_t dk = ipow(d, k, &ov);
        if (ov || dk <= 0 || dk > r) {
            if (ov || dk > r) break;
            d++;
            continue;
        }
        while (r % dk == 0) {
            m *= d;
            r /= dk;
        }
        d++;
    }
    *out_m = m;
    *out_r = r;
}

/* Trial-divide |n| into ascending distinct prime factors. Parallel
 * arrays primes[] / exps[] get the factorisation; *count is the number
 * of distinct primes. Returns false (with partial fill discarded) when
 * more than max_primes distinct factors appear -- caller treats this
 * as "give up, leave the input alone". Cost is O(sqrt(n)) trial
 * divisions, matching factor_out_kth_power's loop bound and intended
 * for the int64-bounded residues that path produces. */
static bool power_factor_int64(int64_t n,
                               int64_t* primes, int64_t* exps,
                               int max_primes, int* count) {
    *count = 0;
    if (n < 0) n = -n;
    if (n <= 1) return true;
    int64_t d = 2;
    while (d <= n / d) {
        if (n % d == 0) {
            if (*count >= max_primes) return false;
            primes[*count] = d;
            int64_t e = 0;
            while (n % d == 0) { n /= d; e++; }
            exps[*count] = e;
            (*count)++;
        }
        d++;
    }
    if (n > 1) {
        if (*count >= max_primes) return false;
        primes[*count] = n;
        exps[*count] = 1;
        (*count)++;
    }
    return true;
}

/* Canonicalise Power[r, b_rem/q] for positive integer r, integer
 * b_rem > 0 and integer q > 1 by splitting r into its prime
 * factorisation and grouping primes by effective exponent.
 *
 * Algorithm: factor r = ∏ p_i^a_i. For each prime, compute
 *   num_i = a_i * b_rem,   int_i = num_i / q,   d_i = num_i mod q
 * (0 <= d_i < q). The coefficient absorbs ∏ p_i^int_i; the radical
 * part groups primes by their reduced (d_i/q) fraction and emits one
 * Power[base_group, d_i/q] per group, where base_group is the product
 * of primes in that group.
 *
 * Triggers only when splitting strictly adds information --
 * i.e. some int_i != 0 (a perfect q-th power inside r contributes a
 * rational coefficient) or the reduced (d_i, q) fractions are not all
 * identical (the input r^(b_rem/q) compactly hides differing per-prime
 * exponents). For uniform exponents with no integer part (e.g.
 * 6^(1/3), 30^(1/3)), returns NULL so the caller keeps the unsplit
 * Power[r, b_rem/q] form -- matching Mathematica's preference for
 * the shorter representation when no information is gained.
 *
 * Examples:
 *   r=18, b_rem=1, q=3 -> 2^(1/3) * 3^(2/3)
 *   r=12, b_rem=1, q=3 -> 2^(2/3) * 3^(1/3)
 *   r=60, b_rem=1, q=3 -> 2^(2/3) * 15^(1/3)   (3,5 share eff 1/3)
 *   r=50, b_rem=2, q=3 -> 5 * 2^(2/3) * 5^(1/3) (5^(4/3) splits int*)
 *   r=6,  b_rem=1, q=3 -> NULL                  (all effs equal)
 *   r=30, b_rem=1, q=3 -> NULL                  (all effs equal)
 *
 * Restricted to b_rem > 0; the negative-numerator case stays in the
 * existing Power[res_base, b_rem/q] form (see the integer-rat block
 * below) to avoid having to reason about (-1)^(int_i) sign aggregation
 * for negative base contexts. */
static Expr* power_split_residue(int64_t r, int64_t b_rem, int64_t q) {
    if (r <= 1 || b_rem <= 0 || q <= 1) return NULL;
    enum { MAX_P = 20 };
    int64_t primes[MAX_P], exps[MAX_P];
    int npr = 0;
    if (!power_factor_int64(r, primes, exps, MAX_P, &npr)) return NULL;
    if (npr <= 1) return NULL;

    int64_t int_part[MAX_P], eff_num[MAX_P], eff_den[MAX_P];
    for (int i = 0; i < npr; i++) {
        __int128_t prod = (__int128_t)exps[i] * (__int128_t)b_rem;
        if (prod > INT64_MAX) return NULL;
        int64_t num = (int64_t)prod;
        int_part[i] = num / q;
        int64_t rem = num - int_part[i] * q;
        if (rem == 0) {
            eff_num[i] = 0;
            eff_den[i] = 1;
        } else {
            int64_t g = gcd(rem, q);
            eff_num[i] = rem / g;
            eff_den[i] = q / g;
        }
    }

    bool any_int = false, uniform_eff = true;
    for (int i = 0; i < npr; i++) {
        if (int_part[i] != 0) any_int = true;
        if (eff_num[i] != eff_num[0] || eff_den[i] != eff_den[0]) uniform_eff = false;
    }
    if (!any_int && uniform_eff) return NULL;

    /* Coefficient = ∏ p_i^int_part[i] in GMP for bigint promotion. */
    mpz_t coeff_z;
    mpz_init_set_ui(coeff_z, 1);
    for (int i = 0; i < npr; i++) {
        if (int_part[i] > 0) {
            mpz_t t; mpz_init_set_si(t, primes[i]);
            mpz_pow_ui(t, t, (unsigned long)int_part[i]);
            mpz_mul(coeff_z, coeff_z, t);
            mpz_clear(t);
        }
    }

    Expr** out = malloc(sizeof(Expr*) * (npr + 1));
    int out_count = 0;
    if (mpz_cmp_ui(coeff_z, 1) != 0) {
        out[out_count++] = expr_bigint_normalize(expr_new_bigint_from_mpz(coeff_z));
    }
    mpz_clear(coeff_z);

    bool used[MAX_P] = {false};
    for (int i = 0; i < npr; i++) {
        if (used[i]) continue;
        used[i] = true;
        if (eff_num[i] == 0) continue;
        mpz_t base_z;
        mpz_init_set_si(base_z, primes[i]);
        for (int j = i + 1; j < npr; j++) {
            if (!used[j] && eff_num[j] == eff_num[i] && eff_den[j] == eff_den[i]) {
                mpz_mul_si(base_z, base_z, primes[j]);
                used[j] = true;
            }
        }
        Expr* base_e = expr_bigint_normalize(expr_new_bigint_from_mpz(base_z));
        mpz_clear(base_z);
        Expr* exp_e = make_rational(eff_num[i], eff_den[i]);
        out[out_count++] = expr_new_function(expr_new_symbol("Power"),
                                             (Expr*[]){base_e, exp_e}, 2);
    }

    if (out_count == 0) { free(out); return expr_new_integer(1); }
    if (out_count == 1) { Expr* r_ = out[0]; free(out); return r_; }
    Expr* result = expr_new_function(expr_new_symbol("Times"), out, out_count);
    free(out);
    return result;
}

/* Returns true when Power[f, p/q] is guaranteed to evaluate to a pure
 * rational coefficient with no remaining q-th-power radical. Used to
 * gate distribution of Power[Times[positives], p/q]: we distribute only
 * when at least one factor cleanly extracts to a rational (the
 * remaining radicals on other factors then rejoin canonically). This
 * matches Mathematica's behaviour for Sqrt[4 Pi] -> 2 Sqrt[Pi] (4 is a
 * perfect square) while leaving (4 Pi)^(2/3) alone (4 has no cube
 * factor). The predicate intentionally rejects Power[B, E] factors --
 * even though composition always works, distributing on Power-only
 * factors with no clean rational extraction would split forms like
 * Sqrt[2^a 3^a] into Sqrt[2^a] Sqrt[3^a], which is neither shorter nor
 * more canonical. */
static bool factor_fully_reduces_under_q(Expr* f, int64_t q) {
    if (f->type == EXPR_INTEGER) {
        int64_t n = f->data.integer;
        if (n <= 0) return false;
        if (n == 1) return true;
        int64_t m, r;
        factor_out_kth_power(n, q, &m, &r);
        return (r == 1);
    }
    int64_t n, d;
    if (is_rational(f, &n, &d)) {
        if (n <= 0) return false;
        int64_t m, r;
        factor_out_kth_power(n, q, &m, &r);
        if (r != 1) return false;
        factor_out_kth_power(d, q, &m, &r);
        return (r == 1);
    }
    return false;
}

/* Companion gate for Power[positive_base, p'/q'] factors inside
 * Power[Times[positives], p_out/q_out].  Triggers distribution when
 * composing the outer exponent with the inner exponent introduces no
 * new radical complexity, i.e. the composed exponent's reduced
 * denominator is at most q' (the inner radical's denominator).
 *
 * Composed exponent = (p' * p_out) / (q' * q_out).  After gcd
 * reduction the denominator is (q' * q_out) / gcd(|p' * p_out|,
 * q' * q_out).  Bounded above by q' iff q_out divides |p' * p_out|.
 *
 * Examples (assuming positive base):
 *   Sqrt[Power[2, -2/3]]   -> q_out=2, p_out=1, p'=-2, q'=3
 *     |p' * p_out| = 2, 2 mod 2 = 0 -> composes (-> Power[2, -1/3])
 *   Sqrt[Power[3, 1/2]]    -> q_out=2, p_out=1, p'=1,  q'=2
 *     |p' * p_out| = 1, 1 mod 2 = 1 -> rejects (would introduce 4th root)
 *
 * This lets `Sqrt[k/2^(2/3)]` distribute for any positive `k` (e.g.
 * `k = 1/3` or `k = Pi`), via the clean composition of the Power
 * factor, while preserving Mathematica-faithful behaviour of leaving
 * `Sqrt[2 Sqrt[3]]` and `Sqrt[2 Pi]` alone. */
static bool power_factor_composes_cleanly(Expr* f,
                                          int64_t p_out, int64_t q_out) {
    if (f->type != EXPR_FUNCTION) return false;
    if (!f->data.function.head
        || f->data.function.head->type != EXPR_SYMBOL
        || f->data.function.head->data.symbol != SYM_Power
        || f->data.function.arg_count != 2) return false;
    Expr* b = f->data.function.args[0];
    Expr* e = f->data.function.args[1];
    if (!is_known_positive_pwr(b)) return false;
    int64_t pi, qi;
    if (!is_rational(e, &pi, &qi)) return false;
    if (qi <= 1) return false;
    int64_t abs_p  = pi < 0 ? -pi : pi;
    int64_t abs_po = p_out < 0 ? -p_out : p_out;
    /* Overflow-safe |p' * p_out|; bail to false on overflow. */
    if (abs_p != 0 && abs_po != 0 && abs_po > INT64_MAX / abs_p) return false;
    int64_t prod = abs_p * abs_po;
    return q_out > 0 && (prod % q_out == 0);
}

/* Given positive n_in > 1, decide whether n_in = b^k for some integer
 * b >= 2 and k >= 2. On true, fills b_out with the smallest such b
 * (equivalently the largest such k) and *k_out with that k.
 *
 * Algorithm: mpz_perfect_power_p gives the fast yes/no, then we
 * iteratively factor a prime out of the exponent by checking
 * mpz_root for small primes. Each successful iteration shrinks
 * mpz_sizeinbase(cur, 2) by a factor of >= p, so total work is
 * O(log(n) * log^2(n)) bit operations -- well below the cost of
 * the radical extraction below. Restricted to positive n: negative
 * perfect powers (like -8 = (-2)^3) are NOT unified because
 * (-n)^(p/q) ≠ ((-n)^(1/k))^(kp/q) in general on the principal branch
 * (e.g., (-8)^(1/6) is not (-2)^(1/2)).
 *
 * The prime table covers k_total up to 127 in a single factor; for
 * inputs requiring larger primes (mpz_perfect_power_p says yes but
 * none of these primes works) the function returns whatever k_total
 * was accumulated, which may be less than maximal. This is sound but
 * occasionally suboptimal; in practice every k_total reachable from
 * a 64-bit base is well below 64. */
static bool find_min_perfect_base(const mpz_t n_in, mpz_t b_out, int64_t* k_out) {
    if (mpz_sgn(n_in) <= 0) return false;
    if (mpz_cmp_ui(n_in, 1) <= 0) return false;
    if (!mpz_perfect_power_p(n_in)) return false;

    static const unsigned long primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
        101,103,107,109,113,127
    };
    const size_t n_primes = sizeof(primes) / sizeof(primes[0]);

    mpz_t cur, root;
    mpz_init_set(cur, n_in);
    mpz_init(root);
    int64_t k_total = 1;

    bool reduced = true;
    while (reduced) {
        reduced = false;
        for (size_t i = 0; i < n_primes; i++) {
            unsigned long p = primes[i];
            /* If cur < 2^p then cur cannot be a p-th power of an int >= 2. */
            if (mpz_sizeinbase(cur, 2) <= p) break;
            if (mpz_root(root, cur, p)) {
                /* Guard against k_total overflow (only triggered by absurdly
                 * large exponents -- for safety, stop and return current). */
                if (k_total > INT64_MAX / (int64_t)p) break;
                mpz_set(cur, root);
                k_total *= (int64_t)p;
                reduced = true;
                break;
            }
        }
    }

    bool found = (k_total >= 2);
    if (found) { mpz_set(b_out, cur); *k_out = k_total; }
    mpz_clear(cur);
    mpz_clear(root);
    return found;
}

static Expr* bigint_pow(const Expr* base, int64_t exp) {
    if (exp < 0) return NULL;
    mpz_t b, r;
    expr_to_mpz(base, b);
    mpz_init(r);
    mpz_pow_ui(r, b, (unsigned long)exp);
    mpz_clear(b);
    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
    mpz_clear(r);
    return result;
}

Expr* make_power(Expr* base, Expr* exp) {
    Expr* args[2] = { base, exp };
    return expr_new_function(expr_new_symbol("Power"), args, 2);
}

static bool is_head_call(Expr* e, const char* sym, size_t argc) {
    return head_is(e, sym) && e->data.function.arg_count == argc;
}

/* Cancel Log inside a Power exponent:
 *   Exp[c Log[a]]                  -> a^c                   (base == E)
 *   Power[base, c Log[base, a]]    -> a^c                   (general)
 *
 * Internally Log[b, a] is represented as Log[a] * Log[b]^(-1), so case 2
 * reduces to matching Log[a] * Power[Log[base], -1] among the exponent's
 * Times-factors. We require exactly one Log[a] factor and, when base != E,
 * exactly one matching Power[Log[base], -1] factor; the remaining factors
 * become the coefficient c in the rewritten Power[a, c]. */
static Expr* simplify_exp_log(Expr* base, Expr* exp) {
    bool exp_is_times = exp->type == EXPR_FUNCTION &&
                        exp->data.function.head->type == EXPR_SYMBOL &&
                        exp->data.function.head->data.symbol == SYM_Times;
    size_t nf = exp_is_times ? exp->data.function.arg_count : 1;
    Expr** factors = exp_is_times ? exp->data.function.args : &exp;

    bool base_is_E = (base->type == EXPR_SYMBOL && base->data.symbol == SYM_E);

    int log_idx = -1;
    Expr* a = NULL;
    for (size_t i = 0; i < nf; i++) {
        if (!is_head_call(factors[i], SYM_Log, 1)) continue;
        log_idx = (int)i;
        a = factors[i]->data.function.args[0];
        break;
    }
    if (log_idx < 0) return NULL;

    int inv_log_idx = -1;
    if (!base_is_E) {
        for (size_t i = 0; i < nf; i++) {
            if ((int)i == log_idx) continue;
            Expr* f = factors[i];
            if (!is_head_call(f, SYM_Power, 2)) continue;
            Expr* ib = f->data.function.args[0];
            Expr* ie = f->data.function.args[1];
            if (!(ie->type == EXPR_INTEGER && ie->data.integer == -1)) continue;
            if (!is_head_call(ib, SYM_Log, 1)) continue;
            if (!expr_eq(ib->data.function.args[0], base)) continue;
            inv_log_idx = (int)i;
            break;
        }
        if (inv_log_idx < 0) return NULL;
    }

    Expr* coeff;
    if (!exp_is_times) {
        coeff = expr_new_integer(1);
    } else {
        Expr** rest = malloc(sizeof(Expr*) * nf);
        size_t kept = 0;
        for (size_t i = 0; i < nf; i++) {
            if ((int)i == log_idx) continue;
            if ((int)i == inv_log_idx) continue;
            rest[kept++] = expr_copy(factors[i]);
        }
        if (kept == 0) {
            free(rest);
            coeff = expr_new_integer(1);
        } else if (kept == 1) {
            coeff = rest[0];
            free(rest);
        } else {
            coeff = expr_new_function(expr_new_symbol("Times"), rest, kept);
            free(rest);
        }
    }

    return eval_and_free(expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(a), coeff }, 2));
}

Expr* builtin_power(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t n_args = res->data.function.arg_count;
    if (n_args == 0) return NULL;
    if (n_args == 1) return expr_copy(res->data.function.args[0]);

    if (n_args > 2) {
        Expr** sub_args = malloc(sizeof(Expr*) * (n_args - 1));
        for (size_t i = 0; i < n_args - 1; i++) {
            sub_args[i] = expr_copy(res->data.function.args[i+1]);
        }
        Expr* sub_power = expr_new_function(expr_new_symbol("Power"), sub_args, n_args - 1);
        free(sub_args);
        
        Expr* final_args[2] = { expr_copy(res->data.function.args[0]), sub_power };
        return expr_new_function(expr_new_symbol("Power"), final_args, 2);
    }

    if (n_args == 1) return NULL;

    Expr* base = res->data.function.args[0];
    Expr* exp = res->data.function.args[1];

    /* Infinity / Indeterminate preprocessing.
     *
     * Mathematica semantics for Power (selected):
     *   Indeterminate ^ x, x ^ Indeterminate          -> Indeterminate
     *   1 ^ Infinity, 1 ^ -Infinity, 1 ^ ComplexInfinity -> Indeterminate (msg)
     *   Infinity ^ 0, (-Infinity) ^ 0, ComplexInfinity ^ 0 -> Indeterminate (msg)
     *   0 ^ Infinity                                   -> 0
     *   0 ^ -Infinity                                  -> ComplexInfinity (msg)
     *   0 ^ ComplexInfinity                            -> Indeterminate (msg)
     *   Infinity ^ Infinity                            -> ComplexInfinity
     *   Infinity ^ -Infinity                           -> 0
     *   Infinity ^ n   (n numeric, n > 0)              -> Infinity
     *   Infinity ^ n   (n numeric, n < 0)              -> 0
     *   ComplexInfinity ^ n (n numeric, n > 0)         -> ComplexInfinity
     *   ComplexInfinity ^ n (n numeric, n < 0)         -> 0
     *
     * Note: 0^positive and 0^negative for ordinary exponents are still handled
     * by the existing logic below; only the infinity-flavoured cases are
     * intercepted here.
     */
    {
        bool b_indet = is_indeterminate_sym(base);
        bool e_indet = is_indeterminate_sym(exp);
        if (b_indet || e_indet) return expr_new_symbol("Indeterminate");

        bool b_inf  = is_infinity_sym(base);
        bool b_ninf = is_neg_infinity_form(base);
        bool b_cinf = is_complex_infinity_sym(base);
        bool e_inf  = is_infinity_sym(exp);
        bool e_ninf = is_neg_infinity_form(exp);
        bool e_cinf = is_complex_infinity_sym(exp);

        bool base_is_one  = (base->type == EXPR_INTEGER && base->data.integer == 1);
        bool base_is_zero_lit = (base->type == EXPR_INTEGER && base->data.integer == 0) ||
                                (base->type == EXPR_REAL && base->data.real == 0.0) ||
                                (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) == 0);
        bool exp_is_zero_lit  = (exp->type == EXPR_INTEGER && exp->data.integer == 0) ||
                                (exp->type == EXPR_REAL && exp->data.real == 0.0) ||
                                (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) == 0);

        bool base_is_inf = b_inf || b_ninf || b_cinf;
        bool exp_is_inf  = e_inf || e_ninf || e_cinf;

        if (base_is_one && exp_is_inf) {
            const char* what = e_inf ? "Infinity" : (e_ninf ? "-Infinity" : "ComplexInfinity");
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression 1^%s encountered.\n", what);
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_inf && exp_is_zero_lit) {
            const char* what = b_inf ? "Infinity" : (b_ninf ? "-Infinity" : "ComplexInfinity");
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression %s^0 encountered.\n", what);
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_zero_lit && e_cinf) {
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression 0^ComplexInfinity encountered.\n");
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_zero_lit && e_inf) {
            return expr_new_integer(0);
        }

        if (base_is_zero_lit && e_ninf) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            return expr_new_symbol("ComplexInfinity");
        }

        if (b_inf) {
            if (e_inf)  return expr_new_symbol("ComplexInfinity");
            if (e_ninf) return expr_new_integer(0);
            if (e_cinf) return expr_new_symbol("Indeterminate");
            int es = expr_numeric_sign(exp);
            if (es > 0) return expr_new_symbol("Infinity");
            if (es < 0) return expr_new_integer(0);
        }

        if (b_cinf) {
            if (e_inf)  return expr_new_symbol("ComplexInfinity");
            if (e_ninf) return expr_new_integer(0);
            int es = expr_numeric_sign(exp);
            if (es > 0) return expr_new_symbol("ComplexInfinity");
            if (es < 0) return expr_new_integer(0);
        }
    }

    /* 0^0 is indeterminate. Check before the generic Power[_,0] -> 1
     * identity below so the special case wins. Catches integer 0, real 0.0,
     * and Power[0, 0.0] / Power[0.0, 0] mixed forms. */
    {
        bool b_zero_lit = (base->type == EXPR_INTEGER && base->data.integer == 0)
                       || (base->type == EXPR_REAL    && base->data.real == 0.0);
        bool e_zero_lit = (exp->type == EXPR_INTEGER && exp->data.integer == 0)
                       || (exp->type == EXPR_REAL    && exp->data.real == 0.0);
        if (b_zero_lit && e_zero_lit) {
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Power::indet: Indeterminate expression 0^0 encountered.\n");
            return expr_new_symbol("Indeterminate");
        }
    }

    if (exp->type == EXPR_INTEGER && exp->data.integer == 0) return expr_new_integer(1);
    if (exp->type == EXPR_INTEGER && exp->data.integer == 1) return expr_copy(base);
    if (base->type == EXPR_INTEGER && base->data.integer == 1) return expr_new_integer(1);

    bool base_is_zero = false;
    if (base->type == EXPR_INTEGER && base->data.integer == 0) base_is_zero = true;
    if (base->type == EXPR_REAL && base->data.real == 0.0) base_is_zero = true;

    bool exp_is_negative = false;
    bool exp_is_positive = false;
    if (exp->type == EXPR_INTEGER) {
        if (exp->data.integer < 0) exp_is_negative = true;
        else if (exp->data.integer > 0) exp_is_positive = true;
    }
    if (exp->type == EXPR_REAL) {
        if (exp->data.real < 0.0) exp_is_negative = true;
        else if (exp->data.real > 0.0) exp_is_positive = true;
    }
    int64_t rn, rd;
    if (is_rational(exp, &rn, &rd)) {
        if (rn < 0) exp_is_negative = true;
        else if (rn > 0) exp_is_positive = true;
    }

    if (base_is_zero && exp_is_negative) {
        if (!arith_warnings_muted())
            fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
        return expr_new_symbol("ComplexInfinity");
    }
    /* 0^positive -> 0 (includes 0^(1/2) = Sqrt[0] = 0, 0^0.5 = 0, etc.).
     * Preserves the base's numeric type: real 0.0 for real exponents or a
     * real-base zero, integer 0 otherwise. */
    if (base_is_zero && exp_is_positive) {
        if (base->type == EXPR_REAL || exp->type == EXPR_REAL) return expr_new_real(0.0);
        return expr_new_integer(0);
    }

    Expr *re_b = NULL, *im_b = NULL, *re_e = NULL, *im_e = NULL;
    bool base_comp = is_complex(base, &re_b, &im_b);
    bool exp_comp = is_complex(exp, &re_e, &im_e);

    /* Symbolic principal-value evaluation for Power[Complex[0, ±1], n/2]
     * (n in {±1, ±3} after mod-4 reduction).  Without these, Mathilda
     * leaves Sqrt[I], Sqrt[-I], (I)^(3/2), (-I)^(3/2) as opaque held
     * Powers, which propagate through intrat output as
     *   Log[Sqrt[I] + x] / (I)^(3/2)
     * and downstream Simplify cannot canonicalise.  The principal
     * values are
     *   Sqrt[I]      = (1 + I) / Sqrt[2]
     *   Sqrt[-I]     = (1 - I) / Sqrt[2]
     *   (I)^(3/2)    = (-1 + I) / Sqrt[2]
     *   (-I)^(3/2)   = (-1 - I) / Sqrt[2]
     * (taking the standard principal branch arg = atan2(im, re)). */
    if (base_comp && !exp_comp
        && re_b->type == EXPR_INTEGER && re_b->data.integer == 0
        && im_b->type == EXPR_INTEGER
        && (im_b->data.integer == 1 || im_b->data.integer == -1)) {
        int64_t en = 0, ed = 1;
        bool exp_is_rational =
            (exp->type == EXPR_INTEGER) ? (en = exp->data.integer, ed = 1, true)
            : is_rational(exp, &en, &ed);
        if (exp_is_rational && ed == 2 && im_b->data.integer == 1) {
            /* I^(n/2):  reduce n mod 4. n=1: (1+I)/Sqrt[2]; n=3: (-1+I)/Sqrt[2];
             * n=-1: (1-I)/Sqrt[2]; n=-3: (-1-I)/Sqrt[2].  Even n folds to
             * I^(n/2) = (I^n)^(1/1) which the integer path already
             * handles. */
            int64_t r = ((en % 4) + 4) % 4;  /* normalise to {0, 1, 2, 3} */
            if (r == 1 || r == 3) {
                Expr* sgn_re = expr_new_integer((r == 1) ? 1 : -1);
                Expr* sgn_im = expr_new_integer(1);
                Expr* sqrt2_inv = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){expr_new_integer(2),
                              expr_new_function(expr_new_symbol("Rational"),
                                  (Expr*[]){expr_new_integer(-1),
                                            expr_new_integer(2)}, 2)}, 2);
                sqrt2_inv = eval_and_free(sqrt2_inv);
                Expr* z = make_complex(sgn_re, sgn_im);
                return eval_and_free(internal_times(
                    (Expr*[]){z, sqrt2_inv}, 2));
            }
        }
        if (exp_is_rational && ed == 2 && im_b->data.integer == -1) {
            /* (-I)^(n/2) = conjugate(I^(n/2)) for real coefficients —
             * mirrors the I case with imaginary sign flipped. */
            int64_t r = ((en % 4) + 4) % 4;
            if (r == 1 || r == 3) {
                Expr* sgn_re = expr_new_integer((r == 1) ? 1 : -1);
                Expr* sgn_im = expr_new_integer(-1);
                Expr* sqrt2_inv = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){expr_new_integer(2),
                              expr_new_function(expr_new_symbol("Rational"),
                                  (Expr*[]){expr_new_integer(-1),
                                            expr_new_integer(2)}, 2)}, 2);
                sqrt2_inv = eval_and_free(sqrt2_inv);
                Expr* z = make_complex(sgn_re, sgn_im);
                return eval_and_free(internal_times(
                    (Expr*[]){z, sqrt2_inv}, 2));
            }
        }
    }

    if (base_comp || exp_comp ||
        base->type == EXPR_REAL || exp->type == EXPR_REAL ||
        base->type == EXPR_BIGINT || exp->type == EXPR_BIGINT
#ifdef USE_MPFR
        || base->type == EXPR_MPFR || exp->type == EXPR_MPFR
#endif
        ) {
        bool base_num = expr_is_numeric_like(base) || base_comp;
        bool exp_num  = expr_is_numeric_like(exp)  || exp_comp;

        if (base_num && exp_num) {
#ifdef USE_MPFR
            /* MPFR real-valued power: fast path when both sides are real
             * and at least one carries MPFR precision. */
            if (!base_comp && !exp_comp && numeric_any_mpfr(base, exp)) {
                /* Only take the real MPFR path when the result is real:
                 * base > 0 OR exponent is an integer. */
                Expr *b_re_check = base;
                bool base_nonneg = (b_re_check->type == EXPR_INTEGER && b_re_check->data.integer >= 0)
                                   || (b_re_check->type == EXPR_REAL && b_re_check->data.real >= 0.0)
                                   || (b_re_check->type == EXPR_BIGINT && mpz_sgn(b_re_check->data.bigint) >= 0)
                                   || (b_re_check->type == EXPR_MPFR && mpfr_sgn(b_re_check->data.mpfr) >= 0);
                int64_t en, ed;
                bool exp_is_int = (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT
                                   || (is_rational(exp, &en, &ed) && ed == 1));
                if (base_nonneg || exp_is_int) {
                    Expr* r = numeric_mpfr_pow(base, exp, 0);
                    if (r) return r;
                }
            }
            /* MPFR complex-valued power: handles the cases that the real
             * MPFR path skipped — negative base with fractional exponent
             * (e.g. (-1)^(1/4)), complex base, or complex exponent
             * (e.g. E^(I Pi/4)). Without this the cpow fallback below
             * would coerce MPFR operands to zero and yield NaN. */
            if (numeric_any_mpfr(base, exp)) {
                Expr* r = numeric_mpfr_complex_pow(base, exp, 0);
                if (r) return r;
            }
#endif
            bool has_real = (base->type == EXPR_REAL || exp->type == EXPR_REAL ||
                            (base_comp && (re_b->type == EXPR_REAL || im_b->type == EXPR_REAL)) ||
                            (exp_comp && (re_e->type == EXPR_REAL || im_e->type == EXPR_REAL)));
#ifdef USE_MPFR
            if (!has_real) {
                /* Treat an MPFR operand as "real" for the double-complex
                 * fallback path too. */
                has_real = (base->type == EXPR_MPFR || exp->type == EXPR_MPFR);
            }
#endif
            if (has_real) {
                double vbase_re = 0, vbase_im = 0, vexp_re = 0, vexp_im = 0;
                int64_t n, d;
                
                if (base_comp) {
                    vbase_re = (re_b->type == EXPR_REAL) ? re_b->data.real : (re_b->type == EXPR_INTEGER) ? (double)re_b->data.integer : (re_b->type == EXPR_BIGINT) ? mpz_get_d(re_b->data.bigint) : 0;
                    vbase_im = (im_b->type == EXPR_REAL) ? im_b->data.real : (im_b->type == EXPR_INTEGER) ? (double)im_b->data.integer : (im_b->type == EXPR_BIGINT) ? mpz_get_d(im_b->data.bigint) : 0;
                    if (is_rational(re_b, &n, &d)) vbase_re = (double)n / d;
                    if (is_rational(im_b, &n, &d)) vbase_im = (double)n / d;
                } else {
                    vbase_re = (base->type == EXPR_REAL) ? base->data.real : (base->type == EXPR_INTEGER) ? (double)base->data.integer : (base->type == EXPR_BIGINT) ? mpz_get_d(base->data.bigint) : 0;
                    if (is_rational(base, &n, &d)) vbase_re = (double)n / d;
                }
                
                if (exp_comp) {
                    vexp_re = (re_e->type == EXPR_REAL) ? re_e->data.real : (re_e->type == EXPR_INTEGER) ? (double)re_e->data.integer : (re_e->type == EXPR_BIGINT) ? mpz_get_d(re_e->data.bigint) : 0;
                    vexp_im = (im_e->type == EXPR_REAL) ? im_e->data.real : (im_e->type == EXPR_INTEGER) ? (double)im_e->data.integer : (im_e->type == EXPR_BIGINT) ? mpz_get_d(im_e->data.bigint) : 0;
                    if (is_rational(re_e, &n, &d)) vexp_re = (double)n / d;
                    if (is_rational(im_e, &n, &d)) vexp_im = (double)n / d;
                } else {
                    vexp_re = (exp->type == EXPR_REAL) ? exp->data.real : (exp->type == EXPR_INTEGER) ? (double)exp->data.integer : (exp->type == EXPR_BIGINT) ? mpz_get_d(exp->data.bigint) : 0;
                    if (is_rational(exp, &n, &d)) vexp_re = (double)n / d;
                }

                double complex z = vbase_re + vbase_im * I;
                double complex w = vexp_re + vexp_im * I;
                
                if (vbase_re < 0.0 && vbase_im == 0.0 && vexp_im == 0.0 && floor(vexp_re) != vexp_re) {
                    double r = pow(-vbase_re, vexp_re);
                    double theta = vexp_re * 3.14159265358979323846;
                    Expr* c_re = expr_new_real(r * cos(theta));
                    Expr* c_im = expr_new_real(r * sin(theta));
                    Expr* comp = make_complex(c_re, c_im);
                    return comp;
                }
                
                double complex result = cpow(z, w);
                double r_re = creal(result);
                double r_im = cimag(result);
                
                if (r_im == 0.0 || (fabs(r_im) < 1e-14 && fabs(r_re) > 1e-14)) {
                    return expr_new_real(r_re);
                }
                
                Expr* c_re = expr_new_real(r_re);
                Expr* c_im = expr_new_real(r_im);
                Expr* comp = make_complex(c_re, c_im);
                return comp;
            }
        }
    }
    
    if (base_comp && exp->type == EXPR_INTEGER) {
        int64_t e = exp->data.integer;
        if (e < 0 && e > -1000) {
            // e < 0, compute (base^-e)^-1
            Expr* pos_pow = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), expr_new_integer(-e)}, 2);
            Expr* denom = evaluate(pos_pow);
            expr_free(pos_pow);
            
            // Now we have denom = a + b I. We need 1 / denom.
            // We can emit: Divide[1, denom], which is: Times[1, Power[denom, -1]] -> WAIT that loops!
            // No, we need to construct the conjugate division manually!
            Expr* re; Expr* im;
            if (is_complex(denom, &re, &im)) {
                // conj = a - b I
                Expr* conj = make_complex(expr_copy(re), expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(im)}, 2));
                // mag_sq = a^2 + b^2
                Expr* a2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(re), expr_new_integer(2)}, 2);
                Expr* b2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(im), expr_new_integer(2)}, 2);
                Expr* mag_sq = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){a2, b2}, 2);
                
                // Result = conj / mag_sq
                Expr* res_ast = expr_new_function(expr_new_symbol("Divide"), (Expr*[]){conj, mag_sq}, 2);
                Expr* result = evaluate(res_ast);
                expr_free(res_ast);
                expr_free(denom);
                return result;
            } else {
                // If it evaluated to Real or Integer, just do Power[denom, -1] natively
                Expr* res_ast = expr_new_function(expr_new_symbol("Divide"), (Expr*[]){expr_new_integer(1), expr_copy(denom)}, 2);
                Expr* result = evaluate(res_ast);
                expr_free(res_ast);
                expr_free(denom);
                return result;
            }
        }
        if (e >= 0 && e < 1000) { 
            Expr** prod_args = malloc(sizeof(Expr*) * (size_t)e);
            for(int64_t i=0; i<e; i++) prod_args[i] = expr_copy(base);
            Expr* prod = expr_new_function(expr_new_symbol("Times"), prod_args, (size_t)e);
            free(prod_args);
            return prod;
        }
    }

    if (base->type == EXPR_BIGINT && exp->type == EXPR_INTEGER) {
        int64_t e = exp->data.integer;
        if (e > 0) {
            return bigint_pow(base, e);
        }
        if (e < 0) {
            Expr* denom = bigint_pow(base, -e);
            if (!denom) return NULL;
            Expr* r_args[2] = { expr_new_integer(1), denom };
            return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
        }
    }

    if (base->type == EXPR_INTEGER && exp->type == EXPR_INTEGER) {
        int64_t b = base->data.integer;
        int64_t e = exp->data.integer;
        if (e > 0) {
            bool overflow = false;
            int64_t res_val = ipow(b, e, &overflow);
            if (overflow) return bigint_pow(base, e);
            return expr_new_integer(res_val);
        }
        if (e < 0) {
            bool overflow = false;
            int64_t res_val = ipow(b, -e, &overflow);
            if (overflow) {
                Expr* denom = bigint_pow(base, -e);
                if (!denom) return NULL;
                Expr* r_args[2] = { expr_new_integer(1), denom };
                return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
            }
            return make_rational(1, res_val);
        }
    }

    int64_t bn, bd, e;
    if (is_rational(base, &bn, &bd) && exp->type == EXPR_INTEGER) {
        e = exp->data.integer;
        if (e == 0) return expr_new_integer(1);
        if (e > 0) {
            bool overflow = false;
            int64_t res_n = ipow(bn, e, &overflow);
            int64_t res_d = ipow(bd, e, &overflow);
            if (!overflow) return make_rational(res_n, res_d);
        } else {
            bool overflow = false;
            int64_t res_n = ipow(bn, -e, &overflow);
            int64_t res_d = ipow(bd, -e, &overflow);
            if (!overflow) return make_rational(res_d, res_n);
        }
    }

    /* Bigint-aware Rational-base / integer-exponent path.  Triggers
     * when the int64 is_rational above fails because either the
     * Rational components or the exponent are BigInt -- without this,
     * Power[Rational[bignum, bignum], -1] (a common shape produced by
     * symbolic and bignum-driven linear-algebra routines) was leaked
     * unsimplified as Times[..., Power[Rational[..,..],-1]], breaking
     * downstream is_zero_poly / Together-based simplification. */
    if (base->type == EXPR_FUNCTION &&
        base->data.function.head->type == EXPR_SYMBOL &&
        base->data.function.head->data.symbol == SYM_Rational &&
        base->data.function.arg_count == 2 &&
        expr_is_integer_like(base->data.function.args[0]) &&
        expr_is_integer_like(base->data.function.args[1]) &&
        (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT)) {
        Expr* num_expr = base->data.function.args[0];
        Expr* den_expr = base->data.function.args[1];

        /* Exponent must fit in a long for mpz_pow_ui.  Out-of-range
         * |exponent| > 2^31 falls through to the symbolic path, which
         * is fine -- nobody asks for Power[Rational, 10^10] in
         * practice. */
        long e_abs = 0;
        bool e_neg = false;
        if (exp->type == EXPR_INTEGER) {
            int64_t ev = exp->data.integer;
            if (ev == 0) return expr_new_integer(1);
            e_neg = ev < 0;
            int64_t magnitude = e_neg ? -ev : ev;
            if (magnitude > 1000000000) goto rational_bigint_giveup;
            e_abs = (long)magnitude;
        } else { /* EXPR_BIGINT */
            if (mpz_sgn(exp->data.bigint) == 0) return expr_new_integer(1);
            if (!mpz_fits_slong_p(exp->data.bigint)) goto rational_bigint_giveup;
            long ev = mpz_get_si(exp->data.bigint);
            e_neg = ev < 0;
            e_abs = e_neg ? -ev : ev;
        }

        mpz_t bnum, bden, rnum, rden;
        mpz_init(bnum); mpz_init(bden);
        expr_to_mpz(num_expr, bnum);
        expr_to_mpz(den_expr, bden);
        mpz_init(rnum); mpz_init(rden);
        mpz_pow_ui(rnum, bnum, (unsigned long)e_abs);
        mpz_pow_ui(rden, bden, (unsigned long)e_abs);
        mpz_clear(bnum); mpz_clear(bden);

        /* Negative exponent flips numerator and denominator.  Restore
         * the canonical sign-on-numerator invariant make_rational
         * enforces:  denominator must be positive. */
        if (e_neg) {
            mpz_swap(rnum, rden);
        }
        if (mpz_sgn(rden) < 0) {
            mpz_neg(rnum, rnum);
            mpz_neg(rden, rden);
        }

        /* Cancel common factor. */
        mpz_t g; mpz_init(g);
        mpz_gcd(g, rnum, rden);
        if (mpz_cmp_ui(g, 1) != 0) {
            mpz_divexact(rnum, rnum, g);
            mpz_divexact(rden, rden, g);
        }
        mpz_clear(g);

        /* If denominator is 1, return the plain (possibly demoted) integer. */
        if (mpz_cmp_ui(rden, 1) == 0) {
            Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(rnum));
            mpz_clear(rnum); mpz_clear(rden);
            return out;
        }

        Expr* num_out = expr_bigint_normalize(expr_new_bigint_from_mpz(rnum));
        Expr* den_out = expr_bigint_normalize(expr_new_bigint_from_mpz(rden));
        mpz_clear(rnum); mpz_clear(rden);
        Expr* r_args[2] = { num_out, den_out };
        return expr_new_function(expr_new_symbol("Rational"), r_args, 2);

        rational_bigint_giveup:;
    }

    /* Rational-base, rational-exponent canonicalisation.
     *
     *   Power[Rational[n, d], p/q]   (q > 1, d > 1, gcd(p, q) = 1)
     *      ->   Power[n, p/q] * Power[d, -p/q]
     *
     * The split routes each integer piece through the existing
     * integer-base + rational-exponent path below, which extracts the
     * q-th-power coefficient and leaves a residue Power[r, b/q] with
     * r free of q-th-power factors. Triggering is gated on at least
     * one piece simplifying:
     *
     *   (a) |n| == 1   -- numerator collapses, and the result is a
     *       pure d^(-p/q) which routes cleanly through the integer
     *       path (find_min_perfect_base can further reduce d there);
     *   (b) m_n > 1    -- numerator has a perfect q-th-power factor;
     *   (c) m_d > 1    -- denominator has a perfect q-th-power factor.
     *
     * Without any of these triggers we leave Power[Rational[n,d], p/q]
     * alone, matching Mathematica's conservative behaviour for cases
     * like Sqrt[2/3] and (4/9)^(2/3) where no rational coefficient can
     * be extracted on either side.
     *
     * Soundness for negative numerator: (n/d)^(p/q) =
     * Power[n, p/q] * Power[d, -p/q] holds on the principal branch
     * because d is strictly positive (make_rational keeps d > 0). A
     * negative n routes through the integer-base negative-base
     * handling: q == 2 pulls out an I, odd q absorbs the sign, even
     * q >= 4 leaves Power[n, p/q] unevaluated (still progress -- the
     * Rational has been split into its integer components).
     *
     * Examples:
     *   (1/54)^(2/3)  -> 1 * 54^(-2/3) -> 1/9 * 2^(-2/3)
     *   (8/27)^(2/3)  -> 8^(2/3) * 27^(-2/3) -> 4 * 1/9 -> 4/9
     *   Sqrt[4/9]     -> Sqrt[4] * 9^(-1/2) -> 2 * 1/3 -> 2/3
     *   Sqrt[1/4]     -> 1 * 4^(-1/2) -> 1/2
     *   (-1/8)^(1/3)  -> (-1)^(1/3) * 8^(-1/3) -> (-1)^(1/3) / 2
     *   (1/4)^(2/3)   -> 4^(-2/3) -> 1/2 * 2^(-1/3) (via perfect-power)
     */
    {
        int64_t pp_rb, qq_rb;
        if (is_rational(base, &bn, &bd) && bd > 1 &&
            is_rational(exp, &pp_rb, &qq_rb) && qq_rb > 1) {
            int64_t abs_bn = (bn < 0) ? -bn : bn;
            int64_t m_n_rb, r_n_rb, m_d_rb, r_d_rb;
            factor_out_kth_power(abs_bn, qq_rb, &m_n_rb, &r_n_rb);
            factor_out_kth_power(bd,     qq_rb, &m_d_rb, &r_d_rb);
            (void)r_n_rb; (void)r_d_rb;
            bool trigger = (abs_bn == 1) || (m_n_rb > 1) || (m_d_rb > 1);
            if (trigger) {
                Expr* num_pow = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(bn), expr_copy(exp) }, 2);
                Expr* neg_exp = make_rational(-pp_rb, qq_rb);
                Expr* den_pow = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(bd), neg_exp }, 2);
                Expr* result = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ num_pow, den_pow }, 2);
                return eval_and_free(result);
            }
        }
    }

    if (exp->type == EXPR_INTEGER && base->type == EXPR_FUNCTION && base->data.function.head->data.symbol == SYM_Times) {
        size_t bc = base->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * bc);
        for (size_t i = 0; i < bc; i++) {
            Expr* p_args[2] = { expr_copy(base->data.function.args[i]), expr_copy(exp) };
            new_args[i] = expr_new_function(expr_new_symbol("Power"), p_args, 2);
        }
        Expr* result = expr_new_function(expr_new_symbol("Times"), new_args, bc);
        free(new_args);
        return result;
    }

    /* Rational-exponent Times-base distribution for the pure-imaginary
     * case: Power[Times[..., Complex[0, k_i], ..., positive_factors],
     * p/q] -> Power[-1, k_sum * p / (2q)] * Power[positive_part, p/q].
     *
     * Triggers only when the base contains at least one factor of the
     * form Complex[0, k] (a pure imaginary) and EVERY other factor is
     * known strictly positive (positive integer/rational, Pi, E). Without
     * the positivity guard we would silently cross a branch cut for
     * Sqrt[I x] (x of unknown sign).
     *
     * Examples:
     *   (I Pi)^(1/2)    -> (-1)^(1/4) Sqrt[Pi]
     *   (-I Pi)^(1/2)   -> (-1)^(-1/4) Sqrt[Pi]   == (-1)^(7/4) Sqrt[Pi]
     *   Sqrt[I 2 Pi]    -> (-1)^(1/4) Sqrt[2 Pi]
     *   Sqrt[I x]       -> unevaluated (x sign unknown)
     */
    if (base->type == EXPR_FUNCTION &&
        base->data.function.head->type == EXPR_SYMBOL &&
        base->data.function.head->data.symbol == SYM_Times) {
        int64_t pp_d, qq_d;
        if (is_rational(exp, &pp_d, &qq_d) && qq_d > 1) {
            size_t bc = base->data.function.arg_count;
            int64_t imag_count = 0;       /* signed count of I units */
            bool ok = true;
            bool saw_imag = false;
            Expr** pos_factors = malloc(sizeof(Expr*) * bc);
            size_t pos_count = 0;
            mpz_t coef_num, coef_den;
            mpz_init_set_ui(coef_num, 1);
            mpz_init_set_ui(coef_den, 1);
            for (size_t i = 0; i < bc; i++) {
                Expr* f = base->data.function.args[i];
                Expr *re = NULL, *im = NULL;
                if (is_complex(f, &re, &im)) {
                    /* Pure imaginary requires re == 0. */
                    bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0)
                                || (re->type == EXPR_REAL    && re->data.real == 0.0);
                    if (!re_zero) { ok = false; break; }
                    int64_t kn, kd;
                    if (im->type == EXPR_INTEGER) { kn = im->data.integer; kd = 1; }
                    else if (is_rational(im, &kn, &kd)) {}
                    else { ok = false; break; }
                    if (kn == 0) { ok = false; break; }
                    if (kn > 0) imag_count += 1;
                    else        { imag_count -= 1; kn = -kn; }
                    /* accumulate |k| into coefficient */
                    mpz_t tn, td;
                    mpz_init_set_si(tn, kn);
                    mpz_init_set_si(td, kd > 0 ? kd : -kd);
                    mpz_mul(coef_num, coef_num, tn);
                    mpz_mul(coef_den, coef_den, td);
                    mpz_clear(tn); mpz_clear(td);
                    saw_imag = true;
                } else {
                    bool pos = is_known_positive_pwr(f);
                    if (!pos) { ok = false; break; }
                    pos_factors[pos_count++] = expr_copy(f);
                }
            }
            if (ok && saw_imag) {
                /* coef_part = (coef_num / coef_den)^(p/q) */
                Expr* coef_part = NULL;
                if (mpz_cmp_ui(coef_num, 1) != 0 || mpz_cmp_ui(coef_den, 1) != 0) {
                    Expr* num_e = expr_bigint_normalize(expr_new_bigint_from_mpz(coef_num));
                    Expr* coef_base;
                    if (mpz_cmp_ui(coef_den, 1) == 0) {
                        coef_base = num_e;
                    } else {
                        Expr* den_e = expr_bigint_normalize(expr_new_bigint_from_mpz(coef_den));
                        coef_base = expr_new_function(expr_new_symbol("Rational"),
                                                     (Expr*[]){num_e, den_e}, 2);
                    }
                    coef_part = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ coef_base, expr_copy(exp) }, 2);
                }
                /* pos_part = Power[Times[pos_factors...], p/q]. The
                 * recursive Power evaluator will further canonicalise
                 * (Pi^(1/2) -> Sqrt[Pi], etc.). */
                Expr* pos_part = NULL;
                if (pos_count == 1) {
                    pos_part = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ pos_factors[0], expr_copy(exp) }, 2);
                } else if (pos_count > 1) {
                    Expr* pos_times = expr_new_function(expr_new_symbol("Times"),
                                                        pos_factors, pos_count);
                    pos_part = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ pos_times, expr_copy(exp) }, 2);
                } else {
                    free(pos_factors);
                    pos_factors = NULL;
                }
                /* sign_part = Power[-1, imag_count * p / (2 * q)]. The
                 * exponent is built as Rational[imag_count*pp_d, 2*qq_d]
                 * and reduced/normalised by make_rational. */
                Expr* sign_exp;
                {
                    /* careful with int64 overflow: imag_count * pp_d. */
                    __int128_t top = (__int128_t)imag_count * (__int128_t)pp_d;
                    __int128_t bot = (__int128_t)2 * (__int128_t)qq_d;
                    if (top >= INT64_MIN && top <= INT64_MAX &&
                        bot >= INT64_MIN && bot <= INT64_MAX) {
                        sign_exp = make_rational((int64_t)top, (int64_t)bot);
                    } else {
                        /* overflow guard: bail */
                        if (pos_part) expr_free(pos_part);
                        if (coef_part) expr_free(coef_part);
                        if (pos_factors) {
                            for (size_t i = 0; i < pos_count; i++) expr_free(pos_factors[i]);
                            free(pos_factors);
                        }
                        mpz_clear(coef_num); mpz_clear(coef_den);
                        goto rat_imag_fallthrough;
                    }
                }
                Expr* sign_part = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(-1), sign_exp }, 2);
                /* assemble Times[sign_part, coef_part?, pos_part?] */
                Expr* parts[3];
                size_t pcount = 0;
                parts[pcount++] = sign_part;
                if (coef_part) parts[pcount++] = coef_part;
                if (pos_part)  parts[pcount++] = pos_part;
                mpz_clear(coef_num); mpz_clear(coef_den);
                /* free pos_factors array (Expr* contents now owned by
                 * pos_times -> pos_part). */
                if (pos_factors) free(pos_factors);
                if (pcount == 1) return parts[0];
                Expr** result_args = malloc(sizeof(Expr*) * pcount);
                for (size_t i = 0; i < pcount; i++) result_args[i] = parts[i];
                Expr* result = expr_new_function(expr_new_symbol("Times"),
                                                 result_args, pcount);
                free(result_args);
                return eval_and_free(result);
            }
            /* clean up the pos_factors copies if we bailed */
            for (size_t i = 0; i < pos_count; i++) expr_free(pos_factors[i]);
            free(pos_factors);
            mpz_clear(coef_num); mpz_clear(coef_den);
        }
    }
rat_imag_fallthrough: ;

    /* Power[Times[positive_factors], p/q] distribution.
     *
     * Distributes the outer Power over every Times factor when:
     *   (a) every factor is known positive (is_known_positive_pwr), and
     *   (b) at least one factor either:
     *       (i) is an integer or rational that PERFECTLY reduces to a
     *           rational coefficient under Power[_, p/q] (residue 1
     *           on both numerator and denominator), OR
     *       (ii) is a Power[positive, p'/q'] whose composition with
     *            the outer p/q introduces no new radical complexity
     *            (composed reduced denominator <= q'), per
     *            power_factor_composes_cleanly.
     *
     * The gating in (b) keeps Mathematica-faithful behaviour:
     *
     *   Sqrt[4 Pi]              -> 2 Sqrt[Pi]      (4 = 2^2, r=1)
     *   Sqrt[(1/9) 2^(-2/3)]    -> 1/3 * 2^(-1/3)  (9 = 3^2, r=1)
     *   (Pi/8)^(2/3)            -> Pi^(2/3) / 4    (8 = 2^3, r=1)
     *
     * vs. leaving (4 Pi)^(2/3) and Sqrt[2 Pi] alone (no factor has a
     * clean rational extraction). Pure-imaginary factors are handled
     * by the earlier rat_imag block, which fires before this one and
     * is gated on saw_imag; falling through here means no imaginary
     * factor is present.
     *
     * After distribution, each resulting Power[factor, p/q] is
     * re-evaluated. Power-typed factors hit the Power-of-Power
     * composition rule below (positive base lets exponents merge
     * without crossing a branch cut). Rational/integer factors hit
     * their dedicated paths. This is the canonical-form companion of
     * the rational-base distribution above: that one factors a single
     * Rational base; this one factors a Times of positives that
     * arises naturally as the output of a prior simplification (e.g.
     * Sqrt[(1/54)^(2/3)] after (1/54)^(2/3) reduces to
     * 1/9 * 2^(-2/3)).
     */
    if (base->type == EXPR_FUNCTION && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Times) {
        int64_t pp_pt, qq_pt;
        if (is_rational(exp, &pp_pt, &qq_pt) && qq_pt > 1) {
            size_t bc = base->data.function.arg_count;
            bool all_pos = true;
            bool has_full_reducer = false;
            for (size_t i = 0; i < bc; i++) {
                Expr* f = base->data.function.args[i];
                if (!is_known_positive_pwr(f)) { all_pos = false; break; }
                if (factor_fully_reduces_under_q(f, qq_pt)
                    || power_factor_composes_cleanly(f, pp_pt, qq_pt)) {
                    has_full_reducer = true;
                }
            }
            if (all_pos && has_full_reducer) {
                Expr** new_args = malloc(sizeof(Expr*) * bc);
                for (size_t i = 0; i < bc; i++) {
                    new_args[i] = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(base->data.function.args[i]),
                                   expr_copy(exp) }, 2);
                }
                Expr* result = expr_new_function(expr_new_symbol("Times"),
                                                 new_args, bc);
                free(new_args);
                return eval_and_free(result);
            }
        }
    }

    /* Power-of-Power composition: Power[Power[B, E_inner], exp] ->
     * Power[B, E_inner * exp].
     *
     * Sound when (a) exp is an integer (the identity holds for any B
     * since integer powers don't cross branch cuts), or (b) B is known
     * positive (the principal-branch identity (B^E_inner)^exp =
     * B^(E_inner * exp) holds without crossing branch cuts). Without
     * the positivity guard for non-integer exp we would silently break
     * Sqrt[(-1)^2] = 1 vs (-1)^(2 * 1/2) = -1.
     *
     * The positive-base extension is what reduces Sqrt[Power[2, -2/3]]
     * to Power[2, -1/3] -- needed for the inner radical that emerges
     * from Sqrt[(1/54)^(2/3)] after the rational-base distribution
     * splits it into 1/9 * Power[2, -2/3].
     */
    if (base->type == EXPR_FUNCTION && base->data.function.head->type == EXPR_SYMBOL
        && base->data.function.head->data.symbol == SYM_Power
        && base->data.function.arg_count == 2) {
        Expr* inner_base = base->data.function.args[0];
        Expr* inner_exp  = base->data.function.args[1];
        bool can_compose = (exp->type == EXPR_INTEGER)
                        || is_known_positive_pwr(inner_base);
        if (can_compose) {
            Expr* t_args[2] = { expr_copy(inner_exp), expr_copy(exp) };
            Expr* new_exp = expr_new_function(expr_new_symbol("Times"), t_args, 2);
            Expr* p_args[2] = { expr_copy(inner_base), new_exp };
            return expr_new_function(expr_new_symbol("Power"), p_args, 2);
        }
    }

    {
        Expr* simp = simplify_exp_log(base, exp);
        if (simp) return simp;
    }

    /* Perfect-power base unification: Power[b^k, p/q] -> Power[b, k*p/q]
     * with k maximal (b smallest). Run before integer-part extraction so
     * the residue under the radical is always over the smallest possible
     * base. Examples:
     *   4^(2/3)   -> 2^(4/3) -> 2 * 2^(1/3)
     *   9^(1/3)   -> 3^(2/3)
     *   8^(5/3)   -> 2^5 = 32
     *   1024^(1/2) -> 2^5 = 32
     * Restricted to positive bases (negative perfect powers don't
     * commute with arbitrary p/q on the principal branch -- see
     * find_min_perfect_base for the full argument). */
    {
        int64_t pp_, qq_;
        /* Accept both EXPR_INTEGER (>= 4, the threshold below which
         * find_min_perfect_base cannot make progress: 0,1 are degenerate
         * and 2,3 are primes) and EXPR_BIGINT (always >= 4 here since any
         * BigInt that survived normalisation does not fit in int64). The
         * downstream find_min_perfect_base wants the full mpz_t, so the
         * conversion goes through expr_to_mpz uniformly. */
        bool base_is_pp_candidate =
            (base->type == EXPR_INTEGER && base->data.integer >= 4)
            || (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) > 0);
        if (base_is_pp_candidate && is_rational(exp, &pp_, &qq_) && qq_ > 1) {
            mpz_t n_z, b_z;
            expr_to_mpz(base, n_z);
            mpz_init(b_z);
            int64_t k_total = 0;
            if (find_min_perfect_base(n_z, b_z, &k_total)) {
                /* New exponent = k_total * pp_ / qq_; check int64 overflow. */
                __int128_t kp = (__int128_t)k_total * (__int128_t)pp_;
                if (kp <= INT64_MAX && kp >= INT64_MIN) {
                    int64_t new_p = (int64_t)kp;
                    Expr* new_base_e;
                    if (mpz_fits_slong_p(b_z)) {
                        new_base_e = expr_new_integer((int64_t)mpz_get_si(b_z));
                    } else {
                        new_base_e = expr_bigint_normalize(expr_new_bigint_from_mpz(b_z));
                    }
                    /* make_rational reduces and collapses to integer when
                     * the gcd makes the denominator 1. */
                    Expr* new_exp_e = make_rational(new_p, qq_);
                    Expr* new_pow = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ new_base_e, new_exp_e }, 2);
                    mpz_clear(b_z);
                    mpz_clear(n_z);
                    /* Recurse via evaluate so the reduced expression flows
                     * through the integer-part extraction below. The new
                     * base is strictly smaller, so termination is guaranteed
                     * (it will not be a perfect power again). evaluate()
                     * copies its input; the caller must free new_pow.
                     * Per the BuiltinFunc convention, we do NOT free res. */
                    Expr* result = evaluate(new_pow);
                    expr_free(new_pow);
                    return result;
                }
            }
            mpz_clear(b_z);
            mpz_clear(n_z);
        }
    }

    /* Integer base, rational exponent p/q with q > 1.
     *
     * For positive base n, decompose n^(p/q) into an integer/rational
     * coefficient and a residual radical:
     *
     *   1. Integer-part extraction (only for p >= q; p < q including all
     *      negative p stays in the residual). Write p = a*q + b with
     *      0 <= b < q so n^(p/q) = n^a * n^(b/q).
     *   2. Perfect-q-th-power reduction. Write n = m^q * r with r free of
     *      q-th-power factors, so n^(b/q) = m^b * r^(b/q).
     *
     * Combined: n^(p/q) = (n^a * m^b) * r^(b/q).
     *
     * Negative b (when 0 < -p < q, i.e., we kept a=0 and b=p<0) yields a
     * rational coefficient n^a / m^|b| and a Power[r, b/q] residue with a
     * negative exponent — the surrounding Times canonicalizer fuses
     * those into a single Power if a same-base coefficient is present.
     *
     * We compute the coefficient in GMP to handle bigint promotion
     * automatically; the residue is at most a single Power node. No
     * recursive evaluator round-trips. */
    int64_t p, q;
    /* Negative BigInt base with even q == 2: route through the existing
     * Sqrt[-n] = I^p * Sqrt[|n|] identity by recursing on the positive
     * BigInt. The int64-only path below cannot fit |n| when it overflows;
     * this specialised branch is the BigInt equivalent of the q == 2
     * arm at the head of the EXPR_INTEGER block. */
    if (base->type == EXPR_BIGINT
        && mpz_sgn(base->data.bigint) < 0
        && is_rational(exp, &p, &q) && q == 2) {
        mpz_t pos;
        mpz_init(pos);
        mpz_neg(pos, base->data.bigint);
        Expr* pos_base = expr_bigint_normalize(expr_new_bigint_from_mpz(pos));
        mpz_clear(pos);

        int64_t pmod = ((p % 4) + 4) % 4;
        int i_sign = (pmod == 1) ? 1 : -1;
        Expr* i_val = expr_new_function(expr_new_symbol("Complex"),
            (Expr*[]){expr_new_integer(0), expr_new_integer(i_sign)}, 2);
        Expr* pos_power_args[2] = { pos_base, expr_copy(exp) };
        Expr* pos_power = expr_new_function(expr_new_symbol("Power"), pos_power_args, 2);
        Expr* rest = evaluate(pos_power);
        expr_free(pos_power);

        bool rest_is_one = (rest->type == EXPR_INTEGER && rest->data.integer == 1);
        if (rest_is_one) {
            expr_free(rest);
            return i_val;
        }
        Expr* t_args[2] = { i_val, rest };
        return expr_new_function(expr_new_symbol("Times"), t_args, 2);
    }

    if (base->type == EXPR_INTEGER && is_rational(exp, &p, &q) && q > 1) {
        int64_t n = base->data.integer;
        bool n_negative = (n < 0);
        if (n_negative) {
            /* q==2: (-n)^(p/2) = I^p * |n|^(p/2) on the principal branch.
             * p is odd in canonical form (gcd(p, 2) = 1), so I^p reduces to
             *   p mod 4 == 1  ->  I
             *   p mod 4 == 3  ->  -I
             * (Negative p is normalised by ((p % 4) + 4) % 4.) Covers both
             * Sqrt[-n] = I Sqrt[n] (p == 1) and the higher-power cases like
             * (-1)^(3/2) = -I, (-12)^(3/2) = -24 I Sqrt[3]. */
            if (q == 2) {
                int64_t pmod = ((p % 4) + 4) % 4;
                int i_sign = (pmod == 1) ? 1 : -1;     /* p odd => pmod in {1,3} */
                Expr* i_val = expr_new_function(expr_new_symbol("Complex"),
                    (Expr*[]){expr_new_integer(0), expr_new_integer(i_sign)}, 2);
                Expr* pos_base = expr_new_integer(-n);
                Expr* tmp_p_args[2] = { pos_base, expr_copy(exp) };
                Expr* tmp_power = expr_new_function(expr_new_symbol("Power"), tmp_p_args, 2);
                Expr* rest = builtin_power(tmp_power);
                if (!rest) rest = tmp_power; else expr_free(tmp_power);
                bool rest_is_one = (rest->type == EXPR_INTEGER
                                    && rest->data.integer == 1);
                if (rest_is_one) {
                    expr_free(rest);
                    return i_val;
                }
                Expr* t_args[2] = { i_val, rest };
                return expr_new_function(expr_new_symbol("Times"), t_args, 2);
            }
            /* Even q >= 4: for base == -1 we can still do integer-part
             * extraction -- (-1)^(p/q) = (-1)^a_int * (-1)^(b_rem/q) with
             * |b_rem| < q -- so e.g. (-1)^(5/4) -> -(-1)^(1/4). For other
             * negative bases the principal-branch form has no clean
             * canonical representation here yet; leave it unevaluated. */
            if (q % 2 == 0 && n != -1) return NULL;
            /* Odd q (any negative base), or base == -1 with even q: fall
             * through.  The transformation is
             *   (-n)^(p/q) = (-n)^a * (-n)^(b/q)        [where p = a*q + b]
             *              = (-1)^a * n^a * (-1)^(b/q) * n^(b/q)
             *              = (-1)^a * (m^b * n^a) * (-r)^(b/q)
             * which canonicalizes to a (signed) integer/rational coefficient
             * times Power[-r, b/q]. The sign comes from (-1)^a. */
        }
        int64_t abs_n = n_negative ? -n : n;

        int64_t a_int = 0;
        int64_t b_rem = p;
        /* Integer-part extraction: write p = a_int*q + b_rem with |b_rem|<q
         * and b_rem same-sign as p. C99 `/` truncates toward zero, which
         * gives exactly the canonical Mathematica decomposition for both
         * signs (e.g. -5/3 -> a_int=-1, b_rem=-2, so 2^(-5/3) = 2^-1 *
         * 2^(-2/3)). Skip when |p|<q (no integer part to extract; the
         * perfect-q-th-power reduction below still runs). */
        if (p >= q || p <= -q) {
            a_int = p / q;
            b_rem = p - a_int * q;
        }

        int64_t m, r;
        factor_out_kth_power(abs_n, q, &m, &r);

        if (a_int != 0 || m > 1) {
            /* Build coeff = abs_n^a_int * m^b_rem in GMP, distributing each
             * factor to either numerator or denominator by its sign. Both
             * a_int and b_rem may be negative; the GCD reduction below
             * normalises into a canonical Rational[num, den] (or plain
             * integer when den==1). */
            mpz_t num_z, den_z;
            mpz_init_set_ui(num_z, 1);
            mpz_init_set_ui(den_z, 1);
            if (a_int != 0) {
                mpz_t t; mpz_init_set_si(t, abs_n);
                mpz_pow_ui(t, t, (unsigned long)(a_int > 0 ? a_int : -a_int));
                if (a_int > 0) mpz_mul(num_z, num_z, t);
                else            mpz_mul(den_z, den_z, t);
                mpz_clear(t);
            }
            if (b_rem != 0 && m > 1) {
                mpz_t t; mpz_init_set_si(t, m);
                mpz_pow_ui(t, t, (unsigned long)(b_rem > 0 ? b_rem : -b_rem));
                if (b_rem > 0) mpz_mul(num_z, num_z, t);
                else            mpz_mul(den_z, den_z, t);
                mpz_clear(t);
            }
            mpz_t g_z; mpz_init(g_z);
            mpz_gcd(g_z, num_z, den_z);
            if (mpz_cmp_ui(g_z, 1) > 0) {
                mpz_divexact(num_z, num_z, g_z);
                mpz_divexact(den_z, den_z, g_z);
            }
            mpz_clear(g_z);

            Expr* coeff;
            if (mpz_cmp_ui(den_z, 1) == 0) {
                coeff = expr_bigint_normalize(expr_new_bigint_from_mpz(num_z));
            } else {
                Expr* num_e = expr_bigint_normalize(expr_new_bigint_from_mpz(num_z));
                Expr* den_e = expr_bigint_normalize(expr_new_bigint_from_mpz(den_z));
                coeff = expr_new_function(expr_new_symbol("Rational"),
                                          (Expr*[]){num_e, den_e}, 2);
            }
            mpz_clear(num_z); mpz_clear(den_z);

            /* For negative n, integer-part a contributes (-1)^a to the
             * coefficient. Even a -> +; odd a -> negate. */
            if (n_negative && (a_int % 2 != 0)) {
                coeff = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), coeff }, 2);
                coeff = eval_and_free(coeff);
            }

            /* residue = Power[res_base, b_rem/q] where res_base is
             * +r for positive n and -r for negative n (the (-1)^(b/q)
             * factor merges with r^(b/q) into (-r)^(b/q)). When b_rem == 0
             * (only reachable if q == 1, excluded here) or res_base == 1
             * the residue is 1. For negative n with r == 1, res_base == -1
             * and we keep Power[-1, b/q] as the residue (Mathematica form
             * "(-1)^(b/q)"). */
            int64_t res_base = n_negative ? -r : r;
            Expr* residue;
            if (b_rem == 0 || res_base == 1) {
                residue = expr_new_integer(1);
            } else {
                /* p/q is in lowest terms by construction, so gcd(b_rem, q)
                 * divides gcd(p, q) = 1. make_rational still normalises. */
                Expr* new_exp = make_rational(b_rem, q);
                residue = expr_new_function(expr_new_symbol("Power"),
                                            (Expr*[]){expr_new_integer(res_base), new_exp}, 2);
            }

            bool coeff_is_one = (coeff->type == EXPR_INTEGER && coeff->data.integer == 1);
            bool residue_is_one = (residue->type == EXPR_INTEGER && residue->data.integer == 1);
            if (residue_is_one) { expr_free(residue); return coeff; }
            if (coeff_is_one) { expr_free(coeff); return residue; }
            return expr_new_function(expr_new_symbol("Times"),
                                     (Expr*[]){coeff, residue}, 2);
        }

        /* Heterogeneous-prime residue split. Reached when a_int == 0
         * AND m == 1, i.e. abs_n has no perfect q-th-power factor and
         * the exponent's |numerator| < q so no integer part to extract
         * -- the existing block above bails. If abs_n still factors
         * into multiple distinct primes with non-uniform per-prime
         * effective exponents (a_i * b_rem mod q), split into a
         * canonical product of distinct-prime powers. This is what
         * upgrades 18^(1/3) from unevaluated to 2^(1/3) * 3^(2/3),
         * 12^(1/3) to 2^(2/3) * 3^(1/3), etc. For uniform-exponent
         * residues like 6^(1/3) and 30^(1/3), power_split_residue
         * returns NULL and the input form survives.
         *
         * Restricted to positive base and positive b_rem; mixed-sign
         * cases keep the existing residue shape. The b_rem > 0
         * variants of integer-extracted residues (e.g. inside the
         * 18^(5/3) = 18 * 18^(2/3) construction above) are picked up
         * on the recursive evaluator pass: Power[18, 2/3] re-enters
         * builtin_power and hits this same path. */
        if (!n_negative && b_rem > 0 && r > 1) {
            Expr* split = power_split_residue(r, b_rem, q);
            if (split) return split;
        }
    }

    /* Trig / hyperbolic reciprocal naming: Power[Cos[x], -k] -> Sec[x]^k, etc.
     * The parser collapses 1/Cos[x] to Power[Cos[x], -1] without going through
     * Times, so the Times-level canonicalizer never sees these solo cases. */
    if (exp->type == EXPR_INTEGER) {
        Expr* tc = trig_canon_power(base, exp->data.integer);
        if (tc) return tc;
    }

    return NULL;
}

Expr* builtin_sqrt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* half = make_rational(1, 2);
    Expr* p_args[2] = { expr_copy(arg), half };
    return expr_new_function(expr_new_symbol("Power"), p_args, 2);
}
