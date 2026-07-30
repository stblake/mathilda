/*
 * random.c - Random number generation builtins for Mathilda
 *
 * RandomInteger[{imin, imax}]  - pseudorandom integer in [imin, imax]
 * RandomInteger[imax]          - pseudorandom integer in [0, imax]
 * RandomInteger[]              - pseudorandomly gives 0 or 1
 * RandomInteger[range, n]      - list of n pseudorandom integers
 * RandomInteger[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomReal[]                 - pseudorandom real in [0, 1)
 * RandomReal[xmax]             - pseudorandom real in [0, xmax)
 * RandomReal[{xmin, xmax}]    - pseudorandom real in [xmin, xmax)
 * RandomReal[range, n]         - list of n pseudorandom reals
 * RandomReal[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomComplex[]              - pseudorandom complex in rectangle [0,1]x[0,1]
 * RandomComplex[zmax]          - pseudorandom complex in rectangle [0,zmax]
 * RandomComplex[{zmin, zmax}]  - pseudorandom complex in rectangle [zmin,zmax]
 * RandomComplex[range, n]      - list of n pseudorandom complex numbers
 * RandomComplex[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomChoice[list]           - pseudorandom choice from list
 * RandomChoice[list, n]        - list of n pseudorandom choices
 * RandomChoice[list, {n1, n2, ...}] - n1*n2*... array of choices
 * RandomChoice[wts->elems]     - weighted pseudorandom choice
 * RandomChoice[wts->elems, n]  - list of n weighted choices
 * RandomChoice[wts->elems, {n1,...}] - array of weighted choices
 *
 * SeedRandom[n]                - seeds the RNG with integer n
 * SeedRandom[]                 - reseeds from system entropy
 *
 * Uses GMP's random state (gmp_randstate_t) for bignum support.
 */

#include "random.h"
#include "checked_int.h"   /* ci_add_i64 */
#include "arithmetic.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include "numeric.h"
#include "ndarray.h"
#include "pack.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Global GMP random state. Retained for everything that genuinely needs
 * arbitrary precision: bignum integer ranges and the MPFR draws. */
static gmp_randstate_t g_rand_state;
static int g_rand_initialized = 0;

/* --------------------------------------------------------------------------
 * Machine-precision generator: xoshiro256++
 * --------------------------------------------------------------------------
 *
 * WHY A SECOND GENERATOR. Every machine-precision draw used to go through GMP:
 * random_uniform_01 did mpz_init x2, recomputed 2^53 with mpz_ui_pow_ui, called
 * mpz_urandomm, then mpz_clear x2 -- two heap allocations and a BIGNUM MODULAR
 * draw to produce 53 bits. Measured at 117 ns per double, which was essentially
 * the entire cost of RandomReal: RandomReal[{0,1}, 10^7] took 1.22 s, and 10^7 x
 * 117 ns is 1.17 s of it.
 *
 * Measured alternatives, same machine, 10^7 draws:
 *     mpz per call (what this replaced)   117.2 ns
 *     hoisted mpz scratch + modulus        18.3 ns   (6.4x)
 *     gmp_urandomb_ui, no mpz               9.0 ns   (13x)
 *     xoshiro256++                          2.4 ns   (49x)
 *
 * xoshiro256++ is the standard choice for exactly this job: 256 bits of state,
 * period 2^256-1, passes BigCrush, and one call is a handful of shifts and xors.
 * It is NOT cryptographic, which RandomReal never claimed to be.
 *
 * GMP KEEPS the bignum and MPFR paths, so arbitrary-precision draws are
 * unchanged. Both generators are seeded together and saved/restored together
 * (see seed_all and random_push_seed), so SeedRandom stays a single coherent
 * operation over both.
 * -------------------------------------------------------------------------- */

static uint64_t g_xs[4];

/* splitmix64: the reference seeder for xoshiro. One 64-bit seed expands to the
 * 256-bit state with good avalanche, so adjacent seeds give unrelated streams --
 * seeding the state directly from a small integer would leave most of it zero. */
static uint64_t splitmix64(uint64_t* x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void xs_seed(uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; i++) g_xs[i] = splitmix64(&sm);
    /* The all-zero state is xoshiro's one fixed point; splitmix64 cannot produce
     * it from any seed, but a corrupted restore could. */
    if (!(g_xs[0] | g_xs[1] | g_xs[2] | g_xs[3])) g_xs[0] = 0x9E3779B97F4A7C15ULL;
}

static inline uint64_t xs_rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t xs_next(void) {
    uint64_t r = xs_rotl(g_xs[0] + g_xs[3], 23) + g_xs[0];
    uint64_t t = g_xs[1] << 17;
    g_xs[2] ^= g_xs[0];
    g_xs[3] ^= g_xs[1];
    g_xs[1] ^= g_xs[2];
    g_xs[0] ^= g_xs[3];
    g_xs[2] ^= t;
    g_xs[3] = xs_rotl(g_xs[3], 45);
    return r;
}

/*
 * Uniform in [0, m), with no modulo bias and -- in the overwhelmingly common
 * case -- no division either.
 *
 * The obvious rejection form computes `threshold = 2^64 mod m` and then
 * `r % m`: TWO 64-bit divisions per draw, each ~20-40 cycles. That made
 * RandomInteger ~3x slower per element than RandomReal, which is only a shift
 * and a multiply. Two ways out, both taken:
 *
 *   - A power-of-two range is a mask. No division, and no rejection either,
 *     since 2^64 is an exact multiple. RandomInteger[1] and RandomInteger[{0,
 *     255}] and friends land here.
 *   - Otherwise Lemire's nearly-divisionless method: one 128-bit multiply, and
 *     the high half is already uniform on [0, m). The division survives only
 *     inside a branch entered with probability m/2^64 -- for any range a person
 *     writes, never.
 *
 * __uint128_t is a compiler extension, so it is guarded on __SIZEOF_INT128__
 * with the plain double-division form as the fallback. Both are exactly
 * uniform; they differ only in speed and in which value a given draw yields, so
 * a build without 128-bit integers has a different stream. That is acceptable
 * for the same reason SeedRandom does not promise a stream across versions.
 */
static inline uint64_t xs_below(uint64_t m) {
    if (m <= 1) return 0;
    if ((m & (m - 1)) == 0) return xs_next() & (m - 1);    /* power of two */
#if defined(__SIZEOF_INT128__)
    __uint128_t prod = (__uint128_t)xs_next() * (__uint128_t)m;
    uint64_t low = (uint64_t)prod;
    if (low < m) {                                  /* p = m/2^64: essentially never */
        uint64_t threshold = (uint64_t)(-(int64_t)m) % m;   /* == 2^64 mod m */
        while (low < threshold) {
            prod = (__uint128_t)xs_next() * (__uint128_t)m;
            low = (uint64_t)prod;
        }
    }
    return (uint64_t)(prod >> 64);
#else
    uint64_t threshold = (uint64_t)(-(int64_t)m) % m;
    uint64_t r;
    do { r = xs_next(); } while (r < threshold);
    return r % m;
#endif
}

/* Seed both generators from one value, so SeedRandom is coherent across the
 * machine and arbitrary-precision paths. */
static void seed_all(uint64_t seed) {
    gmp_randseed_ui(g_rand_state, (unsigned long)seed);
    xs_seed(seed);
}

/* Ensure RNG is initialized */
static void ensure_rand_init(void) {
    if (!g_rand_initialized) {
        gmp_randinit_mt(g_rand_state);
        g_rand_initialized = 1;
        seed_all((uint64_t)time(NULL) ^ (uint64_t)clock());
    }
}

/* Saved-state stack for random_push_seed/random_pop_seed. Depth 4 is ample:
 * the only caller (PossibleZeroQ's sampler) does not recurse into the RNG. */
#define RANDOM_SEED_STACK_DEPTH 4
static gmp_randstate_t g_rand_saved[RANDOM_SEED_STACK_DEPTH];
/* The xoshiro state is saved alongside it: a push/pop that restored only GMP
 * would leave the machine-precision stream advanced by whatever ran in between,
 * so PossibleZeroQ's sampler would perturb the user's RandomReal sequence --
 * exactly what random.h promises it does not do. */
static uint64_t g_xs_saved[RANDOM_SEED_STACK_DEPTH][4];
static int g_rand_saved_depth = 0;

void random_push_seed(uint64_t seed) {
    ensure_rand_init();
    if (g_rand_saved_depth >= RANDOM_SEED_STACK_DEPTH) {
        /* Stack full (unexpected): reseed without saving so behaviour stays
         * deterministic; the matching pop is a no-op (see guard below). */
        seed_all(seed);
        g_rand_saved_depth++;
        return;
    }
    gmp_randinit_set(g_rand_saved[g_rand_saved_depth], g_rand_state);
    for (int i = 0; i < 4; i++) g_xs_saved[g_rand_saved_depth][i] = g_xs[i];
    g_rand_saved_depth++;
    seed_all(seed);
}

void random_pop_seed(void) {
    if (g_rand_saved_depth <= 0) return;            /* unbalanced pop */
    g_rand_saved_depth--;
    if (g_rand_saved_depth >= RANDOM_SEED_STACK_DEPTH) return;  /* overflow push */
    gmp_randclear(g_rand_state);
    gmp_randinit_set(g_rand_state, g_rand_saved[g_rand_saved_depth]);
    gmp_randclear(g_rand_saved[g_rand_saved_depth]);
    for (int i = 0; i < 4; i++) g_xs[i] = g_xs_saved[g_rand_saved_depth][i];
}

/*
 * Internal sampling draw: uniform in [lo, hi], from the GMP stream.
 *
 * Deliberately NOT the xoshiro path that RandomInteger now uses. The zero-test
 * sampler (src/zero_test.c) draws its Schwartz-Zippel points through this, and
 * those points are part of a DECISION PROCEDURE, not a user-visible random
 * sequence: PossibleZeroQ seeds them from the expression's structural hash so a
 * verdict is a pure function of the input.
 *
 * Sharing the user-facing generator meant that speeding up RandomInteger moved
 * those sample points -- and the integrator's answers with them.
 * `Integrate[E^(Log[x]^2), x]` is non-elementary and must stay unintegrated;
 * with a different point set the zero test reached a different verdict and the
 * integrator claimed to solve it. That is a real pre-existing fragility (the
 * verdict should not be that sensitive), but it is not something an RNG
 * optimisation gets to change, so the sampler keeps the generator it was tuned
 * against and the coupling is removed.
 */
int64_t random_internal_int_range(int64_t lo, int64_t hi) {
    ensure_rand_init();
    if (hi < lo) return lo;
    static mpz_t span, draw;
    static int ready = 0;
    if (!ready) { mpz_init(span); mpz_init(draw); ready = 1; }
    mpz_set_si(span, (long)hi);
    mpz_sub_ui(span, span, (unsigned long)0);
    {   /* span = hi - lo + 1, computed in mpz so hi-lo cannot overflow */
        mpz_t lo_z; mpz_init_set_si(lo_z, (long)lo);
        mpz_sub(span, span, lo_z);
        mpz_add_ui(span, span, 1);
        mpz_urandomm(draw, g_rand_state, span);
        mpz_add(draw, draw, lo_z);
        mpz_clear(lo_z);
    }
    return (int64_t)mpz_get_si(draw);
}

/*
 * Generate a single random integer in [imin, imax] (inclusive).
 * imin and imax are given as mpz_t values.
 * Returns a new Expr* (EXPR_INTEGER or EXPR_BIGINT, normalized).
 */
static Expr* random_integer_range(const mpz_t imin, const mpz_t imax) {
    ensure_rand_init();

    /* Scratch reused across calls rather than mpz_init/mpz_clear per draw: this
     * is called once per element of RandomInteger[range, n], and two heap
     * allocations per element dominated the actual work. */
    static mpz_t range_size, result;
    static int scratch_ready = 0;
    if (!scratch_ready) { mpz_init(range_size); mpz_init(result); scratch_ready = 1; }

    mpz_sub(range_size, imax, imin);
    mpz_add_ui(range_size, range_size, 1);

    if (mpz_sgn(range_size) <= 0) return NULL;   /* empty or invalid range */

    /* Machine-width range AND a machine-width base: draw from xoshiro and build
     * the Expr straight from an int64. No mpz is touched at all.
     *
     * Building through expr_bigint_normalize(expr_new_bigint_from_mpz(...)) is
     * what kept RandomInteger slow after the generator was replaced -- an
     * mpz-backed Expr allocated and then normalised back down to a machine
     * integer, once per element. The draw is ~3 ns; that round trip was ~170. */
    if (mpz_fits_ulong_p(range_size) && mpz_fits_slong_p(imin)) {
        uint64_t span = (uint64_t)mpz_get_ui(range_size);
        int64_t  base = (int64_t)mpz_get_si(imin);
        int64_t  draw = (int64_t)xs_below(span);
        int64_t  out;
        /* imin + draw can still leave int64 even when both fit; fall through to
         * the exact path rather than wrap. */
        if (!ci_add_i64(base, draw, &out)) return expr_new_integer(out);
    }

    /* Wide range (or a base past int64): GMP draws it uniformly. */
    if (mpz_fits_ulong_p(range_size)) {
        mpz_set_ui(result, (unsigned long)xs_below((uint64_t)mpz_get_ui(range_size)));
        mpz_add(result, result, imin);
        return expr_bigint_normalize(expr_new_bigint_from_mpz(result));
    }

    /* Genuine bignum range: GMP is the only thing that can draw it uniformly. */
    mpz_urandomm(result, g_rand_state, range_size);
    mpz_add(result, result, imin);
    return expr_bigint_normalize(expr_new_bigint_from_mpz(result));
}

/*
 * Parse a range argument into imin, imax (as mpz_t).
 * Supports:
 *   - Integer or BigInt atom: range is [0, val]
 *   - List[imin, imax]: range is [imin, imax]
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_range(Expr* arg, mpz_t imin, mpz_t imax) {
    if (arg->type == EXPR_INTEGER) {
        mpz_set_si(imin, 0);
        mpz_set_si(imax, arg->data.integer);
        return true;
    }
    if (arg->type == EXPR_BIGINT) {
        mpz_set_si(imin, 0);
        mpz_set(imax, arg->data.bigint);
        return true;
    }
    /* Check for List[imin, imax] */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List &&
        arg->data.function.arg_count == 2) {

        Expr* lo = arg->data.function.args[0];
        Expr* hi = arg->data.function.args[1];

        if (!expr_is_integer_like(lo) || !expr_is_integer_like(hi))
            return false;

        expr_to_mpz(lo, imin);
        expr_to_mpz(hi, imax);
        return true;
    }
    return false;
}

/*
 * Build a multi-dimensional array of random integers.
 * dims is a list expression {n1, n2, ...}, dim_idx is current dimension.
 */
static Expr* random_array(const mpz_t imin, const mpz_t imax,
                          Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_integer_range(imin, imax);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        } else {
            elems[i] = random_array(imin, imax, dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * builtin_randominteger - implements RandomInteger[...]
 *
 * Forms:
 *   RandomInteger[]              -> 0 or 1
 *   RandomInteger[imax]          -> random in [0, imax]
 *   RandomInteger[{imin, imax}]  -> random in [imin, imax]
 *   RandomInteger[range, n]      -> flat list of n values
 *   RandomInteger[range, {n1, n2, ...}] -> nested array
 */
Expr* builtin_randominteger(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* RandomInteger[] -> 0 or 1 */
    if (argc == 0) {
        mpz_t lo, hi;
        mpz_init_set_si(lo, 0);
        mpz_init_set_si(hi, 1);
        Expr* result = random_integer_range(lo, hi);
        mpz_clear(lo);
        mpz_clear(hi);
        return result;
    }

    /* RandomInteger[range] or RandomInteger[range, dims] */
    if (argc == 1 || argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        mpz_t imin, imax;
        mpz_init(imin);
        mpz_init(imax);

        if (!parse_range(range_arg, imin, imax)) {
            mpz_clear(imin);
            mpz_clear(imax);
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (argc == 1) {
            /* Single random integer */
            Expr* result = random_integer_range(imin, imax);
            mpz_clear(imin);
            mpz_clear(imax);
            if (!result) return NULL;
            return result;
        }

        /* argc == 2: RandomInteger[range, n] or RandomInteger[range, {n1, n2, ...}] */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;

            /* Straight into a packed int64 buffer when the whole range is
             * machine-width, mirroring what RandomReal does. Every value is an
             * exact machine Integer by construction, so there is nothing to
             * decide per element -- and this skips an Expr allocation and a
             * pack_offer re-sniff per element, which was the difference between
             * 100 ns and 3 ns each. A wide range still builds the List below and
             * offers it afterwards, because random_integer_range can return a
             * BigInt and pack_offer must then decline for the whole list. */
            if (mpz_fits_slong_p(imin) && mpz_fits_slong_p(imax)) {
                int64_t lo = (int64_t)mpz_get_si(imin);
                int64_t hi = (int64_t)mpz_get_si(imax);
                int64_t span_m1;
                if (hi >= lo && !ci_sub_i64(hi, lo, &span_m1)) {
                    uint64_t span = (uint64_t)span_m1 + 1u;   /* inclusive */
                    int64_t* buf = NULL;
                    Expr* packed = ndbuild_open_i64((int64_t)n, &buf);
                    if (packed) {
                        /* The span is fixed for the whole array, so the shape of
                         * the draw is decided ONCE rather than re-tested per
                         * element. xs_below re-checks `m <= 1` and the
                         * power-of-two case on every call, which is nothing on
                         * its own and measurable across 10^7 of them. */
                        if ((span & (span - 1)) == 0) {
                            uint64_t mask = span - 1;    /* exact, no rejection */
                            for (size_t i = 0; i < n; i++)
                                buf[i] = lo + (int64_t)(xs_next() & mask);
                        } else {
                            for (size_t i = 0; i < n; i++)
                                buf[i] = lo + (int64_t)xs_below(span);
                        }
                        mpz_clear(imin);
                        mpz_clear(imax);
                        return packed;
                    }
                }
            }

            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) {
                mpz_clear(imin);
                mpz_clear(imax);
                return NULL;
            }
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_integer_range(imin, imax);
                if (!elems[i]) {
                    for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                    free(elems);
                    mpz_clear(imin);
                    mpz_clear(imax);
                    return NULL;
                }
            }
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
            free(elems);
            mpz_clear(imin);
            mpz_clear(imax);
            /* Not direct-built: random_integer_range goes through GMP and can
             * return a BigInt for a wide range, which pack_offer then declines
             * for the whole list. Offering after the fact keeps that decision in
             * one place. */
            return pack_offer(list);
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            /* Multi-dimensional array */
            /* Validate all dimensions are non-negative integers */
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    mpz_clear(imin);
                    mpz_clear(imax);
                    return NULL;
                }
            }
            Expr* result = random_array(imin, imax, dim_arg, 0);
            mpz_clear(imin);
            mpz_clear(imax);
            if (!result) return NULL;
            return pack_offer(result);
        }

        mpz_clear(imin);
        mpz_clear(imax);
        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * Read a plain numeric literal into a double. Handles the exact leaf forms
 * (Integer, Real, BigInt, Rational[n,d]) directly, without any numeric
 * evaluation. Returns true on success, false otherwise.
 */
static bool literal_to_real(Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) {
        *out = (double)e->data.integer;
        return true;
    }
    if (e->type == EXPR_REAL) {
        *out = e->data.real;
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        *out = mpz_get_d(e->data.bigint);
        return true;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        return true;
    }
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        *out = (double)n / (double)d;
        return true;
    }
    return false;
}

/*
 * Convert a numeric expression to a double value.
 *
 * Fast path: exact literals via literal_to_real. Slow path: any other
 * expression that N[] can reduce to a machine number — symbolic constants
 * and constant arithmetic such as Pi, -Pi, Sqrt[2], or E/2 — so that range
 * bounds like RandomReal[{-Pi, Pi}] evaluate. Returns false only when the
 * expression is genuinely non-numeric (e.g. a free symbol x).
 */
static bool expr_to_real(Expr* e, double* out) {
    if (literal_to_real(e, out)) return true;

    Expr* num = numericalize(e, numeric_machine_spec());
    if (!num) return false;
    bool ok = literal_to_real(num, out);
    expr_free(num);
    return ok;
}

/*
 * Generate a single uniform random real in [0, 1).
 * Uses 53 bits of randomness from GMP for full double precision.
 */
static double random_uniform_01(void) {
    ensure_rand_init();
    /* Top 53 bits scaled by 2^-53: uniform on [0, 1) with a full double mantissa,
     * the same 53 bits of randomness the GMP version delivered, at 2.4 ns instead
     * of 117 ns. The high bits are used because xoshiro's low bits are the
     * weakest -- that is the documented way to take a double from this family. */
    return (double)(xs_next() >> 11) * 0x1.0p-53;
}

/*
 * Generate a single random real in [xmin, xmax).
 */
static double random_real_value(double xmin, double xmax) {
    double u = random_uniform_01();
    return xmin + u * (xmax - xmin);
}

static Expr* random_real_range(double xmin, double xmax) {
    return expr_new_real(random_real_value(xmin, xmax));
}

/*
 * Parse a real-valued range argument into xmin, xmax.
 * Supports:
 *   - Numeric atom: range is [0, val]
 *   - List[xmin, xmax]: range is [xmin, xmax]
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_real_range(Expr* arg, double* xmin, double* xmax) {
    double val;
    if (expr_to_real(arg, &val)) {
        *xmin = 0.0;
        *xmax = val;
        return true;
    }
    /* Check for List[xmin, xmax] */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List &&
        arg->data.function.arg_count == 2) {

        Expr* lo = arg->data.function.args[0];
        Expr* hi = arg->data.function.args[1];

        if (!expr_to_real(lo, xmin) || !expr_to_real(hi, xmax))
            return false;

        return true;
    }
    return false;
}

/*
 * Build a multi-dimensional array of random reals.
 * dims is a list expression {n1, n2, ...}, dim_idx is current dimension.
 */
static Expr* random_real_array(double xmin, double xmax,
                               Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_real_range(xmin, xmax);
        } else {
            elems[i] = random_real_array(xmin, xmax, dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * Detect a trailing option of the form Rule[WorkingPrecision, n].
 *
 * If arg is `WorkingPrecision -> MachinePrecision` (or omitted), *is_machine
 * is set to true and *bits to 0. If arg is `WorkingPrecision -> digits`
 * where `digits` is a positive numeric literal, *is_machine is set to false
 * and *bits is set to the corresponding MPFR precision (in binary bits).
 * Returns true on a successful match (machine or MPFR), false otherwise.
 *
 * Without USE_MPFR support, a numeric digit count still parses but degrades
 * to machine precision so existing programs keep working — a one-shot
 * stderr warning communicates the loss.
 */
static bool parse_working_precision_opt(Expr* arg, long* bits, bool* is_machine) {
    if (!arg || arg->type != EXPR_FUNCTION) return false;
    Expr* head = arg->data.function.head;
    if (head->type != EXPR_SYMBOL || head->data.symbol.name != SYM_Rule) return false;
    if (arg->data.function.arg_count != 2) return false;
    Expr* key = arg->data.function.args[0];
    if (key->type != EXPR_SYMBOL || key->data.symbol.name != SYM_WorkingPrecision) return false;

    Expr* val = arg->data.function.args[1];
    if (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_MachinePrecision) {
        *bits = 0;
        *is_machine = true;
        return true;
    }

    double digits = 0.0;
    int64_t rn, rd;
    if (val->type == EXPR_INTEGER) digits = (double)val->data.integer;
    else if (val->type == EXPR_REAL) digits = val->data.real;
    else if (is_rational(val, &rn, &rd)) digits = (double)rn / (double)rd;
    else return false;
    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    if (digits <= NUMERIC_MACHINE_PRECISION_DIGITS) {
        /* Mathematica treats requests at or below machine precision as
         * machine: the doubles path is exact enough. */
        *bits = 0;
        *is_machine = true;
    } else {
        *bits = numeric_digits_to_bits(digits);
        *is_machine = false;
    }
    return true;
#else
    static bool warned = false;
    if (!warned) {
        fprintf(stderr,
                "WorkingPrecision: arbitrary precision unavailable "
                "(USE_MPFR=0); using machine precision.\n");
        warned = true;
    }
    *bits = 0;
    *is_machine = true;
    return true;
#endif
}

#ifdef USE_MPFR
/*
 * Extract an MPFR approximation of `e` into `(re, im)`, numericalizing
 * symbolic-but-numeric bounds (Pi, -Pi, Sqrt[2], ...) when the direct
 * get_approx_mpfr extraction fails. `re` and `im` must already be
 * initialized to the target precision. Returns false only when `e` is
 * genuinely non-numeric.
 */
static bool random_get_approx_mpfr(Expr* e, mpfr_t re, mpfr_t im) {
    if (get_approx_mpfr(e, re, im, NULL)) return true;

    NumericSpec spec;
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = (long)mpfr_get_prec(re);
    spec.preserve_inexact = false;
    Expr* num = numericalize(e, spec);
    if (!num) return false;
    bool ok = get_approx_mpfr(num, re, im, NULL);
    expr_free(num);
    return ok;
}

/*
 * Generate a single MPFR random real in [xmin, xmax) at the given precision.
 * xmin and xmax must already be initialized MPFR values.
 */
static Expr* random_real_range_mpfr(const mpfr_t xmin, const mpfr_t xmax, long bits) {
    ensure_rand_init();
    mpfr_t u, out;
    mpfr_init2(u, bits);
    mpfr_init2(out, bits);
    mpfr_urandomb(u, g_rand_state);          /* uniform in [0, 1) */
    mpfr_sub(out, xmax, xmin, MPFR_RNDN);
    mpfr_mul(out, out, u, MPFR_RNDN);
    mpfr_add(out, out, xmin, MPFR_RNDN);
    mpfr_clear(u);
    return expr_new_mpfr_move(out);
}

/*
 * Parse a real-valued range argument into mpfr_t xmin, xmax. The output
 * MPFR values must already be initialized to the desired precision.
 *
 * Supports the same forms as parse_real_range, but bound coercion goes
 * through get_approx_mpfr so exact rationals are preserved at full
 * working precision.
 */
static bool parse_real_range_mpfr(Expr* arg, mpfr_t xmin, mpfr_t xmax) {
    long bits = mpfr_get_prec(xmin);
    mpfr_t scratch_im;
    mpfr_init2(scratch_im, bits);

    bool ok;
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List &&
        arg->data.function.arg_count == 2) {
        ok = random_get_approx_mpfr(arg->data.function.args[0], xmin, scratch_im)
             && mpfr_zero_p(scratch_im);
        if (ok) {
            ok = random_get_approx_mpfr(arg->data.function.args[1], xmax, scratch_im)
                 && mpfr_zero_p(scratch_im);
        }
    } else {
        /* Single value: range is [0, val]. */
        mpfr_set_zero(xmin, +1);
        ok = random_get_approx_mpfr(arg, xmax, scratch_im)
             && mpfr_zero_p(scratch_im);
    }

    mpfr_clear(scratch_im);
    return ok;
}

/*
 * Build a multi-dimensional MPFR array of random reals.
 */
static Expr* random_real_array_mpfr(const mpfr_t xmin, const mpfr_t xmax,
                                    long bits, Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count) return NULL;
    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0) return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);
    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_real_range_mpfr(xmin, xmax, bits);
        } else {
            elems[i] = random_real_array_mpfr(xmin, xmax, bits, dims, dim_idx + 1);
        }
        if (!elems[i]) {
            for (size_t j = 0; j < i; j++) expr_free(elems[j]);
            free(elems);
            return NULL;
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * MPFR-precision RandomReal[...] body. effective_argc excludes the trailing
 * WorkingPrecision option.
 */
static Expr* randomreal_mpfr(Expr* res, size_t effective_argc, long bits) {
    mpfr_t xmin, xmax;
    mpfr_init2(xmin, bits);
    mpfr_init2(xmax, bits);

    if (effective_argc == 0) {
        mpfr_set_zero(xmin, +1);
        mpfr_set_si(xmax, 1, MPFR_RNDN);
        Expr* r = random_real_range_mpfr(xmin, xmax, bits);
        mpfr_clear(xmin); mpfr_clear(xmax);
        return r;
    }

    if (effective_argc == 1 || effective_argc == 2) {
        Expr* range_arg = res->data.function.args[0];
        if (!parse_real_range_mpfr(range_arg, xmin, xmax)) {
            mpfr_clear(xmin); mpfr_clear(xmax);
            return NULL;
        }

        if (effective_argc == 1) {
            Expr* r = random_real_range_mpfr(xmin, xmax, bits);
            mpfr_clear(xmin); mpfr_clear(xmax);
            return r;
        }

        Expr* dim_arg = res->data.function.args[1];
        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) { mpfr_clear(xmin); mpfr_clear(xmax); return NULL; }
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_real_range_mpfr(xmin, xmax, bits);
            }
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
            free(elems);
            mpfr_clear(xmin); mpfr_clear(xmax);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    mpfr_clear(xmin); mpfr_clear(xmax);
                    return NULL;
                }
            }
            Expr* result = random_real_array_mpfr(xmin, xmax, bits, dim_arg, 0);
            mpfr_clear(xmin); mpfr_clear(xmax);
            return result;
        }

        mpfr_clear(xmin); mpfr_clear(xmax);
        return NULL;
    }

    mpfr_clear(xmin); mpfr_clear(xmax);
    return NULL;
}
#endif /* USE_MPFR */

/* Machine-precision RandomReal[...] body. effective_argc excludes any trailing
 * WorkingPrecision option. */
static Expr* randomreal_machine(Expr* res, size_t effective_argc) {
    /* RandomReal[] -> random real in [0, 1) */
    if (effective_argc == 0) {
        return random_real_range(0.0, 1.0);
    }

    /* RandomReal[range] or RandomReal[range, dims] */
    if (effective_argc == 1 || effective_argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        double xmin, xmax;
        if (!parse_real_range(range_arg, &xmin, &xmax)) {
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (effective_argc == 1) {
            return random_real_range(xmin, xmax);
        }

        /* effective_argc == 2: RandomReal[range, n] or RandomReal[range, {n1, n2, ...}] */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            /* Straight into a packed buffer when it is big enough to be worth
             * it: the values are machine reals by construction, so there is
             * nothing to decide per element. */
            double* buf = NULL;
            Expr* packed = ndbuild_open_f64((int64_t)n, &buf);
            if (packed) {
                for (size_t i = 0; i < n; i++) buf[i] = random_real_value(xmin, xmax);
                return packed;
            }
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) return NULL;
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_real_range(xmin, xmax);
            }
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
            free(elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            /* Multi-dimensional array */
            size_t rank = dim_arg->data.function.arg_count;
            for (size_t i = 0; i < rank; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    return NULL;
                }
            }
            if (rank >= 1 && rank < NDARRAY_MAX_RANK) {
                int64_t dims[NDARRAY_MAX_RANK];
                int64_t total = 1;
                for (size_t i = 0; i < rank; i++) {
                    dims[i] = dim_arg->data.function.args[i]->data.integer;
                    total *= dims[i];
                }
                void* vbuf = NULL;
                Expr* packed = ndbuild_open((int)rank, dims, NDT_FLOAT64, &vbuf);
                if (packed) {
                    double* buf = (double*)vbuf;
                    /* Row-major fill order is the same order random_real_array
                     * visits the leaves in, so a seeded stream gives the same
                     * array either way. */
                    for (int64_t i = 0; i < total; i++)
                        buf[i] = random_real_value(xmin, xmax);
                    return packed;
                }
            }
            Expr* result = random_real_array(xmin, xmax, dim_arg, 0);
            if (!result) return NULL;
            return result;
        }

        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * builtin_randomreal - implements RandomReal[...]
 *
 * Forms:
 *   RandomReal[]                                     -> random real in [0, 1)
 *   RandomReal[xmax]                                 -> random real in [0, xmax)
 *   RandomReal[{xmin, xmax}]                         -> random real in [xmin, xmax)
 *   RandomReal[range, n]                             -> flat list of n values
 *   RandomReal[range, {n1, n2, ...}]                 -> nested array
 *   RandomReal[spec, WorkingPrecision -> n]          -> n-digit MPFR result(s)
 *
 * WorkingPrecision may be passed as the last argument of any of the
 * positional forms. A digit count > MachinePrecision triggers MPFR-backed
 * generation; MachinePrecision (or no option) takes the doubles path.
 */
Expr* builtin_randomreal(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    long bits = 0;
    bool is_machine = true;
    if (argc > 0) {
        Expr* last = res->data.function.args[argc - 1];
        if (parse_working_precision_opt(last, &bits, &is_machine)) {
            argc--; /* Peel off the option from the positional count. */
        }
    }

#ifdef USE_MPFR
    if (!is_machine) return randomreal_mpfr(res, argc, bits);
#else
    (void)bits;
#endif
    return randomreal_machine(res, argc);
}

/*
 * builtin_seedrandom - implements SeedRandom[n] and SeedRandom[]
 */
Expr* builtin_seedrandom(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        /* Re-seed from system entropy */
        ensure_rand_init();
        seed_all((uint64_t)time(NULL) ^ (uint64_t)clock());
        return expr_new_symbol(SYM_Null);
    }

    if (argc == 1) {
        Expr* arg = res->data.function.args[0];
        if (arg->type == EXPR_INTEGER) {
            ensure_rand_init();
            seed_all((uint64_t)arg->data.integer);
            return expr_new_symbol(SYM_Null);
        }
        if (arg->type == EXPR_BIGINT) {
            ensure_rand_init();
            gmp_randseed(g_rand_state, arg->data.bigint);
            /* Derive the machine seed from the same bignum, so SeedRandom[huge]
             * is as reproducible as SeedRandom[small]. */
            xs_seed((uint64_t)mpz_get_ui(arg->data.bigint) ^
                    ((uint64_t)mpz_sizeinbase(arg->data.bigint, 2) << 32));
            return expr_new_symbol(SYM_Null);
        }
        return NULL; /* Non-integer seed */
    }

    return NULL; /* Too many arguments */
}

/*
 * Generate a single random complex number with real part in [re_min, re_max)
 * and imaginary part in [im_min, im_max).
 * Returns a new Complex[re, im] expression.
 */
static Expr* random_complex_range(double re_min, double re_max,
                                   double im_min, double im_max) {
    double u1 = random_uniform_01();
    double u2 = random_uniform_01();
    double re = re_min + u1 * (re_max - re_min);
    double im = im_min + u2 * (im_max - im_min);
    return make_complex(expr_new_real(re), expr_new_real(im));
}

/*
 * Convert a numeric expression to its real and imaginary double parts.
 *
 * Handles a Complex[] atom, a plain real literal (imag part 0), and — via
 * numericalize — symbolic-but-numeric corners such as `-Pi - I` (a Plus that
 * only collapses to a Complex once Pi is numeric) or Exp[I Pi/4]. Returns
 * false only when the expression is genuinely non-numeric.
 */
static bool expr_to_complex_parts(Expr* e, double* re, double* im) {
    Expr *rp, *ip;
    if (is_complex(e, &rp, &ip))
        return expr_to_real(rp, re) && expr_to_real(ip, im);
    if (expr_to_real(e, re)) {
        *im = 0.0;
        return true;
    }
    /* Slow path: reduce symbolic expressions like `-Pi - I` to a number. */
    Expr* num = numericalize(e, numeric_machine_spec());
    if (!num) return false;
    bool ok = false;
    if (is_complex(num, &rp, &ip))
        ok = expr_to_real(rp, re) && expr_to_real(ip, im);
    else if (literal_to_real(num, re)) {
        *im = 0.0;
        ok = true;
    }
    expr_free(num);
    return ok;
}

/*
 * Parse a complex-valued range argument into corner coordinates.
 * A complex range is a rectangle in the complex plane defined by two corners.
 * Supports:
 *   - Complex atom: range is [0+0i, val]
 *   - Real/Integer atom: range is [0+0i, val+0i]
 *   - List[zmin, zmax]: range is the rectangle with corners zmin, zmax
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_complex_range(Expr* arg,
                                double* re_min, double* re_max,
                                double* im_min, double* im_max) {
    /* Single value: upper corner, lower corner is origin. */
    if (arg->type != EXPR_FUNCTION ||
        arg->data.function.head->type != EXPR_SYMBOL ||
        arg->data.function.head->data.symbol.name != SYM_List ||
        arg->data.function.arg_count != 2) {
        *re_min = 0.0; *im_min = 0.0;
        return expr_to_complex_parts(arg, re_max, im_max);
    }

    /* List[zmin, zmax] rectangle. */
    return expr_to_complex_parts(arg->data.function.args[0], re_min, im_min)
        && expr_to_complex_parts(arg->data.function.args[1], re_max, im_max);
}

/*
 * Build a multi-dimensional array of random complex numbers.
 */
static Expr* random_complex_array(double re_min, double re_max,
                                   double im_min, double im_max,
                                   Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_complex_range(re_min, re_max, im_min, im_max);
        } else {
            elems[i] = random_complex_array(re_min, re_max, im_min, im_max,
                                            dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

#ifdef USE_MPFR
/*
 * Generate a single MPFR random complex number whose real part is uniform
 * in [re_min, re_max) and whose imaginary part is uniform in [im_min, im_max),
 * each at the given working precision. All bound mpfr_t inputs must be
 * initialized.
 */
static Expr* random_complex_range_mpfr(const mpfr_t re_min, const mpfr_t re_max,
                                       const mpfr_t im_min, const mpfr_t im_max,
                                       long bits) {
    ensure_rand_init();
    mpfr_t u1, u2, re, im;
    mpfr_init2(u1, bits);
    mpfr_init2(u2, bits);
    mpfr_init2(re, bits);
    mpfr_init2(im, bits);
    mpfr_urandomb(u1, g_rand_state);
    mpfr_urandomb(u2, g_rand_state);

    mpfr_sub(re, re_max, re_min, MPFR_RNDN);
    mpfr_mul(re, re, u1, MPFR_RNDN);
    mpfr_add(re, re, re_min, MPFR_RNDN);

    mpfr_sub(im, im_max, im_min, MPFR_RNDN);
    mpfr_mul(im, im, u2, MPFR_RNDN);
    mpfr_add(im, im, im_min, MPFR_RNDN);

    mpfr_clear(u1);
    mpfr_clear(u2);
    return make_complex(expr_new_mpfr_move(re), expr_new_mpfr_move(im));
}

/*
 * Parse a complex-valued range argument into MPFR-precision rectangle
 * corners. The corner mpfr_t values must already be initialized at the
 * desired working precision; this routine does not change their precision.
 *
 * Supports the same forms as parse_complex_range but uses get_approx_mpfr
 * so that integer / rational / Complex bounds keep full working precision.
 */
static bool parse_complex_range_mpfr(Expr* arg,
                                     mpfr_t re_min, mpfr_t re_max,
                                     mpfr_t im_min, mpfr_t im_max) {
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List &&
        arg->data.function.arg_count == 2) {
        if (!random_get_approx_mpfr(arg->data.function.args[0], re_min, im_min))
            return false;
        if (!random_get_approx_mpfr(arg->data.function.args[1], re_max, im_max))
            return false;
        return true;
    }

    /* Single value: rectangle is from origin to arg. */
    mpfr_set_zero(re_min, +1);
    mpfr_set_zero(im_min, +1);
    return random_get_approx_mpfr(arg, re_max, im_max);
}

/*
 * Build a multi-dimensional MPFR array of random complex numbers.
 */
static Expr* random_complex_array_mpfr(const mpfr_t re_min, const mpfr_t re_max,
                                       const mpfr_t im_min, const mpfr_t im_max,
                                       long bits, Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count) return NULL;
    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0) return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);
    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_complex_range_mpfr(re_min, re_max, im_min, im_max, bits);
        } else {
            elems[i] = random_complex_array_mpfr(re_min, re_max, im_min, im_max,
                                                  bits, dims, dim_idx + 1);
        }
        if (!elems[i]) {
            for (size_t j = 0; j < i; j++) expr_free(elems[j]);
            free(elems);
            return NULL;
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * MPFR-precision RandomComplex[...] body. effective_argc excludes any trailing
 * WorkingPrecision option.
 */
static Expr* randomcomplex_mpfr(Expr* res, size_t effective_argc, long bits) {
    mpfr_t re_min, re_max, im_min, im_max;
    mpfr_init2(re_min, bits);
    mpfr_init2(re_max, bits);
    mpfr_init2(im_min, bits);
    mpfr_init2(im_max, bits);

    if (effective_argc == 0) {
        mpfr_set_zero(re_min, +1);
        mpfr_set_zero(im_min, +1);
        mpfr_set_si(re_max, 1, MPFR_RNDN);
        mpfr_set_si(im_max, 1, MPFR_RNDN);
        Expr* r = random_complex_range_mpfr(re_min, re_max, im_min, im_max, bits);
        mpfr_clear(re_min); mpfr_clear(re_max);
        mpfr_clear(im_min); mpfr_clear(im_max);
        return r;
    }

    if (effective_argc == 1 || effective_argc == 2) {
        Expr* range_arg = res->data.function.args[0];
        if (!parse_complex_range_mpfr(range_arg, re_min, re_max, im_min, im_max)) {
            mpfr_clear(re_min); mpfr_clear(re_max);
            mpfr_clear(im_min); mpfr_clear(im_max);
            return NULL;
        }

        if (effective_argc == 1) {
            Expr* r = random_complex_range_mpfr(re_min, re_max, im_min, im_max, bits);
            mpfr_clear(re_min); mpfr_clear(re_max);
            mpfr_clear(im_min); mpfr_clear(im_max);
            return r;
        }

        Expr* dim_arg = res->data.function.args[1];
        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) {
                mpfr_clear(re_min); mpfr_clear(re_max);
                mpfr_clear(im_min); mpfr_clear(im_max);
                return NULL;
            }
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_complex_range_mpfr(re_min, re_max, im_min, im_max, bits);
            }
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
            free(elems);
            mpfr_clear(re_min); mpfr_clear(re_max);
            mpfr_clear(im_min); mpfr_clear(im_max);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    mpfr_clear(re_min); mpfr_clear(re_max);
                    mpfr_clear(im_min); mpfr_clear(im_max);
                    return NULL;
                }
            }
            Expr* result = random_complex_array_mpfr(re_min, re_max, im_min, im_max,
                                                      bits, dim_arg, 0);
            mpfr_clear(re_min); mpfr_clear(re_max);
            mpfr_clear(im_min); mpfr_clear(im_max);
            return result;
        }

        mpfr_clear(re_min); mpfr_clear(re_max);
        mpfr_clear(im_min); mpfr_clear(im_max);
        return NULL;
    }

    mpfr_clear(re_min); mpfr_clear(re_max);
    mpfr_clear(im_min); mpfr_clear(im_max);
    return NULL;
}
#endif /* USE_MPFR */

/* Machine-precision RandomComplex[...] body. effective_argc excludes any
 * trailing WorkingPrecision option. */
static Expr* randomcomplex_machine(Expr* res, size_t effective_argc) {
    /* RandomComplex[] -> random complex in [0,1]+[0,1]i */
    if (effective_argc == 0) {
        return random_complex_range(0.0, 1.0, 0.0, 1.0);
    }

    /* RandomComplex[range] or RandomComplex[range, dims] */
    if (effective_argc == 1 || effective_argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        double re_min, re_max, im_min, im_max;
        if (!parse_complex_range(range_arg, &re_min, &re_max, &im_min, &im_max)) {
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (effective_argc == 1) {
            return random_complex_range(re_min, re_max, im_min, im_max);
        }

        /* effective_argc == 2 */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) return NULL;
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_complex_range(re_min, re_max, im_min, im_max);
            }
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
            free(elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            /* Multi-dimensional array */
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    return NULL;
                }
            }
            Expr* result = random_complex_array(re_min, re_max, im_min, im_max,
                                                dim_arg, 0);
            if (!result) return NULL;
            return result;
        }

        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * builtin_randomcomplex - implements RandomComplex[...]
 *
 * Forms:
 *   RandomComplex[]                                  -> random complex in [0,1]+[0,1]i
 *   RandomComplex[zmax]                              -> random complex in rectangle [0, zmax]
 *   RandomComplex[{zmin, zmax}]                      -> random complex in rectangle [zmin, zmax]
 *   RandomComplex[range, n]                          -> flat list of n values
 *   RandomComplex[range, {n1, n2, ...}]              -> nested array
 *   RandomComplex[spec, WorkingPrecision -> n]       -> n-digit MPFR result(s)
 *
 * WorkingPrecision may be passed as the last argument of any of the
 * positional forms. A digit count > MachinePrecision triggers MPFR-backed
 * generation; MachinePrecision (or no option) takes the doubles path.
 */
Expr* builtin_randomcomplex(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    long bits = 0;
    bool is_machine = true;
    if (argc > 0) {
        Expr* last = res->data.function.args[argc - 1];
        if (parse_working_precision_opt(last, &bits, &is_machine)) {
            argc--;
        }
    }

#ifdef USE_MPFR
    if (!is_machine) return randomcomplex_mpfr(res, argc, bits);
#else
    (void)bits;
#endif
    return randomcomplex_machine(res, argc);
}

/*
 * Pick a random index in [0, n) using uniform distribution.
 */
static size_t random_index(size_t n) {
    ensure_rand_init();
    return (size_t)xs_below((uint64_t)n);   /* rejection: no modulo bias */
}

/*
 * Pick a weighted random index given an array of cumulative weights.
 * cum_weights[i] is the sum of weights[0..i].
 * total is cum_weights[count-1].
 * Uses binary search for efficiency.
 */
static size_t weighted_random_index(double* cum_weights, size_t count, double total) {
    double u = random_uniform_01() * total;
    /* Binary search for the first index where cum_weights[i] > u */
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (cum_weights[mid] <= u)
            lo = mid + 1;
        else
            hi = mid;
    }
    /* Clamp to valid range */
    if (lo >= count) lo = count - 1;
    return lo;
}

/*
 * Build a multi-dimensional array of random choices from a list.
 * Uses uniform selection (each element equally likely).
 */
static Expr* random_choice_array(Expr** choices, size_t choice_count,
                                 Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            size_t idx = random_index(choice_count);
            elems[i] = expr_copy(choices[idx]);
        } else {
            elems[i] = random_choice_array(choices, choice_count,
                                           dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * Build a multi-dimensional array of weighted random choices.
 */
static Expr* weighted_choice_array(Expr** choices, size_t choice_count,
                                   double* cum_weights, double total,
                                   Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            size_t idx = weighted_random_index(cum_weights, choice_count, total);
            elems[i] = expr_copy(choices[idx]);
        } else {
            elems[i] = weighted_choice_array(choices, choice_count,
                                             cum_weights, total,
                                             dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * Check if an expression is Rule[wlist, elist] (i.e., wlist -> elist).
 * If so, set *wlist and *elist to the left and right sides.
 */
static bool is_rule(Expr* e, Expr** wlist, Expr** elist) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rule &&
        e->data.function.arg_count == 2) {
        *wlist = e->data.function.args[0];
        *elist = e->data.function.args[1];
        return true;
    }
    return false;
}

/*
 * Check if an expression is a List with at least one element.
 * If so, set *args and *count.
 */
static bool is_nonempty_list(Expr* e, Expr*** args, size_t* count) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_List &&
        e->data.function.arg_count > 0) {
        *args = e->data.function.args;
        *count = e->data.function.arg_count;
        return true;
    }
    return false;
}

/*
 * builtin_randomchoice - implements RandomChoice[...]
 *
 * Forms:
 *   RandomChoice[{e1, e2, ...}]                     -> uniform random choice
 *   RandomChoice[list, n]                            -> list of n uniform choices
 *   RandomChoice[list, {n1, n2, ...}]                -> array of uniform choices
 *   RandomChoice[{w1,w2,...}->{e1,e2,...}]            -> weighted random choice
 *   RandomChoice[wlist->elist, n]                    -> list of n weighted choices
 *   RandomChoice[wlist->elist, {n1, n2, ...}]        -> array of weighted choices
 */
Expr* builtin_randomchoice(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 1 || argc > 2) return NULL;

    Expr* first_arg = res->data.function.args[0];
    Expr *wlist_expr, *elist_expr;

    /* Check if first arg is Rule[wlist, elist] (weighted form) */
    if (is_rule(first_arg, &wlist_expr, &elist_expr)) {
        /* Weighted random choice */
        Expr** weights_args;
        size_t weights_count;
        Expr** elems_args;
        size_t elems_count;

        if (!is_nonempty_list(wlist_expr, &weights_args, &weights_count))
            return NULL;
        if (!is_nonempty_list(elist_expr, &elems_args, &elems_count))
            return NULL;
        if (weights_count != elems_count)
            return NULL;

        /* Convert weights to doubles and compute cumulative weights */
        double* cum_weights = malloc(sizeof(double) * weights_count);
        if (!cum_weights) return NULL;

        double total = 0.0;
        bool valid = true;
        for (size_t i = 0; i < weights_count; i++) {
            double w;
            if (!expr_to_real(weights_args[i], &w) || w < 0.0) {
                valid = false;
                break;
            }
            total += w;
            cum_weights[i] = total;
        }

        if (!valid || total <= 0.0) {
            free(cum_weights);
            return NULL;
        }

        if (argc == 1) {
            /* Single weighted choice */
            size_t idx = weighted_random_index(cum_weights, elems_count, total);
            free(cum_weights);
            return expr_copy(elems_args[idx]);
        }

        /* argc == 2 */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            size_t n = (size_t)dim_arg->data.integer;
            Expr** result_elems = malloc(sizeof(Expr*) * n);
            if (!result_elems) { free(cum_weights); return NULL; }
            for (size_t i = 0; i < n; i++) {
                size_t idx = weighted_random_index(cum_weights, elems_count, total);
                result_elems[i] = expr_copy(elems_args[idx]);
            }
            free(cum_weights);
            Expr* list = expr_new_function(expr_new_symbol(SYM_List), result_elems, n);
            free(result_elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol.name == SYM_List) {
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    free(cum_weights);
                    return NULL;
                }
            }
            Expr* result = weighted_choice_array(elems_args, elems_count,
                                                 cum_weights, total,
                                                 dim_arg, 0);
            free(cum_weights);
            if (!result) return NULL;
            return result;
        }

        free(cum_weights);
        return NULL;
    }

    /* Uniform random choice: first arg must be a non-empty list */
    Expr** choices;
    size_t choice_count;
    if (!is_nonempty_list(first_arg, &choices, &choice_count))
        return NULL;

    if (argc == 1) {
        /* Single uniform choice */
        size_t idx = random_index(choice_count);
        return expr_copy(choices[idx]);
    }

    /* argc == 2 */
    Expr* dim_arg = res->data.function.args[1];

    if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
        size_t n = (size_t)dim_arg->data.integer;
        Expr** elems = malloc(sizeof(Expr*) * n);
        if (!elems) return NULL;
        for (size_t i = 0; i < n; i++) {
            size_t idx = random_index(choice_count);
            elems[i] = expr_copy(choices[idx]);
        }
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
        free(elems);
        return list;
    }

    if (dim_arg->type == EXPR_FUNCTION &&
        dim_arg->data.function.head->type == EXPR_SYMBOL &&
        dim_arg->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
            Expr* d = dim_arg->data.function.args[i];
            if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                return NULL;
            }
        }
        Expr* result = random_choice_array(choices, choice_count, dim_arg, 0);
        if (!result) return NULL;
        return result;
    }

    return NULL;
}

/*
 * Fisher-Yates partial shuffle: select n elements from indices [0, total)
 * without replacement. Returns an array of n selected indices.
 * Caller must free the returned array.
 */
static size_t* fisher_yates_sample(size_t total, size_t n) {
    size_t* indices = malloc(sizeof(size_t) * total);
    if (!indices) return NULL;
    for (size_t i = 0; i < total; i++) indices[i] = i;

    /* Partial shuffle: only need first n elements */
    for (size_t i = 0; i < n; i++) {
        size_t remaining = total - i;
        size_t j = i + random_index(remaining);
        /* Swap indices[i] and indices[j] */
        size_t tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    /* Return only the first n indices */
    size_t* result = malloc(sizeof(size_t) * n);
    if (!result) { free(indices); return NULL; }
    for (size_t i = 0; i < n; i++) result[i] = indices[i];
    free(indices);
    return result;
}

/*
 * Weighted sampling without replacement.
 * Selects n indices from [0, count) with probabilities proportional to weights.
 * Uses the "remove and renormalize" approach.
 * Returns an array of n selected indices. Caller must free.
 */
static size_t* weighted_sample_without_replacement(double* weights, size_t count, size_t n) {
    /* Make a mutable copy of weights */
    double* w = malloc(sizeof(double) * count);
    if (!w) return NULL;
    for (size_t i = 0; i < count; i++) w[i] = weights[i];

    size_t* result = malloc(sizeof(size_t) * n);
    if (!result) { free(w); return NULL; }

    for (size_t k = 0; k < n; k++) {
        /* Compute total remaining weight */
        double total = 0.0;
        for (size_t i = 0; i < count; i++) total += w[i];
        if (total <= 0.0) { free(w); free(result); return NULL; }

        /* Build cumulative weights for active elements */
        double u = random_uniform_01() * total;
        double cum = 0.0;
        size_t chosen = count - 1; /* fallback */
        for (size_t i = 0; i < count; i++) {
            if (w[i] <= 0.0) continue;
            cum += w[i];
            if (cum > u) {
                chosen = i;
                break;
            }
        }

        result[k] = chosen;
        w[chosen] = 0.0; /* Remove from future selection */
    }

    free(w);
    return result;
}

/*
 * Check if an expression is UpTo[n] and extract n.
 * Returns true if it matches, with *val set to n.
 */
static bool is_upto(Expr* e, int64_t* val) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_UpTo &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER) {
        *val = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

/*
 * builtin_randomsample - implements RandomSample[...]
 *
 * Forms:
 *   RandomSample[{e1, e2, ...}]                        -> random permutation
 *   RandomSample[{e1, e2, ...}, n]                     -> sample of n elements (no replacement)
 *   RandomSample[{e1, e2, ...}, UpTo[n]]               -> sample of min(n, len) elements
 *   RandomSample[{w1,...}->{e1,...}]                    -> weighted random permutation
 *   RandomSample[{w1,...}->{e1,...}, n]                 -> weighted sample of n elements
 *   RandomSample[{w1,...}->{e1,...}, UpTo[n]]           -> weighted sample of min(n, len)
 */
Expr* builtin_randomsample(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 1 || argc > 2) return NULL;

    Expr* first_arg = res->data.function.args[0];
    Expr *wlist_expr, *elist_expr;

    /* Determine the sample size */
    size_t sample_n;
    bool have_n = false;

    /* Check if first arg is Rule[wlist, elist] (weighted form) */
    if (is_rule(first_arg, &wlist_expr, &elist_expr)) {
        /* Weighted random sample */
        Expr** weights_args;
        size_t weights_count;
        Expr** elems_args;
        size_t elems_count;

        if (!is_nonempty_list(wlist_expr, &weights_args, &weights_count))
            return NULL;
        if (!is_nonempty_list(elist_expr, &elems_args, &elems_count))
            return NULL;
        if (weights_count != elems_count)
            return NULL;

        /* Parse weights to doubles */
        double* weights = malloc(sizeof(double) * weights_count);
        if (!weights) return NULL;

        bool valid = true;
        for (size_t i = 0; i < weights_count; i++) {
            if (!expr_to_real(weights_args[i], &weights[i]) || weights[i] < 0.0) {
                valid = false;
                break;
            }
        }
        if (!valid) { free(weights); return NULL; }

        if (argc == 2) {
            Expr* n_arg = res->data.function.args[1];
            int64_t upto_val;
            if (n_arg->type == EXPR_INTEGER && n_arg->data.integer >= 0) {
                sample_n = (size_t)n_arg->data.integer;
                if (sample_n > elems_count) {
                    free(weights);
                    return NULL; /* Cannot sample more than available */
                }
                have_n = true;
            } else if (is_upto(n_arg, &upto_val) && upto_val >= 0) {
                sample_n = (size_t)upto_val;
                if (sample_n > elems_count) sample_n = elems_count;
                have_n = true;
            } else {
                free(weights);
                return NULL;
            }
        } else {
            /* Full permutation */
            sample_n = elems_count;
            have_n = true;
        }

        size_t* selected = weighted_sample_without_replacement(weights, elems_count, sample_n);
        free(weights);
        if (!selected) return NULL;

        Expr** result_elems = malloc(sizeof(Expr*) * sample_n);
        if (!result_elems) { free(selected); return NULL; }
        for (size_t i = 0; i < sample_n; i++) {
            result_elems[i] = expr_copy(elems_args[selected[i]]);
        }
        free(selected);

        Expr* list = expr_new_function(expr_new_symbol(SYM_List), result_elems, sample_n);
        free(result_elems);
        return list;
    }

    /* Uniform random sample: first arg must be a non-empty list */
    Expr** choices;
    size_t choice_count;
    if (!is_nonempty_list(first_arg, &choices, &choice_count))
        return NULL;

    if (argc == 2) {
        Expr* n_arg = res->data.function.args[1];
        int64_t upto_val;
        if (n_arg->type == EXPR_INTEGER && n_arg->data.integer >= 0) {
            sample_n = (size_t)n_arg->data.integer;
            if (sample_n > choice_count) return NULL; /* Can't sample more than available */
            have_n = true;
        } else if (is_upto(n_arg, &upto_val) && upto_val >= 0) {
            sample_n = (size_t)upto_val;
            if (sample_n > choice_count) sample_n = choice_count;
            have_n = true;
        } else {
            return NULL;
        }
    } else {
        /* Full permutation */
        sample_n = choice_count;
        have_n = true;
    }

    (void)have_n;

    size_t* selected = fisher_yates_sample(choice_count, sample_n);
    if (!selected) return NULL;

    Expr** result_elems = malloc(sizeof(Expr*) * sample_n);
    if (!result_elems) { free(selected); return NULL; }
    for (size_t i = 0; i < sample_n; i++) {
        result_elems[i] = expr_copy(choices[selected[i]]);
    }
    free(selected);

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), result_elems, sample_n);
    free(result_elems);
    return list;
}

void random_init(void) {
    symtab_add_builtin("RandomInteger", builtin_randominteger);
    symtab_get_def("RandomInteger")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomReal", builtin_randomreal);
    symtab_get_def("RandomReal")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("SeedRandom", builtin_seedrandom);
    symtab_get_def("SeedRandom")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomComplex", builtin_randomcomplex);
    symtab_get_def("RandomComplex")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomChoice", builtin_randomchoice);
    symtab_get_def("RandomChoice")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomSample", builtin_randomsample);
    symtab_get_def("RandomSample")->attributes |= ATTR_PROTECTED;
}
