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
 * RandomPrime[{imin, imax}]    - pseudorandom prime in [imin, imax], uniform over primes
 * RandomPrime[imax]            - pseudorandom prime in [2, imax]
 * RandomPrime[range, n]        - list of n pseudorandom primes
 * RandomPrime[range, {n1, ...}] - n1*n2*... array of primes
 *
 * SeedRandom[n]                - seeds the RNG with integer n
 * SeedRandom[]                 - reseeds from system entropy
 *
 * Uses GMP's random state (gmp_randstate_t) for bignum support.
 */

#include "random.h"
#include "dist.h"   /* ml_dist_reset_cache -- SeedRandom must clear the Gaussian cache */
#include "checked_int.h"   /* ci_add_i64 */
#include "arithmetic.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include "numeric.h"
#include "ndarray.h"
#include "pack.h"
#include "facint.h"   /* mth_is_prime_u64/_mpz -- RandomPrime shares PrimeQ's test */
#include "message.h"
#include "common.h"   /* builtin_arg_error */
#include "print.h"    /* expr_to_string, for intp/posdim message text */
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

/* All zero until xs_seed runs, which xs_next() guarantees for itself -- see the
 * comment there. Left uninitialised deliberately: the all-zero state is
 * xoshiro's fixed point and so is a LOUD failure (a hang), which is what a
 * missing seed should be. A plausible-looking default would turn the same
 * omission into a silently identical stream on every run. */
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

static void ensure_rand_init(void);   /* seeds BOTH generators; defined below */

/*
 * THE GENERATOR SEEDS ITSELF. This is not a convenience -- it closes a bug
 * class. Every per-element helper (random_integer_range, random_uniform_01)
 * called ensure_rand_init() for itself, so a new caller that drew from xoshiro
 * DIRECTLY inherited none of it: RandomInteger's packed producer did exactly
 * that, and in a fresh process the state was still all zero. All zero is
 * xoshiro's one fixed point, so xs_next() returned 0 forever --
 * RandomInteger[{1,1024}, 300] answered 300 copies of 1, and
 * RandomInteger[{1,10}, 300] never returned at all, spinning in xs_below()'s
 * rejection loop.
 *
 * The check is free at the rate this is called: RandomReal has always paid this
 * same branch per element (random_uniform_01 line ~600) and still draws at
 * 2.3 ns. A predictable branch on a hot global costs nothing next to the shifts
 * and xors below, and it is the only placement a future caller cannot bypass.
 */
static inline uint64_t xs_next(void) {
    ensure_rand_init();
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

/*
 * Batch draws: the SAME stream as n calls of xs_below(m) / random_uniform_01(),
 * element for element, so a seeded RandomChoice answers exactly as before --
 * but with the state in four locals instead of the globals. Through xs_next()
 * every draw paid an out-of-line ensure_rand_init() call and a reload and
 * store of g_xs (the caller's memcpy into the result buffer may alias them),
 * ~6 ns per element; here a draw is the xoshiro step and one 128-bit multiply,
 * ~1 ns. The state is written back once at the end.
 */
static void xs_below_fill(uint64_t m, uint64_t* out, size_t n) {
    ensure_rand_init();
    uint64_t s0 = g_xs[0], s1 = g_xs[1], s2 = g_xs[2], s3 = g_xs[3];
#define XS_STEP(r) do {                                                        \
        (r) = xs_rotl(s0 + s3, 23) + s0;                                       \
        uint64_t t_ = s1 << 17;                                                \
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t_;                      \
        s3 = xs_rotl(s3, 45);                                                  \
    } while (0)
    if (m <= 1) {
        for (size_t i = 0; i < n; i++) out[i] = 0;    /* xs_below(1) draws nothing */
    } else if ((m & (m - 1)) == 0) {
        for (size_t i = 0; i < n; i++) { uint64_t r; XS_STEP(r); out[i] = r & (m - 1); }
    } else {
#if defined(__SIZEOF_INT128__)
        uint64_t threshold = (uint64_t)(-(int64_t)m) % m;   /* == 2^64 mod m */
        for (size_t i = 0; i < n; i++) {
            uint64_t r; XS_STEP(r);
            __uint128_t prod = (__uint128_t)r * (__uint128_t)m;
            while ((uint64_t)prod < threshold) {            /* p < m/2^64: rare */
                XS_STEP(r);
                prod = (__uint128_t)r * (__uint128_t)m;
            }
            out[i] = (uint64_t)(prod >> 64);
        }
#else
        uint64_t threshold = (uint64_t)(-(int64_t)m) % m;
        for (size_t i = 0; i < n; i++) {
            uint64_t r;
            do { XS_STEP(r); } while (r < threshold);
            out[i] = r % m;
        }
#endif
    }
    g_xs[0] = s0; g_xs[1] = s1; g_xs[2] = s2; g_xs[3] = s3;
}

/* n uniforms on [0, 1): the random_uniform_01() stream, batched. */
static void xs_u01_fill(double* out, size_t n) {
    ensure_rand_init();
    uint64_t s0 = g_xs[0], s1 = g_xs[1], s2 = g_xs[2], s3 = g_xs[3];
    for (size_t i = 0; i < n; i++) {
        uint64_t r; XS_STEP(r);
        out[i] = (double)(r >> 11) * 0x1.0p-53;
    }
    g_xs[0] = s0; g_xs[1] = s1; g_xs[2] = s2; g_xs[3] = s3;
#undef XS_STEP
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
 *
 * EXPORTED (see random.h) because src/ml needs deviates that draw from THIS stream:
 * a distribution sampler that used its own generator would ignore SeedRandom, so
 * RandomVariate would not be reproducible while RandomReal was. Note that
 * random_internal_int_range is deliberately NOT the path for that -- it is reserved
 * for internal decision procedures that must not perturb the user-visible stream,
 * which is the opposite requirement.
 */
double random_uniform_01(void) {
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
    /* Reseeding must also drop any cached Gaussian deviate held by src/ml's
     * Box-Muller sampler. Without this, SeedRandom leaves a spare deviate in hand from
     * the PREVIOUS stream, so the first RandomVariate after reseeding differs between
     * two otherwise identical runs -- reproducibility would hold for RandomReal and
     * silently fail for RandomVariate, which is the worst kind of half-working. */
    ml_dist_reset_cache();

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

        /* A genuinely-complex rectangle (imaginary range not identically zero)
         * fills a packed complex64 buffer directly, just as RandomReal fills a
         * float64 one: the values are machine complex by construction, so there
         * is nothing to decide per element, and the common HermiteH[n, cArray]
         * / Sin[cArray] workflow then runs on the buffer instead of threading
         * boxed. The re-then-im draw order matches random_complex_range, so a
         * seeded stream is unchanged. A degenerate zero-imaginary range keeps the
         * old Real-valued List path (make_complex folds im=0 back to a Real). */
        bool complex_valued = (im_min != 0.0 || im_max != 0.0);

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            if (complex_valued) {
                void* vbuf = NULL;
                int64_t n64 = (int64_t)n;
                Expr* packed = ndbuild_open(1, &n64, NDT_COMPLEX64, &vbuf);
                if (packed) {
                    for (size_t i = 0; i < n; i++) {
                        double u1 = random_uniform_01();
                        double u2 = random_uniform_01();
                        ndt_set(vbuf, i, NDT_COMPLEX64,
                                re_min + u1 * (re_max - re_min),
                                im_min + u2 * (im_max - im_min));
                    }
                    return packed;
                }
            }
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
            size_t rank = dim_arg->data.function.arg_count;
            for (size_t i = 0; i < rank; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    return NULL;
                }
            }
            if (complex_valued && rank >= 1 && rank < NDARRAY_MAX_RANK) {
                int64_t dims[NDARRAY_MAX_RANK];
                int64_t total = 1;
                for (size_t i = 0; i < rank; i++) {
                    dims[i] = dim_arg->data.function.args[i]->data.integer;
                    total *= dims[i];
                }
                void* vbuf = NULL;
                Expr* packed = ndbuild_open((int)rank, dims, NDT_COMPLEX64, &vbuf);
                if (packed) {
                    /* Row-major fill order matches random_complex_array's leaf
                     * visitation, so a seeded stream gives the same array. */
                    for (int64_t i = 0; i < total; i++) {
                        double u1 = random_uniform_01();
                        double u2 = random_uniform_01();
                        ndt_set(vbuf, (size_t)i, NDT_COMPLEX64,
                                re_min + u1 * (re_max - re_min),
                                im_min + u2 * (im_max - im_min));
                    }
                    return packed;
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
static size_t* fisher_yates_sample(size_t total, size_t n);   /* defined below */
static size_t  random_index(size_t n);                        /* defined below */
static size_t  cum_search(const double* cum_weights, size_t count, double u);
static bool    is_upto(Expr* e, int64_t* val);                /* defined below */

/* ---------------------------------------------------------------------------
 *  NDArray buffer paths for RandomSample / RandomChoice.
 *
 *  Both heads had no fast path on either surface: a visible NDArray came back
 *  unevaluated, and a packed List was materialised into one Expr per element so
 *  that indices into it could be gathered.
 *
 *  THE DRAW SEQUENCE IS THE SAME ONE. Each of these calls exactly the routine
 *  the List path calls -- fisher_yates_sample for the sample, random_index per
 *  element for the choice -- and differs only in gathering int64/double words
 *  instead of copying Exprs. That is deliberate and it is the whole design
 *  constraint here: a producer that consumed the generator differently would
 *  make SeedRandom[1]; RandomSample[v] give one permutation for a packed v and
 *  another for the identical plain v, which is a surface disagreement in the
 *  one place a user cannot see it coming.
 *
 *  Returns NULL for anything outside the fast domain; the caller degrades.
 * ------------------------------------------------------------------------- */

/* Gather `n` elements of `a` at `idx` into a fresh buffer shaped like `a`. */
static Expr* nd_gather_indices(const Expr* a, const size_t* idx, size_t n) {
    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    const char* src = (const char*)a->data.ndarray.data;
    char* out = malloc(esz * (n ? n : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < n; i++)
        memcpy(out + i * esz, src + idx[i] * esz, esz);
    int64_t odims[1];
    odims[0] = (int64_t)n;
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

/* RandomSample[a] / RandomSample[a, n] over a rank-1 buffer. */
static Expr* nd_random_sample(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* a = res->data.function.args[0];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    size_t total = (size_t)a->data.ndarray.dims[0];
    if (total == 0) return NULL;

    size_t want = total;
    if (argc == 2) {
        Expr* n_arg = res->data.function.args[1];
        int64_t upto_val;
        if (n_arg->type == EXPR_INTEGER && n_arg->data.integer >= 0) {
            if ((uint64_t)n_arg->data.integer > (uint64_t)total) return NULL;
            want = (size_t)n_arg->data.integer;
        } else if (is_upto(n_arg, &upto_val) && upto_val >= 0) {
            want = ((uint64_t)upto_val > (uint64_t)total) ? total : (size_t)upto_val;
        } else {
            return NULL;
        }
    }
    size_t* sel = fisher_yates_sample(total, want);
    if (!sel) return NULL;
    Expr* out = nd_gather_indices(a, sel, want);
    free(sel);
    return out;
}

/* RandomChoice[a, n] over a rank-1 buffer -- n independent draws WITH
 * replacement. The bare RandomChoice[a] returns one element, which is a scalar
 * and has no buffer to build, so it is left to the List path. */
static Expr* nd_random_choice(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* n_arg = res->data.function.args[1];
    if (!is_ndarray(a) || a->data.ndarray.rank != 1) return NULL;
    if (n_arg->type != EXPR_INTEGER || n_arg->data.integer < 0) return NULL;
    size_t total = (size_t)a->data.ndarray.dims[0];
    if (total == 0) return NULL;
    size_t n = (size_t)n_arg->data.integer;

    NDType dt = a->data.ndarray.dtype;
    size_t esz = ndt_elem_size(dt);
    const char* src = (const char*)a->data.ndarray.data;
    char* out = malloc(esz * (n ? n : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < n; i++)
        memcpy(out + i * esz, src + random_index(total) * esz, esz);
    int64_t odims[1];
    odims[0] = (int64_t)n;
    return expr_new_ndarray_like(a, 1, odims, out, dt);
}

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
    return cum_search(cum_weights, count, random_uniform_01() * total);
}

/* The first index whose cumulative weight exceeds u (clamped to count - 1). */
static size_t cum_search(const double* cum_weights, size_t count, double u) {
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
 * Direct packed build for RandomChoice over a plain List of machine numbers.
 *
 * The choice list itself is almost always tiny ({1., 2.}, {-1, 1}) -- far under
 * the packing threshold, so it arrives as a plain List and never reaches
 * nd_random_choice. The RESULT, though, is as long as the user asks, and
 * building it as n boxed copies left every consumer off the buffer path:
 * Log[RandomChoice[{1., 2.}, 10^6]] cost 320 ms against 7 ms for the same
 * values packed. When every choice is a machine Real (or every one a machine
 * Integer) the answer is a uniform machine array, so it is written straight
 * into a buffer of `dims` (row-major).
 *
 * The draws are the SAME sequence the boxed builders make -- one
 * random_index / weighted_random_index call per element, in row-major order,
 * which is the order random_choice_array's depth-first recursion visits the
 * leaves -- so a seeded run gives the same values whichever path it takes.
 * Returns NULL, having drawn nothing, when the choices are mixed or non-machine,
 * the shape is degenerate, or the result is under the packing threshold
 * (ndbuild_open declines); the caller then builds the List as before.
 * `cum_weights` NULL means uniform.
 */
/* The shape of a rectangular nest of Lists whose leaves are all machine
 * numbers of ONE type (*t, EXPR_REAL or EXPR_INTEGER; set from the first leaf
 * when *tset is false).  A non-List is a rank-0 leaf.  False on ragged, mixed, or
 * non-machine leaves, or a nest deeper than max_rank. */
static bool choice_cell_shape(const Expr* e, int64_t* shp, int* r, int max_rank,
                              ExprType* t, bool* tset) {
    if (e->type == EXPR_REAL || e->type == EXPR_INTEGER) {
        if (!*tset) { *t = e->type; *tset = true; }
        *r = 0;
        return e->type == *t;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List ||
        e->data.function.arg_count == 0 || max_rank < 1)
        return false;
    int64_t sub[NDARRAY_MAX_RANK]; int sr = 0;
    if (!choice_cell_shape(e->data.function.args[0], sub, &sr, max_rank - 1, t, tset)) return false;
    for (size_t i = 1; i < e->data.function.arg_count; i++) {
        int64_t s2[NDARRAY_MAX_RANK]; int r2 = 0;
        if (!choice_cell_shape(e->data.function.args[i], s2, &r2, max_rank - 1, t, tset) ||
            r2 != sr || memcmp(s2, sub, sizeof(int64_t) * (size_t)sr) != 0)
            return false;
    }
    shp[0] = (int64_t)e->data.function.arg_count;
    memcpy(shp + 1, sub, sizeof(int64_t) * (size_t)sr);
    *r = sr + 1;
    return true;
}

/* Row-major leaves of a nest that choice_cell_shape accepted. */
static void choice_cell_flatten(const Expr* e, double* d, int64_t* z, size_t* at) {
    if (e->type == EXPR_REAL)    { d[(*at)++] = e->data.real; return; }
    if (e->type == EXPR_INTEGER) { z[(*at)++] = e->data.integer; return; }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        choice_cell_flatten(e->data.function.args[i], d, z, at);
}

static Expr* choice_buffer(Expr** choices, size_t choice_count,
                           const int64_t* dims, int rank,
                           double* cum_weights, double total) {
    if (rank < 1 || choice_count == 0) return NULL;
    /* Every choice a machine number, or every one a rectangular array of them
     * of one shape (RandomChoice[w -> {{3, 1}, {3, -1}}, n] is an n x 2
     * array): one element type, one cell shape. */
    ExprType t0 = EXPR_INTEGER; bool tset = false;
    int64_t cshape[NDARRAY_MAX_RANK]; int crank = 0;
    if (!choice_cell_shape(choices[0], cshape, &crank, NDARRAY_MAX_RANK - 1 - rank, &t0, &tset))
        return NULL;
    for (size_t i = 1; i < choice_count; i++) {
        int64_t s2[NDARRAY_MAX_RANK]; int r2 = 0;
        if (!choice_cell_shape(choices[i], s2, &r2, NDARRAY_MAX_RANK - 1 - rank, &t0, &tset) ||
            r2 != crank || memcmp(s2, cshape, sizeof(int64_t) * (size_t)crank) != 0)
            return NULL;                                /* lossless: one head, one shape */
    }
    size_t cell = 1;
    for (int d = 0; d < crank; d++) cell *= (size_t)cshape[d];
    size_t n = 1;
    for (int d = 0; d < rank; d++) {
        if (dims[d] < 1) return NULL;
        if ((uint64_t)dims[d] > (uint64_t)(SIZE_MAX / 8) / n) return NULL;
        n *= (size_t)dims[d];
    }
    if (cell > 1 && n > (SIZE_MAX / 8) / cell) return NULL;
    int64_t odims[NDARRAY_MAX_RANK];
    memcpy(odims, dims, sizeof(int64_t) * (size_t)rank);
    memcpy(odims + rank, cshape, sizeof(int64_t) * (size_t)crank);
    NDType dt = (t0 == EXPR_REAL) ? NDT_FLOAT64 : NDT_INT64;
    void* buf = NULL;
    Expr* out = ndbuild_open(rank + crank, odims, dt, &buf);
    if (!out) return NULL;
    /* The choices' leaves, flattened once: cell values per choice. */
    double* cd = NULL; int64_t* cz = NULL;
    if (dt == NDT_FLOAT64) cd = malloc(sizeof(double) * choice_count * cell);
    else                   cz = malloc(sizeof(int64_t) * choice_count * cell);
    if (!cd && !cz) { expr_free(out); return NULL; }
    for (size_t i = 0; i < choice_count; i++) {
        size_t at = i * cell;
        choice_cell_flatten(choices[i], cd, cz, &at);
    }
    /* Draw the indices a chunk at a time (the stream of n random_index /
     * weighted_random_index calls, in order), then gather. A scalar cell is a
     * plain assignment; a row cell a memcpy. */
    enum { CHUNK = 1024 };
    uint64_t idx[CHUNK];
    double   uu[CHUNK];
    for (size_t b = 0; b < n; b += CHUNK) {
        size_t bn = (n - b < CHUNK) ? n - b : CHUNK;
        if (cum_weights) {
            xs_u01_fill(uu, bn);
            for (size_t j = 0; j < bn; j++)
                idx[j] = cum_search(cum_weights, choice_count, uu[j] * total);
        } else {
            xs_below_fill((uint64_t)choice_count, idx, bn);
        }
        if (cell == 1 && dt == NDT_FLOAT64) {
            double* o = (double*)buf + b;
            for (size_t j = 0; j < bn; j++) o[j] = cd[idx[j]];
        } else if (cell == 1) {
            int64_t* o = (int64_t*)buf + b;
            for (size_t j = 0; j < bn; j++) o[j] = cz[idx[j]];
        } else {
            for (size_t j = 0; j < bn; j++) {
                size_t i = b + j, k = (size_t)idx[j];
                if (dt == NDT_FLOAT64)
                    memcpy((double*)buf + i * cell, cd + k * cell, sizeof(double) * cell);
                else
                    memcpy((int64_t*)buf + i * cell, cz + k * cell, sizeof(int64_t) * cell);
            }
        }
    }
    free(cd); free(cz);
    return out;
}

/* Read a RandomChoice shape argument (n or {n1, ..., nk}, all machine
 * non-negative integers) into dims. False for any other shape. */
static bool choice_dims(const Expr* dim_arg, int64_t* dims, int* rank) {
    if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
        dims[0] = dim_arg->data.integer;
        *rank = 1;
        return true;
    }
    if (dim_arg->type != EXPR_FUNCTION ||
        dim_arg->data.function.head->type != EXPR_SYMBOL ||
        dim_arg->data.function.head->data.symbol.name != SYM_List) return false;
    size_t k = dim_arg->data.function.arg_count;
    if (k < 1 || k >= NDARRAY_MAX_RANK) return false;
    for (size_t i = 0; i < k; i++) {
        const Expr* d = dim_arg->data.function.args[i];
        if (d->type != EXPR_INTEGER || d->data.integer < 0) return false;
        dims[i] = d->data.integer;
    }
    *rank = (int)k;
    return true;
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
    if (is_ndarray(first_arg)) {
        Expr* nd = nd_random_choice(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
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

        {   /* machine-number choices: write the buffer directly */
            int64_t cdims[NDARRAY_MAX_RANK];
            int crank = 0;
            if (choice_dims(dim_arg, cdims, &crank)) {
                Expr* packed = choice_buffer(elems_args, elems_count, cdims, crank,
                                             cum_weights, total);
                if (packed) { free(cum_weights); return packed; }
            }
        }

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

    {   /* machine-number choices: write the buffer directly */
        int64_t cdims[NDARRAY_MAX_RANK];
        int crank = 0;
        if (choice_dims(dim_arg, cdims, &crank)) {
            Expr* packed = choice_buffer(choices, choice_count, cdims, crank, NULL, 0.0);
            if (packed) return packed;
        }
    }

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

/* Slot map for fisher_yates_sparse: open addressing on the slot index; an
 * absent key holds its own index (the implicit identity array). */
typedef struct { size_t* key; size_t* val; char* used; size_t mask; } FySlots;

static size_t* fy_slot(FySlots* m, size_t k) {
    size_t h = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 17) & m->mask;
    while (m->used[h] && m->key[h] != k) h = (h + 1) & m->mask;
    if (!m->used[h]) { m->used[h] = 1; m->key[h] = k; m->val[h] = k; }
    return &m->val[h];
}

static size_t* fisher_yates_sparse(size_t total, size_t n) {
    size_t cap = 16;
    while (cap < 4 * (n + 1)) cap <<= 1;      /* <= 2n keys: load <= 1/2 */
    FySlots m;
    m.key = malloc(cap * sizeof(size_t));
    m.val = malloc(cap * sizeof(size_t));
    m.used = calloc(cap, 1);
    m.mask = cap - 1;
    size_t* result = malloc(sizeof(size_t) * (n ? n : 1));
    if (!m.key || !m.val || !m.used || !result) {
        free(m.key); free(m.val); free(m.used); free(result);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) {
        size_t j = i + random_index(total - i);
        size_t* a = fy_slot(&m, i);
        size_t* b = fy_slot(&m, j);           /* no growth, so a stays valid */
        size_t tmp = *a; *a = *b; *b = tmp;
        result[i] = *a;
    }
    free(m.key); free(m.val); free(m.used);
    return result;
}

/*
 * Fisher-Yates partial shuffle: select n elements from indices [0, total)
 * without replacement. Returns an array of n selected indices.
 * Caller must free the returned array.
 */
static size_t* fisher_yates_sample(size_t total, size_t n) {
    /* A sparse sample (total much larger than n) runs the same partial shuffle
     * over a hash map of the displaced slots instead of an identity array of
     * length total: identical random_index draws, identical result, O(n)
     * memory. */
    if (total / 4 > n) return fisher_yates_sparse(total, n);
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

size_t* random_sample_indices(size_t total, size_t n) {
    if (n > total) return NULL;
    return fisher_yates_sample(total, n);
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
    if (is_ndarray(first_arg)) {
        Expr* nd = nd_random_sample(res);
        return nd ? nd : ndarray_delist_and_reeval(res);
    }
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

/* --------------------------------------------------------------------------
 * RandomPrime
 * --------------------------------------------------------------------------
 *
 * UNIFORM OVER PRIMES, BY REJECTION. A draw picks uniformly from the candidate
 * set C = {2 if 2 is in range} U {odd k in range} and rejects composites. Every
 * prime in [lo, hi] occurs in C exactly once, so the accepted value is exactly
 * uniform over the primes in the interval -- the Mathematica contract -- and
 * skipping the even numbers halves the expected number of draws to ~ln(N)/2.
 *
 * Primality is PrimeQ's own test (mth_is_prime_u64 / mth_is_prime_mpz in
 * facint.c), so RandomPrime can never return a number PrimeQ rejects.
 *
 * Three sampler shapes, chosen once per call and reused for every element:
 *
 *   RP_TABLE  span <= RP_TABLE_SPAN, machine-width. The primes are enumerated
 *             once and each draw is a single xs_below(count). This is also the
 *             exact noprime test for the only spans that can be prime-free in
 *             machine range, and it cannot loop.
 *   RP_U64    wider machine-width span: xoshiro rejection. No existence scan is
 *             needed -- the largest prime gap below 2^64 is 1550 (after
 *             18361375334787046697), well under RP_TABLE_SPAN, so a span this
 *             wide always contains a prime.
 *   RP_MPZ    a bound past int64: GMP rejection, after an upward scan from lo
 *             that finds one prime (as a witness only -- returning it would bias
 *             the draw) or proves the interval empty. Bignum gaps can exceed any
 *             fixed span, so here existence must be checked.
 * -------------------------------------------------------------------------- */

#define RP_TABLE_SPAN 4096

typedef enum { RP_TABLE, RP_U64, RP_MPZ } RpKind;

typedef struct {
    RpKind    kind;
    /* RP_TABLE */
    uint64_t* table;
    uint64_t  tcount;
    /* RP_U64: candidates are {2}? then odd_lo, odd_lo+2, ..., n_cand total */
    int       has_two;
    uint64_t  odd_lo;
    uint64_t  n_cand;
    /* RP_MPZ: same layout in mpz; scratch hoisted out of the draw loop */
    mpz_t     z_odd_lo, z_n_cand, z_draw;
    bool      z_ready;
} RpSampler;

static void rp_sampler_clear(RpSampler* S) {
    free(S->table);
    S->table = NULL;
    if (S->z_ready) {
        mpz_clear(S->z_odd_lo);
        mpz_clear(S->z_n_cand);
        mpz_clear(S->z_draw);
        S->z_ready = false;
    }
}

/* Build the sampler for [lo, hi] (lo >= 1 and lo <= hi, as parsed). Returns
 * false iff the interval holds no prime; S then needs no clear. */
static bool rp_sampler_init(RpSampler* S, const mpz_t lo_in, const mpz_t hi) {
    memset(S, 0, sizeof *S);
    ensure_rand_init();

    mpz_t lo;
    mpz_init_set(lo, lo_in);
    if (mpz_cmp_ui(lo, 2) < 0) mpz_set_ui(lo, 2);        /* 1 is not prime */
    if (mpz_cmp(hi, lo) < 0) { mpz_clear(lo); return false; }

    if (mpz_sizeinbase(hi, 2) <= 63) {
        /* Machine-width: every candidate fits int64, and so packs. Read out
         * with mpz_export, not mpz_get_si: long is 32 bits on LLP64. */
        uint64_t l = 0, h = 0;
        mpz_export(&l, NULL, -1, sizeof l, 0, 0, lo);
        mpz_export(&h, NULL, -1, sizeof h, 0, 0, hi);
        mpz_clear(lo);

        S->has_two = (l == 2);
        uint64_t ol = (l <= 3) ? 3 : (l | 1u);
        uint64_t oh = (h & 1u) ? h : h - 1;
        uint64_t n_odd = (oh >= ol) ? (oh - ol) / 2 + 1 : 0;
        S->odd_lo = ol;
        S->n_cand = n_odd + (uint64_t)S->has_two;

        if (h - l + 1 <= RP_TABLE_SPAN) {
            S->kind = RP_TABLE;
            S->table = malloc(sizeof(uint64_t) * (size_t)(S->n_cand ? S->n_cand : 1));
            if (!S->table) return false;
            if (S->has_two) S->table[S->tcount++] = 2;
            for (uint64_t i = 0; i < n_odd; i++) {
                uint64_t k = ol + 2 * i;
                if (mth_is_prime_u64(k)) S->table[S->tcount++] = k;
            }
            if (S->tcount == 0) { rp_sampler_clear(S); return false; }
            return true;
        }
        S->kind = RP_U64;
        return true;
    }

    /* Bignum bound. */
    S->kind = RP_MPZ;
    mpz_init(S->z_odd_lo);
    mpz_init(S->z_n_cand);
    mpz_init(S->z_draw);
    S->z_ready = true;

    S->has_two = (mpz_cmp_ui(lo, 2) == 0);
    if (mpz_cmp_ui(lo, 3) <= 0) mpz_set_ui(S->z_odd_lo, 3);
    else { mpz_set(S->z_odd_lo, lo); if (mpz_even_p(lo)) mpz_add_ui(S->z_odd_lo, lo, 1); }
    /* n_odd = (oh - ol)/2 + 1 with oh the largest odd <= hi; hi > 2^63 here,
     * so oh >= ol always. */
    mpz_t oh;
    mpz_init_set(oh, hi);
    if (mpz_even_p(oh)) mpz_sub_ui(oh, oh, 1);
    if (mpz_cmp(oh, S->z_odd_lo) >= 0) {
        mpz_sub(S->z_n_cand, oh, S->z_odd_lo);
        mpz_fdiv_q_2exp(S->z_n_cand, S->z_n_cand, 1);
        mpz_add_ui(S->z_n_cand, S->z_n_cand, 1);
    } else {
        mpz_set_ui(S->z_n_cand, 0);
    }
    mpz_add_ui(S->z_n_cand, S->z_n_cand, (unsigned long)S->has_two);

    /* Existence: 2 is a witness on its own; otherwise scan odd k upward. */
    bool found = S->has_two;
    if (!found) {
        mpz_set(S->z_draw, S->z_odd_lo);
        while (mpz_cmp(S->z_draw, oh) <= 0) {
            if (mth_is_prime_mpz(S->z_draw)) { found = true; break; }
            mpz_add_ui(S->z_draw, S->z_draw, 2);
        }
    }
    mpz_clear(oh);
    mpz_clear(lo);
    if (!found) { rp_sampler_clear(S); return false; }
    return true;
}

/* One prime from a machine-width sampler (RP_TABLE / RP_U64). */
static int64_t rp_draw_i64(RpSampler* S) {
    if (S->kind == RP_TABLE) return (int64_t)S->table[xs_below(S->tcount)];
    for (;;) {
        uint64_t j = xs_below(S->n_cand);
        uint64_t k;
        if (S->has_two) { if (j == 0) return 2; k = S->odd_lo + 2 * (j - 1); }
        else k = S->odd_lo + 2 * j;
        if (mth_is_prime_u64(k)) return (int64_t)k;
    }
}

/* One prime as an Expr, for any sampler kind. */
static Expr* rp_draw_expr(RpSampler* S) {
    if (S->kind != RP_MPZ) return expr_new_integer(rp_draw_i64(S));
    for (;;) {
        mpz_urandomm(S->z_draw, g_rand_state, S->z_n_cand);
        if (S->has_two) {
            if (mpz_sgn(S->z_draw) == 0) return expr_new_integer(2);
            mpz_sub_ui(S->z_draw, S->z_draw, 1);
        }
        mpz_mul_2exp(S->z_draw, S->z_draw, 1);
        mpz_add(S->z_draw, S->z_draw, S->z_odd_lo);
        if (mth_is_prime_mpz(S->z_draw))
            return expr_bigint_normalize(expr_new_bigint_from_mpz(S->z_draw));
    }
}

/* RandomPrime::intp naming `e`. Always returns false, for one-line use. */
static bool rp_intp(Expr* e) {
    char* s = expr_to_string(e);
    mth_message("RandomPrime", "intp", "%s is not a positive integer.", s ? s : "");
    free(s);
    return false;
}

static bool rp_is_positive_integer(Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer >= 1) ||
           (e->type == EXPR_BIGINT && mpz_sgn(e->data.bigint) > 0);
}

/* Parse the range argument into lo <= hi (both >= 1), emitting intp/prmrng
 * on failure. A reversed pair is accepted, as in Mathematica. */
static bool rp_parse_range(Expr* arg, mpz_t lo, mpz_t hi) {
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List) {
        if (arg->data.function.arg_count != 2) {
            char* s = expr_to_string(arg);
            mth_message("RandomPrime", "prmrng",
                        "First argument %s is not a positive integer or a list of two positive integers.",
                        s ? s : "");
            free(s);
            return false;
        }
        Expr* a = arg->data.function.args[0];
        Expr* b = arg->data.function.args[1];
        if (!rp_is_positive_integer(a)) return rp_intp(a);
        if (!rp_is_positive_integer(b)) return rp_intp(b);
        expr_to_mpz(a, lo);
        expr_to_mpz(b, hi);
        if (mpz_cmp(lo, hi) > 0) mpz_swap(lo, hi);
        return true;
    }
    if (!rp_is_positive_integer(arg)) return rp_intp(arg);
    mpz_set_ui(lo, 1);
    expr_to_mpz(arg, hi);
    return true;
}

/* Parse the dimension spec: n >= 1, or a list of such (possibly empty, which
 * Mathematica answers with a one-element list). Fills dims[0..*rank) and the
 * element total; emits posdim on failure. */
static bool rp_parse_dims(Expr* d, int64_t* dims, int* rank, int64_t* total) {
    bool ok = true;
    if (d->type == EXPR_INTEGER) {
        ok = d->data.integer >= 1;
        if (ok) { dims[0] = d->data.integer; *rank = 1; *total = dims[0]; }
    } else if (d->type == EXPR_FUNCTION &&
               d->data.function.head->type == EXPR_SYMBOL &&
               d->data.function.head->data.symbol.name == SYM_List &&
               d->data.function.arg_count < NDARRAY_MAX_RANK) {
        size_t r = d->data.function.arg_count;
        int64_t t = 1;
        for (size_t i = 0; i < r && ok; i++) {
            Expr* di = d->data.function.args[i];
            ok = di->type == EXPR_INTEGER && di->data.integer >= 1 &&
                 !ci_mul_i64(t, di->data.integer, &t);
            if (ok) dims[i] = di->data.integer;
        }
        if (ok && r == 0) { dims[0] = 1; r = 1; }
        if (ok) { *rank = (int)r; *total = t; }
    } else {
        ok = false;
    }
    if (!ok) {
        char* s = expr_to_string(d);
        mth_message("RandomPrime", "posdim",
                    "The dimensions parameter %s is expected to be a positive integer or a list of positive integers.",
                    s ? s : "");
        free(s);
    }
    return ok;
}

/* Nest a flat row-major Expr* array into Lists per dims; consumes flat[]. */
static Expr* rp_nest(Expr** flat, const int64_t* dims, int rank, int level, size_t* pos) {
    size_t n = (size_t)dims[level];
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;
    for (size_t i = 0; i < n; i++)
        elems[i] = (level == rank - 1) ? flat[(*pos)++]
                                       : rp_nest(flat, dims, rank, level + 1, pos);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), elems, n);
    free(elems);
    return list;
}

/*
 * builtin_randomprime - implements RandomPrime[...]
 *
 *   RandomPrime[imax]                  -> prime in [2, imax]
 *   RandomPrime[{imin, imax}]          -> prime in [imin, imax]
 *   RandomPrime[range, n]              -> list of n primes
 *   RandomPrime[range, {n1, n2, ...}]  -> n1 x n2 x ... array of primes
 *
 * A machine-width range writes straight into a packed int64 buffer, like
 * RandomInteger's producer; a bignum range builds Lists and offers them.
 */
Expr* builtin_randomprime(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return builtin_arg_error("RandomPrime", argc, 1, 2);

    mpz_t lo, hi;
    mpz_init(lo);
    mpz_init(hi);
    Expr* out = NULL;
    int64_t dims[NDARRAY_MAX_RANK];
    int rank = 0;
    int64_t total = 1;
    RpSampler S;

    if (!rp_parse_range(res->data.function.args[0], lo, hi)) goto done;
    if (argc == 2 && !rp_parse_dims(res->data.function.args[1], dims, &rank, &total))
        goto done;
    if (!rp_sampler_init(&S, lo, hi)) {
        mth_message("RandomPrime", "noprime", "There are no primes in the specified interval.");
        goto done;
    }

    if (argc == 1) {
        out = rp_draw_expr(&S);
    } else {
        void* buf = NULL;
        if (S.kind != RP_MPZ &&
            (out = ndbuild_open(rank, dims, NDT_INT64, &buf)) != NULL) {
            int64_t* b = (int64_t*)buf;
            for (int64_t i = 0; i < total; i++) b[i] = rp_draw_i64(&S);
        } else {
            Expr** flat = malloc(sizeof(Expr*) * (size_t)total);
            if (flat) {
                for (int64_t i = 0; i < total; i++) flat[i] = rp_draw_expr(&S);
                size_t pos = 0;
                out = rp_nest(flat, dims, rank, 0, &pos);
                free(flat);
                if (out) out = pack_offer(out);
            }
        }
    }
    rp_sampler_clear(&S);

done:
    mpz_clear(lo);
    mpz_clear(hi);
    return out;
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

    symtab_add_builtin("RandomPrime", builtin_randomprime);
    symtab_get_def("RandomPrime")->attributes |= ATTR_PROTECTED;
}
