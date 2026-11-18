/*
 * fi_internal.h — shared internals of the FactorInteger subsystem.
 *
 * Layout (one concern per file, all under src/numbertheory/factorinteger/):
 *
 *   fi_core.c      FiFactors list, FiCtx, deterministic RNG, the recursive
 *                  driver fi_factor() and the Automatic tier table.
 *   fi_builtin.c   FactorInteger[...] and every FactorInteger`Method[...]
 *                  entry point: argument/option parsing, messages, Gaussian
 *                  integers, output construction, registration + docstrings.
 *   fi_u64.c       word-size engine: Montgomery arithmetic, deterministic
 *                  primality, rho / SQUFOF / Lehman / micro-ECM on uint64_t.
 *   fi_u128.c      two-word Montgomery arithmetic and Pollard rho for
 *                  2^64 <= n < 2^127.
 *   fi_trial.c     batched trial division (one bignum pass per word-sized
 *                  product of primes) and the shared small-prime table.
 *   fi_rho.c       Pollard rho (Brent) on bignums.
 *   fi_fermat.c    Fermat (residue-sieved) and Lehman.
 *   fi_squfof.c    Shanks' square forms, racing multipliers.
 *   fi_ecm.c       ECM, Pollard p-1, Williams p+1 (GMP-ECM, with native
 *                  fallbacks when the build has no GMP-ECM).
 *   fi_qs.c        self-initialising quadratic sieve (FLINT qsieve or native).
 *   fi_cfrac.c     continued-fraction factoring (Morrison–Brillhart).
 *   fi_dixon.c     Dixon's random-squares method.
 *   fi_smooth.c    threaded trial division of CFRAC / Dixon candidates.
 *   fi_gf2.c       relation store, relation logs, GF(2) nullspace solver
 *                  shared by CFRAC / Dixon / QS.
 *   fi_par.c       thread regions and work claiming.
 *   fi_rbd.c       Blake rational base descent: options, variants, the
 *                  threaded base search and the hit cache (Mathilda glue).
 *   fi_rbd.h       its reentrant GMP/FLINT core, implemented in
 *                  fi_rbd_descent.c (iterated / exact descents, batched
 *                  gcds), fi_rbd_bases.c (Stern-Brocot and sum-ordered
 *                  bases), fi_rbd_lattice.c (Howgrave-Graham detector) and
 *                  fi_rbd_lll.c (FLINT fmpz_lll or a native integral LLL).
 *
 * A *splitter* takes an odd composite n that is not a perfect power and
 * either writes a non-trivial divisor 1 < d < n and returns 1, or returns 0
 * when its budget is exhausted.  Splitters never loop unboundedly: every one
 * carries an iteration/curve/relation budget, overridable per call through
 * the method suboptions in FiCtx.
 */
#ifndef MATHILDA_FI_INTERNAL_H
#define MATHILDA_FI_INTERNAL_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>
#include "expr.h"   /* mth_mpz_*_u64: GMP _ui calls are 32-bit on LLP64 */

/* ------------------------------------------------------------------ */
/* Methods                                                             */
/* ------------------------------------------------------------------ */

typedef enum {
    FI_AUTOMATIC = 0,
    FI_TRIAL,
    FI_RHO,
    FI_FERMAT,
    FI_LEHMAN,
    FI_SQUFOF,
    FI_CFRAC,
    FI_DIXON,
    FI_QS,
    FI_PM1,
    FI_PP1,
    FI_ECM,
    FI_RBD,
    FI_NMETHODS
} FiMethod;

/* ------------------------------------------------------------------ */
/* Factor list                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t   p;      /* prime (or, for an incomplete result, a composite) */
    int64_t e;      /* exponent; negative for denominator primes */
} FiEntry;

typedef struct {
    FiEntry* v;
    size_t   n, cap;
} FiFactors;

void fi_factors_init(FiFactors* f);
void fi_factors_clear(FiFactors* f);
/* Add p^e, merging with an existing entry for the same p. */
void fi_factors_add(FiFactors* f, const mpz_t p, int64_t e);
void fi_factors_add_ui(FiFactors* f, uint64_t p, int64_t e);
/* Sort ascending by p. */
void fi_factors_sort(FiFactors* f);

/* ------------------------------------------------------------------ */
/* Context                                                             */
/* ------------------------------------------------------------------ */

/* BlakeRationalBaseDescent variants ("Variant" suboption / alias heads). */
typedef enum {
    FI_RBD_EXACT = 0,     /* scaled exact descent, batched gcds (default)  */
    FI_RBD_CLASSIC,       /* the paper's R1-R4 reference descent           */
    FI_RBD_ALLHITS,       /* record every hit, refine into coprime pieces  */
    FI_RBD_LATTICE,       /* Howgrave-Graham detector in the depth band    */
    FI_RBD_TWOPASS,       /* exact pass, then lattice pass                 */
    FI_RBD_COPPERSMITH,   /* known a/b, c, n: one detector call            */
    FI_RBD_NVARIANTS
} FiRbdVariant;

/* "Order" of the enumerated bases (0, the default, is Weighted). */
enum {
    FI_RBD_ORDER_WEIGHTED = 0, /* ascending a / log2(a/b)                    */
    FI_RBD_ORDER_SUM,          /* a + b ascending ([RBD] Appendix A)         */
    FI_RBD_ORDER_SB,           /* Stern-Brocot breadth-first                 */
    FI_RBD_ORDER_HD            /* Stern-Brocot in doubling height bands      */
};

/* BlakeRationalBaseDescent suboptions; zero / unset means "default". */
typedef struct {
    int      variant;
    size_t   nbases;      /* "Base" -> a/b or {a1/b1, ...} (a > b, coprime) */
    mpz_t*   ba;
    mpz_t*   bb;
    int      order;       /* "Order": FI_RBD_ORDER_*                        */
    uint64_t height;      /* "Height": numerator bound of enumerated bases */
    uint64_t max_bases;   /* "Bases": enumeration budget                   */
    int      cmax_set, cmax_inf;
    mpz_t    cmax;        /* "MaxCoefficient" (Infinity: no deferral)      */
    int      mult_set;
    mpz_t    mult;        /* "Multiplier"                                  */
    double   beta;        /* "Beta": divisor size N^beta (lattice)         */
    int      tol_set;
    mpz_t    tol;         /* "Tolerance" X (lattice)                       */
    int      t, d;        /* "Multiplicity", "Dimension" (lattice)         */
    int      coef_set;
    mpz_t    coef;        /* "Coefficient" c (Coppersmith)                 */
    uint64_t expo;        /* "Exponent" n (Coppersmith)                    */
    int      range_set;
    mpq_t    lo, hi;      /* "BaseRange" (lattice pass grid)               */
    mpz_t*   cache;       /* coprime divisors found earlier in this call   */
    size_t   ncache;
    uint64_t ncalls;      /* splitter calls so far (0: the input itself)   */
} FiRbdOpts;

void fi_rbd_opts_init(FiRbdOpts* o);
void fi_rbd_opts_clear(FiRbdOpts* o);
/* Append a/b to the "Base" list (given as p/q with p, q > 0, p != q; it is
 * reduced, and inverted when below 1); reset_bases empties the list. */
void fi_rbd_opts_add_base(FiRbdOpts* o, const mpz_t p, const mpz_t q);
void fi_rbd_opts_reset_bases(FiRbdOpts* o);

typedef struct {
    FiMethod method;

    /* FactorInteger[n, k]: at most k entries in the answer (the last may be
     * an unfactored remainder).  0 = complete factorisation. */
    int64_t  k_limit;

    /* Method suboptions; 0 means "the method's default". */
    uint64_t bound;       /* TrialDivision "Bound"                       */
    uint64_t iterations;  /* Rho/Fermat/Lehman/SQUFOF/CFRAC/Dixon budget */
    uint64_t starts;      /* PollardRho "Starts"                         */
    uint64_t curves;      /* ECM "Curves"                                */
    double   B1, B2;      /* ECM / p-1 / p+1 stage bounds                */
    uint64_t sigma;       /* ECM "Sigma" (first curve; then sigma+1, ...)*/
    FiRbdOpts rbd;        /* BlakeRationalBaseDescent suboptions + cache */

    /* State. */
    int      incomplete;  /* some composite was returned unfactored      */
    uint64_t rng;         /* xorshift64* state, seeded from n            */
    int      warned;      /* set: suppress nofac (internal callers)      */
} FiCtx;

void     fi_ctx_init(FiCtx* c, FiMethod m);
void     fi_ctx_clear(FiCtx* c);
void     fi_ctx_seed(FiCtx* c, const mpz_t n);
uint64_t fi_rand(FiCtx* c);

/* Probable-prime test used throughout (BPSW; deterministic below 2^64). */
int fi_is_prime(const mpz_t n);

/* Factor n >= 1 completely (subject to c->k_limit) with c->method, adding
 * prime powers to out.  A composite the method could not split is added
 * with exponent equal to its multiplicity and c->incomplete is set. */
void fi_factor(const mpz_t n, FiFactors* out, FiCtx* c);

/* Complete factorisation of |n| (Automatic, silent), sorted; returns 0 if a
 * composite survived the search (the caller must then decline). */
int fi_factor_complete_mpz(const mpz_t n, FiFactors* out);

/* Named method on a machine word (allocation-free); -1 if the method has no
 * word splitter.  Fills sorted p[]/e[] (room for 64). */
int fi_word_method_array(uint64_t n, FiCtx* c, uint64_t* p, int64_t* e);

/* ------------------------------------------------------------------ */
/* Word-size engine (fi_u64.c)                                         */
/* ------------------------------------------------------------------ */

int      fi_u64_is_prime(uint64_t n);
/* Complete factorisation of n >= 1 (multiplicity mult). */
void     fi_u64_factor(uint64_t n, int64_t mult, FiFactors* out);
/* Allocation-free: sorted prime powers of n into p[]/e[] (room for 64). */
int      fi_u64_factor_array(uint64_t n, uint64_t* p, int64_t* e);
/* Splitters on odd composite n; return a non-trivial divisor or 0. */
uint64_t fi_u64_rho(uint64_t n, uint64_t c, uint64_t max_iters);
uint64_t fi_u64_squfof(uint64_t n, uint64_t max_iters);
uint64_t fi_u64_lehman(uint64_t n);
uint64_t fi_u64_fermat(uint64_t n, uint64_t max_iters);
uint64_t fi_u64_ecm(uint64_t n, uint64_t* rng, int max_curves);
uint64_t fi_u64_pm1(uint64_t n, uint64_t B1);
uint64_t fi_u64_pp1(uint64_t n, uint64_t B1, uint64_t seed);
uint64_t fi_u64_gcd(uint64_t a, uint64_t b);
uint64_t fi_u64_isqrt(uint64_t n);
/* Prime bitmap up to *limit (65536), built once. */
const unsigned char* fi_prime_bitmap(uint32_t* limit);
/* Lowest-level helpers shared with fi_u128.c. */
void     fi_mul64(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo);

/* ------------------------------------------------------------------ */
/* Two-word engine (fi_u128.c)                                         */
/* ------------------------------------------------------------------ */

/* Pollard–Brent rho on 2^64 <= n < 2^126 using two-limb Montgomery
 * arithmetic.  Returns 1 and sets d on success. */
int fi_u128_rho(mpz_t d, const mpz_t n, uint64_t c, uint64_t max_iters);
/* Two-limb micro-ECM ("tinyecm") for 2^64 <= n < 2^126: up to max_curves
 * Suyama curves at stage-1 bound B1 (B2 = 25 B1).  Returns 1 with d set. */
int fi_u128_ecm(mpz_t d, const mpz_t n, uint64_t* rng, uint32_t B1, int max_curves);

/* ------------------------------------------------------------------ */
/* Trial division (fi_trial.c)                                         */
/* ------------------------------------------------------------------ */

/* Primes below FI_SMALL_PRIME_LIMIT, generated once. */
#define FI_SMALL_PRIME_LIMIT (1u << 20)
const uint32_t* fi_small_primes(size_t* count);

/* Divide out every prime p <= bound from n (in place), adding p^e * mult to
 * out.  Returns the number of distinct primes found. */
size_t fi_trial_divide(mpz_t n, uint64_t bound, int64_t mult, FiFactors* out);

/* ------------------------------------------------------------------ */
/* Bignum splitters                                                    */
/* ------------------------------------------------------------------ */

int fi_split_rho(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters, uint64_t starts);
int fi_split_fermat(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_lehman(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_squfof(mpz_t d, const mpz_t n, FiCtx* c, uint64_t max_iters);
int fi_split_cfrac(mpz_t d, const mpz_t n, FiCtx* c);
int fi_split_dixon(mpz_t d, const mpz_t n, FiCtx* c);
int fi_split_rbd(mpz_t d, const mpz_t n, FiCtx* c);

int fi_split_qs(mpz_t d, const mpz_t n, FiCtx* c);
/* Knuth–Schroeppel multiplier for CFRAC / QS. */
unsigned fi_ks_multiplier(const mpz_t n);

/* ECM family.  B1/curves of 0 take the method's defaults. */
int fi_split_ecm(mpz_t d, const mpz_t n, FiCtx* c, double B1, uint64_t curves);
int fi_split_pm1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2);
int fi_split_pp1(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2);
int fi_split_pp1_seeds(mpz_t d, const mpz_t n, FiCtx* c, double B1, double B2, int nseeds);
/* Automatic's ECM tier: run curves at the B1 level that targets factors of
 * `digits` digits; returns 1 with d set on success. */
int fi_ecm_level(mpz_t d, const mpz_t n, FiCtx* c, int digits);
/* Expected curve count and B1 for a digit level (GMP-ECM table). */
void fi_ecm_level_params(int digits, double* B1, uint64_t* curves);
/* 1 when this build has GMP-ECM. */
int fi_have_gmp_ecm(void);

/* ------------------------------------------------------------------ */
/* GF(2) linear algebra (fi_gf2.c)                                     */
/* ------------------------------------------------------------------ */

/* Sparse relation matrix over GF(2): `nrows` relations, each a list of
 * column indices (factor-base primes with odd exponent).  Finds up to
 * `max_deps` linear dependencies among the rows; each dependency is a bit
 * vector over rows, returned in deps[k * words(nrows) ...].  Returns the
 * number of dependencies found. */
typedef struct {
    uint32_t  nrows, ncols;
    uint32_t* start;     /* row r's columns are col[start[r] .. start[r+1]) */
    uint32_t* col;
} FiSparseMat;

size_t fi_gf2_dependencies(const FiSparseMat* m, uint64_t* deps, size_t max_deps);
static inline size_t fi_gf2_words(size_t nrows) { return (nrows + 63) / 64; }

/* Relation store for the congruence-of-squares methods.  Column 0 is the
 * sign (-1); column c >= 1 is the prime fb[c].  A relation records
 * X^2 == prod fb[col]^exp * L (mod n), col sorted ascending; L = 1 for a
 * full relation, else one large prime (partials pair up automatically). */
typedef struct FiRels FiRels;
FiRels* fi_rels_new(const uint32_t* fb, size_t ncols);
void    fi_rels_free(FiRels* R);
int     fi_rels_add(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                    size_t cnt, uint64_t L, const mpz_t n);
/* Two large primes (each 1 or a prime outside the base). */
int     fi_rels_add2(FiRels* R, const mpz_t X, const uint32_t* col, const uint32_t* exps,
                     size_t cnt, uint64_t L1, uint64_t L2, const mpz_t n);
size_t  fi_rels_count(const FiRels* R);
size_t  fi_rels_partials(const FiRels* R);
/* Columns that occur with an odd exponent in some full relation: once the
 * relation count exceeds this, dependencies must exist. */
size_t  fi_rels_active(const FiRels* R);
/* Linear algebra + square roots; returns 1 with a proper divisor in d. */
int     fi_rels_solve(FiRels* R, const mpz_t n, mpz_t d);

/* Relation log: relations found by a worker thread, kept in the order found
 * and replayed by the main thread (the relation store is not thread-safe).
 * `tag` orders records across workers (CFRAC / Dixon: the candidate index);
 * the (col, exp) lists live in a shared pool. */
typedef struct {
    mpz_t    X;
    uint64_t tag, L1, L2;
    uint32_t off, cnt;
} FiLogRec;

typedef struct {
    FiLogRec* r;
    size_t    n, ninit, cap;   /* r[0..ninit) have X initialised (reused) */
    uint32_t* col;
    uint32_t* exps;
    size_t    npool, cappool;
} FiRelLog;

void fi_rellog_push(FiRelLog* lg, uint64_t tag, const mpz_t X, const uint32_t* col,
                    const uint32_t* exps, size_t cnt, uint64_t L1, uint64_t L2);
/* Add every record to R in log order, then empty the log. */
void fi_rellog_replay(FiRelLog* lg, FiRels* R, const mpz_t n);
/* Visit the records of logs[0..nlogs) merged by ascending tag (each log must
 * already be ascending); stops at the first nonzero visit() and returns it.
 * Empties every log either way. */
int  fi_rellog_merge(FiRelLog* logs, int nlogs,
                     int (*visit)(void* user, const FiLogRec* r, const uint32_t* col,
                                  const uint32_t* exps),
                     void* user);
void fi_rellog_free(FiRelLog* lg);

/* ------------------------------------------------------------------ */
/* Smooth-relation collection for CFRAC and Dixon (fi_smooth.c)        */
/* ------------------------------------------------------------------ */

/* Factor v (> 0) over fb[1..nfb] (fb[0] is the sign column), appending the
 * sorted (col, exp) list after the cnt entries already present; *L is the
 * large-prime cofactor (1 when smooth).  Returns 1 for a full relation or a
 * partial with L <= lp_max.  Early abort: after the first quarter of the base,
 * give up while the cofactor still has more than abort_bits bits. */
int fi_fb_factor(mpz_t v, const uint32_t* fb, size_t nfb, uint32_t* col, uint32_t* exps,
                 size_t* cnt, uint64_t* L, uint64_t lp_max, size_t abort_bits);

/* One region trial-divides candidates [0, count) on `nthreads` threads,
 * which claim chunks of `chunk` indices.  prep() writes candidate idx as
 * X^2 == v (mod n) (and may seed col/exps, e.g. the sign column); returning 0
 * skips it.  Hits go to per-thread logs tagged by index, and
 * fi_smooth_replay() hands them to the main thread in index order, so the
 * relation stream is the same whatever the thread count.  side(), if set,
 * runs on worker 0 before it starts claiming (CFRAC generates its next batch
 * there).  prep() and side() run off the main thread: they may only touch
 * their own state and the scratch they are given. */
typedef struct FiSmooth FiSmooth;
typedef int  (*FiSmoothPrep)(void* user, size_t idx, mpz_t X, mpz_t v, mpz_t tmp,
                             uint32_t* col, uint32_t* exps, size_t* cnt);
FiSmooth* fi_smooth_new(int nthreads, size_t max_nfb);
void      fi_smooth_run(FiSmooth* S, size_t count, size_t chunk, FiSmoothPrep prep,
                        void (*side)(void*), void* user, const uint32_t* fb, size_t nfb,
                        uint64_t lp_max, size_t abort_bits);
int       fi_smooth_replay(FiSmooth* S,
                           int (*visit)(void* user, const FiLogRec* r, const uint32_t* col,
                                        const uint32_t* exps),
                           void* user);
void      fi_smooth_free(FiSmooth* S);

/* ------------------------------------------------------------------ */
/* Messages                                                            */
/* ------------------------------------------------------------------ */

/* FactorInteger::nofac for composite m (unless c->warned). */
void fi_warn_incomplete(FiCtx* c, const mpz_t m);

/* ------------------------------------------------------------------ */
/* Parallel regions (fi_par.c)                                         */
/* ------------------------------------------------------------------ */

#define FI_MAXTHREADS 16
/* Worker threads available (1 without MATHILDA_THREADS). */
int  fi_nthreads(void);
/* Run fn on n jobs laid out stride bytes apart, one thread each, inside a
 * TimeConstrained-deferred region (see fi_par.c). */
void fi_run_parallel(void (*fn)(void*), void* jobs, size_t stride, int n);
/* Nonzero once an enclosing TimeConstrained deadline has passed.  Safe on a
 * worker thread (the main thread is parked in fi_run_parallel's joins and
 * writes no deadline state meanwhile): a long region polls it -- or hands it
 * to GMP-ECM as stop_asap -- and winds down, and the main thread then takes
 * the timeout as soon as the region ends. */
int  fi_stop_requested(void);

/* Work claiming inside a region: indices [0, end) handed out `chunk` at a
 * time under a lock (a plain cursor without MATHILDA_THREADS). */
typedef struct {
    size_t next, end;
    void*  lock;
} FiClaim;
void fi_claim_init(FiClaim* q);
void fi_claim_reset(FiClaim* q, size_t end);
/* Claim [*lo, *hi); returns 0 once the range is exhausted. */
int  fi_claim_next(FiClaim* q, size_t chunk, size_t* lo, size_t* hi);
void fi_claim_destroy(FiClaim* q);
/* Lower the end of the range to at most `end` (a search that has found a
 * hit at index `end` stops handing out later indices); and read it. */
void   fi_claim_cap(FiClaim* q, size_t end);
size_t fi_claim_end(FiClaim* q);

/* ------------------------------------------------------------------ */
/* Bit helpers (GCC builtins where available, portable otherwise)      */
/* ------------------------------------------------------------------ */

static inline int fi_ctz64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_ctzll(x);
#else
    int c = 0;
    while (!(x & 1)) { x >>= 1; c++; }
    return c;
#endif
}

static inline int fi_parity64(uint64_t x) {
#if defined(__GNUC__)
    return __builtin_parityll(x);
#else
    x ^= x >> 32; x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1;
    return (int)(x & 1);
#endif
}

/* ------------------------------------------------------------------ */
/* Tuning hooks                                                        */
/* ------------------------------------------------------------------ */

/* Build with -DFI_TUNING to let environment variables override the tuned
 * constants (FI_QS_FB, FI_QS_BLK, FI_QS_LP, FI_QS_SLACK, FI_QS_SMALL,
 * FI_QS_MED, FI_TINY_B1, FI_TINY_CURVES, FI_DENSE_MAX) for parameter
 * sweeps.  Compiled out otherwise. */
#ifdef FI_TUNING
#include <stdlib.h>
#define FI_TUNE_U(name, var) do { const char* s_ = getenv(name); if (s_) (var) = (unsigned)atol(s_); } while (0)
#define FI_TUNE_D(name, var) do { const char* s_ = getenv(name); if (s_) (var) = atof(s_); } while (0)
#else
#define FI_TUNE_U(name, var) ((void)0)
#define FI_TUNE_D(name, var) ((void)0)
#endif

#endif /* MATHILDA_FI_INTERNAL_H */
