/*
 * fi_core.c — factor list, context, and the recursive driver.
 *
 * fi_factor() runs a work stack of (composite, multiplicity, stage) items.
 * Every method shares the same pipeline:
 *
 *   1. strip powers of two (and, for Automatic / TrialDivision, every prime
 *      up to the trial bound, once, at the top);
 *   2. for each work item: word-sized → the u64 engine (Automatic only; a
 *      named method keeps its own algorithm at every size), prime → output,
 *      perfect power → its primitive root with the multiplicity scaled;
 *   3. otherwise hand it to the method's splitter, and push both parts.
 *
 * A composite the splitter cannot crack is emitted with its multiplicity,
 * FactorInteger::nofac is issued once, and c->incomplete is set.
 *
 * Automatic's splitter is a tier table, resumed per item through `stage`:
 * cheap Pollard rho, one Pollard p-1, then ECM in rising digit levels, with a
 * switch to the quadratic sieve once the ECM effort spent reaches 2/7 of the
 * composite's digit count (yafu's crossover rule) and the composite is small
 * enough for QS to be the faster finisher.  A part produced by a split keeps
 * its parent's stage: the parent's search already ran over it.
 */
#include "fi_internal.h"
#include "term/term_interrupt.h"   /* split_auto stops on a Ctrl-C */

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Factor list                                                         */
/* ------------------------------------------------------------------ */

void fi_factors_init(FiFactors* f) { f->v = NULL; f->n = f->cap = 0; }

void fi_factors_clear(FiFactors* f) {
    for (size_t i = 0; i < f->n; i++) mpz_clear(f->v[i].p);
    free(f->v);
    f->v = NULL;
    f->n = f->cap = 0;
}

static FiEntry* grow(FiFactors* f) {
    if (f->n == f->cap) {
        f->cap = f->cap ? 2 * f->cap : 16;
        f->v = realloc(f->v, f->cap * sizeof(FiEntry));
    }
    return &f->v[f->n++];
}

void fi_factors_add(FiFactors* f, const mpz_t p, int64_t e) {
    for (size_t i = 0; i < f->n; i++)
        if (mpz_cmp(f->v[i].p, p) == 0) { f->v[i].e += e; return; }
    FiEntry* x = grow(f);
    mpz_init_set(x->p, p);
    x->e = e;
}

void fi_factors_add_ui(FiFactors* f, uint64_t p, int64_t e) {
    for (size_t i = 0; i < f->n; i++)
        if (mpz_cmp_ui(f->v[i].p, (unsigned long)p) == 0) { f->v[i].e += e; return; }
    FiEntry* x = grow(f);
    mpz_init_set_ui(x->p, (unsigned long)p);
    x->e = e;
}

static int cmp_entry(const void* a, const void* b) {
    return mpz_cmp(((const FiEntry*)a)->p, ((const FiEntry*)b)->p);
}
void fi_factors_sort(FiFactors* f) {
    if (f->n > 1) qsort(f->v, f->n, sizeof(FiEntry), cmp_entry);
}

/* ------------------------------------------------------------------ */
/* Context                                                             */
/* ------------------------------------------------------------------ */

void fi_ctx_init(FiCtx* c, FiMethod m) {
    memset(c, 0, sizeof *c);
    c->method = m;
    fi_rbd_opts_init(&c->rbd);
    c->rng = UINT64_C(0x9E3779B97F4A7C15);
}

void fi_ctx_clear(FiCtx* c) {
    fi_rbd_opts_clear(&c->rbd);
}

/* The randomised methods (rho starts, ECM sigmas, QS) draw from a generator
 * seeded by n itself, so a given input always takes the same path: the
 * answer, the time and any nofac depend on n alone, never on what the
 * process factored before (MATHILDA_DIVERGENCES F9). */
void fi_ctx_seed(FiCtx* c, const mpz_t n) {
    uint64_t h = UINT64_C(0x9E3779B97F4A7C15);
    size_t sz = mpz_size(n);
    for (size_t i = 0; i < sz; i++) {
        h ^= (uint64_t)mpz_getlimbn(n, (mp_size_t)i);
        h *= UINT64_C(0xff51afd7ed558ccd);
        h ^= h >> 33;
    }
    c->rng = h ? h : 1;
}

uint64_t fi_rand(FiCtx* c) {
    uint64_t x = c->rng;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    c->rng = x;
    return x * UINT64_C(0x2545F4914F6CDD1D);
}

int fi_is_prime(const mpz_t n) {
    if (mpz_sizeinbase(n, 2) <= 64) return fi_u64_is_prime(mpz_get_ui(n));
    return mpz_probab_prime_p(n, 25) != 0;
}

/* ------------------------------------------------------------------ */
/* Work stack                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t   m;
    int64_t mult;
    int     stage;      /* Automatic: next tier to run */
    int     qs_tried;   /* Automatic: QS already attempted on this item */
} Work;

typedef struct { Work* v; size_t n, cap; } WorkStack;

static void push(WorkStack* s, const mpz_t m, int64_t mult, int stage) {
    if (s->n == s->cap) {
        s->cap = s->cap ? 2 * s->cap : 16;
        s->v = realloc(s->v, s->cap * sizeof(Work));
    }
    Work* w = &s->v[s->n++];
    mpz_init_set(w->m, m);
    w->mult = mult;
    w->stage = stage;
    w->qs_tried = 0;
}

/* Reduce m to its primitive root: m = b^e with b not a perfect power.
 * Returns e (1 when m is not a perfect power).  Restarting the exponent scan
 * after every hit accumulates the full exponent (p^6 -> p via 2 then 3). */
static int64_t primitive_root(mpz_t b, const mpz_t m) {
    mpz_set(b, m);
    if (!mpz_perfect_power_p(b)) return 1;
    int64_t e = 1;
    mpz_t r;
    mpz_init(r);
    for (int reduced = 1; reduced; ) {
        reduced = 0;
        unsigned long emax = (unsigned long)mpz_sizeinbase(b, 2);
        for (unsigned long k = 2; k <= emax; k++) {
            if (mpz_root(r, b, k)) {
                mpz_set(b, r);
                e *= (int64_t)k;
                reduced = 1;
                break;
            }
        }
    }
    mpz_clear(r);
    return e;
}

/* ------------------------------------------------------------------ */
/* Named-method splitters                                              */
/* ------------------------------------------------------------------ */

#define DEF(v, d) ((v) ? (v) : (d))

static int split_named(mpz_t d, const mpz_t m, FiCtx* c) {
    int word = mpz_sizeinbase(m, 2) <= 64;
    uint64_t w = word ? mpz_get_ui(m) : 0, r = 0;
    switch (c->method) {
    case FI_RHO:
        if (word) {
            uint64_t starts = DEF(c->starts, 16), it = DEF(c->iterations, UINT64_C(1) << 26);
            for (uint64_t k = 0; k < starts && !r; k++) r = fi_u64_rho(w, 1 + k, it);
            break;
        }
        return fi_split_rho(d, m, c, c->iterations, c->starts);
    case FI_FERMAT:
        if (word && !(w >> 62)) { r = fi_u64_fermat(w, DEF(c->iterations, UINT64_C(1) << 24)); break; }
        return fi_split_fermat(d, m, c, c->iterations);
    case FI_LEHMAN:
        if (word && !(w >> 42)) { r = fi_u64_lehman(w); break; }
        return fi_split_lehman(d, m, c, c->iterations);
    case FI_SQUFOF:
        if (word && !(w >> 62)) { r = fi_u64_squfof(w, c->iterations); if (r) break; }
        return fi_split_squfof(d, m, c, c->iterations);
    case FI_CFRAC:  return fi_split_cfrac(d, m, c);
    case FI_DIXON:  return fi_split_dixon(d, m, c);
    case FI_QS:     return fi_split_qs(d, m, c);
    case FI_PM1:    return fi_split_pm1(d, m, c, c->B1, c->B2);
    case FI_PP1:    return fi_split_pp1(d, m, c, c->B1, c->B2);
    case FI_ECM:    return fi_split_ecm(d, m, c, c->B1, c->curves);
    case FI_RBD:    return fi_split_rbd(d, m, c);
    default:        return 0;
    }
    if (r > 1 && r < w) { mpz_set_ui(d, (unsigned long)r); return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Named methods on words                                              */
/* ------------------------------------------------------------------ */

/* A named method with a word-size splitter factors a word input entirely in
 * word arithmetic (no mpz, no heap): what is measured is the method. */
static int has_word_splitter(FiMethod m) {
    return m == FI_RHO || m == FI_FERMAT || m == FI_LEHMAN || m == FI_SQUFOF ||
           m == FI_PM1 || m == FI_PP1 || m == FI_ECM;
}

/* Split odd composite word n (not a perfect power) with c->method; 0 when
 * the method gives up or does not apply at this size. */
static uint64_t word_split(uint64_t n, FiCtx* c) {
    uint64_t d = 0;
    switch (c->method) {
    case FI_RHO: {
        uint64_t starts = DEF(c->starts, 16), it = DEF(c->iterations, UINT64_C(1) << 26);
        for (uint64_t k = 0; k < starts && !d; k++) d = fi_u64_rho(n, 1 + k, it);
        break;
    }
    case FI_FERMAT:
    case FI_LEHMAN: {
        int fits = c->method == FI_FERMAT ? !(n >> 62) : !(n >> 42);
        if (fits) {
            d = c->method == FI_FERMAT ? fi_u64_fermat(n, DEF(c->iterations, UINT64_C(1) << 24))
                                       : fi_u64_lehman(n);
        } else {
            /* Past the word variant's range: the budgeted bignum version. */
            mpz_t m, f;
            mpz_init_set_ui(m, (unsigned long)n);
            mpz_init(f);
            int ok = c->method == FI_FERMAT ? fi_split_fermat(f, m, c, c->iterations)
                                            : fi_split_lehman(f, m, c, c->iterations);
            if (ok) d = mpz_get_ui(f);
            mpz_clears(m, f, NULL);
        }
        break;
    }
    case FI_SQUFOF:
        if (!(n >> 62)) d = fi_u64_squfof(n, c->iterations);
        if (!d) {
            /* kN overflows a word for the larger multipliers: two-word forms. */
            mpz_t m, f;
            mpz_init_set_ui(m, (unsigned long)n);
            mpz_init(f);
            if (fi_split_squfof(f, m, c, c->iterations)) d = mpz_get_ui(f);
            mpz_clears(m, f, NULL);
        }
        break;
    case FI_PM1:    d = fi_u64_pm1(n, c->B1 > 0 ? (uint64_t)c->B1 : 0); break;
    case FI_PP1:
        for (int s = 0; s < 3 && !d; s++) d = fi_u64_pp1(n, c->B1 > 0 ? (uint64_t)c->B1 : 0, fi_rand(c));
        break;
    case FI_ECM:    d = fi_u64_ecm(n, &c->rng, c->curves ? (int)c->curves : 256); break;
    default: break;
    }
    return (d > 1 && d < n) ? d : 0;
}

/* x^e as a word, or 0 on overflow. */
static uint64_t wpow(uint64_t x, int e) {
    uint64_t r = 1;
    while (e) {
        uint64_t h, l;
        if (e & 1) { fi_mul64(r, x, &h, &l); if (h) return 0; r = l; }
        e >>= 1;
        if (e) { fi_mul64(x, x, &h, &l); if (h) return 0; x = l; }
    }
    return r;
}

/* n = b^e with b not a perfect power; returns e (1 if none) and sets *b.
 * Only prime exponents up to log2(n) are tried, each by a rounded floating
 * root checked exactly; a hit recurses on b, accumulating the exponent. */
static int word_root(uint64_t n, uint64_t* b) {
    static const int PE[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
    int bits = 0;
    for (uint64_t t = n; t; t >>= 1) bits++;
    for (size_t i = 0; i < sizeof PE / sizeof PE[0] && PE[i] < bits; i++) {
        int e = PE[i];
        uint64_t r = (uint64_t)(pow((double)n, 1.0 / e) + 0.5);
        for (uint64_t x = r > 2 ? r - 1 : 2; x <= r + 1; x++) {
            if (wpow(x, e) == n) {
                int e2 = word_root(x, b);
                return e * e2;
            }
        }
    }
    *b = n;
    return 1;
}

/* Sink for word results: prime powers, plus at most a few unfactored
 * composites (flagged). */
typedef struct { uint64_t p[64]; int64_t e[64]; int n; } WordOut;

static void wo_add(WordOut* w, uint64_t p, int64_t e) {
    for (int i = 0; i < w->n; i++) if (w->p[i] == p) { w->e[i] += e; return; }
    w->p[w->n] = p; w->e[w->n++] = e;
}

static void word_factor_method_w(uint64_t n, int64_t mult, WordOut* out, FiCtx* c) {
    if (n <= 1) return;
    if (fi_u64_is_prime(n)) { wo_add(out, n, mult); return; }
    uint64_t b;
    int e = word_root(n, &b);
    if (e > 1) { word_factor_method_w(b, mult * e, out, c); return; }
    if (c->method == FI_PM1 || c->method == FI_PP1 || c->method == FI_ECM) {
        size_t np;
        const uint32_t* P = fi_small_primes(&np);
        for (size_t i = 1; i < np && P[i] < 1024; i++)
            if (n % P[i] == 0 && n != P[i]) {
                word_factor_method_w(P[i], mult, out, c);
                word_factor_method_w(n / P[i], mult, out, c);
                return;
            }
    }
    uint64_t d = word_split(n, c);
    if (!d) {
        mpz_t m;
        mpz_init_set_ui(m, (unsigned long)n);
        c->incomplete = 1;
        fi_warn_incomplete(c, m);
        mpz_clear(m);
        wo_add(out, n, mult);
        return;
    }
    word_factor_method_w(d, mult, out, c);
    word_factor_method_w(n / d, mult, out, c);
}

/* FactorInteger`Method[n] for a machine integer n >= 2: sorted prime
 * powers into p[]/e[] (room for 64); returns the count, or -1 when the
 * method has no word splitter (the caller takes the general path). */
int fi_word_method_array(uint64_t n, FiCtx* c, uint64_t* p, int64_t* e) {
    if (!has_word_splitter(c->method)) return -1;
    WordOut w;
    w.n = 0;
    c->rng = n * UINT64_C(0x9E3779B97F4A7C15) | 1;
    int z = 0;
    while (!(n & 1)) { n >>= 1; z++; }
    if (z) wo_add(&w, 2, z);
    word_factor_method_w(n, 1, &w, c);
    for (int i = 1; i < w.n; i++) {
        uint64_t pp = w.p[i];
        int64_t ee = w.e[i];
        int j = i;
        while (j > 0 && w.p[j - 1] > pp) { w.p[j] = w.p[j - 1]; w.e[j] = w.e[j - 1]; j--; }
        w.p[j] = pp; w.e[j] = ee;
    }
    for (int i = 0; i < w.n; i++) { p[i] = w.p[i]; e[i] = w.e[i]; }
    return w.n;
}

static void word_factor_method(uint64_t n, int64_t mult, FiFactors* out, FiCtx* c, int* gave_up) {
    if (n <= 1) return;
    if (fi_u64_is_prime(n)) { fi_factors_add_ui(out, n, mult); return; }
    uint64_t b;
    int e = word_root(n, &b);
    if (e > 1) { word_factor_method(b, mult * e, out, c, gave_up); return; }
    if (c->method == FI_PM1 || c->method == FI_PP1 || c->method == FI_ECM) {
        /* Same tiny-prime pre-step as the bignum front ends. */
        size_t np;
        const uint32_t* P = fi_small_primes(&np);
        for (size_t i = 1; i < np && P[i] < 1024; i++)
            if (n % P[i] == 0 && n != P[i]) {
                word_factor_method(P[i], mult, out, c, gave_up);
                word_factor_method(n / P[i], mult, out, c, gave_up);
                return;
            }
    }
    uint64_t d = word_split(n, c);
    if (!d) {
        mpz_t m;
        mpz_init_set_ui(m, (unsigned long)n);
        fi_factors_add(out, m, mult);
        c->incomplete = 1;
        fi_warn_incomplete(c, m);
        mpz_clear(m);
        *gave_up = 1;
        return;
    }
    word_factor_method(d, mult, out, c, gave_up);
    word_factor_method(n / d, mult, out, c, gave_up);
}

/* ------------------------------------------------------------------ */
/* Automatic tier table                                                */
/* ------------------------------------------------------------------ */

/* Trial-division bound applied once at the top of an Automatic run. */
#define AUTO_TRIAL_BOUND   (1u << 16)
/* Cheap rho budget: finds factors up to ~2^24 before p-1 / ECM. */
#define AUTO_RHO_ITERS     (1u << 13)
/* Cheap Fermat probe: catches two factors within ~n^(1/4) * 2^7 of each
 * other (a common construction, and the one rho/ECM/QS are worst at). */
#define AUTO_FERMAT_ITERS  (1u << 14)
/* One Pollard p-1 and one Williams p+1 pass (GMP-ECM picks stage-2 B2). */
#define AUTO_PM1_B1        150000.0
#define AUTO_PP1_B1        100000.0
/* QS is the finisher for composites up to this many digits. */
#define AUTO_QS_MAX_DIGITS 100
/* ECM digit levels: 15, 20, ..., stopping after AUTO_ECM_MAX. */
#define AUTO_ECM_FIRST     15
#define AUTO_ECM_STEP      5
#define AUTO_ECM_MAX       65

enum { ST_RHO = 0, ST_TINY, ST_FERMAT, ST_PM1, ST_PP1, ST_ECM0 };

/* Run tier `s` once on m.  Stage s >= ST_ECM0 is ECM level
 * AUTO_ECM_FIRST + (s - ST_ECM0) * STEP.  The p-+1 bounds scale with the
 * size of m so a small composite is not charged a large-number search. */
static int run_tier(int s, mpz_t d, const mpz_t m, size_t bits, int digits, FiCtx* c) {
    switch (s) {
    case ST_RHO:
        return bits < 126 ? fi_u128_rho(d, m, 1, AUTO_RHO_ITERS)
                          : fi_split_rho(d, m, c, AUTO_RHO_ITERS, 1);
    case ST_TINY: {
        /* Two-limb micro-ECM, a light pass: 8 curves at B1 = 150 catch the
         * 9..12 digit factors rho misses for ~0.25 ms.  Measured: a heavier
         * pass (B1 = 300..1000) loses to QS on balanced 96-bit semiprimes. */
        if (bits > 126) return 0;
        uint32_t B1 = 150;
        int curves = 8;
        FI_TUNE_U("FI_TINY_B1", B1);
        FI_TUNE_U("FI_TINY_CURVES", curves);
        return fi_u128_ecm(d, m, &c->rng, B1, curves);
    }
    case ST_FERMAT: return fi_split_fermat(d, m, c, AUTO_FERMAT_ITERS);
    case ST_PM1: {
        double b1 = 2000.0 * digits;
        return fi_split_pm1(d, m, c, b1 < AUTO_PM1_B1 ? b1 : AUTO_PM1_B1, 0);
    }
    case ST_PP1: {
        /* Each seed succeeds with probability 1/2 (it tests p+1 or p-1
         * depending on a Legendre symbol), so composites big enough that ECM
         * or QS will cost seconds (70+ digits) get three seeds; smaller ones
         * one probe, which costs less than the parallel QS it might save. */
        double b1 = 1500.0 * digits;
        return fi_split_pp1_seeds(d, m, c, b1 < AUTO_PP1_B1 ? b1 : AUTO_PP1_B1, 0, digits >= 70 ? 3 : 1);
    }
    default:        return fi_ecm_level(d, m, c, AUTO_ECM_FIRST + (s - ST_ECM0) * AUTO_ECM_STEP);
    }
}

/* Search effort spent once tiers [0, s) have run, as a digit count: the
 * size of factor they would have found with good probability.  Calibrated
 * so that QS (whose cost at 40 digits is about one p-1 pass at B1 = 10^5)
 * takes over directly after the Fermat probe below ~45 digits, after p-1
 * below ~52, after p+1 below ~56, and after ECM level L below 3.5 L. */
static int effort_digits(int s) {
    switch (s) {
    case ST_RHO:    return 0;
    case ST_TINY:   return 7;
    case ST_FERMAT: return 7;
    case ST_PM1:    return 13;
    case ST_PP1:    return 15;
    case ST_ECM0:   return 16;
    default: {
        int done = AUTO_ECM_FIRST + (s - ST_ECM0 - 1) * AUTO_ECM_STEP;
        return done > 16 ? done : 16;
    }
    }
}

/* Returns 1 with d set, 0 when the whole table is exhausted.  The stage only
 * advances when a tier FAILS: after a split both parts re-enter at the tier
 * that succeeded, so a composite rich in small factors keeps peeling them
 * off with cheap rho instead of falling through to ECM / QS.  Before every
 * tier past Fermat, QS takes over once the effort spent reaches 2/7 of the
 * composite's digit count (yafu's crossover): beyond that point a factor is
 * more cheaply found by sieving than by searching further. */
static int split_auto(mpz_t d, Work* w, FiCtx* c) {
    const mpz_t* mp = (const mpz_t*)&w->m;
    size_t bits = mpz_sizeinbase(*mp, 2);
    int digits = (int)(bits * 0.30103) + 1;
    for (;;) {
        /* An interactive interrupt: stop climbing the tiers. Each tier would
         * otherwise run at least one ECM stage 2 (which never polls) before
         * declining, so a Ctrl-C took seconds. The number is left incomplete
         * and the evaluator discards the result. A TimeConstrained deadline
         * keeps its own path (the jump after the parallel join). */
        if (mth_interrupt_should_stop()) return 0;
        int s = w->stage;
        if (s >= ST_PM1 && digits <= AUTO_QS_MAX_DIGITS && !w->qs_tried &&
            7 * effort_digits(s) >= 2 * digits) {
            w->qs_tried = 1;
            if (fi_split_qs(d, *mp, c)) return 1;
            /* QS failed (should not happen): keep searching with ECM. */
        }
        if (s >= ST_ECM0 && AUTO_ECM_FIRST + (s - ST_ECM0) * AUTO_ECM_STEP > AUTO_ECM_MAX) return 0;
        if (run_tier(s, d, *mp, bits, digits, c)) return 1;
        w->stage++;
    }
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

void fi_warn_incomplete(FiCtx* c, const mpz_t m);

static size_t distinct_primes(const FiFactors* out) { return out->n; }

void fi_factor(const mpz_t n_in, FiFactors* out, FiCtx* c) {
    if (mpz_cmp_ui(n_in, 1) <= 0) return;
    fi_ctx_seed(c, n_in);

    /* FactorInteger[n, 1]: the answer is n itself. */
    if (c->k_limit == 1) { fi_factors_add(out, n_in, 1); return; }

    mpz_t n;
    mpz_init_set(n, n_in);

    /* Shared pre-step: powers of two; Automatic and TrialDivision also take
     * every prime up to their bound here, once. */
    if (c->method == FI_AUTOMATIC && mpz_sizeinbase(n, 2) > 64) {
        fi_trial_divide(n, AUTO_TRIAL_BOUND, 1, out);
    } else if (c->method == FI_TRIAL) {
        uint64_t B = c->bound;
        if (!B) {
            /* Default: complete trial division where it is cheap (sqrt(n)
             * below 10^7), else primes up to 10^7. */
            mpz_t s; mpz_init(s); mpz_sqrt(s, n);
            B = (mpz_cmp_ui(s, 10000000u) < 0) ? mpz_get_ui(s) : 10000000u;
            mpz_clear(s);
        }
        if (B > (UINT64_C(1) << 40)) B = UINT64_C(1) << 40;
        fi_trial_divide(n, B, 1, out);
        if (mpz_cmp_ui(n, 1) > 0) {
            if (fi_is_prime(n)) fi_factors_add(out, n, 1);
            else {
                mpz_t b; mpz_init(b);
                int64_t e = primitive_root(b, n);
                if (e > 1 && fi_is_prime(b)) fi_factors_add(out, b, e);
                else { fi_factors_add(out, n, 1); c->incomplete = 1; fi_warn_incomplete(c, n); }
                mpz_clear(b);
            }
        }
        mpz_clear(n);
        return;
    } else {
        mp_bitcnt_t z = mpz_scan1(n, 0);
        if (z) { mpz_tdiv_q_2exp(n, n, z); fi_factors_add_ui(out, 2, (int64_t)z); }
        if (c->k_limit == 0 && has_word_splitter(c->method) && mpz_sizeinbase(n, 2) <= 64) {
            int gave_up = 0;
            word_factor_method(mpz_get_ui(n), 1, out, c, &gave_up);
            mpz_clear(n);
            return;
        }
    }

    WorkStack st = { NULL, 0, 0 };
    if (mpz_cmp_ui(n, 1) > 0) push(&st, n, 1, 0);
    mpz_clear(n);

    mpz_t d, q, b;
    mpz_inits(d, q, b, NULL);
    mpz_t rest;
    mpz_init_set_ui(rest, 1);
    int truncated = 0;

    while (st.n) {
        Work w = st.v[--st.n];            /* take ownership of w.m */

        if (c->k_limit > 1 && distinct_primes(out) + 1 >= (size_t)c->k_limit &&
            !(st.n == 0 && !truncated && fi_is_prime(w.m))) {
            /* Partial factorisation: everything still pending becomes the
             * final entry, the unfactored remainder (a lone prime is still
             * emitted as itself).  The caller merges any overshoot. */
            mpz_t pw; mpz_init(pw);
            mpz_pow_ui(pw, w.m, (unsigned long)w.mult);
            mpz_mul(rest, rest, pw);
            mpz_clear(pw);
            mpz_clear(w.m);
            truncated = 1;
            continue;
        }

        if (c->method == FI_AUTOMATIC && mpz_sizeinbase(w.m, 2) <= 64) {
            fi_u64_factor(mpz_get_ui(w.m), w.mult, out);
            mpz_clear(w.m);
            continue;
        }
        if (fi_is_prime(w.m)) { fi_factors_add(out, w.m, w.mult); mpz_clear(w.m); continue; }
        int64_t e = primitive_root(b, w.m);
        if (e > 1) {
            push(&st, b, w.mult * e, w.stage);
            mpz_clear(w.m);
            continue;
        }

        int r;
        if (c->method == FI_AUTOMATIC) r = split_auto(d, &w, c);
        else r = split_named(d, w.m, c);

        if (r == 1 && mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, w.m) < 0 && mpz_divisible_p(w.m, d)) {
            mpz_divexact(q, w.m, d);
            /* Push the larger part first so the smaller is processed first:
             * small factors come out early, which is what FactorInteger[n,k]
             * wants. */
            if (mpz_cmp(d, q) < 0) { push(&st, q, w.mult, w.stage); push(&st, d, w.mult, w.stage); }
            else                   { push(&st, d, w.mult, w.stage); push(&st, q, w.mult, w.stage); }
        } else {
            fi_factors_add(out, w.m, w.mult);
            c->incomplete = 1;
            fi_warn_incomplete(c, w.m);
        }
        mpz_clear(w.m);
    }
    if (truncated && mpz_cmp_ui(rest, 1) > 0) fi_factors_add(out, rest, 1);
    free(st.v);
    mpz_clears(d, q, b, rest, NULL);
}
