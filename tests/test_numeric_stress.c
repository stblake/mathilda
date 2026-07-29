/* Stress sweep: machine-precision N[] must agree with high-precision N[].
 *
 * The invariant
 *
 *     N[f[x]]  ==  N[f[x], STRESS_DIGITS]   (to within STRESS_RELTOL)
 *
 * holds for every numeric builtin and every exact argument x, and it is
 * exactly what the leaf-rounding bug violated: N[Sin[3141592653589793238]]
 * answered -0.641653 while N[Sin[3141592653589793238], 30] answered
 * -0.446315..., because machine mode rounded the argument to a double
 * before Sin ever saw it. It needs no external reference table, so it
 * applies uniformly across every elementary and special function — which
 * is what makes sweeping all of them tractable.
 *
 * This is a soundness sweep, not an accuracy audit: the tolerance is loose
 * enough that ordinary differences between a machine kernel and its MPFR
 * counterpart pass, and tight enough that the failure mode above (relative
 * errors of order 1) cannot. tests/test_numeric_largearg.c is the tight
 * one — it checks against correctly-rounded MPFR oracles.
 *
 * Combinations neither path can answer (Exp past MPFR's exponent range,
 * special functions with no large-argument implementation) are recorded in
 * kKnownGaps below and asserted to STAY there: an unlisted gap fails the
 * run, and so does a listed gap that no longer occurs. The list cannot
 * silently grow, and it cannot silently rot either.
 *
 * Run in the foreground; it is chatty on failure by design.
 */

/* clock_gettime / fork / waitpid are POSIX, which glibc hides under a
 * strict-C99 dialect. Must come before any #include to have an effect.
 * 200809L rather than 199309L: the older level predates C99 and Darwin's
 * headers then hide snprintf along with it. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#ifdef USE_MPFR
#include <mpfr.h>

/* Digits the reference side is computed at. Well above machine, low enough
 * that the slower special functions still finish. */
#define STRESS_DIGITS 30
#define STRESS_BITS   128
#define STRESS_RELTOL 1e-9

/* test_utils.h arms alarm(60) in a constructor. This sweep is ~1500
 * evaluations, a good share of them 30-digit special functions, so give it
 * a real budget. Still bounded: a hang must not become a hung CI job. */
#define STRESS_TIMEOUT_SECONDS 900

/* Per-case budget. Each case runs in a forked child so that a hang or a
 * crash is a recorded outcome instead of the end of the run — this suite
 * exists to probe hostile inputs, and one bad combination must not cost
 * the coverage of the other eight hundred. N[Erfi[2^53], 30] is the case
 * that motivated it: its series does not terminate at large arguments.
 * 10 s is ~6x the slowest case that does finish (AiryBi at -10^25, 1.8 s),
 * so it separates "slow" from "never" without much idle waiting. */
#define CASE_TIMEOUT_SECONDS 10

/* ------------------------------------------------------------------------
 *  Arguments
 *
 *  The classes that matter: exactly representable as a double (control —
 *  these were always right), not exactly representable (the bug), and
 *  ordinary small values (the fast path, which must not regress).
 * ---------------------------------------------------------------------- */

#define A_BIG     (1u << 0)   /* |x| large: skip where that overflows */
#define A_NEG     (1u << 1)   /* negative */
#define A_UNIT    (1u << 2)   /* |x| <= 1 */
#define A_SMALL   (1u << 3)   /* |x| small but > 1 */

typedef struct {
    const char* label;    /* used in the gap key, so keep it stable */
    const char* text;     /* what actually gets substituted */
    unsigned    klass;
} StressArg;

static const StressArg kArgs[] = {
    /* exactly representable — the control */
    { "2^53",       "9007199254740992",                          A_BIG },
    { "10^20",      "100000000000000000000",                     A_BIG },
    /* not exactly representable — the bug */
    { "2^53+1",     "9007199254740993",                          A_BIG },
    { "pi-digits",  "3141592653589793238",                       A_BIG },
    { "10^25",      "10000000000000000000000000",                A_BIG },
    { "10^40",      "10000000000000000000000000000000000000000", A_BIG },
    { "-10^25",     "-10000000000000000000000000",     A_BIG | A_NEG },
    /* rationals: never exactly representable, and the amplification is
     * driven by the magnitude, not by the size of numerator and
     * denominator */
    { "10^20/3",    "100000000000000000000/3",                   A_BIG },
    { "22/7",       "22/7",                                    A_SMALL },
    { "1/3",        "1/3",                                      A_UNIT },
    /* ordinary controls */
    { "2",          "2",                                       A_SMALL },
    { "-3/2",       "-3/2",                            A_SMALL | A_NEG },
    { "1/2",        "1/2",                                      A_UNIT },
};
static const size_t kArgCount = sizeof(kArgs) / sizeof(kArgs[0]);

/* ------------------------------------------------------------------------
 *  Functions
 * ---------------------------------------------------------------------- */

#define F_ALL     0u
#define F_POS     (1u << 0)   /* argument must be > 0 */
#define F_NOBIG   (1u << 1)   /* argument is an index/order: keep it small */
#define F_UNIT    (1u << 2)   /* argument must satisfy |x| <= 1 */

typedef struct {
    const char* label;    /* gap key; keep stable */
    const char* form;     /* printf template, one %s */
    unsigned    flags;
} StressFn;

static const StressFn kFns[] = {
    /* --- trigonometric --- */
    { "Sin",        "Sin[%s]",        F_ALL },
    { "Cos",        "Cos[%s]",        F_ALL },
    { "Tan",        "Tan[%s]",        F_ALL },
    { "Cot",        "Cot[%s]",        F_ALL },
    { "Sec",        "Sec[%s]",        F_ALL },
    { "Csc",        "Csc[%s]",        F_ALL },
    /* --- inverse trigonometric (complex off the principal interval) --- */
    { "ArcSin",     "ArcSin[%s]",     F_ALL },
    { "ArcCos",     "ArcCos[%s]",     F_ALL },
    { "ArcTan",     "ArcTan[%s]",     F_ALL },
    { "ArcCot",     "ArcCot[%s]",     F_ALL },
    { "ArcSec",     "ArcSec[%s]",     F_ALL },
    { "ArcCsc",     "ArcCsc[%s]",     F_ALL },
    /* --- hyperbolic --- */
    { "Sinh",       "Sinh[%s]",       F_ALL },
    { "Cosh",       "Cosh[%s]",       F_ALL },
    { "Tanh",       "Tanh[%s]",       F_ALL },
    { "Coth",       "Coth[%s]",       F_ALL },
    { "Sech",       "Sech[%s]",       F_ALL },
    { "Csch",       "Csch[%s]",       F_ALL },
    { "ArcSinh",    "ArcSinh[%s]",    F_ALL },
    { "ArcCosh",    "ArcCosh[%s]",    F_ALL },
    { "ArcTanh",    "ArcTanh[%s]",    F_ALL },
    { "ArcCoth",    "ArcCoth[%s]",    F_ALL },
    /* --- exponential / algebraic --- */
    { "Log",        "Log[%s]",        F_ALL },
    { "Exp",        "Exp[%s]",        F_ALL },
    { "Sqrt",       "Sqrt[%s]",       F_ALL },
    { "Abs",        "Abs[%s]",        F_ALL },
    { "Sign",       "Sign[%s]",       F_ALL },
    { "Sinc",       "Sinc[%s]",       F_ALL },
    /* --- gamma family --- */
    { "Gamma",      "Gamma[%s]",      F_ALL },
    { "LogGamma",   "LogGamma[%s]",   F_POS },
    { "PolyGamma",  "PolyGamma[%s]",  F_POS },
    { "PolyGamma1", "PolyGamma[1, %s]", F_POS },
    { "Pochhammer", "Pochhammer[%s, 2]", F_ALL },
    { "Beta",       "Beta[%s, 2]",    F_POS },
    /* --- zeta family --- */
    { "Zeta",       "Zeta[%s]",       F_ALL },
    { "HurwitzZeta", "HurwitzZeta[2, %s]", F_POS },
    { "PolyLog",    "PolyLog[2, %s]", F_ALL },
    { "LerchPhi",   "LerchPhi[1/2, 2, %s]", F_POS },
    { "HarmonicNumber", "HarmonicNumber[%s]", F_ALL },
    /* --- error function family --- */
    { "Erf",        "Erf[%s]",        F_ALL },
    { "Erfc",       "Erfc[%s]",       F_ALL },
    { "Erfi",       "Erfi[%s]",       F_ALL },
    { "InverseErf", "InverseErf[%s]", F_UNIT },
    { "FresnelC",   "FresnelC[%s]",   F_ALL },
    { "FresnelS",   "FresnelS[%s]",   F_ALL },
    /* --- integral functions --- */
    { "SinIntegral",  "SinIntegral[%s]",  F_ALL },
    { "CosIntegral",  "CosIntegral[%s]",  F_POS },
    { "SinhIntegral", "SinhIntegral[%s]", F_ALL },
    { "CoshIntegral", "CoshIntegral[%s]", F_POS },
    { "ExpIntegralEi", "ExpIntegralEi[%s]", F_ALL },
    { "LogIntegral",  "LogIntegral[%s]",  F_POS },
    /* --- Airy / Bessel / Legendre --- */
    { "AiryAi",     "AiryAi[%s]",     F_ALL },
    { "AiryBi",     "AiryBi[%s]",     F_ALL },
    { "AiryAiPrime", "AiryAiPrime[%s]", F_ALL },
    { "BesselJ0",   "BesselJ[0, %s]", F_ALL },
    { "BesselY0",   "BesselY[0, %s]", F_POS },
    { "BesselI0",   "BesselI[0, %s]", F_ALL },
    { "BesselK0",   "BesselK[0, %s]", F_POS },
    { "BesselJfrac", "BesselJ[1/3, %s]", F_POS },
    { "LegendreP2", "LegendreP[2, %s]", F_ALL },
    /* --- Lambert W and friends --- */
    { "ProductLog", "ProductLog[%s]", F_POS },
    /* --- index-taking: keep the argument small --- */
    { "BernoulliB", "BernoulliB[%s]", F_NOBIG },
    { "EulerE",     "EulerE[%s]",     F_NOBIG },
    { "Fibonacci",  "Fibonacci[%s]",  F_NOBIG },
    { "Binomial",   "Binomial[%s, 2]", F_NOBIG },
    { "Factorial",  "Factorial[%s]",  F_NOBIG },
};
static const size_t kFnCount = sizeof(kFns) / sizeof(kFns[0]);

/* ------------------------------------------------------------------------
 *  Known gaps
 *
 *  "<function label>@<argument label>" for combinations where BOTH paths
 *  decline in the same way. Each carries the reason. Anything not listed
 *  here that shows up is a failure, and anything listed here that no
 *  longer shows up is a failure too.
 * ---------------------------------------------------------------------- */

typedef struct {
    const char* key;
    const char* why;
} KnownGap;

static const KnownGap kKnownGaps[] = {
#include "test_numeric_stress_gaps.inc"
};
static const size_t kKnownGapCount = sizeof(kKnownGaps) / sizeof(kKnownGaps[0]);
static bool gap_seen[sizeof(kKnownGaps) / sizeof(kKnownGaps[0])];

static int find_gap(const char* key) {
    for (size_t i = 0; i < kKnownGapCount; ++i)
        if (strcmp(kKnownGaps[i].key, key) == 0) return (int)i;
    return -1;
}

/* ------------------------------------------------------------------------
 *  Result extraction
 * ---------------------------------------------------------------------- */

typedef enum {
    RES_NUMBER,      /* a finite real or complex value */
    RES_DEGENERATE,  /* Infinity / NaN / ComplexInfinity / Indeterminate */
    RES_SYMBOLIC     /* left unevaluated */
} ResKind;

static bool scalar_to_mpfr(const Expr* e, mpfr_t out, bool* finite) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); break;
        case EXPR_BIGINT:  mpfr_set_z(out, e->data.bigint, MPFR_RNDN);         break;
        case EXPR_REAL:    mpfr_set_d(out, e->data.real, MPFR_RNDN);           break;
        case EXPR_MPFR:    mpfr_set(out, e->data.mpfr, MPFR_RNDN);             break;
        case EXPR_FUNCTION: {
            /* Rational[n, d] survives when a builtin returns an exact value. */
            const Expr* h = e->data.function.head;
            if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 2
                && strcmp(h->data.symbol.name, "Rational") == 0) {
                mpfr_t d;
                mpfr_init2(d, mpfr_get_prec(out));
                bool ok = scalar_to_mpfr(e->data.function.args[0], out, NULL)
                       && scalar_to_mpfr(e->data.function.args[1], d, NULL);
                if (ok) mpfr_div(out, out, d, MPFR_RNDN);
                mpfr_clear(d);
                if (!ok) return false;
                break;
            }
            return false;
        }
        default: return false;
    }
    if (finite) *finite = mpfr_number_p(out) != 0;
    return true;
}

static ResKind classify(const Expr* e, mpfr_t re, mpfr_t im) {
    mpfr_set_zero(re, 1);
    mpfr_set_zero(im, 1);
    if (!e) return RES_SYMBOLIC;

    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        if (strcmp(n, "Infinity") == 0 || strcmp(n, "ComplexInfinity") == 0
            || strcmp(n, "Indeterminate") == 0 || strcmp(n, "Underflow") == 0
            || strcmp(n, "Overflow") == 0) {
            return RES_DEGENERATE;
        }
        return RES_SYMBOLIC;
    }

    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h && h->type == EXPR_SYMBOL && e->data.function.arg_count == 2
            && strcmp(h->data.symbol.name, "Complex") == 0) {
            bool fr = true, fi = true;
            if (!scalar_to_mpfr(e->data.function.args[0], re, &fr)) return RES_SYMBOLIC;
            if (!scalar_to_mpfr(e->data.function.args[1], im, &fi)) return RES_SYMBOLIC;
            return (fr && fi) ? RES_NUMBER : RES_DEGENERATE;
        }
        if (h && h->type == EXPR_SYMBOL
            && strcmp(h->data.symbol.name, "DirectedInfinity") == 0) {
            return RES_DEGENERATE;
        }
    }

    bool fin = true;
    if (!scalar_to_mpfr(e, re, &fin)) return RES_SYMBOLIC;
    return fin ? RES_NUMBER : RES_DEGENERATE;
}

/* ------------------------------------------------------------------------
 *  The sweep
 * ---------------------------------------------------------------------- */

static int cases = 0, checks = 0, failures = 0, gaps = 0, skipped = 0;
static int unlisted_gaps = 0;

/* A sweep this wide will eventually meet a combination that takes minutes.
 * Without this the only symptom is a silent SIGALRM with no indication of
 * which case was running, so record the case in progress and have the
 * timeout name it. Async-signal-safe: write(2) on a buffer already filled. */
static char current_case[256] = "<none>";

static void on_timeout(int sig) {
    (void)sig;
    static const char msg[] = "\nFAIL: timed out in case: ";
    ssize_t ignored;
    ignored = write(STDERR_FILENO, msg, sizeof(msg) - 1);
    ignored = write(STDERR_FILENO, current_case, strlen(current_case));
    ignored = write(STDERR_FILENO, "\n", 1);
    (void)ignored;
    _exit(1);
}

/* Wall time, not CPU time: the work happens in a forked child, so clock()
 * in the parent would report zero for every case. */
static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Slowest cases, reported at the end. A case that has quietly become slow
 * is worth knowing about before it becomes a timeout. */
#define SLOW_TOP 12
static struct { char key[64]; double secs; } slowest[SLOW_TOP];

static void note_time(const char* key, double secs) {
    size_t worst = 0;
    for (size_t i = 1; i < SLOW_TOP; ++i)
        if (slowest[i].secs < slowest[worst].secs) worst = i;
    if (secs > slowest[worst].secs) {
        snprintf(slowest[worst].key, sizeof(slowest[worst].key), "%s", key);
        slowest[worst].secs = secs;
    }
}

static void report(const char* key, const char* what, const char* detail) {
    fprintf(stderr, "FAIL: %s [%s]\n  %s\n", what, key, detail);
    failures++;
}

static Expr* eval_str(const char* input) {
    Expr* parsed = parse_expression(input);
    if (!parsed) return NULL;
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    return r;
}

static bool arg_allowed(const StressFn* f, const StressArg* a) {
    if ((f->flags & F_POS)   && (a->klass & A_NEG))    return false;
    if ((f->flags & F_NOBIG) && (a->klass & A_BIG))    return false;
    if ((f->flags & F_NOBIG) && !(a->klass & A_SMALL)) return false;
    if ((f->flags & F_UNIT)  && !(a->klass & A_UNIT))  return false;
    return true;
}

/* What one case concluded. The first character of the child's reply. */
#define OUT_OK      'O'   /* both paths agree */
#define OUT_MISMATCH 'M'  /* both computed, and they disagree */
#define OUT_GAP     'G'   /* neither path produced a finite number */
#define OUT_SPLIT   'S'   /* exactly one path computed: a real divergence */

/* Runs entirely in the forked child. Writes "<outcome>|<detail>" to `fd`. */
static void decide_case(const char* lo, const char* hi, int fd) {
    Expr* elo = eval_str(lo);
    Expr* ehi = eval_str(hi);
    char* slo = elo ? expr_to_string(elo) : NULL;
    char* shi = ehi ? expr_to_string(ehi) : NULL;
    char reply[1024];

    mpfr_t lre, lim, hre, him;
    mpfr_inits2(STRESS_BITS, lre, lim, hre, him, (mpfr_ptr)0);
    ResKind klo = classify(elo, lre, lim);
    ResKind khi = classify(ehi, hre, him);

    if (klo != RES_NUMBER || khi != RES_NUMBER) {
        /* Both declining is a coverage gap even when they decline in
         * different words (unevaluated vs ComplexInfinity). Only one side
         * answering is a genuine divergence between the two paths. */
        char kind = (klo == RES_NUMBER || khi == RES_NUMBER) ? OUT_SPLIT : OUT_GAP;
        snprintf(reply, sizeof(reply), "%c|N[..] = %s ; N[.., %d] = %s",
                 kind, slo ? slo : "(null)", STRESS_DIGITS, shi ? shi : "(null)");
    } else {
        /* scale = |Re(want)| + |Im(want)|, so a component that is
         * legitimately zero is judged against the size of the value as a
         * whole rather than against itself. */
        mpfr_t scale, diff, t;
        mpfr_inits2(STRESS_BITS, scale, diff, t, (mpfr_ptr)0);
        mpfr_abs(scale, hre, MPFR_RNDN);
        mpfr_abs(t, him, MPFR_RNDN);
        mpfr_add(scale, scale, t, MPFR_RNDN);
        mpfr_mul_d(scale, scale, STRESS_RELTOL, MPFR_RNDN);

        bool ok = true;
        mpfr_sub(diff, lre, hre, MPFR_RNDN);
        mpfr_abs(diff, diff, MPFR_RNDN);
        if (mpfr_greater_p(diff, scale)) ok = false;
        mpfr_sub(t, lim, him, MPFR_RNDN);
        mpfr_abs(t, t, MPFR_RNDN);
        if (mpfr_greater_p(t, scale)) ok = false;
        mpfr_clears(scale, diff, t, (mpfr_ptr)0);

        snprintf(reply, sizeof(reply), "%c|N[..] = %s ; N[.., %d] = %s",
                 ok ? OUT_OK : OUT_MISMATCH,
                 slo ? slo : "(null)", STRESS_DIGITS, shi ? shi : "(null)");
    }

    ssize_t ignored = write(fd, reply, strlen(reply));
    (void)ignored;
    mpfr_clears(lre, lim, hre, him, (mpfr_ptr)0);
    free(slo);
    free(shi);
    expr_free(elo);
    expr_free(ehi);
}

static void handle_gap(const char* key, const char* what, const char* detail) {
    int idx = find_gap(key);
    if (idx < 0) {
        fprintf(stderr,
                "FAIL: unlisted gap [%s]\n    %s\n    %s\n"
                "    Add it to test_numeric_stress_gaps.inc if intended.\n",
                key, what, detail);
        unlisted_gaps++;
        failures++;
    } else {
        gap_seen[idx] = true;
        gaps++;
    }
}

static void run_case(const StressFn* f, const StressArg* a) {
    char key[128], body[256], lo[300], hi[320];
    snprintf(key,  sizeof(key),  "%s@%s", f->label, a->label);
    snprintf(body, sizeof(body), f->form, a->text);
    snprintf(lo,   sizeof(lo),   "N[%s]", body);
    snprintf(hi,   sizeof(hi),   "N[%s, %d]", body, STRESS_DIGITS);
    cases++;
    snprintf(current_case, sizeof(current_case), "%s  (%s)", key, hi);

    int fds[2];
    ASSERT_MSG(pipe(fds) == 0, "pipe() failed for %s", key);

    double t0 = now_seconds();
    fflush(NULL);   /* or the child re-flushes our buffered output at _exit */
    pid_t pid = fork();
    ASSERT_MSG(pid >= 0, "fork() failed for %s", key);
    if (pid == 0) {
        close(fds[0]);
        signal(SIGALRM, SIG_DFL);       /* so the parent sees the signal status */
        alarm(CASE_TIMEOUT_SECONDS);
        decide_case(lo, hi, fds[1]);
        close(fds[1]);
        _exit(0);
    }
    close(fds[1]);

    char buf[1024];
    size_t got = 0;
    for (;;) {
        ssize_t n = read(fds[0], buf + got, sizeof(buf) - 1 - got);
        if (n > 0)      got += (size_t)n;
        else if (n == 0) break;
        else if (errno != EINTR) break;
    }
    buf[got] = '\0';
    close(fds[0]);

    int status = 0;
    while (waitpid(pid, &status, 0) < 0 && errno == EINTR) { }
    note_time(key, now_seconds() - t0);

    const char* detail = (got > 2) ? buf + 2 : "(no output)";

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0 || got < 2) {
        /* Killed, timed out, or died before answering. Still an outcome,
         * still has to be a listed gap — silence here would read as
         * coverage. */
        char why[256];
        if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM) {
            snprintf(why, sizeof(why), "TIMED OUT after %ds  (%s)",
                     CASE_TIMEOUT_SECONDS, hi);
        } else if (WIFSIGNALED(status)) {
            snprintf(why, sizeof(why), "died on signal %d  (%s)",
                     WTERMSIG(status), hi);
        } else {
            snprintf(why, sizeof(why), "exited %d without answering  (%s)",
                     WEXITSTATUS(status), hi);
        }
        handle_gap(key, lo, why);
        return;
    }

    switch (buf[0]) {
        case OUT_OK:
            checks++;
            break;
        case OUT_MISMATCH: {
            char msg[1200];
            snprintf(msg, sizeof(msg),
                     "machine N disagrees with %d-digit N\n    %s",
                     STRESS_DIGITS, detail);
            checks++;
            report(key, lo, msg);
            break;
        }
        case OUT_SPLIT: {
            char msg[1200];
            snprintf(msg, sizeof(msg),
                     "one path computed a number and the other did not — "
                     "they are the same computation at different precision\n"
                     "    %s", detail);
            report(key, lo, msg);
            break;
        }
        default:
            handle_gap(key, lo, detail);
            break;
    }
}

static void test_machine_matches_high_precision(void) {
    for (size_t i = 0; i < kFnCount; ++i) {
        for (size_t j = 0; j < kArgCount; ++j) {
            if (!arg_allowed(&kFns[i], &kArgs[j])) { skipped++; continue; }
            run_case(&kFns[i], &kArgs[j]);
        }
    }
}

/* Every recorded gap must still be a gap. A stale entry means the
 * limitation was fixed and the list is now lying about coverage. */
static void test_no_stale_gaps(void) {
    for (size_t i = 0; i < kKnownGapCount; ++i) {
        if (!gap_seen[i]) {
            fprintf(stderr,
                    "FAIL: stale gap [%s] — it now evaluates on both paths.\n"
                    "    Remove it from test_numeric_stress_gaps.inc (%s).\n",
                    kKnownGaps[i].key, kKnownGaps[i].why);
            failures++;
        }
    }
}

/* Identities that must hold at machine precision for arguments a double
 * cannot represent. These fail loudly under the leaf-rounding bug because
 * each side reduces a differently-rounded argument. */
static void test_identities_at_large_arguments(void) {
    static const char* ids[] = {
        "Abs[Sin[10^25]^2 + Cos[10^25]^2 - 1] < 10^-12",
        "Abs[Tan[10^25] - Sin[10^25]/Cos[10^25]] < 10^-12",
        "Abs[Cot[10^25] Tan[10^25] - 1] < 10^-12",
        "Abs[Sec[10^25] Cos[10^25] - 1] < 10^-12",
        "Abs[Csc[10^25] Sin[10^25] - 1] < 10^-12",
        "Abs[Sin[2 (2^53+1)] - 2 Sin[2^53+1] Cos[2^53+1]] < 10^-12",
        "Abs[Cosh[22/7]^2 - Sinh[22/7]^2 - 1] < 10^-12",
        "Abs[Exp[Log[10^25]] / 10^25 - 1] < 10^-12",
        "Abs[Sqrt[10^40]^2 / 10^40 - 1] < 10^-12",
        "Abs[Log[10^25 10^40] - Log[10^25] - Log[10^40]] < 10^-12",
        "Abs[Sinc[3141592653589793238] 3141592653589793238"
        " - Sin[3141592653589793238]] < 10^-12",
    };
    for (size_t i = 0; i < sizeof(ids) / sizeof(ids[0]); ++i) {
        Expr* r = eval_str(ids[i]);
        char* s = r ? expr_to_string(r) : NULL;
        checks++;
        if (!s || strcmp(s, "True") != 0) {
            fprintf(stderr, "FAIL: identity did not hold at machine precision\n"
                            "    %s\n    -> %s\n", ids[i], s ? s : "(null)");
            failures++;
        }
        free(s);
        expr_free(r);
    }
}

/* The whole point of the ladder's control rows: arguments that ARE exactly
 * representable must still take the plain machine path. Their answers were
 * already right, and routing everything through MPFR instead of just the
 * lossy cases would be a silent 5x on ordinary work. */
static void test_representable_arguments_are_unchanged(void) {
    static const struct { const char* input; const char* expect; } fixed[] = {
        { "N[Sin[10^20]]", "-0.645251" },
        { "N[Tan[10^20]]", "-0.844602" },
        { "N[Sin[2^53]]",  "-0.848926" },
        { "N[Sin[1/3]]",   "0.327195"  },
        { "N[Cos[1/2]]",   "0.877583"  },
        { "N[1/3]",        "0.333333"  },
        { "N[Sqrt[2]]",    "1.41421"   },
    };
    for (size_t i = 0; i < sizeof(fixed) / sizeof(fixed[0]); ++i) {
        Expr* r = eval_str(fixed[i].input);
        char* s = r ? expr_to_string(r) : NULL;
        checks++;
        if (!s || strcmp(s, fixed[i].expect) != 0) {
            fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                    fixed[i].input, fixed[i].expect, s ? s : "(null)");
            failures++;
        }
        free(s);
        expr_free(r);
    }
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_MPFR
    signal(SIGALRM, on_timeout);
    alarm(STRESS_TIMEOUT_SECONDS);

    TEST(test_machine_matches_high_precision);
    TEST(test_no_stale_gaps);
    TEST(test_identities_at_large_arguments);
    TEST(test_representable_arguments_are_unchanged);

    printf("\nslowest cases:\n");
    for (size_t i = 0; i < SLOW_TOP; ++i) {
        if (slowest[i].secs > 0.05)
            printf("  %-28s %6.2f s\n", slowest[i].key, slowest[i].secs);
    }

    printf("\n%d cases, %d numeric checks, %d known gaps, %d skipped, "
           "%d failures\n", cases, checks, gaps, skipped, failures);
    if (unlisted_gaps) {
        fprintf(stderr, "%d unlisted gap(s) — see the FAIL lines above.\n",
                unlisted_gaps);
    }
    if (failures) {
        fprintf(stderr, "STRESS FAILED\n");
        return 1;
    }
#else
    printf("USE_MPFR is off: the high-precision reference side is "
           "unavailable, nothing to sweep.\n");
#endif

    printf("All numeric_stress tests passed.\n");
    return 0;
}
