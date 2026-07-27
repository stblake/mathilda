/* Compile[] coverage audit over every NumericFunction head.
 *
 * WHY THIS IS A TEST AND NOT A REPORT.  When `compile_expr` cannot handle a
 * head it returns NULL, the caller quietly interprets, and the result is still
 * correct — just an order of magnitude slower.  A coverage gap is therefore
 * invisible: it looks exactly like a working fast path.  The compilable subset
 * is a cliff, not a slope, and one head just outside it costs the WHOLE body.
 *
 * So coverage is asserted here rather than described somewhere.  Every symbol
 * carrying ATTR_NUMERICFUNCTION is probed by actually compiling `Head[x]` (and
 * `Head[x, y]`), which is the only question that matters: not "is there a kernel
 * registered" but "does a body containing this head compile".
 *
 * KNOWN_GAPS lists the heads that legitimately do not compile yet, with the
 * reason.  The test fails if a head OUTSIDE that list stops compiling (a
 * regression) or if a head INSIDE it starts compiling (the list is stale and
 * should shrink).  Both directions matter: a silently stale exception list is
 * how a coverage audit stops meaning anything.
 *
 * Run: ./compile_coverage_tests */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "sym_intern.h"
#include "symtab.h"
#include "attr.h"
#include "compile/compile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

/* Heads that carry NumericFunction but are not expected to compile yet.
 *
 * Two distinct reasons, and only the first is a coverage gap worth closing:
 *
 *  (a) NO MACHINE KERNEL.  The module computes in MPFR at >= 53 bits and rounds
 *      at the end (see special_functions/zeta.c), so there is no double-precision
 *      code path to lift into the kernel registry — closing these means writing
 *      genuinely new double/double-complex implementations, each needing a parity
 *      test against the MPFR path.  Registered as explicit degrade sentinels in
 *      ndkernels.c, or absent from the registry entirely.
 *
 *  (b) NOT A SCALAR MACHINE FUNCTION.  Exact integer semantics that a double
 *      cannot represent faithfully (GCD, Binomial on big arguments), or a result
 *      that is not one number (ReIm, QuotientRemainder return lists).  These
 *      should stay uncompiled; they are listed so the audit stays honest about
 *      why, not because they are pending work. */
typedef struct { const char* name; const char* why; } Gap;
static const Gap KNOWN_GAPS[] = {
    /* (a) no machine kernel yet — real coverage gaps.
     *
     * The exponential-integral family used to live here; expint_machine.c now
     * provides real double implementations for it, so the ones remaining are
     * those that still need genuinely new numerics rather than a wrapper. */
    { "Zeta",                 "no double kernel (MPFR-only module)" },
    { "HurwitzZeta",          "no double kernel (MPFR-only module)" },
    { "PolyLog",              "no double kernel (MPFR-only module)" },
    { "LerchPhi",             "no double kernel (MPFR-only module)" },
    { "PolyGamma",            "no double kernel (MPFR-only module)" },
    { "Erfi",                 "no double kernel (MPFR-only module)" },
    { "FresnelS",             "no double kernel (MPFR-only module)" },
    { "FresnelC",             "no double kernel (MPFR-only module)" },
    { "AiryAi",               "no double kernel (MPFR-only module)" },
    { "AiryBi",               "no double kernel (MPFR-only module)" },
    { "AiryAiPrime",          "no double kernel (MPFR-only module)" },
    { "AiryBiPrime",          "no double kernel (MPFR-only module)" },
    { "ProductLog",           "no double kernel (MPFR-only module)" },
    { "BesselI",              "no double kernel (MPFR-only module)" },
    { "BesselK",              "no double kernel (MPFR-only module)" },
    { "BesselJZero",          "no double kernel (root-finding module)" },
    { "LegendreP",            "no double kernel (MPFR-only module)" },
    { "BarnesG",              "no double kernel" },
    { "Hyperfactorial",       "no double kernel" },
    { "HarmonicNumber",       "no double kernel" },
    { "Pochhammer",           "no double kernel" },
    { "QPochhammer",          "no double kernel" },
    { "Hypergeometric0F1",    "no double kernel" },
    { "Hypergeometric1F1",    "no double kernel" },
    { "Hypergeometric2F1",    "no double kernel" },
    { "HypergeometricPFQ",    "no double kernel" },
    { "UnitStep",             "no double kernel registered" },
    { "Clip",                 "n-ary; kernel registry is unary/binary only" },
    { "Rescale",              "n-ary; kernel registry is unary/binary only" },
    { "Fibonacci",            "no double kernel" },
    { "LucasL",               "no double kernel" },
    { "Factorial2",           "no double kernel" },
    { "FactorialPower",       "no double kernel" },
    /* (b) deliberately not machine-scalar functions */
    { "GCD",                  "exact integer semantics; not a double kernel" },
    { "LCM",                  "exact integer semantics; not a double kernel" },
    { "Binomial",             "exact integer semantics; not a double kernel" },
    { "DigitSum",             "exact integer semantics; not a double kernel" },
    { "ReIm",                 "returns a list, not one machine number" },
    { "QuotientRemainder",    "returns a list, not one machine number" },
};
static const size_t NGAPS = sizeof KNOWN_GAPS / sizeof KNOWN_GAPS[0];

static const char* gap_reason(const char* nm) {
    for (size_t i = 0; i < NGAPS; i++)
        if (strcmp(KNOWN_GAPS[i].name, nm) == 0) return KNOWN_GAPS[i].why;
    return NULL;
}

/* Collected NumericFunction heads. */
#define MAXH 512
static const char* heads[MAXH];
static int nheads = 0;

static void collect(const char* name, SymbolDef* def, void* user) {
    (void)user;
    if (nheads >= MAXH) return;
    if (!(get_attributes_def(def) & ATTR_NUMERICFUNCTION)) return;
    heads[nheads++] = name;
}

static int cmp_name(const void* a, const void* b) {
    return strcmp(*(const char* const*)a, *(const char* const*)b);
}

/* Does a body containing this head compile, at any arity from 1 to 3? */
static bool head_compiles(const char* h) {
    static const char* an[3] = { "x", "y", "z" };
    const char* inm[3];
    const CompileType RRR[3] = { CT_REAL, CT_REAL, CT_REAL };
    for (int i = 0; i < 3; i++) inm[i] = intern_symbol(an[i]);

    for (int arity = 1; arity <= 3; arity++) {
        char buf[256];
        int p = snprintf(buf, sizeof buf, "%s[x", h);
        for (int i = 1; i < arity; i++) p += snprintf(buf + p, sizeof buf - (size_t)p, ", %s", an[i]);
        snprintf(buf + p, sizeof buf - (size_t)p, "]");

        /* The RAW parse tree: evaluating first would let the interpreter rewrite
         * the very head under test (Sqrt[x] -> x^(1/2), Abs of a symbol, ...) and
         * the audit would be measuring something else. */
        Expr* b = parse_expression(buf);
        if (!b) continue;
        CompiledProgram* prog = compile_expr(b, inm, RRR, (size_t)arity);
        expr_free(b);
        if (prog) { compiled_free(prog); return true; }
    }
    return false;
}

int main(void) {
    core_init();
    symtab_for_each(collect, NULL);
    qsort(heads, (size_t)nheads, sizeof heads[0], cmp_name);

    int covered = 0, gaps = 0, stale = 0, regressed = 0;
    for (int i = 0; i < nheads; i++) {
        const char* h = heads[i];
        bool ok = head_compiles(h);
        const char* why = gap_reason(h);
        if (ok && why) {
            printf("FAIL: %-22s now COMPILES but is still listed as a known gap (%s)\n"
                   "      -> remove it from KNOWN_GAPS in this file\n", h, why);
            failures++; stale++;
        } else if (ok) {
            covered++;
        } else if (why) {
            gaps++;
        } else {
            printf("FAIL: %-22s does NOT compile and is not a known gap\n"
                   "      -> a body using it silently drops to the interpreter\n", h);
            failures++; regressed++;
        }
    }

    printf("\n%d NumericFunction heads: %d compile, %d known gaps"
           " (%.0f%% covered)\n",
           nheads, covered, gaps,
           nheads ? 100.0 * covered / nheads : 0.0);
    if (stale)     printf("%d stale exception(s): coverage improved, update KNOWN_GAPS.\n", stale);
    if (regressed) printf("%d head(s) REGRESSED out of the compilable subset.\n", regressed);
    if (!failures) printf("Compile coverage audit passed.\n");
    return failures ? 1 : 0;
}
