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

/* Heads that carry NumericFunction but do not compile.
 *
 * As of the M5 kernel work this list contains NO pending numerics — every head
 * here is one that SHOULD stay uncompiled, for one of two reasons:
 *
 *  (a) THE INTERPRETER ITSELF DECLINES on real arguments.  `BesselJZero[2., 1.]`,
 *      `BarnesG[3.5]`, `Hyperfactorial[3.5]`, `Factorial2[5.5]` and
 *      `FactorialPower[5.5, 2.]` all come back UNEVALUATED.  A machine kernel
 *      would make the compiled path ANSWER where the interpreter does not —
 *      the one divergence the engine forbids, and a far worse bug than a bail.
 *      (Check this before writing numerics for anything on this list:
 *      `Binomial[5.5, 2.]` DOES evaluate to 12.375, which is why Binomial is
 *      covered and these are not.)
 *
 *  (b) NOT A SCALAR MACHINE FUNCTION.  Exact integer semantics a double cannot
 *      represent faithfully (GCD, LCM, DigitSum), or a result that is not one
 *      number (ReIm, QuotientRemainder return lists).
 *
 * So an entry appearing here is a statement about the FUNCTION, not about
 * unfinished work.  If a genuinely-pending kernel is ever added back to this
 * list, say so explicitly — otherwise the distinction rots. */
typedef struct { const char* name; const char* why; } Gap;
static const Gap KNOWN_GAPS[] = {
    /* (a) the interpreter leaves these unevaluated on machine reals. */
    { "BesselJZero",          "interpreter leaves BesselJZero[n,k] unevaluated" },
    { "BarnesG",              "interpreter leaves BarnesG[real] unevaluated" },
    { "Hyperfactorial",       "interpreter leaves Hyperfactorial[real] unevaluated" },
    { "Factorial2",           "interpreter leaves Factorial2[real] unevaluated" },
    { "FactorialPower",       "interpreter leaves FactorialPower[real,real] unevaluated" },
    /* (b) deliberately not machine-scalar functions. */
    { "GCD",                  "exact integer semantics; not a double kernel" },
    { "LCM",                  "exact integer semantics; not a double kernel" },
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

/* A few heads cannot be probed by the generic `Head[x, y, ...]` form because
 * their real signature is not a flat list of scalars — Clip and Rescale take a
 * LIST of bounds.  Probing them with the wrong shape would report a coverage gap
 * that is really a defect in the probe. */
typedef struct { const char* name; const char* probe; } Probe;
static const Probe PROBES[] = {
    { "Clip",              "Clip[x, {1., 3.}]" },
    { "Rescale",           "Rescale[x, {0., 4.}]" },
    { "HypergeometricPFQ", "HypergeometricPFQ[{1.5}, {2.5}, x]" },
};
static const char* probe_for(const char* nm) {
    for (size_t i = 0; i < sizeof PROBES / sizeof PROBES[0]; i++)
        if (strcmp(PROBES[i].name, nm) == 0) return PROBES[i].probe;
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

/* Does a body containing this head compile, at any arity from 1 to 4?
 * Four, not three: Hypergeometric2F1[a, b, c, z] is 4-ary, and probing only to
 * three reported it as a coverage gap when the kernel was there all along. */
#define PROBE_MAXARITY 4
static bool head_compiles(const char* h) {
    static const char* an[PROBE_MAXARITY] = { "x", "y", "z", "w" };
    const char* inm[PROBE_MAXARITY];
    const CompileType RRR[PROBE_MAXARITY] = { CT_REAL, CT_REAL, CT_REAL, CT_REAL };
    for (int i = 0; i < PROBE_MAXARITY; i++) inm[i] = intern_symbol(an[i]);

    const char* override = probe_for(h);
    if (override) {
        Expr* b = parse_expression(override);
        if (!b) return false;
        CompiledProgram* prog = compile_expr(b, inm, RRR, 1);
        expr_free(b);
        if (prog) { compiled_free(prog); return true; }
        return false;
    }

    for (int arity = 1; arity <= PROBE_MAXARITY; arity++) {
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
