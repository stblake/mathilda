/* test_parallelmixedtower.c
 *
 * Tests for the ParallelMixedTower integration method
 * (Integrate`ParallelMixedTower, Method -> "ParallelMixedTower") -- a parallel
 * (Risch-Norman) integrator over a simple radical in a mixed transcendental
 * tower, ported from src/internal/mixed/ParallelMixed.m -- together with the
 * three Mathilda-core pieces the port needed:
 *   - the Message subsystem (Quiet / Check / Message),
 *   - association element assignment (assoc[key] = val, nested),
 *   - single-position Part assignment of a list value (m[[2]] = {1,2,3}).
 *
 * Correctness of an antiderivative is asserted by the universal predicate
 * Simplify[D[Integrate[f, x], x] - f] === 0 rather than by a fixed output
 * string, so the tests survive surface-form changes.  The method is lazily
 * loaded from disk on first use; the resolver is CWD-independent.
 *
 * Assertions use ASSERT_STR_EQ (a hard exit(1) on mismatch), NOT
 * assert_eval_eq, whose libc assert() is a no-op under -DNDEBUG -- a mismatch
 * there would print FAIL yet let the suite report success.
 */

/* dup/dup2/fileno (stderr capture for the Integrate::nonelem warning test) are
 * POSIX; expose them under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "print.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>   /* dup, dup2, fileno */

/* Evaluate `input`, hard-assert its printed form equals `expected`. */
static void assert_eval(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    ASSERT_STR_EQ(str, expected);
    free(str);
    expr_free(evaluated);
}

/* Evaluate `input` with stderr redirected to a temp file, and report whether
 * `needle` appears in what was written there.  The Integrate::nonelem warning is
 * a raw fprintf(stderr, ...) (deliberately not routed through Quiet[]/Check[],
 * matching RischTranscendental), so this is the only way to assert it fired.
 * dup/dup2 to a tmpfile (not freopen to /dev/tty) so the restore works in CI. */
static int eval_stderr_contains(const char* input, const char* needle) {
    fflush(stderr);
    int saved = dup(fileno(stderr));
    ASSERT(saved != -1);
    FILE* cap = tmpfile();
    ASSERT(cap != NULL);
    dup2(fileno(cap), fileno(stderr));

    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    expr_free(evaluated);

    fflush(stderr);
    dup2(saved, fileno(stderr));   /* restore the real stderr */
    close(saved);

    rewind(cap);
    char buf[8192];
    size_t n = fread(buf, 1, sizeof(buf) - 1, cap);
    buf[n] = '\0';
    fclose(cap);
    return strstr(buf, needle) != NULL;
}

/* ---------------------------------------------------------- Message subsystem */
static void test_message_subsystem(void) {
    /* Quiet evaluates and returns the value, suppressing the diagnostic. */
    assert_eval("Quiet[1/0]", "ComplexInfinity");
    assert_eval("Quiet[Sqrt[4]]", "2");
    /* Check returns failexpr iff a message fired during the first argument. */
    assert_eval("Check[2 + 2, bad]", "4");
    assert_eval("Check[1/0, bad]", "bad");
    assert_eval("Quiet[Check[1/0, FAILED]]", "FAILED");
    assert_eval("Check[Sqrt[-4], bad]", "2*I");   /* no message -> value */
    /* An explicit Message is caught by an enclosing Check. */
    assert_eval("foo::bar = \"m\"; Check[Message[foo::bar]; 99, CAUGHT]", "CAUGHT");
}

/* --------------------------------------------- assoc[key]=val and Part fixes */
static void test_core_assignment_fixes(void) {
    /* Association element assignment on a symbol holding an Association. */
    assert_eval("a1 = <||>; a1[\"k\"] = 5; a1[\"k\"]", "5");
    assert_eval("a2 = <||>; a2[foo[1, 2]] = 7; KeyExistsQ[a2, foo[1, 2]]", "True");
    /* Nested element assignment into an existing sub-association. */
    assert_eval("a3 = <|\"o\" -> <||>|>; a3[\"o\", \"i\"] = 9; a3[\"o\"][\"i\"]", "9");
    /* An ordinary DownValue definition is untouched (symbol has no value). */
    assert_eval("Clear[g1]; g1[1] = 5; g1[1]", "5");
    /* Single-position Part assignment stores the WHOLE list, not its head. */
    assert_eval("m1 = {0, 0, 0}; m1[[2]] = {1, 2, 3}; m1", "{0, {1, 2, 3}, 0}");
    /* Multi-position Part assignment still distributes element-wise. */
    assert_eval("m2 = {0, 0, 0}; m2[[{1, 3}]] = {9, 8}; m2", "{9, 0, 8}");
    assert_eval("m3 = {1, 2, 3, 4}; m3[[2 ;; 3]] = {20, 30}; m3", "{1, 20, 30, 4}");
}

/* ---------------------------------------------------- the integration method */
static void test_method_transcendental(void) {
    /* Qualified-symbol surface, cold (lazy load happens here). */
    assert_eval(
        "Simplify[D[Integrate`ParallelMixedTower[Log[x], x], x] - Log[x]]", "0");
    /* Method-option surface. */
    assert_eval(
        "Simplify[D[Integrate[x Log[x], x, Method -> \"ParallelMixedTower\"], x]"
        " - x Log[x]]", "0");
    assert_eval(
        "Simplify[D[Integrate`ParallelMixedTower[1/(x Log[x]), x], x]"
        " - 1/(x Log[x])]", "0");
}

static void test_method_radical(void) {
    /* A simple radical y^2 = x^2 + 1: the flagship non-transcendental case. */
    assert_eval(
        "Simplify[D[Integrate`ParallelMixedTower[1/Sqrt[x^2 + 1], x], x]"
        " - 1/Sqrt[x^2 + 1]]", "0");
    assert_eval(
        "Simplify[D[Integrate`ParallelMixedTower[x/Sqrt[x^2 + 1], x], x]"
        " - x/Sqrt[x^2 + 1]]", "0");
}

static void test_method_split_specials(void) {
    /* Sqrt[Tan[x]] flattens to Integrate[2 u^2/(1 + u^4), u] with u = Sqrt[Tan[x]];
     * the tower special 1 + u^4 must be split into its four linear factors over
     * the algebraic closure (SplitSpecials) before the ansatz is solvable.
     * Regression for the nested-Function closure bug that left `splittable`
     * False (Mathilda's Function does not close over an enclosing Function's
     * parameter), which disabled the split and reported "not elementary".
     * Branch-sensitive over the Root objects, so verified numerically -- and the
     * head is checked too, since D[unevaluated Integrate] returns the integrand
     * and would let a decline pass the residual test vacuously. */
    assert_eval(
        "r = Integrate[Sqrt[Tan[x]], x, Method -> \"ParallelMixedTower\"];"
        " {Head[r] =!= Integrate,"
        "  Abs[N[(D[r, x] - Sqrt[Tan[x]]) /. x -> 1/2, 25]] < 10^-15}",
        "{True, True}");
}

static void test_method_declines_cleanly(void) {
    /* Genuinely non-elementary: the Method form must DECLINE, leaving Integrate
     * unevaluated so the cascade / caller sees no answer (not a wrong one). */
    assert_eval(
        "Head[Integrate[Exp[x^2], x, Method -> \"ParallelMixedTower\"]]", "Integrate");
    /* The cascade must not spend this stage on a pure rational function -- the
     * gate skips it, and Integrate falls through to BronsteinRational. */
    assert_eval("Integrate[1/(1 + x^2), x]", "ArcTan[x]");
}

static void test_method_certifies_nonelementary(void) {
    /* 1/(x Log[x + Sqrt[x^2+1]]) is provably non-elementary: the residue at a
     * normal prime is Sqrt[1+x^2]/x, outside the constant field (paper Thm 9.2(a),
     * the "not elementary in one line" example).  The .m worker returns the raw
     * {"not elementary", ...} certificate on the qualified-symbol surface. */
    assert_eval(
        "Head[Integrate`ParallelMixedTower[1/(x Log[x + Sqrt[x^2 + 1]]), x]]",
        "List");
    assert_eval(
        "Integrate`ParallelMixedTower[1/(x Log[x + Sqrt[x^2 + 1]]), x][[1]]",
        "\"not elementary\"");   /* a String prints with its quotes */
    /* Through the Method surface the certificate is (as with RischTranscendental)
     * reported and the integral left unevaluated -- never a wrong answer. */
    assert_eval(
        "Head[Integrate[1/(x Log[x + Sqrt[x^2 + 1]]), x, "
        "Method -> \"ParallelMixedTower\"]]",
        "Integrate");
}

static void test_nonelem_warning_emitted(void) {
    /* A PROVED certificate must issue Integrate::nonelem, exactly as
     * RischTranscendental does for its own field decision. */
    ASSERT(eval_stderr_contains(
        "Integrate[1/(x Log[x + Sqrt[x^2 + 1]]), x, Method -> \"ParallelMixedTower\"]",
        "Integrate::nonelem"));
    /* An inconclusive {"failed", ...} give-up (here Exp[x^2], which this stage
     * cannot certify) must stay SILENT: failure of the parallel method proves
     * nothing about elementarity. */
    ASSERT(!eval_stderr_contains(
        "Integrate[Exp[x^2], x, Method -> \"ParallelMixedTower\"]",
        "Integrate::nonelem"));
}

/* --------------------------------------------- soundness/robustness regressions
 * Guards for the defects the sympy-comparison stress test surfaced (see
 * tasks/parallelmixedtower_stress.md).  Each was a wrong answer, a false
 * certificate, an internal leak, or a hang on the strict method surface. */
static void test_soundness_fixes(void) {
    /* Fix 1 (N6): Sqrt[Log[x]] is non-elementary; the worker MUST decline with a
     * {"failed", ...} tuple, not return 0.  Root cause was OptionValue[
     * "SpecialExponent"] leaking (Options[iPIM] undefined) -> symbolic bound ->
     * 0-equation trivial solve -> a bogus 0 antiderivative. */
    assert_eval(
        "MatchQ[Integrate`ParallelMixedTower[Sqrt[Log[x]], x], {\"failed\", ___}]",
        "True");
    /* Fix 1 (F2): the same leak baked an unevaluated OptionValue[...] into the
     * answer for Tan[Sqrt[x]]/Sqrt[x].  Now: no OptionValue in the result, and it
     * differentiates back to the integrand. */
    assert_eval(
        "r = Integrate`ParallelMixedTower[Tan[Sqrt[x]]/Sqrt[x], x];"
        " {FreeQ[r, OptionValue],"
        "  Abs[N[(D[r, x] - Tan[Sqrt[x]]/Sqrt[x]) /. x -> 1/2, 25]] < 10^-15}",
        "{True, True}");
    /* Fix 2 (T2): 1/(x (Log[x]^2+1)) was a FALSE non-elementary certificate
     * (its integral is ArcTan[Log[x]]).  Root cause: a ragged Transpose[{pts,
     * taus}] left unevaluated in RealisePoints, so unequal residues over the
     * log-tower prime were never realised.  Now it solves (not a List). */
    assert_eval(
        "r = Integrate`ParallelMixedTower[1/(x (Log[x]^2 + 1)), x];"
        " {Head[r] =!= List,"
        "  Abs[N[(D[r, x] - 1/(x (Log[x]^2 + 1))) /. x -> 2, 25]] < 10^-15}",
        "{True, True}");
    /* Fix 2 (T10): 1/(x^4-1) -- same false-certificate class over x^2+1. */
    assert_eval(
        "r = Integrate`ParallelMixedTower[1/(x^4 - 1), x];"
        " Abs[N[(D[r, x] - 1/(x^4 - 1)) /. x -> 2, 25]] < 10^-15",
        "True");
    /* Fix 3 (A11): x^3 ArcSin[x]/Sqrt[1-x^4] previously returned a WRONG
     * antiderivative (bad quartic-realisation branch).  The hard verify-or-decline
     * gate now differentiates the answer at generic complex points and demotes a
     * non-verifying result to a clean {"failed", ...} decline. */
    assert_eval(
        "MatchQ[Integrate`ParallelMixedTower[x^3 ArcSin[x]/Sqrt[1 - x^4], x], {\"failed\", ___}]",
        "True");
    /* Fix 5 (R3, R4): simplest genus-0 conic + one rational pole -- were declines,
     * fixed by the same Transpose->Thread correction on the curve branch. */
    assert_eval(
        "r = Integrate`ParallelMixedTower[1/((x + 1) Sqrt[x^2 + x + 1]), x];"
        " Abs[N[(D[r, x] - 1/((x + 1) Sqrt[x^2 + x + 1])) /. x -> 2, 25]] < 10^-15",
        "True");
    assert_eval(
        "r = Integrate`ParallelMixedTower[Sqrt[x^2 + 1]/x, x];"
        " Abs[N[(D[r, x] - Sqrt[x^2 + 1]/x) /. x -> 2, 25]] < 10^-15",
        "True");
}

void test_parallelmixedtower(void) {
    symtab_init();
    core_init();
    /* Production-like environment + chdir to the repo root so the lazy loader
     * resolves src/internal/mixed/ParallelMixed.m. */
    test_load_init_m();

    TEST(test_message_subsystem);
    TEST(test_core_assignment_fixes);
    TEST(test_method_transcendental);
    TEST(test_method_radical);
    TEST(test_method_split_specials);
    TEST(test_method_declines_cleanly);
    TEST(test_method_certifies_nonelementary);
    TEST(test_nonelem_warning_emitted);
    TEST(test_soundness_fixes);

    printf("All ParallelMixedTower tests passed!\n");
}

int main(void) {
    test_parallelmixedtower();
    return 0;
}
