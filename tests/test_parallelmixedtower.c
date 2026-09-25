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
    /* Nested radical Sqrt[x + Sqrt[x]] = Sqrt[Sqrt[x] (1 + Sqrt[x])]: the outer
     * radicand is a genus-0 conic in u = Sqrt[x], so it must stay FUSED as the
     * simple radical y^2 = u^2 + u rather than split (branch-unsafe) into the
     * x^(1/4) Sqrt[1 + Sqrt[x]] form that is a correct antiderivative only on
     * x > 0.  Guards both properties: the answer is free of the split fourth root
     * x^(1/4), and D[r] equals the integrand at a COMPLEX point off the positive
     * reals -- the global-fidelity check the real-only verify gate cannot make,
     * where the old split form differed from the integrand. */
    assert_eval(
        "r = Integrate[Sqrt[x + Sqrt[x]], x, Method -> \"ParallelMixedTower\"];"
        " {Head[r] =!= Integrate, FreeQ[r, Power[x, 1/4]],"
        "  Abs[N[(D[r, x] - Sqrt[x + Sqrt[x]]) /. x -> 2 + I, 25]] < 10^-15}",
        "{True, True, True}");
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
    /* Fix 3 (A11): x^3 ArcSin[x]/Sqrt[1-x^4] used to return a WRONG antiderivative
     * (the y-branch of the quartic realisation was the negative of the one the
     * integrand pair was decomposed on), so the gate demoted it to a decline.  The
     * gate now RESOLVES that branch numerically (accepts -surf when D[surf] == -f
     * on the real domain), so A11 solves and verifies. */
    assert_eval(
        "r = Integrate`ParallelMixedTower[x^3 ArcSin[x]/Sqrt[1 - x^4], x];"
        " Abs[N[(D[r, x] - x^3 ArcSin[x]/Sqrt[1 - x^4]) /. x -> 1/3, 25]] < 10^-15",
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

/* ------------------- non-torsion certificate at finite places (2026-09-22) */
static void test_nontorsion_divisor_certificate(void) {
    /* The hyperbolic parity rule (integrands odd in Sinh: u = Cosh[x] with
     * sinh^2 = cosh^2 - 1).  Before the fix the tower used sinh^2 = 1 + cosh^2
     * and Sinh[x]^3 came back as Cosh[x] + Cosh[x]^3/3. */
    assert_eval(
        "Simplify[D[Integrate`ParallelMixedTower[Sinh[x]^3, x], x] - Sinh[x]^3]",
        "0");
    assert_eval(
        "r = Integrate`ParallelMixedTower[Coth[x] Sqrt[Cosh[x]], x];"
        " Abs[N[(D[r, x] - Coth[x] Sqrt[Cosh[x]]) /. x -> 3/2, 25]] < 10^-15",
        "True");
    /* Proposition 9.4 for a residue divisor at finite places (Part II, example
     * 10.20): the class of (1, Sqrt[2]) - (1, -Sqrt[2]) on y^2 = x^5 + 1 has
     * orders 25 and 29 in the Jacobians over GF(7) and GF(17) -- incompatible
     * with any finite order, so the divisor is not torsion.  The coordinates
     * are given over Q on the power basis of theta = Sqrt[2] (minimal
     * polynomial -2 + z^2, lowest first).  The package is loaded by the
     * earlier tests; the first call below makes that explicit. */
    assert_eval(
        "Integrate`ParallelMixedTower[Log[x], x];"
        " ParallelMixed`Private`NontorsionDivisorCertificate[x^5 + 1, x,"
        "  {{{{1}, {0, 1}}, 1}, {{{1}, {0, -1}}, -1}}, {-2, 0, 1}]",
        "{True, {{7, 25}, {17, 29}}}");
    /* A torsion class must NOT be certified: (0, 1) - (0, -1) on y^2 = x^5 + 1
     * is 5-torsion (div(y - 1) = 5 (0,1) - 5 oo), so every prime reports 5. */
    assert_eval(
        "ParallelMixed`Private`NontorsionDivisorCertificate[x^5 + 1, x,"
        "  {{{{0}, {1}}, 1}, {{{0}, {-1}}, -1}}, None, 4]",
        "{False, {{3, 5}, {7, 5}, {11, 5}, {13, 5}}}");
    /* The certificate for [oo+ - oo-] (Cohen's -72 variant) was dead in
     * Mathilda: its prime loop never ran (MATHILDA_DIVERGENCES.md A1) and the
     * continued fraction over GF(p) never terminated (A3).  Orders 3, 13, 7, 21
     * modulo 7, 11, 13, 17 are pairwise incompatible. */
    assert_eval(
        "ParallelMixed`Private`NontorsionCertificate[x^4 + 10 x^2 - 96 x - 72, x]",
        "{True, {{7, 3}, {11, 13}, {13, 7}, {17, 21}}}");
    /* Coth[x]/(1 + Sech[x]^5)^(3/2) itself: u = Sqrt[Cosh[x]] over y^2 = 1 + u^10
     * (genus 4); the certificate needs ~7000 Jacobian operations at p = 17 and
     * 41, which Mathematica does in 3 s and Mathilda's expression-level GF(p)
     * arithmetic does not fit into $ParallelMixedTimeBudget (MATHILDA_DIVERGENCES.md
     * C).  Until the core primitives are fixed the honest outcome is the clean
     * decline; the certificate is the other admissible answer.  TIGHTEN to the
     * certificate alone once the Modulus primitives are native. */
    assert_eval(
        "MatchQ[Integrate`ParallelMixedTower[Coth[x]/(1 + Sech[x]^5)^(3/2), x],"
        " {\"not elementary\", \"residue divisor not torsion: reduction mod p\","
        "  {{17, 29}, {41, 155}}, _} | {\"failed\", \"time budget exceeded\"}]",
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
    TEST(test_nontorsion_divisor_certificate);

    printf("All ParallelMixedTower tests passed!\n");
}

int main(void) {
    test_parallelmixedtower();
    return 0;
}
