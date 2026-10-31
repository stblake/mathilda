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
 * `needle` appears in what was written there.  Used to assert a diagnostic
 * PRINTS on a bare (non-Quiet) call -- the print is what Quiet[] suppresses and
 * a bare call still shows.  (Integrate::nonelem and the rest now route through
 * the mth_message funnel, so Quiet[] silences them and Check[] catches them; see
 * test_message_subsystem.)  dup/dup2 to a tmpfile (not freopen to /dev/tty) so
 * the restore works in CI. */
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

    /* --- The migrated diagnostics all support Quiet[] and Check[] now. Every
     * one of these was a raw fprintf(stderr, ...) that escaped both before the
     * message-routing migration; a representative head per class. --- */
    /* common.c builtin_arg_error (the argx/argt/... class shared by most heads). */
    assert_eval("Check[Fourier[], CAUGHT]", "CAUGHT");
    assert_eval("Check[IntegerDigits[1, 2, 3, 4], CAUGHT]", "CAUGHT");
    /* power.c / matpow.c -- previously Check-blind / escaped both. */
    assert_eval("Check[Power[0, -2], CAUGHT]", "CAUGHT");
    assert_eval("Check[MatrixPower[{{1, 2}, {3, 4}}, 1/2], CAUGHT]", "CAUGHT");
    /* linalg / eval / regex representatives. */
    assert_eval("Check[Det[{{1, 2, 3}, {4, 5, 6}}], CAUGHT]", "CAUGHT");
    assert_eval("Check[Sin = 5, CAUGHT]", "CAUGHT");   /* Set::wrsym (Protected) */

    /* Quiet[] silences the print at these sites (no message-driven failure). */
    assert_eval("Quiet[MatrixPower[{{1, 2}, {3, 4}}, 1/2]]",
                "MatrixPower[{{1, 2}, {3, 4}}, 1/2]");
    ASSERT(!eval_stderr_contains("Quiet[Det[{{1, 2, 3}, {4, 5, 6}}]]", "Det::"));

    /* CRUCIAL: an INTERNAL probe mute (g_arith_warnings_muted) is invisible to
     * Check, not merely silent -- Limit probes 0*ComplexInfinity while computing
     * Sin[x]/x -> 1, firing Infinity::indet internally, and Check must NOT catch
     * that sampling noise (it did, briefly, when the fix over-noted). */
    assert_eval("Check[Limit[Sin[x]/x, x -> 0], CAUGHT]", "1");
    assert_eval("Check[Limit[1/x, x -> 0], CAUGHT]", "ComplexInfinity");
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

/* The real form of the logarithmic part (src/internal/mixed/logrewrite.m, Rioboo):
 * conjugate pairs of logarithms come back as real logarithms, arctangents and
 * hyperbolic arctangents, so a real integrand gets an answer free of I and, where
 * the pairs cover the whole logarithmic part, free of Log as well. */
static void test_method_real_form(void) {
    /* a conjugate pair over a transcendental generator: ArcTan[Log[x]] */
    assert_eval("Integrate`ParallelMixedTower[1/(x (Log[x]^2 + 1)), x]", "ArcTan[Log[x]]");
    /* over the radical y^2 = x: the arctangent of the radical itself */
    assert_eval("Integrate`ParallelMixedTower[Sqrt[x]/(x + 1), x]", "2 Sqrt[x] - 2 ArcTan[Sqrt[x]]");
    /* a pair conjugate under Sqrt[2] -> -Sqrt[2]: a hyperbolic arctangent, no Log, no I */
    assert_eval(
        "r = Integrate`ParallelMixedTower[Sin[x]/(1 + Sin[x]^2), x];"
        " {FreeQ[r, Complex], FreeQ[r, Log], ! FreeQ[r, ArcTanh],"
        "  Abs[N[(D[r, x] - Sin[x]/(1 + Sin[x]^2)) /. x -> 2/3, 25]] < 10^-15}",
        "{True, True, True, True}");
    /* the arctangent argument is a polynomial in the radical (YQuot): the surface layer
     * rewrites the positive powers of the radical only, and a radical in a denominator
     * was read on the principal branch -- wrong wherever Tan[x] < 0 (Charlwood P5) */
    /* Since v0.6093 the answer must hold for ALL complex x.  The Tan[x] surface
     * reads Sqrt[Cos^4 + Cos^2 + 1] as Cos^2 Sqrt[3 + 3 Tan^2 + Tan^4], true only
     * where Cos[x]^2 >= 0 (the real line), so the method now declines; whatever
     * it returns must verify for all complex x (a declined List passes). */
    assert_eval(
        "f = Cos[x]^2/Sqrt[Cos[x]^4 + Cos[x]^2 + 1]; r = Quiet[Integrate`ParallelMixedTower[f, x]];"
        " If[ListQ[r], First[r] =!= \"not elementary\","
        "   FreeQ[r, Complex] && TrueQ[PossibleZeroQ[D[r, x] - f, Assumptions -> Element[x, Complexes]]]]",
        "True");
    /* two nested-radical fields in one session (Charlwood P4 then A40): the module
     * canonicalises with CanRaw, not with the field detour of Can, whose back-conversion
     * once returned a Dot[{}, Inverse[{}], {}] coefficient here */
    assert_eval(
        "Integrate`ParallelMixedTower[Log[x Sqrt[x^2 + 1] + 1], x];"
        " f = ArcTan[x Sqrt[1 - x^2]]; r = Integrate`ParallelMixedTower[f, x];"
        " {FreeQ[r, Complex | Dot | Inverse], Abs[N[(D[r, x] - f) /. x -> 1/2, 25]] < 10^-15}",
        "{True, True}");
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

/* False certificates through Series: at an order-9 (resp. 7) pole the deep
 * residues came back {0, 0, 0, 0} from a Series that stopped short of the
 * e^-1 term, and the holomorphic-remainder argument then "proved" these
 * elementary integrands non-elementary.  (1-x^2) Sqrt[1+x^2] on the reals. */
static void test_no_false_certificate_series_precision(void) {
    assert_eval("Quiet[ParallelMixed`ParallelIntegrateMixed[Sqrt[1 - x] Sqrt[x + 1] Sqrt[1 - x^4], x]][[1]] =!= \"not elementary\"",
                "True");
    assert_eval("Quiet[ParallelMixed`ParallelIntegrateMixed[Sqrt[x + 1] Sqrt[1 - x^4]/Sqrt[1 - x], x]][[1]] =!= \"not elementary\"",
                "True");
}

/* Run LAST: it replaces a package helper.  A $RecursionLimit overflow does not
 * abort -- the evaluator Holds the call and degrades every enclosing loop to one
 * rewrite -- so nothing produced while one fired may be returned as an answer
 * or a certificate.  Force a real overflow inside the 2-arg entry and require a
 * failure tuple from the package, the qualified symbol and the Method form. */
/* A root nested on a conic: BuildTower parametrises the conic BEFORE it
 * normalises the outer root.  Splitting a generator factor off the outer base
 * first (Sqrt[x (x + y)] -> Sqrt[x] Sqrt[x + y]) flattened x = u^2 and made the
 * conic a quartic in u, so the outer root demanded a second radical and the
 * tower was never built.  The Euler step also takes a parametric square root of
 * the leading (or constant) coefficient: Sqrt[a^2] -> a, Sqrt[b^2] -> b. */
static void test_euler_before_root_split(void) {
    const char *cases[] = {
        "Sqrt[x^2 + x Sqrt[x^2 - 1]]/(x Sqrt[x^2 - 1])",      /* DDI potpourri p8 (v) */
        "Sqrt[Sqrt[a^2 x^2 + b] + a x]",                      /* nested 15 (i): Sqrt[a^2] */
        "Sqrt[Sqrt[a x^2 + b^2] + b]",                        /* nested 14 (i): point (0, b) */
        "Sqrt[a x + Sqrt[a^2 x^2 + b x + c]]/((b x + c) Sqrt[a^2 x^2 + b x + c])",
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Integrate[f, x, Method -> \"ParallelMixedTower\"]}, "
                 "FreeQ[F, Integrate] && TrueQ[PossibleZeroQ[D[F, x] - f, "
                 "Assumptions -> Element[x, Complexes]]]]]", cases[i]);
        assert_eval(buf, "True");
    }
}

/* The u = x^k pre-pass integrates g(u) in a Unique["u"], which is Temporary;
 * the acceptance gate rejected every answer holding a Temporary symbol (an
 * escaped Module local), the integration variable included, so the sub-answer
 * was always dropped and the full curve y^2 = 14 x^4 + 31 x^2 + 15 (genus 1)
 * declined on bounds. */
static void test_xpower_subcall_variable_accepted(void) {
    assert_eval("With[{f = x Sqrt[15 + 31 x^2 + 14 x^4]/(5 + 7 x^2)}, "
                "With[{F = Integrate[f, x, Method -> \"ParallelMixedTower\"]}, "
                "FreeQ[F, Integrate] && TrueQ[PossibleZeroQ[D[F, x] - f, "
                "Assumptions -> Element[x, Complexes]]]]]", "True");
}

/* Several flat radicands (a radical compositum) that BuildTower cannot fold
 * into one radical: the integrand is split into Kummer pieces g M, M a radical
 * monomial, one per class modulo C(x)^*, and each piece is integrated through
 * one radical P (M = kappa P, kappa a locally constant root of unity). */
static void test_kummer_split_flat_radicals(void) {
    const char *cases[] = {
        "Sqrt[2 x^2 + 3]/Sqrt[5/x^2 + 7]",        /* DDI chebychev 5 (i): one class */
        "x^(3/2) Sqrt[x + 1] + Sqrt[x]/x",         /* two classes */
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Integrate[f, x, Method -> \"ParallelMixedTower\"]}, "
                 "FreeQ[F, Integrate] && TrueQ[PossibleZeroQ[D[F, x] - f, "
                 "Assumptions -> Element[x, Complexes]]]]]", cases[i]);
        assert_eval(buf, "True");
    }
    /* a piece proven non-elementary certifies the whole integrand: here the one
     * piece is x^2 (x^6 - 1)^(-2/3), i.e. (1/3) Int (u^2 - 1)^(-2/3) du, which
     * fails all three Chebyshev conditions */
    assert_eval("First[Quiet[ParallelMixed`ParallelIntegrateMixed["
                "(x^3 - 1)^(1/3) (x^3 + 1)^(1/3) x^2/(x^6 - 1), x]]]", "\"not elementary\"");
    /* two different monomials in one class (they differ by a sign that changes
     * across the plane) are not split */
    assert_eval("ParallelMixed`Private`KummerSplit[x/(Sqrt[x - 1] Sqrt[x + 1]) + Sqrt[x^2 - 1], x]",
                "None");
}

/* DDI multiple_radicals_0 (i): the answer through the Euler parameter of
 * Sqrt[x^2 - a] was -f on parts of x^2 < a when a, b, c, d > 0 -- a real-axis
 * cut in 2 sign orthants of 32 -- and PossibleZeroQ over complex x missed it
 * (about one real draw per orthant).  The complex-domain screen now visits
 * each orthant 16 times, half of them real: the gate declines, or answers
 * with an antiderivative right at those points. */
static void test_complex_screen_catches_real_slice_cut(void) {
    assert_eval("With[{f = (x^4 - a b)/((d x^4 - (a d + b d - c) x^2 + a b d) Sqrt[x^2 - a] Sqrt[x^2 - b])}, "
                "With[{F = Integrate[f, x, Method -> \"ParallelMixedTower\"]}, "
                "!FreeQ[F, Integrate] || AllTrue[{1/5, -1/2, -1}, "
                "Abs[N[(D[F, x] - f) /. {a -> 3/2, b -> 5/3, c -> 7/4, d -> 9/5} /. x -> #, 30]] < 10^-15 &]]]",
                "True");
}

/* PolynomialLCM[1, b/a + x^3] stayed unevaluated, and the linear solve over
 * Q(a, b) read it as a nonzero residual: five DDI cube-root integrals with
 * parameters declined "internal: solution does not verify". */
static void test_parametric_lcm_in_linear_solve(void) {
    const char *cases[] = {
        "(a + b x Sqrt[x])^(1/3)",               /* DDI nested 3 (ii) */
        "1/(a x^3 + b)^(1/3)",                   /* DDI chebychev 3 (v) */
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
                 "With[{f = %s}, With[{F = Integrate[f, x, Method -> \"ParallelMixedTower\"]}, "
                 "FreeQ[F, Integrate] && TrueQ[PossibleZeroQ[D[F, x] - f, "
                 "Assumptions -> Element[x, Complexes]]]]]", cases[i]);
        assert_eval(buf, "True");
    }
}

static void test_recursion_overflow_declines(void) {
    assert_eval("ParallelMixed`Private`boom[n_] := 1 + ParallelMixed`Private`boom[n + 1]; "
                "ClearAll[ParallelMixed`Private`VanishingDenominator]; "
                "ParallelMixed`Private`VanishingDenominator[___] := ParallelMixed`Private`boom[0]; 0", "0");
    /* the direct package call: a failure tuple (the reason depends on where the
     * degraded evaluation stops -- since v0.6095 the u = x^2 sub-call overflows
     * first and the 45 s backstop ends the run); never an answer */
    assert_eval("First[Quiet[ParallelMixed`ParallelIntegrateMixed[Sqrt[x^2 + 1]/x, x]]]", "\"failed\"");
    assert_eval("Quiet[Integrate`ParallelMixedTower[Sqrt[x^2 + 1]/x, x]]",
                "{\"failed\", \"internal error: $RecursionLimit exceeded\"}");
    assert_eval("Head[Quiet[Integrate[Sqrt[x^2 + 1]/x, x, Method -> \"ParallelMixedTower\"]]]",
                "Integrate");
}

static void test_method_declines_cleanly(void) {
    /* Genuinely non-elementary: the Method form must DECLINE, leaving Integrate
     * unevaluated so the cascade / caller sees no answer (not a wrong one). */
    assert_eval(
        "Head[Integrate[Exp[x^2], x, Method -> \"ParallelMixedTower\"]]", "Integrate");
    /* The cascade must not spend this stage on a pure rational function -- the
     * gate skips it, and Integrate falls through to BronsteinRational. */
    assert_eval("Integrate[1/(1 + x^2), x]", "ArcTan[x]");
    /* A second radical over a conic y^2 = 2 g^2 + 3 with an irrational leading
     * square root AND no rational point.  The Euler step's decline sat in a
     * Return inside a With, which exits only the With: the step went on with its
     * derivation unassigned, overflowed $RecursionLimit, and the method returned
     * the unassigned local ParallelMixed`Private`r$NNN as the antiderivative --
     * also through the default Integrate cascade (56 cases of the DDI corpus). */
    /* Since v0.6099 the Kummer split integrates it (one class, through
     * x Sqrt[(2 x^2 + 3) (7 x^2 + 5)]/(7 x^2 + 5)); the guard is that the
     * answer is a verified antiderivative, never a leaked local. */
    assert_eval(
        "With[{f = Sqrt[2 x^2 + 3]/Sqrt[5/x^2 + 7]}, With[{F = Quiet[ParallelMixed`ParallelIntegrateMixed[f, x]]}, "
        "!ListQ[F] && FreeQ[F, _Symbol?(MemberQ[Attributes[#], Temporary] &)] && "
        "TrueQ[PossibleZeroQ[D[F, x] - f, Assumptions -> Element[x, Complexes]]]]]",
        "True");
    assert_eval(
        "FreeQ[Quiet[Integrate[Sqrt[2 x^2 + 3]/Sqrt[5/x^2 + 7], x, Method -> \"ParallelMixedTower\"]], "
        "_Symbol?(MemberQ[Attributes[#], Temporary] &) | Integrate]",
        "True");
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
    /* Exp[x^2] now CERTIFIES and so must WARN.  Until v0.240 it was this test's
     * "stays silent" exemplar, because the holomorphic-remainder certificate was
     * hard-gated to a tower carrying a curve; on a curve-free monomial tower the
     * bounds {2,1} are proved, the residual is residue-free and the system is
     * inconsistent, which is Proposition 9.2(b) itself. */
    ASSERT(eval_stderr_contains(
        "Integrate[Exp[x^2], x, Method -> \"ParallelMixedTower\"]",
        "Integrate::nonelem"));
    /* An inconclusive {"failed", ...} give-up must still stay SILENT: failure of
     * the parallel method proves nothing about elementarity.  Sqrt[Log[x]]
     * declines with {"failed", "no solution within bounds", {4, 3}} at every
     * rung of the retry ladder, so the certificate is never reached. */
    ASSERT(!eval_stderr_contains(
        "Integrate[Sqrt[Log[x]], x, Method -> \"ParallelMixedTower\"]",
        "Integrate::nonelem"));
    /* Integrate::nonelem now routes through the funnel: Quiet[] suppresses the
     * print (it was a deliberate raw fprintf before this migration). */
    ASSERT(!eval_stderr_contains(
        "Quiet[Integrate[1/(x Log[x + Sqrt[x^2 + 1]]), x, "
        "Method -> \"ParallelMixedTower\"]]",
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

/* ------------------------------ no false certificate on a curve-free tower
 * The certificate of Proposition 9.2(b) is no longer hard-gated to a tower
 * carrying a curve (that gate was belt-and-braces over the T2/T10 residue-
 * realisation bug, fixed at its root in v0.168/v0.175).  Every integrand here
 * is ELEMENTARY and lives on a curve-free tower, so a {"not elementary", ...}
 * from any of them is a FALSE certificate -- the one failure mode worse than
 * declining.  Asserted as "solves", which is stronger than "does not certify":
 * the certificate branch is entered only when the ansatz system is
 * inconsistent, so a solution proves the branch was never reached. */
static void test_no_false_certificate_curve_free(void) {
    static const struct { const char* f; const char* pt; } cases[] = {
        /* the two historical false certificates the retired guard was added for */
        { "1/(x (Log[x]^2 + 1))",  "2"   },
        { "1/(x^4 - 1)",           "2"   },
        { "1/(x^4 + 1)",           "2"   },
        /* reach sol === {} with ALL BOUNDS PROVED and a residue-free residual in
         * the base run; only the split of the degree-4 special 1 + u^4 over Fbar
         * saves them, so they are the tightest guard on the retired conjunct */
        { "Sqrt[Tan[x]]",          "1/2" },
        { "Tan[x]^(3/2)",          "1/2" },
        { "1/Sqrt[Tan[x]]",        "1/2" },
        /* the elementary twins of the newly-certifying family, on the identical
         * towers with the identical specials */
        { "(1 + 2 x^2) Exp[x^2]",  "3/5" },
        { "(1 - 2 x^2) Exp[-x^2]", "3/5" },
        { "(1 - 3 x^3) Exp[-x^3]", "3/5" },
        { "(x - 1) Exp[x]/x^2",    "7/5" },
        { "Exp[1/x]/x^2",          "7/5" },
        /* algebraically DEPENDENT generators (E^x with E^(2x), Log[x] with
         * Log[x^2]).  Every residue argument over such a tower is vacuous, and
         * all four certified falsely until "StructureTheorem" became the default
         * for ParallelIntegrateMixed.  They are the reason that default moved. */
        { "Exp[2 x]/(1 + Exp[x])",        "7/5" },
        { "Exp[x]/(1 + Exp[x] + Exp[2 x])", "7/5" },
        { "1/(Exp[x] - Exp[-x])",         "7/5" },
    };
    /* Log[x^2]/Log[x] is 2 only for Re x > 0: 2 x is not an antiderivative for
     * all complex x, so since v0.6093 (PossibleZeroQ gate) the method DECLINES --
     * the point of this test is that it must never CERTIFY non-elementarity. */
    assert_eval("First[Quiet[Integrate`ParallelMixedTower[Log[x^2]/Log[x], x]]] =!= \"not elementary\"", "True");
    for (size_t i = 0; i < sizeof cases / sizeof *cases; i++) {
        char buf[768];
        snprintf(buf, sizeof buf,
                 "Module[{r = Integrate`ParallelMixedTower[%s, x]},"
                 " Head[r] =!= List && "
                 " Abs[N[(D[r, x] - (%s)) /. x -> %s, 30]] < 10^-20]",
                 cases[i].f, cases[i].f, cases[i].pt);
        Expr* parsed = parse_expression(buf);
        ASSERT(parsed != NULL);
        Expr* evaluated = evaluate(parsed);
        expr_free(parsed);
        char* s = expr_to_string(evaluated);
        ASSERT_MSG(strcmp(s, "True") == 0,
                   "curve-free elementary integrand did not solve (a List here is "
                   "a FALSE non-elementarity certificate): %s -> %s",
                   cases[i].f, s);
        free(s);
        expr_free(evaluated);
    }
}

/* ------------- the certificate that the retired guard used to suppress
 * Genuinely non-elementary, curve-free: proved bounds, a complete logand set
 * and a verified residue-free residual, which is Proposition 9.2(b) whether or
 * not the tower carries a curve.  Worth ~50 cases of the 312-case corpus. */
static void test_certifies_curve_free_nonelementary(void) {
    static const char* nonelem[] = {
        "Exp[-x^2]", "Exp[x^2]", "Exp[-x^3]", "Sqrt[x] Exp[x]", "Exp[1/x]",
    };
    for (size_t i = 0; i < sizeof nonelem / sizeof *nonelem; i++) {
        char buf[512];
        snprintf(buf, sizeof buf,
                 "Module[{r = Integrate`ParallelMixedTower[%s, x]},"
                 " ListQ[r] && r[[1]] === \"not elementary\"]", nonelem[i]);
        Expr* parsed = parse_expression(buf);
        ASSERT(parsed != NULL);
        Expr* evaluated = evaluate(parsed);
        expr_free(parsed);
        char* s = expr_to_string(evaluated);
        ASSERT_MSG(strcmp(s, "True") == 0,
                   "curve-free non-elementary integrand was not certified: %s -> %s",
                   nonelem[i], s);
        free(s);
        expr_free(evaluated);
    }
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
    TEST(test_method_real_form);
    TEST(test_method_split_specials);
    TEST(test_method_declines_cleanly);
    TEST(test_no_false_certificate_series_precision);
    TEST(test_method_certifies_nonelementary);
    TEST(test_nonelem_warning_emitted);
    TEST(test_soundness_fixes);
    TEST(test_no_false_certificate_curve_free);
    TEST(test_certifies_curve_free_nonelementary);
    TEST(test_nontorsion_divisor_certificate);
    TEST(test_euler_before_root_split);
    TEST(test_xpower_subcall_variable_accepted);
    TEST(test_kummer_split_flat_radicals);
    TEST(test_complex_screen_catches_real_slice_cut);
    TEST(test_parametric_lcm_in_linear_solve);
    TEST(test_recursion_overflow_declines);   /* last: redefines a package helper */

    printf("All ParallelMixedTower tests passed!\n");
}

int main(void) {
    test_parallelmixedtower();
    return 0;
}
