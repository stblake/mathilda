/*
 * test_refine.c — Refine[expr, assum].
 *
 * Covers every example from the Refine specification: assumption-driven
 * rewrites (Sqrt/Power, Log, trig, Abs/Sign/Re, Floor/Ceiling/FractionalPart/
 * Mod), domain and predicate decisions (Element, equations, inequalities via
 * the Reduce/CAD entailment), the option/scoping semantics shared with
 * Simplify/PossibleZeroQ, and the argument-count errors.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "message.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static Expr* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* e = evaluate(parsed);
    expr_free(parsed);
    return e;
}

/* Assert eval(input) prints exactly as `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = eval_str(input);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

static void test_power_sqrt(void) {
    check("Refine[Sqrt[x^2], x > 0]", "x");
    check("Refine[Sqrt[x^2], Element[x, Reals]]", "Abs[x]");
    check("Refine[Sqrt[x^2], x < 0]", "-x");
    check("Refine[Sqrt[x^2 y^2], x > 0 && y < 0]", "-x y");
    check("Refine[(x^3)^(1/3), x >= 0]", "x");
    check("Refine[(a^b)^c, -1 < b < 1]", "a^(b c)");
    check("Refine[a^p b^p, a > 0 && b > 0]", "(a b)^p");
}

static void test_log(void) {
    check("Refine[Log[x], x < 0]", "I Pi + Log[-x]");
    check("Refine[Log[x^p], x > 0 && Element[p, Reals]]", "p Log[x]");
}

static void test_trig(void) {
    check("Refine[Sin[k Pi], Element[k, Integers]]", "0");
    check("Refine[Cos[x + k Pi], Element[k, Integers]]", "(-1)^k Cos[x]");
    check("Refine[ArcTan[Tan[x]], -Pi/2 < Re[x] < Pi/2]", "x");
}

static void test_sign_abs_re(void) {
    check("Refine[Sign[x^2 - x y + y^2 + 1], Element[x | y, Reals]]", "1");
    check("Refine[Re[a + b I], Element[a | b, Reals]]", "a");
    check("Refine[Abs[x], x > 0]", "x");
    check("Refine[Conjugate[x], Element[x, Reals]]", "x");
    /* Arg of a provably-signed quantity (prov_pos descends Times/Power, so a^2,
     * 2a, 1/a^2 count as positive under a > 0); i*(positive) -> Pi/2. */
    check("Refine[Arg[a], a > 0]", "0");
    check("Refine[Arg[a^2], a > 0]", "0");
    check("Refine[Arg[2 a], a > 0]", "0");
    check("Refine[Arg[I a], a > 0]", "1/2 Pi");
    check("Refine[Arg[2 I a], a > 0]", "1/2 Pi");
}

static void test_rounding_mod(void) {
    check("Refine[Floor[2 a + 1], Element[a, Integers]]", "1 + 2 a");
    check("Refine[Ceiling[x], 2 < x <= 3]", "3");
    check("Refine[FractionalPart[a], a < 0 && Mod[a, 1] == 1/3]", "-2/3");
    check("Refine[Mod[a, 4], Element[(a + 3)/4, Integers]]", "1");
}

static void test_element(void) {
    /* Direct fact. */
    check("Refine[Element[k, Reals], Element[k, Integers]]", "True");
    /* Compound-expression domain inference. */
    check("Refine[Element[(2 x + x^p)/(x Gamma[x + 2]), Reals], x > 0 && p > 0]", "True");
    check("Refine[Element[2 k^3 Floor[x]^k, Integers], "
          "Element[k, Integers] && k > 0 && Element[x, Reals]]", "True");
    /* Algebraic-in-inequality implies real. */
    check("Refine[Element[x, Reals], x^2 < 1]", "True");
}

static void test_predicates(void) {
    check("Refine[a^2 - b^2 + 1 == 0, a + b == 0]", "False");
    check("Refine[a^2 - a b + b^2 >= 0, Element[a | b, Reals]]", "True");
    check("Refine[(x - 1)^2 + (y - 2)^2 < 3/2, x^2 + y^2 <= 1]", "False");
    check("Refine[-1 < x < 1, x^2 < 1]", "True");
}

static void test_options_scoping(void) {
    /* Assumptions given both positionally and as an option. */
    check("Refine[Cos[k Pi]^m, Element[k, Integers], Assumptions -> Mod[m, 2] == 0]", "1");
    /* Assuming propagates via $Assumptions. */
    check("Assuming[x > 0, Refine[Sqrt[x^2]]]", "x");
    /* Positional assumption is combined with $Assumptions. */
    check("Assuming[x > 0, Refine[Sqrt[x^2 y^2], y < 0]]", "-x y");
    /* Assumptions -> option prevents Refine from using $Assumptions. */
    check("Assuming[x > 0, Refine[Sqrt[x^2 y^2], Assumptions -> y < 0]]", "-Sqrt[x^2] y");
    /* No assumptions -> identity. */
    check("Refine[Sqrt[x^2]]", "Sqrt[x^2]");
    check("Refine[Sqrt[x^2], True]", "Sqrt[x^2]");
}

static void test_errors_edges(void) {
    /* 0 positional args -> message + unevaluated. */
    check("Refine[]", "Refine[]");
    /* > 2 positional args -> message + unevaluated. */
    check("Refine[a, b, c]", "Refine[a, b, c]");
    /* A TimeConstraint option does not count as a positional assumption. */
    check("Refine[Sqrt[x^2], x > 0, TimeConstraint -> 1]", "x");
    /* Options exist. */
    check("MemberQ[Options[Refine][[All, 1]], Assumptions]", "True");
    check("MemberQ[Options[Refine][[All, 1]], TimeConstraint]", "True");
}

/* Regressions from the first-principles stress corpus (v0.203). Each guards a
 * gap the corpus surfaced; kept in SameQ form so printer-form drift never
 * breaks them. */
static void test_stress_fixes(void) {
    /* Equal predicate under an equality fact: the assumption-aware zero test's
     * Schwartz-Zippel sampler used to ignore coupling equalities and wrongly
     * return False. Now decided by equality substitution; still False under an
     * inequality fact (which does NOT force a==b). */
    check("Refine[a == b, a - b == 0]", "True");
    check("Refine[a == b, a == b]", "True");
    check("Refine[a^2 == b^2, a == b]", "True");
    check("Refine[a == b, a > b]", "False");
    /* PossibleZeroQ shares the fixed zero test. */
    check("PossibleZeroQ[a - b, Assumptions -> a - b == 0]", "True");
    check("PossibleZeroQ[a, Assumptions -> a > 0]", "False");   /* sign FALSE preserved */

    /* Sign domains as queried membership. */
    check("Refine[Element[x, Positive], x > 0]", "True");
    check("Refine[Element[x, Negative], x < 0]", "True");
    check("Refine[Element[x, NonNegative], x >= 0]", "True");
    check("Refine[Element[k^2, NonNegative], Element[k, Reals]]", "True");

    /* NonPositive rewrites (x <= 0 => -x), not just strict x < 0. */
    check("Refine[Sqrt[x^2], x <= 0]", "-x");
    check("Refine[Abs[x], x <= 0]", "-x");

    /* Arg under a sign fact. */
    check("Refine[Arg[x], x > 0]", "0");
    check("Refine[Arg[x], x < 0]", "Pi");

    /* Abs of a complex expression with real parts -> magnitude. */
    check("Refine[Abs[a + b I], Element[a | b, Reals]] === Sqrt[a^2 + b^2]", "True");

    /* Log identities under reality / positivity. */
    check("Refine[Log[x^2], Element[x, Reals]] === 2 Log[Abs[x]]", "True");
    check("Refine[Log[E^x], Element[x, Reals]]", "x");
    check("Refine[Log[a b], a > 0 && b > 0] === Log[a] + Log[b]", "True");

    /* Many assumed symbols but few used: the per-symbol rule synthesis no
     * longer overflows its buffer and drops every rule. */
    check("Refine[Abs[a1], a1 > 0 && a2 > 0 && a3 > 0 && a4 > 0 && a5 > 0 && "
          "a6 > 0 && a7 > 0 && a8 > 0 && a9 > 0 && a10 > 0 && a11 > 0 && a12 > 0 && "
          "a13 > 0 && a14 > 0 && a15 > 0 && a16 > 0 && a17 > 0 && a18 > 0 && "
          "a19 > 0 && a20 > 0]", "a1");
}


/* Evaluate `input`, assert its printed form is `expected` AND that no message
 * fired while it ran -- the direct guard against an internal probe's
 * diagnostic reaching an enclosing Check[] (which would then return its
 * failure branch for a successful computation). */
static void check_quiet(const char* input, const char* expected) {
    unsigned long fired0 = mth_msg_fired_count();
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n", input, expected, s);
    ASSERT_STR_EQ(s, expected);
    if (mth_msg_fired_count() != fired0)
        fprintf(stderr, "FAIL: %s fired a message\n", input);
    ASSERT(mth_msg_fired_count() == fired0);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* Re / Im of provably-real terms: Re[Log[b]] is Log[b] under b > 0 (not the
 * sign-agnostic ComplexExpand form Log[b^2]/2); Re distributes over Plus. */
static void test_refine_re_of_real_terms(void) {
    check_quiet("Refine[Re[Log[b]], b > 0]", "Log[b]");
    check_quiet("Refine[Im[Log[b]], b > 0]", "0");
    check_quiet("Refine[Re[-I Pi - Log[a] + Log[b]], a > 0 && b > 0]", "-Log[a] + Log[b]");
    check_quiet("Refine[Im[-I Pi - Log[a] + Log[b]], a > 0 && b > 0]", "-Pi");
    /* No reality proof: left to ComplexExpand / unevaluated as before. */
    check_quiet("Refine[Re[Log[b]], b < 0]", "1/2 Log[b^2]");
}

/* A DirectedInfinity whose direction has a decidable sign under the facts is
 * the canonical +-Infinity (or DirectedInfinity[+-I]); an undecidable one is
 * left alone.  DirectedInfinity[1/(-1 + s)] under s > 1 is the Limit answer
 * this canonicalises. */
static void test_directed_infinity(void) {
    check("Refine[DirectedInfinity[1/(s - 1)], s > 1]", "Infinity");
    check("Refine[DirectedInfinity[s - 1], s > 1]", "Infinity");
    check("Refine[DirectedInfinity[1 - s], s > 1]", "-Infinity");
    check("Refine[DirectedInfinity[a], a < 0]", "-Infinity");
    check("Refine[DirectedInfinity[a], b > 0]", "DirectedInfinity[a]");
    check("Refine[DirectedInfinity[s - 1], s > 0]", "DirectedInfinity[-1 + s]");
    check("Refine[{x, DirectedInfinity[c]}, c > 0]", "{x, Infinity}");
}

/* The assumption rules are built as text with the user's symbols pasted in;
 * a user symbol named like one of their pattern variables (r, p, rest, base,
 * c, m, k) used to be captured by it: Refine[1/Sqrt[r^2], r > 0] was -2 and
 * Refine[Log[p^3], p > 0] was 3 Log[3]. */
static void test_pattern_name_capture(void) {
    check("Refine[1/Sqrt[r^2], r > 0]", "1/r");
    check("Refine[Log[p^3], p > 0]", "3 Log[p]");
    check("Refine[Sqrt[rest^2 y^2], rest > 0]", "rest Sqrt[y^2]");
    check("Refine[Sqrt[base^2], base > 0]", "base");
    check("Refine[(-1)^(2 m), Element[m, Integers]]", "1");
    check("Simplify[c^2/Sqrt[c^2], c > 0]", "c");
    check("Refine[Sqrt[r^2 rest], r > 0]", "r Sqrt[rest]");
}

/* The Reduce-backed sign proof is vacuous for a non-real argument (over the
 * reals both p > 0 and p <= 0 are unsatisfiable): Arg[-a + I a] was 0 under
 * a > 0, which tripled the PV of 1/(x^3 - a^3).  Arg of a positive multiple of a
 * numeric constant is that constant's Arg. */
static void test_refine_arg_complex(void) {
    check("Refine[Arg[-a + I a], a > 0]", "3/4 Pi");
    check("Refine[Arg[1/2 (-a - Sqrt[-3 a^2])], a > 0]", "-2/3 Pi");
    check("Refine[Arg[a + I], a > 0]", "Arg[I + a]");
    check("Refine[Abs[-a + I a], a > 0]", "Sqrt[2] a");
    check("Refine[Arg[-a], a > 0]", "Pi");
    check("Refine[Abs[x - 1], x > 2]", "-1 + x");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_refine_arg_complex);
    test_power_sqrt();
    test_log();
    test_trig();
    test_sign_abs_re();
    test_rounding_mod();
    test_element();
    test_predicates();
    test_options_scoping();
    test_errors_edges();
    test_stress_fixes();
    TEST(test_refine_re_of_real_terms);
    TEST(test_pattern_name_capture);
    test_directed_infinity();
    printf("All Refine tests passed.\n");
    return 0;
}
