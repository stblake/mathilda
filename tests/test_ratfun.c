/*
 * test_ratfun.c -- the factored rational-function engine (src/poly/ratfun.c):
 * Together, Cancel and Apart over Q, Q(params) and GF(p).
 *
 * The forms are those of a user-level call, which renders the factored form
 * (an internal call from another builtin gets the expanded one; see
 * test_internal_form).  Known differences kept on purpose: Apart's Modulus
 * coefficients are residues (6/(3 + x), not 1/(6 (3 + x))); a Gaussian
 * content is rational, not a Gaussian-integer gcd; Apart also takes Extension.
 *
 * Three kinds of check:
 *   - exact output forms for the cases whose form is settled;
 *   - regressions found while the engine was written (a hang, an incomplete
 *     reduction), each with the input that exposed it;
 *   - a seeded stress: random rational functions whose Together / Cancel /
 *     Apart must equal the input at random rational points and come back in
 *     lowest terms.
 */
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    expr_free(e);
    expr_free(r);
    return s;
}

static int g_failures = 0;

/* Reports every mismatch (not only the first) and fails the run at the end. */
static void assert_form(const char* input, const char* expected) {
    char* s = eval_str(input);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        free(s);
        g_failures++;
        return;
    }
    printf("PASS: %s -> %s\n", input, expected);
    free(s);
}

static void assert_true(const char* input) { assert_form(input, "True"); }

void test_together_forms(void) {
    assert_form("Together[1/x + 1/(x + 1)]", "(1 + 2 x)/(x (1 + x))");
    assert_form("Together[x/2 + y/3]", "1/6 (3 x + 2 y)");
    assert_form("Together[1/(2 x - 2) + 1/(3 x + 3)]", "(1 + 5 x)/(6 (-1 + x) (1 + x))");
    assert_form("Together[a/(x y) + b/(x z)]", "(b y + a z)/(x y z)");
    assert_form("Together[1/(1 + E^x) + Log[x]/(E^(2 x) - 1)]", "(-1 + E^x + Log[x])/((-1 + E^x) (1 + E^x))");
    /* the generator is oriented positively: 1/E^x is a denominator */
    assert_form("Together[1 + E^(-x)]", "E^(-x) (1 + E^x)");
}

void test_cancel_forms(void) {
    assert_form("Cancel[(x^2 - 1)/(x - 1)]", "1 + x");
    assert_form("Cancel[(x^4 - 1)/(x^6 - 1)]", "(1 + x^2)/(1 + x^2 + x^4)");
    assert_form("Cancel[(x^2 - 1)/(x^3 - 1) + (x^2 - 4)/(x - 2)]", "2 + x + (1 + x)/(1 + x + x^2)");
    /* a term with nothing to cancel is left as it is (an E^(...) exponent
     * written 1/(1 - x^2) must not become -1/(-1 + x^2)) */
    assert_form("Cancel[-1 + E^x + 1/x + 1/(1 - x^2)]", "-1 + E^x + 1/x + 1/(1 - x^2)");
    /* a nested sum of fractions is left to the classical path (which
     * combines it; the factored engine would leave 1/(1 + 1/x) alone) */
    assert_form("Cancel[1/(1 + 1/x)]", "x/(1 + x)");
}

void test_regressions(void) {
    /* the monomial split of (3x^3 - 2x^4) left a second (2x - 3) unmerged and
     * the single-term pass looped on the equal pair */
    assert_true("TimeConstrained[Together[(x + 4)/((3 x^3 - 2 x^4)^2 (3 - 2 x))], 10] =!= $Aborted");
    /* a monomial factor x^3 kept whole was treated as coprime to x */
    assert_form("Together[(6 x^4 + 2 x^3)/((5 x - x^3 + 5 x^2)^2 (3 x + 1)^2)]",
                "(2 x)/((1 + 3 x) (-5 - 5 x + x^2)^2)");
}

/* The v0.6074-75 fast paths: the inner sums of a continued fraction skip the
 * squarefree split, linear denominators skip the coprimality and cancellation
 * gcds, and the factored answer is built canonical without an evaluation (an
 * answer equal to its own re-evaluation, and to the evaluated build under
 * MATHILDA_RATFUN_DIRECT_CHECK).  The depth-3 denominator's form is not pinned:
 * a reference splits it against the deeper 1 + x as (1 + x)(3 + 2 x + x^2). */
void test_fast_paths(void) {
    assert_true("t = Together[Nest[x + 1/(1 + #) &, x, 3]]; "
                "Denominator[t] === 3 + 5 x + 3 x^2 + x^3 || "
                "Denominator[t] === (1 + x) (3 + 2 x + x^2)");
    assert_true("t = Together[Nest[x + 1/(1 + #) &, x, 24]]; "
                "(t - Nest[x + 1/(1 + #) &, x, 24] /. x -> 7/3) === 0 && "
                "Exponent[Denominator[t], x] == 24");
    assert_form("Together[Sum[1/(x - i), {i, 1, 4}] - 1/(x - 4)]",
                "(11 - 12 x + 3 x^2)/((-3 + x) (-2 + x) (-1 + x))");
    assert_form("Together[x/(x - 1) - 1/(x - 1) + 1/(x + 2)]", "(3 + x)/(2 + x)");
    assert_true("t = Together[Sum[(i^2 + 1)/(x - i), {i, 1, 9}]]; t === ToExpression[ToString[t, InputForm]]");
    assert_true("t = Together[a/(x + y) + b/(x - y) + 3/(2 x + 1)]; t === ToExpression[ToString[t, InputForm]]");
}

void test_apart_forms(void) {
    assert_form("Apart[1/((x + a) (x + b)), x]", "1/((-a + b) (a + x)) + 1/((a - b) (b + x))");
    assert_form("Apart[x^3/(x - a), x]", "a^2 + a x + x^2 + a^3/(-a + x)");
    assert_form("Apart[x^2/(a (x - 1)), x]", "1/a + 1/(a (-1 + x)) + x/a");
    assert_form("Apart[(x + a)/((x - 1) (x - 2) a), x]", "(2 + a)/(a (-2 + x)) + (-1 - a)/(a (-1 + x))");
    assert_form("Apart[(x^2 + 1)/((x - 1)^3 (x + 2)^2), x]",
                "2/(9 (-1 + x)^3) + 2/(27 (-1 + x)^2) + 1/(27 (-1 + x)) - 5/(27 (2 + x)^2) - 1/(27 (2 + x))");
    assert_form("Apart[1/(x^2 - 1)]", "1/(2 (-1 + x)) - 1/(2 (1 + x))");
    /* five parameters: the classical elimination took 2 s, this is native */
    assert_true("Together[Apart[1/((x + a) (x + b) (x + c) (x + d) (x + f)), x] - "
                "1/((x + a) (x + b) (x + c) (x + d) (x + f))] === 0");
}

void test_modulus(void) {
    assert_form("Together[1/(x - 1) + 1/(x + 1), Modulus -> 3]", "(2 x)/((1 + x) (2 + x))");
    assert_form("Cancel[(x^2 - 1)/(x - 1), Modulus -> 2]", "1 + x");
    assert_form("Cancel[(x^7 - x)/(x^2 - x), Modulus -> 7]", "1 + x + x^2 + x^3 + x^4 + x^5");
    assert_form("Apart[1/(x^3 - 1), x, Modulus -> 7]", "6/(3 + x) + 3/(5 + x) + 5/(6 + x)");
    assert_form("Together[1/(x - 1) + 1/(x + 1), Modulus -> 0]", "(2 x)/((-1 + x) (1 + x))");
    /* a composite modulus: Together::modp, the call stays unevaluated */
    assert_form("Quiet[Together[1/(x - 1) + 1/(x + 1), Modulus -> 4]]",
                "Together[1/(-1 + x) + 1/(1 + x), Modulus -> 4]");
    assert_form("Check[Together[1/x, Modulus -> 4], failed]", "failed");
}

/* The same call from inside another builtin keeps the expanded form its
 * internal callers read: Simplify's Together, Sum's, Apart's classical path. */
void test_internal_form(void) {
    assert_form("Together[1/(x - 1) + 1/(x + 1)]", "(2 x)/((-1 + x) (1 + x))");
    assert_form("Map[Together, {1/(x - 1) + 1/(x + 1)}]", "{(2 x)/((-1 + x) (1 + x))}");
    assert_form("Module[{t = Together[1/(x - 1) + 1/(x + 1)]}, t]", "(2 x)/((-1 + x) (1 + x))");
    assert_form("Hold[Together[1/(x - 1) + 1/(x + 1)]] // ReleaseHold", "(2 x)/((-1 + x) (1 + x))");
    assert_form("LerchPhi[2, -2, a]", "-6 + 4 a - a^2");
}

void test_options(void) {
    assert_form("Options[Together]", "{Extension -> None, Modulus -> 0, Trig -> False}");
    assert_form("Options[Apart]", "{Extension -> None, Modulus -> 0, Trig -> False}");
    /* the default settings are accepted (and dropped) */
    assert_form("Together[1/x + 1/(x + 1), Trig -> False]", "(1 + 2 x)/(x (1 + x))");
    assert_form("Apart[1/(x^2 - 1), x, Trig -> False, Modulus -> 0]", "1/(2 (-1 + x)) - 1/(2 (1 + x))");
    assert_form("Check[Cancel[x/(x - Sqrt[2]), Modulus -> 3, Extension -> Sqrt[2]], failed]", "failed");
}

void test_gaussian(void) {
    assert_form("Together[1/(x - I) + 1/(x + I)]", "(2 x)/((-I + x) (I + x))");
    assert_form("Cancel[(x^2 + 1)/(x^2 - 2 I x - 1)]", "(I + x)/(-I + x)");
    /* integral Gaussian coefficients, not a monic denominator */
    assert_form("Together[(-1 - I)/(x + 3/5 + I/5)]", "(-5 - 5*I)/(3 + I + 5 x)");
    assert_form("Apart[1/((x - I) (x + 2 I)), x]", "(-1/3*I)/(-I + x) + (1/3*I)/(2*I + x)");
    /* Gaussian input: the denominator is factored over Q(i) */
    assert_form("Apart[I/(x^2 + 1), x]", "1/(2 (-I + x)) - 1/(2 (I + x))");
    /* a degree-62 gcd over Q(Sqrt[2]): the modular field gcd, 4 s -> ms */
    assert_true("TimeConstrained[Cancel[Expand[(x^2 - 2) (x^60 + 3 x + 1)]/"
                "Expand[(x - Sqrt[2]) (x^60 - 5 x^7 + 2)], Extension -> Automatic], 10] === "
                "Cancel[(Sqrt[2] + x) (x^60 + 3 x + 1)/(x^60 - 5 x^7 + 2)]");
}

/* Random rational functions (sums of fractions, products, nested sums) in x
 * and up to two parameters: each answer equals the input at a random rational
 * point, and Together's comes back in lowest terms. */
void test_stress(void) {
    assert_true(
        "SeedRandom[20261028];"
        "rnd[vs_] := Total[Table[RandomInteger[{-4, 4}] Times @@ (vs^RandomInteger[{0, 2}, Length[vs]]), {4}]];"
        "pts = {x -> 37/101, a -> 29/83, b -> -53/97};"
        "ok = True;"
        "Do["
        "  vs = Take[{x, a, b}, RandomInteger[{1, 3}]];"
        "  f = Total[Table[rnd[vs]/(rnd[vs] (rnd[vs])^RandomInteger[{1, 2}]), {RandomInteger[{1, 4}]}]]"
        "      + 1/(rnd[vs] + 1/rnd[vs]);"
        "  If[Quiet[Check[f /. pts, $Failed]] === $Failed || !FreeQ[f /. pts, ComplexInfinity | Indeterminate], Continue[]];"
        "  t = Together[f]; c = Cancel[f];"
        "  p = If[Length[vs] > 2, f, Apart[f, x]];"   /* over Q(a, b) these get large */
        "  ok = ok && (t - f /. pts) === 0 && (c - f /. pts) === 0 && (p - f /. pts) === 0 &&"
        "       FreeQ[PolynomialGCD[Numerator[t], Denominator[t]], x | a | b],"
        "  {200}];"
        "ok");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running ratfun tests...\n");
    TEST(test_together_forms);
    TEST(test_cancel_forms);
    TEST(test_regressions);
    TEST(test_fast_paths);
    TEST(test_apart_forms);
    TEST(test_modulus);
    TEST(test_options);
    TEST(test_internal_form);
    TEST(test_gaussian);
    TEST(test_stress);
    if (g_failures) {
        printf("%d ratfun test(s) FAILED\n", g_failures);
        return 1;
    }
    printf("All ratfun tests passed!\n");
    return 0;
}
