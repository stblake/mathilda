/* Builtin-audit regressions: solutions of equations (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Root#1: Root[f, k] numbers the deg f roots of f counted with multiplicity. */
static void test_root_multiplicity(void) {
    audit_chk("Root#1", "N[Root[#^2 (# - 1) &, 2]]", "0.0");
    audit_chk("Root#1", "N[Root[#^2 (# - 1) &, 3]]", "1.0");
    audit_true("Root#1", "Abs[N[Root[(#^2 - 2)^2 &, 3]] - Sqrt[2]] < 10^-12");
    audit_true("Root#1",
               "Max[Abs[Table[N[Root[(#^2 - 2)^2 &, k]], {k, 4}] - "
               "{-Sqrt[2], -Sqrt[2], Sqrt[2], Sqrt[2]}]] < 10^-12");
    audit_true("Root#1", "Abs[N[Root[(# - 1)^3 (# + 2) &, 4], 30] - 1] < 10^-25");
    audit_chk("Root#1", "Check[N[Root[#^2 (# - 1) &, 4]], fired, Root::indx]", "fired");
    /* squarefree polynomials are unchanged */
    audit_true("Root#1", "Abs[N[Root[#^3 - 2 &, 1]] - 2^(1/3)] < 10^-12");
}

/* Solve#1: fusing Log[x] + Log[x - 1] into Log[x (x - 1)] admitted x -> -1,
 * where the left side is Log[2] + 2 Pi I. */
static void test_solve_log_fusion(void) {
    audit_chk("Solve#1", "Solve[Log[x] + Log[x - 1] == Log[2], x]", "{{x -> 2}}");
    audit_true("Solve#1",
               "PossibleZeroQ[(Log[x] + Log[x - 1] - Log[2]) /. "
               "First[Solve[Log[x] + Log[x - 1] == Log[2], x]]]");
    audit_chk("Solve#1", "Solve[Log[x] == 2, x]", "{{x -> E^2}}");
}

/* FindInstance#1 / #2: under Modulus -> p an equation holds in Z/pZ. */
static void test_findinstance_modulus(void) {
    audit_chk("FindInstance#1", "FindInstance[x^2 == 2, x, Modulus -> 7]", "{{x -> 3}}");
    audit_chk("FindInstance#1", "FindInstance[2 x == 1, x, Modulus -> 7]", "{{x -> 4}}");
    audit_chk("FindInstance#1", "FindInstance[x^2 == 3, x, Modulus -> 7]", "{}");   /* non-residue */
    audit_true("FindInstance#1",
               "Mod[(x^2 - 2) /. First[FindInstance[x^2 == 2, x, Modulus -> 7]], 7] == 0");
    audit_chk("FindInstance#2", "FindInstance[x^2 == 2, x, Integers, Modulus -> 7]", "{{x -> 3}}");
    audit_chk("FindInstance#2", "FindInstance[x^2 == 4, x]", "{{x -> -2}}");
}

/* Exists#1 / Resolve#1: an inexact polynomial at its own machine root is an
 * approximate zero, not a sign. */
static void test_reduce_inexact_roots(void) {
    audit_chk("Exists#1", "Resolve[Exists[x, x^2 == 2.0], Reals]", "True");
    audit_chk("Resolve#1", "Resolve[Exists[x, x^2 == 2.0], Reals]", "True");
    audit_chk("Resolve#1", "Resolve[Exists[x, x^2 == 4.0], Reals]", "True");
    audit_chk("Resolve#1", "Resolve[Exists[x, x^2 == -2.0], Reals]", "False");
    audit_true("Exists#1",
               "Module[{r = Reduce[x^2 == 2.0, x, Reals]}, Head[r] === Or && "
               "Abs[r[[2, 2]] - Sqrt[2]] < 10^-12 && Abs[r[[1, 2]] + Sqrt[2]] < 10^-12]");
    audit_true("Exists#1",
               "Module[{r = Reduce[x^3 == 2.5, x, Reals]}, Head[r] === Equal && "
               "Abs[r[[2]] - 2.5^(1/3)] < 10^-12]");
    audit_chk("Exists#1", "Reduce[x^2 + 1.0 == 0, x, Reals]", "False");
    audit_true("Exists#1",   /* all four real roots of an inexact quartic */
               "Length[Reduce[x^4 - 5.2 x^2 + 4.1 == 0, x, Reals]] == 4");
    audit_chk("Exists#1", "Reduce[x^2 == 4, x, Reals]", "x == -2 || x == 2");
}

/* SolveAlways#1: Sin[x], Log[x], ... are functions of x, not parameters. */
static void test_solvealways_transcendental(void) {
    audit_chk("SolveAlways#1", "SolveAlways[a Sin[x] + b == 0, x]", "{{a -> 0, b -> 0}}");
    audit_chk("SolveAlways#1", "SolveAlways[a Log[x] + b == 0, x]", "{{a -> 0, b -> 0}}");
    audit_chk("SolveAlways#1", "SolveAlways[a Sin[x] + b Cos[x] == 0, x]", "{{a -> 0, b -> 0}}");
    audit_chk("SolveAlways#1", "SolveAlways[a Cos[x]^2 + b Sin[x]^2 == 1, x]", "{{a -> 1, b -> 1}}");
    audit_chk("SolveAlways#1", "SolveAlways[a Exp[x] + b == 2 Exp[x] + 3, x]", "{{a -> 2, b -> 3}}");
    audit_chk("SolveAlways#1", "SolveAlways[a Sin[x] + b == Sin[x] + 2, x]", "{{a -> 1, b -> 2}}");
    /* no parameters: a true identity is {{}}, a false one {} */
    audit_chk("SolveAlways#1", "SolveAlways[Sin[x]^2 + Cos[x]^2 == 1, x]", "{{}}");
    audit_chk("SolveAlways#1", "SolveAlways[Sin[2 x] == 2 Sin[x] Cos[x], x]", "{{}}");
    audit_chk("SolveAlways#1", "SolveAlways[Sin[x] == Sin[x] + 1, x]", "{}");
    /* outside the exp-log family: stays unevaluated rather than solving for Sqrt[x] */
    audit_chk("SolveAlways#1", "SolveAlways[a Sqrt[x] + b == 0, x]", "SolveAlways[b + a Sqrt[x] == 0, x]");
    /* polynomial identities unchanged */
    audit_chk("SolveAlways#1", "SolveAlways[a x^2 + b x + c == 0, x]", "{{a -> 0, b -> 0, c -> 0}}");
    audit_chk("SolveAlways#1",
              "SolveAlways[a/(x - 1) + b/(x - 2) == (3 x - 5)/((x - 1) (x - 2)), x]",
              "{{a -> 2, b -> 1}}");
}

/* ToRadicals#1: the inexact Cardano branch cancelled to 0. and every root was 0. */
static void test_toradicals_inexact_cubic(void) {
    audit_true("ToRadicals#1", "Abs[ToRadicals[Root[#^3 - 2.0 &, 1]] - 2^(1/3)] < 10^-14");
    audit_chk("ToRadicals#1", "Head[ToRadicals[Root[#^3 - 2.0 &, 1]]]", "Real");
    audit_true("ToRadicals#1",
               "Max[Abs[Table[ToRadicals[Root[#^3 - 2.0 &, k]] - N[Root[#^3 - 2 &, k]], {k, 3}]]] < 10^-12");
    audit_true("ToRadicals#1",
               "Max[Abs[Table[ToRadicals[Root[#^3 + 2.0 &, k]] - N[Root[#^3 + 2 &, k]], {k, 3}]]] < 10^-12");
    audit_true("ToRadicals#1",
               "Max[Abs[Table[(ToRadicals[Root[#^3 - 3. # + 1 &, k]])^3 - "
               "3 ToRadicals[Root[#^3 - 3. # + 1 &, k]] + 1, {k, 3}]]] < 10^-12");
    /* exact radicals unchanged */
    audit_true("ToRadicals#1", "PossibleZeroQ[ToRadicals[Root[#^3 - 2 &, 1]] - 2^(1/3)]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* CountRoots#1: the zero polynomial vanishes everywhere */
    audit_true("CountRoots#1", "CountRoots[0, x] === Infinity");
    audit_true("CountRoots#1", "CountRoots[0, {x, 0, 1}] === Infinity");
    audit_true("CountRoots#1", "CountRoots[x^2 - 1, x] === 2");
    audit_true("CountRoots#1", "CountRoots[5, x] === 0");
    /* NotElement#3: NotElement of an infinity */
    audit_true("NotElement#3", "NotElement[Infinity, Reals] === True");
    audit_true("NotElement#3", "NotElement[3, Reals] === False");
}

int main(void) {
    audit_init();
    TEST(test_root_multiplicity);
    TEST(test_solve_log_fusion);
    TEST(test_findinstance_modulus);
    TEST(test_reduce_inexact_roots);
    TEST(test_solvealways_transcendental);
    TEST(test_toradicals_inexact_cubic);
    TEST(test_s11_special_values);
    return audit_finish("audit_solutions_of_equations");
}
