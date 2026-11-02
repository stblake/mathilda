/* test_feyncalc_round12.c
 *
 * Kernel fixes found by bisecting FeynCalc's slowest test files (v0.6145):
 *   - SetOptions on a list of symbols, and with a list of rules
 *     (FeynCalc: SetOptions[{B1, B00, B11}, BReduce -> False] did nothing);
 *   - Together treats a user-head node with compound arguments, F[q + 1], as
 *     one opaque variable on the fast path (same answer, ~30x faster);
 *   - applying a pure function shares the unchanged parts of its body, so a
 *     large bound value keeps its fixed-point stamps (FeynCalc's Isolate).
 * Expected values computed by the reference kernel (InputForm strings).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

static void eq_if(const char* e, const char* expected) {
    char buf[2048], quoted[2048];
    snprintf(buf, sizeof buf, "ToString[%s, InputForm]", e);
    snprintf(quoted, sizeof quoted, "\"%s\"", expected);   /* a String prints quoted */
    assert_eval_eq(buf, quoted, 0);
}

static void test_setoptions_lists(void) {
    assert_eval_eq("Options[sog] = {opt -> True, o2 -> 1}; Options[soh] = {opt -> 3}; "
                   "sog[x_, OptionsPattern[]] := sored[x] /; OptionValue[opt]; Null", "Null", 0);
    eq_if("SetOptions[{sog, soh}, opt -> False]", "{{opt -> False, o2 -> 1}, {opt -> False}}");
    eq_if("{Options[sog], sog[2]}", "{{opt -> False, o2 -> 1}, sog[2]}");
    eq_if("SetOptions[sog, {opt -> 5, o2 -> 7}]", "{opt -> 5, o2 -> 7}");
    eq_if("SetOptions[sog, opt -> 6, o2 -> 8]", "{opt -> 6, o2 -> 8}");
}

static void test_together_compound_argument_atoms(void) {
    eq_if("Together[1/F[q + 1] + 1/F[q + 2]]", "(F[1 + q] + F[2 + q])/(F[1 + q]*F[2 + q])");
    eq_if("Together[F[q + 1]/(x - 1) - F[q + 1]/(x + 1)]", "(2*F[1 + q])/((-1 + x)*(1 + x))");
    eq_if("Together[F[1/a + 1/b]/(x^2 - 1) + 1/(x + 1)]",
          "(-1 + x + F[a^(-1) + b^(-1)])/((-1 + x)*(1 + x))");
    /* the cofactor sum NPointTo4Point's check builds: 0.28 s before */
    assert_eval_eq("tgG = Table[s[Min[i, j], Max[i, j]], {i, 4}, {j, 4}] + m IdentityMatrix[4]; "
                   "tgd = Expand[Det[tgG]]; "
                   "tgn = Sum[F[q + i] Expand[Det[Drop[tgG, {i}, {i}]] tgd], {i, 4}]; Null", "Null", 0);
    clock_t t0 = clock();
    assert_eval_eq("Expand[Together[tgn/tgd] - "
                   "Sum[F[q + i] Det[Drop[tgG, {i}, {i}]], {i, 4}]] === 0", "True", 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 0.2, "Together over F[q + i] atoms took %.2f s", secs);
}

static void test_pure_function_shares_body(void) {
    assert_eval_eq("Function[x, (x + #) &][a] === ((a + #) &)", "True", 0);
    eq_if("(#1 + #2 &)[a, b]", "a + b");
    eq_if("Function[{x, y}, f[x, g[2]]][1, 2]", "f[1, g[2]]");
    eq_if("(f[#, ##2] &)[1, 2, 3]", "f[1, 2, 3]");
    assert_eval_eq("pfex = Total[Table[Sqrt[TZ[i] + 2 ZZ[i]]/HH[ToString[i] <> ToString[i + 2]], "
                   "{i, 1, 10000}]]; pff[x_] := Select[Range[10], (# === x) &]; Null", "Null", 0);
    clock_t t0 = clock();
    assert_eval_eq("pff[pfex]", "{}", 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 0.15, "Select with a 10,000-term bound value took %.2f s", secs);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_setoptions_lists);
    TEST(test_together_compound_argument_atoms);
    TEST(test_pure_function_shares_body);
    printf("All FeynCalc round-12 kernel tests passed.\n");
    return 0;
}
