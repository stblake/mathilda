/* test_experimental_diractr.c
 *
 * Experimental`FeynCalc`DiracTr (src/experimental/feyncalc/diractr.c), the
 * native Dirac trace.  Its contract is exact: the value is SameQ -- term
 * order included -- to that of the recursive reference definition
 * (Woxi issue #2142), which every test here evaluates alongside it:
 *
 *     diracTr[] := 4;
 *     diracTr[a_, b_] := 4 g[a, b];
 *     diracTr[a_, rest___] := Module[{list = {rest}, n},
 *       Sum[(-1)^(m + 1) * Expand[g[a, list[[m]]] * diracTr @@ Delete[list, m]],
 *         {m, 1, Length[list]}]];
 *
 * Covered: the small cases written out; distinct arguments in and out of
 * canonical order, atoms and compound arguments; repeated arguments (powers,
 * merged and cancelling terms) including a seeded random sweep; the issue's
 * n = 14 counts; that the fixed-point stamp hides nothing (re-evaluation,
 * later definitions of g); the fallback when g, an argument or the
 * arithmetic carries rules; the size guard and its message.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

#define DT "Experimental`FeynCalc`DiracTr"

/* Evaluate for effect. */
static void run(const char* src) {
    Expr* p = parse_expression(src);
    if (!p) { fprintf(stderr, "FAIL: does not parse: %s\n", src); abort(); }
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

static void is(const char* lhs, const char* rhs) {
    char buf[8192];
    snprintf(buf, sizeof buf, "(%s) === (%s)", lhs, rhs);
    assert_eval_eq(buf, "True", 0);
}

/* DiracTr[args] === diracTr[args] (the reference recursion) */
static void same_as_ref(const char* args) {
    char l[4096], r[4096];
    snprintf(l, sizeof l, DT "[%s]", args);
    snprintf(r, sizeof r, "diracTr[%s]", args);
    is(l, r);
}

static void define_reference(void) {
    run("diracTr[] := 4; diracTr[a_, b_] := 4 g[a, b]; "
        "diracTr[a_, rest___] := Module[{list = {rest}, n}, "
        "Sum[(-1)^(m + 1) * Expand[g[a, list[[m]]] * diracTr @@ Delete[list, m]], "
        "{m, 1, Length[list]}]]");
}

static void test_small_cases(void) {
    is(DT "[]", "4");
    is(DT "[a]", "0");
    is(DT "[a, b, c]", "0");
    is(DT "[a, b, c, d, e]", "0");
    is(DT "[a, b]", "4 g[a, b]");
    is(DT "[b, a]", "4 g[b, a]");
    is(DT "[a, b, c, d]", "4 g[a, d] g[b, c] - 4 g[a, c] g[b, d] + 4 g[a, b] g[c, d]");
    assert_eval_eq("InputForm[" DT "[a, b, c, d]]",
                   "4*g[a, d]*g[b, c] - 4*g[a, c]*g[b, d] + 4*g[a, b]*g[c, d]", 0);
    is("FullForm[" DT "[a, b]]", "FullForm[Times[4, g[a, b]]]");
    for (int k = 0; k <= 5; k++) {
        char buf[64];
        snprintf(buf, sizeof buf, "Array[x, %d]", k);
        char l[128], r[128];
        snprintf(l, sizeof l, DT " @@ %s", buf);
        snprintf(r, sizeof r, "diracTr @@ %s", buf);
        is(l, r);
    }
}

static void test_distinct_against_reference(void) {
    for (int k = 2; k <= 12; k += 2) {
        char l[128], r[128];
        snprintf(l, sizeof l, DT " @@ Array[x, %d]", k);
        snprintf(r, sizeof r, "diracTr @@ Array[x, %d]", k);
        is(l, r);
    }
    /* out of canonical order: per-term factor sort and radix order matter */
    same_as_ref("z, b, a, y");
    same_as_ref("d, c, b, a, f, e");
    same_as_ref("x[3], x[1], x[4], x[2], x[6], x[5]");
    same_as_ref("x[8], x[1], x[7], x[2], x[6], x[3], x[5], x[4]");
    is(DT " @@ Reverse[Array[x, 10]]", "diracTr @@ Reverse[Array[x, 10]]");
    /* atoms and compound arguments of every kind g can hold */
    same_as_ref("1, p + q, f[x], a");
    same_as_ref("mu, 2, nu, 1");
    same_as_ref("\"s\", k, 3/2, p^2, Sin[t], {1, 2}");
    same_as_ref("p1 + p2, mu, p1 - p2, nu, -k, rho");
    same_as_ref("x[1, 2], x[1], x[2, 1], x[0]");
}

static void test_repeated_against_reference(void) {
    same_as_ref("a, a");
    same_as_ref("a, a, b, b");
    same_as_ref("a, b, a, b");
    same_as_ref("a, b, b, a");
    same_as_ref("a, a, a, a");
    same_as_ref("mu, nu, mu, nu, rho, rho");
    same_as_ref("mu, nu, rho, mu, nu, rho");
    same_as_ref("a, a, a, a, a, a");
    same_as_ref("a, b, a, b, a, b, a, b");
    same_as_ref("p, mu, p, mu, q, nu, q, nu");
    is(DT "[a, a, b, b]", "4 g[a, a] g[b, b]");
    is(DT "[a, a, a, a]", "4 g[a, a]^2");
    is(DT "[a, b, a, b]", "4 g[a, b]^2 + 4 g[a, b] g[b, a] - 4 g[a, a] g[b, b]");
    /* seeded random multisets, n <= 10, over a small alphabet */
    run("SeedRandom[2142]; dtCases = Flatten[Table[RandomChoice[{a, b, c, d}, k], "
        "{k, {2, 4, 6, 8, 10}}, {12}], 1];");
    assert_eval_eq("And @@ ((" DT " @@ # === diracTr @@ #) & /@ dtCases)", "True", 0);
    run("SeedRandom[7]; dtCases = Table[RandomChoice[{x[1], x[2], mu, 1}, 8], {20}];");
    assert_eval_eq("And @@ ((" DT " @@ # === diracTr @@ #) & /@ dtCases)", "True", 0);
}

static void test_issue_counts(void) {
    /* Woxi #2142: n = 14 has 13!! = 135135 terms and 91 distinct g */
    run("dt14 = " DT " @@ Array[x, 14];");
    assert_eval_eq("Length[dt14]", "135135", 0);
    assert_eval_eq("Length[Union[Cases[dt14, _g, Infinity]]]", "91", 0);
    assert_eval_eq("Count[dt14, Times[4, __]]", "67568", 0);
    assert_eval_eq("Count[dt14, Times[-4, __]]", "67567", 0);
    /* the signs sum to 1 (the Pfaffian of the all-ones antisymmetric matrix) */
    assert_eval_eq("dt14 /. g[_, _] -> 1", "4", 0);
    assert_eval_eq("Length[" DT " @@ Array[x, 10]]", "945", 0);
    run("Clear[dt14]");
}

static void test_canonical_and_stamp(void) {
    /* the stamped value is what re-canonicalising it gives */
    run("dt8 = " DT " @@ Reverse[Array[x, 8]];");
    is("dt8", "Plus @@ (List @@ dt8)");
    is("dt8", "Expand[dt8]");
    is("dt8", "Sort[dt8]");
    run("dtr = " DT "[a, b, a, b, c, c];");
    is("dtr", "Plus @@ (List @@ dtr)");
    is("dtr", "Expand[dtr]");
    /* a definition made afterwards reaches the stored value */
    run("dt4 = " DT "[a, b, c, d];");
    run("g[a, b] = 0;");
    is("dt4", "4 g[a, d] g[b, c] - 4 g[a, c] g[b, d]");
    run("Clear[g]");
    is("dt4", "4 g[a, d] g[b, c] - 4 g[a, c] g[b, d] + 4 g[a, b] g[c, d]");
    run("Clear[dt8, dtr, dt4]");
}

static void test_fallback_when_g_has_rules(void) {
    run("g[i_, i_] := dim");
    same_as_ref("mu, mu");
    same_as_ref("mu, mu, nu, nu");
    same_as_ref("mu, nu, mu, nu");
    same_as_ref("a, b, c, d, e, f");
    run("Clear[g]");

    run("SetAttributes[g, Orderless]");
    same_as_ref("a, b, a, b");
    same_as_ref("d, c, b, a");
    is(DT "[b, a]", "4 g[a, b]");
    run("ClearAttributes[g, Orderless]");

    run("g[u_, v_] := u + v");   /* g evaluates to a sum: Expand must apply */
    same_as_ref("a, b, c, d");
    same_as_ref("x[1], x[2], x[3], x[4], x[5], x[6]");
    run("Clear[g]");

    run("x /: g[x[1], y_] := y");  /* an UpValue on an argument's head */
    same_as_ref("x[1], x[2], x[3], x[4]");
    run("Clear[x]");

    /* everything plain again: the fast path agrees once more */
    same_as_ref("x[1], x[2], x[3], x[4], x[5], x[6]");
}

static void test_guard_and_attributes(void) {
    /* 25!! > 2^32 terms: message, unevaluated */
    assert_eval_eq("Check[" DT " @@ Range[26], caught]", "caught", 0);
    assert_eval_eq("Quiet[Head[" DT " @@ Range[26]]]", DT, 0);
    /* odd n is 0 at any size, with no message */
    assert_eval_eq("Check[" DT " @@ Range[101], caught]", "0", 0);
    assert_eval_eq("Attributes[" DT "]", "{Protected}", 0);
    { SymbolDef* d = symtab_get_def(DT); ASSERT(d && d->docstring && strlen(d->docstring) > 200); }
}

/* eval.c: x = rhs skips its O(tree) "does rhs mention x?" walk when rhs is
 * still a reusable fixed point after the binding (r = <trace> walked all
 * 135135 terms).  The cases where the walk mattered must be unchanged. */
static void test_set_of_a_fixed_point(void) {
    assert_eval_eq("sn1 = 1; sn1 = sn1 + 1; sn1", "2", 0);
    assert_eval_eq("sz = sg2[sw]; sw = 5; sz", "sg2[5]", 0);
    assert_eval_eq("sp = sf[sq]; sq = 3; sp = sp; sp", "sf[3]", 0);
    assert_eval_eq("sl = {sa, sb}; sl = Append[sl, sl]; sl", "{sa, sb, {sa, sb}}", 0);
    assert_eval_eq("Quiet[sx = sx]", "sx", 0);
    assert_eval_eq("sk = 0; sv = Table[sk++, {3}]; {sv, sk}", "{{0, 1, 2}, 3}", 0);
    /* the value assigned is the trace itself, stamps intact */
    assert_eval_eq("st = " DT " @@ Array[x, 10]; {Length[st], st === (" DT " @@ Array[x, 10])}",
                   "{945, True}", 0);
    run("Clear[sn1, sz, sw, sp, sq, sl, sx, sk, sv, st]");
}

int main(void) {
    symtab_init();
    core_init();
    define_reference();
    TEST(test_small_cases);
    TEST(test_distinct_against_reference);
    TEST(test_repeated_against_reference);
    TEST(test_issue_counts);
    TEST(test_canonical_and_stamp);
    TEST(test_fallback_when_g_has_rules);
    TEST(test_guard_and_attributes);
    TEST(test_set_of_a_fixed_point);
    printf("All Experimental`FeynCalc`DiracTr tests passed.\n");
    return 0;
}
