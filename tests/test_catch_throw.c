/* Unit tests for Catch / Throw (Wolfram-Language exception control flow).
 *
 * Covers the full three-argument semantics: sentinel propagation through
 * ordinary function application, iteration/functional builtins, tag matching,
 * nested inner/outer catch, first-throw-wins, the f[value,tag] handler form,
 * and uncaught-throw handling. See src/eval.c (arg-loop short-circuit +
 * eval_report_uncaught_throw) and src/funcprog.c (builtin_catch/builtin_throw).
 *
 * NOTE: the uncaught-throw tests deliberately emit "Throw::nocatch" on stderr;
 * that is expected output, not a failure. */
#include "print.h"
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ----- Basic catch of a throw inside CompoundExpression ----- */
void test_basic() {
    assert_eval_eq("Catch[a; b; Throw[c]; d; e]", "c", 0);
    assert_eval_eq("Catch[42]", "42", 0);          /* no throw -> body value */
    assert_eval_eq("Catch[1 + 2]", "3", 0);
    assert_eval_eq("Catch[{a, Throw[b], c}]", "b", 0);
}

/* ----- Throw fired inside ordinary function application (Plus args) ----- */
void test_arbitrary_application() {
    assert_eval_eq("f[x_] := If[x > 10, Throw[overflow], x!]", "Null", 0);
    assert_eval_eq("Catch[f[2] + f[11]]", "overflow", 0);   /* 2nd arg throws */
    assert_eval_eq("Catch[f[2] + f[3]]", "8", 0);           /* no throw: 2!+3! */
    assert_eval_eq("Clear[f]", "Null", 0);
    /* Throw surfacing through a ReplaceAll-injected subexpression. */
    assert_eval_eq("Catch[a^2 + b^2 + c^2 /. b :> Throw[bbb]]", "bbb", 0);
}

/* ----- Exits from iteration / functional builtins ----- */
void test_iteration_exits() {
    assert_eval_eq("Catch[Do[If[i == 3, Throw[i]], {i, 10}]]", "3", 0);
    assert_eval_eq("Catch[Map[If[# > 2, Throw[#], #] &, {1, 2, 3, 4}]]", "3", 0);
    assert_eval_eq("Catch[Scan[If[# > 2, Throw[#]] &, {1, 2, 3, 4}]]", "3", 0);
    assert_eval_eq("Catch[Table[If[i == 4, Throw[i], i], {i, 10}]]", "4", 0);
    assert_eval_eq("Catch[Sum[If[i == 5, Throw[i], i], {i, 10}]]", "5", 0);
    assert_eval_eq("Catch[Nest[If[# > 8, Throw[#], # + 3] &, 0, 100]]", "9", 0);
    assert_eval_eq("Catch[Fold[If[#1 > 6, Throw[#1], #1 + #2] &, 0, {1,2,3,4,5}]]", "10", 0);
    assert_eval_eq("Catch[NestWhile[If[# > 5, Throw[#], # + 2] &, 0, True &]]", "6", 0);
    assert_eval_eq("Catch[FixedPoint[If[# > 5, Throw[#], # + 1] &, 0]]", "6", 0);
    assert_eval_eq("Catch[Which[False, a, True, Throw[w]]]", "w", 0);
    /* Anonymous-function map over a list (Throw need not be lexically inside). */
    assert_eval_eq("Catch[If[# < 0, Throw[#]] & /@ {1, 2, 0, -1, 5, 6}]", "-1", 0);
}

/* ----- Tag matching (form is often a pattern; tag re-evaluated) ----- */
void test_tags() {
    assert_eval_eq("Catch[Throw[a, u], u]", "a", 0);   /* tag matches form */
    assert_eval_eq("Catch[Throw[a, u], _]", "a", 0);   /* pattern form */
    /* Inner form does not match tag -> propagates to the outer Catch. */
    assert_eval_eq("Catch[Catch[Throw[v, outer], inner], outer]", "v", 0);
    /* Inner form matches -> outer never sees it. */
    assert_eval_eq("Catch[Catch[Throw[v, inner], inner], outer]", "v", 0);
    /* A tagless Throw is not caught by a form-Catch; the 1-arg outer gets it. */
    assert_eval_eq("Catch[Catch[Throw[v], t]]", "v", 0);
    /* Nested inner/outer catch by tag from the spec. */
    assert_eval_eq("g[x_] := If[x > 10, Throw[overflow], x!]", "Null", 0);
    assert_eval_eq("Catch[g[Catch[Throw[a, u], v]], u]", "a", 0);
    assert_eval_eq("Clear[g]", "Null", 0);
}

/* ----- First evaluated throw wins ----- */
void test_first_throw_wins() {
    assert_eval_eq("Catch[{Throw[a], Throw[b], Throw[c]}]", "a", 0);
    assert_eval_eq("Catch[Throw /@ {a, b, c}]", "a", 0);
}

/* ----- Module-local tag ----- */
void test_module_local_tag() {
    assert_eval_eq("Module[{u}, Catch[Throw[a, u], u]]", "a", 0);
}

/* ----- Three-argument handler form: Catch[expr, form, f] -> f[value, tag] ----- */
void test_handler_form() {
    assert_eval_eq("Catch[Throw[v, tg], tg, {#1, #2} &]", "{v, tg}", 0);
    assert_eval_eq("Catch[10, _, ff]", "10", 0);  /* no throw: handler unused */
    /* Full exception-handler dispatch from the spec. */
    assert_eval_eq("h[x_] := Which[x < 0, Throw[x, error[negative]], "
                   "x == 0, Throw[x, error[zero]], True, 1/Sqrt[x]]", "Null", 0);
    assert_eval_eq("hh[x_] := Catch[h[x], error[_], Function[{value, tag}, "
                   "tag /. {error[negative] :> Indeterminate, "
                   "error[zero] :> Infinity, _ :> Throw[value, tag]}]]", "Null", 0);
    assert_eval_eq("hh /@ {-1, 0, 1}", "{Indeterminate, Infinity, 1}", 0);
    assert_eval_eq("Clear[h]", "Null", 0);
    assert_eval_eq("Clear[hh]", "Null", 0);
}

/* ----- Uncaught throw handling (emits Throw::nocatch on stderr) ----- */
void test_uncaught() {
    assert_eval_eq("Throw[7]", "Hold[Throw[7]]", 0);
    assert_eval_eq("Throw[7, k, {#1, #2} &]", "{7, k}", 0);   /* 3-arg -> f[v,t] */
    assert_eval_eq("Catch[Throw[a, x], y]", "Hold[Throw[a, x]]", 0);  /* tag miss */
}

int main() {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_arbitrary_application);
    TEST(test_iteration_exits);
    TEST(test_tags);
    TEST(test_first_throw_wins);
    TEST(test_module_local_tag);
    TEST(test_handler_form);
    TEST(test_uncaught);

    printf("All catch/throw tests passed!\n");
    symtab_clear();
    return 0;
}
