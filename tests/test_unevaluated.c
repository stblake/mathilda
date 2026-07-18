#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "attr.h"

/* Unevaluated has attributes {HoldAllComplete, Protected} */
void test_unevaluated_attributes() {
    assert_eval_eq("Attributes[Unevaluated]", "{HoldAllComplete, Protected}", 0);
}

/* Unevaluated evaluates to itself (HoldAllComplete prevents inner evaluation) */
void test_unevaluated_self() {
    assert_eval_eq("Unevaluated[1+2]", "Unevaluated[1 + 2]", 0);
    assert_eval_eq("Unevaluated[Sequence[a, b]]", "Unevaluated[Sequence[a, b]]", 0);
    assert_eval_eq("Unevaluated[Evaluate[1+2]]", "Unevaluated[Evaluate[1 + 2]]", 0);
}

/* f[Unevaluated[expr]] strips the wrapper and passes expr unevaluated */
void test_unevaluated_length_plus() {
    /* Explicit Plus[...] form: arg_count is the number of summands. */
    assert_eval_eq("Length[Unevaluated[Plus[5, 6, 7, 8]]]", "4", 0);
    assert_eval_eq("Length[Unevaluated[Plus[a, b, c]]]", "3", 0);
    /* Infix parses n-ary, so Length matches the number of summands/factors.
     * This also regresses the parser: `a+b+c` must NOT become Plus[Plus[a,b],c]. */
    assert_eval_eq("Length[Unevaluated[5 + 6 + 7 + 8]]", "4", 0);
    assert_eval_eq("Length[Unevaluated[1 + 1/0 + 0/0]]", "3", 0);
    assert_eval_eq("Length[Unevaluated[1 + 1/0 + 0/0 + x]]", "4", 0);
    assert_eval_eq("Length[Unevaluated[1 + 1/0 + 0/0 + x + 1]]", "5", 0);
    assert_eval_eq("Length[Unevaluated[a * b * c * d]]", "4", 0);
    /* Mix of + and -: each term is a top-level summand. */
    assert_eval_eq("Length[Unevaluated[a + b - c + d]]", "4", 0);
}

/* Sequences directly inside Unevaluated are NOT flattened */
void test_unevaluated_sequence_not_flattened() {
    assert_eval_eq("Length[Unevaluated[Sequence[a, b]]]", "2", 0);
    assert_eval_eq("Length[Unevaluated[Sequence[a, b, c, d]]]", "4", 0);
}

/* Unevaluated stops Evaluate inside a held function: Evaluate forces the arg,
 * Unevaluated keeps 1+2 from evaluating (result is 1+2, not 3), and the wrapper
 * is stripped per WMA semantics. */
void test_unevaluated_stops_evaluate() {
    assert_eval_eq("Hold[Evaluate[Unevaluated[1+2]]]", "Hold[1 + 2]", 0);
}

/* Unevaluated inside a HoldAll function is stripped even in held positions
 * (WMA semantics); the exposed content stays unevaluated because HoldAll holds it. */
void test_unevaluated_in_holdall() {
    assert_eval_eq("ClearAll[f]; SetAttributes[f, HoldAll]; f[Unevaluated[1+2]]",
                   "f[1 + 2]", 0);
}

/* Unevaluated is stripped in held AND unheld positions (only HoldAllComplete
 * preserves it). The exposed content stays unevaluated in held positions. */
void test_unevaluated_mixed_hold() {
    assert_eval_eq("Hold[Unevaluated[1+2]]", "Hold[1 + 2]", 0);
    /* HoldFirst on mySym: arg 1 is held (1+2 stays unevaluated) but the
     * Unevaluated wrapper is still stripped; arg 2 evaluates normally. */
    assert_eval_eq("ClearAll[mySym]; Attributes[mySym]; SetAttributes[mySym, HoldFirst]; mySym[Unevaluated[1+2], 3+4]",
                   "mySym[1 + 2, 7]", 0);
}

/* Nested Unevaluated: each call strips one layer */
void test_unevaluated_nested() {
    assert_eval_eq("Length[Unevaluated[Unevaluated[Plus[a, b, c]]]]", "1", 0);
}

/* Head itself does not get stripped when Unevaluated is an argument to a
 * symbolic (non-builtin) function with no hold attributes */
void test_unevaluated_symbolic_head() {
    /* After stripping, h receives Plus[a, b] which evaluates normally (printed as "a + b"). */
    assert_eval_eq("ClearAll[h]; h[Unevaluated[Plus[a, b]]]", "h[a + b]", 0);
}

/* HoldComplete prevents Unevaluated stripping (HoldAllComplete semantics) */
void test_holdcomplete_preserves_unevaluated() {
    assert_eval_eq("HoldComplete[Unevaluated[1+2]]", "HoldComplete[Unevaluated[1 + 2]]", 0);
    assert_eval_eq("HoldComplete[Unevaluated[Sequence[a, b]]]",
                   "HoldComplete[Unevaluated[Sequence[a, b]]]", 0);
}

/* ==============================================================
 * HoldComplete / HoldAllComplete evaluator semantics
 * ============================================================== */

/* Attributes */
void test_holdcomplete_attributes() {
    assert_eval_eq("Attributes[HoldComplete]", "{HoldAllComplete, Protected}", 0);
}

/* HoldComplete does not evaluate its arguments */
void test_holdcomplete_no_eval() {
    assert_eval_eq("HoldComplete[1+1]", "HoldComplete[1 + 1]", 0);
    /* Plus[1, 2, 3] is held but the printer still formats it as infix */
    assert_eval_eq("HoldComplete[Plus[1, 2, 3]]", "HoldComplete[1 + 2 + 3]", 0);
}

/* Evaluate inside HoldComplete is NOT honoured */
void test_holdcomplete_blocks_evaluate() {
    assert_eval_eq("HoldComplete[Evaluate[1+2]]", "HoldComplete[Evaluate[1 + 2]]", 0);
    assert_eval_eq("HoldComplete[1+1, Evaluate[1+2], Sequence[3, 4]]",
                   "HoldComplete[1 + 1, Evaluate[1 + 2], Sequence[3, 4]]", 0);
}

/* Sequence is NOT flattened inside HoldComplete */
void test_holdcomplete_blocks_sequence_flatten() {
    assert_eval_eq("HoldComplete[Sequence[a, b]]", "HoldComplete[Sequence[a, b]]", 0);
    assert_eval_eq("HoldComplete[Sequence[a, b], 1+2]",
                   "HoldComplete[Sequence[a, b], 1 + 2]", 0);
}

/* Substitution via ReplaceAll still happens inside HoldComplete */
void test_holdcomplete_replace_all() {
    assert_eval_eq("HoldComplete[f[1+2]] /. f[x_] :> g[x]",
                   "HoldComplete[g[1 + 2]]", 0);
    assert_eval_eq("HoldComplete[Sequence[a, b]] /. Sequence -> List",
                   "HoldComplete[{a, b}]", 0);
}

/* ReleaseHold removes exactly one layer of HoldComplete */
void test_holdcomplete_releasehold() {
    assert_eval_eq("ReleaseHold[HoldComplete[1+2]]", "3", 0);
    assert_eval_eq("ReleaseHold[HoldComplete[Sequence[1, 2]]]", "Sequence[1, 2]", 0);
}

/* ==============================================================
 * HoldPattern
 * ============================================================== */

void test_holdpattern_attributes() {
    assert_eval_eq("Attributes[HoldPattern]", "{HoldAll, Protected}", 0);
}

/* The HoldAll attribute keeps the pattern unevaluated; matching unwraps it. */
void test_holdpattern_basic_rule() {
    assert_eval_eq("a + b /. HoldPattern[_+_] -> 0", "0", 0);
    assert_eval_eq("{a + b, x * y} /. HoldPattern[_+_] -> 99", "{99, x y}", 0);
}

/* Cases with HoldPattern keeps the Rule structure intact. Without HoldPattern,
 * the pattern `a -> _` would be evaluated as a rule itself, not a match spec. */
void test_holdpattern_cases() {
    assert_eval_eq("Cases[{a -> b, c -> d}, HoldPattern[a -> _]]", "{a -> b}", 0);
}

void test_holdpattern_matchq() {
    assert_eval_eq("MatchQ[a + b, HoldPattern[_+_]]", "True", 0);
    assert_eval_eq("MatchQ[a, HoldPattern[_+_]]", "False", 0);
}

/* ==============================================================
 * Hold comparison (confirms HoldAll semantics still work)
 * ============================================================== */

void test_hold_basic() {
    assert_eval_eq("Hold[1+2]", "Hold[1 + 2]", 0);
    /* Hold has HoldAll, so Evaluate overrides */
    assert_eval_eq("Hold[Evaluate[1+2]]", "Hold[3]", 0);
    /* Hold DOES flatten Sequence (HoldAll is not HoldAllComplete) */
    assert_eval_eq("Hold[Sequence[a, b]]", "Hold[a, b]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_unevaluated_attributes);
    TEST(test_unevaluated_self);
    TEST(test_unevaluated_length_plus);
    TEST(test_unevaluated_sequence_not_flattened);
    TEST(test_unevaluated_stops_evaluate);
    TEST(test_unevaluated_in_holdall);
    TEST(test_unevaluated_mixed_hold);
    TEST(test_unevaluated_nested);
    TEST(test_unevaluated_symbolic_head);
    TEST(test_holdcomplete_preserves_unevaluated);

    TEST(test_holdcomplete_attributes);
    TEST(test_holdcomplete_no_eval);
    TEST(test_holdcomplete_blocks_evaluate);
    TEST(test_holdcomplete_blocks_sequence_flatten);
    TEST(test_holdcomplete_replace_all);
    TEST(test_holdcomplete_releasehold);

    TEST(test_holdpattern_attributes);
    TEST(test_holdpattern_basic_rule);
    TEST(test_holdpattern_cases);
    TEST(test_holdpattern_matchq);

    TEST(test_hold_basic);

    printf("All Unevaluated / HoldComplete tests passed!\n");
    return 0;
}
