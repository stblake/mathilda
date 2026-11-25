/*
 * test_abort.c -- Abort[], CheckAbort, AbortProtect and $Aborted.
 *
 * Abort[] is an in-flight sentinel (like Throw) that unwinds every enclosing
 * computation on the ordinary return path; only CheckAbort or the top level
 * (-> $Aborted) stop it.  Each control construct and every site that passes a
 * Throw upward is checked to pass an abort upward too, and the expected values
 * are the reference language's.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>

static void test_basic(void) {
    assert_eval_eq("{Attributes[Abort], Attributes[CheckAbort], Attributes[AbortProtect]}",
                   "{{Protected}, {HoldAll, Protected}, {HoldAll, Protected}}", 0);
    assert_eval_eq("Abort[]", "$Aborted", 0);                       /* top level */
    assert_eval_eq("CheckAbort[1 + 1, \"ab\"]", "2", 0);
    assert_eval_eq("CheckAbort[Abort[], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Abort[]; Print[\"never\"], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("Head[CheckAbort[Abort[], $Aborted]]", "Symbol", 0);
    assert_eval_eq("CheckAbort[f[Abort[]], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[1 + Abort[] + 2, \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[ReleaseHold[Hold[Abort[]]], \"ab\"]", "\"ab\"", 0);
    /* Held, it is just an expression. */
    assert_eval_eq("Hold[Abort[]]", "Hold[Abort[]]", 0);
    assert_eval_eq("Abort[1]", "Abort[1]", 0);
}

static void test_through_constructs(void) {
    assert_eval_eq("CheckAbort[Map[If[# == 2, Abort[], #] &, {1, 2, 3}], \"aborted\"]", "\"aborted\"", 0);
    assert_eval_eq("ab1[x_] := (If[x > 2, Abort[]]; x); CheckAbort[ab1 /@ {1, 2, 3, 4}, \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Table[If[i == 3, Abort[], i], {i, 5}], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Do[If[i == 3, Abort[]], {i, 5}], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[k = 0; While[True, k++; If[k == 4, Abort[]]], {\"ab\", k}]", "{\"ab\", 4}", 0);
    assert_eval_eq("CheckAbort[For[j = 0, j < 10, j++, If[j == 2, Abort[]]], {\"ab\", j}]", "{\"ab\", 2}", 0);
    assert_eval_eq("CheckAbort[Nest[If[# > 3, Abort[], # + 1] &, 0, 10], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Fold[If[#2 == 3, Abort[], #1 + #2] &, 0, Range[5]], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[NestWhile[If[# > 3, Abort[], # + 1] &, 0, True &], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[FixedPoint[If[# > 3, Abort[], # + 1] &, 0], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Select[{1, 2, 3}, If[# == 2, Abort[], True] &], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Scan[If[# == 2, Abort[]] &, {1, 2, 3}], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Which[Abort[], 1, True, 2], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Switch[1, 1, Abort[]], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Reap[Sow[1]; Abort[]; Sow[2]], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Module[{t}, Abort[]], \"ab\"]", "\"ab\"", 0);
    /* Scoping is unwound before failexpr runs. */
    assert_eval_eq("x = 1; CheckAbort[Block[{x = 2}, Abort[]; x], {\"ab\", x}]", "{\"ab\", 1}", 0);
    assert_eval_eq("ab2[0] := Abort[]; ab2[n_] := ab2[n - 1] + 1; CheckAbort[ab2[50], \"deep\"]", "\"deep\"", 0);
}

static void test_not_caught_by_others(void) {
    /* Only CheckAbort catches an abort ... */
    assert_eval_eq("CheckAbort[Catch[Abort[]; 5], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Catch[Abort[], _], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Check[Abort[], \"msg\"], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[Quiet[Abort[]], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[TimeConstrained[Abort[], 5, \"tc\"], \"ab\"]", "\"ab\"", 0);
    /* ... and CheckAbort does not catch a Throw. */
    assert_eval_eq("Catch[CheckAbort[Throw[7], \"ab\"]]", "7", 0);
    assert_eval_eq("Catch[CheckAbort[Throw[7, t], \"ab\"], t]", "7", 0);
}

static void test_nesting(void) {
    assert_eval_eq("CheckAbort[CheckAbort[Abort[], \"inner\"]; 9, \"outer\"]", "9", 0);
    assert_eval_eq("CheckAbort[CheckAbort[Abort[], Abort[]], \"outer\"]", "\"outer\"", 0);
    assert_eval_eq("CheckAbort[{CheckAbort[Abort[], 1], CheckAbort[2, 3]}, \"outer\"]", "{1, 2}", 0);
}

static void test_abortprotect(void) {
    assert_eval_eq("AbortProtect[1 + 1]", "2", 0);
    assert_eval_eq("CheckAbort[AbortProtect[Abort[]; \"after\"], \"ab\"]", "\"ab\"", 0);
    /* The protected body runs to the end; the abort fires after it. */
    assert_eval_eq("ap = 0; CheckAbort[AbortProtect[Abort[]; ap = 1]; ap = 2, {\"ab\", ap}]", "{\"ab\", 1}", 0);
    /* Nested: deferred to the outermost AbortProtect. */
    assert_eval_eq("aq = {}; CheckAbort[AbortProtect[AbortProtect[Abort[]; AppendTo[aq, 1]]; AppendTo[aq, 2]]; "
                   "AppendTo[aq, 3], aq]", "{1, 2}", 0);
    /* No abort requested: nothing fires afterwards. */
    assert_eval_eq("CheckAbort[AbortProtect[5]; 6, \"ab\"]", "6", 0);
    assert_eval_eq("AbortProtect[Abort[]]", "$Aborted", 0);
}

/* CheckAbort's failexpr is abort-protected (the reference: an abort raised in it
 * waits until it ends, then propagates), and PropagateAborts -> True passes the
 * caught abort on after failexpr. */
static void test_checkabort_failexpr(void) {
    assert_eval_eq("cf = {}; r = CheckAbort[CheckAbort[Abort[], Abort[]; AppendTo[cf, \"after\"]], \"o\"]; {r, cf}",
                   "{\"o\", {\"after\"}}", 0);
    assert_eval_eq("CheckAbort[CheckAbort[Abort[], \"in\", PropagateAborts -> True], \"out\"]", "\"out\"", 0);
    assert_eval_eq("CheckAbort[CheckAbort[Abort[], \"in\", PropagateAborts -> False], \"out\"]", "\"in\"", 0);
    assert_eval_eq("CheckAbort[CheckAbort[Abort[], \"in\", PropagateAborts -> Automatic], \"out\"]", "\"in\"", 0);
    assert_eval_eq("cf = {}; CheckAbort[CheckAbort[Abort[], AppendTo[cf, 1], PropagateAborts -> True]; "
                   "AppendTo[cf, 2], cf]", "{1}", 0);
    /* No abort: PropagateAborts changes nothing. */
    assert_eval_eq("CheckAbort[5, 6, PropagateAborts -> True]", "5", 0);
    assert_eval_eq("Options[CheckAbort]", "{PropagateAborts -> Automatic}", 0);
    /* A malformed third argument leaves the call alone. */
    assert_eval_eq("CheckAbort[5, 6, 7]", "CheckAbort[5, 6, 7]", 0);
    assert_eval_eq("CheckAbort[5, 6, PropagateAborts -> 3]", "CheckAbort[5, 6, PropagateAborts -> 3]", 0);
}

static void test_get_and_files(void) {
    char path[128], src[512];
    snprintf(path, sizeof path, "/tmp/Mathilda_abort_%d.m", (int)getpid());
    FILE* fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("abA = 1;\nAbort[];\nabB = 2;\n", fp);
    fclose(fp);
    /* The abort stops the read and leaves Get: the statement around Get aborts. */
    snprintf(src, sizeof src, "r = CheckAbort[Get[\"%s\"], \"caught\"]; {r, abA, abB}", path);
    assert_eval_eq(src, "{\"caught\", 1, abB}", 0);
    snprintf(src, sizeof src, "Get[\"%s\"]", path);
    assert_eval_eq(src, "$Aborted", 0);
    unlink(path);
    /* ToExpression's statements stop at the abort too. */
    assert_eval_eq("CheckAbort[ToExpression[\"abC = 1\\nAbort[]\\nabD = 2\"], \"caught\"]; {abC, abD}",
                   "{1, abD}", 0);
}

/* A Throw or Abort raised inside a user callback -- a predicate, comparator,
 * key function, pattern test or guard -- leaves the builtin that called it,
 * instead of being read as "not True" and dropped.  Each row is run with U
 * replaced by Throw[99] (expect 99 from Catch) and by Abort[] (expect "ab" from
 * CheckAbort). */
static void test_unwind_through_callbacks(void) {
    static const char* const rows[][2] = {
        { "Select", "Select[{1, 2, 3}, If[# == 2, U, True] &]" },
        { "Count", "Count[{1, 2, 3}, _?(If[# == 2, U, True] &)]" },
        { "Cases", "Cases[{1, 2, 3}, _?(If[# == 2, U, True] &)]" },
        { "DeleteCases", "DeleteCases[{1, 2, 3}, _?(If[# == 2, U, True] &)]" },
        { "Position", "Position[{1, 2, 3}, _?(If[# == 2, U, True] &)]" },
        { "MemberQ", "MemberQ[{3, 2, 1}, _?(If[# == 2, U, False] &)]" },
        { "FreeQ", "FreeQ[{3, 2, 1}, _?(If[# == 2, U, False] &)]" },
        { "MatchQ", "MatchQ[2, _?(If[# == 2, U, True] &)]" },
        { "Condition", "MatchQ[2, x_ /; If[x == 2, U, True]]" },
        { "ReplaceAll", "{1, 2, 3} /. x_ /; If[x == 2, U, True] :> 0" },
        { "DownValue guard", "ClearAll[cbg]; cbg[x_] /; If[x == 2, U, True] := x; cbg[2]" },
        { "Sort", "Sort[{3, 2, 1}, If[#1 == 2 || #2 == 2, U, True] &]" },
        { "SortBy", "SortBy[{3, 2, 1}, If[# == 2, U, #] &]" },
        { "GatherBy", "GatherBy[{1, 2, 3}, If[# == 2, U, #] &]" },
        { "SplitBy", "SplitBy[{1, 2, 3}, If[# == 2, U, #] &]" },
        { "Split", "Split[{1, 2, 3}, If[#1 == 2 || #2 == 2, U, True] &]" },
        { "TakeWhile", "TakeWhile[{1, 2, 3}, If[# == 2, U, True] &]" },
        { "LengthWhile", "LengthWhile[{1, 2, 3}, If[# == 2, U, True] &]" },
        { "DeleteDuplicates", "DeleteDuplicates[{1, 2, 3}, If[#1 == 2 || #2 == 2, U, True] &]" },
        { "Union", "Union[{1, 2, 3}, SameTest -> (If[#1 == 2 || #2 == 2, U, True] &)]" },
        { "MaximalBy", "MaximalBy[{1, 2, 3}, If[# == 2, U, #] &]" },
        { "MinimalBy", "MinimalBy[{1, 2, 3}, If[# == 2, U, #] &]" },
        { "Tally", "Tally[{1, 2, 3}, If[#1 == 2 || #2 == 2, U, True] &]" },
        { "KeySelect", "KeySelect[<|1 -> a, 2 -> b|>, If[# == 2, U, True] &]" },
        { "GroupBy", "GroupBy[{1, 2, 3}, If[# == 2, U, #] &]" },
    };
    char src[1024], expr[768];
    for (size_t i = 0; i < sizeof rows / sizeof rows[0]; i++) {
        const char* t = rows[i][1];
        const char* u = strstr(t, "U");
        ASSERT(u != NULL);
        /* Throw */
        snprintf(expr, sizeof expr, "%.*sThrow[99]%s", (int)(u - t), t, u + 1);
        snprintf(src, sizeof src, "Catch[%s]", expr);
        printf("  [%s] throw\n", rows[i][0]);
        assert_eval_eq(src, "99", 0);
        /* Abort */
        snprintf(expr, sizeof expr, "%.*sAbort[]%s", (int)(u - t), t, u + 1);
        snprintf(src, sizeof src, "CheckAbort[%s, \"ab\"]", expr);
        printf("  [%s] abort\n", rows[i][0]);
        assert_eval_eq(src, "\"ab\"", 0);
    }
    /* A caught Throw inside a callback is consumed, not re-raised. */
    assert_eval_eq("Select[{1, 2, 3}, Catch[If[# == 2, Throw[False], True]] &]", "{1, 3}", 0);
    assert_eval_eq("Cases[{1, 2, 3}, _?(CheckAbort[If[# == 2, Abort[], True], False] &)]", "{1, 3}", 0);
    /* Speculative internal evaluation is not a callback: Sum's symbolic probe of
     * the summand does not leak the Throw that its numeric loop never reaches. */
    assert_eval_eq("Catch[Sum[If[i == 5, Throw[i], i], {i, 3}]]", "6", 0);
    assert_eval_eq("Catch[Sum[If[i == 5, Throw[i], i], {i, 10}]]", "5", 0);
}

/* TimeConstrained's timeout is the ordinary Abort[] (cooperative: it unwinds
 * through every frame's cleanup) owned by its TimeConstrained.  Each expected
 * value is the reference kernel's. */
static void test_timeconstrained_abort(void) {
    /* CheckAbort runs failexpr, and (Automatic) the timeout then goes on. */
    assert_eval_eq("TimeConstrained[CheckAbort[While[True], \"ca\"], 0.2, \"tc\"]", "\"tc\"", 0);
    assert_eval_eq("TimeConstrained[CheckAbort[While[True], \"ca\"]; \"after\", 0.2, \"tc\"]", "\"tc\"", 0);
    assert_eval_eq("tl = {}; {TimeConstrained[CheckAbort[While[True], AppendTo[tl, \"f\"]; \"f\", "
                   "PropagateAborts -> True], 0.2, \"tc\"], tl}", "{\"tc\", {\"f\"}}", 0);
    /* PropagateAborts -> False absorbs it, and the frame never fires again. */
    assert_eval_eq("tl = {}; {TimeConstrained[CheckAbort[While[True], AppendTo[tl, \"f\"], "
                   "PropagateAborts -> False]; AppendTo[tl, \"after\"]; \"body\", 0.2, \"tc\"], tl}",
                   "{\"body\", {\"f\", \"after\"}}", 0);
    assert_eval_eq("{TimeConstrained[CheckAbort[While[True], 0, PropagateAborts -> False]; "
                   "t0 = AbsoluteTime[]; While[AbsoluteTime[] - t0 < 0.5]; \"spun\", 0.2, \"tc\"]}",
                   "{\"spun\"}", 0);
    /* failexpr of a CheckAbort inside runs to completion. */
    assert_eval_eq("tl = {}; {TimeConstrained[CheckAbort[While[True], t0 = AbsoluteTime[]; "
                   "While[AbsoluteTime[] - t0 < 0.4]; AppendTo[tl, \"f\"]], 0.2, \"tc\"], tl}",
                   "{\"tc\", {\"f\"}}", 0);
    /* AbortProtect inside the body runs to completion; an enclosing one does not
     * stop the timeout, and the code around it continues. */
    assert_eval_eq("tl = {}; {TimeConstrained[AbortProtect[t0 = AbsoluteTime[]; While[AbsoluteTime[] - t0 < 0.4]; "
                   "AppendTo[tl, 1]]; AppendTo[tl, 2], 0.2, \"tc\"], tl}", "{\"tc\", {1}}", 0);
    assert_eval_eq("tl = {}; {AbortProtect[TimeConstrained[While[True], 0.2, AppendTo[tl, \"tc\"]; \"tc\"]; "
                   "AppendTo[tl, \"past\"]], tl}", "{{\"tc\", \"past\"}, {\"tc\", \"past\"}}", 0);
    assert_eval_eq("TimeConstrained[AbortProtect[Pause[0.4]; \"p\"], 0.2, \"tc\"]", "\"tc\"", 0);
    /* Not caught outside, by Catch or by Check; a plain Abort[] is not consumed. */
    assert_eval_eq("CheckAbort[TimeConstrained[While[True], 0.2], \"outer\"]", "$Aborted", 0);
    assert_eval_eq("TimeConstrained[Catch[While[True]], 0.2, \"tc\"]", "\"tc\"", 0);
    assert_eval_eq("TimeConstrained[Check[While[True], \"chk\"], 0.2, \"tc\"]", "\"tc\"", 0);
    assert_eval_eq("CheckAbort[TimeConstrained[Abort[], 5, \"tc\"], \"ca\"]", "\"ca\"", 0);
    /* Nesting: the expired budget answers. */
    assert_eval_eq("TimeConstrained[TimeConstrained[While[True], 10, \"inner\"], 0.2, \"outer\"]", "\"outer\"", 0);
    assert_eval_eq("TimeConstrained[TimeConstrained[While[True], 0.2, \"inner\"], 10, \"outer\"]", "\"inner\"", 0);
    assert_eval_eq("Table[TimeConstrained[While[True], 0.05, k], {k, 3}]", "{1, 2, 3}", 0);
    /* Blocking and scoped bodies. */
    assert_eval_eq("TimeConstrained[Pause[5], 0.2, \"tc\"]", "\"tc\"", 0);
    assert_eval_eq("x = 1; {TimeConstrained[Block[{x = 2}, While[True]], 0.2, \"tc\"], x}", "{\"tc\", 1}", 0);
    assert_eval_eq("TimeConstrained[Module[{s = 0}, While[True, s++]], 0.2, \"tc\"]", "\"tc\"", 0);
    /* Afterwards nothing is left armed. */
    assert_eval_eq("TimeConstrained[While[True], 0.1, 0]; t0 = AbsoluteTime[]; Pause[0.3]; "
                   "AbsoluteTime[] - t0 >= 0.29", "True", 0);
    assert_eval_eq("CheckAbort[Abort[], \"ab\"]", "\"ab\"", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_basic);
    TEST(test_through_constructs);
    TEST(test_not_caught_by_others);
    TEST(test_nesting);
    TEST(test_abortprotect);
    TEST(test_checkabort_failexpr);
    TEST(test_timeconstrained_abort);
    TEST(test_get_and_files);
    TEST(test_unwind_through_callbacks);
    printf("All abort tests passed.\n");
    return 0;
}
