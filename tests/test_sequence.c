#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "attr.h"

/* ==============================================================
 * Sequence: automatic argument splicing
 * ============================================================== */

/* Sequence[...] splices into an enclosing function's argument list. */
void test_sequence_basic_splice() {
    assert_eval_eq("f[a, Sequence[b, c], d]", "f[a, b, c, d]", 0);
    assert_eval_eq("g[Sequence[1, 2], Sequence[3, 4]]", "g[1, 2, 3, 4]", 0);
    /* Works in list literals too (List is an ordinary head). */
    assert_eval_eq("{a, Sequence[b, c], Sequence[d, e]}", "{a, b, c, d, e}", 0);
}

/* A single-argument Sequence acts like the identity (regression: a lone
 * Sequence[x] changes structure without changing arg_count, so the
 * flatten_sequences fast-path must not skip it). */
void test_sequence_single_arg() {
    assert_eval_eq("{a, Sequence[b], c, Identity[d]}", "{a, b, c, d}", 0);
    assert_eval_eq("f[Sequence[x]]", "f[x]", 0);
    assert_eval_eq("f[a, Sequence[b], c]", "f[a, b, c]", 0);
}

/* An empty Sequence[] evaporates. */
void test_sequence_empty() {
    assert_eval_eq("{Sequence[], a}", "{a}", 0);
    assert_eval_eq("f[Sequence[], a]", "f[a]", 0);
    assert_eval_eq("{Sequence[]}", "{}", 0);
    assert_eval_eq("{a, Sequence[], b}", "{a, b}", 0);
}

/* A bare Sequence (no enclosing function) stays as a Sequence object, and so
 * does one stored in an OwnValue -- splicing only happens at a call site. */
void test_sequence_bare_and_ownvalue() {
    assert_eval_eq("Sequence[a, b, c]", "Sequence[a, b, c]", 0);
    assert_eval_eq("uu = Sequence[a, b, c]; uu", "Sequence[a, b, c]", 0);
    assert_eval_eq("uu = Sequence[a, b, c]; {uu, uu, uu}",
                   "{a, b, c, a, b, c, a, b, c}", 0);
}

/* Replacing with a Sequence splices at each replaced position. */
void test_sequence_replace_splice() {
    assert_eval_eq("{u, u, u} /. u -> Sequence[a, b, c]",
                   "{a, b, c, a, b, c, a, b, c}", 0);
    assert_eval_eq("{a, b, f[x, y], g[w], f[z, y]} /. f -> Sequence",
                   "{a, b, x, y, g[w], z, y}", 0);
    /* List -> Sequence completely flattens one level of sublists. */
    assert_eval_eq("f[{{a, b}, {c, d}, {a}}] /. List -> Sequence",
                   "f[a, b, c, d, a]", 0);
}

/* A __ (BlankSequence) match is a Sequence object when substituted alone. */
void test_sequence_from_blanksequence() {
    assert_eval_eq("f[a, b, c] /. f[x__] -> x", "Sequence[a, b, c]", 1);
    /* Spliced into an enclosing head on the RHS. */
    assert_eval_eq("f[a, b, c] /. f[x__] -> gg[x]", "gg[a, b, c]", 0);
}

/* SlotSequence (##) produces a Sequence object. */
void test_sequence_slotsequence() {
    assert_eval_eq("##&[a, b, c]", "Sequence[a, b, c]", 1);
    assert_eval_eq("f[##2]&[a, b, c, d]", "f[b, c, d]", 0);
}

/* ==============================================================
 * SequenceHold / HoldAllComplete gating
 * ============================================================== */

/* HoldAll (without SequenceHold) does NOT prevent splicing. */
void test_sequence_holdall_still_splices() {
    assert_eval_eq("Hold[Sequence[a, b]]", "Hold[a, b]", 0);
    assert_eval_eq("Hold[Sequence[]]", "Hold[]", 0);
    assert_eval_eq("Hold[a, Sequence[b, c]]", "Hold[a, b, c]", 0);
    assert_eval_eq("Attributes[Hold]", "{HoldAll, Protected}", 0);
}

/* HoldAllComplete implies SequenceHold: Sequence is preserved. */
void test_sequence_holdallcomplete_preserves() {
    assert_eval_eq("HoldComplete[Sequence[]]", "HoldComplete[Sequence[]]", 0);
    assert_eval_eq("HoldComplete[Sequence[a, b]]", "HoldComplete[Sequence[a, b]]", 0);
    assert_eval_eq("HoldComplete[a, Sequence[b, c]]",
                   "HoldComplete[a, Sequence[b, c]]", 0);
    assert_eval_eq("Attributes[HoldComplete]", "{HoldAllComplete, Protected}", 0);
}

/* A user head with the SequenceHold attribute suppresses splicing -- this is
 * the point of the attribute-driven gate. Use fresh symbols per assertion so
 * attributes set in the shared symbol table don't leak between checks. */
void test_sequence_user_sequencehold() {
    assert_eval_eq("SetAttributes[shA, SequenceHold]; shA[a, Sequence[b, c]]",
                   "shA[a, Sequence[b, c]]", 0);
    /* HoldAll alone is not enough; SequenceHold must also be present. */
    assert_eval_eq("SetAttributes[shB, {HoldAll, SequenceHold}]; shB[a, Sequence[b, c]]",
                   "shB[a, Sequence[b, c]]", 0);
    /* Without SequenceHold, a HoldAll user head still splices. */
    assert_eval_eq("SetAttributes[shC, HoldAll]; shC[a, Sequence[b, c]]",
                   "shC[a, b, c]", 0);
}

/* Assignment and rule heads are SequenceHold, so a Sequence can be returned as
 * a result and spliced only at the eventual call site. */
void test_sequence_assignment_returns_sequence() {
    assert_eval_eq("splice[x_] := Sequence[x, x, x]; {a, splice[b], c}",
                   "{a, b, b, b, c}", 0);
    assert_eval_eq("{f[1], g[2], hh[3]} /. g[x_] :> Sequence[x, x]",
                   "{f[1], 2, 2, hh[3]}", 0);
    /* Rule/RuleDelayed keep the literal Sequence on their RHS. */
    assert_eval_eq("u -> Sequence[a, b]", "u -> Sequence[a, b]", 0);
    assert_eval_eq("u :> Sequence[a, b]", "u :> Sequence[a, b]", 0);
}

/* Ordinary functions with special input forms still splice. */
void test_sequence_special_input_forms() {
    assert_eval_eq("a == Sequence[b, c]", "a == b == c", 0);
}

/* Sequence directly inside Unevaluated is preserved (Unevaluated is
 * HoldAllComplete); regression guard for the splice/Unevaluated ordering. */
void test_sequence_unevaluated_preserved() {
    assert_eval_eq("Length[Unevaluated[Sequence[a, b]]]", "2", 0);
    assert_eval_eq("Length[Unevaluated[Sequence[a, b, c, d]]]", "4", 0);
}

/* ==============================================================
 * Attributes and metadata
 * ============================================================== */

void test_sequence_attributes() {
    assert_eval_eq("Attributes[Sequence]", "{Protected}", 0);
    assert_eval_eq("Attributes[SequenceHold]", "{Protected}", 0);
    /* Assignment/replacement heads now carry SequenceHold. */
    assert_eval_eq("Attributes[Set]", "{HoldFirst, Protected, SequenceHold}", 0);
    assert_eval_eq("Attributes[SetDelayed]", "{HoldAll, Protected, SequenceHold}", 0);
    assert_eval_eq("Attributes[Rule]", "{Protected, SequenceHold}", 0);
    assert_eval_eq("Attributes[RuleDelayed]", "{HoldRest, Protected, SequenceHold}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_sequence_basic_splice);
    TEST(test_sequence_single_arg);
    TEST(test_sequence_empty);
    TEST(test_sequence_bare_and_ownvalue);
    TEST(test_sequence_replace_splice);
    TEST(test_sequence_from_blanksequence);
    TEST(test_sequence_slotsequence);

    TEST(test_sequence_holdall_still_splices);
    TEST(test_sequence_holdallcomplete_preserves);
    TEST(test_sequence_user_sequencehold);
    TEST(test_sequence_assignment_returns_sequence);
    TEST(test_sequence_special_input_forms);
    TEST(test_sequence_unevaluated_preserved);

    TEST(test_sequence_attributes);

    printf("All Sequence tests passed!\n");
    return 0;
}
