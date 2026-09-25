/* Tests for named slots (#name), Minus, generic operator (curried) forms,
 * Lookup's lazy default and KeyMemberQ/KeyFreeQ key patterns.
 *
 * Every expected value below was checked against Mathematica 15
 * (wolframscript); where Mathilda's printer spells a result differently
 * (e.g. `2 x` for `2*x`) only the spelling differs, not the value. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "opform.h"
#include <stdio.h>

static int g_checks = 0;

static void check(const char* in, const char* out) {
    g_checks++;
    assert_eval_eq(in, out, 0);
}

static void check_full(const char* in, const char* out) {
    g_checks++;
    assert_eval_eq(in, out, 1);
}

/* ---------------- Named slots ---------------- */

static void test_named_slot_parse(void) {
    check_full("Hold[#name]", "Hold[Slot[\"name\"]]");
    check_full("Hold[#a1b]", "Hold[Slot[\"a1b\"]]");
    check_full("Hold[#a$b]", "Hold[Slot[\"a$b\"]]");
    check_full("Hold[#$a]", "Hold[Slot[\"$a\"]]");
    check_full("Hold[#a`b]", "Hold[Slot[\"a`b\"]]");
    check_full("Hold[#\"quoted\"]", "Hold[Slot[\"quoted\"]]");
    check_full("Hold[#\"a b\"]", "Hold[Slot[\"a b\"]]");
    /* Unchanged forms. */
    check_full("Hold[#]", "Hold[Slot[1]]");
    check_full("Hold[#1]", "Hold[Slot[1]]");
    check_full("Hold[#12]", "Hold[Slot[12]]");
    check_full("Hold[##]", "Hold[SlotSequence[1]]");
    check_full("Hold[##2]", "Hold[SlotSequence[2]]");
    /* #1a is Slot[1]*a and ##a is SlotSequence[1]*a, as in WL. */
    check_full("Hold[#1a]", "Hold[Times[Slot[1], a]]");
    check_full("Hold[##a]", "Hold[Times[SlotSequence[1], a]]");
    /* Adjacent tokens. */
    check_full("Hold[#a_b]", "Hold[Times[Slot[\"a\"], Blank[b]]]");
    check_full("Hold[2#a]", "Hold[Times[2, Slot[\"a\"]]]");
    check_full("Hold[#a^2]", "Hold[Power[Slot[\"a\"], 2]]");
    check_full("Hold[#a[[1]]]", "Hold[Part[Slot[\"a\"], 1]]");
    check_full("Hold[x#]", "Hold[Times[x, Slot[1]]]");
    check_full("Hold[a#b]", "Hold[Times[a, Slot[\"b\"]]]");
    /* '#' inside strings and comments is not a slot. */
    check("StringLength[\"a#b\"]", "3");
    check("(* #a *) 5", "5");
}

static void test_named_slot_apply(void) {
    check("#a &[<|\"a\" -> 1|>]", "1");
    check("#a + #b &[<|\"a\" -> 1, \"b\" -> 2|>]", "3");
    check("Select[{<|\"a\" -> 1|>, <|\"a\" -> 3|>}, #a > 2 &]", "{<|\"a\" -> 3|>}");
    check("SortBy[{<|\"a\" -> 3|>, <|\"a\" -> 1|>}, #a &]", "{<|\"a\" -> 1|>, <|\"a\" -> 3|>}");
    check("#a & /@ {<|\"a\" -> 1|>, <|\"a\" -> 2|>}", "{1, 2}");
    check("#a[[1]] &[<|\"a\" -> {5, 6}|>]", "5");
    check("#\"a b\" &[<|\"a b\" -> 7|>]", "7");
    check("Slot[\"a\"] &[<|\"a\" -> 3|>]", "3");
    check("#a &[<|\"a\" :> 1 + 1|>]", "2");
    check("#a &[<|\"a\" -> 1|>, 7]", "1");
    check("#2 &[<|\"a\" -> 1|>, 7]", "7");
    check("{#a, #} &[<|\"a\" -> 1|>]", "{1, <|\"a\" -> 1|>}");
    check("2#a &[<|\"a\" -> 4|>]", "8");
    check("#a^2 &[<|\"a\" -> 4|>]", "16");
    /* A symbol key is not the string key "a". */
    check("#x &[<|x -> 1, \"x\" -> 2|>]", "2");
    /* Unfillable named slots stay in place (with a message), as in WL. */
    check("#a &[<|\"b\" -> 3|>]", "#a");
    check("#a + #b &[<|\"a\" -> 1|>]", "1 + #b");
    check("#a &[5]", "#a");
    check("#a &[]", "#a");
}

static void test_named_slot_print(void) {
    check("Slot[\"a\"]", "#a");
    check("Slot[\"a b\"]", "#\"a b\"");
    check("Slot[\"1a\"]", "#\"1a\"");
    check_full("Slot[\"a\"]", "Slot[\"a\"]");
    check("Hold[#a + #b &]", "Hold[#a + #b &]");
    check("ToExpression[ToString[Hold[#a + #b &], InputForm]] === Hold[#a + #b &]", "True");
}

/* ---------------- Minus ---------------- */

static void test_minus(void) {
    check("Minus[3]", "-3");
    check("Minus[x]", "-x");
    check_full("Minus[x]", "Times[-1, x]");
    check("Minus[{1, 2}]", "{-1, -2}");
    check("Minus[2.5]", "-2.5");
    check("Minus[1/2]", "-1/2");
    check("Minus[I]", "-I");
    check("Minus[0]", "0");
    check("Minus[Infinity]", "-Infinity");
    check("Minus[x, y]", "Minus[x, y]");
    check("Attributes[Minus]", "{Listable, NumericFunction, Protected}");
    check("SortBy[{3, 1, 2}, Minus]", "{3, 2, 1}");
    check("SortBy[<|\"a\" -> 3, \"b\" -> 1, \"c\" -> 2|>, Minus]",
          "<|\"a\" -> 3, \"c\" -> 2, \"b\" -> 1|>");
    check("KeySortBy[<|1 -> 1, 2 -> 2|>, Minus]", "<|2 -> 2, 1 -> 1|>");
    check("Minus[NDArray[{1., 2.}]]", "NDArray[{-1.0, -2.0}]");
}

/* ---------------- Operator forms ---------------- */

static void test_opform_assoc_heads(void) {
    const char* A = "a0 = <|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>;";
    char buf[512];
#define CA(expr, out) do { snprintf(buf, sizeof buf, "%s %s", A, expr); check(buf, out); } while (0)
    CA("Lookup[\"a\"][a0]", "1");
    CA("Lookup[{\"a\", \"c\"}][a0]", "{1, 3}");
    CA("Lookup[\"z\"][a0]", "Missing[\"KeyAbsent\", \"z\"]");
    CA("Lookup[\"a\"][{a0, a0}]", "{1, 1}");
    CA("KeyTake[{\"a\", \"b\"}][a0]", "<|\"a\" -> 1, \"b\" -> 2|>");
    CA("KeyTake[\"a\"][a0]", "<|\"a\" -> 1|>");
    CA("KeyTake[\"a\"][{a0, a0}]", "{<|\"a\" -> 1|>, <|\"a\" -> 1|>}");
    CA("KeyDrop[\"a\"][a0]", "<|\"b\" -> 2, \"c\" -> 3|>");
    CA("KeyDrop[{\"a\", \"b\"}][a0]", "<|\"c\" -> 3|>");
    CA("KeySelect[# != \"b\" &][a0]", "<|\"a\" -> 1, \"c\" -> 3|>");
    CA("KeyMap[f][a0]", "<|f[\"a\"] -> 1, f[\"b\"] -> 2, f[\"c\"] -> 3|>");
    CA("KeyExistsQ[\"a\"][a0]", "True");
    CA("KeyExistsQ[\"a\"][{\"a\" -> 1}]", "True");
    CA("KeyMemberQ[\"a\"][a0]", "True");
    CA("KeyFreeQ[\"a\"][a0]", "False");
    CA("KeyValueMap[f][a0]", "{f[\"a\", 1], f[\"b\", 2], f[\"c\", 3]}");
    CA("KeySortBy[Minus][<|1 -> 1, 2 -> 2|>]", "<|2 -> 2, 1 -> 1|>");
    CA("Select[OddQ][a0]", "<|\"a\" -> 1, \"c\" -> 3|>");
    CA("Select[# > 1 &][a0]", "<|\"b\" -> 2, \"c\" -> 3|>");
    CA("Apply[f][a0]", "f[1, 2, 3]");
    CA("Map[f][a0]", "<|\"a\" -> f[1], \"b\" -> f[2], \"c\" -> f[3]|>");
    CA("GroupBy[OddQ][a0]", "<|True -> <|\"a\" -> 1, \"c\" -> 3|>, False -> <|\"b\" -> 2|>|>");
    CA("Merge[f][{a0, a0}]", "<|\"a\" -> f[{1, 1}], \"b\" -> f[{2, 2}], \"c\" -> f[{3, 3}]|>");
    CA("AllTrue[Positive][a0]", "True");
    CA("SelectFirst[EvenQ][a0]", "2");
    CA("Append[\"d\" -> 4][a0]", "<|\"a\" -> 1, \"b\" -> 2, \"c\" -> 3, \"d\" -> 4|>");
    CA("Prepend[\"d\" -> 4][a0]", "<|\"d\" -> 4, \"a\" -> 1, \"b\" -> 2, \"c\" -> 3|>");
    CA("ReplacePart[\"a\" -> 9][a0]", "<|\"a\" -> 9, \"b\" -> 2, \"c\" -> 3|>");
    CA("Insert[\"z\" -> 5, 2][a0]", "<|\"a\" -> 1, \"z\" -> 5, \"b\" -> 2, \"c\" -> 3|>");
    CA("DeleteDuplicatesBy[OddQ][a0]", "<|\"a\" -> 1, \"b\" -> 2|>");
    CA("TakeLargestBy[Identity, 2][a0]", "<|\"c\" -> 3, \"b\" -> 2|>");
    CA("SortBy[Minus][a0]", "<|\"c\" -> 3, \"b\" -> 2, \"a\" -> 1|>");
#undef CA
}

static void test_opform_list_heads(void) {
    const char* L = "l0 = {3, 1, 4, 1, 5, 9, 2, 6};";
    char buf[512];
#define CL(expr, out) do { snprintf(buf, sizeof buf, "%s %s", L, expr); check(buf, out); } while (0)
    CL("Select[OddQ][l0]", "{3, 1, 1, 5, 9}");
    CL("Apply[f][{1, 2}]", "f[1, 2]");
    CL("Map[f][{1, 2}]", "{f[1], f[2]}");
    CL("Cases[_?OddQ][l0]", "{3, 1, 1, 5, 9}");
    CL("DeleteCases[1][l0]", "{3, 4, 5, 9, 2, 6}");
    CL("SortBy[Minus][l0]", "{9, 6, 5, 4, 3, 2, 1, 1}");
    CL("GroupBy[OddQ][l0]", "<|True -> {3, 1, 1, 5, 9}, False -> {4, 2, 6}|>");
    CL("CountsBy[OddQ][l0]", "<|True -> 5, False -> 3|>");
    CL("AssociationMap[f][{1, 2}]", "<|1 -> f[1], 2 -> f[2]|>");
    CL("Merge[Total][{<|1 -> 2|>, <|1 -> 3|>}]", "<|1 -> 5|>");
    CL("SelectFirst[EvenQ][l0]", "4");
    CL("AllTrue[EvenQ][l0]", "False");
    CL("AnyTrue[EvenQ][l0]", "True");
    CL("NoneTrue[EvenQ][l0]", "False");
    CL("Count[1][l0]", "2");
    CL("FreeQ[1][l0]", "False");
    CL("MemberQ[1][l0]", "True");
    CL("ReplaceAll[1 -> x][l0]", "{3, x, 4, x, 5, 9, 2, 6}");
    CL("Replace[1 -> x][1]", "x");
    CL("Append[x][l0]", "{3, 1, 4, 1, 5, 9, 2, 6, x}");
    CL("Prepend[x][l0]", "{x, 3, 1, 4, 1, 5, 9, 2, 6}");
    CL("Insert[x, 2][l0]", "{3, x, 1, 4, 1, 5, 9, 2, 6}");
    CL("Insert[x, {{1}, {2}}][l0]", "{x, 3, x, 1, 4, 1, 5, 9, 2, 6}");
    CL("Delete[2][l0]", "{3, 4, 1, 5, 9, 2, 6}");
    CL("ReplacePart[2 -> x][l0]", "{3, x, 4, 1, 5, 9, 2, 6}");
    CL("MapAt[f, 2][l0]", "{3, f[1], 4, 1, 5, 9, 2, 6}");
    CL("DeleteDuplicatesBy[OddQ][l0]", "{3, 4}");
    CL("MaximalBy[Minus][l0]", "{1, 1}");
    CL("MinimalBy[Minus][l0]", "{9}");
    CL("TakeLargestBy[Minus, 2][l0]", "{1, 1}");
    CL("TakeSmallestBy[Minus, 2][l0]", "{9, 6}");
    CL("TakeLargest[2][l0]", "{9, 6}");
    CL("TakeSmallest[2][l0]", "{1, 1}");
    CL("Extract[{1, 2}][{{1, 2}}]", "2");
    CL("Lookup[1 + 1][<|2 -> 3|>]", "3");
#undef CL
}

static void test_opform_negative(void) {
    /* Operator forms stay inert until applied. */
    check("Select[OddQ]", "Select[OddQ]");
    check("Lookup[1 + 1]", "Lookup[1 + 1]");
    /* Wrong operator arity / wrong number of data arguments: unevaluated. */
    check("Select[OddQ, 1][{1, 2, 3}]", "Select[OddQ, 1][{1, 2, 3}]");
    check("Lookup[\"z\", 0][<|\"a\" -> 1|>]", "Lookup[\"z\", 0][<|\"a\" -> 1|>]");
    check("Select[OddQ][{1}, {2}]", "Select[OddQ][{1}, {2}]");
    check("KeyTake[\"a\", 1][<|\"a\" -> 1|>]", "KeyTake[\"a\", 1][<|\"a\" -> 1|>]");
    /* Heads WL gives no operator form to are left alone. */
    check("Drop[1][{1, 2}]", "Drop[1][{1, 2}]");
    /* A form whose rewritten call does not evaluate keeps its curried shape. */
    check("KeyValueMap[f][{1}]", "KeyValueMap[f][{1}]");
    /* User-defined curried heads are untouched. */
    check("ff[1][2]", "ff[1][2]");
    /* Composes with Map / Query-style pipelines. */
    check("Map[Lookup[\"a\"], {<|\"a\" -> 1|>, <|\"a\" -> 2|>}]", "{1, 2}");
    check("Select[{<|\"a\" -> 1|>, <|\"a\" -> 3|>}, KeyMemberQ[\"a\"]]",
          "{<|\"a\" -> 1|>, <|\"a\" -> 3|>}");
}

static void test_opform_register_api(void) {
    /* The registry is open to other modules: register a Mathilda-side head. */
    opform_register("myHeadS4", 1, 1, 0);
    check("myHeadS4[x_, y_] := {x, y}; myHeadS4[2][1]", "{1, 2}");
    opform_register("myHeadS4b", 1, 2, 1);
    check("myHeadS4b[x_, y_Integer, z___] := {x, y, z}; {myHeadS4b[f][1], myHeadS4b[f, g][1]}", "{{f, 1}, {f, 1, g}}");
}

/* Operator forms of the second-tier heads in src/assoc_ops.c (Mathematica 15). */
static void test_opform_assoc_ops_heads(void) {
    check("Discard[EvenQ][{1, 2, 3}]", "{1, 3}");
    check("Discard[EvenQ][<|a -> 1, b -> 2|>]", "<|a -> 1|>");
    check("CountDistinctBy[Mod[#, 2] &][{1, 2, 3, 5}]", "2");
    check("AssociationComap[{f, g}][x]", "<|f -> f[x], g -> g[x]|>");
}

/* ---------------- Lookup laziness / KeyMemberQ patterns ---------------- */

static void test_lookup_lazy_default(void) {
    check("Attributes[Lookup]", "{HoldAll, Protected}");
    check("cnt = 0; Lookup[<|x -> 1|>, x, cnt++; 0]; cnt", "0");
    check("cnt = 0; Lookup[<|x -> 1|>, y, cnt++; 0]", "0");
    check("cnt", "1");
    check("cnt = 0; {Lookup[<|x -> 1|>, {x, y}, cnt++; 0], cnt}", "{{1, 0}, 1}");
    check("cnt = 0; {Lookup[{<|x -> 1|>, <|y -> 2|>}, x, cnt++; 0], cnt}", "{{1, 0}, 1}");
    check("cnt = 0; {Lookup[<|x -> 1|>, Key[x], cnt++; 0], cnt}", "{1, 0}");
    check("aa0 = <|1 -> 2|>; Lookup[aa0, 0 + 1]", "2");
    check("Lookup[{x -> 1}, z]", "Missing[\"KeyAbsent\", z]");
}

static void test_keymemberq_patterns(void) {
    check("KeyMemberQ[<|x -> 1, y -> 2|>, _]", "True");
    check("KeyFreeQ[<|x -> 1, y -> 2|>, _]", "False");
    check("KeyMemberQ[<|1 -> 1|>, _Integer]", "True");
    check("KeyFreeQ[<|x -> 1|>, _Integer]", "True");
    check("KeyMemberQ[<|\"a\" -> 1|>, Except[\"a\"]]", "False");
    check("KeyMemberQ[<|\"a\" -> 1, \"b\" -> 2|>, Except[\"a\"]]", "True");
    check("KeyMemberQ[{x -> 1}, _]", "True");
    check("KeyMemberQ[<|x -> 1|>, x | y]", "True");
    check("KeyMemberQ[<||>, _]", "False");
    /* KeyExistsQ does NOT treat its key as a pattern. */
    check("KeyExistsQ[<|x -> 1|>, _]", "False");
    check("KeyMemberQ[_][<|x -> 1|>]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_named_slot_parse);
    TEST(test_named_slot_apply);
    TEST(test_named_slot_print);
    TEST(test_minus);
    TEST(test_opform_assoc_heads);
    TEST(test_opform_list_heads);
    TEST(test_opform_negative);
    TEST(test_opform_register_api);
    TEST(test_opform_assoc_ops_heads);
    TEST(test_lookup_lazy_default);
    TEST(test_keymemberq_patterns);

    printf("All assoc-forms tests passed (%d assertions).\n", g_checks);
    return 0;
}
