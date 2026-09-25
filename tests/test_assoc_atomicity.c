/* Association atomicity and pattern semantics.
 *
 * In Mathematica an evaluated Association is an ATOM: AtomQ gives True, the
 * level-based functions see its values only (never the Rule wrappers or the
 * keys), and ReplaceAll never rewrites a key. Every expected output below was
 * produced by Mathematica 15.0 (wolframscript) unless the comment says the
 * behaviour is a deliberate Mathilda difference. The shared machinery lives in
 * src/assoc_struct.h. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------- AtomQ / AssociationQ / construction ---------- */

void test_atomq() {
    assert_eval_eq("AtomQ[<|a -> 1|>]", "True", 0);
    assert_eval_eq("AtomQ[<||>]", "True", 0);
    assert_eval_eq("AtomQ[<|a :> 1 + 1|>]", "True", 0);
    /* A malformed node is an ordinary expression, not an atom. */
    assert_eval_eq("AtomQ[Association[1, 2]]", "False", 0);
    assert_eval_eq("MatchQ[<|a -> 1|>, _?AtomQ]", "True", 0);
}

void test_associationq_wellformed_only() {
    assert_eval_eq("AssociationQ[<|a -> 1|>]", "True", 0);
    assert_eval_eq("AssociationQ[<|a :> 1|>]", "True", 0);
    assert_eval_eq("AssociationQ[<||>]", "True", 0);
    assert_eval_eq("AssociationQ[Association[1]]", "False", 0);
    assert_eval_eq("AssociationQ[Association[a -> 1, 2]]", "False", 0);
    assert_eval_eq("AssociationQ[Association[f[a -> 1]]]", "False", 0);
    assert_eval_eq("AssociationQ[{a -> 1}]", "False", 0);
    assert_eval_eq("Length[Association[1, 2]]", "2", 0);
}

/* ---------- Depth / LeafCount / Length ---------- */

void test_depth() {
    assert_eval_eq("Depth[<|a -> 1|>]", "2", 0);
    assert_eval_eq("Depth[<|a -> <|b -> 1|>|>]", "3", 0);
    assert_eval_eq("Depth[<||>]", "2", 0);
    assert_eval_eq("Depth[<|a -> {}|>]", "3", 0);
    assert_eval_eq("Depth[{<|a -> {1}|>}]", "4", 0);
    assert_eval_eq("Depth[<|a -> 1|>, Heads -> True]", "2", 0);
    /* The empty-compound rule of Depth[] itself. */
    assert_eval_eq("Depth[{}]", "2", 0);
    assert_eval_eq("Depth[{{}}]", "3", 0);
}

void test_leafcount_length() {
    assert_eval_eq("LeafCount[<|a -> 1|>]", "2", 0);
    assert_eval_eq("LeafCount[<|a -> 1, b -> f[x]|>]", "4", 0);
    assert_eval_eq("LeafCount[<|f[x] -> 1|>]", "2", 0);
    assert_eval_eq("LeafCount[<||>]", "1", 0);
    assert_eval_eq("LeafCount[<|a -> 1|>, Heads -> False]", "1", 0);
    assert_eval_eq("Length[<|a -> 1, b -> 2|>]", "2", 0);
}

/* ---------- Level ---------- */

void test_level() {
    assert_eval_eq("Level[<|a -> 1, b -> 2|>, {1}]", "{1, 2}", 0);
    assert_eval_eq("Level[<|a -> f[1], b -> 2|>, {2}]", "{1}", 0);
    assert_eval_eq("Level[<|a -> f[1], b -> 2|>, Infinity]", "{1, f[1], 2}", 0);
    assert_eval_eq("Level[{<|a -> f[1]|>}, Infinity]", "{1, f[1], <|a -> f[1]|>}", 0);
    assert_eval_eq("Level[<|a -> 1|>, {-1}]", "{1}", 0);
    assert_eval_eq("Level[<|a -> 1|>, {-2}]", "{<|a -> 1|>}", 0);
    assert_eval_eq("Level[<|a -> {1}|>, -1]", "{1, {1}}", 0);
    assert_eval_eq("Level[<|a -> 1|>, {1}, Heads -> True]", "{Association, 1}", 0);
    assert_eval_eq("Level[<|a -> 1|>, {0, 1}, Heads -> True]",
                   "{Association, 1, <|a -> 1|>}", 0);
    assert_eval_eq("Level[<|x -> 1|>, {0, Infinity}]", "{1, <|x -> 1|>}", 0);
}

/* ---------- FreeQ / OrderedQ / MemberQ / Count ---------- */

void test_freeq() {
    assert_eval_eq("FreeQ[<|a -> 1|>, a]", "True", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, 1]", "False", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, a -> 1]", "True", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, _Rule]", "True", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, Association]", "False", 0);
    assert_eval_eq("FreeQ[<|a -> 1|>, _Association]", "False", 0);
    assert_eval_eq("FreeQ[{<|a -> 1|>}, a]", "True", 0);
    assert_eval_eq("FreeQ[<|a -> x|>, x]", "False", 0);
    assert_eval_eq("FreeQ[{<|a -> 1|>}, 1, {1}]", "True", 0);
    assert_eval_eq("FreeQ[{<|a -> 1|>}, 1, {2}]", "False", 0);
}

void test_orderedq() {
    assert_eval_eq("OrderedQ[<|a -> 2, b -> 1|>]", "False", 0);
    assert_eval_eq("OrderedQ[<|b -> 1, a -> 2|>]", "True", 0);
    assert_eval_eq("OrderedQ[<||>]", "True", 0);
}

void test_memberq_count() {
    assert_eval_eq("MemberQ[<|a -> 1|>, 1]", "True", 0);
    assert_eval_eq("MemberQ[<|x -> 1|>, x]", "False", 0);
    assert_eval_eq("MemberQ[<|x -> 1|>, x -> 1]", "False", 0);
    assert_eval_eq("MemberQ[{<|a -> 1|>}, 1, Infinity]", "True", 0);
    assert_eval_eq("MemberQ[{<|a -> 1|>}, a, Infinity]", "False", 0);
    assert_eval_eq("Count[<|a -> 1, b -> 1|>, 1]", "2", 0);
    assert_eval_eq("Count[{<|a -> 1|>}, 1, Infinity]", "1", 0);
    assert_eval_eq("Count[{<|x -> 1|>}, x, Infinity]", "0", 0);
    assert_eval_eq("Count[<|a -> 1, b -> 1|>, 1, {0, Infinity}]", "2", 0);
}

/* ---------- Position / Cases / DeleteCases ---------- */

void test_position() {
    assert_eval_eq("Position[<|a -> 1, b -> 2|>, 2]", "{{Key[b]}}", 0);
    assert_eval_eq("Position[<|a -> 1, b -> 2|>, a]", "{}", 0);
    assert_eval_eq("Position[<|a -> 1, b -> {1}|>, 1]", "{{Key[a]}, {Key[b], 1}}", 0);
    assert_eval_eq("Position[{<|a -> 1|>}, 1]", "{{1, Key[a]}}", 0);
    assert_eval_eq("Position[{<|a -> <|b -> 1|>|>}, 1]", "{{1, Key[a], Key[b]}}", 0);
    assert_eval_eq("Position[<|a -> 1|>, _]", "{{0}, {Key[a]}, {}}", 0);
    assert_eval_eq("Position[<|a -> 1|>, _, Heads -> False]", "{{Key[a]}, {}}", 0);
    assert_eval_eq("Position[<|a -> 1|>, 1, {1}]", "{{Key[a]}}", 0);
    assert_eval_eq("FirstPosition[{<|a -> 1, b -> 2|>}, 2]", "{1, Key[b]}", 0);
}

void test_cases() {
    assert_eval_eq("Cases[<|a -> 1, b -> x|>, _Integer]", "{1}", 0);
    assert_eval_eq("Cases[<|a -> 1, b -> x|>, _Rule]", "{}", 0);
    assert_eval_eq("Cases[<|a -> {1}, b -> x|>, _Integer, Infinity]", "{1}", 0);
    assert_eval_eq("Cases[<|a -> {1}, b -> 2|>, _Integer, {2}]", "{1}", 0);
    assert_eval_eq("Cases[<|a -> 1|>, _, {0}]", "{<|a -> 1|>}", 0);
    assert_eval_eq("Cases[{<|a -> 1|>}, _Rule, Infinity]", "{}", 0);
    assert_eval_eq("Cases[{<|a -> 1|>}, _Integer, Infinity]", "{1}", 0);
    assert_eval_eq("Cases[<|a -> 1|>, _, Heads -> True]", "{Association, 1}", 0);
}

void test_delete_cases() {
    assert_eval_eq("DeleteCases[<|a -> 1, b -> 2|>, 1]", "<|b -> 2|>", 0);
    assert_eval_eq("DeleteCases[<|a -> 1, b -> 2|>, a -> 1]", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("DeleteCases[{<|a -> 1, b -> 2|>}, 1, Infinity]", "{<|b -> 2|>}", 0);
    assert_eval_eq("DeleteCases[{<|a -> 1, b -> 2|>}, a, Infinity]",
                   "{<|a -> 1, b -> 2|>}", 0);
    assert_eval_eq("DeleteCases[<|a -> {1, 2}, b -> 2|>, 1, {2}]", "<|a -> {2}, b -> 2|>", 0);
}

/* ---------- ReplaceAll / ReplaceRepeated / Replace ---------- */

void test_replace_all_keys_untouched() {
    assert_eval_eq("<|a -> 1, b -> 2|> /. b -> a", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("Keys[<|a -> 1, b -> 2|> /. b -> a]", "{a, b}", 0);
    assert_eval_eq("<|1 -> 2|> /. 1 -> 0", "<|1 -> 2|>", 0);
    assert_eval_eq("<|x -> x|> /. x -> 3", "<|x -> 3|>", 0);
    assert_eval_eq("<|x -> 1, y -> 2|> /. x -> z", "<|x -> 1, y -> 2|>", 0);
    assert_eval_eq("<|a -> <|b -> 1|>|> /. b -> c", "<|a -> <|b -> 1|>|>", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> //. b -> a", "<|a -> 1, b -> 2|>", 0);
}

void test_replace_all_entries_untouched() {
    assert_eval_eq("<|a -> 1, b -> 2|> /. (a -> 1) -> c", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("<|a -> 1|> /. (a -> x_) :> x", "<|a -> 1|>", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> /. (1 -> x_) :> x", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> /. Rule -> List", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> /. {a -> _} :> z", "<|a -> 1, b -> 2|>", 0);
}

void test_replace_all_values_and_whole() {
    assert_eval_eq("<|a -> 1, b -> 2|> /. 1 -> 5", "<|a -> 5, b -> 2|>", 0);
    assert_eval_eq("<|a -> 1, b -> <|c -> 2|>|> /. 2 -> 3", "<|a -> 1, b -> <|c -> 3|>|>", 0);
    assert_eval_eq("<|a -> 1, a -> 2|> /. 2 -> 3", "<|a -> 3|>", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> //. 1 -> 2", "<|a -> 2, b -> 2|>", 0);
    assert_eval_eq("<|a :> 1 + 1|> /. 1 -> 3", "<|a :> 3 + 3|>", 0);
    assert_eval_eq("<|a -> 1|> /. _Association -> z", "z", 0);
    assert_eval_eq("<|a -> 1|> /. <|a -> 1|> -> z", "z", 0);
    assert_eval_eq("{<|a -> 1|>} /. <|a -> x_|> :> x", "{1}", 0);
    assert_eval_eq("{<|a -> 1|>} /. _Association -> z", "{z}", 0);
    /* The head is still offered to the rules. */
    assert_eval_eq("<|a -> 1|> /. Association -> g", "g[a -> 1]", 0);
    assert_eval_eq("<|a -> 1|> /. Association -> List", "{a -> 1}", 0);
    assert_eval_eq("Head[<|a -> 1|> /. 1 -> Sequence[]]", "Association", 0);
    /* Deliberate difference: Mathematica does not re-evaluate the values of
     * a replaced association (st19 gives <|"a" -> 3, "b" -> 3^2|>). Mathilda
     * has no "already evaluated" mark on the rebuilt node, so the evaluator's
     * fixed point reduces the new values, the way it does for a List. */
    assert_eval_eq("<|\"a\" -> x, \"b\" -> x^2|> /. x -> 3", "<|\"a\" -> 3, \"b\" -> 9|>", 0);
}

void test_replace_levels() {
    assert_eval_eq("Replace[<|a -> 1, b -> 2|>, 1 -> 5]", "<|a -> 1, b -> 2|>", 0);
    assert_eval_eq("Replace[<|a -> 1, b -> 2|>, 1 -> 5, {1}]", "<|a -> 5, b -> 2|>", 0);
    assert_eval_eq("Replace[<|a -> 1, b -> 2|>, 1 -> 5, All]", "<|a -> 5, b -> 2|>", 0);
    assert_eval_eq("Replace[<|a -> {1}, b -> 2|>, 1 -> 5, {2}]", "<|a -> {5}, b -> 2|>", 0);
    assert_eval_eq("Replace[<|a -> 1|>, a -> b, {0, Infinity}]", "<|a -> 1|>", 0);
    assert_eval_eq("Replace[<|a -> 1|>, _ -> z, {-1}]", "<|a -> z|>", 0);
    assert_eval_eq("Replace[<|a -> {1}|>, {1} -> z, {-2}]", "<|a -> z|>", 0);
    assert_eval_eq("Replace[{<|a -> 1|>}, 1 -> 2, {2}]", "{<|a -> 2|>}", 0);
    assert_eval_eq("Replace[<|a -> 1|>, <|a -> 1|> -> 2, {0}]", "2", 0);
    assert_eval_eq("Replace[<|a -> 1, b -> 2|>, _Integer :> 0, {1}]", "<|a -> 0, b -> 0|>", 0);
    assert_eval_eq("Replace[<|a -> 1|>, Association -> z, {1}, Heads -> True]", "z[a -> 1]", 0);
}

/* ---------- An association used as a rule set ---------- */

void test_association_as_rules() {
    assert_eval_eq("ReplaceAll[{\"a\", \"b\"}, <|\"a\" -> 1|>]", "{1, \"b\"}", 0);
    assert_eval_eq("Replace[\"a\", <|\"a\" -> 1|>]", "1", 0);
    assert_eval_eq("ReplaceList[\"a\", <|\"a\" -> 1|>]", "{1}", 0);
    assert_eval_eq("ReplaceList[a, <|a -> 1|>, 0]", "{}", 0);
    assert_eval_eq("ReplaceRepeated[{\"a\", \"b\"}, <|\"a\" -> \"b\", \"b\" -> 1|>]", "{1, 1}", 0);
    assert_eval_eq("x /. <|x -> y, y -> z|>", "y", 0);
    assert_eval_eq("x //. <|x -> y, y -> z|>", "z", 0);
    assert_eval_eq("{a -> 1} /. <|a -> 2|>", "{2 -> 1}", 0);
    assert_eval_eq("Replace[{a, b}, <|a -> 1|>, {1}]", "{1, b}", 0);
    assert_eval_eq("{a, b} /. <||>", "{a, b}", 0);
    /* Keys match literally, never as patterns. */
    assert_eval_eq("x /. <|x_ -> 1|>", "x", 0);
    assert_eval_eq("3 /. <|_Integer -> 1|>", "3", 0);
    /* A list of associations is a list of rule sets. */
    assert_eval_eq("a /. {<|a -> 1|>, <|a -> 2|>}", "{1, 2}", 0);
    assert_eval_eq("Replace[a, {<|a -> 1|>, <|a -> 2|>}]", "{1, 2}", 0);
}

/* ---------- Map / MapAll / MapIndexed / Apply / Scan ---------- */

void test_map_levels() {
    assert_eval_eq("Map[f, <|a -> 1|>]", "<|a -> f[1]|>", 0);
    assert_eval_eq("Map[f, <|a -> 1|>, {1}]", "<|a -> f[1]|>", 0);
    assert_eval_eq("Map[f, <|a -> 1|>, {2}]", "<|a -> 1|>", 0);
    assert_eval_eq("Map[f, <|a -> {1}|>, {2}]", "<|a -> {f[1]}|>", 0);
    assert_eval_eq("Map[f, <|a -> {1}|>, {0, 2}]", "f[<|a -> f[{f[1]}]|>]", 0);
    assert_eval_eq("Map[f, <|a -> {1}|>, 2]", "<|a -> f[{f[1]}]|>", 0);
    assert_eval_eq("Map[f, <|a -> {1}|>, Infinity]", "<|a -> f[{f[1]}]|>", 0);
    assert_eval_eq("Map[f, {<|a -> 1|>}, {2}]", "{<|a -> f[1]|>}", 0);
    assert_eval_eq("Map[f, <|a -> 1|>, {-1}]", "<|a -> f[1]|>", 0);
    assert_eval_eq("Map[f, <|a -> g[1]|>, -1]", "<|a -> f[g[f[1]]]|>", 0);
    assert_eval_eq("Map[f, <|a -> 1, b -> 2|>, {0}]", "f[<|a -> 1, b -> 2|>]", 0);
    assert_eval_eq("Map[f, <|a -> 1|>, Heads -> True]", "f[Association][f[1]]", 0);
    /* RuleDelayed stays delayed. */
    assert_eval_eq("Map[f, <|a :> 1 + 1|>]", "<|a :> f[1 + 1]|>", 0);
    /* A malformed association maps over its arguments. */
    assert_eval_eq("Map[f, Association[1, 2], {1}]", "Association[f[1], f[2]]", 0);
    assert_eval_eq("AssociationQ[Map[f, <|a -> 1|>, {1}]]", "True", 0);
}

void test_mapall_mapindexed() {
    assert_eval_eq("MapAll[f, <|a -> {1}|>]", "f[<|a -> f[{f[1]}]|>]", 0);
    assert_eval_eq("MapIndexed[f, <|a -> 1, b -> 2|>]",
                   "<|a -> f[1, {Key[a]}], b -> f[2, {Key[b]}]|>", 0);
    assert_eval_eq("MapIndexed[f, {<|a -> 1|>}, {2}]", "{<|a -> f[1, {1, Key[a]}]|>}", 0);
    assert_eval_eq("MapIndexed[f, <|a -> 1|>, {0, 1}]", "f[<|a -> f[1, {Key[a]}]|>, {}]", 0);
}

void test_apply_levels() {
    assert_eval_eq("Apply[f, <|a -> 1|>]", "f[1]", 0);
    assert_eval_eq("Apply[f, <|a -> 1|>, {0}]", "f[1]", 0);
    assert_eval_eq("Apply[f, <|a -> 1|>, {1}]", "<|a -> 1|>", 0);
    assert_eval_eq("Apply[f, <|a -> {1}|>, {1}]", "<|a -> f[1]|>", 0);
    assert_eval_eq("Apply[List, <|a -> g[1]|>, {1}]", "<|a -> {1}|>", 0);
    assert_eval_eq("Apply[f, {<|a -> {1}|>}, {2}]", "{<|a -> f[1]|>}", 0);
    assert_eval_eq("Apply[f, {<|a -> {1}|>}, {1}]", "{f[{1}]}", 0);
    assert_eval_eq("Apply[f, {<|a -> 1|>}]", "f[<|a -> 1|>]", 0);
    assert_eval_eq("Apply[f, <|a -> {g[1]}|>, {2}]", "<|a -> {f[1]}|>", 0);
    /* A negative level is a depth bound: -1 is levels {1, -1}. */
    assert_eval_eq("Apply[f, <|a -> g[1]|>, -1]", "<|a -> f[1]|>", 0);
    assert_eval_eq("Apply[f, {g[1], h[g[2]]}, -1]", "{f[1], f[f[2]]}", 0);
    assert_eval_eq("Apply[f, <|a -> h[g[2]]|>, {-2}]", "<|a -> h[f[2]]|>", 0);
    assert_eval_eq("MapIndexed[f, <|a -> g[1]|>, -1]",
                   "<|a -> f[g[f[1, {Key[a], 1}]], {Key[a]}]|>", 0);
}

void test_scan_levels() {
    assert_eval_eq("Reap[Scan[Sow, <|a -> 1, b -> 2|>]]", "{Null, {{1, 2}}}", 0);
    assert_eval_eq("Reap[Scan[Sow, <|a -> {1}, b -> 2|>, {2}]]", "{Null, {{1}}}", 0);
    assert_eval_eq("Reap[Scan[Sow, <|a -> {1}, b -> 2|>, Infinity]]",
                   "{Null, {{1, {1}, 2}}}", 0);
    assert_eval_eq("Reap[Scan[Sow, {<|a -> {x}|>}, Infinity]]",
                   "{Null, {{x, {x}, <|a -> {x}|>}}}", 0);
    assert_eval_eq("Reap[Scan[Sow, <|a -> 1|>, Heads -> True]]",
                   "{Null, {{Association, 1}}}", 0);
}

/* ---------- Equal / Unequal ---------- */

void test_equal_unequal() {
    assert_eval_eq("<|a -> 1, b -> 2|> == <|b -> 2, a -> 1|>", "False", 0);
    assert_eval_eq("<|a -> 1|> == <|a -> 1.|>", "True", 0);
    assert_eval_eq("<|a -> 1|> == <|a -> 1|>", "True", 0);
    assert_eval_eq("<|a -> 1|> == <|a -> 2|>", "False", 0);
    assert_eval_eq("<|a -> x|> == <|b -> x|>", "False", 0);
    assert_eval_eq("<|a -> 1|> == <|a -> 1, b -> 2|>", "False", 0);
    assert_eval_eq("<|1 -> a|> == <|1. -> a|>", "False", 0);
    assert_eval_eq("<|a :> 1|> == <|a -> 1|>", "False", 0);
    assert_eval_eq("<||> == <||>", "True", 0);
    assert_eval_eq("<|a -> x|> == <|a -> y|>", "<|a -> x|> == <|a -> y|>", 0);
    assert_eval_eq("<|a -> 1|> == 1", "<|a -> 1|> == 1", 0);
    assert_eval_eq("<|a -> 1|> == <|a -> 1|> == <|a -> 1.|>", "True", 0);
    assert_eval_eq("Equal[<|a -> 2|>, <|a -> 1|>, <|a -> x|>]", "False", 0);
    assert_eval_eq("<|a -> 1|> != <|a -> 2|>", "True", 0);
    assert_eval_eq("<|a -> 1|> != <|a -> 1.|>", "False", 0);
    assert_eval_eq("<|a -> 1, b -> 2|> != <|b -> 2, a -> 1|>", "True", 0);
    assert_eval_eq("<|a -> 1|> != <|a -> x|>", "<|a -> 1|> != <|a -> x|>", 0);
    assert_eval_eq("Unequal[<|a -> 1|>, <|a -> 1|>, <|a -> 2|>]", "False", 0);
}

/* ---------- Pattern matching still sees associations ---------- */

static void define(const char* src) {
    struct Expr* p = parse_expression(src);
    struct Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

void test_patterns_still_match() {
    assert_eval_eq("Cases[{<|a -> 1|>, 2, <|b -> 2|>}, _Association]",
                   "{<|a -> 1|>, <|b -> 2|>}", 0);
    assert_eval_eq("Cases[{<|a -> 1|>, {<|b -> 2|>}}, _Association, Infinity]",
                   "{<|a -> 1|>, <|b -> 2|>}", 0);
    assert_eval_eq("Count[{<|a -> 1|>, <|b -> 2|>}, _Association]", "2", 0);
    assert_eval_eq("MatchQ[<|\"a\" -> 1|>, <|\"a\" -> x_|>]", "True", 0);
    assert_eval_eq("MatchQ[<|\"a\" -> 1|>, _Association]", "True", 0);
    assert_eval_eq("MatchQ[<|a -> 1, b -> 2|>, <|a -> _, b -> _|>]", "True", 0);
    assert_eval_eq("MatchQ[<|a -> 1, b -> 2|>, <|b -> _, a -> _|>]", "False", 0);
    assert_eval_eq("MatchQ[<|\"a\" -> 1, \"b\" -> 2|>, KeyValuePattern[\"a\" -> _]]", "True", 0);
    assert_eval_eq("Cases[{<|a -> 1|>}, KeyValuePattern[a -> _]]", "{<|a -> 1|>}", 0);
    assert_eval_eq("Cases[{<|\"a\" -> 1|>, 3, <|\"b\" -> 2|>}, <|\"a\" -> x_|> :> x]", "{1}", 0);
    define("s1ff[x_Association] := Keys[x]");
    assert_eval_eq("s1ff[<|q -> 1|>]", "{q}", 0);
    assert_eval_eq("s1ff[{q -> 1}]", "s1ff[{q -> 1}]", 0);
    define("s1gg[<|\"k\" -> v_|>] := v + 1");
    assert_eval_eq("s1gg[<|\"k\" -> 41|>]", "42", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_atomq);
    TEST(test_associationq_wellformed_only);
    TEST(test_depth);
    TEST(test_leafcount_length);
    TEST(test_level);
    TEST(test_freeq);
    TEST(test_orderedq);
    TEST(test_memberq_count);
    TEST(test_position);
    TEST(test_cases);
    TEST(test_delete_cases);
    TEST(test_replace_all_keys_untouched);
    TEST(test_replace_all_entries_untouched);
    TEST(test_replace_all_values_and_whole);
    TEST(test_replace_levels);
    TEST(test_association_as_rules);
    TEST(test_map_levels);
    TEST(test_mapall_mapindexed);
    TEST(test_apply_levels);
    TEST(test_scan_levels);
    TEST(test_equal_unequal);
    TEST(test_patterns_still_match);

    printf("All Association atomicity tests passed.\n");
    return 0;
}
