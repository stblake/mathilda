/* Tests for in-place Association mutation (stream S2 of the Association
 * overhaul): a[k] = v, the compound assignments (+=, -=, *=, /=, ++, --) on
 * association elements, a[k] =., KeyDropFrom, Delete[a, Key[k]], nested and
 * curried element assignment, AssociateTo, RuleDelayed preservation, the
 * index/scan consistency on malformed nodes -- and, above all, aliasing: an
 * in-place write must never be visible through another holder of the same
 * association (b = a; a[k] = v leaves b alone), in every mutation form.
 *
 * Every expected value below was checked against Mathematica 15.
 *
 * Each CHECK drives the whole pipeline (parse -> evaluate -> print). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "assoc.h"
#include "assoc_index.h"
#include <stdio.h>

static int g_checks = 0;

#define CHECK(in, out) do { assert_eval_eq((in), (out), 0); g_checks++; } while (0)
#define CHECK_FF(in, out) do { assert_eval_eq((in), (out), 1); g_checks++; } while (0)
#define CHECK_TRUE(cond) do { if (!(cond)) { \
        fprintf(stderr, "FAIL: %s (line %d)\n", #cond, __LINE__); abort(); } \
        g_checks++; } while (0)

static void run(const char* src) {
    struct Expr* p = parse_expression(src);
    assert(p != NULL);
    struct Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* The association stored as `name`'s OwnValue (borrowed), or NULL. */
static Expr* stored(const char* name) {
    Expr** slot = assoc_symbol_slot(name);
    return slot ? *slot : NULL;
}

/* ---------- compound assignment on elements ---------- */

void test_compound_ops_keyed() {
    run("a = <|\"x\" -> 1, \"y\" -> 2|>");
    CHECK("a[\"x\"] += 5", "6");
    CHECK("a", "<|\"x\" -> 6, \"y\" -> 2|>");
    CHECK("a[\"x\"] -= 1", "5");
    CHECK("a[\"x\"] *= 3", "15");
    CHECK("a[\"x\"] /= 2", "15/2");
    CHECK("a[\"x\"]++", "15/2");
    CHECK("a[\"x\"]", "17/2");
    CHECK("a[\"x\"]--", "17/2");
    CHECK("++a[\"x\"]", "17/2");
    CHECK("--a[\"x\"]", "15/2");
    CHECK("a", "<|\"x\" -> 15/2, \"y\" -> 2|>");
}

void test_compound_ops_part_forms() {
    run("a = <|\"x\" -> 1, \"y\" -> 2|>");
    CHECK("a[[Key[\"y\"]]] += 10", "12");
    CHECK("a[[\"y\"]] += 10", "22");
    CHECK("a[[1]] += 100", "101");
    CHECK("a[[Key[\"y\"]]]++", "22");
    CHECK("a", "<|\"x\" -> 101, \"y\" -> 23|>");
    CHECK("a[[\"x\"]] *= 10; a[[\"x\"]] /= 4; a[[\"x\"]] -= 4; a", "<|\"x\" -> 497/2, \"y\" -> 23|>");
    CHECK("--a[[\"y\"]]", "22");
}

void test_compound_missing_key() {
    /* Mathematica does not reject a missing key: the value is Missing[...]
     * and the stored result is the unevaluated sum. */
    run("a = <|\"x\" -> 1|>");
    CHECK("a[\"new\"] += 1", "1 + Missing[\"KeyAbsent\", \"new\"]");
    CHECK("a[\"n2\"]++", "Missing[\"KeyAbsent\", \"n2\"]");
    CHECK("Keys[a]", "{\"x\", \"new\", \"n2\"}");
}

void test_counting_idiom() {
    CHECK("c = <||>; Scan[If[KeyExistsQ[c, #], c[#] += 1, c[#] = 1] &, {1, 2, 1, 3, 1, 2}]; c",
          "<|1 -> 3, 2 -> 2, 3 -> 1|>");
    CHECK("a = <|\"c\" -> 0|>; Do[a[\"c\"] += i, {i, 10}]; a", "<|\"c\" -> 55|>");
}

void test_appendto_element() {
    CHECK("b = <|\"l\" -> {1}|>; AppendTo[b[\"l\"], 2]", "{1, 2}");
    CHECK("b", "<|\"l\" -> {1, 2}|>");
    CHECK("AppendTo[b[[\"l\"]], 5]; PrependTo[b[\"l\"], 0]; b", "<|\"l\" -> {0, 1, 2, 5}|>");
}

/* ---------- deletion ---------- */

void test_unset_element() {
    run("a = <|\"x\" -> 1, \"y\" -> 2|>");
    CHECK("a[\"x\"] =.", "Null");
    CHECK("a", "<|\"y\" -> 2|>");
    CHECK("a[\"nokey\"] =.; a", "<|\"y\" -> 2|>");
    CHECK("n = <|\"p\" -> <|\"q\" -> 1, \"r\" -> 2|>|>; n[\"p\"][\"q\"] =.; n", "<|\"p\" -> <|\"r\" -> 2|>|>");
    CHECK("n = <|\"p\" -> <|\"q\" -> 1, \"r\" -> 2|>|>; n[\"p\", \"q\"] =.; n", "<|\"p\" -> <|\"r\" -> 2|>|>");
}

void test_keydropfrom() {
    CHECK("c = <|1 -> 2, 3 -> 4|>; KeyDropFrom[c, 1]", "<|3 -> 4|>");
    CHECK("c", "<|3 -> 4|>");
    CHECK("c = <|1 -> 2, 3 -> 4, 5 -> 6|>; KeyDropFrom[c, {1, 5, 7}]", "<|3 -> 4|>");
    CHECK("KeyDropFrom[c, 99]", "<|3 -> 4|>");
    CHECK("KeyDropFrom[c, {}]", "<|3 -> 4|>");
    CHECK("a = <|\"x\" -> 1|>; KeyDropFrom[a, Key[\"x\"]]", "<||>");
    CHECK("AssociationQ[a]", "True");
    CHECK("Attributes[KeyDropFrom]", "{HoldFirst, Protected}");
    CHECK("KeyDropFrom[{1}, 1]", "KeyDropFrom[{1}, 1]");
    CHECK("KeyDropFrom[unsetsym, 1]", "KeyDropFrom[unsetsym, 1]");
    CHECK_FF("a = <|\"x\" :> 1 + 1, \"y\" -> 2|>; KeyDropFrom[a, \"y\"]",
             "Association[RuleDelayed[\"x\", Plus[1, 1]]]");
    /* The many-key path (one filtering pass). */
    CHECK("a = AssociationThread[Range[20], Range[20]]; KeyDropFrom[a, Range[2, 20]]", "<|1 -> 1|>");
    CHECK("a = AssociationThread[Range[20], Range[20]]; KeyDropFrom[a, Key /@ Range[3, 20]]", "<|1 -> 1, 2 -> 2|>");
    /* Deleting then re-adding keeps the index consistent. */
    CHECK("a = AssociationThread[Range[10], Range[10]]; Do[KeyDropFrom[a, i], {i, 1, 10, 2}]; "
          "a[3] = 33; a[11] = 11; {a, a[4], a[3], a[11], KeyExistsQ[a, 5]}",
          "{<|2 -> 2, 4 -> 4, 6 -> 6, 8 -> 8, 10 -> 10, 3 -> 33, 11 -> 11|>, 4, 33, 11, False}");
}

void test_delete_key() {
    run("d = <|\"a\" -> 1, \"b\" -> 2|>");
    CHECK("Delete[d, Key[\"a\"]]", "<|\"b\" -> 2|>");
    CHECK("Delete[d, {Key[\"a\"]}]", "<|\"b\" -> 2|>");
    CHECK("Delete[d, \"a\"]", "<|\"b\" -> 2|>");
    CHECK("Delete[d, {{Key[\"a\"]}, {Key[\"b\"]}}]", "<||>");
    CHECK("Delete[d, Key[\"zz\"]]", "<|\"a\" -> 1, \"b\" -> 2|>");
    CHECK("Delete[d, -1]", "<|\"a\" -> 1|>");
    CHECK("Delete[<|\"a\" -> <|\"b\" -> 1, \"c\" -> 2|>|>, {\"a\", \"b\"}]", "<|\"a\" -> <|\"c\" -> 2|>|>");
    CHECK("d", "<|\"a\" -> 1, \"b\" -> 2|>");
}

/* ---------- nested and curried assignment ---------- */

void test_nested_assignment() {
    run("n = <|\"p\" -> <|\"q\" -> 1|>|>");
    CHECK("n[\"p\"][\"q\"] = 5; n", "<|\"p\" -> <|\"q\" -> 5|>|>");
    CHECK("n[\"p\"][\"r\"] = 6; n", "<|\"p\" -> <|\"q\" -> 5, \"r\" -> 6|>|>");
    CHECK("n[[\"p\", \"s\"]] = 7; n", "<|\"p\" -> <|\"q\" -> 5, \"r\" -> 6, \"s\" -> 7|>|>");
    CHECK("n[\"p\", \"t\"] = 8; Keys[n[\"p\"]]", "{\"q\", \"r\", \"s\", \"t\"}");
    CHECK("n[\"p\", \"q\"] += 1; n[\"p\", \"q\"]", "6");
    CHECK("n[[\"p\", \"q\"]]++; n[\"p\"][\"q\"]", "7");
    CHECK("n[\"p\"][\"q\"] *= 2; n[\"p\", \"q\"]", "14");
    /* A missing intermediate key: Set::kval, no change, Set yields the rhs. */
    CHECK("n[\"zz\", \"q\"] = 1", "1");
    CHECK("KeyExistsQ[n, \"zz\"]", "False");
    /* The Part form stores under the missing key instead (Mathematica 15). */
    CHECK("n[[\"zz\", \"q\"]] = 1", "1");
    CHECK("{KeyExistsQ[n, \"zz\"], n[\"zz\"]}", "{True, 1}");
    CHECK("KeyDropFrom[n, \"zz\"]; Keys[n]", "{\"p\"}");
    /* Through a non-association value. */
    CHECK("m = <|\"p\" -> 5|>; m[\"p\"][\"k\"] = 1; m", "<|\"p\" -> 5|>");
    CHECK("m = <|\"x\" -> {1, 2}|>; m[[\"x\", 1]] = 7; m", "<|\"x\" -> {7, 2}|>");
    CHECK("AssociateTo[n[\"p\"], \"u\" -> 9]; n[\"p\", \"u\"]", "9");
}

void test_part_assignment() {
    run("a = <|\"x\" -> 1|>");
    CHECK("a[[Key[\"new\"]]] = 5; a", "<|\"x\" -> 1, \"new\" -> 5|>");
    CHECK("a[[\"new2\"]] = 6; a", "<|\"x\" -> 1, \"new\" -> 5, \"new2\" -> 6|>");
    CHECK("a[[9]] = 1", "1");                         /* Set::partw, unchanged */
    CHECK("Length[a]", "3");
    CHECK("a[[-1]] = 60; a[[\"new2\"]]", "60");
    CHECK("a[[All]] = 0; a", "<|\"x\" -> 0, \"new\" -> 0, \"new2\" -> 0|>");
    CHECK("b = <|1 -> \"one\"|>; b[1] = \"uno\"; b[[1]] = \"eins\"; b", "<|1 -> \"eins\"|>");
    /* Part assignment works on the STORED value, as in Mathematica. */
    CHECK("x = <|\"a\" -> yy|>; yy = 5; x[\"b\"] = 1; OwnValues[x]",
          "{HoldPattern[x] :> <|\"a\" -> yy, \"b\" -> 1|>}");
}

/* ---------- RuleDelayed ---------- */

void test_ruledelayed() {
    CHECK_FF("<|\"a\" -> 1, \"a\" :> 2|>", "Association[RuleDelayed[\"a\", 2]]");
    CHECK_FF("Association[{p :> 1 + 1}]", "Association[RuleDelayed[p, Plus[1, 1]]]");
    CHECK_FF("<|p :> 1 + 1, q -> 2, p :> 2 + 2|>",
             "Association[RuleDelayed[p, Plus[2, 2]], Rule[q, 2]]");
    /* Set replaces a delayed entry by an immediate one; := makes one delayed. */
    CHECK_FF("a = <|\"x\" :> 1 + 1|>; a[\"x\"] = 5; a", "Association[Rule[\"x\", 5]]");
    CHECK_FF("a = <|\"x\" :> 1 + 1|>; a[\"x\"] += 5; a", "Association[Rule[\"x\", 7]]");
    CHECK_FF("a = <|\"x\" -> 1|>; a[\"q\"] := 7; a",
             "Association[Rule[\"x\", 1], RuleDelayed[\"q\", 7]]");
    CHECK_FF("a = <|\"x\" -> 1|>; a[\"x\"] := 1 + 1; a", "Association[RuleDelayed[\"x\", Plus[1, 1]]]");
    CHECK_FF("a = <||>; AssociateTo[a, \"k\" :> 1 + 1]", "Association[RuleDelayed[\"k\", Plus[1, 1]]]");
}

/* ---------- AssociateTo ---------- */

void test_associateto() {
    CHECK("a = <|1 -> 2|>; AssociateTo[a, <|3 -> 4, 1 -> 9|>]", "<|1 -> 9, 3 -> 4|>");
    CHECK("a = <|1 -> 2|>; AssociateTo[a, {5 -> 6, 7 -> 8}]", "<|1 -> 2, 5 -> 6, 7 -> 8|>");
    CHECK("a", "<|1 -> 2, 5 -> 6, 7 -> 8|>");
    CHECK("q = 5; AssociateTo[q, 1 -> 2]", "AssociateTo[q, 1 -> 2]");
    CHECK("AssociateTo[unsetsym, 1 -> 2]", "AssociateTo[unsetsym, 1 -> 2]");
    CHECK("Module[{b = <|\"k\" -> 1|>}, b[\"k\"] += 1; b[\"j\"] = 5; b]", "<|\"k\" -> 2, \"j\" -> 5|>");
}

/* ---------- aliasing: a write is never visible through another holder ---------- */

void test_alias_every_form() {
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"] = 2; {a, b}", "{<|\"x\" -> 2|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"y\"] = 2; {a, b}", "{<|\"x\" -> 1, \"y\" -> 2|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; AssociateTo[a, \"y\" -> 3]; {a, b}",
          "{<|\"x\" -> 1, \"y\" -> 3|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"]++; {a, b}", "{<|\"x\" -> 2|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"] += 5; {a, b}", "{<|\"x\" -> 6|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"] *= 5; {a, b}", "{<|\"x\" -> 5|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; KeyDropFrom[a, \"x\"]; {a, b}", "{<||>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"] =.; {a, b}", "{<||>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[[\"x\"]] = 9; {a, b}", "{<|\"x\" -> 9|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[[1]] = 9; {a, b}", "{<|\"x\" -> 9|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; b = a; a[\"x\"] := 3; {a, b}", "{<|\"x\" :> 3|>, <|\"x\" -> 1|>}");
    CHECK("a = <|\"x\" -> {1}|>; b = a; AppendTo[a[\"x\"], 2]; {a, b}",
          "{<|\"x\" -> {1, 2}|>, <|\"x\" -> {1}|>}");
}

void test_alias_nested() {
    CHECK("n = <|\"p\" -> <|\"q\" -> 1|>|>; m = n; n[\"p\", \"q\"] = 9; {n, m}",
          "{<|\"p\" -> <|\"q\" -> 9|>|>, <|\"p\" -> <|\"q\" -> 1|>|>}");
    /* The inner association is shared with another variable: still untouched. */
    CHECK("n = <|\"p\" -> <|\"q\" -> 1|>|>; inner = n[\"p\"]; n[\"p\", \"q\"] = 9; {n, inner}",
          "{<|\"p\" -> <|\"q\" -> 9|>|>, <|\"q\" -> 1|>}");
    CHECK("n = <|\"p\" -> <|\"q\" -> 1|>|>; inner = n[\"p\"]; n[\"p\"][\"z\"] = 2; {n, inner}",
          "{<|\"p\" -> <|\"q\" -> 1, \"z\" -> 2|>|>, <|\"q\" -> 1|>}");
    CHECK("n = <|\"p\" -> <|\"q\" -> 1|>|>; inner = n[\"p\"]; n[\"p\", \"q\"] =.; {n, inner}",
          "{<|\"p\" -> <||>|>, <|\"q\" -> 1|>}");
    /* Two associations sharing an ENTRY (built from the same rule). */
    CHECK("r = \"k\" -> 1; a = <|r|>; b = <|r, \"j\" -> 2|>; a[\"k\"] = 5; {a, b, r}",
          "{<|\"k\" -> 5|>, <|\"k\" -> 1, \"j\" -> 2|>, \"k\" -> 1}");
}

void test_alias_loops_and_module() {
    CHECK("a = <|1 -> 1|>; b = a; Do[a[i] = i, {i, 2, 4}]; {a, b}",
          "{<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, <|1 -> 1|>}");
    CHECK("a = <|1 -> 1|>; b = a; Do[AssociateTo[a, i -> i], {i, 2, 4}]; {a, b}",
          "{<|1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4|>, <|1 -> 1|>}");
    CHECK("a = <|\"x\" -> 1|>; Module[{}, b = a; Do[a[\"x\"] += 1, {3}]]; {a, b}",
          "{<|\"x\" -> 4|>, <|\"x\" -> 1|>}");
    /* A snapshot taken inside the loop keeps its own state. */
    CHECK("a = <||>; snaps = Table[a[i] = i; a, {i, 3}]",
          "{<|1 -> 1|>, <|1 -> 1, 2 -> 2|>, <|1 -> 1, 2 -> 2, 3 -> 3|>}");
    CHECK("Module[{u = <|\"k\" -> 0|>, v}, v = u; Do[u[\"k\"]++; u[i] = i, {i, 3}]; {u, v}]",
          "{<|\"k\" -> 3, 1 -> 1, 2 -> 2, 3 -> 3|>, <|\"k\" -> 0|>}");
    CHECK("Module[{u = <|1 -> 1, 2 -> 2|>, v}, v = u; Do[KeyDropFrom[u, i], {i, 2}]; {u, v}]",
          "{<||>, <|1 -> 1, 2 -> 2|>}");
    /* A list holding the association, then written through. */
    CHECK("a = <|\"x\" -> 1|>; l = {a, a}; a[\"x\"] = 2; {a, l}",
          "{<|\"x\" -> 2|>, {<|\"x\" -> 1|>, <|\"x\" -> 1|>}}");
    /* Out-of-loop value from a function returning the association. */
    CHECK("a = <|\"x\" -> 1|>; f[] := a; g = f[]; a[\"x\"] = 3; {a, g}",
          "{<|\"x\" -> 3|>, <|\"x\" -> 1|>}");
}

/* ---------- fixed-point bookkeeping: never serve a stale value ---------- */

void test_no_stale_values() {
    /* A symbol value that changes later is still re-evaluated on read. */
    CHECK("Clear[w]; a = <||>; a[\"k\"] = w; w = 7; a[\"k\"]", "7");
    CHECK("a = <|\"k\" -> 1|>; x1 = a; a[\"k\"] = 2; {x1[\"k\"], a[\"k\"], Lookup[a, \"k\"]}", "{1, 2, 2}");
    CHECK("a = <|1 -> True|>; Do[a[i] = True, {i, 2, 5}]; {Length[a], a[5], KeyExistsQ[a, 6]}",
          "{5, True, False}");
    CHECK("a = <||>; Do[a[i] = i^2; If[a[i] != i^2, Print[\"stale\"]], {i, 50}]; Total[a]", "42925");
    CHECK("a = <|\"k\" -> 1|>; Hash[a] == Hash[<|\"k\" -> 1|>]", "True");
    CHECK("a[\"k\"] = 2; Hash[a] == Hash[<|\"k\" -> 2|>]", "True");
    CHECK("a === <|\"k\" -> 2|>", "True");
}

/* ---------- index / scan consistency on malformed nodes ---------- */

void test_malformed_consistency() {
    CHECK("Lookup[Association[f[1, 2]], 1]", "Lookup[Association[f[1, 2]], 1]");
    CHECK("KeyExistsQ[Association[f[1, 2]], 1]", "False");
}

/* ---------- C-level invariants ---------- */

void test_inplace_invariants() {
    run("inv = <||>; Do[inv[i] = i, {i, 1000}]");
    Expr* a = stored("inv");
    CHECK_TRUE(a != NULL);
    CHECK_TRUE(a->refcount == 1);                       /* updated in place */
    CHECK_TRUE(a->data.function.arg_count == 1000);
    CHECK_TRUE(a->data.function.index != NULL);         /* index maintained */
    CHECK_TRUE(assoc_index_matches(a->data.function.index, a->data.function.args, 1000));
    CHECK_TRUE(eval_node_is_ground(a));                 /* still a GROUND fixed point */
    Expr* before = a;
    run("Do[inv[i] = 2 i, {i, 1000}]; Do[inv[i] =., {i, 1, 1000, 2}]");
    a = stored("inv");
    CHECK_TRUE(a == before);                            /* same node throughout */
    CHECK_TRUE(a->data.function.arg_count == 500);
    CHECK_TRUE(assoc_index_matches(a->data.function.index, a->data.function.args, 500));
    CHECK("{inv[2], inv[1000], KeyExistsQ[inv, 999], Total[inv]}", "{4, 2000, False, 501000}");
    /* An alias forces exactly one copy; the alias keeps the old node. */
    run("inv2 = inv; inv[2] = -1");
    CHECK_TRUE(stored("inv") != stored("inv2"));
    CHECK_TRUE(stored("inv2") == before);
    CHECK("{inv[2], inv2[2]}", "{-1, 4}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_compound_ops_keyed);
    TEST(test_compound_ops_part_forms);
    TEST(test_compound_missing_key);
    TEST(test_counting_idiom);
    TEST(test_appendto_element);
    TEST(test_unset_element);
    TEST(test_keydropfrom);
    TEST(test_delete_key);
    TEST(test_nested_assignment);
    TEST(test_part_assignment);
    TEST(test_ruledelayed);
    TEST(test_associateto);
    TEST(test_alias_every_form);
    TEST(test_alias_nested);
    TEST(test_alias_loops_and_module);
    TEST(test_no_stale_values);
    TEST(test_malformed_consistency);
    TEST(test_inplace_invariants);

    printf("All Association mutation tests passed (%d assertions).\n", g_checks);
    return 0;
}
