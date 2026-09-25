/* Tests for src/assoc_ops.c: the second-tier Association heads
 * (KeyIntersection, KeyComplement, MissingQ, Discard, CountDistinct[By],
 * SubsetQ, JoinAcross, ApplyTo, AssociationComap, Splice), the extended forms
 * of Keys/Values/KeyUnion/AssociationMap/PositionIndex/Merge/Association/
 * Normal/DeleteMissing/Lookup/Transpose, and Listable threading over
 * association values.
 *
 * Every expected string was produced by Mathematica 15 for the same input
 * (ToString[expr, InputForm]) and matches Mathilda byte for byte, except where
 * a comment says otherwise. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "expr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_checks = 0;

/* Evaluate ToString[(input), InputForm] and compare the resulting string. */
static void CHECK(const char* input, const char* expected) {
    size_t n = strlen(input) + 32;
    char* wrapped = malloc(n);
    snprintf(wrapped, n, "ToString[(%s), InputForm]", input);
    struct Expr* parsed = parse_expression(wrapped);
    assert(parsed != NULL);
    struct Expr* r = evaluate(parsed);
    expr_free(parsed);
    const char* got = (r && r->type == EXPR_STRING) ? r->data.string : "<not a string>";
    if (strcmp(got, expected) != 0)
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, got);
    assert(strcmp(got, expected) == 0);
    expr_free(r);
    free(wrapped);
    g_checks++;
}

static void test_threading(void) {
    CHECK("<|a->1,b->2|> + 1",
          "<|a -> 2, b -> 3|>");
    CHECK("2 <|a->1,b->2|>",
          "<|a -> 2, b -> 4|>");
    CHECK("<|a->1,b->2|> + <|a->10,b->20|>",
          "<|a -> 11, b -> 22|>");
    CHECK("<|a->1,b->2|> + <|b->10,a->20|>",
          "<|a -> 1, b -> 2|> + <|b -> 10, a -> 20|>");
    CHECK("<|a->1|> + <|a->1,b->2|>",
          "<|a -> 1|> + <|a -> 1, b -> 2|>");
    CHECK("Sqrt[<|a->4,b->9|>]",
          "<|a -> 2, b -> 3|>");
    CHECK("Sin[<|a->0|>]",
          "<|a -> 0|>");
    CHECK("<|a->2|>^3",
          "<|a -> 8|>");
    CHECK("2^<|a->2,b->3|>",
          "<|a -> 4, b -> 8|>");
    CHECK("<|a->1,b->2|> + {1,2}",
          "{<|a -> 2, b -> 3|>, <|a -> 3, b -> 4|>}");
    CHECK("{<|a->1|>,<|a->2|>} + 1",
          "{<|a -> 2|>, <|a -> 3|>}");
    CHECK("<|a->{1,2}|> + 1",
          "<|a -> {2, 3}|>");
    CHECK("<|a->1,b->2|> + x",
          "<|a -> 1 + x, b -> 2 + x|>");
    CHECK("<|a->1,b->2|> - <|a->1,b->1|>",
          "<|a -> 0, b -> 1|>");
    CHECK("<|a->1,b->2|> / 2",
          "<|a -> 1/2, b -> 1|>");
    CHECK("<|a->1,b->2|>*<|a->3,b->4|>",
          "<|a -> 3, b -> 8|>");
    CHECK("<||> + 1",
          "<||>");
    CHECK("-<|a->1|>",
          "<|a -> -1|>");
    CHECK("1/<|a->2|>",
          "<|a -> 1/2|>");
    CHECK("N[<|a->1/2|>]",
          "<|a -> 0.5|>");
    CHECK("Abs[<|a->-2|>]",
          "<|a -> 2|>");
    CHECK("Total[{<|a->1,b->2|>,<|a->3,b->4|>}]",
          "<|a -> 4, b -> 6|>");
    CHECK("Mean[{<|a->1,b->2|>,<|a->3,b->4|>}]",
          "<|a -> 2, b -> 3|>");
    CHECK("Total[<|a->1,b->2|>]",
          "3");
    CHECK("Mean[<|a->1,b->2|>]",
          "3/2");
    CHECK("Mean[{<|a->1|>,<|a->2|>, <|a->4|>}]",
          "<|a -> 7/3|>");
    CHECK("Total[<|a->{1,2},b->{3,4}|>]",
          "{4, 6}");
    CHECK("lf[<|a->1,b->2|>]",
          "lf[<|a -> 1, b -> 2|>]");
    CHECK("SetAttributes[lf2, Listable]; lf2[<|a->1|>, {1,2}]",
          "{<|a -> lf2[1, 1]|>, <|a -> lf2[1, 2]|>}");
    CHECK("StringLength[<|a->\"xx\"|>]",
          "<|a -> 2|>");
    CHECK("Floor[<|a->1.5|>]",
          "<|a -> 1|>");
    CHECK("Mod[<|a->5|>, 3]",
          "<|a -> 2|>");
    /* A RuleDelayed entry stays delayed. (Mathematica sorts the held Plus as
     * Plus[1, Plus[1, 1]]; the argument order inside the held value differs.) */
    CHECK("FullForm[<|a:>1+1, b->2|> + 1]",
          "Association[RuleDelayed[a, Plus[Plus[1, 1], 1]], Rule[b, 3]]");
    CHECK("Options[JoinAcross]",
          "{KeyCollisionFunction -> Left}");
}

static void test_listable_regression(void) {
    CHECK("{1,2} + {3,4}",
          "{4, 6}");
    CHECK("Sin[{0, x}]",
          "{0, Sin[x]}");
    CHECK("f[{1,2}]",
          "f[{1, 2}]");
    /* Calls with neither a List nor an Association argument are untouched. */
    CHECK("Plus[a, b, 2]", "2 + a + b");
    CHECK("Sin[Pi/6]", "1/2");
}

/* Packed and visible-NDArray arguments: association threading carries a
 * buffer value through, Splice accepts a packed list, and the structural heads
 * answer the same for a visible NDArray as for the List it represents. (No
 * Mathematica counterpart for NDArray; the expected values follow the List.) */
static void test_packed_and_ndarray(void) {
    CHECK("<|a->Range[3]|> + 1", "<|a -> {2, 3, 4}|>");
    CHECK("Range[2] + <|a->1|>", "{<|a -> 2|>, <|a -> 3|>}");
    CHECK("Sin[<|a -> NDArray[{0.5}]|>]", "<|a -> NDArray[{0.479426}]|>");
    CHECK("{0, Splice[Range[3]]}", "{0, 1, 2, 3}");
    CHECK("CountDistinct[Range[10]]", "10");
    CHECK("CountDistinct[NDArray[{1.5, 2.5, 1.5}]]", "2");
    CHECK("Discard[NDArray[{1,2,3,4}, DataType -> \"int64\"], EvenQ]", "{1, 3}");
    CHECK("SubsetQ[NDArray[{1,2,3}, DataType -> \"int64\"], {1}]", "True");
    CHECK("CountDistinctBy[NDArray[{1,2,3}, DataType -> \"int64\"], EvenQ]", "2");
}

static void test_keyintersection(void) {
    CHECK("KeyIntersection[{<|a->1,b->2|>,<|b->3,c->4|>}]",
          "{<|b -> 2|>, <|b -> 3|>}");
    CHECK("KeyIntersection[{<|a->1,b->2|>,<|b->3,a->4|>}]",
          "{<|a -> 1, b -> 2|>, <|a -> 4, b -> 3|>}");
    CHECK("KeyIntersection[{}]",
          "{}");
    CHECK("KeyIntersection[{<|a->1,b->2|>, {b->5}}]",
          "{<|b -> 2|>, <|b -> 5|>}");
    CHECK("KeyIntersection[<|a->1|>]",
          "KeyIntersection[<|a -> 1|>]");
}

static void test_keycomplement(void) {
    CHECK("KeyComplement[{<|a->1,b->2|>,<|b->3|>}]",
          "<|a -> 1|>");
    CHECK("KeyComplement[{<|a->1,b->2,c->3|>,<|b->3|>,<|c->1|>}]",
          "<|a -> 1|>");
    CHECK("KeyComplement[{a->1,b->2}]",
          "<|a -> 1|>");
    CHECK("KeyComplement[{}]",
          "KeyComplement[{}]");
}

static void test_missingq(void) {
    CHECK("MissingQ[Missing[]]",
          "True");
    CHECK("MissingQ[Missing[\"KeyAbsent\",x]]",
          "True");
    CHECK("MissingQ[3]",
          "False");
    CHECK("MissingQ[{Missing[]}]",
          "False");
    CHECK("MissingQ[Missing]",
          "False");
}

static void test_discard(void) {
    CHECK("Discard[{1,2,3,4}, EvenQ]",
          "{1, 3}");
    CHECK("Discard[<|a->1,b->2|>, EvenQ]",
          "<|a -> 1|>");
    CHECK("Discard[{1,2,3,4,5,6},EvenQ,2]",
          "{1, 3, 5, 6}");
    CHECK("Discard[{2,4,3}, EvenQ, Infinity]",
          "{3}");
    CHECK("Discard[f[1,2,3], OddQ]",
          "f[2]");
    CHECK("Discard[{1,2,3}, EvenQ, -1]",
          "Discard[{1, 2, 3}, EvenQ, -1]");
    CHECK("Discard[x, EvenQ]",
          "Discard[x, EvenQ]");
}

static void test_countdistinct(void) {
    CHECK("CountDistinct[{1,2,1,3}]",
          "3");
    CHECK("CountDistinct[<|a->1,b->1,c->2|>]",
          "2");
    CHECK("CountDistinct[{}]",
          "0");
    CHECK("CountDistinct[{1,1.,2}]",
          "3");
    CHECK("CountDistinct[f[1,1,2]]",
          "2");
    CHECK("CountDistinct[x]",
          "CountDistinct[x]");
    CHECK("CountDistinctBy[{1,2,3,4},EvenQ]",
          "2");
    CHECK("CountDistinctBy[<|a->1,b->2,c->3|>,OddQ]",
          "2");
}

static void test_subsetq(void) {
    CHECK("SubsetQ[{1,2,3},{1,2}]",
          "True");
    CHECK("SubsetQ[{1,2},{1,3}]",
          "False");
    CHECK("SubsetQ[<|a->1,b->2|>,<|a->1|>]",
          "True");
    CHECK("SubsetQ[<|a->1,b->2|>,<|a->2|>]",
          "True");
    CHECK("SubsetQ[<|a->1,b->2|>,{1}]",
          "True");
    CHECK("SubsetQ[{1,2},{}]",
          "True");
    CHECK("SubsetQ[{1,1,2},{1,1,1}]",
          "True");
    CHECK("SubsetQ[f[1,2],f[1]]",
          "True");
    CHECK("SubsetQ[f[1],g[1]]",
          "SubsetQ[f[1], g[1]]");
}

static void test_splice(void) {
    CHECK("{1, Splice[{2,3}], 4}",
          "{1, 2, 3, 4}");
    CHECK("<|a->1, Splice[{b->2,c->3}]|>",
          "<|a -> 1, b -> 2, c -> 3|>");
    CHECK("f[1, Splice[{2,3}]]",
          "f[1, Splice[{2, 3}]]");
    CHECK("g[Splice[{1,2},_]]",
          "g[1, 2]");
    CHECK("h[Splice[{1,2},h]]",
          "h[1, 2]");
    CHECK("k[Splice[{1,2},Except[h]]]",
          "k[1, 2]");
    CHECK("h[Splice[{1,2},Except[h]]]",
          "h[Splice[{1, 2}, Except[h]]]");
    CHECK("{Splice[{}]}",
          "{}");
    CHECK("Table[Splice[{i,i}], {i,2}]",
          "{1, 1, 2, 2}");
    CHECK("Map[Splice[{#,#}]&, {1,2}]",
          "{1, 1, 2, 2}");
    CHECK("Attributes[Splice]",
          "{Protected}");
}

static void test_keys_values(void) {
    CHECK("Keys[<|a->1,b->2|>, f]",
          "{f[a], f[b]}");
    CHECK("Values[<|a->1,b->2|>, f]",
          "{f[1], f[2]}");
    CHECK("Keys[{<|a->1|>,<|b->2|>}]",
          "{{a}, {b}}");
    CHECK("Values[{<|a->1|>,<|b->2|>}]",
          "{{1}, {2}}");
    CHECK("Keys[{a->1,b->2}, f]",
          "{f[a], f[b]}");
    CHECK("Keys[{{a->1},{b->2}}]",
          "{{a}, {b}}");
    CHECK("Keys[{a->1, <|b->2|>}]",
          "{a, {b}}");
    CHECK("Values[{<|a->1|>, {<|b->2|>}}]",
          "{{1}, {{2}}}");
    CHECK("Keys[{<|a->1|>,<|b->2|>}, f]",
          "{{f[a]}, {f[b]}}");
    CHECK("Keys[{}, f]",
          "{}");
    CHECK("Keys[{<|a->1|>, 3}]",
          "Keys[{<|a -> 1|>, 3}]");
}

static void test_keyunion(void) {
    CHECK("KeyUnion[{<|a->1|>,<|b->2|>}, 0&]",
          "{<|a -> 1, b -> 0|>, <|a -> 0, b -> 2|>}");
    CHECK("KeyUnion[{<|a->1|>,<|b->2|>}, f]",
          "{<|a -> 1, b -> f[b]|>, <|a -> f[a], b -> 2|>}");
    CHECK("KeyUnion[{{a->1},{b->2}}]",
          "{<|a -> 1, b -> Missing[\"KeyAbsent\", b]|>, <|a -> Missing[\"KeyAbsent\", a], b -> 2|>}");
    CHECK("KeyUnion[{<|a->1|>}, f]",
          "{<|a -> 1|>}");
}

static void test_associationmap(void) {
    CHECK("AssociationMap[Reverse, <|a->1,b->2|>]",
          "<|1 -> a, 2 -> b|>");
    CHECK("AssociationMap[#[[1]]->#[[2]]+1&, <|a->1,b->2|>]",
          "<|a -> 2, b -> 3|>");
    CHECK("AssociationMap[{#,#}&, <|a->1|>]",
          "<|a -> 1|>");
    CHECK("AssociationMap[Nothing&, <|a->1|>]",
          "<||>");
    CHECK("AssociationMap[<|b->1|>&, <|a->1, c->2|>]",
          "<|b -> 1|>");
    CHECK("AssociationMap[f, {a,b}]",
          "<|a -> f[a], b -> f[b]|>");
}

static void test_positionindex(void) {
    CHECK("PositionIndex[<|a->x,b->y,c->x|>]",
          "<|x -> {a, c}, y -> {b}|>");
    CHECK("PositionIndex[<|a->1,b->1|>]",
          "<|1 -> {a, b}|>");
    CHECK("PositionIndex[{a,b,a}]",
          "<|a -> {1, 3}, b -> {2}|>");
}

static void test_merge(void) {
    CHECK("Merge[{a->1, a->2}, Total]",
          "<|a -> 3|>");
    CHECK("Merge[{a->1, b->2, a->3}, f]",
          "<|a -> f[{1, 3}], b -> f[{2}]|>");
    CHECK("Merge[{{a->1}, {a->2}}, Total]",
          "<|a -> 3|>");
    CHECK("Merge[{<|a->1|>, a->2}, Total]",
          "<|a -> 3|>");
    CHECK("Merge[{{a->1, a->2}}, f]",
          "<|a -> f[{1, 2}]|>");
    CHECK("Merge[{3}, f]",
          "Merge[{3}, f]");
}

static void test_association(void) {
    CHECK("Association[{{a->1},{b->2}}]",
          "<|a -> 1, b -> 2|>");
    CHECK("Association[{a->1,{b->2,{c->3}}}]",
          "<|a -> 1, b -> 2, c -> 3|>");
    CHECK("Association[{<|a->1|>, {b->2}}]",
          "<|a -> 1, b -> 2|>");
    CHECK("Association[{{}}]",
          "<||>");
    CHECK("Association[{{a->1}, {a->2}}]",
          "<|a -> 2|>");
}

static void test_normal(void) {
    CHECK("Normal[{<|a->1|>, f[<|b->2|>]}]",
          "{{a -> 1}, f[{b -> 2}]}");
    CHECK("Normal[f[g[<|b->2|>]]]",
          "f[g[{b -> 2}]]");
    CHECK("Normal[<|a->{<|b->1|>}|>]",
          "{a -> {<|b -> 1|>}}");
    CHECK("Normal[f[<|a-><|b->1|>|>]]",
          "f[{a -> <|b -> 1|>}]");
    CHECK("Normal[{<|a->1|>, f[<|b->2|>]}, Association]",
          "{{a -> 1}, f[{b -> 2}]}");
    CHECK("Normal[f[SeriesData[x,0,{1,1},0,2,1]]]",
          "f[1 + x]");
}

static void test_counts_tally(void) {
    CHECK("Tally[<|x->a,y->b,z->a|>]",
          "{{a, 2}, {b, 1}}");
    CHECK("Counts[<|x->a,y->b,z->a|>]",
          "<|a -> 2, b -> 1|>");
    CHECK("Tally[Counts[{a,b,a}]]",
          "{{2, 1}, {1, 1}}");
    CHECK("Counts[Tally[{a,b,a}]]",
          "<|{a, 2} -> 1, {b, 1} -> 1|>");
}

static void test_deletemissing(void) {
    CHECK("DeleteMissing[{1, {Missing[], 2}}, 2]",
          "{1, {2}}");
    CHECK("DeleteMissing[{<|a->1,b->Missing[]|>, Missing[]}, 2]",
          "{<|a -> 1|>}");
    CHECK("DeleteMissing[{{Missing[]}, {1, Missing[]}}, 1, 1]",
          "{}");
    CHECK("DeleteMissing[{{{Missing[]}}, 1}, 3, 2]",
          "{{{}}, 1}");
    CHECK("DeleteMissing[{Missing[], {Missing[]}}, Infinity]",
          "{{}}");
    CHECK("DeleteMissing[<|a-><|b->Missing[]|>|>, 2]",
          "<|a -> <||>|>");
    CHECK("DeleteMissing[{{1, {Missing[]}}, {2}}, 1, 2]",
          "{{2}}");
    CHECK("DeleteMissing[{1, Missing[], 2}]",
          "{1, 2}");
    CHECK("DeleteMissing[3]",
          "DeleteMissing[3]");
}

static void test_lookup(void) {
    CHECK("Lookup[{<|a->1|>,{a->3},{}}, a, 0]",
          "{1, 3, 0}");
    CHECK("Lookup[{{a->1}, {b->2}}, a]",
          "{1, Missing[\"KeyAbsent\", a]}");
    CHECK("Lookup[{<|a->1|>,<|b->2|>}, {a,b}, 0]",
          "{{1, 0}, {0, 2}}");
    CHECK("Lookup[{a->1, b->2}, b]",
          "2");
}

static void test_joinacross(void) {
    CHECK("JoinAcross[{<|a->1,b->x|>,<|a->2,b->y|>}, {<|a->1,c->p|>,<|a->3,c->q|>}, Key[a]]",
          "{<|a -> 1, b -> x, c -> p|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>,<|a->2,b->y|>}, {<|a->1,c->p|>,<|a->3,c->q|>}, a, \"Left\"]",
          "{<|a -> 1, b -> x, c -> p|>, <|a -> 2, b -> y, c -> Missing[\"Unmatched\"]|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>,<|a->2,b->y|>}, {<|a->1,c->p|>,<|a->3,c->q|>}, a, \"Right\"]",
          "{<|a -> 1, b -> x, c -> p|>, <|a -> 3, b -> Missing[\"Unmatched\"], c -> q|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>,<|a->2,b->y|>}, {<|a->1,c->p|>,<|a->3,c->q|>}, a, \"Outer\"]",
          "{<|a -> 1, b -> x, c -> p|>, <|a -> 2, b -> y, c -> Missing[\"Unmatched\"]|>, <|a -> 3, b -> Missing[\"Unmatched\"], c -> q|>}");
    CHECK("JoinAcross[{<|a->1,b->x,d->1|>}, {<|a->1,c->p,d->1|>,<|a->1,c->r,d->2|>}, {a,d}]",
          "{<|a -> 1, b -> x, d -> 1, c -> p|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>}, {<|e->1,c->p|>}, Key[a]->Key[e]]",
          "{<|a -> 1, b -> x, e -> 1, c -> p|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>, <|a->1,b->z|>}, {<|a->1,c->p|>,<|a->1,c->r|>}, a]",
          "{<|a -> 1, b -> x, c -> p|>, <|a -> 1, b -> x, c -> r|>, <|a -> 1, b -> z, c -> p|>, <|a -> 1, b -> z, c -> r|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>, <|a->2|>}, {<|a->3,c->p|>, <|a->4, d->1|>}, a, \"Outer\"]",
          "{<|a -> 1, b -> x, c -> Missing[\"NotAvailable\"], d -> Missing[\"NotAvailable\"]|>, <|a -> 2, b -> Missing[\"NotAvailable\"], c -> Missing[\"NotAvailable\"], d -> Missing[\"NotAvailable\"]|>, <|a -> 3, b -> Missing[\"NotAvailable\"], c -> p, d -> Missing[\"NotAvailable\"]|>, <|a -> 4, b -> Missing[\"NotAvailable\"], c -> Missing[\"NotAvailable\"], d -> 1|>}");
    CHECK("JoinAcross[{<|a->2,b->3|>, <|a->3, y->1|>}, {<|a->1,c->1|>}, a, \"Left\"]",
          "{<|a -> 2, b -> 3, y -> Missing[\"NotAvailable\"], c -> Missing[\"Unmatched\"]|>, <|a -> 3, b -> Missing[\"NotAvailable\"], y -> 1, c -> Missing[\"Unmatched\"]|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>}, {<|e->2,c->p|>}, Key[a]->Key[e], \"Outer\"]",
          "{<|a -> 1, b -> x, e -> Missing[\"Unmatched\"], c -> Missing[\"Unmatched\"]|>, <|a -> Missing[\"Unmatched\"], b -> Missing[\"Unmatched\"], e -> 2, c -> p|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>}, {<|a->1,b->p, c->1|>}, a, KeyCollisionFunction -> Right]",
          "{<|a -> 1, b -> p, c -> 1|>}");
    CHECK("JoinAcross[{<|a->1,b->x|>}, {<|a->1,b->p, c->1|>}, a, KeyCollisionFunction -> Function[k, {k[1], k[2]}]]",
          "{<|a -> 1, b[1] -> x, b[2] -> p, c -> 1|>}");
    CHECK("JoinAcross[{<|a->2,b->3|>}, {<|a->1,b->9|>}, a, \"Right\"]",
          "{<|a -> 1, b -> Missing[\"Unmatched\"]|>}");
    CHECK("JoinAcross[{<|a->1|>}, {<|a->1|>}, a, \"Foo\"]",
          "JoinAcross[{<|a -> 1|>}, {<|a -> 1|>}, a, \"Foo\"]");
    CHECK("JoinAcross[{a->1}, {<|a->1|>}, a]",
          "JoinAcross[{a -> 1}, {<|a -> 1|>}, a]");
}

static void test_transpose(void) {
    CHECK("Transpose[<|a-><|x->1,y->2|>,b-><|x->3,y->4|>|>]",
          "<|x -> <|a -> 1, b -> 3|>, y -> <|a -> 2, b -> 4|>|>");
    CHECK("Transpose[{<|x->1,y->2|>,<|x->3,y->4|>}]",
          "{<|x -> 1, y -> 2|>, <|x -> 3, y -> 4|>}");
    CHECK("Transpose[<|a-><|x->1,y->2|>,b-><|y->4,x->3|>|>]",
          "Transpose[<|a -> <|x -> 1, y -> 2|>, b -> <|y -> 4, x -> 3|>|>]");
    CHECK("Transpose[<|a-><||>, b-><||>|>]",
          "<||>");
    CHECK("Transpose[{{1,2},{3,4}}]",
          "{{1, 3}, {2, 4}}");
    CHECK("Transpose[<|a->{1,2}|>]",
          "Transpose[<|a -> {1, 2}|>]");
}

static void test_applyto(void) {
    CHECK("s = {1,2,3}; ApplyTo[s[[2]], f]; s",
          "{1, f[2], 3}");
    CHECK("as = <|a->1|>; ApplyTo[as[a], f]; as",
          "<|a -> f[1]|>");
    CHECK("u = 5; {ApplyTo[u, #^2&], u}",
          "{25, 25}");
    CHECK("ApplyTo[3, f]",
          "ApplyTo[3, f]");
    CHECK("ApplyTo[w, f]",
          "ApplyTo[w, f]");
    CHECK("z = <|a->1|>; ApplyTo[z, Append[#, b->2]&]; z",
          "<|a -> 1, b -> 2|>");
}

static void test_associationcomap(void) {
    CHECK("AssociationComap[{f,g}, 2]",
          "<|f -> f[2], g -> g[2]|>");
    CHECK("AssociationComap[{}, 2]",
          "<||>");
    CHECK("AssociationComap[{f,f}, 2]",
          "<|f -> f[2]|>");
    CHECK("AssociationComap[f, 2]",
          "AssociationComap[f, 2]");
}

static void test_attributes(void) {
    CHECK("Attributes[ApplyTo]",
          "{HoldFirst, Protected}");
    CHECK("Attributes[Key]",
          "{Protected}");
    /* Mathematica adds ReadProtected, which Mathilda does not implement. */
    CHECK("Attributes[Missing]",
          "{Protected}");
    CHECK("Attributes[JoinAcross]",
          "{Protected}");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_threading);
    TEST(test_listable_regression);
    TEST(test_packed_and_ndarray);
    TEST(test_keyintersection);
    TEST(test_keycomplement);
    TEST(test_missingq);
    TEST(test_discard);
    TEST(test_countdistinct);
    TEST(test_subsetq);
    TEST(test_splice);
    TEST(test_keys_values);
    TEST(test_keyunion);
    TEST(test_associationmap);
    TEST(test_positionindex);
    TEST(test_merge);
    TEST(test_association);
    TEST(test_normal);
    TEST(test_counts_tally);
    TEST(test_deletemissing);
    TEST(test_lookup);
    TEST(test_joinacross);
    TEST(test_transpose);
    TEST(test_applyto);
    TEST(test_associationcomap);
    TEST(test_attributes);

    printf("All assoc_ops tests passed (%d checks).\n", g_checks);
    return 0;
}
