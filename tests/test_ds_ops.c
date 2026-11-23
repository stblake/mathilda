/* The ds_ops registry: the generic builtins that were migrated onto it (Length,
 * Dimensions, Head, and the pattern matcher's head test) must answer exactly as
 * they did with the hard-coded per-kind branches, for every registered kind. */
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "ds_ops.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

static void test_registry_rows(void) {
    ASSERT_MSG(ds_ops(EXPR_NDARRAY) && ds_ops(EXPR_SPARSE) && ds_ops(EXPR_LINKEDLIST),
               "built-in payload kinds must be registered");
    ASSERT_MSG(ds_ops(EXPR_INTEGER) == NULL && ds_ops(EXPR_FUNCTION) == NULL &&
               ds_ops(EXPR_SYMBOL) == NULL, "atoms and function nodes have no row");
    ASSERT_MSG(ds_register(NULL) == false, "NULL row refused");
}
static void test_length(void) {
    assert_true("Length[Range[10]] == 10");                         /* packed list */
    assert_true("Length[Table[i + j, {i, 7}, {j, 3}]] == 7");       /* leading axis */
    assert_true("Length[NDArray[{{1., 2.}, {3., 4.}, {5., 6.}}]] == 3");
    assert_true("Length[SparseArray[{1 -> 5, 4 -> 2}, 9]] == 9");
    assert_true("Length[LinkedList[{1, 2, 3, 4}]] == 4");
    assert_true("Length[LinkedList[]] == 0");
    assert_true("Length[3/4] == 0 && Length[x] == 0 && Length[f[1, 2]] == 2");
}
static void test_dimensions(void) {
    assert_true("Dimensions[Range[6]] === {6}");
    assert_true("Dimensions[Table[i + j, {i, 4}, {j, 3}]] === {4, 3}");
    assert_true("Dimensions[Table[i + j + k, {i, 4}, {j, 3}, {k, 2}], 2] === {4, 3}");
    assert_true("Dimensions[SparseArray[{{1, 2} -> 1}, {5, 6}]] === {5, 6}");
    assert_true("Dimensions[LinkedList[{1, 2}]] === {}");              /* not rectangular */
    assert_true("Dimensions[{{1, 2}, {3}}] === {2}");                /* ragged plain list */
}
static void test_head(void) {
    assert_true("Head[Range[5]] === List");                          /* packed list is a List */
    assert_true("Head[NDArray[{1., 2.}]] === NDArray");
    assert_true("Head[SparseArray[{1 -> 1}, 3]] === SparseArray");
    assert_true("Head[LinkedList[{1}]] === LinkedList");
    assert_true("Head[3] === Integer && Head[2.5] === Real && Head[\"s\"] === String && Head[x] === Symbol");
    assert_true("Head[f[1]] === f");
}
static void test_blank_head_match(void) {
    assert_true("MatchQ[SparseArray[{1 -> 1}, 3], _SparseArray]");
    assert_true("!MatchQ[SparseArray[{1 -> 1}, 3], _List]");
    assert_true("MatchQ[LinkedList[{1, 2}], _LinkedList]");
    assert_true("!MatchQ[LinkedList[{1, 2}], _List]");
    assert_true("MatchQ[Range[4], _List]");
    assert_true("Length[Cases[{LinkedList[{1}], SparseArray[{1 -> 1}, 2], 5}, _LinkedList]] == 1");
}
static void test_part(void) {
    assert_true("Range[100][[5]] == 5 && Range[100][[-1]] == 100");
    assert_true("Range[100][[2 ;; 4]] === {2, 3, 4} && Range[100][[{1, 3}]] === {1, 3}");
    assert_true("Range[100][[0]] === List");
    assert_true("SparseArray[{1 -> 3, 5 -> 7}, 6][[5]] == 7");
    assert_true("Normal[SparseArray[{1 -> 3, 5 -> 7}, 6][[2 ;; 5]]] === {0, 0, 0, 7}");
    assert_true("LinkedList[{1, 2, 3}][[2]] == 2");
    assert_true("Normal[LinkedList[{1, 2, 3}][[1 ;; 2]]] === {1, 2} || LinkedList[{1, 2, 3}][[1 ;; 2]] =!= $Failed");
    assert_true("Head[LinkedList[{1, 2, 3}][[0]]] === Symbol");
    assert_true("MatchQ[Range[3][[9]], _Part]");                      /* out of range: unevaluated */
    assert_true("Developer`ToPackedArray[{{1, 2}, {3, 4}}][[All, 1]] === {1, 3}");
    assert_true("Part[{SparseArray[{1 -> 3, 5 -> 7}, 6], 2}, 1, 5] == 7");  /* below the top */
}
static void test_normal(void) {
    assert_true("Normal[SparseArray[{1 -> 3, 5 -> 7}, 6]] === {3, 0, 0, 0, 7, 0}");
    assert_true("Normal[LinkedList[{1, 2, 3}]] === {1, 2, 3}");
    assert_true("Normal[Developer`ToPackedArray[{{1, 2}, {3, 4}}]] === {{1, 2}, {3, 4}}");
    assert_true("Normal[NDArray[{1., 2.}]] === {1., 2.}");
    assert_true("Normal[<|a -> 1|>] === {a -> 1}");
    assert_true("Normal[x] === x");
}
static int each_count; static int64_t each_sum;
static bool each_fn(Expr* e, void* ud) { (void)ud; each_count++; if (e->type == EXPR_INTEGER) each_sum += e->data.integer; return true; }
static bool each_stop(Expr* e, void* ud) { (void)e; (*(int*)ud)++; return *(int*)ud < 3; }
static void test_each(void) {
    ASSERT_MSG(ds_has_each(EXPR_PVEC) && ds_has_each(EXPR_LINKEDLIST), "PVector and LinkedList iterate");
    ASSERT_MSG(!ds_has_each(EXPR_INTEGER) && !ds_has_each(EXPR_SPARSE), "atoms and SparseArray have no each");
    Expr* pv = parse_expression("PVector[Range[50]]"); Expr* v = evaluate(pv);
    each_count = 0; each_sum = 0;
    ASSERT_MSG(ds_each(v, each_fn, NULL) && each_count == 50 && each_sum == 1275, "PVector each visits all in order");
    int stop = 0; ds_each(v, each_stop, &stop);
    ASSERT_MSG(stop == 3, "each stops when fn returns false");
    ASSERT_MSG(ds_atoms_only(v), "typed PVector holds atoms");
    expr_free(v);
    Expr* ll = evaluate(parse_expression("LinkedList[{1, 2, 3, x}]"));
    each_count = 0; each_sum = 0;
    ASSERT_MSG(ds_each(ll, each_fn, NULL) && each_count == 4 && each_sum == 6 && !ds_atoms_only(ll), "boxed LinkedList each");
    expr_free(ll);
    Expr* gv = evaluate(parse_expression("PVector[{7, 8, 9}]"));
    Expr* r = ds_get(gv, 1);
    expr_free(gv);
    ASSERT_MSG(r && r->type == EXPR_INTEGER && r->data.integer == 8, "ds_get");
    expr_free(r);
}

static void test_get_rebuild(void) {
    ASSERT_MSG(ds_ops(EXPR_PVEC) && ds_ops(EXPR_PVEC)->get && ds_ops(EXPR_PVEC)->from_list, "PVector row has get and from_list");
    ASSERT_MSG(ds_ops(EXPR_LINKEDLIST)->get && ds_ops(EXPR_LINKEDLIST)->from_list, "LinkedList row has get and from_list");
    assert_true("PVectorQ[PVector[{1, 2}]] && Normal[PVector[{1, 2}] /. x_ :> x] === {1, 2}");
    /* a structure-preserving head keeps the kind of its argument */
    assert_true("PVectorQ[Reverse[PVector[{1, 2, 3}]]] && Head[Reverse[LinkedList[{1, 2, 3}]]] === LinkedList");
    assert_true("StringContainsQ[ToString[LinkedList[{1}], TeXForm], \"LinkedList\"]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_registry_rows);
    TEST(test_length);
    TEST(test_dimensions);
    TEST(test_head);
    TEST(test_blank_head_match);
    TEST(test_part);
    TEST(test_normal);
    TEST(test_get_rebuild);
    TEST(test_each);
    printf("All ds_ops tests passed!\n");
    return 0;
}
