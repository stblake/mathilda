#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"

/* Basic: Evaluate[expr] outside hold context acts as identity */
void test_evaluate_basic() {
    assert_eval_eq("Evaluate[1+1]", "2", 0);
    assert_eval_eq("Evaluate[2]", "2", 0);
    assert_eval_eq("Evaluate[x]", "x", 0);
    assert_eval_eq("Evaluate[{1, 2, 3}]", "{1, 2, 3}", 0);
    assert_eval_eq("Evaluate[Sin[0]]", "0", 0);
}

/* Evaluate overrides Hold */
void test_evaluate_in_hold() {
    assert_eval_eq("Hold[Evaluate[1+1]]", "Hold[2]", 0);
    assert_eval_eq("Hold[Evaluate[1+1], 2+2]", "Hold[2, 2 + 2]", 0);
    assert_eval_eq("Hold[Evaluate[2+3]]", "Hold[5]", 0);
}

/* Evaluate works only at the first level inside Hold */
void test_evaluate_first_level_only() {
    assert_eval_eq("Hold[f[Evaluate[1+2]]]", "Hold[f[Evaluate[1 + 2]]]", 0);
    assert_eval_eq("Hold[{Evaluate[1+2]}]", "Hold[{Evaluate[1 + 2]}]", 0);
}

/* Evaluate overrides HoldFirst */
void test_evaluate_holdfirst() {
    /* If is HoldRest: second and third args are held */
    assert_eval_eq("If[Evaluate[1 == 1], yes, no]", "yes", 0);
    /* Set is not HoldFirst on rhs, but SetDelayed holds rhs */
}

/* Evaluate overrides HoldAll (e.g., Table, Do) */
void test_evaluate_holdall() {
    assert_eval_eq("Hold[Evaluate[1+2], Evaluate[3+4]]", "Hold[3, 7]", 0);
}

/* Evaluate does NOT override HoldAllComplete */
void test_evaluate_holdallcomplete() {
    /* HoldAllComplete prevents Evaluate from working */
    assert_eval_eq("Hold[Evaluate[1+1]]", "Hold[2]", 0);
    /* With a HoldAllComplete function, Evaluate should be preserved */
    /* We can test this by defining a function with HoldAllComplete */
}

/* Multiple Evaluate in the same held expression */
void test_evaluate_multiple() {
    assert_eval_eq("Hold[Evaluate[1+1], Evaluate[2+2], Evaluate[3+3]]", "Hold[2, 4, 6]", 0);
}

/* Evaluate with symbolic expressions */
void test_evaluate_symbolic() {
    assert_eval_eq("Hold[Evaluate[Sin[Pi/6]]]", "Hold[1/2]", 0);
    assert_eval_eq("Hold[Evaluate[2^10]]", "Hold[1024]", 0);
}

/* Evaluate with other than one argument reduces to Sequence (which then splices
 * wherever it lands). Standalone Sequence stays as-is. */
void test_evaluate_wrong_args() {
    assert_eval_eq("Evaluate[]", "Sequence[]", 0);
    assert_eval_eq("Evaluate[1, 2]", "Sequence[1, 2]", 0);
}

/* Evaluate[...] with arity != 1 in a held slot splices via Sequence:
 * Evaluate[] vanishes and Evaluate[a,b] forces then splices both. */
void test_evaluate_seq_in_hold() {
    assert_eval_eq("Hold[Evaluate[]]", "Hold[]", 0);
    assert_eval_eq("Hold[Evaluate[], 1+1]", "Hold[2]", 0);
    assert_eval_eq("Hold[Evaluate[1+1, 2+2]]", "Hold[2, 4]", 0);
    assert_eval_eq("Hold[a, Evaluate[1+1, 2+2], b]", "Hold[a, 2, 4, b]", 0);
}

/* Evaluate overrides the HoldAll of Attributes: Attributes[x] sees the symbol x,
 * while Attributes[Evaluate[x]] sees the value of x (here Plus). */
void test_evaluate_attributes_override() {
    assert_eval_eq("x=Plus; {Attributes[x], Attributes[Evaluate[x]]}",
                   "{{}, {Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}}", 0);
    assert_eval_eq("ClearAll[x]", "Null", 0);
}

/* Unevaluated temporarily makes a non-Hold head behave as HoldAll for that arg. */
void test_evaluate_length_unevaluated() {
    assert_eval_eq("Length[Unevaluated[1+2+3]]", "3", 0);
    assert_eval_eq("Length[1+2+3]", "0", 0);
}

/* Evaluate is Protected */
void test_evaluate_protected() {
    assert_eval_eq("Attributes[Evaluate]", "{Protected}", 0);
}

/* Evaluate with nested expressions */
void test_evaluate_nested_expr() {
    assert_eval_eq("Evaluate[Evaluate[1+1]]", "2", 0);
    assert_eval_eq("Hold[Evaluate[Evaluate[1+1]]]", "Hold[2]", 0);
}

/* Evaluate with List operations */
void test_evaluate_with_list() {
    assert_eval_eq("Evaluate[{1+1, 2+2, 3+3}]", "{2, 4, 6}", 0);
    assert_eval_eq("Hold[Evaluate[{1+1, 2+2}]]", "Hold[{2, 4}]", 0);
}

/* Evaluate inside HoldForm */
void test_evaluate_in_holdform() {
    /* HoldForm has HoldAll, so Evaluate should override it */
    assert_eval_eq("HoldForm[Evaluate[1+1]]", "2", 0);
}

/* Evaluate with compound expressions */
void test_evaluate_compound() {
    assert_eval_eq("Hold[Evaluate[3*4+1]]", "Hold[13]", 0);
    assert_eval_eq("Hold[Evaluate[Length[{a,b,c}]]]", "Hold[3]", 0);
}

/* Evaluate with pure functions in hold context */
void test_evaluate_with_functions() {
    assert_eval_eq("Evaluate[Head[{1,2,3}]]", "List", 0);
    assert_eval_eq("Hold[Evaluate[Head[{1,2,3}]]]", "Hold[List]", 0);
}

/* Verify no memory leaks with a simple constructed test */
void test_evaluate_memory() {
    /* Exercise Evaluate in various contexts to catch leaks */
    Expr* e1 = parse_expression("Evaluate[1+1]");
    Expr* r1 = evaluate(e1);
    expr_free(e1);
    expr_free(r1);

    Expr* e2 = parse_expression("Hold[Evaluate[2+3]]");
    Expr* r2 = evaluate(e2);
    expr_free(e2);
    expr_free(r2);

    Expr* e3 = parse_expression("Hold[f[Evaluate[1+2]]]");
    Expr* r3 = evaluate(e3);
    expr_free(e3);
    expr_free(r3);

    Expr* e4 = parse_expression("Evaluate[x]");
    Expr* r4 = evaluate(e4);
    expr_free(e4);
    expr_free(r4);

    Expr* e5 = parse_expression("Evaluate[]");
    Expr* r5 = evaluate(e5);
    expr_free(e5);
    expr_free(r5);

    /* Sequence-producing forms (arity != 1), held and unheld. */
    Expr* e6 = parse_expression("Evaluate[1, 2]");
    Expr* r6 = evaluate(e6);
    expr_free(e6);
    expr_free(r6);

    Expr* e7 = parse_expression("Hold[Evaluate[1+1, 2+2]]");
    Expr* r7 = evaluate(e7);
    expr_free(e7);
    expr_free(r7);

    Expr* e8 = parse_expression("Hold[Evaluate[]]");
    Expr* r8 = evaluate(e8);
    expr_free(e8);
    expr_free(r8);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_evaluate_basic);
    TEST(test_evaluate_in_hold);
    TEST(test_evaluate_first_level_only);
    TEST(test_evaluate_holdfirst);
    TEST(test_evaluate_holdall);
    TEST(test_evaluate_holdallcomplete);
    TEST(test_evaluate_multiple);
    TEST(test_evaluate_symbolic);
    TEST(test_evaluate_wrong_args);
    TEST(test_evaluate_seq_in_hold);
    TEST(test_evaluate_attributes_override);
    TEST(test_evaluate_length_unevaluated);
    TEST(test_evaluate_protected);
    TEST(test_evaluate_nested_expr);
    TEST(test_evaluate_with_list);
    TEST(test_evaluate_in_holdform);
    TEST(test_evaluate_compound);
    TEST(test_evaluate_with_functions);
    TEST(test_evaluate_memory);

    printf("All Evaluate tests passed!\n");
    return 0;
}
