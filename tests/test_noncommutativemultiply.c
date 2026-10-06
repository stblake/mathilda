/*
 * test_noncommutativemultiply.c — unit tests for NonCommutativeMultiply (**)
 *
 * NonCommutativeMultiply is an associative (Flat) but non-commutative
 * generalized multiplication with attributes {Flat, OneIdentity, Protected}
 * and NO automatic simplification beyond flattening. Expected behavior and
 * expected strings are taken from the Mathematica reference and verified
 * against the implementation.
 */
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "attr.h"      /* get_attributes, ATTR_FLAT/ONEIDENTITY/PROTECTED */

/* ---- 1. Parsing: a**b**c builds the n-ary NonCommutativeMultiply[a,b,c] ---- */
void test_ncm_parse(void) {
    /* Result of evaluation, in FullForm. */
    assert_eval_eq("a ** b ** c", "NonCommutativeMultiply[a, b, c]", 1);
    /* Under Hold (no evaluation) the parser already produces the flat n-ary
     * form, exactly as Mathematica's parser does. */
    assert_eval_eq("Hold[a ** b ** c]", "Hold[NonCommutativeMultiply[a, b, c]]", 1);
    /* Two-argument base case. */
    assert_eval_eq("a ** b", "NonCommutativeMultiply[a, b]", 1);
}

/* ---- 2. Flat / associativity: a**(b**c) and (a**b)**c flatten alike ---- */
void test_ncm_flat(void) {
    assert_eval_eq("a ** (b ** c)", "NonCommutativeMultiply[a, b, c]", 1);
    assert_eval_eq("(a ** b) ** c", "NonCommutativeMultiply[a, b, c]", 1);
    /* Nested head form flattens via the Flat attribute (no builtin involved). */
    assert_eval_eq("NonCommutativeMultiply[a, NonCommutativeMultiply[b, c]]",
                   "NonCommutativeMultiply[a, b, c]", 1);
    assert_eval_eq("NonCommutativeMultiply[NonCommutativeMultiply[a, b], c]",
                   "NonCommutativeMultiply[a, b, c]", 1);
    /* Operations are associative. */
    assert_eval_eq("a ** (b ** c) == (a ** b) ** c", "True", 0);
}

/* ---- 3. Non-commutative: args are NOT reordered ---- */
void test_ncm_noncommutative(void) {
    assert_eval_eq("b ** a", "NonCommutativeMultiply[b, a]", 1);
    /* a**b and b**a are structurally distinct. */
    assert_eval_eq("(a ** b) === (b ** a)", "False", 0);
    /* a**b == b**a stays unevaluated (compare to a*b == b*a, which is True). */
    assert_eval_eq("a ** b == b ** a", "a ** b == b ** a", 0);
    assert_eval_eq("{a*b == b*a, a**b == b**a}", "{True, a ** b == b ** a}", 0);
}

/* ---- 4. One argument stays unevaluated (unlike Times[a] -> a) ---- */
void test_ncm_one_arg(void) {
    assert_eval_eq("NonCommutativeMultiply[a]", "NonCommutativeMultiply[a]", 1);
    /* OutputForm of the one-arg form is the bracketed head form, not infix. */
    assert_eval_eq("NonCommutativeMultiply[a]", "NonCommutativeMultiply[a]", 0);
}

/* ---- 5. No automatic simplification rules ---- */
void test_ncm_no_simplify(void) {
    assert_eval_eq("0 ** a", "NonCommutativeMultiply[0, a]", 1);
    assert_eval_eq("1 ** a", "NonCommutativeMultiply[1, a]", 1);
    assert_eval_eq("{0 ** a, 1 ** a}", "{0 ** a, 1 ** a}", 0);
    /* Expand / Simplify / FullSimplify do not operate on NCM expressions.
     * (These three invoke the Simplify family, whose small residual leak is
     * pre-existing and unrelated to NonCommutativeMultiply — the pure NCM
     * paths valgrind clean at the core_init baseline.) */
    assert_eval_eq("Expand[(a + b) ** c - a*c]", "-a c + (a + b) ** c", 0);
    assert_eval_eq("Simplify[(a + b) ** c - a*c]", "-a c + (a + b) ** c", 0);
    assert_eval_eq("FullSimplify[(a + b) ** c - a*c]", "-a c + (a + b) ** c", 0);
}

/* ---- 6. Printing: infix with precedence-aware parenthesisation ---- */
void test_ncm_print(void) {
    assert_eval_eq("a ** b ** c", "a ** b ** c", 0);
    /* Plus (looser) is parenthesised as an NCM argument. */
    assert_eval_eq("(a + b) ** c", "(a + b) ** c", 0);
    /* NCM (tighter than Plus) is NOT parenthesised inside Plus. */
    assert_eval_eq("a ** b + c", "a ** b + c", 0);
}

/* ---- 7. Precedence: Plus < Times < Dot < NCM < Power ---- */
void test_ncm_precedence(void) {
    /* Times binds looser: a (b**c). */
    assert_eval_eq("a b ** c", "Times[a, NonCommutativeMultiply[b, c]]", 1);
    /* Dot binds looser: a . (b**c). */
    assert_eval_eq("a . b ** c", "Dot[a, NonCommutativeMultiply[b, c]]", 1);
    /* Power binds tighter: a ** (b^c). */
    assert_eval_eq("a ** b ^ c", "NonCommutativeMultiply[a, Power[b, c]]", 1);
}

/* ---- 8. Attributes {Flat, OneIdentity, Protected} ---- */
void test_ncm_attributes(void) {
    /* Language-level (printed alphabetically). */
    assert_eval_eq("Attributes[NonCommutativeMultiply]",
                   "{Flat, OneIdentity, Protected}", 0);
    /* C-API, bit level. */
    SymbolDef* def = symtab_get_def("NonCommutativeMultiply");
    ASSERT(def != NULL);
    uint32_t at = get_attributes("NonCommutativeMultiply");
    ASSERT((at & ATTR_FLAT)        != 0);
    ASSERT((at & ATTR_ONEIDENTITY) != 0);
    ASSERT((at & ATTR_PROTECTED)   != 0);
}

/* ---- 9. Pattern-matcher integration (Flat + OneIdentity in patterns) ---- */
void test_ncm_patterns(void) {
    /* A plain two-slot pattern matches the base case. */
    Expr* d1 = parse_expression("ncmH[p_ ** q_] := ncmHH[p, q]");
    Expr* r1 = evaluate(d1); expr_free(d1); if (r1) expr_free(r1);
    assert_eval_eq("ncmH[a ** b]", "ncmHH[a, b]", 0);
    /* On three factors, Flat grouping fills the single q_ slot with b**c. */
    assert_eval_eq("ncmH[a ** b ** c]", "ncmHH[a, b ** c]", 0);

    /* A BlankSequence in the first slot (the DOperator[L1__ ** L2_, ...] shape). */
    Expr* d2 = parse_expression("ncmG[x__ ** y_] := ncmNC[{x}, y]");
    Expr* r2 = evaluate(d2); expr_free(d2); if (r2) expr_free(r2);
    assert_eval_eq("ncmG[a ** b ** c]", "ncmNC[{a}, b ** c]", 0);
}

/* ---- 10. Memory exercise (run under valgrind --leak-check=full) ---- */
void test_ncm_memory(void) {
    const char* inputs[] = {
        "a ** b ** c",
        "a ** (b ** c)",
        "(a ** b) ** c",
        "NonCommutativeMultiply[a]",
        "NonCommutativeMultiply[a, NonCommutativeMultiply[b, c]]",
        "0 ** a",
        "1 ** a",
        "b ** a",
        "(a + b) ** c - a*c",
        "a ** b == b ** a",
        NULL
    };
    for (int i = 0; inputs[i]; i++) {
        Expr* e = parse_expression(inputs[i]);
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        expr_free(e);
        if (r) expr_free(r);
    }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_ncm_parse);
    TEST(test_ncm_flat);
    TEST(test_ncm_noncommutative);
    TEST(test_ncm_one_arg);
    TEST(test_ncm_no_simplify);
    TEST(test_ncm_print);
    TEST(test_ncm_precedence);
    TEST(test_ncm_attributes);
    TEST(test_ncm_patterns);
    TEST(test_ncm_memory);

    printf("All NonCommutativeMultiply tests passed!\n");
    return 0;
}
