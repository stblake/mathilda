#include "expand.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

void test_expand() {
    run_test("Expand[(x+3)(x+2)]", "Plus[6, Times[5, x], Power[x, 2]]");
    run_test("Expand[(x+y)^2 (x-y)^2]", "Plus[Power[x, 4], Times[-2, Times[Power[x, 2], Power[y, 2]]], Power[y, 4]]");
    run_test("Expand[(x+y)^5 - (x^5+y^5)]", "Plus[Times[5, Times[Power[x, 4], y]], Times[10, Times[Power[x, 3], Power[y, 2]]], Times[10, Times[Power[x, 2], Power[y, 3]]], Times[5, Times[x, Power[y, 4]]]]");
    run_test("Expand[(f[Subscript[x, 1]] + f[Subscript[x, 2]])^3]", "Plus[Power[f[Subscript[x, 1]], 3], Power[f[Subscript[x, 2]], 3], Times[3, Times[Power[f[Subscript[x, 1]], 2], f[Subscript[x, 2]]]], Times[3, Times[Power[f[Subscript[x, 2]], 2], f[Subscript[x, 1]]]]]");
    run_test("Expand[(Sin[x] + Cos[x])^2]", "Plus[Power[Cos[x], 2], Power[Sin[x], 2], Times[2, Times[Cos[x], Sin[x]]]]");
    run_test("Expand[Sqrt[(1+x)^2]]", "Power[Power[Plus[1, x], 2], Rational[1, 2]]");
    run_test("Expand[1 < (x+y)^2 < 2]", "Inequality[1, Less, Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]], Less, 2]");
    
    // Pattern-based Expand
    run_test("Expand[(x+1)^2 + (y+1)^2, x]", "Plus[1, Times[2, x], Power[x, 2], Power[Plus[1, y], 2]]");
    run_test("Expand[(x+1)^2 + (y+1)^2, y]", "Plus[1, Power[Plus[1, x], 2], Times[2, y], Power[y, 2]]");
    run_test("Expand[(x+1)^2 + (y+1)^2, _Symbol]", "Plus[2, Times[2, x], Power[x, 2], Times[2, y], Power[y, 2]]");
    run_test("Expand[(x^2+1)^2 + (y+1)^2, x^2]", "Plus[1, Times[2, Power[x, 2]], Power[x, 4], Power[Plus[1, y], 2]]");

    // Two-arg Expand across a product: pattern-free factors stay unexpanded and
    // ride along as an atomic coefficient (regression for the (a+b) distribution bug).
    run_test("Expand[(a+b)(x+y)^2, x]", "Plus[Times[Plus[a, b], Power[x, 2]], Times[2, Times[Plus[a, b], x, y]], Times[Plus[a, b], Power[y, 2]]]");
    run_test("Expand[(x+1)^2 (y+1)^2, x]", "Plus[Power[Plus[1, y], 2], Times[2, Times[x, Power[Plus[1, y], 2]]], Times[Power[x, 2], Power[Plus[1, y], 2]]]");
    // A pattern-free numeric coefficient still distributes.
    run_test("Expand[3 (x+y), x]", "Plus[Times[3, x], Times[3, y]]");
}

/* FLINT-accelerated large polynomial expansion and Overflow[] gating. */
void test_expand_flint_overflow() {
    /* Dense univariate power: the multinomial term count C(344,3) ~ 6.7M would
     * once have refused this, but the Newton-box estimate sees a degree-5797
     * univariate (<= 5798 terms) and FLINT expands it. Assert structural facts
     * rather than the 5758-term fullform. */
    run_test("Length[Expand[(1 + x + 5 x^3 + 8 x^17)^341]]", "5758");
    run_test("Exponent[Expand[(1 + x + 5 x^3 + 8 x^17)^341], x]", "5797");
    run_test("Coefficient[Expand[(1 + x + 5 x^3 + 8 x^17)^341], x, 0]", "1");

    /* FLINT result matches the classical binomial expansion exactly. */
    run_test("Expand[(x + 2)^60] === Sum[Binomial[60, k] x^k 2^(60 - k), {k, 0, 60}]", "True");

    /* Genuine high-dimensional blow-up: never silently declined, yields
     * Overflow[], which propagates through Plus. */
    run_test("Expand[(a + b + c + d + e + f + g)^60]", "Overflow[]");
    run_test("Expand[(a + b + c + d + e + f + g)^60 + x^2]", "Overflow[]");

    /* Mixed poly/non-poly: the polynomial part expands, the transcendental
     * kernel is preserved. */
    run_test("Expand[(1 + x)^2 + Sin[y]^2]",
             "Plus[1, Times[2, x], Power[x, 2], Power[Sin[y], 2]]");

    /* Partitioned Times: a non-polynomial factor forces the whole-product FLINT
     * decline, but the polynomial sub-product is still expanded via FLINT and
     * only the transcendental factor is distributed generically. Verify the
     * result equals the reference (order-independent), including a large
     * intermediate that would be slow to distribute term-by-term. */
    run_test("Expand[Log[x] (1 + y)^20 (1 - y)^20 - Log[x] (1 - y^2)^20] === 0", "True");
    run_test("Expand[Sin[x] (p + q) (r + s)] === Expand[Sin[x] Expand[(p + q) (r + s)]]", "True");
    run_test("Expand[Sin[x] (p + q) (r + s)]",
             "Plus[Times[p, r, Sin[x]], Times[q, r, Sin[x]], Times[p, s, Sin[x]], Times[q, s, Sin[x]]]");
}

/* ExpandAll: expand products and integer powers in EVERY part of expr. */
void test_expand_all() {
    /* Expands inside a transcendental argument AND expands the denominator,
     * where a top-level Expand reaches neither. */
    run_test("ExpandAll[1/(1+x)^3 + Sin[(1+x)^3]]",
             "Plus[Power[Plus[1, Times[3, x], Times[3, Power[x, 2]], Power[x, 3]], -1], "
             "Sin[Plus[1, Times[3, x], Times[3, Power[x, 2]], Power[x, 3]]]]");
    /* Contrast: plain Expand leaves both untouched. */
    run_test("Expand[1/(1+x)^3 + Sin[(1+x)^3]]",
             "Plus[Power[Plus[1, x], -3], Sin[Power[Plus[1, x], 3]]]");

    /* Numerator and denominator of a rational both expand, and the numerator
     * distributes over the expanded denominator. */
    run_test("ExpandAll[(x+z)^2/(x+y)^2]",
             "Plus[Times[Power[x, 2], Power[Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]], -1]], "
             "Times[2, Times[x, Power[Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]], -1], z]], "
             "Times[Power[Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]], -1], Power[z, 2]]]");

    /* Expand inside an exponent. */
    run_test("ExpandAll[E^(I a (t-b))]",
             "Power[E, Plus[Times[Complex[0, -1], Times[a, b]], Times[Complex[0, 1], Times[a, t]]]]");

    /* Go into compound (non-symbol) heads. */
    run_test("ExpandAll[((1+a) (1+b))[x]]",
             "Plus[1, a, b, Times[a, b]][x]");

    /* Two-arg ExpandAll: expand only subexpressions containing x. (y+z)^2 has
     * no x and stays factored. */
    run_test("ExpandAll[(f[(x+y)^2] + g[(y+z)^2])^2, x]",
             "Plus[Power[f[Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]]], 2], "
             "Power[g[Power[Plus[y, z], 2]], 2], "
             "Times[2, Times[f[Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2]]], g[Power[Plus[y, z], 2]]]]]");

    /* Threads over lists, expanding a numerator power and a denominator power. */
    run_test("ExpandAll[{(1+x)^2, 1/(1+y)^2}]",
             "List[Plus[1, Times[2, x], Power[x, 2]], "
             "Power[Plus[1, Times[2, y], Power[y, 2]], -1]]");

    /* Threads over an equation, expanding both sides. */
    run_test("ExpandAll[(1+x)^2 == 1/(1+y)^2]",
             "Equal[Plus[1, Times[2, x], Power[x, 2]], "
             "Power[Plus[1, Times[2, y], Power[y, 2]], -1]]");

    /* Nested denominator inside a transcendental head. */
    run_test("ExpandAll[Log[1/(a+b)^2]]",
             "Log[Power[Plus[Power[a, 2], Times[2, Times[a, b]], Power[b, 2]], -1]]");

    /* Atoms and already-expanded input are returned unchanged. */
    run_test("ExpandAll[x]", "x");
    run_test("ExpandAll[a + b + c]", "Plus[a, b, c]");
    run_test("ExpandAll[7]", "7");

    /* Pattern that matches nothing: whole expression left untouched. */
    run_test("ExpandAll[(1+x)^2, y]", "Power[Plus[1, x], 2]");

    /* Deeply nested: product inside a sum inside a function argument. */
    run_test("ExpandAll[f[(a+b)(c+d)]]",
             "f[Plus[Times[a, c], Times[b, c], Times[a, d], Times[b, d]]]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_expand);
    TEST(test_expand_flint_overflow);
    TEST(test_expand_all);

    printf("All expand tests passed!\n");
    return 0;
}
