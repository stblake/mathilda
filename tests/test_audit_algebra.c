/* Builtin-audit regressions: algebra (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* AlgebraicIntegerQ#1 (high): past to_qqbar's degree cap (120) the decline
 * became a definite False. Ring closure now decides sums, products and
 * positive (rational) powers of integers outright -- 2^(1/k) is a root of the
 * monic x^k - 2 -- and a cap decline is undecided (unevaluated), never False. */
static void test_algebraicintegerq_high_degree(void) {
    audit_true("AlgebraicIntegerQ#1",
               "Table[AlgebraicIntegerQ[2^(1/k)], {k, {50, 100, 200, 300, 500}}] === "
               "{True, True, True, True, True}");
    audit_true("AlgebraicIntegerQ#1",
               "Table[AlgebraicIntegerQ[Sum[Sqrt[Prime[k]], {k, n}]], {n, 4, 8}] === "
               "{True, True, True, True, True}");
    audit_true("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[2^(1/1000) + 1]");
    /* A capped non-closure case is left undecided, not answered False. */
    audit_true("AlgebraicIntegerQ#1",
               "Head[AlgebraicIntegerQ[1/2^(1/200)]] === AlgebraicIntegerQ");
    /* Decisions under the cap are unchanged. */
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[Sqrt[2]/2]", "False");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[(1 + Sqrt[5])/2]", "True");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[1/2]", "False");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[Pi]", "False");
}

/* IrreduciblePolynomialQ#1 (high): Extension -> All on a multivariate input
 * factored over Q(i) only. Now exact: irreducible over Q, then one Factor over
 * Q(t), t a root of a squarefree specialisation (a smooth point fixes the one
 * absolute factor through it, so that factor is defined over Q(t)). */
static void test_irreduciblepolynomialq_absolute(void) {
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 - 2 y^2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + x y + y^2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^3 - 2 y^3, Extension -> All]", "False");
    /* x^2 y^2 - 2 = (x y - Sqrt[2])(x y + Sqrt[2]) */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 y^2 - 2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^2, Extension -> All]", "False");
    /* Absolutely irreducible curves and a rank-3 quadric stay True. */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^3, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^4 + y^4 + 1, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^2 + z^2, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x y - 1, Extension -> All]", "True");
    /* Over Q (no extension) is unchanged. */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 - 2 y^2]", "True");
}

int main(void) {
    audit_init();
    TEST(test_algebraicintegerq_high_degree);
    TEST(test_irreduciblepolynomialq_absolute);
    return audit_finish("audit_algebra");
}
