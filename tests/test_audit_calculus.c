/* Builtin-audit regressions: calculus (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* D#1: a number or numeric constant is not a differentiation variable. */
static void test_d_numeric_variable(void) {
    audit_chk("D#1", "Check[D[x^2, 2], fired, D::ivar]", "fired");
    audit_chk("D#1", "Quiet[D[x^2, 2]]", "D[x^2, 2]");
    audit_chk("D#1", "Quiet[D[Pi^x, Pi]]", "D[Pi^x, Pi]");
    audit_chk("D#1", "Quiet[D[x^2, I]]", "D[x^2, I]");
    audit_chk("D#1", "Quiet[D[x^2, x, 0]]", "D[x^2, x, 0]");
    audit_chk("D#1", "Quiet[D[x, 1.5]]", "D[x, 1.5]");
    audit_chk("D#1", "Quiet[D[x^2, {2, 1}]]", "D[x^2, {2, 1}]");
    /* ordinary specs unchanged */
    audit_chk("D#1", "D[x^2, {x, 0}]", "x^2");
    audit_chk("D#1", "D[x^2 y, x, y]", "2 x");
    audit_chk("D#1", "D[x^3, {x, 2}]", "6 x");
}

/* Dt#1: Constants -> {a} was read as a differentiation variable. */
static void test_dt_constants_option(void) {
    audit_chk("Dt#1", "Dt[x^2, x, Constants -> {a}]", "2 x");
    audit_chk("Dt#1", "Dt[a x^2, x, Constants -> {a}]", "2 a x");
    audit_chk("Dt#1", "Dt[a x, Constants -> a]", "a Dt[x]");
    audit_chk("Dt#1", "Dt[a x^2 + y, x, Constants -> {a}]", "2 a x + Dt[y, x]");
    audit_chk("Dt#1", "Dt[x^2 y, x]", "2 x y + x^2 Dt[y, x]");
}

/* Integrate#2: an indefinite variable among several specs. */
static void test_integrate_mixed_specs(void) {
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x, x, y] - x^2 y/2]");
    audit_true("Integrate#2", "PossibleZeroQ[D[Integrate[x, x, y], x, y] - x]");
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x y, {x, 0, 1}, y] - y^2/4]");
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x y, x, {y, 0, 2}] - x^2]");
    audit_true("Integrate#2", "PossibleZeroQ[D[Integrate[Sin[x] Cos[y], x, y], x, y] - Sin[x] Cos[y]]");
    audit_chk("Integrate#2", "Integrate[x y, {x, 0, 1}, {y, 0, 2}]", "1");
}

/* Residue#1: a logarithmic branch point has no residue. */
static void test_residue_branch_point(void) {
    audit_chk("Residue#1", "Residue[Log[z]/z, {z, 0}]", "Residue[Log[z]/z, {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[z], {z, 0}]", "Residue[Log[z], {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[z]/z^2, {z, 0}]", "Residue[Log[z]/z^2, {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[1 + z]/z^2, {z, 0}]", "1");
    audit_chk("Residue#1", "Residue[Exp[z]/z^2, {z, 0}]", "1");
    audit_chk("Residue#1", "Residue[1/(z^2 + 1), {z, I}]", "-1/2*I");
}

/* LaplaceTransform#1: the transform integrates from 0-, so DiracDelta[t] counts in full. */
static void test_laplace_delta_origin(void) {
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[t], t, s]", "1");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[2 t], t, s]", "1/2");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[Cos[t] DiracDelta[t], t, s]", "1");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[E^(3 t) DiracDelta[t], t, s]", "1");
    audit_true("LaplaceTransform#1",
               "PossibleZeroQ[LaplaceTransform[3 DiracDelta[t] + t, t, s] - (3 + 1/s^2)]");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[t - 2], t, s]", "E^(-2 s)");
}

/* Maximize#1: x == I has no real solution. */
static void test_maximize_complex_equality(void) {
    audit_chk("Maximize#1", "Check[Maximize[{-x^2, x == I}, x], fired, Maximize::infeas]", "fired");
    audit_chk("Maximize#1", "Quiet[Maximize[{-x^2, x == I}, x]]", "{-Infinity, {x -> Indeterminate}}");
    audit_chk("Maximize#1", "Quiet[Minimize[{x^2, x == I}, x]]", "{Infinity, {x -> Indeterminate}}");
    audit_chk("Maximize#1", "Maximize[{-x^2, x == 1}, x]", "{-1, {x -> 1}}");
    audit_chk("Maximize#1", "Minimize[{x^2 + y^2, x == 1}, {x, y}]", "{1, {x -> 1, y -> 0}}");
}

/* PolynomialQuotientRemainder#1: p == Q q + R must hold. */
static void test_pqr_nonpolynomial_terms(void) {
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[1/x + x^2, x + 1, x]", "{-1 + x, 1 + 1/x}");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[Sqrt[x] + x^2, x + 1, x]", "{-1 + x, 1 + Sqrt[x]}");
    audit_true("PolynomialQuotientRemainder#1",
               "Module[{qr = PolynomialQuotientRemainder[1/x + x^2, x + 1, x]}, "
               "PossibleZeroQ[qr[[1]] (x + 1) + qr[[2]] - (1/x + x^2)]]");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[x^2 + 1, x + 1, x]", "{-1 + x, 2}");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[a x^2 + 1, x + 1, x]", "{-a + a x, 1 + a}");
}

int main(void) {
    audit_init();
    TEST(test_d_numeric_variable);
    TEST(test_dt_constants_option);
    TEST(test_integrate_mixed_specs);
    TEST(test_residue_branch_point);
    TEST(test_laplace_delta_origin);
    TEST(test_maximize_complex_equality);
    TEST(test_pqr_nonpolynomial_terms);
    return audit_finish("audit_calculus");
}
