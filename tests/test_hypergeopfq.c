/* Tests for HypergeometricPFQ (src/hypergeopfq.c) and its Hypergeometric0F1/
 * 1F1/2F1 convenience heads, plus the z-derivative rule in deriv.c. */
#include "test_utils.h"
#include "core.h"
#include <stdio.h>

/* ---- structural rules ------------------------------------------------ */

static void test_origin(void) {
    /* pFq[a, b, 0] == 1 */
    assert_eval_eq("HypergeometricPFQ[{a1,a2,a3},{b1,b2,b3},0]", "1", 0);
    assert_eval_eq("HypergeometricPFQ[{},{},0]", "1", 0);
}

static void test_exact_stays_symbolic(void) {
    /* All-exact, non-reducible, non-terminating: stays unevaluated. */
    assert_eval_eq("HypergeometricPFQ[{1,1},{3,3,3},2]",
                   "HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2]", 0);
}

static void test_cancellation(void) {
    assert_eval_eq("HypergeometricPFQ[{a,b,c},{a,d,e},z]",
                   "HypergeometricPFQ[{b, c}, {d, e}, z]", 0);
    /* A non-positive integer common element must NOT cancel -- it makes the
     * series terminate to a polynomial instead. */
    assert_eval_eq("HypergeometricPFQ[{-2,a,b,c},{-2,d,e},z]",
        "1 + (a b c z)/(d e) + 1/2 (a (1 + a) b (1 + b) c (1 + c) z^2)/(d (1 + d) e (1 + e))",
        0);
}

static void test_termination(void) {
    /* Upper non-positive integer -> finite polynomial (symbolic z). */
    assert_eval_eq("HypergeometricPFQ[{-3,a},{b,c},x]",
        "1 - 3 (a x)/(b c) + 3 (a (1 + a) x^2)/(b (1 + b) c (1 + c)) - "
        "(a (1 + a) (2 + a) x^3)/(b (1 + b) (2 + b) c (1 + c) (2 + c))", 0);
    /* Degree-0 termination. */
    assert_eval_eq("HypergeometricPFQ[{0},{b},x]", "1", 0);
}

static void test_threading(void) {
    assert_eval_eq("HypergeometricPFQ[{1,2,3,4},{5,6,7},{0.1,0.3,0.5}]",
                   "{1.01164, 1.03627, 1.06296}", 0);
}

/* ---- symbolic reductions to elementary functions -------------------- */

static void test_reductions(void) {
    assert_eval_eq("HypergeometricPFQ[{},{},z]", "E^z", 0);            /* 0F0 */
    assert_eval_eq("HypergeometricPFQ[{a},{},z]", "(1 - z)^(-a)", 0);  /* 1F0 */
    assert_eval_eq("HypergeometricPFQ[{1},{},z]", "1/(1 - z)", 0);     /* geometric */
    assert_eval_eq("HypergeometricPFQ[{1,1},{2},z]", "-Log[1 - z]/z", 0); /* 2F1 */
    assert_eval_eq("HypergeometricPFQ[{},{1/2},z]", "Cosh[2 Sqrt[z]]", 0); /* 0F1 */
    assert_eval_eq("HypergeometricPFQ[{},{3/2},z]",
                   "(1/2 Sinh[2 Sqrt[z]])/Sqrt[z]", 0);                 /* 0F1 */
}

/* ---- convenience heads ---------------------------------------------- */

static void test_convenience_heads(void) {
    /* 0F1[b,z] == pFq[{},{b},z];  with b==1/2 reduces to Cosh. */
    assert_eval_eq("Hypergeometric0F1[1/2, z]", "Cosh[2 Sqrt[z]]", 0);
    /* 1F1[a,b,z] == pFq[{a},{b},z]; 1F1[1,2,z] == (E^z - 1)/z form via series,
     * but here just confirm it routes to pFq and stays symbolic for symbolic z. */
    assert_eval_eq("Hypergeometric1F1[a, b, 0]", "1", 0);
    /* 2F1[1,1,2,z] == -Log[1-z]/z. */
    assert_eval_eq("Hypergeometric2F1[1, 1, 2, z]", "-Log[1 - z]/z", 0);
}

/* ---- numeric: machine, complex, MPFR -------------------------------- */

static void test_numeric_machine(void) {
    assert_eval_eq("HypergeometricPFQ[{1,1},{3,3,3},2.]", "1.07893", 0);
    /* 0F0 numeric agrees with Exp. */
    assert_eval_eq("HypergeometricPFQ[{},{},2.]", "7.38906", 0);
}

static void test_numeric_complex(void) {
    assert_eval_startswith("HypergeometricPFQ[{I,I,I},{2,2,2},-1.0 I]", "0.870032");
}

static void test_numeric_mpfr(void) {
    /* 50-digit value from the Wolfram reference. */
    assert_eval_startswith("N[HypergeometricPFQ[{1,1,1},{3/2,3/2,3/2},10],50]",
                           "530.191888273625904388559616854440877927");
    /* 35-digit request: prefix must match. */
    assert_eval_startswith("N[HypergeometricPFQ[{1,1},{3,3,3},2],35]",
                           "1.0789");
}

static void test_convergence_gate(void) {
    /* p == q+1 with |z| >= 1 : not summed directly (stays unevaluated). */
    assert_eval_eq("HypergeometricPFQ[{10,10},{50},2.]",
                   "HypergeometricPFQ[{10, 10}, {50}, 2.0]", 0);
    /* p == q+1 with |z| < 1 : converges. */
    assert_eval_startswith("HypergeometricPFQ[{1,2},{3},0.5]", "1.54518");
}

/* ---- derivative ----------------------------------------------------- */

static void test_derivative(void) {
    assert_eval_eq("D[HypergeometricPFQ[{a1,a2},{b1,b2,b3},x],x]",
        "(a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)",
        0);
    /* Constant w.r.t. the differentiation variable -> 0. */
    assert_eval_eq("D[HypergeometricPFQ[{a},{b},y],x]", "0", 0);
}

int main(void) {
    core_init();
    TEST(test_origin);
    TEST(test_exact_stays_symbolic);
    TEST(test_cancellation);
    TEST(test_termination);
    TEST(test_threading);
    TEST(test_reductions);
    TEST(test_convenience_heads);
    TEST(test_numeric_machine);
    TEST(test_numeric_complex);
    TEST(test_numeric_mpfr);
    TEST(test_convergence_gate);
    TEST(test_derivative);
    printf("All HypergeometricPFQ tests passed.\n");
    return 0;
}
