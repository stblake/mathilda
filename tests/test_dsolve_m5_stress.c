/*
 * test_dsolve_m5_stress.c — randomized/parametrized stress tests for the M5
 * DSolve methods (Kovacic, Frobenius/PowerSeries).
 *
 * Kovacic is exercised by a FORWARD GENERATOR: pick a rational logarithmic
 * derivative omega, set r = omega' + omega^2, and feed y'' - r y == 0.  By
 * construction z = Exp[Integrate[omega]] is a Liouvillian (Case 1) solution, so
 * DSolve`Kovacic must solve it, and we verify by back-substitution.  This tests
 * the pure-polynomial (Exp), genuine-pole, and apparent-singularity paths.
 *
 * Frobenius is exercised over polynomial-coefficient ordinary-point equations:
 * the pinned power series must back-substitute to O[x]^k.
 *
 * As with the numeric solver stress tests, we assert the VALIDITY of whatever is
 * returned (Head===List then residual ~ 0), never a fixed solution form.  Keep
 * the whole run well under the 60 s test-harness alarm.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}

#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* Assert DSolve`Kovacic solves `eqn` (Head===List) and the residual `resid`
 * back-substitutes to zero.  The Head guard defeats a vacuous PossibleZeroQ on a
 * declined solve (where [[1]] stays symbolic). */
static void kovacic_ok(const char* eqn, const char* resid) {
    char buf[768];
    snprintf(buf, sizeof(buf), "Head[DSolve`Kovacic[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[(%s) /. DSolve`Kovacic[%s, y, x][[1]]]", resid, eqn);
    ASSERT_TRUE(buf);
}

/* ---- Kovacic Case 1 forward generator ----
 * Each family fixes omega and feeds the EXPANDED r = omega' + omega^2 (a clean
 * rational function), so a Liouvillian solution z = Exp[Integrate[omega]] is
 * guaranteed and the back-substitution residual stays in reduced form. */

/* omega = a x + b  ->  r = a^2 x^2 + 2ab x + (a + b^2)  (polynomial r, Exp sol). */
static void t_kov_stress_poly(void) {
    int as[] = {1, 2, 3};
    int bs[] = {-2, 0, 1, 2};
    for (size_t ai = 0; ai < 3; ai++)
        for (size_t bi = 0; bi < 4; bi++) {
            int a = as[ai], b = bs[bi];
            int A = a * a, B = 2 * a * b, C = a + b * b;
            char eqn[256], res[256];
            snprintf(eqn, sizeof(eqn), "y''[x] - (%d x^2 + %d x + %d) y[x] == 0", A, B, C);
            snprintf(res, sizeof(res), "D[y[x],{x,2}] - (%d x^2 + %d x + %d) y[x]", A, B, C);
            kovacic_ok(eqn, res);
        }
}

/* omega = a x + 1/x  ->  r = a^2 x^2 + 3a  (apparent singularity: pole cancels). */
static void t_kov_stress_apparent(void) {
    int as[] = {1, 2, 3};
    for (size_t ai = 0; ai < 3; ai++) {
        int a = as[ai];
        int A = a * a, C = 3 * a;
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn), "y''[x] - (%d x^2 + %d) y[x] == 0", A, C);
        snprintf(res, sizeof(res), "D[y[x],{x,2}] - (%d x^2 + %d) y[x]", A, C);
        kovacic_ok(eqn, res);
    }
}

/* omega = a x + c/x, c not in {0,1}  ->  r = a^2 x^2 + (a + 2ac) + (c^2-c)/x^2
 * (genuine order-2 pole). */
static void t_kov_stress_pole(void) {
    int as[] = {1, 2};
    int cs[] = {2, 3, -2};
    for (size_t ai = 0; ai < 2; ai++)
        for (size_t ci = 0; ci < 3; ci++) {
            int a = as[ai], c = cs[ci];
            int A = a * a, C0 = a + 2 * a * c, Cp = c * c - c;
            char eqn[256], res[256];
            snprintf(eqn, sizeof(eqn),
                     "y''[x] - (%d x^2 + %d + (%d)/x^2) y[x] == 0", A, C0, Cp);
            snprintf(res, sizeof(res),
                     "D[y[x],{x,2}] - (%d x^2 + %d + (%d)/x^2) y[x]", A, C0, Cp);
            kovacic_ok(eqn, res);
        }
}

/* omega = b/x only (Euler, roots 0 and 1-b) -> r = (b^2 - b)/x^2. */
static void t_kov_stress_euler(void) {
    int bs[] = {2, 3, -1, -2};
    for (size_t bi = 0; bi < 4; bi++) {
        int b = bs[bi];
        char eqn[256], res[256];
        snprintf(eqn, sizeof(eqn), "y''[x] - ((%d)/x^2) y[x] == 0", b * b - b);
        snprintf(res, sizeof(res), "D[y[x],{x,2}] - ((%d)/x^2) y[x]", b * b - b);
        kovacic_ok(eqn, res);
    }
}

/* ---- Frobenius / PowerSeries: ordinary-point polynomial-coefficient families ---- */

/* Pinned power series about 0; the truncated residual must be O[x]^k. */
static void frob_ordinary_ok(const char* pcoef, const char* qcoef) {
    char eqn[320], buf[640];
    snprintf(eqn, sizeof(eqn), "y''[x] + (%s) y'[x] + (%s) y[x] == 0", pcoef, qcoef);
    snprintf(buf, sizeof(buf), "Head[DSolve`PowerSeries[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[Module[{b = y[x] /. DSolve`PowerSeries[%s, y, x][[1]]}, "
             "D[b,{x,2}] + (%s) D[b,x] + (%s) b]]", eqn, pcoef, qcoef);
    ASSERT_TRUE(buf);
}

static void t_frob_stress_ordinary(void) {
    const char* ps[] = {"0", "x", "1 + x", "x^2", "2 x"};
    const char* qs[] = {"1", "x", "1 + x^2", "-x", "3"};
    for (size_t i = 0; i < 5; i++)
        for (size_t j = 0; j < 5; j++)
            frob_ordinary_ok(ps[i], qs[j]);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(t_kov_stress_poly);
    TEST(t_kov_stress_apparent);
    TEST(t_kov_stress_pole);
    TEST(t_kov_stress_euler);
    TEST(t_frob_stress_ordinary);

    printf("\nAll DSolve M5 stress tests passed.\n");
    return 0;
}
