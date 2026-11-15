/*
 * test_solve_exponential.c
 *
 * Exact exponential equations in Solve / Reduce (src/solve/solve_explog.c and
 * its glue in src/solve/solveinv.c).
 *
 * The motivating bug: Solve[4^100 + 4^100 == 8^x, x] answered
 * Log[2^201]/Log[8] + 2 Pi I C[1]/Log[8] instead of 67 + 2 Pi I C[1]/Log[8].
 * Over Reals it returned the same complex family, and over Integers it
 * returned {}, a wrong answer. Reduce over Integers inherited that as False.
 *
 * Every expected answer here was derived by hand and checked by
 * back-substitution with PossibleZeroQ (see test_backsubstitution). An exact
 * InputForm string pins the canonical output shape.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "solve/solve_explog.h"
#include "test_utils.h"

/* Evaluate `input` (wrapped in ToString[InputForm[...]] and Quiet) and
 * compare it with `expected`. */
static void run_test(const char* input, const char* expected) {
    size_t n = strlen(input) + 64;
    char* buf = (char*)malloc(n);
    snprintf(buf, n, "Quiet[ToString[InputForm[%s]]]", input);
    Expr* e = parse_expression(buf);
    free(buf);
    if (!e) { printf("FAIL: failed to parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    const char* got = res->type == EXPR_STRING ? res->data.string : "<not a string>";
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        expr_free(res); expr_free(e); ASSERT(0); return;
    }
    printf("PASS: %s\n", input);
    expr_free(res); expr_free(e);
}

/* ------------------------------------------------------------------ *
 *  The multiplicative basis (C level).                                *
 * ------------------------------------------------------------------ */

static void basis_check(const char* const* vals, size_t nv, const char* expect) {
    mpz_t* v = (mpz_t*)malloc(sizeof(mpz_t) * nv);
    for (size_t i = 0; i < nv; i++) mpz_init_set_str(v[i], vals[i], 10);
    ExplogBasis B;
    explog_basis_build((const mpz_t*)v, nv, &B);
    char got[512] = "";
    for (size_t k = 0; k < B.n; k++) {
        char* s = mpz_get_str(NULL, 10, B.p[k]);
        if (k) strcat(got, ",");
        strncat(got, s, sizeof(got) - strlen(got) - 1);
        free(s);
    }
    /* Every input is an exact product of basis powers. */
    int64_t* e = (int64_t*)calloc(B.n ? B.n : 1, sizeof(int64_t));
    for (size_t i = 0; i < nv; i++) ASSERT(explog_basis_exponents(&B, v[i], e));
    free(e);
    if (strcmp(got, expect) != 0) {
        printf("FAIL: basis expected {%s} got {%s}\n", expect, got);
        ASSERT(0);
    } else {
        printf("PASS: basis {%s}\n", got);
    }
    explog_basis_clear(&B);
    for (size_t i = 0; i < nv; i++) mpz_clear(v[i]);
    free(v);
}

static void test_basis(void) {
    { const char* v[] = { "8", "3213876088517980551083924184682325205044405987565585670602752" };
      basis_check(v, 2, "2"); }                                  /* 8, 2^201 */
    { const char* v[] = { "4", "8" };        basis_check(v, 2, "2"); }
    { const char* v[] = { "12", "18" };      basis_check(v, 2, "2,3"); }
    { const char* v[] = { "6", "10", "15" }; basis_check(v, 3, "2,3,5"); }
    { const char* v[] = { "36", "216" };     basis_check(v, 2, "6"); }
    { const char* v[] = { "1", "1" };        basis_check(v, 2, ""); }
    { const char* v[] = { "1000", "100" };   basis_check(v, 2, "10"); }

    /* Exponents, including a 10^6-bit power (no factoring: mpz_remove). */
    ExplogBasis B;
    mpz_t v[2];
    mpz_init_set_ui(v[0], 16);
    mpz_init(v[1]);
    mpz_ui_pow_ui(v[1], 2, 1000000);
    explog_basis_build((const mpz_t*)v, 2, &B);
    ASSERT(B.n == 1 && mpz_cmp_ui(B.p[0], 2) == 0);
    int64_t e = 0;
    ASSERT(explog_basis_exponents(&B, v[1], &e) && e == 1000000);
    mpz_set_ui(v[0], 3);
    ASSERT(!explog_basis_exponents(&B, v[0], &e));             /* 3 not in it */
    explog_basis_clear(&B);
    mpz_clears(v[0], v[1], NULL);
    printf("PASS: basis exponents\n");
}

/* ------------------------------------------------------------------ *
 *  The reported case, every domain.                                   *
 * ------------------------------------------------------------------ */

static void test_reported(void) {
    run_test("Solve[4^100 + 4^100 == 8^x, x]",
        "{{x -> ConditionalExpression[67 + ((2*I)*Pi*C[1])/Log[8], Element[C[1], Integers]]}}");
    run_test("Solve[4^100 + 4^100 == 8^x, x, Reals]", "{{x -> 67}}");
    run_test("Solve[4^100 + 4^100 == 8^x, x, Integers]", "{{x -> 67}}");
    run_test("Solve[4^100 + 4^100 == 8^x, x, Rationals]", "{{x -> 67}}");
    run_test("Reduce[4^100 + 4^100 == 8^x, x]",
        "Element[C[1], Integers] && x == 67 + ((2*I)*Pi*C[1])/Log[8]");
    run_test("Reduce[4^100 + 4^100 == 8^x, x, Reals]", "x == 67");
    run_test("Reduce[4^100 + 4^100 == 8^x, x, Integers]", "x == 67");
}

/* ------------------------------------------------------------------ *
 *  One exponential.                                                   *
 * ------------------------------------------------------------------ */

static void test_single_base(void) {
    run_test("Solve[2^x == 8, x]",
        "{{x -> ConditionalExpression[3 + ((2*I)*Pi*C[1])/Log[2], Element[C[1], Integers]]}}");
    run_test("Solve[2^x == 8, x, Reals]", "{{x -> 3}}");
    run_test("Solve[3^x == 81, x, Integers]", "{{x -> 4}}");
    run_test("Solve[4^x == 8, x]",
        "{{x -> ConditionalExpression[3/2 + ((2*I)*Pi*C[1])/Log[4], Element[C[1], Integers]]}}");
    run_test("Solve[4^x == 8, x, Integers]", "{}");
    run_test("Solve[4^x == 8, x, Rationals]", "{{x -> 3/2}}");
    run_test("Solve[(1/8)^x == 32, x, Reals]", "{{x -> -5/3}}");
    run_test("Solve[(2/3)^x == 8/27, x, Reals]", "{{x -> 3}}");
    run_test("Solve[Sqrt[2]^x == 8, x, Reals]", "{{x -> 6}}");
    run_test("Solve[E^x == 5, x, Reals]", "{{x -> Log[5]}}");
    run_test("Solve[2^x == -8, x, Reals]", "{}");
    run_test("Solve[2^x == 0, x]", "{}");
    run_test("Solve[2^x == 3, x, Integers]", "{}");
    run_test("Solve[2^x == 3, x, Reals]", "{{x -> Log[3]/Log[2]}}");
    run_test("Solve[2^x == a, x, Reals]",
        "{{x -> ConditionalExpression[Log[a]/Log[2], a > 0]}}");
    run_test("Solve[16^x == 2^(10^6), x, Integers]", "{{x -> 250000}}");
    run_test("Solve[8^x == 2^201, x, Rationals]", "{{x -> 67}}");
    run_test("Solve[Pi^x == Pi^2, x, Reals]", "{{x -> 2}}");
    run_test("Solve[2^x == Pi, x, Integers]", "{}");
    run_test("Solve[2^x == 3 - 2 Sqrt[3], x, Reals]", "{}");
    run_test("Solve[2^(a x) == 8, x, Reals]", "{{x -> 3/a}}");
    run_test("Solve[2^(1/x) == 8, x, Reals]", "{{x -> 1/3}}");
}

/* ------------------------------------------------------------------ *
 *  Several exponentials (multiplicative basis of rank 1 or more).     *
 * ------------------------------------------------------------------ */

static void test_several(void) {
    run_test("Solve[4^x == 8^(x + 1), x]",
        "{{x -> ConditionalExpression[-3 + ((2*I)*Pi*C[1])/Log[2], Element[C[1], Integers]]}}");
    run_test("Solve[9^x == 27^(x - 1), x, Reals]", "{{x -> 3}}");
    run_test("Solve[2^(3 x) == 4^(x + 5), x, Integers]", "{{x -> 10}}");
    run_test("Solve[1000^x == 10^(x + 4), x]",
        "{{x -> ConditionalExpression[2 + ((2*I)*Pi*C[1])/Log[100], Element[C[1], Integers]]}}");
    run_test("Solve[2^x == 3^x, x, Reals]", "{{x -> 0}}");
    run_test("Solve[2^x == 3^x, x]",
        "{{x -> ConditionalExpression[((2*I)*Pi*C[1])/Log[3/2], Element[C[1], Integers]]}}");
    run_test("Solve[12^x == 18^x, x, Reals]", "{{x -> 0}}");
    run_test("Solve[2^x 3^x == 36, x, Reals]", "{{x -> 2}}");
    run_test("Solve[2^x 5^x == 1000, x, Integers]", "{{x -> 3}}");
    run_test("Solve[3 2^x == 5^x, x, Reals]", "{{x -> Log[3]/Log[5/2]}}");
    run_test("Solve[(2^60 3^40)^x == 2^120 3^80, x, Reals]", "{{x -> 2}}");
    run_test("Solve[2^x == 2 2^x, x]", "{}");
    run_test("Solve[2^x == -2^x, x]", "{}");
    /* A negative base: the real root sits at C[1] = 1, not at C[1] = 0. */
    run_test("Solve[(-2)^x == 4, x, Reals]", "{{x -> 2}}");
}

/* ------------------------------------------------------------------ *
 *  Nonlinear exponents and polynomials in one exponential.            *
 * ------------------------------------------------------------------ */

static void test_nonlinear_and_kernels(void) {
    run_test("Solve[2^(x^2) == 16, x, Reals]", "{{x -> -2}, {x -> 2}}");
    run_test("Solve[2^(x^2) == 16, x, Integers]", "{{x -> -2}, {x -> 2}}");
    run_test("Solve[2^(x^2 - 1) == 8, x, Reals]", "{{x -> -2}, {x -> 2}}");
    run_test("Solve[4^x - 3 2^x + 2 == 0, x, Reals]", "{{x -> 0}, {x -> 1}}");
    run_test("Solve[4^x - 3 2^x + 2 == 0, x]",
        "{{x -> ConditionalExpression[((2*I)*Pi*C[1])/Log[2], Element[C[1], Integers]]}, "
        "{x -> ConditionalExpression[1 + ((2*I)*Pi*C[1])/Log[2], Element[C[1], Integers]]}}");
    run_test("Solve[9^x - 3^(x + 1) - 54 == 0, x, Reals]", "{{x -> 2}}");
    run_test("Solve[2^x + 2^(x + 1) == 24, x, Reals]", "{{x -> 3}}");
    run_test("Solve[(3^x)^2 == 81, x, Reals]", "{{x -> 2}}");
    run_test("Solve[E^(2 x) - 3 E^x + 2 == 0, x, Reals]", "{{x -> 0}, {x -> Log[2]}}");
}

/* ------------------------------------------------------------------ *
 *  Inequality constraints make the problem real.                      *
 * ------------------------------------------------------------------ */

static void test_constrained(void) {
    run_test("Solve[2^x == 8 && x > 0, x]", "{{x -> 3}}");
    run_test("Solve[2^x == 8 && x > 5, x]", "{}");
    run_test("Solve[2^x == 3 && x > 0, x]", "{{x -> Log[3]/Log[2]}}");
    run_test("Solve[4^x - 3 2^x + 2 == 0 && x > 1/2, x]", "{{x -> 1}}");
    run_test("Reduce[2^x == 8 && x > 0, x]", "x == 3");
}

/* ------------------------------------------------------------------ *
 *  Reduce inherits every domain through Solve; inexact input.         *
 * ------------------------------------------------------------------ */

static void test_reduce_and_inexact(void) {
    /* was False: the real root x = 2 sits at C[1] = 1 */
    run_test("Reduce[(-2)^x == 4, x, Reals]", "x == 2");
    run_test("Reduce[2^x == 3^x, x, Reals]", "x == 0");
    run_test("Reduce[2^x == a, x, Reals]", "a > 0 && x == Log[a]/Log[2]");
    run_test("Reduce[4^x - 3 2^x + 2 == 0, x, Reals]", "x == 0 || x == 1");
    run_test("Reduce[3^x == 81, x, Integers]", "x == 4");
    /* the generated parameter stays an exact name under numericalisation */
    run_test("Solve[Exp[x] == 2.0, x]",
        "{{x -> ConditionalExpression[0.6931471805599453 + (0.0 + 6.283185307179586*I)*C[1], "
        "Element[C[1], Integers]]}}");
    run_test("Solve[2^x == 8.0, x, Reals]", "{{x -> 3.0}}");
}

/* ------------------------------------------------------------------ *
 *  Back-substitution: every answer above solves its equation.         *
 * ------------------------------------------------------------------ */

static void test_backsubstitution(void) {
    /* For each row: strip the ConditionalExpression and check the residual at
     * C[1] = -1, 0, 2 with PossibleZeroQ. */
    Expr* def = parse_expression(
        "chk[eq_, x_, dom___] := Module[{s = Solve[eq, x, dom], r = eq[[1]] - eq[[2]]},"
        "  ListQ[s] && Length[s] > 0 && And @@ Flatten[Table["
        "    PossibleZeroQ[r /. (row /. ConditionalExpression[v_, _] :> v) /. C[1] -> k],"
        "    {row, s}, {k, {-1, 0, 2}}]]]");
    expr_free(evaluate(def));
    expr_free(def);
    const char* eqs[] = {
        "chk[4^100 + 4^100 == 8^x, x]", "chk[4^x == 8^(x + 1), x]",
        "chk[2^x == 3^x, x]", "chk[1000^x == 10^(x + 4), x]",
        "chk[3 2^x == 5^x, x, Reals]", "chk[2^(x^2) == 16, x]",
        "chk[4^x - 3 2^x + 2 == 0, x]", "chk[9^x - 3^(x + 1) - 54 == 0, x]",
        "chk[(-2)^x == 4, x, Reals]", "chk[2^x 3^x == 36, x]",
        "chk[Pi^x == Pi^2, x]", "chk[(2^60 3^40)^x == 2^120 3^80, x]",
    };
    for (size_t i = 0; i < sizeof(eqs) / sizeof(eqs[0]); i++) run_test(eqs[i], "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basis);
    TEST(test_reported);
    TEST(test_single_base);
    TEST(test_several);
    TEST(test_nonlinear_and_kernels);
    TEST(test_constrained);
    TEST(test_reduce_and_inexact);
    TEST(test_backsubstitution);

    printf("All solve-exponential tests passed!\n");
    return 0;
}
