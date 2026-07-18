#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate `input` and assert its FullForm equals `expected`. */
static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "UnitVector %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* Evaluate a boolean expression and assert it reduces to True. */
static void run_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, "True") == 0, "expected True for %s, got %s", input, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- One-argument (2-D) form ---------- */

static void test_uv_1d_k1(void) { run_full("UnitVector[1]", "List[1, 0]"); }
static void test_uv_1d_k2(void) { run_full("UnitVector[2]", "List[0, 1]"); }

/* ---------- Two-argument (n-D) form ---------- */

static void test_uv_3d_k1(void) { run_full("UnitVector[3, 1]", "List[1, 0, 0]"); }
static void test_uv_3d_k2(void) { run_full("UnitVector[3, 2]", "List[0, 1, 0]"); }
static void test_uv_3d_k3(void) { run_full("UnitVector[3, 3]", "List[0, 0, 1]"); }
static void test_uv_1x1(void)   { run_full("UnitVector[1, 1]", "List[1]"); }

/* ---------- Large vector + structural invariants ---------- */

static void test_uv_large_length(void) {
    run_true("Length[UnitVector[100, 4]] == 100");
}
static void test_uv_large_position(void) {
    run_true("Position[UnitVector[100, 4], 1] == {{4}}");
}
static void test_uv_large_zero_count(void) {
    run_true("100 - Count[UnitVector[100, 4], 0] == 1");
}
static void test_uv_invariants_arbitrary(void) {
    /* The request's invariants, on a fixed (n, k). */
    run_true("Length[UnitVector[37, 20]] == 37");
    run_true("Position[UnitVector[37, 20], 1] == {{20}}");
    run_true("37 - Count[UnitVector[37, 20], 0] == 1");
}

/* ---------- Numeric components ---------- */

static void test_uv_n_threads(void) {
    /* N threads over the exact list; components become machine reals. */
    run_full("N[UnitVector[2]]", "List[0.0, 1.0]");
}
static void test_uv_wp_machine(void) {
    run_full("UnitVector[2, WorkingPrecision -> MachinePrecision]",
             "List[0.0, 1.0]");
    run_true("Precision[UnitVector[2, WorkingPrecision -> MachinePrecision][[2]]] "
             "== MachinePrecision");
}
static void test_uv_wp_infinity_exact(void) {
    run_full("UnitVector[3, 2, WorkingPrecision -> Infinity]", "List[0, 1, 0]");
}
static void test_uv_wp_mpfr(void) {
    /* A digit count above machine precision yields high-precision reals. */
    run_true("Precision[UnitVector[3, 2, WorkingPrecision -> 30][[2]]] > 20");
    run_true("Precision[UnitVector[3, 2, WorkingPrecision -> 30][[1]]] > 20");
}

/* ---------- Errors and unevaluated cases ---------- */

static void test_uv_no_args(void) {
    /* UnitVector::argt printed to stderr; call stays unevaluated. */
    run_full("UnitVector[]", "UnitVector[]");
}
static void test_uv_too_many(void) {
    /* UnitVector::nonopt printed to stderr; call stays unevaluated. */
    run_full("UnitVector[1, 2, 3]", "UnitVector[1, 2, 3]");
}
static void test_uv_out_of_range_high(void) {
    run_full("UnitVector[2, 3]", "UnitVector[2, 3]");
}
static void test_uv_out_of_range_1arg(void) {
    /* UnitVector[3] is UnitVector[2, 3]: k=3 out of range in 2 dims. */
    run_full("UnitVector[3]", "UnitVector[3]");
}
static void test_uv_zero_k(void) {
    run_full("UnitVector[3, 0]", "UnitVector[3, 0]");
}
static void test_uv_symbolic(void) {
    run_full("UnitVector[x]", "UnitVector[x]");
    run_full("UnitVector[n, k]", "UnitVector[n, k]");
}

int main(void) {
    symtab_init();
    core_init();

    test_uv_1d_k1();
    test_uv_1d_k2();
    test_uv_3d_k1();
    test_uv_3d_k2();
    test_uv_3d_k3();
    test_uv_1x1();
    test_uv_large_length();
    test_uv_large_position();
    test_uv_large_zero_count();
    test_uv_invariants_arbitrary();
    test_uv_n_threads();
    test_uv_wp_machine();
    test_uv_wp_infinity_exact();
    test_uv_wp_mpfr();
    test_uv_no_args();
    test_uv_too_many();
    test_uv_out_of_range_high();
    test_uv_out_of_range_1arg();
    test_uv_zero_k();
    test_uv_symbolic();

    symtab_clear();
    printf("All UnitVector tests passed.\n");
    return 0;
}
