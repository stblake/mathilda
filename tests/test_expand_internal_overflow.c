/* test_expand_internal_overflow.c -- Overflow[] is a user-level answer only.
 *
 * The Expand builtin answers Overflow[] when the estimated expansion passes the
 * memory ceiling (EXPAND_OVERFLOW_CAP).  For the user that is the documented
 * answer.  A call made from inside another builtin needs a value instead:
 * Overflow[] entered the caller's arithmetic, Overflow[]/Overflow[] is
 * Indeterminate, and Simplify's trig-rational pass (src/simp/trigrat.c), whose
 * denominator rationalisation multiplied a 530k-leaf denominator by its
 * conjugate, returned Indeterminate as "the simplest form" of a perfectly
 * regular algebraic function.  An internal Expand that would overflow now
 * returns its subject unexpanded (the same value).
 *
 *   1. user level: Expand[(a + ... + j)^30] is still Overflow[];
 *   2. inside a builtin (builtin depth > 1): the same call returns its subject;
 *   3. the reported case: Simplify of D[F, x], F the GeneralisedGoursat
 *      antiderivative of a cube-root integrand over Q(2^(1/3), 3^(1/6),
 *      37^(1/3), ...), stored in data/radical_field/
 *      goursat_cbrt_overflow_derivative.m (argv[1]), is no longer
 *      Indeterminate: it equals the integrand at two complex points and is no
 *      larger.
 *
 * Uses ASSERT_MSG (never libc assert: Release defines NDEBUG). */
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

static char data_path[1200] = "data/radical_field/goursat_cbrt_overflow_derivative.m";

static char* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_inputform_display(res);
    expr_free(res);
    expr_free(e);
    return s;
}

static void run_test(const char* input, const char* expected) {
    char* s = eval_str(input);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %.300s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %.300s", input, expected, s);
    free(s);
}

static void test_user_level_overflow(void) {
    /* C(39, 9) ~ 2.1e8 terms: past the ceiling, answered by the estimate alone */
    run_test("Expand[(a + b + c + d + e + f + g + h + i + j)^30]", "Overflow[]");
}

static void test_internal_level_returns_subject(void) {
    /* as if called from inside another builtin */
    int* depth = eval_builtin_depth_counter();
    *depth += 2;
    Expr* e = parse_expression("Expand[(a + b + c + d + e + f + g + h + i + j)^30]");
    Expr* r = evaluate(e);
    *depth -= 2;
    char* s = expr_to_string_inputform_display(r);
    const char* want = "(a + b + c + d + e + f + g + h + i + j)^30";
    if (strcmp(s, want) != 0) printf("FAIL: internal Expand\n  expected: %s\n  got:      %.300s\n", want, s);
    ASSERT_MSG(strcmp(s, want) == 0, "internal Expand: expected %s, got %.300s", want, s);
    free(s);
    expr_free(r); expr_free(e);
    /* an internal Expand within the ceiling still expands */
    *depth += 2;
    e = parse_expression("Expand[(a + b)^2]");
    r = evaluate(e);
    *depth -= 2;
    s = expr_to_string_inputform_display(r);
    ASSERT_MSG(strcmp(s, "a^2 + 2*a*b + b^2") == 0, "internal Expand[(a + b)^2]: got %s", s);
    free(s);
    expr_free(r); expr_free(e);
}

static void test_reported_case(void) {
    char buf[1600];
    /* f is the integrand of benchmark 265's G-s7-cbrt-symmetric-001 */
    snprintf(buf, sizeof buf,
             "ovd = Get[\"%s\"]; "
             "ovf = (-1512 + 6545*x - 8346*x^2 + 3619*x^3 - 2170*x^4 + 369*x^5 - 128*x^6)/"
             "((-2 + 2*x - 2*x^2)^(1/3)*(216 - 580*x + 324*x^2 + 108*x^3 - 548*x^4 + "
             "324*x^5 - 108*x^6 + 32*x^7)); LeafCount[ovd]", data_path);
    run_test(buf, "1152");
    run_test("ovs = Simplify[ovd]; FreeQ[ovs, Indeterminate | Overflow]", "True");
    run_test("Max[Abs[N[(ovs - ovf) /. x -> #, 30] & /@ {N[3/7 + 2 I/5, 40], N[13/5 - I/7, 40]}]] < 10^-20",
             "True");
    run_test("LeafCount[ovs] <= LeafCount[ovf]", "True");
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);
    printf("Starting expand_internal_overflow_tests\n");
    if (argc > 1) snprintf(data_path, sizeof data_path, "%s", argv[1]);
    symtab_init();
    core_init();
    TEST(test_user_level_overflow);
    TEST(test_internal_level_returns_subject);
    TEST(test_reported_case);
    printf("All expand_internal_overflow tests passed!\n");
    return 0;
}
