/* test_algebraicnumberpolynomial.c — AlgebraicNumberPolynomial[a, x].
 *
 * Covers the number passthrough (integer / bigint / rational), the polynomial
 * build from an AlgebraicNumber object (basic, rational coefficients, dropped
 * zero coefficient, higher degree), Listable threading over both arguments, the
 * defining round-trip (poly /. x -> theta == a, checked without a numeric oracle
 * via the RootReduce zero test), the addition-of-algebraic-numbers-as-polynomials
 * identity, and the decline paths (the AlgebraicNumberPolynomial::naobj message
 * for a non-AlgebraicNumber argument, and wrong argument count).
 *
 * The AlgebraicNumber objects it consumes require FLINT (the qqbar engine) to be
 * built, so the whole suite SKIPs cleanly when FLINT is off.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "flint_bridge.h"

#include "test_utils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static Expr* eval_str(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* e = evaluate(parsed);
    expr_free(parsed);
    return e;
}

/* Assert the short-form printed value of `input` equals `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = eval_str(input);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* Value-preservation without a numeric oracle: assert RootReduce[expr] is "0". */
static void check_zero(const char* expr) {
    char buf[1024];
    snprintf(buf, sizeof buf, "RootReduce[%s]", expr);
    Expr* e = eval_str(buf);
    char* s = expr_to_string(e);
    if (strcmp(s, "0") != 0) {
        fprintf(stderr, "FAIL(nonzero): RootReduce[%s] = %s (expected 0)\n", expr, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* Capture stderr while `input` is evaluated; return the collected text (heap,
 * caller frees) and write the printed result into *out (also heap). */
static char* eval_capturing_stderr(const char* input, char** out) {
    const char* path = "/tmp/mathilda_anp_stderr.log";
    fflush(stderr);
    if (!freopen(path, "w+", stderr)) { if (out) *out = NULL; return NULL; }
    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out) *out = expr_to_string(e);
    expr_free(p);
    expr_free(e);
    fflush(stderr);
    freopen("/dev/tty", "w", stderr);

    FILE* f = fopen(path, "r");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    remove(path);
    return buf;
}

/* ---------------- Number passthrough --------------------------------------- */

static void test_number_passthrough(void) {
    check("AlgebraicNumberPolynomial[2, x]", "2");
    check("AlgebraicNumberPolynomial[-3, x]", "-3");
    check("AlgebraicNumberPolynomial[1/2, x]", "1/2");
    check("AlgebraicNumberPolynomial[100!, x] === 100!", "True");
    /* An integer input ignores the variable entirely. */
    check("AlgebraicNumberPolynomial[0, x]", "0");
}

/* ---------------- Polynomial from an AlgebraicNumber object ----------------- */

static void test_polynomial_build(void) {
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x]",
          "1 + 2 x");
    /* Rational coefficient. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 1/2}], x]",
          "1 + 1/2 x");
    /* Degree-3 generator, zero middle coefficient is dropped. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[2^(1/3), {1, 0, 3}], x]",
          "1 + 3 x^2");
    /* Degree-4 generator (algebraic integer), coefficients preserved. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,2,3,4}], x]",
          "1 + 2 x + 3 x^2 + 4 x^3");
    /* A variable other than x is honoured. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], y]",
          "1 + 2 y");
    /* An over-length list folds inside AlgebraicNumber first: {1,0,3} over
     * Sqrt[2] collapses to the rational 7, so the polynomial is that constant. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 0, 3}], x]", "7");
}

/* ---------------- Listable threading --------------------------------------- */

static void test_listable(void) {
    check("AlgebraicNumberPolynomial[{2, AlgebraicNumber[Sqrt[2], {1, 2}]}, x]",
          "{2, 1 + 2 x}");
    /* Threads over the variable argument too. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], {x, y}]",
          "{1 + 2 x, 1 + 2 y}");
    check("Attributes[AlgebraicNumberPolynomial]", "{Listable, Protected}");
}

/* ---------------- Defining round-trip (poly /. x -> theta == a) ------------- */

static void test_round_trip(void) {
    /* AlgebraicNumber is by definition a polynomial in its generator. */
    check_zero("(AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2+Sqrt[3]], {1,2,3,4}], x] "
               "/. x -> Sqrt[2+Sqrt[3]]) - AlgebraicNumber[Sqrt[2+Sqrt[3]], {1,2,3,4}]");
    check_zero("(AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 1/2}], x] "
               "/. x -> Sqrt[2]) - AlgebraicNumber[Sqrt[2], {1, 1/2}]");
    check_zero("(AlgebraicNumberPolynomial[AlgebraicNumber[Root[#^3+#+1&,3], {1,2,1}], x] "
               "/. x -> Root[#^3+#+1&,3]) - AlgebraicNumber[Root[#^3+#+1&,3], {1,2,1}]");
}

/* ---------------- Addition of algebraic numbers via polynomials ------------- */

static void test_addition_via_polynomials(void) {
    /* Same generator: adding the polynomials matches adding the numbers. */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,2,3,4}], x] + "
          "AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,-2,1,1}], x]",
          "2 + 4 x^2 + 5 x^3");
    check_zero("((AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,2,3,4}], x] + "
               "AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,-2,1,1}], x]) "
               "/. x -> Sqrt[2]+Sqrt[3]) - "
               "(AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,2,3,4}] + "
               "AlgebraicNumber[Sqrt[2]+Sqrt[3], {1,-2,1,1}])");
}

/* ---------------- Declines and the naobj message --------------------------- */

static void test_declines_and_message(void) {
    /* A non-AlgebraicNumber, non-rational argument: message + unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("AlgebraicNumberPolynomial[Sqrt[2], x]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberPolynomial::naobj") != NULL,
               "expected AlgebraicNumberPolynomial::naobj, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberPolynomial[Sqrt[2], x]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* Symbolic argument also declines with a message. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberPolynomial[y^2, z]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberPolynomial::naobj") != NULL,
               "expected AlgebraicNumberPolynomial::naobj, got: %s", err);
    free(err);
    free(result);

    /* Wrong argument count: silently unevaluated (no message). */
    check("AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}]]",
          "AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}]]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping "
               "AlgebraicNumberPolynomial tests.\n");
        return 0;
    }
    test_number_passthrough();
    test_polynomial_build();
    test_listable();
    test_round_trip();
    test_addition_via_polynomials();
    test_declines_and_message();
    printf("test_algebraicnumberpolynomial: all passed\n");
    return 0;
}
