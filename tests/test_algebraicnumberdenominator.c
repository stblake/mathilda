/* test_algebraicnumberdenominator.c — AlgebraicNumberDenominator[a].
 *
 * AlgebraicNumberDenominator[a] is the smallest positive integer n with n a an
 * algebraic integer.  Covers rationals/integers, radicals, complex algebraic
 * numbers, Root objects, AlgebraicNumber objects (including the cases where the
 * answer is a proper divisor of the min-poly leading coefficient — the reason the
 * implementation does a per-prime valuation rather than reading that coefficient),
 * Listable threading, the defining integrality invariant
 * (AlgebraicIntegerQ[AlgebraicNumberDenominator[a] a] == True), the decline path
 * (the AlgebraicNumberDenominator::nalg message on a non-algebraic argument), and
 * wrong arity.
 *
 * The algebraic numbers it consumes require FLINT (the qqbar engine), so the whole
 * suite SKIPs cleanly when FLINT is off.
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

/* Capture stderr while `input` is evaluated; return the collected text (heap,
 * caller frees) and write the printed result into *out (also heap). */
static char* eval_capturing_stderr(const char* input, char** out) {
    const char* path = "/tmp/mathilda_andenom_stderr.log";
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

/* ---------------- Rationals and integers ----------------------------------- */

static void test_rationals(void) {
    check("AlgebraicNumberDenominator[1/3]", "3");
    check("AlgebraicNumberDenominator[3]", "1");
    check("AlgebraicNumberDenominator[1/2]", "2");
    check("AlgebraicNumberDenominator[-1/6]", "6");
    check("AlgebraicNumberDenominator[0]", "1");
    check("AlgebraicNumberDenominator[5]", "1");
    check("AlgebraicNumberDenominator[7/10]", "10");
    /* Denominator of a rational is its (reduced) denominator. */
    check("AlgebraicNumberDenominator[6/4]", "2");
}

/* ---------------- Real radicals -------------------------------------------- */

static void test_radicals(void) {
    check("AlgebraicNumberDenominator[1/Sqrt[3]]", "3");
    check("AlgebraicNumberDenominator[(1 + Sqrt[5])/2]", "1");   /* golden ratio: an algebraic integer */
    check("AlgebraicNumberDenominator[Sqrt[2]]", "1");
    check("AlgebraicNumberDenominator[1/Sqrt[2]]", "2");
    check("AlgebraicNumberDenominator[1/Sqrt[Sqrt[2] + 3]]", "7");
    check("AlgebraicNumberDenominator[1/(1 + Sqrt[3])]", "2");
    check("AlgebraicNumberDenominator[Sqrt[2] + Sqrt[3]]", "1");
    check("AlgebraicNumberDenominator[2^(1/3)/2]", "2");
}

/* ---------------- Complex algebraic numbers -------------------------------- */

static void test_complex(void) {
    check("AlgebraicNumberDenominator[(1 + 3 I)^(-1/3)]", "10");
    check("AlgebraicNumberDenominator[I]", "1");
    check("AlgebraicNumberDenominator[1/(1 + I)]", "2");
    check("AlgebraicNumberDenominator[1/Sqrt[1 + I]]", "2");
}

/* ---------------- Root objects --------------------------------------------- */

static void test_root_objects(void) {
    check("AlgebraicNumberDenominator[Root[5 - 6 #1 + 3 #1^3 &, 1]]", "3");
    check("AlgebraicNumberDenominator[Root[#1^2 - 2 &, 1]]", "1");  /* +/- Sqrt[2] */
}

/* ---------------- AlgebraicNumber objects ---------------------------------- */
/* The divergence cases: the answer is a PROPER divisor of the primitive minimal
 * polynomial's leading coefficient, so a naive qqbar_denominator would be wrong. */

static void test_algebraicnumber_objects(void) {
    /* 1/5 + Sqrt[2]: primitive min-poly 25 x^2 - 10 x - 49 (leading 25), denom 5. */
    check("AlgebraicNumberDenominator[AlgebraicNumber[Sqrt[2], {1/5, 1}]]", "5");
    /* (1 + Sqrt[2])/3: primitive min-poly 9 x^2 - 6 x - 1 (leading 9), denom 3. */
    check("AlgebraicNumberDenominator[AlgebraicNumber[Sqrt[2], {1/3, 1/3}]]", "3");
    /* Pure Sqrt[2] (an algebraic integer). */
    check("AlgebraicNumberDenominator[AlgebraicNumber[Sqrt[2], {0, 1}]]", "1");
    /* Also accepts the object built via ToNumberField. */
    check("AlgebraicNumberDenominator[ToNumberField[1/(1 + Sqrt[3]), Sqrt[3]]]", "2");
}

/* ---------------- Listable threading --------------------------------------- */

static void test_listable(void) {
    check("AlgebraicNumberDenominator[{Sqrt[2], 1/Sqrt[2], 1/3}]", "{1, 2, 3}");
    check("AlgebraicNumberDenominator[{1/2, 1/3}]", "{2, 3}");
    check("Attributes[AlgebraicNumberDenominator]", "{Listable, Protected}");
}

/* ---------------- Integrality invariant ------------------------------------ */
/* By definition n a is an algebraic integer; and where a itself is not. */

static void test_integrality_invariant(void) {
    check("AlgebraicIntegerQ[1/Sqrt[1 + I]]", "False");
    check("AlgebraicIntegerQ[AlgebraicNumberDenominator[1/Sqrt[1 + I]] / Sqrt[1 + I]]", "True");
    check("AlgebraicIntegerQ[AlgebraicNumberDenominator[1/Sqrt[3]] / Sqrt[3]]", "True");
    check("AlgebraicIntegerQ[AlgebraicNumberDenominator[(1 + 3 I)^(-1/3)] (1 + 3 I)^(-1/3)]", "True");
    check("AlgebraicIntegerQ[AlgebraicNumberDenominator[AlgebraicNumber[Sqrt[2], {1/5, 1}]] "
          "AlgebraicNumber[Sqrt[2], {1/5, 1}]]", "True");
    /* n is minimal: (n-1) a for the 1/Sqrt[3] case is NOT an algebraic integer. */
    check("AlgebraicIntegerQ[2 / Sqrt[3]]", "False");
}

/* ---------------- Declines and the nalg message ---------------------------- */

static void test_declines_and_message(void) {
    /* A non-algebraic argument: message + unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("AlgebraicNumberDenominator[Pi]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberDenominator::nalg") != NULL,
               "expected AlgebraicNumberDenominator::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberDenominator[Pi]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* A free symbol also declines with a message. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberDenominator[x]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberDenominator::nalg") != NULL,
               "expected AlgebraicNumberDenominator::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberDenominator[x]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* Log[2] is not an explicit algebraic number either. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberDenominator[Log[2]]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberDenominator::nalg") != NULL,
               "expected AlgebraicNumberDenominator::nalg, got: %s", err);
    free(err);
    free(result);

    /* Wrong argument count: silently unevaluated (no message). */
    check("AlgebraicNumberDenominator[1, 2]", "AlgebraicNumberDenominator[1, 2]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping "
               "AlgebraicNumberDenominator tests.\n");
        return 0;
    }
    test_rationals();
    test_radicals();
    test_complex();
    test_root_objects();
    test_algebraicnumber_objects();
    test_listable();
    test_integrality_invariant();
    test_declines_and_message();
    printf("test_algebraicnumberdenominator: all passed\n");
    return 0;
}
