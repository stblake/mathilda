/* test_algebraicnumbertrace.c — AlgebraicNumberTrace[a].
 *
 * AlgebraicNumberTrace[a] is the field trace of the algebraic number a: the sum
 * of the roots of a's minimal polynomial, i.e. -(coeff of x^{deg-1}) / (leading
 * coeff).  AlgebraicNumberTrace[a, Extension -> theta] is the relative trace over
 * Q(theta), equal to the absolute trace scaled by the tower index [Q(theta):Q(a)]
 * (defined only when a lies in Q(theta)).
 *
 * Covers integers/rationals, real and complex radicals, roots of unity, Root and
 * AlgebraicNumber objects, GoldenRatio (the sole algebraic named constant, handled
 * by the shared qqbar converter), Listable threading with a trailing Extension
 * option, the relative-trace cases (including additivity over a fixed field and a
 * compositum built by ToNumberField), and the decline paths
 * (AlgebraicNumberTrace::nalg on a non-algebraic argument, AlgebraicNumberTrace::ext
 * on an argument not in the extension, silent unevaluation on wrong arity).
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
    const char* path = "/tmp/mathilda_antrace_stderr.log";
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

/* ---------------- Integers and rationals ----------------------------------- */

static void test_rationals(void) {
    check("AlgebraicNumberTrace[2]", "2");
    check("AlgebraicNumberTrace[-2/3]", "-2/3");
    check("AlgebraicNumberTrace[0]", "0");
    check("AlgebraicNumberTrace[1]", "1");
    check("AlgebraicNumberTrace[-5]", "-5");
    check("AlgebraicNumberTrace[7/10]", "7/10");
    /* Trace of a rational is the rational itself (minimal polynomial q x - p). */
    check("AlgebraicNumberTrace[6/4]", "3/2");
}

/* ---------------- Real radicals and GoldenRatio ---------------------------- */

static void test_radicals(void) {
    check("AlgebraicNumberTrace[Sqrt[2]]", "0");             /* x^2 - 2 */
    check("AlgebraicNumberTrace[2 Sqrt[2]]", "0");           /* x^2 - 8 */
    check("AlgebraicNumberTrace[5 + Sqrt[2]]", "10");        /* x^2 - 10 x + 23 */
    check("AlgebraicNumberTrace[9 + Sqrt[10]]", "18");       /* x^2 - 18 x + 71 */
    check("AlgebraicNumberTrace[GoldenRatio]", "1");         /* x^2 - x - 1 */
    check("AlgebraicNumberTrace[(1 + Sqrt[5])/2]", "1");     /* same, spelled out */
    check("AlgebraicNumberTrace[(1 + Sqrt[2])/2]", "1");     /* 4 x^2 - 4 x - 1 */
    check("AlgebraicNumberTrace[1/Sqrt[Sqrt[2] + 3]]", "0"); /* 7 x^4 - 6 x^2 + 1 */
    check("AlgebraicNumberTrace[Sqrt[1 + Sqrt[2]]]", "0");   /* x^4 - 2 x^2 - 1 */
    check("AlgebraicNumberTrace[Sqrt[2] + Sqrt[3]]", "0");   /* x^4 - 10 x^2 + 1 */
}

/* ---------------- Complex algebraic numbers / roots of unity --------------- */

static void test_complex(void) {
    check("AlgebraicNumberTrace[1 + I]", "2");               /* x^2 - 2 x + 2 */
    check("AlgebraicNumberTrace[I]", "0");                   /* x^2 + 1 */
    check("AlgebraicNumberTrace[E^(Pi I/8)]", "0");          /* Phi_16 = x^8 + 1 */
    check("AlgebraicNumberTrace[E^(2 Pi I)/3]", "1/3");      /* E^(2 Pi I) = 1 */
}

/* ---------------- Root and AlgebraicNumber objects ------------------------- */

static void test_objects(void) {
    check("AlgebraicNumberTrace[Root[-1 + #1 + #1^2 + #1^3 + #1^4 &, 1]]", "-1");
    check("AlgebraicNumberTrace[Root[#1^4 + 11 #1^3 + #1^2 + #1 + 1 &, 1]]", "-11");
    /* 1 + 2 Sqrt[2] I: minimal polynomial x^2 - 2 x + 9. */
    check("AlgebraicNumberTrace[AlgebraicNumber[Sqrt[2] I, {1, 2}]]", "2");
    check("AlgebraicNumberTrace[AlgebraicNumber[Sqrt[2], {0, 1}]]", "0");
}

/* ---------------- Listable threading --------------------------------------- */

static void test_listable(void) {
    check("AlgebraicNumberTrace[{5 + Sqrt[2], E^(Pi I/8)}]", "{10, 0}");
    check("AlgebraicNumberTrace[{2 Sqrt[2], E^(Pi I/8), 1 + I}]", "{0, 0, 2}");
    check("AlgebraicNumberTrace[{2, -2/3}]", "{2, -2/3}");
    check("Attributes[AlgebraicNumberTrace]", "{Listable, Protected}");
    check("Options[AlgebraicNumberTrace]", "{Extension -> None}");
}

/* ---------------- Relative trace (Extension option) ------------------------ */

static void test_extension(void) {
    /* Rational a in Q(theta): trace = a * [Q(theta):Q]. */
    check("AlgebraicNumberTrace[5, Extension -> Sqrt[2]]", "10");   /* 5 * 2 */
    check("AlgebraicNumberTrace[2, Extension -> Sqrt[5]]", "4");    /* 2 * 2 */
    check("AlgebraicNumberTrace[3, Extension -> Sqrt[2] + Sqrt[3]]", "12"); /* 3 * 4 */
    /* Trace of the generator itself: absolute trace, scaled by n/d = 1. */
    check("AlgebraicNumberTrace[Sqrt[5], Extension -> Sqrt[5]]", "0");
    check("AlgebraicNumberTrace[(1 + Sqrt[2])/2, Extension -> Sqrt[2]]", "1");
    /* Sqrt[2] in Q(zeta_8): d=2, n=4, absolute trace 0, so 0. */
    check("AlgebraicNumberTrace[Sqrt[2], Extension -> E^(Pi I/4)]", "0");
    /* Listable threading with a trailing Extension option. */
    check("AlgebraicNumberTrace[{5, Sqrt[5]}, Extension -> Sqrt[5]]", "{10, 0}");
    /* Extension -> None is the absolute trace (the default). */
    check("AlgebraicNumberTrace[5 + Sqrt[2], Extension -> None]", "10");
}

/* ---------------- Additivity over a fixed field ---------------------------- */
/* The relative trace Tr_{Q(theta)/Q} is additive: Tr(a + b) = Tr(a) + Tr(b). */

static void test_additive(void) {
    check("AlgebraicNumberTrace[5 + (1 + Sqrt[2])/2, Extension -> Sqrt[2]] == "
          "AlgebraicNumberTrace[5, Extension -> Sqrt[2]] + "
          "AlgebraicNumberTrace[(1 + Sqrt[2])/2, Extension -> Sqrt[2]]", "True");
    /* The spec's Total[] form of the same identity. */
    check("AlgebraicNumberTrace[5 + (1 + Sqrt[2])/2, Extension -> Sqrt[2]] == "
          "Total[AlgebraicNumberTrace[{5, (1 + Sqrt[2])/2}, Extension -> Sqrt[2]]]",
          "True");
    /* E^(2 Pi I/3) (degree 2) in the degree-6 compositum Q(2^(1/3), zeta_3):
     * trace = (6/2) * (-1) = -3. */
    check("AlgebraicNumberTrace[E^(2 Pi I/3), "
          "Extension -> ToNumberField[{2^(1/3), E^(2 Pi I/3)}, All][[1, 1]]]", "-3");
}

/* ---------------- Declines and messages ------------------------------------ */

static void test_declines_and_messages(void) {
    /* Non-algebraic argument: AlgebraicNumberTrace::nalg + unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("AlgebraicNumberTrace[Pi]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberTrace::nalg") != NULL,
               "expected AlgebraicNumberTrace::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberTrace[Pi]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* A free symbol declines with a message. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberTrace[x]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberTrace::nalg") != NULL,
               "expected AlgebraicNumberTrace::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberTrace[x]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* An argument not in the extension field: AlgebraicNumberTrace::ext. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberTrace[Sqrt[3], Extension -> Sqrt[2]]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberTrace::ext") != NULL,
               "expected AlgebraicNumberTrace::ext, got: %s", err);
    ASSERT_MSG(result &&
               strcmp(result, "AlgebraicNumberTrace[Sqrt[3], Extension -> Sqrt[2]]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* Wrong positional-argument count: silently unevaluated (no message). */
    check("AlgebraicNumberTrace[1, 2]", "AlgebraicNumberTrace[1, 2]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping "
               "AlgebraicNumberTrace tests.\n");
        return 0;
    }
    test_rationals();
    test_radicals();
    test_complex();
    test_objects();
    test_listable();
    test_extension();
    test_additive();
    test_declines_and_messages();
    printf("test_algebraicnumbertrace: all passed\n");
    return 0;
}
