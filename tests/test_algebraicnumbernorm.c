/* test_algebraicnumbernorm.c — AlgebraicNumberNorm[a].
 *
 * AlgebraicNumberNorm[a] is the field norm of the algebraic number a: the product
 * of the roots of a's minimal polynomial, i.e. (-1)^deg times the monic-over-Q
 * constant term.  AlgebraicNumberNorm[a, Extension -> theta] is the relative norm
 * over Q(theta), equal to the absolute norm raised to the tower index [Q(theta):Q(a)]
 * (defined only when a lies in Q(theta)).
 *
 * Covers integers/rationals, real and complex radicals, roots of unity, Root and
 * AlgebraicNumber objects, GoldenRatio (the sole algebraic named constant, handled
 * by the shared qqbar converter), Listable threading with a trailing Extension
 * option, the relative-norm cases (including multiplicativity over a fixed field
 * and a compositum built by ToNumberField), and the decline paths
 * (AlgebraicNumberNorm::nalg on a non-algebraic argument, AlgebraicNumberNorm::ext
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
    const char* path = "/tmp/mathilda_annorm_stderr.log";
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
    check("AlgebraicNumberNorm[2]", "2");
    check("AlgebraicNumberNorm[-2/3]", "-2/3");
    check("AlgebraicNumberNorm[0]", "0");
    check("AlgebraicNumberNorm[1]", "1");
    check("AlgebraicNumberNorm[-5]", "-5");
    check("AlgebraicNumberNorm[7/10]", "7/10");
    /* Norm of a rational is the rational itself (minimal polynomial q x - p). */
    check("AlgebraicNumberNorm[6/4]", "3/2");
}

/* ---------------- Real radicals and GoldenRatio ---------------------------- */

static void test_radicals(void) {
    check("AlgebraicNumberNorm[Sqrt[2]]", "-2");             /* x^2 - 2 */
    check("AlgebraicNumberNorm[2 Sqrt[2]]", "-8");           /* x^2 - 8 */
    check("AlgebraicNumberNorm[GoldenRatio]", "-1");         /* x^2 - x - 1 */
    check("AlgebraicNumberNorm[(1 + Sqrt[5])/2]", "-1");     /* same, spelled out */
    check("AlgebraicNumberNorm[1/Sqrt[Sqrt[2] + 3]]", "1/7");/* 7 x^4 - 6 x^2 + 1 */
    check("AlgebraicNumberNorm[9 + Sqrt[10]]", "71");        /* x^2 - 18 x + 71 */
    check("AlgebraicNumberNorm[Sqrt[2] + Sqrt[3]]", "1");    /* x^4 - 10 x^2 + 1 */
}

/* ---------------- Complex algebraic numbers / roots of unity --------------- */

static void test_complex(void) {
    check("AlgebraicNumberNorm[1 + I]", "2");                /* x^2 - 2 x + 2 */
    check("AlgebraicNumberNorm[I]", "1");                    /* x^2 + 1 */
    check("AlgebraicNumberNorm[E^(Pi I/8)]", "1");           /* Phi_16 = x^8 + 1 */
    check("AlgebraicNumberNorm[E^(2 Pi I)/3]", "1/3");       /* E^(2 Pi I) = 1 */
}

/* ---------------- Root and AlgebraicNumber objects ------------------------- */

static void test_objects(void) {
    check("AlgebraicNumberNorm[Root[-1 + #1 + #1^2 + #1^3 + #1^4 &, 1]]", "-1");
    /* 1 + 2 Sqrt[2] I: minimal polynomial x^2 - 2 x + 9. */
    check("AlgebraicNumberNorm[AlgebraicNumber[Sqrt[2] I, {1, 2}]]", "9");
    check("AlgebraicNumberNorm[AlgebraicNumber[Sqrt[2], {0, 1}]]", "-2");
}

/* ---------------- Listable threading --------------------------------------- */

static void test_listable(void) {
    check("AlgebraicNumberNorm[{2 Sqrt[2], E^(Pi I/8), 1 + I}]", "{-8, 1, 2}");
    check("AlgebraicNumberNorm[{2, -2/3}]", "{2, -2/3}");
    check("Attributes[AlgebraicNumberNorm]", "{Listable, Protected}");
    check("Options[AlgebraicNumberNorm]", "{Extension -> None}");
}

/* ---------------- Relative norm (Extension option) ------------------------- */

static void test_extension(void) {
    /* Sqrt[2] in Q(zeta_8): d=2, n=4, (-2)^(4/2) = 4. */
    check("AlgebraicNumberNorm[Sqrt[2], Extension -> E^(Pi I/4)]", "4");
    /* Norm over Q(Sqrt[5]): a rational r has norm r^2. */
    check("AlgebraicNumberNorm[2, Extension -> Sqrt[5]]", "4");
    check("AlgebraicNumberNorm[Sqrt[5], Extension -> Sqrt[5]]", "-5");
    check("AlgebraicNumberNorm[2 Sqrt[5], Extension -> Sqrt[5]]", "-20");
    /* Listable threading with a trailing Extension option. */
    check("AlgebraicNumberNorm[{2, Sqrt[5]}, Extension -> Sqrt[5]]", "{4, -5}");
    /* Compositum Q(Sqrt[3], Sqrt[-5]) via ToNumberField; norm of Sqrt[3]+Sqrt[-5]. */
    check("AlgebraicNumberNorm[Sqrt[3] + Sqrt[-5], "
          "Extension -> ToNumberField[{Sqrt[3], Sqrt[-5]}, All][[1, 1]]]", "64");
    /* Extension -> None is the absolute norm (the default). */
    check("AlgebraicNumberNorm[Sqrt[2], Extension -> None]", "-2");
}

/* ---------------- Multiplicativity over a fixed field ---------------------- */
/* The relative norm N_{Q(theta)/Q} is multiplicative: N(a b) = N(a) N(b). */

static void test_multiplicative(void) {
    /* -20 == 4 * (-5) over Q(Sqrt[5]). */
    check("AlgebraicNumberNorm[2 Sqrt[5], Extension -> Sqrt[5]] == "
          "AlgebraicNumberNorm[2, Extension -> Sqrt[5]] "
          "AlgebraicNumberNorm[Sqrt[5], Extension -> Sqrt[5]]", "True");
    /* Absolute norm is multiplicative within a single field too. */
    check("AlgebraicNumberNorm[(1 + Sqrt[5])/2, Extension -> Sqrt[5]] "
          "AlgebraicNumberNorm[(1 - Sqrt[5])/2, Extension -> Sqrt[5]] == "
          "AlgebraicNumberNorm[-1, Extension -> Sqrt[5]]", "True");
}

/* ---------------- Declines and messages ------------------------------------ */

static void test_declines_and_messages(void) {
    /* Non-algebraic argument: AlgebraicNumberNorm::nalg + unevaluated. */
    char* result = NULL;
    char* err = eval_capturing_stderr("AlgebraicNumberNorm[Pi]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberNorm::nalg") != NULL,
               "expected AlgebraicNumberNorm::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberNorm[Pi]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* A free symbol declines with a message. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberNorm[x]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberNorm::nalg") != NULL,
               "expected AlgebraicNumberNorm::nalg, got: %s", err);
    ASSERT_MSG(result && strcmp(result, "AlgebraicNumberNorm[x]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* An argument not in the extension field: AlgebraicNumberNorm::ext. */
    result = NULL;
    err = eval_capturing_stderr("AlgebraicNumberNorm[Sqrt[3], Extension -> Sqrt[2]]", &result);
    ASSERT(err != NULL);
    ASSERT_MSG(strstr(err, "AlgebraicNumberNorm::ext") != NULL,
               "expected AlgebraicNumberNorm::ext, got: %s", err);
    ASSERT_MSG(result &&
               strcmp(result, "AlgebraicNumberNorm[Sqrt[3], Extension -> Sqrt[2]]") == 0,
               "expected the call to stay unevaluated, got: %s", result ? result : "(null)");
    free(err);
    free(result);

    /* Wrong positional-argument count: silently unevaluated (no message). */
    check("AlgebraicNumberNorm[1, 2]", "AlgebraicNumberNorm[1, 2]");
}

int main(void) {
    symtab_init();
    core_init();
    if (!flint_bridge_available()) {
        printf("FLINT not compiled in (USE_FLINT off); skipping "
               "AlgebraicNumberNorm tests.\n");
        return 0;
    }
    test_rationals();
    test_radicals();
    test_complex();
    test_objects();
    test_listable();
    test_extension();
    test_multiplicative();
    test_declines_and_messages();
    printf("test_algebraicnumbernorm: all passed\n");
    return 0;
}
