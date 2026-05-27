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

/* Capture stderr while `input` is parsed + evaluated. Returns the
 * collected stderr text as a heap string (caller frees) and writes the
 * printed result into *out_result_str (also heap-allocated). Same
 * pattern as tests/test_integer_exponent.c. */
static char* eval_capturing_stderr(const char* input, char** out_result_str) {
    const char* path = "/tmp/mathilda_log_stderr.log";
    fflush(stderr);
    if (!freopen(path, "w+", stderr)) {
        if (out_result_str) *out_result_str = NULL;
        return NULL;
    }

    Expr* p = parse_expression(input);
    Expr* e = evaluate(p);
    if (out_result_str) *out_result_str = expr_to_string(e);
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

/* Log[] and Log[a, b, c, ...] should emit `Log::argt` and leave the
 * call unevaluated, matching Mathematica's documented behaviour. */
void test_log_argt(void) {
    struct {
        const char* input;
        size_t argc;
    } cases[] = {
        {"Log[]",            0},
        {"Log[1, 2, 3]",     3},
        {"Log[2, 3, 4, 5]",  4},
        {NULL, 0}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        char* result = NULL;
        char* err = eval_capturing_stderr(cases[i].input, &result);
        ASSERT(result != NULL);
        ASSERT_MSG(strcmp(result, cases[i].input) == 0,
                   "argt %s: expected call unchanged, got %s",
                   cases[i].input, result);
        ASSERT_MSG(err && strstr(err, "Log::argt") != NULL,
                   "argt %s: expected Log::argt diagnostic, got: %s",
                   cases[i].input, err ? err : "(null)");
        char needle[64];
        snprintf(needle, sizeof needle, "called with %zu argument%s",
                 cases[i].argc, cases[i].argc == 1 ? "" : "s");
        ASSERT_MSG(strstr(err, needle) != NULL,
                   "argt %s: missing arg-count phrase '%s' in: %s",
                   cases[i].input, needle, err);
        ASSERT(strstr(err, "1 or 2 arguments are expected") != NULL);
        free(result);
        free(err);
    }
}

void test_logexp_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Log[E]", "1"},
        {"Log[1]", "0"},
        {"Log[0]", "-Infinity"},
        {"Log[0.0]", "Indeterminate"},
        {"Log[-0.0]", "Indeterminate"},
        {"Log[Infinity]", "Infinity"},
        {"Log[2, 8]", "3"},
        {"Log[2, 0]", "-Infinity"},
        {"Log[10, 0]", "-Infinity"},
        {"Log[E, 0]", "-Infinity"},
        {"Log[Pi, 0]", "-Infinity"},
        {"Log[1/2, 0]", "Infinity"},
        {"Log[2, 0.0]", "Indeterminate"},
        {"Log[Pi, 0.0]", "Indeterminate"},
        {"Log[b, b]", "1"},
        {"Exp[0]", "1"},
        {"Exp[-Infinity]", "0"},
        {"Exp[Infinity]", "Infinity"},
        {"Log[10, 10]", "1"},
        {"Log[E, x]", "Log[x]"},
        {"Exp[0.0]", "1.0"},
        {"Log[-1.0]", "0.0 + 3.14159*I"},
        {"Log[-1]", "I Pi"},
        {"Log[-5]", "Log[5] + I Pi"},
        {"Exp[Log[x]]", "x"},
        {"Exp[b Log[a]]", "a^b"},
        {"Exp[2 Log[x]]", "x^2"},
        {"Exp[3 Log[2]]", "8"},
        {"3^Log[3, x]", "x"},
        {"b^Log[b, a]", "a"},
        {"Power[base, b Log[base, a]]", "a^b"},
        {"10^(3 Log[10, x])", "x^3"},
        {"Log[E^4]", "4"},
        {"Log[E^(1/3)]", "1/3"},
        {"Log[E^(-2)]", "-2"},
        {"Log[E^x]", "Log[E^x]"},
        {"Log[2, 2^(1/3)]", "1/3"},
        {"Log[E, E^4]", "4"},
        {"Table[Exp[I * n * Pi / 2], {n, 0, 4}]", "{1, I, -1, -I, 1}"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

/* Phase 3: Exp[Complex[MPFR, MPFR]] and Log[Complex[MPFR, MPFR]] must
 * preserve MPFR precision end-to-end rather than coerce through libc's
 * cexp / clog. Regressions cover:
 *
 *   - Euler's identity Exp[I N[Pi, 50]] == -1 + 0 I, with the residual
 *     imaginary part at the MPFR roundoff floor (~1e-50, not ~1e-15).
 *   - Log of a complex MPFR matches the closed-form
 *     Log[Sqrt[2]] + I Pi/4 at 50 leading digits.
 *   - Log of a negative MPFR real promotes through the new complex
 *     branch (helper extracts the imag-zero MPFR into (re, 0); log of
 *     a negative becomes log(|re|) + I Pi at MPFR precision).
 *
 * Precision is verified via a Precision[...] round-trip — same idiom
 * used by the Abs / Arg / Sign tests in test_core.c. */
void test_exp_log_mpfr_complex(void) {
    /* Exp[Complex[N[1, 50], N[Pi, 50]]] = exp(1) (cos pi + i sin pi)
     * = -e + ~0 i. Both components are MPFR so the Complex head
     * survives Mathilda's `Complex[0, x] -> x I` auto-collapse — the
     * collapse fires when the real component is an exact Integer or
     * Real 0, never for an MPFR (which is "inexact" and not equal to
     * the exact 0). */
    Expr* e1 = parse_expression("Exp[Complex[N[1, 50], N[Pi, 50]]]");
    Expr* r1 = evaluate(e1);
    char* s1 = expr_to_string(r1);
    /* Leading "-2.71828..." matches exp(1)*cos(Pi) = -e at MPFR
     * precision; the imag residual is at the working-precision floor. */
    ASSERT_MSG(strncmp(s1, "-2.7182818284590452353602874713526624977572",
                       42) == 0,
               "Exp[1 + I Pi] (50 digits): expected leading -2.71828..., got %s",
               s1);
    free(s1); expr_free(e1); expr_free(r1);

    Expr* e1p = parse_expression("Precision[Exp[Complex[N[1, 50], N[Pi, 50]]]]");
    Expr* r1p = evaluate(e1p);
    char* s1p = expr_to_string(r1p);
    ASSERT_MSG(strncmp(s1p, "50.", 3) == 0,
               "Precision[Exp[1 + I Pi]] (50 digits): expected 50.*, got %s",
               s1p);
    free(s1p); expr_free(e1p); expr_free(r1p);

    /* Log[Complex[N[1, 80], N[1, 80]]] matches Log[Sqrt[2]] + I Pi/4
     * at 80 digits. Asserts the first 50 digits of the real part. */
    Expr* e2 = parse_expression("Log[Complex[N[1, 80], N[1, 80]]]");
    Expr* r2 = evaluate(e2);
    char* s2 = expr_to_string(r2);
    ASSERT_MSG(strncmp(s2, "0.34657359027997265470861606072908828403775006718012",
                       50) == 0,
               "Log[1+i] (80 digits): expected Log[Sqrt[2]] +..., got %s",
               s2);
    free(s2); expr_free(e2); expr_free(r2);

    Expr* e2p = parse_expression("Precision[Log[Complex[N[1, 80], N[1, 80]]]]");
    Expr* r2p = evaluate(e2p);
    char* s2p = expr_to_string(r2p);
    ASSERT_MSG(strncmp(s2p, "80.", 3) == 0,
               "Precision[Log[1+i]] (80 digits): expected 80.*, got %s",
               s2p);
    free(s2p); expr_free(e2p); expr_free(r2p);

    /* Log[N[-2, 50]] = Log[2] + I Pi at 50 digits — the existing real
     * MPFR Log path rejected negatives; the new complex branch handles
     * them via mpfr_complex_log = log(hypot) + I atan2. */
    Expr* e3 = parse_expression("Log[N[-2, 50]]");
    Expr* r3 = evaluate(e3);
    char* s3 = expr_to_string(r3);
    ASSERT_MSG(strncmp(s3, "0.6931471805599453094172321214581765680755001343602",
                       49) == 0,
               "Log[-2] (50 digits): expected Log[2] + I Pi..., got %s", s3);
    free(s3); expr_free(e3); expr_free(r3);

    Expr* e3p = parse_expression("Precision[Log[N[-2, 50]]]");
    Expr* r3p = evaluate(e3p);
    char* s3p = expr_to_string(r3p);
    ASSERT_MSG(strncmp(s3p, "50.", 3) == 0,
               "Precision[Log[-2]] (50 digits): expected 50.*, got %s", s3p);
    free(s3p); expr_free(e3p); expr_free(r3p);

    /* Pure real MPFR Exp still takes the real-only path: result is a
     * single EXPR_MPFR (no "I" in the printed form, no Complex[
     * wrapper). */
    Expr* e4p = parse_expression("Exp[N[1, 50]]");
    Expr* r4p = evaluate(e4p);
    char* s4p = expr_to_string(r4p);
    ASSERT_MSG(strncmp(s4p, "2.7182818284590452353602874713526624977572",
                       42) == 0,
               "Exp[N[1, 50]] (50 digits): expected leading 2.71828..., got %s",
               s4p);
    ASSERT_MSG(strchr(s4p, 'I') == NULL,
               "Exp[N[1, 50]]: expected pure real result, got %s", s4p);
    free(s4p); expr_free(e4p); expr_free(r4p);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_logexp_forward);
    TEST(test_log_argt);
    TEST(test_exp_log_mpfr_complex);

    return 0;
}
