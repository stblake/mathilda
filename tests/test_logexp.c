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

int main() {
    symtab_init();
    core_init();
    
    TEST(test_logexp_forward);
    TEST(test_log_argt);

    return 0;
}
