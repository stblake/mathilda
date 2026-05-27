
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

void test_trig_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[0]", "0"},
        {"Sin[Pi/2]", "1"},
        {"Sin[Pi/6]", "Rational[1, 2]"},
        {"Sin[Pi/4]", "Power[2, Rational[-1, 2]]"},
        {"Sin[Pi/3]", "Times[Rational[1, 2], Power[3, Rational[1, 2]]]"},
        {"Cos[0]", "1"},
        {"Cos[Pi/2]", "0"},
        {"Cos[Pi/3]", "Rational[1, 2]"},
        {"Tan[Pi/4]", "1"},
        {"Cot[Pi/4]", "1"},
        {"Sec[0]", "1"},
        {"Csc[Pi/2]", "1"},
        // Exact values for d=12
        {"Sin[Pi/12]", "Times[Rational[1, 4], Plus[Power[6, Rational[1, 2]], Times[-1, Power[2, Rational[1, 2]]]]]"},
        {"Cos[Pi/12]", "Times[Rational[1, 4], Plus[Power[2, Rational[1, 2]], Power[6, Rational[1, 2]]]]"},
        // Exact values for d=10, 5
        {"Sin[Pi/10]", "Times[Rational[1, 4], Plus[-1, Power[5, Rational[1, 2]]]]"},
        {"Cos[Pi/5]", "Times[Rational[1, 4], Plus[1, Power[5, Rational[1, 2]]]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_inverse() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"ArcSin[0]", "0"},
        {"ArcSin[1]", "Times[Rational[1, 2], Pi]"},
        {"ArcSin[1/2]", "Times[Rational[1, 6], Pi]"},
        {"ArcCos[0]", "Times[Rational[1, 2], Pi]"},
        {"ArcCos[1]", "0"},
        {"ArcTan[0]", "0"},
        {"ArcTan[1]", "Times[Rational[1, 4], Pi]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        printf("Testing inverse: %s\n", cases[i].input);
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

void test_trig_forward_of_inverse() {
    /* f[f_inv[x]] = x for each direct trig / inverse-trig pair. Uses
     * InputForm-style strings since the simplification strips the inverse
     * and leaves the inner argument untouched. */
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sin[ArcSin[x]]", "x"},
        {"Cos[ArcCos[y]]", "y"},
        {"Tan[ArcTan[a]]", "a"},
        {"Cot[ArcCot[b]]", "b"},
        {"Sec[ArcSec[c]]", "c"},
        {"Csc[ArcCsc[d]]", "d"},
        /* Composite arguments survive as-is */
        {"Sin[ArcSin[x^2 + 1]]", "1 + x^2"},
        {"Cos[ArcCos[2 y]]", "2 y"},
        /* ArcTan[x, y] two-arg form must NOT be stripped: our rule guards
         * on arg_count==1, so Tan[ArcTan[x, y]] stays unevaluated here
         * rather than collapsing via the (wrong) single-argument rule. */
        {"Tan[ArcTan[3, 4]]", "Tan[ArcTan[3, 4]]"},
        /* Opposite direction is NOT folded: ArcSin[Sin[x]] stays put */
        {"ArcSin[Sin[x]]", "ArcSin[Sin[x]]"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0,
                   "Forward-of-inverse %s: expected %s, got %s",
                   cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(res);
        expr_free(e);
    }
}

/* Phase 4: Sin / Cos / Tan / Cot / Sec / Csc on Complex[MPFR, MPFR] must
 * preserve MPFR precision rather than coerce to double via csin/ccos/etc.
 *
 * The asserted prefixes are the first ~40 digits of the closed-form
 * value computed at 50-digit precision; precision is verified via a
 * Precision[...] round-trip (which itself produces an MPFR with the
 * expected magnitude, ~50.27 decimal digits for a 50-digit-input
 * result). */
void test_trig_mpfr_complex(void) {
    struct {
        const char* input;
        const char* expected_prefix;
        size_t prefix_len;
    } cases[] = {
        /* sin(1+i) = sin(1)cosh(1) + i cos(1)sinh(1)
         *          ~ 1.29846 + 0.63496 i */
        {"Sin[Complex[N[1, 50], N[1, 50]]]",
         "1.29845758141597729482604236580781562031343656163522", 40},
        /* cos(1+i) = cos(1)cosh(1) - i sin(1)sinh(1) */
        {"Cos[Complex[N[1, 50], N[1, 50]]]",
         "0.83373002513114904888388539433509447980987478520962", 40},
        /* tan(1+i) ~ 0.27175 + 1.08392 i */
        {"Tan[Complex[N[1, 50], N[1, 50]]]",
         "0.27175258531951171652884372249858892070946411146178", 40},
        /* cot(1+i) = 1/tan(1+i) ~ 0.21762 - 0.86801 i */
        {"Cot[Complex[N[1, 50], N[1, 50]]]",
         "0.21762156185440268136513424360523807352075436916785", 40},
        /* sec(1+i) = 1/cos(1+i) ~ 0.49834 + 0.59108 i */
        {"Sec[Complex[N[1, 50], N[1, 50]]]",
         "0.49833703055518678521380589177216953443287793247109", 40},
        /* csc(1+i) = 1/sin(1+i) ~ 0.62152 - 0.30393 i */
        {"Csc[Complex[N[1, 50], N[1, 50]]]",
         "0.62151801717042842123490780585592014816751214181073", 40},
        {NULL, NULL, 0}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* r = evaluate(e);
        char* s = expr_to_string(r);
        ASSERT_MSG(strncmp(s, cases[i].expected_prefix, cases[i].prefix_len) == 0,
                   "%s: expected leading %s, got %s",
                   cases[i].input, cases[i].expected_prefix, s);
        free(s);
        expr_free(e);
        expr_free(r);
    }

    /* Spot-check precision is preserved for Sin (representative). */
    Expr* ep = parse_expression(
        "Precision[Sin[Complex[N[1, 50], N[1, 50]]]]");
    Expr* rp = evaluate(ep);
    char* sp = expr_to_string(rp);
    ASSERT_MSG(strncmp(sp, "50.", 3) == 0,
               "Precision[Sin[1+i]] (50 digits): expected 50.*, got %s", sp);
    free(sp); expr_free(ep); expr_free(rp);

    /* Real MPFR input still takes the real path (no I in result). */
    Expr* eR = parse_expression("Sin[N[1, 50]]");
    Expr* rR = evaluate(eR);
    char* sR = expr_to_string(rR);
    ASSERT_MSG(strchr(sR, 'I') == NULL,
               "Sin[N[1, 50]]: expected pure real, got %s", sR);
    ASSERT_MSG(strncmp(sR, "0.8414709848078965066525023216302989996225",
                       42) == 0,
               "Sin[N[1, 50]] (50 digits): expected leading 0.84147..., got %s",
               sR);
    free(sR); expr_free(eR); expr_free(rR);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_trig_forward);
    TEST(test_trig_inverse);
    TEST(test_trig_forward_of_inverse);
    TEST(test_trig_mpfr_complex);

    return 0;
}
