
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

void test_hyperbolic_forward() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sinh[0]", "0"},
        {"Cosh[0]", "1"},
        {"Tanh[0]", "0"},
        {"Sinh[Infinity]", "Infinity"},
        {"Cosh[Infinity]", "Infinity"},
        {"Sinh[-Infinity]", "Times[-1, Infinity]"},
        {"Tanh[Infinity]", "1"},
        {"Tanh[-Infinity]", "-1"},
        {"Coth[Infinity]", "1"},
        {"Sech[Infinity]", "0"},
        {"Csch[Infinity]", "0"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

void test_hyperbolic_inverse() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"ArcSinh[0]", "0"},
        {"ArcCosh[1]", "0"},
        {"ArcTanh[0]", "0"},
        {"ArcCoth[Infinity]", "0"},
        {"ArcSech[1]", "0"},
        {"ArcSinh[Infinity]", "Infinity"},
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Inverse %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

void test_hyperbolic_forward_of_inverse() {
    /* f[f_inv[x]] = x for each direct hyperbolic / inverse-hyperbolic
     * pair. Identity holds over the complex numbers since each ArcXh is a
     * right inverse of Xh by construction. */
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        {"Sinh[ArcSinh[x]]", "x"},
        {"Cosh[ArcCosh[y]]", "y"},
        {"Tanh[ArcTanh[z]]", "z"},
        {"Coth[ArcCoth[w]]", "w"},
        {"Sech[ArcSech[u]]", "u"},
        {"Csch[ArcCsch[v]]", "v"},
        {"Sinh[ArcSinh[x^2 - 1]]", "-1 + x^2"},
        /* Opposite direction is NOT folded */
        {"ArcSinh[Sinh[x]]", "ArcSinh[Sinh[x]]"},
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
        expr_free(e);
        expr_free(res);
    }
}

/* Phase 5: forward hyperbolic functions on Complex[MPFR, MPFR] must
 * stay at MPFR precision. The identity sinh(z) = -i sin(i z) means
 * the closed-form values are the trig values shifted by a 90-degree
 * coordinate swap: sinh(a+bi) = sinh a cos b + i cosh a sin b, etc. */
void test_hyperbolic_mpfr_complex(void) {
    struct {
        const char* input;
        const char* expected_prefix;
        size_t prefix_len;
    } cases[] = {
        /* sinh(1+i) ~ 0.63496 + 1.29846 i */
        {"Sinh[Complex[N[1, 50], N[1, 50]]]",
         "0.63496391478473610825508220299150978151708195141938", 40},
        /* cosh(1+i) ~ 0.83373 + 0.98890 i */
        {"Cosh[Complex[N[1, 50], N[1, 50]]]",
         "0.83373002513114904888388539433509447980987478520963", 40},
        /* tanh(1+i) ~ 1.08392 + 0.27175 i */
        {"Tanh[Complex[N[1, 50], N[1, 50]]]",
         "1.08392332733869454347575206121197172134496752747539", 40},
        /* coth(1+i) = 1/tanh(1+i) ~ 0.86801 - 0.21762 i */
        {"Coth[Complex[N[1, 50], N[1, 50]]]",
         "0.86801414289592494863584920891627388827343874994608", 40},
        /* sech(1+i) = 1/cosh(1+i) ~ 0.49834 - 0.59108 i */
        {"Sech[Complex[N[1, 50], N[1, 50]]]",
         "0.49833703055518678521380589177216953443287793247109", 40},
        /* csch(1+i) = 1/sinh(1+i) ~ 0.30393 - 0.62152 i */
        {"Csch[Complex[N[1, 50], N[1, 50]]]",
         "0.30393100162842645033448560450970327348872988190146", 40},
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

    Expr* ep = parse_expression(
        "Precision[Sinh[Complex[N[1, 50], N[1, 50]]]]");
    Expr* rp = evaluate(ep);
    char* sp = expr_to_string(rp);
    ASSERT_MSG(strncmp(sp, "50.", 3) == 0,
               "Precision[Sinh[1+i]] (50 digits): expected 50.*, got %s", sp);
    free(sp); expr_free(ep); expr_free(rp);

    /* Pure real MPFR stays pure real (no I in the printed result). */
    Expr* eR = parse_expression("Cosh[N[1, 50]]");
    Expr* rR = evaluate(eR);
    char* sR = expr_to_string(rR);
    ASSERT_MSG(strchr(sR, 'I') == NULL,
               "Cosh[N[1, 50]]: expected pure real, got %s", sR);
    free(sR); expr_free(eR); expr_free(rR);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hyperbolic_forward);
    TEST(test_hyperbolic_inverse);
    TEST(test_hyperbolic_forward_of_inverse);
    TEST(test_hyperbolic_mpfr_complex);

    return 0;
}
