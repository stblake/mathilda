/* Tests for the error function family: Erf[z] and Erf[z0, z1].
 *
 * Covers exact special values (0, +-1, the +-I directed infinities),
 * symbolic odd-symmetry, machine real & complex numerics, arbitrary-precision
 * (MPFR) reals and complex, precision tracking, the generalized two-argument
 * form, derivatives, Series, symbolic fall-through, Listable threading, and
 * attributes. */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- numeric helpers ------------------------------------------------ */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.12g, got %.12g", input, expected, v);
}

static void assert_complex_close(const char* input, double er, double ei, double tol) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_FUNCTION &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
               r->data.function.arg_count == 2,
               "%s: expected Complex[..], got something else", input);
    Expr* re = r->data.function.args[0];
    Expr* im = r->data.function.args[1];
    ASSERT(re->type == EXPR_REAL && im->type == EXPR_REAL);
    ASSERT_MSG(fabs(re->data.real - er) <= tol && fabs(im->data.real - ei) <= tol,
               "%s: expected %.8g %+.8g I, got %.8g %+.8g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact special values ------------------------------------------- */

void test_erf_exact() {
    assert_eval_eq("Erf[0]", "0", 0);
    assert_eval_eq("Erf[Infinity]", "1", 0);
    assert_eval_eq("Erf[-Infinity]", "-1", 0);
    assert_eval_eq("Erf[ComplexInfinity]", "ComplexInfinity", 0);
    assert_eval_eq("Erf[Indeterminate]", "Indeterminate", 0);
    /* Directed imaginary infinities. */
    assert_eval_eq("Erf[I Infinity]", "DirectedInfinity[Complex[0, 1]]", 1);
    assert_eval_eq("Erf[-I Infinity]", "DirectedInfinity[Complex[0, -1]]", 1);
}

/* ---- symbolic fall-through & odd symmetry --------------------------- */

void test_erf_symbolic() {
    assert_eval_eq("Erf[x]", "Erf[x]", 0);
    /* Odd: Erf[-x] = -Erf[x], Erf[-2 x] = -Erf[2 x]. */
    assert_eval_eq("Erf[-x]", "-Erf[x]", 0);
    assert_eval_eq("Erf[-2 x]", "-Erf[2 x]", 0);
    /* A bare symbolic constant stays put. */
    assert_eval_eq("Erf[a + b]", "Erf[a + b]", 0);
}

/* ---- machine real --------------------------------------------------- */

void test_erf_machine_real() {
    assert_close("Erf[0.95]", 0.8208908072732779, 1e-12);
    assert_close("Erf[1.5]", 0.9661051464753107, 1e-12);
    assert_close("Erf[-1.0]", -0.8427007929497149, 1e-12);
    assert_close("Erf[3.0]", 0.9999779095030014, 1e-12);
    assert_close("Erf[0.0]", 0.0, 0.0);
    /* Odd at machine precision. */
    assert_close("Erf[-2.5]", -erf(2.5), 1e-12);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_erf_arbitrary_precision() {
    /* Reference digits from mpmath (prefix kept short of the last-digit
     * rounding boundary). N[Erf[3/2], 50]. */
    assert_eval_startswith("N[Erf[3/2], 50]",
        "0.9661051464753107270669762616459478586814104");
    /* Erf[1/2]. */
    assert_eval_startswith("N[Erf[1/2], 30]", "0.5204998778130465376827466538");
    /* The precision of the output tracks the precision of the input. */
    assert_eval_startswith("Erf[1.500000000000000000000000`28]",
                           "0.966105146475310727066976");
    /* Required MPFR coverage case (numeric-builtins-cover-MPFR lesson). */
    assert_eval_startswith("N[Erf[1], 35]",
        "0.842700792949714869341220635082");
}

/* ---- machine complex ------------------------------------------------ */

void test_erf_machine_complex() {
    /* Erf[1.5 - I] = 1.0784 + 0.0279637 I (documentation; full value mpmath). */
    assert_complex_close("Erf[1.5 - I]", 1.0783992074989335, 0.0279637112386558, 1e-10);
    /* Erf[I] = I erfi(1) = 1.6504257587975428 I (pure imaginary). */
    assert_complex_close("Erf[1.0 I]", 0.0, 1.6504257587975428, 1e-10);
    /* erf(1 + i). */
    assert_complex_close("Erf[1.0 + 1.0 I]", 1.3161512816979477, 0.1904534692378347, 1e-10);
    /* Conjugate symmetry: Erf[conj z] = conj Erf[z]. */
    assert_complex_close("Erf[0.5 - 1.0 I]", 1.2048475583142180, -1.0244008816084459, 1e-10);
}

/* ---- arbitrary-precision complex ------------------------------------ */

void test_erf_arbitrary_complex() {
    /* Erf[I] to 30 digits is exactly I erfi(1) (mpmath reference). */
    /* Erf[I] = I Erfi[1]: pure imaginary, the real part an exact 0 */
    assert_eval_startswith("N[Erf[I], 30]",
        "1.65042575879754287602533772956*I");
    /* Erf[1/2 + I], mpmath reference. */
    assert_eval_startswith("N[Erf[1/2 + I], 30]",
        "1.20484755831421800270211268210");
}

/* ---- |Im z| beyond the machine kernel: extended exponent, not Overflow -- */

void test_erf_large_complex() {
    /* The machine kernel declines |Im z| > 26 and |erf| ~ 1.7e385 exceeds the
     * double range; Erf returned Overflow[] while Erfc/Erfi did not. mpmath:
     * erf(3+30i) = -1.24897783266730786e385 - 1.14353681145367294e385 i. */
    assert_eval_eq("Abs[Erf[3. + 30. I]/(-1.24897783266730786 10^385 - "
                   "1.14353681145367294 10^385 I) - 1] < 10^-13", "True", 0);
    /* |z| = 100 (arb kernel; the series cap used to run out): mpmath
     * erf(0.5+100i) = -1.97618773074101295e4340 + 3.32714786217951006e4340 i. */
    assert_eval_eq("Abs[Erf[0.5 + 100. I]/(-1.97618773074101295 10^4340 + "
                   "3.32714786217951006 10^4340 I) - 1] < 10^-13", "True", 0);
    assert_eval_eq("Abs[Erfi[100. + 0.5 I]/(3.32714786217951006 10^4340 - "
                   "1.97618773074101295 10^4340 I) - 1] < 10^-13", "True", 0);
    /* Moderate |z| through arb still lands on machine Reals. */
    assert_complex_close("Erf[3. + 20. I]", 8.15453402888295e167, 1.60442324870228e168, 1e154);
    assert_eval_eq("MachineNumberQ[Re[Erf[3. + 20. I]]]", "True", 0);
}

/* ---- two-argument generalized form ---------------------------------- */

void test_erf_two_arg() {
    /* Erf[1.5, 2] = erf(2) - erf(1.5) = 0.0292171. */
    assert_close("Erf[1.5, 2]", erf(2.0) - erf(1.5), 1e-9);
    assert_close("Erf[2.0, 3.0]", erf(3.0) - erf(2.0), 1e-9);
    /* Erf[-Infinity, Infinity] = 1 - (-1) = 2. */
    assert_eval_eq("Erf[-Infinity, Infinity]", "2", 0);
    /* Exact / symbolic two-arg stays unevaluated. */
    assert_eval_eq("Erf[2, 3]", "Erf[2, 3]", 0);
    assert_eval_eq("Erf[a, b]", "Erf[a, b]", 0);
}

/* ---- derivatives & series ------------------------------------------- */

void test_erf_derivatives() {
    assert_eval_eq("D[Erf[x], x]", "2/(E^x^2 Sqrt[Pi])", 0);
    /* Chain rule. */
    assert_eval_eq("D[Erf[x^2], x]", "(4 x)/(E^x^4 Sqrt[Pi])", 0);
    /* Series at the origin (generic D-based Taylor fallback). */
    assert_eval_startswith("Series[Erf[x], {x, 0, 3}]", "(2 x)/Sqrt[Pi] - (2 x^3)/(3 Sqrt[Pi])");
}

/* ---- Listable & attributes ------------------------------------------ */

void test_erf_listable() {
    assert_eval_eq("Erf[{0}]", "{0}", 0);
    assert_eval_eq("Erf[{}]", "{}", 0);
    /* Numeric threading. */
    Expr* e = parse_expression("Erf[{0.5, 1.0, 1.5}]");
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT(r->type == EXPR_FUNCTION &&
           strcmp(r->data.function.head->data.symbol.name, "List") == 0 &&
           r->data.function.arg_count == 3);
    double exp0[3] = { erf(0.5), erf(1.0), erf(1.5) };
    for (int i = 0; i < 3; i++) {
        Expr* el = r->data.function.args[i];
        ASSERT(el->type == EXPR_REAL);
        ASSERT_MSG(fabs(el->data.real - exp0[i]) <= 1e-12,
                   "Erf list element %d: expected %.12g got %.12g",
                   i, exp0[i], el->data.real);
    }
    expr_free(r);
}

void test_erf_attributes() {
    SymbolDef* d = symtab_get_def("Erf");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "Erf must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "Erf must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "Erf must be Protected");
}

/* ---- arity errors --------------------------------------------------- */

void test_erf_arity() {
    /* Wrong arity stays unevaluated (argt diagnostic to stderr). */
    assert_eval_eq("Erf[]", "Erf[]", 0);
    assert_eval_eq("Erf[1, 2, 3]", "Erf[1, 2, 3]", 0);
}

/* The machine-precision value comes from the shared raw-double kernel
 * (sf_machine_error.c) on every surface: scalar == unpacked List == packed
 * List == NDArray == Compile (scalar and rank-1), bit for bit, for machine-real and machine-complex arguments. */
void test_erf_machine_surfaces() {
    assert_eval_eq(
        "Module[{xs = (SeedRandom[1]; RandomReal[{-4, 4}, 200]), sc, cf, cfa}, sc = Erf /@ xs;"
        " cf = Compile[{{x, _Real}}, Erf[x]]; cfa = Compile[{{v, _Real, 1}}, Erf[v]];"
        " {Max[Abs[sc - Erf[FromPackedArray[xs]]]], Max[Abs[sc - Erf[ToPackedArray[xs]]]],"
        "  Max[Abs[sc - FromNDArray[Erf[ToNDArray[xs]]]]], Max[Abs[sc - (cf /@ xs)]],"
        "  Max[Abs[sc - cfa[xs]]]}]",
        "{0.0, 0.0, 0.0, 0.0, 0.0}", 0);
    assert_eval_eq(
        "Module[{xs = (SeedRandom[2]; RandomComplex[{-3 - 3 I, 3 + 3 I}, 200]), sc, cf, cfa},"
        " sc = Erf /@ xs; cf = Compile[{{x, _Complex}}, Erf[x]]; cfa = Compile[{{v, _Complex, 1}}, Erf[v]];"
        " {Max[Abs[sc - Erf[FromPackedArray[xs]]]], Max[Abs[sc - Erf[ToPackedArray[xs]]]],"
        "  Max[Abs[sc - FromNDArray[Erf[ToNDArray[xs]]]]], Max[Abs[sc - (cf /@ xs)]],"
        "  Max[Abs[sc - cfa[xs]]]}]",
        "{0.0, 0.0, 0.0, 0.0, 0.0}", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_erf_exact);
    TEST(test_erf_symbolic);
    TEST(test_erf_machine_real);
    TEST(test_erf_arbitrary_precision);
    TEST(test_erf_machine_complex);
    TEST(test_erf_arbitrary_complex);
    TEST(test_erf_large_complex);
    TEST(test_erf_two_arg);
    TEST(test_erf_derivatives);
    TEST(test_erf_listable);
    TEST(test_erf_attributes);
    TEST(test_erf_arity);

    TEST(test_erf_machine_surfaces);

    printf("All Erf tests passed.\n");
    return 0;
}
