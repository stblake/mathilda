/* Tests for the Lerch transcendent LerchPhi[z, s, a].
 *
 * Covers the exact reductions (z = 0 -> a^-s, s = 0 -> 1/(1-z), z = 1 ->
 * Zeta[s,a], z = -1, positive integer a -> a PolyLog form, negative integer s
 * -> a rational function of z), machine real & complex numerics with
 * cross-consistency against PolyLog and Zeta, arbitrary-precision (MPFR) reals
 * and complexes, the DoublyInfinite and IncludeSingularTerm options, the
 * |z| > 1 symbolic fall-through, the derivative rules, Listable threading,
 * attributes, and argument / option diagnostics. */

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
               strcmp(r->data.function.head->data.symbol, "Complex") == 0 &&
               r->data.function.arg_count == 2,
               "%s: expected Complex[..], got something else", input);
    Expr* re = r->data.function.args[0];
    Expr* im = r->data.function.args[1];
    ASSERT(re->type == EXPR_REAL && im->type == EXPR_REAL);
    ASSERT_MSG(fabs(re->data.real - er) <= tol && fabs(im->data.real - ei) <= tol,
               "%s: expected %.9g %+.9g I, got %.9g %+.9g I",
               input, er, ei, re->data.real, im->data.real);
    expr_free(r);
}

/* ---- exact reductions ----------------------------------------------- */

void test_lp_trivial_reductions() {
    /* z = 0: only the k = 0 term, a^-s. */
    assert_eval_eq("LerchPhi[0, 2, 3]", "1/9", 0);
    assert_eval_eq("LerchPhi[0, 0, a]", "1", 0);
    assert_eval_eq("LerchPhi[0, s, a]", "a^(-s)", 0);
    /* s = 0: geometric sum 1/(1-z), independent of a. */
    assert_eval_eq("LerchPhi[z, 0, a]", "1/(1 - z)", 0);
    assert_eval_eq("LerchPhi[1/2, 0, 7]", "2", 0);
}

void test_lp_reduce_to_polylog() {
    /* z LerchPhi[z, s, 1] = PolyLog[s, z]. */
    assert_eval_eq("LerchPhi[z, s, 1]", "PolyLog[s, z]/z", 0);
    /* Positive integer a > 1: shift down to the PolyLog series. */
    assert_eval_eq("LerchPhi[z, s, 2]", "(-z + PolyLog[s, z])/z^2", 0);
}

void test_lp_reduce_to_zeta() {
    /* z = 1: LerchPhi[1, s, a] = Zeta[s, a]. */
    assert_eval_eq("LerchPhi[1, s, a]", "Zeta[s, a]", 0);
    assert_eval_eq("LerchPhi[1, 3, 1]", "Zeta[3]", 0);
    /* z = -1: 2^-s (Zeta[s, a/2] - Zeta[s, (a+1)/2]). */
    assert_eval_eq("LerchPhi[-1, s, 1/2]",
                   "2^(-s) (-Zeta[s, 3/4] + Zeta[s, 1/4])", 0);
}

void test_lp_nonpositive_integer_a() {
    /* a = 0: LerchPhi[z, s, 0] = PolyLog[s, z]. */
    assert_eval_eq("LerchPhi[z, s, 0]", "PolyLog[s, z]", 0);
    assert_eval_eq("PolyLog[s, z] == LerchPhi[z, s, 0]", "True", 0);
    /* Value at the origin: the singular k = 0 term is dropped. */
    assert_eval_eq("LerchPhi[0, s, 0]", "0", 0);
    assert_eval_eq("LerchPhi[0, 0, 0]", "0", 0);
    /* LerchPhi[1/2, 1, 0] = PolyLog[1, 1/2] = -Log[1/2] (= Log[2]). */
    assert_eval_eq("LerchPhi[1/2, 1, 0]", "-Log[1/2]", 0);
    /* Negative integer a = -1: z PolyLog[s,z] + (-1)^-s. */
    assert_eval_eq("LerchPhi[z, s, -1]", "(-1)^(-s) + z PolyLog[s, z]", 0);
}

void test_lp_minus_one_integer_a() {
    /* z = -1 with integer a routes through PolyLog (the two-Zeta form is
     * indeterminate at integer s). LerchPhi[-1, 1, 1] = Log[2]. */
    assert_eval_eq("LerchPhi[-1, 1, 1]", "Log[2]", 0);
    /* The reported divergent table: {Log[2], 1, ComplexInfinity, -I Pi/2}. */
    assert_eval_eq("Table[LerchPhi[z, 1, 1], {z, -1, 2}]",
                   "{Log[2], 1, ComplexInfinity, (-1/2*I) Pi}", 0);
}

void test_lp_series_at_zero() {
    /* Series at z = 0 with symbolic s, a: coefficient ((k+a)^2)^(-s/2). */
    assert_eval_eq("Normal[Series[LerchPhi[z, s, a], {z, 0, 3}]]",
                   "(a^2)^(-1/2 s) + ((1 + a)^2)^(-1/2 s) z + "
                   "((2 + a)^2)^(-1/2 s) z^2 + ((3 + a)^2)^(-1/2 s) z^3", 0);
    /* Reducible case (a positive integer) goes through PolyLog and matches. */
    assert_eval_eq("Normal[Series[LerchPhi[z, 1, 2], {z, 0, 4}]]",
                   "1/2 + 1/3 z + 1/4 z^2 + 1/5 z^3 + 1/6 z^4", 0);
}

void test_lp_negative_integer_s() {
    /* LerchPhi[z, -n, a] = (z d/dz + a)^n [1/(1-z)], a rational function. */
    assert_eval_eq("LerchPhi[2, -1, a]", "2 - a", 0);
    assert_eval_eq("LerchPhi[3, -1, a]", "1/4 (3 - 2 a)", 0);
    assert_eval_eq("LerchPhi[2, -2, a]", "2 (-3 + a) + a (2 - a)", 0);
    /* Rational z: LerchPhi[1/3, -1, 2] = Sum (k+2)(1/3)^k = 15/4. */
    assert_eval_eq("LerchPhi[1/3, -1, 2]", "15/4", 0);
}

/* ---- machine-precision numerics + cross-consistency ----------------- */

void test_lp_machine_real() {
    assert_close("LerchPhi[0.5, 3, 2.5]", 0.0794983119688099407, 1e-12);
    assert_close("LerchPhi[0.4, 2, 3]", 0.145046476114, 1e-9);
    /* z LerchPhi[z, s, 1] - PolyLog[s, z] = 0 (series vs PolyLog kernel). */
    assert_close("0.4 LerchPhi[0.4, 3, 1] - PolyLog[3, 0.4]", 0.0, 1e-12);
    /* LerchPhi[1, s, a] - Zeta[s, a] = 0 (series vs Zeta kernel). */
    assert_close("LerchPhi[1, 3.0, 2.5] - Zeta[3.0, 2.5]", 0.0, 1e-12);
}

void test_lp_machine_complex() {
    assert_complex_close("LerchPhi[0.3 + 0.2 I, 2, 1.5]",
                         0.495505453148633019, 0.0444652909672175351, 1e-9);
}

/* ---- arbitrary precision (MPFR) ------------------------------------- */

void test_lp_arbitrary_precision() {
    assert_eval_startswith("N[LerchPhi[1/2, 2, 5/2], 30]",
                           "0.219693113434910235649039949138");
    /* Cross-check against PolyLog at arbitrary precision (a = 1 reduction). */
    assert_eval_startswith("N[LerchPhi[3/10, 2, 1], 25]",
                           "1.08709836691825356510011");
    /* Inexact MPFR inputs (z and a carry 40-digit precision). */
    assert_eval_startswith("LerchPhi[0.5`40, 3, 2.5`40]",
                           "0.07949831196880994065085794");
    /* Irrational value to 20 digits. */
    assert_eval_startswith("N[LerchPhi[1/2, 3, 1/2], 20]",
                           "8.168022726140350633");
}

void test_lp_arbitrary_complex() {
    assert_eval_startswith("N[LerchPhi[3/10 + 2/10 I, 2, 3/2], 20]",
                           "0.49550545314863301");
}

/* ---- |z| > 1 numeric continuation (Lerch/Erdelyi expansion) --------- */

void test_lp_large_z_continuation() {
    /* Off the branch cut, non-integer s: Erdelyi continuation gives the value. */
    assert_complex_close("LerchPhi[5. + I, I, I + 2]", -0.581502, 0.384767, 1e-5);
    /* Independently cross-checked against PolyLog: z LerchPhi[z,s,1] = PolyLog[s,z]
     * for |z| > 1 (PolyLog has its own continuation), to ~machine epsilon. */
    assert_close("Abs[(2.+0.5 I) LerchPhi[2.+0.5 I, 1.5, 1] - PolyLog[1.5, 2.+0.5 I]]",
                 0.0, 1e-12);
    assert_close("Abs[(2.+0.5 I) LerchPhi[2.+0.5 I, 0.7+0.3 I, 1] "
                 "- PolyLog[0.7+0.3 I, 2.+0.5 I]]", 0.0, 1e-12);
}

/* ---- |z| > 1 stays symbolic where no continuation is implemented ----- */

void test_lp_large_z_symbolic() {
    /* Integer s on the branch cut needs the logarithmic confluent form. */
    assert_eval_eq("LerchPhi[2, 3, -1.5]", "LerchPhi[2, 3, -1.5]", 0);
    assert_eval_eq("LerchPhi[2.0, 3, 0.5]", "LerchPhi[2.0, 3, 0.5]", 0);
}

/* ---- options -------------------------------------------------------- */

void test_lp_doubly_infinite() {
    /* Symmetric case (a = 1/2) just doubles: 2 Zeta[2,1/2] = Pi^2. */
    assert_close("LerchPhi[1, 2, 0.5, DoublyInfinite -> True]", 9.86960440109, 1e-7);
    /* Default option value (False) leaves the single sum. */
    assert_close("LerchPhi[1, 2, 0.5, DoublyInfinite -> False]", 4.93480220054, 1e-7);
}

void test_lp_include_singular_term() {
    /* IncludeSingularTerm at a non-positive integer a -> the k = -a term blows
     * up for any s. */
    assert_eval_eq("LerchPhi[1/2, 3, -3, IncludeSingularTerm -> True]",
                   "ComplexInfinity", 0);
    assert_eval_eq("LerchPhi[2, 3, -4, IncludeSingularTerm -> True]",
                   "ComplexInfinity", 0);
}

/* ---- derivatives ---------------------------------------------------- */

void test_lp_derivatives() {
    /* d/dz Phi = (LerchPhi[z,-1+s,a] - a LerchPhi[z,s,a]) / z. */
    assert_eval_eq("D[LerchPhi[z, s, a], z]",
                   "(-a LerchPhi[z, s, a] + LerchPhi[z, -1 + s, a])/z", 0);
    /* d/da Phi = -s LerchPhi[z, 1+s, a]. */
    assert_eval_eq("D[LerchPhi[z, s, a], a]", "-s LerchPhi[z, 1 + s, a]", 0);
    /* d/ds has no elementary closed form: generic partial. */
    assert_eval_eq("D[LerchPhi[z, s, a], s]",
                   "Derivative[0, 1, 0][LerchPhi][z, s, a]", 0);
    /* Chain rule through the first argument. */
    assert_eval_eq("D[LerchPhi[x^2, s, a], x]",
                   "(2 (-a LerchPhi[x^2, s, a] + LerchPhi[x^2, -1 + s, a]))/x", 0);
    /* Chain rule through the third argument. */
    assert_eval_eq("D[LerchPhi[z, s, x^3], x]",
                   "-3 s x^2 LerchPhi[z, 1 + s, x^3]", 0);
    /* Higher a-derivatives follow the rising-factorial pattern. */
    assert_eval_eq("Table[D[LerchPhi[z, s, a], {a, k}], {k, 1, 2}]",
                   "{-s LerchPhi[z, 1 + s, a], s (1 + s) LerchPhi[z, 2 + s, a]}", 0);
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_lp_symbolic() {
    assert_eval_eq("LerchPhi[z, s, a]", "LerchPhi[z, s, a]", 0);
    assert_eval_eq("LerchPhi[z, 3, a]", "LerchPhi[z, 3, a]", 0);
}

/* ---- Listable & attributes ------------------------------------------ */

void test_lp_listable() {
    assert_eval_eq("LerchPhi[{1, 0}, 2, 1]", "{1/6 Pi^2, 1}", 0);
    assert_eval_eq("LerchPhi[z, 0, {1, 2}]", "{1/(1 - z), 1/(1 - z)}", 0);
}

void test_lp_attributes() {
    SymbolDef* d = symtab_get_def("LerchPhi");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "LerchPhi must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "LerchPhi must be NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "LerchPhi must be Protected");
}

void test_lp_argcount() {
    /* Fewer than 3 arguments -> LerchPhi::argrx, stays unevaluated. */
    assert_eval_eq("LerchPhi[]", "LerchPhi[]", 0);
    assert_eval_eq("LerchPhi[1, 2]", "LerchPhi[1, 2]", 0);
    /* A non-option beyond position 3 -> LerchPhi::nonopt, stays unevaluated. */
    assert_eval_eq("LerchPhi[1, 2, 3, 4, 5]", "LerchPhi[1, 2, 3, 4, 5]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_lp_trivial_reductions);
    TEST(test_lp_reduce_to_polylog);
    TEST(test_lp_reduce_to_zeta);
    TEST(test_lp_nonpositive_integer_a);
    TEST(test_lp_minus_one_integer_a);
    TEST(test_lp_series_at_zero);
    TEST(test_lp_negative_integer_s);
    TEST(test_lp_machine_real);
    TEST(test_lp_machine_complex);
    TEST(test_lp_arbitrary_precision);
    TEST(test_lp_arbitrary_complex);
    TEST(test_lp_large_z_continuation);
    TEST(test_lp_large_z_symbolic);
    TEST(test_lp_doubly_infinite);
    TEST(test_lp_include_singular_term);
    TEST(test_lp_derivatives);
    TEST(test_lp_symbolic);
    TEST(test_lp_listable);
    TEST(test_lp_attributes);
    TEST(test_lp_argcount);

    printf("All LerchPhi tests passed.\n");
    return 0;
}
