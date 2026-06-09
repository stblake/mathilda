/* Tests for PolyLog: the polylogarithm PolyLog[n, z] (and the accepted-but-
 * symbolic Nielsen form PolyLog[n, p, z]).
 *
 * Covers the exact closed forms (order <= 1 logarithmic/rational, negative
 * integer Eulerian rationals, PolyLog[n, +-1] in terms of Zeta, and the famous
 * PolyLog[2, 1/2] / PolyLog[3, 1/2] values), Listable threading (including the
 * empty list), machine and arbitrary-precision (MPFR) real and complex
 * numerics, symbolic fall-through, the derivative rules, argument-count
 * diagnostics, and attributes. */

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

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- numeric helpers ------------------------------------------------- */

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
               "%s: expected %.10g, got %.10g", input, expected, v);
}

/* ---- exact closed forms: order <= 1 --------------------------------- */

void test_polylog_low_order() {
    /* PolyLog[1, z] = -Log[1 - z];  PolyLog[0, z] = z/(1 - z). */
    assert_eval_eq("PolyLog[1, z]", "-Log[1 - z]", 0);
    assert_eval_eq("PolyLog[0, z]", "z/(1 - z)", 0);
    /* PolyLog[n, 0] = 0 for any n. */
    assert_eval_eq("PolyLog[2, 0]", "0", 0);
    assert_eval_eq("PolyLog[n, 0]", "0", 0);
    assert_eval_eq("PolyLog[-3, 0]", "0", 0);
}

void test_polylog_negative_order() {
    /* Negative integer orders are Eulerian-number rational functions. */
    assert_eval_eq("PolyLog[-1, z]", "z/(1 - z)^2", 0);
    assert_eval_eq("PolyLog[-2, z]", "(z + z^2)/(1 - z)^3", 0);
    assert_eval_eq("PolyLog[-3, z]", "(z + 4 z^2 + z^3)/(1 - z)^4", 0);
    assert_eval_eq("PolyLog[-4, z]", "(z + 11 z^2 + 11 z^3 + z^4)/(1 - z)^5", 0);
    /* Exact numeric argument reduces to an exact rational: Li_{-1}(2) = 2. */
    assert_eval_eq("PolyLog[-1, 2]", "2", 0);
}

/* ---- exact closed forms: order >= 2 special arguments --------------- */

void test_polylog_at_one() {
    /* PolyLog[n, 1] = Zeta[n] for integer n >= 2. */
    assert_eval_eq("PolyLog[2, 1]", "1/6 Pi^2", 0);
    assert_eval_eq("PolyLog[3, 1]", "Zeta[3]", 0);
    assert_eval_eq("PolyLog[4, 1]", "1/90 Pi^4", 0);
    assert_eval_eq("PolyLog[6, 1]", "1/945 Pi^6", 0);
}

void test_polylog_at_minus_one() {
    /* PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n] = -eta(n). */
    assert_eval_eq("PolyLog[2, -1]", "-1/12 Pi^2", 0);
    assert_eval_eq("PolyLog[4, -1]", "-7/720 Pi^4", 0);
}

void test_polylog_half() {
    /* The famous PolyLog[2, 1/2] and PolyLog[3, 1/2] closed forms. */
    assert_eval_eq("PolyLog[2, 1/2]", "-1/2 Log[2]^2 + 1/12 Pi^2", 0);
    assert_eval_eq("PolyLog[3, 1/2]",
                   "1/6 Log[2]^3 - 1/12 Log[2] Pi^2 + 7/8 Zeta[3]", 0);
    /* Numerical agreement with the closed form. */
    assert_close("N[PolyLog[2, 1/2]]", 0.5822405264650125, 1e-12);
}

/* ---- Listable threading --------------------------------------------- */

void test_polylog_listable() {
    assert_eval_eq("PolyLog[{1, 0}, z]", "{-Log[1 - z], z/(1 - z)}", 0);
    assert_eval_eq("PolyLog[{2, 4}, -1]", "{-1/12 Pi^2, -7/720 Pi^4}", 0);
}

void test_polylog_empty_list() {
    /* Listable threading over an empty list yields an empty list. */
    assert_eval_eq("PolyLog[{}, x]", "{}", 0);
    assert_eval_eq("PolyLog[2, {}]", "{}", 0);
}

/* ---- numeric evaluation: machine precision -------------------------- */

void test_polylog_machine_real() {
    /* Direct real series (|z| < 1). */
    assert_close("PolyLog[2, 0.9]", 1.2997147230049587, 1e-12);
    assert_close("PolyLog[2, 0.3]", 0.3261295100754761, 1e-12);
    assert_close("PolyLog[2, -0.5]", -0.4484142069236462, 1e-12);
    /* Order-0 closed form, evaluated numerically: 5/(1-5) = -1.25. */
    assert_close("PolyLog[0, 5.0]", -1.25, 1e-12);
    /* Inexact order is treated as the integer it rounds to. */
    assert_close("PolyLog[2.0, 0.5]", 0.5822405264650125, 1e-12);
}

void test_polylog_branch_cut() {
    /* Real z > 1 sits on the branch cut [1, Infinity); the value is taken
     * continuous from below (Mathematica's convention): for the dilogarithm
     * Im Li_2(x) = -Pi Log[x] and Re Li_2(x) = Pi^2/3 - Log[x]^2/2. */
    assert_close("Re[PolyLog[2, 2.0]]", 2.4674011002723395, 1e-9);
    assert_close("Im[PolyLog[2, 2.0]]", -M_PI * log(2.0), 1e-9);
    assert_close("Im[PolyLog[2, 5.0]]", -M_PI * log(5.0), 1e-8);
    /* Trilogarithm on the cut: Im Li_3(x) = -Pi Log[x]^2/2. */
    assert_close("Im[PolyLog[3, 2.0]]", -M_PI * log(2.0) * log(2.0) / 2.0, 1e-9);
    /* The branch point z = 1 is finite and real: Li_2(1) = Pi^2/6. */
    assert_close("PolyLog[2, 1.0]", M_PI * M_PI / 6.0, 1e-9);
    assert_close("PolyLog[3, 1.0]", 1.2020569031595942, 1e-9); /* Zeta[3] */
}

void test_polylog_machine_complex() {
    /* Complex argument, |z| > 1, complex order: zeta expansion + reflection. */
    assert_close("Re[PolyLog[.2 + I, .5 - I]]", -0.08985259, 1e-6);
    assert_close("Im[PolyLog[.2 + I, .5 - I]]", -0.59586506, 1e-6);
    /* Integer order, complex argument with |z| > 1/2. */
    assert_close("Re[PolyLog[2, 0.5 + 0.3 I]]", 0.5309619, 1e-6);
    assert_close("Im[PolyLog[2, 0.5 + 0.3 I]]", 0.4034368, 1e-6);
}

/* ---- numeric evaluation: arbitrary precision ------------------------ */

void test_polylog_arbitrary_precision() {
    /* PolyLog[1, 1/3] = -Log[2/3]; 50-digit numericalisation. */
    assert_eval_startswith("N[PolyLog[1, 1/3], 50]",
                           "0.40546510810816438197801311546434913657199042346249");
    /* Direct real series at 100-digit precision. */
    assert_eval_startswith(
        "PolyLog[3.0`100, 0.5`100]",
        "0.5372131936080402009406232255949658266704024993403781706897619307");
    /* Non-integer order, high precision (the prompt's timing example value). */
    assert_eval_startswith(
        "PolyLog[2/3, 1/3`100]",
        "0.428036823753647804105765196060585302286899263746351143200978600558");
}

/* ---- symbolic fall-through ------------------------------------------ */

void test_polylog_symbolic() {
    /* No closed form / no inexact operand -> stays unevaluated. */
    assert_eval_eq("PolyLog[2, 1/3]", "PolyLog[2, 1/3]", 0);
    assert_eval_eq("PolyLog[5, 1/2]", "PolyLog[5, 1/2]", 0);
    assert_eval_eq("PolyLog[2, x]",   "PolyLog[2, x]",   0);
    assert_eval_eq("PolyLog[n, z]",   "PolyLog[n, z]",   0);
    /* Nielsen 3-argument form is accepted but left symbolic. */
    assert_eval_eq("PolyLog[2, 3, x]", "PolyLog[2, 3, x]", 0);
}

/* ---- derivatives ---------------------------------------------------- */

void test_polylog_derivative() {
    /* d/dz PolyLog[n, z] = PolyLog[n-1, z]/z. */
    assert_eval_eq("D[PolyLog[n, z], z]", "PolyLog[-1 + n, z]/z", 0);
    /* PolyLog[1, x] folds, giving -Log[1-x]/x. */
    assert_eval_eq("D[PolyLog[2, x], x]", "-Log[1 - x]/x", 0);
    /* Chain rule through an inner function. */
    assert_eval_eq("D[PolyLog[3, x^2], x]", "(2 PolyLog[2, x^2])/x", 0);
    /* Order derivative has no elementary form. */
    assert_eval_eq("D[PolyLog[s, z], s]", "Derivative[1, 0][PolyLog][s, z]", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_polylog_argcount() {
    /* Wrong argument counts emit PolyLog::argt and stay unevaluated. */
    assert_eval_eq("PolyLog[]",  "PolyLog[]",  0);
    assert_eval_eq("PolyLog[2]", "PolyLog[2]", 0);
    assert_eval_eq("PolyLog[1, 2, 3, 4]", "PolyLog[1, 2, 3, 4]", 0);
}

void test_polylog_attributes() {
    SymbolDef* d = symtab_get_def("PolyLog");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "PolyLog must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "PolyLog must be a NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "PolyLog must be Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_polylog_low_order);
    TEST(test_polylog_negative_order);
    TEST(test_polylog_at_one);
    TEST(test_polylog_at_minus_one);
    TEST(test_polylog_half);
    TEST(test_polylog_listable);
    TEST(test_polylog_empty_list);
    TEST(test_polylog_machine_real);
    TEST(test_polylog_branch_cut);
    TEST(test_polylog_machine_complex);
    TEST(test_polylog_arbitrary_precision);
    TEST(test_polylog_symbolic);
    TEST(test_polylog_derivative);
    TEST(test_polylog_argcount);
    TEST(test_polylog_attributes);

    printf("All PolyLog tests passed.\n");
    return 0;
}
