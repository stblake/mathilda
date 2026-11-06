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
    assert_eval_eq("PolyLog[2, 1/2]", "1/12 Pi^2 - 1/2 Log[2]^2", 0);
    assert_eval_eq("PolyLog[3, 1/2]",
                   "-1/12 Pi^2 Log[2] + 1/6 Log[2]^3 + 7/8 Zeta[3]", 0);
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
    /* Nielsen: d/dz S_{n,p}(z) = S_{n-1,p}(z)/z; the orders stay inert. */
    assert_eval_eq("D[PolyLog[n, p, z], z]", "PolyLog[-1 + n, p, z]/z", 0);
    assert_eval_eq("D[PolyLog[n, p, z], n]", "Derivative[1, 0, 0][PolyLog][n, p, z]", 0);
    assert_eval_eq("D[PolyLog[n, p, z], p]", "Derivative[0, 1, 0][PolyLog][n, p, z]", 0);
    assert_eval_eq("D[PolyLog[n, p, y], z]", "0", 0);
}

/* ---- Nielsen machine values: a real value carries no imaginary dust --- */

void test_polylog_nielsen_real_value() {
    /* PolyLog[0, 2, 2.] = Log(-1)^2/2 = -Pi^2/2, exactly real: an integer
     * power of the complex logarithm used to go through cpow and leave
     * -6e-16 I behind. The value must be a Real, on every surface. */
    assert_eval_eq("Head[PolyLog[0, 2, 2.]]", "Real", 0);
    assert_eval_eq("Abs[PolyLog[0, 2, 2.] + Pi^2/2] < 10^-15", "True", 0);
    assert_eval_eq("With[{v = {2., 3., 0.5, -1.5}}, "
                   "(PolyLog[0, 2, #] & /@ v) === PolyLog[0, 2, v] === "
                   "Normal[PolyLog[0, 2, NDArray[v]]] === "
                   "Compile[{{x, _Real, 1}}, PolyLog[0, 2, x]][v]]", "True", 0);
}

/* ---- diagnostics & attributes --------------------------------------- */

void test_polylog_argcount() {
    /* Wrong argument counts emit PolyLog::argt and stay unevaluated. */
    assert_eval_eq("PolyLog[]",  "PolyLog[]",  0);
    assert_eval_eq("PolyLog[2]", "PolyLog[2]", 0);
    assert_eval_eq("PolyLog[1, 2, 3, 4]", "PolyLog[1, 2, 3, 4]", 0);
}

/* ---- FunctionExpand at roots of unity -------------------------------- */

static void check30_eq(const char* a, const char* b) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Abs[N[(%s) - (%s), 40]] < 10^-30 Max[1, Abs[N[%s, 40]]]", a, b, a);
    assert_eval_eq(buf, "True", 0);
}

void test_polylog_roots_of_unity() {
    /* closed forms (the reference system's values) */
    assert_eval_eq("FunctionExpand[PolyLog[3, I] + PolyLog[3, -I]]", "-3/16 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[3, I]]", "(1/32*I) Pi^3 - 3/32 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[2, I] - PolyLog[2, -I]]", "(2*I) Catalan", 0);
    assert_eval_eq("FunctionExpand[PolyLog[2, E^(I Pi/3)] + PolyLog[2, E^(-I Pi/3)]]", "1/18 Pi^2", 0);
    assert_eval_eq("FunctionExpand[PolyLog[3, E^(2 Pi I/3)]]", "(2/81*I) Pi^3 - 4/9 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[5, -I]]", "(-5/1536*I) Pi^5 - 15/512 Zeta[5]", 0);
    assert_eval_eq("FunctionExpand[Re[PolyLog[3, E^(I Pi/3)]]]", "1/3 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[Re[PolyLog[5, E^(I Pi/3)]]]", "25/54 Zeta[5]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[3, I] + PolyLog[3, -I] - Zeta[3]/4]", "-7/16 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[3, (-1)^(2/3)]]", "(2/81*I) Pi^3 - 4/9 Zeta[3]", 0);
    /* a non-root and a symbolic argument are left alone */
    assert_eval_eq("FunctionExpand[PolyLog[3, 2 I]]", "PolyLog[3, 2*I]", 0);
    assert_eval_eq("FunctionExpand[PolyLog[3, z]]", "PolyLog[3, z]", 0);

    /* every (order, root) pair: the expansion against the MPFR polylog */
    const char* roots[] = { "I", "-I", "(-1)^(1/3)", "E^(I Pi/3)", "E^(2 I Pi/3)", "E^(-I Pi/4)",
                            "E^(3 I Pi/4)", "E^(I Pi/6)", "E^(-5 I Pi/6)", "E^(2 I Pi/5)",
                            "E^(I Pi/8)", "-I E^(-I Pi/4)", "E^(7 I Pi/12)" };
    for (int n = 2; n <= 7; n++)
        for (size_t k = 0; k < sizeof roots / sizeof roots[0]; k++) {
            char a[256], b[256];
            snprintf(a, sizeof a, "FunctionExpand[PolyLog[%d, %s]]", n, roots[k]);
            snprintf(b, sizeof b, "PolyLog[%d, N[%s, 50]]", n, roots[k]);
            check30_eq(a, b);
        }

    /* Log at roots of unity, and 1 + root */
    assert_eval_eq("FunctionExpand[Log[E^(I Pi/4)]]", "(1/4*I) Pi", 0);
    assert_eval_eq("FunctionExpand[Log[1 - E^(I Pi/4)]]", "(-3/8*I) Pi + Log[2 Cos[3/8 Pi]]", 0);
    const char* logs[] = { "Log[1 - E^(I Pi/4)]", "Log[1 + I]", "Log[1 - I E^(-I Pi/4)]",
                           "Log[1 + E^(5 I Pi/6)]", "Log[E^(-2 I Pi/3)]" };
    for (size_t k = 0; k < sizeof logs / sizeof logs[0]; k++) {
        char a[256];
        snprintf(a, sizeof a, "FunctionExpand[%s]", logs[k]);
        check30_eq(a, logs[k]);
    }

    /* PolyGamma at 1/2, 1/4, 3/4 */
    assert_eval_eq("FunctionExpand[PolyGamma[1, 1/4]]", "8 Catalan + Pi^2", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[2, 1/4]]", "-2 Pi^3 - 56 Zeta[3]", 0);
    assert_eval_eq("FunctionExpand[PolyGamma[3, 1/4]]", "PolyGamma[3, 1/4]", 0);
    const char* pg[] = { "PolyGamma[1, 1/4]", "PolyGamma[1, 3/4]", "PolyGamma[2, 3/4]",
                         "PolyGamma[4, 1/4]", "PolyGamma[3, 1/2]", "PolyGamma[2, 1/2]" };
    for (size_t k = 0; k < sizeof pg / sizeof pg[0]; k++) {
        char a[256];
        snprintf(a, sizeof a, "FunctionExpand[%s]", pg[k]);
        check30_eq(a, pg[k]);
    }
}

void test_polylog_attributes() {
    SymbolDef* d = symtab_get_def("PolyLog");
    ASSERT(d != NULL);
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "PolyLog must be Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0,
               "PolyLog must be a NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "PolyLog must be Protected");
}

/* ---- real integer order just below z = 1 (the Jonquiere kernel) ------ *
 * Checked against identities whose other terms the kernel does not compute
 * (arguments below 1/2, or negative): Euler's reflection for Li_2, the
 * three-term identity for Li_3, and the defining series for Li_5. */
static void check60_zero(const char* expr) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Abs[%s] < 10^-60", expr);
    assert_eval_eq(buf, "True", 0);
}

void test_polylog_near_one_real() {
    /* Li_2(z) + Li_2(1 - z) = Pi^2/6 - Log[z] Log[1 - z] */
    const char* zs[] = { "3/4", "9/10", "1 - 10^-5", "1 - 10^-30", "1 - 10^-60" };
    for (size_t k = 0; k < sizeof zs / sizeof zs[0]; k++) {
        char e[1024];
        snprintf(e, sizeof e,
                 "With[{z = N[%s, 80]}, PolyLog[2, z] + PolyLog[2, 1 - z] - Pi^2/6 + Log[z] Log[1 - z]]",
                 zs[k]);
        check60_zero(e);
    }
    /* Li_3(z) + Li_3(1-z) + Li_3(1-1/z) = Zeta[3] + Log[z]^3/6 + Pi^2/6 Log[z]
     *                                     - Log[z]^2 Log[1-z]/2 */
    const char* z3[] = { "3/4", "2/3", "99/100", "1 - 10^-30" };
    for (size_t k = 0; k < sizeof z3 / sizeof z3[0]; k++) {
        char e[1024];
        snprintf(e, sizeof e,
                 "With[{z = N[%s, 80]}, PolyLog[3, z] + PolyLog[3, 1 - z] + PolyLog[3, 1 - 1/z]"
                 " - Zeta[3] - Log[z]^3/6 - Pi^2/6 Log[z] + Log[z]^2 Log[1 - z]/2]", z3[k]);
        check60_zero(e);
    }
    /* a z below 1 that rounds to the double 1.0 stays real */
    assert_eval_eq("Head[PolyLog[3, N[1 - 10^-30, 60]]]", "Real", 0);
    /* higher order, against the defining series at z = 3/4 */
    check60_zero("With[{z = N[3/4, 80]}, PolyLog[5, z] - Sum[z^k/k^5, {k, 1, 700}]]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_polylog_low_order);
    TEST(test_polylog_near_one_real);
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
    TEST(test_polylog_nielsen_real_value);
    TEST(test_polylog_argcount);
    TEST(test_polylog_attributes);
    TEST(test_polylog_roots_of_unity);

    printf("All PolyLog tests passed.\n");
    return 0;
}
