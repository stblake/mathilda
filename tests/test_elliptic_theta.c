/* Tests for EllipticTheta[a, u, q] and EllipticThetaPrime[a, u, q] (nome q).
 *
 * Reference values: mpmath jtheta at 500 digits (at 130 digits mpmath itself is
 * wrong for |q| > 0.99; the reference system is wrong there in machine
 * precision too: EllipticThetaPrime[4, 0.7, 0.99] is 5.23882*10^-30, it prints
 * 6.7*10^-16).
 *
 * Covers: exact rules (q = 0, the ConditionalExpression zeros, the two-argument
 * form, the etype message); 64 mpmath points (all four functions and their
 * derivatives; real and complex u and q, q < 0, |q| up to 0.999, large Re u);
 * 40 digits; identities (Jacobi's theta_2^4 + theta_4^4 = theta_3^4,
 * theta_1' = theta_2 theta_3 theta_4 at u = 0, the derivative); D; buffers and
 * Compile[]; attributes.
 */

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

/* {a, prime, Re u, Im u, Re q, Im q, Re value, Im value} */
static const double TH_REF[][8] = {
    { 1, 0, 0.7, 0.0, 0.3, 0.0, 8.3817877516948841156e-1, 0.0 },
    { 1, 0, 0.7, 0.0, -0.3, 0.0, 5.9268189576897985563e-1, 5.9268189576897985563e-1 },
    { 1, 0, 0.2, 0.0, 0.95, 0.0, 9.628674993302892376e-16, 0.0 },
    { 1, 0, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, 1.3418652789216800557e-38, 0.0 },
    { 1, 0, 1.0, 2.0, 0.1, 0.2, -5.0822781712892283433, -8.9905075780131412153 },
    { 1, 0, -5.5, 0.3, 0.6, -0.2, 1.0610999119167477497, 1.4198153112826004472 },
    { 1, 0, 12.3, 0.0, 0.5, 0.0, -1.664900470088892959e-1, 0.0 },
    { 1, 0, 0.4, 0.0, 0.999, 0.0, -1.3205565549165733917e-504, 0.0 },
    { 1, 1, 0.7, 0.0, 0.3, 0.0, 1.3287987781514223629, 0.0 },
    { 1, 1, 0.7, 0.0, -0.3, 0.0, 9.3960262686326952012e-1, 9.3960262686326952012e-1 },
    { 1, 1, 0.2, 0.0, 0.95, 0.0, 5.1464630887963758018e-14, 0.0 },
    { 1, 1, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, -5.1346207883035655171e-36, 0.0 },
    { 1, 1, 1.0, 2.0, 0.1, 0.2, -3.5826608531121697913e1, 2.5873127636212258237e1 },
    { 1, 1, -5.5, 0.3, 0.6, -0.2, 5.329059972207752557, -2.5071690231140959204e-1 },
    { 1, 1, 12.3, 0.0, 0.5, 0.0, 7.7480845587370419593e-1, 0.0 },
    { 1, 1, 0.4, 0.0, 0.999, 0.0, -6.1145499407382972142e-555, 0.0 },
    { 2, 0, 0.7, 0.0, 0.3, 0.0, 1.0638297480113753269, 0.0 },
    { 2, 0, 0.7, 0.0, -0.3, 0.0, 7.5224122884681956463e-1, 7.5224122884681956463e-1 },
    { 2, 0, 0.2, 0.0, 0.95, 0.0, 3.5881384757923322245, 0.0 },
    { 2, 0, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, -1.1467619909505779629e-209, 0.0 },
    { 2, 0, 1.0, 2.0, 0.1, 0.2, 1.5314645573110735506e1, -8.2492102893370679754 },
    { 2, 0, -5.5, 0.3, 0.6, -0.2, 6.0942605744028123042e-1, -3.0249763618128106977e-1 },
    { 2, 0, 12.3, 0.0, 0.5, 0.0, 1.9217765406393842817, 0.0 },
    { 2, 0, 0.4, 0.0, 0.999, 0.0, 1.9774183727120943408e-68, 0.0 },
    { 2, 1, 0.7, 0.0, 0.3, 0.0, -1.2966279695092942712, 0.0 },
    { 2, 1, 0.7, 0.0, -0.3, 0.0, -9.1685442991616595945e-1, -9.1685442991616595945e-1 },
    { 2, 1, 0.2, 0.0, 0.95, 0.0, -2.7981345465408135035e1, 0.0 },
    { 2, 1, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, -1.0139394404675525e-206, 0.0 },
    { 2, 1, 1.0, 2.0, 0.1, 0.2, -1.4805396900647858568e1, -3.6972461201055683763e1 },
    { 2, 1, -5.5, 0.3, 0.6, -0.2, -1.5834157365370173568, 1.3066224375060785647 },
    { 2, 1, 12.3, 0.0, 0.5, 0.0, 1.4771716068429048517, 0.0 },
    { 2, 1, 0.4, 0.0, 0.999, 0.0, -1.5811435989267420825e-65, 0.0 },
    { 3, 0, 0.7, 0.0, 0.3, 0.0, 1.0866969908910608074, 0.0 },
    { 3, 0, 0.7, 0.0, -0.3, 0.0, 8.8277501862550187038e-1, 0.0 },
    { 3, 0, 0.2, 0.0, 0.95, 0.0, 3.5881384757923322245, 0.0 },
    { 3, 0, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, 1.1467619909505779629e-209, 0.0 },
    { 3, 0, 1.0, 2.0, 0.1, 0.2, 1.5276292346400704274e1, -6.5782876212174785236 },
    { 3, 0, -5.5, 0.3, 0.6, -0.2, 6.1253130611063327944e-1, -2.964466767234686738e-1 },
    { 3, 0, 12.3, 0.0, 0.5, 0.0, 1.9218049188744731206, 0.0 },
    { 3, 0, 0.4, 0.0, 0.999, 0.0, 1.9774183727120943408e-68, 0.0 },
    { 3, 1, 0.7, 0.0, 0.3, 0.0, -1.2040410019278929798, 0.0 },
    { 3, 1, 0.7, 0.0, -0.3, 0.0, 1.1606266246243559595, 0.0 },
    { 3, 1, 0.2, 0.0, 0.95, 0.0, -2.7981345465408135035e1, 0.0 },
    { 3, 1, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, 1.0139394404675525e-206, 0.0 },
    { 3, 1, 1.0, 2.0, 0.1, 0.2, -7.6303275525638157473, -4.1504750268962926599e1 },
    { 3, 1, -5.5, 0.3, 0.6, -0.2, -1.5293113541569350405, 1.3268561700739621174 },
    { 3, 1, 12.3, 0.0, 0.5, 0.0, 1.4769402569120903634, 0.0 },
    { 3, 1, 0.4, 0.0, 0.999, 0.0, -1.5811435989267420825e-65, 0.0 },
    { 4, 0, 0.7, 0.0, 0.3, 0.0, 8.8277501862550187038e-1, 0.0 },
    { 4, 0, 0.7, 0.0, -0.3, 0.0, 1.0866969908910608074, 0.0 },
    { 4, 0, 0.2, 0.0, 0.95, 0.0, 9.6286749937442530557e-16, 0.0 },
    { 4, 0, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, 1.3418652789216800557e-38, 0.0 },
    { 4, 0, 1.0, 2.0, 0.1, 0.2, 2.7970219138731325331e-1, 1.277241978519190867e1 },
    { 4, 0, -5.5, 0.3, 0.6, -0.2, 1.060716456685581388, 1.4197830962638043215 },
    { 4, 0, 12.3, 0.0, 0.5, 0.0, 1.9918317805198452266e-1, 0.0 },
    { 4, 0, 0.4, 0.0, 0.999, 0.0, -1.2232156860700653449e-504, 0.0 },
    { 4, 1, 0.7, 0.0, 0.3, 0.0, 1.1606266246243559595, 0.0 },
    { 4, 1, 0.7, 0.0, -0.3, 0.0, -1.2040410019278929798, 0.0 },
    { 4, 1, 0.2, 0.0, 0.95, 0.0, 5.146463088491634264e-14, 0.0 },
    { 4, 1, 2.045260087875171, 0.0, 0.9975231757436976, 0.0, -5.1346207883035655171e-36, 0.0 },
    { 4, 1, 1.0, 2.0, 0.1, 0.2, 3.2406051941863748581e1, -1.2716573902520487504e1 },
    { 4, 1, -5.5, 0.3, 0.6, -0.2, 5.3320626894776750359, -2.5199243634345679384e-1 },
    { 4, 1, 12.3, 0.0, 0.5, 0.0, -6.0150402290854109125e-1, 0.0 },
    { 4, 1, 0.4, 0.0, 0.999, 0.0, -8.1199657557525649098e-554, 0.0 },
};
static const char* const TH_MP = "1.31510723839933155413074436959379114712554717`45";

static Expr* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static void eval_num(const char* input, double* re, double* im) {
    Expr* r = eval_str(input);
    if (r->type == EXPR_REAL) { *re = r->data.real; *im = 0.0; expr_free(r); return; }
    ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
               r->data.function.head->type == EXPR_SYMBOL &&
               strcmp(r->data.function.head->data.symbol.name, "Complex") == 0 &&
               r->data.function.args[0]->type == EXPR_REAL &&
               r->data.function.args[1]->type == EXPR_REAL,
               "%s: expected a machine number", input);
    *re = r->data.function.args[0]->data.real;
    *im = r->data.function.args[1]->data.real;
    expr_free(r);
}

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

static const char* R(double v, char* buf) {
    char t[40];
    snprintf(t, sizeof t, "%.17g", v);
    char* e = strchr(t, 'e');
    if (e) { *e = 0; snprintf(buf, 64, "%s%s*^%s", t, strchr(t, '.') ? "" : ".", e + 1); }
    else snprintf(buf, 64, "%s%s", t, strchr(t, '.') ? "" : ".");
    return buf;
}

#define TH_TOL 2.9e-14     /* 128 ulp */

void test_th_exact() {
    assert_eval_eq("EllipticTheta[1, u, 0]", "0", 0);
    assert_eval_eq("EllipticTheta[3, u, 0]", "1", 0);
    assert_eval_eq("EllipticThetaPrime[3, u, 0]", "0", 0);
    assert_true("EllipticTheta[1, 0, q] === ConditionalExpression[0, Abs[q] < 1]");
    assert_true("EllipticThetaPrime[3, 0, q] === ConditionalExpression[0, Abs[q] < 1]");
    assert_true("EllipticTheta[1, q] === ConditionalExpression[0, Abs[q] < 1]");
    assert_eval_eq("EllipticTheta[3, q]", "EllipticTheta[3, q]", 0);
    assert_eval_eq("EllipticTheta[2, 0, q]", "EllipticTheta[2, 0, q]", 0);
    assert_true("EllipticTheta[2, 0.3] == EllipticTheta[2, 0, 0.3]");
    assert_true("NumericQ[EllipticTheta[3, 1, 1/2]]");
    test_mute_messages();
    assert_eval_eq("Check[EllipticTheta[5, u, q], \"failed\"]", "\"failed\"", 0);
}

void test_th_machine() {
    char buf[400], b0[64], b1[64], b2[64], b3[64];
    for (size_t i = 0; i < sizeof(TH_REF) / sizeof(TH_REF[0]); i++) {
        const double* c = TH_REF[i];
        char us[160], qs[160];
        if (c[3] == 0.0) snprintf(us, sizeof us, "%s", R(c[2], b0));
        else snprintf(us, sizeof us, "Complex[%s, %s]", R(c[2], b0), R(c[3], b1));
        if (c[5] == 0.0) snprintf(qs, sizeof qs, "%s", R(c[4], b2));
        else snprintf(qs, sizeof qs, "Complex[%s, %s]", R(c[4], b2), R(c[5], b3));
        snprintf(buf, sizeof buf, "%s[%d, %s, %s]", c[1] != 0.0 ? "EllipticThetaPrime" : "EllipticTheta",
                 (int)c[0], us, qs);
        double re, im;
        eval_num(buf, &re, &im);
        ASSERT_MSG(hypot(re - c[6], im - c[7]) <= TH_TOL * hypot(c[6], c[7]),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", buf, c[6], c[7], re, im);
    }
    assert_true("MachineNumberQ[EllipticTheta[3, 0.7, 0.3]]");
    assert_true("Head[EllipticTheta[2, 0.7, -0.3]] === Complex");
}

void test_th_arbitrary() {
    char buf[300];
    snprintf(buf, sizeof buf, "With[{v = N[EllipticTheta[3, 1/3, 1/5], 40], r = %s}, Abs[v - r] < 10^-38 Abs[r] && Precision[v] >= 39]", TH_MP);
    assert_true(buf);
    /* a value far below 2^-53: relative accuracy, not absolute */
    assert_true("Abs[EllipticThetaPrime[3, SetPrecision[2.045260087875171, 30], SetPrecision[0.9975231757436976, 30]]/1.0139394404675525`16*^-206 - 1] < 10^-14");
}

void test_th_identities() {
    assert_true("With[{q = 0.37}, Abs[EllipticTheta[2, 0, q]^4 + EllipticTheta[4, 0, q]^4 - EllipticTheta[3, 0, q]^4] < 10^-14]");
    assert_true("With[{q = 0.61`30}, Abs[EllipticThetaPrime[1, 0, q] - EllipticTheta[2, 0, q] EllipticTheta[3, 0, q] EllipticTheta[4, 0, q]] < 10^-27]");
    assert_true("Abs[(D[EllipticTheta[4, u, 0.4], u] /. u -> 0.9) - (EllipticTheta[4, 0.9 + 10.^-6, 0.4] - EllipticTheta[4, 0.9 - 10.^-6, 0.4])/(2 10.^-6)] < 10^-8");
    /* the period: theta_1(u + pi) = -theta_1(u), theta_3(u + pi) = theta_3(u) */
    assert_true("Abs[EllipticTheta[1, 0.4 + Pi, 0.3] + EllipticTheta[1, 0.4, 0.3]] < 10^-15 && Abs[EllipticTheta[3, 0.4 + 7 Pi, 0.3] - EllipticTheta[3, 0.4, 0.3]] < 10^-14");
}

void test_th_derivative() {
    assert_true("D[EllipticTheta[3, u, q], u] === EllipticThetaPrime[3, u, q]");
    assert_true("D[EllipticThetaPrime[3, u, q], u] === Derivative[0, 1, 0][EllipticThetaPrime][3, u, q]");
    assert_true("D[EllipticTheta[3, u, q], q] === Derivative[0, 0, 1][EllipticTheta][3, u, q]");
    assert_true("D[EllipticTheta[1, x^2, q], x] === 2 x EllipticThetaPrime[1, x^2, q]");
}

void test_th_buffers_compile() {
    for (int a = 1; a <= 4; a++) {
        char buf[300];
        snprintf(buf, sizeof buf, "With[{v = Range[-20., 20., 0.05]}, PackedArrayQ[EllipticTheta[%d, v, 0.45]] && EllipticTheta[%d, v, 0.45] === Table[EllipticTheta[%d, x, 0.45], {x, FromPackedArray[v]}]]", a, a, a);
        assert_true(buf);
        snprintf(buf, sizeof buf, "With[{v = Range[-3., 3., 0.05]}, EllipticThetaPrime[%d, v, 0.97] === Table[EllipticThetaPrime[%d, x, 0.97], {x, FromPackedArray[v]}]]", a, a);
        assert_true(buf);
    }
    assert_true("With[{v = Range[-0.9, 0.9, 0.05]}, EllipticTheta[2, 0.5, v] === Table[EllipticTheta[2, 0.5, q], {q, FromPackedArray[v]}]]");
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, EllipticTheta[3, x, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, EllipticThetaPrime[1, v, 0.3]][[1]]", "\"Compiled\" -> True", 0);
    assert_true("With[{f = Compile[{{v, _Real, 1}}, EllipticTheta[4, v, 0.6]], v = Range[-5., 5., 0.01]}, f[v] === EllipticTheta[4, v, 0.6]]");
}

void test_th_attributes() {
    assert_eval_eq("Attributes[EllipticTheta]", "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    assert_eval_eq("Attributes[EllipticThetaPrime]", "{Listable, NHoldFirst, NumericFunction, Protected}", 0);
    SymbolDef* d = symtab_get_def("EllipticTheta");
    ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, "EllipticTheta[a, u, q]"));
}

int main() {
    symtab_init();
    core_init();
    TEST(test_th_exact);
    TEST(test_th_machine);
    TEST(test_th_arbitrary);
    TEST(test_th_identities);
    TEST(test_th_derivative);
    TEST(test_th_buffers_compile);
    TEST(test_th_attributes);
    printf("All EllipticTheta tests passed.\n");
    return 0;
}
