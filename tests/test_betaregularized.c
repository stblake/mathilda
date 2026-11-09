/* Tests for BetaRegularized[z, a, b] and BetaRegularized[z0, z1, a, b].
 *
 * Reference values: mpmath betainc(a, b, 0, z, regularized=True) at 40 digits.
 * (The reference system prints BetaRegularized[0.3, 200., 300.] as
 * 1.0496985243291939*10^-6; the true value is 1.04969852432929467*10^-6.)
 *
 * Covers: exact rules (the non-positive-integer parameters, a = 1, b = 1, the
 * terminating sum, unevaluated forms); machine values in [0, 1] over small and
 * large parameters, complex and out-of-range arguments, negative parameters;
 * the four-argument form; 40 digits; identities (symmetry, recurrence); D in z,
 * a, b against central differences; buffers, NDArray, Compile[]; attributes.
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

static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* {call, Re value, Im value} */
static const struct { const char* call; double re, im; } BR_REF[] = {
    { "BetaRegularized[0.3, 2.5, 3.7]", 0.3190031752843086985, 0.0 },
    { "BetaRegularized[0.3, 200., 300.]", 1.04969852432929467e-6, 0.0 },
    { "BetaRegularized[0.4, 200., 300.]", 0.50242861631993199569, 0.0 },
    { "BetaRegularized[0.7, 30., 20.]", 0.93001290669949251928, 0.0 },
    { "BetaRegularized[1.*^-5, 0.1, 5.]", 0.38685187798145923291, 0.0 },
    { "BetaRegularized[0.05, 0.01, 0.02]", 0.64752169527360224769, 0.0 },
    { "BetaRegularized[0.999, 0.5, 1000.]", 1.0, 0.0 },
    { "BetaRegularized[0.3, -1.5, 3.]", 0.81397657851462184236, 0.0 },
    { "BetaRegularized[0.3, 2., -1.5]", 0.060891806951547848217, 0.0 },
    { "BetaRegularized[-0.5, 2.5, 3.]", 0.0, 2.5798348969852944836 },
    { "BetaRegularized[2., 3.5, 4.2]", -91.218410477698483004, -67.000597077103194468 },
    { "BetaRegularized[Complex[0.3, 0.2], 2.5, 3.7]", 0.27851084586224895152, 0.43938610289759732041 },
    { "BetaRegularized[0.2, 0.5, 2., 3.]", 0.50669999999999998295, 0.0 },
};

#define BR_TOL 1.5e-14

void test_br_exact() {
    assert_true("BetaRegularized[z, a, b] === BetaRegularized[z, a, b]");
    assert_true("BetaRegularized[0, a, b] === BetaRegularized[0, a, b]");
    assert_true("BetaRegularized[z, 1, b] === 1 - (1 - z)^b");
    assert_true("BetaRegularized[z, a, 1] === z^a - 0^a");
    assert_true("BetaRegularized[1/2, 2, 3] === 11/16");
    assert_true("BetaRegularized[1/5, 1/2, 2, 3] === 5067/10000");
    assert_true("BetaRegularized[1/3, 1/2, 1/2] === BetaRegularized[1/3, 1/2, 1/2]");
    assert_true("BetaRegularized[z, -1, 3] === 1 && BetaRegularized[z, 3, -1] === 0");
    assert_true("BetaRegularized[0.3, 0, 3.] === 1. && BetaRegularized[0.3, 2., 0] === 0.");
    assert_true("BetaRegularized[0., 2., 3.] === 0. && BetaRegularized[1., 2., 3.] === 1.");
    assert_true("BetaRegularized[2., 3., 4.] == -48.");
    assert_true("NumericQ[BetaRegularized[1/3, 1/2, 1/2]]");
    test_mute_messages();
    assert_eval_eq("Check[BetaRegularized[0.5, 1.], \"argb\"]", "\"argb\"", 0);
}

void test_br_machine() {
    for (size_t i = 0; i < sizeof(BR_REF) / sizeof(BR_REF[0]); i++) {
        Expr* e = parse_expression(BR_REF[i].call);
        Expr* r = evaluate(e);
        expr_free(e);
        double re = 0, im = 0;
        if (r->type == EXPR_REAL) re = r->data.real;
        else {
            ASSERT_MSG(r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
                       r->data.function.args[0]->type == EXPR_REAL, "%s: not a number", BR_REF[i].call);
            re = r->data.function.args[0]->data.real;
            im = r->data.function.args[1]->data.real;
        }
        expr_free(r);
        ASSERT_MSG(hypot(re - BR_REF[i].re, im - BR_REF[i].im) <= BR_TOL * hypot(BR_REF[i].re, BR_REF[i].im),
                   "%s: expected %.17g + %.17g I, got %.17g + %.17g I", BR_REF[i].call,
                   BR_REF[i].re, BR_REF[i].im, re, im);
    }
    assert_true("MachineNumberQ[BetaRegularized[0.3, 2.5, 3.7]]");
    assert_true("Head[BetaRegularized[2., 3.5, 4.2]] === Complex");
}

void test_br_arbitrary() {
    assert_true("With[{v = N[BetaRegularized[3/10, 5/2, 37/10], 40]}, Abs[v - 0.3190031752843087000427729070649337445919`45] < 10^-38 && Precision[v] >= 39]");
    assert_true("Abs[BetaRegularized[0.3`30, 200, 300] - 1.04969852432929467`17*^-6] < 10^-20");
    assert_true("Precision[BetaRegularized[0.2`30, 0.5`30, 2, 3]] >= 29");
}

void test_br_identities() {
    /* symmetry I_z(a, b) = 1 - I_(1-z)(b, a) and the recurrence
     * I_z(a, b) = I_z(a + 1, b) + z^a (1-z)^b / (a Beta[a, b]) */
    assert_true("Abs[BetaRegularized[0.37, 2.2, 5.1] + BetaRegularized[0.63, 5.1, 2.2] - 1] < 10^-15");
    assert_true("Abs[BetaRegularized[0.37, 2.2, 5.1] - BetaRegularized[0.37, 3.2, 5.1] - 0.37^2.2 0.63^5.1/(2.2 Beta[2.2, 5.1])] < 10^-15");
    assert_true("Abs[BetaRegularized[0.1, 0.6, 3., 4.] - (BetaRegularized[0.6, 3., 4.] - BetaRegularized[0.1, 3., 4.])] < 10^-15");
}

void test_br_derivative() {
    assert_true("D[BetaRegularized[z, a, b], z] === (1 - z)^(b - 1) z^(a - 1)/Beta[a, b]");
    assert_true("D[BetaRegularized[z0, z1, a, b], z1] === (1 - z1)^(b - 1) z1^(a - 1)/Beta[a, b]");
    assert_true("With[{pt = {z -> 3/10, a -> 5/2, b -> 37/10}, h = 10^-20}, "
                "Abs[N[(D[BetaRegularized[z, a, b], a] /. pt) - (BetaRegularized[3/10, 5/2 + h, 37/10] - BetaRegularized[3/10, 5/2 - h, 37/10])/(2 h), 30]] < 10^-30 && "
                "Abs[N[(D[BetaRegularized[z, a, b], b] /. pt) - (BetaRegularized[3/10, 5/2, 37/10 + h] - BetaRegularized[3/10, 5/2, 37/10 - h])/(2 h), 30]] < 10^-30]");
    assert_true("With[{pt = {z0 -> 1/10, z1 -> 3/10, a -> 5/2, b -> 37/10}, h = 10^-20}, "
                "Abs[N[(D[BetaRegularized[z0, z1, a, b], a] /. pt) - (BetaRegularized[1/10, 3/10, 5/2 + h, 37/10] - BetaRegularized[1/10, 3/10, 5/2 - h, 37/10])/(2 h), 30]] < 10^-30]");
}

void test_br_buffers_compile() {
    assert_true("With[{w = Range[0.001, 0.999, 0.001]}, PackedArrayQ[BetaRegularized[w, 2.5, 3.7]] && BetaRegularized[w, 2.5, 3.7] === Table[BetaRegularized[x, 2.5, 3.7], {x, FromPackedArray[w]}]]");
    assert_true("With[{w = Range[0.3, 0.5, 0.001]}, BetaRegularized[w, 200., 300.] === Table[BetaRegularized[x, 200., 300.], {x, FromPackedArray[w]}]]");
    assert_true("With[{w = Range[0.1, 5., 0.1]}, BetaRegularized[0.4, w, 3.] === Table[BetaRegularized[0.4, x, 3.], {x, FromPackedArray[w]}]]");
    assert_true("With[{w = Range[-1., 2., 0.01]}, BetaRegularized[w, 2.5, 3.5] === Table[BetaRegularized[x, 2.5, 3.5], {x, FromPackedArray[w]}]]");
    assert_true("With[{w = Range[0.001, 0.999, 0.001]}, Normal[BetaRegularized[NDArray[w], 2.5, 3.7]] === BetaRegularized[w, 2.5, 3.7]]");
    assert_true("With[{w = Range[0.001, 0.999, 0.001]}, BetaRegularized[0.1, w, 2.5, 3.7] === Table[BetaRegularized[0.1, x, 2.5, 3.7], {x, FromPackedArray[w]}]]");
    assert_true("With[{f = Compile[{{x, _Real, 1}}, BetaRegularized[x, 2.5, 3.7]], w = Range[0.001, 0.999, 0.001]}, f[w] === BetaRegularized[w, 2.5, 3.7]]");
    assert_true("With[{f = Compile[{{x, _Real, 1}}, BetaRegularized[0.1, x, 2.5, 3.7]], w = Range[0.001, 0.999, 0.001]}, f[w] === BetaRegularized[0.1, w, 2.5, 3.7]]");
    assert_true("With[{f = Compile[{{x, _Real}}, BetaRegularized[x, 2.5, 3.7]]}, f[0.3] === BetaRegularized[0.3, 2.5, 3.7]]");
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, BetaRegularized[x, 2.5, 3.7]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, BetaRegularized[v, 2.5, 3.7]][[1]]", "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, BetaRegularized[0.1, x, 2.5, 3.7]][[1]]", "\"Compiled\" -> True", 0);
}

void test_br_attributes() {
    assert_eval_eq("Attributes[BetaRegularized]", "{Listable, NumericFunction, Protected}", 0);
    SymbolDef* d = symtab_get_def("BetaRegularized");
    ASSERT(d != NULL && d->docstring != NULL && strstr(d->docstring, "BetaRegularized[z, a, b]"));
}

int main() {
    symtab_init();
    core_init();
    TEST(test_br_exact);
    TEST(test_br_machine);
    TEST(test_br_arbitrary);
    TEST(test_br_identities);
    TEST(test_br_derivative);
    TEST(test_br_buffers_compile);
    TEST(test_br_attributes);
    printf("All BetaRegularized tests passed.\n");
    return 0;
}
