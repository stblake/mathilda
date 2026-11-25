/* Tests for the hypergeometric-family machine kernels (sf_machine_hyper.c) and
 * their wiring: Hypergeometric0F1 / 1F1 / 2F1 / PFQ, HypergeometricU,
 * WhittakerM / WhittakerW.
 *
 * Covers:
 *   - spot values against mpmath (30 digits, at the double parameters -- what a
 *     machine kernel actually receives), real and complex, across every regime:
 *     the plain series, the double-double re-sum (0F1 at negative z), Miller's
 *     recurrence (0F1 at z << 0), the I_nu asymptotic (0F1 at z >> 0), Kummer
 *     and the DLMF 13.7.2 asymptotic (1F1), Pfaff / 1-z / 1/(1-z) with the
 *     degenerate logarithmic cases, Gauss's sum, the e^(i pi/3) band and the
 *     ODE continuation (2F1, including the cut z > 1 as the limit from below);
 *   - declines: poles, and exact calls staying exact;
 *   - five-surface identity, BIT-FOR-BIT: scalar == unpacked List == packed
 *     List == NDArray == Compile[], real and complex;
 *   - the raw kernels directly (vector driver == scalar kernel).
 *
 * The full ulp sweep (>= 1000 points per regime) is tools/sf_accuracy/hyper.py;
 * this file pins behaviour and catches regressions.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "ndarray.h"
#include "arithmetic.h"
#include "special_functions/sf_machine.h"
#include "special_functions/sf_machine_hyper.h"
#include "test_utils.h"

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers --------------------------------------------------------- */

/* Pull machine values out of a result: a Real, a Complex[Real, Real], or a List
 * / NDArray of those (nested lists flattened in order). Returns the count. */
static int collect(const Expr* e, double* re, double* im, int cap) {
    if (!e || cap <= 0) return 0;
    if (e->type == EXPR_REAL) { re[0] = e->data.real; im[0] = 0.0; return 1; }
    Expr *a, *b;
    if (is_complex((Expr*)e, &a, &b) && a->type == EXPR_REAL && b->type == EXPR_REAL) {
        re[0] = a->data.real; im[0] = b->data.real; return 1;
    }
    if (is_ndarray(e)) {
        Expr* l = ndarray_to_nested_list(e);
        int n = collect(l, re, im, cap);
        expr_free(l);
        return n;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "List") == 0) {
        int n = 0;
        for (size_t i = 0; i < e->data.function.arg_count && n < cap; i++)
            n += collect(e->data.function.args[i], re + n, im + n, cap - n);
        return n;
    }
    return 0;
}

static int eval_vals(const char* input, double* re, double* im, int cap) {
    Expr* p = parse_expression(input);
    ASSERT_MSG(p != NULL, "parse failed: %s", input);
    Expr* r = evaluate(p);
    expr_free(p);
    int n = collect(r, re, im, cap);
    expr_free(r);
    return n;
}

/* A machine literal for x: 17 significant digits plus the backtick, with *^
 * for an exponent (a bare 17-digit literal would parse as arbitrary precision). */
static void mlit(char* buf, size_t sz, double x) {
    char t[64];
    snprintf(t, sizeof t, "%.17g", x);
    char* e = strchr(t, 'e');
    if (e) { *e = '\0'; snprintf(buf, sz, "%s`*^%d", t, atoi(e + 1)); }
    else   snprintf(buf, sz, "%s`", t);
}

/* ---- spot values vs mpmath ------------------------------------------ */

typedef struct { const char* expr; double re, im; int cplx; } Ref;
static const Ref REFS[] = {
    { "Hypergeometric0F1[2, -4.2`]", -0.050175796896052394, 0.0, 0 },
    { "Hypergeometric0F1[2, -0.5`]", 0.7699866217445036, 0.0, 0 },
    { "Hypergeometric0F1[2, 1.5`]", 1.9627864279361782, 0.0, 0 },
    { "Hypergeometric0F1[2, 3.03`]", 3.5057046836944665, 0.0, 0 },
    { "Hypergeometric0F1[2, -300.0`]", 0.004989342033853385, 0.0, 0 },
    { "Hypergeometric0F1[2, 1000.0`]", 4.623747766424143e+24, 0.0, 0 },
    { "Hypergeometric0F1[3/10, -25.1`]", -1.4098059894471846, 0.0, 0 },
    { "Hypergeometric1F1[1, 2, -4.5`]", 0.21975355632483504, 0.0, 0 },
    { "Hypergeometric1F1[1, 2, 0.001`]", 1.0005001667083417, 0.0, 0 },
    { "Hypergeometric1F1[1, 2, 2.96`]", 6.181747214696878, 0.0, 0 },
    { "Hypergeometric1F1[0.3`, 1.7`, -31.0`]", 0.364102303386958, 0.0, 0 },
    { "Hypergeometric1F1[0.3`, 1.7`, -10.0`]", 0.5068221814769915, 0.0, 0 },
    { "Hypergeometric1F1[0.3`, 1.7`, 15.0`]", 24128.873694952043, 0.0, 0 },
    { "Hypergeometric1F1[0.3`, 1.7`, 45.0`]", 5.2608229323052184e+16, 0.0, 0 },
    { "Hypergeometric1F1[-3, 2, 24.08`]", -326.9760213333332, 0.0, 0 },
    { "Hypergeometric1F1[2.5`, 0.7`, -52.9`]", 2.3037691409941514e-05, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, -0.9`]", 0.6373978119200129, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, -0.0145`]", 0.9904372534375129, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, 0.5`]", 1.5451774444795625, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, 0.89`]", 3.326031847468049, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, 0.977`]", 5.856838426182938, 0.0, 0 },
    { "Hypergeometric2F1[1, 2, 3, -12.0`]", 0.13104237003525643, 0.0, 0 },
    { "Hypergeometric2F1[1.5`, 2.25`, 3.1`, 0.95`]", 16.1487147925464, 0.0, 0 },
    { "Hypergeometric2F1[1.5`, 2.25`, 3.1`, -22.5`]", 0.019152367458216545, 0.0, 0 },
    { "Hypergeometric2F1[0.5`, 0.5`, 1, -15.35`]", 0.44265897120283754, 0.0, 0 },
    { "Hypergeometric2F1[1, 0.5`, 3, 1.`]", 1.3333333333333333, 0.0, 0 },
    { "HypergeometricPFQ[{1, 2}, {3, 4}, -0.96`]", 0.8603962347900486, 0.0, 0 },
    { "HypergeometricPFQ[{1, 2}, {3, 4}, 2.99`]", 1.856801395200482, 0.0, 0 },
    { "HypergeometricPFQ[{1, 2}, {3, 4}, 41.9`]", 6442915131440.583, 0.0, 0 },
    { "HypergeometricPFQ[{0.5`, 1, 1.5`}, {2, 2.5`}, 0.8356`]", 1.1950551569664463, 0.0, 0 },
    { "HypergeometricU[1, 2, 0.412`]", 2.4271844660194177, 0.0, 0 },
    { "HypergeometricU[0.7`, 1.3`, 36.5`]", 0.08001180119651184, 0.0, 0 },
    { "HypergeometricU[1.5`, 3, 5.68`]", 0.08278164594343934, 0.0, 0 },
    { "WhittakerM[1, 1/2, 6.87`]", 0.22138856829326842, 0.0, 0 },
    { "WhittakerW[1, 1/2, 6.87`]", 0.22138856829326842, 0.0, 0 },
    { "WhittakerM[0.3`, 0.7`, 4.25`]", 6.385747176724958, 0.0, 0 },
    { "Hypergeometric0F1[2, Complex[-0.85`, 3.0`]]", 0.05176266859664211, 0.9623708444099933, 1 },
    { "Hypergeometric1F1[1, 2, Complex[-2.93`, 1.61`]]", 0.2703798111745964, 0.1303602342293562, 1 },
    { "Hypergeometric1F1[0.3`, 1.7`, Complex[2.36`, -4.25`]]", 0.22639357431652352, -0.5680206613592619, 1 },
    { "Hypergeometric2F1[1, 2, 3, Complex[0.5`, 0.8`]]", 0.863611955538838, 0.683321901354877, 1 },
    { "Hypergeometric2F1[1, 2, 3, Complex[0.8`, 0.8`]]", 0.8215900994813007, 0.948701186865637, 1 },
    { "Hypergeometric2F1[1, 2, 3, Complex[0.0187`, -0.297`]]", 0.9691291242316813, -0.19320242306480656, 1 },
    { "Hypergeometric2F1[1, 2, 3, Complex[3.74`, -2.15`]]", -0.23944746590368174, -0.4808631080717877, 1 },
    { "Hypergeometric2F1[1, 2, 3, Complex[0.5`, 0.866`]]", 0.8138177743517941, 0.6848532561641409, 1 },
    { "Hypergeometric2F1[1, 2, 3, 2.5`]", -0.9297488345946125, -1.0053096491487339, 1 },
    { "HypergeometricPFQ[{1, 2}, {3, 4}, Complex[1.97`, 0.587`]]", 1.438384415080832, 0.18616580841185426, 1 },
    { "HypergeometricU[0.7`, 1.3`, Complex[3.5`, 6.7`]]", 0.176857997591449, -0.15948789097946955, 1 },
    { "WhittakerM[0.3`, 0.7`, Complex[2.27`, -3.82`]]", 0.31452911261205085, -3.012736014166196, 1 },
    { "WhittakerW[0.3`, 0.7`, Complex[1.82`, 7.04`]]", -0.7439794188829774, 0.02944006870254798, 1 },
};

void test_hyper_spot_values() {
    for (size_t i = 0; i < sizeof(REFS) / sizeof(REFS[0]); i++) {
        const Ref* r = &REFS[i];
        double re[2], im[2];
        int n = eval_vals(r->expr, re, im, 2);
        ASSERT_MSG(n == 1, "%s: not a machine number", r->expr);
        /* Relative to the magnitude, floored near a zero of an oscillatory
         * function (0F1 at negative z) by the local envelope (~0.2 there). */
        double mag = hypot(r->re, r->im);
        if (mag < 0.2) mag = 0.2;
        double err = hypot(re[0] - r->re, im[0] - r->im);
        ASSERT_MSG(err <= 1.0e-14 * mag, "%s: got %.17g%+.17gi, want %.17g%+.17gi",
                   r->expr, re[0], im[0], r->re, r->im);
        if (!r->cplx) ASSERT_MSG(im[0] == 0.0, "%s: real input gave a complex value", r->expr);
    }
}

/* ---- declines and exactness ----------------------------------------- */

void test_hyper_declines() {
    /* Exact calls stay exact / symbolic. */
    assert_eval_eq("Hypergeometric1F1[1, 2, 1]", "-1 + E", 0);
    assert_eval_eq("Hypergeometric1F1[-3, 2, x]", "1 - 3/2 x + 1/2 x^2 - 1/24 x^3", 0);
    /* Gauss's sum diverges at c - a - b <= 0: positive terms, so Infinity, as
     * the exact call (audit Hypergeometric2F1#3). */
    assert_eval_eq("Hypergeometric2F1[1, 2, 3, 1.]", "Infinity", 0);
    /* Raw kernel declines: lower-parameter poles, divergent 2F1 at z = 1. */
    double v;
    ASSERT(!sf_machine_hyper1f1(0.5, -2.0, 1.0, &v));
    ASSERT(!sf_machine_hyper2f1(0.5, 0.7, 0.0, 0.3, &v));
    ASSERT(!sf_machine_hyper0f1(-3.0, 1.0, &v));
    ASSERT(!sf_machine_hyper2f1(1.0, 2.0, 3.0, 2.5, &v));    /* real kernel: the cut */
    /* ...but a terminating upper parameter ends the series before the pole. */
    ASSERT(sf_machine_hyper1f1(-2.0, -3.0, 1.5, &v));
    ASSERT_MSG(fabs(v - (1.0 + (-2.0) / (-3.0) * 1.5 + (-2.0) * (-1.0) / ((-3.0) * (-2.0)) * 1.5 * 1.5 / 2.0)) < 1e-15,
               "1F1(-2;-3;1.5) = %.17g", v);
    double rr, ri;
    /* Complex parameters outside |z| <= 0.9 are not transformed: decline. */
    ASSERT(!sf_machine_hyper2f1_c(1.0, 0.5, 2.0, 0.0, 3.0, 0.0, 2.0, 1.0, &rr, &ri));
}

/* ---- five surfaces, bit for bit --------------------------------------- */

/* For `head` (a format with one %s for the argument), compare the scalar
 * evaluation at each point with the unpacked List, the packed List, the
 * NDArray and the Compile[] (argument spec `spec`) results -- exactly. */
static void check_surfaces(const char* head, const char* spec,
                           const double* xr, const double* xi, int n, int cplx) {
    char lits[16][96];
    for (int i = 0; i < n; i++) {
        if (cplx) {
            char a[48], b[48];
            mlit(a, sizeof a, xr[i]); mlit(b, sizeof b, xi[i]);
            snprintf(lits[i], sizeof lits[i], "Complex[%s, %s]", a, b);
        } else {
            mlit(lits[i], sizeof lits[i], xr[i]);
        }
    }
    char list[2048] = "{";
    for (int i = 0; i < n; i++) {
        strcat(list, lits[i]);
        strcat(list, i + 1 < n ? ", " : "}");
    }
    double sr[16], si[16];
    for (int i = 0; i < n; i++) {
        char e[512], c[1024];
        snprintf(e, sizeof e, head, lits[i]);
        ASSERT_MSG(eval_vals(e, &sr[i], &si[i], 1) == 1, "scalar %s", e);
        /* Compile[] */
        char body[256];
        snprintf(body, sizeof body, head, "x");
        snprintf(c, sizeof c, "Compile[{{x, %s}}, %s][%s]", spec, body, lits[i]);
        double cr, ci;
        ASSERT_MSG(eval_vals(c, &cr, &ci, 1) == 1, "compile %s", c);
        ASSERT_MSG(cr == sr[i] && ci == si[i], "%s: Compile %.17g%+.17gi != scalar %.17g%+.17gi",
                   c, cr, ci, sr[i], si[i]);
    }
    const char* wraps[3] = { "FromPackedArray[%s]", "ToPackedArray[%s]", "NDArray[%s]" };
    for (int w = 0; w < 3; w++) {
        char arg[2100], e[2600];
        snprintf(arg, sizeof arg, wraps[w], list);
        snprintf(e, sizeof e, head, arg);
        double vr[16], vi[16];
        int m = eval_vals(e, vr, vi, 16);
        ASSERT_MSG(m == n, "%s: %d values", e, m);
        for (int i = 0; i < n; i++)
            ASSERT_MSG(vr[i] == sr[i] && vi[i] == si[i],
                       "%s [%d]: %.17g%+.17gi != scalar %.17g%+.17gi", e, i, vr[i], vi[i],
                       sr[i], si[i]);
    }
}

void test_hyper_surfaces_real() {
    double x0[] = { -4.2, -1.1, 0.3, 2.7, 4.9 };
    check_surfaces("Hypergeometric0F1[2, %s]", "_Real", x0, NULL, 5, 0);
    double x1[] = { -4.5, -0.7, 0.001, 2.96, 4.4 };
    check_surfaces("Hypergeometric1F1[1, 2, %s]", "_Real", x1, NULL, 5, 0);
    check_surfaces("Hypergeometric1F1[3/10, 17/10, %s]", "_Real", x1, NULL, 5, 0);
    double x2[] = { -0.9, -0.31, 0.12, 0.6, 0.89 };
    check_surfaces("Hypergeometric2F1[1, 2, 3, %s]", "_Real", x2, NULL, 5, 0);
    check_surfaces("Hypergeometric2F1[3/2, 9/4, 31/10, %s]", "_Real", x2, NULL, 5, 0);
    double xp[] = { -2.9, -0.96, 0.4, 2.99 };
    check_surfaces("HypergeometricPFQ[{1, 2}, {3, 4}, %s]", "_Real", xp, NULL, 4, 0);
    double xu[] = { 0.25, 1.5, 4.0, 7.9 };
    check_surfaces("HypergeometricU[1, 2, %s]", "_Real", xu, NULL, 4, 0);
    check_surfaces("WhittakerM[1, 1/2, %s]", "_Real", xu, NULL, 4, 0);
    check_surfaces("WhittakerW[1, 1/2, %s]", "_Real", xu, NULL, 4, 0);
    check_surfaces("WhittakerM[3/10, 7/10, %s]", "_Real", xu, NULL, 4, 0);
}

void test_hyper_surfaces_complex() {
    double re[] = { -2.9, -0.4, 0.5, 1.8 }, im[] = { 1.6, -2.2, 0.8, -0.3 };
    check_surfaces("Hypergeometric0F1[2, %s]", "_Complex", re, im, 4, 1);
    check_surfaces("Hypergeometric1F1[1, 2, %s]", "_Complex", re, im, 4, 1);
    check_surfaces("Hypergeometric1F1[3/10, 17/10, %s]", "_Complex", re, im, 4, 1);
    double r2[] = { -0.7, 0.5, 0.8, 0.02 }, i2[] = { 0.6, 0.8, 0.8, -0.3 };
    check_surfaces("Hypergeometric2F1[1, 2, 3, %s]", "_Complex", r2, i2, 4, 1);
    check_surfaces("HypergeometricPFQ[{1, 2}, {3, 4}, %s]", "_Complex", re, im, 4, 1);
    double r3[] = { 0.6, 1.5, 3.9 }, i3[] = { -1.9, 0.4, 1.2 };
    check_surfaces("HypergeometricU[1, 2, %s]", "_Complex", r3, i3, 3, 1);
    check_surfaces("WhittakerM[1, 1/2, %s]", "_Complex", r3, i3, 3, 1);
    check_surfaces("WhittakerW[1, 1/2, %s]", "_Complex", r3, i3, 3, 1);
}

/* The vector driver (constant parameters over a buffer, 2F1 coefficient
 * tables) agrees with the scalar kernel exactly, across the regimes. */
void test_hyper_vector_driver() {
    enum { N = 2000 };
    static double x[N], y[N];
    for (int i = 0; i < N; i++) x[i] = -12.0 + 12.98 * (double)i / (double)(N - 1);
    const double sets[][3] = { { 1, 2, 3 }, { 1.5, 2.25, 3.1 }, { 0.5, 0.5, 1 },
                               { -3, 0.7, 1.9 }, { 20, 20, 1.5 } };
    for (size_t s = 0; s < sizeof sets / sizeof sets[0]; s++) {
        const double* p = sets[s];
        double ab[2] = { p[0], p[1] };
        size_t done = sf_machine_pfq_vec(ab, 2, &p[2], 1, x, y, N);
        for (size_t i = 0; i < done; i++) {
            double v;
            ASSERT(sf_machine_hyper2f1(p[0], p[1], p[2], x[i], &v));
            ASSERT_MSG(v == y[i], "2F1(%g,%g;%g;%.17g): vec %.17g != scalar %.17g",
                       p[0], p[1], p[2], x[i], y[i], v);
        }
        if (done < N) {   /* the first decline must also decline in scalar */
            double v;
            ASSERT(!sf_machine_hyper2f1(p[0], p[1], p[2], x[done], &v));
        }
    }
}

/* Raw complex kernel: real arguments on the real axis below 1 agree with the
 * real kernel to a few ulp, and the cut z > 1 is the limit from below. */
void test_hyper_branch() {
    double rr, ri, v;
    ASSERT(sf_machine_hyper2f1(1, 1, 2, 0.5, &v));
    ASSERT(sf_machine_hyper2f1_c(1, 0, 1, 0, 2, 0, 0.5, 0.0, &rr, &ri));
    ASSERT_MSG(fabs(rr - v) <= 4e-16 * fabs(v) && ri == 0.0, "2F1 axis: %.17g vs %.17g", rr, v);
    /* 2F1(1,1;2;2) = -log(1-2)/2 = -i pi/2 from below. */
    ASSERT(sf_machine_hyper2f1_c(1, 0, 1, 0, 2, 0, 2.0, -0.0, &rr, &ri));
    ASSERT_MSG(fabs(rr) < 1e-15 && fabs(ri + M_PI / 2) < 1e-15,
               "2F1(1,1;2;2-i0) = %.17g%+.17gi", rr, ri);
}

void test_hyper_attributes() {
    const char* heads[] = { "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1",
                            "HypergeometricPFQ", "HypergeometricU", "WhittakerM",
                            "WhittakerW" };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        SymbolDef* d = symtab_get_def(heads[i]);
        ASSERT_MSG(d && (d->attributes & ATTR_NUMERICFUNCTION), "%s attributes", heads[i]);
        ASSERT_MSG(d->docstring != NULL, "%s has no docstring", heads[i]);
    }
}

int main() {
    symtab_init();
    core_init();

    TEST(test_hyper_spot_values);
    TEST(test_hyper_declines);
    TEST(test_hyper_surfaces_real);
    TEST(test_hyper_surfaces_complex);
    TEST(test_hyper_vector_driver);
    TEST(test_hyper_branch);
    TEST(test_hyper_attributes);

    printf("All hypergeometric machine-kernel tests passed.\n");
    return 0;
}
