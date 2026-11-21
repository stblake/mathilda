/* Tests for the zeta-family machine kernels (src/special_functions/
 * sf_machine_zeta.c) behind Zeta, HurwitzZeta / two-argument Zeta, LerchPhi
 * and PolyLog at machine precision.
 *
 * Covers:
 *   - spot values against mpmath (50 digits), real and complex, each regime of
 *     each kernel: tabulated integers, Euler-Maclaurin, the direct Dirichlet
 *     sum (s >= 20), the functional equation (s < 1/2, real and complex);
 *     the dilogarithm's Bernoulli series and its 1-z / 1/z maps, the cut
 *     convention (real x > 1 from below), the Jonquiere / inversion branches
 *     of PolyLog for other orders; real and complex LerchPhi series;
 *   - regime boundaries (s = 1/2, s = 20, x = 1/2, the EM shift W);
 *   - declines that must fall through to the arbitrary-precision path with
 *     the old answer (poles, Re a <= 0, complex order, |z| >= 1);
 *   - the five-surface bit-match: scalar == unpacked List == packed List ==
 *     NDArray == Compile[] (scalar and rank-1), real and complex.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "sym_names.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate `input` to a machine number; returns false if the result is not a
 * machine Real or Complex[Real, Real]. */
static bool eval_num(const char* input, double* re, double* im) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    bool ok = false;
    if (r->type == EXPR_REAL) {
        *re = r->data.real; *im = 0.0; ok = true;
    } else if (r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
               r->data.function.head->type == EXPR_SYMBOL &&
               r->data.function.head->data.symbol.name == SYM_Complex &&
               r->data.function.args[0]->type == EXPR_REAL &&
               r->data.function.args[1]->type == EXPR_REAL) {
        *re = r->data.function.args[0]->data.real;
        *im = r->data.function.args[1]->data.real;
        ok = true;
    }
    expr_free(r);
    return ok;
}

typedef struct { const char* in; double re, im; } Spot;

/* Relative error (modulus) <= 8 ulp of the reference. */
static void check_spots(const Spot* t, size_t n) {
    for (size_t i = 0; i < n; i++) {
        double re, im;
        ASSERT_MSG(eval_num(t[i].in, &re, &im), "%s: not a machine number", t[i].in);
        double err = hypot(re - t[i].re, im - t[i].im);
        double mag = hypot(t[i].re, t[i].im);
        ASSERT_MSG(err <= 8.0 * 2.220446049250313e-16 * mag,
                   "%s: expected %.17g%+.17gI, got %.17g%+.17gI (rel %.3g)",
                   t[i].in, t[i].re, t[i].im, re, im, err / mag);
    }
}

void test_zeta_machine_spots() {
    static const Spot t[] = {
        { "Zeta[2.5]", 1.3414872572509171798, 0.0 },
        { "Zeta[2.]", 1.6449340668482264365, 0.0 },
        { "Zeta[-3.5]", 0.0044410113354794319585, 0.0 },
        { "Zeta[0.25]", -0.81327840526189165652, 0.0 },
        { "Zeta[19.99]", 1.0000009605985443501, 0.0 },
        { "Zeta[20.01]", 1.0000009473713779418, 0.0 },
        { "Zeta[0.4999]", -1.4599623242214536584, 0.0 },
        { "Zeta[0.5001]", -1.4607468534812965067, 0.0 },
        { "Zeta[-31.02148751594208]", 488760602.97341230981, 0.0 },
        { "Zeta[Complex[2.,1.]]", 1.1503557032549026717, -0.43753086591960788112 },
        { "Zeta[Complex[0.3,5.]]", 0.67564899811602329843, 0.25414478655467744161 },
        { "Zeta[Complex[-2.5,3.]]", 0.068763679033646481628, 0.13398028393783442697 },
    };
    check_spots(t, sizeof t / sizeof t[0]);
    /* Trivial zero and the exact table entries. */
    assert_eval_eq("Zeta[-4.]", "0.0", 0);
    assert_eval_eq("Zeta[0.]", "-0.5", 0);
}

void test_hurwitz_machine_spots() {
    static const Spot t[] = {
        { "HurwitzZeta[2.5,0.7]", 2.9028675777573466296, 0.0 },
        { "Zeta[2.5,0.7]", 2.9028675777573466296, 0.0 },
        { "HurwitzZeta[3.5,9.6]", 0.001594251119821890448, 0.0 },     /* EM shift N = 1 */
        { "HurwitzZeta[3.5,10.4]", 0.0012923057248894297732, 0.0 },   /* no shift */
        { "HurwitzZeta[Complex[2.5,1.],Complex[0.7,-0.3]]",
          0.63507590433651115197, 1.1053756047127642544 },
    };
    check_spots(t, sizeof t / sizeof t[0]);
    /* a = 1 is the Riemann kernel, bit for bit. */
    assert_eval_eq("HurwitzZeta[3.7, 1.] === Zeta[3.7]", "True", 0);
}

void test_lerchphi_machine_spots() {
    static const Spot t[] = {
        { "LerchPhi[0.5,2.,0.7]", 2.2616091544648781387, 0.0 },
        { "LerchPhi[Complex[0.3,0.4],2.5,1.2]", 0.66791526031994653202, 0.069120226233962508687 },
        { "LerchPhi[Complex[-0.5,0.5],Complex[2.,0.5],Complex[1.5,-0.5]]",
          0.24513406428052630723, 0.18290800662703176113 },
    };
    check_spots(t, sizeof t / sizeof t[0]);
}

void test_polylog_machine_spots() {
    static const Spot t[] = {
        { "PolyLog[2,0.5]", 0.5822405264650125059, 0.0 },
        { "PolyLog[2,0.5000001]", 0.58224066509445468198, 0.0 },
        { "PolyLog[2,-1.]", -0.82246703342411321824, 0.0 },
        { "PolyLog[2,-3.]", -1.9393754207667089531, 0.0 },
        { "PolyLog[2,0.9]", 1.299714723004958782, 0.0 },
        { "PolyLog[2,1.*^-8]", 1.000000002500000032e-8, 0.0 },
        { "PolyLog[2,2.]", 2.4674011002723396547, -2.1775860903036021305 },
        { "PolyLog[2,Complex[0.3,0.6]]", 0.20173298164594456678, 0.6702970657586272647 },
        { "PolyLog[2,Complex[3.,-4.]]", -0.60480701206119998349, -3.7336195322943858856 },
        { "PolyLog[2,Complex[0.9,0.3]]", 1.1049863515242157246, 0.61705302808486198385 },
        { "PolyLog[3,0.7]", 0.78006393425766150446, 0.0 },
        { "PolyLog[3,Complex[5.,2.]]", 3.028439989431801062, 4.8747802827637167147 },
        { "PolyLog[2.5,2.]", 2.789660332382777144, -1.3638037005393527944 },
        { "PolyLog[2.5,Complex[0.6,0.5]]", 0.58863942850249106103, 0.63142973911595969942 },
        { "PolyLog[2.5,-3.]", -2.1627007120020566623, 0.0 },
    };
    check_spots(t, sizeof t / sizeof t[0]);
    /* PolyLog[1, x > 1] = -Log[1 - x], from below the cut: -I Pi. */
    double re, im;
    ASSERT(eval_num("PolyLog[1, 2.]", &re, &im));
    ASSERT(fabs(re) < 1e-15 && fabs(im + 3.141592653589793) < 1e-15);
    /* x = 1 is pi^2/6 exactly as the double rounds it. */
    assert_eval_eq("PolyLog[2, 1.] === Zeta[2.]", "True", 0);
    /* The exact zero argument stays exact. */
    assert_eval_eq("PolyLog[2., 0]", "0", 0);
}

/* A decline must fall through to the arbitrary-precision path. */
void test_zeta_family_declines() {
    assert_eval_eq("Zeta[1.]", "ComplexInfinity", 0);
    assert_eval_eq("HurwitzZeta[1., 0.5]", "ComplexInfinity", 0);
    double re, im;
    /* s < 1/2: the kernel declines, the MPFR Euler-Maclaurin answers. */
    ASSERT(eval_num("HurwitzZeta[0.3, 0.5]", &re, &im));
    ASSERT(fabs(re - (-0.20908381885368533)) < 1e-14);
    /* Re a <= 0: declined, still numeric. */
    ASSERT(eval_num("HurwitzZeta[2., -0.5]", &re, &im));
    ASSERT(fabs(re - 8.934802200544679) < 1e-13);
    /* Complex order: declined by the kernel, MPFR path. */
    ASSERT(eval_num("PolyLog[Complex[2.,1.], 0.5]", &re, &im));
    ASSERT(fabs(re - 0.55485354535332032) < 1e-13 && fabs(im - (-0.058049561791697155)) < 1e-13);
    /* LerchPhi at z = 1 is the Hurwitz kernel. */
    ASSERT(eval_num("LerchPhi[1., 2., 0.7]", &re, &im));
    ASSERT(fabs(re - 2.8340491566946109) < 1e-14);
}

/* Re s >= 20, complex: the Dirichlet series on the machine (v0.6662).  Every
 * |s| >= 10^4 used to decline to the arbitrary-precision path -- 1.7 s for
 * Zeta[1.*^6 + 5. I], and FindRoot's iterate near Re s = 10^189 never
 * returned.  References at 40 digits. */
void test_zeta_large_real_part_complex() {
    double re, im;
    ASSERT(eval_num("Zeta[20. + 3. I]", &re, &im));
    ASSERT(fabs(re - 0.999999535282062328) < 4e-16 && fabs(im - (-8.32899148367607721e-7)) < 4e-16);
    ASSERT(eval_num("Zeta[25.5 - 7. I]", &re, &im));
    ASSERT(fabs(re - 1.00000000293327626) < 4e-16 && fabs(im - (-2.08676231340347526e-8)) < 4e-16);
    ASSERT(eval_num("Zeta[1000. + 5. I]", &re, &im));
    ASSERT(re == 1.0 && fabs(im - 2.97241446425992725e-302) < 1e-315);
    /* far out: n^-Re s underflows, so the answer is exactly 1, at once */
    assert_eval_eq("Zeta[1.06*^189 + 1.1*^188 I]", "1.0", 0);
    assert_eval_eq("First[AbsoluteTiming[Do[Zeta[1.*^6 + 5. I], {200}]]] < 0.5", "True", 0);
}

/* scalar == unpacked List == packed List == NDArray == Compile[], bit for bit. */
void test_zeta_family_surfaces() {
    static const char* real_cases[] = {
        /* name, body over x (scalar) */
        "Zeta[#]&", "HurwitzZeta[2.5, #]&", "HurwitzZeta[#, 1.7]&",
        "PolyLog[2, #]&", "PolyLog[3, #]&", "PolyLog[2.5, #]&",
        "LerchPhi[#, 2.5, 1.2]&",
    };
    static const char* real_data[] = {
        "{1.5, 2.25, 3.7, 6.1, 8.9}",                 /* Zeta, Hurwitz s */
        "{0.55, 1.3, 2.75, 3.9, 9.3}",                /* Hurwitz a */
        "{1.6, 2.4, 4.5, 7.5, 21.5}",
        "{-0.85, -0.3, 0.2, 0.6, 0.88}",              /* PolyLog / LerchPhi z */
        "{-0.95, -0.5, 0.1, 0.55, 0.97}",
        "{-0.9, -0.2, 0.3, 0.45, 0.8}",
        "{-0.8, -0.1, 0.2, 0.5, 0.85}",
    };
    char buf[2048];
    for (size_t i = 0; i < sizeof real_cases / sizeof real_cases[0]; i++) {
        snprintf(buf, sizeof buf,
            "Module[{f = %s, xs = %s, u, p, n, cs, cv}, u = FromPackedArray[xs];"
            " p = ToPackedArray[xs]; n = ToNDArray[xs];"
            " cs = Compile[{{x, _Real}}, Evaluate[f[x]]];"
            " cv = Compile[{{v, _Real, 1}}, Evaluate[f[v]]];"
            " Map[f, u] === f[u] === f[p] === Normal[f[n]] === Map[cs, u] === cv[p]]",
            real_cases[i], real_data[i]);
        assert_eval_eq(buf, "True", 0);
    }
    static const char* cplx_cases[] = {
        "Zeta[#]&", "HurwitzZeta[Complex[2.5, 0.5], #]&", "PolyLog[2, #]&",
        "PolyLog[3, #]&", "LerchPhi[#, 2.5, 1.2]&",
    };
    static const char* cplx_data[] = {
        "{Complex[1.5, -2.9], Complex[2.2, 0.4], Complex[4.9, 2.7], Complex[0.2, 7.5]}",
        "{Complex[0.6, -0.9], Complex[1.4, 0.2], Complex[2.9, 0.95], Complex[8.5, -3.]}",
        "{Complex[-0.7, 0.7], Complex[0.6, -0.2], Complex[0.79, 0.5], Complex[3., -9.]}",
        "{Complex[-0.7, 0.7], Complex[0.6, -0.2], Complex[1.5, 0.5], Complex[12., -9.]}",
        "{Complex[-0.6, 0.55], Complex[0.5, -0.2], Complex[0.1, 0.6], Complex[0.59, 0.59]}",
    };
    for (size_t i = 0; i < sizeof cplx_cases / sizeof cplx_cases[0]; i++) {
        snprintf(buf, sizeof buf,
            "Module[{f = %s, xs = %s, u, p, n, cs, cv}, u = FromPackedArray[xs];"
            " p = ToPackedArray[xs]; n = ToNDArray[xs];"
            " cs = Compile[{{x, _Complex}}, Evaluate[f[x]]];"
            " cv = Compile[{{v, _Complex, 1}}, Evaluate[f[v]]];"
            " Map[f, u] === f[u] === f[p] === Normal[f[n]] === Map[cs, u] === cv[p]]",
            cplx_cases[i], cplx_data[i]);
        assert_eval_eq(buf, "True", 0);
    }
    /* All five lower in Compile[] at scalar and rank-1 shape, real and complex. */
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{x, _Real}}, Zeta[x] + HurwitzZeta[x, 2.] + PolyLog[2, x] + LerchPhi[x, 2., 1.5]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Real, 1}}, Zeta[v] + PolyLog[2, v]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Zeta[z] + HurwitzZeta[z, 2.] + PolyLog[2, z] + LerchPhi[z, 2., 1.5]]", "True", 0);
    assert_eval_eq("\"Compiled\" /. CompileDiagnostics[{{v, _Complex, 1}}, Zeta[v] + PolyLog[2, v]]", "True", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_zeta_machine_spots);
    TEST(test_hurwitz_machine_spots);
    TEST(test_lerchphi_machine_spots);
    TEST(test_polylog_machine_spots);
    TEST(test_zeta_family_declines);
    TEST(test_zeta_large_real_part_complex);
    TEST(test_zeta_family_surfaces);

    printf("All zeta-family machine-kernel tests passed.\n");
    return 0;
}
