/* Tests for the Legendre elliptic integrals: EllipticK, EllipticF, EllipticE,
 * EllipticPi.
 *
 * The convention is the trap here, so it is asserted directly: the parameter
 * argument is m = k^2, NOT the modulus k. EllipticK[1/2] is 1.8540746773...,
 * which is K at parameter 1/2; read as a modulus it would be 1.6857503548...,
 * the value at m = 1/4. Two of the numeric tests below are exactly that pair,
 * so a convention slip fails rather than merely shifting answers.
 *
 * Reference values are from mpmath at 35 working digits, cross-checked against
 * direct numerical quadrature of the defining integrals (not recalled
 * constants -- the first draft of this file carried two wrong "Mathematica"
 * values from memory, and the quadrature is what caught them).
 *
 * Coverage: exact reductions; exact arguments staying symbolic; machine and
 * 25-digit numerics for all six arities; complex phi (which
 * EllipticF[ArcSin[z], m] with |z| > 1 produces routinely); EllipticPi past the
 * pole at Sin[t]^2 == 1/n, where the value is the Cauchy principal value;
 * derivative rules verified as residuals rather than as printed forms; the
 * packed/NDArray surfaces agreeing with the List path; attributes and arity.
 *
 * The last group is the one that matters most: the seven Mathematica reference
 * antiderivatives of the elliptic block of the ParallelMixedSpecial stress
 * corpus (cases 100-114), differentiated back to their integrands. That single
 * check exercises the numerics, the branch placement and the derivative rules
 * together, against answers produced by a different CAS.
 *
 * Assertions use ASSERT / ASSERT_MSG / ASSERT_STR_EQ (hard exit(1)), NOT
 * assert_eval_eq, whose libc assert() is a no-op under -DNDEBUG.
 */

#include "attr.h"
#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers --------------------------------------------------------- */

/* Evaluate `input` and return its printed form (caller frees). */
static char* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

static void assert_eval(const char* input, const char* expected) {
    char* s = eval_str(input);
    ASSERT_STR_EQ(s, expected);
    free(s);
}

/* Assert that `input` evaluates to the symbol True. Used for every numeric
 * claim, so the comparison happens in the kernel at its own precision rather
 * than being squeezed through a C double. */
static void assert_true(const char* input) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s: expected True, got %s", input, s);
    free(s);
}

/* |expr - ref| < tol, with expr evaluated at `digits` precision. */
static void assert_num(const char* expr, const char* ref, int digits, const char* tol) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Abs[N[%s, %d] - %s] < %s", expr, digits, ref, tol);
    assert_true(buf);
}

/* ---- exact reductions ----------------------------------------------- */

static void test_exact_reductions(void) {
    assert_eval("EllipticK[0]", "1/2 Pi");
    assert_eval("EllipticK[1]", "ComplexInfinity");
    assert_eval("EllipticE[0]", "1/2 Pi");
    assert_eval("EllipticE[1]", "1");

    assert_eval("EllipticF[0, m]", "0");
    assert_eval("EllipticF[phi, 0]", "phi");
    assert_eval("EllipticF[Pi/2, m]", "EllipticK[m]");

    assert_eval("EllipticE[0, m]", "0");
    assert_eval("EllipticE[phi, 0]", "phi");
    assert_eval("EllipticE[Pi/2, m]", "EllipticE[m]");
    assert_eval("EllipticE[phi, 1]", "Sin[phi]");

    assert_eval("EllipticPi[0, m]", "EllipticK[m]");
    assert_eval("EllipticPi[0, phi, m]", "EllipticF[phi, m]");
    assert_eval("EllipticPi[n, 0, m]", "0");
    assert_eval("EllipticPi[n, Pi/2, m]", "EllipticPi[n, m]");
}

/* An exact, non-special argument stays symbolic, as in the Wolfram Language --
 * only inexact input or an explicit N[...] evaluates. */
static void test_exact_stays_symbolic(void) {
    assert_eval("EllipticK[1/2]", "EllipticK[1/2]");
    assert_eval("EllipticF[1/3, 1/2]", "EllipticF[1/3, 1/2]");
    assert_eval("EllipticE[1/3, 1/2]", "EllipticE[1/3, 1/2]");
    assert_eval("EllipticPi[1/3, 1/2]", "EllipticPi[1/3, 1/2]");
    assert_eval("EllipticPi[1/3, 1/4, 1/2]", "EllipticPi[1/3, 1/4, 1/2]");
    assert_eval("EllipticK[m]", "EllipticK[m]");
}

/* ---- numerics -------------------------------------------------------- */

/* PARAMETER, not modulus. If the argument were read as the modulus k,
 * EllipticK[1/2] would come back as 1.68575..., which is the value at m = 1/4.
 * Both are asserted here, so the pair pins the convention rather than merely
 * shifting every answer by a consistent amount. */
static void test_complete(void) {
    assert_num("EllipticK[1/2]",  "1.854074677301371918433850",   25, "10^-22");
    assert_num("EllipticK[1/4]",  "1.685750354812596042871204",   25, "10^-22");
    assert_num("EllipticE[1/2]",  "1.350643881047675502520175",   25, "10^-22");
    assert_num("EllipticE[1/4]",  "1.467462209339427155459795",   25, "10^-22");
    assert_num("EllipticPi[1/3, 1/2]", "2.310794990754281815442315", 25, "10^-22");
}

static void test_incomplete(void) {
    assert_num("EllipticF[1, 1/2]", "1.083216772845168750444132",  25, "10^-22");
    assert_num("EllipticE[1, 1/2]", "0.9273298836244400669659042", 25, "10^-22");
    assert_num("EllipticPi[1/3, 1, 1/2]", "1.206801275668137959983858", 25, "10^-22");
    /* phi past the principal strip: the quasi-period must be applied.
     * F(phi + k Pi | m) = F(phi | m) + 2 k K(m), so F(1 + Pi | 1/2) is
     * F(1 | 1/2) + 2 K(1/2). */
    assert_true("Abs[N[EllipticF[1 + Pi, 1/2], 30] "
                "- N[EllipticF[1, 1/2] + 2 EllipticK[1/2], 30]] < 10^-25");
    /* odd in phi */
    assert_true("Abs[N[EllipticF[-1, 1/2], 30] + N[EllipticF[1, 1/2], 30]] < 10^-25");
}

/* A complex phi, which the corpus produces whenever ArcSin's argument exceeds
 * 1, and a complete K past m = 1 -- both genuinely off the real axis. */
static void test_complex_branch(void) {
    assert_true("Abs[N[EllipticF[ArcSin[Sqrt[2] Sqrt[1/(1 + 1/2)]], 1/2], 30] "
                "- (1.85407467730137191843385 - 0.8260178762492451854562394 I)] < 10^-22");
    assert_true("Abs[N[EllipticK[3/2], 30] "
                "- (1.656638170236594166448468 "
                "- 1.415737208425956198892166 I)] < 10^-22");
}

/* EllipticPi with n > 1: the path crosses the pole at Sin[t]^2 == 1/n, so the
 * value is the Cauchy principal value and genuinely complex. The machine
 * Carlson kernels have no principal-value R_J and must decline here, leaving
 * the FLINT/Arb path to answer -- so this also guards that fall-through. */
static void test_principal_value(void) {
    assert_true("Abs[N[EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[1/(1 + 7/5)]], 1/2], 30] "
                "- (0.9877396997285220802095515 - 2.720699046351326775891117 I)] < 10^-22");
}

/* ---- derivatives ----------------------------------------------------- */

/* The phi-derivatives ARE the integrands, which is what makes a numeric verify
 * of an antiderivative built from these kernels close. Asserted as printed
 * forms (they are canonical and small) plus a residual against a central
 * difference, so a sign slip cannot hide behind a plausible-looking form. */
static void test_derivative_forms(void) {
    assert_eval("D[EllipticF[p, m], p]", "1/Sqrt[1 - m Sin[p]^2]");
    assert_eval("D[EllipticE[p, m], p]", "Sqrt[1 - m Sin[p]^2]");
    assert_eval("D[EllipticPi[n, p, m], p]",
                "1/(Sqrt[1 - m Sin[p]^2] (1 - n Sin[p]^2))");
    /* Deliberately inert rather than guessed -- see the comment in deriv.c. */
    assert_eval("D[EllipticF[p, m], m]", "Derivative[0, 1][EllipticF][p, m]");
    assert_eval("D[EllipticPi[n, m], n]", "Derivative[1, 0][EllipticPi][n, m]");
}

static void test_derivative_residuals(void) {
    /* Central difference at 30 digits, h = 10^-15: agreement to ~10^-14 is the
     * finite-difference floor, not the derivative's accuracy. */
    const char* cases[] = {
        "EllipticF[p, 1/2]",            "p",  "3/5",
        "EllipticE[p, 1/2]",            "p",  "3/5",
        "EllipticPi[1/3, p, 1/2]",      "p",  "3/5",
        "EllipticK[mm]",                "mm", "2/5",
        "EllipticE[mm]",                "mm", "2/5",
        "EllipticE[3/5, mm]",           "mm", "2/5",
        /* the chain rule through the composite arguments the corpus builds */
        "EllipticF[ArcSin[Sqrt[2] Sqrt[1/(1 + xx)]], 1/2]",       "xx", "7/5",
        "EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[1/(1 + xx)]], 1/2]", "xx", "7/5",
        NULL
    };
    char buf[2048];
    for (int i = 0; cases[i]; i += 3) {
        snprintf(buf, sizeof(buf),
            "Module[{h = 10^-15, f = %s, v = %s, v0 = %s},"
            " Abs[N[D[f, v] /. v -> v0, 30]"
            " - N[((f /. v -> v0 + h) - (f /. v -> v0 - h))/(2 h), 30]] < 10^-12]",
            cases[i], cases[i + 1], cases[i + 2]);
        assert_true(buf);
    }
}

/* ---- the corpus elliptic block, end to end -------------------------- */

/* The seven Mathematica reference antiderivatives of cases 100-114 of the
 * ParallelMixedSpecial stress corpus, differentiated back to their integrands
 * at the corpus sample points. This is the check that ties the numerics, the
 * branch placement and the derivative rules together against another CAS's
 * answers; each residual is 0 or ~10^-29 in practice. */
static void test_corpus_residuals(void) {
    struct { const char* label; const char* ans; const char* f; } c[] = {
        { "#100", "-(Sqrt[2] EllipticF[ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2])",
                  "1/Sqrt[x^3 - x]" },
        { "#101", "(2 Sqrt[-x + x^3])/(1 + x) "
                  "+ 2 Sqrt[2] EllipticE[ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2] "
                  "- Sqrt[2] EllipticF[ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2]",
                  "x/Sqrt[x^3 - x]" },
        { "#102", "(Sqrt[2] EllipticF[ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2])/3 "
                  "- (Sqrt[2] EllipticPi[3/2, ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2])/3",
                  "1/((x - 2) Sqrt[x^3 - x])" },
        { "#103", "EllipticF[ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2]/(2 Sqrt[2]) "
                  "- EllipticPi[2, ArcSin[Sqrt[2] Sqrt[(1 + x)^(-1)]], 1/2]/(2 Sqrt[2])",
                  "1/((x - 3) Sqrt[x^3 - x])" },
        { "#105", "(2 Sqrt[-4 x + x^3])/(2 + x) "
                  "+ 4 EllipticE[ArcSin[2 Sqrt[(2 + x)^(-1)]], 1/2] "
                  "- 2 EllipticF[ArcSin[2 Sqrt[(2 + x)^(-1)]], 1/2]",
                  "x/Sqrt[x^3 - 4 x]" },
        { "#113", "(-2 Sqrt[-E^x + E^(3 x)])/E^x + (4 Sqrt[-E^x + E^(3 x)])/(1 + E^x) "
                  "+ 4 Sqrt[2] EllipticE[ArcSin[Sqrt[2] Sqrt[(1 + E^x)^(-1)]], 1/2] "
                  "- 2 Sqrt[2] EllipticF[ArcSin[Sqrt[2] Sqrt[(1 + E^x)^(-1)]], 1/2]",
                  "Sqrt[Exp[3 x] - Exp[x]] Exp[-x]" },
        { "#114", "-EllipticF[ArcSin[2 Sqrt[(2 + E^x)^(-1)]], 1/2]",
                  "Exp[x]/Sqrt[Exp[3 x] - 4 Exp[x]]" },
    };
    char buf[4096];
    for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++) {
        snprintf(buf, sizeof(buf),
            "Max[Table[Abs[N[(D[%s, x] - (%s)) /. x -> p, 30]], "
            "{p, {7/5, 9/4, 13/5}}]] < 10^-20", c[i].ans, c[i].f);
        char* s = eval_str(buf);
        ASSERT_MSG(strcmp(s, "True") == 0,
                   "corpus case %s: D[answer] - integrand did not vanish (%s)",
                   c[i].label, s);
        free(s);
    }
}

/* ---- packed / NDArray surfaces -------------------------------------- */

/* The plain List, the packed List and the visible NDArray must agree. The
 * machine Carlson kernels serve the buffer; where they decline (m > 1, or
 * 1 - m Sin[phi]^2 < 0) the array is abandoned and the List path answers
 * through FLINT, which is slower but never wrong -- the last pair checks that
 * fall-through rather than only the fast case. */
static void test_nd_surfaces(void) {
    /* Compare through Normal on the array side: a visible NDArray and a plain
     * List do not combine element-wise under Listable threading (the List
     * threads and each element meets the whole array), so `list == ndarray`
     * stays unevaluated and a subtraction builds an outer product. Normal
     * materialises the buffer to the List it presents as, which is the
     * comparison this test actually wants. */
    assert_true("Normal[EllipticK[NDArray[{0.25, 0.5, 0.75}]]] "
                "== EllipticK[{0.25, 0.5, 0.75}]");
    assert_true("Normal[EllipticE[NDArray[{0.25, 0.5, 0.75}]]] "
                "== EllipticE[{0.25, 0.5, 0.75}]");
    assert_true("Normal[EllipticF[NDArray[{0.3, 0.6, 0.9}], 0.5]] "
                "== EllipticF[{0.3, 0.6, 0.9}, 0.5]");
    assert_true("Normal[EllipticE[NDArray[{0.3, 0.6, 0.9}], 0.5]] "
                "== EllipticE[{0.3, 0.6, 0.9}, 0.5]");
    /* a packed List (Range-produced) agrees with the literal List */
    assert_true("EllipticK[N[Range[1, 3]/4]] == EllipticK[{0.25, 0.5, 0.75}]");
    /* the buffer element must equal the scalar path to machine precision */
    assert_true("Abs[First[Normal[EllipticF[NDArray[{0.6}], 0.5]]] "
                "- N[EllipticF[3/5, 1/2]]] < 10^-14");
    /* m > 1: the machine Carlson kernel declines (the value is complex), the
     * array is abandoned and the List path answers through FLINT. Both
     * representations must still agree. */
    assert_true("Normal[EllipticK[NDArray[{1.5}]]] == EllipticK[{1.5}]");
}

/* ---- attributes and arity ------------------------------------------- */

static void test_attributes(void) {
    assert_eval("Attributes[EllipticK]", "{Listable, NumericFunction, Protected}");
    assert_eval("Attributes[EllipticF]", "{Listable, NumericFunction, Protected}");
    assert_eval("Attributes[EllipticE]", "{Listable, NumericFunction, Protected}");
    assert_eval("Attributes[EllipticPi]", "{Listable, NumericFunction, Protected}");
}

static void test_arity(void) {
    /* Wrong arity stays unevaluated, with an argx / argt diagnostic to stderr. */
    assert_eval("EllipticK[1, 2]", "EllipticK[1, 2]");
    assert_eval("EllipticF[1]", "EllipticF[1]");
    assert_eval("EllipticE[]", "EllipticE[]");
    assert_eval("EllipticE[1, 2, 3]", "EllipticE[1, 2, 3]");
    assert_eval("EllipticPi[1]", "EllipticPi[1]");
    assert_eval("EllipticPi[1, 2, 3, 4]", "EllipticPi[1, 2, 3, 4]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_exact_reductions);
    TEST(test_exact_stays_symbolic);
    TEST(test_complete);
    TEST(test_incomplete);
    TEST(test_complex_branch);
    TEST(test_principal_value);
    TEST(test_derivative_forms);
    TEST(test_derivative_residuals);
    TEST(test_corpus_residuals);
    TEST(test_nd_surfaces);
    TEST(test_attributes);
    TEST(test_arity);

    printf("All elliptic integral tests passed.\n");
    return 0;
}
