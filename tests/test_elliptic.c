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

    assert_eval("EllipticPi[0, m]", "EllipticK[m]");
    assert_eval("EllipticPi[0, phi, m]", "EllipticF[phi, m]");
    assert_eval("EllipticPi[n, 0, m]", "0");
    assert_eval("EllipticPi[n, Pi/2, m]", "EllipticPi[n, m]");
}

/* E(phi | 1) = Int_0^phi |Cos[t]| dt, which is Sin[phi] ONLY on the principal
 * strip |phi| <= Pi/2. This file used to assert the unconditional
 * `EllipticE[phi, 1] -> Sin[phi]`, which made EllipticE[2, 1] answer
 * Sin[2] = 0.909297 where the true value is 2 - Sin[2] = 1.090703 -- and
 * jumped 0.18 against its own neighbour at m = 1 - 10^-18. The rule now fires
 * only where it is provably valid, so a symbolic amplitude stays symbolic
 * rather than wrong, and the numeric path (which continues correctly off the
 * strip) answers when asked. */
static void test_m_one_continuation(void) {
    assert_eval("EllipticE[1/3, 1]", "Sin[1/3]");      /* inside the strip */
    assert_eval("EllipticE[Pi/2, 1]", "1");            /* at the edge */
    assert_eval("EllipticE[phi, 1]", "EllipticE[phi, 1]");   /* undecidable */
    assert_eval("EllipticE[2, 1]", "EllipticE[2, 1]"); /* outside: not Sin[2] */
    /* and off the strip the numeric answer is the continuation, not Sin */
    assert_true("Abs[N[EllipticE[2, 1], 20] - (2 - Sin[2])] < 10^-18");
    assert_true("Abs[N[EllipticE[2, 1], 20] - N[EllipticE[2, 1 - 10^-25], 20]] < 10^-20");
}

/* Poles answer ComplexInfinity, in every spelling of the argument.
 * EllipticK[SetPrecision[1, 30]] used to come back as an unevaluated
 * EllipticK[1.0]: ell_is_one did not look at EXPR_MPFR, so the pole check was
 * skipped, Arb was handed an exact 1, and its non-finite ball became NULL.
 * EllipticPi[1., 1/2] was unevaluated for the same reason, with no pole check at
 * all. The INCOMPLETE third kind has no pole at n = 1 and must keep answering. */
static void test_poles(void) {
    assert_eval("EllipticK[1]", "ComplexInfinity");
    assert_eval("EllipticK[1.]", "ComplexInfinity");
    assert_eval("EllipticK[SetPrecision[1, 30]]", "ComplexInfinity");
    assert_eval("EllipticPi[1, 1/2]", "ComplexInfinity");
    assert_eval("EllipticPi[1., 0.5]", "ComplexInfinity");
    /* Pi[1-e | 1/2] grows as 1/Sqrt[e], so the complete form really diverges */
    assert_true("N[EllipticPi[1 - 1/10^8, 1/2], 20] > 22000");
    /* but the incomplete form at n = 1 is finite (quadrature: 1.73199154202) */
    assert_true("Abs[N[EllipticPi[1, 1, 1/2], 20] - 1.73199154202] < 10^-10");
}

/* An inexact argument must give an inexact answer. These reductions used to run
 * BEFORE the numeric path, so EllipticK[0.] returned the exact Pi/2 (head
 * Times) and -- worst -- EllipticPi[0., 1/2] returned the SYMBOLIC
 * EllipticK[1/2], i.e. a numeric call came back un-numeric. The rest of the
 * system is consistent about this (Sin[0.] is 0., Gamma[1.] is 1.). */
static void test_inexact_contagion(void) {
    assert_true("Head[EllipticK[0.]] === Real");
    assert_true("Head[EllipticE[0.]] === Real");
    assert_true("Head[EllipticE[1.]] === Real");
    assert_true("Head[EllipticE[0., 0.5]] === Real");
    assert_true("Head[EllipticPi[0., 1/2]] === Real");
    assert_true("Abs[EllipticK[0.] - N[Pi/2]] < 10^-15");
    assert_true("Abs[EllipticPi[0., 1/2] - N[EllipticK[1/2], 20]] < 10^-15");
    /* the exact spellings are untouched */
    assert_eval("EllipticK[0]", "1/2 Pi");
    assert_eval("EllipticE[1]", "1");
    /* and a reduction Arb cannot read (the exact Pi/2 limit) still fires, then
     * evaluates -- the numeric-first ordering must not strand it */
    assert_true("Abs[EllipticE[Pi/2, 0.5] - N[EllipticE[1/2], 20]] < 10^-15");
    assert_true("Abs[EllipticF[Pi/2, 0.5] - N[EllipticK[1/2], 20]] < 10^-15");
}

/* An exact, non-special argument stays symbolic, as in the Wolfram Language --
 * only inexact input or an explicit N[...] evaluates. */
static void test_exact_stays_symbolic(void) {
    /* 1/2 and -1 are NOT in this group: both are singular values with a closed
     * form, asserted in test_closed_forms below. 1/3 has none. */
    assert_eval("EllipticK[1/3]", "EllipticK[1/3]");
    assert_eval("EllipticF[1/3, 1/2]", "EllipticF[1/3, 1/2]");
    assert_eval("EllipticE[1/3, 1/2]", "EllipticE[1/3, 1/2]");
    assert_eval("EllipticPi[1/3, 1/2]", "EllipticPi[1/3, 1/2]");
    assert_eval("EllipticPi[1/3, 1/4, 1/2]", "EllipticPi[1/3, 1/4, 1/2]");
    assert_eval("EllipticK[m]", "EllipticK[m]");
}

/* The exact values Mathematica gives and this file did not. Each closed form is
 * pinned BOTH as a printed form (so a change is deliberate) and numerically
 * against the 30-digit value, so a transcription slip in the expression cannot
 * pass. The numeric references are mpmath's. */
static void test_closed_forms(void) {
    /* K's two singular values. E has no closed form at either, and must not
     * acquire a guessed one. */
    assert_true("Abs[N[EllipticK[-1], 30] "
                "- 1.31102877714605990523241979494556] < 10^-28");
    assert_true("Abs[N[EllipticK[1/2], 30] "
                "- 1.85407467730137191843385034719526] < 10^-28");
    assert_true("FreeQ[EllipticK[-1], EllipticK]");
    assert_true("FreeQ[EllipticK[1/2], EllipticK]");
    assert_eval("EllipticE[1/2]", "EllipticE[1/2]");
    assert_eval("EllipticE[-1]", "EllipticE[-1]");

    /* the third kind at m = 0, both arities */
    assert_eval("EllipticPi[n, 0]", "(1/2 Pi)/Sqrt[1 - n]");
    assert_true("Abs[N[EllipticPi[3/10, 0], 25] - N[Pi/(2 Sqrt[1 - 3/10]), 25]] < 10^-23");
    assert_eval("EllipticPi[n, phi, 0]", "ArcTanh[Sqrt[-1 + n] Tan[phi]]/Sqrt[-1 + n]");

    /* the value at infinity, in every argument slot */
    assert_eval("EllipticK[Infinity]", "0");
    assert_eval("EllipticE[Infinity]", "ComplexInfinity");
    assert_eval("EllipticE[ComplexInfinity]", "ComplexInfinity");
    assert_eval("EllipticF[phi, Infinity]", "0");
    assert_eval("EllipticE[Infinity, Infinity]", "ComplexInfinity");
    assert_eval("EllipticPi[Infinity, m]", "0");
    assert_eval("EllipticPi[n, Infinity]", "0");

    /* odd in the amplitude: each integrand is even in t */
    assert_eval("EllipticF[-phi, m]", "-EllipticF[phi, m]");
    assert_eval("EllipticE[-phi, m]", "-EllipticE[phi, m]");
    assert_eval("EllipticPi[n, -phi, m]", "-EllipticPi[n, phi, m]");
    assert_eval("EllipticF[-2 x, m]", "-EllipticF[2 x, m]");
    /* superficial negativity only, as the trig heads do: -x-y does not fold */
    assert_eval("EllipticF[-x - y, m]", "EllipticF[-x - y, m]");
    /* and the rule is right, not just tidy */
    assert_true("Abs[N[EllipticF[-3/10, 1/2], 20] + N[EllipticF[3/10, 1/2], 20]] < 10^-18");
    assert_true("Abs[N[EllipticPi[1/3, -3/10, 1/2], 20] "
                "+ N[EllipticPi[1/3, 3/10, 1/2], 20]] < 10^-18");
}

/* Series at m = 0 for the complete integrals.
 *
 * This needs a DEDICATED kernel and cannot be had any other way: the generic
 * Taylor-via-D path evaluates d/dm K at m = 0, which is (Pi/2 - Pi/2)/0, i.e.
 * Times[0, ComplexInfinity] -> Indeterminate, and series.c's has_infinity then
 * abandons the whole expansion -- which is why Series[EllipticK[m], {m, 0, 1}]
 * did not previously emit even the leading Pi/2. Coefficients are
 * a_0 = 1, a_k = a_{k-1} ((2k-1)/(2k))^2, with K = (Pi/2) sum a_k m^k and
 * E = (Pi/2) sum a_k m^k/(1-2k); both match Mathematica term for term. */
static void test_series_at_zero(void) {
    assert_eval("Series[EllipticK[m], {m, 0, 3}]",
                "1/2 Pi + 1/8 Pi m + 9/128 Pi m^2 + 25/512 Pi m^3 + O[m]^4");
    assert_eval("Series[EllipticE[m], {m, 0, 3}]",
                "1/2 Pi + -1/8 Pi m + -3/128 Pi m^2 + -5/512 Pi m^3 + O[m]^4");
    /* the leading term alone, which used to come back unevaluated */
    assert_eval("Series[EllipticK[m], {m, 0, 1}]", "1/2 Pi + 1/8 Pi m + O[m]^2");
    /* composition through an inner series */
    assert_eval("Series[EllipticK[2 x^2], {x, 0, 4}]",
                "1/2 Pi + 1/4 Pi x^2 + 9/32 Pi x^4 + O[x]^5");
    /* an expansion about a regular point must still take Taylor-via-D */
    assert_true("FreeQ[Series[EllipticK[m], {m, 1/2, 2}], Series]");
    /* and the incomplete forms, which already worked, are untouched */
    assert_eval("Series[EllipticF[phi, m], {phi, 0, 3}]",
                "phi + 1/6 m phi^3 + O[phi]^4");
    assert_eval("Series[EllipticE[phi, m], {phi, 0, 3}]",
                "phi + -1/6 m phi^3 + O[phi]^4");
    /* the series really is the function: truncation agrees numerically */
    /* a five-term truncation at m = 1/10; the next term is (Pi/2)(3969/65536)
     * m^5 = 9.5*10^-7, so 10^-5 is the honest bound for this many terms */
    assert_true("Abs[(Pi/2)(1 + 1/4 (1/10) + 9/64 (1/10)^2 + 25/256 (1/10)^3 "
                "+ 1225/16384 (1/10)^4) - N[EllipticK[1/10], 20]] < 10^-5");
}

/* Interval[] threading. Three different mechanisms are at work and all three
 * are pinned, because each failed differently before:
 *
 *  - the COMPLETE forms need bespoke monotone rows (K increasing, E decreasing
 *    below m = 1): their derivatives reintroduce themselves, so the generic
 *    certifier can never bottom out;
 *  - the amplitude slot of F / E / Pi needed Power[Interval, NEGATIVE rational]
 *    to thread, since the phi-derivative is built as Power[rad, -1/2] -- that
 *    gap is why EllipticE[Interval, m] worked and EllipticF[Interval, m] did
 *    not, and it was a defect in power.c, not here;
 *  - the m slot declines, and must do so CLEANLY. With the self-referential
 *    derivative now in closed form, the certifier's descent would otherwise
 *    leave the evaluator unable to reach a fixed point, and the user saw
 *    `$IterationLimit exceeded`. */
static void test_interval(void) {
    assert_eval("EllipticK[Interval[{0.4, 0.5}]]", "Interval[{1.77752, 1.85407}]");
    assert_eval("EllipticE[Interval[{0.4, 0.5}]]", "Interval[{1.35064, 1.39939}]");
    /* amplitude slot */
    assert_true("Head[EllipticF[Interval[{1/4, 1/2}], 1/2]] === Interval");
    assert_true("Head[EllipticE[Interval[{1/4, 1/2}], 1/2]] === Interval");
    assert_true("Head[EllipticPi[1/2, Interval[{1/4, 1/2}], 1/2]] === Interval");
    /* the general Power fix that unblocked them */
    assert_eval("Power[Interval[{1/4, 1/2}], -1/2]", "Interval[{Sqrt[2], 2}]");
    assert_eval("Power[Interval[{-1, 1}], -1/2]", "1/Sqrt[Interval[{-1, 1}]]");
    /* beyond m = 1 the complete forms are complex: must NOT thread */
    assert_eval("EllipticK[Interval[{1.2, 1.5}]]", "EllipticK[Interval[{1.2, 1.5}]]");
    assert_eval("EllipticK[Interval[{0.5, 1.5}]]", "EllipticK[Interval[{0.5, 1.5}]]");
    /* m slot: a clean symbolic decline, not a non-terminating evaluation */
    assert_eval("EllipticF[0.2, Interval[{0.3, 0.4}]]",
                "EllipticF[0.2, Interval[{0.3, 0.4}]]");
    assert_eval("EllipticPi[Interval[{0.2, 0.3}], 0.5]",
                "EllipticPi[Interval[{0.2, 0.3}], 0.5]");
    /* the enclosures are sound: K is increasing, so the endpoints are the bounds */
    assert_true("Abs[First[First[EllipticK[Interval[{0.4, 0.5}]]]] - N[EllipticK[2/5]]] < 10^-5");
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
    /* The m-derivative of F and both derivatives of the COMPLETE Pi are now
     * closed forms, each checked against a central difference (below). The
     * n- and m-derivatives of the INCOMPLETE Pi stay deliberately inert: a fit
     * over the natural candidate basis does not recover them, so guessing is
     * exactly the failure mode deriv.c's note warns about. */
    assert_true("FreeQ[D[EllipticF[p, m], m], Derivative]");
    assert_true("FreeQ[D[EllipticPi[n, m], n], Derivative]");
    assert_true("FreeQ[D[EllipticPi[n, m], m], Derivative]");
    assert_eval("D[EllipticPi[n, p, m], n]", "Derivative[1, 0, 0][EllipticPi][n, p, m]");
    assert_eval("D[EllipticPi[n, p, m], m]", "Derivative[0, 0, 1][EllipticPi][n, p, m]");
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
        /* the three rules that were inert placeholders until they were checked */
        "EllipticF[3/5, mm]",           "mm", "2/5",
        "EllipticPi[nn, 1/2]",          "nn", "3/10",
        "EllipticPi[3/10, mm]",         "mm", "1/2",
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

    /* EllipticPi has no machine kernel, so it is not packed_aware and a VISIBLE
     * NDArray arrives at the builtin untouched. All three argument positions
     * used to leave it UNEVALUATED, which CLAUDE.md counts as a wrong answer;
     * the builtin now delists and re-evaluates. */
    assert_true("Normal[EllipticPi[NDArray[{0.2, 0.5}], 0.5]] "
                "== EllipticPi[{0.2, 0.5}, 0.5]");
    assert_true("Normal[EllipticPi[0.5, NDArray[{0.2, 0.5}]]] "
                "== EllipticPi[0.5, {0.2, 0.5}]");
    assert_true("Normal[EllipticPi[0.5, NDArray[{0.2, 0.5}], 0.25]] "
                "== EllipticPi[0.5, {0.2, 0.5}, 0.25]");
    assert_true("Head[EllipticPi[NDArray[{0.2, 0.5}], 0.5]] =!= EllipticPi");

    /* The binary kernels decline where 1 - m Sin[phi]^2 < 0 as well, not only
     * on the unary m > 1 route; that fall-through had no coverage. */
    assert_true("Normal[EllipticF[NDArray[{1.4}], 2.5]] == EllipticF[{1.4}, 2.5]");
    assert_true("Normal[EllipticE[NDArray[{1.4}], 2.5]] == EllipticE[{1.4}, 2.5]");
}

/* EllipticPi had NO machine kernel at all: ~40 us per element against SciPy's
 * 316 ns for the same Carlson composition, and `Compiled -> False` at both
 * arities. The standing reason was that n > 1 needs a Cauchy principal value --
 * but it does not: the value there is genuinely COMPLEX
 * (Pi[3/2 | 1/2] = -0.456720313453 - 2.72069904635 I), so what a double kernel
 * owes is a DECLINE, exactly as K/E/F already decline outside their real
 * domains. With Carlson R_J (and R_C, which R_J needs) the complete form runs at
 * 41 ns/elt and 2 ulp.
 *
 * Note what must NOT break: `packed_aware` is a property of the SYMBOL, not of
 * one arity, so registering the two-argument kernel stops the transparency gate
 * materialising packed Lists for the THREE-argument form as well. The ND
 * element-wise layer tops out at arity 2, so that form has no buffer path and
 * relies on the builtin's delist fallback -- which the last two assertions pin. */
static void test_ellipticpi_kernel(void) {
    /* dyadic grid, so no input skew: this measures the kernel */
    assert_true("Max[Abs[Normal[EllipticPi[NDArray[N[Range[1, 255]/256]], 0.5]] "
                "/ Table[N[EllipticPi[k/256, 1/2], 30], {k, 1, 255}] - 1]] < 2*10^-15");
    assert_true("Max[Abs[Normal[EllipticPi[0.5, NDArray[N[Range[1, 255]/256]], 0.25]] "
                "/ Table[N[EllipticPi[1/2, k/256, 1/4], 30], {k, 1, 255}] - 1]] < 2*10^-15");
    /* n >= 1: the kernel declines, the buffer is abandoned, Arb answers complex,
     * and the two representations still agree */
    assert_true("Abs[EllipticPi[1.5, 0.5] - N[EllipticPi[3/2, 1/2], 20]] < 10^-15");
    assert_true("Normal[EllipticPi[NDArray[{1.5, 0.5}], 0.5]] "
                "== EllipticPi[{1.5, 0.5}, 0.5]");
    /* the quasi-period, which the incomplete kernel applies off the strip */
    assert_true("Abs[N[EllipticPi[1/3, 7/10 + Pi, 1/2], 25] "
                "- (N[EllipticPi[1/3, 7/10, 1/2], 25] + 2 N[EllipticPi[1/3, 1/2], 25])] "
                "< 10^-22");
    /* a packed List at arity 3 must still answer, element for element, even
     * though the symbol is now packed_aware and the gate no longer materialises */
    assert_true("EllipticPi[0.5, N[Range[1, 300]/250], 0.25] "
                "== Table[EllipticPi[0.5, N[k/250], 0.25], {k, 1, 300}]");
    assert_true("Head[EllipticPi[0.5, N[Range[1, 300]/250], 0.25]] === List");
}

/* Compile[] lowering is driven generically off the ND kernel registry, so these
 * are not merely compilability checks -- they pin that the compiled VALUE agrees
 * with the interpreter, which nothing asserted for any elliptic head before. The
 * three-argument EllipticPi lowers through OP_KERNN from the N-ary kernel, which
 * the NDArray element-wise layer itself cannot use. */
static void test_compile_values(void) {
    /* Part[.., 1, 2] is the RHS of the first rule, i.e. the Compiled flag. The
     * obvious `[[1]] === (Compiled -> True)` does NOT work: the diagnostics
     * carry their own `Compiled` symbol and a literal one written here does not
     * compare equal to it. */
    assert_true("Part[CompileDiagnostics[{{x, _Real}}, EllipticK[x]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{v, _Real, 1}}, EllipticK[v]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{x, _Real}}, EllipticE[x]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{v, _Real, 1}}, EllipticE[v]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{x, _Real}}, EllipticF[x, 0.5]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{v, _Real, 1}}, EllipticF[v, 0.5]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{x, _Real}}, EllipticPi[0.5, x]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{v, _Real, 1}}, EllipticPi[0.5, v]], 1, 2] === True");
    assert_true("Part[CompileDiagnostics[{{x, _Real}}, EllipticPi[0.5, x, 0.25]], 1, 2] === True");

    assert_true("Abs[Compile[{{x, _Real}}, EllipticK[x]][0.3] - EllipticK[0.3]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticE[x]][0.3] - EllipticE[0.3]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticF[x, 0.5]][0.7] - EllipticF[0.7, 0.5]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticE[x, 0.5]][0.7] - EllipticE[0.7, 0.5]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticPi[0.5, x]][0.3] - EllipticPi[0.5, 0.3]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticPi[0.5, x, 0.25]][0.7] "
                "- EllipticPi[0.5, 0.7, 0.25]] < 10^-15");
    /* rank-1 agrees elementwise with the interpreter */
    assert_true("Compile[{{v, _Real, 1}}, EllipticK[v]][{0.2, 0.5, 0.8}] "
                "== EllipticK[{0.2, 0.5, 0.8}]");
    /* Outside the machine domain (m > 1, where the kernel declines because the
     * value is complex) the compiled function does NOT return NaN: the VM bails
     * and the interpreter answers, so the compiled and interpreted values are
     * the SAME complex number. Pinned because it is the behaviour a caller
     * depends on and the opposite of what a "machine-domain only" reading of
     * the compile contract would predict. */
    assert_true("Head[Compile[{{x, _Real}}, EllipticK[x]][1.5]] === Complex");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticK[x]][1.5] - EllipticK[1.5]] < 10^-15");
    assert_true("Abs[Compile[{{x, _Real}}, EllipticPi[1.5, x]][0.5] "
                "- EllipticPi[1.5, 0.5]] < 10^-15");
}

/* The buffer and the scalar path must agree to MACHINE precision, not merely to
 * the 1.42e-14 relative tolerance that `Equal` applies to machine reals -- which
 * is why the assertions above passed while the buffer was 61-106 ulp out.
 *
 * Two separate defects were hiding under that tolerance. The duplication
 * stopped at a hardcoded 0.01 (with a `#define EC_RTOL 1e-16` beside it that was
 * never read), buying ~1e-12 where the comment claimed double precision; and
 * `a = 1 - Sin[r]^2` instead of `Cos[r]^2` lost EIGHT significant digits near
 * r = Pi/2, where the subtraction's absolute error is the size of the answer. */
static void test_buffer_matches_scalar(void) {
    /* The grid is k/256, and the denominator being a POWER OF TWO is the whole
     * point: k/256 is exactly representable, so `N[...]` on the buffer side and
     * the exact Rational on the reference side are the same number, and the
     * comparison measures the kernel. On a k/200 grid they are not: the double
     * differs from the decimal by ~1e-17, dK/dm reaches 100 at k = 199, and the
     * resulting 6.7e-16 of pure INPUT SKEW looks exactly like kernel error.
     * (Measured: 4.4e-16 on k/256 against 6.7e-16 on k/200, same kernel.)
     * Bounds below are the measured worst case with headroom -- 2, 7 and 2.5
     * ulp -- and they still catch the 61/106/70 ulp this used to deliver. */
    assert_true("Max[Abs[Normal[EllipticK[NDArray[N[Range[1, 255]/256]]]] "
                "/ Table[N[EllipticK[k/256], 30], {k, 1, 255}] - 1]] < 2*10^-15");
    assert_true("Max[Abs[Normal[EllipticE[NDArray[N[Range[1, 255]/256]]]] "
                "/ Table[N[EllipticE[k/256], 30], {k, 1, 255}] - 1]] < 4*10^-15");
    assert_true("Max[Abs[Normal[EllipticF[NDArray[N[Range[1, 255]/256]], 0.5]] "
                "/ Table[N[EllipticF[k/256, 1/2], 30], {k, 1, 255}] - 1]] < 2*10^-15");
    /* the Pi/2 approach, where 1 - Sin^2 used to be wrong by 2.7*10^-8 */
    assert_true("Max[Table[Abs[First[Normal[EllipticF[NDArray[{N[Pi/2] - 1.0*10^-k}], 0.99]]] "
                "/ N[EllipticF[Rationalize[N[Pi/2] - 1.0*10^-k, 0], 99/100], 30] - 1], "
                "{k, 2, 12}]] < 5*10^-15");
    assert_true("Max[Table[Abs[First[Normal[EllipticE[NDArray[{N[Pi/2] - 1.0*10^-k}], 0.99]]] "
                "/ N[EllipticE[Rationalize[N[Pi/2] - 1.0*10^-k, 0], 99/100], 30] - 1], "
                "{k, 2, 12}]] < 5*10^-15");
    /* Rationalize[.., 0] above is the exact value of the double, so these two
     * carry no input skew either. Worst measured: 1.7e-15 (F), 1.1e-15 (E),
     * against 2.7e-08 (F) and 9.8e-10 (E) before the Cos^2 fix. */
}

/* N[expr, d] must deliver d GOOD digits, and the two places it can lose them
 * are independent.
 *
 * The ARGUMENT's precision is chosen by numeric_plan_working_spec
 * (src/numeric.c), which used to raise the working precision only when an input
 * LEAF was lossy -- blind to how much the function itself amplifies, which
 * depends on conditioning, not on spelling. With an exact rational input, so
 * nothing was lossy and nothing was raised, N[EllipticK[1 - 10^-17], 20] gave
 * six good digits of twenty, and N[Zeta[1 + 10^-20], 20] gave 2^66 + 1.
 *
 * The KERNEL's own losses are a different question, answered by the accuracy
 * ladder in src/flint_num_bridge.c: EllipticPi[-10^30, 1/2] at 30 digits asks
 * Arb for 100 bits and gets a ball holding 78, so the last six digits were
 * garbage before the ladder retried. */
static void test_precision_is_honest(void) {
    /* argument side: ill-conditioned, exact input, every digit must be right */
    assert_true("Abs[N[EllipticK[99999999999999999/100000000000000000], 20] "
                "- 20.95826765156927898288] < 10^-18");
    assert_true("Abs[N[EllipticF[2, 999999999999999999/1000000000000000000], 25] "
                "- 42.69566795256993004057617] < 10^-22");
    assert_true("Abs[N[Zeta[1 + 1/10^20], 40] "
                "- 100000000000000000000.5772156649015328606072] < 10^-18");
    /* kernel side: the ladder's own regression case */
    assert_true("Abs[N[EllipticPi[-10^30, 1/2], 30] "
                "/ (1570796326794897122662117945335/10^45) - 1] < 10^-29");
    /* a zero-valued result must not climb the ladder looking for relative
     * accuracy it can never have -- it must come straight back */
    assert_eval("FLINT`Zeta[-2.]", "0.0");
    /* and the well-conditioned path is unchanged */
    assert_true("Abs[N[EllipticK[2/5], 40] "
                "- 1.777519371491253323502990072871915020298] < 10^-39");
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
    TEST(test_closed_forms);
    TEST(test_series_at_zero);
    TEST(test_interval);
    TEST(test_m_one_continuation);
    TEST(test_poles);
    TEST(test_inexact_contagion);
    TEST(test_complete);
    TEST(test_incomplete);
    TEST(test_complex_branch);
    TEST(test_principal_value);
    TEST(test_derivative_forms);
    TEST(test_derivative_residuals);
    TEST(test_corpus_residuals);
    TEST(test_nd_surfaces);
    TEST(test_buffer_matches_scalar);
    TEST(test_ellipticpi_kernel);
    TEST(test_compile_values);
    TEST(test_precision_is_honest);
    TEST(test_attributes);
    TEST(test_arity);

    printf("All elliptic integral tests passed.\n");
    return 0;
}
