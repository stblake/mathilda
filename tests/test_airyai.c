/* Tests for AiryAi: the Airy function Ai(z) and its minimal derivative head
 * AiryAiPrime.
 *
 * Covers the exact value at 0 and limits at +-Infinity; machine- and
 * arbitrary-precision (MPFR) real and complex numerics across the Maclaurin
 * (small |z|) and asymptotic (large |z|, incl. the negative-real-axis
 * connection-formula sector) regimes; an independent mpfr_ai oracle on the
 * real axis; magnitude-overflow promotion for large imaginary arguments; list
 * threading; the derivative rules and Taylor series at 0; the asymptotic
 * Series at Infinity; argument-count diagnostics; and attributes. */

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

#ifdef USE_MPFR
#include <mpfr.h>
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
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* ---- exact value, limits -------------------------------------------- */

void test_airyai_exact() {
    /* AiryAi[0] = 1/(3^(2/3) Gamma[2/3]). */
    assert_eval_eq("AiryAi[0]", "1/(3^(2/3) Gamma[2/3])", 0);
    /* Ai decays to 0 at both real infinities. */
    assert_eval_eq("AiryAi[Infinity]", "0", 0);
    assert_eval_eq("AiryAi[-Infinity]", "0", 0);
    /* Exact non-zero arguments stay symbolic (numericalize only under N). */
    assert_eval_eq("AiryAi[2]", "AiryAi[2]", 0);
    assert_eval_eq("AiryAi[x]", "AiryAi[x]", 0);
}

/* ---- machine-precision real ----------------------------------------- */

void test_airyai_machine_real() {
    /* Spec-documented values (loose tol; the mpfr_ai oracle below pins these
     * tightly and independently). */
    assert_close("AiryAi[0.0]",  0.3550280538878172, 1e-7);  /* Ai(0) exact */
    assert_close("AiryAi[1.2]",  0.106126, 1e-5);
    assert_close("AiryAi[1.8]",  0.0470362, 1e-5);
    assert_close("AiryAi[-2.0]", 0.2274074282, 1e-7);
    /* Negative axis: oscillatory, exercises the asymptotic connection path.
     * These two values are cross-checked against the high-precision (Maclaurin)
     * evaluation of the same points. */
    assert_close("AiryAi[-50.0]", -0.1618814236123209, 1e-9);
    /* Large positive: deep in the asymptotic regime (tiny but accurate). */
    assert_close("AiryAi[40.0]", 6.365742658552915e-75, 1e-81);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_airyai_arbitrary_real() {
    /* N[AiryAi[2], 50] -- matches the reference to 50 digits. */
    assert_eval_startswith("N[AiryAi[2], 50]",
        "0.0349241304232743791353220807918076097610602");
    /* Precision tracks the input. */
    assert_eval_startswith("AiryAi[4.8`40]",
        "0.00017032552328643494849005822667939733412");
    assert_eval_startswith("N[AiryAi[1/2], 60]",
        "0.231693606480833489769125254509921739618386475");
}

/* ---- mpfr_ai oracle on the real axis -------------------------------- */
/* Independently compute Ai(x) with MPFR's mpfr_ai over a grid spanning both
 * the Maclaurin (|x| small) and asymptotic (|x| large, both signs) regimes,
 * and check the builtin agrees. This is a strong correctness check that does
 * not depend on any hand-tabulated values. */
void test_airyai_mpfr_oracle() {
#ifdef USE_MPFR
    mpfr_t ax, ref;
    mpfr_init2(ax, 160);
    mpfr_init2(ref, 160);
    /* Half-integer grid from -18 to 18: spans the Maclaurin regime (|x| small)
     * and the asymptotic / connection regime (|x| large, both signs), and
     * covers fractional points the spec uses (1.2, 1.5, 1.8, ...). */
    for (int t = -36; t <= 36; t++) {
        double x = t * 0.5;
        char in[64];
        snprintf(in, sizeof(in), "AiryAi[%.1f]", x);
        double got = eval_real(in);
        mpfr_set_d(ax, x, MPFR_RNDN);
        mpfr_ai(ref, ax, MPFR_RNDN);
        double want = mpfr_get_d(ref, MPFR_RNDN);
        ASSERT_MSG(fabs(got - want) <= 1e-10 * (1.0 + fabs(want)),
                   "AiryAi[%.1f]: builtin %.14g vs mpfr_ai %.14g", x, got, want);
    }
    mpfr_clear(ax);
    mpfr_clear(ref);
#endif
}

/* ---- machine and arbitrary-precision complex ------------------------ */

void test_airyai_complex() {
    /* Maclaurin regime: AiryAi[2.5 + I] = -0.00191209 - 0.0180329 I (spec). */
    assert_close("Re[AiryAi[2.5 + I]]", -0.00191209, 1e-6);
    assert_close("Im[AiryAi[2.5 + I]]", -0.0180329, 1e-6);
    /* Asymptotic connection regime (large negative real part): values
     * cross-checked against the high-precision (Maclaurin) evaluation. */
    assert_close("Re[AiryAi[-30.0 + 5.0 I]]", 64309778642.37744, 1.0);
    assert_close("Im[AiryAi[-30.0 + 5.0 I]]", 72694657726.70228, 1.0);
    /* High precision tracks: real part of a moderate complex argument. */
    assert_eval_startswith("Re[N[AiryAi[1 + I], 30]]", "0.0");
}

/* ---- magnitude-overflow promotion ----------------------------------- */
/* Large imaginary argument: |Ai(150 I)| ~ 1e374 overflows the double range,
 * so the machine-input result is promoted to MPFR rather than becoming Inf. */
void test_airyai_overflow_promotion() {
    Expr* e = parse_expression("AiryAi[150.0 I]");
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    /* Must be a finite Complex (contains 'e+374' magnitude, not 'inf'). */
    ASSERT_MSG(strstr(s, "inf") == NULL && strstr(s, "Inf") == NULL,
               "AiryAi[150.0 I] overflowed to %s", s);
    ASSERT_MSG(strstr(s, "e+374") != NULL,
               "AiryAi[150.0 I] expected ~1e374 magnitude, got %s", s);
    free(s);
    expr_free(r);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_airyai_listable() {
    assert_eval_eq("AiryAi[{a, b}]", "{AiryAi[a], AiryAi[b]}", 0);
    assert_eval_eq("AiryAi[{}]", "{}", 0);
    assert_close("AiryAi[{1.2, 1.8}][[1]]", 0.106126, 1e-5);
    assert_close("AiryAi[{1.2, 1.8}][[2]]", 0.0470362, 1e-5);
}

/* ---- derivatives and Taylor series at 0 ----------------------------- */

void test_airyai_derivatives() {
    assert_eval_eq("D[AiryAi[x], x]", "AiryAiPrime[x]", 0);
    /* Higher derivatives reduce via Ai'' = x Ai. */
    assert_eval_eq("D[AiryAi[x], {x, 2}]", "x AiryAi[x]", 0);
    assert_eval_eq("D[AiryAiPrime[x], x]", "x AiryAi[x]", 0);
    /* AiryAiPrime[0] closed form (needed by the Taylor series). */
    assert_eval_eq("AiryAiPrime[0]", "-1/(3^(1/3) Gamma[1/3])", 0);
    /* Taylor series at 0 (closed-form coefficients). */
    assert_eval_eq("Series[AiryAi[x], {x, 0, 4}]",
        "1/(3^(2/3) Gamma[2/3]) + -1/(3^(1/3) Gamma[1/3]) x + "
        "1/6/(3^(2/3) Gamma[2/3]) x^3 + -1/12/(3^(1/3) Gamma[1/3]) x^4 + O[x]^5", 0);
}

/* ---- asymptotic series at Infinity ---------------------------------- */

void test_airyai_series_infinity() {
    assert_eval_eq("Series[AiryAi[x], {x, Infinity, 2}]",
        "E^(-2/3 x^(3/2)) (1/2/Sqrt[Pi] (1/x)^(1/4) + "
        "-5/96/Sqrt[Pi] (1/x)^(7/4) + O[1/x]^(13/4))", 0);
    assert_eval_eq("Series[AiryAi[x], {x, Infinity, 1}]",
        "E^(-2/3 x^(3/2)) (1/2/Sqrt[Pi] (1/x)^(1/4) + O[1/x]^(7/4))", 0);
}

/* ==================================================================== */
/* AiryAiPrime[z] = Ai'(z): full numeric evaluator                       */
/* ==================================================================== */

/* ---- exact values and limits ---------------------------------------- */

void test_airyaiprime_special() {
    /* AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3]). */
    assert_eval_eq("AiryAiPrime[0]", "-1/(3^(1/3) Gamma[1/3])", 0);
    /* Ai' decays to 0 at +Infinity (recessive solution). */
    assert_eval_eq("AiryAiPrime[Infinity]", "0", 0);
    /* At -Infinity Ai' oscillates with growing amplitude: no limit, so the
     * call is deliberately left unevaluated (unlike AiryAi[-Infinity] = 0). */
    assert_eval_eq("AiryAiPrime[-Infinity]", "AiryAiPrime[-Infinity]", 0);
    /* Exact non-zero / symbolic arguments stay unevaluated. */
    assert_eval_eq("AiryAiPrime[2]", "AiryAiPrime[2]", 0);
    assert_eval_eq("AiryAiPrime[y]", "AiryAiPrime[y]", 0);
}

/* ---- machine-precision real ----------------------------------------- */

void test_airyaiprime_machine_real() {
    /* Spec value. */
    assert_close("AiryAiPrime[0.5]", -0.224911, 1e-5);
    /* Maclaurin regime, both signs. */
    assert_close("AiryAiPrime[1.2]", -0.132785, 1e-5);
    assert_close("AiryAiPrime[1.8]", -0.0685248, 1e-6);
    /* Negative axis: oscillatory, exercises the asymptotic connection path. */
    assert_close("AiryAiPrime[-10.0]", 0.99626504413, 1e-7);
    /* Large positive: deep asymptotic regime (tiny but accurate). */
    assert_close("AiryAiPrime[20.0]", -7.586391625748e-27, 1e-33);
}

/* ---- independent mpfr_ai finite-difference oracle ------------------- */
/* MPFR provides Ai (mpfr_ai) but not Ai'. Reconstruct Ai'(x) independently
 * via a 4th-order central difference of mpfr_ai evaluated at high precision,
 *   f'(x) ~ (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h),  h = 1e-2,
 * and check the builtin agrees across the Maclaurin and asymptotic (both
 * signs) regimes. This does not depend on any Mathilda code. */
void test_airyaiprime_mpfr_oracle() {
#ifdef USE_MPFR
    const double h = 1e-2;
    mpfr_t ax, f1, f2, f3, f4, deriv, tmp;
    mpfr_init2(ax, 200); mpfr_init2(f1, 200); mpfr_init2(f2, 200);
    mpfr_init2(f3, 200); mpfr_init2(f4, 200);
    mpfr_init2(deriv, 200); mpfr_init2(tmp, 200);
    for (int t = -24; t <= 24; t++) {
        double x = t * 0.5;
        char in[64];
        snprintf(in, sizeof(in), "AiryAiPrime[%.1f]", x);
        double got = eval_real(in);

        mpfr_set_d(ax, x + 2 * h, MPFR_RNDN); mpfr_ai(f1, ax, MPFR_RNDN);
        mpfr_set_d(ax, x + h,     MPFR_RNDN); mpfr_ai(f2, ax, MPFR_RNDN);
        mpfr_set_d(ax, x - h,     MPFR_RNDN); mpfr_ai(f3, ax, MPFR_RNDN);
        mpfr_set_d(ax, x - 2 * h, MPFR_RNDN); mpfr_ai(f4, ax, MPFR_RNDN);
        /* deriv = (-f1 + 8 f2 - 8 f3 + f4) / (12 h). */
        mpfr_mul_ui(tmp, f2, 8, MPFR_RNDN);
        mpfr_sub(deriv, tmp, f1, MPFR_RNDN);
        mpfr_mul_ui(tmp, f3, 8, MPFR_RNDN);
        mpfr_sub(deriv, deriv, tmp, MPFR_RNDN);
        mpfr_add(deriv, deriv, f4, MPFR_RNDN);
        mpfr_div_d(deriv, deriv, 12 * h, MPFR_RNDN);
        double want = mpfr_get_d(deriv, MPFR_RNDN);

        ASSERT_MSG(fabs(got - want) <= 1e-6 * (1.0 + fabs(want)),
                   "AiryAiPrime[%.1f]: builtin %.14g vs mpfr_ai FD %.14g",
                   x, got, want);
    }
    mpfr_clear(ax); mpfr_clear(f1); mpfr_clear(f2); mpfr_clear(f3);
    mpfr_clear(f4); mpfr_clear(deriv); mpfr_clear(tmp);
#endif
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_airyaiprime_arbitrary_real() {
    /* N[AiryAiPrime[5/2], 50] -- matches the reference. (Prefix stops short of
     * the last digit, which N leaves unrounded vs the reference's 50-digit
     * rounding: true value continues ...3813577 07..., rounded to ...35771.) */
    assert_eval_startswith("N[AiryAiPrime[5/2], 50]",
        "-0.0262508810359032303648954962972325094463178381357");
    /* Precision tracks the input (0.5 to 20 significant digits). */
    assert_eval_startswith("AiryAiPrime[0.5`20]",
        "-0.224910532664683893136");
    /* High precision on a rational argument. */
    assert_eval_startswith("N[AiryAiPrime[1/2], 60]",
        "-0.224910532664683893135996990328583214825029635610892837136285");
}

/* ---- machine and arbitrary-precision complex ------------------------ */

void test_airyaiprime_complex() {
    /* Maclaurin regime: AiryAiPrime[2.5 + I] = -0.00187921 + 0.0310276 I. */
    assert_close("Re[AiryAiPrime[2.5 + I]]", -0.00187921, 1e-6);
    assert_close("Im[AiryAiPrime[2.5 + I]]", 0.0310276, 1e-6);
    /* High precision tracks (30-digit complex). */
    assert_eval_startswith("Re[AiryAiPrime[2.5`30 + I]]",
        "-0.00187920860963515417533502004628");
    assert_eval_startswith("Im[AiryAiPrime[2.5`30 + I]]",
        "0.0310276242841243496825433750570");
}

/* ---- magnitude-overflow promotion ----------------------------------- */
/* Large imaginary argument: |Ai'(150 I)| ~ 1e376 overflows the double range,
 * so the machine-input result is promoted to MPFR rather than becoming Inf. */
void test_airyaiprime_overflow_promotion() {
    Expr* e = parse_expression("AiryAiPrime[150.0 I]");
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strstr(s, "inf") == NULL && strstr(s, "Inf") == NULL,
               "AiryAiPrime[150.0 I] overflowed to %s", s);
    ASSERT_MSG(strstr(s, "e+376") != NULL,
               "AiryAiPrime[150.0 I] expected ~1e376 magnitude, got %s", s);
    free(s);
    expr_free(r);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_airyaiprime_listable() {
    assert_eval_eq("AiryAiPrime[{a, b}]", "{AiryAiPrime[a], AiryAiPrime[b]}", 0);
    assert_eval_eq("AiryAiPrime[{}]", "{}", 0);
    assert_close("AiryAiPrime[{1.2, 1.8}][[1]]", -0.132785, 1e-5);
    assert_close("AiryAiPrime[{1.2, 1.8}][[2]]", -0.0685248, 1e-6);
}

/* ---- Series (Taylor at 0 and asymptotic at Infinity) ---------------- */

void test_airyaiprime_series() {
    /* Taylor at 0 (closed-form coefficients, via the generic Taylor-via-D
     * path on AiryAiPrime[0] and the Ai'' = z Ai derivative rules). */
    assert_eval_eq("Series[AiryAiPrime[x], {x, 0, 4}]",
        "-1/(3^(1/3) Gamma[1/3]) + 1/2/(3^(2/3) Gamma[2/3]) x^2 + "
        "-1/3/(3^(1/3) Gamma[1/3]) x^3 + O[x]^5", 0);
    /* Asymptotic at Infinity (DLMF 9.7.6): leading -1/(2 Sqrt[Pi]) x^(1/4),
     * coefficients -7/96/Sqrt[Pi], 455/9216/Sqrt[Pi]. */
    assert_eval_eq("Series[AiryAiPrime[x], {x, Infinity, 2}]",
        "E^(-2/3 x^(3/2)) ((-1/2/Sqrt[Pi])/(1/x)^(1/4) + "
        "-7/96/Sqrt[Pi] (1/x)^(5/4) + O[1/x]^(11/4))", 0);
    assert_eval_eq("Series[AiryAiPrime[x], {x, Infinity, 3}]",
        "E^(-2/3 x^(3/2)) ((-1/2/Sqrt[Pi])/(1/x)^(1/4) + "
        "-7/96/Sqrt[Pi] (1/x)^(5/4) + 455/9216/Sqrt[Pi] (1/x)^(11/4) + "
        "O[1/x]^(17/4))", 0);
}

/* ---- derivative rule ------------------------------------------------ */

void test_airyaiprime_derivative() {
    /* D[AiryAiPrime[x], x] = x AiryAi[x]  (from Ai'' = z Ai). */
    assert_eval_eq("D[AiryAiPrime[x], x]", "x AiryAi[x]", 0);
}

/* ---- argument-count diagnostics ------------------------------------- */

void test_airyaiprime_argcount() {
    assert_eval_eq("AiryAiPrime[]", "AiryAiPrime[]", 0);
    assert_eval_eq("AiryAiPrime[1, 2, 3]", "AiryAiPrime[1, 2, 3]", 0);
}

/* ---- argument-count diagnostics ------------------------------------- */

void test_airyai_argcount() {
    /* Wrong arity stays unevaluated (a diagnostic is printed to stderr). */
    assert_eval_eq("AiryAi[]", "AiryAi[]", 0);
    assert_eval_eq("AiryAi[1, 2, 3]", "AiryAi[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_airyai_attributes() {
    SymbolDef* d = symtab_get_def("AiryAi");
    ASSERT_MSG(d != NULL, "AiryAi not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "AiryAi not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "AiryAi not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "AiryAi not Protected");
    ASSERT_MSG((d->attributes & ATTR_READPROTECTED) != 0, "AiryAi not ReadProtected");

    SymbolDef* p = symtab_get_def("AiryAiPrime");
    ASSERT_MSG(p != NULL, "AiryAiPrime not registered");
    ASSERT_MSG((p->attributes & ATTR_LISTABLE) != 0, "AiryAiPrime not Listable");
    ASSERT_MSG((p->attributes & ATTR_NUMERICFUNCTION) != 0, "AiryAiPrime not NumericFunction");
    ASSERT_MSG((p->attributes & ATTR_PROTECTED) != 0, "AiryAiPrime not Protected");
    ASSERT_MSG((p->attributes & ATTR_READPROTECTED) != 0, "AiryAiPrime not ReadProtected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_airyai_exact);
    TEST(test_airyai_machine_real);
    TEST(test_airyai_arbitrary_real);
    TEST(test_airyai_mpfr_oracle);
    TEST(test_airyai_complex);
    TEST(test_airyai_overflow_promotion);
    TEST(test_airyai_listable);
    TEST(test_airyai_derivatives);
    TEST(test_airyai_series_infinity);
    TEST(test_airyai_argcount);
    TEST(test_airyai_attributes);

    /* AiryAiPrime: full numeric evaluator. */
    TEST(test_airyaiprime_special);
    TEST(test_airyaiprime_machine_real);
    TEST(test_airyaiprime_mpfr_oracle);
    TEST(test_airyaiprime_arbitrary_real);
    TEST(test_airyaiprime_complex);
    TEST(test_airyaiprime_overflow_promotion);
    TEST(test_airyaiprime_listable);
    TEST(test_airyaiprime_series);
    TEST(test_airyaiprime_derivative);
    TEST(test_airyaiprime_argcount);

    printf("All AiryAi tests passed.\n");
    return 0;
}
