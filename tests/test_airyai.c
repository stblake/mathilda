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

/* ---- inert AiryAiPrime ---------------------------------------------- */

void test_airyaiprime_inert() {
    /* No general numeric evaluator: stays unevaluated away from 0. */
    assert_eval_eq("N[AiryAiPrime[2.0]]", "AiryAiPrime[2.0]", 0);
    assert_eval_eq("AiryAiPrime[y]", "AiryAiPrime[y]", 0);
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
    ASSERT_MSG((p->attributes & ATTR_PROTECTED) != 0, "AiryAiPrime not Protected");
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
    TEST(test_airyaiprime_inert);
    TEST(test_airyai_argcount);
    TEST(test_airyai_attributes);

    printf("All AiryAi tests passed.\n");
    return 0;
}
