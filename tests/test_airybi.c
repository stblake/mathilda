/* Tests for AiryBi: the Airy function Bi(z) and its minimal derivative head
 * AiryBiPrime.
 *
 * Covers the exact value at 0 and limits at +-Infinity; machine- and
 * arbitrary-precision (MPFR) real and complex numerics across the Maclaurin
 * (small |z|), dominant-asymptotic (large |z|, |arg z| < pi/3), and
 * connection-to-Ai (large |z|, incl. the oscillatory negative-real-axis)
 * regimes; an independent oracle built from the exact DLMF 9.2.10 relation
 *   Bi(z) = e^{i pi/6} Ai(z e^{2 pi i/3}) + e^{-i pi/6} Ai(z e^{-2 pi i/3})
 * which ties Bi to the (mpfr_ai-verified) AiryAi builtin across all sectors
 * -- MPFR provides no Bi, so this is the strongest available numerical check;
 * magnitude-overflow promotion for large arguments; list threading; the
 * derivative rules and Taylor series at 0; the asymptotic Series at Infinity;
 * argument-count diagnostics; and attributes. */

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

/* Oracle: the exact DLMF 9.2.10 connection relation expresses Bi entirely in
 * terms of Ai *values* (no derivatives), so the relative residual
 *   |AiryBi[z] - (e^{i pi/6} AiryAi[z e^{2pi i/3}] + e^{-i pi/6} AiryAi[z e^{-2pi i/3}])|
 *   / (1 + |AiryBi[z]|)
 * must vanish to working precision for every z. AiryAi is itself pinned to
 * mpfr_ai, so this independently verifies AiryBi in every argument sector. */
static void assert_xcheck(const char* zexpr, double tol) {
    /* Use an L1 metric on the real/imag parts so the result is guaranteed to be
     * a machine Real (Abs of a machine Complex can retain a tiny imaginary
     * residue here, which is irrelevant to the cross-check). */
    char buf[768];
    snprintf(buf, sizeof buf,
        "(Abs[Re[AiryBi[%s] - (Exp[I Pi/6] AiryAi[(%s) Exp[2 Pi I/3]] "
        "+ Exp[-I Pi/6] AiryAi[(%s) Exp[-2 Pi I/3]])]] + "
        "Abs[Im[AiryBi[%s] - (Exp[I Pi/6] AiryAi[(%s) Exp[2 Pi I/3]] "
        "+ Exp[-I Pi/6] AiryAi[(%s) Exp[-2 Pi I/3]])]]) / "
        "(1 + Abs[Re[AiryBi[%s]]] + Abs[Im[AiryBi[%s]]])",
        zexpr, zexpr, zexpr, zexpr, zexpr, zexpr, zexpr, zexpr);
    double v = eval_real(buf);
    ASSERT_MSG(v <= tol, "AiryBi cross-check at %s: relative residual %.3e > %.3e",
               zexpr, v, tol);
}

/* ---- exact value, limits -------------------------------------------- */

void test_airybi_exact() {
    /* AiryBi[0] = 1/(3^(1/6) Gamma[2/3]). */
    assert_eval_eq("AiryBi[0]", "1/(3^(1/6) Gamma[2/3])", 0);
    /* Bi grows to Infinity at +Infinity, decays to 0 at -Infinity. */
    assert_eval_eq("AiryBi[Infinity]", "Infinity", 0);
    assert_eval_eq("AiryBi[-Infinity]", "0", 0);
    assert_eval_eq("AiryBi[Indeterminate]", "Indeterminate", 0);
    /* Exact non-zero arguments stay symbolic (numericalize only under N). */
    assert_eval_eq("AiryBi[2]", "AiryBi[2]", 0);
    assert_eval_eq("AiryBi[x]", "AiryBi[x]", 0);
}

/* ---- machine-precision real ----------------------------------------- */

void test_airybi_machine_real() {
    assert_close("AiryBi[0.0]",  0.6149266274460007, 1e-7);  /* Bi(0) exact */
    assert_close("AiryBi[1.2]",  1.421133, 1e-5);
    assert_close("AiryBi[1.5]",  1.878941, 1e-5);
    assert_close("AiryBi[1.8]",  2.595869, 1e-5);
    /* Negative axis: oscillatory, exercises the connection-to-Ai path. */
    assert_close("AiryBi[-2.0]", -0.4123026, 1e-6);
    assert_close("AiryBi[-5.0]", -0.1383691, 1e-6);
    /* Large positive: deep in the dominant asymptotic regime (huge but exact). */
    assert_close("AiryBi[40.0]", 3.953139e72, 1e66);
}

/* ---- arbitrary-precision real --------------------------------------- */

void test_airybi_arbitrary_real() {
    /* N[AiryBi[2], 50] -- matches the reference to ~49 digits. */
    assert_eval_startswith("N[AiryBi[2], 50]",
        "3.2980949999782147102806044252234524220039759634");
    /* Precision tracks the input (`30 -> ~29 digits). */
    assert_eval_startswith("AiryBi[1.8`30]",
        "2.5958693567439062900598912548");
    /* High precision via N on an exact rational. */
    assert_eval_startswith("N[AiryBi[1/2], 60]",
        "0.854277043103155493300048798795243180856787400504739626934");
}

/* ---- Ai-connection oracle across all sectors ------------------------ */
/* The relation holds for every z; sweep the real axis (Maclaurin + dominant
 * + connection regimes, both signs) and a spread of complex points. */
void test_airybi_oracle() {
    for (int t = -36; t <= 36; t++) {
        double x = t * 0.5;
        char z[32];
        snprintf(z, sizeof z, "%.1f", x);
        assert_xcheck(z, 1e-8);
    }
    const char* cpts[] = {
        "2.0 + 3.0 I", "-30.0 + 5.0 I", "8.0 I", "-8.0 + 0.5 I",
        "5.0 - 2.0 I", "-15.0 - 4.0 I", "1.0 + 1.0 I", "20.0 I",
        "0.3 + 0.2 I", "-1.0 - 1.0 I"
    };
    for (size_t i = 0; i < sizeof cpts / sizeof cpts[0]; i++)
        assert_xcheck(cpts[i], 1e-8);
}

/* ---- machine and arbitrary-precision complex ------------------------ */

void test_airybi_complex() {
    /* AiryBi[2.5 + I] = 0.512544 + 5.335 I (spec). */
    assert_close("Re[AiryBi[2.5 + I]]", 0.512544, 1e-5);
    assert_close("Im[AiryBi[2.5 + I]]", 5.335, 1e-3);
    /* High precision tracks for a moderate complex argument (cross-checked by
     * the oracle above; here we just confirm the MPFR path produces digits). */
    assert_eval_startswith("Head[N[AiryBi[1 + I], 30]]", "Complex");
}

/* ---- magnitude-overflow promotion ----------------------------------- */
/* AiryBi[1000.0] ~ 5.4e9154 overflows the double range, so the machine-input
 * result is promoted to MPFR rather than becoming Inf (MachineNumberQ False). */
void test_airybi_overflow_promotion() {
    Expr* e = parse_expression("AiryBi[1000.0]");
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strstr(s, "inf") == NULL && strstr(s, "Inf") == NULL,
               "AiryBi[1000.0] overflowed to %s", s);
    ASSERT_MSG(strstr(s, "e+9154") != NULL,
               "AiryBi[1000.0] expected ~5e9154 magnitude, got %s", s);
    free(s);
    expr_free(r);
}

/* ---- list threading (Listable) -------------------------------------- */

void test_airybi_listable() {
    assert_eval_eq("AiryBi[{a, b}]", "{AiryBi[a], AiryBi[b]}", 0);
    assert_eval_eq("AiryBi[{}]", "{}", 0);
    assert_close("AiryBi[{1.2, 1.8}][[1]]", 1.421133, 1e-5);
    assert_close("AiryBi[{1.2, 1.8}][[2]]", 2.595869, 1e-5);
}

/* ---- derivatives and Taylor series at 0 ----------------------------- */

void test_airybi_derivatives() {
    assert_eval_eq("D[AiryBi[x], x]", "AiryBiPrime[x]", 0);
    /* Higher derivatives reduce via Bi'' = x Bi. */
    assert_eval_eq("D[AiryBi[x], {x, 2}]", "x AiryBi[x]", 0);
    assert_eval_eq("D[AiryBiPrime[x], x]", "x AiryBi[x]", 0);
    /* AiryBiPrime[0] closed form (needed by the Taylor series). */
    assert_eval_eq("AiryBiPrime[0]", "3^(1/6)/Gamma[1/3]", 0);
    /* Taylor series at 0 (closed-form coefficients). */
    assert_eval_eq("Series[AiryBi[x], {x, 0, 4}]",
        "1/(3^(1/6) Gamma[2/3]) + 3^(1/6)/Gamma[1/3] x + "
        "1/6/(3^(1/6) Gamma[2/3]) x^3 + 1/4/(3^(5/6) Gamma[1/3]) x^4 + O[x]^5", 0);
}

/* ---- asymptotic series at Infinity ---------------------------------- */

void test_airybi_series_infinity() {
    assert_eval_eq("Series[AiryBi[x], {x, Infinity, 2}]",
        "E^(2/3 x^(3/2)) ((1/x)^(1/4)/Sqrt[Pi] + "
        "5/48/Sqrt[Pi] (1/x)^(7/4) + O[1/x]^(13/4))", 0);
    assert_eval_eq("Series[AiryBi[x], {x, Infinity, 1}]",
        "((1/x)^(1/4)/Sqrt[Pi] + O[1/x]^(7/4)) E^(2/3 x^(3/2))", 0);
}

/* ---- inert AiryBiPrime ---------------------------------------------- */

void test_airybiprime_inert() {
    /* No general numeric evaluator: stays unevaluated away from 0. */
    assert_eval_eq("N[AiryBiPrime[2.0]]", "AiryBiPrime[2.0]", 0);
    assert_eval_eq("AiryBiPrime[y]", "AiryBiPrime[y]", 0);
}

/* ---- argument-count diagnostics ------------------------------------- */

void test_airybi_argcount() {
    /* Wrong arity stays unevaluated (a diagnostic is printed to stderr). */
    assert_eval_eq("AiryBi[]", "AiryBi[]", 0);
    assert_eval_eq("AiryBi[1, 2, 3]", "AiryBi[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_airybi_attributes() {
    SymbolDef* d = symtab_get_def("AiryBi");
    ASSERT_MSG(d != NULL, "AiryBi not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "AiryBi not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "AiryBi not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "AiryBi not Protected");
    ASSERT_MSG((d->attributes & ATTR_READPROTECTED) != 0, "AiryBi not ReadProtected");

    SymbolDef* p = symtab_get_def("AiryBiPrime");
    ASSERT_MSG(p != NULL, "AiryBiPrime not registered");
    ASSERT_MSG((p->attributes & ATTR_PROTECTED) != 0, "AiryBiPrime not Protected");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_airybi_exact);
    TEST(test_airybi_machine_real);
    TEST(test_airybi_arbitrary_real);
    TEST(test_airybi_oracle);
    TEST(test_airybi_complex);
    TEST(test_airybi_overflow_promotion);
    TEST(test_airybi_listable);
    TEST(test_airybi_derivatives);
    TEST(test_airybi_series_infinity);
    TEST(test_airybiprime_inert);
    TEST(test_airybi_argcount);
    TEST(test_airybi_attributes);

    printf("All AiryBi tests passed.\n");
    return 0;
}
