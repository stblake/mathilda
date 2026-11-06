/* Tests for HypergeometricPFQ (src/hypergeopfq.c) and its Hypergeometric0F1/
 * 1F1/2F1 convenience heads, plus the z-derivative rule in deriv.c. */
#include "test_utils.h"
#include "core.h"
#include <stdio.h>

/* ---- structural rules ------------------------------------------------ */

static void test_origin(void) {
    /* pFq[a, b, 0] == 1 */
    assert_eval_eq("HypergeometricPFQ[{a1,a2,a3},{b1,b2,b3},0]", "1", 0);
    assert_eval_eq("HypergeometricPFQ[{},{},0]", "1", 0);
}

static void test_exact_stays_symbolic(void) {
    /* All-exact, non-reducible, non-terminating: stays unevaluated. */
    assert_eval_eq("HypergeometricPFQ[{1,1},{3,3,3},2]",
                   "HypergeometricPFQ[{1, 1}, {3, 3, 3}, 2]", 0);
}

static void test_cancellation(void) {
    assert_eval_eq("HypergeometricPFQ[{a,b,c},{a,d,e},z]",
                   "HypergeometricPFQ[{b, c}, {d, e}, z]", 0);
    /* A non-positive integer common element must NOT cancel -- it makes the
     * series terminate to a polynomial instead. */
    assert_eval_eq("HypergeometricPFQ[{-2,a,b,c},{-2,d,e},z]",
        "1 + (a b c z)/(d e) + (a (1 + a) b (1 + b) c (1 + c) z^2)/(2 d (1 + d) e (1 + e))",
        0);
}

static void test_termination(void) {
    /* Upper non-positive integer -> finite polynomial (symbolic z). */
    assert_eval_eq("HypergeometricPFQ[{-3,a},{b,c},x]",
        "1 - (3 a x)/(b c) + (3 a (1 + a) x^2)/(b (1 + b) c (1 + c)) - "
        "(a (1 + a) (2 + a) x^3)/(b (1 + b) (2 + b) c (1 + c) (2 + c))", 0);
    /* Degree-0 termination. */
    assert_eval_eq("HypergeometricPFQ[{0},{b},x]", "1", 0);
}

static void test_threading(void) {
    assert_eval_eq("HypergeometricPFQ[{1,2,3,4},{5,6,7},{0.1,0.3,0.5}]",
                   "{1.01164, 1.03627, 1.06296}", 0);
}

/* ---- symbolic reductions to elementary functions -------------------- */

static void test_reductions(void) {
    assert_eval_eq("HypergeometricPFQ[{},{},z]", "E^z", 0);            /* 0F0 */
    assert_eval_eq("HypergeometricPFQ[{a},{},z]", "(1 - z)^(-a)", 0);  /* 1F0 */
    assert_eval_eq("HypergeometricPFQ[{1},{},z]", "1/(1 - z)", 0);     /* geometric */
    assert_eval_eq("HypergeometricPFQ[{1,1},{2},z]", "-Log[1 - z]/z", 0); /* 2F1 */
    assert_eval_eq("HypergeometricPFQ[{},{1/2},z]", "Cosh[2 Sqrt[z]]", 0); /* 0F1 */
    assert_eval_eq("HypergeometricPFQ[{},{3/2},z]",
                   "Sinh[2 Sqrt[z]]/(2 Sqrt[z])", 0);                 /* 0F1 */
    /* 2F1 central-binomial / arcsin family (drives Sum[.../Binomial[2k,k]]). */
    assert_eval_eq("HypergeometricPFQ[{1,1},{1/2},z]",
                   "1/(1 - z) + (Sqrt[z] ArcSin[Sqrt[z]])/(1 - z)^(3/2)", 0);
    assert_eval_eq("HypergeometricPFQ[{1,1},{3/2},z]",
                   "ArcSin[Sqrt[z]]/(Sqrt[1 - z] Sqrt[z])", 0);
    assert_eval_eq("HypergeometricPFQ[{2,1},{3/2},z]",
                   "1/(2 (1 - z)) + ArcSin[Sqrt[z]]/(2 (1 - z)^(3/2) Sqrt[z])", 0);
    /* Very-well-poised Ramanujan 4F3 (central-binomial-cubed 1/Pi class). */
    assert_eval_eq("HypergeometricPFQ[{1/2,1/2,1/2,5/4},{1/4,1,1},-1]", "2/Pi", 0);
    /* A nearby but non-matching 4F3 must NOT reduce. */
    assert_eval_eq("HypergeometricPFQ[{1/2,1/2,1/2,7/4},{3/4,1,1},-1]",
                   "HypergeometricPFQ[{1/2, 1/2, 1/2, 7/4}, {3/4, 1, 1}, -1]", 0);
}

/* ---- named closed forms: error function, incomplete gamma, inverse trig,
 * Gauss. Each identity is checked to 30 digits against the direct series
 * (MPFR parameters keep the exact-parameter rewrite from firing, so the
 * reference side really is the series). ------------------------------- */

static void check30(const char* closed, const char* series) {
    char buf[2048];
    snprintf(buf, sizeof buf, "Abs[N[(%s) - (%s), 40]] < 10^-30 Max[1, Abs[N[%s, 40]]]",
             closed, series, closed);
    assert_eval_eq(buf, "True", 0);
}

static void test_named_closed_forms(void) {
    assert_eval_eq("HypergeometricPFQ[{1/2},{3/2},-z]", "(Sqrt[Pi] Erf[Sqrt[z]])/(2 Sqrt[z])", 0);
    assert_eval_eq("HypergeometricPFQ[{1/2},{3/2},z]", "(Sqrt[Pi] Erfi[Sqrt[z]])/(2 Sqrt[z])", 0);
    assert_eval_eq("HypergeometricPFQ[{-1/2},{1/2},z]", "E^z - Sqrt[Pi] Sqrt[z] Erfi[Sqrt[z]]", 0);
    assert_eval_eq("HypergeometricPFQ[{1},{3/2},-z]",
                   "(Sqrt[Pi] Erfi[Sqrt[z]])/(2 E^z Sqrt[z])", 0);
    assert_eval_eq("HypergeometricPFQ[{2},{3},z]", "(2 (E^z - (-1 + E^z)/z))/z", 0);
    /* a neither half-integer nor a positive integer: declined (the incomplete
     * gamma form has no numerics on its cut, the series does) */
    assert_eval_eq("HypergeometricPFQ[{1/3},{4/3},z]", "HypergeometricPFQ[{1/3}, {4/3}, z]", 0);
    assert_eval_eq("HypergeometricPFQ[{1},{4/3},z]", "HypergeometricPFQ[{1}, {4/3}, z]", 0);
    assert_eval_eq("HypergeometricPFQ[{1/2,1},{3/2},z]", "ArcTanh[Sqrt[z]]/Sqrt[z]", 0);
    assert_eval_eq("HypergeometricPFQ[{1/2,1},{3/2},-z]", "ArcTan[Sqrt[z]]/Sqrt[z]", 0);
    assert_eval_eq("HypergeometricPFQ[{1/2,1/2},{3/2},z]", "ArcSin[Sqrt[z]]/Sqrt[z]", 0);
    assert_eval_eq("HypergeometricPFQ[{1/2,1/2},{3/2},-z]", "ArcSinh[Sqrt[z]]/Sqrt[z]", 0);
    /* symbolic Gauss parameters stay, as in the reference system */
    assert_eval_eq("HypergeometricPFQ[{a,b},{c},1]", "HypergeometricPFQ[{a, b}, {c}, 1]", 0);
    assert_eval_eq("HypergeometricPFQ[{1/3,1/2},{2},1]", "(2 Gamma[7/6])/(Sqrt[Pi] Gamma[5/3])", 0);

    const char* pts[] = { "7/3", "-7/3", "3/10", "-3/10", "2/5 + 3/4 I" };
    for (size_t k = 0; k < sizeof pts / sizeof pts[0]; k++) {
        char c[512], s[512];
        const char* z = pts[k];
        const char* fam[][2] = {
            { "HypergeometricPFQ[{1/2},{3/2},%s]",   "HypergeometricPFQ[{N[1/2,50]},{N[3/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{-1/2},{1/2},%s]",  "HypergeometricPFQ[{N[-1/2,50]},{N[1/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{-5/2},{-3/2},%s]", "HypergeometricPFQ[{N[-5/2,50]},{N[-3/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{5/2},{7/2},%s]",   "HypergeometricPFQ[{N[5/2,50]},{N[7/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{3},{4},%s]",       "HypergeometricPFQ[{N[3,50]},{N[4,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{1},{3/2},%s]",     "HypergeometricPFQ[{N[1,50]},{N[3/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{1},{7/2},%s]",     "HypergeometricPFQ[{N[1,50]},{N[7/2,50]},N[%s,50]]" },
            { "HypergeometricPFQ[{1},{4},%s]",       "HypergeometricPFQ[{N[1,50]},{N[4,50]},N[%s,50]]" },
        };
        for (size_t f = 0; f < sizeof fam / sizeof fam[0]; f++) {
            snprintf(c, sizeof c, fam[f][0], z);
            snprintf(s, sizeof s, fam[f][1], z);
            check30(c, s);
        }
        /* the 2F1 forms inside the unit disc (the series' own domain) */
        if (k >= 2) {
            snprintf(c, sizeof c, "HypergeometricPFQ[{1/2,1},{3/2},%s]", z);
            snprintf(s, sizeof s, "HypergeometricPFQ[{N[1/2,50],N[1,50]},{N[3/2,50]},N[%s,50]]", z);
            check30(c, s);
            snprintf(c, sizeof c, "HypergeometricPFQ[{1/2,1/2},{3/2},%s]", z);
            snprintf(s, sizeof s, "HypergeometricPFQ[{N[1/2,50],N[1/2,50]},{N[3/2,50]},N[%s,50]]", z);
            check30(c, s);
        }
    }
    /* Gauss at z = 1, against reference-system values (35 digits) */
    check30("HypergeometricPFQ[{1/3,1/2},{2},1]", "1.1595952669639283657699920515700208819`40");
    check30("HypergeometricPFQ[{1/4,2/3},{3/2},1]", "1.3241485971871386627937369335599509558`40");
    check30("HypergeometricPFQ[{-1/2},{1/2},7/3]", "-2.9326738896023747980081444006243659817`40");
}

/* 2F1 at arbitrary precision outside the unit disc (off the cut). */
static void test_2f1_mpfr_continuation(void) {
    assert_eval_startswith("N[HypergeometricPFQ[{2/3,11/12},{5/3},-1],20]", "0.76419793513910566523");
    assert_eval_startswith("N[HypergeometricPFQ[{2/3,11/12},{5/3},-2],20]", "0.64424468650536985019");
}

/* ---- convenience heads ---------------------------------------------- */

static void test_convenience_heads(void) {
    /* 0F1[b,z] == pFq[{},{b},z];  with b==1/2 reduces to Cosh. */
    assert_eval_eq("Hypergeometric0F1[1/2, z]", "Cosh[2 Sqrt[z]]", 0);
    /* 1F1[a,b,z] == pFq[{a},{b},z]; 1F1[1,2,z] == (E^z - 1)/z form via series,
     * but here just confirm it routes to pFq and stays symbolic for symbolic z. */
    assert_eval_eq("Hypergeometric1F1[a, b, 0]", "1", 0);
    /* 2F1[1,1,2,z] == -Log[1-z]/z. */
    assert_eval_eq("Hypergeometric2F1[1, 1, 2, z]", "-Log[1 - z]/z", 0);
}

/* ---- numeric: machine, complex, MPFR -------------------------------- */

static void test_numeric_machine(void) {
    assert_eval_eq("HypergeometricPFQ[{1,1},{3,3,3},2.]", "1.07893", 0);
    /* 0F0 numeric agrees with Exp. */
    assert_eval_eq("HypergeometricPFQ[{},{},2.]", "7.38906", 0);
}

static void test_numeric_complex(void) {
    assert_eval_startswith("HypergeometricPFQ[{I,I,I},{2,2,2},-1.0 I]", "0.870032");
}

static void test_numeric_mpfr(void) {
    /* 50-digit value from the Wolfram reference. */
    assert_eval_startswith("N[HypergeometricPFQ[{1,1,1},{3/2,3/2,3/2},10],50]",
                           "530.191888273625904388559616854440877927");
    /* 35-digit request: prefix must match. */
    assert_eval_startswith("N[HypergeometricPFQ[{1,1},{3,3,3},2],35]",
                           "1.0789");
}

static void test_convergence_gate(void) {
    /* p == q+1 with |z| >= 1 : not summed directly (stays unevaluated)... */
    assert_eval_eq("HypergeometricPFQ[{1,2,3},{4,5},2.]",
                   "HypergeometricPFQ[{1, 2, 3}, {4, 5}, 2.0]", 0);
    /* ...except 2F1, which is continued analytically to the whole plane; on
     * the cut z > 1 it is the limit from below (mpmath:
     * -1705.757331635544 - 356.7499860560238 I). */
    assert_eval_startswith("HypergeometricPFQ[{10,10},{50},2.]", "-1705.76 - 356.75");
    /* p == q+1 with |z| < 1 : converges. */
    assert_eval_startswith("HypergeometricPFQ[{1,2},{3},0.5]", "1.54518");
}

/* ---- derivative ----------------------------------------------------- */

static void test_derivative(void) {
    assert_eval_eq("D[HypergeometricPFQ[{a1,a2},{b1,b2,b3},x],x]",
        "(a1 a2 HypergeometricPFQ[{1 + a1, 1 + a2}, {1 + b1, 1 + b2, 1 + b3}, x])/(b1 b2 b3)",
        0);
    /* Constant w.r.t. the differentiation variable -> 0. */
    assert_eval_eq("D[HypergeometricPFQ[{a},{b},y],x]", "0", 0);
}

int main(void) {
    core_init();
    TEST(test_origin);
    TEST(test_exact_stays_symbolic);
    TEST(test_cancellation);
    TEST(test_termination);
    TEST(test_threading);
    TEST(test_reductions);
    TEST(test_named_closed_forms);
    TEST(test_2f1_mpfr_continuation);
    TEST(test_convenience_heads);
    TEST(test_numeric_machine);
    TEST(test_numeric_complex);
    TEST(test_numeric_mpfr);
    TEST(test_convergence_gate);
    TEST(test_derivative);
    printf("All HypergeometricPFQ tests passed.\n");
    return 0;
}
