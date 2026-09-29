/*
 * test_dsolve_m60_stress.c — anti-overfit stress families for M60
 * (DSolve`GeneralizedAiry: the n-th order pure-power potential y^(n) == A x^m y,
 *  and DSolve`OperatorFactor at order 2 + its adjoint/left-factor peel).
 *
 * Every family is a FORWARD GENERATOR: the equation is built from parameters whose
 * closed form is guaranteed by construction, so passing cannot be an accident of a
 * hand-picked example.
 *
 *   F1  pure power      y^(n) == a x^m y, n = 3,4 over an (m, a) grid.  The solution
 *                       set is x^j 0F_{n-1}(...) and is verified by NUMERIC back-
 *                       substitution, since the pFq residual is undecidable by
 *                       zero_test.
 *   F2  gauge           D[x^k y, {x, n}] == a x^(m+k) y, i.e. the pure-power equation
 *                       conjugated by x^(-k) (k = 1 is corpus 2.1.2-241/604).  This
 *                       exercises the depression pre-pass, not the recogniser alone.
 *   F3  symbolic        y''' == a x^b y with SYMBOLIC a and b: the corpus shape
 *                       (2.1.2-595).  Solved symbolically, then instantiated and
 *                       verified -- the M16 trap, where a symbolic-parameter residual
 *                       numericizes to all-NaN and the correct branch is discarded.
 *   F4  order-2 factor  (D - r2)(D - r1) over a rational r grid: the answer must be a
 *                       CLOSED FORM, not a SeriesData.  Before M60, OperatorFactor was
 *                       gated to order >= 3, so a rational first-order right factor
 *                       that Kovacic declined fell through to a truncated Frobenius
 *                       series.
 *   F5  order-3 factor  (D - r3)(D - r2)(D - r1): three constants, residual ~ 0.
 *   F6  adjoint         L = (D + s) o Q with Q IRREDUCIBLE (Airy, Kovacic case 2):
 *                       DSolve`DFactor must report the order-2 right factor and the
 *                       first-order LEFT factor, and the factorization is RECONSTRUCTED
 *                       and compared against the original operator -- so this checks
 *                       the adjoint and left-division algebra, not just a shape.
 *   F7  forced          F4/F5 with polynomial forcing carried through the peel.
 *   F8  declines        non-power potential, degenerate (logarithmic) exponent pattern,
 *                       nonlinear: Head =!= List, bounded.
 *
 * Each family asserts Head[sol] === List FIRST, so a declining method cannot pass
 * vacuously.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}
static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}
#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

/* Solve `eqn` and require `resid` to back-substitute to a numeric ~0 at several
 * sample points, with the first `nc` generated constants instantiated at distinct
 * non-round rationals.  `head` is "DSolve" or a pinned DSolve`<Method>. */
static void solve_ok_x(const char* head, const char* eqn, const char* resid, int nc,
                       bool closed_form) {
    char cs[192] = "";
    const char* vals[] = { "6/10", "13/10", "7/10", "9/10", "11/10", "3/10" };
    for (int k = 1; k <= nc; k++) {
        char one[48];
        snprintf(one, sizeof(one), "%sC[%d] -> %s", (k > 1 ? ", " : ""), k, vals[k - 1]);
        strncat(cs, one, sizeof(cs) - strlen(cs) - 1);
    }
    char buf[1700];
    snprintf(buf, sizeof(buf),
        "With[{sol = %s[%s, y, x]}, Head[sol] === List && Length[sol] >= 1 && %s"
        "Module[{r = (%s) /. sol[[1]] /. {%s}}, "
        "  Max[Table[Abs[N[r /. x -> xv, 25]], {xv, {6/10, 11/10, 17/10}}]] < 10^-8]]",
        head, eqn, (closed_form ? "FreeQ[sol, SeriesData] && " : ""), resid, cs);
    ASSERT_TRUE(buf);
}
static void solve_ok(const char* head, const char* eqn, const char* resid, int nc) {
    solve_ok_x(head, eqn, resid, nc, false);
}
/* As solve_ok, but ALSO requires a genuine closed form (no truncated SeriesData) —
 * checked inside the same solve so the equation is not solved twice. */
static void solve_closed_ok(const char* eqn, const char* resid, int nc) {
    solve_ok_x("DSolve", eqn, resid, nc, true);
}

/* ---- F1: the pure-power potential y^(n) == a x^m y ---------------------- */
static void t_m60_pure_power(void) {
    const int ms3[] = { 1, 2, 3 };
    const int as3[] = { 1, -2 };
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 2; j++) {
            char eqn[160], res[160];
            snprintf(eqn, sizeof(eqn), "y'''[x] == (%d) x^%d y[x]", as3[j], ms3[i]);
            snprintf(res, sizeof(res), "y'''[x] - (%d) x^%d y[x]", as3[j], ms3[i]);
            solve_ok("DSolve", eqn, res, 3);
        }
    /* order 4 */
    for (int m = 1; m <= 2; m++) {
        char eqn[160], res[160];
        snprintf(eqn, sizeof(eqn), "y''''[x] == 3 x^%d y[x]", m);
        snprintf(res, sizeof(res), "y''''[x] - 3 x^%d y[x]", m);
        solve_ok("DSolve", eqn, res, 4);
    }
    /* the pinned method owns m == 0 too (p = n != 0), where the automatic cascade
     * hands the constant-coefficient case to LinearConstantCoefficients first */
    solve_ok("DSolve`GeneralizedAiry", "y'''[x] + y[x] == 0", "y'''[x] + y[x]", 3);
}

/* ---- F2: the depression gauge, D[x^k y, {x,n}] == a x^(m+k) y ------------ */
static void t_m60_gauge_family(void) {
    for (int k = 1; k <= 3; k++)
        for (int m = 1; m <= 2; m++) {
            char eqn[200], res[200];
            snprintf(eqn, sizeof(eqn), "D[x^%d y[x], {x, 3}] == 2 x^%d y[x]", k, m + k);
            snprintf(res, sizeof(res), "D[x^%d y[x], {x, 3}] - 2 x^%d y[x]", k, m + k);
            solve_ok("DSolve", eqn, res, 3);
        }
    /* 2.1.2-241 / -604 exactly: k = 1, n = 3, m = 1 expands to x y''' + 3 y'' == a x^2 y */
    solve_ok("DSolve", "x y'''[x] + 3 y''[x] - x^2 y[x] == 0",
                       "x y'''[x] + 3 y''[x] - x^2 y[x]", 3);
}

/* ---- F3: symbolic coefficient AND symbolic exponent (2.1.2-595) --------- */
static void t_m60_symbolic_exponent(void) {
    ASSERT_TRUE("With[{s = DSolve[y'''[x] - a x^b y[x] == 0, y, x]}, "
                "Head[s] === List && Length[s] >= 1 && !FreeQ[s, HypergeometricPFQ]]");
    /* instantiate the parameters and verify numerically on the original equation */
    ASSERT_TRUE("With[{s = DSolve[y'''[x] - a x^b y[x] == 0, y, x]}, "
                "Head[s] === List && Module[{f = (y /. s[[1]]) /. {a -> 3/2, b -> 5/2}}, "
                "Max[Table[Abs[N[(D[f[x], {x, 3}] - (3/2) x^(5/2) f[x]) /. "
                "  {C[1] -> 6/10, C[2] -> 13/10, C[3] -> 7/10} /. x -> xv, 25]], "
                "  {xv, {6/10, 11/10, 17/10}}]] < 10^-8]]");
}

/* ---- F4: order-2 operator factoring must beat the series fallback ------- */
static void t_m60_order2_factorable(void) {
    /* L = (D - r2)(D - r1), r rational: written as the composition so the generator,
     * not the author, decides the expanded coefficients. */
    const char* r1[] = { "2/x",   "1/x",      "3/x",       "1/(1+x)" };
    const char* r2[] = { "1/x",   "2/(1+x)",  "1/x + 1",   "2/x"     };
    for (size_t i = 0; i < 4; i++) {
        char eqn[320], res[320];
        snprintf(eqn, sizeof(eqn),
            "D[y'[x] - (%s) y[x], x] - (%s) (y'[x] - (%s) y[x]) == 0", r1[i], r2[i], r1[i]);
        snprintf(res, sizeof(res),
            "D[y'[x] - (%s) y[x], x] - (%s) (y'[x] - (%s) y[x])", r1[i], r2[i], r1[i]);
        /* the answer must also be a genuine closed form, not a truncated series */
        solve_closed_ok(eqn, res, 2);
    }
    /* the exact 2.1.2-253 quotient that used to come back as O[x]^6 */
    solve_closed_ok("x^2 (1+x) y''[x] + 2 x (2+x) y'[x] + 2 y[x] == 0",
                    "x^2 (1+x) y''[x] + 2 x (2+x) y'[x] + 2 y[x]", 2);
}

/* ---- F5: order-3 fully factorable operators ----------------------------- */
static void t_m60_order3_factorable(void) {
    const char* r1[] = { "2/x", "1/x",     "2/x"     };
    const char* r2[] = { "1/x", "2/(1+x)", "1/x + 1" };
    const char* r3[] = { "3/x", "1/x",     "1/x"     };
    for (size_t i = 0; i < 3; i++) {
        char inner[200], eqn[560], res[560];
        snprintf(inner, sizeof(inner),
            "(D[y'[x] - (%s) y[x], x] - (%s) (y'[x] - (%s) y[x]))", r1[i], r2[i], r1[i]);
        snprintf(eqn, sizeof(eqn), "D[%s, x] - (%s) %s == 0", inner, r3[i], inner);
        snprintf(res, sizeof(res), "D[%s, x] - (%s) %s", inner, r3[i], inner);
        solve_ok("DSolve", eqn, res, 3);
    }
    /* 2.1.2-250 and -253 themselves */
    solve_ok("DSolve", "-12 y[x] + 3 (2 x^2+1) y''[x] + x (x^2+1) y'''[x] == 0",
                       "-12 y[x] + 3 (2 x^2+1) y''[x] + x (x^2+1) y'''[x]", 3);
    solve_ok("DSolve",
        "-4 (1+3 x) y[x] + 2 x (2+5 x) y'[x] - 2 x^2 (2 x+1) y''[x] + x^3 (x+1) y'''[x] == 0",
        "-4 (1+3 x) y[x] + 2 x (2+5 x) y'[x] - 2 x^2 (2 x+1) y''[x] + x^3 (x+1) y'''[x]", 3);
}

/* ---- F6: the adjoint / left-factor peel (Beke order-(n-1) right factors) -
 * L = (D + s) o Q with Q irreducible, so no first-order RIGHT factor exists and the
 * only way in is through the adjoint.  The returned factorization is reconstructed
 * and compared against the original operator. */
static void t_m60_adjoint_left_factor(void) {
    const char* qs[] = {
        "D[y[x], {x, 2}] - x y[x]",                             /* Airy: irreducible   */
        "D[y[x], {x, 2}] - (x/4 + 5/(16 x^2)) y[x]",            /* Kovacic case 2      */
        "D[y[x], {x, 2}] - (1 + 1/x) y[x]",                     /* irreducible         */
    };
    const char* ss[] = { "1/x", "2/x", "1/(1+x)" };
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++) {
            char buf[1100];
            snprintf(buf, sizeof(buf),
                "Module[{lhs = D[%s, x] + (%s) (%s), f, q, sf, qy}, "
                "  f = DSolve`DFactor[lhs == 0, y[x], x]; "
                "  Head[f] === List && Length[f] == 2 && !FreeQ[f[[1]], Dx^2] && "
                "  FreeQ[f[[2]], Dx^2] && "
                "  (q = f[[1]]; sf = Cancel[Together[f[[2]] - Dx]]; "
                "   qy = Sum[Coefficient[q, Dx, k] D[y[x], {x, k}], {k, 0, 2}]; "
                "   TrueQ[Cancel[Together[D[qy, x] + sf qy - lhs]] == 0])]",
                qs[i], ss[j], qs[i]);
            ASSERT_TRUE(buf);
        }
    /* an operator with NO first-order factor either way stays irreducible (one inert
     * remainder), so the adjoint search cannot invent a factorization */
    ASSERT_TRUE("Module[{f = DSolve`DFactor[y'''[x] + x y[x] == 0, y[x], x]}, "
                "Head[f] === List && Length[f] == 1 && !FreeQ[f[[1]], Dx^3]]");
}

/* ---- F7: forcing carried through the peel ------------------------------- */
static void t_m60_forced(void) {
    const char* g[] = { "x", "x^2 + 1", "1" };
    for (size_t i = 0; i < 3; i++) {
        char eqn[380], res[380];
        snprintf(eqn, sizeof(eqn),
            "D[y'[x] - (2/x) y[x], x] - (1/x) (y'[x] - (2/x) y[x]) == %s", g[i]);
        snprintf(res, sizeof(res),
            "D[y'[x] - (2/x) y[x], x] - (1/x) (y'[x] - (2/x) y[x]) - (%s)", g[i]);
        solve_ok("DSolve", eqn, res, 2);
    }
    /* order 3, forced */
    solve_ok("DSolve",
        "-12 y[x] + 3 (2 x^2+1) y''[x] + x (x^2+1) y'''[x] == x",
        "-12 y[x] + 3 (2 x^2+1) y''[x] + x (x^2+1) y'''[x] - x", 3);
}

/* ---- F8: bounded declines (never a wrong answer) ------------------------ */
static void t_m60_declines(void) {
    /* non-power potential */
    ASSERT_TRUE("Head[DSolve`GeneralizedAiry[y'''[x] - Sin[x] y[x] == 0, y, x]] "
                "=!= List");
    /* degenerate lower parameter: exponents differ by a multiple of p (logarithmic
     * Frobenius case) -- 2.1.2-1204 and the 4th-order x^2 y'''' == a y */
    ASSERT_TRUE("Head[DSolve`GeneralizedAiry[x y'''[x] + y[x] == 0, y, x]] =!= List");
    ASSERT_TRUE("Head[DSolve`GeneralizedAiry[x^2 y''''[x] - 3 y[x] == 0, y, x]] =!= List");
    /* order too low for the recogniser (the 2nd-order rows own these) */
    ASSERT_TRUE("Head[DSolve`GeneralizedAiry[y''[x] - x y[x] == 0, y, x]] =!= List");
    /* nonlinear */
    ASSERT_TRUE("Head[DSolve`GeneralizedAiry[y'''[x] == y[x]^2, y, x]] =!= List");
    ASSERT_TRUE("Head[DSolve`OperatorFactor[y'''[x] == y[x]^2, y, x]] =!= List");
    /* an irreducible operator: OperatorFactor declines both ways */
    ASSERT_TRUE("Head[DSolve`OperatorFactor[y''[x] - x y[x] == 0, y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m60_pure_power);
    TEST(t_m60_gauge_family);
    TEST(t_m60_symbolic_exponent);
    TEST(t_m60_order2_factorable);
    TEST(t_m60_order3_factorable);
    TEST(t_m60_adjoint_left_factor);
    TEST(t_m60_forced);
    TEST(t_m60_declines);

    printf("All DSolve M60 stress tests passed.\n");
    return 0;
}
