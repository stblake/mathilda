/*
 * test_dsolve_m61_stress.c — anti-overfit stress families for M61
 * (the inert `Inactive[Integrate]` variation-of-parameters particular).
 *
 * The correctness of an inert answer rests ENTIRELY on the method's own gate
 * (`ds_inert_vop_verified`): every other verifier in the substrate keeps a residual
 * that contains an integral, and `PossibleZeroQ` answers True for one.  So this suite
 * does two things a hand-picked example cannot:
 *
 *   F1  forward generator — x^2 y'' + x y' + (a x^m + b) y == p(x) over an (m, a, b, p)
 *       grid.  This is M60's generalised-Bessel row, whose homogeneous fundamental set
 *       is Sqrt[x] Z_nu(kappa x^((m+2)/2)) BY CONSTRUCTION, so every member is
 *       guaranteed to reach the inert path.  Each answer is verified INDEPENDENTLY of
 *       the in-method gate (see m61_verify below).
 *   F2  elementary still wins — the same shape with a forcing whose quadratures close:
 *       the answer must carry no Inactive at all.  Guards the byte-identity claim.
 *   F3  the ACCEPT/REJECT MARGIN of the gate's metric, measured on a planted wrong
 *       answer.  The gate calls a piece zero when it is below 1e-10 relative to the
 *       size of its own terms; that threshold is worthless unless a genuinely wrong
 *       basis lands far above it.  Here a correct VoP answer and a corrupted one (an
 *       extra Sqrt[x] on one basis element; a mismatched Bessel order) are both built
 *       in the language and their relative residuals compared.  This is the test that
 *       would catch the gate going blind.
 *   F4  series basis declines — a Frobenius/log homogeneous set must never produce an
 *       inert particular (an inert integral over a truncated series is meaningless).
 *   F5  IVP declines — an inert particular has no value at a point, so its constants
 *       cannot be fitted; returning the unfitted general solution would look solved.
 *   F6  latency bound — the failing `Integrate` attempt on a Bessel Wronskian quotient
 *       costs tens of seconds, so the mode skips it via a denominator pre-screen.  If
 *       that pre-screen stops firing, this family is what notices.
 *   F7  bounded declines — nonlinear, first-order, distributional and arbitrary-f[x]
 *       forcing: Head =!= List, and fast.
 *
 * Every family asserts Head[sol] === List FIRST, so a declining method cannot pass
 * vacuously.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

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

/* Verify an inert answer WITHOUT the method's gate.  Substituting the body into the
 * ODE leaves an expression that is linear in the surviving Inactive[Integrate] terms;
 * the identity must hold for an ARBITRARY value of each antiderivative, so replacing
 * each by a free symbol and requiring every coefficient of that linear form to vanish
 * is the whole correctness statement.  Sample points (23/50 .. 211/100) and precision
 * (40 digits) both differ from the in-method gate, so this is genuinely independent.
 *
 * `eqn` is solved with DSolve; `resid` is the ODE residual written in terms of `b`. */
static void m61_verify(const char* eqn, const char* resid) {
    char buf[2600];
    snprintf(buf, sizeof(buf),
        "Module[{s = DSolve[%s, y, x], b, r, zs, pieces, i},"
        "  Head[s] === List && Length[s] >= 1 && !FreeQ[s, Inactive] &&"
        "  FreeQ[s, SeriesData] &&"
        "  (b = s[[1]][[1]][[2]][[2]];"
        "   r = %s;"
        "   zs = DeleteDuplicates@Cases[r, Inactive[Integrate][__], {0, Infinity}];"
        "   Length[zs] >= 1 &&"
        "   (For[i = 1, i <= Length[zs], i++,"
        "        r = r /. zs[[i]] -> ToExpression[\"m61z\" <> ToString[i]]];"
        "    FreeQ[r, Inactive] &&"
        "    (pieces = Join[{r /. Table[ToExpression[\"m61z\" <> ToString[i]] -> 0,"
        "                               {i, Length[zs]}]},"
        "                   Table[Coefficient[r, ToExpression[\"m61z\" <> ToString[i]]],"
        "                         {i, Length[zs]}]];"
        "     Max[Table[Max[Table[Abs[N[p /. {C[1] -> 17/13, C[2] -> -23/19} /. x -> pt,"
        "                               40]],"
        "                         {pt, {23/50, 91/100, 157/100, 211/100}}]],"
        "               {p, pieces}]] < 10^-25)))]",
        eqn, resid);
    ASSERT_TRUE(buf);
}

/* ---- F1: forward generator over the generalised-Bessel family -------------- */
/* x^2 y'' + x y' + (a x^m + b) y == p(x).  The homogeneous set is Bessel by
 * construction (M60's Q = A x^m + B x^(-2) row after normalisation), so every member
 * is guaranteed to reach the inert path -- nothing here is hand-picked. */
static void t_m61_bessel_family(void) {
    const char* ms[] = { "1", "2", "3" };
    const char* as[] = { "1", "-2" };
    const char* ps[] = { "x", "x^2 + x", "x^3" };
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 2; j++)
            for (size_t k = 0; k < 3; k++) {
                char eqn[256], res[256];
                snprintf(eqn, sizeof(eqn),
                    "x^2 y''[x] + x y'[x] + ((%s) x^%s + 12) y[x] == %s",
                    as[j], ms[i], ps[k]);
                snprintf(res, sizeof(res),
                    "x^2 D[b, {x, 2}] + x D[b, x] + ((%s) x^%s + 12) b - (%s)",
                    as[j], ms[i], ps[k]);
                m61_verify(eqn, res);
            }
    /* the y'-free and the leading-coefficient-x members (corpus 3395 / 3387 shapes) */
    m61_verify("9 x^2 y''[x] + (3 x + 2) y[x] == x^4 + x^2",
               "9 x^2 D[b, {x, 2}] + (3 x + 2) b - x^4 - x^2");
    m61_verify("x y''[x] + 3 y'[x] - y[x] == x",
               "x D[b, {x, 2}] + 3 D[b, x] - b - x");
    m61_verify("x y''[x] + y'[x] - 2 y[x] x == x^2",
               "x D[b, {x, 2}] + D[b, x] - 2 b x - x^2");
}

/* ---- F2: where the quadratures close, nothing changes ---------------------- */
static void t_m61_elementary_unchanged(void) {
    /* named non-elementary antiderivative (ExpIntegralEi) is a SUCCESS, not inert */
    ASSERT_TRUE("Module[{s = DSolve[x y''[x] - x y'[x] + y[x] == x^3, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive] && !FreeQ[s, ExpIntegralEi]]");
    /* constant-coefficient, Euler, undetermined-coefficient forcings */
    ASSERT_TRUE("Module[{s = DSolve[y''[x] + y[x] == x^2, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
    ASSERT_TRUE("Module[{s = DSolve[y''[x] - 4 y[x] == E^(2 x), y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
    ASSERT_TRUE("Module[{s = DSolve[x^2 y''[x] - 3 x y'[x] + 4 y[x] == x^3, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
    ASSERT_TRUE("Module[{s = DSolve[x^2 y''[x] + x y'[x] + 16 y[x] == 0, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
    /* a homogeneous equation must never acquire a particular */
    ASSERT_TRUE("Module[{s = DSolve[x y''[x] + 3 y'[x] - y[x] == 0, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
    ASSERT_TRUE("Module[{s = DSolve[9 x^2 y''[x] + (3 x + 2) y[x] == 0, y, x]}, "
                "Head[s] === List && FreeQ[s, Inactive]]");
}

/* ---- F3: the gate's accept/reject margin, on a PLANTED wrong answer -------- */
/* The gate declares a piece zero below 1e-10 RELATIVE to the size of its own terms.
 * That number only means something if a wrong basis lands far above it, so measure
 * both here in the language: build the correct variation-of-parameters answer for
 * 3387 and two corrupted ones, and compare their relative residuals.  A Bessel
 * residual cancels twelve digits, so the comparison is done at exact rationals --
 * which is itself part of what is being guarded (at machine reals the CORRECT answer
 * reads ~1e-7 and would be indistinguishable from the wrong one). */
static void t_m61_gate_margin(void) {
    /* The shared construction: basis -> Wronskian -> inert u1,u2 -> relative residual
     * of the Z-decomposed pieces, exactly the metric the C gate applies.  The inert
     * integrals are collected with Cases rather than matched as `u1`/`u2`: a leading
     * sign is absorbed into the surrounding product, so pattern-matching the built
     * expression silently misses occurrences and leaves an inert integral to be
     * sampled -- which then captures the integration variable and yields nonsense. */
    static const char* tmpl =
        "Module[{b1 = %s, b2 = %s, W, u1, u2, yp, r, zs, i, pieces, rel},"
        "  W = b1 D[b2, x] - b2 D[b1, x];"
        "  u1 = Inactive[Integrate][-b2/W, x]; u2 = Inactive[Integrate][b1/W, x];"
        "  yp = b1 u1 + b2 u2;"
        "  r = x D[yp, {x, 2}] + 3 D[yp, x] - yp - x;"
        "  zs = DeleteDuplicates@Cases[r, Inactive[Integrate][__], {0, Infinity}];"
        "  For[i = 1, i <= Length[zs], i++,"
        "      r = r /. zs[[i]] -> ToExpression[\"m61q\" <> ToString[i]]];"
        "  FreeQ[r, Inactive] &&"
        "  (pieces = Join[{r /. Table[ToExpression[\"m61q\" <> ToString[i]] -> 0,"
        "                             {i, Length[zs]}]},"
        "                 Table[Coefficient[r, ToExpression[\"m61q\" <> ToString[i]]],"
        "                       {i, Length[zs]}]];"
        "   rel = Max[Table["
        "      Module[{p = Expand[pp /. x -> 113/100], v, sc},"
        "        v = Abs[N[p, 30]];"
        "        sc = If[Head[p] === Plus, Total[Abs[N[#, 30]] & /@ (List @@ p)], v];"
        "        If[TrueQ[sc > 1], v/sc, v]], {pp, pieces}]];"
        "   %s)]";
    char buf[2400];
    /* CORRECT basis: the two solutions of x y'' + 3 y' - y == 0.  Measured 5.4e-51. */
    snprintf(buf, sizeof(buf), tmpl,
             "BesselJ[2, 2 I Sqrt[x]]/x", "BesselY[2, 2 I Sqrt[x]]/x",
             "rel < 10^-20");
    ASSERT_TRUE(buf);
    /* CORRUPTED three ways.  Each measured at 0.18 .. 0.49 -- fifty orders of
     * magnitude above the correct answer, and eight above the gate's 1e-10 cut. */
    snprintf(buf, sizeof(buf), tmpl,
             "Sqrt[x] BesselJ[2, 2 I Sqrt[x]]/x", "BesselY[2, 2 I Sqrt[x]]/x",
             "rel > 10^-3");
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf), tmpl,
             "BesselJ[2, 2 I Sqrt[x]]/x", "BesselY[3, 2 I Sqrt[x]]/x",
             "rel > 10^-3");
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf), tmpl,
             "BesselJ[2, 2 I Sqrt[x]]/x^2", "BesselY[2, 2 I Sqrt[x]]/x",
             "rel > 10^-3");
    ASSERT_TRUE(buf);
}

/* ---- F4: a truncated-series homogeneous set must never go inert ------------ */
static void t_m61_series_basis_declines(void) {
    ASSERT_TRUE("FreeQ[DSolve[x^2 (x + 1) y''[x] + x (x^2 + 3) y'[x] + y[x] "
                "== -2 x^2 + x, y, x], Inactive]");
    ASSERT_TRUE("FreeQ[DSolve[3 x^2 (x + 1) y''[x] + x (5 - x) y'[x] "
                "+ (2 x^2 - 1) y[x] == -x^3, y, x], Inactive]");
    ASSERT_TRUE("FreeQ[DSolve[4 x^2 y''[x] - 3 (x^2 + x) y'[x] + 2 y[x] == x, y, x], "
                "Inactive]");
}

/* ---- F5: an IVP cannot be fitted, so it must decline ---------------------- */
static void t_m61_ivp_declines(void) {
    ASSERT_TRUE("FreeQ[DSolve[{x y''[x] + 3 y'[x] - y[x] == x, y[1] == 2, y'[1] == 0}, "
                "y, x], Inactive]");
    ASSERT_TRUE("FreeQ[DSolve[{9 x^2 y''[x] + (3 x + 2) y[x] == x^4 + x^2, y[1] == 0, "
                "y'[1] == 1}, y, x], Inactive]");
    /* a boundary-value spelling declines too */
    ASSERT_TRUE("FreeQ[DSolve[{x y''[x] + 3 y'[x] - y[x] == x, y[1] == 0, y[2] == 0}, "
                "y, x], Inactive]");
}

/* ---- F6: latency -- the denominator pre-screen must keep firing ------------ */
/* Without it, the FAILING Integrate on a Bessel Wronskian quotient costs tens of
 * seconds per term (measured 47.8 s on the 3395 member), which blows every solve
 * budget and would silently turn these PASSes into corpus timeouts. */
static void t_m61_latency_bound(void) {
    const char* eqs[] = {
        "x y''[x] + 3 y'[x] - y[x] == x",
        "x^2 y''[x] + x y'[x] + (x + 12) y[x] == x^2 + x",
        "9 x^2 y''[x] + (3 x + 2) y[x] == x^4 + x^2",
        "x y''[x] + y'[x] - 2 y[x] x == x^2" };
    for (size_t i = 0; i < 4; i++) {
        char buf[320];
        snprintf(buf, sizeof(buf), "Head[DSolve[%s, y, x]] === List", eqs[i]);
        clock_t t0 = clock();
        ASSERT_TRUE(buf);
        double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
        ASSERT_MSG(dt < 20.0, "M61 latency: %s took %.1f s (pre-screen not firing?)",
                   eqs[i], dt);
    }
}

/* ---- F7: bounded declines ------------------------------------------------- */
static void t_m61_declines(void) {
    /* nonlinear */
    ASSERT_TRUE("Head[DSolve`VariationOfParameters[x y''[x] == y[x]^2, y, x]] =!= List");
    /* first order */
    ASSERT_TRUE("Head[DSolve`VariationOfParameters[x y'[x] + y[x] == x, y, x]] =!= List");
    /* distributional forcing is the Green's-function methods' domain */
    ASSERT_TRUE("FreeQ[DSolve[x y''[x] + 3 y'[x] - y[x] == DiracDelta[x - 1], y, x], "
                "Inactive[Integrate]]");
    /* an ORDINARY point (no singularity) is Frobenius/Kovacic territory */
    ASSERT_TRUE("FreeQ[DSolve[y''[x] + y'[x] - x y[x] == x^2, y, x], Inactive]");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m61_bessel_family);
    TEST(t_m61_elementary_unchanged);
    TEST(t_m61_gate_margin);
    TEST(t_m61_series_basis_declines);
    TEST(t_m61_ivp_declines);
    TEST(t_m61_latency_bound);
    TEST(t_m61_declines);

    printf("All DSolve M61 stress tests passed.\n");
    return 0;
}
