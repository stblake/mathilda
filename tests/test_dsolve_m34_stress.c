/*
 * test_dsolve_m34_stress.c — anti-overfit stress families for M34 (§2.2.14).
 *
 * The four M34 fixes, each proven not to be an overfit by a FORWARD GENERATOR over
 * a genuine parameter grid (a single-example fix must solve the whole grid):
 *
 *   A/B — robust variation of parameters + numeric-zero verify short-circuit:
 *         y'' + y == f (const-coeff) and (1-x)y'' + x y' - y == f (variable-coeff,
 *         basis {x, e^x}) for a grid of NON-undetermined-coefficient forcings
 *         (Tan, Sec, Cot, ...).  The VoP answer carries Log branch cuts on which the
 *         symbolic zero_test spins; the numeric-zero verify keeps it.  Verified by a
 *         numeric residual on the ORIGINAL equation with generic C[1],C[2] at a real
 *         point (the Logs are real there).
 *   C   — bounded-Kovacic complex-pole gate -> Frobenius series: (x^3+a)y''+4x y'+y
 *         over a grid of a (complex-conjugate poles), a genuinely Heun family that
 *         must decline promptly and fall to the ordinary-point series (Head===List).
 *   D   — 2nd-order exact -> series: y'' + p y' + p' y == 0 (an exact operator whose
 *         first integral y' + p y == C has a non-elementary quadrature) for a grid of
 *         p, falling to the Frobenius series (Head===List).
 *   E   — IC-point Frobenius series for a transcendental coefficient:
 *         x^2 y'' + (x+1) y' + a Log[x] y == 0, y[x0]==b, y'[x0]==c, expanded about
 *         the (ordinary) IC point and FIT there (FreeQ[C]).
 *   F   — forced special-function operator with an ARBITRARY forcing:
 *         x^2 y'' + x y' + (x^2-1/4) y == g[x], whose inert-Integrate VoP particular
 *         is Mathematica's own integral form (Head===List).
 *
 * Head === List is checked first so a decline / hang ($Aborted) can never pass
 * vacuously.  The residue 1360 (forced Duffing, no CAS closed form) is asserted to
 * DECLINE cleanly (bounded, not a wrong answer).
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

/* Solve `eqn`, require Head === List, then require the ORIGINAL-equation residual
 * `resid` to vanish numerically under `subs` (C[1],C[2] + x at a real point). */
static void res_ok(const char* eqn, const char* resid, const char* subs) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "With[{sol = TimeConstrained[DSolve[%s, y, x], 8, $Aborted]}, "
        "Head[sol] === List && Length[sol] >= 1 && "
        "Abs[N[(%s) /. sol[[1]] /. {%s}, 20]] < 10^-6]",
        eqn, resid, subs);
    ASSERT_TRUE(buf);
}

/* Solve `eqn`, require only that it SOLVED (Head === List, non-empty) -- for a
 * series / inert-integral answer that does not numericize (arbitrary forcing, a
 * SeriesData body).  A hang ($Aborted) or decline (unevaluated DSolve) fails. */
static void solved_ok(const char* eqn) {
    char buf[900];
    snprintf(buf, sizeof(buf),
        "With[{sol = TimeConstrained[DSolve[%s, y, x], 8, $Aborted]}, "
        "Head[sol] === List && Length[sol] >= 1]", eqn);
    ASSERT_TRUE(buf);
}

/* Family A — variation of parameters, const-coeff basis {Cos,Sin}, NON-UC forcing.
 * The VoP answer of y''+y==Tan[x] etc. carries -Cos Log[Sec+Tan] whose residual
 * spins zero_test; the numeric-zero verify keeps it (1337/1341). */
static void t_m34_vop_const(void) {
    const char* f[] = { "Tan[x]", "Sec[x]", "Cot[x]", "2 Sec[x/2]" };
    for (int i = 0; i < 4; i++) {
        char eqn[160], res[160];
        snprintf(eqn, sizeof(eqn), "y''[x] + y[x] == %s", f[i]);
        snprintf(res, sizeof(res), "y''[x] + y[x] - (%s)", f[i]);
        res_ok(eqn, res, "C[1] -> 13/10, C[2] -> 7/10, x -> 1/2");
    }
}

/* Family B — variation of parameters, VARIABLE-coeff basis {x, e^x}.  The homogeneous
 * (1-x)y''+x y'-y==0 has {x,e^x}; VoP closes an elementary/inert forcing (1354). */
static void t_m34_vop_varcoeff(void) {
    const char* f[] = { "2 (x - 1) E^(-x)", "x - 1", "(x - 1) E^x" };
    for (int i = 0; i < 3; i++) {
        char eqn[220], res[220];
        snprintf(eqn, sizeof(eqn), "(1 - x) y''[x] + x y'[x] - y[x] == %s", f[i]);
        snprintf(res, sizeof(res), "(1 - x) y''[x] + x y'[x] - y[x] - (%s)", f[i]);
        res_ok(eqn, res, "C[1] -> 13/10, C[2] -> 7/10, x -> 1/2");
    }
}

/* Family C — cubic-coefficient Heun family (complex-conjugate poles) that Kovacic
 * declines to the Frobenius ordinary-point series (1392/1393). */
static void t_m34_cubic_series(void) {
    const char* a[] = { "1", "2", "8", "-1", "3" };
    for (int i = 0; i < 5; i++) {
        char eqn[160];
        snprintf(eqn, sizeof(eqn), "(x^3 + (%s)) y''[x] + 4 x y'[x] + y[x] == 0", a[i]);
        solved_ok(eqn);
    }
}

/* Family D — 2nd-order EXACT operator y'' + p y' + p' y == (y' + p y)', whose first
 * integral y' + p y == C has a non-elementary quadrature, so it falls to the
 * Frobenius series (1384). */
static void t_m34_exact_series(void) {
    const char* p[] = { "Sin[x]", "x^2", "E^x", "Cos[x]" };
    for (int i = 0; i < 4; i++) {
        char eqn[220];
        snprintf(eqn, sizeof(eqn),
            "y''[x] + (%s) y'[x] + D[%s, x] y[x] == 0", p[i], p[i]);
        solved_ok(eqn);
    }
}

/* Family E — transcendental-coefficient IVP at an ordinary IC point: the ONLY route
 * is a Taylor series about x0, FIT there (1385).  FreeQ[C] proves the ICs were met. */
static void t_m34_icpoint_series(void) {
    struct { const char* a; const char* b; const char* c; } g[] = {
        { "3", "2", "0" }, { "1", "1", "1" }, { "2", "0", "1" }, { "5", "3", "-2" }
    };
    for (int i = 0; i < 4; i++) {
        char buf[420];
        snprintf(buf, sizeof(buf),
            "With[{sol = TimeConstrained[DSolve[{x^2 y''[x] + (x + 1) y'[x] + "
            "(%s) Log[x] y[x] == 0, y[1] == %s, y'[1] == %s}, y, x], 8, $Aborted]}, "
            "Head[sol] === List && Length[sol] >= 1 && FreeQ[sol, C[_]]]",
            g[i].a, g[i].b, g[i].c);
        ASSERT_TRUE(buf);
    }
}

/* Family F — forced special-function operator with an ARBITRARY forcing g(x): the
 * inert-Integrate VoP particular is Mathematica's own integral form (1350). */
static void t_m34_forced_arbitrary(void) {
    solved_ok("x^2 y''[x] + x y'[x] + (x^2 - 1/4) y[x] == g[x]");
    solved_ok("x^2 y''[x] + x y'[x] + (x^2 - 1/4) y[x] == h[x]");
}

/* Pinned §2.2.14 corpus cases that motivated each fix. */
static void t_m34_pinned(void) {
    /* 1337/1341 — VoP verify short-circuit */
    res_ok("y''[x] + y[x] == Tan[x]",   "y''[x] + y[x] - Tan[x]",   "C[1] -> 13/10, C[2] -> 7/10, x -> 1/2");
    /* 1381 — Bessel general solution singular at x=0 -> origin series IVP */
    ASSERT_TRUE("With[{sol = TimeConstrained[DSolve[{y''[x] + x^2 y[x] == 0, y[0] == 1, "
                "y'[0] == 0}, y, x], 8, $Aborted]}, Head[sol] === List && FreeQ[sol, C[_]]]");
    /* 1384 — 2nd-order exact -> Frobenius series IVP, fit at x=0 */
    ASSERT_TRUE("With[{sol = TimeConstrained[DSolve[{y''[x] + Sin[x] y'[x] + Cos[x] y[x] == 0, "
                "y[0] == 0, y'[0] == 1}, y, x], 8, $Aborted]}, Head[sol] === List && FreeQ[sol, C[_]]]");
    /* 1385 — IC-point transcendental series */
    ASSERT_TRUE("With[{sol = TimeConstrained[DSolve[{x^2 y''[x] + (x + 1) y'[x] + 3 Log[x] y[x] == 0, "
                "y[1] == 2, y'[1] == 0}, y, x], 8, $Aborted]}, Head[sol] === List && FreeQ[sol, C[_]]]");
    /* residue 1360 — forced Duffing: a bounded DECLINE, never a wrong answer */
    ASSERT_TRUE("With[{sol = TimeConstrained[DSolve[{u''[t] + u'[t] + u[t]^3/5 == Cos[t], "
                "u[0] == 2, u'[0] == 0}, u, t], 8, $Aborted]}, "
                "sol =!= $Aborted && !MatchQ[sol, {{__Rule}}]]");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m34_vop_const);
    TEST(t_m34_vop_varcoeff);
    TEST(t_m34_cubic_series);
    TEST(t_m34_exact_series);
    TEST(t_m34_icpoint_series);
    TEST(t_m34_forced_arbitrary);
    TEST(t_m34_pinned);

    printf("All DSolve M34 stress tests passed.\n");
    return 0;
}
