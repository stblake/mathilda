/*
 * test_dsolve_m20_stress.c — anti-overfit stress families for M20
 * (DSolve`PolynomialShiftSubstitution: y' == -phi'/c + g(x) (phi(x) + c y)^p, the
 * x-dependent-shift generalisation of FirstOrderSubstitution, reducing via
 * u = phi(x) + c y to the separable u' == c g(x) u^p and returning the implicit
 * first integral u^(1-p)/(1-p) - Integrate[c g,x] == C[1]).
 *
 * Forward generator over a genuine (phi, c, g, p) grid: from chosen data it builds
 * the ODE whose closed form is guaranteed, so a fix to one example is proven not to
 * be an overfit by requiring the whole grid to solve.  The output is IMPLICIT, so
 * verification is the implicit-function rule: with the solution G(x,y) == C[1], the
 * residual D[G,x] + D[G,y[x]]*RHS must vanish (branch-safe — the radical cancels).
 * Head === List and an Equal relation are checked first so a decline can't pass.
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

/* Build  y' == rhs, solve, require an implicit Equal solution, then verify the
 * implicit-function residual  (dG/dx with y' -> rhs) == 0 numerically at `subs`.
 * D[g,x] is the TOTAL x-derivative (= g_x + g_y y'), so substituting the ODE
 * y' -> rhs into it yields g_x + g_y rhs, which vanishes along the solution
 * (branch-safe — the radical cancels). */
static void implicit_ok(const char* rhs, const char* subs) {
    char buf[2400];
    snprintf(buf, sizeof(buf),
        "With[{sol = DSolve[y'[x] == %s, y, x]}, "
        "Head[sol] === List && Length[sol] >= 1 && Head[sol[[1, 1]]] === Equal && "
        "With[{g = sol[[1, 1, 1]] - sol[[1, 1, 2]]}, "
        "Abs[N[(D[g, x] /. y'[x] -> (%s)) /. {%s}, 20]] < 10^-6]]",
        rhs, rhs, subs);
    ASSERT_TRUE(buf);
}

/* Family A — pure radical p = 1/2 over a (phi, c, g) grid.  R = -phi'/c. */
static void t_m20_sqrt_grid(void) {
    /* {phi, c, g, R = -phi'/c, x-val, y-val} */
    struct { const char* phi; const char* c; const char* g; const char* R;
             const char* xv; const char* yv; } T[] = {
        { "x^3",       "-6", "1 + x^2 + x^3", "x^2/2",       "6/5", "1/10" },
        { "x^4",       "8",  "x^3/(1 + x)",   "-x^3/2",      "6/5", "1/10" },
        { "x^2 + x",   "3",  "x^2",           "-(2 x + 1)/3","7/5", "1/20" },
        { "x^3 + 2 x", "-2", "1 + x",         "(3 x^2 + 2)/2","6/5","1/10" },
    };
    for (int i = 0; i < 4; i++) {
        char rhs[300], subs[120];
        snprintf(rhs, sizeof(rhs), "(%s) + (%s) ((%s) + (%s) y[x])^(1/2)",
                 T[i].R, T[i].g, T[i].phi, T[i].c);
        snprintf(subs, sizeof(subs), "x -> %s, y[x] -> %s", T[i].xv, T[i].yv);
        implicit_ok(rhs, subs);
    }
}

/* Family B — general fractional power p over a grid (cube / two-thirds roots). */
static void t_m20_power_grid(void) {
    struct { const char* phi; const char* c; const char* g; const char* R;
             const char* p; const char* xv; const char* yv; } T[] = {
        { "x^2", "-3", "1 + x",  "2 x/3",   "1/3", "6/5", "1/20" },
        { "x^3", "6",  "x",      "-x^2/2",  "2/3", "6/5", "1/10" },
        { "x",   "2",  "1 + x^2","-1/2",    "1/3", "7/5", "1/20" },
    };
    for (int i = 0; i < 3; i++) {
        char rhs[300], subs[120];
        snprintf(rhs, sizeof(rhs), "(%s) + (%s) ((%s) + (%s) y[x])^(%s)",
                 T[i].R, T[i].g, T[i].phi, T[i].c, T[i].p);
        snprintf(subs, sizeof(subs), "x -> %s, y[x] -> %s", T[i].xv, T[i].yv);
        implicit_ok(rhs, subs);
    }
}

/* Pinned: corpus flagships 402/371 solve; an equation whose radical base is
 * NONLINEAR in y (y^2 inside) is not the method's class -> not solved by it. */
static void t_m20_pinned(void) {
    /* 2.1.2-402 */
    ASSERT_TRUE("With[{sol = DSolve`PolynomialShiftSubstitution["
                "y'[x] == x^2/2 + (1 + x^2 + x^3) Sqrt[x^3 - 6 y[x]], y, x]}, "
                "Head[sol] === List && Length[sol] >= 1 && Head[sol[[1, 1]]] === Equal]");
    /* 2.1.2-371 (symbolic parameter a) */
    ASSERT_TRUE("Head[DSolve`PolynomialShiftSubstitution["
                "y'[x] == -(1/2) Sqrt[a] x^3 (Sqrt[a] + Sqrt[a] x - 2 Sqrt[a x^4 + 8 y[x]])/(1 + x), y, x]] === List");
    /* base Sqrt[x + y^2] is nonlinear in y -> declines (leaves unevaluated) */
    ASSERT_TRUE("Head[DSolve`PolynomialShiftSubstitution["
                "y'[x] == Sqrt[x + y[x]^2], y, x]] =!= List");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m20_sqrt_grid);
    TEST(t_m20_power_grid);
    TEST(t_m20_pinned);

    printf("All DSolve M20 stress tests passed.\n");
    return 0;
}
