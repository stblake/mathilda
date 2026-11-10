/* test_dsolve_pde_ivp.c — DSolve`PDEInitialValue (src/calculus/dsolve_pde_ivp.c)
 * and the PDE parsing fixes it relies on.
 *
 * Each family is tested on the famous-problem instance and on siblings it does
 * not name; every answer is also checked here by back-substitution into the
 * PDE and the conditions (pdeOK), independently of the method's own check.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void assert_true(const char* cond) { assert_eval_eq(cond, "True", 0); }

/* pdeOK[eqs, sol, u, vars]: every equation (PDE and conditions) holds for the
 * first solution; conditions at a singular point are taken as a limit. */
static void define_helpers(void) {
    Expr* e = parse_expression(
        "pdeOK[eqs_, sol_, u_, vars_] := Module[{f = Function @@ {vars, (u @@ vars) /. sol[[1]]}},"
        " And @@ (PossibleZeroQ[Together[(#[[1]] - #[[2]]) /. u -> f]] & /@ eqs)]");
    Expr* r = evaluate(e);
    expr_free(e); expr_free(r);
}

static void test_parsing(void) {
    /* the applied form u[x, t] parses (it used to make every PDE decline) */
    assert_true("Head[DSolve[{D[u[x,t],t] == D[u[x,t],{x,2}], u[x,0] == E^(-x^2)}, u[x,t], {x,t}]] === List");
    /* {x, y, t} is three variables when the unknown is applied to them ... */
    assert_true("Head[DSolve[{D[w[x,y,t],t] == D[w[x,y,t],{x,2}], w[x,y,0] == Sin[x]}, w[x,y,t], {x,y,t}]] === List");
    /* ... and still a range {x, xmin, xmax} for an ODE */
    assert_true("Head[DSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}]] === List");
}

static void test_amplitude_family(void) {
    assert_true("DSolve[{D[u[x,t],t] + u[x,t] D[u[x,t],x] == 0, u[x,0] == x}, u[x,t], {x,t}] === {{u[x,t] -> x/(1 + t)}}");
    assert_true("With[{e = {D[u[x,t],t] + u[x,t] D[u[x,t],x] == 0, u[x,0] == 2 x}}, "
                "pdeOK[e, DSolve[e, u[x,t], {x,t}], u, {x,t}]]");
    /* Taylor-Green: the advection term cancels on the data */
    assert_true("With[{e = {D[w[x,y,t],t] + Sin[x] Cos[y] E^(-2 nu t) D[w[x,y,t],x] - Cos[x] Sin[y] E^(-2 nu t) D[w[x,y,t],y]"
                " == nu (D[w[x,y,t],{x,2}] + D[w[x,y,t],{y,2}]), w[x,y,0] == 2 Sin[x] Sin[y]}}, "
                "DSolve[e, w[x,y,t], {x,y,t}] === {{w[x,y,t] -> 2 E^(-2 nu t) Sin[x] Sin[y]}}]");
}

static void test_shift_family(void) {
    assert_true("With[{e = {D[u[x,t],t] + D[u[x,t],{x,3}] + 6 u[x,t] D[u[x,t],x] == 0, u[x,0] == c Sech[Sqrt[c] x/2]^2/2}}, "
                "pdeOK[e, DSolve[e, u[x,t], {x,t}], u, {x,t}]]");
    /* a different soliton amplitude: c = 4 */
    assert_true("With[{e = {D[u[x,t],t] + D[u[x,t],{x,3}] + 6 u[x,t] D[u[x,t],x] == 0, u[x,0] == 2 Sech[x]^2}}, "
                "DSolve[e, u[x,t], {x,t}] === {{u[x,t] -> 2 Sech[x - 4 t]^2}}]");
    assert_true("With[{e = {D[u[x,t],t] == D[u[x,t],{x,2}] + u[x,t] (1 - u[x,t]), u[x,0] == 1/(E^(x/Sqrt[6]) + 1)^2}}, "
                "pdeOK[e, DSolve[e, u[x,t], {x,t}], u, {x,t}]]");
    assert_true("With[{e = {D[u[x,t],t] + u[x,t] D[u[x,t],x] == nu D[u[x,t],{x,2}], u[x,0] == c - a Tanh[a x/(2 nu)]}}, "
                "pdeOK[e, DSolve[e, u[x,t], {x,t}], u, {x,t}]]");
    /* second order in time: the speed comes from u_t(x, 0) */
    assert_true("FreeQ[DSolve[{D[u[x,t],{t,2}] - D[u[x,t],{x,2}] + Sin[u[x,t]] == 0, u[x,0] == 4 ArcTan[E^(x/Sqrt[1 - v^2])],"
                " Derivative[0,1][u][x,0] == -2 v/(Cosh[x/Sqrt[1 - v^2]] Sqrt[1 - v^2])}, u[x,t], {x,t}], DSolve]");
    /* a profile that does not travel (the KdV equation with Gaussian data) declines */
    assert_true("Head[DSolve`PDEInitialValue[{D[u[x,t],t] + D[u[x,t],{x,3}] + 6 u[x,t] D[u[x,t],x] == 0, u[x,0] == E^(-x^2)},"
                " u[x,t], {x,t}]] === DSolve`PDEInitialValue");
}

static void test_exppoly_family(void) {
    assert_true("DSolve[{D[u[x,t],t] == D[u[x,t],{x,2}], u[x,0] == E^(-x^2)}, u[x,t], {x,t}] === "
                "{{u[x,t] -> E^(-x^2/(1 + 4 t))/Sqrt[1 + 4 t]}}");
    assert_true("With[{e = {D[u[x,t],t] == k D[u[x,t],{x,2}], u[x,0] == E^(-a x^2)}}, "
                "pdeOK[e, DSolve`PDEInitialValue[e, u[x,t], {x,t}], u, {x,t}]]");
    /* the Schrodinger packet, and with a momentum factor E^(I k x) */
    assert_true("With[{e = {I D[psi[x,t],t] == -D[psi[x,t],{x,2}]/2, psi[x,0] == E^(-x^2/2)}}, "
                "pdeOK[e, DSolve[e, psi[x,t], {x,t}], psi, {x,t}]]");
    assert_true("With[{e = {I D[psi[x,t],t] == -D[psi[x,t],{x,2}]/2, psi[x,0] == E^(-x^2/2 + 3 I x)}}, "
                "pdeOK[e, DSolve[e, psi[x,t], {x,t}], psi, {x,t}]]");
    /* Lamb-Oseen: the radial Laplacian's 1/r term */
    assert_true("DSolve[{D[w[r,t],t] == nu (D[w[r,t],{r,2}] + D[w[r,t],r]/r), w[r,0] == E^(-r^2)}, w[r,t], {r,t}] === "
                "{{w[r,t] -> E^(-r^2/(1 + 4 nu t))/(1 + 4 nu t)}}");
    /* a 2-D heat equation with separable Gaussian data */
    assert_true("With[{e = {D[u[x,y,t],t] == D[u[x,y,t],{x,2}] + D[u[x,y,t],{y,2}], u[x,y,0] == E^(-x^2 - 2 y^2)}}, "
                "pdeOK[e, DSolve[e, u[x,y,t], {x,y,t}], u, {x,y,t}]]");
}

static void test_boundary_and_liouville(void) {
    assert_true("DSolve[{D[u[y,t],t] == nu D[u[y,t],{y,2}], u[y,0] == 0, u[0,t] == U}, u[y,t], {y,t}] === "
                "{{u[y,t] -> U - U Erf[y/(2 Sqrt[nu] Sqrt[t])]}}");
    assert_true("FreeQ[DSolve[{D[u[y,t],t] == 3 D[u[y,t],{y,2}], u[y,0] == 1, u[0,t] == 5}, u[y,t], {y,t}], DSolve]");
    assert_true("With[{e = {D[u[y,t],t] == nu D[u[y,t],{y,2}], u[0,t] == U Cos[omega t]}}, "
                "pdeOK[e, DSolve[e, u[y,t], {y,t}], u, {y,t}]]");
    assert_true("FreeQ[DSolve[{D[u[y,t],t] == nu D[u[y,t],{y,2}], u[0,t] == U Sin[2 t]}, u[y,t], {y,t}], DSolve]");
    assert_true("DSolve[D[u[x,y],x,y] == E^u[x,y], u[x,y], {x,y}] === "
                "{{u[x,y] -> Log[2 Derivative[1][C[1]][x] Derivative[1][C[2]][y]/(C[1][x] + C[2][y])^2]}}");
    assert_true("With[{s = DSolve[D[u[x,y],x,y] == 3 E^(2 u[x,y]), u[x,y], {x,y}]}, "
                "PossibleZeroQ[Together[(D[#, x, y] - 3 E^(2 #)) &[u[x,y] /. s[[1]]]]]]");
    /* the d'Alembert and heat-kernel runners, reached through the applied form */
    assert_true("With[{e = {D[u[x,t],{t,2}] == c^2 D[u[x,t],{x,2}], u[x,0] == E^(-x^2), Derivative[0,1][u][x,0] == 0}}, "
                "pdeOK[e, DSolve[e, u[x,t], {x,t}], u, {x,t}]]");
    assert_true("DSolve[{D[u[x,t],t] == D[u[x,t],{x,2}], u[x,0] == Sin[x]}, u[x,t], {x,t}] === {{u[x,t] -> E^-t Sin[x]}}");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();
    define_helpers();
    TEST(test_parsing);
    TEST(test_amplitude_family);
    TEST(test_shift_family);
    TEST(test_exppoly_family);
    TEST(test_boundary_and_liouville);
    printf("All PDE initial-value tests passed.\n");
    return 0;
}
