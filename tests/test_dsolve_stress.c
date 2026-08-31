/*
 * test_dsolve_stress.c — parametrized/forward-generator stress tests spanning the
 * DSolve cascade methods that were previously exercised only by hand-written point
 * cases (Kovacic and Frobenius have their own file, test_dsolve_m5_stress.c).
 *
 * Each family is a FORWARD GENERATOR: the equation is built from parameters whose
 * closed-form solution is guaranteed (a chosen root set, a potential function, an
 * integrating-factor pair, ...), so DSolve must solve it and we verify by
 * back-substitution.  As with every DSolve stress test we assert the VALIDITY of
 * whatever is returned — Head===List first (a declined solve leaves [[1]] symbolic
 * and would let PossibleZeroQ pass vacuously), then the residual ~ 0 — never a
 * fixed solution form.  The whole run stays far under the 60 s harness alarm.
 *
 * Where a method registers a backtick builtin we PIN it (so the family really
 * exercises that method, not a neighbour that happens to catch the same form);
 * the systems and PDE methods have no backtick builtin and are driven through the
 * automatic DSolve dispatch.
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

/* ---- generic solve-and-verify helpers ---- */

/* pinned-method: Head[method[eqn,y,x]]===List, then residual back-substitutes. */
static void method_ok(const char* method, const char* eqn, const char* resid) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Head[%s[%s, y, x]] === List", method, eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[(%s) /. %s[%s, y, x][[1]]]", resid, method, eqn);
    ASSERT_TRUE(buf);
}

/* ---- per-family generators (each builds eqn+residual then calls method_ok) ---- */

/* Target-solution generator: choose the coefficient p and a target solution yt,
 * then q = yt' + p yt so that yt is a particular solution and the integrating-
 * factor integral int(mu q) = mu yt is elementary BY CONSTRUCTION — this keeps
 * every (p, yt) pair in-domain (a raw (p, q) grid can hit a non-elementary
 * integral, e.g. y' + x y == Sin[x] needs int Sin[x] Exp[x^2/2]). */
static void lin1_ok(const char* p, const char* yt) {
    char eqn[320], res[320];
    snprintf(eqn, sizeof(eqn),
             "y'[x] + (%s) y[x] == (D[%s, x] + (%s) (%s))", p, yt, p, yt);
    snprintf(res, sizeof(res),
             "y'[x] + (%s) y[x] - (D[%s, x] + (%s) (%s))", p, yt, p, yt);
    method_ok("DSolve`LinearFirstOrder", eqn, res);
}

static void sep_ok(const char* g, const char* h) {
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y'[x] == (%s) (%s)", g, h);
    snprintf(res, sizeof(res), "y'[x] - (%s) (%s)", g, h);
    method_ok("DSolve`Separable", eqn, res);
}

static void bern_ok(const char* p, const char* q, const char* n) {
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y'[x] + (%s) y[x] == (%s) y[x]^(%s)", p, q, n);
    snprintf(res, sizeof(res), "y'[x] + (%s) y[x] - (%s) y[x]^(%s)", p, q, n);
    method_ok("DSolve`Bernoulli", eqn, res);
}

static void hom_ok(const char* rhs) {
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y'[x] == (%s)", rhs);
    snprintf(res, sizeof(res), "y'[x] - (%s)", rhs);
    method_ok("DSolve`Homogeneous", eqn, res);
}

/* transcendental homogeneous -> implicit first integral G(x,y[x]) == C[1] (an
 * Equal branch, so it is verified by implicit differentiation, not /. rule). */
static void hom_implicit_ok(const char* rhs) {
    char buf[768];
    snprintf(buf, sizeof(buf), "Head[DSolve[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{eq = DSolve[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]] - eq[[2]], x] /. y'[x] -> (%s)]]", rhs, rhs);
    ASSERT_TRUE(buf);
}

/* Exact forward generator: pick a potential phi(x,y); M = phi_x, N = phi_y are an
 * exact pair by construction; the equation is M + N y' == 0 with y -> y[x]. */
static void exact_ok(const char* phi) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "With[{M = D[%s, x], N = D[%s, y]}, "
        "Head[DSolve`Exact[(M /. y -> y[x]) + (N /. y -> y[x]) y'[x] == 0, y, x]] === List]",
        phi, phi);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "With[{M = D[%s, x], N = D[%s, y]}, "
        "PossibleZeroQ[((M /. y -> y[x]) + (N /. y -> y[x]) y'[x]) /. "
        "DSolve`Exact[(M /. y -> y[x]) + (N /. y -> y[x]) y'[x] == 0, y, x][[1]]]]",
        phi, phi);
    ASSERT_TRUE(buf);
}

/* constant-coefficient generators (from a chosen spectrum). */
static void cc_ok(int c1, int c0) {  /* y'' - c1 y' + c0 y == 0 */
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y''[x] - (%d) y'[x] + (%d) y[x] == 0", c1, c0);
    snprintf(res, sizeof(res), "y''[x] - (%d) y'[x] + (%d) y[x]", c1, c0);
    method_ok("DSolve`LinearConstantCoefficients", eqn, res);
}
static void cc_inh_ok(const char* forcing) {  /* roots 1,2 + forcing */
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y''[x] - 3 y'[x] + 2 y[x] == (%s)", forcing);
    snprintf(res, sizeof(res), "y''[x] - 3 y'[x] + 2 y[x] - (%s)", forcing);
    method_ok("DSolve`LinearConstantCoefficients", eqn, res);
}

/* Euler-Cauchy from indicial roots r1,r2: x^2 y'' + (1-r1-r2) x y' + r1 r2 y == 0. */
static void euler_ok(int r1, int r2) {
    int c1 = 1 - r1 - r2, c0 = r1 * r2;
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "x^2 y''[x] + (%d) x y'[x] + (%d) y[x] == 0", c1, c0);
    snprintf(res, sizeof(res), "x^2 y''[x] + (%d) x y'[x] + (%d) y[x]", c1, c0);
    method_ok("DSolve`EulerCauchy", eqn, res);
}

/* reduction of order, 2nd-order missing y: y'' == f(x) y'. */
static void ro_ok(const char* f) {
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y''[x] == (%s) y'[x]", f);
    snprintf(res, sizeof(res), "y''[x] - (%s) y'[x]", f);
    method_ok("DSolve`ReductionOfOrder", eqn, res);
}

/* reduction of order, general missing-y RHS y'' == F(x, y'). */
static void ro_full_ok(const char* F) {
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y''[x] == %s", F);
    snprintf(res, sizeof(res), "y''[x] - (%s)", F);
    method_ok("DSolve`ReductionOfOrder", eqn, res);
}

/* Riccati from a chosen spectrum r1,r2 (q2==1): the linearisation
 * u'' - (r1+r2) u' + r1 r2 u == 0 has roots r1,r2, so the Riccati equation
 * y' == r1 r2 + (r1+r2) y + y^2 has an ELEMENTARY (exp-ratio) solution — the
 * back-substituted residual is decidable. */
static void riccati_ok(int r1, int r2) {
    int q1 = r1 + r2, q0 = r1 * r2;
    char eqn[256], res[256];
    snprintf(eqn, sizeof(eqn), "y'[x] == (%d) + (%d) y[x] + y[x]^2", q0, q1);
    snprintf(res, sizeof(res), "y'[x] - ((%d) + (%d) y[x] + y[x]^2)", q0, q1);
    method_ok("DSolve`Riccati", eqn, res);
}

/* 2x2 constant-coefficient system y'=A y (auto dispatch — no backtick builtin). */
static void sys2_ok(int a, int b, int c, int d) {
    char sys[256], buf[1024];
    snprintf(sys, sizeof(sys),
        "{y'[x] == (%d) y[x] + (%d) z[x], z'[x] == (%d) y[x] + (%d) z[x]}", a, b, c, d);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, {y, z}, x]] === List", sys);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "And @@ (PossibleZeroQ /@ ({y'[x] - ((%d) y[x] + (%d) z[x]), "
        "z'[x] - ((%d) y[x] + (%d) z[x])} /. DSolve[%s, {y, z}, x][[1]]))",
        a, b, c, d, sys);
    ASSERT_TRUE(buf);
}

/* triangular (possibly variable-coefficient) 2-function system. */
static void tri_ok(const char* ry, const char* rz) {
    char sys[256], buf[1024];
    snprintf(sys, sizeof(sys), "{y'[x] == %s, z'[x] == %s}", ry, rz);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, {y, z}, x]] === List", sys);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "And @@ (PossibleZeroQ /@ ({y'[x] - (%s), z'[x] - (%s)} /. "
        "DSolve[%s, {y, z}, x][[1]]))", ry, rz, sys);
    ASSERT_TRUE(buf);
}

/* first-order linear PDE a u_x + b u_y == f (auto dispatch); verify against a
 * concrete arbitrary function C[1][z_] :> Sin[z] (see the ODE-file PDE tests). */
static void pde_ok(int a, int b, const char* f) {
    char eqn[256], buf[1024];
    snprintf(eqn, sizeof(eqn),
        "(%d) D[u[x,y],x] + (%d) D[u[x,y],y] == (%s)", a, b, f);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, u, {x, y}]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "With[{uc = (u[x,y] /. DSolve[%s, u, {x, y}][[1]]) /. C[1][z_] :> Sin[z]}, "
        "PossibleZeroQ[(%d) D[uc,x] + (%d) D[uc,y] - (%s)]]", eqn, a, b, f);
    ASSERT_TRUE(buf);
}

/* ---------------------- families ---------------------- */

static void t_stress_linear1(void) {
    const char* ps[] = {"1", "x", "1/x", "2"};
    const char* yts[] = {"x", "x^2", "Sin[x]", "Exp[x]"};
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 4; j++)
            lin1_ok(ps[i], yts[j]);
}

static void t_stress_separable(void) {
    const char* gs[] = {"1", "x"};
    const char* hs[] = {"y[x]", "y[x]^2", "1 + y[x]^2"};
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 3; j++)
            sep_ok(gs[i], hs[j]);
}

static void t_stress_bernoulli(void) {
    /* (p, q, n) triples spanning n = 2, 3 with constant/variable p, q. */
    bern_ok("1", "1", "2");
    bern_ok("1/x", "1", "2");
    bern_ok("-1", "x", "3");
    bern_ok("1", "1", "3");
}

static void t_stress_homogeneous(void) {
    /* Homogeneous inverts the pure-log (real-root) rational family via the
     * exponentiate-and-clear-radicals fallback (homog_exp_log_invert); the
     * transcendental (ArcTan) subset has no explicit inverse and still declines,
     * so this is a curated set of confirmed in-domain forms, not a dense grid. */
    const char* fs[] = {
        "(x - y[x])/(x + y[x])",
        "(2 x - y[x])/(x + y[x])",
        "y[x]/x + (y[x]/x)^2",
        "y[x]/x + x/y[x]",
        "(x + 2 y[x])/(2 x + y[x])",   /* algebraic fallback: cubic Root branches */
        "(2 x + y[x])/(x + 2 y[x])",
        "(x + 3 y[x])/(3 x + y[x])",
        "(2 x + 3 y[x])/(3 x + 2 y[x])"
    };
    for (size_t i = 0; i < 8; i++) hom_ok(fs[i]);
}

static void t_stress_homogeneous_implicit(void) {
    /* transcendental (ArcTan) rational-homogeneous family: no explicit inverse,
     * returned as the implicit first integral and verified by implicit diff. */
    const char* fs[] = {
        "(x + y[x])/(x - y[x])",
        "(x + y[x])/(2 x + y[x])",
        "(3 x + y[x])/(x + 2 y[x])",
        "(x + y[x])/(x - 2 y[x])",
        "x/(2 x + y[x])"
    };
    for (size_t i = 0; i < 5; i++) hom_implicit_ok(fs[i]);
}

static void t_stress_exact(void) {
    const char* phis[] = {
        "x^2 y", "x^2 + x y + y^2", "x y^2", "x^3 + x y", "x^2 - y^2", "x^3 y"
    };
    for (size_t i = 0; i < 6; i++) exact_ok(phis[i]);
}

static void t_stress_constcoeff(void) {
    /* real distinct roots (r1,r2) -> c1=r1+r2, c0=r1 r2 */
    int rr[][2] = {{1, 2}, {-1, 3}, {0, 1}, {-2, -3}};
    for (size_t i = 0; i < 4; i++)
        cc_ok(rr[i][0] + rr[i][1], rr[i][0] * rr[i][1]);
    /* complex pairs a +- b i -> c1=2a, c0=a^2+b^2 */
    int cx[][2] = {{0, 2}, {1, 1}, {-1, 3}};
    for (size_t i = 0; i < 3; i++)
        cc_ok(2 * cx[i][0], cx[i][0] * cx[i][0] + cx[i][1] * cx[i][1]);
    /* repeated root r -> c1=2r, c0=r^2 */
    int rep[] = {1, 2, -1};
    for (size_t i = 0; i < 3; i++) cc_ok(2 * rep[i], rep[i] * rep[i]);
    /* inhomogeneous (variation of parameters) */
    cc_inh_ok("1");
    cc_inh_ok("x");
    cc_inh_ok("Exp[x]");
}

static void t_stress_euler(void) {
    int rr[][2] = {{1, 2}, {-1, -1}, {2, 3}, {1, -1}};
    for (size_t i = 0; i < 4; i++) euler_ok(rr[i][0], rr[i][1]);
}

static void t_stress_reduce_order(void) {
    const char* fs[] = {"1/x", "2", "-1", "3", "2/x"};
    for (size_t i = 0; i < 5; i++) ro_ok(fs[i]);
}

static void t_stress_reduce_order_riccati(void) {
    /* autonomous a + b (y')^2 -> Tan/Tanh, and Riccati-in-p c x (y')^2 -> 1/(C+..):
     * the first-order sub-solve is exact but the antiderivative is left in a form
     * only PossibleZeroQ can confirm (guard fallback in dsolve_reduce_order.c). */
    const char* fs[] = {
        "1 + y'[x]^2", "1 - y'[x]^2", "2 + y'[x]^2", "1 + 2 y'[x]^2",
        "-2 x y'[x]^2", "-x y'[x]^2", "y'[x]^2/x"
    };
    for (size_t i = 0; i < 7; i++) ro_full_ok(fs[i]);
}

static void t_stress_riccati(void) {
    /* spectra r1,r2: distinct signs/magnitudes + a repeated root, all giving an
     * elementary exp-ratio solution the residual PossibleZeroQ can decide. */
    int rr[][2] = {
        {1, 2}, {-1, -2}, {2, 3}, {1, -1}, {-1, 3}, {1, -3}, {2, 2}
    };
    for (size_t i = 0; i < 7; i++) riccati_ok(rr[i][0], rr[i][1]);
}

static void t_stress_systems(void) {
    /* (a,b,c,d): distinct-real, complex, defective, defective+singular. */
    int m[][4] = {
        {2, 1, 1, 2},   /* eigenvalues 1, 3 */
        {1, -2, 1, -1}, /* eigenvalues +- i */
        {1, -1, 1, 3},  /* eigenvalue 2 doubled, defective */
        {0, 0, -1, 0}   /* defective AND singular */
    };
    for (size_t i = 0; i < 4; i++) sys2_ok(m[i][0], m[i][1], m[i][2], m[i][3]);
}

static void t_stress_triangular(void) {
    tri_ok("2 y[x]", "y[x] + z[x]");        /* constant, triangular */
    tri_ok("y[x]/x", "z[x]/x + y[x]");      /* variable coefficient */
    tri_ok("y[x]/x", "y[x]");               /* variable, z decoupled from own value */
}

static void t_stress_pde(void) {
    pde_ok(1, 1, "0");
    pde_ok(2, 3, "0");
    pde_ok(1, 2, "y");
    pde_ok(3, 5, "x");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(t_stress_linear1);
    TEST(t_stress_separable);
    TEST(t_stress_bernoulli);
    TEST(t_stress_homogeneous);
    TEST(t_stress_homogeneous_implicit);
    TEST(t_stress_exact);
    TEST(t_stress_constcoeff);
    TEST(t_stress_euler);
    TEST(t_stress_reduce_order);
    TEST(t_stress_reduce_order_riccati);
    TEST(t_stress_riccati);
    TEST(t_stress_systems);
    TEST(t_stress_triangular);
    TEST(t_stress_pde);

    printf("\nAll DSolve stress tests passed.\n");
    return 0;
}
