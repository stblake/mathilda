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

/* Higher-order exact linear, built as the total derivative L[y] = d/dx(b1 y'+b0 y):
 * exact by construction.  (b1,b0) are chosen so the reduced first-order solve
 * b1 y' + b0 y == C[2] is elementary; pin DSolve`ExactODE (in AUTO some of these
 * are also Euler-Cauchy and would be claimed there). */
static void exact_ode_ok(const char* b1, const char* b0) {
    char eqn[512], res[512];
    snprintf(eqn, sizeof(eqn), "D[(%s) y'[x] + (%s) y[x], x] == 0", b1, b0);
    snprintf(res, sizeof(res), "D[(%s) y'[x] + (%s) y[x], x]", b1, b0);
    method_ok("DSolve`ExactODE", eqn, res);
}

/* OperatorFactor: compose first-order factors (D - r_i) into a monic operator L,
 * pin DSolve`OperatorFactor, verify Head===List then back-substitution.  Factors
 * are given as a Mathilda list string; choose them (shifted-Euler / all-constant)
 * so every reduction integral is elementary by construction. */
static void operfactor_ok(const char* rs) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Module[{rs=%s, w=y[x], op, s}, Do[w=D[w,x]-rs[[k]] w,{k,Length[rs]}]; op=Expand[w]; "
        "s=DSolve`OperatorFactor[op==0, y, x]; "
        "Head[s]===List && PossibleZeroQ[op /. s[[1]]]]", rs);
    ASSERT_TRUE(buf);
}

/* DFactor: factor the composed operator and reconstruct it from the returned
 * factors (apply the Dx-factors, rightmost first, to a concrete test function). */
static void dfactor_ok(const char* rs) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Module[{rs=%s, w, op, optf, fs, recon, tf=Exp[x]+x^2+x^5}, "
        "w=y[x]; Do[w=D[w,x]-rs[[k]] w,{k,Length[rs]}]; op=Expand[w]; "
        "w=tf; Do[w=D[w,x]-rs[[k]] w,{k,Length[rs]}]; optf=w; "     /* reference: L[tf] */
        "fs = DSolve`DFactor[op==0, y[x], x]; "
        "recon = Fold[Function[{u,f}, D[u,x] + (f /. Dx->0) u], tf, fs]; "
        "Head[fs]===List && Length[fs]==Length[rs] && PossibleZeroQ[recon - optf]]", rs);
    ASSERT_TRUE(buf);
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

/* Lagrange/d'Alembert forward generator: y == x phi(y') + psi(y') built from
 * phi(p), psi(p) (in the marker symbol p). Auto dispatch → parametric output
 * {x->Function[{t},X], y->Function[{t},Y]}; verify by substituting x->X(t),
 * y->Y(t), y'->Y'(t)/X'(t) into the original equation residual. phi=2p (and
 * 3p/2) with polynomial psi keep the integrating-factor integral elementary. */
static void lagrange_ok(const char* phi, const char* psi) {
    char eqn[512], buf[1500];
    snprintf(eqn, sizeof(eqn),
        "y[x] == x (%s /. p -> y'[x]) + (%s /. p -> y'[x])", phi, psi);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "Module[{s = DSolve[%s, y, x][[1]], X, Y, yp}, "
        "X = (x /. s)[t]; Y = (y /. s)[t]; yp = D[Y,t]/D[X,t]; "
        "PossibleZeroQ[Y - (X (%s /. p -> yp) + (%s /. p -> yp))]]",
        eqn, phi, psi);
    ASSERT_TRUE(buf);
}

/* implicit first-integral verify (auto dispatch): DSolve returns {{G==C[1]}} and
 * D[G,x] with y' -> rhs vanishes (Chini/Abel, homogeneous log-spirals). */
static void impl_ok(const char* rhs) {
    char buf[1500];
    snprintf(buf, sizeof(buf), "Head[DSolve[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{eq = DSolve[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> (%s)]]", rhs, rhs);
    ASSERT_TRUE(buf);
}

/* Lie `bivariate` forward generator: every omega = y/x + A(y/x)/x has the
 * genuinely-degree-2 point symmetry (xi = x^2, eta = x y).  For the A(u) fed
 * below the degree-1 (affine `linear`) and one-variable (`abaco1_simple`)
 * determining systems are BOTH trivial (verified in the REPL: deg-1 NullSpace
 * empty, every abaco1_simple ratio depends on both variables), so `bivariate` is
 * the only heuristic that can solve them — a declined solve leaves [[1,1]]
 * non-Equal and fails, so the pass is never vacuous.  The pinned builtin returns
 * the implicit first integral, verified by implicit differentiation. */
static void lie_bivariate_ok(const char* A_of_u) {
    char rhs[512], buf[1400];
    snprintf(rhs, sizeof(rhs), "y[x]/x + (%s /. u -> y[x]/x)/x", A_of_u);
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{r = %s, eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> r]]", rhs, rhs);
    ASSERT_TRUE(buf);
}

/* Lie `abaco1_product` forward generator (Cheb-Terrab & Roche 1998, §4.1): every
 * omega = 2 x y/(x^2 + P(y)) has the rational-but-non-polynomial point symmetry
 * (xi = y/x, eta = 0) — from the invariant family f_x/(g f + J) with f = x^2/2,
 * g = 1/y, so the symmetry is the same for any P (= 2 y J(y)).  The y^4 forms of P
 * below break the scaling symmetry AND admit no polynomial symmetry of degree <= 3
 * (verified in the REPL: the deg-1/2/3 determining NullSpace is empty and every
 * abaco1_simple ratio depends on both variables), so `linear`/`bivariate`/
 * `abaco1_simple` all decline — abaco1_product is the only heuristic that can
 * solve them, and a declined solve leaves [[1,1]] non-Equal and fails, so the pass
 * is never vacuous.  The pinned builtin returns the implicit first integral,
 * verified by implicit differentiation. */
static void lie_product_ok(const char* P_of_y) {
    char rhs[512], buf[1400];
    snprintf(rhs, sizeof(rhs), "2 x y[x]/(x^2 + (%s))", P_of_y);
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{r = %s, eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> r]]", rhs, rhs);
    ASSERT_TRUE(buf);
}

/* Chini reducible by construction: from (f, n, B, C) set g = B - f'/((n-1)f),
 * h = C f^(-1/(n-1)) so the scaling y = f^(-1/(n-1)) u gives u' = u^n + B u + C —
 * B, C are the chosen constants, so the autonomous reduction always succeeds. */
static void chini_ok(const char* f, const char* n, const char* B, const char* C) {
    char rhs[700];
    snprintf(rhs, sizeof(rhs),
        "(%s) y[x]^(%s) + ((%s) - D[%s,x]/(((%s)-1) (%s))) y[x] + ((%s) (%s)^(-1/((%s)-1)))",
        f, n, B, f, n, f, C, f, n);
    impl_ok(rhs);
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

/* scalar-factor variable-coefficient coupled system  Y' == f(x) B Y,
 * B = {{a,b},{c,d}} constant (b,c != 0 so it is genuinely non-triangular).
 * DSolve`LinearSystemVarCoeff: tau = Integrate[f,x], Phi = e^{B tau}. */
static void sys2_varcoeff_ok(const char* f, int a, int b, int c, int d) {
    char sys[320], buf[1024];
    snprintf(sys, sizeof(sys),
        "{y'[x] == (%s)((%d) y[x] + (%d) z[x]), z'[x] == (%s)((%d) y[x] + (%d) z[x])}",
        f, a, b, f, c, d);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, {y, z}, x]] === List", sys);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "And @@ (PossibleZeroQ /@ ({y'[x] - (%s)((%d) y[x] + (%d) z[x]), "
        "z'[x] - (%s)((%d) y[x] + (%d) z[x])} /. DSolve[%s, {y, z}, x][[1]]))",
        f, a, b, f, c, d, sys);
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

/* PDEQuasilinear (Lagrange) — variable-coefficient SEMILINEAR class Pc u_x +
 * Qc u_y == Rc (Pc, Qc free of u so pde1 declines; Rc affine in u): the general
 * solution is EXPLICIT with an arbitrary function of the characteristic
 * invariant.  Guard Head === List first, then back-substitute with C[1] pinned
 * to Sin (u[x,y] in Rc is rewritten to the solution too). */
static void pdequasi_semilinear_ok(const char* Pc, const char* Qc, const char* Rc) {
    char eqn[320], buf[1200];
    snprintf(eqn, sizeof(eqn),
        "(%s) D[u[x,y],x] + (%s) D[u[x,y],y] == (%s)", Pc, Qc, Rc);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, u, {x, y}]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "With[{uc = (u[x,y] /. DSolve[%s, u, {x, y}][[1]]) /. C[1][z_] :> Sin[z]}, "
        "PossibleZeroQ[(((%s) D[uc,x] + (%s) D[uc,y] - (%s)) /. u[x,y] -> uc)]]",
        eqn, Pc, Qc, Rc);
    ASSERT_TRUE(buf);
}

/* PDEQuasilinear — genuinely quasilinear CONSERVATION law Pc u_x + Qc u_y == 0
 * (Pc or Qc depends on u): the solution is the IMPLICIT relation G == C[1][u].
 * Guard the branch is an Equal, then verify by implicit differentiation with
 * C[1] pinned to #^2 (u_vi = -Psi_vi / Psi_U, Psi = G - U^2). */
static void pdequasi_conservation_ok(const char* Pc, const char* Qc) {
    char eqn[320], buf[1200];
    snprintf(eqn, sizeof(eqn),
        "(%s) D[u[x,y],x] + (%s) D[u[x,y],y] == 0", Pc, Qc);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, u, {x, y}][[1,1]]] === Equal", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "Module[{rel, Psi, ux, uy}, "
        "rel = DSolve[%s, u, {x,y}][[1,1]] /. C[1] -> (#^2 &); "
        "Psi = (rel[[1]] - rel[[2]]) /. u[x,y] -> U; "
        "ux = -D[Psi,x]/D[Psi,U]; uy = -D[Psi,y]/D[Psi,U]; "
        "PossibleZeroQ[((%s) ux + (%s) uy) /. u[x,y] -> U]]",
        eqn, Pc, Qc);
    ASSERT_TRUE(buf);
}

/* PDEClairaut — u == x u_x + y u_y + f(u_x, u_y) with a NONLINEAR f (given in the
 * placeholder derivatives p = u_x, q = u_y); the complete integral is
 * C[1] x + C[2] y + f(C[1],C[2]).  Guard Head === List, then verify the complete
 * integral back-substitutes (bare constants C[1], C[2] survive verification). */
static void pdeclairaut_ok(const char* f) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
        "Head[DSolve[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] + "
        "((%s) /. {p -> D[u[x,y],x], q -> D[u[x,y],y]}), u, {x,y}]] === List", f);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "With[{uc = u[x,y] /. DSolve[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] + "
        "((%s) /. {p -> D[u[x,y],x], q -> D[u[x,y],y]}), u, {x,y}][[1]]}, "
        "PossibleZeroQ[uc - (x D[uc,x] + y D[uc,y] + ((%s) /. {p -> D[uc,x], q -> D[uc,y]}))]]",
        f, f);
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

static void t_stress_exactode(void) {
    const char* pairs[][2] = {
        {"x", "2"}, {"x", "3"}, {"x", "0"},
        {"x^2", "x"}, {"x^2", "2 x"}, {"1", "1"},
    };
    for (size_t i = 0; i < 6; i++) exact_ode_ok(pairs[i][0], pairs[i][1]);
    /* inhomogeneous: the forcing is folded into the first integral */
    method_ok("DSolve`ExactODE", "x y''[x] + 3 y'[x] == x", "x y''[x] + 3 y'[x] - x");
}

static void t_stress_operfactor(void) {
    /* all-constant spectra (order 3 & 4) */
    operfactor_ok("{1, 2, 3}");
    operfactor_ok("{-1, 1, 2}");
    operfactor_ok("{1, 2, 3, 4}");
    /* shifted-Euler (shared pole x-b, b != 0 — EulerCauchy declines these) */
    operfactor_ok("{1/(x-1), 2/(x-1), 4/(x-1)}");
    operfactor_ok("{1/(x-2), 3/(x-2), 5/(x-2)}");
    operfactor_ok("{2/(x+1), 3/(x+1), 4/(x+1), 6/(x+1)}");
    /* resonant repeated factor (secular Log in the basis) */
    operfactor_ok("{1/(x-1), 1/(x-1), 3/(x-1)}");
    /* Euler at 0 (claimed by EulerCauchy in AUTO, still solved when pinned) */
    operfactor_ok("{1/x, 2/x, 4/x}");
}

static void t_stress_dfactor(void) {
    dfactor_ok("{1, 2, 3}");
    dfactor_ok("{1, 2, 3, 4}");
    dfactor_ok("{1/(x-1), 2/(x-1), 4/(x-1)}");
    dfactor_ok("{2/(x+1), 3/(x+1), 5/(x+1)}");
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

static void t_stress_lagrange(void) {
    /* phi=2p (mu=t^2, rational parametric — fully decidable residual) */
    const char* psis[] = { "p^2", "p^3", "p^4", "p^2 + p", "p^3 - p^2", "2 p^2 - p" };
    for (size_t i = 0; i < 6; i++) lagrange_ok("2 p", psis[i]);
    /* phi=3p/2 (mu=t^3, rational) */
    lagrange_ok("3 p/2", "p^2");
    lagrange_ok("3 p/2", "p^3");
}

static void t_stress_chini(void) {
    /* (f, n, B, C) reducible by construction -> implicit first integral */
    chini_ok("x^2", "3", "0", "1");
    chini_ok("x^2", "3", "1", "-2");
    chini_ok("x^3", "4", "0", "1");
    chini_ok("x^4", "3", "0", "1");
    /* Abel: the f=x^2 Chini shifted by z = y + 1 (introduces the y^2 term) */
    impl_ok("x^2 y[x]^3 + 3 x^2 y[x]^2 + (3 x^2 - 1/x) y[x] + x^2");
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

static void t_stress_linsys_varcoeff(void) {
    /* scalar factors f(x) with elementary antiderivative (Log[x], x^2/2, ...) x
     * genuinely-coupled constant matrices B (b,c != 0): distinct-real, complex,
     * distinct-real spectra.  Each is A(x)=f(x)B, non-triangular, variable-coeff. */
    const char* fs[] = { "1/x", "x", "2/x" };
    int m[][4] = {
        {2, 1, 1, 2},   /* eigenvalues 1, 3 */
        {0, 1, -1, 0},  /* eigenvalues +- i (complex spectrum -> Cos/Sin[tau]) */
        {1, 2, 2, 1}    /* eigenvalues 3, -1 */
    };
    for (size_t k = 0; k < 3; k++)
        for (size_t i = 0; i < 3; i++)
            sys2_varcoeff_ok(fs[k], m[i][0], m[i][1], m[i][2], m[i][3]);
}

static void t_stress_pde(void) {
    pde_ok(1, 1, "0");
    pde_ok(2, 3, "0");
    pde_ok(1, 2, "y");
    pde_ok(3, 5, "x");
}

/* M6: first-order nonlinear PDEs — quasilinear (Lagrange) + Clairaut. */
static void t_stress_pde_quasilinear(void) {
    /* variable-coefficient semilinear (pde1 declines; quasilinear owns them) */
    pdequasi_semilinear_ok("x", "y", "u[x,y]");    /* x C[1][y/x]          */
    pdequasi_semilinear_ok("x", "-y", "0");        /* C[1][x y]            */
    pdequasi_semilinear_ok("y", "x", "0");         /* C[1][y^2 - x^2]      */
    pdequasi_semilinear_ok("2 x", "1", "0");       /* C[1][y - Log[x]/2]   */
    pdequasi_semilinear_ok("x", "y", "2 u[x,y]");  /* c = 2                */
    pdequasi_semilinear_ok("1", "x", "y");         /* y along characteristic */
    /* genuinely quasilinear conservation laws (implicit) */
    pdequasi_conservation_ok("u[x,y]", "1");       /* inviscid Burgers     */
    pdequasi_conservation_ok("u[x,y]", "-1");
    pdequasi_conservation_ok("1", "u[x,y]");
    pdequasi_conservation_ok("u[x,y]", "x");
}

static void t_stress_pde_clairaut(void) {
    pdeclairaut_ok("p q");          /* u_x u_y             */
    pdeclairaut_ok("p^2 + q^2");    /* u_x^2 + u_y^2       */
    pdeclairaut_ok("p^2");          /* u_x^2               */
    pdeclairaut_ok("p q + p");      /* mixed nonlinear     */
    pdeclairaut_ok("p^2 + q");      /* nonlinear + linear  */
}

/* Factorable: (y' - r1 y)(y' - r2 y) == 0 factors into two linear ODEs; EVERY
 * branch must back-substitute (not just [[1]]), so verify with And @@ Map. */
static void factorable_ok(const char* r1, const char* r2) {
    char eqn[512], buf[1024];
    snprintf(eqn, sizeof(eqn),
             "(y'[x] - (%s) y[x]) (y'[x] - (%s) y[x])", r1, r2);
    snprintf(buf, sizeof(buf), "Head[DSolve[%s == 0, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[(%s) /. #] &, DSolve[%s == 0, y, x]]", eqn, eqn);
    ASSERT_TRUE(buf);
}
static void t_stress_factorable(void) {
    const char* rs[][2] = { {"1","-1"}, {"2","3"}, {"1","0"}, {"-2","5"}, {"1/2","-3"} };
    for (size_t i = 0; i < sizeof(rs)/sizeof(rs[0]); i++) factorable_ok(rs[i][0], rs[i][1]);
}

/* NthAlgebraic: (y')^2 == c (c > 0 constant) -> y' == +/- Sqrt[c], both branches
 * free of y (Quadrature).  Every branch back-substitutes. */
static void nthalg_ok(const char* c) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`NthAlgebraic[(y'[x])^2 == %s, y, x]] === List", c);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[((y'[x])^2 - (%s)) /. #] &, "
             "DSolve`NthAlgebraic[(y'[x])^2 == %s, y, x]]", c, c);
    ASSERT_TRUE(buf);
}
static void t_stress_nth_algebraic(void) {
    const char* cs[] = { "2", "9", "5", "1/4" };
    for (size_t i = 0; i < sizeof(cs)/sizeof(cs[0]); i++) nthalg_ok(cs[i]);
    /* the flagship y-dependent case: (y')^2 == 4 y -> (x + C)^2 */
    char buf[512];
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[((y'[x])^2 - 4 y[x]) /. #] &, "
             "DSolve`NthAlgebraic[(y'[x])^2 == 4 y[x], y, x]]");
    ASSERT_TRUE(buf);
}

/* LinearCoefficients: y' == (a1 x+b1 y+c1)/(a2 x+b2 y+c2), det != 0, real-root
 * (explicit) cases -> every branch back-substitutes to the cleared equation.
 * (Real roots require (a2-b1)^2 + 4 a1 b2 >= 0; the curated tuples all satisfy it.) */
static void lincoeff_ok(int a1,int b1,int c1,int a2,int b2,int c2) {
    char eqn[512], resid[512], buf[1200];
    snprintf(eqn, sizeof(eqn),
             "y'[x] == (%d x + %d y[x] + %d)/(%d x + %d y[x] + %d)", a1,b1,c1,a2,b2,c2);
    snprintf(resid, sizeof(resid),
             "y'[x] (%d x + %d y[x] + %d) - (%d x + %d y[x] + %d)", a2,b2,c2,a1,b1,c1);
    snprintf(buf, sizeof(buf), "Head[DSolve`LinearCoefficients[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[(%s) /. #] &, DSolve`LinearCoefficients[%s, y, x]]",
             resid, eqn);
    ASSERT_TRUE(buf);
}
static void t_stress_lincoeff(void) {
    lincoeff_ok(1,2,-4, 2,1,-5);   /* lines meet at (2,1) */
    lincoeff_ok(2,3,-1, 3,2, 2);
    lincoeff_ok(1,4,-1, 4,1,-1);
}

/* AlmostLinear: build 2 y y' + P(x) y^2 - Q(x) == 0 (u = y^2, u' + P u == Q);
 * every explicit branch back-substitutes. */
static void almostlinear_ok(const char* p, const char* q) {
    char eqn[512], resid[512], buf[1200];
    snprintf(eqn, sizeof(eqn), "2 y[x] y'[x] + (%s) y[x]^2 - (%s) == 0", p, q);
    snprintf(resid, sizeof(resid), "2 y[x] y'[x] + (%s) y[x]^2 - (%s)", p, q);
    snprintf(buf, sizeof(buf), "Head[DSolve`AlmostLinear[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[(%s) /. #] &, DSolve`AlmostLinear[%s, y, x]]", resid, eqn);
    ASSERT_TRUE(buf);
}
static void t_stress_almostlinear(void) {
    const char* ps[] = { "1", "2" };
    const char* qs[] = { "x", "1", "Exp[x]" };
    for (size_t i = 0; i < 2; i++) for (size_t j = 0; j < 3; j++) almostlinear_ok(ps[i], qs[j]);
}

/* SeparableReduced: y' == y^2/(1 + c x y) has r = x y/(1 + c x y) = w/(1+c w),
 * w = x y (n=1); returns the implicit first integral, verified by implicit diff. */
static void sepreduced_ok(const char* c) {
    char rhs[256], buf[1200];
    snprintf(rhs, sizeof(rhs), "y[x]^2/(1 + (%s) x y[x])", c);
    snprintf(buf, sizeof(buf),
             "Head[DSolve`SeparableReduced[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[Module[{eq = DSolve`SeparableReduced[y'[x] == %s, y, x][[1,1]]}, "
             "D[eq[[1]] - eq[[2]], x] /. y'[x] -> (%s)]]", rhs, rhs);
    ASSERT_TRUE(buf);
}
static void t_stress_sepreduced(void) {
    const char* cs[] = { "1", "2", "3" };
    for (size_t i = 0; i < 3; i++) sepreduced_ok(cs[i]);
}

/* Liouville: y'' + g(y)(y')^2 + h(x) y' == 0 for chosen g(y), h(x) with elementary
 * EG = Integrate[Exp[Integrate[g,y]],y] and EH; every explicit branch verifies. */
static void liouville_ok(const char* g, const char* h) {
    char eqn[512], resid[512], buf[1200];
    snprintf(eqn, sizeof(eqn), "y''[x] + (%s)(y'[x])^2 + (%s) y'[x] == 0", g, h);
    snprintf(resid, sizeof(resid), "y''[x] + (%s)(y'[x])^2 + (%s) y'[x]", g, h);
    snprintf(buf, sizeof(buf), "Head[DSolve`Liouville[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[(%s) /. #] &, DSolve`Liouville[%s, y, x]]", resid, eqn);
    ASSERT_TRUE(buf);
}
static void t_stress_liouville(void) {
    liouville_ok("1/y[x]", "1/x");
    liouville_ok("1/y[x]", "1");
    liouville_ok("2/y[x]", "1/x");
}

/* UndeterminedCoefficients: y'' + c1 y' + c0 y == f for UC forcing f (some cases
 * resonant); the single-branch solution back-substitutes. */
static void undetcoeff_ok(const char* c1, const char* c0, const char* f) {
    char eqn[512], resid[512], buf[1200];
    snprintf(eqn, sizeof(eqn), "y''[x] + (%s) y'[x] + (%s) y[x] == %s", c1, c0, f);
    snprintf(resid, sizeof(resid), "y''[x] + (%s) y'[x] + (%s) y[x] - (%s)", c1, c0, f);
    snprintf(buf, sizeof(buf), "Head[DSolve`UndeterminedCoefficients[%s, y, x]] === List", eqn);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[(%s) /. DSolve`UndeterminedCoefficients[%s, y, x][[1]]]", resid, eqn);
    ASSERT_TRUE(buf);
}
static void t_stress_undetcoeff(void) {
    const char* ops[][2] = { {"0","1"}, {"0","-1"}, {"3","2"}, {"-2","1"} };
    const char* fs[] = { "x", "x^2", "Exp[x]", "Sin[x]", "Cos[2 x]", "x Exp[x]" };
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 6; j++)
            undetcoeff_ok(ops[i][0], ops[i][1], fs[j]);
}

/* FirstOrderPowerSeries: y' == F for F analytic at 0; the truncated residual is
 * O[x]^N, so Normal[residual] == 0. */
static void fops_ok(const char* rhs) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`FirstOrderPowerSeries[y'[x] == %s, y, x]] === List", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[Normal[(y'[x] - (%s)) /. "
             "DSolve`FirstOrderPowerSeries[y'[x] == %s, y, x][[1]]]]", rhs, rhs);
    ASSERT_TRUE(buf);
}
static void t_stress_fops(void) {
    const char* rhss[] = { "x + y[x]", "x + y[x]^2", "x^2 + y[x]^3",
                           "Sin[x] + y[x]", "x y[x] + 1", "y[x]^2 - x" };
    for (size_t i = 0; i < sizeof(rhss)/sizeof(rhss[0]); i++) fops_ok(rhss[i]);
}
static void t_stress_lie_bivariate(void) {
    const char* As[] = { "u^2 + u", "u^2 - 1", "u^2 + 2", "u^2 - u",
                         "u^2 + u + 1", "2 u^2 - 3", "u^3 + u", "u^2 + 3 u" };
    for (size_t i = 0; i < sizeof(As)/sizeof(As[0]); i++) lie_bivariate_ok(As[i]);
}
static void t_stress_lie_product(void) {
    const char* Ps[] = { "2 y[x]^4 + 2", "2 y[x]^4 + 2 y[x]", "2 y[x]^4 - 3",
                         "2 y[x]^4 + y[x] + 1", "2 y[x]^4 - 2 y[x]", "3 y[x]^4 + 1",
                         "2 y[x]^4 + 3 y[x] - 1" };
    for (size_t i = 0; i < sizeof(Ps)/sizeof(Ps[0]); i++) lie_product_ok(Ps[i]);
}
/* abaco2_similar (§4.3): [F(x), H(x)] symmetry.  y' == (a x + b y + c)^p with p a
 * non-integer power is the clean isolating family — omega is irrational, so all the
 * rational heuristics (abaco1_simple/linear/abaco1_product) decline and abaco2_similar
 * fires (Q = omega_y/omega_yy, T = Q_x/Q_y free of y).  Verified implicit-diff residual. */
static void lie_similar_ok(const char* rhs) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{r = %s, eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> r]]", rhs, rhs);
    ASSERT_TRUE(buf);
}
static void t_stress_lie_similar(void) {
    const char* rhss[] = { "Sqrt[x + y[x]]", "Sqrt[2 x + y[x]]", "Sqrt[x + 2 y[x]]",
                           "Sqrt[x + y[x] + 1]", "Sqrt[3 x - y[x]]", "Sqrt[x - y[x] + 2]",
                           "(x + y[x])^(1/3)", "Sqrt[4 x + 3 y[x] + 2]" };
    for (size_t i = 0; i < sizeof(rhss)/sizeof(rhss[0]); i++) lie_similar_ok(rhss[i]);
}
/* function_sum (§4.2): the additive symmetry [F(x)+G(y), 0].  Each omega is a
 * pre-constructed member of the invariant family (F=1/x/2/x, G=y/2y with J=0/1/y..),
 * whose 1/omega is transcendental (Log) so the rational heuristics decline and
 * function_sum solves it (attribution-verified).  Implicit-diff residual == 0. */
static void lie_fsum_ok(const char* rhs) {
    char buf[1400];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{r = %s, eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> r]]", rhs, rhs);
    ASSERT_TRUE(buf);
}
static void t_stress_lie_function_sum(void) {
    const char* rhss[] = {
        "(x y[x]^3)/(-1 + x y[x] + x^2 y[x]^2 - 2 Log[1 + x y[x]] - 2 x y[x] Log[1 + x y[x]])",
        "(x y[x]^3)/(-1 + y[x]^2 + x y[x]^3 + x y[x] + x^2 y[x]^2 - 2 Log[1 + x y[x]] - 2 x y[x] Log[1 + x y[x]])",
        "(x y[x]^3)/(-1 + x y[x]^2 + x y[x] + x^2 y[x]^2 + y[x] - 2 Log[1 + x y[x]] - 2 x y[x] Log[1 + x y[x]])",
        "(x y[x]^3)/(-4 + 2 x y[x] + x^2 y[x]^2 - 8 Log[2 + x y[x]] - 4 x y[x] Log[2 + x y[x]])",
        "(4 x y[x]^3)/(-1 + 2 x y[x] + 4 x^2 y[x]^2 - 4 x y[x] Log[1 + 2 x y[x]] - 2 Log[1 + 2 x y[x]])" };
    for (size_t i = 0; i < sizeof(rhss)/sizeof(rhss[0]); i++) lie_fsum_ok(rhss[i]);
}
/* abaco2_unique_unknown (§4.4.1): [F(x),G(y)]/[G(y),F(x)] from a non-integer power of
 * both variables.  omega = (c x/y)(a x^2 + b y^2 + d)^p has symmetry [1/x, -1/y]; the
 * irrational power makes the rational heuristics and abaco2_similar decline
 * (attribution-verified).  Implicit-diff residual == 0. */
static void lie_unique_ok(const char* rhs) {
    char buf[1200];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    ASSERT_TRUE(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{r = %s, eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]]-eq[[2]], x] /. y'[x] -> r]]", rhs, rhs);
    ASSERT_TRUE(buf);
}
static void t_stress_lie_unique_unknown(void) {
    const char* rhss[] = {
        "(x/y[x]) (x^2 + y[x]^2)^(1/3)",   "(x/y[x]) Sqrt[x^2 + y[x]^2]",
        "(x/y[x]) (2 x^2 + y[x]^2)^(1/3)", "(x/y[x]) (x^2 + 3 y[x]^2)^(1/3)",
        "(2 x/y[x]) (x^2 + y[x]^2)^(1/3)", "(x/y[x]) (x^2 + y[x]^2 + 1)^(1/3)" };
    for (size_t i = 0; i < sizeof(rhss)/sizeof(rhss[0]); i++) lie_unique_ok(rhss[i]);
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();   /* match production: deriv.m rules + CRC integral tables */

    TEST(t_stress_factorable);
    TEST(t_stress_nth_algebraic);
    TEST(t_stress_lincoeff);
    TEST(t_stress_almostlinear);
    TEST(t_stress_sepreduced);
    TEST(t_stress_liouville);
    TEST(t_stress_undetcoeff);
    TEST(t_stress_fops);
    TEST(t_stress_lie_bivariate);
    TEST(t_stress_lie_product);
    TEST(t_stress_lie_similar);
    TEST(t_stress_lie_function_sum);
    TEST(t_stress_lie_unique_unknown);
    TEST(t_stress_linear1);
    TEST(t_stress_separable);
    TEST(t_stress_bernoulli);
    TEST(t_stress_homogeneous);
    TEST(t_stress_homogeneous_implicit);
    TEST(t_stress_exact);
    TEST(t_stress_constcoeff);
    TEST(t_stress_euler);
    TEST(t_stress_exactode);
    TEST(t_stress_operfactor);
    TEST(t_stress_dfactor);
    TEST(t_stress_reduce_order);
    TEST(t_stress_reduce_order_riccati);
    TEST(t_stress_riccati);
    TEST(t_stress_lagrange);
    TEST(t_stress_chini);
    TEST(t_stress_systems);
    TEST(t_stress_triangular);
    TEST(t_stress_linsys_varcoeff);
    TEST(t_stress_pde);
    TEST(t_stress_pde_quasilinear);
    TEST(t_stress_pde_clairaut);

    printf("\nAll DSolve stress tests passed.\n");
    return 0;
}
