/* integrate_elliptic.c
 *
 * Definite integrals Int_a^b R(x, Sqrt[P(x)]) dx, P a squarefree cubic or
 * quartic, in closed form through EllipticF / EllipticE / EllipticPi /
 * EllipticK (parameter m = k^2) plus elementary terms.  See
 * integrate_elliptic.h for the overview; the pipeline, stage by stage:
 *
 * 1. NORMALISE.  Every Power[g, k/2] (g a rational function of x) is written
 *    on the interval as (kappa S/D)^k Sqrt[P]^k, where P > 0 is the common
 *    squarefree kernel of all the radicands (decline when two radicands have
 *    different kernels), S the square part and D the denominator of g, and
 *    kappa a constant (+-1, +-I, a square root of a ratio of contents).  The
 *    integrand is then rational in x and y = Sqrt[P]; rationalising the
 *    denominator gives f = A(x) + W(x)/y.  A is integrated by Integrate
 *    (rational), W/Sqrt[P] by the rest of the pipeline.
 *
 * 2. EVEN FORM.  P = pi Q1 Q2 with Q1, Q2 real binary quadratic forms whose
 *    root pairs do not interlace (a cubic has the point Infinity as its fourth
 *    root: one Q is linear).  The two roots mu, nu of the Jacobian of
 *    (Q1, Q2) are the common harmonic pair; the Moebius map
 *    x = (al t + be)/(ga t + de) with x(0) = mu, x(Infinity) = nu makes both
 *    quadratics even: (ga t + de)^2 Q_i(x(t)) = A_i + B_i t^2, and
 *    dx/Sqrt[P] = (al de - be ga) dt/Sqrt[Phat(t)],
 *    Phat = pi (A1 + B1 t^2)(A2 + B2 t^2).  The interval is split so that nu
 *    is never inside a piece (t stays finite).
 *
 * 3. HERMITE.  W's partial fractions in x map term by term to partial
 *    fractions in t (closed-form binomial expansions; the pole of x^n is the
 *    image of x = Infinity).  Hermite reduction against Phat removes every
 *    pole of order >= 2 (and every pole at a root of Phat) and every power
 *    t^n, n >= 3, leaving T0 = Int dt/Sqrt[Phat], T1 = Int t dt/.., T2 =
 *    Int t^2 dt/.., J_tau = Int dt/((t - tau) Sqrt[Phat]) and an algebraic
 *    part c(t) Sqrt[Phat] evaluated at the endpoints.
 *
 * 4. CARLSON.  On a t-range of one sign, s = t^2: the even parts (T0, T2 and
 *    tau/(t^2 - tau^2) of J_tau) become Int g(s) ds / Sqrt[s L1 L2] with
 *    L_i = A_i + B_i s REAL linear factors, i.e. the cubic case of DLMF
 *    19.29.4 (R_F), 19.29.7 (R_D) and 19.29.8 (R_J + R_C); the odd parts
 *    (T1, t/(t^2 - tau^2)) become quadratic radicals (R_C).  The image of the
 *    pole's mirror point -tau may lie in the range: the even and odd parts
 *    then each carry a removable-in-sum simple pole, and both are taken as
 *    Cauchy principal values -- 19.29.8 is valid as a PV with the sign of the
 *    R_C term fixed by sign(X1 X2 X3 (s1 - q) + Y1 Y2 Y3 (s2 - q)) (checked
 *    against quadrature on random configurations), and the odd part by the
 *    explicit logarithm.
 *
 * 5. LEGENDRE.  R_F, R_D, R_J with three distinct positive arguments
 *    u1 < u2 < u3 go to F, E, Pi with lambda = u3 - u1, Sin[phi]^2 =
 *    lambda/u3, m = (u3 - u2)/lambda (DLMF 19.25.5-19.25.14: R_D's special
 *    argument may be the largest, middle or smallest, three formulas); a
 *    negative R_J parameter (principal value) goes through DLMF 19.20.14;
 *    R_C is ArcTan / ArcTanh (19.2.17-19.2.20).
 *
 * 6. VERIFY.  N[result, 30] against a 90-digit double-exponential quadrature
 *    of f (ell_reference); relative disagreement above 10^-20 declines.
 *
 * All exact arithmetic is done by the evaluator through small templates
 * (`ell_x`), so the radicals stay exact; numeric values are used only to make
 * DECISIONS (orderings, signs, zero tests at 60 digits), and the final
 * numeric verification guards every one of them.
 */

#include "integrate_elliptic.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* arith_warnings_mute_push/pop */
#include "message.h"
#include "sym_names.h"
#include "core.h"         /* tc_register_unwind_counter */

#include <math.h>
#include <stdio.h>
#include "print.h"
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define ELL_MAXK   10     /* highest pole order handled            */
#define ELL_MAXP   20     /* polynomial coefficient slots (degree) */
#define ELL_MAXPOL 12     /* distinct poles                         */
#define ELL_MAXALG 64     /* algebraic-part terms                   */

static int ell_depth = 0;

/* MATHILDA_ELL_DEBUG=1 traces the decline points on stderr (development aid;
 * not a user-facing diagnostic). */
static void ell_dbg(const char* m) { if (getenv("MATHILDA_ELL_DEBUG")) fprintf(stderr, "[ell] %s\n", m); }
#define ELL_DBG(msg) do { if (getenv("MATHILDA_ELL_DEBUG")) fputs("[ell] " msg "\n", stderr); } while (0)

/* -------------------------------------------------------------------------
 * Evaluator plumbing.
 * ---------------------------------------------------------------------- */

/* Parsed template cache (keyed by the template's address: every template is
 * a string literal).  The cached bodies live for the session. */
typedef struct { const char* key; Expr* body; } EllTmpl;
static EllTmpl ell_cache[160];
static int ell_ncache = 0;

static Expr* ell_body(const char* t) {
    for (int i = 0; i < ell_ncache; i++)
        if (ell_cache[i].key == t) return ell_cache[i].body;
    Expr* b = parse_expression(t);
    if (!b) return NULL;
    if (ell_ncache < (int)(sizeof ell_cache / sizeof ell_cache[0])) {
        ell_cache[ell_ncache].key = t;
        ell_cache[ell_ncache].body = b;
        ell_ncache++;
        return b;
    }
    return b;   /* cache full: leaks one parse per call; never reached */
}

/* Evaluate the template `t` (a pure-function body in #1..#n) applied to the
 * borrowed arguments.  Returns an owned result (never NULL unless parsing
 * fails). */
static Expr* ell_x(const char* t, int n, ...) {
    Expr* body = ell_body(t);
    if (!body) return NULL;
    Expr* args[12];
    va_list ap;
    va_start(ap, n);
    for (int i = 0; i < n && i < 12; i++) args[i] = expr_copy(va_arg(ap, Expr*));
    va_end(ap);
    Expr* fn = expr_new_function(expr_new_symbol("Function"),
                                 (Expr*[]){ expr_copy(body) }, 1);
    Expr* call = expr_new_function(fn, n ? args : NULL, (size_t)n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static bool ell_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool ell_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && ell_is_sym(e->data.function.head, name);
}
static bool ell_true(Expr* e) {            /* consumes e */
    bool r = ell_is_sym(e, "True");
    expr_free(e);
    return r;
}
static size_t ell_len(const Expr* e) {
    return ell_head_is(e, "List") ? e->data.function.arg_count : 0;
}
static Expr* ell_arg(const Expr* e, size_t i) { return e->data.function.args[i]; }

/* Machine value of an exact numeric expression; NAN when it is not a real
 * number (|Im| relatively above 1e-12). */
static double ell_num(Expr* e) {
    Expr* r = ell_x("N[N[{Re[#1], Im[#1]}, 40]]", 1, e);
    double v = NAN;
    if (ell_len(r) == 2) {
        Expr* re = ell_arg(r, 0);
        Expr* im = ell_arg(r, 1);
        double a = re->type == EXPR_REAL ? re->data.real
                 : re->type == EXPR_INTEGER ? (double)re->data.integer : NAN;
        double b = im->type == EXPR_REAL ? im->data.real
                 : im->type == EXPR_INTEGER ? (double)im->data.integer : NAN;
        if (!isnan(a) && !isnan(b) && fabs(b) <= 1e-12 * (1.0 + fabs(a))) v = a;
    }
    expr_free(r);
    return v;
}

/* Exact zero test for an algebraic number, decided at 60 digits. */
static bool ell_zero(Expr* e) {
    return ell_true(ell_x("#1 === 0 || TrueQ[Abs[N[#1, 60]] < 10^-40]", 1, e));
}

/* Owned arithmetic on borrowed operands. */
static Expr* E_add(Expr* a, Expr* b) { return ell_x("#1 + #2", 2, a, b); }
static Expr* E_sub(Expr* a, Expr* b) { return ell_x("#1 - #2", 2, a, b); }
static Expr* E_int(long v)           { return expr_new_integer((int64_t)v); }

/* acc += v (acc may be NULL = 0); consumes v. */
static void acc_add(Expr** acc, Expr* v) {
    if (!v) return;
    if (!*acc) { *acc = v; return; }
    Expr* s = E_add(*acc, v);
    expr_free(*acc); expr_free(v);
    *acc = s;
}

static bool contains_sym(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == x->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_sym(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_sym(e->data.function.args[i], x)) return true;
    return false;
}

static double binom(int n, int k) {
    double r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return r;
}
static Expr* E_binom(int n, int k) { return E_int((long)(binom(n, k) + 0.5)); }

/* -------------------------------------------------------------------------
 * Points of the extended real line.
 * ---------------------------------------------------------------------- */
typedef struct { Expr* v; int inf; double d; } Pt;   /* inf: -1, 0, +1 */

static bool pt_from(Expr* e, Pt* p) {
    p->v = NULL; p->inf = 0; p->d = 0;
    if (ell_true(ell_x("#1 === Infinity", 1, e)))  { p->inf = 1;  p->d = INFINITY;  return true; }
    if (ell_true(ell_x("#1 === -Infinity", 1, e))) { p->inf = -1; p->d = -INFINITY; return true; }
    double d = ell_num(e);
    if (isnan(d)) return false;
    p->v = expr_copy(e); p->d = d;
    return true;
}
static void pt_free(Pt* p) { if (p->v) expr_free(p->v); p->v = NULL; }

/* -------------------------------------------------------------------------
 * Carlson -> Legendre (exact).  Arguments are borrowed, positive (except the
 * R_J parameter / R_C second argument, which may be negative: principal
 * value), and the result is owned.  NULL on a degenerate configuration.
 * ---------------------------------------------------------------------- */

typedef struct { Expr* u[3]; double d[3]; } Sorted3;

static bool sort3(Expr* x, Expr* y, Expr* z, Sorted3* s) {
    Expr* in[3] = { x, y, z };
    double dv[3];
    for (int i = 0; i < 3; i++) {
        dv[i] = ell_num(in[i]);
        if (isnan(dv[i]) || dv[i] < -1e-30) return false;
    }
    int idx[3] = { 0, 1, 2 };
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            if (dv[idx[j]] < dv[idx[i]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    for (int i = 0; i < 3; i++) { s->u[i] = in[idx[i]]; s->d[i] = dv[idx[i]]; }
    /* three distinct arguments: the Carlson arguments of a squarefree P are
     * always distinct (U_ab^2 - U_ac^2 = d_ad d_bc != 0) */
    for (int i = 0; i < 2; i++) {
        Expr* df = E_sub(s->u[i + 1], s->u[i]);
        bool z = ell_zero(df);
        expr_free(df);
        if (z) return false;
    }
    return true;
}

static Expr* leg_RC(Expr* x, Expr* y) {
    double xd = ell_num(x), yd = ell_num(y);
    if (isnan(xd) || isnan(yd) || xd < -1e-30) return NULL;
    if (yd < 0 || ell_zero(y)) return NULL;   /* PV: see leg_RCpv */
    Expr* d = E_sub(x, y);
    Expr* r;
    if (ell_zero(d))           r = ell_x("1/Sqrt[#1]", 1, x);
    else if (ell_zero(x))      r = ell_x("Pi/(2 Sqrt[#1])", 1, y);
    else if (xd < yd)          r = ell_x("ArcTan[Sqrt[Integrate`Elliptic`Private`simp[(#2 - #1)/#1]]]/Sqrt[Integrate`Elliptic`Private`simp[#2 - #1]]", 2, x, y);
    else                       r = ell_x("ArcTanh[Sqrt[Integrate`Elliptic`Private`simp[(#1 - #2)/#1]]]/Sqrt[Integrate`Elliptic`Private`simp[#1 - #2]]", 2, x, y);
    expr_free(d);
    return r;
}

/* R_C with a principal-value branch: wraps leg_RC. */
static Expr* leg_RCpv(Expr* x, Expr* y) {
    double yd = ell_num(y);
    if (isnan(yd)) return NULL;
    if (yd > 0) return leg_RC(x, y);
    Expr* xm = ell_x("Integrate`Elliptic`Private`simp[#1 - #2]", 2, x, y);
    Expr* my = ell_x("-#1", 1, y);
    Expr* inner = leg_RC(xm, my);
    Expr* r = inner ? ell_x("Sqrt[Integrate`Elliptic`Private`simp[#1/#2]] #3", 3, x, xm, inner) : NULL;
    expr_free(xm); expr_free(my); if (inner) expr_free(inner);
    return r;
}

/* Legendre parameters for sorted u1 < u2 < u3: lambda, m, phi. */
typedef struct { Expr* lam; Expr* m; Expr* phi; bool complete; } LegPar;

static void legpar(const Sorted3* s, LegPar* L) {
    L->lam = ell_x("Integrate`Elliptic`Private`simp[#2 - #1]", 2, s->u[0], s->u[2]);
    L->m = ell_x("Integrate`Elliptic`Private`simp[(#3 - #2)/(#3 - #1)]", 3, s->u[0], s->u[1], s->u[2]);
    L->complete = ell_zero(s->u[0]);
    L->phi = L->complete ? ell_x("Pi/2", 0)
                         : ell_x("ArcSin[Sqrt[Integrate`Elliptic`Private`simp[(#2 - #1)/#2]]]", 2, s->u[0], s->u[2]);
}
static void legpar_free(LegPar* L) { expr_free(L->lam); expr_free(L->m); expr_free(L->phi); }

static Expr* leg_F(LegPar* L) {
    return L->complete ? ell_x("EllipticK[#1]", 1, L->m) : ell_x("EllipticF[#1, #2]", 2, L->phi, L->m);
}
static Expr* leg_E(LegPar* L) {
    return L->complete ? ell_x("EllipticE[#1]", 1, L->m) : ell_x("EllipticE[#1, #2]", 2, L->phi, L->m);
}

static Expr* leg_RF(Expr* x, Expr* y, Expr* z) {
    Sorted3 s;
    if (!sort3(x, y, z, &s)) return NULL;
    LegPar L; legpar(&s, &L);
    Expr* F = leg_F(&L);
    Expr* r = ell_x("#1/Sqrt[#2]", 2, F, L.lam);
    expr_free(F); legpar_free(&L);
    return r;
}

static Expr* leg_RD(Expr* x, Expr* y, Expr* z) {
    Sorted3 s;
    if (!sort3(x, y, z, &s)) return NULL;
    LegPar L; legpar(&s, &L);
    Expr* F = leg_F(&L);
    Expr* E = leg_E(&L);
    Expr* r;
    if (s.u[2] == z)        /* special argument is the largest */
        r = ell_x("3 (#1 - #2)/(#3 #4^(3/2))", 4, F, E, L.m, L.lam);
    else if (s.u[1] == z)   /* middle: DLMF 19.25.11 */
        r = ell_x("3 (#2 - (1 - #3) #1 - #3 Sqrt[#5 #4/(#7 #6)])/(#3 (1 - #3) #4^(3/2))",
                  7, F, E, L.m, L.lam, s.u[0], s.u[1], s.u[2]);
    else                    /* smallest: DLMF 19.25.12 */
        r = ell_x("3 (Sqrt[#5 #3/(#6 #4)] - #1)/((1 - #2) #3^(3/2))",
                  6, E, L.m, L.lam, s.u[0], s.u[1], s.u[2]);
    expr_free(F); expr_free(E); legpar_free(&L);
    return r;
}

static Expr* leg_RJ(Expr* x, Expr* y, Expr* z, Expr* p) {
    double pd = ell_num(p);
    if (isnan(pd) || ell_zero(p)) return NULL;
    Sorted3 s;
    if (!sort3(x, y, z, &s)) return NULL;
    if (pd < 0) {
        /* DLMF 19.20.14 (principal value): with q = -p and
         * pp = (z (x + y + q) - x y)/(z + q) > 0,
         * (q + z) R_J(x,y,z,-q) = (pp - z) R_J(x,y,z,pp) - 3 R_F(x,y,z)
         *     + 3 Sqrt[x y z/(x y + pp q)] R_C(x y + pp q, pp q). */
        Expr* q  = ell_x("-#1", 1, p);
        Expr* pp = ell_x("Integrate`Elliptic`Private`simp[(#3 (#1 + #2 + #4) - #1 #2)/(#3 + #4)]", 4, s.u[0], s.u[1], s.u[2], q);
        Expr* c1 = ell_x("Integrate`Elliptic`Private`simp[#1 #2 + #3 #4]", 4, s.u[0], s.u[1], pp, q);
        Expr* c2 = ell_x("Integrate`Elliptic`Private`simp[#1 #2]", 2, pp, q);
        Expr* rjp = NULL;
        Expr* dz = E_sub(pp, s.u[2]);
        bool ppz = ell_zero(dz);
        if (!ppz) rjp = leg_RJ(s.u[0], s.u[1], s.u[2], pp);
        Expr* rf = leg_RF(s.u[0], s.u[1], s.u[2]);
        Expr* rc = leg_RC(c1, c2);
        Expr* r = NULL;
        if ((ppz || rjp) && rf && rc)
            r = ppz ? ell_x("(-3 #1 + 3 Sqrt[#2 #3 #4/#5] #6)/(#7 + #4)", 7,
                            rf, s.u[0], s.u[1], s.u[2], c1, rc, q)
                    : ell_x("(#1 #2 - 3 #3 + 3 Sqrt[#4 #5 #6/#7] #8)/(#9 + #6)", 9,
                            dz, rjp, rf, s.u[0], s.u[1], s.u[2], c1, rc, q);
        expr_free(q); expr_free(pp); expr_free(c1); expr_free(c2); expr_free(dz);
        if (rjp) expr_free(rjp);
        if (rf) expr_free(rf);
        if (rc) expr_free(rc);
        return r;
    }
    Expr* dz = E_sub(p, s.u[2]);
    if (ell_zero(dz)) { expr_free(dz); return leg_RD(s.u[0], s.u[1], s.u[2]); }
    expr_free(dz);
    for (int i = 0; i < 2; i++) {          /* p equal to another argument */
        Expr* d = E_sub(p, s.u[i]);
        bool z0 = ell_zero(d);
        expr_free(d);
        if (z0) return NULL;
    }
    LegPar L; legpar(&s, &L);
    Expr* n = ell_x("Integrate`Elliptic`Private`simp[(#1 - #2)/#3]", 3, s.u[2], p, L.lam);
    Expr* F = leg_F(&L);
    Expr* P = L.complete ? ell_x("EllipticPi[#1, #2]", 2, n, L.m)
                         : ell_x("EllipticPi[#1, #2, #3]", 3, n, L.phi, L.m);
    Expr* r = ell_x("3 (#1 - #2)/(#3 #4^(3/2))", 4, P, F, n, L.lam);
    expr_free(n); expr_free(F); expr_free(P); legpar_free(&L);
    return r;
}

/* -------------------------------------------------------------------------
 * s-space integrals.  The radicand is Cpi * s^e * L1 * L2 with e = 1 (even
 * part, cubic) or e = 0 (odd part, quadratic), L_i = A_i + B_i s, and the
 * integral runs over [s1, s2] (oriented: s1 > s2 allowed).
 *   kind 'F': Int ds/Sqrt,  'G': Int s ds/Sqrt (cubic only),
 *   'H': Int ds/((s - q) Sqrt).
 * ---------------------------------------------------------------------- */
typedef struct {
    Expr* pi;        /* leading constant of Phat */
    Expr* A[2];
    Expr* B[2];
} EvenQ;

static Expr* carlson_U(Expr* X[4], Expr* Y[4], int i, int j, Expr* dxy) {
    int k = -1, l = -1;
    for (int m = 0; m < 4; m++) {
        if (m == i || m == j) continue;
        if (k < 0) k = m; else l = m;
    }
    return ell_x("(#1 #2 #3 #4 + #5 #6 #7 #8)/#9", 9,
                 X[i], X[j], Y[k], Y[l], Y[i], Y[j], X[k], X[l], dxy);
}

static Expr* s_integral(const EvenQ* Q, int cubic, char kind,
                        Expr* s1in, Expr* s2in, Expr* q) {
    double d1 = ell_num(s1in), d2 = ell_num(s2in);
    if (isnan(d1) || isnan(d2)) return NULL;
    Expr* s1 = s1in; Expr* s2 = s2in;
    int sign = 1;
    if (d1 > d2) { s1 = s2in; s2 = s1in; sign = -1; double t = d1; d1 = d2; d2 = t; }
    if (d2 - d1 <= 0) return E_int(0);
    double mid = 0.5 * (d1 + d2);
    /* orient the linear factors positive on the range */
    Expr* fa[4]; Expr* fb[4];
    int nf = 0;
    if (cubic) { fa[nf] = E_int(0); fb[nf] = E_int(1); nf++; }
    int sg = 1;
    for (int i = 0; i < 2; i++) {
        double a = ell_num(Q->A[i]), b = ell_num(Q->B[i]);
        if (isnan(a) || isnan(b)) { for (int k = 0; k < nf; k++) { expr_free(fa[k]); expr_free(fb[k]); } return NULL; }
        int o = (a + b * mid > 0) ? 1 : -1;
        sg *= o;
        fa[nf] = o > 0 ? expr_copy(Q->A[i]) : ell_x("-#1", 1, Q->A[i]);
        fb[nf] = o > 0 ? expr_copy(Q->B[i]) : ell_x("-#1", 1, Q->B[i]);
        nf++;
    }
    Expr* Cp = sg > 0 ? expr_copy(Q->pi) : ell_x("-#1", 1, Q->pi);
    Expr* result = NULL;
    double cpd = ell_num(Cp);
    if (isnan(cpd) || cpd <= 0) { goto done; }
    Expr* dxy = E_sub(s2, s1);
    double qd = q ? ell_num(q) : 0.0;
    if (q && isnan(qd)) { expr_free(dxy); goto done; }
    if (cubic) {
        fa[3] = E_int(1); fb[3] = E_int(0); nf = 4;
        Expr* X[4]; Expr* Y[4];
        for (int i = 0; i < 4; i++) {
            X[i] = ell_x("Sqrt[#1 + #2 #3]", 3, fa[i], fb[i], s2);
            Y[i] = ell_x("Sqrt[#1 + #2 #3]", 3, fa[i], fb[i], s1);
        }
        Expr* U12 = carlson_U(X, Y, 0, 1, dxy);
        Expr* U13 = carlson_U(X, Y, 0, 2, dxy);
        Expr* U14 = carlson_U(X, Y, 0, 3, dxy);
        Expr* u12 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, U12);
        Expr* u13 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, U13);
        Expr* u14 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, U14);
        Expr* body = NULL;
        if (kind == 'F') {
            Expr* rf = leg_RF(u12, u13, u14);
            if (rf) { body = ell_x("2 #1", 1, rf); expr_free(rf); }
        } else if (kind == 'G') {
            /* 19.29.7 with delta = 4 (the constant factor):
             *   I(e_al - e_4) = Int (a_al + b_al s) ds/Sqrt
             *     = 2/3 d_al,be d_al,ga R_D(U_al,be^2, U_al,ga^2, U_al,4^2) + 2 X_al Y_al/U_al,4,
             * needing U_al,4 != 0 (it vanishes when both limits are roots of
             * the right pair; at most one U vanishes, DLMF 19.29.5).  alpha = 1
             * (the factor s) gives Int s ds/Sqrt directly; alpha = 2, 3 give
             * a_al I0 + b_al Int s ds/Sqrt.  U_24 = U_13, U_34 = U_12. */
            Expr* Ual4[3] = { U14, U13, U12 };
            Expr* ual4[3] = { u14, u13, u12 };
            Expr* uab[3][2] = { { u12, u13 }, { u12, u14 }, { u13, u14 } };
            int al = -1;
            for (int i = 0; i < 3 && al < 0; i++) if (!ell_zero(ual4[i])) al = i;
            Expr* rd = al >= 0 ? leg_RD(uab[al][0], uab[al][1], ual4[al]) : NULL;
            Expr* rf = (rd && al > 0) ? leg_RF(u12, u13, u14) : NULL;
            if (rd && (al == 0 || rf)) {
                int be = al == 0 ? 1 : 0, ga = al == 2 ? 1 : 2;
                Expr* g = ell_x("2/3 (#1 #4 - #2 #3) (#1 #6 - #5 #3) #7 + 2 #8 #9/#10",
                                10, fa[al], fa[be], fb[al], fb[be], fa[ga], fb[ga], rd,
                                X[al], Y[al], Ual4[al]);
                body = al == 0 ? expr_copy(g)
                               : ell_x("(#1 - 2 #2 #3)/#4", 4, g, fa[al], rf, fb[al]);
                expr_free(g);
            }
            if (rd) expr_free(rd);
            if (rf) expr_free(rf);
        } else {
            /* 19.29.8 with alpha = 4 (constant factor), fifth factor s - q */
            Expr* a5 = ell_x("-#1", 1, q);
            Expr* b5 = E_int(1);
            Expr* d5[3];
            for (int i = 0; i < 3; i++)          /* d_i5 = a_i b5 - a5 b_i */
                d5[i] = ell_x("#1 #2 - #3 #4", 4, fa[i], b5, a5, fb[i]);
            /* U45^2 = U41^2 - d42 d43 d15/d45, d4j = b_j, d45 = b5 */
            Expr* u45 = ell_x("Integrate`Elliptic`Private`simp[#1 - #2 #3 #4/#5]", 5, u14, fb[1], fb[2], d5[0], b5);
            Expr* Q2 = ell_x("Integrate`Elliptic`Private`simp[(#1 + #2 #3) (#1 + #2 #4) #5]", 5, a5, b5, s2, s1, u45);
            Expr* S2 = ell_x("Integrate`Elliptic`Private`simp[#1 + #2 #3 #4/#5]", 5, Q2, d5[0], d5[1], d5[2], b5);
            Expr* sv = ell_x("#1 #2 #3 (#7 + #8 #10) + #4 #5 #6 (#7 + #8 #9)", 10,
                             X[0], X[1], X[2], Y[0], Y[1], Y[2], a5, b5, s2, s1);
            double svd = ell_num(sv);
            expr_free(sv);
            Expr* rj = leg_RJ(u12, u13, u14, u45);
            Expr* rc = leg_RCpv(S2, Q2);
            /* both endpoints at roots: S^2 = 0 and the R_C term vanishes
             * (R_C(0, y) is 0 as a principal value, Pi/(2 Sqrt[y]) otherwise
             * -- but then Q^2 > 0 means the pole is outside, sv > 0) */
            bool s2zero = ell_zero(S2);
            if (rj && s2zero && ell_num(Q2) < 0) {
                body = ell_x("2/3 #1 #2 #3/#4 #5", 5, fb[0], fb[1], fb[2], b5, rj);
            } else if (rj && rc && !isnan(svd) && fabs(svd) > 1e-25) {
                body = svd > 0
                    ? ell_x("2/3 #1 #2 #3/#4 #5 + 2 #6", 6, fb[0], fb[1], fb[2], b5, rj, rc)
                    : ell_x("2/3 #1 #2 #3/#4 #5 - 2 #6", 6, fb[0], fb[1], fb[2], b5, rj, rc);
            } else if (rj && rc) {
                /* sign(S) undetermined (both limits at roots, S = 0 while the
                 * R_C term is Pi/(2 Q) != 0): the two candidates differ by
                 * 2 Pi/(2 Q); a machine quadrature of the s-integral (smooth
                 * after s = m - h Cos[th]) picks the right one, and the final
                 * 30-digit verification still guards the result */
                Expr* bp = ell_x("2/3 #1 #2 #3/#4 #5 + 2 #6", 6, fb[0], fb[1], fb[2], b5, rj, rc);
                Expr* bm = ell_x("2/3 #1 #2 #3/#4 #5 - 2 #6", 6, fb[0], fb[1], fb[2], b5, rj, rc);
                Expr* nq = ell_x("With[{m = (#1 + #2)/2, h = (#2 - #1)/2}, NIntegrate[Evaluate[(1/((Integrate`Elliptic`Private`u - #3) "
                                 "Sqrt[Integrate`Elliptic`Private`u (#4 + #5 Integrate`Elliptic`Private`u) (#6 + #7 Integrate`Elliptic`Private`u)]) "
                                 "/. Integrate`Elliptic`Private`u -> m - h Cos[Integrate`Elliptic`Private`th]) h Sin[Integrate`Elliptic`Private`th]], "
                                 "{Integrate`Elliptic`Private`th, 0, Pi/2, Pi}]]", 7, s1, s2, q, fa[1], fb[1], fa[2], fb[2]);
                double nv = ell_num(nq), vp = ell_num(bp), vm = ell_num(bm);
                expr_free(nq);
                double ep = fabs(vp - nv), em = fabs(vm - nv);
                if (!isnan(nv) && !isnan(vp) && !isnan(vm) && fmin(ep, em) < 1e-6 * (1 + fabs(nv)) &&
                    fmax(ep, em) > 1e-3 * (1 + fabs(nv))) {
                    body = ep < em ? bp : bm;
                    expr_free(ep < em ? bm : bp);
                } else { expr_free(bp); expr_free(bm); }
            }
            if (rj) expr_free(rj);
            if (rc) expr_free(rc);
            for (int i = 0; i < 3; i++) expr_free(d5[i]);
            expr_free(a5); expr_free(b5); expr_free(u45); expr_free(Q2); expr_free(S2);
        }
        for (int i = 0; i < 4; i++) { expr_free(X[i]); expr_free(Y[i]); }
        expr_free(U12); expr_free(U13); expr_free(U14);
        expr_free(u12); expr_free(u13); expr_free(u14);
        if (body) {
            result = sign > 0 ? ell_x("#1/Sqrt[#2]", 2, body, Cp) : ell_x("-#1/Sqrt[#2]", 2, body, Cp);
            expr_free(body);
        }
    } else {
        Expr* X[2]; Expr* Y[2];
        for (int i = 0; i < 2; i++) {
            X[i] = ell_x("Sqrt[#1 + #2 #3]", 3, fa[i], fb[i], s2);
            Y[i] = ell_x("Sqrt[#1 + #2 #3]", 3, fa[i], fb[i], s1);
        }
        Expr* V = ell_x("(#1 #4 + #3 #2)/#5", 5, X[0], X[1], Y[0], Y[1], dxy);
        Expr* v2 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, V);
        Expr* body = NULL;
        if (kind == 'F') {
            Expr* U = ell_x("(#1 #2 + #3 #4)/#5", 5, X[0], X[1], Y[0], Y[1], dxy);
            Expr* u2 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, U);
            Expr* rc = leg_RC(u2, v2);
            if (rc) { body = ell_x("2 #1", 1, rc); expr_free(rc); }
            expr_free(U); expr_free(u2);
        } else if (qd > d1 && qd < d2) {
            /* principal value: -(Log|N(s2)| - Log|N(s1)|)/Sqrt[R(q)],
             * N(s) = (2 R(q) + R'(q) (s - q) + 2 Sqrt[R(q) R(s)])/(s - q) */
            Expr* Rq = ell_x("Integrate`Elliptic`Private`simp[(#1 + #2 #5) (#3 + #4 #5)]", 5, fa[0], fb[0], fa[1], fb[1], q);
            Expr* Rp = ell_x("Integrate`Elliptic`Private`simp[#2 (#3 + #4 #5) + #4 (#1 + #2 #5)]", 5, fa[0], fb[0], fa[1], fb[1], q);
            Expr* N2 = ell_x("(2 #1 + #2 (#3 - #4) + 2 Sqrt[#1] #5 #6)/(#3 - #4)", 6, Rq, Rp, s2, q, X[0], X[1]);
            Expr* N1 = ell_x("(2 #1 + #2 (#3 - #4) + 2 Sqrt[#1] #5 #6)/(#3 - #4)", 6, Rq, Rp, s1, q, Y[0], Y[1]);
            double n2 = ell_num(N2), n1 = ell_num(N1);
            if (!isnan(n1) && !isnan(n2) && n1 != 0 && n2 != 0)
            {
                Expr* g2 = E_int(n2 > 0 ? 1 : -1);
                Expr* g1 = E_int(n1 > 0 ? 1 : -1);
                body = ell_x("-(Log[#1 #2] - Log[#3 #4])/Sqrt[#5]", 5, g2, N2, g1, N1, Rq);
                expr_free(g2); expr_free(g1);
            }
            expr_free(Rq); expr_free(Rp); expr_free(N2); expr_free(N1);
        } else {
            /* 19.29.8 with two constant factors: 2 R_C(Q^2 + d15 d25, Q^2),
             * fifth factor oriented positive on the range */
            int o = qd < d1 ? 1 : -1;            /* s - q > 0 iff q < s1 */
            Expr* a5 = o > 0 ? ell_x("-#1", 1, q) : expr_copy(q);
            Expr* b5 = E_int(o);
            Expr* Q2 = ell_x("Integrate`Elliptic`Private`simp[(#1 + #2 #3) (#1 + #2 #4) #5]", 5, a5, b5, s2, s1, v2);
            Expr* c = ell_x("Integrate`Elliptic`Private`simp[#5 + (#1 #3 - #2 #4) (#6 #3 - #2 #7)]", 7,
                            fa[0], a5, b5, fb[0], Q2, fa[1], fb[1]);
            Expr* rc = leg_RC(c, Q2);
            if (rc) { body = ell_x("2 #1 #2", 2, b5, rc); expr_free(rc); }
            expr_free(a5); expr_free(b5); expr_free(Q2); expr_free(c);
        }
        for (int i = 0; i < 2; i++) { expr_free(X[i]); expr_free(Y[i]); }
        expr_free(V); expr_free(v2);
        if (body) {
            result = sign > 0 ? ell_x("#1/Sqrt[#2]", 2, body, Cp) : ell_x("-#1/Sqrt[#2]", 2, body, Cp);
            expr_free(body);
        }
    }
    expr_free(dxy);
done:
    for (int k = 0; k < nf; k++) { expr_free(fa[k]); expr_free(fb[k]); }
    expr_free(Cp);
    return result;
}

/* -------------------------------------------------------------------------
 * t-space: partial fractions, Hermite reduction, evaluation.
 * ---------------------------------------------------------------------- */
typedef struct {
    Expr* tau;
    Expr* L[ELL_MAXK + 1];   /* L[k]: coefficient of (t - tau)^-k */
} TPole;

typedef struct {
    int kind;                /* 0: coef t^m Sqrt[Phat];  1: coef Sqrt[Phat]/(t-tau)^m */
    int m;
    Expr* tau;
    Expr* coef;
} AlgT;

typedef struct {
    Expr* poly[ELL_MAXP];
    TPole pole[ELL_MAXPOL];
    int npole;
    AlgT alg[ELL_MAXALG];
    int nalg;
    Expr* c[5];              /* Phat = c0 + c2 t^2 + c4 t^4 (c1 = c3 = 0) */
    bool ok;
} TSys;

static void tsys_init(TSys* T) {
    memset(T, 0, sizeof *T);
    T->ok = true;
}
static void tsys_free(TSys* T) {
    for (int i = 0; i < ELL_MAXP; i++) if (T->poly[i]) expr_free(T->poly[i]);
    for (int p = 0; p < T->npole; p++) {
        expr_free(T->pole[p].tau);
        for (int k = 0; k <= ELL_MAXK; k++) if (T->pole[p].L[k]) expr_free(T->pole[p].L[k]);
    }
    for (int i = 0; i < T->nalg; i++) {
        expr_free(T->alg[i].coef);
        if (T->alg[i].tau) expr_free(T->alg[i].tau);
    }
    for (int i = 0; i < 5; i++) if (T->c[i]) expr_free(T->c[i]);
}

static void t_poly_add(TSys* T, int i, Expr* v) {          /* consumes v */
    if (i < 0 || i >= ELL_MAXP) { T->ok = false; expr_free(v); return; }
    acc_add(&T->poly[i], v);
}

/* index of pole tau (created when absent); -1 on overflow */
static int t_pole_index(TSys* T, Expr* tau) {
    for (int p = 0; p < T->npole; p++) {
        if (expr_eq(T->pole[p].tau, tau)) return p;
        Expr* d = E_sub(T->pole[p].tau, tau);
        bool z = ell_zero(d);
        expr_free(d);
        if (z) return p;
    }
    if (T->npole >= ELL_MAXPOL) return -1;
    TPole* P = &T->pole[T->npole];
    memset(P, 0, sizeof *P);
    P->tau = expr_copy(tau);
    return T->npole++;
}
static void t_pole_add(TSys* T, int p, int k, Expr* v) {   /* consumes v */
    if (p < 0 || k < 1 || k > ELL_MAXK) { T->ok = false; expr_free(v); return; }
    acc_add(&T->pole[p].L[k], v);
}

/* c (t - tau)^j added to the polynomial part; c borrowed */
static void t_poly_shift(TSys* T, Expr* tau, int j, Expr* c) {
    for (int i = 0; i <= j; i++) {
        Expr* bj = E_binom(j, i);
        Expr* ex = E_int(j - i);
        Expr* v = ell_x("#1 #2 If[#4 === 0, 1, (-#3)^#4]", 4, c, bj, tau, ex);
        expr_free(bj); expr_free(ex);
        t_poly_add(T, i, v);
    }
}

static void t_alg_add(TSys* T, int kind, int m, Expr* tau, Expr* coef) {   /* consumes coef */
    if (T->nalg >= ELL_MAXALG) { T->ok = false; expr_free(coef); return; }
    AlgT* a = &T->alg[T->nalg++];
    a->kind = kind; a->m = m; a->tau = tau ? expr_copy(tau) : NULL; a->coef = coef;
}

/* Taylor coefficients of Phat at tau */
static void phat_taylor(TSys* T, Expr* tau, Expr* e[5]) {
    e[0] = ell_x("Integrate`Elliptic`Private`simp[#1 + #2 #4^2 + #3 #4^4]", 4, T->c[0], T->c[2], T->c[4], tau);
    e[1] = ell_x("Integrate`Elliptic`Private`simp[2 #1 #3 + 4 #2 #3^3]", 3, T->c[2], T->c[4], tau);
    e[2] = ell_x("Integrate`Elliptic`Private`simp[#1 + 6 #2 #3^2]", 3, T->c[2], T->c[4], tau);
    e[3] = ell_x("Integrate`Elliptic`Private`simp[4 #1 #2]", 2, T->c[4], tau);
    e[4] = expr_copy(T->c[4]);
}

/* Hermite reduction against Phat (see the file header).  A pole at a point
 * where Phat != 0 loses its orders >= 2:
 *   d/dt[Sqrt[Phat]/(t-tau)^(k-1)] Sqrt[Phat]
 *       = sum_j e_j (j/2 - (k-1)) (t-tau)^(j-k);
 * a pole at a simple root of Phat (e_0 = 0) loses every order:
 *   d/dt[Sqrt[Phat]/(t-tau)^k] Sqrt[Phat] = sum_{j>=1} e_j (j/2 - k) (t-tau)^(j-k-1);
 * and t^n, n >= 3, drops through
 *   d/dt[t^m Sqrt[Phat]] Sqrt[Phat] = sum_j c_j (m + j/2) t^(m+j-1). */
static void t_hermite(TSys* T) {
    for (int p = 0; p < T->npole && T->ok; p++) {
        TPole* P = &T->pole[p];
        Expr* e[5];
        phat_taylor(T, P->tau, e);
        bool root = ell_zero(e[0]);
        int K = 0;
        for (int k = ELL_MAXK; k >= 1; k--)
            if (P->L[k] && !ell_zero(P->L[k])) { K = k; break; }
        int lowest = root ? 1 : 2;
        for (int k = K; k >= lowest && T->ok; k--) {
            if (!P->L[k]) continue;
            Expr* kE = E_int(k);
            Expr* coef = root ? ell_x("Together[#1/(#3 (1/2 - #2))]", 3, P->L[k], kE, e[1])
                              : ell_x("Together[#1/(-(#2 - 1) #3)]", 3, P->L[k], kE, e[0]);
            expr_free(P->L[k]); P->L[k] = NULL;
            int jlo = root ? 2 : 1;
            for (int j = jlo; j <= 4; j++) {
                int ex = root ? j - k - 1 : j - k;
                Expr* jE = E_int(j);
                Expr* neg = root ? ell_x("-#1 #2 (#3/2 - #4)", 4, coef, e[j], jE, kE)
                                 : ell_x("-#1 #2 (#3/2 - (#4 - 1))", 4, coef, e[j], jE, kE);
                expr_free(jE);
                if (ex < 0) t_pole_add(T, p, -ex, neg);
                else { t_poly_shift(T, P->tau, ex, neg); expr_free(neg); }
            }
            expr_free(kE);
            t_alg_add(T, 1, root ? k : k - 1, P->tau, coef);
        }
        if (root && P->L[1]) { expr_free(P->L[1]); P->L[1] = NULL; }
        for (int j = 0; j < 5; j++) expr_free(e[j]);
    }
    for (int n = ELL_MAXP - 1; n >= 3 && T->ok; n--) {
        if (!T->poly[n]) continue;
        if (ell_zero(T->poly[n])) { expr_free(T->poly[n]); T->poly[n] = NULL; continue; }
        int m = n - 3;
        Expr* mE = E_int(m);
        Expr* coef = ell_x("Together[#1/(#2 (#3 + 2))]", 3, T->poly[n], T->c[4], mE);
        expr_free(T->poly[n]); T->poly[n] = NULL;
        for (int j = 0; j <= 2; j += 2) {
            int idx = m + j - 1;
            if (idx < 0) continue;
            Expr* jE = E_int(j);
            Expr* v = ell_x("-#1 #2 (#3 + #4/2)", 4, coef, T->c[j], mE, jE);
            expr_free(jE);
            t_poly_add(T, idx, v);
        }
        expr_free(mE);
        t_alg_add(T, 0, m, NULL, coef);
    }
}

/* value of the algebraic part at t (borrowed); NULL when singular */
static Expr* t_alg_value(TSys* T, Expr* t) {
    Expr* sq = ell_x("Sqrt[Integrate`Elliptic`Private`simp[#1 + #2 #4^2 + #3 #4^4]]", 4, T->c[0], T->c[2], T->c[4], t);
    Expr* acc = NULL;
    for (int i = 0; i < T->nalg; i++) {
        AlgT* a = &T->alg[i];
        Expr* mE = E_int(a->m);
        Expr* v;
        if (a->kind == 0) v = ell_x("#1 If[#3 === 0, 1, #2^#3] #4", 4, a->coef, t, mE, sq);
        else {
            Expr* d = E_sub(t, a->tau);
            if (ell_zero(d)) {
                expr_free(d); expr_free(mE); expr_free(sq);
                if (acc) expr_free(acc);
                return NULL;
            }
            v = ell_x("#1 #4/#2^#3", 4, a->coef, d, mE, sq);
            expr_free(d);
        }
        expr_free(mE);
        acc_add(&acc, v);
    }
    expr_free(sq);
    return acc ? acc : E_int(0);
}

/* Basis integrals over a one-signed t segment [u1, u2] (oriented):
 * with s = t^2 (dt = +-ds/(2 Sqrt[s]) for t >< 0, t dt = ds/2),
 *   T0 -> +-1/2 Int ds/Sqrt[s L1 L2]          T1 -> 1/2 Int ds/Sqrt[L1 L2]
 *   T2 -> +-1/2 Int s ds/Sqrt[s L1 L2]
 *   J_tau = Int (t + tau) dt/((t^2 - tau^2) Sqrt)
 *         -> +-tau/2 Int ds/((s - tau^2) Sqrt[s L1 L2]) + 1/2 Int ds/((s - tau^2) Sqrt[L1 L2]). */
static Expr* t_segment(TSys* T, const EvenQ* Q, Expr* u1, Expr* u2) {
    double a = ell_num(u1), b = ell_num(u2);
    if (isnan(a) || isnan(b)) return NULL;
    if (a == b) return E_int(0);
    bool neg = (a + b) < 0;
    Expr* s1 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, u1);
    Expr* s2 = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, u2);
    double s1d = a * a, s2d = b * b;
    Expr* es = E_int(neg ? -1 : 1);
    Expr* acc = NULL;
    bool ok = true;
    if (T->poly[0] && !ell_zero(T->poly[0])) {
        Expr* v = s_integral(Q, 1, 'F', s1, s2, NULL);
        if (!v) (ell_dbg("T0 failed"), ok = false);
        else { acc_add(&acc, ell_x("#1 #2 #3/2", 3, es, T->poly[0], v)); expr_free(v); }
    }
    if (ok && T->poly[1] && !ell_zero(T->poly[1])) {
        Expr* v = s_integral(Q, 0, 'F', s1, s2, NULL);
        if (!v) (ell_dbg("T1 failed"), ok = false);
        else { acc_add(&acc, ell_x("#1 #2/2", 2, T->poly[1], v)); expr_free(v); }
    }
    if (ok && T->poly[2] && !ell_zero(T->poly[2])) {
        Expr* v = s_integral(Q, 1, 'G', s1, s2, NULL);
        if (!v) (ell_dbg("T2 failed"), ok = false);
        else { acc_add(&acc, ell_x("#1 #2 #3/2", 3, es, T->poly[2], v)); expr_free(v); }
    }
    for (int p = 0; ok && p < T->npole; p++) {
        TPole* P = &T->pole[p];
        if (!P->L[1] || ell_zero(P->L[1])) continue;
        Expr* q = ell_x("Integrate`Elliptic`Private`simp[#1^2]", 1, P->tau);
        double qd = ell_num(q);
        /* the mirror pole on a segment endpoint: both parts diverge */
        if (isnan(qd) || fabs(qd - s1d) < 1e-12 * (1 + qd) || fabs(qd - s2d) < 1e-12 * (1 + qd)) {
            expr_free(q); (ell_dbg("mirror pole on endpoint"), ok = false); break;
        }
        if (!ell_zero(P->tau)) {
            Expr* v = s_integral(Q, 1, 'H', s1, s2, q);
            if (!v) (ell_dbg("J even failed"), ok = false);
            else { acc_add(&acc, ell_x("#1 #2 #3 #4/2", 4, es, P->L[1], P->tau, v)); expr_free(v); }
        }
        if (ok) {
            Expr* v = s_integral(Q, 0, 'H', s1, s2, q);
            if (!v) (ell_dbg("J odd failed"), ok = false);
            else { acc_add(&acc, ell_x("#1 #2/2", 2, P->L[1], v)); expr_free(v); }
        }
        expr_free(q);
    }
    expr_free(s1); expr_free(s2); expr_free(es);
    if (!ok) { if (acc) expr_free(acc); return NULL; }
    return acc ? acc : E_int(0);
}

/* -------------------------------------------------------------------------
 * The integrand in x: W = sum wpoly[n] x^n + sum c/(x - p)^k, and the
 * factorisation P = pi Q1 Q2.
 * ---------------------------------------------------------------------- */
typedef struct {
    Expr* wpoly[ELL_MAXP];
    Expr* pp[ELL_MAXPOL];
    Expr* pc[ELL_MAXPOL][ELL_MAXK + 1];
    double pd[ELL_MAXPOL];
    int np;
    Expr* Q[2][3];          /* Q_i = Q[i][0] x^2 + Q[i][1] x + Q[i][2] */
    Expr* pi;               /* P = pi Q1 Q2 */
} XSys;

static void xsys_free(XSys* S) {
    for (int i = 0; i < ELL_MAXP; i++) if (S->wpoly[i]) expr_free(S->wpoly[i]);
    for (int p = 0; p < S->np; p++) {
        expr_free(S->pp[p]);
        for (int k = 0; k <= ELL_MAXK; k++) if (S->pc[p][k]) expr_free(S->pc[p][k]);
    }
    for (int i = 0; i < 2; i++) for (int j = 0; j < 3; j++) if (S->Q[i][j]) expr_free(S->Q[i][j]);
    if (S->pi) expr_free(S->pi);
}

static bool x_pole_add(XSys* S, Expr* p, int k, Expr* c) {    /* c consumed */
    if (k < 1 || k > ELL_MAXK) { expr_free(c); return false; }
    int idx = -1;
    for (int i = 0; i < S->np; i++) {
        Expr* d = E_sub(S->pp[i], p);
        bool z = ell_zero(d);
        expr_free(d);
        if (z) { idx = i; break; }
    }
    if (idx < 0) {
        if (S->np >= ELL_MAXPOL) { expr_free(c); return false; }
        idx = S->np++;
        S->pp[idx] = expr_copy(p);
        S->pd[idx] = ell_num(p);
        for (int j = 0; j <= ELL_MAXK; j++) S->pc[idx][j] = NULL;
        if (isnan(S->pd[idx])) { expr_free(c); return false; }  /* complex pole */
    }
    acc_add(&S->pc[idx][k], c);
    return true;
}

/* Partial fractions of W over the reals (linear and real-rooted quadratic
 * denominators). */
static bool x_partial(XSys* S, Expr* W, Expr* x) {
    Expr* ap = ell_x("With[{ap = Apart[#1, #2]}, If[Head[ap] === Plus, List @@ ap, {ap}]]", 2, W, x);
    bool ok = ell_len(ap) > 0;
    for (size_t i = 0; ok && i < ell_len(ap); i++) {
        Expr* T = ell_arg(ap, i);
        Expr* den = ell_x("Denominator[#1]", 1, T);
        if (!contains_sym(den, x)) {
            Expr* cl = ell_x("CoefficientList[#1, #2]", 2, T, x);
            size_t n = ell_len(cl);
            if (n == 0 || n > ELL_MAXP) ok = false;
            for (size_t j = 0; ok && j < n; j++)
                acc_add(&S->wpoly[j], expr_copy(ell_arg(cl, j)));
            expr_free(cl);
            expr_free(den);
            continue;
        }
        /* the single x-dependent factor of the denominator */
        Expr* fb = ell_x("Select[FactorList[#1], Function[{fq}, !FreeQ[fq, #2]]]", 2, den, x);
        expr_free(den);
        if (ell_len(fb) != 1 || ell_len(ell_arg(fb, 0)) != 2) { expr_free(fb); ok = false; break; }
        Expr* base = ell_arg(ell_arg(fb, 0), 0);
        Expr* kx = ell_arg(ell_arg(fb, 0), 1);
        int k = kx->type == EXPR_INTEGER ? (int)kx->data.integer : 0;
        Expr* cl = ell_x("CoefficientList[#1, #2]", 2, base, x);
        size_t deg = ell_len(cl) ? ell_len(cl) - 1 : 0;
        if (deg == 1 && k >= 1) {
            Expr* p = ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, ell_arg(cl, 0), ell_arg(cl, 1));
            Expr* kE = E_int(k);
            Expr* c = ell_x("Integrate`Elliptic`Private`simp[Together[#1 (#2 - #3)^#4]]", 4, T, x, p, kE);
            expr_free(kE);
            if (contains_sym(c, x)) { expr_free(c); ok = false; }
            else ok = x_pole_add(S, p, k, c);
            expr_free(p);
        } else if (deg == 2 && k == 1) {
            Expr* disc = ell_x("#2^2 - 4 #1 #3", 3, ell_arg(cl, 0), ell_arg(cl, 1), ell_arg(cl, 2));
            double dd = ell_num(disc);
            if (isnan(dd) || dd <= 0) ok = false;          /* complex poles */
            for (int sgn = -1; ok && sgn <= 1; sgn += 2) {
                Expr* sE = E_int(sgn);
                Expr* r = ell_x("Integrate`Elliptic`Private`simp[(-#1 + #3 Sqrt[#2])/(2 #4)]", 4,
                                ell_arg(cl, 1), disc, sE, ell_arg(cl, 2));
                Expr* c = ell_x("Integrate`Elliptic`Private`simp[(Numerator[#1]/D[Denominator[#1], #2]) /. #2 -> #3]",
                                3, T, x, r);
                expr_free(sE);
                if (contains_sym(c, x)) { expr_free(c); ok = false; }
                else ok = x_pole_add(S, r, 1, c);
                expr_free(r);
            }
            expr_free(disc);
        } else ok = false;
        expr_free(cl);
        expr_free(fb);
    }
    expr_free(ap);
    return ok;
}

typedef struct { Expr* r; double d; int parent; } RRoot;

/* the monic quadratic form with roots r1, r2 (the rational parent when both
 * come from the same rational quadratic factor) */
static Expr* form_pair(RRoot* r1, RRoot* r2, Expr* par[][3]) {
    if (r1->parent >= 0 && r1->parent == r2->parent)
        return ell_x("{#1, #2, #3}", 3, par[r1->parent][0], par[r1->parent][1], par[r1->parent][2]);
    return ell_x("{1, Integrate`Elliptic`Private`simp[-(#1 + #2)], Integrate`Elliptic`Private`simp[#1 #2]}", 2, r1->r, r2->r);
}

/* Real factorisation P = pi Q1 Q2 into two binary quadratic forms whose root
 * pairs do not interlace on the projective line (a cubic's fourth root is
 * Infinity, giving a linear form). */
static bool x_pair(XSys* S, Expr* P, Expr* x, int pick, int* nopt) {
    Expr* fl = ell_x("FactorList[#1]", 1, P);
    RRoot rr[4]; int nr = 0;
    Expr* cq[2][3]; int nc = 0;
    Expr* par[4][3]; int npar = 0;
    bool ok = ell_len(fl) > 0;
    for (size_t i = 0; ok && i < ell_len(fl); i++) {
        Expr* f = ell_arg(ell_arg(fl, i), 0);
        Expr* m = ell_arg(ell_arg(fl, i), 1);
        if (!contains_sym(f, x)) continue;
        if (!(m->type == EXPR_INTEGER && m->data.integer == 1)) { ok = false; break; }
        Expr* cl = ell_x("CoefficientList[#1, #2]", 2, f, x);
        size_t deg = ell_len(cl) ? ell_len(cl) - 1 : 0;
        if (deg == 1 && nr < 4) {
            rr[nr].r = ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, ell_arg(cl, 0), ell_arg(cl, 1));
            rr[nr].d = ell_num(rr[nr].r); rr[nr].parent = -1; nr++;
        } else if (deg == 2) {
            Expr *c0 = ell_arg(cl, 0), *c1 = ell_arg(cl, 1), *c2 = ell_arg(cl, 2);
            Expr* disc = ell_x("#2^2 - 4 #1 #3", 3, c0, c1, c2);
            double dd = ell_num(disc);
            if (isnan(dd)) ok = false;
            else if (dd < 0) {
                if (nc >= 2) ok = false;
                else {
                    cq[nc][0] = E_int(1);
                    cq[nc][1] = ell_x("#1/#2", 2, c1, c2);
                    cq[nc][2] = ell_x("#1/#2", 2, c0, c2);
                    nc++;
                }
            } else if (nr <= 2 && npar < 4) {
                par[npar][0] = E_int(1);
                par[npar][1] = ell_x("#1/#2", 2, c1, c2);
                par[npar][2] = ell_x("#1/#2", 2, c0, c2);
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    Expr* sE = E_int(sgn);
                    rr[nr].r = ell_x("Integrate`Elliptic`Private`simp[(-#1 + #3 Sqrt[#2])/(2 #4)]", 4, c1, disc, sE, c2);
                    expr_free(sE);
                    rr[nr].d = ell_num(rr[nr].r); rr[nr].parent = npar; nr++;
                }
                npar++;
            } else ok = false;
            expr_free(disc);
        } else if ((deg == 3 || deg == 4) && nr == 0 && nc == 0) {
            /* an irreducible cubic or quartic factor: its roots in radicals
             * (decline on Root objects -- casus irreducibilis, or a quartic
             * Solve leaves implicit), real ones kept, complex ones paired
             * with their conjugates into real quadratic forms
             * x^2 - 2 Re[r] x + Abs[r]^2. */
            Expr* rs = ell_x("With[{rl = #2 /. Solve[#1 == 0, #2]}, If[FreeQ[rl, Root] && "
                             "Length[rl] == Exponent[#1, #2], "
                             "Map[Function[{u}, {u, N[N[Re[u], 40]], N[N[Im[u], 40]]}], rl], {}]]", 2, f, x);
            size_t nrs = ell_len(rs);
            if (nrs != deg) ok = false;
            bool used[4] = { false, false, false, false };
            for (size_t i = 0; ok && i < nrs; i++) {
                if (used[i]) continue;
                Expr* ri = ell_arg(rs, i);
                Expr* imE = ell_arg(ri, 2);
                double im = imE->type == EXPR_REAL ? imE->data.real : (imE->type == EXPR_INTEGER ? (double)imE->data.integer : NAN);
                Expr* reE = ell_arg(ri, 1);
                double re = reE->type == EXPR_REAL ? reE->data.real : (reE->type == EXPR_INTEGER ? (double)reE->data.integer : NAN);
                if (isnan(im) || isnan(re)) { ok = false; break; }
                used[i] = true;
                if (fabs(im) < 1e-20 * (1 + fabs(re))) {
                    if (nr >= 4) { ok = false; break; }
                    rr[nr].r = ell_x("Integrate`Elliptic`Private`simp[ComplexExpand[Re[#1]]]", 1, ell_arg(ri, 0));
                    rr[nr].d = re; rr[nr].parent = -1;
                    if (!ell_true(ell_x("FreeQ[#1, Complex]", 1, rr[nr].r))) { expr_free(rr[nr].r); ok = false; break; }
                    nr++;
                } else {
                    /* the conjugate partner */
                    size_t j;
                    for (j = i + 1; j < nrs; j++) {
                        if (used[j]) continue;
                        Expr* rj = ell_arg(rs, j);
                        Expr* a2 = ell_arg(rj, 1); Expr* b2 = ell_arg(rj, 2);
                        if (a2->type == EXPR_REAL && b2->type == EXPR_REAL &&
                            fabs(a2->data.real - re) < 1e-12 * (1 + fabs(re)) &&
                            fabs(b2->data.real + im) < 1e-12 * (1 + fabs(im))) break;
                    }
                    if (j >= nrs || nc >= 2) { ok = false; break; }
                    used[j] = true;
                    cq[nc][0] = E_int(1);
                    cq[nc][1] = ell_x("Integrate`Elliptic`Private`simp[-2 ComplexExpand[Re[#1]]]", 1, ell_arg(ri, 0));
                    cq[nc][2] = ell_x("Integrate`Elliptic`Private`simp[ComplexExpand[Re[#1]^2 + Im[#1]^2]]", 1, ell_arg(ri, 0));
                    nc++;
                    if (!ell_true(ell_x("FreeQ[{#1, #2}, Complex]", 2, cq[nc - 1][1], cq[nc - 1][2]))) { ok = false; break; }
                }
            }
            expr_free(rs);
        } else ok = false;                      /* irreducible quartic, casus irreducibilis */
        expr_free(cl);
    }
    expr_free(fl);
    for (int i = 0; i < nr; i++) if (isnan(rr[i].d)) ok = false;
    int deg = nr + 2 * nc;
    if (ok && (deg < 3 || deg > 4)) ok = false;
    Expr* forms[2] = { NULL, NULL };
    *nopt = (nr >= 3) ? 2 : 1;
    if (ok) {
        for (int i = 0; i < nr; i++)
            for (int j = i + 1; j < nr; j++)
                if (rr[j].d < rr[i].d) { RRoot t = rr[i]; rr[i] = rr[j]; rr[j] = t; }
        int nf = 0;
        for (int i = 0; i < nc; i++)
            forms[nf++] = ell_x("{#1, #2, #3}", 3, cq[i][0], cq[i][1], cq[i][2]);
        if (nr == 2) forms[nf++] = form_pair(&rr[0], &rr[1], par);
        else if (nr == 1) forms[nf++] = ell_x("{0, 1, -#1}", 1, rr[0].r);
        else if (nr == 4 || nr == 3) {
            /* two non-interlacing pairings (the caller tries both and keeps
             * the simpler closed form) */
            Expr* opt[2][2];
            opt[0][0] = form_pair(&rr[0], &rr[1], par);
            opt[1][0] = form_pair(&rr[1], &rr[2], par);
            if (nr == 4) {
                opt[0][1] = form_pair(&rr[2], &rr[3], par);
                opt[1][1] = form_pair(&rr[3], &rr[0], par);
            } else {
                opt[0][1] = ell_x("{0, 1, -#1}", 1, rr[2].r);
                opt[1][1] = ell_x("{0, 1, -#1}", 1, rr[0].r);
            }
            if (pick < 0 || pick > 1) pick = 0;
            forms[0] = expr_copy(opt[pick][0]);
            forms[1] = expr_copy(opt[pick][1]);
            for (int o = 0; o < 2; o++) { expr_free(opt[o][0]); expr_free(opt[o][1]); }
        }
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 3; j++) S->Q[i][j] = expr_copy(ell_arg(forms[i], (size_t)j));
        S->pi = ell_x("Integrate`Elliptic`Private`simp[Together[#1/((#2 #8^2 + #3 #8 + #4) (#5 #8^2 + #6 #8 + #7))]]", 8,
                      P, S->Q[0][0], S->Q[0][1], S->Q[0][2], S->Q[1][0], S->Q[1][1], S->Q[1][2], x);
        if (contains_sym(S->pi, x)) ok = false;
    }
    for (int i = 0; i < 2; i++) if (forms[i]) expr_free(forms[i]);
    for (int i = 0; i < nr; i++) expr_free(rr[i].r);
    for (int i = 0; i < nc; i++) for (int j = 0; j < 3; j++) expr_free(cq[i][j]);
    for (int i = 0; i < npar; i++) for (int j = 0; j < 3; j++) expr_free(par[i][j]);
    return ok;
}

/* -------------------------------------------------------------------------
 * One piece of the interval, mapped by x = (al t + be)/(ga t + de) with
 * x(0) = m0 and x(Infinity) = minf; minf is not in the closed piece.
 * ---------------------------------------------------------------------- */
static Expr* x_piece(XSys* S, Pt* lo, Pt* hi, Pt* m0, Pt* minf) {
    Expr *al, *be, *ga, *de;
    if (m0->inf)        { al = expr_copy(minf->v); be = E_int(1); ga = E_int(1); de = E_int(0); }
    else if (minf->inf) { al = E_int(1); be = expr_copy(m0->v); ga = E_int(0); de = E_int(1); }
    else                { al = expr_copy(minf->v); be = expr_copy(m0->v); ga = E_int(1); de = E_int(1); }
    bool ga0 = ell_zero(ga);
    Expr* Dl = ell_x("Integrate`Elliptic`Private`simp[#1 #4 - #2 #3]", 4, al, be, ga, de);
    Expr* result = NULL;
    TSys T; tsys_init(&T);
    EvenQ Q; memset(&Q, 0, sizeof Q);
    Expr* tA = NULL; Expr* tB = NULL;
    bool ok = true;

    /* the even quartic */
    for (int i = 0; ok && i < 2; i++) {
        Expr* c1 = ell_x("Integrate`Elliptic`Private`simp[2 #1 #4 #5 + #2 (#4 #7 + #5 #6) + 2 #3 #6 #7]", 7,
                         S->Q[i][0], S->Q[i][1], S->Q[i][2], al, be, ga, de);
        if (!ell_zero(c1)) (ell_dbg("not even"), ok = false);
        expr_free(c1);
        Q.A[i] = ell_x("Integrate`Elliptic`Private`simp[#1 #5^2 + #2 #5 #7 + #3 #7^2]", 7,
                       S->Q[i][0], S->Q[i][1], S->Q[i][2], al, be, ga, de);
        Q.B[i] = ell_x("Integrate`Elliptic`Private`simp[#1 #4^2 + #2 #4 #6 + #3 #6^2]", 6,
                       S->Q[i][0], S->Q[i][1], S->Q[i][2], al, (Expr*)be, ga);
    }
    Q.pi = expr_copy(S->pi);
    if (!ok) goto out;
    T.c[0] = ell_x("Integrate`Elliptic`Private`simp[#1 #2 #3]", 3, Q.pi, Q.A[0], Q.A[1]);
    T.c[1] = E_int(0);
    T.c[2] = ell_x("Integrate`Elliptic`Private`simp[#1 (#2 #5 + #3 #4)]", 5, Q.pi, Q.A[0], Q.A[1], Q.B[0], Q.B[1]);
    T.c[3] = E_int(0);
    T.c[4] = ell_x("Integrate`Elliptic`Private`simp[#1 #2 #3]", 3, Q.pi, Q.B[0], Q.B[1]);
    if (ell_zero(T.c[4]) || ell_zero(T.c[0])) { (ell_dbg("Phat degenerate"), ok = false); goto out; }

    /* t endpoints */
    {
        Pt* ends[2] = { lo, hi };
        Expr** tt[2] = { &tA, &tB };
        for (int i = 0; i < 2; i++) {
            if (ends[i]->inf) {
                if (ga0) { ok = false; break; }
                *tt[i] = ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, de, ga);
            } else
                *tt[i] = ell_x("Integrate`Elliptic`Private`simp[(#2 - #3 #5)/(#4 #5 - #1)]", 5, al, be, de, ga, ends[i]->v);
        }
        if (!ok) goto out;
    }

    /* W -> t partial fractions (times the Jacobian Dl) */
    for (int n = 0; ok && n < ELL_MAXP; n++) {
        if (!S->wpoly[n] || ell_zero(S->wpoly[n])) continue;
        Expr* nE = E_int(n);
        if (ga0) {
            Expr* cl = ell_x("CoefficientList[Expand[#1 #2 ((#3 Integrate`Elliptic`Private`t + #4)/#5)^#6], Integrate`Elliptic`Private`t]",
                             6, S->wpoly[n], Dl, al, be, de, nE);
            for (size_t i = 0; i < ell_len(cl); i++) t_poly_add(&T, (int)i, expr_copy(ell_arg(cl, i)));
            expr_free(cl);
        } else {
            Expr* tinf = ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, de, ga);
            int pidx = -1;
            for (int j = 0; j <= n; j++) {
                Expr* jE = E_int(j);
                Expr* bj = E_binom(n, j);
                Expr* v = ell_x("#1 #2 #3 If[#6 === #7, 1, (#4/#5)^(#6 - #7)] (-#2/#5^2)^#7", 7,
                                S->wpoly[n], Dl, bj, al, ga, nE, jE);
                expr_free(jE); expr_free(bj);
                if (j == 0) t_poly_add(&T, 0, v);
                else {
                    if (pidx < 0) pidx = t_pole_index(&T, tinf);
                    t_pole_add(&T, pidx, j, v);
                }
            }
            expr_free(tinf);
        }
        expr_free(nE);
    }
    for (int p = 0; ok && p < S->np; p++) {
        Expr* g = ell_x("Integrate`Elliptic`Private`simp[#1 - #2 #3]", 3, al, S->pp[p], ga);
        Expr* h = ell_x("Integrate`Elliptic`Private`simp[#1 - #2 #3]", 3, be, S->pp[p], de);
        bool g0 = ell_zero(g);
        Expr* tau = g0 ? NULL : ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, h, g);
        Expr* e = g0 ? NULL : ell_x("Integrate`Elliptic`Private`simp[#1 #2 + #3]", 3, ga, tau, de);
        int pidx = g0 ? -1 : t_pole_index(&T, tau);
        for (int k = 1; ok && k <= ELL_MAXK; k++) {
            Expr* c = S->pc[p][k];
            if (!c) continue;
            Expr* kE = E_int(k);
            if (g0) {        /* the pole sits at x(Infinity): polynomial in t */
                Expr* cl = ell_x("CoefficientList[Expand[#1 #2 (#3 Integrate`Elliptic`Private`t + #4)^#5/#6^#5], Integrate`Elliptic`Private`t]",
                                 6, c, Dl, ga, de, kE, h);
                for (size_t i = 0; i < ell_len(cl); i++) t_poly_add(&T, (int)i, expr_copy(ell_arg(cl, i)));
                expr_free(cl);
            } else {
                for (int j = 0; j <= k; j++) {
                    Expr* jE = E_int(j);
                    Expr* bj = E_binom(k, j);
                    Expr* v = ell_x("#1 #2 #3 If[#5 === 0, 1, #4^#5] If[#7 === #5, 1, #6^(#7 - #5)]/#8^#7", 8,
                                    c, Dl, bj, ga, jE, e, kE, g);
                    expr_free(jE); expr_free(bj);
                    if (j == k) t_poly_add(&T, 0, v);
                    else t_pole_add(&T, pidx, k - j, v);
                }
            }
            expr_free(kE);
        }
        expr_free(g); expr_free(h);
        if (tau) expr_free(tau);
        if (e) expr_free(e);
    }
    if (!T.ok) { (ell_dbg("t-transform overflow"), ok = false); goto out; }
    /* light canonicalisation of the coefficients */
    for (int i = 0; i < ELL_MAXP; i++)
        if (T.poly[i]) { Expr* s = ell_x("Integrate`Elliptic`Private`simp[#1]", 1, T.poly[i]); expr_free(T.poly[i]); T.poly[i] = s; }
    for (int p = 0; p < T.npole; p++)
        for (int k = 1; k <= ELL_MAXK; k++)
            if (T.pole[p].L[k]) { Expr* s = ell_x("Integrate`Elliptic`Private`simp[#1]", 1, T.pole[p].L[k]); expr_free(T.pole[p].L[k]); T.pole[p].L[k] = s; }

    t_hermite(&T);
    if (!T.ok) { (ell_dbg("hermite overflow"), ok = false); goto out; }
    for (int i = 0; i < 3; i++)
        if (T.poly[i]) { Expr* s = ell_x("Integrate`Elliptic`Private`simp[#1]", 1, T.poly[i]); expr_free(T.poly[i]); T.poly[i] = s; }
    for (int p = 0; p < T.npole; p++)
        if (T.pole[p].L[1]) { Expr* s = ell_x("Integrate`Elliptic`Private`simp[#1]", 1, T.pole[p].L[1]); expr_free(T.pole[p].L[1]); T.pole[p].L[1] = s; }

    {
        Expr* vb = t_alg_value(&T, tB);
        Expr* va = vb ? t_alg_value(&T, tA) : NULL;
        if (!va) { if (vb) expr_free(vb); (ell_dbg("alg singular"), ok = false); goto out; }
        result = E_sub(vb, va);
        expr_free(vb); expr_free(va);
        double a = ell_num(tA), b = ell_num(tB);
        if (isnan(a) || isnan(b)) { ok = false; goto out; }
        Expr* zero = E_int(0);
        Expr* seg;
        if (a * b < 0) {
            Expr* s1 = t_segment(&T, &Q, tA, zero);
            Expr* s2 = s1 ? t_segment(&T, &Q, zero, tB) : NULL;
            seg = (s1 && s2) ? E_add(s1, s2) : NULL;
            if (s1) expr_free(s1);
            if (s2) expr_free(s2);
        } else seg = t_segment(&T, &Q, tA, tB);
        expr_free(zero);
        if (!seg) { (ell_dbg("segment failed"), ok = false); goto out; }
        acc_add(&result, seg);
    }
out:
    if (!ok && result) { expr_free(result); result = NULL; }
    tsys_free(&T);
    for (int i = 0; i < 2; i++) { if (Q.A[i]) expr_free(Q.A[i]); if (Q.B[i]) expr_free(Q.B[i]); }
    if (Q.pi) expr_free(Q.pi);
    if (tA) expr_free(tA);
    if (tB) expr_free(tB);
    expr_free(al); expr_free(be); expr_free(ga); expr_free(de); expr_free(Dl);
    return result;
}

static bool pt_inside(Pt* p, Pt* lo, Pt* hi) {
    if (p->inf) return lo->inf || hi->inf;
    double tol = 1e-12 * (1 + fabs(p->d));
    return (lo->inf || p->d >= lo->d - tol) && (hi->inf || p->d <= hi->d + tol);
}

/* A simple rational strictly inside (lo, hi) between u and v (both inside
 * the closed piece, at least one finite). */
static Expr* split_point(Pt* lo, Pt* hi, Pt* u, Pt* v) {
    if (!u->inf && !v->inf) {
        /* not the midpoint: that is the mirror image of x = Infinity under
         * the involution t -> -t, where the even/odd split of a pole at
         * Infinity would put a principal-value pole on the segment end */
        double m = u->d + 0.4 * (v->d - u->d), w = fabs(u->d - v->d) / 1000;
        Expr* mE = expr_new_real(m);
        Expr* wE = expr_new_real(w);
        Expr* r = ell_x("Rationalize[#1, #2]", 2, mE, wE);
        expr_free(mE); expr_free(wE);
        return r;
    }
    Pt* f = u->inf ? v : u;
    if (lo->inf) return E_int((long)floor(f->d) - 1);
    if (hi->inf) return E_int((long)ceil(f->d) + 1);
    return NULL;
}

static Expr* x_interval(XSys* S, Pt* lo, Pt* hi, Pt* mu, Pt* nu, int depth) {
    if (depth > 4) return NULL;
    if (!pt_inside(nu, lo, hi)) return x_piece(S, lo, hi, mu, nu);
    if (!pt_inside(mu, lo, hi)) return x_piece(S, lo, hi, nu, mu);
    Expr* c = split_point(lo, hi, mu, nu);
    if (!c) return NULL;
    Pt pc;
    if (!pt_from(c, &pc)) { expr_free(c); return NULL; }
    expr_free(c);
    Expr* r1 = x_interval(S, lo, &pc, mu, nu, depth + 1);
    Expr* r2 = r1 ? x_interval(S, &pc, hi, mu, nu, depth + 1) : NULL;
    pt_free(&pc);
    Expr* r = (r1 && r2) ? E_add(r1, r2) : NULL;
    if (r1) expr_free(r1);
    if (r2) expr_free(r2);
    return r;
}

/* Int_lo^hi W/Sqrt[P] dx for the (prepared) XSys. */
static Expr* x_integral(XSys* S, Pt* lo, Pt* hi) {
    /* the common harmonic pair: roots of the Jacobian of (Q1, Q2) */
    Expr* J2 = ell_x("Integrate`Elliptic`Private`simp[#1 #5 - #4 #2]", 6, S->Q[0][0], S->Q[0][1], S->Q[0][2],
                     S->Q[1][0], S->Q[1][1], S->Q[1][2]);
    Expr* J1 = ell_x("Integrate`Elliptic`Private`simp[2 (#1 #6 - #4 #3)]", 6, S->Q[0][0], S->Q[0][1], S->Q[0][2],
                     S->Q[1][0], S->Q[1][1], S->Q[1][2]);
    Expr* J0 = ell_x("Integrate`Elliptic`Private`simp[#2 #6 - #5 #3]", 6, S->Q[0][0], S->Q[0][1], S->Q[0][2],
                     S->Q[1][0], S->Q[1][1], S->Q[1][2]);
    Pt mu = { NULL, 0, 0 }, nu = { NULL, 0, 0 };
    Expr* r = NULL;
    bool ok = true;
    if (ell_zero(J2)) {
        if (ell_zero(J1)) ok = false;
        else {
            Expr* m = ell_x("Integrate`Elliptic`Private`simp[-#1/#2]", 2, J0, J1);
            ok = pt_from(m, &mu);
            expr_free(m);
            nu.inf = 1; nu.d = INFINITY;
        }
    } else {
        Expr* disc = ell_x("Integrate`Elliptic`Private`simp[#1^2 - 4 #2 #3]", 3, J1, J2, J0);
        double dd = ell_num(disc);
        if (isnan(dd) || dd <= 0) ok = false;
        else {
            Expr* m1 = ell_x("Integrate`Elliptic`Private`simp[(-#1 - Sqrt[#2])/(2 #3)]", 3, J1, disc, J2);
            Expr* m2 = ell_x("Integrate`Elliptic`Private`simp[(-#1 + Sqrt[#2])/(2 #3)]", 3, J1, disc, J2);
            ok = pt_from(m1, &mu) && pt_from(m2, &nu);
            expr_free(m1); expr_free(m2);
        }
        expr_free(disc);
    }
    if (ok) {
        if (lo->inf && hi->inf) {
            /* the whole line: split at a point that is neither mu nor nu */
            Expr* c = E_int(0);
            for (int k = 0; k < 6; k++) {
                Pt t; pt_from(c, &t);
                bool bad = (!mu.inf && fabs(mu.d - t.d) < 1e-9) || (!nu.inf && fabs(nu.d - t.d) < 1e-9);
                pt_free(&t);
                if (!bad) break;
                expr_free(c); c = E_int(k + 1);
            }
            Pt pc; pt_from(c, &pc);
            expr_free(c);
            Expr* r1 = x_interval(S, lo, &pc, &mu, &nu, 0);
            Expr* r2 = r1 ? x_interval(S, &pc, hi, &mu, &nu, 0) : NULL;
            if (r1 && r2) r = E_add(r1, r2);
            if (r1) expr_free(r1);
            if (r2) expr_free(r2);
            pt_free(&pc);
        } else r = x_interval(S, lo, hi, &mu, &nu, 0);
    }
    pt_free(&mu); pt_free(&nu);
    expr_free(J2); expr_free(J1); expr_free(J0);
    return r;
}

/* -------------------------------------------------------------------------
 * Normalisation of the integrand to A + W/Sqrt[P].
 * ---------------------------------------------------------------------- */
#define ELL_MAXB 8

static bool is_half_power(const Expr* e, const Expr* x, int* k) {
    if (!ell_head_is(e, "Power") || e->data.function.arg_count != 2) return false;
    const Expr* ex = e->data.function.args[1];
    if (!ell_head_is(ex, "Rational") || ex->data.function.arg_count != 2) return false;
    const Expr* n = ex->data.function.args[0];
    const Expr* d = ex->data.function.args[1];
    if (n->type != EXPR_INTEGER || d->type != EXPR_INTEGER || d->data.integer != 2) return false;
    if (!contains_sym(e->data.function.args[0], x)) return false;
    *k = (int)n->data.integer;
    return true;
}

static bool collect_bases(const Expr* e, const Expr* x, Expr** B, int* nb) {
    int k;
    if (is_half_power(e, x, &k)) {
        Expr* g = e->data.function.args[0];
        for (int i = 0; i < *nb; i++) if (expr_eq(B[i], g)) return true;
        if (*nb >= ELL_MAXB) return false;
        B[(*nb)++] = g;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return true;
    if (!collect_bases(e->data.function.head, x, B, nb)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!collect_bases(e->data.function.args[i], x, B, nb)) return false;
    return true;
}

/* Copy of e with Power[B[i], k/2] -> (coef[i] y)^k. */
static Expr* subst_bases(Expr* e, const Expr* x, Expr** B, Expr** coef, int nb, Expr* y) {
    int k;
    if (is_half_power(e, x, &k)) {
        for (int i = 0; i < nb; i++)
            if (expr_eq(B[i], e->data.function.args[0])) {
                Expr* kE = E_int(k);
                Expr* r = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_function(expr_new_symbol("Times"),
                                   (Expr*[]){ expr_copy(coef[i]), expr_copy(y) }, 2), kE }, 2);
                return r;
            }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr** args = n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++)
        args[i] = subst_bases(e->data.function.args[i], x, B, coef, nb, y);
    Expr* r = expr_new_function(subst_bases(e->data.function.head, x, B, coef, nb, y), args, n);
    free(args);
    return r;
}

/* true when the polynomial p has a real root strictly inside (lo, hi) */
static bool poly_root_inside(Expr* p, Expr* x, Pt* lo, Pt* hi) {
    if (!contains_sym(p, x)) return false;
    Expr* rs = ell_x("#2 /. NSolve[#1 == 0, #2]", 2, p, x);
    bool in = !ell_head_is(rs, "List");
    for (size_t i = 0; !in && i < ell_len(rs); i++) {
        Expr* r = ell_arg(rs, i);
        double d = ell_num(r);
        if (isnan(d)) continue;
        double tol = 1e-9 * (1 + fabs(d));
        if ((lo->inf || d > lo->d + tol) && (hi->inf || d < hi->d - tol)) in = true;
    }
    expr_free(rs);
    return in;
}

/* f -> A + W/Sqrt[P], P > 0 on the interval.  Outputs owned. */
static bool ell_normalise(Expr* f, Expr* x, Pt* lo, Pt* hi, Expr* xmid,
                          Expr** Aout, Expr** Wout, Expr** Pout) {
    Expr* B[ELL_MAXB]; int nb = 0;
    if (!collect_bases(f, x, B, &nb) || nb == 0) return false;
    Expr* info[ELL_MAXB];
    for (int i = 0; i < nb; i++)
        info[i] = ell_x("Module[{tg = Together[#1], fl}, fl = FactorSquareFreeList[Expand[Numerator[tg] Denominator[tg]]]; "
                        "{Denominator[tg], Times @@ Map[Function[{pq}, pq[[1]]^Mod[pq[[2]], 2]], fl], "
                        "Times @@ Map[Function[{pq}, pq[[1]]^Quotient[pq[[2]], 2]], fl]}]", 1, B[i]);
    bool ok = true;
    for (int i = 0; i < nb; i++) if (ell_len(info[i]) != 3) ok = false;
    Expr* coef[ELL_MAXB] = { NULL };
    Expr* P = NULL;
    if (ok) {
        Expr* K1 = ell_arg(info[0], 1);
        Expr* k1m = ell_x("#1 /. #2 -> #3", 3, K1, x, xmid);
        double sK = ell_num(k1m);
        expr_free(k1m);
        if (isnan(sK) || sK == 0) ok = false;
        else {
            Expr* sKE = E_int(sK > 0 ? 1 : -1);
            P = ell_x("Expand[#1 #2]", 2, sKE, K1);
            for (int i = 0; ok && i < nb; i++) {
                Expr* Dn = ell_arg(info[i], 0);
                Expr* Ki = ell_arg(info[i], 1);
                Expr* Si = ell_arg(info[i], 2);
                if (poly_root_inside(Dn, x, lo, hi) || poly_root_inside(Si, x, lo, hi)) { ok = false; break; }
                Expr* lam = ell_x("Together[#1/#2]", 2, Ki, K1);
                if (contains_sym(lam, x)) { expr_free(lam); ok = false; break; }
                Expr* rho = ell_x("Sign[(#1/#2) /. #3 -> #4]", 4, Si, Dn, x, xmid);
                coef[i] = ell_x("Integrate`Elliptic`Private`simp[#1 Sqrt[#2 #3] #4/#5]", 5, rho, lam, sKE, Si, Dn);
                expr_free(rho); expr_free(lam);
            }
            expr_free(sKE);
            if (ok && poly_root_inside(P, x, lo, hi)) ok = false;
        }
    }
    Expr* A = NULL; Expr* W = NULL;
    if (ok) {
        Expr* deg = ell_x("Exponent[#1, #2]", 2, P, x);
        if (!(deg->type == EXPR_INTEGER && (deg->data.integer == 3 || deg->data.integer == 4))) ok = false;
        expr_free(deg);
    }
    if (ok) {
        Expr* y = expr_new_symbol("Integrate`Elliptic`Private`y");
        Expr* F = subst_bases(f, x, B, coef, nb, y);
        Expr* Fe = evaluate(F);
        expr_free(F);
        /* (n0 + n1 y)/(d0 + d1 y) -> (e0 + e1 y), y^2 = P */
        Expr* r = ell_x("Module[{tg = Together[#1], nn, dd, n0, n1, d0, d1, den}, "
                        "nn = PolynomialRemainder[Expand[Numerator[tg]], #3^2 - #2, #3]; "
                        "dd = PolynomialRemainder[Expand[Denominator[tg]], #3^2 - #2, #3]; "
                        "n0 = Coefficient[nn, #3, 0]; n1 = Coefficient[nn, #3, 1]; "
                        "d0 = Coefficient[dd, #3, 0]; d1 = Coefficient[dd, #3, 1]; "
                        "den = Expand[d0^2 - d1^2 #2]; "
                        "{Together[(n0 d0 - n1 d1 #2)/den], Together[(n1 d0 - n0 d1) #2/den]}]",
                        3, Fe, P, y);
        expr_free(Fe);
        if (ell_len(r) == 2) {
            A = expr_copy(ell_arg(r, 0));
            W = expr_copy(ell_arg(r, 1));
            if (contains_sym(A, y) || contains_sym(W, y)) ok = false;
            Expr* rat = ell_x("PolynomialQ[Numerator[#1], #3] && PolynomialQ[Denominator[#1], #3] && "
                              "PolynomialQ[Numerator[#2], #3] && PolynomialQ[Denominator[#2], #3]", 3, A, W, x);
            if (!ell_true(rat)) ok = false;
            if (ok && ell_zero(W)) ok = false;
        } else ok = false;
        expr_free(r);
        expr_free(y);
    }
    for (int i = 0; i < nb; i++) { expr_free(info[i]); if (coef[i]) expr_free(coef[i]); }
    if (!ok) {
        if (A) expr_free(A);
        if (W) expr_free(W);
        if (P) expr_free(P);
        return false;
    }
    *Aout = A; *Wout = W; *Pout = P;
    return true;
}

/* Moderate clean-up of a closed form: Simplify, then the purely algebraic
 * terms (no elliptic / inverse-trig / log head) are combined and reduced to a
 * single radical when RootReduce can do so without a Root object. */
static Expr* ell_finish(Expr* r) {
    return ell_x(
        "Module[{s = If[LeafCount[#1] <= 600, TimeConstrained[Simplify[#1], 6, #1], #1], tl, alg, rest, ar}, "
        "s = If[FreeQ[s, Root], s, #1]; "
        "tl = If[Head[s] === Plus, List @@ s, {s}]; "
        "alg = Select[tl, Function[{u}, FreeQ[u, EllipticF | EllipticE | EllipticK | EllipticPi | "
        "ArcTan | ArcTanh | ArcSin | Log | Pi]]]; "
        "rest = Complement[tl, alg]; "
        "ar = TimeConstrained[RootReduce[Plus @@ alg], 4, Plus @@ alg]; "
        "If[FreeQ[ar, Root] && LeafCount[ar] < LeafCount[Plus @@ alg], ar + Plus @@ rest, s]]", 1, r);
}

/* -------------------------------------------------------------------------
 * Driver.
 * ---------------------------------------------------------------------- */
static Expr* ell_core(Expr* f, Expr* x, Pt* lo, Pt* hi) {
    Expr* xmid;
    if (!lo->inf && !hi->inf) xmid = ell_x("(#1 + #2)/2", 2, lo->v, hi->v);
    else if (!lo->inf)        xmid = ell_x("#1 + 1", 1, lo->v);
    else if (!hi->inf)        xmid = ell_x("#1 - 1", 1, hi->v);
    else                      xmid = E_int(0);
    Expr *A = NULL, *W = NULL, *P = NULL;
    Expr* result = NULL;
    if (!ell_normalise(f, x, lo, hi, xmid, &A, &W, &P)) {
        ELL_DBG("normalise declined");
        expr_free(xmid);
        return NULL;
    }
    int nopt = 1;
    long best_cost = 0;
    for (int pick = 0; pick < nopt; pick++) {
        XSys S; memset(&S, 0, sizeof S);
        Expr* r = NULL;
        bool ok = x_pair(&S, P, x, pick, &nopt);
        if (!ok) ELL_DBG("real pairing declined");
        if (ok && !(ok = x_partial(&S, W, x))) ELL_DBG("partial fractions declined");
        for (int p = 0; ok && p < S.np; p++)      /* poles must lie off [lo, hi] */
            if ((lo->inf || S.pd[p] >= lo->d - 1e-12) && (hi->inf || S.pd[p] <= hi->d + 1e-12)) ok = false;
        if (ok) {
            r = x_integral(&S, lo, hi);
            if (!r) ELL_DBG("reduction declined");
        }
        xsys_free(&S);
        if (!r) continue;
        Expr* rs = ell_finish(r);
        expr_free(r);
        r = rs;
        Expr* c = ell_x("LeafCount[#1]", 1, r);
        long cost = c->type == EXPR_INTEGER ? (long)c->data.integer : 1000000L;
        expr_free(c);
        if (!result || cost < best_cost) {
            if (result) expr_free(result);
            result = r; best_cost = cost;
        } else expr_free(r);
    }
    if (result && !ell_zero(A)) {
        Expr* lv = lo->inf ? ell_x("-Infinity", 0) : expr_copy(lo->v);
        Expr* hv = hi->inf ? ell_x("Infinity", 0) : expr_copy(hi->v);
        Expr* ai = ell_x("Integrate[#1, {#2, #3, #4}]", 4, A, x, lv, hv);
        expr_free(lv); expr_free(hv);
        if (ell_head_is(ai, "Integrate") || !ell_true(ell_x("NumericQ[#1]", 1, ai))) {
            expr_free(result); result = NULL;
        } else acc_add(&result, expr_copy(ai));
        expr_free(ai);
    }
    expr_free(A); expr_free(W); expr_free(P); expr_free(xmid);
    return result;
}

/* The reference value of Int_lo^hi f dx (lo < hi) by double-exponential
 * quadrature at 90 digits (tanh-sinh on a finite interval, exp-sinh on a
 * half-line, the whole line split at 0), at step sizes 1/32 and 1/64; NULL
 * unless the two agree to 10^-12 (halving h roughly squares the error of a
 * double-exponential rule, so the finer value is then good to ~10^-24).
 * Self-contained on purpose: NIntegrate at WorkingPrecision 40 was found to return 5.0111 for
 * Int_-1^3 x^2/Sqrt[(x+1)(3-x)(x^2-2x+2)] (= 5.2704), to drop to machine
 * precision when given intermediate points, and to lose digits after a
 * trigonometric change of variables -- none of which a 10^-20 check can
 * absorb.  The square-root endpoint singularities of these integrands need
 * no special handling: the nodes approach the ends double-exponentially and
 * are placed as a + (b - a) q with q computed directly (no cancellation). */
static Expr* ell_ts(Expr* f, Expr* x, Pt* lo, Pt* hi, Expr* h) {
    if (!lo->inf && !hi->inf)
        return ell_x("Integrate`Elliptic`Private`tsF[#1, #2, #3, #4, #5]", 5, f, x, lo->v, hi->v, h);
    if (!lo->inf)
        return ell_x("Integrate`Elliptic`Private`tsH[#1, #2, #3, 1, #4]", 4, f, x, lo->v, h);
    if (!hi->inf)
        return ell_x("Integrate`Elliptic`Private`tsH[#1, #2, #3, -1, #4]", 4, f, x, hi->v, h);
    return ell_x("Integrate`Elliptic`Private`tsH[#1, #2, 0, 1, #3] + "
                 "Integrate`Elliptic`Private`tsH[#1, #2, 0, -1, #3]", 3, f, x, h);
}

static Expr* ell_reference(Expr* f, Expr* x, Pt* lo, Pt* hi) {
    Expr* h1 = ell_x("1/32", 0);
    Expr* h2 = ell_x("1/64", 0);
    Expr* r1 = ell_ts(f, x, lo, hi, h1);
    Expr* r2 = ell_ts(f, x, lo, hi, h2);
    bool agree = ell_true(ell_x("NumericQ[#1] && NumericQ[#2] && Abs[#1 - #2] < 10^-12 (1 + Abs[#2])", 2, r1, r2));
    expr_free(h1); expr_free(h2); expr_free(r1);
    if (!agree) { expr_free(r2); return NULL; }
    return r2;
}

/* N[r, 30] against the reference to a relative 10^-20. */
static bool ell_verify(Expr* r, Expr* f, Expr* x, Pt* lo, Pt* hi, bool fwd) {
    Expr* ref = ell_reference(f, x, lo, hi);
    if (!ref) return false;
    Expr* sg = E_int(fwd ? 1 : -1);
    Expr* ok = ell_x("With[{v = N[#1, 30], w = #3 #2}, NumericQ[w] && NumericQ[v] && "
                     "Abs[v - w] <= 10^-20 (1 + Abs[w])]", 3, r, ref, sg);
    expr_free(sg);
    expr_free(ref);
    return ell_true(ok);
}

/* The complete integrals as quarter-period integrals of a trigonometric
 * quadratic: with Q(x) = A + B Sin[x]^2 (Cos^2 = 1 - Sin^2 first),
 *     Int_0^{Pi/2} Q^(-1/2) dx = EllipticK[-B/A]/Sqrt[A],
 *     Int_0^{Pi/2} Q^(1/2)  dx = Sqrt[A] EllipticE[-B/A],
 * for A > 0 (the substitution-free definition of K and E in the parameter
 * m = -B/A).  Q is even and Pi-periodic in x, so [0, Pi] doubles and
 * [0, 2 Pi] quadruples.  Gauss's lemniscate AGM integral
 * 1/Sqrt[Cos^2 + 2 Sin^2] is K(-1); the ellipse perimeter Sqrt[1 - m Sin^2]
 * is E(m).  The terms of Q are classified natively (constant, c Sin^2,
 * c Cos^2); A > 0 is asked of Greater, then of Simplify under $Assumptions. */
static Expr* tq_eval(Expr* e) { Expr* r = evaluate(e); expr_free(e); return r; }
static Expr* tq_fn2(const char* h, Expr* x, Expr* y) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ x, y }, 2);
}
static bool tq_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}
/* e == Trig[x]^2 for the named trig head. */
static bool tq_trig_sq(const Expr* e, const Expr* x, const char* trig) {
    return tq_head(e, SYM_Power) && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER && e->data.function.args[1]->data.integer == 2
        && tq_head(e->data.function.args[0], trig) && e->data.function.args[0]->data.function.arg_count == 1
        && expr_eq(e->data.function.args[0]->data.function.args[0], x);
}
/* b / Pi as a small integer multiple of 1/2 (1 -> Pi/2, 2 -> Pi, 4 -> 2 Pi). */
static int tq_period_multiple(const Expr* b) {
    static const int64_t num[3] = { 1, 1, 2 }, den[3] = { 2, 1, 1 };
    static const int k[3] = { 1, 2, 4 };
    for (int i = 0; i < 3; i++) {
        Expr* ref = tq_eval(tq_fn2(SYM_Times, make_rational(num[i], den[i]), expr_new_symbol(SYM_Pi)));
        bool m = expr_eq(ref, b);
        expr_free(ref);
        if (m) return k[i];
    }
    return 0;
}

static Expr* tq_core(Expr* f, Expr* x, Expr* a, Expr* b);

/* c (A + B Sin^2)^(+-1/2), c free of x: the constant comes out. */
Expr* integrate_elliptic_trig_quadratic(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!tq_head(f, SYM_Times)) return tq_core(f, x, a, b);
    size_t n = f->data.function.arg_count, iv = n, nv = 0;
    for (size_t i = 0; i < n; i++)
        if (contains_sym(f->data.function.args[i], x)) { iv = i; nv++; }
    if (nv != 1) return NULL;
    Expr* r = tq_core(f->data.function.args[iv], x, a, b);
    if (!r) return NULL;
    Expr** v = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) v[i] = i == iv ? r : expr_copy(f->data.function.args[i]);
    Expr* out = tq_eval(expr_new_function(expr_new_symbol(SYM_Times), v, n));
    free(v);
    return out;
}

static Expr* tq_core(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    int k = tq_period_multiple(b);
    if (!k || !tq_head(f, SYM_Power) || f->data.function.arg_count != 2) return NULL;
    const Expr* ex = f->data.function.args[1];
    Expr* half = make_rational(1, 2), *mhalf = make_rational(-1, 2);
    bool inv = expr_eq(ex, mhalf), dir = expr_eq(ex, half);
    expr_free(half); expr_free(mhalf);
    if (!inv && !dir) return NULL;
    /* base = A + B Sin^2, reading Cos^2 as 1 - Sin^2: classify the terms. */
    Expr* q = tq_eval(expr_new_function(expr_new_symbol(SYM_Expand),
                      (Expr*[]){ expr_copy(f->data.function.args[0]) }, 1));
    size_t nt = tq_head(q, SYM_Plus) ? q->data.function.arg_count : 1;
    Expr* const* t = tq_head(q, SYM_Plus) ? q->data.function.args : &q;
    Expr** A = malloc(sizeof(Expr*) * 2 * nt);
    Expr** B = malloc(sizeof(Expr*) * 2 * nt);
    size_t na = 0, nb = 0;
    bool ok = true, any = false;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* ti = t[i];
        if (!contains_sym(ti, x)) { A[na++] = expr_copy((Expr*)ti); continue; }
        size_t nf = tq_head(ti, SYM_Times) ? ti->data.function.arg_count : 1;
        Expr* const* fa = tq_head(ti, SYM_Times) ? ti->data.function.args : (Expr* const*)&ti;
        int which = 0;                  /* 1 = Sin^2, 2 = Cos^2 */
        Expr** c = malloc(sizeof(Expr*) * nf);
        size_t nc = 0;
        for (size_t j = 0; j < nf && ok; j++) {
            if (!contains_sym(fa[j], x)) { c[nc++] = expr_copy(fa[j]); continue; }
            int w = tq_trig_sq(fa[j], x, SYM_Sin) ? 1 : tq_trig_sq(fa[j], x, SYM_Cos) ? 2 : 0;
            if (!w || which) ok = false; else which = w;
        }
        if (ok) {
            Expr* co = expr_new_function(expr_new_symbol(SYM_Times), c, nc);
            any = true;
            if (which == 1) B[nb++] = co;
            else {                      /* c Cos^2 = c - c Sin^2 */
                A[na++] = expr_copy(co);
                B[nb++] = tq_fn2(SYM_Times, expr_new_integer(-1), co);
            }
        } else {
            for (size_t j = 0; j < nc; j++) expr_free(c[j]);
        }
        free(c);
    }
    expr_free(q);
    Expr* out = NULL;
    if (ok && any) {
        Expr* Av = tq_eval(expr_new_function(expr_new_symbol(SYM_Plus), A, na));
        Expr* Bv = tq_eval(expr_new_function(expr_new_symbol(SYM_Plus), B, nb));
        bool pos = ell_true(tq_eval(tq_fn2(SYM_Greater, expr_copy(Av), expr_new_integer(0))))
                || ell_true(tq_eval(expr_new_function(expr_new_symbol(SYM_Simplify),
                       (Expr*[]){ tq_fn2(SYM_Greater, expr_copy(Av), expr_new_integer(0)) }, 1)));
        bool bz = Bv->type == EXPR_INTEGER && Bv->data.integer == 0;
        if (pos && !bz) {
            /* m = -B/A */
            Expr* m = tq_fn2(SYM_Times, tq_fn2(SYM_Times, expr_new_integer(-1), Bv),
                             tq_fn2(SYM_Power, expr_copy(Av), expr_new_integer(-1)));
            Bv = NULL;
            Expr* sq = tq_fn2(SYM_Power, Av, inv ? make_rational(-1, 2) : make_rational(1, 2));
            Av = NULL;
            Expr* ell = expr_new_function(expr_new_symbol(inv ? SYM_EllipticK : SYM_EllipticE),
                                          (Expr*[]){ m }, 1);
            out = tq_eval(expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_new_integer(k), sq, ell }, 3));
        }
        if (Av) expr_free(Av);
        if (Bv) expr_free(Bv);
    } else {
        for (size_t i = 0; i < na; i++) expr_free(A[i]);
        for (size_t i = 0; i < nb; i++) expr_free(B[i]);
    }
    free(A); free(B);
    return out;
}

/* ---- indefinite: incomplete elliptic integrals ----------------------------
 *
 *   Int (A + B Sin[k x]^2)^(-1/2) dx = EllipticF[k x, -B/A]/(k Sqrt[A])
 *   Int (A + B Sin[k x]^2)^(1/2)  dx = Sqrt[A] EllipticE[k x, -B/A]/k
 *   Int (A + B Cos[k x])^(-1/2)   dx = 2 EllipticF[k x/2, 2B/(A+B)]/(k Sqrt[A+B])
 *   Int (A + B Cos[k x])^(1/2)    dx = 2 Sqrt[A+B] EllipticE[k x/2, 2B/(A+B)]/k
 *
 * (Cos^2 is read as 1 - Sin^2; Cos u = 1 - 2 Sin[u/2]^2.)  Each is the
 * definition of F and E in the parameter m, d/dphi F = (1 - m Sin^2)^(-1/2),
 * so the antiderivative differentiates back on the principal branch (A > 0).
 * The pendulum's first integral Int (C + 2 Cos[th])^(-1/2) dth closes here
 * instead of grinding through the cascade (60 s, unevaluated). */

/* e == g[k x] (or g[x]) with k free of x: sets *k (owned, or NULL for 1). */
static bool tq_lin_arg(const Expr* e, const Expr* x, Expr** k) {
    *k = NULL;
    if (expr_eq(e, x)) return true;
    if (tq_head(e, SYM_Times) && e->data.function.arg_count == 2
        && expr_eq(e->data.function.args[1], x) && !contains_sym(e->data.function.args[0], x)) {
        *k = expr_copy(e->data.function.args[0]);
        return true;
    }
    return false;
}

/* f is Trig[u]^pw with u a linear argument; returns the argument's k via *k
 * and the argument itself via *u (borrowed). */
static bool tq_trig_pow(const Expr* e, const Expr* x, const char* trig, int64_t pw,
                        const Expr** u) {
    (void)x;
    const Expr* t = e;
    if (pw != 1) {
        if (!(tq_head(e, SYM_Power) && e->data.function.arg_count == 2
              && e->data.function.args[1]->type == EXPR_INTEGER
              && e->data.function.args[1]->data.integer == pw)) return false;
        t = e->data.function.args[0];
    }
    if (!tq_head(t, trig) || t->data.function.arg_count != 1) return false;
    *u = t->data.function.args[0];
    return true;
}

/* Cos[u] (A + B Cos[u])^(-1/2) = ((A + B Cos[u])^(1/2) - A (A + B Cos[u])^(-1/2))/B:
 * the first-degree numerator reduces to the two kernel forms below.  Without
 * it the term Cos[u]/Sqrt[-Cos[u]] of G&R 2.619.2 fell to a later stage whose
 * answer held on one sign interval of Cos only.  NULL unless f has exactly
 * that shape (constant factors allowed) and both pieces close. */
static Expr* tq_cos_times_invsqrt(Expr* f, Expr* x) {
    if (!tq_head(f, SYM_Times)) return NULL;
    size_t n = f->data.function.arg_count;
    const Expr *cosf = NULL, *pw = NULL;
    size_t nx = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* fi = f->data.function.args[i];
        if (!contains_sym(fi, x)) continue;
        nx++;
        if (tq_head(fi, SYM_Cos) && fi->data.function.arg_count == 1) cosf = fi;
        else if (tq_head(fi, SYM_Power) && fi->data.function.arg_count == 2) pw = fi;
    }
    if (nx != 2 || !cosf || !pw) return NULL;
    Expr* mhalf = make_rational(-1, 2);
    bool inv = expr_eq(pw->data.function.args[1], mhalf);
    expr_free(mhalf);
    if (!inv) return NULL;
    Expr* Q = pw->data.function.args[0];
    Expr* B = tq_eval(expr_new_function(expr_new_symbol(SYM_Coefficient),
                     (Expr*[]){ expr_copy(Q), expr_copy((Expr*)cosf) }, 2));
    Expr* A = tq_eval(expr_new_function(expr_new_symbol(SYM_Expand),
                     (Expr*[]){ tq_fn2(SYM_Plus, expr_copy(Q),
                                tq_fn2(SYM_Times, expr_new_integer(-1),
                                       tq_fn2(SYM_Times, expr_copy(B), expr_copy((Expr*)cosf)))) }, 1));
    Expr* out = NULL;
    bool ok = B && A && !contains_sym(A, x) && !contains_sym(B, x) &&
              !(B->type == EXPR_INTEGER && B->data.integer == 0);
    if (ok) {
        /* constant cofactor c of Cos[u] Q^(-1/2) */
        Expr** cs = malloc(sizeof(Expr*) * (n + 1));
        size_t w = 0;
        for (size_t i = 0; i < n; i++)
            if (!contains_sym(f->data.function.args[i], x)) cs[w++] = expr_copy(f->data.function.args[i]);
        cs[w++] = tq_fn2(SYM_Power, expr_copy(B), expr_new_integer(-1));
        Expr* c = tq_eval(expr_new_function(expr_new_symbol(SYM_Times), cs, w));
        free(cs);
        Expr* t1 = tq_eval(tq_fn2(SYM_Times, expr_copy(c),
                                  tq_fn2(SYM_Power, expr_copy(Q), make_rational(1, 2))));
        Expr* t2 = tq_eval(tq_fn2(SYM_Times,
                                  tq_fn2(SYM_Times, expr_new_integer(-1), tq_fn2(SYM_Times, c, expr_copy(A))),
                                  tq_fn2(SYM_Power, expr_copy(Q), make_rational(-1, 2))));
        bool t2zero = t2 && t2->type == EXPR_INTEGER && t2->data.integer == 0;
        Expr* I1 = integrate_elliptic_trig_indefinite(t1, x);
        Expr* I2 = !I1 ? NULL : t2zero ? expr_new_integer(0)
                                       : integrate_elliptic_trig_indefinite(t2, x);
        if (I1 && I2) out = tq_eval(tq_fn2(SYM_Plus, I1, I2));
        else { if (I1) expr_free(I1); if (I2) expr_free(I2); }
        expr_free(t1); expr_free(t2);
    }
    if (A) expr_free(A);
    if (B) expr_free(B);
    return out;
}

Expr* integrate_elliptic_trig_indefinite(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;
    {
        Expr* r = tq_cos_times_invsqrt(f, x);
        if (r) return r;
    }
    /* c Q^e: pull a constant factor */
    const Expr* core = f;
    Expr* cfac = NULL;
    if (tq_head(f, SYM_Times)) {
        size_t n = f->data.function.arg_count, iv = n, nv = 0;
        for (size_t i = 0; i < n; i++)
            if (contains_sym(f->data.function.args[i], x)) { iv = i; nv++; }
        if (nv != 1) return NULL;
        core = f->data.function.args[iv];
        Expr** cs = malloc(sizeof(Expr*) * n);
        size_t w = 0;
        for (size_t i = 0; i < n; i++) if (i != iv) cs[w++] = expr_copy(f->data.function.args[i]);
        cfac = expr_new_function(expr_new_symbol(SYM_Times), cs, w);
        free(cs);
    }
    Expr* out = NULL;
    if (!tq_head(core, SYM_Power) || core->data.function.arg_count != 2) goto done;
    Expr* half = make_rational(1, 2), *mhalf = make_rational(-1, 2);
    bool inv = expr_eq(core->data.function.args[1], mhalf);
    bool dir = expr_eq(core->data.function.args[1], half);
    expr_free(half); expr_free(mhalf);
    if (!inv && !dir) goto done;
    {
        Expr* q = tq_eval(expr_new_function(expr_new_symbol(SYM_Expand),
                          (Expr*[]){ expr_copy(core->data.function.args[0]) }, 1));
        size_t nt = tq_head(q, SYM_Plus) ? q->data.function.arg_count : 1;
        Expr* const* t = tq_head(q, SYM_Plus) ? q->data.function.args : &q;
        Expr** A = malloc(sizeof(Expr*) * 2 * nt);
        Expr** B = malloc(sizeof(Expr*) * 2 * nt);
        size_t na = 0, nb = 0;
        int kind = 0;                  /* 1: Sin^2/Cos^2 of u, 2: Cos[u] */
        const Expr* u0 = NULL;
        bool ok = true;
        for (size_t i = 0; i < nt && ok; i++) {
            const Expr* ti = t[i];
            if (!contains_sym(ti, x)) { A[na++] = expr_copy((Expr*)ti); continue; }
            size_t nf = tq_head(ti, SYM_Times) ? ti->data.function.arg_count : 1;
            Expr* const* fa = tq_head(ti, SYM_Times) ? ti->data.function.args : (Expr* const*)&ti;
            int which = 0; const Expr* u = NULL;
            Expr** c = malloc(sizeof(Expr*) * nf);
            size_t nc = 0;
            for (size_t j = 0; j < nf && ok; j++) {
                if (!contains_sym(fa[j], x)) { c[nc++] = expr_copy(fa[j]); continue; }
                int w = 0;
                if (tq_trig_pow(fa[j], x, SYM_Sin, 2, &u)) w = 1;
                else if (tq_trig_pow(fa[j], x, SYM_Cos, 2, &u)) w = 2;
                else if (tq_trig_pow(fa[j], x, SYM_Cos, 1, &u)) w = 3;
                if (!w || which) ok = false; else which = w;
            }
            int k_needed = which == 3 ? 2 : 1;
            if (ok && (kind == 0 || kind == k_needed) && (!u0 || expr_eq(u0, u))) {
                kind = k_needed; u0 = u;
                Expr* co = expr_new_function(expr_new_symbol(SYM_Times), c, nc);
                if (which == 2) {          /* c Cos^2 = c - c Sin^2 */
                    A[na++] = expr_copy(co);
                    B[nb++] = tq_fn2(SYM_Times, expr_new_integer(-1), co);
                } else B[nb++] = co;
            } else {
                for (size_t j = 0; j < nc; j++) expr_free(c[j]);
                ok = false;
            }
            free(c);
        }
        Expr* k = NULL;
        if (ok && kind && u0 && tq_lin_arg(u0, x, &k)) {
            Expr* Av = tq_eval(expr_new_function(expr_new_symbol(SYM_Plus), A, na));
            Expr* Bv = tq_eval(expr_new_function(expr_new_symbol(SYM_Plus), B, nb));
            na = nb = 0;
            Expr *amp, *m, *S, *scale;     /* phi, parameter, Sqrt base, prefactor */
            if (kind == 1) {
                amp = expr_copy((Expr*)u0);
                m = tq_fn2(SYM_Times, tq_fn2(SYM_Times, expr_new_integer(-1), expr_copy(Bv)),
                           tq_fn2(SYM_Power, expr_copy(Av), expr_new_integer(-1)));
                S = expr_copy(Av);
                scale = expr_new_integer(1);
            } else {
                Expr* AB = tq_fn2(SYM_Plus, expr_copy(Av), expr_copy(Bv));
                amp = tq_fn2(SYM_Times, make_rational(1, 2), expr_copy((Expr*)u0));
                m = tq_fn2(SYM_Times, tq_fn2(SYM_Times, expr_new_integer(2), expr_copy(Bv)),
                           tq_fn2(SYM_Power, expr_copy(AB), expr_new_integer(-1)));
                S = AB;
                scale = expr_new_integer(2);
            }
            /* Sqrt[S (1 - m Sin[amp]^2)] = Sqrt[S] Sqrt[1 - m Sin[amp]^2] only
             * for S > 0.  Otherwise (Sqrt[-Cos[x]]: S = -1) the factor that
             * relates the integrand to the elliptic kernel is the LOCALLY
             * CONSTANT ratio Sqrt[base]/Sqrt[base/S] -- +/-I or a real constant
             * on each interval where base keeps its sign -- so the
             * antiderivative is ratio * E (resp. F / ratio), which
             * differentiates back to the integrand on every interval.  The
             * bare Sqrt[S] form was right only where base/S > 0: Integrate[
             * Sqrt[-Cos[x]], x] was 2 I EllipticE[x/2, 2], off by a sign
             * wherever Cos[x] < 0 (G&R 2.619.2, 2.619.11). */
            bool S_pos = ell_true(tq_eval(tq_fn2(SYM_Greater, expr_copy(S), expr_new_integer(0))))
                || ell_true(tq_eval(expr_new_function(expr_new_symbol(SYM_Simplify),
                       (Expr*[]){ tq_fn2(SYM_Greater, expr_copy(S), expr_new_integer(0)) }, 1)));
            Expr* sq;
            if (S_pos) {
                sq = tq_fn2(SYM_Power, S, inv ? make_rational(-1, 2) : make_rational(1, 2));
            } else {
                Expr* base = core->data.function.args[0];
                Expr* nb = tq_eval(expr_new_function(expr_new_symbol(SYM_Together),
                               (Expr*[]){ tq_fn2(SYM_Times, expr_copy(base),
                                                 tq_fn2(SYM_Power, S, expr_new_integer(-1))) }, 1));
                Expr* rb = tq_fn2(SYM_Power, expr_copy(base), make_rational(1, 2));
                Expr* rn = tq_fn2(SYM_Power, nb, make_rational(1, 2));
                sq = inv ? tq_fn2(SYM_Times, rn, tq_fn2(SYM_Power, rb, expr_new_integer(-1)))
                         : tq_fn2(SYM_Times, rb, tq_fn2(SYM_Power, rn, expr_new_integer(-1)));
            }
            Expr* ell = expr_new_function(expr_new_symbol(inv ? SYM_EllipticF : SYM_EllipticE),
                                          (Expr*[]){ amp, m }, 2);
            Expr* kinv = k ? tq_fn2(SYM_Power, k, expr_new_integer(-1)) : expr_new_integer(1);
            out = tq_eval(expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ scale, sq, ell, kinv,
                                     cfac ? expr_copy(cfac) : expr_new_integer(1) }, 5));
            expr_free(Av); expr_free(Bv);
        } else if (k) expr_free(k);
        for (size_t i = 0; i < na; i++) expr_free(A[i]);
        for (size_t i = 0; i < nb; i++) expr_free(B[i]);
        free(A); free(B);
        expr_free(q);
    }
done:
    if (cfac) expr_free(cfac);
    return out;
}

Expr* integrate_elliptic_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    (void)assumptions;
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    {
        Expr* tq = integrate_elliptic_trig_quadratic(f, x, a, b);
        if (tq) return tq;
    }
    if (!contains_sym(f, x) || contains_sym(a, x) || contains_sym(b, x)) return NULL;
    if (ell_depth > 0) return NULL;
    unsigned long fired0 = mth_msg_fired_count();
    mth_msg_suppress_push();
    arith_warnings_mute_push();
    ell_depth++;
    Expr* result = NULL;
    Pt lo, hi;
    bool okp = pt_from(a, &lo);
    if (okp && !pt_from(b, &hi)) { pt_free(&lo); okp = false; }
    if (okp) {
        /* only symbol-free real limits with lo < hi; a reversed interval is
         * the negative of the forward one */
        bool numeric = ell_true(ell_x("(NumericQ[#1] || #1 === -Infinity || #1 === Infinity) && "
                                      "(NumericQ[#2] || #2 === -Infinity || #2 === Infinity) && "
                                      "FreeQ[{#1, #2}, Complex]", 2, a, b));
        bool fwd = lo.d < hi.d;
        if (numeric && lo.d != hi.d && !(lo.inf && lo.inf == hi.inf)) {
            Expr* r = fwd ? ell_core(f, x, &lo, &hi) : ell_core(f, x, &hi, &lo);
            if (r) {
                Expr* s = fwd ? expr_copy(r) : ell_x("-#1", 1, r);
                expr_free(r);
                if (getenv("MATHILDA_ELL_DEBUG")) {
                    char* str = expr_to_string(s);
                    fprintf(stderr, "[ell] candidate: %s\n", str);
                    free(str);
                }
                if (fwd ? ell_verify(s, f, x, &lo, &hi, true) : ell_verify(s, f, x, &hi, &lo, false)) result = s;
                else expr_free(s);
            }
        }
        pt_free(&lo); pt_free(&hi);
    }
    ell_depth--;
    arith_warnings_mute_pop();
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired0);
    return result;
}

/* -------------------------------------------------------------------------
 * Integrate`Elliptic[f, {x, a, b}, Assumptions -> ...]
 * ---------------------------------------------------------------------- */
Expr* builtin_integrate_elliptic(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!ell_head_is(spec, "List") || spec->data.function.arg_count != 3) return NULL;
    Expr* x = spec->data.function.args[0];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if ((ell_head_is(opt, "Rule") || ell_head_is(opt, "RuleDelayed")) &&
            opt->data.function.arg_count == 2 &&
            ell_is_sym(opt->data.function.args[0], "Assumptions")) {
            assumptions = opt->data.function.args[1];
            continue;
        }
        return NULL;
    }
    return integrate_elliptic_try(f, x, spec->data.function.args[1],
                                  spec->data.function.args[2], assumptions);
}

void integrate_elliptic_init(void) {
    tc_register_unwind_counter(&ell_depth);
    /* the coefficient simplifier: Simplify (only on expressions of modest
     * size, and under a 2-second budget -- a timed-out Simplify is both slow
     * and leaky), kept only when it does not trade radicals for Root objects
     * (which neither print well nor compare reliably) */
    Expr* def = parse_expression(
        "Integrate`Elliptic`Private`simp[Integrate`Elliptic`Private`e_] := "
        "With[{Integrate`Elliptic`Private`r = If[LeafCount[Integrate`Elliptic`Private`e] > 250, Integrate`Elliptic`Private`e, TimeConstrained[Simplify[Integrate`Elliptic`Private`e], 2, Integrate`Elliptic`Private`e]]}, "
        "If[FreeQ[Integrate`Elliptic`Private`r, Root], Integrate`Elliptic`Private`r, "
        "Integrate`Elliptic`Private`e]]");
    if (def) { Expr* r = evaluate(def); expr_free(r); expr_free(def); }
    /* the double-exponential reference quadratures (see ell_reference) */
    static const char* const ts_defs[] = {
        "Integrate`Elliptic`Private`tsF[Integrate`Elliptic`Private`f_, Integrate`Elliptic`Private`x_, "
        "Integrate`Elliptic`Private`a_, Integrate`Elliptic`Private`b_, Integrate`Elliptic`Private`h_] := "
        "Module[{k = 0, tot = 0, t, u, q, w, v, go = True}, While[go, "
        "t = k Integrate`Elliptic`Private`h; u = N[Pi/2 Sinh[t], 90]; q = 1/(1 + Exp[2 u]); "
        "w = N[(Integrate`Elliptic`Private`b - Integrate`Elliptic`Private`a) Pi q (1 - q) Cosh[t], 90]; "
        "v = If[k == 0, w N[Integrate`Elliptic`Private`f /. Integrate`Elliptic`Private`x -> "
        "(Integrate`Elliptic`Private`a + Integrate`Elliptic`Private`b)/2, 90], "
        "w (N[Integrate`Elliptic`Private`f /. Integrate`Elliptic`Private`x -> Integrate`Elliptic`Private`a + "
        "(Integrate`Elliptic`Private`b - Integrate`Elliptic`Private`a) q, 90] + "
        "N[Integrate`Elliptic`Private`f /. Integrate`Elliptic`Private`x -> Integrate`Elliptic`Private`b - "
        "(Integrate`Elliptic`Private`b - Integrate`Elliptic`Private`a) q, 90])]; "
        "tot += v; k++; If[q < 10^-70 || k > 400 || (k > 8 && Abs[v] < 10^-45 Abs[tot]), go = False]]; "
        "Integrate`Elliptic`Private`h tot]",
        "Integrate`Elliptic`Private`tsH[Integrate`Elliptic`Private`f_, Integrate`Elliptic`Private`x_, "
        "Integrate`Elliptic`Private`a_, Integrate`Elliptic`Private`sg_, Integrate`Elliptic`Private`h_] := "
        "Module[{tot = 0, k, t, e, v, go, dir}, Do[k = If[dir > 0, 0, -1]; go = True; While[go, "
        "t = k Integrate`Elliptic`Private`h; e = Exp[N[Pi/2 Sinh[t], 90]]; "
        "v = N[Pi/2 Cosh[t], 90] e N[Integrate`Elliptic`Private`f /. Integrate`Elliptic`Private`x -> "
        "Integrate`Elliptic`Private`a + Integrate`Elliptic`Private`sg e, 90]; "
        "tot += v; k += dir; If[e < 10^-70 || e > 10^200 || Abs[k] > 400 || "
        "(Abs[k] > 8 && Abs[v] < 10^-45 Abs[tot]), go = False]], {dir, {1, -1}}]; "
        "Integrate`Elliptic`Private`h tot]"
    };
    for (size_t i = 0; i < sizeof ts_defs / sizeof ts_defs[0]; i++) {
        Expr* d = parse_expression(ts_defs[i]);
        if (d) { Expr* r = evaluate(d); expr_free(r); expr_free(d); }
    }
    symtab_add_builtin("Integrate`Elliptic", builtin_integrate_elliptic);
    symtab_get_def("Integrate`Elliptic")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`Elliptic",
        "Integrate`Elliptic[f, {x, a, b}] evaluates the definite integral of "
        "R(x, Sqrt[P(x)]) over a real interval (limits may be infinite), with "
        "R rational and P a squarefree cubic or quartic whose radicands may be "
        "rational functions sharing one squarefree kernel, in closed form with "
        "EllipticF, EllipticE, EllipticPi, EllipticK and elementary terms "
        "(Carlson's reduction, converted to Legendre form).  The result is "
        "checked numerically; returns unevaluated when it cannot evaluate.  "
        "Accepts Assumptions -> ....");
}
