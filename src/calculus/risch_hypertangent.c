/* risch_hypertangent.c — the hypertangent case (Bronstein §5.10).
 *
 * See risch_hypertangent.h.  IntegrateHypertangentPolynomial reduces p modulo
 * derivatives (PolynomialReduce, since a hypertangent is a nonlinear monomial)
 * to a remainder of degree <= 1, then reads off the log-term coefficient
 * c = coeff(r, t) / (2a), where a = Dt/(t^2+1).
 */

#include "risch_hypertangent.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "risch_field.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Evaluation helpers.                                                 */
/* ------------------------------------------------------------------ */

static Expr* rh_cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* rh_eval_adopt(Expr* call) { Expr* r = evaluate(call); expr_free(call); return r; }
static Expr* rh_fn(const char* head, Expr** args, size_t n) {
    return rh_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rh_call1(const char* head, Expr* a) { return rh_fn(head, (Expr*[]){ a }, 1); }
static Expr* rh_call2(const char* head, Expr* a, Expr* b) { return rh_fn(head, (Expr*[]){ a, b }, 2); }
static Expr* rh_call3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rh_fn(head, (Expr*[]){ a, b, c }, 3);
}
static Expr* rh_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rh_plus(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* rh_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ a, expr_new_integer(n) }, 2);
}
/* t^2 + 1 (evaluated). */
static Expr* rh_t2p1(const Expr* t) {
    return rh_eval_adopt(rh_plus(rh_pow(rh_cp(t), 2), expr_new_integer(1)));
}
/* Dt/(t^2+1) reduced; true (and *eta set, owned) iff it lies in k (free of t). */
static bool rh_eta_in_k(const Expr* Dt, const Expr* t, Expr** eta_out) {
    Expr* t2p1 = rh_t2p1(t);
    Expr* eta = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(t2p1, -1)));  /* adopts t2p1 */
    Expr* fq = rh_fn("FreeQ", (Expr*[]){ rh_cp(eta), rh_cp(t) }, 2);
    bool ok = fq && fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("True");
    expr_free(fq);
    if (!ok) { expr_free(eta); return false; }
    if (eta_out) *eta_out = eta; else expr_free(eta);
    return true;
}
/* True iff Together[e] is the integer 0. */
static bool rh_is_zero(const Expr* e) {
    Expr* z = rh_call1("Together", rh_cp(e));
    bool r = z && z->type == EXPR_INTEGER && z->data.integer == 0;
    expr_free(z);
    return r;
}
/* True iff e is the symbol True. */
static bool rh_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("True");
}
/* Copy the n arguments of a List[...] of arity n into out[0..n); false otherwise. */
static bool rh_list_parts(const Expr* e, size_t n, Expr** out) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != n) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol.name != intern_symbol("List")) return false;
    for (size_t i = 0; i < n; i++) out[i] = expr_copy(e->data.function.args[i]);
    return true;
}

/* ------------------------------------------------------------------ */
/* IntegrateHypertangentPolynomial (§5.10, p.167).                     */
/* ------------------------------------------------------------------ */

bool risch_integrate_hypertangent_poly(const Expr* p, const Expr* t,
                                       const RischDeriv* d, Expr** qo, Expr** co) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    /* a = Dt/(t^2+1) must lie in k (t is hypertangent). */
    Expr* tsq1 = rh_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ rh_pow(rh_cp(t), 2), expr_new_integer(1) }, 2));
    Expr* a = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(rh_cp(tsq1), -1)));
    expr_free(tsq1);
    /* a must be free of t (element of the base field). */
    Expr* fq = rh_fn("FreeQ", (Expr*[]){ rh_cp(a), rh_cp(t) }, 2);
    bool a_in_k = fq && fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("True");
    expr_free(fq);
    if (!a_in_k) { expr_free(a); return false; }

    Expr* q = NULL; Expr* r = NULL;
    if (!risch_field_polynomial_reduce(p, t, d, &q, &r)) { expr_free(a); return false; }

    /* c = coeff(r, t) / (2 a).  (deg_t(r) <= 1 after the reduction.) */
    Expr* r1 = rh_call3("Coefficient", r, rh_cp(t), expr_new_integer(1));  /* adopts r */
    Expr* twoa = rh_eval_adopt(rh_times(expr_new_integer(2), a));          /* adopts a */
    Expr* c = rh_call1("Cancel", rh_times(r1, rh_pow(twoa, -1)));

    *qo = q;
    *co = c;
    return true;
}

/* ------------------------------------------------------------------ */
/* ResidueReduce (Bronstein §5.6, the residue criterion).              */
/* ------------------------------------------------------------------ */

/* For the simple part h in k(t) (proper, squarefree normal denominator in t),
 * compute g2 = sum_j c_j Log(g_j) and beta in {0,1} with h - D[g2] in k[t] when
 * every Rothstein-Trager residue c_j is a constant of the derivation (beta = 1),
 * or beta = 0 (g2 = 0) when a residue is non-constant — in which case h has no
 * elementary integral over k(t) (Thm 5.6.1).  Built on the same LazardRioboo-
 * Trager engine the tower integrator uses (Integrate`TranscendentalLogPart),
 * with the derivation's lower variables {x, t_0, ..., t_{L-1}} as the constant-
 * residue gate.  Owned outputs. */
bool risch_residue_reduce(const Expr* h, const Expr* t, const RischDeriv* d,
                          Expr** g2o, bool* beta) {
    Expr* num = NULL; Expr* den = NULL;
    risch_field_num_den_t(h, t, &num, &den);
    if (risch_field_degree_t(den, t) <= 0) {   /* no simple poles: nothing to log */
        expr_free(num); expr_free(den);
        *g2o = expr_new_integer(0); *beta = true; return true;
    }
    /* Strip the t-free content:  Hden = gcd(den, d den/dt),  rad = den/Hden
     * (squarefree in t),  a = num/Hden. */
    Expr* dden = rh_call2("D", rh_cp(den), rh_cp(t));
    Expr* Hden = rh_call2("PolynomialGCD", rh_cp(den), dden);   /* adopts dden */
    Expr* rad = rh_call1("Cancel", rh_times(rh_cp(den), rh_pow(rh_cp(Hden), -1)));
    Expr* a   = rh_call1("Cancel", rh_times(rh_cp(num), rh_pow(rh_cp(Hden), -1)));
    expr_free(num); expr_free(den); expr_free(Hden);

    Expr* Dd = risch_field_deriv(rad, d);            /* monomial derivation of rad */

    /* Gate list: every derivation variable EXCEPT t (residues must be constants
     * of D, i.e. free of x and the lower monomials). */
    size_t ng = 0;
    for (size_t i = 0; i < d->nvars; i++) if (!expr_eq(d->vars[i], (Expr*)t)) ng++;
    Expr** gv = malloc((ng ? ng : 1) * sizeof(Expr*));
    size_t j = 0;
    for (size_t i = 0; i < d->nvars; i++)
        if (!expr_eq(d->vars[i], (Expr*)t)) gv[j++] = expr_copy(d->vars[i]);
    Expr* gate = expr_new_function(expr_new_symbol("List"), gv, ng);
    free(gv);

    Expr* logpart = rh_fn("Integrate`TranscendentalLogPart",
        (Expr*[]){ a, rad, rh_cp(t), expr_new_symbol("rhZ"), Dd, gate }, 6);
    /* Declined (head unchanged) => a residue is non-constant => beta = 0. */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol.name
             == intern_symbol("Integrate`TranscendentalLogPart")) {
        expr_free(logpart);
        *g2o = expr_new_integer(0); *beta = false; return true;
    }
    *g2o = logpart; *beta = true; return true;
}

/* ------------------------------------------------------------------ */
/* IntegrateHypertangent (Bronstein §5.10, p.172 — the full driver).   */
/* ------------------------------------------------------------------ */

/* 1 - t^2 (evaluated) — the REAL special polynomial for the hyperbolic case
 * (1 - Tanh^2 = Sech^2 > 0), whose log stays real. */
static Expr* rh_1mt2(const Expr* t) {
    return rh_eval_adopt(rh_plus(expr_new_integer(1),
        rh_times(expr_new_integer(-1), rh_pow(rh_cp(t), 2))));
}

/* Shared driver core for the quadratic-special cases (Bronstein §5.10 for the
 * hypertangent t^2+1, and its hyperbolic analogue for 1-t^2).  Composes the four
 * sub-algorithms exactly as the book box does, dispatching the reduced and
 * polynomial steps to the named builtins and using the given special polynomial
 * for the final log term:
 *   (g1, h, r) <- HermiteReduce(f)
 *   (g2, b)    <- ResidueReduce(h);              if b = 0 return (g1 + g2, 0)
 *   p = h - D[g2] + r
 *   (q1, b)    <- <reduced_name>(p);             if b = 0 return (g1+g2+q1, 0)
 *   (q2, c)    <- <poly_name>(p - D[q1])
 *   return (g1+g2+q1+q2 [+ c Log(special) if Dc = 0], Dc = 0). */

/* Normalize a reconciliation remainder (h - D[g2] + r, or p - D[q1]) to its
 * true polynomial-in-t form.  For an irreducible normal pole with real
 * algebraic residues (disc > 0), ResidueReduce emits the residues as conjugate
 * radical pairs a +/- b Sqrt[d]; their contributions to the remainder cancel to
 * a polynomial over Q, but Together treats each Sqrt[d] as opaque and leaves
 * spurious t-poles, so the remainder is not recognised as a polynomial and the
 * downstream polynomial driver spins forever.  RootReduce FIRST reduces the
 * algebraic coefficients to canonical form so the conjugates combine; Cancel
 * [Together] then clears the (now genuinely cancelling) denominators.  On the
 * common rational-residue case (split-linear denominator) there are no radicals
 * and RootReduce is a fast no-op. */
static Expr* rh_polynorm(Expr* e) {   /* adopts e */
    return rh_call1("Cancel", rh_call1("Together", rh_call1("RootReduce", e)));
}

static bool rh_integrate_core(const Expr* f, const Expr* t, const Expr* deriv,
                              const char* reduced_name, const char* poly_name,
                              bool special_plus, Expr** go, bool* beta) {
    Expr* parts[3];
    Expr* H = rh_call3("Risch`HermiteReduce", rh_cp(f), rh_cp(t), rh_cp(deriv));
    bool ok = rh_list_parts(H, 3, parts);
    expr_free(H);
    if (!ok) return false;
    Expr* g1 = parts[0], *h = parts[1], *r = parts[2];

    Expr* rr[2];
    Expr* RR = rh_call3("Risch`ResidueReduce", rh_cp(h), rh_cp(t), rh_cp(deriv));
    ok = rh_list_parts(RR, 2, rr);
    expr_free(RR);
    if (!ok) { expr_free(g1); expr_free(h); expr_free(r); return false; }
    Expr* g2 = rr[0]; bool b1 = rh_is_true(rr[1]); expr_free(rr[1]);

    if (!b1) {                                   /* h non-elementary (Thm 5.6.1) */
        *go = rh_call1("Cancel", rh_plus(g1, g2));   /* adopts g1, g2 */
        expr_free(h); expr_free(r);
        *beta = false; return true;
    }

    /* p = polynomial part = Cancel[Together[RootReduce[h - D[g2] + r]]]
     * (RootReduce so real-algebraic conjugate residues combine — see rh_polynorm). */
    Expr* Dg2 = rh_call2("Risch`Derivation", rh_cp(g2), rh_cp(deriv));
    Expr* p = rh_polynorm(
        rh_plus(rh_plus(h, rh_times(expr_new_integer(-1), Dg2)), r));  /* adopts h, r, Dg2 */

    Expr* ir[2];
    Expr* IR = rh_call3(reduced_name, rh_cp(p), rh_cp(t), rh_cp(deriv));
    ok = rh_list_parts(IR, 2, ir);
    expr_free(IR);
    if (!ok) { expr_free(g1); expr_free(g2); expr_free(p); return false; }
    Expr* q1 = ir[0]; bool b2 = rh_is_true(ir[1]); expr_free(ir[1]);

    if (!b2) {                                   /* reduced part non-elementary */
        *go = rh_call1("Cancel", rh_plus(rh_plus(g1, g2), q1));   /* adopts g1,g2,q1 */
        expr_free(p);
        *beta = false; return true;
    }

    /* p2 = Cancel[Together[RootReduce[p - D[q1]]]] (same algebraic-residue care). */
    Expr* Dq1 = rh_call2("Risch`Derivation", rh_cp(q1), rh_cp(deriv));
    Expr* p2 = rh_polynorm(
        rh_plus(p, rh_times(expr_new_integer(-1), Dq1)));         /* adopts p, Dq1 */

    Expr* ip[2];
    Expr* IP = rh_call3(poly_name, rh_cp(p2), rh_cp(t), rh_cp(deriv));
    expr_free(p2);
    ok = rh_list_parts(IP, 2, ip);
    expr_free(IP);
    if (!ok) { expr_free(g1); expr_free(g2); expr_free(q1); return false; }
    Expr* q2 = ip[0]; Expr* c = ip[1];

    Expr* Dc = rh_call2("Risch`Derivation", rh_cp(c), rh_cp(deriv));
    bool dc_zero = rh_is_zero(Dc);
    expr_free(Dc);

    Expr* base = rh_plus(rh_plus(rh_plus(g1, g2), q1), q2);   /* adopts g1,g2,q1,q2 */
    if (dc_zero) {                               /* elementary: + c Log(special) */
        Expr* special = special_plus ? rh_t2p1(t) : rh_1mt2(t);
        Expr* logterm = rh_times(c, rh_call1("Log", special));  /* adopts c, special */
        *go = rh_call1("Cancel", rh_plus(base, logterm));
        *beta = true;
    } else {                                     /* Dc != 0: non-elementary */
        expr_free(c);
        *go = rh_call1("Cancel", base);
        *beta = false;
    }
    return true;
}

/* IntegrateHypertangent (Bronstein §5.10): the circular tangent case (t^2+1). */
bool risch_integrate_hypertangent(const Expr* f, const Expr* t,
                                  const Expr* deriv, Expr** go, bool* beta) {
    return rh_integrate_core(f, t, deriv,
        "Risch`IntegrateHypertangentReduced",
        "Risch`IntegrateHypertangentPolynomial", true, go, beta);
}

/* ================================================================== */
/* Hyperbolic tangent case:  Dt = eta (t^2 - 1),  t = Tanh(∫eta).       */
/* ================================================================== */
/*
 * The special polynomial is t^2 - 1 = (t-1)(t+1), which SPLITS over k.  As a
 * result the reduced-case pole peel — a coupled system [[0,-2m eta],[-2m eta,0]]
 * (c;d) = (a;b) — DECOUPLES over the reals: with P = c+d, Q = c-d it becomes two
 * independent base Risch differential equations  DP - 2m eta P = a+b  and
 * DQ + 2m eta Q = a-b, so no C(i) is needed (unlike the circular tangent case).
 */

/* t^2 - 1 (evaluated). */
static Expr* rh_t2m1(const Expr* t) {
    return rh_eval_adopt(rh_plus(rh_pow(rh_cp(t), 2), expr_new_integer(-1)));
}
/* FreeQ[e, v]. */
static bool rh_free_of(const Expr* e, const Expr* v) {
    Expr* fq = rh_fn("FreeQ", (Expr*[]){ rh_cp(e), rh_cp(v) }, 2);
    bool r = fq && fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("True");
    expr_free(fq);
    return r;
}
/* The base integration variable (derivative 1) of the derivation. */
static const Expr* rh_base_var(const RischDeriv* d) {
    for (size_t i = 0; i < d->nvars; i++) {
        const Expr* dv = d->dvars[i];
        if (dv && dv->type == EXPR_INTEGER && dv->data.integer == 1) return d->vars[i];
    }
    return NULL;
}
/* Solve D[y] + f y = g over C(x) via Risch`RischDE; NULL when no rational sol. */
static Expr* rh_rischde(const Expr* f, const Expr* g, const Expr* x) {
    Expr* y = rh_call3("Risch`RischDE", rh_cp(f), rh_cp(g), rh_cp(x));
    if (y && y->type == EXPR_FUNCTION && y->data.function.head->type == EXPR_SYMBOL
        && y->data.function.head->data.symbol.name == intern_symbol("Risch`RischDE")) {
        expr_free(y);
        return NULL;
    }
    return y;
}
/* Hyperbolic (a = +1) coupled system  Dc - f2 d = g1,  Dd - f2 c = g2  over k:
 * decouples as P=c+d: DP - f2 P = g1+g2,  Q=c-d: DQ + f2 Q = g1-g2, then
 * c=(P+Q)/2, d=(P-Q)/2.  Returns false if either Risch DE has no solution. */
static bool rh_hyperbolic_coupled(const Expr* f2, const Expr* g1, const Expr* g2,
                                  const Expr* x, Expr** co, Expr** ddo) {
    Expr* negf2 = rh_call1("Cancel", rh_times(expr_new_integer(-1), rh_cp(f2)));
    Expr* sum = rh_call1("Cancel", rh_plus(rh_cp(g1), rh_cp(g2)));
    Expr* dif = rh_call1("Cancel", rh_plus(rh_cp(g1), rh_times(expr_new_integer(-1), rh_cp(g2))));
    Expr* P = rh_rischde(negf2, sum, x);          /* Dy - f2 y = g1 + g2 */
    expr_free(negf2); expr_free(sum);
    if (!P) { expr_free(dif); return false; }
    Expr* Q = rh_rischde(f2, dif, x);             /* Dy + f2 y = g1 - g2 */
    expr_free(dif);
    if (!Q) { expr_free(P); return false; }
    *co  = rh_call1("Cancel", rh_times(rh_plus(rh_cp(P), rh_cp(Q)), rh_pow(expr_new_integer(2), -1)));
    *ddo = rh_call1("Cancel", rh_times(rh_plus(rh_cp(P), rh_times(expr_new_integer(-1), rh_cp(Q))),
                                       rh_pow(expr_new_integer(2), -1)));
    expr_free(P); expr_free(Q);
    return true;
}
/* Multiplicity of (t^2-1) in a polynomial over k[t]. */
static long rh_mult_t2m1(const Expr* poly, const Expr* t) {
    Expr* s = rh_t2m1(t);
    Expr* cur = rh_cp(poly);
    long m = 0;
    for (;;) {
        if (rh_is_zero(cur)) break;
        Expr* q = risch_field_divexact_t(cur, s, t);
        if (!q) break;
        expr_free(cur); cur = q; m++;
        if (m > 4096) break;
    }
    expr_free(cur); expr_free(s);
    return m;
}

/* IntegrateHypertanhReduced: for a hyperbolic-tangent monomial t (Dt = eta(t^2-1))
 * and p in k(t) whose only poles are at t^2-1, peel the poles one multiplicity at
 * a time via the decoupled hyperbolic coupled system, writing q in k(t) with
 * p - D[q] in k[t] (*beta true), or *beta false when p has no elementary integral
 * over k(t).  Owned q. */
static bool rh_integrate_hypertanh_reduced(const Expr* p, const Expr* t,
                                           const RischDeriv* d, Expr** qo, bool* beta) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    if (rh_is_zero(p)) { *qo = expr_new_integer(0); *beta = true; return true; }
    Expr* s = rh_t2m1(t);
    Expr* eta = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(rh_cp(s), -1)));  /* Dt/(t^2-1) */
    if (!rh_free_of(eta, t)) { expr_free(eta); expr_free(s); return false; }

    Expr* num = NULL; Expr* den = NULL;
    risch_field_num_den_t(p, t, &num, &den);
    long m = rh_mult_t2m1(den, t) - rh_mult_t2m1(num, t);
    expr_free(num); expr_free(den);
    if (m <= 0) { expr_free(eta); expr_free(s); *qo = expr_new_integer(0); *beta = true; return true; }

    Expr* h = rh_call1("Cancel", rh_call1("Together", rh_times(rh_pow(rh_cp(s), m), rh_cp(p))));
    Expr* hnum = NULL; Expr* hden = NULL;
    risch_field_num_den_t(h, t, &hnum, &hden);
    bool h_is_poly = (risch_field_degree_t(hden, t) == 0);
    expr_free(hnum); expr_free(hden);
    if (!h_is_poly) { expr_free(h); expr_free(eta); expr_free(s); return false; }

    Expr* qd = NULL; Expr* r = NULL;
    if (!risch_field_divmod_t(h, s, t, &qd, &r)) {
        expr_free(h); expr_free(eta); expr_free(s); return false;
    }
    expr_free(h); expr_free(qd);
    Expr* a = rh_call3("Coefficient", rh_cp(r), rh_cp(t), expr_new_integer(1));
    Expr* b = rh_call3("Coefficient", r, rh_cp(t), expr_new_integer(0));  /* adopts r */

    Expr* f2 = rh_call1("Cancel", rh_times(expr_new_integer(2 * m), rh_cp(eta)));
    const Expr* x = rh_base_var(d);
    if (!x) { expr_free(f2); expr_free(a); expr_free(b); expr_free(eta); expr_free(s); return false; }

    Expr* c = NULL; Expr* dd = NULL;
    bool solved = rh_hyperbolic_coupled(f2, a, b, x, &c, &dd);
    expr_free(f2); expr_free(a); expr_free(b);
    if (!solved) { expr_free(eta); expr_free(s); *qo = expr_new_integer(0); *beta = false; return true; }

    /* q0 = (c t + d)/(t^2-1)^m. */
    Expr* ctd = rh_plus(rh_times(c, rh_cp(t)), dd);   /* adopts c, dd */
    Expr* q0 = rh_call1("Cancel", rh_call1("Together", rh_times(ctd, rh_pow(rh_cp(s), -m))));
    expr_free(eta); expr_free(s);

    Expr* Dq0 = risch_field_deriv(q0, d);
    Expr* pnext = rh_call1("Cancel", rh_call1("Together",
        rh_plus(rh_cp(p), rh_times(expr_new_integer(-1), Dq0))));
    Expr* qrest = NULL; bool brest = false;
    bool ok = rh_integrate_hypertanh_reduced(pnext, t, d, &qrest, &brest);
    expr_free(pnext);
    if (!ok) { expr_free(q0); return false; }
    *qo = rh_call1("Cancel", rh_call1("Together", rh_plus(qrest, q0)));  /* adopts qrest, q0 */
    *beta = brest;
    return true;
}

/* IntegrateHypertanhPolynomial: for p in k[t], write q in k[t] and c in k with
 * p - D[q] - c D(t^2-1)/(t^2-1) in k (the hyperbolic analogue of §5.10 p.167).
 * Since D(t^2-1)/(t^2-1) = 2 eta t with eta = Dt/(t^2-1), c = coeff(r,t)/(2 eta)
 * after PolynomialReduce. */
static bool rh_integrate_hypertanh_poly(const Expr* p, const Expr* t,
                                        const RischDeriv* d, Expr** qo, Expr** co) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    Expr* s = rh_t2m1(t);
    Expr* eta = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(s, -1)));  /* adopts s */
    if (!rh_free_of(eta, t)) { expr_free(eta); return false; }
    Expr* q = NULL; Expr* r = NULL;
    if (!risch_field_polynomial_reduce(p, t, d, &q, &r)) { expr_free(eta); return false; }
    Expr* r1 = rh_call3("Coefficient", r, rh_cp(t), expr_new_integer(1));  /* adopts r */
    Expr* twoeta = rh_eval_adopt(rh_times(expr_new_integer(2), eta));      /* adopts eta */
    Expr* c = rh_call1("Cancel", rh_times(r1, rh_pow(twoeta, -1)));
    *qo = q; *co = c;
    return true;
}

/* IntegrateHypertanh: the full driver for the hyperbolic tangent case (special
 * 1-t^2 for the real log term).  Same composition as IntegrateHypertangent. */
bool risch_integrate_hypertanh(const Expr* f, const Expr* t,
                               const Expr* deriv, Expr** go, bool* beta) {
    return rh_integrate_core(f, t, deriv,
        "Risch`IntegrateHypertanhReduced",
        "Risch`IntegrateHypertanhPolynomial", false, go, beta);
}

/* ------------------------------------------------------------------ */
/* Builtin.                                                            */
/* ------------------------------------------------------------------ */

/* Risch`IntegrateHypertangentPolynomial[p, t, deriv] -> {q, c}. */
static Expr* builtin_risch_int_hypertan_poly(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* q = NULL; Expr* c = NULL;
    bool ok = risch_integrate_hypertangent_poly(p, t, &d, &q, &c);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q, c }, 2);
}

/* Risch`ResidueReduce[h, t, deriv] -> {g2, beta}. */
static Expr* builtin_risch_residue_reduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* h = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* g2 = NULL; bool beta = false;
    bool ok = risch_residue_reduce(h, t, &d, &g2, &beta);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ g2, expr_new_symbol(beta ? "True" : "False") }, 2);
}

/* Risch`IntegrateHypertangent[f, t, deriv] -> {g, beta}. */
static Expr* builtin_risch_integrate_hypertangent(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* f = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    /* Fail fast unless t is a hypertangent monomial (Dt/(t^2+1) in k). */
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    const Expr* Dt = risch_deriv_lookup(&d, t);
    bool hyper = Dt && rh_eta_in_k(Dt, t, NULL);
    risch_deriv_free(&d);
    if (!hyper) return NULL;

    Expr* g = NULL; bool beta = false;
    bool ok = risch_integrate_hypertangent(f, t, res->data.function.args[2], &g, &beta);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ g, expr_new_symbol(beta ? "True" : "False") }, 2);
}

/* eta = Dt/(t^2-1); true iff it lies in k (t is a hyperbolic tangent monomial). */
static bool rh_eta_in_k_hyp(const Expr* Dt, const Expr* t) {
    Expr* s = rh_t2m1(t);
    Expr* eta = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(s, -1)));  /* adopts s */
    bool ok = rh_free_of(eta, t);
    expr_free(eta);
    return ok;
}

/* Risch`IntegrateHypertanhReduced[p, t, deriv] -> {q, beta}. */
static Expr* builtin_risch_int_hypertanh_reduced(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* q = NULL; bool beta = false;
    bool ok = rh_integrate_hypertanh_reduced(p, t, &d, &q, &beta);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ q, expr_new_symbol(beta ? "True" : "False") }, 2);
}

/* Risch`IntegrateHypertanhPolynomial[p, t, deriv] -> {q, c}. */
static Expr* builtin_risch_int_hypertanh_poly(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* q = NULL; Expr* c = NULL;
    bool ok = rh_integrate_hypertanh_poly(p, t, &d, &q, &c);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q, c }, 2);
}

/* Risch`IntegrateHypertanh[f, t, deriv] -> {g, beta}. */
static Expr* builtin_risch_integrate_hypertanh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* f = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    /* Fail fast unless t is a hyperbolic tangent monomial (Dt/(t^2-1) in k). */
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    const Expr* Dt = risch_deriv_lookup(&d, t);
    bool hyper = Dt && rh_eta_in_k_hyp(Dt, t);
    risch_deriv_free(&d);
    if (!hyper) return NULL;
    Expr* g = NULL; bool beta = false;
    bool ok = risch_integrate_hypertanh(f, t, res->data.function.args[2], &g, &beta);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ g, expr_new_symbol(beta ? "True" : "False") }, 2);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rh_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_hypertangent_init(void) {
    rh_install("Risch`IntegrateHypertangentPolynomial", builtin_risch_int_hypertan_poly,
        "Risch`IntegrateHypertangentPolynomial[p, t, deriv] integrates the\n"
        "polynomial p over a hypertangent monomial t (Dt = a (t^2+1)): returns\n"
        "{q, c} with p - D[q] - c D[t^2+1]/(t^2+1) in the base field, so that the\n"
        "integral is q + c Log[t^2+1] plus a base-field integral (Bronstein 5.10).\n"
        "If D[c] is nonzero the integral is not elementary.");
    rh_install("Risch`ResidueReduce", builtin_risch_residue_reduce,
        "Risch`ResidueReduce[h, t, deriv] applies the residue criterion (Bronstein\n"
        "5.6) to the simple part h in k(t): returns {g2, beta} with g2 a sum of\n"
        "logarithms and h - D[g2] in k[t] when beta is True (every Rothstein-Trager\n"
        "residue is a constant of the derivation), or {0, False} when a residue is\n"
        "non-constant, in which case h has no elementary integral over k(t).");
    rh_install("Risch`IntegrateHypertangent", builtin_risch_integrate_hypertangent,
        "Risch`IntegrateHypertangent[f, t, deriv] integrates f in k(t) over a\n"
        "hypertangent monomial t (Dt = a (t^2+1)) directly, keeping tan/arctan real\n"
        "without a complex-exponential rewrite (Bronstein 5.10): returns {g, beta}\n"
        "with f - D[g] in k when beta is True, or f - D[g] having no elementary\n"
        "integral over k(t) when beta is False.  Composes Hermite reduction, the\n"
        "residue criterion, and the reduced and polynomial hypertangent cases.");
    rh_install("Risch`IntegrateHypertanhReduced", builtin_risch_int_hypertanh_reduced,
        "Risch`IntegrateHypertanhReduced[p, t, deriv] peels the t^2-1 poles of a\n"
        "reduced p over a hyperbolic tangent monomial t (Dt = a (t^2-1)) via the\n"
        "decoupled hyperbolic coupled system, returning {q, beta} with p - D[q] in\n"
        "k[t] when beta is True.");
    rh_install("Risch`IntegrateHypertanhPolynomial", builtin_risch_int_hypertanh_poly,
        "Risch`IntegrateHypertanhPolynomial[p, t, deriv] integrates the polynomial p\n"
        "over a hyperbolic tangent monomial t (Dt = a (t^2-1)): returns {q, c} with\n"
        "p - D[q] - c D[t^2-1]/(t^2-1) in the base field.  If D[c] is nonzero the\n"
        "integral is not elementary.");
    rh_install("Risch`IntegrateHypertanh", builtin_risch_integrate_hypertanh,
        "Risch`IntegrateHypertanh[f, t, deriv] integrates f in k(t) over a hyperbolic\n"
        "tangent monomial t (Dt = a (t^2-1)) directly, keeping tanh/artanh real:\n"
        "returns {g, beta} with f - D[g] in k when beta is True.  The special\n"
        "polynomial t^2-1 splits over k, so the reduced case decouples into two real\n"
        "base Risch differential equations (no C(i)).");
}
