/* risch_coupled.c — coupled differential systems + hypertangent reduced case.
 *
 * See risch_coupled.h.  The coupled 2x2 real system (Bronstein Eq. 8.1, a = -1)
 * is solved by the §8.1 reduction to a single Risch DE over C(i)(x): form
 * f = f1 + f2 i, g = g1 + g2 i, solve D y + f y = g with the base-field solver
 * (Risch`RischDE), then split y = y1 + y2 i by algebraic conjugation
 * (y1 = (y + conj y)/2, y2 = (y - conj y)/(2 i)), where conj replaces every
 * Complex[re, im] by Complex[re, -im].  IntegrateHypertangentReduced (§5.10)
 * sits on top, peeling t^2+1 poles via the coupled system Eq. (5.20).
 */

#include "risch_coupled.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "risch_field.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Evaluation helpers (mirroring risch_hypertangent.c).                */
/* ------------------------------------------------------------------ */

static Expr* rc_cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* rc_eval_adopt(Expr* call) { Expr* r = evaluate(call); expr_free(call); return r; }
static Expr* rc_fn(const char* head, Expr** args, size_t n) {
    return rc_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rc_call1(const char* head, Expr* a) { return rc_fn(head, (Expr*[]){ a }, 1); }
static Expr* rc_call2(const char* head, Expr* a, Expr* b) { return rc_fn(head, (Expr*[]){ a, b }, 2); }
static Expr* rc_call3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rc_fn(head, (Expr*[]){ a, b, c }, 3);
}
/* Unevaluated builders (no eval): raw Plus/Times/Power trees. */
static Expr* rc_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rc_plus(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* rc_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ a, expr_new_integer(n) }, 2);
}
/* Normalize a rational function in the tower variables to a single reduced
 * fraction.  Cancel alone does NOT combine a sum of fractions (e.g. the output
 * of the monomial derivation) over a common denominator; Together does, and
 * also cancels — essential for the pole-order valuation in the reduced case to
 * be computed on a genuinely reduced denominator.  adopts a. */
static Expr* rc_cancel(Expr* a) { return rc_call1("Cancel", rc_call1("Together", a)); }

/* The imaginary unit as an evaluated literal Complex[0,1]. */
static Expr* rc_I(void) {
    return rc_fn("Complex", (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
}

/* t^2 + 1 as an (unevaluated) Expr in the monomial variable t. */
static Expr* rc_t2p1(const Expr* t) {
    return rc_plus(rc_pow(rc_cp(t), 2), expr_new_integer(1));
}

/* True iff e is free of the variable v (FreeQ[e, v] == True). */
static bool rc_free_of(const Expr* e, const Expr* v) {
    Expr* fq = rc_fn("FreeQ", (Expr*[]){ rc_cp(e), rc_cp(v) }, 2);
    bool r = fq && fq->type == EXPR_SYMBOL && fq->data.symbol.name == intern_symbol("True");
    expr_free(fq);
    return r;
}

/* True iff e is (rationally) the zero element:  Together[e] == 0. */
static bool rc_is_zero(const Expr* e) {
    Expr* z = rc_call1("Together", rc_cp(e));
    bool r = z && z->type == EXPR_INTEGER && z->data.integer == 0;
    expr_free(z);
    return r;
}

/* ------------------------------------------------------------------ */
/* Base-field Risch DE bridge and Gaussian conjugation.                */
/* ------------------------------------------------------------------ */

/* Solve D[y] + f y = g over C(x) (possibly Gaussian) via the Risch`RischDE
 * builtin.  Returns owned y, or NULL when there is no rational solution (the
 * builtin bubbles back unevaluated). */
static Expr* rc_rischde(const Expr* f, const Expr* g, const Expr* x) {
    Expr* y = rc_call3("Risch`RischDE", rc_cp(f), rc_cp(g), rc_cp(x));
    if (y && y->type == EXPR_FUNCTION && y->data.function.head->type == EXPR_SYMBOL
        && y->data.function.head->data.symbol.name == intern_symbol("Risch`RischDE")) {
        expr_free(y);
        return NULL;                               /* unevaluated => no solution */
    }
    return y;
}

/* Algebraic conjugation: replace every Complex[re, im] in y by Complex[re, -im]
 * (i.e. i -> -i), leaving all real/symbolic parts untouched.  For y in C(i)(x)
 * this returns the Gaussian conjugate.  Owned. */
static Expr* rc_conjugate(const Expr* y) {
    /* Rule:  Complex[a_, b_] :> a + (-1) b I  */
    Expr* pat_a = expr_new_function(expr_new_symbol("Pattern"),
        (Expr*[]){ expr_new_symbol("rc$a"),
                   expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2);
    Expr* pat_b = expr_new_function(expr_new_symbol("Pattern"),
        (Expr*[]){ expr_new_symbol("rc$b"),
                   expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2);
    Expr* lhs = expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ pat_a, pat_b }, 2);
    Expr* rhs = rc_plus(expr_new_symbol("rc$a"),
        expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_new_symbol("rc$b"), rc_I() }, 3));
    Expr* rule = expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){ lhs, rhs }, 2);
    return rc_call2("ReplaceAll", rc_cp(y), rule);   /* adopts rule */
}

/* Split a Gaussian element w = re + im*i (re, im free of the imaginary unit)
 * into its real and imaginary parts via algebraic conjugation:
 *   re = (w + conj w)/2,  im = (w - conj w)/(2 i).
 * Works whether re, im lie in k (constants in t) or in k[t] (w a polynomial in
 * t with Gaussian coefficients — conjugation leaves the real symbol t alone).
 * Both outputs owned. */
static void rc_split_gaussian(const Expr* w, Expr** re, Expr** im) {
    Expr* wc = rc_conjugate(w);
    Expr* sum = rc_plus(rc_cp(w), rc_cp(wc));
    *re = rc_cancel(rc_times(sum, rc_pow(expr_new_integer(2), -1)));
    Expr* diff = rc_plus(rc_cp(w), rc_times(expr_new_integer(-1), wc));  /* adopts wc */
    Expr* twoi = rc_eval_adopt(rc_times(expr_new_integer(2), rc_I()));
    *im = rc_cancel(rc_times(diff, rc_pow(twoi, -1)));                   /* adopts twoi */
}

/* The base integration variable of a derivation d: the one with derivative 1.
 * Borrowed pointer into d, or NULL if none is listed. */
static const Expr* rc_base_var(const RischDeriv* d) {
    for (size_t i = 0; i < d->nvars; i++) {
        const Expr* dv = d->dvars[i];
        if (dv && dv->type == EXPR_INTEGER && dv->data.integer == 1) return d->vars[i];
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* CoupledDESystem (a = -1), base field C(x)  (Bronstein §8.1).         */
/* ------------------------------------------------------------------ */

bool risch_coupled_desystem(const Expr* f1, const Expr* f2,
                            const Expr* g1, const Expr* g2, const Expr* x,
                            Expr** y1o, Expr** y2o) {
    /* f = f1 + f2 i,  g = g1 + g2 i. */
    Expr* f = rc_cancel(rc_plus(rc_cp(f1), rc_times(rc_cp(f2), rc_I())));
    Expr* g = rc_cancel(rc_plus(rc_cp(g1), rc_times(rc_cp(g2), rc_I())));
    Expr* y = rc_rischde(f, g, x);
    expr_free(f); expr_free(g);
    if (!y) return false;                          /* no solution */

    /* y1 = Re y, y2 = Im y (over k = C(x)). */
    rc_split_gaussian(y, y1o, y2o);
    expr_free(y);
    return true;
}

/* ------------------------------------------------------------------ */
/* IntegrateHypertangentReduced (Bronstein §5.10, p.169).              */
/* ------------------------------------------------------------------ */

/* Multiplicity of (t^2+1) in the polynomial poly over k[t]. */
static long rc_mult_t2p1(const Expr* poly, const Expr* t) {
    Expr* t2p1 = rc_eval_adopt(rc_t2p1(t));
    Expr* cur = rc_cp(poly);
    long m = 0;
    for (;;) {
        if (rc_is_zero(cur)) break;                /* 0 is divisible forever */
        Expr* q = risch_field_divexact_t(cur, t2p1, t);
        if (!q) break;
        expr_free(cur); cur = q; m++;              /* q now owned as cur */
        if (m > 4096) break;                       /* runaway guard */
    }
    expr_free(cur); expr_free(t2p1);
    return m;
}

bool risch_integrate_hypertangent_reduced(const Expr* p, const Expr* t,
                                          const RischDeriv* d, Expr** qo, bool* beta) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    /* p == 0 is already a (trivial) polynomial in k[t]: nothing to peel. */
    if (rc_is_zero(p)) { *qo = expr_new_integer(0); *beta = true; return true; }
    /* eta = Dt/(t^2+1) must lie in k (t is a hypertangent monomial). */
    Expr* t2p1 = rc_eval_adopt(rc_t2p1(t));
    Expr* eta = rc_cancel(rc_times(rc_cp(Dt), rc_pow(rc_cp(t2p1), -1)));
    if (!rc_free_of(eta, t)) { expr_free(eta); expr_free(t2p1); return false; }

    /* m = -nu_{t^2+1}(p) = mult in denominator - mult in numerator. */
    Expr* num = NULL; Expr* den = NULL;
    risch_field_num_den_t(p, t, &num, &den);
    long m = rc_mult_t2p1(den, t) - rc_mult_t2p1(num, t);
    expr_free(num); expr_free(den);
    if (m <= 0) {                                  /* no pole at t^2+1: done */
        expr_free(eta); expr_free(t2p1);
        *qo = expr_new_integer(0); *beta = true;
        return true;
    }

    /* h = (t^2+1)^m p  must be a polynomial in k[t] (p is reduced). */
    Expr* h = rc_cancel(rc_times(rc_pow(rc_cp(t2p1), m), rc_cp(p)));
    Expr* hnum = NULL; Expr* hden = NULL;
    risch_field_num_den_t(h, t, &hnum, &hden);
    bool h_is_poly = (risch_field_degree_t(hden, t) == 0);
    expr_free(hnum); expr_free(hden);
    if (!h_is_poly) { expr_free(h); expr_free(eta); expr_free(t2p1); return false; }

    /* PolyDivide(h, t^2+1) -> (qd, r), deg_t(r) < 2, so r = a t + b. */
    Expr* qd = NULL; Expr* r = NULL;
    if (!risch_field_divmod_t(h, t2p1, t, &qd, &r)) {
        expr_free(h); expr_free(eta); expr_free(t2p1); return false;
    }
    expr_free(h); expr_free(qd);
    Expr* a = rc_call3("Coefficient", rc_cp(r), rc_cp(t), expr_new_integer(1));
    Expr* b = rc_call3("Coefficient", r, rc_cp(t), expr_new_integer(0));  /* adopts r */

    /* Coupled system (5.20): (Dc; Dd) + [[0, -2 m eta],[2 m eta, 0]](c; d) = (a; b).
     * In CoupledDESystem's [[f1,-f2],[f2,f1]] form: f1 = 0, f2 = 2 m eta. */
    Expr* f2 = rc_cancel(rc_times(expr_new_integer(2 * m), rc_cp(eta)));
    Expr* zero = expr_new_integer(0);
    /* The base variable of the derivation is the one with derivative 1. */
    const Expr* x = rc_base_var(d);
    if (!x) { expr_free(f2); expr_free(zero); expr_free(a); expr_free(b);
              expr_free(eta); expr_free(t2p1); return false; }

    Expr* c = NULL; Expr* dd = NULL;
    bool solved = risch_coupled_desystem(zero, f2, a, b, x, &c, &dd);
    expr_free(zero); expr_free(f2); expr_free(a); expr_free(b);
    if (!solved) {                                 /* no solution: not elementary */
        expr_free(eta); expr_free(t2p1);
        *qo = expr_new_integer(0); *beta = false;
        return true;
    }

    /* q0 = (c t + d)/(t^2+1)^m. */
    Expr* ct = rc_times(c, rc_cp(t));              /* adopts c */
    Expr* ctd = rc_plus(ct, dd);                   /* adopts dd */
    Expr* q0 = rc_cancel(rc_times(ctd, rc_pow(rc_cp(t2p1), -m)));
    expr_free(eta); expr_free(t2p1);

    /* Recurse on p - D[q0] (pole order strictly drops). */
    Expr* Dq0 = risch_field_deriv(q0, d);
    Expr* pnext = rc_cancel(rc_plus(rc_cp(p), rc_times(expr_new_integer(-1), Dq0)));
    Expr* qrest = NULL; bool beta_rest = false;
    bool ok = risch_integrate_hypertangent_reduced(pnext, t, d, &qrest, &beta_rest);
    expr_free(pnext);
    if (!ok) { expr_free(q0); return false; }

    *qo = rc_cancel(rc_plus(qrest, q0));           /* adopts qrest, q0 */
    *beta = beta_rest;
    return true;
}

/* ------------------------------------------------------------------ */
/* CoupledDECancelTan (Bronstein §8.4, p.265).                          */
/* ------------------------------------------------------------------ */

/* Substitute the hypertangent monomial t -> sqrt(-1) into c in k[t], returning
 * c(sqrt(-1)) in k(sqrt(-1)) = k[i]  (powers t^{2j} collapse via i^2 = -1).
 * Owned. */
static Expr* rc_subst_t_i(const Expr* c, const Expr* t) {
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ rc_cp(t), rc_I() }, 2);
    Expr* r = rc_call2("ReplaceAll", rc_cp(c), rule);     /* adopts rule */
    return rc_cancel(r);                                   /* adopts r */
}

/* CoupledDECancelTan(b0, b2, c1, c2, D, n) — cancellation, tangent case.
 *
 * With t a hypertangent monomial of d (Dt/(t^2+1) = eta in k), b0,b2 in k and
 * c1,c2 in k[t], solves the tangent-cancellation coupled system
 *     (Dq1; Dq2) + [[b0 - n eta t, -b2], [b2, b0 - n eta t]] (q1; q2) = (c1; c2)
 * for q1,q2 in k[t] with deg_t <= n, keeping t REAL: the field only goes complex
 * one level down, inside the base CoupledDESystem over k = C(x).  x is the base
 * variable (Dx = 1), eta = Dt/(t^2+1); both are invariant across the recursion
 * and passed in precomputed.  Returns false for "no solution".  q1,q2 owned. */
static bool rc_cancel_tan(const Expr* b0, const Expr* b2,
                          const Expr* c1, const Expr* c2,
                          const Expr* t, const RischDeriv* d,
                          const Expr* x, const Expr* eta, long n,
                          Expr** q1o, Expr** q2o) {
    /* Base case n = 0: the system is [[b0,-b2],[b2,b0]](q1;q2)=(c1;c2) over k,
     * which is a genuine CoupledDESystem provided c1,c2 lie in k. */
    if (n == 0) {
        if (!rc_free_of(c1, t) || !rc_free_of(c2, t)) return false;
        return risch_coupled_desystem(b0, b2, c1, c2, x, q1o, q2o);
    }

    /* z1 + z2 i = c1(i) + c2(i) i, with z1, z2 in k. */
    Expr* c1i = rc_subst_t_i(c1, t);
    Expr* c2i = rc_subst_t_i(c2, t);
    Expr* w = rc_cancel(rc_plus(c1i, rc_times(c2i, rc_I())));   /* adopts c1i, c2i */
    Expr* z1 = NULL; Expr* z2 = NULL;
    rc_split_gaussian(w, &z1, &z2);
    expr_free(w);

    /* (s1, s2) <- CoupledDESystem(b0, b2 - n eta, z1, z2). */
    Expr* f2s = rc_cancel(rc_plus(rc_cp(b2), rc_times(expr_new_integer(-n), rc_cp(eta))));
    Expr* s1 = NULL; Expr* s2 = NULL;
    bool ok_s = risch_coupled_desystem(b0, f2s, z1, z2, x, &s1, &s2);
    expr_free(f2s);
    if (!ok_s) { expr_free(z1); expr_free(z2); return false; }

    /* Numerator of c:
     *   Nre = c1 - z1 + n eta (s1 t + s2)
     *   Nim = c2 - z2 + n eta (s2 t - s1)
     *   N   = Nre + Nim i     (an exact multiple of t - i, Bronstein p.264). */
    Expr* neta = rc_cancel(rc_times(expr_new_integer(n), rc_cp(eta)));
    Expr* s1t_s2 = rc_plus(rc_times(rc_cp(s1), rc_cp(t)), rc_cp(s2));
    Expr* Nre = rc_plus(rc_plus(rc_cp(c1), rc_times(expr_new_integer(-1), z1)),
                        rc_times(rc_cp(neta), s1t_s2));            /* adopts z1 */
    Expr* s2t_ms1 = rc_plus(rc_times(rc_cp(s2), rc_cp(t)),
                            rc_times(expr_new_integer(-1), rc_cp(s1)));
    Expr* Nim = rc_plus(rc_plus(rc_cp(c2), rc_times(expr_new_integer(-1), z2)),
                        rc_times(neta, s2t_ms1));                  /* adopts z2, neta */
    Expr* N = rc_cancel(rc_plus(Nre, rc_times(Nim, rc_I())));      /* adopts Nre, Nim */

    /* c = N / (t - i) in k(i)[t]; split into d1 + d2 i (d1, d2 in k[t]). */
    Expr* pI = rc_eval_adopt(rc_plus(rc_cp(t), rc_times(expr_new_integer(-1), rc_I())));
    Expr* Q = risch_field_divexact_t(N, pI, t);
    expr_free(N); expr_free(pI);
    if (!Q) { expr_free(s1); expr_free(s2); return false; }       /* not exact */
    Expr* d1 = NULL; Expr* d2 = NULL;
    rc_split_gaussian(Q, &d1, &d2);
    expr_free(Q);

    /* Recurse with b2 -> b2 + eta, n -> n - 1. */
    Expr* b2n = rc_cancel(rc_plus(rc_cp(b2), rc_cp(eta)));
    Expr* h1 = NULL; Expr* h2 = NULL;
    bool ok_h = rc_cancel_tan(b0, b2n, d1, d2, t, d, x, eta, n - 1, &h1, &h2);
    expr_free(b2n); expr_free(d1); expr_free(d2);
    if (!ok_h) { expr_free(s1); expr_free(s2); return false; }

    /* q1 = h1 t + h2 + s1,  q2 = h2 t - h1 + s2. */
    *q1o = rc_cancel(rc_plus(rc_plus(rc_times(rc_cp(h1), rc_cp(t)), rc_cp(h2)), rc_cp(s1)));
    *q2o = rc_cancel(rc_plus(rc_plus(rc_times(rc_cp(h2), rc_cp(t)),
                                     rc_times(expr_new_integer(-1), rc_cp(h1))), rc_cp(s2)));
    expr_free(h1); expr_free(h2); expr_free(s1); expr_free(s2);
    return true;
}

/* ------------------------------------------------------------------ */
/* Builtins.                                                           */
/* ------------------------------------------------------------------ */

/* Risch`CoupledDESystem[f1, f2, g1, g2, x] -> {y1, y2} | unevaluated. */
static Expr* builtin_risch_coupled_desystem(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 5) return NULL;
    const Expr* f1 = res->data.function.args[0];
    const Expr* f2 = res->data.function.args[1];
    const Expr* g1 = res->data.function.args[2];
    const Expr* g2 = res->data.function.args[3];
    const Expr* x  = res->data.function.args[4];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* y1 = NULL; Expr* y2 = NULL;
    if (!risch_coupled_desystem(f1, f2, g1, g2, x, &y1, &y2)) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ y1, y2 }, 2);
}

/* Risch`IntegrateHypertangentReduced[p, t, deriv] -> {q, beta}. */
static Expr* builtin_risch_int_hypertan_reduced(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* q = NULL; bool beta = false;
    bool ok = risch_integrate_hypertangent_reduced(p, t, &d, &q, &beta);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    Expr* bsym = expr_new_symbol(beta ? "True" : "False");
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q, bsym }, 2);
}

/* Risch`CoupledDECancelTan[b0, b2, c1, c2, t, deriv, n] -> {q1, q2} | uneval. */
static Expr* builtin_risch_coupled_cancel_tan(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 7) return NULL;
    Expr** a = res->data.function.args;
    const Expr* b0 = a[0];
    const Expr* b2 = a[1];
    const Expr* c1 = a[2];
    const Expr* c2 = a[3];
    const Expr* t  = a[4];
    const Expr* nn = a[6];
    if (t->type != EXPR_SYMBOL) return NULL;
    if (nn->type != EXPR_INTEGER || nn->data.integer < 0) return NULL;   /* n >= 0 */
    long n = (long) nn->data.integer;

    RischDeriv d;
    if (!risch_deriv_from_rules(a[5], &d)) return NULL;

    /* t must be a hypertangent monomial: eta = Dt/(t^2+1) must lie in k. */
    const Expr* Dt = risch_deriv_lookup(&d, t);
    if (!Dt) { risch_deriv_free(&d); return NULL; }
    Expr* t2p1 = rc_eval_adopt(rc_t2p1(t));
    Expr* eta = rc_cancel(rc_times(rc_cp(Dt), rc_pow(t2p1, -1)));         /* adopts t2p1 */
    if (!rc_free_of(eta, t)) { expr_free(eta); risch_deriv_free(&d); return NULL; }

    const Expr* x = rc_base_var(&d);
    if (!x) { expr_free(eta); risch_deriv_free(&d); return NULL; }

    Expr* q1 = NULL; Expr* q2 = NULL;
    bool ok = rc_cancel_tan(b0, b2, c1, c2, t, &d, x, eta, n, &q1, &q2);
    expr_free(eta); risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q1, q2 }, 2);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rc_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_coupled_init(void) {
    rc_install("Risch`CoupledDESystem", builtin_risch_coupled_desystem,
        "Risch`CoupledDESystem[f1, f2, g1, g2, x] solves the coupled 2x2 real\n"
        "differential system D[y1] + f1 y1 - f2 y2 == g1, D[y2] + f2 y1 + f1 y2\n"
        "== g2 for y1, y2 in C(x), via the Bronstein Chapter 8 reduction to a\n"
        "single Risch differential equation over C(i)(x).  Returns {y1, y2}, or\n"
        "stays unevaluated when the system has no rational solution.");
    rc_install("Risch`CoupledDECancelTan", builtin_risch_coupled_cancel_tan,
        "Risch`CoupledDECancelTan[b0, b2, c1, c2, t, deriv, n] solves the\n"
        "tangent-cancellation coupled system (Bronstein Chapter 8, the\n"
        "hypertangent case) for a hypertangent monomial t (Dt/(t^2+1) = eta in\n"
        "k): D[q1] + (b0 - n eta t) q1 - b2 q2 == c1 and D[q2] + b2 q1 +\n"
        "(b0 - n eta t) q2 == c2, for q1, q2 polynomials in t of degree at most n\n"
        "with b0, b2 in k and c1, c2 in k[t].  Recurses over k keeping t real,\n"
        "solving the base coupled system over k = C(x) via CoupledDESystem.\n"
        "Returns {q1, q2}, or stays unevaluated when there is no such solution.");
    rc_install("Risch`IntegrateHypertangentReduced", builtin_risch_int_hypertan_reduced,
        "Risch`IntegrateHypertangentReduced[p, t, deriv] integrates the reduced\n"
        "part of p over a hypertangent monomial t (Dt = eta (t^2+1)): peels the\n"
        "t^2+1 poles by solving a coupled system per multiplicity, returning\n"
        "{q, beta} with p - D[q] a polynomial in t when beta is True, or beta\n"
        "False when p has no elementary integral over the tangent field.");
}
