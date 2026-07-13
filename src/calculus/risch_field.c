/* risch_field.c — differential-field primitives for the Bronstein Risch algorithm.
 *
 * See risch_field.h for the contract.  Everything here reduces to Mathilda's
 * Expr builtins (which already dispatch to FLINT fast paths for rational
 * operands): the monomial derivation is a chain-rule sum of D[]'s, and the
 * field gcd / exact division in k[t] are PolynomialGCD / Cancel with the
 * pure-x "content" divided out so that k = C(other vars) behaves as a field.
 */

#include "risch_field.h"

#include "expr.h"
#include "eval.h"
#include "sym_intern.h"

#include <stdlib.h>

/* const-accepting deep copy (expr_copy takes a mutable pointer). */
static Expr* rf_cp(const Expr* e) { return expr_copy((Expr*)e); }

/* ------------------------------------------------------------------ */
/* Small evaluation helpers (adopt their argument Exprs, free the call). */
/* ------------------------------------------------------------------ */

static Expr* rf_eval_adopt(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}
static Expr* rf_fn(const char* head, Expr** args, size_t n) {
    return rf_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rf_call1(const char* head, Expr* a) {
    return rf_fn(head, (Expr*[]){ a }, 1);
}
static Expr* rf_call2(const char* head, Expr* a, Expr* b) {
    return rf_fn(head, (Expr*[]){ a, b }, 2);
}
static Expr* rf_call3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rf_fn(head, (Expr*[]){ a, b, c }, 3);
}

/* a * b  (unevaluated node, adopting a and b). */
static Expr* rf_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
/* a^n (unevaluated node, adopting a). */
static Expr* rf_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ a, expr_new_integer(n) }, 2);
}
/* -a (unevaluated node, adopting a). */
static Expr* rf_neg(Expr* a) {
    return expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), a }, 2);
}
/* a + b and a - b (evaluated). */
static Expr* rf_add(Expr* a, Expr* b) {
    return rf_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ a, b }, 2));
}
/* coefficient of t^k in p (an element of k). Owned. */
static Expr* rf_coeff(const Expr* p, const Expr* t, long k) {
    return rf_call3("Coefficient", rf_cp(p), rf_cp(t), expr_new_integer(k));
}

static bool rf_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("True");
}

/* ------------------------------------------------------------------ */
/* Derivation parsing.                                                 */
/* ------------------------------------------------------------------ */

bool risch_deriv_from_rules(const Expr* rules, RischDeriv* out) {
    out->nvars = 0; out->vars = NULL; out->dvars = NULL;
    if (!rules || rules->type != EXPR_FUNCTION) return false;
    const Expr* head = rules->data.function.head;
    if (!(head && head->type == EXPR_SYMBOL && head->data.symbol.name == intern_symbol("List")))
        return false;
    size_t n = rules->data.function.arg_count;
    if (n == 0) return false;
    Expr** vars = malloc(n * sizeof(Expr*));
    Expr** dvars = malloc(n * sizeof(Expr*));
    if (!vars || !dvars) { free(vars); free(dvars); return false; }
    for (size_t i = 0; i < n; i++) {
        const Expr* ri = rules->data.function.args[i];
        if (!ri || ri->type != EXPR_FUNCTION ||
            ri->data.function.head->type != EXPR_SYMBOL ||
            ri->data.function.head->data.symbol.name != intern_symbol("Rule") ||
            ri->data.function.arg_count != 2) {
            free(vars); free(dvars); return false;
        }
        vars[i] = ri->data.function.args[0];   /* borrowed */
        dvars[i] = ri->data.function.args[1];  /* borrowed */
    }
    out->nvars = n; out->vars = vars; out->dvars = dvars;
    return true;
}

void risch_deriv_free(RischDeriv* d) {
    if (!d) return;
    free(d->vars); free(d->dvars);
    d->vars = d->dvars = NULL; d->nvars = 0;
}

/* ------------------------------------------------------------------ */
/* Monomial derivation: D[p] = sum_i dvars[i] * d p / d vars[i].       */
/* ------------------------------------------------------------------ */

Expr* risch_field_deriv(const Expr* p, const RischDeriv* d) {
    size_t n = d->nvars;
    Expr** terms = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        Expr* dp = rf_call2("D", rf_cp(p), rf_cp(d->vars[i]));  /* d p / d var_i */
        terms[i] = rf_times(rf_cp(d->dvars[i]), dp);
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), terms, n);
    free(terms);
    return rf_eval_adopt(sum);
}

/* ------------------------------------------------------------------ */
/* Degree in t.                                                        */
/* ------------------------------------------------------------------ */

long risch_field_degree_t(const Expr* p, const Expr* t) {
    /* Length[CoefficientList[Expand[p], t]] - 1; the empty list (p == 0) maps
     * to -1. */
    Expr* ex = rf_call1("Expand", rf_cp(p));
    Expr* cl = rf_call2("CoefficientList", ex, rf_cp(t));
    long deg = -1;
    if (cl && cl->type == EXPR_FUNCTION &&
        cl->data.function.head->type == EXPR_SYMBOL &&
        cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        deg = (long)cl->data.function.arg_count - 1;
    }
    expr_free(cl);
    return deg;
}

/* Numerator[Together[a]] — clear pure-x denominators so k[t] elements become
 * genuine polynomials in t with polynomial (not rational) coefficients. */
static Expr* rf_clear_field_denoms(const Expr* a) {
    Expr* tg = rf_call1("Together", rf_cp(a));
    return rf_call1("Numerator", tg);  /* adopts tg */
}

/* ------------------------------------------------------------------ */
/* Field gcd in k[t]: monic-in-t gcd over k = C(other vars).          */
/* ------------------------------------------------------------------ */

Expr* risch_field_gcd_t(const Expr* a, const Expr* b, const Expr* t) {
    Expr* ca = rf_clear_field_denoms(a);
    Expr* cb = rf_clear_field_denoms(b);
    /* gcd(a, 0) = a and gcd(0, b) = b — PolynomialGCD does not reliably honour
     * a zero operand, so handle it explicitly (degree -1 marks the zero poly). */
    Expr* g;
    if (risch_field_degree_t(ca, t) < 0)      { g = cb; expr_free(ca); }
    else if (risch_field_degree_t(cb, t) < 0) { g = ca; expr_free(cb); }
    else g = rf_call2("PolynomialGCD", ca, cb);  /* over Q[all vars]; FLINT fast path */
    long dg = risch_field_degree_t(g, t);
    if (dg <= 0) {                  /* shares only a unit of the field k */
        expr_free(g);
        return expr_new_integer(1);
    }
    /* Divide by the leading coefficient in t (an element of k) to make monic. */
    Expr* lc = rf_call3("Coefficient", rf_cp(g), rf_cp(t), expr_new_integer(dg));
    Expr* monic = rf_call1("Cancel", rf_times(g, rf_pow(lc, -1)));  /* adopts g, lc */
    return monic;
}

/* ------------------------------------------------------------------ */
/* Field exact division in k(t): a/b, requiring a polynomial-in-t result. */
/* ------------------------------------------------------------------ */

Expr* risch_field_divexact_t(const Expr* a, const Expr* b, const Expr* t) {
    Expr* q = rf_call1("Cancel", rf_times(rf_cp(a), rf_pow(rf_cp(b), -1)));
    Expr* pq = rf_call2("PolynomialQ", rf_cp(q), rf_cp(t));
    bool ok = rf_is_true(pq);
    expr_free(pq);
    if (!ok) { expr_free(q); return NULL; }
    return q;
}

/* ------------------------------------------------------------------ */
/* Normal / special classification (Def. 3.4.2).                       */
/* ------------------------------------------------------------------ */

bool risch_field_is_normal(const Expr* p, const Expr* t, const RischDeriv* d) {
    Expr* Dp = risch_field_deriv(p, d);
    Expr* g = risch_field_gcd_t(p, Dp, t);
    long dg = risch_field_degree_t(g, t);
    expr_free(Dp); expr_free(g);
    return dg == 0;
}

bool risch_field_is_special(const Expr* p, const Expr* t, const RischDeriv* d) {
    long dp = risch_field_degree_t(p, t);
    if (dp <= 0) return false;   /* units / constants are not special polynomials */
    Expr* Dp = risch_field_deriv(p, d);
    Expr* g = risch_field_gcd_t(p, Dp, t);
    long dg = risch_field_degree_t(g, t);
    expr_free(Dp); expr_free(g);
    return dg == dp;
}

/* ------------------------------------------------------------------ */
/* Numerator / denominator over k(t), denominator made monic in t.     */
/* ------------------------------------------------------------------ */

void risch_field_num_den_t(const Expr* f, const Expr* t, Expr** a, Expr** d) {
    Expr* tg = rf_call1("Together", rf_cp(f));
    Expr* num = rf_call1("Numerator", rf_cp(tg));
    Expr* den = rf_call1("Denominator", tg);   /* adopts tg */
    long dd = risch_field_degree_t(den, t);
    if (dd <= 0) {
        /* denominator carries no t: absorb it into the numerator, d = 1. */
        *a = rf_call1("Cancel", rf_times(num, rf_pow(den, -1)));  /* num/den */
        *d = expr_new_integer(1);
        return;
    }
    /* make the denominator monic in t by dividing out its leading coefficient
     * (an element of k), pushing that unit into the numerator. */
    Expr* lc = rf_coeff(den, t, dd);
    *d = rf_call1("Cancel", rf_times(den, rf_pow(rf_cp(lc), -1)));  /* den/lc, adopts den */
    *a = rf_call1("Cancel", rf_times(num, rf_pow(lc, -1)));         /* num/lc, adopts num, lc */
}

/* ------------------------------------------------------------------ */
/* Polynomial division in k[t] (schoolbook long division; the field    */
/* coefficient arithmetic is handled by Cancel).                       */
/* ------------------------------------------------------------------ */

bool risch_field_divmod_t(const Expr* a, const Expr* b, const Expr* t,
                          Expr** qo, Expr** ro) {
    long db = risch_field_degree_t(b, t);
    if (db < 0) return false;                       /* division by zero */
    Expr* lcb = rf_coeff(b, t, db);
    Expr* q = expr_new_integer(0);
    Expr* r = rf_call1("Expand", rf_cp(a));
    for (;;) {
        long dr = risch_field_degree_t(r, t);
        if (dr < db) break;                         /* done (dr < 0 means r == 0) */
        Expr* lcr = rf_coeff(r, t, dr);
        Expr* ratio = rf_call1("Cancel", rf_times(lcr, rf_pow(rf_cp(lcb), -1)));
        Expr* term = (dr == db) ? ratio
                   : rf_eval_adopt(rf_times(ratio, rf_pow(rf_cp(t), dr - db)));
        q = rf_add(q, rf_cp(term));                 /* q += term */
        Expr* tb = rf_eval_adopt(rf_times(term, rf_cp(b)));  /* term*b, adopts term */
        r = rf_call1("Expand", rf_add(r, rf_neg(tb)));       /* r -= term*b */
    }
    expr_free(lcb);
    *qo = rf_call1("Cancel", q);
    *ro = rf_call1("Cancel", r);
    return true;
}

/* ------------------------------------------------------------------ */
/* Extended gcd in k[t] (iterative; u a + v b = g).                    */
/* ------------------------------------------------------------------ */

void risch_field_xgcd_t(const Expr* a, const Expr* b, const Expr* t,
                        Expr** go, Expr** uo, Expr** vo) {
    Expr* r0 = rf_call1("Expand", rf_cp(a));
    Expr* r1 = rf_call1("Expand", rf_cp(b));
    Expr* u0 = expr_new_integer(1), *u1 = expr_new_integer(0);
    Expr* v0 = expr_new_integer(0), *v1 = expr_new_integer(1);
    while (risch_field_degree_t(r1, t) >= 0) {      /* r1 != 0 */
        Expr* qq = NULL; Expr* rem = NULL;
        risch_field_divmod_t(r0, r1, t, &qq, &rem); /* r0 = qq r1 + rem */
        /* (u0,u1) <- (u1, u0 - qq u1); likewise v */
        Expr* u2 = rf_call1("Cancel",
            rf_add(rf_cp(u0), rf_neg(rf_eval_adopt(rf_times(rf_cp(qq), rf_cp(u1))))));
        Expr* v2 = rf_call1("Cancel",
            rf_add(rf_cp(v0), rf_neg(rf_eval_adopt(rf_times(rf_cp(qq), rf_cp(v1))))));
        expr_free(qq);
        expr_free(r0); r0 = r1; r1 = rem;
        expr_free(u0); u0 = u1; u1 = u2;
        expr_free(v0); v0 = v1; v1 = v2;
    }
    expr_free(r1); expr_free(u1); expr_free(v1);
    *go = r0; *uo = u0; *vo = v0;
}

/* ------------------------------------------------------------------ */
/* Diophantine solve b dn + c ds = r with deg_t(b) < deg_t(ds).        */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/* Derivation lookup and PolynomialReduce (Bronstein §5.4).            */
/* ------------------------------------------------------------------ */

const Expr* risch_deriv_lookup(const RischDeriv* d, const Expr* t) {
    if (!t || t->type != EXPR_SYMBOL) return NULL;
    for (size_t i = 0; i < d->nvars; i++)
        if (d->vars[i]->type == EXPR_SYMBOL && d->vars[i]->data.symbol.name == t->data.symbol.name)
            return d->dvars[i];
    return NULL;
}

bool risch_field_polynomial_reduce(const Expr* p, const Expr* t, const RischDeriv* d,
                                   Expr** qo, Expr** ro) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    long delta = risch_field_degree_t(Dt, t);      /* delta(t) = deg_t(Dt) */
    if (delta < 2) return false;                    /* requires a nonlinear monomial */
    Expr* lambda = rf_coeff(Dt, t, delta);          /* lc_t(Dt) */

    Expr* q = expr_new_integer(0);
    Expr* r = rf_call1("Expand", rf_cp(p));
    for (;;) {
        long dr = risch_field_degree_t(r, t);
        if (dr < delta) break;
        long m = dr - delta + 1;
        /* q0 = (lc_t(r) / (m lambda)) t^m ; D[q0] cancels the leading term of r. */
        Expr* lcr = rf_coeff(r, t, dr);
        Expr* mlam = rf_eval_adopt(rf_times(expr_new_integer(m), rf_cp(lambda)));
        Expr* c = rf_call1("Cancel", rf_times(lcr, rf_pow(mlam, -1)));
        Expr* q0 = rf_eval_adopt(rf_times(c, rf_pow(rf_cp(t), m)));  /* adopts c */
        q = rf_add(q, rf_cp(q0));                                   /* q += q0 */
        Expr* Dq0 = risch_field_deriv(q0, d);
        expr_free(q0);
        r = rf_call1("Expand", rf_add(r, rf_neg(Dq0)));            /* r -= D[q0] */
    }
    expr_free(lambda);
    *qo = rf_call1("Cancel", q);
    *ro = rf_call1("Cancel", r);
    return true;
}

void risch_field_diophantine_t(const Expr* dn, const Expr* ds, const Expr* r,
                               const Expr* t, Expr** bo, Expr** co) {
    Expr *g = NULL, *u = NULL, *v = NULL;
    risch_field_xgcd_t(dn, ds, t, &g, &u, &v);      /* u dn + v ds = g (g a unit) */
    /* normalize to u' dn + v' ds = 1 by dividing through by g */
    Expr* uu = rf_call1("Cancel", rf_times(u, rf_pow(rf_cp(g), -1)));  /* adopts u */
    Expr* vv = rf_call1("Cancel", rf_times(v, rf_pow(rf_cp(g), -1)));  /* adopts v */
    expr_free(g);
    /* particular solution b0 = r u', c0 = r v' */
    Expr* b0 = rf_call1("Expand", rf_eval_adopt(rf_times(rf_cp(r), uu)));
    Expr* c0 = rf_call1("Expand", rf_eval_adopt(rf_times(rf_cp(r), vv)));
    /* reduce b modulo ds: b0 = quo ds + b, then c = c0 + quo dn */
    Expr* quo = NULL; Expr* b = NULL;
    risch_field_divmod_t(b0, ds, t, &quo, &b);
    Expr* c = rf_call1("Cancel",
        rf_add(c0, rf_eval_adopt(rf_times(quo, rf_cp(dn)))));  /* adopts c0, quo */
    expr_free(b0);
    *bo = b; *co = c;
}
