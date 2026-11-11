/* residue.c -- Residue[expr, {z, z0}], the symbolic residue.
 *
 * The residue of f at an isolated singularity z = z0 is the coefficient of
 * (z - z0)^-1 in the Laurent expansion of f. We obtain it directly from the
 * series engine: expand f to order (z - z0)^0 (which always spans the -1 term,
 * however deep the pole), then read the coefficient at exponent -1 out of the
 * resulting SeriesData[z, z0, {coefs}, nmin, nmax, den].
 *
 * A residue is well defined only for an ordinary Laurent expansion (den == 1).
 * A fractional-power (Puiseux) expansion, den > 1, signals a branch point,
 * where the residue is undefined -- we leave the call unevaluated, matching
 * Mathematica (e.g. Residue[1/Sqrt[z], {z, 0}]).
 *
 * Algebraic pole locations.  The series engine decides whether z0 is a pole by
 * evaluating the denominator there and testing it against zero; but for a pole
 * whose location is a SUM of radicals (e.g. z0 = -2 + Sqrt[3], a root of
 * 1 + 4 z + z^2), Denominator(z0) is an expression like
 * 1 + 4 (-2 + Sqrt[3]) + (-2 + Sqrt[3])^2 that does not auto-simplify to 0, so
 * the pole is missed and the residue wrongly comes out 0.  We defeat this by
 * expanding about z0 EXPLICITLY: substitute z -> z0 + w, then Expand the
 * denominator of the result -- polynomial expansion collapses the radical
 * arithmetic (Sqrt[3]^2 -> 3, ...) so the vanishing constant term becomes a
 * literal 0 and the w-factor of the pole is exposed.  Reading the (z-z0)^-1
 * coefficient is then a plain Series-at-0 of the expanded form.
 */

#include "residue.h"
#include "series.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "internal.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "numeric.h"   /* expr_is_numeric_quantity */
#include "rationalize.h" /* internal_contains_inexact */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* Emit the too-few-arguments diagnostic and leave the call unevaluated. */
static Expr* residue_emit_argcount(size_t argc) {
    mth_message("Residue", "argm", "Residue called with %zu argument%s; 2 or more arguments are expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Build and evaluate Series[expr, {z, z0, order}]. expr_new_function copies the
 * args array (memcpy) and adopts the element pointers, so the temporary arrays
 * are freed while their contents live on in the new nodes. Returns the
 * (owned) result of evaluation, which the caller must free. */
static Expr* residue_series(Expr* expr, Expr* z, Expr* z0, int64_t order) {
    Expr** spec_args = calloc(3, sizeof(Expr*));
    if (!spec_args) return NULL;
    spec_args[0] = expr_copy(z);
    spec_args[1] = expr_copy(z0);
    spec_args[2] = expr_new_integer(order);
    Expr* series_spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 3);
    free(spec_args);

    Expr** call_args = calloc(2, sizeof(Expr*));
    if (!call_args) { expr_free(series_spec); return NULL; }
    call_args[0] = expr_copy(expr);
    call_args[1] = series_spec;
    Expr* series_call = expr_new_function(expr_new_symbol("Series"), call_args, 2);
    free(call_args);

    return eval_and_free(series_call);
}

/* Build head[a] and evaluate it, freeing the call.  `a` is consumed. */
static Expr* residue_eval1(const char* head, Expr* a) {
    Expr** args = calloc(1, sizeof(Expr*));
    if (!args) { expr_free(a); return NULL; }
    args[0] = a;
    Expr* call = expr_new_function(expr_new_symbol(head), args, 1);
    free(args);
    return eval_and_free(call);
}

/* Expand[p /. z -> z0 + w]: shift a polynomial `p` in z to the local variable w
 * about z0 and expand it out.  `p` is consumed; `z`, `z0`, `w` borrowed.
 * Returns an owned polynomial in w, or NULL. */
static Expr* residue_shift_poly(Expr* p, Expr* z, Expr* z0, Expr* w) {
    Expr* shift = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(z0), expr_copy(w) }, 2);
    Expr* rule  = expr_new_function(expr_new_symbol(SYM_Rule),
                      (Expr*[]){ expr_copy(z), shift }, 2);
    Expr* repl  = expr_new_function(expr_new_symbol("ReplaceAll"),
                      (Expr*[]){ p, rule }, 2);       /* adopts p */
    Expr* g = eval_and_free(repl);
    if (!g) return NULL;
    return residue_eval1("Expand", g);                /* consumes g */
}

/* Rewrite the rational integrand `f` so that its singularity at z = z0 is
 * expanded about w = 0: returns Expand[P(z0+w)] / Expand[Q(z0+w)] in the fresh
 * variable `w`, where P/Q = Together[f].  The numerator and denominator
 * polynomials are shifted SEPARATELY and each expanded, rather than Together-ing
 * the shifted integrand.  This is what keeps an algebraic pole location (z0 a
 * radical) tractable: the binomial expansion of Q(z0+w) only ever produces
 * powers of the single radical present in z0, so its coefficients stay in one
 * minimal spelling that the polynomial Expand engine can reduce (e.g.
 * ((-1)^(1/4))^4 -> -1).  Together-ing the shifted rational instead re-normalises
 * the algebraic field and can emit a SECOND, independent spelling of the same
 * number (e.g. (-1)^(3/4) beside (-1)^(1/4)); Expand then treats the two
 * spellings as independent generators, the base-power reductions no longer fire,
 * and the downstream series inversion blows up combinatorially.  P and Q are
 * coprime, so shifting them without re-cancelling loses no common w-factor and
 * the w^-1 coefficient -- the residue -- is unchanged.  The expansion also
 * collapses Q's constant term to a literal 0, exposing the pole as a w-factor.
 * Returns an owned Expr* (in `w`) or NULL.  `f`, `z`, `z0`, `w` are borrowed. */
static Expr* residue_shift_form(Expr* f, Expr* z, Expr* z0, Expr* w) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return NULL;
    Expr* P = residue_eval1("Numerator", expr_copy(tog));
    Expr* Q = residue_eval1("Denominator", tog);       /* consumes tog */
    if (!P || !Q) { if (P) expr_free(P); if (Q) expr_free(Q); return NULL; }

    Expr* A = residue_shift_poly(P, z, z0, w);         /* consumes P */
    if (!A) { expr_free(Q); return NULL; }
    Expr* B = residue_shift_poly(Q, z, z0, w);         /* consumes Q */
    if (!B) { expr_free(A); return NULL; }

    /* form = A * B^-1 */
    Expr* bpow = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ B, expr_new_integer(-1) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ A, bpow }, 2);
}

/* True iff PolynomialQ[p, z]. */
static bool residue_polyq(Expr* p, Expr* z) {
    Expr* call = internal_polynomialq(
        (Expr*[]){ expr_copy(p), expr_copy(z) }, 2);
    Expr* v = eval_and_free(call);
    bool ok = v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
    if (v) expr_free(v);
    return ok;
}

/* True iff f is a rational function of z: Together[f] has polynomial numerator
 * AND polynomial denominator in z.  The shift+Expand pole-detection preprocessing
 * is applied ONLY in this case; for transcendental / special-function integrands
 * (Cot, Zeta near its pole, unknown f[z], ...) the series engine's built-in
 * knowledge of the expansion at z0 is preferable, and shifting the argument would
 * defeat it. */
static bool residue_is_rational_in(Expr* f, Expr* z) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return false;
    Expr* num = residue_eval1("Numerator", expr_copy(tog));
    Expr* den = residue_eval1("Denominator", tog);   /* consumes tog */
    if (!num || !den) { if (num) expr_free(num); if (den) expr_free(den); return false; }
    bool ok = residue_polyq(num, z) && residue_polyq(den, z);
    expr_free(num);
    expr_free(den);
    return ok;
}

/* Build head[a, b] and evaluate it, freeing the call.  `a`, `b` consumed. */
static Expr* residue_eval2(const char* head, Expr* a, Expr* b) {
    Expr* call = expr_new_function(expr_new_symbol(head),
                                   (Expr*[]){ a, b }, 2);
    return eval_and_free(call);
}

/* expr /. z -> z0, evaluated.  `expr`, `z`, `z0` borrowed; returns owned. */
static Expr* residue_subst(Expr* expr, Expr* z, Expr* z0) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(z), expr_copy(z0) }, 2);
    return residue_eval2("ReplaceAll", expr_copy(expr), rule);
}

/* True iff PossibleZeroQ[e] is True.  `e` borrowed. */
static bool residue_is_zero(Expr* e) {
    Expr* pz = residue_eval1("PossibleZeroQ", expr_copy(e));
    bool z = pz && pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True;
    if (pz) expr_free(pz);
    return z;
}

/* Simple-pole fast path for a rational integrand.  For f = P/Q with a SIMPLE
 * zero of Q at z0, Res_{z0} f = P(z0)/Q'(z0).  This bypasses the Laurent-series
 * expansion, which for an algebraic pole location (z0 a nested radical such as
 * (-1)^(1/4)) can blow up catastrophically in the generic series inverter
 * (unbounded growth of un-reduced radical coefficients).  Returns an owned
 * residue Expr* on success, or NULL to signal "not a decidable simple pole" so
 * the caller falls back to the series engine.  `f`, `z`, `z0` borrowed. */
static Expr* residue_simple_pole(Expr* f, Expr* z, Expr* z0) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return NULL;
    Expr* P = residue_eval1("Numerator", expr_copy(tog));
    Expr* Q = residue_eval1("Denominator", tog);   /* consumes tog */
    if (!P || !Q) { if (P) expr_free(P); if (Q) expr_free(Q); return NULL; }

    /* z0 must actually be a pole: Q(z0) == 0.  If not provably zero, defer to the
     * series engine (which correctly returns 0 at an analytic point). */
    Expr* Qat = residue_subst(Q, z, z0);
    bool q_zero = Qat && residue_is_zero(Qat);
    if (Qat) expr_free(Qat);
    if (!q_zero) { expr_free(P); expr_free(Q); return NULL; }

    /* Simple pole requires Q'(z0) != 0. */
    Expr* Qp = residue_eval2("D", Q, expr_copy(z));   /* consumes Q */
    if (!Qp) { expr_free(P); return NULL; }
    Expr* Qpat = residue_subst(Qp, z, z0);
    expr_free(Qp);
    if (!Qpat || residue_is_zero(Qpat)) {   /* undecidable or higher-order pole */
        expr_free(P); if (Qpat) expr_free(Qpat); return NULL;
    }

    /* Res = P(z0) / Q'(z0). */
    Expr* Pat = residue_subst(P, z, z0);
    expr_free(P);
    if (!Pat) { expr_free(Qpat); return NULL; }
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ Qpat, expr_new_integer(-1) }, 2);
    return residue_eval2("Times", Pat, inv);
}

/* Read the (svar-spt)^-1 coefficient from a single Series expansion to the given
 * target `order`. Sets *status: 0 = coefficient returned (owned); 1 = analytic
 * point (no principal part), returns NULL; 2 = failure / branch point, NULL;
 * 3 = the O-term lies at or before exponent -1 so the coefficient is not yet
 * resolved at this order, NULL. *pole_m (when non-NULL) receives the pole order
 * -nmin (>=1) or 0 at an analytic point, for the caller's convergence sizing. */
static Expr* residue_coeff_once(Expr* expr, Expr* svar, Expr* spt,
                                int64_t order, int* status, int64_t* pole_m) {
    *status = 2;
    if (pole_m) *pole_m = 0;
    Expr* sd = residue_series(expr, svar, spt, order);
    if (!sd) return NULL;
    if (!is_series_data(sd)) { expr_free(sd); return NULL; }

    Expr** a     = sd->data.function.args;
    Expr* coefs  = a[2];
    Expr* nmin_e = a[3];
    Expr* nmax_e = a[4];
    Expr* den_e  = a[5];
    if (nmin_e->type != EXPR_INTEGER || nmax_e->type != EXPR_INTEGER ||
        den_e->type != EXPR_INTEGER || coefs->type != EXPR_FUNCTION) {
        expr_free(sd);
        return NULL;
    }
    int64_t nmin = nmin_e->data.integer;
    int64_t den  = den_e->data.integer;

    /* Fractional exponents -> branch point -> residue undefined. */
    if (den != 1) { expr_free(sd); return NULL; }

    if (pole_m) *pole_m = (nmin <= -1) ? -nmin : 0;

    int64_t len   = (int64_t)coefs->data.function.arg_count;
    int64_t index = -1 - nmin;
    if (index < 0) { expr_free(sd); *status = 1; return NULL; }  /* analytic */
    if (index < len) {
        Expr* r = expr_copy(coefs->data.function.args[index]);
        expr_free(sd);
        *status = 0;
        return r;
    }
    expr_free(sd);
    *status = 3;   /* O-term <= -1: need a higher order */
    return NULL;
}

/* Resolve the (svar-spt)^-1 coefficient starting from `start_order`, raising the
 * order while the O-term still hides exponent -1. Returns the owned coefficient
 * (status 0), an owned integer 0 (analytic, status 1), or NULL (branch/failure).
 * On success *used is the order actually used and *pole_m the pole order. */
static Expr* residue_coeff_resolved(Expr* expr, Expr* svar, Expr* spt,
                                    int64_t start_order, int64_t* used,
                                    int64_t* pole_m) {
    int64_t order = start_order < 0 ? 0 : start_order;
    for (int attempt = 0; attempt < 10; attempt++) {
        int status; int64_t m = 0;
        Expr* c = residue_coeff_once(expr, svar, spt, order, &status, &m);
        if (pole_m) *pole_m = m;
        if (status == 0) { if (used) *used = order; return c; }
        if (status == 1) { if (used) *used = order; return expr_new_integer(0); }
        if (status == 2) return NULL;
        /* status == 3: O-term still at/below -1; raise and retry. */
        order += 3;
        if (order > 256) return NULL;
    }
    return NULL;
}

/* True iff PossibleZeroQ[a - b] is True.  `a`, `b` borrowed. */
static bool residue_coeffs_agree(Expr* a, Expr* b) {
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Plus),
                     (Expr*[]){ expr_copy(a),
                                expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2) }, 2);
    bool z = residue_is_zero(diff);
    expr_free(diff);
    return z;
}

/* Self-validating Laurent-coefficient extraction. The series engine can return a
 * (svar-spt)^-1 coefficient that is PRESENT but inaccurate: a higher-order pole's
 * regular cofactor is truncated too early, so the principal part drops product-
 * rule cross terms (the symbolic double/triple-pole bug). A padded series still
 * claims validity at that order, so a single read cannot detect the loss. Instead
 * we exploit that the coefficient CONVERGES once the order reaches pole_order-1
 * and is exactly stable thereafter: compute it at successive orders and accept
 * only when two consecutive orders agree. Returns an owned copy of the residue,
 * 0 at an analytic point, or NULL at a branch point / when no series converges. */
static Expr* residue_extract(Expr* expr, Expr* svar, Expr* spt) {
    /* Probe at order 0 to learn the pole order (and short-circuit the analytic
     * and branch-point cases). */
    int st0; int64_t m0 = 0;
    Expr* c0 = residue_coeff_once(expr, svar, spt, 0, &st0, &m0);
    if (st0 == 2) return NULL;                 /* branch / failure */
    if (st0 == 1) return expr_new_integer(0);  /* analytic -> residue 0 */
    if (c0) { expr_free(c0); }                 /* discard: we re-derive with margin */

    int64_t m = (m0 >= 1) ? m0 : 1;
    int64_t used = 0, pm = 0;
    /* Start a safe margin above the convergence threshold (pole_order - 1). */
    Expr* prev = residue_coeff_resolved(expr, svar, spt, m - 1, &used, &pm);
    if (!prev) return NULL;

    int64_t ord = used;
    for (int it = 0; it < 12; it++) {
        int64_t used2 = 0, pm2 = 0;
        Expr* cur = residue_coeff_resolved(expr, svar, spt, ord + 1, &used2, &pm2);
        if (!cur) return prev;                 /* can't go higher: trust prev */
        if (residue_coeffs_agree(prev, cur)) { expr_free(prev); return cur; }
        expr_free(prev);
        prev = cur;
        ord = used2;
        if (ord > 64) return prev;             /* bounded effort */
    }
    return prev;
}

/* Does f have a pole (not an analytic point) at z = z0?  True when the
 * denominator of Together[f] vanishes there. Borrowed args. */
static bool residue_has_pole_at(Expr* f, Expr* z, Expr* z0) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return false;
    Expr* Q = residue_eval1("Denominator", tog);   /* consumes tog */
    if (!Q) return false;
    Expr* Qat = residue_subst(Q, z, z0);
    expr_free(Q);
    bool zero = Qat && residue_is_zero(Qat);
    if (Qat) expr_free(Qat);
    return zero;
}

/* True iff `e` is an Integer, BigInt or Rational. */
static bool residue_is_rational_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Rational;
}

/* True iff `e` is built from z, rational numbers, Plus, Times and integer
 * Powers -- a rational function of z over Q, in any (unexpanded) spelling. */
static bool residue_rational_over_q(const Expr* e, const Expr* z) {
    if (residue_is_rational_number(e)) return true;
    if (e->type == EXPR_SYMBOL) return expr_eq(e, z);
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Power)
        return n == 2 && e->data.function.args[1]->type == EXPR_INTEGER &&
               residue_rational_over_q(e->data.function.args[0], z);
    if (h != SYM_Plus && h != SYM_Times) return false;
    for (size_t i = 0; i < n; i++)
        if (!residue_rational_over_q(e->data.function.args[i], z)) return false;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Essential singularity: Laurent polynomial times Exp[c + a w + b/w]   *
 * ------------------------------------------------------------------ *
 *
 * No Laurent series with a finite principal part exists at an essential
 * singularity, so the series engine cannot read the residue there.  For
 * f = P(w) e^(c + a w + b/w), w = z - z0, P a Laurent polynomial in w, the
 * residue is a convergent double series with a closed form.  Writing
 * e^(a w) = sum a^m w^m/m! and e^(b/w) = sum b^n w^-n/n!, the w^-1 coefficient
 * of w^j e^(a w + b/w) collects the pairs with j + m - n = -1:
 *     R_nu = sum_m a^m b^(m+nu)/(m! (m+nu)!),   nu = j + 1,
 * the series of the modified Bessel function,
 *     nu >= 0:  R_nu = b^nu (a b)^(-nu/2) I_nu(2 Sqrt[a b]),
 *     nu <  0:  R_nu = a^-nu (a b)^(nu/2) I_-nu(2 Sqrt[a b])
 * ((a b)^(-n/2) I_n(2 Sqrt[a b]) = sum (a b)^m/(m! (m+n)!) is single-valued,
 * so the branch of the square root is immaterial).  At a = 0 only m = 0
 * survives: R_nu = b^nu/nu! for nu >= 0, else 0.  The residue is
 * e^c sum_j p_j R_(j+1). */

static bool residue_free_of(const Expr* e, const Expr* w) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) return !expr_eq((Expr*)e, (Expr*)w);
    if (e->type != EXPR_FUNCTION) return true;
    if (!residue_free_of(e->data.function.head, w)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!residue_free_of(e->data.function.args[i], w)) return false;
    return true;
}

static bool residue_head_is(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* Split a monomial term t = coef * w^j (coef free of w, j an integer).
 * Returns the owned coefficient and sets *j, or NULL if t is not of that form. */
static Expr* residue_monomial(Expr* t, Expr* w, int64_t* j) {
    Expr* one[1] = { t };
    Expr** fs = one;
    size_t n = 1;
    if (residue_head_is(t, SYM_Times)) {
        fs = t->data.function.args;
        n = t->data.function.arg_count;
    }
    *j = 0;
    Expr** keep = malloc((n ? n : 1) * sizeof(Expr*));
    if (!keep) return NULL;
    size_t nk = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = fs[i];
        if (expr_eq(f, w)) { *j += 1; continue; }
        if (residue_head_is(f, SYM_Power) && f->data.function.arg_count == 2 &&
            expr_eq(f->data.function.args[0], w) &&
            f->data.function.args[1]->type == EXPR_INTEGER) {
            *j += f->data.function.args[1]->data.integer;
            continue;
        }
        if (!residue_free_of(f, w)) {
            for (size_t k = 0; k < nk; k++) expr_free(keep[k]);
            free(keep);
            return NULL;
        }
        keep[nk++] = expr_copy(f);
    }
    Expr* c = nk == 0 ? expr_new_integer(1)
            : nk == 1 ? keep[0]
            : expr_new_function(expr_new_symbol(SYM_Times), keep, nk);
    free(keep);
    return c;
}

/* R_nu of the comment above, unevaluated (a, b borrowed). */
static Expr* residue_bessel_term(Expr* a, Expr* b, int64_t nu) {
    if (a->type == EXPR_INTEGER && a->data.integer == 0) {
        if (nu < 0) return expr_new_integer(0);
        Expr* fact = expr_new_function(expr_new_symbol(SYM_Factorial),
                         (Expr*[]){ expr_new_integer(nu) }, 1);
        return expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ expr_copy(b), expr_new_integer(nu) }, 2),
                              expr_new_function(expr_new_symbol(SYM_Power),
                                  (Expr*[]){ fact, expr_new_integer(-1) }, 2) }, 2);
    }
    int64_t m = nu >= 0 ? nu : -nu;
    Expr* ab = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
    Expr* lead = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(nu >= 0 ? b : a), expr_new_integer(m) }, 2);
    Expr* half = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_new_integer(-m),
                                expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_new_integer(2), expr_new_integer(-1) }, 2) }, 2);
    Expr* scale = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_copy(ab), half }, 2);
    Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(2),
                               expr_new_function(expr_new_symbol(SYM_Sqrt),
                                   (Expr*[]){ ab }, 1) }, 2);
    Expr* bes = expr_new_function(expr_new_symbol(SYM_BesselI),
                    (Expr*[]){ expr_new_integer(m), arg }, 2);
    return expr_new_function(expr_new_symbol(SYM_Times),
               (Expr*[]){ lead, scale, bes }, 3);
}

/* Split an expanded expression into (coefficient, w-exponent) monomials.
 * Returns the count (0 on a non-monomial term); coefs owned by the caller. */
#define RESIDUE_MAX_TERMS 4096
static size_t residue_split_terms(Expr* e, Expr* w, Expr** coefs, int64_t* exps) {
    Expr* one[1] = { e };
    Expr** ts = one;
    size_t nt = 1;
    if (residue_head_is(e, SYM_Plus)) { ts = e->data.function.args; nt = e->data.function.arg_count; }
    if (nt > RESIDUE_MAX_TERMS) return 0;
    for (size_t i = 0; i < nt; i++) {
        coefs[i] = residue_monomial(ts[i], w, &exps[i]);
        if (!coefs[i]) {
            for (size_t k = 0; k < i; k++) expr_free(coefs[k]);
            return 0;
        }
    }
    return nt;
}

static Expr* residue_exp_essential(Expr* f, Expr* z, Expr* z0) {
    Expr* w = expr_new_symbol("Residue`$w");
    Expr* shift = expr_new_function(expr_new_symbol(SYM_Plus),
                      (Expr*[]){ expr_copy(z0), expr_copy(w) }, 2);
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(z), shift }, 2);
    Expr* fs = residue_eval2("ReplaceAll", expr_copy(f), rule);
    Expr* out = NULL;
    Expr *ex = NULL, *Pe = NULL, *a = NULL, *b = NULL, *c = NULL;
    Expr** coefs = NULL;
    int64_t* exps = NULL;
    if (!fs) goto done;
    {
        /* Separate the exponential factors (E^g, g depending on w) from the rest. */
        Expr* one[1] = { fs };
        Expr** fl = one;
        size_t n = 1;
        if (residue_head_is(fs, SYM_Times)) { fl = fs->data.function.args; n = fs->data.function.arg_count; }
        Expr** rest = malloc(n * sizeof(Expr*));
        Expr** ex_t = malloc(n * sizeof(Expr*));
        size_t nr = 0, ne = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* g = fl[i];
            if (residue_head_is(g, SYM_Power) && g->data.function.arg_count == 2 &&
                g->data.function.args[0]->type == EXPR_SYMBOL &&
                g->data.function.args[0]->data.symbol.name == SYM_E &&
                !residue_free_of(g->data.function.args[1], w))
                ex_t[ne++] = expr_copy(g->data.function.args[1]);
            else
                rest[nr++] = expr_copy(g);
        }
        if (ne == 0) {
            for (size_t i = 0; i < nr; i++) expr_free(rest[i]);
            free(rest); free(ex_t);
            goto done;
        }
        Expr* expo = ne == 1 ? ex_t[0] : expr_new_function(expr_new_symbol(SYM_Plus), ex_t, ne);
        Expr* P = nr == 0 ? expr_new_integer(1)
                : nr == 1 ? rest[0] : expr_new_function(expr_new_symbol(SYM_Times), rest, nr);
        free(rest); free(ex_t);
        ex = residue_eval1("Expand", expo);
        Pe = residue_eval1("Expand", P);
    }
    if (!ex || !Pe) goto done;
    coefs = malloc(RESIDUE_MAX_TERMS * sizeof(Expr*));
    exps = malloc(RESIDUE_MAX_TERMS * sizeof(int64_t));
    if (!coefs || !exps) goto done;
    {
        /* exponent = c + a w + b/w, exactly */
        size_t nt = residue_split_terms(ex, w, coefs, exps);
        if (nt == 0) goto done;
        Expr** part[3];
        size_t np[3] = { 0, 0, 0 };
        bool ok = true;
        for (int k = 0; k < 3; k++) part[k] = malloc(nt * sizeof(Expr*));
        for (size_t i = 0; i < nt; i++) {
            if (exps[i] < -1 || exps[i] > 1) { ok = false; expr_free(coefs[i]); continue; }
            part[exps[i] + 1][np[exps[i] + 1]++] = coefs[i];
        }
        Expr* v[3];
        for (int k = 0; k < 3; k++) {
            if (!ok) { for (size_t i = 0; i < np[k]; i++) expr_free(part[k][i]); v[k] = NULL; }
            else if (np[k] == 0) v[k] = expr_new_integer(0);
            else v[k] = eval_and_free(np[k] == 1 ? part[k][0]
                          : expr_new_function(expr_new_symbol(SYM_Plus), part[k], np[k]));
            free(part[k]);
        }
        b = v[0]; c = v[1]; a = v[2];
        if (!ok || !a || !b || !c) goto done;
        if (b->type == EXPR_INTEGER && b->data.integer == 0) goto done;   /* not essential */
    }
    {
        /* P = sum_j p_j w^j, a Laurent polynomial in w */
        size_t nt = residue_split_terms(Pe, w, coefs, exps);
        if (nt == 0) goto done;
        Expr** acc = malloc(nt * sizeof(Expr*));
        for (size_t i = 0; i < nt; i++)
            acc[i] = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ coefs[i], residue_bessel_term(a, b, exps[i] + 1) }, 2);
        Expr* sum = nt == 1 ? acc[0] : expr_new_function(expr_new_symbol(SYM_Plus), acc, nt);
        free(acc);
        Expr* ec = expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_symbol(SYM_E), expr_copy(c) }, 2);
        out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                              (Expr*[]){ ec, sum }, 2));
    }
done:
    free(coefs);
    free(exps);
    if (fs) expr_free(fs);
    if (ex) expr_free(ex);
    if (Pe) expr_free(Pe);
    if (a) expr_free(a);
    if (b) expr_free(b);
    if (c) expr_free(c);
    expr_free(w);
    return out;
}

/* True iff `e` is an ENTIRE function of z built syntactically from entire
 * pieces: z-free leaves, z, Plus/Times, Power with a non-negative integer
 * exponent or with a z-free base (E^g, c^g), and Exp/Sin/Cos/Sinh/Cosh of an
 * entire argument.  Conservative: any other head is "not known entire". */
static bool residue_entire_in(const Expr* e, const Expr* z) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return true;           /* atom: z or z-free */
    if (residue_free_of(e, z)) return true;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Plus || h == SYM_Times) {
        for (size_t i = 0; i < n; i++)
            if (!residue_entire_in(e->data.function.args[i], z)) return false;
        return true;
    }
    if (h == SYM_Power && n == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* p = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer >= 0)
            return residue_entire_in(b, z);
        return residue_free_of(b, z) && residue_entire_in(p, z);
    }
    if ((h == SYM_Exp || h == SYM_Sin || h == SYM_Cos || h == SYM_Sinh ||
         h == SYM_Cosh) && n == 1)
        return residue_entire_in(e->data.function.args[0], z);
    return false;
}

/* True iff Together[f] = N/Q with Q a polynomial in z and N entire in z (e.g.
 * E^(I a z)/(c^2 + (b^2 + z^2)^2), the Jordan-lemma integrand): every
 * singularity of f is a zero of Q, so at a simple zero z0 of Q the residue is
 * exactly N(z0)/Q'(z0). */
static bool residue_entire_over_poly(Expr* f, Expr* z) {
    Expr* tog = residue_eval1("Together", expr_copy(f));
    if (!tog) return false;
    Expr* num = residue_eval1("Numerator", expr_copy(tog));
    Expr* den = residue_eval1("Denominator", tog);   /* consumes tog */
    bool ok = num && den && residue_polyq(den, z) && residue_entire_in(num, z);
    if (num) expr_free(num);
    if (den) expr_free(den);
    return ok;
}

Expr* residue_compute(Expr* f, Expr* z, Expr* z0) {
    if (!f || !z || !z0 || z->type != EXPR_SYMBOL) return NULL;

    /* A rational function over Q at a rational point: the residue is ONE
     * Laurent coefficient, read off a single expansion about z0.  The series
     * engine expands such an f in exact rational arithmetic (every
     * coefficient rational, every truncation order tracked exactly), so the
     * coefficient it reports is the coefficient -- none of the cross-checking
     * below, which guards symbolic coefficients, and no shift-and-Expand,
     * which exists to expose an ALGEBRAIC pole location.  For p(z)/z^k at 0
     * this is the z^(k-1) coefficient of p, at the cost of one truncated
     * product of p's factors. */
    if (residue_is_rational_number(z0) && residue_rational_over_q(f, z)) {
        int64_t used = 0, pm = 0;
        Expr* r = residue_coeff_resolved(f, z, z0, 0, &used, &pm);
        if (r) return r;
    }

    if (residue_is_rational_in(f, z)) {
        /* Simple-pole fast path: Res = P(z0)/Q'(z0) when Q has a simple zero at
         * z0.  Skips the Laurent-series inversion, which blows up on algebraic
         * pole locations (nested radicals).  NULL falls through to the series
         * engine, which still owns higher-order poles and undecidable cases. */
        Expr* fast = residue_simple_pole(f, z, z0);
        if (fast) return fast;

        /* Rational integrand: expand about z0 with an EXPANDED denominator so an
         * algebraic pole location (z0 a sum of radicals) is exposed as a w-factor. */
        Expr* w = expr_new_symbol("Residue`$w");
        Expr* form = residue_shift_form(f, z, z0, w);
        if (!form) { expr_free(w); return NULL; }
        Expr* zero = expr_new_integer(0);
        Expr* result = residue_extract(form, w, zero);
        expr_free(form);
        expr_free(w);
        expr_free(zero);
        return result;
    }

    /* Entire numerator over a polynomial denominator (the Jordan-lemma
     * integrand E^(I a z) R(z)): at a simple zero of Q the residue is
     * N(z0)/Q'(z0), exactly as in the rational fast path.  The series engine
     * below decides "is z0 a pole" by evaluating Q(z0), which does not
     * auto-simplify to 0 for a nested-radical root (z0 =
     * Sqrt[(-2 b^2 + Sqrt[-4 c^2])/2] of c^2 + (b^2 + z^2)^2) -- it then reads
     * an analytic expansion and reports a WRONG residue 0. */
    bool entire_poly = residue_entire_over_poly(f, z);
    if (entire_poly) {
        Expr* fast = residue_simple_pole(f, z, z0);
        if (fast) return fast;
    }

    /* An essential singularity of e^(b/(z - z0)) type: no Laurent series with
     * a finite principal part exists, so the residue is read in closed form. */
    {
        Expr* ess = residue_exp_essential(f, z, z0);
        if (ess) return ess;
    }

    /* Transcendental / special-function integrand: expand directly about z0 so
     * the series engine can use its knowledge of the function's Laurent series
     * there (e.g. Zeta at 1, Cot / 1/Sin^n at 0, unknown f[z]/z^n). */
    Expr* r = residue_extract(f, z, z0);

    /* Dropped-pole guard. A concrete odd-denominator rational power composed
     * DIRECTLY with a pole at a nonzero z0 (e.g. z^(1/3)/(1+z^2)^2 at z = I)
     * makes the series engine lose the principal part, so residue_extract reads
     * an analytic series and returns 0 -- a silent WRONG answer. When z0 is in
     * fact a pole, re-expand the SHIFTED integrand f /. z -> z0 + w about w = 0:
     * there z^p becomes the analytic binomial (z0+w)^p and the pole sits in w
     * directly, which the series engine handles correctly. */
    if (r && r->type == EXPR_INTEGER && r->data.integer == 0 &&
        residue_has_pole_at(f, z, z0)) {
        Expr* w = expr_new_symbol("Residue`$w");
        Expr* shiftpt = expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_copy(z0), expr_copy(w) }, 2);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                         (Expr*[]){ expr_copy(z), shiftpt }, 2);
        Expr* fs = residue_eval2("ReplaceAll", expr_copy(f), rule);
        Expr* zero = expr_new_integer(0);
        if (fs) {
            Expr* r2 = residue_extract(fs, w, zero);
            expr_free(fs);
            if (r2 && !(r2->type == EXPR_INTEGER && r2->data.integer == 0)) {
                expr_free(r);
                expr_free(w); expr_free(zero);
                return r2;
            }
            if (r2) expr_free(r2);
        }
        /* Entire/polynomial integrand at a higher-order algebraic pole: the
         * shifted integrand still hides the vanishing constant term of
         * Q(z0 + w) behind unexpanded radicals.  Expand numerator and
         * denominator about z0 separately (as for a rational f) so the pole
         * is exposed as an explicit power of w. */
        if (entire_poly) {
            Expr* form = residue_shift_form(f, z, z0, w);
            Expr* r3 = form ? residue_extract(form, w, zero) : NULL;
            if (form) expr_free(form);
            if (r3 && !(r3->type == EXPR_INTEGER && r3->data.integer == 0)) {
                expr_free(r);
                expr_free(w); expr_free(zero);
                return r3;
            }
            if (r3) expr_free(r3);
        }
        expr_free(w); expr_free(zero);
    }
    return r;
}

static size_t residue_leaves(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 1;
    size_t c = residue_leaves(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += residue_leaves(e->data.function.args[i]);
    return c;
}

static bool residue_has_root(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "Root") == 0) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (residue_has_root(e->data.function.args[i])) return true;
    return false;
}

/* True iff `e` has a non-atomic denominator: Power[b, negative] with b a
 * compound expression (a sum or a radical), as in 180/(-18 - 18 I Sqrt[3])^2
 * or 1/(2 Sqrt[3]).  A literal Rational is not a denominator here. */
static bool residue_has_compound_denominator(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* p = e->data.function.args[1];
        bool neg = (p->type == EXPR_INTEGER && p->data.integer < 0) ||
                   (p->type == EXPR_FUNCTION &&
                    p->data.function.head->type == EXPR_SYMBOL &&
                    p->data.function.head->data.symbol.name == SYM_Rational &&
                    p->data.function.args[0]->type == EXPR_INTEGER &&
                    p->data.function.args[0]->data.integer < 0);
        if (neg && b->type == EXPR_FUNCTION) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (residue_has_compound_denominator(e->data.function.args[i])) return true;
    return false;
}

/* Total bit length of the integer leaves of `e`: the tie-break between two
 * equally large spellings, so (54 - 3 I Sqrt[3])/432 yields to
 * (18 - I Sqrt[3])/144.  A bigint counts a flat 64 bits: only the ordering
 * of two small spellings matters. */
static size_t residue_int_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        uint64_t v = e->data.integer < 0 ? (uint64_t)(-(e->data.integer + 1)) + 1
                                         : (uint64_t)e->data.integer;
        size_t b = 0;
        while (v) { b++; v >>= 1; }
        return b;
    }
    if (e->type == EXPR_BIGINT) return 64;
    if (e->type != EXPR_FUNCTION) return 0;
    size_t c = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += residue_int_bits(e->data.function.args[i]);
    return c;
}

/* Is candidate `s` a tidier spelling of `r`?  Strictly fewer leaves; or as
 * few while clearing a compound denominator (180/(-18 - 18 I Sqrt[3])^2 ->
 * (-5 - 5 I Sqrt[3])/72); or as few, introducing no compound denominator,
 * with smaller integers (a common content removed).  A Root object is
 * smaller but no tidier than the radical it replaces. */
static bool residue_tidier(Expr* s, Expr* r) {
    if (!s || residue_has_root(s)) return false;
    size_t ls = residue_leaves(s), lr = residue_leaves(r);
    if (ls != lr) return ls < lr;
    bool dr = residue_has_compound_denominator(r);
    bool ds = residue_has_compound_denominator(s);
    if (dr != ds) return dr;
    return residue_int_bits(s) < residue_int_bits(r);
}

/* An exact algebraic residue read off a Laurent series at an algebraic point
 * keeps the series' unreduced coefficient arithmetic (Residue[z^(1/3)/(1 +
 * z^2)^2, {z, I}] as a three-term sum of I^(k/3) quotients); Simplify, cheap
 * on a numeric constant, collapses it, and Together removes a common content
 * Simplify leaves behind ((54 - 3 I Sqrt[3])/432 -> (18 - I Sqrt[3])/144).
 * The tidiest candidate is kept (residue_tidier).  Consumes r. */
static Expr* residue_tidy_constant(Expr* r) {
    if (!r || r->type != EXPR_FUNCTION || !expr_is_numeric_quantity(r) ||
        internal_contains_inexact(r) || residue_leaves(r) <= 12)
        return r;
    static const char* const heads[] = { "Simplify", "Together" };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        Expr* s = residue_eval1(heads[i], expr_copy(r));
        if (residue_tidier(s, r)) { expr_free(r); r = s; }
        else if (s) expr_free(s);
    }
    return r;
}

Expr* builtin_residue(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return residue_emit_argcount(argc);
    if (argc != 2) return NULL;   /* only the two-argument form is handled */

    Expr* expr = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* The location spec must be List[z, z0] with z a symbol. */
    if (spec->type != EXPR_FUNCTION ||
        spec->data.function.head->type != EXPR_SYMBOL ||
        spec->data.function.head->data.symbol.name != SYM_List ||
        spec->data.function.arg_count != 2)
        return NULL;
    Expr* z  = spec->data.function.args[0];
    Expr* z0 = spec->data.function.args[1];
    if (z->type != EXPR_SYMBOL) return NULL;

    return residue_tidy_constant(residue_compute(expr, z, z0));
}

void residue_init(void) {
    symtab_add_builtin("Residue", builtin_residue);
    symtab_get_def("Residue")->attributes |= ATTR_PROTECTED;
}
