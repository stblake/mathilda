/* ============================================================================
 * limit_osc.c -- Oscillatory normal form for Limit at an infinite point.
 * ============================================================================
 *
 * MOTIVATION
 * ----------
 * Every other layer of limit.c treats an oscillation as an opaque "bounded
 * head": it can squeeze  Sin[g(x)] / x -> 0  because |Sin| <= 1, and it can
 * declare a *bare* Sin[g(x)] indeterminate, but it cannot see through a sum
 * of oscillations. So
 *
 *     Limit[(Cos[x^2]/x^2 - Cos[(x+1)^2]/(x+1)^2) x^3, x -> Infinity]
 *
 * -- asymptotically  2 x Sin[x^2+x+1/2] Sin[x+1/2] + O(1)  -- is out of reach:
 * the two crude |Cos| <= 1 envelopes are each of size x, so the squeeze fails,
 * and there is no single dominant term for the mrv/Series machinery to latch
 * onto. What is needed is a normal form in which distinct oscillations cannot
 * cancel, plus a criterion on that normal form.
 *
 * THE NORMAL FORM
 * ---------------
 * At x -> +/-Infinity write
 *
 *     f(x) = c_0(x) + SUM_{j in S} c_j(x) E^(I theta_j(x))                (NF)
 *
 * where the phases theta_j are pairwise-distinct real functions of x with no
 * constant term, and the amplitudes c_j are oscillation-free. TrigToExp
 * followed by Expand produces exactly this: each Sin/Cos becomes a pair of
 * exponentials, products distribute, and the terms group by exponent. Two
 * canonicalisation steps make the grouping a normal form rather than a
 * rewriting:
 *
 *   1. the exponent of each E-factor is Expand-ed, so E^(I (1+x)^2) and
 *      E^(I (1 + 2x + x^2)) land in the same group;
 *   2. the *constant* part of a phase is absorbed into the amplitude
 *      (E^(I (theta + c)) = E^(I c) E^(I theta)), so Cos[x] and Cos[x+1]
 *      share the phase x and their amplitudes combine.
 *
 * Step 2 also buys the key structural fact used below: distinct phases have a
 * *non-constant* difference.
 *
 * WHY (NF) IS A DECISION FORM
 * ---------------------------
 * Distinct phases are asymptotically orthogonal, so no cancellation can hide
 * between groups. Concretely, with gamma_j := lim |c_j|:
 *
 *  R1 (squeeze). If gamma_j = 0 for every j in S then |f - c_0| <= SUM |c_j|
 *     -> 0, hence lim f = lim c_0. Exact, and needs no hypothesis on the
 *     phases at all.
 *
 *  R3 (dominated oscillation). If SUM_{j in S} |c_j| / |c_0| -> r < 1 and
 *     c_0 -> +/-Infinity then |f - c_0| <= SUM |c_j| keeps the sign of f
 *     equal to that of c_0 and |f| >= (1-r-eps) |c_0| -> Infinity, so
 *     lim f = lim c_0. Note r < 1, not r = 0: x^2 (2 + Cos[x]) has r = 1/2
 *     and diverges, while x^2 (1 + Cos[x]) has r = 1 and does not.
 *
 *  R0 (strictly dominant oscillation; intermediate value theorem). Suppose one
 *     group d (together with its conjugate mate, if any) strictly dominates
 *     every other group, theta_d -> +/-Infinity continuously, arg(c_d) is
 *     bounded, and |c_d| does not tend to 0. Then
 *
 *         f = c_0 + A cos(theta_d + phi) + o(A),   A = |c_d| + |c_mate|,
 *
 *     with phi bounded, so theta_d + phi -> +/-Infinity continuously and by the
 *     IVT there are x_n, y_n -> Infinity with cos(...) = +1 and -1. Provided
 *     c_0 does not swamp A -- either lim c_0 is finite, or |c_0| < A eventually
 *     -- f has two distinct accumulation points (or accumulates at both
 *     +Infinity and -Infinity). No limit.
 *
 *  R2 (Weyl mean / mean-square). Suppose every theta_j (j in S) is a real
 *     polynomial in x of degree >= 1 with numeric coefficients. Every phase
 *     *difference* is then a non-constant polynomial (constant terms were
 *     stripped), so by van der Corput each cross term satisfies
 *     (1/T) INT_0^T c_j conj(c_k) E^(I (theta_j - theta_k)) = o(max_j (1/T)
 *     INT_0^T |c_j|^2) as long as the amplitudes are eventually monotone in
 *     modulus (true of every exp-log amplitude the rewrite produces) and grow
 *     at most polynomially. Hence
 *
 *         (1/T) INT_0^T |f|^2  ->  |lim c_0|^2 + SUM_{j in S} gamma_j^2,
 *         (1/T) INT_0^T f      ->  lim c_0.
 *
 *     If f -> L (finite) then the two Cesaro means converge to L and |L|^2, so
 *     L = lim c_0 and then SUM gamma_j^2 = 0. Contrapositive: if some gamma_j
 *     is non-zero, f has no finite limit. The divergent alternative f ->
 *     +/-Infinity is excluded either because f is bounded (all gamma_j finite
 *     and lim c_0 finite), or -- for real f -- by the window mean
 *     (1/T) INT_T^2T f, which integration by parts bounds by
 *     O(max_j |c_j(2T)| / T^deg theta_j) + lim c_0: bounded whenever
 *     |c_j| = O(x^deg theta_j), and therefore incompatible with f -> Infinity.
 *
 * Rules are tried in the order R1, R3, R0, R2 -- the two that produce a value
 * first, then the two that produce the Indeterminate verdict.
 *
 * WHAT THE MODULE REFUSES
 * -----------------------
 * Abstention (NULL) is the default whenever a hypothesis cannot be *verified*:
 * a phase that neither diverges nor is polynomial, an amplitude that still
 * carries an oscillation (Tan and Sec produce those -- TrigToExp leaves an
 * exponential in a denominator), a sub-limit the caller could not compute, a
 * modulus whose limit is not a decidable number, a realness certificate that
 * does not check out. In particular a symbolic amplitude such as `a Sin[x]`
 * abstains rather than claiming Indeterminate, because a = 0 is not excluded.
 *
 * FINITE LIMIT POINTS
 * -------------------
 * This entry point only accepts +/-Infinity. A finite point is handled by the
 * caller (`layer_oscillatory` in limit.c) through x = a +/- 1/t with
 * t -> +Infinity: an oscillation at a point is an oscillation at infinity in
 * t, with the identical normal form. A two-sided limit runs both sides and
 * requires them to agree.
 *
 * MEMORY
 * ------
 * Every helper returning Expr* returns a freshly-allocated tree owned by the
 * caller; every helper taking Expr* borrows it unless the name says _adopt.
 * The group table owns its phase/amp/sq trees and is released by groups_free.
 * ========================================================================= */

#include "limit_osc.h"

#include "expr.h"
#include "eval.h"
#include "common.h"
#include "sym_names.h"
#include "arithmetic.h"   /* is_infinity_sym */

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Resource guard. Expanding a product of k oscillatory factors costs 2^k
 * terms; past this many *groups* we abstain rather than grind. This bounds
 * work only -- it can turn an answer into an abstention, never into a wrong
 * answer, so it is not a cap on the decision procedure itself. */
#define OSC_MAX_GROUPS 512

/* ---------------------------------------------------------------------- */
/* Tiny builders / predicates (local copies; limit.c's are static)         */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v)     { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* n, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(mk_sym(n), args, 1);
}

static Expr* mk_fn2(const char* n, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(n), args, 2);
}

static Expr* mk_plus(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }
static Expr* mk_times(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }

/* Evaluate and free the source (evaluate() copies what it needs). */
static Expr* simp(Expr* e) {
    if (!e) return NULL;
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static Expr* expand_of(Expr* e_adopt) { return simp(mk_fn1("Expand", e_adopt)); }

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool is_lit_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    return false;
}

static bool is_neg_infinity(const Expr* e) {
    if (!head_is(e, SYM_Times) || e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (!is_infinity_sym(a) && !is_infinity_sym(b)) return false;
    const Expr* other = is_infinity_sym(a) ? b : a;
    if (other->type == EXPR_INTEGER) return other->data.integer < 0;
    if (other->type == EXPR_REAL)    return other->data.real < 0.0;
    return false;
}

/* Any flavour of infinity / undefined anywhere inside `e`. */
static bool has_divergence(const Expr* e) {
    if (!e) return true;
    if (is_infinity_sym((Expr*)e) || is_sym(e, "ComplexInfinity") ||
        is_sym(e, "Indeterminate") || is_sym(e, "Undefined") ||
        head_is(e, SYM_DirectedInfinity) || head_is(e, SYM_Interval)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (has_divergence(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_divergence(e->data.function.args[i])) return true;
    }
    return false;
}

static bool expr_has(const Expr* e, const Expr* target) {
    if (!e) return false;
    if (expr_eq((Expr*)e, (Expr*)target)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (expr_has(e->data.function.head, target)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_has(e->data.function.args[i], target)) return true;
    }
    return false;
}

static bool free_of(const Expr* e, const Expr* x) { return !expr_has(e, x); }

/* Ask the evaluator a yes/no question: does `q` (adopted) evaluate to True? */
static bool ask_true(Expr* q_adopt) {
    Expr* r = simp(q_adopt);
    bool yes = is_sym(r, "True");
    expr_free(r);
    return yes;
}

static bool possible_zero(const Expr* e) {
    return ask_true(mk_fn1("PossibleZeroQ", expr_copy((Expr*)e)));
}

/* PossibleZeroQ on a difference, freeing both operands (which are adopted). */
static bool same_value(Expr* a_adopt, Expr* b_adopt) {
    Expr* d = simp(mk_plus(a_adopt, mk_times(mk_int(-1), b_adopt)));
    bool eq = d && possible_zero(d);
    expr_free(d);
    return eq;
}

static bool numeric_q(const Expr* e) {
    return ask_true(mk_fn1("NumericQ", expr_copy((Expr*)e)));
}

/* A literal number (the only thing a fully-evaluated numeric factor product
 * can be): machine/bignum integer, real, MPFR, Rational[p,q] or Complex[a,b]. */
static bool is_number_leaf(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default: break;
    }
    return head_is(e, SYM_Rational) || head_is(e, SYM_Complex);
}

/* True iff `e` carries an explicit non-zero imaginary part anywhere. */
static bool has_imaginary(const Expr* e) {
    if (!e) return false;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2)
        return !is_lit_zero(e->data.function.args[1]);
    if (e->type == EXPR_FUNCTION) {
        if (has_imaginary(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_imaginary(e->data.function.args[i])) return true;
    }
    return false;
}

/* Complex conjugation *under the standing assumption that every symbol is
 * real* -- which is exactly the setting here, x being a real limit variable.
 * Negating every literal imaginary part is then conjugation, and it is exact
 * (unlike Conjugate[], which stays inert on symbolic arguments). */
static Expr* conj_real(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        return simp(mk_fn2("Complex", expr_copy(e->data.function.args[0]),
                           mk_times(mk_int(-1),
                                    expr_copy(e->data.function.args[1]))));
    }
    size_t n = e->data.function.arg_count;
    Expr** args = n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    if (n && !args) return expr_copy((Expr*)e);
    for (size_t i = 0; i < n; i++) args[i] = conj_real(e->data.function.args[i]);
    Expr* r = expr_new_function(conj_real(e->data.function.head), args, n);
    free(args);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Oscillation detection                                                   */
/* ---------------------------------------------------------------------- */

/* Heads whose value oscillates without decaying to a limit, so that leaving
 * one inside an "amplitude" would break the normal form. Sin/Cos are removed
 * by TrigToExp; Tan/Cot/Sec/Csc survive as exponential *quotients* and are
 * caught here; Bessel J/Y have no exponential rewrite at all. */
static bool is_oscillatory_head(const Expr* e) {
    return head_is(e, SYM_Sin)  || head_is(e, SYM_Cos)  ||
           head_is(e, SYM_Tan)  || head_is(e, SYM_Cot)  ||
           head_is(e, SYM_Sec)  || head_is(e, SYM_Csc)  ||
           head_is(e, SYM_BesselJ) || head_is(e, SYM_BesselY);
}

/* E^z (or Exp[z]); the exponent is reported through z_out. */
static bool is_exponential(const Expr* e, const Expr** z_out) {
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        is_sym(e->data.function.args[0], "E")) {
        if (z_out) *z_out = e->data.function.args[1];
        return true;
    }
    if (head_is(e, SYM_Exp) && e->data.function.arg_count == 1) {
        if (z_out) *z_out = e->data.function.args[0];
        return true;
    }
    return false;
}

/* Does `e` contain an oscillation in x? */
static bool has_oscillation(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_oscillatory_head(e)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_has(e->data.function.args[i], x)) return true;
    }
    const Expr* z = NULL;
    if (is_exponential(e, &z) && has_imaginary(z) && expr_has(z, x))
        return true;
    if (has_oscillation(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_oscillation(e->data.function.args[i], x)) return true;
    return false;
}

/* ---------------------------------------------------------------------- */
/* Group table                                                             */
/* ---------------------------------------------------------------------- */

typedef struct {
    Expr* phase;   /* theta_j: expanded, x-dependent, zero constant term      */
    Expr* amp;     /* c_j                                                     */
    Expr* sq;      /* |c_j|^2 = c_j conj(c_j), expanded                       */
    int   deg;     /* degree of phase as a numeric-coefficient polynomial, or
                    * -1 when it is not one                                   */
    int   mate;    /* index of the group whose phase is -theta_j, else -1     */
} OscGroup;

typedef struct {
    OscGroup* v;
    size_t    n, cap;
    int       zero;   /* index of the theta = 0 group, or -1                  */
} OscTable;

static void groups_free(OscTable* t) {
    for (size_t i = 0; i < t->n; i++) {
        expr_free(t->v[i].phase);
        expr_free(t->v[i].amp);
        expr_free(t->v[i].sq);
    }
    free(t->v);
    t->v = NULL; t->n = t->cap = 0;
}

/* Add `amp` (adopted) to the group keyed by `phase` (adopted on insert). */
static bool groups_add(OscTable* t, Expr* phase, Expr* amp) {
    for (size_t i = 0; i < t->n; i++) {
        if (expr_eq(t->v[i].phase, phase)) {
            t->v[i].amp = mk_plus(t->v[i].amp, amp);
            expr_free(phase);
            return true;
        }
    }
    if (t->n == OSC_MAX_GROUPS) { expr_free(phase); expr_free(amp); return false; }
    if (t->n == t->cap) {
        size_t nc = t->cap ? t->cap * 2 : 8;
        OscGroup* nv = (OscGroup*)realloc(t->v, nc * sizeof(OscGroup));
        if (!nv) { expr_free(phase); expr_free(amp); return false; }
        t->v = nv; t->cap = nc;
    }
    t->v[t->n].phase = phase;
    t->v[t->n].amp   = amp;
    t->v[t->n].sq    = NULL;
    t->v[t->n].deg   = -1;
    t->v[t->n].mate  = -1;
    t->n++;
    return true;
}

/* ---------------------------------------------------------------------- */
/* Term splitting: t  ->  amp * E^(I phase)                                */
/* ---------------------------------------------------------------------- */

/* A growable list of *borrowed* sub-expression pointers. */
typedef struct { Expr** v; size_t n, cap; } ExprList;

static void list_free(ExprList* L) { free(L->v); L->v = NULL; L->n = L->cap = 0; }

static bool list_push(ExprList* L, Expr* e) {
    if (L->n == L->cap) {
        size_t nc = L->cap ? L->cap * 2 : 8;
        Expr** nv = (Expr**)realloc(L->v, nc * sizeof(Expr*));
        if (!nv) return false;
        L->v = nv; L->cap = nc;
    }
    L->v[L->n++] = e;
    return true;
}

/* Collect the `sym`-parts of `e` (Times factors / Plus terms), descending
 * through *nested* same-head nodes. Flatness is an evaluator invariant that
 * an un-re-evaluated builtin result can violate -- TrigToExp, for one,
 * returns Times[c, Times[x, E^(I x)]] -- and a normal form that missed the
 * inner E-factor would silently mis-classify the whole limit. */
static bool collect_parts(Expr* e, const char* sym, ExprList* L) {
    if (head_is(e, sym)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!collect_parts(e->data.function.args[i], sym, L)) return false;
        return true;
    }
    return list_push(L, e);
}

static bool factors_of(Expr* e, ExprList* L) {
    return collect_parts(e, SYM_Times, L);
}

static bool terms_of(Expr* e, ExprList* L) {
    return collect_parts(e, SYM_Plus, L);
}

/* Split `u` into (numeric factor, remaining factor). Both are fresh. */
static void split_numeric_coefficient(Expr* u, Expr** num, Expr** rest) {
    if (is_number_leaf(u)) { *num = expr_copy(u); *rest = mk_int(1); return; }
    ExprList fs = {0};
    if (!factors_of(u, &fs)) { list_free(&fs); *num = mk_int(1); *rest = expr_copy(u); return; }
    Expr* nm = mk_int(1);
    Expr* rs = mk_int(1);
    for (size_t i = 0; i < fs.n; i++) {
        if (is_number_leaf(fs.v[i])) nm = mk_times(nm, expr_copy(fs.v[i]));
        else                         rs = mk_times(rs, expr_copy(fs.v[i]));
    }
    list_free(&fs);
    *num  = simp(nm);
    *rest = simp(rs);
}

/* Real and imaginary parts of an evaluated *literal* number. Borrowed out. */
static void number_re_im(Expr* num, Expr** re, Expr** im) {
    if (head_is(num, SYM_Complex) && num->data.function.arg_count == 2) {
        *re = expr_copy(num->data.function.args[0]);
        *im = expr_copy(num->data.function.args[1]);
    } else {
        *re = expr_copy(num);
        *im = mk_int(0);
    }
}

/* Decompose one expanded product `t` of the TrigToExp/Expand output into an
 * amplitude and a real phase with no constant term.
 *
 * Returns false when `t` is outside the class -- an exponent whose numeric
 * coefficient multiplies something that itself carries an imaginary part, for
 * instance. On success *out_amp and *out_phase are fresh and owned. */
static bool split_term(Expr* t, Expr* x, Expr** out_amp, Expr** out_phase) {
    ExprList fs = {0};
    if (!factors_of(t, &fs)) { list_free(&fs); return false; }

    Expr* amp  = mk_int(1);
    Expr* zsum = mk_int(0);
    for (size_t i = 0; i < fs.n; i++) {
        const Expr* z = NULL;
        if (is_exponential(fs.v[i], &z))
            zsum = mk_plus(zsum, expr_copy((Expr*)z));
        else
            amp = mk_times(amp, expr_copy(fs.v[i]));
    }
    list_free(&fs);
    zsum = expand_of(zsum);
    if (!zsum) { expr_free(amp); return false; }

    /* zsum = rho + I theta with rho, theta real. */
    Expr* rho   = mk_int(0);
    Expr* theta = mk_int(0);
    ExprList us = {0};
    bool ok = terms_of(zsum, &us);
    for (size_t i = 0; i < us.n && ok; i++) {
        Expr *num = NULL, *rest = NULL;
        split_numeric_coefficient(us.v[i], &num, &rest);
        if (!is_number_leaf(num) || has_imaginary(rest)) {
            ok = false;
        } else {
            Expr *re = NULL, *im = NULL;
            number_re_im(num, &re, &im);
            rho   = mk_plus(rho,   mk_times(re, expr_copy(rest)));
            theta = mk_plus(theta, mk_times(im, expr_copy(rest)));
        }
        expr_free(num);
        expr_free(rest);
    }
    list_free(&us);
    expr_free(zsum);
    if (!ok) { expr_free(amp); expr_free(rho); expr_free(theta); return false; }

    rho   = simp(rho);
    theta = expand_of(theta);
    if (!rho || !theta) { expr_free(amp); expr_free(rho); expr_free(theta); return false; }

    /* The real part of the exponent belongs to the amplitude. */
    if (!is_lit_zero(rho)) amp = mk_times(amp, mk_fn2("Power", mk_sym("E"), rho));
    else                   expr_free(rho);

    /* Strip the x-free part of the phase into the amplitude, so that phases
     * differing by a constant share a group and distinct phases differ by a
     * non-constant function. */
    Expr* tc = mk_int(0);
    Expr* tv = mk_int(0);
    ExprList th = {0};
    if (!terms_of(theta, &th)) {
        list_free(&th); expr_free(theta); expr_free(amp);
        expr_free(tc); expr_free(tv); return false;
    }
    for (size_t i = 0; i < th.n; i++) {
        if (free_of(th.v[i], x)) tc = mk_plus(tc, expr_copy(th.v[i]));
        else                     tv = mk_plus(tv, expr_copy(th.v[i]));
    }
    list_free(&th);
    expr_free(theta);
    tc = simp(tc);
    tv = expand_of(tv);
    if (!tc || !tv) { expr_free(amp); expr_free(tc); expr_free(tv); return false; }

    if (!is_lit_zero(tc)) {
        amp = mk_times(amp, mk_fn2("Power", mk_sym("E"),
                                   mk_times(mk_fn2("Complex", mk_int(0), mk_int(1)),
                                            tc)));
    } else {
        expr_free(tc);
    }

    *out_amp   = simp(amp);
    *out_phase = tv;
    if (!*out_amp) { expr_free(tv); return false; }
    return true;
}

/* ---------------------------------------------------------------------- */
/* Building the normal form                                                */
/* ---------------------------------------------------------------------- */

/* Returns false when f is outside the class (caller abstains). */
static bool build_normal_form(Expr* f, Expr* x, OscTable* t) {
    t->v = NULL; t->n = t->cap = 0; t->zero = -1;

    Expr* g = simp(mk_fn1("TrigToExp", expr_copy(f)));
    if (!g) return false;
    g = expand_of(g);
    if (!g) return false;

    ExprList ts = {0};
    bool ok = terms_of(g, &ts);
    for (size_t i = 0; i < ts.n && ok; i++) {
        Expr *amp = NULL, *phase = NULL;
        if (!split_term(ts.v[i], x, &amp, &phase)) { ok = false; break; }
        ok = groups_add(t, phase, amp);
    }
    list_free(&ts);
    expr_free(g);
    if (!ok) { groups_free(t); return false; }

    /* Normalise the accumulated amplitudes. */
    for (size_t i = 0; i < t->n; i++) {
        t->v[i].amp = simp(t->v[i].amp);
        if (!t->v[i].amp) { groups_free(t); return false; }
    }

    /* Every amplitude must be oscillation-free, or the normal form is not
     * complete and no verdict is available. This has to run BEFORE the
     * zero-prune below: PossibleZeroQ samples numerically, so an amplitude
     * like (2 Sin[t])^(t^2) -- which underflows wherever |2 Sin t| < 1 --
     * reads as zero, and dropping it would leave an empty normal form and a
     * confident (wrong) limit. Purity is structural and cannot be fooled. */
    for (size_t i = 0; i < t->n; i++) {
        if (has_oscillation(t->v[i].amp, x) ||
            has_oscillation(t->v[i].phase, x) ||
            has_imaginary(t->v[i].phase)) {
            groups_free(t);
            return false;
        }
    }

    /* Drop the groups whose amplitudes cancelled (Cos[x]^2 + Sin[x]^2 kills
     * both E^(+/-2 I x) groups). Iterating downwards keeps indices stable.
     * A PossibleZeroQ false positive on a *decaying* amplitude is harmless
     * here: dropping a group with |c_j| -> 0 is exactly what rule R1 would
     * have done with it. */
    for (size_t i = t->n; i > 0; i--) {
        size_t k = i - 1;
        if (possible_zero(t->v[k].amp)) {
            expr_free(t->v[k].phase);
            expr_free(t->v[k].amp);
            memmove(&t->v[k], &t->v[k + 1], (t->n - k - 1) * sizeof(OscGroup));
            t->n--;
        }
    }

    /* Locate the non-oscillatory group and the conjugate mates. */
    for (size_t i = 0; i < t->n; i++)
        if (is_lit_zero(t->v[i].phase)) t->zero = (int)i;

    for (size_t i = 0; i < t->n; i++) {
        if ((int)i == t->zero) continue;
        for (size_t k = 0; k < t->n; k++) {
            if (k == i || (int)k == t->zero) continue;
            Expr* s = simp(mk_plus(expr_copy(t->v[i].phase),
                                   expr_copy(t->v[k].phase)));
            bool opposite = is_lit_zero(s);
            expr_free(s);
            if (opposite) { t->v[i].mate = (int)k; break; }
        }
    }

    /* |c|^2 = c conj(c) -- always real, and far friendlier to the limit
     * engine than Abs[], which stays inert on symbolic arguments. */
    for (size_t i = 0; i < t->n; i++) {
        Expr* c = t->v[i].amp;
        t->v[i].sq = expand_of(mk_times(expr_copy(c), conj_real(c)));
        if (!t->v[i].sq) { groups_free(t); return false; }
    }

    /* Phase degree, for the R2 gate: a polynomial in x of degree >= 1 whose
     * coefficients are all numeric. */
    for (size_t i = 0; i < t->n; i++) {
        if ((int)i == t->zero) continue;
        Expr* th = t->v[i].phase;
        if (!ask_true(mk_fn2("PolynomialQ", expr_copy(th), expr_copy(x)))) continue;
        Expr* d = simp(mk_fn2("Exponent", expr_copy(th), expr_copy(x)));
        if (d && d->type == EXPR_INTEGER && d->data.integer >= 1 &&
            d->data.integer < 1024) {
            int deg = (int)d->data.integer;
            bool numeric = true;
            for (int k = 0; k <= deg && numeric; k++) {
                Expr* cargs[3] = { expr_copy(th), expr_copy(x), mk_int(k) };
                Expr* ck = simp(expr_new_function(mk_sym("Coefficient"), cargs, 3));
                numeric = ck && numeric_q(ck);
                expr_free(ck);
            }
            if (numeric) t->v[i].deg = deg;
        }
        expr_free(d);
    }
    return true;
}

/* ---------------------------------------------------------------------- */
/* Classification of the amplitude moduli                                  */
/* ---------------------------------------------------------------------- */

typedef enum { SQ_ZERO, SQ_POSITIVE, SQ_INFINITE, SQ_UNKNOWN } SqClass;

/* lim |c|^2, coarsened to the three cases the rules distinguish. Anything
 * that is not a decidable number (a symbolic parameter, say) is UNKNOWN and
 * makes the caller abstain -- `a Sin[x]` must not be declared Indeterminate,
 * because a = 0 has the limit 0. */
static SqClass classify_sq(Expr* sq, Expr* x, LimitOscSubFn sub, void* ctx) {
    Expr* v = sub(sq, ctx);
    if (!v) return SQ_UNKNOWN;
    SqClass c = SQ_UNKNOWN;
    if (is_lit_zero(v))          c = SQ_ZERO;
    else if (is_infinity_sym(v)) c = SQ_INFINITE;
    else if (!has_divergence(v) && free_of(v, x) && numeric_q(v))
        c = possible_zero(v) ? SQ_ZERO : SQ_POSITIVE;
    expr_free(v);
    return c;
}

/* lim (a / b) as a decision on "a is negligible against b". */
static bool ratio_is_zero(Expr* a, Expr* b, LimitOscSubFn sub, void* ctx) {
    Expr* q = simp(mk_times(expr_copy(a), mk_fn2("Power", expr_copy(b), mk_int(-1))));
    Expr* v = sub(q, ctx);
    expr_free(q);
    bool zero = v && is_lit_zero(v);
    expr_free(v);
    return zero;
}

/* lim (a / b) exists, is a real number and is < bound. */
static bool ratio_below(Expr* a, Expr* b, int64_t bound, Expr* x,
                        LimitOscSubFn sub, void* ctx) {
    Expr* q = simp(mk_times(expr_copy(a), mk_fn2("Power", expr_copy(b), mk_int(-1))));
    Expr* v = sub(q, ctx);
    expr_free(q);
    if (!v) return false;
    bool ok = false;
    if (!has_divergence(v) && free_of(v, x) && numeric_q(v))
        ok = ask_true(mk_fn2("Less", expr_copy(v), mk_int(bound)));
    expr_free(v);
    return ok;
}

/* Is `v` a finite value -- a limit that exists and is not an infinity? */
static bool is_finite_value(Expr* v, Expr* x) {
    return v && !has_divergence(v) && free_of(v, x);
}

/* arg(c) is bounded: c splits as (x-free factor) * (factor with no explicit
 * imaginary part). Sufficient, cheap, and true of every amplitude TrigToExp
 * produces from a real integrand. */
static bool amplitude_arg_bounded(Expr* c, Expr* x) {
    ExprList fs = {0};
    if (!factors_of(c, &fs)) { list_free(&fs); return false; }
    bool ok = true;
    for (size_t i = 0; i < fs.n && ok; i++)
        if (!free_of(fs.v[i], x) && has_imaginary(fs.v[i])) ok = false;
    list_free(&fs);
    return ok;
}

/* The normal form describes a real-valued function: the theta = 0 amplitude
 * is real and every other group is paired with its conjugate. */
static bool normal_form_is_real(OscTable* t) {
    for (size_t i = 0; i < t->n; i++) {
        Expr* c = t->v[i].amp;
        if ((int)i != t->zero && t->v[i].mate < 0) return false;
        Expr* peer = ((int)i == t->zero) ? c : t->v[t->v[i].mate].amp;
        if (!same_value(expr_copy(peer), conj_real(c))) return false;
    }
    return true;
}

/* ---------------------------------------------------------------------- */
/* The layer                                                               */
/* ---------------------------------------------------------------------- */

Expr* limit_oscillatory(Expr* f, Expr* x, Expr* point,
                        LimitOscSubFn sub, void* subctx) {
    if (!f || !x || !point || !sub) return NULL;
    if (!is_infinity_sym(point) && !is_neg_infinity(point)) return NULL;
    if (x->type != EXPR_SYMBOL) return NULL;
    if (!has_oscillation(f, x)) return NULL;

    OscTable t;
    if (!build_normal_form(f, x, &t)) return NULL;

    /* c_0 and its limit. A missing theta = 0 group means c_0 = 0. */
    Expr* c0 = (t.zero >= 0) ? expr_copy(t.v[t.zero].amp) : mk_int(0);
    Expr* L0 = sub(c0, subctx);

    /* Nothing oscillatory survived (Cos[x]^2 + Sin[x]^2, Sin[x + 2 Pi] -
     * Sin[x], ...): the limit is simply that of what is left. */
    bool any_osc = false;
    for (size_t i = 0; i < t.n; i++) if ((int)i != t.zero) any_osc = true;
    if (!any_osc) {
        expr_free(c0);
        groups_free(&t);
        if (L0 && !is_sym(L0, "Indeterminate")) return L0;
        expr_free(L0);
        return NULL;
    }

    /* Classify every oscillatory modulus. */
    SqClass* cls = (SqClass*)malloc(t.n * sizeof(SqClass));
    if (!cls) { expr_free(c0); expr_free(L0); groups_free(&t); return NULL; }
    bool all_known = true, all_zero = true, all_finite = true;
    for (size_t i = 0; i < t.n; i++) {
        if ((int)i == t.zero) { cls[i] = SQ_UNKNOWN; continue; }
        cls[i] = classify_sq(t.v[i].sq, x, sub, subctx);
        if (cls[i] == SQ_UNKNOWN)  all_known  = false;
        if (cls[i] != SQ_ZERO)     all_zero   = false;
        if (cls[i] == SQ_INFINITE) all_finite = false;
    }

    Expr* result = NULL;

    /* ---- R1: every oscillation decays. |f - c_0| <= SUM |c_j| -> 0. ---- */
    if (!result && all_known && all_zero && L0 && free_of(L0, x) &&
        !is_sym(L0, "Indeterminate") && !is_sym(L0, "ComplexInfinity")) {
        result = L0;
        L0 = NULL;
    }

    /* ---- R3: the non-oscillatory part swamps every oscillation. --------
     * |f - c_0| <= SUM_j |c_j|, so when SUM_j |c_j| / |c_0| tends to some
     * r < 1 the sign of f is eventually that of c_0 and |f| >= (1-r-eps)
     * |c_0| -> Infinity. Working through lim |c_j|^2/|c_0|^2 keeps every
     * comparison a plain limit of a ratio; the square roots are then taken
     * on *numbers*, which is why r = 1/2 is decided for x^2 (2 + Cos[x])
     * where a term-by-term "o(c_0)" test would abstain. */
    if (!result && all_known && t.zero >= 0 && L0 &&
        (is_infinity_sym(L0) || is_neg_infinity(L0))) {
        Expr* total = mk_int(0);
        bool ok = true;
        for (size_t i = 0; i < t.n && ok; i++) {
            if ((int)i == t.zero) continue;
            Expr* q = simp(mk_times(expr_copy(t.v[i].sq),
                                    mk_fn2("Power", expr_copy(t.v[t.zero].sq),
                                           mk_int(-1))));
            Expr* v = sub(q, subctx);
            expr_free(q);
            if (!v || has_divergence(v) || !free_of(v, x) || !numeric_q(v)) ok = false;
            else total = mk_plus(total, mk_fn1("Sqrt", expr_copy(v)));
            expr_free(v);
        }
        total = simp(total);
        if (ok && total && ask_true(mk_fn2("Less", expr_copy(total), mk_int(1)))) {
            result = L0;
            L0 = NULL;
        }
        expr_free(total);
    }

    /* ---- R0: a strictly dominant oscillation; IVT gives two limits. ---- */
    if (!result && all_known) {
        for (size_t d = 0; d < t.n && !result; d++) {
            if ((int)d == t.zero) continue;
            if (cls[d] != SQ_POSITIVE && cls[d] != SQ_INFINITE) continue;
            if (!amplitude_arg_bounded(t.v[d].amp, x)) continue;

            /* The phase must sweep continuously to +/-Infinity. */
            Expr* pv = sub(t.v[d].phase, subctx);
            bool phase_diverges = pv && (is_infinity_sym(pv) || is_neg_infinity(pv));
            expr_free(pv);
            if (!phase_diverges) continue;

            /* Strict dominance over every other oscillatory group. */
            int mate = t.v[d].mate;
            bool dominant = true;
            for (size_t k = 0; k < t.n && dominant; k++) {
                if ((int)k == t.zero || k == d || (int)k == mate) continue;
                dominant = ratio_is_zero(t.v[k].sq, t.v[d].sq, sub, subctx);
            }
            if (!dominant) continue;
            /* A mate must be a genuine conjugate, or the "2|c| cos(theta+phi)"
             * reduction is not available. */
            if (mate >= 0 &&
                !same_value(expr_copy(t.v[mate].amp), conj_real(t.v[d].amp)))
                continue;
            /* Unmated *and* unbounded means a complex-valued f whose modulus
             * diverges -- ComplexInfinity, not Indeterminate, is arguably the
             * answer there (x E^(I x)), so abstain rather than pick. An
             * unmated bounded oscillation (E^(I x)) is safe: |f| stays
             * bounded, so no infinity is on the table. */
            if (mate < 0 && cls[d] == SQ_INFINITE) continue;

            /* c_0 must not swamp the oscillation: either it converges (then
             * the two accumulation points differ by 2A regardless), or its
             * modulus stays below the oscillation envelope A = (1 + [mate])|c_d|. */
            bool c0_ok = is_finite_value(L0, x);
            if (!c0_ok && t.zero >= 0)
                c0_ok = ratio_below(t.v[t.zero].sq, t.v[d].sq,
                                    mate >= 0 ? 4 : 1, x, sub, subctx);
            if (!c0_ok && t.zero < 0) c0_ok = true;
            if (!c0_ok) continue;

            result = mk_sym("Indeterminate");
        }
    }

    /* ---- R2: Weyl mean / mean-square over polynomial phases. ----------- */
    if (!result && all_known && !all_zero && is_finite_value(L0, x)) {
        bool poly = true;
        for (size_t i = 0; i < t.n && poly; i++) {
            if ((int)i == t.zero) continue;
            if (t.v[i].deg < 1) poly = false;
        }
        if (poly) {
            /* f bounded: nothing to exclude on the +/-Infinity side. */
            bool no_infinity = all_finite;
            /* Otherwise: real-valued, with |c_j| = O(x^deg theta_j), so the
             * window mean (1/T) INT_T^2T f stays bounded. */
            if (!no_infinity && normal_form_is_real(&t)) {
                no_infinity = true;
                for (size_t i = 0; i < t.n && no_infinity; i++) {
                    if ((int)i == t.zero) continue;
                    Expr* xp = simp(mk_fn2("Power", expr_copy(x),
                                           mk_int(2 * t.v[i].deg)));
                    Expr* q  = simp(mk_times(expr_copy(t.v[i].sq),
                                             mk_fn2("Power", xp, mk_int(-1))));
                    Expr* v  = sub(q, subctx);
                    expr_free(q);
                    no_infinity = v && !has_divergence(v) && numeric_q(v);
                    expr_free(v);
                }
            }
            if (no_infinity) result = mk_sym("Indeterminate");
        }
    }

    free(cls);
    expr_free(c0);
    expr_free(L0);
    groups_free(&t);
    return result;
}
