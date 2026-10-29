/* integrate_newton_leibniz.c
 *
 * Definite integration by the Newton-Leibniz rule (the fundamental theorem
 * of calculus).  See integrate_newton_leibniz.h for the overview.
 *
 * Pipeline for Integrate[f, {x, a, b}]:
 *
 *   1. F = Integrate[f, x]            -- reuse the indefinite cascade.
 *      If F still contains an unevaluated Integrate, bail (NULL): we never
 *      fabricate a value for an integrand we cannot antidifferentiate.
 *
 *   2. Locate the real singular points of f and F strictly inside (a, b)
 *      (poles of the rational part -- denominator roots).  If any pole's
 *      position is undecidable (symbolic coefficients / bounds), bail (NULL)
 *      rather than risk a wrong answer.
 *
 *   3. Split [a, b] at the interior poles into segments [p_i, p_{i+1}] and
 *      form the telescoping sum
 *          Sum_i ( F(p_{i+1}^-) - F(p_i^+) ).
 *      Each boundary is evaluated with the Limit engine: a plain limit at
 *      infinite endpoints, a one-sided limit ("FromBelow"/"FromAbove") at
 *      interior poles, and direct substitution at ordinary finite endpoints
 *      (falling back to a one-sided limit when substitution is singular).
 *      This is where improper integrals acquire their correct finite value
 *      and genuinely divergent ones surface Infinity / ComplexInfinity /
 *      Indeterminate.
 *
 * Branch-form antiderivatives.  Many continuous integrands (e.g. 1/(2+Cos[x]))
 * antidifferentiate to an ALREADY-CONTINUOUS branch form: the Weierstrass
 * integrator emits F = ArcTan[Tan[x/2]...] + K Floor[(x-b)/(2Pi)], where the
 * Floor step exactly cancels the ArcTan jump so F is continuous on (a, b).  For
 * such F, F(b)-F(a) by plain substitution is already the correct value; there
 * is no interior discontinuity to split at (and the Limit engine could not take
 * a one-sided limit through Floor anyway).  We therefore (i) detect poles on the
 * integrand f ONLY -- a continuous f has an antiderivative with no poles of its
 * own -- and (ii) guard the result: antiderivatives that could carry a genuine
 * UNCORRECTED branch jump (step heads, or an inverse-trig node sitting over a
 * pole-bearing trig head) are cross-checked against NIntegrate and bail rather
 * than return a wrong number.  A genuinely divergent integral emits
 * Integrate::idiv and is left unevaluated (Mathematica-faithful).
 */

#include "integrate_newton_leibniz.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "print.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* -------------------------------------------------------------------------
 * Small expression-construction / evaluation helpers.
 * ---------------------------------------------------------------------- */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c }, 3);
}

/* Evaluate `call`, free the call expression, return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True iff `e` is the compound `name[...]` (by head name, not interned ptr). */
static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* True iff any subexpression of `e` is a call with head `name`. */
static bool contains_head(const Expr* e, const char* name) {
    if (!e) return false;
    if (head_name_is(e, name)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_head(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head(e->data.function.args[i], name)) return true;
    return false;
}

/* True iff the symbol `x` occurs anywhere in `e`. */
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == x->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff some subexpression is `heads[i][...]` (any i) AND contains `x`. */
static bool contains_any_head_over_x(const Expr* e, const Expr* x,
                                     const char* const* heads, size_t nh) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < nh; i++)
        if (head_name_is(e, heads[i]) && contains_symbol(e, x)) return true;
    if (contains_any_head_over_x(e->data.function.head, x, heads, nh)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_any_head_over_x(e->data.function.args[i], x, heads, nh))
            return true;
    return false;
}

/* Step / piecewise heads over x: the pieces the Floor continuity-correction is
 * built from, and what the Limit engine cannot see through. */
static bool nl_has_step_head(const Expr* F, const Expr* x) {
    static const char* const H[] = { "Floor", "Ceiling", "Round", "Sign",
                                     "UnitStep", "Piecewise", "Mod", "Quotient" };
    return contains_any_head_over_x(F, x, H, sizeof H / sizeof *H);
}

/* An inverse-circular/hyperbolic node whose subtree contains a pole-bearing
 * trig head over x (e.g. ArcTan[Tan[x/2]/Sqrt[3]]).  This is the fingerprint of
 * a branch antiderivative whose jump may be uncorrected (no accompanying step
 * term), so it must be validated before we trust F(b)-F(a). */
static bool nl_has_arctrig_over_pole(const Expr* e, const Expr* x) {
    static const char* const INV[]  = { "ArcTan", "ArcCot", "ArcTanh", "ArcCoth" };
    static const char* const POLE[] = { "Tan", "Cot", "Sec", "Csc" };
    if (!e || e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < sizeof INV / sizeof *INV; i++)
        if (head_name_is(e, INV[i]) &&
            contains_any_head_over_x(e, x, POLE, sizeof POLE / sizeof *POLE))
            return true;
    if (nl_has_arctrig_over_pole(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_has_arctrig_over_pole(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` contains a Power with a negative exponent whose base depends on
 * x: a denominator that can vanish on the real line. */
static bool nl_has_x_denominator(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        bool neg = (ex->type == EXPR_INTEGER && ex->data.integer < 0) ||
                   (head_name_is(ex, "Rational") && ex->data.function.arg_count == 2 &&
                    ex->data.function.args[0]->type == EXPR_INTEGER &&
                    ex->data.function.args[0]->data.integer < 0);
        if (neg && contains_symbol(e->data.function.args[0], x)) return true;
    }
    if (nl_has_x_denominator(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_has_x_denominator(e->data.function.args[i], x)) return true;
    return false;
}

/* An ArcTan / ArcCot of a quotient over x, e.g. -ArcTan[Sqrt[1 - x^2]/x]: the
 * form a complex logarithm takes in real form.  It jumps by Pi wherever the
 * quotient's denominator crosses zero, so F(b) - F(a) is wrong when that
 * happens inside (a, b) -- the antiderivative is only piecewise continuous.
 * Same fingerprint class as nl_has_arctrig_over_pole: validate before trusting. */
static bool nl_has_arctan_of_quotient(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if ((head_name_is(e, "ArcTan") || head_name_is(e, "ArcCot")) &&
        e->data.function.arg_count == 1 &&
        nl_has_x_denominator(e->data.function.args[0], x))
        return true;
    if (nl_has_arctan_of_quotient(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_has_arctan_of_quotient(e->data.function.args[i], x)) return true;
    return false;
}

/* -------------------------------------------------------------------------
 * Numeric classification of bounds / roots.
 * ---------------------------------------------------------------------- */

/* Extract a machine double from an already-numeric leaf.  Returns false for
 * anything that is not a concrete real number. */
static bool nl_real_double(const Expr* e, double* out) {
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d) && d != 0) {
        *out = (double)n / (double)d;
        return true;
    }
    return false;
}

/* Classify `e` numerically via N[e].
 *   0 = not a real number (symbolic, or genuinely complex);
 *   1 = finite real  (value in *val);
 *   2 = +Infinity;
 *   3 = -Infinity.                                                        */
static int nl_numeric(Expr* e, double* val) {
    if (is_infinity_sym(e))      return 2;
    if (is_neg_infinity_form(e)) return 3;

    Expr* n = eval_take(mk_fn1("N", expr_copy(e)));
    if (!n) return 0;

    int r = 0;
    if (is_infinity_sym(n))           r = 2;
    else if (is_neg_infinity_form(n)) r = 3;
    else {
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex(n, &re, &im)) {
            /* On the real axis only when the imaginary part vanishes. */
            if (expr_numeric_sign(im) == 0 && nl_real_double(re, val)) r = 1;
        } else if (nl_real_double(n, val)) {
            r = 1;
        }
    }
    expr_free(n);
    return r;
}

/* Interval classification of a finite real root `rv` against [a, b]. */
typedef enum { NL_EXTERIOR = 0, NL_INTERIOR = 1, NL_UNKNOWN = 2 } NLClass;

/* Small relative tolerance for endpoint coincidence. */
static bool nl_close(double u, double v) {
    double s = fabs(u) + fabs(v) + 1.0;
    return fabs(u - v) <= 1e-12 * s;
}

/* Is finite `rv` strictly greater than the bound classified by (kind, val)?
 * kind: 1 finite (value in val), 2 +Infinity, 3 -Infinity.  A finite rv is
 * always above -Infinity and never above +Infinity. */
static bool nl_rv_above(double rv, int kind, double val) {
    if (kind == 3) return true;     /* rv > -Infinity */
    if (kind == 2) return false;    /* rv > +Infinity is impossible */
    return (rv > val) && !nl_close(rv, val);
}

/* Is finite `rv` strictly less than the bound (kind, val)? */
static bool nl_rv_below(double rv, int kind, double val) {
    if (kind == 3) return false;    /* rv < -Infinity is impossible */
    if (kind == 2) return true;     /* rv < +Infinity */
    return (rv < val) && !nl_close(rv, val);
}

/* Classify a finite real root `rv` against the integration endpoints a, b.
 * Order-INDEPENDENT: the pole is interior iff it lies strictly between a and b
 * regardless of whether a < b or a > b (reversed limits), so a divergent
 * integral is caught for both orientations rather than silently producing a
 * wrong (branch-form) value for the reversed case. */
static NLClass nl_classify(double rv, Expr* a, Expr* b) {
    double av, bv;
    int ak = nl_numeric(a, &av);
    int bk = nl_numeric(b, &bv);
    if (ak == 0 || bk == 0) return NL_UNKNOWN;   /* symbolic bound */

    /* Strictly between a and b in either orientation. */
    bool between = (nl_rv_above(rv, ak, av) && nl_rv_below(rv, bk, bv)) ||
                   (nl_rv_below(rv, ak, av) && nl_rv_above(rv, bk, bv));
    return between ? NL_INTERIOR : NL_EXTERIOR;
}

/* True iff Refine[pred, as] proves `pred` to the symbol True.  Consumes pred.
 * Quiet: an internal classification probe; a Solve::nongen note from a
 * parametric inequality must not surface from Integrate. */
static bool nl_refine_true(Expr* pred, Expr* as) {
    Expr* inner = as ? mk_fn2("Refine", pred, expr_copy(as)) : mk_fn1("Refine", pred);
    Expr* r = eval_take(mk_fn1("Quiet", inner));
    bool t = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    if (r) expr_free(r);
    return t;
}

/* Assumption-aware classification of a (possibly symbolic) root EXPRESSION `val`
 * against the endpoints a, b, used when a numeric decision is impossible (a
 * symbolic root like x = a, or a symbolic bound) and Assumptions are supplied.
 * INTERIOR iff provably strictly between a and b (either orientation); EXTERIOR
 * iff provably at-or-outside an endpoint (so an integration endpoint x = a is
 * EXTERIOR, not a split point); otherwise UNKNOWN. */
static NLClass nl_classify_symbolic(Expr* val, Expr* a, Expr* b, Expr* as) {
    if (!as) return NL_UNKNOWN;
    Expr* in = mk_fn2("Or",
        mk_fn2("And", mk_fn2("Less", expr_copy(a), expr_copy(val)),
                      mk_fn2("Less", expr_copy(val), expr_copy(b))),
        mk_fn2("And", mk_fn2("Less", expr_copy(b), expr_copy(val)),
                      mk_fn2("Less", expr_copy(val), expr_copy(a))));
    if (nl_refine_true(in, as)) return NL_INTERIOR;
    Expr* ext = mk_fn2("Or",
        mk_fn2("And", mk_fn2("LessEqual", expr_copy(val), expr_copy(a)),
                      mk_fn2("LessEqual", expr_copy(val), expr_copy(b))),
        mk_fn2("And", mk_fn2("GreaterEqual", expr_copy(val), expr_copy(a)),
                      mk_fn2("GreaterEqual", expr_copy(val), expr_copy(b))));
    if (nl_refine_true(ext, as)) return NL_EXTERIOR;
    return NL_UNKNOWN;
}

/* -------------------------------------------------------------------------
 * Interior-pole detection.
 * ---------------------------------------------------------------------- */

typedef struct {
    Expr** pts;   /* owned copies of the interior pole expressions           */
    double* vals; /* their numeric values, used for ordering / dedup         */
    size_t n, cap;
    bool undecidable; /* a pole exists whose interior-ness we cannot decide  */
    /* The set holds a single SYMBOLIC interior pole (x = a, proved interior
     * under the assumptions): its vals[] entry is NaN.  One such pole needs no
     * ordering -- the breakpoints are [lo, pole, hi] -- so it is placeable;
     * any second pole alongside it cannot be ordered and is undecidable. */
    bool symbolic;
} PoleSet;

static void poleset_init(PoleSet* s) {
    s->pts = NULL; s->vals = NULL; s->n = 0; s->cap = 0; s->undecidable = false;
    s->symbolic = false;
}

static void poleset_free(PoleSet* s) {
    for (size_t i = 0; i < s->n; i++) expr_free(s->pts[i]);
    free(s->pts); free(s->vals);
    poleset_init(s);
}

/* Insert `val` (an owned Expr, numeric value `v`) keeping the array sorted
 * ascending and free of duplicates.  Consumes `val` (frees it on dedup). */
static void poleset_insert(PoleSet* s, Expr* val, double v) {
    if (s->symbolic) {                  /* cannot order against a symbolic pole */
        expr_free(val); s->undecidable = true; return;
    }
    for (size_t i = 0; i < s->n; i++) {
        if (nl_close(s->vals[i], v)) { expr_free(val); return; }
    }
    if (s->n == s->cap) {
        size_t nc = s->cap ? s->cap * 2 : 4;
        s->pts  = realloc(s->pts,  nc * sizeof(*s->pts));
        s->vals = realloc(s->vals, nc * sizeof(*s->vals));
        s->cap = nc;
    }
    /* Find sorted insertion point. */
    size_t k = s->n;
    while (k > 0 && s->vals[k - 1] > v) {
        s->pts[k]  = s->pts[k - 1];
        s->vals[k] = s->vals[k - 1];
        k--;
    }
    s->pts[k]  = val;
    s->vals[k] = v;
    s->n++;
}

/* The RHS of the (single) `Rule[x, val]` inside a Solve solution element,
 * or NULL if `el` is not of that shape.  Returns a borrowed pointer. */
static Expr* nl_rule_rhs_for(Expr* el, const Expr* x) {
    /* el may be `List[Rule[x, val], ...]` (univariate Solve) or a bare Rule. */
    if (head_name_is(el, "Rule") && el->data.function.arg_count == 2) {
        Expr* lhs = el->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == x->data.symbol.name)
            return el->data.function.args[1];
        return NULL;
    }
    if (head_name_is(el, "List")) {
        for (size_t i = 0; i < el->data.function.arg_count; i++) {
            Expr* rr = nl_rule_rhs_for(el->data.function.args[i], x);
            if (rr) return rr;
        }
    }
    return NULL;
}

/* The integer-constrained parameter of a Solve condition, e.g. C[1] in
 * Element[C[1], Integers] (recursing through And/Or).  Borrowed; NULL if none. */
static const Expr* nl_integer_param_of(const Expr* cond) {
    if (!cond || cond->type != EXPR_FUNCTION) return NULL;
    if (head_name_is(cond, "Element") && cond->data.function.arg_count == 2) {
        Expr* dom = cond->data.function.args[1];
        if (dom->type == EXPR_SYMBOL && strcmp(dom->data.symbol.name, "Integers") == 0)
            return cond->data.function.args[0];
    }
    for (size_t i = 0; i < cond->data.function.arg_count; i++) {
        const Expr* p = nl_integer_param_of(cond->data.function.args[i]);
        if (p) return p;
    }
    return NULL;
}

/* ReplaceAll[body, param -> val], evaluated.  Consumes val. */
static Expr* nl_subst_param(const Expr* body, const Expr* param, Expr* val) {
    return eval_take(mk_fn2("ReplaceAll", expr_copy((Expr*)body),
                            mk_fn2("Rule", expr_copy((Expr*)param), val)));
}

/* A trig pole solves as a periodic family ConditionalExpression[body(C),
 * Element[C, Integers]] that nl_numeric cannot localize.  Enumerate its concrete
 * members strictly inside (a, b) -- needs numeric bounds and a numeric period --
 * and insert each as an interior pole.  Sets undecidable when it cannot. */
static void nl_enum_periodic_pole(const Expr* body, const Expr* param,
                                  Expr* a, Expr* b, PoleSet* s) {
    double av, bv;
    if (nl_numeric(a, &av) != 1 || nl_numeric(b, &bv) != 1) { s->undecidable = true; return; }
    if (bv < av) { double t = av; av = bv; bv = t; }

    Expr* base = nl_subst_param(body, param, expr_new_integer(0));
    Expr* b1   = nl_subst_param(body, param, expr_new_integer(1));
    double basev, b1v;
    bool okb = base && nl_numeric(base, &basev) == 1;
    bool ok1 = b1   && nl_numeric(b1, &b1v) == 1;
    if (base) expr_free(base);
    if (b1) expr_free(b1);
    if (!okb || !ok1) { s->undecidable = true; return; }

    double period = b1v - basev;
    if (fabs(period) < 1e-12) return;
    double k0 = (av - basev) / period, k1 = (bv - basev) / period;
    double klo = floor(fmin(k0, k1)) - 1.0, khi = ceil(fmax(k0, k1)) + 1.0;
    if (khi - klo > 100000.0) { s->undecidable = true; return; }
    for (long k = (long)klo; k <= (long)khi; k++) {
        double pt = basev + period * (double)k;
        if (pt <= av + 1e-9 || pt >= bv - 1e-9) continue;        /* strict interior */
        Expr* member = nl_subst_param(body, param, expr_new_integer((int64_t)k));
        double mv;
        if (member && nl_numeric(member, &mv) == 1) poleset_insert(s, member, mv);
        else if (member) expr_free(member);
    }
}

/* Solve[g == 0, x, Reals] and add every real root strictly interior to (a, b)
 * to `s` (a concrete root directly, a periodic ConditionalExpression family via
 * nl_enum_periodic_pole).  `g` is borrowed. */
static void nl_add_zero_roots(Expr* g, Expr* x, Expr* a, Expr* b, PoleSet* s, Expr* as) {
    if (!g || !contains_symbol(g, x)) return;
    /* A product vanishes exactly where one of its factors does, so solve each
     * x-dependent factor on its own.  This keeps a never-zero factor from
     * poisoning the whole denominator: Together[E^-t/t^2] has the denominator
     * t^2 E^t, which Solve cannot treat as a unit, so the pole set went
     * undecidable and every improper e^{-t}/t^n integral declined. */
    if (head_name_is(g, "Times")) {
        for (size_t i = 0; i < g->data.function.arg_count; i++)
            nl_add_zero_roots(g->data.function.args[i], x, a, b, s, as);
        return;
    }
    /* E^u is never zero for any finite u: it contributes no pole. */
    if (head_name_is(g, "Power") && g->data.function.arg_count == 2 &&
        g->data.function.args[0]->type == EXPR_SYMBOL &&
        g->data.function.args[0]->data.symbol.name == SYM_E)
        return;
    Expr* eq   = mk_fn2("Equal", expr_copy(g), expr_new_integer(0));
    /* Quiet: Solve is an internal pole probe; its Solve::ifun note on a
     * transcendental denominator (Cos[x] == 0) must not surface from Integrate. */
    Expr* sols = eval_take(mk_fn1("Quiet",
                     mk_fn3("Solve", eq, expr_copy(x), mk_sym("Reals"))));
    if (!sols) { s->undecidable = true; return; }
    if (!head_name_is(sols, "List")) {
        /* Solve could not handle it (e.g. x^b == 0 for symbolic b).  A Power
         * denominator base^exp vanishes exactly where its base does, so solve
         * the base instead rather than failing undecidable. */
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2 &&
            contains_symbol(g->data.function.args[0], x)) {
            expr_free(sols);
            nl_add_zero_roots(g->data.function.args[0], x, a, b, s, as);
            return;
        }
        s->undecidable = true; expr_free(sols); return;
    }

    for (size_t i = 0; i < sols->data.function.arg_count; i++) {
        Expr* val = nl_rule_rhs_for(sols->data.function.args[i], x);
        if (!val) { s->undecidable = true; continue; }

        /* Periodic trig root family (Cos[x] == 0 -> Pi/2 + C[1] Pi): enumerate
         * its concrete interior members rather than failing undecidable. */
        if (head_name_is(val, "ConditionalExpression") &&
            val->data.function.arg_count == 2) {
            const Expr* param = nl_integer_param_of(val->data.function.args[1]);
            Expr* body = val->data.function.args[0];
            if (param && !contains_symbol(body, x))
                nl_enum_periodic_pole(body, param, a, b, s);
            else
                s->undecidable = true;
            continue;
        }

        double rv = 0.0;
        int rk = nl_numeric(val, &rv);
        if (rk == 2 || rk == 3) continue;                  /* infinite root: skip */

        /* Classify numerically when possible; fall back to an assumption-aware
         * symbolic comparison for a symbolic root (x = a) or symbolic bound. */
        NLClass cls = (rk == 1) ? nl_classify(rv, a, b) : NL_UNKNOWN;
        if (cls == NL_UNKNOWN) cls = nl_classify_symbolic(val, a, b, as);

        if (cls == NL_INTERIOR) {
            if (rk == 1) poleset_insert(s, expr_copy(val), rv);
            else if (s->symbolic && s->n == 1 && expr_eq(s->pts[0], val))
                ;                                 /* same symbolic root again */
            else if (s->n == 0 && !s->undecidable) {
                /* The sole interior pole, symbolic (x = a under a > 0): placeable. */
                poleset_insert(s, expr_copy(val), NAN);
                s->symbolic = true;
            } else s->undecidable = true;  /* symbolic root among others: cannot order */
        } else if (cls == NL_UNKNOWN) {
            s->undecidable = true;
        }
        /* NL_EXTERIOR: ignore */
    }
    expr_free(sols);
}

/* Poles contributed by trig heads that are secretly ratios: Tan[u] and Sec[u]
 * blow up where Cos[u] == 0, Cot[u] and Csc[u] where Sin[u] == 0.  These are
 * invisible to Denominator[Together[.]] because the evaluator keeps (and
 * re-canonicalizes Sin/Cos back to) the atomic head, so they must be collected
 * directly from the head.  Walks `e` for such heads with an x-dependent
 * argument. */
static void nl_collect_trig_poles(Expr* e, Expr* x, Expr* a, Expr* b, PoleSet* s, Expr* as) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.arg_count == 1 && contains_symbol(e->data.function.args[0], x)) {
        Expr* u = e->data.function.args[0];
        if (head_name_is(e, "Tan") || head_name_is(e, "Sec")) {
            Expr* den = mk_fn1("Cos", expr_copy(u));
            nl_add_zero_roots(den, x, a, b, s, as);
            expr_free(den);
        } else if (head_name_is(e, "Cot") || head_name_is(e, "Csc")) {
            Expr* den = mk_fn1("Sin", expr_copy(u));
            nl_add_zero_roots(den, x, a, b, s, as);
            expr_free(den);
        }
    }
    nl_collect_trig_poles(e->data.function.head, x, a, b, s, as);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nl_collect_trig_poles(e->data.function.args[i], x, a, b, s, as);
}

/* Collect the interior real poles of `expr` (as a function of `x`) into `s`:
 * the real roots of Denominator[Together[expr]] plus the Cos/Sin zeros under any
 * Tan/Cot/Sec/Csc head, all lying strictly in (a, b).  Sets s->undecidable if a
 * pole exists whose classification cannot be decided. */
static void nl_collect_poles(Expr* expr, Expr* x, Expr* a, Expr* b, PoleSet* s, Expr* as) {
    Expr* tog = eval_take(mk_fn1("Together", expr_copy(expr)));
    if (tog) {
        Expr* den = eval_take(mk_fn1("Denominator", tog));   /* tog consumed */
        if (den) { nl_add_zero_roots(den, x, a, b, s, as); expr_free(den); }
    }
    nl_collect_trig_poles(expr, x, a, b, s, as);
}

/* -------------------------------------------------------------------------
 * Boundary evaluation of the antiderivative.
 * ---------------------------------------------------------------------- */

typedef enum { NL_DIR_NONE, NL_DIR_ABOVE, NL_DIR_BELOW } NLDir;

static bool nl_is_infinite(Expr* e) {
    return is_infinity_sym(e) || is_neg_infinity_form(e) ||
           is_complex_infinity_sym(e);
}

/* True iff any subexpression of `e` is a singular sentinel (Infinity /
 * -Infinity / ComplexInfinity / Indeterminate).  A boundary substitution that
 * merely *contains* one -- e.g. ArcTan[ComplexInfinity/Sqrt[3]] from evaluating
 * a Weierstrass antiderivative exactly on a Tan[x/2] pole -- is not a usable
 * closed form and must fall back to a one-sided limit. */
static bool nl_contains_singular(const Expr* e) {
    if (!e) return false;
    if (nl_is_infinite((Expr*)e) || is_indeterminate_sym((Expr*)e)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (nl_contains_singular(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_contains_singular(e->data.function.args[i])) return true;
    return false;
}

/* True iff `e` contains 0^p with a non-numeric exponent p.  Such a power is an
 * undecided boundary value, not a closed form: x^(1+nu) MeijerG[..., x^2/4]
 * substituted at 0 gives 0^(1+nu) MeijerG[..., 0], a 0 * Infinity product
 * whose true limit depends on how the exponents combine. */
static bool nl_contains_zero_power(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is((Expr*)e, "Power") && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* p = e->data.function.args[1];
        if (b->type == EXPR_INTEGER && b->data.integer == 0 &&
            p->type != EXPR_INTEGER && p->type != EXPR_REAL && p->type != EXPR_BIGINT &&
            !head_name_is((Expr*)p, "Rational"))
            return true;
    }
    if (nl_contains_zero_power(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_contains_zero_power(e->data.function.args[i])) return true;
    return false;
}

/* A boundary value is usable directly only if it is a finite closed form:
 * free of x, and free of any singular sentinel / undecided 0^p / unresolved
 * Limit or Integrate. */
static bool nl_is_finite_closed(Expr* v, const Expr* x) {
    if (!v) return false;
    if (nl_contains_singular(v)) return false;
    if (nl_contains_zero_power(v)) return false;
    if (head_name_is(v, "Interval")) return false;
    if (contains_head(v, "Limit") || contains_head(v, "Integrate")) return false;
    if (contains_symbol(v, x)) return false;
    return true;
}

/* Limit[F, x -> target, Direction -> dir, Assumptions -> as]. */
static Expr* nl_limit(Expr* F, Expr* x, Expr* target, NLDir dir, Expr* as) {
    Expr* rule = mk_fn2("Rule", expr_copy(x), expr_copy(target));
    Expr* call;
    if (dir == NL_DIR_NONE && !as) {
        call = mk_fn2("Limit", expr_copy(F), rule);
    } else {
        /* Build the trailing options: Direction (if directed) and Assumptions. */
        Expr** opts = malloc(3 * sizeof(Expr*));
        size_t no = 0;
        opts[no++] = expr_copy(F);
        opts[no++] = rule;
        if (dir != NL_DIR_NONE)
            opts[no++] = mk_fn2("Rule", mk_sym("Direction"),
                             expr_new_string(dir == NL_DIR_ABOVE ? "FromAbove" : "FromBelow"));
        if (as)
            opts[no++] = mk_fn2("Rule", mk_sym("Assumptions"), expr_copy(as));
        call = expr_new_function(mk_sym("Limit"), opts, no);
        free(opts);
    }
    return eval_take(call);
}

/* Evaluate F at `target`, approaching per `dir`.  Uses direct substitution
 * (Refine'd under `as`) at ordinary finite points and the Limit engine at
 * infinite / singular points (or when substitution is itself singular). */
static Expr* nl_eval_at(Expr* F, Expr* x, Expr* target, NLDir dir, Expr* as) {
    arith_warnings_mute_push();
    Expr* out;
    if (nl_is_infinite(target)) {
        out = nl_limit(F, x, target, NL_DIR_NONE, as);
    } else {
        /* Quiet: a direct boundary substitution can hit a Power::infy on a form
         * like 0^(1-b) that the Refine/Limit below then resolves; it is an
         * internal probe and must not surface from Integrate. */
        Expr* sub = eval_take(mk_fn1("Quiet",
                                 mk_fn2("ReplaceAll", expr_copy(F),
                                        mk_fn2("Rule", expr_copy(x),
                                               expr_copy(target)))));
        /* Under Assumptions, Refine the substituted boundary so a symbolic
         * endpoint closes cleanly: ArcSin[a/Sqrt[a^2]] -> Pi/2 for a>0, and a
         * transient 0^(1-b) -> 0 for 0<b<1 (left unresolved it later leaks a
         * Power::infy).  Refine of an already-clean value is a no-op.  Quiet: an
         * internal boundary-cleanup probe. */
        if (sub && as && !nl_is_infinite(sub) && !contains_symbol(sub, x)) {
            Expr* ref = eval_take(mk_fn1("Quiet", mk_fn2("Refine", sub, expr_copy(as))));  /* consumes sub */
            if (ref) sub = ref;
        }
        if (sub && nl_is_infinite(sub)) {
            /* Direct substitution lands on a bare directed infinity: F itself
             * blows up at this interior finite point (e.g. F = -ArcTanh[x] at
             * x = 1, where ArcTanh[1] = Infinity).  That is a genuine boundary
             * divergence -- keep the sentinel so it propagates to the idiv
             * path -- and it is a strictly better answer than the one-sided
             * Limit, which the engine may leave unevaluated. */
            out = sub;
        } else if (nl_is_finite_closed(sub, x)) {
            out = sub;
        } else {
            if (sub) expr_free(sub);
            out = nl_limit(F, x, target, dir, as);
        }
    }
    arith_warnings_mute_pop();
    return out;
}

/* -------------------------------------------------------------------------
 * Divergence reporting + continuity cross-check.
 * ---------------------------------------------------------------------- */

/* Emit the Mathematica-style Integrate::idiv warning to stderr, naming the
 * integrand and interval.  Shared (declared in the header) so other definite
 * mechanisms can report divergence with the same diagnostic. */
void integrate_emit_idiv(Expr* f, Expr* a, Expr* b) {
    char* fs = expr_to_string(f);
    char* as = expr_to_string(a);
    char* bs = expr_to_string(b);
    mth_message("Integrate", "idiv", "Integral of %s does not converge on {%s, %s}.", fs ? fs : "?", as ? as : "?", bs ? bs : "?");
    free(fs); free(as); free(bs);
}

/* Numerically cross-check a symbolic definite value `V` against NIntegrate.
 * Returns true iff both are finite reals agreeing to a loose relative
 * tolerance -- an uncorrected branch jump is a large (O(Pi/period)) discrepancy
 * far above the tolerance, while quadrature noise sits far below it, so this
 * catches every value-affecting jump without being fooled by NIntegrate's own
 * inaccuracy.  Returns false (reject) if NIntegrate cannot produce a value. */
static bool nl_crosscheck_ok(Expr* f, Expr* x, Expr* a, Expr* b, Expr* V) {
    double vd;
    if (nl_numeric(V, &vd) != 1) {
        /* A real value telescoped from a complex-written antiderivative carries
         * an imaginary rounding residue (~1e-16); accept its real part when the
         * imaginary part is negligible, reject a genuinely complex value. */
        Expr* nv = eval_take(mk_fn1("N", expr_copy(V)));
        Expr *re = NULL, *im = NULL;
        double rv, iv;
        bool ok = nv && is_complex(nv, &re, &im) && nl_real_double(re, &rv) &&
                  nl_real_double(im, &iv) && fabs(iv) <= 1e-9 * (1.0 + fabs(rv));
        if (nv) expr_free(nv);
        if (!ok) return false;                    /* V not a finite (real) value */
        vd = rv;
    }

    arith_warnings_mute_push();
    Expr* spec = mk_fn3("List", expr_copy(x), expr_copy(a), expr_copy(b));
    Expr* ni   = eval_take(mk_fn2("NIntegrate", expr_copy(f), spec)); /* spec consumed */
    arith_warnings_mute_pop();

    double nd;
    int k = ni ? nl_numeric(ni, &nd) : 0;
    if (ni) expr_free(ni);
    if (k != 1) return false;                     /* cannot validate -> reject */
    return fabs(vd - nd) <= 1e-3 * (1.0 + fabs(vd));
}

/* Evaluate `expr` with x := xv (a machine real), classifying like nl_numeric:
 * 0 = not a real number, 1 = finite real (in *out), 2/3 = +/-Infinity. */
static int nl_eval_num_at(Expr* expr, Expr* x, double xv, double* out) {
    Expr* sub = eval_take(mk_fn2("ReplaceAll", expr_copy(expr),
                                 mk_fn2("Rule", expr_copy(x), expr_new_real(xv))));
    if (!sub) return 0;
    int k = nl_numeric(sub, out);
    expr_free(sub);
    return k;
}

/* Numeric resolution of an undecidable pole set: scan Denominator[Together[f]]
 * across (a, b) for a real root.  Solve sometimes returns spurious complex
 * conditional roots for a nowhere-vanishing trig denominator (e.g. 2+Cos[x]),
 * flagging the pole set undecidable when there is in fact no real pole; a direct
 * sign-change / near-tangency scan settles it.  Returns
 *    1 = a real root lies strictly inside (a, b) (genuine interior pole),
 *    0 = no real root found on the grid,
 *   -1 = cannot scan (non-numeric/infinite bounds, or den not real on the grid). */
static int nl_denominator_interior_root_scan(Expr* f, Expr* x, Expr* a, Expr* b) {
    Expr* tog = eval_take(mk_fn1("Together", expr_copy(f)));
    if (!tog) return -1;
    Expr* den = eval_take(mk_fn1("Denominator", tog));   /* tog consumed */
    if (!den) return -1;
    if (!contains_symbol(den, x)) { expr_free(den); return 0; }

    double av, bv;
    if (nl_numeric(a, &av) != 1 || nl_numeric(b, &bv) != 1) { expr_free(den); return -1; }
    if (bv < av) { double t = av; av = bv; bv = t; }

    const int N = 512;
    double span = bv - av;
    double prev = 0.0;
    bool have_prev = false;
    int result = 0;
    for (int i = 0; i <= N && result == 0; i++) {
        double xv = av + span * ((double)i / (double)N);
        double dv;
        int k = nl_eval_num_at(den, x, xv, &dv);
        if (k != 1) { result = -1; break; }        /* den not real -> undecidable */
        bool interior = (i != 0 && i != N);
        if (interior && fabs(dv) < 1e-9) { result = 1; break; }   /* tangency */
        if (have_prev && prev * dv < 0.0)  { result = 1; break; } /* sign change */
        prev = dv; have_prev = true;
    }
    expr_free(den);
    return result;
}

/* True if exponent `e` denotes a non-integer power (a radical): a Rational with
 * denominator != 1, or a non-integer real.  A radical's real branch can kink or
 * jump where its radicand vanishes, which is the hazard the scan below looks
 * for. */
static bool nl_is_fractional_exponent(const Expr* e) {
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) return d != 1 && d != -1;
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        return fabs(v - nearbyint(v)) > 1e-12;
    }
    return false;
}

/* Grid-scan `g` over (av, bv) for an interior real zero (sign change across a
 * subinterval, or a near-tangency at an interior node).  Samples where g is not
 * a finite real are skipped (they break the running sign but are not evidence of
 * a zero). */
static bool nl_expr_interior_zero(Expr* g, const Expr* x, double av, double bv) {
    if (bv < av) { double t = av; av = bv; bv = t; }
    const int N = 256;
    double span = bv - av;
    double prev = 0.0; bool have_prev = false;
    for (int i = 0; i <= N; i++) {
        double xv = av + span * ((double)i / (double)N);
        double gv;
        if (nl_eval_num_at(g, (Expr*)x, xv, &gv) != 1) { have_prev = false; continue; }
        bool interior = (i != 0 && i != N);
        if (interior && fabs(gv) < 1e-9) return true;          /* tangency */
        if (have_prev && prev * gv < 0.0) return true;          /* sign change */
        prev = gv; have_prev = true;
    }
    return false;
}

/* Recursive walk: does `e` contain a fractional-power node Power[g, p] (p a
 * non-integer) whose radicand g vanishes strictly inside (av, bv)? */
static bool nl_radical_hazard_walk(Expr* e, const Expr* x, double av, double bv) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* expo = e->data.function.args[1];
        if (nl_is_fractional_exponent(expo) && contains_symbol(base, x) &&
            nl_expr_interior_zero(base, x, av, bv))
            return true;
    }
    if (nl_radical_hazard_walk(e->data.function.head, x, av, bv)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_radical_hazard_walk(e->data.function.args[i], x, av, bv)) return true;
    return false;
}

/* True when the antiderivative F carries a radical whose radicand vanishes
 * inside the integration interval: the real branch of such a radical (e.g.
 * Sqrt[1 - Cos[2x]] = Sqrt[2] Abs[Sin x], which touches 0 at x = Pi) can kink or
 * flip sign there, so a single-branch F(b) - F(a) is unsafe and must be
 * numerically cross-checked.  Needs finite numeric bounds to scan; returns false
 * for symbolic bounds (which the caller's own undecidable-pole path covers). */
static bool nl_find_radical_hazard(Expr* F, const Expr* x, Expr* a, Expr* b) {
    double av, bv;
    if (nl_numeric(a, &av) != 1 || nl_numeric(b, &bv) != 1) return false;
    return nl_radical_hazard_walk(F, x, av, bv);
}

/* -------------------------------------------------------------------------
 * Core Newton-Leibniz driver.
 * ---------------------------------------------------------------------- */

/* Cauchy principal value over [a, b] given the antiderivative F and the interior
 * pole set: valid only when every interior pole is odd-order, detected by a sign
 * change of f across it (odd -> one-sided divergences cancel; even -> no PV).
 * PV = Re[F(b) - F(a)] with F's real branch: the branch-cut crossing at each odd
 * pole contributes a purely imaginary I*Pi*residue that Re removes, leaving the
 * finite real principal value.  Returns NULL when the gate fails or the value is
 * not a finite real. */
static Expr* nl_principal_value(Expr* F, Expr* f, Expr* x, Expr* a, Expr* b,
                                PoleSet* poles, Expr* as);

/* A known analytic constant that must NOT be instantiated when probing realness. */
static bool nl_is_const_name(const char* n) {
    static const char* const C[] = { "Pi", "E", "I", "EulerGamma", "Infinity",
        "Degree", "GoldenRatio", "Catalan", "ComplexInfinity", "Indeterminate",
        "Glaisher", "Khinchin" };
    for (size_t i = 0; i < sizeof C / sizeof *C; i++)
        if (strcmp(n, C[i]) == 0) return true;
    return false;
}

/* Collect `s -> 2` rules for each distinct non-constant LEAF variable symbol of
 * `e`.  A symbol in head position (Plus, Log, CosIntegral, ...) is a function
 * name, not a variable, so the head is recursed into only when it is itself a
 * function -- a bare-symbol head is skipped. */
static void nl_collect_vars(const Expr* e, Expr*** rules, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (nl_is_const_name(e->data.symbol.name)) return;
        for (size_t i = 0; i < *n; i++)
            if (strcmp((*rules)[i]->data.function.args[0]->data.symbol.name,
                       e->data.symbol.name) == 0) return;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                          *rules = realloc(*rules, *cap * sizeof(Expr*)); }
        (*rules)[(*n)++] = mk_fn2("Rule", expr_copy((Expr*)e), expr_new_integer(2));
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    if (e->data.function.head && e->data.function.head->type == EXPR_FUNCTION)
        nl_collect_vars(e->data.function.head, rules, n, cap);   /* a computed head */
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        nl_collect_vars(e->data.function.args[i], rules, n, cap);
}

/* A principal value computed as Re[F(b)-F(a)] is trustworthy only if it is in
 * fact real: instantiate every non-constant free variable at a sample value and
 * confirm N[.] is a finite real.  Catches the case where Re[] failed to strip a
 * branch-cut's imaginary residue because a special function (CosIntegral[-a]) is
 * treated as real symbolically but evaluates complex -- better to decline
 * (idiv) than return a spurious complex principal value. */
static bool nl_realness_ok(Expr* pvval, Expr* as) {
    Expr** rules = NULL; size_t n = 0, cap = 0;
    nl_collect_vars(pvval, &rules, &n, &cap);
    /* The sample point must satisfy the assumptions: under 0 < a < b the
     * all-equal sample a = b = 2 sits on the boundary a == b, where the genuine
     * PV Log[b - a] - Log[a] is -Infinity, and the probe would reject a correct
     * real value.  Try the all-2 sample first (unchanged when it is admissible),
     * then distinct increasing and decreasing samples, keeping the first the
     * assumptions accept; with none admissible, fall back to all-2. */
    static const int64_t INC[] = { 2, 3, 5, 7, 11, 13, 17, 19 };
    const size_t NI = sizeof INC / sizeof *INC;
    int pick = 0;
    if (as && n > 0 && n <= NI) {
        for (int pat = 0; pat < 3; pat++) {
            Expr** rr = malloc(n * sizeof(Expr*));
            for (size_t k = 0; k < n; k++) {
                int64_t v = pat == 0 ? 2 : pat == 1 ? INC[k] : INC[n - 1 - k];
                rr[k] = mk_fn2("Rule", expr_copy(rules[k]->data.function.args[0]),
                               expr_new_integer(v));
            }
            Expr* rl = expr_new_function(mk_sym("List"), rr, n);
            free(rr);
            Expr* ok = eval_take(mk_fn2("ReplaceAll", expr_copy(as), rl));
            bool admissible = ok && ok->type == EXPR_SYMBOL && ok->data.symbol.name == SYM_True;
            if (ok) expr_free(ok);
            if (admissible) { pick = pat; break; }
        }
    }
    for (size_t k = 0; k < n && pick != 0; k++) {
        int64_t v = pick == 1 ? INC[k] : INC[n - 1 - k];
        expr_free(rules[k]->data.function.args[1]);
        rules[k]->data.function.args[1] = expr_new_integer(v);
    }
    Expr* subbed;
    /* An internal probe: instantiating at a sample can hit a transient
     * Power::infy / Infinity::indet (e.g. 0^(1-b) -> 0^-1 at b=2).  Quiet alone
     * would still let an enclosing Check[] see it, so the arithmetic probe mute
     * brackets it too. */
    arith_warnings_mute_push();
    if (n == 0) {
        subbed = eval_take(mk_fn1("Quiet", mk_fn1("N", expr_copy(pvval))));
    } else {
        Expr* rl = expr_new_function(mk_sym("List"), rules, n);  /* copies the array */
        subbed = eval_take(mk_fn1("Quiet", mk_fn1("N", mk_fn2("ReplaceAll", expr_copy(pvval), rl))));
    }
    arith_warnings_mute_pop();
    free(rules);
    double v;
    bool ok = subbed && nl_numeric(subbed, &v) == 1;   /* finite real only */
    if (subbed) expr_free(subbed);
    return ok;
}

/* True iff `e`, with its free variables instantiated at a sample value, N[.]'s
 * to a finite complex number with a nonzero imaginary part.  This catches a
 * value that is provably complex even when the symbolic Im[.] cannot see it
 * (CosIntegral[-a] is treated as real symbolically but evaluates complex).  A
 * non-numericizing (genuinely symbolic) real result returns false -- it is NOT
 * flagged, so a legitimate symbolic answer is never rejected.  Exposed as
 * integrate_result_is_complex() so the dispatcher's Refine cleanup can avoid
 * INTRODUCING complexity (stripping a Re that was protecting realness). */
bool nl_numeric_is_complex(Expr* e) {
    Expr** rules = NULL; size_t n = 0, cap = 0;
    nl_collect_vars(e, &rules, &n, &cap);
    Expr* subbed;
    /* Quiet: instantiating at a sample can hit a transient Power::infy; probe. */
    if (n == 0) {
        subbed = eval_take(mk_fn1("Quiet", mk_fn1("N", expr_copy(e))));
    } else {
        Expr* rl = expr_new_function(mk_sym("List"), rules, n);
        subbed = eval_take(mk_fn1("Quiet", mk_fn1("N", mk_fn2("ReplaceAll", expr_copy(e), rl))));
    }
    free(rules);
    bool is_cplx = false;
    if (subbed) {
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex(subbed, &re, &im)) {
            double iv;
            is_cplx = nl_real_double(im, &iv) && fabs(iv) > 1e-9;
        }
        expr_free(subbed);
    }
    return is_cplx;
}

static Expr* nl_principal_value(Expr* F, Expr* f, Expr* x, Expr* a, Expr* b,
                                PoleSet* poles, Expr* as) {
    /* Numeric bounds let the odd-order gate size its probe delta from the
     * interval; symbolic bounds (with Assumptions) are allowed too -- the gate
     * is a LOCAL sign-flip test at each (concrete-valued) interior pole, so it
     * only needs a small probe delta, and the endpoint values come from the
     * assumption-aware boundary evaluation below. */
    double av, bv;
    int ak = nl_numeric(a, &av), bk = nl_numeric(b, &bv);
    bool num_bounds = (ak == 1 && bk == 1);
    if (num_bounds && bv < av) { double t = av; av = bv; bv = t; }
    if (!num_bounds && !as) return NULL;   /* symbolic bounds need assumptions */

    /* Odd-order gate: sign of f must flip across every interior pole.  A sole
     * symbolic pole (x = a) has no numeric probe point; it is admitted when it
     * is provably SIMPLE: Limit[(x - p) f, x -> p] is a finite nonzero residue
     * under the assumptions (a simple pole is odd-order, its one-sided
     * divergences cancel).  Higher-order symbolic poles decline. */
    if (poles->symbolic) {
        if (!as || poles->n != 1) return NULL;
        Expr* p = poles->pts[0];
        Expr* g = mk_fn2("Times", mk_fn2("Plus", expr_copy(x),
                             mk_fn2("Times", expr_new_integer(-1), expr_copy(p))),
                         expr_copy(f));
        arith_warnings_mute_push();
        Expr* res = nl_limit(g, x, p, NL_DIR_NONE, as);
        Expr* nz = (res && nl_is_finite_closed(res, x))
            ? eval_take(mk_fn2("Refine", mk_fn2("Unequal", expr_copy(res),
                                                 expr_new_integer(0)), expr_copy(as)))
            : NULL;
        arith_warnings_mute_pop();
        expr_free(g);
        bool simple = nz && nz->type == EXPR_SYMBOL && nz->data.symbol.name == SYM_True;
        if (res) expr_free(res);
        if (nz) expr_free(nz);
        if (!simple) return NULL;
    }
    for (size_t i = 0; i < poles->n && !poles->symbolic; i++) {
        double pv = poles->vals[i];
        double near_other = HUGE_VAL;
        for (size_t j = 0; j < poles->n; j++)
            if (j != i) near_other = fmin(near_other, fabs(pv - poles->vals[j]));
        double near = num_bounds ? fmin(fmin(pv - av, bv - pv), near_other) : near_other;
        if (num_bounds && !(near > 0.0)) return NULL;
        /* delta: a small fraction of the nearest feature, capped at 1e-3 so a
         * symbolic-bound single pole (near == HUGE_VAL) still gets a tight probe. */
        double delta = (near > 0.0 && near < HUGE_VAL) ? fmin(near * 1e-3, 1e-3) : 1e-3;
        double fl, fr;
        if (nl_eval_num_at(f, x, pv - delta, &fl) != 1) return NULL;
        if (nl_eval_num_at(f, x, pv + delta, &fr) != 1) return NULL;
        if (fl == 0.0 || fr == 0.0 || fl * fr > 0.0) return NULL;   /* even order / flat */
    }

    Expr* Fa = nl_eval_at(F, x, a, NL_DIR_ABOVE, as);
    Expr* Fb = nl_eval_at(F, x, b, NL_DIR_BELOW, as);
    if (!Fa || !Fb) { if (Fa) expr_free(Fa); if (Fb) expr_free(Fb); return NULL; }

    Expr* diff = mk_fn2("Plus", Fb, mk_fn2("Times", expr_new_integer(-1), Fa));
    /* Re removes the pole branch-cut's I*Pi residue; Simplify under the
     * assumptions resolves Log[-a] = Log[a] + I Pi for a symbolic bound. */
    Expr* pvval = eval_take(as ? mk_fn1("Quiet", mk_fn2("Simplify",
                                        mk_fn1("Re", expr_copy(diff)), expr_copy(as)))
                               : mk_fn1("Simplify", mk_fn1("Re", expr_copy(diff))));
    if (pvval && nl_is_finite_closed(pvval, x) && nl_realness_ok(pvval, as)) {
        expr_free(diff);
        return pvval;   /* genuinely real (0, Log[b/a], 2 SinhIntegral[a], ...) */
    }
    if (pvval) expr_free(pvval);
    /* Simplify under the assumptions can wrongly STRIP the Re from a term it
     * treats as real but that evaluates complex (CosIntegral[-a]), leaving a
     * spurious complex value.  Fall back to the literal Re[F(b)-F(a)], which
     * forces the real part; accept it only if it is a real finite closed form. */
    Expr* pvre = eval_take(mk_fn1("Re", diff));   /* consumes diff */
    if (pvre && nl_is_finite_closed(pvre, x) && nl_realness_ok(pvre, as)) return pvre;
    if (pvre) expr_free(pvre);
    return NULL;
}

Expr* integrate_newton_leibniz_try(Expr* f, Expr* x, Expr* a, Expr* b,
                                   const char* method) {
    return integrate_newton_leibniz_try_pv(f, x, a, b, method, false, NULL);
}

/* True if the integrand f has a NON-INTEGRABLE singularity at the finite point pt:
 * Limit[(x - pt) f, x -> pt] is nonzero or infinite (a pole of order >= 1).  A
 * removable or integrable singularity (continuous f, Log, x^(-p) with p<1) gives
 * the limit 0.  The test is on the integrand itself (not on the antiderivative's
 * boundary limit, which a nested evaluation context can mis-return as finite), so
 * it is context-independent.  An undecided limit (x survives / Indeterminate) is
 * treated as integrable -- never force a divergence that isn't proven. */
static bool nl_nonintegrable_at(Expr* f, Expr* x, Expr* pt, Expr* as) {
    /* f finite at pt by direct substitution: not a pole.  (The limit test below
     * on a radical endpoint, e.g. b = 2 - Sqrt[3], can return an unsimplified
     * zero whose N value is a rounding residue, which read as a simple pole and
     * declared Integrate[Log[x]/(1 + x^2), {x, 0, 2 - Sqrt[3]}] divergent.) */
    {
        arith_warnings_mute_push();
        Expr* at = eval_take(mk_fn1("Quiet", mk_fn1("N",
                       mk_fn2("ReplaceAll", expr_copy(f),
                              mk_fn2("Rule", expr_copy(x), expr_copy(pt))))));
        arith_warnings_mute_pop();
        double v;
        bool finite = at && nl_numeric(at, &v) == 1 && isfinite(v);
        if (at) expr_free(at);
        if (finite) return false;
    }
    Expr* shift = mk_fn2("Plus", expr_copy(x),
                      mk_fn2("Times", expr_new_integer(-1), expr_copy(pt)));
    Expr* g = mk_fn2("Times", shift, expr_copy(f));
    Expr* lim = nl_limit(g, x, pt, NL_DIR_NONE, as);   /* borrows g */
    expr_free(g);
    if (!lim) return false;
    bool nonint = false;
    if (nl_is_infinite(lim)) {
        nonint = true;                                  /* order > 1 pole */
    } else if (nl_is_finite_closed(lim, x)) {
        double v;
        bool is_zero = (lim->type == EXPR_INTEGER && lim->data.integer == 0) ||
                       (lim->type == EXPR_REAL && lim->data.real == 0.0) ||
                       (nl_numeric(lim, &v) == 1 && fabs(v) < 1e-12);
        nonint = !is_zero;                              /* nonzero finite -> simple pole */
    }
    expr_free(lim);
    return nonint;
}

/* True iff e contains the symbol named `name`. */
static bool contains_symbol_named(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol_named(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol_named(e->data.function.args[i], name)) return true;
    return false;
}

/* Node count of e, stopping once it exceeds cap (returns a value > cap). */
static size_t nl_node_count_capped(const Expr* e, size_t cap) {
    if (!e) return 0;
    size_t n = 1;
    if (e->type != EXPR_FUNCTION) return n;
    n += nl_node_count_capped(e->data.function.head, cap);
    for (size_t i = 0; i < e->data.function.arg_count && n <= cap; i++)
        n += nl_node_count_capped(e->data.function.args[i], cap - (n < cap ? n : cap));
    return n;
}

/* True iff F holds a PolyLog or Log whose argument contains an exponential of
 * x (E^(c x) or E^(I c x)).  The polylogarithm-ladder antiderivatives
 * (src/internal/PolyLogLadder.m) are written in T = E^(I b x): as x runs
 * along the real axis such an argument circles the origin and can cross the
 * branch cut of Log/PolyLog, where F jumps although the integrand is smooth --
 * so the telescoped value must be confirmed by the continuity cross-check. */
static bool nl_arg_has_exp_of_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(e->data.function.args[0]->data.symbol.name, "E") == 0 &&
        contains_symbol(e->data.function.args[1], (Expr*)x)) return true;
    if (head_name_is(e, "Exp") && contains_symbol(e, (Expr*)x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nl_arg_has_exp_of_x(e->data.function.args[i], x)) return true;
    return false;
}

static bool nl_has_exp_polylog(const Expr* F, const Expr* x) {
    if (!F || F->type != EXPR_FUNCTION) return false;
    if ((head_name_is(F, "PolyLog") || head_name_is(F, "Log")) &&
        nl_arg_has_exp_of_x(F, x)) return true;
    for (size_t i = 0; i < F->data.function.arg_count; i++)
        if (nl_has_exp_polylog(F->data.function.args[i], x)) return true;
    return false;
}

Expr* integrate_newton_leibniz_try_pv(Expr* f, Expr* x, Expr* a, Expr* b,
                                      const char* method, bool principal_value,
                                      Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;

    /* 1. Antiderivative via the indefinite cascade. */
    Expr* F;
    if (method) {
        F = eval_take(mk_fn3("Integrate", expr_copy(f), expr_copy(x),
                             mk_fn2("Rule", mk_sym("Method"),
                                    expr_new_string(method))));
    } else {
        F = eval_take(mk_fn2("Integrate", expr_copy(f), expr_copy(x)));
    }
    if (!F) return NULL;
    if (contains_head(F, "Integrate")) { expr_free(F); return NULL; }

    /* 2. Interior poles of the INTEGRAND f only.  A continuous f has an
     * antiderivative that is C^1 (hence pole-free) wherever f is, so the only
     * breakpoints that matter are f's own real poles.  Detecting on F as well
     * would manufacture spurious singularities from its representation (e.g.
     * Cos[x/2] inside a Weierstrass ArcTan[Tan[x/2]]), which are not
     * discontinuities of F -- that is handled by the continuity guard below. */
    PoleSet poles; poleset_init(&poles);
    nl_collect_poles(f, x, a, b, &poles, assumptions);
    bool force_validate = false;
    if (poles.undecidable) {
        /* Solve could not place the poles symbolically.  Try to settle it
         * numerically: if the denominator has no real root in (a, b) the
         * "poles" were spurious (complex conditional roots of a nowhere-zero
         * trig denominator); proceed but force the numeric cross-check since
         * detection was incomplete.  Any real interior root, or an unscannable
         * denominator, is left unevaluated. */
        if (nl_denominator_interior_root_scan(f, x, a, b) != 0) {
            poleset_free(&poles); expr_free(F); return NULL;
        }
        poles.undecidable = false;
        force_validate = true;
    }

    /* 2b. Principal value: with interior poles present, compute the Cauchy PV
     * (gated on odd-order poles) instead of the telescoping sum, which would
     * diverge.  With no interior pole the ordinary FTC path below already gives
     * the correct value, so PV falls through to it. */
    if (principal_value && poles.n > 0) {
        Expr* pv = nl_principal_value(F, f, x, a, b, &poles, assumptions);
        poleset_free(&poles);
        expr_free(F);
        if (!pv) integrate_emit_idiv(f, a, b);   /* even-order pole: no PV */
        return pv;
    }

    /* 2c. Ordinary (non-principal-value) integral with interior poles: a genuine
     * non-integrable interior pole makes it DIVERGE.  The telescoping sum below
     * would silently return the Cauchy principal value -- an odd-order pole's two
     * infinite one-sided limits of F cancel (e.g. Integrate[1/x, {x,-1,1}] -> 0),
     * which is only correct under PrincipalValue -> True.  Detect divergence: if F
     * has a non-finite one-sided limit at any interior pole, emit idiv and decline.
     * A removable / integrable singularity keeps F finite there and falls through
     * to the FTC sum unchanged. */
    if (!principal_value) {
        bool divergent = false;
        /* Interior poles: a genuine non-integrable one makes the ordinary integral
         * diverge; the telescoping sum below would otherwise silently return the
         * Cauchy principal value (an odd pole's infinities cancel -> 0). */
        for (size_t i = 0; i < poles.n && !divergent; i++)
            if (nl_nonintegrable_at(f, x, poles.pts[i], assumptions)) divergent = true;
        /* Finite endpoints: a non-integrable boundary pole (1/x on [0,1]) also
         * diverges.  The telescoping's boundary limit can be mis-returned as finite
         * in a nested context (e.g. the symmetry reduction of 1/x on [-1,1] via the
         * half Integrate[1/x,{x,0,1}]), so test the integrand directly. */
        double ev;
        if (!divergent && nl_numeric(a, &ev) == 1 &&
            nl_nonintegrable_at(f, x, a, assumptions)) divergent = true;
        if (!divergent && nl_numeric(b, &ev) == 1 &&
            nl_nonintegrable_at(f, x, b, assumptions)) divergent = true;
        if (divergent) {
            poleset_free(&poles);
            expr_free(F);
            integrate_emit_idiv(f, a, b);
            return NULL;
        }
    }

    /* 3. Breakpoints [a, poles..., b] and the telescoping sum. */
    size_t nbp = poles.n + 2;
    Expr** bp = malloc(nbp * sizeof(*bp));
    bp[0] = a;                                  /* borrowed */
    for (size_t i = 0; i < poles.n; i++) bp[i + 1] = poles.pts[i]; /* borrowed */
    bp[nbp - 1] = b;                            /* borrowed */

    Expr* total = expr_new_integer(0);
    bool failed = false;
    for (size_t i = 0; i + 1 < nbp && !failed; i++) {
        Expr* lo = nl_eval_at(F, x, bp[i],     NL_DIR_ABOVE, assumptions);
        Expr* hi = nl_eval_at(F, x, bp[i + 1], NL_DIR_BELOW, assumptions);
        if (!lo || !hi) {
            if (lo) expr_free(lo);
            if (hi) expr_free(hi);
            failed = true;
            break;
        }
        /* Quiet: forming F(hi) - F(lo) can transiently evaluate a boundary term
         * like 0^(1-b) that raises Power::infy before the telescoping sum /
         * final Simplify resolves it; it is internal and must not surface. */
        arith_warnings_mute_push();
        Expr* contrib = eval_take(mk_fn1("Quiet", mk_fn2("Plus", hi,
                                         mk_fn2("Times", expr_new_integer(-1), lo))));
        total = eval_take(mk_fn1("Quiet", mk_fn2("Plus", total, contrib)));
        arith_warnings_mute_pop();
    }

    free(bp);
    poleset_free(&poles);

    /* Capture F's continuity-risk fingerprint before releasing it: a step head,
     * an inverse-trig node over a pole-bearing trig head, or a radical whose
     * radicand vanishes inside (a,b) all mean F could carry an uncorrected
     * branch jump and must be cross-checked. */
    bool F_triggers = nl_has_step_head(F, x) || nl_has_arctrig_over_pole(F, x)
                      || nl_has_arctan_of_quotient(F, x)
                      || nl_find_radical_hazard(F, x, a, b)
                      || nl_has_exp_polylog(F, x);
    expr_free(F);

    if (failed || !total) { if (total) expr_free(total); return NULL; }

    /* 4a. Genuine divergence: an interior pole (or a boundary blow-up) leaves a
     * divergence sentinel.  Warn as Mathematica does and leave the integral
     * unevaluated. */
    if (nl_is_infinite(total) || is_indeterminate_sym(total)) {
        integrate_emit_idiv(f, a, b);
        expr_free(total);
        return NULL;
    }

    /* 4b. An unresolved Limit/Integrate, a stray x, an Interval, or a residual
     * singular sentinel buried inside the result means the engine could not
     * decide -- leave the integral unevaluated. */
    if (contains_head(total, "Limit") || contains_head(total, "Integrate") ||
        head_name_is(total, "Interval") || contains_symbol(total, x) ||
        nl_contains_singular(total) || nl_contains_zero_power(total)) {
        expr_free(total);
        return NULL;
    }

    /* 4b'. A real-looking integrand over a real interval cannot have a complex
     * definite value: a nonzero imaginary part means the antiderivative crossed
     * a branch cut at an interior singularity the pole scan missed (e.g. Tan
     * with its pole hidden inside the atomic head gives -Log[Cos] -> -I Pi, or
     * an odd integrand telescoping to CosIntegral[a]-CosIntegral[-a] = -I Pi),
     * so the integral actually diverges (or needs a principal value).  Tested
     * NUMERICALLY (symbolic Im[.] cannot see the CosIntegral branch); emit idiv
     * and decline rather than return the complex artifact. */
    if (!contains_head(f, "Complex") && nl_numeric_is_complex(total)) {
        integrate_emit_idiv(f, a, b);
        expr_free(total);
        return NULL;
    }

    /* 4b''. Real form.  A real integrand's value, telescoped from a complex-
     * written antiderivative (a polylogarithm ladder in E^(I x), Log[1 + I x]),
     * arrives as a combination like I (Log[1 - I] - Log[1 + I]) Pi + PolyLog
     * values whose imaginary parts cancel.  ComplexExpand + Simplify expose
     * the real closed form (x Csc[x] on [0, Pi/2] -> 2 Catalan); it is adopted
     * only if it is free of I and numerically equal to the original.  Gated to
     * small PolyLog-bearing values: Simplify[ComplexExpand[..]] on a large
     * complex-exponential value (DSolve's variation of parameters over
     * Sin[t] + DiracDelta[t - 3 Pi]) ran past the 8 s corpus budget. */
    if (!contains_head(f, "Complex") && contains_head(total, "Complex") &&
        (contains_head(total, "PolyLog") || contains_symbol_named(total, "Catalan") ||
         contains_head(total, "Log")) &&
        nl_node_count_capped(total, 400) <= 400) {
        Expr* ce = eval_take(mk_fn1("Quiet", mk_fn1("Simplify", mk_fn1("ComplexExpand", expr_copy(total)))));
        if (ce && !contains_head(ce, "Complex") && !contains_symbol(ce, x)) {
            Expr* dif = eval_take(mk_fn1("N", mk_fn2("Plus", expr_copy(ce),
                                     mk_fn2("Times", expr_new_integer(-1), expr_copy(total)))));
            double dre;
            bool same = dif && ((nl_numeric(dif, &dre) == 1 && fabs(dre) < 1e-9) ||
                                (dif->type == EXPR_INTEGER && dif->data.integer == 0));
            if (!same && dif && head_name_is(dif, "Complex") && dif->data.function.arg_count == 2) {
                double r0, i0;
                same = nl_numeric(dif->data.function.args[0], &r0) == 1 &&
                       nl_numeric(dif->data.function.args[1], &i0) == 1 &&
                       fabs(r0) < 1e-9 && fabs(i0) < 1e-9;
            }
            if (dif) expr_free(dif);
            if (same) { expr_free(total); total = ce; ce = NULL; }
        }
        if (ce) expr_free(ce);
    }

    /* 4c. Continuity guard.  A finite value is trusted only if a numeric
     * NIntegrate cross-check agrees when either the antiderivative is jump-risk
     * (branch/step form) or pole detection was incomplete (bounds must be
     * numeric or +/-Infinity); otherwise bail rather than risk a wrong number. */
    if (F_triggers || force_validate) {
        double av, bv;
        bool bounds_num = nl_numeric(a, &av) && nl_numeric(b, &bv);
        if (!bounds_num || !nl_crosscheck_ok(f, x, a, b, total)) {
            expr_free(total);
            return NULL;
        }
    }
    return total;
}

/* -------------------------------------------------------------------------
 * Spec parsing + builtins.
 * ---------------------------------------------------------------------- */

/* Extract (x, a, b) from a `{x, a, b}` List spec.  Returns borrowed
 * pointers; false if the spec is malformed. */
static bool nl_parse_spec(Expr* spec, Expr** x, Expr** a, Expr** b) {
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3)
        return false;
    Expr* v = spec->data.function.args[0];
    if (v->type != EXPR_SYMBOL) return false;
    *x = v;
    *a = spec->data.function.args[1];
    *b = spec->data.function.args[2];
    return true;
}

Expr* builtin_integrate_newton_leibniz(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x; Expr* a; Expr* b;
    if (!nl_parse_spec(res->data.function.args[1], &x, &a, &b)) return NULL;
    return integrate_newton_leibniz_try(f, x, a, b, NULL);
}

Expr* builtin_integrate_singular_points(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* x; Expr* a; Expr* b;
    if (!nl_parse_spec(res->data.function.args[1], &x, &a, &b)) return NULL;

    PoleSet poles; poleset_init(&poles);
    nl_collect_poles(expr, x, a, b, &poles, NULL);
    /* The detector returns only the confirmed interior poles; the
     * undecidable flag is what makes the *integrator* bail, not this. */

    Expr** args = poles.n ? malloc(poles.n * sizeof(*args)) : NULL;
    for (size_t i = 0; i < poles.n; i++) args[i] = expr_copy(poles.pts[i]);
    Expr* out = expr_new_function(mk_sym("List"), args, poles.n);
    free(args);
    poleset_free(&poles);
    return out;
}

void integrate_newton_leibniz_init(void) {
    symtab_add_builtin("Integrate`NewtonLeibniz", builtin_integrate_newton_leibniz);
    symtab_get_def("Integrate`NewtonLeibniz")->attributes |=
        ATTR_PROTECTED;
    symtab_set_docstring("Integrate`NewtonLeibniz",
        "Integrate`NewtonLeibniz[f, {x, a, b}] evaluates the definite integral "
        "of f over x from a to b by the fundamental theorem of calculus: it "
        "antidifferentiates f, then forms F(b)-F(a) via the Limit engine, "
        "splitting at interior poles and taking one-sided limits so improper "
        "integrals converge and divergent ones return Infinity / "
        "ComplexInfinity / Indeterminate.  Returns unevaluated when the "
        "antiderivative is unknown or a pole position is undecidable.");

    symtab_add_builtin("Integrate`SingularPoints", builtin_integrate_singular_points);
    symtab_get_def("Integrate`SingularPoints")->attributes |=
        ATTR_PROTECTED;
    symtab_set_docstring("Integrate`SingularPoints",
        "Integrate`SingularPoints[expr, {x, a, b}] returns the sorted list of "
        "real poles of expr (roots of Denominator[Together[expr]]) lying "
        "strictly inside the interval (a, b).");
}
