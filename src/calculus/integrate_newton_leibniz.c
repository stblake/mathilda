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
           strcmp(e->data.function.head->data.symbol, name) == 0;
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
        return e->data.symbol == x->data.symbol;
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

/* -------------------------------------------------------------------------
 * Interior-pole detection.
 * ---------------------------------------------------------------------- */

typedef struct {
    Expr** pts;   /* owned copies of the interior pole expressions           */
    double* vals; /* their numeric values, used for ordering / dedup         */
    size_t n, cap;
    bool undecidable; /* a pole exists whose interior-ness we cannot decide  */
} PoleSet;

static void poleset_init(PoleSet* s) {
    s->pts = NULL; s->vals = NULL; s->n = 0; s->cap = 0; s->undecidable = false;
}

static void poleset_free(PoleSet* s) {
    for (size_t i = 0; i < s->n; i++) expr_free(s->pts[i]);
    free(s->pts); free(s->vals);
    poleset_init(s);
}

/* Insert `val` (an owned Expr, numeric value `v`) keeping the array sorted
 * ascending and free of duplicates.  Consumes `val` (frees it on dedup). */
static void poleset_insert(PoleSet* s, Expr* val, double v) {
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
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == x->data.symbol)
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

/* Collect the interior real poles of `expr` (as a function of `x`) into `s`.
 * A pole is a real root of Denominator[Together[expr]] lying strictly in
 * (a, b).  Sets s->undecidable if a pole exists whose classification cannot
 * be decided. */
static void nl_collect_poles(Expr* expr, Expr* x, Expr* a, Expr* b, PoleSet* s) {
    Expr* tog = eval_take(mk_fn1("Together", expr_copy(expr)));
    if (!tog) return;
    Expr* den = eval_take(mk_fn1("Denominator", tog));   /* tog consumed */
    if (!den) return;

    if (!contains_symbol(den, x)) { expr_free(den); return; }  /* no poles */

    /* Solve[den == 0, x, Reals] -- restrict to real roots up front. */
    Expr* eq    = mk_fn2("Equal", den, expr_new_integer(0));   /* den consumed */
    Expr* solve = mk_fn3("Solve", eq, expr_copy(x), mk_sym("Reals"));
    Expr* sols  = eval_take(solve);
    if (!sols) { s->undecidable = true; return; }

    if (!head_name_is(sols, "List")) {   /* unevaluated Solve / condition */
        s->undecidable = true;
        expr_free(sols);
        return;
    }

    for (size_t i = 0; i < sols->data.function.arg_count; i++) {
        Expr* val = nl_rule_rhs_for(sols->data.function.args[i], x);
        if (!val) { s->undecidable = true; continue; }

        double rv;
        int rk = nl_numeric(val, &rv);
        if (rk == 0) { s->undecidable = true; continue; }  /* symbolic root */
        if (rk != 1) continue;                             /* infinite root: skip */

        NLClass cls = nl_classify(rv, a, b);
        if (cls == NL_INTERIOR)      poleset_insert(s, expr_copy(val), rv);
        else if (cls == NL_UNKNOWN)  s->undecidable = true;
        /* NL_EXTERIOR: ignore */
    }
    expr_free(sols);
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

/* A boundary value is usable directly only if it is a finite closed form:
 * free of x, and free of any singular sentinel / unresolved Limit or Integrate. */
static bool nl_is_finite_closed(Expr* v, const Expr* x) {
    if (!v) return false;
    if (nl_contains_singular(v)) return false;
    if (head_name_is(v, "Interval")) return false;
    if (contains_head(v, "Limit") || contains_head(v, "Integrate")) return false;
    if (contains_symbol(v, x)) return false;
    return true;
}

/* Limit[F, x -> target, Direction -> dir]. */
static Expr* nl_limit(Expr* F, Expr* x, Expr* target, NLDir dir) {
    Expr* rule = mk_fn2("Rule", expr_copy(x), expr_copy(target));
    Expr* call;
    if (dir == NL_DIR_NONE) {
        call = mk_fn2("Limit", expr_copy(F), rule);
    } else {
        Expr* d = mk_fn2("Rule", mk_sym("Direction"),
                         expr_new_string(dir == NL_DIR_ABOVE ? "FromAbove"
                                                             : "FromBelow"));
        call = mk_fn3("Limit", expr_copy(F), rule, d);
    }
    return eval_take(call);
}

/* Evaluate F at `target`, approaching per `dir`.  Uses direct substitution
 * at ordinary finite points and the Limit engine at infinite / singular
 * points (or when substitution is itself singular). */
static Expr* nl_eval_at(Expr* F, Expr* x, Expr* target, NLDir dir) {
    arith_warnings_mute_push();
    Expr* out;
    if (nl_is_infinite(target)) {
        out = nl_limit(F, x, target, NL_DIR_NONE);
    } else {
        Expr* sub = eval_take(mk_fn2("ReplaceAll", expr_copy(F),
                                     mk_fn2("Rule", expr_copy(x),
                                            expr_copy(target))));
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
            out = nl_limit(F, x, target, dir);
        }
    }
    arith_warnings_mute_pop();
    return out;
}

/* -------------------------------------------------------------------------
 * Divergence reporting + continuity cross-check.
 * ---------------------------------------------------------------------- */

/* Emit the Mathematica-style Integrate::idiv warning to stderr, naming the
 * integrand and interval. */
static void nl_emit_idiv(Expr* f, Expr* a, Expr* b) {
    char* fs = expr_to_string(f);
    char* as = expr_to_string(a);
    char* bs = expr_to_string(b);
    fprintf(stderr,
            "Integrate::idiv: Integral of %s does not converge on {%s, %s}.\n",
            fs ? fs : "?", as ? as : "?", bs ? bs : "?");
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
    if (nl_numeric(V, &vd) != 1) return false;   /* V not a finite real */

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

/* -------------------------------------------------------------------------
 * Core Newton-Leibniz driver.
 * ---------------------------------------------------------------------- */

Expr* integrate_newton_leibniz_try(Expr* f, Expr* x, Expr* a, Expr* b,
                                   const char* method) {
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
    nl_collect_poles(f, x, a, b, &poles);
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

    /* 3. Breakpoints [a, poles..., b] and the telescoping sum. */
    size_t nbp = poles.n + 2;
    Expr** bp = malloc(nbp * sizeof(*bp));
    bp[0] = a;                                  /* borrowed */
    for (size_t i = 0; i < poles.n; i++) bp[i + 1] = poles.pts[i]; /* borrowed */
    bp[nbp - 1] = b;                            /* borrowed */

    Expr* total = expr_new_integer(0);
    bool failed = false;
    for (size_t i = 0; i + 1 < nbp && !failed; i++) {
        Expr* lo = nl_eval_at(F, x, bp[i],     NL_DIR_ABOVE);
        Expr* hi = nl_eval_at(F, x, bp[i + 1], NL_DIR_BELOW);
        if (!lo || !hi) {
            if (lo) expr_free(lo);
            if (hi) expr_free(hi);
            failed = true;
            break;
        }
        arith_warnings_mute_push();
        Expr* contrib = eval_take(mk_fn2("Plus", hi,
                                         mk_fn2("Times", expr_new_integer(-1), lo)));
        total = eval_take(mk_fn2("Plus", total, contrib));
        arith_warnings_mute_pop();
    }

    free(bp);
    poleset_free(&poles);

    /* Capture F's continuity-risk fingerprint before releasing it: a step head
     * or an inverse-trig node over a pole-bearing trig head means F could carry
     * an uncorrected branch jump and must be cross-checked. */
    bool F_triggers = nl_has_step_head(F, x) || nl_has_arctrig_over_pole(F, x);
    expr_free(F);

    if (failed || !total) { if (total) expr_free(total); return NULL; }

    /* 4a. Genuine divergence: an interior pole (or a boundary blow-up) leaves a
     * divergence sentinel.  Warn as Mathematica does and leave the integral
     * unevaluated. */
    if (nl_is_infinite(total) || is_indeterminate_sym(total)) {
        nl_emit_idiv(f, a, b);
        expr_free(total);
        return NULL;
    }

    /* 4b. An unresolved Limit/Integrate, a stray x, an Interval, or a residual
     * singular sentinel buried inside the result means the engine could not
     * decide -- leave the integral unevaluated. */
    if (contains_head(total, "Limit") || contains_head(total, "Integrate") ||
        head_name_is(total, "Interval") || contains_symbol(total, x) ||
        nl_contains_singular(total)) {
        expr_free(total);
        return NULL;
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
    nl_collect_poles(expr, x, a, b, &poles);
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
        ATTR_PROTECTED | ATTR_READPROTECTED;
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
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`SingularPoints",
        "Integrate`SingularPoints[expr, {x, a, b}] returns the sorted list of "
        "real poles of expr (roots of Denominator[Together[expr]]) lying "
        "strictly inside the interval (a, b).");
}
