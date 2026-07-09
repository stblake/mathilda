/* integrate_diffunderint.c
 *
 * Definite integration by differentiation under the integral sign (Leibniz
 * rule / "Feynman's trick").  See integrate_diffunderint.h for the overview.
 *
 * Method (unified first-order-ODE framing, Boulnois arXiv:2308.09619):
 * for a definite integral I(p) = Integrate[f(x,p), {x,a,b}] with a free
 * parameter p, we solve the first-order ODE I'(p) = lambda(p) I(p) + M(p):
 *
 *   Stage A (lambda = 0, pure quadrature): differentiate the integrand,
 *     evaluate the (simpler) inner integral J(p) = Integrate[D[f,p], {x,a,b}]
 *     with the existing engine, integrate J(p) back over the parameter, and fix
 *     the constant with an EXACT base value I(p0).
 *
 * Verification is symbolic and correct-by-construction (PossibleZeroQ[D[I,p]-J]
 * plus an exact base).  No NIntegrate anywhere (project rule): the conditional-
 * convergence pitfall (Conrad section 12) is caught for free -- a non-integrable
 * D[f,p] makes the inner Integrate fail to close, so that parameter is skipped.
 */

#include "integrate_diffunderint.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* arith_warnings_mute_push/pop */
#include "sym_names.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifdef DIUI_DEBUG
#include <time.h>
static double diui_ms(void) {
    static clock_t t0 = 0;
    if (!t0) t0 = clock();
    return 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
}
#endif

/* -------------------------------------------------------------------------
 * Recursion guard.  The method recurses into the full Integrate cascade (for
 * the inner integral J and for directly-integrable base values), which can
 * re-enter DiffUnderInt.  Depth 2 is required: e.g. Integrate[Exp[-x^2]
 * Sin[a x]/x, ...] differentiates to a Gaussian cosine transform that is itself
 * a DiffUnderInt target.  The indefinite parameter-integration Integrate[J, p]
 * is definite-only-agnostic and cannot re-enter this method.
 * ---------------------------------------------------------------------- */
static int diui_depth = 0;
#define DIUI_MAX_DEPTH 1

/* -------------------------------------------------------------------------
 * Small expression-construction / evaluation helpers (mirrors the idiom in
 * integrate_newton_leibniz.c).
 * ---------------------------------------------------------------------- */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_int(long v)        { return expr_new_integer((int64_t)v); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c }, 3);
}

/* Small unevaluated arithmetic constructors (consume their args), used to build
 * the closed-form family expressions before a single Simplify at the end. */
static Expr* t_add(Expr* a, Expr* b) { return mk_fn2("Plus", a, b); }
static Expr* t_mul(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* t_neg(Expr* a)          { return mk_fn2("Times", mk_int(-1), a); }
static Expr* t_pow(Expr* a, long n)  { return mk_fn2("Power", a, mk_int(n)); }
static Expr* t_rat(long p, long q)   { return mk_fn2("Rational", mk_int(p), mk_int(q)); }

/* Evaluate `call`, free the call expression, return the (owned) result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* ev1/ev2 build `name[...]` from freshly-owned args and evaluate; the args are
 * consumed. */
static Expr* ev1(const char* name, Expr* a)            { return eval_take(mk_fn1(name, a)); }
static Expr* ev2(const char* name, Expr* a, Expr* b)   { return eval_take(mk_fn2(name, a, b)); }

/* True iff `e` is the compound `name[...]` (by head name). */
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

/* True iff the bare symbol `name` occurs anywhere in `e`. */
static bool contains_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol_name(e->data.function.args[i], name)) return true;
    return false;
}

/* True iff the symbol `x` occurs anywhere in `e`. */
static bool contains_symbol(const Expr* e, const Expr* x) {
    return contains_symbol_name(e, x->data.symbol);
}

/* True iff `e` contains Power[x, p] with p a constant other than 0 or 1 -- a
 * genuinely nonlinear / singular power of x (x^2, x^(-2), Sqrt[x], ...). */
static bool has_nonlinear_x_power(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == x->data.symbol) {
            if (ex->type == EXPR_INTEGER &&
                (ex->data.integer >= 2 || ex->data.integer <= -1)) return true;
            if (ex->type == EXPR_REAL) return true;
            if (ex->type == EXPR_FUNCTION &&
                head_name_is(ex, "Rational")) return true;   /* Sqrt[x] etc. */
        }
    }
    if (has_nonlinear_x_power(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_nonlinear_x_power(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` contains an exponential Exp[g] / E^g whose exponent g is
 * nonlinear or singular in x (a Gaussian e^{-x^2}, e^{-a^2/x^2}, ...).  The
 * indefinite/definite integrator currently HANGS on these forms, so the
 * DiffUnderInt method must decline them up front rather than spawn an inner
 * integral that never returns. */
static bool contains_gaussian_exp(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* arg = NULL;
    if (head_name_is(e, "Exp") && e->data.function.arg_count == 1)
        arg = e->data.function.args[0];
    else if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base->type == EXPR_SYMBOL && strcmp(base->data.symbol, "E") == 0)
            arg = e->data.function.args[1];
    }
    if (arg && contains_symbol(arg, x) && has_nonlinear_x_power(arg, x))
        return true;
    if (contains_gaussian_exp(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_gaussian_exp(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` contains a forward trig / hyperbolic function OF x (Sin[.. x ..],
 * Cos, Tan, ..., Sinh, ...).  The general integrator hangs on trig integrands
 * over finite periods, so such forms must be routed to a family, never the
 * engine. */
static bool has_trig_of_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    static const char* T[] = { "Sin","Cos","Tan","Cot","Sec","Csc",
                               "Sinh","Cosh","Tanh","Coth","Sech","Csch" };
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        for (size_t i = 0; i < sizeof(T)/sizeof(T[0]); i++)
            if (strcmp(h, T[i]) == 0 && contains_symbol(e, x)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_trig_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` contains a radical of x: Power[base, e] with a fractional-CONSTANT
 * exponent and base depending on x (Sqrt[1-x^2], (a^2-x^2)^(1/2), ...).  A
 * symbolic exponent like x^a is NOT a radical and stays engine-safe. */
static bool has_radical_of_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        bool frac = (ex->type == EXPR_REAL) ||
                    (ex->type == EXPR_FUNCTION && head_name_is(ex, "Rational"));
        if (frac && contains_symbol(base, x)) return true;
    }
    if (has_radical_of_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_radical_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* True unless `e` carries a non-finite / undecided marker (infinities,
 * Indeterminate, an unresolved Integrate, ...): the base value and the
 * antiderivative-at-base must be genuine finite closed forms. */
static bool is_finite_value(const Expr* e) {
    if (!e) return false;
    if (contains_head(e, "DirectedInfinity")) return false;
    if (contains_head(e, "Integrate"))        return false;
    static const char* bad[] = { "Indeterminate", "ComplexInfinity", "Infinity",
                                 "Underflow", "Overflow", "Null", "$Aborted" };
    for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++)
        if (contains_symbol_name(e, bad[i])) return false;
    return true;
}

/* Strict symbolic zero test: Simplify[e] reduces to the literal 0.
 *
 * We deliberately do NOT use PossibleZeroQ here.  PossibleZeroQ's shrinkage-
 * trend heuristic false-positives on decaying expressions (e.g.
 * PossibleZeroQ[-E^(-a x)] === True because the magnitude shrinks toward 0 as
 * the sample point grows), which would both mis-skip a valid parameter (a
 * nonzero D[f,p]) and, worse, mis-accept a wrong result during verification.
 * A literal Simplify-to-0 is a genuine proof of the identities involved (the
 * targets are elementary rational/log/exp/inverse-trig forms), so it is both
 * safe for acceptance and correct for the zero-integrand base test. */
static bool is_zero_q(const Expr* e) {
    Expr* s = ev1("Simplify", expr_copy((Expr*)e));
    bool z = s && ((s->type == EXPR_INTEGER && s->data.integer == 0) ||
                   (s->type == EXPR_REAL && s->data.real == 0.0));
    if (s) expr_free(s);
    return z;
}

/* Simplify[e] or Simplify[e, assumptions].  The assumptions (a>0, ...) are what
 * collapse the radical/inverse-hyperbolic forms the family evaluators emit
 * (Sqrt[1/a^2] -> 1/a, ArcCoth[...] -> Log[...]) into clean closed forms; a
 * clean J is also essential for the subsequent parameter-integration to stay
 * fast.  Borrows `e`; `assumptions` may be NULL. */
static Expr* simplify_with(const Expr* e, const Expr* assumptions) {
    if (assumptions)
        return ev2("Simplify", expr_copy((Expr*)e), expr_copy((Expr*)assumptions));
    return ev1("Simplify", expr_copy((Expr*)e));
}

/* is_zero_q with assumptions. */
static bool is_zero_with(const Expr* e, const Expr* assumptions) {
    Expr* s = simplify_with(e, assumptions);
    bool z = s && ((s->type == EXPR_INTEGER && s->data.integer == 0) ||
                   (s->type == EXPR_REAL && s->data.real == 0.0));
    if (s) expr_free(s);
    return z;
}

/* ReplaceAll[e, var -> val], evaluated.  Borrows all three. */
static Expr* subst(const Expr* e, const Expr* var, const Expr* val) {
    return ev2("ReplaceAll", expr_copy((Expr*)e),
               mk_fn2("Rule", expr_copy((Expr*)var), expr_copy((Expr*)val)));
}

/* D[e, var], evaluated.  Borrows both. */
static Expr* deriv(const Expr* e, const Expr* var) {
    return ev2("D", expr_copy((Expr*)e), expr_copy((Expr*)var));
}

/* Definite inner integral Integrate[g, {x,a,b}] (+ Assumptions when present),
 * evaluated by the full engine.  Borrows all. */
static Expr* integrate_definite_of(const Expr* g, const Expr* x, const Expr* a,
                                   const Expr* b, const Expr* assumptions) {
    Expr* spec = mk_fn3("List", expr_copy((Expr*)x), expr_copy((Expr*)a),
                        expr_copy((Expr*)b));
    Expr* call;
    if (assumptions)
        call = mk_fn3("Integrate", expr_copy((Expr*)g), spec,
                      mk_fn2("Rule", mk_sym("Assumptions"),
                             expr_copy((Expr*)assumptions)));
    else
        call = mk_fn2("Integrate", expr_copy((Expr*)g), spec);
    /* Callers gate this to engine-safe integrands only (see inner_definite);
     * TimeConstrained is not used because it does not bound a nested evaluate. */
    return eval_take(call);
}

/* Gaussian parameter back-integration -> Erf (the engine cannot integrate a
 * Gaussian).  Forward-declared here; defined with the Gaussian family below. */
static Expr* integrate_gaussian_param(const Expr* J, const Expr* p);

/* Indefinite integral Integrate[J, p] over the parameter.  Borrows both.  J is a
 * closed form from a family (rational / arctan / log-like), so integrating it
 * over the parameter is elementary and fast -- EXCEPT a Gaussian J = c e^{-k p^2}
 * (from the Gaussian cosine-moment family), whose antiderivative is an Erf the
 * engine does not produce; that case is handled directly. */
static Expr* integrate_over_param(const Expr* J, const Expr* p) {
    Expr* gp = integrate_gaussian_param(J, p);
    if (gp) return gp;
    return ev2("Integrate", expr_copy((Expr*)J), expr_copy((Expr*)p));
}

/* Value of the antiderivative G at p = p0.  Direct substitution when it yields a
 * finite value; otherwise the Limit engine (handles 0, Infinity, removable
 * singularities).  Borrows all. */
static Expr* eval_at_param(const Expr* G, const Expr* p, const Expr* p0) {
    Expr* s = subst(G, p, p0);
    if (s && is_finite_value(s)) {
        Expr* ss = ev1("Simplify", s);
        if (ss && is_finite_value(ss)) return ss;
        if (ss) expr_free(ss);
    } else if (s) {
        expr_free(s);
    }
    Expr* lim = ev2("Limit", expr_copy((Expr*)G),
                    mk_fn2("Rule", expr_copy((Expr*)p), expr_copy((Expr*)p0)));
    if (lim && is_finite_value(lim)) return lim;
    if (lim) expr_free(lim);
    return NULL;
}

/* -------------------------------------------------------------------------
 * Parameter collection + assumptions (minimal, self-contained; the residue
 * method has an equivalent but file-static version).
 * ---------------------------------------------------------------------- */

typedef struct { const char* sym; double lo, hi; } ParamBound;

/* Real machine value of a concrete numeric expression, else false (symbolic). */
static bool numeric_double(const Expr* e, double* out) {
    Expr* n = ev1("N", expr_copy((Expr*)e));
    bool ok = false;
    if (n) {
        if (n->type == EXPR_INTEGER) { *out = (double)n->data.integer; ok = true; }
        else if (n->type == EXPR_REAL) { *out = n->data.real; ok = true; }
        expr_free(n);
    }
    return ok;
}

/* Collect distinct free-parameter symbols of `e` (not x, not I, not a numeric
 * constant like Pi/E) into `pb` (up to `cap`).  Returns the count. */
static size_t collect_params(const Expr* e, const Expr* x, ParamBound* pb,
                             size_t cap, size_t n) {
    if (!e) return n;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == x->data.symbol) return n;
        if (strcmp(e->data.symbol, "I") == 0) return n;
        for (size_t i = 0; i < n; i++)
            if (pb[i].sym == e->data.symbol) return n;
        double tmp;
        if (numeric_double(e, &tmp)) return n;       /* Pi, E, EulerGamma, ... */
        if (n < cap) { pb[n].sym = e->data.symbol; pb[n].lo = -HUGE_VAL;
                       pb[n].hi = HUGE_VAL; n++; }
        return n;
    }
    if (e->type != EXPR_FUNCTION) return n;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n = collect_params(e->data.function.args[i], x, pb, cap, n);
    return n;
}

static void bound_apply(ParamBound* pb, size_t np, const char* var, double c, int dir) {
    for (size_t i = 0; i < np; i++) {
        if (pb[i].sym != var) continue;
        if (dir > 0) { if (c > pb[i].lo) pb[i].lo = c; }
        else         { if (c < pb[i].hi) pb[i].hi = c; }
        return;
    }
}

static void relation_apply(ParamBound* pb, size_t np, const char* op,
                           Expr* L, Expr* R) {
    bool lt = (strcmp(op, "Less") == 0 || strcmp(op, "LessEqual") == 0);
    bool gt = (strcmp(op, "Greater") == 0 || strcmp(op, "GreaterEqual") == 0);
    if (!lt && !gt) return;
    double c;
    if (L->type == EXPR_SYMBOL && numeric_double(R, &c))
        bound_apply(pb, np, L->data.symbol, c, lt ? -1 : +1);
    else if (R->type == EXPR_SYMBOL && numeric_double(L, &c))
        bound_apply(pb, np, R->data.symbol, c, lt ? +1 : -1);
}

/* Tighten parameter bounds from an Assumptions fact (And/List conjunctions,
 * chained Inequality, and the four ordered binary relations). */
static void absorb_fact(ParamBound* pb, size_t np, Expr* fact) {
    if (!fact || fact->type != EXPR_FUNCTION) return;
    if (fact->data.function.head->type != EXPR_SYMBOL) return;
    const char* h = fact->data.function.head->data.symbol;
    size_t ac = fact->data.function.arg_count;
    if (strcmp(h, "And") == 0 || strcmp(h, "List") == 0) {
        for (size_t i = 0; i < ac; i++) absorb_fact(pb, np, fact->data.function.args[i]);
        return;
    }
    if (strcmp(h, "Inequality") == 0 && ac >= 3) {
        for (size_t i = 0; i + 2 < ac; i += 2) {
            Expr* opsym = fact->data.function.args[i + 1];
            if (opsym->type != EXPR_SYMBOL) continue;
            relation_apply(pb, np, opsym->data.symbol,
                           fact->data.function.args[i], fact->data.function.args[i + 2]);
        }
        return;
    }
    if (ac == 2) relation_apply(pb, np, h,
                                fact->data.function.args[0], fact->data.function.args[1]);
}

/* A generic representative strictly inside [lo, hi] (deterministic).  Used only
 * to read off the SIGN of a real part for a convergence gate, never for a
 * numeric answer. */
static double pick_rep(double lo, double hi, size_t idx) {
    static const double SEEDS[] = { 1.3172, 2.7391, 0.6180, 3.1490,
                                    1.4213, 2.2360, 0.8177 };
    const size_t NS = sizeof(SEEDS) / sizeof(SEEDS[0]);
    double seed = SEEDS[idx % NS];
    bool lo_fin = lo > -HUGE_VAL, hi_fin = hi < HUGE_VAL;
    if (lo_fin && hi_fin) return lo + (hi - lo) * 0.5182;
    if (lo_fin)           return lo + seed;
    if (hi_fin)           return hi - seed;
    return seed;
}

/* -------------------------------------------------------------------------
 * Fast inner definite-integral families.
 *
 * The general integrator is too slow / incomplete on the parameter-dependent
 * forms that Feynman's trick produces, so DiffUnderInt evaluates the standard
 * families itself with closed-form formulas, using the engine only for the fast
 * complex *algebra* (never for the integral).  The general engine remains the
 * final fallback (TimeConstrained-bounded).
 * ---------------------------------------------------------------------- */

/* a == literal 0 */
static bool is_zero_expr(const Expr* a) {
    return a && a->type == EXPR_INTEGER && a->data.integer == 0;
}
/* b == +Infinity (symbol Infinity or DirectedInfinity[1]) */
static bool is_pos_inf(const Expr* b) {
    if (!b) return false;
    if (b->type == EXPR_SYMBOL && strcmp(b->data.symbol, "Infinity") == 0) return true;
    if (head_name_is(b, "DirectedInfinity") && b->data.function.arg_count == 1) {
        Expr* d = b->data.function.args[0];
        return d->type == EXPR_INTEGER && d->data.integer == 1;
    }
    return false;
}

/* a == -Infinity (DirectedInfinity[-1] or Times[-1, Infinity]) */
static bool is_neg_inf(const Expr* a) {
    if (!a) return false;
    if (head_name_is(a, "DirectedInfinity") && a->data.function.arg_count == 1) {
        Expr* d = a->data.function.args[0];
        return d->type == EXPR_INTEGER && d->data.integer == -1;
    }
    if (head_name_is(a, "Times") && a->data.function.arg_count == 2) {
        Expr* c = a->data.function.args[0];
        return c->type == EXPR_INTEGER && c->data.integer == -1 &&
               is_pos_inf(a->data.function.args[1]);
    }
    return false;
}

/* Structural degree of a MONOMIAL in x: 0 if x-free, else sum of x-powers over a
 * Times, k for Power[x,k] (k a nonneg integer), 1 for x.  Returns -1 if `e` is
 * not a clean monomial with a nonnegative integer power of x. */
static long x_monomial_degree(const Expr* e, const Expr* x) {
    if (!contains_symbol(e, x)) return 0;
    if (e->type == EXPR_SYMBOL && e->data.symbol == x->data.symbol) return 1;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == x->data.symbol &&
            ex->type == EXPR_INTEGER && ex->data.integer >= 0)
            return (long)ex->data.integer;
        return -1;
    }
    if (head_name_is(e, "Times")) {
        long tot = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            long d = x_monomial_degree(e->data.function.args[i], x);
            if (d < 0) return -1;
            tot += d;
        }
        return tot;
    }
    return -1;
}

/* True iff Re(-alpha) > 0 at a representative point of the parameter box -- the
 * convergence condition for ∫₀^∞ x^n e^{alpha x} dx.  Reading only the sign at a
 * single interior point is sound because the assumptions pin alpha's real part
 * to one connected side. */
static bool alpha_re_neg_value(const Expr* alpha, const ParamBound* pb,
                               size_t np, double* out) {
    Expr* na = ev1("Simplify", mk_fn2("Times", mk_int(-1), expr_copy((Expr*)alpha)));
    for (size_t i = 0; i < np && na; i++) {
        Expr* rep = expr_new_real(pick_rep(pb[i].lo, pb[i].hi, i));
        Expr* sym = mk_sym(pb[i].sym);
        Expr* n2  = subst(na, sym, rep);
        expr_free(sym); expr_free(rep); expr_free(na);
        na = n2;
    }
    if (!na) return false;
    Expr* re = ev1("Re", na);
    bool ok = re && numeric_double(re, out);
    if (re) expr_free(re);
    return ok;
}

/* Strict decay Re(alpha) < 0: convergence for ∫₀^∞ x^n e^{alpha x} dx (n>=0). */
static bool alpha_decays(const Expr* alpha, const ParamBound* pb, size_t np) {
    double v; return alpha_re_neg_value(alpha, pb, np, &v) && v > 1e-9;
}

/* Non-growth Re(alpha) <= 0: allows the pure-oscillatory case (Re(alpha)=0) used
 * by the sinc / x^{-k} forms, whose convergence comes from the regularization
 * rather than exponential decay (the s-integral pole at s=alpha then sits off
 * the positive real s-axis). */
static bool alpha_no_growth(const Expr* alpha, const ParamBound* pb, size_t np) {
    double v; return alpha_re_neg_value(alpha, pb, np, &v) && v > -1e-6;
}

/* Collect the multiplicative leaf factors of a product `T`, flattening nested
 * Times (Expand can leave Times[c, Times[x^-1, E^..]] unflattened).  Returns the
 * count written into out[] (up to cap); factors are borrowed. */
static size_t collect_factors(Expr* T, Expr** out, size_t cap, size_t n) {
    if (head_name_is(T, "Times")) {
        for (size_t i = 0; i < T->data.function.arg_count && n < cap; i++)
            n = collect_factors(T->data.function.args[i], out, cap, n);
        return n;
    }
    if (n < cap) out[n++] = T;
    return n;
}

/* Laplace / Fourier half-line: ∫₀^∞ g dx where, after TrigToExp + Expand, every
 * term is c * x^n * e^{alpha x} with a nonnegative integer n and Re(alpha) < 0.
 * Each term integrates to c * e^{alpha0} * n! / (-alpha)^{n+1} (alpha0 the x-free
 * part of the exponent); conjugate pairs recombine to a real closed form under
 * Simplify.  Returns NULL if the integrand is not of this form or a term fails
 * the decay (convergence) gate.  Covers e.g. ∫₀^∞ e^{-a x} cos(b x) dx. */
static Expr* laplace_halfline(const Expr* g, const Expr* x,
                              const ParamBound* pb, size_t np,
                              const Expr* assumptions) {
    Expr* gg = ev1("Expand", ev1("TrigToExp", expr_copy((Expr*)g)));
    if (!gg) return NULL;
    if (contains_gaussian_exp(gg, x)) { expr_free(gg); return NULL; }

    /* Additive terms. */
    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(gg, "Plus")) {
        nt = gg->data.function.arg_count;
        terms = gg->data.function.args;
    } else { nt = 1; single[0] = gg; terms = single; }

    Expr* total = mk_int(0);
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        /* Separate exponential factors (accumulate their arguments) from the
         * rest of the product. */
        Expr* exparg = mk_int(0);          /* running sum of exp arguments */
        Expr* rest   = mk_int(1);          /* running product of the remainder */
        Expr* facbuf[64];
        size_t nf = collect_factors(T, facbuf, 64, 0);
        Expr** facs = facbuf;
        bool saw_exp = false;
        for (size_t i = 0; i < nf; i++) {
            Expr* F = facs[i];
            const Expr* ea = NULL;
            if (head_name_is(F, "Exp") && F->data.function.arg_count == 1)
                ea = F->data.function.args[0];
            else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                     F->data.function.args[0]->type == EXPR_SYMBOL &&
                     strcmp(F->data.function.args[0]->data.symbol, "E") == 0)
                ea = F->data.function.args[1];
            if (ea) { saw_exp = true;
                      exparg = mk_fn2("Plus", exparg, expr_copy((Expr*)ea)); }
            else    { rest   = mk_fn2("Times", rest, expr_copy(F)); }
        }
        if (!saw_exp) { ok = false; expr_free(exparg); expr_free(rest); break; }

        exparg = ev1("Simplify", exparg);
        Expr* alpha = deriv(exparg, x);                 /* linear coeff in x */
        Expr* d2    = ev1("Simplify", deriv(alpha, x)); /* must vanish (linear) */
        if (!d2 || !is_zero_q(d2)) { if (d2) expr_free(d2);
            ok = false; expr_free(alpha); expr_free(exparg); expr_free(rest); break; }
        expr_free(d2);
        Expr* zero = mk_int(0);
        Expr* off = ev1("Simplify", subst(exparg, x, zero));
        expr_free(zero); expr_free(exparg);

        Expr* restc = ev1("Simplify", rest);
        long n = x_monomial_degree(restc, x);
        Expr* one = mk_int(1);
        Expr* c = (n >= 0) ? subst(restc, x, one) : NULL;
        expr_free(one);
        if (n < 0 || !c || contains_symbol(c, x) || !alpha_decays(alpha, pb, np)) {
            if (c) expr_free(c);
            ok = false; expr_free(alpha); expr_free(off); expr_free(restc); break;
        }
        expr_free(restc);
        /* term integral = c * e^{off} * n! / (-alpha)^{n+1} */
        Expr* neg_alpha = mk_fn2("Times", mk_int(-1), alpha);   /* consumes alpha */
        Expr* ti = mk_fn2("Times",
                     mk_fn2("Times", c, mk_fn1("Exp", off)),
                     mk_fn2("Times", mk_fn1("Factorial", mk_int(n)),
                            mk_fn2("Power", neg_alpha, mk_int(-(n + 1)))));
        total = mk_fn2("Plus", total, ti);
    }
    if (head_name_is(gg, "Plus")) { /* terms borrowed from gg */ }
    expr_free(gg);
    if (!ok) { expr_free(total); return NULL; }

    Expr* res = simplify_with(total, assumptions);
    expr_free(total);
    if (res && is_finite_value(res) && !contains_symbol(res, x)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* Signed x-degree of a monomial (allows negative powers, x^(-1), x^(-2)).
 * Sets *ok = false if `e` is not a clean integer-power monomial in x. */
static long x_deg_signed(const Expr* e, const Expr* x, bool* ok) {
    if (!contains_symbol(e, x)) return 0;
    if (e->type == EXPR_SYMBOL && e->data.symbol == x->data.symbol) return 1;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == x->data.symbol &&
            ex->type == EXPR_INTEGER)
            return (long)ex->data.integer;
        *ok = false; return 0;
    }
    if (head_name_is(e, "Times")) {
        long tot = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            tot += x_deg_signed(e->data.function.args[i], x, ok);
        return tot;
    }
    *ok = false; return 0;
}

static Expr* rational_halfline(const Expr* g, const Expr* x,
                               const ParamBound* pb, size_t np,
                               const Expr* assumptions);   /* forward decl */
static Expr* rational_halfline_general(const Expr* g, const Expr* s,
                                       const ParamBound* pb, size_t np,
                                       const Expr* assumptions);   /* forward decl */

/* Sinc / Frullani half-line: ∫₀^∞ H(x)/x^k dx (k >= 1), for an integrand that
 * after TrigToExp+Expand is a sum of c * x^n * e^{alpha x} terms with some n < 0.
 * Uses 1/x^k = ∫_0^∞ s^{k-1}/(k-1)! e^{-s x} ds to write the value as
 *   ∫_0^∞ (s^{k-1}/(k-1)!) M(s) ds,   M(s) = ∫₀^∞ H(x) e^{-s x} dx = sum of
 *   c (n+k)! / (s - alpha)^{n+k+1},
 * a RATIONAL function of s whose half-line integral the engine evaluates safely
 * (antiderivative is Log/ArcTan).  Covers e.g. ∫₀^∞ Sin[q x]/x dx = Pi/2,
 * ∫₀^∞ e^{-p x} Sin[q x]/x dx = ArcTan[q/p], and the /x^2 forms (#6,#14,#23). */
static Expr* laplace_sinc_halfline(const Expr* g, const Expr* x,
                                   const ParamBound* pb, size_t np,
                                   const Expr* assumptions) {
    Expr* gg = ev1("Expand", ev1("TrigToExp", expr_copy((Expr*)g)));
    if (!gg) return NULL;
    if (contains_gaussian_exp(gg, x)) { expr_free(gg); return NULL; }

    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(gg, "Plus")) { nt = gg->data.function.arg_count; terms = gg->data.function.args; }
    else { nt = 1; single[0] = gg; terms = single; }
    if (nt > 64) { expr_free(gg); return NULL; }

    Expr* C[64]; Expr* AL[64]; long N[64]; size_t m = 0;
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        Expr* exparg = mk_int(0); Expr* rest = mk_int(1);
        Expr* facbuf[64];
        size_t nf = collect_factors(T, facbuf, 64, 0);
        Expr** facs = facbuf;
        bool saw_exp = false;
        for (size_t i = 0; i < nf; i++) {
            Expr* F = facs[i];
            const Expr* ea = NULL;
            if (head_name_is(F, "Exp") && F->data.function.arg_count == 1) ea = F->data.function.args[0];
            else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                     F->data.function.args[0]->type == EXPR_SYMBOL &&
                     strcmp(F->data.function.args[0]->data.symbol, "E") == 0) ea = F->data.function.args[1];
            if (ea) { saw_exp = true; exparg = mk_fn2("Plus", exparg, expr_copy((Expr*)ea)); }
            else    { rest = mk_fn2("Times", rest, expr_copy(F)); }
        }
        if (!saw_exp) { ok = false; expr_free(exparg); expr_free(rest); break; }
        exparg = ev1("Simplify", exparg);
        Expr* alpha = deriv(exparg, x);
        Expr* d2 = ev1("Simplify", deriv(alpha, x));
        bool lin = d2 && is_zero_q(d2); if (d2) expr_free(d2);
        Expr* zero = mk_int(0);
        Expr* off  = ev1("Simplify", subst(exparg, x, zero));
        expr_free(zero); expr_free(exparg);
        Expr* restc = ev1("Simplify", rest);
        bool degok = true;
        long n = x_deg_signed(restc, x, &degok);
        Expr* one = mk_int(1);
        Expr* c = degok ? subst(restc, x, one) : NULL;   /* coefficient c e^{off} */
        expr_free(one); expr_free(restc);
        bool ang = alpha_no_growth(alpha, pb, np);
#ifdef DIUI_DEBUG
        fprintf(stderr, "DIUI:     sinc term: saw_exp lin=%d degok=%d n=%ld c=%d cx=%d ang=%d\n",
                (int)lin, (int)degok, n, c!=NULL, c?(int)contains_symbol(c,x):-1, (int)ang);
#endif
        if (!lin || !degok || !c || contains_symbol(c, x) || !ang) {
            if (c) expr_free(c); expr_free(alpha); expr_free(off); ok = false; break;
        }
        C[m]  = mk_fn2("Times", c, mk_fn1("Exp", off));   /* c * e^{off} */
        AL[m] = alpha; N[m] = n; m++;
    }
    expr_free(gg);
#ifdef DIUI_DEBUG
    fprintf(stderr, "DIUI:   sinc: m=%zu ok=%d\n", m, (int)ok);
#endif
    if (!ok || m == 0) { for (size_t j = 0; j < m; j++) { expr_free(C[j]); expr_free(AL[j]); } return NULL; }

    /* Only the pure k=1 case (all terms x^{-1} e^{alpha x}) is handled.  The
     * value is ∫_0^∞ M(s) ds with M(s) = sum c_j/(s - alpha_j), the Laplace
     * transform of x*integrand -- a REAL rational function of s.  Rather than
     * integrate M via an antiderivative + Limit[..,s->Infinity] (which the engine
     * cannot resolve for a symbolic parameter sign), we hand M to the even-
     * rational half-line evaluator: for the oscillatory case (alpha pure
     * imaginary) M(s) is even in s, e.g. Sin[q x]/x -> M = q/(s^2+q^2),
     * ∫_0^∞ = Pi/2.  Non-even M (a genuine e^{-p x} decay, p>0) is deferred. */
    bool all_m1 = true;
    for (size_t j = 0; j < m; j++) if (N[j] != -1) { all_m1 = false; break; }
    if (!all_m1) { for (size_t j = 0; j < m; j++) { expr_free(C[j]); expr_free(AL[j]); } return NULL; }

    Expr* s = mk_sym("$diuiSig$");
    Expr* M = mk_int(0);
    for (size_t j = 0; j < m; j++) {
        Expr* pole = mk_fn2("Plus", expr_copy(s), mk_fn2("Times", mk_int(-1), AL[j]));
        M = mk_fn2("Plus", M, mk_fn2("Times", C[j], mk_fn2("Power", pole, mk_int(-1))));
    }
    Expr* Ms = simplify_with(M, assumptions);
    expr_free(M);
    /* Oscillatory M (pure imaginary poles) is even in s -> even family; a genuine
     * decay (Re(alpha) < 0) makes M non-even -> the general real-rational family,
     * which returns real ArcTan/Log directly. */
    Expr* res = Ms ? rational_halfline(Ms, s, pb, np, assumptions) : NULL;
    if (!res && Ms) res = rational_halfline_general(Ms, s, pb, np, assumptions);
#ifdef DIUI_DEBUG
    fprintf(stderr, "DIUI:   sinc: rational_halfline(M) -> %s\n", res ? "HIT" : "miss");
#endif
    if (Ms) expr_free(Ms);
    expr_free(s);
    if (res && is_finite_value(res) && !contains_symbol(res, x)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* True iff `e` is > 0 at a representative point of the parameter box (a sign
 * gate, used for a pole parameter d = beta/alpha that must be positive for
 * ∫₀^∞ v^{-1/2}/(v+d)^m dv to converge). */
static bool real_positive(const Expr* e, const ParamBound* pb, size_t np) {
    Expr* s = ev1("Simplify", expr_copy((Expr*)e));
    for (size_t i = 0; i < np && s; i++) {
        Expr* rep = expr_new_real(pick_rep(pb[i].lo, pb[i].hi, i));
        Expr* sym = mk_sym(pb[i].sym);
        Expr* s2  = subst(s, sym, rep);
        expr_free(sym); expr_free(rep); expr_free(s);
        s = s2;
    }
    double v; bool ok = s && numeric_double(s, &v) && v > 1e-9;
    if (s) expr_free(s);
    return ok;
}

/* Even rational half-line: ∫₀^∞ R(x) dx for an even rational R with no real
 * poles (denominator a product of (x^2 + d_k)^{m_k}, d_k > 0) and a >= 2 degree
 * drop.  Substituting v = x^2 gives (1/2) ∫₀^∞ v^{-1/2} R~(v) dv; a partial-
 * fraction split of R~ in v reduces each simple factor to the Beta integral
 * ∫₀^∞ v^{-1/2}/(v+d)^m dv = sqrt(Pi) Gamma(m-1/2)/Gamma(m) d^{1/2-m}.
 * Covers e.g. ∫₀^∞ dx/((1+a^2 x^2)(1+x^2)) = Pi/(2(1+a)). */
static Expr* rational_halfline(const Expr* g, const Expr* x,
                               const ParamBound* pb, size_t np,
                               const Expr* assumptions) {
    /* Evenness: g(x) - g(-x) === 0. */
    Expr* negx = mk_fn2("Times", mk_int(-1), expr_copy((Expr*)x));
    Expr* gneg = subst(g, x, negx);
    expr_free(negx);
    Expr* diff = mk_fn2("Plus", expr_copy((Expr*)g),
                        mk_fn2("Times", mk_int(-1), gneg));
    bool even = is_zero_q(diff);
    expr_free(diff);
    if (!even) return NULL;

    /* R~(v) = Together[g /. x -> Sqrt[v]]; must be rational in v (no x, no
     * fractional power of v). */
    Expr* v = mk_sym("$diuiRatV$");
    Expr* half = mk_fn2("Rational", mk_int(1), mk_int(2));
    Expr* sqrtv = mk_fn2("Power", expr_copy(v), half);
    Expr* Rv = ev1("Together", subst(g, x, sqrtv));
    expr_free(sqrtv);
    if (!Rv || contains_symbol(Rv, x)) { if (Rv) expr_free(Rv); expr_free(v); return NULL; }

    Expr* ap = ev2("Apart", Rv, expr_copy(v));         /* consumes Rv */
    if (!ap) { expr_free(v); return NULL; }

    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(ap, "Plus")) { nt = ap->data.function.arg_count; terms = ap->data.function.args; }
    else { nt = 1; single[0] = ap; terms = single; }

    Expr* total = mk_int(0);
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        Expr* den = ev1("Denominator", expr_copy(T));
        if (!contains_symbol(den, v)) {                 /* polynomial part */
            expr_free(den);
            if (!is_zero_q(T)) { ok = false; break; }   /* divergent */
            continue;                                    /* zero -> skip */
        }
        Expr* num = ev1("Numerator", expr_copy(T));
        /* den = base^m (m>=1). */
        Expr* base; long m;
        if (head_name_is(den, "Power") && den->data.function.arg_count == 2 &&
            den->data.function.args[1]->type == EXPR_INTEGER &&
            den->data.function.args[1]->data.integer >= 1) {
            base = expr_copy(den->data.function.args[0]);
            m = (long)den->data.function.args[1]->data.integer;
        } else { base = expr_copy(den); m = 1; }
        expr_free(den);
        /* base must be linear in v (alpha v + beta), num free of v. */
        Expr* alpha = ev1("Simplify", mk_fn3("Coefficient", expr_copy(base), expr_copy(v), mk_int(1)));
        Expr* beta  = ev1("Simplify", mk_fn3("Coefficient", expr_copy(base), expr_copy(v), mk_int(0)));
        Expr* quad  = ev1("Simplify", mk_fn3("Coefficient", base, expr_copy(v), mk_int(2)));
        bool lin = quad && is_zero_q(quad);
        if (quad) expr_free(quad);
        Expr* d = (alpha && beta) ? ev1("Simplify", mk_fn2("Times", expr_copy(beta),
                        mk_fn2("Power", expr_copy(alpha), mk_int(-1)))) : NULL;
        if (!lin || !alpha || !d || contains_symbol(num, v) ||
            !real_positive(d, pb, np)) {
            if (num) expr_free(num); if (alpha) expr_free(alpha);
            if (beta) expr_free(beta); if (d) expr_free(d);
            ok = false; break;
        }
        expr_free(beta);
        /* ti = (1/2) num alpha^{-m} Pi Binomial[2m-2,m-1] 4^{-(m-1)} d^{(1-2m)/2} */
        Expr* ti = mk_fn2("Times",
            mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)), num),
            mk_fn2("Times",
                mk_fn2("Times", mk_fn2("Power", alpha, mk_int(-m)), mk_sym("Pi")),
                mk_fn2("Times",
                    mk_fn2("Times", mk_fn2("Binomial", mk_int(2*m-2), mk_int(m-1)),
                                    mk_fn2("Power", mk_int(4), mk_int(-(m-1)))),
                    mk_fn2("Power", d, mk_fn2("Rational", mk_int(1-2*m), mk_int(2))))));
        total = mk_fn2("Plus", total, ti);
    }
    expr_free(ap); expr_free(v);
    if (!ok) { expr_free(total); return NULL; }
    Expr* res = simplify_with(total, assumptions);
    expr_free(total);
    if (res && is_finite_value(res) && !contains_symbol(res, x)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* Simplify[Coefficient[e, v, k]]. */
static Expr* coeff_of(const Expr* e, const Expr* v, long k) {
    return ev1("Simplify", mk_fn3("Coefficient", expr_copy((Expr*)e),
                                  expr_copy((Expr*)v), mk_int(k)));
}

/* General real-rational half-line: ∫₀^∞ R(s) ds for a REAL rational R with no
 * poles on [0,∞) and a degree drop >= 2, WITHOUT the evenness restriction of
 * rational_halfline.  Partial-fractions R over s; each simple real factor
 * (s+r)^m (r>0) and irreducible quadratic factor ((s+β0)^2 + γ^2) (m = 1)
 * integrates to a rational / Log / ArcTan boundary value.  The Log(s->∞) pieces
 * must cancel (their coefficients sum to zero, guaranteed by the degree drop);
 * otherwise the integral diverges and we return NULL.  Produces REAL ArcTan/Log
 * output directly (no complex-Log reduction needed) -- this is what unlocks the
 * DECAYING sinc ∫₀^∞ e^{-p x} Sin[q x]/x dx = ArcTan[q/p], whose Laplace image
 * M(s) = q/((s+p)^2+q^2) is a non-even rational the even family declines. */
static Expr* rational_halfline_general(const Expr* g, const Expr* s,
                                       const ParamBound* pb, size_t np,
                                       const Expr* assumptions) {
    Expr* ap = ev2("Apart", expr_copy((Expr*)g), expr_copy((Expr*)s));
    if (!ap) return NULL;

    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(ap, "Plus")) { nt = ap->data.function.arg_count; terms = ap->data.function.args; }
    else { nt = 1; single[0] = ap; terms = single; }

    Expr* total    = mk_int(0);   /* finite boundary value accumulator */
    Expr* logcoeff = mk_int(0);   /* sum of Log(s->∞) coefficients; must vanish */
    bool ok = true;

    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        Expr* den = ev1("Denominator", expr_copy(T));
        if (!contains_symbol(den, s)) {            /* polynomial / constant part */
            expr_free(den);
            if (!is_zero_q(T)) ok = false;         /* nonzero const -> divergent */
            continue;
        }
        Expr* num = ev1("Numerator", expr_copy(T));
        /* den = base^m (base linear or irreducible quadratic in s). */
        Expr* base; long m;
        if (head_name_is(den, "Power") && den->data.function.arg_count == 2 &&
            den->data.function.args[1]->type == EXPR_INTEGER &&
            den->data.function.args[1]->data.integer >= 1) {
            base = ev1("Expand", expr_copy(den->data.function.args[0]));
            m = (long)den->data.function.args[1]->data.integer;
        } else { base = ev1("Expand", expr_copy(den)); m = 1; }
        expr_free(den);

        Expr* a2  = coeff_of(base, s, 2);
        Expr* a1  = coeff_of(base, s, 1);
        Expr* a0  = coeff_of(base, s, 0);
        Expr* B   = coeff_of(num,  s, 1);
        Expr* C   = coeff_of(num,  s, 0);
        Expr* nHi = coeff_of(num,  s, 2);
        expr_free(base); expr_free(num);
        bool numlin = nHi && is_zero_q(nHi) && B && C &&
                      !contains_symbol(B, s) && !contains_symbol(C, s);
        if (nHi) expr_free(nHi);
        bool quad = a2 && !is_zero_q(a2);
        if (!numlin || !a1 || !a0) {
            if (a2) expr_free(a2); if (a1) expr_free(a1); if (a0) expr_free(a0);
            if (B) expr_free(B); if (C) expr_free(C); ok = false; break;
        }

        if (quad) {
            /* base = a2[(s+β0)^2 + γ2], β0 = a1/(2 a2), γ2 = a0/a2 - β0^2 > 0. */
            if (m != 1) { expr_free(a2);expr_free(a1);expr_free(a0);expr_free(B);expr_free(C); ok=false; break; }
            Expr* beta0 = ev1("Simplify", t_mul(a1, t_pow(t_mul(mk_int(2), expr_copy(a2)), -1)));
            Expr* g2    = ev1("Simplify", t_add(t_mul(a0, t_pow(expr_copy(a2), -1)),
                                                t_neg(t_pow(expr_copy(beta0), 2))));
            if (!beta0 || !g2 || !real_positive(g2, pb, np)) {
                if (beta0) expr_free(beta0); if (g2) expr_free(g2);
                expr_free(a2); expr_free(B); expr_free(C); ok = false; break;
            }
            Expr* gam = ev1("Simplify", mk_fn2("Power", expr_copy(g2), t_rat(1, 2)));
            /* logcoeff += B/a2  (Log((s+β0)^2+γ^2) ~ 2 Log s at ∞). */
            logcoeff = t_add(logcoeff, t_mul(expr_copy(B), t_pow(expr_copy(a2), -1)));
            /* finite log at 0: -(B/(2 a2)) Log(β0^2 + γ^2). */
            total = t_add(total,
                t_mul(t_neg(t_mul(expr_copy(B), t_pow(t_mul(mk_int(2), expr_copy(a2)), -1))),
                      mk_fn1("Log", t_add(t_pow(expr_copy(beta0), 2), expr_copy(g2)))));
            /* ArcTan part: coef * [ArcTan((s+β0)/γ)]_0^∞ = coef*(π/2 - ArcTan(β0/γ)),
             * coef = (C - B β0)/(a2 γ).  For β0 > 0 present the boundary value in
             * the reciprocal-free form ArcTan(γ/β0) (= π/2 - ArcTan(β0/γ)); this
             * keeps a decaying-sinc result as ArcTan[q/p] rather than
             * π/2 - ArcTan[p/q], which the parameter back-integration cannot
             * integrate (the engine chokes on ArcTan of a reciprocal argument). */
            Expr* coef = t_mul(t_add(expr_copy(C), t_neg(t_mul(expr_copy(B), expr_copy(beta0)))),
                               t_pow(t_mul(expr_copy(a2), expr_copy(gam)), -1));
            Expr* atval;
            if (real_positive(beta0, pb, np))
                atval = mk_fn1("ArcTan", t_mul(expr_copy(gam), t_pow(expr_copy(beta0), -1)));
            else
                atval = t_add(t_mul(t_rat(1, 2), mk_sym("Pi")),
                              t_neg(mk_fn1("ArcTan", t_mul(expr_copy(beta0), t_pow(expr_copy(gam), -1)))));
            total = t_add(total, t_mul(coef, atval));
            expr_free(a2); expr_free(B); expr_free(C);
            expr_free(beta0); expr_free(g2); expr_free(gam);
        } else {
            /* Simple real factor: (s+r)^m, r = a0/a1 > 0 (pole off [0,∞)); after
             * Apart a linear-power term has a CONSTANT numerator (B == 0). */
            expr_free(a2);
            Expr* r = ev1("Simplify", t_mul(expr_copy(a0), t_pow(expr_copy(a1), -1)));
            if (!is_zero_q(B) || !r || !real_positive(r, pb, np)) {
                expr_free(a1); expr_free(a0); expr_free(B); expr_free(C);
                if (r) expr_free(r); ok = false; break;
            }
            if (m == 1) {
                Expr* co = ev1("Simplify", t_mul(expr_copy(C), t_pow(expr_copy(a1), -1))); /* A/a1 */
                logcoeff = t_add(logcoeff, expr_copy(co));
                total = t_add(total, t_neg(t_mul(co, mk_fn1("Log", expr_copy(r)))));
            } else {
                /* (A/a1^m) r^{1-m}/(m-1). */
                total = t_add(total,
                    t_mul(t_mul(expr_copy(C), t_pow(expr_copy(a1), -m)),
                          t_mul(t_pow(expr_copy(r), 1 - m), t_rat(1, m - 1))));
            }
            expr_free(a1); expr_free(a0); expr_free(B); expr_free(C); expr_free(r);
        }
    }
    expr_free(ap);
    if (ok && !is_zero_with(logcoeff, assumptions)) ok = false;   /* divergent */
    expr_free(logcoeff);
    if (!ok) { expr_free(total); return NULL; }
    Expr* res = simplify_with(total, assumptions);
    expr_free(total);
    if (res && is_finite_value(res) && !contains_symbol(res, s)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* Gaussian moment half-line: ∫₀^∞ c x^n e^{-p x^2} {1 | Cos[q x]} dx (p > 0),
 *   no trig:      (c/2) Γ((n+1)/2) p^{-(n+1)/2}   [n = 0 gives (c/2) Sqrt(π/p)],
 *   cosine, n=0:  (c/2) Sqrt(π/p) e^{-q^2/(4 p)}.
 * The Sin moment is a Dawson/Erfi form and is deliberately declined: the
 * DiffUnderInt targets differentiate the Sin away (∫₀^∞ e^{-x^2} Sin[a x]/x dx
 * differentiates to the cosine moment).  Returns NULL for any other form. */
static Expr* gaussian_halfline(const Expr* g, const Expr* x,
                               const ParamBound* pb, size_t np,
                               const Expr* assumptions) {
    if (!contains_gaussian_exp(g, x)) return NULL;    /* fast reject non-Gaussian */
    Expr* gg = ev1("Expand", expr_copy((Expr*)g));
    if (!gg) return NULL;
    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(gg, "Plus")) { nt = gg->data.function.arg_count; terms = gg->data.function.args; }
    else { nt = 1; single[0] = gg; terms = single; }

    Expr* total = mk_int(0);
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        Expr* facbuf[64];
        size_t nf = collect_factors(T, facbuf, 64, 0);
        Expr* exparg = NULL;             /* the -p x^2 (+off) exponent */
        Expr* qcos   = NULL;             /* q x from a single Cos[q x] factor */
        bool sawsin = false, multi_exp = false, multi_trig = false;
        Expr* rest = mk_int(1);
        for (size_t i = 0; i < nf; i++) {
            Expr* F = facbuf[i];
            const Expr* ea = NULL;
            if (head_name_is(F, "Exp") && F->data.function.arg_count == 1) ea = F->data.function.args[0];
            else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                     F->data.function.args[0]->type == EXPR_SYMBOL &&
                     strcmp(F->data.function.args[0]->data.symbol, "E") == 0) ea = F->data.function.args[1];
            if (ea) { if (exparg) multi_exp = true; else exparg = expr_copy((Expr*)ea); continue; }
            if ((head_name_is(F, "Cos") || head_name_is(F, "Sin")) &&
                F->data.function.arg_count == 1 && contains_symbol(F->data.function.args[0], x)) {
                if (qcos || sawsin) multi_trig = true;
                if (head_name_is(F, "Sin")) sawsin = true;
                else qcos = expr_copy(F->data.function.args[0]);
                continue;
            }
            rest = t_mul(rest, expr_copy(F));
        }
        if (!exparg || multi_exp || multi_trig || sawsin) {
            if (exparg) expr_free(exparg); if (qcos) expr_free(qcos); expr_free(rest);
            ok = false; break;
        }
        /* exponent = -p x^2 + off, no linear term. */
        Expr* p   = ev1("Simplify", t_neg(mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)x), mk_int(2))));
        Expr* lin = ev1("Simplify", mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)x), mk_int(1)));
        Expr* off = ev1("Simplify", mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)x), mk_int(0)));
        expr_free(exparg);
        bool pform = lin && is_zero_q(lin) && p && real_positive(p, pb, np);
        if (lin) expr_free(lin);
        /* rest = c x^n. */
        Expr* restc = ev1("Simplify", rest);
        long n = x_monomial_degree(restc, x);
        Expr* one = mk_int(1);
        Expr* c = (n >= 0) ? subst(restc, x, one) : NULL;
        expr_free(one); expr_free(restc);
        if (!pform || n < 0 || !c || contains_symbol(c, x)) {
            if (c) expr_free(c); if (p) expr_free(p); if (off) expr_free(off);
            if (qcos) expr_free(qcos); ok = false; break;
        }
        Expr* ce = t_mul(c, mk_fn1("Exp", off));   /* c e^{off} */
        Expr* ti = NULL;
        if (qcos) {
            /* q x must be pure linear: q = Coefficient[.,x,1], no constant part. */
            Expr* q  = ev1("Simplify", mk_fn3("Coefficient", expr_copy(qcos), expr_copy((Expr*)x), mk_int(1)));
            Expr* q0 = ev1("Simplify", mk_fn3("Coefficient", qcos, expr_copy((Expr*)x), mk_int(0)));
            bool qok = n == 0 && q && q0 && is_zero_q(q0) && !contains_symbol(q, x);
            if (q0) expr_free(q0);
            if (!qok) { if (q) expr_free(q); expr_free(ce); expr_free(p); ok = false; break; }
            /* (c e^{off}/2) Sqrt(π/p) e^{-q^2/(4 p)} */
            ti = t_mul(t_mul(t_rat(1, 2), ce),
                   t_mul(mk_fn2("Power", t_mul(mk_sym("Pi"), t_pow(expr_copy(p), -1)), t_rat(1, 2)),
                         mk_fn1("Exp", t_mul(t_rat(-1, 4), t_mul(t_pow(q, 2), t_pow(expr_copy(p), -1))))));
            expr_free(p);
        } else {
            /* (c e^{off}/2) Γ((n+1)/2) p^{-(n+1)/2} */
            ti = t_mul(t_mul(t_rat(1, 2), ce),
                   t_mul(mk_fn1("Gamma", t_rat(n + 1, 2)),
                         mk_fn2("Power", p, t_rat(-(n + 1), 2))));
        }
        total = t_add(total, ti);
    }
    expr_free(gg);
    if (!ok) { expr_free(total); return NULL; }
    Expr* res = simplify_with(total, assumptions);
    expr_free(total);
    if (res && is_finite_value(res) && !contains_symbol(res, x)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* Gaussian parameter back-integration: ∫ c e^{-k p^2} dp = c (1/2) Sqrt(π/k) Erf(Sqrt(k) p),
 * for a J that is a sum of pure Gaussian-in-p terms (k a positive real constant).
 * The engine does not produce this Erf antiderivative, so DiffUnderInt supplies
 * it directly; NULL for anything non-Gaussian (the engine handles those). */
static Expr* integrate_gaussian_param(const Expr* J, const Expr* p) {
    if (!contains_gaussian_exp(J, p)) return NULL;
    Expr* jj = ev1("Expand", expr_copy((Expr*)J));
    if (!jj) return NULL;
    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(jj, "Plus")) { nt = jj->data.function.arg_count; terms = jj->data.function.args; }
    else { nt = 1; single[0] = jj; terms = single; }

    Expr* total = mk_int(0);
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* T = terms[t];
        Expr* facbuf[64];
        size_t nf = collect_factors(T, facbuf, 64, 0);
        Expr* exparg = NULL; bool multi_exp = false;
        Expr* rest = mk_int(1);
        for (size_t i = 0; i < nf; i++) {
            Expr* F = facbuf[i];
            const Expr* ea = NULL;
            if (head_name_is(F, "Exp") && F->data.function.arg_count == 1) ea = F->data.function.args[0];
            else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                     F->data.function.args[0]->type == EXPR_SYMBOL &&
                     strcmp(F->data.function.args[0]->data.symbol, "E") == 0) ea = F->data.function.args[1];
            if (ea) { if (exparg) multi_exp = true; else exparg = expr_copy((Expr*)ea); continue; }
            rest = t_mul(rest, expr_copy(F));
        }
        if (!exparg || multi_exp) { if (exparg) expr_free(exparg); expr_free(rest); ok = false; break; }
        /* exponent = -k p^2 + off, no linear term, k a positive constant. */
        Expr* k   = ev1("Simplify", t_neg(mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)p), mk_int(2))));
        Expr* lin = ev1("Simplify", mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)p), mk_int(1)));
        Expr* off = ev1("Simplify", mk_fn3("Coefficient", expr_copy(exparg), expr_copy((Expr*)p), mk_int(0)));
        expr_free(exparg);
        double kv;
        bool kok = lin && is_zero_q(lin) && k && numeric_double(k, &kv) && kv > 1e-12;
        if (lin) expr_free(lin);
        /* rest must be a p-free constant c (poly x Gaussian antiderivatives are
         * not elementary in general). */
        Expr* c = ev1("Simplify", rest);
        if (!kok || !c || contains_symbol(c, p)) {
            if (c) expr_free(c); if (k) expr_free(k); if (off) expr_free(off);
            ok = false; break;
        }
        Expr* ce = t_mul(c, mk_fn1("Exp", off));
        /* c e^{off} (1/2) Sqrt(π/k) Erf(Sqrt(k) p) */
        Expr* ti = t_mul(t_mul(t_rat(1, 2), ce),
                     t_mul(mk_fn2("Power", t_mul(mk_sym("Pi"), t_pow(expr_copy(k), -1)), t_rat(1, 2)),
                           mk_fn1("Erf", t_mul(mk_fn2("Power", k, t_rat(1, 2)), expr_copy((Expr*)p)))));
        total = t_add(total, ti);
    }
    expr_free(jj);
    if (!ok) { expr_free(total); return NULL; }
    Expr* res = ev1("Simplify", total);
    if (res && is_finite_value(res)) return res;
    if (res) expr_free(res);
    return NULL;
}

/* Region-aware inner definite integral.
 *
 * The general integrator HANGS (uninterruptibly from inside this builtin --
 * TimeConstrained does NOT bound a nested evaluate) on the half-line
 * exponential/trig forms Feynman's trick produces, so we never hand those to it.
 *   - half-line {0, +Inf}: only the fast closed-form families are tried; if none
 *     matches, return NULL (that parameter/base is abandoned, no hang);
 *   - finite region: the engine is used only for provably-safe integrands (no
 *     trig-of-x, no radical-of-x, no Gaussian), which it evaluates quickly
 *     (e.g. Integrate[x^a, {x,0,1}]); other finite forms need a family. */
static Expr* inner_definite(const Expr* g, const Expr* x, const Expr* a,
                            const Expr* b, const Expr* assumptions,
                            const ParamBound* pb, size_t np) {
    if (is_zero_expr(a) && is_pos_inf(b)) {
        Expr* r = laplace_halfline(g, x, pb, np, assumptions);
        if (!r) r = laplace_sinc_halfline(g, x, pb, np, assumptions);
        if (!r) r = rational_halfline(g, x, pb, np, assumptions);
        if (!r) r = rational_halfline_general(g, x, pb, np, assumptions);
        if (!r) r = gaussian_halfline(g, x, pb, np, assumptions);
#ifdef DIUI_DEBUG
        fprintf(stderr, "DIUI:   [%.0fms] half-line family -> %s\n",
                diui_ms(), r ? "HIT" : "miss");
#endif
        return r;                                /* families only -- no engine */
    }
    if (has_trig_of_x(g, x) || has_radical_of_x(g, x) || contains_gaussian_exp(g, x))
        return NULL;                             /* would hang the engine */
    return integrate_definite_of(g, x, a, b, assumptions);
}

/* -------------------------------------------------------------------------
 * Base-point search: find (p0, I0) with an EXACT known integral value I(p0).
 * Ordered by trustworthiness: zero-integrand bases first (I0 = 0 exactly), then
 * a directly-integrable base computed by the engine.  On success returns owned
 * *out_p0, *out_i0.
 * ---------------------------------------------------------------------- */
static bool find_base(const Expr* f, const Expr* x, const Expr* a, const Expr* b,
                      const Expr* assumptions, const Expr* p,
                      const ParamBound* pb, size_t np,
                      Expr** out_p0, Expr** out_i0) {
    /* 1. other-parameter zero base: f|_{p->q} identically 0 in x. */
    for (size_t i = 0; i < np; i++) {
        if (pb[i].sym == p->data.symbol) continue;
        Expr* q = mk_sym(pb[i].sym);
        Expr* fq = subst(f, p, q);
        bool zero = fq && is_zero_q(fq);
        if (fq) expr_free(fq);
        if (zero) { *out_p0 = q; *out_i0 = mk_int(0); return true; }
        expr_free(q);
    }
    /* 2. numeric zero base: f|_{p->c} identically 0 in x, c in {0, 1, -1}. */
    static const long CS[] = { 0, 1, -1 };
    for (size_t i = 0; i < sizeof(CS)/sizeof(CS[0]); i++) {
        Expr* c = mk_int(CS[i]);
        Expr* fc = subst(f, p, c);
        bool zero = fc && is_zero_q(fc);
        if (fc) expr_free(fc);
        if (zero) { *out_p0 = c; *out_i0 = mk_int(0); return true; }
        expr_free(c);
    }
    /* 3. directly-integrable base: I(p0) = Integrate[f|_{p->p0}, {x,a,b}] closes.
     *    Only the cheap numeric anchor p0 = 0 is tried: substituting one
     *    parameter for another tends to CREATE a harder (or engine-hanging)
     *    integral, so it is deliberately excluded here (the zero-integrand
     *    other-parameter base above already covers the q-collapse cases). */
    Expr* cands[2]; size_t nc = 0;
    cands[nc++] = mk_int(0);
    bool found = false;
    for (size_t i = 0; i < nc && !found; i++) {
        Expr* fp0 = subst(f, p, cands[i]);
        Expr* I0  = fp0 ? inner_definite(fp0, x, a, b, assumptions, pb, np) : NULL;
        if (fp0) expr_free(fp0);
        if (I0) { Expr* s = simplify_with(I0, assumptions); expr_free(I0); I0 = s; }
        if (I0 && !contains_head(I0, "Integrate") && is_finite_value(I0)) {
            *out_p0 = expr_copy(cands[i]); *out_i0 = I0; found = true;
        } else if (I0) {
            expr_free(I0);
        }
    }
    for (size_t i = 0; i < nc; i++) expr_free(cands[i]);
    return found;
}

/* Output cleanup.  Simplify canonicalises c*Log[w] into the contracted
 * Log[w^c] form (e.g. -(1/2)Log[u] + (1/2)Log[v] -> Log[1/Sqrt[u]] + Log[Sqrt[v]]),
 * which reads poorly for the Frullani/Laplace family results.  1-arg PowerExpand
 * pulls those powers back out; we keep the expanded form ONLY when it is provably
 * equal to the verified answer (Simplify[clean - I] === 0 under the assumptions),
 * so a branch-changing PowerExpand can never corrupt a correct result.  Borrows
 * `I`; returns an owned cleaned copy, or NULL to keep `I` as-is. */
static Expr* diui_finalize(const Expr* I, const Expr* assumptions) {
    Expr* clean = ev1("PowerExpand", expr_copy((Expr*)I));
    if (!clean) return NULL;
    Expr* diff = t_add(expr_copy(clean), t_neg(expr_copy((Expr*)I)));
    bool same = is_zero_with(diff, assumptions);
    expr_free(diff);
    if (same && is_finite_value(clean)) return clean;
    expr_free(clean);
    return NULL;
}

/* -------------------------------------------------------------------------
 * Stage A -- pure quadrature (lambda = 0).
 * ---------------------------------------------------------------------- */
static Expr* stage_quadrature(const Expr* f, const Expr* x, const Expr* a,
                              const Expr* b, const Expr* assumptions,
                              const ParamBound* pb, size_t np) {
    for (size_t pi = 0; pi < np; pi++) {
        Expr* p = mk_sym(pb[pi].sym);
#ifdef DIUI_DEBUG
        fprintf(stderr, "DIUI: [%.0fms] try param %s\n", diui_ms(), pb[pi].sym);
#endif

        /* g = Simplify[D[f, p]].  Skip if f is independent of p. */
        Expr* g = ev1("Simplify", deriv(f, p));
        if (!g || is_zero_q(g)) {
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI:   g is zero / null -> skip\n");
#endif
            if (g) expr_free(g); expr_free(p); continue; }

        /* J(p) = Integrate[g, {x,a,b}].  Must close to a finite closed form.
         * Simplify collapses the FTC boundary artifacts (e.g. the lower-limit
         * `0^(1+a)` term in Integrate[x^a,{x,0,1}]) that would otherwise stop
         * the parameter-integration from closing. */
        Expr* J = inner_definite(g, x, a, b, assumptions, pb, np);
        expr_free(g);
        if (J) { Expr* Js = simplify_with(J, assumptions); expr_free(J); J = Js; }
        if (!J || contains_head(J, "Integrate") || !is_finite_value(J)) {
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI:   [%.0fms] J did not close / non-finite -> skip\n", diui_ms());
#endif
            if (J) expr_free(J); expr_free(p); continue;
        }

        /* G(p) = Integrate[J, p]  (antiderivative over the parameter). */
        Expr* G = integrate_over_param(J, p);
        if (!G || contains_head(G, "Integrate") || !is_finite_value(G)) {
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI:   [%.0fms] G did not close -> skip\n", diui_ms());
#endif
            if (G) expr_free(G); expr_free(J); expr_free(p); continue;
        }

        /* Base point with an exact known value I(p0). */
        Expr* p0 = NULL; Expr* I0 = NULL;
        if (!find_base(f, x, a, b, assumptions, p, pb, np, &p0, &I0)) {
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI:   [%.0fms] no base found -> skip\n", diui_ms());
#endif
            expr_free(G); expr_free(J); expr_free(p); continue;
        }

        /* I(p) = G(p) - G(p0) + I(p0). */
        Expr* Gp0 = eval_at_param(G, p, p0);
        if (!Gp0) {
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI:   G(p0) not finite -> skip\n");
#endif
            expr_free(p0); expr_free(I0); expr_free(G); expr_free(J); expr_free(p); continue; }

        Expr* sum = mk_fn2("Plus", expr_copy(G),
                           mk_fn2("Plus",
                                  mk_fn2("Times", mk_int(-1), Gp0),
                                  expr_copy(I0)));
        Expr* I = simplify_with(sum, assumptions);
        expr_free(sum);
        expr_free(p0); expr_free(I0); expr_free(G);

        /* Symbolic verification: D[I,p] - J === 0 (under the assumptions). */
        bool ok = false;
        if (I && is_finite_value(I)) {
            Expr* chk = mk_fn2("Plus", deriv(I, p),
                               mk_fn2("Times", mk_int(-1), expr_copy(J)));
            ok = is_zero_with(chk, assumptions);
            expr_free(chk);
        }
#ifdef DIUI_DEBUG
        fprintf(stderr, "DIUI:   [%.0fms] verify D[I,p]-J==0 : %s\n", diui_ms(), ok ? "PASS" : "FAIL");
#endif
        expr_free(J); expr_free(p);
        if (ok) {
            Expr* cl = diui_finalize(I, assumptions);
            if (cl) { expr_free(I); return cl; }
            return I;
        }
        if (I) expr_free(I);
    }
    return NULL;
}

/* -------------------------------------------------------------------------
 * Public entry point.
 * ---------------------------------------------------------------------- */
/* Whole-line divergence gate.  Integrate[f, {x,-Inf,Inf}] with a real pole of f
 * on the axis diverges (only a principal value exists), so DiffUnderInt cannot
 * help -- and pursuing it drives a non-terminating escalation (Integrate[x^k f],
 * plus Exp->Cosh/Sinh rewrites that swell the Risch-Norman linear system to
 * minutes).  Detect a real pole by a numeric sign-change / zero scan of the
 * x-only denominator: Solve is unreliable on transcendental denominators (for
 * Exp[x]-1 vs Exp[x]+1 it returns the same complex family under Reals).  Gated to
 * a parameter-free denominator; conservative (returns false) otherwise, and the
 * numerator being nonzero at an Exp-tower root is implicit (a genuine pole). */
static bool whole_line_divergent_pole(const Expr* f, const Expr* x) {
    Expr* T = ev1("Together", expr_copy((Expr*)f));
    if (!T) return false;
    Expr* den = ev1("Denominator", T);                /* consumes T */
    if (!den) return false;
    if (!contains_symbol(den, x)) { expr_free(den); return false; }
    ParamBound pb[4];
    if (collect_params(den, x, pb, 4, 0) > 0) { expr_free(den); return false; }
    bool pole = false, have_prev = false; double prev = 0.0;
    for (int k = -40; k <= 40 && !pole; k++) {
        Expr* at = eval_take(mk_fn2("ReplaceAll", expr_copy(den),
                       mk_fn2("Rule", expr_copy((Expr*)x), mk_int(k))));
        double v;
        bool got = at && numeric_double(at, &v);
        if (at) expr_free(at);
        if (!got) { have_prev = false; continue; }
        if (v == 0.0) { pole = true; break; }
        if (have_prev && ((prev < 0.0) != (v < 0.0))) pole = true;   /* sign change */
        prev = v; have_prev = true;
    }
    expr_free(den);
    return pole;
}

Expr* integrate_diffunderint_try(Expr* f, Expr* x, Expr* a, Expr* b,
                                 Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;          /* no x: not our business */
    if (diui_depth >= DIUI_MAX_DEPTH) return NULL;
    /* Divergent whole-line integrand (real axis pole): decline before the
     * escalation.  Only the pure two-sided-infinite case; half-line / finite
     * integrals keep their behaviour (a boundary singularity is not interior). */
    if (is_neg_inf(a) && is_pos_inf(b) && whole_line_divergent_pole(f, x))
        return NULL;
    /* Gaussian integrands (Exp[-p x^2], ...) are now handled: the half-line inner
     * integrals go only to the closed-form families (never the general engine,
     * which hangs on Exp[nonlinear-in-x]), gaussian_halfline supplies the moment
     * integrals, and integrate_gaussian_param supplies the Erf back-integration.
     * A finite-region Gaussian is still declined inside inner_definite (it would
     * hand the engine a hanging form), so such f simply finds no closing
     * parameter and returns unevaluated -- no hang. */

    /* Free parameters of the integrand (besides x). */
    ParamBound pb[16];
    size_t np = collect_params(f, x, pb, 16, 0);
    if (np == 0) return NULL;                         /* no parameter to vary */
    if (assumptions) absorb_fact(pb, np, assumptions);

    /* The method probes many divergent candidate sub-expressions (Limit at a
     * pole, 1/0 in a degenerate family M(s), ...).  Those arithmetic warnings
     * are spurious search noise -- mute them like Limit does, keeping the return
     * values (ComplexInfinity/Indeterminate) that gate the search unchanged. */
    diui_depth++;
    arith_warnings_mute_push();
    Expr* result = stage_quadrature(f, x, a, b, assumptions, pb, np);
    arith_warnings_mute_pop();
    diui_depth--;
    return result;
}

/* -------------------------------------------------------------------------
 * `Integrate`DiffUnderInt[f, {x,a,b}]` (optionally Assumptions -> ...) builtin.
 * ---------------------------------------------------------------------- */
Expr* builtin_integrate_diffunderint(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f    = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3)
        return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol == SYM_Rule ||
             opt->data.function.head->data.symbol == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                strcmp(lhs->data.symbol, "Assumptions") == 0) {
                assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;                                  /* unknown trailing arg */
    }
    return integrate_diffunderint_try(f, x, a, b, assumptions);
}

void integrate_diffunderint_init(void) {
    symtab_add_builtin("Integrate`DiffUnderInt", builtin_integrate_diffunderint);
    symtab_get_def("Integrate`DiffUnderInt")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`DiffUnderInt",
        "Integrate`DiffUnderInt[f, {x, a, b}] evaluates a parameter-dependent "
        "definite integral by differentiation under the integral sign (the "
        "Leibniz rule / Feynman's trick): it differentiates the integrand with "
        "respect to a free parameter, evaluates the resulting simpler definite "
        "integral, integrates back over the parameter, and fixes the constant "
        "from an exact base value.  Accepts Assumptions -> ... constraining the "
        "parameters.  Returns unevaluated when no parameter yields a closed "
        "form.");
}
