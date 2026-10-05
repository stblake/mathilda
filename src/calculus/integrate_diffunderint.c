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

/* True iff the bare symbol `name` occurs anywhere in `e`. */
static bool contains_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol_name(e->data.function.args[i], name)) return true;
    return false;
}

/* True iff the symbol `x` occurs anywhere in `e`. */
static bool contains_symbol(const Expr* e, const Expr* x) {
    return contains_symbol_name(e, x->data.symbol.name);
}

/* True iff `e` contains Power[x, p] with p a constant other than 0 or 1 -- a
 * genuinely nonlinear / singular power of x (x^2, x^(-2), Sqrt[x], ...). */
static bool has_nonlinear_x_power(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == x->data.symbol.name) {
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
        if (base->type == EXPR_SYMBOL && strcmp(base->data.symbol.name, "E") == 0)
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
        const char* h = e->data.function.head->data.symbol.name;
        for (size_t i = 0; i < sizeof(T)/sizeof(T[0]); i++)
            if (strcmp(h, T[i]) == 0 && contains_symbol(e, x)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_trig_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e`'s head is a forward trig / hyperbolic function. */
static bool is_trig_head_sym(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION ||
        e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    static const char* T[] = { "Sin","Cos","Tan","Cot","Sec","Csc",
                               "Sinh","Cosh","Tanh","Coth","Sech","Csch" };
    for (size_t i = 0; i < sizeof(T)/sizeof(T[0]); i++)
        if (strcmp(h, T[i]) == 0) return true;
    return false;
}

/* True iff some multiplicative term of `g` carries TWO OR MORE forward
 * trig/hyperbolic factors OF x (Sin[b x] Cos[c x], Sin[b x] Sin[c x], ...),
 * where a Power[trig, k] counts as a SINGLE factor.  Such a product drives the
 * Laplace/sinc complex-exponential path into a tangled Log/ArcTanh form that
 * carries a spurious imaginary part (or an Indeterminate) for some real
 * parameter signs; TrigReduce first rewrites the product as a sum of single
 * trig-of-x terms the families integrate cleanly.  A single trig factor raised
 * to a power (Sin[a x]^2, ...) is deliberately NOT a match -- those already
 * close through the existing paths and must stay untouched. */
static bool has_trig_product_of_x(const Expr* g, const Expr* x) {
    if (!g || g->type != EXPR_FUNCTION) return false;
    if (head_name_is(g, "Plus")) {
        for (size_t i = 0; i < g->data.function.arg_count; i++)
            if (has_trig_product_of_x(g->data.function.args[i], x)) return true;
        return false;
    }
    if (head_name_is(g, "Times")) {
        int c = 0;
        for (size_t i = 0; i < g->data.function.arg_count; i++) {
            Expr* f = g->data.function.args[i];
            Expr* base = f;
            if (head_name_is(f, "Power") && f->data.function.arg_count == 2)
                base = f->data.function.args[0];
            if (is_trig_head_sym(base) && contains_symbol(base, x)) c++;
        }
        return c >= 2;
    }
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

/* True iff `e` is a RATIONAL function of x: x appears only inside Plus/Times and
 * as the base of an INTEGER-exponent Power.  A transcendental kernel of x
 * (Exp[x] = Power[E, x], Log[x], Sin[x], ...) or an x-dependent / non-integer
 * exponent (a^x, x^(s-1), Sqrt[x]) is NOT rational.  The rational half-line
 * families below feed g to Apart[g, x], which is only meaningful -- and only
 * terminates cheaply -- for a rational g; a symbolic-exponent or exp-geometric
 * integrand (Mellin/Ramanujan territory) otherwise drives Apart into an
 * expensive non-terminating rewrite.  Conservative: an x-dependent head we do
 * not model returns false (decline), never a spurious accept. */
static bool is_rational_in_x(const Expr* e, const Expr* x) {
    if (!e || !contains_symbol(e, x)) return true;         /* x-free coefficient */
    if (e->type == EXPR_SYMBOL) return true;               /* the bare x itself */
    if (e->type != EXPR_FUNCTION) return true;
    if (head_name_is(e, "Plus") || head_name_is(e, "Times")) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!is_rational_in_x(e->data.function.args[i], x)) return false;
        return true;
    }
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (contains_symbol(ex, x)) return false;          /* a^x, x^x */
        if (ex->type != EXPR_INTEGER) return false;        /* x^(s-1), Sqrt[x] */
        return is_rational_in_x(base, x);
    }
    return false;                                          /* Exp/Log/Sin/... of x */
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

/* simplify_with that CONSUMES `e` (convenience for the finite-domain closers,
 * which build a throwaway tree and immediately simplify it). */
static Expr* simplify_take(Expr* e, const Expr* assumptions) {
    Expr* r = simplify_with(e, assumptions);
    expr_free(e);
    return r;
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
    /* Normalize before handing J to the engine: pull Logs out of radicals /
     * reciprocal powers (Log[1/u^(1/4)] -> -(1/4) Log[u]).  The indefinite engine
     * grinds UNINTERRUPTIBLY on a Log[1/poly^rational] back-integrand (it enters an
     * algebraic / radical integration path), yet the Log/ArcTan/rational spelling
     * PowerExpand produces integrates instantly.  PowerExpand can shift a branch,
     * but the caller re-verifies D[I,p]-J===0, so a branch slip is rejected, never
     * returned. */
    Expr* Jn = ev1("PowerExpand", expr_copy((Expr*)J));
    if (!Jn) Jn = expr_copy((Expr*)J);
    Expr* Je = ev1("Expand", Jn);                 /* consumes Jn */
    if (!Je) Je = expr_copy((Expr*)J);
    /* Bounded gate, mirroring inner_definite's families-only discipline: the
     * parameter back-integrand must be free of radical / trig / Gaussian OF p --
     * those drive the engine into the same uninterruptible grind the inner-integral
     * families avoid (TimeConstrained does not bound a nested evaluate).  Log /
     * ArcTan / rational / polynomial in p are elementary, fast and terminating. */
    if (has_radical_of_x(Je, p) || has_trig_of_x(Je, p) ||
        contains_gaussian_exp(Je, p)) {
        expr_free(Je);
        return NULL;                              /* abandon this parameter, no hang */
    }
    return ev2("Integrate", Je, expr_copy((Expr*)p));   /* consumes Je */
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
        if (e->data.symbol.name == x->data.symbol.name) return n;
        if (strcmp(e->data.symbol.name, "I") == 0) return n;
        for (size_t i = 0; i < n; i++)
            if (pb[i].sym == e->data.symbol.name) return n;
        double tmp;
        if (numeric_double(e, &tmp)) return n;       /* Pi, E, EulerGamma, ... */
        if (n < cap) { pb[n].sym = e->data.symbol.name; pb[n].lo = -HUGE_VAL;
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

/* Re[sym] / Im[sym] -> sym; any other expression is returned unchanged.  Lets an
 * assumption such as `Re[a] > 0` register a bound on the bare symbol `a` (the
 * bound is only ever used for SIGN gating, so the Re/Im distinction is immaterial
 * for the real-parameter families served here). */
static Expr* strip_re_im(Expr* e) {
    if ((head_name_is(e, "Re") || head_name_is(e, "Im")) &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_SYMBOL)
        return e->data.function.args[0];
    return e;
}

static void relation_apply(ParamBound* pb, size_t np, const char* op,
                           Expr* L, Expr* R) {
    bool lt = (strcmp(op, "Less") == 0 || strcmp(op, "LessEqual") == 0);
    bool gt = (strcmp(op, "Greater") == 0 || strcmp(op, "GreaterEqual") == 0);
    if (!lt && !gt) return;
    L = strip_re_im(L);
    R = strip_re_im(R);
    double c;
    if (L->type == EXPR_SYMBOL && numeric_double(R, &c))
        bound_apply(pb, np, L->data.symbol.name, c, lt ? -1 : +1);
    else if (R->type == EXPR_SYMBOL && numeric_double(L, &c))
        bound_apply(pb, np, R->data.symbol.name, c, lt ? +1 : -1);
}

/* Tighten parameter bounds from an Assumptions fact (And/List conjunctions,
 * chained Inequality, and the four ordered binary relations). */
static void absorb_fact(ParamBound* pb, size_t np, Expr* fact) {
    if (!fact || fact->type != EXPR_FUNCTION) return;
    if (fact->data.function.head->type != EXPR_SYMBOL) return;
    const char* h = fact->data.function.head->data.symbol.name;
    size_t ac = fact->data.function.arg_count;
    if (strcmp(h, "And") == 0 || strcmp(h, "List") == 0) {
        for (size_t i = 0; i < ac; i++) absorb_fact(pb, np, fact->data.function.args[i]);
        return;
    }
    if (strcmp(h, "Inequality") == 0 && ac >= 3) {
        for (size_t i = 0; i + 2 < ac; i += 2) {
            Expr* opsym = fact->data.function.args[i + 1];
            if (opsym->type != EXPR_SYMBOL) continue;
            relation_apply(pb, np, opsym->data.symbol.name,
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
    if (b->type == EXPR_SYMBOL && strcmp(b->data.symbol.name, "Infinity") == 0) return true;
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
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == x->data.symbol.name) return 1;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == x->data.symbol.name &&
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
                     strcmp(F->data.function.args[0]->data.symbol.name, "E") == 0)
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
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == x->data.symbol.name) return 1;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == x->data.symbol.name &&
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
                     strcmp(F->data.function.args[0]->data.symbol.name, "E") == 0) ea = F->data.function.args[1];
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
            if (c) { expr_free(c); } expr_free(alpha); expr_free(off); ok = false; break;
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
    if (!is_rational_in_x(g, x)) return NULL;   /* not our family (see is_rational_in_x) */
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
            if (num) { expr_free(num); } if (alpha) { expr_free(alpha); }
            if (beta) { expr_free(beta); } if (d) { expr_free(d); }
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
    if (!is_rational_in_x(g, s)) return NULL;   /* rational-in-s only; else Apart hangs */
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
            if (a2) { expr_free(a2); } if (a1) { expr_free(a1); } if (a0) { expr_free(a0); }
            if (B) { expr_free(B); } if (C) { expr_free(C); } ok = false; break;
        }

        if (quad) {
            /* base = a2[(s+β0)^2 + γ2], β0 = a1/(2 a2), γ2 = a0/a2 - β0^2 > 0. */
            if (m != 1) { expr_free(a2);expr_free(a1);expr_free(a0);expr_free(B);expr_free(C); ok=false; break; }
            Expr* beta0 = ev1("Simplify", t_mul(a1, t_pow(t_mul(mk_int(2), expr_copy(a2)), -1)));
            Expr* g2    = ev1("Simplify", t_add(t_mul(a0, t_pow(expr_copy(a2), -1)),
                                                t_neg(t_pow(expr_copy(beta0), 2))));
            if (!beta0 || !g2 || !real_positive(g2, pb, np)) {
                if (beta0) { expr_free(beta0); } if (g2) { expr_free(g2); }
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
                if (r) { expr_free(r); } ok = false; break;
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
                     strcmp(F->data.function.args[0]->data.symbol.name, "E") == 0) ea = F->data.function.args[1];
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
            if (exparg) { expr_free(exparg); } if (qcos) { expr_free(qcos); } expr_free(rest);
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
            if (c) { expr_free(c); } if (p) { expr_free(p); } if (off) { expr_free(off); }
            if (qcos) { expr_free(qcos); } ok = false; break;
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
                     strcmp(F->data.function.args[0]->data.symbol.name, "E") == 0) ea = F->data.function.args[1];
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
            if (c) { expr_free(c); } if (k) { expr_free(k); } if (off) { expr_free(off); }
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
        /* Collapse a PRODUCT of trig-of-x into a sum of single trig terms
         * (Sin[b x] Sin[c x] -> (Cos[(b-c)x] - Cos[(b+c)x])/2) before the
         * families run: the complex-exponential path otherwise yields a tangled
         * Log[b-c]/ArcTanh form that is non-real (or Indeterminate) for some real
         * parameter signs.  Left untouched when no such product is present. */
        Expr* gred = NULL;
        const Expr* gg = g;
        if (has_trig_product_of_x(g, x)) {
            gred = ev1("TrigReduce", expr_copy((Expr*)g));
            if (gred) gg = gred;
        }
        Expr* r = laplace_halfline(gg, x, pb, np, assumptions);
        if (!r) r = laplace_sinc_halfline(gg, x, pb, np, assumptions);
        if (!r) r = rational_halfline(gg, x, pb, np, assumptions);
        if (!r) r = rational_halfline_general(gg, x, pb, np, assumptions);
        if (!r) r = gaussian_halfline(gg, x, pb, np, assumptions);
        if (gred) expr_free(gred);
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
                      bool zero_base_only,
                      Expr** out_p0, Expr** out_i0) {
    /* 1. other-parameter zero base: f|_{p->q} identically 0 in x. */
    for (size_t i = 0; i < np; i++) {
        if (pb[i].sym == p->data.symbol.name) continue;
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
    /* A zero-integrand base is EXACT (I0 = 0, sign-independent).  The computed
     * base below is sign-sensitive (the family renders it for one sign branch)
     * and the D[I,p]-J verification cannot catch a wrong constant base, so a
     * caller that has another parameter to try prefers to skip it here. */
    if (zero_base_only) return false;
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

/* =========================================================================
 * Finite-domain Feynman families.
 *
 * The half-line families above all rely on a PRE-EXISTING parameter and an
 * engine-safe (rational/exp) inner integral.  Three important finite-domain
 * families fail both assumptions: they are numeric (no parameter to vary) or
 * their differentiated inner integral is trig/radical on a finite interval,
 * which the general engine cannot do (it hangs, or -- for the Form-B radical
 * integral -- returns a WRONG value).  So these closers (a) canonicalise the
 * integrand with a change of variables, (b) INTRODUCE an artificial Feynman
 * parameter, and (c) supply the inner integral in closed form themselves,
 * never routing it through the engine.
 *
 *   Family 1 (power-log):    Log[1 + c x^p] / (x Sqrt[1 - x^(2p)])  on {0,1}
 *                            -> (Pi^2/8 - ArcCos[c]^2/2) / p
 *   Family 2 (secant-radical): Sec[2x] Log[1 + c Sqrt[1 - Tan[x]^2]] on {0,Pi/4}
 *                            -> Pi^2/8 - ArcCos[c]^2/2
 *   Family 3 (tangent-power):  Csc[2x]^2 Log[1 + Tan[x]^a] on {0,Pi/4}
 *                            -> (Pi Csc[Pi/a] - a)/4
 *
 * Families 1 and 2 share the same inner integral ArcCos[q]/Sqrt[1-q^2] and the
 * same K(c) = Pi^2/8 - ArcCos[c]^2/2 shape; family 3 is a direct Beta/digamma
 * evaluation with a rational-anchor (a0 = 3) self-check.
 * ====================================================================== */

/* Exponent q with `term` = (x-free) * x^q for a single monomial in x; NULL if
 * `term` is not such a pure monomial.  Owned Integer/Rational/Real. */
static Expr* monomial_x_exponent(const Expr* term, const Expr* x) {
    if (!term) return NULL;
    if (term->type == EXPR_SYMBOL)
        return (term->data.symbol.name == x->data.symbol.name) ? mk_int(1) : NULL;
    if (head_name_is(term, "Power") && term->data.function.arg_count == 2) {
        Expr* base = term->data.function.args[0];
        Expr* ex   = term->data.function.args[1];
        if (base->type == EXPR_SYMBOL &&
            base->data.symbol.name == x->data.symbol.name && !contains_symbol(ex, x))
            return expr_copy(ex);
        return NULL;
    }
    if (head_name_is(term, "Times")) {
        Expr* found = NULL;
        for (size_t i = 0; i < term->data.function.arg_count; i++) {
            Expr* fac = term->data.function.args[i];
            if (!contains_symbol(fac, x)) continue;
            if (found) { expr_free(found); return NULL; }     /* two x-factors */
            found = monomial_x_exponent(fac, x);
            if (!found) return NULL;
        }
        return found;                                         /* NULL if none */
    }
    return NULL;
}

/* True if `e` contains Power[var, exp] with a NON-INTEGER exponent (fractional
 * constant like t^(5/2), or symbolic like t^a) of the bare integration variable.
 * The ArcCos family (power-log, secant-radical) never does -- its bare-variable
 * powers are all integers (u, u^-1, u^2), with fractional powers only of the
 * radicand (1-var^2)^(1/2) -- whereas the tangent-power family is exactly
 * var^(non-integer), and Simplify of a radical times var^(non-integer) HANGS.
 * stage_finite_feynman gates on this to stay off family-3 integrands. */
static bool has_noninteger_var_power(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL &&
            base->data.symbol.name == var->data.symbol.name &&
            ex->type != EXPR_INTEGER)
            return true;
    }
    if (has_noninteger_var_power(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_noninteger_var_power(e->data.function.args[i], var)) return true;
    return false;
}

/* First subexpression of the form Log[Plus[1, ...]] (borrowed), else NULL. */
static const Expr* find_log1plus(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (head_name_is(e, "Log") && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        if (head_name_is(arg, "Plus"))
            for (size_t i = 0; i < arg->data.function.arg_count; i++) {
                Expr* t = arg->data.function.args[i];
                if (t->type == EXPR_INTEGER && t->data.integer == 1) return e;
            }
    }
    const Expr* r = find_log1plus(e->data.function.head);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if ((r = find_log1plus(e->data.function.args[i]))) return r;
    return NULL;
}

/* x-exponent of the (single) non-1 summand of a Log[Plus[1, W]] node, requiring
 * W to be a pure monomial c*x^q.  Owned exponent or NULL. */
static Expr* log1plus_monomial_exponent(const Expr* lognode, const Expr* x) {
    Expr* arg = lognode->data.function.args[0];               /* Plus[...] */
    Expr* W = NULL; int nonone = 0;
    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* trm = arg->data.function.args[i];
        if (trm->type == EXPR_INTEGER && trm->data.integer == 1) continue;
        W = trm; nonone++;
    }
    if (nonone != 1) return NULL;
    return monomial_x_exponent(W, x);
}

/* Rewrite the single Log[Plus[1, W]] node to Log[Plus[1, Times[p, W]]].  Owned
 * copy; *count receives the number of such nodes rewritten. */
static Expr* feyn_rewrite(const Expr* e, const Expr* p, int* count) {
    if (!e) return NULL;
    if (head_name_is(e, "Log") && e->data.function.arg_count == 1) {
        Expr* arg = e->data.function.args[0];
        if (head_name_is(arg, "Plus")) {
            size_t ac = arg->data.function.arg_count; int one_idx = -1;
            for (size_t i = 0; i < ac; i++) {
                Expr* t = arg->data.function.args[i];
                if (t->type == EXPR_INTEGER && t->data.integer == 1) { one_idx = (int)i; break; }
            }
            if (one_idx >= 0) {
                (*count)++;
                size_t rem = ac - 1; Expr* W = NULL;
                if (rem == 1) {
                    for (size_t i = 0; i < ac; i++)
                        if ((int)i != one_idx) { W = expr_copy(arg->data.function.args[i]); break; }
                } else {
                    Expr** ws = malloc(rem * sizeof(Expr*)); size_t k = 0;
                    for (size_t i = 0; i < ac; i++)
                        if ((int)i != one_idx) ws[k++] = expr_copy(arg->data.function.args[i]);
                    W = expr_new_function(mk_sym("Plus"), ws, rem); free(ws);
                }
                return mk_fn1("Log", t_add(mk_int(1), t_mul(expr_copy((Expr*)p), W)));
            }
        }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head = feyn_rewrite(e->data.function.head, p, count);
    size_t n = e->data.function.arg_count;
    Expr** args = n ? malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) args[i] = feyn_rewrite(e->data.function.args[i], p, count);
    Expr* r = expr_new_function(head, args, n);
    if (args) free(args);
    return r;
}

/* Introduce the Feynman parameter: exactly one Log[1+W] -> Log[1+p W]. */
static Expr* insert_feynman_param(const Expr* canon, const Expr* p) {
    int count = 0;
    Expr* r = feyn_rewrite(canon, p, &count);
    if (count != 1) { if (r) expr_free(r); return NULL; }
    return r;
}

/* Simplify[e /. var -> val]; consumes `val`, borrows e/var. */
static Expr* at_point(const Expr* e, const Expr* var, Expr* val, const Expr* assum) {
    Expr* s = subst(e, var, val);
    expr_free(val);
    if (!s) return NULL;
    return simplify_take(s, assum);
}

/* Small algebraic builders. */
static Expr* one_plus_sq(const Expr* v)  { return t_add(mk_int(1), mk_fn2("Power", expr_copy((Expr*)v), mk_int(2))); }
static Expr* one_minus_sq(const Expr* v) { return t_add(mk_int(1), t_neg(mk_fn2("Power", expr_copy((Expr*)v), mk_int(2)))); }
/* Sqrt[1 - v^2]. */
static Expr* sqrt_1msq(const Expr* v)    { return mk_fn2("Power", one_minus_sq(v), t_rat(1, 2)); }
/* Times[2, x] (the "2x" argument of the trig-of-2x rewrite rules). */
static Expr* two_x_times(const Expr* x)  { return mk_fn2("Times", mk_int(2), expr_copy((Expr*)x)); }

/* A * ArcCos[q] * (1 - q^2)^(-1/2); consumes A and q. */
static Expr* mk_arccos_over_sqrt(Expr* A, Expr* q) {
    Expr* invs = mk_fn2("Power", t_add(mk_int(1), t_neg(mk_fn2("Power", expr_copy(q), mk_int(2)))),
                        t_rat(-1, 2));
    return t_mul(A, t_mul(mk_fn1("ArcCos", q), invs));
}

/* -(A/(2 alpha)) * ArcCos[q]^2; consumes A, alpha, q. */
static Expr* mk_G(Expr* A, Expr* alpha, Expr* q) {
    Expr* frac = t_mul(A, mk_fn2("Power", t_mul(mk_int(2), alpha), mk_int(-1)));
    Expr* ac2  = mk_fn2("Power", mk_fn1("ArcCos", q), mk_int(2));
    return t_neg(t_mul(frac, ac2));
}

/* (Pi Csc[Pi/e] - e)/4; borrows e. */
static Expr* emit_tanpow_form(const Expr* e) {
    Expr* csc = mk_fn1("Csc", t_mul(mk_sym("Pi"), mk_fn2("Power", expr_copy((Expr*)e), mk_int(-1))));
    Expr* inner = t_add(t_mul(mk_sym("Pi"), csc), t_neg(expr_copy((Expr*)e)));
    return t_mul(t_rat(1, 4), inner);
}

/* (1/4)(1 + t^2) t^(-2) Log[1 + t^expo]; borrows t, expo. */
static Expr* mk_tanpow_canon(const Expr* t, const Expr* expo) {
    Expr* logarg = t_add(mk_int(1), mk_fn2("Power", expr_copy((Expr*)t), expr_copy((Expr*)expo)));
    Expr* body = t_mul(one_plus_sq(t),
                       t_mul(mk_fn2("Power", expr_copy((Expr*)t), mk_int(-2)),
                             mk_fn1("Log", logarg)));
    return t_mul(t_rat(1, 4), body);
}

/* Family 1 normalizer: c(x) Log[1 + c x^p]/(x Sqrt[1 - x^(2p)]) on {0,1} with the
 * power substitution u = x^p, yielding (1/p) Log[1 + c u]/(u Sqrt[1-u^2]) on {0,1}.
 * On success returns the canonical integrand and *out_var = u; else NULL. */
static Expr* normalize_power_sub(const Expr* f, const Expr* x, const Expr* a,
                                 const Expr* b, Expr** out_var) {
    if (!(is_zero_expr(a) && b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;
    const Expr* lg = find_log1plus(f);
    if (!lg) return NULL;
    Expr* p = log1plus_monomial_exponent(lg, x);              /* the power p */
    if (!p) return NULL;
    /* p must be a positive real/rational constant. */
    {
        double pv;
        if (!numeric_double(p, &pv) || pv <= 0.0) { expr_free(p); return NULL; }
    }
    /* Residue check: f * x * Sqrt[1 - x^(2p)] / Log[1 + c x^p] must be x-free. */
    Expr* x2p  = mk_fn2("Power", expr_copy((Expr*)x), mk_fn2("Times", mk_int(2), expr_copy(p)));
    Expr* sqrt = mk_fn2("Power", t_add(mk_int(1), t_neg(x2p)), t_rat(1, 2));
    Expr* resid = simplify_take(
        t_mul(expr_copy((Expr*)f),
              t_mul(expr_copy((Expr*)x),
                    t_mul(sqrt, mk_fn2("Power", expr_copy((Expr*)lg), mk_int(-1))))),
        NULL);
    bool match = resid && !contains_symbol(resid, x) && is_finite_value(resid) && !is_zero_q(resid);
    if (resid) expr_free(resid);
    if (!match) { expr_free(p); return NULL; }
    /* Apply u = x^p : Integrate[f,{x,0,1}] = (1/p) Integrate[f(u^(1/p)) u^(1/p-1),{u,0,1}]. */
    Expr* u = mk_sym("$diuiSubU$");
    Expr* usub = mk_fn2("Power", expr_copy(u), mk_fn2("Power", expr_copy(p), mk_int(-1)));
    Expr* fsub = subst(f, x, usub); expr_free(usub);
    Expr* jac  = t_mul(mk_fn2("Power", expr_copy(p), mk_int(-1)),
                       mk_fn2("Power", expr_copy(u),
                              t_add(mk_fn2("Power", expr_copy(p), mk_int(-1)), mk_int(-1))));
    /* PowerExpand first: the substitution x = u^(1/p) creates nested powers like
     * (u^2)^(1/2) (for p = 1/2) that stay as Sqrt[u^2] (= Abs[u]) without a
     * positivity assumption; the integration variable u is on (0,1), so
     * PowerExpand's assume-positive is exactly correct and gives Log[1+u]. */
    Expr* canon = simplify_take(ev1("PowerExpand", t_mul(fsub, jac)), NULL);
    expr_free(p);
    if (!canon || contains_symbol(canon, x) || !contains_symbol(canon, u)) {
        if (canon) expr_free(canon);
        expr_free(u);
        return NULL;
    }
    *out_var = u;
    return canon;
}

/* Rational-trig-of-2x normalizer for {0, Pi/4}: t = Tan[x] (dx = dt/(1+t^2),
 * t : 0 -> 1).  Uses explicit rewrite rules (ReplaceAll[x -> ArcTan[t]] is not
 * reliable for Sec[2x]).  On success returns the canonical integrand in t on
 * {0,1} and *out_var = t; NULL if any x-dependence is left unresolved. */
static Expr* normalize_tan_half(const Expr* f, const Expr* x, const Expr* a,
                                const Expr* b, Expr** out_var) {
    if (!is_zero_expr(a)) return NULL;
    {   /* b == Pi/4 */
        Expr* d = t_add(expr_copy((Expr*)b), t_neg(t_mul(t_rat(1, 4), mk_sym("Pi"))));
        bool isq = is_zero_q(d); expr_free(d);
        if (!isq) return NULL;
    }
    Expr* T = mk_sym("$diuiTanT$");
    Expr* r_sec = mk_fn2("Rule", mk_fn1("Sec", two_x_times(x)),
                         t_mul(one_plus_sq(T), mk_fn2("Power", one_minus_sq(T), mk_int(-1))));
    Expr* r_csc = mk_fn2("Rule", mk_fn1("Csc", two_x_times(x)),
                         t_mul(one_plus_sq(T), mk_fn2("Power", t_mul(mk_int(2), expr_copy(T)), mk_int(-1))));
    Expr* r_cos = mk_fn2("Rule", mk_fn1("Cos", two_x_times(x)),
                         t_mul(one_minus_sq(T), mk_fn2("Power", one_plus_sq(T), mk_int(-1))));
    Expr* r_sin = mk_fn2("Rule", mk_fn1("Sin", two_x_times(x)),
                         t_mul(t_mul(mk_int(2), expr_copy(T)), mk_fn2("Power", one_plus_sq(T), mk_int(-1))));
    Expr* r_tan = mk_fn2("Rule", mk_fn1("Tan", expr_copy((Expr*)x)), expr_copy(T));
    Expr* rules[5] = { r_sec, r_csc, r_cos, r_sin, r_tan };
    Expr* rulelist = expr_new_function(mk_sym("List"), rules, 5);
    Expr* sub = ev2("ReplaceAll", expr_copy((Expr*)f), rulelist);
    if (!sub) { expr_free(T); return NULL; }
    Expr* canon = simplify_take(t_mul(sub, mk_fn2("Power", one_plus_sq(T), mk_int(-1))), NULL);
    if (!canon || contains_symbol(canon, x) || !contains_symbol(canon, T)) {
        if (canon) expr_free(canon);
        expr_free(T);
        return NULL;
    }
    *out_var = T;
    return canon;
}

/* Closed-form inner integrals on {0,1} of the ArcCos family, EMITTED (never
 * delegated: the engine cannot do Form A at a symbolic parameter and is
 * numerically WRONG on Form B).  Both evaluate to A * ArcCos[q]/Sqrt[1-q^2]:
 *   Form A:  A / ((1 + q u) Sqrt[1-u^2])
 *   Form B:  A / ((1 + q Sqrt[1-u^2]) Sqrt[1-u^2])
 * The candidate (A, q) is read off by point evaluation and then VERIFIED by
 * reconstructing g and Simplifying to 0, so a non-matching g is rejected.
 * On success sets owned *out_A, *out_q and returns true. */
static bool inner_arccos_family(const Expr* g, const Expr* var,
                                const Expr* assumptions, Expr** out_A, Expr** out_q) {
    Expr* s = simplify_take(t_mul(expr_copy((Expr*)g), sqrt_1msq(var)), assumptions);
    if (!s) return false;
    Expr* r = simplify_take(mk_fn2("Power", s, mk_int(-1)), assumptions);  /* r = 1/(g Sqrt) */
    if (!r) return false;
    Expr* r0 = at_point(r, var, mk_int(0), assumptions);
    Expr* r1 = at_point(r, var, mk_int(1), assumptions);
    expr_free(r);
    if (!r0 || !r1) { if (r0) expr_free(r0); if (r1) expr_free(r1); return false; }

    bool ok = false; Expr* A = NULL; Expr* q = NULL;

    /* Form A: r = (1 + q var)/A ; r|0 = 1/A, r|1 = (1+q)/A. */
    {
        Expr* Aa = simplify_take(mk_fn2("Power", expr_copy(r0), mk_int(-1)), assumptions);
        Expr* qa = simplify_take(t_mul(t_add(expr_copy(r1), t_neg(expr_copy(r0))),
                                       mk_fn2("Power", expr_copy(r0), mk_int(-1))), assumptions);
        if (Aa && qa && is_finite_value(Aa) && is_finite_value(qa) &&
            !contains_symbol(Aa, var) && !contains_symbol(qa, var)) {
            Expr* denom = t_mul(t_add(mk_int(1), t_mul(expr_copy(qa), expr_copy((Expr*)var))),
                                sqrt_1msq(var));
            Expr* recon = t_mul(expr_copy(Aa), mk_fn2("Power", denom, mk_int(-1)));
            if (is_zero_with(t_add(expr_copy((Expr*)g), t_neg(recon)), assumptions)) {
                A = expr_copy(Aa);
                q = expr_copy(qa);
                ok = true;
            }
        }
        if (Aa) expr_free(Aa);
        if (qa) expr_free(qa);
    }
    /* Form B: r = (1 + q Sqrt[1-var^2])/A ; r|1 = 1/A, r|0 = (1+q)/A. */
    if (!ok) {
        Expr* Ab = simplify_take(mk_fn2("Power", expr_copy(r1), mk_int(-1)), assumptions);
        Expr* qb = simplify_take(t_add(t_mul(expr_copy(Ab), expr_copy(r0)), mk_int(-1)), assumptions);
        if (Ab && qb && is_finite_value(Ab) && is_finite_value(qb) &&
            !contains_symbol(Ab, var) && !contains_symbol(qb, var)) {
            Expr* denom = t_mul(t_add(mk_int(1), t_mul(expr_copy(qb), sqrt_1msq(var))),
                                sqrt_1msq(var));
            Expr* recon = t_mul(expr_copy(Ab), mk_fn2("Power", denom, mk_int(-1)));
            if (is_zero_with(t_add(expr_copy((Expr*)g), t_neg(recon)), assumptions)) {
                A = expr_copy(Ab);
                q = expr_copy(qb);
                ok = true;
            }
        }
        if (Ab) expr_free(Ab);
        if (qb) expr_free(qb);
    }
    expr_free(r0);
    expr_free(r1);
    if (ok) { *out_A = A; *out_q = q; return true; }
    if (A) expr_free(A);
    if (q) expr_free(q);
    return false;
}

/* Families 1 and 2: normalize -> introduce parameter -> differentiate -> ArcCos
 * inner (emitted) -> closed-form back-integration -> verify D[I,p]-J===0 -> I(1). */
static Expr* stage_finite_feynman(const Expr* f, const Expr* x, const Expr* a,
                                  const Expr* b, const Expr* assumptions) {
    Expr* var = NULL;
    Expr* canon = normalize_power_sub(f, x, a, b, &var);
    if (!canon) { if (var) { expr_free(var); var = NULL; }
                  canon = normalize_tan_half(f, x, a, b, &var); }
    if (!canon) return NULL;
    /* A non-integer power of the bare integration variable (t^(5/2), t^a) is the
     * tangent-power family (handled by stage_tangent_power); Simplify of the
     * differentiated radical inner integral would HANG on it. */
    if (has_noninteger_var_power(canon, var)) {
        expr_free(canon);
        if (var) expr_free(var);
        return NULL;
    }

    Expr* p = mk_sym("$diuiFeyn$");
    Expr* fp = insert_feynman_param(canon, p);
    Expr* g = NULL, *A = NULL, *q = NULL, *alpha = NULL, *J = NULL, *G = NULL,
         *I = NULL, *result = NULL;
    if (!fp) goto done;

    /* Base I(p=0)=0: the parameterised integrand must vanish at p=0 (Log[1+0]=0). */
    { Expr* z = mk_int(0); Expr* fp0 = subst(fp, p, z); expr_free(z);
      bool zq = fp0 && is_zero_q(fp0); if (fp0) expr_free(fp0); if (!zq) goto done; }

    g = simplify_take(deriv(fp, p), assumptions);
    if (!g || is_zero_q(g)) goto done;
    if (!inner_arccos_family(g, var, assumptions, &A, &q)) goto done;

    alpha = simplify_take(deriv(q, p), assumptions);           /* dq/dp, must be const */
    if (!alpha || contains_symbol(alpha, p) || contains_symbol(alpha, var) ||
        !is_finite_value(alpha) || is_zero_q(alpha)) goto done;

    J = mk_arccos_over_sqrt(expr_copy(A), expr_copy(q));
    G = mk_G(expr_copy(A), expr_copy(alpha), expr_copy(q));
    { Expr* z = mk_int(0); Expr* G0 = eval_at_param(G, p, z); expr_free(z);
      if (!G0) goto done;
      I = simplify_with(t_add(expr_copy(G), t_neg(G0)), assumptions); }   /* t_neg consumes G0 */
    if (!I || !is_finite_value(I)) goto done;

    /* Verify D[I,p] - J === 0 (catches a back-integration error). */
    { Expr* chk = t_add(deriv(I, p), t_neg(expr_copy(J)));
      bool v = is_zero_with(chk, assumptions); expr_free(chk);
      if (!v) goto done; }

    /* Evaluate at the artificial parameter p = 1. */
    { Expr* one = mk_int(1); Expr* Ival = eval_at_param(I, p, one); expr_free(one);
      if (!Ival) goto done;
      Expr* cl = diui_finalize(Ival, assumptions);
      if (cl) { expr_free(Ival); result = cl; } else result = Ival; }

done:
    if (canon) expr_free(canon);
    if (var) expr_free(var);
    if (p) expr_free(p);
    if (fp) expr_free(fp);
    if (g) expr_free(g);
    if (A) expr_free(A);
    if (q) expr_free(q);
    if (alpha) expr_free(alpha);
    if (J) expr_free(J);
    if (G) expr_free(G);
    if (I) expr_free(I);
    return result;
}

/* Family 3: Csc[2x]^2 Log[1 + Tan[x]^a] on {0,Pi/4}, a direct Beta/digamma
 * evaluation (not a Feynman loop; its differentiated inner integral has no clean
 * closed form).  Normalizes via t = Tan[x] to (1/4)(1+t^2)/t^2 Log[1+t^a] and
 * emits (Pi Csc[Pi/a] - a)/4 (the digamma reflection integral of 1/(u+1)).
 * Correct-by-construction: certified by the value-independent Csc identity and by
 * the exact rational anchor a0 = 3 (the only engine-safe/Simplify-closing one). */
static Expr* stage_tangent_power(const Expr* f, const Expr* x, const Expr* a,
                                 const Expr* b, const Expr* assumptions) {
    Expr* t = NULL;
    Expr* canon = normalize_tan_half(f, x, a, b, &t);
    if (!canon) return NULL;
    Expr* expo = NULL, *target = NULL, *F = NULL, *V = NULL, *result = NULL;

    const Expr* lg = find_log1plus(canon);
    if (!lg) goto done;
    expo = log1plus_monomial_exponent(lg, t);
    if (!expo) goto done;

    /* Recognizer: canon == (1/4)(1+t^2)/t^2 Log[1 + t^expo]. */
    target = mk_tanpow_canon(t, expo);
    { Expr* diff = t_add(expr_copy(canon), t_neg(expr_copy(target)));
      bool m = is_zero_with(diff, assumptions); expr_free(diff);
      if (!m) goto done; }

    /* Value-independent certification of the weight transform. */
    { Expr* lhs = mk_fn2("Power", mk_fn1("Csc", two_x_times(x)), mk_int(2));   /* Csc[2x]^2 */
      Expr* tanx = mk_fn1("Tan", expr_copy((Expr*)x));
      Expr* num  = mk_fn2("Power", t_add(mk_int(1), mk_fn2("Power", expr_copy(tanx), mk_int(2))), mk_int(2));
      Expr* den  = t_mul(mk_int(4), mk_fn2("Power", tanx, mk_int(2)));
      Expr* rhs  = t_mul(num, mk_fn2("Power", den, mk_int(-1)));
      bool idok = is_zero_q(t_add(lhs, t_neg(rhs)));
      if (!idok) goto done; }

    /* Emit the closed form. */
    F = simplify_with(emit_tanpow_form(expo), assumptions);
    if (!F || !is_finite_value(F)) goto done;

    /* Anchor a0 = 3: V = Integrate[(1+t^2)/t^2 Log[1+t^3], {t,0,1}] is an
     * engine-safe Log integrand (no radical/trig), and the nested Integrate
     * re-enters at diui_depth >= DIUI_MAX_DEPTH so DiffUnderInt declines and the
     * ordinary definite engine evaluates it (no recursion, no hang). */
    { Expr* ig = t_mul(one_plus_sq(t),
                       t_mul(mk_fn2("Power", expr_copy(t), mk_int(-2)),
                             mk_fn1("Log", t_add(mk_int(1), mk_fn2("Power", expr_copy(t), mk_int(3))))));
      Expr* spec = mk_fn3("List", expr_copy(t), mk_int(0), mk_int(1));
      V = eval_take(mk_fn2("Integrate", ig, spec)); }
    if (!V || contains_head(V, "Integrate") || !is_finite_value(V)) goto done;
    { Expr* three = mk_int(3); Expr* F3 = emit_tanpow_form(three); expr_free(three);
      Expr* lhs = t_mul(t_rat(1, 4), expr_copy(V));
      bool anchor = is_zero_with(t_add(lhs, t_neg(F3)), assumptions);
      if (!anchor) goto done; }

    result = F; F = NULL;

done:
    if (canon) expr_free(canon);
    if (t) expr_free(t);
    if (expo) expr_free(expo);
    if (target) expr_free(target);
    if (F) expr_free(F);
    if (V) expr_free(V);
    return result;
}

/* =========================================================================
 * Stage B -- first-order linear ODE in the parameter (lambda != 0).
 *
 * When Stage A's inner integral J = Integrate[D[f,p], {x,a,b}] does not close to
 * an x-free form, J may still be expressible AS a multiple of the original
 * integral: J = lambda(p) I + M(p), where I = Integrate[f, {x,a,b}] and lambda,
 * M are x-free.  This holds exactly when, by integration by parts in x,
 *     D[f,p] - lambda(p) f = D[U,x]
 * for some U in the span of f and g's transcendental atoms and an x-free
 * lambda(p); then J = lambda I + [U]_a^b.  That is the first-order linear ODE
 *     I'(p) = lambda(p) I(p) + M(p),   M = [U]_a^b,
 * solved by the integrating factor mu = Exp[-Integral lambda dp].  The classic
 * Feynman Gaussian Integrate[Exp[-a^2 x^2] Cos[b x], {x,0,Inf}] closes this way
 * (differentiating in b gives J = -(b/(2a^2)) I, so I = I(0) Exp[-b^2/(4a^2)]).
 * ---------------------------------------------------------------------- */

/* A maximal transcendental-of-x leaf: Exp/Log/trig/hyperbolic/inverse-trig of an
 * argument containing x, or Exp spelled Power[E, arg(x)].  These are treated as
 * independent indeterminates when matching the IBP identity. */
static bool is_primitive_atom(const Expr* F, const Expr* x) {
    if (!F || F->type != EXPR_FUNCTION || !contains_symbol(F, x)) return false;
    if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
        F->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(F->data.function.args[0]->data.symbol.name, "E") == 0)
        return true;                                   /* Exp[arg] as Power[E,arg] */
    if (F->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = F->data.function.head->data.symbol.name;
    static const char* H[] = { "Exp","Log","Sin","Cos","Tan","Cot","Sec","Csc",
                               "Sinh","Cosh","Tanh","Coth","Sech","Csch",
                               "ArcTan","ArcSin","ArcCos","ArcSinh","ArcCosh","ArcTanh" };
    for (size_t i = 0; i < sizeof(H)/sizeof(H[0]); i++)
        if (strcmp(h, H[i]) == 0) return true;
    return false;
}

/* Append each DISTINCT primitive atom of `e` (borrowed pointers) to out[]. */
static void collect_primitive_atoms(const Expr* e, const Expr* x,
                                    Expr** out, size_t* n, size_t cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (is_primitive_atom(e, x)) {
        for (size_t i = 0; i < *n; i++) if (expr_eq(out[i], e)) return;
        if (*n < cap) out[(*n)++] = (Expr*)e;
        return;                                        /* do not recurse into an atom */
    }
    collect_primitive_atoms(e->data.function.head, x, out, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_primitive_atoms(e->data.function.args[i], x, out, n, cap);
}

/* The atom-product shape of a single multiplicative term: the product of its
 * primitive-atom factors (and Power[atom,k] factors), x-monomial and x-free
 * coefficients stripped.  Owned, canonicalised by evaluation; NULL if the term
 * carries no atom. */
static Expr* atom_shape_of_term(Expr* term, const Expr* x) {
    Expr* facbuf[64];
    size_t nf = collect_factors(term, facbuf, 64, 0);
    Expr* shape = NULL;
    for (size_t i = 0; i < nf; i++) {
        Expr* F = facbuf[i];
        bool keep = is_primitive_atom(F, x);           /* Exp=Power[E,arg], Sin, Cos, ... */
        if (!keep && head_name_is(F, "Power") && F->data.function.arg_count == 2)
            keep = is_primitive_atom(F->data.function.args[0], x);  /* Sin[..]^2, ... */
        if (keep) {
            Expr* fc = expr_copy(F);
            shape = shape ? t_mul(shape, fc) : fc;
        }
    }
    if (!shape) return NULL;
    return eval_take(shape);                            /* canonical factor order */
}

/* Append the DISTINCT atom-product shapes of e's additive terms (owned). */
static void collect_term_shapes(const Expr* e, const Expr* x,
                                Expr** out, size_t* n, size_t cap) {
    Expr* ee = ev1("Expand", expr_copy((Expr*)e));
    if (!ee) return;
    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(ee, "Plus")) { nt = ee->data.function.arg_count; terms = ee->data.function.args; }
    else { nt = 1; single[0] = ee; terms = single; }
    for (size_t t = 0; t < nt; t++) {
        Expr* sh = atom_shape_of_term(terms[t], x);
        if (!sh) continue;
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq(out[i], sh)) { dup = true; break; }
        if (dup || *n >= cap) { expr_free(sh); continue; }
        out[(*n)++] = sh;
    }
    expr_free(ee);
}

/* [U(x,p)]_{x=a}^{x=b}, via direct substitution or Limit at an improper bound.
 * Returns NULL (decline) if either boundary is non-finite or still x-dependent --
 * the IBP identity J = lambda I + [U] is only valid when the boundary converges. */
static Expr* boundary_value(const Expr* U, const Expr* x, const Expr* a,
                            const Expr* b, const Expr* assumptions) {
    Expr* Ub = eval_at_param(U, x, b);
    if (!Ub) return NULL;
    Expr* Ua = eval_at_param(U, x, a);
    if (!Ua) { expr_free(Ub); return NULL; }
    Expr* diff = mk_fn2("Plus", Ub, t_neg(Ua));        /* consumes Ub, Ua */
    Expr* M = simplify_with(diff, assumptions);
    expr_free(diff);
    if (M && is_finite_value(M) && !contains_symbol(M, x)) return M;
    if (M) expr_free(M);
    return NULL;
}

/* True iff `e` contains a negative power of x (1/x^k, Laurent).  The IBP ansatz
 * does not model Laurent integrands (CoefficientList has no negative degrees);
 * such a form is left for the self-similar recognizer. */
static bool has_inv_x_power(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (contains_symbol(base, x)) {
            if (ex->type == EXPR_INTEGER && ex->data.integer < 0) return true;
            if (ex->type == EXPR_REAL && ex->data.real < 0) return true;
            if (ex->type == EXPR_FUNCTION && head_name_is(ex, "Rational") &&
                ex->data.function.args[0]->type == EXPR_INTEGER &&
                ex->data.function.args[0]->data.integer < 0) return true;
        }
    }
    if (has_inv_x_power(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_inv_x_power(e->data.function.args[i], x)) return true;
    return false;
}

/* Detect J = lambda(p) I + M(p) by the IBP ansatz.  On success returns true with
 * owned *out_lambda, *out_M. */
static bool detect_linear_ode(const Expr* f, const Expr* g, const Expr* x,
                              const Expr* a, const Expr* b,
                              const Expr* assumptions,
                              Expr** out_lambda, Expr** out_M) {
    enum { CAP = 16 };
    /* Laurent integrand (1/x^k): the polynomial coefficient-matching below cannot
     * represent it -- leave it for the self-similar recognizer. */
    if (has_inv_x_power(f, x) || has_inv_x_power(g, x)) return false;
    Expr* shapes[CAP]; size_t nsh = 0;
    collect_term_shapes(f, x, shapes, &nsh, CAP);
    collect_term_shapes(g, x, shapes, &nsh, CAP);
    bool ok = false;
    Expr* U = NULL; Expr* lam = NULL;
    if (nsh == 0 || nsh >= CAP) goto done;

    /* U = sum_k c_k shapes[k];  lam a fresh coefficient. */
    for (size_t k = 0; k < nsh; k++) {
        char nm[32]; snprintf(nm, sizeof(nm), "DIUIc%zu", k);
        Expr* term = t_mul(mk_sym(nm), expr_copy(shapes[k]));
        U = U ? t_add(U, term) : term;
    }
    lam = mk_sym("DIUIlam");

    /* R = Expand[g - lam f - D[U,x]], then atoms -> fresh Y indeterminates. */
    {
    Expr* dU = deriv(U, x);
    Expr* R = ev1("Expand",
                  t_add(expr_copy((Expr*)g),
                        t_add(t_neg(t_mul(expr_copy(lam), expr_copy((Expr*)f))),
                              t_neg(dU))));
    /* Collect the primitive atoms of R itself: D[U,x] introduces derivative atoms
     * (Sin from Cos, ...) not present in f/g, and an atom left unsubstituted would
     * carry x into CoefficientList's coefficients.  Collecting from R is complete
     * by construction. */
    Expr* atoms[CAP]; size_t nat = 0;                  /* borrowed from R */
    if (R) collect_primitive_atoms(R, x, atoms, &nat, CAP);
    if (!R || nat == 0 || nat >= CAP) { if (R) expr_free(R); goto done; }
    Expr* rules[CAP];
    for (size_t i = 0; i < nat; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "DIUIy%zu", i);
        rules[i] = mk_fn2("Rule", expr_copy(atoms[i]), mk_sym(nm));
    }
    Expr* rulelist = expr_new_function(mk_sym("List"), rules, nat);
    Expr* R2 = ev2("ReplaceAll", R, rulelist);
    /* Any transcendental-of-x surviving substitution means the basis is
     * incomplete -> CoefficientList would mishandle x; decline. */
    if (R2 && (has_trig_of_x(R2, x) || has_radical_of_x(R2, x) ||
               contains_gaussian_exp(R2, x))) { expr_free(R2); R2 = NULL; }

    /* Coefficient equations of R2 in {x, Y_i}: each must vanish. */
    Expr* vars[CAP + 1]; vars[0] = expr_copy((Expr*)x);
    for (size_t i = 0; i < nat; i++) { char nm[32]; snprintf(nm, sizeof(nm), "DIUIy%zu", i); vars[i+1] = mk_sym(nm); }
    Expr* varlist = expr_new_function(mk_sym("List"), vars, nat + 1);
    Expr* CL = R2 ? ev2("CoefficientList", R2, varlist) : NULL;
    if (!R2) expr_free(varlist);
    Expr* flat = CL ? ev1("Flatten", CL) : NULL;
    Expr* eqsrc = flat ? ev2("DeleteCases", flat, mk_int(0)) : NULL;

    if (eqsrc && head_name_is(eqsrc, "List") && eqsrc->data.function.arg_count > 0) {
        size_t ne = eqsrc->data.function.arg_count;
        Expr** eqitems = (Expr**)malloc(ne * sizeof(Expr*));
        for (size_t i = 0; i < ne; i++)
            eqitems[i] = mk_fn2("Equal", expr_copy(eqsrc->data.function.args[i]), mk_int(0));
        Expr* eqlist = expr_new_function(mk_sym("List"), eqitems, ne);
        free(eqitems);
        Expr** sv = (Expr**)malloc((nsh + 1) * sizeof(Expr*));
        for (size_t k = 0; k < nsh; k++) { char nm[32]; snprintf(nm, sizeof(nm), "DIUIc%zu", k); sv[k] = mk_sym(nm); }
        sv[nsh] = expr_copy(lam);
        Expr* svlist = expr_new_function(mk_sym("List"), sv, nsh + 1);
        free(sv);
        Expr* sol = ev2("Solve", eqlist, svlist);
        if (sol && head_name_is(sol, "List") && sol->data.function.arg_count > 0) {
            Expr* first = sol->data.function.args[0];              /* borrowed */
            Expr* lamr = eval_take(mk_fn2("ReplaceAll", expr_copy(lam), expr_copy(first)));
            Expr* Ur   = eval_take(mk_fn2("ReplaceAll", expr_copy(U),   expr_copy(first)));
            Expr* lamv = lamr ? simplify_with(lamr, assumptions) : NULL;
            Expr* Us   = Ur   ? simplify_with(Ur, assumptions)   : NULL;
            if (lamr) expr_free(lamr);
            if (Ur)   expr_free(Ur);
            /* lambda must be x-free, finite; reject a still-symbolic (c_k-laden)
             * underdetermined solution. */
            bool good = lamv && Us && is_finite_value(lamv) && !contains_symbol(lamv, x) &&
                        !contains_head(lamv, "DIUIlam") && !contains_head(Us, "DIUIlam");
            /* Verify the IBP identity literally (atoms may not be independent). */
            if (good) {
                Expr* idc = t_add(expr_copy((Expr*)g),
                                  t_add(t_neg(t_mul(expr_copy(lamv), expr_copy((Expr*)f))),
                                        t_neg(deriv(Us, x))));
                good = is_zero_with(idc, assumptions);
                expr_free(idc);
            }
            if (good) {
                Expr* M = boundary_value(Us, x, a, b, assumptions);
                if (M) { *out_lambda = lamv; *out_M = M; lamv = NULL; ok = true; }
            }
            if (lamv) expr_free(lamv);
            if (Us) expr_free(Us);
        }
        if (sol) expr_free(sol);
        /* eqlist, svlist were consumed by ev2("Solve", ...) -- do not free them. */
    }
    if (eqsrc) expr_free(eqsrc);
    }

done:
    for (size_t i = 0; i < nsh; i++) expr_free(shapes[i]);
    if (U) expr_free(U);
    if (lam) expr_free(lam);
    return ok;
}

/* Self-similar Gaussian recognizer: f = Exp[c0 + c2 x^2 + cm2/x^2] with c2, cm2
 * both negative.  Under the reciprocal substitution x -> k/x, k = Sqrt[cm2/c2],
 * the exponent is invariant, so Integrate[g, {x,0,Inf}] maps to lambda I with
 * lambda = (g|_{x->k/x}) (k/x^2) / f (x-free when it fires).  M = 0 (the kernel
 * decays at both ends).  Closes Integrate[Exp[-a^2 x^2 - b^2/x^2], {x,0,Inf}]:
 * differentiating in b gives lambda = -2a, so I = I(0) Exp[-2 a b]. */
static bool detect_selfsimilar_ode(const Expr* f, const Expr* g, const Expr* x,
                                   const Expr* assumptions, Expr** out_lambda) {
    const Expr* arg = NULL;
    if (head_name_is(f, "Exp") && f->data.function.arg_count == 1)
        arg = f->data.function.args[0];
    else if (head_name_is(f, "Power") && f->data.function.arg_count == 2 &&
             f->data.function.args[0]->type == EXPR_SYMBOL &&
             strcmp(f->data.function.args[0]->data.symbol.name, "E") == 0)
        arg = f->data.function.args[1];
    if (!arg || !contains_symbol(arg, x)) return false;

    Expr* c2  = ev1("Simplify", mk_fn3("Coefficient", expr_copy((Expr*)arg), expr_copy((Expr*)x), mk_int(2)));
    Expr* ax2 = ev1("Expand", t_mul(expr_copy((Expr*)arg), t_pow(expr_copy((Expr*)x), 2)));
    Expr* cm2 = ev1("Simplify", mk_fn3("Coefficient", ax2, expr_copy((Expr*)x), mk_int(0)));
    bool ok = false;
    /* c0 = arg - c2 x^2 - cm2/x^2; the shape {2,0,-2} is valid iff c0 is x-free
     * (Coefficient[arg,x,0] cannot be used -- it lumps the 1/x^2 term into x^0). */
    Expr* c0 = simplify_with(
        t_add(expr_copy((Expr*)arg),
              t_neg(t_add(t_mul(expr_copy(c2), t_pow(expr_copy((Expr*)x), 2)),
                          t_mul(expr_copy(cm2), t_pow(expr_copy((Expr*)x), -2))))),
        assumptions);
    bool shape_ok = c0 && !contains_symbol(c0, x) &&
                    !contains_symbol(c2, x) && !contains_symbol(cm2, x);
    if (shape_ok) {
        /* k = Sqrt[cm2/c2] = Sqrt[beta/alpha] (both negative -> ratio positive). */
        Expr* k = ev1("Simplify",
                      mk_fn2("Power", t_mul(expr_copy(cm2), t_pow(expr_copy(c2), -1)), t_rat(1, 2)));
        if (k && !contains_symbol(k, x)) {
            Expr* xsub = t_mul(expr_copy(k), t_pow(expr_copy((Expr*)x), -1));   /* k/x */
            Expr* gsub = subst(g, x, xsub); expr_free(xsub);
            /* J = Integrate[g(k/x) (k/x^2), {x,0,Inf}] = lambda I. */
            Expr* lam = NULL;
            if (gsub) {
                Expr* h = t_mul(gsub, t_mul(expr_copy(k), t_pow(expr_copy((Expr*)x), -2)));
                Expr* ratio = t_mul(h, t_pow(expr_copy((Expr*)f), -1));
                lam = simplify_with(ratio, assumptions);
                expr_free(ratio);
            }
            if (lam && !contains_symbol(lam, x) && is_finite_value(lam) && !is_zero_q(lam)) {
                *out_lambda = lam; ok = true;
            } else if (lam) { expr_free(lam); }
        }
        if (k) expr_free(k);
    }
    expr_free(c2); expr_free(c0); expr_free(cm2);
    return ok;
}

/* Solve I'(p) = lambda(p) I(p) + M(p) with I(p0) = I0, by integrating factor.
 * All parameter antiderivatives route through integrate_over_param (bounded), so
 * no unbounded engine call is formed.  Returns owned I, or NULL. */
static Expr* solve_linear_ode(const Expr* lambda, const Expr* M, const Expr* p,
                              const Expr* p0, const Expr* I0,
                              const Expr* assumptions) {
    Expr* L = integrate_over_param(lambda, p);          /* Integral lambda dp */
    if (!L || contains_head(L, "Integrate") || !is_finite_value(L)) {
        if (L) expr_free(L);
        return NULL;
    }
    Expr* L0 = eval_at_param(L, p, p0);
    if (!L0) { expr_free(L); return NULL; }
    Expr* I;
    if (is_zero_q(M)) {
        /* Homogeneous: I = I0 Exp[L - L0]. */
        I = t_mul(expr_copy((Expr*)I0),
                  mk_fn1("Exp", mk_fn2("Plus", expr_copy(L), t_neg(expr_copy(L0)))));
    } else {
        /* Inhomogeneous: I = Exp[L](N - N0) + Exp[L - L0] I0, N = Integral Exp[-L] M dp. */
        Expr* muM = t_mul(mk_fn1("Exp", t_neg(expr_copy(L))), expr_copy((Expr*)M));
        Expr* N = integrate_over_param(muM, p); expr_free(muM);
        if (!N || contains_head(N, "Integrate") || !is_finite_value(N)) {
            if (N) expr_free(N);
            expr_free(L); expr_free(L0);
            return NULL;
        }
        Expr* N0 = eval_at_param(N, p, p0);
        if (!N0) { expr_free(N); expr_free(L); expr_free(L0); return NULL; }
        Expr* term1 = t_mul(mk_fn1("Exp", expr_copy(L)), mk_fn2("Plus", N, t_neg(N0)));
        Expr* term2 = t_mul(mk_fn1("Exp", mk_fn2("Plus", expr_copy(L), t_neg(expr_copy(L0)))),
                            expr_copy((Expr*)I0));
        I = t_add(term1, term2);
    }
    expr_free(L); expr_free(L0);
    Expr* Is = simplify_with(I, assumptions); expr_free(I);
    if (Is && is_finite_value(Is)) return Is;
    if (Is) expr_free(Is);
    return NULL;
}

/* Stage B driver for one parameter: detect the linear ODE, pin the base, solve,
 * and verify D[I,p] - (lambda I + M) === 0. */
static Expr* stage_linear_ode(const Expr* f, const Expr* g, const Expr* x,
                              const Expr* a, const Expr* b, const Expr* assumptions,
                              const Expr* p, const ParamBound* pb, size_t np) {
    Expr* lam = NULL; Expr* M = NULL;
    if (!detect_linear_ode(f, g, x, a, b, assumptions, &lam, &M)) {
        /* Fall back to the self-similar Gaussian recognizer (Exp[-a^2 x^2 - b^2/x^2]),
         * which IBP-in-x cannot find; it yields lambda directly, with M = 0. */
        if (detect_selfsimilar_ode(f, g, x, assumptions, &lam)) M = mk_int(0);
        else return NULL;
    }

    Expr* p0 = NULL; Expr* I0 = NULL;
    if (!find_base(f, x, a, b, assumptions, p, pb, np, false, &p0, &I0)) {
        expr_free(lam); expr_free(M); return NULL;
    }
    Expr* I = solve_linear_ode(lam, M, p, p0, I0, assumptions);
    expr_free(p0); expr_free(I0);
    if (!I) { expr_free(lam); expr_free(M); return NULL; }

    /* D[I,p] - (lambda I + M) === 0 under the assumptions. */
    Expr* chk = t_add(deriv(I, p),
                      t_neg(t_add(t_mul(expr_copy(lam), expr_copy(I)), expr_copy(M))));
    bool ok = is_zero_with(chk, assumptions);
    expr_free(chk); expr_free(lam); expr_free(M);
#ifdef DIUI_DEBUG
    fprintf(stderr, "DIUI:   [%.0fms] Stage B verify : %s\n", diui_ms(), ok ? "PASS" : "FAIL");
#endif
    if (!ok) { expr_free(I); return NULL; }
    Expr* cl = diui_finalize(I, assumptions);
    if (cl) { expr_free(I); return cl; }
    return I;
}

/* -------------------------------------------------------------------------
 * Stage A -- pure quadrature (lambda = 0).
 * ---------------------------------------------------------------------- */
/* One parameter attempt.  `zero_base_only` restricts find_base to an EXACT
 * (sign-independent) zero base; see the two-pass driver below.  Returns an owned
 * verified closed form, or NULL to move on. */
static Expr* stage_quadrature_param(const Expr* f, const Expr* x, const Expr* a,
                                    const Expr* b, const Expr* assumptions,
                                    const Expr* p, const ParamBound* pb, size_t np,
                                    bool zero_base_only) {
    /* g = Simplify[D[f, p]].  Skip if f is independent of p. */
    Expr* g = ev1("Simplify", deriv(f, p));
    if (!g || is_zero_q(g)) { if (g) expr_free(g); return NULL; }

    /* J(p) = Integrate[g, {x,a,b}].  Must close to a finite closed form.
     * Simplify collapses the FTC boundary artifacts (e.g. the lower-limit
     * `0^(1+a)` term in Integrate[x^a,{x,0,1}]) that would otherwise stop the
     * parameter-integration from closing. */
    Expr* J = inner_definite(g, x, a, b, assumptions, pb, np);
    if (J) { Expr* Js = simplify_with(J, assumptions); expr_free(J); J = Js; }
    if (!J || contains_head(J, "Integrate") || !is_finite_value(J)) {
        if (J) expr_free(J);
        /* Stage A's inner integral did not close.  Try Stage B: J = lambda I + M,
         * a first-order linear ODE in the parameter (e.g. Exp[-a^2 x^2] Cos[b x]).
         * Only in the second (computed-base) pass, to avoid redundant work. */
        Expr* r = zero_base_only ? NULL
                 : stage_linear_ode(f, g, x, a, b, assumptions, p, pb, np);
        expr_free(g);
        return r;
    }
    expr_free(g);

    /* G(p) = Integrate[J, p]  (antiderivative over the parameter). */
    Expr* G = integrate_over_param(J, p);
    if (!G || contains_head(G, "Integrate") || !is_finite_value(G)) {
        if (G) expr_free(G);
        expr_free(J);
        return NULL;
    }

    /* Base point with an exact known value I(p0). */
    Expr* p0 = NULL; Expr* I0 = NULL;
    if (!find_base(f, x, a, b, assumptions, p, pb, np, zero_base_only, &p0, &I0)) {
        expr_free(G); expr_free(J); return NULL;
    }

    /* I(p) = G(p) - G(p0) + I(p0). */
    Expr* Gp0 = eval_at_param(G, p, p0);
    if (!Gp0) { expr_free(p0); expr_free(I0); expr_free(G); expr_free(J); return NULL; }

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
    fprintf(stderr, "DIUI:   [%.0fms] verify (zb=%d) D[I,p]-J==0 : %s\n",
            diui_ms(), zero_base_only, ok ? "PASS" : "FAIL");
#endif
    expr_free(J);
    if (ok) {
        Expr* cl = diui_finalize(I, assumptions);
        if (cl) { expr_free(I); return cl; }
        return I;
    }
    if (I) expr_free(I);
    return NULL;
}

static Expr* stage_quadrature(const Expr* f, const Expr* x, const Expr* a,
                              const Expr* b, const Expr* assumptions,
                              const ParamBound* pb, size_t np) {
    /* Two passes over the parameters.  Pass 1 accepts ONLY a parameter with an
     * exact zero base (I0 = 0, sign-independent); pass 2 allows the sign-sensitive
     * computed base.  A parameter whose integrand vanishes at the base point
     * reconstructs the integral cleanly for every parameter sign, whereas a
     * computed base is rendered for one sign branch -- and the D[I,p]-J check
     * cannot catch a wrong constant base, so a zero-base parameter, when one
     * exists, must win.  (E.g. Exp[-a x] Sin[b x] Sin[c x]/x: the b/c paths have
     * the exact base b=0 -> 0, while the a-path's base Integrate[Sin[b x] Sin[c
     * x]/x] is only valid for b>c.) */
    for (int pass = 0; pass < 2; pass++) {
        bool zero_base_only = (pass == 0);
        for (size_t pi = 0; pi < np; pi++) {
            Expr* p = mk_sym(pb[pi].sym);
#ifdef DIUI_DEBUG
            fprintf(stderr, "DIUI: [%.0fms] try param %s (zb=%d)\n",
                    diui_ms(), pb[pi].sym, zero_base_only);
#endif
            Expr* r = stage_quadrature_param(f, x, a, b, assumptions, p, pb, np,
                                             zero_base_only);
            expr_free(p);
            if (r) return r;
        }
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

    /* Free parameters of the integrand (besides x).  np == 0 is NO LONGER a
     * decline: the finite-domain closers below introduce their own artificial
     * parameter, so a purely numeric integrand (families 1 and 2) is handled. */
    ParamBound pb[16];
    size_t np = collect_params(f, x, pb, 16, 0);
    if (assumptions) absorb_fact(pb, np, assumptions);

    /* The method probes many divergent candidate sub-expressions (Limit at a
     * pole, 1/0 in a degenerate family M(s), ...).  Those arithmetic warnings
     * are spurious search noise -- mute them like Limit does, keeping the return
     * values (ComplexInfinity/Indeterminate) that gate the search unchanged. */
    diui_depth++;
    arith_warnings_mute_push();
    /* Finite-domain closers: they canonicalise a trig/radical integrand with a
     * change of variables, introduce an artificial Feynman parameter (or evaluate
     * directly), and supply the inner integral in closed form -- never routing a
     * trig/radical inner integral through the general engine.  Each is a cheap
     * recognizer that declines instantly on a non-match, so the existing half-line
     * paths are unaffected.  stage_tangent_power runs FIRST because it owns the
     * var^(non-integer) family and never forms the radical product that would
     * hang stage_finite_feynman's Simplify on such an integrand.  stage_quadrature
     * (which needs a genuine pre-existing parameter and an engine-safe inner
     * integral) runs last. */
    Expr* result = stage_tangent_power(f, x, a, b, assumptions);
    if (!result) result = stage_finite_feynman(f, x, a, b, assumptions);
    if (!result && np > 0)
        result = stage_quadrature(f, x, a, b, assumptions, pb, np);
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
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                strcmp(lhs->data.symbol.name, "Assumptions") == 0) {
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
        ATTR_PROTECTED;
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
