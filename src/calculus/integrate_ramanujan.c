/* integrate_ramanujan.c
 *
 * Half-line definite integration Integrate[f, {x, 0, Infinity}] by the
 * Mellin-transform / Ramanujan Master Theorem method.  See
 * integrate_ramanujan.h for the overview.
 *
 * Strategy (Layer 1 -- table + operational rules):
 *   Expand the integrand into additive terms.  Each term is decomposed into
 *       C * x^rho * kernel(x)
 *   where C is x-free, x^rho is a pure power (rho a number or a symbol), and
 *   kernel is a single transcendental/algebraic factor matched against a table
 *   of proven base Mellin transforms M_f(s) = Integrate[x^(s-1) f, {x,0,Inf}]:
 *
 *       Exp[c x]              -> Gamma(s) (-c)^(-s)            0 < Re s
 *       Exp[c x^2]            -> (1/2) (-c)^(-s/2) Gamma(s/2)  0 < Re s
 *       (p + q x^m)^(-a)      -> (1/m) p^(s/m-a) q^(-s/m)
 *                                 Beta(s/m, a - s/m)          0 < Re s < m Re a
 *       Cos[lam x]           -> Gamma(s) Cos(pi s/2) lam^(-s) 0 < Re s < 1
 *       Sin[lam x]           -> Gamma(s) Sin(pi s/2) lam^(-s) -1 < Re s < 1
 *       BesselJ[nu, lam x]   -> 2^(s-1) lam^(-s)
 *                                 Gamma((nu+s)/2)/Gamma((nu-s)/2+1)
 *
 *   The power prefactor sets s = rho + 1; a linear internal scaling folds into
 *   the transform (the lam^(-s) / (-c)^(-s) factors).  Each application is
 *   gated on its convergence strip (checked by Simplify: numerically for a
 *   numeric s, or against the supplied Assumptions for a symbolic s), so the
 *   result is unconditionally correct.  A removable 0*Infinity coincidence at a
 *   numeric s (e.g. Sin[x]/x at s=0) is resolved by a Limit fallback.
 *
 * Terms are summed independently: each must converge on (0, Infinity) on its
 * own (absolute convergence of each term => the sum of integrals equals the
 * integral of the sum).  A product of two transcendental kernels is a Mellin
 * convolution (Meijer-G territory) and is out of scope -> NULL.
 *
 * Verification is symbolic and correct-by-construction: every table identity is
 * a theorem and the strip gate guarantees convergence.  No NIntegrate anywhere
 * (project rule).
 */

#include "integrate_ramanujan.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "parse.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* -------------------------------------------------------------------------
 * Small expression-construction / evaluation helpers (mirrors the idiom in
 * integrate_diffunderint.c / integrate_newton_leibniz.c).
 * ---------------------------------------------------------------------- */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_int(long v)        { return expr_new_integer((int64_t)v); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}

/* Evaluate `call`, free the call expression, return the (owned) result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}
static Expr* ev1(const char* name, Expr* a)          { return eval_take(mk_fn1(name, a)); }
static Expr* ev2(const char* name, Expr* a, Expr* b) { return eval_take(mk_fn2(name, a, b)); }

/* Const-friendly deep copy: recognizers hold borrowed const subexpressions. */
static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

/* Compact algebra builders (each consumes its Expr* arguments). */
static Expr* Gamma_(Expr* z)          { return mk_fn1("Gamma", z); }
static Expr* Beta_(Expr* a, Expr* b)  { return mk_fn2("Beta", a, b); }
static Expr* Pw(Expr* b, Expr* e)     { return mk_fn2("Power", b, e); }
static Expr* Tms(Expr* a, Expr* b)    { return mk_fn2("Times", a, b); }
static Expr* Tms3(Expr* a, Expr* b, Expr* c) { return Tms(a, Tms(b, c)); }
static Expr* Pls(Expr* a, Expr* b)    { return mk_fn2("Plus", a, b); }
static Expr* Neg(Expr* e)             { return Tms(mk_int(-1), e); }
static Expr* Rat(long p, long q)      { return mk_fn2("Rational", mk_int(p), mk_int(q)); }
static Expr* Half(Expr* e)            { return Tms(Rat(1, 2), e); }   /* e/2 */
static Expr* ExpE(Expr* z)            { return Pw(mk_sym("E"), z); }
static Expr* Gt(Expr* a, Expr* b)     { return mk_fn2("Greater", a, b); }
static Expr* Lt(Expr* a, Expr* b)     { return mk_fn2("Less", a, b); }
static Expr* And2(Expr* a, Expr* b)   { return mk_fn2("And", a, b); }

static Expr* simp(Expr* e) { return ev1("Simplify", e); }
static Expr* simp2(Expr* e, Expr* as) {
    return as ? ev2("Simplify", e, cp(as)) : ev1("Simplify", e);
}
/* Structural size, for choosing the tidier of two equivalent forms. */
static long leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    long n = 1 + leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n += leaf_count(e->data.function.args[i]);
    return n;
}

/* Final cleanup.  The pFq reductions can leave Gamma-ratio artifacts such as
 * Gamma[1+a]/Gamma[a]; the bounded, terminating rewrite
 * Gamma[z+n] -> Pochhammer[z,n] Gamma[z] exposes them for Simplify to cancel.
 * But that same rewrite un-consolidates an isolated Gamma[1+nu-s] into a bulkier
 * product when there is nothing to cancel.  So simplify both with and without
 * the rewrite and keep whichever is structurally smaller -- the cancellation
 * cases shrink, the isolated-Gamma cases keep their compact form.  (FullSimplify
 * would subsume this but can hang on symbolic-exponent forms like (1-s)^(-nu).) */
static Expr* fullsimp2(Expr* e, Expr* as) {
    Expr* plain = as ? ev2("Simplify", cp(e), cp(as)) : ev1("Simplify", cp(e));
    static const char* GRULE =
        "Gamma[z_ + n_Integer] :> Pochhammer[z, n] Gamma[z]";
    Expr* rule = parse_expression(GRULE);
    Expr* norm = NULL;
    if (rule) {
        Expr* repl = ev2("ReplaceRepeated", cp(e), rule);
        norm = as ? ev2("Simplify", repl, cp(as)) : ev1("Simplify", repl);
    }
    expr_free(e);
    if (norm && (!plain || leaf_count(norm) < leaf_count(plain))) {
        if (plain) expr_free(plain);
        return norm;
    }
    if (norm) expr_free(norm);
    return plain;
}

/* True iff `e` is the compound `name[...]` (by head name). */
static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
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
static bool contains_symbol(const Expr* e, const Expr* x) {
    return x && x->type == EXPR_SYMBOL && contains_symbol_name(e, x->data.symbol.name);
}

/* Non-consuming: Simplify[e] is x-free. */
static bool free_of_x_now(const Expr* e, const Expr* x) {
    Expr* v = simp(cp((Expr*)e));
    bool f = v && !contains_symbol(v, x);
    if (v) expr_free(v);
    return f;
}
/* Non-consuming: Simplify[e] is the literal 0. */
static bool is_zero_now(const Expr* e) {
    Expr* v = simp(cp((Expr*)e));
    bool z = v && ((v->type == EXPR_INTEGER && v->data.integer == 0) ||
                   (v->type == EXPR_REAL && v->data.real == 0.0));
    if (v) expr_free(v);
    return z;
}
/* Consuming: Simplify[pred] (optionally with assumptions) is the symbol True. */
static bool prove_true(Expr* pred, Expr* as) {
    Expr* v = as ? ev2("Simplify", pred, cp(as)) : ev1("Simplify", pred);
    bool t = v && v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "True") == 0;
    if (v) expr_free(v);
    return t;
}

/* D[e, x] (non-consuming inputs; owned result). */
static Expr* dx(const Expr* e, const Expr* x) {
    return ev2("D", cp((Expr*)e), cp((Expr*)x));
}

/* A finite closed form: non-NULL and free of the infinite/indeterminate
 * sentinels the evaluator produces for divergent Gamma/Power evaluations. */
static bool is_finite_value(const Expr* e) {
    if (!e) return false;
    if (contains_symbol_name(e, "Indeterminate")) return false;
    if (contains_symbol_name(e, "ComplexInfinity")) return false;
    if (contains_symbol_name(e, "Infinity")) return false;   /* covers +/-Infinity */
    if (head_name_is(e, "DirectedInfinity")) return false;
    /* An unevaluated transform (a leftover Limit / SeriesCoefficient) is not a
     * closed value -- never return it as if finite. */
    if (contains_symbol_name(e, "Limit")) return false;
    if (contains_symbol_name(e, "SeriesCoefficient")) return false;
    return true;
}

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

/* Collect multiplicative leaf factors of `T`, flattening nested Times.  Returns
 * the count written into out[] (up to cap); factors are borrowed. */
static size_t collect_factors(Expr* T, Expr** out, size_t cap, size_t n) {
    if (head_name_is(T, "Times")) {
        for (size_t i = 0; i < T->data.function.arg_count && n < cap; i++)
            n = collect_factors(T->data.function.args[i], out, cap, n);
        return n;
    }
    if (n < cap) out[n++] = T;
    return n;
}

/* If F is x or Power[x, e] with e x-free, return an owned copy of the exponent;
 * otherwise NULL. */
static Expr* x_power_exponent(const Expr* F, const Expr* x) {
    if (F->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
        strcmp(F->data.symbol.name, x->data.symbol.name) == 0)
        return mk_int(1);
    if (head_name_is(F, "Power") && F->data.function.arg_count == 2) {
        Expr* base = F->data.function.args[0];
        Expr* ex   = F->data.function.args[1];
        if (base->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
            strcmp(base->data.symbol.name, x->data.symbol.name) == 0 &&
            free_of_x_now(ex, x))
            return cp(ex);
    }
    return NULL;
}

/* The exponent argument of an exponential factor: Exp[arg] or Power[E, arg].
 * Returns a borrowed pointer, or NULL. */
static Expr* exp_exponent(const Expr* K) {
    if (head_name_is(K, "Exp") && K->data.function.arg_count == 1)
        return K->data.function.args[0];
    if (head_name_is(K, "Power") && K->data.function.arg_count == 2) {
        Expr* base = K->data.function.args[0];
        if (base->type == EXPR_SYMBOL && strcmp(base->data.symbol.name, "E") == 0)
            return K->data.function.args[1];
    }
    return NULL;
}

/* -------------------------------------------------------------------------
 * Base Mellin-transform recognizers.
 *
 * Each takes the kernel factor K, the variable x, and a "spectral variable"
 * sv (the transform argument s).  On a match it sets *M to M_f(sv) and *P to
 * the convergence predicate in terms of sv, both owned, and returns true.
 * sv is passed in (rather than fixed to s) so the caller can rebuild the
 * transform in a fresh symbol for the removable-point Limit fallback.
 * ---------------------------------------------------------------------- */

/* Exp[c x], Re(c) < 0:  M = e^{c0} Gamma(sv) (-c)^(-sv),  0 < Re sv. */
static bool rec_exp(const Expr* K, const Expr* x, const Expr* sv,
                    Expr** M, Expr** P) {
    Expr* arg = exp_exponent(K);
    if (!arg) return false;
    Expr* c1 = dx(arg, x);                                  /* linear coeff */
    if (!free_of_x_now(c1, x)) { expr_free(c1); return false; }
    Expr* c0 = simp(Pls(cp(arg), Neg(Tms(cp(c1), cp(x)))));
    if (!free_of_x_now(c0, x)) { expr_free(c1); expr_free(c0); return false; }
    *M = Tms3(ExpE(cp(c0)), Gamma_(cp(sv)),
              Pw(Neg(cp(c1)), Neg(cp(sv))));
    /* Convergence Re(c) < 0, phrased as the positive decay rate -c > 0 so weak
     * assumption reasoning (a > 0) discharges it. */
    *P = And2(Gt(cp(sv), mk_int(0)), Gt(Neg(cp(c1)), mk_int(0)));
    expr_free(c1); expr_free(c0);
    return true;
}

/* Exp[c x^2], Re(c) < 0:  M = e^{c0} (1/2) (-c)^(-sv/2) Gamma(sv/2). */
static bool rec_gauss(const Expr* K, const Expr* x, const Expr* sv,
                      Expr** M, Expr** P) {
    Expr* arg = exp_exponent(K);
    if (!arg) return false;
    Expr* c1d = dx(arg, x);
    Expr* c2d = dx(c1d, x);                                 /* = 2 c2 */
    if (!free_of_x_now(c2d, x) || is_zero_now(c2d)) {       /* not a pure quadratic */
        expr_free(c1d); expr_free(c2d); return false;
    }
    /* No linear cross term (a genuine e^{b x} shift would need Erf). */
    Expr* cross = simp(Pls(cp(c1d), Neg(Tms(cp(c2d), cp(x)))));
    bool bad = !is_zero_now(cross);
    expr_free(cross);
    if (bad) { expr_free(c1d); expr_free(c2d); return false; }
    Expr* c2 = simp(Half(cp(c2d)));
    Expr* c0 = simp(Pls(cp(arg),
                        Neg(Tms(cp(c2), Pw(cp(x), mk_int(2))))));
    if (!free_of_x_now(c0, x)) {
        expr_free(c1d); expr_free(c2d); expr_free(c2); expr_free(c0); return false;
    }
    *M = Tms(ExpE(cp(c0)),
             Tms3(Rat(1, 2), Pw(Neg(cp(c2)), Neg(Half(cp(sv)))),
                  Gamma_(Half(cp(sv)))));
    /* Convergence Re(c) < 0, as the positive rate -c > 0 (see rec_exp). */
    *P = And2(Gt(cp(sv), mk_int(0)), Gt(Neg(cp(c2)), mk_int(0)));
    expr_free(c1d); expr_free(c2d); expr_free(c2); expr_free(c0);
    return true;
}

/* (p + q x^m)^(-a), p,q > 0, m > 0:
 *   M = (1/m) p^(sv/m - a) q^(-sv/m) Beta(sv/m, a - sv/m),  0 < Re sv < m Re a. */
static bool rec_algebraic(const Expr* K, const Expr* x, const Expr* sv,
                          Expr** M, Expr** P) {
    if (!head_name_is(K, "Power") || K->data.function.arg_count != 2) return false;
    Expr* base = K->data.function.args[0];
    Expr* exp  = K->data.function.args[1];
    if (!free_of_x_now(exp, x)) return false;
    if (!head_name_is(base, "Plus")) return false;         /* need a constant term */
    Expr* a = simp(Neg(cp(exp)));                   /* a = -exp */

    /* Split base = (x-free part p) + (single x-part). */
    Expr* p = mk_int(0);
    Expr* xpart = NULL;
    for (size_t i = 0; i < base->data.function.arg_count; i++) {
        Expr* ai = base->data.function.args[i];
        if (free_of_x_now(ai, x)) { p = Pls(p, cp(ai)); continue; }
        if (xpart) { expr_free(p); expr_free(a); return false; }   /* two x-parts */
        xpart = ai;                                        /* borrow */
    }
    p = simp(p);
    if (!xpart || is_zero_now(p)) { expr_free(p); expr_free(a); return false; }

    /* Parse xpart = q x^m. */
    Expr* facs[32];
    size_t nf = collect_factors(xpart, facs, 32, 0);
    Expr* q = mk_int(1);
    Expr* m = NULL;
    for (size_t i = 0; i < nf; i++) {
        Expr* F = facs[i];
        if (free_of_x_now(F, x)) { q = Tms(q, cp(F)); continue; }
        Expr* e = x_power_exponent(F, x);
        if (!e || m) { if (e) expr_free(e); expr_free(q); if (m) expr_free(m);
                       expr_free(p); expr_free(a); return false; }
        m = e;
    }
    q = simp(q);
    if (!m || !prove_true(Gt(cp(m), mk_int(0)), NULL)) {
        if (m) { expr_free(m); } expr_free(q); expr_free(p); expr_free(a); return false;
    }

    Expr* sm = simp(Tms(cp(sv), Pw(cp(m), mk_int(-1))));   /* s/m */
    *M = Tms(Pw(cp(m), mk_int(-1)),
             Tms3(Pw(cp(p), Pls(cp(sm), Neg(cp(a)))),
                  Pw(cp(q), Neg(cp(sm))),
                  Beta_(cp(sm), Pls(cp(a), Neg(cp(sm))))));
    *P = And2(And2(Gt(cp(p), mk_int(0)), Gt(cp(q), mk_int(0))),
              And2(Gt(cp(sm), mk_int(0)), Lt(cp(sm), cp(a))));
    expr_free(p); expr_free(a); expr_free(q); expr_free(m); expr_free(sm);
    return true;
}

/* Cos[lam x] / Sin[lam x], lam > 0:
 *   Sin: M = pi / (2 Cos(pi sv/2) Gamma(1-sv)) lam^(-sv),  -1 < Re sv < 1
 *   Cos: M = pi / (2 Sin(pi sv/2) Gamma(1-sv)) lam^(-sv),   0 < Re sv < 1
 * These are the reflection-formula rewrites of Gamma(sv) Sin/Cos(pi sv/2):
 * using Gamma(sv) = pi/(Sin(pi sv) Gamma(1-sv)) and Sin(pi sv) =
 * 2 Sin(pi sv/2) Cos(pi sv/2), the transform becomes a form with NO removable
 * 0*Infinity coincidence (the Sin[x]/x value pi/2 falls out directly at sv=0,
 * where Gamma(sv) Sin(pi sv/2) would need a limit the engine cannot take). */
static bool rec_trig(const Expr* K, const Expr* x, const Expr* sv,
                     Expr** M, Expr** P) {
    bool is_sin = head_name_is(K, "Sin");
    bool is_cos = head_name_is(K, "Cos");
    if ((!is_sin && !is_cos) || K->data.function.arg_count != 1) return false;
    Expr* arg = K->data.function.args[0];
    Expr* lam = dx(arg, x);
    if (!free_of_x_now(lam, x)) { expr_free(lam); return false; }
    Expr* cross = simp(Pls(cp(arg), Neg(Tms(cp(lam), cp(x)))));
    bool bad = !is_zero_now(cross);                        /* no constant phase */
    expr_free(cross);
    if (bad) { expr_free(lam); return false; }
    Expr* trig = is_sin                                    /* denominator trig */
        ? mk_fn1("Cos", Half(Tms(mk_sym("Pi"), cp(sv))))
        : mk_fn1("Sin", Half(Tms(mk_sym("Pi"), cp(sv))));
    Expr* denom = Tms3(mk_int(2), trig, Gamma_(Pls(mk_int(1), Neg(cp(sv)))));
    *M = Tms3(mk_sym("Pi"), Pw(denom, mk_int(-1)), Pw(cp(lam), Neg(cp(sv))));
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), mk_int(is_sin ? -1 : 0)),
                   Lt(cp(sv), mk_int(1))));
    expr_free(lam);
    return true;
}

/* BesselJ[nu, lam x], lam > 0:
 *   M = 2^(sv-1) lam^(-sv) Gamma((nu+sv)/2)/Gamma((nu-sv)/2+1);
 *   -Re nu < Re sv < 3/2. */
static bool rec_bessel(const Expr* K, const Expr* x, const Expr* sv,
                       Expr** M, Expr** P) {
    if (!head_name_is(K, "BesselJ") || K->data.function.arg_count != 2) return false;
    Expr* nu  = K->data.function.args[0];
    Expr* arg = K->data.function.args[1];
    if (!free_of_x_now(nu, x)) return false;
    Expr* lam = dx(arg, x);
    if (!free_of_x_now(lam, x)) { expr_free(lam); return false; }
    Expr* cross = simp(Pls(cp(arg), Neg(Tms(cp(lam), cp(x)))));
    bool bad = !is_zero_now(cross);
    expr_free(cross);
    if (bad) { expr_free(lam); return false; }
    *M = Tms3(Pw(mk_int(2), Pls(cp(sv), mk_int(-1))),
              Pw(cp(lam), Neg(cp(sv))),
              Tms(Gamma_(Half(Pls(cp(nu), cp(sv)))),
                  Pw(Gamma_(Pls(Half(Pls(cp(nu), Neg(cp(sv)))),
                                mk_int(1))),
                     mk_int(-1))));
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), Neg(cp(nu))),
                   Lt(cp(sv), Rat(3, 2))));
    expr_free(lam);
    return true;
}

/* Log[1 + lam x], lam > 0:  M = pi/(sv Sin(pi sv)) lam^(-sv),  -1 < Re sv < 0.
 * (From Integrate[x^s/(1+x)] = -pi/Sin(pi s) by parts.) */
static bool rec_log(const Expr* K, const Expr* x, const Expr* sv,
                    Expr** M, Expr** P) {
    if (!head_name_is(K, "Log") || K->data.function.arg_count != 1) return false;
    Expr* arg = K->data.function.args[0];
    if (!head_name_is(arg, "Plus")) return false;             /* need 1 + lam x */
    Expr* c = mk_int(0);
    Expr* xpart = NULL;
    for (size_t i = 0; i < arg->data.function.arg_count; i++) {
        Expr* ai = arg->data.function.args[i];
        if (free_of_x_now(ai, x)) { c = Pls(c, cp(ai)); continue; }
        if (xpart) { expr_free(c); return false; }
        xpart = ai;
    }
    c = simp(c);
    bool c_is_one = c->type == EXPR_INTEGER && c->data.integer == 1;
    expr_free(c);
    if (!xpart || !c_is_one) return false;                    /* constant term must be 1 */
    Expr* lam = dx(xpart, x);
    if (!free_of_x_now(lam, x)) { expr_free(lam); return false; }
    Expr* cross = simp(Pls(cp(xpart), Neg(Tms(cp(lam), cp(x)))));
    bool bad = !is_zero_now(cross);                           /* xpart == lam x, linear */
    expr_free(cross);
    if (bad) { expr_free(lam); return false; }
    Expr* denom = Tms(cp(sv), mk_fn1("Sin", Tms(mk_sym("Pi"), cp(sv))));
    *M = Tms3(mk_sym("Pi"), Pw(denom, mk_int(-1)), Pw(cp(lam), Neg(cp(sv))));
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), mk_int(-1)), Lt(cp(sv), mk_int(0))));
    expr_free(lam);
    return true;
}

/* ArcTan[lam x], lam > 0:  M = -pi/(2 sv Cos(pi sv/2)) lam^(-sv),  -1 < Re sv < 0. */
static bool rec_arctan(const Expr* K, const Expr* x, const Expr* sv,
                       Expr** M, Expr** P) {
    if (!head_name_is(K, "ArcTan") || K->data.function.arg_count != 1) return false;
    Expr* arg = K->data.function.args[0];
    Expr* lam = dx(arg, x);
    if (!free_of_x_now(lam, x)) { expr_free(lam); return false; }
    Expr* cross = simp(Pls(cp(arg), Neg(Tms(cp(lam), cp(x)))));
    bool bad = !is_zero_now(cross);                           /* no constant phase */
    expr_free(cross);
    if (bad) { expr_free(lam); return false; }
    Expr* denom = Tms3(mk_int(2), cp(sv),
                       mk_fn1("Cos", Half(Tms(mk_sym("Pi"), cp(sv)))));
    *M = Tms3(mk_int(-1),
              Tms(mk_sym("Pi"), Pw(denom, mk_int(-1))),
              Pw(cp(lam), Neg(cp(sv))));
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), mk_int(-1)), Lt(cp(sv), mk_int(0))));
    expr_free(lam);
    return true;
}

/* HypergeometricPFQ[{a1..ap}, {b1..bq}, -lam x], lam > 0, p <= q+1:
 *   M = (prod_j Gamma(b_j) / prod_i Gamma(a_i)) Gamma(sv)
 *        (prod_i Gamma(a_i - sv) / prod_j Gamma(b_j - sv)) lam^(-sv),
 *   0 < Re sv < min_i Re a_i.
 * This is the Ramanujan Master Theorem applied to the pFq power series:
 *   pFq(-lam x) = sum (-1)^k phi(k) (lam x)^k / k!,  phi(k) = prod (a_i)_k / prod (b_j)_k,
 *   giving  Integrate[x^(sv-1) pFq] = Gamma(sv) phi(-sv) lam^(-sv). */
static bool rec_pfq(const Expr* K, const Expr* x, const Expr* sv,
                    Expr** M, Expr** P) {
    if (!head_name_is(K, "HypergeometricPFQ") || K->data.function.arg_count != 3)
        return false;
    Expr* as  = K->data.function.args[0];
    Expr* bs  = K->data.function.args[1];
    Expr* arg = K->data.function.args[2];
    if (!head_name_is(as, "List") || !head_name_is(bs, "List")) return false;
    if (as->data.function.arg_count > bs->data.function.arg_count + 1) return false;
    if (!free_of_x_now(as, x) || !free_of_x_now(bs, x)) return false;
    Expr* c1 = dx(arg, x);
    if (!free_of_x_now(c1, x)) { expr_free(c1); return false; }
    Expr* cross = simp(Pls(cp(arg), Neg(Tms(cp(c1), cp(x)))));
    bool bad = !is_zero_now(cross);                           /* arg == c1 x, linear */
    expr_free(cross);
    if (bad) { expr_free(c1); return false; }
    Expr* lam = simp(Neg(cp(c1)));                            /* arg = -lam x, lam > 0 */
    expr_free(c1);

    Expr* prod = Gamma_(cp(sv));                              /* Gamma(sv) */
    Expr* smin = NULL;                                        /* strip: sv < a_i, all i */
    for (size_t i = 0; i < as->data.function.arg_count; i++) {
        Expr* ai = as->data.function.args[i];
        prod = Tms(prod, Gamma_(Pls(cp(ai), Neg(cp(sv)))));   /* Gamma(a_i - sv) */
        prod = Tms(prod, Pw(Gamma_(cp(ai)), mk_int(-1)));     /* / Gamma(a_i)    */
        Expr* ci = Lt(cp(sv), cp(ai));
        smin = smin ? And2(smin, ci) : ci;
    }
    for (size_t j = 0; j < bs->data.function.arg_count; j++) {
        Expr* bj = bs->data.function.args[j];
        prod = Tms(prod, Gamma_(cp(bj)));                     /* Gamma(b_j)      */
        prod = Tms(prod, Pw(Gamma_(Pls(cp(bj), Neg(cp(sv)))), mk_int(-1)));
    }                                                         /* / Gamma(b_j-sv) */
    *M = Tms(prod, Pw(cp(lam), Neg(cp(sv))));
    Expr* pos = Gt(cp(sv), mk_int(0));
    *P = And2(Gt(cp(lam), mk_int(0)), smin ? And2(pos, smin) : pos);
    expr_free(lam);
    return true;
}

/* PolyLog[nu, -lam x], lam > 0:  M = pi (-sv)^(-nu) lam^(-sv) / Sin(pi sv),
 *   -1 < Re sv < 0.
 * From the Fermi-Dirac / Bose integral representation
 *   -Li_nu(-x) = (1/Gamma(nu)) Integrate[t^(nu-1) x/(e^t+x), {t,0,Inf}]
 * with Integrate[x^s/(beta+x), {x,0,Inf}] = -beta^s pi/Sin(pi s).  Consistency:
 * at nu=1, -Li_1(-x)=Log(1+x) and the transform reduces to rec_log's. */
static bool rec_polylog(const Expr* K, const Expr* x, const Expr* sv,
                        Expr** M, Expr** P) {
    if (!head_name_is(K, "PolyLog") || K->data.function.arg_count != 2) return false;
    Expr* nu  = K->data.function.args[0];
    Expr* arg = K->data.function.args[1];
    if (!free_of_x_now(nu, x)) return false;
    Expr* c1 = dx(arg, x);
    if (!free_of_x_now(c1, x)) { expr_free(c1); return false; }
    Expr* cross = simp(Pls(cp(arg), Neg(Tms(cp(c1), cp(x)))));
    bool bad = !is_zero_now(cross);                           /* arg == c1 x, linear */
    expr_free(cross);
    if (bad) { expr_free(c1); return false; }
    Expr* lam = simp(Neg(cp(c1)));                            /* arg = -lam x, lam > 0 */
    expr_free(c1);
    Expr* denom = mk_fn1("Sin", Tms(mk_sym("Pi"), cp(sv)));
    *M = Tms(mk_sym("Pi"),
             Tms(Pw(Neg(cp(sv)), Neg(cp(nu))),
                 Tms(Pw(cp(lam), Neg(cp(sv))), Pw(denom, mk_int(-1)))));
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), mk_int(-1)), Lt(cp(sv), mk_int(0))));
    expr_free(lam);
    return true;
}

/* -------------------------------------------------------------------------
 * Interval-arithmetic bound prover -- the assumption gate for a SYMBOLIC
 * fugacity in the exponential-geometric kernel below.
 *
 * The built-in assumption engine (Simplify[pred, as]) discharges only a
 * predicate that syntactically matches an assumption atom: from `0 < z < 1` it
 * proves neither `1/z > 0` (reciprocal sign) nor `-1 <= -z <= 1` (strict ->
 * nonstrict + transitivity).  So no reformulation of prove_true closes the
 * general Bose integral 1/(z^-1 e^x - 1) with a symbolic z.
 *
 * Instead we prove the two sign / interval conditions on A_eff and gamma' by a
 * small, sound interval evaluation over the parameter box read off the
 * Assumptions.  Interval arithmetic yields a guaranteed ENCLOSURE of the true
 * value set, so a proof is never a false positive; a dependency (a symbol
 * occurring more than once) only widens the enclosure, which can at worst make
 * the gate decline -- it can never wrongly accept.  Endpoint openness is
 * tracked (and propagated conservatively) so a strict bound `gamma' > -1` --
 * needed to exclude the Bose boundary z = 1 where the convergence strip changes
 * -- can be certified from an open assumption box.
 * ---------------------------------------------------------------------- */

typedef struct { double lo, hi; bool lo_open, hi_open; } Iv;

/* Machine value of a concrete numeric expression (Integer/Real/Rational/...),
 * else false for a symbolic one. */
static bool num_dbl(const Expr* e, double* out) {
    Expr* n = ev1("N", cp(e));
    bool ok = false;
    if (n) {
        if (n->type == EXPR_INTEGER) { *out = (double)n->data.integer; ok = true; }
        else if (n->type == EXPR_REAL) { *out = n->data.real; ok = true; }
        expr_free(n);
    }
    return ok;
}

/* Tighten `iv` by one binary relation.  `op` is a comparison head; the relation
 * is (sym op c) when sym_on_left, else (c op sym). */
static void iv_tighten(Iv* iv, const char* op, double c, bool sym_on_left) {
    bool lt = (strcmp(op, "Less") == 0 || strcmp(op, "LessEqual") == 0);
    bool gt = (strcmp(op, "Greater") == 0 || strcmp(op, "GreaterEqual") == 0);
    bool strict = (strcmp(op, "Less") == 0 || strcmp(op, "Greater") == 0);
    if (!lt && !gt) return;
    bool upper = sym_on_left ? lt : gt;              /* sym < c => c bounds above */
    if (upper) {
        if (c < iv->hi) { iv->hi = c; iv->hi_open = strict; }
        else if (c == iv->hi && strict) iv->hi_open = true;
    } else {
        if (c > iv->lo) { iv->lo = c; iv->lo_open = strict; }
        else if (c == iv->lo && strict) iv->lo_open = true;
    }
}

/* Read the tightest box [lo, hi] the Assumptions impose on the bare symbol
 * `sym` (And/List conjunctions, chained Inequality, the four binary relations).
 * Returns the box; an absent bound stays +-inf (open at infinity). */
static Iv sym_bound(const Expr* fact, const char* sym) {
    Iv iv = { -HUGE_VAL, HUGE_VAL, true, true };
    if (!fact || fact->type != EXPR_FUNCTION) return iv;
    const Expr* head = fact->data.function.head;
    if (head->type != EXPR_SYMBOL) return iv;
    const char* h = head->data.symbol.name;
    size_t ac = fact->data.function.arg_count;
    if (strcmp(h, "And") == 0 || strcmp(h, "List") == 0) {
        for (size_t i = 0; i < ac; i++) {
            Iv sub = sym_bound(fact->data.function.args[i], sym);
            if (sub.lo > iv.lo) { iv.lo = sub.lo; iv.lo_open = sub.lo_open; }
            else if (sub.lo == iv.lo && sub.lo_open) iv.lo_open = true;
            if (sub.hi < iv.hi) { iv.hi = sub.hi; iv.hi_open = sub.hi_open; }
            else if (sub.hi == iv.hi && sub.hi_open) iv.hi_open = true;
        }
        return iv;
    }
    if (strcmp(h, "Inequality") == 0 && ac >= 3) {
        for (size_t i = 0; i + 2 < ac; i += 2) {
            Expr* L = fact->data.function.args[i];
            Expr* opE = fact->data.function.args[i + 1];
            Expr* R = fact->data.function.args[i + 2];
            double c;
            if (opE->type != EXPR_SYMBOL) continue;
            if (L->type == EXPR_SYMBOL && strcmp(L->data.symbol.name, sym) == 0 && num_dbl(R, &c))
                iv_tighten(&iv, opE->data.symbol.name, c, true);
            else if (R->type == EXPR_SYMBOL && strcmp(R->data.symbol.name, sym) == 0 && num_dbl(L, &c))
                iv_tighten(&iv, opE->data.symbol.name, c, false);
        }
        return iv;
    }
    if (ac == 2) {
        Expr* L = fact->data.function.args[0];
        Expr* R = fact->data.function.args[1];
        double c;
        if (L->type == EXPR_SYMBOL && strcmp(L->data.symbol.name, sym) == 0 && num_dbl(R, &c))
            iv_tighten(&iv, h, c, true);
        else if (R->type == EXPR_SYMBOL && strcmp(R->data.symbol.name, sym) == 0 && num_dbl(L, &c))
            iv_tighten(&iv, h, c, false);
    }
    return iv;
}

static bool iv_finite(double v) { return v > -HUGE_VAL && v < HUGE_VAL; }

/* c = a * b.  When one operand is a nonzero exact point, openness carries
 * exactly; otherwise (both genuine intervals) endpoints are marked CLOSED --
 * a sound over-approximation of attainment (never claims a strict bound it
 * cannot guarantee). */
static bool iv_mul(Iv a, Iv b, Iv* out) {
    bool a_pt = (a.lo == a.hi && !a.lo_open && !a.hi_open);
    bool b_pt = (b.lo == b.hi && !b.lo_open && !b.hi_open);
    if (a_pt || b_pt) {
        Iv pt = a_pt ? a : b, iv = a_pt ? b : a;
        double k = pt.lo;
        if (k == 0.0) { out->lo = out->hi = 0.0; out->lo_open = out->hi_open = false; return true; }
        if (k > 0) { out->lo = k * iv.lo; out->hi = k * iv.hi;
                     out->lo_open = iv.lo_open; out->hi_open = iv.hi_open; }
        else       { out->lo = k * iv.hi; out->hi = k * iv.lo;
                     out->lo_open = iv.hi_open; out->hi_open = iv.lo_open; }
        return true;
    }
    if (!iv_finite(a.lo) || !iv_finite(a.hi) || !iv_finite(b.lo) || !iv_finite(b.hi))
        return false;
    double c[4] = { a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi };
    double mn = c[0], mx = c[0];
    for (int i = 1; i < 4; i++) { if (c[i] < mn) mn = c[i]; if (c[i] > mx) mx = c[i]; }
    out->lo = mn; out->hi = mx; out->lo_open = out->hi_open = false;
    return true;
}

/* c = a + b.  The sum-min is attained only when both mins are attained, so an
 * open endpoint on either side makes the sum endpoint open (sound even under a
 * shared symbol: a > a.lo forces a + b > a.lo + b.lo). */
static bool iv_add(Iv a, Iv b, Iv* out) {
    out->lo = a.lo + b.lo; out->hi = a.hi + b.hi;
    out->lo_open = a.lo_open || b.lo_open;
    out->hi_open = a.hi_open || b.hi_open;
    return true;
}

/* out = 1 / r.  Requires r sign-definite away from 0; 1/x is decreasing on each
 * sign, so 1/[lo,hi] = [1/hi, 1/lo] with openness swapped (0 -> +-inf). */
static bool iv_recip(Iv r, Iv* out) {
    bool pos = (r.lo > 0.0) || (r.lo == 0.0 && r.lo_open);
    bool neg = (r.hi < 0.0) || (r.hi == 0.0 && r.hi_open);
    if (!pos && !neg) return false;                  /* straddles 0 */
    out->lo = (r.hi == 0.0) ? -HUGE_VAL : (iv_finite(r.hi) ? 1.0 / r.hi : 0.0);
    out->hi = (r.lo == 0.0) ?  HUGE_VAL : (iv_finite(r.lo) ? 1.0 / r.lo : 0.0);
    out->lo_open = r.hi_open;
    out->hi_open = r.lo_open;
    return true;
}

/* Sound interval enclosure of `e` over the Assumptions box.  Supports numeric
 * constants, bounded symbols, Plus, Times, and Power[base, +-1]; declines
 * (returns false) on anything else -- which only forfeits a proof, never
 * fabricates one. */
static bool iv_eval(const Expr* e, const Expr* as, Iv* out) {
    double v;
    if (num_dbl(e, &v)) { out->lo = out->hi = v; out->lo_open = out->hi_open = false; return true; }
    if (e->type == EXPR_SYMBOL) {
        Iv b = sym_bound(as, e->data.symbol.name);
        if (!iv_finite(b.lo) && !iv_finite(b.hi)) return false;    /* unbounded */
        *out = b;
        return true;
    }
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (ex->type != EXPR_INTEGER) return false;
        long n = (long)ex->data.integer;
        Iv b;
        if (!iv_eval(base, as, &b)) return false;
        if (n == 1) { *out = b; return true; }
        if (n == -1) return iv_recip(b, out);
        return false;                                              /* |n| >= 2: decline */
    }
    if (head_name_is(e, "Times")) {
        Iv acc = { 1.0, 1.0, false, false };
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Iv f;
            if (!iv_eval(e->data.function.args[i], as, &f)) return false;
            if (!iv_mul(acc, f, &acc)) return false;
        }
        *out = acc; return true;
    }
    if (head_name_is(e, "Plus")) {
        Iv acc = { 0.0, 0.0, false, false };
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Iv t;
            if (!iv_eval(e->data.function.args[i], as, &t)) return false;
            if (!iv_add(acc, t, &acc)) return false;
        }
        *out = acc; return true;
    }
    return false;
}

/* Prove `e > 0` everywhere on the Assumptions box (lower enclosure end > 0). */
static bool iv_prove_pos(const Expr* e, const Expr* as) {
    Iv v;
    return iv_eval(e, as, &v) && ((v.lo > 0.0) || (v.lo == 0.0 && v.lo_open));
}

/* Prove the fugacity gate -1 < e <= 1 on the box: the lower bound is STRICT
 * (exclude the Bose boundary gamma' = -1, whose convergence strip differs),
 * the upper is nonstrict (the Fermi boundary gamma' = 1 keeps Re s > 0). */
static bool iv_prove_fugacity(const Expr* e, const Expr* as) {
    Iv v;
    if (!iv_eval(e, as, &v)) return false;
    bool lower = (v.lo > -1.0) || (v.lo == -1.0 && v.lo_open);
    return lower && (v.hi <= 1.0);
}

/* A_eff > 0: engine proof, else the interval gate (symbolic fugacity). */
static bool prove_pos(const Expr* e, Expr* as) {
    return prove_true(Gt(cp(e), mk_int(0)), as) || iv_prove_pos(e, as);
}

/* 1/(A e^(c x) + gamma), A e^(c0) > 0, c > 0, gamma real with -1 <= gamma' <= 1
 * (gamma' = gamma / A_eff, A_eff = A e^(c0)):  the exponential-geometric kernel of
 * the Bose-Einstein / Fermi-Dirac / general-fugacity integrals.  With u = e^(-c x),
 *   1/(e^(c x)+g') = u/(1+g' u) = (-1/g') sum_{j>=1} (-g')^j u^j    (|g'| e^(-c x) < 1)
 * and Integrate[x^(sv-1) u^j] = Gamma(sv) (j c)^(-sv), so term-by-term
 *   M = (1/A_eff) Gamma(sv) c^(-sv) (-1/g') PolyLog(sv, -g').
 * Special cases:  g' = -1 (Bose) -> Gamma(sv) c^(-sv) Zeta(sv), strip Re sv > 1;
 *                 g' = +1 (Fermi) -> Gamma(sv) c^(-sv) (-PolyLog(sv,-1)) = ... eta(sv),
 *                 strip Re sv > 0 (the -PolyLog form stays finite at sv=1, unlike
 *                 (1-2^(1-sv)) Zeta(sv) which is 0*Infinity there).
 * No interior pole on (0,Inf) for |g'| <= 1 (e^(c x) = -g' has x = ln(-g')/c <= 0);
 * at x=0 the denominator is A_eff(1+g'), zero only for g' = -1, hence the extra
 * Re sv > 1 there.  The keyhole/reflection identity PolyLog(sv,1) = Zeta(sv) is
 * used to land the Bose case on the canonical Zeta head (Simplify does not do this
 * for a symbolic sv). */
static bool rec_expgeom(const Expr* K, const Expr* x, const Expr* sv,
                        Expr* assumptions, Expr** M, Expr** P) {
    if (!head_name_is(K, "Power") || K->data.function.arg_count != 2) return false;
    Expr* base = K->data.function.args[0];
    Expr* ex   = K->data.function.args[1];
    if (ex->type != EXPR_INTEGER || ex->data.integer != -1) return false;  /* K = 1/base */
    if (!head_name_is(base, "Plus")) return false;

    /* Split base into x-free constant gamma and exactly one x-dependent term. */
    Expr* gamma = mk_int(0);
    Expr* expterm = NULL;                              /* borrowed A e^(...) summand */
    for (size_t i = 0; i < base->data.function.arg_count; i++) {
        Expr* ai = base->data.function.args[i];
        if (free_of_x_now(ai, x)) { gamma = Pls(gamma, cp(ai)); continue; }
        if (expterm) { expr_free(gamma); return false; }        /* two x-terms */
        expterm = ai;
    }
    gamma = simp(gamma);
    if (!expterm) { expr_free(gamma); return false; }

    /* expterm = A * Exp[arg], A x-free, arg = c x + c0 linear in x. */
    Expr* facs[32];
    size_t nf = collect_factors((Expr*)expterm, facs, 32, 0);
    Expr* A = mk_int(1);
    Expr* earg = NULL;
    for (size_t i = 0; i < nf; i++) {
        Expr* F = facs[i];
        if (free_of_x_now(F, x)) { A = Tms(A, cp(F)); continue; }
        Expr* a = exp_exponent(F);                              /* Exp[a] / Power[E,a] */
        if (!a || earg) { expr_free(A); if (earg) expr_free(earg); expr_free(gamma); return false; }
        earg = cp(a);
    }
    if (!earg) { expr_free(A); expr_free(gamma); return false; }

    Expr* c1 = dx(earg, x);                                     /* rate c = d(arg)/dx */
    if (!free_of_x_now(c1, x)) { expr_free(c1); expr_free(earg); expr_free(A); expr_free(gamma); return false; }
    Expr* c0 = simp(Pls(cp(earg), Neg(Tms(cp(c1), cp(x)))));    /* constant offset */
    expr_free(earg);
    if (!free_of_x_now(c0, x)) { expr_free(c0); expr_free(c1); expr_free(A); expr_free(gamma); return false; }

    Expr* Aeff = simp(Tms(cp(A), ExpE(cp(c0))));               /* A_eff = A e^(c0) */
    expr_free(A); expr_free(c0);
    if (!prove_pos(c1, assumptions) || !prove_pos(Aeff, assumptions)) {
        expr_free(c1); expr_free(Aeff); expr_free(gamma); return false;
    }

    Expr* gp = simp(Tms(cp(gamma), Pw(cp(Aeff), mk_int(-1))));  /* gamma' = gamma/A_eff */
    expr_free(gamma);
    /* Real gamma' with -1 <= gamma' <= 1 (a complex gamma' leaves the inequalities
     * unevaluated, so prove_true is False -> declined).  A concrete gamma' is
     * discharged by the engine; a SYMBOLIC fugacity is admitted by the interval
     * gate when the Assumptions bound it into (-1, 1] (e.g. 0 < z < 1 proves the
     * general Bose integral 1/(z^-1 e^x - 1) -> Gamma[s] PolyLog[s, z]). */
    if (!prove_true(And2(mk_fn2("GreaterEqual", cp(gp), mk_int(-1)),
                         mk_fn2("LessEqual", cp(gp), mk_int(1))), assumptions) &&
        !iv_prove_fugacity(gp, assumptions)) {
        expr_free(c1); expr_free(Aeff); expr_free(gp); return false;
    }

    Expr* gpp1 = simp(Pls(cp(gp), mk_int(1)));
    bool bose = is_zero_now(gpp1);                              /* gamma' == -1 */
    expr_free(gpp1);

    /* closed = (-1/gamma') PolyLog(sv, -gamma'), canonicalised to Zeta for Bose. */
    Expr* closed = bose
        ? mk_fn1("Zeta", cp(sv))
        : Tms(Tms(mk_int(-1), Pw(cp(gp), mk_int(-1))),
              mk_fn2("PolyLog", cp(sv), Neg(cp(gp))));
    *M = Tms(Pw(cp(Aeff), mk_int(-1)),
             Tms3(Gamma_(cp(sv)), Pw(cp(c1), Neg(cp(sv))), closed));
    /* Strip: Re sv > 0, tightened to Re sv > 1 for the Bose case (denominator
     * zero at x=0).  sv > 1 already implies sv > 0, so the Bose bound is stated
     * alone -- carrying the redundant sv > 0 would survive Simplify (which does
     * not detect the subsumption) and wrongly claim validity on 0 < sv <= 1. */
    Expr* strip = Gt(cp(sv), mk_int(bose ? 1 : 0));
    *P = And2(Gt(cp(c1), mk_int(0)), strip);
    expr_free(c1); expr_free(Aeff); expr_free(gp);
    return true;
}

/* ---- monomial internal substitution K = g(x^k) ------------------------- */

/* Record exponent `e` (an owned Simplify'd copy) into out[] if not present. */
static void expo_add(Expr** out, size_t* n, size_t cap, const Expr* e) {
    Expr* es = simp(cp(e));
    for (size_t i = 0; i < *n; i++)
        if (expr_eq(out[i], es)) { expr_free(es); return; }
    if (*n < cap) out[(*n)++] = es; else expr_free(es);
}
/* Collect the distinct exponents with which x occurs (bare x counts as 1).
 * Returns false if x occurs under an exponent that itself depends on x. */
static bool expo_walk(const Expr* e, const Expr* x, Expr** out, size_t* n, size_t cap) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) {
        if (contains_symbol(e, x)) {
            Expr* one = mk_int(1); expo_add(out, n, cap, one); expr_free(one);
        }
        return true;
    }
    if (e->type != EXPR_FUNCTION) return true;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
            strcmp(base->data.symbol.name, x->data.symbol.name) == 0) {
            if (!free_of_x_now(ex, x)) return false;
            expo_add(out, n, cap, ex);
            return true;                                      /* don't recurse into x^ex */
        }
    }
    if (!expo_walk(e->data.function.head, x, out, n, cap)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!expo_walk(e->data.function.args[i], x, out, n, cap)) return false;
    return true;
}
/* If x occurs in K solely as x^k for one constant k, return an owned copy of k. */
static Expr* single_distinct_x_exponent(const Expr* K, const Expr* x) {
    Expr* out[8]; size_t n = 0;
    bool ok = expo_walk(K, x, out, &n, 8);
    Expr* r = (ok && n == 1) ? out[0] : NULL;
    for (size_t i = 0; i < n; i++) if (out[i] != r) expr_free(out[i]);
    return r;
}

/* Try every base recognizer at spectral variable sv. */
static bool try_recognizers(const Expr* K, const Expr* x, const Expr* sv,
                            Expr* assumptions, Expr** M, Expr** P) {
    if (rec_exp(K, x, sv, M, P))       return true;
    if (rec_gauss(K, x, sv, M, P))     return true;
    if (rec_algebraic(K, x, sv, M, P)) return true;
    if (rec_trig(K, x, sv, M, P))      return true;
    if (rec_bessel(K, x, sv, M, P))    return true;
    if (rec_log(K, x, sv, M, P))       return true;
    if (rec_arctan(K, x, sv, M, P))    return true;
    if (rec_pfq(K, x, sv, M, P))       return true;
    if (rec_polylog(K, x, sv, M, P))   return true;
    if (rec_expgeom(K, x, sv, assumptions, M, P)) return true;
    return false;
}

/* Dispatch a kernel to a recognizer, transparently handling a monomial internal
 * substitution K = g(x^k), k != 1.  With y = x^k,
 *   Integrate[x^(sv-1) g(x^k)] = (1/k) Integrate[y^(sv/k-1) g(y)],
 * so the transform is the base transform evaluated at sv/k, scaled by 1/k, and
 * the convergence strip transfers verbatim (an exact change of variables). */
static bool dispatch_kernel(const Expr* K, const Expr* x, const Expr* sv,
                            Expr* assumptions, Expr** M, Expr** P) {
    *M = *P = NULL;
    if (try_recognizers(K, x, sv, assumptions, M, P)) return true;

    Expr* k = single_distinct_x_exponent(K, x);
    if (!k) return false;
    if (k->type == EXPR_INTEGER && k->data.integer == 1) { expr_free(k); return false; }

    Expr* rule   = mk_fn2("Rule", Pw(cp(x), cp(k)), cp(x));   /* x^k -> x */
    Expr* Klin   = ev2("ReplaceAll", cp(K), rule);
    Expr* sv_eff = simp(Tms(cp(sv), Pw(cp(k), mk_int(-1))));  /* sv / k */
    bool ok = Klin && try_recognizers(Klin, x, sv_eff, assumptions, M, P);
    if (ok) *M = Tms(Pw(cp(k), mk_int(-1)), *M);             /* x (1/k) */
    if (Klin) expr_free(Klin);
    expr_free(sv_eff);
    expr_free(k);
    return ok;
}

/* If F is Log[x] or Power[Log[x], k] (k a positive integer) whose argument is
 * exactly the variable x, return the log weight k (>= 1); else 0.  Only the bare
 * Log[x] qualifies -- Log[1+x] / Log[a x] are ordinary kernels, not weights. */
static long log_x_weight(const Expr* F, const Expr* x) {
    const Expr* logf;
    long k;
    if (head_name_is(F, "Log") && F->data.function.arg_count == 1) { logf = F; k = 1; }
    else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
             head_name_is(F->data.function.args[0], "Log") &&
             F->data.function.args[0]->data.function.arg_count == 1) {
        Expr* e = F->data.function.args[1];
        if (e->type != EXPR_INTEGER || e->data.integer <= 0) return 0;
        logf = F->data.function.args[0];
        k = (long)e->data.integer;
    } else return 0;
    Expr* arg = logf->data.function.args[0];
    return (arg->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
            strcmp(arg->data.symbol.name, x->data.symbol.name) == 0) ? k : 0;
}

/* Decompose `term` into C (x-free), rho (net power of x), a Log[x] weight kw
 * (the total power of a bare Log[x] factor), and the list of x-dependent
 * non-(x-power, non-Log[x]) kernel factors (borrowed pointers into `term`).
 * Returns false only on kernel-list overflow. */
static bool split_term(Expr* term, const Expr* x, Expr** C_out, Expr** rho_out,
                       Expr** kernels, size_t* nk, size_t cap, long* kw_out) {
    Expr* facs[64];
    size_t nf = collect_factors(term, facs, 64, 0);
    Expr* C = mk_int(1);
    Expr* rho = mk_int(0);
    long kw = 0;
    *nk = 0;
    for (size_t i = 0; i < nf; i++) {
        Expr* F = facs[i];
        if (!contains_symbol(F, x)) { C = Tms(C, cp(F)); continue; }
        Expr* e = x_power_exponent(F, x);
        if (e) { rho = Pls(rho, e); continue; }
        long lw = log_x_weight(F, x);
        if (lw > 0) { kw += lw; continue; }
        if (*nk >= cap) { expr_free(C); expr_free(rho); return false; }
        kernels[(*nk)++] = F;
    }
    *C_out = simp(C);
    *rho_out = simp(rho);
    *kw_out = kw;
    return true;
}

/* If `base` is exactly 1 + lam x (constant term 1, a single linear x-term),
 * return lam (owned); else NULL. */
static Expr* parse_unit_plus_linear(const Expr* base, const Expr* x) {
    if (!head_name_is(base, "Plus")) return NULL;
    Expr* c = mk_int(0);
    Expr* xpart = NULL;
    for (size_t i = 0; i < base->data.function.arg_count; i++) {
        Expr* ai = base->data.function.args[i];
        if (free_of_x_now(ai, x)) { c = Pls(c, cp(ai)); continue; }
        if (xpart) { expr_free(c); return NULL; }
        xpart = ai;
    }
    c = simp(c);
    bool one = c->type == EXPR_INTEGER && c->data.integer == 1;
    expr_free(c);
    if (!xpart || !one) return NULL;
    Expr* lam = dx(xpart, x);
    if (!free_of_x_now(lam, x)) { expr_free(lam); return NULL; }
    Expr* cross = simp(Pls(cp(xpart), Neg(Tms(cp(lam), cp(x)))));
    bool bad = !is_zero_now(cross);
    expr_free(cross);
    if (bad) { expr_free(lam); return NULL; }
    return lam;
}

/* Set *base (borrowed) on first use, else require the candidate to equal it. */
static bool base_ok(Expr** base, Expr* cand) {
    if (!*base) { *base = cand; return true; }
    Expr* d = simp(Pls(cp(*base), Neg(cp(cand))));
    bool same = is_zero_now(d);
    expr_free(d);
    return same;
}

/* Tier F -- parametric differentiation.  Kernel factors forming
 *   Log[1+lam x]^n * (1+lam x)^(-w0)   (n >= 1)
 * equal (-1)^n d^n/dw^n[(1+lam x)^(-w)] at w=w0, so
 *   M = (-1)^n d^n/dw^n[ lam^(-sv) Beta(sv, w-sv) ] |_{w=w0},
 * with strip -n < Re sv < w0 (w0 = 0 => upper bound 0).  The Beta factor is
 * built in the manifestly-w=0-regular form lam^(-sv) Gamma(sv) w Gamma(w-sv) /
 * Gamma(1+w) (= lam^(-sv) Gamma(sv) Gamma(w-sv)/Gamma(w)), so the w->w0 limit is
 * a plain substitution even at w0=0.  Consumes nothing; borrows kernels[]. */
static bool rec_logpow(Expr** kernels, size_t nk, const Expr* x, const Expr* sv,
                       Expr** M, Expr** P) {
    long n = 0;
    Expr* base = NULL;                       /* shared 1+lam x (borrowed) */
    Expr* w0 = mk_int(0);                     /* sum of exponent offsets */
    for (size_t i = 0; i < nk; i++) {
        Expr* F = kernels[i];
        if (head_name_is(F, "Log") && F->data.function.arg_count == 1) {
            if (!base_ok(&base, F->data.function.args[0])) { expr_free(w0); return false; }
            n += 1;
        } else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                   head_name_is(F->data.function.args[0], "Log") &&
                   F->data.function.args[0]->data.function.arg_count == 1) {
            Expr* mexp = F->data.function.args[1];
            if (mexp->type != EXPR_INTEGER || mexp->data.integer <= 0) { expr_free(w0); return false; }
            if (!base_ok(&base, F->data.function.args[0]->data.function.args[0])) { expr_free(w0); return false; }
            n += (long)mexp->data.integer;
        } else if (head_name_is(F, "Power") && F->data.function.arg_count == 2 &&
                   head_name_is(F->data.function.args[0], "Plus")) {
            Expr* e = F->data.function.args[1];
            if (!free_of_x_now(e, x)) { expr_free(w0); return false; }
            if (!base_ok(&base, F->data.function.args[0])) { expr_free(w0); return false; }
            w0 = simp(Pls(w0, Neg(cp(e))));   /* (1+lam x)^e contributes -e to w0 */
        } else { expr_free(w0); return false; }
    }
    if (n < 1 || !base) { expr_free(w0); return false; }
    Expr* lam = parse_unit_plus_linear(base, x);
    if (!lam) { expr_free(w0); return false; }

    Expr* w = mk_sym("RmtParamW");
    Expr* T = Tms(Pw(cp(lam), Neg(cp(sv))),
                  Tms(Gamma_(cp(sv)),
                      Tms(cp(w),
                          Tms(Gamma_(Pls(cp(w), Neg(cp(sv)))),
                              Pw(Gamma_(Pls(mk_int(1), cp(w))), mk_int(-1))))));
    Expr* Dn = T;
    for (long i = 0; i < n; i++) Dn = ev2("D", Dn, cp(w));
    Expr* atw0 = ev2("ReplaceAll", Dn, mk_fn2("Rule", cp(w), cp(w0)));
    *M = Tms(mk_int((n % 2) ? -1 : 1), atw0);
    *P = And2(Gt(cp(lam), mk_int(0)),
              And2(Gt(cp(sv), mk_int(-n)), Lt(cp(sv), cp(w0))));
    expr_free(lam); expr_free(w); expr_free(w0);
    return true;
}

/* Build (M, P) for a whole term's kernel list: single-kernel base recognizers
 * (with the monomial wrapper), else the parametric-differentiation family. */
static bool dispatch_term(Expr** kernels, size_t nk, const Expr* x,
                          const Expr* sv, Expr* assumptions, Expr** M, Expr** P) {
    *M = *P = NULL;
    if (nk == 1 && dispatch_kernel(kernels[0], x, sv, assumptions, M, P)) return true;
    if (rec_logpow(kernels, nk, x, sv, M, P)) return true;
    return false;
}

/* True iff `e` is the bare symbol `name`. */
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* Value of Integrate[term, {x, 0, Infinity}], or NULL if the term is out of
 * scope or provably divergent.  On success the convergence strip, reduced under
 * `assumptions`, is returned via *cond_out: NULL if it is discharged (proved
 * True), else the residual predicate to carry as a ConditionalExpression. */
static Expr* mellin_term(Expr* term, const Expr* x, Expr* assumptions,
                         Expr** cond_out) {
    *cond_out = NULL;
    Expr *C = NULL, *rho = NULL;
    Expr* kernels[16]; size_t nk = 0; long kw = 0;
    if (!split_term(term, x, &C, &rho, kernels, &nk, 16, &kw)) return NULL;
    if (nk == 0) { expr_free(C); expr_free(rho); return NULL; }  /* pure (log)power diverges */

    Expr* s = simp(Pls(cp(rho), mk_int(1)));             /* s = rho + 1 */
    expr_free(rho);

    Expr *M = NULL, *P = NULL;
    if (kw == 0) {
        if (!dispatch_term(kernels, nk, x, s, assumptions, &M, &P)) {
            expr_free(C); expr_free(s); return NULL;
        }
    } else {
        /* x^rho Log[x]^kw R(x): since d/ds x^(s-1) = x^(s-1) Log x, a Log[x]^kw
         * weight is the kw-th s-derivative of the base transform M_R(s).  Build
         * M_R in a fresh spectral symbol, differentiate kw times, then set s.
         * Log^kw is dominated by x^(+-eps), so R's OPEN strip carries unchanged. */
        Expr* sv = mk_sym("RmtLogSpectral");
        if (!dispatch_term(kernels, nk, x, sv, assumptions, &M, &P)) {
            expr_free(sv); expr_free(C); expr_free(s); return NULL;
        }
        for (long i = 0; i < kw; i++) M = ev2("D", M, cp(sv));
        M = ev2("ReplaceAll", M, mk_fn2("Rule", cp(sv), cp(s)));
        P = ev2("ReplaceAll", P, mk_fn2("Rule", cp(sv), cp(s)));
        expr_free(sv);
    }
    /* Reduce the strip against the assumptions.  True -> discharged; False ->
     * provably divergent (decline); otherwise carry it as a residual. */
    Expr* Ps = simp2(cp(P), assumptions);
    expr_free(P);
    if (sym_is(Ps, "False")) {
        expr_free(Ps); expr_free(C); expr_free(s); expr_free(M); return NULL;
    }

    Expr* val = simp2(Tms(cp(C), cp(M)), assumptions);
    expr_free(C); expr_free(s); expr_free(M);
    if (!is_finite_value(val)) {                           /* boundary/divergent */
        if (val) expr_free(val);
        expr_free(Ps);
        return NULL;
    }
    if (sym_is(Ps, "True")) expr_free(Ps);                 /* strip discharged */
    else *cond_out = Ps;                                   /* residual condition */
    return val;
}

/* Reduction pre-pass (applied before Expand): rewrite special-function and
 * product kernels that are hypergeometric but not stored as HypergeometricPFQ
 * into pFq form, so rec_pfq (with the monomial wrapper for x^k arguments) closes
 * them.
 *
 * Crucially this runs on the *whole* integrand before Expand distributes, so a
 * cancellation kernel such as the lower incomplete gamma
 *   Gamma[a] - Gamma[a, x] = z^a/a 1F1(a; a+1; -z)
 * is recognised as one convergent kernel -- if Expand split it first, each half
 * (Gamma[a] x^(...) and -Gamma[a,x] x^(...)) would be individually divergent.
 *
 * The Bessel-square rule turns a genuine product of two transcendental kernels
 * (a Mellin convolution) into a single 1F2 via
 *   J_nu(z)^2 = (z/2)^(2nu)/Gamma(nu+1)^2 1F2(nu+1/2; 2nu+1, nu+1; -z^2),
 * so the convolution closes through rec_pfq rather than a Barnes integral.
 * Each rule is a proven identity (all verified numerically in the tests). */
static Expr* reduce_to_hypergeometric(const Expr* f) {
    /* Cheap guard: only pay the parse/replace cost when a reducible head is
     * present. */
    if (!contains_symbol_name(f, "Erf") && !contains_symbol_name(f, "Erfc") &&
        !contains_symbol_name(f, "Gamma") && !contains_symbol_name(f, "BesselJ"))
        return cp(f);
    static const char* RULES =
        "{ Erf[u_] :> (2/Sqrt[Pi]) u HypergeometricPFQ[{1/2}, {3/2}, -u^2], "
        "  Erfc[u_] :> 1 - (2/Sqrt[Pi]) u HypergeometricPFQ[{1/2}, {3/2}, -u^2], "
        "  Gamma[a_] - Gamma[a_, z_] :> z^a/a HypergeometricPFQ[{a}, {a+1}, -z], "
        "  BesselJ[nu_, z_]^2 :> (z/2)^(2 nu)/Gamma[nu+1]^2 "
        "      HypergeometricPFQ[{nu+1/2}, {2 nu+1, nu+1}, -z^2] }";
    Expr* rules = parse_expression(RULES);
    if (!rules) return cp(f);
    return ev2("ReplaceRepeated", cp(f), rules);
}

/* -------------------------------------------------------------------------
 * Frullani integrals -- Integrate[(f(a x) - f(b x))/x, {x, 0, Infinity}]
 *                         = (f(0+) - f(Infinity)) Log(b/a),   a, b > 0.
 *
 * A whole-integrand pre-pass (BEFORE Expand splits the sum): each half of a
 * Frullani integrand is individually divergent, so the term-by-term Mellin path
 * cannot see it -- it must be recognised as a single unit.  The scale ratio
 * rho = b/a is read off the two terms structurally and the pairing is verified
 * by the exact identity (t1 /. x -> rho x) + t2 == 0; the two boundary values
 * are the finite limits of f.  Correct-by-construction (Frullani's theorem).
 * ---------------------------------------------------------------------- */

/* The x-scale k of a Frullani term t = C f(k x): C, k x-free, k != 0, with a
 * single x-dependent factor that is Exp[k x] or a unary g[k x] whose argument is
 * exactly k x (no constant offset).  Returns k (owned) or NULL. */
static Expr* frullani_scale(const Expr* t, const Expr* x) {
    Expr* facs[32];
    size_t nf = collect_factors((Expr*)t, facs, 32, 0);
    Expr* F = NULL;                                    /* the sole x-dependent factor */
    for (size_t i = 0; i < nf; i++) {
        if (free_of_x_now(facs[i], x)) continue;
        if (F) return NULL;                            /* more than one x-factor */
        F = facs[i];
    }
    if (!F) return NULL;
    Expr* arg = exp_exponent(F);                       /* Exp[..] / Power[E,..] */
    if (!arg) {
        if (F->type == EXPR_FUNCTION && F->data.function.arg_count == 1)
            arg = F->data.function.args[0];            /* unary g[arg] */
        else return NULL;
    }
    Expr* k = dx(arg, x);
    if (!free_of_x_now(k, x) || is_zero_now(k)) { expr_free(k); return NULL; }
    Expr* off = simp(Pls(cp(arg), Neg(Tms(cp(k), cp(x)))));   /* arg - k x */
    bool linear = is_zero_now(off);
    expr_free(off);
    if (!linear) { expr_free(k); return NULL; }
    return k;
}

/* Integrate[(f(a x) - f(b x))/x, {x, 0, Infinity}].  Borrows f, x, assumptions;
 * returns the owned value or NULL (out of scope / boundary limits not finite). */
static Expr* frullani_try(Expr* f, const Expr* x, Expr* assumptions) {
    Expr* G = ev1("Expand", Tms(cp((Expr*)x), cp(f)));        /* x f = f(a x) - f(b x) */
    if (!G || !head_name_is(G, "Plus") || G->data.function.arg_count != 2) {
        if (G) { expr_free(G); } return NULL;
    }
    Expr* t1 = G->data.function.args[0];
    Expr* t2 = G->data.function.args[1];
    Expr* k1 = frullani_scale(t1, x);
    Expr* k2 = frullani_scale(t2, x);
    Expr* result = NULL;
    if (k1 && k2) {
        Expr* rho = simp(Tms(cp(k2), Pw(cp(k1), mk_int(-1))));    /* b/a */
        Expr* rm1 = simp(Pls(cp(rho), mk_int(-1)));
        /* rho > 0  <=>  k1, k2 have the same sign.  Prove each scale's sign
         * separately in the ">0" orientation: Simplify discharges a linear sign
         * against the assumptions, but not the ratio b/a > 0 directly (nor a
         * "-a < 0" form, which it normalises to "a > 0" without re-discharging). */
        bool both_pos = prove_true(Gt(cp(k1), mk_int(0)), assumptions) &&
                        prove_true(Gt(cp(k2), mk_int(0)), assumptions);
        bool both_neg = prove_true(Gt(Neg(cp(k1)), mk_int(0)), assumptions) &&
                        prove_true(Gt(Neg(cp(k2)), mk_int(0)), assumptions);
        bool ok = (both_pos || both_neg) && !is_zero_now(rm1);
        expr_free(rm1);
        if (ok) {
            /* Verify the pairing: (t1 /. x -> rho x) + t2 == 0 (either split order
             * gives the same value, so one check suffices). */
            Expr* subst = ev2("ReplaceAll", cp(t1),
                              mk_fn2("Rule", cp((Expr*)x), Tms(cp(rho), cp((Expr*)x))));
            Expr* chk = simp(Pls(subst, cp(t2)));
            bool identity = is_zero_now(chk);
            expr_free(chk);
            if (identity) {
                Expr* f0   = ev2("Limit", cp(t1), mk_fn2("Rule", cp((Expr*)x), mk_int(0)));
                Expr* finf = ev2("Limit", cp(t1),
                                 mk_fn2("Rule", cp((Expr*)x), mk_sym("Infinity")));
                if (is_finite_value(f0) && is_finite_value(finf) &&
                    !contains_symbol(f0, x) && !contains_symbol(finf, x)) {
                    result = simp2(Tms(Pls(cp(f0), Neg(cp(finf))),
                                       mk_fn1("Log", cp(rho))), assumptions);
                    if (!is_finite_value(result) || contains_symbol(result, x)) {
                        if (result) { expr_free(result); } result = NULL;
                    }
                }
                if (f0) expr_free(f0);
                if (finf) expr_free(finf);
            }
        }
        expr_free(rho);
    }
    if (k1) expr_free(k1);
    if (k2) expr_free(k2);
    expr_free(G);
    return result;
}

/* -------------------------------------------------------------------------
 * Sin[r x]^k / x^m on [0, Inf) -- the ssp family (Wang pp. 86-87).
 *
 * Reduce the power with TrigReduce (Sin^k -> a sum of 1 / cos(c x) / sin(c x))
 * and apply the analytically-continued Mellin transforms of cos / sin,
 *   Int_0^Inf x^(s-1) cos(c x) dx = Gamma(s) cos(pi s/2) c^(-s),
 *   Int_0^Inf x^(s-1) sin(c x) dx = Gamma(s) sin(pi s/2) c^(-s),   c > 0,
 * at s = 1 - m.  Gamma(s) has a pole at the non-positive integer s, cancelled
 * by the trig zero exactly in the convergent regime; the reflection forms
 *   Gamma(s) cos(pi s/2) = pi / (2 sin(pi s/2) Gamma(1-s)),
 *   Gamma(s) sin(pi s/2) = pi / (2 cos(pi s/2) Gamma(1-s))
 * make that finite value manifest (and blow up to Infinity for a divergent
 * integrand, which the finiteness gate then rejects).  The scaleless constant
 * term of TrigReduce (c = 0) has no Mellin transform and is dropped.
 * ---------------------------------------------------------------------- */

/* Coefficient c of x in a purely-linear argument c*x (c free of x, and
 * arg - c*x simplifies to 0).  Returns NULL if arg is not of that form. */
static Expr* linear_coeff(const Expr* arg, const Expr* x) {
    Expr* c = simp(mk_fn2("Divide", cp((Expr*)arg), cp((Expr*)x)));
    if (!c || contains_symbol(c, x)) { if (c) expr_free(c); return NULL; }
    Expr* resid = simp(mk_fn2("Subtract", cp((Expr*)arg), Tms(cp(c), cp((Expr*)x))));
    bool zero = resid && ((resid->type == EXPR_INTEGER && resid->data.integer == 0) ||
                          (resid->type == EXPR_REAL && resid->data.real == 0.0));
    if (resid) expr_free(resid);
    if (!zero) { expr_free(c); return NULL; }
    return c;
}

/* Mellin contribution of one TrigReduce term `coeff * trig(c x)` at s = 1 - m,
 * or the constant/zero term (returns 0), or NULL if the term is not of the
 * expected form.  `sExpr` is the (owned-by-caller) exponent 1 - m. */
static Expr* sinpow_term(Expr* term, const Expr* x, const Expr* sExpr) {
    /* Split the term into an x-free coefficient and the trig kernel. */
    size_t nf; Expr** fs; Expr* one[1];
    if (head_name_is(term, "Times")) { nf = term->data.function.arg_count; fs = term->data.function.args; }
    else { nf = 1; one[0] = term; fs = one; }

    Expr* coeff = mk_int(1);
    Expr* kernel = NULL;                 /* the Cos[..]/Sin[..] factor */
    for (size_t i = 0; i < nf; i++) {
        Expr* g = fs[i];
        if (!contains_symbol(g, x)) { coeff = Tms(coeff, cp(g)); continue; }
        if ((head_name_is(g, "Cos") || head_name_is(g, "Sin")) &&
            g->data.function.arg_count == 1 && !kernel) {
            kernel = g; continue;
        }
        expr_free(coeff); return NULL;   /* x-dependent non-trig factor */
    }
    if (!kernel) { expr_free(coeff); return mk_int(0); }   /* scaleless constant -> drop */

    Expr* c = linear_coeff(kernel->data.function.args[0], x);
    if (!c) { expr_free(coeff); return NULL; }

    /* Reflection-form bracket, finite in the convergent regime. */
    bool is_cos = head_name_is(kernel, "Cos");
    Expr* halfPis = Tms(Rat(1, 2), Tms(mk_sym("Pi"), cp(sExpr)));   /* pi s / 2 */
    Expr* trig_denom = is_cos ? mk_fn1("Sin", cp(halfPis)) : mk_fn1("Cos", cp(halfPis));
    expr_free(halfPis);
    Expr* bracket = mk_fn2("Divide", mk_sym("Pi"),
        Tms(mk_int(2), Tms(trig_denom, Gamma_(mk_fn2("Subtract", mk_int(1), cp(sExpr))))));
    /* coeff * bracket * c^(-s). */
    Expr* term_val = Tms3(coeff, bracket, Pw(c, Neg(cp(sExpr))));
    return term_val;
}

Expr* integrate_sinpowmono_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!is_zero_expr(a) || !is_pos_inf(b)) return NULL;
    if (!contains_symbol(f, x)) return NULL;

    /* Parse f = C * Sin[r x]^k * x^(-m), C free of x, r > 0, k>=1, m>=1 ints. */
    size_t nf; Expr** fs; Expr* one[1];
    if (head_name_is(f, "Times")) { nf = f->data.function.arg_count; fs = f->data.function.args; }
    else { nf = 1; one[0] = f; fs = one; }

    Expr* Cc = mk_int(1);
    Expr* rarg = NULL;         /* Sin argument */
    long k = 0, m = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* g = fs[i];
        if (!contains_symbol(g, x)) { Cc = Tms(Cc, cp(g)); continue; }
        /* Sin[arg] or Sin[arg]^k */
        if (head_name_is(g, "Sin") && g->data.function.arg_count == 1 && !rarg) {
            rarg = g->data.function.args[0]; k = 1; continue;
        }
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
            Expr* base = g->data.function.args[0];
            Expr* e    = g->data.function.args[1];
            if (head_name_is(base, "Sin") && base->data.function.arg_count == 1 &&
                e->type == EXPR_INTEGER && e->data.integer >= 1 && !rarg) {
                rarg = base->data.function.args[0]; k = e->data.integer; continue;
            }
            if (base->type == EXPR_SYMBOL && base->data.symbol.name == x->data.symbol.name &&
                e->type == EXPR_INTEGER && e->data.integer < 0 && m == 0) {
                m = -e->data.integer; continue;
            }
        }
        ok = false;
    }
    if (!ok || !rarg || k < 1 || m < 1) { expr_free(Cc); return NULL; }

    Expr* r = linear_coeff(rarg, x);
    if (!r) { expr_free(Cc); return NULL; }
    /* Need r > 0 for c^(-s) to be real/well-defined. */
    if (!prove_pos(r, assumptions)) { expr_free(Cc); expr_free(r); return NULL; }
    expr_free(r);

    Expr* sExpr = mk_int(1 - m);                       /* s = 1 - m */
    /* G = Expand[TrigReduce[Sin[r x]^k]] -- TrigReduce returns a *factored*
     * form like (1/2)(1 - Cos[2x]); Expand distributes it into a sum of
     * coefficient*harmonic terms the loop below can split. */
    Expr* sink = Pw(mk_fn1("Sin", cp(rarg)), mk_int(k));
    Expr* G = ev1("Expand", mk_fn1("TrigReduce", sink));
    if (!G) { expr_free(Cc); expr_free(sExpr); return NULL; }

    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(G, "Plus")) { nt = G->data.function.arg_count; terms = G->data.function.args; }
    else { nt = 1; single[0] = G; terms = single; }

    Expr* total = mk_int(0);
    bool good = true;
    for (size_t t = 0; t < nt && good; t++) {
        Expr* v = sinpow_term(terms[t], x, sExpr);
        if (!v) { good = false; break; }
        total = Pls(total, v);
    }
    expr_free(G); expr_free(sExpr);
    if (!good) { expr_free(total); expr_free(Cc); return NULL; }

    Expr* value = fullsimp2(Tms(Cc, total), assumptions);
    if (!value || !is_finite_value(value) || contains_symbol(value, x)) {
        if (value) expr_free(value);
        return NULL;
    }
    return value;
}

/* -------------------------------------------------------------------------
 * R(x) Log[x]^n on [0, Inf) with R a proper rational whose poles are all on the
 * negative real axis (real linear factors) -- the log*rat / logquad0 family for
 * multi-pole R that the single-generator (a + x^m)^(-p) Mellin path does not
 * cover.  Uses the Mellin transform of a shifted pole,
 *   Int_0^Inf x^(s-1)/(x+a)^p dx = a^(s-p) B(s, p-s)   (a > 0, 0 < Re s < p),
 * summed over the partial fractions of R to give M(s); then
 *   Int_0^Inf R(x) Log[x]^n dx = d^n/ds^n M(s) |_{s=1}.
 * Each per-term Mellin converges on 0 < Re s < p even where the log-weighted
 * integral's partial fractions individually diverge, so M(s) is analytic and
 * the derivative-at-1 is the finite answer.
 * ---------------------------------------------------------------------- */

/* Parse one Apart term `c * (x + a)^(-p)` (a > 0 real, p >= 1 int, c free of x).
 * On success returns the Mellin contribution c * a^(s-p) * Beta[s, p-s] (with
 * the caller's symbolic `s`), else NULL. */
static Expr* ratlog_term(Expr* term, const Expr* x, const Expr* s, Expr* as) {
    size_t nf; Expr** fs; Expr* one[1];
    if (head_name_is(term, "Times")) { nf = term->data.function.arg_count; fs = term->data.function.args; }
    else { nf = 1; one[0] = term; fs = one; }

    Expr* c = mk_int(1);
    Expr* a = NULL; long p = 0;
    for (size_t i = 0; i < nf; i++) {
        Expr* g = fs[i];
        if (!contains_symbol(g, x)) { c = Tms(c, cp(g)); continue; }
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2 &&
            g->data.function.args[1]->type == EXPR_INTEGER &&
            g->data.function.args[1]->data.integer < 0 && !a) {
            Expr* base = g->data.function.args[0];
            /* base must be x + a with unit x-coefficient. */
            Expr* slope = simp(mk_fn2("D", cp(base), cp((Expr*)x)));
            bool unit = slope && slope->type == EXPR_INTEGER && slope->data.integer == 1;
            if (slope) expr_free(slope);
            if (!unit) { expr_free(c); return NULL; }
            Expr* aval = simp(mk_fn2("Subtract", cp(base), cp((Expr*)x)));   /* a = base - x */
            if (!aval || contains_symbol(aval, x)) { if (aval) expr_free(aval); expr_free(c); return NULL; }
            a = aval; p = -g->data.function.args[1]->data.integer;
            continue;
        }
        expr_free(c); return NULL;                 /* polynomial part or unexpected factor */
    }
    if (!a) { expr_free(c); return NULL; }          /* pole-free term -> divergent */
    if (!prove_pos(a, as)) { expr_free(c); expr_free(a); return NULL; }

    /* c * a^(s-p) * Beta[s, p-s].  For a simple pole (p = 1) Beta[s, 1-s] is
     * singular at s = 1 and the Limit engine cannot resolve the derivative
     * there, so use the reflection form Beta[s, 1-s] = Pi/Sin[Pi s], which it
     * handles; for p >= 2 the Beta is finite at s = 1 and needs no rewrite. */
    Expr* sp = mk_int(p);
    Expr* kernel = (p == 1)
        ? mk_fn2("Divide", mk_sym("Pi"), mk_fn1("Sin", Tms(mk_sym("Pi"), cp(s))))
        : Beta_(cp(s), mk_fn2("Subtract", cp(sp), cp(s)));
    Expr* contrib = Tms3(c, Pw(a, mk_fn2("Subtract", cp(s), cp(sp))), kernel);
    expr_free(sp);
    return contrib;
}

Expr* integrate_ratlogpow_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!is_zero_expr(a) || !is_pos_inf(b)) return NULL;
    if (!contains_symbol(f, x)) return NULL;

    /* Separate the Log[x]^n weight from the rational part R. */
    size_t nf; Expr** fs; Expr* one[1];
    if (head_name_is(f, "Times")) { nf = f->data.function.arg_count; fs = f->data.function.args; }
    else { nf = 1; one[0] = f; fs = one; }

    long n = 0;
    Expr* R = mk_int(1);
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* g = fs[i];
        if (head_name_is(g, "Log") && g->data.function.arg_count == 1 &&
            g->data.function.args[0]->type == EXPR_SYMBOL &&
            g->data.function.args[0]->data.symbol.name == x->data.symbol.name && n == 0) {
            n = 1; continue;
        }
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
            Expr* base = g->data.function.args[0];
            Expr* e    = g->data.function.args[1];
            if (head_name_is(base, "Log") && base->data.function.arg_count == 1 &&
                base->data.function.args[0]->type == EXPR_SYMBOL &&
                base->data.function.args[0]->data.symbol.name == x->data.symbol.name &&
                e->type == EXPR_INTEGER && e->data.integer >= 1 && n == 0) {
                n = e->data.integer; continue;
            }
        }
        if (contains_symbol_name(g, "Log")) { ok = false; break; }   /* other Log -> out of scope */
        R = Tms(R, cp(g));
    }
    if (!ok || n < 1) { expr_free(R); return NULL; }

    /* Partial-fraction R over x (Together first so a product of factors becomes
     * a single fraction Apart can decompose). */
    Expr* Rtog = ev1("Together", R);          /* consumes R */
    Expr* Ap = Rtog ? ev2("Apart", Rtog, cp((Expr*)x)) : NULL;
    if (!Ap) return NULL;

    Expr* s = mk_sym("Integrate`RL`s");
    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(Ap, "Plus")) { nt = Ap->data.function.arg_count; terms = Ap->data.function.args; }
    else { nt = 1; single[0] = Ap; terms = single; }

    Expr* M = mk_int(0);
    bool good = true;
    for (size_t t = 0; t < nt && good; t++) {
        Expr* v = ratlog_term(terms[t], x, s, assumptions);
        if (!v) { good = false; break; }
        M = Pls(M, v);
    }
    expr_free(Ap);
    if (!good) { expr_free(M); expr_free(s); return NULL; }

    /* Int R Log^n = M^(n)(1) = n! * [coefficient of (s-1)^n in M(s)].  The
     * SeriesCoefficient handles the removable singularity of M at s = 1 (where
     * the individual reflection factors Pi/Sin[Pi s] blow up but the residue sum
     * stays analytic) more robustly than a repeated-derivative Limit. */
    Expr* spec = expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ cp(s), mk_int(1), mk_int(n) }, 3);
    Expr* coef = eval_take(mk_fn2("SeriesCoefficient", M, spec));   /* consumes M */
    Expr* ats1 = Tms(mk_fn1("Factorial", mk_int(n)), coef);
    expr_free(s);

    Expr* value = fullsimp2(ats1, assumptions);
    if (!value || !is_finite_value(value) || contains_symbol_name(value, "Integrate`RL`s")) {
        if (value) expr_free(value);
        return NULL;
    }
    return value;
}

/* -------------------------------------------------------------------------
 * Public entry point.
 * ---------------------------------------------------------------------- */
Expr* integrate_ramanujan_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!is_zero_expr(a) || !is_pos_inf(b)) return NULL;   /* half-line only */
    if (!contains_symbol(f, x)) return NULL;

    /* Frullani pre-pass: (f(a x) - f(b x))/x must match as a whole (each half is
     * individually divergent, so the term-by-term Mellin path below cannot). */
    Expr* fru = frullani_try(f, x, assumptions);
    if (fru) return fru;

    /* Sin[r x]^k / x^m pre-pass (the ssp family): the constant TrigReduce term is
     * scaleless and the individual cos/sin terms straddle a Gamma pole, so this
     * too must be summed as a whole rather than through the term loop below. */
    Expr* ssp = integrate_sinpowmono_try(f, x, a, b, assumptions);
    if (ssp) return ssp;

    /* R(x) Log[x]^n pre-pass (log*rat): the partial fractions of R Log^n are
     * individually divergent, so the Mellin-derivative M^(n)(1) is summed as a
     * whole rather than through the term loop. */
    Expr* rlp = integrate_ratlogpow_try(f, x, a, b, assumptions);
    if (rlp) return rlp;

    Expr* fr = reduce_to_hypergeometric(f);
    if (!fr) return NULL;
    Expr* g = ev1("Expand", fr);
    if (!g) return NULL;

    size_t nt; Expr** terms; Expr* single[1];
    if (head_name_is(g, "Plus")) {
        nt = g->data.function.arg_count;
        terms = g->data.function.args;
    } else { nt = 1; single[0] = g; terms = single; }

    Expr* total = mk_int(0);
    Expr* cond  = NULL;                                     /* And of residual strips */
    bool ok = true;
    for (size_t t = 0; t < nt && ok; t++) {
        Expr* c = NULL;
        Expr* v = mellin_term(terms[t], x, assumptions, &c);
        if (!v) { ok = false; if (c) expr_free(c); break; }
        total = Pls(total, v);
        if (c) cond = cond ? And2(cond, c) : c;
    }
    expr_free(g);
    if (!ok) { expr_free(total); if (cond) expr_free(cond); return NULL; }

    Expr* res = fullsimp2(total, assumptions);
    if (!res || !is_finite_value(res) || contains_symbol(res, x)) {
        if (res) expr_free(res);
        if (cond) expr_free(cond);
        return NULL;
    }
    if (!cond) return res;                                  /* strip fully discharged */

    /* Carry the residual strip as a ConditionalExpression (it collapses to the
     * bare value when the assumptions later prove it, to Undefined if refuted). */
    Expr* cs = simp2(cond, assumptions);
    if (sym_is(cs, "False")) { expr_free(cs); expr_free(res); return NULL; }
    if (sym_is(cs, "True"))  { expr_free(cs); return res; }
    return eval_take(mk_fn2("ConditionalExpression", res, cs));
}

/* -------------------------------------------------------------------------
 * `Integrate`RamanujanMasterTheorem[f, {x,0,Infinity}]` (optionally
 * Assumptions -> ...) builtin.
 * ---------------------------------------------------------------------- */
Expr* builtin_integrate_ramanujan(Expr* res) {
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
        return NULL;                                       /* unknown trailing arg */
    }
    return integrate_ramanujan_try(f, x, a, b, assumptions);
}

/* Shared spec parse for the half-line builtins: head[f, {x,0,Inf}, Assumptions].
 * Returns false on malformed input; on success sets the borrowed pointers. */
static bool half_line_parse(Expr* res, Expr** f, Expr** x, Expr** a, Expr** b,
                            Expr** assumptions) {
    if (res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return false;
    *f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return false;
    *x = spec->data.function.args[0];
    *a = spec->data.function.args[1];
    *b = spec->data.function.args[2];
    if ((*x)->type != EXPR_SYMBOL) return false;
    *assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol.name, "Assumptions") == 0) {
                *assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return false;
    }
    return true;
}

Expr* builtin_integrate_sinpowmono(Expr* res) {
    Expr *f, *x, *a, *b, *as;
    if (!half_line_parse(res, &f, &x, &a, &b, &as)) return NULL;
    return integrate_sinpowmono_try(f, x, a, b, as);
}

Expr* builtin_integrate_oscpower(Expr* res) {
    Expr *f, *x, *a, *b, *as;
    if (!half_line_parse(res, &f, &x, &a, &b, &as)) return NULL;
    /* Cos/Sin[b x^n] is closed by the Mellin engine's monomial-substitution
     * path; this builtin is the named entry point / Method for that family. */
    return integrate_ramanujan_try(f, x, a, b, as);
}

Expr* builtin_integrate_ratlogpow(Expr* res) {
    Expr *f, *x, *a, *b, *as;
    if (!half_line_parse(res, &f, &x, &a, &b, &as)) return NULL;
    return integrate_ratlogpow_try(f, x, a, b, as);
}

void integrate_ramanujan_init(void) {
    symtab_add_builtin("Integrate`SinPowerMonomial", builtin_integrate_sinpowmono);
    symtab_get_def("Integrate`SinPowerMonomial")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`SinPowerMonomial",
        "Integrate`SinPowerMonomial[f, {x, 0, Infinity}] evaluates Sin[r x]^k / "
        "x^m over the half line (the ssp family): TrigReduce lowers the power to "
        "a sum of cos/sin harmonics and the analytically-continued Mellin "
        "transforms of cos and sin (reflection form, so the Gamma pole cancels "
        "against the trig zero) are summed.  Returns unevaluated when the "
        "integrand is not of this form or the result is not finite (divergent).");

    symtab_add_builtin("Integrate`RationalLog", builtin_integrate_ratlogpow);
    symtab_get_def("Integrate`RationalLog")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`RationalLog",
        "Integrate`RationalLog[f, {x, 0, Infinity}] evaluates R(x) Log[x]^n over "
        "the half line for a proper rational R whose poles all lie on the "
        "negative real axis: the Mellin transform M(s) = Integrate[x^(s-1) R(x)] "
        "is summed over the partial fractions of R (each pole x = -a contributing "
        "a^(s-p) Beta[s, p-s]) and the integral is the n-th s-derivative of M at "
        "s = 1.  Returns unevaluated when R has a positive-real / complex pole, is "
        "improper, or the result is not finite.");

    symtab_add_builtin("Integrate`OscillatoryPower", builtin_integrate_oscpower);
    symtab_get_def("Integrate`OscillatoryPower")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`OscillatoryPower",
        "Integrate`OscillatoryPower[f, {x, 0, Infinity}] evaluates the Fresnel-"
        "type oscillatory integrals Cos[b x^n] and Sin[b x^n] over the half line "
        "to Gamma[1/n]/(n b^(1/n)) {Cos,Sin}[Pi/(2 n)] (b > 0, n > 1).  It is the "
        "named entry point for this family, which the Mellin / Ramanujan engine "
        "closes via its monomial-substitution path.  Returns unevaluated "
        "otherwise.");

    symtab_add_builtin("Integrate`RamanujanMasterTheorem",
                       builtin_integrate_ramanujan);
    symtab_get_def("Integrate`RamanujanMasterTheorem")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`RamanujanMasterTheorem",
        "Integrate`RamanujanMasterTheorem[f, {x, 0, Infinity}] evaluates a "
        "half-line definite integral by the Mellin-transform / Ramanujan Master "
        "Theorem method: the integrand is split term-by-term into a power of x "
        "times a base kernel whose Mellin transform is known in closed form -- "
        "exponential, Gaussian, algebraic binomial, Cos, Sin, ArcTan, Log[1+x], "
        "BesselJ, HypergeometricPFQ, PolyLog, or the exponential-geometric kernel "
        "1/(e^(c x)+g) -> Gamma PolyLog (Bose-Einstein 1/(e^x-1) -> Gamma Zeta, "
        "Fermi-Dirac 1/(e^x+1) -> Gamma eta) -- with a monomial substitution for "
        "x^k arguments (e.g. Sin[Sqrt[x]]), a Log[x]^k weight handled as the k-th "
        "s-derivative of the base transform, and parametric differentiation for "
        "Log[1+x]-power kernels.  A whole-integrand Frullani pre-pass evaluates "
        "(f(a x)-f(b x))/x -> (f(0)-f(Infinity)) Log(b/a).  Erf, the lower "
        "incomplete gamma, and BesselJ^2 are reduced to pFq form first.  Each "
        "transform is gated on its convergence strip; when Assumptions do not "
        "prove the strip the result is a ConditionalExpression stating it.  "
        "Returns unevaluated when the integrand is out of scope or provably "
        "divergent.");
}
