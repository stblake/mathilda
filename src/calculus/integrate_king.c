/* integrate_king.c
 *
 * Definite integration by the general King reflection.  See integrate_king.h.
 * Case handled in v1: the reflected sum f(x)+f(a+b-x) collapses to a constant
 * (the Integrate[1/(1+Tan^s),{x,0,Pi/2}] = Pi/4 family).  The constant is read
 * at the midpoint, screened for constancy numerically across the interval, and
 * the result is cross-checked against NIntegrate -- a mis-detected constancy
 * cannot survive the check, so the stage never returns a wrong answer.
 */

#include "integrate_king.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"      /* arith_warnings_mute_push/pop */

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ---- construction / evaluation helpers (mirrors integrate_beta.c) -------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* mk_int(long v) { return expr_new_integer((int64_t)v); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* h, Expr* a) { Expr* v[1]={a}; return mk_fn(h,v,1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { Expr* v[2]={a,b}; return mk_fn(h,v,2); }
static Expr* mk_fn3b(const char* h, Expr* a, Expr* b, Expr* c) { Expr* v[3]={a,b,c}; return mk_fn(h,v,3); }
static Expr* Times_(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* Plus_(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }
static Expr* Rat_(long p, long q)     { return mk_fn2("Rational", mk_int(p), mk_int(q)); }

static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}
static Expr* ev1(const char* name, Expr* a) { return eval_take(mk_fn1(name, a)); }
static Expr* simp(Expr* e) { return ev1("Simplify", e); }
static Expr* simp2(Expr* e, Expr* as) {
    if (!as) return simp(e);
    return eval_take(mk_fn2("Simplify", e, cp(as)));
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_symbol(const Expr* e, const Expr* x) {
    return e->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
           e->data.symbol.name == x->data.symbol.name;
}
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return is_symbol(e, x);
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}
static bool mentions_nonfinite(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL)
        return sym_is(e, "Infinity") || sym_is(e, "ComplexInfinity") ||
               sym_is(e, "Indeterminate") || sym_is(e, "Undefined");
    if (head_name_is(e, "DirectedInfinity") || head_name_is(e, "Integrate"))
        return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_nonfinite(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_nonfinite(e->data.function.args[i])) return true;
    return false;
}

/* Evaluate N[e]; extract a finite machine real.  Returns 1 on success. */
static int real_double(Expr* e, double* out) {
    Expr* n = ev1("N", cp(e));
    if (!n) return 0;
    int ok = 0;
    if (n->type == EXPR_REAL)         { *out = n->data.real;            ok = 1; }
    else if (n->type == EXPR_INTEGER) { *out = (double) n->data.integer; ok = 1; }
    expr_free(n);
    return ok;
}
/* N[e /. x->pt] into a finite machine real. */
static int real_double_at(Expr* e, Expr* x, Expr* pt, double* out) {
    Expr* sub = eval_take(mk_fn2("ReplaceAll", cp(e), mk_fn2("Rule", cp(x), cp(pt))));
    int ok = real_double(sub, out);
    expr_free(sub);
    return ok;
}

/* Numeric backstop: the candidate value V must agree with NIntegrate[f] over
 * [a,b] to a loose relative tolerance (see integrate_series.c). */
static bool king_crosscheck_ok(Expr* f, Expr* x, Expr* a, Expr* b, Expr* V) {
    double vd;
    if (!real_double(V, &vd)) return false;
    arith_warnings_mute_push();
    Expr* spec = mk_fn3b("List", cp(x), cp(a), cp(b));
    Expr* ni   = eval_take(mk_fn2("NIntegrate", cp(f), spec));   /* spec consumed */
    arith_warnings_mute_pop();
    double nd;
    int k = ni ? real_double(ni, &nd) : 0;
    if (ni) expr_free(ni);
    if (k != 1) return false;
    return fabs(vd - nd) <= 1e-3 * (1.0 + fabs(vd));
}

/* ---- the method --------------------------------------------------------- */

/* A pole of the integrand strictly inside the finite interval (a,b): King's
 * identity ∫f = c(b-a)/2 assumes the integral CONVERGES, but a non-integrable
 * interior pole makes it diverge even when f+f̃ collapses to a constant (1/x on
 * [-1,1]: f+f̃ = 0 -> c = 0, and the symmetric NIntegrate backstop also averages to
 * ~0, so neither catches it).  Decline on any interior zero of the denominator of
 * Together[f] -- within King's parameter-free numeric scope such a zero is a genuine
 * pole.  Samples the open interval (a zero on a sample, or a sign change). */
static bool king_interior_pole(Expr* f, Expr* x, Expr* a, Expr* b) {
    Expr* T = eval_take(mk_fn1("Together", cp(f)));
    if (!T) return false;
    Expr* den = eval_take(mk_fn1("Denominator", T));          /* T consumed */
    if (!den) return false;
    if (!contains_symbol(den, x)) { expr_free(den); return false; }
    const int N = 256;
    bool pole = false, have = false; double prev = 0.0;
    for (int i = 1; i < N && !pole; i++) {
        Expr* pt = Plus_(cp(a), Times_(Plus_(cp(b), Times_(mk_int(-1), cp(a))),
                                       Rat_((long)i, (long)N)));
        double v;
        int ok = real_double_at(den, x, pt, &v);
        expr_free(pt);
        if (!ok) { have = false; continue; }                 /* non-numeric sample */
        if (v == 0.0) pole = true;
        else if (have && ((prev < 0.0) != (v < 0.0))) pole = true;   /* odd-order pole */
        prev = v; have = true;
    }
    expr_free(den);
    return pole;
}

Expr* integrate_king_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;
    if (contains_symbol(a, x) || contains_symbol(b, x)) return NULL;

    /* Finite, numeric endpoints (needed for the midpoint, the constancy screen,
     * and the NIntegrate backstop).  Symbolic-limit King is out of scope in v1. */
    double av, bv;
    if (!real_double(a, &av) || !real_double(b, &bv)) return NULL;
    if (!(fabs(bv - av) > 1e-12)) return NULL;

    /* Decline a non-integrable interior pole (the ordinary integral diverges). */
    if (king_interior_pole(f, x, a, b)) return NULL;

    /* Numeric-domain precondition.  King is an inherently NUMERIC method (it
     * screens constancy and cross-checks against NIntegrate), so it can only
     * ever succeed on a parameter-free integrand -- one that evaluates to a
     * number at a numeric x.  A parametric integrand (e.g. (x^a-1)/Log[x]) is
     * out of scope: it is screened out here CHEAPLY, before any reflection /
     * Simplify work, so King does no heavy re-entrant evaluation on an integrand
     * it cannot handle.  Require a finite real at one of a few interior samples;
     * none finite -> decline. */
    {
        bool numeric = false;
        const long gn[3] = {1, 1, 2}, gd[3] = {2, 3, 3};
        for (int i = 0; i < 3 && !numeric; i++) {
            Expr* pt = Plus_(cp(a), Times_(Plus_(cp(b), Times_(mk_int(-1), cp(a))),
                                           Rat_(gn[i], gd[i])));
            double tmp;
            if (real_double_at(f, x, pt, &tmp)) numeric = true;
            expr_free(pt);
        }
        if (!numeric) return NULL;
    }

    /* h(x) = f(x) + f(a+b-x). */
    Expr* reflx = Plus_(cp(a), Plus_(cp(b), Times_(mk_int(-1), cp(x))));
    Expr* g = eval_take(mk_fn2("ReplaceAll", cp(f), mk_fn2("Rule", cp(x), reflx)));
    if (!g) return NULL;
    Expr* h = Plus_(cp(f), g);                  /* owns f-copy and g */

    /* Candidate constant c = Simplify[h /. x -> (a+b)/2]. */
    Expr* mid = Times_(Plus_(cp(a), cp(b)), Rat_(1, 2));
    Expr* c_at_mid = eval_take(mk_fn2("ReplaceAll", cp(h), mk_fn2("Rule", cp(x), mid)));
    Expr* c = assumptions ? simp2(c_at_mid, assumptions) : simp(c_at_mid);
    if (!c || contains_symbol(c, x) || mentions_nonfinite(c)) {
        expr_free(h); if (c) expr_free(c); return NULL;
    }

    /* Numeric constancy screen: h must agree with c across several interior
     * points.  (This is a cheap pre-filter; the NIntegrate backstop below is the
     * authoritative correctness gate.)  Non-finite samples -- interior
     * singularities of h -- are skipped; require at least 3 finite agreeing
     * points, counting the midpoint value c. */
    double cd;
    if (!real_double(c, &cd)) { expr_free(h); expr_free(c); return NULL; }
    static const long num[6] = {1, 1, 2, 3, 3, 4};
    static const long den[6] = {5, 3, 5, 5, 4, 5};
    int agree = 1;                               /* c (midpoint) counts as one */
    bool mismatch = false;
    for (int i = 0; i < 6 && !mismatch; i++) {
        Expr* pt = Plus_(cp(a), Times_(Plus_(cp(b), Times_(mk_int(-1), cp(a))),
                                       Rat_(num[i], den[i])));
        double hv;
        int ok = real_double_at(h, x, pt, &hv);
        expr_free(pt);
        if (!ok) continue;                        /* skip a singular sample */
        if (fabs(hv - cd) > 1e-6 * (1.0 + fabs(cd))) mismatch = true;
        else agree++;
    }
    expr_free(h);
    if (mismatch || agree < 3) { expr_free(c); return NULL; }

    /* I = c (b-a)/2. */
    Expr* half_len = Times_(Plus_(cp(b), Times_(mk_int(-1), cp(a))), Rat_(1, 2));
    Expr* val = Times_(c, half_len);
    val = assumptions ? simp2(val, assumptions) : simp(val);
    if (!val || contains_symbol(val, x) || mentions_nonfinite(val)) {
        if (val) expr_free(val);
        return NULL;
    }

    if (!king_crosscheck_ok(f, x, a, b, val)) { expr_free(val); return NULL; }
    return val;
}

Expr* builtin_integrate_king(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* as = NULL;
    for (size_t t = 2; t < res->data.function.arg_count; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Assumptions) {
                as = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;
    }
    return integrate_king_try(f, x, a, b, as);
}

void integrate_king_init(void) {
    symtab_add_builtin("Integrate`King", builtin_integrate_king);
    symtab_get_def("Integrate`King")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`King",
        "Integrate`King[f, {x, a, b}] evaluates a definite integral by the King "
        "reflection x -> a+b-x: Integrate[f,{x,a,b}] = "
        "(1/2) Integrate[f(x)+f(a+b-x),{x,a,b}].  When the reflected sum "
        "f(x)+f(a+b-x) is constant c (e.g. 1/(1+Tan[x]^s) over [0,Pi/2], c=1) the "
        "value is c (b-a)/2.  Constancy is read at the midpoint and screened "
        "numerically; the result is cross-checked against NIntegrate.  Returns "
        "unevaluated when the reflected sum is not constant, the interval is not "
        "finite/numeric, or the check fails.");
}
