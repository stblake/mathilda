/* integrate_feynparam.c
 *
 * Definite integration of the parameter / Frullani-logarithm family:
 *   Integrate[(Sum c_i x^(e_i))/Log[x], {x,0,1}] = Sum c_i Log[e_i+1]
 * (Sum c_i = 0, each e_i+1 > 0).  See integrate_feynparam.h.  The inner
 * x-power integral is taken directly as 1/(e+1), so the method stays clear of
 * the Beta/Gamma form the simplifier cannot reduce and does no long re-entrant
 * Integrate.  Correct-by-construction + an NIntegrate backstop.
 */

#include "integrate_feynparam.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ---- helpers (mirror integrate_series.c / integrate_king.c) ------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* mk_int(long v) { return expr_new_integer((int64_t)v); }
static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* h, Expr* a) { Expr* v[1]={a}; return mk_fn(h,v,1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { Expr* v[2]={a,b}; return mk_fn(h,v,2); }
static Expr* mk_fn3b(const char* h, Expr* a, Expr* b, Expr* c) { Expr* v[3]={a,b,c}; return mk_fn(h,v,3); }
static Expr* Pow_(Expr* a, Expr* b)   { return mk_fn2("Power", a, b); }
static Expr* Times_(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* Plus_(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }

static Expr* eval_take(Expr* call) { Expr* r = evaluate(call); expr_free(call); return r; }
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
/* True iff `e` contains the subexpression Log[x] (arg exactly x). */
static bool contains_log_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Log") && e->data.function.arg_count == 1 &&
        is_symbol(e->data.function.args[0], x)) return true;
    if (contains_log_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_log_x(e->data.function.args[i], x)) return true;
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
static bool is_zero_simp_(Expr* e) {   /* consumes e */
    Expr* s = simp(e);
    bool z = (s->type == EXPR_INTEGER && s->data.integer == 0);
    expr_free(s);
    return z;
}
static bool prove_pos(Expr* val, Expr* as) {   /* is val > 0 provable? consumes val */
    Expr* pred = mk_fn2("Greater", val, mk_int(0));
    Expr* r = as ? eval_take(mk_fn2("Refine", pred, cp(as))) : ev1("Refine", pred);
    bool t = sym_is(r, "True");
    expr_free(r);
    return t;
}
static int real_double(Expr* e, double* out) {
    Expr* n = ev1("N", cp(e));
    if (!n) return 0;
    int ok = 0;
    if (n->type == EXPR_REAL)         { *out = n->data.real;            ok = 1; }
    else if (n->type == EXPR_INTEGER) { *out = (double) n->data.integer; ok = 1; }
    expr_free(n);
    return ok;
}

/* Collect the free parameter symbols of `e` (symbols in argument position that
 * are not x and not a known constant) into out[], up to cap; returns the count
 * (or cap+1 to signal overflow -- caller declines). */
static bool is_known_const(const char* s) {
    static const char* K[] = {"Pi","E","EulerGamma","Catalan","Degree","GoldenRatio",
        "I","Infinity","ComplexInfinity","Indeterminate","Glaisher","Khinchin", NULL};
    for (int i = 0; K[i]; i++) if (strcmp(s, K[i]) == 0) return true;
    return false;
}
static void collect_params(const Expr* e, const Expr* x, Expr** out, size_t* n, size_t cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (is_symbol(e, x) || is_known_const(e->data.symbol.name)) return;
        for (size_t i = 0; i < *n && i < cap; i++)
            if (out[i]->data.symbol.name == e->data.symbol.name) return;  /* dup */
        if (*n < cap) out[*n] = (Expr*)e;
        (*n)++;   /* may exceed cap -> caller sees overflow */
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    /* arguments only -- a head symbol (Log, Power, Plus, f) is not a parameter */
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_params(e->data.function.args[i], x, out, n, cap);
}

/* ---- the method --------------------------------------------------------- */

Expr* integrate_feynparam_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    /* Interval exactly [0,1] in v1. */
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;
    if (!contains_symbol(f, x)) return NULL;
    /* Needs a 1/Log[x] weight -- cheap structural gate. */
    if (!contains_log_x(f, x)) return NULL;

    /* g = f * Log[x]; the Log[x] denominator must cancel cleanly. */
    Expr* g = eval_take(Times_(cp(f), mk_fn1("Log", cp(x))));
    if (!g || contains_log_x(g, x) || mentions_nonfinite(g)) {
        if (g) expr_free(g);
        return NULL;
    }

    /* Split g into summands; each must be c_i * x^(e_i) with c_i, e_i free of x. */
    Expr** terms; size_t nterms; Expr* one[1];
    if (head_name_is(g, "Plus")) { terms = g->data.function.args; nterms = g->data.function.arg_count; }
    else { one[0] = g; terms = one; nterms = 1; }
    if (nterms == 0 || nterms > 64) { expr_free(g); return NULL; }

    Expr* csum = mk_int(0);                 /* running Sum c_i (must be 0) */
    Expr* result = mk_int(0);               /* running Sum c_i Log[e_i+1] */
    bool ok = true;

    for (size_t i = 0; i < nterms && ok; i++) {
        Expr* s = terms[i];
        /* e_i = Exponent[s, x]; must be free of x. */
        Expr* ei = eval_take(mk_fn2("Exponent", cp(s), cp(x)));
        if (!ei || contains_symbol(ei, x) || mentions_nonfinite(ei)) { if (ei) expr_free(ei); ok = false; break; }
        /* c_i = Simplify[s * x^(-e_i)]; must be free of x. */
        Expr* ci = simp(Times_(cp(s), Pow_(cp(x), Times_(mk_int(-1), cp(ei)))));
        if (!ci || contains_symbol(ci, x) || mentions_nonfinite(ci)) {
            if (ci) expr_free(ci);
            expr_free(ei); ok = false; break;
        }
        /* Require e_i + 1 > 0 (convergence of the x-power integral at x=0). */
        if (!prove_pos(Plus_(cp(ei), mk_int(1)), assumptions)) {
            expr_free(ci); expr_free(ei); ok = false; break;
        }
        csum = Plus_(csum, cp(ci));
        /* += c_i Log[e_i + 1] */
        result = Plus_(result, Times_(ci, mk_fn1("Log", Plus_(ei, mk_int(1)))));
    }
    expr_free(g);
    if (!ok) { expr_free(csum); expr_free(result); return NULL; }

    /* Convergence at x = 1 requires Sum c_i = 0. */
    if (!is_zero_simp_(csum)) { expr_free(result); return NULL; }   /* csum consumed */

    Expr* val = assumptions ? simp2(result, assumptions) : simp(result);
    if (!val || contains_symbol(val, x) || mentions_nonfinite(val)) {
        if (val) expr_free(val);
        return NULL;
    }

    /* Numeric backstop: instantiate any free parameters at an in-domain sample
     * and compare N[val] with NIntegrate[f]. */
    Expr* params[8]; size_t np = 0;
    collect_params(f, x, params, &np, 8);
    if (np > 8) { expr_free(val); return NULL; }

    static const long cand[2][8] = {
        { 2, 3, 5, 7, 11, 13, 17, 19 },
        { 3, 5, 7, 11, 13, 17, 19, 23 },
    };
    bool verified = false;
    for (int set = 0; set < 2 && !verified; set++) {
        /* Build substitution rules params[i] -> cand[set][i]. */
        Expr* rules[8];
        for (size_t i = 0; i < np; i++)
            rules[i] = mk_fn2("Rule", cp(params[i]), mk_int(cand[set][i]));
        Expr* sub = (np > 0) ? mk_fn("List", rules, np) : NULL;

        /* Assumptions must hold at this sample (skip when none). */
        bool assum_ok = true;
        if (assumptions && sub) {
            Expr* chk = eval_take(mk_fn1("Refine",
                            eval_take(mk_fn2("ReplaceAll", cp(assumptions), cp(sub)))));
            assum_ok = sym_is(chk, "True");
            expr_free(chk);
        }
        if (!assum_ok) { if (sub) expr_free(sub); continue; }

        Expr* fn  = sub ? eval_take(mk_fn2("ReplaceAll", cp(f), cp(sub))) : cp(f);
        Expr* vn  = sub ? eval_take(mk_fn2("ReplaceAll", cp(val), cp(sub))) : cp(val);
        if (sub) expr_free(sub);

        double vd;
        int vok = real_double(vn, &vd);
        expr_free(vn);
        if (!vok) { expr_free(fn); continue; }

        arith_warnings_mute_push();
        Expr* ni = eval_take(mk_fn2("NIntegrate", fn,
                        mk_fn3b("List", cp(x), mk_int(0), mk_int(1))));   /* fn consumed */
        arith_warnings_mute_pop();
        double nd;
        int nok = ni ? real_double(ni, &nd) : 0;
        if (ni) expr_free(ni);
        if (nok && fabs(vd - nd) <= 1e-3 * (1.0 + fabs(vd))) verified = true;
    }

    if (!verified) { expr_free(val); return NULL; }
    return val;
}

Expr* builtin_integrate_feynparam(Expr* res) {
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
    return integrate_feynparam_try(f, x, a, b, as);
}

void integrate_feynparam_init(void) {
    symtab_add_builtin("Integrate`FeynmanParameter", builtin_integrate_feynparam);
    symtab_get_def("Integrate`FeynmanParameter")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`FeynmanParameter",
        "Integrate`FeynmanParameter[f, {x, 0, 1}] evaluates the parameter / "
        "Frullani-logarithm family Integrate[(Sum c_i x^(e_i))/Log[x], {x,0,1}] = "
        "Sum c_i Log[e_i+1] (differentiation under the integral sign w.r.t. the "
        "exponent, with the inner x-power integral taken as 1/(e+1)).  Requires "
        "Sum c_i = 0 (so the 1/Log[x] singularity at x=1 is integrable) and each "
        "e_i+1 > 0 (proved from Assumptions).  E.g. (x^a-1)/Log[x] -> Log[a+1], "
        "(x^a-x^b)/Log[x] -> Log[(a+1)/(b+1)].  Numerically cross-checked against "
        "NIntegrate; returns unevaluated otherwise.");
}
