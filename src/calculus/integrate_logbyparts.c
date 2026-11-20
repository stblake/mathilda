/* integrate_logbyparts.c — Log[g] times a trig/hyperbolic kernel -> one
 * integration by parts.
 *
 * Recognizes a product  c Log[g(x)] K(x)  with exactly one Log factor, where K
 * is free of Log and carries at least one trig or hyperbolic kernel of x, and
 * integrates it by a single integration by parts with u = Log[g], dv = K dx:
 *
 *     INT Log[g] K dx  =  Log[g] V  -  INT V (g'/g) dx,    V = INT K dx.
 *
 * On the family this targets — Log[x] Sin[a x], Log[x] Cos[a x], and the
 * polynomial/constant-weighted variants — V is an elementary trig antiderivative
 * and V (g'/g) = (trig)/x closes to SinIntegral / CosIntegral, so the whole
 * integral is a clean Log*trig + Si/Ci form.
 *
 * Why this exists as its own recognizer.  A Log*trig integrand is outside every
 * elementary stage of the cascade and the general search does not merely decline,
 * it GRINDS: ParallelMixedSpecial (the last, special-function stage) closes
 * Integrate[Log[x] Sin[x], x] correctly but in ~3.8 s, and the sibling
 * Integrate[x^3 Sin[x] Log[x]^2, x] in ~12 s.  DSolve[y'' + 4y == Log[x], y, x]
 * needs exactly the pair INT Log[x] Sin[2x] and INT Log[x] Cos[2x] for its
 * variation-of-parameters particular, so the two slow searches were the whole
 * ~9.4 s cost of that ODE (over the corpus harness's 8 s wall).  This stage runs
 * AFTER the cheap elementary stages (the cascade's `if (!result)` short-circuit
 * means it only ever sees integrands those stages all declined) and BEFORE the
 * ParallelMixedTower / ParallelMixedSpecial tail, so it costs the slow path
 * nothing it was not already going to pay, and turns the Log*trig closes fast.
 *
 * Correctness rests on an exact Simplify diff-back (the one acceptance test, as
 * in integrate_gammapower.c), so a mis-recognition can only decline, never emit a
 * wrong closed form.  Its recursive sub-integrals run with g_integrate_no_special
 * raised, so a Log*trig case whose IBP residual is itself non-elementary declines
 * promptly instead of paying a second ParallelMixedSpecial search.
 */

#include "integrate_logbyparts.h"
#include "integrate.h"       /* g_integrate_no_special */

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

static Expr* lb_sym(const char* s) { return expr_new_symbol(s); }
static Expr* lb_fn1(const char* h, Expr* a) {
    return expr_new_function(lb_sym(h), (Expr*[]){ a }, 1);
}
static Expr* lb_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(lb_sym(h), (Expr*[]){ a, b }, 2);
}
/* evaluate() BORROWS its argument, so free the input node. */
static Expr* lb_ev(Expr* node) { Expr* r = evaluate(node); expr_free(node); return r; }
static Expr* lb_ev2(const char* h, Expr* a, Expr* b) { return lb_ev(lb_fn2(h, a, b)); }

static bool lb_free_of(const Expr* e, const Expr* x) {
    if (!e) return true;
    if (expr_eq((Expr*)e, (Expr*)x)) return false;
    if (e->type != EXPR_FUNCTION) return true;
    if (!lb_free_of(e->data.function.head, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!lb_free_of(e->data.function.args[i], x)) return false;
    return true;
}

/* True when some node of e is a function with head symbol named `name`. */
static bool lb_has_head(const Expr* e, const char* name) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const char* iname = intern_symbol(name);
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == iname)
        return true;
    if (lb_has_head(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lb_has_head(e->data.function.args[i], name)) return true;
    return false;
}

/* True when e contains a trig/hyperbolic kernel with an x-dependent argument. */
static bool lb_has_trig_of_x(const Expr* e, const Expr* x) {
    static const char* const heads[] = {
        "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
        "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch"
    };
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        for (size_t k = 0; k < sizeof(heads) / sizeof(heads[0]); k++) {
            if (h == intern_symbol(heads[k])) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    if (!lb_free_of(e->data.function.args[i], x)) return true;
            }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lb_has_trig_of_x(e->data.function.args[i], x)) return true;
    return lb_has_trig_of_x(e->data.function.head, x);
}

/* Is e a one-argument Log[g]? */
static bool lb_is_log1(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1;
}

/* Split f = c Log[g] K into the Log argument g (fresh) and the cofactor K
 * (fresh, the product of all non-Log factors).  Requires EXACTLY one Log factor
 * (a bare one-argument Log[g], not Log[g]^k), g x-dependent and free of Log, and
 * every other factor free of Log.  Returns false otherwise. */
static bool lb_split(Expr* f, Expr* x, Expr** pg, Expr** pK) {
    const char* s_times = intern_symbol("Times");
    size_t nf = 1;
    Expr* one[1] = { f };
    Expr** fac = one;
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == s_times) {
        fac = f->data.function.args;
        nf  = f->data.function.arg_count;
    }

    Expr* g = NULL;
    Expr** keep = malloc(nf * sizeof(Expr*));
    size_t nkeep = 0;
    for (size_t i = 0; i < nf; i++) {
        Expr* t = fac[i];
        if (!g && lb_is_log1(t)) {
            Expr* arg = t->data.function.args[0];
            if (!lb_free_of(arg, x) && !lb_has_head(arg, "Log")) {
                g = expr_copy(arg);
                continue;
            }
        }
        /* a second Log anywhere in a non-Log factor disqualifies the split */
        if (lb_has_head(t, "Log")) {
            if (g) expr_free(g);
            for (size_t j = 0; j < nkeep; j++) expr_free(keep[j]);
            free(keep);
            return false;
        }
        keep[nkeep++] = expr_copy(t);
    }

    if (!g || nkeep == 0) {
        if (g) expr_free(g);
        for (size_t j = 0; j < nkeep; j++) expr_free(keep[j]);
        free(keep);
        return false;
    }
    *pg = g;
    *pK = (nkeep == 1) ? keep[0]
        : lb_ev(expr_new_function(lb_sym("Times"), keep, nkeep));
    free(keep);
    return true;
}

Expr* integrate_logbyparts_try(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;

    Expr* g = NULL; Expr* K = NULL;
    if (!lb_split(f, x, &g, &K)) return NULL;

    /* Gate: confine to the class the slow stages grind on — a trig/hyperbolic
     * kernel of x in the cofactor.  Log*polynomial, Log*exp and Log*rational all
     * already close cheaply (earlier stages, so they never reach here anyway),
     * and this keeps the recursive sub-integrals in the elementary-trig regime. */
    if (!lb_has_trig_of_x(K, x)) { expr_free(g); expr_free(K); return NULL; }

    /* V = INT K dx and W = INT V (g'/g) dx, both run with the heavy tail stages
     * suppressed: on this class they close at a cheap stage, and a case where
     * they do not must decline here rather than pay ParallelMixedSpecial twice. */
    g_integrate_no_special++;
    Expr* V = lb_ev2("Integrate", expr_copy(K), expr_copy(x));
    bool v_ok = V && !lb_has_head(V, "Integrate") && !lb_has_head(V, "Log")
             && !lb_has_head(V, "Inactive");
    Expr* W = NULL;
    if (v_ok) {
        Expr* gp    = lb_ev2("D", lb_fn1("Log", expr_copy(g)), expr_copy(x));
        Expr* resid = lb_ev2("Times", expr_copy(V), gp);
        W = lb_ev2("Integrate", resid, expr_copy(x));
    }
    g_integrate_no_special--;

    bool w_ok = v_ok && W && !lb_has_head(W, "Integrate") && !lb_has_head(W, "Inactive");
    if (!w_ok) { expr_free(g); expr_free(K); expr_free(V); expr_free(W); return NULL; }

    /* result = Log[g] V - W */
    Expr* result = lb_ev2("Subtract",
                       lb_ev2("Times", lb_fn1("Log", g), V),   /* consumes g, V */
                       W);
    expr_free(K);

    /* Diff-back is the only acceptance test.  D[SinIntegral[z]] prints as Sinc[z]
     * rather than Sin[z]/z, which Simplify does not expand on its own, so rewrite
     * it before the comparison; this is a spelling identity (Sinc[z] == Sin[z]/z),
     * not a weakening of the certificate. */
    Expr* d    = lb_ev2("D", expr_copy(result), expr_copy(x));
    Expr* raw  = lb_ev2("Subtract", d, expr_copy(f));
    Expr* rule = parse_expression("Sinc[a_] :> Sin[a]/a");
    Expr* repl = lb_ev2("ReplaceAll", raw, rule);
    Expr* diff = lb_ev(lb_fn1("Simplify", expr_copy(repl)));
    bool ok = (diff && diff->type == EXPR_INTEGER && diff->data.integer == 0);
    expr_free(diff);
    /* Simplify leaves Cosh/Sinh beside E^x uncombined (a hyperbolic residual of a
     * PolyLog answer, Log[1 + E^x] Tanh[x]^2): the zero test decides */
    if (!ok) {
        Expr* z = lb_ev(lb_fn1("PossibleZeroQ", expr_copy(repl)));
        ok = z && z->type == EXPR_SYMBOL && z->data.symbol.name == intern_symbol("True");
        if (z) expr_free(z);
    }
    expr_free(repl);
    if (!ok) { expr_free(result); return NULL; }
    return result;
}

/* ------------------------------------------------------------------------ */
/* Special-function kernels by parts (Campaign D7).                          */
/*                                                                          */
/*   INT S K dx  =  S V  -  INT V S' dx,    V = INT K dx elementary,          */
/*                                                                          */
/* for exactly one factor S in {ExpIntegralEi, Erf, Erfi, Erfc, LogIntegral}  */
/* of an x-dependent argument, K free of those heads.  S' is elementary      */
/* (E^u u'/u, 2 E^(-u^2) u'/Sqrt[Pi], u'/Log[u], ...), so V S' is elementary  */
/* and the cascade closes INT V S' with another special function where one   */
/* exists: INT E^x Ei[x] = E^x Ei[x] - Ei[2 x], INT Log[x] li(x) =           */
/* (x Log[x] - x) li(x) - x^2/2 + li(x^2)... The Risch decision over the     */
/* tower with S as a primitive monomial proves such integrands non-elementary */
/* (or leaves them undecided) and has no closed form to give; FriCAS answers  */
/* them this way.  Accepted only on PossibleZeroQ[D[F] - f].                 */
/* ------------------------------------------------------------------------ */

static bool sb_is_special(const Expr* e, const Expr* x) {
    static const char* const heads[] = { "ExpIntegralEi", "Erf", "Erfi", "Erfc", "LogIntegral" };
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    for (size_t k = 0; k < sizeof(heads) / sizeof(heads[0]); k++)
        if (h == intern_symbol(heads[k])) return !lb_free_of(e->data.function.args[0], x);
    return false;
}
static bool sb_has_special(const Expr* e) {
    static const char* const heads[] = { "ExpIntegralEi", "Erf", "Erfi", "Erfc", "LogIntegral",
                                         "ExpIntegralE", "SinIntegral", "CosIntegral",
                                         "SinhIntegral", "CoshIntegral", "PolyLog", "Gamma" };
    for (size_t k = 0; k < sizeof(heads) / sizeof(heads[0]); k++)
        if (lb_has_head(e, heads[k])) return true;
    return false;
}

Expr* integrate_specialbyparts_try(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;
    const char* s_times = intern_symbol("Times");
    size_t nf = 1;
    Expr* one[1] = { f };
    Expr** fac = one;
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == s_times) {
        fac = f->data.function.args;
        nf  = f->data.function.arg_count;
    }
    Expr* S = NULL;
    Expr** keep = malloc(nf * sizeof(Expr*));
    size_t nkeep = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* t = fac[i];
        if (!S && sb_is_special(t, x) && !sb_has_special(t->data.function.args[0])) { S = t; continue; }
        if (sb_has_special(t)) ok = false;   /* a second special factor, or S^2 */
        else keep[nkeep++] = t;
    }
    if (!ok || !S) { free(keep); return NULL; }
    /* keep[] holds pointers borrowed from f: the product gets its own copies */
    Expr* K;
    if (nkeep == 0) K = expr_new_integer(1);
    else if (nkeep == 1) K = expr_copy(keep[0]);
    else {
        for (size_t j = 0; j < nkeep; j++) keep[j] = expr_copy(keep[j]);
        K = lb_ev(expr_new_function(lb_sym("Times"), keep, nkeep));
    }
    free(keep);
    /* both sub-integrals with the heavy tail stages suppressed (as LogByParts):
     * a residual that does not close cheaply declines here */
    g_integrate_no_special++;
    Expr* V = lb_ev2("Integrate", K, expr_copy(x));   /* consumes K */
    bool v_ok = V && !lb_has_head(V, "Integrate") && !lb_has_head(V, "Inactive") && !sb_has_special(V);
    if (!v_ok) { g_integrate_no_special--; if (V) expr_free(V); return NULL; }
    Expr* Sp = lb_ev2("D", expr_copy(S), expr_copy(x));
    Expr* W = lb_ev2("Integrate", lb_ev2("Times", expr_copy(V), Sp), expr_copy(x));
    g_integrate_no_special--;
    if (!W || lb_has_head(W, "Integrate") || lb_has_head(W, "Inactive")) {
        if (W) expr_free(W);
        expr_free(V);
        return NULL;
    }
    Expr* result = lb_ev2("Subtract", lb_ev2("Times", expr_copy(S), V), W);   /* consumes V, W */
    Expr* d = lb_ev2("Subtract", lb_ev2("D", expr_copy(result), expr_copy(x)), expr_copy(f));
    Expr* z = lb_ev(lb_fn1("PossibleZeroQ", d));
    bool zok = z && z->type == EXPR_SYMBOL && z->data.symbol.name == intern_symbol("True");
    if (z) expr_free(z);
    if (!zok) { expr_free(result); return NULL; }
    return result;
}
