/* intrat.c — rational-function integration package (`Integrate``).
 *
 * Phase 1 of docs/notes/plans/INTEGRATE_PLAN.md.  Implements the scaffolding,
 * polynomial integration, Mack's linear Hermite reduction, the
 * derivative-recognition fast path (`c*D'/D^k`), and the four
 * one-line helpers (Content / Primitive / Monic / LeadingCoefficient)
 * used by the later phases.  Everything is exposed under the
 * `Integrate`` context so each stage is REPL-testable in isolation.
 *
 * Algorithm references:
 *   - Bronstein, *Symbolic Integration I*, Springer 2nd ed., 2004.
 *   - Trager, *Algebraic Factoring and Rational Function
 *     Integration*, ACM SYMSAC 1976.
 * The Mathematica baseline being ported is `IntegrateRational.m` in
 * the repo root — line numbers in the comments below refer to it.
 *
 * Memory contract: every public BuiltinFunc follows the Mathilda
 * convention — caller (evaluator) owns `res`; the function returns a
 * freshly-allocated Expr* on success or NULL on failure.
 */

#include "intrat.h"
#include "parse.h"
#include "intrat_internal.h"
#include "intrat_subres.h"
#include "intrat_alg.h"
#include "intsimp.h"
#include "common.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "poly.h"
#include "expand.h"
#include "sym_names.h"
#include "print.h"
#include "options.h"
#include "rationalize.h"
#include "root.h"
#include "core.h"
#include "flint_qqbar.h"
#ifdef USE_FLINT
#include "flint_bridge_internal.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* Trace scaffolding (matches Mathematica `debugprint`).               */
/*                                                                    */
/* When `Integrate``$Verbose` evaluates to True the trace helper       */
/* prints an `IN ...` line on entry / `OUT ...` line on exit of every  */
/* public helper.  Cheap when off (single OwnValue lookup).            */
/* ------------------------------------------------------------------ */

static bool intrat_trace_enabled(void) {
    SymbolDef* def = symtab_lookup("Integrate`$Verbose");
    if (!def || !def->own_values) return false;
    /* The OwnValue is `Integrate`$Verbose -> True/False`.  Walk it. */
    Rule* r = def->own_values;
    while (r) {
        if (r->replacement && r->replacement->type == EXPR_SYMBOL
            && r->replacement->data.symbol.name == SYM_True) return true;
        r = r->next;
    }
    return false;
}

static void intrat_trace(const char* fn, const char* dir, Expr* e) {
    if (!intrat_trace_enabled()) return;
    char* s = expr_to_string(e);
    fprintf(stderr, "[Integrate`] %s %s %s\n", fn, dir, s ? s : "(null)");
    if (s) free(s);
}

/* ------------------------------------------------------------------ */
/* Small wrappers around Mathilda internal_* / evaluate idioms.        */
/*                                                                    */
/* All of these hand the result back as a fresh Expr* (the caller     */
/* owns it).  None of them free their input arguments — every helper  */
/* makes its own copies.                                              */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/* Algorithm choices: Method -> {"BronsteinRational", "Hermite" -> ...,  */
/* "LogPart" -> ..., "RealForm" -> ...} and the same rules as trailing  */
/* options of Integrate`BronsteinRational.  An installed setting is      */
/* valid for the top-level evaluation that installed it only (a         */
/* TimeConstrained unwind that skips the restore cannot leak it into     */
/* the next input).                                                      */
/*                                                                      */
/* The Automatic defaults were chosen by measurement on benchmark 266    */
/* (benchmarks/266-rational-integration/REPORT.md, "Algorithm grid").    */
/* ------------------------------------------------------------------ */
#define INTRAT_DEFAULT_HERMITE  INTRAT_HERMITE_LINEAR
#define INTRAT_DEFAULT_LOGPART  INTRAT_LOGPART_RT

static IntratMethodOptions g_ir_opts = { -1, -1, -1 };
static uint64_t g_ir_opts_stamp = 0;

void intrat_method_options_default(IntratMethodOptions* o) {
    o->hermite = o->logpart = o->realform = -1;
}

static IntratMethodOptions ir_opts(void) {
    IntratMethodOptions o;
    if (g_ir_opts_stamp == eval_toplevel_id()) o = g_ir_opts;
    else intrat_method_options_default(&o);
    return o;
}

void intrat_method_options_install(const IntratMethodOptions* o, IntratMethodOptions* saved) {
    *saved = ir_opts();
    g_ir_opts = *o;
    g_ir_opts_stamp = eval_toplevel_id();
}

void intrat_method_options_restore(const IntratMethodOptions* saved) {
    g_ir_opts = *saved;
    g_ir_opts_stamp = eval_toplevel_id();
}

int intrat_method_suboption(const Expr* rule, IntratMethodOptions* o) {
    if (!rule || rule->type != EXPR_FUNCTION || rule->data.function.arg_count != 2
        || rule->data.function.head->type != EXPR_SYMBOL
        || (rule->data.function.head->data.symbol.name != SYM_Rule
            && rule->data.function.head->data.symbol.name != SYM_RuleDelayed)) return 0;
    const Expr* k = rule->data.function.args[0];
    const Expr* v = rule->data.function.args[1];
    if (k->type != EXPR_STRING) return 0;
    const char* vs = v->type == EXPR_STRING ? v->data.string : NULL;
    bool vauto = v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic;
    bool vnone = v->type == EXPR_SYMBOL && strcmp(v->data.symbol.name, "None") == 0;
    if (strcmp(k->data.string, "Hermite") == 0) {
        if (vauto) { o->hermite = -1; return 1; }
        if (!vs) return 0;
        if (strcmp(vs, "Linear") == 0 || strcmp(vs, "Mack") == 0) { o->hermite = INTRAT_HERMITE_LINEAR; return 1; }
        if (strcmp(vs, "Quadratic") == 0) { o->hermite = INTRAT_HERMITE_QUADRATIC; return 1; }
        if (strcmp(vs, "HorowitzOstrogradsky") == 0) { o->hermite = INTRAT_HERMITE_HOROWITZ; return 1; }
        return 0;
    }
    if (strcmp(k->data.string, "LogPart") == 0) {
        if (vauto) { o->logpart = -1; return 1; }
        if (!vs) return 0;
        if (strcmp(vs, "LazardRiobooTrager") == 0) { o->logpart = INTRAT_LOGPART_LRT; return 1; }
        if (strcmp(vs, "RothsteinTrager") == 0) { o->logpart = INTRAT_LOGPART_RT; return 1; }
        if (strcmp(vs, "Czichowski") == 0) { o->logpart = INTRAT_LOGPART_CZICHOWSKI; return 1; }
        return 0;
    }
    if (strcmp(k->data.string, "RealForm") == 0) {
        if (vauto) { o->realform = -1; return 1; }
        if (vnone) { o->realform = INTRAT_REALFORM_NONE; return 1; }
        if (!vs) return 0;
        if (strcmp(vs, "Rioboo") == 0) { o->realform = INTRAT_REALFORM_RIOBOO; return 1; }
        if (strcmp(vs, "LogToReal") == 0) { o->realform = INTRAT_REALFORM_LOGTOREAL; return 1; }
        return 0;
    }
    return 0;
}

static Expr* call_eval(const char* head, Expr** args, size_t argc) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, argc);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Expand then free the input.  expr_expand borrows its argument and
 * returns a fresh tree, so any caller that owns a freshly-built temp
 * (e.g. internal_times(...) on the spot) must free it after the
 * expand or the temp leaks.  This wraps the common idiom. */
static inline Expr* expand_and_free(Expr* e) {
    if (!e) return NULL;
    Expr* r = expr_expand(e);
    expr_free(e);
    return r;
}

/* canonic[expr] = Cancel[Together[expr]].  Phase 1 keeps this
 * stripped down — the full Mathematica `canonic` adds a RootReduce /
 * ToRadicals tail that we'll bring online in later phases.
 *
 * The full Mathematica canonic also switches to
 * `Extension -> Automatic` only when the input contains a radical
 * (intsimp_has_radical); the Mathilda port stays plain-rational for
 * Phase 1 since our Cancel/Together don't yet honour Extension. */
Expr* intrat_canonic(Expr* e) {
    Expr* tg = internal_together((Expr*[]){expr_copy(e)}, 1);
    Expr* result = internal_cancel((Expr*[]){tg}, 1);
    return result;
}

static Expr* intrat_d(Expr* e, Expr* x) {
    Expr* args[2] = {expr_copy(e), expr_copy(x)};
    return call_eval("D", args, 2);
}

static Expr* intrat_numerator(Expr* e) {
    return internal_numerator((Expr*[]){expr_copy(e)}, 1);
}

Expr* intrat_denominator(Expr* e) {
    return internal_denominator((Expr*[]){expr_copy(e)}, 1);
}

static Expr* intrat_polyq(Expr* a, Expr* b, Expr* x) {
    return internal_polynomialquotient(
        (Expr*[]){expr_copy(a), expr_copy(b), expr_copy(x)}, 3);
}

/* PolynomialRemainder[a, b, var] over the named variable.  Used by
 * the Phase 2 LRT pipeline to reduce mod Q_i(t). */
static Expr* intrat_polyrem_t(Expr* a, Expr* b, Expr* var) {
    return internal_polynomialremainder(
        (Expr*[]){expr_copy(a), expr_copy(b), expr_copy(var)}, 3);
}

/* PolynomialQuotientRemainder, returning the {q, r} components separately. */
static bool intrat_polyqr(Expr* a, Expr* b, Expr* x, Expr** out_q, Expr** out_r) {
    Expr* args[3] = {expr_copy(a), expr_copy(b), expr_copy(x)};
    Expr* call = expr_new_function(
        expr_new_symbol(SYM_PolynomialQuotientRemainder), args, 3);
    Expr* result = evaluate(call);
    expr_free(call);
    if (!result || result->type != EXPR_FUNCTION
        || result->data.function.head->type != EXPR_SYMBOL
        || result->data.function.head->data.symbol.name != SYM_List
        || result->data.function.arg_count != 2) {
        if (result) expr_free(result);
        *out_q = NULL; *out_r = NULL;
        return false;
    }
    *out_q = expr_copy(result->data.function.args[0]);
    *out_r = expr_copy(result->data.function.args[1]);
    expr_free(result);
    return true;
}

/* exquo[a, b, x] = quotient asserting zero remainder.  Bronstein's
 * "exact quotient" — IntegrateRational.m:1238-1248.  On a non-zero
 * remainder we emit a soft warning (matching the reference) and
 * return the quotient anyway; downstream consumers will detect the
 * symbolic remainder during their own canonicalisation. */
static Expr* intrat_exquo(Expr* a, Expr* b, Expr* x) {
    Expr *q = NULL, *r = NULL;
    if (!intrat_polyqr(a, b, x, &q, &r)) return NULL;
    /* A nonzero remainder means a caller's divisibility premise failed (an
     * undecided zero test over parameters or radicals).  The quotient used
     * to be returned anyway, after a raw stderr line that bypassed Quiet and
     * Check, and Hermite reduction / the S_i correction carried on with wrong
     * data.  Decline: every caller treats NULL as "cannot integrate". */
    bool exact = !r || is_zero_poly(r);
    if (r) expr_free(r);
    if (!exact) { if (q) expr_free(q); return NULL; }
    return q;
}

static bool intrat_polyq_test(Expr* expr, Expr* var) {
    Expr* args[2] = {expr_copy(expr), expr_copy(var)};
    Expr* call = internal_polynomialq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
    if (val) expr_free(val);
    return ok;
}

bool intrat_freeq_test(Expr* expr, Expr* var) {
    Expr* args[2] = {expr_copy(expr), expr_copy(var)};
    Expr* call = internal_freeq(args, 2);
    Expr* val  = evaluate(call);
    expr_free(call);
    bool ok = (val && val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
    if (val) expr_free(val);
    return ok;
}

/* Monic GCD in K(coeffs)[x].  Mathematica's HermiteReduce baseline
 * uses First[PolynomialExtendedGCD[a, b, x]] which returns the GCD
 * normalised so leading_coefficient(g, x) == 1.  Using the content-
 * normalised PolynomialGCD instead breaks the algorithm: when D has
 * non-trivial integer content (e.g. d = 4 x^2 (5-4x)^2) the dbar
 * term silently scales by content(d), the dstar / dbarhat / dtil
 * triplet downstream gets out of sync, and the f = D[g, x] + h
 * identity fails by a constant factor.
 *
 * Implementation: compute PolynomialGCD, then divide by the leading
 * coefficient in x.  PolynomialExtendedGCD itself is correct here but
 * orders of magnitude slower on hard inputs because it also tracks
 * Bezout coefficients we never use.
 *
 * Result is expanded — every downstream consumer feeds the gcd into
 * get_degree_poly, which sees through `Power[var, n]` only, not
 * `Power[base, n]` for multi-term `base`. */
/* True when e carries a radical of a SYMBOLIC base free of x, such as
 * Sqrt[b^2 + c^2]: PolynomialGCD treats it as an independent variable, so a
 * gcd that holds only through the relation Sqrt[p]^2 = p is missed. */
static bool intrat_has_symbolic_radical(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* r = e->data.function.args[1];
        bool frac_exp = r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
             r->data.function.head->data.symbol.name == SYM_Rational;
        if (frac_exp && b->type != EXPR_INTEGER && b->type != EXPR_BIGINT &&
            !(b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL &&
              b->data.function.head->data.symbol.name == SYM_Rational) &&
            intrat_freeq_test((Expr*)b, (Expr*)x))
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_has_symbolic_radical(e->data.function.args[i], x)) return true;
    return false;
}

/* gcd(a, b) in K[x] by Euclid's algorithm, every remainder put over a common
 * denominator (Together), so the evaluator's Sqrt[p]^2 -> p reduction decides
 * which remainders vanish.  The gcd of (Sqrt[b^2+c^2] - b) u^2 + 2 c u +
 * Sqrt[b^2+c^2] + b and its derivative is the linear factor -- the quadratic is
 * a square in the field -- where PolynomialGCD answered 1, Hermite reduction
 * then left the "squarefree" quadratic to Rothstein-Trager, whose resultant had
 * a zero discriminant: Integrate[1/(Sqrt[b^2 + c^2] + b Cos[x] + c Sin[x]), x]
 * came out 0 after 1/0 and 0 ComplexInfinity messages.  NULL when a remainder's
 * degree fails to drop (an undecided zero test): the caller keeps PolynomialGCD. */
static Expr* intrat_euclid_gcd(Expr* a, Expr* b, Expr* x) {
    Expr* u = expr_copy(a);
    Expr* v = expr_copy(b);
    for (int it = 0; it < 64; it++) {
        int dv = get_degree_poly(v, x);
        if (is_zero_poly(v)) { expr_free(v); return u; }
        if (dv <= 0) { expr_free(u); expr_free(v); return expr_new_integer(1); }
        Expr* r = intrat_polyrem_t(u, v, x);
        Expr* rt = r ? internal_together((Expr*[]){r}, 1) : NULL;
        Expr* re = rt ? expr_expand(rt) : NULL;
        if (rt) expr_free(rt);
        if (!re) { expr_free(u); expr_free(v); return NULL; }
        Expr* rn = internal_numerator((Expr*[]){expr_copy(re)}, 1);
        Expr* rne = rn ? expr_expand(rn) : NULL;
        if (rn) expr_free(rn);
        expr_free(re);
        if (!rne) { expr_free(u); expr_free(v); return NULL; }
        if (!is_zero_poly(rne) && get_degree_poly(rne, x) >= dv) {
            expr_free(rne); expr_free(u); expr_free(v); return NULL;
        }
        expr_free(u);
        u = v;
        v = rne;     /* the numerator: a scalar multiple of the remainder */
    }
    expr_free(u); expr_free(v);
    return NULL;
}

static Expr* intrat_polygcd_monic(Expr* a, Expr* b, Expr* x) {
    Expr* g = NULL;
    if (intrat_has_symbolic_radical(a, x) || intrat_has_symbolic_radical(b, x))
        g = intrat_euclid_gcd(a, b, x);
    if (!g) g = internal_polynomialgcd(
        (Expr*[]){expr_copy(a), expr_copy(b)}, 2);
    Expr* g_e = expr_expand(g);
    expr_free(g);
    if (!g_e || is_zero_poly(g_e)) return g_e;
    int deg = get_degree_poly(g_e, x);
    if (deg <= 0) {
        /* Degree-0 gcd: normalise to 1 to keep the algorithm scale-
         * invariant (and so PolynomialQuotient[d, gcd] returns d/lc
         * faithfully). */
        expr_free(g_e);
        return expr_new_integer(1);
    }
    Expr* lc = get_coeff(g_e, x, deg);
    if (!lc) return g_e;
    /* If lc is already 1, nothing to do. */
    if (lc->type == EXPR_INTEGER && lc->data.integer == 1) {
        expr_free(lc);
        return g_e;
    }
    Expr* div = internal_divide((Expr*[]){g_e, lc}, 2);
    Expr* expanded = expr_expand(div);
    expr_free(div);
    return expanded;
}

/* lc[p, x] = leading coefficient of p in x.  When deg=0 the polynomial
 * is the constant itself. */
static Expr* intrat_lc(Expr* p, Expr* x) {
    int deg = get_degree_poly(p, x);
    if (deg <= 0) return expr_copy(p);
    return get_coeff(p, x, deg);
}

/* content[p, x] = PolynomialGCD @@ CoefficientList[p, x]. */
static Expr* intrat_content(Expr* p, Expr* x) {
    Expr* clist = internal_coefficientlist(
        (Expr*[]){expr_copy(p), expr_copy(x)}, 2);
    Expr* cl = evaluate(clist);
    expr_free(clist);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        return expr_new_integer(1);
    }
    size_t n = cl->data.function.arg_count;
    if (n == 0) { expr_free(cl); return expr_new_integer(0); }
    /* Fold via PolynomialGCD.  The single-variable form
     * PolynomialGCD[a, b, c, ...] handles the multi-arg form. */
    Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) args[i] = expr_copy(cl->data.function.args[i]);
    expr_free(cl);
    Expr* g = internal_polynomialgcd(args, n);
    /* internal_polynomialgcd (via internal_call_impl/expr_new_function)
     * memcpy's the slot pointers into its own backing store; the caller-
     * owned `args` buffer is no longer needed and must be released. */
    free(args);
    /* When coefficients carry parametric fractions (e.g. 1/b^3 factors
     * after Apart on a parametric integrand), our PolynomialGCD bails
     * and returns the call held.  Treat that as content == 1 — still
     * a correct (non-strict) primitive, and crucially keeps held GCD
     * calls from polluting downstream S polynomials. */
    if (head_is(g, SYM_PolynomialGCD)) {
        expr_free(g);
        g = expr_new_integer(1);
    }
    return g;
}

/* primitive[p, x] = p / content[p, x], expanded.  When p is free of x
 * the content equals p, so the result is 1; we shortcut to `p` to
 * match Mathematica's primitive[p, x] /; FreeQ[p, x] := p case. */
static Expr* intrat_primitive(Expr* p, Expr* x) {
    if (intrat_freeq_test(p, x)) return expr_copy(p);
    Expr* c = intrat_content(p, x);
    if (!c || (c->type == EXPR_INTEGER && c->data.integer == 1)) {
        if (c) expr_free(c);
        return expr_copy(p);
    }
    if (expr_is_numeric_like(c)) {
        Expr* div = internal_divide(
            (Expr*[]){expr_copy(p), c}, 2);
        Expr* expanded = expr_expand(div);
        expr_free(div);
        return expanded;
    }
    /* Polynomial content (in the other variables, e.g. the residue variable
     * t of a subresultant).  Expand[p/c] does NOT divide: it distributes 1/c
     * over the terms and leaves each one an uncancelled fraction, so the
     * "primitive part" comes back a rational function in t.  Downstream
     * (PrimitivePart_mod, LogToReal) then substitutes algebraic residue roots
     * into those denominators, and Together/PolynomialGCD over the resulting
     * nested-radical field stalls (issue #133: (x^2+x+1)/(x^4+1) hung).  The
     * content divides every coefficient exactly, so divide coefficient-wise
     * with Cancel and rebuild the polynomial in x. */
    Expr* cl = eval_and_free(internal_coefficientlist(
        (Expr*[]){expr_copy(p), expr_copy(x)}, 2));
    if (!cl || !head_is(cl, SYM_List)) {
        if (cl) expr_free(cl);
        expr_free(c);
        return expr_copy(p);
    }
    size_t n = cl->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t k = 0; k < n; k++) {
        Expr* q = eval_and_free(internal_cancel((Expr*[]){
            internal_divide((Expr*[]){
                expr_copy(cl->data.function.args[k]), expr_copy(c)}, 2)}, 1));
        Expr* qe = expr_expand(q); expr_free(q);
        terms[k] = internal_times((Expr*[]){qe,
            internal_power((Expr*[]){expr_copy(x),
                expr_new_integer((int64_t)k)}, 2)}, 2);
    }
    expr_free(cl);
    expr_free(c);
    Expr* sum = internal_plus(terms, n);
    free(terms);
    Expr* out = expr_expand(sum);
    expr_free(sum);
    return out;
}

/* monic[p, x] = p / lc[p, x], expanded.  When p is free of x we
 * follow Mathematica's convention and return 1. */
static Expr* intrat_monic_internal(Expr* p, Expr* x) {
    if (intrat_freeq_test(p, x)) return expr_new_integer(1);
    Expr* lc = intrat_lc(p, x);
    if (!lc) return expr_copy(p);
    if (lc->type == EXPR_INTEGER && lc->data.integer == 1) {
        expr_free(lc);
        return expr_expand(p);
    }
    Expr* div = internal_divide((Expr*[]){expr_copy(p), lc}, 2);
    Expr* result = expr_expand(div);
    expr_free(div);
    /* Cancel any rational coefficients that appear after the divide. */
    Expr* cancelled = internal_cancel((Expr*[]){result}, 1);
    return cancelled;
}

/* ------------------------------------------------------------------ */
/* ExtendedEuclidean (Diophantine).                                    */
/*                                                                    */
/* Given a, b, c with c in the ideal (a, b), return r, s' such that    */
/*    r*a + s'*b = c   and   deg(r) < deg(b).                          */
/*                                                                    */
/* Direct port of IntegrateRational.m:1203-1210.                       */
/* ------------------------------------------------------------------ */

static bool intrat_extended_euclidean(Expr* a, Expr* b, Expr* c, Expr* x,
                                      Expr** out_r, Expr** out_sp) {
    *out_r = *out_sp = NULL;

    Expr* eg_args[3] = { expr_copy(a), expr_copy(b), expr_copy(x) };
    Expr* eg = internal_polynomialextendedgcd(eg_args, 3);
    Expr* eg_eval = evaluate(eg);
    expr_free(eg);
    if (!eg_eval || eg_eval->type != EXPR_FUNCTION
        || eg_eval->data.function.arg_count != 2) {
        if (eg_eval) expr_free(eg_eval);
        return false;
    }
    Expr* g = eg_eval->data.function.args[0];
    Expr* coefs = eg_eval->data.function.args[1];
    if (coefs->type != EXPR_FUNCTION
        || coefs->data.function.head->type != EXPR_SYMBOL
        || coefs->data.function.head->data.symbol.name != SYM_List
        || coefs->data.function.arg_count != 2) {
        expr_free(eg_eval);
        return false;
    }

    Expr* g_copy = expr_copy(g);
    Expr* t_coef = expr_copy(coefs->data.function.args[0]);
    Expr* s_coef = expr_copy(coefs->data.function.args[1]);
    expr_free(eg_eval);

    /* q = c / g via PolynomialQuotient (g divides c). */
    Expr* q = intrat_polyq(c, g_copy, x);
    expr_free(g_copy);
    if (!q) {
        expr_free(t_coef); expr_free(s_coef);
        return false;
    }

    /* t = q*t_coef, s = q*s_coef, both expanded. */
    Expr* qt_raw = internal_times((Expr*[]){expr_copy(q), t_coef}, 2);
    Expr* t_new  = expr_expand(qt_raw); expr_free(qt_raw);
    Expr* qs_raw = internal_times((Expr*[]){expr_copy(q), s_coef}, 2);
    Expr* s_new  = expr_expand(qs_raw); expr_free(qs_raw);
    expr_free(q);

    /* {q2, r} = PolynomialQuotientRemainder[t_new, b, x] */
    Expr *q2 = NULL, *r = NULL;
    if (!intrat_polyqr(t_new, b, x, &q2, &r)) {
        expr_free(t_new); expr_free(s_new);
        return false;
    }
    expr_free(t_new);

    /* sp = Together[s_new + q2*a]. */
    Expr* q2a = internal_times((Expr*[]){q2, expr_copy(a)}, 2);
    Expr* sum = internal_plus((Expr*[]){s_new, q2a}, 2);
    Expr* sp = internal_together((Expr*[]){sum}, 1);

    *out_r = r;
    *out_sp = sp;
    return true;
}

/* AXIOM-style primitivePart(p, q) generalised for the LRT pipeline.
 * Mirrors SUBRESP.primitivePart in intrf.spad.pamphlet:48-52:
 *   r := extEuc(lc(p, var_p), q, 1, var_q).coef1     // cofactor mod q
 *   if r carries unsimplified radicals, Simplify it (otherwise the
 *     downstream PolynomialRemainder step on r*p over an algebraic
 *     extension fails to terminate -- the radical-coefficient
 *     coefficients don't reduce in the underlying algebraic field
 *     and the pseudo-division loop diverges)
 *   rem := PolynomialRemainder(r * p, q, var_q)
 *   return PrimitivePart(rem, var_p)
 *
 * Returns owned Expr* on success, NULL on extEuc / division failure.
 * Borrows p, q, var_p, var_q (caller retains ownership of all four).
 *
 * This consolidates two parallel patches in intrat_int_rational_log_part
 * (the i = deg_d case and the general i case).  Applying the radical-
 * Simplify guard uniformly is defensive -- costs nothing on rational-
 * coefficient inputs because intsimp_has_radical short-circuits. */
static Expr* intrat_primitive_part_mod(Expr* p, Expr* q,
                                       Expr* var_p, Expr* var_q) {
    Expr* lc_p = intrat_lc(p, var_p);
    Expr* one = expr_new_integer(1);
    Expr *r = NULL, *_unused = NULL;
    bool ok = intrat_extended_euclidean(lc_p, q, one, var_q, &r, &_unused);
    expr_free(lc_p);
    expr_free(one);
    if (!ok) return NULL;
    if (_unused) expr_free(_unused);

    r = intsimp_simplify_if_radical(r);

    Expr* rp_raw = internal_times(
        (Expr*[]){r, expr_copy(p)}, 2);
    Expr* rp = expr_expand(rp_raw); expr_free(rp_raw);
    Expr* rem = intrat_polyrem_t(rp, q, var_q);
    expr_free(rp);
    if (!rem) return NULL;
    /* S is the argument of a Log, so it is only defined up to a constant
     * multiple.  We take the integer-content-primitive part (rather than the
     * var_p-monic form): clearing rational denominators in the K(var_q)
     * coefficient field keeps the downstream real-log expansion (LogToReal /
     * LogToAtan over the residue-root algebraic field) on integer-coefficient
     * polynomials, which is markedly faster and avoids the fractional-coefficient
     * blow-up that stalls high-degree cases such as 1/(x^5 + 1). */
    Expr* out = intrat_primitive(rem, var_p);
    expr_free(rem);
    return out;
}

/* ------------------------------------------------------------------ */
/* HermiteReduce (Mack's linear version).                              */
/*                                                                    */
/* Direct port of IntegrateRational.m:1303-1323.  Returns {g, h} with  */
/* f == D[g, x] + h and Denominator[h] squarefree.                     */
/* ------------------------------------------------------------------ */

static Expr* intrat_hermite_reduce(Expr* f, Expr* x) {
    intrat_trace("HermiteReduce", "IN", f);

    /* Numerator / Denominator may return unexpanded results — every
     * downstream `get_degree_poly` call needs the polynomial in
     * expanded form to see the leading term, so we expand once
     * upfront. */
    Expr* a_raw = intrat_numerator(f);
    Expr* a = expr_expand(a_raw); expr_free(a_raw);
    Expr* d_raw = intrat_denominator(f);
    Expr* d = expr_expand(d_raw); expr_free(d_raw);

    /* Over Q: the native (FLINT) reduction in the selected variant
     * (Hermite -> "Linear" | "Quadratic" | "HorowitzOstrogradsky").  The
     * Expr version below (Mack's linear) serves every other coefficient
     * domain. */
    {
        IntratMethodOptions o = ir_opts();
        int hv = o.hermite >= 0 ? o.hermite : INTRAT_DEFAULT_HERMITE;
        Expr* nat = intrat_alg_hermite(a, d, x, hv);
        if (nat) {
            expr_free(a); expr_free(d);
            intrat_trace("HermiteReduce", "OUT", nat);
            return nat;
        }
    }
    Expr* g = expr_new_integer(0);

    /* dbar = gcd(d, D[d, x]) — monic gcd via PolynomialExtendedGCD. */
    Expr* dprime = intrat_d(d, x);
    Expr* dprime_e = expr_expand(dprime); expr_free(dprime);
    Expr* dbar = intrat_polygcd_monic(d, dprime_e, x);
    expr_free(dprime_e);
    if (!dbar) { expr_free(a); expr_free(d); expr_free(g); return NULL; }
    /* dstar = d / dbar.  intrat_exquo is NULL when PolynomialQuotientRemainder
     * did not come back as a pair -- under an abort, a TimeConstrained timeout
     * mid-reduction -- and each NULL below was handed on to internal_divide,
     * which dereferenced it (integrate_risch_transcendental_tests, SIGSEGV in
     * builtin_divide under load, 2026-11-20).  Decline instead. */
    Expr* dstar = intrat_exquo(d, dbar, x);
    expr_free(d);
    if (!dstar) { expr_free(a); expr_free(dbar); expr_free(g); return NULL; }

    while (true) {
        int deg_dbar = get_degree_poly(dbar, x);
        if (deg_dbar <= 0) break;

        /* dbartwo = monic-gcd(dbar, D[dbar, x]) */
        Expr* dbar_prime = intrat_d(dbar, x);
        Expr* dbar_prime_e = expr_expand(dbar_prime); expr_free(dbar_prime);
        Expr* dbartwo = intrat_polygcd_monic(dbar, dbar_prime_e, x);
        expr_free(dbar_prime_e);
        if (!dbartwo) {
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }

        /* dbarhat = dbar / dbartwo */
        Expr* dbarhat = intrat_exquo(dbar, dbartwo, x);
        if (!dbarhat) {
            expr_free(dbartwo); expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }

        /* dtil = -exquo(dstar * D[dbar, x], dbar, x) */
        Expr* dbar_prime2 = intrat_d(dbar, x);
        Expr* dt_raw = internal_times((Expr*[]){expr_copy(dstar), dbar_prime2}, 2);
        Expr* numerator_dtil = expr_expand(dt_raw); expr_free(dt_raw);
        Expr* exq = intrat_exquo(numerator_dtil, dbar, x);
        expr_free(numerator_dtil);
        if (!exq) {
            expr_free(dbarhat); expr_free(dbartwo);
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }
        Expr* dtil_raw = internal_times((Expr*[]){expr_new_integer(-1), exq}, 2);
        Expr* dtil = expr_expand(dtil_raw); expr_free(dtil_raw);

        /* {b_coef, c_coef} = ExtendedEuclidean(dtil, dbarhat, a, x) */
        Expr *b_coef = NULL, *c_coef = NULL;
        if (!intrat_extended_euclidean(dtil, dbarhat, a, x, &b_coef, &c_coef)) {
            expr_free(dtil); expr_free(dbarhat); expr_free(dbartwo);
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }

        /* a := c - exquo(D[b, x] * dstar, dbarhat, x) */
        Expr* b_prime = intrat_d(b_coef, x);
        Expr* bp_raw = internal_times((Expr*[]){b_prime, expr_copy(dstar)}, 2);
        Expr* bp_dstar = expr_expand(bp_raw); expr_free(bp_raw);
        Expr* corr = intrat_exquo(bp_dstar, dbarhat, x);
        expr_free(bp_dstar);
        if (!corr) {
            expr_free(b_coef); expr_free(c_coef); expr_free(dtil); expr_free(dbarhat); expr_free(dbartwo);
            expr_free(a); expr_free(dbar); expr_free(dstar); expr_free(g);
            return NULL;
        }
        Expr* neg_corr = internal_times(
            (Expr*[]){expr_new_integer(-1), corr}, 2);
        Expr* new_a = internal_plus((Expr*[]){c_coef, neg_corr}, 2);
        new_a = eval_and_free(new_a);
        expr_free(a);
        a = new_a;

        /* g := g + canonic(b/dbar) */
        Expr* bd = internal_divide(
            (Expr*[]){b_coef, expr_copy(dbar)}, 2);
        Expr* bd_can = intrat_canonic(bd);
        expr_free(bd);
        Expr* new_g = internal_plus((Expr*[]){g, bd_can}, 2);
        new_g = eval_and_free(new_g);
        g = new_g;

        /* dbar := dbartwo */
        expr_free(dbar); expr_free(dbarhat); expr_free(dtil);
        dbar = dbartwo;
    }
    expr_free(dbar);

    /* h = canonic(a/dstar) */
    Expr* h_raw = internal_divide(
        (Expr*[]){a, dstar}, 2);
    Expr* h = intrat_canonic(h_raw);
    expr_free(h_raw);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){g, h}, 2);
    intrat_trace("HermiteReduce", "OUT", result);
    return result;
}

/* ------------------------------------------------------------------ */
/* IntegratePolynomial — IntegrateRational.m:1153-1162.                */
/* ------------------------------------------------------------------ */

static Expr* intrat_integrate_polynomial(Expr* poly, Expr* x) {
    intrat_trace("IntegratePolynomial", "IN", poly);

    /* CoefficientList[Expand[poly], x] gives {a_0, a_1, ..., a_n}. */
    Expr* expanded = expr_expand(poly);
    Expr* cl_call = internal_coefficientlist(
        (Expr*[]){expanded, expr_copy(x)}, 2);
    Expr* cl = evaluate(cl_call);
    expr_free(cl_call);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        /* Fall back to the trivial constant case: a polynomial that is
         * a pure constant has CoefficientList -> {const}. */
        return NULL;
    }
    size_t n = cl->data.function.arg_count;
    if (n == 0) { expr_free(cl); return expr_new_integer(0); }

    /* Build Σ a_k x^(k+1) / (k+1). */
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        Expr* a = cl->data.function.args[k];
        Expr* x_pow = internal_power(
            (Expr*[]){expr_copy(x), expr_new_integer((int64_t)(k + 1))}, 2);
        Expr* div = internal_divide(
            (Expr*[]){x_pow, expr_new_integer((int64_t)(k + 1))}, 2);
        terms[k] = internal_times(
            (Expr*[]){expr_copy(a), div}, 2);
    }
    expr_free(cl);
    Expr* sum = internal_plus(terms, n);
    free(terms);
    Expr* sum_eval = eval_and_free(sum);

    intrat_trace("IntegratePolynomial", "OUT", sum_eval);
    return sum_eval;
}

/* ------------------------------------------------------------------ */
/* Pre-Hermite derivative-recognition fast path.                        */
/*                                                                    */
/* If f = num/den with FactorSquareFree[den] = c0 * pol^k for a single */
/* squarefree pol, and num == c * D[pol, x] for some constant c, then  */
/*   k == 1: ∫ f dx = (c/c0) Log[pol]                                  */
/*   k >= 2: ∫ f dx = -(c/c0) / ((k-1) * pol^(k-1))                    */
/*                                                                    */
/* Returns NULL if the structure does not match — the caller falls     */
/* back to the standard Hermite + LRT pipeline.                        */
/* ------------------------------------------------------------------ */

/* Walk a FactorSquareFree[Times[..., pol^k, ...]] result and decide
 * whether it has the shape `c0 * pol^k` for a single non-constant
 * squarefree polynomial `pol` in `x`.  Returns true on match and
 * fills `*c0_out`, `*pol_out`, `*k_out`. */
static bool extract_factored_denominator(Expr* fs, Expr* x,
                                         Expr** c0_out, Expr** pol_out,
                                         int* k_out) {
    Expr* fs_expanded = fs;
    /* Times-headed: fs = c0 * factor1 * factor2 * ... */
    if (fs_expanded->type == EXPR_FUNCTION
        && fs_expanded->data.function.head->type == EXPR_SYMBOL
        && fs_expanded->data.function.head->data.symbol.name == SYM_Times) {
        Expr* c0 = expr_new_integer(1);
        Expr* pol = NULL;
        int k = 0;
        for (size_t i = 0; i < fs_expanded->data.function.arg_count; i++) {
            Expr* arg = fs_expanded->data.function.args[i];
            if (intrat_freeq_test(arg, x)) {
                Expr* prod = internal_times(
                    (Expr*[]){c0, expr_copy(arg)}, 2);
                c0 = eval_and_free(prod);
                continue;
            }
            /* Non-constant factor: must be unique and Power[pol, k] form. */
            if (pol) { expr_free(c0); return false; }
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Power
                && arg->data.function.arg_count == 2
                && arg->data.function.args[1]->type == EXPR_INTEGER
                && arg->data.function.args[1]->data.integer >= 1) {
                pol = expr_copy(arg->data.function.args[0]);
                k = (int)arg->data.function.args[1]->data.integer;
            } else {
                pol = expr_copy(arg);
                k = 1;
            }
        }
        if (!pol) { expr_free(c0); return false; }
        *c0_out = c0; *pol_out = pol; *k_out = k;
        return true;
    }
    /* Single Power[pol, k]. */
    if (fs_expanded->type == EXPR_FUNCTION
        && fs_expanded->data.function.head->type == EXPR_SYMBOL
        && fs_expanded->data.function.head->data.symbol.name == SYM_Power
        && fs_expanded->data.function.arg_count == 2
        && fs_expanded->data.function.args[1]->type == EXPR_INTEGER
        && fs_expanded->data.function.args[1]->data.integer >= 1
        && !intrat_freeq_test(fs_expanded->data.function.args[0], x)) {
        *c0_out = expr_new_integer(1);
        *pol_out = expr_copy(fs_expanded->data.function.args[0]);
        *k_out = (int)fs_expanded->data.function.args[1]->data.integer;
        return true;
    }
    /* Bare polynomial. */
    if (!intrat_freeq_test(fs_expanded, x)) {
        *c0_out = expr_new_integer(1);
        *pol_out = expr_copy(fs_expanded);
        *k_out = 1;
        return true;
    }
    return false;
}

/* Returns the closed-form integral if the derivative-recognition
 * pattern matches; NULL otherwise. */
static Expr* intrat_derivative_recognition(Expr* num, Expr* den, Expr* x) {
    /* Factor the denominator into (constant) * (single squarefree pol)^k. */
    Expr* fs = internal_factorsquarefree(
        (Expr*[]){expr_copy(den)}, 1);
    Expr* fs_eval = evaluate(fs);
    expr_free(fs);
    if (!fs_eval) return NULL;

    Expr *c0 = NULL, *pol = NULL;
    int k = 0;
    bool ok = extract_factored_denominator(fs_eval, x, &c0, &pol, &k);
    expr_free(fs_eval);
    if (!ok) return NULL;

    /* Compute c = num / (c0 * D[pol, x]).  If FreeQ[c, x] we have the
     * derivative-recognition form. */
    Expr* dpol = intrat_d(pol, x);
    Expr* c0_dpol_raw = internal_times((Expr*[]){expr_copy(c0), dpol}, 2);
    Expr* c0_dpol = expr_expand(c0_dpol_raw); expr_free(c0_dpol_raw);
    Expr* ratio = internal_divide(
        (Expr*[]){expr_copy(num), c0_dpol}, 2);
    Expr* c_const = intrat_canonic(ratio);
    expr_free(ratio);

    if (!intrat_freeq_test(c_const, x)) {
        expr_free(c0); expr_free(pol); expr_free(c_const);
        return NULL;
    }

    Expr* result;
    if (k == 1) {
        /* ∫ c * D[pol]/pol dx = c * Log[pol] */
        Expr* log_pol = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){pol}, 1);
        result = internal_times(
            (Expr*[]){c_const, log_pol}, 2);
    } else {
        /* ∫ c * D[pol]/pol^k dx = -c / ((k-1) * pol^(k-1)) */
        Expr* pol_pow = internal_power(
            (Expr*[]){pol, expr_new_integer((int64_t)(k - 1))}, 2);
        Expr* denom = internal_times(
            (Expr*[]){expr_new_integer((int64_t)(k - 1)), pol_pow}, 2);
        Expr* div = internal_divide(
            (Expr*[]){c_const, denom}, 2);
        result = internal_times(
            (Expr*[]){expr_new_integer(-1), div}, 2);
    }
    expr_free(c0);
    return eval_and_free(result);
}

/* ====================================================================
 * Phase 2 — Lazard-Rioboo-Trager log part
 * ====================================================================
 *
 * Implements `IntRationalLogPart[A/D, x, t]` (Bronstein, Symbolic
 * Integration I, p. 51) for square-free denominators D in K[x],
 * along with the supporting helpers SquareFree, ExtractConstants,
 * ApartList, and a structural pre-detector for the cyclotomic /
 * n-th-root denominator family that appears repeatedly in Bronstein
 * and the BronsteinRational test suite.
 *
 * Output form (RootSum -> False, the default for Phase 2's structural
 * tests):
 *   list of {Q_i(t), S_i(t, x)} pairs encoding
 *      ∫ A/D dx == Σ_i RootSum[Q_i(t)==0, t Log[S_i(t, x)]].
 *
 * Phase 4's LogToReal will consume the same (Q_i, S_i) pairs and
 * convert each one into real Log + ArcTan form when Q_i factors
 * tractably.  Until then, the Phase 2 builtin returns the symbolic
 * RootSum sum (RootSum -> True) or the raw pair list (RootSum ->
 * False).
 */

/* ----- SquareFree: {poly, multiplicity} list ----- */

/* Bucket FactorSquareFree[p]'s output into a list of (poly,
 * multiplicity) pairs, indexed densely by multiplicity from 1..max.
 * If the input has no factor of multiplicity i, the i-th entry is
 * {1, i}.  Mirrors the BronsteinRational baseline's `SquareFree[]`
 * helper at IntegrateRational.m:1474-1487.
 *
 * Implementation note: Mathilda's FactorSquareFree returns the
 * factorisation as `Times[c, Power[p1, k1], Power[p2, k2], ...]`
 * (or a single Power[]/atomic factor).  We walk that structure,
 * grouping factors by multiplicity. */
/* True iff e carries a radical constant: Power[c, r] with c a number (no
 * symbol) and r a non-integer rational -- Sqrt[2], 3^(1/3), ... */
static bool intrat_has_radical_constant(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* r = e->data.function.args[1];
        bool num_base = b->type == EXPR_INTEGER || b->type == EXPR_BIGINT ||
            (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL &&
             b->data.function.head->data.symbol.name == SYM_Rational);
        bool frac_exp = r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
             r->data.function.head->data.symbol.name == SYM_Rational;
        if (num_base && frac_exp) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_has_radical_constant(e->data.function.args[i])) return true;
    return false;
}

static Expr* intrat_squarefree_list(Expr* p) {
    /* Run FactorSquareFree first.  Over a number field (radical constants in
     * the coefficients) the plain call treats Sqrt[2] as an independent
     * indeterminate and misses repeated factors: the Rothstein-Trager
     * resultant (99 - 70 Sqrt[2]) (16 t^2 - 1)^2 of
     * (x^2 - 3 + 2 Sqrt[2])/(x^4 + (2 - 4 Sqrt[2]) x^2 + 17 - 12 Sqrt[2])
     * came back as one "squarefree" quartic, and LogToReal divided by its zero
     * discriminant (Integrate answered Indeterminate).  Extension -> Automatic
     * factors over the field the coefficients generate. */
    Expr* fs;
    if (intrat_has_radical_constant(p)) {
        Expr* opt = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol("Extension"), expr_new_symbol(SYM_Automatic) }, 2);
        fs = expr_new_function(expr_new_symbol("FactorSquareFree"),
            (Expr*[]){ expr_copy(p), opt }, 2);
    } else {
        fs = internal_factorsquarefree((Expr*[]){expr_copy(p)}, 1);
    }
    Expr* fs_eval = evaluate(fs);
    expr_free(fs);
    if (!fs_eval) return NULL;

    /* Walk the factor structure, accumulating (factor, mult) pairs. */
    typedef struct { Expr* poly; int mult; } Pair;
    Pair* pairs = (Pair*)malloc(sizeof(Pair) * 8);
    size_t pcap = 8, pcount = 0;

#define PUSH_PAIR(POLY, MULT) do { \
        if (pcount >= pcap) { pcap *= 2; \
            pairs = (Pair*)realloc(pairs, sizeof(Pair) * pcap); } \
        pairs[pcount].poly = (POLY); \
        pairs[pcount].mult = (MULT); \
        pcount++; \
    } while (0)

    if (fs_eval->type == EXPR_FUNCTION
        && fs_eval->data.function.head->type == EXPR_SYMBOL
        && fs_eval->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < fs_eval->data.function.arg_count; i++) {
            Expr* arg = fs_eval->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Power
                && arg->data.function.arg_count == 2
                && arg->data.function.args[1]->type == EXPR_INTEGER) {
                int k = (int)arg->data.function.args[1]->data.integer;
                if (k >= 1) {
                    PUSH_PAIR(expr_copy(arg->data.function.args[0]), k);
                    continue;
                }
            }
            PUSH_PAIR(expr_copy(arg), 1);
        }
    } else if (fs_eval->type == EXPR_FUNCTION
        && fs_eval->data.function.head->type == EXPR_SYMBOL
        && fs_eval->data.function.head->data.symbol.name == SYM_Power
        && fs_eval->data.function.arg_count == 2
        && fs_eval->data.function.args[1]->type == EXPR_INTEGER) {
        int k = (int)fs_eval->data.function.args[1]->data.integer;
        if (k >= 1) PUSH_PAIR(expr_copy(fs_eval->data.function.args[0]), k);
        else PUSH_PAIR(expr_copy(fs_eval), 1);
    } else {
        /* Single bare factor — multiplicity 1. */
        PUSH_PAIR(expr_copy(fs_eval), 1);
    }
    expr_free(fs_eval);

#undef PUSH_PAIR

    /* Determine max multiplicity. */
    int max_mult = 0;
    for (size_t i = 0; i < pcount; i++) {
        if (pairs[i].mult > max_mult) max_mult = pairs[i].mult;
    }
    if (max_mult == 0) {
        /* Constant input — return {{c, 1}}. */
        Expr* pair = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){expr_copy(p), expr_new_integer(1)}, 2);
        Expr* list = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){pair}, 1);
        for (size_t i = 0; i < pcount; i++) expr_free(pairs[i].poly);
        free(pairs);
        return list;
    }

    /* Build output list of length max_mult.  Index i (1-based) has
     * the product of all factors of multiplicity i, or 1 if none. */
    Expr** out_args = (Expr**)malloc(sizeof(Expr*) * max_mult);
    for (int m = 1; m <= max_mult; m++) {
        /* Collect all pairs with this multiplicity. */
        Expr** factors = NULL;
        size_t fcount = 0;
        for (size_t i = 0; i < pcount; i++) {
            if (pairs[i].mult == m) {
                factors = (Expr**)realloc(factors, sizeof(Expr*) * (fcount + 1));
                factors[fcount++] = expr_copy(pairs[i].poly);
            }
        }
        Expr* poly_for_m;
        if (fcount == 0) {
            poly_for_m = expr_new_integer(1);
        } else if (fcount == 1) {
            poly_for_m = expr_expand(factors[0]);
            expr_free(factors[0]);
        } else {
            Expr* prod = internal_times(factors, fcount);
            poly_for_m = expr_expand(prod);
            expr_free(prod);
        }
        if (factors) free(factors);
        out_args[m - 1] = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){poly_for_m, expr_new_integer(m)}, 2);
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), out_args, max_mult);
    free(out_args);

    for (size_t i = 0; i < pcount; i++) expr_free(pairs[i].poly);
    free(pairs);
    return list;
}

/* Extract the i-th element of a Mathematica-style List expression.
 * Returns a fresh copy.  i is 1-based to match Mathematica.  Returns
 * NULL if `lst` is not a List or i is out of range. */
static Expr* list_get(Expr* lst, size_t i) {
    if (!lst || lst->type != EXPR_FUNCTION
        || lst->data.function.head->type != EXPR_SYMBOL
        || lst->data.function.head->data.symbol.name != SYM_List
        || i == 0 || i > lst->data.function.arg_count) return NULL;
    return expr_copy(lst->data.function.args[i - 1]);
}

/* ----- ExtractConstants: pull constant prefactors out of f ----- */

/* Split p into (numeric prefactor, x-bearing residue).  Walk a Times
 * head; every factor that is FreeQ[_, x] absorbs into the prefactor,
 * everything else stays in the residue.  When the input is a single
 * factor we either return (1, p) or (p, 1).  Mirrors the part of
 * IntegrateRational.m:1013-1023 that wraps FactorSquareFreeList +
 * Power@@first, but we work directly on Times args to avoid an extra
 * factorisation pass. */
static void intrat_split_constant(Expr* p, Expr* x,
                                  Expr** out_const, Expr** out_residue) {
    if (intrat_freeq_test(p, x)) {
        *out_const = expr_copy(p);
        *out_residue = expr_new_integer(1);
        return;
    }
    if (p->type != EXPR_FUNCTION
        || p->data.function.head->type != EXPR_SYMBOL
        || p->data.function.head->data.symbol.name != SYM_Times) {
        *out_const = expr_new_integer(1);
        *out_residue = expr_copy(p);
        return;
    }
    /* Times-headed: split args. */
    size_t n = p->data.function.arg_count;
    Expr** consts = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** resid  = (Expr**)malloc(sizeof(Expr*) * n);
    size_t cn = 0, rn = 0;
    for (size_t i = 0; i < n; i++) {
        if (intrat_freeq_test(p->data.function.args[i], x)) {
            consts[cn++] = expr_copy(p->data.function.args[i]);
        } else {
            resid[rn++] = expr_copy(p->data.function.args[i]);
        }
    }
    *out_const  = (cn == 0) ? expr_new_integer(1)
                : (cn == 1) ? consts[0]
                : eval_and_free(internal_times(consts, cn));
    *out_residue = (rn == 0) ? expr_new_integer(1)
                : (rn == 1) ? resid[0]
                : eval_and_free(internal_times(resid, rn));
    free(consts); free(resid);
}

/* extractConstants[f, x] = {numconst/denconst, num_residue/den_residue}.
 * Returns a List[const, simplified] expression. */
static Expr* intrat_extract_constants(Expr* f, Expr* x) {
    Expr* num = intrat_numerator(f);
    Expr* den = intrat_denominator(f);
    Expr* num_e = eval_and_free(num);
    Expr* den_e = eval_and_free(den);
    /* Apply FactorSquareFree once to surface the constant prefactor. */
    Expr* num_fs = internal_factorsquarefree((Expr*[]){num_e}, 1);
    num_fs = eval_and_free(num_fs);
    Expr* den_fs = internal_factorsquarefree((Expr*[]){den_e}, 1);
    den_fs = eval_and_free(den_fs);

    Expr *nc = NULL, *nr = NULL, *dc = NULL, *dr = NULL;
    intrat_split_constant(num_fs, x, &nc, &nr);
    intrat_split_constant(den_fs, x, &dc, &dr);
    expr_free(num_fs); expr_free(den_fs);

    Expr* const_part = internal_divide(
        (Expr*[]){nc, dc}, 2);
    const_part = eval_and_free(const_part);
    Expr* residue = internal_divide(
        (Expr*[]){nr, dr}, 2);
    residue = eval_and_free(residue);

    return expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){const_part, residue}, 2);
}

/* ----- ApartList: partial-fraction expansion as a list ----- */

/* apartList[f, x] = list of partial-fraction terms.  Phase 2 keeps
 * the lightweight implementation: call Apart[f, x], and split the
 * resulting Plus into a list of summands.  When Apart returns a
 * single non-Plus term we wrap it in a singleton list. */
static Expr* intrat_apart_list(Expr* f, Expr* x, const Expr* alpha) {
    /* Build Apart[f, x, Extension -> alpha?]. */
    size_t argc = (alpha ? 3 : 2);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * argc);
    args[0] = expr_copy(f);
    args[1] = expr_copy(x);
    if (alpha) {
        args[2] = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){expr_new_symbol(SYM_Extension), expr_copy((Expr*)alpha)}, 2);
    }
    Expr* call = expr_new_function(expr_new_symbol(SYM_Apart), args, argc);
    Expr* ap = evaluate(call);
    expr_free(call);
    free(args);

    if (!ap) return NULL;

    if (ap->type == EXPR_FUNCTION
        && ap->data.function.head->type == EXPR_SYMBOL
        && ap->data.function.head->data.symbol.name == SYM_Plus) {
        size_t n = ap->data.function.arg_count;
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) terms[i] = expr_copy(ap->data.function.args[i]);
        expr_free(ap);
        /* expr_new_function memcpy's the args buffer into its own
         * backing store, so the caller's `terms` buffer must be freed. */
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), terms, n);
        free(terms);
        return list;
    }
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ap}, 1);
}

/* ----- Cyclotomic / nth-root pre-detection ----- */

/* nthrootQ[p, x]: True iff p has the form `a x^n + b` with a, b free
 * of x.  Caught by Mathematica's MatchQ pattern at
 * IntegrateRational.m:1093.  The integer n is returned via *n_out.
 * Reserved for Phase 4's LogToReal closer (n-th-root-of-unity case);
 * unused in Phase 2 but kept here so the helper does not need to be
 * resurrected when Phase 4 lands. */
MATHILDA_MAYBE_UNUSED
static bool intrat_nthroot_q(Expr* p, Expr* x, int* n_out, Expr** a_out, Expr** b_out) {
    /* Run Collect[p, x] to gather the polynomial in canonical form. */
    Expr* coll = expr_new_function(expr_new_symbol(SYM_Collect),
        (Expr*[]){expr_copy(p), expr_copy(x)}, 2);
    Expr* c = evaluate(coll);
    expr_free(coll);
    if (!c) return false;

    /* Walk a Plus-headed result. */
    Expr* a = NULL;
    Expr* b = NULL;
    int n = 0;

    Expr** terms;
    size_t nterms;
    if (c->type == EXPR_FUNCTION
        && c->data.function.head->type == EXPR_SYMBOL
        && c->data.function.head->data.symbol.name == SYM_Plus) {
        terms = c->data.function.args;
        nterms = c->data.function.arg_count;
    } else {
        terms = &c;
        nterms = 1;
    }

    for (size_t i = 0; i < nterms; i++) {
        Expr* term = terms[i];
        if (intrat_freeq_test(term, x)) {
            if (b) { expr_free(b); expr_free(c); if (a) expr_free(a); return false; }
            b = expr_copy(term);
            continue;
        }
        /* term must be c * x^n with c free of x (or x^n alone, c=1). */
        Expr* coeff = NULL;
        int deg = 0;
        if (term->type == EXPR_FUNCTION
            && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol.name == SYM_Power
            && term->data.function.arg_count == 2
            && expr_eq(term->data.function.args[0], x)
            && term->data.function.args[1]->type == EXPR_INTEGER) {
            deg = (int)term->data.function.args[1]->data.integer;
            coeff = expr_new_integer(1);
        } else if (expr_eq(term, x)) {
            deg = 1; coeff = expr_new_integer(1);
        } else if (term->type == EXPR_FUNCTION
            && term->data.function.head->type == EXPR_SYMBOL
            && term->data.function.head->data.symbol.name == SYM_Times) {
            /* coeff * x^n form. */
            Expr** factors = term->data.function.args;
            size_t nf = term->data.function.arg_count;
            Expr** ccoll = (Expr**)malloc(sizeof(Expr*) * nf);
            size_t ccn = 0;
            for (size_t k = 0; k < nf; k++) {
                Expr* f = factors[k];
                if (intrat_freeq_test(f, x)) {
                    ccoll[ccn++] = expr_copy(f);
                } else if (deg == 0 && f->type == EXPR_FUNCTION
                    && f->data.function.head->type == EXPR_SYMBOL
                    && f->data.function.head->data.symbol.name == SYM_Power
                    && f->data.function.arg_count == 2
                    && expr_eq(f->data.function.args[0], x)
                    && f->data.function.args[1]->type == EXPR_INTEGER) {
                    deg = (int)f->data.function.args[1]->data.integer;
                } else if (deg == 0 && expr_eq(f, x)) {
                    deg = 1;
                } else {
                    deg = -1; break;
                }
            }
            if (deg <= 0) {
                for (size_t kk = 0; kk < ccn; kk++) expr_free(ccoll[kk]);
                free(ccoll);
                if (a) expr_free(a);
                if (b) expr_free(b);
                expr_free(c);
                return false;
            }
            coeff = (ccn == 0) ? expr_new_integer(1)
                  : (ccn == 1) ? ccoll[0]
                  : eval_and_free(internal_times(ccoll, ccn));
            free(ccoll);
        } else {
            if (a) expr_free(a);
            if (b) expr_free(b);
            expr_free(c);
            return false;
        }
        if (a) { expr_free(a); expr_free(b ? b : NULL); expr_free(c);
                 expr_free(coeff); return false; }
        a = coeff;
        n = deg;
    }

    expr_free(c);
    if (!a || n <= 1) {
        if (a) expr_free(a);
        if (b) expr_free(b);
        return false;
    }
    *n_out = n;
    *a_out = a;
    *b_out = b ? b : expr_new_integer(0);
    return true;
}

/* ----- IntRationalLogPart: the LRT log-part computation ----- */

/* Find the unique element of the prs chain whose degree in x equals
 * `target_deg`.  Returns a fresh copy or NULL if none exists. */
static Expr* find_prs_at_degree(Expr* prs, Expr* x, int target_deg) {
    if (prs->type != EXPR_FUNCTION
        || prs->data.function.head->type != EXPR_SYMBOL
        || prs->data.function.head->data.symbol.name != SYM_List) return NULL;
    for (size_t i = 0; i < prs->data.function.arg_count; i++) {
        Expr* el = prs->data.function.args[i];
        Expr* el_e = expr_expand(el);
        int d = get_degree_poly(el_e, x);
        if (d == target_deg) {
            return el_e;
        }
        expr_free(el_e);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Simple-root expansion (Phase 8d-bonus).                             */
/* ------------------------------------------------------------------ */

/* Tiny tree builders used only by the radical-form root extraction.
 * Each takes ownership of every Expr* arg and returns a fresh tree
 * the caller owns.  All paths route through eval_and_free at the
 * end so the outer Plus/Times can fold trivial identities. */
static Expr* nlp_neg(Expr* a) {
    return internal_times((Expr*[]){ expr_new_integer(-1), a }, 2);
}
static Expr* nlp_sqrt(Expr* a) {
    /* Power[a, 1/2] — Mathilda's Power evaluator turns this into the
     * canonical Sqrt[..] surface form during printing. */
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ a,
                   expr_new_function(expr_new_symbol(SYM_Rational),
                       (Expr*[]){ expr_new_integer(1),
                                  expr_new_integer(2) }, 2) },
        2);
}
/* Forward declarations for helpers defined later in this file but
 * referenced by the palindromic-quartic expander immediately below. */
static Expr* subst_t(Expr* e, Expr* t, Expr* sub_expr);
static void split_re_im(Expr* p, Expr** re_out, Expr** im_out);

static Expr* nlp_sub_body(Expr* body, Expr* bvar, Expr* value) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy(bvar), expr_copy(value) }, 2);
    Expr* sub = internal_replace_all(
        (Expr*[]){ expr_copy(body), rule }, 2);
    return eval_and_free(sub);
}

/* True iff every variable in `poly` is structurally equal to `bvar`.
 * Phase 8d-bonus only fires on non-parametric polynomials: when poly
 * has parameters (e.g. b + a t^4), the radical-formula substitution
 * produces forms whose D[..., x] - integrand reduces to 0 only under
 * a stronger algebraic-extension simplifier than Mathilda currently
 * has.  Bailing out for parametric cases preserves the held RootSum
 * form, which the user can verify formally and which avoids
 * apparent-regression noise on the corpus. */
static bool poly_only_uses(Expr* poly, Expr* bvar) {
    Expr* vars_call = expr_new_function(expr_new_symbol(SYM_Variables),
        (Expr*[]){ expr_copy(poly) }, 1);
    Expr* vars = evaluate(vars_call);
    expr_free(vars_call);
    if (!vars || vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) {
        if (vars) expr_free(vars);
        return false;
    }
    bool ok = true;
    for (size_t i = 0; i < vars->data.function.arg_count; i++) {
        if (!expr_eq(vars->data.function.args[i], bvar)) { ok = false; break; }
    }
    expr_free(vars);
    return ok;
}

/* Build the real form of
 *   RootSum[ Function[t, d_t],
 *            Function[t, a_t · Log[x − t] / dd_t] ]
 * for a palindromic quartic d_t = c0 t^4 + c1 t^3 + c2 t^2 + c1 t + c0
 * (c0 == c4, c1 == c3).
 *
 * Strategy: factor d_t over R via the u = t + 1/t substitution
 *   c0 u^2 + c1 u + (c2 − 2 c0) = 0  ⇒  u = u_±,
 * giving two complex-conjugate root pairs of d_t:
 *   α_k = u_k/2 + i v_k,   v_k = Sqrt[(4 − u_k^2)/4].
 * Each conjugate pair contributes
 *   2 Re[ (a(α_k)/d'(α_k)) · Log[x − α_k] ]
 *     = C_re · Log[(x − u_k/2)^2 + v_k^2]
 *     + 2 C_im · ArcTan[v_k/(x − u_k/2)]
 * where C = a(α_k)/d'(α_k) = C_re + i C_im.
 *
 * Avoids LogToReal/LogToAtan entirely (those bog down on the nested-
 * radical S polynomials from the scaled-palindromic LRT path).
 * Returns NULL when d_t isn't palindromic, when u_± aren't real, or
 * when any sign-gate is undecidable. */
static Expr* expand_palindromic_quartic_real(
    Expr* a_t, Expr* d_t, Expr* dd_t, Expr* t, Expr* x)
{
    /* Detect palindromic quartic in t. */
    int deg = get_degree_poly(d_t, t);
    if (deg != 4) return NULL;
    Expr* c0 = get_coeff(d_t, t, 0);
    Expr* c1 = get_coeff(d_t, t, 1);
    Expr* c2 = get_coeff(d_t, t, 2);
    Expr* c3 = get_coeff(d_t, t, 3);
    Expr* c4 = get_coeff(d_t, t, 4);

    Expr* c0mc4 = internal_subtract(
        (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2);
    Expr* c1mc3 = internal_subtract(
        (Expr*[]){expr_copy(c1), expr_copy(c3)}, 2);
    bool palindromic = intsimp_zero_q(c0mc4) && intsimp_zero_q(c1mc3);
    expr_free(c0mc4); expr_free(c1mc3);

    bool c0_nonzero = !(c0->type == EXPR_INTEGER && c0->data.integer == 0);
    if (!palindromic || !c0_nonzero) {
        expr_free(c0); expr_free(c1); expr_free(c2);
        expr_free(c3); expr_free(c4);
        return NULL;
    }

    /* Outer disc Du = c1^2 − 4 c0 (c2 − 2 c0).  Soundness gate:
     * Du >= 0 so u_± are real. */
    Expr* c1sq = expand_and_free(internal_power(
        (Expr*[]){expr_copy(c1), expr_new_integer(2)}, 2));
    Expr* two_c0 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(c0)}, 2));
    Expr* c2m2c0 = eval_and_free(internal_subtract(
        (Expr*[]){expr_copy(c2), two_c0}, 2));
    Expr* four_term = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(4),
                  expr_copy(c0),
                  c2m2c0}, 3));
    Expr* neg_four_term = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), four_term}, 2));
    Expr* Du = eval_and_free(internal_plus(
        (Expr*[]){c1sq, neg_four_term}, 2));
    Expr* Du_exp = expr_expand(Du); expr_free(Du);
    int Du_sign = intsimp_sign_pos_assumption(Du_exp);
    if (Du_sign == 0) Du_sign = intsimp_numeric_sign(Du_exp);
    if (Du_sign < 0) {
        expr_free(Du_exp);
        expr_free(c0); expr_free(c1); expr_free(c2);
        expr_free(c3); expr_free(c4);
        return NULL;
    }
    Expr* sqrt_Du = intsimp_pos_sqrt(Du_exp);
    expr_free(Du_exp);

    Expr* neg_c1 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(c1)}, 2));
    Expr* two_c0_b = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(c0)}, 2));

    Expr* u_num_p = eval_and_free(internal_plus(
        (Expr*[]){expr_copy(neg_c1), expr_copy(sqrt_Du)}, 2));
    Expr* u_num_m = eval_and_free(internal_subtract(
        (Expr*[]){neg_c1, sqrt_Du}, 2));
    Expr* u_plus = eval_and_free(internal_divide(
        (Expr*[]){u_num_p, expr_copy(two_c0_b)}, 2));
    Expr* u_minus = eval_and_free(internal_divide(
        (Expr*[]){u_num_m, two_c0_b}, 2));

    expr_free(c0); expr_free(c1); expr_free(c2);
    expr_free(c3); expr_free(c4);

    /* For each u, compute the conjugate-pair real contribution. */
    Expr* us[2] = {u_plus, u_minus};
    Expr* contributions[2] = {NULL, NULL};
    bool ok = true;

    for (int k = 0; k < 2 && ok; k++) {
        Expr* u = us[k];

        /* v_sq = (4 - u^2)/4.  Must be > 0 for v_im real. */
        Expr* u_sq = expand_and_free(internal_power(
            (Expr*[]){expr_copy(u), expr_new_integer(2)}, 2));
        Expr* neg_u_sq = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), u_sq}, 2));
        Expr* four_m_usq = eval_and_free(internal_plus(
            (Expr*[]){expr_new_integer(4), neg_u_sq}, 2));
        Expr* v_sq = eval_and_free(internal_divide(
            (Expr*[]){four_m_usq, expr_new_integer(4)}, 2));
        Expr* v_sq_exp = expr_expand(v_sq); expr_free(v_sq);

        int v_sq_sign = intsimp_sign_pos_assumption(v_sq_exp);
        if (v_sq_sign == 0) v_sq_sign = intsimp_numeric_sign(v_sq_exp);
        if (v_sq_sign <= 0) {
            expr_free(v_sq_exp);
            ok = false; break;
        }

        Expr* v_im = intsimp_pos_sqrt(v_sq_exp);

        Expr* u_half = eval_and_free(internal_divide(
            (Expr*[]){expr_copy(u), expr_new_integer(2)}, 2));

        /* α = u_half + I · v_im. */
        Expr* I_unit = expr_new_function(expr_new_symbol(SYM_Complex),
            (Expr*[]){expr_new_integer(0), expr_new_integer(1)}, 2);
        Expr* I_v = eval_and_free(internal_times(
            (Expr*[]){I_unit, expr_copy(v_im)}, 2));
        Expr* alpha = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(u_half), I_v}, 2));

        /* a(α) and d'(α) — split each into real + imag. */
        Expr* a_at_alpha = subst_t(a_t, t, alpha);
        Expr* dp_at_alpha = subst_t(dd_t, t, alpha);
        expr_free(alpha);

        Expr *Aa = NULL, *Ba = NULL;
        split_re_im(a_at_alpha, &Aa, &Ba);
        expr_free(a_at_alpha);
        Expr *Ad = NULL, *Bd = NULL;
        split_re_im(dp_at_alpha, &Ad, &Bd);
        expr_free(dp_at_alpha);

        /* c = (Aa + i Ba) / (Ad + i Bd)
         *   = (Aa Ad + Ba Bd + i(Ba Ad − Aa Bd)) / (Ad^2 + Bd^2). */
        Expr* AaAd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Aa), expr_copy(Ad)}, 2));
        Expr* BaBd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Ba), expr_copy(Bd)}, 2));
        Expr* re_num = eval_and_free(internal_plus(
            (Expr*[]){AaAd, BaBd}, 2));
        Expr* BaAd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Ba), expr_copy(Ad)}, 2));
        Expr* AaBd = expand_and_free(internal_times(
            (Expr*[]){expr_copy(Aa), expr_copy(Bd)}, 2));
        Expr* im_num = eval_and_free(internal_subtract(
            (Expr*[]){BaAd, AaBd}, 2));
        Expr* Ad_sq = expand_and_free(internal_power(
            (Expr*[]){Ad, expr_new_integer(2)}, 2));
        Expr* Bd_sq = expand_and_free(internal_power(
            (Expr*[]){Bd, expr_new_integer(2)}, 2));
        Expr* denom = eval_and_free(internal_plus(
            (Expr*[]){Ad_sq, Bd_sq}, 2));
        Expr* C_re = eval_and_free(internal_divide(
            (Expr*[]){re_num, expr_copy(denom)}, 2));
        Expr* C_im = eval_and_free(internal_divide(
            (Expr*[]){im_num, denom}, 2));
        expr_free(Aa); expr_free(Ba);

        /* L = (x − u_half)^2 + v_sq.
         * T = v_im / (x − u_half). */
        Expr* neg_uh = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(u_half)}, 2));
        Expr* x_m_uh = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(x), neg_uh}, 2));
        Expr* x_m_uh_sq = expand_and_free(internal_power(
            (Expr*[]){expr_copy(x_m_uh), expr_new_integer(2)}, 2));
        Expr* L = eval_and_free(internal_plus(
            (Expr*[]){x_m_uh_sq, v_sq_exp}, 2));
        Expr* T = eval_and_free(internal_divide(
            (Expr*[]){expr_copy(v_im), expr_copy(x_m_uh)}, 2));
        expr_free(x_m_uh); expr_free(v_im); expr_free(u_half);

        Expr* log_L = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){L}, 1);
        Expr* atan_T = expr_new_function(expr_new_symbol(SYM_ArcTan),
            (Expr*[]){T}, 1);

        Expr* term_log = eval_and_free(internal_times(
            (Expr*[]){expr_copy(C_re), log_L}, 2));
        Expr* term_atan = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(2),
                      expr_copy(C_im), atan_T}, 3));
        expr_free(C_re); expr_free(C_im);

        contributions[k] = eval_and_free(internal_plus(
            (Expr*[]){term_log, term_atan}, 2));
    }

    expr_free(u_plus); expr_free(u_minus);

    if (!ok) {
        if (contributions[0]) expr_free(contributions[0]);
        if (contributions[1]) expr_free(contributions[1]);
        return NULL;
    }

    Expr* total = eval_and_free(internal_plus(contributions, 2));
    return total;
}

/* When `poly` (a univariate polynomial in `bvar`) is one of the
 * structural shapes whose roots have a closed-form radical formula,
 * expand `RootSum[Function[bvar, poly], Function[bvar, body]]` into
 * an explicit Plus[body[α₁], body[α₂], ...] and return it.  Returns
 * NULL when the polynomial does not match any of the supported
 * shapes — currently:
 *
 *   degree 1                                         (1 root)
 *   degree 2                                         (2 roots)
 *   degree 4 with only constant / x^2 / x^4 terms    (4 roots)
 *
 * Gated on poly being free of parameters other than bvar — see
 * poly_only_uses above for the rationale.
 *
 * Once Solve / ToRadicals lands, this routine is superseded by a
 * call to Solve[poly == 0, bvar, Reals] followed by the same body
 * substitution loop. */
static bool intrat_contains_complex(const Expr* e);

/* True when e contains a RootSum (an unexpanded log sum).  Borrowed. */
static bool intrat_contains_rootsum(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_RootSum) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_contains_rootsum(e->data.function.args[i])) return true;
    return false;
}

static Expr* expand_simple_rootsum(Expr* poly, Expr* bvar, Expr* body) {
    /* Parameters are allowed only for a complex-coefficient linear or
     * quadratic factor: that is the factor LogToReal now declines (it is valid
     * only over a real field), and the root formula holds for any
     * coefficients, so 1/(x^2 + I a) keeps a closed form.  Real parametric
     * factors keep their RootSum, as before. */
    if (!poly_only_uses(poly, bvar)
        && (!intrat_contains_complex(poly) || get_degree_poly(poly, bvar) > 2))
        return NULL;

    Expr* cl_call = expr_new_function(expr_new_symbol(SYM_CoefficientList),
        (Expr*[]){ expr_copy(poly), expr_copy(bvar) }, 2);
    Expr* cl = evaluate(cl_call);
    expr_free(cl_call);
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List) {
        if (cl) expr_free(cl);
        return NULL;
    }

    size_t n = cl->data.function.arg_count;
    Expr* result = NULL;

    if (n == 2) {
        /* c0 + c1 t = 0  =>  t = -c0 / c1. */
        Expr* c0 = expr_copy(cl->data.function.args[0]);
        Expr* c1 = expr_copy(cl->data.function.args[1]);
        Expr* root_raw = internal_divide(
            (Expr*[]){ nlp_neg(c0), c1 }, 2);
        Expr* root = eval_and_free(root_raw);
        result = nlp_sub_body(body, bvar, root);
        expr_free(root);
    } else if (n == 3) {
        /* c0 + c1 t + c2 t^2 = 0
         * t = (-c1 ± Sqrt[c1² - 4 c0 c2]) / (2 c2). */
        Expr* c0 = expr_copy(cl->data.function.args[0]);
        Expr* c1 = expr_copy(cl->data.function.args[1]);
        Expr* c2 = expr_copy(cl->data.function.args[2]);

        Expr* disc_raw = internal_plus((Expr*[]){
            internal_times((Expr*[]){ expr_copy(c1), expr_copy(c1) }, 2),
            nlp_neg(internal_times((Expr*[]){
                expr_new_integer(4),
                expr_copy(c0), expr_copy(c2) }, 3))
        }, 2);
        Expr* disc = eval_and_free(disc_raw);
        Expr* sqd  = eval_and_free(nlp_sqrt(disc));

        Expr* two_c2 = eval_and_free(internal_times(
            (Expr*[]){ expr_new_integer(2), expr_copy(c2) }, 2));

        Expr* num_p = eval_and_free(internal_plus((Expr*[]){
            nlp_neg(expr_copy(c1)), expr_copy(sqd) }, 2));
        Expr* num_m = eval_and_free(internal_plus((Expr*[]){
            nlp_neg(expr_copy(c1)), nlp_neg(expr_copy(sqd)) }, 2));
        expr_free(sqd); expr_free(c0); expr_free(c1); expr_free(c2);

        Expr* root_p = eval_and_free(internal_divide(
            (Expr*[]){ num_p, expr_copy(two_c2) }, 2));
        Expr* root_m = eval_and_free(internal_divide(
            (Expr*[]){ num_m, two_c2 }, 2));

        Expr* t1 = nlp_sub_body(body, bvar, root_p);
        Expr* t2 = nlp_sub_body(body, bvar, root_m);
        expr_free(root_p); expr_free(root_m);

        result = eval_and_free(internal_plus((Expr*[]){ t1, t2 }, 2));
    } else if (n == 5) {
        /* Biquadratic: only c0, c2, c4 are non-zero.
         *   c0 + c2 t² + c4 t⁴ = 0
         *   u = (-c2 ± Sqrt[c2² - 4 c0 c4]) / (2 c4),  t = ±Sqrt[u]. */
        Expr* c1 = cl->data.function.args[1];
        Expr* c3 = cl->data.function.args[3];
        bool zero_c1 = (c1->type == EXPR_INTEGER && c1->data.integer == 0);
        bool zero_c3 = (c3->type == EXPR_INTEGER && c3->data.integer == 0);
        if (zero_c1 && zero_c3) {
            Expr* c0 = expr_copy(cl->data.function.args[0]);
            Expr* c2 = expr_copy(cl->data.function.args[2]);
            Expr* c4 = expr_copy(cl->data.function.args[4]);

            Expr* disc_raw = internal_plus((Expr*[]){
                internal_times((Expr*[]){ expr_copy(c2), expr_copy(c2) }, 2),
                nlp_neg(internal_times((Expr*[]){
                    expr_new_integer(4),
                    expr_copy(c0), expr_copy(c4) }, 3))
            }, 2);
            Expr* disc = eval_and_free(disc_raw);
            Expr* sqd  = eval_and_free(nlp_sqrt(disc));

            Expr* two_c4 = eval_and_free(internal_times(
                (Expr*[]){ expr_new_integer(2), expr_copy(c4) }, 2));

            Expr* up_num = eval_and_free(internal_plus((Expr*[]){
                nlp_neg(expr_copy(c2)), expr_copy(sqd) }, 2));
            Expr* um_num = eval_and_free(internal_plus((Expr*[]){
                nlp_neg(expr_copy(c2)), nlp_neg(expr_copy(sqd)) }, 2));
            expr_free(sqd); expr_free(c0); expr_free(c2); expr_free(c4);

            Expr* u_p = eval_and_free(internal_divide(
                (Expr*[]){ up_num, expr_copy(two_c4) }, 2));
            Expr* u_m = eval_and_free(internal_divide(
                (Expr*[]){ um_num, two_c4 }, 2));

            Expr* sp = eval_and_free(nlp_sqrt(u_p));
            Expr* sm = eval_and_free(nlp_sqrt(u_m));

            Expr* t1 = nlp_sub_body(body, bvar, sp);
            /* nlp_sub_body borrows its value: the negated roots are freed here */
            Expr* nsp = eval_and_free(nlp_neg(expr_copy(sp)));
            Expr* nsm = eval_and_free(nlp_neg(expr_copy(sm)));
            Expr* t2 = nlp_sub_body(body, bvar, nsp);
            Expr* t3 = nlp_sub_body(body, bvar, sm);
            Expr* t4 = nlp_sub_body(body, bvar, nsm);
            expr_free(sp); expr_free(sm); expr_free(nsp); expr_free(nsm);

            result = eval_and_free(internal_plus(
                (Expr*[]){ t1, t2, t3, t4 }, 4));
        }
    }

    expr_free(cl);
    return result;
}

/* ------------------------------------------------------------------ */
/* NaiveLogPart — RootSum fallback for the log part.                   */
/* ------------------------------------------------------------------ */

/* Direct port of IntegrateRational.m:1116-1124.
 *
 *   NaiveLogPart[f, x] = Σ_{α : d(α) = 0}  a(α) Log(x - α) / d'(α)
 *
 * for f = a/d.  Returned in held-symbolic RootSum form
 *
 *   RootSum[ Function[t, d(t)],
 *            Function[t, a(t) Log(x - t) / d'(t)] ]
 *
 * with `t` a fresh bound variable distinct from the outer integration
 * variable x.  The result is the universal fallback when LogToReal
 * cannot close the log part to a real elementary expression — every
 * proper rational integrand admits this representation, and the
 * D[RootSum] rule wired in src/deriv.c differentiates it back to f.
 *
 * This routine deliberately does NOT call ToRadicals on the
 * biquadratic case (cf. IntegrateRational.m:1121).  Algebraic-extension
 * radical closure is the responsibility of Solve / ToRadicals which we
 * have not yet implemented; until then RootSum is the honest answer. */
static Expr* intrat_naive_log_part(Expr* f, Expr* x) {
    intrat_trace("NaiveLogPart", "IN", f);

    /* a/d split.  We canonicalise f via Together so a and d are
     * polynomial in x even when f arrived as a sum of partial-fraction
     * pieces or with rational coefficients. */
    Expr* tg = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(tg);
    expr_free(tg);
    if (!combined) return NULL;

    Expr* a = eval_and_free(intrat_numerator(combined));
    Expr* d = eval_and_free(intrat_denominator(combined));
    expr_free(combined);

    /* d'(x), expanded so polynomial divisions downstream see it in
     * normal form. */
    Expr* dprime = intrat_d(d, x);
    Expr* dd = expr_expand(dprime); expr_free(dprime);

    /* Fresh bound variable.  The Private subcontext keeps it well clear
     * of any user symbol the caller might have set globally — using a
     * bare `t` would be unsafe because ReplaceAll evaluates the rule's
     * RHS, so a user-bound t = 5 would corrupt the substitution. */
    Expr* t = expr_new_symbol("Integrate`Private`t");

    /* Substitute t for x in a, d, dd via ReplaceAll. */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy(x), expr_copy(t) }, 2);

    Expr* a_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(a), expr_copy(rule) }, 2));
    Expr* d_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(d), expr_copy(rule) }, 2));
    Expr* dd_t = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(dd), expr_copy(rule) }, 2));
    expr_free(rule); expr_free(a); expr_free(d); expr_free(dd);

    /* Phase D2 — palindromic-quartic real-form expansion.  Builds the
     * conjugate-pair contributions directly, sidestepping the LogToReal
     * / LogToAtan chain that overwhelms on nested-radical S
     * polynomials.  Try this before constructing the held RootSum body
     * so we never produce a Function-headed RootSum for inputs whose
     * denominator is a palindromic quartic. */
    {
        Expr* palindromic_real = expand_palindromic_quartic_real(
            a_t, d_t, dd_t, t, x);
        if (palindromic_real) {
            expr_free(t); expr_free(a_t); expr_free(d_t); expr_free(dd_t);
            intrat_trace("NaiveLogPart", "OUT (palindromic-real)",
                         palindromic_real);
            return palindromic_real;
        }
    }

    /* Build the body: a_t * Log[x - t] / dd_t. */
    Expr* neg_t = internal_times((Expr*[]){
        expr_new_integer(-1),
        expr_copy(t)
    }, 2);
    Expr* x_minus_t = internal_plus(
        (Expr*[]){ expr_copy(x), neg_t }, 2);
    Expr* x_minus_t_e = eval_and_free(x_minus_t);

    Expr* log_term = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){ x_minus_t_e }, 1);

    Expr* numer = internal_times(
        (Expr*[]){ a_t, log_term }, 2);
    Expr* numer_e = eval_and_free(numer);

    Expr* body = internal_divide(
        (Expr*[]){ numer_e, dd_t }, 2);
    Expr* body_e = eval_and_free(body);

    /* Phase 8d-bonus finishing pass: if the polynomial d(t) has a
     * closed-form radical-formula root set (linear, quadratic, or
     * biquadratic), expand the RootSum to an explicit Plus over the
     * 1/2/4 roots and evaluate the body at each.  Otherwise fall
     * back to the held RootSum form.  Once Solve / ToRadicals lands
     * the simple-shape detector is replaced by Solve. */
    Expr* expanded = expand_simple_rootsum(d_t, t, body_e);
    if (expanded) {
        expr_free(t); expr_free(d_t); expr_free(body_e);
        intrat_trace("NaiveLogPart", "OUT (expanded)", expanded);
        return expanded;
    }

    /* Wrap in RootSum[Function[t, d_t], Function[t, body_e]]. */
    Expr* result = root_make_rootsum(t, d_t, body_e);
    intrat_trace("NaiveLogPart", "OUT", result);
    return result;
}

/* Divide p by the polynomial GCD of its NONZERO coefficients viewed as a
 * polynomial in cvar (the content of p in cvar).  internal PolynomialGCD
 * bails to a held call when handed a zero argument (PolynomialGCD[x^2, 0,
 * 4 x^2] returns 1, not x^2), so zero coefficients are filtered first.
 * Returns Cancel[p / content] when the content is a non-trivial polynomial,
 * else a copy of p.  Borrows p; caller owns the result.  Used by the
 * transcendental log-part path to strip the x-only content that the
 * denominator-clearing scaling introduces into the Rothstein-Trager
 * resultant, exposing the (constant-residue) polynomial in the residue
 * variable. */
static Expr* intrat_strip_content_var(Expr* p, Expr* cvar) {
    Expr* cl = eval_and_free(internal_coefficientlist(
        (Expr*[]){expr_copy(p), expr_copy(cvar)}, 2));
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List
        || cl->data.function.arg_count == 0) {
        if (cl) expr_free(cl);
        return expr_copy(p);
    }
    size_t n = cl->data.function.arg_count;
    Expr** nz = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nnz = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* c = cl->data.function.args[i];
        if (!is_zero_poly(c)) nz[nnz++] = expr_copy(c);
    }
    expr_free(cl);
    if (nnz == 0) { free(nz); return expr_copy(p); }
    Expr* content = (nnz == 1) ? nz[0] : internal_polynomialgcd(nz, nnz);
    free(nz);
    if (head_is(content, SYM_PolynomialGCD)
        || (content->type == EXPR_INTEGER
            && (content->data.integer == 1 || content->data.integer == -1))) {
        expr_free(content);
        return expr_copy(p);
    }
    Expr* div = internal_divide((Expr*[]){expr_copy(p), content}, 2);
    return eval_and_free(internal_cancel((Expr*[]){div}, 1));
}

/* True iff polynomial p is free of every variable named in `xgate` (a single
 * symbol, or a List of symbols) — i.e. the roots of p are constants of the
 * derivation.  Used by the residue criterion's κ_D split. */
static bool intrat_freeq_all(Expr* p, Expr* xgate) {
    if (xgate->type == EXPR_FUNCTION
        && xgate->data.function.head->type == EXPR_SYMBOL
        && xgate->data.function.head->data.symbol.name == SYM_List) {
        for (size_t gi = 0; gi < xgate->data.function.arg_count; gi++)
            if (!intrat_freeq_test(p, xgate->data.function.args[gi]))
                return false;
        return true;
    }
    return intrat_freeq_test(p, xgate);
}

/* Recursively walk a Factor[] result, routing each irreducible factor (base of
 * a Power, or a bare polynomial in `t`) into the constant-residue product *Qs
 * (free of the gate vars) or the non-constant product *Qn.  Numeric / t-free
 * factors are dropped.  *Qs, *Qn are running products (owned), each starting at
 * integer 1; internal_times adopts them, so they are consumed and replaced. */
static void intrat_route_factor(Expr* f, Expr* t, Expr* xgate,
                                Expr** Qs, Expr** Qn) {
    if (f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL) {
        const char* h = f->data.function.head->data.symbol.name;
        if (h == SYM_Times) {
            for (size_t k = 0; k < f->data.function.arg_count; k++)
                intrat_route_factor(f->data.function.args[k], t, xgate, Qs, Qn);
            return;
        }
        if (h == SYM_Power) {   /* squarefree bucket => exponent 1: route base */
            intrat_route_factor(f->data.function.args[0], t, xgate, Qs, Qn);
            return;
        }
    }
    if (get_degree_poly(f, t) <= 0) return;   /* numeric / t-free: drop */
    Expr** dst = intrat_freeq_all(f, xgate) ? Qs : Qn;
    *dst = expand_and_free(internal_times((Expr*[]){*dst, expr_copy(f)}, 2));
}

/* Split a squarefree bucket polynomial Qi (in `t`) into its constant-residue
 * sub-product (returned) and, multiplied into *Qn_acc, its non-constant-residue
 * sub-product.  Returns the constant product (owned; integer 1 if none). */
static Expr* intrat_bucket_const_part(Expr* Qi, Expr* t, Expr* xgate,
                                      Expr** Qn_acc) {
    Expr* factored = eval_and_free(internal_factor((Expr*[]){expr_copy(Qi)}, 1));
    Expr* Qs = expr_new_integer(1);
    if (factored) { intrat_route_factor(factored, t, xgate, &Qs, Qn_acc); expr_free(factored); }
    return Qs;
}

/* IntRationalLogPart core. Returns either the {Q, S} pair list
 * (root_sum=false) or the symbolic RootSum sum (root_sum=true).
 * On algorithmic failure returns NULL.
 *
 * `a`, `d` are the (owned) numerator / denominator in the main variable `x`
 * with squarefree `d`; `t` is the residue variable.  Two optional hooks
 * generalise the pure-rational computation to a transcendental monomial
 * extension (the recursive-Risch Lazard-Rioboo-Trager log part):
 *   - `dprime_override != NULL` injects the monomial derivation `D(d)` in
 *     place of the ordinary `d/dx` (here `x` is the monomial variable);
 *   - `xgate != NULL` requires the residues (roots of the resultant) to be
 *     constants of the derivation (free of `xgate`).  Bronstein's residue
 *     criterion (Thm 5.6.1) splits the resultant r = r_s · r_n by whether each
 *     factor's residues are constant: the constant part r_s yields elementary
 *     logs, the non-constant part r_n does not.  When `remainder_out != NULL`
 *     and the split is genuinely mixed, the elementary r_s logs are returned AND
 *     *remainder_out is set to the (owned) r_n rational-function remainder to be
 *     reported unintegrated (partial log part).  When `remainder_out == NULL`
 *     (caller refuses partials) a mixed resultant declines (NULL), as does the
 *     all-non-constant case; the all-constant case is unchanged.
 * With both NULL this is the original pure-rational log part.
 *
 * Direct port of IntegrateRational.m:751-801. */
/* Root-mode log part over Q(Root) (defined with the Q(alpha) helpers). */
typedef struct IrRootMode IrRootMode;
static Expr* intrat_log_part_core_root(Expr* a, Expr* d, Expr* x, Expr* t,
                                       bool root_sum, bool* handled);

static Expr* intrat_log_part_core(Expr* a, Expr* d, Expr* x, Expr* t,
                                  bool root_sum, Expr* dprime_override,
                                  Expr* xgate, Expr** remainder_out, bool decide) {
    if (remainder_out) *remainder_out = NULL;
    /* Coefficients in Q(alpha) for one Root object alpha: the whole log
     * part over that number field (see intrat_log_part_core_root). */
    if (!dprime_override && !xgate) {
        bool handled = false;
        Expr* rr = intrat_log_part_core_root(a, d, x, t, root_sum, &handled);
        if (handled) return rr;
    }
    /* prs = SubresultantPolynomialRemainders[d, a - t*D(d), x] */
    Expr* dprime_e;
    if (dprime_override) {
        dprime_e = expr_expand(dprime_override);
    } else {
        Expr* dprime = intrat_d(d, x);
        dprime_e = expr_expand(dprime); expr_free(dprime);
    }
    Expr* tdprime = internal_times(
        (Expr*[]){expr_copy(t), dprime_e}, 2);
    Expr* atdp_raw = internal_subtract(
        (Expr*[]){expr_copy(a), tdprime}, 2);
    Expr* atdp;
    if (dprime_override) {
        /* the monomial derivation D(d) can carry x-denominators (e.g. the
         * u'/u of a logarithmic kernel); clear them so the resultant / PRS
         * operate on genuine polynomials.  The x-only denominator scales the
         * resultant by a unit that the x-content strip below removes. */
        Expr* tg = eval_and_free(internal_together((Expr*[]){atdp_raw}, 1));
        Expr* nn = intrat_numerator(tg);
        atdp = expr_expand(nn); expr_free(nn); expr_free(tg);
    } else {
        atdp = expr_expand(atdp_raw); expr_free(atdp_raw);
    }

    /* The subresultant PRS of d and a - t D(d) in Z[t][x], natively (FLINT)
     * for rational coefficients; the generic Expr pseudo-remainder chain
     * otherwise (parameters, radicals, the transcendental callers).  The
     * generic chain has exponential coefficient growth: 19 s for
     * Czichowski's degree-14 example, against 10 ms natively. */
    Expr* prs = intrat_subres_chain(d, atdp, x, t);
    if (!prs) prs = intrat_subres_chain_mv(d, atdp, x, t);
    if (!prs) {
        Expr* prs_args[3] = { expr_copy(d), expr_copy(atdp), expr_copy(x) };
        Expr* prs_call = expr_new_function(
            expr_new_symbol(SYM_SubresultantPolynomialRemainders), prs_args, 3);
        prs = evaluate(prs_call);
        expr_free(prs_call);
    }
    if (!prs) {
        expr_free(a); expr_free(d); expr_free(atdp); return NULL;
    }

    /* resultant = primitive[Resultant[d, a - t*D[d,x], x], t]
     *
     * NOTE: this Resultant[] call is NOT redundant with the PRS chain
     * computed above, despite both being subresultant-PRS-flavoured
     * computations.  Mathilda's SubresultantPolynomialRemainders is a
     * pseudo-remainder chain (poly.c:1637-1642), not the Lazard-scaled
     * subresultant chain.  The chain's primitive-part-in-t and
     * degree-in-x are content-invariant and so are correct substrates
     * for the S_i lookups below (i >= 1); but the chain's degree-0-in-x
     * entry is NOT in general the resultant, even after primitive in
     * t.  Concrete witness: x/(1 + x^8) has resultant (1 + 4096 t^4)^2
     * but PRS-bottom primitive in t is t (1 + 4096 t^4), structurally
     * different. */
    /* The residue map is a/D(d) = P(x)/W(x) where atdp = P - t*W is linear in t.
     * Snapshot it before the Resultant call adopts `atdp` — the κ_D split's
     * remainder reconstruction resultant needs it (partial log part). */
    Expr* atdp_keep = xgate ? expr_copy(atdp) : NULL;
    Expr* res_args[3] = { expr_copy(d), atdp, expr_copy(x) };
    Expr* res_call = expr_new_function(
        expr_new_symbol(SYM_Resultant), res_args, 3);
    Expr* resultant_raw = evaluate(res_call);
    expr_free(res_call);
    if (!resultant_raw) {
        expr_free(a); expr_free(d); expr_free(prs);
        if (atdp_keep) expr_free(atdp_keep);
        return NULL;
    }
    Expr* resultant = intrat_primitive(resultant_raw, t);
    expr_free(resultant_raw);

    /* κ_D split state (partial log part).  In `partial_mode` the resultant has a
     * genuine mix of constant- and non-constant-residue factors: emit only the
     * constant-residue (r_s) logs and report the r_n part via *remainder_out. */
    bool partial_mode = false;
    Expr* remainder_local = NULL;

    if (xgate) {
        /* Strip the lower-field content the denominator-clearing scaling leaves,
         * then classify the residues (roots of the resultant) by whether they
         * are constants of the derivation (free of every gate variable).  A
         * single symbol gates a single-kernel extension; a `List` gates a tower
         * proper part.  Bronstein's residue criterion (Thm 5.6.1) keeps the
         * elementary constant-residue logs even when non-constant residues are
         * present; the fully-constant case is unchanged. */
        Expr* stripped = intrat_strip_content_var(resultant, t);
        expr_free(resultant);
        resultant = stripped;
        bool free_all = intrat_freeq_all(resultant, xgate);
        if (!free_all && decide) {
            /* Decision mode: ANY non-constant residue proves this simple part has no
             * elementary integral (Bronstein residue criterion, Thm 5.6.1(ii)).  We
             * only need the VERDICT, not the elementary constant-residue logs or a
             * reconstructed remainder — return an inert marker the recursive-Risch
             * caller (rt_field_lrt_logpart) decodes into a NONELEMENTARY verdict. */
            expr_free(resultant); expr_free(a); expr_free(d); expr_free(prs);
            if (atdp_keep) expr_free(atdp_keep);
            return expr_new_symbol("Integrate`$NonConstantResidue");
        }
        if (!free_all) {
            /* Mixed or fully non-constant.  A caller that cannot represent a
             * partial result (remainder_out == NULL — e.g. the exponential
             * coupled path) declines exactly as before.  Otherwise defer to the
             * per-factor partition after the squarefree split below. */
            if (!remainder_out) {
                expr_free(resultant); expr_free(a); expr_free(d); expr_free(prs);
                if (atdp_keep) expr_free(atdp_keep);
                return NULL;
            }
            partial_mode = true;
        }
    }

    /* Q = SquareFree[resultant]: list of {poly_i, mult_i}. */
    Expr* Q = intrat_squarefree_list(resultant);
    expr_free(resultant);
    if (!Q) {
        expr_free(a); expr_free(d); expr_free(prs);
        if (atdp_keep) expr_free(atdp_keep);
        return NULL;
    }

    int deg_d = get_degree_poly(d, x);
    size_t nQ = Q->data.function.arg_count;

    /* Partial κ_D partition: per bucket, the constant-residue sub-product
     * Qconst[i] (used in the output instead of the full bucket) and the global
     * non-constant product Q_n_poly (whose residues drive the remainder).  Only
     * built in partial_mode; NULL entries mean "no constant residues in bucket". */
    Expr** Qconst = NULL;
    if (partial_mode) {
        Qconst = (Expr**)malloc(sizeof(Expr*) * nQ);
        Expr* Q_n_poly = expr_new_integer(1);
        bool any_const = false;
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
            Expr* Qi = list_get(Q->data.function.args[i_idx], 1);
            if (!Qi || get_degree_poly(Qi, t) <= 0) { if (Qi) expr_free(Qi); Qconst[i_idx] = NULL; continue; }
            Qconst[i_idx] = intrat_bucket_const_part(Qi, t, xgate, &Q_n_poly);
            expr_free(Qi);
            if (Qconst[i_idx] && get_degree_poly(Qconst[i_idx], t) > 0) any_const = true;
        }
        if (!any_const) {
            /* No elementary logs to emit — decline (also prevents an infinite
             * Integrate[remainder] regress).  Decision mode never reaches here: it
             * short-circuits at the non-constant-residue detection above. */
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst); expr_free(Q_n_poly); expr_free(Q);
            expr_free(a); expr_free(d); expr_free(prs);
            if (atdp_keep) expr_free(atdp_keep);
            return NULL;
        }
        /* Remainder = the d_n-part of a/d, where d_n = gcd(d, Res_t(Q_n, atdp))
         * collects the poles whose residue is a root of Q_n (non-constant).
         * atdp = P - t*W is linear in t, so Res_t(Q_n, atdp) = W^deg * Q_n(P/W)
         * (a poly in x); its gcd with d picks exactly those poles. */
        Expr* resn_call = expr_new_function(expr_new_symbol(SYM_Resultant),
            (Expr*[]){ Q_n_poly, atdp_keep, expr_copy(t) }, 3);   /* adopts both */
        atdp_keep = NULL;
        Expr* Res_n = eval_and_free(resn_call);
        Expr* d_n = Res_n ? intrat_polygcd_monic(d, Res_n, x) : NULL;
        if (Res_n) expr_free(Res_n);
        Expr* d_s = d_n ? intrat_exquo(d, d_n, x) : NULL;
        Expr* rr = NULL; Expr* ss = NULL;
        bool ok = d_n && d_s
            && intrat_extended_euclidean(d_n, d_s, a, x, &rr, &ss);  /* rr d_n + ss d_s = a */
        if (rr) expr_free(rr);
        if (ok) {
            /* a/(d_n d_s) = rr/d_s + ss/d_n; the non-constant (d_n) part is ss/d_n. */
            remainder_local = eval_and_free(internal_together((Expr*[]){
                internal_divide((Expr*[]){ ss, expr_copy(d_n) }, 2) }, 1));
        } else if (ss) { expr_free(ss); }
        if (d_n) expr_free(d_n);
        if (d_s) expr_free(d_s);
        if (!remainder_local) {
            /* Could not reconstruct the remainder — decline rather than emit a
             * partial result whose residual identity we cannot certify. */
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst); expr_free(Q);
            expr_free(a); expr_free(d); expr_free(prs);
            return NULL;
        }
    }
    if (atdp_keep) { expr_free(atdp_keep); atdp_keep = NULL; }

    /* Compute S[i] for each i.  Where Q[[i,1]] has non-trivial t-
     * degree, we derive a polynomial in (t, x). */
    Expr** S = (Expr**)malloc(sizeof(Expr*) * nQ);
    for (size_t i = 0; i < nQ; i++) S[i] = NULL;

    for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
        int i = (int)(i_idx + 1);
        Expr* Qi_pair = Q->data.function.args[i_idx];
        Expr* Qi = list_get(Qi_pair, 1);
        if (!Qi) continue;

        if (get_degree_poly(Qi, t) <= 0) {
            expr_free(Qi);
            continue;
        }

        if (i == deg_d) {
            /* Special case (Geddes-Czapor-Labahn correction
             * IntegrateRational.m:765-768): use d itself with lc[d, x]
             * as the leading-coefficient witness in the Diophantine. */
            S[i_idx] = intrat_primitive_part_mod(d, Qi, x, t);
        } else {
            Expr* s_raw = find_prs_at_degree(prs, x, i);
            if (!s_raw) { expr_free(Qi); continue; }
            Expr* s = intrat_primitive(s_raw, x);
            expr_free(s_raw);

            /* A = SquareFree[primitive[lc[s, x], t]]. */
            Expr* lc_s = intrat_lc(s, x);
            Expr* lc_s_prim = intrat_primitive(lc_s, t);
            expr_free(lc_s);
            Expr* A = intrat_squarefree_list(lc_s_prim);
            expr_free(lc_s_prim);
            if (!A) { expr_free(Qi); expr_free(s); continue; }

            /* Divide common factors between A_j and Qi out of s. */
            for (size_t j = 0; j < A->data.function.arg_count; j++) {
                Expr* A_j_pair = A->data.function.args[j];
                Expr* A_jpoly = list_get(A_j_pair, 1);
                if (!A_jpoly) continue;
                Expr* g_args[3] = {A_jpoly, expr_copy(Qi), expr_copy(t)};
                /* Use PolynomialExtendedGCD-derived monic gcd to
                 * avoid the content-scaling issue we hit in Phase 1. */
                Expr* g = intrat_polygcd_monic(g_args[0], g_args[1], g_args[2]);
                expr_free(A_jpoly); expr_free(g_args[1]); expr_free(g_args[2]);
                if (!g || (g->type == EXPR_INTEGER && g->data.integer == 1)) {
                    if (g) expr_free(g);
                    continue;
                }
                /* s := exquo[s, g^j, t] (j is 1-based index = j+1). */
                int jpow = (int)(j + 1);
                Expr* gj = expand_and_free(internal_power(
                    (Expr*[]){g, expr_new_integer((int64_t)jpow)}, 2));
                Expr* sn = intrat_exquo(s, gj, t);
                expr_free(gj);
                if (sn) { expr_free(s); s = sn; }
            }
            expr_free(A);

            /* s := PrimitivePart_mod(s, Qi) -- see helper definition. */
            S[i_idx] = intrat_primitive_part_mod(s, Qi, x, t);
            expr_free(s);
        }
        expr_free(Qi);
    }

    /* Every Q_i of positive degree has an S_i (Bronstein §2.5: the
     * subresultant of x-degree i exists for each nontrivial squarefree
     * factor of the resultant).  A failed lookup or normalisation above used
     * to `continue`, and the output loops below then skipped that Q_i: its
     * whole RootSum vanished and the antiderivative was silently wrong.
     * Decline instead. */
    {
        bool incomplete = false;
        for (size_t i_idx = 0; i_idx < nQ && !incomplete; i_idx++) {
            if (S[i_idx]) continue;
            Expr* Qi = partial_mode ? (Qconst[i_idx] ? expr_copy(Qconst[i_idx]) : NULL)
                                    : list_get(Q->data.function.args[i_idx], 1);
            if (Qi && get_degree_poly(Qi, t) > 0) incomplete = true;
            if (Qi) expr_free(Qi);
        }
        if (incomplete) {
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
            free(S);
            if (partial_mode) {
                for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
                free(Qconst);
                if (remainder_local) expr_free(remainder_local);
            }
            expr_free(Q); expr_free(prs);
            expr_free(a); expr_free(d);
            return NULL;
        }
    }
    expr_free(prs);
    expr_free(a); expr_free(d);

    /* Build output. */
    if (root_sum) {
        /* Sum of RootSum heads. */
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * nQ);
        size_t nterms = 0;
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
            if (!S[i_idx]) continue;
            /* In partial_mode restrict the RootSum to the constant-residue
             * sub-product of the bucket (S[i] mod Qi agrees on its sub-roots). */
            Expr* Qi;
            if (partial_mode) {
                if (!Qconst[i_idx] || get_degree_poly(Qconst[i_idx], t) <= 0) continue;
                Qi = expr_copy(Qconst[i_idx]);
            } else {
                Qi = list_get(Q->data.function.args[i_idx], 1);
                if (!Qi) continue;
                if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }
            }

            /* RootSum[Function[t, Qi], Function[t, t Log[S[i]]]].  Canonicalise
             * the held Function bodies up front (Function is HoldAll, so the
             * evaluator will not normalise them afterwards): sort the resultant
             * factor Qi and fold the t*Log[S] product so the RootSum prints in
             * normal form. */
            Qi = eval_and_free(Qi);
            Expr* func1 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){expr_copy(t), Qi}, 2);
            Expr* logS = expr_new_function(expr_new_symbol(SYM_Log),
                (Expr*[]){expr_copy(S[i_idx])}, 1);
            Expr* tlog = eval_and_free(internal_times(
                (Expr*[]){expr_copy(t), logS}, 2));
            Expr* func2 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){expr_copy(t), tlog}, 2);
            terms[nterms++] = expr_new_function(expr_new_symbol(SYM_RootSum),
                (Expr*[]){func1, func2}, 2);
        }
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
        free(S);
        expr_free(Q);
        if (partial_mode) {
            for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
            free(Qconst);
            *remainder_out = remainder_local;
        }
        Expr* result;
        if (nterms == 0) { free(terms); result = expr_new_integer(0); }
        else if (nterms == 1) { result = terms[0]; free(terms); }
        else { result = internal_plus(terms, nterms); free(terms); result = eval_and_free(result); }
        intrat_trace("IntRationalLogPart", "OUT", result);
        return result;
    }

    /* RootSum -> False: list of {Q_i, S_i} pairs. */
    Expr** pairs = (Expr**)malloc(sizeof(Expr*) * nQ);
    size_t npairs = 0;
    for (size_t i_idx = 0; i_idx < nQ; i_idx++) {
        if (!S[i_idx]) continue;
        Expr* Qi;
        if (partial_mode) {
            if (!Qconst[i_idx] || get_degree_poly(Qconst[i_idx], t) <= 0) continue;
            Qi = expr_copy(Qconst[i_idx]);
        } else {
            Qi = list_get(Q->data.function.args[i_idx], 1);
            if (!Qi) continue;
            if (get_degree_poly(Qi, t) <= 0) { expr_free(Qi); continue; }
        }
        pairs[npairs++] = expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){Qi, expr_copy(S[i_idx])}, 2);
    }
    for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (S[i_idx]) expr_free(S[i_idx]);
    free(S);
    expr_free(Q);
    if (partial_mode) {
        for (size_t i_idx = 0; i_idx < nQ; i_idx++) if (Qconst[i_idx]) expr_free(Qconst[i_idx]);
        free(Qconst);
        *remainder_out = remainder_local;
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), pairs, npairs);
    free(pairs);
    intrat_trace("IntRationalLogPart", "OUT", result);
    return result;
}

/* Pure-rational log part: Together-canonicalise f, split into a/d, and run
 * the core with the ordinary derivation (no transcendental hooks). */
static Expr* intrat_int_rational_log_part(Expr* f, Expr* x, Expr* t,
                                          bool root_sum) {
    intrat_trace("IntRationalLogPart", "IN", f);
    /* Together-canonicalise f so the numerator and denominator are clean
     * polynomials with no internal fractions in parametric coefficients. */
    Expr* fc = eval_and_free(internal_together((Expr*[]){expr_copy(f)}, 1));
    Expr* a = eval_and_free(intrat_numerator(fc));
    Expr* d_raw = intrat_denominator(fc);
    Expr* d = expr_expand(d_raw); expr_free(d_raw);
    expr_free(fc);
    return intrat_log_part_core(a, d, x, t, root_sum, NULL, NULL, NULL, false);
}

/* Defined below (Phase 4 consumer); forward-declared for the transcendental
 * log part which hands its {Q_i, S_i} pairs straight to Rioboo's LogToReal. */
static Expr* intrat_log_to_real_pairs(Expr* pair_list, Expr* x, Expr* t);

/* True when e contains a non-real number (Complex[_, _]). */
static bool intrat_has_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex) return true;
    if (intrat_has_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_has_complex(e->data.function.args[i])) return true;
    return false;
}

/* Is the residue-factor set {Q_i} of a pair list {{Q_i, S_i}, ...} closed
 * under complex conjugation?  Every Q_i with a non-real coefficient must have
 * a Q_j with Q_j / conj(Q_i) free of the residue variable (equal up to a
 * constant); conj maps Complex[a, b] to Complex[a, -b] (the parameters and
 * the residue variable are real symbols here).  Real factors pass. */
static bool intrat_pairs_conjugate_closed(Expr* pairs) {
    if (!pairs || pairs->type != EXPR_FUNCTION) return false;
    size_t n = pairs->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        Expr* pi = pairs->data.function.args[i];
        if (pi->type != EXPR_FUNCTION || pi->data.function.arg_count != 2) return false;
        Expr* Qi = pi->data.function.args[0];
        if (!intrat_has_complex(Qi)) continue;
        Expr* rule = parse_expression("Complex[intratCa_, intratCb_] :> Complex[intratCa, -intratCb]");
        Expr* cq = eval_and_free(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){expr_copy(Qi), rule}, 2));
        bool found = false;
        for (size_t j = 0; j < n && !found; j++) {
            Expr* pj = pairs->data.function.args[j];
            if (pj->type != EXPR_FUNCTION || pj->data.function.arg_count != 2) continue;
            Expr* ratio = eval_and_free(internal_together((Expr*[]){
                internal_divide((Expr*[]){expr_copy(pj->data.function.args[0]), expr_copy(cq)}, 2)}, 1));
            /* constant: no variable of Q_i (the residue variable) survives */
            Expr* vars = eval_and_free(expr_new_function(expr_new_symbol("Variables"),
                (Expr*[]){expr_copy(Qi)}, 1));
            bool free = true;
            if (vars && vars->type == EXPR_FUNCTION)
                for (size_t k = 0; k < vars->data.function.arg_count && free; k++) {
                    Expr* fq = eval_and_free(expr_new_function(expr_new_symbol("FreeQ"),
                        (Expr*[]){expr_copy(ratio), expr_copy(vars->data.function.args[k])}, 2));
                    free = fq && fq->type == EXPR_SYMBOL && fq->data.symbol.name == SYM_True;
                    if (fq) expr_free(fq);
                }
            if (vars) expr_free(vars);
            if (ratio) expr_free(ratio);
            found = free;
        }
        expr_free(cq);
        if (!found) return false;
    }
    return true;
}

/* Sum_{Q(r) = 0} r Log[S(r)] written out when Q has degree 1 or 2 in z;
 * NULL otherwise (or on a malformed coefficient list).  All borrowed. */
static Expr* intrat_expand_small_root_sum(Expr* Q, Expr* S, Expr* z) {
    Expr* cl = eval_and_free(internal_coefficientlist(
        (Expr*[]){expr_expand(Q), expr_copy(z)}, 2));
    if (!cl || cl->type != EXPR_FUNCTION
        || cl->data.function.head->type != EXPR_SYMBOL
        || cl->data.function.head->data.symbol.name != SYM_List
        || cl->data.function.arg_count < 2 || cl->data.function.arg_count > 3) {
        if (cl) expr_free(cl);
        return NULL;
    }
    Expr** c = cl->data.function.args;
    Expr* roots[2]; size_t nr = 0;
    if (cl->data.function.arg_count == 2) {
        /* c0 + c1 z: r = -c0/c1 */
        roots[nr++] = eval_and_free(internal_times((Expr*[]){expr_new_integer(-1),
            internal_divide((Expr*[]){expr_copy(c[0]), expr_copy(c[1])}, 2)}, 2));
    } else {
        /* c0 + c1 z + c2 z^2: r = (-c1 +- Sqrt[c1^2 - 4 c0 c2])/(2 c2) */
        Expr* disc = internal_plus((Expr*[]){
            internal_power((Expr*[]){expr_copy(c[1]), expr_new_integer(2)}, 2),
            internal_times((Expr*[]){expr_new_integer(-4), expr_copy(c[0]), expr_copy(c[2])}, 3)}, 2);
        Expr* sq = eval_and_free(expr_new_function(expr_new_symbol("Sqrt"), (Expr*[]){disc}, 1));
        for (int sg = 1; sg >= -1; sg -= 2)
            roots[nr++] = eval_and_free(internal_divide((Expr*[]){
                internal_plus((Expr*[]){
                    internal_times((Expr*[]){expr_new_integer(-1), expr_copy(c[1])}, 2),
                    internal_times((Expr*[]){expr_new_integer(sg), expr_copy(sq)}, 2)}, 2),
                internal_times((Expr*[]){expr_new_integer(2), expr_copy(c[2])}, 2)}, 2));
        expr_free(sq);
    }
    expr_free(cl);
    Expr* terms[2];
    for (size_t i = 0; i < nr; i++) {
        Expr* rule = expr_new_function(expr_new_symbol("Rule"), (Expr*[]){expr_copy(z), expr_copy(roots[i])}, 2);
        Expr* Sr = eval_and_free(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){expr_copy(S), rule}, 2));
        terms[i] = internal_times((Expr*[]){roots[i],
            expr_new_function(expr_new_symbol(SYM_Log), (Expr*[]){Sr}, 1)}, 2);
    }
    if (nr == 1) return eval_and_free(terms[0]);
    return eval_and_free(internal_plus(terms, nr));
}

/* {{Q_1, S_1}, ...} (the RootSum -> False pair list of intrat_log_part_core,
 * Q_i in the residue variable z, S_i in z and the monomial) as
 *     Sum_i RootSum[Function[z, Q_i], Function[z, z Log[S_i]]]
 * -- the log part with its residues left algebraic.  NULL on a malformed
 * or empty list.  `pairs` and `z` are borrowed. */
static Expr* intrat_root_sum_pairs(Expr* pairs, Expr* z) {
    if (!pairs || pairs->type != EXPR_FUNCTION) return NULL;
    size_t n = pairs->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* pr = pairs->data.function.args[i];
        if (pr->type != EXPR_FUNCTION || pr->data.function.arg_count != 2) {
            for (size_t k = 0; k < nt; k++) expr_free(terms[k]);
            free(terms);
            return NULL;
        }
        Expr* Qi = eval_and_free(expr_copy(pr->data.function.args[0]));
        /* Degree 1 or 2 in z: the roots are explicit, so write the sum out,
         * Sum_{Q(r) = 0} r Log[S(r)] (the quadratic formula; Sqrt of a
         * complex discriminant is fine -- the log part is an identity in the
         * roots, whatever they are).  Degree >= 3 stays a RootSum. */
        Expr* expl = intrat_expand_small_root_sum(Qi, pr->data.function.args[1], z);
        if (expl) { expr_free(Qi); terms[nt++] = expl; continue; }
        Expr* func1 = expr_new_function(expr_new_symbol(SYM_Function),
            (Expr*[]){expr_copy(z), Qi}, 2);
        Expr* logS = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){expr_copy(pr->data.function.args[1])}, 1);
        Expr* zlog = eval_and_free(internal_times((Expr*[]){expr_copy(z), logS}, 2));
        Expr* func2 = expr_new_function(expr_new_symbol(SYM_Function),
            (Expr*[]){expr_copy(z), zlog}, 2);
        terms[nt++] = expr_new_function(expr_new_symbol(SYM_RootSum),
            (Expr*[]){func1, func2}, 2);
    }
    /* No pair: the core dropped every factor (S_i undetermined, e.g. a
     * residue factor over Q(i)) -- a failure, never a zero log part. */
    if (nt == 0) { free(terms); return NULL; }
    Expr* result;
    if (nt == 1) { result = terms[0]; free(terms); }
    else { result = internal_plus(terms, nt); free(terms); result = eval_and_free(result); }
    return result;
}

/* Transcendental (recursive-Risch) Lazard-Rioboo-Trager log part.
 *
 * Integrates the PROPER, SQUAREFREE-denominator fraction  a(tau)/d(tau)  of a
 * single transcendental monomial extension tau (= Log[u] or E^u over Q(x)),
 * whose log part carries ALGEBRAIC (non-rational) residues — the class the
 * single-constant-per-factor SolveAlways ansatz cannot express.  `Dd` is the
 * monomial derivation D(d) = d/dx(d) + (D tau) d/dtau(d), computed by the
 * caller.  The residues are the roots of the Rothstein-Trager resultant
 * Res_tau(a - z Dd, d) (z = the residue variable); each squarefree factor
 * Q_i(z) with its log argument S_i(z, tau) is converted to real Log + ArcTan
 * form by Rioboo's LogToReal.  Returns the antiderivative as a function of
 * `tau` (the caller substitutes tau -> kernel), or NULL when the integral is
 * not elementary in this form / a factor exceeds LogToReal's bounded scope.
 * `a`, `d`, `Dd`, `tau`, `z`, `xreal` are all borrowed.
 *
 * When `remainder_out != NULL`, the residue criterion may return a PARTIAL log
 * part (Bronstein Thm 5.6.1): the elementary constant-residue logs, with the
 * non-constant-residue rational remainder written to *remainder_out (owned, a
 * function of `tau`) for the caller to report unintegrated.  *remainder_out is
 * NULL when the whole log part is elementary. */
static Expr* intrat_transcendental_log_part(Expr* a, Expr* d, Expr* tau,
                                            Expr* z, Expr* Dd, Expr* xreal,
                                            Expr** remainder_out, bool decide) {
    if (remainder_out) *remainder_out = NULL;
    Expr* rem = NULL;
    Expr* pairs = intrat_log_part_core(expr_copy(a), expr_copy(d), tau, z,
                                       /*root_sum=*/false, Dd, xreal, &rem, decide);
    if (!pairs) { if (rem) expr_free(rem); return NULL; }
    /* Decision mode short-circuit: a non-constant residue yields the inert marker
     * Integrate`$NonConstantResidue (no pair list) — propagate it up verbatim for
     * the recursive-Risch caller to decode into a NONELEMENTARY verdict. */
    if (pairs->type == EXPR_SYMBOL
        && strcmp(pairs->data.symbol.name, "Integrate`$NonConstantResidue") == 0) {
        if (rem) expr_free(rem);
        return pairs;
    }
    Expr* real = intrat_log_to_real_pairs(pairs, tau, z);
    /* Rioboo's LogToReal needs the residues closed under complex conjugation:
     * every residue factor Q_i with a non-real coefficient must meet its
     * conjugate among the factors.  The trig callers (tau = E^(I x)) satisfy
     * it, and need the real reconstruction.  A factor without its conjugate
     * (1 + 4 I z^2, of 1/(u (I + Log[u]^2))) does not, and LogToReal returned
     * a wrong 0 there; v0.6593 declined every complex pair instead, which sent
     * Integrate[Log[61/36 - 5/3 Cos[2 x]] Sin[x] Sin[9 x], {x, 0, Pi}] from 5 s
     * past 100 s (fixed v0.6599).  Not closed: the RootSum below. */
    if (real && !intrat_pairs_conjugate_closed(pairs)) {
        expr_free(real);
        real = NULL;
    }
    if (!real) {
        /* LogToReal declined (a residue factor of degree >= 3), or its
         * precondition failed (above): the residues are algebraic, and Bronstein §5.6 allows
         * them in the algebraic closure, so the log part is the RootSum
         *     Sum_i RootSum[Q_i(z), z Log[S_i(z, tau)]]
         * (review F7: 1/(u (1 + Log[u] + Log[u]^3)) was a 27-68 s decline;
         * the rational integrator already answers such cases this way). */
        real = intrat_root_sum_pairs(pairs, z);
    }
    expr_free(pairs);
    if (!real) { if (rem) expr_free(rem); return NULL; }
    if (remainder_out) *remainder_out = rem;
    else if (rem) expr_free(rem);   /* caller refuses partial: shouldn't happen */
    return real;
}

/* ====================================================================
 * Phase 3 — LogToAtan (Rioboo's recursive conversion).
 * ====================================================================
 *
 * Given A, B in K[x] with B != 0, return a sum of ArcTan of polynomials
 * in K[x] such that
 *    D[f, x] = D[I Log[(A + I B)/(A - I B)], x].
 *
 * Direct port of IntegrateRational.m:1529-1561.  Used by Phase 4's
 * LogToReal: each pair (A, B) of complex Re/Im polynomial parts feeds
 * LogToAtan to produce the corresponding real ArcTan summand.
 */

static Expr* intrat_log_to_atan(Expr* a, Expr* b, Expr* x) {
    intrat_trace("LogToAtan", "IN", a);

    /* If b is non-zero and free of x, return 0 (when a free of x) or
     * 2 ArcTan[a/b] otherwise. */
    bool b_free = intrat_freeq_test(b, x);
    if (!is_zero_poly(b) && b_free) {
        if (intrat_freeq_test(a, x)) return expr_new_integer(0);
        Expr* ratio = internal_divide((Expr*[]){expr_copy(a), expr_copy(b)}, 2);
        Expr* canon = intrat_canonic(ratio); expr_free(ratio);
        Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
            (Expr*[]){canon}, 1);
        Expr* res = internal_times((Expr*[]){expr_new_integer(2), arctan}, 2);
        return eval_and_free(res);
    }

    /* A := Collect[a, x] (we use intrat_canonic for the closest
     * Mathilda equivalent — it folds rational coefficients to canonical
     * form without RootReduce).  The .m baseline applies simproot to
     * each coefficient, but Phase 3 deliberately stays minimal so the
     * tests can still match against the un-simplified output. */
    Expr* A = expr_expand(a);
    Expr* B = expr_expand(b);

    /* If A is divisible by B in K[x], LogToAtan reduces to a single
     * 2 ArcTan[A/B] term (the recursion would otherwise terminate in
     * one step). */
    if (!is_zero_poly(B)) {
        Expr* rem = intrat_polyrem_t(A, B, x);
        if (rem && is_zero_poly(rem)) {
            expr_free(rem);
            Expr* q = intrat_polyq(A, B, x);
            if (intrat_freeq_test(q, x)) {
                expr_free(q); expr_free(A); expr_free(B);
                return expr_new_integer(0);
            }
            Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
                (Expr*[]){q}, 1);
            Expr* res = internal_times((Expr*[]){expr_new_integer(2), arctan}, 2);
            expr_free(A); expr_free(B);
            return eval_and_free(res);
        }
        if (rem) expr_free(rem);
    }

    /* If deg(A) < deg(B), swap (recurse with -B, A). */
    int dA = get_degree_poly(A, x);
    int dB = get_degree_poly(B, x);
    if (dA < dB) {
        Expr* negB = expand_and_free(internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(B)}, 2));
        Expr* recurse = intrat_log_to_atan(negB, A, x);
        expr_free(negB); expr_free(A); expr_free(B);
        return recurse;
    }

    /* {g, {d, c}} = PolynomialExtendedGCD[B, -A, x]   so B*d - A*c = g. */
    Expr* negA = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(A)}, 2));
    Expr* eg_args[3] = { expr_copy(B), negA, expr_copy(x) };
    Expr* eg = internal_polynomialextendedgcd(eg_args, 3);
    Expr* eg_eval = evaluate(eg);
    expr_free(eg);
    if (!eg_eval || eg_eval->type != EXPR_FUNCTION
        || eg_eval->data.function.arg_count != 2) {
        if (eg_eval) expr_free(eg_eval);
        expr_free(A); expr_free(B);
        return NULL;
    }
    Expr* g_poly = expr_copy(eg_eval->data.function.args[0]);
    Expr* coefs = eg_eval->data.function.args[1];
    if (coefs->type != EXPR_FUNCTION || coefs->data.function.arg_count != 2) {
        expr_free(eg_eval); expr_free(g_poly);
        expr_free(A); expr_free(B); return NULL;
    }
    Expr* d_coef = expr_copy(coefs->data.function.args[0]);
    Expr* c_coef = expr_copy(coefs->data.function.args[1]);
    expr_free(eg_eval);

    /* v = collectnum[(A*d + B*c)/g, x] // Together. */
    Expr* Ad_raw = internal_times(
        (Expr*[]){expr_copy(A), expr_copy(d_coef)}, 2);
    Expr* Ad = expr_expand(Ad_raw); expr_free(Ad_raw);
    Expr* Bc_raw = internal_times(
        (Expr*[]){expr_copy(B), expr_copy(c_coef)}, 2);
    Expr* Bc = expr_expand(Bc_raw); expr_free(Bc_raw);
    Expr* sum_raw = internal_plus(
        (Expr*[]){Ad, Bc}, 2);
    Expr* div = internal_divide(
        (Expr*[]){sum_raw, g_poly}, 2);
    Expr* v = intrat_canonic(div); expr_free(div);

    expr_free(A); expr_free(B);

    /* If v is free of x, return LogToAtan[d, c, x] (no ArcTan
     * contribution at this level). */
    if (intrat_freeq_test(v, x)) {
        expr_free(v);
        Expr* recurse = intrat_log_to_atan(d_coef, c_coef, x);
        expr_free(d_coef); expr_free(c_coef);
        return recurse;
    }

    /* Else: 2 ArcTan[v] + LogToAtan[d, c, x]. */
    Expr* arctan = expr_new_function(expr_new_symbol(SYM_ArcTan),
        (Expr*[]){v}, 1);
    Expr* head_term = internal_times(
        (Expr*[]){expr_new_integer(2), arctan}, 2);
    Expr* recurse = intrat_log_to_atan(d_coef, c_coef, x);
    expr_free(d_coef); expr_free(c_coef);
    if (!recurse) recurse = expr_new_integer(0);
    Expr* sum = internal_plus(
        (Expr*[]){head_term, recurse}, 2);
    Expr* result = eval_and_free(sum);
    intrat_trace("LogToAtan", "OUT", result);
    return result;
}

/* ----- Shared helpers used by both the Phase 2 linear-Q closer and
 *       the Phase 4 LogToReal dispatcher. ----- */

/* Extract the (t - α) form of a linear polynomial in `t` and return
 * the rational root α.  Returns NULL when the input is not linear in
 * t.  Fresh Expr*. */
static Expr* extract_linear_root(Expr* poly, Expr* t) {
    int deg = get_degree_poly(poly, t);
    if (deg != 1) return NULL;
    Expr* a_coef = get_coeff(poly, t, 1);
    Expr* b_coef = get_coeff(poly, t, 0);
    Expr* neg_b = internal_times(
        (Expr*[]){expr_new_integer(-1), b_coef}, 2);
    Expr* alpha = internal_divide((Expr*[]){neg_b, a_coef}, 2);
    return eval_and_free(alpha);
}

/* Produce a single Log term: α * Log[S(α, x)]. */
static Expr* build_log_term(Expr* alpha, Expr* S, Expr* t) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_copy(t), expr_copy(alpha)}, 2);
    Expr* substituted = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){expr_copy(S), rule}, 2);
    Expr* Sat = evaluate(substituted);
    expr_free(substituted);
    Expr* logS = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){Sat}, 1);
    Expr* term = internal_times(
        (Expr*[]){expr_copy(alpha), logS}, 2);
    return eval_and_free(term);
}

/* ====================================================================
 * Phase 4 — LogToReal (Rioboo's complex-to-real conversion).
 * ====================================================================
 *
 * Given r in K[t] and s in K[t,x], produce a real function f satisfying
 *   D[f, x] == D[RootSum[r(α)==0, α Log[s(α, x)]], x].
 *
 * Phase 4 implements the bounded-Solve subset called out in the plan:
 * each factor of r over Q is dispatched as
 *   linear (t - α): contribute α Log[s(α, x)];
 *   quadratic a t^2 + b t + c with disc >= 0:
 *       two real roots (-b ± sqrt(disc))/(2a) treated as two linear cases;
 *   quadratic a t^2 + b t + c with disc < 0:
 *       complex conjugate pair u ± I v with u = -b/(2a),
 *       v = sqrt(-disc)/(2a); contribute
 *           u Log[A^2 + B^2] + v LogToAtan[A, B, x]
 *       where A + I B = s(u + I v, x) split into real / imaginary parts.
 *
 * Higher-degree factors (or symbolic discriminants we can't decide)
 * fall through to the symbolic RootSum form so the result still
 * differentiates correctly via the residue formula.
 */

/* Walk an expanded `term` and split it into (real_part, imag_part).
 * Recognises Complex[a, b] either as a bare factor or as a factor
 * inside Times.  Anything else is treated as fully real. */
static void split_term_re_im(Expr* term, Expr** re_out, Expr** im_out) {
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol.name == SYM_Complex
        && term->data.function.arg_count == 2) {
        *re_out = expr_copy(term->data.function.args[0]);
        *im_out = expr_copy(term->data.function.args[1]);
        return;
    }
    if (term->type == EXPR_FUNCTION
        && term->data.function.head->type == EXPR_SYMBOL
        && term->data.function.head->data.symbol.name == SYM_Times) {
        /* Look for a Complex factor among the Times args. */
        Expr* complex_factor = NULL;
        size_t complex_idx = 0;
        size_t n = term->data.function.arg_count;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = term->data.function.args[i];
            if (arg->type == EXPR_FUNCTION
                && arg->data.function.head->type == EXPR_SYMBOL
                && arg->data.function.head->data.symbol.name == SYM_Complex) {
                complex_factor = arg; complex_idx = i; break;
            }
        }
        if (!complex_factor) {
            *re_out = expr_copy(term);
            *im_out = expr_new_integer(0);
            return;
        }
        /* Build "rest" by dropping the complex factor. */
        size_t k = (n == 0) ? 0 : n - 1;
        Expr* rest;
        if (k == 0) rest = expr_new_integer(1);
        else if (k == 1) {
            rest = expr_copy(term->data.function.args[complex_idx == 0 ? 1 : 0]);
        } else {
            Expr** rest_args = (Expr**)malloc(sizeof(Expr*) * k);
            size_t kk = 0;
            for (size_t i = 0; i < n; i++) {
                if (i == complex_idx) continue;
                rest_args[kk++] = expr_copy(term->data.function.args[i]);
            }
            rest = eval_and_free(internal_times(rest_args, kk));
            free(rest_args);
        }
        Expr* a_part = expr_copy(complex_factor->data.function.args[0]);
        Expr* b_part = expr_copy(complex_factor->data.function.args[1]);
        *re_out = eval_and_free(internal_times(
            (Expr*[]){a_part, expr_copy(rest)}, 2));
        *im_out = eval_and_free(internal_times(
            (Expr*[]){b_part, rest}, 2));
        return;
    }
    /* No Complex anywhere — purely real. */
    *re_out = expr_copy(term);
    *im_out = expr_new_integer(0);
}

/* Expand p, then walk its Plus head, splitting each summand. */
static void split_re_im(Expr* p, Expr** re_out, Expr** im_out) {
    Expr* expanded = expr_expand(p);
    if (!(expanded->type == EXPR_FUNCTION
        && expanded->data.function.head->type == EXPR_SYMBOL
        && expanded->data.function.head->data.symbol.name == SYM_Plus)) {
        split_term_re_im(expanded, re_out, im_out);
        expr_free(expanded);
        return;
    }
    size_t n = expanded->data.function.arg_count;
    Expr** re_terms = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** im_terms = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nre = 0, nim = 0;
    for (size_t i = 0; i < n; i++) {
        Expr *re_t = NULL, *im_t = NULL;
        split_term_re_im(expanded->data.function.args[i], &re_t, &im_t);
        if (re_t && !is_zero_poly(re_t)) re_terms[nre++] = re_t;
        else if (re_t) expr_free(re_t);
        if (im_t && !is_zero_poly(im_t)) im_terms[nim++] = im_t;
        else if (im_t) expr_free(im_t);
    }
    *re_out = (nre == 0) ? expr_new_integer(0)
            : (nre == 1) ? re_terms[0]
            : eval_and_free(internal_plus(re_terms, nre));
    *im_out = (nim == 0) ? expr_new_integer(0)
            : (nim == 1) ? im_terms[0]
            : eval_and_free(internal_plus(im_terms, nim));
    free(re_terms); free(im_terms);
    expr_free(expanded);
}

/* Substitute t -> sub_expr inside e (full ReplaceAll). */
static Expr* subst_t(Expr* e, Expr* t, Expr* sub_expr) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_copy(t), expr_copy(sub_expr)}, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){expr_copy(e), rule}, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

/* ------------------------------------------------------------------ */
/* Arithmetic in Q(alpha) for a single Root-object coefficient.        */
/*                                                                    */
/* The pipeline's polynomial arithmetic (Together, resultants,         */
/* PolynomialExtendedGCD, ...) treats a Root object alpha = Root[f, k] */
/* as an opaque symbol, so nothing ever reduces alpha^n modulo f.  An  */
/* integrand that is a rational function of x over Q(alpha) therefore  */
/* came back with alpha^6..alpha^14 and ~10^160 coefficients that      */
/* cancel catastrophically under any fixed-precision N (the            */
/* Root[-1 - 4 #1 + #1^3 &, 2] Goursat repro), and the sign of a        */
/* discriminant in Q(alpha) was never decided (intrat_constant_sign).  */
/*                                                                    */
/* intrat_reduce_root_coefficients rewrites an expression over         */
/* K = Q(alpha) or K = Q(alpha)(Sqrt[d]) -- d in Q(alpha) the one       */
/* x-free square root the log part introduces for a split quadratic --  */
/* by exact identities only:                                           */
/*   - every maximal subexpression rational in (x, alpha, Sqrt[d])      */
/*     becomes N/D with N, D reduced modulo alpha's minimal polynomial  */
/*     and Sqrt[d]^2 = d; a D free of x is inverted in K, a D in x is   */
/*     made monic in x;                                                */
/*   - the summands of a Plus that share their non-K part are merged,   */
/*     their K coefficients summed and reduced.                        */
/* Only a single distinct Root object with a rational univariate       */
/* defining polynomial is handled; anything else is returned as is.    */
/* ------------------------------------------------------------------ */

static bool ir_is_root_object(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Root
        && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_FUNCTION
        && e->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.function.head->data.symbol.name
               == SYM_Function
        && e->data.function.args[1]->type == EXPR_INTEGER;
}

/* Record the distinct Root objects of `e` (stops at two). */
static void ir_collect_roots(const Expr* e, const Expr** found, int* n) {
    if (!e || e->type != EXPR_FUNCTION || *n > 1) return;
    if (ir_is_root_object(e)) {
        if (*n == 0) { found[0] = e; *n = 1; }
        else if (!expr_eq((Expr*)found[0], (Expr*)e)) *n = 2;
        return;
    }
    ir_collect_roots(e->data.function.head, found, n);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ir_collect_roots(e->data.function.args[i], found, n);
}

static bool ir_head_is(const Expr* e, const char* h) {
    return e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}

/* A rational expression: numbers and symbols under Plus, Times and
 * integer Power. */
static bool ir_is_rational(const Expr* e) {
    switch (e->type) {
    case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_SYMBOL: return true;
    case EXPR_FUNCTION: break;
    default: return false;
    }
    if (ir_head_is(e, SYM_Rational) || ir_head_is(e, SYM_Complex)) return true;
    if (ir_head_is(e, SYM_Power))
        return e->data.function.arg_count == 2
            && e->data.function.args[1]->type == EXPR_INTEGER
            && ir_is_rational(e->data.function.args[0]);
    if (!ir_head_is(e, SYM_Plus) && !ir_head_is(e, SYM_Times)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!ir_is_rational(e->data.function.args[i])) return false;
    return true;
}

typedef struct {
    Expr* z;   /* placeholder for alpha */
    Expr* p;   /* minimal (defining) polynomial of alpha, in z */
    Expr* x;   /* integration variable, or NULL */
    Expr* w;   /* placeholder for Sqrt[d], or NULL */
    Expr* d;   /* the radicand d, in z (when w is set) */
    Expr* params;  /* root / algebraic mode: List of the parameter symbols, or NULL */
} IrField;

/* True when every symbol of the rational expression `e` is a field
 * generator (z or w): `e` is then an element of K. */
static bool ir_in_field(const Expr* e, const IrField* F) {
    if (e->type == EXPR_SYMBOL)
        return expr_eq((Expr*)e, F->z) || (F->w && expr_eq((Expr*)e, F->w));
    if (e->type != EXPR_FUNCTION) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!ir_in_field(e->data.function.args[i], F)) return false;
    return true;
}

/* Normal form of a polynomial: Sqrt[d]^2 -> d, then modulo p(z).
 * Consumes e. */
static Expr* ir_rem(Expr* e, const IrField* F) {
    Expr* ex = expand_and_free(e);
    if (F->w && !intrat_freeq_test(ex, F->w)) {
        Expr* w2 = eval_and_free(internal_subtract((Expr*[]){
            internal_power((Expr*[]){expr_copy(F->w), expr_new_integer(2)}, 2),
            expr_copy(F->d)}, 2));
        Expr* r = eval_and_free(internal_polynomialremainder(
            (Expr*[]){ex, w2, expr_copy(F->w)}, 3));
        ex = expand_and_free(r);
    }
    return eval_and_free(internal_polynomialremainder(
        (Expr*[]){ex, expr_copy(F->p), expr_copy(F->z)}, 3));
}

/* 1/c for c in Q(z) (free of w) as a reduced polynomial, or NULL. */
static Expr* ir_inverse_z(Expr* c, const IrField* F) {
    Expr* eg = eval_and_free(internal_polynomialextendedgcd(
        (Expr*[]){expr_copy(c), expr_copy(F->p), expr_copy(F->z)}, 3));
    Expr* inv = NULL;
    if (eg && ir_head_is(eg, SYM_List) && eg->data.function.arg_count == 2) {
        Expr* g = eg->data.function.args[0];
        Expr* uv = eg->data.function.args[1];
        if (intrat_freeq_test(g, F->z) && !is_zero_poly(g) && ir_in_field(g, F)
            && ir_head_is(uv, SYM_List) && uv->data.function.arg_count == 2)
            inv = eval_and_free(internal_divide(
                (Expr*[]){expr_copy(uv->data.function.args[0]), expr_copy(g)}, 2));
    }
    if (eg) expr_free(eg);
    return inv;
}

/* 1/c for c in K, reduced, or NULL (c borrowed). */
static Expr* ir_inverse(Expr* c, const IrField* F) {
    if (!F->w || intrat_freeq_test(c, F->w)) return ir_inverse_z(c, F);
    /* 1/c = conj(c) / (c conj(c)), conj: Sqrt[d] -> -Sqrt[d]. */
    Expr* neg_w = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(F->w)}, 2));
    Expr* conj = subst_t(c, F->w, neg_w);
    expr_free(neg_w);
    Expr* norm = ir_rem(internal_times((Expr*[]){expr_copy(c), expr_copy(conj)}, 2), F);
    Expr* inv = NULL;
    if (norm && intrat_freeq_test(norm, F->w)) {
        Expr* ni = ir_inverse_z(norm, F);
        if (ni) {
            inv = ir_rem(internal_times((Expr*[]){conj, ni}, 2), F);
            conj = NULL;  /* consumed */
        }
    }
    if (conj) expr_free(conj);
    if (norm) expr_free(norm);
    return inv;
}

/* A reduced n0 + n1 w kept as n0 + w n1 with n1 unexpanded, so the
 * printed form carries Sqrt[d] once (consumes n). */
static Expr* ir_collect_w(Expr* n, const IrField* F) {
    if (!F->w || intrat_freeq_test(n, F->w)) return n;
    Expr* n0 = get_coeff(n, F->w, 0);
    Expr* n1 = get_coeff(n, F->w, 1);
    expr_free(n);
    return eval_and_free(internal_plus((Expr*[]){n0,
        internal_times((Expr*[]){expr_copy(F->w), n1}, 2)}, 2));
}

/* Split the rational content off the polynomial *p (FactorTermsList):
 * *p becomes the primitive part; returns the content. */
static Expr* ir_content(Expr** p) {
    Expr* ftl = call_eval("FactorTermsList", (Expr*[]){expr_copy(*p)}, 1);
    Expr* c = NULL;
    if (ftl && ir_head_is(ftl, SYM_List) && ftl->data.function.arg_count == 2) {
        Expr* c0 = ftl->data.function.args[0];
        if (c0->type == EXPR_INTEGER || c0->type == EXPR_BIGINT
            || ir_head_is(c0, SYM_Rational)) {
            c = expr_copy(c0);
            expr_free(*p);
            *p = expr_copy(ftl->data.function.args[1]);
        }
    }
    if (ftl) expr_free(ftl);
    return c ? c : expr_new_integer(1);
}

/* Reduce a rational expression in (x, z, w) to N/D over K. */
static Expr* ir_reduce_rational(Expr* e, const IrField* F) {
    Expr* tg = eval_and_free(internal_together((Expr*[]){expr_copy(e)}, 1));
    Expr* N = ir_rem(eval_and_free(intrat_numerator(tg)), F);
    Expr* D = ir_rem(eval_and_free(intrat_denominator(tg)), F);
    expr_free(tg);
    bool has_x = F->x && !intrat_freeq_test(D, F->x);
    if (has_x || !(D->type == EXPR_INTEGER || D->type == EXPR_BIGINT
                   || ir_head_is(D, SYM_Rational))) {
        /* Scale by the inverse of D itself (D free of x) or of its
         * leading coefficient in x, when that is an element of K. */
        Expr* lc = NULL;
        if (!has_x) lc = expr_copy(D);
        else {
            int dx = get_degree_poly(D, F->x);
            if (dx > 0) lc = get_coeff(D, F->x, dx);
        }
        Expr* inv = (lc && ir_in_field(lc, F)) ? ir_inverse(lc, F) : NULL;
        if (lc) expr_free(lc);
        if (inv) {
            N = ir_rem(internal_times((Expr*[]){N, expr_copy(inv)}, 2), F);
            if (has_x) D = ir_rem(internal_times((Expr*[]){D, inv}, 2), F);
            else { expr_free(D); expr_free(inv); D = expr_new_integer(1); }
        }
    }
    /* Rational content in front: c (N'/D'), so a Log[c p] downstream
     * sheds its constant and a kernel argument is primitive. */
    Expr* cN = ir_content(&N);
    Expr* cD = ir_content(&D);
    N = ir_collect_w(N, F);
    D = ir_collect_w(D, F);
    Expr* c = eval_and_free(internal_divide((Expr*[]){cN, cD}, 2));
    return eval_and_free(internal_times((Expr*[]){c,
        eval_and_free(internal_divide((Expr*[]){N, D}, 2))}, 2));
}

/* An element of K: rational, free of x, every symbol a generator. */
static bool ir_is_field_element(const Expr* e, const IrField* F) {
    return ir_is_rational(e) && ir_in_field(e, F)
        && (!F->x || intrat_freeq_test((Expr*)e, F->x));
}

/* Split a summand into its K coefficient and the rest (both fresh). */
static void ir_split_term(Expr* t, const IrField* F, Expr** coef, Expr** rest) {
    if (ir_is_field_element(t, F)) {
        *coef = expr_copy(t); *rest = expr_new_integer(1); return;
    }
    if (!ir_head_is(t, SYM_Times)) {
        *coef = expr_new_integer(1); *rest = expr_copy(t); return;
    }
    size_t n = t->data.function.arg_count;
    Expr** cs = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** rs = (Expr**)malloc(sizeof(Expr*) * n);
    size_t nc = 0, nr = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = t->data.function.args[i];
        if (ir_is_field_element(f, F)) cs[nc++] = expr_copy(f);
        else rs[nr++] = expr_copy(f);
    }
    *coef = nc ? eval_and_free(internal_times(cs, nc)) : expr_new_integer(1);
    *rest = nr ? eval_and_free(internal_times(rs, nr)) : expr_new_integer(1);
    free(cs); free(rs);
}

static Expr* ir_reduce_walk(Expr* e, const IrField* F);

/* Merge the summands of a Plus (args already reduced) by their non-K
 * part, reducing each summed K coefficient. */
static Expr* ir_merge_plus(Expr** args, size_t n, const IrField* F) {
    Expr** rests = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** coefs = (Expr**)malloc(sizeof(Expr*) * n);
    size_t ng = 0;
    for (size_t i = 0; i < n; i++) {
        Expr *c = NULL, *r = NULL;
        ir_split_term(args[i], F, &c, &r);
        size_t g = 0;
        while (g < ng && !expr_eq(rests[g], r)) g++;
        if (g < ng) {
            coefs[g] = eval_and_free(internal_plus((Expr*[]){coefs[g], c}, 2));
            expr_free(r);
        } else {
            rests[ng] = r; coefs[ng] = c; ng++;
        }
    }
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (ng ? ng : 1));
    for (size_t g = 0; g < ng; g++) {
        Expr* c = ir_reduce_rational(coefs[g], F);
        expr_free(coefs[g]);
        terms[g] = internal_times((Expr*[]){c, rests[g]}, 2);
    }
    Expr* out = ng ? eval_and_free(internal_plus(terms, ng)) : expr_new_integer(0);
    free(terms); free(rests); free(coefs);
    return out;
}

static Expr* ir_reduce_walk(Expr* e, const IrField* F) {
    if (intrat_freeq_test(e, F->z) && (!F->w || intrat_freeq_test(e, F->w)))
        return expr_copy(e);
    if (ir_is_rational(e)) return ir_reduce_rational(e, F);
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = ir_reduce_walk(e->data.function.args[i], F);
    Expr* out;
    if (ir_head_is(e, SYM_Plus)) {
        out = ir_merge_plus(args, n, F);
        for (size_t i = 0; i < n; i++) expr_free(args[i]);
    } else {
        out = eval_and_free(expr_new_function(
            expr_copy(e->data.function.head), args, n));
    }
    free(args);
    return out;
}

/* Set up Q(alpha) for the single Root object of `e`: alpha (borrowed
 * into *alpha), the placeholder z and its defining polynomial p(z).
 * Returns false when `e` has no Root object, more than one, or one
 * whose defining polynomial is not a rational univariate polynomial. */
static bool ir_field_setup(const Expr* e, const Expr** alpha, IrField* F) {
    const Expr* found[2] = { NULL, NULL };
    int n = 0;
    ir_collect_roots(e, found, &n);
    if (n != 1) return false;
    F->z = expr_new_symbol("Integrate`Private`rz$");
    F->x = NULL; F->w = NULL; F->d = NULL; F->params = NULL;
    Expr* app = expr_new_function(expr_copy(found[0]->data.function.args[0]),
        (Expr*[]){expr_copy(F->z)}, 1);
    F->p = expand_and_free(eval_and_free(app));
    if (!F->p || !ir_is_rational(F->p) || !ir_in_field(F->p, F)
        || !intrat_polyq_test(F->p, F->z) || get_degree_poly(F->p, F->z) < 2) {
        if (F->p) expr_free(F->p);
        expr_free(F->z);
        return false;
    }
    *alpha = found[0];
    return true;
}

static void ir_field_free(IrField* F) {
    expr_free(F->z); expr_free(F->p);
    if (F->params) { expr_free(F->params); F->params = NULL; }
    if (F->w) expr_free(F->w);
    if (F->d) expr_free(F->d);
}

/* Find the x-free square roots Power[d, +-1/2] with d in Q(z) in `e`;
 * *d_out gets the radicand when they all share one (stops at two). */
static void ir_collect_sqrt(const Expr* e, const IrField* F,
                            const Expr** d_out, int* n) {
    if (!e || e->type != EXPR_FUNCTION || *n > 1) return;
    if (ir_head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* ex = e->data.function.args[1];
        const Expr* b = e->data.function.args[0];
        if (ir_head_is(ex, SYM_Rational) && ex->data.function.arg_count == 2
            && ex->data.function.args[1]->type == EXPR_INTEGER
            && ex->data.function.args[1]->data.integer == 2
            && ex->data.function.args[0]->type == EXPR_INTEGER
            && (ex->data.function.args[0]->data.integer == 1
                || ex->data.function.args[0]->data.integer == -1)
            && !intrat_freeq_test((Expr*)b, F->z)
            && ir_is_rational(b) && ir_in_field(b, F)) {
            if (*n == 0) { *d_out = b; *n = 1; }
            else if (!expr_eq((Expr*)*d_out, (Expr*)b)) *n = 2;
            return;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ir_collect_sqrt(e->data.function.args[i], F, d_out, n);
}

/* ReplaceAll[e, {lhs1 -> rhs1, lhs2 -> rhs2}] (all borrowed). */
static Expr* ir_subst2(Expr* e, Expr* l1, Expr* r1, Expr* l2, Expr* r2) {
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){expr_copy(l1), expr_copy(r1)}, 2),
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){expr_copy(l2), expr_copy(r2)}, 2)}, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){expr_copy(e), rules}, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Reduce every K coefficient of `e` (borrowed) -- see the block comment
 * above; `x` (may be NULL) is the integration variable.  Returns a fresh
 * tree (a copy of `e` when nothing applies, or when the rewrite would be
 * larger). */
static Expr* intrat_reduce_root_coefficients(Expr* e, Expr* x) {
    const Expr* alpha = NULL;
    IrField F;
    if (!e || !ir_field_setup(e, &alpha, &F)) return e ? expr_copy(e) : NULL;
    F.x = x;
    Expr* ez = subst_t(e, (Expr*)alpha, F.z);
    /* One shared x-free Sqrt[d], d in Q(z), becomes the generator w. */
    const Expr* d = NULL;
    int nd = 0;
    ir_collect_sqrt(ez, &F, &d, &nd);
    if (nd == 1) {
        F.d = expr_copy((Expr*)d);
        F.w = expr_new_symbol("Integrate`Private`rw$");
        Expr* half = eval_and_free(internal_divide(
            (Expr*[]){expr_new_integer(1), expr_new_integer(2)}, 2));
        Expr* mhalf = eval_and_free(internal_divide(
            (Expr*[]){expr_new_integer(-1), expr_new_integer(2)}, 2));
        Expr* sq = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){expr_copy(F.d), half}, 2);
        Expr* isq = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){expr_copy(F.d), mhalf}, 2);
        Expr* w_over_d = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
            expr_copy(F.w), expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){expr_copy(F.d), expr_new_integer(-1)}, 2)}, 2);
        Expr* ezw = ir_subst2(ez, sq, F.w, isq, w_over_d);
        expr_free(sq); expr_free(isq); expr_free(w_over_d);
        expr_free(ez);
        ez = ezw;
    }
    Expr* red = ir_reduce_walk(ez, &F);
    expr_free(ez);
    Expr* out;
    if (F.w) {
        Expr* sqd = expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(F.d)}, 1);
        Expr* back = subst_t(red, F.w, sqd);
        expr_free(sqd);
        expr_free(red);
        red = back;
    }
    out = subst_t(red, F.z, (Expr*)alpha);
    expr_free(red);
    ir_field_free(&F);
    if (!out || leaf_count_internal(out, false) > leaf_count_internal(e, false)) {
        if (out) expr_free(out);
        return expr_copy(e);
    }
    return out;
}

/* Sign of an x-free constant: the positive-symbol rule first, then,
 * for a closed algebraic number in a single Root object (the
 * discriminant of a quadratic over Q(alpha)), the exact Sign[] of its
 * reduced form.  Before this the discriminant of such a quadratic was
 * always "unknown", so a quadratic with two real roots took the ArcTan
 * branch with an imaginary Sqrt[-disc]: an antiderivative with complex
 * values and branch jumps on the real line. */
static int intrat_constant_sign(Expr* c) {
    int s = intsimp_sign_pos_assumption(c);
    if (s != 0) return s;
    const Expr* alpha = NULL;
    IrField F;
    if (!ir_field_setup(c, &alpha, &F)) return 0;
    Expr* cz = subst_t(c, (Expr*)alpha, F.z);
    bool closed = ir_is_rational(cz) && ir_in_field(cz, &F)
        && !intrat_contains_complex(cz);
    expr_free(cz);
    ir_field_free(&F);
    if (!closed) return 0;
    Expr* red = intrat_reduce_root_coefficients(c, NULL);
    Expr* sg = call_eval(SYM_Sign, (Expr*[]){red}, 1);
    if (sg && sg->type == EXPR_INTEGER) {
        if (sg->data.integer > 0) s = 1;
        else if (sg->data.integer < 0) s = -1;
    }
    if (sg) expr_free(sg);
    return s;
}

/* ------------------------------------------------------------------ */
/* The Lazard-Rioboo-Trager log part over K = Q(alpha), alpha a Root.  */
/*                                                                    */
/* intrat_log_part_core treats every coefficient as an element of a   */
/* rational function field: a Root object alpha is an indeterminate,   */
/* nothing reduces alpha^n modulo its minimal polynomial, and every    */
/* gcd / Bezout cofactor in the residue variable t is computed over    */
/* Q(alpha)(t) with alpha transcendental.  Over Q(Root[837/4 + #^3 &,   */
/* 1]) the Rothstein-Trager resultant of a 92-leaf integrand came back */
/* with alpha^14 and 26-figure coefficients (it is 262144 t^6 - 1 once */
/* reduced), and the inverse of lc(S_1) modulo Q_1 -- a quintic in t   */
/* against a sextic -- did not finish in 60 s as a PolynomialExtended- */
/* GCD over Q(alpha, t): the rational functions in alpha of the        */
/* Euclidean chain swell without bound.  The answer lives in K[t]: in  */
/* root mode alpha is the placeholder z, every polynomial is reduced   */
/* modulo p(z), and every operation in t that depends on the field     */
/* structure (gcd, exact quotient, inverse modulo Q_i, the squarefree  */
/* split of lc(s_i)) runs natively over the number field               */
/* (FLINT`PolynomialExtendedGCD / FLINT`PolynomialQuotientRemainder    */
/* with Extension -> {z, p}).  The Q_i must be squarefree and pairwise */
/* coprime over K, which the generic FactorSquareFree over Q[z, t]     */
/* does not promise, so that is certified over K too (decline else).   */
/* ------------------------------------------------------------------ */

struct IrRootMode {
    IrField F;
    Expr* alpha;     /* the Root object (owned copy); NULL in algebraic mode */
    /* algebraic mode (irk_begin_alg): z is a primitive element theta of the
     * field of the coefficients' radicals / I / Root objects, and zb[j] is
     * theta^j written in the caller's own atoms (the product basis), so a
     * reduced polynomial in z goes back as sum_j c_j zb[j] */
    Expr** zb;
    int nzb;
};

/* Extension -> {z, p(z)} (fresh). */
static Expr* irk_ext_opt(const IrField* F) {
    return expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){
        expr_new_symbol(SYM_Extension),
        expr_new_function(expr_new_symbol(SYM_List),
            (Expr*[]){ expr_copy(F->z), expr_copy(F->p) }, 2) }, 2);
}

/* "Parameters" -> {x, params...} (fresh); NULL when there are none. */
static Expr* irk_params_opt(Expr* x, const IrField* F) {
    size_t np = F->params ? F->params->data.function.arg_count : 0;
    if (!x && np == 0) return NULL;
    Expr** ps = (Expr**)malloc(sizeof(Expr*) * (np + 1));
    size_t n = 0;
    if (x) ps[n++] = expr_copy(x);
    for (size_t i = 0; i < np; i++) ps[n++] = expr_copy(F->params->data.function.args[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Rule), (Expr*[]){
        expr_new_string("Parameters"),
        expr_new_function(expr_new_symbol(SYM_List), ps, n) }, 2);
    free(ps);
    return r;
}

/* The symbols of a and d other than x and z (the parameters of Q(params)):
 * every FLINT kernel call names them, so the field arithmetic runs over
 * K(params) instead of declining. */
static void irk_set_params(IrField* F, Expr* a, Expr* d, Expr* x) {
    size_t vc = 0, vcap = 8;
    Expr** vs = (Expr**)malloc(sizeof(Expr*) * vcap);
    collect_variables(a, &vs, &vc, &vcap);
    collect_variables(d, &vs, &vc, &vcap);
    Expr** keep = (Expr**)malloc(sizeof(Expr*) * (vc ? vc : 1));
    size_t nk = 0;
    for (size_t i = 0; i < vc; i++) {
        bool dup = false;
        for (size_t j = 0; j < nk && !dup; j++) dup = expr_eq(keep[j], vs[i]);
        if (vs[i]->type == EXPR_SYMBOL && !dup && !expr_eq(vs[i], x) && !expr_eq(vs[i], F->z))
            keep[nk++] = expr_copy(vs[i]);
    }
    for (size_t i = 0; i < vc; i++) expr_free(vs[i]);
    free(vs);
    if (F->params) { expr_free(F->params); F->params = NULL; }
    if (nk) F->params = expr_new_function(expr_new_symbol(SYM_List), keep, nk);
    free(keep);
}

/* FLINT`head[a, b, t, Extension -> {z, p}(, "Parameters" -> {x})]: the
 * evaluated List, or NULL when the kernel declined (left unevaluated,
 * e.g. FLINT compiled out, a reducible p, a coefficient outside K). */
static Expr* irk_call(const char* head, Expr* a, Expr* b, Expr* t,
                      const IrField* F, Expr* x) {
    Expr* args[5] = { expr_copy(a), expr_copy(b), expr_copy(t), irk_ext_opt(F), NULL };
    size_t n = 4;
    Expr* po = irk_params_opt(x, F);
    if (po) args[n++] = po;
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(head), args, n));
    if (r && ir_head_is(r, SYM_List) && r->data.function.arg_count == 2) return r;
    if (r) expr_free(r);
    return NULL;
}

/* Monic gcd over K of two polynomials in t free of x, or NULL. */
static Expr* irk_gcd(Expr* a, Expr* b, Expr* t, const IrField* F) {
    Expr* r = irk_call(SYM_FLINT_PolynomialExtendedGCD, a, b, t, F, NULL);
    if (!r) return NULL;
    Expr* g = eval_and_free(expr_copy(r->data.function.args[0]));
    expr_free(r);
    return g;
}

/* Division in K[t] (the coefficients may carry x, read as a parameter):
 * quotient and remainder, evaluated; false when the kernel declined. */
static bool irk_divrem(Expr* a, Expr* b, Expr* t, const IrField* F, Expr* x,
                       Expr** q, Expr** rem) {
    *q = *rem = NULL;
    Expr* r = irk_call(SYM_FLINT_PolynomialQuotientRemainder, a, b, t, F, x);
    if (!r) return false;
    *q = eval_and_free(expr_copy(r->data.function.args[0]));
    *rem = eval_and_free(expr_copy(r->data.function.args[1]));
    expr_free(r);
    return true;
}

static bool irk_is_one(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* Exact quotient a / b in K[t]; NULL unless the remainder is 0. */
static Expr* irk_exquo(Expr* a, Expr* b, Expr* t, const IrField* F, Expr* x) {
    Expr *q, *r;
    if (!irk_divrem(a, b, t, F, x, &q, &r)) return NULL;
    bool exact = is_zero_poly(r);
    expr_free(r);
    if (!exact) { expr_free(q); return NULL; }
    return q;
}

/* Yun's squarefree decomposition over K of p in K[t] (free of x): a
 * List {A_1, A_2, ...}, A_j the product of the factors of multiplicity
 * j (an A_j may be 1), or NULL. */
static Expr* irk_squarefree(Expr* p, Expr* t, const IrField* F) {
    Expr* dp = expand_and_free(intrat_d(p, t));
    if (is_zero_poly(dp)) {                  /* p constant in t */
        expr_free(dp);
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    Expr* c = irk_gcd(p, dp, t, F);
    expr_free(dp);
    if (!c) return NULL;
    Expr* w = irk_exquo(p, c, t, F, NULL);
    if (!w) { expr_free(c); return NULL; }
    Expr** out = NULL;
    size_t n = 0, cap = 0;
    bool ok = true;
    while (get_degree_poly(w, t) > 0) {
        Expr* y = irk_gcd(w, c, t, F);
        Expr* aj = y ? irk_exquo(w, y, t, F, NULL) : NULL;
        Expr* cn = y ? irk_exquo(c, y, t, F, NULL) : NULL;
        if (!y || !aj || !cn) {
            if (y) expr_free(y);
            if (aj) expr_free(aj);
            if (cn) expr_free(cn);
            ok = false;
            break;
        }
        if (n == cap) { cap = cap ? 2 * cap : 4; out = (Expr**)realloc(out, sizeof(Expr*) * cap); }
        out[n++] = aj;
        expr_free(w); w = y;
        expr_free(c); c = cn;
    }
    expr_free(w); expr_free(c);
    if (!ok) {
        for (size_t i = 0; i < n; i++) expr_free(out[i]);
        free(out);
        return NULL;
    }
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);
    return lst;
}

/* Enter root mode for the log part of a/d: one Root object alpha with a
 * rational univariate defining polynomial among the coefficients.  On
 * success *a and *d are rewritten over z (alpha -> z, reduced modulo
 * p) and the mode is returned; NULL leaves *a and *d untouched. */
static IrRootMode* irk_begin(Expr** a, Expr** d, Expr* x) {
    Expr* ad = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_copy(*a), expr_copy(*d) }, 2);
    const Expr* alpha = NULL;
    IrField F;
    if (!ir_field_setup(ad, &alpha, &F)) { expr_free(ad); return NULL; }
    IrRootMode* M = (IrRootMode*)malloc(sizeof(IrRootMode));
    M->F = F;
    M->alpha = expr_copy((Expr*)alpha);
    M->zb = NULL; M->nzb = 0;
    expr_free(ad);
    /* The kernels need p irreducible (the minimal polynomial of a
     * canonical Root object) and FLINT present; probe once, on a
     * coprime pair (the multi-modular xgcd is slow on proportional
     * inputs -- 70 ms for gcd(x, x)). */
    Expr* one = expr_new_integer(1);
    Expr* g = irk_gcd(x, one, x, &M->F);
    expr_free(one);
    if (!g) {
        expr_free(M->alpha);
        ir_field_free(&M->F);
        free(M);
        return NULL;
    }
    expr_free(g);
    M->F.x = x;
    Expr* az = subst_t(*a, M->alpha, M->F.z);
    Expr* dz = subst_t(*d, M->alpha, M->F.z);
    expr_free(*a); expr_free(*d);
    *a = ir_rem(az, &M->F);
    *d = ir_rem(dz, &M->F);
    irk_set_params(&M->F, *a, *d, x);
    return M;
}

#ifdef USE_FLINT
/* AlgebraicNumber[theta, {c0, ..}] -> c0 + c1 z + ... (every other node
 * copied). */
static Expr* irk_algnum_to_z(const Expr* e, Expr* z) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (ir_head_is(e, SYM_AlgebraicNumber) && e->data.function.arg_count == 2 &&
        ir_head_is(e->data.function.args[1], SYM_List)) {
        const Expr* co = e->data.function.args[1];
        size_t k = co->data.function.arg_count;
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * (k ? k : 1));
        for (size_t i = 0; i < k; i++)
            terms[i] = internal_times((Expr*[]){ expr_copy(co->data.function.args[i]),
                internal_power((Expr*[]){ expr_copy(z), expr_new_integer((int64_t)i) }, 2) }, 2);
        Expr* r = k ? internal_plus(terms, k) : expr_new_integer(0);
        free(terms);
        return r;
    }
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = irk_algnum_to_z(e->data.function.args[i], z);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}

static int irk_cmp_exprs(const void* a, const void* b) {
    return expr_compare(*(Expr* const*)a, *(Expr* const*)b);
}
#endif

/* Enter algebraic mode for the log part of a/d: the coefficients carry
 * constant algebraic atoms (I, Sqrt[3], 2^(1/3), (-1)^(1/3), Root objects)
 * spanning K = Q(theta).  Each atom is replaced by its image, a polynomial in
 * z standing for theta; the answer is rendered back over the product basis of
 * the atoms.  The generic core treats every atom as an independent
 * indeterminate that nothing reduces, so its subresultant chain swells
 * without bound (a degree-12 denominator over Q(I Sqrt[3]) ran past a
 * minute).  NULL leaves *a and *d untouched. */
static IrRootMode* irk_begin_alg(Expr** a, Expr** d, Expr* x) {
#ifdef USE_FLINT
    Expr** atoms = NULL;
    size_t na = 0, cap = 0;
    if (!fb_collect_atoms(*a, &atoms, &na, &cap) || !fb_collect_atoms(*d, &atoms, &na, &cap) ||
        na == 0 || na > FG_RENDER_MAX_ATOMS) {
        for (size_t i = 0; i < na; i++) expr_free(atoms[i]);
        free(atoms);
        return NULL;
    }
    if (na > 1) qsort(atoms, na, sizeof(Expr*), irk_cmp_exprs);
    IrRootMode* M = NULL;
    Expr* imgs = NULL;
    Expr* mp = NULL;
    {
        const Expr* ci = fb_field_images(atoms, na);
        if (ci) imgs = expr_copy((Expr*)ci);       /* the cache slot may be recycled */
    }
    const Expr* theta = NULL;
    bool ok = imgs != NULL;
    for (size_t i = 0; i < na && ok; i++) {
        const Expr* im = imgs->data.function.args[i];
        ok = ir_head_is(im, SYM_AlgebraicNumber) && im->data.function.arg_count == 2;
        if (!ok) break;
        if (!theta) theta = im->data.function.args[0];
        else ok = expr_eq((Expr*)theta, im->data.function.args[0]);
    }
    if (ok) {
        mp = flint_qqbar_gen_minpoly_coeffs(theta);
        ok = mp && ir_head_is(mp, SYM_List) && mp->data.function.arg_count >= 3;
    }
    if (ok) {
        slong n = (slong)mp->data.function.arg_count - 1;
        FgRender R;
        fb_render_init(&R, atoms, imgs->data.function.args, na, theta, n);
        ok = R.ready != 0;
        M = (IrRootMode*)calloc(1, sizeof(IrRootMode));
        M->F.z = expr_new_symbol("Integrate`Private`rz$");
        Expr** pt = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
        for (slong j = 0; j <= n; j++)
            pt[j] = internal_times((Expr*[]){ expr_copy(mp->data.function.args[j]),
                internal_power((Expr*[]){ expr_copy(M->F.z), expr_new_integer((int64_t)j) }, 2) }, 2);
        M->F.p = expand_and_free(eval_and_free(internal_plus(pt, (size_t)(n + 1))));
        free(pt);
        M->zb = (Expr**)calloc((size_t)n, sizeof(Expr*));
        M->nzb = (int)n;
        for (slong j = 0; j < n && ok; j++) {
            Expr** co = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
            for (slong i = 0; i < n; i++) co[i] = expr_new_integer(i == j ? 1 : 0);
            Expr* an = expr_new_function(expr_new_symbol(SYM_AlgebraicNumber), (Expr*[]){
                expr_copy((Expr*)theta), expr_new_function(expr_new_symbol(SYM_List), co, (size_t)n) }, 2);
            free(co);
            Expr* rendered = fb_render_tree(an, &R);
            expr_free(an);
            M->zb[j] = rendered ? eval_and_free(rendered) : NULL;
            ok = M->zb[j] != NULL && intrat_freeq_test(M->zb[j], x);
        }
        fb_render_clear(&R);
        if (ok) {
            /* the kernels need FLINT and p irreducible: probe once */
            Expr* one = expr_new_integer(1);
            Expr* g = irk_gcd(x, one, x, &M->F);
            expr_free(one);
            ok = g != NULL;
            if (g) expr_free(g);
        }
    }
    if (ok) {
        Expr* aa = fb_subst_atoms(*a, atoms, imgs->data.function.args, na);
        Expr* dd = fb_subst_atoms(*d, atoms, imgs->data.function.args, na);
        Expr* az = aa ? eval_and_free(irk_algnum_to_z(aa, M->F.z)) : NULL;
        Expr* dz = dd ? eval_and_free(irk_algnum_to_z(dd, M->F.z)) : NULL;
        if (aa) expr_free(aa);
        if (dd) expr_free(dd);
        ok = az && dz;
        if (ok) {
            M->F.x = x;
            expr_free(*a); expr_free(*d);
            *a = ir_rem(az, &M->F);
            *d = ir_rem(dz, &M->F);
            irk_set_params(&M->F, *a, *d, x);
        } else {
            if (az) expr_free(az);
            if (dz) expr_free(dz);
        }
    }
    if (imgs) expr_free(imgs);
    if (mp) expr_free(mp);
    for (size_t i = 0; i < na; i++) expr_free(atoms[i]);
    free(atoms);
    if (!ok && M) {
        M->F.x = NULL;
        if (M->zb) {
            for (int j = 0; j < M->nzb; j++) if (M->zb[j]) expr_free(M->zb[j]);
            free(M->zb);
        }
        ir_field_free(&M->F);
        free(M);
        M = NULL;
    }
    return M;
#else
    (void)a; (void)d; (void)x;
    return NULL;
#endif
}

static void irk_end(IrRootMode* M) {
    if (!M) return;
    if (M->alpha) expr_free(M->alpha);
    if (M->zb) {
        for (int j = 0; j < M->nzb; j++) expr_free(M->zb[j]);
        free(M->zb);
    }
    M->F.x = NULL;
    ir_field_free(&M->F);
    free(M);
}

/* Normal form modulo p(z) (consumes e). */
static Expr* irk_reduce(IrRootMode* M, Expr* e) {
    return ir_rem(e, &M->F);
}

/* z -> alpha (consumes e); in algebraic mode z^j -> zb[j] after reduction
 * modulo p(z), so the answer is written in the caller's radicals. */
static Expr* irk_back(IrRootMode* M, Expr* e) {
    if (M->zb) {
        Expr* r = ir_rem(e, &M->F);
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)M->nzb);
        size_t nt = 0;
        for (int j = 0; j < M->nzb; j++) {
            Expr* c = get_coeff(r, M->F.z, (int)j);
            if (!c) continue;
            if (is_zero_poly(c)) { expr_free(c); continue; }
            terms[nt++] = internal_times((Expr*[]){ c, expr_copy(M->zb[j]) }, 2);
        }
        expr_free(r);
        Expr* out = nt == 0 ? expr_new_integer(0) : internal_plus(terms, nt);
        free(terms);
        return expand_and_free(eval_and_free(out));
    }
    Expr* r = subst_t(e, M->F.z, M->alpha);
    expr_free(e);
    return r;
}

static bool irk_factors_ok(IrRootMode* M, Expr* Q, Expr* t);

/* The residue factors {{Q_1, 1}, {Q_2, 2}, ...} reduced over K (a fresh
 * list; Q is borrowed and not touched -- expr_copy shares nodes), or
 * NULL unless every non-constant Q_i is squarefree over K and the Q_i
 * are pairwise coprime over K (the LRT precondition). */
static Expr* irk_certified_factors(IrRootMode* M, Expr* Qin, Expr* t) {
    size_t n = Qin->data.function.arg_count;
    Expr** prs = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr* pr = Qin->data.function.args[i];
        if (!ir_head_is(pr, SYM_List) || pr->data.function.arg_count != 2) {
            for (size_t k = 0; k < i; k++) expr_free(prs[k]);
            free(prs);
            return NULL;
        }
        prs[i] = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
            irk_reduce(M, expr_copy(pr->data.function.args[0])),
            expr_copy(pr->data.function.args[1]) }, 2);
    }
    Expr* Q = expr_new_function(expr_new_symbol(SYM_List), prs, n);
    free(prs);
    if (!irk_factors_ok(M, Q, t)) { expr_free(Q); return NULL; }
    return Q;
}

static bool irk_factors_ok(IrRootMode* M, Expr* Q, Expr* t) {
    size_t n = Q->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        Expr* Qi = Q->data.function.args[i]->data.function.args[0];
        if (get_degree_poly(Qi, t) <= 0) continue;
        Expr* dQ = expand_and_free(intrat_d(Qi, t));
        Expr* g = irk_gcd(Qi, dQ, t, &M->F);
        expr_free(dQ);
        bool sf = g && irk_is_one(g);
        if (g) expr_free(g);
        if (!sf) return false;
        for (size_t j = i + 1; j < n; j++) {
            Expr* Qj = Q->data.function.args[j]->data.function.args[0];
            if (get_degree_poly(Qj, t) <= 0) continue;
            Expr* h = irk_gcd(Qi, Qj, t, &M->F);
            bool cop = h && irk_is_one(h);
            if (h) expr_free(h);
            if (!cop) return false;
        }
    }
    return true;
}

/* PrimitivePart_mod over K: pp_x(inv(lc_x(p)) p mod Q), the inverse
 * taken in K[t]/(Q).  NULL on failure (lc not invertible modulo Q). */
static Expr* irk_primitive_part_mod(IrRootMode* M, Expr* p, Expr* Q,
                                    Expr* x, Expr* t) {
    Expr* lc = intrat_lc(p, x);
    Expr *q0, *lc_r;
    bool ok = irk_divrem(lc, Q, t, &M->F, NULL, &q0, &lc_r);
    expr_free(lc);
    if (!ok) return NULL;
    expr_free(q0);
    Expr* eg = irk_call(SYM_FLINT_PolynomialExtendedGCD, lc_r, Q, t, &M->F, NULL);
    expr_free(lc_r);
    if (!eg) return NULL;
    Expr* g = eg->data.function.args[0];
    Expr* uv = eg->data.function.args[1];
    if (!irk_is_one(g) || !ir_head_is(uv, SYM_List) || uv->data.function.arg_count != 2) {
        expr_free(eg);
        return NULL;
    }
    Expr* inv = eval_and_free(expr_copy(uv->data.function.args[0]));
    expr_free(eg);
    Expr* rp = expand_and_free(internal_times((Expr*[]){ inv, expr_copy(p) }, 2));
    Expr *q1, *rem;
    ok = irk_divrem(rp, Q, t, &M->F, x, &q1, &rem);
    expr_free(rp);
    if (!ok) return NULL;
    expr_free(q1);
    Expr* out = intrat_primitive(rem, x);
    expr_free(rem);
    return out;
}

/* S_i of the LRT over K for the residue factor Qi (squarefree over K):
 * the subresultant s of degree i in x (d itself when i == deg d), its
 * factors in common with Qi divided out with their multiplicity in
 * lc_x(s), then PrimitivePart_mod.  NULL on failure. */
static Expr* irk_lrt_S(IrRootMode* M, Expr* prs, Expr* d, Expr* Qi,
                       int i, int deg_d, Expr* x, Expr* t) {
    if (i == deg_d) return irk_primitive_part_mod(M, d, Qi, x, t);
    Expr* s_raw = find_prs_at_degree(prs, x, i);
    if (!s_raw) return NULL;
    Expr* s_red = irk_reduce(M, s_raw);
    Expr* s = intrat_primitive(s_red, x);
    expr_free(s_red);
    Expr* lc_s = intrat_lc(s, x);
    Expr* A = irk_squarefree(lc_s, t, &M->F);
    expr_free(lc_s);
    if (!A) { expr_free(s); return NULL; }
    for (size_t j = 0; j < A->data.function.arg_count; j++) {
        Expr* Aj = A->data.function.args[j];
        if (get_degree_poly(Aj, t) <= 0) continue;
        Expr* g = irk_gcd(Aj, Qi, t, &M->F);
        if (!g) { expr_free(A); expr_free(s); return NULL; }
        if (get_degree_poly(g, t) <= 0) { expr_free(g); continue; }
        Expr* gj = expand_and_free(internal_power(
            (Expr*[]){ g, expr_new_integer((int64_t)(j + 1)) }, 2));
        Expr* sn = irk_exquo(s, gj, t, &M->F, x);
        expr_free(gj);
        if (!sn) { expr_free(A); expr_free(s); return NULL; }
        expr_free(s);
        s = sn;
    }
    expr_free(A);
    Expr* S = irk_primitive_part_mod(M, s, Qi, x, t);
    expr_free(s);
    return S;
}

/* The pure-rational log part of a/d over K = Q(alpha) (root mode, see
 * above).  Same output contract as intrat_log_part_core (a RootSum sum,
 * or the {{Q_i, S_i}, ...} pair list), with Q_i and S_i reduced over K.
 * *handled is false -- a and d untouched, the caller runs the generic
 * core -- when the coefficients do not span such a K, a kernel declines
 * or the residue factors are not squarefree / coprime over K.  When
 * *handled is true a and d have been consumed. */
static Expr* intrat_log_part_core_root(Expr* a, Expr* d, Expr* x, Expr* t,
                                       bool root_sum, bool* handled) {
    *handled = false;
    Expr* az = expr_copy(a);
    Expr* dz = expr_copy(d);
    IrRootMode* M = irk_begin(&az, &dz, x);
    /* not one Root object: the field of the radicals, I and Root objects.
     * The generic chain is cheap for a short one and the field kernels cost a
     * fixed ~0.1 s per call sequence over a degree-4 field (a quartic over
     * Q(I, Sqrt[3]): 0.47 s here, 30 ms generic), so the field is entered for
     * the cases the generic chain cannot carry: a long chain (deg_x d >= 5:
     * a degree-12 denominator over Q(I Sqrt[3]) ran past a minute) or
     * parameters beside the radicals (u/(a + b u^3) over Q(k)(I Sqrt[3]):
     * 0.56 s generic, 0.19 s here). */
    if (!M) {
        size_t nv = 0, vcap = 8;
        Expr** vs = (Expr**)malloc(sizeof(Expr*) * vcap);
        collect_variables(az, &vs, &nv, &vcap);
        collect_variables(dz, &vs, &nv, &vcap);
        bool params = false;
        for (size_t i = 0; i < nv; i++) {
            if (vs[i]->type == EXPR_SYMBOL && !expr_eq(vs[i], x)) params = true;
            expr_free(vs[i]);
        }
        free(vs);
        if (params || get_degree_poly(dz, x) >= 5) M = irk_begin_alg(&az, &dz, x);
    }
    if (!M) { expr_free(az); expr_free(dz); return NULL; }

    Expr* dprime = expand_and_free(intrat_d(dz, x));
    Expr* atdp = irk_reduce(M, internal_subtract((Expr*[]){ expr_copy(az),
        internal_times((Expr*[]){ expr_copy(t), dprime }, 2) }, 2));
    /* the genuine subresultant chain over Q[z, t][x] (polynomial size), the
     * pseudo-remainder chain only when that declines */
    Expr* prs = intrat_subres_chain_mv(dz, atdp, x, t);
    if (!prs) prs = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_SubresultantPolynomialRemainders),
        (Expr*[]){ expr_copy(dz), expr_copy(atdp), expr_copy(x) }, 3));
    Expr* res_raw = eval_and_free(expr_new_function(expr_new_symbol(SYM_Resultant),
        (Expr*[]){ expr_copy(dz), atdp, expr_copy(x) }, 3));

    /* Q: the squarefree split of the resultant as computed (alpha
     * unreduced); Qred: the same factors reduced over K, certified, and
     * the ones every computation in t uses.  The output takes whichever
     * spelling of each factor is smaller: the two are equal in K[t], and
     * 256 alpha^6 t^4 - 1 is a smaller radicand than its reduced form. */
    Expr* Q = NULL;
    Expr* Qred = NULL;
    Expr** S = NULL;
    size_t nQ = 0;
    bool ok = prs && res_raw && ir_head_is(prs, SYM_List);
    if (ok) {
        Expr* resultant = intrat_primitive(res_raw, t);
        expr_free(res_raw);
        res_raw = NULL;
        Q = intrat_squarefree_list(resultant);
        ok = Q && ir_head_is(Q, SYM_List);
        if (ok) {
            Qred = irk_certified_factors(M, Q, t);
            ok = Qred != NULL;
        }
        /* The squarefree split over Q[z, t] reads z as an indeterminate, so
         * after reduction modulo p(z) a factor may repeat (or two share a
         * root) over K: split the reduced resultant over K itself (Yun, with
         * the field gcd) instead of giving up. */
        if (!ok) {
            if (Q) { expr_free(Q); Q = NULL; }
            Expr* rr = irk_reduce(M, expr_copy(resultant));
            Expr* A = irk_squarefree(rr, t, &M->F);
            expr_free(rr);
            if (A && A->data.function.arg_count > 0) {
                size_t na = A->data.function.arg_count;
                Expr** pr = (Expr**)malloc(sizeof(Expr*) * na);
                for (size_t j = 0; j < na; j++)
                    pr[j] = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
                        expr_copy(A->data.function.args[j]), expr_new_integer((int64_t)(j + 1)) }, 2);
                Q = expr_new_function(expr_new_symbol(SYM_List), pr, na);
                free(pr);
                Qred = irk_certified_factors(M, Q, t);
                ok = Qred != NULL;
            }
            if (A) expr_free(A);
        }
        expr_free(resultant);
    }
    if (ok) {
        int deg_d = get_degree_poly(dz, x);
        nQ = Q->data.function.arg_count;
        S = (Expr**)calloc(nQ ? nQ : 1, sizeof(Expr*));
        for (size_t k = 0; k < nQ && ok; k++) {
            Expr* Qi = Qred->data.function.args[k]->data.function.args[0];
            if (get_degree_poly(Qi, t) <= 0) continue;
            S[k] = irk_lrt_S(M, prs, dz, Qi, (int)(k + 1), deg_d, x, t);
            if (!S[k]) ok = false;
        }
    }
    if (res_raw) expr_free(res_raw);
    if (prs) expr_free(prs);
    expr_free(az); expr_free(dz);
    if (!ok) {
        if (S) { for (size_t k = 0; k < nQ; k++) if (S[k]) expr_free(S[k]); free(S); }
        if (Q) expr_free(Q);
        if (Qred) expr_free(Qred);
        irk_end(M);
        return NULL;
    }

    /* z -> alpha, and the output in the generic core's two shapes. */
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (nQ ? nQ : 1));
    size_t nt = 0;
    for (size_t k = 0; k < nQ; k++) {
        if (!S[k]) continue;
        Expr* Qi = eval_and_free(irk_back(M,
            expr_copy(Qred->data.function.args[k]->data.function.args[0])));
        Expr* Qo = eval_and_free(irk_back(M,
            expr_copy(Q->data.function.args[k]->data.function.args[0])));
        if (leaf_count_internal(Qo, false) < leaf_count_internal(Qi, false)) {
            expr_free(Qi);
            Qi = Qo;
        } else {
            expr_free(Qo);
        }
        Expr* Sk = eval_and_free(irk_back(M, S[k]));
        S[k] = NULL;
        if (root_sum) {
            Expr* func1 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){ expr_copy(t), Qi }, 2);
            Expr* tlog = eval_and_free(internal_times((Expr*[]){ expr_copy(t),
                expr_new_function(expr_new_symbol(SYM_Log), (Expr*[]){ Sk }, 1) }, 2));
            Expr* func2 = expr_new_function(expr_new_symbol(SYM_Function),
                (Expr*[]){ expr_copy(t), tlog }, 2);
            terms[nt++] = expr_new_function(expr_new_symbol(SYM_RootSum),
                (Expr*[]){ func1, func2 }, 2);
        } else {
            terms[nt++] = expr_new_function(expr_new_symbol(SYM_List),
                (Expr*[]){ Qi, Sk }, 2);
        }
    }
    free(S);
    expr_free(Q);
    expr_free(Qred);
    irk_end(M);
    Expr* result;
    if (!root_sum) {
        result = expr_new_function(expr_new_symbol(SYM_List), terms, nt);
    } else if (nt == 0) {
        result = expr_new_integer(0);
    } else if (nt == 1) {
        result = terms[0];
    } else {
        result = eval_and_free(internal_plus(terms, nt));
    }
    free(terms);
    expr_free(a); expr_free(d);
    *handled = true;
    intrat_trace("IntRationalLogPart", "OUT (over Q(Root))", result);
    return result;
}

/* Sign tests, canonic-zero, and the positive-symbol Sqrt simplifier
 * used by logtoreal_dispatch live in intsimp.c — see intsimp.h. */

/* True when `e` contains an explicit Complex[..] number (e.g. I). */
static bool intrat_contains_complex(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex) return true;
    if (intrat_contains_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_contains_complex(e->data.function.args[i])) return true;
    return false;
}

/* Try to dispatch a single quadratic factor `Q(t) = a t^2 + b t + c`
 * to its real-form contribution.  Returns NULL on failure (caller
 * falls back to symbolic RootSum). */
static Expr* logtoreal_quadratic_core(Expr* a, Expr* b, Expr* c,
                                      Expr* s, Expr* x, Expr* t,
                                      bool over_root);

/* Over Q(alpha) for a Root object alpha the quadratic is first made
 * monic, with b/a and c/a reduced modulo alpha's minimal polynomial,
 * so the discriminant and the residues stay small. */
static Expr* logtoreal_quadratic(Expr* a, Expr* b, Expr* c,
                                 Expr* s, Expr* x, Expr* t) {
    Expr* abc = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){expr_copy(a), expr_copy(b), expr_copy(c)}, 3);
    const Expr* alpha = NULL;
    IrField F;
    bool over_root = ir_field_setup(abc, &alpha, &F);
    expr_free(abc);
    if (!over_root) return logtoreal_quadratic_core(a, b, c, s, x, t, false);
    ir_field_free(&F);
    Expr* bq = eval_and_free(internal_divide((Expr*[]){expr_copy(b), expr_copy(a)}, 2));
    Expr* cq = eval_and_free(internal_divide((Expr*[]){expr_copy(c), expr_copy(a)}, 2));
    Expr* bm = intrat_reduce_root_coefficients(bq, NULL);
    Expr* cm = intrat_reduce_root_coefficients(cq, NULL);
    expr_free(bq); expr_free(cq);
    Expr* one = expr_new_integer(1);
    Expr* r = logtoreal_quadratic_core(one, bm, cm, s, x, t, true);
    expr_free(one); expr_free(bm); expr_free(cm);
    return r;
}

static Expr* logtoreal_quadratic_core(Expr* a, Expr* b, Expr* c,
                                      Expr* s, Expr* x, Expr* t,
                                      bool over_root) {
    /* discriminant = b^2 - 4 a c. */
    Expr* b2 = expand_and_free(internal_power((Expr*[]){expr_copy(b), expr_new_integer(2)}, 2));
    Expr* fourac = expand_and_free(
        internal_times((Expr*[]){expr_new_integer(4), expr_copy(a), expr_copy(c)}, 3));
    Expr* neg_fourac = internal_times(
        (Expr*[]){expr_new_integer(-1), fourac}, 2);
    Expr* disc = eval_and_free(internal_plus((Expr*[]){b2, neg_fourac}, 2));

    /* Decide the sign of disc.  Rational cases give a definitive
     * sign; symbolic cases are decided under the positive-symbol
     * assumption (matching IntegrateRational.m's `Refine[#>0]&` filter
     * that gates the parametric path).  When the sign cannot be
     * determined (parametric quadratic with mixed-sign discriminant
     * such as 1/(a x^2 + b x + c)), default to the negative-disc /
     * ArcTan branch with a held Sqrt[-disc].  That formula is the
     * principal antiderivative for both signs of the discriminant
     * (ArcTan continues analytically through imaginary argument),
     * matching Mathematica's parametric-quadratic output. */
    {   /* Over Q(alpha) for a Root object alpha: reduce modulo its
         * minimal polynomial, so the sign test and Sqrt see the small
         * normal form. */
        if (over_root) {
            Expr* dr = intrat_reduce_root_coefficients(disc, x);
            expr_free(disc);
            /* Pull the rational content out, so Sqrt[disc] takes its
             * square part. */
            disc = call_eval("FactorTerms", (Expr*[]){dr}, 1);
        }
    }
    int disc_sign = intrat_constant_sign(disc);

    Expr* two_a = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(2), expr_copy(a)}, 2));

    if (disc_sign > 0) {
        /* Two real roots.  α± = (-b ± Sqrt[disc])/(2 a). */
        Expr* sqrt_d = expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(disc)}, 1);
        sqrt_d = eval_and_free(sqrt_d);
        Expr* neg_b = internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(b)}, 2);
        neg_b = eval_and_free(neg_b);
        Expr* num_plus = eval_and_free(internal_plus(
            (Expr*[]){expr_copy(neg_b), expr_copy(sqrt_d)}, 2));
        Expr* num_minus = eval_and_free(internal_subtract(
            (Expr*[]){neg_b, sqrt_d}, 2));
        Expr* alpha_plus = eval_and_free(internal_divide(
            (Expr*[]){num_plus, expr_copy(two_a)}, 2));
        Expr* alpha_minus = eval_and_free(internal_divide(
            (Expr*[]){num_minus, two_a}, 2));

        /* Each contributes α * Log[s(α, x)]. */
        Expr* s_at_plus  = subst_t(s, t, alpha_plus);
        Expr* s_at_minus = subst_t(s, t, alpha_minus);
        Expr* log_plus = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){s_at_plus}, 1);
        Expr* log_minus = expr_new_function(expr_new_symbol(SYM_Log),
            (Expr*[]){s_at_minus}, 1);
        Expr* term_plus  = internal_times((Expr*[]){alpha_plus, log_plus}, 2);
        Expr* term_minus = internal_times((Expr*[]){alpha_minus, log_minus}, 2);
        expr_free(disc);
        Expr* sum = internal_plus(
            (Expr*[]){term_plus, term_minus}, 2);
        return eval_and_free(sum);
    }

    /* Complex conjugate pair: u = -b/(2a), v = Sqrt[-disc]/(2a). */
    Expr* neg_b2 = eval_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(b)}, 2));
    Expr* u_root = eval_and_free(internal_divide(
        (Expr*[]){neg_b2, expr_copy(two_a)}, 2));
    /* Distribute -1 across `disc` before passing to the Sqrt
     * simplifier.  Without an Expand, Mathilda keeps `-1 * disc` as
     * Times[-1, Plus[...]]; intsimp_pos_sqrt then walks the Times
     * factors and pulls Sqrt[-1] = I out of the `-1`, ending in an
     * imaginary v_root that collapses logtoreal_quadratic to 0. */
    Expr* neg_disc_raw = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), disc}, 2));
    /* Try the positive-symbol Sqrt simplifier first (e.g.
     * Sqrt[4 a^2 b^2] -> 2 a b); fall back to the held Sqrt[..] when
     * it can't extract a clean radical form. */
    Expr* sqrt_neg_disc = intsimp_pos_sqrt(neg_disc_raw);
    /* A sign-undecided monomial such as -16 Sin[k] has a negative
     * numeric factor that intsimp_pos_sqrt turns into Sqrt[-16] = 4 I,
     * making v_root imaginary and collapsing the result to 0.  Keep
     * the radical held instead: v_root must stay formally real. */
    if (intrat_contains_complex(sqrt_neg_disc)) {
        expr_free(sqrt_neg_disc);
        sqrt_neg_disc = expr_new_function(expr_new_symbol(SYM_Sqrt),
            (Expr*[]){expr_copy(neg_disc_raw)}, 1);
    }
    expr_free(neg_disc_raw);
    Expr* v_root = eval_and_free(internal_divide(
        (Expr*[]){sqrt_neg_disc, two_a}, 2));

    /* Substitute t -> u_root + I*v_root in s and split into A + I B. */
    Expr* iv = internal_times((Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Complex),
            (Expr*[]){expr_new_integer(0), expr_new_integer(1)}, 2),
        expr_copy(v_root)
    }, 2);
    iv = eval_and_free(iv);
    Expr* substituted_t = eval_and_free(internal_plus(
        (Expr*[]){expr_copy(u_root), iv}, 2));
    Expr* s_complex = subst_t(s, t, substituted_t);
    expr_free(substituted_t);

    Expr *A = NULL, *B = NULL;
    split_re_im(s_complex, &A, &B);
    expr_free(s_complex);

    /* Term 1: u_root * Log[A^2 + B^2]. */
    Expr* A2 = expand_and_free(internal_power(
        (Expr*[]){expr_copy(A), expr_new_integer(2)}, 2));
    Expr* B2 = expand_and_free(internal_power(
        (Expr*[]){expr_copy(B), expr_new_integer(2)}, 2));
    Expr* mod2 = eval_and_free(internal_plus((Expr*[]){A2, B2}, 2));
    Expr* mod2_can = intrat_canonic(mod2); expr_free(mod2);
    Expr* logmod = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){mod2_can}, 1);
    Expr* term1 = internal_times((Expr*[]){expr_copy(u_root), logmod}, 2);
    term1 = eval_and_free(term1);

    /* Term 2: v_root * LogToAtan[A, B, x]. */
    Expr* atan_part = intrat_log_to_atan(A, B, x);
    if (!atan_part) atan_part = expr_new_integer(0);
    Expr* term2 = internal_times((Expr*[]){expr_copy(v_root), atan_part}, 2);
    term2 = eval_and_free(term2);

    expr_free(A); expr_free(B); expr_free(u_root); expr_free(v_root);

    Expr* sum = internal_plus((Expr*[]){term1, term2}, 2);
    return eval_and_free(sum);
}

/* Handle a sparse nth-root depressed polynomial c_n t^n + c_0 (all
 * intermediate coefficients zero) for n >= 3.  Factor over R via the
 * standard cyclotomic decomposition:
 *
 *   roots are r e^{i θ_k} with r = |c_0 / c_n|^(1/n) and angles
 *     q := -c_0/c_n > 0:  θ_k = 2 k π / n   (k = 0..n-1)
 *     q < 0:              θ_k = (2 k + 1) π / n
 *
 * Real roots correspond to angles 0 and π.  Other roots come in
 * complex-conjugate pairs (k, n−k) for q > 0 or (k, n−1−k) for q < 0.
 * Each pair contributes a quadratic factor
 *     t^2 − 2 r cos(θ_k) t + r^2
 * whose discriminant −4 r^2 sin^2 θ_k is strictly negative, so
 * logtoreal_quadratic always lands in its ArcTan branch.
 *
 * Generalises the deg-3 and deg-4-Sophie-Germain branches that used
 * to live in logtoreal_dispatch.  Catches 1/(b ± a x^n) for all
 * n ≥ 3 once the LRT pipeline reduces to its sparse Q-in-t form.
 *
 * Returns the accumulated real form, or NULL when:
 *   • the polynomial isn't sparse (some intermediate c_k ≠ 0),
 *   • c_0 = 0 (caller should have factored t out upstream), or
 *   • the sign of −c_0/c_n cannot be decided symbolically or numerically.
 */
static Expr* logtoreal_nthroot_sparse(Expr* base, int deg,
                                      Expr* s, Expr* x, Expr* t) {
    if (deg < 3) return NULL;

    for (int k = 1; k < deg; k++) {
        Expr* ck = get_coeff(base, t, k);
        bool zero = (ck->type == EXPR_INTEGER && ck->data.integer == 0);
        expr_free(ck);
        if (!zero) return NULL;
    }

    Expr* c0 = get_coeff(base, t, 0);
    Expr* cn = get_coeff(base, t, deg);
    if (c0->type == EXPR_INTEGER && c0->data.integer == 0) {
        expr_free(c0); expr_free(cn);
        return NULL;
    }

    /* q := −c_0 / c_n. */
    Expr* neg_c0 = expand_and_free(internal_times(
        (Expr*[]){expr_new_integer(-1), expr_copy(c0)}, 2));
    Expr* q = eval_and_free(internal_divide(
        (Expr*[]){neg_c0, expr_copy(cn)}, 2));
    expr_free(c0); expr_free(cn);

    int q_sign = intsimp_sign_pos_assumption(q);
    if (q_sign == 0) q_sign = intsimp_numeric_sign(q);
    if (q_sign == 0) { expr_free(q); return NULL; }

    bool q_pos = (q_sign > 0);

    /* r = |q|^(1/n) using the principal real n-th root. */
    Expr* abs_q = q_pos
        ? q
        : eval_and_free(internal_times(
              (Expr*[]){expr_new_integer(-1), q}, 2));
    Expr* inv_n = expr_new_function(
        expr_new_symbol(SYM_Rational),
        (Expr*[]){expr_new_integer(1), expr_new_integer(deg)}, 2);
    inv_n = eval_and_free(inv_n);
    Expr* r = eval_and_free(internal_power(
        (Expr*[]){abs_q, inv_n}, 2));

    size_t cap = (size_t)deg / 2 + 4;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * cap);
    size_t nterms = 0;
    bool fail = false;

    /* --- Linear factors for real roots. -------------------------- */
    if (q_pos) {
        /* (t − r) — always present when q > 0. */
        Expr* lin = build_log_term(r, s, t);
        if (!lin) { fail = true; goto done; }
        terms[nterms++] = lin;
        if (deg % 2 == 0) {
            /* (t + r) when n is even. */
            Expr* neg_r = eval_and_free(internal_times(
                (Expr*[]){expr_new_integer(-1), expr_copy(r)}, 2));
            Expr* lin2 = build_log_term(neg_r, s, t);
            expr_free(neg_r);
            if (!lin2) { fail = true; goto done; }
            terms[nterms++] = lin2;
        }
    } else {
        /* q < 0: real root only when n is odd, at t = −r. */
        if (deg % 2 == 1) {
            Expr* neg_r = eval_and_free(internal_times(
                (Expr*[]){expr_new_integer(-1), expr_copy(r)}, 2));
            Expr* lin = build_log_term(neg_r, s, t);
            expr_free(neg_r);
            if (!lin) { fail = true; goto done; }
            terms[nterms++] = lin;
        }
    }

    /* --- Quadratic factors for complex-conjugate pairs. ---------- */
    int k_start, k_end;
    if (q_pos) {
        k_start = 1;
        k_end = (deg - 1) / 2;          /* floor((n-1)/2) */
    } else {
        k_start = 0;
        k_end = (deg % 2 == 1)
            ? (deg - 3) / 2              /* odd n: skip k=(n-1)/2 (θ=π real root) */
            : (deg - 2) / 2;             /* even n: floor((n-2)/2) */
    }

    Expr* r_sq = expand_and_free(internal_power(
        (Expr*[]){expr_copy(r), expr_new_integer(2)}, 2));

    for (int k = k_start; k <= k_end; k++) {
        int angle_num_k = q_pos ? (2 * k) : (2 * k + 1);
        Expr* angle_num = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(angle_num_k),
                      expr_new_symbol(SYM_Pi)}, 2));
        Expr* angle = eval_and_free(internal_divide(
            (Expr*[]){angle_num, expr_new_integer(deg)}, 2));
        Expr* cos_v = expr_new_function(expr_new_symbol(SYM_Cos),
            (Expr*[]){angle}, 1);
        cos_v = eval_and_free(cos_v);

        Expr* b_quad = eval_and_free(internal_times(
            (Expr*[]){expr_new_integer(-2),
                      expr_copy(r),
                      cos_v}, 3));

        Expr* one_a = expr_new_integer(1);
        Expr* r_sq_copy = expr_copy(r_sq);
        Expr* quad = logtoreal_quadratic(
            one_a,
            b_quad,
            r_sq_copy,
            s, x, t);
        expr_free(one_a);
        expr_free(b_quad);
        expr_free(r_sq_copy);
        if (!quad) { fail = true; break; }
        terms[nterms++] = quad;
    }

    expr_free(r_sq);
    expr_free(r);

done:
    if (fail || nterms == 0) {
        for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
        free(terms);
        return NULL;
    }
    if (nterms == 1) {
        Expr* result = terms[0];
        free(terms);
        return result;
    }
    Expr* sum = internal_plus(terms, nterms);
    free(terms);
    return eval_and_free(sum);
}

/* Walk a Factor[Q]'s output and dispatch every factor.  Returns the
 * accumulated real form, or NULL if any factor exceeds the bounded-
 * Solve scope (degree > 2 in t). */
static Expr* logtoreal_dispatch(Expr* factored, Expr* s, Expr* x, Expr* t) {
    Expr** factors;
    size_t nfactors;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol.name == SYM_Times) {
        factors = factored->data.function.args;
        nfactors = factored->data.function.arg_count;
    } else {
        factors = &factored;
        nfactors = 1;
    }
    Expr** acc = (Expr**)malloc(sizeof(Expr*) * (nfactors * 2));
    size_t acc_n = 0;

    for (size_t i = 0; i < nfactors; i++) {
        Expr* fac = factors[i];
        Expr* base = fac;
        int multiplicity = 1;
        if (fac->type == EXPR_FUNCTION
            && fac->data.function.head->type == EXPR_SYMBOL
            && fac->data.function.head->data.symbol.name == SYM_Power
            && fac->data.function.arg_count == 2
            && fac->data.function.args[1]->type == EXPR_INTEGER
            && fac->data.function.args[1]->data.integer >= 1) {
            base = fac->data.function.args[0];
            multiplicity = (int)fac->data.function.args[1]->data.integer;
        }
        int deg = get_degree_poly(base, t);
        if (deg == 0) continue;  /* numeric coefficient */

        Expr* contribution = NULL;
        if (deg == 1) {
            Expr* alpha = extract_linear_root(base, t);
            if (!alpha) goto fail;
            contribution = build_log_term(alpha, s, t);
            expr_free(alpha);
        } else if (deg == 2) {
            Expr* a = get_coeff(base, t, 2);
            Expr* b = get_coeff(base, t, 1);
            Expr* c = get_coeff(base, t, 0);
            contribution = logtoreal_quadratic(a, b, c, s, x, t);
            expr_free(a); expr_free(b); expr_free(c);
            if (!contribution) goto fail;
        } else if (deg >= 3) {
            /* General sparse nth-root path: c_n t^n + c_0 (all
             * intermediate coefficients zero), dispatched via the
             * cyclotomic factorisation of t^n − q.  Subsumes the
             * previous deg-3 nth-root branch and the deg-4
             * Sophie-Germain shortcut; catches 1/(b ± a x^n) for
             * every n ≥ 3.  When sparse fails for deg == 4 we fall
             * through to the deg-4 biquadratic / palindromic
             * specialisations below. */
            contribution = logtoreal_nthroot_sparse(base, deg, s, x, t);
            if (!contribution && deg != 4) goto fail;
        }
        if (!contribution && deg == 4) {
            /* Non-sparse quartic: try biquadratic (c1 = c3 = 0,
             * c2 ≠ 0) or palindromic specialisations.  The fully
             * sparse case (c1 = c2 = c3 = 0) was already handled by
             * logtoreal_nthroot_sparse above; if it didn't fire, the
             * sign of c0/c4 was undecidable and there's nothing more
             * to do at this layer. */
            Expr* c0 = get_coeff(base, t, 0);
            Expr* c1 = get_coeff(base, t, 1);
            Expr* c2 = get_coeff(base, t, 2);
            Expr* c3 = get_coeff(base, t, 3);
            Expr* c4 = get_coeff(base, t, 4);

            bool c1z = (c1->type == EXPR_INTEGER && c1->data.integer == 0);
            bool c2z = (c2->type == EXPR_INTEGER && c2->data.integer == 0);
            bool c3z = (c3->type == EXPR_INTEGER && c3->data.integer == 0);

            if (c1z && c2z && c3z) {
                /* Fully sparse but sparse handler bailed (sign of
                 * c0/c4 undecidable both symbolically and numerically).
                 * No further branch can do better. */
                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);
                goto fail;
            } else if (c1z && c3z) {
                /* Genuine biquadratic in t: c4 t^4 + c2 t^2 + c0
                 * (c2 != 0).  Substitute u = t^2 to land in
                 *   c4 u^2 + c2 u + c0 = 0
                 * with roots u_± = (-c2 ± Sqrt[c2^2 - 4 c0 c4]) / (2 c4).
                 * Q(t) factors as c4 (t^2 - u_+)(t^2 - u_-).
                 *
                 * Soundness gate: only fire when c2^2 - 4 c0 c4 is
                 * provably positive under the positive-symbol
                 * assumption.  Otherwise u_± would be a complex
                 * conjugate pair, the (t^2 - u_±) factors would have
                 * complex coefficients, and logtoreal_quadratic
                 * (which expects a real-coefficient quadratic in t)
                 * would silently miscompute -- typically collapsing
                 * to 0 because the imaginary parts cancel between
                 * the two sub-factor calls.  Catches the parametric
                 * family 1/(a x^4 + b x^2 + c) under the
                 * b^2 > 4 a c assumption (which is always true in
                 * the positive-symbol model) while leaving cases
                 * with negative or sign-unknown inner discriminant
                 * to fall through to NaiveLogPart's held RootSum
                 * (where they were before this branch existed). */
                Expr* c2sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(c2), expr_new_integer(2)}, 2));
                Expr* fourc0c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(4), expr_copy(c0), expr_copy(c4)}, 3));
                Expr* neg_fourc0c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), fourc0c4}, 2));
                Expr* discsq = eval_and_free(internal_plus(
                    (Expr*[]){c2sq, neg_fourc0c4}, 2));
                int inner_disc_sign = intsimp_sign_pos_assumption(discsq);
                if (inner_disc_sign <= 0) {
                    /* Inner disc not provably positive.  Before giving
                     * up, try the Sophie-Germain-with-c2 factorisation
                     * over R:
                     *   c4 t^4 + c2 t^2 + c0
                     *      = c4 (t^2 + α t + β)(t^2 − α t + β)
                     * with β = Sqrt[c0/c4], α^2 = 2 β − c2/c4.
                     * Real when (a) c0/c4 > 0 (so β is real) and
                     * (b) 2 β − c2/c4 > 0 (so α is real).  Both gates
                     * fire under the positive-symbol assumption, and
                     * algebraically: c2^2 − 4 c0 c4 < 0 ⇒
                     * |c2| < 2 Sqrt[c0 c4] ⇒ 2 Sqrt[c0/c4] − c2/c4 =
                     * (2 Sqrt[c0 c4] − c2) / c4 > 0 when c4 > 0,
                     * so the gate is automatic in the negative-inner-
                     * disc regime we land in here.
                     *
                     * Each resulting real quadratic has its own
                     * discriminant α^2 − 4 β < 0, which routes through
                     * logtoreal_quadratic's ArcTan branch.  This is
                     * the path Mathematica's BronsteinRational takes
                     * for inputs like (-1+x^2)/(1-2x^2+2x^4) whose
                     * Q-in-t is 2 t^4 − 2 t^2 + 1 (inner disc −4). */
                    expr_free(discsq);

                    Expr* ratio = eval_and_free(internal_divide(
                        (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2));
                    bool ratio_pos = (intsimp_sign_pos_assumption(ratio) > 0);

                    if (!ratio_pos) {
                        expr_free(ratio);
                        expr_free(c0); expr_free(c1); expr_free(c2);
                        expr_free(c3); expr_free(c4);
                        goto fail;
                    }

                    Expr* beta = intsimp_pos_sqrt(ratio); /* Sqrt[c0/c4] */
                    expr_free(ratio);

                    Expr* c2_over_c4 = eval_and_free(internal_divide(
                        (Expr*[]){expr_copy(c2), expr_copy(c4)}, 2));
                    Expr* neg_c2_over_c4 = expand_and_free(internal_times(
                        (Expr*[]){expr_new_integer(-1), c2_over_c4}, 2));

                    /* α^2 = 2 β − c2/c4. */
                    Expr* two_beta = eval_and_free(internal_times(
                        (Expr*[]){expr_new_integer(2), expr_copy(beta)}, 2));
                    Expr* alpha_sq = eval_and_free(internal_plus(
                        (Expr*[]){two_beta, neg_c2_over_c4}, 2));
                    Expr* alpha_sq_expanded = expr_expand(alpha_sq);
                    expr_free(alpha_sq);

                    if (intsimp_sign_pos_assumption(alpha_sq_expanded) <= 0) {
                        /* α^2 not provably positive — bail. */
                        expr_free(alpha_sq_expanded);
                        expr_free(beta);
                        expr_free(c0); expr_free(c1); expr_free(c2);
                        expr_free(c3); expr_free(c4);
                        goto fail;
                    }

                    Expr* alpha = intsimp_pos_sqrt(alpha_sq_expanded);
                    expr_free(alpha_sq_expanded);
                    Expr* neg_alpha = eval_and_free(internal_times(
                        (Expr*[]){expr_new_integer(-1),
                                  expr_copy(alpha)}, 2));

                    Expr* one_p1 = expr_new_integer(1);
                    Expr* beta_copy = expr_copy(beta);
                    Expr* part1 = logtoreal_quadratic(
                        one_p1,
                        alpha,
                        beta_copy,
                        s, x, t);
                    expr_free(one_p1); expr_free(beta_copy);

                    Expr* one_p2 = expr_new_integer(1);
                    Expr* part2 = logtoreal_quadratic(
                        one_p2,
                        neg_alpha,
                        beta,
                        s, x, t);
                    expr_free(one_p2);
                    expr_free(alpha); expr_free(neg_alpha); expr_free(beta);

                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);

                    if (part1 && part2) {
                        contribution = eval_and_free(internal_plus(
                            (Expr*[]){part1, part2}, 2));
                        for (int m = 0; m < multiplicity; m++) {
                            acc[acc_n++] = (m == 0) ? contribution : expr_copy(contribution);
                        }
                        continue; /* next factor */
                    }
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
                Expr* sqrt_disc = intsimp_pos_sqrt(discsq);
                expr_free(discsq);

                Expr* neg_c2 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), expr_copy(c2)}, 2));
                Expr* two_c4 = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));

                Expr* num_plus = eval_and_free(internal_plus(
                    (Expr*[]){expr_copy(neg_c2), expr_copy(sqrt_disc)}, 2));
                Expr* num_minus = eval_and_free(internal_subtract(
                    (Expr*[]){neg_c2, sqrt_disc}, 2));
                Expr* u_plus = eval_and_free(internal_divide(
                    (Expr*[]){num_plus, expr_copy(two_c4)}, 2));
                Expr* u_minus = eval_and_free(internal_divide(
                    (Expr*[]){num_minus, two_c4}, 2));

                Expr* neg_u_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), u_plus}, 2));
                Expr* neg_u_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), u_minus}, 2));

                Expr* one_a1 = expr_new_integer(1);
                Expr* zero_b1 = expr_new_integer(0);
                Expr* part1 = logtoreal_quadratic(
                    one_a1,
                    zero_b1,
                    neg_u_plus,
                    s, x, t);
                expr_free(one_a1); expr_free(zero_b1); expr_free(neg_u_plus);

                Expr* one_a2 = expr_new_integer(1);
                Expr* zero_b2 = expr_new_integer(0);
                Expr* part2 = logtoreal_quadratic(
                    one_a2,
                    zero_b2,
                    neg_u_minus,
                    s, x, t);
                expr_free(one_a2); expr_free(zero_b2); expr_free(neg_u_minus);

                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);

                if (part1 && part2) {
                    contribution = eval_and_free(internal_plus(
                        (Expr*[]){part1, part2}, 2));
                } else {
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
            } else {
                /* Scaled-palindromic quartic:
                 *   Q(t) = c4 t^4 + c3 t^3 + c2 t^2 + c1 t + c0
                 * with c4 c1^2 = c3^2 c0 (the symmetry condition that
                 * generalises the c0=c4, c1=c3 palindromic family).
                 * Under u = r t with r^2 = c3/c1 = Sqrt[c4/c0] the
                 * polynomial r^n Q(u/r) is palindromic in u:
                 *   c4 u^4 + (c3 r) u^3 + (c2 r^2) u^2 + (c1 r^3) u
                 *   + (c0 r^4)   with c0 r^4 = c4 and c1 r^3 = c3 r.
                 * Substitute v = u + 1/u (so u^2 + 1/u^2 = v^2 − 2):
                 *   c4 v^2 + (c3 r) v + (c2 r^2 − 2 c4) = 0.
                 * Solve for v_±; each (u^2 − v_± u + 1) translates to
                 *   r^2 t^2 − r v_± t + 1
                 * which is dispatched through logtoreal_quadratic.
                 *
                 * Catches Q-in-t for 1/(x^5+1) (Q = 625 t^4 + 125 t^3
                 * + 25 t^2 + 5 t + 1 — palindromic in u = 5t) and the
                 * scaled cyclotomic-quartic family.  When the symmetry
                 * condition fails or any sign-gate is undecidable we
                 * bail to RootSum. */

                /* Symmetry check: c4 c1^2 == c3^2 c0  AND  c0 != 0.
                 * General scaled-palindromic detection.  We further
                 * restrict to r = 1 (pure palindromic) below, because
                 * the scaled cases produce nested-radical S(t/r, x)
                 * polynomials that overwhelm the downstream LogToAtan
                 * GCD machinery.  Scaled-palindromic-with-clean-r is a
                 * future extension. */
                Expr* sym_lhs = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(c4),
                              expand_and_free(internal_power(
                                  (Expr*[]){expr_copy(c1),
                                            expr_new_integer(2)}, 2))}, 2));
                Expr* sym_rhs = expand_and_free(internal_times(
                    (Expr*[]){expand_and_free(internal_power(
                                  (Expr*[]){expr_copy(c3),
                                            expr_new_integer(2)}, 2)),
                              expr_copy(c0)}, 2));
                Expr* sym_diff = internal_subtract(
                    (Expr*[]){sym_lhs, sym_rhs}, 2);
                bool symmetric = intsimp_zero_q(sym_diff);
                expr_free(sym_diff);

                bool c0_nonzero = !(c0->type == EXPR_INTEGER && c0->data.integer == 0);
                bool c1_nonzero = !(c1->type == EXPR_INTEGER && c1->data.integer == 0);

                if (!symmetric || !c0_nonzero || !c1_nonzero) {
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }

                /* r^2 = c3 / c1.  Restrict to r = 1 (c1 == c3 AND
                 * c0 == c4) so the resulting LogToReal substitution
                 * stays radical-clean. */
                Expr* c1mc3 = internal_subtract(
                    (Expr*[]){expr_copy(c1), expr_copy(c3)}, 2);
                Expr* c0mc4 = internal_subtract(
                    (Expr*[]){expr_copy(c0), expr_copy(c4)}, 2);
                bool r_is_one = intsimp_zero_q(c1mc3) && intsimp_zero_q(c0mc4);
                expr_free(c1mc3); expr_free(c0mc4);
                if (!r_is_one) {
                    /* Scaled-palindromic with r != 1: bail (the
                     * substituted S(t/r, x) blows up downstream). */
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }
                Expr* r_sq = expr_new_integer(1);
                Expr* r    = expr_new_integer(1);

                /* v-quadratic: c4 v^2 + (c3 r) v + (c2 r^2 − 2 c4) = 0.
                 * disc_v = (c3 r)^2 − 4 c4 (c2 r^2 − 2 c4)
                 *        = c3^2 r^2 − 4 c2 c4 r^2 + 8 c4^2.
                 * Soundness gate: disc_v >= 0 so v_± are real. */
                Expr* c3r = eval_and_free(internal_times(
                    (Expr*[]){expr_copy(c3), expr_copy(r)}, 2));
                Expr* c3r_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(c3r), expr_new_integer(2)}, 2));
                Expr* two_c4 = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));
                Expr* c2_rsq = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(c2), expr_copy(r_sq)}, 2));
                Expr* v_const = eval_and_free(internal_subtract(
                    (Expr*[]){c2_rsq, two_c4}, 2));
                Expr* four_c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(4), expr_copy(c4),
                              expr_copy(v_const)}, 3));
                Expr* neg_four_c4 = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), four_c4}, 2));
                Expr* disc_v = eval_and_free(internal_plus(
                    (Expr*[]){c3r_sq, neg_four_c4}, 2));
                Expr* disc_v_exp = expr_expand(disc_v); expr_free(disc_v);
                int disc_v_sign = intsimp_sign_pos_assumption(disc_v_exp);
                if (disc_v_sign == 0) disc_v_sign = intsimp_numeric_sign(disc_v_exp);

                if (disc_v_sign < 0) {
                    /* v complex — bail. */
                    expr_free(v_const); expr_free(disc_v_exp);
                    expr_free(c3r); expr_free(r); expr_free(r_sq);
                    expr_free(c0); expr_free(c1); expr_free(c2);
                    expr_free(c3); expr_free(c4);
                    goto fail;
                }
                Expr* sqrt_disc_v = intsimp_pos_sqrt(disc_v_exp);
                expr_free(disc_v_exp);

                /* v_± = (−c3 r ± Sqrt[disc_v]) / (2 c4). */
                Expr* neg_c3r = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), expr_copy(c3r)}, 2));
                Expr* v_num_p = eval_and_free(internal_plus(
                    (Expr*[]){expr_copy(neg_c3r), expr_copy(sqrt_disc_v)}, 2));
                Expr* v_num_m = eval_and_free(internal_subtract(
                    (Expr*[]){neg_c3r, sqrt_disc_v}, 2));
                Expr* two_c4_b = eval_and_free(internal_times(
                    (Expr*[]){expr_new_integer(2), expr_copy(c4)}, 2));
                Expr* v_plus = eval_and_free(internal_divide(
                    (Expr*[]){v_num_p, expr_copy(two_c4_b)}, 2));
                Expr* v_minus = eval_and_free(internal_divide(
                    (Expr*[]){v_num_m, two_c4_b}, 2));

                expr_free(v_const); expr_free(c3r);
                expr_free(c0); expr_free(c1); expr_free(c2);
                expr_free(c3); expr_free(c4);

                /* Each v_± yields the quadratic in t:
                 *   r^2 t^2 − (r v_±) t + 1
                 * Discriminant Δ_t = r^2 v^2 − 4 r^2 = r^2 (v^2 − 4).
                 * Verify Δ_t < 0 so logtoreal_quadratic's ArcTan
                 * branch fires.  When the symbolic walker is
                 * inconclusive we fall back to numerical N[]. */
                Expr* v_plus_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(v_plus), expr_new_integer(2)}, 2));
                Expr* v_plus_sq_m4 = eval_and_free(internal_plus(
                    (Expr*[]){v_plus_sq, expr_new_integer(-4)}, 2));
                Expr* v_plus_check_arg = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r_sq), expr_copy(v_plus_sq_m4)}, 2));
                expr_free(v_plus_sq_m4);
                int v_plus_check = intsimp_sign_pos_assumption(v_plus_check_arg);
                if (v_plus_check == 0) v_plus_check = intsimp_numeric_sign(v_plus_check_arg);
                expr_free(v_plus_check_arg);

                Expr* v_minus_sq = expand_and_free(internal_power(
                    (Expr*[]){expr_copy(v_minus), expr_new_integer(2)}, 2));
                Expr* v_minus_sq_m4 = eval_and_free(internal_plus(
                    (Expr*[]){v_minus_sq, expr_new_integer(-4)}, 2));
                Expr* v_minus_check_arg = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r_sq), expr_copy(v_minus_sq_m4)}, 2));
                expr_free(v_minus_sq_m4);
                int v_minus_check = intsimp_sign_pos_assumption(v_minus_check_arg);
                if (v_minus_check == 0) v_minus_check = intsimp_numeric_sign(v_minus_check_arg);
                expr_free(v_minus_check_arg);

                if (v_plus_check >= 0 || v_minus_check >= 0) {
                    /* One of the quadratic factors would have a non-
                     * negative discriminant.  Bail. */
                    expr_free(v_plus); expr_free(v_minus);
                    expr_free(r); expr_free(r_sq);
                    goto fail;
                }

                /* Build (r^2) t^2 + (−r v) t + 1 for each v_±. */
                Expr* rv_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r), expr_copy(v_plus)}, 2));
                Expr* neg_rv_plus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), rv_plus}, 2));
                Expr* rv_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_copy(r), expr_copy(v_minus)}, 2));
                Expr* neg_rv_minus = expand_and_free(internal_times(
                    (Expr*[]){expr_new_integer(-1), rv_minus}, 2));

                Expr* r_sq_c1 = expr_copy(r_sq);
                Expr* one_c1 = expr_new_integer(1);
                Expr* part1 = logtoreal_quadratic(
                    r_sq_c1,
                    neg_rv_plus,
                    one_c1,
                    s, x, t);
                expr_free(r_sq_c1); expr_free(neg_rv_plus); expr_free(one_c1);

                Expr* r_sq_c2 = expr_copy(r_sq);
                Expr* one_c2 = expr_new_integer(1);
                Expr* part2 = logtoreal_quadratic(
                    r_sq_c2,
                    neg_rv_minus,
                    one_c2,
                    s, x, t);
                expr_free(r_sq_c2); expr_free(neg_rv_minus); expr_free(one_c2);

                expr_free(v_plus); expr_free(v_minus);
                expr_free(r); expr_free(r_sq);

                if (part1 && part2) {
                    contribution = eval_and_free(internal_plus(
                        (Expr*[]){part1, part2}, 2));
                } else {
                    if (part1) expr_free(part1);
                    if (part2) expr_free(part2);
                    goto fail;
                }
            }
        }
        if (!contribution) goto fail;
        for (int m = 0; m < multiplicity; m++) {
            acc[acc_n++] = (m == 0) ? contribution : expr_copy(contribution);
        }
    }
    if (acc_n == 0) { free(acc); return expr_new_integer(0); }
    if (acc_n == 1) { Expr* r = acc[0]; free(acc); return r; }
    Expr* sum = internal_plus(acc, acc_n);
    free(acc);
    return eval_and_free(sum);

fail:
    for (size_t k = 0; k < acc_n; k++) expr_free(acc[k]);
    free(acc);
    return NULL;
}

/* LogToReal[r, s, x, t]: factor r over Q, then dispatch each linear /
 * quadratic factor through logtoreal_quadratic.  Returns NULL when r
 * has any irreducible factor of degree >= 3 over Q — the caller is
 * expected to fall back to NaiveLogPart in that case.  We deliberately
 * do NOT guess algebraic extensions here; algebraic-extension closure
 * is the responsibility of Solve (in solve.c, forthcoming). */
static Expr* intrat_log_to_real(Expr* r, Expr* s, Expr* x, Expr* t) {
    intrat_trace("LogToReal", "IN", r);


    Expr* factored = internal_factor((Expr*[]){expr_copy(r)}, 1);
    Expr* factored_e = evaluate(factored);
    expr_free(factored);
    Expr* result = logtoreal_dispatch(factored_e, s, x, t);
    expr_free(factored_e);
    if (result) { intrat_trace("LogToReal", "OUT", result); return result; }

    intrat_trace("LogToReal", "OUT (failed)", r);
    return NULL;
}

/* ----- Linear-Q closer (Phase 2's contribution to closed forms) -----
 *
 * Try to evaluate the {Q_i, S_i} pair list to a closed-form Log sum
 * by factoring each Q_i over Q[t] (or Q[t, alpha]) and contributing
 * α Log[S_i(α, x)] for every linear factor (t - α).
 *
 * Phase 2 closes the case where every linear factor lies over Q —
 * which covers `1/((x-a)(x-b))`, `1/(x^n - c^n)` for rational c, and
 * the rational-roots fragment of any LRT problem.  Quadratic factors
 * with negative or symbolic discriminants stay unhandled here; Phase
 * 4's LogToReal will pick them up.
 */

/* Walk a Factor[poly]'s Times-headed result and append every linear
 * (t - α) factor's contribution to *terms.  Returns false if any
 * factor is non-linear in t with non-trivial t-degree (caller falls
 * back to the symbolic form). */
static bool collect_linear_factors(Expr* factored, Expr* S, Expr* t,
                                    Expr*** terms, size_t* nterms,
                                    size_t* tcap) {
    Expr** factors;
    size_t nfactors;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol.name == SYM_Times) {
        factors = factored->data.function.args;
        nfactors = factored->data.function.arg_count;
    } else {
        factors = &factored;
        nfactors = 1;
    }

    for (size_t i = 0; i < nfactors; i++) {
        Expr* fac = factors[i];

        /* (t - α)^n: bind α once, contribute n times.  Power[base, k]
         * with base linear in t. */
        Expr* base = fac;
        int multiplicity = 1;
        if (fac->type == EXPR_FUNCTION
            && fac->data.function.head->type == EXPR_SYMBOL
            && fac->data.function.head->data.symbol.name == SYM_Power
            && fac->data.function.arg_count == 2
            && fac->data.function.args[1]->type == EXPR_INTEGER
            && fac->data.function.args[1]->data.integer >= 1) {
            base = fac->data.function.args[0];
            multiplicity = (int)fac->data.function.args[1]->data.integer;
        }

        int deg = get_degree_poly(base, t);
        if (deg == 0) continue;  /* constant factor */
        if (deg != 1) return false;

        Expr* alpha = extract_linear_root(base, t);
        if (!alpha) return false;
        Expr* contribution = build_log_term(alpha, S, t);
        expr_free(alpha);
        for (int m = 0; m < multiplicity; m++) {
            if (*nterms >= *tcap) {
                *tcap *= 2;
                *terms = (Expr**)realloc(*terms, sizeof(Expr*) * (*tcap));
            }
            (*terms)[(*nterms)++] = (m == 0) ? contribution : expr_copy(contribution);
        }
    }
    return true;
}

/* Try to evaluate the {Q, S} pair list to a closed-form Log sum
 * by factoring each Q_i over Q.  Returns NULL when any Q_i has
 * factors of degree > 1 (Phase 4 territory). */
static Expr* intrat_linear_q_closer(Expr* pair_list, Expr* x, Expr* t) {
    (void)x;
    if (!pair_list || pair_list->type != EXPR_FUNCTION
        || pair_list->data.function.head->type != EXPR_SYMBOL
        || pair_list->data.function.head->data.symbol.name != SYM_List) return NULL;
    size_t n = pair_list->data.function.arg_count;
    /* Empty pair list = LRT producer couldn't decompose.  See the
     * matching guard in intrat_log_to_real_pairs above for why returning
     * Integer[0] would be a silent wrong-answer. */
    if (n == 0) return NULL;

    size_t tcap = 8, nterms = 0;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * tcap);

    for (size_t i = 0; i < n; i++) {
        Expr* pair = pair_list->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) {
            for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
            free(terms);
            if (Qi) expr_free(Qi);
            if (Si) expr_free(Si);
            return NULL;
        }

        Expr* factored = internal_factor((Expr*[]){expr_copy(Qi)}, 1);
        Expr* factored_e = evaluate(factored);
        expr_free(factored);

        bool ok = collect_linear_factors(factored_e, Si, t, &terms, &nterms, &tcap);
        expr_free(factored_e);
        expr_free(Qi); expr_free(Si);
        if (!ok) {
            for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
            free(terms);
            return NULL;
        }
    }

    Expr* sum;
    if (nterms == 0) { sum = expr_new_integer(0); free(terms); }
    else if (nterms == 1) { sum = terms[0]; free(terms); }
    else { sum = internal_plus(terms, nterms); free(terms); sum = eval_and_free(sum); }
    return sum;
}

/* ====================================================================
 * Phase 6 + 7 — output simplification of the resulting integral
 * ====================================================================
 *
 * Log pairing (c Log[A] ± c Log[B] -> Log[A·B] / Log[A/B] / ArcTanh),
 * Log[c·p] -> Log[p] when c is free of x, Plus distribution, and
 * ArcTan / ArcTanh sign normalisation all live in intsimp.c — see
 * intsimp.h for the public surface.
 */


/* ====================================================================
 * Phase 5 — IntegrateRealRationalFunction-style per-summand loop.
 * ====================================================================
 *
 * Apart-split a squarefree-denominator residual h, ExtractConstants
 * each summand to keep coefficients out of the subresultant chain,
 * and process each piece independently.  If *every* piece closes via
 * the derivative-recognition / LRT / LogToReal pipeline we return
 * the sum; otherwise NULL so the caller can attempt the whole-h LRT
 * path (or surface the integral as unevaluated).
 *
 * Direct port of IntegrateRational.m:534-587, minus the `NaiveLogPart`
 * fallback (we prefer to surface the integral as unevaluated rather
 * than emit a symbolic RootSum that Mathilda can't yet differentiate
 * back to the integrand).
 */

/* LRT consumer: apply Rioboo's intrat_log_to_real to every {Q_i, S_i}
 * pair and sum the results.  Returns the closed real-elementary form
 * on success, or NULL if any pair's Q_i has an irreducible factor of
 * degree >= 5 (intrat_log_to_real's failure mode) or any list entry
 * is malformed.  Caller retains ownership of `pair_list`. */
static Expr* intrat_log_to_real_pairs(Expr* pair_list, Expr* x, Expr* t) {
    if (!pair_list || pair_list->type != EXPR_FUNCTION
        || pair_list->data.function.head->type != EXPR_SYMBOL
        || pair_list->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }
    size_t npairs = pair_list->data.function.arg_count;
    /* Empty pair list means the LRT producer couldn't decompose the
     * input -- typically the denominator was irreducible over the
     * working ring (e.g. x^16 - x^8 - 1, an irreducible biquadratic-
     * in-x^8 with sqrt-5 roots).  Returning Integer[0] here would
     * make the caller treat the integral as identically zero --
     * silent wrong-answer.  Return NULL so try_lrt_close falls
     * through to NaiveLogPart's universal RootSum form. */
    if (npairs == 0) return NULL;
    Expr** lr_terms = (Expr**)malloc(sizeof(Expr*) * npairs);
    size_t lr_count = 0;

    for (size_t i = 0; i < npairs; i++) {
        Expr* pair = pair_list->data.function.args[i];
        Expr* Qi = list_get(pair, 1);
        Expr* Si = list_get(pair, 2);
        if (!Qi || !Si) {
            if (Qi) expr_free(Qi);
            if (Si) expr_free(Si);
            for (size_t k = 0; k < lr_count; k++) expr_free(lr_terms[k]);
            free(lr_terms);
            return NULL;
        }
        Expr* term = intrat_log_to_real(Qi, Si, x, t);
        expr_free(Qi); expr_free(Si);
        if (!term) {
            for (size_t k = 0; k < lr_count; k++) expr_free(lr_terms[k]);
            free(lr_terms);
            return NULL;
        }
        lr_terms[lr_count++] = term;
    }

    Expr* result;
    if (lr_count == 0) {
        result = expr_new_integer(0);
    } else if (lr_count == 1) {
        result = lr_terms[0];
        lr_terms[0] = NULL;  /* prevent double-free */
    } else {
        Expr* sum = internal_plus(lr_terms, lr_count);
        memset(lr_terms, 0, sizeof(Expr*) * lr_count);
        result = eval_and_free(sum);
    }
    free(lr_terms);
    return result;
}

/* Try to close `num/den` via the LRT-then-real-form pipeline.
 *
 * Architecture (mirrors AXIOM's INTRAT/IRRF separation in
 * intrf.spad.pamphlet:776-778, plus irexpand):
 *   - Producer:  intrat_int_rational_log_part -- single LRT call per
 *                try_lrt_close, emits the {Q_i, S_i} pair list.
 *   - LRT consumers (uniform `(pairs, x, t) -> Expr*|NULL` shape):
 *       * intrat_log_to_real_pairs  -- Rioboo real-elementary form
 *       * intrat_linear_q_closer    -- every Q_i splits over Q
 *   - Non-LRT fallback:  intrat_naive_log_part(num/den) -- Lagrange
 *     RootSum form, universal so try_lrt_close never returns NULL
 *     for a well-formed proper rational integrand.
 *
 * Closure preference matches IntegrateRational.m:560-572. */
/* Function[z, body] -> Function[body /. z -> #1] throughout e, so a RootSum
 * built over the residue variable prints as the reference's `... &` form.
 * Borrows e and z; owned result. */
static Expr* intrat_slot_functions(Expr* e, Expr* z) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Function
        && e->data.function.arg_count == 2 && expr_eq(e->data.function.args[0], z)) {
        Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot), (Expr*[]){ expr_new_integer(1) }, 1);
        Expr* body = subst_t(e->data.function.args[1], z, slot);
        expr_free(slot);
        return expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ body }, 1);
    }
    size_t n = e->data.function.arg_count;
    Expr** a = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = intrat_slot_functions(e->data.function.args[i], z);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return eval_and_free(r);
}

static Expr* try_lrt_close(Expr* num, Expr* den, Expr* x) {
    if (is_zero_poly(num)) return expr_new_integer(0);

    Expr* h = internal_divide(
        (Expr*[]){expr_copy(num), expr_copy(den)}, 2);
    Expr* h_eval = eval_and_free(h);
    Expr* t_sym = expr_new_symbol("Integrate`Private`tt$");

    /* The log part: Rothstein-Trager or Czichowski natively over Q when
     * selected (LogPart -> ...), else -- and on any decline -- the
     * Lazard-Rioboo-Trager producer, which also serves parameters,
     * radicals and Root coefficients. */
    IntratMethodOptions iro = ir_opts();
    int lp = iro.logpart >= 0 ? iro.logpart : INTRAT_DEFAULT_LOGPART;
    int rf = iro.realform;
    Expr* pairs = lp != INTRAT_LOGPART_LRT ? intrat_alg_logpart(num, den, x, t_sym, lp) : NULL;
    /* Automatic: Rothstein-Trager declines when an irreducible residue factor
     * has degree > 8 (its gcd over Q(a) then costs more); there Czichowski's
     * lex basis is the fastest of the three on benchmark 266 (rootsum, stress:
     * 4-5x faster than RT or LRT), and LRT remains the general fallback. */
    if (!pairs && iro.logpart < 0)
        pairs = intrat_alg_logpart(num, den, x, t_sym, INTRAT_LOGPART_CZICHOWSKI);
    if (!pairs) pairs = intrat_int_rational_log_part(h_eval, x, t_sym, false);
    if (pairs) {
        /* Coefficients in Q(alpha) for a Root object alpha: reduce the
         * resultant factors Q_i(t) and S_i(t, x) modulo alpha's minimal
         * polynomial before LogToReal does arithmetic on them. */
        Expr* pr = intrat_reduce_root_coefficients(pairs, x);
        expr_free(pairs);
        pairs = pr;
    }

    Expr* result = NULL;
    if (pairs) {
        /* LRT consumer 1: Rioboo LogToReal on every pair.  It turns the log
         * sum over a conjugate-closed set of residues into real ArcTan/Log
         * terms, which is valid only over a REAL coefficient field: the
         * conjugate of a residue must be a residue and S(conj a, x) =
         * conj S(a, x).  With I among the integrand's coefficients its
         * conjugate-pair branch cancelled against itself --
         * Integrate[1/(x^2 + I), x] was 0, and 1/(x^2 + I x + 1), 1/(x^3 + I)
         * were nonzero but wrong (53 of 158 random Gaussian-coefficient
         * integrands).  Such an integrand takes the plain complex log sum
         * (see the !real_field branch below).  The
         * transcendental callers (t = E^(I x) for trig integrands) keep
         * LogToReal: their I is part of the real reconstruction. */
        /* RealForm -> None: the formal sum, residues left algebraic. */
        if (rf == INTRAT_REALFORM_NONE) {
            result = intrat_root_sum_pairs(pairs, t_sym);
            if (result) {
                Expr* sl = intrat_slot_functions(result, t_sym);
                expr_free(result);
                result = sl;
            }
        }
        bool real_field = !intrat_contains_complex(h_eval);
        /* A non-real coefficient field: no real form exists (the residues are
         * not closed under conjugation), whatever RealForm asks for.  The
         * answer is the complex log sum over the residues of the pair list:
         * first the linear closer (every Q_i split into rational or Gaussian
         * linear factors, one Log each -- the form the Risch callers read
         * back), else the sum written out where every Q_i has degree <= 2
         * (the roots of the quadratic formula summed in pairs: one value
         * under either branch of the Sqrt).  NaiveLogPart over the ROOTS OF
         * THE DENOMINATOR, the route these integrands took before, is the
         * same function at up to 8x the size (u/(u^4 + (1 + I) u^2 + 2): 379
         * leaves against 46), and its bulkier residues sent the
         * algebraic-constant pass into casus-irreducibilis radicals (1515
         * leaves in 1.2 s, against 358 in 0.2 s, for the rational part of
         * ((4 + 2 Sqrt[2]) u^2)/((3 + 2 Sqrt[2]) (I + 5 u^2 + 5 I u^4 + u^6))).
         * Comparing the two before the output passes does not work: the
         * denominator-root form is sometimes smaller raw, and that repro then
         * came out at 668 leaves where the residue form compacts to 358.  A
         * Q_i of degree >= 3 keeps NaiveLogPart, which may still close it
         * (1/(u^4 + I) in symmetric radicals). */
        if (!result && !real_field && rf != INTRAT_REALFORM_NONE)
            result = intrat_linear_q_closer(pairs, x, t_sym);
        if (!result && !real_field) {
            Expr* cs = intrat_root_sum_pairs(pairs, t_sym);
            if (cs && !intrat_contains_rootsum(cs)) result = cs;
            else if (cs) expr_free(cs);
        }
        /* RealForm -> "LogToReal": Bronstein's general LogToReal (exact roots
         * of every Q_i, Rioboo's LogToAtan), Root objects allowed. */
        if (!result && real_field && rf == INTRAT_REALFORM_LOGTOREAL)
            result = intrat_alg_log_to_real(pairs, x, t_sym, true);
        /* "Rioboo" and Automatic: the shape dispatch (quadratic factors in
         * closed form, sparse and biquadratic quartics, ...).  Automatic does
         * not add the general LogToReal: on benchmark 266 it was never smaller
         * than this dispatch's answer (or the RootSum the reference keeps) and
         * only cost time. */
        if (!result && real_field && rf != INTRAT_REALFORM_NONE)
            result = intrat_log_to_real_pairs(pairs, x, t_sym);
        /* LRT consumer 2: linear-Q closer (factor each Q_i over Q); a
         * non-real field has tried it above. */
        if (!result && real_field && rf != INTRAT_REALFORM_NONE)
            result = intrat_linear_q_closer(pairs, x, t_sym);
        expr_free(pairs);
    }
    expr_free(t_sym);

    /* Non-LRT universal fallback (Lagrange RootSum on full a/d). */
    if (!result) result = intrat_naive_log_part(h_eval, x);

    expr_free(h_eval);
    return result;
}

/* Process h piece-by-piece via Apart, ExtractConstants, and the
 * derivative-recognition / LRT pipeline.  Returns the closed sum or
 * NULL if any piece can't be closed. */
static Expr* intrat_integrate_summands(Expr* h, Expr* x) {
    Expr* pieces = intrat_apart_list(h, x, NULL);
    if (!pieces || pieces->type != EXPR_FUNCTION
        || pieces->data.function.head->type != EXPR_SYMBOL
        || pieces->data.function.head->data.symbol.name != SYM_List) {
        if (pieces) expr_free(pieces);
        return NULL;
    }
    size_t n = pieces->data.function.arg_count;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t nterms = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* piece = pieces->data.function.args[i];

        /* Skip pure polynomial summands — they belong to the Q part. */
        if (intrat_polyq_test(piece, x)) {
            Expr* poly_int = intrat_integrate_polynomial(piece, x);
            if (!poly_int) goto fail;
            terms[nterms++] = poly_int;
            continue;
        }

        /* ExtractConstants — pull scalar prefactors out of num/den. */
        Expr* ec = intrat_extract_constants(piece, x);
        Expr* c_const = list_get(ec, 1);
        Expr* residue = list_get(ec, 2);
        expr_free(ec);
        if (!c_const || !residue) {
            if (c_const) expr_free(c_const);
            if (residue) expr_free(residue);
            goto fail;
        }

        Expr* num = intrat_numerator(residue);
        Expr* den = intrat_denominator(residue);
        Expr* num_e = eval_and_free(num);
        Expr* den_e = eval_and_free(den);
        expr_free(residue);

        /* Derivative-recognition fast path. */
        Expr* fast = intrat_derivative_recognition(num_e, den_e, x);
        Expr* piece_int = NULL;
        if (fast) piece_int = fast;
        else      piece_int = try_lrt_close(num_e, den_e, x);
        expr_free(num_e); expr_free(den_e);

        if (!piece_int) {
            expr_free(c_const);
            goto fail;
        }
        Expr* scaled = internal_times(
            (Expr*[]){c_const, piece_int}, 2);
        terms[nterms++] = eval_and_free(scaled);
    }

    expr_free(pieces);
    if (nterms == 0) { free(terms); return expr_new_integer(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    Expr* sum = internal_plus(terms, nterms);
    free(terms);
    return eval_and_free(sum);

fail:
    expr_free(pieces);
    for (size_t k = 0; k < nterms; k++) expr_free(terms[k]);
    free(terms);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* BronsteinRational top-level (Phase 1 skeleton + Phase 2 hookup).    */
/* ------------------------------------------------------------------ */

/* True when e contains a Root object anywhere. */
static bool intrat_tree_has_root(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Root) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (intrat_tree_has_root(e->data.function.args[i])) return true;
    return false;
}

static Expr* intrat_integrate_rational(Expr* f, Expr* x) {
    intrat_trace("BronsteinRational", "IN", f);

    /* Pure-polynomial fast path. */
    if (intrat_polyq_test(f, x)) {
        return intrat_integrate_polynomial(f, x);
    }

    /* Combine into a single fraction num/den. */
    Expr* tg = internal_together((Expr*[]){expr_copy(f)}, 1);
    Expr* combined = evaluate(tg);
    expr_free(tg);
    if (!combined) return NULL;

    Expr* num = intrat_numerator(combined);
    Expr* num_eval = eval_and_free(num);
    Expr* den = intrat_denominator(combined);
    Expr* den_eval = eval_and_free(den);
    expr_free(combined);

    if (!intrat_polyq_test(num_eval, x) || !intrat_polyq_test(den_eval, x)) {
        expr_free(num_eval); expr_free(den_eval);
        return NULL;
    }

    /* Polynomial part: Q + R/den. */
    Expr *Q = NULL, *R = NULL;
    if (!intrat_polyqr(num_eval, den_eval, x, &Q, &R)) {
        expr_free(num_eval); expr_free(den_eval);
        return NULL;
    }
    expr_free(num_eval);

    Expr* poly_int = expr_new_integer(0);
    if (!is_zero_poly(Q)) {
        Expr* tmp = intrat_integrate_polynomial(Q, x);
        if (tmp) { expr_free(poly_int); poly_int = tmp; }
    }
    expr_free(Q);

    if (is_zero_poly(R)) {
        expr_free(R); expr_free(den_eval);
        Expr* result = eval_and_free(poly_int);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Try the derivative-recognition fast path on R/den. */
    Expr* fast = intrat_derivative_recognition(R, den_eval, x);
    if (fast) {
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, fast}, 2);
        Expr* result = eval_and_free(sum);
        expr_free(R); expr_free(den_eval);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* HermiteReduce R/den. */
    Expr* proper = internal_divide(
        (Expr*[]){R, den_eval}, 2);
    Expr* hr = intrat_hermite_reduce(proper, x);
    expr_free(proper);
    if (!hr || hr->type != EXPR_FUNCTION || hr->data.function.arg_count != 2) {
        if (hr) expr_free(hr);
        expr_free(poly_int);
        return NULL;
    }
    Expr* g = expr_copy(hr->data.function.args[0]);
    Expr* h = expr_copy(hr->data.function.args[1]);
    expr_free(hr);

    /* The Hermite rational part, in its most compact form (and kept out of
     * the Expand-based log-part passes below, which split its numerators
     * over their denominators). */
    {
        Expr* gc = intsimp_compact_rational(g, x);
        if (gc) { expr_free(g); g = gc; }
    }

    /* If h is zero, we have a closed form. */
    if (is_zero_poly(h)) {
        expr_free(h);
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, g}, 2);
        Expr* result = eval_and_free(sum);
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Phase 5 — IntegrateRealRationalFunction-style per-summand loop.
     *
     * We process the squarefree-denominator residual h piece-by-piece:
     * Apart-split it, ExtractConstants to keep coefficients out of the
     * subresultant chain, then run derivative-recognition / LRT /
     * LogToReal on each piece.  If every piece closes we sum them up
     * with the Hermite g and the polynomial part; otherwise we fall
     * back to running LRT on h whole (Phase 4 path).
     */
    Expr* per_summand = intrat_integrate_summands(h, x);
    if (per_summand) {
        expr_free(h);
        /* g joins after the log-part passes (see intsimp_compact_rational). */
        Expr* sum = internal_plus(
            (Expr*[]){poly_int, per_summand}, 2);
        Expr* sum_eval = eval_and_free(sum);
        /* Phase A-1: distribute scalar Times-over-Plus so per-summand
         * `c_k * piece_int_k` accumulators flatten into a single
         * Plus[...] of `c_k · term` summands.  Without this, the
         * subsequent log-pairing rule (intsimp_log_to_arctanh) cannot
         * see across the Times boundary.
         * Phase A-2: strip Log[c · p] -> Log[p] when FreeQ[c, x] so
         * downstream Log[A] + Log[B] merging sees the cleaned
         * arguments. */
        /* Coefficients in Q(alpha) for a Root object alpha: reduce
         * them modulo its minimal polynomial before the Log-pairing
         * and constant-stripping passes read the kernel arguments. */
        Expr* sum_red = intrat_reduce_root_coefficients(sum_eval, x);
        expr_free(sum_eval);
        Expr* distributed = intsimp_distribute_plus(sum_red);
        expr_free(sum_red);
        Expr* stripped = intsimp_strip_log_constants(distributed, x);
        expr_free(distributed);
        /* Phase 6 final pass: combine c Log[A] ± c Log[B] terms.
         * Phase 7 final pass: normalise ArcTan/ArcTanh argument signs.*/
        Expr* simplified = intsimp_log_to_arctanh(stripped, x);
        expr_free(stripped);
        /* Re-distribute and re-strip in case the pairing rule
         * introduces new Times[scalar, Plus[…]] or Log[c · p] shapes. */
        Expr* simplified2 = intsimp_distribute_plus(simplified);
        expr_free(simplified);
        Expr* simplified3 = intsimp_strip_log_constants(simplified2, x);
        expr_free(simplified2);
        /* Normalise the algebraic constants (Sqrt[2], ...) the log part
         * leaves unreduced: kernel arguments get RootReduce-normal
         * coefficients and terms sharing a kernel merge. */
        Expr* reduced = intsimp_reduce_algebraic_constants(simplified3, x);
        expr_free(simplified3);
        Expr* normalized = intsimp_normalize_inverse_trig_signs(reduced);
        expr_free(reduced);
        /* Same for a Root object: coefficients reduced modulo its
         * minimal polynomial, kernel arguments monic in x. */
        Expr* result = intrat_reduce_root_coefficients(normalized, x);
        if (!expr_eq(result, normalized)) {
            /* the reduction pulls rational content to the front of
             * Log arguments: drop it as the earlier passes do */
            Expr* stripped_r = intsimp_strip_log_constants(result, x);
            expr_free(result);
            result = stripped_r;
        }
        expr_free(normalized);
        if (!is_zero_poly(g)) {
            /* g over Q(alpha) for a Root alpha gets the same reduction the
             * log part had above */
            if (intrat_tree_has_root(g)) {
                Expr* gr = intrat_reduce_root_coefficients(g, x);
                if (gr) { expr_free(g); g = gr; }
            }
            result = eval_and_free(internal_plus((Expr*[]){ g, result }, 2));
        } else {
            expr_free(g);
        }
        intrat_trace("BronsteinRational", "OUT", result);
        return result;
    }

    /* Beyond Phase 5's combined scope — leave unevaluated. */
    expr_free(h); expr_free(g); expr_free(poly_int);
    intrat_trace("BronsteinRational", "OUT (unresolved)", f);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* BuiltinFunc wrappers.  Each follows the Mathilda convention: caller  */
/* owns res; return NULL on failure.                                   */
/* ------------------------------------------------------------------ */

/* Recognise a trailing option (Rule[lhs, rhs] / RuleDelayed[...]).
 * Returns true if `opt` is a Rule with one of the BronsteinRational
 * option names on the lhs. */
static bool is_intrat_option(Expr* opt) {
    if (!opt || opt->type != EXPR_FUNCTION
        || opt->data.function.head->type != EXPR_SYMBOL
        || (opt->data.function.head->data.symbol.name != SYM_Rule
            && opt->data.function.head->data.symbol.name != SYM_RuleDelayed)
        || opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    if (lhs->type == EXPR_SYMBOL) {
        return lhs->data.symbol.name == SYM_Extension;
    }
    if (lhs->type == EXPR_STRING) {
        const char* s = lhs->data.string;
        return strcmp(s, "PFD") == 0 || strcmp(s, "LogToArcTan") == 0
            || strcmp(s, "Radicals") == 0 || strcmp(s, "Hermite") == 0
            || strcmp(s, "LogPart") == 0 || strcmp(s, "RealForm") == 0;
    }
    return false;
}

/* Recursively scan for any EXPR_REAL leaf.  Mirrors the same guard in
 * src/integrate.c; here it is silent (no Integrate::inexact message)
 * so the user-visible diagnostic comes from the public Integrate head. */
static bool intrat_contains_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_FUNCTION) {
        if (intrat_contains_real(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (intrat_contains_real(e->data.function.args[i])) return true;
        }
    }
    return false;
}

Expr* builtin_intrat_integraterational(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Strip recognised trailing options.  Phase 7 stores option
     * values via extract_bool_option / extract_extension_option but
     * keeps them advisory for now — the algorithmic path uses the
     * defaults. */
    size_t argc = res->data.function.arg_count;
    IntratMethodOptions mo;
    intrat_method_options_default(&mo);
    bool have_mo = false;
    while (argc > 0 && is_intrat_option(res->data.function.args[argc - 1])) {
        if (intrat_method_suboption(res->data.function.args[argc - 1], &mo)) have_mo = true;
        argc--;
    }
    if (argc != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    IntratMethodOptions saved_mo;
    if (have_mo) intrat_method_options_install(&mo, &saved_mo);
    /* Try to coerce inexact integrands to exact rationals.  Silent
     * here — the public Integrate head owns the user-visible
     * diagnostic.  If Rationalize cannot eliminate every EXPR_REAL,
     * bubble back unevaluated. */
    Expr* coerced = NULL;
    if (intrat_contains_real(f)) {
        coerced = internal_rationalize_expr(f, 0.0, RATIONALIZE_DEFAULT);
        if (!coerced || intrat_contains_real(coerced)) {
            if (coerced) expr_free(coerced);
            if (have_mo) intrat_method_options_restore(&saved_mo);
            return NULL;
        }
        f = coerced;
    }
    Expr* result = intrat_integrate_rational(f, x);
    if (coerced) expr_free(coerced);
    if (have_mo) intrat_method_options_restore(&saved_mo);
    return result;
}

Expr* builtin_intrat_integratepolynomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    if (!intrat_polyq_test(f, x)) return NULL;
    return intrat_integrate_polynomial(f, x);
}

Expr* builtin_intrat_hermitereduce(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_hermite_reduce(f, x);
}

Expr* builtin_intrat_helpers_content(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    if (is_zero_poly(p)) return expr_new_integer(0);
    return intrat_content(p, x);
}

Expr* builtin_intrat_helpers_primitive(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_primitive(p, x);
}

Expr* builtin_intrat_helpers_monic(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_monic_internal(p, x);
}

Expr* builtin_intrat_helpers_leadingcoefficient(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_lc(p, x);
}

/* Phase 2 builtins. */

Expr* builtin_intrat_intrationallogpart(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    /* Strip a trailing RootSum -> True/False option. */
    size_t argc = res->data.function.arg_count;
    bool root_sum = false;
    if (argc >= 1) {
        Expr* last = res->data.function.args[argc - 1];
        if (last->type == EXPR_FUNCTION
            && last->data.function.head->type == EXPR_SYMBOL
            && (last->data.function.head->data.symbol.name == SYM_Rule
                || last->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && last->data.function.arg_count == 2
            && last->data.function.args[0]->type == EXPR_SYMBOL
            && last->data.function.args[0]->data.symbol.name == SYM_RootSum) {
            Expr* val = last->data.function.args[1];
            root_sum = (val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_True);
            argc--;
        }
    }
    if (argc != 3) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    Expr* t = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    return intrat_int_rational_log_part(f, x, t, root_sum);
}

Expr* builtin_intrat_helpers_squarefree(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return intrat_squarefree_list(res->data.function.args[0]);
}

Expr* builtin_intrat_helpers_extractconstants(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_extract_constants(f, x);
}

Expr* builtin_intrat_helpers_apartlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    const Expr* alpha = extract_extension_option(res, &argc);
    if (argc != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_apart_list(f, x, alpha);
}

Expr* builtin_intrat_logtoatan(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_log_to_atan(a, b, x);
}

Expr* builtin_intrat_logtoreal(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* r = res->data.function.args[0];
    Expr* s = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* t = res->data.function.args[3];
    if (x->type != EXPR_SYMBOL || t->type != EXPR_SYMBOL) return NULL;
    return intrat_log_to_real(r, s, x, t);
}

Expr* builtin_intrat_transcendental_log_part(Expr* res) {
    if (res->type != EXPR_FUNCTION
        || (res->data.function.arg_count != 6 && res->data.function.arg_count != 7))
        return NULL;
    Expr* a    = res->data.function.args[0];
    Expr* d    = res->data.function.args[1];
    Expr* tau  = res->data.function.args[2];
    Expr* z    = res->data.function.args[3];
    Expr* Dd   = res->data.function.args[4];
    Expr* xreal = res->data.function.args[5];
    /* Optional 7th argument is the "decide" marker (any expression): when present
     * the residue criterion, on an all-non-constant-residue resultant, returns
     * Integrate`PartialLogPart[0, a/d] instead of declining — so a caller seeking a
     * non-elementarity VERDICT (never re-integrating the remainder) learns the whole
     * simple part is non-elementary (Thm 5.6.1(ii)).  Absent for the integrator. */
    bool decide = (res->data.function.arg_count == 7);
    if (tau->type != EXPR_SYMBOL || z->type != EXPR_SYMBOL)
        return NULL;
    /* `xreal` is the residue-constant gate: a single symbol (single-kernel
     * extension) or a List of symbols (tower proper part — free of x and every
     * lower tower variable). */
    bool gate_ok = (xreal->type == EXPR_SYMBOL)
        || (xreal->type == EXPR_FUNCTION
            && xreal->data.function.head->type == EXPR_SYMBOL
            && xreal->data.function.head->data.symbol.name == SYM_List);
    if (!gate_ok) return NULL;
    Expr* rem = NULL;
    Expr* logs = intrat_transcendental_log_part(a, d, tau, z, Dd, xreal, &rem, decide);
    if (!logs) { if (rem) expr_free(rem); return NULL; }
    if (rem) {
        /* Partial log part: wrap the elementary logs and the unintegrated
         * non-constant-residue remainder in an inert head the recursive-Risch
         * caller (rt_field_lrt_logpart) decodes.  The head has no definition, so
         * it survives evaluation unchanged. */
        return expr_new_function(expr_new_symbol("Integrate`PartialLogPart"),
                                 (Expr*[]){ logs, rem }, 2);
    }
    return logs;
}

Expr* builtin_intrat_logtoarctanh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* e = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intsimp_log_to_arctanh(e, x);
}

Expr* builtin_intrat_naive_log_part(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return intrat_naive_log_part(f, x);
}

/* Integrate`RothsteinTragerResultant[num, den, z, x]
 *   = Resultant[num - z D[den, x], den, x],
 * the Rothstein-Trager resultant whose roots in z are the residues of num/den.
 * This is the named parametric-resultant primitive the LRT log part uses and that
 * Cherry's special-function argument generators reuse (CHERRY_DESIGN.md §3.3): the
 * log / exponential-integral (Ei) arguments are recovered from its roots. Built on
 * the existing `Resultant` builtin — no new resultant machinery. */
static Expr* builtin_intrat_rt_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];
    Expr* z   = res->data.function.args[2];
    Expr* x   = res->data.function.args[3];
    Expr* dd = eval_and_free(expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy(den), expr_copy(x) }, 2));            /* D[den, x] */
    if (!dd) return NULL;
    Expr* arg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(num),
            expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_copy(z), dd }, 3) }, 2)); /* num - z D[den] */
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Resultant),
        (Expr*[]){ arg, expr_copy(den), expr_copy(x) }, 3));       /* Resultant[.,den,x] */
    return r;
}

/* Integrate`ExpIntegralEiResultant[g1, p, q, alpha, x]
 *   = Resultant[g1, p + alpha q, x],
 * the parametric resultant whose roots in the parameter alpha are Cherry's
 * ei-argument constants (the u_i = f + alpha_i of the P1 generator, resolving
 * g~1/g1 = Sum c_i u_i'/u_i; Cherry 1989 §4.2 / CHERRY_PLAN.md §3.2).  Distinct
 * from RothsteinTragerResultant (which resultants against num - z D[den]); here
 * the parameter enters linearly through p + alpha q with q the exponent
 * denominator.  A thin, testable sibling built on the existing Resultant. */
static Expr* builtin_intrat_ei_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 5) return NULL;
    Expr* g1    = res->data.function.args[0];
    Expr* p     = res->data.function.args[1];
    Expr* q     = res->data.function.args[2];
    Expr* alpha = res->data.function.args[3];
    Expr* x     = res->data.function.args[4];
    Expr* arg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(p),
            expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(alpha), expr_copy(q) }, 2) }, 2));  /* p + alpha q */
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Resultant),
        (Expr*[]){ expr_copy(g1), arg, expr_copy(x) }, 3));             /* Resultant[g1,.,x] */
    return r;
}

/* ------------------------------------------------------------------ */
/* Init.                                                              */
/* ------------------------------------------------------------------ */

static void install(const char* name, Expr* (*fn)(Expr*), const char* docstring) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    if (docstring) symtab_set_docstring(name, docstring);
}

void intrat_init(void) {
    /* Package symbols live under `Integrate``. Names must be the fully
     * qualified strings — the symtab is keyed by interned name and the
     * parser routes `Integrate`Foo` directly to that key. */
    install("Integrate`BronsteinRational",
            builtin_intrat_integraterational,
            "Integrate`BronsteinRational[f, x, opts] integrates a rational function f\n"
            "of x (Bronstein, Symbolic Integration I, ch. 2); a non-rational f stays\n"
            "unevaluated. Options, also given as Method -> {\"BronsteinRational\", ...}:\n"
            "\"Hermite\" -> \"Linear\" | \"Quadratic\" | \"HorowitzOstrogradsky\";\n"
            "\"LogPart\" -> \"LazardRiobooTrager\" | \"RothsteinTrager\" | \"Czichowski\";\n"
            "\"RealForm\" -> \"Rioboo\" | \"LogToReal\" | None.  Each defaults to Automatic.");

    install("Integrate`IntegratePolynomial",
            builtin_intrat_integratepolynomial,
            "Integrate`IntegratePolynomial[p, x] integrates a polynomial p in x\n"
            "term-by-term using the rule a x^n -> a x^(n+1)/(n+1). Returns\n"
            "unevaluated for non-polynomial inputs.");

    install("Integrate`HermiteReduce",
            builtin_intrat_hermitereduce,
            "Integrate`HermiteReduce[f, x] returns {g, h} such that\n"
            "f == D[g, x] + h with Denominator[h] squarefree (Mack's linear\n"
            "version of Hermite reduction; Bronstein, Symbolic Integration I,\n"
            "p. 44).");

    install("Integrate`Helpers`Content",
            builtin_intrat_helpers_content,
            "Integrate`Helpers`Content[p, x] gives PolynomialGCD applied to\n"
            "the coefficient list of p as a polynomial in x.");

    install("Integrate`Helpers`Primitive",
            builtin_intrat_helpers_primitive,
            "Integrate`Helpers`Primitive[p, x] is p divided by Content[p, x],\n"
            "the primitive-part operation in the polynomial ring K[x].");

    install("Integrate`Helpers`Monic",
            builtin_intrat_helpers_monic,
            "Integrate`Helpers`Monic[p, x] divides p by its leading coefficient\n"
            "in x, returning a monic polynomial. For p free of x returns 1.");

    install("Integrate`Helpers`LeadingCoefficient",
            builtin_intrat_helpers_leadingcoefficient,
            "Integrate`Helpers`LeadingCoefficient[p, x] gives the coefficient\n"
            "of the highest power of x in p; for p free of x returns p.");

    /* Phase 2 — Lazard-Rioboo-Trager log part. */
    install("Integrate`IntRationalLogPart",
            builtin_intrat_intrationallogpart,
            "Integrate`IntRationalLogPart[A/D, x, t] computes the logarithmic\n"
            "part of an integral of a rational function with squarefree\n"
            "denominator D, returning a list of {Q_i(t), S_i(t,x)} pairs\n"
            "(default) or a sum of RootSum heads (option RootSum -> True).\n"
            "Bronstein, Symbolic Integration I, p. 51 (Lazard-Rioboo-Trager).");

    install("Integrate`Helpers`SquareFree",
            builtin_intrat_helpers_squarefree,
            "Integrate`Helpers`SquareFree[p] returns a list of\n"
            "{factor, multiplicity} pairs indexed densely by multiplicity\n"
            "1..max, with the i-th entry holding the product of all\n"
            "squarefree factors of multiplicity i (or {1, i} if none).");

    install("Integrate`Helpers`ExtractConstants",
            builtin_intrat_helpers_extractconstants,
            "Integrate`Helpers`ExtractConstants[f, x] returns {const, residue}\n"
            "with const free of x and residue holding the remaining\n"
            "polynomial-in-x parts of Numerator[f] and Denominator[f].");

    install("Integrate`Helpers`ApartList",
            builtin_intrat_helpers_apartlist,
            "Integrate`Helpers`ApartList[f, x] returns Apart[f, x]'s output\n"
            "as a List of summands. Accepts a trailing Extension -> alpha\n"
            "rule that is forwarded to the underlying Apart call.");

    /* Phase 3 — Rioboo's recursive LogToAtan. */
    install("Integrate`LogToAtan",
            builtin_intrat_logtoatan,
            "Integrate`LogToAtan[A, B, x] returns a sum of ArcTan of polynomials\n"
            "in K[x] whose derivative equals D[I Log[(A + I B)/(A - I B)], x].\n"
            "Used by LogToReal (Phase 4) to convert complex log pairs into real\n"
            "ArcTan summands. Bronstein, Symbolic Integration I, p. 63 (Rioboo).");

    /* Phase 6 — Log post-processing. */
    install("Integrate`LogToArcTanh",
            builtin_intrat_logtoarctanh,
            "Integrate`LogToArcTanh[expr, x] combines pairs of Log terms in\n"
            "expr's top-level Plus head: c Log[A] + c Log[B] -> c Log[A*B];\n"
            "c Log[A] - c Log[B] -> 2 c ArcTanh[(B-A)/(B+A)] when the\n"
            "ArcTanh argument is rational in x, else c Log[A/B] when A/B is\n"
            "rational in x. Beautifies the rational integrator's output\n"
            "without affecting differentiation; cf. IntegrateRational.m\n"
            "lines 1722-1761.");

    /* Phase 4 — Rioboo's LogToReal. */
    install("Integrate`LogToReal",
            builtin_intrat_logtoreal,
            "Integrate`LogToReal[r, s, x, t] converts the symbolic\n"
            "RootSum[r(α)==0, α Log[s(α, x)]] into a real Log + ArcTan form\n"
            "by factoring r over Q and dispatching each linear / quadratic\n"
            "factor (positive discriminant -> two real Logs; negative\n"
            "discriminant -> Log + LogToAtan-derived ArcTan). Returns the\n"
            "call unevaluated when r has factors of degree > 2 in t (the\n"
            "bounded-Solve scope of Phase 4).");

    /* Recursive-Risch transcendental LRT log part (consumed by
     * Integrate`RischTranscendental; not a user-facing surface form). */
    install("Integrate`RothsteinTragerResultant",
            builtin_intrat_rt_resultant,
            "Integrate`RothsteinTragerResultant[num, den, z, x] = Resultant[num - z D[den,x],\n"
            "den, x], the Rothstein-Trager resultant whose roots in z are the residues of\n"
            "num/den. Argument generator for the log / exponential-integral (Ei) parts of\n"
            "Cherry's special-function integration.");

    install("Integrate`ExpIntegralEiResultant",
            builtin_intrat_ei_resultant,
            "Integrate`ExpIntegralEiResultant[g1, p, q, a, x] = Resultant[g1, p + a q, x],\n"
            "the parametric resultant whose roots in the parameter a are the ExpIntegralEi\n"
            "argument constants alpha_i of Cherry's P1 generator (u_i = f + alpha_i, f = p/q).\n"
            "Built on Resultant; the ei sibling of RothsteinTragerResultant.");

    install("Integrate`TranscendentalLogPart",
            builtin_intrat_transcendental_log_part,
            "Integrate`TranscendentalLogPart[a, d, tau, z, Dd, g] computes the\n"
            "logarithmic part of Integrate[a/d] for a transcendental monomial\n"
            "tau (a/d proper with squarefree d in tau), where Dd is the monomial\n"
            "derivation D(d), z is the residue variable, and g is the residue-\n"
            "constant gate: a symbol (single kernel) or a List of symbols (tower\n"
            "proper part, gating x and every lower tower variable). Returns the\n"
            "real Log + ArcTan form as a function of tau (Rothstein-Trager\n"
            "resultant + Rioboo LogToReal), or unevaluated when the integral is\n"
            "not elementary in this form.");

    /* Phase 8b — NaiveLogPart RootSum fallback. */
    install("Integrate`NaiveLogPart",
            builtin_intrat_naive_log_part,
            "Integrate`NaiveLogPart[f, x] returns the held-symbolic\n"
            "RootSum form of the logarithmic part of Integrate[f, x]\n"
            "for a rational function f = a/d:\n"
            "  RootSum[Function[t, d(t)],\n"
            "          Function[t, a(t) Log(x - t) / d'(t)]].\n"
            "Used as the universal fallback by the rational integrator\n"
            "whenever LogToReal cannot close the log part to a real\n"
            "elementary expression.  The result differentiates back to f\n"
            "via the D[RootSum, x] rule wired in src/deriv.c.\n"
            "Direct port of IntegrateRational.m:1116-1124.");

    /* Trace flag: Integrate`$Verbose. Default False. */
    {
        Expr* pat = expr_new_symbol("Integrate`$Verbose");
        Expr* val = expr_new_symbol(SYM_False);
        symtab_add_own_value("Integrate`$Verbose", pat, val);
        expr_free(pat); expr_free(val);
    }

    /* RootSum is registered in src/root.c (Phase 8b-prereq) — the
     * canonical home of held algebraic-root machinery.  Nothing to
     * do here. */
}
