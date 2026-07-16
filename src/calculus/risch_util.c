/* risch_util.c — shared low-level helpers for the transcendental Risch modules.
 *
 * Small host-grounded eval / predicate / polynomial primitives (the `rt_` core)
 * shared across integrate_risch_transcendental.c and its tower / RDE / case
 * modules.  See risch_util.h for the public surface and the ownership notes on
 * each wrapper.
 */

#include "risch_util.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "sym_intern.h"
#include "simp_trigexp_zero.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

/* ================================================================== */
/* Small evaluation helpers.                                          */
/* ================================================================== */

/* True iff `e` is the unevaluated call `name[...]`. */
bool rt_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == intern_symbol(name);
}

/* True iff `e` is a nonzero integer constant; writes its value to *out. */
static bool rt_is_int_const(const Expr* e, long* out) {
    if (e && e->type == EXPR_INTEGER && e->data.integer != 0) {
        *out = (long)e->data.integer; return true;
    }
    return false;
}

/* True iff `e` is a nonzero rational number constant (Integer, BigInt, or
 * Rational[p, q]) — used to group multiplicatively commensurate exponents. */
bool rt_is_rat_const(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer != 0;
    if (e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == intern_symbol("Rational");
}

/* Numerator / denominator (den > 0) of a rational-constant Expr (Integer or
 * Rational[p, q]).  Returns false if `e` is not a small rational constant
 * (BigInt rationals are out of range here and decline).  Used to synthesize a
 * commensurability-class primitive exponent from non-integer rational ratios. */
static bool rt_rat_num_den(const Expr* e, long* num, long* den) {
    long i;
    if (rt_is_int_const(e, &i)) { *num = i; *den = 1; return true; }
    if (e && e->type == EXPR_FUNCTION && rt_head_is(e, "Rational") &&
        e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p && p->type == EXPR_INTEGER && q && q->type == EXPR_INTEGER &&
            q->data.integer != 0) {
            long a = (long)p->data.integer, b = (long)q->data.integer;
            if (b < 0) { a = -a; b = -b; }
            *num = a; *den = b; return true;
        }
    }
    return false;
}

static long rt_gcd_l(long a, long b) {
    a = labs(a); b = labs(b);
    while (b) { long t = a % b; a = b; b = t; }
    return a;
}
static long rt_lcm_l(long a, long b) {
    if (a == 0 || b == 0) return 0;
    return labs(a / rt_gcd_l(a, b) * b);
}

/* Build `head[args...]` (adopting the owned `args` element pointers) and
 * evaluate it, freeing the constructed call.  Returns evaluate()'s result. */
Expr* rt_eval_call(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Evaluate `e`, taking ownership: frees `e` and returns evaluate()'s
 * result (evaluate itself does not consume its argument). */
Expr* rt_eval_own(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* ================================================================== */
/* Generic builtin-call and predicate helpers (grounded in the host).  */
/* ================================================================== */

/* NOTE ON CORRECTNESS.  This is a decision procedure: every branch is
 * correct by construction and is NOT checked by differentiating the
 * result.  Each case fires only behind an exact structural certificate
 * (Cancel of a logarithmic derivative to a polynomial, an exactly linear
 * denominator, a genuinely x-free cofactor, ...) that already proves the
 * closed form it emits.                                                */

Expr* rt_eval1(const char* head, Expr* a) {
    return rt_eval_call(head, (Expr*[]){ a }, 1);
}
Expr* rt_eval2(const char* head, Expr* a, Expr* b) {
    return rt_eval_call(head, (Expr*[]){ a, b }, 2);
}
Expr* rt_eval3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rt_eval_call(head, (Expr*[]){ a, b, c }, 3);
}

/* Given a commensurability class of x-dependent exponents `ws[0..nw-1]` (every
 * ratio ws[j]/ws[0] a rational constant), synthesize the class *primitive*
 * exponent p = ws[0] / M (M = lcm of the ratio denominators) so that every
 * E^(ws[j]) = (E^p)^(kof[j]) with kof[j] a nonzero integer.  This generalizes
 * the earlier "some member is an integer multiple of every other" test — which
 * only handled integer ratios (E^(2u) = (E^u)^2) — to genuine rational ratios
 * (E^(x/2), E^(x/3) → primitive E^(x/6), k = 3, 2).  Returns an owned,
 * Cancel-normalized primitive Expr* and fills kof[], or NULL if any ratio is
 * not a small rational constant (genuinely independent kernels — a tower /
 * sum, handled elsewhere) or on arithmetic overflow. */
Expr* rt_class_primitive(Expr** ws, size_t nw, long* kof) {
    if (nw == 0) return NULL;
    long* nums = malloc(nw * sizeof(long));
    long* dens = malloc(nw * sizeof(long));
    long M = 1;
    bool ok = true;
    for (size_t j = 0; j < nw && ok; j++) {
        Expr* ratio = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(ws[j]), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(ws[0]), expr_new_integer(-1) }, 2) }, 2));
        long a, b;
        if (ratio && rt_rat_num_den(ratio, &a, &b) && a != 0) {
            nums[j] = a; dens[j] = b; M = rt_lcm_l(M, b);
        } else ok = false;
        if (ratio) expr_free(ratio);
    }
    Expr* prim = NULL;
    if (ok && M > 0) {
        for (size_t j = 0; j < nw; j++) kof[j] = nums[j] * (M / dens[j]);
        prim = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(ws[0]), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(M), expr_new_integer(-1) }, 2) }, 2));
    }
    free(nums); free(dens);
    return prim;
}

bool rt_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("True");
}
/* FreeQ[e, x] */
bool rt_free_of_x(Expr* e, Expr* x) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* FreeQ[e, head] — true iff the symbol `head` occurs nowhere in e. */
bool rt_free_of_head(Expr* e, const char* head) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_new_symbol(head));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* PolynomialQ[e, x] */
bool rt_is_poly(Expr* e, Expr* x) {
    Expr* r = rt_eval2("PolynomialQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* Degree of a polynomial in x via Length[CoefficientList[e, x]] - 1
 * (Mathilda has no Exponent builtin).  Returns -1 if CoefficientList does
 * not reduce to a List.  Call only after rt_is_poly. */
long rt_degree(Expr* e, Expr* x) {
    Expr* cl = rt_eval2("CoefficientList", expr_copy(e), expr_copy(x));
    long d = -1;
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List"))
        d = (long)cl->data.function.arg_count - 1;
    if (cl) expr_free(cl);
    return d;
}
/* Coefficient[e, x, k] */
Expr* rt_coeff(Expr* e, Expr* x, long k) {
    return rt_eval3("Coefficient", expr_copy(e), expr_copy(x),
                    expr_new_integer(k));
}
/* Together[e] === 0 — an exact zero test for rational functions of the
 * field kernels (Together is exact and, unlike Simplify, cheap; it is only
 * ever applied here to constants and small rational cofactors). */
/* Is `e` identically zero?  A structural recursion handles the shapes that
 * dominate the Risch path — a single fraction Times[num, Power[den, -1]] or a
 * polynomial — WITHOUT combining over a common denominator: a product is zero
 * iff some factor is zero, and a reciprocal factor Power[g, negative] = 1/g is
 * never zero, so its (possibly degree-2n) denominator is never examined. Only a
 * Plus (whose terms might cancel) or an opaque head falls back to the exact
 * Together test, and then only on the smaller numerator polynomials. This keeps
 * the zero-test off the big rational denominators that made it cost ~0.5 s. */
bool rt_is_zero(Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer == 0;
        case EXPR_REAL:    return e->data.real == 0.0;
        case EXPR_BIGINT:  return false;   /* normalized: |value| >= 2^63 != 0 */
        case EXPR_SYMBOL:
        case EXPR_STRING:  return false;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL) {
        if (rt_head_is(e, "Times")) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* fac = e->data.function.args[i];
                if (rt_head_is(fac, "Power") && fac->data.function.arg_count == 2) {
                    Expr* ex = fac->data.function.args[1];
                    if (ex->type == EXPR_INTEGER && ex->data.integer < 0) continue;
                    if (rt_head_is(ex, "Rational") && ex->data.function.arg_count == 2 &&
                        ex->data.function.args[0]->type == EXPR_INTEGER &&
                        ex->data.function.args[0]->data.integer < 0) continue;
                }
                if (rt_is_zero(fac)) return true;
            }
            return false;
        }
        if (rt_head_is(e, "Power") && e->data.function.arg_count == 2) {
            Expr* ex = e->data.function.args[1];
            if (ex->type == EXPR_INTEGER && ex->data.integer < 0) return false;
            if (rt_head_is(ex, "Rational") && ex->data.function.arg_count == 2 &&
                ex->data.function.args[0]->type == EXPR_INTEGER &&
                ex->data.function.args[0]->data.integer < 0) return false;
            return rt_is_zero(e->data.function.args[0]);
        }
    }
    Expr* s = rt_eval1("Together", expr_copy(e));
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* Parse `tmpl`, substitute the named placeholder symbols with `vals`
 * (borrowed; copied in), evaluate, and return the result (or NULL). */
Expr* rt_template(const char* tmpl, const char** names,
                         Expr** vals, size_t n) {
    Expr* t = parse_expression(tmpl);
    if (!t) return NULL;
    Expr** rules = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_symbol(names[i]), expr_copy(vals[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, n);
    free(rules);
    Expr* ra = expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ t, rl }, 2);
    Expr* out = evaluate(ra);
    expr_free(ra);
    return out;
}
/* Find (borrowed) the argument u of the first Log[u] subexpression of `e`
 * whose argument depends on x. */
Expr* rt_find_log_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    Expr* r = rt_find_log_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_log_of_x(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

/* Find (borrowed) the exponent u of the first exponential kernel of `e`
 * whose exponent depends on x, matching either Exp[u] or E^u (Power[E,u]). */
Expr* rt_find_exp_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == intern_symbol("E")
            && !rt_free_of_x(ex, x))
            return ex;
    }
    Expr* r = rt_find_exp_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_exp_of_x(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

/* Diff-back verifier shared across the transcendental Risch modules: is
 * Simplify[D[result,x] - f] identically 0?  Fast exact path first
 * (Together[TrigToExp[...]] and the rational-trig point-grid), Simplify last.
 * Correct-by-construction cases skip this; search-based ones gate on it. */
bool rt_verify_antideriv(Expr* result, Expr* f, Expr* x) {
    Expr* d = rt_eval2("D", expr_copy(result), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    /* Fast EXACT zero-test before the (occasionally hang-prone) Simplify: TrigToExp
     * collapses the circular/hyperbolic kernels to exponentials, and Together over
     * that single rational function of E^(...) reduces a genuinely-zero difference
     * to 0 at once — a sound certificate.  This makes the verification of a
     * trig/hyperbolic antiderivative decidable in O(rational).  NB: TrigToExp must
     * come FIRST — a plain Together on the raw hyperbolic form itself blows up doing
     * GCD over the opaque Cosh/Sinh kernels (the Sech[x] / Sec[x]^3 hang):
     * Simplify[D[2 ArcTan[Cosh x + Sinh x]] - Sech x] and Together of it both hang,
     * yet Together[TrigToExp[...]] = 0 immediately. */
    Expr* t2 = rt_eval1("Together", rt_eval1("TrigToExp", expr_copy(diff)));
    if (t2 && rt_is_zero(t2)) { expr_free(t2); expr_free(diff); return true; }
    if (t2) expr_free(t2);
    /* Exact single-kernel rational-trig decision (Together over the raw Laurent
     * TrigToExp form above leaves E^(-i x) poles uncombined for a Log/ArcTan-of-Tan
     * diff-back — e.g. the §5.10 hypertangent answer for a rational function of
     * Tan[x] — so it reads non-zero; and Simplify then mis-reduces the Sqrt-wrapped
     * form to a spurious residual).  trigexp_rational_is_zero clears the kernel by
     * an EXACT rational point-grid (a proof, never a numeric sample), so a genuine
     * rational-in-Tan antiderivative certifies here in sub-ms.  Rigorous both ways;
     * declines (UNKNOWN) to Simplify on multi-kernel / opaque-log forms. */
    TrigExpZeroResult tez = trigexp_rational_is_zero(diff);
    if (tez == TRIGEXP_ZERO_TRUE) { expr_free(diff); return true; }
    /* Residual (algebraic / log) forms the fast test leaves unreduced fall back to
     * Simplify (unchanged from the original gate). */
    Expr* s = rt_eval1("Simplify", diff);
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* A single-extension case is valid only when its kernel's defining function `u`
 * (a Log argument or an exp exponent) is a rational function of x ALONE — i.e.
 * free of any other exp/log of x.  A NESTED kernel (e.g. u = E^x for the outer
 * kernel of E^(E^x)) is a two-extension tower: the single-kernel derivation
 * Dt = u' theta would carry the unsubstituted inner kernel, and SolveAlways would
 * treat it as a free parameter and certify a WRONG residue.  Such integrands must
 * be left to the tower cases (rt_log_tower_case / rt_exp_tower_case). */
bool rt_kernel_simple(Expr* u, Expr* x) {
    return rt_find_exp_of_x(u, x) == NULL && rt_find_log_of_x(u, x) == NULL;
}
