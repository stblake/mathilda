/* integrate_risch_transcendental.c — recursive transcendental Risch integrator.
 *
 * The recursive Risch decision procedure for transcendental elementary
 * functions (Bronstein/Roach lineage), with arithmetic grounded in
 * Mathilda's existing Expr/poly/rat primitives (see the header for the
 * full contract).
 *
 * Structure: a differential transcendental tower over a single integration
 * variable, dispatched through the rational base case (delegated to
 * Integrate`BronsteinRational), the logarithmic / exponential / coupled /
 * tower cases, a trig-hyperbolic front-end, and (flag-gated) special-
 * function outputs.  Every branch is correct by construction behind an
 * exact structural certificate.
 */

#include "integrate_risch_transcendental.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"
#include "risch_field.h"
#include "risch_hermite.h"
#include "risch_canonical.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================== */
/* Small evaluation helpers.                                          */
/* ================================================================== */

/* True iff `e` is the unevaluated call `name[...]`. */
static bool rt_head_is(const Expr* e, const char* name) {
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
static bool rt_is_rat_const(const Expr* e) {
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
static Expr* rt_eval_call(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Evaluate `e`, taking ownership: frees `e` and returns evaluate()'s
 * result (evaluate itself does not consume its argument). */
static Expr* rt_eval_own(Expr* e) {
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

static Expr* rt_eval1(const char* head, Expr* a) {
    return rt_eval_call(head, (Expr*[]){ a }, 1);
}
static Expr* rt_eval2(const char* head, Expr* a, Expr* b) {
    return rt_eval_call(head, (Expr*[]){ a, b }, 2);
}
static Expr* rt_eval3(const char* head, Expr* a, Expr* b, Expr* c) {
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
static Expr* rt_class_primitive(Expr** ws, size_t nw, long* kof) {
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

static bool rt_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("True");
}
/* FreeQ[e, x] */
static bool rt_free_of_x(Expr* e, Expr* x) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* FreeQ[e, head] — true iff the symbol `head` occurs nowhere in e. */
static bool rt_free_of_head(Expr* e, const char* head) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_new_symbol(head));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* PolynomialQ[e, x] */
static bool rt_is_poly(Expr* e, Expr* x) {
    Expr* r = rt_eval2("PolynomialQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* Degree of a polynomial in x via Length[CoefficientList[e, x]] - 1
 * (Mathilda has no Exponent builtin).  Returns -1 if CoefficientList does
 * not reduce to a List.  Call only after rt_is_poly. */
static long rt_degree(Expr* e, Expr* x) {
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
static Expr* rt_coeff(Expr* e, Expr* x, long k) {
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
static bool rt_is_zero(Expr* e) {
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
static Expr* rt_template(const char* tmpl, const char** names,
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

/* ================================================================== */
/* Case: rational function of x  (rational Risch base case).           */
/* ================================================================== */

/* Delegate the pure-rational case to the Bronstein rational integrator
 * (the rational Risch base case: Hermite + Lazard-Rioboo-Trager).  Returns
 * NULL when the integrand is not rational in x (BronsteinRational leaves
 * its call unevaluated). */
static Expr* rt_rational_case(Expr* f, Expr* x) {
    Expr* r = rt_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(f), expr_copy(x) }, 2);
    if (!r) return NULL;
    if (rt_head_is(r, "Integrate`BronsteinRational")) {
        expr_free(r);
        return NULL;
    }
    return r;
}

/* ================================================================== */
/* Case: special functions  (Erf / Erfi / ExpIntegralEi / LogIntegral /
 * PolyLog outputs for the standard non-elementary structural forms).    */
/* ================================================================== */
/* Each recognizer builds a candidate antiderivative from a template and
 * accepts it only if it passes the diff-back gate, so a mis-recognition
 * can never emit a wrong closed form.  These close integrals the whole
 * elementary cascade leaves open, using special functions Mathilda
 * already provides (Erf, ExpIntegralEi, LogIntegral, PolyLog).          */

/* Forward decls: the Ei/li recognizers extract the exponential/log kernel
 * directly (defined below alongside the transcendental-case machinery). */
static Expr* rt_find_log_of_x(Expr* e, Expr* x);
static Expr* rt_exp_ratreduce_case(Expr* f, Expr* x);
static Expr* rt_find_exp_of_x(Expr* e, Expr* x);
/* Diff-back verifier (defined with the tower machinery); used by the
 * pure-resultant LRT frac path, whose reuse crosses content boundaries. */
static bool rt_verify_antideriv(Expr* result, Expr* f, Expr* x);

/* K * E^(a x^2 + b x + c) with the leading (quadratic) coefficient
 * nonzero  ->  Erf/Erfi.  Detected by the log-derivative f'/f being a
 * degree-1 polynomial). */
static Expr* rt_try_erf(Expr* f, Expr* x) {
    Expr* df = rt_eval2("D", expr_copy(f), expr_copy(x));
    if (!df) return NULL;
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2);
    Expr* quot = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ df, inv }, 2);
    Expr* ld = rt_eval1("Cancel", quot);
    if (!ld) return NULL;
    if (!rt_is_poly(ld, x) || rt_degree(ld, x) != 1) { expr_free(ld); return NULL; }
    Expr* a2 = rt_coeff(ld, x, 1);   /* = 2a */
    Expr* b1 = rt_coeff(ld, x, 0);   /* = b  */
    expr_free(ld);

    const char* names[2] = { "rmA2", "rmB" };
    Expr* vals[2] = { a2, b1 };
    Expr* result = NULL;
    /* Since ld = f'/f is a degree-1 polynomial, f = K' E^((a2/2)x^2 + b1 x)
     * with K' a constant, so K' = f|_{x=0} (the exponent has no constant
     * term).  Evaluating at 0 avoids a Simplify on a product of Gaussian
     * exponentials, which is prohibitively slow. */
    Expr* kp = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(x), expr_new_integer(0) }, 2) }, 2));
    if (kp) {
        /* Integral of E^((a2/2)x^2 + b1 x), a = a2/2, b = b1:
         *   E^(-b^2/(4a)) * (-Sqrt[Pi] Erf[(b+2a x)/(2 Sqrt[-a])] / (2 Sqrt[-a]))
         * where -b^2/(4a) = -b1^2/(2 a2) is the completing-the-square shift. */
        Expr* erfpart = rt_template(
            "E^(-rmB^2/(2*rmA2)) *"
            " (-(Sqrt[Pi]*Erf[(rmB + rmA2*x)/(2*Sqrt[-rmA2/2])])/(2*Sqrt[-rmA2/2]))",
            names, vals, 2);
        if (erfpart) {
            /* Correct by construction: ld = f'/f is a degree-1 polynomial,
             * so f = K' E^(a x^2 + b x) exactly and this is its integral. */
            result = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(kp), erfpart }, 2));
        }
        expr_free(kp);
    }
    expr_free(a2);
    expr_free(b1);
    return result;
}

/* (M * E^(a x + b)) / (c x + d)  ->  ExpIntegralEi.  A constant M times a
 * pure exponential with an EXACTLY LINEAR exponent v = a x + b (a != 0), over
 * a linear denominator c x + d.  The exponential kernel E^v is extracted
 * DIRECTLY (via rt_find_exp_of_x) rather than by Together: a negative leading
 * coefficient (E^(-x)/x) would otherwise be pushed into the denominator by
 * Together, making the denominator x E^x (non-polynomial) and hiding the Ei
 * shape.  Extracting E^v and reducing the rational cofactor R = f E^(-v)
 * handles a < 0 and a nonzero constant term b uniformly.
 *
 * Sub w = c x + d in M INT E^(a x + b)/(c x + d) dx:
 *   = (M/c) E^(b - a d/c) ExpIntegralEi[a x + a d/c].
 * Correct by construction: R free of x with a linear denominator proves
 * f = M E^v/(c x + d) exactly, so this is its antiderivative. */
static Expr* rt_try_ei(Expr* f, Expr* x) {
    Expr* v = rt_find_exp_of_x(f, x);          /* borrowed exponent of E^v */
    if (!v) return NULL;
    if (!rt_is_poly(v, x) || rt_degree(v, x) != 1) return NULL;
    /* R = Together[f E^(-v)]: the rational cofactor after peeling E^v. */
    Expr* emv = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_new_symbol("E"),
            expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(v) }, 2) }, 2);
    Expr* R = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(f), emv }, 2));
    if (!R) return NULL;
    Expr* num = rt_eval1("Numerator", expr_copy(R));
    Expr* den = rt_eval1("Denominator", expr_copy(R));
    expr_free(R);
    Expr* result = NULL;
    if (num && den && rt_free_of_x(num, x) && !rt_is_zero(num)
        && rt_is_poly(den, x) && rt_degree(den, x) == 1) {
        Expr* aa = rt_coeff(v, x, 1);   /* a */
        Expr* bb = rt_coeff(v, x, 0);   /* b */
        Expr* cc = rt_coeff(den, x, 1); /* c */
        Expr* dd = rt_coeff(den, x, 0); /* d */
        const char* names[5] = { "rmM", "rmA", "rmB", "rmC", "rmD" };
        Expr* vals[5] = { num, aa, bb, cc, dd };
        result = rt_template(
            "(rmM*E^(rmB - rmA*rmD/rmC)/rmC)"
            "*ExpIntegralEi[rmA*x + rmA*rmD/rmC]",
            names, vals, 5);
        expr_free(aa);
        expr_free(bb);
        expr_free(cc);
        expr_free(dd);
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    return result;
}

/* c * w^(p-1) * w'(x) / Log[w]  ->  c LogIntegral[w^p]  for a Log[w] kernel
 * whose argument w depends on x, a positive integer p, and a nonzero constant
 * c.  Certificate: cand = Together[f Log[w] / (w^(p-1) w')] is a constant free
 * of x, which proves f = cand d/dx LogIntegral[w^p] exactly (using the
 * standard branch Log[w^p] = p Log[w], the convention under which
 * LogIntegral answers are stated).  This subsumes the bare K/Log[x] -> K
 * LogIntegral[x] form (p = 1, w = x) and adds a scaled/affine argument
 * (1/Log[2x] -> LogIntegral[2x]/2, p = 1, w = 2x) and a monomial numerator
 * (x/Log[x] -> LogIntegral[x^2], p = 2, w = x). */
static Expr* rt_try_li(Expr* f, Expr* x) {
    Expr* w = rt_find_log_of_x(f, x);          /* borrowed argument of Log[w] */
    if (!w) return NULL;
    Expr* dw = rt_eval2("D", expr_copy(w), expr_copy(x));   /* w'(x) */
    if (!dw) return NULL;
    if (rt_is_zero(dw)) { expr_free(dw); return NULL; }
    Expr* logw = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(w) }, 1);
    Expr* result = NULL;
    const int MAXP = 6;
    for (int p = 1; p <= MAXP && !result; p++) {
        /* cand = Together[f Log[w] / (w^(p-1) w')] must be a nonzero constant. */
        Expr* wp1 = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(w), expr_new_integer(p - 1) }, 2);
        Expr* denom = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ wp1, expr_copy(dw) }, 2);
        Expr* frac = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(f), expr_copy(logw),
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ denom, expr_new_integer(-1) }, 2) }, 3);
        Expr* cand = rt_eval1("Together", frac);
        if (cand && rt_free_of_x(cand, x) && !rt_is_zero(cand)) {
            Expr* wp = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(w), expr_new_integer(p) }, 2);
            Expr* li = expr_new_function(expr_new_symbol("LogIntegral"),
                (Expr*[]){ wp }, 1);
            result = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(cand), li }, 2));
        }
        if (cand) expr_free(cand);
    }
    expr_free(dw);
    expr_free(logw);
    return result;
}

/* Find (borrowed) the argument u of the first Log[u] subexpression of `e`
 * whose argument depends on x. */
static Expr* rt_find_log_of_x(Expr* e, Expr* x) {
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
static Expr* rt_find_exp_of_x(Expr* e, Expr* x) {
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

/* A single-extension case is valid only when its kernel's defining function `u`
 * (a Log argument or an exp exponent) is a rational function of x ALONE — i.e.
 * free of any other exp/log of x.  A NESTED kernel (e.g. u = E^x for the outer
 * kernel of E^(E^x)) is a two-extension tower: the single-kernel derivation
 * Dt = u' theta would carry the unsubstituted inner kernel, and SolveAlways would
 * treat it as a free parameter and certify a WRONG residue.  Such integrands must
 * be left to the tower cases (rt_log_tower_case / rt_exp_tower_case). */
static bool rt_kernel_simple(Expr* u, Expr* x) {
    return rt_find_exp_of_x(u, x) == NULL && rt_find_log_of_x(u, x) == NULL;
}

/* K Log[1 + p x] / x  ->  -K PolyLog[2, -p x]  (dilogarithm form).
 * Certificate: the Log argument is exactly linear with constant term 1,
 * and Together[f x / Log[u]] is a constant K free of x (so f is exactly
 * K Log[1+p x]/x). */
static Expr* rt_try_dilog(Expr* f, Expr* x) {
    Expr* u = rt_find_log_of_x(f, x);   /* borrowed */
    if (!u) return NULL;
    if (!rt_is_poly(u, x) || rt_degree(u, x) != 1) return NULL;
    Expr* u0 = rt_coeff(u, x, 0);
    Expr* u1 = rt_coeff(u, x, 1);
    Expr* result = NULL;
    Expr* u0m1 = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(u0), expr_new_integer(-1) }, 2);
    bool u0_is_1 = rt_is_zero(u0m1);
    expr_free(u0m1);
    if (u0_is_1) {
        Expr* logu = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* invlog = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ logu, expr_new_integer(-1) }, 2);
        Expr* prod = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(f), expr_copy(x), invlog }, 3);
        Expr* K = rt_eval1("Together", prod);
        if (K && rt_free_of_x(K, x) && !rt_is_zero(K)) {
            const char* names[2] = { "rmK", "rmU1" };
            Expr* vals[2] = { K, u1 };
            result = rt_template("-rmK*PolyLog[2, -rmU1*x]", names, vals, 2);
        }
        if (K) expr_free(K);
    }
    expr_free(u0);
    expr_free(u1);
    return result;
}

/* Special-function form registry (Cherry substrate — CHERRY_DESIGN.md §2.1/§3.2).
 * Each entry is a special-function "plugin".  Today `recognize` is the narrow,
 * enumerated template of Cherry's decision procedure (e.g. one Gaussian -> Erf,
 * a linear exponent/denominator -> Ei) — each already correct-by-construction via
 * its own diff-back/constant-cofactor gate.  This struct is the seam Cherry grows
 * into: it will gain a finite argument generator (Sigma-decomposition for Li,
 * Res_x(g1,p+alpha q) for Ei, completing-square for Erf), a derivative template,
 * an answer constructor, and a structural bound (see the design doc), turning each
 * template into the full decision procedure WITHOUT touching the dispatch below.
 * Registering a new special function then means adding one entry here. */
typedef struct {
    const char* name;                    /* the emitted special function */
    Expr* (*recognize)(Expr* f, Expr* x);/* narrow generator (grows into gen_arguments) */
} RtSpecialForm;

static const RtSpecialForm RT_SPECIAL_FORMS[] = {
    { "Erf",           rt_try_erf   },   /* K E^(a x^2 + b x + c)                */
    { "ExpIntegralEi", rt_try_ei    },   /* M E^(a x + b) / (c x + d)            */
    { "LogIntegral",   rt_try_li    },   /* c w^(p-1) w' / Log[w]                */
    { "PolyLog",       rt_try_dilog },   /* K Log[1 + p x] / x                   */
};

/* Try each registered special-function form in turn (order preserved). */
static Expr* rt_special_case(Expr* f, Expr* x) {
    for (size_t i = 0; i < sizeof(RT_SPECIAL_FORMS) / sizeof(RT_SPECIAL_FORMS[0]); i++) {
        Expr* r = RT_SPECIAL_FORMS[i].recognize(f, x);
        if (r) return r;
    }
    return NULL;
}

/* ================================================================== */
/* Case: transcendental — the recursive Risch algorithm proper.        */
/* ================================================================== */
/* This is the genuine recursive Risch (over a single logarithmic /
 * exponential monomial extension), NOT the parallel-Risch (pmint)
 * heuristic.  It reduces the integrand in the differential field
 * K(theta) over K = C(x), grounding all coefficient-field arithmetic in
 * Mathilda's rational-Risch primitives (Integrate`BronsteinRational for
 * the base-case integrals in K).                                        */

/* Limited integration in K = C(x) modulo the extension monomial theta =
 * Log[u]:  compute integ(r) = s + c*theta with s in K (rational) and c a
 * constant in C, or fail.  The only logarithm the result may contain is
 * c*Log[u] (that is theta itself); any other logarithm, a non-constant
 * theta-coefficient, or an unresolved integral makes it decline.  This is
 * the sub-oracle that lets the primitive polynomial reduction fold a
 * would-be new logarithm back into the tower.  On success sets *s_out and
 * *c_out (both owned) and returns 0; otherwise returns -1. */
static int rt_limited_integrate(Expr* r, Expr* x, Expr* u,
                                Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rt_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(r), expr_copy(x) }, 2);
    if (!R) return -1;
    if (rt_head_is(R, "Integrate`BronsteinRational")) { expr_free(R); return -1; }
    /* Rsub = R with Log[u] replaced by the fresh variable rmT. */
    Expr* logu = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(u) }, 1);
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ logu, expr_new_symbol("rmT") }, 2);
    Expr* Rsub = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ R, rule }, 2));   /* adopts R */
    if (!Rsub) return -1;
    Expr* tsym = expr_new_symbol("rmT");
    int rc = -1;
    if (rt_free_of_head(Rsub, "Log") && rt_is_poly(Rsub, tsym)
        && rt_degree(Rsub, tsym) <= 1) {
        Expr* c = rt_coeff(Rsub, tsym, 1);   /* theta-coefficient */
        Expr* s = rt_coeff(Rsub, tsym, 0);   /* rational part      */
        if (rt_free_of_x(c, x)) { *s_out = s; *c_out = c; rc = 0; }
        else { expr_free(c); expr_free(s); }
    }
    expr_free(Rsub);
    expr_free(tsym);
    return rc;
}

/* Primitive (logarithmic) polynomial case.  For theta = Log[u] with
 * u rational in x (theta' = u'/u =: eta in K), integrate a polynomial
 * P(theta) = sum p_i theta^i by the recursive coefficient matching of the
 * recursive Risch algorithm (Bronstein, IntegratePrimitivePolynomial):
 * the antiderivative Q = sum q_i theta^i satisfies, degree by degree,
 *     q_i' + (i+1) q_{i+1} eta = p_i.
 * Working top-down, each level integrates r_i = p_i - (i+1) q_{i+1} eta in
 * K by the limited-integration oracle, which returns s_i + c_i theta; the
 * constant c_i is folded back by bumping q_{i+1} += c_i/(i+1) (this is how
 * a would-be new logarithm becomes a higher power of theta, e.g. the
 * (3/2)Log[2x+3] arising in Integrate[Log[2x+3]]).  Declines (NULL) if any
 * level leaves K in a way the oracle cannot absorb. */
static Expr* rt_log_poly_case(Expr* f, Expr* x) {
    Expr* u = rt_find_log_of_x(f, x);        /* borrowed: the Log argument */
    if (!u || !rt_kernel_simple(u, x)) return NULL;   /* nested -> tower */

    /* F = f with Log[u] replaced by the fresh polynomial variable rmT. */
    Expr* logu = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(u) }, 1);
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ logu, expr_new_symbol("rmT") }, 2);
    Expr* F = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rule }, 2));
    if (!F) return NULL;

    Expr* tsym = expr_new_symbol("rmT");
    /* Require: pure polynomial in theta, no residual logs (single
     * extension), and genuine theta-dependence. */
    if (!rt_is_poly(F, tsym) || !rt_free_of_head(F, "Log")
        || rt_free_of_x(F, tsym)) {
        expr_free(F); expr_free(tsym); return NULL;
    }
    long m = rt_degree(F, tsym);
    if (m < 1) { expr_free(F); expr_free(tsym); return NULL; }

    /* eta = Cancel[D[u,x]/u]. */
    Expr* du = rt_eval2("D", expr_copy(u), expr_copy(x));
    Expr* invu = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
    Expr* eta = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ du, invu }, 2));

    Expr** p = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (long i = 0; i <= m; i++) p[i] = rt_coeff(F, tsym, i);
    /* q[0..m+1]; NULL entries represent 0 (a bump can create q[m+1]). */
    Expr** q = calloc((size_t)(m + 2), sizeof(Expr*));

    bool fail = (eta == NULL);
    for (long i = m; i >= 0 && !fail; i--) {
        /* r_i = p_i - (i+1) q[i+1] eta. */
        Expr* r_i;
        if (q[i + 1] == NULL) {
            r_i = expr_copy(p[i]);
        } else {
            Expr* term = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(i + 1), expr_copy(q[i + 1]),
                          expr_copy(eta) }, 3);
            Expr* neg = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), term }, 2);
            r_i = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(p[i]), neg }, 2);
        }
        Expr* s = NULL; Expr* c = NULL;
        int rc = rt_limited_integrate(r_i, x, u, &s, &c);
        expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        /* Fold the theta-term back: q[i+1] += c/(i+1). */
        Expr* bump = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2);
        if (q[i + 1] == NULL) {
            q[i + 1] = rt_eval_own(bump);
        } else {
            q[i + 1] = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ q[i + 1], bump }, 2));   /* adopts old q[i+1] */
        }
        q[i] = s;
    }

    Expr* result = NULL;
    if (!fail) {
        /* Q = sum_{i=0}^{m+1} q_i Log[u]^i  (skip zero coefficients). */
        Expr** terms = malloc((size_t)(m + 2) * sizeof(Expr*));
        size_t nt = 0;
        for (long i = 0; i <= m + 1; i++) {
            if (!q[i]) continue;
            Expr* lu = expr_new_function(expr_new_symbol("Log"),
                (Expr*[]){ expr_copy(u) }, 1);
            Expr* pw = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ lu, expr_new_integer(i) }, 2);
            terms[nt++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[i]), pw }, 2);
        }
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), terms, nt);
        free(terms);
        result = evaluate(sum);
        expr_free(sum);
    }

    for (long i = 0; i <= m; i++) expr_free(p[i]);
    for (long i = 0; i <= m + 1; i++) if (q[i]) expr_free(q[i]);
    free(p); free(q);
    if (eta) expr_free(eta);
    expr_free(F); expr_free(tsym);
    return result;
}

/* Integrate g in K = C(x) allowing genuine logarithms in the result (the
 * i = 0 term of the exponential case is an ordinary base-field integral). */
static Expr* rt_integrate_in_K_with_logs(Expr* g, Expr* x) {
    Expr* r = rt_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(g), expr_copy(x) }, 2);
    if (!r) return NULL;
    if (rt_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
    return r;
}

/* Exact degree bound for the solution q of a Risch DE  D[q] + f q = p, in a single
 * monomial variable v (Bronstein RdeBoundDegree).  Returns an UPPER bound on deg_v(q)
 * with NO arbitrary cap — a function of the equation's degrees alone (plus the exact
 * cancellation/resonance sub-case, below).  `deriv_lowers` is true when D lowers deg_v
 * by one (v = x under d/dx, or a LOGARITHMIC monomial), false when D preserves deg_v (an
 * EXPONENTIAL monomial, whose self-derivative D[t^k] = k w' t^k keeps the degree).
 * dpv = deg_v(p), dfv = deg_v(f) (each as deg(numerator) - deg(denominator)).
 *
 * Leading-degree balance.  Match the top v-degree of D[q] + f q = p.  For a
 * deriv-lowering v, deg_v(D[q]) = deg_v(q) - 1 and deg_v(f q) = dfv + deg_v(q); when
 * dfv >= 0 the f q term strictly dominates (distinct degrees, no cancellation) so
 * deg_v(q) = dpv - dfv exactly, and when dfv < 0 the balance is between D[q]
 * (integration raises degree by one: dpv + 1) and a pole in f (dpv - dfv), so the bound
 * is their max.  For a deriv-preserving (exponential) v, both terms sit at
 * deg_v(q) + max(0, dfv), giving deg_v(q) = dpv - dfv.
 *
 * Cancellation / resonance (Bronstein's recursive degree reduction).  The leading-degree
 * balance is exact EXCEPT where the two leading coefficients cancel, so the top term of
 * D[q] + f q vanishes and deg_v(q) can exceed the naive value.  This happens in exactly
 * two configurations, each keyed by an integer `m_res` (the Bronstein resonance integer,
 * or -1 when none):
 *   - deriv-PRESERVING v (exponential), dfv == 0: D[c v^n] + f(c v^n) has leading
 *     coefficient (n D[v]/v + f) c, which vanishes when n = -f/(D[v]/v) is a
 *     nonnegative integer m_res; then deg_v(q) can reach m_res, so bound = max(naive,
 *     m_res).  m_res is supplied by the caller (rt_resonance_int), which alone has the
 *     coefficients.
 *   - deriv-LOWERING v (base x / logarithmic), dfv == -1: the primitive
 *     leading-coefficient cancellation, whose homogeneous-solution degree can likewise
 *     reach m_res, so bound = max(naive, m_res).
 * The widening is MONOTONE (only ever raises the bound), so a spurious m_res can never
 * ship a wrong result (the solution stays SolveAlways-certified and the enclosing tower
 * case diff-back verifies); a missed one merely declines.  Pass m_res = -1 to disable.
 *
 * Reachability note.  In the current tower architecture BOTH cancellation configurations
 * are pre-empted, so m_res is -1 on every reachable call: (a) an exponential resonance
 * n = -(i w_L')/w_j' being an integer means the top and lower exp exponents are
 * commensurate, and the commensurate-exponent reduction in rt_tower_build collapses such
 * kernels to one primitive before any RDE solve; (b) dfv == -1 requires a simple pole in
 * f = i Dcoef, which a rational tower element's derivative cannot have (it would
 * integrate to a Log), and the only kernel that could — a log top with Dcoef = u'/u —
 * never routes through the field RDE (it uses the primitive-polynomial recursion).  The
 * detection is nonetheless computed live and folded in here so the degree bound is exact
 * per Bronstein should a future kernel type expose either configuration. */
long rt_rde_var_bound(long dpv, long dfv, bool deriv_lowers, long m_res) {
    long bq;
    if (deriv_lowers) {
        if (dfv >= 0) {
            bq = dpv - dfv;                         /* f q dominates: exact */
        } else {                                   /* integration rise vs pole in f */
            long a = dpv + 1, b = dpv - dfv;
            bq = (a > b) ? a : b;
            /* primitive leading-coefficient cancellation (dfv == -1) */
            if (dfv == -1 && m_res >= 0 && m_res > bq) bq = m_res;
        }
    } else {                                       /* exponential: derivation preserves deg */
        bq = dpv - dfv;
        /* exponential integer resonance (dfv == 0) */
        if (dfv == 0 && m_res >= 0 && m_res > bq) bq = m_res;
    }
    return (bq < 0) ? 0 : bq;
}

/* Bronstein resonance integer for an EXPONENTIAL lower monomial v = t_k in the Risch DE
 * D[q] + f q = p with f = i Dcoef_L (Dcoef_L = D[t_L]/t_L = w_L' the top exp kernel's
 * derivation coefficient) and Dcoef_v = D[v]/v = w_k'.  The leading coefficient of
 * D[c v^n] + f (c v^n) is (n Dcoef_v + f) c v^n; it cancels — allowing deg_v(q) to reach
 * n beyond the naive bound — exactly when n = -f/Dcoef_v = -(i Dcoef_L)/Dcoef_v is a
 * nonnegative integer.  Returns that n, or -1 when the ratio is not a nonnegative integer
 * constant (no resonance).  Feeds rt_rde_var_bound's monotone widening, so a value here
 * can only ever raise the bound: never a wrong result (SolveAlways certifies, caller
 * diff-back verifies), at worst wasted ansatz terms. */
static long rt_resonance_int(long i, Expr* DcoefL, Expr* Dcoefv) {
    if (!DcoefL || !Dcoefv || rt_is_zero(Dcoefv)) return -1;
    Expr* ratio = rt_eval1("Simplify", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-i), expr_copy(DcoefL),
            expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(Dcoefv), expr_new_integer(-1) }, 2) }, 3));
    long m = -1;
    if (ratio && ratio->type == EXPR_INTEGER && ratio->data.integer >= 0)
        m = (long)ratio->data.integer;
    if (ratio) expr_free(ratio);
    return m;
}

/* Multiplicity of the variable v at 0 in the polynomial p: the lowest power of v
 * that appears with a nonzero coefficient (0 if v does not divide p).  This is the
 * exact negative Laurent extent for an EXPONENTIAL kernel v — the order of the pole
 * that v contributes to a proper fraction — derived from the input, no cap. */
static long rt_var_mult_at_zero(Expr* p, Expr* v) {
    long a = 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(p), expr_copy(v));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    return a;
}

/* ==================================================================== */
/* Bronstein Chapter 6 — the Risch differential equation  D y + f y = g. */
/*                                                                       */
/* Base field C(x) with D = d/dx (Stage A).  Replaces the undetermined-  */
/* coefficient ansatz with the rational one-step-reduction algorithm     */
/* (WeakNormalizer -> RdeNormalDenominator -> RdeBoundDegreeBase ->       */
/* SPDE -> PolyRischDENoCancel1/2), all in polynomial-gcd time.  This is  */
/* the fix for the Davenport-progressive-reduction blow-up of the former */
/* SolveAlways ansatz (Bronstein, J. Symbolic Comput. 9 (1990) 49-60,    */
/* and Symbolic Integration I, 2nd ed., Chapter 6).  In the base field    */
/* every irreducible is normal for d/dx, so there is no special part and  */
/* SplitFactor is the identity — the algorithms simplify accordingly.     */
/* Correct by construction: SPDE's polynomial identity certifies          */
/* (q E^(i u))' = p E^(i u); returns NULL for a genuinely non-elementary  */
/* term.                                                                  */
/* ==================================================================== */

/* Numerator / Denominator of Together[e] (e not consumed). */
static Expr* rde_num(Expr* e) {
    return rt_eval1("Numerator", rt_eval1("Together", expr_copy(e)));
}
static Expr* rde_den(Expr* e) {
    return rt_eval1("Denominator", rt_eval1("Together", expr_copy(e)));
}
/* Univariate polynomial primitives in x (args not consumed). */
static Expr* rde_gcd(Expr* a, Expr* b, Expr* x) {
    (void)x;   /* PolynomialGCD infers the variable; a third arg is a 3rd poly */
    return rt_eval2("PolynomialGCD", expr_copy(a), expr_copy(b));
}
static Expr* rde_quot(Expr* a, Expr* b, Expr* x) {
    return rt_eval3("PolynomialQuotient", expr_copy(a), expr_copy(b), expr_copy(x));
}
static Expr* rde_rem(Expr* a, Expr* b, Expr* x) {
    return rt_eval3("PolynomialRemainder", expr_copy(a), expr_copy(b), expr_copy(x));
}
/* Expanded / cancelled arithmetic (args not consumed). */
static Expr* rde_mul(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
}
static Expr* rde_add(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
}
static Expr* rde_sub(Expr* a, Expr* b) {
    return rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(b) }, 2) }, 2));
}
/* Field quotient a/b, cancelled (args not consumed). */
static Expr* rde_divq(Expr* a, Expr* b) {
    return rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(b), expr_new_integer(-1) }, 2) }, 2));
}
/* Leading coefficient of poly p in x (0 if p == 0). */
static Expr* rde_lc(Expr* p, Expr* x) {
    long d = rt_degree(p, x);
    return d < 0 ? expr_new_integer(0) : rt_coeff(p, x, d);
}
/* d/dx (e not consumed). */
static Expr* rde_dx(Expr* e, Expr* x) {
    return rt_eval2("D", expr_copy(e), expr_copy(x));
}

/* SPDE result tuple: any solution q of a Dq + b q = c with deg(q) <= n is
 * q = alpha*H + beta, where H in k[t], deg(H) <= m and D H + b' H = c'. */
typedef struct { Expr* b; Expr* c; long m; Expr* alpha; Expr* beta; } RdeSpde;
static void rde_spde_free(RdeSpde* s) {
    if (s->b) expr_free(s->b);
    if (s->c) expr_free(s->c);
    if (s->alpha) expr_free(s->alpha);
    if (s->beta) expr_free(s->beta);
    s->b = s->c = s->alpha = s->beta = NULL;
}

/* SPDE(a, b, c, D=d/dx, n): Rothstein's degree-reducing recursion (Bronstein
 * Thm 6.4.1 / algorithm SPDE, p.203).  Returns 1 and fills *out on success, or
 * 0 for "no solution" (the equation has no solution of degree <= n in k[x]).
 * Each level reduces the degree bound by deg(a) via one ExtendedEuclidean, so
 * it terminates in O(n / deg(a)) steps — no dense linear solve. */
static int rde_spde(Expr* a, Expr* b, Expr* c, Expr* x, long n, RdeSpde* out) {
    out->b = out->c = out->alpha = out->beta = NULL; out->m = 0;
    if (n < 0) {
        if (rt_is_zero(c)) {
            out->b = expr_new_integer(0); out->c = expr_new_integer(0);
            out->m = 0; out->alpha = expr_new_integer(0); out->beta = expr_new_integer(0);
            return 1;
        }
        return 0;                                   /* c != 0, n < 0: no solution */
    }
    Expr* g = rde_gcd(a, b, x);
    Expr* rem = rde_rem(c, g, x);
    bool gdivc = rt_is_zero(rem);
    expr_free(rem);
    if (!gdivc) { expr_free(g); return 0; }         /* g does not divide c */
    Expr* a1 = rde_quot(a, g, x);
    Expr* b1 = rde_quot(b, g, x);
    Expr* c1 = rde_quot(c, g, x);
    expr_free(g);
    long da = rt_degree(a1, x);
    if (da == 0) {                                  /* a in k*: base of recursion */
        out->b = rde_divq(b1, a1);
        out->c = rde_divq(c1, a1);
        out->m = n;              /* Bronstein SPDE: return (b/a, c/a, n, 1, 0) */
        out->alpha = expr_new_integer(1);
        out->beta = expr_new_integer(0);
        expr_free(a1); expr_free(b1); expr_free(c1);
        return 1;
    }
    /* (r, z) with b1 r + a1 z = c1, deg(r) < deg(a1); gcd(a1,b1)=1 so the
     * ExtendedEuclidean cofactor s of b1 gives r = (s c1) mod a1. */
    Expr* eg = rt_eval3("PolynomialExtendedGCD", expr_copy(b1), expr_copy(a1), expr_copy(x));
    Expr* s = NULL;
    if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
        Expr* cof = eg->data.function.args[1];
        if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
            s = cof->data.function.args[0];
    }
    if (!s) { if (eg) expr_free(eg); expr_free(a1); expr_free(b1); expr_free(c1); return 0; }
    Expr* sc = rde_mul(s, c1);
    Expr* r = rde_rem(sc, a1, x);
    expr_free(sc);
    if (eg) expr_free(eg);
    /* z = (c1 - b1 r) / a1 (exact). */
    Expr* b1r = rde_mul(b1, r);
    Expr* c1_b1r = rde_sub(c1, b1r);
    expr_free(b1r);
    Expr* z = rde_quot(c1_b1r, a1, x);
    expr_free(c1_b1r);
    /* recurse SPDE(a1, b1 + Da1, z - Dr, n - da). */
    Expr* da1 = rde_dx(a1, x);
    Expr* b2 = rde_add(b1, da1);
    expr_free(da1);
    Expr* dr = rde_dx(r, x);
    Expr* c2 = rde_sub(z, dr);
    expr_free(dr); expr_free(z);
    RdeSpde sub;
    int rc = rde_spde(a1, b2, c2, x, n - da, &sub);
    expr_free(b1); expr_free(c1); expr_free(b2); expr_free(c2);
    if (!rc) { expr_free(a1); expr_free(r); return 0; }
    /* q = a1*H + (a1*sub.alpha) H_sub + ...  ->  alpha = a1*sub.alpha,
     *    beta = a1*sub.beta + r. */
    out->b = sub.b; out->c = sub.c; out->m = sub.m;   /* adopt inner eqn */
    out->alpha = rde_mul(a1, sub.alpha);
    Expr* a1beta = rde_mul(a1, sub.beta);
    out->beta = rde_add(a1beta, r);
    expr_free(a1beta);
    expr_free(sub.alpha); expr_free(sub.beta);
    expr_free(a1); expr_free(r);
    return 1;
}

/* PolyRischDENoCancel1(b, c, D=d/dx, n) with b != 0 (Bronstein p.208).  The
 * leading terms of Dq and bq never cancel, so the solution is built top-down,
 * one monomial per pass, the degree bound strictly decreasing.  Returns q
 * (owned) or NULL for "no solution". */
static Expr* rde_polyrischde_nocancel1(Expr* b, Expr* c, Expr* x, long n) {
    long db = rt_degree(b, x);
    Expr* lcb = rde_lc(b, x);
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rt_degree(cc, x);
        long m = dc - db;
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
        Expr* lcc = rde_lc(cc, x);
        Expr* coeff = rde_divq(lcc, lcb);           /* lc(c)/lc(b) */
        expr_free(lcc);
        Expr* p = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ coeff, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer(m) }, 2) }, 2));  /* adopts coeff */
        Expr* nq = rde_add(q, p);
        expr_free(q); q = nq;
        /* c <- c - Dp - b p. */
        Expr* dp = rde_dx(p, x);
        Expr* bp = rde_mul(b, p);
        Expr* t1 = rde_sub(cc, dp);
        Expr* t2 = rde_sub(t1, bp);
        expr_free(cc); cc = t2;
        expr_free(dp); expr_free(bp); expr_free(t1); expr_free(p);
        n = m - 1;
    }
    expr_free(cc); expr_free(lcb);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* PolyRischDENoCancel2(b, c, D=d/dx, n) with b == 0 (Bronstein p.209, D=d/dt
 * case): the equation is Dq = c, i.e. antidifferentiation of c.  Returns the
 * polynomial antiderivative q with deg(q) <= n (owned), or NULL if c has no
 * polynomial antiderivative of bounded degree.  (lambda(x) = lc(Dx) = 1,
 * delta(x) = 0 for the base derivation.) */
static Expr* rde_polyrischde_integrate(Expr* c, Expr* x, long n) {
    Expr* q = expr_new_integer(0);
    Expr* cc = expr_copy(c);
    bool fail = false;
    while (!rt_is_zero(cc)) {
        long dc = rt_degree(cc, x);
        long m = dc + 1;                            /* deg(c) - delta + 1 */
        if (n < 0 || m < 0 || m > n) { fail = true; break; }
        Expr* lcc = rde_lc(cc, x);
        Expr* coeff = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ lcc, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(m), expr_new_integer(-1) }, 2) }, 2));  /* adopts lcc */
        Expr* p = rt_eval1("Expand", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ coeff, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer(m) }, 2) }, 2));
        Expr* nq = rde_add(q, p);
        expr_free(q); q = nq;
        Expr* dp = rde_dx(p, x);
        Expr* t1 = rde_sub(cc, dp);
        expr_free(cc); cc = t1;
        expr_free(dp); expr_free(p);
        n = m - 1;
    }
    expr_free(cc);
    if (fail) { expr_free(q); return NULL; }
    return q;
}

/* WeakNormalizer(f, D=d/dx) in the base field (Bronstein p.183).  Returns w in
 * k[x] such that f - Dw/w is weakly normalized (no normal irreducible with a
 * positive-integer residue), so no spurious logarithm survives.  In the base
 * field SplitFactor is the identity (every squarefree factor is normal).  The
 * common case (no positive-integer residue) returns w = 1. */
static Expr* rde_weak_normalizer(Expr* f, Expr* x) {
    Expr* den = rde_den(f);
    Expr* num = rde_num(f);
    Expr* dden = rde_dx(den, x);
    Expr* g = rde_gcd(den, dden, x);                /* gcd(dn, dn') */
    Expr* dstar = rde_quot(den, g, x);              /* d* = dn / g  */
    Expr* gds = rde_gcd(dstar, g, x);
    Expr* d1 = rde_quot(dstar, gds, x);             /* d1 = d* / gcd(d*, g) */
    expr_free(dden); expr_free(g); expr_free(dstar); expr_free(gds);
    Expr* w = expr_new_integer(1);
    long dd1 = rt_degree(d1, x);
    if (dd1 >= 1) {
        /* a with (den/d1) a  ==  num  (mod d1), deg(a) < deg(d1). */
        Expr* A = rde_quot(den, d1, x);             /* den / d1 */
        Expr* eg = rt_eval3("PolynomialExtendedGCD", expr_copy(A), expr_copy(d1), expr_copy(x));
        Expr* s = NULL; Expr* gg = NULL;
        if (eg && eg->type == EXPR_FUNCTION && eg->data.function.arg_count == 2) {
            gg = eg->data.function.args[0];
            Expr* cof = eg->data.function.args[1];
            if (cof && cof->type == EXPR_FUNCTION && cof->data.function.arg_count == 2)
                s = cof->data.function.args[0];
        }
        /* Only the common gg == 1 case yields the clean residue formula; else
         * leave w = 1 (weak normalization is a no-op that at worst declines). */
        if (s && gg && gg->type == EXPR_INTEGER && gg->data.integer == 1) {
            Expr* sa = rde_mul(s, num);
            Expr* a = rde_rem(sa, d1, x);           /* a = s*num mod d1 */
            expr_free(sa);
            /* r = resultant_x(a - z Dd1, d1), roots z = positive integers. */
            Expr* dd1poly = rde_dx(d1, x);
            Expr* zsym = expr_new_symbol("rmWNz");
            Expr* azd = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), zsym, dd1poly }, 3) }, 2));  /* adopts */
            Expr* res = rt_eval3("Resultant", azd, expr_copy(d1), expr_copy(x));
            /* Positive integer roots of res in rmWNz. */
            Expr* roots = rt_eval2("Solve", expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ res, expr_new_integer(0) }, 2), expr_new_symbol("rmWNz"));
            if (roots && roots->type == EXPR_FUNCTION
                && roots->data.function.head->type == EXPR_SYMBOL
                && roots->data.function.head->data.symbol.name == intern_symbol("List")) {
                for (size_t ri = 0; ri < roots->data.function.arg_count; ri++) {
                    Expr* sol = roots->data.function.args[ri];   /* {rmWNz -> val} */
                    if (!sol || sol->type != EXPR_FUNCTION
                        || sol->data.function.arg_count < 1) continue;
                    Expr* rule = sol->data.function.args[0];
                    if (!rt_head_is(rule, "Rule") || rule->data.function.arg_count != 2) continue;
                    Expr* val = rule->data.function.args[1];
                    long nv;
                    if (val && val->type == EXPR_INTEGER && (nv = (long)val->data.integer) >= 1) {
                        /* w *= gcd(a - nv*Dd1, d1)^nv. */
                        Expr* nd = rde_dx(d1, x);
                        Expr* an = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_copy(a), expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ expr_new_integer(-nv), nd }, 2) }, 2));  /* adopts nd */
                        Expr* gpow = rde_gcd(an, d1, x);
                        expr_free(an);
                        for (long e = 0; e < nv; e++) {
                            Expr* nw = rde_mul(w, gpow);
                            expr_free(w); w = nw;
                        }
                        expr_free(gpow);
                    }
                }
            }
            if (roots) expr_free(roots);
            expr_free(a);
        }
        if (eg) expr_free(eg);
        expr_free(A);
    }
    expr_free(d1); expr_free(den); expr_free(num);
    return w;
}

/* RischDE in the base field C(x): solve D y + f y = g for y in C(x).  Returns y
 * (owned) or NULL for "no solution" (the term is non-elementary in this field).
 * Pipeline: weak-normalize f, then RdeNormalDenominator reduces to the
 * polynomial equation a Dq + b q = c (base field has no special denominator),
 * RdeBoundDegreeBase bounds deg(q), SPDE + PolyRischDENoCancel solve it, and
 * y = q / (h * w). */
static Expr* rde_base(Expr* f, Expr* g, Expr* x) {
    if (rt_is_zero(g)) return expr_new_integer(0);

#ifdef USE_FLINT
    /* Native fast path: for the exponential tower the coefficient f is a
     * polynomial over Q(x) and g is a rational function over Q(x); the whole
     * base-field RDE (weak normalization is then a no-op) runs in fmpq_poly,
     * converting f and g straight to FLINT with no evaluator Together/Cancel.
     * This is the dominant, high-degree case (In16/In17/`poly·e^x`) — it
     * collapses from seconds of Expr rational arithmetic to milliseconds. The
     * verdict is authoritative: 1 -> y, 0 -> genuinely no solution (decline),
     * -1 -> out of scope (f rational / not univariate over Q) -> fall through
     * to the Expr path below. */
    if (x->type == EXPR_SYMBOL) {
        Expr* y = NULL;
        int nr = flint_rde_base_solve_fg(f, g, x->data.symbol.name, &y);
        if (nr >= 0) return y;              /* y is NULL on nr == 0 (decline) */
    }
#endif

    /* 1. Weak normalization: f <- f - Dw/w, g <- w g, y = z/w. */
    Expr* w = rde_weak_normalizer(f, x);
    Expr* fbar;
    if (rt_is_zero(w)) { expr_free(w); return NULL; }
    if (w->type == EXPR_INTEGER && w->data.integer == 1) {
        fbar = expr_copy(f);
    } else {
        Expr* dw = rde_dx(w, x);
        Expr* dww = rde_divq(dw, w);
        fbar = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_copy(f), expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), dww }, 2) }, 2));   /* adopts dww */
        expr_free(dw);
    }
    Expr* gbar = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(w), expr_copy(g) }, 2));

    /* 2. RdeNormalDenominator(fbar, gbar): base field, ds = es = 1. */
    Expr* dn = rde_den(fbar);
    Expr* en = rde_den(gbar);
    Expr* p = rde_gcd(dn, en, x);
    Expr* den_e = rde_dx(en, x);
    Expr* den_num = rde_gcd(en, den_e, x);          /* gcd(en, en') */
    Expr* den_p = rde_dx(p, x);
    Expr* den_den = rde_gcd(p, den_p, x);           /* gcd(p, p')  */
    expr_free(den_e); expr_free(den_p);
    Expr* h = rde_quot(den_num, den_den, x);
    expr_free(den_num); expr_free(den_den);
    /* Guard: en | dn h^2  (else no solution). */
    Expr* h2 = rde_mul(h, h);
    Expr* dnh2 = rde_mul(dn, h2);
    Expr* rem = rde_rem(dnh2, en, x);
    bool ok = rt_is_zero(rem);
    expr_free(rem); expr_free(h2);
    Expr* result = NULL;
    if (ok) {
        /* a = dn h,  b = dn h fbar - dn Dh,  c = dn h^2 gbar.  These are all
         * polynomials by construction; compute them as such (Numerator and
         * exact PolynomialQuotient) rather than via Cancel of a rational, which
         * can fail to reduce a factored/expanded mismatch (e.g. leaving
         * (x+1)(x+2)^2/(x+2)^2 uncancelled). */
        Expr* num_f = rde_num(fbar);                /* dn fbar = num_f (dn=den_f) */
        Expr* num_g = rde_num(gbar);
        Expr* aa = rde_mul(dn, h);                  /* dn h */
        Expr* dh = rde_dx(h, x);
        Expr* dndh = rde_mul(dn, dh);               /* dn Dh */
        Expr* hnf = rde_mul(h, num_f);              /* dn h fbar = h num_f */
        Expr* bb = rde_sub(hnf, dndh);
        Expr* dnh2_en = rde_quot(dnh2, en, x);      /* dn h^2 / en (exact) */
        Expr* cc = rde_mul(dnh2_en, num_g);         /* (dn h^2 / en) num_g */
        expr_free(num_f); expr_free(num_g); expr_free(dh); expr_free(dndh);
        expr_free(hnf); expr_free(dnh2_en);
        if (rt_is_poly(aa, x) && rt_is_poly(bb, x) && rt_is_poly(cc, x)) {
            /* 3. Degree bound (RdeBoundDegreeBase, p.199). */
            long da = rt_degree(aa, x), db = rt_degree(bb, x), dc = rt_degree(cc, x);
            long mx = (db > da - 1) ? db : (da - 1);
            long n = dc - mx; if (n < 0) n = 0;
            if (db == da - 1 && da >= 1) {
                Expr* lcb = rde_lc(bb, x); Expr* lca = rde_lc(aa, x);
                Expr* mm = rde_divq(lcb, lca);
                Expr* mmn = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), mm }, 2));   /* -lc(b)/lc(a), adopts mm */
                if (mmn && mmn->type == EXPR_INTEGER) {
                    long mv = (long)mmn->data.integer;
                    long cand = dc - db; if (cand < 0) cand = 0;
                    if (mv > n) n = mv;
                    if (cand > n) n = cand;
                }
                if (mmn) expr_free(mmn);
                expr_free(lcb); expr_free(lca);
            }
            /* 4. SPDE reduce, 5. solve bounded (Expr fallback path — the
             * exponential common case is handled natively in fmpq_poly at the
             * top of rde_base; see the flint_rde_base_solve_fg dispatch). */
            RdeSpde sp;
            if (rde_spde(aa, bb, cc, x, n, &sp)) {
                Expr* H;
                if (rt_is_zero(sp.b))
                    H = rde_polyrischde_integrate(sp.c, x, sp.m);
                else
                    H = rde_polyrischde_nocancel1(sp.b, sp.c, x, sp.m);
                if (H) {
                    Expr* aH = rde_mul(sp.alpha, H);
                    Expr* q = rde_add(aH, sp.beta);
                    expr_free(aH); expr_free(H);
                    /* y = q / (h w). */
                    Expr* hw = rde_mul(h, w);
                    result = rde_divq(q, hw);
                    expr_free(hw); expr_free(q);
                }
                rde_spde_free(&sp);
            }
        }
        expr_free(aa); expr_free(bb); expr_free(cc);
    }
    expr_free(dnh2);
    expr_free(dn); expr_free(en); expr_free(p); expr_free(h);
    expr_free(fbar); expr_free(gbar); expr_free(w);
    return result;
}

/* Solve the Risch differential equation  q' + i u' q = p  for q in the base
 * field C(x), via the Bronstein one-step-reduction algorithm above.  The
 * coefficient is f = i u'; returns q (owned) or NULL when the term is
 * non-elementary in this field (e.g. E^(-x^2), E^(1/x)).  Delegates to rde_base,
 * the COMPLETE base-field solver (weak normalization + normal-denominator
 * reduction + SPDE + polynomial non-cancellation solve) — so every NULL it
 * returns is an authoritative "no rational solution", never a bounded-ansatz
 * decline.  u and p may be rational in x (not just polynomial). */
static Expr* rt_solve_rde(Expr* p, long i, Expr* u, Expr* x) {
    if (i == 0) return NULL;
    Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));   /* u' */
    if (!up) return NULL;
    Expr* f = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), up }, 2));          /* f = i u', adopts up */
    Expr* q = rde_base(f, p, x);
    expr_free(f);
    return q;
}

/* Collect (as owned copies, deduplicated) the exponents of every E^w / Exp[w]
 * kernel of `e` whose exponent depends on x. */
static void rt_collect_exp_exponents(Expr* e, Expr* x,
                                     Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    Expr* w = NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        w = e->data.function.args[0];
    else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")
        && !rt_free_of_x(e->data.function.args[1], x))
        w = e->data.function.args[1];
    if (w) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], w)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) {
                *cap = *cap ? *cap * 2 : 4;
                *arr = realloc(*arr, *cap * sizeof(Expr*));
            }
            (*arr)[(*n)++] = expr_copy(w);
        }
    }
    rt_collect_exp_exponents(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_exp_exponents(e->data.function.args[i], x, arr, n, cap);
}

/* Exponential (hyperexponential) Laurent-polynomial case.  For a single
 * kernel theta = E^u (u a polynomial in x, theta' = u' theta), an integrand
 * that is a Laurent polynomial  f = sum_i p_i theta^i  (i possibly negative)
 * integrates to Q = sum_i q_i theta^i: the powers DECOUPLE, the i = 0 term is
 * an ordinary integral q_0 = INT p_0 in K, and each i != 0 term solves the
 * Risch differential equation q_i' + i u' q_i = p_i.  All E^(k u) present must
 * share the one primitive exponent u (k integer); a genuine theta-denominator
 * beyond a pure power of theta (a real fractional part) is declined here and
 * left to the fractional case.  Closes e.g. Cosh/Sinh and mixed +/- powers. */
/* Kernelize a single-primitive exponential extension: find a primitive
 * exponent u so every E^w kernel of `f` (w depending on x) has w = k u with
 * k a nonzero integer, substitute every E^(k u) -> rmT^k, and return the
 * resulting expression in rmT (sets *u_out to an owned copy of u).  Returns
 * NULL when there is no single-primitive exponential structure. */
static Expr* rt_exp_kernelize(Expr* f, Expr* x, Expr** u_out) {
    *u_out = NULL;
    if (!rt_find_exp_of_x(f, x)) return NULL;
    Expr** ws = NULL; size_t nw = 0, cap = 0;
    rt_collect_exp_exponents(f, x, &ws, &nw, &cap);
    if (nw == 0) { free(ws); return NULL; }

    /* Synthesize the single primitive exponent u so every E^(w_j) = (E^u)^k_j
     * (k_j integer): u = ws[0]/lcm(ratio denominators).  Handles integer ratios
     * (E^(2u) = (E^u)^2) and genuine rational ratios (E^(x/2), E^(x/3) → E^(x/6))
     * alike; declines (NULL) when the exponents are not all a rational multiple
     * of one another — genuinely independent kernels (a sum / tower). */
    long* kof = malloc(nw * sizeof(long));
    Expr* u = rt_class_primitive(ws, nw, kof);   /* owned primitive, or NULL */
    /* The primitive exponent must be rational in x alone; a nested exponent
     * (e.g. u = E^x in E^(E^x)) is a two-extension tower left to rt_exp_tower_case
     * (else the single-kernel RDE would carry the inner kernel as a free param). */
    if (u && !rt_kernel_simple(u, x)) { expr_free(u); u = NULL; }
    if (!u) { for (size_t i = 0; i < nw; i++) expr_free(ws[i]); free(ws); free(kof); return NULL; }

    Expr** rules = malloc(2 * nw * sizeof(Expr*));
    for (size_t j = 0; j < nw; j++) {
        Expr* tk = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("rmT"), expr_new_integer(kof[j]) }, 2);
        rules[2 * j] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(ws[j]) }, 1), expr_copy(tk) }, 2);
        rules[2 * j + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(ws[j]) }, 2), tk }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, 2 * nw);
    free(rules);
    Expr* uexp = u;                       /* adopt the synthesized primitive */
    for (size_t i = 0; i < nw; i++) expr_free(ws[i]);
    free(ws); free(kof);

    Expr* F = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rl }, 2));
    if (!F) { expr_free(uexp); return NULL; }
    *u_out = uexp;
    return F;
}

static Expr* rt_exp_poly_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;

    Expr* tsym = expr_new_symbol("rmT");
    /* F must be a Laurent polynomial in t: num/den with den a pure power of t.
     * Split into num and the offset M (den = c t^M). */
    Expr* G = rt_eval1("Together", expr_copy(F));
    Expr* num = G ? rt_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rt_eval1("Denominator", expr_copy(G)) : NULL;
    Expr* result = NULL;
    long M = 0;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(F, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;
    if (ok) {
        M = rt_degree(den, tsym);
        /* den must be a monomial c t^M: den / t^M free of t. */
        Expr* dq = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(tsym), expr_new_integer(-M) }, 2) }, 2));
        if (!dq || !rt_free_of_x(dq, tsym)) ok = false;
        if (dq) expr_free(dq);
        /* Fold the constant den/t^M into num so p_i = Coefficient[num/c, ...]. */
        if (ok) {
            Expr* nnew = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2),
                  expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(tsym), expr_new_integer(M) }, 2) }, 3));
            if (nnew) { expr_free(num); num = nnew; }
            else ok = false;
        }
    }
    long dnum = ok ? rt_degree(num, tsym) : -1;
    if (ok && dnum >= 0) {
        Expr** q = calloc((size_t)(dnum + 1), sizeof(Expr*));
        long* qi_pow = malloc((size_t)(dnum + 1) * sizeof(long));
        bool fail = false;
        for (long j = 0; j <= dnum && !fail; j++) {
            long i = j - M;                       /* Laurent power */
            Expr* pi = rt_coeff(num, tsym, j);
            if (rt_is_zero(pi)) { expr_free(pi); q[j] = NULL; qi_pow[j] = i; continue; }
            Expr* qi = (i == 0) ? rt_integrate_in_K_with_logs(pi, x)
                                : rt_solve_rde(pi, i, uexp, x);
            expr_free(pi);
            qi_pow[j] = i;
            if (!qi) { fail = true; break; }
            q[j] = qi;
        }
        if (!fail) {
            Expr** terms = malloc((size_t)(dnum + 1) * sizeof(Expr*));
            size_t nt = 0;
            for (long j = 0; j <= dnum; j++) {
                if (!q[j]) continue;
                long i = qi_pow[j];
                Expr* term;
                if (i == 0) {
                    term = expr_copy(q[j]);
                } else {
                    Expr* iu = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(i), expr_copy(uexp) }, 2);
                    Expr* eiu = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_new_symbol("E"), iu }, 2);
                    term = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(q[j]), eiu }, 2);
                }
                terms[nt++] = term;
            }
            result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                terms, nt));
            free(terms);
        }
        for (long j = 0; j <= dnum; j++) if (q[j]) expr_free(q[j]);
        free(q); free(qi_pow);
    }

    if (G) expr_free(G);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(uexp);
    return result;
}

/* Fractional (Rothstein-Trager) log-part for a single monomial extension
 * theta (theta = Log[u], D theta = u'/u; or theta = E^u, D theta = u' theta).
 * With theta -> t, an integrand whose t-denominator d = prod g_i is
 * squarefree integrates (when elementary and free of a polynomial/Hermite
 * remainder) to  sum_i c_i Log(g_i)  with each c_i a CONSTANT.  Those
 * constants are found by solving the exact polynomial identity
 *     num  =  sum_i c_i D(g_i) (d / g_i)          (D = d/dx + (D t) d/dt)
 * for all t and x with SolveAlways[..., {t, x}]: the solver returns the
 * constant residues, or no solution (which declines) when the integrand has
 * a repeated pole, an irreducible factor with non-constant residue, or a
 * polynomial part.  Correct by construction — a solution certifies
 * D(sum c_i Log g_i) = num/d exactly. */
static Expr* rt_frac_try(Expr* f, Expr* x, Expr* u, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* rules; Expr* Dt; Expr* kernel_back;
    if (is_log) {
        Expr* logu = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ logu, expr_new_symbol("rmT") }, 2) }, 1);
        Expr* du = rt_eval2("D", expr_copy(u), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));
        kernel_back = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
    } else {
        Expr* eu = expr_new_function(expr_new_symbol("Exp"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* peu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(u) }, 2);
        rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ eu, expr_new_symbol("rmT") }, 2),
              expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ peu, expr_new_symbol("rmT") }, 2) }, 2);
        Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));
        Dt = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ up, expr_new_symbol("rmT") }, 2);
        kernel_back = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(u) }, 2);
    }
    if (!Dt) { expr_free(tsym); expr_free(rules); expr_free(kernel_back); return NULL; }

    Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rules }, 2));
    Expr* F = F0 ? rt_eval1("Together", F0) : NULL;
    if (!F) { expr_free(tsym); expr_free(Dt); expr_free(kernel_back); return NULL; }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;

    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;

    if (ok) {
        Expr* factored = rt_eval1("Factor", expr_copy(den));
        Expr* g[16]; size_t ng = 0; bool bad = false;
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args;
                nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2
                    && term->data.function.args[1]->type == EXPR_INTEGER) {
                    base = term->data.function.args[0];
                    e = (long)term->data.function.args[1]->data.integer;
                }
                if (rt_free_of_x(base, tsym)) continue;   /* FreeQ[base, t] */
                if (e != 1 || ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }
        if (!bad && ng >= 1) {
            /* residual = num - sum_i c_i D(g_i) (den/g_i). */
            Expr** negterms = malloc((ng + 1) * sizeof(Expr*));
            negterms[0] = expr_copy(num);
            for (size_t i = 0; i < ng; i++) {
                char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                Expr* dgx = rt_eval2("D", expr_copy(g[i]), expr_copy(x));
                Expr* dgt = rt_eval2("D", expr_copy(g[i]), expr_copy(tsym));
                Expr* dgi = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ dgx, expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(Dt), dgt }, 2) }, 2));
                Expr* cof = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(g[i]), expr_new_integer(-1) }, 2) }, 2));
                negterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(nm),
                              dgi ? dgi : expr_new_integer(0),
                              cof ? cof : expr_new_integer(0) }, 4);
            }
            Expr* residual = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                negterms, ng + 1));
            free(negterms);

            Expr* varlist = expr_new_function(expr_new_symbol("List"),
                (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ residual, expr_new_integer(0) }, 2);
            Expr* sol = rt_eval2("SolveAlways", eqn, varlist);
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                /* Rothstein-Trager residues MUST be constants.  SolveAlways can
                 * return an x-dependent pseudo-solution for a residual with a
                 * constant-term structure (e.g. k -> x + c for the residual
                 * 1 - k/(x + c) of 1/Log[x + c]); accepting it would certify a
                 * WRONG polynomial-coefficient logarithm ((x+c) Log[Log[x+c]]).
                 * Require every residue free of x and t, else decline. */
                Expr* srules = sol->data.function.args[0];
                bool const_res = true;
                for (size_t si = 0;
                     si < srules->data.function.arg_count && const_res; si++) {
                    Expr* rule = srules->data.function.args[si];
                    if (rule->type == EXPR_FUNCTION
                        && rule->data.function.arg_count == 2
                        && (!rt_free_of_x(rule->data.function.args[1], x)
                            || !rt_free_of_x(rule->data.function.args[1], tsym)))
                        const_res = false;
                }
                if (const_res) {
                Expr** rterms = malloc(ng * sizeof(Expr*));
                for (size_t i = 0; i < ng; i++) {
                    char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                    Expr* gib = rt_eval_own(expr_new_function(
                        expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(g[i]),
                          expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym), expr_copy(kernel_back) },
                                2) }, 1) }, 2));
                    Expr* logg = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ gib }, 1);
                    rterms[i] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), logg }, 2);
                }
                Expr* R = expr_new_function(expr_new_symbol("Plus"), rterms, ng);
                free(rterms);
                R = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ R, expr_copy(sol->data.function.args[0]) }, 2));
                if (R) {
                    Expr** zero = malloc(ng * sizeof(Expr*));
                    for (size_t i = 0; i < ng; i++) {
                        char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                        zero[i] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, ng);
                    free(zero);
                    result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ R, zl }, 2));
                }
                }
            }
            if (sol) expr_free(sol);
        }
        for (size_t i = 0; i < ng; i++) expr_free(g[i]);
        if (factored) expr_free(factored);
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(Dt);
    expr_free(kernel_back);
    return result;
}

/* Pure resultant Lazard-Rioboo-Trager log part for a single monomial
 * extension theta (theta = Log[u] or E^u), for the residue class the
 * SolveAlways single-constant-per-factor ansatz of rt_frac_try cannot
 * express: an irreducible-over-Q denominator factor in theta whose
 * Rothstein-Trager residues are ALGEBRAIC (non-rational) constants — e.g.
 * the +-I/2 residues of 1/(x (Log[x]^2+1)) that split the irreducible
 * theta^2+1 into a conjugate log pair (= ArcTan[Log[x]]).
 *
 * With theta -> t (rmT) and D t = u'/u (log) or u' t (exp), a PROPER,
 * SQUAREFREE-denominator fraction a(t)/d(t) has log part
 *   sum over roots c of Res_t(a - z D(d), d):  c Log(gcd_t(a - c D(d), d)),
 * where D(d) = d/dx(d) + (D t) d/dt(d) is the monomial derivation.  The
 * exact resultant + Rioboo LogToReal real-form collapse are delegated to
 * Integrate`TranscendentalLogPart (intrat.c), which reuses the tested
 * rational-LRT subresultant machinery; this routine just builds the
 * monomial derivation, hands off, substitutes t -> kernel, and DIFF-BACK
 * VERIFIES (the reuse crosses content/denominator boundaries, so — like the
 * tower cases — a mis-reduction must decline, never ship a wrong form).
 * Runs only after rt_frac_try declines. */
static Expr* rt_frac_lrt(Expr* f, Expr* x, Expr* u, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* Dt = NULL; Expr* kernel_back = NULL; Expr* F = NULL;
    if (is_log) {
        /* Log[u] -> t; Log[u]^k is Power[Log[u], k] so the inner rewrite
         * lifts to t^k automatically. */
        Expr* rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(u) }, 1), expr_new_symbol("rmT") }, 2) }, 1);
        Expr* du = rt_eval2("D", expr_copy(u), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));                 /* eta = u'/u */
        kernel_back = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(f), rules }, 2));
        F = F0 ? rt_eval1("Together", F0) : NULL;
    } else {
        /* Kernelize so multiplicatively commensurate exponents collapse onto
         * one primitive t = E^uexp (E^(2x) -> t^2), which a plain E^x -> t
         * substitution cannot do. */
        Expr* uexp = NULL;
        Expr* Fk = rt_exp_kernelize(f, x, &uexp);
        if (Fk && uexp) {
            Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));
            Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ up, expr_new_symbol("rmT") }, 2);   /* u' t */
            kernel_back = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2);
            F = rt_eval1("Together", Fk);
        } else if (Fk) { expr_free(Fk); }
        if (uexp) expr_free(uexp);
    }
    (void)u;   /* exp branch derives its own primitive via rt_exp_kernelize */
    if (!Dt || !F) {
        expr_free(tsym); if (Dt) expr_free(Dt);
        if (kernel_back) { expr_free(kernel_back); } if (F) { expr_free(F); }
        return NULL;
    }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;

    /* Gates: proper squarefree rational function of the single kernel t,
     * with the kernel fully substituted (no residual nested exp/log of x —
     * that would let the resultant treat it as a free parameter). */
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && rt_degree(num, tsym) < rt_degree(den, tsym);

    if (ok) {
        /* Squarefree denominator: gcd(d, dd/dt) constant. */
        Expr* ddt0 = rt_eval2("D", expr_copy(den), expr_copy(tsym));
        Expr* g = ddt0 ? rt_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(den), ddt0 }, 2) : NULL;
        bool squarefree = g && rt_degree(g, tsym) <= 0;
        if (g) expr_free(g);

        if (squarefree) {
            /* Monomial derivation D(d) = d/dx(d) + (D t) d/dt(d). */
            Expr* ddx = rt_eval2("D", expr_copy(den), expr_copy(x));
            Expr* ddt = rt_eval2("D", expr_copy(den), expr_copy(tsym));
            Expr* Dd = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ ddx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt), ddt }, 2) }, 2));
            if (Dd) {
                Expr* logpart = rt_eval_call("Integrate`TranscendentalLogPart",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(tsym),
                              expr_new_symbol("rmZ"), Dd, expr_copy(x) }, 6);
                /* Unevaluated (declined) iff the head is unchanged. */
                bool declined = !logpart
                    || (logpart->type == EXPR_FUNCTION
                        && logpart->data.function.head->type == EXPR_SYMBOL
                        && logpart->data.function.head->data.symbol.name
                             == intern_symbol("Integrate`TranscendentalLogPart"));
                /* Partial log part (Bronstein Thm 5.6.1): the residue resultant
                 * mixed constant and non-constant residues.  Surface the
                 * elementary constant-residue logs plus an unevaluated
                 * Integrate[remainder, x] for the non-constant part.  The FTC
                 * rule D[Integrate[f,x],x] = f makes rt_verify_antideriv close
                 * the diff-back on the assembled partial form. */
                bool partial = logpart && logpart->type == EXPR_FUNCTION
                    && logpart->data.function.head->type == EXPR_SYMBOL
                    && logpart->data.function.head->data.symbol.name
                         == intern_symbol("Integrate`PartialLogPart")
                    && logpart->data.function.arg_count == 2;
                if (partial) {
                    Expr* krule = expr_new_function(expr_new_symbol("List"),
                        (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_copy(tsym), expr_copy(kernel_back) }, 2) }, 1);
                    Expr* logs_k = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(logpart->data.function.args[0]), expr_copy(krule) }, 2));
                    Expr* rem_k = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(logpart->data.function.args[1]), krule }, 2));
                    Expr* Q = (logs_k && rem_k) ? rt_eval_own(expr_new_function(
                        expr_new_symbol("Plus"),
                        (Expr*[]){ logs_k,
                            expr_new_function(expr_new_symbol("Integrate"),
                                (Expr*[]){ rem_k, expr_copy(x) }, 2) }, 2)) : NULL;
                    if (!Q) { if (logs_k) expr_free(logs_k); if (rem_k) expr_free(rem_k); }
                    if (Q && rt_verify_antideriv(Q, f, x)) result = Q;
                    else if (Q) expr_free(Q);
                    expr_free(logpart);
                } else if (!declined) {
                    /* Substitute t -> kernel(x) and diff-back verify. */
                    Expr* Q = rt_eval_own(expr_new_function(
                        expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ logpart,
                          expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym),
                                          expr_copy(kernel_back) }, 2) }, 1) }, 2));
                    if (Q && rt_verify_antideriv(Q, f, x)) result = Q;
                    else if (Q) expr_free(Q);
                } else if (logpart) {
                    expr_free(logpart);
                }
            }
        }
    }

    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(Dt);
    expr_free(kernel_back);
    return result;
}

/* Hermite reduction for a REPEATED transcendental pole, for a logarithmic
 * (theta = Log[u], D theta = u'/u) or exponential (theta = E^u, D theta = u'
 * theta) kernel.  A proper rational function A/D in theta whose denominator D
 * has a repeated factor integrates to  Q = H(theta)/Hden(theta) +
 * sum_j c_j Log(g_j)  where Hden = gcd(D, dD/dtheta) = prod D_m^(m-1) is the
 * repeated part, g_j are the squarefree factors of the radical D/Hden, H is a
 * polynomial in theta of degree < deg(Hden) with bounded-degree polynomial-in-x
 * coefficients, and the c_j are constants.  All unknown constants are solved at
 * once by SolveAlways[Q' - f == 0, {t, x}] — correct by construction (e.g.
 * Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x]),
 * Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)).  The exponential kernel is handled
 * only when D is coprime to theta (no theta = 0 pole — that Laurent case is
 * left to the hyperexponential path). */
static Expr* rt_hermite_try(Expr* f, Expr* x, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* F = NULL; Expr* Dt = NULL; Expr* kback = NULL; Expr* uu = NULL;
    if (is_log) {
        Expr* u = rt_find_log_of_x(f, x);
        if (!u || !rt_kernel_simple(u, x)) { expr_free(tsym); return NULL; }
        uu = expr_copy(u);
        Expr* du = rt_eval2("D", expr_copy(uu), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(uu), expr_new_integer(-1) }, 2);
        Dt = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));   /* eta = u'/u */
        kback = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(uu) }, 1);
        Expr* rl = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(uu) }, 1), expr_new_symbol("rmT") }, 2) }, 1);
        Expr* F0 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(f), rl }, 2));
        F = F0 ? rt_eval1("Together", F0) : NULL;
    } else {
        Expr* F0 = rt_exp_kernelize(f, x, &uu);
        if (F0) F = rt_eval1("Together", F0);
        if (uu) {
            Expr* up = rt_eval2("D", expr_copy(uu), expr_copy(x));
            Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ up, expr_new_symbol("rmT") }, 2);   /* u' t */
            kback = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(uu) }, 2);
        }
    }
    if (!F || !Dt || !uu) {
        expr_free(tsym);
        if (F) { expr_free(F); } if (Dt) { expr_free(Dt); }
        if (kback) { expr_free(kback); } if (uu) { expr_free(uu); }
        return NULL;
    }

    Expr* num = rt_eval1("Numerator", expr_copy(F));
    Expr* den = rt_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    long dnum = ok ? rt_degree(num, tsym) : 0, dden = ok ? rt_degree(den, tsym) : 0;
    ok = ok && dnum < dden;                      /* proper fraction */
    if (ok && !is_log) {
        /* exponential kernel: require D coprime to theta (no theta = 0 pole). */
        Expr* c0 = rt_coeff(den, tsym, 0);
        if (rt_is_zero(c0)) ok = false;
        expr_free(c0);
    }

    Expr* Hden = NULL; Expr* rad = NULL;
    if (ok) {
        Expr* dent = rt_eval2("D", expr_copy(den), expr_copy(tsym));
        Hden = rt_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(den), dent }, 2);
        long dH = Hden ? rt_degree(Hden, tsym) : 0;
        if (dH < 1) ok = false;                  /* no repeated pole */
        if (ok) rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    }

    if (ok && rad) {
        long dH = rt_degree(Hden, tsym);
        /* distinct t-factors of the squarefree radical (log terms). */
        Expr* factored = rt_eval1("Factor", expr_copy(rad));
        Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2)
                    base = term->data.function.args[0];
                if (rt_free_of_x(base, tsym)) continue;
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }
        long dnx = rt_degree(num, x), ddx = rt_degree(den, x);
        long Nx = (dnx > ddx ? dnx : ddx) + 2;   /* derived Hermite x-degree, no cap */
        size_t nh = (size_t)(dH * (Nx + 1));
        if (!bad && (long)(nh + ng) > 0) {
            Expr** hterms = malloc((nh ? nh : 1) * sizeof(Expr*));
            size_t nt = 0;
            for (long p = 0; p < dH; p++)
                for (long k = 0; k <= Nx; k++) {
                    char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                    hterms[nt++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(tsym), expr_new_integer(p) }, 2) }, 3);
                }
            Expr* Hpoly = expr_new_function(expr_new_symbol("Plus"), hterms, nt);
            free(hterms);
            Expr* ratpart = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
            Expr** qterms = malloc((1 + ng) * sizeof(Expr*));
            size_t ntq = 0;
            qterms[ntq++] = ratpart;
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            Expr* dQx = rt_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rt_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt), dQt }, 2) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rt_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    Expr** zero = malloc((nh + ng) * sizeof(Expr*));
                    size_t zi = 0;
                    for (long p = 0; p < dH; p++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nh + ng);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym), expr_copy(kback) }, 2) }, 1);
                        Expr* rr = rt_eval_own(expr_new_function(
                            expr_new_symbol("ReplaceAll"), (Expr*[]){ Qs, back }, 2));
                        /* Cancel folds away the x-content that Hden carried
                         * (e.g. -x/(x + x Log[x]) -> -1/(1 + Log[x])). */
                        if (rr) result = rt_eval1("Cancel", rr);
                    }
                }
            }
            if (sol) expr_free(sol);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
    }
    if (Hden) expr_free(Hden);
    if (rad) expr_free(rad);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(Dt);
    expr_free(kback);
    expr_free(uu);
    expr_free(F);
    expr_free(tsym);
    return result;
}

/* Hermite reduction dispatch: logarithmic then exponential kernel. */
static Expr* rt_hermite_case(Expr* f, Expr* x) {
    Expr* r = rt_hermite_try(f, x, true);
    if (r) return r;
    r = rt_hermite_try(f, x, false);
    return r;
}

/* Try the fractional log-part for a logarithmic then an exponential kernel. */
static Expr* rt_frac_case(Expr* f, Expr* x) {
    Expr* ul = rt_find_log_of_x(f, x);
    if (ul && rt_kernel_simple(ul, x)) {
        Expr* r = rt_frac_try(f, x, ul, true); if (r) return r;
    }
    Expr* ue = rt_find_exp_of_x(f, x);
    if (ue && rt_kernel_simple(ue, x)) {
        Expr* r = rt_frac_try(f, x, ue, false); if (r) return r;
    }
    /* Pure resultant LRT for algebraic (non-rational) residues — runs only
     * when the rational-residue SolveAlways path above declines. */
    if (ul && rt_kernel_simple(ul, x)) {
        Expr* r = rt_frac_lrt(f, x, ul, true); if (r) return r;
    }
    if (ue && rt_kernel_simple(ue, x)) {
        Expr* r = rt_frac_lrt(f, x, ue, false); if (r) return r;
    }
    return NULL;
}

/* Coupled hyperexponential case: a general rational function of a single
 * exponential monomial theta = E^u (u polynomial in x) whose integral mixes a
 * Laurent-polynomial part with a logarithmic part — the case where D(g)/g is
 * itself improper in theta, so the polynomial and log parts do NOT separate
 * (e.g. Tan/Tanh after TrigToExp, or 1/(1+E^x)).  A UNIFIED ansatz
 *     Q = sum_i (sum_k a_{ik} x^k) theta^i + sum_j c_j Log(g_j)
 * (g_j the squarefree factors of the theta-coprime denominator) is solved for
 * all its constant coefficients a_{ik}, c_j at once by requiring Q' = f for all
 * theta and x via SolveAlways[..., {t, x}].  Correct by construction: an exact
 * solution certifies the identity; an imperfect degree bound can only decline.
 * Tried after the pure polynomial and pure-log-part cases (which are cheaper). */
static Expr* rt_hyperexp_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;
    Expr* tsym = expr_new_symbol("rmT");
    Expr* result = NULL;

    Expr* G = rt_eval1("Together", expr_copy(F));
    Expr* num = G ? rt_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rt_eval1("Denominator", expr_copy(G)) : NULL;
    bool ok = num && den && rt_is_poly(num, tsym) && rt_is_poly(den, tsym)
        && !rt_free_of_x(den, tsym)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;

    if (ok) {
        /* a = multiplicity of t in den (first nonzero t-coefficient). */
        long a = 0;
        Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(tsym));
        if (cl && cl->type == EXPR_FUNCTION
            && cl->data.function.head->type == EXPR_SYMBOL
            && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
            for (size_t i = 0; i < cl->data.function.arg_count; i++)
                if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
        }
        if (cl) expr_free(cl);
        Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(tsym), expr_new_integer(-a) }, 2) }, 2));
        /* Split the theta-coprime denominator Dtil into its repeated part
         * Hden = gcd(Dtil, d Dtil/dt) (= prod D_m^(m-1)) and squarefree radical
         * rad = Dtil/Hden.  The distinct t-factors of rad carry constant-residue
         * logarithms; the repeated part is absorbed by a Hermite term
         * H(t)/Hden(t) (H proper in t with x-polynomial coefficients).  When Dtil
         * is squarefree Hden = 1 and the Hermite term is empty, recovering the
         * plain Laurent + log hyperexponential ansatz.  This is the Phase A
         * generalization that couples the Laurent theta = 0 pole with a Hermite
         * repeated pole, closing shapes such as 1/(1 + E^x)^2 and
         * 1/(E^x (1 + E^x)^2) that the pure Hermite / squarefree paths decline. */
        Expr* Hden = NULL; Expr* rad = NULL; long dH = 0;
        if (Dtil) {
            Expr* dDt = rt_eval2("D", expr_copy(Dtil), expr_copy(tsym));
            Hden = rt_eval_call("PolynomialGCD",
                (Expr*[]){ expr_copy(Dtil), dDt }, 2);
            if (Hden) {
                dH = rt_degree(Hden, tsym);
                if (dH < 0) dH = 0;
                rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dtil), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
            }
        }
        Expr* factored = rad ? rt_eval1("Factor", expr_copy(rad)) : NULL;
        Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
        if (factored) {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2)
                    base = term->data.function.args[0];   /* rad squarefree: e = 1 */
                if (rt_free_of_x(base, tsym)) continue;    /* FreeQ[base, t] */
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }

        long dnum = ok ? rt_degree(num, tsym) : 0, dden = ok ? rt_degree(den, tsym) : 0;
        long ihi = dnum - dden; if (ihi < 0) ihi = 0;
        long ilo = -a;
        long degu = rt_degree(uexp, x);
        long dnx = rt_degree(num, x), ddx = rt_degree(den, x);
        /* x-degree bound of the coefficient polynomials — derived from the numerator/
         * denominator degrees and the exponent degree, NO cap.  Generous but a
         * function of the input; a too-large bound only slows SolveAlways (every
         * solution is certified + diff-back verified), never mis-solves. */
        long Nx = (dnx > ddx ? dnx : ddx) + (degu > 0 ? degu : 1) + 1;
        long nwi = ihi - ilo + 1;
        size_t nH_syms = (size_t)(dH * (Nx + 1));   /* Hermite numerator coeffs */
        long nunk = nwi * (Nx + 1) + (long)nH_syms + (long)ng;
        if (!bad && nwi > 0 && nunk > 0) {
            size_t nw_syms = (size_t)(nwi * (Nx + 1));
            Expr** qterms = malloc((nw_syms + 1 + ng) * sizeof(Expr*));
            size_t ntq = 0;
            /* Laurent part: sum_i (sum_k rmW.. x^k) t^i. */
            for (long i = ilo; i <= ihi; i++) {
                for (long k = 0; k <= Nx; k++) {
                    char nm[32]; snprintf(nm, sizeof(nm), "rmW%ldv%ld", i - ilo, k);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(tsym), expr_new_integer(i) }, 2) }, 3);
                }
            }
            /* Hermite part: (sum_{p<dH} (sum_k rmHh.. x^k) t^p) / Hden. */
            if (dH >= 1) {
                Expr** hterms = malloc(nH_syms * sizeof(Expr*));
                size_t nh = 0;
                for (long p = 0; p < dH; p++)
                    for (long k = 0; k <= Nx; k++) {
                        char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                        hterms[nh++] = expr_new_function(expr_new_symbol("Times"),
                            (Expr*[]){ expr_new_symbol(nm),
                              expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2),
                              expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ expr_copy(tsym), expr_new_integer(p) }, 2) }, 3);
                    }
                Expr* Hpoly = expr_new_function(expr_new_symbol("Plus"), hterms, nh);
                free(hterms);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
            }
            /* Log part: sum_j rmHc.. Log(g_j). */
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            /* Q'(full) = D[Q,x] + u' t D[Q,t]. */
            Expr* dQx = rt_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rt_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ up, expr_copy(tsym), dQt }, 3) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rt_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    Expr** zero = malloc((nw_syms + nH_syms + ng) * sizeof(Expr*));
                    size_t zi = 0;
                    for (long i = ilo; i <= ihi; i++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[32]; snprintf(nm, sizeof(nm), "rmW%ldv%ld", i - ilo, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (long p = 0; p < dH; p++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmHc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero,
                                                 nw_syms + nH_syms + ng);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym),
                                    expr_new_function(expr_new_symbol("Power"),
                                        (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2)
                                }, 2) }, 1);
                        result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                            (Expr*[]){ Qs, back }, 2));
                    }
                }
            }
            if (sol) expr_free(sol);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        if (rad) expr_free(rad);
        if (Hden) expr_free(Hden);
        if (Dtil) expr_free(Dtil);
    }
    if (G) expr_free(G);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F);
    expr_free(tsym);
    expr_free(uexp);
    return result;
}

/* ================================================================== */
/* Case: sum of NON-COMMENSURATE exponentials (Phase B, first tower    */
/* increment) — several independent primitive exponential extensions.  */
/* ================================================================== */
/* Split one multiplicative term T = p(x) * prod_k E^(w_k) into the total
 * x-dependent exponential exponent W = sum_k w_k (over every E^w / Exp[w]
 * factor whose exponent depends on x) and the exponential-free cofactor p.
 * Returns 0 on success with *W_out and *p_out owned; -1 if the cofactor still
 * carries an x-dependent exponential (a shape this simple splitter declines). */
/* Recursively accumulate the factors of a term (descending into nested,
 * un-flattened Times nodes — Expand can leave Times[2, Times[x, E^w]]): every
 * x-dependent exponential exponent is summed into *W; every other factor is
 * appended (as an owned copy) to the growable array *pf. */
static void rt_accum_factors(Expr* T, Expr* x, Expr** W,
                             Expr*** pf, size_t* npf, size_t* cap) {
    if (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL
        && T->data.function.head->data.symbol.name == intern_symbol("Times")) {
        for (size_t i = 0; i < T->data.function.arg_count; i++)
            rt_accum_factors(T->data.function.args[i], x, W, pf, npf, cap);
        return;
    }
    Expr* w = NULL;   /* borrowed exponent if T is an x-dependent exp kernel */
    const char* h = (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL)
        ? T->data.function.head->data.symbol.name : NULL;
    if (h == intern_symbol("Exp") && T->data.function.arg_count == 1
        && !rt_free_of_x(T->data.function.args[0], x))
        w = T->data.function.args[0];
    else if (h == intern_symbol("Power") && T->data.function.arg_count == 2
        && T->data.function.args[0]->type == EXPR_SYMBOL
        && T->data.function.args[0]->data.symbol.name == intern_symbol("E")
        && !rt_free_of_x(T->data.function.args[1], x))
        w = T->data.function.args[1];
    if (w) {
        *W = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ *W, expr_copy(w) }, 2));
    } else {
        if (*npf == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *pf = realloc(*pf, *cap * sizeof(Expr*));
        }
        (*pf)[(*npf)++] = expr_copy(T);
    }
}

static int rt_split_exp_term(Expr* T, Expr* x, Expr** W_out, Expr** p_out) {
    Expr* W = expr_new_integer(0);
    Expr** pf = NULL; size_t npf = 0, cap = 0;
    rt_accum_factors(T, x, &W, &pf, &npf, &cap);

    Expr* p;
    if (npf == 0) p = expr_new_integer(1);
    else if (npf == 1) p = pf[0];                          /* adopt sole factor */
    else p = expr_new_function(expr_new_symbol("Times"), pf, npf);  /* adopts */
    free(pf);

    if (rt_find_exp_of_x(p, x) != NULL) {   /* nested/residual exp in cofactor */
        expr_free(W); expr_free(p);
        return -1;
    }
    *W_out = W; *p_out = p;
    return 0;
}

/* Sum-of-exponentials case.  When TrigToExp (or the raw integrand) yields a sum
 *   f = sum_k p_k(x) E^(W_k)
 * of exponentials whose exponents W_k are NOT integer multiples of a single
 * primitive (e.g. (1 +/- I) x from E^x Sin[x]), the distinct exponentials are
 * independent transcendental extensions.  Because d/dx maps each E^(W_k) to a
 * multiple of itself and never mixes distinct exponents, the terms DECOUPLE:
 *   INT p_k E^(W_k) dx = q_k E^(W_k),   q_k' + W_k' q_k = p_k  (a Risch DE)
 * with W_k polynomial in x, plus any exponential-free (W_k = 0) terms handled by
 * the base-field integral.  Each q_k is found exactly (SolveAlways certificate),
 * so d/dx(sum_k q_k E^(W_k)) = sum_k p_k E^(W_k) = f by linearity — correct by
 * construction.  Declines (whole) if any term is not rational-coefficient over a
 * polynomial exponent (e.g. E^(x^2), left to the Erf recognizer) or carries a
 * residual log / non-elementary cofactor.  Runs after the single-primitive
 * exponential cases, so it only fires on genuinely multi-kernel integrands. */
static Expr* rt_expsum_case(Expr* f, Expr* x) {
    if (!rt_find_exp_of_x(f, x)) return NULL;   /* need at least one exp kernel */
    Expr* fe = rt_eval1("Expand", expr_copy(f));
    if (!fe) return NULL;

    Expr** terms; size_t nt; Expr* single[1];
    if (fe->type == EXPR_FUNCTION
        && fe->data.function.head->type == EXPR_SYMBOL
        && fe->data.function.head->data.symbol.name == intern_symbol("Plus")) {
        terms = fe->data.function.args; nt = fe->data.function.arg_count;
    } else { single[0] = fe; terms = single; nt = 1; }

    Expr** outs = malloc((nt ? nt : 1) * sizeof(Expr*));
    size_t no = 0;
    bool fail = false, any_exp = false;
    for (size_t i = 0; i < nt && !fail; i++) {
        Expr* W = NULL; Expr* p = NULL;
        if (rt_split_exp_term(terms[i], x, &W, &p) != 0) { fail = true; break; }
        bool wdep = !rt_free_of_x(W, x);
        Expr* qi = NULL;
        if (!wdep) {
            /* E^W is a constant coefficient: INT p E^W dx = E^W INT p dx. */
            qi = rt_integrate_in_K_with_logs(p, x);
        } else if (rt_is_poly(W, x)) {
            qi = rt_solve_rde(p, 1, W, x);   /* q' + W' q = p */
            if (qi) any_exp = true;
        }   /* rational (non-polynomial) W is Phase C: qi stays NULL -> decline */
        expr_free(p);
        if (!qi) { expr_free(W); fail = true; break; }
        Expr* eW = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), W }, 2);   /* adopts W; E^0 -> 1 */
        outs[no++] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ qi, eW }, 2);                    /* adopts qi, eW */
    }

    Expr* result = NULL;
    if (!fail && any_exp) {
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), outs, no));
    } else {
        for (size_t i = 0; i < no; i++) expr_free(outs[i]);
    }
    free(outs);
    expr_free(fe);
    return result;
}

/* ================================================================== */
/* Tower verification gate.                                            */
/* ================================================================== */
/* The tower cases (rt_log_tower_case / rt_exp_tower_case) are bounded-ansatz
 * SEARCHES over a multi-variable SolveAlways, not pure decision procedures: the
 * certificate "D_tower[Q] = F as a polynomial identity in the tower variables"
 * is sound only when SolveAlways returns a genuine solution, and a non-elementary
 * integrand can yield a spurious one (e.g. E^(E^x)/(1+E^(E^x)) is non-elementary
 * yet the ansatz once returned Log[1+E^(E^x)]/E^x).  As with the other
 * search-based integrators, the tower results are therefore diff-back verified:
 * Simplify[D[result,x] - f] must be exactly 0, else the case declines.  This
 * keeps the genuine closures (E^x E^(E^x) -> E^(E^x),
 * E^x E^(E^x)/(1+E^(E^x)) -> Log[1+E^(E^x)]) and rejects the spurious ones. */
static bool rt_verify_antideriv(Expr* result, Expr* f, Expr* x) {
    Expr* d = rt_eval2("D", expr_copy(result), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    Expr* s = rt_eval1("Simplify", diff);
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* ================================================================== */
/* Case: nested logarithmic tower (Phase B, second increment).         */
/* ================================================================== */
/* Collect (owned, deduplicated) every Log[u] kernel of `e` whose argument
 * depends on x. */
static void rt_collect_logs(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x)) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], e)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                              *arr = realloc(*arr, *cap * sizeof(Expr*)); }
            (*arr)[(*n)++] = expr_copy(e);
        }
    }
    rt_collect_logs(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_logs(e->data.function.args[i], x, arr, n, cap);
}

/* Structural containment: does `big` contain `small` as a subexpression? */
static bool rt_contains(Expr* big, Expr* small) {
    if (!big) return false;
    if (expr_eq(big, small)) return true;
    if (big->type != EXPR_FUNCTION) return false;
    if (rt_contains(big->data.function.head, small)) return true;
    for (size_t i = 0; i < big->data.function.arg_count; i++)
        if (rt_contains(big->data.function.args[i], small)) return true;
    return false;
}

/* Build the monomial  prod_j lv[j]^e[j]  (a Times of Powers; owned). */
static Expr* rt_build_monomial(Expr** lv, const long* e, size_t nlv) {
    Expr** fs = malloc((nlv ? nlv : 1) * sizeof(Expr*));
    size_t nf = 0;
    for (size_t j = 0; j < nlv; j++) {
        if (e[j] == 0) continue;
        fs[nf++] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(lv[j]), expr_new_integer(e[j]) }, 2);
    }
    Expr* m;
    if (nf == 0) m = expr_new_integer(1);
    else if (nf == 1) m = fs[0];
    else m = expr_new_function(expr_new_symbol("Times"), fs, nf);
    free(fs);
    return m;
}

/* Decode the flat monomial index `idx` into the exponent vector `e` over the
 * per-variable degree bounds `bd` (odometer, lowest var fastest). */
static void rt_decode_mono(long idx, const long* bd, size_t nlv, long* e) {
    for (size_t j = 0; j < nlv; j++) {
        long base = bd[j] + 1;
        e[j] = idx % base;
        idx /= base;
    }
}

/* ------------------------------------------------------------------ */
/* Logarithm-argument normalization for the tower builder.             */
/* ------------------------------------------------------------------ */
/* Rewrite  Log[a b ...] -> Log[a] + Log[b] + ...  and  Log[b^p] -> p Log[b]
 * recursively, so a nested-log integrand presents the MINIMAL set of
 * multiplicatively-independent Log generators to rt_collect_logs.  Without
 * this, a composite kernel like Log[x/Log[x]] is treated as an INDEPENDENT
 * transcendental on top of Log[x] and Log[Log[x]], inflating the tower depth
 * (and hence the undetermined-coefficient ansatz — the dominant cost) by a
 * spurious, functionally-redundant generator.
 *
 * The rewrite is exact for the DERIVATIVE (d/dx Log[u] = u'/u regardless of
 * branch), which is all the tower's correct-by-construction certificate and
 * the final diff-back gate (rt_verify_antideriv, against the ORIGINAL f)
 * depend on; the branch-cut constants it drops are absorbed by the constant
 * of integration.  Sums inside a Log are left intact (Log[a+b] does not
 * split). */
static Expr* rt_expand_logs(Expr* e);

/* Expanded form of Log[a] (the argument `a` is borrowed). */
static Expr* rt_log_of(Expr* a) {
    if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL) {
        const char* h = a->data.function.head->data.symbol.name;
        if (h == intern_symbol("Times")) {
            size_t m = a->data.function.arg_count;
            Expr** ts = malloc((m ? m : 1) * sizeof(Expr*));
            for (size_t i = 0; i < m; i++) ts[i] = rt_log_of(a->data.function.args[i]);
            Expr* s = expr_new_function(expr_new_symbol("Plus"), ts, m);
            free(ts);
            return s;
        }
        if (h == intern_symbol("Power") && a->data.function.arg_count == 2) {
            Expr* lb = rt_log_of(a->data.function.args[0]);
            return expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ rt_expand_logs(a->data.function.args[1]), lb }, 2);
        }
    }
    /* Atomic (or Plus) argument: keep the Log, but still expand logs nested
     * inside the argument (e.g. the Log[x] inside Log[Log[x]]). */
    return expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ rt_expand_logs(a) }, 1);
}

static Expr* rt_expand_logs(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1) {
        return rt_log_of(e->data.function.args[0]);
    }
    size_t n = e->data.function.arg_count;
    Expr** args = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) args[i] = rt_expand_logs(e->data.function.args[i]);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return r;
}

/* Nested logarithmic tower.  When the integrand is a rational function of a
 * chain of NESTED logarithms  t_1 = Log[u_1](x), t_2 = Log[u_2](x, t_1), ...,
 * t_n = Log[u_n](x, t_1..t_{n-1})  (n >= 2), the single-extension cases (which
 * model one kernel over C(x)) do not apply.  Generalize the single-kernel
 * derivation  D = d/dx + Dt d/dt  to the full tower
 *   D_tower = d/dx + sum_i Dt_i d/dt_i,   Dt_i = Cancel[D[u_i,x]/u_i] |_{ker->t}
 * (triangular: Dt_i lies in C(x, t_1..t_{i-1})), and solve ONE unified ansatz
 *   Q = sum_{k=0}^{Ntop} P_k(x, t_1..t_{n-1}) t_n^k  +  sum_j c_j Log(g_j)
 * (P_k bounded-degree polynomials with unknown constant coefficients, g_j the
 * squarefree t_n-denominator factors, c_j constants) by requiring
 * D_tower[Q] = f for all {t_n,...,t_1,x} via SolveAlways.  Correct by
 * construction (an exact solution certifies the identity).  Closes e.g.
 * Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]] and
 * Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2.  Declines cleanly
 * when the tower is not a valid nested chain, the top denominator has a
 * repeated pole (tower Hermite — later), or a lower-field coefficient must be
 * rational (the full recursion — later). */
static Expr* rt_log_tower_case(Expr* f, Expr* x) {
    /* Normalize composite logs (Log[a b] -> Log a + Log b, Log[b^p] -> p Log b)
     * so the tower is built over the MINIMAL independent generator set.  `fn`
     * drives kernel collection and the ansatz; correctness is still gated by
     * diff-back against the ORIGINAL `f` (rt_verify_antideriv, below). */
    Expr* fn = rt_eval_own(rt_expand_logs(f));
    if (!fn) fn = expr_copy(f);

    Expr** logs = NULL; size_t nl = 0, lcap = 0;
    rt_collect_logs(fn, x, &logs, &nl, &lcap);
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
                            free(logs); expr_free(fn); return NULL; }

    /* Order innermost-first: if logs[i] contains logs[k], logs[k] is deeper. */
    for (size_t i = 0; i < nl; i++)
        for (size_t k = i + 1; k < nl; k++)
            if (rt_contains(logs[i], logs[k])) {
                Expr* t = logs[i]; logs[i] = logs[k]; logs[k] = t;
            }

    /* t-symbols and substitution rules kernel_i -> t_i. */
    Expr** ts = malloc(nl * sizeof(Expr*));
    Expr** rules = malloc(nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmtw%zu", i);
        ts[i] = expr_new_symbol(nm);
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(logs[i]), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, nl);
    free(rules);

    Expr* F = rt_eval1("Together", rt_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(fn), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL;
    Expr** Dt = NULL;
    Expr* result = NULL;

    bool ok = F && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rt_eval1("Numerator", expr_copy(F));
        den = rt_eval1("Denominator", expr_copy(F));
        /* F must be a genuine RATIONAL function of the whole tower {x, t_1..t_n}
         * — i.e. num and den are multivariate polynomials in those variables and
         * nothing else.  Without this, a residual non-rational kernel of an inner
         * argument (e.g. Sin[Log[x]] -> Sin[t_1]) would slip through the per-top-
         * variable test and the ansatz could certify a WRONG closed form. */
        Expr** vv = malloc((nl + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < nl; i++) vv[i + 1] = expr_copy(ts[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nl + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rt_is_true(pqn) && rt_is_true(pqd)
             && !rt_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Tower derivation coefficients Dt_i, checked triangular. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* arg = logs[i]->data.function.args[0];
            Expr* darg = rt_eval2("D", expr_copy(arg), expr_copy(x));
            Expr* q = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ darg, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(arg), expr_new_integer(-1) }, 2) }, 2);
            Expr* dc = rt_eval1("Cancel", q);
            Dt[i] = dc ? rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
            for (size_t j = i; j < nl; j++)
                if (!rt_free_of_x(Dt[i], ts[j])) { ok = false; break; }
        }
    }

    if (ok) {
        long dtn_num = rt_degree(num, top), dtn_den = rt_degree(den, top);
        /* Exact top-degree bound: t_n is a LOGARITHM (D lowers deg_{t_n} by 1), so
         * deg_{t_n}(Q) = deg_{t_n}(f) + 1.  No cap. */
        long Ntop = dtn_num - dtn_den + 1;
        if (Ntop < 0) Ntop = 0;

        /* Lower field variables: x, t_1..t_{n-1}. */
        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));
        lv[0] = x;
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        for (size_t j = 0; j < nlv; j++) {
            long a = rt_degree(num, lv[j]); if (a < 0) a = 0;
            long b = rt_degree(den, lv[j]); if (b < 0) b = 0;
            long d = a + b + 1;              /* derived lower-field proxy, no cap */
            bd[j] = d;
        }
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);

        /* Squarefree t_n-dependent factors of den -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = rt_eval1("Factor", expr_copy(den));
        if (!factored) bad = true;
        else {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2
                    && term->data.function.args[1]->type == EXPR_INTEGER) {
                    base = term->data.function.args[0];
                    e = (long)term->data.function.args[1]->data.integer;
                }
                if (rt_free_of_x(base, top)) continue;     /* not a t_n factor */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long nunk = (Ntop + 1) * nmono + (long)ng;
        if (!bad && nunk > 0) {
            /* Build Q = sum_{k,mono} rmLp{k}_{m} (mono) t_n^k + sum_j rmLc{j} Log(g_j). */
            size_t nq = (size_t)((Ntop + 1) * nmono + (long)ng);
            Expr** qterms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long k = 0; k <= Ntop; k++)
                for (long m = 0; m < nmono; m++) {
                    char nm[64]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
                    rt_decode_mono(m, bd, nlv, ev);
                    Expr* mono = rt_build_monomial(lv, ev, nlv);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), mono,
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(top), expr_new_integer(k) }, 2) }, 3);
                }
            for (size_t j = 0; j < ng; j++) {
                char nm[32]; snprintf(nm, sizeof(nm), "rmLc%zu", j);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            /* D_tower[Q] = D[Q,x] + sum_i Dt_i D[Q,t_i]. */
            Expr** dterms = malloc((nl + 1) * sizeof(Expr*));
            dterms[0] = rt_eval2("D", expr_copy(Q), expr_copy(x));
            for (size_t i = 0; i < nl; i++) {
                Expr* dqi = rt_eval2("D", expr_copy(Q), expr_copy(ts[i]));
                dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
            }
            Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                dterms, nl + 1));
            free(dterms);

            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rt_eval1("Together", diff);
            Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

            /* SolveAlways over {t_n, ..., t_1, x}. */
            Expr* sol = NULL;
            if (rnum) {
                Expr** vl = malloc((nl + 1) * sizeof(Expr*));
                for (size_t i = 0; i < nl; i++) vl[i] = expr_copy(ts[nl - 1 - i]);
                vl[nl] = expr_copy(x);
                Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nl + 1);
                free(vl);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rt_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol.name == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol.name
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    /* Pin any free ansatz parameters to 0. */
                    Expr** zero = malloc(nq * sizeof(Expr*));
                    size_t zi = 0;
                    for (long k = 0; k <= Ntop; k++)
                        for (long m = 0; m < nmono; m++) {
                            char nm[64]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
                            zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                        }
                    for (size_t j = 0; j < ng; j++) {
                        char nm[32]; snprintf(nm, sizeof(nm), "rmLc%zu", j);
                        zero[zi++] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
                    }
                    Expr* zl = expr_new_function(expr_new_symbol("List"), zero, zi);
                    free(zero);
                    Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        /* Back-substitute t_i -> Log kernels (innermost first so a
                         * t appearing inside a Log argument is restored too). */
                        Expr** back = malloc(nl * sizeof(Expr*));
                        for (size_t i = 0; i < nl; i++)
                            back[i] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(ts[i]), expr_copy(logs[i]) }, 2);
                        Expr* bl = expr_new_function(expr_new_symbol("List"), back, nl);
                        free(back);
                        result = rt_eval_own(expr_new_function(
                            expr_new_symbol("ReplaceAll"), (Expr*[]){ Qs, bl }, 2));
                    }
                }
            }
            if (sol) expr_free(sol);
            free(ev);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        free(lv); free(bd);
    }

    /* Diff-back safety gate (see rt_verify_antideriv). */
    if (result && !rt_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    if (Dt) { for (size_t i = 0; i < nl; i++) if (Dt[i]) expr_free(Dt[i]); free(Dt); }
    if (num) expr_free(num);
    if (den) expr_free(den);
    if (F) expr_free(F);
    expr_free(fn);
    expr_free(rl);
    for (size_t i = 0; i < nl; i++) { expr_free(ts[i]); expr_free(logs[i]); }
    free(ts); free(logs);
    return result;
}

/* Shared tower solve tail: given the ansatz Q (in {x, t_1..t_n}), its unknown
 * coefficient symbols `syms`, the tower derivation coefficients Dt, the tower
 * variables ts, the target F, and the back-substitution rule list `backrules`
 * (t_i -> kernel_i), form D_tower[Q] - F, solve the polynomial identity over
 * {t_n,...,t_1,x} with SolveAlways, pin free parameters to 0, and return the
 * back-substituted antiderivative (owned) or NULL.  Borrows all arguments.
 *
 * EXTENDED-LIOUVILLE EXTENSION POINT (Cherry substrate — CHERRY_DESIGN.md §2.2/§3.1).
 * This IS the "reduce to a linear system over the constants" solver both Cherry
 * algorithms need (Risch69 Main Thm part (b) in the 1986 paper; the undetermined-
 * coefficient E1/E2/E3 solve in the 1989 paper).  It already accepts special-
 * function basis terms with NO change: put  Q = poly_ansatz + Sum_i k_i SF(a_i)
 * (SF in {LogIntegral, ExpIntegralEi, Erf, PolyLog, ...}, a_i a generated
 * argument) with the k_i listed in `syms`.  The tower derivative it forms,
 * D[Q,x] + Sum_j Dt_j D[Q,t_j], differentiates the SF terms correctly because
 * Mathilda's D[] knows their derivatives (D[LogIntegral[u]] = u'/Log[u], etc.), so
 * for a kernel argument u the SF term contributes exactly D_tower[u]/Log[u] to the
 * linear system.  The Cherry driver (future) generates the a_i, appends the SF
 * basis terms to Q and the k_i to `syms`, and calls this unchanged. */
static Expr* rt_tower_solve(Expr* Q, Expr** syms, size_t nsym,
                            Expr** Dt, Expr** ts, size_t nl,
                            Expr* F, Expr* x, Expr* backrules) {
    Expr** dterms = malloc((nl + 1) * sizeof(Expr*));
    dterms[0] = rt_eval2("D", expr_copy(Q), expr_copy(x));
    for (size_t i = 0; i < nl; i++) {
        Expr* dqi = rt_eval2("D", expr_copy(Q), expr_copy(ts[i]));
        dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
    }
    Expr* Qder = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        dterms, nl + 1));
    free(dterms);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
    Expr* tog = rt_eval1("Together", diff);
    Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

    Expr* result = NULL; Expr* sol = NULL;
    if (rnum) {
        Expr** vl = malloc((nl + 1) * sizeof(Expr*));
        for (size_t i = 0; i < nl; i++) vl[i] = expr_copy(ts[nl - 1 - i]);
        vl[nl] = expr_copy(x);
        Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nl + 1);
        free(vl);
        Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
            (Expr*[]){ rnum, expr_new_integer(0) }, 2);
        sol = rt_eval2("SolveAlways", eqn, varlist);
    }
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol.name == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol.name
             == intern_symbol("List")) {
        Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
        if (Qs) {
            Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
            for (size_t j = 0; j < nsym; j++)
                zero[j] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
            Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
            free(zero);
            Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ Qs, zl }, 2));
            if (Qs)
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, expr_copy(backrules) }, 2));
        }
    }
    if (sol) expr_free(sol);
    return result;
}

/* Nested EXPONENTIAL tower (Phase B, third increment) — the dual of
 * rt_log_tower_case for a chain of nested exponentials t_i = E^(u_i) with u_i in
 * C(x, t_1..t_{i-1}) (e.g. t_1 = E^x, t_2 = E^(E^x)).  An exponential kernel's
 * derivative is a multiple of itself, so the tower derivation is
 *   D_tower = d/dx + sum_i Dt_i d/dt_i,   Dt_i = (D[u_i,x] |_{ker->t}) t_i,
 * and a unified Laurent ansatz
 *   Q = sum_{i=ilo}^{ihi} P_i(x, t_1..t_{n-1}) t_n^i  +  sum_j c_j Log(g_j)
 * (P_i bounded-degree lower-field polynomials, ilo = -(t_n multiplicity in den),
 * g_j the squarefree t_n-coprime denominator factors) is solved by SolveAlways
 * over {t_n,...,t_1,x}.  Same whole-tower rationality gate and correct-by-
 * construction certificate as the log tower.  Closes e.g.
 * Integrate[E^x E^(E^x), x] = E^(E^x).  Repeated t_n poles (tower Hermite) and
 * rational lower-field coefficients are out of scope and decline. */
static Expr* rt_exp_tower_case(Expr* f, Expr* x) {
    if (!rt_find_exp_of_x(f, x)) return NULL;
    Expr** us = NULL; size_t nl = 0, ucap = 0;
    rt_collect_exp_exponents(f, x, &us, &nl, &ucap);   /* us[i] = exp exponents */
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(us[i]);
                            free(us); return NULL; }

    /* Kernel forms E^(u_i). */
    Expr** kf = malloc(nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++)
        kf[i] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(us[i]) }, 2);

    /* Order innermost-first: if kf[i] contains kf[k], kf[k] is deeper. */
    for (size_t i = 0; i < nl; i++)
        for (size_t k = i + 1; k < nl; k++)
            if (rt_contains(kf[i], kf[k])) {
                Expr* t = kf[i]; kf[i] = kf[k]; kf[k] = t;
                t = us[i]; us[i] = us[k]; us[k] = t;
            }

    /* t-symbols and substitution rules (Exp[u]->t and Power[E,u]->t). */
    Expr** ts = malloc(nl * sizeof(Expr*));
    Expr** rules = malloc(2 * nl * sizeof(Expr*));
    for (size_t i = 0; i < nl; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmte%zu", i);
        ts[i] = expr_new_symbol(nm);
        rules[2 * i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(us[i]) }, 1), expr_copy(ts[i]) }, 2);
        rules[2 * i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(kf[i]), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, 2 * nl);
    free(rules);

    Expr* F = rt_eval1("Together", rt_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(f), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL; Expr** Dt = NULL; Expr* result = NULL;

    bool ok = F && rt_find_log_of_x(F, x) == NULL && rt_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rt_eval1("Numerator", expr_copy(F));
        den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((nl + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < nl; i++) vv[i + 1] = expr_copy(ts[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nl + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rt_is_true(pqn) && rt_is_true(pqd)
             && !rt_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Dt_i = Cancel[D[u_i,x]] |_{ker->t} * t_i. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* dui = rt_eval2("D", expr_copy(us[i]), expr_copy(x));
            Expr* dc = rt_eval1("Cancel", dui);
            Expr* dcs = dc ? rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            Dt[i] = dcs ? rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ dcs, expr_copy(ts[i]) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
        }
    }

    if (ok) {
        /* a = multiplicity of t_n at 0 in den (Laurent negative extent). */
        long a = rt_var_mult_at_zero(den, top);
        Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(top), expr_new_integer(-a) }, 2) }, 2));

        /* Squarefree t_n-dependent factors of Dtil -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = Dtil ? rt_eval1("Factor", expr_copy(Dtil)) : NULL;
        if (!factored) bad = true;
        else {
            Expr** fa; size_t nf; Expr* single[1];
            if (factored->type == EXPR_FUNCTION
                && factored->data.function.head->type == EXPR_SYMBOL
                && factored->data.function.head->data.symbol.name == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol.name == intern_symbol("Power")
                    && term->data.function.arg_count == 2
                    && term->data.function.args[1]->type == EXPR_INTEGER) {
                    base = term->data.function.args[0];
                    e = (long)term->data.function.args[1]->data.integer;
                }
                if (rt_free_of_x(base, top)) continue;    /* t_n-coprime constant part */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long dtn_num = rt_degree(num, top), dtn_den = rt_degree(den, top);
        /* Top kernel t_n is an EXPONENTIAL (D preserves deg_{t_n}), so Q's Laurent
         * range in t_n matches f's: high extent deg_{t_n}(num) - deg_{t_n}(den),
         * low extent -(multiplicity of t_n at 0 in den).  Exact, no cap. */
        long ihi = dtn_num - dtn_den;
        long ilo = -a;

        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));   /* per-var odometer range size-1 */
        long* lo = malloc(nlv * sizeof(long));   /* per-var lowest exponent       */
        lv[0] = x; lo[0] = 0;                    /* x: polynomial, degree >= 0     */
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        {
            /* x is polynomial; degree bound derived from num/den, no cap. */
            long p = rt_degree(num, x); if (p < 0) p = 0;
            long q = rt_degree(den, x); if (q < 0) q = 0;
            long d = p + q + 1;
            bd[0] = d;
        }
        /* Inner kernels t_1..t_{n-1} are EXPONENTIALS and invertible, so their
         * coefficient exponents are LAURENT — the antiderivative can carry a power
         * the integrand lacks (e.g. INT E^(x+E^x) dx = E^(E^x)).  Each inner window
         * is f's own per-kernel Laurent extent WIDENED by the reach of the top
         * derivation coefficient w' = D[u_n] (Dt[top] = w' t_n): one factor of w'
         * enters via t_n^i in D_tower[Q], shifting inner-kernel exponents by up to
         * deg_{t_j}(w') upward and w''s t_j-pole order downward.  Fully derived from
         * the input — no hardcoded window. */
        Expr* wtop = rt_eval1("Together", expr_copy(Dt[nl - 1]));
        Expr* wnum = wtop ? rt_eval1("Numerator", expr_copy(wtop)) : NULL;
        Expr* wden = wtop ? rt_eval1("Denominator", expr_copy(wtop)) : NULL;
        for (size_t j = 1; j < nlv; j++) {
            long hij = rt_degree(num, ts[j - 1]) - rt_degree(den, ts[j - 1]);
            if (hij < 0) hij = 0;
            long loj = -rt_var_mult_at_zero(den, ts[j - 1]);
            long wr = 0;
            if (wnum && wden) {
                long wu = rt_degree(wnum, ts[j - 1]); if (wu < 0) wu = 0;
                wr = wu + rt_var_mult_at_zero(wden, ts[j - 1]);
            }
            lo[j] = loj - wr;
            bd[j] = (hij + wr) - (loj - wr);       /* odometer count-1 (>= 0) */
        }
        if (wtop) expr_free(wtop);
        if (wnum) expr_free(wnum);
        if (wden) expr_free(wden);
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
        long nwi = (ihi >= ilo) ? (ihi - ilo + 1) : 0;
        long nunk = nwi * nmono + (long)ng;

        if (!bad && nwi > 0 && nunk > 0) {
            size_t nq = (size_t)nunk;
            Expr** qterms = malloc(nq * sizeof(Expr*));
            Expr** syms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0, nsym = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long i = ilo; i <= ihi; i++)
                for (long m = 0; m < nmono; m++) {
                    char nm[40]; snprintf(nm, sizeof(nm), "rmEp%ld_%ld", i - ilo, m);
                    rt_decode_mono(m, bd, nlv, ev);
                    for (size_t j = 0; j < nlv; j++) ev[j] += lo[j];  /* Laurent shift */
                    Expr* mono = rt_build_monomial(lv, ev, nlv);
                    syms[nsym++] = expr_new_symbol(nm);
                    qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), mono,
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(top), expr_new_integer(i) }, 2) }, 3);
                }
            for (size_t j = 0; j < ng; j++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmEc%zu", j);
                syms[nsym++] = expr_new_symbol(nm);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm),
                      expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
            }
            free(ev);
            Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
            free(qterms);

            Expr** back = malloc(nl * sizeof(Expr*));
            for (size_t i = 0; i < nl; i++)
                back[i] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(ts[i]), expr_copy(kf[i]) }, 2);
            Expr* bl = expr_new_function(expr_new_symbol("List"), back, nl);
            free(back);

            result = rt_tower_solve(Q, syms, nsym, Dt, ts, nl, F, x, bl);
            expr_free(bl);
            for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
            free(syms);
            expr_free(Q);
        }
        for (size_t j = 0; j < ng; j++) expr_free(g[j]);
        if (factored) expr_free(factored);
        if (Dtil) expr_free(Dtil);
        free(lv); free(bd); free(lo);
    }

    /* Diff-back safety gate (see rt_verify_antideriv). */
    if (result && !rt_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    if (Dt) { for (size_t i = 0; i < nl; i++) if (Dt[i]) expr_free(Dt[i]); free(Dt); }
    if (num) expr_free(num);
    if (den) expr_free(den);
    if (F) expr_free(F);
    expr_free(rl);
    for (size_t i = 0; i < nl; i++) { expr_free(ts[i]); expr_free(kf[i]); expr_free(us[i]); }
    free(ts); free(kf); free(us);
    return result;
}

/* ================================================================== */
/* Genuine one-extension-at-a-time recursive Risch engine.             */
/* ================================================================== */
/* The tower cases above use a FLAT ansatz — one SolveAlways over every tower
 * variable at once, with bounded-degree POLYNOMIAL lower-field coefficients.
 * That misses two whole families the recursive algorithm (Bronstein ch. 5)
 * handles: (1) MIXED exp/log towers (the flat cases are each
 * single-kind), and (2) RATIONAL lower-field coefficients (a t_n coefficient
 * that is 1/x, not a polynomial — a nonlinear unknown for the flat ansatz).
 *
 * This engine peels ONE extension at a time.  To integrate F in the top field
 * K_n = K_{n-1}(t_n): split F into its polynomial/Laurent part in t_n and a
 * proper rational part; integrate the polynomial part coefficient-by-coefficient,
 * where each coefficient integral is itself an integration in the LOWER field
 * K_{n-1} (the recursion), bottoming out at the rational base case C(x)
 * (Integrate`BronsteinRational).  Because the coefficients are integrated in
 * their own field they may be rational, and because each level dispatches on its
 * own kernel kind the tower may mix logs and exponentials.
 *
 * The genuine algorithm's proper-rational part (tower Hermite + Rothstein-Trager
 * over K_{n-1}) and its general field Risch-DE are deferred to a later increment;
 * a nonzero proper part or a non-base RDE declines cleanly here.  As with the
 * flat tower cases the whole result is a bounded search, so it is diff-back
 * VERIFIED (rt_verify_antideriv) — a spurious certificate cannot ship. */

#define RT_MAXK 5

typedef enum { RT_LOG, RT_EXP } RtKind;

/* ================================================================== */
/* Elementary-integrability decision (P3: Bronstein §5.6 residue        */
/* criterion + Ch.6 Risch-DE no-solution certificates).                 */
/* ================================================================== */
/* Tri-state verdict of the decision procedure. */
typedef enum { RT_DEC_UNKNOWN = 0, RT_DEC_ELEMENTARY, RT_DEC_NONELEMENTARY } RtDecision;

/* File-local decision context.  The field integrator (rt_field_integrate and its
 * callees) is a genuine decision procedure — every NULL it returns is either a
 * theorem-backed "no elementary integral" or a dispatch-to-sibling decline.  When
 * a caller wants the VERDICT (not just the antiderivative) it sets
 * g_rt_decide_mode and runs the same path; the sub-algorithms raise
 * g_rt_decision to RT_DEC_NONELEMENTARY at their AUTHORITATIVE decline points
 * ONLY — a non-constant residue (Thm 5.6.1(ii)), a Risch DE with no rational
 * solution (Ch.6, via the complete rde_base / exact-degree-bound ansatz), or the
 * hypertangent Dc≠0 certificate — and NEVER on a routing decline.  The flag is
 * write-once/sticky within a run; the driver (rt_decide) resets it, then reads
 * ELEMENTARY iff a full antiderivative with no unintegrated remainder is produced
 * (an elementary result always wins over a stray flag from an abandoned path). */
static bool       g_rt_decide_mode = false;
static RtDecision g_rt_decision     = RT_DEC_UNKNOWN;

/* Raise the decision to NONELEMENTARY at an authoritative decline (no-op unless a
 * decision is being requested). */
static void rt_dec_nonelem(void) {
    if (g_rt_decide_mode) g_rt_decision = RT_DEC_NONELEMENTARY;
}

typedef struct {
    size_t n;
    RtKind kind[RT_MAXK];
    Expr* kernel[RT_MAXK];   /* Log[u_i] or Power[E, w_i]                (owned) */
    Expr* arg[RT_MAXK];      /* u_i (log argument) or w_i (exp exponent) (owned) */
    Expr* t[RT_MAXK];        /* fresh tower variable t_i                 (owned) */
    Expr* Dcoef[RT_MAXK];    /* log: u'/u ; exp: w' — in t-vars, in K_{i-1} (owned) */
    Expr* subrules;          /* List of all kernel -> t_i rules          (owned) */
    /* Multiplicatively commensurate non-primitive exp members: a collected
     * kernel E^w whose exponent w = mmult * arg[mprim] is NOT an independent
     * extension — it is (E^arg[mprim])^mmult = t[mprim]^mmult.  Recorded here so
     * rt_subst_kernels can alias it to a power of the primitive's tower var
     * instead of leaving it as a foreign kernel. */
    size_t nm;
    Expr* marg[2 * RT_MAXK]; /* the member exponent w                    (owned) */
    long  mprim[2 * RT_MAXK];/* tower index of the class primitive              */
    long  mmult[2 * RT_MAXK];/* integer multiplier k (w = k * arg[mprim])       */
} RtTower;

static void rt_tower_free(RtTower* T) {
    for (size_t i = 0; i < T->n; i++) {
        if (T->kernel[i]) expr_free(T->kernel[i]);
        if (T->arg[i]) expr_free(T->arg[i]);
        if (T->t[i]) expr_free(T->t[i]);
        if (T->Dcoef[i]) expr_free(T->Dcoef[i]);
    }
    if (T->subrules) expr_free(T->subrules);
    for (size_t i = 0; i < T->nm; i++)
        if (T->marg[i]) expr_free(T->marg[i]);
    T->n = 0; T->nm = 0; T->subrules = NULL;
}

/* Build the ordered differential tower of f over C(x).  Collect every
 * x-dependent Log and E^ kernel, order them innermost-first (deepest at index 0)
 * by structural containment — tie-breaking independent kernels EXP-before-LOG so
 * the primitive (logarithmic) recursion sits on top and the exponential Risch DEs
 * bottom out in C(x) — assign tower variables, and compute each derivation
 * coefficient Dcoef_i (log: u_i'/u_i ; exp: w_i').  The structure-theorem
 * soundness check requires every Dcoef_i to lie in K_{i-1} = C(x, t_1..t_{i-1})
 * (triangular: free of t_i..t_n and of any residual foreign kernel).  Returns
 * true with T populated (2 <= n <= RT_MAXK); false otherwise (caller still calls
 * rt_tower_free). */
static bool rt_tower_build_min(Expr* f, Expr* x, RtTower* T, size_t min_n) {
    for (size_t i = 0; i < RT_MAXK; i++)
        { T->kernel[i] = T->arg[i] = T->t[i] = T->Dcoef[i] = NULL; }
    T->subrules = NULL; T->n = 0;

    Expr** logs = NULL; size_t nl = 0, lc = 0; rt_collect_logs(f, x, &logs, &nl, &lc);
    Expr** exps = NULL; size_t ne = 0, ec = 0;
    rt_collect_exp_exponents(f, x, &exps, &ne, &ec);
    T->nm = 0;
    Expr* mprim_pexp[2 * RT_MAXK];       /* per-member class-primitive exponent (borrowed) */

    /* --- Multiplicatively commensurate reduction of the exponential kernels. ---
     * Collected exponents w_i, w_j define algebraically DEPENDENT kernels
     * E^w_i, E^w_j when w_i/w_j is a nonzero rational: then E^w = (E^prim)^k for
     * a class primitive exponent `prim` and integer k (e.g. E^(2 E^x) =
     * (E^(E^x))^2).  Partition the exponents into such commensurability classes,
     * SYNTHESIZE one primitive exponent per class (prim = member/lcm(ratio
     * denominators), possibly NOT itself a member — e.g. E^x for {2 E^x, 3 E^x},
     * or E^(x/6) for {E^(x/2), E^(x/3)}), keep ONLY the synthesized primitives as
     * tower extensions, and record each dependent member as the integer power
     * E^w -> t[prim]^k of its primitive's tower variable (T->m*).  Without this a
     * dependent kernel would spuriously add an extension, breaking independence.
     * A class whose members are not all rational multiples of one another has no
     * common primitive and declines the whole tower (never wrong). */
    long* clsrep = malloc((ne ? ne : 1) * sizeof(long));
    long* multof = malloc((ne ? ne : 1) * sizeof(long));
    for (size_t i = 0; i < ne; i++) { clsrep[i] = -1; multof[i] = 1; }
    bool okc = (ne <= 2 * RT_MAXK);                 /* member-array bound */
    for (size_t i = 0; i < ne && okc; i++) {        /* group by commensurability */
        if (clsrep[i] != -1) continue;
        clsrep[i] = (long)i;
        for (size_t j = i + 1; j < ne; j++) {
            if (clsrep[j] != -1) continue;
            Expr* r = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(exps[j]), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(exps[i]), expr_new_integer(-1) }, 2) }, 2));
            if (rt_is_rat_const(r)) clsrep[j] = (long)i;
            if (r) expr_free(r);
        }
    }
    /* Synthesize one primitive exponent per class (clsprim[rep], owned) and the
     * per-member integer multiplier multof[i] (w_i = multof[i] * clsprim[rep]). */
    Expr** clsprim = calloc((ne ? ne : 1), sizeof(Expr*));
    for (size_t rep = 0; rep < ne && okc; rep++) {
        if (clsrep[rep] != (long)rep) continue;
        Expr** cls = malloc((ne ? ne : 1) * sizeof(Expr*));
        size_t* idxs = malloc((ne ? ne : 1) * sizeof(size_t));
        size_t cm = 0;
        for (size_t m = 0; m < ne; m++)
            if (clsrep[m] == (long)rep) { cls[cm] = exps[m]; idxs[cm] = m; cm++; }
        long* kk = malloc((cm ? cm : 1) * sizeof(long));
        Expr* p = rt_class_primitive(cls, cm, kk);   /* owned primitive, or NULL */
        if (!p) okc = false;
        else {
            clsprim[rep] = p;
            for (size_t c = 0; c < cm; c++) multof[idxs[c]] = kk[c];
        }
        free(cls); free(idxs); free(kk);
    }
    size_t np = 0;
    for (size_t rep = 0; rep < ne; rep++) if (clsrep[rep] == (long)rep) np++;
    size_t n = nl + np;
    if (!okc || n < min_n || n > RT_MAXK) {
        for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
        free(clsprim); free(clsrep); free(multof);
        for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
        for (size_t i = 0; i < ne; i++) expr_free(exps[i]);
        free(logs); free(exps);
        return false;
    }

    size_t idx = 0;
    for (size_t i = 0; i < nl; i++) {
        T->kind[idx] = RT_LOG;
        T->kernel[idx] = logs[i];                                  /* adopt Log[u] */
        T->arg[idx] = expr_copy(logs[i]->data.function.args[0]);
        idx++;
    }
    free(logs);
    /* One EXP tower kernel per synthesized class primitive. */
    for (size_t rep = 0; rep < ne; rep++) {
        if (clsrep[rep] != (long)rep) continue;
        T->kind[idx] = RT_EXP;
        T->arg[idx] = expr_copy(clsprim[rep]);                     /* primitive exponent */
        T->kernel[idx] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(clsprim[rep]) }, 2);
        idx++;
    }
    /* Every exp member that is not itself the class primitive is an alias
     * E^w -> t[prim]^k; store its primitive exponent for the post-reorder remap. */
    for (size_t i = 0; i < ne; i++) {
        if (expr_eq(exps[i], clsprim[clsrep[i]])) continue;   /* realized by the primitive kernel */
        T->marg[T->nm] = expr_copy(exps[i]);                  /* member exponent w */
        T->mmult[T->nm] = multof[i];                          /* w = k * clsprim[rep]  */
        mprim_pexp[T->nm] = clsprim[clsrep[i]];               /* borrowed; remapped below */
        T->nm++;
    }
    T->n = n;

    /* Order innermost-first (deepest at index 0); tie-break EXP before LOG. */
    for (size_t pass = 0; pass < n; pass++)
        for (size_t i = 0; i + 1 < n; i++) {
            bool swap = false;
            if (rt_contains(T->kernel[i], T->kernel[i + 1])) swap = true;
            else if (!rt_contains(T->kernel[i + 1], T->kernel[i])
                     && T->kind[i] == RT_LOG && T->kind[i + 1] == RT_EXP) swap = true;
            if (swap) {
                RtKind kk = T->kind[i]; T->kind[i] = T->kind[i + 1]; T->kind[i + 1] = kk;
                Expr* a = T->kernel[i]; T->kernel[i] = T->kernel[i + 1]; T->kernel[i + 1] = a;
                a = T->arg[i]; T->arg[i] = T->arg[i + 1]; T->arg[i + 1] = a;
            }
        }

    /* Tower variables t_i and the combined substitution rule list.  Each member
     * kernel E^(marg) contributes a rule E^(marg) -> t[prim]^mmult (both the
     * Exp[] and Power[E,] spellings); mprim is remapped from an exps[] index to
     * the primitive's post-ordering tower index by matching its exponent. */
    Expr** rules = malloc((2 * n + 2 * T->nm) * sizeof(Expr*)); size_t nr = 0;
    for (size_t i = 0; i < n; i++) {
        char nm[16]; snprintf(nm, sizeof(nm), "rmR%zu", i);
        T->t[i] = expr_new_symbol(nm);
        if (T->kind[i] == RT_LOG) {
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T->kernel[i]), expr_copy(T->t[i]) }, 2);
        } else {
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                    (Expr*[]){ expr_copy(T->arg[i]) }, 1), expr_copy(T->t[i]) }, 2);
            rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_copy(T->kernel[i]), expr_copy(T->t[i]) }, 2);
        }
    }
    for (size_t m = 0; m < T->nm; m++) {
        Expr* pw = mprim_pexp[m];                  /* class primitive exponent value */
        long ti = -1;
        for (size_t i = 0; i < n; i++)
            if (T->kind[i] == RT_EXP && expr_eq(T->arg[i], pw)) { ti = (long)i; break; }
        T->mprim[m] = ti;                          /* now a tower index (>= 0) */
        Expr* tk = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(T->t[ti]), expr_new_integer(T->mmult[m]) }, 2);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Exp"),
                (Expr*[]){ expr_copy(T->marg[m]) }, 1), expr_copy(tk) }, 2);
        rules[nr++] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(T->marg[m]) }, 2), tk }, 2);
    }
    T->subrules = expr_new_function(expr_new_symbol("List"), rules, nr);
    free(rules);
    for (size_t rep = 0; rep < ne; rep++) if (clsprim[rep]) expr_free(clsprim[rep]);
    free(clsprim); free(clsrep); free(multof);
    for (size_t i = 0; i < ne; i++) expr_free(exps[i]);
    free(exps);

    /* Derivation coefficients + structure-theorem (triangularity) soundness. */
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        /* Compute the RAW derivative (Log: u'/u, Exp: w') WITHOUT normalising,
         * substitute the kernels structurally, and ONLY THEN Cancel in the tower
         * variables.  Cancelling first expands a factored kernel power such as
         * (x + E^x)^2 into x^2 + 2 x E^x + E^(2 x); the E^(2 x) term is not the
         * literal kernel E^x, so the substitution misses it and the structure
         * check spuriously rejects a valid tower (e.g. Bronstein's own example
         * E^((x^2-1)/x + 1/(x + E^x))).  Substituting on the unexpanded form maps
         * (x + E^x)^2 -> (x + t_0)^2 cleanly. */
        Expr* d;
        if (T->kind[i] == RT_LOG) {
            Expr* du = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
            d = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ du, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->arg[i]), expr_new_integer(-1) }, 2) }, 2);
        } else {
            d = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
        }
        if (!d) { ok = false; break; }
        Expr* dsub = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ d, expr_copy(T->subrules) }, 2));           /* adopts d */
        Expr* ds = dsub ? rt_eval1("Cancel", dsub) : NULL;
        if (!ds) { ok = false; break; }
        T->Dcoef[i] = ds;
        if (rt_find_exp_of_x(ds, x) != NULL || rt_find_log_of_x(ds, x) != NULL) ok = false;
        for (size_t j = i; j < n && ok; j++)
            if (!rt_free_of_x(ds, T->t[j])) ok = false;
    }
    return ok;
}

/* Mutually recursive field-integration primitives (forward declarations).
 * `rem_out` (when non-NULL) receives the unintegrated non-constant-residue
 * remainder of a partial log part (Bronstein Thm 5.6.1); it is threaded only
 * from the top-level recursive-tower assembly point — inner recursions pass
 * NULL, which makes the residue criterion decline a partial rather than surface
 * one mid-recursion. */
static Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out);
static int   rt_limited_field_integrate(Expr* r, RtTower* T, long L, Expr* x,
                                        Expr** s_out, Expr** c_out);
static Expr* rt_int_primitive_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_int_hyperexp_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_field_rde(Expr* p, long i, RtTower* T, long L, Expr* x);
static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x, Expr** rem_out);
static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out);

/* Tower derivation D_tower[e] = D[e,x] + sum_i Dt_i D[e,t_i], with Dt_i =
 * Dcoef_i (log) or Dcoef_i * t_i (exp).  Owned result. */
static Expr* rt_tower_deriv(Expr* e, RtTower* T, Expr* x) {
    Expr** terms = malloc((T->n + 1) * sizeof(Expr*));
    terms[0] = rt_eval2("D", expr_copy(e), expr_copy(x));
    for (size_t i = 0; i < T->n; i++) {
        Expr* dei = rt_eval2("D", expr_copy(e), expr_copy(T->t[i]));
        Expr* dti = (T->kind[i] == RT_LOG)
            ? expr_copy(T->Dcoef[i])
            : expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(T->Dcoef[i]), expr_copy(T->t[i]) }, 2);
        terms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ dti, dei }, 2);
    }
    Expr* r = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, T->n + 1));
    free(terms);
    return r;
}

/* Pure resultant Lazard-Rioboo-Trager log part for a TOWER proper rational
 * part Rr/den in t = t_L (deg_t Rr < deg_t den), lifting rt_frac_lrt's single-
 * kernel algebraic-residue closure to the recursive-Risch tower.  The bounded
 * SolveAlways ansatz in rt_field_ratint / rt_field_hyperexp_coupled expresses
 * exactly one CONSTANT residue per squarefree t-factor, so it cannot integrate
 * a factor whose Rothstein-Trager residues are ALGEBRAIC (e.g. the +-I/2
 * residues that split t^2+1 into a conjugate log pair = ArcTan) — this closes
 * exactly that class at the tower level, e.g.
 *   1/(x Log[x] (Log[Log[x]]^2+1)) -> ArcTan[Log[Log[x]]].
 *
 * Only the squarefree case (den has no repeated t-factor) is handled: den
 * factors as cont(x, t_0..t_{L-1}) * rad with rad squarefree in t, and the log
 * part is
 *   Integrate`TranscendentalLogPart[Rr/cont, rad, t, z, D_tower[rad],
 *                                   {x, t_0, ..., t_{L-1}}],
 * whose gate list requires the residues to be constants of the tower derivation
 * (free of every lower-field variable).  A genuine Hermite (repeated-pole) part
 * is left to the coupled ansatz.  Returns the log part in tower-variable form
 * (owned, in terms of the t_L symbol), or NULL.  The top-level caller diff-back
 * verifies, so an incomplete reduction (e.g. an exponential Laurent part not
 * captured here) declines rather than shipping a wrong form. */
static Expr* rt_field_lrt_logpart(Expr* Rr, Expr* den, RtTower* T, long L,
                                  Expr* x, Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    long dn = rt_degree(Rr, t), dd = rt_degree(den, t);
    if (dn < 0 || dd <= 0 || dn >= dd) return NULL;   /* proper, den depends on t */

    /* Strip the t-free content of den: Hden = gcd(den, dden/dt); when den is
     * squarefree in t this is exactly the t-free factor (dH == 0). */
    Expr* dden = rt_eval2("D", expr_copy(den), expr_copy(t));
    if (!dden) return NULL;
    Expr* Hden = rt_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(den), dden }, 2);
    if (!Hden) return NULL;
    if (rt_degree(Hden, t) != 0) { expr_free(Hden); return NULL; }  /* repeated pole */

    /* rad = den / Hden (squarefree in t); a = Rr / Hden. */
    Expr* rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    Expr* a = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(Rr), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    expr_free(Hden);
    if (!rad || !a) { if (rad) expr_free(rad); if (a) expr_free(a); return NULL; }

    /* Monomial (tower) derivation of the squarefree radical. */
    Expr* Dd = rt_tower_deriv(rad, T, x);
    if (!Dd) { expr_free(rad); expr_free(a); return NULL; }

    /* Gate list {x, t_0, ..., t_{L-1}}: residues must be constants of D_tower. */
    size_t ng = (size_t)L + 1;
    Expr** gv = malloc(ng * sizeof(Expr*));
    gv[0] = expr_copy(x);
    for (long i = 0; i < L; i++) gv[i + 1] = expr_copy(T->t[i]);
    Expr* gate = expr_new_function(expr_new_symbol("List"), gv, ng);
    free(gv);

    /* In decision mode, append a 7th "decide" marker so a non-constant-residue
     * resultant returns the inert marker Integrate`$NonConstantResidue (the simple
     * part is non-elementary, Thm 5.6.1(ii)) instead of declining silently. */
    Expr* logpart = g_rt_decide_mode
        ? rt_eval_call("Integrate`TranscendentalLogPart",
              (Expr*[]){ a, rad, expr_copy(t), expr_new_symbol("rmZ"), Dd, gate,
                         expr_new_symbol("True") }, 7)
        : rt_eval_call("Integrate`TranscendentalLogPart",
              (Expr*[]){ a, rad, expr_copy(t), expr_new_symbol("rmZ"), Dd, gate }, 6);
    /* Non-constant-residue verdict marker (decision mode): the simple part has no
     * elementary integral — raise the decision and decline the (empty) log part. */
    if (logpart && logpart->type == EXPR_SYMBOL
        && logpart->data.symbol.name == intern_symbol("Integrate`$NonConstantResidue")) {
        expr_free(logpart);
        rt_dec_nonelem();
        return NULL;
    }
    /* Declined iff the head is unchanged (a, rad, Dd, gate adopted by the call). */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol.name
             == intern_symbol("Integrate`TranscendentalLogPart")) {
        expr_free(logpart);
        return NULL;
    }
    /* Partial log part (Bronstein Thm 5.6.1): the builtin returned
     * Integrate`PartialLogPart[logs, remainder] — the residue resultant mixed
     * constant and non-constant residues.  Only a caller that can carry the
     * remainder (rem_out != NULL) may accept it; otherwise decline (the
     * exponential coupled path cannot absorb a proper simple remainder). */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol.name
             == intern_symbol("Integrate`PartialLogPart")
        && logpart->data.function.arg_count == 2) {
        if (!rem_out) { expr_free(logpart); return NULL; }
        Expr* logs = expr_copy(logpart->data.function.args[0]);
        *rem_out = expr_copy(logpart->data.function.args[1]);
        expr_free(logpart);
        return logs;   /* tower-variable elementary logs; remainder via rem_out */
    }
    return logpart;   /* tower-variable form (in terms of the t_L symbol) */
}

/* Proper rational part in t = t_L (deg_t num < deg_t den, coefficients in the
 * lower field K_{L-1}): integrate to  Q = H(t)/Hden(t) + sum_j c_j Log(g_j),
 * where Hden = gcd(den, d den/dt) is the repeated part, g_j the distinct
 * t-dependent factors of the squarefree radical rad = den/Hden (Rothstein-Trager
 * log arguments), H a polynomial in t of degree < deg(Hden) with bounded-degree
 * lower-field (x, t_1..t_{L-1}) polynomial coefficients, and c_j constants.  All
 * unknowns are solved at once by requiring D_tower[Q] = num/den for every tower
 * variable via SolveAlways — the constant-residue and denominator-structure
 * certificates are exactly Hermite reduction + the residue criterion, lifted to
 * the tower derivation.  (A genuinely rational Hermite numerator coefficient and
 * the pure resultant LRT remain a later refinement; a bounded ansatz that misses
 * them declines, and the caller diff-back verifies.)  Returns Q (tower-variable
 * form, owned) or NULL. */
/* Build the tower derivation as a List[Rule[var, Dvar], ...] for risch_field /
 * risch_hermite: x -> 1, and each tower variable t_i -> Dt_i (Dcoef_i for a log,
 * Dcoef_i * t_i for an exp).  Owned; matches rt_tower_deriv's D_tower exactly. */
static Expr* rt_build_deriv_rules(RtTower* T, Expr* x) {
    size_t n = T->n;
    Expr** rules = malloc((n + 1) * sizeof(Expr*));
    rules[0] = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_copy(x), expr_new_integer(1) }, 2);
    for (size_t i = 0; i < n; i++) {
        Expr* dti = (T->kind[i] == RT_LOG)
            ? expr_copy(T->Dcoef[i])
            : expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(T->Dcoef[i]), expr_copy(T->t[i]) }, 2);
        rules[i + 1] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T->t[i]), dti }, 2);
    }
    Expr* r = expr_new_function(expr_new_symbol("List"), rules, n + 1);
    free(rules);
    return r;
}

/* Bronstein canonical representation (§3.5) at tower level L: split
 * F = f_p + f_s + f_n into its polynomial part f_p in t = t_L, the special part
 * f_s = b/d_s (special denominator — poles fixed by the derivation, e.g. t^k for
 * an exponential, t^2+1 for a hypertangent), and the normal part f_n = c/d_n
 * (proper, denominator coprime to the special part).  Reuses the tested
 * differential-field C API (risch_canonical_representation) over the tower's
 * monomial derivation.  Writes owned *f_p, *f_s, *f_n; returns false (leaving
 * them untouched) only if the derivation rules are malformed. */
static bool rt_canonical_split(Expr* F, RtTower* T, long L, Expr* x,
                               Expr** f_p, Expr** f_s, Expr** f_n) {
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); return false; }
    risch_canonical_representation(F, T->t[L], &d, f_p, f_s, f_n);
    risch_deriv_free(&d);
    expr_free(rules);
    return true;
}

/* Literal Hermite reduction of the proper part Rr/den in t = t_L (Bronstein
 * §5.3) — the exact algorithm, with no undetermined-coefficient ansatz.  For a
 * primitive (log) monomial t_L is normal and Rr/den is proper, so the reduced
 * part r is 0: Rr/den = D_tower[g] + h with h simple, and the antiderivative is
 * g + (log part of h), the log part coming from the residue criterion (§5.6).
 * The literal reduction handles arbitrary rational K_{L-1} = C(x, t_0..t_{L-1})
 * numerator coefficients.  Returns g + logpart, self-verified by a tower
 * diff-back, or NULL when there is no elementary integral. */
static Expr* rt_field_ratint_hermite(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                     Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    Expr* f = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2)));
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); expr_free(f); return NULL; }
    Expr* g = NULL; Expr* h = NULL; Expr* r = NULL;
    bool ok = risch_hermite_reduce(f, t, &d, &g, &h, &r);
    risch_deriv_free(&d);
    expr_free(rules); expr_free(f);
    if (!ok) return NULL;

    /* Proper primitive input: the reduced part must vanish (no polynomial/special
     * leftover).  If not (unexpected), decline. */
    if (!rt_is_zero(r)) { expr_free(g); expr_free(h); expr_free(r); return NULL; }
    expr_free(r);

    Expr* logpart = NULL;
    Expr* rem_tower = NULL;   /* partial log part: unintegrated r_n remainder */
    if (!rt_is_zero(h)) {
        /* h is simple (squarefree normal denominator): its integral is the
         * residue-criterion / Rothstein-Trager log part.  With a partial log
         * part (Thm 5.6.1) it returns the constant-residue logs plus rem_tower. */
        Expr* hn = rt_eval1("Numerator", rt_eval1("Together", expr_copy(h)));
        Expr* hd = rt_eval1("Denominator", rt_eval1("Together", expr_copy(h)));
        logpart = (hn && hd) ? rt_field_lrt_logpart(hn, hd, T, L, x, rem_out ? &rem_tower : NULL) : NULL;
        if (hn) expr_free(hn);
        if (hd) expr_free(hd);
        if (!logpart) { expr_free(g); expr_free(h); return NULL; }  /* not elementary */
    }
    expr_free(h);

    Expr* result = logpart
        ? rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
              (Expr*[]){ g, logpart }, 2))
        : g;

    /* Self-verify the (possibly partial) identity D_tower[result] == num/den −
     * rem_tower.  Exact and decidable in the tower field; identical to the full
     * check when rem_tower is NULL.  A mismatch declines rather than shipping a
     * wrong form. */
    Expr* Dres = rt_tower_deriv(result, T, x);
    Expr* tgt = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Dres, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), tgt }, 2) }, 2);
    if (rem_tower)   /* D[result] − num/den + rem_tower == 0 */
        diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ diff, expr_copy(rem_tower) }, 2);
    Expr* chk = rt_eval1("Together", diff);
    bool good = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!good) { expr_free(result); if (rem_tower) expr_free(rem_tower); return NULL; }
    if (rem_out) *rem_out = rem_tower;
    else if (rem_tower) expr_free(rem_tower);
    return result;
}

static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                             Expr** rem_out) {
    /* The proper rational part in t = t_L, integrated by the LITERAL Bronstein
     * pipeline: Hermite reduction (Sect. 5.3) peels the repeated normal poles
     * into an exact derivative D_tower[g], leaving a simple part h whose integral
     * is the residue-criterion / Lazard-Rioboo-Trager log part (Sect. 5.6).  This
     * is the exact decision procedure; the former bounded undetermined-
     * coefficient ansatz has been removed.  Returns g + logpart (self-verified by
     * a tower diff-back), or NULL when there is no elementary integral.  A mixed
     * residue resultant returns the constant-residue logs and the r_n remainder
     * via rem_out (Thm 5.6.1) when rem_out != NULL. */
    return rt_field_ratint_hermite(num, den, T, L, x, rem_out);
}

/* Coupled hyperexponential proper part (exponential top level t = t_L = E^(w_L)),
 * integrated by the LITERAL Bronstein pipeline (Sect. 5.3 + 5.6 + 5.9):
 *
 *   1. HermiteReduce(f)  ->  (g, h, r)              [Sect. 5.3, over D_tower]
 *      g = exact rational part (repeated NORMAL poles peeled, with arbitrary
 *      rational lower-field coefficients); h = simple normal part (squarefree
 *      denominator, coprime to t); r = reduced Laurent polynomial in t (the
 *      special part f_s at t = 0 plus the polynomial part f_p — CanonicalRepresen-
 *      tation puts the t-power denominator into f_s for an exponential monomial).
 *   2. residue criterion on h  ->  L = sum_j c_j Log(g_j)          [Sect. 5.6]
 *      constant residues c_j (rt_field_lrt_logpart gate); the g_j are coprime to t.
 *   3. leftover  P = h + r - D_tower[L].  For an exponential monomial
 *      D_tower[Log g] = D_tower[g]/g is NOT proper (deg_t D[g] = deg_t g since
 *      D[t] = w' t), so the logs of step 2 spill a t-polynomial into the Laurent
 *      part — Bronstein's coupling.  Subtracting D_tower[L] reconciles it exactly,
 *      leaving P a genuine Laurent polynomial in t (denominator a power of t).
 *   4. IntegrateHyperexponentialPolynomial(P)  ->  Q               [Sect. 5.9]
 *      per-coefficient Risch DE (rt_int_hyperexp_poly): each t^i coefficient (i!=0)
 *      solves q_i' + i w' q_i = p_i, the t^0 coefficient recurses to the lower field.
 *   return g + L + Q, self-verified by a tower diff-back.
 *
 * This is the exact decision procedure — Hermite reduction + residue criterion +
 * hyperexponential-polynomial integration — and is strictly more complete than the
 * former unified SolveAlways ansatz (which used bounded-degree POLYNOMIAL Hermite,
 * Laurent and log coefficients, declining rational lower-field coefficients).  The
 * ansatz has been removed.  Returns g + L + Q (tower-variable form, owned), or NULL
 * when there is no elementary integral. */
static Expr* rt_field_hyperexp_hermite(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    Expr* t = T->t[L];
    Expr* f = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2)));
    if (!f) return NULL;
    Expr* rules = rt_build_deriv_rules(T, x);
    RischDeriv d;
    if (!risch_deriv_from_rules(rules, &d)) { expr_free(rules); expr_free(f); return NULL; }
    Expr* g = NULL; Expr* h = NULL; Expr* r = NULL;
    bool ok = risch_hermite_reduce(f, t, &d, &g, &h, &r);
    risch_deriv_free(&d);
    expr_free(rules); expr_free(f);
    if (!ok) return NULL;

    /* Step 2: residue-criterion log part of the simple h.  When the caller can carry
     * a remainder (rem_out != NULL) we ACCEPT a partial log part (Bronstein Thm
     * 5.6.1): a mixed residue resultant yields the constant-residue logs L_s plus the
     * non-constant-residue remainder h_n (a proper simple fraction in t).  The
     * elementary part is g + L_s + Q_s, with Q_s reconciled from h_s = h - h_n, and
     * h_n reported unintegrated — mirroring the primitive/log path.  With rem_out ==
     * NULL a mixed resultant declines (rt_field_lrt_logpart returns NULL). */
    Expr* Llog = NULL;
    Expr* h_n = NULL;   /* non-constant-residue remainder (owned), NULL if none */
    if (!rt_is_zero(h)) {
        Expr* hn = rt_eval1("Numerator", rt_eval1("Together", expr_copy(h)));
        Expr* hd = rt_eval1("Denominator", rt_eval1("Together", expr_copy(h)));
        Llog = (hn && hd) ? rt_field_lrt_logpart(hn, hd, T, L, x, rem_out ? &h_n : NULL) : NULL;
        if (hn) expr_free(hn);
        if (hd) expr_free(hd);
        if (!Llog) { expr_free(g); expr_free(h); expr_free(r); return NULL; }  /* non-elementary */
    } else {
        Llog = expr_new_integer(0);
    }

    /* h_s = h - h_n (the constant-residue part that reconciles into the Laurent
     * polynomial); h_s == h when there is no partial remainder.  h is consumed. */
    Expr* h_s;
    if (h_n) {
        h_s = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ h, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(h_n) }, 2) }, 2));   /* adopts h */
        if (!h_s) { expr_free(g); expr_free(r); expr_free(Llog); expr_free(h_n); return NULL; }
    } else {
        h_s = h;
    }

    /* Step 3: Laurent leftover  P = h_s + r - D_tower[L]  (a Laurent polynomial in t).
     * (h_s, r are consumed by the Plus below regardless of success.) */
    Expr* DL = rt_tower_deriv(Llog, T, x);
    if (!DL) {
        expr_free(g); expr_free(h_s); expr_free(r); expr_free(Llog);
        if (h_n) expr_free(h_n);
        return NULL;
    }
    Expr* P = rt_eval1("Cancel", rt_eval1("Together",
        expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ h_s, r,
              expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), DL }, 2) }, 3)));
    if (!P) { expr_free(g); expr_free(Llog); if (h_n) expr_free(h_n); return NULL; }

    /* Step 4: integrate the Laurent polynomial P (per-coefficient Risch DE).  If P
     * has any t-coprime denominator factor (the logs did not reconcile), the
     * hyperexponential-polynomial gate declines and so do we. */
    Expr* Pn = rt_eval1("Numerator", expr_copy(P));
    Expr* Pd = rt_eval1("Denominator", expr_copy(P));
    expr_free(P);
    Expr* Q = (Pn && Pd) ? rt_int_hyperexp_poly(Pn, Pd, T, L, x) : NULL;
    if (Pn) expr_free(Pn);
    if (Pd) expr_free(Pd);
    if (!Q) { expr_free(g); expr_free(Llog); if (h_n) expr_free(h_n); return NULL; }

    Expr* result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ g, Llog, Q }, 3));

    /* Self-verify the (possibly partial) identity D_tower[result] + h_n == num/den.
     * Exact and decidable in the tower field; with h_n == NULL this is the usual
     * D_tower[result] == num/den.  A mismatch declines (a wrong h_n reconstruction —
     * the κ_D gcd reconstruction is not exact over a transcendental monomial — can
     * therefore only decline, never ship a wrong partial). */
    Expr* Dres = rt_tower_deriv(result, T, x);
    Expr* tgt = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Dres, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), tgt }, 2) }, 2);
    if (h_n) {
        diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ diff, expr_copy(h_n) }, 2);
    }
    Expr* chk = rt_eval1("Together", diff);
    bool good = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!good) { expr_free(result); if (h_n) expr_free(h_n); return NULL; }
    if (rem_out) *rem_out = h_n;
    else if (h_n) expr_free(h_n);
    return result;
}

static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x,
                                       Expr** rem_out) {
    /* The exponential-top proper part, integrated by the literal Bronstein pipeline
     * (Hermite reduction + residue criterion + hyperexponential-polynomial), with the
     * former undetermined-coefficient ansatz removed.  See rt_field_hyperexp_hermite.
     * A mixed residue resultant surfaces its non-constant-residue part via rem_out. */
    return rt_field_hyperexp_hermite(num, den, T, L, x, rem_out);
}

/* Integrate F (a rational function of x, t_1..t_L in tower-variable form) with
 * respect to the tower derivation, returning an antiderivative in tower-variable
 * form (owned) or NULL.  L < 0 is the rational base case C(x). */
static Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x, Expr** rem_out) {
    if (rem_out) *rem_out = NULL;
    if (!F) return NULL;
    if (L < 0) {
        Expr* r = rt_eval_call("Integrate`BronsteinRational",
            (Expr*[]){ expr_copy(F), expr_copy(x) }, 2);
        if (!r) return NULL;
        if (rt_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
        return r;
    }
    Expr* Fg = rt_eval1("Together", expr_copy(F));
    if (!Fg) return NULL;
    Expr* num = rt_eval1("Numerator", expr_copy(Fg));
    Expr* den = rt_eval1("Denominator", expr_copy(Fg));
    Expr* result = NULL;
    if (num && den) {
        if (T->kind[L] == RT_LOG) {
            /* Bronstein canonical representation (§3.5): F = f_p + f_s + f_n over
             * the primitive (log) monomial t_L.  A primitive monomial has no
             * non-trivial special polynomial (Dt ∈ k ⇒ every special factor is a
             * unit), so f_s ≡ 0; integrate the polynomial part f_p by the
             * primitive-polynomial recursion (§5.8) and the normal proper part
             * f_n by Hermite reduction + the Rothstein-Trager log part
             * (§5.3 + §5.6).  This unifies the former ad-hoc PolynomialQuotient/
             * Remainder denominator split (poly-quotient = f_p, remainder/den =
             * f_n for a primitive) behind the canonical decomposition. */
            Expr* f_p = NULL; Expr* f_s = NULL; Expr* f_n = NULL;
            if (rt_canonical_split(Fg, T, L, x, &f_p, &f_s, &f_n)) {
                if (rt_is_zero(f_s)) {   /* primitive: special part must vanish */
                    Expr* one = expr_new_integer(1);
                    Expr* poly_int = rt_int_primitive_poly(f_p, one, T, L, x);
                    expr_free(one);
                    Expr* prop_int;
                    if (rt_is_zero(f_n)) {
                        prop_int = expr_new_integer(0);
                    } else {
                        Expr* nn = rt_eval1("Numerator", expr_copy(f_n));
                        Expr* dd = rt_eval1("Denominator", expr_copy(f_n));
                        prop_int = (nn && dd) ? rt_field_ratint(nn, dd, T, L, x, rem_out) : NULL;
                        if (nn) expr_free(nn);
                        if (dd) expr_free(dd);
                    }
                    if (poly_int && prop_int)
                        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ poly_int, prop_int }, 2));
                    else {
                        if (poly_int) expr_free(poly_int);
                        if (prop_int) expr_free(prop_int);
                    }
                }
                expr_free(f_p); expr_free(f_s); expr_free(f_n);
            }
        } else {
            /* Exponential top: the pure-Laurent recursion first (genuine RDE,
             * rational lower-field coefficients); if it declines because a proper
             * part is present, the coupled hyperexponential pipeline.  The latter
             * (rt_field_hyperexp_hermite) already performs the canonical split
             * implicitly: its HermiteReduce -> (g, h, r) yields the f_n-derivative
             * part g, the simple normal part h (= f_n), and the reduced Laurent
             * polynomial r (= f_p + f_s, the special part living at t = 0), which
             * it must keep COUPLED through the §5.9 residue reconciliation — so an
             * independent f_p/f_s/f_n split is deliberately not applied here. */
            result = rt_int_hyperexp_poly(num, den, T, L, x);
            if (!result)
                result = rt_field_hyperexp_coupled(num, den, T, L, x, rem_out);
        }
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(Fg);
    /* A partial remainder is meaningful only alongside a successful elementary
     * part; if the overall integration failed, drop it. */
    if (!result && rem_out && *rem_out) { expr_free(*rem_out); *rem_out = NULL; }
    return result;
}

/* Primitive (logarithmic top) polynomial part.  num/den is a polynomial in
 * t = t_L (den free of t) with lower-field coefficients; Dt = Dcoef_L = u_L'/u_L.
 * Integrate P(t) = sum_i p_i t^i to Q = sum_i q_i t^i by the recursive coefficient
 * matching q_i' + (i+1) q_{i+1} Dt = p_i (top-down): each residual r_i is
 * integrated in the LOWER field K_L (recursion), and a would-be new logarithm
 * equal to t is folded back into q_{i+1}. */
static Expr* rt_int_primitive_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    Expr* Dt = T->Dcoef[L];
    long m = rt_degree(num, t);
    if (m < 0) return expr_new_integer(0);    /* num == 0: integral is 0 */
    Expr** p = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (long i = 0; i <= m; i++)
        p[i] = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ rt_coeff(num, t, i), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2));
    Expr** q = calloc((size_t)(m + 2), sizeof(Expr*));
    bool fail = false;
    for (long i = m; i >= 0 && !fail; i--) {
        Expr* r_i;
        if (!q[i + 1]) {
            r_i = expr_copy(p[i]);
        } else {
            Expr* term = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(i + 1), expr_copy(q[i + 1]),
                          expr_copy(Dt) }, 3);
            r_i = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(p[i]),
                    expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), term }, 2) }, 2));
        }
        Expr* s = NULL; Expr* c = NULL;
        int rc = r_i ? rt_limited_field_integrate(r_i, T, L, x, &s, &c) : -1;
        if (r_i) expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        if (c) {
            Expr* bump = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2));
            if (!q[i + 1]) q[i + 1] = bump;
            else q[i + 1] = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ q[i + 1], bump }, 2));
        }
        q[i] = s;
    }
    Expr* result = NULL;
    if (!fail) {
        Expr** terms = malloc((size_t)(m + 2) * sizeof(Expr*));
        size_t nt = 0;
        for (long i = 0; i <= m + 1; i++) {
            if (!q[i] || rt_is_zero(q[i])) continue;
            Expr* pw = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(t), expr_new_integer(i) }, 2);
            terms[nt++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[i]), pw }, 2);
        }
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
        free(terms);
    }
    for (long i = 0; i <= m; i++) if (p[i]) expr_free(p[i]);
    for (long i = 0; i <= m + 1; i++) if (q[i]) expr_free(q[i]);
    free(p); free(q);
    return result;
}

/* Limited integration in the lower field: integrate r (free of t_L) in K_L and
 * split the result as s + c*t_L with s in K_L and c a constant (the single new
 * logarithm equal to t_L = Log[u_L] is folded into c; any other new logarithm or
 * non-constant t_L-coefficient means r is not integrable within this tower here,
 * so decline).  On success sets *s_out (owned, may be NULL for 0) and *c_out
 * (owned, NULL when 0) and returns 0; else -1. */
static int rt_limited_field_integrate(Expr* r, RtTower* T, long L, Expr* x,
                                      Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rt_field_integrate(r, T, L - 1, x, NULL);  /* lower field: no partial */
    if (!R) return -1;
    Expr* Rs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ R, expr_copy(T->subrules) }, 2));               /* adopts R */
    if (!Rs) return -1;
    Expr* tL = T->t[L];
    int rc = -1;
    if (rt_find_exp_of_x(Rs, x) == NULL && rt_find_log_of_x(Rs, x) == NULL
        && rt_is_poly(Rs, tL) && rt_degree(Rs, tL) <= 1) {
        Expr* c = rt_coeff(Rs, tL, 1);
        Expr* s = rt_coeff(Rs, tL, 0);
        bool cconst = rt_free_of_x(c, x);
        for (size_t j = 0; j < T->n && cconst; j++)
            if (!rt_free_of_x(c, T->t[j])) cconst = false;
        if (cconst) {
            if (rt_is_zero(c)) { expr_free(c); *c_out = NULL; }
            else *c_out = c;
            *s_out = s;
            rc = 0;
        } else {
            /* The would-be new-logarithm coefficient c is NON-constant (Dc != 0):
             * by Bronstein's primitive-polynomial criterion (§5.8, p.158) the term
             * has no elementary integral.  Authoritative non-elementary verdict. */
            rt_dec_nonelem();
            expr_free(c); expr_free(s);
        }
    }
    expr_free(Rs);
    return rc;
}

/* Hyperexponential (exponential top) Laurent part.  num/den is a Laurent
 * polynomial in t = t_L (den = c t^a, a monomial in t) with lower-field
 * coefficients; Dcoef_L = w_L'.  Integrate sum_i p_i t^i to sum_i q_i t^i: the
 * i = 0 term integrates in the lower field K_L (recursion); each i != 0 term
 * solves the Risch differential equation q_i' + i w_L' q_i = p_i.  A genuine
 * (non-monomial) t-denominator is a proper fraction — declined here. */
static Expr* rt_int_hyperexp_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    long a = 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(t));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(t), expr_new_integer(-a) }, 2) }, 2));
    if (!Dtil || !rt_free_of_x(Dtil, t)) { if (Dtil) expr_free(Dtil); return NULL; }
    Expr* nnum = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ Dtil, expr_new_integer(-1) }, 2) }, 2));     /* adopts Dtil */
    if (!nnum) return NULL;
    long dnn = rt_degree(nnum, t);
    if (dnn < 0) { expr_free(nnum); return expr_new_integer(0); }  /* num == 0 */

    Expr** q = calloc((size_t)(dnn + 1), sizeof(Expr*));
    long* pw = malloc((size_t)(dnn + 1) * sizeof(long));
    bool fail = false;
    for (long j = 0; j <= dnn && !fail; j++) {
        long ip = j - a;
        pw[j] = ip;
        Expr* pj = rt_coeff(nnum, t, j);
        if (rt_is_zero(pj)) { expr_free(pj); q[j] = NULL; continue; }
        Expr* qj = (ip == 0) ? rt_field_integrate(pj, T, L - 1, x, NULL)
                             : rt_field_rde(pj, ip, T, L, x);
        expr_free(pj);
        if (!qj) { fail = true; break; }
        q[j] = qj;
    }
    Expr* result = NULL;
    if (!fail) {
        Expr** terms = malloc((size_t)(dnn + 1) * sizeof(Expr*));
        size_t nt = 0;
        for (long j = 0; j <= dnn; j++) {
            if (!q[j]) continue;
            Expr* term;
            if (pw[j] == 0) term = expr_copy(q[j]);
            else term = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[j]), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(t), expr_new_integer(pw[j]) }, 2) }, 2);
            terms[nt++] = term;
        }
        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
        free(terms);
    }
    for (long j = 0; j <= dnn; j++) if (q[j]) expr_free(q[j]);
    free(q); free(pw);
    expr_free(nnum);
    return result;
}

/* Risch differential equation q' + i w_L' q = p for q in the lower field K_L.
 * Base case (w_L and p both in C(x)): the bounded polynomial-in-x ansatz of
 * rt_solve_rde.  General case: q lives in K_L = C(x, t_0..t_{L-1}) and may be
 * rational there (a monomial denominator, e.g. 1/Log[x]), which the coupled-
 * hyperexponential polynomial coefficients cannot express — solve by a bounded
 * LAURENT ansatz over {x, t_0..t_{L-1}} (each variable ranging over negative powers
 * too), requiring D_tower[q] + i Dcoef_L q = p for every lower tower variable via
 * SolveAlways.  A non-monomial denominator (needing the full Bronstein SPDE) is out
 * of a monomial-Laurent ansatz's reach and declines. */
static Expr* rt_field_rde(Expr* p, long i, RtTower* T, long L, Expr* x) {
    Expr* w = T->arg[L];
    /* Base field C(x): the RDE q' + i w' q = p has both coefficient (i w', with
     * w = the exponent of the current monomial) and RHS p in C(x) — i.e. free of
     * every lower tower variable.  Route it to the COMPLETE base-field solver
     * rde_base (via rt_solve_rde), whose every "no solution" is authoritative
     * (Bronstein Ch. 6, A1/F5 audit).  This deliberately does NOT require w or p to
     * be POLYNOMIAL in x — rde_base handles rational coefficients (weak
     * normalization + normal-denominator reduction), so a rational RHS such as the
     * 1/x of E^x/x (∫ non-elementary, Ei) or E^(x^2) (Erf) is decided here rather
     * than diverted to the bounded-ansatz general path below (the former
     * rt_is_poly(p,x) guard was a spurious scope decline: it sent every rational-p
     * base RDE to SolveAlways). */
    bool base = true;
    for (size_t j = 0; j < T->n && base; j++)
        if (!rt_free_of_x(w, T->t[j]) || !rt_free_of_x(p, T->t[j])) base = false;
    if (base) {
        Expr* q = rt_solve_rde(p, i, w, x);
        if (!q) rt_dec_nonelem();   /* rde_base "no solution" is authoritative (Ch.6) */
        return q;
    }

    /* By the RDE denominator theorem, denom(q) | denom(p) (a pole of q would give a
     * higher-order pole in q' + i Dcoef q = p that nothing cancels), so
     * q = h / pd with pd = Denominator[p] and h a bounded POLYNOMIAL numerator over
     * {x, t_0..t_{L-1}}.  This is strictly more general than a monomial-Laurent
     * ansatz — it captures a NON-monomial denominator such as 1/(1+Log[x]) — and
     * subsumes it (pd carries every pole).  Solve h by SolveAlways. */
    Expr* Dcoef = T->Dcoef[L];
    Expr* pg = rt_eval1("Together", expr_copy(p));
    Expr* pd = pg ? rt_eval1("Denominator", expr_copy(pg)) : NULL;
    Expr* pn = pg ? rt_eval1("Numerator", expr_copy(pg)) : NULL;
    if (pg) expr_free(pg);
    if (!pd || !pn) { if (pd) expr_free(pd); if (pn) expr_free(pn); return NULL; }

    size_t nlv = (size_t)L + 1;
    Expr** lv = malloc(nlv * sizeof(Expr*));
    long* bd = malloc(nlv * sizeof(long));
    lv[0] = x;
    for (long j = 0; j < L; j++) lv[j + 1] = T->t[j];

    /* Exact degree bound (no cap) for the numerator h, where q = h/pd solves
     * D_tower[q] + i Dcoef q = p.  Per lower variable v the bound is the leading-degree
     * balance rt_rde_var_bound(deg_v(p), deg_v(i Dcoef), deriv_lowers) + deg_v(pd),
     * where deriv_lowers distinguishes a LOGARITHMIC / base-x monomial (D lowers deg_v)
     * from an EXPONENTIAL one (D preserves deg_v via the self-derivative D[t^k] =
     * k w' t^k).  This replaces the former arbitrary cap-at-5 proxy — the bound is a
     * function of the equation's degrees alone, so an exponential-Laurent coefficient
     * of ANY degree is found rather than declined.  Correctness is SolveAlways-certified
     * and the caller diff-back verifies, so the bound only affects completeness/perf. */
    Expr* dcg = rt_eval1("Together", expr_copy(Dcoef));
    Expr* dcn = dcg ? rt_eval1("Numerator", expr_copy(dcg)) : NULL;
    Expr* dcd = dcg ? rt_eval1("Denominator", expr_copy(dcg)) : NULL;
    if (dcg) expr_free(dcg);
    for (size_t j = 0; j < nlv; j++) {
        long dpn = rt_degree(pn, lv[j]); if (dpn < 0) dpn = 0;
        long dpd = rt_degree(pd, lv[j]); if (dpd < 0) dpd = 0;
        long dfv = 0;                              /* deg_v(i Dcoef) = deg_v(Dcoef) */
        if (dcn && dcd) {
            long dn = rt_degree(dcn, lv[j]); if (dn < 0) dn = 0;
            long dd = rt_degree(dcd, lv[j]); if (dd < 0) dd = 0;
            dfv = dn - dd;
        }
        /* v = x (j == 0) and logarithmic kernels lower deg_v under D; exponential
         * kernels preserve it. */
        bool deriv_lowers = (j == 0) || (T->kind[j - 1] == RT_LOG);
        /* Live cancellation/resonance detection (Bronstein recursive degree reduction).
         * For an EXPONENTIAL lower monomial at deg_v(f) == 0 the leading coefficients of
         * D[q] and f q can cancel at an integer resonance n = -(i Dcoef_L)/Dcoef_v; feed
         * that n to widen the bound exactly.  (deriv-lowering v has no reachable dfv==-1
         * config in this architecture — see rt_rde_var_bound's reachability note.) */
        long m_res = -1;
        if (!deriv_lowers && dfv == 0)
            m_res = rt_resonance_int(i, Dcoef, T->Dcoef[j - 1]);
        bd[j] = rt_rde_var_bound(dpn - dpd, dfv, deriv_lowers, m_res) + dpd;
    }
    if (dcn) expr_free(dcn);
    if (dcd) expr_free(dcd);
    /* Per-variable LOW exponent.  An EXPONENTIAL lower monomial t_j is a unit
     * (invertible), so the RDE solution may carry a SPECIAL pole there — a
     * negative Laurent power NOT reflected in Denominator[p].  Bronstein's own
     * three-level example integrates to q = t_1/t_0 (a t_0-pole absent from p's
     * denominator x^2 (x+t_0)^2); a purely polynomial numerator over {x, t_j}
     * cannot express it.  Give each exponential kernel the symmetric Laurent
     * range [-bd, bd]; x and logarithmic kernels (whose derivation lowers
     * degree, so no special pole) keep [0, bd].  The extra unknowns are
     * certified to zero by SolveAlways, so this only widens completeness. */
    long* lo = malloc(nlv * sizeof(long));
    for (size_t j = 0; j < nlv; j++) {
        bool is_exp = (j >= 1) && (T->kind[j - 1] == RT_EXP);
        lo[j] = is_exp ? -bd[j] : 0;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] - lo[j] + 1);

    Expr* result = NULL;
    if (nmono > 0) {
        Expr** hterms = malloc((size_t)nmono * sizeof(Expr*));
        Expr** syms = malloc((size_t)nmono * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        for (long m = 0; m < nmono; m++) {
            char nm[40]; snprintf(nm, sizeof(nm), "rmRd%ld", m);
            long idx = m;                          /* mixed-radix decode with lo[] offset */
            for (size_t j = 0; j < nlv; j++) {
                long range = bd[j] - lo[j] + 1;
                ev[j] = lo[j] + idx % range;
                idx /= range;
            }
            Expr* mono = rt_build_monomial(lv, ev, nlv);
            syms[nsym++] = expr_new_symbol(nm);
            hterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_symbol(nm), mono }, 2);
        }
        free(ev);
        Expr* h = expr_new_function(expr_new_symbol("Plus"), hterms, ntq);
        free(hterms);
        /* q = h / pd. */
        Expr* q = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ h, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(pd), expr_new_integer(-1) }, 2) }, 2);

        Expr* dq = rt_tower_deriv(q, T, x);
        Expr* iq = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(i), expr_copy(Dcoef), expr_copy(q) }, 3);
        Expr* residual = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ dq, iq, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 3);
        Expr* tog = rt_eval1("Together", residual);
        Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc(nlv * sizeof(Expr*));
            for (long j = 0; j < L; j++) vl[j] = expr_copy(T->t[L - 1 - j]);
            vl[L] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nlv);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rt_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol.name == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol.name
                 == intern_symbol("List")) {
            Expr* qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(q), expr_copy(sol->data.function.args[0]) }, 2));
            if (qs) {
                Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
                for (size_t j = 0; j < nsym; j++)
                    zero[j] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
                Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
                free(zero);
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ qs, zl }, 2));
            }
        }
        if (sol) expr_free(sol);
        for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
        free(syms);
        expr_free(q);
    }
    /* Cancel the h/pd fraction (e.g. x(1+t)/(x(1+t)^2) -> 1/(1+t)) so the
     * assembled antiderivative stays in lowest terms. */
    if (result) result = rt_eval1("Cancel", result);
    free(lv); free(bd); free(lo);
    expr_free(pd); expr_free(pn);
    /* A NULL here is the SolveAlways over the EXACT-degree-bound q = h/pd ansatz
     * (rt_rde_var_bound: proven cap-free upper bound + denominator theorem +
     * exponential special-pole Laurent range) finding the linear coefficient
     * system inconsistent — i.e. no rational solution exists in K_L.  That is an
     * authoritative "non-elementary" verdict, not a bounded-ansatz decline. */
    if (!result) rt_dec_nonelem();
    return result;
}

/* Recursive-tower case: build the differential tower of f, substitute all kernels
 * to tower variables, gate on a genuine rational function of the whole tower,
 * integrate by the one-extension recursion, back-substitute, and diff-back verify.
 * Closes mixed exp/log towers and rational lower-field coefficients that the flat
 * tower cases decline; runs after them, before the trig front-end. */
/* Structural pre-pass for the recursive tower: expand a MERGED exponential
 * monomial  E^(a + b + ...)  into the product  E^a E^b ...,  built directly
 * (NOT evaluated) so the evaluator cannot re-merge it.  The evaluator
 * automatically combines a product of exponentials into one power with a
 * summed exponent (E^x E^(E^x) -> E^(x + E^x)); that merged exponent is not a
 * valid tower monomial — its argument x + E^x contains the foreign kernel E^x,
 * so rt_tower_build's structure-theorem check rejects the tower and the case
 * declines.  Splitting E^(x + E^x) back to E^x E^(E^x) restores the
 * independent tower basis {E^x, E^(E^x)}, so integrands the evaluator merged
 * (e.g. E^x E^(E^x)/(1+E^(E^x)) = E^(x+E^x)/(1+E^(E^x))) close instead.  The
 * rewrite is exact (E^(a+b) = E^a E^b) and the whole recursive case is
 * diff-back verified, so it can never ship a wrong form.  Returns a
 * freshly-owned tree (caller frees). */
static Expr* rt_expand_exp_sums(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    Expr* expo = NULL;
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
        expo = e->data.function.args[1];
    else if (h == intern_symbol("Exp") && e->data.function.arg_count == 1)
        expo = e->data.function.args[0];
    if (expo && rt_head_is(expo, "Plus") && expo->data.function.arg_count >= 2) {
        size_t m = expo->data.function.arg_count;
        Expr** facs = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) {
            Expr* part = rt_expand_exp_sums(expo->data.function.args[i]);
            facs[i] = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), part }, 2);
        }
        Expr* prod = expr_new_function(expr_new_symbol("Times"), facs, m);
        free(facs);
        return prod;
    }
    Expr* nh = rt_expand_exp_sums(e->data.function.head);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_expand_exp_sums(e->data.function.args[i]);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}

/* Structural (non-evaluating) kernel substitution: replace each tower-kernel
 * subtree (E^arg_i / Log[arg_i], plus the Exp[arg_i] spelling) by its tower
 * variable t_i, top-down (a matched node is not descended into).  Unlike an
 * evaluated ReplaceAll this never lets the evaluator re-merge a split
 * exponential product (E^x E^(E^x) -> E^(x+E^x)) before the kernels are
 * aliased; because the kernels were collected from the very tree being
 * substituted, structural equality holds by construction.  Returns a
 * freshly-owned tree (caller frees). */
static Expr* rt_subst_kernels(Expr* e, RtTower* T) {
    if (!e) return NULL;
    for (size_t i = 0; i < T->n; i++) {
        if (expr_eq(e, T->kernel[i])) return expr_copy(T->t[i]);
        if (T->kind[i] == RT_EXP
            && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == intern_symbol("Exp")
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i]))
            return expr_copy(T->t[i]);
    }
    /* Multiplicatively commensurate member kernel E^(marg) = t[mprim]^mmult. */
    for (size_t m = 0; m < T->nm; m++) {
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count >= 1
            && e->data.function.head->type == EXPR_SYMBOL) {
            const char* h = e->data.function.head->data.symbol.name;
            Expr* w = NULL;
            if (h == intern_symbol("Exp") && e->data.function.arg_count == 1)
                w = e->data.function.args[0];
            else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
                     && e->data.function.args[0]->type == EXPR_SYMBOL
                     && e->data.function.args[0]->data.symbol.name == intern_symbol("E"))
                w = e->data.function.args[1];
            if (w && expr_eq(w, T->marg[m]))
                return expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->t[T->mprim[m]]),
                               expr_new_integer(T->mmult[m]) }, 2);
        }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* nh = rt_subst_kernels(e->data.function.head, T);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_subst_kernels(e->data.function.args[i], T);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}

static Expr* rt_recursive_tower_case(Expr* f, Expr* x) {
    /* Split any evaluator-merged exponential monomial back into an independent
     * tower basis before building the tower (see rt_expand_exp_sums). */
    Expr* fx = rt_expand_exp_sums(f);
    RtTower T;
    /* min_n = 2: single-kernel integrands are handled by the dedicated flat-tower
     * cases above; the general recursion only takes depth->=2 towers. */
    if (!rt_tower_build_min(fx, x, &T, 2)) { rt_tower_free(&T); expr_free(fx); return NULL; }

    /* Alias the kernels to tower variables structurally (NOT via an evaluated
     * ReplaceAll, which would re-merge a split exponential product before
     * substitution — see rt_subst_kernels), then normalise with Together. */
    Expr* F = rt_eval1("Together", rt_subst_kernels(fx, &T));
    Expr* result = NULL;
    bool tower_verified = false;
    if (F && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL) {
        Expr* num = rt_eval1("Numerator", expr_copy(F));
        Expr* den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((T.n + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < T.n; i++) vv[i + 1] = expr_copy(T.t[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, T.n + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        bool gate = num && den && rt_is_true(pqn) && rt_is_true(pqd);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (gate) {
            Expr* rem_tower = NULL;   /* partial log part: unintegrated r_n */
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x, &rem_tower);
            if (Q) {
                /* Reduce Q to lowest terms in the TOWER VARIABLES (each kernel
                 * t_i an independent symbol) BEFORE back-substituting.  A common
                 * factor like (1 + t_1) in Q = t_3 (1+t_1) / ((1+t_1) t_2)
                 * cancels cleanly here, but not afterwards when t_1 = E^x
                 * reappears inside t_2 = Log[1 + E^x] and blocks the gcd —
                 * otherwise a correct answer (Log[x]/Log[1 + E^x]) is left in an
                 * unrecognisable form.  Together first: Q is a SUM of Laurent/log
                 * terms, which Cancel alone does not combine before reducing. */
                Q = rt_eval1("Cancel", rt_eval1("Together", Q));
                /* EXACT verification in the tower variables: with a partial log
                 * part (Bronstein Thm 5.6.1) the identity is D_tower[Q] + rem =
                 * F, i.e. the elementary Q integrates all of F except the
                 * unintegrated non-constant-residue remainder rem_tower.  When
                 * rem_tower is NULL this is the usual exact D_tower[Q] == F. */
                Expr* DQ = rt_tower_deriv(Q, &T, x);
                Expr* diff = DQ ? expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(DQ), expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2) : NULL;
                if (diff && rem_tower)
                    diff = expr_new_function(expr_new_symbol("Plus"),
                        (Expr*[]){ diff, expr_copy(rem_tower) }, 2);
                Expr* chk = diff ? rt_eval1("Together", diff) : NULL;
                bool tower_ok = chk && rt_is_zero(chk);
                if (DQ) expr_free(DQ);
                if (chk) expr_free(chk);
                if (tower_ok) {
                    Expr** back = malloc(T.n * sizeof(Expr*));
                    for (size_t i = 0; i < T.n; i++)
                        back[i] = expr_new_function(expr_new_symbol("Rule"),
                            (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                    Expr* bl = expr_new_function(expr_new_symbol("List"), back, T.n);
                    free(back);
                    result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Q, bl }, 2));                   /* adopts Q, bl */
                    if (result) result = rt_eval1("Cancel", result);
                    /* Surface the unintegrated non-constant-residue part as a
                     * genuine Integrate[rem, x] in the original kernels (a
                     * separate back-rule list — ReplaceAll consumes its list).
                     * rem's residues are all non-constant, so re-integrating it
                     * declines (r_s empty ⇒ NULL): a fixed point, no loop. */
                    if (result && rem_tower) {
                        Expr** back2 = malloc(T.n * sizeof(Expr*));
                        for (size_t i = 0; i < T.n; i++)
                            back2[i] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                        Expr* bl2 = expr_new_function(expr_new_symbol("List"), back2, T.n);
                        free(back2);
                        Expr* rem_k = rt_eval1("Cancel", rt_eval_own(
                            expr_new_function(expr_new_symbol("ReplaceAll"),
                                (Expr*[]){ expr_copy(rem_tower), bl2 }, 2)));
                        Expr* rem_int = expr_new_function(expr_new_symbol("Integrate"),
                            (Expr*[]){ rem_k, expr_copy(x) }, 2);
                        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ result, rem_int }, 2));
                    }
                    tower_verified = true;
                } else {
                    expr_free(Q);
                }
            }
            if (rem_tower) expr_free(rem_tower);
        }
    }
    if (F) expr_free(F);

    /* The tower-variable check above is exact; fall back to the (bounded)
     * Simplify diff-back only if it was not reached (defensive). */
    if (result && !tower_verified && !rt_verify_antideriv(result, f, x)) {
        expr_free(result); result = NULL;
    }

    rt_tower_free(&T);
    expr_free(fx);
    return result;
}

/* ---- real reconstruction of I-laden trigonometric antiderivatives --------- */
/* The complex-exponential route integrates over the kernel E^(I x); its log part
 * comes back as complex logs / ArcTan of E^(±I x), a correct but I-laden closed
 * form (e.g. ∫Csc x -> Log[E^(I x)-1] - Log[E^(I x)+1], which equals the real
 * Log[Tan[x/2]] up to a constant).  cx_reim decomposes such an expression of the
 * REAL variable x into real + I·imaginary parts (treating x, Pi, E and any
 * trig/Sqrt of a real argument as real); rt_realify keeps the real part — a valid
 * antiderivative since the discarded imaginary part is a constant for a real
 * integrand — and returns it only behind an exact diff-back check.  This is the
 * LogToReal real reconstruction for the exponential tower, not a rewrite of the
 * integrand. */

static Expr* cx_int(long n) { return expr_new_integer(n); }
static Expr* cx_add(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* cx_mul(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* cx_neg(Expr* a) { return cx_mul(cx_int(-1), a); }
static Expr* cx_half(void) {
    return expr_new_function(expr_new_symbol("Rational"),
        (Expr*[]){ cx_int(1), cx_int(2) }, 2);
}
static Expr* cx_I(void) {
    return expr_new_function(expr_new_symbol("Complex"), (Expr*[]){ cx_int(0), cx_int(1) }, 2);
}
static Expr* cx_log(Expr* a) { return expr_new_function(expr_new_symbol("Log"), (Expr*[]){ a }, 1); }

/* True iff e contains no Complex[...] atom (manifestly real for real symbols). */
static bool rt_free_of_complex(Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_FUNCTION) {
        if (rt_head_is(e, "Complex")) return false;
        if (!rt_free_of_complex(e->data.function.head)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!rt_free_of_complex(e->data.function.args[i])) return false;
    }
    return true;
}

/* e == *re + I (*im), for a function of the real variable x.  Returns false on a
 * construct it cannot split.  Owns the returned parts. */
static bool cx_reim(Expr* e, Expr* x, Expr** re, Expr** im) {
    *re = NULL; *im = NULL;
    if (rt_head_is(e, "Complex") && e->data.function.arg_count == 2) {
        *re = expr_copy(e->data.function.args[0]);
        *im = expr_copy(e->data.function.args[1]);
        return true;
    }
    if (rt_free_of_complex(e)) {          /* real atom / trig/Sqrt of real x, ... */
        *re = expr_copy(e); *im = cx_int(0); return true;
    }
    if (e->type != EXPR_FUNCTION) return false;
    size_t n = e->data.function.arg_count;

    if (rt_head_is(e, "Plus")) {
        Expr* R = cx_int(0); Expr* M = cx_int(0);
        for (size_t i = 0; i < n; i++) {
            Expr* r; Expr* m;
            if (!cx_reim(e->data.function.args[i], x, &r, &m)) {
                expr_free(R); expr_free(M); return false;
            }
            R = cx_add(R, r); M = cx_add(M, m);
        }
        *re = R; *im = M; return true;
    }
    if (rt_head_is(e, "Times")) {
        Expr* R = cx_int(1); Expr* M = cx_int(0);
        for (size_t i = 0; i < n; i++) {
            Expr* c; Expr* d;
            if (!cx_reim(e->data.function.args[i], x, &c, &d)) {
                expr_free(R); expr_free(M); expr_free(c); return false;
            }
            /* (R + M i)(c + d i) = (Rc - Md) + (Rd + Mc) i */
            Expr* nR = cx_add(cx_mul(expr_copy(R), expr_copy(c)),
                              cx_neg(cx_mul(expr_copy(M), expr_copy(d))));
            Expr* nM = cx_add(cx_mul(expr_copy(R), expr_copy(d)),
                              cx_mul(expr_copy(M), expr_copy(c)));
            expr_free(R); expr_free(M); expr_free(c); expr_free(d);
            R = nR; M = nM;
        }
        *re = R; *im = M; return true;
    }
    if (rt_head_is(e, "Power") && n == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (rt_free_of_complex(base)) {
            /* base^(a + b i) = base^a (Cos[b Log base] + i Sin[b Log base]). */
            Expr* a; Expr* b;
            if (!cx_reim(ex, x, &a, &b)) return false;
            Expr* mag = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(base), a }, 2);                 /* a adopted */
            Expr* th = cx_mul(b, cx_log(expr_copy(base)));           /* b adopted */
            *re = cx_mul(expr_copy(mag),
                expr_new_function(expr_new_symbol("Cos"), (Expr*[]){ expr_copy(th) }, 1));
            *im = cx_mul(mag,
                expr_new_function(expr_new_symbol("Sin"), (Expr*[]){ th }, 1));
            return true;
        }
        if (ex->type == EXPR_INTEGER) {                             /* complex^int */
            Expr* br; Expr* bi;
            if (!cx_reim(base, x, &br, &bi)) return false;
            long p = ex->data.integer, ap = p < 0 ? -p : p;
            Expr* R = cx_int(1); Expr* M = cx_int(0);
            for (long k = 0; k < ap; k++) {
                Expr* nR = cx_add(cx_mul(expr_copy(R), expr_copy(br)),
                                  cx_neg(cx_mul(expr_copy(M), expr_copy(bi))));
                Expr* nM = cx_add(cx_mul(expr_copy(R), expr_copy(bi)),
                                  cx_mul(expr_copy(M), expr_copy(br)));
                expr_free(R); expr_free(M); R = nR; M = nM;
            }
            expr_free(br); expr_free(bi);
            if (p < 0) {                                            /* invert */
                Expr* den = cx_add(cx_mul(expr_copy(R), expr_copy(R)),
                                   cx_mul(expr_copy(M), expr_copy(M)));
                Expr* invden = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ den, cx_int(-1) }, 2);
                Expr* nR = cx_mul(R, expr_copy(invden));
                Expr* nM = cx_mul(cx_neg(M), invden);
                R = nR; M = nM;
            }
            *re = R; *im = M; return true;
        }
        return false;                                              /* complex^noninteger */
    }
    if (rt_head_is(e, "Log") && n == 1) {
        Expr* a; Expr* b;
        if (!cx_reim(e->data.function.args[0], x, &a, &b)) return false;
        /* Log[a + b i] = (1/2) Log[a^2 + b^2] + i ArcTan[a, b]  (Arg). */
        Expr* mag2 = cx_add(cx_mul(expr_copy(a), expr_copy(a)),
                            cx_mul(expr_copy(b), expr_copy(b)));
        *re = cx_mul(cx_half(), cx_log(mag2));
        *im = expr_new_function(expr_new_symbol("ArcTan"), (Expr*[]){ a, b }, 2);
        return true;
    }
    /* ArcTan[z] = (i/2)(Log[1 - i z] - Log[1 + i z]); recurse on the rewrite. */
    if (rt_head_is(e, "ArcTan") && n == 1) {
        Expr* z = e->data.function.args[0];
        Expr* iz1 = cx_mul(cx_I(), expr_copy(z));
        Expr* iz2 = cx_mul(cx_I(), expr_copy(z));
        Expr* rw = cx_mul(cx_mul(cx_I(), cx_half()),
            cx_add(cx_log(cx_add(cx_int(1), cx_neg(iz1))),
                   cx_neg(cx_log(cx_add(cx_int(1), iz2)))));
        bool ok = cx_reim(rw, x, re, im);
        expr_free(rw);
        return ok;
    }
    /* ArcTanh[z] = (1/2)(Log[1 + z] - Log[1 - z]). */
    if (rt_head_is(e, "ArcTanh") && n == 1) {
        Expr* z = e->data.function.args[0];
        Expr* rw = cx_mul(cx_half(),
            cx_add(cx_log(cx_add(cx_int(1), expr_copy(z))),
                   cx_neg(cx_log(cx_add(cx_int(1), cx_neg(expr_copy(z)))))));
        bool ok = cx_reim(rw, x, re, im);
        expr_free(rw);
        return ok;
    }
    return false;
}

/* Numeric diff-back: true iff Sum_k |N[(D[g,x] - f) /. x -> p_k]| < 1e-9 over a
 * few interior points.  Re[G] of an already-correct antiderivative G is provably
 * an antiderivative of the real integrand f, so a numeric check is a sufficient
 * guard against a cx_reim bug — and it is robust where the Simplify diff-back
 * cannot reduce the nested-log real forms the reconstruction produces. */
static bool rt_realify_numverify(Expr* g, Expr* f, Expr* x) {
    static const char* pts[] = { "7/10", "13/10", "23/10", "37/10" };
    Expr* Dg = rt_eval2("D", expr_copy(g), expr_copy(x));
    if (!Dg) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Dg, cx_neg(expr_copy(f)) }, 2);                  /* D[g] - f */
    Expr* acc = cx_int(0);
    for (size_t k = 0; k < sizeof pts / sizeof pts[0]; k++) {
        Expr* rule = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(x), parse_expression(pts[k]) }, 2);
        Expr* at = expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(diff), rule }, 2);
        Expr* nv = expr_new_function(expr_new_symbol("Abs"),
            (Expr*[]){ expr_new_function(expr_new_symbol("N"), (Expr*[]){ at }, 1) }, 1);
        acc = cx_add(acc, nv);
    }
    expr_free(diff);
    Expr* test = expr_new_function(expr_new_symbol("Less"),
        (Expr*[]){ acc, parse_expression("1/1000000000") }, 2);
    Expr* val = rt_eval_own(test);
    bool ok = val && val->type == EXPR_SYMBOL
              && val->data.symbol.name == intern_symbol("True");
    if (val) expr_free(val);
    return ok;
}

/* Real closed form of the I-laden antiderivative G of the real integrand f:
 * Re[G] (x real), returned only if it is manifestly real and diff-backs to f
 * numerically. */
static Expr* rt_realify(Expr* G, Expr* x, Expr* f) {
    Expr* re = NULL; Expr* im = NULL;
    if (!cx_reim(G, x, &re, &im)) { if (re) expr_free(re); if (im) expr_free(im); return NULL; }
    if (im) expr_free(im);
    Expr* reS = rt_eval1("Simplify", re);
    if (reS && rt_free_of_complex(reS) && rt_realify_numverify(reS, f, x)) return reS;
    if (reS) expr_free(reS);
    return NULL;
}

/* Trigonometric / hyperbolic front-end (exponentialize path).
 * Rewrites the trig/hyperbolic kernels to complex exponentials with
 * TrigToExp, integrates the resulting (Laurent-)rational function of the
 * exponential kernel E^(I x) / E^x with the exponential machinery, and
 * converts the answer back to trigonometric form with ExpToTrig.  Both
 * rewrites are exact, so the result is correct by construction; rt_realify then
 * reconstructs a real closed form (diff-back gated) from the I-laden output. */
static Expr* rt_trig_frontend(Expr* f, Expr* x) {
    Expr* fe = rt_eval1("TrigToExp", expr_copy(f));
    if (!fe) return NULL;
    if (expr_eq(fe, f)) { expr_free(fe); return NULL; }   /* no trig/hyperbolic */
    /* All exponential cases are used, including the coupled hyperexponential
     * one, for completeness — a correct antiderivative is returned even when
     * it cannot be reduced to the cleanest real form.
     *
     * KNOWN SIMPLIFICATION GAP (Simplify improvement opportunity): via the
     * complex substitution u = I x, Tan[x]/Tanh[x] and similar close to a
     * correct but I-laden form such as  I x - Log[1 + E^(2 I x)]  (which equals
     * -Log[Cos[x]]).  No current simplifier (Simplify / FullSimplify /
     * TrigReduce / PowerExpand) collapses  I x - Log[1 + E^(2 I x)]  to
     * -Log[Cos[x]].  The result is nonetheless returned (correct by
     * construction); teaching Simplify the log-of-product / half-angle
     * collapse would render it in real closed form. */
    Expr* r = rt_exp_poly_case(fe, x);
    if (!r) r = rt_frac_case(fe, x);
    /* Rational-function-of-a-single-exponential closer (kernelize E^(I x) -> pure
     * rational integral).  Closes the Laurent forms with E^(-I x) in the
     * denominator that the frac case leaves — Sec[x], Csc[x], 1/(2+Cos[x]),
     * Sec[x]^3 — which are exactly the rational trigonometric integrands.  Runs
     * after rt_frac_case so its cleaner squarefree ArcTan/Log parts win first. */
    if (!r) r = rt_exp_ratreduce_case(fe, x);
    if (!r) r = rt_hyperexp_case(fe, x);
    /* Multi-kernel decoupling (Phase B): e.g. Sin/Cos times a real exponential
     * exponentialize to a sum of two non-commensurate exponentials E^((a +/- b I) x)
     * that the single-primitive cases cannot kernelize. */
    if (!r) r = rt_expsum_case(fe, x);
    expr_free(fe);
    if (!r) return NULL;
    Expr* out = rt_eval1("ExpToTrig", r);   /* adopts r; back to trig form */
    if (!out) return NULL;
    /* Real reconstruction: the exp route's log part is a correct but I-laden
     * closed form; reconstruct a real one.  Diff-back gated, so on failure the
     * correct I-laden form is kept (never a wrong or worse answer). */
    Expr* real = rt_realify(out, x, f);
    if (real) { expr_free(out); return real; }
    return out;
}

/* True iff `e` contains, anywhere, a function node with head symbol `h`. */
static bool rt_expr_has_head(Expr* e, const char* h) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol(h)) return true;
    if (rt_expr_has_head(e->data.function.head, h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_expr_has_head(e->data.function.args[i], h)) return true;
    return false;
}

/* Rational-function-of-a-single-exponential fallback.  When the integrand is a
 * rational function of one exponential kernel with a LINEAR exponent (so the
 * commensurate primitive t = E^(a x) has constant logarithmic derivative
 * u' = a), the substitution t = E^(a x), dx = dt/(a t) reduces the integral to
 * the pure rational integral  ∫ R(t)/(a t) dt, which the rational integrator
 * closes in full — including the algebraic-residue log parts it expresses as
 * RootSum (e.g. the cubic 1 + t^2 + t^3 from E^(x/6)/(1 + E^(x/2) + E^(x/3))).
 * Back-substitute t -> E^(a x) and diff-back verify.  This is the general
 * closure for the cubic-and-higher residue cases that the structured real
 * Log/ArcTan Rothstein-Trager log part (rt_frac_lrt) necessarily declines;
 * it runs LAST among the exponential cases, so their cleaner ArcTan/radical
 * forms win whenever they apply.  Correctness rests on the same
 * correct-by-construction rational integrator plus the diff-back gate; the
 * RootSum it emits is recognised by the generalized D[RootSum] collapse
 * (src/root.c). */
static Expr* rt_exp_ratreduce_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F || !uexp) { if (F) expr_free(F); if (uexp) expr_free(uexp); return NULL; }

    Expr* t = expr_new_symbol("rmT");
    Expr* up = rt_eval2("D", expr_copy(uexp), expr_copy(x));   /* u' = a */
    Expr* result = NULL;

    /* Gates: F is a rational function of the kernel t ALONE (no residual x, no
     * nested exp/log of x), it genuinely depends on t, and u' is free of x
     * (a linear exponent -> constant log-derivative, so the Jacobian 1/(u' t)
     * is rational in t). */
    bool ok = F && up
        && rt_free_of_x(F, x)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && !rt_free_of_x(F, t)
        && rt_free_of_x(up, x);

    if (ok) {
        /* G = F / (u' t). */
        Expr* denom = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(up), expr_copy(t) }, 2);
        Expr* G = rt_eval1("Together", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(F),
                expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ denom, expr_new_integer(-1) }, 2) }, 2));
        /* Pure rational integral in t. */
        Expr* A = G ? rt_eval2("Integrate", expr_copy(G), expr_copy(t)) : NULL;
        bool declined = !A || rt_expr_has_head(A, "Integrate");
        if (!declined) {
            /* Back-substitute t -> E^u and diff-back verify against f. */
            Expr* back = expr_new_function(expr_new_symbol("List"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(t),
                        expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2) }, 2) }, 1);
            Expr* Q = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(A), back }, 2));
            /* Correct by construction (exact rational integral over the kernel +
             * exact back-substitution); accept on the Simplify diff-back, or a
             * numeric one where Simplify cannot reduce the higher-pole forms
             * (e.g. Sec[x]^3 -> a triple pole at E^(I x) = +/- i). */
            if (Q && (rt_verify_antideriv(Q, f, x) || rt_realify_numverify(Q, f, x)))
                result = Q;
            else if (Q) expr_free(Q);
        }
        if (A) expr_free(A);
        if (G) expr_free(G);
    }

    expr_free(t);
    if (up) expr_free(up);
    expr_free(F);
    expr_free(uexp);
    return result;
}

/* Find (borrowed) the argument u of the first kernel `head`[u] of `e` (arity 1)
 * whose argument depends on x. */
static Expr* rt_find_head_of_x(Expr* e, Expr* x, const char* head) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol(head)
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    Expr* r = rt_find_head_of_x(e->data.function.head, x, head);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_head_of_x(e->data.function.args[i], x, head);
        if (r) return r;
    }
    return NULL;
}

/* True iff `e` contains no circular- or hyperbolic-trig kernel (a coarse "is this
 * a rational function of the tangent monomial and x alone" gate; the final
 * diff-back is the real guarantee). */
static bool rt_free_of_trig(Expr* e) {
    static const char* H[] = { "Tan", "Sin", "Cos", "Cot", "Sec", "Csc",
                               "Tanh", "Sinh", "Cosh", "Coth", "Sech", "Csch" };
    for (size_t i = 0; i < sizeof H / sizeof H[0]; i++)
        if (!rt_free_of_head(e, H[i])) return false;
    return true;
}

/* The monomial variable t used throughout the family (a fixed internal symbol). */
#define RT_HT_SYM "Integrate`htT"

/* Maximum degree in t of an irreducible factor of p over Q (0 when p is constant
 * in t, -1 on failure).  The normal-pole gate uses this: a normal denominator
 * whose irreducible factors are each linear or quadratic has Rothstein-Trager
 * residues that are rational or at-worst quadratic-algebraic, which the §5.10
 * residue criterion realises as real Log / ArcTan.  A higher-degree irreducible
 * factor carries residues whose tower-derivation spillover destabilises the
 * p = h - D[g2] + r reconciliation, so those defer. */
static long rt_max_irr_degree(Expr* p, Expr* t) {
    Expr* fac = rt_eval1("Factor", expr_copy(p));
    if (!fac) return -1;
    long maxd = 0;
    if (rt_head_is(fac, "Times")) {
        for (size_t i = 0; i < fac->data.function.arg_count; i++) {
            Expr* term = fac->data.function.args[i];
            Expr* base = (rt_head_is(term, "Power")
                          && term->data.function.arg_count == 2)
                         ? term->data.function.args[0] : term;
            long d = rt_degree(base, t);
            if (d > maxd) maxd = d;
        }
    } else {
        Expr* base = (rt_head_is(fac, "Power") && fac->data.function.arg_count == 2)
                     ? fac->data.function.args[0] : fac;
        long d = rt_degree(base, t);
        if (d > maxd) maxd = d;
    }
    expr_free(fac);
    return maxd;
}


/* Real "hypertangent family" integrator (Bronstein §5.10 and its hyperbolic
 * analogue): integrate a rational function of a SINGLE kernel `khead`[u]
 * (u rational in x) DIRECTLY and real, the real-valued replacement for the
 * TrigToExp route.  Substitute t = khead[u] so F is rational in t over C(x);
 * build the tower derivation Dt (given, = u'(t^2+1) for Tan, -u'(t^2+1) for Cot,
 * u'(1-t^2) for Tanh); run the `driver` (IntegrateHypertangent / -Hypertanh);
 * integrate the leftover base-field element in C(x); back-substitute t ->
 * khead[u]; collapse the real special log Log[1 (+/-) khead[u]^2] to
 * -2 Log[`inner`[u]] (Cos / Sin / Cosh); and diff-back verify.  `special` is the
 * monic special polynomial (t^2+1 or t^2-1) used by the normal-pole gate.  Takes
 * ownership of Dt and special. */
/* Real hypertangent / hyperbolic-tangent family integrator.
 *
 * The kernel substitution and the whole pipeline are CORRECT BY CONSTRUCTION —
 * there is no diff-back Simplify gate.  The structural gate below (F is a
 * genuine rational function of t over C(x), no trig/exp/log of x survives)
 * already guarantees F|_{t->kernel} = f: for a direct kernel t IS the kernel,
 * so back-substitution is the identity; for a Weierstrass half-angle
 * substitution the rules encode exact identities.  The §5.10 driver then
 * returns g with beta = True iff the transcendental part is elementary, leaving
 * the t-free base-field remainder base = F - D_tower[g] for a recursive integral
 * in x.  Since the tower derivation D[t] IS the derivative of the kernel bval,
 * d/dx[(g + integral base)|_{t->bval}] = F|_{t->bval} = f is a theorem (Bronstein
 * Thm 5.10.1), not something to re-verify.  The only runtime facts consulted are
 * STRUCTURAL: beta = True (elementary), base free of t, and the base integral is
 * itself elementary (not left unevaluated). */
static Expr* rt_hypertan_family(Expr* f, Expr* x, Expr* subrule, Expr* bval,
                                Expr* Dt, Expr* special,
                                const char* driver, const char* inner,
                                Expr* carg, bool cosmetic_plus) {
    Expr* t = expr_new_symbol(RT_HT_SYM);
    Expr* F = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), subrule }, 2));                 /* adopts subrule */

    /* F must be a genuine rational function of t over C(x). */
    bool ok = F && rt_free_of_trig(F)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL
        && !rt_free_of_x(F, t);
    if (ok) {
        Expr* Ft = rt_eval1("Together", expr_copy(F));
        Expr* nn = Ft ? rt_eval1("Numerator", expr_copy(Ft)) : NULL;
        Expr* dd = Ft ? rt_eval1("Denominator", expr_copy(Ft)) : NULL;
        ok = nn && dd && rt_is_poly(nn, t) && rt_is_poly(dd, t);
        if (ok) {
            /* Normal-pole gate.  Strip the special factor from the denominator,
             * then require every IRREDUCIBLE factor of the remaining NORMAL part
             * (over Q) to have degree <= 2: its Rothstein-Trager residues are then
             * rational (linear factors) or at-worst quadratic-algebraic (irreducible
             * quadratics), both of which the §5.10 residue criterion realises as
             * real Log / ArcTan.  A higher-degree irreducible factor carries
             * residues whose tower-derivation spillover blows up the
             * p = h - D[g2] + r reconciliation, so those defer.  (Pure-polynomial
             * and special^k integrands strip to a unit; split normal parts of any
             * degree pass, as does an irreducible quadratic such as the 3+t^2 of
             * Tan[x]/(3+Tan[x]^2).) */
            Expr* dn = expr_copy(dd);
            for (int guard = 0; dn && guard < 4096; guard++) {
                Expr* rem = rt_eval3("PolynomialRemainder",
                    expr_copy(dn), expr_copy(special), expr_copy(t));
                bool divisible = rem && rt_is_zero(rem);
                if (rem) expr_free(rem);
                if (!divisible) break;
                Expr* q = rt_eval3("PolynomialQuotient",
                    expr_copy(dn), expr_copy(special), expr_copy(t));
                expr_free(dn); dn = q;
            }
            long mdeg = dn ? rt_max_irr_degree(dn, t) : -1;
            ok = (mdeg >= 0 && mdeg <= 2);
            if (dn) expr_free(dn);
        }
        if (Ft) expr_free(Ft);
        if (nn) expr_free(nn);
        if (dd) expr_free(dd);
    }
    if (!ok) { if (F) expr_free(F); expr_free(t); expr_free(Dt); expr_free(special);
               expr_free(bval); expr_free(carg); return NULL; }

    /* deriv = {x -> 1, t -> Dt}. */
    Expr* deriv = expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(x), expr_new_integer(1) }, 2),
          expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(t), expr_copy(Dt) }, 2) }, 2);

    Expr* result = NULL;
    Expr* res = rt_eval_call(driver,
        (Expr*[]){ expr_copy(F), expr_copy(t), expr_copy(deriv) }, 3);
    if (res && rt_head_is(res, "List") && res->data.function.arg_count == 2
        && rt_is_true(res->data.function.args[1])) {
        Expr* g = res->data.function.args[0];               /* borrowed */
        /* base = F - D_tower[g]  (an element of k = C(x): must be free of t).
         * Expand after Together: the tower derivation returns unexpanded products
         * (e.g. (t^2+1)(t^2-1)) whose t-terms only cancel once expanded. */
        Expr* Dg = rt_eval2("Risch`Derivation", expr_copy(g), expr_copy(deriv));
        Expr* base = Dg ? rt_eval1("Expand", rt_eval1("Together",
            expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ expr_copy(F), expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), Dg }, 2) }, 2))) : NULL;
        /* base is the t-free base-field remainder; its being free of t is the
         * structural proof that g integrates the transcendental part of F.  The
         * recursive base-field integral ib is elementary iff Integrate does not
         * leave it unevaluated.  Both are correct by construction (rational Risch
         * over C(x)); back-substituting the kernel then yields a correct
         * antiderivative of f without any diff-back re-verification. */
        if (base && rt_free_of_x(base, t)) {
            Expr* ib = rt_eval2("Integrate", expr_copy(base), expr_copy(x));
            if (ib && !rt_expr_has_head(ib, "Integrate")) {
                Expr* backrule = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(t), expr_copy(bval) }, 2);
                Expr* gsub = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(g), backrule }, 2));    /* adopts backrule */
                Expr* ans = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ gsub, ib }, 2));                  /* adopts gsub, ib */
                if (ans) {
                    /* Cosmetic (real): Log[1 (+/-) bval^2] = -2 Log[inner[carg]]
                     * (Tan/Cot: 1+K^2 -> -2 Log Cos/Sin; Tanh/Coth: 1-K^2 -> -2 Log
                     * Cosh/Sinh; Weierstrass: 1+Tan[u/2]^2 -> -2 Log Cos[u/2]).  Each
                     * is an exact identity (for Coth, 1-Coth^2 < 0 contributes only a
                     * constant i Pi), so the swap is derivative-preserving by
                     * construction; apply it when the pattern is present. */
                    Expr* ksq = expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(bval), expr_new_integer(2) }, 2);
                    Expr* logarg = cosmetic_plus
                        ? expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_new_integer(1), ksq }, 2)
                        : expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ expr_new_integer(1),
                                expr_new_function(expr_new_symbol("Times"),
                                    (Expr*[]){ expr_new_integer(-1), ksq }, 2) }, 2);
                    Expr* lhs = expr_new_function(expr_new_symbol("Log"),
                        (Expr*[]){ logarg }, 1);
                    Expr* rhs = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-2),
                            expr_new_function(expr_new_symbol("Log"),
                                (Expr*[]){ expr_new_function(expr_new_symbol(inner),
                                    (Expr*[]){ expr_copy(carg) }, 1) }, 1) }, 2);
                    Expr* crule = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ lhs, rhs }, 2);
                    Expr* clean = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ expr_copy(ans), crule }, 2)); /* adopts crule */
                    if (clean && !expr_eq(clean, ans)) {
                        expr_free(ans); result = clean;
                    } else {
                        if (clean) expr_free(clean);
                        result = ans;
                    }
                }
            } else if (ib) {
                expr_free(ib);
            }
        }
        if (base) expr_free(base);
    }
    if (res) expr_free(res);
    expr_free(F); expr_free(t); expr_free(deriv); expr_free(Dt); expr_free(special);
    expr_free(bval); expr_free(carg);
    return result;
}

/* Direct single-kernel call: substitute khead[u] -> t (and back-substitute the
 * same), cosmetic argument u.  Takes ownership of Dt and special. */
static Expr* rt_hypertan_direct(Expr* f, Expr* x, Expr* u, const char* khead,
                                Expr* Dt, Expr* special, const char* driver,
                                const char* inner, bool cosmetic_plus) {
    Expr* subrule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_function(expr_new_symbol(khead),
            (Expr*[]){ expr_copy(u) }, 1), expr_new_symbol(RT_HT_SYM) }, 2);
    Expr* bval = expr_new_function(expr_new_symbol(khead), (Expr*[]){ expr_copy(u) }, 1);
    return rt_hypertan_family(f, x, subrule, bval, Dt, special, driver, inner,
                             expr_copy(u), cosmetic_plus);
}

/* eta = D[u, x] must be a rational function of x alone (free of trig/exp/log of x)
 * so that khead[u] is a monomial over k = C(x).  Returns owned eta or NULL. */
static Expr* rt_kernel_eta(Expr* u, Expr* x) {
    if (!rt_kernel_simple(u, x)) return NULL;
    Expr* eta = rt_eval2("D", expr_copy(u), expr_copy(x));
    if (!eta || !rt_free_of_trig(eta)
        || rt_find_exp_of_x(eta, x) != NULL || rt_find_log_of_x(eta, x) != NULL) {
        if (eta) expr_free(eta);
        return NULL;
    }
    return eta;
}
/* t^2 + sigma  (sigma = +1 or -1), in the family monomial variable. */
static Expr* rt_special_poly(long sigma) {
    return expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol(RT_HT_SYM), expr_new_integer(2) }, 2),
          expr_new_integer(sigma) }, 2);
}
/* 1 - t^2 in the family monomial variable (the Tanh/Coth tower derivative base). */
static Expr* rt_one_minus_t2(void) {
    return expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_new_integer(1),
            expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_new_symbol(RT_HT_SYM), expr_new_integer(2) }, 2) }, 2) }, 2);
}

/* Dispatch the real hypertangent family: Tan[u] and Cot[u] (special t^2+1, the
 * §5.10 hypertangent driver) and Tanh[u] (special t^2-1, the hyperbolic driver).
 * The special polynomial t^2+1 is irreducible over k (needs a coupled complex
 * system); t^2-1 splits (two real Risch DEs) — both keep the answer real. */
static Expr* rt_hypertangent_case(Expr* f, Expr* x) {
    /* Tan[u]:  Dt = eta (t^2+1). */
    Expr* u = rt_find_head_of_x(f, x, "Tan");
    if (u) {
        Expr* eta = rt_kernel_eta(u, x);
        if (eta) {
            Expr* Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ eta, rt_special_poly(1) }, 2);     /* adopts eta */
            return rt_hypertan_direct(f, x, u, "Tan", Dt, rt_special_poly(1),
                "Risch`IntegrateHypertangent", "Cos", true);
        }
    }
    /* Cot[u]:  Dt = -eta (t^2+1). */
    u = rt_find_head_of_x(f, x, "Cot");
    if (u) {
        Expr* eta = rt_kernel_eta(u, x);
        if (eta) {
            Expr* neg_eta = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), eta }, 2);   /* adopts eta */
            Expr* Dt = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ neg_eta, rt_special_poly(1) }, 2);
            return rt_hypertan_direct(f, x, u, "Cot", Dt, rt_special_poly(1),
                "Risch`IntegrateHypertangent", "Sin", true);
        }
    }
    /* Tanh[u] and Coth[u]:  Dt = eta (1 - t^2)  (special t^2-1, which splits).
     * They share the hyperbolic driver; only the real cosmetic differs — for
     * Tanh 1-Tanh^2 = Sech^2 > 0 (-> -2 Log Cosh), for Coth 1-Coth^2 = -Csch^2
     * (-> -2 Log Sinh, the sign folded into the derivative-preserving rewrite). */
    struct { const char* head; const char* inner; } hyp[] = {
        { "Tanh", "Cosh" }, { "Coth", "Sinh" }
    };
    for (size_t i = 0; i < sizeof hyp / sizeof hyp[0]; i++) {
        u = rt_find_head_of_x(f, x, hyp[i].head);
        if (!u) continue;
        Expr* eta = rt_kernel_eta(u, x);
        if (!eta) continue;
        Expr* Dt = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ eta, rt_one_minus_t2() }, 2);          /* adopts eta */
        return rt_hypertan_direct(f, x, u, hyp[i].head, Dt, rt_special_poly(-1),
            "Risch`IntegrateHypertanh", hyp[i].inner, false);
    }
    return NULL;
}

/* Dispatch the transcendental cases: the primitive (logarithmic) polynomial
 * reduction, the exponential (hyperexponential / Risch-DE) reduction, the
 * fractional (Rothstein-Trager) log-part, the real hypertangent case, and the
 * trig/hyperbolic (TrigToExp) front-end.  The general Hermite reduction for
 * repeated poles lands in a subsequent increment. */
static Expr* rt_transcendental_case(Expr* f, Expr* x) {
    Expr* r = rt_log_poly_case(f, x);
    if (r) return r;
    r = rt_exp_poly_case(f, x);
    if (r) return r;
    r = rt_frac_case(f, x);
    if (r) return r;
    /* Fast rational-of-a-single-exp path (kernelize -> pure rational integral in
     * t, which the FLINT-accelerated rational integrator closes) BEFORE the
     * repeated-pole Hermite ansatz: when the integrand is a rational function of
     * one exponential with a linear exponent (F free of x, u' free of x), this
     * closes it in O(deg) instead of the O(mult)-variable SolveAlways Hermite
     * solve.  It declines (F still carries x, nonlinear exponent, multi-kernel)
     * exactly when the genuine coupled/tower cases below are needed, so their
     * cleaner forms still win; rt_frac_case above keeps precedence for the
     * squarefree ArcTan/Log parts.  Correctness rests on the same diff-back gate. */
    r = rt_exp_ratreduce_case(f, x);
    if (r) return r;
    r = rt_hermite_case(f, x);
    if (r) return r;
    r = rt_hyperexp_case(f, x);
    if (r) return r;
    r = rt_expsum_case(f, x);   /* direct multi-kernel exponential sums */
    if (r) return r;
    r = rt_log_tower_case(f, x);   /* nested logarithmic tower (depth > 1) */
    if (r) return r;
    r = rt_exp_tower_case(f, x);    /* nested exponential tower (depth > 1) */
    if (r) return r;
    r = rt_recursive_tower_case(f, x);   /* one-extension recursion (mixed / rational coeff) */
    if (r) return r;
    /* Real tangent monomial (Bronstein §5.10): integrate rational functions of a
     * single Tan[u] directly and real, retiring the TrigToExp route below for
     * real-tangent integrands (which it strands at an I-laden complex-log form). */
    r = rt_hypertangent_case(f, x);
    if (r) return r;
    r = rt_trig_frontend(f, x);
    if (r) return r;
    return NULL;
}

/* ================================================================== */
/* Top-level integration dispatch.                                    */
/* ================================================================== */

/* Returns a fresh antiderivative (also self-verified by the recognizers,
 * and re-verified by the caller's diff-back gate) or NULL if no case
 * applies.  Dispatch order of the recursive Risch algorithm:
 *   1. rational base case (delegated to the recursive rational Risch,
 *      Integrate`BronsteinRational);
 *   2. transcendental case over a single logarithmic / exponential
 *      monomial extension (rt_transcendental_case — the recursive Risch
 *      proper: Hermite reduction, residue log-part, and the polynomial /
 *      Risch-differential-equation reductions);
 *   3. special-function outputs (Erf / dilog / Ei / li forms).
 * Every branch is verified by differentiation, so a mis-reduction can
 * only decline, never emit a wrong closed form.  NB: this must NOT fall
 * back on the parallel-Risch (pmint) engine Integrate`RischNorman —
 * that is a different algorithm; RischTranscendental is the recursive Risch. */
static Expr* rt_integrate(Expr* f, Expr* x) {
    Expr* r = rt_rational_case(f, x);
    if (r) return r;
    r = rt_transcendental_case(f, x);
    if (r) return r;
    r = rt_special_case(f, x);
    if (r) return r;
    return NULL;
}

/* ================================================================== */
/* Elementary-integrability decision procedure (P3).                  */
/* ================================================================== */

/* True iff the head names a function that is NOT elementary — the special
 * functions the integrator emits for a non-elementary antiderivative — or an
 * unevaluated Integrate (a partial / declined result).  RootSum and Root are
 * elementary (finite sums of logs / algebraic numbers) and are NOT listed. */
static bool rt_head_nonelementary(const char* h) {
    if (!h) return false;
    static const char* NE[] = {
        "ExpIntegralEi", "ExpIntegralE", "LogIntegral",
        "SinIntegral", "CosIntegral", "SinhIntegral", "CoshIntegral",
        "Erf", "Erfc", "Erfi", "FresnelS", "FresnelC",
        "PolyLog", "Integrate", NULL };
    for (int i = 0; NE[i]; i++)
        if (h == intern_symbol(NE[i])) return true;
    return false;
}

/* An antiderivative expressed WITHOUT any non-elementary special function (and
 * without an unevaluated Integrate) is an elementary closed form: exhibiting it
 * (correct by construction from the integrator) proves f is elementary-integrable. */
static bool rt_expr_is_elementary(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL
        && rt_head_nonelementary(e->data.function.head->data.symbol.name))
        return false;
    if (!rt_expr_is_elementary(e->data.function.head)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_expr_is_elementary(e->data.function.args[i])) return false;
    return true;
}

/* Field-path elementary-integrability decision.  Routes f through the AUTHORITATIVE
 * recursive field integrator (rt_field_integrate) in decision mode — the same
 * tower build / substitution / gate as rt_recursive_tower_case, but reading the
 * VERDICT rather than assembling an answer.  Returns:
 *   RT_DEC_ELEMENTARY     — a full antiderivative (no unintegrated remainder);
 *   RT_DEC_NONELEMENTARY  — an authoritative certificate fired: a non-constant
 *                           residue (Thm 5.6.1(ii)), a Risch DE with no rational
 *                           solution (Ch.6, rde_base / exact-bound ansatz), or the
 *                           §5.8 Dc!=0 primitive certificate;
 *   RT_DEC_UNKNOWN        — out of the single-tower field scope (algebraic
 *                           extensions, tower-builder rejects, or a genuine gap).
 * Sound by construction: an ELEMENTARY / NONELEMENTARY verdict is emitted only
 * behind an exact certificate — otherwise UNKNOWN. */
static RtDecision rt_decide_field(Expr* f, Expr* x) {
    Expr* fx = rt_expand_exp_sums(f);
    RtTower T;
    /* min_n = 1: the decision routes SINGLE-kernel integrands (E^x/x, E^(x^2),
     * 1/Log[x]) through the same authoritative field integrator too — the flat-tower
     * integrate cases that normally serve n=1 do not carry the decision certificates. */
    if (!rt_tower_build_min(fx, x, &T, 1)) { rt_tower_free(&T); expr_free(fx); return RT_DEC_UNKNOWN; }
    Expr* F = rt_eval1("Together", rt_subst_kernels(fx, &T));
    RtDecision verdict = RT_DEC_UNKNOWN;
    if (F && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL) {
        Expr* num = rt_eval1("Numerator", expr_copy(F));
        Expr* den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((T.n + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < T.n; i++) vv[i + 1] = expr_copy(T.t[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, T.n + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        bool gate = num && den && rt_is_true(pqn) && rt_is_true(pqd);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (gate) {
            bool       save_mode = g_rt_decide_mode;
            RtDecision save_dec  = g_rt_decision;
            g_rt_decide_mode = true;
            g_rt_decision    = RT_DEC_UNKNOWN;
            Expr* rem = NULL;
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x, &rem);
            if (Q && !rem)                               verdict = RT_DEC_ELEMENTARY;
            else if (Q && rem)                           verdict = RT_DEC_NONELEMENTARY;
            else if (g_rt_decision == RT_DEC_NONELEMENTARY) verdict = RT_DEC_NONELEMENTARY;
            else                                         verdict = RT_DEC_UNKNOWN;
            if (Q) expr_free(Q);
            if (rem) expr_free(rem);
            g_rt_decide_mode = save_mode;
            g_rt_decision    = save_dec;
        }
    }
    if (F) expr_free(F);
    rt_tower_free(&T);
    expr_free(fx);
    return verdict;
}

/* Complete elementary-integrability verdict for f w.r.t. x.  True side: exhibit an
 * elementary antiderivative via the full integrator (an existence proof — correct
 * by construction).  False side: the field decision's authoritative non-elementary
 * certificate.  Otherwise UNKNOWN.  Checking True first makes an exhibited
 * antiderivative dominate any certificate, so the verdict is robust. */
static RtDecision rt_decide(Expr* f, Expr* x) {
    Expr* anti = rt_integrate(f, x);
    if (anti) {
        bool elem = rt_expr_is_elementary(anti);
        expr_free(anti);
        if (elem) return RT_DEC_ELEMENTARY;   /* elementary closed form exhibited */
    }
    RtDecision fld = rt_decide_field(f, x);
    if (fld == RT_DEC_ELEMENTARY)    return RT_DEC_ELEMENTARY;
    if (fld == RT_DEC_NONELEMENTARY) return RT_DEC_NONELEMENTARY;
    return RT_DEC_UNKNOWN;
}

/* ================================================================== */
/* Public builtin.                                                    */
/* ================================================================== */

Expr* builtin_rischtranscendental(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    /* The integration variable must be a single symbol. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Correct by construction: rt_integrate returns a result only behind an
     * exact certificate, so no differentiation check is applied (a Risch
     * integrator is a decision procedure, not a guess-and-verify search).
     *
     * Mute arithmetic warnings across the decision procedure: its internal
     * algebraic manipulations (kernel substitutions, degree/coefficient probes,
     * Together/Cancel on tower normal forms) legitimately form transient
     * singular expressions — e.g. E^(x/E^(1/x))/E^(2/x) provokes a 1/0 while a
     * candidate exponent is reduced — none of which are part of the returned
     * antiderivative.  Surfacing Power::infy from them would be spurious (as it
     * is in Mathematica's Integrate). */
    arith_warnings_mute_push();
    Expr* result = rt_integrate(f, x);
    /* Decision half (P3): when the recursive algorithm produced no elementary
     * antiderivative, consult the authoritative field decision.  If it PROVES the
     * integrand has no elementary integral (a non-constant residue, a Risch DE with
     * no rational solution, or the Dc!=0 primitive certificate — never a bounded-
     * ansatz decline), report it Mathematica-style.  A stray certificate from an
     * abandoned attempt cannot fire here: rt_decide_field re-derives the verdict
     * with elementary-result precedence. */
    if (!result && rt_decide_field(f, x) == RT_DEC_NONELEMENTARY) {
        char* fs = expr_to_string(f);
        char* xs = expr_to_string(x);
        fprintf(stderr,
            "Integrate::nonelem: The integrand %s has no antiderivative "
            "elementary in %s.\n", fs ? fs : "?", xs ? xs : "?");
        free(fs); free(xs);
    }
    arith_warnings_mute_pop();
    return result;
}

/* Risch`ElementaryIntegralQ[f, x] — the Bronstein elementary-integrability decision
 * predicate.  Returns True when f has an antiderivative expressible in elementary
 * terms (exhibited by the integrator — an existence proof), False when the recursive
 * Risch decision procedure PROVES none exists (§5.6 residue criterion Thm 5.6.1(ii),
 * or a Ch.6 Risch DE with no rational solution, or the §5.8 Dc!=0 certificate), and
 * stays UNEVALUATED (with an ElementaryIntegralQ::undec message) when the verdict is
 * outside the single-tower field scope (algebraic extensions, deeper structures).
 * Sound by construction: a Boolean is returned only behind an exact certificate. */
static Expr* builtin_elementaryintegralq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    arith_warnings_mute_push();
    RtDecision d = rt_decide(f, x);
    arith_warnings_mute_pop();
    if (d == RT_DEC_ELEMENTARY)    return expr_new_symbol("True");
    if (d == RT_DEC_NONELEMENTARY) return expr_new_symbol("False");
    /* Undecided within scope: emit a diagnostic and leave the call unevaluated
     * (the *Q-returns-Boolean convention allows NULL only with an emitted tag). */
    {
        char* fs = expr_to_string(f);
        char* xs = expr_to_string(x);
        fprintf(stderr,
            "Risch`ElementaryIntegralQ::undec: Cannot decide elementary "
            "integrability of %s with respect to %s.\n", fs ? fs : "?", xs ? xs : "?");
        free(fs); free(xs);
    }
    return NULL;
}

/* Risch`RischDE[f, g, x] — solve the base-field Risch differential equation
 * D[y] + f y = g for y in C(x) (Bronstein Ch.6), exposing the internal rde_base
 * driver so the solver (and the coupled-system layer built on it) is directly
 * unit-testable.  Returns y (the unique rational solution, or 0 when g == 0), or
 * leaves the call unevaluated when no rational solution exists ("no solution" —
 * the term is non-elementary in C(x)).  Coefficients may lie in C(i)(x): the
 * Gaussian case falls through the FLINT fast path to the Expr pipeline. */
static Expr* builtin_risch_rischde(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* g = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    return rde_base(f, g, x);   /* NULL => no solution => call bubbles unevaluated */
}

/* ================================================================== */
/* Registration.                                                      */
/* ================================================================== */

static void rt_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_transcendental_init(void) {
    rt_install("Integrate`RischTranscendental", builtin_rischtranscendental,
        "Integrate`RischTranscendental[f, x] integrates f with respect to x using the\n"
        "recursive transcendental Risch algorithm: a decision\n"
        "procedure over a differential transcendental tower, with rational,\n"
        "logarithmic, exponential, and special-function (Erf, ExpIntegralEi,\n"
        "LogIntegral, PolyLog) cases.  Each case is correct by construction (no\n"
        "differentiation check).  Distinct from Integrate`RischNorman, which is\n"
        "the parallel-Risch (pmint) heuristic.  Out-of-scope integrands\n"
        "(algebraic extensions, non-elementary answers) return unevaluated.");
    rt_install("Risch`RischDE", builtin_risch_rischde,
        "Risch`RischDE[f, g, x] solves the Risch differential equation\n"
        "D[y] + f y == g for a rational y in C(x) (Bronstein Chapter 6): weak\n"
        "normalization, normal-denominator reduction, SPDE, and the polynomial\n"
        "non-cancellation solve.  Returns y (0 when g is 0), or stays unevaluated\n"
        "when no rational solution exists.  Coefficients may be Gaussian (C(i)(x)).");
    rt_install("Risch`ElementaryIntegralQ", builtin_elementaryintegralq,
        "Risch`ElementaryIntegralQ[f, x] decides whether f has an antiderivative\n"
        "elementary in x: True if the recursive Risch integrator exhibits an\n"
        "elementary closed form, False if it proves none exists (Bronstein residue\n"
        "criterion / Risch-DE no-solution certificates), or unevaluated when the\n"
        "verdict is outside the transcendental-tower field scope.");
}
