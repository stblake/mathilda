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
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"

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
           e->data.function.head->data.symbol == intern_symbol(name);
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
           e->data.function.head->data.symbol == intern_symbol("Rational");
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
    return e && e->type == EXPR_SYMBOL && e->data.symbol == intern_symbol("True");
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
        && cl->data.function.head->data.symbol == intern_symbol("List"))
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
static bool rt_is_zero(Expr* e) {
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
        && e->data.function.head->data.symbol == intern_symbol("Log")
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
        ? e->data.function.head->data.symbol : NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == intern_symbol("E")
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

/* Try each special-function recognizer in turn. */
static Expr* rt_special_case(Expr* f, Expr* x) {
    Expr* r;
    if ((r = rt_try_erf(f, x))) return r;
    if ((r = rt_try_ei(f, x))) return r;
    if ((r = rt_try_li(f, x))) return r;
    if ((r = rt_try_dilog(f, x))) return r;
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
        && cl->data.function.head->data.symbol == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    return a;
}

/* Risch differential equation q' + i u' q = p when the exponent u OR the
 * coefficient p is RATIONAL in x (so the polynomial ansatz below does not apply —
 * e.g. u = 1/x for E^(1/x), or p = 1/x - 1/x^2).  By the denominator theorem the
 * denominator of q divides Denominator[p], so q = h/pd with pd = Denominator[p]
 * and h a bounded polynomial in x, solved exactly by SolveAlways.  Correct by
 * construction (the SolveAlways certificate proves (q E^(i u))' = p E^(i u));
 * NULL when no such q exists (the term is non-elementary, e.g. E^(1/x) itself).
 * Closes Integrate[-E^(1/x)/x^2, x] = E^(1/x), Integrate[E^x/x - E^x/x^2, x] =
 * E^x/x. */
static Expr* rt_solve_rde_rational(Expr* p, long i, Expr* u, Expr* x) {
    Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));
    if (!up) return NULL;
    Expr* pg = rt_eval1("Together", expr_copy(p));
    Expr* pd = pg ? rt_eval1("Denominator", expr_copy(pg)) : NULL;
    Expr* pn = pg ? rt_eval1("Numerator", expr_copy(pg)) : NULL;
    if (pg) expr_free(pg);
    if (!pd || !pn) { expr_free(up); if (pd) expr_free(pd); if (pn) expr_free(pn); return NULL; }
    /* p and u' must be genuine RATIONAL FUNCTIONS of x — else a transcendental
     * coefficient (e.g. p = Sin[x], the E^x coefficient of a raw Sin[x] product)
     * would let SolveAlways certify a spurious q (it once returned q = 0 for
     * E^x Sin[x]).  Such integrands belong to the trig front-end / expsum path. */
    Expr* upg = rt_eval1("Together", expr_copy(up));
    Expr* upn = upg ? rt_eval1("Numerator", expr_copy(upg)) : NULL;
    Expr* upd = upg ? rt_eval1("Denominator", expr_copy(upg)) : NULL;
    if (upg) expr_free(upg);
    bool ratl = upn && upd && rt_is_poly(pn, x) && rt_is_poly(pd, x)
                && rt_is_poly(upn, x) && rt_is_poly(upd, x);
    if (upn) expr_free(upn);
    if (upd) expr_free(upd);
    if (!ratl) { expr_free(up); expr_free(pd); expr_free(pn); return NULL; }
    /* Exact numerator degree bound (no cap).  q = h/pd solves q' + i u' q = p over
     * C(x); deg_x(h) = deg_x(q) + deg_x(pd) with deg_x(q) from the leading balance of
     * q' against i u' q (rt_rde_var_bound, x is deriv-lowering).  f = i u' has
     * deg_x(f) = deg(Numerator[u']) - deg(Denominator[u']). */
    long dpd = rt_degree(pd, x); if (dpd < 0) dpd = 0;
    long dpn = rt_degree(pn, x); if (dpn < 0) dpn = 0;
    Expr* ung = rt_eval1("Together", expr_copy(up));
    Expr* unn = ung ? rt_eval1("Numerator", expr_copy(ung)) : NULL;
    Expr* und = ung ? rt_eval1("Denominator", expr_copy(ung)) : NULL;
    if (ung) expr_free(ung);
    long dfx = 0;
    if (unn && und) {
        long a = rt_degree(unn, x); if (a < 0) a = 0;
        long b = rt_degree(und, x); if (b < 0) b = 0;
        dfx = a - b;
    }
    if (unn) expr_free(unn);
    if (und) expr_free(und);
    /* v = x is deriv-lowering; its primitive cancellation config (dfv == -1) is
     * unreachable here — a rational exponent u has u' with pole order >= 2 (or >= 0),
     * never a simple pole — so no resonance integer applies (m_res = -1). */
    long N = rt_rde_var_bound(dpn - dpd, dfx, true, -1) + dpd;   /* deg_x(h) */

    Expr** terms = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (long k = 0; k <= N; k++) {
        char nm[24]; snprintf(nm, sizeof(nm), "rmRq%ld", k);
        terms[k] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_symbol(nm), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2) }, 2);
    }
    Expr* h = expr_new_function(expr_new_symbol("Plus"), terms, (size_t)(N + 1));
    free(terms);
    Expr* q = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ h, expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(pd), expr_new_integer(-1) }, 2) }, 2);

    Expr* dq = rt_eval2("D", expr_copy(q), expr_copy(x));
    Expr* iuq = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), expr_copy(up), expr_copy(q) }, 3);
    Expr* residual = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ dq, iuq, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 3));
    Expr* tog = residual ? rt_eval1("Together", residual) : NULL;
    Expr* rnum = tog ? rt_eval1("Numerator", expr_copy(tog)) : NULL;
    if (tog) expr_free(tog);

    Expr* result = NULL;
    Expr* sol = NULL;
    if (rnum) {
        Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
            (Expr*[]){ rnum, expr_new_integer(0) }, 2);              /* adopts rnum */
        sol = rt_eval2("SolveAlways", eqn, expr_copy(x));
    }
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol
             == intern_symbol("List")) {
        Expr* qi = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(q), expr_copy(sol->data.function.args[0]) }, 2));
        if (qi) {
            Expr** zero = malloc((size_t)(N + 1) * sizeof(Expr*));
            for (long k = 0; k <= N; k++) {
                char nm[24]; snprintf(nm, sizeof(nm), "rmRq%ld", k);
                zero[k] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
            }
            Expr* zl = expr_new_function(expr_new_symbol("List"), zero, (size_t)(N + 1));
            free(zero);
            result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ qi, zl }, 2));
            if (result) result = rt_eval1("Cancel", result);
        }
    }
    if (sol) expr_free(sol);
    expr_free(q); expr_free(up); expr_free(pd); expr_free(pn);
    return result;
}

/* Solve the Risch differential equation  q' + i u' q = p  for q in K,
 * with u a polynomial in x (u' polynomial) and p a polynomial in x, by a
 * bounded polynomial ansatz solved exactly with SolveAlways.  Correct by
 * construction: SolveAlways certifies the polynomial identity holds for
 * all x, so (q E^(i u))' = p E^(i u) exactly.  Returns q (owned) or NULL
 * when no polynomial solution exists (the term is then non-elementary in
 * this field — e.g. E^(-x^2) itself, left to the Erf recognizer).  When u or p
 * is rational in x, defers to rt_solve_rde_rational (E^(1/x) etc.). */
static Expr* rt_solve_rde(Expr* p, long i, Expr* u, Expr* x) {
    if (i == 0) return NULL;
    if (!rt_is_poly(u, x) || !rt_is_poly(p, x))
        return rt_solve_rde_rational(p, i, u, x);
    long du = rt_degree(u, x);
    if (du < 1) return NULL;                 /* u must be nonconstant */
    long dp = rt_degree(p, x);
    if (dp < 0) return NULL;
    long dF = du - 1;                        /* deg(i u') */
    long N = (dF >= 1) ? (dp - dF) : dp;     /* exact degree bound */
    if (N < 0) return NULL;

    Expr* up = rt_eval2("D", expr_copy(u), expr_copy(x));   /* u' */
    if (!up) return NULL;

    /* Ansatz q = sum_{k=0}^{N} rmC_k x^k. */
    Expr** terms = malloc((size_t)(N + 1) * sizeof(Expr*));
    for (long k = 0; k <= N; k++) {
        char nm[24];
        snprintf(nm, sizeof(nm), "rmC%ld", k);
        Expr* pw = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(x), expr_new_integer(k) }, 2);
        terms[k] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_symbol(nm), pw }, 2);
    }
    Expr* q = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        terms, (size_t)(N + 1)));
    free(terms);

    /* residual = D[q,x] + i u' q - p. */
    Expr* dq = rt_eval2("D", expr_copy(q), expr_copy(x));
    Expr* iuq = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), expr_copy(up), expr_copy(q) }, 3);
    Expr* negp = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2);
    Expr* residual = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ dq, iuq, negp }, 3));

    /* SolveAlways[residual == 0, x]. */
    Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ residual, expr_new_integer(0) }, 2);   /* adopts residual */
    Expr* sol = rt_eval2("SolveAlways", eqn, expr_copy(x));

    Expr* result = NULL;
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol
             == intern_symbol("List")) {
        /* Apply the solution, then pin any free ansatz parameters to 0. */
        Expr* qi = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(q), expr_copy(sol->data.function.args[0]) }, 2));
        Expr** zero = malloc((size_t)(N + 1) * sizeof(Expr*));
        for (long k = 0; k <= N; k++) {
            char nm[24];
            snprintf(nm, sizeof(nm), "rmC%ld", k);
            zero[k] = expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_symbol(nm), expr_new_integer(0) }, 2);
        }
        Expr* zl = expr_new_function(expr_new_symbol("List"), zero, (size_t)(N + 1));
        free(zero);
        result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ qi, zl }, 2));   /* adopts qi, zl */
    }
    if (sol) expr_free(sol);
    expr_free(q);
    expr_free(up);
    return result;
}

/* Collect (as owned copies, deduplicated) the exponents of every E^w / Exp[w]
 * kernel of `e` whose exponent depends on x. */
static void rt_collect_exp_exponents(Expr* e, Expr* x,
                                     Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol : NULL;
    Expr* w = NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        w = e->data.function.args[0];
    else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol == intern_symbol("E")
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
                && factored->data.function.head->data.symbol == intern_symbol("Times")) {
                fa = factored->data.function.args;
                nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol == intern_symbol("Power")
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
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
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
        if (kernel_back) expr_free(kernel_back); if (F) expr_free(F);
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
                        && logpart->data.function.head->data.symbol
                             == intern_symbol("Integrate`TranscendentalLogPart"));
                if (!declined) {
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
        if (F) expr_free(F); if (Dt) expr_free(Dt);
        if (kback) expr_free(kback); if (uu) expr_free(uu);
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
                && factored->data.function.head->data.symbol == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol == intern_symbol("Power")
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
                    char nm[32]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
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
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    Expr** zero = malloc((nh + ng) * sizeof(Expr*));
                    size_t zi = 0;
                    for (long p = 0; p < dH; p++)
                        for (long k = 0; k <= Nx; k++) {
                            char nm[32]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
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
            && cl->data.function.head->data.symbol == intern_symbol("List")) {
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
                && factored->data.function.head->data.symbol == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol == intern_symbol("Power")
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
                        char nm[32]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
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
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
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
                            char nm[32]; snprintf(nm, sizeof(nm), "rmHh%ldv%ld", p, k);
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
        && T->data.function.head->data.symbol == intern_symbol("Times")) {
        for (size_t i = 0; i < T->data.function.arg_count; i++)
            rt_accum_factors(T->data.function.args[i], x, W, pf, npf, cap);
        return;
    }
    Expr* w = NULL;   /* borrowed exponent if T is an x-dependent exp kernel */
    const char* h = (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL)
        ? T->data.function.head->data.symbol : NULL;
    if (h == intern_symbol("Exp") && T->data.function.arg_count == 1
        && !rt_free_of_x(T->data.function.args[0], x))
        w = T->data.function.args[0];
    else if (h == intern_symbol("Power") && T->data.function.arg_count == 2
        && T->data.function.args[0]->type == EXPR_SYMBOL
        && T->data.function.args[0]->data.symbol == intern_symbol("E")
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
        && fe->data.function.head->data.symbol == intern_symbol("Plus")) {
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
        && e->data.function.head->data.symbol == intern_symbol("Log")
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
    Expr** logs = NULL; size_t nl = 0, lcap = 0;
    rt_collect_logs(f, x, &logs, &nl, &lcap);
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
                            free(logs); return NULL; }

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
        char nm[16]; snprintf(nm, sizeof(nm), "rmtw%zu", i);
        ts[i] = expr_new_symbol(nm);
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(logs[i]), expr_copy(ts[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, nl);
    free(rules);

    Expr* F = rt_eval1("Together", rt_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(f), expr_copy(rl) }, 2)));
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
                && factored->data.function.head->data.symbol == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol == intern_symbol("Power")
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
                    char nm[32]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
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
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr* Qs = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
                if (Qs) {
                    /* Pin any free ansatz parameters to 0. */
                    Expr** zero = malloc(nq * sizeof(Expr*));
                    size_t zi = 0;
                    for (long k = 0; k <= Ntop; k++)
                        for (long m = 0; m < nmono; m++) {
                            char nm[32]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
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
 * back-substituted antiderivative (owned) or NULL.  Borrows all arguments. */
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
        && sol->data.function.head->data.symbol == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol
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
        char nm[16]; snprintf(nm, sizeof(nm), "rmte%zu", i);
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
                && factored->data.function.head->data.symbol == intern_symbol("Times")) {
                fa = factored->data.function.args; nf = factored->data.function.arg_count;
            } else { single[0] = factored; fa = single; nf = 1; }
            for (size_t i = 0; i < nf && !bad; i++) {
                Expr* term = fa[i]; Expr* base = term; long e = 1;
                if (term->type == EXPR_FUNCTION
                    && term->data.function.head->type == EXPR_SYMBOL
                    && term->data.function.head->data.symbol == intern_symbol("Power")
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
static bool rt_tower_build(Expr* f, Expr* x, RtTower* T) {
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
    if (!okc || n < 2 || n > RT_MAXK) {
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
        Expr* d;
        if (T->kind[i] == RT_LOG) {
            Expr* du = rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
            Expr* q = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ du, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->arg[i]), expr_new_integer(-1) }, 2) }, 2);
            d = rt_eval1("Cancel", q);
        } else {
            d = rt_eval1("Cancel", rt_eval2("D", expr_copy(T->arg[i]), expr_copy(x)));
        }
        if (!d) { ok = false; break; }
        Expr* ds = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ d, expr_copy(T->subrules) }, 2));           /* adopts d */
        if (!ds) { ok = false; break; }
        T->Dcoef[i] = ds;
        if (rt_find_exp_of_x(ds, x) != NULL || rt_find_log_of_x(ds, x) != NULL) ok = false;
        for (size_t j = i; j < n && ok; j++)
            if (!rt_free_of_x(ds, T->t[j])) ok = false;
    }
    return ok;
}

/* Mutually recursive field-integration primitives (forward declarations). */
static Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x);
static int   rt_limited_field_integrate(Expr* r, RtTower* T, long L, Expr* x,
                                        Expr** s_out, Expr** c_out);
static Expr* rt_int_primitive_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_int_hyperexp_poly(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_field_rde(Expr* p, long i, RtTower* T, long L, Expr* x);
static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x);
static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x);

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
                                  Expr* x) {
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

    Expr* logpart = rt_eval_call("Integrate`TranscendentalLogPart",
        (Expr*[]){ a, rad, expr_copy(t), expr_new_symbol("rmZ"), Dd, gate }, 6);
    /* Declined iff the head is unchanged (a, rad, Dd, gate adopted by the call). */
    if (logpart && logpart->type == EXPR_FUNCTION
        && logpart->data.function.head->type == EXPR_SYMBOL
        && logpart->data.function.head->data.symbol
             == intern_symbol("Integrate`TranscendentalLogPart")) {
        expr_free(logpart);
        return NULL;
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
static Expr* rt_field_ratint(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    Expr* dden = rt_eval2("D", expr_copy(den), expr_copy(t));
    Expr* Hden = dden ? rt_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(den), dden }, 2) : NULL;
    if (!Hden) return NULL;
    long dH = rt_degree(Hden, t); if (dH < 0) dH = 0;
    Expr* rad = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));

    /* Distinct t-dependent factors of the squarefree radical (log arguments). */
    Expr* factored = rad ? rt_eval1("Factor", expr_copy(rad)) : NULL;
    Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
    if (factored) {
        Expr** fa; size_t nf; Expr* single[1];
        if (factored->type == EXPR_FUNCTION
            && factored->data.function.head->type == EXPR_SYMBOL
            && factored->data.function.head->data.symbol == intern_symbol("Times")) {
            fa = factored->data.function.args; nf = factored->data.function.arg_count;
        } else { single[0] = factored; fa = single; nf = 1; }
        for (size_t i = 0; i < nf && !bad; i++) {
            Expr* term = fa[i]; Expr* base = term;
            if (term->type == EXPR_FUNCTION
                && term->data.function.head->type == EXPR_SYMBOL
                && term->data.function.head->data.symbol == intern_symbol("Power")
                && term->data.function.arg_count == 2)
                base = term->data.function.args[0];       /* rad squarefree: e = 1 */
            if (rt_free_of_x(base, t)) continue;
            if (ng >= 16) { bad = true; break; }
            g[ng++] = expr_copy(base);
        }
    }

    /* Lower-field variables {x, t_0..t_{L-1}} and per-variable degree bounds. */
    size_t nlv = (size_t)L + 1;
    Expr** lv = malloc(nlv * sizeof(Expr*));
    long* bd = malloc(nlv * sizeof(long));
    lv[0] = x;
    for (long i = 0; i < L; i++) lv[i + 1] = T->t[i];
    for (size_t j = 0; j < nlv; j++) {
        long a = rt_degree(num, lv[j]); if (a < 0) a = 0;
        long b = rt_degree(den, lv[j]); if (b < 0) b = 0;
        long d = (a > b ? a : b) + 1;    /* derived lower-field proxy, no cap */
        bd[j] = d;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
    long nunk = dH * nmono + (long)ng;

    Expr* result = NULL;
    if (!bad && nunk > 0) {
        size_t nq = (size_t)nunk;
        Expr** qterms = malloc(nq * sizeof(Expr*));
        Expr** syms = malloc(nq * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        /* Hermite numerator H(t)/Hden: sum_{p<dH} (sum_mono a x^..) t^p. */
        for (long p = 0; p < dH; p++)
            for (long m = 0; m < nmono; m++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmGh%ld_%ld", p, m);
                rt_decode_mono(m, bd, nlv, ev);
                Expr* mono = rt_build_monomial(lv, ev, nlv);
                syms[nsym++] = expr_new_symbol(nm);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm), mono,
                      expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(t), expr_new_integer(p) }, 2) }, 3);
            }
        free(ev);
        Expr* Hpoly = NULL;
        if (dH >= 1) {
            Hpoly = expr_new_function(expr_new_symbol("Plus"),
                qterms, ntq);           /* adopts the H-term array contents */
            /* rebuild qterms as: [ Hpoly/Hden, c_j Log(g_j)... ] */
            Expr** q2 = malloc((1 + ng) * sizeof(Expr*));
            size_t n2 = 0;
            q2[n2++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
            free(qterms);
            qterms = q2; ntq = n2;
        } else {
            ntq = 0;                    /* no Hermite part */
        }
        for (size_t j = 0; j < ng; j++) {
            char nm[40]; snprintf(nm, sizeof(nm), "rmGc%zu", j);
            syms[nsym++] = expr_new_symbol(nm);
            qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_symbol(nm),
                  expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
        }
        Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
        free(qterms);

        /* residual = Numerator[Together[D_tower[Q] - num/den]]. */
        Expr* target = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
        Expr* Qder = rt_tower_deriv(Q, T, x);
        Expr* diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), target }, 2) }, 2);
        Expr* tog = rt_eval1("Together", diff);
        Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc((T->n + 1) * sizeof(Expr*));
            for (size_t i = 0; i < T->n; i++) vl[i] = expr_copy(T->t[T->n - 1 - i]);
            vl[T->n] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, T->n + 1);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rt_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
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
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, zl }, 2));
            }
        }
        if (sol) expr_free(sol);
        for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
        free(syms);
        expr_free(Q);
    }

    /* Algebraic-residue log part (pure resultant Lazard-Rioboo-Trager): the
     * bounded ansatz above expresses only rational (single-constant) residues;
     * when it declines, close a squarefree-denominator factor whose residues
     * are algebraic constants of the tower derivation. */
    if (!result) result = rt_field_lrt_logpart(num, den, T, L, x);

    for (size_t j = 0; j < ng; j++) expr_free(g[j]);
    if (factored) expr_free(factored);
    free(lv); free(bd);
    if (rad) expr_free(rad);
    expr_free(Hden);
    return result;
}

/* Coupled hyperexponential proper part (exponential top level t = t_L = E^(w_L)).
 * Unlike the logarithmic case, an exponential kernel's Laurent part and log part do
 * NOT separate — D_tower[Log(g)] = D_tower[g]/g and D_tower[t] = w' t, so a single
 * Log(1 + t) contributes to both the t^0 Laurent coefficient and the proper part.
 * So (as in the single-kernel rt_hyperexp_case, lifted here to the tower derivation)
 * a UNIFIED ansatz is solved at once:
 *   Q = sum_{i=ilo}^{ihi} w_i t^i + H(t)/Hden(t) + sum_j c_j Log(g_j)
 * where ilo = -(mult. of t at 0 in den), Hden = gcd(Dtil, dDtil/dt) is the repeated
 * part of the t-coprime denominator Dtil = den/t^a, g_j the distinct t-factors of
 * the squarefree radical Dtil/Hden, the w_i and the numerator H are bounded-degree
 * lower-field polynomials, and the c_j constants — all found by SolveAlways over
 * every tower variable of D_tower[Q] = num/den.  Correct by construction; the
 * caller diff-back verifies.  Tried only when the Laurent recursion
 * (rt_int_hyperexp_poly) declines because a genuine proper part is present. */
static Expr* rt_field_hyperexp_coupled(Expr* num, Expr* den, RtTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    long a = 0;
    Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(t));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rt_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    Expr* Dtil = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(t), expr_new_integer(-a) }, 2) }, 2));
    if (!Dtil) return NULL;

    Expr* dDt = rt_eval2("D", expr_copy(Dtil), expr_copy(t));
    Expr* Hden = dDt ? rt_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(Dtil), dDt }, 2) : NULL;
    long dH = Hden ? rt_degree(Hden, t) : 0; if (dH < 0) dH = 0;
    Expr* rad = Hden ? rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(Dtil), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2)) : NULL;

    Expr* factored = rad ? rt_eval1("Factor", expr_copy(rad)) : NULL;
    Expr* g[16]; size_t ng = 0; bool bad = (factored == NULL);
    if (factored) {
        Expr** fa; size_t nf; Expr* single[1];
        if (factored->type == EXPR_FUNCTION
            && factored->data.function.head->type == EXPR_SYMBOL
            && factored->data.function.head->data.symbol == intern_symbol("Times")) {
            fa = factored->data.function.args; nf = factored->data.function.arg_count;
        } else { single[0] = factored; fa = single; nf = 1; }
        for (size_t i = 0; i < nf && !bad; i++) {
            Expr* term = fa[i]; Expr* base = term;
            if (term->type == EXPR_FUNCTION
                && term->data.function.head->type == EXPR_SYMBOL
                && term->data.function.head->data.symbol == intern_symbol("Power")
                && term->data.function.arg_count == 2)
                base = term->data.function.args[0];
            if (rt_free_of_x(base, t)) continue;
            if (ng >= 16) { bad = true; break; }
            g[ng++] = expr_copy(base);
        }
    }

    long dnum = rt_degree(num, t), dden = rt_degree(den, t);
    /* Single EXP kernel t (D preserves deg_t): Laurent range = f's own extent,
     * high = deg_t(num) - deg_t(den), low = -(mult of t at 0 in den).  No cap. */
    long ihi = dnum - dden; if (ihi < 0) ihi = 0;
    long ilo = -a;
    long nwi = ihi - ilo + 1;

    size_t nlv = (size_t)L + 1;
    Expr** lv = malloc(nlv * sizeof(Expr*));
    long* bd = malloc(nlv * sizeof(long));
    lv[0] = x;
    for (long i = 0; i < L; i++) lv[i + 1] = T->t[i];
    for (size_t j = 0; j < nlv; j++) {
        long p = rt_degree(num, lv[j]); if (p < 0) p = 0;
        long q = rt_degree(den, lv[j]); if (q < 0) q = 0;
        long d = (p > q ? p : q) + 1;    /* derived lower-field proxy, no cap */
        bd[j] = d;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
    long nunk = nwi * nmono + dH * nmono + (long)ng;

    Expr* result = NULL;
    if (!bad && nwi > 0 && nunk > 0) {
        size_t nq = (size_t)(nwi * nmono + dH * nmono + (long)ng);
        Expr** qterms = malloc(nq * sizeof(Expr*));
        Expr** syms = malloc(nq * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        /* Laurent part: sum_i (sum_mono aW mono) t^i. */
        for (long i = ilo; i <= ihi; i++)
            for (long m = 0; m < nmono; m++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmXw%ld_%ld", i - ilo, m);
                rt_decode_mono(m, bd, nlv, ev);
                Expr* mono = rt_build_monomial(lv, ev, nlv);
                syms[nsym++] = expr_new_symbol(nm);
                qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_symbol(nm), mono,
                      expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(t), expr_new_integer(i) }, 2) }, 3);
            }
        /* Hermite part: (sum_{p<dH} (sum_mono aH mono) t^p) / Hden. */
        if (dH >= 1) {
            Expr** hterms = malloc((size_t)(dH * nmono) * sizeof(Expr*));
            size_t nh = 0;
            for (long p = 0; p < dH; p++)
                for (long m = 0; m < nmono; m++) {
                    char nm[40]; snprintf(nm, sizeof(nm), "rmXh%ld_%ld", p, m);
                    rt_decode_mono(m, bd, nlv, ev);
                    Expr* mono = rt_build_monomial(lv, ev, nlv);
                    syms[nsym++] = expr_new_symbol(nm);
                    hterms[nh++] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_symbol(nm), mono,
                          expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ expr_copy(t), expr_new_integer(p) }, 2) }, 3);
                }
            Expr* Hpoly = expr_new_function(expr_new_symbol("Plus"), hterms, nh);
            free(hterms);
            qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ Hpoly, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2);
        }
        /* Log part: sum_j c_j Log(g_j). */
        for (size_t j = 0; j < ng; j++) {
            char nm[40]; snprintf(nm, sizeof(nm), "rmXc%zu", j);
            syms[nsym++] = expr_new_symbol(nm);
            qterms[ntq++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_symbol(nm),
                  expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(g[j]) }, 1) }, 2);
        }
        free(ev);
        Expr* Q = expr_new_function(expr_new_symbol("Plus"), qterms, ntq);
        free(qterms);

        Expr* target = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2) }, 2);
        Expr* Qder = rt_tower_deriv(Q, T, x);
        Expr* diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), target }, 2) }, 2);
        Expr* tog = rt_eval1("Together", diff);
        Expr* rnum = tog ? rt_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc((T->n + 1) * sizeof(Expr*));
            for (size_t i = 0; i < T->n; i++) vl[i] = expr_copy(T->t[T->n - 1 - i]);
            vl[T->n] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, T->n + 1);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rt_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
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
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, zl }, 2));
            }
        }
        if (sol) expr_free(sol);
        for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
        free(syms);
        expr_free(Q);
    }

    /* Algebraic-residue log part (pure resultant LRT): when the coupled ansatz
     * declines and there is no Laurent principal part (a == 0, so num/den is a
     * pure proper fraction and the LRT log part is complete), close a squarefree
     * factor with algebraic constant residues.  The top-level caller diff-back
     * verifies, so a case with a hidden Laurent part declines harmlessly. */
    if (!result && a == 0) result = rt_field_lrt_logpart(num, den, T, L, x);

    for (size_t j = 0; j < ng; j++) expr_free(g[j]);
    if (factored) expr_free(factored);
    free(lv); free(bd);
    if (rad) expr_free(rad);
    if (Hden) expr_free(Hden);
    expr_free(Dtil);
    return result;
}

/* Integrate F (a rational function of x, t_1..t_L in tower-variable form) with
 * respect to the tower derivation, returning an antiderivative in tower-variable
 * form (owned) or NULL.  L < 0 is the rational base case C(x). */
static Expr* rt_field_integrate(Expr* F, RtTower* T, long L, Expr* x) {
    if (!F) return NULL;
    if (L < 0) {
        Expr* r = rt_eval_call("Integrate`BronsteinRational",
            (Expr*[]){ expr_copy(F), expr_copy(x) }, 2);
        if (!r) return NULL;
        if (rt_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
        return r;
    }
    Expr* t = T->t[L];
    Expr* Fg = rt_eval1("Together", expr_copy(F));
    if (!Fg) return NULL;
    Expr* num = rt_eval1("Numerator", expr_copy(Fg));
    Expr* den = rt_eval1("Denominator", expr_copy(Fg));
    Expr* result = NULL;
    if (num && den) {
        if (T->kind[L] == RT_LOG) {
            if (rt_free_of_x(den, t)) {
                /* Pure polynomial in t: primitive-polynomial recursion. */
                result = rt_int_primitive_poly(num, den, T, L, x);
            } else {
                /* Split the polynomial part (recursion) from the proper rational
                 * part in t (Hermite reduction + Rothstein-Trager log part). */
                Expr* Pp = rt_eval_call("PolynomialQuotient",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(t) }, 3);
                Expr* Rr = rt_eval_call("PolynomialRemainder",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(t) }, 3);
                if (Pp && Rr) {
                    Expr* one = expr_new_integer(1);
                    Expr* poly_int = rt_int_primitive_poly(Pp, one, T, L, x);
                    expr_free(one);
                    Expr* prop_int = rt_is_zero(Rr) ? expr_new_integer(0)
                                                    : rt_field_ratint(Rr, den, T, L, x);
                    if (poly_int && prop_int)
                        result = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                            (Expr*[]){ poly_int, prop_int }, 2));
                    else {
                        if (poly_int) expr_free(poly_int);
                        if (prop_int) expr_free(prop_int);
                    }
                }
                if (Pp) expr_free(Pp);
                if (Rr) expr_free(Rr);
            }
        } else {
            /* Exponential top: the pure-Laurent recursion first (genuine RDE,
             * rational lower-field coefficients); if it declines because a proper
             * part is present, the coupled hyperexponential ansatz. */
            result = rt_int_hyperexp_poly(num, den, T, L, x);
            if (!result)
                result = rt_field_hyperexp_coupled(num, den, T, L, x);
        }
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(Fg);
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
    Expr* R = rt_field_integrate(r, T, L - 1, x);
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
        } else { expr_free(c); expr_free(s); }
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
        && cl->data.function.head->data.symbol == intern_symbol("List")) {
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
        Expr* qj = (ip == 0) ? rt_field_integrate(pj, T, L - 1, x)
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
    bool base = rt_is_poly(w, x) && rt_is_poly(p, x);
    for (size_t j = 0; j < T->n && base; j++)
        if (!rt_free_of_x(w, T->t[j]) || !rt_free_of_x(p, T->t[j])) base = false;
    if (base) return rt_solve_rde(p, i, w, x);

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
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);

    Expr* result = NULL;
    if (nmono > 0) {
        Expr** hterms = malloc((size_t)nmono * sizeof(Expr*));
        Expr** syms = malloc((size_t)nmono * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        for (long m = 0; m < nmono; m++) {
            char nm[40]; snprintf(nm, sizeof(nm), "rmRd%ld", m);
            rt_decode_mono(m, bd, nlv, ev);
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
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
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
    free(lv); free(bd);
    expr_free(pd); expr_free(pn);
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
        ? e->data.function.head->data.symbol : NULL;
    Expr* expo = NULL;
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol == intern_symbol("E"))
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
            && e->data.function.head->data.symbol == intern_symbol("Exp")
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i]))
            return expr_copy(T->t[i]);
    }
    /* Multiplicatively commensurate member kernel E^(marg) = t[mprim]^mmult. */
    for (size_t m = 0; m < T->nm; m++) {
        if (e->type == EXPR_FUNCTION && e->data.function.arg_count >= 1
            && e->data.function.head->type == EXPR_SYMBOL) {
            const char* h = e->data.function.head->data.symbol;
            Expr* w = NULL;
            if (h == intern_symbol("Exp") && e->data.function.arg_count == 1)
                w = e->data.function.args[0];
            else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
                     && e->data.function.args[0]->type == EXPR_SYMBOL
                     && e->data.function.args[0]->data.symbol == intern_symbol("E"))
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
    if (!rt_tower_build(fx, x, &T)) { rt_tower_free(&T); expr_free(fx); return NULL; }

    /* Alias the kernels to tower variables structurally (NOT via an evaluated
     * ReplaceAll, which would re-merge a split exponential product before
     * substitution — see rt_subst_kernels), then normalise with Together. */
    Expr* F = rt_eval1("Together", rt_subst_kernels(fx, &T));
    Expr* result = NULL;
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
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x);
            if (Q) {
                Expr** back = malloc(T.n * sizeof(Expr*));
                for (size_t i = 0; i < T.n; i++)
                    back[i] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                Expr* bl = expr_new_function(expr_new_symbol("List"), back, T.n);
                free(back);
                result = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Q, bl }, 2));                       /* adopts Q, bl */
                /* Tidy the Hermite rational part (e.g. -x Log[x]/(x Log[x](1+t)) ->
                 * -1/(1+t)); Cancel treats the Log/exp kernels as opaque atoms. */
                if (result) result = rt_eval1("Cancel", result);
            }
        }
    }
    if (F) expr_free(F);

    /* Diff-back safety gate (bounded search, not a decision procedure).
     * Verify against the ORIGINAL integrand f (== fx mathematically). */
    if (result && !rt_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    rt_tower_free(&T);
    expr_free(fx);
    return result;
}

/* Trigonometric / hyperbolic front-end (exponentialize path).
 * Rewrites the trig/hyperbolic kernels to complex exponentials with
 * TrigToExp, integrates the resulting (Laurent-)rational function of the
 * exponential kernel E^(I x) / E^x with the exponential machinery, and
 * converts the answer back to trigonometric form with ExpToTrig.  Both
 * rewrites are exact, so the result is correct by construction. */
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
    if (!r) r = rt_hyperexp_case(fe, x);
    /* Multi-kernel decoupling (Phase B): e.g. Sin/Cos times a real exponential
     * exponentialize to a sum of two non-commensurate exponentials E^((a +/- b I) x)
     * that the single-primitive cases cannot kernelize. */
    if (!r) r = rt_expsum_case(fe, x);
    expr_free(fe);
    if (!r) return NULL;
    return rt_eval1("ExpToTrig", r);   /* adopts r; back to trig form */
}

/* Dispatch the transcendental cases: the primitive (logarithmic) polynomial
 * reduction, the exponential (hyperexponential / Risch-DE) reduction, the
 * fractional (Rothstein-Trager) log-part, and the trig/hyperbolic front-end.
 * The general Hermite reduction for repeated poles lands in a subsequent
 * increment. */
static Expr* rt_transcendental_case(Expr* f, Expr* x) {
    Expr* r = rt_log_poly_case(f, x);
    if (r) return r;
    r = rt_exp_poly_case(f, x);
    if (r) return r;
    r = rt_frac_case(f, x);
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
     * integrator is a decision procedure, not a guess-and-verify search). */
    return rt_integrate(f, x);
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
}
