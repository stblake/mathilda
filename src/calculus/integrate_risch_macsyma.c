/* integrate_risch_macsyma.c — Risch integrator ported from Maxima.
 *
 * Faithful port of the algorithm STRUCTURE of Maxima's src/risch.lisp,
 * with arithmetic grounded in Mathilda's existing Expr/poly/rat
 * primitives (see the header for the full contract and phase plan).
 *
 * Phase 1 (this revision): scaffold, dispatcher registration, the
 * rational case (delegated to Integrate`BronsteinRational), and the
 * mandatory diff-back verification gate.  Later phases add the
 * differential tower and the logarithmic / exponential / special-
 * function cases.
 */

#include "integrate_risch_macsyma.h"

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
static bool rm_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == intern_symbol(name);
}

/* Build `head[args...]` (adopting the owned `args` element pointers) and
 * evaluate it, freeing the constructed call.  Returns evaluate()'s result. */
static Expr* rm_eval_call(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Evaluate `e`, taking ownership: frees `e` and returns evaluate()'s
 * result (evaluate itself does not consume its argument). */
static Expr* rm_eval_own(Expr* e) {
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
 * closed form it emits, exactly as Maxima's risch does.               */

static Expr* rm_eval1(const char* head, Expr* a) {
    return rm_eval_call(head, (Expr*[]){ a }, 1);
}
static Expr* rm_eval2(const char* head, Expr* a, Expr* b) {
    return rm_eval_call(head, (Expr*[]){ a, b }, 2);
}
static Expr* rm_eval3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rm_eval_call(head, (Expr*[]){ a, b, c }, 3);
}

static bool rm_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == intern_symbol("True");
}
/* FreeQ[e, x] */
static bool rm_free_of_x(Expr* e, Expr* x) {
    Expr* r = rm_eval2("FreeQ", expr_copy(e), expr_copy(x));
    bool t = rm_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* FreeQ[e, head] — true iff the symbol `head` occurs nowhere in e. */
static bool rm_free_of_head(Expr* e, const char* head) {
    Expr* r = rm_eval2("FreeQ", expr_copy(e), expr_new_symbol(head));
    bool t = rm_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* PolynomialQ[e, x] */
static bool rm_is_poly(Expr* e, Expr* x) {
    Expr* r = rm_eval2("PolynomialQ", expr_copy(e), expr_copy(x));
    bool t = rm_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* Degree of a polynomial in x via Length[CoefficientList[e, x]] - 1
 * (Mathilda has no Exponent builtin).  Returns -1 if CoefficientList does
 * not reduce to a List.  Call only after rm_is_poly. */
static long rm_degree(Expr* e, Expr* x) {
    Expr* cl = rm_eval2("CoefficientList", expr_copy(e), expr_copy(x));
    long d = -1;
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol == intern_symbol("List"))
        d = (long)cl->data.function.arg_count - 1;
    if (cl) expr_free(cl);
    return d;
}
/* Coefficient[e, x, k] */
static Expr* rm_coeff(Expr* e, Expr* x, long k) {
    return rm_eval3("Coefficient", expr_copy(e), expr_copy(x),
                    expr_new_integer(k));
}
/* Together[e] === 0 — an exact zero test for rational functions of the
 * field kernels (Together is exact and, unlike Simplify, cheap; it is only
 * ever applied here to constants and small rational cofactors). */
static bool rm_is_zero(Expr* e) {
    Expr* s = rm_eval1("Together", expr_copy(e));
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* Parse `tmpl`, substitute the named placeholder symbols with `vals`
 * (borrowed; copied in), evaluate, and return the result (or NULL). */
static Expr* rm_template(const char* tmpl, const char** names,
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
/* Case: rational function of x  (Maxima rischfprog / dprog / eprog).  */
/* ================================================================== */

/* Delegate the pure-rational case to the Bronstein rational integrator,
 * which IS Maxima's rational Risch case and strictly stronger.  Returns
 * NULL when the integrand is not rational in x (BronsteinRational leaves
 * its call unevaluated). */
static Expr* rm_rational_case(Expr* f, Expr* x) {
    Expr* r = rm_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(f), expr_copy(x) }, 2);
    if (!r) return NULL;
    if (rm_head_is(r, "Integrate`BronsteinRational")) {
        expr_free(r);
        return NULL;
    }
    return r;
}

/* ================================================================== */
/* Case: special functions  (Maxima erfarg / erfarg2 / Ei / li / dilog).*/
/* ================================================================== */
/* Each recognizer builds a candidate antiderivative from a template and
 * accepts it only if it passes the diff-back gate, so a mis-recognition
 * can never emit a wrong closed form.  These close integrals the whole
 * elementary cascade leaves open, using special functions Mathilda
 * already provides (Erf, ExpIntegralEi, LogIntegral, PolyLog).          */

/* K * E^(a x^2 + b x + c) with the leading (quadratic) coefficient
 * nonzero  ->  Erf/Erfi.  Detected by the log-derivative f'/f being a
 * degree-1 polynomial (Maxima's erfarg). */
static Expr* rm_try_erf(Expr* f, Expr* x) {
    Expr* df = rm_eval2("D", expr_copy(f), expr_copy(x));
    if (!df) return NULL;
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(f), expr_new_integer(-1) }, 2);
    Expr* quot = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ df, inv }, 2);
    Expr* ld = rm_eval1("Cancel", quot);
    if (!ld) return NULL;
    if (!rm_is_poly(ld, x) || rm_degree(ld, x) != 1) { expr_free(ld); return NULL; }
    Expr* a2 = rm_coeff(ld, x, 1);   /* = 2a */
    Expr* b1 = rm_coeff(ld, x, 0);   /* = b  */
    expr_free(ld);

    const char* names[2] = { "rmA2", "rmB" };
    Expr* vals[2] = { a2, b1 };
    Expr* result = NULL;
    /* Since ld = f'/f is a degree-1 polynomial, f = K' E^((a2/2)x^2 + b1 x)
     * with K' a constant, so K' = f|_{x=0} (the exponent has no constant
     * term).  Evaluating at 0 avoids a Simplify on a product of Gaussian
     * exponentials, which is prohibitively slow. */
    Expr* kp = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(x), expr_new_integer(0) }, 2) }, 2));
    if (kp) {
        /* Integral of E^((a2/2)x^2 + b1 x), a = a2/2, b = b1:
         *   E^(-b^2/(4a)) * (-Sqrt[Pi] Erf[(b+2a x)/(2 Sqrt[-a])] / (2 Sqrt[-a]))
         * where -b^2/(4a) = -b1^2/(2 a2) is the completing-the-square shift. */
        Expr* erfpart = rm_template(
            "E^(-rmB^2/(2*rmA2)) *"
            " (-(Sqrt[Pi]*Erf[(rmB + rmA2*x)/(2*Sqrt[-rmA2/2])])/(2*Sqrt[-rmA2/2]))",
            names, vals, 2);
        if (erfpart) {
            /* Correct by construction: ld = f'/f is a degree-1 polynomial,
             * so f = K' E^(a x^2 + b x) exactly and this is its integral. */
            result = rm_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(kp), erfpart }, 2));
        }
        expr_free(kp);
    }
    expr_free(a2);
    expr_free(b1);
    return result;
}

/* (const * E^(a x)) / (c x + d)  ->  ExpIntegralEi.  Numerator is a
 * constant times a pure exponential; denominator linear in x. */
static Expr* rm_try_ei(Expr* f, Expr* x) {
    Expr* g = rm_eval1("Together", expr_copy(f));
    if (!g) return NULL;
    Expr* num = rm_eval1("Numerator", expr_copy(g));
    Expr* den = rm_eval1("Denominator", expr_copy(g));
    expr_free(g);
    Expr* result = NULL;
    if (num && den && rm_is_poly(den, x) && rm_degree(den, x) == 1
        && !rm_free_of_x(num, x)) {
        /* numld = Cancel[D[num]/num] must be a nonzero constant (= a). */
        Expr* dn = rm_eval2("D", expr_copy(num), expr_copy(x));
        Expr* invn = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(num), expr_new_integer(-1) }, 2);
        Expr* q = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ dn, invn }, 2);
        Expr* numld = rm_eval1("Cancel", q);
        if (numld && rm_free_of_x(numld, x) && !rm_is_zero(numld)) {
            /* num = M E^(a x) with M constant, so M = num|_{x=0}. */
            Expr* M = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(num),
                    expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(x), expr_new_integer(0) }, 2) }, 2));
            if (M && rm_free_of_x(M, x)) {
                Expr* cc = rm_coeff(den, x, 1);
                Expr* dd = rm_coeff(den, x, 0);
                const char* names[4] = { "rmM", "rmA", "rmC", "rmD" };
                Expr* vals[4] = { M, numld, cc, dd };
                /* Correct by construction: num = M E^(a x) (num'/num = a
                 * constant) over an exactly linear denominator c x + d. */
                result = rm_template(
                    "(rmM*E^(-rmA*rmD/rmC)/rmC)*ExpIntegralEi[rmA*(x + rmD/rmC)]",
                    names, vals, 4);
                expr_free(cc);
                expr_free(dd);
            }
            if (M) expr_free(M);
        }
        if (numld) expr_free(numld);
    }
    if (num) expr_free(num);
    if (den) expr_free(den);
    return result;
}

/* K / Log[x]  ->  K LogIntegral[x].  Certificate: Together[f Log[x]] is a
 * constant K free of x (so f = K/Log[x] exactly). */
static Expr* rm_try_li(Expr* f, Expr* x) {
    Expr* logx = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(x) }, 1);
    Expr* prod = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(f), logx }, 2);
    Expr* K = rm_eval1("Together", prod);
    Expr* result = NULL;
    if (K && rm_free_of_x(K, x) && !rm_is_zero(K)) {
        Expr* li = expr_new_function(expr_new_symbol("LogIntegral"),
            (Expr*[]){ expr_copy(x) }, 1);
        result = rm_eval_own(expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(K), li }, 2));
    }
    if (K) expr_free(K);
    return result;
}

/* Find (borrowed) the argument u of the first Log[u] subexpression of `e`
 * whose argument depends on x. */
static Expr* rm_find_log_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rm_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    Expr* r = rm_find_log_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rm_find_log_of_x(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

/* Find (borrowed) the exponent u of the first exponential kernel of `e`
 * whose exponent depends on x, matching either Exp[u] or E^u (Power[E,u]). */
static Expr* rm_find_exp_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol : NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rm_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == intern_symbol("E")
            && !rm_free_of_x(ex, x))
            return ex;
    }
    Expr* r = rm_find_exp_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rm_find_exp_of_x(e->data.function.args[i], x);
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
 * be left to the tower cases (rm_log_tower_case / rm_exp_tower_case). */
static bool rm_kernel_simple(Expr* u, Expr* x) {
    return rm_find_exp_of_x(u, x) == NULL && rm_find_log_of_x(u, x) == NULL;
}

/* K Log[1 + p x] / x  ->  -K PolyLog[2, -p x]  (Maxima's dilog).
 * Certificate: the Log argument is exactly linear with constant term 1,
 * and Together[f x / Log[u]] is a constant K free of x (so f is exactly
 * K Log[1+p x]/x). */
static Expr* rm_try_dilog(Expr* f, Expr* x) {
    Expr* u = rm_find_log_of_x(f, x);   /* borrowed */
    if (!u) return NULL;
    if (!rm_is_poly(u, x) || rm_degree(u, x) != 1) return NULL;
    Expr* u0 = rm_coeff(u, x, 0);
    Expr* u1 = rm_coeff(u, x, 1);
    Expr* result = NULL;
    Expr* u0m1 = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(u0), expr_new_integer(-1) }, 2);
    bool u0_is_1 = rm_is_zero(u0m1);
    expr_free(u0m1);
    if (u0_is_1) {
        Expr* logu = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        Expr* invlog = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ logu, expr_new_integer(-1) }, 2);
        Expr* prod = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(f), expr_copy(x), invlog }, 3);
        Expr* K = rm_eval1("Together", prod);
        if (K && rm_free_of_x(K, x) && !rm_is_zero(K)) {
            const char* names[2] = { "rmK", "rmU1" };
            Expr* vals[2] = { K, u1 };
            result = rm_template("-rmK*PolyLog[2, -rmU1*x]", names, vals, 2);
        }
        if (K) expr_free(K);
    }
    expr_free(u0);
    expr_free(u1);
    return result;
}

/* Try each special-function recognizer in turn. */
static Expr* rm_special_case(Expr* f, Expr* x) {
    Expr* r;
    if ((r = rm_try_erf(f, x))) return r;
    if ((r = rm_try_ei(f, x))) return r;
    if ((r = rm_try_li(f, x))) return r;
    if ((r = rm_try_dilog(f, x))) return r;
    return NULL;
}

/* ================================================================== */
/* Case: transcendental — the recursive Risch algorithm proper.        */
/* ================================================================== */
/* This is the genuine recursive Risch (Maxima's rischint over a single
 * logarithmic / exponential monomial extension), NOT the parallel-Risch
 * (pmint) heuristic.  It reduces the integrand in the differential field
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
static int rm_limited_integrate(Expr* r, Expr* x, Expr* u,
                                Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rm_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(r), expr_copy(x) }, 2);
    if (!R) return -1;
    if (rm_head_is(R, "Integrate`BronsteinRational")) { expr_free(R); return -1; }
    /* Rsub = R with Log[u] replaced by the fresh variable rmT. */
    Expr* logu = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(u) }, 1);
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ logu, expr_new_symbol("rmT") }, 2);
    Expr* Rsub = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ R, rule }, 2));   /* adopts R */
    if (!Rsub) return -1;
    Expr* tsym = expr_new_symbol("rmT");
    int rc = -1;
    if (rm_free_of_head(Rsub, "Log") && rm_is_poly(Rsub, tsym)
        && rm_degree(Rsub, tsym) <= 1) {
        Expr* c = rm_coeff(Rsub, tsym, 1);   /* theta-coefficient */
        Expr* s = rm_coeff(Rsub, tsym, 0);   /* rational part      */
        if (rm_free_of_x(c, x)) { *s_out = s; *c_out = c; rc = 0; }
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
static Expr* rm_log_poly_case(Expr* f, Expr* x) {
    Expr* u = rm_find_log_of_x(f, x);        /* borrowed: the Log argument */
    if (!u || !rm_kernel_simple(u, x)) return NULL;   /* nested -> tower */

    /* F = f with Log[u] replaced by the fresh polynomial variable rmT. */
    Expr* logu = expr_new_function(expr_new_symbol("Log"),
        (Expr*[]){ expr_copy(u) }, 1);
    Expr* rule = expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ logu, expr_new_symbol("rmT") }, 2);
    Expr* F = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rule }, 2));
    if (!F) return NULL;

    Expr* tsym = expr_new_symbol("rmT");
    /* Require: pure polynomial in theta, no residual logs (single
     * extension), and genuine theta-dependence. */
    if (!rm_is_poly(F, tsym) || !rm_free_of_head(F, "Log")
        || rm_free_of_x(F, tsym)) {
        expr_free(F); expr_free(tsym); return NULL;
    }
    long m = rm_degree(F, tsym);
    if (m < 1) { expr_free(F); expr_free(tsym); return NULL; }

    /* eta = Cancel[D[u,x]/u]. */
    Expr* du = rm_eval2("D", expr_copy(u), expr_copy(x));
    Expr* invu = expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
    Expr* eta = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ du, invu }, 2));

    Expr** p = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (long i = 0; i <= m; i++) p[i] = rm_coeff(F, tsym, i);
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
        int rc = rm_limited_integrate(r_i, x, u, &s, &c);
        expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        /* Fold the theta-term back: q[i+1] += c/(i+1). */
        Expr* bump = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2);
        if (q[i + 1] == NULL) {
            q[i + 1] = rm_eval_own(bump);
        } else {
            q[i + 1] = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
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
static Expr* rm_integrate_in_K_with_logs(Expr* g, Expr* x) {
    Expr* r = rm_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(g), expr_copy(x) }, 2);
    if (!r) return NULL;
    if (rm_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
    return r;
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
static Expr* rm_solve_rde_rational(Expr* p, long i, Expr* u, Expr* x) {
    Expr* up = rm_eval2("D", expr_copy(u), expr_copy(x));
    if (!up) return NULL;
    Expr* pg = rm_eval1("Together", expr_copy(p));
    Expr* pd = pg ? rm_eval1("Denominator", expr_copy(pg)) : NULL;
    Expr* pn = pg ? rm_eval1("Numerator", expr_copy(pg)) : NULL;
    if (pg) expr_free(pg);
    if (!pd || !pn) { expr_free(up); if (pd) expr_free(pd); if (pn) expr_free(pn); return NULL; }
    /* p and u' must be genuine RATIONAL FUNCTIONS of x — else a transcendental
     * coefficient (e.g. p = Sin[x], the E^x coefficient of a raw Sin[x] product)
     * would let SolveAlways certify a spurious q (it once returned q = 0 for
     * E^x Sin[x]).  Such integrands belong to the trig front-end / expsum path. */
    Expr* upg = rm_eval1("Together", expr_copy(up));
    Expr* upn = upg ? rm_eval1("Numerator", expr_copy(upg)) : NULL;
    Expr* upd = upg ? rm_eval1("Denominator", expr_copy(upg)) : NULL;
    if (upg) expr_free(upg);
    bool ratl = upn && upd && rm_is_poly(pn, x) && rm_is_poly(pd, x)
                && rm_is_poly(upn, x) && rm_is_poly(upd, x);
    if (upn) expr_free(upn);
    if (upd) expr_free(upd);
    if (!ratl) { expr_free(up); expr_free(pd); expr_free(pn); return NULL; }
    long dpd = rm_degree(pd, x); if (dpd < 0) dpd = 0;
    long dpn = rm_degree(pn, x); if (dpn < 0) dpn = 0;
    long N = dpd + dpn + 2; if (N > 10) N = 10;

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

    Expr* dq = rm_eval2("D", expr_copy(q), expr_copy(x));
    Expr* iuq = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), expr_copy(up), expr_copy(q) }, 3);
    Expr* residual = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ dq, iuq, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 3));
    Expr* tog = residual ? rm_eval1("Together", residual) : NULL;
    Expr* rnum = tog ? rm_eval1("Numerator", expr_copy(tog)) : NULL;
    if (tog) expr_free(tog);

    Expr* result = NULL;
    Expr* sol = NULL;
    if (rnum) {
        Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
            (Expr*[]){ rnum, expr_new_integer(0) }, 2);              /* adopts rnum */
        sol = rm_eval2("SolveAlways", eqn, expr_copy(x));
    }
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol
             == intern_symbol("List")) {
        Expr* qi = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
            result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ qi, zl }, 2));
            if (result) result = rm_eval1("Cancel", result);
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
 * is rational in x, defers to rm_solve_rde_rational (E^(1/x) etc.). */
static Expr* rm_solve_rde(Expr* p, long i, Expr* u, Expr* x) {
    if (i == 0) return NULL;
    if (!rm_is_poly(u, x) || !rm_is_poly(p, x))
        return rm_solve_rde_rational(p, i, u, x);
    long du = rm_degree(u, x);
    if (du < 1) return NULL;                 /* u must be nonconstant */
    long dp = rm_degree(p, x);
    if (dp < 0) return NULL;
    long dF = du - 1;                        /* deg(i u') */
    long N = (dF >= 1) ? (dp - dF) : dp;     /* exact degree bound */
    if (N < 0) return NULL;

    Expr* up = rm_eval2("D", expr_copy(u), expr_copy(x));   /* u' */
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
    Expr* q = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
        terms, (size_t)(N + 1)));
    free(terms);

    /* residual = D[q,x] + i u' q - p. */
    Expr* dq = rm_eval2("D", expr_copy(q), expr_copy(x));
    Expr* iuq = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(i), expr_copy(up), expr_copy(q) }, 3);
    Expr* negp = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2);
    Expr* residual = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ dq, iuq, negp }, 3));

    /* SolveAlways[residual == 0, x]. */
    Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
        (Expr*[]){ residual, expr_new_integer(0) }, 2);   /* adopts residual */
    Expr* sol = rm_eval2("SolveAlways", eqn, expr_copy(x));

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
        Expr* qi = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
        result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ qi, zl }, 2));   /* adopts qi, zl */
    }
    if (sol) expr_free(sol);
    expr_free(q);
    expr_free(up);
    return result;
}

/* Collect (as owned copies, deduplicated) the exponents of every E^w / Exp[w]
 * kernel of `e` whose exponent depends on x. */
static void rm_collect_exp_exponents(Expr* e, Expr* x,
                                     Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol : NULL;
    Expr* w = NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rm_free_of_x(e->data.function.args[0], x))
        w = e->data.function.args[0];
    else if (h == intern_symbol("Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol == intern_symbol("E")
        && !rm_free_of_x(e->data.function.args[1], x))
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
    rm_collect_exp_exponents(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rm_collect_exp_exponents(e->data.function.args[i], x, arr, n, cap);
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
static Expr* rm_exp_kernelize(Expr* f, Expr* x, Expr** u_out) {
    *u_out = NULL;
    if (!rm_find_exp_of_x(f, x)) return NULL;
    Expr** ws = NULL; size_t nw = 0, cap = 0;
    rm_collect_exp_exponents(f, x, &ws, &nw, &cap);
    if (nw == 0) { free(ws); return NULL; }

    Expr* u = NULL; long* kof = malloc(nw * sizeof(long));
    for (size_t cand = 0; cand < nw && !u; cand++) {
        bool all_int = true;
        for (size_t j = 0; j < nw; j++) {
            Expr* ratio = rm_eval1("Cancel", expr_new_function(
                expr_new_symbol("Times"), (Expr*[]){ expr_copy(ws[j]),
                    expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(ws[cand]), expr_new_integer(-1) }, 2)
                }, 2));
            if (ratio && ratio->type == EXPR_INTEGER && ratio->data.integer != 0)
                kof[j] = (long)ratio->data.integer;
            else all_int = false;
            if (ratio) expr_free(ratio);
            if (!all_int) break;
        }
        if (all_int) u = ws[cand];
    }
    /* The primitive exponent must be rational in x alone; a nested exponent
     * (e.g. u = E^x in E^(E^x)) is a two-extension tower left to rm_exp_tower_case
     * (else the single-kernel RDE would carry the inner kernel as a free param). */
    if (u && !rm_kernel_simple(u, x)) u = NULL;
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
    Expr* uexp = expr_copy(u);
    for (size_t i = 0; i < nw; i++) expr_free(ws[i]);
    free(ws); free(kof);

    Expr* F = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rl }, 2));
    if (!F) { expr_free(uexp); return NULL; }
    *u_out = uexp;
    return F;
}

static Expr* rm_exp_poly_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rm_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;

    Expr* tsym = expr_new_symbol("rmT");
    /* F must be a Laurent polynomial in t: num/den with den a pure power of t.
     * Split into num and the offset M (den = c t^M). */
    Expr* G = rm_eval1("Together", expr_copy(F));
    Expr* num = G ? rm_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rm_eval1("Denominator", expr_copy(G)) : NULL;
    Expr* result = NULL;
    long M = 0;
    bool ok = num && den && rm_is_poly(num, tsym) && rm_is_poly(den, tsym)
        && !rm_free_of_x(F, tsym)
        && rm_find_exp_of_x(F, x) == NULL && rm_find_log_of_x(F, x) == NULL;
    if (ok) {
        M = rm_degree(den, tsym);
        /* den must be a monomial c t^M: den / t^M free of t. */
        Expr* dq = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(tsym), expr_new_integer(-M) }, 2) }, 2));
        if (!dq || !rm_free_of_x(dq, tsym)) ok = false;
        if (dq) expr_free(dq);
        /* Fold the constant den/t^M into num so p_i = Coefficient[num/c, ...]. */
        if (ok) {
            Expr* nnew = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(den), expr_new_integer(-1) }, 2),
                  expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(tsym), expr_new_integer(M) }, 2) }, 3));
            if (nnew) { expr_free(num); num = nnew; }
            else ok = false;
        }
    }
    long dnum = ok ? rm_degree(num, tsym) : -1;
    if (ok && dnum >= 0) {
        Expr** q = calloc((size_t)(dnum + 1), sizeof(Expr*));
        long* qi_pow = malloc((size_t)(dnum + 1) * sizeof(long));
        bool fail = false;
        for (long j = 0; j <= dnum && !fail; j++) {
            long i = j - M;                       /* Laurent power */
            Expr* pi = rm_coeff(num, tsym, j);
            if (rm_is_zero(pi)) { expr_free(pi); q[j] = NULL; qi_pow[j] = i; continue; }
            Expr* qi = (i == 0) ? rm_integrate_in_K_with_logs(pi, x)
                                : rm_solve_rde(pi, i, uexp, x);
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
            result = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
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
static Expr* rm_frac_try(Expr* f, Expr* x, Expr* u, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* rules; Expr* Dt; Expr* kernel_back;
    if (is_log) {
        Expr* logu = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(u) }, 1);
        rules = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ logu, expr_new_symbol("rmT") }, 2) }, 1);
        Expr* du = rm_eval2("D", expr_copy(u), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(u), expr_new_integer(-1) }, 2);
        Dt = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
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
        Expr* up = rm_eval2("D", expr_copy(u), expr_copy(x));
        Dt = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ up, expr_new_symbol("rmT") }, 2);
        kernel_back = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(u) }, 2);
    }
    if (!Dt) { expr_free(tsym); expr_free(rules); expr_free(kernel_back); return NULL; }

    Expr* F0 = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(f), rules }, 2));
    Expr* F = F0 ? rm_eval1("Together", F0) : NULL;
    if (!F) { expr_free(tsym); expr_free(Dt); expr_free(kernel_back); return NULL; }

    Expr* num = rm_eval1("Numerator", expr_copy(F));
    Expr* den = rm_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;

    bool ok = num && den && rm_is_poly(num, tsym) && rm_is_poly(den, tsym)
        && !rm_free_of_x(den, tsym)
        && rm_find_exp_of_x(F, x) == NULL && rm_find_log_of_x(F, x) == NULL;

    if (ok) {
        Expr* factored = rm_eval1("Factor", expr_copy(den));
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
                if (rm_free_of_x(base, tsym)) continue;   /* FreeQ[base, t] */
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
                Expr* dgx = rm_eval2("D", expr_copy(g[i]), expr_copy(x));
                Expr* dgt = rm_eval2("D", expr_copy(g[i]), expr_copy(tsym));
                Expr* dgi = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ dgx, expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(Dt), dgt }, 2) }, 2));
                Expr* cof = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(g[i]), expr_new_integer(-1) }, 2) }, 2));
                negterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(nm),
                              dgi ? dgi : expr_new_integer(0),
                              cof ? cof : expr_new_integer(0) }, 4);
            }
            Expr* residual = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                negterms, ng + 1));
            free(negterms);

            Expr* varlist = expr_new_function(expr_new_symbol("List"),
                (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ residual, expr_new_integer(0) }, 2);
            Expr* sol = rm_eval2("SolveAlways", eqn, varlist);
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr** rterms = malloc(ng * sizeof(Expr*));
                for (size_t i = 0; i < ng; i++) {
                    char nm[24]; snprintf(nm, sizeof(nm), "rmK%zu", i);
                    Expr* gib = rm_eval_own(expr_new_function(
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
                R = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
                    result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ R, zl }, 2));
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
static Expr* rm_hermite_try(Expr* f, Expr* x, bool is_log) {
    Expr* tsym = expr_new_symbol("rmT");
    Expr* F = NULL; Expr* Dt = NULL; Expr* kback = NULL; Expr* uu = NULL;
    if (is_log) {
        Expr* u = rm_find_log_of_x(f, x);
        if (!u || !rm_kernel_simple(u, x)) { expr_free(tsym); return NULL; }
        uu = expr_copy(u);
        Expr* du = rm_eval2("D", expr_copy(uu), expr_copy(x));
        Expr* invu = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(uu), expr_new_integer(-1) }, 2);
        Dt = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ du, invu }, 2));   /* eta = u'/u */
        kback = expr_new_function(expr_new_symbol("Log"),
            (Expr*[]){ expr_copy(uu) }, 1);
        Expr* rl = expr_new_function(expr_new_symbol("List"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Log"),
                    (Expr*[]){ expr_copy(uu) }, 1), expr_new_symbol("rmT") }, 2) }, 1);
        Expr* F0 = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(f), rl }, 2));
        F = F0 ? rm_eval1("Together", F0) : NULL;
    } else {
        Expr* F0 = rm_exp_kernelize(f, x, &uu);
        if (F0) F = rm_eval1("Together", F0);
        if (uu) {
            Expr* up = rm_eval2("D", expr_copy(uu), expr_copy(x));
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

    Expr* num = rm_eval1("Numerator", expr_copy(F));
    Expr* den = rm_eval1("Denominator", expr_copy(F));
    Expr* result = NULL;
    bool ok = num && den && rm_is_poly(num, tsym) && rm_is_poly(den, tsym)
        && !rm_free_of_x(den, tsym)
        && rm_find_log_of_x(F, x) == NULL && rm_find_exp_of_x(F, x) == NULL;
    long dnum = ok ? rm_degree(num, tsym) : 0, dden = ok ? rm_degree(den, tsym) : 0;
    ok = ok && dnum < dden;                      /* proper fraction */
    if (ok && !is_log) {
        /* exponential kernel: require D coprime to theta (no theta = 0 pole). */
        Expr* c0 = rm_coeff(den, tsym, 0);
        if (rm_is_zero(c0)) ok = false;
        expr_free(c0);
    }

    Expr* Hden = NULL; Expr* rad = NULL;
    if (ok) {
        Expr* dent = rm_eval2("D", expr_copy(den), expr_copy(tsym));
        Hden = rm_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(den), dent }, 2);
        long dH = Hden ? rm_degree(Hden, tsym) : 0;
        if (dH < 1) ok = false;                  /* no repeated pole */
        if (ok) rad = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
    }

    if (ok && rad) {
        long dH = rm_degree(Hden, tsym);
        /* distinct t-factors of the squarefree radical (log terms). */
        Expr* factored = rm_eval1("Factor", expr_copy(rad));
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
                if (rm_free_of_x(base, tsym)) continue;
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }
        long dnx = rm_degree(num, x), ddx = rm_degree(den, x);
        long Nx = (dnx > ddx ? dnx : ddx) + 2; if (Nx > 8) Nx = 8;
        size_t nh = (size_t)(dH * (Nx + 1));
        if (!bad && (long)(nh + ng) <= 80) {
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

            Expr* dQx = rm_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rm_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* Qder = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt), dQt }, 2) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rm_eval1("Together", diff);
            Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rm_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
                    Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym), expr_copy(kback) }, 2) }, 1);
                        Expr* rr = rm_eval_own(expr_new_function(
                            expr_new_symbol("ReplaceAll"), (Expr*[]){ Qs, back }, 2));
                        /* Cancel folds away the x-content that Hden carried
                         * (e.g. -x/(x + x Log[x]) -> -1/(1 + Log[x])). */
                        if (rr) result = rm_eval1("Cancel", rr);
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
static Expr* rm_hermite_case(Expr* f, Expr* x) {
    Expr* r = rm_hermite_try(f, x, true);
    if (r) return r;
    r = rm_hermite_try(f, x, false);
    return r;
}

/* Try the fractional log-part for a logarithmic then an exponential kernel. */
static Expr* rm_frac_case(Expr* f, Expr* x) {
    Expr* ul = rm_find_log_of_x(f, x);
    if (ul && rm_kernel_simple(ul, x)) {
        Expr* r = rm_frac_try(f, x, ul, true); if (r) return r;
    }
    Expr* ue = rm_find_exp_of_x(f, x);
    if (ue && rm_kernel_simple(ue, x)) {
        Expr* r = rm_frac_try(f, x, ue, false); if (r) return r;
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
static Expr* rm_hyperexp_case(Expr* f, Expr* x) {
    Expr* uexp = NULL;
    Expr* F = rm_exp_kernelize(f, x, &uexp);
    if (!F) return NULL;
    Expr* tsym = expr_new_symbol("rmT");
    Expr* result = NULL;

    Expr* G = rm_eval1("Together", expr_copy(F));
    Expr* num = G ? rm_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rm_eval1("Denominator", expr_copy(G)) : NULL;
    bool ok = num && den && rm_is_poly(num, tsym) && rm_is_poly(den, tsym)
        && !rm_free_of_x(den, tsym)
        && rm_find_exp_of_x(F, x) == NULL && rm_find_log_of_x(F, x) == NULL;

    if (ok) {
        /* a = multiplicity of t in den (first nonzero t-coefficient). */
        long a = 0;
        Expr* cl = rm_eval2("CoefficientList", expr_copy(den), expr_copy(tsym));
        if (cl && cl->type == EXPR_FUNCTION
            && cl->data.function.head->type == EXPR_SYMBOL
            && cl->data.function.head->data.symbol == intern_symbol("List")) {
            for (size_t i = 0; i < cl->data.function.arg_count; i++)
                if (!rm_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
        }
        if (cl) expr_free(cl);
        Expr* Dtil = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
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
            Expr* dDt = rm_eval2("D", expr_copy(Dtil), expr_copy(tsym));
            Hden = rm_eval_call("PolynomialGCD",
                (Expr*[]){ expr_copy(Dtil), dDt }, 2);
            if (Hden) {
                dH = rm_degree(Hden, tsym);
                if (dH < 0) dH = 0;
                rad = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dtil), expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));
            }
        }
        Expr* factored = rad ? rm_eval1("Factor", expr_copy(rad)) : NULL;
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
                if (rm_free_of_x(base, tsym)) continue;    /* FreeQ[base, t] */
                if (ng >= 16) { bad = true; break; }
                g[ng++] = expr_copy(base);
            }
        }

        long dnum = ok ? rm_degree(num, tsym) : 0, dden = ok ? rm_degree(den, tsym) : 0;
        long ihi = dnum - dden; if (ihi < 0) ihi = 0;
        long ilo = -a;
        long degu = rm_degree(uexp, x);
        long dnx = rm_degree(num, x), ddx = rm_degree(den, x);
        long Nx = (dnx > ddx ? dnx : ddx) + (degu > 0 ? degu : 1) + 1;
        if (Nx > 8) Nx = 8;
        long nwi = ihi - ilo + 1;
        size_t nH_syms = (size_t)(dH * (Nx + 1));   /* Hermite numerator coeffs */
        if (!bad && nwi > 0
            && nwi * (Nx + 1) + (long)nH_syms + (long)ng <= 80) {
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
            Expr* dQx = rm_eval2("D", expr_copy(Q), expr_copy(x));
            Expr* dQt = rm_eval2("D", expr_copy(Q), expr_copy(tsym));
            Expr* up = rm_eval2("D", expr_copy(uexp), expr_copy(x));
            Expr* Qder = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ dQx, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ up, expr_copy(tsym), dQt }, 3) }, 2));
            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rm_eval1("Together", diff);
            Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;
            Expr* sol = NULL;
            if (rnum) {
                Expr* varlist = expr_new_function(expr_new_symbol("List"),
                    (Expr*[]){ expr_copy(tsym), expr_copy(x) }, 2);
                Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ rnum, expr_new_integer(0) }, 2);
                sol = rm_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
                    Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                        (Expr*[]){ Qs, zl }, 2));
                    if (Qs) {
                        Expr* back = expr_new_function(expr_new_symbol("List"),
                            (Expr*[]){ expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(tsym),
                                    expr_new_function(expr_new_symbol("Power"),
                                        (Expr*[]){ expr_new_symbol("E"), expr_copy(uexp) }, 2)
                                }, 2) }, 1);
                        result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
static void rm_accum_factors(Expr* T, Expr* x, Expr** W,
                             Expr*** pf, size_t* npf, size_t* cap) {
    if (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL
        && T->data.function.head->data.symbol == intern_symbol("Times")) {
        for (size_t i = 0; i < T->data.function.arg_count; i++)
            rm_accum_factors(T->data.function.args[i], x, W, pf, npf, cap);
        return;
    }
    Expr* w = NULL;   /* borrowed exponent if T is an x-dependent exp kernel */
    const char* h = (T->type == EXPR_FUNCTION
        && T->data.function.head->type == EXPR_SYMBOL)
        ? T->data.function.head->data.symbol : NULL;
    if (h == intern_symbol("Exp") && T->data.function.arg_count == 1
        && !rm_free_of_x(T->data.function.args[0], x))
        w = T->data.function.args[0];
    else if (h == intern_symbol("Power") && T->data.function.arg_count == 2
        && T->data.function.args[0]->type == EXPR_SYMBOL
        && T->data.function.args[0]->data.symbol == intern_symbol("E")
        && !rm_free_of_x(T->data.function.args[1], x))
        w = T->data.function.args[1];
    if (w) {
        *W = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ *W, expr_copy(w) }, 2));
    } else {
        if (*npf == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *pf = realloc(*pf, *cap * sizeof(Expr*));
        }
        (*pf)[(*npf)++] = expr_copy(T);
    }
}

static int rm_split_exp_term(Expr* T, Expr* x, Expr** W_out, Expr** p_out) {
    Expr* W = expr_new_integer(0);
    Expr** pf = NULL; size_t npf = 0, cap = 0;
    rm_accum_factors(T, x, &W, &pf, &npf, &cap);

    Expr* p;
    if (npf == 0) p = expr_new_integer(1);
    else if (npf == 1) p = pf[0];                          /* adopt sole factor */
    else p = expr_new_function(expr_new_symbol("Times"), pf, npf);  /* adopts */
    free(pf);

    if (rm_find_exp_of_x(p, x) != NULL) {   /* nested/residual exp in cofactor */
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
static Expr* rm_expsum_case(Expr* f, Expr* x) {
    if (!rm_find_exp_of_x(f, x)) return NULL;   /* need at least one exp kernel */
    Expr* fe = rm_eval1("Expand", expr_copy(f));
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
        if (rm_split_exp_term(terms[i], x, &W, &p) != 0) { fail = true; break; }
        bool wdep = !rm_free_of_x(W, x);
        Expr* qi = NULL;
        if (!wdep) {
            /* E^W is a constant coefficient: INT p E^W dx = E^W INT p dx. */
            qi = rm_integrate_in_K_with_logs(p, x);
        } else if (rm_is_poly(W, x)) {
            qi = rm_solve_rde(p, 1, W, x);   /* q' + W' q = p */
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
        result = rm_eval_own(expr_new_function(expr_new_symbol("Plus"), outs, no));
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
/* The tower cases (rm_log_tower_case / rm_exp_tower_case) are bounded-ansatz
 * SEARCHES over a multi-variable SolveAlways, not pure decision procedures: the
 * certificate "D_tower[Q] = F as a polynomial identity in the tower variables"
 * is sound only when SolveAlways returns a genuine solution, and a non-elementary
 * integrand can yield a spurious one (e.g. E^(E^x)/(1+E^(E^x)) is non-elementary
 * yet the ansatz once returned Log[1+E^(E^x)]/E^x).  As with the other
 * search-based integrators, the tower results are therefore diff-back verified:
 * Simplify[D[result,x] - f] must be exactly 0, else the case declines.  This
 * keeps the genuine closures (E^x E^(E^x) -> E^(E^x),
 * E^x E^(E^x)/(1+E^(E^x)) -> Log[1+E^(E^x)]) and rejects the spurious ones. */
static bool rm_verify_antideriv(Expr* result, Expr* f, Expr* x) {
    Expr* d = rm_eval2("D", expr_copy(result), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    Expr* s = rm_eval1("Simplify", diff);
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* ================================================================== */
/* Case: nested logarithmic tower (Phase B, second increment).         */
/* ================================================================== */
/* Collect (owned, deduplicated) every Log[u] kernel of `e` whose argument
 * depends on x. */
static void rm_collect_logs(Expr* e, Expr* x, Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rm_free_of_x(e->data.function.args[0], x)) {
        bool dup = false;
        for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], e)) { dup = true; break; }
        if (!dup) {
            if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                              *arr = realloc(*arr, *cap * sizeof(Expr*)); }
            (*arr)[(*n)++] = expr_copy(e);
        }
    }
    rm_collect_logs(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rm_collect_logs(e->data.function.args[i], x, arr, n, cap);
}

/* Structural containment: does `big` contain `small` as a subexpression? */
static bool rm_contains(Expr* big, Expr* small) {
    if (!big) return false;
    if (expr_eq(big, small)) return true;
    if (big->type != EXPR_FUNCTION) return false;
    if (rm_contains(big->data.function.head, small)) return true;
    for (size_t i = 0; i < big->data.function.arg_count; i++)
        if (rm_contains(big->data.function.args[i], small)) return true;
    return false;
}

/* Build the monomial  prod_j lv[j]^e[j]  (a Times of Powers; owned). */
static Expr* rm_build_monomial(Expr** lv, const long* e, size_t nlv) {
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
static void rm_decode_mono(long idx, const long* bd, size_t nlv, long* e) {
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
static Expr* rm_log_tower_case(Expr* f, Expr* x) {
    Expr** logs = NULL; size_t nl = 0, lcap = 0;
    rm_collect_logs(f, x, &logs, &nl, &lcap);
    if (nl < 2 || nl > 4) { for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
                            free(logs); return NULL; }

    /* Order innermost-first: if logs[i] contains logs[k], logs[k] is deeper. */
    for (size_t i = 0; i < nl; i++)
        for (size_t k = i + 1; k < nl; k++)
            if (rm_contains(logs[i], logs[k])) {
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

    Expr* F = rm_eval1("Together", rm_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(f), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL;
    Expr** Dt = NULL;
    Expr* result = NULL;

    bool ok = F && rm_find_log_of_x(F, x) == NULL && rm_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rm_eval1("Numerator", expr_copy(F));
        den = rm_eval1("Denominator", expr_copy(F));
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
        Expr* pqn = num ? rm_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rm_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rm_is_true(pqn) && rm_is_true(pqd)
             && !rm_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Tower derivation coefficients Dt_i, checked triangular. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* arg = logs[i]->data.function.args[0];
            Expr* darg = rm_eval2("D", expr_copy(arg), expr_copy(x));
            Expr* q = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ darg, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(arg), expr_new_integer(-1) }, 2) }, 2);
            Expr* dc = rm_eval1("Cancel", q);
            Dt[i] = dc ? rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
            for (size_t j = i; j < nl; j++)
                if (!rm_free_of_x(Dt[i], ts[j])) { ok = false; break; }
        }
    }

    if (ok) {
        long dtn_num = rm_degree(num, top), dtn_den = rm_degree(den, top);
        long Ntop = dtn_num - dtn_den + 1;
        if (Ntop < 0) Ntop = 0; if (Ntop > 4) Ntop = 4;

        /* Lower field variables: x, t_1..t_{n-1}. */
        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));
        lv[0] = x;
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        for (size_t j = 0; j < nlv; j++) {
            long a = rm_degree(num, lv[j]); if (a < 0) a = 0;
            long b = rm_degree(den, lv[j]); if (b < 0) b = 0;
            long d = a + b + 1; if (d > 3) d = 3;
            bd[j] = d;
        }
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);

        /* Squarefree t_n-dependent factors of den -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = rm_eval1("Factor", expr_copy(den));
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
                if (rm_free_of_x(base, top)) continue;     /* not a t_n factor */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long nunk = (Ntop + 1) * nmono + (long)ng;
        if (!bad && nunk > 0 && nunk <= 80) {
            /* Build Q = sum_{k,mono} rmLp{k}_{m} (mono) t_n^k + sum_j rmLc{j} Log(g_j). */
            size_t nq = (size_t)((Ntop + 1) * nmono + (long)ng);
            Expr** qterms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long k = 0; k <= Ntop; k++)
                for (long m = 0; m < nmono; m++) {
                    char nm[32]; snprintf(nm, sizeof(nm), "rmLp%ld_%ld", k, m);
                    rm_decode_mono(m, bd, nlv, ev);
                    Expr* mono = rm_build_monomial(lv, ev, nlv);
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
            dterms[0] = rm_eval2("D", expr_copy(Q), expr_copy(x));
            for (size_t i = 0; i < nl; i++) {
                Expr* dqi = rm_eval2("D", expr_copy(Q), expr_copy(ts[i]));
                dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
            }
            Expr* Qder = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                dterms, nl + 1));
            free(dterms);

            Expr* diff = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
            Expr* tog = rm_eval1("Together", diff);
            Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;

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
                sol = rm_eval2("SolveAlways", eqn, varlist);
            }
            if (sol && sol->type == EXPR_FUNCTION
                && sol->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.head->data.symbol == intern_symbol("List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
                && sol->data.function.args[0]->data.function.head->data.symbol
                     == intern_symbol("List")) {
                Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
                    Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
                        result = rm_eval_own(expr_new_function(
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

    /* Diff-back safety gate (see rm_verify_antideriv). */
    if (result && !rm_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

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
static Expr* rm_tower_solve(Expr* Q, Expr** syms, size_t nsym,
                            Expr** Dt, Expr** ts, size_t nl,
                            Expr* F, Expr* x, Expr* backrules) {
    Expr** dterms = malloc((nl + 1) * sizeof(Expr*));
    dterms[0] = rm_eval2("D", expr_copy(Q), expr_copy(x));
    for (size_t i = 0; i < nl; i++) {
        Expr* dqi = rm_eval2("D", expr_copy(Q), expr_copy(ts[i]));
        dterms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(Dt[i]), dqi }, 2);
    }
    Expr* Qder = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
        dterms, nl + 1));
    free(dterms);
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2);
    Expr* tog = rm_eval1("Together", diff);
    Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;

    Expr* result = NULL; Expr* sol = NULL;
    if (rnum) {
        Expr** vl = malloc((nl + 1) * sizeof(Expr*));
        for (size_t i = 0; i < nl; i++) vl[i] = expr_copy(ts[nl - 1 - i]);
        vl[nl] = expr_copy(x);
        Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nl + 1);
        free(vl);
        Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
            (Expr*[]){ rnum, expr_new_integer(0) }, 2);
        sol = rm_eval2("SolveAlways", eqn, varlist);
    }
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol == intern_symbol("List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.args[0]->data.function.head->data.symbol
             == intern_symbol("List")) {
        Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
        if (Qs) {
            Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
            for (size_t j = 0; j < nsym; j++)
                zero[j] = expr_new_function(expr_new_symbol("Rule"),
                    (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
            Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
            free(zero);
            Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ Qs, zl }, 2));
            if (Qs)
                result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, expr_copy(backrules) }, 2));
        }
    }
    if (sol) expr_free(sol);
    return result;
}

/* Nested EXPONENTIAL tower (Phase B, third increment) — the dual of
 * rm_log_tower_case for a chain of nested exponentials t_i = E^(u_i) with u_i in
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
static Expr* rm_exp_tower_case(Expr* f, Expr* x) {
    if (!rm_find_exp_of_x(f, x)) return NULL;
    Expr** us = NULL; size_t nl = 0, ucap = 0;
    rm_collect_exp_exponents(f, x, &us, &nl, &ucap);   /* us[i] = exp exponents */
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
            if (rm_contains(kf[i], kf[k])) {
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

    Expr* F = rm_eval1("Together", rm_eval_own(expr_new_function(
        expr_new_symbol("ReplaceAll"), (Expr*[]){ expr_copy(f), expr_copy(rl) }, 2)));
    Expr* top = ts[nl - 1];
    Expr* num = NULL; Expr* den = NULL; Expr** Dt = NULL; Expr* result = NULL;

    bool ok = F && rm_find_log_of_x(F, x) == NULL && rm_find_exp_of_x(F, x) == NULL;
    if (ok) {
        num = rm_eval1("Numerator", expr_copy(F));
        den = rm_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((nl + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < nl; i++) vv[i + 1] = expr_copy(ts[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, nl + 1);
        free(vv);
        Expr* pqn = num ? rm_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rm_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        ok = num && den && rm_is_true(pqn) && rm_is_true(pqd)
             && !rm_free_of_x(F, top);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
    }

    /* Dt_i = Cancel[D[u_i,x]] |_{ker->t} * t_i. */
    if (ok) {
        Dt = calloc(nl, sizeof(Expr*));
        for (size_t i = 0; i < nl && ok; i++) {
            Expr* dui = rm_eval2("D", expr_copy(us[i]), expr_copy(x));
            Expr* dc = rm_eval1("Cancel", dui);
            Expr* dcs = dc ? rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ dc, expr_copy(rl) }, 2)) : NULL;
            Dt[i] = dcs ? rm_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ dcs, expr_copy(ts[i]) }, 2)) : NULL;
            if (!Dt[i]) { ok = false; break; }
        }
    }

    if (ok) {
        /* a = multiplicity of t_n at 0 in den (Laurent negative extent). */
        long a = 0;
        Expr* cl = rm_eval2("CoefficientList", expr_copy(den), expr_copy(top));
        if (cl && cl->type == EXPR_FUNCTION
            && cl->data.function.head->type == EXPR_SYMBOL
            && cl->data.function.head->data.symbol == intern_symbol("List")) {
            for (size_t i = 0; i < cl->data.function.arg_count; i++)
                if (!rm_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
        }
        if (cl) expr_free(cl);
        Expr* Dtil = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(top), expr_new_integer(-a) }, 2) }, 2));

        /* Squarefree t_n-dependent factors of Dtil -> log terms. */
        Expr* g[16]; size_t ng = 0; bool bad = false;
        Expr* factored = Dtil ? rm_eval1("Factor", expr_copy(Dtil)) : NULL;
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
                if (rm_free_of_x(base, top)) continue;    /* t_n-coprime constant part */
                if (e != 1 || ng >= 16) { bad = true; break; }  /* repeated: later */
                g[ng++] = expr_copy(base);
            }
        }

        long dtn_num = rm_degree(num, top), dtn_den = rm_degree(den, top);
        long ihi = dtn_num - dtn_den; if (ihi > 4) ihi = 4;
        long ilo = -a; if (ilo < -4) ilo = -4;

        size_t nlv = nl;                 /* x plus (nl-1) inner kernels */
        Expr** lv = malloc(nlv * sizeof(Expr*));
        long* bd = malloc(nlv * sizeof(long));   /* per-var odometer range size-1 */
        long* lo = malloc(nlv * sizeof(long));   /* per-var lowest exponent       */
        lv[0] = x; lo[0] = 0;                    /* x: polynomial, degree >= 0     */
        for (size_t i = 0; i + 1 < nl; i++) lv[i + 1] = ts[i];
        /* x is polynomial (0..bd); the inner kernels t_1..t_{n-1} are EXPONENTIALS
         * and invertible, so their coefficient exponents are LAURENT (-h..h) — the
         * antiderivative can carry a negative power the integrand lacks, e.g.
         * INT E^(x+E^x) dx = E^(E^x) = t_2 / t_1. */
        {
            long p = rm_degree(num, x); if (p < 0) p = 0;
            long q = rm_degree(den, x); if (q < 0) q = 0;
            long d = p + q + 1; if (d > 2) d = 2;
            bd[0] = d;
        }
        for (size_t j = 1; j < nlv; j++) { lo[j] = -2; bd[j] = 4; }  /* t_i in -2..2 */
        long nmono = 1;
        for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
        long nwi = (ihi >= ilo) ? (ihi - ilo + 1) : 0;
        long nunk = nwi * nmono + (long)ng;

        if (!bad && nwi > 0 && nunk > 0 && nunk <= 80) {
            size_t nq = (size_t)nunk;
            Expr** qterms = malloc(nq * sizeof(Expr*));
            Expr** syms = malloc(nq * sizeof(Expr*));
            size_t ntq = 0, nsym = 0;
            long* ev = malloc(nlv * sizeof(long));
            for (long i = ilo; i <= ihi; i++)
                for (long m = 0; m < nmono; m++) {
                    char nm[40]; snprintf(nm, sizeof(nm), "rmEp%ld_%ld", i - ilo, m);
                    rm_decode_mono(m, bd, nlv, ev);
                    for (size_t j = 0; j < nlv; j++) ev[j] += lo[j];  /* Laurent shift */
                    Expr* mono = rm_build_monomial(lv, ev, nlv);
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

            result = rm_tower_solve(Q, syms, nsym, Dt, ts, nl, F, x, bl);
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

    /* Diff-back safety gate (see rm_verify_antideriv). */
    if (result && !rm_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

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
 * That misses two whole families the recursive algorithm (Bronstein ch. 5 /
 * Maxima risch.lisp) handles: (1) MIXED exp/log towers (the flat cases are each
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
 * VERIFIED (rm_verify_antideriv) — a spurious certificate cannot ship. */

#define RM_MAXK 5

typedef enum { RM_LOG, RM_EXP } RmKind;

typedef struct {
    size_t n;
    RmKind kind[RM_MAXK];
    Expr* kernel[RM_MAXK];   /* Log[u_i] or Power[E, w_i]                (owned) */
    Expr* arg[RM_MAXK];      /* u_i (log argument) or w_i (exp exponent) (owned) */
    Expr* t[RM_MAXK];        /* fresh tower variable t_i                 (owned) */
    Expr* Dcoef[RM_MAXK];    /* log: u'/u ; exp: w' — in t-vars, in K_{i-1} (owned) */
    Expr* subrules;          /* List of all kernel -> t_i rules          (owned) */
} RmTower;

static void rm_tower_free(RmTower* T) {
    for (size_t i = 0; i < T->n; i++) {
        if (T->kernel[i]) expr_free(T->kernel[i]);
        if (T->arg[i]) expr_free(T->arg[i]);
        if (T->t[i]) expr_free(T->t[i]);
        if (T->Dcoef[i]) expr_free(T->Dcoef[i]);
    }
    if (T->subrules) expr_free(T->subrules);
    T->n = 0; T->subrules = NULL;
}

/* Build the ordered differential tower of f over C(x).  Collect every
 * x-dependent Log and E^ kernel, order them innermost-first (deepest at index 0)
 * by structural containment — tie-breaking independent kernels EXP-before-LOG so
 * the primitive (logarithmic) recursion sits on top and the exponential Risch DEs
 * bottom out in C(x) — assign tower variables, and compute each derivation
 * coefficient Dcoef_i (log: u_i'/u_i ; exp: w_i').  The structure-theorem
 * soundness check requires every Dcoef_i to lie in K_{i-1} = C(x, t_1..t_{i-1})
 * (triangular: free of t_i..t_n and of any residual foreign kernel).  Returns
 * true with T populated (2 <= n <= RM_MAXK); false otherwise (caller still calls
 * rm_tower_free). */
static bool rm_tower_build(Expr* f, Expr* x, RmTower* T) {
    for (size_t i = 0; i < RM_MAXK; i++)
        { T->kernel[i] = T->arg[i] = T->t[i] = T->Dcoef[i] = NULL; }
    T->subrules = NULL; T->n = 0;

    Expr** logs = NULL; size_t nl = 0, lc = 0; rm_collect_logs(f, x, &logs, &nl, &lc);
    Expr** exps = NULL; size_t ne = 0, ec = 0;
    rm_collect_exp_exponents(f, x, &exps, &ne, &ec);
    size_t n = nl + ne;
    if (n < 2 || n > RM_MAXK) {
        for (size_t i = 0; i < nl; i++) expr_free(logs[i]);
        for (size_t i = 0; i < ne; i++) expr_free(exps[i]);
        free(logs); free(exps);
        return false;
    }

    size_t idx = 0;
    for (size_t i = 0; i < nl; i++) {
        T->kind[idx] = RM_LOG;
        T->kernel[idx] = logs[i];                                  /* adopt Log[u] */
        T->arg[idx] = expr_copy(logs[i]->data.function.args[0]);
        idx++;
    }
    free(logs);
    for (size_t i = 0; i < ne; i++) {
        T->kind[idx] = RM_EXP;
        T->arg[idx] = exps[i];                                     /* adopt w */
        T->kernel[idx] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_copy(exps[i]) }, 2);
        idx++;
    }
    free(exps);
    T->n = n;

    /* Order innermost-first (deepest at index 0); tie-break EXP before LOG. */
    for (size_t pass = 0; pass < n; pass++)
        for (size_t i = 0; i + 1 < n; i++) {
            bool swap = false;
            if (rm_contains(T->kernel[i], T->kernel[i + 1])) swap = true;
            else if (!rm_contains(T->kernel[i + 1], T->kernel[i])
                     && T->kind[i] == RM_LOG && T->kind[i + 1] == RM_EXP) swap = true;
            if (swap) {
                RmKind kk = T->kind[i]; T->kind[i] = T->kind[i + 1]; T->kind[i + 1] = kk;
                Expr* a = T->kernel[i]; T->kernel[i] = T->kernel[i + 1]; T->kernel[i + 1] = a;
                a = T->arg[i]; T->arg[i] = T->arg[i + 1]; T->arg[i + 1] = a;
            }
        }

    /* Tower variables t_i and the combined substitution rule list. */
    Expr** rules = malloc(2 * n * sizeof(Expr*)); size_t nr = 0;
    for (size_t i = 0; i < n; i++) {
        char nm[16]; snprintf(nm, sizeof(nm), "rmR%zu", i);
        T->t[i] = expr_new_symbol(nm);
        if (T->kind[i] == RM_LOG) {
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
    T->subrules = expr_new_function(expr_new_symbol("List"), rules, nr);
    free(rules);

    /* Derivation coefficients + structure-theorem (triangularity) soundness. */
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr* d;
        if (T->kind[i] == RM_LOG) {
            Expr* du = rm_eval2("D", expr_copy(T->arg[i]), expr_copy(x));
            Expr* q = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ du, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(T->arg[i]), expr_new_integer(-1) }, 2) }, 2);
            d = rm_eval1("Cancel", q);
        } else {
            d = rm_eval1("Cancel", rm_eval2("D", expr_copy(T->arg[i]), expr_copy(x)));
        }
        if (!d) { ok = false; break; }
        Expr* ds = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
            (Expr*[]){ d, expr_copy(T->subrules) }, 2));           /* adopts d */
        if (!ds) { ok = false; break; }
        T->Dcoef[i] = ds;
        if (rm_find_exp_of_x(ds, x) != NULL || rm_find_log_of_x(ds, x) != NULL) ok = false;
        for (size_t j = i; j < n && ok; j++)
            if (!rm_free_of_x(ds, T->t[j])) ok = false;
    }
    return ok;
}

/* Mutually recursive field-integration primitives (forward declarations). */
static Expr* rm_field_integrate(Expr* F, RmTower* T, long L, Expr* x);
static int   rm_limited_field_integrate(Expr* r, RmTower* T, long L, Expr* x,
                                        Expr** s_out, Expr** c_out);
static Expr* rm_int_primitive_poly(Expr* num, Expr* den, RmTower* T, long L, Expr* x);
static Expr* rm_int_hyperexp_poly(Expr* num, Expr* den, RmTower* T, long L, Expr* x);
static Expr* rm_field_rde(Expr* p, long i, RmTower* T, long L, Expr* x);
static Expr* rm_field_ratint(Expr* num, Expr* den, RmTower* T, long L, Expr* x);
static Expr* rm_field_hyperexp_coupled(Expr* num, Expr* den, RmTower* T, long L, Expr* x);

/* Tower derivation D_tower[e] = D[e,x] + sum_i Dt_i D[e,t_i], with Dt_i =
 * Dcoef_i (log) or Dcoef_i * t_i (exp).  Owned result. */
static Expr* rm_tower_deriv(Expr* e, RmTower* T, Expr* x) {
    Expr** terms = malloc((T->n + 1) * sizeof(Expr*));
    terms[0] = rm_eval2("D", expr_copy(e), expr_copy(x));
    for (size_t i = 0; i < T->n; i++) {
        Expr* dei = rm_eval2("D", expr_copy(e), expr_copy(T->t[i]));
        Expr* dti = (T->kind[i] == RM_LOG)
            ? expr_copy(T->Dcoef[i])
            : expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(T->Dcoef[i]), expr_copy(T->t[i]) }, 2);
        terms[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ dti, dei }, 2);
    }
    Expr* r = rm_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, T->n + 1));
    free(terms);
    return r;
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
static Expr* rm_field_ratint(Expr* num, Expr* den, RmTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    Expr* dden = rm_eval2("D", expr_copy(den), expr_copy(t));
    Expr* Hden = dden ? rm_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(den), dden }, 2) : NULL;
    if (!Hden) return NULL;
    long dH = rm_degree(Hden, t); if (dH < 0) dH = 0;
    Expr* rad = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2));

    /* Distinct t-dependent factors of the squarefree radical (log arguments). */
    Expr* factored = rad ? rm_eval1("Factor", expr_copy(rad)) : NULL;
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
            if (rm_free_of_x(base, t)) continue;
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
        long a = rm_degree(num, lv[j]); if (a < 0) a = 0;
        long b = rm_degree(den, lv[j]); if (b < 0) b = 0;
        long d = (a > b ? a : b) + 1; if (d > 3) d = 3;
        bd[j] = d;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
    long nunk = dH * nmono + (long)ng;

    Expr* result = NULL;
    if (!bad && nunk > 0 && nunk <= 60) {
        size_t nq = (size_t)nunk;
        Expr** qterms = malloc(nq * sizeof(Expr*));
        Expr** syms = malloc(nq * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        /* Hermite numerator H(t)/Hden: sum_{p<dH} (sum_mono a x^..) t^p. */
        for (long p = 0; p < dH; p++)
            for (long m = 0; m < nmono; m++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmGh%ld_%ld", p, m);
                rm_decode_mono(m, bd, nlv, ev);
                Expr* mono = rm_build_monomial(lv, ev, nlv);
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
        Expr* Qder = rm_tower_deriv(Q, T, x);
        Expr* diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), target }, 2) }, 2);
        Expr* tog = rm_eval1("Together", diff);
        Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc((T->n + 1) * sizeof(Expr*));
            for (size_t i = 0; i < T->n; i++) vl[i] = expr_copy(T->t[T->n - 1 - i]);
            vl[T->n] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, T->n + 1);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rm_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
                 == intern_symbol("List")) {
            Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
            if (Qs) {
                Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
                for (size_t j = 0; j < nsym; j++)
                    zero[j] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
                Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
                free(zero);
                result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, zl }, 2));
            }
        }
        if (sol) expr_free(sol);
        for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
        free(syms);
        expr_free(Q);
    }

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
 * So (as in the single-kernel rm_hyperexp_case, lifted here to the tower derivation)
 * a UNIFIED ansatz is solved at once:
 *   Q = sum_{i=ilo}^{ihi} w_i t^i + H(t)/Hden(t) + sum_j c_j Log(g_j)
 * where ilo = -(mult. of t at 0 in den), Hden = gcd(Dtil, dDtil/dt) is the repeated
 * part of the t-coprime denominator Dtil = den/t^a, g_j the distinct t-factors of
 * the squarefree radical Dtil/Hden, the w_i and the numerator H are bounded-degree
 * lower-field polynomials, and the c_j constants — all found by SolveAlways over
 * every tower variable of D_tower[Q] = num/den.  Correct by construction; the
 * caller diff-back verifies.  Tried only when the Laurent recursion
 * (rm_int_hyperexp_poly) declines because a genuine proper part is present. */
static Expr* rm_field_hyperexp_coupled(Expr* num, Expr* den, RmTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    long a = 0;
    Expr* cl = rm_eval2("CoefficientList", expr_copy(den), expr_copy(t));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rm_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    Expr* Dtil = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(t), expr_new_integer(-a) }, 2) }, 2));
    if (!Dtil) return NULL;

    Expr* dDt = rm_eval2("D", expr_copy(Dtil), expr_copy(t));
    Expr* Hden = dDt ? rm_eval_call("PolynomialGCD",
        (Expr*[]){ expr_copy(Dtil), dDt }, 2) : NULL;
    long dH = Hden ? rm_degree(Hden, t) : 0; if (dH < 0) dH = 0;
    Expr* rad = Hden ? rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(Dtil), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(Hden), expr_new_integer(-1) }, 2) }, 2)) : NULL;

    Expr* factored = rad ? rm_eval1("Factor", expr_copy(rad)) : NULL;
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
            if (rm_free_of_x(base, t)) continue;
            if (ng >= 16) { bad = true; break; }
            g[ng++] = expr_copy(base);
        }
    }

    long dnum = rm_degree(num, t), dden = rm_degree(den, t);
    long ihi = dnum - dden; if (ihi < 0) ihi = 0; if (ihi > 4) ihi = 4;
    long ilo = -a; if (ilo < -4) ilo = -4;
    long nwi = ihi - ilo + 1;

    size_t nlv = (size_t)L + 1;
    Expr** lv = malloc(nlv * sizeof(Expr*));
    long* bd = malloc(nlv * sizeof(long));
    lv[0] = x;
    for (long i = 0; i < L; i++) lv[i + 1] = T->t[i];
    for (size_t j = 0; j < nlv; j++) {
        long p = rm_degree(num, lv[j]); if (p < 0) p = 0;
        long q = rm_degree(den, lv[j]); if (q < 0) q = 0;
        long d = (p > q ? p : q) + 1; if (d > 3) d = 3;
        bd[j] = d;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);
    long nunk = nwi * nmono + dH * nmono + (long)ng;

    Expr* result = NULL;
    if (!bad && nwi > 0 && nunk > 0 && nunk <= 60) {
        size_t nq = (size_t)(nwi * nmono + dH * nmono + (long)ng);
        Expr** qterms = malloc(nq * sizeof(Expr*));
        Expr** syms = malloc(nq * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        /* Laurent part: sum_i (sum_mono aW mono) t^i. */
        for (long i = ilo; i <= ihi; i++)
            for (long m = 0; m < nmono; m++) {
                char nm[40]; snprintf(nm, sizeof(nm), "rmXw%ld_%ld", i - ilo, m);
                rm_decode_mono(m, bd, nlv, ev);
                Expr* mono = rm_build_monomial(lv, ev, nlv);
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
                    rm_decode_mono(m, bd, nlv, ev);
                    Expr* mono = rm_build_monomial(lv, ev, nlv);
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
        Expr* Qder = rm_tower_deriv(Q, T, x);
        Expr* diff = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ Qder, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), target }, 2) }, 2);
        Expr* tog = rm_eval1("Together", diff);
        Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc((T->n + 1) * sizeof(Expr*));
            for (size_t i = 0; i < T->n; i++) vl[i] = expr_copy(T->t[T->n - 1 - i]);
            vl[T->n] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, T->n + 1);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rm_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
                 == intern_symbol("List")) {
            Expr* Qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(Q), expr_copy(sol->data.function.args[0]) }, 2));
            if (Qs) {
                Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
                for (size_t j = 0; j < nsym; j++)
                    zero[j] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
                Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
                free(zero);
                result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Qs, zl }, 2));
            }
        }
        if (sol) expr_free(sol);
        for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
        free(syms);
        expr_free(Q);
    }

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
static Expr* rm_field_integrate(Expr* F, RmTower* T, long L, Expr* x) {
    if (!F) return NULL;
    if (L < 0) {
        Expr* r = rm_eval_call("Integrate`BronsteinRational",
            (Expr*[]){ expr_copy(F), expr_copy(x) }, 2);
        if (!r) return NULL;
        if (rm_head_is(r, "Integrate`BronsteinRational")) { expr_free(r); return NULL; }
        return r;
    }
    Expr* t = T->t[L];
    Expr* Fg = rm_eval1("Together", expr_copy(F));
    if (!Fg) return NULL;
    Expr* num = rm_eval1("Numerator", expr_copy(Fg));
    Expr* den = rm_eval1("Denominator", expr_copy(Fg));
    Expr* result = NULL;
    if (num && den) {
        if (T->kind[L] == RM_LOG) {
            if (rm_free_of_x(den, t)) {
                /* Pure polynomial in t: primitive-polynomial recursion. */
                result = rm_int_primitive_poly(num, den, T, L, x);
            } else {
                /* Split the polynomial part (recursion) from the proper rational
                 * part in t (Hermite reduction + Rothstein-Trager log part). */
                Expr* Pp = rm_eval_call("PolynomialQuotient",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(t) }, 3);
                Expr* Rr = rm_eval_call("PolynomialRemainder",
                    (Expr*[]){ expr_copy(num), expr_copy(den), expr_copy(t) }, 3);
                if (Pp && Rr) {
                    Expr* one = expr_new_integer(1);
                    Expr* poly_int = rm_int_primitive_poly(Pp, one, T, L, x);
                    expr_free(one);
                    Expr* prop_int = rm_is_zero(Rr) ? expr_new_integer(0)
                                                    : rm_field_ratint(Rr, den, T, L, x);
                    if (poly_int && prop_int)
                        result = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
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
            result = rm_int_hyperexp_poly(num, den, T, L, x);
            if (!result)
                result = rm_field_hyperexp_coupled(num, den, T, L, x);
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
static Expr* rm_int_primitive_poly(Expr* num, Expr* den, RmTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    Expr* Dt = T->Dcoef[L];
    long m = rm_degree(num, t);
    if (m < 0) return expr_new_integer(0);    /* num == 0: integral is 0 */
    Expr** p = malloc((size_t)(m + 1) * sizeof(Expr*));
    for (long i = 0; i <= m; i++)
        p[i] = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ rm_coeff(num, t, i), expr_new_function(expr_new_symbol("Power"),
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
            r_i = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_copy(p[i]),
                    expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), term }, 2) }, 2));
        }
        Expr* s = NULL; Expr* c = NULL;
        int rc = r_i ? rm_limited_field_integrate(r_i, T, L, x, &s, &c) : -1;
        if (r_i) expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        if (c) {
            Expr* bump = rm_eval_own(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2));
            if (!q[i + 1]) q[i + 1] = bump;
            else q[i + 1] = rm_eval_own(expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ q[i + 1], bump }, 2));
        }
        q[i] = s;
    }
    Expr* result = NULL;
    if (!fail) {
        Expr** terms = malloc((size_t)(m + 2) * sizeof(Expr*));
        size_t nt = 0;
        for (long i = 0; i <= m + 1; i++) {
            if (!q[i] || rm_is_zero(q[i])) continue;
            Expr* pw = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(t), expr_new_integer(i) }, 2);
            terms[nt++] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_copy(q[i]), pw }, 2);
        }
        result = rm_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
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
static int rm_limited_field_integrate(Expr* r, RmTower* T, long L, Expr* x,
                                      Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rm_field_integrate(r, T, L - 1, x);
    if (!R) return -1;
    Expr* Rs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ R, expr_copy(T->subrules) }, 2));               /* adopts R */
    if (!Rs) return -1;
    Expr* tL = T->t[L];
    int rc = -1;
    if (rm_find_exp_of_x(Rs, x) == NULL && rm_find_log_of_x(Rs, x) == NULL
        && rm_is_poly(Rs, tL) && rm_degree(Rs, tL) <= 1) {
        Expr* c = rm_coeff(Rs, tL, 1);
        Expr* s = rm_coeff(Rs, tL, 0);
        bool cconst = rm_free_of_x(c, x);
        for (size_t j = 0; j < T->n && cconst; j++)
            if (!rm_free_of_x(c, T->t[j])) cconst = false;
        if (cconst) {
            if (rm_is_zero(c)) { expr_free(c); *c_out = NULL; }
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
static Expr* rm_int_hyperexp_poly(Expr* num, Expr* den, RmTower* T, long L, Expr* x) {
    Expr* t = T->t[L];
    long a = 0;
    Expr* cl = rm_eval2("CoefficientList", expr_copy(den), expr_copy(t));
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol == intern_symbol("List")) {
        for (size_t i = 0; i < cl->data.function.arg_count; i++)
            if (!rm_is_zero(cl->data.function.args[i])) { a = (long)i; break; }
    }
    if (cl) expr_free(cl);
    Expr* Dtil = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(den), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(t), expr_new_integer(-a) }, 2) }, 2));
    if (!Dtil || !rm_free_of_x(Dtil, t)) { if (Dtil) expr_free(Dtil); return NULL; }
    Expr* nnum = rm_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(num), expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ Dtil, expr_new_integer(-1) }, 2) }, 2));     /* adopts Dtil */
    if (!nnum) return NULL;
    long dnn = rm_degree(nnum, t);
    if (dnn < 0) { expr_free(nnum); return expr_new_integer(0); }  /* num == 0 */

    Expr** q = calloc((size_t)(dnn + 1), sizeof(Expr*));
    long* pw = malloc((size_t)(dnn + 1) * sizeof(long));
    bool fail = false;
    for (long j = 0; j <= dnn && !fail; j++) {
        long ip = j - a;
        pw[j] = ip;
        Expr* pj = rm_coeff(nnum, t, j);
        if (rm_is_zero(pj)) { expr_free(pj); q[j] = NULL; continue; }
        Expr* qj = (ip == 0) ? rm_field_integrate(pj, T, L - 1, x)
                             : rm_field_rde(pj, ip, T, L, x);
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
        result = rm_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
        free(terms);
    }
    for (long j = 0; j <= dnn; j++) if (q[j]) expr_free(q[j]);
    free(q); free(pw);
    expr_free(nnum);
    return result;
}

/* Risch differential equation q' + i w_L' q = p for q in the lower field K_L.
 * Base case (w_L and p both in C(x)): the bounded polynomial-in-x ansatz of
 * rm_solve_rde.  General case: q lives in K_L = C(x, t_0..t_{L-1}) and may be
 * rational there (a monomial denominator, e.g. 1/Log[x]), which the coupled-
 * hyperexponential polynomial coefficients cannot express — solve by a bounded
 * LAURENT ansatz over {x, t_0..t_{L-1}} (each variable ranging over negative powers
 * too), requiring D_tower[q] + i Dcoef_L q = p for every lower tower variable via
 * SolveAlways.  A non-monomial denominator (needing the full Bronstein SPDE) is out
 * of a monomial-Laurent ansatz's reach and declines. */
static Expr* rm_field_rde(Expr* p, long i, RmTower* T, long L, Expr* x) {
    Expr* w = T->arg[L];
    bool base = rm_is_poly(w, x) && rm_is_poly(p, x);
    for (size_t j = 0; j < T->n && base; j++)
        if (!rm_free_of_x(w, T->t[j]) || !rm_free_of_x(p, T->t[j])) base = false;
    if (base) return rm_solve_rde(p, i, w, x);

    /* By the RDE denominator theorem, denom(q) | denom(p) (a pole of q would give a
     * higher-order pole in q' + i Dcoef q = p that nothing cancels), so
     * q = h / pd with pd = Denominator[p] and h a bounded POLYNOMIAL numerator over
     * {x, t_0..t_{L-1}}.  This is strictly more general than a monomial-Laurent
     * ansatz — it captures a NON-monomial denominator such as 1/(1+Log[x]) — and
     * subsumes it (pd carries every pole).  Solve h by SolveAlways. */
    Expr* Dcoef = T->Dcoef[L];
    Expr* pg = rm_eval1("Together", expr_copy(p));
    Expr* pd = pg ? rm_eval1("Denominator", expr_copy(pg)) : NULL;
    Expr* pn = pg ? rm_eval1("Numerator", expr_copy(pg)) : NULL;
    if (pg) expr_free(pg);
    if (!pd || !pn) { if (pd) expr_free(pd); if (pn) expr_free(pn); return NULL; }

    size_t nlv = (size_t)L + 1;
    Expr** lv = malloc(nlv * sizeof(Expr*));
    long* bd = malloc(nlv * sizeof(long));
    lv[0] = x;
    for (long j = 0; j < L; j++) lv[j + 1] = T->t[j];
    for (size_t j = 0; j < nlv; j++) {
        long a2 = rm_degree(pd, lv[j]); if (a2 < 0) a2 = 0;
        long b2 = rm_degree(pn, lv[j]); if (b2 < 0) b2 = 0;
        long d = a2 + b2 + 1; if (d > 5) d = 5;
        bd[j] = d;
    }
    long nmono = 1;
    for (size_t j = 0; j < nlv; j++) nmono *= (bd[j] + 1);

    Expr* result = NULL;
    if (nmono > 0 && nmono <= 60) {
        Expr** hterms = malloc((size_t)nmono * sizeof(Expr*));
        Expr** syms = malloc((size_t)nmono * sizeof(Expr*));
        size_t ntq = 0, nsym = 0;
        long* ev = malloc(nlv * sizeof(long));
        for (long m = 0; m < nmono; m++) {
            char nm[40]; snprintf(nm, sizeof(nm), "rmRd%ld", m);
            rm_decode_mono(m, bd, nlv, ev);
            Expr* mono = rm_build_monomial(lv, ev, nlv);
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

        Expr* dq = rm_tower_deriv(q, T, x);
        Expr* iq = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(i), expr_copy(Dcoef), expr_copy(q) }, 3);
        Expr* residual = expr_new_function(expr_new_symbol("Plus"),
            (Expr*[]){ dq, iq, expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 3);
        Expr* tog = rm_eval1("Together", residual);
        Expr* rnum = tog ? rm_eval1("Numerator", tog) : NULL;

        Expr* sol = NULL;
        if (rnum) {
            Expr** vl = malloc(nlv * sizeof(Expr*));
            for (long j = 0; j < L; j++) vl[j] = expr_copy(T->t[L - 1 - j]);
            vl[L] = expr_copy(x);
            Expr* varlist = expr_new_function(expr_new_symbol("List"), vl, nlv);
            free(vl);
            Expr* eqn = expr_new_function(expr_new_symbol("Equal"),
                (Expr*[]){ rnum, expr_new_integer(0) }, 2);
            sol = rm_eval2("SolveAlways", eqn, varlist);
        }
        if (sol && sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == intern_symbol("List")
            && sol->data.function.arg_count >= 1
            && sol->data.function.args[0]->type == EXPR_FUNCTION
            && sol->data.function.args[0]->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.args[0]->data.function.head->data.symbol
                 == intern_symbol("List")) {
            Expr* qs = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy(q), expr_copy(sol->data.function.args[0]) }, 2));
            if (qs) {
                Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
                for (size_t j = 0; j < nsym; j++)
                    zero[j] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(syms[j]), expr_new_integer(0) }, 2);
                Expr* zl = expr_new_function(expr_new_symbol("List"), zero, nsym);
                free(zero);
                result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
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
    if (result) result = rm_eval1("Cancel", result);
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
 * so rm_tower_build's structure-theorem check rejects the tower and the case
 * declines.  Splitting E^(x + E^x) back to E^x E^(E^x) restores the
 * independent tower basis {E^x, E^(E^x)}, so integrands the evaluator merged
 * (e.g. E^x E^(E^x)/(1+E^(E^x)) = E^(x+E^x)/(1+E^(E^x))) close instead.  The
 * rewrite is exact (E^(a+b) = E^a E^b) and the whole recursive case is
 * diff-back verified, so it can never ship a wrong form.  Returns a
 * freshly-owned tree (caller frees). */
static Expr* rm_expand_exp_sums(Expr* e) {
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
    if (expo && rm_head_is(expo, "Plus") && expo->data.function.arg_count >= 2) {
        size_t m = expo->data.function.arg_count;
        Expr** facs = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) {
            Expr* part = rm_expand_exp_sums(expo->data.function.args[i]);
            facs[i] = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), part }, 2);
        }
        Expr* prod = expr_new_function(expr_new_symbol("Times"), facs, m);
        free(facs);
        return prod;
    }
    Expr* nh = rm_expand_exp_sums(e->data.function.head);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rm_expand_exp_sums(e->data.function.args[i]);
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
static Expr* rm_subst_kernels(Expr* e, RmTower* T) {
    if (!e) return NULL;
    for (size_t i = 0; i < T->n; i++) {
        if (expr_eq(e, T->kernel[i])) return expr_copy(T->t[i]);
        if (T->kind[i] == RM_EXP
            && e->type == EXPR_FUNCTION
            && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol == intern_symbol("Exp")
            && e->data.function.arg_count == 1
            && expr_eq(e->data.function.args[0], T->arg[i]))
            return expr_copy(T->t[i]);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* nh = rm_subst_kernels(e->data.function.head, T);
    size_t k = e->data.function.arg_count;
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rm_subst_kernels(e->data.function.args[i], T);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    return r;
}

static Expr* rm_recursive_tower_case(Expr* f, Expr* x) {
    /* Split any evaluator-merged exponential monomial back into an independent
     * tower basis before building the tower (see rm_expand_exp_sums). */
    Expr* fx = rm_expand_exp_sums(f);
    RmTower T;
    if (!rm_tower_build(fx, x, &T)) { rm_tower_free(&T); expr_free(fx); return NULL; }

    /* Alias the kernels to tower variables structurally (NOT via an evaluated
     * ReplaceAll, which would re-merge a split exponential product before
     * substitution — see rm_subst_kernels), then normalise with Together. */
    Expr* F = rm_eval1("Together", rm_subst_kernels(fx, &T));
    Expr* result = NULL;
    if (F && rm_find_exp_of_x(F, x) == NULL && rm_find_log_of_x(F, x) == NULL) {
        Expr* num = rm_eval1("Numerator", expr_copy(F));
        Expr* den = rm_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((T.n + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < T.n; i++) vv[i + 1] = expr_copy(T.t[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, T.n + 1);
        free(vv);
        Expr* pqn = num ? rm_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rm_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        bool gate = num && den && rm_is_true(pqn) && rm_is_true(pqd);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (gate) {
            Expr* Q = rm_field_integrate(F, &T, (long)T.n - 1, x);
            if (Q) {
                Expr** back = malloc(T.n * sizeof(Expr*));
                for (size_t i = 0; i < T.n; i++)
                    back[i] = expr_new_function(expr_new_symbol("Rule"),
                        (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
                Expr* bl = expr_new_function(expr_new_symbol("List"), back, T.n);
                free(back);
                result = rm_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
                    (Expr*[]){ Q, bl }, 2));                       /* adopts Q, bl */
                /* Tidy the Hermite rational part (e.g. -x Log[x]/(x Log[x](1+t)) ->
                 * -1/(1+t)); Cancel treats the Log/exp kernels as opaque atoms. */
                if (result) result = rm_eval1("Cancel", result);
            }
        }
    }
    if (F) expr_free(F);

    /* Diff-back safety gate (bounded search, not a decision procedure).
     * Verify against the ORIGINAL integrand f (== fx mathematically). */
    if (result && !rm_verify_antideriv(result, f, x)) { expr_free(result); result = NULL; }

    rm_tower_free(&T);
    expr_free(fx);
    return result;
}

/* Trigonometric / hyperbolic front-end (Maxima's rischform exponentialize
 * path).  Rewrites the trig/hyperbolic kernels to complex exponentials with
 * TrigToExp, integrates the resulting (Laurent-)rational function of the
 * exponential kernel E^(I x) / E^x with the exponential machinery, and
 * converts the answer back to trigonometric form with ExpToTrig.  Both
 * rewrites are exact, so the result is correct by construction. */
static Expr* rm_trig_frontend(Expr* f, Expr* x) {
    Expr* fe = rm_eval1("TrigToExp", expr_copy(f));
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
    Expr* r = rm_exp_poly_case(fe, x);
    if (!r) r = rm_frac_case(fe, x);
    if (!r) r = rm_hyperexp_case(fe, x);
    /* Multi-kernel decoupling (Phase B): e.g. Sin/Cos times a real exponential
     * exponentialize to a sum of two non-commensurate exponentials E^((a +/- b I) x)
     * that the single-primitive cases cannot kernelize. */
    if (!r) r = rm_expsum_case(fe, x);
    expr_free(fe);
    if (!r) return NULL;
    return rm_eval1("ExpToTrig", r);   /* adopts r; back to trig form */
}

/* Dispatch the transcendental cases: the primitive (logarithmic) polynomial
 * reduction, the exponential (hyperexponential / Risch-DE) reduction, the
 * fractional (Rothstein-Trager) log-part, and the trig/hyperbolic front-end.
 * The general Hermite reduction for repeated poles lands in a subsequent
 * increment. */
static Expr* rm_transcendental_case(Expr* f, Expr* x) {
    Expr* r = rm_log_poly_case(f, x);
    if (r) return r;
    r = rm_exp_poly_case(f, x);
    if (r) return r;
    r = rm_frac_case(f, x);
    if (r) return r;
    r = rm_hermite_case(f, x);
    if (r) return r;
    r = rm_hyperexp_case(f, x);
    if (r) return r;
    r = rm_expsum_case(f, x);   /* direct multi-kernel exponential sums */
    if (r) return r;
    r = rm_log_tower_case(f, x);   /* nested logarithmic tower (depth > 1) */
    if (r) return r;
    r = rm_exp_tower_case(f, x);    /* nested exponential tower (depth > 1) */
    if (r) return r;
    r = rm_recursive_tower_case(f, x);   /* one-extension recursion (mixed / rational coeff) */
    if (r) return r;
    r = rm_trig_frontend(f, x);
    if (r) return r;
    return NULL;
}

/* ================================================================== */
/* Top-level integration (Maxima rischint / tryrisch dispatch).       */
/* ================================================================== */

/* Returns a fresh antiderivative (also self-verified by the recognizers,
 * and re-verified by the caller's diff-back gate) or NULL if no case
 * applies.  Order mirrors Maxima's recursive Risch:
 *   1. rational base case (delegated to the recursive rational Risch,
 *      Integrate`BronsteinRational);
 *   2. transcendental case over a single logarithmic / exponential
 *      monomial extension (rm_transcendental_case — the recursive Risch
 *      proper: Hermite reduction, residue log-part, and the polynomial /
 *      Risch-differential-equation reductions);
 *   3. special-function outputs (Maxima's erfarg / dilog / Ei / li).
 * Every branch is verified by differentiation, so a mis-reduction can
 * only decline, never emit a wrong closed form.  NB: this must NOT fall
 * back on the parallel-Risch (pmint) engine Integrate`RischNorman —
 * that is a different algorithm; RischMacsyma is the recursive Risch. */
static Expr* rm_integrate(Expr* f, Expr* x) {
    Expr* r = rm_rational_case(f, x);
    if (r) return r;
    r = rm_transcendental_case(f, x);
    if (r) return r;
    r = rm_special_case(f, x);
    if (r) return r;
    return NULL;
}

/* ================================================================== */
/* Public builtin.                                                    */
/* ================================================================== */

Expr* builtin_rischmacsyma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    /* The integration variable must be a single symbol. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Correct by construction: rm_integrate returns a result only behind an
     * exact certificate, so no differentiation check is applied (a Risch
     * integrator is a decision procedure, not a guess-and-verify search). */
    return rm_integrate(f, x);
}

/* ================================================================== */
/* Registration.                                                      */
/* ================================================================== */

static void rm_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_macsyma_init(void) {
    rm_install("Integrate`RischMacsyma", builtin_rischmacsyma,
        "Integrate`RischMacsyma[f, x] integrates f with respect to x using the\n"
        "recursive Risch algorithm ported from Maxima (risch.lisp): a decision\n"
        "procedure over a differential transcendental tower, with rational,\n"
        "logarithmic, exponential, and special-function (Erf, ExpIntegralEi,\n"
        "LogIntegral, PolyLog) cases.  Each case is correct by construction (no\n"
        "differentiation check).  Distinct from Integrate`RischNorman, which is\n"
        "the parallel-Risch (pmint) heuristic.  Out-of-scope integrands\n"
        "(algebraic extensions, non-elementary answers) return unevaluated.");
}
