/* risch_special.c — special-function outputs & base cases for the Risch integrator.
 *
 * Pure-rational delegation (Integrate`BronsteinRational), the Erf/Erfi,
 * ExpIntegralEi, LogIntegral and dilogarithm recognisers, and the
 * logarithmic-polynomial case.  Tried ahead of the tower machinery.  See
 * risch_special.h.
 */

#include "risch_special.h"
#include "risch_util.h"
#include "cherry_ei.h"
#include "cherry_li.h"
#include "cherry_dilog.h"

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

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================== */
/* Case: rational function of x  (rational Risch base case).           */
/* ================================================================== */

/* Delegate the pure-rational case to the Bronstein rational integrator
 * (the rational Risch base case: Hermite + Lazard-Rioboo-Trager).  Returns
 * NULL when the integrand is not rational in x (BronsteinRational leaves
 * its call unevaluated). */
Expr* rt_rational_case(Expr* f, Expr* x) {
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
    { "Erf",           rt_try_erf     }, /* K E^(a x^2 + b x + c)                */
    { "ExpIntegralEi", rt_try_ei      }, /* M E^(a x + b) / (c x + d)   [fast path] */
    { "ExpIntegralEi", rt_cherry_ei   }, /* g E^f, g,f in C(x): ei + erf (Cherry 1989) */
    { "LogIntegral",   rt_try_li      }, /* c w^(p-1) w' / Log[w]      [fast path] */
    { "LogIntegral",   rt_cherry_li   }, /* multi-li over C(x,Log[w])  (Cherry 1986) */
    { "PolyLog",       rt_try_dilog   }, /* K Log[1 + p x] / x        [fast path] */
    { "PolyLog",       rt_cherry_dilog }, /* R Log[w] -> LogLog + PolyLog[2] (Cherry) */
};

/* Try each registered special-function form in turn (order preserved). */
Expr* rt_special_case(Expr* f, Expr* x) {
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
 * *c_out (both owned) and returns 0; otherwise returns -1.
 *
 * `is_bottom` selects the theta^0 (constant) coefficient of the antiderivative,
 * which — unlike the higher theta^i coefficients, which must be rational so the
 * result stays a polynomial in theta over K — may be ANY elementary function of
 * x.  Bronstein's IntegratePrimitivePolynomial (§5.8) integrates the lowest term
 * by a full LimitedIntegrate in K, so a genuine new logarithm or an ArcTan from
 * an irreducible-quadratic log argument (e.g. the Sqrt[3] ArcTan of
 * Integrate[(x^5-1) Log[x^2-x+1], x]) is a legitimate part of the answer, not a
 * decline.  At the bottom level we therefore accept the complete lower-field
 * antiderivative wholesale (s = R, c = 0).  The caller diff-back verifies. */
static int rt_limited_integrate(Expr* r, Expr* x, Expr* u, int is_bottom,
                                Expr** s_out, Expr** c_out) {
    *s_out = NULL; *c_out = NULL;
    Expr* R = rt_eval_call("Integrate`BronsteinRational",
        (Expr*[]){ expr_copy(r), expr_copy(x) }, 2);
    if (!R) return -1;
    if (rt_head_is(R, "Integrate`BronsteinRational")) { expr_free(R); return -1; }
    if (is_bottom) { *s_out = R; *c_out = NULL; return 0; }
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
Expr* rt_log_poly_case(Expr* f, Expr* x) {
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
        int rc = rt_limited_integrate(r_i, x, u, i == 0, &s, &c);
        expr_free(r_i);
        if (rc != 0) { fail = true; break; }
        /* Fold the theta-term back: q[i+1] += c/(i+1).  The bottom level returns
         * c == NULL (the whole antiderivative is kept in s = q[0]); nothing to fold. */
        if (c) {
            Expr* bump = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ c, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_integer(i + 1), expr_new_integer(-1) }, 2) }, 2);
            if (q[i + 1] == NULL) {
                q[i + 1] = rt_eval_own(bump);
            } else {
                q[i + 1] = rt_eval_own(expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ q[i + 1], bump }, 2));   /* adopts old q[i+1] */
            }
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
    /* The recursion is correct by construction, but the bottom level now admits an
     * arbitrary elementary lower-field antiderivative (new logs / ArcTan); diff-back
     * verify so any surprise from the base integrator can only decline, never ship
     * a wrong form. */
    if (result && !rt_verify_antideriv(result, f, x)) {
        expr_free(result); result = NULL;
    }
    return result;
}
