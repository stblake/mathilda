/*
 * functionexpand.c -- FunctionExpand for Mathilda.
 *
 *   FunctionExpand[expr]
 *   FunctionExpand[expr, assum]
 *
 * expands special functions in expr into simpler ones, bottom-up. The rules:
 *
 *   MeijerG   -> named functions (meijer_g.c, meijerg_function_expand): the
 *                BesselK / HypergeometricU shapes in z or 1/z, else the named
 *                Slater expansion, else for p == q, m + n <= p a Piecewise of
 *                the inner and outer expansions. A log case is left alone.
 *
 *   PolyLog[n, w], integer n >= 2, w = E^(I Pi r) an exact root of unity
 *                (r rational, written with I, -1, (-1)^r, E^(I Pi r) and
 *                products of them):  Li_n(w) = R + I S, with x = r/2 mod 1,
 *       - the "Bernoulli part" (R for even n, S for odd n) always closed:
 *             Li_n(e^{2 pi i x}) + (-1)^n Li_n(e^{-2 pi i x})
 *                 = -(2 pi i)^n B_n(x) / n!            (0 <= x <= 1);
 *       - the "Clausen part" (R for odd n, S for even n) as a multiple of
 *         Zeta[n] when the cosine values are rational (q = den x in
 *         {1, 2, 3, 4, 6}: Ramanujan sums,  sum_k c_q(k)/k^n =
 *         Zeta[n] sum_{d|q} mu(q/d) d^{1-n}, shared by the phi(q) <= 2
 *         primitive p); S = +-Catalan for n = 2, q = 4; else by the Hurwitz
 *         decomposition  q^{-n} sum_{j=1..q} f(2 pi j p/q) zeta(n, j/q)  with
 *         zeta(n, a) = (-1)^n PolyGamma[n-1, a]/(n-1)!.
 *     Conjugate pairs therefore cancel their imaginary parts exactly:
 *     PolyLog[3, I] + PolyLog[3, -I] -> -3 Zeta[3]/16.
 *
 *   PolyLog[2, 2] = Pi^2/4 - I Pi Log[2] and PolyLog[3, 2] = 7/8 Zeta[3] +
 *                Pi^2/4 Log[2] - I Pi Log[2]^2/2 (inversion formula; PolyLog's
 *                own branch, Im < 0 on the cut).
 *
 *   PolyLog[2, z], z an exact number one anharmonic map (1 - z, 1/z,
 *                z/(z - 1), ...) away from a root of unity, 1/2 or 2:
 *                Li_2((1 + I)/2) = 5 Pi^2/96 - Log[2]^2/8 + I (Catalan -
 *                Pi Log[2]/8), Li_2(1 + I), Li_2(1/2 - I/2), ...
 *
 *   Log[w] -> I Pi r and Log[1 + w] -> Log[2 Cos[Pi r/2]] + I Pi r/2
 *                (-1 < r < 1) for the same roots of unity w.
 *
 * Everything else is left as it is. The assumptions argument is accepted for
 * compatibility and not yet used.
 *
 * fe_reduce_definite_result (functionexpand.h) applies the same machinery as
 * the final pass of a definite Integrate.
 *
 * Memory contract: the builtin takes ownership of res but must not free it.
 */

#include "functionexpand.h"
#include "polygamma.h"   /* polygamma_gauss_digamma */
#include "meijer_g.h"
#include "hrep.h"           /* hrep_pfq_expr: Roach's Formula algorithm */
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "arithmetic.h"         /* is_rational */
#include "sym_names.h"
#include "numeric.h"            /* numericalize_machine_unplanned */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Small helpers                                                      */
/* ------------------------------------------------------------------ */

static bool fe_is_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

/* `fn`, the source of a pure function, applied to the (copied) arguments and
 * evaluated. Owned result, or NULL if the template does not parse. */
static Expr* fe_apply(const char* fn, Expr* const* args, size_t n) {
    Expr* f = parse_expression(fn);
    if (!f) return NULL;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    Expr* call = expr_new_function(f, a, n);
    free(a);
    return eval_and_free(call);
}

static int64_t fe_gcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* Bound on the denominators carried; far beyond any root of unity a closed
 * form is wanted for, and it keeps every product below int64 range. */
#define FE_MAX_DEN 100000

/* r += n/d, reduced; false on a denominator past FE_MAX_DEN. */
static bool fe_rat_add(int64_t* rn, int64_t* rd, int64_t n, int64_t d) {
    if (d <= 0 || d > FE_MAX_DEN || *rd > FE_MAX_DEN) return false;
    int64_t g = fe_gcd(*rd, d);
    int64_t den = (*rd / g) * d;
    if (den > FE_MAX_DEN) return false;
    int64_t num = *rn * (d / g) + n * (*rd / g);
    int64_t h = fe_gcd(num, den);
    if (h == 0) h = 1;
    *rn = num / h; *rd = den / h;
    return true;
}

/* An exact rational leaf n/d (Integer or Rational). */
static bool fe_rational(const Expr* e, int64_t* n, int64_t* d) {
    return is_rational(e, n, d) && *d > 0;
}

/* Accumulate the angle (in units of Pi) of one factor of a root of unity. */
static bool fe_root_factor(const Expr* f, int64_t* rn, int64_t* rd) {
    int64_t n, d;
    if (f->type == EXPR_INTEGER) {
        if (f->data.integer == 1) return true;
        if (f->data.integer == -1) return fe_rat_add(rn, rd, 1, 1);
        return false;
    }
    if (fe_is_head(f, SYM_Complex) && f->data.function.arg_count == 2) {
        const Expr* re = f->data.function.args[0];
        const Expr* im = f->data.function.args[1];
        if (re->type != EXPR_INTEGER || re->data.integer != 0 ||
            im->type != EXPR_INTEGER) return false;
        if (im->data.integer == 1)  return fe_rat_add(rn, rd, 1, 2);
        if (im->data.integer == -1) return fe_rat_add(rn, rd, -1, 2);
        return false;
    }
    if (fe_is_head(f, SYM_Power) && f->data.function.arg_count == 2) {
        const Expr* b = f->data.function.args[0];
        const Expr* x = f->data.function.args[1];
        /* (-1)^(n/d) */
        if (b->type == EXPR_INTEGER && b->data.integer == -1 && fe_rational(x, &n, &d))
            return fe_rat_add(rn, rd, n, d);
        /* E^(Complex[0, n/d] Pi) */
        if (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_E &&
            fe_is_head(x, SYM_Times) && x->data.function.arg_count == 2) {
            const Expr* c = x->data.function.args[0];
            const Expr* pi = x->data.function.args[1];
            if (pi->type == EXPR_SYMBOL && pi->data.symbol.name == SYM_Pi &&
                fe_is_head(c, SYM_Complex) && c->data.function.arg_count == 2 &&
                c->data.function.args[0]->type == EXPR_INTEGER &&
                c->data.function.args[0]->data.integer == 0 &&
                fe_rational(c->data.function.args[1], &n, &d))
                return fe_rat_add(rn, rd, n, d);
        }
        return false;
    }
    return false;
}

/* w an exact root of unity E^(I Pi r): r = *rn / *rd reduced into (-1, 1]. */
static bool fe_root_of_unity(const Expr* w, int64_t* rn, int64_t* rd) {
    *rn = 0; *rd = 1;
    if (fe_is_head(w, SYM_Times)) {
        for (size_t i = 0; i < w->data.function.arg_count; i++)
            if (!fe_root_factor(w->data.function.args[i], rn, rd)) return false;
    } else if (!fe_root_factor(w, rn, rd)) {
        return false;
    }
    /* Reduce modulo 2 into (-1, 1]. */
    int64_t two = 2 * *rd;
    *rn %= two;
    if (*rn <= -*rd) *rn += two;
    if (*rn > *rd) *rn -= two;
    return true;
}

/* Rule families of the walk: MeijerG; PolyLog/Log at roots of unity (with
 * the PolyGamma closed forms); the PolyGamma closed forms alone. */
enum { FE_MEIJER = 1, FE_ROOTS = 2, FE_PGAMMA = 4, FE_PFQ = 8, FE_ALL = 15 };
static Expr* fe_walk(const Expr* e, bool* changed, int which);
static Expr* fe_expand(const Expr* e);
static Expr* fe_log_reciprocal(const Expr* e, bool* changed);

/* ------------------------------------------------------------------ */
/*  PolyGamma at 1/2, 1/4, 3/4                                         */
/* ------------------------------------------------------------------ */

/* PolyGamma[m, x] = (-1)^(m+1) m! zeta(m+1, x), integer m >= 1, with the
 * Hurwitz zeta closed where it is: n = m + 1,
 *     zeta(n, 1/2)       = (2^n - 1) Zeta[n],
 *     zeta(n, 1/4 | 3/4) = ((4^n - 2^n) Zeta[n] +- 4^n beta(n)) / 2,
 * the Dirichlet beta closed for n = 2 (Catalan) and odd n = 2k+1
 * (beta = (-1)^k EulerE[2k] Pi^n / (4^(k+1) (2k)!)); even n > 2 stays.
 * FunctionExpand[PolyGamma[1, 1/4]] = Pi^2 + 8 Catalan. */
static Expr* fe_polygamma_closed(const Expr* e) {
    if (!fe_is_head(e, SYM_PolyGamma) || e->data.function.arg_count != 2) return NULL;
    const Expr* me = e->data.function.args[0];
    if (me->type != EXPR_INTEGER || me->data.integer < 0 || me->data.integer > 200) return NULL;
    int64_t xn, xd;
    if (!fe_rational(e->data.function.args[1], &xn, &xd)) return NULL;
    /* Digamma at a positive non-integer rational: Gauss's theorem (the
     * reference: FunctionExpand[PolyGamma[0, 1/3]] =
     * -EulerGamma - Pi/(2 Sqrt[3]) - 3 Log[3]/2). */
    if (me->data.integer == 0)
        return (xn > 0 && xd > 1) ? polygamma_gauss_digamma(xn, xd) : NULL;
    int64_t n = me->data.integer + 1;
    Expr* M = expr_new_integer(me->data.integer);
    Expr* out = NULL;
    if (xd == 1 && xn >= 1 && xn <= 10000) {
        /* Positive integer k, every order (zeta(odd) included):
         * n! ((-1)^(n+1) (Zeta[n+1] - H_{k-1}^(n+1))), n! outside as in the
         * reference: FunctionExpand[PolyGamma[2, 3]] = 2 (9/8 - Zeta[3]). */
        Expr* args[2] = { M, expr_new_integer(xn) };
        out = fe_apply("(#1! ((-1)^(#1 + 1) (Zeta[#1 + 1] - HarmonicNumber[#2 - 1, #1 + 1])))&",
                       args, 2);
        expr_free(args[1]);
        expr_free(M);
        return out;                     /* already in the reference's form */
    }
    if (xd == 2 && xn > 1 && xn <= 20001) {
        /* Positive half-integer k + 1/2 past 1/2: (2^(n) - 1) Zeta[n] less the
         * recurrence terms.  Odd orders keep n! outside (6 (-16 + Pi^4/6)),
         * even orders (zeta(odd)) come out expanded (448/27 - 14 Zeta[3]), as
         * in the reference. */
        Expr* args[2] = { M, expr_new_integer((xn - 1) / 2) };
        out = fe_apply(me->data.integer % 2 == 1
            ? "(#1! ((2^(#1 + 1) - 1) Zeta[#1 + 1] - Sum[(jj + 1/2)^(-#1 - 1), {jj, 0, #2 - 1}]))&"
            : "Expand[-#1! ((2^(#1 + 1) - 1) Zeta[#1 + 1] - Sum[(jj + 1/2)^(-#1 - 1), {jj, 0, #2 - 1}])]&",
            args, 2);
        expr_free(args[1]);
        expr_free(M);
        return out;
    }
    if (xn == 1 && xd == 2) {
        out = fe_apply("((-1)^(#1 + 1) #1! (2^(#1 + 1) - 1) Zeta[#1 + 1])&", &M, 1);
    } else if (xd == 4 && (xn == 1 || xn == 3) && (n == 2 || n % 2 == 1)) {
        Expr* args[2] = { M, expr_new_integer(xn == 1 ? 1 : -1) };
        out = fe_apply(n == 2
            ? "((-1)^(#1 + 1) #1! ((4^(#1+1) - 2^(#1+1)) Zeta[#1+1] + #2 4^(#1+1) Catalan)/2)&"
            : "((-1)^(#1 + 1) #1! ((4^(#1+1) - 2^(#1+1)) Zeta[#1+1] + #2 4^(#1+1)"
              " (-1)^(#1/2) EulerE[#1] Pi^(#1+1) / (4^(#1/2 + 1) #1!))/2)&",
            args, 2);
        expr_free(args[1]);
    }
    expr_free(M);
    if (out) { Expr* x = fe_expand(out); expr_free(out); out = x; }
    return out;
}

/* ------------------------------------------------------------------ */
/*  PolyLog at roots of unity                                          */
/* ------------------------------------------------------------------ */

/* Clausen part of Li_n(e^{2 pi i p/q}) -- the cosine sum (odd n) or the sine
 * sum (even n) -- as an evaluated, owned Expr. */
static Expr* fe_clausen_part(int64_t n, int64_t p, int64_t q) {
    Expr* N = expr_new_integer(n);
    Expr* P = expr_new_integer(p);
    Expr* Q = expr_new_integer(q);
    Expr* args[3] = { N, P, Q };
    Expr* out = NULL;
    bool odd = (n % 2) != 0;
    if (odd && (q == 1 || q == 2 || q == 3 || q == 4 || q == 6)) {
        /* Rational cosines: the Ramanujan-sum closed form (p is primitive
         * mod q and the phi(q) <= 2 primitive residues share it). */
        out = fe_apply("(Zeta[#1] Sum[MoebiusMu[#3/dd] dd^(1 - #1), {dd, Divisors[#3]}]"
                       " / EulerPhi[#3])&", args, 3);
    } else if (!odd && (q == 1 || q == 2)) {
        out = expr_new_integer(0);
    } else if (!odd && n == 2 && q == 4) {
        /* Im Li_2(+-I) = +-Catalan (DirichletBeta[2]). */
        out = fe_apply("(Sin[2 Pi #2/#3] Catalan)&", args, 3);
    } else {
        out = fe_apply(odd
            ? "(#3^(-#1) Sum[Cos[2 Pi jj #2/#3] If[jj == #3, Zeta[#1],"
              " (-1)^#1 PolyGamma[#1 - 1, jj/#3]/(#1 - 1)!], {jj, 1, #3}])&"
            : "(#3^(-#1) Sum[Sin[2 Pi jj #2/#3] If[jj == #3, 0,"
              " (-1)^#1 PolyGamma[#1 - 1, jj/#3]/(#1 - 1)!], {jj, 1, #3}])&",
            args, 3);
    }
    expr_free(N); expr_free(P); expr_free(Q);
    return out;
}

/* PolyLog[n, w] -> R + I S for integer n >= 2 and an exact root of unity w;
 * NULL (nothing owned) when the call is not of that form. */
static Expr* fe_polylog_root(const Expr* e) {
    if (!fe_is_head(e, SYM_PolyLog) || e->data.function.arg_count != 2) return NULL;
    const Expr* ne = e->data.function.args[0];
    if (ne->type != EXPR_INTEGER || ne->data.integer < 2 || ne->data.integer > 200) return NULL;
    int64_t rn, rd;
    if (!fe_root_of_unity(e->data.function.args[1], &rn, &rd)) return NULL;
    int64_t n = ne->data.integer;
    /* x = r/2 mod 1 in [0, 1): p/q. */
    int64_t p = rn, q = 2 * rd;
    if (p < 0) p += q;
    int64_t g = fe_gcd(p, q);
    if (g == 0) g = 1;
    p /= g; q /= g;
    if (p == 0) q = 1;
    Expr* N = expr_new_integer(n);
    Expr* X = make_rational(p, q);
    Expr* args[2] = { N, X };
    bool odd = (n % 2) != 0;
    Expr* bern = fe_apply(odd
        ? "((-1)^((#1 + 1)/2) (2 Pi)^#1 BernoulliB[#1, #2] / (2 #1!))&"
        : "((-1)^(#1/2 + 1) (2 Pi)^#1 BernoulliB[#1, #2] / (2 #1!))&", args, 2);
    expr_free(N); expr_free(X);
    /* The Clausen part from the representative p <= q/2: cosine sums are even
     * in p -> q - p and sine sums odd, so the conjugate root w-bar yields the
     * SAME expression (times -1 for a sine sum) and a conjugate pair cancels
     * structurally, not merely numerically. */
    Expr* cl;
    if (2 * p > q) {
        cl = fe_clausen_part(n, q - p, q);
        if (cl && !odd) {
            Expr* m1 = expr_new_integer(-1);
            Expr* ng = fe_apply("(-#1)&", &cl, 1);
            expr_free(m1); expr_free(cl); cl = ng;
        }
    } else {
        cl = fe_clausen_part(n, p, q);
    }
    if (!bern || !cl) { expr_free(bern); expr_free(cl); return NULL; }
    /* even n: R = bern, S = cl;  odd n: R = cl, S = bern. */
    Expr* re = odd ? cl : bern;
    Expr* im = odd ? bern : cl;
    Expr* ri[2] = { re, im };
    Expr* out = fe_apply("(#1 + I #2)&", ri, 2);
    expr_free(re); expr_free(im);
    /* The Hurwitz terms at 1/2, 1/4, 3/4 have closed forms of their own. */
    if (out) {
        bool ch = false;
        Expr* w2 = fe_walk(out, &ch, FE_PGAMMA);
        expr_free(out);
        out = w2;
    }
    return out;
}

/* PolyLog[2, 2] and PolyLog[3, 2] by the inversion formula
 * Li_n(z) + (-1)^n Li_n(1/z) = -(2 Pi I)^n/n! B_n(1/2 + Log[-z]/(2 Pi I)) at
 * z = 2, with Li_n(1/2) closed for n = 2, 3 and Log[-2] = Log[2] + I Pi (the
 * branch PolyLog itself evaluates on: the value approached from below the cut):
 *     PolyLog[2, 2] = Pi^2/4 - I Pi Log[2],
 *     PolyLog[3, 2] = 7/8 Zeta[3] + Pi^2/4 Log[2] - I Pi Log[2]^2/2.
 * NULL (nothing owned) for any other call. */
static Expr* fe_polylog_two(const Expr* e) {
    if (!fe_is_head(e, SYM_PolyLog) || e->data.function.arg_count != 2) return NULL;
    const Expr* ne = e->data.function.args[0];
    const Expr* z = e->data.function.args[1];
    if (ne->type != EXPR_INTEGER || z->type != EXPR_INTEGER || z->data.integer != 2)
        return NULL;
    if (ne->data.integer == 2)
        return fe_apply("(Pi^2/4 - I Pi Log[2])&", NULL, 0);
    if (ne->data.integer == 3)
        return fe_apply("(7/8 Zeta[3] + Pi^2/4 Log[2] - I Pi Log[2]^2/2)&", NULL, 0);
    return NULL;
}

/* PolyLog[2, z] at an exact Gaussian rational z whose orbit under the six
 * anharmonic maps z -> {1 - z, 1/z, z/(z - 1), 1/(1 - z), (z - 1)/z} meets a
 * closed point w: a root of unity (fe_polylog_root), 1/2 or 2.  Li_2(z) is
 * then Li_2(w) plus elementary terms, by the reflection, Landen and inversion
 * formulas and their compositions:
 *     w = 1 - z:        Pi^2/6 - Log[z] Log[1 - z] - Li_2(w)
 *     w = z/(z - 1):    -Log[1 - z]^2/2 - Li_2(w)
 *     w = 1/z:          -Pi^2/6 - Log[-z]^2/2 - Li_2(w)
 *     w = 1/(1 - z):    Pi^2/3 - Log[z] Log[1 - z] + Log[z - 1]^2/2 + Li_2(w)
 *     w = (z - 1)/z:    -Pi^2/3 - Log[-z]^2/2 + Log[1/z] Log[1 - 1/z] + Li_2(w)
 * Each holds off the branch cuts of the logs and dilogarithms involved; rather
 * than track the cuts case by case, the candidate (its Logs split by
 * ComplexExpand into real and imaginary parts) is accepted only when it agrees
 * with PolyLog[2, z] to 30 digits -- a wrong branch is off by a multiple of
 * 2 Pi I Log or Pi^2, never close.  Li_2((1 + I)/2) = 5 Pi^2/96 - Log[2]^2/8
 * + I (Catalan - Pi Log[2]/8), Li_2(1 + I) = Pi^2/16 + I (Catalan + Pi Log[2]/4).
 * NULL (nothing owned) when no image is closed. */
static int fe_li2_depth = 0;

static bool fe_li2_closed_point(const Expr* w) {
    int64_t rn, rd, n, d;
    if (fe_root_of_unity(w, &rn, &rd)) return true;
    if (fe_rational(w, &n, &d) && ((n == 1 && d == 2) || (n == 2 && d == 1))) return true;
    return false;
}

/* An exact numeric constant (rationals, Gaussian rationals, radicals, roots
 * of unity): NumericQ with infinite precision. */
static bool fe_exact_numeric(const Expr* z) {
    int64_t n, d;
    if (fe_rational(z, &n, &d)) return true;
    if (z->type != EXPR_FUNCTION) return false;
    Expr* zc = expr_copy((Expr*)z);
    Expr* ok = fe_apply("(NumericQ[#1] && Precision[#1] === Infinity)&", &zc, 1);
    expr_free(zc);
    bool r = ok && ok->type == EXPR_SYMBOL && ok->data.symbol.name == SYM_True;
    expr_free(ok);
    return r;
}

static Expr* fe_polylog2_anharmonic(const Expr* e) {
    if (!fe_is_head(e, SYM_PolyLog) || e->data.function.arg_count != 2) return NULL;
    const Expr* ne = e->data.function.args[0];
    const Expr* z = e->data.function.args[1];
    if (ne->type != EXPR_INTEGER || ne->data.integer != 2) return NULL;
    if (fe_li2_depth > 0) return NULL;
    if (fe_li2_closed_point(z) || !fe_exact_numeric(z)) return NULL;
    static const char* const maps[5] = {
        "(1 - #1)&", "(#1/(#1 - 1))&", "(1/#1)&", "(1/(1 - #1))&", "((#1 - 1)/#1)&" };
    static const char* const forms[5] = {
        "(Pi^2/6 - Log[#1] Log[1 - #1] - PolyLog[2, #2])&",
        "(-Log[1 - #1]^2/2 - PolyLog[2, #2])&",
        "(-Pi^2/6 - Log[-#1]^2/2 - PolyLog[2, #2])&",
        "(Pi^2/3 - Log[#1] Log[1 - #1] + Log[#1 - 1]^2/2 + PolyLog[2, #2])&",
        "(-Pi^2/3 - Log[-#1]^2/2 + Log[1/#1] Log[1 - 1/#1] + PolyLog[2, #2])&" };
    Expr* zc = expr_copy((Expr*)z);
    Expr* out = NULL;
    fe_li2_depth++;
    for (int k = 0; k < 5 && !out; k++) {
        /* z = 0 or 1 make some maps singular; those z never get here (they
         * evaluate), and 1/0 would only yield ComplexInfinity, not closed. */
        Expr* w = fe_apply(maps[k], &zc, 1);
        if (!w || !fe_li2_closed_point(w)) { expr_free(w); continue; }
        Expr* args[2] = { zc, w };
        Expr* cand = fe_apply(forms[k], args, 2);
        expr_free(w);
        if (!cand) continue;
        bool ch = false;
        Expr* red = fe_walk(cand, &ch, FE_ROOTS);
        expr_free(cand);
        Expr* t0 = fe_apply("(ComplexExpand[#1])&", &red, 1);
        expr_free(red);
        if (!t0) continue;
        bool lr = false;                      /* Log[1/2] -> -Log[2] */
        Expr* t1 = fe_log_reciprocal(t0, &lr);
        expr_free(t0);
        Expr* tidy = fe_expand(t1);
        expr_free(t1);
        if (!tidy) continue;
        Expr* chk[2] = { tidy, zc };
        /* both sides at 40 digits, then compared: N of the (cancelling)
         * difference itself would chase 40 significant digits of ~0 */
        Expr* ok = fe_apply("(FreeQ[#1, PolyLog] && "
                            "TrueQ[Abs[N[#1, 40] - N[PolyLog[2, #2], 40]] < 10^-30])&", chk, 2);
        bool good = ok && ok->type == EXPR_SYMBOL && ok->data.symbol.name == SYM_True;
        expr_free(ok);
        if (good) out = tidy; else expr_free(tidy);
    }
    fe_li2_depth--;
    expr_free(zc);
    return out;
}

/* Log[w] -> I Pi r, Log[1 + w] -> Log[2 Cos[Pi r/2]] + I Pi r/2 (|r| < 1). */
static Expr* fe_log_root(const Expr* e) {
    if (!fe_is_head(e, SYM_Log) || e->data.function.arg_count != 1) return NULL;
    const Expr* c = e->data.function.args[0];
    int64_t rn, rd;
    if (fe_is_head(c, SYM_Plus) && c->data.function.arg_count == 2 &&
        c->data.function.args[0]->type == EXPR_INTEGER &&
        c->data.function.args[0]->data.integer == 1 &&
        fe_root_of_unity(c->data.function.args[1], &rn, &rd)) {
        if (rn == rd) return NULL;               /* 1 + (-1): Log[0] */
        Expr* R = make_rational(rn, rd);
        Expr* out = fe_apply("(Log[2 Cos[Pi #1/2]] + I Pi #1/2)&", &R, 1);
        expr_free(R);
        return out;
    }
    if (fe_root_of_unity(c, &rn, &rd) && rn != 0) {
        Expr* R = make_rational(rn, rd);
        Expr* out = fe_apply("(I Pi #1)&", &R, 1);
        expr_free(R);
        return out;
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  The walk                                                           */
/* ------------------------------------------------------------------ */


/* HypergeometricPFQ[{a..}, {b..}, z] with every parameter put over one
 * denominator by Together, evaluated, so that parameters which agree only
 * up to that form cancel (1/2 (2 + 2 v) against 1 + v); NULL when no
 * parameter changes. Owned. */
static Expr* fe_pfq_together(const Expr* e) {
    if (e->data.function.arg_count != 3) return NULL;
    bool any = false;
    Expr* lists[2];
    for (int k = 0; k < 2; k++) {
        const Expr* L = e->data.function.args[k];
        if (!fe_is_head(L, SYM_List)) {
            for (int j = 0; j < k; j++) expr_free(lists[j]);
            return NULL;
        }
        size_t n = L->data.function.arg_count;
        Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) {
            const Expr* par = L->data.function.args[i];
            v[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                                                   (Expr*[]){ expr_copy((Expr*)par) }, 1));
            if (!expr_eq(v[i], par)) any = true;
        }
        lists[k] = expr_new_function(expr_new_symbol(SYM_List), v, n);
        free(v);
    }
    if (!any) { expr_free(lists[0]); expr_free(lists[1]); return NULL; }
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_HypergeometricPFQ),
        (Expr*[]){ lists[0], lists[1], expr_copy(e->data.function.args[2]) }, 3));
}

/* ParabolicCylinderD[n, z] at an exact integer order (|n| <= 200):
 *   n >= 0:  2^(-n/2) E^(-z^2/4) HermiteH[n, z/Sqrt[2]]            (DLMF 12.7.2)
 *   n <  0:  D_0 = E^(-z^2/4), D_-1 = E^(z^2/4) Sqrt[Pi/2] Erfc[z/Sqrt[2]]
 *            (DLMF 12.7.5) and D_{k-1} = (z D_k - D_{k+1})/k downwards
 *            (the recurrence 12.8.1). NULL for any other node. */
static Expr* fe_pcfd_integer(const Expr* e) {
    if (!fe_is_head(e, SYM_ParabolicCylinderD) || e->data.function.arg_count != 2) return NULL;
    const Expr* nu = e->data.function.args[0];
    if (nu->type != EXPR_INTEGER || nu->data.integer > 200 || nu->data.integer < -200)
        return NULL;
    Expr* const* args = e->data.function.args;
    if (nu->data.integer >= 0)
        return fe_apply("(2^(-#1/2) Exp[-#2^2/4] HermiteH[#1, #2/Sqrt[2]])&", args, 2);
    return fe_apply(
        "Module[{d0 = Exp[-#2^2/4], d1 = Exp[#2^2/4] Sqrt[Pi/2] Erfc[#2/Sqrt[2]], t},"
        " Do[t = Expand[(#2 d1 - d0)/k]; d0 = d1; d1 = t, {k, -1, #1 + 1, -1}]; d1]&",
        args, 2);
}

/* A copy of e with every expandable node replaced, bottom-up; sets *changed
 * when anything was rewritten. `which` selects the rule families. Returns an
 * owned Expr. */
/* ---- native builders for the point tables (arguments consumed) ---- */
static Expr* fp_fn(const char* h, int n, ...) {
    Expr* a[8];
    va_list ap;
    va_start(ap, n);
    for (int i = 0; i < n; i++) a[i] = va_arg(ap, Expr*);
    va_end(ap);
    return expr_new_function(expr_new_symbol(h), a, (size_t)n);
}
static Expr* fp_i(int64_t v) { return expr_new_integer(v); }
static Expr* fp_q(int64_t n, int64_t d) { return make_rational(n, d); }
static Expr* fp_pow(Expr* b, Expr* e) { return fp_fn(SYM_Power, 2, b, e); }
static Expr* fp_sqrt(Expr* b) { return fp_pow(b, fp_q(1, 2)); }
static Expr* fp_pi(void) { return expr_new_symbol(SYM_Pi); }
static Expr* fp_eval(Expr* e) { Expr* r = evaluate(e); expr_free(e); return r; }
/* phi = (1 + Sqrt[5])/2 */
static Expr* fp_phi(void) {
    return fp_fn(SYM_Times, 2, fp_q(1, 2), fp_fn(SYM_Plus, 2, fp_i(1), fp_sqrt(fp_i(5))));
}

/* The value of an exact numeric z as a double (real part), or NAN. */
static double fp_value(const Expr* z) {
    Expr* n = numericalize_machine_unplanned(z);
    double v = NAN;
    if (n && n->type == EXPR_REAL) v = n->data.real;
    else if (n && n->type == EXPR_INTEGER) v = (double)n->data.integer;
    if (n) expr_free(n);
    return v;
}

/* A table entry: the point's double value (the screen), and builders for the
 * exact point and the closed form. */
typedef struct { double at; Expr* (*point)(void); Expr* (*value)(void); } FePoint;

/* The closed form of the first entry whose point equals z exactly.  The
 * double screen leaves at most one candidate; RootReduce[z - point] === 0
 * is then the exact test. */
static Expr* fe_point_table(const Expr* z, const FePoint* tab, size_t n) {
    double v = fp_value(z);
    if (v != v) return NULL;
    for (size_t i = 0; i < n; i++) {
        if (fabs(v - tab[i].at) > 1e-12 * (1 + fabs(tab[i].at))) continue;
        Expr* d = fp_eval(fp_fn("RootReduce", 1, fp_fn(SYM_Plus, 2, expr_copy((Expr*)z),
                                 fp_fn(SYM_Times, 2, fp_i(-1), tab[i].point()))));
        bool hit = d && d->type == EXPR_INTEGER && d->data.integer == 0;
        if (d) expr_free(d);
        if (hit) return fp_eval(tab[i].value());
    }
    return NULL;
}

/* PolyLog[2, z] at the golden-ratio points (Landen): with phi = (1+Sqrt[5])/2,
 *   Li2(phi^-2) = Pi^2/15 - Log[phi]^2,     Li2(phi^-1) = Pi^2/10 - Log[phi]^2,
 *   Li2(-phi^-1) = -Pi^2/15 + Log[phi]^2/2, Li2(-phi) = -Pi^2/10 - Log[phi]^2,
 * the only algebraic points besides 0, +-1, 1/2 where Li2 is known in closed
 * form (each verified to 30 digits against the series). */
static Expr* fp_lg2(void) { return fp_pow(fp_fn(SYM_Log, 1, fp_phi()), fp_i(2)); }
static Expr* fp_pi2(int64_t n, int64_t d) { return fp_fn(SYM_Times, 2, fp_q(n, d), fp_pow(fp_pi(), fp_i(2))); }
static Expr* fg_p1(void) { return fp_pow(fp_phi(), fp_i(-2)); }
static Expr* fg_v1(void) { return fp_fn(SYM_Plus, 2, fp_pi2(1, 15), fp_fn(SYM_Times, 2, fp_i(-1), fp_lg2())); }
static Expr* fg_p2(void) { return fp_pow(fp_phi(), fp_i(-1)); }
static Expr* fg_v2(void) { return fp_fn(SYM_Plus, 2, fp_pi2(1, 10), fp_fn(SYM_Times, 2, fp_i(-1), fp_lg2())); }
static Expr* fg_p3(void) { return fp_fn(SYM_Times, 2, fp_i(-1), fp_pow(fp_phi(), fp_i(-1))); }
static Expr* fg_v3(void) { return fp_fn(SYM_Plus, 2, fp_pi2(-1, 15), fp_fn(SYM_Times, 2, fp_q(1, 2), fp_lg2())); }
static Expr* fg_p4(void) { return fp_fn(SYM_Times, 2, fp_i(-1), fp_phi()); }
static Expr* fg_v4(void) { return fp_fn(SYM_Plus, 2, fp_pi2(-1, 10), fp_fn(SYM_Times, 2, fp_i(-1), fp_lg2())); }

static Expr* fe_polylog_golden(const Expr* e) {
    if (!fe_is_head(e, SYM_PolyLog) || e->data.function.arg_count != 2) return NULL;
    const Expr* ne = e->data.function.args[0];
    if (ne->type != EXPR_INTEGER || ne->data.integer != 2) return NULL;
    const Expr* z = e->data.function.args[1];
    if (!fe_is_head((Expr*)z, SYM_Plus) && !fe_is_head((Expr*)z, SYM_Times)) return NULL;
    const double phi = 1.6180339887498948482;
    const FePoint tab[4] = {
        { 1 / (phi * phi), fg_p1, fg_v1 }, { 1 / phi, fg_p2, fg_v2 },
        { -1 / phi, fg_p3, fg_v3 },        { -phi, fg_p4, fg_v4 } };
    return fe_point_table(z, tab, 4);
}

/* EllipticK at the singular moduli k_r^2 with class number h(-4r) = 1 that are
 * not automatic (K(1/2), K(-1) evaluate on their own), and their complements
 * K(1 - k_r^2) = Sqrt[r] K(k_r^2) (Chowla-Selberg; each verified to 30 digits):
 *   r = 3: m = (2 - Sqrt[3])/4,  K = 3^(1/4) Gamma[1/3]^3/(2^(7/3) Pi)
 *   r = 2: m = (Sqrt[2] - 1)^2,  K = Sqrt[Sqrt[2] + 1] Gamma[1/8] Gamma[3/8]/(2^(13/4) Sqrt[Pi]) */
static Expr* fk_k3(void) {
    return fp_fn(SYM_Times, 4, fp_pow(fp_i(3), fp_q(1, 4)),
                 fp_pow(fp_fn(SYM_Gamma, 1, fp_q(1, 3)), fp_i(3)),
                 fp_pow(fp_i(2), fp_q(-7, 3)), fp_pow(fp_pi(), fp_i(-1)));
}
static Expr* fk_k2(void) {
    return fp_fn(SYM_Times, 5, fp_sqrt(fp_fn(SYM_Plus, 2, fp_sqrt(fp_i(2)), fp_i(1))),
                 fp_fn(SYM_Gamma, 1, fp_q(1, 8)), fp_fn(SYM_Gamma, 1, fp_q(3, 8)),
                 fp_pow(fp_i(2), fp_q(-13, 4)), fp_pow(fp_pi(), fp_q(-1, 2)));
}
static Expr* fk_s21(void) {     /* (Sqrt[2] - 1)^2 */
    return fp_pow(fp_fn(SYM_Plus, 2, fp_sqrt(fp_i(2)), fp_i(-1)), fp_i(2));
}
static Expr* fk_p1(void) { return fp_fn(SYM_Times, 2, fp_q(1, 4), fp_fn(SYM_Plus, 2, fp_i(2), fp_fn(SYM_Times, 2, fp_i(-1), fp_sqrt(fp_i(3))))); }
static Expr* fk_p2(void) { return fp_fn(SYM_Times, 2, fp_q(1, 4), fp_fn(SYM_Plus, 2, fp_i(2), fp_sqrt(fp_i(3)))); }
static Expr* fk_v2(void) { return fp_fn(SYM_Times, 2, fp_sqrt(fp_i(3)), fk_k3()); }
static Expr* fk_p4(void) { return fp_fn(SYM_Plus, 2, fp_i(1), fp_fn(SYM_Times, 2, fp_i(-1), fk_s21())); }
static Expr* fk_v4(void) { return fp_fn(SYM_Times, 2, fp_sqrt(fp_i(2)), fk_k2()); }

static Expr* fe_elliptick_singular(const Expr* e) {
    if (!fe_is_head(e, SYM_EllipticK) || e->data.function.arg_count != 1) return NULL;
    const Expr* m = e->data.function.args[0];
    if (!fe_is_head((Expr*)m, SYM_Plus) && !fe_is_head((Expr*)m, SYM_Times)
        && !fe_is_head((Expr*)m, SYM_Power)) return NULL;
    const double s3 = 1.7320508075688772935, s21 = 0.17157287525380990239;
    const FePoint tab[4] = {
        { (2 - s3) / 4, fk_p1, fk_k3 }, { (2 + s3) / 4, fk_p2, fk_v2 },
        { s21, fk_s21, fk_k2 },          { 1 - s21, fk_p4, fk_v4 } };
    return fe_point_table(m, tab, 4);
}

/* ---- Cos at rational multiples of Pi with a Fermat-prime denominator ------
 *
 * Gauss's construction.  For a prime p with p - 1 = 2^k, the Galois group of
 * Q(zeta_p) is cyclic of order 2^k, generated by zeta -> zeta^g (g a primitive
 * root).  The Gaussian periods of the subgroup of index 2^j,
 *     eta_{j,i} = sum over s == i (mod 2^j) of zeta^(g^s),
 * halve at every level: the two children a, b of a level-(j-1) period are the
 * roots of a quadratic, eta_{a,b} = (eta_parent +- Sqrt[Delta])/2, where
 * Delta = (eta_a - eta_b)^2 is fixed by the parent subgroup and so an INTEGER
 * combination of level-(j-1) periods -- read off by counting exponent sums
 * (no product of radicals is ever formed, so the radicands stay linear in
 * the known periods).  The sign is the numeric order of the two children.
 * At level k-1 the period is zeta^t + zeta^-t = 2 Cos[2 Pi t/p]:
 *   Cos[2 Pi/17] = (-1 + Sqrt[17] + Sqrt[34 - 2 Sqrt[17]]
 *                   + 2 Sqrt[17 + 3 Sqrt[17] - ...])/16.
 * Cos[a Pi/p] reduces to Cos[2 Pi t/p] by parity; Cos[a Pi/(2 p)] by the
 * half-angle formula.  p = 3, 5 are evaluated by Cos itself; p <= 257 here. */
static bool fe_fermat_prime(int64_t p) {
    if (p < 17 || p > 257) return false;
    for (int64_t d = 2; d * d <= p; d++) if (p % d == 0) return false;
    int64_t q = p - 1;
    return (q & (q - 1)) == 0;
}

static Expr* fe_g2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* fe_ev(Expr* e) { Expr* r = evaluate(e); expr_free(e); return r; }
/* Expand[e]: distribute the rational halves so like radicals combine. */
static Expr* fe_gexpand(Expr* e) {
    return fe_ev(expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){ e }, 1));
}

/* Sqrt[delta] with delta's rational coefficients cleared by a square:
 * Sqrt[17/2 - Sqrt[17]/2] -> Sqrt[34 - 2 Sqrt[17]]/2.  Consumes delta. */
static int64_t fe_lcm(int64_t a, int64_t b) {
    int64_t x = a, y = b; while (y) { int64_t t = x % y; x = y; y = t; }
    return a / x * b;
}
static Expr* fe_sqrt_integral(Expr* delta) {
    int64_t L = 1;
    size_t nt = fe_is_head(delta, SYM_Plus) ? delta->data.function.arg_count : 1;
    Expr* const* t = fe_is_head(delta, SYM_Plus) ? delta->data.function.args : &delta;
    for (size_t i = 0; i < nt && L < 1000000; i++) {
        const Expr* c = t[i];
        if (fe_is_head((Expr*)c, SYM_Times) && c->data.function.arg_count >= 1) c = c->data.function.args[0];
        int64_t nn, dd;
        if (is_rational(c, &nn, &dd) && dd > 1) L = fe_lcm(L, dd);
    }
    /* smallest q with L | q^2 */
    int64_t q = 1, r = L;
    for (int64_t pr = 2; pr * pr <= r || r > 1; pr++) {
        if (pr * pr > r) pr = r;
        int v = 0; while (r % pr == 0) { r /= pr; v++; }
        for (int k = 0; k < (v + 1) / 2; k++) q *= pr;
        if (r == 1) break;
    }
    if (q == 1) return fe_g2(SYM_Power, delta, make_rational(1, 2));
    Expr* scaled = fe_gexpand(fe_g2(SYM_Times, expr_new_integer(q * q), delta));
    return fe_g2(SYM_Times, make_rational(1, q), fe_g2(SYM_Power, scaled, make_rational(1, 2)));
}

/* 2 Cos[2 Pi t/p] in radicals (t not == 0 mod p). */
static Expr* fe_gauss_2cos(int64_t p, int64_t t) {
    int k = 0; while (((int64_t)1 << k) < p - 1) k++;
    int n = (int)(p - 1);
    int* pw = malloc(sizeof(int) * (size_t)n);      /* g^s mod p */
    int* lg = malloc(sizeof(int) * (size_t)p);      /* log_g      */
    int g = 2;
    for (;; g++) {                                   /* primitive root */
        int64_t x = 1; bool prim = true;
        for (int s = 1; s < n; s++) { x = x * g % p; if (x == 1) { prim = false; break; } }
        if (prim) break;
    }
    { int64_t x = 1; for (int s = 0; s < n; s++) { pw[s] = (int)x; lg[x] = s; x = x * g % p; } }
    /* periods by level: E[j][i], val[j][i] */
    Expr*** E = malloc(sizeof(Expr**) * (size_t)k);
    double** V = malloc(sizeof(double*) * (size_t)k);
    for (int j = 0; j < k; j++) {
        int np = 1 << j;
        E[j] = calloc((size_t)np, sizeof(Expr*));
        V[j] = calloc((size_t)np, sizeof(double));
        for (int s = 0; s < n; s++) V[j][s % np] += cos(2.0 * M_PI * pw[s] / (double)p);
    }
    E[0][0] = expr_new_integer(-1);
    int64_t* cnt = malloc(sizeof(int64_t) * (size_t)p);
    bool ok = true;
    for (int j = 1; j < k && ok; j++) {
        int half = 1 << (j - 1), np = 1 << j;
        for (int i = 0; i < half && ok; i++) {
            int a = i, b = i + half;
            /* Delta = (eta_a - eta_b)^2 = sum s_x s_y zeta^(x+y) */
            memset(cnt, 0, sizeof(int64_t) * (size_t)p);
            for (int sx = 0; sx < n; sx++) {
                int cx = sx % np; if (cx != a && cx != b) continue;
                int sgx = cx == a ? 1 : -1;
                for (int sy = 0; sy < n; sy++) {
                    int cy = sy % np; if (cy != a && cy != b) continue;
                    int sgy = cy == a ? 1 : -1;
                    cnt[(pw[sx] + pw[sy]) % p] += sgx * sgy;
                }
            }
            /* group by level-(j-1) class; counts must agree within a class */
            int64_t* cl = malloc(sizeof(int64_t) * (size_t)half);
            bool* seen = calloc((size_t)half, sizeof(bool));
            for (int r = 1; r < p && ok; r++) {
                int c = lg[r] % half;
                if (!seen[c]) { seen[c] = true; cl[c] = cnt[r]; }
                else if (cl[c] != cnt[r]) ok = false;
            }
            if (ok) {
                Expr** tm = malloc(sizeof(Expr*) * (size_t)(half + 1));
                int w = 0;
                tm[w++] = expr_new_integer(cnt[0]);
                for (int c = 0; c < half; c++)
                    if (cl[c]) tm[w++] = fe_g2(SYM_Times, expr_new_integer(cl[c]), expr_copy(E[j - 1][c]));
                Expr* delta = fe_gexpand(expr_new_function(expr_new_symbol(SYM_Plus), tm, (size_t)w));
                free(tm);
                Expr* sq = fe_sqrt_integral(delta);
                int sg = V[j][a] >= V[j][b] ? 1 : -1;
                for (int ch = 0; ch < 2; ch++) {
                    int sign = ch == 0 ? sg : -sg;
                    Expr* v = fe_g2(SYM_Times, make_rational(1, 2),
                        fe_g2(SYM_Plus, expr_copy(E[j - 1][i]),
                              fe_g2(SYM_Times, expr_new_integer(sign), expr_copy(sq))));
                    E[j][ch == 0 ? a : b] = fe_gexpand(v);
                }
                expr_free(sq);
            }
            free(cl); free(seen);
        }
    }
    Expr* out = NULL;
    if (ok) {
        int idx = lg[((t % p) + p) % p] % (1 << (k - 1));
        out = expr_copy(E[k - 1][idx]);
    }
    for (int j = 0; j < k; j++) {
        for (int i = 0; i < (1 << j); i++) if (E[j][i]) expr_free(E[j][i]);
        free(E[j]); free(V[j]);
    }
    free(E); free(V); free(cnt); free(pw); free(lg);
    return out;
}

static Expr* fe_cos_gauss(const Expr* e) {
    if (!fe_is_head(e, SYM_Cos) || e->data.function.arg_count != 1) return NULL;
    const Expr* arg = e->data.function.args[0];
    /* arg == (a/b) Pi */
    int64_t a, b;
    if (!(fe_is_head((Expr*)arg, SYM_Times) && arg->data.function.arg_count == 2
          && arg->data.function.args[1]->type == EXPR_SYMBOL
          && arg->data.function.args[1]->data.symbol.name == SYM_Pi
          && is_rational(arg->data.function.args[0], &a, &b))) return NULL;
    bool half = false;
    int64_t p = b;
    if (b % 2 == 0 && fe_fermat_prime(b / 2)) { p = b / 2; half = true; }
    if (!fe_fermat_prime(p)) return NULL;
    /* Cos[a' Pi/p] with a' = a (or the doubled angle for the half-angle case) */
    int64_t ap = a;                     /* angle a' Pi / p */
    int64_t t; int sign = 1;
    ap %= 2 * p; if (ap < 0) ap += 2 * p;
    if (ap % 2 == 0) t = ap / 2;
    else { t = (ap - p) / 2; sign = -1; }     /* cos(x) = -cos(x - Pi) */
    if (((t % p) + p) % p == 0) return NULL;
    Expr* c2 = fe_gauss_2cos(p, t);          /* 2 Cos[a' Pi/p] * sign */
    if (!c2) return NULL;
    Expr* c = fe_gexpand(fe_g2(SYM_Times, make_rational(sign, 2), c2));
    if (!half) return c;
    /* Cos[a Pi/(2p)] = s Sqrt[(1 + Cos[a Pi/p])/2] */
    double ang = M_PI * (double)a / (double)b;
    int s = cos(ang) >= 0 ? 1 : -1;
    return fe_ev(fe_g2(SYM_Times, expr_new_integer(s),
        fe_g2(SYM_Power, fe_g2(SYM_Times, make_rational(1, 2), fe_g2(SYM_Plus, expr_new_integer(1), c)),
              make_rational(1, 2))));
}

static Expr* fe_walk(const Expr* e, bool* changed, int which) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        bool c = false;
        a[i] = fe_walk(e->data.function.args[i], &c, which);
        any |= c;
    }
    Expr* node = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (any) {
        *changed = true;
        node = eval_and_free(node);
    }
    Expr* s = NULL;
    if ((which & FE_MEIJER) && fe_is_head(node, SYM_MeijerG)) {
        /* The G walk first: it handles the log case and names the result
         * directly (K, Y, Ci, Gamma(a, z), ...). Then the Slater route. */
        if (which & FE_PFQ) s = hrep_meijerg_expr(node);
        if (s) { expr_free(node); *changed = true; return s; }
        s = meijerg_function_expand(node);
        /* The Slater expansion leaves pFq terms. Each goes through Roach's
         * Formula algorithm (Roach '97: a G function without the log case
         * is a sum of pFq), so a MeijerG reduces to named functions. */
        if (s && (which & FE_PFQ)) {
            bool c2 = false;
            Expr* s2 = fe_walk(s, &c2, FE_PFQ);
            if (c2) { expr_free(s); s = eval_and_free(s2); }
            else expr_free(s2);
        }
    }
    if (!s && (which & FE_PFQ) && fe_is_head(node, SYM_HypergeometricPFQ)) {
        Expr* t = fe_pfq_together(node);
        if (t) {
            expr_free(node); node = t; *changed = true;
        }
        if (fe_is_head(node, SYM_HypergeometricPFQ)) s = hrep_pfq_expr(node);
    }
    if (!s && (which & FE_ROOTS)) {
        s = fe_polylog_root(node);
        if (!s) s = fe_polylog_two(node);
        if (!s) s = fe_polylog2_anharmonic(node);
        if (!s) s = fe_polylog_golden(node);
        if (!s) s = fe_cos_gauss(node);
        if (!s) s = fe_elliptick_singular(node);
        if (!s) s = fe_log_root(node);
        if (!s) s = fe_pcfd_integer(node);
    }
    if (!s && (which & (FE_ROOTS | FE_PGAMMA))) s = fe_polygamma_closed(node);
    if (s) { expr_free(node); node = s; *changed = true; }
    return node;
}

/* Split a term c P (P a sum, c the product of the other factors) into its
 * sum factor P (borrowed) and coefficient c (owned); NULL when the term has
 * not exactly one sum factor. */
static const Expr* fe_term_sum(const Expr* t, Expr** coef) {
    if (!fe_is_head(t, SYM_Times)) return NULL;
    const Expr* sum = NULL;
    size_t n = t->data.function.arg_count;
    Expr** rest = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* f = t->data.function.args[i];
        if (fe_is_head(f, SYM_Plus)) {
            if (sum) { for (size_t j = 0; j < k; j++) expr_free(rest[j]); free(rest); return NULL; }
            sum = f;
        } else {
            rest[k++] = expr_copy((Expr*)f);
        }
    }
    if (!sum) { for (size_t j = 0; j < k; j++) expr_free(rest[j]); free(rest); return NULL; }
    *coef = k == 0 ? expr_new_integer(1)
          : k == 1 ? rest[0]
          : expr_new_function(expr_new_symbol(SYM_Times), rest, k);
    free(rest);
    return sum;
}

/* Drop the pairs c P + c Q of a sum whose sum factors cancel, P + Q = 0 --
 * the conjugate Clausen parts of PolyLog at a root of unity are built so
 * they do, but -(c P) is c (-P) (the sign distributes over P), which Plus
 * does not collect against c P. Only the cancelling pairs go; every other
 * term keeps its form. Consumes `sum`. */
static Expr* fe_cancel_pairs(Expr* sum) {
    if (!fe_is_head(sum, SYM_Plus)) return sum;
    size_t n = sum->data.function.arg_count;
    bool* gone = calloc(n ? n : 1, sizeof(bool));
    bool any = false;
    for (size_t i = 0; i < n; i++) {
        if (gone[i]) continue;
        Expr* ci = NULL;
        const Expr* pi = fe_term_sum(sum->data.function.args[i], &ci);
        if (!pi) continue;
        for (size_t j = i + 1; j < n && !gone[i]; j++) {
            if (gone[j]) continue;
            Expr* cj = NULL;
            const Expr* pj = fe_term_sum(sum->data.function.args[j], &cj);
            if (!pj) continue;
            if (expr_eq(ci, cj)) {
                Expr* z = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
                    (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_copy((Expr*)pi), expr_copy((Expr*)pj) }, 2) }, 1));
                if (z && z->type == EXPR_INTEGER && z->data.integer == 0) {
                    gone[i] = gone[j] = true;
                    any = true;
                }
                if (z) expr_free(z);
            }
            expr_free(cj);
        }
        expr_free(ci);
    }
    if (!any) { free(gone); return sum; }
    Expr** keep = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        if (!gone[i]) keep[m++] = expr_copy(sum->data.function.args[i]);
    free(gone);
    expr_free(sum);
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), keep, m));
    free(keep);
    return out;
}

Expr* builtin_function_expand(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na != 1 && na != 2) return NULL;
    bool changed = false;
    Expr* out = fe_walk(res->data.function.args[0], &changed, FE_ALL);
    /* Unchanged: hand back the (already evaluated) argument itself. */
    if (out && changed) out = fe_cancel_pairs(eval_and_free(out));
    return out;
}

/* ------------------------------------------------------------------ */
/*  Definite-integral finalisation                                     */
/* ------------------------------------------------------------------ */

static size_t fe_count_head(const Expr* e, const char* h) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    size_t c = fe_is_head(e, h) ? 1 : 0;
    c += fe_count_head(e->data.function.head, h);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += fe_count_head(e->data.function.args[i], h);
    return c;
}

static size_t fe_leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    size_t c = fe_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += fe_leaf_count(e->data.function.args[i]);
    return c;
}

/* Refine[e, assum], evaluated and owned. */
static Expr* fe_refine(const Expr* e, const Expr* assum) {
    Expr* call = expr_new_function(expr_new_symbol(SYM_Refine),
                    (Expr*[]){ expr_copy((Expr*)e), expr_copy((Expr*)assum) }, 2);
    return eval_and_free(call);
}

/* A copy of e with the CONDITION of every Piecewise pair Refined on its own
 * under assum (and the Piecewise re-evaluated, dropping False branches and
 * collapsing to a True one). Refine of a whole Piecewise rewrites its
 * conditions but does not decide them -- Abs[b^2/a^2] > 1 became b^2/a^2 > 1
 * under b > a > 0 -- while the standalone condition is decided. Owned. */
static Expr* fe_refine_piecewise(const Expr* e, const Expr* assum) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = fe_refine_piecewise(e->data.function.args[i], assum);
    Expr* node = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (fe_is_head(node, SYM_Piecewise) && node->data.function.arg_count >= 1 &&
        fe_is_head(node->data.function.args[0], SYM_List)) {
        Expr* pairs = node->data.function.args[0];
        for (size_t i = 0; i < pairs->data.function.arg_count; i++) {
            Expr* pr = pairs->data.function.args[i];
            if (!fe_is_head(pr, SYM_List) || pr->data.function.arg_count != 2) continue;
            Expr* rc = fe_refine(pr->data.function.args[1], assum);
            expr_free(pr->data.function.args[1]);
            pr->data.function.args[1] = rc;
        }
    }
    return eval_and_free(node);
}

/* Expand[e], evaluated and owned: collects the like terms a conjugate pair's
 * expansions leave inside separate products. */
static Expr* fe_expand(const Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
                            (Expr*[]){ expr_copy((Expr*)e) }, 1));
}

/* A copy of `e` with every Log[1/n] (integer n >= 2) written -Log[n]; sets
 * *changed when one was rewritten.  The same number, one Log fewer to merge:
 * an endpoint value Log[1 - 1/(1 + x)] at x = 1 otherwise stays apart from the
 * Log[2] of PolyLog[2, 1/2] and PolyLog[3, 1/2]. */
static Expr* fe_log_reciprocal(const Expr* e, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    int64_t n = 0, d = 0;
    if (fe_is_head(e, SYM_Log) && e->data.function.arg_count == 1 &&
        is_rational(e->data.function.args[0], &n, &d) && n == 1 && d >= 2) {
        *changed = true;
        Expr* lg = expr_new_function(expr_new_symbol(SYM_Log),
                        (Expr*[]){ expr_new_integer(d) }, 1);
        return expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), lg }, 2);
    }
    size_t na = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (na ? na : 1));
    for (size_t i = 0; i < na; i++)
        args[i] = fe_log_reciprocal(e->data.function.args[i], changed);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, na);
    free(args);
    return out;
}

/* Expand[e] with every Log of an exact Gaussian rational split into real and
 * imaginary parts (Log[1 + I] -> Log[2]/2 + I Pi/4) and every Log[1/n]
 * written -Log[n]: endpoint values of a ladder in x +- I come in conjugate
 * pairs whose imaginary parts then cancel, and Log[1/2] merges with Log[2].
 * Owned result (never NULL unless evaluation failed). */
static Expr* fe_tidy_logs(const Expr* e) {
    Expr* ec = expr_copy((Expr*)e);
    Expr* g = fe_apply("(#1 /. Log[c_Complex /; Precision[c] === Infinity] :> "
                       "ComplexExpand[Log[c]])&", &ec, 1);
    expr_free(ec);
    bool changed = false;
    Expr* c0 = fe_log_reciprocal(g ? g : e, &changed);
    expr_free(g);
    Expr* c = fe_expand(c0);
    expr_free(c0);
    return c;
}

/* Conjugate pairs of trilogarithms, c (Li_3(z) + Li_3(conj z)) = 2 c Re Li_3(z),
 * closed when z or 1/z lies on the line Re = 1/2 and the third point of
 * Landen's identity is closed (a root of unity):
 *     Li_3(w) + Li_3(1 - w) + Li_3(1 - 1/w)
 *         = Zeta[3] + Log[w]^3/6 + Pi^2/6 Log[w] - Log[w]^2 Log[1 - w]/2,
 * with 1 - w = conj w on that line; for Re 1/z = 1/2 the inversion formula
 *     Li_3(z) = Li_3(1/z) - Pi^2/6 Log[-z] - Log[-z]^3/6
 * first moves z there.  Li_3(1 + I) + Li_3(1 - I) = Pi^2 Log[2]/16 +
 * 35 Zeta[3]/32 -- the pair a ladder in 1 +- I x leaves at x = 1
 * (Int_0^1 Log[1 + x^2]^2/x = Zeta[3]/8).  The real closed form is accepted
 * only when it agrees with the pair to 30 digits.  Owned result: the
 * rewritten r, or a copy of r when no pair closed. */
static Expr* fe_li3_conjugate_pairs(const Expr* r) {
    Expr* rc = expr_copy((Expr*)r);
    Expr* zs = fe_apply("(Union[Cases[#1, PolyLog[3, z_ /; NumericQ[z] && "
                        "Precision[z] === Infinity && TrueQ[Im[z] > 0]] :> z, {0, Infinity}]])&",
                        &rc, 1);
    if (!zs || !fe_is_head(zs, SYM_List) || zs->data.function.arg_count == 0) {
        expr_free(zs); return rc;
    }
    for (size_t i = 0; i < zs->data.function.arg_count; i++) {
        Expr* z = zs->data.function.args[i];
        Expr* cand = fe_apply(
            "(With[{z = #1, e3 = Function[u, Log[u]^3/6 + Pi^2/6 Log[u] - Log[u]^2 Log[1 - u]/2]},"
            " Which[Re[z] == 1/2, Zeta[3] + e3[z] - PolyLog[3, 1 - 1/z],"
            "  Re[1/z] == 1/2, Zeta[3] + e3[1/z] - PolyLog[3, 1 - z]"
            "    - (Pi^2/6 Log[-z] + Log[-z]^3/6) - (Pi^2/6 Log[-Conjugate[z]] + Log[-Conjugate[z]]^3/6),"
            "  True, $Failed]])&", &z, 1);
        if (!cand || cand->type == EXPR_SYMBOL) { expr_free(cand); continue; }
        bool ch = false;
        Expr* red = fe_walk(cand, &ch, FE_ROOTS);
        expr_free(cand);
        Expr* re = fe_apply("(Expand[ComplexExpand[Re[#1 /. Log[c_Complex /; Precision[c] === Infinity] :> "
                            "ComplexExpand[Log[c]]]]])&", &red, 1);
        expr_free(red);
        if (!re) continue;
        Expr* args[3] = { rc, z, re };
        Expr* nr = fe_apply(
            "(With[{c1 = Coefficient[Expand[#1], PolyLog[3, #2]],"
            "       c2 = Coefficient[Expand[#1], PolyLog[3, Conjugate[#2]]]},"
            " If[c1 =!= 0 && FreeQ[c1, PolyLog] && Expand[c1 - c2] === 0 && FreeQ[#3, PolyLog] &&"
            "    TrueQ[Abs[N[#3, 40] - N[PolyLog[3, #2] + PolyLog[3, Conjugate[#2]], 40]] < 10^-30],"
            "  Expand[Expand[#1] - c1 (PolyLog[3, #2] + PolyLog[3, Conjugate[#2]]) + c1 #3], $Failed]])&",
            args, 3);
        expr_free(re);
        if (nr && nr->type != EXPR_SYMBOL) { expr_free(rc); rc = nr; }
        else expr_free(nr);
    }
    expr_free(zs);
    return rc;
}

Expr* fe_reduce_definite_result(Expr* r, const Expr* assumptions) {
    if (!r) return r;
    /* (1) MeijerG -> named functions. Accepted only when it leaves no MeijerG,
     * HypergeometricPFQ or Piecewise behind that the original did not have:
     * a reduction that trades a MeijerG for an unresolved Piecewise (the
     * p == q, m + n <= p case with an undecided |z| vs 1) or for raw pFq's is
     * not a closed form. The Piecewise is resolved by Refine under the call's
     * assumptions when they decide |z| against 1. */
    if (fe_count_head(r, SYM_MeijerG) > 0) {
        bool changed = false;
        Expr* c = fe_walk(r, &changed, FE_MEIJER | FE_PFQ);
        if (changed && assumptions &&
            fe_count_head(c, SYM_Piecewise) > fe_count_head(r, SYM_Piecewise)) {
            Expr* rc = fe_refine_piecewise(c, assumptions);
            expr_free(c); c = rc;
        }
        if (changed &&
            fe_count_head(c, SYM_MeijerG) < fe_count_head(r, SYM_MeijerG) &&
            fe_count_head(c, SYM_HypergeometricPFQ) <= fe_count_head(r, SYM_HypergeometricPFQ) &&
            fe_count_head(c, SYM_Piecewise) <= fe_count_head(r, SYM_Piecewise)) {
            expr_free(r); r = c;
        } else {
            expr_free(c);
        }
    }
    /* (1b) HypergeometricPFQ -> a closed form by Roach's Formula algorithm,
     * accepted when it leaves strictly fewer pFq and no new MeijerG. */
    if (fe_count_head(r, SYM_HypergeometricPFQ) > 0) {
        bool changed = false;
        Expr* c = fe_walk(r, &changed, FE_PFQ);
        if (changed &&
            fe_count_head(c, SYM_HypergeometricPFQ) < fe_count_head(r, SYM_HypergeometricPFQ) &&
            fe_count_head(c, SYM_MeijerG) <= fe_count_head(r, SYM_MeijerG)) {
            expr_free(r); r = c;
        } else {
            expr_free(c);
        }
    }
    /* (2) PolyLog / Log at roots of unity. These are always valid but not
     * always simpler (Log[1 + I] -> Log[2]/2 + I Pi/4), so accepted only when
     * the expanded result is strictly smaller and has introduced no PolyGamma
     * -- i.e. when the parts of conjugate pairs have actually cancelled:
     * PolyLog[3, I] + PolyLog[3, -I] - Zeta[3]/4 -> -7 Zeta[3]/16. */
    if (fe_count_head(r, SYM_PolyLog) + fe_count_head(r, SYM_Log) > 0) {
        bool changed = false;
        Expr* rp = fe_li3_conjugate_pairs(r);
        if (!expr_eq(rp, r)) changed = true;
        Expr* c = fe_walk(rp, &changed, FE_ROOTS);
        expr_free(rp);
        if (changed) { Expr* x = fe_tidy_logs(c); expr_free(c); c = x; }
        /* compared with r tidied the same way, so that the logs a closed
         * dilogarithm brings (Log[1/2], Log[1 + I]) are not counted against it */
        Expr* rt = changed ? fe_tidy_logs(r) : NULL;
        size_t rsz = fe_leaf_count(r);
        if (rt && fe_leaf_count(rt) < rsz) rsz = fe_leaf_count(rt);
        expr_free(rt);
        if (changed && c && fe_leaf_count(c) < rsz &&
            fe_count_head(c, SYM_PolyGamma) <= fe_count_head(r, SYM_PolyGamma)) {
            expr_free(r); r = c;
        } else {
            expr_free(c);
        }
    }
    /* (3) A Newton-Leibniz value assembled from Log and PolyLog endpoint values
     * (Log[1/2] from Log[1 - 1/(1 + x)], PolyLog[3, 1/2] = 7/8 Zeta[3] - ...)
     * comes back as a sum of unexpanded groups.  Log[1/n] -> -Log[n] and one
     * Expand merge the like terms; accepted only when the result is strictly
     * smaller, so a value that is not a sum of such groups keeps its form. */
    if (fe_count_head(r, SYM_Log) > 0) {
        Expr* c = fe_tidy_logs(r);
        if (c && fe_leaf_count(c) < fe_leaf_count(r)) {
            expr_free(r); r = c;
        } else {
            expr_free(c);
        }
    }
    return r;
}

void functionexpand_init(void) {
    symtab_add_builtin("FunctionExpand", builtin_function_expand);
    symtab_get_def("FunctionExpand")->attributes |= ATTR_PROTECTED;
}
