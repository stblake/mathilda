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
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

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

/* A copy of e with every expandable node replaced, bottom-up; sets *changed
 * when anything was rewritten. `which` selects the rule families. Returns an
 * owned Expr. */
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
        if (!s) s = fe_log_root(node);
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
        Expr* c = fe_walk(r, &changed, FE_ROOTS);
        if (changed) { Expr* x = fe_expand(c); expr_free(c); c = x; }
        if (changed && fe_leaf_count(c) < fe_leaf_count(r) &&
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
        bool changed = false;
        Expr* c0 = fe_log_reciprocal(r, &changed);
        Expr* c = fe_expand(c0);
        expr_free(c0);
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
