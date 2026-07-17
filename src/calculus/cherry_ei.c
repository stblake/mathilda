/* cherry_ei.c — Cherry's rational exponential integral: ExpIntegralEi + Erf(i).
 *
 * G. W. Cherry, "An Analysis of the Rational Exponential Integral" (SIAM J.
 * Comput. 18(5), 1989), restricted to the base field C(x): integrands
 * gamma = g E^f with g, f in C(x).  Emits BOTH the exponential-integral (ei) and
 * error-function (erf) terms of the extended Liouville form (Thm 2.3), solved as
 * one linear system.  See CHERRY_PLAN.md §3-4 and risch_special.c.
 *
 * The extended-Liouville form (Thm 2.3) is
 *   INT g E^f dx = y E^f + Sum_i c_i E^(-alpha_i) ei(u_i),   u_i = f + alpha_i,
 * whose derivative gives the rational MATCHING IDENTITY the engine solves:
 *   g = y' + f' y + Sum_i c_i * f'/(f + alpha_i).                          (*)
 * Everything in (*) is rational in x — the exponential kernel E^f divides out —
 * so this is a plain linear system over the constants c_i and the coefficients
 * of y, closed by SolveAlways (the 1989 undetermined-coefficient E1/E2/E3 solve,
 * NOT Risch's differential equation).
 *
 * Argument generation:
 *   P1 (Cherry §4.2): the ei arguments f + alpha_i come from the resultant
 *       h(alpha) = Res_x(g1, p + alpha q)      (f = p/q reduced, q monic;
 *       g1 = the part of den(g) coprime to q).  Its roots are the alpha_i.
 *   P2 (Cherry §4.2, case p in F): a single q-side term u = f itself (alpha = 0)
 *       when f = p/q is a proper fraction with constant numerator.
 * The candidate set is finite (bounded by deg h); over-inclusion is harmless
 * because SolveAlways pins an unused c_i to zero and the diff-back gate rejects
 * any spurious survivor.
 *
 * Scope of THIS engine (first Cherry increment):
 *   - base field C(x) only (a single, x-rational exponential kernel);
 *   - RATIONAL ei-argument constants alpha_i.  Algebraic alpha (Ex 5.3 sqrt2,
 *     p.894 e^x/(x^2-2)) need RootReduce/qqbar and are a later phase (§7); the
 *     engine DECLINES cleanly (NULL) on an algebraic root rather than emit a
 *     partial answer.
 *   - erf terms (perfect-square q) are engine C2, not here.
 * Correct by construction: the emitted answer is accepted only after the exact
 * diff-back gate rt_verify_antideriv, and rejected if it still contains an
 * unevaluated Integrate head (see the driver-side guard below).
 */

#include "cherry_ei.h"
#include "risch_util.h"
#include "risch_singleext.h"

#include "expr.h"
#include "eval.h"
#include "print.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ---- small Expr builders (local, to keep the algorithm readable) ---- */

static Expr* mk_sym(const char* s)        { return expr_new_symbol(s); }
static Expr* mk_int(long n)               { return expr_new_integer(n); }

static Expr* mk_pow(Expr* b, Expr* e) {
    return expr_new_function(mk_sym("Power"), (Expr*[]){ b, e }, 2);
}
static Expr* mk_neg(Expr* a) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ mk_int(-1), a }, 2);
}
static Expr* mk_plus2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Plus"), (Expr*[]){ a, b }, 2);
}
static Expr* mk_times2(Expr* a, Expr* b) {
    return expr_new_function(mk_sym("Times"), (Expr*[]){ a, b }, 2);
}

/* Product of the degree-1 (rational-root) irreducible factors of d, at their
 * multiplicities.  The elementary part y may have a pole only at a RATIONAL root
 * of the integrand; an algebraic-root factor (irreducible over Q of degree >= 2)
 * stays OUT of y's denominator (y is polynomial there).  This matches Cherry's
 * structure AND avoids a shared Together/SolveAlways-over-Q(sqrt d) stack
 * overflow that fires when an algebraic ei pole (x -+ sqrt d) is reduced against
 * a y-denominator carrying the same reducible factor x^2 - d. */
/* Accumulate a single Factor[] factor into *acc if it is a degree-1 (rational
 * root) factor, possibly raised to a power:  base^e with deg_x(base)==1, or a
 * bare linear.  Higher-degree (irreducible over Q) and constant factors are
 * skipped. */
static void accum_linear_factor(Expr* f, Expr* x, Expr** acc) {
    if (!f) return;
    if (f->type == EXPR_FUNCTION && rt_head_is(f, "Power")
        && f->data.function.arg_count == 2) {
        Expr* base = f->data.function.args[0];
        Expr* exn  = f->data.function.args[1];
        if (rt_degree(base, x) == 1)
            *acc = rt_eval_own(mk_times2(*acc, mk_pow(expr_copy(base), expr_copy(exn))));
        return;
    }
    if (rt_degree(f, x) == 1)
        *acc = rt_eval_own(mk_times2(*acc, expr_copy(f)));
}

static Expr* rational_pole_denominator(Expr* d, Expr* x) {
    Expr* fac = rt_eval1("Factor", expr_copy(d));
    Expr* acc = mk_int(1);
    if (fac && fac->type == EXPR_FUNCTION && rt_head_is(fac, "Times")) {
        for (size_t i = 0; i < fac->data.function.arg_count; i++)
            accum_linear_factor(fac->data.function.args[i], x, &acc);
    } else if (fac) {
        accum_linear_factor(fac, x, &acc);
    }
    if (fac) expr_free(fac);
    return acc ? acc : mk_int(1);
}

/* The part of d coprime to q:  strip every factor d shares with q. */
static Expr* poly_coprime_part(Expr* d, Expr* q, Expr* x) {
    Expr* cur = expr_copy(d);
    for (int guard = 0; guard < 64; guard++) {
        Expr* g = rt_eval_call("PolynomialGCD",
            (Expr*[]){ expr_copy(cur), expr_copy(q) }, 2);
        if (!g) break;
        if (rt_degree(g, x) <= 0) { expr_free(g); break; }
        Expr* nxt = rt_eval_call("PolynomialQuotient",
            (Expr*[]){ expr_copy(cur), g, expr_copy(x) }, 3);   /* adopts g */
        if (!nxt) break;
        expr_free(cur);
        cur = nxt;
    }
    return cur;
}

/* Admit a candidate constant iff it generates an extension of degree <= 2 over Q
 * — rational, or a real OR COMPLEX quadratic (sqrt d, +-I, the (1 +- i sqrt3)/2 of
 * Cherry d12).  The complex-conjugate ExpIntegralEi/Erfi diff-back is tractable
 * over a quadratic extension but blows up for higher degree (x^4-2, cube roots),
 * so degree >= 3 is deferred (declines cleanly, never a wrong form — the diff-back
 * still certifies).  A SYMBOLIC constant (sqrt a with a a parameter) is not
 * algebraic over Q, so MinimalPolynomial declines to a non-polynomial form and it
 * is admitted (the diff-back certifies). */
static bool real_admissible(Expr* val) {
    Expr* z = mk_sym("chei$mpz");
    Expr* mp = rt_eval2("MinimalPolynomial", expr_copy(val), expr_copy(z));
    bool accept = true;                       /* symbolic / undecided -> admit */
    if (mp && rt_is_poly(mp, z)) accept = (rt_degree(mp, z) <= 2);
    if (mp) expr_free(mp);
    expr_free(z);
    return accept;
}

/* True iff val has a numerically-nonzero imaginary part (a genuine complex
 * number).  Symbolic values (Im not numerically decidable) -> false. */
static bool numeric_complex(Expr* val) {
    Expr* lt = rt_eval2("Less",
        rt_eval1("Abs", rt_eval1("Im", rt_eval2("N", expr_copy(val), mk_int(30)))),
        mk_pow(mk_int(10), mk_int(-15)));
    bool cplx = lt && lt->type == EXPR_SYMBOL
        && lt->data.symbol.name == intern_symbol("False");
    if (lt) expr_free(lt);
    return cplx;
}

/* True iff e is the $Failed symbol (a declined builtin result). */
static bool is_failed(const Expr* e) {
    return e && e->type == EXPR_SYMBOL
        && e->data.symbol.name == intern_symbol("$Failed");
}

/* Perfect-square beta finder for erf (Cherry 1989 §3).  With q = s^2 (s a
 * polynomial), find the constants beta such that p + beta q = r^2 for a
 * polynomial r; the error-function argument is then u~ = r/s (u~^2 = f + beta).
 * By Thm 3.2 there are at most two.  Candidates: the roots (in b) of the
 * discriminant Res_x(p + b q, d/dx(p + b q)) — where p + b q acquires a repeated
 * factor — together with b = 0 (p itself a square).  Each is validated by
 * PolynomialSqrt plus an exact Expand[r^2 - (p + beta q)] == 0 certificate, and
 * real-gated.  Fills betas[]/rs[] (owned), returns the count. */
static size_t gen_beta_candidates(Expr* p, Expr* q, Expr* x,
                                  Expr** betas, Expr** rs, size_t cap) {
    size_t n = 0;
    Expr* cands[64]; size_t nc = 0;
    cands[nc++] = mk_int(0);                                   /* b = 0 */

    Expr* b = mk_sym("chei$b");
    Expr* pbq = mk_plus2(expr_copy(p), mk_times2(expr_copy(b), expr_copy(q)));
    Expr* disc = rt_eval_call("Resultant",
        (Expr*[]){ expr_copy(pbq), rt_eval2("D", expr_copy(pbq), expr_copy(x)),
                   expr_copy(x) }, 3);
    if (disc && !rt_free_of_x(disc, b)) {
        Expr* eqn = expr_new_function(mk_sym("Equal"),
            (Expr*[]){ expr_copy(disc), mk_int(0) }, 2);
        Expr* sols = rt_eval2("Solve", eqn, expr_copy(b));
        if (sols && sols->type == EXPR_FUNCTION && rt_head_is(sols, "List")) {
            for (size_t i = 0; i < sols->data.function.arg_count && nc < 64; i++) {
                Expr* rule = sols->data.function.args[i];
                if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                    || rule->data.function.arg_count != 1) continue;
                Expr* r = rule->data.function.args[0];
                if (r->type != EXPR_FUNCTION || !rt_head_is(r, "Rule")) continue;
                Expr* val = r->data.function.args[1];
                if (!rt_free_of_x(val, x) || !rt_free_of_x(val, b)) continue;
                cands[nc++] = expr_copy(val);
            }
        }
        if (sols) expr_free(sols);
    }
    if (disc) expr_free(disc);
    expr_free(pbq);
    expr_free(b);

    for (size_t i = 0; i < nc; i++) {
        Expr* bv = cands[i];
        bool dup = false;
        for (size_t k = 0; k < n; k++) if (expr_eq(betas[k], bv)) { dup = true; break; }
        if (dup || !real_admissible(bv) || numeric_complex(bv)) continue; /* erf: real beta */
        Expr* pbqv = rt_eval1("Expand",
            mk_plus2(expr_copy(p), mk_times2(expr_copy(bv), expr_copy(q))));
        Expr* r = rt_eval2("PolynomialSqrt", expr_copy(pbqv), expr_copy(x));
        bool ok = r && !is_failed(r) && rt_is_poly(r, x);
        if (ok) {                                             /* exact certificate */
            Expr* chk = rt_eval1("Expand",
                mk_plus2(mk_pow(expr_copy(r), mk_int(2)), mk_neg(expr_copy(pbqv))));
            ok = chk && rt_is_zero(chk);
            if (chk) expr_free(chk);
        }
        if (ok && n < cap) { betas[n] = expr_copy(bv); rs[n] = expr_copy(r); n++; }
        if (pbqv) expr_free(pbqv);
        if (r) expr_free(r);
    }
    for (size_t i = 0; i < nc; i++) expr_free(cands[i]);
    return n;
}

/* Collect the finite set of candidate ei-argument constants alpha:
 *   P1: roots (in z) of Res_x(g1, p + z q) — RATIONAL or ALGEBRAIC.  Cherry's
 *       theorems assume an algebraically-closed constant field, so sqrt2, +-I,
 *       and Root[] objects are all admissible ei-argument constants (Ex 5.3
 *       sqrt2, p.894 e^x/(x^2-2), d12 complex roots).  The matching-identity
 *       SolveAlways solves the coefficients over the extension the alpha_i
 *       generate, and the diff-back gate rejects any candidate that fails to
 *       verify — so accepting an algebraic root is always sound.
 *   P2: alpha = 0 when q is non-trivial and p is a constant (u = f itself).
 * Fills alphas[] (owned) up to cap, returns the count. */
static size_t gen_alpha_candidates(Expr* g1, Expr* p, Expr* q, Expr* x,
                                   Expr** alphas, size_t cap) {
    size_t n = 0;
    long degq = rt_degree(q, x);
    long degp = rt_degree(p, x);

    /* P1: resultant in the parameter z = chei$z. */
    if (rt_degree(g1, x) >= 1) {
        Expr* z = mk_sym("chei$z");
        Expr* pzq = mk_plus2(expr_copy(p), mk_times2(expr_copy(z), expr_copy(q)));
        Expr* res = rt_eval_call("Resultant",
            (Expr*[]){ expr_copy(g1), pzq, expr_copy(x) }, 3);   /* adopts pzq */
        if (res && !rt_free_of_x(res, z)) {
            Expr* eqn = expr_new_function(mk_sym("Equal"),
                (Expr*[]){ expr_copy(res), mk_int(0) }, 2);
            Expr* sols = rt_eval2("Solve", eqn, expr_copy(z));
            if (sols && sols->type == EXPR_FUNCTION
                && rt_head_is(sols, "List")) {
                for (size_t i = 0; i < sols->data.function.arg_count; i++) {
                    Expr* rule = sols->data.function.args[i];   /* {z -> val} */
                    if (rule->type != EXPR_FUNCTION || !rt_head_is(rule, "List")
                        || rule->data.function.arg_count != 1) continue;
                    Expr* r = rule->data.function.args[0];
                    if (r->type != EXPR_FUNCTION || !rt_head_is(r, "Rule")) continue;
                    Expr* val = r->data.function.args[1];
                    /* a genuine constant root (free of x and of the parameter z) */
                    if (!rt_free_of_x(val, x) || !rt_free_of_x(val, z)) continue;
                    /* Degree gate (<= 2): real or complex-quadratic. */
                    if (!real_admissible(val)) continue;
                    /* A complex root is tractable only as the SOLE conjugate pair:
                     * f a polynomial (q = 1, so no P2 term) and g1 an irreducible
                     * quadratic (deg 2).  Any mixing of complex ei poles with other
                     * poles (a P2 term, a rational factor) blows up Together over
                     * Q(i) — defer those. */
                    if (numeric_complex(val) && (degq != 0 || rt_degree(g1, x) != 2))
                        continue;
                    /* dedup */
                    bool dup = false;
                    for (size_t k = 0; k < n; k++)
                        if (expr_eq(alphas[k], val)) { dup = true; break; }
                    if (!dup && n < cap) alphas[n++] = expr_copy(val);
                }
            }
            if (sols) expr_free(sols);
        }
        if (res) expr_free(res);
        expr_free(z);
    }

    /* P2 (case p in F): proper fraction with constant numerator -> u = f, alpha = 0. */
    if (degq >= 1 && degp == 0) {
        bool dup = false;
        for (size_t k = 0; k < n; k++)
            if (rt_is_zero(alphas[k])) { dup = true; break; }
        if (!dup && n < cap) alphas[n++] = mk_int(0);
    }
    return n;
}

/* A1 number-field fallback for a LONE complex-conjugate ei pair (Cherry C = C-bar).
 * a0, a1 are conjugate complex constants (roots of the P1 resultant of an irreducible
 * quadratic pole).  Mathilda's Solve/Together do not reliably solve the coefficient
 * system over Q(i sqrt d): Together fails to cancel the two complex-linear factors
 * (x + a_i) back against the real quadratic (x + a0)(x + a1), so the residual is
 * over-determined and Solve then declines.  This solves the system over Q instead, by
 * the number-field basis method: write a = center +- chs with chs the algebraic
 * generator (chs^2 = disc, disc = ((a0-a1)/2)^2 a negative rational), carry the
 * conjugate coefficients as c0 = cza + czb*chs, c1 = cza - czb*chs and y with rational
 * unknowns, reduce the residual modulo chs^2 - disc, split it into the {1, chs} Q-basis
 * (both parts must vanish), and Solve over Q — where Solve is reliable.  The generator
 * is then back-substituted (chs -> (a0-a1)/2) and the answer diff-back verified.
 * Borrows all arguments; returns a fresh, verified antiderivative or NULL. */
static Expr* rt_cherry_ei_conjpair_nf(Expr* f, Expr* v, Expr* g, Expr* a0, Expr* a1,
                                      Expr* sden, long Ny, Expr* x) {
    Expr* rhalf = expr_new_function(mk_sym("Rational"), (Expr*[]){ mk_int(1), mk_int(2) }, 2);
    /* center = (a0+a1)/2 is the RATIONAL real part; Simplify (not Together) is needed to
     * cancel the conjugate imaginary parts (Together leaves e.g. 1/2(1/2(2-2I√2)+...)). */
    Expr* center = rt_eval1("Simplify",
        mk_times2(mk_plus2(expr_copy(a0), expr_copy(a1)), expr_copy(rhalf)));   /* (a0+a1)/2 */
    Expr* hd = rt_eval1("Together",
        mk_times2(mk_plus2(expr_copy(a0), mk_neg(expr_copy(a1))), rhalf));      /* (a0-a1)/2 = chs value */
    Expr* disc = rt_eval1("Simplify", mk_pow(expr_copy(hd), mk_int(2)));        /* hd^2 (neg rational) */
    if (!center || !hd || !disc || !rt_free_of_x(center, x) || !rt_free_of_x(disc, x)) {
        if (center) expr_free(center); if (hd) expr_free(hd); if (disc) expr_free(disc);
        return NULL;
    }
    Expr* chs = mk_sym("chei$s");
    long ny = Ny < 0 ? 0 : Ny;

    /* unknowns: y_0..y_ny (rational), cza, czb */
    size_t nun = (size_t)(ny + 1) + 2;
    Expr** unk = malloc(nun * sizeof(Expr*));
    size_t ui = 0;
    Expr** yterms = malloc((size_t)(ny + 1) * sizeof(Expr*));
    for (long j = 0; j <= ny; j++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chei$y%ld", j);
        Expr* yj = mk_sym(nm);
        unk[ui++] = expr_copy(yj);
        yterms[j] = mk_times2(yj, mk_pow(expr_copy(x), mk_int(j)));
    }
    Expr* yans = mk_times2(expr_new_function(mk_sym("Plus"), yterms, (size_t)(ny + 1)),
                           mk_pow(expr_copy(sden), mk_int(-1)));
    free(yterms);
    Expr* cza = mk_sym("chei$cza"); Expr* czb = mk_sym("chei$czb");
    unk[ui++] = expr_copy(cza); unk[ui++] = expr_copy(czb);
    Expr* c0e = mk_plus2(expr_copy(cza), mk_times2(expr_copy(czb), expr_copy(chs)));       /* cza + czb chs */
    Expr* c1e = mk_plus2(expr_copy(cza), mk_neg(mk_times2(expr_copy(czb), expr_copy(chs))));/* cza - czb chs */

    /* rhs = y' + v' y + c0 v'/(v + (center+chs)) + c1 v'/(v + (center-chs)) */
    Expr* vderiv = rt_eval2("D", expr_copy(v), expr_copy(x));
    Expr* yder = rt_eval2("D", expr_copy(yans), expr_copy(x));
    Expr* a0s = mk_plus2(expr_copy(center), expr_copy(chs));
    Expr* a1s = mk_plus2(expr_copy(center), mk_neg(expr_copy(chs)));
    Expr* t1 = mk_times2(expr_copy(vderiv), expr_copy(yans));
    Expr* t2 = mk_times2(c0e, mk_times2(expr_copy(vderiv),
                    mk_pow(mk_plus2(expr_copy(v), a0s), mk_int(-1))));   /* adopts c0e, a0s */
    Expr* t3 = mk_times2(c1e, mk_times2(expr_copy(vderiv),
                    mk_pow(mk_plus2(expr_copy(v), a1s), mk_int(-1))));   /* adopts c1e, a1s */
    expr_free(vderiv);
    Expr* rhs = expr_new_function(mk_sym("Plus"), (Expr*[]){ yder, t1, t2, t3 }, 4);
    Expr* rnum = rt_eval1("Numerator", rt_eval1("Together",
        mk_plus2(expr_copy(g), mk_neg(rhs))));

    Expr* result = NULL;
    if (rnum) {
        Expr* minp = mk_plus2(mk_pow(expr_copy(chs), mk_int(2)), mk_neg(expr_copy(disc))); /* chs^2 - disc */
        Expr* rd = rt_eval_call("PolynomialRemainder",
            (Expr*[]){ rt_eval1("Expand", rnum), minp, expr_copy(chs) }, 3);   /* adopts rnum, minp */
        if (rd) {
            Expr* A = rt_eval_call("Coefficient",
                (Expr*[]){ expr_copy(rd), expr_copy(chs), mk_int(0) }, 3);     /* chs^0 part */
            Expr* B = rt_eval_call("Coefficient",
                (Expr*[]){ rd, expr_copy(chs), mk_int(1) }, 3);                /* adopts rd; chs^1 part */
            Expr* clA = A ? rt_eval2("CoefficientList", A, expr_copy(x)) : NULL;
            Expr* clB = B ? rt_eval2("CoefficientList", B, expr_copy(x)) : NULL;
            Expr* allc = (clA && clB) ? rt_eval2("Join", clA, clB) : NULL;     /* adopts clA, clB */
            if (!allc) { if (clA) expr_free(clA); if (clB) expr_free(clB); }
            Expr* eqs = allc ? rt_eval1("Thread", expr_new_function(mk_sym("Equal"),
                (Expr*[]){ allc, mk_int(0) }, 2)) : NULL;                      /* adopts allc */
            Expr** ul = malloc(nun * sizeof(Expr*));
            for (size_t jj = 0; jj < nun; jj++) ul[jj] = expr_copy(unk[jj]);
            Expr* unklist = expr_new_function(mk_sym("List"), ul, nun);
            free(ul);
            Expr* sol = eqs ? rt_eval2("Solve", eqs, unklist) : (expr_free(unklist), (Expr*)NULL);
            bool solved = sol && sol->type == EXPR_FUNCTION && rt_head_is(sol, "List")
                && sol->data.function.arg_count >= 1
                && sol->data.function.args[0]->type == EXPR_FUNCTION
                && rt_head_is(sol->data.function.args[0], "List");
            if (solved) {
                Expr* rules = sol->data.function.args[0];
                /* apply the solution, then back-substitute chs -> hd (the real value) */
                Expr* chs2hd = expr_new_function(mk_sym("List"),
                    (Expr*[]){ expr_new_function(mk_sym("Rule"),
                        (Expr*[]){ expr_copy(chs), expr_copy(hd) }, 2) }, 1);
                Expr* c0v = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                    (Expr*[]){ rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                        (Expr*[]){ mk_plus2(expr_copy(cza), mk_times2(expr_copy(czb), expr_copy(chs))),
                                   expr_copy(rules) }, 2)), expr_copy(chs2hd) }, 2));
                Expr* c1v = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                    (Expr*[]){ rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                        (Expr*[]){ mk_plus2(expr_copy(cza), mk_neg(mk_times2(expr_copy(czb), expr_copy(chs)))),
                                   expr_copy(rules) }, 2)), expr_copy(chs2hd) }, 2));
                Expr* yv = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                    (Expr*[]){ expr_copy(yans), expr_copy(rules) }, 2));
                expr_free(chs2hd);
                /* answer = y E^v + c0 E^(-a0) ei(v+a0) + c1 E^(-a1) ei(v+a1) */
                Expr* ei0 = expr_new_function(mk_sym("ExpIntegralEi"),
                    (Expr*[]){ mk_plus2(expr_copy(v), expr_copy(a0)) }, 1);
                Expr* ei1 = expr_new_function(mk_sym("ExpIntegralEi"),
                    (Expr*[]){ mk_plus2(expr_copy(v), expr_copy(a1)) }, 1);
                Expr* term0 = mk_times2(yv, mk_pow(mk_sym("E"), expr_copy(v)));
                Expr* term1 = mk_times2(c0v, mk_times2(mk_pow(mk_sym("E"), mk_neg(expr_copy(a0))), ei0));
                Expr* term2 = mk_times2(c1v, mk_times2(mk_pow(mk_sym("E"), mk_neg(expr_copy(a1))), ei1));
                Expr* Q = rt_eval_own(expr_new_function(mk_sym("Plus"),
                    (Expr*[]){ term0, term1, term2 }, 3));
                if (Q && rt_free_of_head(Q, "Integrate") && rt_verify_antideriv(Q, f, x))
                    result = Q;
                else if (Q) expr_free(Q);
            }
            if (sol) expr_free(sol);
        }
    }

    for (size_t jj = 0; jj < nun; jj++) expr_free(unk[jj]);
    free(unk);
    expr_free(yans); expr_free(cza); expr_free(czb); expr_free(chs);
    expr_free(center); expr_free(hd); expr_free(disc);
    return result;
}

Expr* rt_cherry_ei(Expr* f, Expr* x) {
    /* 1. The exponential kernel E^v (v = the Cherry exponent f).  Base field:
     *    v must be a rational function of x alone (no nested exp/log). */
    Expr* v = rt_find_exp_of_x(f, x);          /* borrowed */
    if (!v || !rt_kernel_simple(v, x)) return NULL;

    /* 1a. Constant-offset reduction.  If the exponent v has a nonzero x-free
     *     polynomial part c (v = c + proper-fraction), factor E^v = E^c E^(v-c): the
     *     constant E^c folds into the cofactor.  Without this, the constant term
     *     inflates deg(p) and defeats the P2 (deg p == 0) recognition — e.g.
     *     E^(1/x+2) = E^2 E^(1/x) would decline.  Recurse on f E^(-c) (whose exponent
     *     v-c has no constant term) and scale by E^c.  Fires ONLY when the polynomial
     *     part of v is a nonzero CONSTANT (degree 0), so a genuine polynomial exponent
     *     part (E^(x+1), E^((x^4+a)/x^2)) is untouched. */
    {
        Expr* vt = rt_eval1("Together", expr_copy(v));
        Expr* pv = vt ? rt_eval1("Numerator", expr_copy(vt)) : NULL;
        Expr* qv = vt ? rt_eval1("Denominator", expr_copy(vt)) : NULL;
        if (vt) expr_free(vt);
        Expr* vc = NULL;
        if (pv && qv && rt_is_poly(pv, x) && rt_is_poly(qv, x)) {
            Expr* quo = rt_eval_call("PolynomialQuotient",
                (Expr*[]){ expr_copy(pv), expr_copy(qv), expr_copy(x) }, 3);
            if (quo && rt_is_poly(quo, x) && rt_degree(quo, x) == 0 && !rt_is_zero(quo))
                vc = expr_copy(quo);
            if (quo) expr_free(quo);
        }
        if (pv) expr_free(pv);
        if (qv) expr_free(qv);
        if (vc) {
            Expr* emc = mk_pow(mk_sym("E"), mk_neg(expr_copy(vc)));
            Expr* fstr = rt_eval1("Together", mk_times2(expr_copy(f), emc));  /* f E^(-c) */
            Expr* sub = fstr ? rt_cherry_ei(fstr, x) : NULL;
            if (fstr) expr_free(fstr);
            if (sub) {
                Expr* ec = mk_pow(mk_sym("E"), vc);           /* adopts vc */
                return rt_eval_own(mk_times2(ec, sub));       /* E^c * INT f E^(-c) */
            }
            expr_free(vc);   /* recursion declined — fall through to the normal path */
        }
    }

    /* g = Together[f E^(-v)] — the rational cofactor after peeling E^v. */
    Expr* emv = mk_pow(mk_sym("E"), mk_neg(expr_copy(v)));
    Expr* g = rt_eval1("Together", mk_times2(expr_copy(f), emv));
    if (!g) return NULL;
    if (rt_find_exp_of_x(g, x) || rt_find_log_of_x(g, x)) { expr_free(g); return NULL; }

    Expr* gnum = rt_eval1("Numerator", expr_copy(g));
    Expr* gden = rt_eval1("Denominator", expr_copy(g));
    if (!gnum || !gden || !rt_is_poly(gnum, x) || !rt_is_poly(gden, x)) {
        if (gnum) expr_free(gnum);
        if (gden) expr_free(gden);
        expr_free(g);
        return NULL;
    }

    /* 2. v = p/q reduced; make q monic (so P1/P2 read leading coefficients cleanly). */
    Expr* vt = rt_eval1("Together", expr_copy(v));
    Expr* p  = rt_eval1("Numerator", expr_copy(vt));
    Expr* q  = rt_eval1("Denominator", expr_copy(vt));
    expr_free(vt);
    if (!p || !q || !rt_is_poly(p, x) || !rt_is_poly(q, x)) {
        if (p) expr_free(p);
        if (q) expr_free(q);
        expr_free(gnum); expr_free(gden); expr_free(g);
        return NULL;
    }
    long degq = rt_degree(q, x);
    if (degq >= 1) {
        Expr* lcq = rt_coeff(q, x, degq);
        Expr* invlc = mk_pow(lcq, mk_int(-1));                 /* adopts lcq */
        Expr* qn = rt_eval1("Expand", mk_times2(q, expr_copy(invlc)));   /* adopts q */
        Expr* pn = rt_eval1("Expand", mk_times2(p, invlc));    /* adopts p, invlc */
        q = qn ? qn : q; p = pn ? pn : p;
    }

    /* 3. g1 = the factor of den(g) coprime to q (the ei-pole carrier for P1). */
    Expr* g1 = poly_coprime_part(gden, q, x);

    /* 4. Generate the finite candidate ei-argument constants alpha_i (rational
     *    or real-algebraic — Cherry's algebraically-closed constant field). */
    Expr* alphas[64];
    size_t m = gen_alpha_candidates(g1, p, q, x, alphas, 64);
    expr_free(g1);

    /* 4b. erf candidates: when q = s^2 is a perfect square (s a polynomial),
     *     completing the square yields <=2 error-function arguments u~_j = r_j/s
     *     with r_j^2 = p + beta_j q (Cherry 1989 §3, Thm 3.2). */
    Expr* s = rt_eval2("PolynomialSqrt", expr_copy(q), expr_copy(x));
    bool q_square = s && !is_failed(s) && rt_is_poly(s, x);
    Expr* betas[8]; Expr* rs[8]; size_t me = 0;
    if (q_square) me = gen_beta_candidates(p, q, x, betas, rs, 8);
    if (s && !q_square) { expr_free(s); s = NULL; }

    /* No ei AND no erf candidate -> not a special-function problem: a pure
     * elementary exponential integral belongs to the transcendental case that
     * already ran ahead of the special layer.  Decline (and avoid a spurious
     * y-only solve, slow on a non-elementary integrand like E^x/(x^2+1) once its
     * complex roots are gated out). */
    if (m == 0 && me == 0) {
        for (size_t i = 0; i < m; i++) expr_free(alphas[i]);
        if (s) expr_free(s);
        expr_free(gnum); expr_free(gden); expr_free(g);
        expr_free(p); expr_free(q);
        return NULL;
    }

    /* 5. Build the y-ansatz  y = Y(x)/s,  s = den(g) * q  at FULL multiplicity
     *    (NOT its radical): the elementary part carries a pole of multiplicity
     *    m-1 where g has multiplicity m (the ei/log-derivative term removes one
     *    power), so den(y) is generally NOT squarefree — e.g.
     *    INT e^x/(x+1)^3 = (1/2) e^-1 ei(x+1) - e^x (x+2)/(2 (x+1)^2), den(y) =
     *    (x+1)^2.  Using the full denominator (times q for the exponent's poles)
     *    covers den(y) | den(g) q; the excess powers are pinned to zero by
     *    SolveAlways, and the diff-back gate guarantees soundness regardless.
     *    A generous DERIVED degree bound on the numerator Y. */
    Expr* dq = rt_eval1("Expand", mk_times2(expr_copy(gden), expr_copy(q)));
    Expr* sden = rational_pole_denominator(dq, x);
    /* the Y-degree bound must still cover a polynomial elementary part over the
     * FULL denominator degree (algebraic poles push their contribution into Y). */
    long Ny = rt_degree(gnum, x) + rt_degree(dq, x) + degq + 2;
    expr_free(dq);
    if (Ny < 0) Ny = 0;

    /* A1 (complex C = C-bar): when an ei-argument constant is COMPLEX (a lone
     * conjugate pair over Q(i)/Q(i sqrt3) — e.g. d12 (x^2+1)e^x/(x^2+x+1) over
     * Q(i sqrt3)), the coefficient system is solved over that number field.  The
     * generous real-case bound above inflates Y to ~6 unknowns, and Solve over the
     * algebraic field does NOT reduce that larger mixed system (whereas the tight
     * system solves cleanly — verified).  For an admitted complex candidate the
     * only pole is the irreducible quadratic g1 (the numeric-complex gate requires
     * q = 1, deg g1 = 2, so sden has no rational factor), hence y is a pure
     * polynomial whose degree is the polynomial part of g plus a small margin.
     * Tighten Ny to keep the number-field Solve tractable; the exact diff-back gate
     * still guarantees soundness, so a too-small bound can only decline. */
    bool has_complex = false;
    for (size_t i = 0; i < m; i++)
        if (numeric_complex(alphas[i])) { has_complex = true; break; }
    if (has_complex) {
        long poly_part = rt_degree(gnum, x) - rt_degree(gden, x);
        if (poly_part < 0) poly_part = 0;
        long Ny_tight = poly_part + rt_degree(sden, x) + 1;
        if (Ny_tight >= 0 && Ny_tight < Ny) Ny = Ny_tight;
    }

    size_t nsym = (size_t)(Ny + 1) + m + me;
    Expr** syms = malloc((nsym ? nsym : 1) * sizeof(Expr*));
    size_t si = 0;

    /* Y = sum_{j=0}^{Ny} a_j x^j */
    Expr** yterms = malloc((size_t)(Ny + 1) * sizeof(Expr*));
    for (long j = 0; j <= Ny; j++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chei$a%ld", j);
        Expr* a = mk_sym(nm);
        syms[si++] = expr_copy(a);
        yterms[j] = mk_times2(a, mk_pow(expr_copy(x), mk_int(j)));
    }
    Expr* Ypoly = expr_new_function(mk_sym("Plus"), yterms, (size_t)(Ny + 1));
    free(yterms);
    Expr* yansatz = mk_times2(Ypoly, mk_pow(expr_copy(sden), mk_int(-1)));

    Expr* vderiv = rt_eval2("D", expr_copy(v), expr_copy(x));
    Expr* yderiv = rt_eval2("D", expr_copy(yansatz), expr_copy(x));

    /* rhs = y' + v' y + sum_i c_i v'/(v+alpha_i)
     *                 + sum_j k_j (2/Sqrt[Pi]) (r_j/s)'
     * (the e^(-beta) answer prefactor cancels e^(u~^2)=e^(f+beta), so the erf
     *  contribution to the divided-by-e^f identity stays rational in x.) */
    size_t nrhs = 2 + m + me;
    Expr** rhs = malloc(nrhs * sizeof(Expr*));
    rhs[0] = yderiv;
    rhs[1] = mk_times2(expr_copy(vderiv), expr_copy(yansatz));
    Expr** csyms = malloc((m ? m : 1) * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chei$c%zu", i);
        Expr* c = mk_sym(nm);
        csyms[i] = expr_copy(c);
        syms[si++] = expr_copy(c);
        Expr* vpa = mk_plus2(expr_copy(v), expr_copy(alphas[i]));
        rhs[2 + i] = mk_times2(c, mk_times2(expr_copy(vderiv),
                                            mk_pow(vpa, mk_int(-1))));
    }
    expr_free(vderiv);
    Expr** ksyms = malloc((me ? me : 1) * sizeof(Expr*));
    for (size_t j = 0; j < me; j++) {
        char nm[32]; snprintf(nm, sizeof(nm), "chei$k%zu", j);
        Expr* k = mk_sym(nm);
        ksyms[j] = expr_copy(k);
        syms[si++] = expr_copy(k);
        Expr* us = rt_eval1("Together",
            mk_times2(expr_copy(rs[j]), mk_pow(expr_copy(s), mk_int(-1))));  /* r_j/s */
        Expr* dus = rt_eval2("D", us, expr_copy(x));                         /* adopts us */
        Expr* twopi = mk_times2(mk_int(2), mk_pow(mk_sym("Pi"),
            expr_new_function(mk_sym("Rational"), (Expr*[]){ mk_int(-1), mk_int(2) }, 2)));
        rhs[2 + m + j] = mk_times2(k, mk_times2(twopi, dus));
    }
    Expr* rhs_sum = expr_new_function(mk_sym("Plus"), rhs, nrhs);
    free(rhs);

    /* resid = g - rhs;  numerator of Together must vanish identically in x. */
    Expr* resid = mk_plus2(expr_copy(g), mk_neg(rhs_sum));
    Expr* rnum = rt_eval1("Numerator", rt_eval1("Together", resid));
    Expr* sol = NULL;
    if (rnum) {
        /* Every x-coefficient of the numerator must vanish.  Solve for OUR
         * unknowns EXPLICITLY (the ansatz syms) rather than via SolveAlways —
         * SolveAlways would also try to solve for any symbolic PARAMETER in the
         * integrand (e.g. the `a` in E^((x^4+a)/x^2)), which corrupts the system. */
        Expr* clist = rt_eval2("CoefficientList", rnum, expr_copy(x));   /* adopts rnum */
        Expr* eqs = rt_eval1("Thread", expr_new_function(mk_sym("Equal"),
            (Expr*[]){ clist, mk_int(0) }, 2));                          /* {c_i == 0} */
        Expr** ul = malloc((nsym ? nsym : 1) * sizeof(Expr*));
        for (size_t jj = 0; jj < nsym; jj++) ul[jj] = expr_copy(syms[jj]);
        Expr* unklist = expr_new_function(mk_sym("List"), ul, nsym);
        free(ul);
        sol = rt_eval2("Solve", eqs, unklist);
    }

    /* 6. Assemble the answer from the solution and diff-back verify. */
    Expr* result = NULL;
    bool solved = sol && sol->type == EXPR_FUNCTION && rt_head_is(sol, "List")
        && sol->data.function.arg_count >= 1
        && sol->data.function.args[0]->type == EXPR_FUNCTION
        && rt_head_is(sol->data.function.args[0], "List");
    if (solved) {
        Expr* rules = sol->data.function.args[0];

        /* answer = y E^v + sum_i c_i E^(-alpha_i) ExpIntegralEi[v+alpha_i]
         *                + sum_j k_j E^(-beta_j) Erfi[r_j/s]                */
        size_t nans = 1 + m + me;
        Expr** ans = malloc(nans * sizeof(Expr*));
        ans[0] = mk_times2(expr_copy(yansatz), mk_pow(mk_sym("E"), expr_copy(v)));
        for (size_t i = 0; i < m; i++) {
            Expr* ei = expr_new_function(mk_sym("ExpIntegralEi"),
                (Expr*[]){ mk_plus2(expr_copy(v), expr_copy(alphas[i])) }, 1);
            Expr* wt = mk_pow(mk_sym("E"), mk_neg(expr_copy(alphas[i])));
            ans[1 + i] = mk_times2(expr_copy(csyms[i]), mk_times2(wt, ei));
        }
        for (size_t j = 0; j < me; j++) {
            Expr* us = rt_eval1("Together",
                mk_times2(expr_copy(rs[j]), mk_pow(expr_copy(s), mk_int(-1))));
            Expr* erfi = expr_new_function(mk_sym("Erfi"), (Expr*[]){ us }, 1);  /* adopts us */
            Expr* wt = mk_pow(mk_sym("E"), mk_neg(expr_copy(betas[j])));
            ans[1 + m + j] = mk_times2(expr_copy(ksyms[j]), mk_times2(wt, erfi));
        }
        Expr* Q = expr_new_function(mk_sym("Plus"), ans, nans);
        free(ans);

        /* substitute the solution, then pin any free unknown to 0 */
        Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
            (Expr*[]){ Q, expr_copy(rules) }, 2));
        if (Q) {
            Expr** zero = malloc((nsym ? nsym : 1) * sizeof(Expr*));
            for (size_t jj = 0; jj < nsym; jj++)
                zero[jj] = expr_new_function(mk_sym("Rule"),
                    (Expr*[]){ expr_copy(syms[jj]), mk_int(0) }, 2);
            Expr* zl = expr_new_function(mk_sym("List"), zero, nsym);
            free(zero);
            Q = rt_eval_own(expr_new_function(mk_sym("ReplaceAll"),
                (Expr*[]){ Q, zl }, 2));
        }

        /* Driver-side emission gate (CHERRY_PLAN.md §1.3 hazard 2 / A4 F5): a
         * Cherry answer is a COMPLETE antiderivative — reject any candidate that
         * still hides an unevaluated Integrate head (which the FTC rule would
         * otherwise let slip through the diff-back), then the exact gate. */
        if (Q && rt_free_of_head(Q, "Integrate") && rt_verify_antideriv(Q, f, x))
            result = Q;
        else if (Q) expr_free(Q);
    }
    if (sol) expr_free(sol);

    /* A1 number-field fallback: a lone complex-conjugate ei pair (m == 2, no erf)
     * that the primary Solve over Q(i sqrt d) could not crack — solve it over Q by
     * the {1, chs} basis method instead (rt_cherry_ei_conjpair_nf).  Fires ONLY when
     * the direct path already failed, so every case the direct Solve closes stays
     * byte-identical; the fallback newly closes e.g. E^x/(x^2+2x+3) (Q(i sqrt2)). */
    if (!result && m == 2 && me == 0
        && numeric_complex(alphas[0]) && numeric_complex(alphas[1]))
        result = rt_cherry_ei_conjpair_nf(f, v, g, alphas[0], alphas[1], sden, Ny, x);

    for (size_t i = 0; i < m; i++) { expr_free(alphas[i]); expr_free(csyms[i]); }
    free(csyms);
    for (size_t j = 0; j < me; j++) { expr_free(betas[j]); expr_free(rs[j]); expr_free(ksyms[j]); }
    free(ksyms);
    for (size_t j = 0; j < nsym; j++) expr_free(syms[j]);
    free(syms);
    if (s) expr_free(s);
    expr_free(yansatz); expr_free(sden);
    expr_free(gnum); expr_free(gden); expr_free(g);
    expr_free(p); expr_free(q);
    return result;
}

/* ================================================================== */
/* Multi-term exponential (Cherry Thm 5.4 case b, flat level).         */
/* ================================================================== */

/* f rational in a SINGLE exponential kernel E^w (w rational in x) whose Laurent
 * expansion has several essential terms  Sum_i p_i(x) E^(i w).  The single-shape
 * rt_cherry_ei cannot peel these together — its cofactor g = f E^(-v) must be free
 * of exp, which fails as soon as two commensurate exponentials (E^w, E^(2w), …)
 * coexist (e.g. (E^x + E^(2x))/(x-1)).  Cherry's Thm 5.4 case b integrates the term
 * by term: Laurent-split in t = E^w, then integrate each essential term
 * p_i E^(i w) with rt_cherry_ei (recovering its ei/erf part) and the t^0 term with
 * the rational integrator, and sum.  Reached only after rt_cherry_ei declined, so a
 * single-term integrand never lands here.  Diff-back verified as a whole; a mis-split
 * can only decline. */
Expr* rt_cherry_exp_multiterm(Expr* f, Expr* x) {
    Expr* v = rt_find_exp_of_x(f, x);
    if (!v || !rt_kernel_simple(v, x)) return NULL;

    /* Kernelize to t = rmT = E^uexp (commensurate exponents E^(k w) -> t^k). */
    Expr* uexp = NULL;
    Expr* F = rt_exp_kernelize(f, x, &uexp);
    if (!F || !uexp) { if (F) expr_free(F); if (uexp) expr_free(uexp); return NULL; }
    Expr* t = mk_sym("rmT");

    Expr* G   = rt_eval1("Together", expr_copy(F));
    Expr* num = G ? rt_eval1("Numerator", expr_copy(G)) : NULL;
    Expr* den = G ? rt_eval1("Denominator", expr_copy(G)) : NULL;
    Expr* result = NULL;

    /* Require a PURE Laurent form: den a monomial c*t^a in t (a t-dependent
     * denominator pole is the Hermite/log hyperexponential case handled elsewhere),
     * num a polynomial in t, and no residual exp/log after kernelization. */
    bool shape = num && den && rt_is_poly(num, t) && rt_is_poly(den, t)
        && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL;
    long a = -1; Expr* lc = NULL;
    if (shape) {
        Expr* cl = rt_eval2("CoefficientList", expr_copy(den), expr_copy(t));
        if (cl && cl->type == EXPR_FUNCTION && rt_head_is(cl, "List")) {
            size_t nz = 0;
            for (size_t i = 0; i < cl->data.function.arg_count; i++)
                if (!rt_is_zero(cl->data.function.args[i])) {
                    if (a < 0) a = (long)i;
                    nz++;
                }
            if (nz == 1 && a >= 0) lc = expr_copy(cl->data.function.args[a]);  /* monomial */
        }
        if (cl) expr_free(cl);
    }

    if (lc) {
        long dnum = rt_degree(num, t);
        Expr** pieces = malloc((size_t)(dnum + 1) * sizeof(Expr*));
        size_t np = 0;
        bool okall = true;
        bool any_ei = false;   /* require >=2 essential terms (else rt_cherry_ei would have fired) */
        int nterms = 0;
        for (long j = 0; j <= dnum && okall; j++) {
            Expr* c = rt_coeff(num, t, j);
            if (rt_is_zero(c)) { expr_free(c); continue; }
            long ip = j - a;
            Expr* pi = rt_eval1("Cancel", mk_times2(c, mk_pow(expr_copy(lc), mk_int(-1)))); /* p_i(x) */
            nterms++;
            Expr* piece = NULL;
            if (ip == 0) {
                /* the t^0 term integrates in C(x) */
                Expr* r = rt_eval_call("Integrate`BronsteinRational",
                    (Expr*[]){ expr_copy(pi), expr_copy(x) }, 2);
                if (r && !rt_head_is(r, "Integrate`BronsteinRational")) piece = r;
                else if (r) expr_free(r);
            } else {
                /* the essential term p_i E^(ip uexp) -> ei/erf via rt_cherry_ei */
                Expr* sub = mk_times2(expr_copy(pi),
                    mk_pow(mk_sym("E"), mk_times2(mk_int(ip), expr_copy(uexp))));
                Expr* sub_e = rt_eval1("Together", sub);   /* normalize E^(ip w) */
                piece = sub_e ? rt_cherry_ei(sub_e, x) : NULL;
                if (sub_e) expr_free(sub_e);
                if (piece) any_ei = true;
            }
            expr_free(pi);
            if (!piece) { okall = false; }
            else pieces[np++] = piece;
        }
        /* Only emit a genuine MULTI-term answer with at least one ei/erf piece — a
         * single essential term is rt_cherry_ei's job (already tried and declined for
         * a reason we must not paper over). */
        if (okall && any_ei && nterms >= 2 && np > 0) {
            Expr* sum = expr_new_function(mk_sym("Plus"), pieces, np);
            Expr* cand = rt_eval_own(sum);
            if (cand && rt_free_of_head(cand, "Integrate") && rt_verify_antideriv(cand, f, x))
                result = cand;
            else if (cand) expr_free(cand);
        } else {
            for (size_t i = 0; i < np; i++) expr_free(pieces[i]);
        }
        free(pieces);
        expr_free(lc);
    }

    if (G) expr_free(G);
    if (num) expr_free(num);
    if (den) expr_free(den);
    expr_free(F); expr_free(uexp); expr_free(t);
    return result;
}
