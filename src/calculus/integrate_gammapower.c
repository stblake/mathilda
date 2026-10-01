/* integrate_gammapower.c — symbolic power times an exponential -> incomplete Gamma.
 *
 * Recognizes  K x^p E^(a x^m)  where the exponent p is SYMBOLIC (carries a symbol
 * rather than being a number), m a positive integer, and K, a free of x, and
 * integrates it as
 *
 *     INT x^p E^(a x^m) dx  =  -(1/m) (-a)^(-s) Gamma[s, -a x^m],   s = (p+1)/m,
 *
 * which differentiates back exactly:  d/dx Gamma[s, u] = -u^(s-1) E^(-u) u', so
 * with u = -a x^m the prefactors cancel to x^(ms-1) E^(a x^m) = x^p E^(a x^m)
 * under the usual (-a x^m)^(s-1) = (-a)^(s-1) x^(m(s-1)) branch convention (the
 * same one Mathematica's answer for this family carries).
 *
 * Why this exists as its own recognizer, and why it is gated to a SYMBOLIC
 * exponent.  With p symbolic the integrand is outside every elementary stage of
 * the cascade -- there is no antiderivative in the integrand's own field -- and
 * the general stages do not merely decline, they SEARCH: `Integrate[x^n E^(-x), x]`
 * cost a measured 12.9 s to come back unevaluated, and 12.4 s for the Gaussian
 * sibling `Integrate[x^n E^(-x^2), x]`.  Two of those are the whole cost of
 * `DSolve[y''-y == x^n]`, whose variation-of-parameters particular is exactly this
 * pair of integrals.  A NUMERIC exponent is deliberately left alone: for a
 * non-negative integer p the answer is elementary and the existing stages give it,
 * and for a numeric non-integer they give the cleaner Erf/Gamma-free form, so
 * firing here would be a quality regression rather than a gain.
 *
 * The candidate is accepted only after an exact Simplify diff-back, so a
 * mis-recognition (a wrong branch, a mis-read exponent) can never emit a wrong
 * closed form -- it declines and the cascade continues.
 */

#include "integrate_gammapower.h"

#include "expr.h"
#include "eval.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>

static Expr* gp_sym(const char* s) { return expr_new_symbol(s); }
static Expr* gp_int(long n)        { return expr_new_integer(n); }
static Expr* gp_fn1(const char* h, Expr* a) {
    return expr_new_function(gp_sym(h), (Expr*[]){ a }, 1);
}
static Expr* gp_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(gp_sym(h), (Expr*[]){ a, b }, 2);
}
/* evaluate() BORROWS its argument, so free the input node. */
static Expr* gp_ev(Expr* node) { Expr* r = evaluate(node); expr_free(node); return r; }
static Expr* gp_ev2(const char* h, Expr* a, Expr* b) { return gp_ev(gp_fn2(h, a, b)); }

static bool gp_free_of(const Expr* e, const Expr* x) {
    if (!e) return true;
    if (expr_eq((Expr*)e, (Expr*)x)) return false;
    if (e->type != EXPR_FUNCTION) return true;
    if (!gp_free_of(e->data.function.head, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!gp_free_of(e->data.function.args[i], x)) return false;
    return true;
}

/* True when e carries a symbol other than x.  Numeric WRAPPER heads do not count:
 * Rational[1,2] and Complex[0,1] are EXPR_FUNCTIONs whose head is a symbol, so a
 * bare structural scan reads `Sqrt[x]`'s exponent 1/2 as "symbolic" and diverts a
 * numeric-exponent integrand that the later stages answer better.  The NumberQ
 * test below is the actual gate; this is the complementary "there is a parameter
 * in there" requirement. */
static bool gp_has_foreign_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return !expr_eq((Expr*)e, (Expr*)x);
    if (e->type != EXPR_FUNCTION) return false;
    if (gp_has_foreign_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (gp_has_foreign_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* NumberQ[e]: a literal Integer / Rational / Real / Complex.  Pi and Sqrt[2] are
 * NOT numbers by this test, which is what we want — a symbolic-constant exponent
 * has no elementary answer either, so it belongs to this recognizer. */
static bool gp_is_number(const Expr* e) {
    Expr* q = gp_ev(gp_fn1("NumberQ", expr_copy((Expr*)e)));
    bool r = (q && q->type == EXPR_SYMBOL && q->data.symbol.name == intern_symbol("True"));
    expr_free(q);
    return r;
}

/* Split a product into the factors matching x^p (symbolic p), E^(a x^m), and the
 * rest.  Returns false unless exactly one of each was found.  `*pexp`, `*acoef`
 * are fresh; `*m` is the integer power of x in the exponential's argument. */
static bool gp_split(Expr* f, Expr* x, Expr** pexp, Expr** acoef, long* m,
                     Expr** rest) {
    const char* s_times = intern_symbol("Times");
    const char* s_power = intern_symbol("Power");
    size_t nf = 1;
    Expr* one[1] = { f };
    Expr** fac = one;
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == s_times) {
        fac = f->data.function.args;
        nf  = f->data.function.arg_count;
    }

    Expr* pe = NULL; Expr* ac = NULL; long mm = 0;
    Expr** keep = malloc(nf * sizeof(Expr*));
    size_t nkeep = 0;
    for (size_t i = 0; i < nf; i++) {
        Expr* t = fac[i];
        if (t->type == EXPR_FUNCTION && t->data.function.head->type == EXPR_SYMBOL
            && t->data.function.head->data.symbol.name == s_power
            && t->data.function.arg_count == 2) {
            Expr* base = t->data.function.args[0];
            Expr* ex   = t->data.function.args[1];
            /* x^p with p symbolic and free of x */
            if (!pe && expr_eq(base, x) && gp_free_of(ex, x)
                && gp_has_foreign_symbol(ex, x) && !gp_is_number(ex)) {
                pe = expr_copy(ex);
                continue;
            }
            /* E^(a x^m): base is the symbol E, exponent is a monomial in x */
            if (!ac && base->type == EXPR_SYMBOL
                && base->data.symbol.name == intern_symbol("E")) {
                /* exponent = coefficient * x^m, m a positive integer */
                Expr* arg = ex;
                Expr* co = NULL; long mv = 0;
                if (expr_eq(arg, x)) { co = gp_int(1); mv = 1; }
                else if (arg->type == EXPR_FUNCTION
                         && arg->data.function.head->type == EXPR_SYMBOL
                         && arg->data.function.head->data.symbol.name == s_power
                         && arg->data.function.arg_count == 2
                         && expr_eq(arg->data.function.args[0], x)
                         && arg->data.function.args[1]->type == EXPR_INTEGER
                         && arg->data.function.args[1]->data.integer > 0) {
                    co = gp_int(1);
                    mv = (long)arg->data.function.args[1]->data.integer;
                } else if (arg->type == EXPR_FUNCTION
                           && arg->data.function.head->type == EXPR_SYMBOL
                           && arg->data.function.head->data.symbol.name == s_times) {
                    /* c * x^m (c free of x), in any factor order */
                    size_t na = arg->data.function.arg_count;
                    Expr** cf = malloc(na * sizeof(Expr*));
                    size_t ncf = 0;
                    for (size_t j = 0; j < na; j++) {
                        Expr* u = arg->data.function.args[j];
                        if (!mv && expr_eq(u, x)) { mv = 1; continue; }
                        if (!mv && u->type == EXPR_FUNCTION
                            && u->data.function.head->type == EXPR_SYMBOL
                            && u->data.function.head->data.symbol.name == s_power
                            && u->data.function.arg_count == 2
                            && expr_eq(u->data.function.args[0], x)
                            && u->data.function.args[1]->type == EXPR_INTEGER
                            && u->data.function.args[1]->data.integer > 0) {
                            mv = (long)u->data.function.args[1]->data.integer;
                            continue;
                        }
                        if (!gp_free_of(u, x)) { mv = 0; break; }
                        cf[ncf++] = expr_copy(u);
                    }
                    if (mv > 0) {
                        co = (ncf == 0) ? gp_int(1)
                           : (ncf == 1) ? cf[0]
                           : gp_ev(expr_new_function(gp_sym("Times"), cf, ncf));
                    } else {
                        for (size_t j = 0; j < ncf; j++) expr_free(cf[j]);
                    }
                    free(cf);
                }
                if (co && mv > 0) { ac = co; mm = mv; continue; }
                if (co) expr_free(co);
            }
        }
        if (!gp_free_of(t, x)) {       /* an x-dependent factor we cannot absorb */
            if (pe) expr_free(pe);
            if (ac) expr_free(ac);
            for (size_t j = 0; j < nkeep; j++) expr_free(keep[j]);
            free(keep);
            return false;
        }
        keep[nkeep++] = expr_copy(t);
    }

    if (!pe || !ac) {
        if (pe) expr_free(pe);
        if (ac) expr_free(ac);
        for (size_t j = 0; j < nkeep; j++) expr_free(keep[j]);
        free(keep);
        return false;
    }
    *pexp = pe; *acoef = ac; *m = mm;
    *rest = (nkeep == 0) ? gp_int(1)
          : (nkeep == 1) ? keep[0]
          : gp_ev(expr_new_function(gp_sym("Times"), keep, nkeep));
    free(keep);
    return true;
}

Expr* integrate_gammapower_try(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;

    Expr* p = NULL; Expr* a = NULL; Expr* K = NULL; long m = 0;
    if (!gp_split(f, x, &p, &a, &m, &K)) return NULL;

    /* s = (p+1)/m,  z = -a x^m,  result = -K/m (-a)^(-s) Gamma[s, z]. */
    Expr* s = gp_ev2("Divide", gp_ev2("Plus", p, gp_int(1)), gp_int(m));
    Expr* na = gp_ev(gp_fn1("Minus", expr_copy(a)));
    Expr* z  = gp_ev2("Times", expr_copy(na),
                      (m == 1) ? expr_copy(x) : gp_ev2("Power", expr_copy(x), gp_int(m)));
    Expr* pref = gp_ev2("Power", na, gp_ev(gp_fn1("Minus", expr_copy(s))));
    Expr* gam  = gp_ev2("Gamma", expr_copy(s), z);
    Expr* body = gp_ev2("Times",
                     gp_ev2("Times", gp_ev2("Divide", gp_ev(gp_fn1("Minus", K)), gp_int(m)),
                            pref),
                     gam);
    expr_free(a); expr_free(s);

    /* Diff-back is the only acceptance criterion.  Two passes, because the closed
     * form is stated in the principal-branch convention
     * (-a x^m)^k == (-a)^k x^(m k), which Simplify does not apply on its own:
     *   - a < 0, m == 1 verifies exactly (nothing to re-branch), the common case;
     *   - a > 0 (so -a < 0) or m > 1 leaves (-x)^k / (x^2)^k in the derivative, and
     *     only PowerExpand — i.e. adopting the very convention the formula is
     *     written in — collapses it.
     * The second pass therefore quotients out exactly that convention and nothing
     * else: it is still a differentiation certificate, not a weaker test, and a
     * genuinely mis-read exponent or coefficient fails both passes and declines. */
    Expr* d = gp_ev2("D", expr_copy(body), expr_copy(x));
    Expr* raw = gp_ev2("Subtract", d, expr_copy(f));
    Expr* diff = gp_ev(gp_fn1("Simplify", expr_copy(raw)));
    bool ok = (diff && diff->type == EXPR_INTEGER && diff->data.integer == 0);
    expr_free(diff);
    if (!ok) {
        Expr* pe2 = gp_ev(gp_fn1("Simplify", gp_ev(gp_fn1("PowerExpand", expr_copy(raw)))));
        ok = (pe2 && pe2->type == EXPR_INTEGER && pe2->data.integer == 0);
        expr_free(pe2);
    }
    expr_free(raw);
    if (!ok) { expr_free(body); return NULL; }
    return body;
}
