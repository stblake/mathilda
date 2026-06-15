/* integrate_linratiorad.c
 *
 * Integration of rational functions of x and radicals of a common *linear
 * fractional* (Möbius) argument -- the "linear ratio radicals" substitution.
 * See integrate_linratiorad.h for the high-level description.
 *
 * ---------------------------------------------------------------------
 * Algorithm
 * ---------------------------------------------------------------------
 * Identify the integrand as
 *
 *     R[x, ((a x + b)/(c x + d))^(n1/m1), ((a x + b)/(c x + d))^(n2/m2), ...]
 *
 * with R rational, a, b, c, d free of x (a d - b c != 0), and a single shared
 * linear-fractional base r = (a x + b)/(c x + d):
 *
 *   1. Scan f for every Power[base, p/q] with |q| > 1 whose base depends on x.
 *      All such bases must be structurally equal; m = LCM of the |q|.  Reject
 *      on distinct bases.  (No polynomiality is required here -- the ratio base
 *      `Times[a x + b, Power[c x + d, -1]]` is not a polynomial.  That is also
 *      what cleanly partitions this method from LinearRadicals / Quadratic-
 *      Radicals, whose scans demand a polynomial base.)
 *   2. Canonicalise the shared base with Cancel[Together[.]] and read off the
 *      Möbius coefficients: a, b from the numerator, c, d from the denominator.
 *      Require numerator degree <= 1, denominator degree == 1 (a genuine
 *      ratio -- a constant denominator is LinearRadicals' job) and
 *      a d - b c != 0.
 *   3. Substitute the radicals to a fresh symbol u via poly_subst_radical_to_gen
 *      (shared with LinearRadicals / the algebraic-factoring path):
 *      r^(p/q) -> u^(p m / q), bare r -> u^m.  The result exy = R(x, u) must be
 *      a rational function of {x, u}.
 *   4. Form the rationalised integrand
 *          exy · m (a d - b c) u^(m-1)/(a - c u^m)^2  /.  x -> (d u^m - b)/(a - c u^m),
 *      Cancel[Together[.]] it, and recurse: Integrate[., u].
 *   5. Back-substitute u -> ((a x + b)/(c x + d))^(1/m).  The Möbius
 *      substitution is an exact bijection introducing no branch issues, so the
 *      antiderivative is correct by construction once the rational sub-integral
 *      closes -- no differentiate-back verification is performed (it is
 *      unnecessary and, for symbolic parameters, prohibitively expensive).
 *      Mirrors integrate_linrad.c.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take` consumes its
 * argument.  We never expr_free(res).  The small helper set mirrors
 * integrate_linrad.c (each integration-method file keeps a private copy).
 */

#include "integrate_linratiorad.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "poly.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirrors integrate_linrad.c)       */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(name), args, 2);
}

static Expr* mk_fn3(const char* name, Expr* a, Expr* b, Expr* c) {
    Expr* args[3] = { a, b, c };
    return expr_new_function(expr_new_symbol(name), args, 3);
}

/* -b c, built as Times[-1, b, c]; consumes b and c. */
static Expr* mk_neg_prod(Expr* b, Expr* c) {
    return mk_fn3("Times", mk_int(-1), b, c);
}

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True if `f` contains no subexpression structurally equal to `x`. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_unintegrated(e->data.function.args[i])) return true;
    }
    return false;
}

/* Cancel[Together[e]], consuming `e`. */
static Expr* cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* ReplaceAll[expr, from -> to]; consumes `expr`; returns owned (evaluated). */
static Expr* replace_one(Expr* expr, const Expr* from, const Expr* to) {
    Expr* rule = mk_fn2("Rule", expr_copy((Expr*)from), expr_copy((Expr*)to));
    Expr* call = internal_replace_all((Expr*[]){ expr, rule }, 2);
    return eval_take(call);
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2("Integrate", expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* True if `e` is the integer/real literal 0. */
static bool is_numeric_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    return false;
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define LRR_MAX_DEPTH 8
static int lrr_depth = 0;
static unsigned long lrr_sym_counter = 0;

/* ---------------------------------------------------------------------- */
/* Identification: the shared linear-fractional base and its radical LCM  */
/* ---------------------------------------------------------------------- */

/* Walk `e` collecting every Power[base, p/q] (|q| > 1) whose base depends on x.
 * All such bases must be structurally equal; m = running LCM of the radical
 * denominators.  Polynomiality is *not* required -- the base is a ratio.
 *  *first_base : owned copy of the shared base (NULL until found).
 *  *m          : running LCM of the radical denominators.
 *  *found      : count flag (0 = none yet).
 *  *ok         : cleared on any rejecting condition. */
static void lrr_scan(Expr* e, Expr* x, Expr** first_base,
                     int64_t* m, int* found, bool* ok) {
    if (!*ok || !e) return;
    if (e->type != EXPR_FUNCTION) return;

    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && (q > 1 || q < -1)) {
            if (!expr_free_of(base, x)) {
                int64_t qa = q < 0 ? -q : q;
                if (*found == 0) {
                    *first_base = expr_copy(base);
                    *m = qa;
                    *found = 1;
                } else if (!expr_eq(base, *first_base)) {
                    *ok = false;   /* distinct radical bases */
                    return;
                } else {
                    *m = lcm(*m, qa);
                }
            }
            /* base free of x: a constant radical -- ignored. */
        }
    }

    lrr_scan(e->data.function.head, x, first_base, m, found, ok);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        lrr_scan(e->data.function.args[i], x, first_base, m, found, ok);
    }
}

/* Read the Möbius coefficients a, b, c, d of base = (a x + b)/(c x + d).
 * Returns true (filling owned a,b,c,d) when the canonicalised base is a genuine
 * linear-fractional function of x: numerator degree <= 1, denominator degree
 * == 1, and a d - b c != 0.  On false nothing is allocated.  Borrows base/x. */
static bool lrr_mobius_coeffs(Expr* base, Expr* x,
                              Expr** a, Expr** b, Expr** c, Expr** d) {
    bool ok = false;
    Expr* tb   = cancel_together(expr_copy(base));     /* canonical ratio */
    Expr* num  = NULL, *den = NULL, *nume = NULL, *dene = NULL;
    Expr* aa = NULL, *bb = NULL, *cc = NULL, *dd = NULL, *det = NULL;
    if (!tb) goto done;

    num = eval_take(mk_fn1("Numerator",   expr_copy(tb)));
    den = eval_take(mk_fn1("Denominator", expr_copy(tb)));
    if (!num || !den) goto done;
    nume = eval_take(internal_expand((Expr*[]){ expr_copy(num) }, 1));
    dene = eval_take(internal_expand((Expr*[]){ expr_copy(den) }, 1));
    if (!nume || !dene) goto done;

    /* Both numerator and denominator must be polynomials in x, with the
     * denominator genuinely linear (degree 1) -- a constant denominator means
     * a bare linear base, which belongs to Integrate`LinearRadicals. */
    {
        Expr* vars[1] = { x };
        if (!is_polynomial(nume, vars, 1) || !is_polynomial(dene, vars, 1))
            goto done;
        if (get_degree_poly(nume, x) > 1) goto done;
        if (get_degree_poly(dene, x) != 1) goto done;
    }

    aa = get_coeff(nume, x, 1);
    bb = get_coeff(nume, x, 0);
    cc = get_coeff(dene, x, 1);
    dd = get_coeff(dene, x, 0);
    if (!aa || !bb || !cc || !dd) goto done;

    /* a d - b c != 0 (degenerate Möbius maps to a constant). */
    det = cancel_together(mk_fn2("Plus",
        mk_fn2("Times", expr_copy(aa), expr_copy(dd)),
        mk_neg_prod(expr_copy(bb), expr_copy(cc))));
    if (det && is_numeric_zero(det)) goto done;

    *a = aa; *b = bb; *c = cc; *d = dd;
    aa = bb = cc = dd = NULL;     /* ownership transferred */
    ok = true;

done:
    if (tb)   expr_free(tb);
    if (num)  expr_free(num);
    if (den)  expr_free(den);
    if (nume) expr_free(nume);
    if (dene) expr_free(dene);
    if (aa)   expr_free(aa);
    if (bb)   expr_free(bb);
    if (cc)   expr_free(cc);
    if (dd)   expr_free(dd);
    if (det)  expr_free(det);
    return ok;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

static Expr* lrr_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (expr_free_of(f, x))       return NULL;   /* nothing to integrate in x */
    if (lrr_depth >= LRR_MAX_DEPTH) return NULL;

    /* --- 1. Identify the shared base and the radical LCM m. --- */
    Expr* base = NULL;
    int64_t m = 1;
    int found = 0;
    bool ok = true;
    lrr_scan(f, x, &base, &m, &found, &ok);
    if (!ok || !found || m < 2) { if (base) expr_free(base); return NULL; }

    /* --- 2. Möbius coefficients of base = (a x + b)/(c x + d). --- */
    Expr* a = NULL, *b = NULL, *c = NULL, *d = NULL;
    if (!lrr_mobius_coeffs(base, x, &a, &b, &c, &d)) { expr_free(base); return NULL; }

    /* Fresh substitution variable u. */
    char uname[96];
    snprintf(uname, sizeof(uname), "Integrate`LinearRatioRadicals`u$%lu",
             lrr_sym_counter++);
    Expr* u = expr_new_symbol(uname);

    Expr* result = NULL;
    Expr* exy = NULL, *integrand = NULL, *G = NULL;

    /* --- 3. Substitute radicals -> u and require a rational function. --- */
    {
        Expr* one = mk_int(1);
        exy = poly_subst_radical_to_gen(f, base, one, m, uname);
        expr_free(one);
    }
    if (exy) exy = cancel_together(exy);
    if (!exy) goto done;
    {
        Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(exy)));
        Expr* den = eval_take(mk_fn1("Denominator", expr_copy(exy)));
        Expr* vars2[2] = { x, u };
        bool rat_ok = num && den &&
                      is_polynomial(num, vars2, 2) &&
                      is_polynomial(den, vars2, 2);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (!rat_ok) goto done;
    }

    /* --- 4. exy · m (a d - b c) u^(m-1)/(a - c u^m)^2  /.  x -> X(u). --- */
    {
        /* aceum = a - c u^m ; jacobian denominator (a - c u^m)^2. */
        Expr* det = mk_fn2("Plus",
            mk_fn2("Times", expr_copy(a), expr_copy(d)),
            mk_neg_prod(expr_copy(b), expr_copy(c)));        /* a d - b c */
        Expr* aceum = mk_fn2("Plus", expr_copy(a),
            mk_fn3("Times", mk_int(-1), expr_copy(c),
                   mk_fn2("Power", expr_copy(u), mk_int(m))));
        Expr* jac = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
            mk_int(m),
            det,                                            /* consumes det */
            mk_fn2("Power", expr_copy(u), mk_int(m - 1)),
            mk_fn2("Power", mk_fn2("Power", expr_copy(aceum), mk_int(2)),
                   mk_int(-1)),                             /* (a - c u^m)^-2 */
            exy                                             /* consumes exy */
        }, 5);
        exy = NULL;
        /* x -> (d u^m - b)/(a - c u^m). */
        Expr* sub = mk_fn2("Times",
            mk_fn2("Plus",
                mk_fn2("Times", expr_copy(d), mk_fn2("Power", expr_copy(u), mk_int(m))),
                mk_fn2("Times", mk_int(-1), expr_copy(b))),
            mk_fn2("Power", expr_copy(aceum), mk_int(-1)));
        integrand = replace_one(jac, x, sub);               /* consumes jac */
        expr_free(sub);
        expr_free(aceum);
        if (integrand) integrand = cancel_together(integrand);
    }
    if (!integrand) goto done;

    /* --- 5. Recurse and back-substitute. --- */
    lrr_depth++;
    G = integrate_in(integrand, u);
    lrr_depth--;
    if (!G) goto done;

    {
        /* u -> ((a x + b)/(c x + d))^(1/m), i.e. base^(1/m). */
        Expr* back = mk_fn2("Power", expr_copy(base), make_rational(1, m));
        result = replace_one(G, u, back);                   /* consumes G */
        G = NULL;
        expr_free(back);
    }

done:
    expr_free(base);
    expr_free(a);
    expr_free(b);
    expr_free(c);
    expr_free(d);
    expr_free(u);
    if (exy)       expr_free(exy);
    if (integrand) expr_free(integrand);
    if (G)         expr_free(G);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_linratiorad_try(Expr* f, Expr* x) {
    return lrr_core(f, x);
}

Expr* builtin_integrate_linratiorad(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return lrr_core(f, x);
}

void integrate_linratiorad_init(void) {
    symtab_add_builtin("Integrate`LinearRatioRadicals", builtin_integrate_linratiorad);
    symtab_get_def("Integrate`LinearRatioRadicals")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`LinearRatioRadicals",
        "Integrate`LinearRatioRadicals[f, x] integrates a rational function of x\n"
        "and radicals ((a x + b)/(c x + d))^(m/n) sharing one linear-fractional\n"
        "argument. It substitutes u = ((a x + b)/(c x + d))^(1/n) (n = LCM of the\n"
        "radical denominators), inverts the Mobius map to reduce f to a rational\n"
        "function of u, integrates that, and back-substitutes. The substitution is\n"
        "exact, so the result is not re-verified. Strict: returns unevaluated when\n"
        "f is not of this form or the reduced integral does not close.");
}
