/* integrate_quadrad.c
 *
 * Integration of rational functions of x and square roots of a common quadratic
 * argument -- the classical Euler substitution.  See integrate_quadrad.h for the
 * high-level description.
 *
 * ---------------------------------------------------------------------
 * Algorithm
 * ---------------------------------------------------------------------
 * Identify the integrand as
 *
 *     R[x, (a x^2 + b x + c)^(n1/2), (a x^2 + b x + c)^(n2/2), ...]
 *
 * with R rational, a, b, c free of x, every radical a square root (denominator
 * 2) and a single shared degree-2 radicand rad = a x^2 + b x + c:
 *
 *   1. Scan f for every Power[base, p/q] whose x-dependent base is a degree-2
 *      polynomial.  Each such radical must have q == 2 and all bases must be
 *      structurally equal.  Reject otherwise.
 *   2. Substitute the radicals to a fresh symbol y via poly_subst_radical_to_gen
 *      (rad^(p/2) -> y^p, bare rad -> y^2).  exy = R(x, y) must be rational in
 *      {x, y}.
 *   3. Choose exactly ONE Euler substitution so the auxiliary constants stay
 *      real (see qr_select):
 *        a > 0                       Euler I    : y = Sqrt[a] x + u
 *        a < 0 and b^2-4ac > 0       Euler III  : y = (x - alpha) u
 *        a symbolic                  Euler I    (best-effort, a > 0 branch)
 *      Each builds x = X(u), the radical image y = ysub(u), and dx = X'(u) du.
 *   4. Form the rationalised integrand exy * dx /. {y -> ysub, x -> X}, Cancel/
 *      Together it, and recurse: Integrate[., u].
 *   5. Back-substitute u -> U(x).  The Euler substitution is an exact bijection
 *      on the relevant domain, so the antiderivative is correct by construction
 *      once the rational sub-integral closes -- no differentiate-back check is
 *      performed (it is unnecessary and, for symbolic parameters, prohibitively
 *      expensive).  Mirrors integrate_linrad.c.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take` consumes its
 * argument.  We never expr_free(res).  The small helper set mirrors
 * integrate_linrad.c (each integration-method file keeps a private copy).
 */

#include "integrate_quadrad.h"

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

/* ReplaceAll[expr, {rule1, rule2}] -- simultaneous list of rules.  Consumes
 * `expr`, `rule1`, `rule2`; returns owned (evaluated). */
static Expr* replace_rules2(Expr* expr, Expr* rule1, Expr* rule2) {
    Expr* rules = expr_new_function(expr_new_symbol("List"),
                      (Expr*[]){ rule1, rule2 }, 2);
    Expr* call = internal_replace_all((Expr*[]){ expr, rules }, 2);
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

/* True if `e` has at least one free variable (Variables[e] =!= {}); used for the
 * symbolic-leading-coefficient branch.  Borrows e. */
static bool has_variables(Expr* e) {
    Expr* v = eval_take(internal_variables((Expr*[]){ expr_copy(e) }, 1));
    bool nonempty = v && v->type == EXPR_FUNCTION && head_is(v, SYM_List)
                    && v->data.function.arg_count > 0;
    if (v) expr_free(v);
    return nonempty;
}

static bool is_int_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define QR_MAX_DEPTH 8
static int qr_depth = 0;
static unsigned long qr_sym_counter = 0;

/* ---------------------------------------------------------------------- */
/* Identification: the shared degree-2 radicand                           */
/* ---------------------------------------------------------------------- */

/* Walk `e` collecting every Power[base, p/q] whose x-dependent base is a
 * degree-2 polynomial in x.  Each such radical must be a square root (q == 2)
 * and all bases must be structurally equal.
 *  *first_rad : owned copy of the shared radicand (NULL until found).
 *  *found     : count flag (0 = none yet).
 *  *ok        : cleared on any rejecting condition. */
static void qr_scan(Expr* e, Expr* x, Expr** first_rad, int* found, bool* ok) {
    if (!*ok || !e) return;
    if (e->type != EXPR_FUNCTION) return;

    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q)) {
            int64_t qa = q < 0 ? -q : q;
            if (qa > 1 && !expr_free_of(base, x)) {
                /* A genuine radical whose argument depends on x. */
                Expr* vars[1] = { x };
                if (is_polynomial(base, vars, 1) && get_degree_poly(base, x) == 2) {
                    if (qa != 2) {
                        *ok = false;     /* non-square root of the quadratic */
                        return;
                    }
                    if (*found == 0) {
                        *first_rad = expr_copy(base);
                        *found = 1;
                    } else if (!expr_eq(base, *first_rad)) {
                        *ok = false;     /* distinct quadratic radicands */
                        return;
                    }
                }
                /* Radical of a non-degree-2 base (linear, cubic, ...): left in
                 * place -- the downstream rationality check rejects it. */
            }
        }
    }

    qr_scan(e->data.function.head, x, first_rad, found, ok);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        qr_scan(e->data.function.args[i], x, first_rad, found, ok);
    }
}

/* Sign of the discriminant b^2 - 4 a c: +1, -1, or 0 (unknown). */
static int qr_disc_sign(Expr* a, Expr* b, Expr* c) {
    Expr* disc = eval_take(mk_fn2("Plus",
        mk_fn2("Power", expr_copy(b), mk_int(2)),
        mk_fn3("Times", mk_int(-4), expr_copy(a), expr_copy(c))));
    int s = disc ? expr_numeric_sign(disc) : 0;
    if (disc) expr_free(disc);
    return s;
}

/* ---------------------------------------------------------------------- */
/* Euler substitutions: build X(u), ysub(u), dx(u), U(x)                  */
/* ---------------------------------------------------------------------- */

typedef struct {
    Expr* X;     /* x as a function of u                       */
    Expr* ysub;  /* the radical Sqrt[rad] as a function of u   */
    Expr* dx;    /* dx/du                                      */
    Expr* U;     /* u as a function of x (back-substitution)   */
} QrSub;

/* Euler's first substitution (a > 0, or symbolic a): Sqrt[rad] = Sqrt[a] x + u.
 *   X  = (u^2 - c)/(b - 2 Sqrt[a] u)
 *   dx = 2(-Sqrt[a] u^2 + b u - Sqrt[a] c)/(b - 2 Sqrt[a] u)^2
 *   U  = Sqrt[rad] - Sqrt[a] x
 * All coefficient arguments are borrowed. */
static QrSub qr_euler1(Expr* a, Expr* b, Expr* c, Expr* rad, Expr* x, Expr* u) {
    QrSub s;
    Expr* sa  = mk_fn1("Sqrt", expr_copy(a));
    Expr* den = mk_fn2("Plus", expr_copy(b),
                    mk_fn3("Times", mk_int(-2), expr_copy(sa), expr_copy(u)));

    s.X = mk_fn2("Times",
            mk_fn2("Plus", mk_fn2("Power", expr_copy(u), mk_int(2)),
                           mk_fn2("Times", mk_int(-1), expr_copy(c))),
            mk_fn2("Power", expr_copy(den), mk_int(-1)));

    s.ysub = mk_fn2("Plus", mk_fn2("Times", expr_copy(sa), expr_copy(s.X)),
                            expr_copy(u));

    {
        Expr* dxnum = mk_fn3("Plus",
            mk_fn3("Times", mk_int(-1), expr_copy(sa),
                   mk_fn2("Power", expr_copy(u), mk_int(2))),
            mk_fn2("Times", expr_copy(b), expr_copy(u)),
            mk_fn3("Times", mk_int(-1), expr_copy(sa), expr_copy(c)));
        s.dx = mk_fn3("Times", mk_int(2), dxnum,
            mk_fn2("Power", mk_fn2("Power", expr_copy(den), mk_int(2)), mk_int(-1)));
    }

    s.U = mk_fn2("Plus", mk_fn1("Sqrt", expr_copy(rad)),
            mk_fn3("Times", mk_int(-1), expr_copy(sa), expr_copy(x)));

    expr_free(sa);
    expr_free(den);
    return s;
}

/* Euler's third substitution (a < 0, b^2-4ac > 0): Sqrt[rad] = (x - alpha) u,
 * with alpha, beta the real roots of rad.
 *   X  = (a beta - alpha u^2)/(a - u^2)
 *   dx = 2 u a (beta - alpha)/(u^2 - a)^2
 *   U  = Sqrt[rad]/(x - alpha)
 * All coefficient arguments are borrowed. */
static QrSub qr_euler3(Expr* a, Expr* b, Expr* c, Expr* rad, Expr* x, Expr* u) {
    QrSub s;
    Expr* disc  = mk_fn2("Plus", mk_fn2("Power", expr_copy(b), mk_int(2)),
                      mk_fn3("Times", mk_int(-4), expr_copy(a), expr_copy(c)));
    Expr* sdisc = mk_fn1("Sqrt", disc);                 /* consumes disc */
    Expr* two_a = mk_fn2("Times", mk_int(2), expr_copy(a));

    Expr* alpha = mk_fn2("Times",
        mk_fn2("Plus", mk_fn2("Times", mk_int(-1), expr_copy(b)),
                       mk_fn2("Times", mk_int(-1), expr_copy(sdisc))),
        mk_fn2("Power", expr_copy(two_a), mk_int(-1)));
    Expr* beta = mk_fn2("Times",
        mk_fn2("Plus", mk_fn2("Times", mk_int(-1), expr_copy(b)), expr_copy(sdisc)),
        mk_fn2("Power", expr_copy(two_a), mk_int(-1)));

    {
        Expr* xnum = mk_fn2("Plus", mk_fn2("Times", expr_copy(a), expr_copy(beta)),
            mk_fn3("Times", mk_int(-1), expr_copy(alpha),
                   mk_fn2("Power", expr_copy(u), mk_int(2))));
        Expr* xden = mk_fn2("Plus", expr_copy(a),
            mk_fn2("Times", mk_int(-1), mk_fn2("Power", expr_copy(u), mk_int(2))));
        s.X = mk_fn2("Times", xnum, mk_fn2("Power", expr_copy(xden), mk_int(-1)));
        expr_free(xden);
    }

    s.ysub = mk_fn2("Times",
        mk_fn2("Plus", expr_copy(s.X), mk_fn2("Times", mk_int(-1), expr_copy(alpha))),
        expr_copy(u));

    {
        Expr* dxnum = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_int(2), expr_copy(u), expr_copy(a),
            mk_fn2("Plus", expr_copy(beta), mk_fn2("Times", mk_int(-1), expr_copy(alpha)))
        }, 4);
        Expr* dxden = mk_fn2("Power",
            mk_fn2("Plus", mk_fn2("Power", expr_copy(u), mk_int(2)),
                           mk_fn2("Times", mk_int(-1), expr_copy(a))), mk_int(2));
        s.dx = mk_fn2("Times", dxnum, mk_fn2("Power", expr_copy(dxden), mk_int(-1)));
        expr_free(dxden);
    }

    s.U = mk_fn2("Times", mk_fn1("Sqrt", expr_copy(rad)),
        mk_fn2("Power", mk_fn2("Plus", expr_copy(x),
            mk_fn2("Times", mk_int(-1), expr_copy(alpha))), mk_int(-1)));

    expr_free(sdisc);
    expr_free(two_a);
    expr_free(alpha);
    expr_free(beta);
    return s;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

static Expr* qr_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (expr_free_of(f, x))      return NULL;   /* nothing to integrate in x */
    if (qr_depth >= QR_MAX_DEPTH) return NULL;

    Expr* rad = NULL, *rade = NULL;
    Expr* a = NULL, *b = NULL, *c = NULL;
    Expr* exy = NULL, *y = NULL, *u = NULL, *su = NULL;
    Expr* result = NULL;

    /* --- 1. Identify the shared degree-2 radicand. --- */
    int found = 0;
    bool ok = true;
    qr_scan(f, x, &rad, &found, &ok);
    if (!ok || !found || !rad) goto done;

    /* a x^2 + b x + c coefficients (a != 0 since deg == 2). */
    rade = eval_take(internal_expand((Expr*[]){ expr_copy(rad) }, 1));
    if (!rade) goto done;
    a = get_coeff(rade, x, 2);
    b = get_coeff(rade, x, 1);
    c = get_coeff(rade, x, 0);
    if (!a || !b || !c) goto done;
    if (is_int_zero(b) && is_int_zero(c)) goto done;   /* degenerate Sqrt[a x^2] */

    /* --- 3. Choose exactly one real-valued Euler substitution. --- */
    int branch;
    {
        int sa = expr_numeric_sign(a);
        if (sa > 0) {
            branch = 1;
        } else if (sa < 0 && qr_disc_sign(a, b, c) > 0) {
            branch = 3;
        } else if (has_variables(a)) {
            branch = 1;
        } else {
            goto done;   /* no real branch (e.g. a < 0, disc <= 0) */
        }
    }

    /* --- 2. Substitute radicals -> y and require a rational function. --- */
    char yname[96];
    snprintf(yname, sizeof yname, "Integrate`QuadraticRadicals`y$%lu", qr_sym_counter++);
    y = expr_new_symbol(yname);
    {
        Expr* one = mk_int(1);
        exy = poly_subst_radical_to_gen(f, rad, one, 2, yname);
        expr_free(one);
    }
    if (exy) exy = cancel_together(exy);
    if (!exy) goto done;
    {
        Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(exy)));
        Expr* den = eval_take(mk_fn1("Denominator", expr_copy(exy)));
        Expr* vars2[2] = { x, y };
        bool rat_ok = num && den &&
                      is_polynomial(num, vars2, 2) &&
                      is_polynomial(den, vars2, 2);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (!rat_ok) goto done;
    }

    /* Fresh substitution variable u. */
    char uname[96];
    snprintf(uname, sizeof uname, "Integrate`QuadraticRadicals`u$%lu", qr_sym_counter++);
    u = expr_new_symbol(uname);

    /* --- 4. Rationalise: exy * dx /. {y -> ysub, x -> X}, then Cancel. --- */
    {
        QrSub s = (branch == 1) ? qr_euler1(a, b, c, rad, x, u)
                                : qr_euler3(a, b, c, rad, x, u);
        su = s.U;                                       /* freed at done */
        Expr* prod  = mk_fn2("Times", expr_copy(exy), s.dx);   /* consumes s.dx */
        Expr* ruleY = mk_fn2("Rule", expr_copy(y), s.ysub);    /* consumes s.ysub */
        Expr* ruleX = mk_fn2("Rule", expr_copy(x), s.X);       /* consumes s.X */
        Expr* integrand = replace_rules2(prod, ruleY, ruleX);  /* consumes all */
        if (integrand) integrand = cancel_together(integrand);

        /* --- 5. Recurse and back-substitute. --- */
        if (integrand) {
            qr_depth++;
            Expr* G = integrate_in(integrand, u);
            qr_depth--;
            expr_free(integrand);
            if (G) result = replace_one(G, u, su);     /* consumes G; borrows su */
        }
    }

done:
    if (rad)  expr_free(rad);
    if (rade) expr_free(rade);
    if (a)    expr_free(a);
    if (b)    expr_free(b);
    if (c)    expr_free(c);
    if (exy)  expr_free(exy);
    if (y)    expr_free(y);
    if (u)    expr_free(u);
    if (su)   expr_free(su);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_quadrad_try(Expr* f, Expr* x) {
    return qr_core(f, x);
}

Expr* builtin_integrate_quadrad(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return qr_core(f, x);
}

void integrate_quadrad_init(void) {
    symtab_add_builtin("Integrate`QuadraticRadicals", builtin_integrate_quadrad);
    symtab_get_def("Integrate`QuadraticRadicals")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`QuadraticRadicals",
        "Integrate`QuadraticRadicals[f, x] integrates a rational function of x and\n"
        "square roots (a x^2 + b x + c)^(m/2) sharing one quadratic argument. It\n"
        "applies a single real-valued Euler substitution chosen by the sign of a\n"
        "(and, when a < 0, of the discriminant), reduces f to a rational function\n"
        "of u, integrates that, and back-substitutes. The substitution is exact,\n"
        "so the result is not re-verified. Strict: returns unevaluated when f is\n"
        "not of this form or the reduced integral does not close.");
}
