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

/* Small algebraic builders.  Every one ADOPTS its Expr* arguments (expr_new_function
 * takes ownership of the pointer array), matching the nested-mk_fn style below. */
static Expr* mk_neg(Expr* e)     { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_inv(Expr* e)     { return mk_fn2("Power", e, mk_int(-1)); }
static Expr* mk_pow_int(Expr* base, int64_t n) { return mk_fn2("Power", base, mk_int(n)); }
static Expr* mk_sqrt(Expr* e)    { return mk_fn1("Sqrt", e); }

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
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List),
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
        Expr* dxnum = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
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
/* Closed-form path (preferred): partial fractions + standard reductions. */
/*                                                                        */
/* The Euler substitution below always closes but emits Log-of-radical    */
/* antiderivatives.  For the common shape rational(x)/Sqrt[rad] this path  */
/* instead writes the answer with ArcTan / ArcTanh / ArcSinh / ArcSin /    */
/* Log directly (matching Mathematica): split the integrand as            */
/*     f = E(x) + H(x)/Sqrt[rad]        (E, H rational in x)               */
/* integrate E as a rational function, and reduce H by partial fractions   */
/* into the standard building blocks                                       */
/*     I_n  = Integrate[x^n / Sqrt[rad], x]      (reduction down to I0)    */
/*     J(p) = Integrate[1/((x-p) Sqrt[rad]), x]  (ArcTanh / ArcTan)        */
/* Each block is an exact identity (differentiates back to its integrand), */
/* so the assembled sum needs no verification.  Any shape it cannot fully  */
/* decompose (higher-order poles, irreducible-quadratic denominators, a    */
/* non-closing rational part) returns NULL and qr_core falls back to the   */
/* Euler substitution -- nothing that closes today stops closing.          */
/* ---------------------------------------------------------------------- */

/* Best-effort numeric sign of a parameter expression: exact test first,
 * then a probe substituting each free variable with a distinct small positive
 * rational and taking the sign of N[.].  The choice it feeds (ArcTanh vs ArcTan)
 * is cosmetic -- both differentiate back to the same integrand -- so a wrong or
 * unknown (0) guess is harmless. */
static int qr_sign_probe(Expr* e) {
    int s = expr_numeric_sign(e);
    if (s != 0) return s;
    Expr* vars = eval_take(mk_fn1("Variables", expr_copy(e)));
    if (!vars || !head_is(vars, SYM_List)) { if (vars) expr_free(vars); return 0; }
    static const int64_t pr[] = { 2, 3, 5, 7, 11, 13, 17, 19 };
    Expr* subd = expr_copy(e);
    for (size_t i = 0; i < vars->data.function.arg_count && subd; i++) {
        Expr* rule = mk_fn2("Rule", expr_copy(vars->data.function.args[i]),
                                    make_rational(pr[i % 8], pr[(i + 3) % 8]));
        subd = eval_take(internal_replace_all((Expr*[]){ subd, rule }, 2));
    }
    expr_free(vars);
    if (!subd) return 0;
    Expr* nn = eval_take(mk_fn2("N", subd, mk_int(30)));
    int r = nn ? expr_numeric_sign(nn) : 0;
    if (nn) expr_free(nn);
    return r;
}

/* Fold a polynomial in y (= Sqrt[rad], so y^2 -> rad) to its even (parity 0) or
 * odd (parity 1) part as a polynomial in x: sum over matching powers of
 * Coefficient[poly, y, pw] * rad^((pw-parity)/2).  Borrows poly,y,rad. */
static Expr* qr_fold_y(Expr* poly, Expr* y, Expr* rad, int parity) {
    int d = get_degree_poly(poly, y);
    Expr* acc = mk_int(0);
    for (int pw = parity; pw >= 0 && pw <= d; pw += 2) {
        Expr* ck = get_coeff(poly, y, pw);
        if (!ck) continue;
        acc = mk_fn2("Plus", acc,
            mk_fn2("Times", ck, mk_pow_int(expr_copy(rad), (pw - parity) / 2)));
    }
    return cancel_together(acc);
}

/* Split f = E(x) + H(x)/Sqrt[rad] with E,H rational in x, via y = Sqrt[rad] and
 * the conjugate over y^2 = rad.  Sets *E_out,*H_out (owned) and returns true;
 * false on failure.  Borrows f,rad,x. */
static bool qr_split(Expr* f, Expr* rad, Expr* x, Expr** E_out, Expr** H_out) {
    char yname[96];
    snprintf(yname, sizeof yname, "Integrate`QuadraticRadicals`g$%lu", qr_sym_counter++);
    Expr* y = expr_new_symbol(yname);
    Expr* one = mk_int(1);
    Expr* exy = poly_subst_radical_to_gen(f, rad, one, 2, yname);
    expr_free(one);
    if (exy) exy = cancel_together(exy);
    if (!exy) { expr_free(y); return false; }

    Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(exy)));
    Expr* den = eval_take(mk_fn1("Denominator", expr_copy(exy)));
    expr_free(exy);
    Expr* vars2[2] = { x, y };
    if (!num || !den || !is_polynomial(num, vars2, 2) || !is_polynomial(den, vars2, 2)) {
        if (num) expr_free(num); if (den) expr_free(den); expr_free(y); return false;
    }

    Expr* n0 = qr_fold_y(num, y, rad, 0), *n1 = qr_fold_y(num, y, rad, 1);
    Expr* d0 = qr_fold_y(den, y, rad, 0), *d1 = qr_fold_y(den, y, rad, 1);
    expr_free(num); expr_free(den); expr_free(y);
    if (!n0 || !n1 || !d0 || !d1) {
        if (n0) expr_free(n0); if (n1) expr_free(n1);
        if (d0) expr_free(d0); if (d1) expr_free(d1); return false;
    }

    /* denom = d0^2 - rad d1^2 (the norm of the conjugate). */
    Expr* denom = cancel_together(mk_fn2("Plus", mk_pow_int(expr_copy(d0), 2),
        mk_neg(mk_fn2("Times", expr_copy(rad), mk_pow_int(expr_copy(d1), 2)))));
    /* E = (n0 d0 - rad n1 d1)/denom. */
    Expr* E = cancel_together(mk_fn2("Times",
        mk_fn2("Plus", mk_fn2("Times", expr_copy(n0), expr_copy(d0)),
            mk_neg(mk_fn3("Times", expr_copy(rad), expr_copy(n1), expr_copy(d1)))),
        mk_inv(expr_copy(denom))));
    /* H = rad (n1 d0 - n0 d1)/denom. */
    Expr* H = cancel_together(mk_fn3("Times", expr_copy(rad),
        mk_fn2("Plus", mk_fn2("Times", expr_copy(n1), expr_copy(d0)),
            mk_neg(mk_fn2("Times", expr_copy(n0), expr_copy(d1)))),
        mk_inv(expr_copy(denom))));
    expr_free(n0); expr_free(n1); expr_free(d0); expr_free(d1); expr_free(denom);
    if (!E || !H) { if (E) expr_free(E); if (H) expr_free(H); return false; }
    *E_out = E; *H_out = H; return true;
}

/* I0 = Integrate[1/Sqrt[rad], x], rad = a x^2 + b x + c.  ArcSinh / ArcSin / Log
 * by the signs of a and the discriminant 4ac-b^2 (Log form for symbolic a).
 * Borrows a,b,c,rad,x; returns owned. */
static Expr* qr_I0(Expr* a, Expr* b, Expr* c, Expr* rad, Expr* x) {
    int sa = expr_numeric_sign(a);
    Expr* disc = eval_take(mk_fn2("Plus",
        mk_fn3("Times", mk_int(4), expr_copy(a), expr_copy(c)),
        mk_neg(mk_pow_int(expr_copy(b), 2))));         /* 4 a c - b^2 */
    int sd = disc ? expr_numeric_sign(disc) : 0;
    Expr* lin = mk_fn2("Plus", mk_fn3("Times", mk_int(2), expr_copy(a), expr_copy(x)),
                               expr_copy(b));           /* 2 a x + b */
    /* Form selection so the derivative stays Simplify-reducible AND the branch is
     * stable under the differentiate-back sampler:
     *  - numeric a with a PERFECT-SQUARE |disc|: the pretty ArcSinh/ArcSin;
     *  - otherwise (numeric a with non-square disc, or symbolic a): the Log form
     *    below, whose derivative cancels algebraically without a leftover
     *    Sqrt[rational*rad] that Simplify cannot pull apart, and which -- unlike
     *    ArcSinh[.../Sqrt[disc]] when disc<0 -- has no spurious branch jump at the
     *    numeric sample points the Goursat verifier probes.  (The symbolic-a
     *    pure-polynomial case, where the Log form's symbolic Simplify is costly,
     *    is declined in qr_direct so it routes to the Euler substitution.) */
    Expr* res = NULL;
    int64_t rn, rd;
    if (sa > 0 && sd > 0) {
        Expr* sq = eval_take(mk_sqrt(expr_copy(disc)));           /* Sqrt[4ac-b^2] */
        if (sq && is_rational(sq, &rn, &rd))
            res = mk_fn2("Times", mk_inv(mk_sqrt(expr_copy(a))),
                mk_fn1("ArcSinh", mk_fn2("Times", expr_copy(lin), mk_inv(sq))));
        else if (sq) expr_free(sq);
    } else if (sa < 0 && sd < 0) {
        Expr* sq = eval_take(mk_sqrt(mk_neg(expr_copy(disc))));   /* Sqrt[b^2-4ac] */
        if (sq && is_rational(sq, &rn, &rd))
            res = mk_neg(mk_fn2("Times", mk_inv(mk_sqrt(mk_neg(expr_copy(a)))),
                mk_fn1("ArcSin", mk_fn2("Times", expr_copy(lin), mk_inv(sq)))));
        else if (sq) expr_free(sq);
    }
    if (disc) expr_free(disc);
    if (!res) {   /* (1/Sqrt[a]) Log[2 a x + b + 2 Sqrt[a] Sqrt[rad]] */
        res = mk_fn2("Times", mk_inv(mk_sqrt(expr_copy(a))),
            mk_fn1("Log", mk_fn2("Plus", expr_copy(lin),
                mk_fn3("Times", mk_int(2), mk_sqrt(expr_copy(a)), mk_sqrt(expr_copy(rad))))));
    }
    expr_free(lin);
    return eval_take(res);
}

/* I_n = Integrate[x^n/Sqrt[rad], x], n >= 0, via the standard reduction
 *   I_n = x^{n-1} Sqrt[rad]/(n a) - ((2n-1)b)/(2 n a) I_{n-1}
 *                                 - ((n-1)c)/(n a)   I_{n-2}.
 * Borrows a,b,c,rad,x; returns owned. */
static Expr* qr_In(int n, Expr* a, Expr* b, Expr* c, Expr* rad, Expr* x) {
    if (n <= 0) return qr_I0(a, b, c, rad, x);
    Expr* Im1 = qr_In(n - 1, a, b, c, rad, x);
    Expr* na  = mk_fn2("Times", mk_int(n), expr_copy(a));
    Expr* term1 = mk_fn2("Times",
        mk_fn2("Times", mk_pow_int(expr_copy(x), n - 1), mk_sqrt(expr_copy(rad))),
        mk_inv(expr_copy(na)));
    Expr* term2 = mk_fn2("Times",
        mk_fn3("Times", mk_int(-(2 * n - 1)), expr_copy(b),
               mk_inv(mk_fn2("Times", mk_int(2), expr_copy(na)))), Im1);
    Expr* sum = mk_fn2("Plus", term1, term2);
    if (n >= 2) {
        Expr* Im2 = qr_In(n - 2, a, b, c, rad, x);
        sum = mk_fn2("Plus", sum, mk_fn2("Times",
            mk_fn3("Times", mk_int(-(n - 1)), expr_copy(c), mk_inv(expr_copy(na))), Im2));
    }
    expr_free(na);
    return eval_take(sum);
}

/* J(p) = Integrate[1/((x-p) Sqrt[rad]), x], A = rad(p):
 *   A>0 : -(1/Sqrt[A])  ArcTanh[(2A + (2 a p + b)(x-p))/(2 Sqrt[A]  Sqrt[rad])]
 *   A<0 :  (1/Sqrt[-A]) ArcTan [(2A + (2 a p + b)(x-p))/(2 Sqrt[-A] Sqrt[rad])]
 * A == 0 (pole at a branch point) -> NULL.  Borrows p,a,b,c,rad,x. */
static Expr* qr_Jp(Expr* p, Expr* a, Expr* b, Expr* c, Expr* rad, Expr* x) {
    Expr* A = eval_take(mk_fn3("Plus",
        mk_fn2("Times", expr_copy(a), mk_pow_int(expr_copy(p), 2)),
        mk_fn2("Times", expr_copy(b), expr_copy(p)), expr_copy(c)));
    if (!A) return NULL;
    int sA = qr_sign_probe(A);
    if (sA == 0) {
        Expr* z = eval_take(mk_fn1("PossibleZeroQ", expr_copy(A)));
        bool zero = z && z->type == EXPR_SYMBOL && z->data.symbol == SYM_True;
        if (z) expr_free(z);
        if (zero) { expr_free(A); return NULL; }   /* pole on a branch point */
        sA = 1;                                     /* indeterminate -> ArcTanh */
    }
    Expr* B   = mk_fn2("Plus", mk_fn3("Times", mk_int(2), expr_copy(a), expr_copy(p)),
                               expr_copy(b));                        /* 2 a p + b */
    Expr* xmp = mk_fn2("Plus", expr_copy(x), mk_neg(expr_copy(p))); /* x - p */
    Expr* numArg = mk_fn2("Plus", mk_fn2("Times", mk_int(2), expr_copy(A)),
                                  mk_fn2("Times", B, xmp));          /* 2A + B(x-p) */
    Expr* res;
    if (sA > 0) {
        Expr* sq = mk_sqrt(expr_copy(A));
        res = mk_neg(mk_fn2("Times", mk_inv(expr_copy(sq)),
            mk_fn1("ArcTanh", mk_fn2("Times", numArg,
                mk_inv(mk_fn3("Times", mk_int(2), sq, mk_sqrt(expr_copy(rad))))))));
    } else {
        Expr* sq = mk_sqrt(mk_neg(expr_copy(A)));
        res = mk_fn2("Times", mk_inv(expr_copy(sq)),
            mk_fn1("ArcTan", mk_fn2("Times", numArg,
                mk_inv(mk_fn3("Times", mk_int(2), sq, mk_sqrt(expr_copy(rad)))))));
    }
    expr_free(A);
    return eval_take(res);
}

/* Closed-form driver.  Borrows f,x,rad,a,b,c; returns owned antiderivative or
 * NULL (caller then tries the Euler substitution). */
static Expr* qr_direct(Expr* f, Expr* x, Expr* rad, Expr* a, Expr* b, Expr* c) {
    Expr* E = NULL, *H = NULL;
    if (!qr_split(f, rad, x, &E, &H)) return NULL;

    Expr* total = mk_int(0);
    bool okflag = true;
    bool has_pole = false;   /* did any 1/((x-p)Sqrt[rad]) term (-> ArcTan/ArcTanh) appear? */

    /* Rational (non-radical) part E. */
    if (!is_int_zero(E)) {
        if (expr_free_of(E, x)) {
            total = mk_fn2("Plus", total, mk_fn2("Times", expr_copy(E), expr_copy(x)));
        } else {
            qr_depth++;
            Expr* gE = integrate_in(E, x);
            qr_depth--;
            if (!gE) okflag = false;
            else total = mk_fn2("Plus", total, gE);
        }
    }
    expr_free(E);
    if (!okflag) { expr_free(H); expr_free(total); return NULL; }

    /* Radical part: partial-fraction H, classify each term. */
    Expr* Hp = eval_take(mk_fn2("Apart", H, expr_copy(x)));   /* consumes H */
    if (!Hp) { expr_free(total); return NULL; }

    Expr* poly_acc = mk_int(0);
    Expr* one_term[1] = { Hp };
    size_t nterms; Expr** terms;
    if (head_is(Hp, SYM_Plus)) { nterms = Hp->data.function.arg_count; terms = Hp->data.function.args; }
    else                       { nterms = 1; terms = one_term; }

    for (size_t i = 0; i < nterms && okflag; i++) {
        Expr* t = terms[i];
        Expr* dent = eval_take(mk_fn1("Expand",
                        eval_take(mk_fn1("Denominator", expr_copy(t)))));
        if (!dent) { okflag = false; break; }
        if (expr_free_of(dent, x)) {          /* polynomial-in-x term */
            poly_acc = mk_fn2("Plus", poly_acc, expr_copy(t));
            expr_free(dent);
            continue;
        }
        int deg = get_degree_poly(dent, x);
        if (deg != 1) { expr_free(dent); okflag = false; break; }   /* higher/irreducible pole */
        Expr* A1 = get_coeff(dent, x, 1);
        Expr* A0 = get_coeff(dent, x, 0);
        expr_free(dent);
        Expr* numt = eval_take(mk_fn1("Numerator", expr_copy(t)));
        if (!A1 || !A0 || !numt || !expr_free_of(numt, x)) {
            if (A1) expr_free(A1); if (A0) expr_free(A0); if (numt) expr_free(numt);
            okflag = false; break;
        }
        /* term = numt/(A1 x + A0) = (numt/A1)/(x - p), p = -A0/A1. */
        Expr* p     = cancel_together(mk_fn2("Times", mk_neg(A0), mk_inv(expr_copy(A1))));
        Expr* coeff = cancel_together(mk_fn2("Times", numt, mk_inv(A1)));
        Expr* J = p ? qr_Jp(p, a, b, c, rad, x) : NULL;
        if (p) expr_free(p);
        if (!J) { if (coeff) expr_free(coeff); okflag = false; break; }
        has_pole = true;
        total = mk_fn2("Plus", total, mk_fn2("Times", coeff, J));
    }
    expr_free(Hp);
    if (!okflag) { expr_free(poly_acc); expr_free(total); return NULL; }

    /* Polynomial part: sum_n Coefficient[poly,x,n] * I_n. */
    poly_acc = eval_take(mk_fn1("Expand", poly_acc));
    if (poly_acc && !is_int_zero(poly_acc)) {
        int pd = get_degree_poly(poly_acc, x);
        if (pd < 0) okflag = false;
        for (int n = 0; n <= pd && okflag; n++) {
            Expr* qn = get_coeff(poly_acc, x, n);
            if (!qn || is_int_zero(qn)) { if (qn) expr_free(qn); continue; }
            Expr* In = qr_In(n, a, b, c, rad, x);
            if (!In) { expr_free(qn); okflag = false; break; }
            total = mk_fn2("Plus", total, mk_fn2("Times", qn, In));
        }
    }
    if (poly_acc) expr_free(poly_acc);
    if (!okflag) { expr_free(total); return NULL; }

    /* A pure polynomial-over-Sqrt[rad] integrand with a SYMBOLIC leading
     * coefficient (no ArcTan/ArcTanh pole term) gains nothing here over the Euler
     * substitution -- and its Log-form I_n derivative makes a symbolic Simplify
     * (the differentiate-back check outside the parametric Goursat path) very
     * costly.  Decline so qr_core falls back to Euler; the parametric Goursat
     * genus-0 pieces always carry a pole, so they keep this closed form. */
    if (!has_pole && expr_numeric_sign(a) == 0) { expr_free(total); return NULL; }

    return eval_take(total);
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

    /* --- Preferred closed-form path (ArcTan/ArcTanh/ArcSinh/ArcSin/Log). ---
     * Emits clean antiderivatives directly; falls through to the Euler
     * substitution below when it cannot fully decompose the integrand. */
    result = qr_direct(f, x, rad, a, b, c);
    if (result) goto done;

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
