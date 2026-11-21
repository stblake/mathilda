/* integrate_fresnel.c — Fresnel integrals of a Gaussian-phase trig integrand.
 *
 * Recognizes  K Sin[Q]  and  K Cos[Q]  with Q = a x^2 + b x + c a quadratic in x
 * (a != 0) and K free of x, and integrates by completing the square: with
 * u = x + b/(2a) and d = c - b^2/(4a),
 *   INT Sin[a u^2 + d] du = Cos[d] IS + Sin[d] IC,
 *   INT Cos[a u^2 + d] du = Cos[d] IC - Sin[d] IS,
 * where  IS = Sqrt[Pi/(2a)] FresnelS[Sqrt[2a/Pi] u],
 *        IC = Sqrt[Pi/(2a)] FresnelC[Sqrt[2a/Pi] u]
 * (Mathilda's FresnelS[z]' = Sin[Pi z^2/2], FresnelC[z]' = Cos[Pi z^2/2]).  This
 * is the trigonometric sibling of the K E^(a x^2+b x+c) -> Erf recognizer.  The
 * candidate is accepted only after an exact Simplify diff-back, so a
 * mis-recognition can never emit a wrong closed form.
 */

#include "integrate_fresnel.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stddef.h>

static Expr* mk_sym(const char* s)  { return expr_new_symbol(s); }
static Expr* mk_int(long n)         { return expr_new_integer(n); }
static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2);
}
/* evaluate() BORROWS its argument (does not free it), so free the input node. */
static Expr* evf(Expr* node)                        { Expr* r = evaluate(node);
                                                      expr_free(node); return r; }
static Expr* ev1(const char* h, Expr* a)            { return evf(mk_fn1(h, a)); }
static Expr* ev2(const char* h, Expr* a, Expr* b)   { return evf(mk_fn2(h, a, b)); }

static bool free_of(Expr* e, Expr* x) {
    Expr* r = ev2("FreeQ", expr_copy(e), expr_copy(x));
    bool y = r && r->type == EXPR_SYMBOL && r->data.symbol.name == intern_symbol("True");
    if (r) expr_free(r);
    return y;
}

/* Is e = Sin[Q] or Cos[Q] with Q a degree-2 polynomial in x?  On success returns
 * the head name ("Sin"/"Cos") and borrows *Q = the argument. */
static const char* trig_quad(Expr* e, Expr* x, Expr** Q) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return NULL;
    if (e->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != intern_symbol("Sin") && h != intern_symbol("Cos")) return NULL;
    Expr* arg = e->data.function.args[0];
    Expr* pq = ev2("PolynomialQ", expr_copy(arg), expr_copy(x));
    bool poly = pq && pq->type == EXPR_SYMBOL && pq->data.symbol.name == intern_symbol("True");
    if (pq) expr_free(pq);
    if (!poly) return NULL;
    *Q = arg;
    return h;
}

/* Is z the exact integer 0? */
static bool is_zero(Expr* z) {
    return z && z->type == EXPR_INTEGER && z->data.integer == 0;
}

static Expr* integrate_fresnel_kernel(Expr* f, Expr* x);

/* ---- the power family (v0.6662) -------------------------------------------
 *
 * Two shapes reduce to the kernel above by parts:
 *   (P) K u^(2m) Trig[a u^2 + c]       (m a nonzero integer, a != 0)
 *   (H) K x^p    Trig[a x + c]         (p a half-integer; x = u^2 gives
 *                                       2 K u^(2p+1) Trig[a u^2 + c], an (P))
 * With S = Sin[a u^2 + c], C = Cos[a u^2 + c] and J_m^T = INT u^(2m) T du:
 *   m > 0:  J_m^S = -u^(2m-1) C/(2a) + (2m-1)/(2a) J_(m-1)^C
 *           J_m^C =  u^(2m-1) S/(2a) - (2m-1)/(2a) J_(m-1)^S
 *   m < 0:  J_m^S =  u^(2m+1) S/(2m+1) - 2a/(2m+1) J_(m+1)^C
 *           J_m^C =  u^(2m+1) C/(2m+1) + 2a/(2m+1) J_(m+1)^S
 * down (or up) to J_0, the Fresnel kernel.  For (H) the answer is put back
 * with u = Sqrt[x], exact on the principal branch (x^(k+1/2) == Sqrt[x]^(2k+1)).
 * Before this the whole family -- Sqrt[t] Sin[t], t^(3/2) Cos[3 t],
 * Sin[t]/Sqrt[t] -- fell through to ParallelMixedSpecial, which spent ~11 s
 * declining each one before RoachMeijerG closed it, and DSolve's separable
 * 2.2.23-2299 (Sqrt[t] Sin[t] y + y' == 0) ran out of its quadrature budget.
 * The result passes the same exact diff-back gate as the kernel's. */
#define FRESNEL_MAX_POWER 12

/* J_m of u^(2m) Trig[a u^2 + c]; trig is "Sin" or "Cos".  Owned result. */
static Expr* fresnel_power(int m, const char* trig, Expr* a, Expr* c, Expr* u) {
    if (m == 0) {
        Expr* arg = ev2("Plus", ev2("Times", expr_copy(a), mk_fn2("Power", expr_copy(u), mk_int(2))),
                        expr_copy(c));
        Expr* k = mk_fn1(trig, arg);
        Expr* r = integrate_fresnel_kernel(k, u);
        expr_free(k);
        return r;
    }
    bool is_sin = trig[0] == 'S';
    Expr* arg = ev2("Plus", ev2("Times", expr_copy(a), mk_fn2("Power", expr_copy(u), mk_int(2))),
                    expr_copy(c));
    const char* other = is_sin ? "Cos" : "Sin";
    if (m > 0) {
        Expr* rest = fresnel_power(m - 1, other, a, c, u);
        if (!rest) { expr_free(arg); return NULL; }
        /* sign: Sin -> -u^(2m-1) C/(2a) + (2m-1)/(2a) J^C ; Cos -> +u^(2m-1) S/(2a) - ... */
        Expr* boundary = ev2("Times", mk_fn2("Power", expr_copy(u), mk_int(2 * m - 1)),
                             ev2("Times", mk_fn1(other, expr_copy(arg)),
                                 mk_fn2("Power", ev2("Times", mk_int(2), expr_copy(a)), mk_int(-1))));
        Expr* coef = ev2("Times", mk_int(2 * m - 1),
                         mk_fn2("Power", ev2("Times", mk_int(2), expr_copy(a)), mk_int(-1)));
        Expr* sgn = mk_int(is_sin ? -1 : 1);
        Expr* r = ev2("Plus", ev2("Times", expr_copy(sgn), boundary),
                      ev2("Times", ev2("Times", ev1("Minus", sgn), coef), rest));
        expr_free(arg);
        return r;
    }
    /* m < 0 */
    Expr* rest = fresnel_power(m + 1, other, a, c, u);
    if (!rest) { expr_free(arg); return NULL; }
    Expr* boundary = ev2("Times", mk_fn2("Power", expr_copy(u), mk_int(2 * m + 1)),
                         ev2("Times", mk_fn1(trig, expr_copy(arg)),
                             mk_fn2("Power", mk_int(2 * m + 1), mk_int(-1))));
    Expr* coef = ev2("Times", ev2("Times", mk_int(2), expr_copy(a)),
                     mk_fn2("Power", mk_int(2 * m + 1), mk_int(-1)));
    Expr* r = ev2("Plus", boundary,
                  ev2("Times", ev2("Times", mk_int(is_sin ? -1 : 1), coef), rest));
    expr_free(arg);
    return r;
}

/* Is e the exact integer 0? (forward: defined below) */
static bool is_zero(Expr* z);

/* Split f = K * x^p * Trig[L]: K free of x, p rational (1 if absent), Trig a
 * Sin/Cos whose argument L contains x.  Borrows *trig_e; *K, *p owned. */
static bool fresnel_split(Expr* f, Expr* x, Expr** K, Expr** p, Expr** trig_e) {
    Expr* const* fs = &f;
    size_t nf = 1;
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == intern_symbol("Times")) {
        fs = f->data.function.args; nf = f->data.function.arg_count;
    }
    Expr* tr = NULL; Expr* pw = NULL; Expr* kk = mk_int(1);
    for (size_t i = 0; i < nf; i++) {
        Expr* g = fs[i];
        if (free_of(g, x)) { kk = ev2("Times", kk, expr_copy(g)); continue; }
        bool is_trig = g->type == EXPR_FUNCTION && g->data.function.arg_count == 1
            && g->data.function.head->type == EXPR_SYMBOL
            && (g->data.function.head->data.symbol.name == intern_symbol("Sin")
                || g->data.function.head->data.symbol.name == intern_symbol("Cos"));
        if (is_trig && !tr) { tr = g; continue; }
        /* x itself, or x^q with q rational */
        if (!pw && g->type == EXPR_SYMBOL && expr_eq(g, x)) { pw = mk_int(1); continue; }
        if (!pw && g->type == EXPR_FUNCTION && g->data.function.arg_count == 2
            && g->data.function.head->type == EXPR_SYMBOL
            && g->data.function.head->data.symbol.name == intern_symbol("Power")
            && expr_eq(g->data.function.args[0], x)) {
            Expr* q = g->data.function.args[1];
            bool rat = q->type == EXPR_INTEGER || (q->type == EXPR_FUNCTION
                && q->data.function.head->type == EXPR_SYMBOL
                && q->data.function.head->data.symbol.name == intern_symbol("Rational"));
            if (rat) { pw = expr_copy(q); continue; }
        }
        expr_free(kk); if (pw) expr_free(pw);
        return false;
    }
    if (!tr) { expr_free(kk); if (pw) expr_free(pw); return false; }
    *K = kk; *p = pw ? pw : mk_int(0); *trig_e = tr;
    return true;
}

/* Coefficients of L = a x^d + c (exactly, a != 0, no other x terms). */
static bool fresnel_coeffs(Expr* L, Expr* x, int d, Expr** a, Expr** c) {
    Expr* A = evf(expr_new_function(mk_sym("Coefficient"),
                  (Expr*[]){ expr_copy(L), expr_copy(x), mk_int(d) }, 3));
    Expr* C0 = evf(expr_new_function(mk_sym("Coefficient"),
                  (Expr*[]){ expr_copy(L), expr_copy(x), mk_int(0) }, 3));
    Expr* recon = ev2("Plus", ev2("Times", expr_copy(A), mk_fn2("Power", expr_copy(x), mk_int(d))),
                      expr_copy(C0));
    Expr* rem = ev1("Expand", ev2("Plus", expr_copy(L), ev2("Times", mk_int(-1), recon)));
    bool ok = A && !is_zero(A) && free_of(A, x) && free_of(C0, x) && is_zero(rem);
    if (rem) expr_free(rem);
    if (!ok) { expr_free(A); expr_free(C0); return false; }
    *a = A; *c = C0;
    return true;
}

/* The power family above; NULL when f is not one of its shapes. */
static Expr* integrate_fresnel_power(Expr* f, Expr* x) {
    Expr *K = NULL, *p = NULL, *tr = NULL;
    if (!fresnel_split(f, x, &K, &p, &tr)) return NULL;
    const char* trig = tr->data.function.head->data.symbol.name == intern_symbol("Sin") ? "Sin" : "Cos";
    Expr* L = tr->data.function.args[0];
    Expr *a = NULL, *c = NULL, *res = NULL;
    /* p = num/den */
    int64_t num = 0, den = 1;
    if (p->type == EXPR_INTEGER) num = p->data.integer;
    else if (p->data.function.args[0]->type == EXPR_INTEGER && p->data.function.args[1]->type == EXPR_INTEGER) {
        num = p->data.function.args[0]->data.integer; den = p->data.function.args[1]->data.integer;
    } else den = 0;
    if (den == 1 && num != 0 && num % 2 == 0 && num / 2 >= -FRESNEL_MAX_POWER
        && num / 2 <= FRESNEL_MAX_POWER && fresnel_coeffs(L, x, 2, &a, &c)) {
        /* (P): x^(2m) Trig[a x^2 + c] */
        Expr* J = fresnel_power((int)(num / 2), trig, a, c, x);
        if (J) res = ev2("Times", expr_copy(K), J);
    } else if (den == 2 && fresnel_coeffs(L, x, 1, &a, &c)) {
        /* (H): x^(num/2) Trig[a x + c]; u^(2p+1) = u^(num+1), num odd */
        int64_t m2 = num + 1;                      /* even */
        if (m2 / 2 >= -FRESNEL_MAX_POWER && m2 / 2 <= FRESNEL_MAX_POWER) {
            Expr* u = mk_sym("Integrate`fresnelU");
            Expr* J = fresnel_power((int)(m2 / 2), trig, a, c, u);
            if (J) {
                Expr* back = evf(expr_new_function(mk_sym("ReplaceAll"),
                    (Expr*[]){ J, expr_new_function(mk_sym("Rule"),
                        (Expr*[]){ expr_copy(u), mk_fn1("Sqrt", expr_copy(x)) }, 2) }, 2));
                res = ev2("Times", ev2("Times", mk_int(2), expr_copy(K)), back);
            }
            expr_free(u);
        }
    }
    if (a) expr_free(a);
    if (c) expr_free(c);
    expr_free(K); expr_free(p);
    if (!res) return NULL;
    /* the exact diff-back gate, as for the kernel */
    Expr* chk = ev1("Simplify", ev2("Plus", ev2("D", expr_copy(res), expr_copy(x)),
                                    ev2("Times", mk_int(-1), expr_copy(f))));
    bool ok = chk && chk->type == EXPR_INTEGER && chk->data.integer == 0;
    if (chk) expr_free(chk);
    if (!ok) { expr_free(res); return NULL; }
    return res;
}

Expr* integrate_fresnel_try(Expr* f, Expr* x) {
    Expr* r = integrate_fresnel_kernel(f, x);
    if (r) return r;
    return integrate_fresnel_power(f, x);
}

static Expr* integrate_fresnel_kernel(Expr* f, Expr* x) {
    /* Locate the Sin[Q]/Cos[Q] kernel and the x-free prefactor K. */
    const char* trig = NULL;
    Expr* Q = NULL;                 /* borrowed */
    Expr* Kexpr = NULL;             /* owned: the constant prefactor */

    if ((trig = trig_quad(f, x, &Q))) {
        Kexpr = mk_int(1);
    } else if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
               && f->data.function.head->data.symbol.name == intern_symbol("Times")) {
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            Expr* fac = f->data.function.args[i];
            const char* h = trig_quad(fac, x, &Q);
            if (!h) continue;
            /* K = f / fac  (must be free of x). */
            Expr* k = ev2("Times", expr_copy(f),
                          expr_new_function(mk_sym("Power"),
                              (Expr*[]){ expr_copy(fac), mk_int(-1) }, 2));
            if (free_of(k, x)) { trig = h; Kexpr = k; break; }
            expr_free(k);
            Q = NULL;
        }
    }
    if (!trig || !Kexpr) { if (Kexpr) expr_free(Kexpr); return NULL; }

    /* Q = a x^2 + b x + c  via Coefficient[Q, x, k]. */
    Expr* a = evf(expr_new_function(mk_sym("Coefficient"),
            (Expr*[]){ expr_copy(Q), expr_copy(x), mk_int(2) }, 3));
    Expr* b = evf(expr_new_function(mk_sym("Coefficient"),
            (Expr*[]){ expr_copy(Q), expr_copy(x), mk_int(1) }, 3));
    Expr* c = evf(expr_new_function(mk_sym("Coefficient"),
            (Expr*[]){ expr_copy(Q), expr_copy(x), mk_int(0) }, 3));

    /* Require degree EXACTLY 2: a != 0 and Q == a x^2 + b x + c. */
    Expr* recon = expr_new_function(mk_sym("Plus"),
        (Expr*[]){ expr_new_function(mk_sym("Times"),
                     (Expr*[]){ expr_copy(a), mk_fn2("Power", expr_copy(x), mk_int(2)) }, 2),
                   expr_new_function(mk_sym("Times"),
                     (Expr*[]){ expr_copy(b), expr_copy(x) }, 2),
                   expr_copy(c) }, 3);
    Expr* rem = ev1("Expand", expr_new_function(mk_sym("Plus"),
        (Expr*[]){ expr_copy(Q), expr_new_function(mk_sym("Times"),
                     (Expr*[]){ mk_int(-1), recon }, 2) }, 2));   /* Q - recon */
    bool quad = is_zero(a) ? false : is_zero(rem);
    if (rem) expr_free(rem);
    if (!quad) {
        expr_free(a); expr_free(b); expr_free(c); expr_free(Kexpr);
        return NULL;
    }

    /* h = -b/(2a),  u = x - h = x + b/(2a),  d = c - b^2/(4a). */
    /* frX stands for the integration variable: the template used to spell it
     * as the literal symbol x, so for any other variable the candidate was in
     * x, its t-derivative 0, the diff-back gate rejected it and the stage
     * declined -- Integrate[Cos[t^2], t] was Integrate::nonelem while
     * Integrate[Cos[x^2], x] was FresnelC (since v0.15x, July 2026). */
    const char* names[5] = { "frA", "frB", "frC", "frK", "frX" };
    Expr* vals[5] = { a, b, c, Kexpr, x };
    /* Build the answer via a template evaluated with the coefficients. */
    Expr* rules[5];
    for (int i = 0; i < 5; i++)
        rules[i] = expr_new_function(mk_sym("Rule"),
            (Expr*[]){ mk_sym(names[i]), expr_copy(vals[i]) }, 2);
    Expr* rulelist = expr_new_function(mk_sym("List"), rules, 5);

    /* u, d, IS, IC as symbolic templates in frA,frB,frC (substituted below). */
    const char* tmpl = (trig[0] == 'S')   /* Sin vs Cos */
        ? "frK (Cos[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelS[Sqrt[2 frA/Pi](frX + frB/(2 frA))]"
          "   + Sin[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelC[Sqrt[2 frA/Pi](frX + frB/(2 frA))])"
        : "frK (Cos[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelC[Sqrt[2 frA/Pi](frX + frB/(2 frA))]"
          "   - Sin[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelS[Sqrt[2 frA/Pi](frX + frB/(2 frA))])";
    Expr* result = NULL;
    Expr* parsed = parse_expression(tmpl);
    if (!parsed) {
        expr_free(rulelist);
    } else {
        Expr* subst = evf(expr_new_function(mk_sym("ReplaceAll"),
            (Expr*[]){ parsed, rulelist }, 2));   /* consumes parsed, rulelist */
        if (subst) {
            /* exact diff-back gate */
            Expr* diff = expr_new_function(mk_sym("Plus"),
                (Expr*[]){ ev2("D", expr_copy(subst), expr_copy(x)),
                           expr_new_function(mk_sym("Times"),
                               (Expr*[]){ mk_int(-1), expr_copy(f) }, 2) }, 2);
            Expr* chk = ev1("Simplify", diff);
            if (chk && chk->type == EXPR_INTEGER && chk->data.integer == 0)
                result = subst;
            else if (subst)
                expr_free(subst);
            if (chk) expr_free(chk);
        }
    }

    expr_free(a); expr_free(b); expr_free(c); expr_free(Kexpr);
    return result;
}
