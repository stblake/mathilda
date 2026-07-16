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

Expr* integrate_fresnel_try(Expr* f, Expr* x) {
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
    const char* names[4] = { "frA", "frB", "frC", "frK" };
    Expr* vals[4] = { a, b, c, Kexpr };
    /* Build the answer via a template evaluated with the coefficients. */
    Expr* rules[4];
    for (int i = 0; i < 4; i++)
        rules[i] = expr_new_function(mk_sym("Rule"),
            (Expr*[]){ mk_sym(names[i]), expr_copy(vals[i]) }, 2);
    Expr* rulelist = expr_new_function(mk_sym("List"), rules, 4);

    /* u, d, IS, IC as symbolic templates in frA,frB,frC (substituted below). */
    const char* tmpl = (trig[0] == 'S')   /* Sin vs Cos */
        ? "frK (Cos[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelS[Sqrt[2 frA/Pi](x + frB/(2 frA))]"
          "   + Sin[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelC[Sqrt[2 frA/Pi](x + frB/(2 frA))])"
        : "frK (Cos[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelC[Sqrt[2 frA/Pi](x + frB/(2 frA))]"
          "   - Sin[frC - frB^2/(4 frA)] Sqrt[Pi/(2 frA)] FresnelS[Sqrt[2 frA/Pi](x + frB/(2 frA))])";
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
