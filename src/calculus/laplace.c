/* laplace.c — LaplaceTransform[f, t, s]: the one-sided Laplace transform
 *     F(s) = Int_0^Infinity f(t) E^(-s t) dt.
 *
 * Algorithm (each step exact; the integral is the only analysis).  The
 * factor classification is native C; Expand, Integrate and Simplify are
 * called as builtins:
 *   0. Cosh[u] and Sinh[u], u linear in t, are written as exponentials, so
 *      they are shifts of step 2; the transform is recombined by Simplify.
 *   1. Lists thread; sums are transformed term by term (linearity); factors
 *      free of t come out, also of a transform that stays held.
 *   2. Exponential shift: E^(a t + b) g(t) -> E^b G(s - a), G = L[g].  This
 *      keeps the convergence region of the remaining integral at Re s > 0
 *      whatever the exponential rate a.
 *   3. A t-free remainder is c/s; a derivative f^(n)(t) of an undefined f is
 *      s^n F(s) - s^(n-1) f(0) - ... - f^(n-1)(0), F held.
 *   4. Otherwise the defining integral, under the user's assumptions and
 *      s > 0 (the abscissa of convergence of the shifted g is 0 for every
 *      transform of exponential order 0 that Integrate closes).  A result
 *      ConditionalExpression[F, cond] (Sin[b t]: s > Abs[Im[b]]) is F,
 *      simplified under cond: the transform is the analytic function the
 *      integral defines on its region of convergence.  If Integrate leaves
 *      the integral unevaluated, so is the transform.
 * LaplaceTransform[BesselJ[0, t], t, s] = 1/Sqrt[1 + s^2] (Lipschitz),
 * L[t^2 E^(3 t)] = 2/(s - 3)^3, L[Sin[a t]] = a/(a^2 + s^2).
 */
#include "eval.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "simp_sf_engine.h" /* sf_contains_head */
#include "sym_names.h"
#include "common.h"     /* head_is */
#include <stdlib.h>
#include <string.h>

static Expr* lt_fn(const char* h, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(h), a, n);
}
static Expr* lt_fn2(const char* h, Expr* x, Expr* y) { return lt_fn(h, (Expr*[]){ x, y }, 2); }
static Expr* lt_eval(Expr* e) { Expr* r = evaluate(e); expr_free(e); return r; }

static bool lt_has(const Expr* e, const Expr* t) {
    if (expr_eq(e, t)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (lt_has(e->data.function.head, t)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lt_has(e->data.function.args[i], t)) return true;
    return false;
}

/* p (evaluated, Expand-ed) == p0 + p1 t with p0, p1 free of t and p1 != 0. */
static bool lt_linear(const Expr* p, const Expr* t, Expr** p0, Expr** p1) {
    size_t nt = head_is(p, SYM_Plus) ? p->data.function.arg_count : 1;
    Expr* const* tm = head_is(p, SYM_Plus) ? p->data.function.args : (Expr* const*)&p;
    Expr** a = malloc(sizeof(Expr*) * nt);
    Expr** b = malloc(sizeof(Expr*) * nt);
    size_t na = 0, nb = 0;
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* x = tm[i];
        if (!lt_has(x, t)) { b[nb++] = expr_copy((Expr*)x); continue; }
        if (expr_eq(x, t)) { a[na++] = expr_new_integer(1); continue; }
        if (!head_is(x, SYM_Times)) { ok = false; break; }
        size_t nf = x->data.function.arg_count, nr = 0, nv = 0;
        Expr** r = malloc(sizeof(Expr*) * nf);
        for (size_t j = 0; j < nf; j++) {
            const Expr* y = x->data.function.args[j];
            if (expr_eq(y, t)) nv++;
            else if (lt_has(y, t)) nv += 2;
            else r[nr++] = expr_copy((Expr*)y);
        }
        if (nv == 1) a[na++] = lt_fn(SYM_Times, r, nr);
        else { for (size_t j = 0; j < nr; j++) expr_free(r[j]); ok = false; }
        free(r);
    }
    if (ok && na > 0) {
        *p1 = lt_eval(lt_fn(SYM_Plus, a, na));
        *p0 = lt_eval(lt_fn(SYM_Plus, b, nb));
    } else {
        for (size_t i = 0; i < na; i++) expr_free(a[i]);
        for (size_t i = 0; i < nb; i++) expr_free(b[i]);
        ok = false;
    }
    free(a); free(b);
    return ok;
}

static Expr* lt_transform(const Expr* f, const Expr* t, const Expr* s, const Expr* as);

/* The transform of g, or LaplaceTransform[g, t, s, Assumptions -> as] held. */
static Expr* lt_or_held(const Expr* g, const Expr* t, const Expr* s, const Expr* as) {
    Expr* v = lt_transform(g, t, s, as);
    if (v) return v;
    /* c h(t) with c free of t: c LaplaceTransform[h] (2 f[t] -> 2 L[f[t]]) */
    if (head_is(g, SYM_Times)) {
        size_t n = g->data.function.arg_count, nc = 0, nr = 0;
        Expr** cf = malloc(sizeof(Expr*) * n);
        Expr** rf = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            Expr* x = g->data.function.args[i];
            if (lt_has(x, t)) rf[nr++] = expr_copy(x); else cf[nc++] = expr_copy(x);
        }
        if (nc > 0 && nr > 0) {
            Expr* h = lt_eval(lt_fn(SYM_Times, rf, nr));
            Expr* held = lt_or_held(h, t, s, as);
            expr_free(h);
            free(rf);
            Expr* c = lt_fn(SYM_Times, cf, nc);
            free(cf);
            return lt_eval(lt_fn2(SYM_Times, c, held));
        }
        for (size_t i = 0; i < nc; i++) expr_free(cf[i]);
        for (size_t i = 0; i < nr; i++) expr_free(rf[i]);
        free(cf); free(rf);
    }
    if (as->type == EXPR_SYMBOL && as->data.symbol.name == SYM_True)
        return lt_fn(SYM_LaplaceTransform, (Expr*[]){ expr_copy((Expr*)g), expr_copy((Expr*)t),
                                                      expr_copy((Expr*)s) }, 3);
    return lt_fn(SYM_LaplaceTransform, (Expr*[]){ expr_copy((Expr*)g), expr_copy((Expr*)t),
        expr_copy((Expr*)s), lt_fn2(SYM_Rule, expr_new_symbol(SYM_Assumptions),
                                    expr_copy((Expr*)as)) }, 4);
}

/* The transform of f^(n)(t), f an undefined function (Derivative[n][f][t]):
 *     s^n F(s) - Sum[s^(n - 1 - k) f^(k)(0), {k, 0, n - 1}],
 * F = LaplaceTransform[f[t], t, s] held (from integrating by parts n times;
 * the boundary terms at Infinity vanish on the region of convergence).  NULL
 * when r is not of that form. */
static Expr* lt_derivative(const Expr* r, const Expr* t, const Expr* s) {
    if (r->type != EXPR_FUNCTION || r->data.function.arg_count != 1
        || !expr_eq(r->data.function.args[0], t))
        return NULL;
    const Expr* dh = r->data.function.head;              /* Derivative[n][f] */
    if (dh->type != EXPR_FUNCTION || dh->data.function.arg_count != 1
        || dh->data.function.head->type != EXPR_FUNCTION)
        return NULL;
    const Expr* d = dh->data.function.head;              /* Derivative[n] */
    const Expr* fn = dh->data.function.args[0];
    if (!head_is(d, SYM_Derivative) || d->data.function.arg_count != 1
        || d->data.function.args[0]->type != EXPR_INTEGER
        || fn->type != EXPR_SYMBOL || lt_has(fn, t))
        return NULL;
    int64_t n = d->data.function.args[0]->data.integer;
    if (n < 1 || n > 64) return NULL;
    Expr* ft = lt_fn(fn->data.symbol.name, (Expr*[]){ expr_copy((Expr*)t) }, 1);
    Expr* F = lt_fn(SYM_LaplaceTransform, (Expr*[]){ ft, expr_copy((Expr*)t),
                                                     expr_copy((Expr*)s) }, 3);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(n + 1));
    terms[0] = lt_fn2(SYM_Times, lt_fn2(SYM_Power, expr_copy((Expr*)s), expr_new_integer(n)), F);
    for (int64_t k = 0; k < n; k++) {
        Expr* fk = k == 0 ? expr_new_symbol(fn->data.symbol.name)
                          : expr_new_function(lt_fn(SYM_Derivative,
                                (Expr*[]){ expr_new_integer(k) }, 1),
                                (Expr*[]){ expr_new_symbol(fn->data.symbol.name) }, 1);
        Expr* f0 = expr_new_function(fk, (Expr*[]){ expr_new_integer(0) }, 1);
        terms[k + 1] = lt_fn(SYM_Times, (Expr*[]){ expr_new_integer(-1),
            lt_fn2(SYM_Power, expr_copy((Expr*)s), expr_new_integer(n - 1 - k)), f0 }, 3);
    }
    Expr* out = lt_eval(lt_fn(SYM_Plus, terms, (size_t)(n + 1)));
    free(terms);
    return out;
}

/* ConditionalExpression[F, cond] -> F, simplified under cond (where it was
 * derived: s/Sqrt[(b^2 + s^2)^2] under b > 0 is s/(b^2 + s^2)); anything
 * else is returned as is.  Consumes v. */
static Expr* lt_generic_value(Expr* v) {
    if (v && head_is(v, SYM_ConditionalExpression) && v->data.function.arg_count == 2) {
        Expr* r = lt_eval(lt_fn2(SYM_Simplify, expr_copy(v->data.function.args[0]),
                                 expr_copy(v->data.function.args[1])));
        expr_free(v);
        return r;
    }
    return v;
}

/* e with every Cosh[u] / Sinh[u], u linear in t, written as (E^u +- E^-u)/2:
 * the exponentials are then shifts (step 2), so L[Cosh[a t] g] =
 * (G(s - a) + G(s + a))/2 for any g Integrate closes on its own -- the
 * defining integral of Cosh[2 t] E^(-s t) itself did not close.  *hit is set
 * when anything was rewritten.  Returns a new expression. */
static Expr* lt_hyper_to_exp(const Expr* e, const Expr* t, bool* hit) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if ((head_is(e, SYM_Cosh) || head_is(e, SYM_Sinh)) && e->data.function.arg_count == 1
        && lt_has(e->data.function.args[0], t)) {
        Expr* u = lt_eval(lt_fn(SYM_Expand, (Expr*[]){ expr_copy(e->data.function.args[0]) }, 1));
        Expr *u0 = NULL, *u1 = NULL;
        if (u && lt_linear(u, t, &u0, &u1)) {
            expr_free(u0); expr_free(u1);
            *hit = true;
            Expr* ep = lt_fn2(SYM_Power, expr_new_symbol(SYM_E), expr_copy(u));
            Expr* em = lt_fn2(SYM_Power, expr_new_symbol(SYM_E),
                              lt_fn2(SYM_Times, expr_new_integer(-1), u));
            if (head_is(e, SYM_Sinh)) em = lt_fn2(SYM_Times, expr_new_integer(-1), em);
            return lt_fn2(SYM_Times, lt_fn2(SYM_Power, expr_new_integer(2), expr_new_integer(-1)),
                          lt_fn2(SYM_Plus, ep, em));
        }
        if (u) expr_free(u);
        return expr_copy((Expr*)e);
    }
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = lt_hyper_to_exp(e->data.function.args[i], t, hit);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return out;
}

static Expr* lt_transform(const Expr* f, const Expr* t, const Expr* s, const Expr* as) {
    bool hyper = false;
    Expr* fx = lt_hyper_to_exp(f, t, &hyper);
    if (hyper) {
        /* the transform of the exponential form, recombined over s */
        Expr* v = lt_transform(fx, t, s, as);
        expr_free(fx);
        if (v && !sf_contains_head(v, SYM_LaplaceTransform))
            v = lt_eval(lt_fn(SYM_Simplify, (Expr*[]){ v }, 1));
        return v;
    }
    expr_free(fx);
    Expr* g = lt_eval(lt_fn(SYM_Expand, (Expr*[]){ expr_copy((Expr*)f) }, 1));
    if (!g) return NULL;
    Expr* out = NULL;
    if (head_is(g, SYM_List) || head_is(g, SYM_Plus)) {
        /* 1. linearity */
        size_t n = g->data.function.arg_count;
        Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) v[i] = lt_or_held(g->data.function.args[i], t, s, as);
        out = lt_eval(lt_fn(head_is(g, SYM_List) ? SYM_List : SYM_Plus, v, n));
        free(v);
        expr_free(g);
        return out;
    }
    size_t nf = head_is(g, SYM_Times) ? g->data.function.arg_count : 1;
    Expr* const* fac = head_is(g, SYM_Times) ? g->data.function.args : &g;
    Expr** cf = malloc(sizeof(Expr*) * (nf + 1));
    Expr** rf = malloc(sizeof(Expr*) * (nf + 1));
    Expr** pf = malloc(sizeof(Expr*) * (nf + 1));
    size_t nc = 0, nr = 0, np = 0;
    for (size_t i = 0; i < nf; i++) {
        const Expr* x = fac[i];
        Expr *q0 = NULL, *q1 = NULL;
        if (!lt_has(x, t)) { cf[nc++] = expr_copy((Expr*)x); continue; }
        if (head_is(x, SYM_Power) && x->data.function.arg_count == 2
            && x->data.function.args[0]->type == EXPR_SYMBOL
            && x->data.function.args[0]->data.symbol.name == SYM_E
            && lt_linear(x->data.function.args[1], t, &q0, &q1)) {
            expr_free(q0); expr_free(q1);
            pf[np++] = expr_copy(x->data.function.args[1]);
            continue;
        }
        rf[nr++] = expr_copy((Expr*)x);
    }
    expr_free(g);
    Expr* c = lt_eval(lt_fn(SYM_Times, cf, nc));
    Expr* r = lt_eval(lt_fn(SYM_Times, rf, nr));
    Expr* p = lt_eval(lt_fn(SYM_Expand, (Expr*[]){ lt_fn(SYM_Plus, pf, np) }, 1));
    free(cf); free(rf); free(pf);
    Expr *p0 = NULL, *p1 = NULL, *v_der = NULL;
    if (np > 0 && lt_linear(p, t, &p0, &p1)) {
        /* 2. exponential shift: E^(p1 t + p0) g -> E^p0 G(s - p1) */
        Expr* v = lt_transform(r, t, s, as);
        if (v && !sf_contains_head(v, SYM_LaplaceTransform)) {
            Expr* sh = lt_fn2(SYM_Plus, expr_copy((Expr*)s),
                              lt_fn2(SYM_Times, expr_new_integer(-1), expr_copy(p1)));
            Expr* vs = lt_eval(lt_fn2(SYM_ReplaceAll, v, lt_fn2(SYM_Rule, expr_copy((Expr*)s), sh)));
            out = lt_eval(lt_fn(SYM_Times, (Expr*[]){ expr_copy(c),
                lt_fn2(SYM_Power, expr_new_symbol(SYM_E), expr_copy(p0)), vs }, 3));
        } else if (v) expr_free(v);
        expr_free(p0); expr_free(p1);
    } else if ((v_der = lt_derivative(r, t, s)) != NULL) {
        /* 3a. a derivative of an undefined function */
        out = lt_eval(lt_fn2(SYM_Times, expr_copy(c), v_der));
    } else if (r->type == EXPR_INTEGER && r->data.integer == 1) {
        /* 3. a t-free remainder */
        out = lt_eval(lt_fn2(SYM_Times, expr_copy(c),
                             lt_fn2(SYM_Power, expr_copy((Expr*)s), expr_new_integer(-1))));
    } else {
        /* 4. the defining integral under s > 0 */
        Expr* ker = lt_fn2(SYM_Times, expr_copy(r), lt_fn2(SYM_Power, expr_new_symbol(SYM_E),
                     lt_fn(SYM_Times, (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)s),
                                                 expr_copy((Expr*)t) }, 3)));
        Expr* rng = lt_fn(SYM_List, (Expr*[]){ expr_copy((Expr*)t), expr_new_integer(0),
                                                expr_new_symbol(SYM_Infinity) }, 3);
        Expr* asm_ = lt_fn2(SYM_And, expr_copy((Expr*)as),
                            lt_fn2(SYM_Greater, expr_copy((Expr*)s), expr_new_integer(0)));
        Expr* v = lt_eval(lt_fn(SYM_Integrate, (Expr*[]){ ker, rng,
                     lt_fn2(SYM_Rule, expr_new_symbol(SYM_Assumptions), asm_) }, 3));
        /* A ConditionalExpression[F, cond] is F where the integral converges
         * (Sin[b t]: s > Abs[Im[b]]; t^n: n > -1).  The transform is the
         * analytic function F, defined by the integral on its region of
         * convergence, so F is the transform: its generic value. */
        v = lt_generic_value(v);
        if (v && !sf_contains_head(v, SYM_Integrate)
              && !sf_contains_head(v, SYM_ConditionalExpression))
            out = lt_eval(lt_fn2(SYM_Times, expr_copy(c),
                                 lt_fn(SYM_Simplify, (Expr*[]){ v }, 1)));
        else if (v) expr_free(v);
    }
    expr_free(c); expr_free(r); expr_free(p);
    return out;
}

/* LaplaceTransform[f, t, s] or [f, t, s, Assumptions -> a]. */
static Expr* builtin_laplace_transform(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n != 3 && n != 4) return NULL;
    Expr* t = res->data.function.args[1];
    Expr* s = res->data.function.args[2];
    if (t->type != EXPR_SYMBOL || s->type != EXPR_SYMBOL) return NULL;
    Expr* as = NULL;
    if (n == 4) {
        Expr* o = res->data.function.args[3];
        if (!(head_is(o, SYM_Rule) && o->data.function.arg_count == 2
              && o->data.function.args[0]->type == EXPR_SYMBOL
              && o->data.function.args[0]->data.symbol.name == SYM_Assumptions))
            return NULL;
        as = o->data.function.args[1];
    }
    Expr* tru = as ? NULL : expr_new_symbol(SYM_True);
    Expr* r = lt_transform(res->data.function.args[0], t, s, as ? as : tru);
    if (tru) expr_free(tru);
    return r;
}

void laplace_init(void) {
    symtab_add_builtin(SYM_LaplaceTransform, builtin_laplace_transform);
    symtab_get_def(SYM_LaplaceTransform)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_LaplaceTransform,
        "LaplaceTransform[f, t, s] gives the one-sided Laplace transform "
        "Integrate[f E^(-s t), {t, 0, Infinity}] of f as a function of s. "
        "Linear; threads over lists; an exponential factor E^(a t) shifts s -> s - a. "
        "LaplaceTransform[f, t, s, Assumptions -> a] uses the assumptions a.");
}
