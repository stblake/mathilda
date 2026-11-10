/* laplace.c — LaplaceTransform[f, t, s]: the one-sided Laplace transform
 *     F(s) = Int_0^Infinity f(t) E^(-s t) dt.
 *
 * Algorithm (each step exact; the integral is the only analysis).  The
 * factor classification is native C; Expand, Integrate and Simplify are
 * called as builtins:
 *   1. Lists thread; sums are transformed term by term (linearity); factors
 *      free of t come out.
 *   2. Exponential shift: E^(a t + b) g(t) -> E^b G(s - a), G = L[g].  This
 *      keeps the convergence region of the remaining integral at Re s > 0
 *      whatever the exponential rate a.
 *   3. A t-free remainder is c/s.
 *   4. Otherwise the defining integral, under the user's assumptions and
 *      s > 0 (the abscissa of convergence of the shifted g is 0 for every
 *      transform of exponential order 0 that Integrate closes).  If Integrate
 *      leaves the integral unevaluated, so is the transform.
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
    if (as->type == EXPR_SYMBOL && as->data.symbol.name == SYM_True)
        return lt_fn(SYM_LaplaceTransform, (Expr*[]){ expr_copy((Expr*)g), expr_copy((Expr*)t),
                                                      expr_copy((Expr*)s) }, 3);
    return lt_fn(SYM_LaplaceTransform, (Expr*[]){ expr_copy((Expr*)g), expr_copy((Expr*)t),
        expr_copy((Expr*)s), lt_fn2(SYM_Rule, expr_new_symbol(SYM_Assumptions),
                                    expr_copy((Expr*)as)) }, 4);
}

static Expr* lt_transform(const Expr* f, const Expr* t, const Expr* s, const Expr* as) {
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
    Expr *p0 = NULL, *p1 = NULL;
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
