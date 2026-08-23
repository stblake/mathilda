/*
 * reduce_form.c
 *
 * The DNF logical normal-form algebra for `Reduce` (see reduce_form.h and
 * REDUCE_PLAN.md).  This file owns the RAtom / RConj / RForm data structures,
 * their memory management, the And/Or/Not combinators, the constant-atom
 * simplifier, and emission back to a Mathematica-form Expr.  Atom
 * canonicalisation and input parsing live in reduce_atom.c.
 */
#include "reduce_form.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  RAtom                                                              *
 * ------------------------------------------------------------------ */

void ratom_free_members(RAtom* a) {
    if (!a) return;
    if (a->poly) { expr_free(a->poly); a->poly = NULL; }
    if (a->elem_dom) { expr_free(a->elem_dom); a->elem_dom = NULL; }
    if (a->display) { expr_free(a->display); a->display = NULL; }
}

RAtom ratom_dup(const RAtom* a) {
    RAtom r = *a;
    r.poly = a->poly ? expr_copy(a->poly) : NULL;
    r.elem_dom = a->elem_dom ? expr_copy(a->elem_dom) : NULL;
    r.display = a->display ? expr_copy(a->display) : NULL;
    return r;
}

bool ratom_eq(const RAtom* x, const RAtom* y) {
    if (x->rel != y->rel) return false;
    if (!expr_eq(x->poly, y->poly)) return false;
    if (x->rel == R_ELEM) {
        if (!x->elem_dom || !y->elem_dom) return x->elem_dom == y->elem_dom;
        return expr_eq(x->elem_dom, y->elem_dom);
    }
    return true;
}

RAtom ratom_negate(const RAtom* a, bool* ok) {
    RAtom r;
    r.elem_dom = NULL;
    r.display = NULL;
    r.nonconst_denom = a->nonconst_denom;   /* negation shares the denominator */
    r.main_var = a->main_var;
    r.deg_main = a->deg_main;
    r.is_linear = a->is_linear;
    switch (a->rel) {
        case R_EQ: r.rel = R_NE; r.poly = expr_copy(a->poly); break;
        case R_NE: r.rel = R_EQ; r.poly = expr_copy(a->poly); break;
        case R_LT: /* !(p<0) == (p>=0) == (-p<=0) */
            r.rel = R_LE;
            r.poly = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_copy(a->poly) }, 2));
            break;
        case R_LE: /* !(p<=0) == (p>0) == (-p<0) */
            r.rel = R_LT;
            r.poly = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_copy(a->poly) }, 2));
            break;
        default:   /* R_ELEM negation not supported yet */
            r.rel = a->rel;
            r.poly = a->poly ? expr_copy(a->poly) : NULL;
            if (ok) *ok = false;
            break;
    }
    return r;
}

/* ------------------------------------------------------------------ *
 *  RConj                                                              *
 * ------------------------------------------------------------------ */

RConj* rconj_new(void) {
    RConj* c = calloc(1, sizeof(RConj));
    return c;
}

void rconj_free(RConj* c) {
    if (!c) return;
    for (int i = 0; i < c->n; i++) ratom_free_members(&c->a[i]);
    free(c->a);
    free(c);
}

void rconj_push(RConj* c, RAtom a) {
    /* structural dedup: drop an atom already present */
    for (int i = 0; i < c->n; i++) {
        if (ratom_eq(&c->a[i], &a)) { ratom_free_members(&a); return; }
    }
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 4;
        c->a = realloc(c->a, (size_t)c->cap * sizeof(RAtom));
    }
    c->a[c->n++] = a;
}

RConj* rconj_dup(const RConj* c) {
    RConj* r = rconj_new();
    r->is_false = c->is_false;
    if (c->n) {
        r->cap = c->n;
        r->a = malloc((size_t)r->cap * sizeof(RAtom));
        for (int i = 0; i < c->n; i++) r->a[r->n++] = ratom_dup(&c->a[i]);
    }
    return r;
}

/* Merge a deep copy of every atom of `src` into `dst`. */
static void rconj_merge_into(RConj* dst, const RConj* src) {
    if (src->is_false) dst->is_false = true;
    for (int i = 0; i < src->n; i++) rconj_push(dst, ratom_dup(&src->a[i]));
}

/* ------------------------------------------------------------------ *
 *  RForm                                                              *
 * ------------------------------------------------------------------ */

RForm* rform_new(void) {
    return calloc(1, sizeof(RForm));
}

void rform_free(RForm* f) {
    if (!f) return;
    for (int i = 0; i < f->n; i++) rconj_free(f->c[i]);
    free(f->c);
    free(f);
}

RForm* rform_true(void) {
    RForm* f = rform_new();
    f->is_true = true;
    return f;
}

RForm* rform_false(void) {
    return rform_new();   /* n == 0, is_true == false */
}

void rform_add_conj(RForm* f, RConj* c) {
    if (f->n == f->cap) {
        f->cap = f->cap ? f->cap * 2 : 4;
        f->c = realloc(f->c, (size_t)f->cap * sizeof(RConj*));
    }
    f->c[f->n++] = c;
}

RForm* rform_single_atom(RAtom a) {
    RForm* f = rform_new();
    RConj* c = rconj_new();
    rconj_push(c, a);
    rform_add_conj(f, c);
    return f;
}

RForm* rform_or(RForm* a, RForm* b) {
    if (a->is_true) { rform_free(b); return a; }
    if (b->is_true) { rform_free(a); return b; }
    /* concatenate conjunction lists, transferring ownership of b's conjs */
    for (int i = 0; i < b->n; i++) rform_add_conj(a, b->c[i]);
    b->n = 0;
    rform_free(b);
    return a;
}

RForm* rform_and(RForm* a, RForm* b) {
    if (a->is_true) { rform_free(a); return b; }
    if (b->is_true) { rform_free(b); return a; }
    if (a->n == 0) { rform_free(a); rform_free(b); return rform_false(); }
    if (b->n == 0) { rform_free(a); rform_free(b); return rform_false(); }
    /* distributive product: (∨ Ai) ∧ (∨ Bj) = ∨_{i,j} (Ai ∧ Bj) */
    RForm* out = rform_new();
    for (int i = 0; i < a->n; i++) {
        for (int j = 0; j < b->n; j++) {
            RConj* nc = rconj_dup(a->c[i]);
            rconj_merge_into(nc, b->c[j]);
            rform_add_conj(out, nc);
        }
    }
    rform_free(a);
    rform_free(b);
    return out;
}

void rform_simplify(RForm* f, Expr** vars, int nv) {
    (void)vars; (void)nv;   /* Phase 0: only the constant-atom pass */
    if (!f || f->is_true) return;

    int w = 0;   /* write cursor for surviving conjunctions */
    for (int i = 0; i < f->n; i++) {
        RConj* c = f->c[i];
        if (c->is_false) { rconj_free(c); continue; }

        /* Decide constant atoms in place. */
        int cw = 0;
        for (int k = 0; k < c->n; k++) {
            int v = reduce_atom_decide_constant(&c->a[k]);
            if (v == 1) {            /* trivially true: drop the atom */
                ratom_free_members(&c->a[k]);
            } else if (v == 0) {     /* trivially false: kill the conjunction */
                c->is_false = true;
                /* free remaining atoms including this one */
                for (int m = k; m < c->n; m++) ratom_free_members(&c->a[m]);
                /* free atoms already kept */
                for (int m = 0; m < cw; m++) ratom_free_members(&c->a[m]);
                c->n = 0;
                break;
            } else {                 /* undecided: keep */
                c->a[cw++] = c->a[k];
            }
        }
        if (c->is_false) { rconj_free(c); continue; }
        c->n = cw;

        if (c->n == 0) {             /* empty conjunction == True */
            rconj_free(c);
            for (int j = 0; j < w; j++) rconj_free(f->c[j]);
            f->n = 0;
            f->is_true = true;
            return;
        }
        f->c[w++] = c;
    }
    f->n = w;

    /* Dedup identical conjunctions (cheap O(n^2), n small in practice). */
    for (int i = 0; i < f->n; i++) {
        for (int j = i + 1; j < f->n; j++) {
            RConj* ci = f->c[i];
            RConj* cj = f->c[j];
            if (ci->n != cj->n) continue;
            bool same = true;
            for (int k = 0; k < ci->n && same; k++) {
                bool found = false;
                for (int m = 0; m < cj->n; m++)
                    if (ratom_eq(&ci->a[k], &cj->a[m])) { found = true; break; }
                if (!found) same = false;
            }
            if (same) {
                rconj_free(cj);
                for (int m = j; m < f->n - 1; m++) f->c[m] = f->c[m + 1];
                f->n--;
                j--;
            }
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Emission back to Expr                                              *
 * ------------------------------------------------------------------ */

/* Build the Expr for one atom:  poly REL 0   (or Element[poly, dom]).  A solved
 * atom carrying a `display` override emits that verbatim. */
static Expr* ratom_to_expr(const RAtom* a) {
    if (a->display) return expr_copy(a->display);
    if (a->rel == R_ELEM) {
        return expr_new_function(
            expr_new_symbol(SYM_Element),
            (Expr*[]){ expr_copy(a->poly), expr_copy(a->elem_dom) }, 2);
    }
    const char* head;
    switch (a->rel) {
        case R_EQ: head = SYM_Equal;     break;
        case R_NE: head = SYM_Unequal;   break;
        case R_LT: head = SYM_Less;      break;
        case R_LE: head = SYM_LessEqual; break;
        default:   head = SYM_Equal;     break;   /* unreachable */
    }
    return expr_new_function(
        expr_new_symbol(head),
        (Expr*[]){ expr_copy(a->poly), expr_new_integer(0) }, 2);
}

/* Build the Expr for one conjunction. */
static Expr* rconj_to_expr(const RConj* c) {
    if (c->n == 1) return ratom_to_expr(&c->a[0]);
    Expr** args = malloc((size_t)c->n * sizeof(Expr*));
    for (int i = 0; i < c->n; i++) args[i] = ratom_to_expr(&c->a[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_And), args, (size_t)c->n);
    free(args);
    return out;
}

Expr* rform_to_expr(const RForm* f, Expr** vars, int nv) {
    (void)vars; (void)nv;
    if (f->is_true) return expr_new_symbol(SYM_True);
    if (f->n == 0)  return expr_new_symbol(SYM_False);

    Expr* built;
    if (f->n == 1) {
        built = rconj_to_expr(f->c[0]);
    } else {
        Expr** args = malloc((size_t)f->n * sizeof(Expr*));
        for (int i = 0; i < f->n; i++) args[i] = rconj_to_expr(f->c[i]);
        built = expr_new_function(expr_new_symbol(SYM_Or), args, (size_t)f->n);
        free(args);
    }
    /* One evaluate() pass flattens nested And/Or and tidies the relations. */
    return eval_and_free(built);
}
