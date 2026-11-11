/* algfield.c -- Together / Cancel over the number field of the constant
 * algebraic coefficients.  See algfield.h for the contract.
 *
 * The field engine (algfield_normalize), in four steps:
 *
 *   1. Collect the algebraic atoms.  Read through Plus / Times / Power (the
 *      structure Together combines through), every radical b^(p/q) with a
 *      constant algebraic base contributes the generator b^(1/q), every
 *      Complex the generator I, every constant Root object itself.  Bases are
 *      searched too, so Sqrt[1 + Sqrt[2]] contributes Sqrt[2] as well.
 *
 *   2. Build K = Q(generators) once.  A nested generator whose exact canonical
 *      form (flint_qqbar_canonical) is a radical expression over Q -- a
 *      denestable one, Sqrt[3/4 - 1/Sqrt[2]] = (-1 + Sqrt[2])/2 -- is not a
 *      field generator at all: its atoms join the list and it is located in K
 *      afterwards.  The rest go to one primitive-element construction
 *      (flint_qqbar_to_number_field_common).  Every atom is now an
 *      AlgebraicNumber[theta, {..}].
 *
 *   3. Substitute and combine.  The tree with AlgebraicNumber coefficients is
 *      evaluated (the evaluator's same-field arithmetic folds the constants)
 *      and handed to Together / Cancel, whose AlgebraicNumber paths reduce
 *      over K[vars].
 *
 *   4. Map back.  K has a Q-basis of radical monomials in the generators: for
 *      the generators in order (fewest nesting levels, then smallest), the
 *      basis B of the field generated so far grows by B g, B g^2, ... until a
 *      block is dependent (that power is the degree of g over the field so
 *      far).  An atom that lies in the field of the earlier ones adds nothing,
 *      which is what makes the answer denested.  Each AlgebraicNumber of the
 *      result is rewritten in that basis through the inverse of the d x d
 *      change-of-basis matrix.
 *
 * Exact throughout (qqbar / antic / fmpq), so no numeric check is needed. */

#include "algfield.h"

#include "sym_names.h"
#include "eval.h"
#include "common.h"
#include "core.h"
#include "arithmetic.h"
#include "flint_qqbar.h"
#include "poly.h"
#include "rat_internal.h"

#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_vec.h>
#endif

/* Above this field degree the basis and the antic arithmetic stop being
 * cheap; decline and leave the input to the generator model. */
#define AF_DEGREE_CAP 32
/* Distinct algebraic generators collected before declining. */
#define AF_MAX_GENS   12

/* ------------------------------------------------------------------ */
/*  Small structural helpers                                           */
/* ------------------------------------------------------------------ */

static bool af_is_ratnum(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
           (head_is(e, SYM_Rational) && e->data.function.arg_count == 2);
}

/* `e` is Rational[p, q] with machine-size p, q and q >= 2. */
static bool af_rat_exponent(const Expr* e, int64_t* p, int64_t* q) {
    if (!head_is(e, SYM_Rational) || e->data.function.arg_count != 2) return false;
    const Expr* n = e->data.function.args[0];
    const Expr* d = e->data.function.args[1];
    if (n->type != EXPR_INTEGER || d->type != EXPR_INTEGER) return false;
    if (d->data.integer < 2) return false;
    *p = n->data.integer;
    *q = d->data.integer;
    return true;
}

static bool af_is_int_exponent(const Expr* e) {
    return e->type == EXPR_INTEGER;
}

#ifdef USE_FLINT
/* The number-field path's helpers (af_constant_radical, af_depth) are compiled
 * only with FLINT: an unused static is an error under
 * -Werror=unused-function in the no-FLINT builds (Windows, USE_MPFR=0). */

/* A radical b^(p/q) whose base is a constant algebraic number. */
static bool af_constant_radical(const Expr* e, int64_t* p, int64_t* q) {
    if (!head_is(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    if (!af_rat_exponent(e->data.function.args[1], p, q)) return false;
    const Expr* b = e->data.function.args[0];
    if (af_is_ratnum(b)) return true;
    return b->type == EXPR_FUNCTION && flint_qqbar_is_constant_algebraic(b);
}

/* Nesting depth: the number of radical layers above the innermost rational
 * radicand (Sqrt[2] and I are 0, Sqrt[1 + Sqrt[2]] is 1). */
static int af_depth(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return 0;
    int64_t p, q;
    if (af_constant_radical(e, &p, &q)) {
        const Expr* b = e->data.function.args[0];
        return af_is_ratnum(b) ? 0 : 1 + af_depth(b);
    }
    int m = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        int d = af_depth(e->data.function.args[i]);
        if (d > m) m = d;
    }
    return m;
}
#endif /* USE_FLINT */

/* Does `e` hold a radical (a power with a non-integer rational exponent) or a
 * Root object? */
static bool af_has_radical(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Root)) return true;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2 &&
        head_is(e->data.function.args[1], SYM_Rational))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (af_has_radical(e->data.function.args[i])) return true;
    return false;
}

/* The base of a NESTED constant radical: a constant algebraic number that
 * itself carries a radical.  A Gaussian base (I^(1/3), (1 + I)^(1/2)) is a
 * root of unity times a radical of a rational and keeps the paths it had. */
static bool af_nested_base(const Expr* b) {
    return b->type == EXPR_FUNCTION && !af_is_ratnum(b) && af_has_radical(b) &&
           flint_qqbar_is_constant_algebraic(b);
}

static bool af_gate_walk(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (af_gate_walk(e->data.function.args[i])) return true;
        return false;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        int64_t p, q;
        if (af_rat_exponent(e->data.function.args[1], &p, &q)) return af_nested_base(b);
        if (af_is_int_exponent(e->data.function.args[1])) return af_gate_walk(b);
    }
    return false;
}

bool algfield_has_nested_constant_radical(const Expr* e) {
    return e && af_gate_walk(e);
}

/* ------------------------------------------------------------------ */
/*  Default mode: nested radicals as opaque generators                 */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr**   base;   /* owned radicand b                         */
    int64_t* m;      /* lcm of the exponent denominators         */
    Expr**   sym;    /* owned fresh generator symbol g = b^(1/m) */
    size_t   n, cap;
} AfKern;

static void af_kern_free(AfKern* k) {
    for (size_t i = 0; i < k->n; i++) {
        expr_free(k->base[i]);
        if (k->sym[i]) expr_free(k->sym[i]);
    }
    free(k->base); free(k->m); free(k->sym);
}

static int64_t af_gcd64(int64_t a, int64_t b) {
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
}

static long af_kern_find(const AfKern* k, const Expr* b) {
    for (size_t i = 0; i < k->n; i++)
        if (expr_eq(k->base[i], b)) return (long)i;
    return -1;
}

/* Gather the outermost nested constant radicals (the gate's atoms). */
static bool af_kern_gather(const Expr* e, AfKern* k) {
    if (e->type != EXPR_FUNCTION) return true;
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!af_kern_gather(e->data.function.args[i], k)) return false;
        return true;
    }
    if (!head_is(e, SYM_Power) || e->data.function.arg_count != 2) return true;
    const Expr* b = e->data.function.args[0];
    int64_t p, q;
    if (af_rat_exponent(e->data.function.args[1], &p, &q)) {
        if (!af_nested_base(b)) return true;
        long at = af_kern_find(k, b);
        if (at >= 0) {
            int64_t g = af_gcd64(k->m[at], q);
            if (k->m[at] / g > (INT64_C(1) << 20) / q) return false;   /* lcm too large */
            k->m[at] = k->m[at] / g * q;
            return true;
        }
        if (k->n == k->cap) {
            size_t nc = k->cap ? 2 * k->cap : 4;
            k->base = realloc(k->base, nc * sizeof(Expr*));
            k->m    = realloc(k->m,    nc * sizeof(int64_t));
            k->sym  = realloc(k->sym,  nc * sizeof(Expr*));
            k->cap = nc;
        }
        k->base[k->n] = expr_copy((Expr*)b);
        k->m[k->n] = q;
        k->sym[k->n] = NULL;
        k->n++;
        return k->n <= AF_MAX_GENS;
    }
    if (af_is_int_exponent(e->data.function.args[1])) return af_kern_gather(b, k);
    return true;
}

/* b^(p/q) -> g^(p m/q) along the same structure. */
static Expr* af_kern_forward(const Expr* e, const AfKern* k) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    bool ring = head_is(e, SYM_Plus) || head_is(e, SYM_Times);
    if (!ring && head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        int64_t p, q;
        if (af_rat_exponent(e->data.function.args[1], &p, &q)) {
            long at = af_kern_find(k, e->data.function.args[0]);
            if (at < 0) return expr_copy((Expr*)e);
            return expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_copy(k->sym[at]),
                           expr_new_integer(p * (k->m[at] / q)) }, 2);
        }
        if (!af_is_int_exponent(e->data.function.args[1])) return expr_copy((Expr*)e);
        ring = true;
    }
    if (!ring) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = af_kern_forward(e->data.function.args[i], k);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

/* g -> b^(1/m) everywhere (the generators are fresh, so a plain walk). */
static Expr* af_kern_back(const Expr* e, const AfKern* k) {
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < k->n; i++)
            if (expr_eq((Expr*)e, k->sym[i]))
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(k->base[i]), make_rational(1, k->m[i]) }, 2);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = af_kern_back(e->data.function.args[i], k);
    Expr* r = expr_new_function(af_kern_back(e->data.function.head, k), a, n);
    free(a);
    return r;
}

/* compute(head[x]) with x consumed; NULL when the builtin declines. */
static Expr* af_call1(Expr* (*compute)(Expr*), const char* head, Expr* x) {
    Expr* call = expr_new_function(expr_new_symbol(head), (Expr*[]){ x }, 1);
    Expr* r = compute(call);
    expr_free(call);
    return r;
}

static bool af_is_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Each nested radical's base in its Together form, along the gate's
 * structure.  A radicand that combines to a radical monomial
 * (-1/(9 2^(2/3)) + 2/9 2^(1/3) = 1/(3 2^(2/3))) is then not nested at all. */
static Expr* af_normalize_bases(const Expr* e, bool* changed) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    bool ring = head_is(e, SYM_Plus) || head_is(e, SYM_Times);
    if (!ring && head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        int64_t p, q;
        if (af_rat_exponent(e->data.function.args[1], &p, &q)) {
            if (!af_nested_base(b)) return expr_copy((Expr*)e);
            Expr* nb = eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
                (Expr*[]){ expr_copy((Expr*)b) }, 1));
            if (expr_eq(nb, b)) { expr_free(nb); return expr_copy((Expr*)e); }
            *changed = true;
            return expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ nb, expr_copy(e->data.function.args[1]) }, 2);
        }
        ring = af_is_int_exponent(e->data.function.args[1]);
    }
    if (!ring) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = af_normalize_bases(e->data.function.args[i], changed);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

static Expr* af_kernel_model_core(const Expr* arg, const char* head,
                                  Expr* (*compute)(Expr*));

Expr* algfield_kernel_model(const Expr* call, const char* head,
                            Expr* (*compute)(Expr*)) {
    if (!call || call->type != EXPR_FUNCTION || call->data.function.arg_count != 1)
        return NULL;
    const Expr* arg = call->data.function.args[0];
    bool changed = false;
    Expr* nb = af_normalize_bases(arg, &changed);
    if (!changed) {
        expr_free(nb);
        return af_kernel_model_core(arg, head, compute);
    }
    nb = eval_and_free(nb);
    Expr* out;
    if (algfield_has_nested_constant_radical(nb)) {
        out = af_kernel_model_core(nb, head, compute);
        if (!out) { out = nb; nb = NULL; }
    } else {
        /* nothing nested is left: the ordinary cascade handles it */
        out = af_call1(compute, head, expr_copy(nb));
        if (!out) { out = nb; nb = NULL; }
    }
    if (nb) expr_free(nb);
    return out;
}

static Expr* af_kernel_model_core(const Expr* arg, const char* head,
                                  Expr* (*compute)(Expr*)) {
    AfKern k = {0};
    if (!af_kern_gather(arg, &k) || k.n == 0) { af_kern_free(&k); return NULL; }

    /* fresh generator names, distinct from every symbol of the input and from
     * one another (each is checked against the input with the previous ones
     * already substituted) */
    Expr* sub = expr_copy((Expr*)arg);
    for (size_t i = 0; i < k.n; i++) {
        char* nm = poly_make_fresh_gen(sub);
        k.sym[i] = expr_new_symbol(nm);
        free(nm);
        AfKern one = k; one.n = i + 1;
        Expr* s2 = af_kern_forward(sub, &one);
        expr_free(sub);
        sub = s2;
    }
    sub = eval_and_free(sub);

    Expr* r = af_call1(compute, head, sub);
    if (!r) { af_kern_free(&k); return NULL; }

    /* g^m = b: reduce numerator and denominator below degree m in each g.
     * Value-preserving, and what keeps the answer compact; for Cancel only a
     * single fraction is reduced (a sum stays a sum). */
    if (head == SYM_Together || !head_is(r, SYM_Plus)) {
        Expr* num; Expr* den;
        extract_num_den(r, &num, &den);
        for (size_t i = 0; i < k.n; i++) {
            Expr* rel = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(k.sym[i]), expr_new_integer(k.m[i]) }, 2),
                           expr_new_function(expr_new_symbol(SYM_Times),
                               (Expr*[]){ expr_new_integer(-1), expr_copy(k.base[i]) }, 2) }, 2);
            rel = eval_and_free(rel);
            num = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialRemainder),
                (Expr*[]){ num, expr_copy(rel), expr_copy(k.sym[i]) }, 3));
            den = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialRemainder),
                (Expr*[]){ den, rel, expr_copy(k.sym[i]) }, 3));
        }
        if (af_is_zero(den)) {
            /* cannot happen for a well-defined input (the denominator's value
             * is unchanged by the reduction); keep the unreduced answer */
            expr_free(num); expr_free(den);
        } else {
            Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2));
            Expr* r2 = af_call1(compute, head, q);
            if (r2) { expr_free(r); r = r2; }
        }
    }
    Expr* out = eval_and_free(af_kern_back(r, &k));
    expr_free(r);
    af_kern_free(&k);
    return out;
}

/* ------------------------------------------------------------------ */
/*  The field engine                                                   */
/* ------------------------------------------------------------------ */

#ifdef USE_FLINT

typedef struct {
    Expr** g;        /* owned, evaluated generator expressions          */
    int*   depth;
    bool*  derived;  /* denestable: located in K after it is built      */
    Expr** an;       /* owned AlgebraicNumber (or rational) of g[i]      */
    size_t n, cap;
} AfGens;

static void af_gens_free(AfGens* G) {
    for (size_t i = 0; i < G->n; i++) {
        expr_free(G->g[i]);
        if (G->an[i]) expr_free(G->an[i]);
    }
    free(G->g); free(G->depth); free(G->derived); free(G->an);
}

/* Takes ownership of `gen`.  False once the generator budget is exceeded. */
static bool af_add_gen(AfGens* G, Expr* gen) {
    for (size_t i = 0; i < G->n; i++)
        if (expr_eq(G->g[i], gen)) { expr_free(gen); return true; }
    if (G->n == AF_MAX_GENS) { expr_free(gen); return false; }
    if (G->n == G->cap) {
        size_t nc = G->cap ? 2 * G->cap : 4;
        G->g       = realloc(G->g,       nc * sizeof(Expr*));
        G->depth   = realloc(G->depth,   nc * sizeof(int));
        G->derived = realloc(G->derived, nc * sizeof(bool));
        G->an      = realloc(G->an,      nc * sizeof(Expr*));
        G->cap = nc;
    }
    G->g[G->n] = gen;
    G->depth[G->n] = af_depth(gen);
    G->derived[G->n] = false;
    G->an[G->n] = NULL;
    G->n++;
    return true;
}

static Expr* af_imag_unit(void) {
    return make_complex(expr_new_integer(0), expr_new_integer(1));
}

static Expr* af_root_gen(const Expr* b, int64_t q) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)b), make_rational(1, q) }, 2));
}

static bool af_exact_num(const Expr* e) {
    return af_is_ratnum(e);
}

/* Step 1.  False declines the engine (an inexact or AlgebraicNumber
 * coefficient, an operator it does not read, too many generators). */
static bool af_collect(const Expr* e, AfGens* G) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_SYMBOL: return true;
        case EXPR_FUNCTION: break;
        default: return false;
    }
    if (head_is(e, SYM_Rational)) return true;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        if (!af_exact_num(e->data.function.args[0]) || !af_exact_num(e->data.function.args[1]))
            return false;
        return af_add_gen(G, af_imag_unit());
    }
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!af_collect(e->data.function.args[i], G)) return false;
        return true;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        int64_t p, q;
        if (af_is_int_exponent(e->data.function.args[1])) return af_collect(b, G);
        if (af_constant_radical(e, &p, &q)) {
            if (!af_add_gen(G, af_root_gen(b, q))) return false;
            return af_collect(b, G);
        }
        return false;
    }
    if (head_is(e, SYM_Root) && flint_qqbar_is_constant_algebraic(e))
        return af_add_gen(G, expr_copy((Expr*)e));
    return false;
}

static bool af_has_head(const Expr* e, const char* h) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (af_has_head(e->data.function.args[i], h)) return true;
    return false;
}

/* generator order: fewest nesting levels, then smallest, then canonical */
static int af_gen_cmp(const AfGens* G, size_t i, size_t j) {
    if (G->depth[i] != G->depth[j]) return G->depth[i] < G->depth[j] ? -1 : 1;
    int64_t li = leaf_count_internal(G->g[i], true), lj = leaf_count_internal(G->g[j], true);
    if (li != lj) return li < lj ? -1 : 1;
    return expr_compare(G->g[i], G->g[j]);
}

static void af_gens_sort(AfGens* G) {
    for (size_t i = 1; i < G->n; i++)          /* n <= AF_MAX_GENS: insertion */
        for (size_t j = i; j > 0 && af_gen_cmp(G, j - 1, j) > 0; j--) {
            Expr* t = G->g[j]; G->g[j] = G->g[j - 1]; G->g[j - 1] = t;
            int d = G->depth[j]; G->depth[j] = G->depth[j - 1]; G->depth[j - 1] = d;
            bool v = G->derived[j]; G->derived[j] = G->derived[j - 1]; G->derived[j - 1] = v;
            Expr* a = G->an[j]; G->an[j] = G->an[j - 1]; G->an[j - 1] = a;
        }
}

static bool af_is_algnum(const Expr* e) {
    return head_is(e, SYM_AlgebraicNumber) && e->data.function.arg_count == 2 &&
           head_is(e->data.function.args[1], SYM_List);
}

/* Step 2.  On success *theta is the shared generator (NULL when K = Q) and
 * *d the field degree. */
static bool af_build_field(AfGens* G, Expr** theta, slong* d) {
    *theta = NULL; *d = 1;
    /* denestable nested generators: their radical atoms replace them as
     * field generators */
    size_t n0 = G->n;
    for (size_t i = 0; i < n0; i++) {
        if (G->depth[i] == 0) continue;
        Expr* c = flint_qqbar_canonical(G->g[i], QQBAR_METHOD_AUTOMATIC);
        if (!c) continue;
        if (!af_has_head(c, SYM_Root) && af_depth(c) == 0) {
            AfGens tmp = {0};
            bool ok = af_collect(c, &tmp);
            if (ok) {
                G->derived[i] = true;
                for (size_t j = 0; j < tmp.n && ok; j++) {
                    ok = af_add_gen(G, tmp.g[j]);
                    tmp.g[j] = NULL;
                }
            }
            for (size_t j = 0; j < tmp.n; j++) if (tmp.g[j]) expr_free(tmp.g[j]);
            tmp.n = 0;
            af_gens_free(&tmp);
            if (!ok) { expr_free(c); return false; }
        }
        expr_free(c);
    }
    af_gens_sort(G);

    const Expr** F = malloc(G->n * sizeof(Expr*));
    size_t nf = 0;
    for (size_t i = 0; i < G->n; i++) if (!G->derived[i]) F[nf++] = G->g[i];
    bool ok = true;
    if (nf > 0) {
        Expr* L = flint_qqbar_to_number_field_common(F, nf, 0);
        if (!L || !head_is(L, SYM_List) || L->data.function.arg_count != nf) ok = false;
        size_t k = 0;
        for (size_t i = 0; i < G->n && ok; i++) {
            if (G->derived[i]) continue;
            G->an[i] = expr_copy(L->data.function.args[k++]);
            if (af_is_algnum(G->an[i])) {
                const Expr* t = G->an[i]->data.function.args[0];
                if (!*theta) *theta = expr_copy((Expr*)t);
                else if (!expr_eq(*theta, t)) ok = false;
            } else if (!af_is_ratnum(G->an[i])) ok = false;
        }
        if (L) expr_free(L);
    }
    free(F);
    for (size_t i = 0; i < G->n && ok; i++) {
        if (!G->derived[i]) continue;
        G->an[i] = *theta ? flint_qqbar_to_number_field(G->g[i], *theta)
                          : flint_qqbar_canonical(G->g[i], QQBAR_METHOD_AUTOMATIC);
        if (!G->an[i]) { ok = false; break; }
        if (af_is_algnum(G->an[i])) {
            if (!*theta || !expr_eq(*theta, G->an[i]->data.function.args[0])) ok = false;
        } else if (!af_is_ratnum(G->an[i])) ok = false;
    }
    if (ok && *theta) {
        /* the degree: the length of any power-basis coefficient list */
        for (size_t i = 0; i < G->n; i++)
            if (af_is_algnum(G->an[i])) {
                slong m = (slong)G->an[i]->data.function.args[1]->data.function.arg_count;
                if (m > *d) *d = m;
            }
        Expr* mp = flint_qqbar_gen_minpoly_coeffs(*theta);
        if (mp && head_is(mp, SYM_List) && mp->data.function.arg_count >= 2)
            *d = (slong)mp->data.function.arg_count - 1;
        if (mp) expr_free(mp);
    }
    if (ok && *d > AF_DEGREE_CAP) ok = false;
    if (!ok && *theta) { expr_free(*theta); *theta = NULL; }
    return ok;
}

static long af_gen_index(const AfGens* G, const Expr* gen) {
    for (size_t i = 0; i < G->n; i++) if (expr_eq(G->g[i], gen)) return (long)i;
    return -1;
}

/* Step 3's substitution, along af_collect's structure. */
static Expr* af_subst(const Expr* e, const AfGens* G, bool* ok) {
    if (e->type != EXPR_FUNCTION || head_is(e, SYM_Rational)) return expr_copy((Expr*)e);
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        Expr* ig = af_imag_unit();
        long at = af_gen_index(G, ig);
        expr_free(ig);
        if (at < 0) { *ok = false; return expr_copy((Expr*)e); }
        return expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(e->data.function.args[0]),
                       expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_copy(e->data.function.args[1]),
                                      expr_copy(G->an[at]) }, 2) }, 2);
    }
    if (head_is(e, SYM_Root)) {
        long at = af_gen_index(G, e);
        if (at < 0) { *ok = false; return expr_copy((Expr*)e); }
        return expr_copy(G->an[at]);
    }
    int64_t p, q;
    if (af_constant_radical(e, &p, &q)) {
        Expr* gen = af_root_gen(e->data.function.args[0], q);
        long at = af_gen_index(G, gen);
        expr_free(gen);
        if (at < 0) { *ok = false; return expr_copy((Expr*)e); }
        return expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(G->an[at]), expr_new_integer(p) }, 2);
    }
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = af_subst(e->data.function.args[i], G, ok);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

static bool af_fmpq_from_expr(const Expr* e, fmpq_t out) {
    if (e->type == EXPR_INTEGER) { fmpq_set_si(out, e->data.integer, 1); return true; }
    if (e->type == EXPR_BIGINT) {
        fmpz_set_mpz(fmpq_numref(out), e->data.bigint);
        fmpz_one(fmpq_denref(out));
        return true;
    }
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) {
        fmpq_t n, d; fmpq_init(n); fmpq_init(d);
        bool ok = af_fmpq_from_expr(e->data.function.args[0], n) &&
                  af_fmpq_from_expr(e->data.function.args[1], d) && !fmpq_is_zero(d);
        if (ok) fmpq_div(out, n, d);
        fmpq_clear(n); fmpq_clear(d);
        return ok;
    }
    return false;
}

static Expr* af_expr_from_fmpq(const fmpq_t x) {
    mpz_t n, d; mpz_init(n); mpz_init(d);
    fmpz_get_mpz(n, fmpq_numref(x));
    fmpz_get_mpz(d, fmpq_denref(x));
    Expr* r = make_rational_mpz(n, d);
    mpz_clear(n); mpz_clear(d);
    return r;
}

/* The coordinates of x (an AlgebraicNumber over theta, or a rational) in the
 * power basis of theta: row `row` of `M`. */
static bool af_vec(const Expr* x, const Expr* theta, slong d, fmpq_mat_t M, slong row) {
    for (slong j = 0; j < d; j++) fmpq_zero(fmpq_mat_entry(M, row, j));
    if (af_is_ratnum(x)) return af_fmpq_from_expr(x, fmpq_mat_entry(M, row, 0));
    if (!af_is_algnum(x) || !theta || !expr_eq(x->data.function.args[0], theta)) return false;
    const Expr* L = x->data.function.args[1];
    if ((slong)L->data.function.arg_count > d) return false;
    for (size_t j = 0; j < L->data.function.arg_count; j++)
        if (!af_fmpq_from_expr(L->data.function.args[j], fmpq_mat_entry(M, row, (slong)j)))
            return false;
    return true;
}

typedef struct {
    slong      d;
    Expr**     mono;   /* owned basis monomials (radical expressions) */
    fmpq_mat_t inv;    /* coordinates -> basis coefficients           */
    const Expr* theta;
} AfBasis;

static Expr* af_times(Expr* a, Expr* b) {           /* consumes both */
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ a, b }, 2));
}

/* Step 4's basis.  `bv` (d x d) collects the power-basis vectors of the
 * accepted monomials; `ban` their AlgebraicNumbers. */
static bool af_basis_build(const AfGens* G, const Expr* theta, slong d, AfBasis* B) {
    B->d = d; B->theta = theta;
    B->mono = calloc((size_t)d, sizeof(Expr*));
    Expr** ban = calloc((size_t)d, sizeof(Expr*));
    fmpq_mat_t bv; fmpq_mat_init(bv, d, d);
    slong nb = 1;
    B->mono[0] = expr_new_integer(1);
    ban[0] = expr_new_integer(1);
    bool ok = af_vec(ban[0], theta, d, bv, 0);
    for (size_t gi = 0; gi < G->n && ok && nb < d; gi++) {
        slong nb0 = nb;
        Expr* gpow = expr_new_integer(1);        /* g^j as an AlgebraicNumber */
        for (int64_t j = 1; j < 64 && nb + nb0 <= d && ok; j++) {
            gpow = af_times(gpow, expr_copy(G->an[gi]));
            fmpq_mat_t t; fmpq_mat_init(t, nb + nb0, d);
            for (slong r = 0; r < nb; r++)
                for (slong c = 0; c < d; c++)
                    fmpq_set(fmpq_mat_entry(t, r, c), fmpq_mat_entry(bv, r, c));
            Expr** cand = malloc((size_t)nb0 * sizeof(Expr*));
            for (slong k = 0; k < nb0; k++) {
                cand[k] = af_times(expr_copy(ban[k]), expr_copy(gpow));
                if (ok && !af_vec(cand[k], theta, d, t, nb + k)) ok = false;
            }
            bool grow = false;
            if (ok) {
                fmpq_mat_t rr; fmpq_mat_init(rr, nb + nb0, d);
                grow = fmpq_mat_rref(rr, t) == nb + nb0;        /* returns the rank */
                fmpq_mat_clear(rr);
            }
            if (grow) {
                Expr* gj = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(G->g[gi]), expr_new_integer(j) }, 2));
                for (slong k = 0; k < nb0; k++) {
                    for (slong c = 0; c < d; c++)
                        fmpq_set(fmpq_mat_entry(bv, nb + k, c), fmpq_mat_entry(t, nb + k, c));
                    ban[nb + k] = cand[k];
                    B->mono[nb + k] = af_times(expr_copy(B->mono[k]), expr_copy(gj));
                }
                expr_free(gj);
                nb += nb0;
            } else {
                for (slong k = 0; k < nb0; k++) expr_free(cand[k]);
            }
            free(cand);
            fmpq_mat_clear(t);
            if (!grow) break;
        }
        expr_free(gpow);
    }
    if (ok && nb != d) ok = false;
    fmpq_mat_init(B->inv, d, d);
    if (ok) ok = fmpq_mat_inv(B->inv, bv) != 0;
    for (slong k = 0; k < d; k++) if (ban[k]) expr_free(ban[k]);
    free(ban);
    fmpq_mat_clear(bv);
    return ok;
}

static void af_basis_free(AfBasis* B) {
    if (B->mono) for (slong k = 0; k < B->d; k++) if (B->mono[k]) expr_free(B->mono[k]);
    free(B->mono);
    fmpq_mat_clear(B->inv);
}

/* AlgebraicNumber -> sum of basis monomials; anything else rebuilt. */
static Expr* af_unmap(const Expr* e, const AfBasis* B, bool* ok) {
    if (af_is_algnum(e)) {
        slong d = B->d;
        fmpq_mat_t v; fmpq_mat_init(v, 1, d);
        Expr* out = NULL;
        if (af_vec(e, B->theta, d, v, 0)) {
            fmpq_mat_t y; fmpq_mat_init(y, 1, d);
            fmpq_mat_mul(y, v, B->inv);
            Expr** terms = malloc((size_t)d * sizeof(Expr*));
            size_t nt = 0;
            for (slong k = 0; k < d; k++) {
                if (fmpq_is_zero(fmpq_mat_entry(y, 0, k))) continue;
                terms[nt++] = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ af_expr_from_fmpq(fmpq_mat_entry(y, 0, k)),
                               expr_copy(B->mono[k]) }, 2);
            }
            out = nt == 0 ? expr_new_integer(0)
                          : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
            free(terms);
            fmpq_mat_clear(y);
        } else {
            *ok = false;
        }
        fmpq_mat_clear(v);
        return out ? out : expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = af_unmap(e->data.function.args[i], B, ok);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

/* ---- presentation: the answer over K in the plain normal form ----------- */

/* One term of an expanded polynomial over K: its coefficient in the radical
 * basis (y, d entries) and its monomial in the variables. */
typedef struct {
    fmpq*  y;
    Expr*  mono;
    slong  deg;
} AfTerm;

static void af_terms_free(AfTerm* t, size_t n, slong d) {
    for (size_t i = 0; i < n; i++) {
        _fmpq_vec_clear(t[i].y, d);
        expr_free(t[i].mono);
    }
    free(t);
}

/* Split an expanded polynomial with rational / AlgebraicNumber coefficients
 * into terms.  *lead (optional) receives the coefficient of the leading term
 * (highest total degree, then the largest monomial) as an Expr over K. */
static AfTerm* af_split_poly(const Expr* P, const AfBasis* B, size_t* n_out,
                             Expr** lead, bool* ok) {
    bool sum = head_is(P, SYM_Plus);
    size_t nt = sum ? P->data.function.arg_count : 1;
    AfTerm* T = calloc(nt, sizeof(AfTerm));
    slong d = B->d;
    fmpq_mat_t v, yv; fmpq_mat_init(v, 1, d); fmpq_mat_init(yv, 1, d);
    long best = -1;
    Expr* best_coef = NULL;
    size_t made = 0;
    for (size_t i = 0; i < nt && *ok; i++) {
        const Expr* t = sum ? P->data.function.args[i] : P;
        bool prod = head_is(t, SYM_Times);
        size_t nf = prod ? t->data.function.arg_count : 1;
        fmpq_t r; fmpq_init(r); fmpq_one(r);
        const Expr* an = NULL;
        Expr** mf = malloc((nf ? nf : 1) * sizeof(Expr*));
        size_t nm = 0;
        slong deg = 0;
        for (size_t j = 0; j < nf && *ok; j++) {
            const Expr* f = prod ? t->data.function.args[j] : t;
            if (af_is_ratnum(f)) {
                fmpq_t q; fmpq_init(q);
                if (af_fmpq_from_expr(f, q)) fmpq_mul(r, r, q); else *ok = false;
                fmpq_clear(q);
            } else if (af_is_algnum(f)) {
                if (an) *ok = false;               /* the evaluator folds these */
                an = f;
            } else {
                mf[nm++] = expr_copy((Expr*)f);
                if (head_is(f, SYM_Power) && f->data.function.arg_count == 2 &&
                    f->data.function.args[1]->type == EXPR_INTEGER)
                    deg += (slong)f->data.function.args[1]->data.integer;
                else
                    deg += 1;
            }
        }
        if (*ok) {
            if (an) *ok = af_vec(an, B->theta, d, v, 0);
            else {
                for (slong k = 0; k < d; k++) fmpq_zero(fmpq_mat_entry(v, 0, k));
                fmpq_one(fmpq_mat_entry(v, 0, 0));
            }
        }
        Expr* mono = nm == 0 ? expr_new_integer(1)
                   : nm == 1 ? mf[0]
                   : expr_new_function(expr_new_symbol(SYM_Times), mf, nm);
        free(mf);
        T[made].y = _fmpq_vec_init(d);
        T[made].mono = mono;
        T[made].deg = deg;
        made++;
        if (*ok) {
            for (slong k = 0; k < d; k++)
                fmpq_mul(fmpq_mat_entry(v, 0, k), fmpq_mat_entry(v, 0, k), r);
            fmpq_mat_mul(yv, v, B->inv);
            for (slong k = 0; k < d; k++) fmpq_set(T[made - 1].y + k, fmpq_mat_entry(yv, 0, k));
            if (lead && (best < 0 || deg > T[best].deg ||
                         (deg == T[best].deg && expr_compare(mono, T[best].mono) > 0))) {
                best = (long)(made - 1);
                if (best_coef) expr_free(best_coef);
                Expr* rr = af_expr_from_fmpq(r);
                best_coef = an ? af_times(rr, expr_copy((Expr*)an)) : rr;
            }
        }
        fmpq_clear(r);
    }
    fmpq_mat_clear(v); fmpq_mat_clear(yv);
    *n_out = made;
    if (lead) *lead = best_coef;
    else if (best_coef) expr_free(best_coef);
    return T;
}

static Expr* af_expand(Expr* e) {                    /* consumes e */
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){ e }, 1));
}

/* sum_t sum_k (s y_tk) mono_k mono_t, evaluated */
static Expr* af_terms_expr(const AfTerm* T, size_t n, const AfBasis* B, const fmpq_t s) {
    size_t m = 0;
    Expr** a = malloc((n * (size_t)B->d + 1) * sizeof(Expr*));
    fmpq_t c; fmpq_init(c);
    for (size_t i = 0; i < n; i++)
        for (slong k = 0; k < B->d; k++) {
            if (fmpq_is_zero(T[i].y + k)) continue;
            fmpq_mul(c, T[i].y + k, s);
            a[m++] = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ af_expr_from_fmpq(c), expr_copy(B->mono[k]), expr_copy(T[i].mono) }, 3);
        }
    fmpq_clear(c);
    Expr* r = m == 0 ? expr_new_integer(0) : expr_new_function(expr_new_symbol(SYM_Plus), a, m);
    free(a);
    return eval_and_free(r);
}

/* Joint content of the two term lists: s = lcm(denominators)/gcd(numerators),
 * so s * every coefficient is an integer and together they have gcd 1. */
static void af_joint_scale(const AfTerm* TN, size_t nn, const AfTerm* TD, size_t nd,
                           slong d, fmpq_t s) {
    fmpz_t L, G, t; fmpz_init(L); fmpz_init(G); fmpz_init(t);
    fmpz_one(L);
    for (int pass = 0; pass < 2; pass++) {
        const AfTerm* T = pass ? TD : TN;
        size_t n = pass ? nd : nn;
        for (size_t i = 0; i < n; i++)
            for (slong k = 0; k < d; k++)
                if (!fmpq_is_zero(T[i].y + k)) fmpz_lcm(L, L, fmpq_denref(T[i].y + k));
    }
    for (int pass = 0; pass < 2; pass++) {
        const AfTerm* T = pass ? TD : TN;
        size_t n = pass ? nd : nn;
        for (size_t i = 0; i < n; i++)
            for (slong k = 0; k < d; k++) {
                if (fmpq_is_zero(T[i].y + k)) continue;
                fmpz_divexact(t, L, fmpq_denref(T[i].y + k));
                fmpz_mul(t, t, fmpq_numref(T[i].y + k));
                fmpz_gcd(G, G, t);
            }
    }
    if (fmpz_is_zero(G)) fmpz_one(G);
    fmpq_set_fmpz_frac(s, L, G);
    fmpz_clear(L); fmpz_clear(G); fmpz_clear(t);
}

/* r (over K, AlgebraicNumber coefficients) -> N/D: the denominator made monic
 * over K, every coefficient written in the radical basis, and the rational
 * content of N and D cleared jointly (integer coefficients with gcd 1, the
 * leading denominator coefficient positive). */
static Expr* af_present(const Expr* r, const AfBasis* B, bool* ok) {
    Expr* num; Expr* den;
    extract_num_den((Expr*)r, &num, &den);
    num = af_expand(num);
    den = af_expand(den);
    size_t nd = 0, nn = 0;
    Expr* lc = NULL;
    AfTerm* TD = af_split_poly(den, B, &nd, &lc, ok);
    af_terms_free(TD, nd, B->d);
    Expr* out = NULL;
    if (*ok && lc) {
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ lc, expr_new_integer(-1) }, 2));
        lc = NULL;
        num = af_expand(af_times(num, expr_copy(inv)));
        den = af_expand(af_times(den, inv));
        AfTerm* TN = af_split_poly(num, B, &nn, NULL, ok);
        TD = af_split_poly(den, B, &nd, NULL, ok);
        if (*ok) {
            fmpq_t s; fmpq_init(s);
            af_joint_scale(TN, nn, TD, nd, B->d, s);
            Expr* N = af_terms_expr(TN, nn, B, s);
            Expr* D = af_terms_expr(TD, nd, B, s);
            fmpq_clear(s);
            if (af_is_zero(D)) { expr_free(N); expr_free(D); *ok = false; }
            else out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ N, expr_new_function(expr_new_symbol(SYM_Power),
                                          (Expr*[]){ D, expr_new_integer(-1) }, 2) }, 2));
        }
        af_terms_free(TN, nn, B->d);
        af_terms_free(TD, nd, B->d);
    } else {
        *ok = false;
    }
    if (lc) expr_free(lc);
    expr_free(num); expr_free(den);
    return out;
}

/* A Plus with a denominator in some term: Cancel keeps such a sum a sum. */
static bool af_sum_of_fractions(const Expr* r) {
    if (!head_is(r, SYM_Plus)) return false;
    for (size_t i = 0; i < r->data.function.arg_count; i++) {
        Expr* n; Expr* d;
        extract_num_den(r->data.function.args[i], &n, &d);
        bool frac = !af_is_ratnum(d) && !af_is_algnum(d);
        expr_free(n); expr_free(d);
        if (frac) return true;
    }
    return false;
}

Expr* algfield_normalize(const Expr* e, const char* head) {
    if (!e || (head != SYM_Together && head != SYM_Cancel)) return NULL;
    AfGens G = {0};
    if (!af_collect(e, &G) || G.n == 0) { af_gens_free(&G); return NULL; }

    Expr* theta = NULL;
    slong d = 1;
    if (!af_build_field(&G, &theta, &d)) { af_gens_free(&G); return NULL; }

    Expr* out = NULL;
    bool ok = true;
    Expr* ek = af_subst(e, &G, &ok);
    if (ok) {
        ek = eval_and_free(ek);
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(head), (Expr*[]){ ek }, 1));
        AfBasis B = {0};
        if (af_basis_build(&G, theta, d, &B)) {
            if (head == SYM_Cancel && af_sum_of_fractions(r)) {
                Expr* m = af_unmap(r, &B, &ok);       /* Cancel keeps the sum */
                if (ok) out = eval_and_free(m);
                else expr_free(m);
            } else {
                out = af_present(r, &B, &ok);
                if (!ok && out) { expr_free(out); out = NULL; }
            }
        }
        af_basis_free(&B);
        expr_free(r);
    } else {
        expr_free(ek);
    }
    if (theta) expr_free(theta);
    af_gens_free(&G);
    return out;
}

#else  /* !USE_FLINT */

Expr* algfield_normalize(const Expr* e, const char* head) { (void)e; (void)head; return NULL; }

#endif
