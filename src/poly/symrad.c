/* symrad.c -- Together / Cancel over Q(variables)(Sqrt[b]), b a polynomial.
 * See symrad.h for the contract and the normal form. */

#include "symrad.h"

#include "sym_names.h"
#include "eval.h"
#include "arithmetic.h"   /* is_rational */
#include "numeric.h"      /* expr_is_numeric_quantity */
#include "core.h"         /* leaf_count_internal */
#include "flint_qqbar.h"  /* flint_qqbar_to_number_field_common */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static bool sr_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == h;
}

static bool sr_is_ratnum(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
           (sr_head(e, SYM_Rational) && e->data.function.arg_count == 2);
}

static Expr* sr_fn1(const char* h, Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1));
}
static Expr* sr_fn2(const char* h, Expr* a, Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2));
}
static Expr* sr_times(Expr* a, Expr* b) { return sr_fn2(SYM_Times, a, b); }
static Expr* sr_plus(Expr* a, Expr* b)  { return sr_fn2(SYM_Plus, a, b); }
static Expr* sr_minus(Expr* a, Expr* b) {
    return sr_plus(a, sr_times(expr_new_integer(-1), b));
}
static Expr* sr_pow(Expr* a, int64_t n) { return sr_fn2(SYM_Power, a, expr_new_integer(n)); }

static bool sr_is_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* ---- the scan --------------------------------------------------------- */

#define SR_MAX_CONST 12
typedef struct {
    const Expr* base;      /* the one symbolic radicand (borrowed)        */
    bool  bad;
    bool  constalg;        /* a constant algebraic number occurs          */
    Expr* gens[SR_MAX_CONST];   /* their generators r^(1/q) and I (owned)  */
    size_t ngens;
} SrScan;

static void sr_scan_clear(SrScan* S) {
    for (size_t i = 0; i < S->ngens; i++) expr_free(S->gens[i]);
    S->ngens = 0;
}

/* Record a constant generator (takes ownership). */
static void sr_note_gen(SrScan* S, Expr* g) {
    S->constalg = true;
    for (size_t i = 0; i < S->ngens; i++)
        if (expr_eq(S->gens[i], g)) { expr_free(g); return; }
    if (S->ngens == SR_MAX_CONST) { expr_free(g); S->bad = true; return; }
    S->gens[S->ngens++] = g;
}

/* A polynomial over Q in genuine variables: symbols that are not numeric
 * constants, rationals, Plus, Times, non-negative integer powers. */
static bool sr_poly(const Expr* e, bool* has_var) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return true;
        case EXPR_SYMBOL:
            if (expr_is_numeric_quantity(e)) return false;
            *has_var = true;
            return true;
        case EXPR_FUNCTION: break;
        default: return false;
    }
    if (sr_head(e, SYM_Rational)) return sr_is_ratnum(e);
    if (sr_head(e, SYM_Plus) || sr_head(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!sr_poly(e->data.function.args[i], has_var)) return false;
        return true;
    }
    if (sr_head(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* x = e->data.function.args[1];
        return x->type == EXPR_INTEGER && x->data.integer >= 0 &&
               sr_poly(e->data.function.args[0], has_var);
    }
    return false;
}

static void sr_scan(const Expr* e, SrScan* S) {
    if (S->bad) return;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_SYMBOL: return;
        case EXPR_FUNCTION: break;
        default: S->bad = true; return;               /* inexact */
    }
    if (sr_head(e, SYM_Rational)) return;
    if (sr_head(e, SYM_Complex) && e->data.function.arg_count == 2) {
        if (!sr_is_ratnum(e->data.function.args[0]) || !sr_is_ratnum(e->data.function.args[1]))
            S->bad = true;
        else sr_note_gen(S, expr_new_function(expr_new_symbol(SYM_Complex),
                                (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2));
        return;
    }
    if (sr_head(e, SYM_Plus) || sr_head(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) sr_scan(e->data.function.args[i], S);
        return;
    }
    if (sr_head(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER) { sr_scan(b, S); return; }
        int64_t p, q;
        if (!is_rational(x, &p, &q) || q < 2) { S->bad = true; return; }
        if (sr_is_ratnum(b)) {                              /* 2^(1/3), (-1)^(1/5) */
            sr_note_gen(S, sr_fn2(SYM_Power, expr_copy((Expr*)b),
                                  expr_new_function(expr_new_symbol(SYM_Rational),
                                      (Expr*[]){ expr_new_integer(1), expr_new_integer(q) }, 2)));
            return;
        }
        bool has_var = false;
        if (q != 2 || !sr_poly(b, &has_var) || !has_var || b->type == EXPR_SYMBOL) {
            S->bad = true;
            return;
        }
        if (!S->base) S->base = b;
        else if (!expr_eq(S->base, (Expr*)b)) S->bad = true;    /* two radicands */
        return;
    }
    S->bad = true;                                       /* another head */
}

/* ---- the pieces ------------------------------------------------------- */

static bool sr_mentions(const Expr* e, const char* name) {
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sr_mentions(e->data.function.args[i], name)) return true;
    return false;
}

/* b^(p/2) -> b^((p-1)/2) s, everything else rebuilt. */
static Expr* sr_subst(const Expr* e, const Expr* b, const char* s) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (sr_head(e, SYM_Power) && e->data.function.arg_count == 2 &&
        expr_eq(e->data.function.args[0], (Expr*)b)) {
        int64_t p, q;
        if (is_rational(e->data.function.args[1], &p, &q) && q == 2) {
            /* p odd: b^(p/2) = b^((p-1)/2) b^(1/2) */
            int64_t m = (p - 1) / 2;                   /* exact: p is odd */
            return expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy((Expr*)b), expr_new_integer(m) }, 2),
                           expr_new_symbol(s) }, 2);
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = sr_subst(e->data.function.args[i], b, s);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

/* p (a polynomial in s) modulo s^2 - b:  *c0 + *c1 s. */
static bool sr_reduce(const Expr* p, const Expr* b, const char* s, Expr** c0, Expr** c1) {
    Expr* cl = sr_fn2(SYM_CoefficientList, expr_copy((Expr*)p), expr_new_symbol(s));
    if (!sr_head(cl, SYM_List)) { expr_free(cl); return false; }
    Expr* e0 = expr_new_integer(0);
    Expr* e1 = expr_new_integer(0);
    for (size_t j = 0; j < cl->data.function.arg_count; j++) {
        Expr* t = expr_copy(cl->data.function.args[j]);
        if (sr_is_zero(t)) { expr_free(t); continue; }
        if (j >= 2) t = sr_times(t, sr_pow(expr_copy((Expr*)b), (int64_t)(j / 2)));
        if (j % 2 == 0) e0 = sr_plus(e0, t); else e1 = sr_plus(e1, t);
    }
    expr_free(cl);
    *c0 = e0; *c1 = e1;
    return true;
}

/* Is the squarefree part of b a constant, and if so is it a rational square?
 * Reads FactorSquareFreeList[b] = {{c, 1}, {p1, e1}, ...}. */
static bool sr_sqfree_part(const Expr* b, bool* constant, bool* square, Expr** cpart) {
    Expr* L = sr_fn1("FactorSquareFreeList", expr_copy((Expr*)b));
    bool ok = sr_head(L, SYM_List) && L->data.function.arg_count >= 1;
    *constant = true; *square = false;
    Expr* c = NULL;
    for (size_t i = 0; ok && i < L->data.function.arg_count; i++) {
        const Expr* pe = L->data.function.args[i];
        if (!sr_head(pe, SYM_List) || pe->data.function.arg_count != 2 ||
            pe->data.function.args[1]->type != EXPR_INTEGER) { ok = false; break; }
        const Expr* f = pe->data.function.args[0];
        int64_t ex = pe->data.function.args[1]->data.integer;
        if (sr_is_ratnum(f)) { if (ex % 2) c = expr_copy((Expr*)f); continue; }
        if (ex % 2) *constant = false;
    }
    if (ok && *constant) {
        /* a constant squarefree part c: a square iff Sqrt[c] is rational */
        Expr* r = c ? sr_fn1(SYM_Sqrt, expr_copy(c)) : expr_new_integer(1);
        *square = sr_is_ratnum(r);
        expr_free(r);
    }
    if (cpart && ok && *constant) *cpart = c ? c : expr_new_integer(1);
    else if (c) expr_free(c);
    expr_free(L);
    return ok;
}

/* [Q(gens) : Q], or 0 when the field cannot be built. */
static int64_t sr_field_degree(Expr* const* gens, size_t n) {
    if (n == 0) return 1;
    Expr* L = flint_qqbar_to_number_field_common((const Expr* const*)gens, n, 0);
    if (!L || !sr_head(L, SYM_List)) { if (L) expr_free(L); return 0; }
    int64_t d = 1;
    for (size_t i = 0; i < L->data.function.arg_count; i++) {
        const Expr* im = L->data.function.args[i];
        if (sr_head(im, SYM_AlgebraicNumber) && im->data.function.arg_count == 2 &&
            sr_head(im->data.function.args[1], SYM_List))
            d = (int64_t)im->data.function.args[1]->data.function.arg_count;
    }
    expr_free(L);
    return d;
}

static Expr* sr_canon(Expr* x, bool constalg) {          /* consumes x */
    if (!constalg) return sr_fn1(SYM_Together, x);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Together),
        (Expr*[]){ x, expr_new_function(expr_new_symbol(SYM_Rule),
                          (Expr*[]){ expr_new_symbol(SYM_Extension),
                                     expr_new_symbol(SYM_Automatic) }, 2) }, 2));
}

static Expr* sr_subst_sym(const Expr* e, const char* s, const Expr* val) {
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, s) == 0) return expr_copy((Expr*)val);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = sr_subst_sym(e->data.function.args[i], s, val);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return r;
}

/* (n0 + n1 s)/(d0 + d1 s) as A + B s with A, B in lowest terms over Q(vars)
 * (multiplying the conjugate d0 - d1 s through when d1 != 0).  False when the
 * norm d0^2 - d1^2 b vanishes (not a field).  *A and *B are fresh. */
static bool sr_coords(const Expr* n0, const Expr* n1, const Expr* d0, const Expr* d1,
                      const Expr* b, bool constalg, Expr** A, Expr** B) {
    Expr *a, *c;
    if (sr_is_zero(d1)) {
        a = sr_times(expr_copy((Expr*)n0), sr_pow(expr_copy((Expr*)d0), -1));
        c = sr_times(expr_copy((Expr*)n1), sr_pow(expr_copy((Expr*)d0), -1));
    } else {
        /* ((n0 d0 - n1 d1 b) + (n1 d0 - n0 d1) s)/(d0^2 - d1^2 b) */
        Expr* nm = sr_fn1("Expand", sr_minus(sr_pow(expr_copy((Expr*)d0), 2),
                                sr_times(sr_pow(expr_copy((Expr*)d1), 2), expr_copy((Expr*)b))));
        if (sr_is_zero(nm)) { expr_free(nm); return false; }
        Expr* ni = sr_pow(nm, -1);
        a = sr_times(sr_minus(sr_times(expr_copy((Expr*)n0), expr_copy((Expr*)d0)),
                              sr_times(sr_times(expr_copy((Expr*)n1), expr_copy((Expr*)d1)),
                                       expr_copy((Expr*)b))),
                     expr_copy(ni));
        c = sr_times(sr_minus(sr_times(expr_copy((Expr*)n1), expr_copy((Expr*)d0)),
                              sr_times(expr_copy((Expr*)n0), expr_copy((Expr*)d1))),
                     ni);
    }
    *A = sr_canon(a, constalg);
    *B = sr_canon(c, constalg);
    return true;
}

/* A + B s as one fraction (or its reciprocal), s -> Sqrt[b]. */
static Expr* sr_render(const Expr* A, const Expr* B, const char* s, const Expr* rt, bool invert) {
    Expr* x = sr_fn1(SYM_Together, sr_plus(expr_copy((Expr*)A),
                                           sr_times(expr_copy((Expr*)B), expr_new_symbol(s))));
    if (invert) x = sr_fn1(SYM_Together, sr_pow(x, -1));
    Expr* out = eval_and_free(sr_subst_sym(x, s, rt));
    expr_free(x);
    return out;
}

/* The normal form of one expression whose radicand is `b` (scan passed).
 *
 * The value is A + B Sqrt[b] with A, B in lowest terms over Q(vars) -- a
 * canonical form.  It is written as the shorter of two spellings that are
 * each a function of the value alone, so the choice is canonical too:
 * A + B Sqrt[b] as one fraction, and 1/(A' + B' Sqrt[b]) from the coordinates
 * of the reciprocal.  The second keeps a radical denominator that has nothing
 * to cancel, as Together leaves 1/(x + Sqrt[p + q]) (and 1/(x - Sqrt[2])) as
 * written, and it is the cancellation of (x - Sqrt[b])/(x^2 - b) to
 * 1/(x + Sqrt[b]). */
static Expr* sr_normal(const Expr* e, const Expr* b, bool constalg) {
    char s[32];
    int ctr = 0;
    do snprintf(s, sizeof s, "$symrad%d$", ctr++); while (sr_mentions(e, s) && ctr < 100);
    Expr* r = sr_fn1(SYM_Together, sr_subst(e, b, s));
    Expr* nn = sr_fn1(SYM_Numerator, expr_copy(r));
    Expr* dd = sr_fn1(SYM_Denominator, r);
    Expr *n0 = NULL, *n1 = NULL, *d0 = NULL, *d1 = NULL;
    Expr* out = NULL;
    Expr* rt = sr_fn2(SYM_Power, expr_copy((Expr*)b),
                      expr_new_function(expr_new_symbol(SYM_Rational),
                          (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2));
    Expr *A = NULL, *B = NULL;
    if (sr_reduce(nn, b, s, &n0, &n1) && sr_reduce(dd, b, s, &d0, &d1) &&
        sr_coords(n0, n1, d0, d1, b, constalg, &A, &B)) {
        if (sr_is_zero(A) && sr_is_zero(B)) {
            out = expr_new_integer(0);
        } else {
            out = sr_render(A, B, s, rt, false);
            Expr *A2 = NULL, *B2 = NULL;
            if (sr_coords(d0, d1, n0, n1, b, constalg, &A2, &B2)) {
                Expr* inv = sr_render(A2, B2, s, rt, true);
                if (leaf_count_internal(inv, true) < leaf_count_internal(out, true)) {
                    expr_free(out); out = inv;
                } else expr_free(inv);
                expr_free(A2); expr_free(B2);
            }
        }
    }
    if (A) expr_free(A);
    if (B) expr_free(B);
    expr_free(rt);
    if (n0) expr_free(n0);
    if (n1) expr_free(n1);
    if (d0) expr_free(d0);
    if (d1) expr_free(d1);
    expr_free(nn); expr_free(dd);
    return out;
}

Expr* symrad_normalize(const Expr* e, const char* head) {
    if (!e || (head != SYM_Together && head != SYM_Cancel)) return NULL;
    SrScan S; memset(&S, 0, sizeof S);
    sr_scan(e, &S);
    if (S.bad || !S.base) { sr_scan_clear(&S); return NULL; }
    bool constant = false, square = false;
    Expr* cpart = NULL;
    if (!sr_sqfree_part(S.base, &constant, &square, &cpart) || (constant && square)) {
        /* b a square (Sqrt[(k - 1)^2]): s^2 - b splits, no field */
        if (cpart) expr_free(cpart);
        sr_scan_clear(&S);
        return NULL;
    }
    if (constant && S.ngens > 0) {
        /* a constant squarefree part c (-3 (k - 1)^2): Sqrt[b]/(k - 1) squares
         * to c, so the ring is a field only if Sqrt[c] is not already in the
         * field of the constant algebraic numbers (it is in Q((-1)^(1/3)), and
         * there s^2 - b splits) */
        bool field = false;
        if (S.ngens < SR_MAX_CONST) {
            int64_t d0 = sr_field_degree(S.gens, S.ngens);
            S.gens[S.ngens++] = sr_fn1(SYM_Sqrt, expr_copy(cpart));
            int64_t d1 = sr_field_degree(S.gens, S.ngens);
            field = d0 > 0 && d1 == 2 * d0;
        }
        if (!field) { expr_free(cpart); sr_scan_clear(&S); return NULL; }
    }
    if (cpart) expr_free(cpart);
    sr_scan_clear(&S);
    Expr* b = expr_copy((Expr*)S.base);
    Expr* out = NULL;
    if (head == SYM_Cancel && sr_head(e, SYM_Plus)) {
        /* Cancel keeps a sum a sum: each term over the field */
        size_t n = e->data.function.arg_count;
        Expr** parts = malloc(n * sizeof(Expr*));
        bool ok = true;
        for (size_t i = 0; i < n && ok; i++) {
            SrScan T; memset(&T, 0, sizeof T);
            sr_scan(e->data.function.args[i], &T);
            sr_scan_clear(&T);
            parts[i] = T.base ? sr_normal(e->data.function.args[i], b, S.constalg)
                              : sr_canon(expr_copy(e->data.function.args[i]), S.constalg);
            if (!parts[i]) { ok = false; for (size_t j = 0; j < i; j++) expr_free(parts[j]); }
        }
        if (ok) out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), parts, n));
        free(parts);
    } else {
        out = sr_normal(e, b, S.constalg);
    }
    expr_free(b);
    return out;
}
