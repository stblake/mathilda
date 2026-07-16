/*
 * minpoly.c
 * ---------
 * Implementation of MinimalPolynomial (see minpoly.h for the user-facing
 * description).
 *
 * Pipeline for MinimalPolynomial[s, x]:
 *
 *   1. Atom walk (mp_walk): recursively rewrite the algebraic number s into a
 *      "value expression" V written in fresh auxiliary symbols, recording for
 *      each auxiliary symbol t_i a polynomial defining relation p_i (in t_i and
 *      earlier auxiliaries).  Every relation is kept polynomial — negative
 *      powers become reciprocal variables (D*w - 1) so no fractions appear.
 *   2. Build g = (x - V) and eliminate each t_i (highest index first) by
 *      g <- Resultant[g, p_i, t_i].  The introduction order guarantees t_i is
 *      present in g (in V or in a later relation already substituted in) when
 *      it is eliminated, and that each relation references only earlier
 *      auxiliaries — so the chain terminates in a univariate polynomial G(x).
 *   3. Clear denominators (Numerator[Together[G]]), make primitive over Z and
 *      factor.  Numerically evaluate s to high precision and pick the unique
 *      irreducible factor that vanishes at s.
 *   4. Return that factor, primitive with positive leading coefficient.
 *
 * Extension -> a uses the tower law: if s in Q(a) then the characteristic
 * polynomial of s over Q(a) is m_s(x)^([Q(a):Q] / [Q(s):Q]).  Membership is
 * checked via the primitive-element degree [Q(a,s):Q] == [Q(a):Q].
 *
 * Memory: the builtin takes ownership of res and returns a fresh Expr* or NULL
 * (never frees res).  The internal_* helpers consume their argument Exprs and
 * return an owned result; numericalize does not consume its input.
 */
#include "minpoly.h"
#include "sym_names.h"
#include "expr.h"
#include "internal.h"
#include "symtab.h"
#include "attr.h"
#include "numeric.h"
#include "zupoly.h"
#include "poly.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ---------------------------------------------------------------------- */
/*  Small structural helpers                                              */
/* ---------------------------------------------------------------------- */

static bool head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* True iff a symbol named `name` appears anywhere in e (head or args). */
static bool mp_has_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type == EXPR_FUNCTION) {
        if (mp_has_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (mp_has_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* Replace every Slot[1] in e with the symbol `var`; returns a fresh copy. */
static Expr* mp_subst_slot(const Expr* e, const char* var) {
    if (!e) return NULL;
    if (head_is(e, "Slot") && e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[0]->data.integer == 1) {
        return expr_new_symbol(var);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = mp_subst_slot(e->data.function.head, var);
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = mp_subst_slot(e->data.function.args[i], var);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* Replace every symbol named `from` in e with symbol `to`; fresh copy. */
static Expr* mp_subst_symbol(const Expr* e, const char* from, const char* to) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, from) == 0)
        return expr_new_symbol(to);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = mp_subst_symbol(e->data.function.head, from, to);
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = mp_subst_symbol(e->data.function.args[i], from, to);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* Replace every symbol named `from` in e with Slot[1]; fresh copy. */
static Expr* mp_subst_symbol_with_slot(const Expr* e, const char* from) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, from) == 0)
        return expr_new_function(expr_new_symbol(SYM_Slot),
                                 (Expr*[]){ expr_new_integer(1) }, 1);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = mp_subst_symbol_with_slot(e->data.function.head, from);
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = mp_subst_symbol_with_slot(e->data.function.args[i], from);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* ---------------------------------------------------------------------- */
/*  Auxiliary-variable context for the atom walk                          */
/* ---------------------------------------------------------------------- */

typedef struct {
    const Expr* orig;   /* the input s (for fresh-name collision avoidance) */
    char**      vars;   /* auxiliary symbol names, in introduction order    */
    Expr**      rels;   /* defining relation for each (owned)               */
    size_t      n, cap;
    const char* imag;   /* shared imaginary-unit var name (borrowed), or 0  */
    bool        failed; /* a non-algebraic atom was encountered             */
} AtomCtx;

static void mp_ctx_free(AtomCtx* c) {
    for (size_t i = 0; i < c->n; i++) {
        free(c->vars[i]);
        expr_free(c->rels[i]);
    }
    free(c->vars);
    free(c->rels);
}

/* A fresh symbol name not occurring in `orig` and not already allocated.
 * Returns a heap string the caller is responsible for (it is handed to
 * mp_store, which records it for later freeing). */
static char* mp_fresh_name(AtomCtx* c) {
    char buf[64];
    for (int k = 0; ; k++) {
        snprintf(buf, sizeof(buf), "$mpg%d$", k);
        if (mp_has_symbol(c->orig, buf)) continue;
        bool used = false;
        for (size_t i = 0; i < c->n; i++)
            if (strcmp(c->vars[i], buf) == 0) { used = true; break; }
        if (!used) break;
    }
    size_t len = strlen(buf);
    char* out = malloc(len + 1);
    memcpy(out, buf, len + 1);
    return out;
}

static void mp_store(AtomCtx* c, char* name, Expr* rel) {
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 8;
        c->vars = realloc(c->vars, sizeof(char*) * c->cap);
        c->rels = realloc(c->rels, sizeof(Expr*) * c->cap);
    }
    c->vars[c->n] = name;
    c->rels[c->n] = rel;
    c->n++;
}

/* True iff any auxiliary variable currently appears in e. */
static bool mp_has_auxvar(const Expr* e, const AtomCtx* c) {
    for (size_t i = 0; i < c->n; i++)
        if (mp_has_symbol(e, c->vars[i])) return true;
    return false;
}

/* The shared imaginary unit: introduced once with relation i^2 + 1. */
static const char* mp_imag(AtomCtx* c) {
    if (c->imag) return c->imag;
    char* nm = mp_fresh_name(c);
    Expr* rel = internal_plus((Expr*[]){
        internal_power((Expr*[]){ expr_new_symbol(nm), expr_new_integer(2) }, 2),
        expr_new_integer(1) }, 2);
    mp_store(c, nm, rel);
    c->imag = nm; /* borrowed alias into c->vars */
    return nm;
}

/* ---------------------------------------------------------------------- */
/*  Atom walk                                                             */
/* ---------------------------------------------------------------------- */

static Expr* mp_walk(const Expr* e, AtomCtx* c);

/* Introduce t with t^q = B^p (B already a value expression, consumed). */
static Expr* mp_make_radical(Expr* B, int64_t p, int64_t q, AtomCtx* c) {
    char* nm = mp_fresh_name(c);
    Expr* Bp = (p == 1) ? B
        : internal_power((Expr*[]){ B, expr_new_integer(p) }, 2);
    Expr* rel = internal_subtract((Expr*[]){
        internal_power((Expr*[]){ expr_new_symbol(nm), expr_new_integer(q) }, 2),
        Bp }, 2);
    mp_store(c, nm, rel);
    return expr_new_symbol(nm);
}

/* Reciprocal of a value expression (consumed).  Pure-numeric values reciprocate
 * directly; otherwise a reciprocal variable w with relation val*w - 1. */
static Expr* mp_reciprocal(Expr* val, AtomCtx* c) {
    if (!val) return NULL;
    if (!mp_has_auxvar(val, c))
        return internal_power((Expr*[]){ val, expr_new_integer(-1) }, 2);
    char* nm = mp_fresh_name(c);
    Expr* rel = internal_subtract((Expr*[]){
        internal_times((Expr*[]){ val, expr_new_symbol(nm) }, 2),
        expr_new_integer(1) }, 2);
    mp_store(c, nm, rel);
    return expr_new_symbol(nm);
}

/* Recognise Power[E, exp] as a root of unity e^(I Pi r): on success sets
 * *q_out to the (positive) denominator q of r = p/q so that the relation
 * t^(2q) - 1 vanishes at the value.  Returns false otherwise. */
static bool mp_root_of_unity_order(const Expr* exp, long* q_out) {
    const Expr* fac_single[1] = { exp };
    Expr* const* fac;
    size_t nf;
    if (head_is(exp, "Times")) {
        fac = exp->data.function.args;
        nf  = exp->data.function.arg_count;
    } else {
        fac = (Expr* const*)fac_single;
        nf  = 1;
    }
    mpq_t rat, t;
    mpq_init(rat); mpq_set_ui(rat, 1, 1);
    mpq_init(t);
    int pi = 0, im = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        const Expr* f = fac[i];
        if (f->type == EXPR_SYMBOL && strcmp(f->data.symbol.name, "Pi") == 0) {
            pi++;
        } else if (f->type == EXPR_INTEGER) {
            mpq_set_si(t, f->data.integer, 1); mpq_mul(rat, rat, t);
        } else if (head_is(f, "Rational") &&
                   f->data.function.args[0]->type == EXPR_INTEGER &&
                   f->data.function.args[1]->type == EXPR_INTEGER) {
            mpq_set_si(t, f->data.function.args[0]->data.integer,
                          (unsigned long)f->data.function.args[1]->data.integer);
            mpq_canonicalize(t); mpq_mul(rat, rat, t);
        } else if (head_is(f, "Complex") &&
                   f->data.function.args[0]->type == EXPR_INTEGER &&
                   f->data.function.args[0]->data.integer == 0) {
            const Expr* iv = f->data.function.args[1];
            if (iv->type == EXPR_INTEGER) {
                mpq_set_si(t, iv->data.integer, 1); mpq_mul(rat, rat, t);
            } else if (head_is(iv, "Rational") &&
                       iv->data.function.args[0]->type == EXPR_INTEGER &&
                       iv->data.function.args[1]->type == EXPR_INTEGER) {
                mpq_set_si(t, iv->data.function.args[0]->data.integer,
                              (unsigned long)iv->data.function.args[1]->data.integer);
                mpq_canonicalize(t); mpq_mul(rat, rat, t);
            } else ok = false;
            im++;
        } else {
            ok = false;
        }
    }
    bool matched = ok && pi == 1 && (im % 2 == 1);
    if (matched) {
        mpq_canonicalize(rat);
        *q_out = mpz_get_si(mpq_denref(rat));
        if (*q_out <= 0) matched = false;
    }
    mpq_clear(rat); mpq_clear(t);
    return matched;
}

static Expr* mp_walk(const Expr* e, AtomCtx* c) {
    if (c->failed || !e) { c->failed = true; return NULL; }

    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return expr_copy((Expr*)e);
        case EXPR_FUNCTION:
            break;                 /* handled below */
        default:                   /* Real, MPFR, String, bare Symbol */
            c->failed = true;
            return NULL;
    }

    const Expr* head = e->data.function.head;
    if (head->type != EXPR_SYMBOL) { c->failed = true; return NULL; }
    const char* h = head->data.symbol.name;
    size_t na = e->data.function.arg_count;
    Expr* const* A = e->data.function.args;

    if (strcmp(h, "Rational") == 0) return expr_copy((Expr*)e);

    if (strcmp(h, "Complex") == 0 && na == 2) {
        /* a + b*I, with a, b exact rationals. */
        if (!expr_is_numeric_like(A[0]) || !expr_is_numeric_like(A[1])) {
            c->failed = true; return NULL;
        }
        const char* inm = mp_imag(c);
        return internal_plus((Expr*[]){
            expr_copy((Expr*)A[0]),
            internal_times((Expr*[]){
                expr_copy((Expr*)A[1]), expr_new_symbol(inm) }, 2) }, 2);
    }

    if ((strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) && na >= 1) {
        Expr** vals = malloc(sizeof(Expr*) * na);
        for (size_t i = 0; i < na; i++) {
            vals[i] = mp_walk(A[i], c);
            if (c->failed) {
                for (size_t j = 0; j < i; j++) expr_free(vals[j]);
                free(vals);
                return NULL;
            }
        }
        Expr* r = (h[0] == 'P') ? internal_plus(vals, na)
                                : internal_times(vals, na);
        free(vals);
        return r;
    }

    if (strcmp(h, "Sqrt") == 0 && na == 1) {
        Expr* B = mp_walk(A[0], c);
        if (c->failed) return NULL;
        return mp_make_radical(B, 1, 2, c);
    }

    if (strcmp(h, "Power") == 0 && na == 2) {
        const Expr* base = A[0];
        const Expr* exp  = A[1];

        if (base->type == EXPR_SYMBOL && strcmp(base->data.symbol.name, "E") == 0) {
            long q;
            if (mp_root_of_unity_order(exp, &q)) {
                char* nm = mp_fresh_name(c);
                Expr* rel = internal_subtract((Expr*[]){
                    internal_power((Expr*[]){ expr_new_symbol(nm),
                        expr_new_integer(2 * (int64_t)q) }, 2),
                    expr_new_integer(1) }, 2);
                mp_store(c, nm, rel);
                return expr_new_symbol(nm);
            }
            c->failed = true; return NULL;   /* e^(transcendental) */
        }

        if (exp->type == EXPR_INTEGER) {
            int64_t n = exp->data.integer;
            if (n == 0) return expr_new_integer(1);
            if (n > 0) {
                Expr* B = mp_walk(base, c);
                if (c->failed) return NULL;
                return internal_power((Expr*[]){ B, expr_new_integer(n) }, 2);
            }
            /* n < 0: reciprocal of the positive power. */
            Expr* pos = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ expr_copy((Expr*)base), expr_new_integer(-n) }, 2);
            Expr* val = mp_walk(pos, c);
            expr_free(pos);
            if (c->failed) return NULL;
            return mp_reciprocal(val, c);
        }

        if (head_is(exp, "Rational") &&
            exp->data.function.args[0]->type == EXPR_INTEGER &&
            exp->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t p = exp->data.function.args[0]->data.integer;
            int64_t q = exp->data.function.args[1]->data.integer;
            if (q < 0) { p = -p; q = -q; }
            if (p < 0) {
                Expr* pos = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy((Expr*)base),
                        expr_new_function(expr_new_symbol(SYM_Rational),
                            (Expr*[]){ expr_new_integer(-p),
                                       expr_new_integer(q) }, 2) }, 2);
                Expr* val = mp_walk(pos, c);
                expr_free(pos);
                if (c->failed) return NULL;
                return mp_reciprocal(val, c);
            }
            Expr* B = mp_walk(base, c);
            if (c->failed) return NULL;
            return mp_make_radical(B, p, q, c);
        }

        c->failed = true; return NULL;
    }

    if (strcmp(h, "Root") == 0 && na >= 1) {
        const Expr* fn = A[0];
        if (!head_is(fn, "Function")) { c->failed = true; return NULL; }
        size_t fnn = fn->data.function.arg_count;
        char* nm = mp_fresh_name(c);
        Expr* rel = NULL;
        if (fnn == 1) {
            rel = mp_subst_slot(fn->data.function.args[0], nm);
        } else if (fnn == 2 &&
                   fn->data.function.args[0]->type == EXPR_SYMBOL) {
            rel = mp_subst_symbol(fn->data.function.args[1],
                                  fn->data.function.args[0]->data.symbol.name, nm);
        }
        if (!rel) { free(nm); c->failed = true; return NULL; }
        mp_store(c, nm, rel);
        return expr_new_symbol(nm);
    }

    c->failed = true;
    return NULL;
}

/* ---------------------------------------------------------------------- */
/*  Polynomial post-processing                                            */
/* ---------------------------------------------------------------------- */

/* Degree of an integer-coefficient polynomial in var; -1 if not integer. */
static int mp_degree(const Expr* poly, const Expr* var) {
    ZUPoly* p = expr_to_zupoly(poly, var);
    if (!p) return -1;
    int d = p->deg;
    zupoly_free(p);
    return d;
}

/* Primitive part with positive leading coefficient (integer poly), or NULL. */
static Expr* mp_primitivize(const Expr* poly, const Expr* var) {
    ZUPoly* p = expr_to_zupoly(poly, var);
    if (!p) return NULL;
    if (p->deg < 0) { zupoly_free(p); return NULL; }
    ZUPoly* pp = zupoly_primitive_part(p);
    zupoly_free(p);
    const mpz_t* lc = zupoly_getcoef(pp, pp->deg);
    if (lc && mpz_sgn(*lc) < 0) {
        ZUPoly* ng = zupoly_neg(pp);
        zupoly_free(pp);
        pp = ng;
    }
    Expr* out = zupoly_to_expr(pp, var);
    zupoly_free(pp);
    return out;
}

/* Collect distinct non-constant irreducible factors of a Factor[] result. */
static void mp_collect(const Expr* e, const Expr* var,
                       Expr*** out, size_t* n, size_t* cap) {
    if (head_is(e, "Times")) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            mp_collect(e->data.function.args[i], var, out, n, cap);
        return;
    }
    if (head_is(e, "Power") && e->data.function.arg_count == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        mp_collect(e->data.function.args[0], var, out, n, cap);
        return;
    }
    if (!mp_has_symbol(e, var->data.symbol.name)) return;   /* numeric content */
    for (size_t i = 0; i < *n; i++)
        if (expr_eq((*out)[i], e)) return;             /* already have it */
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 4;
        *out = realloc(*out, sizeof(Expr*) * (*cap));
    }
    (*out)[(*n)++] = expr_copy((Expr*)e);
}

/* Read a real numeric Expr as a double (large sentinel if non-numeric). */
static double mp_to_double(const Expr* e) {
    if (!e) return 1e300;
    switch (e->type) {
        case EXPR_INTEGER: return (double)e->data.integer;
        case EXPR_REAL:    return e->data.real;
        case EXPR_BIGINT:  return mpz_get_d(e->data.bigint);
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
        default: break;
    }
    if (head_is(e, "Rational") &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        return (double)e->data.function.args[0]->data.integer /
               (double)e->data.function.args[1]->data.integer;
    }
    return 1e300;
}

/* |f(value)| as a double, value substituted for var, evaluated under spec. */
static double mp_eval_abs(const Expr* f, const Expr* var,
                          const Expr* value, NumericSpec spec) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy((Expr*)var), expr_copy((Expr*)value) }, 2);
    Expr* sub = internal_replace_all((Expr*[]){ expr_copy((Expr*)f), rule }, 2);
    Expr* num = numericalize(sub, spec);
    expr_free(sub);
    Expr* ab = internal_abs((Expr*[]){ num }, 1);
    double d = mp_to_double(ab);
    expr_free(ab);
    return d;
}

/* ---------------------------------------------------------------------- */
/*  Core: minimal polynomial of s over Q, in `var`                        */
/* ---------------------------------------------------------------------- */

static Expr* mp_core(const Expr* s, const Expr* var) {
    AtomCtx c;
    memset(&c, 0, sizeof(c));
    c.orig = s;

    Expr* V = mp_walk(s, &c);
    if (c.failed || !V) {
        if (V) expr_free(V);
        mp_ctx_free(&c);
        return NULL;
    }

    /* g = x - V, then eliminate auxiliaries highest-index first. */
    Expr* g = internal_subtract((Expr*[]){ expr_copy((Expr*)var), V }, 2);
    bool bail = false;
    for (size_t ii = c.n; ii-- > 0; ) {
        if (!mp_has_symbol(g, c.vars[ii])) continue;
        g = internal_resultant((Expr*[]){
            g, expr_copy(c.rels[ii]), expr_new_symbol(c.vars[ii]) }, 3);
        if (g->type == EXPR_INTEGER && g->data.integer == 0) { bail = true; break; }
    }
    if (!bail)
        for (size_t i = 0; i < c.n; i++)
            if (mp_has_symbol(g, c.vars[i])) { bail = true; break; }
    if (bail) { expr_free(g); mp_ctx_free(&c); return NULL; }

    /* Univariate G(x); clear denominators and make it a Z-polynomial. */
    Expr* G = internal_expand((Expr*[]){ g }, 1);
    G = internal_together((Expr*[]){ G }, 1);
    G = internal_numerator((Expr*[]){ G }, 1);
    G = internal_expand((Expr*[]){ G }, 1);

    /* Factor and gather candidate irreducible factors. */
    Expr* fac = internal_factor((Expr*[]){ expr_copy(G) }, 1);
    Expr** cands = NULL;
    size_t ncand = 0, capc = 0;
    mp_collect(fac, var, &cands, &ncand, &capc);
    expr_free(fac);
    if (ncand == 0) {                       /* irreducible already, or constant */
        cands = malloc(sizeof(Expr*));
        cands[0] = expr_copy(G);
        ncand = 1;
    }

    /* Pick the factor that vanishes at s (high-precision numeric test). */
    NumericSpec spec;
#ifdef USE_MPFR
    spec.mode = NUMERIC_MODE_MPFR;
    spec.bits = numeric_digits_to_bits(80);
#else
    spec = numeric_machine_spec();
#endif
    Expr* Ns = numericalize(s, spec);
    int best = -1;
    double bestmag = 1e300;
    for (size_t i = 0; i < ncand; i++) {
        double m = mp_eval_abs(cands[i], var, Ns, spec);
        if (m < bestmag) { bestmag = m; best = (int)i; }
    }
    expr_free(Ns);

    Expr* chosen = (best >= 0) ? expr_copy(cands[best]) : NULL;
    for (size_t i = 0; i < ncand; i++) expr_free(cands[i]);
    free(cands);
    expr_free(G);
    mp_ctx_free(&c);

    /* eps: with high precision a true root is essentially 0; non-roots are
     * O(1).  A loose ceiling guards against a non-algebraic slip-through. */
    if (!chosen || bestmag > 1e-3) {
        if (chosen) expr_free(chosen);
        return NULL;
    }

    Expr* prim = mp_primitivize(chosen, var);
    if (prim) { expr_free(chosen); return prim; }
    return chosen;
}

/* ---------------------------------------------------------------------- */
/*  Extension -> a : characteristic polynomial of s in Q(a) over Q(a)     */
/* ---------------------------------------------------------------------- */

/* [Q(a,s):Q] estimated as the max minimal-polynomial degree of a + r*s over
 * a few multipliers r (primitive-element theorem); -1 on failure. */
static int mp_compositum_degree(const Expr* a, const Expr* s, const Expr* var) {
    int best = -1;
    for (int r = 1; r <= 2; r++) {
        Expr* gamma = internal_plus((Expr*[]){
            expr_copy((Expr*)a),
            internal_times((Expr*[]){ expr_new_integer(r),
                                      expr_copy((Expr*)s) }, 2) }, 2);
        Expr* mg = mp_core(gamma, var);
        expr_free(gamma);
        if (mg) {
            int dd = mp_degree(mg, var);
            expr_free(mg);
            if (dd > best) best = dd;
        }
    }
    return best;
}

static Expr* mp_extension(const Expr* s, const Expr* x, const Expr* a) {
    Expr* ms = mp_core(s, x);
    if (!ms) return NULL;
    int e = mp_degree(ms, x);
    Expr* ma = mp_core(a, x);
    if (!ma) { expr_free(ms); return NULL; }
    int d = mp_degree(ma, x);
    expr_free(ma);

    int comp = mp_compositum_degree(a, s, x);
    if (e <= 0 || d <= 0 || comp != d || d % e != 0) {
        fprintf(stderr, "MinimalPolynomial::ext: the second argument is not "
                        "an element of the field extension generated by the "
                        "Extension option.\n");
        expr_free(ms);
        return NULL;
    }
    int k = d / e;
    if (k == 1) return ms;

    Expr* pw  = internal_power((Expr*[]){ ms, expr_new_integer(k) }, 2);
    Expr* ex  = internal_expand((Expr*[]){ pw }, 1);
    Expr* prim = mp_primitivize(ex, x);
    if (prim) { expr_free(ex); return prim; }
    return ex;
}

/* ---------------------------------------------------------------------- */
/*  Builtin dispatch                                                      */
/* ---------------------------------------------------------------------- */

Expr* builtin_minimalpolynomial(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr* const* args = res->data.function.args;

    /* Separate positional args from an Extension -> a option. */
    Expr* ext = NULL;
    Expr* pos[3];
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        const Expr* a = args[i];
        if ((head_is(a, "Rule") || head_is(a, "RuleDelayed")) &&
            a->data.function.arg_count == 2 &&
            a->data.function.args[0]->type == EXPR_SYMBOL &&
            strcmp(a->data.function.args[0]->data.symbol.name, "Extension") == 0) {
            ext = a->data.function.args[1];
        } else if (npos < 3) {
            pos[npos++] = args[i];
        } else {
            npos++;     /* too many; falls through to arity error */
        }
    }

    if (npos == 1 && !ext) {
        /* Pure-function form: MinimalPolynomial[s]. */
        const Expr* s = pos[0];
        char* nm = poly_make_fresh_gen((Expr*)s);
        Expr* var = expr_new_symbol(nm);
        Expr* p = mp_core(s, var);
        expr_free(var);
        if (!p) { free(nm); return NULL; }
        Expr* body = mp_subst_symbol_with_slot(p, nm);
        expr_free(p);
        free(nm);
        return expr_new_function(expr_new_symbol(SYM_Function),
                                 (Expr*[]){ body }, 1);
    }

    if (npos == 2) {
        const Expr* s = pos[0];
        const Expr* x = pos[1];
        if (x->type != EXPR_SYMBOL) {
            fprintf(stderr, "MinimalPolynomial::ivar: the second argument %s "
                            "is not a valid variable.\n", "");
            return NULL;
        }
        if (ext) return mp_extension(s, x, ext);
        return mp_core(s, x);
    }

    fprintf(stderr, "MinimalPolynomial::argt: MinimalPolynomial called with an "
                    "unsupported number of arguments; 1 or 2 (plus an optional "
                    "Extension rule) are expected.\n");
    return NULL;
}

void minpoly_init(void) {
    symtab_add_builtin("MinimalPolynomial", builtin_minimalpolynomial);
    symtab_get_def("MinimalPolynomial")->attributes |=
        ATTR_PROTECTED | ATTR_LISTABLE;
}
