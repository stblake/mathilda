/*
 * print_lean.c -- LeanForm renderer and the strict Lean translator
 * (see print_lean.h).
 *
 * One recursive printer, build(), writes the Lean text of a node into a
 * growable buffer and returns the Lean precedence of what it wrote; put()
 * parenthesises that text when the context needs a tighter binding.  The
 * display form and the strict translator differ only in the context flags:
 * in strict mode an untranslatable node records the first refusal (reason +
 * node) and the caller discards the text.
 *
 * Rendering rules (one line each; print_lean.h has the precedence table):
 *   - integers as decimal; a negative literal is a negation, so it is
 *     parenthesised as an operand: 2 * (-3);
 *   - Rational[p, q] -> (p / q), or (p / q : T) in strict mode with a number
 *     type; refused under Int;
 *   - Plus renders negative terms as subtraction (a - 2 * b) and, when the
 *     canonical order puts a negative term first, moves the first
 *     non-negative term to the front (x - 3, not -3 + x);
 *   - Times keeps a rational coefficient as a factor ((3 / 4) * x), moves
 *     Power[_, -k] factors to a denominator (x / (y * z)), and prints a
 *     negative coefficient as a negation of the rest: -(2 * x), -x;
 *   - Power[a, n] (n >= 0 an integer) -> a ^ n; Power[a, -n] -> 1 / a ^ n;
 *     Power[a, 1/2] -> Real.sqrt a; Power[a, k] with k a symbol whose Lean
 *     name is in the pass-through set -> a ^ k (the generator vouches that k
 *     is a Nat-typed binder, so `^` is the monoid power); any other exponent
 *     only in display form;
 *   - relations = != < <= > >=, chains and Inequality[...] become
 *     conjunctions; And Or Not Implies Equivalent -> /\ \/ not -> <->;
 *   - ForAll / Exists -> binders (typed with the number type, if any);
 *     Function[x, b] -> fun x => b; List -> [a, b];
 *   - Log Exp Sin Cos Tan ArcTan Sinh Cosh, Pi, E, I -> Real.* / Complex.I
 *     (refused when core_only);
 *   - any other head f[x, y] -> the application f x y (strict mode: only for
 *     heads in the pass-through set);
 *   - symbols: a plain identifier passes through; anything else (a context
 *     mark, $, a Lean keyword) is escaped as <<name>> (French quotes).  The
 *     strict translator refuses a name that contains a French quote or a
 *     control character, two distinct symbols with the same Lean text, and --
 *     with LeanFormOpts.closed -- any symbol not bound by an enclosing
 *     quantifier / Function and not in the pass-through set (a free name such
 *     as `default` would silently resolve to a Lean global).
 */
#include "print_lean.h"

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "context.h"
#include "message.h"
#include "ndarray.h"
#include "print.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ */
/* Output buffer                                                       */
/* ------------------------------------------------------------------ */

typedef struct { char* s; size_t n, cap; } LBuf;

static void lb_putn(LBuf* b, const char* s, size_t n) {
    if (b->n + n + 1 > b->cap) {
        size_t nc = b->cap ? b->cap : 64;
        while (nc < b->n + n + 1) nc *= 2;
        char* p = realloc(b->s, nc);
        if (!p) abort();
        b->s = p;
        b->cap = nc;
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = '\0';
}

static void lb_puts(LBuf* b, const char* s) { lb_putn(b, s, strlen(s)); }

static char* lb_take(LBuf* b) {
    if (!b->s) { char* e = malloc(1); if (!e) abort(); e[0] = '\0'; return e; }
    return b->s;
}

/* ------------------------------------------------------------------ */
/* Context and precedences                                             */
/* ------------------------------------------------------------------ */

enum {
    P_BINDER = 0,     /* fun / forall / exists: parenthesised unless top level */
    P_IFF = 20,
    P_IMP = 25,
    P_OR = 30,
    P_AND = 35,
    P_NOT = 40,
    P_REL = 50,
    P_PLUS = 65,
    P_NEG = 66,       /* -x and negative literals: parenthesised as operands */
    P_TIMES = 70,
    P_POW = 75,
    P_APP = 1024,
    P_ATOM = 2000
};

#define FMT_DEPTH_LIMIT 64

typedef struct {
    bool strict;
    LeanFormOpts o;
    const char* err;      /* first refusal (strict mode) */
    const Expr* bad;
    int fmt_depth;
    /* Strict mode only: the binder variables in scope (interned names) ... */
    const char** bound;
    size_t n_bound, cap_bound;
    /* ... and every identifier written so far, as {Mathilda name, Lean text},
     * so two distinct symbols never print as the same Lean identifier. */
    const char** seen_name;
    char** seen_text;
    size_t n_seen, cap_seen;
} LCtx;

static void ctx_release(LCtx* c) {
    for (size_t i = 0; i < c->n_seen; i++) free(c->seen_text[i]);
    free((void*)c->seen_name);
    free(c->seen_text);
    free((void*)c->bound);
}

static int build(LCtx* c, const Expr* e, LBuf* b);

static int refuse(LCtx* c, const Expr* e, const char* why, LBuf* b) {
    if (!c->err) { c->err = why; c->bad = e; }
    lb_puts(b, "sorry");
    return P_ATOM;
}

/* e at a position that needs binding strength `need`. */
static void put(LCtx* c, const Expr* e, int need, LBuf* b) {
    LBuf t = {0};
    int p = build(c, e, &t);
    if (p < need) lb_puts(b, "(");
    if (t.s) lb_putn(b, t.s, t.n);
    if (p < need) lb_puts(b, ")");
    free(t.s);
}

/* ------------------------------------------------------------------ */
/* Predicates on numbers                                               */
/* ------------------------------------------------------------------ */

static const char* head_name(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    return h->type == EXPR_SYMBOL ? h->data.symbol.name : NULL;
}

static bool head_is(const Expr* e, const char* sym, size_t argc) {
    return head_name(e) == sym && e->data.function.arg_count == argc;
}

static bool is_int(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

static bool is_rational(const Expr* e) {
    return head_is(e, SYM_Rational, 2) && is_int(e->data.function.args[0]) &&
           is_int(e->data.function.args[1]);
}

static bool is_real_atom(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool is_number(const Expr* e) {
    return is_int(e) || is_rational(e) || is_real_atom(e) || head_is(e, SYM_Complex, 2);
}

static int int_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (e->data.integer > 0) - (e->data.integer < 0);
    return mpz_sgn(e->data.bigint);
}

/* Sign of a real number (0 for a complex number or a non-number). */
static int real_sign(const Expr* e) {
    if (is_int(e)) return int_sign(e);
    if (is_rational(e)) return int_sign(e->data.function.args[0]);
    if (e->type == EXPR_REAL) return (e->data.real > 0) - (e->data.real < 0);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_sgn(e->data.mpfr);
#endif
    return 0;
}

static bool is_int_value(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* -e for a real number, as a fresh node. */
static Expr* num_negate(const Expr* e) {
    if (e->type == EXPR_INTEGER && e->data.integer != INT64_MIN)
        return expr_new_integer(-e->data.integer);
    if (is_int(e)) {
        mpz_t z; mpz_init(z); expr_to_mpz((Expr*)e, z); mpz_neg(z, z);
        Expr* r = expr_new_bigint_from_mpz(z); mpz_clear(z);
        return r;
    }
    if (e->type == EXPR_REAL) return expr_new_real(-e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        Expr* r = expr_new_mpfr_bits(mpfr_get_prec(e->data.mpfr));
        mpfr_neg(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
        return r;
    }
#endif
    if (is_rational(e)) {
        Expr* a[2] = { num_negate(e->data.function.args[0]),
                       expr_copy(e->data.function.args[1]) };
        return expr_new_function(expr_new_symbol(SYM_Rational), a, 2);
    }
    return expr_copy((Expr*)e);
}

/* A term that Plus prints as a subtraction: a negative real number, or a
 * product whose leading numeric coefficient is negative. */
static bool is_negative_term(const Expr* t) {
    if (real_sign(t) < 0) return true;
    if (head_name(t) == SYM_Times && t->data.function.arg_count >= 2)
        return real_sign(t->data.function.args[0]) < 0;
    return false;
}

/* -t for a negative term, as a fresh node (Times[-1, x] -> x). */
static Expr* negate_term(const Expr* t) {
    if (is_number(t)) return num_negate(t);
    size_t n = t->data.function.arg_count;
    Expr* c = num_negate(t->data.function.args[0]);
    bool unit = is_int_value(c, 1);
    if (unit && n == 2) { expr_free(c); return expr_copy(t->data.function.args[1]); }
    size_t m = unit ? n - 1 : n;
    Expr** a = malloc(m * sizeof *a);
    if (!a) abort();
    size_t k = 0;
    if (unit) expr_free(c); else a[k++] = c;
    for (size_t i = 1; i < n; i++) a[k++] = expr_copy(t->data.function.args[i]);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), a, m);
    free(a);
    return r;
}

/* An expression built only from numbers and arithmetic heads: a relation
 * between such sides gets a type ascription so numerals do not default to Nat. */
static bool is_closed_numeric(const Expr* e) {
    if (is_number(e)) return true;
    const char* h = head_name(e);
    if (h != SYM_Plus && h != SYM_Times && h != SYM_Power) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!is_closed_numeric(e->data.function.args[i])) return false;
    return true;
}

/* ------------------------------------------------------------------ */
/* Atoms                                                               */
/* ------------------------------------------------------------------ */

static void put_int(const Expr* e, LBuf* b) {
    if (e->type == EXPR_INTEGER) {
        char t[32];
        snprintf(t, sizeof t, "%" PRId64, e->data.integer);
        lb_puts(b, t);
        return;
    }
    char* s = malloc(mpz_sizeinbase(e->data.bigint, 10) + 2);
    if (!s) abort();
    mpz_get_str(s, 10, e->data.bigint);
    lb_puts(b, s);
    free(s);
}

static void put_double(double d, LBuf* b) {
    char t[40];
    snprintf(t, sizeof t, "%.17g", d);
    lb_puts(b, t);
    if (!strpbrk(t, ".eEn")) lb_puts(b, ".0");   /* 100 -> 100.0; n: nan / inf */
}

static const char* num_type_name(LeanNumType t) {
    switch (t) {
    case LEAN_NUM_INT: return "Int";
    case LEAN_NUM_RAT: return "Rat";
    case LEAN_NUM_REAL: return "Real";
    default: return NULL;
    }
}

static bool division_allowed(const LCtx* c) {
    if (!c->strict) return true;
    return c->o.num != LEAN_NUM_INT && c->o.allow_division;
}

static const char* division_refusal(const LCtx* c) {
    return c->o.num == LEAN_NUM_INT ? "division over Int" : "division not allowed";
}

/* A real number (integer, rational, machine or MPFR real). */
static int build_number(LCtx* c, const Expr* e, LBuf* b) {
    if (real_sign(e) < 0) {
        Expr* m = num_negate(e);
        lb_puts(b, "-");
        put(c, m, P_ATOM, b);
        expr_free(m);
        return P_NEG;
    }
    if (is_int(e)) { put_int(e, b); return P_ATOM; }
    if (is_rational(e)) {
        if (!division_allowed(c)) return refuse(c, e, division_refusal(c), b);
        const char* ty = c->strict ? num_type_name(c->o.num) : NULL;
        lb_puts(b, "(");
        put_int(e->data.function.args[0], b);
        lb_puts(b, " / ");
        put_int(e->data.function.args[1], b);
        if (ty) { lb_puts(b, " : "); lb_puts(b, ty); }
        lb_puts(b, ")");
        return P_ATOM;
    }
    if (c->strict) return refuse(c, e, "inexact number; Rationalize it first", b);
    if (e->type == EXPR_REAL) { put_double(e->data.real, b); return P_ATOM; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) { put_double(mpfr_get_d(e->data.mpfr, MPFR_RNDN), b); return P_ATOM; }
#endif
    return P_ATOM;
}

static const char* const LEAN_KEYWORDS[] = {
    "abbrev", "at", "attribute", "axiom", "break", "by", "calc", "catch", "class",
    "continue", "decreasing_by", "def", "deriving", "do", "else", "end", "example",
    "exists", "export", "extends", "finally", "for", "forall", "fun", "have", "if",
    "import", "in", "inductive", "infix", "infixl", "infixr", "instance", "let",
    "local", "macro", "match", "mut", "mutual", "namespace", "noncomputable",
    "nomatch", "nofun", "notation", "obtain", "open", "partial", "postfix",
    "prefix", "private", "protected", "rec", "return", "section", "set_option",
    "show", "structure", "suffices", "syntax", "termination_by", "then", "theorem",
    "try", "unless", "universe", "unsafe", "variable", "where", "with", "Prop",
    "Sort", "Type", NULL
};

static bool is_plain_ident(const char* s) {
    if (!s || !((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z'))) return false;
    for (const char* p = s; *p; p++)
        if (!((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z') || (*p >= '0' && *p <= '9')))
            return false;
    for (size_t i = 0; LEAN_KEYWORDS[i]; i++)
        if (strcmp(s, LEAN_KEYWORDS[i]) == 0) return false;
    return true;
}

/* A Mathilda symbol name as a Lean identifier: plain, or escaped <<name>>. */
static void put_ident(const char* name, LBuf* b) {
    const char* nm = context_display_name(name);
    if (is_plain_ident(nm)) { lb_puts(b, nm); return; }
    lb_puts(b, "\xC2\xAB");          /* U+00AB << */
    lb_puts(b, nm);
    lb_puts(b, "\xC2\xBB");          /* U+00BB >> */
}

/* Strict mode: may `name` be written as an identifier?  Refused when its
 * escaped form could break out of the French quotes (a << or >> or a control
 * character inside the name), or when another symbol already printed as the
 * same Lean text (two Mathilda symbols must stay two Lean variables). */
static const char* ident_problem(LCtx* c, const char* name) {
    if (!c->strict) return NULL;
    const char* nm = context_display_name(name);
    for (const unsigned char* p = (const unsigned char*)nm; *p; p++) {
        if (*p < 0x20 || *p == 0x7F) return "symbol name with a control character";
        if (p[0] == 0xC2 && (p[1] == 0xAB || p[1] == 0xBB))
            return "symbol name containing a French quote";
    }
    LBuf t = {0};
    put_ident(name, &t);
    const char* txt = t.s ? t.s : "";
    for (size_t i = 0; i < c->n_seen; i++) {
        if (strcmp(c->seen_text[i], txt) == 0) {
            free(t.s);
            return c->seen_name[i] == name ? NULL
                                           : "two symbols with the same Lean name";
        }
    }
    if (c->n_seen == c->cap_seen) {
        size_t nc = c->cap_seen ? 2 * c->cap_seen : 16;
        const char** a = realloc((void*)c->seen_name, nc * sizeof *a);
        if (!a) abort();
        c->seen_name = a;
        char** tx = realloc(c->seen_text, nc * sizeof *tx);
        if (!tx) abort();
        c->seen_text = tx;
        c->cap_seen = nc;
    }
    c->seen_name[c->n_seen] = name;
    c->seen_text[c->n_seen] = t.s ? t.s : lb_take(&t);
    c->n_seen++;
    return NULL;
}

static bool is_bound(const LCtx* c, const char* name) {
    for (size_t i = c->n_bound; i > 0; i--)
        if (c->bound[i - 1] == name) return true;
    return false;
}

static void push_bound(LCtx* c, const char* name) {
    if (c->n_bound == c->cap_bound) {
        size_t nc = c->cap_bound ? 2 * c->cap_bound : 8;
        const char** a = realloc((void*)c->bound, nc * sizeof *a);
        if (!a) abort();
        c->bound = a;
        c->cap_bound = nc;
    }
    c->bound[c->n_bound++] = name;
}

/* Push the binder symbols of v (x or {x, y}); the number pushed. */
static size_t push_binders(LCtx* c, const Expr* v) {
    if (v->type == EXPR_SYMBOL) { push_bound(c, v->data.symbol.name); return 1; }
    size_t k = 0;
    if (head_name(v) == SYM_List)
        for (size_t i = 0; i < v->data.function.arg_count; i++)
            if (v->data.function.args[i]->type == EXPR_SYMBOL) {
                push_bound(c, v->data.function.args[i]->data.symbol.name);
                k++;
            }
    return k;
}

static const char* pass_head(const LCtx* c, const char* h, char* buf, size_t bufsz);

static int build_symbol(LCtx* c, const Expr* e, LBuf* b) {
    const char* nm = e->data.symbol.name;
    if (nm == SYM_True) { lb_puts(b, "True"); return P_ATOM; }
    if (nm == SYM_False) { lb_puts(b, "False"); return P_ATOM; }
    if (nm == SYM_Pi || nm == SYM_E || nm == SYM_I) {
        if (c->strict && c->o.core_only)
            return refuse(c, e, "real/complex constant needs Mathlib (Tier 1)", b);
        if (nm == SYM_Pi) { lb_puts(b, "Real.pi"); return P_ATOM; }
        if (nm == SYM_I) { lb_puts(b, "Complex.I"); return P_ATOM; }
        lb_puts(b, "Real.exp 1");
        return P_APP;
    }
    if (c->strict && (nm == SYM_Infinity || nm == SYM_ComplexInfinity || nm == SYM_Indeterminate))
        return refuse(c, e, "infinite or indeterminate quantity", b);
    if (c->strict && c->o.closed && !is_bound(c, nm)) {
        char buf[256];
        if (!pass_head(c, nm, buf, sizeof buf))
            return refuse(c, e, "free variable (not bound by a quantifier)", b);
    }
    const char* why = ident_problem(c, nm);
    if (why) return refuse(c, e, why, b);
    put_ident(nm, b);
    return P_ATOM;
}

static void put_string_literal(const char* s, LBuf* b) {
    lb_puts(b, "\"");
    for (const char* p = s; *p; p++) {
        switch (*p) {
        case '"': lb_puts(b, "\\\""); break;
        case '\\': lb_puts(b, "\\\\"); break;
        case '\n': lb_puts(b, "\\n"); break;
        case '\t': lb_puts(b, "\\t"); break;
        case '\r': lb_puts(b, "\\r"); break;
        default: lb_putn(b, p, 1);
        }
    }
    lb_puts(b, "\"");
}

/* ------------------------------------------------------------------ */
/* Arithmetic                                                          */
/* ------------------------------------------------------------------ */

static int build_plus(LCtx* c, const Expr* e, LBuf* b) {
    size_t n = e->data.function.arg_count;
    Expr* const* t = e->data.function.args;
    /* Lead with the first non-negative term: x - 3 rather than -3 + x. */
    size_t lead = 0;
    if (is_negative_term(t[0]))
        for (size_t i = 1; i < n; i++)
            if (!is_negative_term(t[i])) { lead = i; break; }
    put(c, t[lead], P_PLUS, b);
    for (size_t i = 0; i < n; i++) {
        if (i == lead) continue;
        if (is_negative_term(t[i])) {
            Expr* m = negate_term(t[i]);
            lb_puts(b, " - ");
            put(c, m, P_PLUS + 1, b);
            expr_free(m);
        } else {
            lb_puts(b, " + ");
            put(c, t[i], P_PLUS + 1, b);
        }
    }
    return P_PLUS;
}

/* Power[_, negative integer] or Power[_, -1/2]: a denominator factor. */
static bool is_denominator_factor(const Expr* f) {
    if (!head_is(f, SYM_Power, 2)) return false;
    const Expr* x = f->data.function.args[1];
    if (is_int(x)) return int_sign(x) < 0;
    return is_rational(x) && is_int_value(x->data.function.args[0], -1) &&
           is_int_value(x->data.function.args[1], 2);
}

/* The denominator factor b^k of b^-k (or Sqrt[b] of b^(-1/2)), fresh. */
static Expr* invert_factor(const Expr* f) {
    Expr* a[2] = { expr_copy(f->data.function.args[0]), num_negate(f->data.function.args[1]) };
    if (is_int_value(a[1], 1)) { expr_free(a[1]); return a[0]; }
    return expr_new_function(expr_new_symbol(SYM_Power), a, 2);
}

/* f1 * f2 * ... (left-associative), the factors given. */
static void put_product(LCtx* c, Expr* const* fs, size_t n, LBuf* b) {
    for (size_t i = 0; i < n; i++) {
        if (i) lb_puts(b, " * ");
        put(c, fs[i], i ? P_TIMES + 1 : P_TIMES, b);
    }
}

static int build_times(LCtx* c, const Expr* e, LBuf* b) {
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    if (n == 1) return build(c, a[0], b);
    if (real_sign(a[0]) < 0) {
        Expr* m = negate_term(e);
        lb_puts(b, "-");
        put(c, m, P_ATOM, b);
        expr_free(m);
        return P_NEG;
    }
    Expr** num = malloc(n * sizeof *num);
    Expr** den = malloc(n * sizeof *den);
    if (!num || !den) abort();
    size_t nn = 0, nd = 0;
    for (size_t i = 0; i < n; i++) {
        if (i == 0 && is_int_value(a[i], 1)) continue;
        if (is_denominator_factor(a[i])) den[nd++] = invert_factor(a[i]);
        else num[nn++] = (Expr*)a[i];
    }
    int prec = P_TIMES;
    if (nd && !division_allowed(c)) {
        prec = refuse(c, e, division_refusal(c), b);
    } else {
        if (nn) put_product(c, num, nn, b);
        else lb_puts(b, "1");
        if (nd) {
            lb_puts(b, " / ");
            if (nd == 1) put(c, den[0], P_TIMES + 1, b);
            else { lb_puts(b, "("); put_product(c, den, nd, b); lb_puts(b, ")"); }
        }
    }
    for (size_t i = 0; i < nd; i++) expr_free(den[i]);
    free(num);
    free(den);
    return prec;
}

/* A symbol exponent the generator declared Nat-typed: its Lean name (context
 * marks as dots, as for pass-through heads) is in the pass-through set. */
static bool nat_exponent_symbol(const LCtx* c, const Expr* x) {
    char buf[256];
    return x->type == EXPR_SYMBOL && pass_head(c, x->data.symbol.name, buf, sizeof buf) != NULL;
}

static int build_power(LCtx* c, const Expr* e, LBuf* b) {
    const Expr* base = e->data.function.args[0];
    const Expr* x = e->data.function.args[1];
    bool half = is_rational(x) && is_int_value(x->data.function.args[1], 2) &&
                (is_int_value(x->data.function.args[0], 1) ||
                 is_int_value(x->data.function.args[0], -1));
    if (half) {
        if (c->strict && c->o.core_only) return refuse(c, e, "radical needs Mathlib (Tier 1)", b);
        bool inv = is_int_value(x->data.function.args[0], -1);
        if (inv) {
            if (!division_allowed(c)) return refuse(c, e, division_refusal(c), b);
            lb_puts(b, "1 / ");
        }
        lb_puts(b, "Real.sqrt ");
        put(c, base, P_ATOM, b);
        return inv ? P_TIMES : P_APP;
    }
    if (is_int(x) && int_sign(x) < 0) {
        if (!division_allowed(c)) return refuse(c, e, division_refusal(c), b);
        Expr* d = invert_factor(e);
        lb_puts(b, "1 / ");
        put(c, d, P_TIMES + 1, b);
        expr_free(d);
        return P_TIMES;
    }
    if (c->strict && !is_int(x) && !nat_exponent_symbol(c, x))
        return refuse(c, e, "non-integer exponent", b);
    put(c, base, P_POW + 1, b);
    lb_puts(b, " ^ ");
    put(c, x, P_POW, b);
    return P_POW;
}

/* ------------------------------------------------------------------ */
/* Relations, logic, binders                                           */
/* ------------------------------------------------------------------ */

static const char* rel_op(const char* h) {
    if (h == SYM_Equal) return " = ";
    if (h == SYM_Unequal) return " \xE2\x89\xA0 ";        /* U+2260 */
    if (h == SYM_Less) return " < ";
    if (h == SYM_LessEqual) return " \xE2\x89\xA4 ";      /* U+2264 */
    if (h == SYM_Greater) return " > ";
    if (h == SYM_GreaterEqual) return " \xE2\x89\xA5 ";   /* U+2265 */
    return NULL;
}

/* lhs op rhs; a closed numeric relation ascribes its left side in strict mode. */
static void put_relation(LCtx* c, const Expr* l, const char* op, const Expr* r, LBuf* b) {
    const char* ty = c->strict ? num_type_name(c->o.num) : NULL;
    if (ty && is_closed_numeric(l) && is_closed_numeric(r)) {
        lb_puts(b, "(");
        put(c, l, P_BINDER, b);
        lb_puts(b, " : ");
        lb_puts(b, ty);
        lb_puts(b, ")");
    } else {
        put(c, l, P_REL + 1, b);
    }
    lb_puts(b, op);
    put(c, r, P_REL + 1, b);
}

static int build_relation(LCtx* c, const Expr* e, const char* op, LBuf* b) {
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    if (n == 2) { put_relation(c, a[0], op, a[1], b); return P_REL; }
    /* Unequal[a, b, c]: pairwise distinct; the others chain. */
    bool first = true;
    for (size_t i = 0; i + 1 < n; i++) {
        size_t jmax = head_name(e) == SYM_Unequal ? n : i + 2;
        for (size_t j = i + 1; j < jmax; j++) {
            if (!first) lb_puts(b, " \xE2\x88\xA7 ");        /* U+2227 /\ */
            first = false;
            LBuf t = {0};
            put_relation(c, a[i], op, a[j], &t);
            lb_puts(b, t.s);
            free(t.s);
        }
    }
    return P_AND;
}

/* Inequality[a, Less, b, LessEqual, c, ...] as a conjunction of relations. */
static int build_inequality(LCtx* c, const Expr* e, LBuf* b) {
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    for (size_t i = 0; i + 2 < n; i += 2) {
        const char* op = a[i + 1]->type == EXPR_SYMBOL ? rel_op(a[i + 1]->data.symbol.name) : NULL;
        if (!op) return refuse(c, e, "malformed Inequality", b);
        if (i) lb_puts(b, " \xE2\x88\xA7 ");
        put_relation(c, a[i], op, a[i + 2], b);
    }
    return n > 3 ? P_AND : P_REL;
}

static int build_nary(LCtx* c, const Expr* e, const char* op, int prec, LBuf* b) {
    size_t n = e->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (i) lb_puts(b, op);
        put(c, e->data.function.args[i], prec + 1, b);
    }
    return prec;
}

/* The binder variables of ForAll / Exists / Function: x or {x, y}. */
static bool put_binder_vars(LCtx* c, const Expr* v, LBuf* b) {
    if (v->type == EXPR_SYMBOL) {
        if (ident_problem(c, v->data.symbol.name)) return false;
        put_ident(v->data.symbol.name, b);
        return true;
    }
    if (head_name(v) != SYM_List || v->data.function.arg_count == 0) return false;
    for (size_t i = 0; i < v->data.function.arg_count; i++) {
        const Expr* s = v->data.function.args[i];
        if (s->type != EXPR_SYMBOL) return false;
        if (ident_problem(c, s->data.symbol.name)) return false;
        if (i) lb_puts(b, " ");
        put_ident(s->data.symbol.name, b);
    }
    return true;
}

/* ForAll[x, body], ForAll[x, cond, body] (and Exists). */
static int build_quantifier(LCtx* c, const Expr* e, bool forall, LBuf* b) {
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;
    LBuf vars = {0};
    if (!put_binder_vars(c, a[0], &vars)) {
        free(vars.s);
        return refuse(c, e, "binder is not a symbol or a list of symbols", b);
    }
    lb_puts(b, forall ? "\xE2\x88\x80 " : "\xE2\x88\x83 ");   /* U+2200 / U+2203 */
    lb_puts(b, vars.s);
    free(vars.s);
    const char* ty = num_type_name(c->o.num);
    if (ty) { lb_puts(b, " : "); lb_puts(b, ty); }
    lb_puts(b, ", ");
    size_t pushed = push_binders(c, a[0]);
    if (n == 3) {
        put(c, a[1], forall ? P_IMP + 1 : P_AND + 1, b);
        lb_puts(b, forall ? " \xE2\x86\x92 " : " \xE2\x88\xA7 ");   /* -> or /\ */
        put(c, a[2], forall ? P_IMP : P_AND + 1, b);
    } else {
        put(c, a[1], P_BINDER, b);
    }
    c->n_bound -= pushed;
    return P_BINDER;
}

static int build_function(LCtx* c, const Expr* e, LBuf* b) {
    Expr* const* a = e->data.function.args;
    if (e->data.function.arg_count != 2) {
        if (c->strict) return refuse(c, e, "pure function with Slot", b);
        return -1;    /* display: fall back to the application form */
    }
    LBuf vars = {0};
    if (!put_binder_vars(c, a[0], &vars)) {
        free(vars.s);
        if (c->strict) return refuse(c, e, "Function parameters are not symbols", b);
        return -1;
    }
    lb_puts(b, "fun ");
    lb_puts(b, vars.s);
    free(vars.s);
    lb_puts(b, " => ");
    size_t pushed = push_binders(c, a[0]);
    put(c, a[1], P_BINDER, b);
    c->n_bound -= pushed;
    return P_BINDER;
}

/* ------------------------------------------------------------------ */
/* Heads                                                               */
/* ------------------------------------------------------------------ */

/* Elementary transcendental heads and their Mathlib names. */
static const char* transcendental_name(const char* h, size_t n) {
    if (n == 1) {
        if (h == SYM_Log) return "Real.log";
        if (h == SYM_Exp) return "Real.exp";
        if (h == SYM_Sin) return "Real.sin";
        if (h == SYM_Cos) return "Real.cos";
        if (h == SYM_Tan) return "Real.tan";
        if (h == SYM_ArcTan) return "Real.arctan";
        if (h == SYM_Sinh) return "Real.sinh";
        if (h == SYM_Cosh) return "Real.cosh";
        if (h == SYM_Sqrt) return "Real.sqrt";
    }
    if (n == 2 && h == SYM_Log) return "Real.logb";
    return NULL;
}

/* The pass-through Lean name of a strict-mode head: its display name with
 * context marks as dots (Mathilda`matMul -> Mathilda.matMul), when listed. */
static const char* pass_head(const LCtx* c, const char* h, char* buf, size_t bufsz) {
    const char* nm = context_display_name(h);
    size_t L = strlen(nm);
    if (L + 1 > bufsz) return NULL;
    for (size_t i = 0; i <= L; i++) buf[i] = nm[i] == '`' ? '.' : nm[i];
    for (size_t i = 0; i < c->o.n_pass_heads; i++)
        if (c->o.pass_heads[i] && strcmp(c->o.pass_heads[i], buf) == 0) return buf;
    return NULL;
}

/* head arg1 arg2 ... ; `name` overrides the printed head. */
static int build_application(LCtx* c, const Expr* e, const char* name, LBuf* b) {
    size_t n = e->data.function.arg_count;
    if (name) lb_puts(b, name);
    else put(c, e->data.function.head, P_APP + 1, b);
    for (size_t i = 0; i < n; i++) {
        lb_puts(b, " ");
        put(c, e->data.function.args[i], P_ATOM, b);
    }
    return n ? P_APP : P_ATOM;
}

static int build_list(LCtx* c, const Expr* e, LBuf* b) {
    lb_puts(b, "[");
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (i) lb_puts(b, ", ");
        put(c, e->data.function.args[i], P_BINDER, b);
    }
    lb_puts(b, "]");
    return P_ATOM;
}

static int build_function_node(LCtx* c, const Expr* e, LBuf* b) {
    const char* h = head_name(e);
    size_t n = e->data.function.arg_count;
    Expr* const* a = e->data.function.args;

    if (h) {
        if ((h == SYM_HoldForm || h == SYM_LeanForm) && n == 1) return build(c, a[0], b);
        if (is_number(e)) {
            if (head_is(e, SYM_Complex, 2)) {
                if (c->strict && c->o.core_only) return refuse(c, e, "complex number needs Mathlib (Tier 1)", b);
                if (c->strict && (!is_int(a[0]) && !is_rational(a[0])))
                    return refuse(c, e, "inexact number; Rationalize it first", b);
                lb_puts(b, "(");
                put(c, a[0], P_PLUS, b);
                lb_puts(b, " + ");
                put(c, a[1], P_TIMES, b);
                lb_puts(b, " * Complex.I)");
                return P_ATOM;
            }
            return build_number(c, e, b);
        }
        if (h == SYM_Plus && n >= 2) return build_plus(c, e, b);
        if (h == SYM_Times && n >= 1) return build_times(c, e, b);
        if (h == SYM_Power && n == 2) return build_power(c, e, b);
        if (h == SYM_List) return build_list(c, e, b);
        const char* op = rel_op(h);
        if (op && n >= 2) return build_relation(c, e, op, b);
        if (h == SYM_Inequality && n >= 3 && n % 2 == 1) return build_inequality(c, e, b);
        if (h == SYM_And && n >= 2) return build_nary(c, e, " \xE2\x88\xA7 ", P_AND, b);
        if (h == SYM_Or && n >= 2) return build_nary(c, e, " \xE2\x88\xA8 ", P_OR, b);
        if (h == SYM_Not && n == 1) {
            lb_puts(b, "\xC2\xAC");                     /* U+00AC not */
            put(c, a[0], P_NOT, b);
            return P_NOT;
        }
        if (h == SYM_Implies && n == 2) {
            put(c, a[0], P_IMP + 1, b);
            lb_puts(b, " \xE2\x86\x92 ");               /* U+2192 -> */
            put(c, a[1], P_IMP, b);
            return P_IMP;
        }
        if (h == SYM_Equivalent && n >= 2) {
            if (n == 2) {
                put(c, a[0], P_IFF + 1, b);
                lb_puts(b, " \xE2\x86\x94 ");           /* U+2194 <-> */
                put(c, a[1], P_IFF + 1, b);
                return P_IFF;
            }
            for (size_t i = 0; i + 1 < n; i++) {
                if (i) lb_puts(b, " \xE2\x88\xA7 ");
                lb_puts(b, "(");
                put(c, a[i], P_IFF + 1, b);
                lb_puts(b, " \xE2\x86\x94 ");
                put(c, a[i + 1], P_IFF + 1, b);
                lb_puts(b, ")");
            }
            return P_AND;
        }
        if ((h == SYM_ForAll || h == SYM_Exists) && (n == 2 || n == 3))
            return build_quantifier(c, e, h == SYM_ForAll, b);
        if (h == SYM_Function) {
            LBuf t = {0};
            int p = build_function(c, e, &t);
            if (p >= 0) { if (t.s) lb_putn(b, t.s, t.n); free(t.s); return p; }
            free(t.s);
        }
        if (h == SYM_Slot && c->strict) return refuse(c, e, "pure function with Slot", b);
        const char* tn = transcendental_name(h, n);
        if (tn) {
            if (c->strict && c->o.core_only)
                return refuse(c, e, "transcendental function needs Mathlib (Tier 1)", b);
            return build_application(c, e, tn, b);
        }
        if (c->strict) {
            char buf[256];
            const char* pn = pass_head(c, h, buf, sizeof buf);
            if (!pn) return refuse(c, e, "no Lean counterpart for this head", b);
            return build_application(c, e, pn, b);
        }
    } else if (c->strict) {
        return refuse(c, e, "compound head", b);
    }
    return build_application(c, e, NULL, b);
}

static int build(LCtx* c, const Expr* e, LBuf* b) {
    if (!e) { lb_puts(b, "Null"); return P_ATOM; }

    /* Format[f[...], LeanForm] := rhs (display form only). */
    if (!c->strict && (e->type == EXPR_FUNCTION || e->type == EXPR_SYMBOL) &&
        c->fmt_depth < FMT_DEPTH_LIMIT) {
        Expr* shown = print_form_format_rewrite((Expr*)e, SYM_LeanForm);
        if (shown) {
            c->fmt_depth++;
            int p = build(c, shown, b);
            c->fmt_depth--;
            expr_free(shown);
            return p;
        }
    }

    switch (e->type) {
    case EXPR_INTEGER:
    case EXPR_BIGINT:
    case EXPR_REAL:
#ifdef USE_MPFR
    case EXPR_MPFR:
#endif
        return build_number(c, e, b);
    case EXPR_STRING:
        put_string_literal(e->data.string, b);
        return P_ATOM;
    case EXPR_SYMBOL:
        return build_symbol(c, e, b);
    case EXPR_FUNCTION:
        return build_function_node(c, e, b);
    case EXPR_NDARRAY: {
        Expr* l = ndarray_to_nested_list((Expr*)e);
        int p = l ? build(c, l, b) : refuse(c, e, "unsupported array", b);
        if (l) expr_free(l);
        return p;
    }
    default: {
        if (c->strict) return refuse(c, e, "unsupported atom", b);
        char* s = expr_to_string_outputform((Expr*)e);
        lb_puts(b, s ? s : "");
        free(s);
        return P_ATOM;
    }
    }
}

/* ------------------------------------------------------------------ */
/* Entry points                                                        */
/* ------------------------------------------------------------------ */

char* lean_form_string(const Expr* e) {
    LCtx c;
    memset(&c, 0, sizeof c);
    LBuf b = {0};
    build(&c, e, &b);
    ctx_release(&c);
    return lb_take(&b);
}

char* lean_form_string_ex(const Expr* e, const LeanFormOpts* o, const char** err,
                          const Expr** bad) {
    LCtx c;
    memset(&c, 0, sizeof c);
    c.strict = true;
    if (o) c.o = *o;
    if (c.o.num == LEAN_NUM_INT) c.o.allow_division = 0;
    LBuf b = {0};
    build(&c, e, &b);
    ctx_release(&c);
    if (err) *err = c.err;
    if (bad) *bad = c.bad;
    if (c.err) { free(b.s); return NULL; }
    return lb_take(&b);
}

Expr* builtin_lean_tolean(Expr* res) {
    if (res->type != EXPR_FUNCTION ||
        (res->data.function.arg_count != 4 && res->data.function.arg_count != 5))
        return NULL;
    Expr* const* a = res->data.function.args;
    LeanFormOpts o;
    memset(&o, 0, sizeof o);
    o.allow_division = 1;
    if (res->data.function.arg_count == 5) {
        if (a[4]->type != EXPR_SYMBOL ||
            (a[4]->data.symbol.name != SYM_True && a[4]->data.symbol.name != SYM_False))
            return NULL;
        o.closed = a[4]->data.symbol.name == SYM_True;
    }
    if (a[1]->type == EXPR_STRING) {
        const char* t = a[1]->data.string;
        if (strcmp(t, "Int") == 0) o.num = LEAN_NUM_INT;
        else if (strcmp(t, "Rat") == 0) o.num = LEAN_NUM_RAT;
        else if (strcmp(t, "Real") == 0) o.num = LEAN_NUM_REAL;
        else return NULL;
    } else if (!(a[1]->type == EXPR_SYMBOL && a[1]->data.symbol.name == SYM_None)) {
        return NULL;
    }
    if (a[2]->type != EXPR_SYMBOL ||
        (a[2]->data.symbol.name != SYM_True && a[2]->data.symbol.name != SYM_False))
        return NULL;
    o.core_only = a[2]->data.symbol.name == SYM_True;
    if (head_name(a[3]) != SYM_List) return NULL;
    size_t np = a[3]->data.function.arg_count;
    const char** pass = np ? malloc(np * sizeof *pass) : NULL;
    if (np && !pass) abort();
    for (size_t i = 0; i < np; i++) {
        const Expr* s = a[3]->data.function.args[i];
        if (s->type != EXPR_STRING) { free((void*)pass); return NULL; }
        pass[i] = s->data.string;
    }
    o.pass_heads = pass;
    o.n_pass_heads = np;

    const char* err = NULL;
    const Expr* bad = NULL;
    char* s = lean_form_string_ex(a[0], &o, &err, &bad);
    free((void*)pass);
    if (!s) {
        const Expr* shown = bad;
        if (shown && head_is(shown, SYM_HoldForm, 1)) shown = shown->data.function.args[0];
        char* txt = shown ? expr_to_string_inputform_display((Expr*)shown) : NULL;
        mth_message("LeanForm", "untrans", "cannot translate %s to Lean (%s).",
                    txt ? txt : "expression", err ? err : "unsupported");
        free(txt);
        return expr_new_symbol(SYM_DollarFailed);
    }
    return expr_new_string_adopt(s);
}
