/*
 * ratfun.c -- factored rational functions over Q: Together, Cancel, Apart.
 *
 * See ratfun.h for the contract.  The engine keeps a rational function as
 *
 *     c * prod_i p_i^e_i          (Rf: c in Q, e_i != 0 integers)
 *
 * where every p_i is a non-constant fmpz_mpoly, PRIMITIVE (integer content 1)
 * and NORMALISED (its leading coefficient in lex order is positive).  The
 * variables are the symbols and transcendental kernels of the input in the
 * evaluator's canonical order; variable 0 is the most significant, so "leading
 * coefficient positive" is the sign convention of the printed answer
 * ((a - b), (x - y), (-1 + x)).  Apart puts its variable first instead, so the
 * factors it divides by have a positive leading coefficient in that variable.
 *
 * Two regimes produce the answer:
 *
 *  - the SUM regime, for a Plus with a denominator anywhere in it (Together's
 *    combination step).  Each term is first cancelled on its own; the
 *    denominator factors of all terms are split into their monomial part and
 *    their squarefree parts, and a gcd-free basis B_1..B_s of those is built.
 *    The common denominator is L = prod B_k^M_k (M_k the largest multiplicity
 *    of B_k in any term), the numerator N = sum_t c_t num_t (L / den_t) is
 *    expanded, and N is cancelled against the squarefree product prod B_k --
 *    one gcd in the usual case that nothing cancels.
 *
 *  - the SINGLE-TERM regime, for a product (the top level of Together on a
 *    non-sum, and every term of Cancel).  The factors are kept as they are; the
 *    denominator factors are split into their monomial content and made
 *    pairwise coprime, and every numerator factor is split against them by gcd
 *    until nothing is shared:  (x^2-1)^2/((x-1)(x+5)) -> (x-1)(x+1)^2/(x+5).
 *
 * The answer is rendered as Times[c, factors...] and evaluated once.  When c is
 * exactly -1 the sign is folded into the single non-monomial numerator factor,
 * or failing that into the single non-monomial denominator factor, of exponent
 * +-1:  1/(1-x) stays (1-x)^-1 rather than -(1/(-1+x)).
 *
 * Kernels: Log[u], the inverse trig/hyperbolic functions and inert user-head
 * atoms (x[3]) are independent generators; the exponentials E^(k u) for integer
 * k become the powers g^k of ONE generator g = E^u, so E^(2x) - 1 factors as
 * (E^x - 1)(E^x + 1).  Anything else (I, Sqrt[2], Sin[x], a^x, Reals) declines.
 */
#include "ratfun.h"

#ifdef USE_FLINT

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_mpoly_factor.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

#include "eval.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "flint_bridge_internal.h"

/* gr_poly_xgcd arrived in FLINT 3.2; the Euclidean variant has the same
 * signature and is older (see flint_bridge.c). */
#if defined(__FLINT_RELEASE) && __FLINT_RELEASE >= 30200
#define RF_GR_POLY_XGCD gr_poly_xgcd
#else
#define RF_GR_POLY_XGCD gr_poly_xgcd_euclidean
#endif

/* ------------------------------------------------------------------ */
/*  Variables: the atoms of the input in canonical order               */
/* ------------------------------------------------------------------ */

typedef struct {
    const Expr** atom;      /* borrowed variable atoms (atom[gidx] is gexp)  */
    size_t n, cap;
    const Expr** xk;        /* distinct exponential kernels E^e (borrowed)   */
    int64_t* xpow;          /* E^e = g^xpow                                  */
    size_t xn, xcap;
    Expr* gexp;             /* owned generator E^u, or NULL                  */
    slong gidx;             /* variable index of g, -1 when none             */
    const fmpz_mpoly_ctx_struct* ctx;   /* cached; never cleared here         */
    int fail;               /* a gcd failed or a division by zero was met    */
    long splits;            /* factors split by a common gcd (rf_single)     */
} RfCtx;

static int rf_push_ptr(const Expr*** v, size_t* n, size_t* cap, const Expr* e) {
    if (*n == *cap) {
        size_t nc = *cap ? *cap * 2 : 8;
        const Expr** nv = realloc((void*)*v, nc * sizeof(*nv));
        if (!nv) return 0;
        *v = nv; *cap = nc;
    }
    (*v)[(*n)++] = e;
    return 1;
}

/* Symbol names are interned, so a pointer compare decides most lookups. */
static int rf_atom_eq(const Expr* a, const Expr* b) {
    if (a == b) return 1;
    if (a->type == EXPR_SYMBOL || b->type == EXPR_SYMBOL)
        return a->type == b->type && (a->data.symbol.name == b->data.symbol.name ||
                                      strcmp(a->data.symbol.name, b->data.symbol.name) == 0);
    return expr_eq(a, b);
}

static slong rf_atom_index(const RfCtx* C, const Expr* e) {
    for (size_t i = 0; i < C->n; i++)
        if ((slong)i != C->gidx && rf_atom_eq(C->atom[i], e)) return (slong)i;
    return -1;
}

static int rf_head_is(const Expr* e, const char* sym) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}

/* Symbols that are not variables of a rational function. */
static int rf_bad_symbol(const Expr* e) {
    const char* s = e->data.symbol.name;
    return strcmp(s, "Infinity") == 0 || strcmp(s, "ComplexInfinity") == 0 ||
           strcmp(s, "Indeterminate") == 0 || strcmp(s, "I") == 0 ||
           strcmp(s, "Null") == 0 || strcmp(s, "True") == 0 || strcmp(s, "False") == 0 ||
           strncmp(s, "$flk", 4) == 0 || strncmp(s, "~flk", 4) == 0;
}

static int rf_add_atom(RfCtx* C, const Expr* e) {
    for (size_t i = 0; i < C->n; i++) if (rf_atom_eq(C->atom[i], e)) return 1;
    return rf_push_ptr(&C->atom, &C->n, &C->cap, e);
}

/* Pass 1: collect the atoms; 0 declines. */
static int rf_scan(const Expr* e, RfCtx* C) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT:
            return 1;
        case EXPR_SYMBOL:
            return rf_bad_symbol(e) ? 0 : rf_add_atom(C, e);
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL) return 0;
            const char* hn = h->data.symbol.name;
            size_t n = e->data.function.arg_count;
            if (hn == SYM_Rational) return 1;
            if (hn == SYM_Plus || hn == SYM_Times) {
                for (size_t i = 0; i < n; i++)
                    if (!rf_scan(e->data.function.args[i], C)) return 0;
                return 1;
            }
            if (hn == SYM_Power) {
                if (n != 2) return 0;
                const Expr* b = e->data.function.args[0];
                const Expr* x = e->data.function.args[1];
                if (x->type == EXPR_INTEGER) return rf_scan(b, C);
                if (b->type == EXPR_SYMBOL && b->data.symbol.name == SYM_E) {
                    for (size_t i = 0; i < C->xn; i++) if (expr_eq(C->xk[i], e)) return 1;
                    return rf_push_ptr(&C->xk, &C->xn, &C->xcap, e);
                }
                return 0;
            }
            if (fb_km_is_generic_kernel_head(hn) || fb_km_is_inert_head_atom(e))
                return rf_add_atom(C, e);
            return 0;
        }
        default:
            return 0;
    }
}

/* evaluate(a/b) when it is a nonzero machine integer (an exponent ratio). */
static int rf_int_ratio(const Expr* a, const Expr* b, int64_t* k) {
    Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_copy((Expr*)a),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy((Expr*)b), expr_new_integer(-1) }, 2) }, 2));
    int ok = q && q->type == EXPR_INTEGER && q->data.integer != 0;
    if (ok) *k = q->data.integer;
    if (q) expr_free(q);
    return ok;
}

/* A syntactically negative quantity: a negative number, a product with a
 * negative numeric factor, or a sum whose first term is negative. */
static int rf_expr_negative(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) < 0;
    if (rf_head_is(e, SYM_Rational) && e->data.function.arg_count == 2)
        return rf_expr_negative(e->data.function.args[0]);
    if ((rf_head_is(e, SYM_Times) || rf_head_is(e, SYM_Plus)) && e->data.function.arg_count > 0)
        return rf_expr_negative(e->data.function.args[0]);
    return 0;
}

static int rf_cmp_atoms(const void* a, const void* b) {
    return expr_compare(*(const Expr* const*)a, *(const Expr* const*)b);
}

/* One ORD_LEX context per variable count, built on first use and kept for the
 * process lifetime (a context is a small read-only descriptor). */
#define RF_CTX_CAP 64
static fmpz_mpoly_ctx_struct* g_rf_ctx[RF_CTX_CAP + 1];

static const fmpz_mpoly_ctx_struct* rf_ctx_for(slong nv) {
    if (nv < 1 || nv > RF_CTX_CAP) return NULL;
    if (!g_rf_ctx[nv]) {
        fmpz_mpoly_ctx_struct* c = malloc(sizeof *c);
        if (!c) return NULL;
        fmpz_mpoly_ctx_init(c, nv, ORD_LEX);
        g_rf_ctx[nv] = c;
    }
    return g_rf_ctx[nv];
}

static void rf_ctx_clear(RfCtx* C) {
    free((void*)C->atom); free((void*)C->xk); free(C->xpow);
    if (C->gexp) expr_free(C->gexp);
    memset(C, 0, sizeof *C);
    C->gidx = -1;
}

/* MATHILDA_TOGETHER_FORM=factored renders the factored form (the denominator a
 * product of its gcd-free basis); the default renders the classical form, an
 * expanded numerator over an expanded denominator with the integer content
 * folded into each, which is what the internal callers read. */
static int rf_factored_form(void) {
    static int f = -1;
    if (f < 0) {
        const char* e = getenv("MATHILDA_TOGETHER_FORM");
        f = (e && strcmp(e, "factored") == 0) ? 1 : 0;
    }
    return f;
}

/* Scan `e`, settle the exponential generator, order the variables.  With
 * `first` (Apart's variable) that atom is moved to index 0. */
static int rf_ctx_build(const Expr* e, RfCtx* C, const Expr* first) {
    memset(C, 0, sizeof *C);
    C->gidx = -1;
    if (!rf_scan(e, C)) { rf_ctx_clear(C); return 0; }
    if (C->xn > 0) {
        /* A fundamental exponent u that every exponent is an integer multiple
         * of, taken from the kernels themselves. */
        const Expr* u = NULL;
        C->xpow = malloc(C->xn * sizeof(int64_t));
        if (!C->xpow) { rf_ctx_clear(C); return 0; }
        for (size_t c = 0; c < C->xn && !u; c++) {
            const Expr* cand = C->xk[c]->data.function.args[1];
            int good = 1;
            for (size_t j = 0; j < C->xn && good; j++)
                good = rf_int_ratio(C->xk[j]->data.function.args[1], cand, &C->xpow[j]);
            if (good) u = cand;
        }
        if (!u) { rf_ctx_clear(C); return 0; }
        /* Orient the generator positively, g = E^x rather than E^(-x), so a
         * 1/E^x is a denominator: Together[1 + E^(-x)] -> (1 + E^x)/E^x. */
        Expr* uu = expr_copy((Expr*)u);
        if (rf_expr_negative(uu)) {
            expr_free(uu);
            uu = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)u) }, 2));
            for (size_t j = 0; j < C->xn; j++) C->xpow[j] = -C->xpow[j];
        }
        C->gexp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_new_symbol(SYM_E), uu }, 2));
        if (!rf_push_ptr(&C->atom, &C->n, &C->cap, C->gexp)) { rf_ctx_clear(C); return 0; }
    }
    if (C->n == 0) { rf_ctx_clear(C); return 0; }
    qsort((void*)C->atom, C->n, sizeof(*C->atom), rf_cmp_atoms);
    if (first) {
        size_t k = C->n;
        for (size_t i = 0; i < C->n; i++) if (rf_atom_eq(C->atom[i], first)) { k = i; break; }
        if (k == C->n) { rf_ctx_clear(C); return 0; }
        const Expr* a = C->atom[k];
        memmove((void*)(C->atom + 1), (void*)C->atom, k * sizeof(*C->atom));
        C->atom[0] = a;
    }
    for (size_t i = 0; i < C->n && C->gexp; i++)
        if (C->atom[i] == C->gexp) C->gidx = (slong)i;
    C->ctx = rf_ctx_for((slong)C->n);
    if (!C->ctx) { rf_ctx_clear(C); return 0; }
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Rf: c * prod p_i^e_i                                               */
/* ------------------------------------------------------------------ */

typedef struct { fmpz_mpoly_struct p; slong e; } RfFac;
/* `reduced`: known to be in lowest terms (set by rf_single and rf_sum, kept
 * by powers and by multiplication with a constant), so a sum need not cancel
 * the term again. */
typedef struct { fmpq_t c; RfFac* f; slong n, cap; int reduced; } Rf;

static void rf_init(Rf* R) { fmpq_init(R->c); fmpq_one(R->c); R->f = NULL; R->n = R->cap = 0; R->reduced = 1; }

static void rf_clear(Rf* R, const RfCtx* C) {
    for (slong i = 0; i < R->n; i++) fmpz_mpoly_clear(&R->f[i].p, C->ctx);
    free(R->f);
    fmpq_clear(R->c);
    R->f = NULL; R->n = R->cap = 0;
}

/* Remove factor i (the last factor moves into slot i). */
static void rf_remove(Rf* R, slong i, const RfCtx* C) {
    fmpz_mpoly_clear(&R->f[i].p, C->ctx);
    R->f[i] = R->f[R->n - 1];
    R->n--;
}

/* c *= z^e for an integer z; 0 on 0^negative. */
static int rf_mul_fmpz_pow(fmpq_t c, const fmpz_t z, slong e) {
    if (e == 0) return 1;
    if (e < 0 && fmpz_is_zero(z)) return 0;
    fmpz_t p; fmpz_init(p);
    fmpz_pow_ui(p, z, (ulong)(e < 0 ? -e : e));
    if (e > 0) fmpq_mul_fmpz(c, c, p);
    else       fmpq_div_fmpz(c, c, p);
    fmpz_clear(p);
    return 1;
}

/* Make p primitive with a positive leading coefficient; return the unit-and-
 * content it was divided by in `g` (g = 0 for p = 0). */
static void rf_normalise(fmpz_mpoly_t p, fmpz_t g, const fmpz_mpoly_ctx_struct* ctx) {
    if (p->length == 0) { fmpz_zero(g); return; }
    _fmpz_vec_content(g, p->coeffs, p->length);
    if (fmpz_sgn(p->coeffs + 0) < 0) fmpz_neg(g, g);
    if (!fmpz_is_one(g)) fmpz_mpoly_scalar_divexact_fmpz(p, p, g, ctx);
}

static int rf_push_gen(Rf* R, slong v, slong e, RfCtx* C);

/* Push p^e; p is consumed (left zero).  p is normalised here, its content and
 * sign going into c; a constant p multiplies c; equal factors merge. */
static int rf_push(Rf* R, fmpz_mpoly_t p, slong e, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    if (e == 0) { fmpz_mpoly_zero(p, ctx); return 1; }
    fmpz_t g; fmpz_init(g);
    if (fmpz_mpoly_is_fmpz(p, ctx)) {
        fmpz_mpoly_get_fmpz(g, p, ctx);
        int ok = rf_mul_fmpz_pow(R->c, g, e);
        fmpz_clear(g);
        fmpz_mpoly_zero(p, ctx);
        if (!ok) C->fail = 1;
        return ok;
    }
    rf_normalise(p, g, ctx);
    if (!fmpz_is_one(g)) rf_mul_fmpz_pow(R->c, g, e);
    fmpz_clear(g);
    if (p->length == 1 && fmpz_mpoly_total_degree_si(p, ctx) > 1) {
        /* a monomial x^a y^b becomes the generator factors x^(a e), y^(b e):
         * every monomial factor is ONE generator (rf_is_mono relies on it) */
        slong nv = (slong)C->n;
        ulong* ex = malloc((size_t)nv * sizeof(ulong));
        if (!ex) { C->fail = 1; return 0; }
        fmpz_mpoly_get_term_exp_ui(ex, p, 0, ctx);
        fmpz_mpoly_zero(p, ctx);
        int ok = 1;
        for (slong v = 0; v < nv && ok; v++)
            if (ex[v]) ok = rf_push_gen(R, v, e * (slong)ex[v], C);
        free(ex);
        return ok;
    }
    for (slong i = 0; i < R->n; i++) {
        if (fmpz_mpoly_equal(&R->f[i].p, p, ctx)) {
            R->f[i].e += e;
            if (R->f[i].e == 0) rf_remove(R, i, C);
            fmpz_mpoly_zero(p, ctx);
            return 1;
        }
    }
    if (R->n == R->cap) {
        slong nc = R->cap ? 2 * R->cap : 4;
        RfFac* nf = realloc(R->f, (size_t)nc * sizeof(RfFac));
        if (!nf) { C->fail = 1; return 0; }
        R->f = nf; R->cap = nc;
    }
    fmpz_mpoly_init(&R->f[R->n].p, ctx);
    fmpz_mpoly_swap(&R->f[R->n].p, p, ctx);
    R->f[R->n].e = e;
    R->n++;
    return 1;
}

/* R *= S (S consumed and reset to 1). */
static int rf_mul_into(Rf* R, Rf* S, RfCtx* C) {
    int ok = 1;
    if (R->n == 0) R->reduced = S->reduced;
    else if (S->n != 0) R->reduced = 0;
    fmpq_mul(R->c, R->c, S->c);
    for (slong i = 0; i < S->n && ok; i++) ok = rf_push(R, &S->f[i].p, S->f[i].e, C);
    rf_clear(S, C);
    rf_init(S);
    return ok;
}

static int rf_pow(Rf* R, slong k, RfCtx* C) {
    if (k < 0) {
        if (fmpq_is_zero(R->c)) { C->fail = 1; return 0; }
        fmpq_inv(R->c, R->c);
    }
    fmpq_pow_si(R->c, R->c, k < 0 ? -k : k);
    for (slong i = 0; i < R->n; i++) R->f[i].e *= k;
    return 1;
}

static int rf_is_poly(const Rf* R) {
    for (slong i = 0; i < R->n; i++) if (R->f[i].e < 0) return 0;
    return 1;
}

/* A single-term factor; normalised, it is one variable (rf_push_gen never
 * pushes a multi-variable monomial). */
static int rf_is_mono(const fmpz_mpoly_struct* p) { return p->length == 1; }

/* prod p_i^e_i over one side: the factors with e > 0 (or e < 0 when `den`). */
static void rf_expand_side(fmpz_mpoly_t out, const Rf* R, int den, const RfCtx* C) {
    fmpz_mpoly_one(out, C->ctx);
    fmpz_mpoly_t t; fmpz_mpoly_init(t, C->ctx);
    for (slong i = 0; i < R->n; i++) {
        slong e = den ? -R->f[i].e : R->f[i].e;
        if (e <= 0) continue;
        if (e == 1) fmpz_mpoly_mul(out, out, &R->f[i].p, C->ctx);
        else {
            fmpz_mpoly_pow_ui(t, &R->f[i].p, (ulong)e, C->ctx);
            fmpz_mpoly_mul(out, out, t, C->ctx);
        }
    }
    fmpz_mpoly_clear(t, C->ctx);
}

static int rf_push_gen(Rf* R, slong v, slong e, RfCtx* C) {
    fmpz_mpoly_t g; fmpz_mpoly_init(g, C->ctx);
    fmpz_mpoly_gen(g, v, C->ctx);
    int ok = rf_push(R, g, e, C);
    fmpz_mpoly_clear(g, C->ctx);
    return ok;
}

/* Split each denominator factor into its monomial content and the rest. */
static int rf_split_den_monomials(Rf* R, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    slong nv = (slong)C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    if (!ex) return 0;
    fmpz_mpoly_t m, q;
    fmpz_mpoly_init(m, ctx); fmpz_mpoly_init(q, ctx);
    int ok = 1;
    for (slong i = 0; i < R->n && ok; i++) {
        if (R->f[i].e >= 0 || rf_is_mono(&R->f[i].p)) continue;
        fmpz_mpoly_term_content(m, &R->f[i].p, ctx);
        if (fmpz_mpoly_is_fmpz(m, ctx)) continue;
        fmpz_mpoly_get_term_exp_ui(ex, m, 0, ctx);
        fmpz_one(m->coeffs + 0);
        if (!fmpz_mpoly_divides(q, &R->f[i].p, m, ctx)) { C->fail = 1; ok = 0; break; }
        slong e = R->f[i].e;
        /* re-push the quotient so it merges with an equal factor: x^3 (2x-3)
         * from (3x^3 - 2x^4) next to a (3 - 2x) of its own */
        rf_remove(R, i, C);
        ok = rf_push(R, q, e, C);
        for (slong v = 0; v < nv && ok; v++)
            if (ex[v]) ok = rf_push_gen(R, v, e * (slong)ex[v], C);
        i = -1;                                    /* the list was reordered */
    }
    fmpz_mpoly_clear(m, ctx); fmpz_mpoly_clear(q, ctx);
    free(ex);
    return ok;
}

/* 1 when gcd(A, B) is non-constant (set into G); a FLINT failure sets fail. */
static int rf_gcd_nontrivial(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                             RfCtx* C) {
    if (rf_is_mono(A) && rf_is_mono(B) && !fmpz_mpoly_equal(A, B, C->ctx)) return 0;
    if (!fmpz_mpoly_gcd(G, A, B, C->ctx)) { C->fail = 1; return 0; }
    return !fmpz_mpoly_is_fmpz(G, C->ctx);
}

/* Replace factor i (p^e) by g^e (p/g)^e. */
static int rf_split_factor(Rf* R, slong i, const fmpz_mpoly_t g, RfCtx* C) {
    C->splits++;
    fmpz_mpoly_t q, gg;
    fmpz_mpoly_init(q, C->ctx); fmpz_mpoly_init(gg, C->ctx);
    int ok = fmpz_mpoly_divides(q, &R->f[i].p, g, C->ctx);
    if (ok) {
        slong e = R->f[i].e;
        fmpz_mpoly_set(gg, g, C->ctx);
        rf_remove(R, i, C);
        ok = rf_push(R, q, e, C) && rf_push(R, gg, e, C);
    } else C->fail = 1;
    fmpz_mpoly_clear(q, C->ctx); fmpz_mpoly_clear(gg, C->ctx);
    return ok;
}

/* Single-term regime (see the file comment).  Each pass finds one pair of
 * factors -- a denominator factor against a numerator factor or a later
 * denominator factor -- with a non-constant gcd g, and splits both by g;
 * rf_push's merge then cancels or combines the copies of g.  Every split
 * lowers a degree, so the loop ends; inputs have few factors. */
static int rf_single(Rf* R, RfCtx* C) {
    if (!rf_split_den_monomials(R, C)) return 0;
    fmpz_mpoly_t g; fmpz_mpoly_init(g, C->ctx);
    fmpz_t u; fmpz_init(u);
    int found = 1, ok = 1;
    while (found && ok && !C->fail) {
        found = 0;
        for (slong i = 0; i < R->n && !found; i++) {
            if (R->f[i].e >= 0) continue;
            for (slong j = 0; j < R->n && !found; j++) {
                if (j == i || (R->f[j].e < 0 && j < i)) continue;
                if (fmpz_mpoly_equal(&R->f[i].p, &R->f[j].p, C->ctx)) {
                    C->splits++;
                    /* the same factor twice: merge (rf_push merges on entry,
                     * this only guards the invariant) */
                    slong a = i > j ? i : j, b = i > j ? j : i;
                    R->f[b].e += R->f[a].e;
                    rf_remove(R, a, C);
                    if (R->f[b].e == 0) rf_remove(R, b, C);
                    found = 1;
                    continue;
                }
                if (!rf_gcd_nontrivial(g, &R->f[i].p, &R->f[j].p, C)) continue;
                rf_normalise(g, u, C->ctx);
                /* split the higher slot first: rf_remove moves the last
                 * factor into the vacated slot, which leaves lower slots alone */
                slong a = i > j ? i : j, b = i > j ? j : i;
                int same_a = fmpz_mpoly_equal(&R->f[a].p, g, C->ctx);
                int same_b = fmpz_mpoly_equal(&R->f[b].p, g, C->ctx);
                if (!same_a) ok = rf_split_factor(R, a, g, C);
                if (ok && !same_b) ok = rf_split_factor(R, b, g, C);
                found = 1;
            }
        }
    }
    fmpz_clear(u);
    fmpz_mpoly_clear(g, C->ctx);
    if (ok && !C->fail) R->reduced = 1;
    return ok && !C->fail;
}

/* ------------------------------------------------------------------ */
/*  The sum regime                                                     */
/* ------------------------------------------------------------------ */

/* The combination step runs over a small polynomial-ring interface so the one
 * algorithm serves both representations: sparse fmpz_mpoly, and dense fmpz_poly
 * when the input has a single variable (several times faster there: dense
 * multiplication, division and the modular gcd instead of heap arithmetic). */
typedef struct RfRing RfRing;
struct RfRing {
    size_t sz;                                    /* bytes per element */
    const fmpz_mpoly_ctx_struct* ctx;             /* the mpoly context */
    void (*init)(void* a, const RfRing* R);
    void (*clear)(void* a, const RfRing* R);
    void (*set)(void* a, const void* b, const RfRing* R);
    void (*one)(void* a, const RfRing* R);
    void (*swap)(void* a, void* b, const RfRing* R);
    int  (*is_const)(const void* a, const RfRing* R);
    int  (*is_zero)(const void* a, const RfRing* R);
    int  (*is_linear)(const void* a, const RfRing* R);   /* total degree 1 */
    int  (*equal)(const void* a, const void* b, const RfRing* R);
    void (*mul)(void* a, const void* b, const void* c, const RfRing* R);
    void (*add)(void* a, const void* b, const void* c, const RfRing* R);
    void (*pow)(void* a, const void* b, ulong e, const RfRing* R);
    void (*smul)(void* a, const void* b, const fmpz_t c, const RfRing* R);
    int  (*gcd)(void* g, const void* a, const void* b, const RfRing* R);
    int  (*divides)(void* q, const void* a, const void* b, const RfRing* R);
    void (*normalise)(void* a, const RfRing* R);
    void (*from_mpoly)(void* a, const fmpz_mpoly_t b, const RfRing* R);
    void (*to_mpoly)(fmpz_mpoly_t a, const void* b, const RfRing* R);
};

/* -- fmpz_mpoly -- */
#define RMP(a) ((fmpz_mpoly_struct*)(a))
#define RMPc(a) ((const fmpz_mpoly_struct*)(a))
static void rm_init(void* a, const RfRing* R) { fmpz_mpoly_init(RMP(a), R->ctx); }
static void rm_clear(void* a, const RfRing* R) { fmpz_mpoly_clear(RMP(a), R->ctx); }
static void rm_set(void* a, const void* b, const RfRing* R) { fmpz_mpoly_set(RMP(a), RMPc(b), R->ctx); }
static void rm_one(void* a, const RfRing* R) { fmpz_mpoly_one(RMP(a), R->ctx); }
static void rm_swap(void* a, void* b, const RfRing* R) { fmpz_mpoly_swap(RMP(a), RMP(b), R->ctx); }
static int rm_is_const(const void* a, const RfRing* R) { return fmpz_mpoly_is_fmpz(RMPc(a), R->ctx); }
static int rm_is_zero(const void* a, const RfRing* R) { return fmpz_mpoly_is_zero(RMPc(a), R->ctx); }
static int rm_is_linear(const void* a, const RfRing* R) { return fmpz_mpoly_total_degree_si(RMPc(a), R->ctx) == 1; }
static int rm_equal(const void* a, const void* b, const RfRing* R) { return fmpz_mpoly_equal(RMPc(a), RMPc(b), R->ctx); }
static void rm_mul(void* a, const void* b, const void* c, const RfRing* R) { fmpz_mpoly_mul(RMP(a), RMPc(b), RMPc(c), R->ctx); }
static void rm_add(void* a, const void* b, const void* c, const RfRing* R) { fmpz_mpoly_add(RMP(a), RMPc(b), RMPc(c), R->ctx); }
static void rm_pow(void* a, const void* b, ulong e, const RfRing* R) { fmpz_mpoly_pow_ui(RMP(a), RMPc(b), e, R->ctx); }
static void rm_smul(void* a, const void* b, const fmpz_t c, const RfRing* R) { fmpz_mpoly_scalar_mul_fmpz(RMP(a), RMPc(b), c, R->ctx); }
static int rm_gcd(void* g, const void* a, const void* b, const RfRing* R) {
    if (RMPc(a)->length == 1 && RMPc(b)->length == 1 && !fmpz_mpoly_equal(RMPc(a), RMPc(b), R->ctx)) {
        fmpz_mpoly_one(RMP(g), R->ctx);       /* distinct generators */
        return 1;
    }
    return fmpz_mpoly_gcd(RMP(g), RMPc(a), RMPc(b), R->ctx);
}
static int rm_divides(void* q, const void* a, const void* b, const RfRing* R) { return fmpz_mpoly_divides(RMP(q), RMPc(a), RMPc(b), R->ctx); }
static void rm_normalise(void* a, const RfRing* R) {
    fmpz_t g; fmpz_init(g);
    rf_normalise(RMP(a), g, R->ctx);
    fmpz_clear(g);
}
static void rm_from(void* a, const fmpz_mpoly_t b, const RfRing* R) { fmpz_mpoly_set(RMP(a), b, R->ctx); }
static void rm_to(fmpz_mpoly_t a, const void* b, const RfRing* R) { fmpz_mpoly_set(a, RMPc(b), R->ctx); }

/* -- fmpz_poly (one variable) -- */
#define RUP(a) ((fmpz_poly_struct*)(a))
#define RUPc(a) ((const fmpz_poly_struct*)(a))
static void ru_init(void* a, const RfRing* R) { (void)R; fmpz_poly_init(RUP(a)); }
static void ru_clear(void* a, const RfRing* R) { (void)R; fmpz_poly_clear(RUP(a)); }
static void ru_set(void* a, const void* b, const RfRing* R) { (void)R; fmpz_poly_set(RUP(a), RUPc(b)); }
static void ru_one(void* a, const RfRing* R) { (void)R; fmpz_poly_one(RUP(a)); }
static void ru_swap(void* a, void* b, const RfRing* R) { (void)R; fmpz_poly_swap(RUP(a), RUP(b)); }
static int ru_is_const(const void* a, const RfRing* R) { (void)R; return fmpz_poly_length(RUPc(a)) <= 1; }
static int ru_is_zero(const void* a, const RfRing* R) { (void)R; return fmpz_poly_is_zero(RUPc(a)); }
static int ru_is_linear(const void* a, const RfRing* R) { (void)R; return fmpz_poly_length(RUPc(a)) == 2; }
static int ru_equal(const void* a, const void* b, const RfRing* R) { (void)R; return fmpz_poly_equal(RUPc(a), RUPc(b)); }
static void ru_mul(void* a, const void* b, const void* c, const RfRing* R) { (void)R; fmpz_poly_mul(RUP(a), RUPc(b), RUPc(c)); }
static void ru_add(void* a, const void* b, const void* c, const RfRing* R) { (void)R; fmpz_poly_add(RUP(a), RUPc(b), RUPc(c)); }
static void ru_pow(void* a, const void* b, ulong e, const RfRing* R) { (void)R; fmpz_poly_pow(RUP(a), RUPc(b), e); }
static void ru_smul(void* a, const void* b, const fmpz_t c, const RfRing* R) { (void)R; fmpz_poly_scalar_mul_fmpz(RUP(a), RUPc(b), c); }
static int ru_gcd(void* g, const void* a, const void* b, const RfRing* R) { (void)R; fmpz_poly_gcd(RUP(g), RUPc(a), RUPc(b)); return 1; }
static int ru_divides(void* q, const void* a, const void* b, const RfRing* R) { (void)R; return fmpz_poly_divides(RUP(q), RUPc(a), RUPc(b)); }
static void ru_normalise(void* a, const RfRing* R) {
    (void)R;
    if (fmpz_poly_is_zero(RUP(a))) return;
    fmpz_t g; fmpz_init(g);
    fmpz_poly_content(g, RUP(a));
    if (fmpz_sgn(fmpz_poly_lead(RUP(a))) < 0) fmpz_neg(g, g);
    if (!fmpz_is_one(g)) fmpz_poly_scalar_divexact_fmpz(RUP(a), RUP(a), g);
    fmpz_clear(g);
}
static void ru_from(void* a, const fmpz_mpoly_t b, const RfRing* R) { fmpz_mpoly_get_fmpz_poly(RUP(a), b, 0, R->ctx); }
static void ru_to(fmpz_mpoly_t a, const void* b, const RfRing* R) { fmpz_mpoly_set_fmpz_poly(a, RUPc(b), 0, R->ctx); }

static void rf_ring_init(RfRing* R, const RfCtx* C) {
    R->ctx = C->ctx;
    if (C->n == 1) {
        R->sz = sizeof(fmpz_poly_struct);
        R->init = ru_init; R->clear = ru_clear; R->set = ru_set; R->one = ru_one;
        R->swap = ru_swap; R->is_const = ru_is_const; R->is_zero = ru_is_zero; R->is_linear = ru_is_linear;
        R->equal = ru_equal; R->mul = ru_mul; R->add = ru_add; R->pow = ru_pow;
        R->smul = ru_smul; R->gcd = ru_gcd; R->divides = ru_divides;
        R->normalise = ru_normalise; R->from_mpoly = ru_from; R->to_mpoly = ru_to;
    } else {
        R->sz = sizeof(fmpz_mpoly_struct);
        R->init = rm_init; R->clear = rm_clear; R->set = rm_set; R->one = rm_one;
        R->swap = rm_swap; R->is_const = rm_is_const; R->is_zero = rm_is_zero; R->is_linear = rm_is_linear;
        R->equal = rm_equal; R->mul = rm_mul; R->add = rm_add; R->pow = rm_pow;
        R->smul = rm_smul; R->gcd = rm_gcd; R->divides = rm_divides;
        R->normalise = rm_normalise; R->from_mpoly = rm_from; R->to_mpoly = rm_to;
    }
}

/* Storage for one element of either ring (aligned for both). */
typedef union { fmpz_mpoly_struct m; fmpz_poly_struct u; } RfElt;

/* A growable array of ring elements. */
typedef struct { char* v; slong n, cap; const RfRing* R; } RfVec;
#define RV(V, i) ((void*)((V)->v + (size_t)(i) * (V)->R->sz))

static void rv_init(RfVec* V, const RfRing* R) { V->v = NULL; V->n = V->cap = 0; V->R = R; }
static slong rv_push(RfVec* V) {
    if (V->n == V->cap) {
        slong nc = V->cap ? 2 * V->cap : 8;
        char* nv = realloc(V->v, (size_t)nc * V->R->sz);
        if (!nv) return -1;
        V->v = nv; V->cap = nc;
    }
    V->R->init(RV(V, V->n), V->R);
    return V->n++;
}
static void rv_clear(RfVec* V) {
    for (slong i = 0; i < V->n; i++) V->R->clear(RV(V, i), V->R);
    free(V->v);
    V->v = NULL; V->n = V->cap = 0;
}

/* A gcd-free basis: elements (squarefree, primitive, normalised) and, for
 * each, the distinct denominator polys D_j it divides. */
typedef struct { RfVec el; slong** own; slong* nown; slong* cown; slong ocap; } RfBasisR;

static void rb_init(RfBasisR* B, const RfRing* R) {
    rv_init(&B->el, R);
    B->own = NULL; B->nown = NULL; B->cown = NULL; B->ocap = 0;
}

static int rb_own_add(RfBasisR* B, slong k, slong j) {
    for (slong i = 0; i < B->nown[k]; i++) if (B->own[k][i] == j) return 1;
    if (B->nown[k] == B->cown[k]) {
        slong nc = B->cown[k] ? 2 * B->cown[k] : 4;
        slong* no = realloc(B->own[k], (size_t)nc * sizeof(slong));
        if (!no) return 0;
        B->own[k] = no; B->cown[k] = nc;
    }
    B->own[k][B->nown[k]++] = j;
    return 1;
}

static slong rb_new(RfBasisR* B) {
    if (B->el.n == B->ocap) {
        slong nc = B->ocap ? 2 * B->ocap : 8;
        slong** o = realloc(B->own, (size_t)nc * sizeof(slong*));
        if (!o) return -1;
        B->own = o;
        slong* n1 = realloc(B->nown, (size_t)nc * sizeof(slong));
        if (!n1) return -1;
        B->nown = n1;
        slong* c1 = realloc(B->cown, (size_t)nc * sizeof(slong));
        if (!c1) return -1;
        B->cown = c1;
        B->ocap = nc;
    }
    slong k = rv_push(&B->el);
    if (k < 0) return -1;
    B->own[k] = NULL; B->nown[k] = B->cown[k] = 0;
    return k;
}

static void rb_clear(RfBasisR* B) {
    for (slong k = 0; k < B->el.n; k++) free(B->own[k]);
    free(B->own); free(B->nown); free(B->cown);
    rv_clear(&B->el);
}

/* Insert the squarefree normalised D (index j).  S (the product of the
 * basis, or NULL) proves "coprime to everything" with one gcd. */
static int rb_insert(RfBasisR* B, const void* D, slong j, void* S, RfCtx* C) {
    const RfRing* R = B->el.R;
    RfElt tq, tg, tr;
    void* q = &tq; void* g = &tg; void* r = &tr;
    R->init(q, R); R->init(g, R); R->init(r, R);
    R->set(q, D, R);
    int ok = 1, coprime = 0;
    if (S && B->el.n > 0) {
        if (!R->gcd(g, q, S, R)) { C->fail = 1; ok = 0; }
        else coprime = R->is_const(g, R);
    }
    slong n0 = B->el.n;
    int qlin = R->is_linear(q, R);
    for (slong k = 0; k < n0 && ok && !coprime && !R->is_const(q, R); k++) {
        /* distinct primitive linear polynomials are coprime */
        if (qlin && R->is_linear(RV(&B->el, k), R) && !R->equal(q, RV(&B->el, k), R)) continue;
        if (!R->gcd(g, q, RV(&B->el, k), R)) { C->fail = 1; ok = 0; break; }
        if (R->is_const(g, R)) continue;
        R->normalise(g, R);
        if (!R->divides(r, q, g, R)) { C->fail = 1; ok = 0; break; }
        R->swap(q, r, R);
        if (!R->equal(g, RV(&B->el, k), R)) {
            /* the element splits: it keeps g, a new element takes p/g */
            if (!R->divides(r, RV(&B->el, k), g, R)) { C->fail = 1; ok = 0; break; }
            slong nk = rb_new(B);
            if (nk < 0) { ok = 0; break; }
            R->normalise(r, R);
            R->swap(RV(&B->el, nk), r, R);
            for (slong o = 0; o < B->nown[k] && ok; o++) ok = rb_own_add(B, nk, B->own[k][o]);
            R->swap(RV(&B->el, k), g, R);
        }
        if (ok) ok = rb_own_add(B, k, j);
    }
    if (ok && !R->is_const(q, R)) {
        R->normalise(q, R);
        slong nk = rb_new(B);
        if (nk < 0) ok = 0;
        else {
            if (S) R->mul(S, S, q, R);
            R->swap(RV(&B->el, nk), q, R);
            ok = rb_own_add(B, nk, j);
        }
    }
    R->clear(q, R); R->clear(g, R); R->clear(r, R);
    return ok;
}

/* Replace every non-squarefree denominator factor by its squarefree parts. */

static int rf_den_squarefree(Rf* R, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    int ok = 1;
    for (slong i = R->n - 1; i >= 0 && ok; i--) {
        if (R->f[i].e >= 0 || rf_is_mono(&R->f[i].p)) continue;
        if (fmpz_mpoly_total_degree_si(&R->f[i].p, ctx) < 2) continue;
        fmpz_mpoly_factor_t sq;
        fmpz_mpoly_factor_init(sq, ctx);
        if (!fmpz_mpoly_factor_squarefree(sq, &R->f[i].p, ctx)) {
            fmpz_mpoly_factor_clear(sq, ctx);
            C->fail = 1;
            return 0;
        }
        if (!(sq->num == 1 && fmpz_is_one(sq->exp + 0))) {
            slong e = R->f[i].e;
            rf_remove(R, i, C);
            ok = rf_mul_fmpz_pow(R->c, sq->constant, e);
            for (slong s = 0; s < sq->num && ok; s++) {
                fmpz_mpoly_t f; fmpz_mpoly_init(f, ctx);
                fmpz_mpoly_set(f, sq->poly + s, ctx);
                ok = rf_push(R, f, e * fmpz_get_si(sq->exp + s), C);
                fmpz_mpoly_clear(f, ctx);
            }
            i = R->n;               /* the factor list was reordered: rescan */
        }
        fmpz_mpoly_factor_clear(sq, ctx);
    }
    return ok;
}

/* Common denominator of rationals: Lc = lcm of the denominators. */
static void rf_lcm_dens(fmpz_t Lc, Rf* T, slong nt) {
    fmpz_one(Lc);
    for (slong t = 0; t < nt; t++)
        if (!fmpq_is_zero(T[t].c)) fmpz_lcm(Lc, Lc, fmpq_denref(T[t].c));
}

/* out = sum of the terms, all polynomial (no denominator factors). */
static int rf_poly_sum(Rf* out, Rf* T, slong nt, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    fmpz_t Lc, k, cf; fmpz_init(Lc); fmpz_init(k); fmpz_init(cf);
    rf_lcm_dens(Lc, T, nt);
    ulong* ex = malloc((size_t)C->n * sizeof(ulong));
    fmpz_mpoly_t acc, P;
    fmpz_mpoly_init(acc, ctx); fmpz_mpoly_init(P, ctx);
    for (slong t = 0; t < nt && ex; t++) {
        if (fmpq_is_zero(T[t].c)) continue;
        fmpz_divexact(k, Lc, fmpq_denref(T[t].c));
        fmpz_mul(k, k, fmpq_numref(T[t].c));
        rf_expand_side(P, &T[t], 0, C);
        /* append the terms; one sort + combine at the end keeps this
         * O(N log N) in the total term count instead of a merge per term */
        for (slong i = 0; i < P->length; i++) {
            fmpz_mul(cf, P->coeffs + i, k);
            fmpz_mpoly_get_term_exp_ui(ex, P, i, ctx);
            fmpz_mpoly_push_term_fmpz_ui(acc, cf, ex, ctx);
        }
    }
    fmpz_mpoly_sort_terms(acc, ctx);
    fmpz_mpoly_combine_like_terms(acc, ctx);
    fmpq_one(out->c);
    fmpq_div_fmpz(out->c, out->c, Lc);
    int ok = ex != NULL;
    if (ok) {
        if (acc->length == 0) fmpq_zero(out->c);
        else ok = rf_push(out, acc, 1, C);
    }
    free(ex);
    fmpz_mpoly_clear(acc, ctx); fmpz_mpoly_clear(P, ctx);
    fmpz_clear(Lc); fmpz_clear(k); fmpz_clear(cf);
    return ok;
}

/* Sum with exactly one fractional term t (already reduced): the common
 * denominator is den_t itself and N = num_t + den_t * (the polynomial terms)
 * is coprime to it -- the polynomial part adds a multiple of den_t to a
 * numerator that shares nothing with it -- so no gcd is needed at all. */
static int rf_sum_one_den(Rf* out, Rf* T, slong nt, slong td, const fmpz_t Lc, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    fmpz_mpoly_t N, D, P;
    fmpz_mpoly_init(N, ctx); fmpz_mpoly_init(D, ctx); fmpz_mpoly_init(P, ctx);
    fmpz_t k; fmpz_init(k);
    rf_expand_side(D, &T[td], 1, C);
    for (slong t = 0; t < nt; t++) {
        if (fmpq_is_zero(T[t].c)) continue;
        fmpz_divexact(k, Lc, fmpq_denref(T[t].c));
        fmpz_mul(k, k, fmpq_numref(T[t].c));
        rf_expand_side(P, &T[t], 0, C);
        if (t != td) fmpz_mpoly_mul(P, P, D, ctx);
        fmpz_mpoly_scalar_mul_fmpz(P, P, k, ctx);
        fmpz_mpoly_add(N, N, P, ctx);
    }
    fmpq_one(out->c);
    fmpq_div_fmpz(out->c, out->c, Lc);
    int ok = rf_push(out, N, 1, C);
    for (slong i = 0; i < T[td].n && ok; i++) {
        if (T[td].f[i].e >= 0) continue;
        fmpz_mpoly_set(P, &T[td].f[i].p, ctx);
        ok = rf_push(out, P, T[td].f[i].e, C);
    }
    out->reduced = 1;
    fmpz_clear(k);
    fmpz_mpoly_clear(N, ctx); fmpz_mpoly_clear(D, ctx); fmpz_mpoly_clear(P, ctx);
    return ok;
}

/* Sum regime: out = sum of the terms (consumed). */
static int rf_sum(Rf* out, Rf* T, slong nt, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    int ok = 1;

    /* 1. each term on its own: cancel, monomial split of the denominator */
    slong nden = 0, td = -1;
    for (slong t = 0; t < nt && ok; t++) {
        if (fmpq_is_zero(T[t].c) || rf_is_poly(&T[t])) continue;
        ok = T[t].reduced ? rf_split_den_monomials(&T[t], C) : rf_single(&T[t], C);
        nden++; td = t;
    }
    if (!ok) return 0;
    fmpz_t Lc; fmpz_init(Lc);
    rf_lcm_dens(Lc, T, nt);
    if (nden == 1) {
        /* the factored form makes the lone denominator squarefree too:
         * Together[1/(x^2 + 2x + 1) + 1] -> (2 + 2x + x^2)/(1 + x)^2 */
        ok = (!rf_factored_form() || rf_den_squarefree(&T[td], C)) &&
             rf_sum_one_den(out, T, nt, td, Lc, C);
        fmpz_clear(Lc);
        return ok && !C->fail;
    }

    /* 2. squarefree parts of the denominator factors */
    for (slong t = 0; t < nt && ok; t++)
        if (!fmpq_is_zero(T[t].c) && !rf_is_poly(&T[t])) ok = rf_den_squarefree(&T[t], C);
    if (!ok) { fmpz_clear(Lc); return 0; }

    /* 3. the distinct denominator polys D_j; per term the (j, m) it uses */
    RfRing R; rf_ring_init(&R, C);
    RfVec Dq; rv_init(&Dq, &R);
    slong** tj = calloc((size_t)nt, sizeof(slong*));
    slong** tm = calloc((size_t)nt, sizeof(slong*));
    slong* tn = calloc((size_t)nt, sizeof(slong));
    fmpz_mpoly_struct** Dsrc = NULL;           /* borrowed: Dq[j] came from here */
    slong ndsrc = 0, cdsrc = 0;
    if (!tj || !tm || !tn) ok = 0;
    for (slong t = 0; t < nt && ok; t++) {
        if (fmpq_is_zero(T[t].c)) continue;
        slong cnt = 0;
        for (slong i = 0; i < T[t].n; i++) if (T[t].f[i].e < 0) cnt++;
        if (!cnt) continue;
        tj[t] = malloc((size_t)cnt * sizeof(slong));
        tm[t] = malloc((size_t)cnt * sizeof(slong));
        if (!tj[t] || !tm[t]) { ok = 0; break; }
        for (slong i = 0; i < T[t].n && ok; i++) {
            if (T[t].f[i].e >= 0) continue;
            slong j = -1;
            for (slong d = 0; d < ndsrc; d++)
                if (fmpz_mpoly_equal(Dsrc[d], &T[t].f[i].p, ctx)) { j = d; break; }
            if (j < 0) {
                if (ndsrc == cdsrc) {
                    slong nc = cdsrc ? 2 * cdsrc : 16;
                    fmpz_mpoly_struct** nd = realloc(Dsrc, (size_t)nc * sizeof(*nd));
                    if (!nd) { ok = 0; break; }
                    Dsrc = nd; cdsrc = nc;
                }
                Dsrc[ndsrc] = &T[t].f[i].p;
                j = ndsrc++;
                slong q = rv_push(&Dq);
                if (q < 0) { ok = 0; break; }
                R.from_mpoly(RV(&Dq, q), &T[t].f[i].p, &R);
            }
            tj[t][tn[t]] = j;
            tm[t][tn[t]] = -T[t].f[i].e;
            tn[t]++;
        }
    }
    free(Dsrc);
    slong nd = Dq.n;

    /* 4. gcd-free basis; S filters the coprime insertions when it is large */
    RfBasisR B; rb_init(&B, &R);
    RfElt tS; void* S = &tS;
    R.init(S, &R); R.one(S, &R);
    int use_S = nd > 6, coprime = 0;
    int small = 1;                      /* the shortcut pays for small factors only */
    for (slong j = 0; j < nd && small; j++)
        if (C->n == 1 && fmpz_poly_degree((fmpz_poly_struct*)RV(&Dq, j)) > 4) small = 0;
    if (nd > 6 && C->n == 1 && small) {
        /* One variable: the D_j are squarefree and distinct, so they are
         * pairwise coprime exactly when their product is squarefree -- one gcd
         * with the derivative instead of a gcd per insertion. */
        fmpz_poly_t Pp, Pd, G;
        fmpz_poly_init(Pp); fmpz_poly_init(Pd); fmpz_poly_init(G);
        fmpz_poly_one(Pp);
        for (slong j = 0; j < nd; j++) fmpz_poly_mul(Pp, Pp, (fmpz_poly_struct*)RV(&Dq, j));
        fmpz_poly_derivative(Pd, Pp);
        fmpz_poly_gcd(G, Pp, Pd);
        coprime = fmpz_poly_length(G) <= 1;
        fmpz_poly_clear(Pp); fmpz_poly_clear(Pd); fmpz_poly_clear(G);
    }
    if (coprime) {
        for (slong j = 0; j < nd && ok; j++) {
            slong k = rb_new(&B);
            if (k < 0) { ok = 0; break; }
            R.set(RV(&B.el, k), RV(&Dq, j), &R);
            ok = rb_own_add(&B, k, j);
        }
    } else {
        for (slong j = 0; j < nd && ok; j++) ok = rb_insert(&B, RV(&Dq, j), j, use_S ? S : NULL, C);
    }

    /* 5. multiplicities: per term its (k, m) list; M[k]; how many terms use k */
    slong nb = B.el.n;
    slong** jk = calloc((size_t)(nd ? nd : 1), sizeof(slong*));
    slong* njk = calloc((size_t)(nd ? nd : 1), sizeof(slong));
    slong* M = calloc((size_t)(nb ? nb : 1), sizeof(slong));
    slong* mt = calloc((size_t)(nb ? nb : 1), sizeof(slong));
    slong* users = calloc((size_t)(nb ? nb : 1), sizeof(slong));
    if (!jk || !njk || !M || !mt || !users) ok = 0;
    for (slong k = 0; k < nb && ok; k++)
        for (slong o = 0; o < B.nown[k]; o++) njk[B.own[k][o]]++;
    for (slong j = 0; j < nd && ok; j++) {
        jk[j] = malloc((size_t)(njk[j] ? njk[j] : 1) * sizeof(slong));
        if (!jk[j]) ok = 0;
        njk[j] = 0;
    }
    for (slong k = 0; k < nb && ok; k++)
        for (slong o = 0; o < B.nown[k]; o++) { slong j = B.own[k][o]; jk[j][njk[j]++] = k; }
    for (slong t = 0; t < nt && ok; t++) {
        for (slong i = 0; i < tn[t]; i++)
            for (slong q = 0; q < njk[tj[t][i]]; q++) mt[jk[tj[t][i]][q]] += tm[t][i];
        for (slong i = 0; i < tn[t]; i++)
            for (slong q = 0; q < njk[tj[t][i]]; q++) {
                slong k = jk[tj[t][i]][q];
                if (mt[k]) { users[k]++; if (mt[k] > M[k]) M[k] = mt[k]; mt[k] = 0; }
            }
    }
    int disjoint = 1;
    for (slong k = 0; k < nb; k++) if (users[k] > 1) disjoint = 0;
    /* single[t]: the one basis element of term t's denominator (-1: none;
     * -2: more than one).  When no term has two, the terms are summed per
     * element first and the per-element sums are disjoint. */
    slong* single = calloc((size_t)(nt ? nt : 1), sizeof(slong));
    slong* sm = calloc((size_t)(nt ? nt : 1), sizeof(slong));
    int grouped = !disjoint && single && sm;
    for (slong t = 0; t < nt && grouped; t++) {
        single[t] = -1;
        for (slong i = 0; i < tn[t]; i++)
            for (slong q = 0; q < njk[tj[t][i]]; q++) {
                slong k = jk[tj[t][i]][q];
                if (single[t] == -1) { single[t] = k; sm[t] = tm[t][i]; }
                else if (single[t] == k) sm[t] += tm[t][i];
                else single[t] = -2;
            }
        if (single[t] == -2) grouped = 0;
    }

    /* 6. numerator.  Leaves (k_t num_t, den_t); disjoint denominators sum by a
     * product tree (P1 Q2 + P2 Q1, Q1 Q2), otherwise against L = prod B^M. */
    RfVec P, Q; rv_init(&P, &R); rv_init(&Q, &R);
    RfElt tN, tt, tL, tc; void* N = &tN; void* tmp = &tt; void* Lp = &tL; void* cof = &tc;
    R.init(N, &R); R.init(tmp, &R); R.init(Lp, &R); R.init(cof, &R);
    fmpz_mpoly_t num; fmpz_mpoly_init(num, ctx);
    fmpz_t kc; fmpz_init(kc);
    for (slong t = 0; t < nt && ok; t++) {
        if (fmpq_is_zero(T[t].c)) continue;
        slong ip = rv_push(&P), iq = rv_push(&Q);
        if (ip < 0 || iq < 0) { ok = 0; break; }
        rf_expand_side(num, &T[t], 0, C);
        R.from_mpoly(RV(&P, ip), num, &R);
        fmpz_divexact(kc, Lc, fmpq_denref(T[t].c));
        fmpz_mul(kc, kc, fmpq_numref(T[t].c));
        R.smul(RV(&P, ip), RV(&P, ip), kc, &R);
        R.one(RV(&Q, iq), &R);
        for (slong i = 0; i < tn[t]; i++)
            for (slong q = 0; q < njk[tj[t][i]]; q++) {
                slong kk = jk[tj[t][i]][q];
                if (tm[t][i] == 1) R.mul(RV(&Q, iq), RV(&Q, iq), RV(&B.el, kk), &R);
                else {
                    R.pow(tmp, RV(&B.el, kk), (ulong)tm[t][i], &R);
                    R.mul(RV(&Q, iq), RV(&Q, iq), tmp, &R);
                }
            }
    }
    if (ok && grouped) {
        /* leaves per element: (sum_t P_t b^(M - m_t), b^M); the polynomial
         * terms (no element) join with Q = 1 */
        RfVec GP, GQ; rv_init(&GP, &R); rv_init(&GQ, &R);
        slong* gi = malloc((size_t)(nb + 1) * sizeof(slong));
        for (slong k = 0; k <= nb; k++) gi[k] = -1;
        slong li = 0;
        for (slong t = 0; t < nt && ok; t++) {
            if (fmpq_is_zero(T[t].c)) continue;
            slong i = li++;                      /* the leaf of term t */
            slong k = single[t] < 0 ? nb : single[t];
            if (gi[k] < 0) {
                gi[k] = rv_push(&GP); slong q = rv_push(&GQ);
                if (gi[k] < 0 || q < 0) { ok = 0; break; }
                if (k == nb) R.one(RV(&GQ, q), &R);
                else R.pow(RV(&GQ, q), RV(&B.el, k), (ulong)M[k], &R);
            }
            if (k < nb && sm[t] < M[k]) {
                R.pow(tmp, RV(&B.el, k), (ulong)(M[k] - sm[t]), &R);
                R.mul(tmp, RV(&P, i), tmp, &R);
                R.add(RV(&GP, gi[k]), RV(&GP, gi[k]), tmp, &R);
            } else R.add(RV(&GP, gi[k]), RV(&GP, gi[k]), RV(&P, i), &R);
        }
        free(gi);
        rv_clear(&P); rv_clear(&Q);
        P = GP; Q = GQ;
    }
    if (ok && grouped) {
        slong n = P.n;
        while (n > 1) {
            slong h = 0;
            for (slong i = 0; i + 1 < n; i += 2, h++) {
                R.mul(tmp, RV(&P, i), RV(&Q, i + 1), &R);
                R.mul(RV(&P, i + 1), RV(&P, i + 1), RV(&Q, i), &R);
                R.add(RV(&P, h), tmp, RV(&P, i + 1), &R);
                R.mul(RV(&Q, h), RV(&Q, i), RV(&Q, i + 1), &R);
            }
            if (n & 1) { R.swap(RV(&P, h), RV(&P, n - 1), &R); R.swap(RV(&Q, h), RV(&Q, n - 1), &R); h++; }
            n = h;
        }
        if (n == 1) R.swap(N, RV(&P, 0), &R);
    } else if (ok) {
        /* A product tree whose denominators are kept as sparse exponent
         * vectors over the basis: merging P1/B^e1 + P2/B^e2 multiplies each
         * side by just the basis powers the other one has more of, so there is
         * neither a division nor a full denominator until the end. */
        slong nl = P.n;
        slong** ek = malloc((size_t)(nl ? nl : 1) * sizeof(slong*));   /* sorted keys */
        slong** em = malloc((size_t)(nl ? nl : 1) * sizeof(slong*));   /* exponents  */
        slong* en = calloc((size_t)(nl ? nl : 1), sizeof(slong));
        if (!ek || !em || !en) ok = 0;
        slong li = 0;
        for (slong t = 0; t < nt && ok; t++) {
            if (fmpq_is_zero(T[t].c)) continue;
            slong cnt = 0;
            for (slong i = 0; i < tn[t]; i++) cnt += njk[tj[t][i]];
            ek[li] = malloc((size_t)(cnt ? cnt : 1) * sizeof(slong));
            em[li] = malloc((size_t)(cnt ? cnt : 1) * sizeof(slong));
            if (!ek[li] || !em[li]) { ok = 0; break; }
            for (slong i = 0; i < tn[t]; i++)
                for (slong q = 0; q < njk[tj[t][i]]; q++) {
                    slong k = jk[tj[t][i]][q], m = tm[t][i], z = en[li];
                    /* insertion into the sorted key list, merging a repeat */
                    slong pos = 0;
                    while (pos < z && ek[li][pos] < k) pos++;
                    if (pos < z && ek[li][pos] == k) { em[li][pos] += m; continue; }
                    for (slong w = z; w > pos; w--) { ek[li][w] = ek[li][w - 1]; em[li][w] = em[li][w - 1]; }
                    ek[li][pos] = k; em[li][pos] = m; en[li]++;
                }
            li++;
        }
        slong n = nl;
        while (ok && n > 1) {
            slong h = 0;
            for (slong i = 0; i + 1 < n && ok; i += 2, h++) {
                slong na = en[i], nb2 = en[i + 1];
                slong* rk = malloc((size_t)(na + nb2 + 1) * sizeof(slong));
                slong* rm = malloc((size_t)(na + nb2 + 1) * sizeof(slong));
                if (!rk || !rm) { free(rk); free(rm); ok = 0; break; }
                R.one(Lp, &R); R.one(cof, &R);          /* Lp: cofactor of side i, cof: of i+1 */
                slong x = 0, y = 0, z = 0;
                while (x < na || y < nb2) {
                    slong k, ma = 0, mb = 0;
                    if (y >= nb2 || (x < na && ek[i][x] < ek[i + 1][y])) { k = ek[i][x]; ma = em[i][x++]; }
                    else if (x >= na || ek[i + 1][y] < ek[i][x]) { k = ek[i + 1][y]; mb = em[i + 1][y++]; }
                    else { k = ek[i][x]; ma = em[i][x++]; mb = em[i + 1][y++]; }
                    slong mx = ma > mb ? ma : mb;
                    if (mx > ma) {
                        if (mx - ma == 1) R.mul(Lp, Lp, RV(&B.el, k), &R);
                        else { R.pow(tmp, RV(&B.el, k), (ulong)(mx - ma), &R); R.mul(Lp, Lp, tmp, &R); }
                    }
                    if (mx > mb) {
                        if (mx - mb == 1) R.mul(cof, cof, RV(&B.el, k), &R);
                        else { R.pow(tmp, RV(&B.el, k), (ulong)(mx - mb), &R); R.mul(cof, cof, tmp, &R); }
                    }
                    rk[z] = k; rm[z] = mx; z++;
                }
                R.mul(tmp, RV(&P, i), Lp, &R);
                R.mul(RV(&P, i + 1), RV(&P, i + 1), cof, &R);
                R.add(RV(&P, h), tmp, RV(&P, i + 1), &R);
                free(ek[i]); free(em[i]); free(ek[i + 1]); free(em[i + 1]);
                ek[h] = rk; em[h] = rm; en[h] = z;
            }
            if (ok && (n & 1)) {
                R.swap(RV(&P, h), RV(&P, n - 1), &R);
                ek[h] = ek[n - 1]; em[h] = em[n - 1]; en[h] = en[n - 1];
                h++;
            }
            n = h;
        }
        if (ok && n == 1) {
            R.swap(N, RV(&P, 0), &R);
            /* the root's exponents are M; elements no term reached stay 0 */
            free(ek[0]); free(em[0]);
        }
        free(ek); free(em); free(en);
    }

    /* 7. cancel N against the denominator, element by element */
    RfVec E; rv_init(&E, &R);
    slong* EM = NULL; slong cem = 0;
    if (ok && !R.is_zero(N, &R)) {
        for (slong b = 0; b < nb && ok; b++) {
            if (M[b] == 0) continue;
            slong i = rv_push(&E);
            if (i < 0 || (i >= cem && !(EM = realloc(EM, (size_t)(cem = 2 * cem + 8) * sizeof(slong))))) { ok = 0; break; }
            R.swap(RV(&E, i), RV(&B.el, b), &R);
            EM[i] = M[b];
        }
        int maybe = ok && !R.is_const(N, &R);
        if (maybe) {
            /* one gcd with the squarefree denominator decides "nothing cancels" */
            R.one(tmp, &R);
            for (slong i = 0; i < E.n; i++) R.mul(tmp, tmp, RV(&E, i), &R);
            if (!R.gcd(cof, N, tmp, &R)) { C->fail = 1; ok = 0; }
            maybe = ok && !R.is_const(cof, &R);
        }
        for (slong i = 0; i < E.n && ok && maybe; i++) {
            /* element b^M: with g = gcd(N, b), b/g is coprime to N (b is
             * squarefree); divide N by g while it can, up to M times, and
             * split g again by whatever of it still divides N */
            while (ok && EM[i] > 0) {
                if (!R.gcd(cof, N, RV(&E, i), &R)) { C->fail = 1; ok = 0; break; }
                if (R.is_const(cof, &R)) break;
                R.normalise(cof, &R);
                if (!R.equal(cof, RV(&E, i), &R)) {
                    if (!R.divides(tmp, RV(&E, i), cof, &R)) { C->fail = 1; ok = 0; break; }
                    slong ni = rv_push(&E);
                    if (ni < 0 || (ni >= cem && !(EM = realloc(EM, (size_t)(cem = 2 * cem + 8) * sizeof(slong))))) { ok = 0; break; }
                    R.normalise(tmp, &R);
                    R.swap(RV(&E, ni), tmp, &R);
                    EM[ni] = EM[i];
                    R.set(RV(&E, i), cof, &R);
                }
                if (!R.divides(tmp, N, RV(&E, i), &R)) { C->fail = 1; ok = 0; break; }
                R.swap(N, tmp, &R);
                EM[i]--;
            }
        }
    }

    /* 8. assemble c * N * prod E_i^-M_i */
    if (ok) {
        fmpq_one(out->c);
        fmpq_div_fmpz(out->c, out->c, Lc);
        if (R.is_zero(N, &R)) fmpq_zero(out->c);
        else {
            R.to_mpoly(num, N, &R);
            ok = rf_push(out, num, 1, C);
            for (slong i = 0; i < E.n && ok; i++) {
                if (EM[i] <= 0) continue;
                R.to_mpoly(num, RV(&E, i), &R);
                ok = rf_push(out, num, -EM[i], C);
            }
        }
        out->reduced = 1;
    }

    free(EM);
    rv_clear(&E); rv_clear(&P); rv_clear(&Q);
    R.clear(N, &R); R.clear(tmp, &R); R.clear(Lp, &R); R.clear(cof, &R); R.clear(S, &R);
    fmpz_mpoly_clear(num, ctx);
    fmpz_clear(kc); fmpz_clear(Lc);
    for (slong j = 0; j < nd; j++) if (jk) free(jk[j]);
    free(jk); free(njk); free(M); free(mt); free(users); free(single); free(sm);
    rb_clear(&B);
    rv_clear(&Dq);
    for (slong t = 0; t < nt; t++) { if (tj) free(tj[t]); if (tm) free(tm[t]); }
    free(tj); free(tm); free(tn);
    return ok && !C->fail;
}

/* ------------------------------------------------------------------ */
/*  Expr -> Rf                                                         */
/* ------------------------------------------------------------------ */

static int rf_from_expr(const Expr* e, Rf* R, RfCtx* C);

/* Exponent vector of a monomial term c * atom^k * ... (c a number, k >= 0);
 * 0 when the term has any other shape. */
static int rf_monomial_term(const Expr* t, const RfCtx* C, ulong* ex, const Expr** coef) {
    for (size_t v = 0; v < C->n; v++) ex[v] = 0;
    *coef = NULL;
    size_t n = 1;
    const Expr* const* fs = &t;
    if (rf_head_is(t, SYM_Times)) { n = t->data.function.arg_count; fs = (const Expr* const*)t->data.function.args; }
    for (size_t i = 0; i < n; i++) {
        const Expr* f = fs[i];
        if (f->type == EXPR_INTEGER || f->type == EXPR_BIGINT || rf_head_is(f, SYM_Rational)) {
            if (*coef || i != 0) return 0;
            *coef = f;
            continue;
        }
        const Expr* b = f;
        int64_t k = 1;
        if (rf_head_is(f, SYM_Power) && f->data.function.arg_count == 2 &&
            f->data.function.args[1]->type == EXPR_INTEGER) {
            b = f->data.function.args[0];
            k = f->data.function.args[1]->data.integer;
            if (k < 0) return 0;
        }
        if (b->type != EXPR_SYMBOL) return 0;      /* kernels take the general path */
        slong v = rf_atom_index(C, b);
        if (v < 0) return 0;
        ex[v] += (ulong)k;
    }
    return 1;
}

/* A polynomial sum of monomial terms straight into one fmpz_mpoly (one sort
 * and combine); 0 when some term is not a monomial term. */
static int rf_plus_direct(const Expr* e, Rf* R, RfCtx* C) {
    size_t nt = e->data.function.arg_count;
    ulong* ex = malloc(C->n * sizeof(ulong));
    if (!ex) return 0;
    fmpz_t Lc; fmpz_init(Lc); fmpz_one(Lc);
    fmpq_t q; fmpq_init(q);
    const Expr* cf;
    int ok = 1;
    for (size_t t = 0; t < nt && ok; t++) {
        ok = rf_monomial_term(e->data.function.args[t], C, ex, &cf);
        if (ok && cf && rf_head_is(cf, SYM_Rational)) {
            ok = fb_expr_to_fmpq(cf, q);
            if (ok) fmpz_lcm(Lc, Lc, fmpq_denref(q));
        }
    }
    if (ok) {
        fmpz_mpoly_t acc; fmpz_mpoly_init2(acc, (slong)nt, C->ctx);
        fmpz_t k; fmpz_init(k);
        for (size_t t = 0; t < nt; t++) {
            rf_monomial_term(e->data.function.args[t], C, ex, &cf);
            if (!cf) fmpz_set(k, Lc);
            else if (cf->type == EXPR_INTEGER) fmpz_mul_si(k, Lc, cf->data.integer);
            else if (cf->type == EXPR_BIGINT) { fmpz_set_mpz(k, cf->data.bigint); fmpz_mul(k, k, Lc); }
            else {
                fb_expr_to_fmpq(cf, q);
                fmpz_divexact(k, Lc, fmpq_denref(q));
                fmpz_mul(k, k, fmpq_numref(q));
            }
            fmpz_mpoly_push_term_fmpz_ui(acc, k, ex, C->ctx);
        }
        fmpz_mpoly_sort_terms(acc, C->ctx);
        fmpz_mpoly_combine_like_terms(acc, C->ctx);
        if (R->n) R->reduced = 0;
        fmpq_div_fmpz(R->c, R->c, Lc);
        if (acc->length == 0) fmpq_zero(R->c);
        else ok = rf_push(R, acc, 1, C);
        fmpz_clear(k);
        fmpz_mpoly_clear(acc, C->ctx);
    }
    fmpq_clear(q); fmpz_clear(Lc);
    free(ex);
    return ok;
}

static int rf_from_plus(const Expr* e, Rf* R, RfCtx* C) {
    if (rf_plus_direct(e, R, C)) return 1;
    if (C->fail) return 0;
    slong nt = (slong)e->data.function.arg_count;
    Rf* T = malloc((size_t)nt * sizeof(Rf));
    if (!T) return 0;
    int ok = 1, poly = 1;
    slong done = 0;
    for (slong t = 0; t < nt && ok; t++) {
        rf_init(&T[t]);
        done = t + 1;
        ok = rf_from_expr(e->data.function.args[t], &T[t], C);
        if (ok && !rf_is_poly(&T[t])) poly = 0;
    }
    Rf S; rf_init(&S);
    if (ok) ok = poly ? rf_poly_sum(&S, T, nt, C) : rf_sum(&S, T, nt, C);
    if (ok) ok = rf_mul_into(R, &S, C);
    rf_clear(&S, C);
    for (slong t = 0; t < done; t++) rf_clear(&T[t], C);
    free(T);
    return ok;
}

/* R *= e. */
static int rf_from_expr(const Expr* e, Rf* R, RfCtx* C) {
    if (C->fail) return 0;
    switch (e->type) {
        case EXPR_INTEGER:
            fmpq_mul_si(R->c, R->c, e->data.integer);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            fmpq_mul_fmpz(R->c, R->c, z);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_SYMBOL: {
            slong v = rf_atom_index(C, e);
            if (R->n) R->reduced = 0;
            return v >= 0 && rf_push_gen(R, v, 1, C);
        }
        case EXPR_FUNCTION: {
            const char* h = e->data.function.head->type == EXPR_SYMBOL
                          ? e->data.function.head->data.symbol.name : NULL;
            size_t n = e->data.function.arg_count;
            if (h == SYM_Rational) {
                fmpq_t q; fmpq_init(q);
                int ok = fb_expr_to_fmpq(e, q);
                if (ok) fmpq_mul(R->c, R->c, q);
                fmpq_clear(q);
                return ok;
            }
            if (h == SYM_Times) {
                for (size_t i = 0; i < n; i++)
                    if (!rf_from_expr(e->data.function.args[i], R, C)) return 0;
                return 1;
            }
            if (h == SYM_Plus) return rf_from_plus(e, R, C);
            if (h == SYM_Power && n == 2) {
                const Expr* b = e->data.function.args[0];
                const Expr* x = e->data.function.args[1];
                if (x->type == EXPR_INTEGER) {
                    Rf S; rf_init(&S);
                    int ok = rf_from_expr(b, &S, C) && rf_pow(&S, (slong)x->data.integer, C) &&
                             rf_mul_into(R, &S, C);
                    rf_clear(&S, C);
                    return ok;
                }
                if (R->n) R->reduced = 0;
                for (size_t i = 0; i < C->xn; i++)
                    if (expr_eq(C->xk[i], e)) return rf_push_gen(R, C->gidx, (slong)C->xpow[i], C);
                return 0;
            }
            slong v = rf_atom_index(C, e);
            if (R->n) R->reduced = 0;
            return v >= 0 && rf_push_gen(R, v, 1, C);
        }
        default:
            return 0;
    }
}

/* ------------------------------------------------------------------ */
/*  Rf -> Expr                                                         */
/* ------------------------------------------------------------------ */

static Expr* rf_fmpz_to_expr(const fmpz_t z) {
    if (fmpz_fits_si(z)) return expr_new_integer((int64_t)fmpz_get_si(z));
    mpz_t m; mpz_init(m);
    fmpz_get_mpz(m, z);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
    mpz_clear(m);
    return r;
}

static Expr* rf_atom_expr(const RfCtx* C, slong v) { return expr_copy((Expr*)C->atom[v]); }

/* The polynomial as an (unevaluated) Plus of Times terms. */
static Expr* rf_poly_to_expr(const fmpz_mpoly_t p, const RfCtx* C) {
    slong len = p->length;
    if (len == 0) return expr_new_integer(0);
    slong nv = (slong)C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    Expr** terms = malloc((size_t)len * sizeof(Expr*));
    Expr** fs = malloc((size_t)(nv + 1) * sizeof(Expr*));
    for (slong i = 0; i < len; i++) {
        fmpz_mpoly_get_term_exp_ui(ex, p, i, C->ctx);
        size_t nf = 0;
        if (!fmpz_is_one(p->coeffs + i)) fs[nf++] = rf_fmpz_to_expr(p->coeffs + i);
        for (slong v = 0; v < nv; v++) {
            if (!ex[v]) continue;
            Expr* a = rf_atom_expr(C, v);
            fs[nf++] = ex[v] == 1 ? a
                : expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ a, expr_new_integer((int64_t)ex[v]) }, 2);
        }
        if (nf == 0) terms[i] = expr_new_integer(1);
        else if (nf == 1) terms[i] = fs[0];
        else terms[i] = expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
    }
    Expr* r = len == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)len);
    free(ex); free(terms); free(fs);
    return r;
}

/* Index of the factor that absorbs a leading -1, or -1: the one non-monomial
 * numerator factor of exponent 1; else, with no non-monomial numerator factor,
 * the one non-monomial denominator factor of exponent -1. */
static slong rf_absorb_slot(const Rf* R) {
    slong nn = 0, nd = 0, in = -1, id = -1;
    for (slong i = 0; i < R->n; i++) {
        if (rf_is_mono(&R->f[i].p)) continue;
        if (R->f[i].e > 0) { nn++; in = i; } else { nd++; id = i; }
    }
    if (nn == 1) return R->f[in].e == 1 ? in : -1;
    if (nn == 0 && nd == 1) return R->f[id].e == -1 ? id : -1;
    return -1;
}

/* Times[c, factors...] (unevaluated).  `extra` (consumed, may be NULL) is one
 * more factor appended after the sign decision. */
static Expr* rf_to_expr_raw(const Rf* R, const RfCtx* C, Expr* extra) {
    if (fmpq_is_zero(R->c)) { if (extra) expr_free(extra); return expr_new_integer(0); }
    slong absorb = -1;
    if (fmpz_equal_si(fmpq_numref(R->c), -1) && fmpz_is_one(fmpq_denref(R->c)))
        absorb = rf_absorb_slot(R);
    Expr** a = malloc((size_t)(R->n + 2) * sizeof(Expr*));
    size_t na = 0;
    if (absorb < 0 && !fmpq_is_one(R->c)) {
        fmpq_t c; fmpq_init(c); fmpq_set(c, R->c);
        a[na++] = fb_fmpq_to_expr(c);
        fmpq_clear(c);
    }
    for (slong i = 0; i < R->n; i++) {
        Expr* b = rf_poly_to_expr(&R->f[i].p, C);
        if (i == absorb)
            b = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), b }, 2));
        a[na++] = R->f[i].e == 1 ? b
            : expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ b, expr_new_integer((int64_t)R->f[i].e) }, 2);
    }
    if (extra) a[na++] = extra;
    Expr* r = na == 0 ? expr_new_integer(1)
            : na == 1 ? a[0]
            : expr_new_function(expr_new_symbol(SYM_Times), a, na);
    free(a);
    return r;
}

/* Classical form: (p * num) / (q * den), both sides expanded. */
static Expr* rf_to_expr_expanded(const Rf* R, const RfCtx* C) {
    if (fmpq_is_zero(R->c)) return expr_new_integer(0);
    fmpz_mpoly_t N, D;
    fmpz_mpoly_init(N, C->ctx); fmpz_mpoly_init(D, C->ctx);
    rf_expand_side(N, R, 0, C);
    rf_expand_side(D, R, 1, C);
    fmpz_mpoly_scalar_mul_fmpz(N, N, fmpq_numref(R->c), C->ctx);
    fmpz_mpoly_scalar_mul_fmpz(D, D, fmpq_denref(R->c), C->ctx);
    /* The classical sign: flip num and den when every term of the
     * denominator's top total degree is negative (5 - x^2 -> x^2 - 5). */
    if (!fmpz_mpoly_is_fmpz(D, C->ctx)) {
        slong top = fmpz_mpoly_total_degree_si(D, C->ctx), nv = (slong)C->n;
        ulong* ex = malloc((size_t)nv * sizeof(ulong));
        int flip = ex != NULL;
        for (slong i = 0; i < D->length && flip; i++) {
            fmpz_mpoly_get_term_exp_ui(ex, D, i, C->ctx);
            slong td = 0;
            for (slong v = 0; v < nv; v++) td += (slong)ex[v];
            if (td == top && fmpz_sgn(D->coeffs + i) > 0) flip = 0;
        }
        free(ex);
        if (flip) { fmpz_mpoly_neg(N, N, C->ctx); fmpz_mpoly_neg(D, D, C->ctx); }
    }
    Expr* ne = rf_poly_to_expr(N, C);
    Expr* r;
    if (fmpz_mpoly_is_one(D, C->ctx)) r = ne;
    else r = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ ne, expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ rf_poly_to_expr(D, C), expr_new_integer(-1) }, 2) }, 2);
    fmpz_mpoly_clear(N, C->ctx); fmpz_mpoly_clear(D, C->ctx);
    return eval_and_free(r);
}

static Expr* rf_to_expr(const Rf* R, const RfCtx* C) {
    if (!rf_factored_form()) return rf_to_expr_expanded(R, C);
    return eval_and_free(rf_to_expr_raw(R, C, NULL));
}

/* ------------------------------------------------------------------ */
/*  Together / Cancel                                                  */
/* ------------------------------------------------------------------ */

/* A denominator to work on: a negative integer power of a non-number, or a
 * Rational coefficient. */
static int rf_has_den(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return 0;
    if (rf_head_is(e, SYM_Rational)) return 1;
    if (rf_head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER && x->data.integer < 0 &&
            (b->type == EXPR_SYMBOL || b->type == EXPR_FUNCTION))
            return 1;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rf_has_den(e->data.function.args[i])) return 1;
    return 0;
}

/* Together/Cancel of one term (a non-Plus, or any expression for Together). */
static Expr* rf_normal_form(const Expr* e, RfCtx* C) {
    Rf R; rf_init(&R);
    Expr* out = NULL;
    if (rf_from_expr(e, &R, C) &&
        (rf_head_is(e, SYM_Plus) ||
         (R.reduced ? rf_split_den_monomials(&R, C) : rf_single(&R, C))) && !C->fail)
        out = rf_to_expr(&R, C);
    rf_clear(&R, C);
    return out;
}

Expr* ratfun_together(const Expr* e) {
    if (!e) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, NULL)) return NULL;
    /* Nothing to combine: a polynomial (or a product of them) is left as it
     * is, factored or expanded -- Together[(1 + x)^2 (2 + x)] stays. */
    Expr* out = rf_has_den(e) ? rf_normal_form(e, &C)
              : rf_factored_form() ? expr_copy((Expr*)e) : NULL;
    rf_ctx_clear(&C);
    return out;
}

/* A sum carrying a denominator strictly inside `e` (below its top node). */
static int rf_nested_sum_den(const Expr* e, int top) {
    if (e->type != EXPR_FUNCTION) return 0;
    if (!top && rf_head_is(e, SYM_Plus)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            const Expr* t = e->data.function.args[i];
            if (rf_has_den(t) && !(t->type == EXPR_FUNCTION && rf_head_is(t, SYM_Rational)))
                return 1;
            /* a Rational coefficient alone is no nested fraction */
            if (t->type == EXPR_FUNCTION && rf_head_is(t, SYM_Times))
                for (size_t k = 0; k < t->data.function.arg_count; k++) {
                    const Expr* u = t->data.function.args[k];
                    if (rf_has_den(u) && !rf_head_is(u, SYM_Rational)) return 1;
                }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rf_nested_sum_den(e->data.function.args[i], 0)) return 1;
    return 0;
}

Expr* ratfun_cancel(const Expr* e) {
    if (!e || !rf_has_den(e)) return NULL;
    /* Cancel does not combine a nested sum of fractions
     * (Cancel[1/(1 + 1/x)] stays); leave those to the classical path. */
    if (rf_head_is(e, SYM_Plus)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (rf_nested_sum_den(e->data.function.args[i], 1)) return NULL;
    } else if (rf_nested_sum_den(e, 1)) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, NULL)) return NULL;
    Expr* out = NULL;
    if (rf_head_is(e, SYM_Plus)) {
        size_t n = e->data.function.arg_count;
        Expr** a = malloc(n * sizeof(Expr*));
        size_t k = 0;
        for (; k < n; k++) {
            /* a term with nothing to cancel is left as it is:
             * Cancel[x + 1/(1 - x^2)] keeps 1/(1 - x^2) */
            const Expr* t = e->data.function.args[k];
            long s0 = C.splits;
            a[k] = rf_has_den(t) ? rf_normal_form(t, &C) : expr_copy((Expr*)t);
            if (!a[k]) break;
            if (C.splits == s0 && rf_has_den(t)) { expr_free(a[k]); a[k] = expr_copy((Expr*)t); }
        }
        if (k == n) out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), a, n));
        else for (size_t i = 0; i < k; i++) expr_free(a[i]);
        free(a);
    } else {
        out = rf_normal_form(e, &C);
    }
    rf_ctx_clear(&C);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Apart                                                              */
/* ------------------------------------------------------------------ */

/* The variable of Apart is variable 0 of the context (rf_ctx_build puts it
 * first), so every factor is normalised to a positive leading coefficient in
 * it:  Apart[x^3/(x - a), x] -> a^2 + a x + x^2 + a^3/(-a + x).  The other
 * variables (parameters) are coefficients: the field is K = Q(params). */

/* An irreducible factor of the denominator in x, with its multiplicity. */
typedef struct { fmpz_mpoly_struct p; slong k; } RfApFac;

/* Expanded polynomial in x with rational coefficients (the classical form of
 * a partial-fraction numerator). */
static Expr* rf_fmpq_poly_expr(const fmpq_poly_t a, const RfCtx* C) {
    slong len = fmpq_poly_length(a);
    if (len == 0) return expr_new_integer(0);
    Expr** t = malloc((size_t)len * sizeof(Expr*));
    size_t n = 0;
    fmpq_t q; fmpq_init(q);
    for (slong i = 0; i < len; i++) {
        fmpq_poly_get_coeff_fmpq(q, a, i);
        if (fmpq_is_zero(q)) continue;
        Expr* c = fb_fmpq_to_expr(q);
        if (i == 0) { t[n++] = c; continue; }
        Expr* xp = i == 1 ? rf_atom_expr(C, 0)
            : expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){ rf_atom_expr(C, 0), expr_new_integer((int64_t)i) }, 2);
        t[n++] = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ c, xp }, 2);
    }
    fmpq_clear(q);
    Expr* r = n == 1 ? t[0] : expr_new_function(expr_new_symbol(SYM_Plus), t, n);
    free(t);
    return r;
}

/* The numerator a(x) over Q as an Rf: content times its primitive part. */
static int rf_from_fmpq_poly(Rf* R, const fmpq_poly_t a, RfCtx* C) {
    fmpz_poly_t z; fmpz_poly_init(z);
    fmpq_poly_get_numerator(z, a);
    fmpz_mpoly_t m; fmpz_mpoly_init(m, C->ctx);
    fmpz_mpoly_set_fmpz_poly(m, z, 0, C->ctx);
    fmpq_one(R->c);
    fmpq_div_fmpz(R->c, R->c, fmpq_poly_denref(a));
    int ok = rf_push(R, m, 1, C);
    fmpz_mpoly_clear(m, C->ctx);
    fmpz_poly_clear(z);
    return ok;
}

/* One term coef / P^j in the output form.  `coef` is consumed. */
static Expr* rf_apart_term(Rf* coef, const fmpz_mpoly_t P, slong j, const RfCtx* C) {
    Expr* den = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ rf_poly_to_expr(P, C), expr_new_integer(-(int64_t)j) }, 2);
    return rf_to_expr_raw(coef, C, den);
}

/* Partial fractions over Q (no parameters): fmpq_poly throughout. */
static Expr* rf_apart_q(const fmpz_mpoly_t Nm, const fmpq_t Kc, RfApFac* F, slong nf, RfCtx* C) {
    fmpz_poly_t z; fmpz_poly_init(z);
    fmpq_poly_t N, D, Qp, Rm, mi, g, s, t, Bi, cur, quot, rem;
    fmpq_poly_struct* P = malloc((size_t)nf * sizeof(fmpq_poly_struct));
    fmpq_poly_struct* Pk = malloc((size_t)nf * sizeof(fmpq_poly_struct));
    fmpq_poly_init(N); fmpq_poly_init(D); fmpq_poly_init(Qp); fmpq_poly_init(Rm);
    fmpq_poly_init(mi); fmpq_poly_init(g); fmpq_poly_init(s); fmpq_poly_init(t);
    fmpq_poly_init(Bi); fmpq_poly_init(cur); fmpq_poly_init(quot); fmpq_poly_init(rem);
    fmpz_mpoly_get_fmpz_poly(z, Nm, 0, C->ctx);
    fmpq_poly_set_fmpz_poly(N, z);
    fmpq_poly_scalar_mul_fmpq(N, N, Kc);
    fmpq_poly_one(D);
    for (slong i = 0; i < nf; i++) {
        fmpq_poly_init(P + i); fmpq_poly_init(Pk + i);
        fmpz_mpoly_get_fmpz_poly(z, &F[i].p, 0, C->ctx);
        fmpq_poly_set_fmpz_poly(P + i, z);
        fmpq_poly_pow(Pk + i, P + i, (ulong)F[i].k);
        fmpq_poly_mul(D, D, Pk + i);
    }
    fmpq_poly_divrem(Qp, Rm, N, D);

    size_t cap = 16, nt = 0;
    Expr** terms = malloc(cap * sizeof(Expr*));
    int factored = rf_factored_form(), ok = 1;
    if (!fmpq_poly_is_zero(Qp)) terms[nt++] = rf_fmpq_poly_expr(Qp, C);
    for (slong i = 0; i < nf && ok; i++) {
        /* B_i = R (m_i^-1 mod P_i^k) mod P_i^k,  m_i = D / P_i^k */
        fmpq_poly_div(mi, D, Pk + i);
        fmpq_poly_xgcd(g, s, t, Pk + i, mi);
        fmpq_poly_mul(Bi, Rm, t);
        fmpq_poly_rem(Bi, Bi, Pk + i);
        /* P_i-adic digits c_0 .. c_{k-1}: B_i/P^k = sum_j c_{k-j} / P^j */
        fmpq_poly_set(cur, Bi);
        slong k = F[i].k;
        for (slong d = 0; d < k && ok; d++) {
            fmpq_poly_divrem(quot, rem, cur, P + i);
            fmpq_poly_swap(cur, quot);
            if (fmpq_poly_is_zero(rem)) continue;
            Expr* term;
            if (factored) {
                Rf co; rf_init(&co);
                ok = rf_from_fmpq_poly(&co, rem, C);
                term = ok ? rf_apart_term(&co, &F[i].p, k - d, C) : NULL;
                rf_clear(&co, C);
            } else {
                Expr* den = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ rf_poly_to_expr(&F[i].p, C), expr_new_integer(-(int64_t)(k - d)) }, 2);
                term = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ rf_fmpq_poly_expr(rem, C), den }, 2);
            }
            if (!term) { ok = 0; break; }
            if (nt == cap) { cap *= 2; terms = realloc(terms, cap * sizeof(Expr*)); }
            terms[nt++] = term;
        }
    }
    Expr* out = NULL;
    if (ok) out = nt == 0 ? expr_new_integer(0)
                : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nt));
    else for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
    free(terms);
    for (slong i = 0; i < nf; i++) { fmpq_poly_clear(P + i); fmpq_poly_clear(Pk + i); }
    free(P); free(Pk);
    fmpq_poly_clear(N); fmpq_poly_clear(D); fmpq_poly_clear(Qp); fmpq_poly_clear(Rm);
    fmpq_poly_clear(mi); fmpq_poly_clear(g); fmpq_poly_clear(s); fmpq_poly_clear(t);
    fmpq_poly_clear(Bi); fmpq_poly_clear(cur); fmpq_poly_clear(quot); fmpq_poly_clear(rem);
    fmpz_poly_clear(z);
    return out;
}

/* Parameter polynomials live in `pctx` (variables 1.. of the main context
 * renumbered 0..). */
static void rf_param_to_main(fmpz_mpoly_t out, const fmpz_mpoly_t p, ulong xpow,
                             const fmpz_mpoly_ctx_struct* pctx, const RfCtx* C, ulong* ex) {
    fmpz_mpoly_zero(out, C->ctx);
    for (slong i = 0; i < p->length; i++) {
        fmpz_mpoly_get_term_exp_ui(ex + 1, p, i, pctx);
        ex[0] = xpow;
        fmpz_mpoly_push_term_fmpz_ui(out, p->coeffs + i, ex, C->ctx);
    }
    fmpz_mpoly_sort_terms(out, C->ctx);
    fmpz_mpoly_combine_like_terms(out, C->ctx);
}

/* fmpz_mpoly (x first) -> gr_poly in x over K = Q(params). */
static int rf_to_kpoly(gr_poly_t A, const fmpz_mpoly_t p, gr_ctx_t K,
                       const fmpz_mpoly_ctx_struct* pctx, const RfCtx* C, ulong* ex) {
    if (gr_poly_zero(A, K) != GR_SUCCESS) return 0;
    slong deg = fmpz_mpoly_degree_si(p, 0, C->ctx);
    if (deg < 0) return 1;
    fmpz_mpoly_struct* cf = malloc((size_t)(deg + 1) * sizeof(fmpz_mpoly_struct));
    if (!cf) return 0;
    for (slong d = 0; d <= deg; d++) fmpz_mpoly_init(cf + d, pctx);
    for (slong i = 0; i < p->length; i++) {
        fmpz_mpoly_get_term_exp_ui(ex, p, i, C->ctx);
        fmpz_mpoly_push_term_fmpz_ui(cf + ex[0], p->coeffs + i, ex + 1, pctx);
    }
    fmpz_mpoly_q_t q; fmpz_mpoly_q_init(q, pctx);
    int ok = 1;
    for (slong d = 0; d <= deg && ok; d++) {
        fmpz_mpoly_sort_terms(cf + d, pctx);
        fmpz_mpoly_combine_like_terms(cf + d, pctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_numref(q), cf + d, pctx);
        fmpz_mpoly_one(fmpz_mpoly_q_denref(q), pctx);
        ok = gr_poly_set_coeff_scalar(A, d, q, K) == GR_SUCCESS;
    }
    fmpz_mpoly_q_clear(q, pctx);
    for (slong d = 0; d <= deg; d++) fmpz_mpoly_clear(cf + d, pctx);
    free(cf);
    return ok;
}

/* R *= the parameter rational function q (numerator kept whole, the
 * denominator factored), times x^xpow. */
static int rf_mul_kelt(Rf* R, const fmpz_mpoly_q_t q, ulong xpow,
                       const fmpz_mpoly_ctx_struct* pctx, RfCtx* C, ulong* ex) {
    fmpz_mpoly_t m; fmpz_mpoly_init(m, C->ctx);
    rf_param_to_main(m, fmpz_mpoly_q_numref(q), xpow, pctx, C, ex);
    int ok = rf_push(R, m, 1, C);
    const fmpz_mpoly_struct* d = fmpz_mpoly_q_denref(q);
    if (ok && !fmpz_mpoly_is_one(d, pctx)) {
        fmpz_mpoly_factor_t f; fmpz_mpoly_factor_init(f, pctx);
        if (!fmpz_mpoly_factor(f, d, pctx)) ok = 0;
        else {
            fmpq_div_fmpz(R->c, R->c, f->constant);
            for (slong i = 0; i < f->num && ok; i++) {
                rf_param_to_main(m, f->poly + i, 0, pctx, C, ex);
                ok = rf_push(R, m, -fmpz_get_si(f->exp + i), C);
            }
        }
        fmpz_mpoly_factor_clear(f, pctx);
    }
    fmpz_mpoly_clear(m, C->ctx);
    return ok;
}

/* a(x) in K[x] -> Rf: n(x, params) / L(params), L the lcm of the coefficient
 * denominators (factored), n the numerator over it (expanded). */
static int rf_from_kpoly(Rf* R, const gr_poly_t a, gr_ctx_t K,
                         const fmpz_mpoly_ctx_struct* pctx, RfCtx* C, ulong* ex) {
    slong len = gr_poly_length(a, K);
    fmpz_mpoly_q_t c, s; fmpz_mpoly_q_init(c, pctx); fmpz_mpoly_q_init(s, pctx);
    fmpz_mpoly_t L, g, n, t, tm;
    fmpz_mpoly_init(L, pctx); fmpz_mpoly_init(g, pctx); fmpz_mpoly_init(t, pctx);
    fmpz_mpoly_init(n, C->ctx); fmpz_mpoly_init(tm, C->ctx);
    fmpz_mpoly_one(L, pctx);
    int ok = 1;
    for (slong i = 0; i < len && ok; i++) {
        ok = gr_poly_get_coeff_scalar(c, a, i, K) == GR_SUCCESS;
        if (!ok || fmpz_mpoly_q_is_zero(c, pctx)) continue;
        const fmpz_mpoly_struct* d = fmpz_mpoly_q_denref(c);
        if (!fmpz_mpoly_gcd(g, L, d, pctx)) { ok = 0; break; }
        if (!fmpz_mpoly_divides(t, d, g, pctx)) { ok = 0; break; }
        fmpz_mpoly_mul(L, L, t, pctx);
    }
    for (slong i = 0; i < len && ok; i++) {
        ok = gr_poly_get_coeff_scalar(c, a, i, K) == GR_SUCCESS;
        if (!ok || fmpz_mpoly_q_is_zero(c, pctx)) continue;
        if (!fmpz_mpoly_divides(t, L, fmpz_mpoly_q_denref(c), pctx)) { ok = 0; break; }
        fmpz_mpoly_mul(t, t, fmpz_mpoly_q_numref(c), pctx);
        rf_param_to_main(tm, t, (ulong)i, pctx, C, ex);
        fmpz_mpoly_add(n, n, tm, C->ctx);
    }
    if (ok) {
        fmpz_mpoly_q_one(s, pctx);
        fmpz_mpoly_swap(fmpz_mpoly_q_denref(s), L, pctx);
        ok = rf_push(R, n, 1, C) && rf_mul_kelt(R, s, 0, pctx, C, ex);
    }
    fmpz_mpoly_q_clear(c, pctx); fmpz_mpoly_q_clear(s, pctx);
    fmpz_mpoly_clear(L, pctx); fmpz_mpoly_clear(g, pctx); fmpz_mpoly_clear(t, pctx);
    fmpz_mpoly_clear(n, C->ctx); fmpz_mpoly_clear(tm, C->ctx);
    return ok;
}

/* Partial fractions over K = Q(params): gr_poly over fmpz_mpoly_q. */
static Expr* rf_apart_k(const fmpz_mpoly_t Nm, const Rf* Kc, RfApFac* F, slong nf, RfCtx* C) {
    slong np = (slong)C->n - 1;
    const fmpz_mpoly_ctx_struct* pctx = rf_ctx_for(np);
    if (!pctx) return NULL;
    ulong* ex = malloc((size_t)C->n * sizeof(ulong));
    if (!ex) return NULL;
    gr_ctx_t K;
    gr_ctx_init_fmpz_mpoly_q(K, np, ORD_LEX);
    gr_poly_t N, D, Qp, Rm, mi, g, s, t, Bi, cur, quot, rem;
    gr_poly_init(N, K); gr_poly_init(D, K); gr_poly_init(Qp, K); gr_poly_init(Rm, K);
    gr_poly_init(mi, K); gr_poly_init(g, K); gr_poly_init(s, K); gr_poly_init(t, K);
    gr_poly_init(Bi, K); gr_poly_init(cur, K); gr_poly_init(quot, K); gr_poly_init(rem, K);
    gr_poly_struct* P = malloc((size_t)nf * sizeof(gr_poly_struct));
    gr_poly_struct* Pk = malloc((size_t)nf * sizeof(gr_poly_struct));
    int ok = P && Pk && ex;
    slong ninit = 0;
    int st = GR_SUCCESS;
    ok = ok && rf_to_kpoly(N, Nm, K, pctx, C, ex);
    if (ok) st |= gr_poly_one(D, K);
    for (slong i = 0; i < nf && ok; i++) {
        gr_poly_init(P + i, K); gr_poly_init(Pk + i, K); ninit = i + 1;
        ok = rf_to_kpoly(P + i, &F[i].p, K, pctx, C, ex);
        if (ok) {
            st |= gr_poly_pow_ui(Pk + i, P + i, (ulong)F[i].k, K);
            st |= gr_poly_mul(D, D, Pk + i, K);
        }
    }
    if (ok) st |= gr_poly_divrem(Qp, Rm, N, D, K);
    ok = ok && st == GR_SUCCESS;

    size_t cap = 16, nt = 0;
    Expr** terms = malloc(cap * sizeof(Expr*));
    gr_ptr c;
    GR_TMP_INIT(c, K);
    /* the polynomial part, power by power: a^2 + a x + x^2 */
    for (slong r = 0; ok && r < gr_poly_length(Qp, K); r++) {
        ok = gr_poly_get_coeff_scalar(c, Qp, r, K) == GR_SUCCESS;
        if (!ok || fmpz_mpoly_q_is_zero((fmpz_mpoly_q_struct*)c, pctx)) continue;
        Rf co; rf_init(&co);
        fmpq_set(co.c, Kc->c);
        for (slong i = 0; i < Kc->n && ok; i++) {
            fmpz_mpoly_t m; fmpz_mpoly_init(m, C->ctx);
            fmpz_mpoly_set(m, &Kc->f[i].p, C->ctx);
            ok = rf_push(&co, m, Kc->f[i].e, C);
            fmpz_mpoly_clear(m, C->ctx);
        }
        ok = ok && rf_mul_kelt(&co, (fmpz_mpoly_q_struct*)c, (ulong)r, pctx, C, ex);
        if (ok) {
            if (nt == cap) { cap *= 2; terms = realloc(terms, cap * sizeof(Expr*)); }
            terms[nt++] = rf_to_expr_raw(&co, C, NULL);
        }
        rf_clear(&co, C);
    }
    for (slong i = 0; i < nf && ok; i++) {
        st = GR_SUCCESS;
        st |= gr_poly_div(mi, D, Pk + i, K);
        st |= RF_GR_POLY_XGCD(g, s, t, Pk + i, mi, K);
        st |= gr_poly_mul(Bi, Rm, t, K);
        st |= gr_poly_rem(Bi, Bi, Pk + i, K);
        st |= gr_poly_set(cur, Bi, K);
        if (st != GR_SUCCESS) { ok = 0; break; }
        slong k = F[i].k;
        for (slong d = 0; d < k && ok; d++) {
            if (gr_poly_divrem(quot, rem, cur, P + i, K) != GR_SUCCESS) { ok = 0; break; }
            gr_poly_swap(cur, quot, K);
            if (gr_poly_length(rem, K) == 0) continue;
            Rf co; rf_init(&co);
            fmpq_set(co.c, Kc->c);
            for (slong q = 0; q < Kc->n && ok; q++) {
                fmpz_mpoly_t m; fmpz_mpoly_init(m, C->ctx);
                fmpz_mpoly_set(m, &Kc->f[q].p, C->ctx);
                ok = rf_push(&co, m, Kc->f[q].e, C);
                fmpz_mpoly_clear(m, C->ctx);
            }
            Rf an; rf_init(&an);
            ok = ok && rf_from_kpoly(&an, rem, K, pctx, C, ex) && rf_mul_into(&co, &an, C);
            rf_clear(&an, C);
            if (ok) {
                if (nt == cap) { cap *= 2; terms = realloc(terms, cap * sizeof(Expr*)); }
                terms[nt++] = rf_apart_term(&co, &F[i].p, k - d, C);
            }
            rf_clear(&co, C);
        }
    }
    GR_TMP_CLEAR(c, K);
    Expr* out = NULL;
    if (ok) out = nt == 0 ? expr_new_integer(0)
                : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nt));
    else for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
    free(terms);
    for (slong i = 0; i < ninit; i++) { gr_poly_clear(P + i, K); gr_poly_clear(Pk + i, K); }
    free(P); free(Pk);
    gr_poly_clear(N, K); gr_poly_clear(D, K); gr_poly_clear(Qp, K); gr_poly_clear(Rm, K);
    gr_poly_clear(mi, K); gr_poly_clear(g, K); gr_poly_clear(s, K); gr_poly_clear(t, K);
    gr_poly_clear(Bi, K); gr_poly_clear(cur, K); gr_poly_clear(quot, K); gr_poly_clear(rem, K);
    gr_ctx_clear(K);
    free(ex);
    return out;
}

Expr* ratfun_apart(const Expr* e, const Expr* x) {
    if (!e || !x || x->type != EXPR_SYMBOL) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, x)) return NULL;
    if (C.gidx == 0 || C.atom[0]->type != EXPR_SYMBOL) { rf_ctx_clear(&C); return NULL; }
    Expr* out = NULL;
    Rf R; rf_init(&R);
    int ok = rf_from_expr(e, &R, &C) && (rf_head_is(e, SYM_Plus) || rf_single(&R, &C)) && !C.fail;
    /* Split: the numerator N (expanded), the x-free constant Kc, and the
     * irreducible x-factors of the denominator with their multiplicities. */
    fmpz_mpoly_t N; fmpz_mpoly_init(N, C.ctx);
    Rf Kc; rf_init(&Kc);
    RfApFac* F = NULL; slong nf = 0, cf = 0;
    if (ok) {
        fmpq_set(Kc.c, R.c);
        rf_expand_side(N, &R, 0, &C);
        fmpz_mpoly_factor_t fa; fmpz_mpoly_factor_init(fa, C.ctx);
        for (slong i = 0; i < R.n && ok; i++) {
            if (R.f[i].e > 0) continue;
            if (fmpz_mpoly_degree_si(&R.f[i].p, 0, C.ctx) == 0) {
                fmpz_mpoly_t m; fmpz_mpoly_init(m, C.ctx);
                fmpz_mpoly_set(m, &R.f[i].p, C.ctx);
                ok = rf_push(&Kc, m, R.f[i].e, &C);
                fmpz_mpoly_clear(m, C.ctx);
                continue;
            }
            if (!fmpz_mpoly_factor(fa, &R.f[i].p, C.ctx)) { ok = 0; break; }
            fmpz_t cz; fmpz_init_set(cz, fa->constant);
            ok = rf_mul_fmpz_pow(Kc.c, cz, R.f[i].e);
            fmpz_clear(cz);
            for (slong j = 0; j < fa->num && ok; j++) {
                slong m = fmpz_get_si(fa->exp + j) * -R.f[i].e;
                fmpz_mpoly_t p; fmpz_mpoly_init(p, C.ctx);
                fmpz_mpoly_set(p, fa->poly + j, C.ctx);
                if (fmpz_mpoly_degree_si(p, 0, C.ctx) == 0) {
                    ok = rf_push(&Kc, p, -m, &C);
                } else {
                    fmpz_t u; fmpz_init(u);
                    rf_normalise(p, u, C.ctx);
                    ok = rf_mul_fmpz_pow(Kc.c, u, -m);
                    fmpz_clear(u);
                    slong k = -1;
                    for (slong q = 0; q < nf; q++) if (fmpz_mpoly_equal(&F[q].p, p, C.ctx)) k = q;
                    if (k >= 0) F[k].k += m;
                    else if (ok) {
                        if (nf == cf) { cf = cf ? 2 * cf : 8; F = realloc(F, (size_t)cf * sizeof(RfApFac)); }
                        fmpz_mpoly_init(&F[nf].p, C.ctx);
                        fmpz_mpoly_swap(&F[nf].p, p, C.ctx);
                        F[nf].k = m;
                        nf++;
                    }
                }
                fmpz_mpoly_clear(p, C.ctx);
            }
        }
        fmpz_mpoly_factor_clear(fa, C.ctx);
    }
    if (ok && nf > 0) {
        if (C.n == 1) {
            out = rf_apart_q(N, Kc.c, F, nf, &C);
        } else {
            out = rf_apart_k(N, &Kc, F, nf, &C);
        }
    }
    for (slong i = 0; i < nf; i++) fmpz_mpoly_clear(&F[i].p, C.ctx);
    free(F);
    rf_clear(&Kc, &C);
    fmpz_mpoly_clear(N, C.ctx);
    rf_clear(&R, &C);
    rf_ctx_clear(&C);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Modulus -> p: Together, Cancel and (univariate) Apart over GF(p)    */
/* ------------------------------------------------------------------ */

/* A fraction over GF(p)[vars] in lowest terms with a monic denominator (the
 * subset of FLINT's fmpz_mod_mpoly_q this needs; that module is newer than
 * the FLINT 3.0 the build supports). */
typedef struct { fmpz_mod_mpoly_struct num, den; } rfmq_struct;
typedef rfmq_struct rfmq_t[1];
#define rfmq_numref(x) (&(x)->num)
#define rfmq_denref(x) (&(x)->den)

static void rfmq_init(rfmq_t x, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_init(&x->num, m); fmpz_mod_mpoly_init(&x->den, m);
    fmpz_mod_mpoly_one(&x->den, m);
}
static void rfmq_clear(rfmq_t x, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_clear(&x->num, m); fmpz_mod_mpoly_clear(&x->den, m);
}
static void rfmq_canonicalise(rfmq_t x, const fmpz_mod_mpoly_ctx_t m) {
    if (fmpz_mod_mpoly_is_zero(&x->num, m)) { fmpz_mod_mpoly_one(&x->den, m); return; }
    if (!fmpz_mod_mpoly_is_fmpz(&x->den, m)) {
        fmpz_mod_mpoly_t g; fmpz_mod_mpoly_init(g, m);
        if (fmpz_mod_mpoly_gcd(g, &x->num, &x->den, m) && !fmpz_mod_mpoly_is_fmpz(g, m)) {
            fmpz_mod_mpoly_divides(&x->num, &x->num, g, m);
            fmpz_mod_mpoly_divides(&x->den, &x->den, g, m);
        }
        fmpz_mod_mpoly_clear(g, m);
    }
    fmpz_t c; fmpz_init(c);
    fmpz_mod_mpoly_get_term_coeff_fmpz(c, &x->den, 0, m);       /* leading */
    if (!fmpz_is_one(c)) {
        fmpz_mod_inv(c, c, m->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz(&x->num, &x->num, c, m);
        fmpz_mod_mpoly_scalar_mul_fmpz(&x->den, &x->den, c, m);
    }
    fmpz_clear(c);
}
static void rfmq_set_si(rfmq_t x, slong v, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_set_si(&x->num, v, m); fmpz_mod_mpoly_one(&x->den, m);
}
static void rfmq_set_fmpz(rfmq_t x, const fmpz_t v, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_set_fmpz(&x->num, v, m); fmpz_mod_mpoly_one(&x->den, m);
}
static void rfmq_set_fmpq(rfmq_t x, const fmpq_t v, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_set_fmpz(&x->num, fmpq_numref(v), m);
    fmpz_mod_mpoly_set_fmpz(&x->den, fmpq_denref(v), m);
    rfmq_canonicalise(x, m);
}
static void rfmq_add(rfmq_t r, const rfmq_t a, const rfmq_t b, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_t n, d, t;
    fmpz_mod_mpoly_init(n, m); fmpz_mod_mpoly_init(d, m); fmpz_mod_mpoly_init(t, m);
    if (fmpz_mod_mpoly_equal(&a->den, &b->den, m)) {
        fmpz_mod_mpoly_add(n, &a->num, &b->num, m);
        fmpz_mod_mpoly_set(d, &a->den, m);
    } else {
        fmpz_mod_mpoly_mul(n, &a->num, &b->den, m);
        fmpz_mod_mpoly_mul(t, &b->num, &a->den, m);
        fmpz_mod_mpoly_add(n, n, t, m);
        fmpz_mod_mpoly_mul(d, &a->den, &b->den, m);
    }
    fmpz_mod_mpoly_swap(&r->num, n, m); fmpz_mod_mpoly_swap(&r->den, d, m);
    fmpz_mod_mpoly_clear(n, m); fmpz_mod_mpoly_clear(d, m); fmpz_mod_mpoly_clear(t, m);
    rfmq_canonicalise(r, m);
}
static void rfmq_mul(rfmq_t r, const rfmq_t a, const rfmq_t b, const fmpz_mod_mpoly_ctx_t m) {
    fmpz_mod_mpoly_t n, d;
    fmpz_mod_mpoly_init(n, m); fmpz_mod_mpoly_init(d, m);
    fmpz_mod_mpoly_mul(n, &a->num, &b->num, m);
    fmpz_mod_mpoly_mul(d, &a->den, &b->den, m);
    fmpz_mod_mpoly_swap(&r->num, n, m); fmpz_mod_mpoly_swap(&r->den, d, m);
    fmpz_mod_mpoly_clear(n, m); fmpz_mod_mpoly_clear(d, m);
    rfmq_canonicalise(r, m);
}
static void rfmq_inv(rfmq_t r, const rfmq_t a, const fmpz_mod_mpoly_ctx_t m) {
    if (r != a) { fmpz_mod_mpoly_set(&r->num, &a->num, m); fmpz_mod_mpoly_set(&r->den, &a->den, m); }
    fmpz_mod_mpoly_swap(&r->num, &r->den, m);
    rfmq_canonicalise(r, m);
}

typedef struct {
    RfCtx* C;
    fmpz_mod_mpoly_ctx_t m;
    fmpz_mod_ctx_t z;
    fmpz_t p;
} RfMod;

static void rfm_init(RfMod* M, RfCtx* C, const fmpz_t p) {
    M->C = C;
    fmpz_init_set(M->p, p);
    fmpz_mod_mpoly_ctx_init(M->m, (slong)C->n, ORD_LEX, p);
    fmpz_mod_ctx_init(M->z, p);
}

static void rfm_clear(RfMod* M) {
    fmpz_mod_mpoly_ctx_clear(M->m);
    fmpz_mod_ctx_clear(M->z);
    fmpz_clear(M->p);
}

/* A rational number modulo p; 0 when its denominator vanishes there. */
static int rfm_number(rfmq_t out, const fmpq_t q, RfMod* M) {
    fmpz_t d; fmpz_init(d);
    fmpz_mod(d, fmpq_denref(q), M->p);
    int ok = !fmpz_is_zero(d);
    if (ok) rfmq_set_fmpq(out, q, M->m);
    fmpz_clear(d);
    return ok;
}

static int rfm_gen_pow(rfmq_t out, slong v, slong k, RfMod* M) {
    fmpz_mod_mpoly_t g; fmpz_mod_mpoly_init(g, M->m);
    fmpz_mod_mpoly_gen(g, v, M->m);
    fmpz_mod_mpoly_pow_ui(rfmq_numref(out), g, (ulong)(k < 0 ? -k : k), M->m);
    fmpz_mod_mpoly_one(rfmq_denref(out), M->m);
    if (k < 0) rfmq_inv(out, out, M->m);
    fmpz_mod_mpoly_clear(g, M->m);
    return 1;
}

static int rfm_from_expr(const Expr* e, rfmq_t out, RfMod* M);

/* A sum: the polynomial terms are pushed into one accumulator (one sort and
 * combine at the end, not a merge per term); the fractional ones are added
 * as fractions. */
static int rfm_from_plus(const Expr* e, rfmq_t out, RfMod* M) {
    size_t n = e->data.function.arg_count;
    slong nv = (slong)M->C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    if (!ex) return 0;
    fmpz_mod_mpoly_t acc; fmpz_mod_mpoly_init(acc, M->m);
    rfmq_t t; rfmq_init(t, M->m);
    fmpz_t c; fmpz_init(c);
    rfmq_set_si(out, 0, M->m);
    int ok = 1, frac = 0;
    for (size_t i = 0; i < n && ok; i++) {
        ok = rfm_from_expr(e->data.function.args[i], t, M);
        if (!ok) break;
        if (fmpz_mod_mpoly_is_one(rfmq_denref(t), M->m)) {
            const fmpz_mod_mpoly_struct* q = rfmq_numref(t);
            for (slong k = 0; k < fmpz_mod_mpoly_length(q, M->m); k++) {
                fmpz_mod_mpoly_get_term_coeff_fmpz(c, q, k, M->m);
                fmpz_mod_mpoly_get_term_exp_ui(ex, q, k, M->m);
                fmpz_mod_mpoly_push_term_fmpz_ui(acc, c, ex, M->m);
            }
        } else {
            rfmq_add(out, out, t, M->m);
            frac = 1;
        }
    }
    if (ok) {
        fmpz_mod_mpoly_sort_terms(acc, M->m);
        fmpz_mod_mpoly_combine_like_terms(acc, M->m);
        if (!frac) {
            fmpz_mod_mpoly_swap(rfmq_numref(out), acc, M->m);
            fmpz_mod_mpoly_one(rfmq_denref(out), M->m);
        } else {
            fmpz_mod_mpoly_swap(rfmq_numref(t), acc, M->m);
            fmpz_mod_mpoly_one(rfmq_denref(t), M->m);
            rfmq_add(out, out, t, M->m);
        }
    }
    fmpz_clear(c);
    rfmq_clear(t, M->m);
    fmpz_mod_mpoly_clear(acc, M->m);
    free(ex);
    return ok;
}

/* Expr -> reduced fraction over GF(p); 0 declines. */
static int rfm_from_expr(const Expr* e, rfmq_t out, RfMod* M) {
    RfCtx* C = M->C;
    switch (e->type) {
        case EXPR_INTEGER:
            rfmq_set_si(out, e->data.integer, M->m);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            rfmq_set_fmpz(out, z, M->m);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_SYMBOL: {
            slong v = rf_atom_index(C, e);
            return v >= 0 && rfm_gen_pow(out, v, 1, M);
        }
        case EXPR_FUNCTION: {
            const char* h = e->data.function.head->type == EXPR_SYMBOL
                          ? e->data.function.head->data.symbol.name : NULL;
            size_t n = e->data.function.arg_count;
            if (h == SYM_Rational) {
                fmpq_t q; fmpq_init(q);
                int ok = fb_expr_to_fmpq(e, q) && rfm_number(out, q, M);
                fmpq_clear(q);
                return ok;
            }
            if (h == SYM_Plus) return rfm_from_plus(e, out, M);
            if (h == SYM_Times) {
                rfmq_t t; rfmq_init(t, M->m);
                rfmq_set_si(out, 1, M->m);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = rfm_from_expr(e->data.function.args[i], t, M);
                    if (!ok) break;
                    rfmq_mul(out, out, t, M->m);
                }
                rfmq_clear(t, M->m);
                return ok;
            }
            if (h == SYM_Power && n == 2) {
                const Expr* b = e->data.function.args[0];
                const Expr* x = e->data.function.args[1];
                if (x->type == EXPR_INTEGER) {
                    rfmq_t base; rfmq_init(base, M->m);
                    int ok = rfm_from_expr(b, base, M);
                    int64_t k = x->data.integer;
                    if (ok && k < 0) {
                        if (fmpz_mod_mpoly_is_zero(rfmq_numref(base), M->m)) ok = 0;
                        else rfmq_inv(base, base, M->m);
                        k = -k;
                    }
                    if (ok) {
                        /* binary powering of the (reduced) fraction */
                        rfmq_set_si(out, 1, M->m);
                        while (k) {
                            if (k & 1) rfmq_mul(out, out, base, M->m);
                            k >>= 1;
                            if (k) rfmq_mul(base, base, base, M->m);
                        }
                    }
                    rfmq_clear(base, M->m);
                    return ok;
                }
                for (size_t i = 0; i < C->xn; i++)
                    if (expr_eq(C->xk[i], e)) return rfm_gen_pow(out, C->gidx, (slong)C->xpow[i], M);
                return 0;
            }
            slong v = rf_atom_index(C, e);
            return v >= 0 && rfm_gen_pow(out, v, 1, M);
        }
        default:
            return 0;
    }
}

static Expr* rfm_poly_to_expr(const fmpz_mod_mpoly_t P, RfMod* M) {
    RfCtx* C = M->C;
    slong len = fmpz_mod_mpoly_length(P, M->m);
    if (len == 0) return expr_new_integer(0);
    slong nv = (slong)C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    Expr** terms = malloc((size_t)len * sizeof(Expr*));
    Expr** fs = malloc((size_t)(nv + 1) * sizeof(Expr*));
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i < len; i++) {
        fmpz_mod_mpoly_get_term_coeff_fmpz(c, P, i, M->m);
        fmpz_mod_mpoly_get_term_exp_ui(ex, P, i, M->m);
        size_t nf = 0;
        if (!fmpz_is_one(c)) fs[nf++] = rf_fmpz_to_expr(c);
        for (slong v = 0; v < nv; v++) {
            if (!ex[v]) continue;
            Expr* a = rf_atom_expr(C, v);
            fs[nf++] = ex[v] == 1 ? a
                : expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ a, expr_new_integer((int64_t)ex[v]) }, 2);
        }
        terms[i] = nf == 0 ? expr_new_integer(1) : nf == 1 ? fs[0]
                 : expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
    }
    fmpz_clear(c);
    Expr* r = len == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)len);
    free(ex); free(terms); free(fs);
    return r;
}

/* num/den (reduced, den monic): the denominator as its irreducible factors in
 * the factored form, expanded otherwise. */
static Expr* rfm_fraction_expr(rfmq_t q, RfMod* M) {
    Expr* num = rfm_poly_to_expr(rfmq_numref(q), M);
    const fmpz_mod_mpoly_struct* d = rfmq_denref(q);
    if (fmpz_mod_mpoly_is_one(d, M->m)) return eval_and_free(num);
    Expr* den = NULL;
    if (rf_factored_form()) {
        fmpz_mod_mpoly_factor_t f; fmpz_mod_mpoly_factor_init(f, M->m);
        /* fully factored: Together[1/(x-1) + 1/(x+1), Modulus -> 3] has the
         * denominator (1 + x)(2 + x) */
        if (fmpz_mod_mpoly_factor(f, d, M->m) && fmpz_is_one(f->constant)) {
            Expr** a = malloc((size_t)(f->num ? f->num : 1) * sizeof(Expr*));
            for (slong i = 0; i < f->num; i++) {
                Expr* b = rfm_poly_to_expr(f->poly + i, M);
                a[i] = fmpz_is_one(f->exp + i) ? b
                    : expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ b, rf_fmpz_to_expr(f->exp + i) }, 2);
            }
            den = f->num == 1 ? a[0] : expr_new_function(expr_new_symbol(SYM_Times), a, (size_t)f->num);
            free(a);
        }
        fmpz_mod_mpoly_factor_clear(f, M->m);
    }
    if (!den) den = rfm_poly_to_expr(d, M);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ num, expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2));
}

/* One variable: lazy fractions of dense fmpz_mod_poly (no gcd until the end,
 * then one), far cheaper than the multivariate canonicalisation per step. */
typedef struct { fmpz_mod_poly_t n, d; } RfmU;

static void rfmu_init(RfmU* a, RfMod* M) { fmpz_mod_poly_init(a->n, M->z); fmpz_mod_poly_init(a->d, M->z); }
static void rfmu_clear(RfmU* a, RfMod* M) { fmpz_mod_poly_clear(a->n, M->z); fmpz_mod_poly_clear(a->d, M->z); }

/* e = c x^k with c an integer (the term shape of an expanded polynomial). */
static int rfmu_monomial(const Expr* e, const RfCtx* C, fmpz_t c, slong* k) {
    const Expr* v = e;
    fmpz_one(c);
    if (rf_head_is(e, SYM_Times) && e->data.function.arg_count == 2) {
        const Expr* a = e->data.function.args[0];
        if (a->type == EXPR_INTEGER) fmpz_set_si(c, a->data.integer);
        else if (a->type == EXPR_BIGINT) fmpz_set_mpz(c, a->data.bigint);
        else return 0;
        v = e->data.function.args[1];
    } else if (e->type == EXPR_INTEGER) { fmpz_set_si(c, e->data.integer); *k = 0; return 1; }
    else if (e->type == EXPR_BIGINT) { fmpz_set_mpz(c, e->data.bigint); *k = 0; return 1; }
    if (v->type == EXPR_SYMBOL && rf_atom_index(C, v) == 0) { *k = 1; return 1; }
    if (rf_head_is(v, SYM_Power) && v->data.function.arg_count == 2 &&
        v->data.function.args[0]->type == EXPR_SYMBOL &&
        rf_atom_index(C, v->data.function.args[0]) == 0 &&
        v->data.function.args[1]->type == EXPR_INTEGER && v->data.function.args[1]->data.integer > 0) {
        *k = (slong)v->data.function.args[1]->data.integer;
        return 1;
    }
    return 0;
}

static int rfmu_from_expr(const Expr* e, RfmU* out, RfMod* M) {
    const fmpz_mod_ctx_struct* z = M->z;
    RfCtx* C = M->C;
    fmpz_mod_poly_one(out->d, z);
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: {
            fmpz_t c; fmpz_init(c);
            if (e->type == EXPR_INTEGER) fmpz_set_si(c, e->data.integer); else fmpz_set_mpz(c, e->data.bigint);
            fmpz_mod_poly_set_fmpz(out->n, c, z);
            fmpz_clear(c);
            return 1;
        }
        case EXPR_SYMBOL:
            if (rf_atom_index(C, e) != 0) return 0;
            fmpz_mod_poly_zero(out->n, z);
            fmpz_mod_poly_set_coeff_ui(out->n, 1, 1, z);
            return 1;
        case EXPR_FUNCTION: {
            const char* h = e->data.function.head->type == EXPR_SYMBOL
                          ? e->data.function.head->data.symbol.name : NULL;
            size_t n = e->data.function.arg_count;
            if (h == SYM_Rational) {
                fmpq_t q; fmpq_init(q);
                int ok = fb_expr_to_fmpq(e, q);
                fmpz_t dd; fmpz_init(dd);
                if (ok) { fmpz_mod(dd, fmpq_denref(q), M->p); ok = !fmpz_is_zero(dd); }
                if (ok) {
                    fmpz_mod_poly_set_fmpz(out->n, fmpq_numref(q), z);
                    fmpz_mod_poly_set_fmpz(out->d, fmpq_denref(q), z);
                }
                fmpz_clear(dd); fmpq_clear(q);
                return ok;
            }
            if (h == SYM_Plus || h == SYM_Times) {
                RfmU t; rfmu_init(&t, M);
                fmpz_mod_poly_t u; fmpz_mod_poly_init(u, z);
                if (h == SYM_Plus) fmpz_mod_poly_zero(out->n, z); else fmpz_mod_poly_one(out->n, z);
                int ok = 1;
                fmpz_t mc, cur; fmpz_init(mc); fmpz_init(cur);
                for (size_t i = 0; i < n && ok; i++) {
                    slong mk;
                    if (h == SYM_Plus && fmpz_mod_poly_is_one(out->d, z) &&
                        rfmu_monomial(e->data.function.args[i], C, mc, &mk)) {
                        /* c x^k straight into its coefficient */
                        fmpz_mod_poly_get_coeff_fmpz(cur, out->n, mk, z);
                        fmpz_mod_set_fmpz(mc, mc, z);
                        fmpz_mod_add(cur, cur, mc, z);
                        fmpz_mod_poly_set_coeff_fmpz(out->n, mk, cur, z);
                        continue;
                    }
                    ok = rfmu_from_expr(e->data.function.args[i], &t, M);
                    if (!ok) break;
                    if (h == SYM_Times) {
                        if (fmpz_mod_poly_length(t.n, z) == 1 && fmpz_mod_poly_is_one(t.d, z))
                            fmpz_mod_poly_scalar_mul_fmpz(out->n, out->n, t.n->coeffs, z);
                        else {
                            fmpz_mod_poly_mul(out->n, out->n, t.n, z);
                            fmpz_mod_poly_mul(out->d, out->d, t.d, z);
                        }
                    } else if (fmpz_mod_poly_equal(out->d, t.d, z)) {
                        fmpz_mod_poly_add(out->n, out->n, t.n, z);
                    } else {
                        fmpz_mod_poly_mul(u, t.n, out->d, z);
                        fmpz_mod_poly_mul(out->n, out->n, t.d, z);
                        fmpz_mod_poly_add(out->n, out->n, u, z);
                        fmpz_mod_poly_mul(out->d, out->d, t.d, z);
                    }
                }
                fmpz_clear(mc); fmpz_clear(cur);
                fmpz_mod_poly_clear(u, z);
                rfmu_clear(&t, M);
                return ok;
            }
            if (h == SYM_Power && n == 2 && e->data.function.args[1]->type == EXPR_INTEGER) {
                int64_t k = e->data.function.args[1]->data.integer;
                const Expr* b = e->data.function.args[0];
                if (b->type == EXPR_SYMBOL && k > 0 && rf_atom_index(C, b) == 0) {
                    fmpz_mod_poly_zero(out->n, z);              /* x^k directly */
                    fmpz_mod_poly_set_coeff_ui(out->n, (slong)k, 1, z);
                    return 1;
                }
                int ok = rfmu_from_expr(b, out, M);
                if (!ok) return 0;
                if (k < 0) {
                    if (fmpz_mod_poly_is_zero(out->n, z)) return 0;
                    fmpz_mod_poly_swap(out->n, out->d, z);
                    k = -k;
                }
                fmpz_mod_poly_pow(out->n, out->n, (ulong)k, z);
                fmpz_mod_poly_pow(out->d, out->d, (ulong)k, z);
                return 1;
            }
            return 0;
        }
        default:
            return 0;
    }
}

static Expr* rfm_normal_form(const Expr* e, RfMod* M) {
    if (M->C->n == 1 && M->C->gidx < 0) {
        RfmU u; rfmu_init(&u, M);
        Expr* out = NULL;
        if (rfmu_from_expr(e, &u, M) && !fmpz_mod_poly_is_zero(u.d, M->z)) {
            /* one gcd, then into the canonical (monic denominator) form */
            fmpz_mod_poly_t g; fmpz_mod_poly_init(g, M->z);
            fmpz_mod_poly_gcd(g, u.n, u.d, M->z);
            if (fmpz_mod_poly_degree(g, M->z) > 0) {
                fmpz_mod_poly_div(u.n, u.n, g, M->z);
                fmpz_mod_poly_div(u.d, u.d, g, M->z);
            }
            fmpz_mod_poly_clear(g, M->z);
            rfmq_t q; rfmq_init(q, M->m);
            fmpz_mod_mpoly_set_fmpz_mod_poly(rfmq_numref(q), u.n, 0, M->m);
            fmpz_mod_mpoly_set_fmpz_mod_poly(rfmq_denref(q), u.d, 0, M->m);
            rfmq_canonicalise(q, M->m);
            out = rfm_fraction_expr(q, M);
            rfmq_clear(q, M->m);
        }
        rfmu_clear(&u, M);
        return out;
    }
    rfmq_t q; rfmq_init(q, M->m);
    Expr* out = rfm_from_expr(e, q, M) ? rfm_fraction_expr(q, M) : NULL;
    rfmq_clear(q, M->m);
    return out;
}

Expr* ratfun_together_mod(const Expr* e, const mpz_t p) {
    if (!e) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, NULL)) return NULL;
    fmpz_t pz; fmpz_init(pz); fmpz_set_mpz(pz, p);
    RfMod M; rfm_init(&M, &C, pz);
    Expr* out = rfm_normal_form(e, &M);
    rfm_clear(&M); fmpz_clear(pz);
    rf_ctx_clear(&C);
    return out;
}

Expr* ratfun_cancel_mod(const Expr* e, const mpz_t p) {
    if (!e) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, NULL)) return NULL;
    fmpz_t pz; fmpz_init(pz); fmpz_set_mpz(pz, p);
    RfMod M; rfm_init(&M, &C, pz);
    Expr* out = NULL;
    if (rf_head_is(e, SYM_Plus)) {
        size_t n = e->data.function.arg_count, k = 0;
        Expr** a = malloc(n * sizeof(Expr*));
        for (; k < n; k++) if (!(a[k] = rfm_normal_form(e->data.function.args[k], &M))) break;
        if (k == n) out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), a, n));
        else for (size_t i = 0; i < k; i++) expr_free(a[i]);
        free(a);
    } else out = rfm_normal_form(e, &M);
    rfm_clear(&M); fmpz_clear(pz);
    rf_ctx_clear(&C);
    return out;
}

static Expr* rfm_upoly_expr(const fmpz_mod_poly_t a, RfMod* M) {
    fmpz_mod_mpoly_t m; fmpz_mod_mpoly_init(m, M->m);
    fmpz_mod_mpoly_set_fmpz_mod_poly(m, a, 0, M->m);
    Expr* r = rfm_poly_to_expr(m, M);
    fmpz_mod_mpoly_clear(m, M->m);
    return r;
}

/* Apart over GF(p), one variable: factor the denominator, CRT, P-adic digits. */
Expr* ratfun_apart_mod(const Expr* e, const Expr* x, const mpz_t p) {
    if (!e || !x || x->type != EXPR_SYMBOL) return NULL;
    RfCtx C;
    if (!rf_ctx_build(e, &C, x)) return NULL;
    if (C.n != 1 || C.gidx == 0) { rf_ctx_clear(&C); return NULL; }
    fmpz_t pz; fmpz_init(pz); fmpz_set_mpz(pz, p);
    RfMod M; rfm_init(&M, &C, pz);
    const fmpz_mod_ctx_struct* z = M.z;
    Expr* out = NULL;
    rfmq_t q; rfmq_init(q, M.m);
    fmpz_mod_poly_t N, D, Qp, R, Pk, m, g, s, t, B, cur, quo, rem;
    fmpz_mod_poly_init(N, z); fmpz_mod_poly_init(D, z); fmpz_mod_poly_init(Qp, z);
    fmpz_mod_poly_init(R, z); fmpz_mod_poly_init(Pk, z); fmpz_mod_poly_init(m, z);
    fmpz_mod_poly_init(g, z); fmpz_mod_poly_init(s, z); fmpz_mod_poly_init(t, z);
    fmpz_mod_poly_init(B, z); fmpz_mod_poly_init(cur, z); fmpz_mod_poly_init(quo, z);
    fmpz_mod_poly_init(rem, z);
    fmpz_mod_poly_factor_t fa; fmpz_mod_poly_factor_init(fa, z);
    if (rfm_from_expr(e, q, &M) &&
        fmpz_mod_mpoly_get_fmpz_mod_poly(N, rfmq_numref(q), 0, M.m) &&
        fmpz_mod_mpoly_get_fmpz_mod_poly(D, rfmq_denref(q), 0, M.m) &&
        fmpz_mod_poly_degree(D, z) > 0) {
        /* D is monic (canonical); R/D with deg R < deg D */
        fmpz_mod_poly_divrem(Qp, R, N, D, z);
        fmpz_mod_poly_factor(fa, D, z);
        size_t cap = 16, nt = 0;
        Expr** terms = malloc(cap * sizeof(Expr*));
        if (!fmpz_mod_poly_is_zero(Qp, z)) terms[nt++] = rfm_upoly_expr(Qp, &M);
        for (slong i = 0; i < fa->num; i++) {
            slong k = fa->exp[i];
            fmpz_mod_poly_pow(Pk, fa->poly + i, (ulong)k, z);
            fmpz_mod_poly_div(m, D, Pk, z);
            fmpz_mod_poly_xgcd(g, s, t, Pk, m, z);       /* s Pk + t m = 1 */
            fmpz_mod_poly_mulmod(B, R, t, Pk, z);
            fmpz_mod_poly_set(cur, B, z);
            for (slong d = 0; d < k; d++) {
                fmpz_mod_poly_divrem(quo, rem, cur, fa->poly + i, z);
                fmpz_mod_poly_swap(cur, quo, z);
                if (fmpz_mod_poly_is_zero(rem, z)) continue;
                Expr* den = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ rfm_upoly_expr(fa->poly + i, &M), expr_new_integer(-(int64_t)(k - d)) }, 2);
                if (nt == cap) { cap *= 2; terms = realloc(terms, cap * sizeof(Expr*)); }
                terms[nt++] = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ rfm_upoly_expr(rem, &M), den }, 2);
            }
        }
        out = nt == 0 ? expr_new_integer(0)
            : eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, nt));
        free(terms);
    }
    fmpz_mod_poly_factor_clear(fa, z);
    fmpz_mod_poly_clear(N, z); fmpz_mod_poly_clear(D, z); fmpz_mod_poly_clear(Qp, z);
    fmpz_mod_poly_clear(R, z); fmpz_mod_poly_clear(Pk, z); fmpz_mod_poly_clear(m, z);
    fmpz_mod_poly_clear(g, z); fmpz_mod_poly_clear(s, z); fmpz_mod_poly_clear(t, z);
    fmpz_mod_poly_clear(B, z); fmpz_mod_poly_clear(cur, z); fmpz_mod_poly_clear(quo, z);
    fmpz_mod_poly_clear(rem, z);
    rfmq_clear(q, M.m);
    rfm_clear(&M); fmpz_clear(pz);
    rf_ctx_clear(&C);
    return out;
}

/* ------------------------------------------------------------------ */
/*  The Modulus option                                                 */
/* ------------------------------------------------------------------ */

int ratfun_modulus_option(const Expr* res, const char* head, mpz_t p, Expr** stripped) {
    *stripped = NULL;
    if (res->type != EXPR_FUNCTION) return 0;
    size_t n = res->data.function.arg_count;
    long at = -1;            /* the Modulus rule */
    long tf = -1;            /* a Trig -> False rule (the default: dropped) */
    for (size_t i = 1; i < n; i++) {
        const Expr* a = res->data.function.args[i];
        if (!(rf_head_is(a, SYM_Rule) && a->data.function.arg_count == 2 &&
              a->data.function.args[0]->type == EXPR_SYMBOL)) continue;
        const char* k = a->data.function.args[0]->data.symbol.name;
        const Expr* v = a->data.function.args[1];
        if (k == SYM_Modulus && at < 0) at = (long)i;
        if (k == SYM_Trig && tf < 0 && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False)
            tf = (long)i;
    }
    int drop_mod = 0;
    if (at >= 0) {
        const Expr* v = res->data.function.args[at]->data.function.args[1];
        if (expr_is_integer_like(v)) {
            expr_to_mpz(v, p);
            mpz_abs(p, p);
            if (mpz_sgn(p) != 0 && mpz_probab_prime_p(p, 30) == 0) {
                char* s = expr_to_string((Expr*)v);
                mth_message(head, "modp", "Value of option Modulus -> %s should be a prime number or zero.",
                            s ? s : "?");
                free(s);
                return -1;
            }
            if (mpz_sgn(p) != 0) {
                for (size_t i = 1; i < n; i++) {
                    const Expr* a = res->data.function.args[i];
                    if ((long)i == at || !rf_head_is(a, SYM_Rule) || a->data.function.arg_count != 2) continue;
                    const Expr* key = a->data.function.args[0];
                    const Expr* val = a->data.function.args[1];
                    if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_Extension &&
                        !(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_None)) {
                        mth_message(head, "mext", "Conflicting option setting. If Modulus is not 0, then "
                                    "Extension should be None.");
                        return -1;
                    }
                }
                return 1;
            }
            drop_mod = 1;                   /* Modulus -> 0 is the default */
        }
    }
    if (!drop_mod && tf < 0) return 0;
    Expr** a = malloc(n * sizeof(Expr*));
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
        if (!((long)i == at && drop_mod) && (long)i != tf) a[k++] = expr_copy(res->data.function.args[i]);
    *stripped = expr_new_function(expr_copy(res->data.function.head), a, k);
    free(a);
    return 0;
}

#else  /* !USE_FLINT */

int ratfun_modulus_option(const Expr* res, const char* head, mpz_t p, Expr** stripped) {
    (void)res; (void)head; (void)p; *stripped = NULL; return 0;
}
Expr* ratfun_together_mod(const Expr* e, const mpz_t p) { (void)e; (void)p; return NULL; }
Expr* ratfun_cancel_mod(const Expr* e, const mpz_t p) { (void)e; (void)p; return NULL; }
Expr* ratfun_apart_mod(const Expr* e, const Expr* x, const mpz_t p) { (void)e; (void)x; (void)p; return NULL; }

Expr* ratfun_together(const Expr* e) { (void)e; return NULL; }
Expr* ratfun_cancel(const Expr* e) { (void)e; return NULL; }
Expr* ratfun_apart(const Expr* e, const Expr* x) { (void)e; (void)x; return NULL; }

#endif /* USE_FLINT */
