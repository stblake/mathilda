/*
 * flint_bridge.c
 * --------------
 * Expr <-> FLINT boundary for the algebraic-extension arithmetic engine.
 * See flint_bridge.h and ALGEBRAIC_EXTENSION_ARITHMETIC_PLAN.md.
 *
 * M1: rational multivariate case R = Q[x_1..x_n] only. When USE_FLINT is
 * undefined the whole module degrades to stubs (available() -> 0, ops -> NULL).
 */
#include "flint_bridge.h"

#include <stdlib.h>
#include <string.h>

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "poly.h"
#include "sym_names.h"

#ifdef USE_FLINT

#include <math.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_mat.h>
#include <flint/mpoly.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/gr_vec.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

/* ------------------------------------------------------------------ */
/*  Variable set: the free symbols that become fmpq_mpoly generators.  */
/* ------------------------------------------------------------------ */

typedef struct {
    char** names;
    size_t count;
    size_t cap;
} VarSet;

static int varset_add(VarSet* vs, const char* name) {
    for (size_t i = 0; i < vs->count; i++)
        if (strcmp(vs->names[i], name) == 0) return 1;
    if (vs->count == vs->cap) {
        size_t nc = vs->cap ? vs->cap * 2 : 8;
        char** np = realloc(vs->names, nc * sizeof(char*));
        if (!np) return 0;
        vs->names = np;
        vs->cap = nc;
    }
    vs->names[vs->count] = mathilda_strdup(name);
    if (!vs->names[vs->count]) return 0;
    vs->count++;
    return 1;
}

static void varset_free(VarSet* vs) {
    for (size_t i = 0; i < vs->count; i++) free(vs->names[i]);
    free(vs->names);
}

static int cmp_str(const void* a, const void* b) {
    return strcmp(*(char* const*)a, *(char* const*)b);
}

static int var_index(const VarSet* vs, const char* name) {
    for (size_t i = 0; i < vs->count; i++)
        if (strcmp(vs->names[i], name) == 0) return (int)i;
    return -1;
}

/* Head symbol name of a function Expr, or NULL. */
static const char* fn_head_name(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) return h->data.symbol;
    return NULL;
}

/*
 * Recognise `e` as a member of Q[x_1..x_n] and collect its variables.
 * Returns 0 (caller bails to the classical path) for anything outside that
 * ring: inexact reals, negative/symbolic powers, unrecognised heads
 * (Sqrt[..], f[..] — those are extension generators, M2/M3 not M1).
 */
static int collect_vars(const Expr* e, VarSet* vs) {
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return 1;
        case EXPR_SYMBOL:
            return varset_add(vs, e->data.symbol);
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) return 0;
            if (strcmp(h, "Rational") == 0) return 1; /* numeric coefficient */
            if (strcmp(h, "Power") == 0) {
                if (e->data.function.arg_count != 2) return 0;
                const Expr* base = e->data.function.args[0];
                const Expr* exp  = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return 0;
                return collect_vars(base, vs);
            }
            if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    if (!collect_vars(e->data.function.args[i], vs)) return 0;
                return 1;
            }
            return 0;
        }
        default:
            return 0; /* EXPR_REAL / EXPR_MPFR / EXPR_STRING: not exact over Q */
    }
}

/* ------------------------------------------------------------------ */
/*  Expr -> fmpq_mpoly                                                 */
/* ------------------------------------------------------------------ */

/* Fill an fmpq from an integer-like Expr (EXPR_INTEGER or EXPR_BIGINT). */
static int fmpz_from_int_expr(fmpz_t out, const Expr* e) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

static int to_mpoly(const Expr* e, fmpq_mpoly_t out,
                    const fmpq_mpoly_ctx_t ctx, const VarSet* vs) {
    switch (e->type) {
        case EXPR_INTEGER:
            fmpq_mpoly_set_si(out, e->data.integer, ctx);
            return 1;
        case EXPR_BIGINT: {
            fmpq_t c; fmpq_init(c);
            fmpz_set_mpz(fmpq_numref(c), e->data.bigint);
            fmpz_one(fmpq_denref(c));
            fmpq_mpoly_set_fmpq(out, c, ctx);
            fmpq_clear(c);
            return 1;
        }
        case EXPR_SYMBOL: {
            int idx = var_index(vs, e->data.symbol);
            if (idx < 0) return 0;
            fmpq_mpoly_gen(out, idx, ctx);
            return 1;
        }
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) return 0;
            size_t n = e->data.function.arg_count;

            if (strcmp(h, "Rational") == 0) {
                if (n != 2) return 0;
                fmpq_t c; fmpq_init(c);
                if (!fmpz_from_int_expr(fmpq_numref(c), e->data.function.args[0]) ||
                    !fmpz_from_int_expr(fmpq_denref(c), e->data.function.args[1])) {
                    fmpq_clear(c);
                    return 0;
                }
                fmpq_canonicalise(c);
                fmpq_mpoly_set_fmpq(out, c, ctx);
                fmpq_clear(c);
                return 1;
            }
            if (strcmp(h, "Plus") == 0) {
                fmpq_mpoly_zero(out, ctx);
                fmpq_mpoly_t t; fmpq_mpoly_init(t, ctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = to_mpoly(e->data.function.args[i], t, ctx, vs);
                    if (ok) fmpq_mpoly_add(out, out, t, ctx);
                }
                fmpq_mpoly_clear(t, ctx);
                return ok;
            }
            if (strcmp(h, "Times") == 0) {
                fmpq_mpoly_one(out, ctx);
                fmpq_mpoly_t t; fmpq_mpoly_init(t, ctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = to_mpoly(e->data.function.args[i], t, ctx, vs);
                    if (ok) fmpq_mpoly_mul(out, out, t, ctx);
                }
                fmpq_mpoly_clear(t, ctx);
                return ok;
            }
            if (strcmp(h, "Power") == 0) {
                if (n != 2) return 0;
                const Expr* exp = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return 0;
                fmpq_mpoly_t b; fmpq_mpoly_init(b, ctx);
                int ok = to_mpoly(e->data.function.args[0], b, ctx, vs);
                if (ok) ok = fmpq_mpoly_pow_ui(out, b, (ulong)exp->data.integer, ctx);
                fmpq_mpoly_clear(b, ctx);
                return ok;
            }
            return 0;
        }
        default:
            return 0;
    }
}

/* ------------------------------------------------------------------ */
/*  fmpq_mpoly -> Expr                                                 */
/* ------------------------------------------------------------------ */

static Expr* expr_from_mpz_local(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

static Expr* expr_from_fmpq_local(const fmpq_t q) {
    mpz_t nz, dz;
    mpz_init(nz);
    mpz_init(dz);
    fmpz_get_mpz(nz, fmpq_numref(q));
    fmpz_get_mpz(dz, fmpq_denref(q));
    Expr* r;
    if (mpz_cmp_si(dz, 1) == 0) {
        r = expr_from_mpz_local(nz);
    } else {
        Expr* args[2] = { expr_from_mpz_local(nz), expr_from_mpz_local(dz) };
        r = expr_new_function(expr_new_symbol("Rational"), args, 2);
    }
    mpz_clear(nz);
    mpz_clear(dz);
    return r;
}

/*
 * Render an fmpq_mpoly as an (unsimplified but correct) Expr tree:
 * Plus[ Times[coeff, var^e, ...], ... ]. The evaluator canonicalises it after
 * the builtin returns, so we need not produce sorted/normal form here.
 */
static Expr* mpoly_to_expr(const fmpq_mpoly_t P, const fmpq_mpoly_ctx_t ctx,
                           const VarSet* vs) {
    slong len = fmpq_mpoly_length(P, ctx);
    if (len == 0) return expr_new_integer(0);

    slong nv = (slong)vs->count;
    ulong* exps = malloc(sizeof(ulong) * (size_t)(nv > 0 ? nv : 1));
    Expr** terms = malloc(sizeof(Expr*) * (size_t)len);
    if (!exps || !terms) { free(exps); free(terms); return NULL; }

    fmpq_t c;
    fmpq_init(c);
    for (slong i = 0; i < len; i++) {
        fmpq_mpoly_get_term_coeff_fmpq(c, P, i, ctx);
        fmpq_mpoly_get_term_exp_ui(exps, P, i, ctx);

        Expr** factors = malloc(sizeof(Expr*) * (size_t)(nv + 1));
        if (!factors) {
            for (slong j = 0; j < i; j++) expr_free(terms[j]);
            fmpq_clear(c);
            free(exps);
            free(terms);
            return NULL;
        }
        size_t nf = 0;
        factors[nf++] = expr_from_fmpq_local(c);
        for (slong v = 0; v < nv; v++) {
            if (exps[v] == 0) continue;
            Expr* var = expr_new_symbol(vs->names[v]);
            if (exps[v] == 1) {
                factors[nf++] = var;
            } else {
                Expr* pa[2] = { var, expr_new_integer((int64_t)exps[v]) };
                factors[nf++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
        }
        terms[i] = (nf == 1) ? factors[0]
                             : expr_new_function(expr_new_symbol("Times"), factors, nf);
        free(factors);
    }
    fmpq_clear(c);
    free(exps);

    Expr* result = (len == 1)
                       ? terms[0]
                       : expr_new_function(expr_new_symbol("Plus"), terms, (size_t)len);
    free(terms);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Public API                                                        */
/* ------------------------------------------------------------------ */

int flint_bridge_available(void) { return 1; }

Expr* flint_multivariate_gcd(const Expr* a, const Expr* b) {
    if (!a || !b) return NULL;

    VarSet vs;
    memset(&vs, 0, sizeof vs);
    if (!collect_vars(a, &vs) || !collect_vars(b, &vs)) { varset_free(&vs); return NULL; }
    if (vs.count == 0) { varset_free(&vs); return NULL; } /* numeric: classical path */
    qsort(vs.names, vs.count, sizeof(char*), cmp_str);

    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, (slong)vs.count, ORD_LEX);

    fmpq_mpoly_t A, B, G;
    fmpq_mpoly_init(A, ctx);
    fmpq_mpoly_init(B, ctx);
    fmpq_mpoly_init(G, ctx);

    Expr* out = NULL;
    if (to_mpoly(a, A, ctx, &vs) && to_mpoly(b, B, ctx, &vs)) {
        if (fmpq_mpoly_gcd(G, A, B, ctx))
            out = mpoly_to_expr(G, ctx, &vs);
    }

    fmpq_mpoly_clear(A, ctx);
    fmpq_mpoly_clear(B, ctx);
    fmpq_mpoly_clear(G, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    varset_free(&vs);
    return out;
}

/* Multivariate exact division a / b over Q[x_1..x_n]. Returns the quotient
 * when b divides a exactly, NULL otherwise (or out of scope). Mirrors
 * flint_multivariate_gcd's recognition/conversion; uses fmpq_mpoly_divides. */
static Expr* flint_multivariate_divexact(const Expr* a, const Expr* b) {
    if (!a || !b) return NULL;

    VarSet vs;
    memset(&vs, 0, sizeof vs);
    if (!collect_vars(a, &vs) || !collect_vars(b, &vs)) { varset_free(&vs); return NULL; }
    if (vs.count == 0) { varset_free(&vs); return NULL; }
    qsort(vs.names, vs.count, sizeof(char*), cmp_str);

    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, (slong)vs.count, ORD_LEX);

    fmpq_mpoly_t A, B, Q;
    fmpq_mpoly_init(A, ctx);
    fmpq_mpoly_init(B, ctx);
    fmpq_mpoly_init(Q, ctx);

    Expr* out = NULL;
    if (to_mpoly(a, A, ctx, &vs) && to_mpoly(b, B, ctx, &vs) &&
        !fmpq_mpoly_is_zero(B, ctx)) {
        if (fmpq_mpoly_divides(Q, A, B, ctx))   /* 1 iff B | A exactly */
            out = mpoly_to_expr(Q, ctx, &vs);
    }

    fmpq_mpoly_clear(A, ctx);
    fmpq_mpoly_clear(B, ctx);
    fmpq_mpoly_clear(Q, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    varset_free(&vs);
    return out;
}

/* ================================================================== */
/*  M2: univariate GCD over a real quadratic number field Q(sqrt d)    */
/* ================================================================== */

/* Detection state gathered over both GCD operands. */
typedef struct {
    const char* xname;   /* the single polynomial variable (NULL until seen) */
    long        d;       /* radicand of Sqrt[d]; 0 = none seen yet           */
    int         multi_x; /* >1 distinct variable seen (disqualifies)         */
    int         multi_d; /* >1 distinct radical seen (disqualifies)          */
    int         other;   /* an out-of-scope construct was seen               */
} NFDetect;

static int is_perfect_square_l(long n) {
    if (n < 0) return 0;
    long r = (long)(0.5 + sqrt((double)n));
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r * r == n;
}

/* True if `e` is Sqrt[k] = Power[k, Rational[1,2]] with k a positive integer;
 * writes k to *out. */
static int is_sqrt_int(const Expr* e, long* out) {
    if (e->type != EXPR_FUNCTION) return 0;
    const char* h = fn_head_name(e);
    if (!h || strcmp(h, "Power") != 0 || e->data.function.arg_count != 2) return 0;
    const Expr* base = e->data.function.args[0];
    const Expr* exp  = e->data.function.args[1];
    if (base->type != EXPR_INTEGER || base->data.integer <= 0) return 0;
    if (exp->type != EXPR_FUNCTION) return 0;
    const char* eh = fn_head_name(exp);
    if (!eh || strcmp(eh, "Rational") != 0 || exp->data.function.arg_count != 2) return 0;
    const Expr* en = exp->data.function.args[0];
    const Expr* ed = exp->data.function.args[1];
    if (en->type == EXPR_INTEGER && en->data.integer == 1 &&
        ed->type == EXPR_INTEGER && ed->data.integer == 2) {
        *out = base->data.integer;
        return 1;
    }
    return 0;
}

static void nf_detect(const Expr* e, NFDetect* st) {
    if (st->other) return;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return;
        case EXPR_SYMBOL:
            if (!st->xname) st->xname = e->data.symbol;
            else if (strcmp(st->xname, e->data.symbol) != 0) st->multi_x = 1;
            return;
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) { st->other = 1; return; }
            if (strcmp(h, "Rational") == 0) return;
            long k;
            if (is_sqrt_int(e, &k)) {
                if (is_perfect_square_l(k)) { st->other = 1; return; }
                if (st->d == 0) st->d = k;
                else if (st->d != k) st->multi_d = 1;
                return;
            }
            if (strcmp(h, "Power") == 0) {
                if (e->data.function.arg_count != 2) { st->other = 1; return; }
                const Expr* exp = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) { st->other = 1; return; }
                nf_detect(e->data.function.args[0], st);
                return;
            }
            if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    nf_detect(e->data.function.args[i], st);
                return;
            }
            st->other = 1;
            return;
        }
        default:
            st->other = 1;
            return;
    }
}

/*
 * Generalised to ANY single algebraic generator s with minimal polynomial m(s)
 * (degree m >= 2): a "generator spec" recognises an Expr atom as a power s^k of
 * the generator (Sqrt[d] -> s^1; a root of unity (-1)^(p/q) -> s^k), and the
 * core reduces coefficients modulo m(s). This serves both Q(sqrt d) (M2) and
 * Q(zeta_n) (cyclotomic), and is the foundation for the tower collapse.
 */
typedef struct {
    const char* xname;                              /* the polynomial variable    */
    /* If `e` is a generator atom, fill `img` with its value as a polynomial in
     * the single abstract generator s and return 1; otherwise return 0. For a
     * simple generator this is a monomial s^k; for a collapsed tower it is the
     * generator's image in the primitive element. */
    int (*atom_image)(const Expr* e, const void* data, fmpq_poly_t img);
    const void* data;
} GenSpec;

/* Expr -> fmpq_mpoly over { x = gen 0, s = gen 1 }, substituting each generator
 * atom by its image polynomial in s. */
static int to_mpoly_gen(const Expr* e, fmpq_mpoly_t out,
                        const fmpq_mpoly_ctx_t ctx, const GenSpec* gs) {
    if (e->type == EXPR_FUNCTION) {
        fmpq_poly_t img;
        fmpq_poly_init(img);
        if (gs->atom_image(e, gs->data, img)) {
            fmpq_mpoly_zero(out, ctx);
            ulong exp[2];
            exp[0] = 0;
            fmpq_t c;
            fmpq_init(c);
            for (slong j = 0; j <= fmpq_poly_degree(img); j++) {
                fmpq_poly_get_coeff_fmpq(c, img, j);
                if (fmpq_is_zero(c)) continue;
                exp[1] = (ulong)j;
                fmpq_mpoly_set_coeff_fmpq_ui(out, c, exp, ctx);  /* s^j coeff */
            }
            fmpq_clear(c);
            fmpq_poly_clear(img);
            return 1;
        }
        fmpq_poly_clear(img);
    }
    switch (e->type) {
        case EXPR_INTEGER:
            fmpq_mpoly_set_si(out, e->data.integer, ctx);
            return 1;
        case EXPR_BIGINT: {
            fmpq_t c; fmpq_init(c);
            fmpz_set_mpz(fmpq_numref(c), e->data.bigint);
            fmpz_one(fmpq_denref(c));
            fmpq_mpoly_set_fmpq(out, c, ctx);
            fmpq_clear(c);
            return 1;
        }
        case EXPR_SYMBOL:
            if (strcmp(e->data.symbol, gs->xname) == 0) { fmpq_mpoly_gen(out, 0, ctx); return 1; }
            return 0;
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) return 0;
            size_t n = e->data.function.arg_count;
            if (strcmp(h, "Rational") == 0) {
                if (n != 2) return 0;
                fmpq_t c; fmpq_init(c);
                if (!fmpz_from_int_expr(fmpq_numref(c), e->data.function.args[0]) ||
                    !fmpz_from_int_expr(fmpq_denref(c), e->data.function.args[1])) {
                    fmpq_clear(c); return 0;
                }
                fmpq_canonicalise(c);
                fmpq_mpoly_set_fmpq(out, c, ctx);
                fmpq_clear(c);
                return 1;
            }
            if (strcmp(h, "Plus") == 0) {
                fmpq_mpoly_zero(out, ctx);
                fmpq_mpoly_t t; fmpq_mpoly_init(t, ctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = to_mpoly_gen(e->data.function.args[i], t, ctx, gs);
                    if (ok) fmpq_mpoly_add(out, out, t, ctx);
                }
                fmpq_mpoly_clear(t, ctx);
                return ok;
            }
            if (strcmp(h, "Times") == 0) {
                fmpq_mpoly_one(out, ctx);
                fmpq_mpoly_t t; fmpq_mpoly_init(t, ctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = to_mpoly_gen(e->data.function.args[i], t, ctx, gs);
                    if (ok) fmpq_mpoly_mul(out, out, t, ctx);
                }
                fmpq_mpoly_clear(t, ctx);
                return ok;
            }
            if (strcmp(h, "Power") == 0) {
                if (n != 2) return 0;
                const Expr* exp = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return 0;
                fmpq_mpoly_t b; fmpq_mpoly_init(b, ctx);
                int ok = to_mpoly_gen(e->data.function.args[0], b, ctx, gs);
                if (ok) ok = fmpq_mpoly_pow_ui(out, b, (ulong)exp->data.integer, ctx);
                fmpq_mpoly_clear(b, ctx);
                return ok;
            }
            return 0;
        }
        default:
            return 0;
    }
}

/* Distribute an fmpq_mpoly in (x=gen0, s=gen1) into a gr_poly in x over Q(s),
 * reducing each x-coefficient (a polynomial in s) modulo the minimal polynomial. */
static int mpoly_gen_to_grpoly(const fmpq_mpoly_t P, const fmpq_mpoly_ctx_t ctx2,
                               const fmpq_poly_t minpoly, gr_ctx_t nf, gr_poly_t out) {
    slong len = fmpq_mpoly_length(P, ctx2);
    ulong exps[2];
    slong dx = 0;
    for (slong i = 0; i < len; i++) {
        fmpq_mpoly_get_term_exp_ui(exps, P, i, ctx2);
        if ((slong)exps[0] > dx) dx = (slong)exps[0];
    }

    gr_ptr alpha, tmp, coef;
    GR_TMP_INIT(alpha, nf);
    GR_TMP_INIT(tmp, nf);
    GR_TMP_INIT(coef, nf);
    fmpq_t c, cc;
    fmpq_init(c);
    fmpq_init(cc);
    fmpq_poly_t spoly, rem;
    fmpq_poly_init(spoly);
    fmpq_poly_init(rem);

    int status = GR_SUCCESS;
    status |= gr_gen(alpha, nf);
    for (slong i = 0; i <= dx; i++) {
        /* gather this x-power's coefficient as a polynomial in s ... */
        fmpq_poly_zero(spoly);
        for (slong t = 0; t < len; t++) {
            fmpq_mpoly_get_term_exp_ui(exps, P, t, ctx2);
            if ((slong)exps[0] != i) continue;
            fmpq_mpoly_get_term_coeff_fmpq(c, P, t, ctx2);
            fmpq_poly_set_coeff_fmpq(spoly, (slong)exps[1], c);
        }
        /* ... reduce mod m(s), then evaluate the remainder at the generator. */
        fmpq_poly_rem(rem, spoly, minpoly);
        status |= gr_zero(coef, nf);
        for (slong j = fmpq_poly_degree(rem); j >= 0; j--) {     /* Horner in s */
            status |= gr_mul(coef, coef, alpha, nf);
            fmpq_poly_get_coeff_fmpq(cc, rem, j);
            status |= gr_set_fmpq(tmp, cc, nf);
            status |= gr_add(coef, coef, tmp, nf);
        }
        status |= gr_poly_set_coeff_scalar(out, i, coef, nf);
    }

    fmpq_poly_clear(spoly);
    fmpq_poly_clear(rem);
    fmpq_clear(c);
    fmpq_clear(cc);
    GR_TMP_CLEAR(alpha, nf);
    GR_TMP_CLEAR(tmp, nf);
    GR_TMP_CLEAR(coef, nf);
    return status == GR_SUCCESS;
}

/* Build sum_j fp[j] * gen^j as an Expr (gen^0 = constant, gen^1 = gen). */
static Expr* gen_coef_expr(const fmpq_poly_t fp, fmpq_t scratch, const Expr* gen) {
    slong dr = fmpq_poly_degree(fp);
    if (dr < 0) return expr_new_integer(0);

    Expr** parts = malloc(sizeof(Expr*) * (size_t)(dr + 1));
    size_t np = 0;
    for (slong j = 0; j <= dr; j++) {
        fmpq_poly_get_coeff_fmpq(scratch, fp, j);
        if (fmpq_is_zero(scratch)) continue;
        Expr* part;
        if (j == 0) {
            part = expr_from_fmpq_local(scratch);
        } else {
            Expr* gp;
            if (j == 1) {
                gp = expr_copy((Expr*)gen);
            } else {
                Expr* pa[2] = { expr_copy((Expr*)gen), expr_new_integer((int64_t)j) };
                gp = expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
            if (fmpq_is_one(scratch)) {
                part = gp;
            } else {
                Expr* ta[2] = { expr_from_fmpq_local(scratch), gp };
                part = expr_new_function(expr_new_symbol("Times"), ta, 2);
            }
        }
        parts[np++] = part;
    }
    Expr* r;
    if (np == 0)       r = expr_new_integer(0);
    else if (np == 1)  r = parts[0];
    else               r = expr_new_function(expr_new_symbol("Plus"), parts, np);
    free(parts);
    return r;
}

/* A coefficient builder turns one field element (in the power basis, as the
 * fmpq_poly from nf_elem_get_fmpq_poly) into an Expr. The simple-generator and
 * collapsed-tower readbacks differ only in this callback. */
typedef Expr* (*CoefBuilder)(const fmpq_poly_t fp, const void* cdata);

/* Simple generator: coefficient = sum_j fp[j] * gen^j. */
static Expr* single_gen_coef(const fmpq_poly_t fp, const void* cdata) {
    const Expr* gen = *(const Expr* const*)cdata;
    fmpq_t scratch;
    fmpq_init(scratch);
    Expr* r = gen_coef_expr(fp, scratch, gen);
    fmpq_clear(scratch);
    return r;
}

/* gr_poly over Q(s) -> Expr (unsimplified Plus[ coef * x^i ]). */
static Expr* grpoly_to_expr(const gr_poly_t G, gr_ctx_t nf, nf_t my_nf,
                            const char* xname, CoefBuilder cb, const void* cdata) {
    slong len = gr_poly_length(G, nf);
    if (len == 0) return expr_new_integer(0);

    Expr** terms = malloc(sizeof(Expr*) * (size_t)len);
    if (!terms) return NULL;

    gr_ptr coef;
    GR_TMP_INIT(coef, nf);
    fmpq_poly_t fp;
    fmpq_poly_init(fp);

    for (slong i = 0; i < len; i++) {
        GR_IGNORE(gr_poly_get_coeff_scalar(coef, G, i, nf));
        nf_elem_get_fmpq_poly(fp, (nf_elem_struct*)coef, my_nf);
        Expr* coefe = cb(fp, cdata);
        if (i == 0) {
            terms[i] = coefe;
        } else {
            Expr* xp;
            if (i == 1) {
                xp = expr_new_symbol(xname);
            } else {
                Expr* xa[2] = { expr_new_symbol(xname), expr_new_integer((int64_t)i) };
                xp = expr_new_function(expr_new_symbol("Power"), xa, 2);
            }
            Expr* ta[2] = { coefe, xp };
            terms[i] = expr_new_function(expr_new_symbol("Times"), ta, 2);
        }
    }

    fmpq_poly_clear(fp);
    GR_TMP_CLEAR(coef, nf);

    Expr* result = (len == 1) ? terms[0]
                              : expr_new_function(expr_new_symbol("Plus"), terms, (size_t)len);
    free(terms);
    return result;
}

typedef enum { ABSOP_GCD, ABSOP_DIVEXACT } AbsOp;

/* Shared engine: a binary operation (GCD, or exact division a/b) over the simple
 * algebraic extension Q(s) defined by `minpoly`, with `gs` mapping input atoms
 * to images in s and `cb` building Expr coefficients from the result. For
 * ABSOP_DIVEXACT, returns NULL when b does not divide a exactly. */
static Expr* absfield_op_core(const Expr* a, const Expr* b, const char* xname,
                              const fmpq_poly_t minpoly, const GenSpec* gs,
                              CoefBuilder cb, const void* cdata, AbsOp op) {
    gr_ctx_t nf;
    gr_ctx_init_nf(nf, minpoly);
    nf_t my_nf;
    nf_init(my_nf, minpoly);

    fmpq_mpoly_ctx_t ctx2;
    fmpq_mpoly_ctx_init(ctx2, 2, ORD_LEX);       /* gen0 = x, gen1 = s */
    fmpq_mpoly_t A, B;
    fmpq_mpoly_init(A, ctx2);
    fmpq_mpoly_init(B, ctx2);
    gr_poly_t GA, GB, GG;
    gr_poly_init(GA, nf);
    gr_poly_init(GB, nf);
    gr_poly_init(GG, nf);

    Expr* out = NULL;
    if (to_mpoly_gen(a, A, ctx2, gs) && to_mpoly_gen(b, B, ctx2, gs) &&
        mpoly_gen_to_grpoly(A, ctx2, minpoly, nf, GA) &&
        mpoly_gen_to_grpoly(B, ctx2, minpoly, nf, GB)) {
        if (op == ABSOP_GCD) {
            if (gr_poly_gcd(GG, GA, GB, nf) == GR_SUCCESS)
                out = grpoly_to_expr(GG, nf, my_nf, xname, cb, cdata);
        } else {  /* ABSOP_DIVEXACT: a / b, only if the remainder vanishes */
            gr_poly_t R;
            gr_poly_init(R, nf);
            if (gr_poly_divrem(GG, R, GA, GB, nf) == GR_SUCCESS &&
                gr_poly_length(R, nf) == 0)
                out = grpoly_to_expr(GG, nf, my_nf, xname, cb, cdata);
            gr_poly_clear(R, nf);
        }
    }

    gr_poly_clear(GA, nf);
    gr_poly_clear(GB, nf);
    gr_poly_clear(GG, nf);
    fmpq_mpoly_clear(A, ctx2);
    fmpq_mpoly_clear(B, ctx2);
    fmpq_mpoly_ctx_clear(ctx2);
    nf_clear(my_nf);
    gr_ctx_clear(nf);
    return out;
}

/* --- Q(sqrt d): single radical generator -------------------------------- */

static int sqrt_atom_image(const Expr* e, const void* data, fmpq_poly_t img) {
    long k;
    if (is_sqrt_int(e, &k) && k == *(const long*)data) {
        fmpq_poly_set_coeff_si(img, 1, 1);   /* img = s */
        return 1;
    }
    return 0;
}

static Expr* numberfield_op(const Expr* a, const Expr* b, AbsOp op) {
    if (!a || !b) return NULL;

    NFDetect st;
    memset(&st, 0, sizeof st);
    nf_detect(a, &st);
    nf_detect(b, &st);
    if (st.other || st.multi_x || st.multi_d) return NULL;
    if (st.d == 0) return NULL;                 /* no radical: rational path's job */
    const char* xname = st.xname ? st.xname : "x";

    fmpq_poly_t minpoly;                         /* y^2 - d */
    fmpq_poly_init(minpoly);
    fmpq_poly_set_coeff_si(minpoly, 2, 1);
    fmpq_poly_set_coeff_si(minpoly, 0, -st.d);

    long d = st.d;
    GenSpec gs = { xname, sqrt_atom_image, &d };
    Expr* sd_args[1] = { expr_new_integer(d) };
    Expr* gen = expr_new_function(expr_new_symbol("Sqrt"), sd_args, 1);

    Expr* out = absfield_op_core(a, b, xname, minpoly, &gs, single_gen_coef, &gen, op);
    expr_free(gen);
    fmpq_poly_clear(minpoly);
    return out;
}

Expr* flint_numberfield_gcd(const Expr* a, const Expr* b) {
    return numberfield_op(a, b, ABSOP_GCD);
}

/* --- Q(zeta_n): single cyclotomic generator ----------------------------- */

static long lcm_l(long a, long b) {
    if (a == 0 || b == 0) return 0;
    long g, x = a, y = b;
    while (y) { g = x % y; x = y; y = g; }
    return (a / x) * b;
}

/* A root of unity Power[-1, Rational[p, q]] -> (p, q); returns 1 on match. */
static int as_root_of_unity(const Expr* e, long* p, long* q) {
    if (e->type != EXPR_FUNCTION) return 0;
    const char* h = fn_head_name(e);
    if (!h || strcmp(h, "Power") != 0 || e->data.function.arg_count != 2) return 0;
    const Expr* base = e->data.function.args[0];
    const Expr* exp  = e->data.function.args[1];
    if (!(base->type == EXPR_INTEGER && base->data.integer == -1)) return 0;
    if (exp->type != EXPR_FUNCTION) return 0;
    const char* eh = fn_head_name(exp);
    if (!eh || strcmp(eh, "Rational") != 0 || exp->data.function.arg_count != 2) return 0;
    const Expr* en = exp->data.function.args[0];
    const Expr* ed = exp->data.function.args[1];
    if (en->type != EXPR_INTEGER || ed->type != EXPR_INTEGER || ed->data.integer < 2) return 0;
    *p = en->data.integer;
    *q = ed->data.integer;
    return 1;
}

typedef struct { long Q; long N; } CycData;

/* (-1)^(p/q) = zeta_N^{p * (N/(2q))}, with N = 2Q and q | Q. */
static int cyc_atom_image(const Expr* e, const void* data, fmpq_poly_t img) {
    const CycData* cd = (const CycData*)data;
    long p, q;
    if (!as_root_of_unity(e, &p, &q)) return 0;
    long k = (p * (cd->Q / q)) % cd->N;
    if (k < 0) k += cd->N;
    fmpq_poly_set_coeff_si(img, k, 1);   /* img = s^k */
    return 1;
}

typedef struct {
    const char* xname;
    long Q;            /* lcm of root-of-unity denominators; 0 = none           */
    int  multi_x;
    int  other;        /* radical / Complex / unsupported head seen             */
} CycDetect;

static void cyc_detect(const Expr* e, CycDetect* st) {
    if (st->other) return;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return;
        case EXPR_SYMBOL:
            if (!st->xname) st->xname = e->data.symbol;
            else if (strcmp(st->xname, e->data.symbol) != 0) st->multi_x = 1;
            return;
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) { st->other = 1; return; }
            long k, p, q;
            if (is_sqrt_int(e, &k)) { st->other = 1; return; }  /* radical: tower, not here */
            if (as_root_of_unity(e, &p, &q)) {
                st->Q = st->Q ? lcm_l(st->Q, q) : q;
                return;
            }
            if (strcmp(h, "Rational") == 0) return;
            if (strcmp(h, "Power") == 0) {
                if (e->data.function.arg_count != 2) { st->other = 1; return; }
                const Expr* exp = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) { st->other = 1; return; }
                cyc_detect(e->data.function.args[0], st);
                return;
            }
            if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    cyc_detect(e->data.function.args[i], st);
                return;
            }
            st->other = 1;
            return;
        }
        default:
            st->other = 1;
            return;
    }
}

static Expr* cyclotomic_op(const Expr* a, const Expr* b, AbsOp op) {
    if (!a || !b) return NULL;

    CycDetect st;
    memset(&st, 0, sizeof st);
    cyc_detect(a, &st);
    cyc_detect(b, &st);
    if (st.other || st.multi_x) return NULL;
    if (st.Q == 0) return NULL;                  /* no root of unity */
    const char* xname = st.xname ? st.xname : "x";

    long N = 2 * st.Q;                           /* generator zeta_N = (-1)^(1/Q) */
    fmpz_poly_t cy;
    fmpz_poly_init(cy);
    fmpz_poly_cyclotomic(cy, (ulong)N);          /* minimal polynomial Phi_N */
    fmpq_poly_t minpoly;
    fmpq_poly_init(minpoly);
    fmpq_poly_set_fmpz_poly(minpoly, cy);
    fmpz_poly_clear(cy);

    CycData cd = { st.Q, N };
    GenSpec gs = { xname, cyc_atom_image, &cd };
    /* gen = (-1)^(1/Q) = Power[-1, Rational[1, Q]] */
    Expr* rargs[2] = { expr_new_integer(1), expr_new_integer(st.Q) };
    Expr* rat = expr_new_function(expr_new_symbol("Rational"), rargs, 2);
    Expr* pargs[2] = { expr_new_integer(-1), rat };
    Expr* gen = expr_new_function(expr_new_symbol("Power"), pargs, 2);

    Expr* out = absfield_op_core(a, b, xname, minpoly, &gs, single_gen_coef, &gen, op);
    expr_free(gen);
    fmpq_poly_clear(minpoly);
    return out;
}

Expr* flint_cyclotomic_gcd(const Expr* a, const Expr* b) {
    return cyclotomic_op(a, b, ABSOP_GCD);
}

/* --- Q(sqrt d_1, ..., sqrt d_r): radical tower via primitive-element collapse - */

#define TOWER_MAXGEN 4    /* up to a 2^4 = 16-dimensional product basis */

typedef struct {
    const char* xname;
    int  r;                      /* number of distinct square-root generators */
    long d[TOWER_MAXGEN];        /* their (squarefree) radicands              */
    int  multi_x;
    int  other;
} TowerDetect;

static void tower_detect(const Expr* e, TowerDetect* st) {
    if (st->other) return;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return;
        case EXPR_SYMBOL:
            if (!st->xname) st->xname = e->data.symbol;
            else if (strcmp(st->xname, e->data.symbol) != 0) st->multi_x = 1;
            return;
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) { st->other = 1; return; }
            long k, p, q;
            if (is_sqrt_int(e, &k)) {
                if (is_perfect_square_l(k)) { st->other = 1; return; }
                for (int i = 0; i < st->r; i++) if (st->d[i] == k) return;
                if (st->r >= TOWER_MAXGEN) { st->other = 1; return; }
                st->d[st->r++] = k;
                return;
            }
            if (as_root_of_unity(e, &p, &q)) { st->other = 1; return; }  /* not pure radical */
            if (strcmp(h, "Rational") == 0) return;
            if (strcmp(h, "Power") == 0) {
                if (e->data.function.arg_count != 2) { st->other = 1; return; }
                const Expr* exp = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER || exp->data.integer < 0) { st->other = 1; return; }
                tower_detect(e->data.function.args[0], st);
                return;
            }
            if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    tower_detect(e->data.function.args[i], st);
                return;
            }
            st->other = 1;
            return;
        }
        default:
            st->other = 1;
            return;
    }
}

/*
 * Collapse Q(sqrt d_1, ..., sqrt d_r) to a single primitive element theta via
 * fully-rational linear algebra: each generator multiplies by a known matrix in
 * the product (bitmask) basis {prod sqrt(d_i)^{e_i}}, so T = sum c_j (mult-by-
 * sqrt d_j) is multiplication by theta = sum c_j sqrt(d_j). The columns
 * {1, theta, ..., theta^{D-1}} (in product coords) form P; if P is invertible,
 * theta is primitive (degree D = 2^r) and:
 *   - minpoly(theta) is read from theta^D = P (Pinv theta^D);
 *   - sqrt(d_i)'s image in the power basis is Pinv applied to its product-basis
 *     unit vector.
 * Fills `minpoly`, `P` (power->product transform, for readback) and `img[i]`
 * (sqrt d_i as a polynomial in theta). Returns 0 if the radicals are not
 * independent (no primitive theta of full degree found).
 */
static int tower_build(const TowerDetect* det, fmpq_poly_t minpoly,
                       fmpq_mat_t P, fmpq_poly_struct* img) {
    int r = det->r;
    slong D = (slong)1 << r;
    static const long cvecs[5][TOWER_MAXGEN] = {
        {1,2,3,4},{1,3,5,7},{2,5,7,11},{1,4,9,16},{3,7,13,17}
    };

    fmpq_mat_t T, Pinv, vcur, vnext, mcol, unit;
    fmpq_mat_init(T, D, D);
    fmpq_mat_init(Pinv, D, D);
    fmpq_mat_init(vcur, D, 1);
    fmpq_mat_init(vnext, D, 1);
    fmpq_mat_init(mcol, D, 1);
    fmpq_mat_init(unit, D, 1);

    int ok = 0;
    for (int attempt = 0; attempt < 5 && !ok; attempt++) {
        const long* c = cvecs[attempt];
        /* T = multiply-by-theta in the product basis */
        for (slong i = 0; i < D; i++)
            for (slong jj = 0; jj < D; jj++)
                fmpq_zero(fmpq_mat_entry(T, i, jj));
        for (slong e = 0; e < D; e++) {
            for (int j = 0; j < r; j++) {
                slong bit = (slong)1 << j, res;
                long mult;
                if (e & bit) { res = e & ~bit; mult = c[j] * det->d[j]; }  /* sqrt d_j squared */
                else         { res = e |  bit; mult = c[j]; }
                fmpq_t add; fmpq_init(add); fmpq_set_si(add, mult, 1);
                fmpq_add(fmpq_mat_entry(T, res, e), fmpq_mat_entry(T, res, e), add);
                fmpq_clear(add);
            }
        }
        /* P columns: v_0 = "1" (mask 0), v_k = T v_{k-1} */
        for (slong i = 0; i < D; i++) fmpq_zero(fmpq_mat_entry(vcur, i, 0));
        fmpq_set_si(fmpq_mat_entry(vcur, 0, 0), 1, 1);
        for (slong i = 0; i < D; i++) fmpq_set(fmpq_mat_entry(P, i, 0), fmpq_mat_entry(vcur, i, 0));
        for (slong k = 1; k < D; k++) {
            fmpq_mat_mul(vnext, T, vcur);
            for (slong i = 0; i < D; i++) fmpq_set(fmpq_mat_entry(P, i, k), fmpq_mat_entry(vnext, i, 0));
            fmpq_mat_set(vcur, vnext);
        }
        if (!fmpq_mat_inv(Pinv, P)) continue;    /* singular: theta not primitive, retry */
        ok = 1;

        /* minpoly: theta^D = sum m_k theta^k, m = Pinv (T v_{D-1}) */
        fmpq_mat_mul(vnext, T, vcur);
        fmpq_mat_mul(mcol, Pinv, vnext);
        fmpq_poly_zero(minpoly);
        fmpq_poly_set_coeff_si(minpoly, D, 1);
        for (slong k = 0; k < D; k++) {
            fmpq_t neg; fmpq_init(neg);
            fmpq_neg(neg, fmpq_mat_entry(mcol, k, 0));
            fmpq_poly_set_coeff_fmpq(minpoly, k, neg);
            fmpq_clear(neg);
        }
        /* images: sqrt(d_i) = Pinv * (its product-basis unit vector) */
        for (int i = 0; i < r; i++) {
            for (slong t = 0; t < D; t++) fmpq_zero(fmpq_mat_entry(unit, t, 0));
            fmpq_set_si(fmpq_mat_entry(unit, (slong)1 << i, 0), 1, 1);
            fmpq_mat_mul(mcol, Pinv, unit);
            fmpq_poly_zero(&img[i]);
            for (slong k = 0; k < D; k++) {
                fmpq* v = fmpq_mat_entry(mcol, k, 0);
                if (!fmpq_is_zero(v)) fmpq_poly_set_coeff_fmpq(&img[i], k, v);
            }
        }
    }

    fmpq_mat_clear(T);
    fmpq_mat_clear(Pinv);
    fmpq_mat_clear(vcur);
    fmpq_mat_clear(vnext);
    fmpq_mat_clear(mcol);
    fmpq_mat_clear(unit);
    return ok;
}

typedef struct { const TowerDetect* det; const fmpq_poly_struct* img; } TowerGen;

static int tower_atom_image(const Expr* e, const void* data, fmpq_poly_t img) {
    const TowerGen* tg = (const TowerGen*)data;
    long k;
    if (!is_sqrt_int(e, &k)) return 0;
    for (int i = 0; i < tg->det->r; i++)
        if (tg->det->d[i] == k) { fmpq_poly_set(img, &tg->img[i]); return 1; }
    return 0;
}

/* The product-basis element for `mask`: product of Sqrt[d_i] over set bits. */
static Expr* tower_basis_expr(slong mask, const TowerDetect* det) {
    Expr* factors[TOWER_MAXGEN];
    size_t nf = 0;
    for (int i = 0; i < det->r; i++) {
        if (!(mask & ((slong)1 << i))) continue;
        Expr* sd[1] = { expr_new_integer(det->d[i]) };
        factors[nf++] = expr_new_function(expr_new_symbol("Sqrt"), sd, 1);
    }
    if (nf == 1) return factors[0];
    return expr_new_function(expr_new_symbol("Times"), factors, nf);
}

typedef struct { const TowerDetect* det; const fmpq_mat_struct* P; } TowerRD;

/* Tower readback: power-basis coords fp -> product coords (P fp) -> Expr. */
static Expr* tower_coef(const fmpq_poly_t fp, const void* cdata) {
    const TowerRD* rd = (const TowerRD*)cdata;
    slong D = (slong)1 << rd->det->r;

    fmpq_mat_t w, pc;
    fmpq_mat_init(w, D, 1);
    fmpq_mat_init(pc, D, 1);
    fmpq_t cj;
    fmpq_init(cj);
    for (slong k = 0; k < D; k++) {
        fmpq_poly_get_coeff_fmpq(cj, fp, k);
        fmpq_set(fmpq_mat_entry(w, k, 0), cj);
    }
    fmpq_mat_mul(pc, rd->P, w);

    Expr** parts = malloc(sizeof(Expr*) * (size_t)D);
    size_t np = 0;
    for (slong mask = 0; mask < D; mask++) {
        fmpq* v = fmpq_mat_entry(pc, mask, 0);
        if (fmpq_is_zero(v)) continue;
        if (mask == 0) {
            parts[np++] = expr_from_fmpq_local(v);
        } else {
            Expr* basis = tower_basis_expr(mask, rd->det);
            if (fmpq_is_one(v)) {
                parts[np++] = basis;
            } else {
                Expr* ta[2] = { expr_from_fmpq_local(v), basis };
                parts[np++] = expr_new_function(expr_new_symbol("Times"), ta, 2);
            }
        }
    }

    fmpq_clear(cj);
    fmpq_mat_clear(w);
    fmpq_mat_clear(pc);

    Expr* r;
    if (np == 0)      r = expr_new_integer(0);
    else if (np == 1) r = parts[0];
    else              r = expr_new_function(expr_new_symbol("Plus"), parts, np);
    free(parts);
    return r;
}

static Expr* tower_op(const Expr* a, const Expr* b, AbsOp op) {
    if (!a || !b) return NULL;

    TowerDetect det;
    memset(&det, 0, sizeof det);
    tower_detect(a, &det);
    tower_detect(b, &det);
    if (det.other || det.multi_x) return NULL;
    if (det.r < 2) return NULL;                  /* 0/1 radical: other paths */
    const char* xname = det.xname ? det.xname : "x";

    slong D = (slong)1 << det.r;
    fmpq_poly_t minpoly;
    fmpq_poly_init(minpoly);
    fmpq_mat_t P;
    fmpq_mat_init(P, D, D);
    fmpq_poly_struct img[TOWER_MAXGEN];
    for (int i = 0; i < det.r; i++) fmpq_poly_init(&img[i]);

    Expr* out = NULL;
    if (tower_build(&det, minpoly, P, img)) {
        TowerGen tg = { &det, img };
        GenSpec gs = { xname, tower_atom_image, &tg };
        TowerRD rd = { &det, P };
        out = absfield_op_core(a, b, xname, minpoly, &gs, tower_coef, &rd, op);
    }

    for (int i = 0; i < det.r; i++) fmpq_poly_clear(&img[i]);
    fmpq_mat_clear(P);
    fmpq_poly_clear(minpoly);
    return out;
}

Expr* flint_tower_gcd(const Expr* a, const Expr* b) {
    return tower_op(a, b, ABSOP_GCD);
}

/* ================================================================== */
/*  M3: parametric single radical  Q(t_1..t_p)(sqrt k), k a SYMBOL      */
/* ================================================================== */
/*
 * The Goursat blocker ring is Q(a, b, k)(sqrt k)[x]: parameters a, b, k and a
 * single radical generator sqrt(k) whose radicand k is itself a parameter.
 * Because sqrt(k) is transcendental over Q(a, b) with k = (sqrt k)^2, this
 * field is *isomorphic to the rational function field* Q(a, b, sqrt k) — there
 * is no genuine algebraic relation to reduce modulo. So GCD/division over it
 * collapses to ORDINARY multivariate GCD over Q after the substitution
 *
 *     sqrt(k) -> S,   k -> S^2
 *
 * (S a fresh transcendental), handled directly by fmpq_mpoly_gcd. The readback
 * S -> sqrt(k) then lets S^2 fold back to k under evaluation. This is the whole
 * `p >= 1, r = 1, deg = 2` regime — no fq_nmod residue fields, no CRT, no outer
 * loop are needed here (those are for constant-coefficient minimal polynomials
 * that define a genuine GF(p^d), a later milestone). It is fully rigorous:
 * fmpq_mpoly_gcd is exact and the substitution is a field isomorphism.
 */

/* Fresh transcendental standing in for sqrt(k). Chosen not to collide with any
 * user symbol (leading '$', internal-marker style). */
#define PS_RADSYM "$flint$sqrtk$"

/* Parse `e` as an integer or Rational[p,q]; returns 1 and writes p,q (q>0). */
static int ps_as_ratio(const Expr* e, long* p, long* q) {
    if (e->type == EXPR_INTEGER) { *p = (long)e->data.integer; *q = 1; return 1; }
    if (e->type == EXPR_FUNCTION) {
        const char* h = fn_head_name(e);
        if (h && strcmp(h, "Rational") == 0 && e->data.function.arg_count == 2) {
            const Expr* n = e->data.function.args[0];
            const Expr* d = e->data.function.args[1];
            if (n->type == EXPR_INTEGER && d->type == EXPR_INTEGER && d->data.integer > 0) {
                *p = (long)n->data.integer; *q = (long)d->data.integer; return 1;
            }
        }
    }
    return 0;
}

typedef struct {
    const char* kname;   /* the single symbolic radicand (NULL until seen) */
    int has_sqrt;        /* at least one Sqrt[k] occurred                  */
    int bad;             /* out-of-scope construct                         */
    int allow_neg_pow;   /* accept Power[..,-n] (rational-function coeffs) */
} PSDetect;

/* Recognise Q[params, x, Sqrt[k]] with a single symbolic radicand k, all of
 * whose powers are integer or half-integer. Any other radical (Sqrt[int],
 * Sqrt[other symbol], cube root, root of unity), negative power, or unknown
 * head disqualifies (bad = 1). */
static void ps_detect(const Expr* e, PSDetect* st) {
    if (st->bad) return;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_SYMBOL:
            return;
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) { st->bad = 1; return; }
            if (strcmp(h, "Rational") == 0) return;
            if (strcmp(h, "Power") == 0) {
                if (e->data.function.arg_count != 2) { st->bad = 1; return; }
                const Expr* base = e->data.function.args[0];
                const Expr* exp  = e->data.function.args[1];
                long p, q;
                if (!ps_as_ratio(exp, &p, &q)) { st->bad = 1; return; }
                if (q == 1) {
                    if (p < 0 && !st->allow_neg_pow) { st->bad = 1; return; }  /* denominator: not a polynomial */
                    ps_detect(base, st);
                    return;
                }
                /* Fractional power: only sqrt of the single symbolic radicand. */
                if (q != 2) { st->bad = 1; return; }
                if (base->type != EXPR_SYMBOL) { st->bad = 1; return; }
                if (!st->kname) st->kname = base->data.symbol;
                else if (strcmp(st->kname, base->data.symbol) != 0) { st->bad = 1; return; }
                st->has_sqrt = 1;
                return;
            }
            if (strcmp(h, "Plus") == 0 || strcmp(h, "Times") == 0) {
                for (size_t i = 0; i < e->data.function.arg_count; i++)
                    ps_detect(e->data.function.args[i], st);
                return;
            }
            st->bad = 1;
            return;
        }
        default:
            st->bad = 1;
            return;
    }
}

/* Substitute the radical away: bare k -> S^2, Power[k, p/2] -> S^p,
 * Power[k, m] (integer m>=0) -> S^(2m). S = PS_RADSYM. Returns a fresh tree. */
static Expr* ps_subst_in(const Expr* e, const char* kname) {
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return expr_copy((Expr*)e);
        case EXPR_SYMBOL:
            if (strcmp(e->data.symbol, kname) == 0) {
                Expr* pa[2] = { expr_new_symbol(PS_RADSYM), expr_new_integer(2) };
                return expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
            return expr_copy((Expr*)e);
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (h && strcmp(h, "Power") == 0 && e->data.function.arg_count == 2 &&
                e->data.function.args[0]->type == EXPR_SYMBOL &&
                strcmp(e->data.function.args[0]->data.symbol, kname) == 0) {
                long p, q;
                if (ps_as_ratio(e->data.function.args[1], &p, &q)) {
                    long s_exp = (q == 1) ? 2 * p : p;   /* k^m -> S^(2m); k^(p/2) -> S^p */
                    if (s_exp == 1) return expr_new_symbol(PS_RADSYM);
                    Expr* pa[2] = { expr_new_symbol(PS_RADSYM), expr_new_integer(s_exp) };
                    return expr_new_function(expr_new_symbol("Power"), pa, 2);
                }
            }
            size_t n = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
            for (size_t i = 0; i < n; i++)
                args[i] = ps_subst_in(e->data.function.args[i], kname);
            Expr* head = ps_subst_in(e->data.function.head, kname);
            Expr* r = expr_new_function(head, args, n);
            free(args);
            return r;
        }
        default:
            return expr_copy((Expr*)e);
    }
}

/* Readback: replace the fresh symbol S by Sqrt[k] = Power[k, 1/2]. The caller
 * evaluates the result, folding S^2 -> k, S^3 -> k Sqrt[k], etc. */
static Expr* ps_subst_out(const Expr* e, const char* kname) {
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
            return expr_copy((Expr*)e);
        case EXPR_SYMBOL:
            if (strcmp(e->data.symbol, PS_RADSYM) == 0) {
                Expr* ra[2] = { expr_new_integer(1), expr_new_integer(2) };
                Expr* half = expr_new_function(expr_new_symbol("Rational"), ra, 2);
                Expr* pa[2] = { expr_new_symbol(kname), half };
                return expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
            return expr_copy((Expr*)e);
        case EXPR_FUNCTION: {
            size_t n = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
            for (size_t i = 0; i < n; i++)
                args[i] = ps_subst_out(e->data.function.args[i], kname);
            Expr* head = ps_subst_out(e->data.function.head, kname);
            Expr* r = expr_new_function(head, args, n);
            free(args);
            return r;
        }
        default:
            return expr_copy((Expr*)e);
    }
}

/* GCD (or exact division) over Q(t_1..t_p)(sqrt k) via the k = S^2 collapse. */
static Expr* parametric_sqrt_op(const Expr* a, const Expr* b, AbsOp op) {
    if (!a || !b) return NULL;
    PSDetect st;
    memset(&st, 0, sizeof st);
    ps_detect(a, &st);
    ps_detect(b, &st);
    if (st.bad || !st.has_sqrt || !st.kname) return NULL;

    Expr* a2 = ps_subst_in(a, st.kname);
    Expr* b2 = ps_subst_in(b, st.kname);

    Expr* g = NULL;
    if (a2 && b2) {
        if (op == ABSOP_GCD)
            g = flint_multivariate_gcd(a2, b2);
        else  /* ABSOP_DIVEXACT: exact quotient a/b over Q[params, S, x] */
            g = flint_multivariate_divexact(a2, b2);
    }
    expr_free(a2);
    expr_free(b2);
    if (!g) return NULL;

    Expr* out = ps_subst_out(g, st.kname);
    expr_free(g);
    return eval_and_free(out);   /* fold S^2 -> k and canonicalise */
}

Expr* flint_parametric_sqrt_gcd(const Expr* a, const Expr* b) {
    return parametric_sqrt_op(a, b, ABSOP_GCD);
}

/* Multivariate resultant Res_var(a, b) over Q[x_1..x_n] (all other variables
 * treated as coefficients). NULL when out of scope or `var` is absent. */
static Expr* flint_multivariate_resultant(const Expr* a, const Expr* b,
                                          const char* varname) {
    if (!a || !b || !varname) return NULL;

    VarSet vs;
    memset(&vs, 0, sizeof vs);
    if (!collect_vars(a, &vs) || !collect_vars(b, &vs)) { varset_free(&vs); return NULL; }
    if (vs.count == 0 || var_index(&vs, varname) < 0) { varset_free(&vs); return NULL; }
    qsort(vs.names, vs.count, sizeof(char*), cmp_str);
    slong vi = (slong)var_index(&vs, varname);

    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, (slong)vs.count, ORD_LEX);
    fmpq_mpoly_t A, B, R;
    fmpq_mpoly_init(A, ctx);
    fmpq_mpoly_init(B, ctx);
    fmpq_mpoly_init(R, ctx);

    Expr* out = NULL;
    if (to_mpoly(a, A, ctx, &vs) && to_mpoly(b, B, ctx, &vs)) {
        if (fmpq_mpoly_resultant(R, A, B, vi, ctx))
            out = mpoly_to_expr(R, ctx, &vs);
    }

    fmpq_mpoly_clear(A, ctx);
    fmpq_mpoly_clear(B, ctx);
    fmpq_mpoly_clear(R, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    varset_free(&vs);
    return out;
}

/*
 * Res_var(a, b) over the parametric radical field Q(t..)(sqrt k). Same collapse
 * as the GCD: sqrt(k) -> S, k -> S^2 makes a, b ordinary polynomials over Q, so
 * FLINT's fmpq_mpoly_resultant computes the (multivariate) resultant fast; then
 * read S back to sqrt(k). This is what makes Rothstein-Trager's
 * Res_x(A - t D', D) over Q(a,b,k)(sqrt k) tractable — the classical symbolic
 * subresultant PRS blows up on the parametric radical coefficients.
 */
Expr* flint_parametric_sqrt_resultant(const Expr* a, const Expr* b,
                                      const Expr* var) {
    if (!a || !b || !var || var->type != EXPR_SYMBOL) return NULL;
    PSDetect st;
    memset(&st, 0, sizeof st);
    ps_detect(a, &st);
    ps_detect(b, &st);
    if (st.bad || !st.has_sqrt || !st.kname) return NULL;
    if (strcmp(var->data.symbol, st.kname) == 0) return NULL;  /* var == radicand */

    Expr* a2 = ps_subst_in(a, st.kname);
    Expr* b2 = ps_subst_in(b, st.kname);
    Expr* r = NULL;
    if (a2 && b2) r = flint_multivariate_resultant(a2, b2, var->data.symbol);
    expr_free(a2);
    expr_free(b2);
    if (!r) return NULL;

    Expr* out = ps_subst_out(r, st.kname);
    expr_free(r);
    return eval_and_free(out);
}

/* Multivariate factorisation over Q[x_1..x_n]. Renders the result as an
 * (unevaluated) Times[ constant, base_1^exp_1, ... ] Expr the caller owns.
 * Returns NULL out of scope. `squarefree` selects squarefree factorisation
 * (group by multiplicity, no irreducible split) instead of full factoring. */
static Expr* flint_multivariate_factor_impl(const Expr* p, int squarefree) {
    if (!p) return NULL;

    VarSet vs;
    memset(&vs, 0, sizeof vs);
    if (!collect_vars(p, &vs)) { varset_free(&vs); return NULL; }
    if (vs.count == 0) { varset_free(&vs); return NULL; }
    qsort(vs.names, vs.count, sizeof(char*), cmp_str);

    /* ORD_DEGLEX so that term 0 of each factor is its highest-total-degree
     * (deglex tie-break) leading term — used below to normalise the factor's
     * sign to Mathematica's positive-leading-coefficient convention. */
    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, (slong)vs.count, ORD_DEGLEX);
    fmpq_mpoly_t P;
    fmpq_mpoly_init(P, ctx);
    fmpq_mpoly_factor_t fac;
    fmpq_mpoly_factor_init(fac, ctx);

    Expr* out = NULL;
    int ok = to_mpoly(p, P, ctx, &vs) &&
             (squarefree ? fmpq_mpoly_factor_squarefree(fac, P, ctx)
                         : fmpq_mpoly_factor(fac, P, ctx));
    if (ok) {
        slong nf = fmpq_mpoly_factor_length(fac, ctx);
        fmpq_t c; fmpq_init(c);
        fmpq_mpoly_factor_get_constant_fmpq(c, fac, ctx);

        /* Slot 0 is reserved for the (possibly sign-adjusted) content, which
         * is filled in after the factor loop so that sign flips can be folded
         * into it.  Polynomial factors occupy slots 1.. */
        size_t nterms = (size_t)nf + 1;
        Expr** factors = malloc(sizeof(Expr*) * nterms);
        factors[0] = NULL;
        size_t nt = 1;

        int sign_flip = 0;      /* parity of factor-sign flips, folded into c */
        fmpq_t lc; fmpq_init(lc);
        fmpq_mpoly_t base; fmpq_mpoly_init(base, ctx);
        for (slong i = 0; i < nf; i++) {
            fmpq_mpoly_factor_get_base(base, fac, i, ctx);
            slong e = fmpq_mpoly_factor_get_exp_si(fac, i, ctx);
            /* Normalise the factor to a positive leading coefficient (highest
             * total degree, deglex tie-break) — Mathematica's convention —
             * folding the discarded -1 into the content (via (-1)^e) so the
             * product is unchanged.  Avoids ugly forms like -(k - u^2). */
            if (fmpq_mpoly_length(base, ctx) > 0) {
                fmpq_mpoly_get_term_coeff_fmpq(lc, base, 0, ctx);
                if (fmpq_sgn(lc) < 0) {
                    fmpq_mpoly_neg(base, base, ctx);
                    if (e & 1) sign_flip ^= 1;
                }
            }
            Expr* be = mpoly_to_expr(base, ctx, &vs);
            if (!be) continue;
            if (e == 1) {
                factors[nt++] = be;
            } else {
                Expr* pa[2] = { be, expr_new_integer((int64_t)e) };
                factors[nt++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
        }
        fmpq_mpoly_clear(base, ctx);
        fmpq_clear(lc);

        if (sign_flip) fmpq_neg(c, c);
        if (!fmpq_is_one(c)) factors[0] = expr_from_fmpq_local(c);
        fmpq_clear(c);

        size_t start = factors[0] ? 0 : 1;   /* skip empty content slot */
        size_t count = nt - start;
        if (count == 0)      out = expr_new_integer(1);
        else if (count == 1) out = factors[start];
        else                 out = expr_new_function(expr_new_symbol("Times"),
                                                      factors + start, count);
        free(factors);
    }

    fmpq_mpoly_factor_clear(fac, ctx);
    fmpq_mpoly_clear(P, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    varset_free(&vs);
    return out;
}

/*
 * Factor (or squarefree-factor) a polynomial over the parametric radical field
 * Q(t..)(sqrt k), k a free symbol, via the sqrt(k) -> S, k -> S^2 collapse: the
 * field is the rational function field Q(t.., sqrt k), so factoring reduces to
 * ORDINARY multivariate factoring over Q (fmpq_mpoly_factor) — which FLINT does,
 * unlike the constant-radicand number-field case (gr_poly_factor is GR_UNABLE).
 * Returns NULL when there is no symbolic radical or for out-of-scope input.
 */
static Expr* parametric_sqrt_factor(const Expr* p, int squarefree) {
    if (!p) return NULL;
    PSDetect st;
    memset(&st, 0, sizeof st);
    ps_detect(p, &st);
    if (st.bad || !st.has_sqrt || !st.kname) return NULL;

    Expr* p2 = ps_subst_in(p, st.kname);
    if (!p2) return NULL;
    Expr* f = flint_multivariate_factor_impl(p2, squarefree);
    expr_free(p2);
    if (!f) return NULL;

    Expr* out = ps_subst_out(f, st.kname);
    expr_free(f);
    return eval_and_free(out);
}

Expr* flint_parametric_sqrt_factor(const Expr* p) {
    return parametric_sqrt_factor(p, 0);
}

Expr* flint_parametric_sqrt_factor_squarefree(const Expr* p) {
    return parametric_sqrt_factor(p, 1);
}

/* ================================================================== */
/*  Univariate polynomials over the rational function field Q(vars),   */
/*  via gr_poly over fmpz_mpoly_q. Used for xgcd/gcd/divrem over        */
/*  Q(params)(sqrt k)[x] (the parametric radical field, sqrt k -> S).   */
/* ================================================================== */

/* Collect EVERY symbol occurring in `e` (unlike collect_vars this does not
 * validate polynomial-ness — coefficients may be rational functions). */
static void collect_all_symbols(const Expr* e, VarSet* vs) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) { varset_add(vs, e->data.symbol); return; }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_all_symbols(e->data.function.args[i], vs);
}

/* Expr (a rational function in the field variables `fvars`, i.e. the input
 * already had sqrt(k) -> S applied) -> fmpz_mpoly_q. Handles +,-(via Times),
 * *, /, integer powers (positive and negative). Returns 1 on success. */
static int expr_to_mpolyq(const Expr* e, fmpz_mpoly_q_t out,
                          const fmpz_mpoly_ctx_t mctx, const VarSet* fvars) {
    switch (e->type) {
        case EXPR_INTEGER:
            fmpz_mpoly_q_set_si(out, e->data.integer, mctx);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            fmpz_mpoly_q_set_fmpz(out, z, mctx);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_SYMBOL: {
            int idx = var_index(fvars, e->data.symbol);
            if (idx < 0) return 0;
            fmpz_mpoly_q_gen(out, (slong)idx, mctx);
            return 1;
        }
        case EXPR_FUNCTION: {
            const char* h = fn_head_name(e);
            if (!h) return 0;
            size_t n = e->data.function.arg_count;
            if (strcmp(h, "Rational") == 0) {
                if (n != 2) return 0;
                fmpq_t c; fmpq_init(c);
                if (!fmpz_from_int_expr(fmpq_numref(c), e->data.function.args[0]) ||
                    !fmpz_from_int_expr(fmpq_denref(c), e->data.function.args[1])) {
                    fmpq_clear(c); return 0;
                }
                fmpq_canonicalise(c);
                fmpz_mpoly_q_set_fmpq(out, c, mctx);
                fmpq_clear(c);
                return 1;
            }
            if (strcmp(h, "Plus") == 0) {
                fmpz_mpoly_q_zero(out, mctx);
                fmpz_mpoly_q_t t; fmpz_mpoly_q_init(t, mctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = expr_to_mpolyq(e->data.function.args[i], t, mctx, fvars);
                    if (ok) fmpz_mpoly_q_add(out, out, t, mctx);
                }
                fmpz_mpoly_q_clear(t, mctx);
                return ok;
            }
            if (strcmp(h, "Times") == 0) {
                fmpz_mpoly_q_one(out, mctx);
                fmpz_mpoly_q_t t; fmpz_mpoly_q_init(t, mctx);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = expr_to_mpolyq(e->data.function.args[i], t, mctx, fvars);
                    if (ok) fmpz_mpoly_q_mul(out, out, t, mctx);
                }
                fmpz_mpoly_q_clear(t, mctx);
                return ok;
            }
            if (strcmp(h, "Power") == 0) {
                if (n != 2) return 0;
                const Expr* base = e->data.function.args[0];
                const Expr* exp  = e->data.function.args[1];
                if (exp->type != EXPR_INTEGER) return 0;
                long m = (long)exp->data.integer;
                fmpz_mpoly_q_t bq; fmpz_mpoly_q_init(bq, mctx);
                int ok = expr_to_mpolyq(base, bq, mctx, fvars);
                if (ok) {
                    fmpz_mpoly_q_one(out, mctx);
                    long am = m < 0 ? -m : m;
                    for (long i = 0; i < am; i++)
                        fmpz_mpoly_q_mul(out, out, bq, mctx);
                    if (m < 0) fmpz_mpoly_q_inv(out, out, mctx);
                }
                fmpz_mpoly_q_clear(bq, mctx);
                return ok;
            }
            return 0;
        }
        default:
            return 0;
    }
}

/* fmpz_mpoly -> Expr (integer coefficients; variable names from `vs`). */
static Expr* fmpz_mpoly_to_expr(const fmpz_mpoly_t P, const fmpz_mpoly_ctx_t ctx,
                                const VarSet* vs) {
    slong len = fmpz_mpoly_length(P, ctx);
    if (len == 0) return expr_new_integer(0);
    slong nv = (slong)vs->count;
    ulong* exps = malloc(sizeof(ulong) * (size_t)(nv > 0 ? nv : 1));
    Expr** terms = malloc(sizeof(Expr*) * (size_t)len);
    fmpz_t c; fmpz_init(c);
    mpz_t z; mpz_init(z);
    for (slong i = 0; i < len; i++) {
        fmpz_mpoly_get_term_coeff_fmpz(c, P, i, ctx);
        fmpz_mpoly_get_term_exp_ui(exps, P, i, ctx);
        fmpz_get_mpz(z, c);
        Expr** factors = malloc(sizeof(Expr*) * (size_t)(nv + 1));
        size_t nf = 0;
        factors[nf++] = expr_from_mpz_local(z);
        for (slong v = 0; v < nv; v++) {
            if (exps[v] == 0) continue;
            Expr* var = expr_new_symbol(vs->names[v]);
            if (exps[v] == 1) factors[nf++] = var;
            else {
                Expr* pa[2] = { var, expr_new_integer((int64_t)exps[v]) };
                factors[nf++] = expr_new_function(expr_new_symbol("Power"), pa, 2);
            }
        }
        terms[i] = (nf == 1) ? factors[0]
                             : expr_new_function(expr_new_symbol("Times"), factors, nf);
        free(factors);
    }
    mpz_clear(z);
    fmpz_clear(c);
    free(exps);
    Expr* r = (len == 1) ? terms[0]
                         : expr_new_function(expr_new_symbol("Plus"), terms, (size_t)len);
    free(terms);
    return r;
}

/* gr_poly (over fmpz_mpoly_q) -> Expr polynomial in `xname`, coefficients
 * rendered as rational functions num/den in the field variables. */
static Expr* grpolyq_to_expr(const gr_poly_t P, gr_ctx_t gctx,
                             const fmpz_mpoly_ctx_t mctx, const VarSet* fvars,
                             const char* xname) {
    slong len = gr_poly_length(P, gctx);
    if (len == 0) return expr_new_integer(0);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)len);
    gr_ptr coeff;
    GR_TMP_INIT(coeff, gctx);
    for (slong i = 0; i < len; i++) {
        GR_IGNORE(gr_poly_get_coeff_scalar(coeff, P, i, gctx));
        fmpz_mpoly_q_struct* q = (fmpz_mpoly_q_struct*)coeff;
        Expr* num = fmpz_mpoly_to_expr(fmpz_mpoly_q_numref(q), mctx, fvars);
        Expr* coeffe;
        if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), mctx)) {
            coeffe = num;
        } else {
            Expr* den = fmpz_mpoly_to_expr(fmpz_mpoly_q_denref(q), mctx, fvars);
            Expr* deninv = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ den, expr_new_integer(-1) }, 2);
            coeffe = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ num, deninv }, 2);
        }
        if (i == 0) {
            terms[i] = coeffe;
        } else {
            Expr* xp;
            if (i == 1) xp = expr_new_symbol(xname);
            else xp = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol(xname), expr_new_integer((int64_t)i) }, 2);
            terms[i] = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ coeffe, xp }, 2);
        }
    }
    GR_TMP_CLEAR(coeff, gctx);
    Expr* r = (len == 1) ? terms[0]
                         : expr_new_function(expr_new_symbol("Plus"), terms, (size_t)len);
    free(terms);
    return r;
}

/* Fill gr_poly A (over gctx) from the Expr polynomial `p` in variable xvar,
 * degree `deg`, coefficients converted through mctx/fvars. Returns 1 on ok. */
static int build_grpolyq(const Expr* p, Expr* xvar, int deg, gr_poly_t A,
                         gr_ctx_t gctx, const fmpz_mpoly_ctx_t mctx,
                         const VarSet* fvars) {
    for (int d = 0; d <= deg; d++) {
        Expr* c = get_coeff((Expr*)p, xvar, d);
        fmpz_mpoly_q_t cq; fmpz_mpoly_q_init(cq, mctx);
        int ok = c ? expr_to_mpolyq(c, cq, mctx, fvars) : 0;
        if (ok) GR_IGNORE(gr_poly_set_coeff_scalar(A, (slong)d, (gr_srcptr)cq, gctx));
        fmpz_mpoly_q_clear(cq, mctx);
        if (c) expr_free(c);
        if (!ok) return 0;
    }
    return 1;
}

/* Readback + sqrt(k) restore + evaluate for one xgcd output polynomial. */
static Expr* fieldpoly_out(const gr_poly_t P, gr_ctx_t gctx,
                           const fmpz_mpoly_ctx_t mctx, const VarSet* fvars,
                           const char* xname, const char* kname) {
    Expr* e = grpolyq_to_expr(P, gctx, mctx, fvars, xname);
    Expr* back = ps_subst_out(e, kname);   /* S -> Sqrt[k] */
    expr_free(e);
    return eval_and_free(back);
}

typedef enum { FIELD_XGCD, FIELD_DIVREM, FIELD_GCD } FieldOp;

/* Shared core: set up Q(params, S)[x] as gr_poly over fmpz_mpoly_q, build a, b,
 * run the requested operation, and render the Mathematica-shaped result:
 *   FIELD_XGCD   -> List[g, List[u, v]]  (u a + v b = g)
 *   FIELD_DIVREM -> List[q, r]           (a = q b + r)
 *   FIELD_GCD    -> g
 * Returns NULL when there is no symbolic radical / var is the radicand /
 * conversion fails. */
static Expr* parametric_field_op(const Expr* a, const Expr* b, const Expr* var,
                                 FieldOp op) {
    if (!a || !b || !var || var->type != EXPR_SYMBOL) return NULL;
    PSDetect st;
    memset(&st, 0, sizeof st);
    st.allow_neg_pow = 1;   /* the field bridge handles rational-function coeffs */
    ps_detect(a, &st);
    ps_detect(b, &st);
    if (st.bad || !st.has_sqrt || !st.kname) return NULL;
    const char* xname = var->data.symbol;
    if (strcmp(xname, st.kname) == 0) return NULL;

    Expr* a2 = ps_subst_in(a, st.kname);
    Expr* b2 = ps_subst_in(b, st.kname);

    /* Field variables = every symbol in a2, b2 except x (the poly variable). */
    VarSet allv; memset(&allv, 0, sizeof allv);
    collect_all_symbols(a2, &allv);
    collect_all_symbols(b2, &allv);
    VarSet fvars; memset(&fvars, 0, sizeof fvars);
    for (size_t i = 0; i < allv.count; i++)
        if (strcmp(allv.names[i], xname) != 0) varset_add(&fvars, allv.names[i]);
    varset_free(&allv);
    if (fvars.count == 0) { varset_free(&fvars); expr_free(a2); expr_free(b2); return NULL; }
    qsort(fvars.names, fvars.count, sizeof(char*), cmp_str);

    Expr* xvar = expr_new_symbol(xname);
    int da = get_degree_poly(a2, xvar);
    int db = get_degree_poly(b2, xvar);

    gr_ctx_t gctx;
    gr_ctx_init_fmpz_mpoly_q(gctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    gr_poly_t A, B, G, S, T;
    gr_poly_init(A, gctx); gr_poly_init(B, gctx); gr_poly_init(G, gctx);
    gr_poly_init(S, gctx); gr_poly_init(T, gctx);

    Expr* result = NULL;
    if (da >= 0 && db >= 0 &&
        build_grpolyq(a2, xvar, da, A, gctx, mctx, &fvars) &&
        build_grpolyq(b2, xvar, db, B, gctx, mctx, &fvars)) {
        if (op == FIELD_XGCD) {
            if (gr_poly_xgcd(G, S, T, A, B, gctx) == GR_SUCCESS) {
                Expr* g = fieldpoly_out(G, gctx, mctx, &fvars, xname, st.kname);
                Expr* u = fieldpoly_out(S, gctx, mctx, &fvars, xname, st.kname);
                Expr* v = fieldpoly_out(T, gctx, mctx, &fvars, xname, st.kname);
                Expr* cof = expr_new_function(expr_new_symbol("List"), (Expr*[]){u, v}, 2);
                result = expr_new_function(expr_new_symbol("List"), (Expr*[]){g, cof}, 2);
            }
        } else if (op == FIELD_DIVREM) {
            /* G := quotient, S := remainder (reuse the buffers). */
            if (gr_poly_divrem(G, S, A, B, gctx) == GR_SUCCESS) {
                Expr* q = fieldpoly_out(G, gctx, mctx, &fvars, xname, st.kname);
                Expr* r = fieldpoly_out(S, gctx, mctx, &fvars, xname, st.kname);
                result = expr_new_function(expr_new_symbol("List"), (Expr*[]){q, r}, 2);
            }
        } else { /* FIELD_GCD */
            if (gr_poly_gcd(G, A, B, gctx) == GR_SUCCESS)
                result = fieldpoly_out(G, gctx, mctx, &fvars, xname, st.kname);
        }
    }

    gr_poly_clear(A, gctx); gr_poly_clear(B, gctx); gr_poly_clear(G, gctx);
    gr_poly_clear(S, gctx); gr_poly_clear(T, gctx);
    fmpz_mpoly_ctx_clear(mctx);
    gr_ctx_clear(gctx);
    expr_free(xvar);
    varset_free(&fvars);
    expr_free(a2);
    expr_free(b2);
    return result;
}

Expr* flint_parametric_field_xgcd(const Expr* a, const Expr* b, const Expr* var) {
    return parametric_field_op(a, b, var, FIELD_XGCD);
}

/* Quotient/remainder over Q(params)(sqrt k): returns List[quotient, remainder]
 * with a = quotient*b + remainder, or NULL when out of scope. */
Expr* flint_parametric_field_divrem(const Expr* a, const Expr* b, const Expr* var) {
    return parametric_field_op(a, b, var, FIELD_DIVREM);
}

/*
 * Cancel/Together a whole rational function e over a parametric radical field
 * Q(t..)(sqrt k) by converting it to an fmpz_mpoly_q element (which stores it in
 * lowest terms num/den automatically) and reading it back. Unlike
 * flint_cancel_fraction this handles a Plus of fractions, nested fractions, and
 * any scalar field expression in one shot — the case the extract-num/den path
 * cannot (it sees den = 1 on a sum and bails, dropping to the slow QA path).
 * Returns the reduced num/den Expr, or NULL when out of scope.
 */
Expr* flint_parametric_field_normalize(const Expr* e) {
    if (!e) return NULL;
    PSDetect st;
    memset(&st, 0, sizeof st);
    st.allow_neg_pow = 1;
    ps_detect(e, &st);
    if (st.bad || !st.has_sqrt || !st.kname) return NULL;

    Expr* e2 = ps_subst_in(e, st.kname);
    VarSet fvars;
    memset(&fvars, 0, sizeof fvars);
    collect_all_symbols(e2, &fvars);
    if (fvars.count == 0) { varset_free(&fvars); expr_free(e2); return NULL; }
    qsort(fvars.names, fvars.count, sizeof(char*), cmp_str);

    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_q_t q;
    fmpz_mpoly_q_init(q, mctx);

    Expr* out = NULL;
    if (expr_to_mpolyq(e2, q, mctx, &fvars)) {
        fmpz_mpoly_q_canonicalise(q, mctx);
        Expr* num = fmpz_mpoly_to_expr(fmpz_mpoly_q_numref(q), mctx, &fvars);
        Expr* r;
        if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), mctx)) {
            r = num;
        } else {
            Expr* den = fmpz_mpoly_to_expr(fmpz_mpoly_q_denref(q), mctx, &fvars);
            Expr* deninv = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ den, expr_new_integer(-1) }, 2);
            r = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ num, deninv }, 2);
        }
        out = eval_and_free(ps_subst_out(r, st.kname));
        expr_free(r);
    }

    fmpz_mpoly_q_clear(q, mctx);
    fmpz_mpoly_ctx_clear(mctx);
    varset_free(&fvars);
    expr_free(e2);
    return out;
}

/* Algebraic-extension GCD only: Q(sqrt d) -> Q(zeta_n) -> radical tower. Returns
 * NULL when no algebraic generator is present (plain Q[x] / parametric), so a
 * caller can keep its existing path for those. This is the entry consumers like
 * PolynomialGCD use — it engages exactly where the classical path under-reduces. */
Expr* flint_extension_gcd(const Expr* a, const Expr* b) {
    Expr* r = flint_numberfield_gcd(a, b);
    if (r) return r;
    r = flint_cyclotomic_gcd(a, b);
    if (r) return r;
    r = flint_tower_gcd(a, b);
    if (r) return r;
    return flint_parametric_sqrt_gcd(a, b);   /* Q(t..)(sqrt k), k a symbol */
}

/* Full chain: plain rational multivariate first, then the extension fields. */
Expr* flint_polynomial_gcd(const Expr* a, const Expr* b) {
    Expr* r = flint_multivariate_gcd(a, b);
    if (r) return r;
    return flint_extension_gcd(a, b);
}

/* Exact division a/b over the detected extension field (Q(sqrt d), Q(zeta_n),
 * radical tower). Returns NULL when no algebraic generator is present or when b
 * does not divide a exactly — so a consumer (Cancel's divide-back) can fall back
 * to its classical path. */
Expr* flint_extension_divexact(const Expr* a, const Expr* b) {
    Expr* r = numberfield_op(a, b, ABSOP_DIVEXACT);
    if (r) return r;
    r = cyclotomic_op(a, b, ABSOP_DIVEXACT);
    if (r) return r;
    r = tower_op(a, b, ABSOP_DIVEXACT);
    if (r) return r;
    return parametric_sqrt_op(a, b, ABSOP_DIVEXACT);   /* Q(t..)(sqrt k) */
}

/* ================================================================== */
/*  Plain multivariate Q[x..] public wrappers (resultant / factor).    */
/*  GCD is already public (flint_multivariate_gcd/flint_polynomial_gcd).*/
/* ================================================================== */

Expr* flint_polynomial_resultant(const Expr* a, const Expr* b, const Expr* var) {
    if (!var || var->type != EXPR_SYMBOL) return NULL;
    return flint_multivariate_resultant(a, b, var->data.symbol);
}

Expr* flint_polynomial_factor(const Expr* p) {
    return flint_multivariate_factor_impl(p, 0);
}

Expr* flint_polynomial_factor_squarefree(const Expr* p) {
    return flint_multivariate_factor_impl(p, 1);
}

/* fmpz -> Expr (Integer or BigInt, normalised). */
static Expr* fz_to_expr(const fmpz_t z) {
    mpz_t m; mpz_init(m);
    fmpz_get_mpz(m, z);
    Expr* e = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
    mpz_clear(m);
    return e;
}

/* Build the (unevaluated) monomial term  coeff * var^deg  (deg >= 0). The
 * coeff==1 / deg in {0,1} simplifications are left to the caller's eval pass. */
static Expr* uni_term(const fmpz_t coeff, const Expr* var, int deg) {
    Expr* c = fz_to_expr(coeff);
    if (deg == 0) return c;
    Expr* pw = (deg == 1)
        ? expr_copy((Expr*)var)
        : expr_new_function(expr_new_symbol("Power"),
              (Expr*[]){ expr_copy((Expr*)var), expr_new_integer(deg) }, 2);
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ c, pw }, 2);
}

/* One primitive fmpz_poly factor -> (unevaluated) Plus Expr in `var`. */
static Expr* uni_poly_to_expr(const fmpz_poly_t f, const Expr* var) {
    slong d = fmpz_poly_degree(f);
    if (d < 0) return expr_new_integer(0);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(d + 1));
    size_t nt = 0;
    fmpz_t c; fmpz_init(c);
    for (slong j = 0; j <= d; j++) {
        fmpz_poly_get_coeff_fmpz(c, f, j);
        if (fmpz_is_zero(c)) continue;
        terms[nt++] = uni_term(c, var, (int)j);
    }
    fmpz_clear(c);
    Expr* out;
    if (nt == 0)      out = expr_new_integer(0);
    else if (nt == 1) out = terms[0];
    else              out = expr_new_function(expr_new_symbol("Plus"), terms, nt);
    free(terms);
    return out;
}

static Expr* flint_univariate_factor(const Expr* P, const Expr* var) {
    if (!P || !var || var->type != EXPR_SYMBOL) return NULL;
    int deg = get_degree_poly((Expr*)P, (Expr*)var);
    if (deg < 1) return NULL;

    fmpz_poly_t G;
    fmpz_poly_init(G);
    int ok = 1;
    for (int i = 0; i <= deg && ok; i++) {
        Expr* c = get_coeff((Expr*)P, (Expr*)var, i);
        if (c->type == EXPR_INTEGER) {
            fmpz_poly_set_coeff_si(G, i, c->data.integer);
        } else if (c->type == EXPR_BIGINT) {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, c->data.bigint);
            fmpz_poly_set_coeff_fmpz(G, i, z);
            fmpz_clear(z);
        } else {
            ok = 0;   /* Rational / symbolic coefficient -> classical path */
        }
        expr_free(c);
    }

    Expr* out = NULL;
    if (ok && fmpz_poly_degree(G) >= 1) {
        /* fmpz_poly_factor stores the signed content in fac->c and the
         * primitive, positive-leading irreducible factors in fac->p[i] with
         * multiplicities fac->exp[i] — exactly the classical convention. */
        fmpz_poly_factor_t fac;
        fmpz_poly_factor_init(fac);
        fmpz_poly_factor(fac, G);

        size_t cap = (size_t)fac->num + 1;
        Expr** factors = malloc(sizeof(Expr*) * cap);
        size_t nt = 0;
        if (!fmpz_is_one(&fac->c))
            factors[nt++] = fz_to_expr(&fac->c);
        for (slong i = 0; i < fac->num; i++) {
            Expr* be = uni_poly_to_expr(fac->p + i, var);
            slong e = fac->exp[i];
            if (e == 1) {
                factors[nt++] = be;
            } else {
                factors[nt++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ be, expr_new_integer((int64_t)e) }, 2);
            }
        }
        if (nt == 0)      out = expr_new_integer(1);
        else if (nt == 1) out = factors[0];
        else              out = expr_new_function(expr_new_symbol("Times"), factors, nt);
        free(factors);

        fmpz_poly_factor_clear(fac);
    }
    fmpz_poly_clear(G);
    if (out) out = eval_and_free(out);
    return out;
}

Expr* flint_univariate_factor_auto(const Expr* p) {
    if (!p) return NULL;
    /* Detect the variable set with FLINT's own scan (bignum-safe, and it bails
     * on denominators / radicals / out-of-scope heads). We only handle the
     * single-variable case; multivariate stays on the classical path. */
    VarSet vs;
    memset(&vs, 0, sizeof vs);
    if (!collect_vars(p, &vs) || vs.count != 1) { varset_free(&vs); return NULL; }
    Expr* var = expr_new_symbol(vs.names[0]);
    varset_free(&vs);
    Expr* out = flint_univariate_factor(p, var);
    expr_free(var);
    return out;
}

/* ================================================================== */
/*  FLINT` context: direct REPL access to the FLINT-backed kernels.    */
/*  Thin wrappers over the bridge functions; each returns NULL (leaving */
/*  FLINT`f[...] unevaluated) when the argument is out of FLINT's scope.*/
/* ================================================================== */

static Expr* builtin_flint_polynomialgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return flint_polynomial_gcd(res->data.function.args[0],
                                res->data.function.args[1]);
}

static Expr* builtin_flint_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    return flint_polynomial_resultant(res->data.function.args[0],
                                      res->data.function.args[1],
                                      res->data.function.args[2]);
}

static Expr* builtin_flint_factor(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return flint_polynomial_factor(res->data.function.args[0]);
}

static Expr* builtin_flint_factorsquarefree(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return flint_polynomial_factor_squarefree(res->data.function.args[0]);
}

void flint_bridge_init(void) {
    /* FLINT keeps thread-local memory pools (fmpz temporaries, prime caches)
     * that it reuses across calls and only releases at thread/program teardown.
     * Without this, valgrind reports them as "possibly lost" pool noise. Release
     * them at normal exit; safe because no FLINT call outlives main(). */
    atexit(flint_cleanup_master);

    symtab_add_builtin(SYM_FLINT_PolynomialGCD, builtin_flint_polynomialgcd);
    symtab_get_def(SYM_FLINT_PolynomialGCD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialGCD,
        "FLINT`PolynomialGCD[a, b] gives the monic greatest common divisor of "
        "the polynomials a and b over the rationals, computed directly via "
        "FLINT (fmpq_mpoly_gcd). Multivariate. Returns unevaluated if an "
        "argument is not a polynomial over Q.");

    symtab_add_builtin(SYM_FLINT_Resultant, builtin_flint_resultant);
    symtab_get_def(SYM_FLINT_Resultant)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Resultant,
        "FLINT`Resultant[a, b, x] gives the resultant of the polynomials a and "
        "b eliminating the variable x, over the rationals, computed directly "
        "via FLINT (fmpq_mpoly_resultant). Other variables are treated as "
        "coefficients. Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_Factor, builtin_flint_factor);
    symtab_get_def(SYM_FLINT_Factor)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Factor,
        "FLINT`Factor[p] gives the irreducible factorisation of the polynomial "
        "p over the rationals, computed directly via FLINT "
        "(fmpq_mpoly_factor), as Times[const, factor^exp, ...]. Multivariate. "
        "Returns unevaluated if p is not a polynomial over Q.");

    symtab_add_builtin(SYM_FLINT_FactorSquareFree, builtin_flint_factorsquarefree);
    symtab_get_def(SYM_FLINT_FactorSquareFree)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_FactorSquareFree,
        "FLINT`FactorSquareFree[p] gives the squarefree factorisation of the "
        "polynomial p over the rationals, computed directly via FLINT "
        "(fmpq_mpoly_factor_squarefree). Returns unevaluated if out of scope.");
}

#else /* !USE_FLINT */

int   flint_bridge_available(void) { return 0; }
Expr* flint_multivariate_gcd(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }
void  flint_bridge_init(void) { /* no FLINT: nothing to register */ }

#endif /* USE_FLINT */
