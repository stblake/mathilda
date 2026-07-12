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
#include <flint/gr_mat.h>

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

/* Shared core: compute the multivariate GCD of a, b over Q[vars]. When
 * `normalize` is false the raw FLINT (monic) GCD is returned; when true the
 * result is rescaled to the primitive-integer, positive-leading associate
 * that Mathilda's classical PolynomialGCD path produces — i.e. Gauss's
 * lemma: content(gcd) = gcd(content a, content b), pp(gcd) = pp(monic gcd).
 * Returns NULL for numeric or non-polynomial input (caller falls back). */
static Expr* flint_multivariate_gcd_core(const Expr* a, const Expr* b,
                                         int normalize) {
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
        if (fmpq_mpoly_gcd(G, A, B, ctx)) {
            if (normalize && !fmpq_mpoly_is_zero(G, ctx)) {
                /* pp(G): divide out G's own (rational) content, leaving a
                 * primitive-integer polynomial with positive leading
                 * coefficient (FLINT's GCD is monic, so content is positive
                 * and the leading term stays positive after division). */
                fmpq_t cg; fmpq_init(cg);
                fmpq_mpoly_content(cg, G, ctx);
                if (!fmpq_is_zero(cg))
                    fmpq_mpoly_scalar_div_fmpq(G, G, cg, ctx);
                fmpq_clear(cg);

                /* content(gcd) = gcd(content A, content B) — reinstates the
                 * integer content the classical path carries. */
                fmpq_t cA, cB, cc; fmpq_init(cA); fmpq_init(cB); fmpq_init(cc);
                fmpq_mpoly_content(cA, A, ctx);
                fmpq_mpoly_content(cB, B, ctx);
                fmpq_gcd(cc, cA, cB);
                if (!fmpq_is_zero(cc))
                    fmpq_mpoly_scalar_mul_fmpq(G, G, cc, ctx);
                fmpq_clear(cA); fmpq_clear(cB); fmpq_clear(cc);
            }
            out = mpoly_to_expr(G, ctx, &vs);
        }
    }

    fmpq_mpoly_clear(A, ctx);
    fmpq_mpoly_clear(B, ctx);
    fmpq_mpoly_clear(G, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    varset_free(&vs);
    return out;
}

Expr* flint_multivariate_gcd(const Expr* a, const Expr* b) {
    return flint_multivariate_gcd_core(a, b, 0);
}

Expr* flint_multivariate_gcd_normalized(const Expr* a, const Expr* b) {
    return flint_multivariate_gcd_core(a, b, 1);
}

/* Multivariate exact division a / b over Q[x_1..x_n]. Returns the quotient
 * when b divides a exactly, NULL otherwise (or out of scope). Mirrors
 * flint_multivariate_gcd's recognition/conversion; uses fmpq_mpoly_divides. */
Expr* flint_multivariate_divexact(const Expr* a, const Expr* b) {
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
/*  Native univariate-over-Q Risch differential-equation solver        */
/* ================================================================== */
/*
 * Bronstein's base-field RDE core (SPDE ladder + PolyRischDENoCancel) executed
 * entirely in FLINT `fmpq_poly` arithmetic. The Expr-tree implementation in
 * integrate_risch_transcendental.c is faithful but routes every one of the
 * ~O(n) recursion levels' polynomial operations back through the generic
 * evaluator (Expand/Plus/Times/PolynomialGCD/CoefficientList/D as builtins) —
 * on a degree-100 integrand that constant-factor overhead dominates (In17:
 * ~17.5 s, ~95% in evaluate_step/builtin_plus/expr_eq/…, not the poly math).
 * Working in packed fmpq_poly collapses the whole ladder to native FLINT
 * kernels; the Expr path remains the fallback for multivariate / tower base
 * fields (where the coefficients carry other symbols and conversion declines).
 */

/* Expr -> fmpq_poly (univariate in `xname`, rational coefficients). Returns 1
 * and fills `out` on success; 0 if `e` is not such a polynomial (out is left
 * zeroed). Mirrors to_mpoly's coefficient handling. */
static int expr_accum_fmpq_poly(const Expr* e, const char* xname, fmpq_poly_t out) {
    switch (e->type) {
        case EXPR_INTEGER:
            fmpq_poly_set_si(out, e->data.integer);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            fmpq_poly_set_fmpz(out, z);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_SYMBOL:
            if (strcmp(e->data.symbol, xname) != 0) return 0;   /* other symbol */
            fmpq_poly_zero(out);
            fmpq_poly_set_coeff_si(out, 1, 1);                  /* out = x */
            return 1;
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
                fmpq_poly_set_fmpq(out, c);
                fmpq_clear(c);
                return 1;
            }
            if (strcmp(h, "Plus") == 0) {
                fmpq_poly_zero(out);
                fmpq_poly_t t; fmpq_poly_init(t);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = expr_accum_fmpq_poly(e->data.function.args[i], xname, t);
                    if (ok) fmpq_poly_add(out, out, t);
                }
                fmpq_poly_clear(t);
                return ok;
            }
            if (strcmp(h, "Times") == 0) {
                fmpq_poly_one(out);
                fmpq_poly_t t; fmpq_poly_init(t);
                int ok = 1;
                for (size_t i = 0; i < n && ok; i++) {
                    ok = expr_accum_fmpq_poly(e->data.function.args[i], xname, t);
                    if (ok) fmpq_poly_mul(out, out, t);
                }
                fmpq_poly_clear(t);
                return ok;
            }
            if (strcmp(h, "Power") == 0) {
                if (n != 2) return 0;
                const Expr* ex = e->data.function.args[1];
                if (ex->type != EXPR_INTEGER || ex->data.integer < 0) return 0;
                fmpq_poly_t b; fmpq_poly_init(b);
                int ok = expr_accum_fmpq_poly(e->data.function.args[0], xname, b);
                if (ok) fmpq_poly_pow(out, b, (ulong)ex->data.integer);
                fmpq_poly_clear(b);
                return ok;
            }
            return 0;
        }
        default:
            return 0;   /* EXPR_REAL / EXPR_MPFR / … : not a rational polynomial */
    }
}

/* num <- num/gcd(num,den), den <- den/gcd(num,den) (both monic-normalised by
 * fmpq_poly). Leaves a coprime pair in lowest terms. */
static void fq_ratfunc_reduce(fmpq_poly_t num, fmpq_poly_t den) {
    if (fmpq_poly_is_zero(num)) { fmpq_poly_set_si(den, 1); return; }
    fmpq_poly_t g; fmpq_poly_init(g);
    fmpq_poly_gcd(g, num, den);
    if (!fmpq_poly_is_one(g)) {
        fmpq_poly_div(num, num, g);
        fmpq_poly_div(den, den, g);
    }
    fmpq_poly_clear(g);
}

/* Expr rational function -> (num/den) as fmpq_poly, univariate in `xname` over
 * Q. Combines Plus of fractions and Power[·,-k] denominators directly in
 * fmpq_poly — no evaluator Together. Returns 1 on success, 0 if `e` contains a
 * construct outside Q(xname) (another symbol, a fractional/ symbolic power, an
 * inexact number). num, den must be pre-initialised; the result is NOT reduced
 * (caller runs fq_ratfunc_reduce if it wants lowest terms). */
static int expr_to_fmpq_ratfunc(const Expr* e, const char* xname,
                                fmpq_poly_t num, fmpq_poly_t den) {
    if (e->type == EXPR_FUNCTION) {
        const char* h = fn_head_name(e);
        if (h && strcmp(h, "Power") == 0 && e->data.function.arg_count == 2) {
            const Expr* ex = e->data.function.args[1];
            if (ex->type == EXPR_INTEGER && ex->data.integer < 0) {
                fmpq_poly_t b; fmpq_poly_init(b);
                int ok = expr_accum_fmpq_poly(e->data.function.args[0], xname, b);
                if (ok) { fmpq_poly_pow(den, b, (ulong)(-ex->data.integer));
                          fmpq_poly_set_si(num, 1); }
                fmpq_poly_clear(b);
                return ok;
            }
        }
        if (h && strcmp(h, "Times") == 0) {
            fmpq_poly_set_si(num, 1); fmpq_poly_set_si(den, 1);
            fmpq_poly_t tn, td; fmpq_poly_init(tn); fmpq_poly_init(td);
            int ok = 1;
            for (size_t i = 0; i < e->data.function.arg_count && ok; i++) {
                ok = expr_to_fmpq_ratfunc(e->data.function.args[i], xname, tn, td);
                if (ok) { fmpq_poly_mul(num, num, tn); fmpq_poly_mul(den, den, td); }
            }
            fmpq_poly_clear(tn); fmpq_poly_clear(td);
            return ok;
        }
        if (h && strcmp(h, "Plus") == 0) {
            fmpq_poly_zero(num); fmpq_poly_set_si(den, 1);
            fmpq_poly_t tn, td, p1, p2; fmpq_poly_init(tn); fmpq_poly_init(td);
            fmpq_poly_init(p1); fmpq_poly_init(p2);
            int ok = 1;
            for (size_t i = 0; i < e->data.function.arg_count && ok; i++) {
                ok = expr_to_fmpq_ratfunc(e->data.function.args[i], xname, tn, td);
                if (ok) {   /* num/den + tn/td = (num*td + tn*den)/(den*td) */
                    fmpq_poly_mul(p1, num, td);
                    fmpq_poly_mul(p2, tn, den);
                    fmpq_poly_add(num, p1, p2);
                    fmpq_poly_mul(den, den, td);
                }
            }
            fmpq_poly_clear(tn); fmpq_poly_clear(td); fmpq_poly_clear(p1); fmpq_poly_clear(p2);
            return ok;
        }
    }
    /* Polynomial / constant leaf (integer, bigint, rational, symbol, x^k, sums
     * and products thereof): denominator 1. */
    if (!expr_accum_fmpq_poly(e, xname, num)) return 0;
    fmpq_poly_set_si(den, 1);
    return 1;
}

/* fmpq_poly -> Expr, as an unsimplified Plus[Times[coeff, x^k], …] tree; the
 * evaluator canonicalises it after the calling builtin returns. */
static Expr* fmpq_poly_to_expr_x(const fmpq_poly_t p, const char* xname) {
    if (fmpq_poly_is_zero(p)) return expr_new_integer(0);
    slong d = fmpq_poly_degree(p);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(d + 1));
    size_t nt = 0;
    fmpq_t c; fmpq_init(c);
    for (slong k = 0; k <= d; k++) {
        fmpq_poly_get_coeff_fmpq(c, p, k);
        if (fmpq_is_zero(c)) continue;
        Expr* coeff = expr_from_fmpq_local(c);
        Expr* term;
        if (k == 0) {
            term = coeff;
        } else {
            Expr* xk = (k == 1)
                ? expr_new_symbol(xname)
                : expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_new_symbol(xname), expr_new_integer((int64_t)k) }, 2);
            term = expr_new_function(expr_new_symbol("Times"), (Expr*[]){ coeff, xk }, 2);
        }
        terms[nt++] = term;
    }
    fmpq_clear(c);
    Expr* r = (nt == 1)
        ? terms[0]
        : expr_new_function(expr_new_symbol("Plus"), terms, nt);
    free(terms);
    return r;
}

/* SPDE(a, b, c, D=d/dx, n): Bronstein Thm 6.4.1, over fmpq_poly. Fills
 * (b_out, c_out, *m_out, alpha_out, beta_out) and returns 1 on success, or 0
 * for "no solution". All outputs must be pre-initialised by the caller. */
static int fq_spde(const fmpq_poly_t a, const fmpq_poly_t b, const fmpq_poly_t c,
                   long n, fmpq_poly_t b_out, fmpq_poly_t c_out, long* m_out,
                   fmpq_poly_t alpha_out, fmpq_poly_t beta_out) {
    if (n < 0) {
        if (fmpq_poly_is_zero(c)) {
            fmpq_poly_zero(b_out); fmpq_poly_zero(c_out); *m_out = 0;
            fmpq_poly_zero(alpha_out); fmpq_poly_zero(beta_out);
            return 1;
        }
        return 0;                                       /* c != 0, n < 0 */
    }
    fmpq_poly_t g, rem, a1, b1, c1;
    fmpq_poly_init(g); fmpq_poly_init(rem);
    fmpq_poly_init(a1); fmpq_poly_init(b1); fmpq_poly_init(c1);
    fmpq_poly_gcd(g, a, b);
    fmpq_poly_rem(rem, c, g);
    if (!fmpq_poly_is_zero(rem)) {                       /* g does not divide c */
        fmpq_poly_clear(g); fmpq_poly_clear(rem);
        fmpq_poly_clear(a1); fmpq_poly_clear(b1); fmpq_poly_clear(c1);
        return 0;
    }
    fmpq_poly_div(a1, a, g);
    fmpq_poly_div(b1, b, g);
    fmpq_poly_div(c1, c, g);
    fmpq_poly_clear(g); fmpq_poly_clear(rem);
    long da = (long)fmpq_poly_degree(a1);
    if (da == 0) {                                      /* a1 in k*: base case */
        fmpq_poly_div(b_out, b1, a1);                   /* b1 / a1 (scalar) */
        fmpq_poly_div(c_out, c1, a1);                   /* c1 / a1 (scalar) */
        *m_out = n;                                     /* SPDE returns (…, n, 1, 0) */
        fmpq_poly_set_si(alpha_out, 1);
        fmpq_poly_zero(beta_out);
        fmpq_poly_clear(a1); fmpq_poly_clear(b1); fmpq_poly_clear(c1);
        return 1;
    }
    /* xgcd: S*b1 + T*a1 = 1 (gcd(a1,b1)=1); s = S. r = (s c1) mod a1. */
    fmpq_poly_t G2, S, T, sc, r;
    fmpq_poly_init(G2); fmpq_poly_init(S); fmpq_poly_init(T);
    fmpq_poly_init(sc); fmpq_poly_init(r);
    fmpq_poly_xgcd(G2, S, T, b1, a1);
    fmpq_poly_mul(sc, S, c1);
    fmpq_poly_rem(r, sc, a1);
    fmpq_poly_clear(G2); fmpq_poly_clear(S); fmpq_poly_clear(T); fmpq_poly_clear(sc);
    /* z = (c1 - b1 r) / a1 (exact). */
    fmpq_poly_t b1r, z, da1, b2, dr, c2;
    fmpq_poly_init(b1r); fmpq_poly_init(z);
    fmpq_poly_mul(b1r, b1, r);
    fmpq_poly_sub(z, c1, b1r);
    fmpq_poly_div(z, z, a1);
    fmpq_poly_clear(b1r);
    /* recurse SPDE(a1, b1 + Da1, z - Dr, n - da). */
    fmpq_poly_init(da1); fmpq_poly_init(b2); fmpq_poly_init(dr); fmpq_poly_init(c2);
    fmpq_poly_derivative(da1, a1);
    fmpq_poly_add(b2, b1, da1);
    fmpq_poly_derivative(dr, r);
    fmpq_poly_sub(c2, z, dr);
    fmpq_poly_clear(b1); fmpq_poly_clear(c1); fmpq_poly_clear(da1);
    fmpq_poly_clear(z); fmpq_poly_clear(dr);

    fmpq_poly_t sub_b, sub_c, sub_alpha, sub_beta;
    fmpq_poly_init(sub_b); fmpq_poly_init(sub_c);
    fmpq_poly_init(sub_alpha); fmpq_poly_init(sub_beta);
    long sub_m = 0;
    int rc = fq_spde(a1, b2, c2, n - da, sub_b, sub_c, &sub_m, sub_alpha, sub_beta);
    fmpq_poly_clear(b2); fmpq_poly_clear(c2);
    if (!rc) {
        fmpq_poly_clear(a1); fmpq_poly_clear(r);
        fmpq_poly_clear(sub_b); fmpq_poly_clear(sub_c);
        fmpq_poly_clear(sub_alpha); fmpq_poly_clear(sub_beta);
        return 0;
    }
    /* adopt inner eqn; alpha = a1*sub.alpha, beta = a1*sub.beta + r. */
    fmpq_poly_set(b_out, sub_b);
    fmpq_poly_set(c_out, sub_c);
    *m_out = sub_m;
    fmpq_poly_mul(alpha_out, a1, sub_alpha);
    fmpq_poly_mul(beta_out, a1, sub_beta);
    fmpq_poly_add(beta_out, beta_out, r);
    fmpq_poly_clear(a1); fmpq_poly_clear(r);
    fmpq_poly_clear(sub_b); fmpq_poly_clear(sub_c);
    fmpq_poly_clear(sub_alpha); fmpq_poly_clear(sub_beta);
    return 1;
}

/* PolyRischDENoCancel1(b, c, D=d/dx, n), b != 0 (Bronstein p.208): leading
 * terms of Dq and bq never cancel, so q is built top-down one monomial per
 * pass. Fills q_out and returns 1, or returns 0 for "no solution". */
static int fq_polyrischde_nocancel1(const fmpq_poly_t b, const fmpq_poly_t c,
                                    long n, fmpq_poly_t q_out) {
    long db = (long)fmpq_poly_degree(b);
    fmpq_t lcb; fmpq_init(lcb);
    fmpq_poly_get_coeff_fmpq(lcb, b, db);
    fmpq_poly_zero(q_out);
    fmpq_poly_t cc, p, dp, bp;
    fmpq_poly_init(cc); fmpq_poly_init(p); fmpq_poly_init(dp); fmpq_poly_init(bp);
    fmpq_poly_set(cc, c);
    int ok = 1;
    fmpq_t lcc, coeff; fmpq_init(lcc); fmpq_init(coeff);
    while (!fmpq_poly_is_zero(cc)) {
        long dc = (long)fmpq_poly_degree(cc);
        long m = dc - db;
        if (n < 0 || m < 0 || m > n) { ok = 0; break; }
        fmpq_poly_get_coeff_fmpq(lcc, cc, dc);
        fmpq_div(coeff, lcc, lcb);                      /* lc(c)/lc(b) */
        fmpq_poly_zero(p);
        fmpq_poly_set_coeff_fmpq(p, m, coeff);          /* p = coeff x^m */
        fmpq_poly_add(q_out, q_out, p);
        fmpq_poly_derivative(dp, p);
        fmpq_poly_mul(bp, b, p);
        fmpq_poly_sub(cc, cc, dp);                      /* c <- c - Dp - b p */
        fmpq_poly_sub(cc, cc, bp);
        n = m - 1;
    }
    fmpq_clear(lcc); fmpq_clear(coeff); fmpq_clear(lcb);
    fmpq_poly_clear(cc); fmpq_poly_clear(p); fmpq_poly_clear(dp); fmpq_poly_clear(bp);
    return ok;
}

/* PolyRischDENoCancel2(b=0, c, D=d/dx, n): the equation is Dq = c, i.e. a
 * bounded-degree polynomial antiderivative of c. Fills q_out, returns 1, or 0
 * when no such q exists (lambda = 1, delta = 0 for the base derivation). */
static int fq_polyrischde_integrate(const fmpq_poly_t c, long n, fmpq_poly_t q_out) {
    fmpq_poly_zero(q_out);
    fmpq_poly_t cc, p, dp;
    fmpq_poly_init(cc); fmpq_poly_init(p); fmpq_poly_init(dp);
    fmpq_poly_set(cc, c);
    int ok = 1;
    fmpq_t lcc, coeff, mq; fmpq_init(lcc); fmpq_init(coeff); fmpq_init(mq);
    while (!fmpq_poly_is_zero(cc)) {
        long dc = (long)fmpq_poly_degree(cc);
        long m = dc + 1;                                /* deg(c) - delta + 1 */
        if (n < 0 || m < 0 || m > n) { ok = 0; break; }
        fmpq_poly_get_coeff_fmpq(lcc, cc, dc);
        fmpq_set_si(mq, m, 1);
        fmpq_div(coeff, lcc, mq);                       /* lc(c)/m */
        fmpq_poly_zero(p);
        fmpq_poly_set_coeff_fmpq(p, m, coeff);
        fmpq_poly_add(q_out, q_out, p);
        fmpq_poly_derivative(dp, p);
        fmpq_poly_sub(cc, cc, dp);
        n = m - 1;
    }
    fmpq_clear(lcc); fmpq_clear(coeff); fmpq_clear(mq);
    fmpq_poly_clear(cc); fmpq_poly_clear(p); fmpq_poly_clear(dp);
    return ok;
}

/* Solve a Dq + b q = c for q in Q[x], deg(q) <= n, via SPDE + PolyRischDE.
 * Fills q_out (fmpq_poly) and returns 1, or returns 0 for "no solution". */
static int fq_rde_solve_ladder(const fmpq_poly_t a, const fmpq_poly_t b,
                               const fmpq_poly_t c, long n, fmpq_poly_t q_out) {
    fmpq_poly_t sb, sc, salpha, sbeta;
    fmpq_poly_init(sb); fmpq_poly_init(sc);
    fmpq_poly_init(salpha); fmpq_poly_init(sbeta);
    long m = 0;
    int result = 0;
    if (fq_spde(a, b, c, n, sb, sc, &m, salpha, sbeta)) {
        fmpq_poly_t H;
        fmpq_poly_init(H);
        int hok = fmpq_poly_is_zero(sb)
            ? fq_polyrischde_integrate(sc, m, H)
            : fq_polyrischde_nocancel1(sb, sc, m, H);
        if (hok) {
            fmpq_poly_mul(q_out, salpha, H);            /* q = alpha H + beta */
            fmpq_poly_add(q_out, q_out, sbeta);
            result = 1;
        }
        fmpq_poly_clear(H);
    }
    fmpq_poly_clear(sb); fmpq_poly_clear(sc);
    fmpq_poly_clear(salpha); fmpq_poly_clear(sbeta);
    return result;
}

/* Native base-field RDE for the exponential tower, univariate over Q.  Solves
 *   Dq + f q = g   for y = q  (D = d/d(xvar), the base derivation), where f is
 * a POLYNOMIAL (f = i·u' for an exp monomial e^u, u a polynomial in xvar, so
 * den(f) = 1 and WeakNormalizer is a no-op, w = 1) and g is a rational function
 * over Q(xvar).  Runs RdeNormalDenominator + RdeBoundDegreeBase + the SPDE
 * ladder + PolyRischDE entirely in fmpq_poly, converting f and g straight from
 * Expr to fmpq_poly (no evaluator Together/Cancel), which is where the seconds
 * of Expr rational arithmetic on a degree-2n denominator are eliminated.
 * Returns:
 *    1  solved — *y_out = the (reduced) rational solution y as an Expr,
 *    0  no solution of bounded degree exists (authoritative decline),
 *   -1  out of scope — g is not univariate over Q, or f is not a polynomial
 *       (a genuinely rational coefficient, i.e. the primitive/log tower where
 *       weak normalization is non-trivial) — caller falls back to Expr. */
int flint_rde_base_solve_fg(const Expr* f, const Expr* g,
                            const char* xvar, Expr** y_out) {
    if (y_out) *y_out = NULL;
    fmpq_poly_t FN, FD, NG, EN;
    fmpq_poly_init(FN); fmpq_poly_init(FD); fmpq_poly_init(NG); fmpq_poly_init(EN);
    if (!expr_to_fmpq_ratfunc(f, xvar, FN, FD) ||
        !expr_to_fmpq_ratfunc(g, xvar, NG, EN)) {
        fmpq_poly_clear(FN); fmpq_poly_clear(FD); fmpq_poly_clear(NG); fmpq_poly_clear(EN);
        return -1;                                      /* not univariate over Q */
    }
    fq_ratfunc_reduce(FN, FD);
    if (!fmpq_poly_is_one(FD)) {                         /* f rational: weak-norm */
        fmpq_poly_clear(FN); fmpq_poly_clear(FD); fmpq_poly_clear(NG); fmpq_poly_clear(EN);
        return -1;                                      /* case not handled here */
    }
    fq_ratfunc_reduce(NG, EN);
    /* NF = f (polynomial), dn = 1, num_g = NG, en = EN.  RdeNormalDenominator
     * with dn = 1:  p = gcd(1, en) = 1,  h = gcd(en, en'). */
    fmpq_poly_t H, enp;
    fmpq_poly_init(H); fmpq_poly_init(enp);
    fmpq_poly_derivative(enp, EN);
    fmpq_poly_gcd(H, EN, enp);                           /* h = gcd(en, en') */
    fmpq_poly_clear(enp);
    /* Guard: en | h^2 (dn = 1). */
    fmpq_poly_t H2, grem;
    fmpq_poly_init(H2); fmpq_poly_init(grem);
    fmpq_poly_mul(H2, H, H);
    fmpq_poly_rem(grem, H2, EN);
    int result;
    if (!fmpq_poly_is_zero(grem)) {
        result = 0;                                     /* en ∤ h^2 */
    } else {
        /* a = h,  b = h·f - Dh,  c = (h^2 / en) num_g. */
        fmpq_poly_t aa, bb, cc, dh, hnf, h2_en;
        fmpq_poly_init(aa); fmpq_poly_init(bb); fmpq_poly_init(cc);
        fmpq_poly_init(dh); fmpq_poly_init(hnf); fmpq_poly_init(h2_en);
        fmpq_poly_set(aa, H);
        fmpq_poly_derivative(dh, H);
        fmpq_poly_mul(hnf, H, FN);
        fmpq_poly_sub(bb, hnf, dh);
        fmpq_poly_div(h2_en, H2, EN);                   /* exact by the guard */
        fmpq_poly_mul(cc, h2_en, NG);
        fmpq_poly_clear(dh); fmpq_poly_clear(hnf); fmpq_poly_clear(h2_en);
        /* RdeBoundDegreeBase (Bronstein p.199), including the b/a resonance. */
        long da = (long)fmpq_poly_degree(aa);
        long db = (long)fmpq_poly_degree(bb);
        long dc = (long)fmpq_poly_degree(cc);
        long mx = (db > da - 1) ? db : (da - 1);
        long n = dc - mx; if (n < 0) n = 0;
        if (db == da - 1 && da >= 1) {
            fmpq_t lcb, lca, mm; fmpq_init(lcb); fmpq_init(lca); fmpq_init(mm);
            fmpq_poly_get_coeff_fmpq(lcb, bb, db);
            fmpq_poly_get_coeff_fmpq(lca, aa, da);
            fmpq_div(mm, lcb, lca);
            fmpq_neg(mm, mm);                           /* -lc(b)/lc(a) */
            if (fmpz_is_one(fmpq_denref(mm)) && fmpz_fits_si(fmpq_numref(mm))) {
                long mv = fmpz_get_si(fmpq_numref(mm));
                long cand = dc - db; if (cand < 0) cand = 0;
                if (mv > n) n = mv;
                if (cand > n) n = cand;
            }
            fmpq_clear(lcb); fmpq_clear(lca); fmpq_clear(mm);
        }
        fmpq_poly_t q;
        fmpq_poly_init(q);
        if (fq_rde_solve_ladder(aa, bb, cc, n, q)) {
            /* y = q / h (w = 1), reduced to lowest terms in fmpq_poly. */
            fq_ratfunc_reduce(q, H);
            if (y_out) {
                Expr* yn = fmpq_poly_to_expr_x(q, xvar);
                if (fmpq_poly_is_one(H)) {
                    *y_out = yn;
                } else {
                    Expr* yd = fmpq_poly_to_expr_x(H, xvar);
                    *y_out = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ yn, expr_new_function(expr_new_symbol("Power"),
                            (Expr*[]){ yd, expr_new_integer(-1) }, 2) }, 2);
                }
            }
            result = 1;
        } else {
            result = 0;
        }
        fmpq_poly_clear(q);
        fmpq_poly_clear(aa); fmpq_poly_clear(bb); fmpq_poly_clear(cc);
    }
    fmpq_poly_clear(H2); fmpq_poly_clear(grem); fmpq_poly_clear(H);
    fmpq_poly_clear(FN); fmpq_poly_clear(FD); fmpq_poly_clear(NG); fmpq_poly_clear(EN);
    return result;
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

/* True if `e` carries an actual denominator to combine — a Power[base, k] with
 * k a negative integer (or negative-numerator Rational) and base a symbol or
 * compound (not a numeric literal, whose negative power is just a rational
 * constant). Only such inputs are combined by Together; a denominator-free
 * product/polynomial is left factored, so flint_rational_together declines it. */
static int expr_has_denominator(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    if (fn_head_name(e) && strcmp(fn_head_name(e), "Power") == 0 &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* ex   = e->data.function.args[1];
        int neg = (ex->type == EXPR_INTEGER && ex->data.integer < 0);
        if (!neg && ex->type == EXPR_FUNCTION && fn_head_name(ex) &&
            strcmp(fn_head_name(ex), "Rational") == 0 && ex->data.function.arg_count == 2 &&
            ex->data.function.args[0]->type == EXPR_INTEGER &&
            ex->data.function.args[0]->data.integer < 0)
            neg = 1;
        if (neg && (base->type == EXPR_SYMBOL || base->type == EXPR_FUNCTION))
            return 1;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_has_denominator(e->data.function.args[i])) return 1;
    return 0;
}

/*
 * Together for a plain rational function over Q: combine into a single fraction
 * in lowest terms via fmpz_mpoly_q (which stores num/den reduced automatically)
 * and read it back. This is the fast path for the common case the algebraic /
 * parametric normalizers above decline — an ordinary rational in symbols, whose
 * classical together_recursive expands and GCDs a degree-2n denominator at
 * O(seconds). Returns the reduced num/den Expr, or NULL when:
 *   - the input has no denominator to combine (Together leaves a product /
 *     polynomial factored — we must NOT expand it), or
 *   - the input is not a plain rational over Q (a symbolic/fractional power such
 *     as a^x or Sqrt[2], or a transcendental kernel — expr_to_mpolyq declines),
 * in which case the caller keeps its classical path. Output form matches the
 * classical Together (expanded, reduced num/den).
 */
Expr* flint_rational_together(const Expr* e) {
    if (!e || !expr_has_denominator(e)) return NULL;
    VarSet fvars;
    memset(&fvars, 0, sizeof fvars);
    collect_all_symbols(e, &fvars);
    if (fvars.count == 0) { varset_free(&fvars); return NULL; }
    qsort(fvars.names, fvars.count, sizeof(char*), cmp_str);

    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_q_t q;
    fmpz_mpoly_q_init(q, mctx);

    Expr* out = NULL;
    if (expr_to_mpolyq(e, q, mctx, &fvars)) {
        fmpz_mpoly_q_canonicalise(q, mctx);
        Expr* num = fmpz_mpoly_to_expr(fmpz_mpoly_q_numref(q), mctx, &fvars);
        if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), mctx)) {
            out = eval_and_free(num);
        } else {
            Expr* den = fmpz_mpoly_to_expr(fmpz_mpoly_q_denref(q), mctx, &fvars);
            Expr* deninv = expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ den, expr_new_integer(-1) }, 2);
            Expr* r = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ num, deninv }, 2);
            out = eval_and_free(r);
        }
    }

    fmpz_mpoly_q_clear(q, mctx);
    fmpz_mpoly_ctx_clear(mctx);
    varset_free(&fvars);
    return out;
}

/* ================================================================== */
/*  Genuine-algebraic parametric tower: reduction over Q(params)[gens] */
/*                                                                     */
/*  Handles rational functions over the field                          */
/*      K = Q(t_1..t_p)(alpha_1, ..., alpha_r)                         */
/*  where each generator alpha_l is a *genuine* algebraic element that  */
/*  the sqrt-only parametric path and the constant-radicand            */
/*  number-field path both reject:                                     */
/*    - a radical Power[base, p/q] of ANY index q >= 2 (cube roots      */
/*      included) whose radicand `base` carries a free symbol (a        */
/*      polynomial / symbol radicand, e.g. (x(1-x)(1-k x))^(1/3),       */
/*      k^(1/3)); and                                                   */
/*    - a root of unity (-1)^(p/q).                                     */
/*                                                                     */
/*  K is modelled as the quotient ring  Q[params, gens] / I  with       */
/*      I = ( gen_l^{Q_l} - radicand_l,  Phi_{2Q}(omega),  ... ).       */
/*  Each minimal polynomial is monic in its own, distinct generator, so */
/*  ordering the generators as the leading LEX variables makes the set  */
/*  a Groebner basis (Buchberger's first criterion: pairwise-coprime    */
/*  leading monomials), and fmpz_mpoly_divrem_ideal returns the         */
/*  canonical normal form. An expression E, written N/D over            */
/*  Q[params,gens], is zero in K iff N reduces to 0 modulo I -- a        */
/*  sound rigorous zero test (also complete when the minimal polynomials */
/*  are irreducible, which holds for a non-perfect-power radicand and    */
/*  for Phi_N). No inversion / primes / CRT are needed for this         */
/*  reduction; it is what lets Together / Cancel / Simplify collapse a   */
/*  verified-zero algebraic identity (e.g. D[Integrate[f],x]-f) to 0    */
/*  by construction, without ever consulting a numeric zero oracle.     */
/* ================================================================== */

#define GENALG_MAXGEN 8

typedef struct {
    int   is_rou;      /* root of unity (base == -1): reduce mod Phi_{2Q}    */
    long  Q;           /* index: gen = base^(1/Q), gen^Q = base              */
    Expr* base_exp;    /* Expand[radicand], owned (NULL for a root of unity) */
    int   genuine;     /* radicand carries a free symbol (the new regime)    */
    char  name[24];    /* fresh generator symbol, never leaks into output    */
} GAGen;

typedef struct {
    GAGen gens[GENALG_MAXGEN];
    int   ngen;
    int   bad;         /* out-of-scope construct (inexact real, > MAXGEN, …) */
} GADetect;

static void ga_free(GADetect* st) {
    for (int i = 0; i < st->ngen; i++)
        if (st->gens[i].base_exp) { expr_free(st->gens[i].base_exp); st->gens[i].base_exp = NULL; }
}

/* Expand[base] -> canonical radicand, so two printed forms of the same
 * polynomial (e.g. x(1-x)(1-k x) and x(-1+x)(-1+k x)) map to one generator. */
static Expr* ga_expand(const Expr* base) {
    Expr* c = expr_copy((Expr*)base);
    Expr* args[1] = { c };
    Expr* call = expr_new_function(expr_new_symbol("Expand"), args, 1);
    return eval_and_free(call);
}

/* Find an existing generator (root of unity: the single omega; algebraic: same
 * expanded radicand) and merge its index via lcm, or add a new one. Takes
 * ownership of base_exp (freed on merge/failure). Returns the generator index. */
static int ga_find_or_add(GADetect* st, int is_rou, long q, Expr* base_exp) {
    for (int i = 0; i < st->ngen; i++) {
        if (is_rou && st->gens[i].is_rou) {
            st->gens[i].Q = lcm_l(st->gens[i].Q, q);
            if (base_exp) expr_free(base_exp);
            return i;
        }
        if (!is_rou && !st->gens[i].is_rou && base_exp && st->gens[i].base_exp &&
            expr_eq(st->gens[i].base_exp, base_exp)) {
            st->gens[i].Q = lcm_l(st->gens[i].Q, q);
            expr_free(base_exp);
            return i;
        }
    }
    if (st->ngen >= GENALG_MAXGEN) { st->bad = 1; if (base_exp) expr_free(base_exp); return -1; }
    int i = st->ngen++;
    st->gens[i].is_rou   = is_rou;
    st->gens[i].Q        = q;
    st->gens[i].base_exp = base_exp;
    st->gens[i].genuine  = 0;
    snprintf(st->gens[i].name, sizeof st->gens[i].name, "$galg%d$", i);
    if (!is_rou && base_exp) {
        VarSet tmp; memset(&tmp, 0, sizeof tmp);
        collect_all_symbols(base_exp, &tmp);
        st->gens[i].genuine = (tmp.count > 0);
        varset_free(&tmp);
    }
    return i;
}

/* Walk `e` collecting the algebraic generators and marking out-of-scope forms. */
static void ga_walk(const Expr* e, GADetect* st) {
    if (st->bad || !e) return;
    if (e->type == EXPR_REAL || e->type == EXPR_MPFR) { st->bad = 1; return; }
    if (e->type != EXPR_FUNCTION) return;   /* integers / symbols / strings: fine */
    const char* h = fn_head_name(e);
    if (h && strcmp(h, "Power") == 0 && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        long p, q;
        if (ps_as_ratio(exp, &p, &q) && q >= 2) {
            if (base->type == EXPR_INTEGER && base->data.integer == -1) {
                ga_find_or_add(st, 1, q, NULL);          /* root of unity */
                return;                                  /* constant base */
            }
            Expr* be = ga_expand(base);
            if (!be) { st->bad = 1; return; }
            ga_find_or_add(st, 0, q, be);                /* takes ownership of be */
            ga_walk(base, st);                           /* nested radicals */
            return;
        }
    }
    ga_walk(e->data.function.head, st);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        ga_walk(e->data.function.args[i], st);
}

/* Cheap pre-scan: any radical Power[_, Rational[p, q>=2]] present? */
static int ga_has_candidate(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 0;
    const char* h = fn_head_name(e);
    if (h && strcmp(h, "Power") == 0 && e->data.function.arg_count == 2) {
        long p, q;
        if (ps_as_ratio(e->data.function.args[1], &p, &q) && q >= 2) return 1;
    }
    if (ga_has_candidate(e->data.function.head)) return 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ga_has_candidate(e->data.function.args[i])) return 1;
    return 0;
}

/* Replace each radical / root of unity by the matching generator raised to the
 * corresponding power: Power[base, p/q] -> gen^(p * Q/q), (-1)^(p/q) -> omega^…. */
static Expr* ga_subst_in(const Expr* e, const GADetect* st) {
    if (e->type == EXPR_FUNCTION) {
        const char* h = fn_head_name(e);
        if (h && strcmp(h, "Power") == 0 && e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            const Expr* exp  = e->data.function.args[1];
            long p, q;
            if (ps_as_ratio(exp, &p, &q) && q >= 2) {
                int gi = -1;
                if (base->type == EXPR_INTEGER && base->data.integer == -1) {
                    for (int i = 0; i < st->ngen; i++)
                        if (st->gens[i].is_rou) { gi = i; break; }
                } else {
                    Expr* be = ga_expand(base);
                    if (be) {
                        for (int i = 0; i < st->ngen; i++)
                            if (!st->gens[i].is_rou && st->gens[i].base_exp &&
                                expr_eq(st->gens[i].base_exp, be)) { gi = i; break; }
                        expr_free(be);
                    }
                }
                if (gi >= 0) {
                    long ge = p * (st->gens[gi].Q / q);
                    Expr* g = expr_new_symbol(st->gens[gi].name);
                    if (ge == 1) return g;
                    Expr* pa[2] = { g, expr_new_integer(ge) };
                    return expr_new_function(expr_new_symbol("Power"), pa, 2);
                }
            }
        }
        size_t n = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) args[i] = ga_subst_in(e->data.function.args[i], st);
        Expr* head = ga_subst_in(e->data.function.head, st);
        Expr* r = expr_new_function(head, args, n);
        free(args);
        return r;
    }
    return expr_copy((Expr*)e);
}

/* Minimal-polynomial divisor gen^Q - radicand as an fmpz_mpoly (monic in gen).
 * Fails (0) when the radicand is not an integer-coefficient polynomial in the
 * field variables (rational coefficient, nested radical, …) — the caller then
 * falls through to its classical path. */
static int ga_build_radical_relation(fmpz_mpoly_t B, slong genvar, long Q,
                                     const Expr* base_exp,
                                     const fmpz_mpoly_ctx_t mctx, const VarSet* fvars) {
    fmpz_mpoly_q_t bq; fmpz_mpoly_q_init(bq, mctx);
    int ok = expr_to_mpolyq(base_exp, bq, mctx, fvars);
    if (ok) { fmpz_mpoly_q_canonicalise(bq, mctx); ok = fmpz_mpoly_is_one(fmpz_mpoly_q_denref(bq), mctx); }
    if (ok) {
        fmpz_mpoly_t gp; fmpz_mpoly_init(gp, mctx);
        fmpz_mpoly_gen(gp, genvar, mctx);
        fmpz_mpoly_pow_ui(gp, gp, (ulong)Q, mctx);
        fmpz_mpoly_sub(B, gp, fmpz_mpoly_q_numref(bq), mctx);
        fmpz_mpoly_clear(gp, mctx);
    }
    fmpz_mpoly_q_clear(bq, mctx);
    return ok;
}

/* Cyclotomic minimal polynomial Phi_N(gen) as an fmpz_mpoly (monic in gen). */
static int ga_build_cyclotomic_relation(fmpz_mpoly_t B, slong genvar, long N,
                                        const fmpz_mpoly_ctx_t mctx) {
    fmpz_poly_t cy; fmpz_poly_init(cy);
    fmpz_poly_cyclotomic(cy, (ulong)N);
    slong d = fmpz_poly_degree(cy);
    fmpz_mpoly_zero(B, mctx);
    fmpz_mpoly_t gp, term; fmpz_mpoly_init(gp, mctx); fmpz_mpoly_init(term, mctx);
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i <= d; i++) {
        fmpz_poly_get_coeff_fmpz(c, cy, i);
        if (fmpz_is_zero(c)) continue;
        fmpz_mpoly_gen(gp, genvar, mctx);
        fmpz_mpoly_pow_ui(gp, gp, (ulong)i, mctx);
        fmpz_mpoly_scalar_mul_fmpz(term, gp, c, mctx);
        fmpz_mpoly_add(B, B, term, mctx);
    }
    fmpz_clear(c); fmpz_mpoly_clear(gp, mctx); fmpz_mpoly_clear(term, mctx);
    fmpz_poly_clear(cy);
    return 1;
}

/*
 * Reduce a whole rational expression `e` over a genuine-algebraic parametric
 * tower K = Q(params)(alpha_1..alpha_r). Returns the integer 0 when `e` is
 * identically zero in K (the reported Goursat "D[F,x]-integrand" class), else
 * NULL — so Cancel/Together keep their classical path for anything this
 * conservative first increment does not fully resolve, guaranteeing it can only
 * add the ability to reach 0, never worsen existing output. Never mutates `e`.
 */
Expr* flint_algebraic_field_normalize(const Expr* e) {
    if (!e || !ga_has_candidate(e)) return NULL;

    GADetect st; memset(&st, 0, sizeof st);
    ga_walk(e, &st);
    int genuine = 0;
    for (int i = 0; i < st.ngen; i++) if (st.gens[i].genuine) genuine = 1;
    if (st.bad || st.ngen == 0 || !genuine) { ga_free(&st); return NULL; }

    Expr* ebar = ga_subst_in(e, &st);

    /* Variable set: generators first (leading LEX), then the parameters. The
     * parameters are the free symbols of the substituted expression AND of every
     * radicand -- a radicand can carry a parameter (e.g. k in k^(1/3)) that the
     * substitution removes from `ebar` but the relation gen^Q - radicand still
     * needs. */
    VarSet fvars; memset(&fvars, 0, sizeof fvars);
    for (int i = 0; i < st.ngen; i++) varset_add(&fvars, st.gens[i].name);
    VarSet all; memset(&all, 0, sizeof all);
    collect_all_symbols(ebar, &all);
    for (int i = 0; i < st.ngen; i++)
        if (st.gens[i].base_exp) collect_all_symbols(st.gens[i].base_exp, &all);
    VarSet params; memset(&params, 0, sizeof params);
    for (size_t i = 0; i < all.count; i++)
        if (var_index(&fvars, all.names[i]) < 0) varset_add(&params, all.names[i]);
    qsort(params.names, params.count, sizeof(char*), cmp_str);
    for (size_t i = 0; i < params.count; i++) varset_add(&fvars, params.names[i]);
    varset_free(&params);
    varset_free(&all);

    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_q_t q; fmpz_mpoly_q_init(q, mctx);

    Expr* out = NULL;
    if (expr_to_mpolyq(ebar, q, mctx, &fvars)) {
        fmpz_mpoly_q_canonicalise(q, mctx);
        fmpz_mpoly_struct* B = malloc(sizeof(fmpz_mpoly_struct) * (size_t)st.ngen);
        fmpz_mpoly_struct* Qb = malloc(sizeof(fmpz_mpoly_struct) * (size_t)st.ngen);
        for (int i = 0; i < st.ngen; i++) { fmpz_mpoly_init(&B[i], mctx); fmpz_mpoly_init(&Qb[i], mctx); }
        int okB = 1;
        for (int i = 0; i < st.ngen && okB; i++) {
            if (st.gens[i].is_rou)
                okB = ga_build_cyclotomic_relation(&B[i], (slong)i, 2 * st.gens[i].Q, mctx);
            else
                okB = ga_build_radical_relation(&B[i], (slong)i, st.gens[i].Q,
                                                st.gens[i].base_exp, mctx, &fvars);
        }
        if (okB) {
            const fmpz_mpoly_struct* Bset[GENALG_MAXGEN];
            fmpz_mpoly_struct* Qset[GENALG_MAXGEN];
            for (int i = 0; i < st.ngen; i++) { Bset[i] = &B[i]; Qset[i] = &Qb[i]; }
            fmpz_mpoly_t R, RD; fmpz_mpoly_init(R, mctx); fmpz_mpoly_init(RD, mctx);
            fmpz_mpoly_divrem_ideal(Qset, R, fmpz_mpoly_q_numref(q),
                                    (fmpz_mpoly_struct* const*)Bset, st.ngen, mctx);
            /* Reduce the denominator too: report 0 only for a genuine 0 with a
             * non-vanishing denominator, never a 0/0 (a denominator that is a
             * non-zero free-variable polynomial but vanishes in the field). */
            if (fmpz_mpoly_is_zero(R, mctx)) {
                fmpz_mpoly_divrem_ideal(Qset, RD, fmpz_mpoly_q_denref(q),
                                        (fmpz_mpoly_struct* const*)Bset, st.ngen, mctx);
                if (!fmpz_mpoly_is_zero(RD, mctx)) out = expr_new_integer(0);
            }
            fmpz_mpoly_clear(R, mctx); fmpz_mpoly_clear(RD, mctx);
        }
        for (int i = 0; i < st.ngen; i++) { fmpz_mpoly_clear(&B[i], mctx); fmpz_mpoly_clear(&Qb[i], mctx); }
        free(B); free(Qb);
    }

    fmpz_mpoly_q_clear(q, mctx);
    fmpz_mpoly_ctx_clear(mctx);
    varset_free(&fvars);
    expr_free(ebar);
    ga_free(&st);
    return out;
}

/* ================================================================== */
/*  Milestone B: full field arithmetic (inversion / canonical form)   */
/*  over the genuine-algebraic tower K = Q(params)(gen_0..gen_{r-1}).  */
/*                                                                     */
/*  K is a Q(params)-vector space with the monomial basis             */
/*     { gen_0^{e_0} ... gen_{r-1}^{e_{r-1}} : 0 <= e_i < d_i },       */
/*  d_i = deg of gen_i's minimal polynomial (radical: Q_i; root of     */
/*  unity: deg Phi_{2Q_i}). Multiplication by a field element D is a   */
/*  Q(params)-linear map M_D on this space (each product reduced mod   */
/*  the ideal I to its normal form via fmpz_mpoly_divrem_ideal). Then  */
/*      (N / D)  =  M_D^{-1} . coords(N)                               */
/*  solved over the field Q(params) (fraction field of Z[params], the  */
/*  gr fmpz_mpoly_q ring) gives the *rationalised* canonical form of   */
/*  N/D: numerator a reduced polynomial in the generators, denominator */
/*  det(M_D) = Norm_{K/Q(params)}(D), free of every generator. This is */
/*  the engine behind RootReduce; it is rigorous (no numeric oracle)   */
/*  and complete when the minimal polynomials are irreducible.         */
/* ================================================================== */

/* Degree of generator i over Q(params): radical -> Q; root of unity ->
 * deg Phi_{2Q} = eulerphi(2Q). */
static long ga_gen_dim(const GAGen* g) {
    if (!g->is_rou) return g->Q;
    fmpz_poly_t cy; fmpz_poly_init(cy);
    fmpz_poly_cyclotomic(cy, (ulong)(2 * g->Q));
    long d = (long)fmpz_poly_degree(cy);
    fmpz_poly_clear(cy);
    return d < 1 ? 1 : d;
}

/* Radical / root-of-unity Expr for generator i (used when reading the
 * canonical form back out): radical -> base^(1/Q), root of unity -> (-1)^(1/Q). */
static Expr* ga_radical_of(const GAGen* g) {
    Expr* base = g->is_rou ? expr_new_integer(-1) : expr_copy(g->base_exp);
    Expr* ex   = expr_new_function(expr_new_symbol("Rational"),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(g->Q) }, 2);
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ base, ex }, 2);
}

/* Replace every generator symbol ($galgN$) in `e` by its radical form. */
static Expr* ga_subst_back(const Expr* e, const GADetect* st) {
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < st->ngen; i++)
            if (strcmp(e->data.symbol, st->gens[i].name) == 0)
                return ga_radical_of(&st->gens[i]);
        return expr_copy((Expr*)e);
    }
    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) args[i] = ga_subst_back(e->data.function.args[i], st);
        Expr* head = ga_subst_back(e->data.function.head, st);
        Expr* r = expr_new_function(head, args, n);
        free(args);
        return r;
    }
    return expr_copy((Expr*)e);
}

/* Accumulate the reduced polynomial R (already normal form mod I) into a
 * coordinate array `coords[total]` of Z[params] polynomials, indexed by the
 * mixed-radix combination of the generator exponents (via `stride`,`dims`). */
static void ga_extract_coords(const fmpz_mpoly_t R, const fmpz_mpoly_ctx_t mctx,
                              int ngen, int nparam, const long* dims,
                              const long* stride, const fmpz_mpoly_ctx_t pctx,
                              fmpz_mpoly_struct* coords) {
    slong len = fmpz_mpoly_length(R, mctx);
    slong nv = fmpz_mpoly_ctx_nvars(mctx);
    ulong* exp = malloc(sizeof(ulong) * (size_t)(nv > 0 ? nv : 1));
    ulong* pe  = malloc(sizeof(ulong) * (size_t)(nparam > 0 ? nparam : 1));
    fmpz_t c; fmpz_init(c);
    for (slong t = 0; t < len; t++) {
        fmpz_mpoly_get_term_coeff_fmpz(c, R, t, mctx);
        fmpz_mpoly_get_term_exp_ui(exp, R, t, mctx);
        long row = 0; int bad = 0;
        for (int i = 0; i < ngen; i++) {
            if ((long)exp[i] >= dims[i]) { bad = 1; break; }
            row += (long)exp[i] * stride[i];
        }
        if (bad) continue;
        for (int p = 0; p < nparam; p++) pe[p] = exp[ngen + p];
        fmpz_mpoly_t pm; fmpz_mpoly_init(pm, pctx);
        fmpz_mpoly_set_coeff_fmpz_ui(pm, c, pe, pctx);
        fmpz_mpoly_add(&coords[row], &coords[row], pm, pctx);
        fmpz_mpoly_clear(pm, pctx);
    }
    fmpz_clear(c);
    free(exp); free(pe);
}

/*
 * Canonical (rationalised) form of a rational expression `e` over the
 * genuine-algebraic tower K: returns N/D reduced so the denominator is free of
 * every radical/root-of-unity generator (denominator rationalised), else NULL
 * when `e` carries no algebraic generator, is out of scope, is a genuine 0/0, or
 * the multiplication map is singular (a reducible minimal polynomial — outside
 * the completeness guarantee). Never mutates `e`. This is the RootReduce engine.
 */
Expr* flint_algebraic_field_canonical(const Expr* e) {
    if (!e || !ga_has_candidate(e)) return NULL;

    GADetect st; memset(&st, 0, sizeof st);
    ga_walk(e, &st);
    if (st.bad || st.ngen == 0) { ga_free(&st); return NULL; }

    long dims[GENALG_MAXGEN], stride[GENALG_MAXGEN], total = 1;
    for (int i = 0; i < st.ngen; i++) {
        dims[i] = ga_gen_dim(&st.gens[i]);
        if (dims[i] < 1) { ga_free(&st); return NULL; }
    }
    for (int i = 0; i < st.ngen; i++) {
        total *= dims[i];
        if (total > 256) { ga_free(&st); return NULL; }   /* dimension cap */
    }
    stride[st.ngen - 1] = 1;
    for (int i = st.ngen - 2; i >= 0; i--) stride[i] = stride[i + 1] * dims[i + 1];

    Expr* ebar = ga_subst_in(e, &st);

    /* fvars = generators first (leading LEX), then sorted parameters. */
    VarSet fvars; memset(&fvars, 0, sizeof fvars);
    for (int i = 0; i < st.ngen; i++) varset_add(&fvars, st.gens[i].name);
    VarSet all; memset(&all, 0, sizeof all);
    collect_all_symbols(ebar, &all);
    for (int i = 0; i < st.ngen; i++)
        if (st.gens[i].base_exp) collect_all_symbols(st.gens[i].base_exp, &all);
    VarSet params; memset(&params, 0, sizeof params);
    for (size_t i = 0; i < all.count; i++)
        if (var_index(&fvars, all.names[i]) < 0) varset_add(&params, all.names[i]);
    qsort(params.names, params.count, sizeof(char*), cmp_str);
    int nparam = (int)params.count;
    for (int i = 0; i < nparam; i++) varset_add(&fvars, params.names[i]);

    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_ctx_t pctx;
    fmpz_mpoly_ctx_init(pctx, (slong)(nparam > 0 ? nparam : 1), ORD_LEX);

    fmpz_mpoly_q_t q; fmpz_mpoly_q_init(q, mctx);
    Expr* out = NULL;

    if (expr_to_mpolyq(ebar, q, mctx, &fvars)) {
        fmpz_mpoly_q_canonicalise(q, mctx);
        fmpz_mpoly_struct* B  = malloc(sizeof(fmpz_mpoly_struct) * (size_t)st.ngen);
        fmpz_mpoly_struct* Qb = malloc(sizeof(fmpz_mpoly_struct) * (size_t)st.ngen);
        for (int i = 0; i < st.ngen; i++) { fmpz_mpoly_init(&B[i], mctx); fmpz_mpoly_init(&Qb[i], mctx); }
        int okB = 1;
        for (int i = 0; i < st.ngen && okB; i++) {
            if (st.gens[i].is_rou)
                okB = ga_build_cyclotomic_relation(&B[i], (slong)i, 2 * st.gens[i].Q, mctx);
            else
                okB = ga_build_radical_relation(&B[i], (slong)i, st.gens[i].Q,
                                                st.gens[i].base_exp, mctx, &fvars);
        }
        if (okB) {
            const fmpz_mpoly_struct* Bset[GENALG_MAXGEN];
            fmpz_mpoly_struct* Qset[GENALG_MAXGEN];
            for (int i = 0; i < st.ngen; i++) { Bset[i] = &B[i]; Qset[i] = &Qb[i]; }

            fmpz_mpoly_t numR, denR;
            fmpz_mpoly_init(numR, mctx); fmpz_mpoly_init(denR, mctx);
            fmpz_mpoly_divrem_ideal(Qset, numR, fmpz_mpoly_q_numref(q),
                                    (fmpz_mpoly_struct* const*)Bset, st.ngen, mctx);
            fmpz_mpoly_divrem_ideal(Qset, denR, fmpz_mpoly_q_denref(q),
                                    (fmpz_mpoly_struct* const*)Bset, st.ngen, mctx);

            if (fmpz_mpoly_is_zero(denR, mctx)) {
                /* genuine 0/0 in the field: leave alone */
            } else if (fmpz_mpoly_is_zero(numR, mctx)) {
                out = expr_new_integer(0);
            } else {
                int den_has_gen = 0;
                for (int i = 0; i < st.ngen; i++)
                    if (fmpz_mpoly_degree_si(denR, (slong)i, mctx) > 0) { den_has_gen = 1; break; }
                if (!den_has_gen) {
                    /* already rationalised: render numR/denR and substitute radicals back. */
                    Expr* ne = fmpz_mpoly_to_expr(numR, mctx, &fvars);
                    Expr* frac;
                    if (fmpz_mpoly_is_one(denR, mctx)) {
                        frac = ne;
                    } else {
                        Expr* de = fmpz_mpoly_to_expr(denR, mctx, &fvars);
                        frac = expr_new_function(expr_new_symbol("Times"),
                            (Expr*[]){ ne, expr_new_function(expr_new_symbol("Power"),
                                (Expr*[]){ de, expr_new_integer(-1) }, 2) }, 2);
                    }
                    Expr* back = ga_subst_back(frac, &st);
                    expr_free(frac);
                    out = eval_and_free(back);
                } else {
                    /* Invert denR in K: solve M_denR . x = coords(numR). */
                    gr_ctx_t fld; gr_ctx_init_fmpz_mpoly_q(fld, (slong)(nparam > 0 ? nparam : 1), ORD_LEX);
                    gr_mat_t M, X, RHS;
                    gr_mat_init(M,   (slong)total, (slong)total, fld);
                    gr_mat_init(X,   (slong)total, 1, fld);
                    gr_mat_init(RHS, (slong)total, 1, fld);

                    fmpz_mpoly_struct* coords = malloc(sizeof(fmpz_mpoly_struct) * (size_t)total);
                    for (long r = 0; r < total; r++) fmpz_mpoly_init(&coords[r], pctx);

                    /* RHS column = coords(numR) */
                    ga_extract_coords(numR, mctx, st.ngen, nparam, dims, stride, pctx, coords);
                    for (long r = 0; r < total; r++) {
                        fmpz_mpoly_q_struct* ent = (fmpz_mpoly_q_struct*)gr_mat_entry_ptr(RHS, r, 0, fld);
                        fmpz_mpoly_set(fmpz_mpoly_q_numref(ent), &coords[r], pctx);
                        fmpz_mpoly_one(fmpz_mpoly_q_denref(ent), pctx);
                        fmpz_mpoly_zero(&coords[r], pctx);
                    }
                    /* Column j of M = coords(denR * basis_j) reduced mod I. */
                    for (long j = 0; j < total; j++) {
                        fmpz_mpoly_t bj; fmpz_mpoly_init(bj, mctx);
                        fmpz_mpoly_one(bj, mctx);
                        for (int i = 0; i < st.ngen; i++) {
                            long ei = (j / stride[i]) % dims[i];
                            if (ei) {
                                fmpz_mpoly_t gp; fmpz_mpoly_init(gp, mctx);
                                fmpz_mpoly_gen(gp, (slong)i, mctx);
                                fmpz_mpoly_pow_ui(gp, gp, (ulong)ei, mctx);
                                fmpz_mpoly_mul(bj, bj, gp, mctx);
                                fmpz_mpoly_clear(gp, mctx);
                            }
                        }
                        fmpz_mpoly_t prod, R; fmpz_mpoly_init(prod, mctx); fmpz_mpoly_init(R, mctx);
                        fmpz_mpoly_mul(prod, denR, bj, mctx);
                        fmpz_mpoly_divrem_ideal(Qset, R, prod,
                                                (fmpz_mpoly_struct* const*)Bset, st.ngen, mctx);
                        for (long r = 0; r < total; r++) fmpz_mpoly_zero(&coords[r], pctx);
                        ga_extract_coords(R, mctx, st.ngen, nparam, dims, stride, pctx, coords);
                        for (long r = 0; r < total; r++) {
                            fmpz_mpoly_q_struct* ent = (fmpz_mpoly_q_struct*)gr_mat_entry_ptr(M, r, j, fld);
                            fmpz_mpoly_set(fmpz_mpoly_q_numref(ent), &coords[r], pctx);
                            fmpz_mpoly_one(fmpz_mpoly_q_denref(ent), pctx);
                        }
                        fmpz_mpoly_clear(bj, mctx); fmpz_mpoly_clear(prod, mctx); fmpz_mpoly_clear(R, mctx);
                    }

                    if (gr_mat_nonsingular_solve(X, M, RHS, fld) == GR_SUCCESS) {
                        VarSet pvars; memset(&pvars, 0, sizeof pvars);
                        for (int p = 0; p < nparam; p++) varset_add(&pvars, params.names[p]);

                        /* Common denominator Dc = lcm of the solution-coordinate
                         * denominators, so the canonical form comes out as ONE
                         * rationalised fraction (a reduced polynomial in the
                         * radicals over a generator-free denominator) rather than
                         * a sum of fractions. lcm(a,b) = a * (b / gcd(a,b)). */
                        fmpz_mpoly_t Dc, gg, tt;
                        fmpz_mpoly_init(Dc, pctx); fmpz_mpoly_init(gg, pctx); fmpz_mpoly_init(tt, pctx);
                        fmpz_mpoly_one(Dc, pctx);
                        for (long j = 0; j < total; j++) {
                            fmpz_mpoly_q_struct* xj = (fmpz_mpoly_q_struct*)gr_mat_entry_ptr(X, j, 0, fld);
                            if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(xj), pctx)) continue;
                            const fmpz_mpoly_struct* dj = fmpz_mpoly_q_denref(xj);
                            if (fmpz_mpoly_is_one(dj, pctx)) continue;
                            fmpz_mpoly_gcd(gg, Dc, dj, pctx);
                            fmpz_mpoly_divides(tt, dj, gg, pctx);   /* tt = dj / gcd */
                            fmpz_mpoly_mul(Dc, Dc, tt, pctx);       /* Dc = lcm(Dc, dj) */
                        }

                        Expr** terms = malloc(sizeof(Expr*) * (size_t)total);
                        size_t nt = 0;
                        for (long j = 0; j < total; j++) {
                            fmpz_mpoly_q_struct* xj = (fmpz_mpoly_q_struct*)gr_mat_entry_ptr(X, j, 0, fld);
                            if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(xj), pctx)) continue;
                            /* scaled numerator coefficient = num_j * (Dc / den_j), exact. */
                            fmpz_mpoly_t scaled; fmpz_mpoly_init(scaled, pctx);
                            fmpz_mpoly_divides(tt, Dc, fmpz_mpoly_q_denref(xj), pctx);
                            fmpz_mpoly_mul(scaled, fmpz_mpoly_q_numref(xj), tt, pctx);
                            Expr* coeffe = fmpz_mpoly_to_expr(scaled, pctx, &pvars);
                            fmpz_mpoly_clear(scaled, pctx);
                            Expr* mon = NULL;
                            for (int i = 0; i < st.ngen; i++) {
                                long ei = (j / stride[i]) % dims[i];
                                if (!ei) continue;
                                Expr* gp = (ei == 1) ? expr_new_symbol(st.gens[i].name)
                                    : expr_new_function(expr_new_symbol("Power"),
                                        (Expr*[]){ expr_new_symbol(st.gens[i].name), expr_new_integer(ei) }, 2);
                                mon = mon ? expr_new_function(expr_new_symbol("Times"), (Expr*[]){ mon, gp }, 2) : gp;
                            }
                            terms[nt++] = mon ? expr_new_function(expr_new_symbol("Times"),
                                                    (Expr*[]){ coeffe, mon }, 2)
                                              : coeffe;
                        }
                        Expr* numer = (nt == 0) ? expr_new_integer(0)
                                   : (nt == 1) ? terms[0]
                                   : expr_new_function(expr_new_symbol("Plus"), terms, nt);
                        free(terms);
                        Expr* frac;
                        if (fmpz_mpoly_is_one(Dc, pctx)) {
                            frac = numer;
                        } else {
                            Expr* denE = fmpz_mpoly_to_expr(Dc, pctx, &pvars);
                            frac = expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){ numer, expr_new_function(expr_new_symbol("Power"),
                                    (Expr*[]){ denE, expr_new_integer(-1) }, 2) }, 2);
                        }
                        Expr* back = ga_subst_back(frac, &st);
                        expr_free(frac);
                        out = eval_and_free(back);
                        fmpz_mpoly_clear(Dc, pctx); fmpz_mpoly_clear(gg, pctx); fmpz_mpoly_clear(tt, pctx);
                        varset_free(&pvars);
                    }

                    for (long r = 0; r < total; r++) fmpz_mpoly_clear(&coords[r], pctx);
                    free(coords);
                    gr_mat_clear(M, fld); gr_mat_clear(X, fld); gr_mat_clear(RHS, fld);
                    gr_ctx_clear(fld);
                }
            }
            fmpz_mpoly_clear(numR, mctx); fmpz_mpoly_clear(denR, mctx);
        }
        for (int i = 0; i < st.ngen; i++) { fmpz_mpoly_clear(&B[i], mctx); fmpz_mpoly_clear(&Qb[i], mctx); }
        free(B); free(Qb);
    }

    fmpz_mpoly_q_clear(q, mctx);
    fmpz_mpoly_ctx_clear(pctx);
    fmpz_mpoly_ctx_clear(mctx);
    varset_free(&fvars);
    varset_free(&params);
    varset_free(&all);
    expr_free(ebar);
    ga_free(&st);
    return out;
}

/*
 * Together/Cancel a whole expression over the tower by combining it into a single
 * fraction with the radicals treated as *free* kernels (generators) — i.e. the
 * Wolfram-faithful behaviour that keeps radicals in the denominator and does NOT
 * use the minimal-polynomial relations (no rationalisation). Each radical
 * Power[base, p/q] becomes gen^p and the whole expression is combined via
 * fmpz_mpoly_q over Q[params, gens, vars]; the fraction field already stores the
 * result in lowest terms, so common factors — including radical ones exposed by
 * the p/q -> gen^p substitution, e.g. x^2 - k^(2/3) = x^2 - gen^2 = (x-gen)(x+gen)
 * — cancel. This is what combines a Plus of radical fractions that the sqrt-only
 * flint_parametric_field_normalize cannot (cube and higher roots).
 *
 * Scope gate: fires only when every radical generator is genuine (its radicand
 * carries a free symbol, so its powers appear as gen powers and combine
 * correctly); a *constant*-radicand radical (Sqrt[2], 2^(1/3), ...) is rejected,
 * because its relation gen^q = const is lost under the free-kernel treatment and
 * the dedicated number-field path must handle it. Roots of unity are allowed
 * (they combine as kernels, matching WL). Returns NULL out of scope. The caller
 * (flint_cancel_fraction) additionally restricts this to Plus inputs so it never
 * pre-empts the relation-aware single-fraction GCD path. Never mutates `e`.
 */
Expr* flint_algebraic_field_together(const Expr* e) {
    if (!e || !ga_has_candidate(e)) return NULL;

    GADetect st; memset(&st, 0, sizeof st);
    ga_walk(e, &st);
    if (st.bad || st.ngen == 0) { ga_free(&st); return NULL; }
    int genuine = 0;
    for (int i = 0; i < st.ngen; i++) {
        if (st.gens[i].is_rou) continue;
        if (st.gens[i].genuine) genuine = 1;
        else { ga_free(&st); return NULL; }   /* constant radical: defer */
    }
    if (!genuine) { ga_free(&st); return NULL; }

    Expr* ebar = ga_subst_in(e, &st);
    VarSet fvars; memset(&fvars, 0, sizeof fvars);
    collect_all_symbols(ebar, &fvars);                       /* gens + params + vars, all free */
    for (int i = 0; i < st.ngen; i++)
        if (st.gens[i].base_exp) collect_all_symbols(st.gens[i].base_exp, &fvars);
    if (fvars.count == 0) { varset_free(&fvars); expr_free(ebar); ga_free(&st); return NULL; }
    qsort(fvars.names, fvars.count, sizeof(char*), cmp_str);

    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_ctx_init(mctx, (slong)fvars.count, ORD_LEX);
    fmpz_mpoly_q_t q; fmpz_mpoly_q_init(q, mctx);

    Expr* out = NULL;
    if (expr_to_mpolyq(ebar, q, mctx, &fvars)) {
        fmpz_mpoly_q_canonicalise(q, mctx);
        Expr* num = fmpz_mpoly_to_expr(fmpz_mpoly_q_numref(q), mctx, &fvars);
        Expr* r;
        if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(q), mctx)) {
            r = num;
        } else {
            Expr* den = fmpz_mpoly_to_expr(fmpz_mpoly_q_denref(q), mctx, &fvars);
            r = expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ num, expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ den, expr_new_integer(-1) }, 2) }, 2);
        }
        Expr* back = ga_subst_back(r, &st);
        expr_free(r);
        out = eval_and_free(back);
    }

    fmpz_mpoly_q_clear(q, mctx);
    fmpz_mpoly_ctx_clear(mctx);
    varset_free(&fvars);
    expr_free(ebar);
    ga_free(&st);
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
Expr* flint_algebraic_field_normalize(const Expr* e) { (void)e; return NULL; }
Expr* flint_algebraic_field_canonical(const Expr* e) { (void)e; return NULL; }
Expr* flint_algebraic_field_together(const Expr* e) { (void)e; return NULL; }
Expr* flint_rational_together(const Expr* e) { (void)e; return NULL; }
int   flint_rde_base_solve_fg(const Expr* f, const Expr* g,
                              const char* xvar, Expr** y_out) {
    (void)f; (void)g; (void)xvar;
    if (y_out) *y_out = NULL;
    return -1;
}
void  flint_bridge_init(void) { /* no FLINT: nothing to register */ }

#endif /* USE_FLINT */
