/*
 * flint_num_bridge.c
 * ------------------
 * FLINT arb/acb-backed numeric evaluation of transcendental special functions.
 * See flint_num_bridge.h. The acb result (a complex ball) is rendered to an
 * Expr through numeric_mpfr_make_complex (the same "Im rounds to 0 -> real
 * leaf" convention the hand-rolled kernels use). Everything is behind
 * USE_FLINT && USE_MPFR (acb -> Expr needs MPFR); the fallback build provides
 * stubs so the module links either way.
 */

#include "flint_num_bridge.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#if defined(USE_FLINT) && defined(USE_MPFR)

#include "numeric_complex.h"   /* numeric_mpfr_make_complex */
#include "numeric.h"           /* numeric_min_inexact_bits  */
#include <gmp.h>
#include <mpfr.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_dirichlet.h>

/* ------------------------------------------------------------------ */
/*  Expr -> acb  (numeric scalars only)                                */
/* ------------------------------------------------------------------ */

/* A real numeric scalar (Integer / BigInt / Rational / Real / MPFR) -> arb.
 * Set exactly (radius 0) so the ball's midpoint carries the input value; the
 * kernel's own error bounds then govern the result. Returns 0 for anything
 * else (symbolic / Complex / compound), letting the caller bail. */
static int scalar_to_arb(const Expr* e, arb_t out, slong prec) {
    switch (e->type) {
        case EXPR_INTEGER:
            arb_set_si(out, e->data.integer);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            arb_set_fmpz(out, z);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_REAL:
            arb_set_d(out, e->data.real);
            return 1;
        case EXPR_MPFR:
            arf_set_mpfr(arb_midref(out), e->data.mpfr);
            mag_zero(arb_radref(out));
            return 1;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h && h->type == EXPR_SYMBOL
                && strcmp(h->data.symbol.name, "Rational") == 0
                && e->data.function.arg_count == 2) {
                const Expr* a = e->data.function.args[0];
                const Expr* b = e->data.function.args[1];
                if ((a->type != EXPR_INTEGER && a->type != EXPR_BIGINT) ||
                    (b->type != EXPR_INTEGER && b->type != EXPR_BIGINT)) return 0;
                fmpq_t q; fmpq_init(q);
                if (a->type == EXPR_INTEGER) fmpz_set_si(fmpq_numref(q), a->data.integer);
                else                          fmpz_set_mpz(fmpq_numref(q), a->data.bigint);
                if (b->type == EXPR_INTEGER) fmpz_set_si(fmpq_denref(q), b->data.integer);
                else                          fmpz_set_mpz(fmpq_denref(q), b->data.bigint);
                fmpq_canonicalise(q);
                arb_set_fmpq(out, q, prec);
                fmpq_clear(q);
                return 1;
            }
            return 0;
        }
        default:
            return 0;
    }
}

/* A numeric scalar, possibly Complex[re, im], -> acb. */
static int expr_to_acb(const Expr* e, acb_t out, slong prec) {
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h && h->type == EXPR_SYMBOL
            && strcmp(h->data.symbol.name, "Complex") == 0
            && e->data.function.arg_count == 2) {
            return scalar_to_arb(e->data.function.args[0], acb_realref(out), prec)
                && scalar_to_arb(e->data.function.args[1], acb_imagref(out), prec);
        }
    }
    arb_zero(acb_imagref(out));
    return scalar_to_arb(e, acb_realref(out), prec);
}

/* Extract an fmpz from an integer-like Expr. */
static int expr_to_fmpz(const Expr* e, fmpz_t out) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  acb -> Expr                                                        */
/* ------------------------------------------------------------------ */

/* Render a finite acb ball to a numeric Expr at `out_bits` precision. NULL if
 * the result is not finite (pole / overflow / indeterminate). */
static Expr* acb_to_expr(const acb_t z, slong out_bits) {
    if (!acb_is_finite(z)) return NULL;
    mpfr_prec_t mp = out_bits > 2 ? (mpfr_prec_t)out_bits : 53;
    mpfr_t re, im;
    mpfr_init2(re, mp);
    mpfr_init2(im, mp);
    arf_get_mpfr(re, arb_midref(acb_realref(z)), MPFR_RNDN);
    arf_get_mpfr(im, arb_midref(acb_imagref(z)), MPFR_RNDN);
    Expr* out = numeric_mpfr_make_complex(re, im);
    mpfr_clear(re);
    mpfr_clear(im);
    return out;
}

/* Target precision from the numeric arguments: the min inexact-bit count over
 * the arguments (an exact arg imposes no constraint), floored at 53. */
static slong pick_out_bits(const Expr* const* args, int n) {
    long best = 0;
    for (int i = 0; i < n; i++) {
        long b = numeric_min_inexact_bits(args[i]);
        if (b > 0 && (best == 0 || b < best)) best = b;
    }
    if (best < 53) best = 53;
    return (slong)best;
}

/* Guard bits added to the working precision above the target. */
#define NB_GUARD 32

/* ------------------------------------------------------------------ */
/*  Kernels                                                            */
/* ------------------------------------------------------------------ */

Expr* flint_num_zeta(const Expr* s) {
    const Expr* args[1] = { s };
    slong outb = pick_out_bits(args, 1);
    slong wp = outb + NB_GUARD;
    acb_t S, R; acb_init(S); acb_init(R);
    Expr* out = NULL;
    if (expr_to_acb(s, S, wp)) {
        acb_dirichlet_zeta(R, S, wp);
        out = acb_to_expr(R, outb);
    }
    acb_clear(S); acb_clear(R);
    return out;
}

Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a) {
    const Expr* args[2] = { s, a };
    slong outb = pick_out_bits(args, 2);
    slong wp = outb + NB_GUARD;
    acb_t S, A, R; acb_init(S); acb_init(A); acb_init(R);
    Expr* out = NULL;
    if (expr_to_acb(s, S, wp) && expr_to_acb(a, A, wp)) {
        acb_dirichlet_hurwitz(R, S, A, wp);
        out = acb_to_expr(R, outb);
    }
    acb_clear(S); acb_clear(A); acb_clear(R);
    return out;
}

Expr* flint_num_polygamma(const Expr* n, const Expr* z) {
    const Expr* args[2] = { n, z };
    slong outb = pick_out_bits(args, 2);
    slong wp = outb + NB_GUARD;
    acb_t N, Z, R; acb_init(N); acb_init(Z); acb_init(R);
    Expr* out = NULL;
    if (expr_to_acb(n, N, wp) && expr_to_acb(z, Z, wp)) {
        acb_polygamma(R, N, Z, wp);
        out = acb_to_expr(R, outb);
    }
    acb_clear(N); acb_clear(Z); acb_clear(R);
    return out;
}

/* StieltjesGamma[n] or StieltjesGamma[n, a]; n a non-negative integer. */
Expr* flint_num_stieltjes(const Expr* n, const Expr* a) {
    fmpz_t N; fmpz_init(N);
    if (!expr_to_fmpz(n, N) || fmpz_sgn(N) < 0) { fmpz_clear(N); return NULL; }
    const Expr* args[1] = { a ? a : n };  /* a governs precision; default a=1 exact */
    slong outb = a ? pick_out_bits(args, 1) : 53;
    slong wp = outb + NB_GUARD;
    acb_t A, R; acb_init(A); acb_init(R);
    Expr* out = NULL;
    int ok = 1;
    if (a) ok = expr_to_acb(a, A, wp);
    else   acb_one(A);
    if (ok) {
        acb_dirichlet_stieltjes(R, N, A, wp);
        out = acb_to_expr(R, outb);
    }
    acb_clear(A); acb_clear(R);
    fmpz_clear(N);
    return out;
}

/* ------------------------------------------------------------------ */
/*  FLINT` context builtins                                            */
/* ------------------------------------------------------------------ */

static Expr* builtin_flint_zeta(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return flint_num_zeta(res->data.function.args[0]);
}

static Expr* builtin_flint_hurwitzzeta(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return flint_num_hurwitz_zeta(res->data.function.args[0],
                                  res->data.function.args[1]);
}

static Expr* builtin_flint_polygamma(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return flint_num_polygamma(res->data.function.args[0],
                               res->data.function.args[1]);
}

static Expr* builtin_flint_stieltjes(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 1) return flint_num_stieltjes(res->data.function.args[0], NULL);
    if (n == 2) return flint_num_stieltjes(res->data.function.args[0],
                                           res->data.function.args[1]);
    return NULL;
}

void flint_num_bridge_init(void) {
    symtab_add_builtin(SYM_FLINT_Zeta, builtin_flint_zeta);
    symtab_get_def(SYM_FLINT_Zeta)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Zeta,
        "FLINT`Zeta[s] gives the numeric value of the Riemann zeta function at "
        "the numeric argument s (real or complex), computed to the precision "
        "of s (machine precision for exact s) via FLINT's rigorous acb "
        "arithmetic (acb_dirichlet_zeta). Unevaluated for symbolic s or at the "
        "pole s = 1.");

    symtab_add_builtin(SYM_FLINT_HurwitzZeta, builtin_flint_hurwitzzeta);
    symtab_get_def(SYM_FLINT_HurwitzZeta)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_HurwitzZeta,
        "FLINT`HurwitzZeta[s, a] gives the numeric value of the Hurwitz zeta "
        "function via FLINT (acb_dirichlet_hurwitz), to the precision of the "
        "arguments. Unevaluated for symbolic arguments or at a pole.");

    symtab_add_builtin(SYM_FLINT_PolyGamma, builtin_flint_polygamma);
    symtab_get_def(SYM_FLINT_PolyGamma)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolyGamma,
        "FLINT`PolyGamma[n, z] gives the numeric value of the n-th derivative "
        "of the digamma function (n = 0 is digamma) via FLINT "
        "(acb_polygamma), to the precision of the arguments. Unevaluated for "
        "symbolic arguments or at a pole.");

    symtab_add_builtin(SYM_FLINT_StieltjesGamma, builtin_flint_stieltjes);
    symtab_get_def(SYM_FLINT_StieltjesGamma)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_StieltjesGamma,
        "FLINT`StieltjesGamma[n] and FLINT`StieltjesGamma[n, a] give the "
        "numeric value of the n-th Stieltjes constant (generalized, at a) for "
        "a non-negative integer n, via FLINT (acb_dirichlet_stieltjes). "
        "Unevaluated for negative or non-integer n.");
}

#else /* !(USE_FLINT && USE_MPFR) */

Expr* flint_num_zeta(const Expr* s) { (void)s; return NULL; }
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a) { (void)s; (void)a; return NULL; }
Expr* flint_num_polygamma(const Expr* n, const Expr* z) { (void)n; (void)z; return NULL; }
Expr* flint_num_stieltjes(const Expr* n, const Expr* a) { (void)n; (void)a; return NULL; }
void  flint_num_bridge_init(void) { /* no FLINT/MPFR: nothing to register */ }

#endif
