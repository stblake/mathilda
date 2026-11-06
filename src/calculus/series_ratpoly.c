/*
 * series_ratpoly.c -- truncated power-series arithmetic over Q on FLINT's
 * fmpq_poly, for the series engine.  See series_ratpoly.h for the contract.
 *
 * The conversions are the only per-coefficient work: a coefficient array is
 * loaded by scaling every numerator to the lcm L of the denominators (one
 * fmpz_poly, then one division by L, so the load is linear rather than the
 * quadratic cost of setting rational coefficients one at a time), and stored
 * back one canonical fmpq at a time.  Everything in between -- mullow,
 * inv_series, pow_trunc, compose_series -- is FLINT's quasi-linear arithmetic.
 */
#include "series_ratpoly.h"
#include "sym_names.h"

#include <stdlib.h>

int srp_is_rational(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return 1;
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Rational) return 0;
    const Expr* n = e->data.function.args[0];
    const Expr* d = e->data.function.args[1];
    return (n->type == EXPR_INTEGER || n->type == EXPR_BIGINT) &&
           (d->type == EXPR_INTEGER || d->type == EXPR_BIGINT);
}

#ifdef USE_FLINT

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include "../poly/flint_bridge_internal.h"   /* fb_expr_to_fmpq / fb_fmpq_to_expr */

static int all_rational(Expr* const* c, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!srp_is_rational(c[i])) return 0;
    return 1;
}

/* c[0..n-1] (all rational, checked by the caller) -> out. */
static void load(Expr* const* c, size_t n, fmpq_poly_t out) {
    fmpq* q = flint_malloc(sizeof(fmpq) * (n ? n : 1));
    fmpz_t L, s;
    fmpz_init_set_ui(L, 1);
    fmpz_init(s);
    for (size_t i = 0; i < n; i++) {
        fmpq_init(q + i);
        fb_expr_to_fmpq(c[i], q + i);
        if (!fmpz_is_one(fmpq_denref(q + i))) fmpz_lcm(L, L, fmpq_denref(q + i));
    }
    fmpz_poly_t P;
    fmpz_poly_init2(P, (slong)n);
    for (size_t i = 0; i < n; i++) {
        if (fmpz_is_zero(fmpq_numref(q + i))) continue;
        fmpz_divexact(s, L, fmpq_denref(q + i));
        fmpz_mul(s, s, fmpq_numref(q + i));
        fmpz_poly_set_coeff_fmpz(P, (slong)i, s);
    }
    fmpq_poly_set_fmpz_poly(out, P);
    if (!fmpz_is_one(L)) fmpq_poly_scalar_div_fmpz(out, out, L);
    fmpz_poly_clear(P);
    for (size_t i = 0; i < n; i++) fmpq_clear(q + i);
    flint_free(q);
    fmpz_clear(L);
    fmpz_clear(s);
}

/* p's coefficients of t^0 .. t^(n-1) -> out (fresh owned Exprs). */
static void store(const fmpq_poly_t p, size_t n, Expr** out) {
    fmpq_t c;
    fmpq_init(c);
    slong len = fmpq_poly_length(p);
    for (size_t i = 0; i < n; i++) {
        if ((slong)i >= len) { out[i] = expr_new_integer(0); continue; }
        fmpq_poly_get_coeff_fmpq(c, p, (slong)i);
        out[i] = fb_fmpq_to_expr(c);
    }
    fmpq_clear(c);
}

int srp_mullow(Expr* const* a, size_t na, Expr* const* b, size_t nb,
               size_t n, Expr** out) {
    if (!all_rational(a, na) || !all_rational(b, nb)) return 0;
    fmpq_poly_t A, B, R;
    fmpq_poly_init(A); fmpq_poly_init(B); fmpq_poly_init(R);
    load(a, na, A);
    load(b, nb, B);
    if (n > 0 && !fmpq_poly_is_zero(A) && !fmpq_poly_is_zero(B))
        fmpq_poly_mullow(R, A, B, (slong)n);
    store(R, n, out);
    fmpq_poly_clear(A); fmpq_poly_clear(B); fmpq_poly_clear(R);
    return 1;
}

int srp_inv_series(Expr* const* a, size_t na, size_t n, Expr** out) {
    if (na == 0 || !all_rational(a, na)) return 0;
    fmpq_poly_t A, R;
    fmpq_poly_init(A); fmpq_poly_init(R);
    load(a, na, A);
    int ok = 0;
    fmpq_t a0;
    fmpq_init(a0);
    fmpq_poly_get_coeff_fmpq(a0, A, 0);
    if (!fmpq_is_zero(a0)) {
        if (n > 0) fmpq_poly_inv_series(R, A, (slong)n);
        store(R, n, out);
        ok = 1;
    }
    fmpq_clear(a0);
    fmpq_poly_clear(A); fmpq_poly_clear(R);
    return ok;
}

int srp_pow_trunc(Expr* const* a, size_t na, uint64_t e, size_t n, Expr** out) {
    if (!all_rational(a, na)) return 0;
    fmpq_poly_t A, R;
    fmpq_poly_init(A); fmpq_poly_init(R);
    load(a, na, A);
    if (n > 0) fmpq_poly_pow_trunc(R, A, (ulong)e, (slong)n);
    store(R, n, out);
    fmpq_poly_clear(A); fmpq_poly_clear(R);
    return 1;
}

int srp_compose_series(Expr* const* k, size_t nk, Expr* const* u, size_t nu,
                       size_t n, Expr** out) {
    if (!all_rational(k, nk) || !all_rational(u, nu)) return 0;
    fmpq_poly_t K, U, R;
    fmpq_poly_init(K); fmpq_poly_init(U); fmpq_poly_init(R);
    load(k, nk, K);
    load(u, nu, U);
    int ok = 0;
    fmpq_t u0;
    fmpq_init(u0);
    fmpq_poly_get_coeff_fmpq(u0, U, 0);
    if (fmpq_is_zero(u0)) {
        if (n > 0) fmpq_poly_compose_series(R, K, U, (slong)n);
        store(R, n, out);
        ok = 1;
    }
    fmpq_clear(u0);
    fmpq_poly_clear(K); fmpq_poly_clear(U); fmpq_poly_clear(R);
    return ok;
}

static Expr* binop(const Expr* a, const Expr* b, int mul) {
    if (!srp_is_rational(a) || !srp_is_rational(b)) return NULL;
    fmpq_t x, y;
    fmpq_init(x); fmpq_init(y);
    fb_expr_to_fmpq(a, x);
    fb_expr_to_fmpq(b, y);
    if (mul) fmpq_mul(x, x, y); else fmpq_add(x, x, y);
    Expr* r = fb_fmpq_to_expr(x);
    fmpq_clear(x); fmpq_clear(y);
    return r;
}

Expr* srp_add2(const Expr* a, const Expr* b) { return binop(a, b, 0); }
Expr* srp_mul2(const Expr* a, const Expr* b) { return binop(a, b, 1); }

#else /* !USE_FLINT */

int srp_mullow(Expr* const* a, size_t na, Expr* const* b, size_t nb,
               size_t n, Expr** out) {
    (void)a; (void)na; (void)b; (void)nb; (void)n; (void)out; return 0;
}
int srp_inv_series(Expr* const* a, size_t na, size_t n, Expr** out) {
    (void)a; (void)na; (void)n; (void)out; return 0;
}
int srp_pow_trunc(Expr* const* a, size_t na, uint64_t e, size_t n, Expr** out) {
    (void)a; (void)na; (void)e; (void)n; (void)out; return 0;
}
int srp_compose_series(Expr* const* k, size_t nk, Expr* const* u, size_t nu,
                       size_t n, Expr** out) {
    (void)k; (void)nk; (void)u; (void)nu; (void)n; (void)out; return 0;
}
Expr* srp_add2(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }
Expr* srp_mul2(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }

#endif /* USE_FLINT */
