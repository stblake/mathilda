/* qafactor.c — Trager algebraic-factoring helpers (Phase G).
 *
 * Phase G3 lives here: norm of f ∈ Q(α)[x] via Resultant_y(P_α, g(x, y)).
 *
 * Strategy: lift QAExt and QAUPoly to picocas Expr*'s in two free
 * symbols (x_name, y_name), then call `internal_resultant` to drop into
 * the Sylvester-matrix routine in poly.c.  This keeps qafactor.c a
 * thin glue layer — all multivariate arithmetic is delegated. */

#include "qafactor.h"
#include "expr.h"
#include "expand.h"
#include "internal.h"
#include <gmp.h>
#include <stdlib.h>

/* ============================== mpq → Expr ============================== */

/* Convert an mpq_t to a picocas Expr.  Returns Integer when the
 * denominator is 1, otherwise Rational[num, den] (num, den as
 * normalised Integer / BigInt nodes).  Handles arbitrary precision. */
static Expr* mpq_to_expr(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);

    Expr* result;
    if (mpz_cmp_ui(den, 1) == 0) {
        result = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
    } else {
        Expr* n = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
        Expr* d = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
        Expr* args[2] = { n, d };
        result = expr_new_function(expr_new_symbol("Rational"), args, 2);
    }

    mpz_clear(num);
    mpz_clear(den);
    return result;
}

/* ====================== Expr-tree construction helpers ================== */

/* Build var^k (k >= 0).  k == 0 → Integer 1, k == 1 → Symbol[var]. */
static Expr* var_power(const char* var, long k) {
    if (k == 0) return expr_new_integer(1);
    if (k == 1) return expr_new_symbol(var);
    Expr* args[2] = { expr_new_symbol(var), expr_new_integer((int64_t)k) };
    return expr_new_function(expr_new_symbol("Power"), args, 2);
}

/* Multiply `coef` (an Expr*, ownership taken) by var^k, with sensible
 * elision when coef is ±1 / k is 0. */
static Expr* term_coef_times_var_power(Expr* coef, const char* var, long k) {
    if (k == 0) return coef;

    Expr* vp = var_power(var, k);

    if (coef->type == EXPR_INTEGER && coef->data.integer == 1) {
        expr_free(coef);
        return vp;
    }
    if (coef->type == EXPR_INTEGER && coef->data.integer == -1) {
        expr_free(coef);
        Expr* args[2] = { expr_new_integer(-1), vp };
        return expr_new_function(expr_new_symbol("Times"), args, 2);
    }

    Expr* args[2] = { coef, vp };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Build a polynomial expression from an mpq_t-array of coefficients
 * coefs[0..n], where coefs[i] is the coefficient of var^i. */
static Expr* mpq_array_to_poly_expr(const mpq_t* coefs,
                                    size_t n_plus_1,
                                    const char* var) {
    /* First pass: count non-zero terms so we can size the Plus arg list. */
    size_t n_terms = 0;
    for (size_t i = 0; i < n_plus_1; i++) {
        if (mpq_sgn(coefs[i]) != 0) n_terms++;
    }
    if (n_terms == 0) return expr_new_integer(0);

    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n_terms);
    size_t out = 0;
    for (size_t i = 0; i < n_plus_1; i++) {
        if (mpq_sgn(coefs[i]) == 0) continue;
        Expr* coef = mpq_to_expr(coefs[i]);
        terms[out++] = term_coef_times_var_power(coef, var, (long)i);
    }

    if (n_terms == 1) {
        Expr* r = terms[0];
        free(terms);
        return r;
    }

    Expr* result = expr_new_function(expr_new_symbol("Plus"), terms, n_terms);
    free(terms);
    return result;
}

/* ============================== Public API ============================== */

Expr* qaext_to_expr(const QAExt* ext, const char* y_name) {
    /* P_α(y) = sum_{i=0..deg} ext->coef[i] * y^i */
    return mpq_array_to_poly_expr(ext->coef, ext->deg + 1, y_name);
}

Expr* qaupoly_to_expr(const QAUPoly* f,
                      const char* x_name,
                      const char* y_name) {
    if (f->deg < 0) return expr_new_integer(0);

    /* Build f(x, α) as a sum of x-monomials whose coefficients are
     * polynomials-in-y derived from each QANum's α-coefficients. */
    size_t n_x_terms = 0;
    for (int i = 0; i <= f->deg; i++) {
        if (!qa_is_zero(f->c[i])) n_x_terms++;
    }
    if (n_x_terms == 0) return expr_new_integer(0);

    Expr** x_terms = (Expr**)malloc(sizeof(Expr*) * n_x_terms);
    size_t out = 0;
    for (int i = 0; i <= f->deg; i++) {
        if (qa_is_zero(f->c[i])) continue;
        /* Coefficient is a polynomial-in-y of α-degree < ext->deg. */
        Expr* y_poly = mpq_array_to_poly_expr(f->c[i]->coef,
                                              f->ext->deg,
                                              y_name);
        x_terms[out++] = term_coef_times_var_power(y_poly, x_name, (long)i);
    }

    Expr* g;
    if (n_x_terms == 1) {
        g = x_terms[0];
    } else {
        g = expr_new_function(expr_new_symbol("Plus"), x_terms, n_x_terms);
    }
    free(x_terms);
    return g;
}

Expr* qaupoly_norm(const QAUPoly* f,
                   const char* x_name,
                   const char* y_name) {
    if (f->deg < 0) return NULL;

    /* Build P_α(y) and g(x, y). */
    Expr* p_alpha = qaext_to_expr(f->ext, y_name);
    Expr* g       = qaupoly_to_expr(f, x_name, y_name);
    Expr* y_var   = expr_new_symbol(y_name);

    /* Resultant_y(P_α, g) — internal_resultant takes ownership of args
     * (they get attached to a temporary Resultant[...] node which is
     * freed inside internal_call_impl regardless of outcome). */
    Expr* args[3] = { p_alpha, g, y_var };
    Expr* r = internal_resultant(args, 3);
    if (!r) return NULL;

    /* The resultant routine already calls Expand internally on the
     * Sylvester determinant, but the fast paths (Times / Power
     * factorisations of P or Q) skip that — re-expand to guarantee a
     * canonical sum-of-monomials form. */
    Expr* expanded = expr_expand(r);
    expr_free(r);
    return expanded;
}
