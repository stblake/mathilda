/*
 * mpoly.h
 * -------
 * Sparse multivariate polynomial over the integers Z with arbitrary-
 * precision (GMP `mpz_t`) coefficients.
 *
 * Used as the substrate for the n-variate Hensel factoring pipeline
 * (Phase F2 of plans/FACTOR_PLAN.md).  Where ZUPoly is "univariate Z[x]"
 * and BPoly is "bivariate Z[x, y]", MPoly is "Z[x_1, ..., x_n]" for
 * arbitrary n.
 *
 * Storage:
 *   - Sparse: each non-zero monomial is stored as (exponent_tuple,
 *     coefficient).  Zero-coefficient terms are never stored.
 *   - Two parallel row-major arrays: `exps` of `n_terms * n_vars`
 *     ints and `coefs` of `n_terms` mpz_t.
 *   - Terms are kept in lex descending order: term i < term j iff
 *     exps[i] > exps[j] under lexicographic comparison (highest
 *     x_0 first, then x_1, etc.).  Term 0 is the leading term.
 *
 * Invariants (after every public operation):
 *   - n_vars >= 0.  (n_vars == 0 represents Z, with n_terms <= 1.)
 *   - n_terms == 0 iff the polynomial is zero.
 *   - All exponents are >= 0.
 *   - All stored coefficients are non-zero.
 *   - Terms are in lex descending order with no duplicate monomials.
 *
 * Memory:
 *   - All MPoly* are heap-allocated; use `mpoly_free` to destroy.
 *   - Functions returning `MPoly*` always return a fresh allocation.
 *   - `mpoly_get_coef` returns a borrowed pointer to an internal
 *     mpz_t (or NULL if the monomial is absent); do NOT free it.
 *
 * Thread safety: the `mpoly_normalize` and `mpoly_mul` paths use a
 * file-static sort context; not safe for concurrent use across
 * threads.  Mathilda is single-threaded so this is acceptable.
 */

#ifndef MPOLY_H
#define MPOLY_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MPoly {
    int     n_vars;
    int*    exps;       /* n_terms * n_vars ints, row-major */
    mpz_t*  coefs;      /* n_terms mpz_t */
    size_t  n_terms;
    size_t  cap;        /* allocated capacity (>= n_terms) */
} MPoly;

/* ---------------------------------------------------------------- */
/*  Construction / destruction                                      */
/* ---------------------------------------------------------------- */

/* Empty (zero) polynomial in n_vars variables. */
MPoly* mpoly_new(int n_vars);

/* Alias for mpoly_new -- legibility at call sites. */
MPoly* mpoly_zero(int n_vars);

/* Constant integer polynomial. */
MPoly* mpoly_from_int(int n_vars, int64_t c);
MPoly* mpoly_from_mpz(int n_vars, const mpz_t c);

/* Single-variable monomial: x_i^k * c.  Useful for constructing
 * test inputs and seed values for the lift. */
MPoly* mpoly_monomial(int n_vars, int var_idx, int k, int64_t c);

/* Deep copy. */
MPoly* mpoly_copy(const MPoly* p);

void   mpoly_free(MPoly* p);

/* Reserve capacity for at least `cap` terms.  Useful before bulk
 * push operations to avoid repeated realloc. */
void   mpoly_reserve(MPoly* p, size_t cap);

/* ---------------------------------------------------------------- */
/*  Term manipulation                                               */
/* ---------------------------------------------------------------- */

/* Append a term (exps borrowed, coef copied) to the END of p without
 * sorting or merging.  After bulk push the caller MUST call
 * `mpoly_normalize` before any other public operation. */
void   mpoly_push_term(MPoly* p, const int* exps, const mpz_t coef);
void   mpoly_push_term_si(MPoly* p, const int* exps, int64_t coef);

/* Set the coefficient of the given monomial.  If the monomial exists,
 * its coefficient is replaced; if not, a new term is inserted in
 * sorted order.  Setting coef = 0 removes the monomial. */
void   mpoly_set_coef(MPoly* p, const int* exps, const mpz_t coef);

/* Look up the coefficient of the given monomial.  Returns a borrowed
 * pointer to the internal mpz_t, or NULL if the monomial is absent. */
const mpz_t* mpoly_get_coef(const MPoly* p, const int* exps);

/* Sort terms into lex descending order, merge duplicates, drop
 * zero-coefficient entries.  Idempotent on already-normal polynomials. */
void   mpoly_normalize(MPoly* p);

/* ---------------------------------------------------------------- */
/*  Predicates and queries                                          */
/* ---------------------------------------------------------------- */

bool   mpoly_is_zero(const MPoly* p);
bool   mpoly_eq(const MPoly* a, const MPoly* b);

/* Maximum power of variable var_idx over all terms.  Returns -1 for
 * the zero polynomial. */
int    mpoly_deg_var(const MPoly* p, int var_idx);

/* Total degree (max sum of exponents over all terms).  -1 for zero. */
int    mpoly_total_deg(const MPoly* p);

/* True iff p does not depend on x_var_idx (every term has exponent 0
 * for that variable).  Equivalent to deg_var(p, var_idx) <= 0.  The
 * zero polynomial counts as constant in every variable. */
bool   mpoly_is_constant_in_var(const MPoly* p, int var_idx);

/* ---------------------------------------------------------------- */
/*  Arithmetic                                                      */
/* ---------------------------------------------------------------- */

MPoly* mpoly_add(const MPoly* a, const MPoly* b);
MPoly* mpoly_sub(const MPoly* a, const MPoly* b);
MPoly* mpoly_neg(const MPoly* a);
MPoly* mpoly_mul(const MPoly* a, const MPoly* b);
MPoly* mpoly_scale(const MPoly* a, const mpz_t k);
MPoly* mpoly_scale_si(const MPoly* a, int64_t k);

/* ---------------------------------------------------------------- */
/*  Substitution and projection                                     */
/* ---------------------------------------------------------------- */

/* Substitute x_var_idx := alpha (an integer) in p.  Returns a new
 * MPoly with the same n_vars; every output term's exponent for
 * var_idx is 0 (the variable is "consumed" but the slot remains).
 * Useful for test-point evaluation in Hensel lifts. */
MPoly* mpoly_subst_var_int(const MPoly* p, int var_idx, int64_t alpha);
MPoly* mpoly_subst_var_mpz(const MPoly* p, int var_idx, const mpz_t alpha);

/* Substitute x_var_idx := x_var_idx + alpha.  Returns a new MPoly
 * with the same n_vars and same monomial structure (just shifted).
 * Used in Hensel lifts to translate the lift point to y = 0. */
MPoly* mpoly_shift_var_int(const MPoly* p, int var_idx, int64_t alpha);

/* Coefficient extraction: returns the polynomial that is the
 * coefficient of x_var_idx^k in p.  The returned polynomial keeps
 * the same n_vars; every term's exponent for var_idx is forced to 0. */
MPoly* mpoly_coef_of_var(const MPoly* p, int var_idx, int k);

/* Leading coefficient in x_var_idx (i.e. coef_of_var at deg_var). */
MPoly* mpoly_lc_var(const MPoly* p, int var_idx);

/* ---------------------------------------------------------------- */
/*  Expr round-trip                                                 */
/* ---------------------------------------------------------------- */

struct Expr;

/* Convert an Expr to MPoly with caller-supplied variable ordering.
 * `vars` is an array of n_vars Expr* (typically symbols); the MPoly
 * variable index i corresponds to vars[i].  Returns NULL if e is not
 * a polynomial in those variables (contains non-integer coefficients,
 * other symbols, or non-polynomial subexpressions like Sin[x]). */
MPoly* expr_to_mpoly(struct Expr* e, struct Expr** vars, int n_vars);

/* Convert MPoly back to Expr.  Returns a Plus[Times[c, vars...], ...]
 * tree (or 0 for the zero polynomial, or a single Times for a single-
 * term polynomial).  Caller owns the returned Expr*. */
struct Expr* mpoly_to_expr(const MPoly* p, struct Expr** vars);

#ifdef __cplusplus
}
#endif

#endif /* MPOLY_H */
