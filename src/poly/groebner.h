/* groebner.h
 *
 * Buchberger Gröbner-basis core for the GroebnerBasis[] builtin.
 *
 * The polynomial substrate is a sparse multivariate polynomial over Q
 * (rational coefficients, GMP `mpq_t`) carried as parallel arrays of
 * exponent tuples and rational coefficients.  This is independent of
 * `src/poly/mpoly.{c,h}` (which is restricted to Z and lex order); we
 * keep it self-contained because Buchberger needs:
 *
 *   - rational coefficients (the default Mathematica coefficient domain),
 *   - a configurable monomial order (lex, grevlex, elimination), and
 *   - re-sorting after every arithmetic op since the leading-monomial
 *     query is order-dependent.
 *
 * Memory: every `GBPoly*` returned by this module is heap-allocated and
 * must be released by `gb_poly_free`.  Internal arrays are realloc-grown
 * on demand; the caller never touches them directly.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

struct Expr;

/* Monomial orders.  The first `elim_pivot` variables are the "elim"
 * block under GB_ORDER_ELIM; otherwise the field is ignored.
 *
 * GB_ORDER_MATRIX compares monomials by a user-supplied integer weight
 * matrix (see GBWeightMatrix): exponent vectors a, b are ranked by the
 * lexicographic comparison of the weight vectors M*a vs M*b (larger
 * weight = larger monomial).  The standard orders are special cases.
 * A GBPoly using this order borrows a `wmat` pointer that the caller
 * keeps alive for the whole computation. */
typedef enum {
    GB_ORDER_LEX = 0,
    GB_ORDER_GREVLEX,
    GB_ORDER_ELIM,
    GB_ORDER_MATRIX
} GBOrder;

/* A weight matrix defining a monomial order on `n_vars` variables.
 * `w` is an `n_rows * n_vars` row-major integer array, caller-owned and
 * outliving every GBPoly that references it.  Construct directly; there
 * is no allocator here (the field arrays are owned by whoever built the
 * struct). */
typedef struct GBWeightMatrix {
    int      n_rows;
    int      n_vars;
    int64_t* w;           /* n_rows * n_vars, row-major; BORROWED by polys */
} GBWeightMatrix;

typedef struct GBPoly {
    int     n_vars;       /* number of variables */
    GBOrder order;        /* current monomial order */
    int     elim_pivot;   /* number of leading elim vars (ELIM order) */
    const GBWeightMatrix* wmat;  /* borrowed; non-NULL iff order==MATRIX */

    int*    exps;         /* n_terms * n_vars, row-major */
    mpq_t*  coefs;        /* n_terms */
    size_t  n_terms;
    size_t  cap;
} GBPoly;

/* ------------------------------------------------------------------ */
/*  Construction / destruction                                         */
/* ------------------------------------------------------------------ */

GBPoly* gb_poly_new(int n_vars, GBOrder order, int elim_pivot);
GBPoly* gb_poly_copy(const GBPoly* p);
void    gb_poly_free(GBPoly* p);

/* Reserve capacity for `cap` terms. */
void    gb_poly_reserve(GBPoly* p, size_t cap);

/* Switch `p` to GB_ORDER_MATRIX with the borrowed weight matrix `wmat`.
 * Does NOT re-sort: the caller must follow up with gb_poly_normalize()
 * (the term order has changed, so the existing term order is stale). */
void    gb_poly_set_wmat(GBPoly* p, const GBWeightMatrix* wmat);

/* True iff the `n_rows * n_vars` integer matrix `w` defines a global
 * monomial (term) order on `n_vars` variables.  Both conditions are
 * required and independent:
 *   (1) rank(w) == n_vars  -- the order is injective on exponent vectors
 *       (otherwise gb_poly_normalize would merge distinct monomials);
 *   (2) for every column, its first non-zero entry (top to bottom) is
 *       positive -- each variable is > 1, so the order is well-founded
 *       (otherwise Buchberger need not terminate). */
bool    gb_wmat_validate(const int64_t* w, int n_rows, int n_vars);

/* Append a term to the END without sorting; caller MUST follow up with
 * gb_poly_normalize() before any other operation. */
void    gb_poly_push_term(GBPoly* p, const int* exps, const mpq_t coef);
void    gb_poly_push_term_si(GBPoly* p, const int* exps, int64_t num, int64_t den);

/* Sort by the active order, merge same-monomial duplicates, drop zero
 * coefficients.  Idempotent. */
void    gb_poly_normalize(GBPoly* p);

/* ------------------------------------------------------------------ */
/*  Predicates and queries                                             */
/* ------------------------------------------------------------------ */

bool         gb_poly_is_zero(const GBPoly* p);
bool         gb_poly_is_constant(const GBPoly* p);
const int*   gb_poly_lm(const GBPoly* p);     /* borrowed; NULL if zero */
const mpq_t* gb_poly_lc(const GBPoly* p);     /* borrowed; NULL if zero */

/* True iff every term of `p` has exponent 0 for every variable in
 * `vars` (an array of `n` variable indices).  Used by the elimination
 * filter at the end of the 3-arg form. */
bool gb_poly_free_of_vars(const GBPoly* p, const int* vars, int n);

/* ------------------------------------------------------------------ */
/*  Arithmetic                                                         */
/* ------------------------------------------------------------------ */

GBPoly* gb_poly_neg(const GBPoly* a);
GBPoly* gb_poly_add(const GBPoly* a, const GBPoly* b);
GBPoly* gb_poly_sub(const GBPoly* a, const GBPoly* b);
GBPoly* gb_poly_scale(const GBPoly* a, const mpq_t c);
GBPoly* gb_poly_mul_by_monomial(const GBPoly* a, const int* exps, const mpq_t c);

/* Scale by the multiplicative inverse of the leading coefficient,
 * yielding a monic polynomial in-place. */
void gb_poly_make_monic(GBPoly* p);

/* ------------------------------------------------------------------ */
/*  Buchberger primitives                                              */
/* ------------------------------------------------------------------ */

/* S-polynomial of two non-zero polynomials.  Returns a fresh GBPoly
 * normalised under the shared order. */
GBPoly* gb_spoly(const GBPoly* f, const GBPoly* g);

/* Multivariate division of `p` by the list `basis[0..n-1]`.  Repeatedly
 * subtracts off `c * x^e * basis[i]` for each term of `p` whose
 * monomial is divisible by some basis leading monomial.  Returns the
 * remainder (a fresh GBPoly; the caller owns it).  The basis polys are
 * not modified. */
GBPoly* gb_reduce(const GBPoly* p, GBPoly* const* basis, size_t n);

/* Like gb_reduce, but also records the division quotients: on return
 *   p == sum_i (*quot_out)[i] * basis[i] + remainder
 * with the remainder being the returned polynomial.  `*quot_out` is a
 * freshly-allocated array of `n` GBPoly* (one quotient per basis
 * element); the caller owns each GBPoly* AND the array, and must release
 * them with gb_poly_free / free.  Used by the Gröbner walk lift step. */
GBPoly* gb_divmod(const GBPoly* p, GBPoly* const* basis, size_t n,
                  GBPoly*** quot_out);

/* The w-initial form of `g`: the sub-polynomial of all terms whose
 * weight  sum_v w[v]*exp[v]  is maximal.  Returns a fresh GBPoly sharing
 * `g`'s order/wmat (caller owns it).  `g` must be non-zero. */
GBPoly* gb_initial_form(const GBPoly* g, const int64_t* w, int n_vars);

/* Run Buchberger's algorithm with Gebauer–Möller pair criteria.
 * Inputs: `F[0..n-1]` -- read-only; the function takes ownership of
 * neither (caller frees).
 * Output: a fresh dynamically-allocated array of `*out_n` GBPoly*
 * forming a Gröbner basis under the common order.  The basis is
 * interreduced and made monic.  Caller owns each GBPoly* in the array
 * AND the array itself; release with `gb_basis_free`. */
GBPoly** gb_buchberger(GBPoly* const* F, size_t n, size_t* out_n);

/* Free a basis array as returned by gb_buchberger. */
void gb_basis_free(GBPoly** basis, size_t n);

/* Turn a generating set that is already a Gröbner basis under G[0]'s
 * order into the canonical reduced basis: interreduce, reduce to a fixed
 * point, make primitive over Z, and sort ascending by leading monomial.
 * Operates in place (the array is only compacted); *nG is updated.  Used
 * by the Gröbner walk to canonicalise its final basis. */
void gb_finalize_basis(GBPoly** G, size_t* nG);

/* ------------------------------------------------------------------ */
/*  Gröbner walk                                                        */
/* ------------------------------------------------------------------ */

/* Compute a Gröbner basis of <F[0..n-1]> under the target order via the
 * Collart–Kalkbrener–Mall Gröbner walk: a cheap DegreeReverseLexico-
 * graphic basis is computed first, then converted to the target order by
 * walking the Gröbner fan.  The result is identical (as a reduced,
 * primitive-over-Z, ascending-LM-sorted basis) to gb_buchberger run
 * directly under the target order; the walk is typically faster for
 * orders (such as Lexicographic) that are expensive to compute directly.
 *
 * `target_order` selects the destination order; `target_wmat` is the
 * weight matrix when `target_order == GB_ORDER_MATRIX` (NULL otherwise).
 * Output ownership matches gb_buchberger (caller frees via
 * gb_basis_free).  On any internal safety-guard trip the function falls
 * back to a direct gb_buchberger run, returning the identical basis. */
GBPoly** gb_groebner_walk(GBPoly* const* F, size_t n,
                          GBOrder target_order,
                          const GBWeightMatrix* target_wmat,
                          size_t* out_n);

/* ------------------------------------------------------------------ */
/*  Expr round-trip                                                    */
/* ------------------------------------------------------------------ */

/* Convert an Expr polynomial-in-`vars` to a GBPoly with the given
 * monomial order.  Handles Integer, BigInt, and Rational coefficients
 * and the standard Plus/Times/Power tree shape.  Returns NULL if `e`
 * contains a non-polynomial substructure in the given variables (any
 * symbol that is not in `vars`, any non-integer Power exponent, any
 * head that is not Plus/Times/Power/Rational/integer literal). */
GBPoly* gb_from_expr(struct Expr* e, struct Expr** vars, int n_vars,
                     GBOrder order, int elim_pivot,
                     const GBWeightMatrix* wmat);

/* Convert a GBPoly back to a Plus[Times[c, var^k, ...], ...] Expr*.
 * Coefficients become Integer/BigInt/Rational as appropriate.  The
 * caller owns the returned Expr*. */
struct Expr* gb_to_expr(const GBPoly* p, struct Expr** vars);

/* Convert a GBPoly to Expr form for the InexactNumbers coefficient
 * domain: the polynomial is made monic and every coefficient is rendered
 * as a `bits`-bit MPFR real.  The caller owns the returned Expr*. */
struct Expr* gb_to_expr_inexact(const GBPoly* p, struct Expr** vars,
                                mpfr_prec_t bits);

#ifdef __cplusplus
}
#endif

#endif /* GROEBNER_H */
