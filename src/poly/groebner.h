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

#ifdef __cplusplus
extern "C" {
#endif

struct Expr;

/* Monomial orders.  The first `elim_pivot` variables are the "elim"
 * block under GB_ORDER_ELIM; otherwise the field is ignored. */
typedef enum {
    GB_ORDER_LEX = 0,
    GB_ORDER_GREVLEX,
    GB_ORDER_ELIM
} GBOrder;

typedef struct GBPoly {
    int     n_vars;       /* number of variables */
    GBOrder order;        /* current monomial order */
    int     elim_pivot;   /* number of leading elim vars (ELIM order) */

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
                     GBOrder order, int elim_pivot);

/* Convert a GBPoly back to a Plus[Times[c, var^k, ...], ...] Expr*.
 * Coefficients become Integer/BigInt/Rational as appropriate.  The
 * caller owns the returned Expr*. */
struct Expr* gb_to_expr(const GBPoly* p, struct Expr** vars);

#ifdef __cplusplus
}
#endif

#endif /* GROEBNER_H */
