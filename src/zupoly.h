/*
 * zupoly.h
 * --------
 * Univariate polynomials over the integers Z with arbitrary-precision
 * coefficients (GMP `mpz_t`).
 *
 * Used as the coefficient type for the bivariate factoring pipeline
 * (`bpoly.h`, `mvfactor.c`): a bivariate polynomial in Z[x, y] is
 * stored as a ZUPoly[] indexed by x-degree, where each entry is a
 * polynomial in y.
 *
 * Why a separate type from facpoly.c's `UPoly`?
 *   - `UPoly` uses 64-bit machine integers and is purpose-built for
 *     modular arithmetic in finite-field algorithms (Cantor-Zassenhaus
 *     EDF/DDF, Hensel lifting modulo p^k where p^k stays well below
 *     2^63).  Its coefficient growth is bounded by the modulus.
 *   - The bivariate Hensel lift operates over Z[y] with no modular
 *     bound, and coefficients can grow exponentially with the lift
 *     count (Mignotte bound ~2^deg * |P|_∞).  Silent 64-bit overflow
 *     here would produce wrong factorisations.  GMP eliminates that
 *     class of bugs.
 *
 * Memory:
 *   - All ZUPoly* are heap-allocated.  Use `zupoly_new` / `zupoly_zero`
 *     / `zupoly_from_int` to construct, `zupoly_free` to destroy.
 *   - Functions returning `ZUPoly*` always return a fresh allocation
 *     (never a borrowed reference).  Caller owns it.
 *   - `zupoly_getcoef` returns a borrowed pointer to the internal
 *     mpz_t and must NOT be freed by the caller.
 *
 * Invariants (held after every public operation):
 *   - `deg >= -1`, with `deg == -1` representing the zero polynomial.
 *   - For `deg >= 0`, `c[deg]` is non-zero (no leading zeros).
 *   - `cap > deg` (room exists for the leading term).
 *   - All slots `c[0..cap-1]` are mpz_init'd; slots `c[deg+1..cap-1]`
 *     hold zero.  This avoids re-init/clear on every coefficient
 *     update at the cost of a small fixed allocation.
 */

#ifndef ZUPOLY_H
#define ZUPOLY_H

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int    deg;    /* highest power with non-zero coefficient; -1 for zero */
    int    cap;    /* allocated coefficient slots; cap > deg */
    mpz_t* c;      /* c[i] is coefficient of x^i, for 0 <= i <= deg.
                    * Slots c[deg+1..cap-1] are zeroed and reusable. */
} ZUPoly;

/* ---------------------------------------------------------------------- */
/*  Construction / destruction                                            */
/* ---------------------------------------------------------------------- */

/* Allocate a zero polynomial with at least `cap` coefficient slots
 * pre-initialised.  cap >= 1; pass 1 for the zero polynomial. */
ZUPoly* zupoly_new(int cap);

/* The zero polynomial. */
ZUPoly* zupoly_zero(void);

/* The constant polynomial p(x) = n. */
ZUPoly* zupoly_from_int(int64_t n);

/* Deep copy. */
ZUPoly* zupoly_copy(const ZUPoly* p);

/* Free.  Safe to pass NULL. */
void zupoly_free(ZUPoly* p);

/* ---------------------------------------------------------------------- */
/*  Coefficient access                                                    */
/* ---------------------------------------------------------------------- */

/* Set the coefficient of x^i to `v`.  Resizes as needed.  After this
 * call, `p->deg >= i` if v != 0 (and may grow), or `p->deg < i` if v
 * was zeroing the leading coefficient. */
void zupoly_setcoef(ZUPoly* p, int i, const mpz_t v);
void zupoly_setcoef_si(ZUPoly* p, int i, int64_t v);

/* Get coefficient of x^i.  For i in [0, deg], returns a borrowed
 * pointer to the internal mpz_t.  For i out of range or for the zero
 * polynomial, returns NULL. */
const mpz_t* zupoly_getcoef(const ZUPoly* p, int i);

/* Re-establish the deg/cap invariants after direct mutation of c[]. */
void zupoly_normalize(ZUPoly* p);

/* ---------------------------------------------------------------------- */
/*  Predicates                                                            */
/* ---------------------------------------------------------------------- */

bool zupoly_is_zero(const ZUPoly* p);
bool zupoly_eq(const ZUPoly* a, const ZUPoly* b);

/* Compare lex by degree first, then by coefficients high-to-low.  The
 * resulting order is well-defined (total) but not arithmetically
 * meaningful -- use only for canonical sorting / equality. */
int zupoly_cmp(const ZUPoly* a, const ZUPoly* b);

/* ---------------------------------------------------------------------- */
/*  Arithmetic                                                            */
/* ---------------------------------------------------------------------- */

ZUPoly* zupoly_add(const ZUPoly* a, const ZUPoly* b);
ZUPoly* zupoly_sub(const ZUPoly* a, const ZUPoly* b);
ZUPoly* zupoly_mul(const ZUPoly* a, const ZUPoly* b);
ZUPoly* zupoly_neg(const ZUPoly* a);

/* a * s where s is an mpz_t scalar. */
ZUPoly* zupoly_scale(const ZUPoly* a, const mpz_t s);
ZUPoly* zupoly_scale_si(const ZUPoly* a, int64_t s);

/* Divide a by b returning (q, r) with a = q*b + r and deg(r) < deg(b).
 * Returns true on success, false if b is zero.  Works correctly only
 * when b is monic OR when the division is exact; for general b over Z
 * use `zupoly_pseudodivrem`.  The caller is responsible for freeing
 * *q and *r. */
bool zupoly_divrem_monic(const ZUPoly* a, const ZUPoly* b,
                         ZUPoly** q, ZUPoly** r);

/* Pseudo-division: lc(b)^(deg(a)-deg(b)+1) * a = q*b + r with
 * deg(r) < deg(b).  Always succeeds for non-zero b over Z (no
 * rational coefficients introduced).  Returns true on success. */
bool zupoly_pseudodivrem(const ZUPoly* a, const ZUPoly* b,
                         ZUPoly** q, ZUPoly** r);

/* Exact division.  Returns a/b on success, NULL if b ∤ a (over Z). */
ZUPoly* zupoly_divexact(const ZUPoly* a, const ZUPoly* b);

/* GCD over Z, computed via subresultant pseudo-remainder sequence.
 * Result has positive leading coefficient and is primitive (content
 * 1).  gcd(0, 0) = 0; gcd(0, b) = primitive_part(b) up to sign. */
ZUPoly* zupoly_gcd(const ZUPoly* a, const ZUPoly* b);

/* Content: positive gcd of all coefficients (an mpz_t).  Sets `out`. */
void zupoly_content(const ZUPoly* p, mpz_t out);

/* Primitive part: p / content(p).  Sign matches sign of lc(p). */
ZUPoly* zupoly_primitive_part(const ZUPoly* p);

/* ---------------------------------------------------------------------- */
/*  Substitution / evaluation                                             */
/* ---------------------------------------------------------------------- */

/* Evaluate p(α) ∈ Z.  Sets `out`. */
void zupoly_eval_si(const ZUPoly* p, int64_t alpha, mpz_t out);
void zupoly_eval(const ZUPoly* p, const mpz_t alpha, mpz_t out);

/* Shift the variable: returns p(x + α). */
ZUPoly* zupoly_shift_si(const ZUPoly* p, int64_t alpha);

/* ---------------------------------------------------------------------- */
/*  Conversion to / from picocas Expr                                     */
/* ---------------------------------------------------------------------- */

struct Expr;

/* Try to interpret `e` as a polynomial in `var` with integer
 * coefficients.  Returns NULL if any coefficient is not an integer
 * (e.g., rational, real, or symbolic non-constant).  Caller owns
 * the returned ZUPoly.
 *
 * `e` is expected to be expanded; if not, expansion happens internally. */
ZUPoly* expr_to_zupoly(const struct Expr* e, const struct Expr* var);

/* Build a fresh Expr representing p(var).  Caller owns the result. */
struct Expr* zupoly_to_expr(const ZUPoly* p, const struct Expr* var);

/* ---------------------------------------------------------------------- */
/*  Debug                                                                 */
/* ---------------------------------------------------------------------- */

/* Print as `c0 + c1*var + c2*var^2 + ...` to stderr. */
void zupoly_print(const ZUPoly* p, const char* var);

/* ---------------------------------------------------------------------- */
/*  Diophantine equation (Bezout-style) over Z[x]                         */
/* ---------------------------------------------------------------------- */

/* Solve  delta_u(x) * v(x) + delta_v(x) * u(x) = e(x)  for delta_u,
 * delta_v ∈ Z[x] with deg(delta_u) < deg(u).
 *
 * Preconditions:
 *   - u and v are monic (lc(u) == lc(v) == 1).
 *   - u and v are coprime over Q[x] (equivalently, gcd(u, v) is a
 *     constant in Q).
 *   - e is integer-coefficient.
 *
 * The implementation runs the extended Euclidean algorithm in Q[x]
 * via internal mpq_t arithmetic.  For monic coprime u, v with
 * integer e, the Diophantine solution is unique in Q[x] subject to
 * the degree constraint, and -- for the inputs that arise in
 * bivariate Hensel lifting after Wang's leading-coefficient
 * correction -- it lies in Z[x].
 *
 * Returns true on success, with *du_out and *dv_out owned by the
 * caller.  Returns false if:
 *   - u or v is not monic / not coprime.
 *   - The solution exists in Q[x] but has non-integer coefficients
 *     (the caller should consider a different lifting strategy).
 *
 * On failure, *du_out and *dv_out are set to NULL. */
bool zupoly_diophantine(const ZUPoly* u, const ZUPoly* v, const ZUPoly* e,
                        ZUPoly** du_out, ZUPoly** dv_out);

#ifdef __cplusplus
}
#endif

#endif /* ZUPOLY_H */
