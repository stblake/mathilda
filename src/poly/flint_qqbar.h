/* flint_qqbar.h — exact algebraic-number canonicalisation for RootReduce.
 *
 * Backs the WL-faithful behaviour of RootReduce on *constant* algebraic
 * numbers (no free symbols): an expression built from integers, rationals,
 * radicals Power[base, p/q], roots of unity (-1)^(p/q), the imaginary unit,
 * and Root[Function[t, poly], k] objects, combined by +, -, *, /, ^, is
 * canonicalised to a single representative — a rational number, a quadratic
 * radical expression, or a Root[Function[minpoly&], k] object.
 *
 * The engine is FLINT's `qqbar` (exact real/complex algebraic numbers via
 * minimal polynomial + isolating enclosure): rigorous, with no numeric zero
 * oracle. When FLINT is not compiled in (USE_FLINT off) every entry point is a
 * graceful no-op returning NULL / -1, matching the rest of the FLINT bridge.
 *
 * Parametric algebraic *functions* (radicals whose radicand carries a free
 * variable, e.g. the Goursat k^(1/3) towers) are NOT handled here — they stay
 * with flint_algebraic_field_canonical in flint_bridge.c. RootReduce dispatches
 * between the two by the presence of free symbols.
 */
#ifndef FLINT_QQBAR_H
#define FLINT_QQBAR_H

#include "expr.h"

/* RootReduce Method selector (see the WL Method -> option). */
typedef enum {
    QQBAR_METHOD_AUTOMATIC = 0,
    QQBAR_METHOD_RECURSIVE,
    QQBAR_METHOD_NUMBERFIELD
} QQBarMethod;

/* True if `e` has no free symbol and every head is one the qqbar converter can
 * handle (a necessary, cheap pre-filter; conversion may still fail on a degree
 * blow-up or an unsupported atom). Integers, rationals, radicals, roots of
 * unity, the imaginary unit, Root objects and the named constant GoldenRatio
 * (the sole algebraic one) qualify; other bare symbols, Pi, E, Log, Sin, ...
 * do not. */
int flint_qqbar_is_constant_algebraic(const Expr* e);

/* Canonicalise the constant algebraic number `e` under the chosen Method.
 * Returns a fresh owned Expr (rational / quadratic radical / Root object) on
 * success, or NULL if `e` is not a constant algebraic number, exceeds the
 * degree cap, or FLINT is unavailable. Never mutates `e`. */
Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod method);

/* Walk `e` and replace every MAXIMAL constant-algebraic subexpression with its
 * flint_qqbar_canonical form (folding combinations of Root / AlgebraicNumber /
 * radicals / roots-of-unity / I / rationals; a vanishing combination collapses
 * to 0). The free-variable structure is left intact — flint_qqbar_canonical
 * declines (NULL) for anything carrying a free symbol, so the walker descends
 * into it. Returns a fresh, *unevaluated* tree (the caller evaluates once, so
 * e.g. Times[0, x] drops); returns a plain copy for an atomic number leaf; NULL
 * only when FLINT is compiled out. Never mutates `e`. Cheap: it reduces numbers,
 * never the multivariate polynomial engine, so it cannot blow up. */
Expr* flint_qqbar_reduce_coeffs(const Expr* e, QQBarMethod method);

/* Exact algebraic equality of two constant algebraic numbers.
 * Returns 1 (equal), 0 (unequal), or -1 (undecided: either operand is not a
 * constant algebraic number, or FLINT is off). */
int flint_qqbar_equal(const Expr* a, const Expr* b);

/* Sign of (a - b) for real constant algebraic numbers: -1, 0, or 1;
 * returns -2 when undecided (non-real, non-algebraic, or FLINT off). */
int flint_qqbar_compare(const Expr* a, const Expr* b);

/* Realness of a constant algebraic number: 1 (a real algebraic number),
 * 0 (a genuinely non-real algebraic number), or -1 (undecided: `e` is not a
 * constant algebraic number, exceeds the degree cap, or FLINT is off).  Lets a
 * Reals-domain filter distinguish "drop this non-real branch" (0) from "cannot
 * decide -- leave the whole statement unevaluated" (-1). */
int flint_qqbar_is_real(const Expr* e);

/* Exact sign of the real part (imag == 0) or the imaginary part (imag != 0)
 * of a constant algebraic number: -1, 0 or 1; -2 when undecided (`e` is not a
 * constant algebraic number, exceeds the degree cap, or FLINT is off).  Backs
 * N's branch-cut side decision (numeric.c, nbc_settle): an exact radicand
 * whose floating-point imaginary part is rounding noise either side of 0. */
int flint_qqbar_part_sign(const Expr* e, int imag);

/* ------------------------------------------------------------------ */
/*  AlgebraicNumber / ToNumberField Expr-level entry points            */
/*  (used by the builtins in algebraicnumber.c / tonumberfield.c).     */
/*  Every one returns a fresh owned Expr on success and NULL to        */
/*  decline — a non-constant-algebraic argument, malformed coeffs, a   */
/*  degree-cap overflow, or FLINT compiled out. None mutate an input.  */
/* ------------------------------------------------------------------ */

/* Canonicalise AlgebraicNumber[gen, coeffs] (coeffs a List of integers /
 * rationals): reduce gen to an algebraic integer and coeffs to the power basis
 * of its minimal polynomial. Returns the reduced AlgebraicNumber[g, {..}], or a
 * plain Integer/Rational when the value is rational. */
Expr* flint_qqbar_algebraic_number(const Expr* gen, const Expr* coeffs);

/* ToNumberField[a, theta]: express the constant algebraic number `a` in the
 * field Q(theta) as AlgebraicNumber[g, {..}] (g the algebraic-integer generator
 * of Q(theta)); returns the rational `a` when Q(theta) = Q. NULL if a is not in
 * Q(theta). */
Expr* flint_qqbar_to_number_field(const Expr* a, const Expr* theta);

/* ToNumberField[{a0..a_{n-1}}] (or `All`): build a common field Q(a0,..) via a
 * primitive element and return List[AlgebraicNumber[g, {..}], ..]. `smallest`
 * is advisory (Automatic and All both use the qqbar primitive element). */
Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n,
                                         int smallest);
/* The same, the numbers read one at a time through get(ctx, i, &owned)
 * (an owned answer is freed after use), or from the array ctx when get is
 * NULL: a LinkedList's slots are read on demand this way; out_ll answers
 * with a LinkedList. */
Expr* flint_qqbar_to_number_field_common_fn(Expr* (*get)(const void*, size_t, bool*),
                                            const void* ctx, size_t n, int smallest, bool out_ll);

/* ToNumberField[x]: express the single algebraic number x as an explicit
 * AlgebraicNumber in Q(x). */
Expr* flint_qqbar_to_number_field_self(const Expr* x);

/* NumberFieldIntegralBasis[a]: a Z-module basis of the ring of integers O_K of
 * K = Q(a), as List[element, ..] with each element an algebraic integer
 * (AlgebraicNumber[g, {..}], or a plain Integer/Rational).  The equation order
 * Z[phi] is enlarged to O_K via Round 2 (Pohst-Zassenhaus) when it is not
 * maximal.  Returns NULL to decline: `a` is not a constant algebraic number,
 * exceeds the degree cap, or O_K cannot be certified (disc will not factor / a
 * prime too large / FLINT compiled out). */
Expr* flint_qqbar_integral_basis(const Expr* a);

/* AlgebraicIntegerQ[x]: 1 (x is an algebraic integer), 0 (it is not), or -1
 * (undecided -- past the qqbar degree cap, or FLINT compiled out).  Ring closure
 * (sums, products, positive powers of integers) decides first; otherwise x is
 * an algebraic integer iff its primitive integer minimal polynomial is monic. */
int flint_qqbar_algebraic_integer_q(const Expr* x);

/* AlgebraicNumberDenominator[x]: the smallest positive integer d such that d*x is
 * an algebraic integer.  Returns 1 and sets *out to that d (a fresh owned Integer
 * or BigInt Expr) on success; 0 when x is not a constant algebraic number; -1 when
 * FLINT is compiled out.  Computed exactly from the minimal polynomial by a
 * per-prime valuation (see flint_qqbar.c) -- not the raw min-poly leading
 * coefficient, which over-counts. */
int flint_qqbar_algebraic_number_denominator(const Expr* x, Expr** out);

/* AlgebraicNumberNorm[a] (theta == NULL): the absolute field norm N_{Q(a)/Q}(a),
 * the product of a's conjugates = (-1)^deg times the monic-over-Q constant term of
 * a's minimal polynomial.  AlgebraicNumberNorm[a, Extension -> theta] (theta != NULL):
 * the relative norm N_{Q(theta)/Q}(a) for a in K = Q(theta), equal to the absolute
 * norm raised to the tower index [K:Q(a)].  Returns 1 with *out set to a fresh owned
 * Integer/Rational; 0 when a (or theta) is not a constant algebraic number; 2 when
 * theta is given but a is not an element of Q(theta); -1 when FLINT is compiled out. */
int flint_qqbar_algebraic_number_norm(const Expr* a, const Expr* theta, Expr** out);

/* AlgebraicNumberTrace[a] (theta == NULL): the absolute field trace Tr_{Q(a)/Q}(a),
 * the sum of a's conjugates = -(coeff of x^{deg-1}) / (leading coeff) of a's minimal
 * polynomial.  AlgebraicNumberTrace[a, Extension -> theta] (theta != NULL): the
 * relative trace Tr_{Q(theta)/Q}(a) for a in K = Q(theta), equal to the absolute
 * trace scaled by the tower index [K:Q(a)].  Returns 1 with *out set to a fresh owned
 * Integer/Rational; 0 when a (or theta) is not a constant algebraic number; 2 when
 * theta is given but a is not an element of Q(theta); -1 when FLINT is compiled out. */
int flint_qqbar_algebraic_number_trace(const Expr* a, const Expr* theta, Expr** out);

/* Same-field arithmetic for the Plus/Times/Power combination pre-passes. `a`
 * (and `b`) are AlgebraicNumber objects; combining happens only when generators
 * are structurally equal (expr_eq). `r` is an Integer/Rational scalar. */
Expr* flint_qqbar_algnum_add(const Expr* a, const Expr* b);
Expr* flint_qqbar_algnum_mul(const Expr* a, const Expr* b);
Expr* flint_qqbar_algnum_pow(const Expr* a, long p);
Expr* flint_qqbar_algnum_add_rational(const Expr* a, const Expr* r);
Expr* flint_qqbar_algnum_scale_rational(const Expr* a, const Expr* r);

/* Minimal polynomial of the algebraic generator `gen`, as a fresh List of its
 * integer coefficients low-to-high {c0, c1, ..., cn}, or NULL when gen is not a
 * recognised algebraic generator (or FLINT is compiled out).  Backs the
 * field-coefficient Expand fast path in flint_bridge.c. */
Expr* flint_qqbar_gen_minpoly_coeffs(const Expr* gen);

/* Exact reduced row echelon form of an r×c row-major matrix (flat[], row-major)
 * whose entries are rationals and AlgebraicNumber[theta, ...] over ONE common
 * theta, via native FLINT antic nf_elem arithmetic.  The RREF is canonical, so
 * the result equals the classical AlgebraicNumber RowReduce for any Method.
 * Returns a List-of-Lists Expr (caller owns) or NULL to fall through to the
 * classical path (no AlgebraicNumber entry, a second distinct generator, a
 * non-field entry, or no FLINT).  Does not take ownership of flat[]. */
Expr* flint_qqbar_nf_mat_rref(Expr** flat, int r, int c);

/* ------------------------------------------------------------------ */
/*  Exact algebra over the number field generated by the input's atoms */
/* ------------------------------------------------------------------ */
/* Entries are constant algebraic expressions -- rationals combined by
 * Plus / Times / integer Power with algebraic atoms (radicals of rationals,
 * I, Root objects, other constant-algebraic leaves).  The atoms
 * are mapped once into a number field K = Q(phi) (primitive element), every
 * operation runs in native nf_elem arithmetic, and results are rendered in a
 * monomial basis of the original atoms (a + b I + c Sqrt[3] + ...).  Each
 * declines (NULL / 0) when an entry is outside that domain -- a free symbol,
 * an inexact number, no algebraic atom at all (the fmpq paths own that), more
 * than 8 atoms, or [K : Q] > 64 -- and the caller keeps its classical path.
 * None take ownership of their inputs. */

/* Solve A.X = B (A r x c, B r x k, row-major).  1: solved, *out = c*k fresh
 * entries (free variables 0; the caller owns the entries and the array);
 * -1: inconsistent; 0: declined. */
int   flint_qqbar_alg_mat_solve(Expr* const* A, int r, int c, Expr* const* B, int k,
                                Expr*** out);
/* Reduced row echelon form, as a List of Lists. */
Expr* flint_qqbar_alg_mat_rref(Expr* const* flat, int r, int c);
/* Determinant of an n x n matrix. */
Expr* flint_qqbar_alg_mat_det(Expr* const* flat, int n);
/* Inverse of an n x n matrix; NULL with *singular = 1 for a singular one. */
Expr* flint_qqbar_alg_mat_inverse(Expr* const* flat, int n, int* singular);

/* An exact zero test over the field generated by the atoms of `exprs` (each
 * a constant algebraic expression, as above).  NULL to decline.
 * flint_qqbar_alg_field_is_zero: 1 when `e` is zero, 0 when it is not, -1
 * when `e` is outside the field (an atom not among those of `exprs`, a free
 * symbol, a division by zero). */
typedef struct FlintAlgField FlintAlgField;
FlintAlgField* flint_qqbar_alg_field_new(Expr* const* exprs, size_t n);
int   flint_qqbar_alg_field_is_zero(FlintAlgField* A, const Expr* e);
void  flint_qqbar_alg_field_free(FlintAlgField* A);

/* Free the module's static memos (qqbar values, Root-object forms, sorted
 * roots, generator minpolys).  Exit-time only; FLINT builds only. */
void flint_qqbar_release_caches(void);

#endif /* FLINT_QQBAR_H */
