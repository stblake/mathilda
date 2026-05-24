/* qafactor.h — Trager algebraic-factoring helpers (Phase G).
 *
 * This module sits on top of qa.{c,h} (Q(α) elements) and qaupoly.{c,h}
 * (Q(α)[x] polynomials) and connects them to Mathilda's existing
 * polynomial machinery (Resultant, Factor, ZUPoly).
 *
 * Phase G3 — Norm via resultant.
 *
 *   Given f(x, α) ∈ Q(α)[x] and minimal polynomial P_α(y) ∈ Q[y], the
 *   field-norm
 *
 *      N(f)(x) = Resultant_y( P_α(y), g(x, y) )
 *
 *   where g(x, y) ∈ Q[x, y] is f with α textually replaced by y, lies
 *   in Q[x].  When f is irreducible over Q(α), N(f) is a power of an
 *   irreducible Q-factor up to (separable) shifts — this is the
 *   foundation of Trager's algorithm.
 *
 * Implementation strategy.  Rather than reinvent multivariate
 * resultants we serialise the QAUPoly to a Mathilda Expr (a polynomial
 * in two free symbols, x_name and y_name) and call the existing
 * `internal_resultant`, which dispatches to the Sylvester-matrix
 * routine in poly.c.  The result comes back as a univariate Expr in
 * x_name, expanded into canonical form. */

#ifndef MATHILDA_QAFACTOR_H
#define MATHILDA_QAFACTOR_H

#include <stdbool.h>

#include "qa.h"
#include "qaupoly.h"

struct Expr;

/* Build P_α(y) as a Mathilda Expr (a polynomial in `y_name`).
 * Caller owns the returned Expr. */
struct Expr* qaext_to_expr(const QAExt* ext, const char* y_name);

/* Build f(x, α) → g(x, y) ∈ Q[x, y] as a Mathilda Expr by textual
 * substitution α → y.  The result is a polynomial in two free
 * symbols `x_name` and `y_name`.  Caller owns the returned Expr. */
struct Expr* qaupoly_to_expr(const QAUPoly* f,
                             const char* x_name,
                             const char* y_name);

/* Compute the field norm
 *
 *     N(f)(x) = Resultant_y(P_α(y), g(x, y))   ∈ Q[x]
 *
 * via Mathilda's own Resultant builtin.  The result is returned as a
 * Mathilda Expr (a univariate polynomial in `x_name`, post-Expand).
 * Returns NULL if f is the zero polynomial.  Caller owns the returned
 * Expr. */
struct Expr* qaupoly_norm(const QAUPoly* f,
                          const char* x_name,
                          const char* y_name);

/* ============================ Phase G4 ============================ */

/* Result of Trager's sqfr_norm.  Caller owns `g` (qaupoly_free) and
 * `R` (expr_free).  `s == -1` signals failure (e.g. f == 0 or the
 * shift loop did not converge within QA_SQFR_NORM_MAX_TRIES). */
typedef struct QASqfrNormResult {
    int   s;
    QAUPoly*     g;
    struct Expr* R;
} QASqfrNormResult;

/* Maximum shifts to try before giving up.  Trager's theorem guarantees
 * finite termination; in practice s ∈ {0, 1, 2} for nearly all inputs.
 * Public for tests / introspection. */
#ifndef QA_SQFR_NORM_MAX_TRIES
#define QA_SQFR_NORM_MAX_TRIES 32
#endif

/* Trager's squarefree-norm: find smallest s ≥ 0 such that
 *
 *     R(x) = Norm( f(x − s α) )
 *
 * is squarefree over Q.  Returns the shift, the shifted polynomial
 * g(x) = f(x − sα) ∈ Q(α)[x], and the squarefree norm R(x) ∈ Q[x].
 *
 * Squarefree-detection is gcd(R, R') over Q[x] via the existing
 * PolynomialGCD builtin: R is squarefree iff that gcd has degree 0
 * in x. */
QASqfrNormResult qa_sqfr_norm(const QAUPoly* f,
                              const char* x_name,
                              const char* y_name);

/* Free the owned fields of a sqfr_norm result.  No-op for failure
 * results.  Does NOT free the struct itself (it is a value type). */
void qa_sqfr_norm_result_free(QASqfrNormResult* r);

/* Trager's algebraic-factoring of squarefree f ∈ Q(α)[x].
 *
 * Returns a freshly-allocated array of `*n_out` monic irreducible
 * factors over Q(α).  The product of these factors equals f up to
 * a leading-coefficient unit (which equals 1 when f itself is monic).
 *
 * For irreducible f returns an array of length 1 containing a copy
 * of f (made monic).
 *
 * Returns NULL with *n_out = 0 on failure (zero input, sqfr_norm
 * non-convergence, or arithmetic failure in qaupoly_gcd).
 *
 * Ownership: caller must qaupoly_free each entry and free() the
 * outer array. */
QAUPoly** qa_alg_factor(const QAUPoly* f,
                        const char* x_name,
                        const char* y_name,
                        int* n_out);

/* ============================ Phase G5 ============================ */

/* Resolve an `Extension -> α` argument to a (QAExt, surface-render-form)
 * pair.
 *
 * Recognised forms:
 *   Sqrt[c]         (rational c, |c| ≤ INT64 / 2)   →  P_α = y² − c
 *   c^(1/n)         (rational c, integer n ≥ 2)     →  P_α = yⁿ − c
 *   I                                                →  P_α = y² + 1
 *
 * On success returns a fresh QAExt (caller owns: qaext_free) and writes
 * a fresh `*render_out` Expr (caller owns: expr_free).  The render-form
 * is what α should look like in user-visible output — typically a copy
 * of `alpha_expr` itself.
 *
 * Returns NULL on unsupported or non-rational inputs; *render_out is
 * left untouched in that case. */
QAExt* qa_resolve_extension(const struct Expr* alpha_expr,
                            struct Expr** render_out);

/* Render a QAUPoly over Q(α) as a Mathilda Expr in `x_name` and
 * `alpha_render` (the surface form for α, e.g. Sqrt[2]).  The result
 * is post-evaluated so that `Sqrt[c]^k` collapses to its canonical
 * form. */
struct Expr* qaupoly_to_expr_alpha(const QAUPoly* f,
                                   const char* x_name,
                                   const struct Expr* alpha_render);

/* Lift a polynomial Expr in `var` (whose coefficients may contain the
 * surface form `alpha_render`, e.g. Sqrt[2]) to a QAUPoly over `ext`.
 * Returns NULL when any coefficient is outside Q(α) (e.g. contains a
 * free symbol other than `var` and `alpha_render`).  Caller owns the
 * returned QAUPoly via qaupoly_free.
 *
 * Used by Extension-aware GCD / LCM / Cancel / Together / Apart paths
 * that need direct access to qaupoly_gcd / qaupoly_divrem rather than
 * the higher-level qa_factor_with_extension wrapper. */
QAUPoly* qa_expr_to_qaupoly(const struct Expr* poly,
                            const struct Expr* var,
                            const struct Expr* alpha_render,
                            const QAExt* ext);

/* Public Mathilda-level wrapper: factor `poly` ∈ Q(α)[x] over Q(α),
 * where `alpha_expr` selects the extension (per qa_resolve_extension)
 * and `var` is the polynomial indeterminate.
 *
 * The input may itself contain `alpha_expr` as a sub-expression (e.g.
 * `Factor[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]`); occurrences
 * are recognised structurally and lifted into the α-component of the
 * QAUPoly coefficients.
 *
 * Returns a fresh Expr — typically `Times[h_1, h_2, ...]` of the monic
 * factors, prefixed by the leading-coefficient unit when it is not 1.
 * Returns NULL on unsupported inputs (non-Q(α) coefficient, sqfr_norm
 * non-convergence, multiple variables, etc.); the caller should fall
 * back to non-extension factoring in that case. */
struct Expr* qa_factor_with_extension(const struct Expr* poly,
                                      const struct Expr* alpha_expr,
                                      const struct Expr* var);

/* ============================ Phase G6 ============================ */

/* QATower — compositum extension Q(α_1, ..., α_n) reduced to a single
 * primitive element γ via Trager §3.  Built iteratively: start with
 * Q(α_1), then absorb α_2, α_3, ..., each step finding the smallest
 * shift s_i ≥ 0 such that γ_i = γ_{i-1} + s_i α_i is a primitive
 * element of the new compositum.
 *
 * Invariants:
 *   - `ext` is the (monic, irreducible-over-Q) minimal polynomial of γ.
 *   - `alphas[i]` is α_i expressed as a Q(γ)-element (a polynomial in γ
 *     with rational coefficients).  All `alphas[i]` share `ext`.
 *   - `alpha_renders[i]` is the user-visible surface form of α_i
 *     (e.g. Sqrt[2], 2^(1/3), I).
 *   - `gamma_render` is the surface form for γ — accumulated as
 *     α_1 + s_2 α_2 + s_3 α_3 + ... + s_n α_n.  Used only for output
 *     rendering (substituted into the internal placeholder symbol). */
typedef struct QATower {
    QAExt* ext;
    int n;
    QANum** alphas;
    struct Expr** alpha_renders;
    struct Expr* gamma_render;
} QATower;

/* Resolve a list of n ≥ 1 algebraic generators to a tower.  Each
 * α_i must be recognised by qa_resolve_extension.  Iteratively builds
 * the compositum via Trager's primitive-element algorithm; returns
 * NULL on non-convergence or if any α_i is unrecognised.  The single-
 * generator case (n == 1) is delegated to `qa_resolve_extension` and
 * wrapped as a one-element tower for symmetry.  Caller owns the
 * returned tower (qa_tower_free). */
QATower* qa_resolve_extension_tower(struct Expr* const* alpha_exprs, int n);
void     qa_tower_free(QATower* t);

/* Public Mathilda-level wrapper for towers.  Parallel to
 * qa_factor_with_extension but accepts a list of α-generators
 * (`alpha_exprs[0..n_alphas-1]`).  For n_alphas == 1 this is exactly
 * qa_factor_with_extension; for n_alphas ≥ 2 it builds the compositum
 * Q(γ), substitutes each α_i in the input with its Q(γ)-representation,
 * and dispatches to the same Trager core as G5.
 *
 * Returns NULL on resolution / shift-search / lift failures (caller
 * falls back to non-extension factoring). */
struct Expr* qa_factor_with_extension_tower(const struct Expr* poly,
                                            struct Expr* const* alpha_exprs,
                                            int n_alphas,
                                            const struct Expr* var);

/* Tower-aware Cancel: combine `arg` into a single fraction over Q(γ),
 * cancel num/den by their qaupoly GCD, and render with γ's surface
 * form (`t->gamma_render`).
 *
 * Steps (mirrors `qa_factor_with_extension_tower`'s substitution pattern):
 *   1. Substitute each α_i in `arg` with its Q(γ) representation
 *      (a polynomial in QA_ALPHA_INTERNAL).
 *   2. Run `Together` (no extension) to combine sum-of-fractions into a
 *      single fraction.
 *   3. Lift numerator and denominator to QAUPoly over Q(γ) and divide
 *      by their `qaupoly_gcd`.
 *   4. Render back through `qaupoly_to_expr_alpha` with γ's surface form.
 *
 * Returns NULL on any lift / GCD failure (caller falls back to the
 * non-tower path).  Caller owns the returned Expr.  Use when
 * `extension_autodetect` returns a tower with n ≥ 2 generators (the
 * single-generator case is handled by the existing
 * `cancel_with_extension` and `together_recursive_ext`). */
struct Expr* qa_cancel_with_tower(const struct Expr* arg, const QATower* t);

/* Multi-α PolynomialGCD over the compositum Q(γ).
 *
 * Substitutes each α_i in every input with its polynomial-in-γ form,
 * lifts each input to QAUPoly[x] over t->ext, computes the iterated
 * `qaupoly_gcd`, and renders the result with γ's surface form
 * (`t->gamma_render`).  Requires the substituted inputs to share
 * exactly one polynomial variable besides γ (the multivariate case is
 * outside the tier-1 scope).  Returns NULL on substitution / lift /
 * GCD failures so the caller can fall back to the no-extension
 * builtin_polynomialgcd path.
 *
 * Used by `builtin_polynomialgcd` when `extension_autodetect_args`
 * returns a tower with n ≥ 2.  The single-generator case is handled
 * by the existing `polynomialgcd_with_extension`. */
struct Expr* qa_polynomialgcd_with_tower(struct Expr* const* argv,
                                         size_t argc,
                                         const QATower* t);

/* Multi-α PolynomialLCM over the compositum Q(γ).  Folds left-to-right
 * via `lcm(a, b) = a * b / gcd(a, b)`, using `qa_polynomialgcd_with_tower`
 * internally.  Returns NULL on the same failure modes. */
struct Expr* qa_polynomiallcm_with_tower(struct Expr* const* argv,
                                         size_t argc,
                                         const QATower* t);

/* Phase D: multivariate tower-based PolynomialGCD/PolynomialLCM via
 * substitute-back through `internal_polynomialgcd` / `internal_polynomiallcm`.
 *
 * Each α_i in the inputs is replaced by its polynomial-in-γ form
 * (using QA_ALPHA_INTERNAL as the γ symbol), the no-extension
 * multivariate builtin is invoked treating γ as a polynomial variable,
 * and the result is mapped back through `t->gamma_render` with Expand
 * + evaluate to canonicalise.
 *
 * Correctness note: the GCD computed this way is the
 * `Q[γ, x, y, ...]`-GCD, which is a (possibly non-maximal) Q(γ)-divisor
 * of both inputs.  When the canonical Q(γ)[x,y,...]-GCD has γ-degree
 * < deg(γ_min) and the inputs are γ-primitive, the two definitions
 * coincide.  Otherwise this returns a result that is mathematically a
 * common divisor / common multiple but not necessarily in canonical
 * form — downstream Cancel / Together passes can reduce further.
 *
 * Used by `builtin_polynomialgcd` and `builtin_polynomiallcm` as a
 * fallback when the univariate `qa_polynomialgcd_with_tower` /
 * `qa_polynomiallcm_with_tower` path bails on having more than one
 * polynomial variable besides γ.  Returns NULL when the substituted
 * form does not simplify (the no-extension call returned unchanged). */
struct Expr* qa_polynomialgcd_with_tower_multivar(struct Expr* const* argv,
                                                  size_t argc,
                                                  const QATower* t);
struct Expr* qa_polynomiallcm_with_tower_multivar(struct Expr* const* argv,
                                                  size_t argc,
                                                  const QATower* t);

/* Predicate: true when any tower generator has a non-integer base, i.e.
 * surfaces as `Sqrt[non_int]` or `Power[non_int, p/q]` (a nested radical
 * like `Sqrt[5 + 2 Sqrt[6]]`).
 *
 * Why callers care: `qa_cancel_with_tower`'s Step 1 substitutes each
 * α_i with its polynomial-in-γ form, but a nested-radical α_i was
 * originally surfaced from a `Power[Plus[...], p/q]` (q ≥ 2) node in
 * the input, and the substitution leaves that `Power` opaque — its
 * base becomes a polynomial in γ but the `Power[base, p/q]` itself
 * remains structurally a `Power` with a non-integer exponent that
 * Step 4's no-extension `Together` cannot combine. The result inflates
 * past the leaf-count gate and is rejected. Skipping the tower path
 * (and the N×fallback) when this returns true avoids the wasted
 * substitution + Together + qaupoly-lift + rejection cycle. */
bool qa_tower_has_nested_radical(const QATower* t);

/* Syntactic predicate on an EXPRESSION (not a tower): true when `e`
 * contains a nested radical (`Sqrt[non_integer]` or
 * `Power[non_integer, p/q]` with q ≥ 2 whose radicand itself contains
 * a radical) AND `e` has no free polynomial variable.
 *
 * Callers use this as a fast O(input-size) prefilter so they can skip
 * the much-more-expensive `extension_autodetect` + tower-build cascade.
 * Two conditions must both hold:
 *   1. Nested radical present — the γ-substitution can't reduce the
 *      outer `Power[non_int, p/q]` form (Step 4 leaves it opaque, the
 *      leaf-count gate rejects the result).
 *   2. No free polynomial variable — qaupoly_gcd over Q(γ)[x] needs a
 *      polynomial variable to find non-trivial common factors; without
 *      one, the tower's cancellation work degenerates to constant
 *      arithmetic in Q(γ) that the input's nested radicals prevent
 *      from completing usefully.
 *
 * When only one of the two holds (e.g. nested radicals with free
 * vars — the `D[Integrate[a x/(x^3+2), x], x]` headline case), the
 * tower path IS needed and returns true would suppress it. */
bool expr_has_nested_radical_radicand(const struct Expr* e);

/* ============================ Phase G9 ============================ */
/* Automatic algebraic-extension detection.                            */

/* Maximum number of distinct rational-base algebraic generators that
 * the auto-detector will collect from a single expression.  More than
 * this and the auto-detector bails (returns NULL), because the tower
 * construction has cost exponential in the generator count and the
 * input is almost certainly outside the intended domain of use. */
#ifndef QA_AUTODETECT_MAX_GENS
#define QA_AUTODETECT_MAX_GENS 4
#endif

/* Walk `e` collecting every algebraic-number generator implicit in the
 * expression and build the corresponding compositum.  Recognised
 * generator shapes (tier 1):
 *
 *   Power[c, p/q]        c integer, |c| ≥ 2, gcd(p, q) = 1, q ≥ 2
 *   Sqrt[c]              c integer, |c| ≥ 2
 *
 * Multiple occurrences sharing the same integer base `c` are merged
 * into a single generator `Power[c, 1/lcm(q_i)]`, so e.g. an expression
 * containing both `Power[2, 1/3]` and `Power[2, 1/2]` resolves to a
 * single generator 2^(1/6) of degree 6.
 *
 * Returns NULL when:
 *   - the expression contains no algebraic generators (caller stays
 *     in Q[x]);
 *   - the expression contains a `Power[u, p/q]` (q > 1) whose base u
 *     is not an integer literal (e.g. a rational, a polynomial, or a
 *     nested radical);
 *   - more than QA_AUTODETECT_MAX_GENS distinct integer bases appear;
 *   - the tower construction in `qa_resolve_extension_tower` fails.
 *
 * Caller owns the returned tower (free with `qa_tower_free`). */
QATower* extension_autodetect(const struct Expr* e);

/* Multi-expression variant of `extension_autodetect`: scan every entry
 * of `args[0..argc-1]` and build the joint compositum.  Equivalent in
 * effect to wrapping the inputs in a single `List[...]` and calling
 * `extension_autodetect`, but skips the temporary Expr.  Used by
 * polynomial builtins that take multiple polynomial arguments
 * (PolynomialGCD, PolynomialLCM, PolynomialQuotient, PolynomialRemainder).
 *
 * Returns NULL on any of the conditions documented for
 * `extension_autodetect`. */
QATower* extension_autodetect_args(struct Expr* const* args, size_t argc);

/* Phase E: single-generator polynomial-radicand Cancel/Together path.
 *
 * Handles inputs containing exactly one distinct radical of the form
 * `Sqrt[poly]` or `Power[poly, 1/q]` where `poly` is a polynomial
 * expression with free symbols.  Substitutes the radical with a fresh
 * symbol, runs `Together`, reduces numerator and denominator modulo
 * `S^q - radicand` via `PolynomialRemainder`, and substitutes back.
 *
 * Used as a fallback by `builtin_cancel` / `builtin_together` when the
 * standard `extension_autodetect` path returns NULL (which it does for
 * polynomial-radicand radicals — the QAExt machinery has Q-coefficient
 * minimal polynomials and cannot represent the resulting Q(params)-
 * coefficient case).
 *
 * Returns NULL when:
 *   - No polynomial-radicand radical is present.
 *   - More than one distinct polynomial-radicand radical is present
 *     (multi-radical inputs need Groebner-basis reasoning that this
 *     path does not attempt — see the Cardano gap documented in the
 *     2026-05 changelog).
 *   - `Power[poly, p/q]` with reduced p != 1.
 *   - The post-reduction result is identical to the input (no
 *     simplification fired) or has inflated past a sanity gate.
 *
 * Caller owns the returned Expr (free with `expr_free`). */
struct Expr* qa_cancel_with_poly_radical(const struct Expr* arg);

#endif
