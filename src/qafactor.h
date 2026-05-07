/* qafactor.h — Trager algebraic-factoring helpers (Phase G).
 *
 * This module sits on top of qa.{c,h} (Q(α) elements) and qaupoly.{c,h}
 * (Q(α)[x] polynomials) and connects them to picocas's existing
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
 * resultants we serialise the QAUPoly to a picocas Expr (a polynomial
 * in two free symbols, x_name and y_name) and call the existing
 * `internal_resultant`, which dispatches to the Sylvester-matrix
 * routine in poly.c.  The result comes back as a univariate Expr in
 * x_name, expanded into canonical form. */

#ifndef PICOCAS_QAFACTOR_H
#define PICOCAS_QAFACTOR_H

#include "qa.h"
#include "qaupoly.h"

struct Expr;

/* Build P_α(y) as a picocas Expr (a polynomial in `y_name`).
 * Caller owns the returned Expr. */
struct Expr* qaext_to_expr(const QAExt* ext, const char* y_name);

/* Build f(x, α) → g(x, y) ∈ Q[x, y] as a picocas Expr by textual
 * substitution α → y.  The result is a polynomial in two free
 * symbols `x_name` and `y_name`.  Caller owns the returned Expr. */
struct Expr* qaupoly_to_expr(const QAUPoly* f,
                             const char* x_name,
                             const char* y_name);

/* Compute the field norm
 *
 *     N(f)(x) = Resultant_y(P_α(y), g(x, y))   ∈ Q[x]
 *
 * via picocas's own Resultant builtin.  The result is returned as a
 * picocas Expr (a univariate polynomial in `x_name`, post-Expand).
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

/* Render a QAUPoly over Q(α) as a picocas Expr in `x_name` and
 * `alpha_render` (the surface form for α, e.g. Sqrt[2]).  The result
 * is post-evaluated so that `Sqrt[c]^k` collapses to its canonical
 * form. */
struct Expr* qaupoly_to_expr_alpha(const QAUPoly* f,
                                   const char* x_name,
                                   const struct Expr* alpha_render);

/* Public picocas-level wrapper: factor `poly` ∈ Q(α)[x] over Q(α),
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

/* Public picocas-level wrapper for towers.  Parallel to
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

#endif
