/*
 * flint_bridge.h
 * --------------
 * Boundary between Mathilda's Expr trees and FLINT's polynomial types, for the
 * fast/rigorous algebraic-extension arithmetic engine described in
 * ALGEBRAIC_EXTENSION_ARITHMETIC_PLAN.md.
 *
 * The whole module is a no-op when Mathilda is built without FLINT
 * (USE_FLINT undefined): the entry points still exist and link, but
 * flint_bridge_available() returns 0 and the conversion/operation routines
 * return NULL, so callers fall back to the classical (slower but rigorous)
 * path. This mirrors the USE_MPFR / USE_LAPACK / USE_GRAPHICS graceful-degrade
 * policy.
 *
 * M1 scope (this file): the rational multivariate case R = Q[x_1..x_n] only,
 * exercised end-to-end through a scaffolding builtin so the bridge can be
 * verified before any real consumer (Cancel/Together/Factor) is re-pointed.
 * Number fields (gr/ANTIC) and the parametric Q(t)(alpha) outer loop arrive in
 * M2/M3.
 */
#ifndef FLINT_BRIDGE_H
#define FLINT_BRIDGE_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 1 if this build links FLINT (>= 3.0), 0 otherwise. */
int flint_bridge_available(void);

/*
 * Multivariate GCD over Q.
 *
 * Recognises `a` and `b` as elements of Q[x_1..x_n] (variables = the free
 * symbols occurring in them; coefficients integer or Rational), computes their
 * GCD with FLINT's fmpq_mpoly_gcd, and renders the (FLINT-monic) result back to
 * a fresh Expr the caller owns.
 *
 * Returns NULL — leaving the caller to use the classical path — when FLINT is
 * absent, when either input is not a polynomial over Q in recognisable
 * variables (e.g. contains an inexact real, a negative/symbolic power, or an
 * unrecognised head), or on any internal failure. Never mutates its arguments.
 */
Expr* flint_multivariate_gcd(const Expr* a, const Expr* b);

/*
 * As flint_multivariate_gcd, but rescales the result to the primitive-integer,
 * positive-leading associate that Mathilda's classical PolynomialGCD path
 * returns (Gauss's lemma: content(gcd) = gcd(content a, content b)). This makes
 * it a transparent drop-in for the two-argument PolynomialGCD that Cancel /
 * Together / Simplify issue, with identical output but far greater speed on
 * multivariate polynomials over Q. Same NULL-fallback contract as above.
 */
Expr* flint_multivariate_gcd_normalized(const Expr* a, const Expr* b);

/* Exact division a / b over Q[x_1..x_n]: the quotient when b divides a exactly,
 * else NULL (inexact, or out of scope — parametric / algebraic coefficients).
 * Fast and overflow-free (fmpq_mpoly_divides); bypasses the classical
 * pseudo-division that is quadratic-plus on dense, large-coefficient inputs. */
Expr* flint_multivariate_divexact(const Expr* a, const Expr* b);

/*
 * Fully expand a polynomial over Q[x_1..x_n] to its distributed normal form
 * (Plus of Times[coeff, x_i^e_i, ...] monomials), via fmpq_mpoly.  This is the
 * accelerator behind `Expand`: converting the tree to fmpq_mpoly does the
 * multiplication and collection in packed FLINT arithmetic — orders of
 * magnitude faster than the generic Expr schoolbook multiply on dense,
 * high-degree inputs (e.g. (1 + x + 5 x^3 + 8 x^17)^341).
 *
 * Returns a fresh Expr the caller owns, or NULL — leaving the caller on the
 * classical path — when FLINT is absent, when `e` is purely numeric (nothing to
 * expand), or when `e` is not a polynomial over Q in recognisable variables (an
 * inexact real, a negative/fractional/symbolic power, or a transcendental head
 * such as Sin/Log).  Never mutates its argument.  The caller is responsible for
 * bounding the result size beforehand; this routine will build whatever the
 * input demands. */
Expr* flint_expand_polynomial(const Expr* e);

/* True (1) when `e` is a polynomial over Q in recognisable variables — i.e. the
 * kind of expression flint_expand_polynomial accepts (built from integers,
 * rationals, symbols, Plus/Times, and non-negative integer Power). Numeric
 * constants qualify. Lets a caller partition a product into the factors FLINT
 * can multiply and the rest. Always 0 without FLINT. */
int flint_is_polynomial_over_q(const Expr* e);

/*
 * Native base-field Risch differential-equation solver for the exponential
 * tower, over Q(xvar): solves  Dq + f q = g  (D = d/d(xvar)) for the solution
 * y, where f is a polynomial (f = i·u' for e^u, u a polynomial in xvar, so
 * WeakNormalizer is trivial) and g is a rational function over Q(xvar). Runs
 * RdeNormalDenominator + RdeBoundDegreeBase + Bronstein's SPDE ladder +
 * PolyRischDENoCancel entirely in fmpq_poly — f and g are converted straight to
 * fmpq_poly, so NONE of the base-field arithmetic touches the generic evaluator
 * (this is where the seconds of Together/Cancel on a degree-2n denominator go).
 * Returns:
 *    1  solved — *y_out = the reduced rational solution y as an Expr,
 *    0  no solution of bounded degree exists (authoritative decline),
 *   -1  out of scope — g not univariate over Q, or f not a polynomial (the
 *       primitive/log tower, where weak normalization is non-trivial) — the
 *       caller falls back to the Expr RDE path.
 * Without FLINT this always returns -1.
 */
int flint_rde_base_solve_fg(const Expr* f, const Expr* g,
                            const char* xvar, Expr** y_out);

/*
 * Univariate GCD over a real quadratic number field Q(sqrt d) (M2, p=0, r=1).
 *
 * Recognises `a`, `b` as univariate polynomials in a single variable x whose
 * coefficients lie in Q(sqrt d) for one square-root generator `Sqrt[d]` (d a
 * positive non-square integer). Computes the GCD rigorously in the field —
 * using the relation d = (sqrt d)^2, so e.g. gcd(x^2 - 2, x - Sqrt[2]) is the
 * true x - Sqrt[2], a split the free-variable / rational view cannot see — via
 * FLINT's generic-ring number-field layer (gr_ctx_init_nf + gr_poly_gcd), and
 * renders the FLINT-monic result back to an Expr with coefficients p + q Sqrt[d].
 *
 * Returns NULL (caller falls back) when there is no radical (the pure-rational
 * case — use flint_multivariate_gcd), when more than one distinct radical or
 * more than one variable appears (towers / multivariate-over-a-field are later
 * milestones), or on any out-of-scope construct. Never mutates its arguments.
 */
Expr* flint_numberfield_gcd(const Expr* a, const Expr* b);

/*
 * Univariate GCD over a cyclotomic field Q(zeta_n) (M2, p=0, r=1).
 *
 * Recognises roots of unity `(-1)^(p/q) = Power[-1, Rational[p,q]]` in the
 * coefficients (one variable), works in Q(zeta_N) with N = 2*lcm(q) and minimal
 * polynomial Phi_N, and renders the monic GCD back with coefficients in powers
 * of the generator (-1)^(1/Q). Returns NULL when there is no root of unity, when
 * a radical also appears (a tower — later work), the imaginary unit `Complex[0,1]`
 * appears (deferred), or for any other out-of-scope construct.
 */
Expr* flint_cyclotomic_gcd(const Expr* a, const Expr* b);

/*
 * Univariate GCD over a tower of square-root extensions Q(√d_1,…,√d_r), r ≥ 2
 * (M2). The tower is collapsed to a single primitive element θ = Σ c_i √d_i by
 * rational linear algebra (each √d_i acts as a known matrix in the product
 * basis), the GCD is computed in the absolute field Q(θ), and the result is read
 * back into the product basis (single radicals Sqrt[∏ d_i]). Returns NULL when
 * fewer than two distinct radicals appear, when the radicals are not independent
 * (no primitive θ of full degree 2^r), or for any out-of-scope construct.
 */
Expr* flint_tower_gcd(const Expr* a, const Expr* b);

/*
 * Dispatch helpers for consumers (PolynomialGCD, Cancel, Together, …).
 *
 * flint_extension_gcd: tries Q(√d) → Q(ζ_n) → tower, returning NULL when no
 *   algebraic generator is present — so a caller keeps its classical path for
 *   plain Q[x] and parametric inputs, and only the under-reduced extension cases
 *   are taken over. This is the recommended consumer entry point.
 * flint_polynomial_gcd: the full chain, plain rational multivariate first.
 *
 * Both return NULL when FLINT is not compiled in. Neither mutates its arguments;
 * the returned Expr is owned by the caller.
 */
Expr* flint_extension_gcd(const Expr* a, const Expr* b);
Expr* flint_polynomial_gcd(const Expr* a, const Expr* b);

/*
 * Plain multivariate polynomial operations over Q[x_1..x_n] (no algebraic
 * extension): the rational analogues of the extension helpers above, exposed
 * for the System` builtins (PolynomialGCD/Resultant/Factor) and the FLINT`
 * context wrappers. Each returns a fresh Expr the caller owns, or NULL when the
 * input is not a polynomial over Q in recognisable variables / out of scope /
 * FLINT absent, so the caller keeps its classical path.
 *
 *   flint_polynomial_resultant: Res_x(a, b) with `var` the elimination
 *     variable (a symbol); all other free symbols are coefficients.
 *   flint_polynomial_factor / _squarefree: irreducible / squarefree
 *     factorisation over Q, rendered as Times[const, base^exp, ...].
 *
 * (GCD is flint_multivariate_gcd / flint_polynomial_gcd above.) These emit
 * FLINT-canonical output (monic factors/GCD over Q); the System`-builtin
 * callers renormalise to Mathematica's primitive-over-Z convention.
 */
Expr* flint_polynomial_resultant(const Expr* a, const Expr* b, const Expr* var);
Expr* flint_polynomial_factor(const Expr* p);
Expr* flint_polynomial_factor_squarefree(const Expr* p);

/*
 * Univariate integer-polynomial factorisation via FLINT (fmpz_poly_factor).
 * Recognises `p` as a polynomial in exactly one variable (detected with FLINT's
 * own bignum-safe variable scan, so it is robust to large coefficients) whose
 * coefficients are all integers (EXPR_INTEGER or EXPR_BIGINT). Returns the
 * factorisation as an evaluated Times[content, factor^exp, ...] Expr the caller
 * owns, in Mathematica's convention: each irreducible factor is primitive with
 * positive leading coefficient (so 4x^2-9 -> (2x-3)(2x+3), NOT the monic
 * 4(x-3/2)(x+3/2) an fmpq factorisation would give) and the signed integer
 * content is emitted as a separate factor. Returns NULL when `p` is not a
 * univariate polynomial over Z (multivariate, rational/symbolic coefficient,
 * a denominator, degree < 1) or FLINT is absent, so the caller keeps its
 * classical path. This is the transparent fast path for Factor[] on univariate
 * integer polynomials — and it fixes a classical-pipeline bug where a bignum
 * coefficient is misclassified as a variable and its term dropped.
 */
Expr* flint_univariate_factor_auto(const Expr* p);

/*
 * Univariate GCD over a parametric radical field Q(t_1..t_p)(sqrt k) where the
 * radicand k is itself a free symbol (parameter), not an integer — the
 * `p >= 1, r = 1, deg = 2` regime, i.e. the Goursat blocker ring Q(a,b,k)(sqrt k).
 * Since sqrt(k) is transcendental over Q(t) with k = (sqrt k)^2, the field is
 * Q(t_1..t_p, sqrt k), so the GCD reduces to an ordinary multivariate GCD over Q
 * after substituting sqrt(k) -> S, k -> S^2 (fmpq_mpoly_gcd), then reading back.
 * Returns NULL when there is no symbolic radical, when a second radical / cube
 * root / root of unity also appears, or for any out-of-scope construct.
 */
Expr* flint_parametric_sqrt_gcd(const Expr* a, const Expr* b);

/*
 * Resultant Res_var(a, b) over a parametric radical field Q(t_1..t_p)(sqrt k),
 * k a free symbol. Via the same sqrt(k) -> S, k -> S^2 collapse used by the GCD,
 * so FLINT's fmpq_mpoly_resultant does the work. This is what makes the
 * Rothstein-Trager resultant Res_x(A - t D', D) over Q(a,b,k)(sqrt k) tractable
 * (the classical symbolic subresultant PRS blows up on the radical
 * coefficients). Returns NULL when there is no symbolic radical, when `var` is
 * the radicand, or for any out-of-scope construct.
 */
Expr* flint_parametric_sqrt_resultant(const Expr* a, const Expr* b,
                                      const Expr* var);

/*
 * Factor / squarefree-factor a polynomial over a parametric radical field
 * Q(t_1..t_p)(sqrt k), k a free symbol. Via the sqrt(k) -> S, k -> S^2 collapse
 * the field is the rational function field Q(t.., sqrt k), so factoring reduces
 * to ordinary multivariate factoring over Q (fmpq_mpoly_factor) — FLINT does
 * this, unlike the constant-radicand number-field case. Returns a Times[...] of
 * factors (owned by the caller), or NULL when no symbolic radical is present or
 * the input is out of scope.
 */
Expr* flint_parametric_sqrt_factor(const Expr* p);
Expr* flint_parametric_sqrt_factor_squarefree(const Expr* p);

/*
 * Extended GCD of two univariate polynomials a, b in `var` whose coefficients
 * lie in a parametric radical field Q(t_1..t_p)(sqrt k) (k a free symbol).
 * Represents Q(t.., sqrt k)[var] as gr_poly over FLINT's rational-function ring
 * fmpz_mpoly_q (sqrt k -> S, so the field is Q(t.., S)); gr_poly_xgcd does the
 * work — fast where the classical Euclidean xgcd over the radical field blows
 * up. Returns the Mathematica-shaped result List[g, List[u, v]] with
 * u*a + v*b == g (g monic over the field), or NULL when there is no symbolic
 * radical / var is the radicand / the input is out of scope.
 */
Expr* flint_parametric_field_xgcd(const Expr* a, const Expr* b, const Expr* var);

/*
 * Quotient/remainder of a divided by b in `var` over a parametric radical field
 * Q(t..)(sqrt k) — same gr_poly-over-fmpz_mpoly_q representation as the xgcd.
 * Returns List[quotient, remainder] with a == quotient*b + remainder, or NULL
 * when there is no symbolic radical / var is the radicand / out of scope.
 */
Expr* flint_parametric_field_divrem(const Expr* a, const Expr* b, const Expr* var);

/*
 * Cancel/Together a whole rational function `e` over a parametric radical field
 * Q(t..)(sqrt k) via fmpz_mpoly_q (stores it in lowest terms automatically).
 * Handles a Plus of fractions / nested fractions in one shot — the case the
 * num/den-extraction path drops to the slow QA path on. Returns reduced num/den,
 * or NULL out of scope.
 */
Expr* flint_parametric_field_normalize(const Expr* e);

/*
 * Together for a plain rational function over Q: combine into a single reduced
 * fraction via fmpz_mpoly_q, the fast path for the ordinary rational-in-symbols
 * case the algebraic/parametric normalizers decline. Returns the (expanded,
 * reduced) num/den Expr — matching the classical Together output — or NULL when
 * the input has no denominator to combine (a product/polynomial Together leaves
 * factored) or is not a plain rational over Q (symbolic/fractional power,
 * transcendental kernel). NULL without FLINT.
 */
Expr* flint_rational_together(const Expr* e);

/*
 * Cancel for a plain rational function over Q: reduce a single fraction's
 * num/den to lowest terms via fmpz_mpoly_q. Like flint_rational_together but
 * with Cancel's stricter gate — fires only when a denominator is present and NO
 * denominator sits inside a Plus (Cancel, unlike Together, leaves a sum of
 * fractions uncombined). Returns the (expanded, reduced) num/den Expr matching
 * classical Cancel, or NULL out of scope. NULL without FLINT.
 */
Expr* flint_rational_cancel(const Expr* e);

/*
 * FLINT-native partial fraction decomposition over Q. Given a proper rational
 * R / (C * prod_i bases[i]^{ks[i]}) with the bases[] the distinct irreducible
 * factors of the denominator (Exprs), ks[] their multiplicities, C the numeric
 * content, and var the fraction variable name, returns the fractional part
 * sum_i sum_{j=1}^{ks[i]} A_{ij}/bases[i]^j (deg A_{ij} < deg bases[i]) as an
 * unsimplified Expr — the caller adds the polynomial part and lets the evaluator
 * canonicalise. Output matches the classical RowReduce Apart exactly. Computed
 * entirely in fmpq_poly (distinct-factor CRT split + p-adic expansion), so it
 * replaces the O(S^2+) symbolic Gaussian elimination for the pure-Q case.
 * Returns NULL when any operand carries another symbol / radical / fractional
 * power (the multivariate case keeps the classical path). NULL without FLINT.
 */
Expr* flint_apart_over_q(const Expr* R, const Expr* const* bases,
                         const int64_t* ks, int m,
                         const Expr* C, const char* var);

/*
 * Exact division a/b over the detected extension field. Returns the quotient, or
 * NULL when there is no algebraic generator or b ∤ a exactly (so Cancel's
 * divide-back can fall through to its classical path). NULL without FLINT.
 */
Expr* flint_extension_divexact(const Expr* a, const Expr* b);

/*
 * Reduce a rational expression over a genuine-algebraic parametric tower
 * K = Q(t_1..t_p)(alpha_1..alpha_r), where each generator alpha_l is a radical
 * Power[base, p/q] of ANY index q>=2 whose radicand carries a free symbol (cube
 * roots, polynomial / symbol radicands like (x(1-x)(1-k x))^(1/3), k^(1/3)) or a
 * root of unity (-1)^(p/q). K is modelled as Q[params,gens]/I with I generated by
 * the monic minimal polynomials (gen^q - radicand, Phi_{2Q}); ordering the
 * generators as the leading LEX variables makes them a Groebner basis, so
 * fmpz_mpoly_divrem_ideal gives the canonical normal form. Returns the integer 0
 * when the expression is identically zero in K — the rigorous route by which
 * Together/Cancel/Simplify collapse a verified-zero algebraic identity (e.g.
 * D[Integrate[f],x]-f) WITHOUT any numeric zero oracle — and NULL otherwise, so
 * callers keep their classical path for everything else. NULL without FLINT.
 * This is the parametric-algebraic regime the sqrt-only path (flint_parametric_*)
 * and the constant-radicand number-field path both reject. Never mutates `e`.
 */
Expr* flint_algebraic_field_normalize(const Expr* e);

/*
 * Canonical (rationalised) form of a rational expression over the same
 * genuine-algebraic tower K as flint_algebraic_field_normalize, computed by
 * inverting the denominator in K: K is a finite-dimensional Q(params)-vector
 * space, multiplication-by-denominator is a linear map, and solving that system
 * over Q(params) (gr fmpz_mpoly_q) yields N/D with the denominator rationalised
 * to Norm_{K/Q(params)}(D) — free of every radical/root-of-unity generator, the
 * numerator a reduced polynomial in the generators. Rigorous and complete when
 * the minimal polynomials are irreducible. This is the engine behind RootReduce.
 * Returns NULL when `e` has no algebraic generator, is out of scope, is a genuine
 * 0/0, or the multiplication map is singular. NULL without FLINT. Never mutates.
 */
Expr* flint_algebraic_field_canonical(const Expr* e);

/*
 * Combine an expression over the tower into a single fraction with the radicals
 * treated as free kernels (WL-faithful Together/Cancel — keeps radicals in the
 * denominator, does NOT rationalise): each Power[base, p/q] becomes gen^p and
 * the whole expression is put over a common denominator in lowest terms via
 * fmpz_mpoly_q, so radical common factors exposed by the p/q -> gen^p mapping
 * cancel. Generalises the sqrt-only flint_parametric_field_normalize to cube and
 * higher roots. Fires only when every radical generator is genuine (symbol/poly
 * radicand; a constant radical like Sqrt[2] is deferred to the number-field
 * path); roots of unity allowed. Returns NULL out of scope / without FLINT.
 */
Expr* flint_algebraic_field_together(const Expr* e);

/*
 * Unified reduction for the ratcanon rewrite (RATCANON_REWRITE_PLAN.md, Phase 3).
 * Reduce `frac` (num/den in generator symbols; transcendental generators are
 * free variables, the `n_alg` algebraic generators `alg_syms[i]` carry the
 * relation `relations[i] == 0`) to WL-faithful lowest terms: combine via
 * fmpz_mpoly_q, reduce num/den modulo the relation ideal (fmpz_mpoly_divrem_ideal,
 * algebraic generators ordered leading), re-canonicalise. Fully reduces the
 * free/transcendental and sum-of-conjugates algebraic cases; does NOT rationalize
 * and does NOT do the pre-formed-fraction field GCD (caller delegates). Returns a
 * fresh Expr in generator symbols, or NULL out of scope / FLINT absent. Never
 * mutates its arguments.
 */
Expr* flint_tower_reduce(const Expr* frac, const char* const* alg_syms,
                         const Expr* const* relations, int n_alg);

/*
 * Callback invoked once per monomial term by flint_linear_system_terms.
 *   base_exp[0..nbase-1] — exponents in the base variables (the row key).
 *   col                  — 0-based unknown column (0..nunk-1), or nunk for
 *                          the unknown-free (constant / RHS) part of the term.
 *   coeff                — the term's rational coefficient (borrowed; valid
 *                          only for the duration of the call).
 */
typedef void (*flint_lsys_term_fn)(const long* base_exp, int nbase,
                                   int col, const mpq_t coeff, void* user);

/*
 * Extract the coefficient structure of a linear system encoded as a single
 * polynomial `equation` that is linear in `unknowns` and polynomial over Q in
 * `vars`. Converts `equation` to one fmpq_mpoly over {vars ∪ unknowns} — FLINT
 * performs the (potentially large) product expansion far faster than a
 * hand-rolled mpq distributor — then invokes `cb` once per resulting monomial
 * term (see flint_lsys_term_fn). This is the fast replacement for a manual
 * distribute-and-accumulate pass over a nested Plus/Times equation tree.
 *
 * All of `vars` and `unknowns` must be distinct symbols. Returns 1 on success
 * (every term delivered to `cb`). Returns 0 WITHOUT invoking `cb` at all when:
 * FLINT is absent; a var/unknown is not a symbol or the two sets overlap;
 * `equation` contains a symbol outside vars∪unknowns, an inexact/irrational
 * coefficient, or a negative/symbolic power (i.e. it is not in Q[vars∪unknowns]);
 * or it is nonlinear in some unknown. The caller then falls back to its own
 * builder. Does not take ownership of any argument.
 */
int flint_linear_system_terms(const Expr* equation,
                              Expr* const* vars, int nvars,
                              Expr* const* unknowns, int nunk,
                              flint_lsys_term_fn cb, void* user);

/* Registers the M1 scaffolding builtin(s). Called from core_init(). */
void flint_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_BRIDGE_H */
