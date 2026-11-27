/* nfrat.h -- Together / Cancel of rational functions whose coefficients are
 * constant algebraic numbers, computed natively over the number field.
 *
 * A rational function in symbols x1..xn whose coefficients are built from
 * constant algebraic atoms -- I, Sqrt[n], n^(1/k), (-1)^(p/q), Root objects,
 * nested radicals -- is an element of K(x1..xn), K = Q(atoms).  The engine:
 *
 *   1. collects the atoms and the variables in one walk, and builds K = Q(theta)
 *      once (the primitive-element images are cached process-wide);
 *   2. converts the expression tree DIRECTLY into a factored rational function
 *      c prod Q_j^f_j / (u prod P_i^e_i) over K[x], with polynomials held in
 *      fmpq_mpoly over (x1..xn, tau) reduced modulo theta's minimal polynomial.
 *      No intermediate ever goes through the evaluator.  Denominator pieces are
 *      kept as written (made monic over K and shared through a piece table), so
 *      a common denominator is an index merge, not a product;
 *   3. cancels every numerator piece against every denominator piece with the
 *      certified modular gcd over K (any number of variables);
 *   4. writes the answer back in the caller's own radicals: the denominator over
 *      a squarefree, gcd-free basis, content in front.  A non-rational constant
 *      is never inverted: a numeric algebraic denominator stays a denominator,
 *      as in the reference (Together[x/(1 + Sqrt[2]) + 1, Extension ->
 *      Automatic] is (1 + x + Sqrt[2] x)/(1 + Sqrt[2])).
 *
 * Cancel differs from Together only at the top: a sum is cancelled term by
 * term, each term reduced to one fraction.
 *
 * Every entry returns the evaluated answer, or NULL to decline (no algebraic
 * atom, no variable, an inexact or AlgebraicNumber coefficient, a head other
 * than Plus / Times / integer Power, a field whose primitive element is not an
 * algebraic integer, an expansion over the size budget, or FLINT compiled out).
 * Inputs are never mutated.
 */
#ifndef NFRAT_H
#define NFRAT_H

#include <stdbool.h>
#include "expr.h"

/* Cheap structural gate: does `e` contain a constant algebraic atom (a Complex,
 * a rational power of a number, a Root object) inside the rational structure?
 * One linear walk; no field arithmetic. */
bool nfrat_has_algebraic(const Expr* e);

/* Does `e` carry Root[f, k]^m with |m| >= deg f in its rational structure?
 * Together / Cancel keep such a power as written (the reference does not
 * reduce it modulo f), so the field engine declines it; the caller reads the
 * Root objects as independent generators instead. */
bool nfrat_has_unreduced_root_power(const Expr* e);

/* A negative power (integer, or a rational power of a number such as
 * 1/Sqrt[2]) inside the rational structure: something to combine. */
bool nfrat_has_denominator(const Expr* e);

/* Number of distinct constant algebraic atoms of `e` (I counted once), stopping
 * at `stop`.  Structural, no field arithmetic. */
int nfrat_count_atoms(const Expr* e, int stop);

/* Together[e] / Cancel[e] over K.  `factored` selects the reference's factored
 * output (user level); 0 gives the expanded numerator / denominator that
 * internal callers read.  `ext` (may be NULL) is an explicit Extension value
 * -- one constant algebraic number or a List of them -- whose generators join
 * the input's, so K = Q(ext, atoms of e). */
Expr* nfrat_together(const Expr* e, int factored, const Expr* ext);
Expr* nfrat_cancel(const Expr* e, int factored, const Expr* ext);

/* Does `e` hold, in its rational structure, a radical Power[B, p/q] whose
 * radicand B has a symbol? */
bool nfrat_has_symbolic_radical(const Expr* e);

/* Together[e] / Cancel[e] (`cancel` != 0) over K for an `e` that also carries
 * radicals of non-constant radicands, each read as a free kernel as the
 * reference does: one generator per radicand as written, B^(p/q) -> g^(p Q/q)
 * (Q the lcm of that radicand's indices), a bare-symbol radicand eliminated
 * (s -> g^Q).  The answer is ONE fraction (Together) over the factored
 * denominator.  NULL to decline: no constant algebraic atom (the generator
 * paths of rat.c handle that), nothing to combine, or nfrat declines. */
Expr* nfrat_radical_kernels(const Expr* e, int factored, const Expr* ext, int cancel);

/* Simplify's canonical form in the radical function field
 *     L = K(x)[h_1..h_r] / (h_i^Q_i - R_i),
 * K the number field of the constant algebraic atoms (Q when there is none),
 * one generator per radicand class (radicands equal up to a positive rational
 * factor; principal-branch exact), denominators rationalised by the norm
 * (Cayley-Hamilton over the traces), the content cancelled over K.  The
 * result is Sum c_m(x) y^m with each c_m in K(x) factored and written in the
 * caller's radicals, a radicand factor absorbed into its radical when that
 * shortens the term.  Integer 0 is a proof that `e` vanishes.  *field (may be
 * NULL) is set to 1 when Kummer theory proves L a field of full degree, so
 * that the form is canonical and a nonzero answer proves `e` nonzero.  NULL
 * to decline: no radical of a non-constant radicand, a nested one, another
 * head, an inexact number, dimension above 64, a zero divisor met while
 * inverting, or the size budget. */
Expr* nfrat_radical_normal_form(const Expr* e, int* field);

#endif /* NFRAT_H */
