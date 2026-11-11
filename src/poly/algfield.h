/* algfield.h -- Together / Cancel over the number field of the constant
 * algebraic coefficients.
 *
 * A rational function whose coefficients are constant algebraic numbers built
 * from NESTED radicals (Sqrt[3/4 - 1/Sqrt[2]], Sqrt[1 + Sqrt[2]], (1 + I Sqrt[3])^(1/3))
 * defeats the generator model the rest of Together uses, where every radical is
 * an indeterminate tied to its radicand.  A nested radical's relation has
 * algebraic coefficients, it may even be reducible (Sqrt[3/4 - 1/Sqrt[2]] is
 * 1/Sqrt[2] - 1/2), and the multivariate gcds over the atoms swell without
 * bound.  This module gives the two answers that class needs:
 *
 *   - algfield_normalize: the exact computation over K = Q(atoms).  Every atom
 *     becomes an AlgebraicNumber over one primitive element (FLINT qqbar), the
 *     combination runs natively over K, and each coefficient of the answer is
 *     written back over a basis of radical monomials in the atoms themselves,
 *     simplest atoms first, so an atom that lies in the field of simpler ones
 *     (a denestable radical) drops out.  This backs Extension -> Automatic.
 *
 *   - algfield_kernel_model: the default (no Extension) reading, run by
 *     Together / Cancel only after their fast paths decline.  Each nested
 *     radical's base first takes its Together form (a radicand that combines to
 *     a radical monomial is not nested after all, and goes back to the ordinary
 *     cascade).  Each remaining nested radical b^(p/q) is an opaque generator
 *     g = b^(1/m), the combination runs with g a fresh symbol, and numerator and
 *     denominator are reduced modulo g^m - b.  No denesting, no
 *     rationalisation, as in the reference.
 *
 * Both are gated by algfield_has_nested_constant_radical, a single linear walk,
 * so ordinary inputs (no nested constant radical) pay nothing else.
 */
#ifndef ALGFIELD_H
#define ALGFIELD_H

#include <stdbool.h>
#include "expr.h"

/* True when `e`, read through Plus / Times / Power (the ring structure
 * Together combines through), contains a radical b^(p/q) (q >= 2) whose base b
 * is a constant algebraic number that itself carries a radical or a Root
 * (a Gaussian base such as I^(1/3) does not count).  Function arguments are
 * opaque and are not searched. */
bool algfield_has_nested_constant_radical(const Expr* e);

/* `head` is SYM_Together or SYM_Cancel.  Reduces `e` over the number field of
 * its constant algebraic atoms (see above) and returns the evaluated result,
 * or NULL to decline: no algebraic atom, an inexact or AlgebraicNumber
 * coefficient, a field of degree above the cap, or FLINT compiled out.  Never
 * mutates `e`. */
Expr* algfield_normalize(const Expr* e, const char* head);

/* Default-mode reading of a nested constant radical: an opaque generator.
 * `call` is the Together[...] / Cancel[...] call (head `head`); `compute` is
 * the builtin's dispatch, called directly on a rewritten call so the answer
 * keeps the caller's user-level form.  Returns the evaluated result or NULL to
 * decline.  Never mutates `call`. */
Expr* algfield_kernel_model(const Expr* call, const char* head,
                            Expr* (*compute)(Expr*));

#endif /* ALGFIELD_H */
