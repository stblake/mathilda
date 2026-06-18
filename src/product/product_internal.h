#ifndef PRODUCT_INTERNAL_H
#define PRODUCT_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "expr.h"

/*
 * Shared helpers for the Product sub-algorithms (Product`Telescoping,
 * Product`Rational, Product`Geometric, Product`QProduct).  Each takes/returns
 * owned Expr* per the usual contract and routes through the evaluator so
 * existing builtins (ReplaceAll, Factor, Together, ...) do the heavy lifting.
 *
 * The direct multiplicative analogue of sum_internal.h: prod_div replaces
 * sum_sub (a / b instead of a - b), and two product-specific helpers
 * (prod_has_symbolic_power, prod_linear_factors) gate the factoring stages.
 */

/* Build head[args...] from `n` owned args, evaluate, free the call. */
Expr* prod_eval(const char* head, Expr** args, size_t n);

/* Evaluate e /. var -> val  (val and e are copied, not consumed). */
Expr* prod_subst(Expr* e, Expr* var, Expr* val);

/* Factor[e]; falls back to a copy of e if Factor leaves itself unevaluated.
 * e is not consumed.  CALLER MUST ENSURE e carries no Power[base, var] (a
 * symbolic exponent makes Factor/Together loop -- see prod_has_symbolic_power). */
Expr* prod_factor(Expr* e);

/* True if e is free of var (FreeQ[e, var]). */
bool prod_free_of(Expr* e, Expr* var);

/* Convenience: integer node. */
Expr* prod_int(int64_t v);

/* a / b  =  Times[a, Power[b, -1]], evaluated  (both copied). */
Expr* prod_div(Expr* a, Expr* b);

/* Parse the (f, i[, imin, imax]) argument shape shared by every stage.
 * Returns true and sets f and var (and, when definite, imin and imax) for
 * head[f,i] (indefinite, argc==2) or head[f,i,imin,imax] (definite, argc==4).
 * The out pointers alias res's args (not copied). */
bool product_stage_args(Expr* res, Expr** f, Expr** var,
                        Expr** imin, Expr** imax, bool* definite);

/* VerifyConvergence option for infinite products (default 1/True): when 0 the
 * Product`Infinite convergence gate is skipped.  Set by the dispatcher around
 * the cascade; read by Product`Infinite. */
extern int g_product_verify_convergence;

/* True if e contains Power[base, p] where p involves var (e.g. a^i, q^k).
 * Such a factor makes Together/Factor loop, so the factoring stages
 * (Telescoping/Rational) bail and leave it to Product`Geometric/QProduct. */
bool prod_has_symbolic_power(Expr* e, Expr* var);

/* Factor a polynomial in `var` into a leading constant and a list of
 * (root, multiplicity) pairs over Q.  On success returns true and fills:
 *   *lead_out   -- owned leading-coefficient Expr (c)
 *   *roots_out  -- owned array of `n_out` owned root Expr* (r_i, so the factor
 *                  is (var - r_i)); caller frees each and the array
 *   *mults_out  -- malloc'd array of `n_out` int multiplicities; caller frees
 *   *all_linear_out -- true iff every irreducible factor was linear over Q
 *                      (false signals an irreducible quadratic+ was present;
 *                      such factors are NOT represented in roots_out)
 * Returns false if e is not polynomial in var (nothing allocated).
 * Built on Factor + get_degree_poly/get_coeff; e must be free of Power[*, var]. */
bool prod_linear_factors(Expr* e, Expr* var,
                         Expr** lead_out,
                         Expr*** roots_out, int** mults_out, size_t* n_out,
                         bool* all_linear_out);

#endif
