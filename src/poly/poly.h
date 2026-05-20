#ifndef POLY_H
#define POLY_H

#include "expr.h"

/* ------------------------------------------------------------------ */
/* Public built-in entry points (called by the evaluator).            */
/* Each follows the Mathilda convention: caller owns `res` (input) and */
/* the builtin returns either NULL (could not evaluate -- res stays   */
/* unchanged) or a freshly allocated Expr* that the caller must free. */
/* ------------------------------------------------------------------ */

Expr* builtin_polynomialq(Expr* res);
Expr* builtin_variables(Expr* res);
Expr* builtin_coefficient(Expr* res);
Expr* builtin_coefficientlist(Expr* res);
Expr* builtin_polynomialgcd(Expr* res);
Expr* builtin_polynomiallcm(Expr* res);
Expr* builtin_polynomialquotient(Expr* res);
Expr* builtin_polynomialremainder(Expr* res);
Expr* builtin_polynomialquotientremainder(Expr* res);
Expr* builtin_subresultantpolynomialremainders(Expr* res);
Expr* builtin_polynomialmod(Expr* res);
Expr* builtin_polynomialextendedgcd(Expr* res);
Expr* builtin_collect(Expr* res);
Expr* builtin_decompose(Expr* res);
Expr* builtin_hornerform(Expr* res);
Expr* builtin_resultant(Expr* res);
Expr* builtin_discriminant(Expr* res);

/* ------------------------------------------------------------------ */
/* Internal polynomial helpers used by neighbouring modules           */
/* (expand.c, facpoly.c, parfrac.c, rat.c, simp.c, ...).              */
/* ------------------------------------------------------------------ */

bool is_polynomial(Expr* e, Expr** vars, size_t var_count);
bool is_zero_poly(Expr* e);
int  get_degree_poly(Expr* e, Expr* var);
Expr* get_coeff(Expr* e, Expr* var, int d);
Expr* exact_poly_div(Expr* A, Expr* B, Expr** vars, size_t var_count);
Expr* poly_gcd_internal(Expr* A, Expr* B, Expr** vars, size_t var_count);
Expr* poly_content(Expr* A, Expr** vars, size_t var_count);
void  collect_variables(Expr* e, Expr*** vars_ptr, size_t* count, size_t* capacity);
int   compare_expr_ptrs(const void* a, const void* b);
bool  contains_any_symbol_from(Expr* expr, Expr* var);

/* ------------------------------------------------------------------ */
/* Algebraic-generator substitution helpers.                          */
/*                                                                    */
/* Detect a (base B, atom A) pair such that the input is a polynomial */
/* / rational function in Power[B, A] (after fixing a common rational */
/* scaling of exponents), substitute Power[B, A] -> g (g a fresh      */
/* symbol), let the caller run the polynomial operation in g, then    */
/* substitute back.  Two flavours fall under this scheme:             */
/*                                                                    */
/*   * Radical case (A = 1): exponents are pure rationals p/q.        */
/*     We pick m = lcm of q's, substitute B -> g^m so that            */
/*     B^(p/q) -> g^(p*m/q) is an integer power.  Triggers when       */
/*     at least one fractional exponent exists (m > 1).               */
/*                                                                    */
/*   * Exponential case (A != 1): exponents are rational multiples    */
/*     of a common atom A, e.g. Power[E, 2x] and Power[E, x] share    */
/*     atom A = x.  We pick m = lcm of the rational scalings'         */
/*     denominators (commonly 1), substitute g = Power[B, A/m] so     */
/*     that Power[B, c*A] -> g^(c*m) is an integer power.  Triggers   */
/*     when at least two matching Power sites exist (i.e. the         */
/*     polynomial in g is non-trivial).                               */
/*                                                                    */
/* Both B and A may be arbitrary expressions; structural identity     */
/* across occurrences is checked via expr_eq.                         */
/* ------------------------------------------------------------------ */

bool  poly_find_radical_gen(Expr* e, Expr** base_out, Expr** atom_out, int64_t* m_out);
Expr* poly_subst_radical_to_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen);
Expr* poly_subst_radical_from_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen);
char* poly_make_fresh_gen(Expr* e);

void poly_init(void);

#endif /* POLY_H */
