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
Expr* builtin_subresultants(Expr* res);
void  subresultants_init(void);
Expr* builtin_subresultantpolynomials(Expr* res);
void  subresultantpolynomials_init(void);
Expr* builtin_monomiallist(Expr* res);
Expr* builtin_coefficientrules(Expr* res);
Expr* builtin_coefficientarrays(Expr* res);
Expr* builtin_distributedtermslist(Expr* res);
Expr* builtin_fromcoefficientrules(Expr* res);
void  monomials_init(void);
Expr* builtin_polynomialreduce(Expr* res);
void  polynomialreduce_init(void);

/* ------------------------------------------------------------------ */
/* Internal polynomial helpers used by neighbouring modules           */
/* (expand.c, facpoly.c, parfrac.c, rat.c, simp.c, ...).              */
/* ------------------------------------------------------------------ */

bool is_polynomial(Expr* e, Expr** vars, size_t var_count);
bool is_zero_poly(Expr* e);
int  get_degree_poly(Expr* e, Expr* var);
Expr* get_coeff(Expr* e, Expr* var, int d);
/* Coefficient of var^n in an already-expanded polynomial (var atomic). */
Expr* get_coeff_expanded(Expr* expanded, Expr* var, int n);
/* Bulk coefficient extraction 0..max_deg in one pass; out_coeffs is a   */
/* freshly malloc'd array of (max_deg+1) Expr* the caller frees. Returns */
/* false if var is not atomic.                                           */
bool get_all_coeffs_expanded(Expr* expanded, Expr* var, int max_deg,
                             Expr*** out_coeffs);
/* Full pseudo-remainder lc(B)^(deg A - deg B + 1) * A mod B (deg A>=B). */
Expr* pseudo_rem_standard(Expr* A, Expr* B, Expr* x);
/* Coefficient-wise division of a polynomial in var by a scalar denom.   */
Expr* poly_divide_by_scalar(Expr* poly, Expr* denom, Expr* var);
/* True if e contains an algebraic number Power[base, p/q] with q > 1.   */
bool subres_has_algebraic(Expr* e);
Expr* exact_poly_div(Expr* A, Expr* B, Expr** vars, size_t var_count);
/* exact_poly_div divides by a numeric radical constant (Sqrt[2], 1 + Sqrt[2])
 * in its coefficient base case, except inside a computation whose input carries
 * a nested radical.  enter(e) opens such a scope when e carries one and returns
 * the bump to pass to release(); Together and Cancel bracket their recursion. */
int  exact_poly_div_nested_enter(const Expr* e);
void exact_poly_div_nested_release(int bump);
/* Turn the numeric-radical unit off (and back on) around a caller that wants
 * the classical division: the exact linear-algebra divider. */
void exact_poly_div_radical_units_suspend(void);
void exact_poly_div_radical_units_resume(void);
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

/* ------------------------------------------------------------------ */
/* Opaque-kernel generators.                                          */
/*                                                                    */
/* collect_variables enrols every non-symbol function atom -- Sin[x], */
/* BesselJ[1/3, x], f[x], E^x -- as an independent polynomial         */
/* generator.  Code that then specialises a SYMBOL generator by       */
/* structural substitution (ReplaceAll[P, x -> 0]) also rewrites the  */
/* x inside those kernels, which is not a specialisation of the       */
/* polynomial at all (BesselJ[-2/3, 0] is ComplexInfinity).  Mapping  */
/* each kernel to a fresh symbol first makes the generators the real  */
/* indeterminates the polynomial algorithms assume.                   */
/*                                                                    */
/* poly_kernels_to_symbols: returns a new tree with every maximal     */
/* kernel generator replaced by a fresh symbol (recording the map),   */
/* or NULL when `e` has no kernel generator.                          */
/* poly_symbols_to_kernels: the inverse (structural; the caller       */
/* evaluates).                                                        */
/* ------------------------------------------------------------------ */
typedef struct {
    Expr** kern;    /* owned kernel atoms                 */
    Expr** sym;     /* owned fresh symbols, parallel      */
    size_t n;
} PolyKernelMap;

Expr* poly_kernels_to_symbols(const Expr* e, PolyKernelMap* map);
Expr* poly_symbols_to_kernels(const Expr* e, const PolyKernelMap* map);
void  poly_kernel_map_free(PolyKernelMap* map);

/* True when `e` binds or differentiates in a variable (an unevaluated Integrate/D/
 * Sum/Product/Limit/Series/Dt, a Function, an Inactive[op][...] form, a Hold
 * wrapper): an opaque atom to every kernel -> generator substitution, which is
 * unsound inside it (Integrate[1/g, x] == x/g for a generator symbol g). */
bool  poly_is_binding_operator(const Expr* e);

void poly_init(void);

#endif /* POLY_H */
