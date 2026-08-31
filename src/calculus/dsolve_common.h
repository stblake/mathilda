/*
 * dsolve_common.h — shared substrate for the DSolve method cascade.
 *
 * A user call DSolve[eqns, funcs, vars, opts] is parsed once into a
 * DSolveProblem (dsolve_parse); every method then reads that struct, produces
 * candidate general-solution bodies, and the substrate verifies them, fits any
 * initial/boundary conditions, and assembles the Mathematica-shaped result
 * {{u -> Function[{x}, ...]}} / {{u[x] -> ...}}.
 *
 * Scope: Phase-1 ordinary differential equations.  The struct carries the
 * fields the current methods need and is deliberately extensible (systems /
 * PDEs add fields, they do not reshape it).  The small derivative matchers are
 * reimplemented here (rather than shared from the numeric NDSolve module, whose
 * copies are file-static) so the symbolic solver has no link dependency on the
 * numeric one.
 */
#ifndef MATHILDA_DSOLVE_COMMON_H
#define MATHILDA_DSOLVE_COMMON_H

#include "../expr.h"
#include <stddef.h>
#include <stdbool.h>

/* One initial/boundary condition, e.g. y''[0]==5. */
typedef struct {
    size_t fi;      /* index into DSolveProblem.fun_names                     */
    int    order;   /* derivative order (0 = value condition)                 */
    Expr*  point;   /* owned: where imposed (the funcapp argument, e.g. 0)    */
    Expr*  value;   /* owned: required value (e.g. 5)                         */
} DSolveCond;

typedef struct {
    /* dependent functions */
    const char** fun_names;   /* interned; length nfun (borrowed pointers)     */
    size_t       nfun;
    int*         max_order;   /* per-function highest derivative order          */

    /* independent variables */
    const char** ind_names;   /* interned; length nind                         */
    size_t       nind;
    Expr*        xmin;        /* owned or NULL (from {x, xmin, xmax})           */
    Expr*        xmax;        /* owned or NULL                                 */

    /* equations, as residuals R == 0 (R = lhs - rhs), derivatives in
     * Derivative[...] form; conditions are split out into `conds`.            */
    Expr**       eq_residuals; /* owned; length neq                            */
    size_t       neq;

    /* initial / boundary conditions */
    DSolveCond*  conds;       /* length ncond                                  */
    size_t       ncond;

    bool         is_pde;      /* nind >= 2                                     */
    bool         applied;     /* arg2 was u[x] (bare-expr result) vs u (Function) */

    /* parsed options */
    const char*  param_head;  /* GeneratedParameters head, default "C"          */
    Expr*        assumptions; /* owned or NULL                                 */
    const char*  method;      /* pinned Method name, or NULL = Automatic         */
    bool         include_singular;
} DSolveProblem;

/* Parse a DSolve/DSolve`Method call `res` into `P`.  Returns false (leaving the
 * call to stay symbolic) on a shape it does not understand.  On true the caller
 * must dsolve_problem_free(P). */
bool  dsolve_parse(Expr* res, DSolveProblem* P);
void  dsolve_problem_free(DSolveProblem* P);

/* A method's core routine: given the parsed problem, return a freshly malloc'd
 * array of *nbranch general-solution bodies (each an expression in the
 * independent-variable symbol and the generated constants C[k]) for the single
 * dependent function, or NULL to decline.  Ownership of the array and its
 * elements passes to the caller. */
typedef Expr** (*DSolveTryFn)(DSolveProblem* P, size_t* nbranch);

/* Run one method end-to-end: call `fn`, drop any branch whose back-substitution
 * is decidably non-zero, fit initial/boundary conditions, and assemble the
 * result.  Returns the {{...}} rule-list, or NULL if the method declined or no
 * branch survived verification. */
Expr* dsolve_run(DSolveProblem* P, DSolveTryFn fn);

/* ---- systems (nfun > 1): a solution is one body per dependent function ---- */
/* A system method: given the parsed problem, produce `bodies[i]` for each
 * dependent function (in x and the generated constants), or NULL to decline. */
typedef Expr** (*DSolveSysFn)(DSolveProblem* P);

/* Assemble {{y1 -> Function[{x},b1], ..., yn -> Function[{x},bn]}} (or the
 * applied u[x] -> b form), renaming C[k] to the GeneratedParameters head.
 * `bodies` is borrowed (length nfun). */
Expr* dsolve_assemble_system(const DSolveProblem* P, Expr** bodies);

/* Fit the system's initial/boundary conditions in place: collect the constants
 * across all bodies, solve the conditions for them, and back-substitute into
 * every body.  No-op when there are no conditions. */
void  dsolve_fit_system(const DSolveProblem* P, Expr** bodies);

/* True unless some equation, under y_j -> Function[{x}, bodies[j]] for all j,
 * back-substitutes to a decidably non-zero residual. */
bool  dsolve_verify_system(const DSolveProblem* P, Expr** bodies);

/* Run one system method end-to-end: call `fn` for the per-function bodies,
 * verify, fit conditions, assemble.  NULL if declined or verification fails. */
Expr* dsolve_run_system(DSolveProblem* P, DSolveSysFn fn);

/* Renumber the generated constants C[1..m] in `body` to C[*offset+1..*offset+m]
 * and advance *offset by m, so per-function solutions in a system do not
 * collide.  `body` is consumed; result owned.  (m <= 0 is a no-op.) */
Expr* dsolve_renumber_constants(Expr* body, int m, int* offset);

/* Extract the solution body for the interned name `fname` from a scalar DSolve
 * result {{fname -> Function[{x}, body]}}.  NULL if absent; `r` borrowed,
 * result owned. */
Expr* dsolve_extract_system_body(Expr* r, const char* fname);

/* Run one PDE method (single function of the nind independent variables): call
 * `fn` for the body (bodies[0]), verify it against the equation, and assemble
 * {{u -> Function[{v1,...,vk}, body]}} (or the applied u[...] -> body form). */
Expr* dsolve_run_pde(DSolveProblem* P, DSolveSysFn fn);

/* Shared REPL entry for a `DSolve`Method[...]` builtin: parse `res`, run `fn`
 * once (strict, no cascade), free the problem.  Returns the result or NULL. */
Expr* dsolve_method_builtin(Expr* res, DSolveTryFn fn);

/* For a single ODE, return F such that (the order-n derivative) == F, where F
 * is an expression in the independent-variable symbol and the dependent
 * funcapps of order < n.  NULL if the equation is not (uniquely) solvable for
 * the top derivative (e.g. it enters nonlinearly).  Caller owns the result. */
Expr* dsolve_solve_top_derivative(DSolveProblem* P, int n);

/* The single first-order ODE's residual with y[x] -> the plain symbol `Yname`
 * and y'[x] -> the plain symbol `Pname` (both interned).  Lets methods that need
 * the form NOT solved for y' (Exact, Clairaut) work algebraically.  Owned. */
Expr* dsolve_algebraic_residual(DSolveProblem* P, const char* Yname, const char* Pname);

/* Solve the linear first-order ODE  v' + Pcoef(x) v == Qcoef(x)  by the
 * integrating factor mu = Exp[Integrate[Pcoef, x]], returning the general body
 * (Integrate[mu Qcoef, x] + C[1]) / mu, or NULL if an integral is not
 * elementary.  Pcoef and Qcoef are consumed. */
Expr* dsolve_linear_factor_solve(Expr* Pcoef, Expr* Qcoef, const char* xvar);

/* Extract the RHS values of the single-variable solutions in a Solve result
 * List[List[Rule[var, val]], ...] whose rule LHS is the symbol `varname`.
 * Returns a malloc'd array of *n owned values, or NULL if none.  Borrows solres. */
Expr** dsolve_extract_solutions(Expr* solres, const char* varname, size_t* n);

/* Distinct roots of a polynomial with per-root multiplicity + Im part, used to
 * build the fundamental set of a linear ODE (constant- or Euler-coefficient).
 * `roots`/`im` are owned arrays of length `ndist`; `total` is the sum of the
 * multiplicities (== degree when the polynomial fully factors). */
typedef struct {
    Expr** roots;   /* distinct roots (owned)                    */
    int*   mult;    /* multiplicity of each                       */
    Expr** im;      /* Im[root] (owned)                          */
    bool*  isreal;  /* Im[root] == 0                             */
    size_t ndist;
    int    total;   /* sum of mult                               */
} DSolveRoots;

/* Solve poly == 0 for `var`, dedupe, and fill multiplicities (derivative test)
 * + Im parts.  Returns false (leaving *out zeroed) if no roots are found. */
bool  dsolve_analyze_roots(const Expr* poly, const char* var, int degree, DSolveRoots* out);
void  dsolve_roots_free(DSolveRoots* r);

/* For a single ODE linear in y and its derivatives, extract the coefficient
 * c_k(x) of y^(k)[x] (k = 0..order) and the forcing g(x), so the equation reads
 * Σ c_k y^(k) == g.  Returns false when the equation is not linear.  On true,
 * *coeffs is a malloc'd array of order+1 owned Exprs, *forcing is owned, and
 * *order is the highest derivative order. */
bool  dsolve_linear_coeffs(DSolveProblem* P, Expr*** coeffs, Expr** forcing, int* order);

/* Particular solution of a linear ODE by variation of parameters over the
 * fundamental set `basis` (length n), forcing `g`, leading coefficient
 * `leadcoef` (the coefficient of y^(n); constant or x-dependent).  Returns the
 * particular solution (Simplify-reduced), or NULL if an integral is not
 * elementary.  `basis`, `g`, `leadcoef` are borrowed. */
Expr* dsolve_variation_of_parameters(Expr** basis, size_t n, const Expr* g,
                                     const Expr* leadcoef, const char* xvar);

/* ---- shared building blocks used by the methods ---- */

/* head[a] / head[a, b] — head is an interned symbol pointer; args consumed. */
Expr* ds_call1(const char* head, Expr* a);
Expr* ds_call2(const char* head, Expr* a, Expr* b);

/* D[e, v] / Integrate[e, v] / Solve[eq, v] from C; args consumed, result owned. */
Expr* ds_d(Expr* e, Expr* v);
Expr* ds_integrate(Expr* e, Expr* v);
Expr* ds_solve(Expr* eq, Expr* v);

/* ReplaceAll[body, from -> to], evaluated; body/from/to consumed. */
Expr* ds_subst(Expr* body, Expr* from, Expr* to);

/* Build the literal y[x] (order 0) or Derivative[m][y][x] (m>=1). */
Expr* ds_make_funcapp(const char* fname, int order, const char* xvar);

/* The generated constant C[k]. */
Expr* ds_const(int k);

/* True iff the interned symbol `name` occurs anywhere in `e` (head or arg). */
bool  ds_contains(const Expr* e, const char* name);

/* True iff `e` does not depend on the variable `var` — tested robustly as
 * D[e, var] == 0, so it holds even when `e` is unsimplified (e.g. a rational
 * form that only cancels the variable after Together).  `var` is a symbol name. */
bool  ds_free_of(const Expr* e, const char* var);

/* Simplify[e]; `e` consumed, result owned.  Used to reduce a value that is known
 * to be constant (e.g. a recovered exponent) to its closed form. */
Expr* ds_simplify(Expr* e);

/* Materialize a packed NDArray list to a plain List (recursively over its
 * direct elements); consumes e, returns owned.  No-op on a plain expression. */
Expr* ds_delist(Expr* e);

/* Zero recognition: strict (proved zero) / proved-nonzero, via zero_test. */
bool  ds_is_zero(const Expr* e);      /* ZERO_TEST_TRUE                        */
bool  ds_is_nonzero(const Expr* e);   /* ZERO_TEST_FALSE                       */

/* True iff `e` still contains the head `Integrate` / `Solve` (i.e. a call the
 * subsolver could not carry out) — used by methods that require a closed form. */
bool  ds_has_head(const Expr* e, const char* head);

#endif /* MATHILDA_DSOLVE_COMMON_H */
