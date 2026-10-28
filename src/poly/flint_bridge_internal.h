/*
 * flint_bridge_internal.h
 * -----------------------
 * Internal surface of flint_bridge.c for sibling FLINT engines in src/poly/
 * (flint_trager.c).  Not a public API: it exposes, under an `fb_` prefix, the
 * Expr <-> fmpq_mpoly conversions and the number-field machinery that
 * flint_field_gcd is built from -- the primitive-element normalisation of a set
 * of algebraic atoms, the tau lift, the certified modular field GCD and the
 * render-back into the caller's own radicals -- so a new engine reuses them
 * instead of duplicating them.
 *
 * Only meaningful when USE_FLINT is defined; the whole header is empty
 * otherwise.
 */
#ifndef FLINT_BRIDGE_INTERNAL_H
#define FLINT_BRIDGE_INTERNAL_H

#ifdef USE_FLINT

#include <stddef.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>

#include "expr.h"

/* The name of the field generator's variable in the tau lift: every field
 * element c0 + c1 theta + ... becomes c0 + c1 tau + ... in Q[x.., tau]. */
#define FIELD_TAU_NAME "$MTHFieldExpandTau"

/* The free symbols that become fmpq_mpoly generators, in generator order. */
typedef struct {
    char** names;
    size_t count;
    size_t cap;
} VarSet;

/* The change of basis from theta-coordinates to the product basis of the
 * caller's own atoms (alpha_1^e1 ... alpha_k^ek, e_i < deg alpha_i).  Member m
 * of the basis has the exponent digits of m in the mixed radix sel_deg[], the
 * first selected atom varying fastest.  ready == 0 when the atoms do not span. */
#define FG_RENDER_MAX_ATOMS 8
typedef struct {
    int      ready;
    slong    n;              /* [K:Q]                                     */
    slong    nprod;          /* product-basis size (== n when usable)     */
    fmpq_mat_t Pinv;         /* theta-coords -> product-basis coords      */
    Expr**   prod;           /* nprod surface Exprs for the basis members */
    size_t   nsel;           /* atoms selected into the basis             */
    slong    sel_deg[FG_RENDER_MAX_ATOMS];
    const Expr* sel_atom[FG_RENDER_MAX_ATOMS];  /* borrowed from the caller */
} FgRender;

/* Cached ORD_LEX context for `nvars` variables (never cleared by the caller). */
const fmpq_mpoly_ctx_struct* fb_ctx_lex(slong nvars);

/* Variable collection (1 on ok; 0 when `e` is not a polynomial over Q in
 * symbols), release, and Expr -> fmpq_mpoly in the order of `vs`. */
int  fb_collect_vars(const Expr* e, VarSet* vs);
void fb_varset_free(VarSet* vs);
int  fb_to_mpoly(const Expr* e, fmpq_mpoly_t out, const fmpq_mpoly_ctx_struct* ctx,
                 const VarSet* vs);

/* Rational-like Expr (Integer, BigInt, Rational) -> fmpq; 1 on ok. */
int  fb_expr_to_fmpq(const Expr* e, fmpq_t out);
Expr* fb_fmpq_to_expr(const fmpq_t q);

/* The distinct non-rational algebraic constants of `e` (appended to *atoms,
 * which the caller frees); 0 when `e` is not a polynomial over such constants. */
int  fb_collect_atoms(const Expr* e, Expr*** atoms, size_t* na, size_t* cap);

/* The common-field images (a List of AlgebraicNumber[theta, {..}], parallel to
 * atoms) -- BORROWED from a process cache; NULL when the field cannot be built. */
const Expr* fb_field_images(Expr* const* atoms, size_t na);

/* Each atom replaced by its image. */
Expr* fb_subst_atoms(const Expr* e, Expr* const* atoms, Expr* const* imgs, size_t na);

/* AlgebraicNumber[theta, {..}] -> polynomial in FIELD_TAU_NAME; a bare
 * Complex[a, b] -> a + b tau when `gaussian` (theta == I), else a decline. */
Expr* fb_field_subst_tau(const Expr* e, int gaussian);

/* M(tau) from the integer coefficient List of theta's minimal polynomial. */
int  fb_field_build_minpoly(fmpq_mpoly_t M, const Expr* mp, const fmpq_mpoly_ctx_struct* ctx,
                            slong nvars, slong tau_idx);

/* A reduced R (tau-degree < n) -> Expr with AlgebraicNumber[theta, ..]
 * coefficients. */
Expr* fb_field_mpoly_to_expr(const fmpq_mpoly_t R, const fmpq_mpoly_ctx_struct* ctx,
                             const VarSet* vs, slong tau_idx, slong n, const Expr* theta);

/* Render-back into the caller's atoms (see FgRender). */
void  fb_render_init(FgRender* R, Expr* const* atoms, Expr* const* imgs, size_t k,
                     const Expr* theta, slong n);
void  fb_render_clear(FgRender* R);
Expr* fb_render_tree(const Expr* e, const FgRender* R);
Expr* fb_expand_algnum(const Expr* e);

/* Certified gcd in Q[gens, tau]/(M) of A and B (both reduced mod M and
 * nonzero; M monic integral, theta an algebraic integer).  G is monic in the
 * LEX-leading gens-monomial.  1 with G set, 0 on decline. */
int  fb_field_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                        const fmpq_mpoly_t M, const fmpq_mpoly_ctx_struct* ctx,
                        const Expr* mp, slong ngens, slong tau_idx, slong n);

/* Does the monic G divide A in Q[gens, tau]/(M)?  When it does and Q is not
 * NULL, Q is set to the quotient reduced mod M. */
int  fb_field_mpoly_divides(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_t G,
                            const fmpq_mpoly_t M, const fmpq_mpoly_ctx_struct* ctx);

/* Kernel classification shared with the Together/Cancel kernel map: an
 * independent transcendental kernel head (Log, the inverse trig/hyperbolic
 * functions), and an inert user-head atom (x[3], Pair[..]) that behaves as an
 * opaque polynomial variable. */
int  fb_km_is_generic_kernel_head(const char* h);
int  fb_km_is_inert_head_atom(const Expr* e);

#endif /* USE_FLINT */
#endif /* FLINT_BRIDGE_INTERNAL_H */
