#ifndef MATHILDA_SIMP_H
#define MATHILDA_SIMP_H

#include "expr.h"

/*
 * simp.c -- Simplify, Assuming, and the assumption-context infrastructure.
 *
 * Simplify implements a heuristic search over a battery of algebraic
 * transforms (Together, Cancel, Expand, Factor, FactorSquareFree, Apart,
 * TrigExpand, TrigFactor, TrigToExp/ExpToTrig round-trip, ...) and returns
 * the candidate of minimum complexity. The default complexity measure
 * matches Mathematica's: total leaf count plus the decimal-digit count of
 * any integer leaves, so e.g. 100 Log[2] is preferred over its expanded
 * Log[2^100] form. A user-supplied ComplexityFunction option overrides the
 * default.
 *
 * Assumption-aware rewrite passes (run as additional seed candidates) use
 * AssumeCtx, a normalized fact set parsed from the user's assumption
 * expression. AssumeCtx is exposed here so future modules (Refine,
 * Integrate, ...) can share it.
 *
 * Assuming[a, body] is implemented by desugaring to
 *   Block[{$Assumptions = $Assumptions && a}, body]
 * and reusing Block's existing scoping machinery; nested Assuming calls
 * therefore compose naturally.
 */

/* Module init -- registers Simplify, Assuming, $Assumptions, Element.
 * Called from core_init(). */
void simp_init(void);

/* Default complexity measure: LeafCount(e) + sum of decimal-digit counts
 * of every integer/bigint leaf. Used by Simplify when no
 * ComplexityFunction is supplied. */
size_t simp_default_complexity(const Expr* e);

/* ----------------------------------------------------------------------- */
/* Assumption context                                                      */
/* ----------------------------------------------------------------------- */

/*
 * A flat list of normalized facts derived from a user assumption
 * expression. And[...] is flattened, list arguments are converted to
 * conjunctions, and uninformative True is dropped. If the assumption set
 * is internally inconsistent (e.g. contains literal False), `inconsistent`
 * is set to true and downstream callers should bail out.
 */
typedef struct AssumeCtx {
    Expr** facts;        /* each fact owned by the ctx; deep-copied */
    size_t count;
    size_t capacity;
    bool   inconsistent;
} AssumeCtx;

/* Build a context from an assumption expression. The expression is not
 * mutated. Returns a heap-owned ctx (never NULL); free with
 * assume_ctx_free. */
AssumeCtx* assume_ctx_from_expr(const Expr* assum);

/* Free a context allocated by assume_ctx_from_expr. */
void assume_ctx_free(AssumeCtx* ctx);

/* Domain queries. Each returns true only if the assumption set provably
 * entails the property under v1's narrow reasoning rules (direct facts
 * plus the obvious lattice Integer => Rational => Real => Complex,
 * x > 0 => x >= 0 => x in Reals, x in Reals && Equal[x,0]==False handled
 * trivially). No transitive chaining across multiple inequality facts. */
bool assume_known_positive(const AssumeCtx* ctx, const Expr* x);
bool assume_known_nonneg  (const AssumeCtx* ctx, const Expr* x);
bool assume_known_negative(const AssumeCtx* ctx, const Expr* x);
bool assume_known_nonpos  (const AssumeCtx* ctx, const Expr* x);
bool assume_known_real    (const AssumeCtx* ctx, const Expr* x);
bool assume_known_integer (const AssumeCtx* ctx, const Expr* x);
bool assume_known_even    (const AssumeCtx* ctx, const Expr* x);

#endif /* MATHILDA_SIMP_H */
