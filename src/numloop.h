/*
 * numloop.h -- Automatic numeric loop fast-path.
 *
 * A tiny "compile the body to double bytecode and iterate without allocating"
 * fast path for the imperative loop builtins (Nest, Do, For, While) when the
 * loop body is numeric-closed arithmetic over machine-real variables. This is
 * what lets `Nest[3.5 # (1-#)&, 1./Pi, 10^6]` run at libm-loop speed instead of
 * re-interpreting (substitute + evaluate + allocate) the body a million times.
 *
 * Each entry point either returns a freshly-owned Expr* result (fast path
 * taken) or NULL to mean "not applicable -- run the interpreter as usual". The
 * fast path is *transparent*: callers try it first and fall back silently.
 *
 * Correctness contract: the fast path fires only when every referenced variable
 * numericalizes to a machine real AND the result is guaranteed inexact. Its
 * double operations then match the interpreter's existing Real path bit-for-bit
 * (Plus/Times coerce to double; Power uses libm pow). Any operand outside the
 * compilable subset, or a non-finite intermediate (overflow / Log of a negative
 * -> the interpreter would go complex), aborts the attempt and returns NULL so
 * the interpreter produces the authoritative answer.
 *
 * "Guaranteed inexact" is per RESULT POSITION, not per call, because every value
 * these heads produce arrives one of two ways:
 *
 *   COMPUTED by the body -- Real if the body carries a Real literal (which
 *     contaminates the arithmetic whatever the inputs were), or if every input
 *     the body reads is already Real.
 *
 *   PASSED THROUGH from an input, unevaluated -- keeps that input's own type.
 *     The seed is passed through at out[[1]] of every *List head, by Nest at
 *     n = 0, and by NestWhile when the test fails immediately; a list element is
 *     passed through by any body that simply returns it (#2 &). The interpreter
 *     answers `NestList[# + 0. &, 1, 3]` with {1, 1., 1., 1.} -- an exact head
 *     on an inexact tail -- and a uniform buffer of doubles cannot express that.
 *
 * So a body carrying a Real literal is enough to license a computed position but
 * NEVER a passed-through one: those require the input itself to be Real. Where
 * an entry point can expose an input unchanged, it demands a Real there, and the
 * mixed exact/inexact case falls to the interpreter rather than being rounded
 * into uniformity.
 */
#ifndef NUMLOOP_H
#define NUMLOOP_H

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

/* Force the numeric fast-path on/off (default on; also settable via the env var
 * MATHILDA_NO_NUMLOOP). Used by the differential test to compare the compiled
 * and interpreted loop paths bit-for-bit. */
void numloop_set_enabled(bool on);

/* Nest[f, x0, n] scalar form: f a pure Function of Slot[1]. Returns the final
 * iterate boxed as a Real, or NULL. */
Expr* numloop_nest(const Expr* f, const Expr* x0, int64_t n);

/* NestList[f, x0, n]: as numloop_nest, returning {x0, f[x0], ..., f^n[x0]} as a
 * List of n+1 Reals. A dense-array seed falls back -- the scalar path's fused
 * in-place array map produces the final array, and keeping n+1 of them is the
 * interpreter's allocation problem, not loop overhead this can remove. */
Expr* numloop_nestlist(const Expr* f, const Expr* x0, int64_t n);

/* Map[f, expr] at the default level {1}: f a numeric pure function (or numeric
 * head) applied to each machine-number element of `expr`. Returns a same-head
 * expression of Reals, or NULL to fall back (non-numeric element, exact result,
 * uncompilable f). */
Expr* numloop_map(const Expr* f, const Expr* expr);

/* Scan[f, expr]: as numloop_map, discarding each result and returning Null. A
 * numeric-closed body has no side effects, so the loop is observably a no-op --
 * it is still run so that a non-finite element hands the call back to the
 * interpreter, keeping that boundary exactly where the interpreter puts it. */
Expr* numloop_scan(const Expr* f, const Expr* expr);

/* Do[body, {n}] count form. The body is a numeric assignment block: a single
 * Set[var, rhs], or a CompoundExpression of such assignments (multiple state
 * variables, sequential within an iteration). Writes each assigned variable's
 * final value back and returns Null, or NULL to fall back. n must be >= 1. */
Expr* numloop_do_count(const Expr* body, int64_t n);

/* Do[body, {i, imin, imax, di}] integer range form: like numloop_do_count but
 * with `var` an integer loop counter the body may read. Returns Null (the
 * iterator stays localised -- its OwnValue is untouched) or NULL. */
Expr* numloop_do_range(const Expr* body, const Expr* var,
                       int64_t imin, int64_t imax, int64_t di);

/* For[start, test, incr, body] canonical numeric counter loop:
 * start = Set[i, <int>], test = i <cmp> <int>, incr = Increment[i],
 * body = Set[x, rhs] numeric-closed in {x, i}. Writes i and x back and returns
 * Null, or NULL to fall back. */
Expr* numloop_for(const Expr* start, const Expr* test,
                  const Expr* incr, const Expr* body);

/* While[test, body]: test a comparison of numeric-closed operands in the body's
 * variable, body = Set[x, rhs]. Writes x back and returns Null, or NULL. */
Expr* numloop_while(const Expr* test, const Expr* body);

/* Fold[f, x0, list]: f a binary Function (Slot[1]=accumulator, Slot[2]=element)
 * over a list of machine numbers. Returns the final accumulator boxed as a
 * Real, or NULL. */
Expr* numloop_fold(const Expr* f, const Expr* x0, const Expr* list);

/* FoldList[f, x0, list]: as numloop_fold, returning the seed and every partial
 * result as a List of Length[list]+1 Reals. */
Expr* numloop_foldlist(const Expr* f, const Expr* x0, const Expr* list);

/* FixedPoint[f, x0]: iterate f from a machine-real seed until the value stops
 * changing (default SameTest). Returns the fixed point boxed as a Real, or
 * NULL (unsupported form / non-convergence -> interpreter). */
Expr* numloop_fixedpoint(const Expr* f, const Expr* x0);

/* FixedPointList[f, x0]: as numloop_fixedpoint, returning every iterate,
 * INCLUDING the repeated final value -- the list ends ..., fp, fp. */
Expr* numloop_fixedpointlist(const Expr* f, const Expr* x0);

/* NestWhile[f, x0, test]: iterate f while the unary predicate `test` (a pure
 * Function whose body is a comparison of numeric-closed operands in Slot[1])
 * holds. Handles the default m = 1 form. Returns the first value that fails the
 * test boxed as a Real, or NULL. */
Expr* numloop_nestwhile(const Expr* f, const Expr* x0, const Expr* test);

/* NestWhileList[f, x0, test]: as numloop_nestwhile, returning the seed and
 * every iterate up to and including the first value that fails the test. */
Expr* numloop_nestwhilelist(const Expr* f, const Expr* x0, const Expr* test);

/* Accumulate[list]: running sums of a list whose every element is a machine
 * Real, as a same-head expression of Reals. The only entry point with no
 * function argument to compile -- it replaces an evaluate(Plus[..]) per
 * element. Returns NULL if any element is not a Real (a mixed exact/inexact
 * list keeps an exact prefix, which a double buffer cannot express). */
Expr* numloop_accumulate(const Expr* list);

#endif /* NUMLOOP_H */
