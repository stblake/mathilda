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
 * numericalizes to a machine real AND the result is guaranteed inexact (the
 * seed is already a machine Real, or the body carries a Real literal that
 * contaminates the arithmetic). Its double operations then match the
 * interpreter's existing Real path bit-for-bit (Plus/Times coerce to double;
 * Power uses libm pow). Any operand outside the compilable subset, or a
 * non-finite intermediate (overflow / Log of a negative -> the interpreter
 * would go complex), aborts the attempt and returns NULL so the interpreter
 * produces the authoritative answer.
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

/* FixedPoint[f, x0]: iterate f from a machine-real seed until the value stops
 * changing (default SameTest). Returns the fixed point boxed as a Real, or
 * NULL (unsupported form / non-convergence -> interpreter). */
Expr* numloop_fixedpoint(const Expr* f, const Expr* x0);

/* NestWhile[f, x0, test]: iterate f while the unary predicate `test` (a pure
 * Function whose body is a comparison of numeric-closed operands in Slot[1])
 * holds. Handles the default m = 1 form. Returns the first value that fails the
 * test boxed as a Real, or NULL. */
Expr* numloop_nestwhile(const Expr* f, const Expr* x0, const Expr* test);

#endif /* NUMLOOP_H */
