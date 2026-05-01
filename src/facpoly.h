#ifndef FACPOLY_H
#define FACPOLY_H

#include "expr.h"

Expr* builtin_factorsquarefree(Expr* res);
Expr* builtin_factor(Expr* res);
Expr* builtin_factorterms(Expr* res);
Expr* builtin_factortermslist(Expr* res);
void facpoly_init(void);

/* Per-call Factor memo.
 *
 * `Simplify` runs the candidate-set search across many transforms,
 * each of which may invoke `Factor` on the running candidates.  Many
 * of these calls are on identical (or structurally equal) inputs --
 * e.g., the same trig polynomial after a no-op transform.  A memo
 * keyed on `expr_hash` + `expr_eq` lets us return the cached result
 * for the second, third, ... identical call.
 *
 * Lifecycle: the caller (typically `builtin_simplify`) uses
 * `factor_memo_push` to install a fresh memo at entry, then
 * `factor_memo_pop` to remove it at exit.  Push/pop are nested:
 * the top-of-stack memo is consulted by `builtin_factor`, and on
 * a miss the result is stored before returning.
 *
 * The memo MUST be scoped to a single Simplify call because Factor
 * results depend on the current `$Assumptions`, on user-defined
 * DownValues, and on numeric option settings -- none of which can
 * change inside a single Simplify call but all of which can change
 * between calls.
 *
 * Memory: the memo owns deep copies of every key and value it
 * caches; both are freed when the memo is popped. */
struct FactorMemo;            /* opaque to callers */
typedef struct FactorMemo FactorMemo;

FactorMemo* factor_memo_new(void);
void factor_memo_free(FactorMemo* m);
void factor_memo_push(FactorMemo* m);
void factor_memo_pop(void);

/* Direct memo access for users beyond builtin_factor (e.g.,
 * transform_trig_roundtrip in simp.c uses the same memo to cache
 * its own results, keyed on `TrigRoundtrip[X]` rather than
 * `Factor[X]` -- the key includes the head, so the cache entries
 * never collide).
 *
 * factor_memo_active(): returns the top-of-stack memo, or NULL if
 *   none is active.  Use this to decide whether to bother building
 *   a memo key.
 *
 * factor_memo_lookup(memo, key): returns a borrowed pointer to the
 *   cached value, or NULL on miss.  Caller must NOT free the result.
 *
 * factor_memo_store(memo, key, value): stores deep copies of both.
 *   The caller retains ownership of key and value.
 *
 * Memo lifetime is bounded by the enclosing factor_memo_push /
 * factor_memo_pop pair (which Simplify wraps around its candidate-
 * set search). */
FactorMemo* factor_memo_active(void);
const Expr* factor_memo_lookup(FactorMemo* m, Expr* key);
void factor_memo_store(FactorMemo* m, Expr* key, Expr* value);

#endif // FACPOLY_H
