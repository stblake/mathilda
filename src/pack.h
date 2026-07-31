#ifndef MATHILDA_PACK_H
#define MATHILDA_PACK_H

/* ---------------------------------------------------------------------------
 * Automatic packed arrays.
 *
 * A PACKED LIST is an ordinary List stored as a dense machine-precision buffer
 * (an EXPR_NDARRAY with present_as == NDA_HEAD_LIST) instead of one Expr node
 * per element. It is invisible: Head is List, it prints as {1., 2., 3.},
 * ListQ/VectorQ/MatrixQ are True, AtomQ is False, patterns match, and
 * expr_eq/expr_hash/expr_compare agree with the plain List element for element.
 * Only NDArrayQ[] tells the two apart.
 *
 * WHY. On a 10^6-element List of Reals, `Length[x]` costs ~30 ms and `x[[7]]`
 * ~32 ms -- not because either is doing work, but because evaluate_step sweeps
 * every argument on every pass (src/eval.c, the per-argument loop), mallocs two
 * 8 MB args arrays, and rebuilds the node. That is ~30 ns per element per
 * evaluation pass, charged every time the list is an argument to anything.
 * Packed, the same operations are 2 us and 82 us. Packing is not an
 * optimisation of list operations; it is the removal of a standing tax on
 * holding a list at all.
 *
 * THE CONTRACT. Packing changes representation and NOTHING else. Not a value,
 * not an element's head, not a printed form, not an ordering. `Range[300]`
 * packed must still answer `Integer` to `Head[Range[300][[1]]]` and
 * `500000500000` (an Integer, promoting past int64 when it has to) to
 * `Total[Range[10^6]]`. Where a packed representation cannot express the exact
 * answer -- an overflow, a Rational, a symbolic result -- the operation
 * abandons and the interpreter re-runs it on the materialised List, which is
 * always right. Slower on that path, never wrong. This is the same rule the
 * Compile[] engine lives under (src/compile/compile_internal.h) and the same
 * one src/numloop.h states for exactness.
 *
 * WHERE THE PIECES LIVE.
 *   expr.h            NDPresentation, and expr_new_ndarray_like vs _raw
 *   expr.c / sort.c   the identity trio: expr_eq, expr_hash, expr_compare
 *   print.c           streaming {...} rendering
 *   eval.c            the transparency gate (a packed list never reaches a head
 *                     that has not opted in)
 *   this file         the packing decision, the producer API, and the builtins
 *
 * See docs/design/packed_arrays.md.
 * ------------------------------------------------------------------------- */

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Minimum TOTAL element count (product of dims, not the leading dimension) for
 * automatic packing. Matches Mathematica's $MinPackedArrayLength default.
 *
 * The performance break-even in this codebase is around n = 2: packing costs
 * ~26 ns/element and saves ~30 ns/element on the very next evaluation pass. So
 * 250 is chosen for BLAST RADIUS, not for cost. Below it lives essentially all
 * symbolic and pattern-matching code, where a 300 ns saving is not worth
 * routing a value through a second representation. Counting total elements
 * rather than length is deliberate: a 3 x 100000 matrix is the single
 * best-value shape in the space and a length test would refuse it. */
#define PACK_MIN_ELEMENTS ((size_t)250)

/* Master switch. Off when MATHILDA_NO_PACK is set in the environment (read
 * once, lazily), or when pack_set_enabled(false) has been called. The
 * differential test drives the C setter to prove that packing on and packing
 * off give bit-identical answers across the corpus. */
bool pack_enabled(void);
void pack_set_enabled(bool on);

/* Set the first time a packed list is built, and never cleared -- a packed list
 * can still be reachable from a symbol long after packing is switched off.
 *
 * The transparency gate in evaluate_step runs on every function node in the
 * system, so it is guarded by this rather than by pack_enabled(): a session
 * that never packs anything pays one predictable load per node instead of an
 * arg scan. Read through pack_any_created(), which inlines to that load.
 * Written only by pack.c. */
extern bool pack_g_any_created;
static inline bool pack_any_created(void) { return pack_g_any_created; }

/* Test-only threshold override; PACK_MIN_ELEMENTS is restored by passing 0.
 * Lets the differential suite exercise packing on small, inspectable values,
 * and check that the threshold itself is not observable (n = 249/250/251 must
 * agree after length normalisation). */
void   pack_set_min_elements(size_t n);
size_t pack_min_elements(void);

/* ---------------------------------------------------------------------------
 * The general path: offer an already-built List for packing.
 *
 * Consumes `list` (owned) and returns either a packed list -- freeing the
 * input -- or `list` unchanged. Never changes a value, so a caller can wrap any
 * List return in it without reasoning about the contents.
 *
 * Declining is cheap: the shape/type probe allocates nothing and bails at the
 * FIRST element that cannot be packed, so a list of symbols costs O(1) and a
 * list with one bad element at the end costs one pointer walk.
 *
 * Deliberately NOT hooked into expr_new_function: that is a constructor called
 * with NULL slots for callers that fill them afterwards, it runs again on every
 * evaluation pass (so a list that fails to pack would pay the probe forever),
 * and it sits on every function node in the system including Plus[a,b]. Like
 * autocompile_new, this is called by the producers that want it.
 * ------------------------------------------------------------------------- */
Expr* pack_offer(Expr* list);

/* pack_offer ignoring the threshold and the master switch: what ToNDArray[]
 * does. Returns `list` unchanged if it is not packable at all. */
Expr* pack_force(Expr* list, bool have_dtype, NDType dtype);

/* Re-pack an operation's List output back into a buffer, for an operation whose
 * INPUT was a packed list (`src`). Consumes `list`; returns it unchanged when it
 * cannot pack.
 *
 * The dtype is SNIFFED from `list`, not taken from `src`, and that is the whole
 * point. Repacking at the source's dtype silently coerced any element the
 * operation introduced: Join[Range[1., 300.], {1}] came back with 1. in the last
 * slot instead of the exact Integer 1, and Riffle[Range[1., 300.], 1] turned
 * every interleaved 1 into 1.. A mixed exact/inexact result now declines and
 * stays an ordinary List, which is the correct value.
 *
 * Only for a PACKED source. A visible NDArray[...] keeps its documented
 * machine-buffer semantics, where coercing an appended exact value is what the
 * user asked for by naming the head. */
Expr* pack_repack_like(const Expr* src, Expr* list);

/* Materialise a packed list (or a visible NDArray) to a nested List. Returns
 * NULL if `e` is not an ndarray at all. Caller owns the result. */
/* Pack every plain numeric List argument of a Listable call up to a buffer, so
 * the head's own NDArray kernel takes the call instead of the buffer being
 * materialised to meet the List. Rewrites `call`'s args in place and sets
 * *changed when it does. All or nothing; see the definition for the four
 * conditions under which it declines and leaves the call to thread.
 *
 * `broadcast_ok` and `int64_ok` are the head's SymbolDef claims. */
void pack_lift_listable_args(Expr* call, bool broadcast_ok, bool int64_ok,
                             bool* changed);

Expr* pack_unpack(const Expr* e);

/* evaluate(), then materialise a packed result. Borrows `e`, exactly like
 * evaluate().
 *
 * THE ONE GAP THE GATE DOES NOT COVER. The gate protects a builtin's ARGUMENTS,
 * because the evaluator puts them there. Nothing protects the value a builtin
 * gets back from its own internal evaluate() call: Median builds Sort[data],
 * evaluates it, and then indexes the result's args -- and once Sort packs, that
 * result is a buffer. The union overlay makes such a site fail loudly (the code
 * tests EXPR_FUNCTION and takes the else branch, so Median[Range[300]] came back
 * UNEVALUATED rather than wrong), but loud is still broken.
 *
 * Use this at any internal evaluate() whose result is then walked structurally.
 * The surface is small and greppable: only a head that can now produce a packed
 * list matters -- Range, Table, ConstantArray, Array, RandomReal, RandomInteger,
 * Sort, Select, NestList, FoldList, Map. */
Expr* pack_eval_plain(Expr* e);

/* ---------------------------------------------------------------------------
 * The direct path: build the buffer, never the Expr nodes.
 *
 * Packing `Range[1., 10^6]` after building 10^6 Expr nodes costs 340 ms + 52 ms;
 * writing 10^6 doubles into a buffer costs under 1 ms. A producer that already
 * knows its shape and dtype must not build the List at all. This is
 * docs/design/compile_state.md's finding in another setting: array speed comes
 * from removing intermediate buffers, not from removing interpretation.
 *
 * Open/fill rather than a per-element callback, because the producers have
 * different loop shapes (Range advances a double, Table runs an AutoCompiled,
 * RandomReal calls an RNG) and a callback would reintroduce one indirect call
 * per element -- exactly the cost being removed.
 *
 * NULL means "not packing" -- disabled, under threshold, or out of memory --
 * and the caller then builds its List exactly as before, so every call site
 * keeps its existing path intact underneath. On success the caller MUST write
 * every element before the node escapes.
 * ------------------------------------------------------------------------- */
Expr* ndbuild_open(int rank, const int64_t* dims, NDType dtype, void** buf_out);

/* Open a buffer for an array DERIVED from `src`, inheriting its presentation and
 * ignoring both the size threshold and pack_enabled().
 *
 * A derived array is not a producer, and the two obey different rules. The
 * interpreter's own array ops derive through expr_new_ndarray_like, so
 * Sin[packedList] is packed at ANY length -- no threshold is re-applied to
 * something that is already a buffer. Applying the producer rule here made
 * Map[#^2 &, ToNDArray[{1., 2., 3., 4.}]] come back as a plain List because four
 * elements are under the threshold, while Sin on the same value stayed packed.
 * pack_enabled() is not consulted for the same reason: the input is already a
 * buffer, so materialising the output would be a downgrade rather than a
 * decline. */
Expr* ndbuild_open_like(const Expr* src, int rank, const int64_t* dims,
                        NDType dtype, void** buf_out);
Expr* ndbuild_open_f64(int64_t n, double** buf_out);
Expr* ndbuild_open_i64(int64_t n, int64_t** buf_out);

/* Abandon a partially-filled buffer: convert the `written` elements already in
 * it to freshly-owned Exprs in `args_out` (caller-allocated, capacity >=
 * `written`) and free the node. Lets a producer start optimistically and still
 * finish on the List path when element k turns out not to be representable. */
void ndbuild_abandon(Expr* node, size_t written, Expr** args_out);

/* Record that the transparency gate is about to materialise `arg` (a packed
 * list) because `head` has not opted in. A no-op unless MATHILDA_PACK_DIAG
 * contains "gate" or "all", in which case an exit summary lists every such head
 * ordered by elements -- a work list of fast paths not taken.
 *
 * Complements tools/check_packed_aware.py, which catches the opposite case: a
 * head that HAS a fast path and is missing from the AWARE list. This one catches
 * a head with NO fast path that real workloads keep feeding buffers to, which no
 * static check can see. */
void pack_gate_report(const Expr* head, const Expr* arg);

/* Register ToNDArray / FromNDArray. Called from core_init. */
void pack_init(void);

#endif /* MATHILDA_PACK_H */
