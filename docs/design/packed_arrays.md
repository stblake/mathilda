# Packed arrays — design

How an ordinary `List` comes to be stored as a dense buffer without anything
being able to tell, and what that costs. User-facing behaviour is in
[`../spec/builtins/packed-arrays.md`](../spec/builtins/packed-arrays.md).

Status: automatic packing is on and it crosses the `Compile[]` boundary (§7). §6
lists what remains — complex packing and the remaining int64 fast paths, both
optimisations rather than correctness items. §8 records the six wrong answers
that switching it on exposed, because each says something about where the
abstraction leaks.

---

## 1. The problem

Mathilda already had a dense machine-precision array (`EXPR_NDARRAY`) with ~83
fast-path dispatch sites, ~80 element-wise kernels, BLAS/LAPACK bridges and a
degrade path. None of it helped ordinary list code, because reaching it meant
typing `NDArray[...]` and accepting a visibly different value.

Measured on a 10^6-element list of Reals:

| | ordinary `List` | packed |
|---|---|---|
| `Length[x]` | 30.6 ms | 2 µs |
| `x[[7]]` | 32.3 ms | 82 µs |
| `Total[x]` | 95.8 ms | 4.0 ms |
| `Sin[x]` | 320 ms | 16.8 ms |

`Length` is the instructive row: 30 ms to answer a question the node header
already knows. The cost is `evaluate_step`'s per-argument sweep — 10^6
`evaluate()` calls, two 8 MB `args` allocations, a rebuilt node — about **30 ns
per element per evaluation pass, charged whenever the list is an argument to
anything**. Packing does not make list operations faster so much as remove a
standing tax on holding a list at all.

`docs/design/compile_state.md` §0d records the same finding from the other side:
a compiled `Map` reads 1.0x over a plain `List` and 277x packed-in/packed-out,
the gap being entirely boundary marshalling.

---

## 2. Decisions

### D1 — one node type, two presentations

`NDArrayData` carries `NDPresentation present_as`:

- `NDA_HEAD_NDARRAY` (0) — the explicit `NDArray[...]`. `Head` is `NDArray`,
  prints wrapped, `AtomQ` is `True`, and a `Part` assignment coerces the
  right-hand side into the buffer: the user asked for a machine buffer.
- `NDA_HEAD_LIST` — a packed list. Indistinguishable from the nested `List`.

Zero is the visible surface, so every construction site predating this keeps its
meaning. A separate `EXPR_PACKEDLIST` type was considered and rejected for cut 1
— it would make every `switch (e->type)` a compile error, which is safer, but at
the cost of touching the 83 `is_ndarray` sites that legitimately want *either*
surface. If the propagation audit (§3) starts producing bugs, that is the
fallback.

`present_as` **is** part of a value's identity, unlike `refcount` or
`last_evaluated_at`: it decides the `Head`, and nothing with different `Head`s is
`SameQ`. So `NDArray[{1.,2.}] === {1.,2.}` is `False` while
`ToNDArray[{1.,2.}] === {1.,2.}` is `True`.

### D2 — a gate in the evaluator, not a retrofit across the tree

`e->data.function.args` is read at ~7100 sites in ~400 files. There is no
accessor layer to hook, and adding one is not a change anyone should attempt.

Instead: **a packed list never reaches a head that has not opted in.**
`evaluate_step` materialises packed arguments for any head without
`SymbolDef.packed_aware`. Every unaware consumer — `Count`, `Level`, `Position`,
`Cases`, `ReplaceAll`, `AtomQ`, `ListQ`, `Insert`, `Append`, user `DownValues`
with structural patterns — is then correct with no code of its own.

This inverts the usual default. Opting a head in is a claim it has been checked;
forgetting to opt one in costs a materialisation, not an answer.

What makes the gate safe rather than merely convenient: `NDArrayData` overlays
`function` in the union so that `arg_count` aliases the `data` heap pointer. A
walker that slipped past the gate and iterated `args[i]` would segfault
essentially immediately, not return a plausible answer. The silent-wrong-answer
class is code that tests `type == EXPR_FUNCTION` and takes the else branch —
exactly what the gate covers.

### D3 — the no-nesting invariant

`List` is deliberately **not** packed-aware. So `{a, packed}` materialises its
packed element as the enclosing `List` is built, and a packed node can never sit
inside a plain `EXPR_FUNCTION` tree.

That is what makes the gate's O(argc) top-level scan *complete*. Without it the
gate would have to walk the whole tree on every unaware head — an unconditional
regression on symbolic workloads, to protect a case that mostly does not arise.

The cost is real and visible: a head whose argument is a list *of* arrays never
sees packed elements. `MapThread[f, {p, q}]` and multi-vector `Transpose` are
the two that matter. Answers are unaffected. Lifting this needs a transitive
"contains a packed node" bit on `Expr` — there is room in the `refcount` word —
plus a measured check that the OR-fold it adds to `expr_new_function`, the
hottest allocator in the system, is affordable.

### D4 — exactness governs everything

An operation over a packed buffer either produces **the exact answer the
interpreter would**, or returns `NULL` so the materialised list path runs
instead. It never yields a `Real` where the interpreter yields an `Integer`, and
never wraps on overflow.

This is `src/compile/autocompile.c`'s doctrine ("the user did not ask for any of
this … would turn a transparent speed-up into a silently wrong answer") and
`src/numloop.h`'s exactness contract, applied to storage instead of to loops.

Concretely: an all-`Integer` list packs to `NDT_INT64` and its elements come back
as exact `Integer`s past 2^53; a mixed `{1, 2.5}` declines, because one dtype
cannot give element 1 an `Integer` head and element 2 a `Real` head.

### D5 — direct construction is where the win is

Packing `Range[1., 10^6]` *after* building 10^6 nodes costs 340 ms + 52 ms;
writing 10^6 doubles costs under 1 ms. Producers that already know their shape
and dtype must never build the list. `ndbuild_open` hands out the raw buffer for
exactly that; open/fill rather than a per-element callback, because a callback
reintroduces one indirect call per element — the cost being removed.

---

## 3. Where the pieces live

| file | role |
|---|---|
| `src/expr.h` | `NDPresentation`; `expr_new_ndarray_raw` vs `expr_new_ndarray_like` |
| `src/expr.c` | `expr_eq`, `expr_hash`, `expr_unshare` |
| `src/sort.c` | `expr_compare` |
| `src/print.c` | streaming `{...}` rendering |
| `src/eval.c` | the transparency gate (step 2.7) |
| `src/symtab.{c,h}` | `packed_aware`, `symtab_set_packed_aware` |
| `src/pack.{c,h}` | the packing decision, producer API, `ToNDArray`/`FromNDArray`, the aware-head list |
| `src/checked_int.h` | overflow-checked integer arithmetic, shared with `Compile[]` |

### The constructor split

`expr_new_ndarray` was renamed `expr_new_ndarray_raw` and joined by
`expr_new_ndarray_like(src, ...)`, which inherits `src`'s presentation. All 48
call sites were then visited — the rename makes that compiler-enforced rather
than grep-enforced, which matters because the failure mode is silent: an array
derived from a packed list with the raw constructor comes back a *visible*
`NDArray`, and transparency breaks at the first operation applied to it.

The same split exists one level up: `ndarray_from_nested_list_like` for the
repack after a materialise-and-reuse detour. `ndstruct_delist_repack` and
`map_try_repack` both go through it — without that, `Map`, `Select`, `Rest`,
`Most` and `Join` all returned visible `NDArray`s from packed input.

### The identity trio

`expr_eq`, `expr_hash` and `expr_compare` sit *below* the evaluator, so the gate
cannot protect them. Each answers "what would `ndarray_to_nested_list` give?"
without calling it — these run on 10^6-element values inside `Association`
lookups and `Union`.

- **`expr_hash`** must be bit-identical to the materialised list's hash, or an
  `Association` keyed by one cannot be found with the other and `Union` splits
  one value in two. FNV-1a composes (a child's hash is an independent run folded
  in with `h ^= child; h *= prime`), which is what makes this computable from
  the buffer.
- **`expr_compare`** puts packed lists in the *List* tier, not the NDArray tier.
  The buffer fast path is valid only when rank and every dim match: element-wise,
  `{2.,0.}` sorts after `{1.,9.,9.}` on its first element but before it on
  length, and list order settles length first. Mismatched shapes materialise.
- **`expr_eq`** compares a packed list against a plain one elementwise, and two
  packed lists without treating dtype as identity — `float32` and `float64` both
  materialise to `Real`s, so they are equal; `int64` and `float64` are not,
  because their heads differ.

`tests/test_packed_list.c` asserts all three against the *actual* materialised
form computed at run time, rather than against hand-written expectations that
could be wrong in the same direction as the code.

---

## 4. The aware-head list

Lives in one reviewable block, `pack_mark_aware_heads()`, rather than as calls
scattered through the modules that own each head. Heads with an element-wise
kernel opt in automatically via the kernel setters.

Deliberately absent, and why:

- `List`, `Rule`, `RuleDelayed`, `Association`, `Hold`, `HoldForm` — these
  *enforce* D3.
- `ListQ`, `VectorQ`, `MatrixQ`, `AtomQ`, `Insert`, `Delete`, `Position`,
  `ReplacePart`, `Append`, `Prepend`, `LeafCount`, `Simplify` — correct by
  omission. Several would be actively wrong if marked aware: `is_listq` is a
  bare head test that reports `False`, `Append` returns unevaluated on an
  ndarray, and `LeafCount` would count a 10^6-element list as one node and skew
  `Simplify`'s complexity metric. Marking these aware is an optimisation with a
  correctness precondition, not a free win.

Two predicates that look like they should have been flipped and were not:

- **`is_atomic` (`src/part.c`)** still reports `true` for a packed list. All 13
  of its call sites go straight on to read `data.function.arg_count`, which
  aliases the payload's `data` pointer — flipping it would turn a heap pointer
  into an element count at 8 of them. Presentation is decided in `expr_head` and
  the `AtomQ` builtin instead.
- **`is_listq` (`src/list/list_common.c`)** likewise: it is called deep inside
  `matrixq.c` immediately before that code walks `args[j]`.

---

## 5. Gate placement

`src/eval.c`, step 2.7 — after the `Unevaluated` strip, before the attribute
block. Forced in both directions:

| must run AFTER | why |
|---|---|
| `flatten_sequences` | `f[Sequence[packed]]` has no packed *argument* until the splice; both passes key on `EXPR_FUNCTION`, so a packed node travels through untouched |
| `Unevaluated` strip | same |

| must run BEFORE | why |
|---|---|
| `Flat` | an unaware associative head must see a real list |
| `Listable` | `has_list_arg` tests for an `EXPR_FUNCTION` headed `List` and returns **false** for a packed node, so an unaware `Listable` head would silently skip threading. Gating first fixes that with no change to `has_list_arg` — while an *aware* head still skips threading, which is what lets its kernel fire |
| `Orderless` | needs `expr_compare`'s packed fork |
| `DownValues` | the matcher cannot descend a buffer |
| the builtin | the point |

Two details that are load-bearing rather than defensive:

- The gate sits **outside** the `head is a symbol` block, so `f[x][y]` is
  covered. A pure `Function`/`&` head is exempt: substitution places the packed
  value into the body where each inner head is gated on the next pass, and
  materialising here would cost `Map[#^2 &, packed]` its fast path.
- Awareness is `packed_aware && !down_values`. `Protected` is opt-in per symbol,
  so even an aware head can carry a user `DownValue` with a structural pattern.

The gate is guarded by `pack_any_created()` — a flag set the first time anything
is packed and never cleared. It runs on every function node in the system, so a
session that never packs pays one predictable load rather than an argument scan.
`bench_eval` shows no measurable change.

---

## 6. Not done yet

1. **The remaining int64 fast paths.** Correctness is settled — every head not on
   the `INT64_OK` list materialises, so the answer is the ordinary list's by
   construction. What is left is speed, and the list is now short:
   - `ndarray_map_unary` / `map_binary` should get exact int64 output buffers for
     the identity-ish heads (`Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign`,
     `Abs`), which would let them back onto `packed_aware` (see §8, defect 5).
     For everything else the right answer is to *bail*: `Sin[{1,2,3}]` is
     symbolic in the interpreter, so computing in float would be wrong, not slow.
   - `ndstruct_delist_repack`'s family — `Join`, `Riffle`, `RotateLeft/Right`,
     `Partition`, `Differences`, `First`/`Last`/`Most`/`Rest` — is exactness-safe
     now that it re-sniffs (§8, defect 3) but has not been re-verified head by
     head against the int64 write paths, so they still materialise.
   - `Chop` must not touch exact values at all.
   `MATHILDA_PACK_DIAG=1` reports every remaining lossy touch of an int64
   buffer, `=abort` aborts at the first so a debugger shows the caller, and
   `nd_int64_lossy_hit` is a stable breakpoint target. Find them by running, not
   by reading.
2. **Complex packing.** Declined for now: a complex element with a zero
   imaginary part materialises as `Complex[re, 0.]`, which the evaluator never
   produces (it folds to a `Real`), so a packed complex list would not
   round-trip. Needs a presentation-aware fold in the element builder.
3. **`bench_pack.c`.** The numbers in the changelog were measured by hand at the
   REPL; they are not yet gated on the `bench_eval` normalized-baseline pattern.
4. **A visible-surface decision, deliberately deferred.**
   `NDArray[{1.,2.,3.}] * {10,20,30}` broadcasts into a 3x3 outer product
   instead of threading elementwise. §8's defect 2 fixed this for packed lists by
   materialising; the visible form was left alone because the right answer there
   is a semantic question (should it come back as a `List` or an `NDArray`?)
   rather than a bug in packing. Likewise `MapIndexed` has never repacked, on
   either surface.

---

## 7. The `Compile[]` boundary

The largest single win in the work, and the one the plan predicted:
`compile_state.md` recorded a compiled `Map` at 1.0x over a plain `List` because
both sides were dominated by marshalling Expr nodes at the boundary. Automatic
packing supplies the packed-in half; this section is the rest.

### The presentation rule

A compiled call must answer with the head the interpreter would give for the same
input. There were two possible inputs and now there are three, so one boolean
("was this argument an NDArray?") had to become a three-way kind:

| argument | result |
|---|---|
| plain `List` | plain `List` — the interpreter threads and returns one |
| packed `List` | packed `List`, **at any size** |
| `NDArray[...]` | `NDArray[...]` |

Joined most-visible-wins across the slots, matching `nd_present_src2` inside the
interpreter's own array ops.

"At any size" is the load-bearing detail, and getting it wrong was the one bug in
this part: a *derived* array inherits presentation with no threshold re-applied
(`Sin[packedList]` is packed however short it is), while a *producer* applies the
threshold. Using the producer rule for a derived result made
`Map[#^2 &, ToNDArray[{1., 2., 3., 4.}]]` come back as a plain List — four
elements being under the threshold — while `Sin` on the same value stayed packed.
Hence `ndbuild_open_like`, which is the derived opener: no threshold, no
`pack_enabled()` consult, presentation inherited.

A body that BUILDS its array has no presentation to inherit. The old code
materialised such a result because "the construct has no packed form"; that stopped
being true when `ConstantArray[1., 300]` started packing, so a built result now
goes through `pack_offer` — the producer rule, threshold and switch included.

Two refusals. A **complex** result never wears the packed presentation, because
`Complex[re, 0.]` does not round-trip (§2 D4's reason, in a new place): the same
call printed `{1. + 0.*I, 4. + 0.*I}` against a plain List's `{1., 4.}`, since the
List path is folded by the evaluator on the way out and a streamed buffer is not.
And an `int64` cast that cannot be exact refuses rather than rounding.

*Safety note.* Writing `present_as` on the result is only sound because the result
is never a borrowed argument: `compile.c` refuses to compile a body whose array
result is not a temporary (`CT_IS_ARRAY(res.type) && !res.tmp` bails). So the node
is freshly owned with refcount 1, and the write cannot reach the caller's value.

### What actually cost the 50x

Not the marshalling. The gate. A `CompiledFunction`'s head is an `EXPR_COMPILED`
with no `SymbolDef`, so it read as *unaware* and the gate materialised the very
buffer the boundary exists to borrow — and then the result was materialised again
on the way out. `f[Range[1., 200000.]]` measured 75x slower than
`f[NDArray[Range[1., 200000.]]]`, two values differing only in `present_as`.

*The general lesson:* an allowlist keyed on a symbol table cannot see a head that
is not a symbol. `EXPR_COMPILED` needed the same explicit exemption the pure
`Function` head already had, and for the same reason.

### The prerequisite: a wrong compiled answer

`Compile[{{u, _Integer, 1}}, u * 2][{1, 2, 3}]` gave `{2., 4., 6.}`. A scalar
operand of an array opcode could only be Real or Complex — `ak_of` had no integer
kind — so the literal went into a real register, `vm_box_scalar` boxed it as a
`Real`, and `ndarray_elementwise` then *correctly* widened the `int64` buffer. No
packing involved; a three-element plain List reproduces it. `AK_INT` fixes it, and
had to be fixed first: returning `int64` buffers from the boundary would have made
it far more visible.

---

## 8. What switching it on exposed

Six wrong answers, all reachable from code that never mentions an array, and
none of them visible by reading the code. They were found by a mechanical sweep:
evaluate the same expression over a packed value and over the identical plain
list (`MATHILDA_NO_PACK=1`), for all 134 aware-head cases and all 166 registered
kernels, and diff. Worth recording because each marks a place where the
transparency abstraction leaks, and the leaks were not where the plan predicted.

**1. The evaluator's own fixed-point test undid the gate.** `evaluate` stops when
`expr_eq(current, next)`, and `expr_eq` is *deliberately* blind to whether a list
is packed — that blindness is the feature that makes `Association` lookups and
`Union` work. So a step whose only effect was materialising a buffer looked like
no progress, and `current` (still packed) was kept over `next` (materialised).
`Table[i j, {i,300}, {j,300}]` came back as a `List` of 300 packed rows:
`Dimensions` `{300}`, and `Total[m, 2]` a list instead of a number. Fixed with a
monotonic counter that the gate bumps and `evaluate`'s loop brackets each step
with; when the two are equal *and* the gate fired, take the newer one. The
counter can only over-report (a nested `evaluate`'s gate bumps it too), and an
over-report costs one pointer choice between two equal values.

*The general lesson:* an invariant enforced by mutation is not safe under a
convergence test that cannot see the mutation.

**2. A `Listable` head holding both a plain list and a buffer.** `has_list_arg`
tests for an `EXPR_FUNCTION` headed `List` and does not see a buffer, so
threading fired on the plain list and treated the buffer as a scalar:
`{1.,2.,3.} * {10,20,30}` is `{10., 40., 90.}`, but with the left side packed it
produced the 3x3 outer product. The gate's original design note anticipated the
*opposite* direction (an unaware Listable head skipping threading) and fixed it
by gating first; this is the mixed case, where being aware is what breaks it. The
gate now materialises for a Listable head that also holds a plain list.
Pre-existing for an explicit `NDArray[...]`; packing is what made it reachable.
Surfaced as a `ListConvolve` reference computation off by 103.

**3. A repack coerced elements the operation had introduced.**
`ndstruct_delist_repack` and `map_try_repack` repacked at the *source's* dtype,
so `Join[Range[1.,300.], {1}]` gave `1.` for the appended exact `1` and
`Riffle[Range[1.,300.], 1]` turned every interleaved `1` into `1.`. Now
re-sniffs, via `pack_repack_like`, and declines on a mixed result. A visible
`NDArray[...]` still coerces, which is what naming that head asks for.

**4. An internal `evaluate()` bypasses the gate entirely.** The gate protects a
builtin's *arguments*, because the evaluator puts them there; nothing protects
the value a builtin gets back from its own `evaluate()` call. `Median` builds
`Sort[data]`, evaluates it, then indexes the result's args — and once `Sort`
packed, `Median[Range[300]]` came back UNEVALUATED. The union overlay did its job
here (the code tests `EXPR_FUNCTION` and takes the else branch, so this failed
loudly rather than silently), but loud is still broken. Fixed with
`pack_eval_plain` at the two `stats.c` sites and a blanket materialise in
`internal_call_impl`, which covers 181 wrappers and ~721 call sites at once.

The audit surface turned out to be tiny, which is the useful part: only a head
that can now *produce* a packed list matters, and grepping for internal
constructions of those heads found three sites in the whole tree.

**5. Six kernels answer with different element HEADS than the list does.**
`Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign` and `Im` are registered as
`real_closed`, so they keep the float64 dtype and write `1.0` where the list
gives the exact Integer `1`. Plus `Clip` (clamps to its *exact* bounds, so
`Clip[{1.,2.,3.}]` is `{1., 1, 1}`) and `Precision`/`Accuracy` (`Listable` in
Mathilda, so the list answers once per element while the buffer path gave a
single scalar). These acquire `packed_aware` *implicitly* from
`symtab_set_ndarray_*_kernel`, so they needed clearing rather than omitting —
hence `symtab_clear_packed_aware` and the `NOT_AWARE` list in `pack.c`.

*The general lesson:* an implicit opt-in is the opposite of what D2 wanted.
Having a machine kernel is not the same claim as "my answer's element heads
match the list's", and the kernel registration conflated them.

**6. A resting expression could still hold a buffer.** An aware `Listable` head
skips threading so its kernel can fire; when no kernel matches the call's arity —
`Mod[buffer]`, a one-argument call on a binary kernel — the node came to rest
un-threaded, as `Mod[{1.,2.,3.}]` where the plain list gives
`{Mod[1.], Mod[2.], Mod[3.]}`. Fixed with a POST-GATE at `evaluate_step`'s final
return: if nothing rewrote the node, every head above has declined the buffer, so
materialise. Converges in one extra pass, and it closes the class rather than
this instance — including the `gg[xx][packed]` hazard the original gate note had
to argue about.

Two more, in the tooling rather than the system: `expr_to_latex` had no arm for a
buffer at all (so every packed result came back with an empty `latex` field over
the NDJSON pipe, as had every visible `NDArray[...]` since it was written), and
`bench_assoc` calibrated on `Total[Range[n]]` — which packing makes ~114x faster,
failing all nine of its operations at once with nothing actually slower. That is
the same landmine Phase 0 had already defused in `bench_eval`, in a file nobody
thought to check.

---

## 9. The int64 gate

Packing infers `NDT_INT64` for an all-`Integer` list, which puts integer buffers
in front of code written when only `Compile[]` could make one. `src/expr.h`
described that containment as load-bearing, and it was: the generic
`ndt_get`/`ndt_set` pair routes through `double`, exact only to 2^53, and most
of the ND layer uses it freely.

Two wrong answers were live as soon as `ToNDArray` existed:

```mathematica
Total[ToNDArray[{1, 2, 3}]]   (* 6.  -- should be the Integer 6 *)
Sin[ToNDArray[{1, 2, 3}]]     (* {0, 0, 0} -- should be {Sin[1], Sin[2], Sin[3]} *)
```

Rather than audit ~75 lossy sites before anything could ship, the gate carries a
second, narrower bit: `SymbolDef.packed_int64_ok`. An int64 packed argument is
materialised for every head *not* on that list, so the answer is the ordinary
list's by construction. The list starts short — the metadata heads that read
rank and dims without touching an element, `Part` and `Normal` which go through
the exact accessors, and the pass-through control heads — and grows only as
sites are given genuine exact int64 paths.

What this buys: the storage win (memory, and the ~30 ns/element/pass evaluator
sweep) applies to integer lists immediately, while every computation on them
stays exact. `Total[Range[10^6]]` over a packed list is `500000500000`, an
Integer, promoting past int64 exactly as the interpreter does.

What it costs: integer arithmetic over a packed list currently runs at ordinary
list speed. That is the backlog in §6.2, and it is a backlog of optimisations,
not of correctness.

---

## 10. The surface after the HPC benchmark (2026-07-30)

Building [`performance.md`](performance.md) — 38 classical HPC kernels, Mathilda
against Wolfram 14.0 — turned out to be the most productive audit of this feature
so far. Every gap below was *correct* and *quiet*: it produced the right answer
on the slow path, broke no test, and would not have been found by reading code.

**The pattern in all of them is the same.** The gate's rule is "materialise for
any head that has not opted in", so a missing opt-in never shows up as a wrong
answer — only as a buffer being unpacked and re-packed around an operation that
could have used it directly. That is a failure mode with no symptom, and the only
instrument that finds it is a timing comparison against a system that does not
have the gap.

| what was not opted in | measured cost |
|---|---|
| the 26 linear-algebra heads | `LinearSolve` of a 90×90 real system **stack-overflowed** (below) |
| `Nest`/`NestList`/`NestWhile`/`FixedPoint` | `Nest[f, u0, 10]` 2.19 s vs 0.017 s for the identical `Do` loop — 118× |
| `RotateLeft`/`Join`/`Partition`/`Differences`/`Riffle`/`Pad` | 30×–237× each; `Differences` at 10⁶ was 834 ms |
| a user helper `jac[u_] := …` | the classical Jacobi stencil, 21.6 s vs 0.137 s — 158× |

The linear-algebra one was not merely slow. Those heads dispatch on
`linalg_call_has_ndarray(res)`, which is tag-based; with the gate materialising
first, that test read *false* and a machine-real solve went down the exact
fraction-free path, one polynomial GCD per pivot, recursing until the stack ran
out. **A missing opt-in is normally a performance bug and was here a crash**,
because the fast path and the correct algorithm were selected by the same test.

### The DownValue exemption

The last of these needed a new idea rather than a list entry. A user symbol has
no `packed_aware` bit and cannot be given one, so `jac[packedArray]` materialised
— rightly in general, since `src/match.c` cannot descend a buffer to match
`f[{a_, b_}]`.

But the shape numerical code actually uses is a *bare* pattern variable: bound,
substituted whole, never looked inside. The gate now exempts a head whose every
DownValue **binds opaquely** — each top-level argument of every rule's LHS is
exactly `Pattern[sym, Blank[]]`, with no head restriction, no `PatternTest`, no
`Condition` and no nested structure (`dv_binds_opaquely` in `src/eval.c`). A
head-restricted `_h` looks at the head, a literal `{a_, b_}` looks inside, a
sequence pattern looks at the length; all still materialise.

This is the first place the gate reasons about *patterns* rather than about
heads, and it is deliberately the narrowest such rule that covers
`f[x_] := body`. Anything it does not recognise keeps today's behaviour exactly.

### What the exactness rule cost, and why it is still right

Four of the new buffer paths introduce an element that was not in the input —
`Riffle`'s separator, `PadLeft`/`PadRight`'s fill. Each may only use the buffer
when that element is exactly representable at the dtype *with a matching head*,
so `PadLeft[{1., 2., 3.}, 5]` — which is `{0, 0, 1., 2., 3.}`, exact zeros beside
Reals — declines and stays a boxed list at 1.19 s for 10⁷ elements.

That is not a Mathilda quirk: Mathematica gives the same mixed list and its own
`PadRight` is likewise unpacked, at 483 ms. Given an explicit Real fill both
systems pack and Mathilda is the faster of the two (18.7 ms vs 24.1 ms, and
237× faster than before this work). The rule costs one benchmark row and buys
the guarantee that packing never changes an element's head.

### Still open

- **No narrowing (float64 → int64) kernel category.** `UnitStep` answers with
  exact Integers, so it cannot have an ordinary real-closed kernel, and has none
  at all — 500 ns/element. That single gap costs Game of Life 658× and the
  vectorised Monte Carlo π 27×, and the same category would let `Floor`,
  `Ceiling`, `Round`, `IntegerPart` and `Sign` off `NOT_AWARE`.
- **`Fourier` does not reach FFTW** for a packed buffer — 40–80× against
  Mathematica, and a constant factor rather than an algorithmic one.
- **Serial buffer ops.** `Reverse`, `RotateLeft`, `Differences`, `Accumulate`
  and `Dot` now sit in a 2–4× band against Mathematica's vectorised versions.
  These are memory-bound; the remaining work is SIMD and threading, not
  algorithms.
