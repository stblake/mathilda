# LinkedList: design

Status: **design, not yet implemented** (2026-11-19). Implementation is phased; see §12.

## 1. Motivation

Mathilda has no cheap way to build a list one element at a time.

| Idiom | Cost of n appends | Limitation |
|---|---|---|
| `AppendTo[s, e]` on a List | O(n²): every call copies the list | quadratic |
| `s = {s, e}` then `Flatten` | O(n) | breaks when the elements are themselves lists |
| `Sow`/`Reap` | O(n) | only inside a `Reap`, so the accumulator cannot be kept, inspected, or shared |

`LinkedList` is a first-class, mutable sequence object:
- appending and prepending are O(1) amortised, and so is removing at either end;
- indexing is O(1);
- it can be held in a variable, passed around and inspected while it is being built;
- it works with the heads people already use on lists, and it runs inside `Compile[]` and auto-compiled loops.

## 2. Design decisions

1. **One public head, `LinkedList`.** There is no family of companion functions (`LinkedListAppend`, `LinkedListPart`, …). The existing list heads are overloaded instead, as SparseArray already does: an atom that behaves like a list.
2. **Closure.** A head that returns a `List` for a `List` returns a `LinkedList` for a `LinkedList`.
3. **Mutation is explicit.** Only the `*To` forms, Part assignment, and the new `Push`/`Pop`/`PopFirst` modify a list in place. Every other head leaves its argument unchanged and returns a new LinkedList, as it would for a List.
4. **Full Part specifications.** `Part` (reading and assigning) and `Delete` accept everything they accept for a List: integers, negative integers, `All`, `Span`s with steps, index lists, and several levels.
5. **Compiled.** LinkedList operations lower in `Compile[]` and in auto-compilation (Table/Map/Sum/… and Do/For/While loops).
6. **Stack and queue operations.** `Push`, `Pop` and `PopFirst` are new heads. Without them, removing an element would mean a copying `Most`/`Rest`, and a stack or queue loop would be O(n²) even though the storage can remove from either end in O(1).

## 3. Surface semantics

### 3.1 The object

```
ll = LinkedList[]                 (* empty *)
ll = LinkedList[{1, 2, 3}]        (* from a List, a packed array or an NDArray *)
ll                                (* prints LinkedList[<3>] *)
```

- **Atom.** `AtomQ[ll]` is True, `Head[ll]` is `LinkedList`, and `MatchQ[ll, _LinkedList]` is True. It prints as `LinkedList[<n>]` in StandardForm, InputForm and FullForm, because a mutable object has no readable literal form.
- **Reference semantics.** `ll2 = ll` aliases the same list: a change through either name is visible through both. Passing `ll` to a function passes the same list, so a helper such as `enqueue[q_, x_] := Push[q, x]` changes its caller's list. A LinkedList is the one shared, mutable value in the system; §3.4 records why.
- **Copy constructor.** `LinkedList[ll]` returns an independent **shallow** copy in O(n), copied chunk by chunk with memcpy, with the same storage kind. Boxed elements are shared by reference count, and a LinkedList stored inside `ll` is still the same object in the copy. There is no deep copy. Any copying head (`ll[[All]]`, `Take[ll, All]`, …) also returns a copy.
- **Equality is identity.** `ll === ll2` is True exactly when they are the same object, and a list compares only with itself. Two different lists with equal elements are not `SameQ`. Compare their contents with `Normal[a] === Normal[b]`.
- **Bad constructor arguments.** `LinkedList[x]`, where `x` is neither a list nor a LinkedList, raises `LinkedList::list` and stays unevaluated. More than one argument raises `LinkedList::argt`.

### 3.2 Overloaded heads

Legend:
- **R** reads the list.
- **M** mutates it in place and returns the list.
- **C** copies: it returns a new LinkedList and leaves the argument alone.

| Head | Kind | Result for a LinkedList argument |
|---|---|---|
| `Length[ll]` | R | Integer |
| `First[ll]`, `Last[ll]` | R | element (an empty list gives the same message as `First[{}]`) |
| `ll[[i]]`, `ll[[-i]]` | R | element; any further levels apply ordinary `Part` to that element |
| `ll[[All]]`, `ll[[i;;j;;k]]`, `ll[[{i1, i2, …}]]` | R | **LinkedList** of the selected elements, with further levels mapped over each one |
| `Normal[ll]` | R | `List`, packed when the storage is typed (§4.3) |
| `MemberQ`, `Count`, `Total` | R | same as for `Normal[ll]` |
| `AppendTo[ll, e]`, `PrependTo[ll, e]` | M | `ll`, O(1) amortised |
| `ll[[spec…]] = v` | M | `v`; the full Part-assignment spec, including spans and lists; `ll[[i, j]] = v` rewrites element `i` copy-on-write |
| `Push[ll, e]` | M | `ll`; same as `AppendTo` |
| `Pop[ll]` | M | the removed **last** element (LIFO), O(1) |
| `PopFirst[ll]` | M | the removed **first** element (FIFO with `Push`), O(1) |
| `Append`, `Prepend`, `Insert`, `ReplacePart` | C | LinkedList |
| `Delete[ll, spec]` | C | LinkedList; full Delete spec: `n`, `-n`, `{n}`, `{{n1}, {n2}}`, `{i, j}` (deletes inside element `i`) |
| `Take`, `Drop`, `Most`, `Rest`, `Reverse` | C | LinkedList |
| `Join[a, b, …]` | C | LinkedList when any argument is a LinkedList (any mix with Lists) |
| `Map[f, ll]` (level 1), `Select[ll, crit]` | C | LinkedList |
| Listable heads (`Sin[ll]`, `ll + 1`, `ll1 * ll2`) | C | LinkedList, threaded as over a List |

**Mutators and their first argument.** The mutators are `HoldFirst` (`AppendTo`, `PrependTo` and `Part` assignment already are; `Push`/`Pop`/`PopFirst` are too). The first argument may be a symbol holding the list, or any expression that evaluates to the object. Because the object is a reference, the symbol is never rebound: the payload it points to changes.

**Plain Lists.** `Push`, `Pop` and `PopFirst` also accept a symbol whose value is a plain List. They rebind it the way `AppendTo` does, which is O(n), so code written for a stack works on either representation.

**Every other head** (phase 6, the "LinkedList everywhere" campaign of 2026-11) treats a LinkedList as the List of its elements, so `H[ll, ...]` answers what `H[Normal[ll], ...]` answers. The table above lists the heads that predate that campaign. The rules that govern the rest are below.

- **No conversion.** A head reads the chunks through `ll_span_at` / `LLCursor`, `ll_elem_view` or `ll_seq_element`. It never builds a List, packed List or NDArray of the elements; only `Normal` does. The counter `Internal`LinkedListConversionCount[]` makes this checkable, and `make check-linkedlist` asserts it for every head.
  - An algorithm that needs random access across the whole sequence (sorting, a median, an FFT) may gather the slots into a raw scratch buffer (`ll_gather`).
  - The machine engines take their operands through the container-free `NDOperand` cores (`src/ndarray.h`), which is how Listable threading reaches the very kernels a packed List does.
- **Bit-identical to `Normal`.** `Normal[ll]` of typed storage is packed, and where the List heads answer a packed and a plain List differently, a LinkedList follows its own `Normal`.
- **Closure.** A head that returns a List for a List returns a LinkedList for a LinkedList.
  - Sub-sequences of the argument (Partition, Split and Gather groups, Subsets, Cases, Select) are LinkedLists.
  - Fixed-shape records (Position index lists, Tally pairs, MinMax, Dimensions) and constructors such as `Table` and `Level` stay Lists.
- **A rank-1 container, not a List.** VectorQ, MatrixQ, ArrayQ, Dimensions, Depth, Level, Position, Cases, FreeQ, Apply, Replace and the iterators see its elements one level down, with no head part. `ListQ[ll]` is False and `{__}` does not match it.
- **Snapshots.** A head that runs user code per element (Map, Select, Scan, Fold, the tests, the iterators `{i, ll}`, a `:>` replacement) reads a shallow `ll_copy` taken first. A callback that changes the list therefore changes neither the walk nor its length, just as for a List: `Do[AppendTo[ll, i], {i, ll}]` runs once per original element.
- **Identity.** ReplaceAll and Replace return the list itself when nothing matched. Any change gives a new LinkedList.
- **Values are stored evaluated.** The evaluator never looks inside a LinkedList, so a walker that builds new elements (Map at a level, MapIndexed, MapAt, Apply below level 0) evaluates each one as it stores it. A Throw or Abort raised there leaves the head instead of being stored.

### 3.3 Examples

```
ll = LinkedList[];
Do[AppendTo[ll, i^2], {i, 5}];
Length[ll]                    (* 5 *)
ll[[2 ;; 4]]                  (* LinkedList[<3>] *)
Normal[ll[[-2 ;;]]]           (* {16, 25} *)
ll[[1]] = 100; Normal[ll]     (* {100, 4, 9, 16, 25} *)
Normal[Delete[ll, {{1}, {3}}]](* {4, 16, 25}; ll unchanged *)

s = LinkedList[]; Push[s, a]; Push[s, b];
Pop[s]                        (* b *)
PopFirst[s]                   (* a *)
Pop[s]                        (* Pop::empty, stays Pop[s] *)
```

### 3.4 Decision: assignment shares, copying is explicit

Three models were considered.

1. **Pointer (chosen).** Assignment, argument passing, `With`, rule replacement and return values all share the one list.
   - Assignment is O(1), and costs are predictable with no hidden copies.
   - Helpers can change their caller's list, which stack, queue and graph-search code relies on.
   - Compile arguments and auto-compiled globals work without write-back.
   - The cost is that aliasing is visible, so it is documented here, in the docstring and in the user docs, and `LinkedList[ll]` is the obvious way to copy.
2. **`Set` copies, nothing else does (rejected).** Values reach other places without `Set`: function arguments, `With`, `{ll, ll}`, Part taken from a container, return values. So `b = ll; Push[b, 1]` would leave `ll` alone while `f[ll]` changed it. That inconsistency is worse than either clean model, and `cur = ll` inside a loop would become O(n).
3. **Value semantics with copy-on-write (rejected).** Mutate in place only when nothing else refers to the list; otherwise copy and rebind.
   - Functions could no longer mutate a list (`q = enqueue[q, x]` everywhere).
   - "Nothing else refers to it" is fragile: `Out[n]`, memo tables and evaluator temporaries hold references, which turns an apparent O(1) append into a silent O(n) copy.
   - If value semantics is ever wanted, it belongs in plain Lists (spare capacity, so `AppendTo` on a List nobody else refers to is amortised O(1)), not in a new type.

## 4. Storage

### 4.1 Layout: a deque of geometric chunk runs

A cons-cell list (one node per element) gives O(1) append but O(n) indexing, one allocation per element, and poor locality. The full Part spec needs indexing, so a LinkedList is stored as **two unrolled linked lists of chunks** whose capacities grow geometrically. Chunk *k* holds `B·2^k` elements, with `B = 8`.

```
            front run (grows leftward on PrependTo)     back run (grows rightward on AppendTo)
 logical:   [ f_{F-1} ... f_1 f_0 ]                     [ b_0 b_1 b_2 ... b_{N-1} ]
 index i:     0 ... F-1                                   F ... F+N-1

 back run:  chunk0 (8) -> chunk1 (16) -> chunk2 (32) -> ...   + dir[64] of chunk pointers
```

- Each run keeps `head` and `tail` pointers, a start offset (O(1) removal at the outer end), a length, and a fixed directory `dir[64]` of its chunk pointers. 64 chunks of doubling size cover any addressable length.
- **Locating element *j* of a run.** It is in chunk `k = floor(log2(j/B + 1))`, computed with one count-leading-zeros instruction, at offset `j − B·(2^k − 1)`. There is no pointer chasing.
- **Logical index *i* (0-based).** If `i < F`, the element is `front[F − 1 − i]`; otherwise it is `back[i − F]`. Both positions are adjusted by their run's start offset.
- **Elements never move.** An append or prepend writes into the run's tail chunk and links a new chunk when that one is full.

### 4.2 Complexity

| Operation | Cost |
|---|---|
| `AppendTo`, `PrependTo`, `Push` | O(1) amortised, with no element copied; one allocation per chunk, so O(log n) in total |
| `Pop`, `PopFirst` | O(1); a chunk that becomes empty is freed |
| `ll[[i]]`, `ll[[i]] = v`, `First`, `Last`, `Length` | O(1) |
| Span, list or `All` selection | O(k) for k selected elements |
| Copying heads | O(n) (they build a new list) |
| `Normal` on typed storage | O(n), using memcpy per chunk |

### 4.3 Typed storage

A list's `kind` is one of `EMPTY`, `I64`, `REAL`, `CPLX` or `BOXED`.
- **Typed kinds.** While every element is a machine number of one kind, chunks hold raw `int64_t`, `double` or (re, im) `double` pairs. This is the same rule that decides whether a List packs. Integers and reals do not mix silently: a real appended to an `I64` list makes it `BOXED`, just as `{1, 2.5}` does not pack to an integer array.
- **Conversion to boxed.** The first element that does not conform converts the list to `BOXED` once, in O(n). Boxed chunks hold `Expr*`, each owned through `expr_copy`, which is a reference-count increment. A list never converts back to a typed kind.
- **What typed storage buys.** `Normal` hands back a packed array built by memcpy; `LinkedList[packed]`, `Join[ll, packed]` and packed right-hand sides of Part assignment copy buffers in directly; Listable threading runs the machine kernels (§7) chunk by chunk; and the compiled code (§8) reads and writes raw slots.

### 4.4 Payload header

The header holds:
- a reference count (the payload's own count, separate from the Expr node's);
- `length` and `kind`;
- the two runs;
- a stable `id`, taken from a global counter, used for hashing and ordering;
- the bytes charged so far.

Every chunk allocation and release goes through `MTH_MEM_CHARGE`/`MTH_MEM_RELEASE`, so `MemoryConstrained` limits apply to lists.

### 4.5 Abort safety

An element is fully written, and for a boxed element its ownership taken, before `length` is increased. An `Abort[]` or `TimeConstrained` interruption part way through an operation therefore leaves a list that is consistent and safe to free. The operations do not allocate after the point of publication.

## 5. Object model

### 5.1 The Expr node

A new atom type `EXPR_LINKEDLIST` with the union member `struct LinkedList* llist`, modelled on `EXPR_COMPILED`:

| Operation | Behaviour |
|---|---|
| `expr_new_linkedlist(ll)` | takes ownership of one payload reference |
| `expr_copy` | reference-count increment on the node (unchanged; every type already works this way) |
| `expr_unshare` | creates a new node that **shares** the payload (retain) |
| `expr_free` | releases the payload |
| `expr_eq` | payload pointer identity |
| `expr_hash` | mixes in the payload `id` |
| `expr_compare` | orders by `id`; atoms sort among the other opaque types |
| `Head` / `_LinkedList` | `SYM_LinkedList` |

**Why identity, not value.** Hashes are cached in nodes, including in parent nodes. If equality and hashing were value-based, every mutation would leave stale hashes on every List or Association that holds the list. With identity semantics a mutation never invalidates anything.

**The rule for mutators.** Never `expr_unshare` before mutating. Copy-on-write is the rule for *values*; for a LinkedList the shared payload *is* the meaning.

### 5.2 Cycles

A mutator could store a LinkedList into itself, either directly (`AppendTo[ll, ll]`) or through a chain of lists (`a` contains `b` and `b` contains `a`). Either way it creates a reference-count cycle that is never freed.
- **The check.** `AppendTo`, `PrependTo`, `Push` and Part assignment scan a boxed value for LinkedList atoms. Typed and plain numeric values are skipped, so the scan costs nothing on the numeric path. For each LinkedList found, the mutator checks whether the target can be reached from it.
- **The response.** On a cycle the call emits `General::llcyc` under the calling head and stays unevaluated.
- **Copying heads** always build a fresh list, so they cannot create a cycle.

### 5.3 Ownership of returned elements

`Pop` and `PopFirst` move the element out of the list. A boxed element changes owner with no copy and no free; a typed element is boxed as a fresh machine number. Part reads return `expr_copy` of the stored element.

## 6. Evaluator interaction

The value of `Length[ll]` changes while the OwnValue of `ll` stays the same node, so no result computed from a LinkedList may be reused from a cache. That means a call with an `EXPR_LINKEDLIST` argument:
- is not memoised as a pure builtin call;
- does not get a dependency-recorded no-fire certificate;
- does not get an INERT or GROUND fixed-point stamp.

The gate tests "any argument is `EXPR_LINKEDLIST`" rather than listing heads, so overloads added later are covered automatically. The same mechanism already exempts volatile builtins such as `RandomReal`. A LinkedList atom is never ground.

The regression test for this is `ll = LinkedList[]; Table[AppendTo[ll, i]; Length[ll], {3}]`, which must give `{1, 2, 3}`, and `Length[ll]` must change between two top-level inputs separated by an append.

## 7. Packed and NDArray surfaces

- **Into typed storage.** `LinkedList[packed | NDArray]`, `Join[ll, packed]` and Part assignment with a packed right-hand side copy the buffer straight into typed chunks.
- **Out of typed storage.** `Normal[ll]` builds an NDArray buffer straight from the chunks and returns it through `pack_offer`, so whether it comes back packed follows the usual rule.
- **Listable threading.** A LinkedList argument threads like a List (every LinkedList argument the same length, the rest scalars) and the result is a LinkedList. There are three paths, in order:
  1. Machine reals through a unary kernel, `Plus` or `Times` go chunk to chunk (`ll_listable_direct`).
  2. Every other registered kernel, and the buffer arithmetic of `Plus`, `Times` and `Power`, runs through the engines' `NDOperand` cores (`ll_listable_engine`). The slots are handed over in place when one ascending span holds them, else gathered once into a raw buffer.
  3. Anything else evaluates `head[a1_i, ...]` element by element from snapshots.

  The kernels are the ones that serve packed arrays and `Compile[]`, so a LinkedList, a packed List and compiled code give bit-identical results.
- **`AWARE` list.** `LinkedList` goes on `AWARE` (and `INT64_OK`) in `src/pack.c`, and `Join`, `Append` and `Part` must be on it too, so packed operands reach the LinkedList code without being unpacked first.
- **Audit probes.** `check-nd-surfaces` and `check-array-exactness` get LinkedList probes over all three input representations.

## 8. Compile[] and auto-compilation

As built in phase 5 (v0.6634). Where this differs from the original plan, the plan's version and the reason for the change are given.

### 8.1 Type

`CT_LLIST + elem` is 44..47, above the arbitrary-precision types (40..42). The base is 44 so that `elem` sits in the low two bits, as for the array and Association encodings. `elem` is Integer, Real or Complex. Elem 0 (`CT_BOOL`) is the argument parser's marker for "element type not given" and never reaches the emitter. Every existing scalar, array, Association and managed-type check rejects these values; `num_common`, `coerce`, `scalar_only`, the inliner and `OP_CALL` reject them explicitly. `CompileDiagnostics` reports `ResultType -> "LinkedList"`.

A register holds the `EXPR_LINKEDLIST` node in `.arr`:
- **owned, in the array bank:** a list the body created, a Module local, or a retained reference to an argument or a folded global. The per-call sweep releases these on every exit path.
- **borrowed, in an argument register:** a list passed in. It is the caller's own object and is never released.

### 8.2 Element type: one per list, never widened

Every value pushed into a compiled list must have the list's element type. The plan joined Int and Real to Real, but that changes the answer: the interpreter keeps `{1, 2.5}` as a mixed, boxed list, and a widened compiled list would give `{1., 2.5}`.

The element type is found as follows:
- **A Module local `b = LinkedList[]`** has no type until something is pushed, so it is guessed, Real first. A push of another machine type records what it wants (`ll_want_*`), and the enclosing Module lowers its body again with that type (`llist_retry`, at most three attempts). If a second, different type is pushed, the list is mixed: the body bails and the interpreter answers.
- **`b = LinkedList[arr]`** takes `arr`'s element type.
- **An argument** is declared `{b, _LinkedList}` or `{b, _LinkedList, _Integer | _Real | _Complex}`, following `{p, _Association, _Real}`. Without the element type, `Compile[]` tries Real, then Integer, then Complex, and keeps the first that compiles. At call time a non-empty list of another kind runs through the interpreter.
- **A folded global** takes its current kind. If it is still empty, the loop runner tries Real, Integer and Complex in turn (`COMPILE_LL_EMPTY_INT`, `COMPILE_LL_EMPTY_CPLX`). A boxed global is not compiled.

### 8.3 Lowered operations

Every opcode calls the shared C helpers the interpreter uses (`ll_push_back_*`, `ll_get_*`, `ll_set_*`, `ll_pop_typed`, `ll_sum_machine`, `ll_append_range`, `ll_to_ndarray`), so the compiled and interpreted results are the same by construction. Every opcode is `K_ARR`:
- The optimiser treats it as opaque: it never moves, removes or CSEs one.
- A loop containing one is left out of LICM, because a list is mutable shared state, so even `Length` must not be hoisted past a push.
- The VM dispatches it out of line (`compile_llist_vm.c`), like the other array ops.
- The strip-mined parallel loops are generated only for fused array code, so a loop holding a list op never fans out across threads.

| Surface | Opcode | Notes |
|---|---|---|
| `LinkedList[]`, `LinkedList[arr]` | `LL_NEW`, `LL_FROMARR` | |
| `b = a` (a local bound to a list) | `LL_RETAIN` | the same list, by reference |
| a folded global | `LL_REF` | program-owned reference (`p->llrefs`) |
| `AppendTo`/`Push`, `PrependTo` | `LL_PUSH` (`LLF_FRONT`) | value is the list |
| `Pop`, `PopFirst` | `LL_POP` | empty: the call fails, and the interpreter answers and emits `Pop::empty` |
| `b[[i]] = v` (one integer index) | `LL_SET` | value is `v` |
| `b[[i]]`, `First`, `Last` | `LL_GET` | out of range: the call fails |
| `Length`, `Total` | `LL_LEN`, `LL_TOTAL` | an Integer Total past int64 fails; a Real or Complex Total is the packed tree, bit for bit |
| `Normal[b]` | `LL_NORMAL` | rank-1 array, built (so it is packed only past the size threshold); empty: fails |
| `Most`, `Rest`, `Reverse`, `Take[b, k]`, `Drop[b, k]`, `Append`, `Prepend`, `Join[b1, b2]` | `LL_STRUCT` + sub-op | new list, chunk span copies |

A compiled function may return a LinkedList. A returned argument is the same object (`LL_RETAIN`), as in the interpreter.

### 8.4 Failure and the undo log

A compiled call that cannot finish is abandoned, and the interpreter runs the body again from the start. Causes include an overflow, an index out of range, an empty `Pop`, or a result the boundary declines (underflow). That was sound only while compiled code had no side effects, and a list the run did not create is a side effect.

- **The undo log.** A program with any list op runs inside an undo frame (`ll_undo_begin` / `_commit` / `_rollback`, `src/linkedlist.c`). While a frame is open, every push, pop and Part write to a list created before the frame is logged, and a rollback reverts them in reverse order.
- **Which lists.** "Created before" is one compare, because list ids are issued in creation order. Lists the run created need no record, since they are thrown away too.
- **Merging.** Consecutive pushes onto the same end of one list merge into a single record, so `Do[AppendTo[ll, x], {10^6}]` logs one record.
- **Nesting.** Frames nest. An inner commit keeps its records for the outer frame, and only the outermost commit discards the log.
- **Where frames open.** `compiled_eval` opens one, and `compiled_function_apply` opens a second around its post-run checks (underflow, degenerate entries), so a result declined after a successful run is also undone.

The test: a compiled function that appends, prepends, writes and pops, then reads past the end, leaves the caller's list exactly as one interpreted run leaves it, with one `Part::partw`.

### 8.5 Who may mutate a list from outside

An undo frame covers one run of one program, and that is enough for a `Compile[]` call and for a whole loop. It is not enough for a per-element auto-compiled program (`Table`, `Map`, `Sum`). When a later element fails there after earlier ones succeeded, the interpreter redoes the whole call, earlier elements included. The plan had expected `Table[AppendTo[ll, f[i]]; …]` to compile; it must not. So:
- **Reads are allowed everywhere.** A folded global list may be read in any auto-compiled body (`Table[ll[[i]]^2, {i, n}]`).
- **Mutating a folded global, or binding one to a local** (which could then be mutated), needs `COMPILE_LL_MUTATE_GLOBALS`. Only the whole-loop runner sets that flag.
- **No pasting.** A `CompiledFunction` with a list parameter is never inlined into, or called from, another program, so it always runs as its own call with its own frame.

### 8.6 Auto-compilation of Do, For and While

`numloop` declines a body that is not a block of numeric assignments. When it does, and the loop names a global whose value is a LinkedList (a depth-bounded scan), `src/iter.c` hands the whole loop to `autocompile_run_closed_stmt`, which:
- compiles the loop as one program with folded globals, `COMPILE_LL_MUTATE_GLOBALS` and `COMPILE_STRICT_JOIN`;
- runs it inside an undo frame;
- returns `Null` on success, or declines, after which the interpreter runs the loop.

`ll = LinkedList[]; Do[AppendTo[ll, N[i]], {i, 10^6}]` takes 37 ms, against about 0.5 s interpreted.

**Strict joins.** `Compile[]` unifies the branches of `If`, `Which`, `Switch` and `Piecewise`, and the values assigned to a local, as the reference's Compile does: `If[EvenQ[i], i, N[i]]` is always Real. The loop runner stands in for the interpreter, whose pushes are observable, so under `COMPILE_STRICT_JOIN` a join of two types bails rather than widening, and the list keeps `{1., 2, 3., 4}`.

**`N`.** `N[x]` on a machine number had no lowering, which sent every `Do[AppendTo[ll, N[i]], …]` to the interpreter. It now lowers: an Integer becomes the Real `(double) x`, which is the interpreter's conversion; a Real or Complex passes through; an Integer array becomes `1. v`.

### 8.7 Out of scope inside Compile

These bail cleanly to the interpreter:
- boxed or nested elements, including lists of lists and lists of LinkedLists;
- `Delete`, `Insert` and `ReplacePart`;
- `Map`, `Select` and Listable threading over a list;
- deeper Part, sequence specs in `Take`/`Drop`, and Part writes other than a single integer index;
- rebinding a list local, `AddTo` on a list element, and `LinkedList[ll]` (a copy).

## 9. Messages

All messages go through `mth_message` (`src/message.h`), so `Quiet` and `Check` work on them.

| Message | Text | When |
|---|---|---|
| `LinkedList::list` | List expected at position 1 in `` `1` ``. | constructor given a non-list |
| `LinkedList::argt` (and similar) | the standard arity messages via `builtin_arg_error` | wrong argument count |
| `General::llcyc` | Storing `` `1` `` would make the linked list contain itself. | a mutator would create a cycle (§5.2) |
| `Pop::empty`, `PopFirst::empty` | Cannot pop from an empty list `` `1` ``. | popping an empty list |
| existing `Part::partw`, `Part::pkspec`, `Delete::partw`, `Take::take`, … | unchanged | bad specs, reported exactly as for the same spec on a List |

Every error leaves the call unevaluated.

## 9a. Symbol registration and docstrings

Each new head, `LinkedList`, `Push`, `Pop` and `PopFirst`, is registered as a full builtin.

1. **Interned name.** Declare `SYM_LinkedList`, `SYM_Push`, `SYM_Pop` and `SYM_PopFirst` in `src/sym_names.h`. Define them in `src/sym_names.c`, assigned in `sym_names_init()` with `intern_symbol(...)`.
2. **Mathilda-only table.** Add all four names to the table in `sym_mark_mathilda_context()` (`src/sym_names.c`), so a package that defines its own `Push` gets its own symbol instead of writing to the protected builtin.
3. **Builtin, attributes, docstring.** `linkedlist_init()` (called from `core_init` in `src/core.c`) registers each one with `symtab_add_builtin`, sets its attributes, and sets its `?Name` text with `symtab_set_docstring`:

| Symbol | Attributes |
|---|---|
| `LinkedList` | `Protected` |
| `Push` | `HoldFirst`, `Protected` |
| `Pop` | `HoldFirst`, `Protected` |
| `PopFirst` | `HoldFirst`, `Protected` |

4. **REPL.** Once registered, the heads are available in the REPL: Tab completion, `?LinkedList`, and `Names["Pop*"]`.

Draft docstrings, in the usual `Name[args]\n\t…` form:

```
LinkedList[]
	Creates an empty mutable linked list. LinkedList[list] creates one holding
	the elements of list (a List, packed array or NDArray); LinkedList[ll]
	returns an independent shallow copy of the linked list ll. Assignment and
	argument passing share the same list rather than copying it: AppendTo, PrependTo, Push, Pop, PopFirst and Part assignment
	change it in place, while Part, Take, Drop, Append, Delete, Map, ... return a
	new LinkedList. Normal[ll] gives the elements as a List.

Push[ll, e]
	Appends e to the end of the linked list ll in place, in O(1), and returns ll.
	Push[s, e] for a symbol s holding a List appends e and resets s.

Pop[ll]
	Removes the last element of the linked list ll in place, in O(1), and
	returns it. Pop[s] for a symbol s holding a List removes its last element.
	An empty list gives Pop::empty and stays unevaluated.

PopFirst[ll]
	Removes the first element of the linked list ll in place, in O(1), and
	returns it. Push with PopFirst gives a first-in, first-out queue.
	An empty list gives PopFirst::empty and stays unevaluated.
```

**Overloaded heads.** These keep their existing registration. Each one's docstring gains a single line noting that it accepts a LinkedList, and what it returns for one. The heads are `Part`, `Length`, `First`, `Last`, `Normal`, `AppendTo`, `PrependTo`, `Append`, `Prepend`, `Insert`, `Delete`, `ReplacePart`, `Take`, `Drop`, `Most`, `Rest`, `Reverse`, `Join`, `Map`, `Select`, `MemberQ`, `Count` and `Total`.

**Documentation.** Each head gets an entry in `docs/spec/builtins/data-structures.md`, and each phase adds a changelog note.

## 10. Implementation map

| File | Role |
|---|---|
| `src/linkedlist.{c,h}` | payload, chunk runs, typed and boxed storage, shared `ll_*` helpers, the constructor, `Push`/`Pop`/`PopFirst`, `linkedlist_init` |
| `src/linkedlist_ops.c` | one entry per overloaded head (`linkedlist_part`, `linkedlist_delete`, …), called from a one-line hook at the top of each builtin |
| `src/expr.{c,h}`, `src/sort.c`, `src/part.c`, `src/match.c`, `src/print*.c`, `src/eval.c`, `src/numeric.c`, `src/precision.c`, … | the new atom type in every type switch |
| `src/sym_names.{c,h}` | `SYM_LinkedList`, `SYM_Push`, `SYM_Pop`, `SYM_PopFirst`; all four entered in the Mathilda-only name table |
| `src/pack.c` | the `AWARE`/`INT64_OK` entries |
| `src/compile/compile.h`, `compile_infer.c`, `compile_emit_llist.c` (new), `compile_internal.h`, `compile_vm.c`, `disasm.c`, `compiled_function.c` | the type, inference, codegen, opcodes, VM, disassembly and `ResultType` name |
| `src/compile/autocompile.c`, `src/iter.c` | the closed-statement auto-compile path for Do, For and While |
| `docs/spec/builtins/data-structures.md` | user documentation |

## 11. Test strategy

- **Differential test against List.** For every overloaded head and every spec form (Part read and assignment, Delete, Take, Drop, …), compare the result for `ll` with the result for the same call on `Normal[ll]`, allowing the `LinkedList`/`List` head difference. That one harness covers the closure rule, full-spec parity and the error messages.
- **Mutation semantics.** Aliasing; a list that converts from typed to boxed part way through; prepends and appends interleaved; `Push`, `Pop` and `PopFirst` round trips, including on an empty list and on a plain List symbol; rejection of a cycle.
- **Cache regression.** The §6 cases.
- **Messages.** Each message is checked under `Quiet` and under `Check`.
- **Compilation.** `CompileDiagnostics` reports `Compiled -> True` for the Module-local and `_LinkedList`-argument forms; compiled results are bit-identical to the interpreter's; an auto-compiled `Do` loop takes the compiled path.
- **Performance.** 10⁶ appends, and 10⁶ push+pop pairs, run in linear time, close to `Reap`/`Sow`.
- **Gates.** `check-c99`, `check-messages`, `check-packed-aware`, `check-nd-surfaces`, `check-array-exactness`, `check-compile-coverage`, `check-fastpath-sweep`, and `ctest -LE nightly`.
- **Leaks.** valgrind on a script that creates, aliases, mutates, nests, rejects a cycle, aborts an append, and runs a compiled loop. It must show no new definite leaks against the startup baseline.

## 12. Phases

Each phase is one commit with a version bump and a tag.

| Phase | Scope |
|---|---|
| 1 | object, storage, constructor, printing, type switches, cache gating |
| 2 | reading heads, full-spec Part reads, `Normal`, packed surfaces |
| 3 | mutators (`AppendTo`, `PrependTo`, `Push`/`Pop`/`PopFirst`, full-spec Part assignment), cycle check |
| 4 | copying heads, full-spec `Delete`, Listable threading with machine kernels |
| 5 | `Compile[]` lowering, then auto-compilation |
| 6 | everywhere a List works: the audit (`make check-linkedlist`) and conversion counter; Map, Select and Listable reading the slots (`NDOperand` cores, `NumFn`); level walkers, functional and threading heads, structure queries, iterators; sorting, sets and grouping; strings and statistics; vector numerics and symbolic mappers (§3.2) |
