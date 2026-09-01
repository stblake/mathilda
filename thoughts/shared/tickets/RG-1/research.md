---
created: 2026-08-30T21:18:27-0400
researcher: Michael Sollami
source_sha: 7701df5b36cab11000efad1f3b2dff094ff2316f
branch: main
repository: mathilda
topic: "How does RandomGraph work today, and what would it take to add the RandomGraph[{n,m}, k] form that returns k random graphs?"
tags: [research, codebase, graph, generators, random, memory]
subsystems: [graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-30
last_updated_by: Michael Sollami
---

# Research: RandomGraph today, and adding the `RandomGraph[{n,m}, k]` form

**Date**: 2026-08-30T21:18:27-0400
**Researcher**: Michael Sollami
**Git Commit**: 7701df5b36cab11000efad1f3b2dff094ff2316f
**Branch**: main
**Repository**: mathilda

## TL;DR

`RandomGraph` is a single-arity builtin (`generators.c:103-134`) that hard-rejects any
second argument, so the `k` form is a ~15-line change following the `RandomInteger`
precedent. Measurement turned up an unrelated, pre-existing leak: each `RandomGraph`
call leaks its whole candidate-edge tree — 364 KB at n=50. Left unfixed by decision;
it bounds how the `k` form should loop.

## Summary

`RandomGraph[{n, m}]` builds all `n(n-1)/2` candidate `UndirectedEdge` expressions,
sampling `m` of them by constructing and evaluating a real `RandomSample[...]` call —
which is why it honors `SeedRandom`. The arity gate is a literal
`arg_count != 1` at `generators.c:104`, so `RandomGraph[{5,4}, 3]` returns unevaluated
today (measured, not inferred).

Adding the flat `k` form is mechanical. Five sibling heads in `src/random.c` already
implement a trailing count argument, and they agree on every convention that matters:
`k = 0` yields `{}`, and a negative, non-integer, or symbolic count returns `NULL`
(unevaluated) with no `Message[]`. There is no shared helper to reuse — each head
duplicates its own validation loop — so a sixth local copy is the idiomatic choice
here, not a refactor.

The one surprise is a measured pre-existing memory leak in the generator family, which
is recorded as a finding below and deliberately not in scope.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] Flat `k` only, or the nested `{k1,k2}` form too? — Flat `RandomGraph[{n,m}, k]`
  only; the nested form is explicitly out of scope.
- [x] Should the O(n²) candidate-materialization be reworked as part of this? — No.
  Reuse the existing single-graph path in a loop; record the cost as a known limit.
- [x] Was the `k` form deliberately deferred, or is there a prior attempt to avoid
  repeating? — Neither. Greenfield: `src/graph/generators.c` has exactly one commit
  (`56035303`), the graph-subsystem bulk add, and the 1-arg form was simply all that
  got written.
- [x] Is the generator-family memory leak in scope? — No. Record it as a finding with
  file:line and the measured numbers; do not fix it.

## Research Review

### Blocking

_None._

### Worth Flagging

_None._

### Resolved

_None._

## Requires Approval

The measured leak below (Finding 4) interacts with the agreed "reuse as-is in a loop"
approach: at n=50 the existing single-graph path leaks ~364 KB per call, so
`RandomGraph[{50,20}, 1000]` would leak ~364 MB. The scope decision to not fix it
stands, but whoever implements the `k` form should know the leak scales with `k` and
may want a documented cap or a note in the docstring.

---

## Research Question

> How does RandomGraph work today, and what would it take to add the RandomGraph[{n,m}, k] form that returns k random graphs?

## Detailed Findings

### 1. Current behavior — measured, not inferred

Built at `7701df5b` (needed `SDKROOT=$(xcrun --show-sdk-path) make -j8`; a bare `make`
fails with `fatal error: stdio.h: No such file or directory` under gcc-16 in this
environment). Actual REPL output:

| Input | Output |
|---|---|
| `RandomGraph[{5, 4}]` | `Graph[<5 vertices, 4 edges>]` |
| `RandomGraph[{5, 4}, 3]` | `RandomGraph[{5, 4}, 3]` — unevaluated |
| `RandomGraph[{5, 4}, 1]` | `RandomGraph[{5, 4}, 1]` — unevaluated |
| `RandomGraph[{5, 11}]` | unevaluated (m > n(n-1)/2, as documented) |
| `SeedRandom[42]` twice | identical `EdgeList` — reproducibility holds |

So the `k` form is not partially present or silently wrong; it is entirely absent, and
the head is otherwise sound.

### 2. Implementation and the exact arity gate

`src/graph/generators.c:103-134`. The gate is two lines:

```c
Expr* builtin_random_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;     /* :104 */
    const Expr* spec = res->data.function.args[0];
    if (!graph_is_list(spec) || spec->data.function.arg_count != 2) return NULL;
```

Then: build every candidate edge (`:113-121`), sample via a constructed
`RandomSample[cand_list, m]` that is really `evaluate()`d (`:124-127`), and assemble
`Graph[List[1..n], sampled]` (`:130-133`).

Two details worth carrying into the change:

- `builtin_random_graph` does **not** use the file's own `make_graph` helper
  (`:36-41`) — it inlines the assembly at `:130-133`. Every other generator uses
  `make_graph`. A `k`-form refactor is the natural moment to extract the
  single-graph body into a static helper that both the 1-arg and `k`-arg paths call.
- It reaches `RandomSample` through the evaluator rather than a C-level RNG helper.
  That is load-bearing for `SeedRandom` and should be preserved, not optimized away.

### 3. Prior art for the trailing count argument — five copies, no shared helper

`src/random.c` implements this exact shape five times over, each with its own
validation loop and its own recursive array builder:

- `builtin_randominteger` — `src/random.c:418-553`, builder `random_array` at `:370-406`
- `randomreal_machine` — `src/random.c:941-1024` (+ MPFR twin at `:842-936`)
- `randomcomplex_machine` — `src/random.c:1394-1449` (+ MPFR twin at `:1276-1389`)
- `builtin_randomchoice` — `src/random.c:1730-1868`
- `builtin_randomsample` — `src/random.c:1967-2033` (flat `n` and `UpTo[n]` only, no
  nested dims)
- `builtin_randomvariate` — `src/ml/dist.c:465-485` (flat `n` only)
- `builtin_constant_array` — `src/list/constant_array.c:38-110`, the same idea
  structured a third way

There is no shared "parse a count-or-dimension spec" function anywhere; `src/random.h`
exports only the `builtin_*` entry points and two RNG primitives. Given flat-`k`-only
scope, writing a sixth local check is consistent with the codebase, and a
generalizing refactor would be out of scope creep.

**The conventions these agree on** — measured against the live binary for all five
heads, not inferred from one:

| Count | `RandomInteger` | `RandomReal` | `RandomComplex` | `RandomChoice` | `RandomSample` |
|---|---|---|---|---|---|
| `0` | `{}` | `{}` | `{}` | `{}` | `{}` |
| `-1` | uneval | uneval | uneval | uneval | uneval |
| `2.5` | uneval | uneval | — | — | — |
| symbolic | uneval | — | — | uneval | uneval |

Unanimous. Worth noting explicitly because `RandomSample` is the head this change
actually calls into, and § 7 below shows it *does* have an extra guard —
`is_nonempty_list` at `src/random.c:1707`. That guard is on the **list** argument, not
the count: `RandomSample[{a,b}, 0]` → `{}` like everyone else, while
`RandomSample[{}, 0]` is unevaluated. So the count convention is uniform, and the
divergence is confined to the empty-candidate-list case that causes the § 7 bug.

`grep` for `Message(`/`symtab_message` in `src/random.c` returns zero matches — the
uniform convention is silent `NULL`, per SPEC §4. `as_count` in
`generators.c:25-28` already implements exactly this predicate (non-negative machine
`EXPR_INTEGER`, else `-1`), so the `k` form should parse `k` with `as_count` and get
the right behavior for free.

Note `k = 1` should return a **one-element list** `{Graph[...]}`, not a bare graph —
that is the Mathematica semantics and matches `RandomInteger[{1,10}, 1]` → `{7}`.

### 4. FINDING (out of scope, do not fix): the graph generators leak, and `RandomGraph` leaks badly

`expr_new_function` **memcpy**s the caller's argument array rather than adopting it
(`src/expr.c:257`):

```c
if (args) memcpy(e->data.function.args, args, sizeof(Expr*) * arg_count);
```

So the caller must `free()` the C array it passed. `int_vertices` (`:43-47`) `calloc`s
a fresh array every call, and **no** generator frees it:

- `generators.c:59` — `make_graph(int_vertices(n), ...)`, array leaked
- `generators.c:73` — same
- `generators.c:100` — same
- `generators.c:131` — same, inside `builtin_random_graph`
- `generators.c:92` — the `PathGraph[{v1,...}]` explicit-vertex branch, which leaks a
  *separately* allocated `verts` array (`calloc` at `:83`) by the same mechanism

`builtin_random_graph` correctly frees its own `cand` array at `:121`; the `edges`
arrays at `:54`, `:69`, `:86`, `:98` are also never freed.

Per-call array count, precisely: `CompleteGraph`, `CycleGraph`, and both `PathGraph`
branches leak **two** arrays each (vertices + edges). `builtin_random_graph` leaks
**one** — it allocates `cand` rather than an `edges` array and does free `cand` at
`:121`, so only `int_vertices(n)` at `:131` is lost.

Measured with `MallocStackLogging=1 leaks --atExit`, 2000 iterations each:

| Expression | Leaks | Bytes |
|---|---|---|
| `RandomSample[Range[1225], 20]` | **0** | 0 |
| `CycleGraph[50]` | 4,000 | 1,792,000 (~1.8 MB) |
| `CompleteGraph[50]` | 4,000 | 21,376,000 (~21 MB) |
| `RandomGraph[{50, 20}]` | **12,265,994** | **727,679,616 (~727 MB)** |

Two distinct leaks, and the second is the serious one:

- **(A) The C arrays** — 2 per call (`int_vertices` + `edges`), affecting all four
  generators. `leaks` names it directly: `ROOT LEAK: <calloc in
  builtin_complete_graph>`, 2000 instances. Small and bounded: O(n) pointers.
- **(B) The candidate-edge *tree*, `RandomGraph` only** — ~6,066 `Expr` nodes leaked
  per call at n=50, i.e. ~364 KB/call. That is the full 1,225-edge candidate list
  (1225 edges × ~5 nodes each ≈ 6,125). The comment at `:127` asserts
  `evaluate(sample_call)` "consumes sample_call", and the measurement says the
  candidate tree it wraps is not in fact reclaimed. `RandomSample` itself is clean
  (0 leaks above), so this is `builtin_random_graph`'s own defect, not the sampler's.

Why this matters for the `k` form specifically: leak (B) is per-call and the agreed
approach is to call the existing path `k` times, so the leak multiplies by `k`.
`RandomGraph[{50,20}, 1000]` leaks ~364 MB. This does not block the feature — it
bounds it, and it is worth a line in the changelog rather than silence.

### 5. What the change touches

Per the `graph.h` "Phase" convention (one builtin per TU, prototypes in `graph.h`,
registration in `graph.c`):

1. `src/graph/generators.c:103-134` — relax the arity gate; extract the single-graph
   body into a static helper; loop it `k` times into a `List`.
2. `src/graph/graph.h:103` — update the trailing signature comment
   (`/* RandomGraph[{n, m}] */`).
3. `src/graph/graph.c:122-124` — extend the docstring to name the `k` form
   (mandatory per CLAUDE.md). Attributes need no change: `ATTR_PROTECTED` only, and
   the graph module sets attributes via the direct
   `symtab_get_def(name)->attributes |= ...` idiom, never `symtab_set_attributes`.
4. `tests/test_graph.c:213-228` — extend `test_random_graph`; style is
   `assert_eval_eq(input, expected_string, fullform_flag)`.
5. `docs/spec/builtins/graphs.md:115-117` and
   `docs/spec/changelog/2026-08-24.md` (Monday of this ISO week) — required by
   CLAUDE.md.

### 6. The packed/NDArray/Compile rule does not apply here

CLAUDE.md requires every new numeric builtin to support packed arrays, NDArray
kernels, and `Compile[]` lowering. Checked, and it does not bite:

- No `Graph`/`RandomGraph`/`CompleteGraph`/`CycleGraph` entry in `src/pack.c`'s
  `AWARE` list (`src/pack.c:489-789`); the only `Graph` match in that file is an
  unrelated comment at `:780`.
- No matches in `src/ndkernels.c`, `src/ndinteger.c`, or `src/compile/`.
- `RandomGraph` returns a `Graph` object, not a machine array — exactly the
  "non-machine object" case CLAUDE.md names as a legitimate exemption.

One caveat for honesty: `PathGraph` *is* on the known-gap ratchet list in
`tools/nd_fastpath_sweep.py:689-690` (`OFF_BUFFER`), which is a backlog list, not an
excused-exemption list. `RandomGraph` appears in neither, so no audit currently
demands anything of it — but there is also no recorded, explicit `EXEMPT` justification
for it anywhere in the tree. Adding one line to the appropriate audit list would make
the exemption deliberate rather than merely unnoticed.

### 7. Adjacent pre-existing bug: `RandomGraph[{0,0}]` and `{1,0}` are unevaluated

Measured:

```
RandomGraph[{0, 0}]  ->  RandomGraph[{0, 0}]   (unevaluated)
RandomGraph[{1, 0}]  ->  RandomGraph[{1, 0}]   (unevaluated)
RandomGraph[{2, 0}]  ->  Graph[<2 vertices, 0 edges>]
```

Root cause isolated: for n ≤ 1, `maxe = 0`, so `cand_list` is the empty `List[]`.
`builtin_randomsample` guards its input with `is_nonempty_list`
(`src/random.c:1707`, used at `:1751-1753`), so `RandomSample[{}, 0]` is itself
unevaluated — confirmed directly — and `graph_is_list(sampled)` then fails at
`generators.c:128`, returning `NULL`. `Graph[{}, {}]` and `CompleteGraph[0]` both
work fine, so this is specifically `RandomGraph`'s dependence on the evaluator-level
sampler. Not in scope, but a `k`-form implementation that early-returns for
`maxe == 0` would incidentally fix it.

## Code References

- `src/graph/generators.c:103-134` — `builtin_random_graph`, the whole current impl
- `src/graph/generators.c:104` — the `arg_count != 1` gate that rejects `k`
- `src/graph/generators.c:25-28` — `as_count`, the count predicate to reuse for `k`
- `src/graph/generators.c:36-41` — `make_graph`, which `builtin_random_graph` bypasses
- `src/graph/generators.c:43-47` — `int_vertices`, fresh `calloc` per call, leaked
- `src/graph/generators.c:121` — `free(cand)`, the one array that *is* freed
- `src/graph/generators.c:127` — `evaluate(sample_call)`, comment claims it consumes
- `src/graph/graph.h:99-103` — Phase 4 generator prototypes and naming convention
- `src/graph/graph.c:120-124` — registration, `ATTR_PROTECTED`, docstring
- `src/core.c:944-945` — `graph_init()` call site
- `src/expr.c:239-262` — `expr_new_function`; `:257` is the memcpy that implies the leak
- `src/random.c:418-553` — `builtin_randominteger`, the closest model for the `k` form
- `src/random.c:1707` — `is_nonempty_list`, cause of the n≤1 bug
- `src/random.c:1967-2033` — `builtin_randomsample`
- `tests/test_graph.c:213-228` — `test_random_graph`
- `tests/test_graph.c:310-311` — test registration
- `docs/spec/builtins/graphs.md:115-117` — the `RandomGraph` spec entry
- `tools/nd_fastpath_sweep.py:689-690` — `PathGraph` on the `OFF_BUFFER` ratchet

## Architecture Insights

- **The graph subsystem mirrors `src/linalg/`**: one builtin per translation unit,
  prototypes grouped by "Phase" banner in `graph.h:16-18`, registration in `graph.c`
  under matching banners. A new arity is a same-file change, not a new TU.
- **Generators delegate randomness to the evaluator, not to C.** Constructing and
  `evaluate()`-ing a real `RandomSample[...]` call is what makes `SeedRandom`
  reproducibility fall out for free. It is also the source of both defects found here
  (the n≤1 unevaluated case and the leaked candidate tree) — the coupling is
  load-bearing and fragile at the same time.
- **"Unevaluated" is the universal error channel.** Zero `Message[]` calls in
  `src/random.c`; `NULL` return per SPEC §4 is the convention for every bad count
  across seven heads in three subsystems. A new `k` form should not be the first head
  to emit a diagnostic.
- **Duplication is the accepted local idiom** for count-spec parsing — five copies in
  one file. Worth knowing so a reviewer does not read a sixth as sloppiness.

## Historical Context (from thoughts/)

Thin. `thoughts/` holds three documents, none touching the graph subsystem:

- `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md`
- `thoughts/shared/research/2026-08-18-unitbox-builtin.md`
- `thoughts/shared/plans/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarks.md`

No subsystem doc exists for `graph` (`thoughts/shared/subsystems/` is empty), so there
were no pre-declared conventions to fold in. Writing one is a reasonable follow-on
given the Phase/one-builtin-per-TU conventions documented above.

Git history is equally thin and confirms greenfield: `src/graph/generators.c` has a
single commit, `56035303` ("Add graph subsystem: construction, queries, matrix views,
generators, algorithms, visualization"). The `RandomGraph` spec text dates from
`docs/spec/changelog/2026-06-29.md:3567-3568`.

## Related Research

None in `thoughts/shared/research/` touches the graph subsystem.
