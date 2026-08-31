---
ticket: RG-2
created: 2026-08-31
reviewer: ais:adversarial-reviewer (dispatched from /ais:verify-implementation)
source_sha: 854419997e4b0d8b35fb401783a17b13db2495c9 + uncommitted working tree
branch: find-vertex-coloring
scope: Phases 1-2 (Phases 3-4 descoped to RG-5 per plan.md:471-491)
result: no HIGH findings — 3 MEDIUM, 2 LOW
---

# Adversarial review — FindVertexColoring (RG-2 Phase 2)

This is a **log**, not a live list. Findings stay here with their resolution recorded.
A finding with no resolution marker is open.

## Reviewed

`src/graph/vertexcoloring.c` (new), `src/graph/graph.{c,h}`, `src/sym_names.{c,h}`,
`tests/test_graph.c`, `tests/test_graph_slow.c` (new), `tests/CMakeLists.txt`.

## Verified correct (context for the findings below)

- **Minimality holds.** Differentially tested against an independent Python brute-force
  chromatic number over 120 random graphs (n=1..13, all densities) — **zero mismatches**.
  The docstring's "MINIMAL" claim is accurate and the K_{2,2} / greedy-upper-bound trap is
  genuinely avoided.
- **Builtin ownership contract is correct on every path.** `res` is never freed;
  `graph_build_adj` takes `const Expr*` and borrows; all six early returns return `NULL`
  without freeing; `elems` is freed after the memcpy in `expr_new_function`.
- **Self-loops and duplicate edges** are rejected upstream by `graph_check`
  (`src/graph/graph_util.c:314,320`) and cannot reach the search.
- **No RG-4-style stranded-static hazard.** `vertexcoloring.c` has no mutable static or
  global state (its four `static` tokens are all functions). Two successive
  `TimeConstrained` aborts each returned `$Aborted` and the third call answered correctly.

---

## MEDIUM — Every `TimeConstrained` abort leaks ~53 KB, and the docstring points users at exactly that path

**Where:** `src/graph/vertexcoloring.c:270` (the `tc_check_deadline()` poll) and
`:365-394` (the head's `a`/`colour` allocations); docstring at `src/graph/graph.c:177-182`

**Failure scenario:** measured —
`SeedRandom[1]; g=RandomGraph[{128,2000}]; Do[TimeConstrained[FindVertexColoring[g],1],{8}]`
moves `MemoryInUse[]` from 19,693,568 to 20,119,552: **~53 KB per abort, unbounded**.
`siglongjmp` unwinds past `graph_adj_free(a)`, `free(colour)`, and the search's
`work`/`best`. A script mapping `TimeConstrained[FindVertexColoring[#],5]&` over a few
thousand graphs leaks hundreds of MB.

**Why it was missed:** `test_graph_slow.c:141-144` and the plan's Risks section both
correctly identify this as "the tree's standing behaviour for every abortable builtin" and
explicitly declare the row is *not* a leak test. That framing is accurate for the few-KB
scratch buffers — but it obscures that the leaked `GraphAdj` is the largest allocation the
head makes, and that the docstring **actively instructs users to use the leaking path**
("wrap it in TimeConstrained -- that is the intended lever"). A generic tree-wide tradeoff
becomes a specific one when a docstring recommends it as the primary usage mode.

---

## MEDIUM — `fvc_search` is exported with fixed-size stack arrays but the vertex cap lives only in the caller

**Where:** `src/graph/vertexcoloring.c:275` (`char seen[FVC_MAX_VERTICES + 2]`) and `:296`
(`char forbid[FVC_MAX_VERTICES + 2]`), against the cap check at `:373`, which is in
`builtin_find_vertex_coloring` and not in `fvc_search`

**Failure scenario:** `fvc_bb` writes `seen[c]` for `c <= used + 1` and `forbid[c]` for
`c <= limit + 1`, where `used` and `limit` are bounded by `best_k <= ub <= n` — the
*graph's* vertex count, not `FVC_MAX_VERTICES`. `fvc_search` never checks `a->n`. A caller
passing a graph with `n > 130` whose DSATUR bound exceeds 129 while the clique bound is
strictly smaller (so the `lb >= ub` short-circuit at `:330` does not fire) overruns both
130-byte stack arrays. `fvc_search` is declared **non-static in `graph.h:157`** with a
header comment inviting cross-TU use, and `FVC_MAX_VERTICES` is *not* visible in
`graph.h`, so no external caller can learn the precondition. `tests/test_graph.c:301` and
`tests/test_graph_slow.c:48` already call it directly.

**Status:** latent, not demonstrated. The reviewer could not construct a triggering graph
within budget (needs χ > 129 with a strictly looser clique bound). Closed by one
`if (a->n > FVC_MAX_VERTICES) return 0;` at the top of `fvc_search`.

---

## MEDIUM — A refusal is indistinguishable from malformed input, silently degrades into arithmetic, and emits no message

**Where:** `src/graph/vertexcoloring.c:363, 366, 373, 376, 383, 387` — six distinct
`return NULL` paths with no `Message[]`; head registration at `src/graph/graph.c:172`

**Failure scenario:** measured — `Length[FindVertexColoring[CompleteGraph[129]]]` returns
**`1`**, not an error, because `Length` of the unevaluated
`FindVertexColoring[Graph[...]]` counts its one argument. So
`Length[FindVertexColoring[g]] == VertexCount[g]` is a plausible check that quietly
returns `False`, and downstream code consuming a "list of colours" gets a scalar. The same
`NULL` covers five conditions a user would act on differently: not a graph (`:366`), over
the 128-vertex cap (`:373`), budget exhausted after ~100 s (`:382`), and two OOM paths
(`:376`, `:387` — plus `fvc_search:337` and `:323`, where a `calloc` failure turns a valid,
answerable graph into a silent refusal). A user who waits 100 seconds and gets their input
back cannot tell "too hard" from "you typed the argument wrong".

**Why it was missed:** the design correctly prioritised never returning a valid-but-larger
colouring, and `NULL` is the subsystem's uniform idiom. But every other `src/graph/` head
refuses in microseconds on malformed input; this one refuses after a minute or two of real
work, which makes the missing diagnostic qualitatively different. AC-9/AC-13/AC-14 all
assert only `Head[...]`, which passes identically for all five causes.

---

## LOW — The riskiest logic in the change — the budget and the abort channel — has zero CI coverage

**Where:** `tests/CMakeLists.txt:983-987` (`EXCLUDE_FROM_ALL`, no `add_test()`)

**Failure scenario:** AC-10b, AC-10d, AC-10e and AC-10f live only in `graph_slow_tests`. A
future change that moves the `++s->steps` counter above the leaf check at
`vertexcoloring.c:249-255`, or that lets `fvc_search` return `s->best_k` on an aborted
search, is caught by nothing on any build or PR. The default suite's only budget-adjacent
assertion is `steps == 0` for `K128`/`C128`, which exercises the short-circuit, not the
budget.

**Why it was missed:** stated plainly in `test_graph_slow.c:19-22` and AC-10d, and the
reasoning (a 60 s `alarm()` in `test_utils.h`; slow tests get deleted) is sound. But the
~2 s `TimeConstrained` row (`test_timeconstrained_aborts_a_slow_instance`, out of the
file's ~2.5 min) could plausibly run on every build and currently does not.

---

## LOW — `fvc_clique_bound`'s complexity comment understates the work by ~two orders of magnitude

**Where:** `src/graph/vertexcoloring.c:190-192` — "O(n^3) worst case, which at n <= 128 is
a few million integer comparisons"

**Failure scenario:** nothing today (`FindVertexColoring[CompleteGraph[128]]` measured at
8.1 ms end-to-end). But the inner `fvc_adjacent` at `:202` is itself O(deg), making the
real cost O(n³·deg) ≈ 10⁸ operations for a dense graph at the cap. The comment is
load-bearing for whoever later considers raising `FVC_MAX_VERTICES`: it makes the bound
look free-and-linear-in-the-cap when it is quartic in the dense case, and unlike `fvc_bb`
it never polls `tc_check_deadline`, so time spent there is not interruptible on hosts where
`SIGPROF` is unreliable.

---

## Could not assess

- **Reachability of the `fvc_search` stack overflow** — no triggering graph constructed.
  Reported as an unenforced precondition on an exported symbol, not a demonstrated crash.
- **Leak-freedom of the normal (non-aborted) path** — no `valgrind` on Apple Silicon, and
  an ASan-instrumented harness would not link against the test object set within budget.
  Code reasoning shows every allocation paired and `MemoryInUse[]` was stable across normal
  calls, but that is weaker than a tool.
- **`docs/spec/builtins/` and the weekly changelog** carry no `FindVertexColoring` entry,
  which CLAUDE.md requires for every new builtin. **Not a finding** — plan.md:478-491
  records an explicit human decision moving documentation (Phase 4) to RG-5. Noted so it is
  not lost if RG-5 slips.
- **No domain checklist configured** — `.claude/GUIDANCE_ROLES.md:108` sets
  `code-review-checklist = not-configured`, so the review ran on the generic list alone.
