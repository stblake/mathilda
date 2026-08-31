**Risk triage: 30/100 — MEDIUM.** Recommended review: one reviewer with domain familiarity
in the touched area. 24 files, 4925 reviewable lines across `origin/main...HEAD`. Drivers:
blast radius (+15, five top-level modules — `.claude`, `docs`, `src`, `tests`, `thoughts`)
and diff size (+15, past the 600-line review-quality threshold). Scored deterministically by
`pr-risk-triage`; routing signal only, it never approves.

**Read that number with one correction.** It was scored twice. Against the **shipping code
alone** — the eight files of commit 1 — it is **15/100, LOW**, 816 lines, two modules. The
jump to MEDIUM comes entirely from the three companion commits: kit configuration, a
decision record, a generated repo map, and ~3,200 lines of markdown ticket artifacts. None
of that is executable. A reviewer should size their attention to the LOW score and read the
rest as supporting material; the triage counts directories and lines and cannot tell prose
from C.

## Summary

Adds `FindVertexColoring[g]`, which returns a **minimal** vertex colouring of a graph as a
list of integers in `VertexList` order — the number of distinct colours equals the chromatic
number, and no edge joins two vertices of equal colour. Minimality is proven by exact
backtracking search seeded by a DSATUR upper bound and a greedy-clique lower bound, not
approximated by a greedy pass.

Exact colouring is NP-hard, so the head is built around one principle: it either answers
correctly or it refuses, and it never returns a valid-but-larger colouring. Two guards
enforce that, and exceeding either returns the expression **unevaluated** — graphs above 128
vertices are refused outright, and the search gives up after 8 million nodes. The budget is
a node count rather than a clock, so the answer does not depend on how fast the host is.
`TimeConstrained` is the intended lever for bounding wall time, and the search polls for the
deadline so it is honoured.

This is Phases 1–2 of the RG-2 plan. Phase 3 (packed return path) and Phase 4
(`docs/spec/builtins/` entry and changelog) are deferred to **RG-5** by recorded decision —
they change the return representation and the docs, neither of which a reviewer needs in
order to assess what was built.

## Changes

**New builtin**
- `src/graph/vertexcoloring.c` (new, 395 lines) — DSATUR upper bound (`fvc_dsatur`), greedy
  clique lower bound (`fvc_clique_bound`), and the exact branch-and-bound search
  (`fvc_search`/`fvc_bb`) with node-budget accounting and `TimeConstrained` deadline polling.
- `src/graph/graph.c` — registers `FindVertexColoring`, sets `ATTR_PROTECTED`, and installs
  a docstring covering the minimality guarantee, both guards, the `TimeConstrained` lever,
  and the fact that cost is driven by **density, not vertex count**.
- `src/graph/graph.h` — prototypes for the three helpers.
- `src/sym_names.{c,h}` — one additive interned name, `SYM_FindVertexColoring`.

**The two bounds are the design, not an optimisation.** Ascending `k = 1..ub` with only an
upper bound meant `CompleteGraph[128]` — 128 vertices, *under* the cap and therefore
accepted — had to refute k=1…127 before answering, which is a hang. The clique lower bound
plus an `lb == ub` short-circuit answers it in **zero search nodes**.

**Tests**
- `tests/test_graph.c` (+172 lines) — AC-1 … AC-18 as expression-level rows, plus
  C-level assertions on the bounds and node counter.
- `tests/test_graph_slow.c` (new, 183 lines) — the four rows that must actually spend the
  budget to prove anything.
- `tests/CMakeLists.txt` — a `graph_slow_tests` target, `EXCLUDE_FROM_ALL` with no
  `add_test()`.

## Testing

**Acceptance criteria: 18 of 18 in-scope rows pass**, each re-run directly against the built
binary rather than read from a checkbox. Full table in
`thoughts/shared/tickets/RG-2/validation.md`. The discriminating rows:

- `AC-18` — a bipartite graph colours in 2. Greedy or DSATUR on an unlucky vertex order
  returns 3, so this is the row that fails first if minimality regresses.
- `AC-8` — non-integer vertices on the path `c–a–b`, asserted position-sensitively
  (`col[[1]] === col[[3]] && col[[1]] =!= col[[2]]`). A weaker "first two differ" check
  passes under both `VertexList` order and sorted order and therefore catches nothing.
- `AC-10c` — `CompleteGraph[128]` answers 128 at **`steps == 0`**, asserted on the node
  counter rather than on elapsed time.

**Suites run**
- `./graph_tests` — all green (17 tests), including `test_vertex_coloring` and
  `test_vertex_coloring_internals`.
- `./graph_slow_tests` — all four green: sparse-at-the-cap answers at 0 nodes / 0.00 s;
  `TimeConstrained` aborts at 2.00 s; a 100-vertex dense instance **answers** at 3,872,135
  nodes / 28.2 s (the row that stops the budget being tightened until correct answers become
  refusals); a 128-vertex dense instance **refuses** at 8,000,001 nodes / 95.6 s.
- `make check-c99`, `make check-packed-aware` — pass.
- Independent differential check: minimality tested against a Python brute-force chromatic
  number over 120 random graphs (n=1..13, all densities) — **zero mismatches**.

### What the verification receipt does and does not cover

`thoughts/shared/receipts/ladder-854419997e4b+dirty.md` reports one green rung, and it is
worth reading precisely rather than as "verified":

- The ladder's `static`, `typecheck` and `integration` rungs are all **not-configured** in
  `.claude/VERIFICATION_LADDER.md` for this repo. They are coverage gaps, not passes.
- The single green rung is `unit`, and in this repo `unit` resolves to
  `make check-c99 && make check-packed-aware` — a **portability gate, not the CMake test
  suite**. A green `unit` rung does not mean the tests pass.
- `graph_tests` and `graph_slow_tests` were run **by hand**. They pass, and no receipt
  vouches for them; nothing in the ladder's current config would have caught them failing.

Both facts belong together. A receipt quoted as more than it measured is worse than no
receipt.

### Pre-existing test failures — a full ctest run, and it found more than the plan recorded

A complete `ctest` ran on this branch: **236 of 243 passed, 7 failed, 553 s**. The plan
(`plan.md:313-334`) documents three. The other four are recorded here rather than left for a
reviewer to trip over, with the honest status of each.

**The three the plan documents.** Verified directly and unpiped on this branch:

```
qrdecomposition_machine_tests -> exit 139 (SIGSEGV, no output at all)
qrdecomposition_mpfr_tests    -> exit 139 (SIGSEGV, no output at all)
plot3d_tests                  -> exit 134 (SIGABRT)
```

**The four the plan does not.** Each re-run individually on an idle machine, because the
original ctest ran concurrently with other test binaries and some of these are
timing-sensitive:

| Test | Alone | What it is |
|---|---|---|
| `primenu_tests` | **exit 0 — passes** | A ctest **flake**, not a failure. It passes in isolation; it failed only under the loaded parallel run. |
| `image_tests` | exit 134 | `ImageCorrelate` with `"NormalizedCrossCorrelation"` locates the patch at `{{11, 9}}`, expected `{{5, 6}}`. |
| `ndarray_linalg_tests` | exit 134 | `Det[NDArray[...]]` returns `Det[{Hold[List][1.0, ...], ...}]` — an `NDArray` unwrapping defect. |
| `bench_pack` | exit 1 | The packed-array performance gate: 6 workloads more than 2.5× baseline (e.g. `Differences int64` at 5.39×). Reproduces on an idle machine, so not a load artifact. |

**Provenance, stated exactly, because these three are weaker claims than the QR ones.**
`image_tests`, `ndarray_linalg_tests` and `bench_pack` are **not proven pre-existing** — no
baseline run was done for them. What can be said is that all three live in subsystems this
branch does not touch: the diff is confined to `src/graph/` plus a purely additive interned
name in `src/sym_names.{c,h}`, and nothing in it modifies `src/pack.c`'s `AWARE`/`INT64_OK`
lists, `src/ndarray.c`, or the image code. `make check-packed-aware` passes. A reviewer who
wants certainty should re-run these three at `origin/main`; that is the check this PR did
not do.

- The QR crash is **inside the QR test's own `extract_matrix`**, on its 2×2 real case —
  `EXC_BAD_ACCESS` at address `0x2`, after the 1×1 case passes. Nothing on that path touches
  `src/graph/`, and this branch's only change outside `src/graph/` is a purely additive
  interned name in `src/sym_names.{c,h}`.
- `plot3d_tests` fails a `Plot3D` mesh assertion (`Expected: 9, Actual: 18`) — unrelated to
  graph colouring by inspection of the failure itself.

**Provenance, stated exactly.** The *current* failures above were verified in this pass. The
*baseline* proof — `git stash -u`, rebuild, and `qrdecomposition_machine_tests` still exits
139 with the changes removed — was performed in an earlier session and is recorded at
`plan.md:313-334`; it was not re-derived here.

Methodological note worth carrying: the QR binaries crash with **no output whatsoever**, and
`./qrdecomposition_machine_tests | tail -3` reports exit 0, because `$?` after a pipeline is
*tail's* status. They must be run unpiped or the segfault reads as a silent pass.

### Manual verification

- [x] `FindVertexColoring[CycleGraph[5]]` → 3 colours, no adjacent pair equal
- [x] `FindVertexColoring[CompleteGraph[129]]` → unevaluated, promptly, no crash
- [x] `MATHILDA_NO_PACK=1` → identical values
- [x] `?FindVertexColoring` names the guarantee, both guards and the density cost model
- [x] Three repeated calls return identical results (search is deterministic)
- [ ] **Consider asking for this to be split** — flagged by risk triage: defect detection
      drops sharply past ~600 lines and this is 816. Counter-argument for the reviewer to
      weigh: 355 of those lines are tests, and the 395-line `vertexcoloring.c` is a single
      self-contained algorithm with no partial state worth landing on its own.
- [ ] **Identify consumers of the changed interfaces** — the triage score counts
      directories, not dependents. Relevant here: `fvc_search`, `fvc_dsatur` and
      `fvc_clique_bound` are newly exported in `graph.h`, so the answer today is "the tests
      only", and finding 2 below is precisely about what happens when that stops being true.
- [ ] **A leak verdict does not exist, deliberately.** `valgrind` has no Apple Silicon
      support, so on this host the criterion is a *wait on Linux CI*, not an unrun check.
      What was run instead: macOS `leaks` on `./graph_tests` reports 760 leaks / 198,368
      bytes with **no frame in `builtin_find_vertex_coloring`, `graph_build_adj` or any
      `fvc_*`** — the whole total attributes to the per-test harness leak every suite in the
      tree shares. Weaker than valgrind: `leaks` sees only what is reachable at exit, so a
      transient leak inside a call that later exits cleanly is invisible to it.

## Notes for reviewers

An adversarial review pass found **no HIGH findings**; three MEDIUM and two LOW are recorded
in full at `thoughts/shared/tickets/RG-2/adversarial.md`. **All five are open** — none is
fixed by this PR. The three worth a reviewer's attention:

1. **`TimeConstrained` abort leaks ~53 KB per abort.** `siglongjmp` unwinds past
   `graph_adj_free(a)`, `free(colour)` and the search's own buffers. Measured: eight aborts
   move `MemoryInUse[]` by ~426 KB. This is pre-existing tree-wide behaviour for every
   abortable builtin (`FactorInteger` does the same) and fixing it properly needs a cleanup
   registry — `docs/design/timeconstrained-abort-channel.md`. What makes it specific here is
   that **the docstring recommends that exact path** as "the intended lever". Note also that
   the plan's own Risks section (`plan.md:338-345`) describes this as "a few KB", which
   understates it by ~15× because it did not account for the `GraphAdj`; the decision to
   accept the leak stands, the recorded magnitude was wrong.
2. **`fvc_search` is exported with a precondition no caller can discover.** It writes
   `char seen[FVC_MAX_VERTICES + 2]` and `forbid[...]`, but the 128-vertex cap is enforced
   in `builtin_find_vertex_coloring`, not in `fvc_search`. The symbol is non-static in
   `graph.h:157` with a comment inviting cross-TU use, and `FVC_MAX_VERTICES` is not visible
   in `graph.h`. Latent rather than demonstrated — no triggering graph was constructed
   (it needs χ > 129 with a strictly looser clique bound). One line at the top of
   `fvc_search` closes it, and arguably belongs in this PR, since this PR is what exported
   the symbol.
3. **A refusal is indistinguishable from malformed input and emits no `Message[]`.** One
   `NULL` covers five conditions a user would act on differently: not-a-graph, over-cap,
   budget-exhausted-after-~100s, and two OOM paths. Consequently
   `Length[FindVertexColoring[CompleteGraph[129]]]` returns **`1`** — `Length` of the
   unevaluated expression — so a plausible `Length[...] == VertexCount[g]` check silently
   returns `False`. Returning `NULL` is the uniform `src/graph/` idiom, but every other head
   there refuses in microseconds on bad input; refusing after a minute of real work with no
   diagnostic is qualitatively different. AC-9/13/14/17 assert only `Head[...]`, which passes
   identically for all five causes, so the acceptance criteria cannot distinguish them either.

**Deliberate omissions, so they do not read as oversights**
- No `docs/spec/builtins/` entry and no changelog note, which `CLAUDE.md` requires for every
  new builtin. Deferred to RG-5 with Phase 4 by recorded decision (`plan.md:471-491`).
- The head returns a **plain `List`** always; `AC-19`/`AC-20` (packedness) are not satisfied.
  Also RG-5. The values are pinned by this PR's tests, which is exactly what makes the
  Phase 3 swap safe to do later — a divergence shows up as a failing row.
- The budget and abort logic — the riskiest code here — has **zero CI coverage**, because
  those rows must spend the full budget and outlast `test_utils.h`'s 60 s `alarm()`. A
  regression there is caught only by someone running `make graph_slow_tests` by hand. The
  ~2 s `TimeConstrained` row is cheap enough to promote into the default suite and probably
  should be.

**What the risk score itself does not assess.** The triage above is deterministic and reads
paths and line counts only. Its own declared blind spots, which are part of the signal
rather than boilerplate:
- **Intent and correctness** — it cannot tell whether the change is right, or a good idea.
- **Semantic weight** — a one-line change to a default can score near zero and still be
  severe. Read against this PR that cuts both ways: 816 lines drove the score, but the
  genuinely dangerous lines are the handful in `fvc_search`'s stack arrays (finding 2).
- **Whether the tests cover this change** — it reads paths, not test bodies. A renamed
  fixture counts exactly as much as a new assertion. The coverage claim in Testing above
  rests on the 18 re-run acceptance rows, not on the score.
- **Blast radius is approximate** — it counts directories, not dependencies.

**Rollback** is a file delete: the head is new and nothing else calls it.

## JIRA Ticket

RG-2
