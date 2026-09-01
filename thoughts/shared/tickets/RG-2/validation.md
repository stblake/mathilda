---
ticket: RG-2
created: 2026-08-31
validator: /ais:check-against-plan (chained from /ais:verify-implementation)
plan: thoughts/shared/tickets/RG-2/plan.md
source_sha: 854419997e4b0d8b35fb401783a17b13db2495c9 + uncommitted working tree
branch: find-vertex-coloring
scope_checked: Phases 1-2 only (Phases 3-4 descoped to RG-5 per plan.md:471-491)
---

# Validation Report: RG-2 — FindVertexColoring

## Implementation Status

| Phase | Status |
|---|---|
| Phase 1: The algorithm, unregistered | Fully implemented — all automated criteria ticked, one manual criterion deliberately unticked (see below) |
| Phase 2: Registration | Fully implemented — all automated and manual criteria ticked |
| Phase 3: The packed producer path | **Not attempted — out of scope.** Descoped to RG-5 by human decision, plan.md:471-491 |
| Phase 4: Documentation | **Not attempted — out of scope.** Descoped to RG-5, same decision |

The scope split is a recorded human decision, not an omission. It is the reason AC-19,
AC-20 and AC-21 are unsatisfied and why `docs/spec/builtins/` and the weekly changelog
carry no `FindVertexColoring` entry despite CLAUDE.md requiring one for every new builtin.

## Automated Verification Results

- **PASS** — `make check-c99`
- **PASS** — `make check-packed-aware`
- **PASS** — full build, `SDKROOT` exported, no new `-Wall -Wextra` diagnostics
- **PASS** — `./graph_tests`, all tests green, including `test_vertex_coloring` and
  `test_vertex_coloring_internals`
- **PASS** — `./graph_slow_tests`, all four rows green:
  - `test_sparse_at_the_cap_is_fast` — `CycleGraph[128]` chi=2 and
    `RandomGraph[{128,200}]` chi=3, both at **0 search nodes, 0.00 s**
  - `test_timeconstrained_aborts_a_slow_instance` — `$Aborted` at 2.00 s
  - `test_budget_does_not_refuse_a_solvable_instance` (AC-10f) — chi=8, 3,872,135 nodes,
    28.2 s, inside budget
  - `test_budget_refuses_rather_than_guessing` (AC-10d) — returns 0 at 8,000,001 nodes,
    95.6 s
- **NOT RUN** — ladder `static`, `typecheck`, `integration` rungs: all declared
  `not-configured` in `.claude/VERIFICATION_LADDER.md`. The only automated rung that ran is
  `unit`, and for this repo that is `check-c99 && check-packed-aware` — a **portability
  gate, not the CMake suite**. The CMake suites above were run by hand and are not in the
  ladder's config, so nothing in the ladder would have caught them failing.
- **PASS** — traceability join (`spec-as-test/scripts/trace.py`): `NOT APPLICABLE`, no
  `prd.md` under the tickets tree. Engineering-only ticket; says nothing about test coverage.

### Acceptance criteria — re-run independently against the built binary

Not read from the plan's checkboxes; each expression evaluated fresh in `./Mathilda`.

| AC | Expected | Actual | Verdict |
|---|---|---|---|
| AC-1 | `5` | `5` | PASS |
| AC-2 | `2` | `2` | PASS |
| AC-3 | `3` | `3` | PASS |
| AC-4 | `2` | `2` | PASS |
| AC-5 | `{1}` | `{1}` | PASS |
| AC-6 | `7` | `7` | PASS |
| AC-7 | `True` | `True` | PASS |
| AC-8 | `{3, True}` | `{3, True}` | PASS — the position-sensitive form, so `VertexList` order is genuinely pinned |
| AC-9 | `FindVertexColoring` | `FindVertexColoring` | PASS |
| AC-10 | `128` | `128` | PASS |
| AC-10b | unevaluated, bounded | unevaluated at 8,000,001 nodes / 95.6 s | PASS (slow target) |
| AC-10c | `128`, `steps == 0` | `128`, `steps == 0` | PASS (default suite) |
| AC-10d | `fvc_search` returns 0 | returns 0 | PASS (slow target) |
| AC-10e | `$Aborted` | `$Aborted` at 2.00 s | PASS (slow target) |
| AC-10f | answered inside budget | chi=8, 3,872,135 nodes | PASS (slow target) |
| AC-11 | `{1}` | `{1}` | PASS |
| AC-12 | `{}` | `{}` | PASS |
| AC-13 | `FindVertexColoring` | `FindVertexColoring` | PASS |
| AC-14 | `FindVertexColoring` | `FindVertexColoring` | PASS |
| AC-15 | `2` | `2` | PASS |
| AC-16 | `2` | `2` | PASS |
| AC-17 | `FindVertexColoring` | `FindVertexColoring` | PASS |
| AC-18 | `2` | `2` | PASS — minimality beats greedy, the discriminating row |
| AC-19 | `is_packed_list(r)` true | head returns a plain `List` always | **NOT SATISFIED — out of scope, RG-5** |
| AC-20 | `!is_packed_list(r)` + `{1,2,3}` | value correct, packedness not implemented | **NOT SATISFIED — out of scope, RG-5** |
| AC-21 | `3` under `MATHILDA_NO_PACK=1` | `3` | PASS (checked despite being deferred — passes trivially, since nothing is packed yet) |

**18 of 18 in-scope rows pass. 2 of 3 deferred rows unsatisfied by design.**

Also verified: `Options[FindVertexColoring]` is `{}`; the docstring names the minimality
guarantee, both guards, the `TimeConstrained` lever and the density-not-count cost model;
three repeated calls on `CycleGraph[7]` return identical results (search is deterministic).

## Code Review Findings

### Matches plan

- Adjacency is the undirected neighbourhood; AC-15 confirms direction is ignored.
- Ownership follows the subsystem norm: `res` is never freed, `graph_build_adj` borrows,
  all six early returns return `NULL` without freeing. Independently audited by the
  adversarial reviewer.
- Both bounds present as specified: `fvc_clique_bound` (greedy clique lower) and DSATUR
  (upper), with the `lb == ub` short-circuit — which is what makes `CompleteGraph[128]`
  answer at zero search nodes, the plan's own stated reason for adding the lower bound.
- The three helpers are non-`static` with prototypes in `graph.h`, and registration stayed
  in Phase 2 — the "Phase 1 registers nothing" decision is intact.
- Slow rows live in `tests/test_graph_slow.c` under `EXCLUDE_FROM_ALL` with no
  `add_test()`, exactly as specified.

### Deviations from plan

1. **The abort leak is ~15x larger than the plan's Risks section says.** plan.md:338-345
   describes "the search's scratch buffers (a few KB, bounded by the 128-vertex cap)".
   Measured: **~53 KB per abort**, because `siglongjmp` also unwinds past
   `graph_adj_free(a)` — the `GraphAdj` is the largest allocation the head makes, and the
   plan's wording does not account for it. The *decision* to accept the leak is unchanged
   and still defensible; the recorded magnitude is wrong and should be corrected before
   this text is reused as the justification.
2. **`fvc_search` is exported with a precondition no caller can discover.** It writes
   `char seen[FVC_MAX_VERTICES + 2]` and `forbid[...]` sized to the cap, but the cap is
   enforced in `builtin_find_vertex_coloring`, not in `fvc_search`. The symbol is non-static
   in `graph.h:157` with a comment inviting cross-TU use, and `FVC_MAX_VERTICES` is not
   visible in `graph.h`. Latent, not demonstrated — the reviewer could not construct a
   triggering graph. The plan specified the export; it did not consider that exporting moves
   the invariant across the module's API boundary.

### Potential issues (from the adversarial pass — full detail in `adversarial.md`)

- MEDIUM: the ~53 KB-per-abort leak on the very path the docstring recommends as "the
  intended lever".
- MEDIUM: unenforced vertex-count precondition on the exported `fvc_search`.
- MEDIUM: a refusal is indistinguishable from malformed input and emits no `Message[]` —
  one `NULL` covers five conditions (not-a-graph, over-cap, budget-exhausted, two OOM
  paths). `Length[FindVertexColoring[CompleteGraph[129]]]` returns **`1`**, so a plausible
  `Length[...] == VertexCount[g]` check silently returns `False`. AC-9/13/14/17 assert only
  `Head[...]`, which passes identically for all five causes — so the plan's own criteria
  cannot distinguish them either.
- LOW: the budget and abort logic — the riskiest code here — has zero CI coverage.
- LOW: `fvc_clique_bound`'s complexity comment understates real cost by ~2 orders of
  magnitude (O(n³·deg), not O(n³)), which matters to whoever later raises the cap.

None is HIGH; none blocks. All are recorded rather than fixed, per this ticket's scope.

## Pre-existing failures — not caused by this change

plan.md:313-334 records three `ctest` failures — `qrdecomposition_machine_tests` and
`qrdecomposition_mpfr_tests` (SEGFAULT) and `plot3d_tests` (Subprocess aborted) — and
proves them pre-existing by re-running at baseline with `git stash -u`, where
`qrdecomposition_machine_tests` still exits 139. The QR crash is inside the QR test's own
`extract_matrix`, on its 2×2 real case, `EXC_BAD_ACCESS` at `0x2`.

**Provenance note:** that baseline stash-and-rebuild was performed by an earlier session
and is being read from the plan, not re-derived here. This validation did not re-run it. A
full `ctest` was started during this pass and had not finished when the report was written.

Methodological note worth carrying (plan.md:331-334): the QR binaries crash with **no
output at all**, and `./qrdecomposition_machine_tests | tail -3` reports exit 0 because
`$?` after a pipeline is *tail's* status. They must be run unpiped or the segfault reads as
a silent pass.

## Manual Testing Required

Already performed this pass:
- [x] `FindVertexColoring[CycleGraph[5]]` → 3 colours, no adjacent pair equal
- [x] `FindVertexColoring[CompleteGraph[129]]` → unevaluated, promptly, no crash
- [x] `MATHILDA_NO_PACK=1` → identical values
- [x] `?FindVertexColoring` docstring content

Outstanding, and not closable on this host:
- [ ] **A leak verdict.** Deliberately unticked in the plan (Phase 1) and still unticked.
      `valgrind` has no Apple Silicon support, so this is a WAIT on Linux CI, not an
      unrun check. macOS `leaks` reports 760 leaks / 198,368 bytes on `./graph_tests` with
      **no frame in `builtin_find_vertex_coloring`, `graph_build_adj` or any `fvc_*`** —
      the total attributes to the per-test harness leak every suite in the tree shares.
      Weaker than valgrind: `leaks` sees only what is reachable at exit, so a transient
      leak inside a call that later exits cleanly is invisible to it.

## Recommendations

1. **Correct the "a few KB" figure** in plan.md:338-345 to the measured ~53 KB and name
   `GraphAdj` as the dominant leaked allocation. Cheap, and it stops the understatement
   being inherited by RG-5 or by `docs/design/timeconstrained-abort-channel.md`.
2. **Add `if (a->n > FVC_MAX_VERTICES) return 0;` at the top of `fvc_search`.** One line,
   closes an unenforced precondition on an exported symbol. Arguably in scope for RG-2
   since RG-2 is what exported it.
3. **File RG-5 after submit**, carrying Phases 3-4, AC-19/AC-20, the `docs/spec/builtins/`
   entry and the changelog note — plus the three MEDIUM findings above if they are not
   fixed here.
4. **Consider promoting the ~2 s `TimeConstrained` row into the default suite.** It is the
   one slow-target row cheap enough to run on every build, and it covers the abort channel,
   which currently has no CI coverage at all.
5. **Wire the ladder's `integration` rung to the CMake suite** in
   `.claude/VERIFICATION_LADDER.md`. Right now `graph_tests` passing is a fact about what
   someone typed by hand, not something any receipt can vouch for.

## Verdict

**MATCHES PLAN**, for the scope RG-2 actually claims. Every in-scope acceptance row was
re-verified independently against the built binary rather than read from a checkbox. The
two deviations are documentation-accuracy and a latent precondition, not behavioural
divergence: the head does what the plan says it does, including refusing rather than
guessing.
