# Worked Exemplar: RPI Loop on a Real Ticket

**What this is.** A full Research → Plan → adversarial Plan-Review → Implement → Verify
pass on a real ticket in a real codebase (Mathilda, ~365 kLoC C99), using the `ais`
(AI-SDLC Starter Kit) plugin, `v8.0.0`. Every artifact quoted below is the actual file, the
actual sub-agent output, or the actual REPL/build output produced during that pass — not a
reconstruction written after the fact to look tidy. Where something below *is* reconstructed
from the session's edit history rather than pulled from a separate committed snapshot, that
is stated explicitly at the point it happens, not smoothed over.

**Why this ticket.** Two other candidate tickets existed in the same session. This one was
picked because its adversarial plan-review pass caught two independently real, verified
defects — an inexact-arithmetic bug and a wrong test-update claim — which is a more
instructive demonstration of what the review step is actually for than a single catch would
be.

**Source repo state**: `mathilda`, commit range `bbcd9bde..52303d9a` on branch `main`
(all commits below are real, in that range, unmodified for this document).

---

## 1. The ticket, before any research

There was no externally-authored ticket file — this ticket was **self-generated** by the
same session, from evidence a prior, already-shipped ticket (edge weights /
`WeightedAdjacencyMatrix`) left behind in its own `## Non-goals` section:

> "No weighted-shortest-path / Dijkstra mode for `FindShortestPath` or `GraphDistance`
> (research doc's resolved Open Question — deferred to a follow-up; real algorithmic scope
> growth, not a few hours)."
> — `thoughts/shared/plans/2026-08-22-graph-edge-weights.md`

The prompt that started this specific ticket, verbatim, from the operator running the
session:

> "take a second, harder ticket in this repo, chosen the same way as the first — from what
> the code and tests show a maintainer would want, sized at a few hours. Run the full loop
> again and append findings to the same journal, numbering continuously."

No ticket description existed beyond that — the research phase's job was to turn "the
obvious next thing a prior ticket's own Non-goals named" into an actual scoped ticket.

---

## 2. Research artifact, as written

Full text of `thoughts/shared/research/2026-08-23-weighted-shortest-path.md`, unedited:

````markdown
---
created: 2026-08-23T04:15:30Z
researcher: Michael Sollami
source_sha: 3d872247
branch: main
repository: mathilda
topic: "Second graph extension: what does the code/tests show a maintainer would want next?"
tags: [research, codebase, graph, shortest-path, dijkstra]
subsystems: [graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-23
last_updated_by: Michael Sollami
---

# Research: Second graph extension — weighted shortest path

**Date**: 2026-08-23T04:15:30Z
**Researcher**: Michael Sollami
**Git Commit**: 3d872247
**Branch**: main
**Repository**: mathilda

## TL;DR
Ticket 1 (edge weights) explicitly deferred making `FindShortestPath`/`GraphDistance`
weight-aware, naming it in its own `Non-goals` as the natural next step. `GraphAdj`
(`graph_util.c`, shared by 8 builtins) stores **no per-edge weight at all** — only
successor/predecessor vertex indices — so real Wolfram Language semantics (both builtins
auto-dispatch to a weighted algorithm when `EdgeWeight` is present) require a real Dijkstra
implementation, not a config flag. This is harder than ticket 1: a new algorithm, not just a
new builtin, and it changes two existing builtins' behavior on weighted graphs rather than
adding new read-only ones. Sized at a few hours given a simple O(V²) array-based Dijkstra
(matching this codebase's own precedent: `VertexConnectivity`'s docstring calls itself
"exact brute-force ... intended for small graphs").

## Summary
`src/graph/shortestpath.c` implements unweighted BFS for both `FindShortestPath[g,s,t]` and
`GraphDistance[g,s,t]`, routed through the shared `GraphAdj` (`graph_build_adj`). Real
Wolfram Language's own `FindShortestPath`/`GraphDistance` automatically use edge weights when
present and fall back to unweighted BFS otherwise — that is the behavior ticket 1's own
research and plan explicitly named as deferred. `GraphAdj` has no weight storage, so this
requires building a small, local weighted-adjacency pass (reusing `graph_resolve_edge_weights`
from ticket 1) rather than touching the shared structure 8 other builtins depend on —
learning directly from ticket 1's plan-reviewer-caught lesson about `graph_build_adj` being a
sensitive shared choke point.

## Open Questions

### Unresolved
_None._

### Resolved
- [x] Does `GraphAdj` already carry weights that a Dijkstra pass could reuse? — No, confirmed
      by reading `graph_util.c`'s `GraphAdj` struct and `graph_build_adj`'s fill loop: only
      `int` successor/predecessor indices, no weight field anywhere.
  - [x] Should this touch the shared `GraphAdj`/`graph_build_adj`? — No: build a local,
      call-scoped weighted adjacency inside `shortestpath.c` instead, to avoid widening the
      blast radius of a structure 8 builtins depend on (direct lesson from ticket 1's
      `plan-reviewer` finding).
- [x] How should non-numeric or negative weights be handled? — Fall back to the existing
      unweighted BFS behavior rather than failing: Dijkstra requires non-negative numeric
      weights to be correct, and a previously-working call should not start returning
      unevaluated just because a graph happens to carry a symbolic or negative weight.
      Documented as an explicit limitation, matching this codebase's existing style
      (`VertexConnectivity`'s own "intended for small graphs" self-limitation).
- [x] Prior attempt or known constraint? — None found in git history; this session's own
      ticket-1 Non-goals is the only prior signal, and it points at doing exactly this.
````

*(Full document continues with Detailed Findings, Code References, Architecture Insights —
see the file itself; truncated here for length. Note the "8 builtins" figure above — this
was carried into the plan and turned out to be wrong; see §4 below.)*

**Cost of this stage**: no research sub-agent was dispatched — the author already held
working context on `src/graph/` from the prior ticket and read `shortestpath.c` directly.
This is a real, load-bearing shortcut, not a null cost: it means the "research" step's
value here is smaller than usual (no independent verification pass), and the exact
figure it carried forward unverified (see §4) is arguably a direct consequence of skipping
that independent check.

---

## 3. Plan, as first written (before adversarial review)

The plan was written directly to `thoughts/shared/plans/2026-08-23-weighted-shortest-path.md`
in one pass. Two sections as they stood **before** `plan-reviewer` read the file — quoted
verbatim from the session's own edit history (the exact `old_string` of the fix applied
later), since no separate git commit captured this intermediate state on its own:

**Overview, as first written:**
> "`GraphAdj` (`graph_util.c`), the structure both builtins currently use, has no weight
> storage and is shared by **6 other builtins** — so this plan builds a separate, call-scoped
> weighted adjacency inside `shortestpath.c` rather than widening `GraphAdj` itself..."

**Phase 1 §2, "Dijkstra", as first written:**
> "**File**: `src/graph/shortestpath.c`
> **Changes**: a static `dijkstra()` mirroring `bfs()`'s signature/shape (fills `parent[]`,
> and a `double dist[]` this time since weights may be non-integer), O(V²) array scan for the
> minimum-unvisited-distance vertex each iteration (no heap, matching `VertexConnectivity`'s
> existing complexity precedent)."

**Components & Files Affected, as first written (excerpt):**
> "Update `test_edge_weights`'s `FindShortestPath`/`GraphDistance` AC-11 lines to the new
> weighted behavior"

The `grill-me` section-contract check
(`skills/grill-me/scripts/check_plan_contract.py`) was run against this draft and failed
once on a word-count overage in `## Decisions` (204/200 words) — fixed by tightening prose,
unrelated to the substantive findings below, and not itself evidence of anything wrong with
the plan's content.

**Cost of this stage**: one `Write` call for the full plan (~310 lines), two contract-check
runs (instant, scripted), no interactive back-and-forth with the operator — the plan-open
`grill-me` question pass was skipped for this ticket on the operator's own explicit
instruction to proceed without waiting. Wall-clock: git commit timestamps bracket
research+plan-writing at **12m 53s** (`3d872247` at `21:14:01` to `bbcd9bde` at
`21:26:54`). Stated plainly rather than smoothed over: this bracket is **not a clean
measurement of this stage alone** — it also contains reading `shortestpath.c` fully and, in
this specific run, an unrelated side investigation (independently verifying an upstream
version-drift claim from a peer session, logged elsewhere as this session's GR-15) that
happened to fall in the same window. The true research+plan-writing time is smaller than
12m 53s; this document reports the honest, checkable bracket rather than an unverifiable
finer split.

---

## 4. Adversarial plan-review — verbatim reviewer output

A fresh sub-agent, briefed with the `ais` kit's own `agents/plan-reviewer.md` rubric (nine
areas, two adversarial lenses — for a *plan* artifact, scope-boundary + testability),
dispatched with instructions to read every file the plan cited directly rather than trust
its citations. Output below, unedited in substance (only HTML-entity-escaping in the
original source restored to plain characters, and the two `[WORTH FLAGGING]` findings'
`Where`/`Lens` metadata lines dropped for length — every `Why it doesn't hold` and `What
would resolve it` sentence is the reviewer's own unaltered wording, not a paraphrase):

> ## Plan review
>
> **Reviewed:** `thoughts/shared/plans/2026-08-23-weighted-shortest-path.md`
> **Rubric areas with findings:** Untestable steps, Acceptance criteria have external ground
> truth, Hidden assumptions
>
> ### Findings
>
> **[BLOCKING] Dijkstra's `double dist[]` has no stated path back to an exact `Expr`, so
> AC-2 (and AC-1/AC-6's implicit numeric correctness) is not achievable as specified**
> - Where: Phase 1 §2 "Dijkstra" ("a `double dist[]` this time since weights may be
>   non-integer") vs. Acceptance Criteria AC-2 (`GraphDistance[...] -> 3`)
> - Lens: both lenses agree — testability (the AC's literal expected value `3`, not `3.`,
>   has no guaranteed path from the stated design) and scope-boundary (a hidden, undecided
>   departure from this codebase's exact-arithmetic convention)
> - Why it doesn't hold: `src/print.c:259` documents that `EXPR_REAL` prints distinctly from
>   `EXPR_INTEGER` (e.g. `3. ` vs `3`), and `tests/test_utils.h`'s `assert_eval_eq` does
>   exact string comparison. The plan's only stated design is a raw `double` accumulator
>   with no mentioned conversion step; a direct `expr_new_real(dist[it])` implementation
>   would print `3.` and fail AC-2's literal `3`. Nowhere in Decisions, Open Questions, or
>   Phase 1 is "how the accumulated distance becomes an `Expr`" ever decided.
> - What would resolve it: an explicit decision — e.g. accumulate distances as `Expr*` via
>   the existing exact-arithmetic add/compare helpers, or accumulate as `double` and add a
>   stated post-pass that converts a whole-valued result back to `EXPR_INTEGER`/
>   `EXPR_RATIONAL` before falling back to `EXPR_REAL`. Either way this belongs in Phase 1
>   §2, not left implicit.
>
> **[BLOCKING] The claim that "test_edge_weights' AC-11 lines for these two builtins must
> change" is false for `FindShortestPath` — only `GraphDistance`'s assertion needs to
> change**
> - Where: Decisions ("Ticket 1's `test_edge_weights` AC-11 lines for these two builtins
>   must change") and Components & Files Affected vs. `tests/test_graph.c:354-359`
> - Lens: testability (acceptance criteria / ground truth)
> - Why it doesn't hold: the AC-11 test graph is `Graph[{1,2,3},{1->2,2->3},
>   EdgeWeight->{5,7}]` — a directed path with exactly one route from vertex 1 to vertex 3
>   (`1->2->3`). Dijkstra and BFS necessarily agree on the *path* here (there is no
>   alternative route to disagree over); only the *distance* changes, from the BFS hop-count
>   `2` to the weighted total `12` (5+7). `FindShortestPath`'s existing assertion `{1, 2, 3}`
>   is already correct under the new behavior and does not need to change.
> - What would resolve it: correct the claim to name only `GraphDistance`'s AC-11 assertion
>   (`2` → `12`) as requiring a value change.
>
> **[WORTH FLAGGING] `graph_weights_usable`'s numeric-type list omits `EXPR_MPFR`, a
> first-class numeric leaf type in the default build, and reinvents rather than reuses the
> codebase's own generic numeric-type helper**
> - Why it doesn't hold: `src/expr.h` guards `EXPR_MPFR` behind `#ifdef USE_MPFR`, but the
>   makefile defaults `USE_MPFR ?= 1`, so it is a live leaf type in the ordinary build.
>   `src/expr.c:412-436` already defines `expr_is_numeric_like()` — the exact "existing
>   numeric-value helper" the plan gestures at — which explicitly includes `EXPR_MPFR`. A
>   weight built from a high-precision real would be judged "not usable" and silently fall
>   back to unweighted BFS.
> - What would resolve it: name `expr_is_numeric_like` (minus its `Complex` branch) as the
>   actual reused helper.
>
> **[WORTH FLAGGING] The recurring "6 other builtins" count is wrong (and internally
> inconsistent with the plan's own file list) — the real number is 5**
> - Why it doesn't hold: verified directly against source and against `tests/test_graph.c`'s
>   own AC-11 block, which exercises exactly 5 other `graph_build_adj`-routed builtins:
>   `ConnectedComponents`, `WeaklyConnectedComponents`, `FindSpanningTree`,
>   `ConnectedGraphQ`, `VertexConnectivity`. `grep` of `graph_build_adj` call sites also
>   shows 6 call sites (not 7 — `shortestpath.c` has a single shared call site inside
>   `resolve()`, used by both builtins) across 4 files (not 5). This repeated miscount is
>   inherited from the (already-implemented) prior ticket's plan, which made the same error.
> - What would resolve it: correct the count to "5 other builtins" everywhere it appears.
>
> ### Could not assess
> - None — all claims checked against live source (`shortestpath.c`, `graph_util.c`,
>   `graph.h`, `expr.h`, `expr.c`, `print.c`, `tests/test_graph.c`,
>   `docs/spec/builtins/graphs.md`) rather than taken on the plan's or research doc's
>   citations.

**Cost of this stage**: one sub-agent dispatch, **10 minutes 27 seconds** of the agent's own
processing time (self-reported by the agent runtime: `duration_ms: 626367`), 32 tool calls,
~147.8k tokens — this is the one number in this document that is a direct tool-reported
figure rather than something computed from commit timestamps, and it is worth flagging that
it does **not** cleanly reconcile with the commit-timestamp gap either side of it (see §8) —
stated here rather than quietly resolved, since a document about honest cost accounting
should not paper over its own inconsistent numbers. This was, at the time, the single most
expensive step in the entire ticket by both of the numbers available for it — genuinely felt
like the most likely candidate to cut if this were being done under time pressure, and would
have been the wrong cut: see below.

---

## 5. What changed in the plan because of these findings

Both BLOCKING findings were fixed before the plan was approved — not deferred, not
overruled. The actual diffs applied (verbatim `new_string` from the session's edit
history):

**Fix 1 — exact-value reconstruction, added to Phase 1 §2:**
> "...using a `double dist[]` **for internal vertex-selection comparisons only** — O(V²)
> array scan for the minimum-unvisited-distance vertex each iteration... This resolves a
> real gap a `plan-reviewer` pass caught in the previous draft: a raw `double` accumulator
> returned directly as `GraphDistance`'s result would print as `EXPR_REAL` (e.g. `12.`)...
> and this codebase treats exact arithmetic as load-bearing throughout. Fix: once
> `dijkstra()` finds the parent chain to `t`, reconstruct the **exact** total by evaluating
> `Plus[w1, ..., wk]` (via `evaluate()`) over the actual `Expr*` weights... giving
> `GraphDistance` an exact `EXPR_INTEGER`/`Rational` result whenever the inputs are exact...
> `FindShortestPath` needs no such reconstruction — it returns the vertex path, not a
> distance value."

**Fix 2 — corrected test-update claim:**
> "Update `test_edge_weights`'s `GraphDistance` AC-11 line only (`"2"` → `"12"`...) —
> `FindShortestPath`'s AC-11 assertion (`{1, 2, 3}`) is unaffected, since that specific test
> graph has only one path from vertex 1 to vertex 3, so BFS and Dijkstra necessarily agree
> on it."

Both findings, and their resolutions, were transcribed into the plan's own `## Plan Review`
section under `### Resolved` (the kit's own convention — a Blocking finding that gets fixed
moves there with a one-line note, rather than being deleted or left sitting under
`### Blocking`, which would still gate `/implement-plan`). The two WORTH FLAGGING findings
(the `EXPR_MPFR` omission, the builtin-count miscount) were fixed the same way, non-blocking
but not ignored either.

**What would have shipped without this review pass, concretely:**
- `GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}],1,4]` would
  have returned `3.` (an `EXPR_REAL`) instead of the exact `3` (`EXPR_INTEGER`) a Wolfram
  Language user would expect from an all-integer-weighted graph — a visible, wrong-type
  result, not just an internal inefficiency.
- The test suite would likely have shipped with an **incorrect** `FindShortestPath`
  assertion "corrected" to some other value where none was needed, on a test graph where the
  correct answer already matched — a follow-the-plan-literally implementer had no reason to
  independently re-derive that the original assertion was already right.

---

## 6. Implementation

Real diff, `git diff bbcd9bde..81bcb7a6 --stat -- src/`:
```
 src/graph/graph.h        |  14 ++
 src/graph/graph_util.c   |  50 ++++++
 src/graph/shortestpath.c | 231 ++++++++++++++++++++++++++++++++++++++------
 3 files changed, 271 insertions(+), 24 deletions(-)
```

Key additions: `graph_weights_usable(g)` and `graph_weight_to_double(w)` in
`graph_util.c` (reusing `expr_is_numeric_like`, per the reviewer's fix); a local `WAdj`
weighted-adjacency struct, `dijkstra()`, and `exact_path_weight()` in `shortestpath.c`
(the exact-value reconstruction from §5); both builtins branch once, at the top, on
`graph_weights_usable(g)`.

Build, exactly as run:
```
$ export SDKROOT=$(xcrun --show-sdk-path)   # this machine's toolchain needed this; see §9
$ make -j$(sysctl -n hw.ncpu)
... EXIT: 0, zero warnings, zero errors
```

**Cost of this stage**: commit timestamps put writing + building + manually verifying the
implementation between `fd666669` (`21:32:27`) and `81bcb7a6` (`21:40:13`), a **7m 46s**
gap — the smallest-caveat number in this document, since no sub-agent or other interleaved
work is known to have fallen inside this specific window.

---

## 7. Acceptance criteria — each one, and how it was actually checked

Every row below was run against the live built `./Mathilda` binary via a `-file` script,
and the actual printed output is quoted, not paraphrased.

| ID | Criterion | Command run | Output |
|---|---|---|---|
| AC-1 | Min-weight path, not min-hop path | `FindShortestPath[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}],1,4]` | `{1, 2, 3, 4}` |
| AC-2 | Min total weight, exact value | `GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4},EdgeWeight->{1,1,1,10}],1,4]` | `3` |
| AC-2 (type) | Result is exact, not `Real` | `Head[GraphDistance[...]]` | `Integer` |
| AC-3 | Unweighted graphs unaffected | `FindShortestPath[CycleGraph[6],1,4]` | `{1, 2, 3, 4}` |
| AC-4 | Symbolic weight falls back to BFS | `FindShortestPath[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{a,7}],1,3]` | `{1, 2, 3}` |
| AC-5 | Negative weight falls back to BFS | `GraphDistance[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{-1,7}],1,3]` | `2` |
| AC-6 | Undirected weighted graph, symmetric | `FindShortestPath[Graph[{1,2,3},{1<->2,2<->3},EdgeWeight->{1,1}],1,3]` | `{1, 2, 3}` |
| AC-7 | Unreachable target | `FindShortestPath[Graph[{1,2,3},{1->2},EdgeWeight->{5}],1,3]` / `GraphDistance[...]` | `{}` / `Infinity` |
| (extra) | Rational weights stay exact | `GraphDistance[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{1/2,1/3}],1,3]` | `5/6` |
| (regression) | Ticket 1's AC-11, corrected | `GraphDistance[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}],1,3]` | `12` (was asserted `2` before this ticket) |

Automated: `make check-c99` and `make check-packed-aware` both exit 0, no new findings.
`tests/build-main/graph_tests` — 17 tests including the new `test_weighted_shortest_path` —
all pass. The repo's own verification-ladder `unit` rung reported FAILED in this run, for a
reason unrelated to this change (a pre-existing, alphabetically-earlier, unrelated flaky
optimization test halts the `for t in *_tests` loop before `graph_tests` ever runs) —
confirmed by running `graph_tests` standalone inside the ladder's own build directory, where
it passes cleanly. Reported here rather than omitted, per the same "don't let a green
summary hide a real gap" principle this whole exercise is about.

**Cost of this stage**: commit timestamps put verification + journal writeup at **3m 54s**
(`81bcb7a6` to `52303d9a`) — the fastest gap of the four, because every check here is either
scripted (`make check-*`) or a short REPL script, with no further design work.

---

## 8. Total cost, honestly — including where the numbers don't add up cleanly

Two independent measurements exist, and they do not fully reconcile — reported here as
found, not smoothed into a single tidy figure:

- **Total wall-clock span across the whole ticket**, by git commit timestamps: `3d872247`
  (`21:14:01`, the commit immediately before this ticket's own work started) to `52303d9a`
  (`21:44:07`, this ticket's verification+writeup commit) = **30m 6s**, across 4 commits.
- **The one direct, tool-reported duration available**: the plan-review sub-agent's own
  `10m 27s` (`duration_ms: 626367`).

| Stage | What the commit history shows | What made it worth it, or not |
|---|---|---|
| Research + plan writing | `3d872247` → `bbcd9bde`, 12m 53s gap — but this window also contains an unrelated side investigation (this session's GR-15), so the true research+plan time is smaller than 12m 53s, not a clean measurement | No sub-agent dispatched — a real shortcut, and the likely source of the one uncaught factual error (the "8 builtins" count) that survived into the plan |
| Adversarial plan-review | `bbcd9bde` → `fd666669`, **5m 33s** gap — this is the number that does not reconcile: the review agent itself reports `10m 27s` of its own processing, longer than the gap between the two commits either side of it. Not resolved further here; reported as an open inconsistency rather than picking whichever number looks better | Caught two real, ship-affecting defects regardless of which duration number is right |
| Implementation | `fd666669` → `81bcb7a6`, 7m 46s — the cleanest bracket in this table, no known interleaved work | Straightforward once the plan's design decisions were actually settled — review moved cost *earlier*, where it's cheaper to fix, not away |
| Verification + writeup | `81bcb7a6` → `52303d9a`, 3m 54s | Fast because every check was scripted or a short REPL script |

Turn counts are not reported here at all: no exact log of tool-call counts per phase was
kept during the run, and giving a number would imply more precision than actually exists.
The honest summary: this ticket took on the order of **half an hour of wall-clock and four
commits**, with a genuine, unresolved discrepancy in exactly how that time split across the
research/plan and review stages — included rather than hidden, because a document arguing
for the value of catching untrustworthy claims should not itself ship one.

**Where this felt like overhead at the time, stated plainly**: the plan-review dispatch is a
10-and-a-half-minute pause with no visible progress from the operator's seat, immediately
after a plan that already looked complete, well-cited, and ready. That is exactly the moment
a real user under time pressure is most likely to skip it — and exactly the moment where, in
this run, skipping it would have shipped a wrong-typed `GraphDistance` result and a
needlessly "corrected" test assertion.

---

## 9. What this exemplar does not hide

**GR-01 touched every phase of this ticket.** No RPI command in this session (`/research-codebase`,
`/create-plan`, `/implement-plan`, `/verify-implementation`) was ever invocable via the
harness's own Skill-tool dispatch — the plugin was installed and reported enabled
(`ais@ais`, `8.0.0`) throughout, but calling it by name returned `Unknown skill`, and a
freshly spawned sub-agent confirmed the same skills were invisible to it too. Every command
this document describes as "run" was instead read from the plugin's cache directory as a
markdown file and followed by hand, including running its referenced Python scripts
directly via `Bash`. This was a real cost in this session, not a hypothetical one — it is
also, per a later report the operator has not independently re-verified, one that (a) did
not reproduce across three fresh cold installs on clean repos, and (b) traces to a
specific, narrower root cause (a lost marketplace registration entry, not a general defect
in the kit's commands themselves) combined with `/reload-plugins` having no
agent-invocable form. Both halves of that sentence matter for reading this exemplar
honestly: the failure was real and had a real cost in this specific run, and it is reported
as not general and not reproduced, rather than either smoothed over or overstated into "the
kit doesn't work."

**The "8 builtins" miscount** in §2's research doc was never independently re-checked by a
second pass in this ticket — it happened to get caught anyway, because the *plan*
(reviewed) restated it and the reviewer verified every claim against live source rather than
trusting the plan's or the research doc's own citations. Had the plan not restated the
count, or had the reviewer trusted the citation instead of re-deriving it, this specific
factual error would likely still be sitting, uncorrected, in a shipped document.
