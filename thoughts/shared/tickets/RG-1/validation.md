---
ticket: RG-1
created: 2026-08-30
source_sha: 7701df5b36cab11000efad1f3b2dff094ff2316f
subsystems: [graph]
type: validation
lifecycle: active
---

# Validation Report: `RandomGraph[{n, m}, k]` (RG-1)

**Checked:** 2026-08-30, working tree on top of HEAD `7701df5b` (the RG-1 work is
**uncommitted** — this command normally reads git history, so every judgment below is
against the worktree instead).
**Plan:** `thoughts/shared/tickets/RG-1/plan.md`
**Adversarial review:** `thoughts/shared/tickets/RG-1/adversarial.md` (HIGH resolved
23:03:49, re-verified)

## Implementation Status

- ✓ **Phase 1: The `k` form** — fully implemented.
- ✓ **Phase 2: Tests, docstring, and docs** — fully implemented; one automated box
  (full suite) remains unrun, carried over from the plan itself.

## Traceability

`trace.py` → **NOT APPLICABLE** — no `prd.md` under the tickets tree. No product
requirements to join, which says nothing about test coverage.

## Automated Verification Results

| Criterion (from the plan) | Result | Evidence |
|---|---|---|
| Build succeeds | ✓ | `make all` clean; `Mathilda` (9,191,912 bytes) newer than `generators.c` |
| `make check-c99` | ✓ | `tools/check_c99_portability.py`, no findings, exit 0 |
| Graph unit tests | ✓ | `graph_tests` rebuilt from the current tree and run: 7 tests incl. `test_random_graph`, ends `All graph tests passed!` |
| No new `-Wall -Wextra` diagnostics | ✓ | clean rebuild of `graph_tests` + `mathilda_common`, no diagnostics |
| **Full suite no worse than before** | ✗ **NOT RUN** | Excluded from this run by the human (466-binary suite). Unchanged from the plan, where this same box was already `NOT RUN`. |

## Acceptance Criteria

18 of 19 verified live against `./mathilda`; all 18 pass.

| AC | Expected | Actual | |
|---|---|---|---|
| AC-1 | `3` | `3` | ✓ |
| AC-2 | `True` | `True` | ✓ |
| AC-3 | `{6}` | `{6}` | ✓ |
| AC-4 | `{5}` | `{5}` | ✓ |
| AC-5 | `1` | `1` | ✓ |
| AC-6 | `{}` | `{}` | ✓ |
| AC-7 | `RandomGraph` | `RandomGraph` | ✓ |
| AC-8 | `RandomGraph` | `RandomGraph` | ✓ |
| AC-9 | `RandomGraph` | `RandomGraph` | ✓ |
| AC-10 | `RandomGraph` | `RandomGraph` | ✓ |
| AC-11 | `True` | `True` | ✓ |
| AC-12 | `True` | `True` | ✓ |
| AC-13 | `{Graph, 5}` | `{Graph, 5}` | ✓ |
| AC-14 | `0` | `0` | ✓ |
| AC-15 | `0` | `0` | ✓ |
| AC-16 | `{0}` | `{0}` | ✓ |
| **AC-17** | RNG stream position identical to the pre-change binary | **NOT VERIFIED** | ✗ |
| AC-18 | `RandomGraph` | `RandomGraph` | ✓ |
| AC-19 | `RandomGraph` | `RandomGraph` | ✓ |

All 19 also appear as `assert_eval_eq` rows (or the seeded `===` check) in
`tests/test_graph.c` except AC-17, which cannot be expressed in-process.

## Code Review Findings

### Matches plan

- `one_random_graph(long n, unsigned long long maxe, long m)` and `vertex_list(n)` exist
  with the planned signatures and the planned inline assembly — `make_graph` correctly
  not used, per the resolved BLOCKING finding in the plan's own review.
- The `maxe == 0` early return is keyed on `maxe` **only**, not on `m == 0` — the exact
  distinction the plan's second BLOCKING finding turned on. `RandomGraph[{5,0}]` keeps its
  `RandomSample[cand, 0]` call.
- Spec validation (`n`, `m`, `maxe`, `m <= maxe`) sits in the dispatcher, so AC-10 fails
  identically at both arities. Confirmed live.
- `k` loop frees accumulated graphs and `gs` on a `NULL` element; no partial list.
- `#include <stdint.h> /* SIZE_MAX */` added at `:24`; header comment line added at `:8`;
  `graph.h:103` signature comment updated; docstring names the `k` form and the memory
  caveat (verified via `Information[RandomGraph]`).
- Non-goals honoured: no packed/NDArray/`Compile[]` surface added (correct — `RandomGraph`
  returns a `Graph`, the "non-machine object" exemption); leak (A) untouched in the other
  four generators; no `k` cap; candidate-tree leak documented, not fixed.

### Deviations from plan — both improvements

1. **An overflow guard the plan did not specify.** `generators.c:174` adds
   `if (n > 2147483647L) return NULL;` **before** the multiply, with a 8-line comment at
   `:166-173`. The plan's Phase 1 code block computed `maxe` in `unsigned long long` and
   bounded the *product* against `SIZE_MAX` — which cannot catch a wrap. This deviation
   fixes a real out-of-bounds heap write (the HIGH in `adversarial.md`) that the plan's
   own code block would have shipped. Correct, and the stronger design.
2. **A stronger AC-18 witness in the tests.** The plan's AC-18 uses `n = 2^62`, which
   wraps to a value *above* the `SIZE_MAX` bound and so was rejected by coincidence.
   `tests/test_graph.c` adds `n = 4294967297` (`2^32 + 1`) at both arities — the witness
   that actually discriminates — and keeps a comment explaining why. Both confirmed live:
   `RandomGraph[{4294967297, 1}]` and `RandomGraph[{2147483648, 1}]` return unevaluated.

### Potential issues

1. **MEDIUM — the changelog asserts a measurement nobody made.**
   `docs/spec/changelog/2026-08-24.md:2515`: "the stream position after
   `RandomGraph[{5, 0}]` is unchanged against the pre-change binary." That is AC-17, and
   AC-17 was never run — no `Mathilda.pre` exists anywhere in the tree, and the plan's own
   manual checkbox for it is unticked. The claim is very likely true (inspection shows no
   RNG consumption added or removed on the `n >= 2` path, and the early return is keyed on
   `maxe`, not `m`), but the changelog states it as a completed comparison. Either run the
   comparison or soften the sentence to what was actually established.
2. **MEDIUM — the changelog's hardening paragraph describes the superseded design.**
   `:2525-2528` credits "`n(n-1)/2` is computed in `unsigned long long` and bounded against
   `SIZE_MAX / sizeof(Expr*)`" — the pre-fix guard — and never mentions the
   `n > 2^31 - 1` pre-multiply bound that is the load-bearing part. A reader auditing the
   overflow story is pointed at the weaker half.
3. **MEDIUM — "both `calloc`s are `NULL`-checked" undercounts the allocation sites.**
   Same paragraph. There are three on this path; `int_vertices` (`generators.c:45-49`)
   still `calloc`s and dereferences without a check. Pre-existing code, but the diff adds
   a caller (`vertex_list`, `:108`) and the changelog advertises full coverage. Carried
   from `adversarial.md`, severity unchanged.

### Resolved since the adversarial review

- **HIGH (overflow wrap → OOB heap write)** — fixed at 23:03:49; re-verified by reading
  `generators.c:174` and by running both witnesses live. Marked resolved in
  `adversarial.md`.
- **LOW (leak claim had no artifact)** — the artifact now exists; see below.

## Leak verification (Phase 1 manual criterion, now run)

Both shapes the plan's revised criterion names, via
`MallocStackLogging=1 leaks --atExit -- ./mathilda -file <script>`:

| Script | Result | Verdict |
|---|---|---|
| `Do[RandomGraph[{1, 0}], {50}]` | **0 leaks for 0 total leaked bytes** | ✓ the `maxe == 0` path is completely leak-free, as the changelog claims |
| `Do[RandomGraph[{50, 20}], {50}]` | 306,600 leaks / 18,169,600 bytes; **0** records naming `builtin_random_graph` or `int_vertices` (51 `ROOT LEAK` records total, none in this head) | ✓ leak (A) is gone for this head; the volume is 363.4 KB/call, matching Finding 4B's documented ~364 KB and attributable to the candidate tree, not a `calloc` root frame |

This is the discrimination the plan's second WORTH FLAGGING finding asked for, and it
holds.

## Manual Testing Required

1. **AC-17 — the one real gap.** Build `7701df5b` to `./Mathilda.pre` (a `git worktree` at
   that SHA avoids disturbing the uncommitted work), then on both binaries:
   - [ ] `(SeedRandom[7]; RandomGraph[{5,0}]; RandomInteger[{1, 10^6}])` — confirm the two
         integers match. This is what proves the early return did not silently extend to
         `m == 0`.
2. **Full suite.**
   - [ ] `cd tests/build && cmake .. && make -j8 && for t in *_tests; do ./$t; done` —
         confirm no regression outside `graph_tests`.
3. **Docs read-through.**
   - [x] Changelog entry is in the current ISO week's file (`2026-08-24.md`, Monday of the
         week containing 2026-08-30) — confirmed.
   - [x] Memory caveat stated in both the docstring and the changelog, with the number
         (364 KB) — confirmed.
   - [ ] `docs/spec/builtins/graphs.md` reads correctly alongside its sibling generator
         entries (human judgment).

## Recommendations

1. Fix the three changelog inaccuracies before committing — items 1–3 above. They are
   prose, not code, and two of them describe verification that did not happen, which is
   the kind of error that survives into the permanent record.
2. Run AC-17 or downgrade the claim. Do not commit the sentence as written.
3. Consider `NULL`-checking `int_vertices` in this diff after all. It is one line, it is on
   the path this ticket touches, and leaving it makes the changelog's coverage claim wrong.
4. Then run the full suite once, and commit the graph work separately from the unrelated
   staged `thoughts/shared/tickets/DEMO-2/*` deletions.

## Verdict

**Implementation matches the plan, and in two places improves on it.** All 18
in-process acceptance criteria pass in a from-scratch `graph_tests` build and live in the
REPL; the leak criterion is now backed by measurement; the adversarial HIGH is genuinely
fixed. Nothing in the code is wrong.

The gaps are in the **record**, not the implementation: three sentences in the changelog
overstate what was verified or describe a superseded guard. Plus two unrun checks — AC-17
(needs the pre-change binary) and the full suite (deliberately excluded from this run).
