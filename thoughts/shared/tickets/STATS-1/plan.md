---
ticket: STATS-1
flow: qrispy
created: 2026-08-27
source_sha: 15b088da
subsystems: []
type: plan
lifecycle: active
status: draft
---

# Quantile & Dispersion Family (Quantile, InterquartileRange, MeanDeviation, MedianDeviation) Implementation Plan

## TL;DR
*(Human, ≤80 words)*
Add four missing Wolfram-core statistics builtins to `src/stats/`: general `Quantile`
(4-parameter form, list-of-qs), `InterquartileRange`, `MeanDeviation`, `MedianDeviation`.
The quantile engine inside `Quartiles` is generalized into a shared helper with one
deliberate semantic correction (integer-h/Ceiling case) that Quartiles' own defaults
never exercise. Risk: regressing Quartiles/Median, and the tests link trap. We'll know
it worked when the new assert_eval_eq tests and the full existing suite are green and
Quartiles' answers are unchanged.

## Overview
*(Human, ≤250 words)*
Mathilda's statistics module has exact-arithmetic Mean/Median/Variance-class builtins,
but the general `Quantile` — the repo's own declared gap (benchmarks/ABSENT.md:30) —
plus `InterquartileRange`, `MeanDeviation`, and `MedianDeviation` are absent (verified
at runtime and by source sweep; research.md).

`Quartiles` (src/stats/quartiles.c:118-176) implements Wolfram's parameterized quantile
formula — h = a+(n+b)q on Sort[data], edge clamps, interpolation weight c + d·frac —
hardcoded to q ∈ {1/4,1/2,3/4}, defaults {{1/2,0},{0,1}}. Generalizing needs one
correction: the loop interpolates A[j] → A[j+1] unconditionally, but Wolfram's
Floor/Ceiling neighbors collapse to s[[h]] at integer h. Quartiles' c=0 masks this
(weight 0); Quantile's defaults (c=1, d=0) hit it on every integer nq — uncorrected, the
engine returns s[[nq+1]], failing AC-1/3/6/7. Phase 1 extracts the loop into
`stats_common.c` with integer h returning s[[h]] — provably identical for Quartiles
(AC-12 + suite) — and adds `Quantile` with defaults {{0,0},{1,0}}:
`Quantile[{1,2,3,4},1/2]` is `2` while `Median` stays `5/2`, matching Mathematica.

Phase 2 adds the three dispersion builtins as thin compositions evaluated through the
evaluator (exact in, exact out): InterquartileRange from the Quartiles parameterization,
MeanDeviation = Mean[Abs[data − Mean[data]]], MedianDeviation = Median[Abs[data −
Median[data]]].

Registration and test wiring land IN Phase 1 (registration in stats.c, the new
`quantile_family_tests` binary, COMMON_SRC additions — the known link trap — and the
fix for a pre-existing bug found in research: `stats_tests` is built but never
`add_test`-registered, so ctest never runs it). Phase 2 extends the same files. Phase 3
adds docstrings, spec examples, the weekly changelog entry, and runs the repo's
fast-path audits and the verification ladder.

## Decisions
*(Human, ≤200 words)*
- **Quantile family over CDF completion** because Quantile is the repo's own tracked gap,
  its engine already exists in-tree, and the CDF vertical would touch the inert
  distribution-head design in src/ml/dist.c plus special functions — larger blast radius
  for a LOW-tier ticket. (assumed — beta test)
- **Extract a shared engine into stats_common.c rather than copy the loop** because two
  divergent copies of interpolation arithmetic is precisely how Quartiles/Quantile drift
  apart; stats_common.c exists for cross-TU helpers. Cost: touching proven Quartiles
  code — covered by regression AC-12 and existing test_stats.c cases.
- **Correctness path only — no AWARE[] entries, ND_REDS, or ndreduce kernels**, but
  NDArray/packed input still answered correctly: no non-AWARE stats head exists as
  precedent, and the pack gate never touches a *visible* NDArray — so each new head
  materializes such input itself via `pack_unpack` (pack.h:150-152) to the exact List
  path. Correct, slower; kernels a follow-on. Omissions declared to the audit tools.
  (assumed — beta test)
- **Fix stats_tests add_test registration in this change** because our verification
  claims "the stats suite ran"; leaving a suite that never runs would false-green V.
- **Decline symbolic/distribution arguments** (return NULL / unevaluated copy per
  existing conventions) rather than partially evaluate.

## Non-goals
*(Human, ≤150 words)*
- No distribution-argument support (`Quantile[NormalDistribution[…], q]`) — the shipped
  distributions are inert heads consumed by RandomVariate/PDF only.
- No CDF/InverseCDF/SurvivalFunction, GeometricMean/HarmonicMean, Mode, TrimmedMean,
  WeightedData — separate verticals.
- No packed-array (AWARE) fast paths, no ndreduce.c kernels, no Compile[] ND_REDS
  entries for the new heads; no INT64_OK changes anywhere.
- No Association support for the new heads (Median has it; adding it here widens scope).
- No benchmarks/ABSENT.md or benchmark-row updates — that file belongs to the numeric
  sweep tooling.
- No changes to Quartiles' observable behavior, including its message texts and its own
  {{1/2,0},{0,1}} defaults.

## Acceptance Criteria
*(Agent, no cap — dense table)*

| ID | Given | When | Then | Input | Expected |
|---|---|---|---|---|---|
| AC-1 | integer list | default Quantile at 1/2 | left-continuous type-1 result, NOT the median | `Quantile[{1, 2, 3, 4}, 1/2]` | `2` |
| AC-2 | integer list | Quantile at 1/4 | first element | `Quantile[{1, 2, 3, 4}, 1/4]` | `1` |
| AC-3 | integer list | list of qs | threads over q | `Quantile[{1, 2, 3, 4}, {1/4, 3/4}]` | `{1, 3}` |
| AC-4a | integer list | q = 1 | clamps to max | `Quantile[{1, 2, 3, 4}, 1]` | `4` |
| AC-4b | integer list | q = 0 | clamps to min | `Quantile[{1, 2, 3, 4}, 0]` | `1` |
| AC-5 | integer list | explicit Quartiles parameters | interpolated median, matches Quartiles[[2]] | `Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]` | `5/2` |
| AC-6 | real list | default Quantile | real result | `Quantile[{1., 2., 3., 4.}, 1/2]` | `2.` |
| AC-7 | unsorted list | default Quantile | sorts first | `Quantile[{3, 1, 4, 2}, 1/2]` | `2` |
| AC-8 | 1..8 | InterquartileRange | q3−q1 under Quartiles parameters | `InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]` | `4` |
| AC-9 | integer list | MeanDeviation | exact rational/integer | `MeanDeviation[{1, 2, 3, 4}]` | `1` |
| AC-10 | integer list | MedianDeviation | exact | `MedianDeviation[{1, 2, 3, 4}]` | `1` |
| AC-11a | symbolic data | Quantile | unevaluated | `Quantile[x, 1/2]` | `Quantile[x, 1/2]` |
| AC-11b | symbolic q | Quantile | unevaluated | `Quantile[{1, 2}, q]` | `Quantile[{1, 2}, q]` |
| AC-12 | regression | Quartiles after engine extraction | unchanged | `Quartiles[{1, 2, 3, 4}]` | `{3/2, 5/2, 7/2}` |
| AC-13 | matrix | Quantile columnwise | per-column result | `Quantile[{{1, 2}, {3, 4}}, 1/2]` | `{1, 2}` |
| AC-14 | mixed exact input | exact arithmetic preserved | rational out | `MeanDeviation[{1/2, 3/2}]` | `1/2` |
| AC-15 | empty list | all four heads | unevaluated, no crash | `Quantile[{}, 1/2]` | `Quantile[{}, 1/2]` |
| AC-16 | q outside [0,1] | Quantile | message + unevaluated | `Quantile[{1, 2, 3}, 2]` | `Quantile[{1, 2, 3}, 2]` |
| AC-17 | 4x3 matrix | InterquartileRange | per-column IQR, incl. the k=3 collision case | `InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]` | `{6, 6, 6}` |
| AC-18 | visible NDArray, int64 | Quantile | materialized to exact path, exact result | `Quantile[NDArray[{1, 2, 3, 4}, DataType->"int64"], 1/2]` | `2` |
| AC-19 | visible NDArray, float64 | MeanDeviation | materialized, real result | `MeanDeviation[NDArray[{1., 2., 3., 4.}, DataType->"float64"]]` | `1.` |

## Open Questions
*(Human, no cap — a question list, not prose)*

### Unresolved
_None._

### Resolved
- [x] Engine extraction vs copy — extract into stats_common.c; AC-12 + existing
  test_stats.c guard the refactor. (assumed — beta test)
- [x] Message tag for out-of-range q — use `Quantile::q100` styled after existing
  `::rectn` raw-printf convention (median.c:53); text approximates Mathematica's, exact
  wording not load-bearing. (assumed — beta test)
- [x] New sym_names entries? — not required: expr_new_symbol takes plain strings
  (stats_common.c:stats_apply_columnwise passes "Median" as a literal); pointer-compare
  reads (SYM_List etc.) all use existing symbols.

## Plan Review
*(Agent, transcribed verbatim from the plan-reviewer pass in step 3 below — never guessed)*

### Blocking
_None._

### Worth Flagging
**[WORTH FLAGGING] InterquartileRange's decline guard cannot distinguish "3-list of scalar quartiles" from "3-column matrix result"**
- Where: Phase 2, interquartilerange.c
- Why: columnwise Quartiles on a k-column matrix returns k per-column triples; k=3 slips a bare 3-element-List guard and computes Quartiles(col3)−Quartiles(col1) — silent wrong answer
- Addressed in this revision (did not need to block): matrix input now recurses via stats_apply_columnwise BEFORE the vector path; vector guard requires 3 SCALAR elements; AC-17 pins the 3-column collision case

### Resolved
**[BLOCKING] The "behavior-preserving extraction + Wolfram defaults" mechanism produces the wrong answer for the plan's own AC-1, AC-3, AC-6, and AC-7**
- Where: Overview/Phase 1 vs AC table; engine at quartiles.c:129-177
- Why: the loop interpolates A[j] → A[j+1] unconditionally; Wolfram uses Floor/Ceiling neighbors, which collapse to s[[h]] at integer h. Quantile's defaults (c=1,d=0) hit this on every integer nq (AC-1: engine would give 3, not 2); Quartiles' c=0 masks it
- Resolved: plan rewritten — the extracted engine returns s[[h]] at integer h (Ceiling upper neighbor), stated as a deliberate semantic correction, provably identical for Quartiles (AC-12 + suite)

**[BLOCKING] Phase 1 and Phase 2 success criteria depend on artifacts that don't exist until Phase 3**
- Where: old Phase 1/2 Automated Verification vs old Phase 3 (registration + all test wiring)
- Why: ctest -R would match nothing; Quantile unregistered in the REPL — every early gate unexecutable
- Resolved: phases restructured — registration, test file, COMMON_SRC additions, quantile_family_tests 3-liner, and the stats_tests add_test fix all land in Phase 1; Phase 2 extends them; old Phase 4 is now Phase 3 and carries the full-suite gate

**[BLOCKING] The no-fast-path decision rests on a precedent that is false in the source, and leaves the visible-NDArray path unaddressed**
- Where: Decisions bullet 3; Phase 4 exemption reason
- Why: Quartiles is AWARE with an ndred kernel (pack.c:553-566, ndreduce.c:550) — no non-AWARE stats head exists; and the pack gate never touches a visible NDArray, so unregistered-fast-path heads would return unevaluated on NDArray input — the repo's wrong-answer class
- Resolved: each new head materializes NDArray/packed input in-head via pack_unpack (pack.h:150-152) to the exact List path; Decisions bullet corrected (no precedent claim); audit exemption reason rewritten to be accurate; AC-18/AC-19 pin the behavior

## Requires Approval
*(Human, ≤100 words)*
_None._ — LOW tier (qualify.md 10/100 projected), no architectural impact, no new
dependencies; scope assumptions are flagged in Open Questions and research.md.

## Architecture Impact
*(Agent, no cap — fixed-shape lines)*
- New services introduced: none
- APIs changed: none (four new builtin heads; no existing head's behavior changes)
- Data crossing a service boundary: none
- New external dependency: none
- Deployment topology change: none

## Subsystems & Dependencies
*(Agent, no cap — fixed-shape lines)*
- Subsystems touched: none (no subsystem docs exist in this repo)
- Interdependencies surfaced: none

## Risks and Rollback
*(Human, ≤200 words)*
_None — standard tier, no architectural impact._ (The Quartiles refactor risk is covered
by AC-12 and the existing suite; rollback is `git revert` of a single commit.)

---

## Current State Analysis
*(Agent, no cap — findings, not prose)*
See research.md Detailed Findings. Load-bearing facts: quartiles.c:83-96 already parses
`{{a,b},{c,d}}`; :99-102 defaults {1/2,0,0,1}; :104-110 Sort via pack_eval_plain (packed
trap); :118-176 the per-q engine; median.c:52-57 message convention; stats.c registration
hub; COMMON_SRC stats block tests/CMakeLists.txt:475-490; stats_tests never add_test'd
(:1654-1655, :3216); main build wildcards src/stats/*.c (makefile:362) so only tests
need manual wiring; docstrings central in info.c:4150+.

## Desired End State
*(Human, ≤150 words)*
`Quantile`, `InterquartileRange`, `MeanDeviation`, `MedianDeviation` evaluate with exact
arithmetic and Wolfram semantics per the AC table; Quartiles/Median behavior unchanged;
all four registered ATTR_PROTECTED in stats_init with terse info.c docstrings; new
`quantile_family_tests` binary registered and green under ctest; `stats_tests` actually
runs under ctest; docs/spec/builtins/statistics.md gains an examples section per head;
weekly changelog entry exists; `make check-c99`, full build, and the audit sweep are
clean (new heads' fast-path omissions declared, not silent).

### Key Discoveries:
- quartiles.c:118-176 is the general engine — Quantile is a parameterization
- tests link via explicit COMMON_SRC list — new TUs must be added there or every test
  binary fails to link (the trap the kit's research flagged)
- stats_tests was never registered with add_test — pre-existing, fixed here

## Components & Files Affected
*(Agent, no cap — dense table)*

| File | Change |
|---|---|
| `src/stats/stats_common.h` | declare `stats_quantile_point()` engine (integer-h → s[[h]]; Ceiling upper neighbor) + shared param-parse helper if useful |
| `src/stats/stats_common.c` | new: engine extracted verbatim-in-behavior from quartiles.c:118-176 |
| `src/stats/quartiles.c:118-176` | replace inline loop with calls to the shared engine; no observable change (c=0 makes integer-h weight 0 — old and new agree) |
| `src/stats/quantile.c` | new TU: builtin_quantile (scalar q, list of qs, optional `{{a,b},{c,d}}`, matrix columnwise, NDArray materialization via pack_unpack, decline rules, `::q100` + `::rectn` messages) |
| `src/stats/interquartilerange.c` | new TU: matrix → stats_apply_columnwise recursion FIRST; vector → Quartiles, guard requires 3 SCALAR elements (stats_is_real_numeric each); NDArray via pack_unpack |
| `src/stats/meandeviation.c` | new TU: Mean[Abs[data − Mean[data]]] through the evaluator |
| `src/stats/mediandeviation.c` | new TU: Median[Abs[data − Median[data]]] through the evaluator |
| `src/stats/stats.h` | four new prototypes |
| `src/stats/stats.c` | four registrations + ATTR_PROTECTED |
| `src/info.c:~4152` | four terse docstrings next to Median/Quartiles |
| `tests/test_quantile_family.c` | new: assert_eval_eq for AC-1..AC-16 |
| `tests/CMakeLists.txt:475-490` | add 4 new TUs to COMMON_SRC stats block |
| `tests/CMakeLists.txt:~1655` | add_test for quantile_family_tests (3-line pattern) AND the missing `add_test(NAME stats_tests COMMAND stats_tests)` |
| `docs/spec/builtins/statistics.md` | reference + examples for the four heads |
| `docs/spec/changelog/2026-08-24.md` | entry for the four builtins + stats_tests fix |

## Core Flow Diagram

```mermaid
flowchart TD
    Q["Quantile[data, q(, params)]"] --> V{data a List of reals?<br/>q numeric in 0..1?}
    V -- "no: symbolic/empty" --> U[return NULL / expr_copy - unevaluated]
    V -- "bad elements / bad q" --> M["printf Head::tag message<br/>return expr_copy(res)"]
    V -- yes --> S["Sort via pack_eval_plain<br/>(quartiles.c:104 pattern)"]
    S --> E["shared engine (stats_common.c):<br/>h = a+(n+b)q; clamp; j=Floor h;<br/>A[j] + (c+d(h-j))(A[j+1]-A[j])"]
    E --> R[exact Expr result]
    QT["Quartiles[data(, params)]"] --> S2[same sort] --> E
    IQR[InterquartileRange] --> QT
```

## Alternatives Considered
*(Human, ≤150 words)*

### CDF completion for Normal/Uniform (PDF exists, CDF absent)
**Rejected because:** touches the inert-head distribution design (src/ml/dist.c strcmp
dispatch) and Erf plumbing; a bigger, less self-contained change than closing the
repo-declared Quantile gap. Second-best candidate; noted for a follow-on ticket.

### Copy the engine into quantile.c instead of extracting
**Rejected because:** two copies of interpolation arithmetic drift; extraction is ~40
lines moved with an existing regression suite covering the donor.

### Implement Quantile in C doubles for speed
**Rejected because:** the subsystem's exactness discipline (mean.c:98-108 documents a
real wrong-answer bug from approximating) — build Expr arithmetic, read doubles only for
clamping/floor, like the existing engine.

## Implementation Approach
*(Human, ≤200 words)*
Phase order minimizes risk to existing behavior: extract-and-regress first (Quartiles
must be provably unchanged before anything new lands on the engine), then the new heads,
then wiring/docs. Every builtin follows the six-step recipe in docs/extending.md:10-131.
All arithmetic goes through the evaluator (eval_and_free of SYM_Plus/SYM_Times trees);
doubles only for index math, mirroring quartiles.c. Each phase ends with the build and
the focused test binary; the full ladder runs at the end (P letter, with --baseline and
--receipt).

## Phase 1: Engine + Quantile + registration + test wiring

### Overview
Extract quartiles.c's per-q loop body into `stats_common.c` as
`Expr* stats_quantile_point(Expr** sorted_args, size_t n, Expr* q, Expr* a, Expr* b, Expr* c, Expr* d)`
(fresh Expr*; Indeterminate for non-numeric h — current behavior; **integer h returns
sorted[[h]]** — the Ceiling correction Quartiles' c=0 defaults cannot see). Rewrite
Quartiles' loop over it. Add `src/stats/quantile.c`, register Quantile in stats.h/stats.c,
create tests/test_quantile_family.c with this phase's AC rows, and do ALL CMake wiring
now (COMMON_SRC additions, quantile_family_tests 3-liner, stats_tests add_test fix) so
every later phase has a runnable gate.

### Changes Required:

#### 1. `src/stats/stats_common.{h,c}` — the engine
Signature above; body is quartiles.c:118-176's loop interior, unchanged in behavior
(same clamps h<=1 / h>=n, same j clamping, same eval_and_free arithmetic).

#### 2. `src/stats/quartiles.c` — call the engine
Replace :118-176 with three `stats_quantile_point(...)` calls (q = 1/4, 1/2, 3/4).

#### 3. `src/stats/quantile.c` — builtin_quantile
- Arity 2-3; NDArray args: none of the fast-path kernels — rely on the pack gate boxing
  (head not AWARE), so plain List arrives.
- data must be List (else expr_copy unevaluated); empty → expr_copy; matrix (first elem
  a List) → columnwise: rebuild Quantile per column carrying q and params (quartiles.c:33-65
  pattern).
- elements must pass stats_is_real_numeric else `Quantile::rectn` message + expr_copy
  (median.c:52-57 pattern).
- q: scalar real-numeric in [0,1] (inclusive; clamps make 0/1 exact min/max) → engine
  with defaults a=0,b=0,c=1,d=0; list of qs → List of per-q results; numeric outside
  [0,1] → `Quantile::q100` message + expr_copy; symbolic q → NULL-equivalent decline
  via expr_copy (stays unevaluated).
- NDArray data (`ndred_call_has_ndarray`): materialize via `pack_unpack(data)` and
  continue on the resulting List — correct exact path, no kernel (AC-18).
- optional third arg `{{a,b},{c,d}}` parsed exactly as quartiles.c:83-96; malformed →
  expr_copy unevaluated.
- Sort once via pack_eval_plain, then per-q engine calls.

### Success Criteria:

#### Automated Verification:
- [x] Build succeeds: `make -j8` (gcc-16, zero new warnings)
- [x] Test config + focused run: `cmake -S tests -B tests/build -DCMAKE_C_COMPILER=gcc-16 && cmake --build tests/build -j8 --target quantile_family_tests stats_tests && ctest --test-dir tests/build -R "quantile_family_tests|stats_tests" --output-on-failure` (AC-1..AC-7, AC-11, AC-12, AC-13, AC-15, AC-16, AC-18 written in this phase; wiring exists as of this phase)
- [x] Static analysis passes: `python3 tools/check_c99_portability.py`

#### Manual Verification:
- [x] `./Mathilda`: Quantile answers match Wolfram on the AC rows; Quartiles unchanged
- [x] Edge cases: singleton list, q=0, q=1, unsorted input
- [ ] No regressions: full existing ctest suite green

**Implementation Note**: auto mode — no human pause available; deviations recorded for
plan-deviation-tracking. (assumed — beta test)

---

## Phase 2: InterquartileRange, MeanDeviation, MedianDeviation

### Overview
Three thin TUs composing existing builtins through the evaluator; registration and
tests extend Phase 1's files.

### Changes Required:

#### 1. `src/stats/interquartilerange.c`
Arity 1. **Matrix input is detected FIRST** (first element a List) and recurses
per-column via `stats_apply_columnwise("InterquartileRange", data)` — never through the
vector path, closing the k=3 collision the reviewer flagged (a 3-column matrix's
Quartiles result is itself a 3-list and would otherwise slip the guard). NDArray input:
pack_unpack first. Vector path: evaluate `Quartiles[data]`; guard requires a 3-element
List whose elements EACH pass stats_is_real_numeric (scalars, not lists); anything else
frees and returns expr_copy(res). Then return eval of `q3 - q1` (AC-8, AC-17).

#### 2. `src/stats/meandeviation.c` / `mediandeviation.c`
Arity 1; List of real-numerics (message `::rectn` else); build
`Mean[Abs[Plus[data, Times[-1, Mean[data]]]]]` (resp. Median) as an Expr tree,
eval_and_free once; if the center (inner Mean/Median eval) fails to reduce to a
non-List numeric/symbolic-free value, decline with expr_copy. Listable arithmetic
threads data − m elementwise (verified in tests AC-9/AC-10/AC-14).

### Success Criteria:

#### Automated Verification:
- [x] Build succeeds: `make -j8`
- [x] Tests pass: `cmake --build tests/build -j8 && ctest --test-dir tests/build -R quantile_family_tests --output-on-failure` (AC-8, AC-9, AC-10, AC-14, AC-17, AC-19 added)
- [x] Static analysis passes: `python3 tools/check_c99_portability.py`

#### Manual Verification:
- [x] Wolfram cross-check of AC-8/9/10/17 values (worked by hand in research)
- [x] Empty/singleton lists decline cleanly

---

## Phase 3: Docs, audits, ladder receipt

### Overview
Close the six-step recipe (step 6) and the fast-path audit obligations; run P.

### Changes Required:

#### 1. `src/info.c` — four terse docstrings (no examples) next to :4150
#### 2. `docs/spec/builtins/statistics.md` — signatures + examples for the four heads
#### 3. `docs/spec/changelog/2026-08-24.md` — entry (file exists; no Mathilda_spec.md row needed)\n#### 3b. Full-suite gate moved here from old Phase 3: `ctest --test-dir tests/build --output-on-failure` all green; `ctest -N` lists BOTH quantile_family_tests AND stats_tests; one fault-injected AC goes red then reverts clean
#### 4. Audits: run `python3 tools/check_packed_aware.py`, `python3 tools/nd_fastpath_sweep.py`, `python3 tools/nd_surface_audit.py`, `python3 tools/compile_coverage.py`, `python3 tools/check_array_exactness.py`; if any flags the new heads, add them to that tool's declared exempt/baseline list with the reason "correctness path only — boxed by the pack gate; Quartiles precedent" (never silent).

### Success Criteria:

#### Automated Verification:
- [ ] Ladder: `python3 <kit>/skills/verification-ladder/scripts/ladder.py --baseline origin/main --receipt`
- [ ] Audit tools exit clean or with declared exemptions only

#### Manual Verification:
- [x] `?Quantile` in the REPL shows the docstring
- [x] Changelog renders coherently

---

## Testing Strategy
*(Human, ≤150 words)*
The AC table is the literal test list — test_quantile_family.c parametrizes it
one assert_eval_eq per row (spec-as-test shape; no PRD upstream, so trace.py will
report NOT APPLICABLE — the honest engineering-only outcome). Regression protection:
AC-12 plus the whole existing suite, now including the resurrected stats_tests.
Fault-injection check in Phase 3 guards against a suite that passes vacuously.

### Edge Cases & Integration Scenarios:
- Singleton list at q=0/1; unsorted input; exact rationals in, rationals out
- Matrix columnwise for Quantile AND InterquartileRange (incl. the 3-column collision
  case, AC-17); deviation heads on matrices decline or thread — pinned by a test either way
- Visible NDArray input on Quantile and MeanDeviation (AC-18/19) — materialized, never
  silently unevaluated

### Manual Testing Steps:
1. `./Mathilda` REPL: run each AC input, compare
2. `ctest -N | grep -c _tests` — count includes the two stats binaries

## Performance Considerations
*(Human, ≤100 words)*
Sort dominates: O(n log n) via existing Sort (packed-aware). The new heads are not on
AWARE[]/ndreduce, so very large packed inputs get boxed — acceptable for the
correctness-first phase and precedented (declared to the audits); fast-path kernels are
a natural follow-on.

## Migration Notes
*(Human, ≤100 words)*
None — additive builtins; no data or config migration.

## Deviations from this plan (recorded during implementation)

Captured for `plan-deviation-tracking`; every one is a change the plan did not
authorize in advance, stated rather than absorbed.

1. **Non-goal violated, deliberately.** "No changes to Quartiles' observable
   behavior" does not hold for a *mixed* exact/inexact list at integer h:
   `Quartiles[{1,2,3.,4,5,6}]` was `{2.0, 3.5, 5}`, is now `{2, 3.5, 5}`. The old
   result came from evaluating `A[j] + 0*(A[j+1]-A[j])`, where the Real neighbour
   contaminated an exact element; Wolfram's Floor/Ceiling definition never
   consults that neighbour at integer h. Kept the corrected behavior, pinned by
   `test_quartiles_mixed_exactness_at_integer_h`, and recorded in the changelog.
2. **Engine gained an exact-selection rule the plan did not describe.** At weight
   exactly 0 or 1 the engine copies the named element instead of computing the
   interpolation. Found by the adversarial pass: with default parameters
   (c=1,d=0) the arithmetic form returned `0.` for
   `Quantile[{-1.0*10^308, 2., 3., 5.}, 3/10]` and turned exact data Real for
   inexact q. HIGH severity; no acceptance criterion had covered that branch.
3. **Scope additions from the adversarial pass**, none planned: parameter-matrix
   head validation and numeric-parameter validation (symbolic params used to
   yield `Indeterminate` or a half-evaluated expression); `N[]` fallback so an
   exact irrational q answers like its `N[]` form; the non-finite element guard
   the plan's Phase 2 text required but the first implementation omitted
   (`MeanDeviation[{1, 2, Infinity}]` returned a half-evaluated expression).
4. **Known, accepted, NOT fixed** (recorded rather than silently left): raw
   `printf` messages bypass `Quiet`/`Off` (matches the existing `median.c`
   convention), and are emitted twice per declining call (evaluator re-entry);
   ragged matrix input recurses to `$RecursionLimit` (pre-existing in
   `Quartiles`, replicated); `::q100` prints twice for NDArray input;
   `compile_coverage.py` exits 1 on 22 pre-existing heads, none of them these.
   Two entries that were on this list have since been CLOSED — see deviation 5.

5. **Follow-up session (same branch, after the ladder receipt at `0bd00b64`).**
   Three of this list's own entries turned out to be reachable wrong answers
   rather than acceptable debt, so they were fixed and the non-goal below was
   violated a second time:
   - *"Complex elements pass `stats_is_real_numeric` (pre-existing helper
     flaw)"* — this is a silent wrong answer, not a cosmetic gap:
     `Median[{1, 2 + I, 3}]` returned `3` and `Quartiles[{1, 2 + I, 3, 4}]`
     returned `{2, 7/2, 3 + I/2}`, a complex quartile, from heads that had just
     printed nothing. Now decided on the imaginary part. **This violates the
     plan's "No changes to Quartiles' observable behavior" non-goal a second
     time**, deliberately, for the same reason as deviation 1: the old output
     was wrong. A `Complex[x, 0]` element, which really occurs at MPFR precision,
     is still accepted. The nested case (`Sqrt[2 + I]`) is NOT closed and stays
     on the known-gap list, in the code comment and in the shipped spec.
   - *"A third copy of the interpolation formula remains in
     `ndreduce.c:568-581`"* — accepted as duplication, but the duplicate was
     carrying a live overflow of its own, so the two surfaces disagreed:
     `Quartiles[NDArray[{-1.0*10^308, 1.0*10^308}, DataType -> "float64"]]` gave
     `inf.0` where the boxed path gave `0.0`. Both now use the convex form and
     both are pinned on the same data. The duplication itself is still there.
   - The interpolation branch of `stats_quantile_point` had the same overflow:
     `Quantile[{-1.0*10^308, 1.0*10^308}, 1/2, {{1/2,0},{0,1}}]` was `Infinity`.
     Fixed for `w` in `[0,1]` only — the adversarial pass on the fix showed that
     applying it unconditionally introduced a NEW wrong answer (`NaN`) for the
     extrapolating weights `{{a,b},{c,d}}` permits, so the historical form is
     kept outside the unit interval and both cases are pinned.
   Also in this session, and NOT in the plan: a leaf fast path in
   `stats_is_real_numeric` (1.15x-2.85x on 200k-element input, measured
   base-vs-head over three interleaved rounds; table in the changelog), and
   seven new test functions covering the above plus three coverage gaps the
   phase-1 set did not reach.

## References
- Original ticket: (mission brief; no tracker)
- Related research: `thoughts/shared/tickets/STATS-1/research.md`
- Similar implementation: `src/stats/quartiles.c:118-176`, `src/stats/median.c:15-70`
- QRISPY receipts: `thoughts/shared/tickets/STATS-1/qualify.md` (Q), trace.md (S, at verify)
