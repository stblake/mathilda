---
type: ladder-receipt
ticket: STATS-1
date: 2026-08-28
flow: follow-up (no new ticket; converge on STATS-1)
sha: working tree at time of run, parent ad5a6fad
baseline: ad5a6fad
supersedes: none — companion to ladder-receipt.md, which covers 0bd00b64
---

# Verification ladder receipt — STATS-1 follow-up

Scope of the change under test: `src/stats/stats_common.c`, `src/ndreduce.c`
(`ndred_quartiles` only), `tests/test_quantile_family.c`, plus the spec and
changelog pages CLAUDE.md requires for a modified builtin, and the plan's own
deviation record.

Rung commands are the four pinned in `.claude/VERIFICATION_LADDER.md`. Verdicts
are per rung, no averaging. UNAVAILABLE and NOT COMPLETE are recorded as
themselves, never as a pass.

| Rung | Command | Verdict |
|---|---|---|
| static | `python3 tools/check_c99_portability.py` | **PASS** (exit 0, no findings) |
| typecheck (= build) | `make -j8` (gcc-16, `-Werror=` set) | **PASS** (exit 0, no warnings on the changed units) |
| unit (focused) | `./quantile_family_tests`, `./stats_tests` | **PASS** — both exit 0; quantile family now 33 test functions (26 before) |
| unit (full suite) | `ctest --test-dir tests/build --output-on-failure` | **COMPLETE — 227/231 passed, 4 failed, 720.8 s, machine idle.** All four classified INHERITED below |
| integration | not-configured | **GAP** — this repo has no integration tier (unchanged) |
| judge | `adversarial-reviewer` (agent), human (absent) | **PARTIAL** — one agent pass ran against the first draft of this change and found one HIGH and five MEDIUM; all were fixed or narrowed, see below. No human judge. No second agent pass on the fixes |

## Repo audits

| Tool | Verdict |
|---|---|
| `tools/check_packed_aware.py` | **PASS** — exit 0, "every head with an NDArray fast path opts in (or is exempt with a reason)" |
| `leaks --atExit -- ./quantile_family_tests` | **PASS** — 0 leaks for 0 total leaked bytes |
| `tools/compile_coverage.py`, `nd_surface_audit.py`, `check_array_exactness.py`, `nd_fastpath_sweep.py` | **NOT RUN** this session. Stated, not smoothed: no claim is made about them |

## The four full-suite failures, classified

Same four as the `0bd00b64` receipt, with identical output signatures.

1. **image_tests** — `ImageCorrelate[..., "NormalizedCrossCorrelation"]` expects
   `{{5, 6}}`, gets `{{11, 9}}`. Byte-identical to the failure recorded at
   `0bd00b64`. → **INHERITED**.
2. **ndarray_linalg_tests** — `Det[NDArray[{{1.,2.,3.},{4.,5.,6.}}]]` returns a
   `Hold[List]`-wrapped form plus a `$RecursionLimit` message. Byte-identical to
   the failure recorded at `0bd00b64`. → **INHERITED**.
3. **plot3d_tests** — expects `9`, gets `18`. Byte-identical. → **INHERITED**.
4. **bench_pack** — the hardware-pinned performance gate. Fails on 6 workloads
   (`user f[v_] over packed` 2.66x, `Jacobi stencil` 3.05x, `Total (10^5)` 2.69x,
   `Sort (10^5)` 2.93x, `Total int64` 10.70x, `Differences int64` 5.68x). This
   gate reports the machine, so it was classified by measurement, not by
   mechanism: each of the six workloads was timed directly against the
   base-commit binary, 3 trials in each of 3 interleaved base/head rounds.

   | workload | HEAD best-of-9 | BASE best-of-9 |
   |---|---|---|
   | `Total[v]` (real, 10^5) | 0.000009 | 0.000009 |
   | `Total[vi]` (int64, 10^5) | 0.000023 | 0.000024 |
   | `Sort[v]` (real, 10^5) | 0.001742 | 0.001690 |
   | `Differences[vi]` (int64, 10^5) | 0.000025 | 0.000025 |
   | `Map[f, v]` (user f over packed) | 0.023413 | 0.023738 |
   | Jacobi-style stencil (200^2, 20) | 0.001137 | 0.001146 |

   Identical within noise on every one. → **INHERITED** — the gate's baseline
   constants do not match this machine, and it fails the same way without this
   change.

**One further failure, and how it was classified.** An earlier full-suite run of
this same tree reported **5** failures, the extra one being `primenu_tests`
(`PrimeNu[2491...1238]` expected `8`, got `7`, after
`FactorInteger::nofac: ... no factor was found within the search bounds`). That
run was made while an adversarial-review agent was running its own `ctest -j8` on
the same machine. Re-run alone on an idle machine the test **passes**, taking
8.18 s where the failing run took 3.41 s — the factorisation is effort-budgeted,
so the failure is contention, not a wrong answer. → **LOAD-INDUCED, not
introduced.** Recorded here rather than dropped, because "it passed the second
time" is exactly the kind of thing that should be written down.

## The judge rung, in detail

The adversarial pass ran against the first draft and is the reason this receipt is
not the same document it would have been an hour earlier. What it found:

- **HIGH** — the interpolation formula has a second live copy in
  `src/ndreduce.c` (`ndred_quartiles`), which still overflowed. The two surfaces
  disagreed on the same data, and the spec page edited in the same change had
  already been rewritten to describe the fixed form. Fixed in both places; both
  now pinned.
- **MEDIUM** — applying the convex form unconditionally introduced a NEW wrong
  answer (`NaN`) for weights outside `[0,1]`, which `{{a,b},{c,d}}` permits.
  Narrowed to `w ∈ [0,1]`; the historical form is kept outside it; both pinned.
- **MEDIUM** — `Complex[x, 0]` is not always normalised away (the MPFR arm
  isn't), so rejecting on the head declined a genuinely real value. Now decided
  on the imaginary part.
- **MEDIUM** — the complex fix is structural and does not see through a numeric
  head, so `Median[{1, Sqrt[2 + I], 3}]` still answers `3`. NOT fixed: narrowed
  the spec text and the code comment to claim only what the code does, and
  recorded it as a known gap in three places.
- **MEDIUM** — one new assertion was vacuous (it passed under the old code too).
  Replaced with the input that actually reproduces the bug.
- **MEDIUM** — the plan's deviation record still listed two of these as
  "Known, accepted, NOT fixed". Updated.
- Found nothing on memory ownership; independently confirmed with `leaks`.

## RECOMMENDATION (ceiling)

The strongest claim these rungs support: **the change builds clean under the
repo's `-Werror` set, passes the static portability gate and the packed-aware
audit, leaks nothing, and the full 231-test suite runs 227 passed / 4 failed on
an idle machine — where all four failures are demonstrated INHERITED, three by
byte-identical output at the base commit and one by direct base-vs-head timings
of all six flagged workloads. Both statistics suites pass. The three wrong
answers this change targets are each reproduced before and pinned after.**

Not "review-ready". One agent adversarial pass ran and its findings were fixed or
explicitly narrowed, but **no second pass ran against the fixes**, and no human
reviewed this work. The `Sqrt[2 + I]` hole is open by decision, not by oversight.
