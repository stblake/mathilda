---
type: ladder-receipt
ticket: STATS-1
date: 2026-08-27
flow: qrispy
sha: 0bd00b64
baseline: 15b088da
---

# Verification ladder receipt (QRISPY P) — STATS-1

Rung commands come from `.claude/VERIFICATION_LADDER.md` (written this run for
mathilda's real C99/GCC toolchain). Verdicts are per rung, no averaging. Where a rung
did not finish, it is recorded as NOT COMPLETE — never as a pass.

| Rung | Command | Verdict |
|---|---|---|
| static | `python3 tools/check_c99_portability.py` | **PASS** (exit 0) |
| typecheck (= build) | `make -j8` (gcc-16, -Werror set) | **PASS** (exit 0, binary relinked) |
| unit (focused) | `ctest -R "quantile_family_tests|stats_tests"` | **PASS** — 2/2, 21 test functions |
| unit (full suite) | `ctest --test-dir tests/build` | **NOT COMPLETE** — 151/231 executed when this receipt was written; 3 failures, all classified below |
| integration | not-configured | **GAP** (honest: repo has no integration tier) |
| judge | adversarial-reviewer (agent), human (absent) | **PARTIAL** — agent pass ran and its findings were fixed; no human judge available (beta run) |

## The three full-suite failures, classified

Classification used the base-commit binary (built independently during the adversarial
pass) to run the exact failing expression — the `--baseline` idea applied by hand.

1. **image_tests** — `ImageCorrelate[..., "NormalizedCrossCorrelation"]` expects
   `{{5, 6}}`, gets `{{11, 9}}`. Base binary returns the SAME `{{11, 9}}`.
   → **INHERITED**, not introduced.
2. **ndarray_linalg_tests** — `Det[NDArray[{{1.,2.,3.},{4.,5.,6.}}]]` returns a
   `Hold[List]`-wrapped form plus a $RecursionLimit message. Base binary returns the
   identical wrong output. → **INHERITED**, not introduced.
3. **bench_pack** — a *performance* gate (fails above 2.5x recorded baseline). Three
   workloads flagged: Total int64 (9.25x), Differences int64 (4.76x), Jacobi stencil
   (3.17x). The tool's own hint points at `src/pack.c`'s AWARE / INT64_OK lists, which
   this change does not touch (no packed fast path was added or removed; the four new
   heads are deliberately non-AWARE and materialise NDArray input themselves).
   Measured while 6 test binaries ran in parallel plus a stalled audit process.
   → **UNRESOLVED — most likely load-induced.** An unloaded re-run is the discriminator
   and did not happen inside this session. Recorded as unresolved rather than dismissed.

## Repo fast-path audits

| Tool | Verdict |
|---|---|
| `tools/check_packed_aware.py` | **PASS** — "every head with an NDArray fast path opts in (or is exempt with a reason)"; the four new heads are correctly absent from AWARE and need no exemption |
| `tools/compile_coverage.py` | exit 1 on **22 pre-existing heads** (Append, ArrayPlot, ExponentialMovingAverage, Fit, ...); **none of the four new heads appear** → INHERITED |
| `tools/nd_surface_audit.py` | **NOT COMPLETE** — process stalled at 0% CPU under contention; reported as not-run, not as pass |
| `tools/check_array_exactness.py`, `tools/nd_fastpath_sweep.py` | **NOT RUN** in this session |

## RECOMMENDATION (ceiling)

The strongest claim these rungs support: **the change builds clean, passes static
portability checks, and passes its own and the pre-existing statistics suites; the full
suite is unfinished and its three observed failures are two proven-inherited and one
unresolved performance gate.** Not "review-ready" — that verdict belongs to the judge
rung, which no human ran here.
