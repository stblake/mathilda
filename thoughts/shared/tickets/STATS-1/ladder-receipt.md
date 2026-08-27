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
| unit (full suite) | `ctest --test-dir tests/build` | **COMPLETE — 227/231 passed (98%), 4 failed, exit 8.** All four classified INHERITED below. `stats_tests` and `quantile_family_tests` both PASSED in the full run |
| integration | not-configured | **GAP** (honest: repo has no integration tier) |
| judge | adversarial-reviewer (agent), human (absent) | **PARTIAL** — agent pass ran and its findings were fixed; no human judge available (beta run) |

## The four full-suite failures, classified

Classification used the base-commit binary (built independently during the adversarial
pass) to run the exact failing expression — the `--baseline` idea applied by hand.

1. **image_tests** — `ImageCorrelate[..., "NormalizedCrossCorrelation"]` expects
   `{{5, 6}}`, gets `{{11, 9}}`. Base binary returns the SAME `{{11, 9}}`.
   → **INHERITED**, not introduced.
2. **ndarray_linalg_tests** — `Det[NDArray[{{1.,2.,3.},{4.,5.,6.}}]]` returns a
   `Hold[List]`-wrapped form plus a $RecursionLimit message. Base binary returns the
   identical wrong output. → **INHERITED**, not introduced.
3. **plot3d_tests** — `Length[Plot3D[x+y, ..., PlotPoints->4, MaxRecursion->0,
   Mesh->None][[1]]]` expects `9`, gets `18`. Base binary returns the identical `18`.
   → **INHERITED**, not introduced.
4. **bench_pack** — a *performance* gate (fails above 2.5x its recorded baseline).
   RESOLVED after the suite finished and the machine went idle. The load hypothesis was
   WRONG: unloaded it fails harder (6 workloads, incl. Total int64 11.40x, Sort 2.87x,
   Jacobi 3.12x), so contention was not the cause. The discriminator that settles it is
   a direct base-vs-head timing of the flagged workloads, 3 trials each:

   | workload | HEAD | BASE (15b088da) |
   |---|---|---|
   | `Total[Range[10^5]]` | 0.000167 / 0.000104 / 0.0001 | 0.000105 / 0.000116 / 0.0001 |
   | `Sort[Table[Mod[k*7919,10^5],{k,10^5}]]` | 0.021417 / 0.021323 / 0.021685 | 0.021974 / 0.020890 / 0.021197 |
   | `Differences[Range[10^5]]` | 0.000167 / 0.000187 / 0.000148 | 0.000173 / 0.000200 / 0.000146 |

   Identical within noise on every one. This change touches no packed fast path
   (`src/pack.c` untouched; the four new heads are deliberately non-AWARE), and the
   evidence agrees. → **INHERITED** — the gate's baseline constants do not match this
   machine, and it fails the same way without this change.

## Repo fast-path audits

| Tool | Verdict |
|---|---|
| `tools/check_packed_aware.py` | **PASS** — "every head with an NDArray fast path opts in (or is exempt with a reason)"; the four new heads are correctly absent from AWARE and need no exemption |
| `tools/compile_coverage.py` | exit 1 on **22 pre-existing heads** (Append, ArrayPlot, ExponentialMovingAverage, Fit, ...); **none of the four new heads appear** → INHERITED |
| `tools/nd_surface_audit.py` | **NOT COMPLETE** — process stalled at 0% CPU under contention; reported as not-run, not as pass |
| `tools/check_array_exactness.py`, `tools/nd_fastpath_sweep.py` | **NOT RUN** in this session |

## RECOMMENDATION (ceiling)

The strongest claim these rungs support: **the change builds clean, passes the static
portability gate, and the full 231-test suite runs at 227 passed / 4 failed — where all
four failures are demonstrated INHERITED, three by identical wrong output at the base
commit and one by identical base-vs-head timings. Both statistics suites pass, one of
which had never executed under ctest before this change.**

Not "review-ready" — that verdict belongs to the judge rung. An agent adversarial pass
ran and its findings were fixed; no human judge reviewed this work.

Outstanding coverage gaps, stated rather than smoothed over: no integration tier exists
in this repo; `nd_surface_audit.py` never completed (stalled, then killed);
`check_array_exactness.py` and `nd_fastpath_sweep.py` were not run.
