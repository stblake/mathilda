# Weekly gap-driven benchmark job: Mathilda vs Python vs Mathematica

Branch: `work/2026-08-03`. Plan: `~/.claude/plans/imperative-scribbling-ember.md`.

## Goal

A repeatable weekly job built from kept `.m`/`.py` file pairs that runs 31
experiments in three systems, joins the rows by label, and emits a ranked report
naming the week's dev work — with *slower* and *absent/incomplete* kept apart so a
missing algorithm is never averaged into a speed number, plus a coverage % of how
many probed heads Mathilda actually has.

## Findings that shaped the build (research, before any code)

- [x] The checked-in `Mathilda` binary was **stale and silently emitted nothing**
      (`-v` printed nothing, exit 0). Rebuilt to `0.024`.
- [x] Hardcoded tool paths in the existing harnesses are wrong on this host:
      `HPC_PYTHON=/usr/local/bin/python3.11` absent (python3 is 3.14.3),
      `wolframscript` at `/usr/local/bin`. `numba 0.65.0` **is** installed, which
      contradicts a standing caveat in `docs/experiments/README.md`.
- [x] `Get` resolves against **cwd** (`$InputFileName`/`DirectoryName` absent), and
      `Get["../harness.m"]` works identically in Mathilda and `wolframscript` —
      so ONE shared prelude, not 60 copies.
- [x] Mathilda exits 0 after errors; stdout is line-buffered on redirect
      (`repl.c:777`). So: classify from parsed rows + stderr, never `$?`; partial
      rows survive a timeout kill.
- [x] **`RandomReal[{}, dims]` returns unevaluated in Mathilda; Mathematica
      supports it.** All 15 WolframMark tests use that spelling — a verbatim port
      would report fictitious ~0.1 ms wins on 8 tests (Fourier, Eigenvalues, Dot,
      Transpose, SVD, LinearSolve, matrix arithmetic, elementary functions).
- [x] Absent heads found so far: `FindFit`, `AccuracyGoal`, `DSolve`, `Reduce`,
      `SparseArray`, `MatrixExp`, `EllipticK`, `Refine`.

## Build

- [x] 1. `benchmarks/harness.m` + `harness.py` — the ONLY copy of bench/check/require
- [x] 2. `benchmarks/data.m` + `data.py` — shared seeded input generators
- [x] 3. `benchmarks/run_all.py` — discover, run, parse, join, classify, rank, render
- [x] 4. `makefile` target `bench-gap` (+ `.PHONY`)
- [x] 5. Group D: `29-graph-ops`, `30-string-ops` — subsystems nothing had timed
      (a WolframMark group was built, then removed: it is a *hardware*
      benchmark, not a CAS one — it earned its keep by exposing the
      `RandomReal[{}, dims]` incompatibility and was retired)
- [x] 6. Group A: 01–10 symbolic vs sympy
- [x] 7. Group B: 11–20 numeric libraries vs scipy
- [x] 8. Group C: 21–28 the array substrate (the eight open roadmap items)
- [x] 9. Run the full job — 202 rows, 19.7 min three-system (8.0 min without
      Mathematica), inside the ~15 min budget for the two-system case
- [x] 10. `benchmarks/README.md` + generated `REPORT.md` / `ABSENT.md`
- [x] 11. Coverage %: 87.4% (160/183 declared heads present)
- [x] 12. `SPEC.md` companion-docs pointer + `docs/spec/changelog/2026-08-03.md`

## Harness errors found and fixed before publishing (each would have been a false finding)

- [x] `bench` 2-arity never dispatched: `HoldRest` held `$BenchReps`, so
      `reps_Integer` could not match. Rows emitted no timing at all.
- [x] `$BenchReps` as a **Module local** made `Table[..., {reps}]` return
      unevaluated (a Mathilda bug in its own right — see the changelog). Rows read
      `0.001 Round[1e+06 Table[0.0, {reps}]]`.
- [x] `data.m` used `Product[...]`, which legitimately returns an **unexpanded
      closed form**; `Exponent` then read 30 instead of 20 and `PolynomialGCD`
      bailed. Looked exactly like a `PolynomialGCD` bug and was not one.
- [x] Fit/interpolation data written `i/100` is an **exact Rational** in Wolfram,
      so Mathilda did exact rational linear algebra against numpy float64.
      Reported 239000× and 70000×; both artifacts. Now `N[]`-pinned.
- [x] `x^5+x+1` antiderivative carries an unresolved root sum that renders
      differently per system, breaking the check (not the timing) → `x^4+1`.
- [x] Convention errors in the Python column, each caught because Mathilda and
      Mathematica **agreed** and numpy/sympy was the outlier: `ArcTan[x,y]` =
      `atan2(y,x)`, `ListConvolve` kernel reversal, Wolfram `Fourier`'s
      `e^{+2πi}` sign and `1/√n` scaling, `Quartiles` quantile definition,
      numpy views vs materialised copies (`.T`, slices, `reshape`).

## Review

**Shipped.** `benchmarks/` — 31 experiments as kept `.m`/`.py` pairs, one shared
`harness.{m,py}` + `data.{m,py}` (not 62 copies), `run_all.py`, `make bench-gap`.
186 cases across four areas, coverage 87.4%, 20.2 min for all three systems /
8.0 min without Mathematica. Outputs: `REPORT.md`, `ABSENT.md`, `history.jsonl`
(source of truth) + `HISTORY.md` (rendered view), `results/<date>.json`.

**Monitoring.** Time-weighted progress bar with an ETA calibrated from the
previous run's per-experiment durations — count-based extrapolation was useless
here (it said 58s at experiment 7 of a 20-minute run) because case costs span
three orders of magnitude. `-v` for per-system lines, `--from-json` to re-render
without re-measuring, `--only`/`--system` write `.partial` files and never touch
the canonical weekly artefacts.

**Findings, in value order.**
1. `Mean`/`Variance` int64 overflow — **fixed this pass**, with a regression test.
2. `Fit` falls off a fast path at 4 basis terms: 3.2 ms → 1035 ms between 3 and 4
   terms at 500 points; 23.16 s vs Mathematica's 0.078 ms at degree 5.
3. Graph accessors are O(n²) and do not cache — 10× vertices → 98.8× time, and
   `EdgeCount` costs the same 5.4 s as `VertexDegree`, so `Graph[…]` is rescanned
   from its edge list on every call. Mathematica: 0.000 ms.
4. `Table[e, {k}]` returns unevaluated when `k` is a `Module` local.
5. `SeriesCoefficient[Series[…], n]` returns unevaluated.
6. `RandomReal[{}, dims]` returns unevaluated (Mathematica accepts it).
7. Build trap: FFTW/FLINT/GMP-ECM absent → `Fourier` was O(n²), 1050.9 ms → 1.0 ms
   at 32768 after installing. `make` alone does not relink; `make clean` required.

**By area** (median ratio, vs Mathematica / vs Python): A symbolic 1.34× / 0.10×,
B numeric libraries 0.44× / 0.98×, C array substrate 2.36× / 1.19×, D uncovered
subsystems 0.66× / 322×. So ~10× ahead of sympy, at parity with scipy and ahead of
Mathematica on the numeric libraries, ~2.4× behind Mathematica on the array
substrate — where the roadmap already pointed — and on the untimed subsystems
*ahead* of Mathematica while hundreds of times behind Python. Splitting graph and
string out of C is what made that last split visible; mixed in, both areas read as
a flat 2.1×/2.3× that said nothing.

**Harness errors caught before publishing** (each would have been a false
finding): `HoldRest` swallowing `$BenchReps`; `Table` with a Module-local count;
`Product[…]` returning an unexpanded closed form (looked exactly like a
`PolynomialGCD` bug); `i/100` being an exact Rational so Mathilda did exact
rational least squares against numpy float64 (reported 239000×, an artifact);
`wolframscript` rendering `ToString[1/2]` as a 3-line 2-D fraction so the parser
read `1`; and several Python-side convention errors, each caught because Mathilda
and Mathematica **agreed** and the Python column was the outlier.

**Not done, deliberately.** `tests/test_stats.c` has 3 pre-existing failures in
`test_central_moment` (print-form assertions that do not hold in the test binary,
which does not load the `.m` bootstrap). Verified pre-existing by stashing this
branch's `src/` changes and reproducing them. Left for a separate change rather
than repaired inside a benchmarking PR.
