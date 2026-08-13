# `benchmarks/` — the weekly gap-driven benchmark job

Thirty experiments, each a kept `.m`/`.py` pair, run in three systems and
joined into a ranked report that names the week's dev work.

The unit of measurement is a **case**: one named kernel, timed in every system
and joined across them by its label. There are 186 of them.

## How to run it

```bash
make bench-gap                     # the whole thing, 3 systems, ~20 min
```

That is the only command you need weekly. It prints a live progress bar with a
calibrated ETA and writes four files. While it runs:

```
30 experiments x mathilda+wolfram+python  (timeout 240s/file)
last run's total for these: 20m12s — that is the ETA basis

  ██████████████░░░░░░░░░░░░  54%  19/30  7m40s    ETA 5m32s   statistics
```

The bar is weighted by **time**, not by experiment count — costs here span three
orders of magnitude, so counting experiments would show 87% while the two
four-minute timeout experiments were still ahead. It replays the previous run's
per-experiment durations, rescaled by how the current run is tracking. The first
run on a fresh checkout has no basis and says so instead of inventing a number.

Useful variations:

```bash
python3 benchmarks/run_all.py --list                  # what exists
python3 benchmarks/run_all.py --only 19,29            # a subset (writes .partial)
python3 benchmarks/run_all.py --system mathilda,python   # skip Mathematica, ~8 min
python3 benchmarks/run_all.py --timeout 600           # per-FILE seconds
python3 benchmarks/run_all.py -v                      # per-system lines, no bar
python3 benchmarks/run_all.py --from-json results/2026-08-04.json   # re-render only
python3 benchmarks/run_all.py --check-labels          # fail on .m/.py label drift
```

`--only` and `--system` produce a **subset**, so they write `REPORT.partial.md` /
`ABSENT.partial.md` / `<date>-partial.json` and leave the canonical weekly files
untouched. They are also excluded from ETA calibration and from `HISTORY.md`. This
matters: the first time a 3-experiment smoke test was run it overwrote the full
suite's report and the timing basis.

`--from-json` re-renders the report from a saved run without re-measuring — use it
after changing report wording or a threshold, rather than paying 20 minutes again.

### Outputs

| file | what it holds |
|---|---|
| `REPORT.md` | trend, ranked gaps, per-area medians, coverage %, build warnings |
| `ABSENT.md` | the absence register — what Mathilda does not have |
| `history.jsonl` | **the source of truth for the arc** — one JSON object per full run, append-only |
| `HISTORY.md` | a rendered *view* of `history.jsonl`. Never edit it; it is regenerated |
| `results/<YYYY-MM-DD>.json` | every raw case, plus the host and version block |

A second full run on the same day **archives** the first as
`results/<date>.1.json` rather than overwriting it — otherwise re-running to
check a fix would destroy the before-picture the trend needs. The trend then
compares against that archive, so same-day iteration works.

`history.jsonl` is data, not prose, because the arc is the thing you will want to
query and plot — not re-parse out of markdown. One object per line means it diffs
cleanly in git and reads without a parser:

```bash
# is group C closing on NumPy?
python3 -c "import json
for l in open('benchmarks/history.jsonl'):
    r = json.loads(l)
    g = r['groups']['C numeric roadmap (numpy)']
    print(r['date'], r['commit'], g['median_vs_python'], g['slower'])"
```

Each record carries the git commit, so a row ties to the code that produced it,
and per-area medians against *each* baseline separately — so "is the symbolic
half closing on Mathematica" is one series, not something to reconstruct from old
reports.

### Reading the report

Top to bottom, it answers four questions in order:

1. **Trend since the last full run** — coverage and class-count deltas, then
   *Regressed* before *Improved*. Regressions are listed first deliberately: a
   benchmark job that only reports wins is a marketing document.
2. **Ranked gaps** — every `SLOWER` case by ratio against the best other system.
   This is the work queue.
3. **Where the slowness is, by area** — median ratio per group against
   Mathematica *and* against Python, side by side. Medians, not means: one
   600000× case would make a mean meaningless.
4. **Absences and check failures** — capability gaps and unverifiable cases,
   never mixed into the speed numbers.

The two median columns are the ones to read against each other. Behind
**Mathematica** is a gap against a *competitor* — same evaluation semantics, same
automatic compilation — so it is an engineering gap. Behind **Python** is a gap
against *the machine*: on the dense linear-algebra areas all three call the same
Accelerate BLAS, so any spread is pure overhead. A case can read acceptably
against one and badly against the other, and that difference is usually the
finding.

A ratio whose baseline ran under 0.005 ms is rendered `>N×` with a footnote, not
as a precise number: 0.001 ms is at clock resolution, so `6097090.0×` would imply
seven significant figures the measurement cannot support. The gap is real; its
magnitude is not precise.

---

## Run this at the end of each week

```bash
make -n 2>&1 | grep "not detected"   # 1. must print NOTHING — see The build gate
make clean && make -j8               #    clean is required: CFLAGS changes do
                                     #    not invalidate .o files
make bench-gap                       # 2. ~20 min, live progress bar + ETA
git diff benchmarks/REPORT.md        # 3. the diff IS the week's progress
```

Then read, in this order:

1. **`REPORT.md` → Trend** — what regressed since last week, before what improved.
2. **`REPORT.md` → Ranked gaps** — the work queue, worst first.
3. **`REPORT.md` → by area** — which subsystem the week should go into.
4. **`ABSENT.md`** — feature requests, not performance bugs. Never a ratio.
5. **`HISTORY.md`** — one line per run. If the arc is flat over a month, the
   weekly reports have been read and not acted on.

Commit all four generated files. The trend section and `HISTORY.md` both work by
diffing against the previous committed run, so skipping the commit costs you the
comparison.

---

## Why this exists along`docs/experiments/`

[`docs/experiments/`](../docs/experiments/) holds twenty experiments as
**narrative**: each one investigated a hypothesis, and the write-ups are the
record. Its `run.sh` prints the three systems' raw text side by side — it does
not parse, join, rank, or classify.
[`comparisons/hpc_bench.py`](../comparisons/hpc_bench.py) does aggregate, but its
kernels are inline in a 1791-line Python file, so there is nothing on disk to
edit or run standalone.

This directory is the **job**: file pairs on disk, a machine-readable contract,
and a ranked report that can be diffed week over week.

It also covers ground the existing corpus does not. All twenty existing
experiments are numeric, against NumPy. Group A here is the first time the
symbolic half of the system — `src/calculus/`, `src/poly/`, `src/simp/`,
`src/sum/`, among the largest subsystems in the tree — has had **any** external
baseline.

---

## The five case classes, and why absence is never a ratio

| class | meaning | drives |
|---|---|---|
| `SLOWER` | every system answered; Mathilda ≥ 1.5× the best other | kernel work — **has a ratio** |
| `AHEAD` | Mathilda fastest, or within 1.5× | nothing; regression guard |
| `ABSENT` | the head does not exist in Mathilda | feature work — **never a ratio** |
| `INCOMPLETE` | head exists but returned unevaluated, errored, or timed out | correctness work |
| `CHECK-FAIL` | the systems disagree on the check value | **timing discarded** |

`SLOWER` and `ABSENT` cases are reported in separate files on purpose. Averaging "we do
not have `DSolve`" into a speed number turns a missing algorithm into a fake 40×,
and sends a week of work at the wrong file.

A `CHECK-FAIL` case prints **no timing at all**. An unverified timing is worse than
a missing one, because the columns may be timing different programs. That is not
hypothetical — see *The RandomReal bug* below.

The 1.5× threshold, not 1.0×, because run-to-run spread on large-working-set cases
is documented at ±20% in `docs/experiments/README.md`; 1.5 sits outside that band.

---

## The build gate — read this before trusting any case

A missing **optional** dependency does not disable a feature in this tree. It
silently substitutes a slower algorithm. The makefile says so and then the
warning scrolls past in a `-j8` build:

```
makefile:295: FFTW not detected; building with USE_FFTW=0
              (Fourier uses a naive O(n^2) fallback)
makefile:256: FLINT >= 3.0 not detected; building with USE_FLINT=0
              (algebraic-extension GCD/Factor use the classical fallback)
makefile:129: GMP-ECM not detected; building with USE_ECM=0
              (advanced factorisation disabled)
```

Measured on this tree, with and without FFTW, same source, same host:

| `Fourier` size | no FFTW | with FFTW |
|---|---|---|
| 8192 | 46.8 ms | — |
| 16384 | 234.4 ms | — |
| 32768 | **1050.9 ms** | **1.0 ms** |

The doubling ratios of ~4.5–5× are the signature of an O(n²) fallback. Reported
without the build block, that case reads *"Mathilda is 1000× slower at FFT"* and
sends a week of work at `src/fourier.c` — when the fix is `brew install fftw`.

`Mathilda -v` cannot tell you this: it lists what **is** compiled in (`GMP`,
`MPFR`, `Raylib`, `Accelerate`, `Readline`) and says nothing about what is
missing. So `run_all.py` re-runs the makefile's own detection on every run and
puts the result at the top of `REPORT.md`.

```bash
brew install fftw flint gmp-ecm          # macOS
sudo apt install libfftw3-dev libflint-dev libecm-dev   # Debian/Ubuntu
make clean && make -j8                   # clean is required
```

**`make` alone is not enough.** Changing `CFLAGS` does not invalidate existing
`.o` files, so a plain rebuild after installing a dependency relinks nothing and
the fallback stays in the binary. This cost a full diagnosis pass to notice.

---

## Layout, and the one shared copy of everything

```
benchmarks/
  harness.m / harness.py    THE ONLY copy of bench / check / require
  data.m    / data.py       shared deterministic input generators
  run_all.py                the job
  NN-slug/slug.m            kernels only, 2 lines of preamble
  NN-slug/slug.py           kernels only, 2 lines of preamble
```

The existing corpus duplicates its reporting helpers into every file ("shared
reporting helpers (identical in every experiment file)"). With 30 experiments
that would be 60 copies to keep in step, so this directory shares them instead:

- **`.m` side**: `Get["../harness.m"]`. Mathilda's `Get` resolves against **cwd**
  (`$InputFileName` and `DirectoryName` do not exist), so `run_all.py` sets cwd to
  the experiment folder. Verified to load identically under `wolframscript`.
- **`.py` side**: a `sys.path` insert off `__file__`, so it is cwd-independent.

`data.m` / `data.py` is the second thing worth sharing: the value gate only means
something if both halves compute over the *same* inputs, and one copy per
language makes that true by construction rather than by two authors agreeing on a
formula.

### Running one experiment by hand

```bash
cd benchmarks/01-integrate-rational       # cwd matters for the .m side
../../Mathilda -file integrate_rational.m
wolframscript  -file integrate_rational.m
python3 integrate_rational.py             # works from anywhere
```

### The emit contract

One tagged tab-separated line per case, identical in all three systems:

```
BENCH   <label>  <milliseconds>
CHECK   <label>  <value>
REQUIRE <head>   present|absent
SKIP    <label>  <head>
```

`<label>` is the join key between the two halves. `run_all.py --check-labels`
fails an experiment whose halves emit different label sets, because silent
case-drift is how a table starts comparing two different programs.

Timing method: one untimed warm-up, then the **minimum** of 3 runs
(`AbsoluteTiming` / `perf_counter`, never `Timing[]`, which sums CPU time across
threads and overstates every threaded path by roughly the core count). Cases that
are already expensive on their own use `benchOnce` — one timed run, no warm-up.

---

## The RandomReal bug — why the value gate is not optional

`RandomReal[{}, dims]` — an empty range, which Mathematica reads as the default
`{0, 1}` — **returns unevaluated in Mathilda**.

This was found by porting Wolfram's WolframMark suite, every one of whose fifteen
tests spells its data that way. A verbatim port hands a symbolic
`RandomReal[{}, {1050, 1050}]` to `Dot`, `Fourier`, `Eigenvalues`, `Transpose`,
`SingularValueDecomposition` and `LinearSolve` — all of which return in ~0.1 ms
having computed nothing. Eight of the fifteen tests would have reported a
fictitious 10–100× win over Mathematica.

Nothing about the timings looks wrong. The only thing that catches it is the
check value, which is why `CHECK-FAIL` outranks every timing in this job.

(The WolframMark experiment itself was removed — see *Why there is no WolframMark
group* below. The bug it exposed is real and outlives it.)

## The 30 experiments

**Group A — symbolic, baseline `sympy`.** The half of the system that had no
external baseline before this directory existed.

| # | experiment | subsystem |
|---|---|---|
| 01 | `integrate-rational` | `src/calculus/intrat.c`, `src/parfrac.c` |
| 02 | `integrate-transcendental` | `src/calculus/risch_*.c` |
| 03 | `simplify` | `src/simp/` |
| 04 | `solve` | `src/solve.c`, `src/poly/solvepoly.c` |
| 05 | `factor-univariate` | `src/poly/facpoly.c`, `zupoly.c` |
| 06 | `factor-multivariate` | `src/poly/mvfactor.c`, `mpoly.c` |
| 07 | `groebner` | `src/poly/groebner*.c` |
| 08 | `series-limit` | `src/calculus/series.c`, `gruntz.c` |
| 09 | `sum-product` | `src/sum/`, `src/product/` |
| 10 | `polynomial-algebra` | `src/poly/subresultants.c`, `src/rat.c` |

**Group B — numerical libraries, baseline `scipy`.**

| # | experiment | subsystem |
|---|---|---|
| 11 | `special-functions` | `src/special_functions/`, `sf_machine.c` |
| 12 | `quadrature` | `src/numerical_calculus/gkadapt.c`, `dequad.c` |
| 13 | `ode-integration` | `src/numerical_calculus/ndsolve*.c` |
| 14 | `root-finding` | `src/numerical_roots/` |
| 15 | `optimization` | `src/findmin.c` |
| 16 | `interpolation` | `src/interp.c` |
| 17 | `fft-and-signal` | `src/fourier.c`, `src/convolutions.c` |
| 18 | `dense-linalg-drivers` | `src/linalg/lapack_bridge.c` |
| 19 | `statistics` | `src/stats/` |
| 20 | `curve-fitting` | `src/fit.c`, `src/linalg/lstsq.c` |

**Group C — the array substrate.** Exactly the eight open roadmap items ranked in
`docs/experiments/README.md`, so the job tracks them continuously instead of once.
These are not algorithms: a gap here is packing, dispatch, or per-call overhead.

| # | experiment | roadmap item |
|---|---|---|
| 21 | `strided-views` | #2 free `Transpose`, sliding window |
| 22 | `small-matrix-lapack` | #3 a 6×6 `Inverse` at 737 µs → ~2 µs |
| 23 | `per-operation-constant` | #7 the ~8 µs floor per array op |
| 24 | `heterogeneous-return` | #1 mixed-dtype return costing 96×/120× |
| 25 | `vector-transcendentals` | #8 Accelerate `vv*` |
| 26 | `structural-ops` | #9 threading `Reverse`/`RotateLeft`/`Part` |
| 27 | `elementwise-binary` | #10 `MapThread[Min]` 7.7 vs 48 GB/s |
| 28 | `span-selector` | #6 `start/step/n` without a position array |

**Group D — uncovered subsystems.** Neither roadmap items nor numeric, which is
why they are not folded into C.

| # | experiment | why | baseline |
|---|---|---|---|
| 29 | `graph-ops` | `src/graph/` — 19 files, never timed by anything | networkx (pure Python) |
| 30 | `string-ops` | `src/strings/` + `regex/` — 28 files, never timed | `re` (compiled C) |

Read D's two median columns apart, not together — they disagree violently and the
disagreement is the finding. Ahead of Mathematica, hundreds of times behind
Python, because the graph accessors are O(n²) and rescan the edge list on every
call while the string operations are fine. A single median over the group would
hide both facts.

**Group E — advanced numerical analysis, baseline `scipy`/`numpy` (one `mpmath`).**
Where group B measures the *basic* numeric surface (1-D smooth quadrature, scalar
root-finding, dense solve/inv/eig, 1-D FFT), group E goes one level deeper into the
same subsystems: matrix decompositions and eigenproblems beyond `Eigenvalues`,
multidimensional and oscillatory quadrature, stiff ODEs and PDEs, nonlinear
systems, high-degree polynomial roots, the discrete cosine/sine transforms, and
regularized least squares — plus arbitrary-precision numerics against `mpmath`,
the one non-`scipy` baseline in the group.

| # | experiment | subsystem | baseline |
|---|---|---|---|
| 53 | `matrix-decompositions` | `src/linalg/{ludecomp,qrdecomp,svdecomp,inv,matrank,nullspace}.c` | `scipy.linalg`/`numpy.linalg` |
| 54 | `eigenproblems` | `src/linalg/eigen_*.c` | `scipy.linalg`/`scipy.sparse.linalg` |
| 55 | `vectorized-special-functions` | `src/ndkernels.c`, `src/special_functions/` | `scipy.special` |
| 56 | `multidim-quadrature` | `src/numerical_calculus/{cubature,oscint,dequad,denint,levincoll}.c` | `scipy.integrate` |
| 57 | `stiff-ode-pde` | `src/numerical_calculus/ndsolve_{implicit,adams,mol,stencil}.c` | `scipy.integrate` |
| 58 | `nonlinear-systems` | `src/numerical_roots/{findroot,nsolve_system}.c` | `scipy.optimize` |
| 59 | `polynomial-roots` | `src/numerical_roots/{nroots,nroots_aberth,nroots_jt}.c` | `numpy` |
| 60 | `dct-dst-transforms` | `src/fourier.c` | `scipy.fft` |
| 61 | `regularized-least-squares` | `src/{fit,linalg/lstsq}.c` | `scipy`/`numpy` |
| 62 | `arbitrary-precision` | `src/precision.c`, `src/*/*_mpfr.c` | `mpmath` |

Group E is authored and verified against **Python 3.11** (numpy/scipy/mpmath, see
`requirements.txt`). Pin the interpreter at run time — the runner does not hardcode
a Python and the group needs no code change to select 3.11:

```bash
/usr/local/bin/python3.11 -m pip install -r benchmarks/requirements.txt
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 53,54,55,56,57,58,59,60,61,62 --check-labels
```

### Why there is no WolframMark group

There was one, and it was removed. WolframMark is a **hardware** benchmark: it
holds Mathematica constant and varies the machine, scoring against a reference
system. It answers *"how fast is this laptop at running Mathematica"*, not *"how
good is this CAS"* — so it is not a coverage measure, its sizes encode
Mathematica's performance profile (the 1.2M-point × 11 DFT is a ~1 s test there
and would have taken roughly forty hours against Mathilda's pre-FFTW O(n²)
fallback), and its aggregate score is meaningless across implementations.

It did earn its keep once: it surfaced that `RandomReal[{}, dims]` returns
unevaluated in Mathilda where Mathematica reads it as the default `{0,1}` range —
which would have produced fictitious ~0.1 ms "wins" on eight of its fifteen tests
while computing nothing. That finding is recorded in the changelog; a guard for it
belongs in `tests/`, not in a timing job.

## Reading a group A case correctly

Symbolic timings are frequently dominated by **algorithm choice, not execution
speed**. A group A case reading 40× behind sympy usually means *sympy has a method
we do not*, which is different work from the numeric findings — nearly all of
which were "the fast path existed and was not taken".

That is what the `SLOWER` / `ABSENT` / `INCOMPLETE` split is for. Always read
group A cases together with `ABSENT.md`.

Note also that sympy and networkx are **pure Python**, so a Mathilda win there is
a weaker claim than a win against scipy, numpy or `re` — which are compiled. Each
`.py` docstring says which kind of baseline it is.

---

## Backlog

Cut only to hold a balanced 10/10/10 across groups A–C. The runner discovers
folders, so adding one is a folder and nothing else.

- `definite-integrate` — definite integrals, contour and residue cases
- `algebraic-numbers` — `RootReduce`, `MinimalPolynomial`, radical denesting
- `pattern-rewriting` — `ReplaceAll` / `ReplaceRepeated` dispatch throughput
- `arbitrary-precision` — `N[…, 1000]` against `mpmath`
