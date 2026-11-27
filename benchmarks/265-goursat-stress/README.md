# 265: the GeneralisedGoursat stress corpora

1183 pseudo-elliptic integrals `F(x)/R(x)^p` (`p` in 1/2, 1/3, 2/3, 1/4, 3/4) from the
GeneralisedGoursat test corpora, `tests/data/generalised_goursat/*.m`. The benchmark
runs three Mathilda Integrate methods on them, against goursat.py (the sympy reference
implementation that `Method -> "GeneralisedGoursat"` ports). Plain sympy `integrate`
is an opt-in column. The package is still being fixed, so its numbers are a baseline.
A rerun is one command.

```bash
python3 tools/goursat/goursat_stress.py                        # THE TIMED RUN: 30 s/case, 6 jobs; quiet machine only
python3 tools/goursat/goursat_stress.py --correctness          # correctness pass: 60 s/case, 4 jobs, any load
python3 tools/goursat/goursat_stress.py --corpora stress_7 --methods GeneralisedGoursat
python3 tools/goursat/goursat_stress.py --families sqrt-klein --fresh
python3 tools/goursat/goursat_stress.py --with-sympy           # add plain sympy integrate (slowest column)
python3 tools/goursat/goursat_stress.py --goursatpy fresh      # re-time goursat.py too (hours)
python3 tools/goursat/goursat_stress.py --retabulate           # REPORT.md from the last saved run
```

**Status (2026-11-27).** Only the correctness pass has been run. A first timed run was
stopped after 25 minutes, because the machine's load average reached 46 on 16 CPUs
with other agents building. Its partial cache is kept in
`results/load-affected-2026-11-27/` with a NOTE, and it is not a baseline. The
current REPORT.md is the correctness pass, labelled as such. Its times are not a
baseline, and its timeouts are provisional. The timed baseline is the first command
above, run on a quiet machine. Before starting it, check `uptime`: the run records
the load it saw, and REPORT.md prints it.

There are two modes:

- **Timed** (the default). `REPORT.md`, `results/<date>.json` and a
  `history.jsonl` line with `"mode": "timed"`. This is the baseline.
- **Correctness** (`--correctness`). The same runs and the same PossibleZeroQ
  verdicts, at 60 s and 4 jobs. Wrong answers and verdicts do not depend on load. It
  writes `REPORT.md` headed "correctness pass; timings pending a quiet machine",
  `results/<date>-correctness.json`, and a `history.jsonl` line with
  `"mode": "correctness"`.

The run is resumable. Every finished run is cached in `results/cache.jsonl`. For
Mathilda the cache key includes the binary's mtime and the newest `src/internal/*.m`
mtime, so editing the package or rebuilding re-runs only Mathilda. A full run writes
`REPORT.md`, `results/<date>.json` and a line of `history.jsonl`. A partial run
(`--corpora`, `--families`, `--limit`, or a subset of `--methods`) writes
`REPORT-partial.md`. `--check` exits 1 on any wrong answer or false certificate.

## The corpus (`corpus.tsv`)

`tools/goursat/gen_stress_corpus.m` generates the corpus. Each case is one `CASE`
row in the integration shootout's format:

    CASE  id  class  family  expect  f  F  points  f(points)  fullform  corpus  label  truth

| field | value |
|---|---|
| class | `G` (a new class, outside the shootout's R/T/A/M/B) |
| family | the corpus's own `Family` (`sqrt-klein`, `cbrt-mixed`, `pseudoelliptic-sqrt`, `chebyshev`, ...) |
| expect | `elem` / `nonelem` / `unknown`, from the corpus's `Truth` (`elementary` / `non-elementary` / `undecided` or none) |
| F | `None`: the corpora carry no reference antiderivative |
| points | three rational points where f is finite (real with positive radicands where possible). They are used only for the systems' own checks |
| corpus, label, truth | where the case comes from. The label is the key into goursat.py's timings |

| corpus | cases | families | truth |
|---|---|---|---|
| unit | 128 | chebyshev (104), unit (23), cubic-field (1) | Chebyshev's theorem; constructed |
| examples | 29 | examples | none |
| literature | 51 | literature | 32 elementary, 8 non-elementary, 11 none |
| welz | 30 | welz | elementary |
| stress_2026, stress_7, stress_11 | 225 each | hermite-exact, cbrt-symmetric, cbrt-mixed, cbrt-nonelementary, cbrt-second-level, sqrt-klein, sqrt-mixed, sqrt-nonelementary, sqrt-order3 (25 each per seed) | constructed |
| pseudoelliptic_sqrt | 152 | pseudoelliptic-sqrt | elementary |
| pseudoelliptic_cbrt | 118 | pseudoelliptic-cbrt | elementary |

The integrand is the corpus's own, with `t` renamed `x`. It is **not** stripped of
constant factors, unlike 264's integrands, because it must stay the integrand that
goursat.py was timed on.

## What runs

| column | call |
|---|---|
| `Mathilda:Automatic` | `Integrate[f, x]` |
| `Mathilda:GoursatAlgebraic` | `Integrate[f, x, Method -> "GoursatAlgebraic"]` (the C method) |
| `Mathilda:GeneralisedGoursat` | `Integrate[f, x, Method -> "GeneralisedGoursat"]` (the `.m` port, `src/internal/goursat/`) |
| `goursat.py` | the sympy reference: `goursat_integrate` / `integrate_algebraic` |
| `sympy` | plain `sympy.integrate` (sympy 1.14, python3.11). Opt-in, with `--with-sympy` |

- **Limits.** Each Mathilda and sympy case runs under a 30 s per-case limit
  (`--timeout`; 60 s with `--correctness`), with `MemoryConstrained` at 4 GB. Cases run in batches of eight per
  process, and each method is warmed up once per process before any case is timed.
- **GeneralisedGoursat's budget.** The stage's own 30 s budget is turned off
  (``Integrate`$GeneralisedGoursatBudget = Infinity``), so the per-case limit is what
  ends a slow case. That way a timeout is recorded as `TIMEOUT` and not as a decline.
- **GeneralisedGoursat's certificate.** The Method option discards the stage's
  non-elementary certificate and declines. So when the Method call declines, the
  stage symbol ``Integrate`GeneralisedGoursat[f, x]`` is asked again, untimed, only to
  read the certificate. Its time is recorded separately as `cert_time`.
- **goursat.py.** By default its verdicts and times are read from
  `tests/data/generalised_goursat/sympy_<corpus>.tsv` (`tools/goursat/bench_sympy.py`).
  Those runs time the reduction only (no verification), used a 120 s limit, and
  check the result numerically. A TSV has no answer text, so these answers are not
  judged by PossibleZeroQ.
- **Re-timing goursat.py.** `--goursatpy fresh` runs goursat.py again through
  `adapters.run_goursatpy` (in `tools/integration_shootout/adapters.py`, next to the
  sympy adapter). It calls `integrate_algebraic(f, x)` on every case, times only that
  call, and sends its answers through the PossibleZeroQ cross-check. It is slow: the
  recorded TSVs sum to 4.2 CPU-hours, 3.9 of them on the three stress corpora.

## Verification

Every answer from Mathilda and sympy is judged by `PossibleZeroQ[D[F, x] - f]` in
Mathilda, using `crosscheck.run_batch` with a 60 s limit:

- zero is **verified** (`ELEM`, or `SPECIAL` if the answer contains a special function);
- nonzero is **wrong**;
- anything else is **unverified**, never verified. That covers undecided, a
  PossibleZeroQ timeout, and a `D[F, x] - f` that does not evaluate numerically.

There is no numeric fallback for the verdict. Each system's own three-point check is
recorded in the raw record (`nums`), but it decides nothing.

An answer that carries `$Failed` or a package-private symbol is wrong outright. A
certificate (`CERT`) on an `expect elem` case is a **false certificate**. Wrong
answers and false certificates are the hard failures, and REPORT.md lists each one.

## Report

REPORT.md gives the following, per method, overall and per family:

- counts of verified, wrong, unverified, certified, uneval, timeout and error;
- the median time of the decided runs and the maximum time;
- the median LeafCount of verified answers;
- the per-case time ratio against goursat.py;
- a per-corpus table comparing GeneralisedGoursat with goursat.py;
- verdicts that disagree with the corpus labels: a certificate on an `elem` case, or
  a verified elementary answer on a `nonelem` case. The report also tabulates how
  each run ends on the `nonelem` cases, and lists the cases where GeneralisedGoursat
  and goursat.py reach opposite verdicts;
- hand-written findings for a run, kept as a `"notes"` string in its
  `results/<run>.json`, which `--retabulate` reprints;
- the worst offenders: the slowest decided cases, the cases goursat.py decides and a
  method does not, and the cases where a method is slowest relative to goursat.py.

## Registration

The 26x integration benchmarks are registered as `make` targets of their own
(`bench-integrate`, `bench-integrate-compact`, ...), not under `make bench-gap`.
`bench-gap` (`benchmarks/run_all.py`) runs only the folders that hold `.m`/`.py`
pairs, and this folder holds none. The matching target, to add to the makefile next
to `bench-integrate-bottom`:

```make
# `make bench-integrate-goursat` -- benchmark 265, the GeneralisedGoursat stress
# corpora: Automatic, GoursatAlgebraic and GeneralisedGoursat against goursat.py
# and sympy (tools/goursat/goursat_stress.py).  Resumable (cache).
bench-integrate-goursat:
	python3 tools/goursat/goursat_stress.py
```
