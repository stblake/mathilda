# Benchmark 88 — Thue-equation stress tests

An adversarial correctness-and-performance harness for the Thue-equation solver
(`Solve[F(x,y) == m && Element[{x,y}, Integers], {x,y}, Integers]`,
`src/solvethue.c`).  The goal is to **break it** — find wrong answers, crashes,
timeouts, and bottlenecks — not to look good.

## The oracle: PARI/GP `thue()`

Every case is cross-checked against **PARI/GP's `thue()`**, the gold-standard
reference: it solves *arbitrary* Thue equations (any degree, any `m`,
non-monogenic fields), unconditionally.  A case where Mathilda returns a finite
set that **differs** from PARI's is a completeness bug (`WRONG`).  This is far
stronger than the three worked examples the solver was written against — it
independently checks the Baker + LLL bound on ~100 equations.

```bash
python3 run.py                       # all cases (needs `gp` on PATH)
python3 run.py --only binom-quartic  # label-substring subset
MATH_TIMEOUT=30 PARI_TIMEOUT=20 python3 run.py
```

Writes `REPORT.md` (ranked, by verdict/family) and `results.json`.  Exits
nonzero on any `WRONG` or `CRASH`.

## The randomized grid (`grid.py`)

`cases.py` is *curated* (hand-picked forms). `grid.py` is its *randomized*
counterpart: a deterministic-seeded generator of random binary forms (degree
3–6, mixed `m` including `|m| != 1`), each cross-checked against PARI `thue()`
through the same runners — hundreds of forms nobody chose. It earns its keep on
the forms Mathilda *solves* (most hard random forms `DECLINE`, which is safe),
where a set differing from PARI is a completeness bug; generation is weighted
toward the solve paths so `CORRECT` is genuinely exercised. Reproducible: a fixed
seed regenerates the identical corpus, so any `WRONG` is a stable, re-runnable
failure. Writes `GRID_REPORT.md` + `grid_results.json`; exits nonzero on
`WRONG`/`CRASH` only.

**Oracle adjudication.** PARI `thue()` is the reference but not infallible — on
some totally-imaginary fields it silently returns an *incomplete* set. So on any
disagreement the grid does not blindly trust PARI: it checks each disputed point
against `F(x,y)==m` directly (soundness needs no box) and classifies
`MATHILDA_WRONG` (a real bug — a spurious point, or a genuine PARI solution
Mathilda lacks) vs `PARI_WRONG` (Mathilda is sound and holds a genuine solution
PARI missed). Only `MATHILDA_WRONG`/`CRASH` fail the run.

```bash
python3 grid.py                  # 300 cases, seed 20260820
python3 grid.py --n 400 --seed 20260820
```

Last run (seed `20260820`, 400 cases): **278 CORRECT / 119 DECLINE / 1 PARI_WRONG
/ 0 WRONG / 0 CRASH**. Every form Mathilda solved matched PARI except one, where
PARI `thue()` was the incomplete side (`x^4-2x^3y+4x^2y^2-3xy^3+y^4 == 5` over
`Q(zeta_5)`: PARI `[]`, true set `{(1,2),(-1,-2)}`, brute-verified).

## The corpus (`cases.py`)

~100 equations chosen to stress every code path and every known weakness:

| family | what it probes |
|---|---|
| binomial cubic `x³−dy³=m` | Delone–Nagell; **fundamental-unit size grows with `d`** |
| simplest / cyclic cubics | monogenic cyclic fields, many solutions (Thomas) |
| binomial quartic `x⁴−dy⁴=m` | the **ℚ-dependent (subfield-unit) degenerate lattice** |
| general quartics | signatures (4,0)/(2,1)/(0,2), biquadratic subfields |
| quintic / sextic / septic | rank ≥ 3–4: unit search + Baker precision |
| non-monogenic | Gate-1 decline (Round-2 not implemented) |
| high `m` (|m|>1) | out of the |m|=1 scope |
| reducible forms | should decline (or reduce) |
| adversarial precision | clustered roots, big coefficients, large regulators |

Each case stores the defining polynomial `P(x)=F(x,1)` as a coeff vector; both
the Mathilda `Solve[...]` and the PARI `thue(thueinit(P,0), m)` are built from
it, so the two halves solve the identical equation.

## Verdicts

`CORRECT` (matches PARI) · `WRONG` (**differs — bug**) · `CRASH` (**bug**) ·
`TIMEOUT` (bottleneck) · `DECLINE` (unevaluated; PARI solved — an honest
coverage gap) · `UNVERIFIED` (PARI could not cross-check in the budget).

## Findings (see `REPORT.md` for the live table)

**Correctness: no bugs.** Of ~100 cases, **48 `CORRECT`, 56 `DECLINE`, 0
`WRONG`, 0 `CRASH`**: every equation Mathilda solves, it solves *completely and
correctly* — its set is identical to PARI's, including the hard ones (the
conductor-7 Thomas cubic's 9 solutions, the biquadratic field's 12, ℚ(2^{1/4})'s
subfield-degenerate lattice, binomial quartics, and even some quintic/sextic/
septic forms). This is strong independent evidence the Baker/de-Weger bound is
right, not merely tuned to the three original targets.

**The dominant gap is the fundamental-unit search (Gate 2).** The unit engine
finds units by an integer-coefficient box search, and a fundamental unit's
coordinates grow with the field regulator:

| field | fundamental unit | regulator | found? |
|---|---|---:|---|
| ℚ(∛2)  | `1 − θ` (coords ≤ 1) | 1.35 | ✓ |
| ℚ(∛6)  | `3θ² − 6θ + 1` (coord 6) | 5.79 | ✓ (adaptive box) |
| ℚ(∛15) | `−12θ² + 30θ − 1` (coord 30) | 9.69 | ✗ decline |
| ℚ(∛41) | 11-digit coordinates | 56.3 | ✗ decline |

A box search is exponential in the coordinate size, so it is the wrong
algorithm for anything but small regulators. The honest behaviour is a **safe
decline** (never a wrong answer). The proper fix is Voronoi's algorithm /
reduction-based unit finding (what PARI's `bnfinit` does); until then a
degree-adaptive box (cubic ≤ 12) catches only the smallest-regulator gaps.

**By-design scope declines (all safe, all confirmed against PARI):**
- **Non-monogenic fields** (Gate 1): many pure quartics `ℚ(d^{1/4})`, some
  cubics (`ℚ(∛17)`, index 3). Needs a Round-2 maximal-order step.
- **|m| ≠ 1**: the unit reduction assumes `β = x−θy` is a unit; general `m`
  needs μ-enumeration over bounded-norm ideal representatives.
- **Reducible forms**: could split into linear/Pell sub-problems; not wired.
- **Degree ≥ 5** and large-regulator quartics: rank/precision out of reach.

**Performance.** Where Mathilda solves, it is fast: cubic solves run in
~3–30 ms and the degenerate quartic in ~7 ms — a full Baker-bound +
LLL-reduction each. The slowest correct solve is the many-solution conductor-7
Thomas cubic (the L-trick enumeration), still well under a second.

## The roadmap this benchmark writes

1. **Voronoi / reduction-based fundamental units** — by far the highest-value
   item: it converts most cubic (and many quartic) `DECLINE`s into solves.
2. **Round-2 maximal order** — unlocks non-monogenic fields.
3. **General `m`** via bounded-norm ideal representatives.
4. **Reducible-form reduction** to the linear/Pell sub-solvers.
