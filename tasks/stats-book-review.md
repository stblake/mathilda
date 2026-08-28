# Statistics builtins — review for book coverage (Phase 0)

Date: 2026-08-28. Purpose: survey the statistical builtins to decide which to
feature in the new book sections (Ch 3 tour + Ch 4 deep dive). **Document-only** —
no code changes. Verified against the built binary via `./Mathilda -file`
(scratchpad `statreview.m`).

## Result: all candidate builtins work; nothing broken. All 37 are `\B{}`-linkable.

## Core `src/stats/` (15) — verified behaviour & output form

| Builtin | Verified example | Output | Notes / featured in |
|---|---|---|---|
| `Mean` | `Mean[michelson20]` | `909` (exact) | exact on int; columnwise on matrices. §4.8.1, §3 tour |
| `Median` | `Median[{26,28,29,30,31,32}]` | `59/2` | exact; robust to outliers. §4.8.1 |
| `Commonest` | `Commonest[michelson20]` | `{980, 1000}` | ties → list. §4.8.1 (mode), §4.8.7 |
| `Variance` | `Variance[michelson20]` | `209180/19` (exact Rational, n−1) | exact-arith showcase. §4.8.2 |
| `StandardDeviation` | `StandardDeviation[michelson20]` | `2 Sqrt[52295/19]`, `N→104.926` | exact surd. §4.8.2 |
| `RootMeanSquare` | `RootMeanSquare[michelson20]` | `2 Sqrt[209185]` | §4.8.2 |
| `Min`/`Max` | `Min/Max[michelson20]` | `650` / `1070` | list-module heads. §4.8.3 |
| `Quartiles` | `Quartiles[michelson20]` | `{850, 940, 980}` | five-number summary. §4.8.3, §3 tour |
| `Moment` | `Moment[sk,2]` | `493/10` | raw moment. §4.8.4 |
| `CentralMoment` | `CentralMoment[sk,3]` | `43452/125` | about the mean, /n. §4.8.4 |
| `Skewness` | `Skewness[sk]` | `38624/907 Sqrt[3/907] ~ 2.449` | §4.8.4 |
| `Kurtosis` | `Kurtosis[sk]` | `6124153/822649 ~ 7.444` | Pearson (not excess). §4.8.4 |
| `Covariance` | `Covariance[{1,3/2},{2,11}]` | `9/4`; matrix form OK | conj on 2nd arg. §4.8.5 |
| `Correlation` | `Correlation[{5,3/4,1},{2,1/2,1}]` | `2 Sqrt[3/13]`; matrix unit diag | §4.8.5 |
| `MovingAverage` | `MovingAverage[ser,3]`, weighted `{1,2,1}` | exact rationals | §4.8.6, §3 tour |
| `MovingMedian` | `MovingMedian[ser,3]` | `{5,5,6,3,6,4}` | §4.8.6 |
| `ExponentialMovingAverage` | `EMA[ser,1/3]` | exact rationals; `0.2`→machine | §4.8.6 |

## Adjacent surface (verified)

- `Tally[cat] → {{a,5},{b,3},{c,2}}`, `Counts[cat] → <|a->5,b->3,c->2|>`,
  `Commonest[cat] → {a}`. §4.8.7.
- Columnwise: `Mean[{{1,10},{2,20},{3,30}}] → {2,20}`, `Variance → {1,100}`. §4.8.
- `Histogram` exists (graphics) — figure only, subject to headless constraint.

## Distribution bridge (ml/dist.c) — verified

- `SeedRandom[42]; s = RandomVariate[NormalDistribution[100,15],1000];`
  `Mean[s] → 99.85`, `StandardDeviation[s] → 14.95` (LLN; reproducible via seed).
- `PDF[NormalDistribution[],{0.,1.,2.}] → {0.398942,0.241971,0.053991}` (correct).
- `PDF[UniformDistribution[{0,1}],{-0.5,0.5,1.5}] → {0.,1.,0.}`.
- `LearnDistribution[cav]` → head `LearnedDistribution`.
- Supported distributions: Normal, Uniform (PDF); LearnDistribution methods
  Multinormal/GaussianMixture/ContingencyTable. **Absent:** CDF, InverseCDF,
  Quantile, Probability, Expectation, Mode, GeometricMean, HarmonicMean, BinCounts
  — do NOT reference these in the book.

## Classical datasets chosen (all typed as literal data in examples/)

- **Michelson (1879)** speed-of-light, first 20 integer deviations — exact
  arithmetic + five-number summary + bimodal Commonest.
- **Cavendish (1798)** 29 density-of-Earth measurements (real) — machine-number
  location/spread.
- **Median-robustness** small set with a −44 outlier.
- **Anscombe's quartet (1973)** — headline covariance/correlation example.
- A short integer series — moving statistics.
- A small categorical list — frequencies.
- `NormalDistribution[100,15]` sample — the distribution bridge.

## Writing order (guarantees verified prose)

1. Write all `.m` example files (below `book/examples/`).
2. `cd book && make examples`.
3. Read `book/generated/**` transcripts.
4. Write `.tex` prose citing the actual printed outputs (qualitative for the RNG
   sample mean; exact for everything deterministic).
