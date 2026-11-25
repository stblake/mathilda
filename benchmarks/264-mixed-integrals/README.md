# 264: a radical on top of a Log/Exp tower

500 more integrals of class M, outside the easy classes. These are the integrals
FriCAS closed none of in benchmarks 262 (216) and 263 (400). As before, the benchmark
tests Mathilda's algorithmic integrators against FriCAS, with Axiom, Maxima and sympy
as comparison columns.

## The corpus

| class | cases | grammar | methods run | integrands |
|---|---|---|---|---|
| M | 500 | M4 | Automatic, PMT, Bronstein | a radical whose radicand holds an exponential of x |

Every integrand is built on a unit of the tower. With `g = P^2 - c E^h`,

    (D[P] - D[h] P/2) / Sqrt[g]   integrates to   Log[P + Sqrt[g]] - h/2,

and no kernel divides it. M4 draws that shape and perturbations of it. The answer is
required to be larger than the integrand in compact form (ratio 1, as for M in 262 and
263). All 500 integrands are distinct, none is in an earlier corpus or generator pool
(256, 257, 259-263), and all are class `none` in `easy.tsv`.

## How it was built

`gen_compact.py --plan 264`:
- **Attempts.** 50 M4 seeds (401-450), 15 cases each, made 4268 attempts.
- **Rejected as easy.** 568 candidates were rejected as degenerate radicands.
- **Kept.** 750 cases were accepted; after de-duplication and the exclusions 605 were
  available, and the first 500 were kept.
- **Reference answers.** All 500 come from Mathilda. FriCAS, Axiom and Maxima
  answered none of the 45 integrands sent to them.

## Verification

Every answer is judged by `PossibleZeroQ[D[F, x] - f]` in Mathilda (`pzq_verify.py`,
v0.6744): zero is verified, nonzero is wrong (an antiderivative must hold for every
complex x), and anything else is unverified. The systems' own three-point checks are
recorded but decide nothing.

## Results (2026-11-26, Mathilda v0.6744, FriCAS 1.3.6)

| run | closed (500) | wrong | how the rest end |
|---|---|---|---|
| Mathilda Automatic | **500** | 0 | |
| Bronstein | **500** | 0 | |
| ParallelMixedTower | 0 | 0 | 493 unevaluated, 7 timeouts |
| FriCAS | 0 | 0 | 500 errors: "implementation incomplete (has polynomial part)" 270, "(constant residues)" 230 |
| Axiom | 0 | 0 | 500 errors: "implementation incomplete (constant residues)" |
| Maxima | 0 | 0 | 500 unevaluated |
| sympy | 0 | 0 | 231 unevaluated, 269 timeouts |

- **Mathilda's Bronstein method closes every case.** The median time is 0.13 s, and
  Automatic's is 0.16 s. Every answer is verified by PossibleZeroQ.
- **No other system closes any.** FriCAS stops in its own integrator with
  "implementation incomplete" on every one.

## Reproducing

```bash
python3 tools/integration_shootout/gen_compact.py --plan 264 --jobs 6
python3 tools/integration_shootout/easy_classify.py --bench benchmarks/264-mixed-integrals
SHOOTOUT_TRUNC=1000000 python3 tools/integration_shootout/shootout.py \
    --bench benchmarks/264-mixed-integrals --systems mathilda,fricas,axiom,maxima,sympy
python3 tools/integration_shootout/pzq_verify.py --bench benchmarks/264-mixed-integrals
python3 tools/integration_shootout/compare_fricas.py
python3 tools/integration_shootout/latex_report.py --bench benchmarks/264-mixed-integrals
```

`report/report.pdf` opens with a page of every integral. It then shows, for each case,
the integral and every system's answer, verbatim and typeset, with its PossibleZeroQ
verdict and, for an error, the system's own message.
