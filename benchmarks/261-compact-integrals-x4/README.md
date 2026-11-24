# 261: compact integrands, integrals at least four times larger

A follow-on to [benchmark 257](../257-compact-integrals/README.md). It tests the
completeness of Mathilda's algorithmic integrators, Bronstein, Trager,
RischTranscendental and ParallelMixedTower, against FriCAS, the other algorithmic
integrator, on integrands from mixed towers: transcendental (Log/Exp) extensions
together with one simple radical extension. Axiom (Tim Daly's distribution, the
system FriCAS forked), Maxima and sympy are comparison columns.

Every integrand f is compact, and its elementary antiderivative F is more than
**four** times larger in compact form: `4 LeafCount[f] < LeafCount[compactForm[F]]`.
257 required only twice.

## The corpus

| class | cases | methods run | integrands |
|---|---|---|---|
| R | 50 | Automatic, PMT, RischTranscendental | rational |
| T | 250 | Automatic, PMT, RischTranscendental | a Log/Exp tower, no radical |
| A | 250 | Automatic, PMT, Trager | one radical, no Log/Exp |
| M | 550 | Automatic, PMT, Bronstein | a radical on top of a Log/Exp tower |
| B | 450 | Automatic, PMT, Bronstein | a radical below a Log/Exp kernel |

The corpus is weighted to the mixed classes M and B: 1000 of its 1550 cases.

## How it was built

`gen_compact.py --plan 261` uses 257's search, with these changes:

- **Ratio 4, strict.** F must be more than four times f, measured in compact form,
  so a padded F cannot pass.
- **Fresh seeds** (101 and up). The integrands of benchmarks 256, 257, 259 and 260
  are excluded.
- **Fair mode.** An integrand Mathilda answers wrongly, or only with special
  functions, goes to FriCAS like one Mathilda leaves unevaluated. A Mathilda miss
  therefore stays in the corpus whenever FriCAS can close it (family
  `c<class>-fricas`). In 257 such integrands were silently dropped.
- The FriCAS acceptance applies the same strict ratio.

Wolfram is never consulted. Cases that neither Mathilda nor FriCAS can close are not
in the corpus: no reference antiderivative exists for them.

## Verification and tallies

Every system first verifies its own answer, at three rational points at 40 digits
(see `tools/integration_shootout/adapters.py`). For FriCAS and Axiom, an exact
`normalize(D(r, x) - f) = 0` also counts as verified when the point check fails
on a branch choice. Those cases are listed separately.

**One rule for "wrong".** An answer that fails its own check, or that its own check
cannot evaluate, is translated into Mathilda (`crosscheck.py`). It is then judged by
`PossibleZeroQ[D[F, x] - f]`, the test the generator applies to Mathilda's answers:
- zero: verified;
- nonzero: **wrong**;
- undecided (60 s), or an answer that does not translate: unverified, never wrong.

A True counts only when `D[F, x] - f` evaluates numerically at two test points.
`PossibleZeroQ` answers True for an expression it cannot evaluate (PZQ-root-1 in
`tasks/todo.md`), so FriCAS's nested `rootOf` is not translated. No answer that
needed the cross-check contains one.

Checks of the cross-check:
- **Translator, both ways.** 40 answers that had passed their own checks (10 per
  system) all come out zero. The same answers multiplied by `1 + x/10^6` all come
  out nonzero.
- **Every verdict against the corpus points.** Each of the 122 answers
  `PossibleZeroQ` found zero also passes the 40-digit check at the three corpus
  points in Mathilda. Each of the 167 it found nonzero fails it at all three.

Elementary and special-function answers are counted apart for every system, by
the same function list (`translate.special_q`).

## Results (2026-11-24, Mathilda v0.6714, FriCAS 1.3.6)

Closed with a verified answer, out of the cases each run applies to:

| run | R (50) | T (250) | A (250) | M (550) | B (450) | all | wrong |
|---|---|---|---|---|---|---|---|
| Mathilda Automatic | 50 | 250 | 250 | 549 | 450 | **1549 / 1550** | 0 |
| RischTranscendental | 49 | 250 | · | · | · | 299 / 300 | 0 |
| Trager | · | · | 250 | · | · | 250 / 250 | 0 |
| Bronstein | · | · | · | 548 | 450 | 998 / 1000 | 0 |
| ParallelMixedTower | 0 | 250 | 240 | 497 | 386 | 1373 / 1550 | 0 |
| FriCAS | 50 | 249 | 243 | 481 | 449 | 1472 / 1550 | 0 |
| Axiom | 49 | 239 | 206 | 383 | 442 | 1319 / 1550 | 156 |
| Maxima | 30 | 108 | 124 | 214 | 421 | 897 / 1550 | 0 |
| sympy | 30 | 173 | 26 | 24 | 107 | 360 / 1550 (41 special) | 11 |

- **Special-function answers.** Neither Mathilda nor FriCAS gives any; every count
  above is elementary except sympy's 41.
- **Mathilda against FriCAS.**
  - Best Mathilda method: Mathilda alone closes 78 cases and FriCAS alone closes none.
  - The dedicated method of each class (RischTranscendental, Trager, Bronstein):
    77 against 2.
  - FriCAS misses 78 cases, and 69 of them are in class M, a radical over a Log/Exp
    tower. Over the whole corpus its misses are 32 timeouts, 43 errors and
    3 unverified answers.
- **What the cross-check changed.**
  - FriCAS: 97 of the 100 answers that failed or could not complete their own check
    are zero. So FriCAS's own point check, not its answers, was wrong on those, at a
    branch of `complexNumeric`.
  - Axiom: 17 of 173 are zero, and 156 are wrong by both checks.
  - sympy: 22 `Piecewise` answers do not translate and stay unverified.
- **Speed** (cases both close, each system's own time): Mathilda's dedicated methods
  have a median of 152 ms against FriCAS's 170 ms, and are faster on 787 of 1470.
  FriCAS has the shorter tail, with a 90th percentile of 220 ms against 538 ms.
- **Mathilda timeouts.** Every one ran the full 30 s, so none came from a poisoned
  session. There are 116 in PMT and 1 each in Automatic (cM-102-291),
  RischTranscendental and Bronstein.

### How the run went

- **Generator.** Building the corpus hit a load-dependent evaluator bug: aborts escape
  their `TimeConstrained` frame (EVAL-abort-2 in `tasks/todo.md`, replay in
  `tests/data/eval_abort2/`). It only made the generator skip attempts, so seeds
  were added with `--extend` until every class reached its quota.
- **First pass.** It ran at a load average of about 90. Every TIMEOUT, ERROR and
  not-fully-verified record (1409) was then dropped from the cache and rerun at low
  load, and the tables are from that second pass.

## Reproducing

```bash
python3 tools/integration_shootout/gen_compact.py --plan 261 --jobs 12
SHOOTOUT_TRUNC=1000000 python3 tools/integration_shootout/shootout.py \
    --bench benchmarks/261-compact-integrals-x4 --systems mathilda,fricas,axiom,maxima,sympy
python3 tools/integration_shootout/compare_fricas.py
python3 tools/integration_shootout/latex_report.py --bench benchmarks/261-compact-integrals-x4
```

`report/report.pdf` shows every answer from every Mathilda method and every other
system, how each was verified, including the cross-check verdict, and, in its
appendix, every integral of the corpus.
