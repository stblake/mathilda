# 262: integrals outside the easy classes

A follow-on to [benchmark 261](../261-compact-integrals-x4/README.md). Most of 257's
and 261's integrands turned out to be **easy**: a rational function, one a
substitution makes rational, or one `DerivativeDivides` closes. Such an integral needs
no decision procedure. `easy_classify.py` (easy.m) classes 3191 of the 4876 cases of
256-261 as easy, and every comparison now leaves them out. 262 is 1000 integrals
generated with the easy classes excluded at the source.

As in 261, it tests the completeness of Mathilda's algorithmic integrators (Bronstein,
Trager, RischTranscendental and ParallelMixedTower) against FriCAS, the other
algorithmic integrator. Axiom (Tim Daly's distribution, the system FriCAS forked),
Maxima and sympy are comparison columns.

## The corpus

| class | cases | methods run | integrands |
|---|---|---|---|
| T | 330 | Automatic, PMT, RischTranscendental | a Log/Exp tower, no radical |
| A | 336 | Automatic, PMT, Trager | one radical, no Log/Exp |
| M | 216 | Automatic, PMT, Bronstein | a radical on top of a Log/Exp tower |
| B | 118 | Automatic, PMT, Bronstein | a radical below a Log/Exp kernel |

There is no R class: every rational integrand is easy.

By family (the system that supplied the reference antiderivative): Mathilda 974
(`c<class>-mathilda`), Axiom 15 (`cA-axiom`), Maxima 11 (`cA-maxima`, `cB-maxima`).

## How it was built

`gen_compact.py --plan 262` uses 261's search, with these changes:

- **Easy integrands are rejected at the source.** Every candidate is classed by
  easy.m before it is accepted. Over 14312 attempts, 1496 candidates were rejected as
  `dd` and 336 as `degenerate`. No rational, linear-radical, Euler or nested-radical
  candidate came up.
- **Ratio 2** (F more than twice f in compact form), as in 257; M runs at ratio 1.
- **The M4 grammar** for most M seeds (M3 for the first four), and fresh seeds (201 and
  up). The integrands of 256, 257, 259, 260 and 261 are excluded.
- **Reference systems.** An integrand Mathilda cannot close goes to FriCAS, Axiom
  and Maxima. An answer from one of them that verifies puts the case in the corpus,
  in that system's family.

B3 repeated itself: the 118 B cases hold 110 distinct integrands. Benchmark 263 adds
the B4 grammar for this reason.

Wolfram is never consulted. Cases that no system can close are not in the corpus,
because no reference antiderivative exists for them.

## Verification and tallies

Every answer of every system is judged by one test, `PossibleZeroQ[D[F, x] - f]` in
Mathilda (`pzq_verify.py`). The other systems' answers are translated into Mathilda
first (`crosscheck.to_wl`).
- zero: verified;
- nonzero: **wrong** (the standard is an antiderivative valid for every complex x);
- anything else (undecided, untranslatable, or an answer text cut when the run was
  recorded): unverified, never wrong.

Each system's own three-point check is still recorded, but decides nothing.

The run itself is unchanged since 2026-11-25. Two corrections were applied to the
saved results afterwards:
- **2026-11-26: re-judged.** Every answer was re-judged with v0.6744's PossibleZeroQ.
- **2026-11-26: errors re-run.** Every FriCAS and Axiom ERROR record was re-run with
  the fixed adapter, which records each error's message (`refresh_errors.py`). Among
  FriCAS's 230 errors were runs that had answered but lost the answer to a fault in the
  adapter's output file; those now count.

## Results (run 2026-11-25 at Mathilda v0.6732; verdicts by PossibleZeroQ, v0.6744)

Closed with a verified answer, out of the cases each run applies to:

| run | T (330) | A (336) | M (216) | B (118) | all | wrong |
|---|---|---|---|---|---|---|
| Mathilda Automatic | 330 | 317 | 216 | 107 | **970 / 1000** | 0 |
| RischTranscendental | 330 | · | · | · | 330 / 330 | 0 |
| Trager | · | 312 | · | · | 312 / 336 | 0 |
| Bronstein | · | · | 216 | 118 | **334 / 334** | 0 |
| ParallelMixedTower | 328 | 271 | 1 | 60 | 660 / 1000 | 0 |
| FriCAS | 330 | 268 | 0 | 48 | 646 / 1000 | 0 |
| Axiom | 330 | 241 | 0 | 0 | 571 / 1000 | 2 |
| Maxima | 327 | 238 | 0 | 73 | 638 / 1000 | 3 |
| sympy | 293 | 228 | 0 | 0 | 521 / 1000 (228 special) | 15 |

- **Mathilda against FriCAS.** Mathilda (any method) closes 351 cases FriCAS does
  not; FriCAS closes 13 Mathilda does not, all in class A (`cA-axiom`: Trager times
  out on them).
  - **M:** FriCAS closes none of the 216; Bronstein closes all.
  - **B:** FriCAS closes 48 of 118; Bronstein closes all.
  - **FriCAS's misses:** 219 errors and 131 timeouts. The errors are all FriCAS's own
    "implementation incomplete": constant residues 109, has polynomial part 107, trace 0
    3.
- **Automatic was below the dedicated method in 14 cases** (11 B, 3 A). Automatic gave
  up at about 20 s, while Bronstein or Trager needs 10-26 s on these sixth-root
  radicals, for example `cB-201-559`, `Log[x]/(x^6 (3 + 2 x^6)^(1/6))`. Fixed in
  v0.6744: an expired cascade budget now resumes the method, and these close. The table
  shows the run at v0.6732. Conversely, Automatic closes 8 A cases where Trager times
  out or declines.
- **Trager** leaves 24 A cases undecided: 21 timeouts and 3 unevaluated.
- **No special-function answers** from Mathilda or FriCAS.
- **Soundness.** No Mathilda method gives a wrong answer. FriCAS has 4 unverified
  answers and none wrong; Axiom 2 wrong, Maxima 3, sympy 15.
- **Speed** (623 cases both close, each system's own time, under heavy machine load):
  - the dedicated methods have a median of 328 ms against FriCAS's 100 ms, and are
    faster on 300 of 623;
  - 154 Mathilda cases take over 1 s. That is the `.m` Bronstein and Trager code on
    the hard cases, and the open performance item this corpus shows.

### How the run went

The first pass shared the machine with other heavy jobs and had 544 TIMEOUTs. Those records were dropped from the cache and rerun; 522 remain,
at a load average of 50-66. Most of them are sympy (224), FriCAS (130) and PMT.
Both systems' times are inflated by the load.

## Reproducing

```bash
python3 tools/integration_shootout/gen_compact.py --plan 262 --jobs 12
python3 tools/integration_shootout/easy_classify.py --bench benchmarks/262-hard-integrals
SHOOTOUT_TRUNC=1000000 python3 tools/integration_shootout/shootout.py \
    --bench benchmarks/262-hard-integrals --systems mathilda,fricas,axiom,maxima,sympy
python3 tools/integration_shootout/compare_fricas.py
python3 tools/integration_shootout/latex_report.py --bench benchmarks/262-hard-integrals
```

`report/report.pdf` shows every answer from every Mathilda method and every other
system, how each was verified, including the cross-check verdict, and, in its
appendix, every integral of the corpus. All 1000 cases are `none` in `easy.tsv`. Long answers are broken over `multline*`
lines, and long integers in the verbatim text wrap.
