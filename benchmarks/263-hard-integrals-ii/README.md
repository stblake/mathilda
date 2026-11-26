# 263: more integrals outside the easy classes

The second corpus generated outside the easy classes, after
[benchmark 262](../262-hard-integrals/README.md). 3191 of the 4876 cases of 256-261
turned out to be easy (rational, rationalised by a substitution, or closed by
`DerivativeDivides`), and every comparison now leaves them out. 263's 1850 cases
replace most of them, so that the comparison rests on integrals that need a decision
procedure: 4535 cases across 256-263.

As in 261 and 262, it tests Mathilda's algorithmic integrators (Bronstein, Trager,
RischTranscendental and ParallelMixedTower) against FriCAS, with Axiom, Maxima and
sympy as comparison columns.

## The corpus

| class | cases | grammar | methods run | integrands |
|---|---|---|---|---|
| T | 500 | T3 | Automatic, PMT, RischTranscendental | a Log/Exp tower, no radical |
| A | 550 | A3 | Automatic, PMT, Trager | one radical, no Log/Exp |
| M | 400 | M4 | Automatic, PMT, Bronstein | a radical on top of a Log/Exp tower |
| B | 400 | B4 | Automatic, PMT, Bronstein | a radical below a Log/Exp kernel |

Every integrand is distinct, and none is in an earlier corpus or generator pool
(256, 257, 259, 260, 261 and 262). All 1850 were class `none` in `easy.tsv` when it was built; after the constant factors were removed (2026-11-26) the corpus holds 1801, of which 37 are `dd` at v0.6760 and left out.

By family (the system that supplied the reference antiderivative): Mathilda 1809,
FriCAS 30 (`cA-fricas`, `cT-fricas`), Maxima 10, Axiom 1.

## How it was built

`gen_compact.py --plan 263`, with 262's search and these changes:

- **B4, a new grammar.** B3 repeated itself in 262 (110 distinct integrands in 118
  cases). B4 builds the derivative of a product holding a radical of a cubic or
  quartic under `Exp` or `Log`. Its base shapes are `D[R E^y]`, `D[R Log[S + y]]`,
  `D[P y E^(k x)]`, `D[y Log[S + c]]`, `D[R y Log[x]]` and
  `D[x E^y] + D[R Log[x + y]]`. They are perturbed by dropping a numerator term or by
  adding a multiple of `R base/(x + c)`. Since F is then no larger than f, B runs at
  ratio 0 (any elementary F). M runs at ratio 1 as in 262; T and A at ratio 2.
- **Fresh seeds** (301 and up) for every class.
- **More A seeds.** A3 repeats itself too: its first 35 seeds gave 467 new distinct
  integrands from 1050 accepted cases. Twenty more seeds (336-355) were run with
  `--extend` to reach 550.
- **Easy integrands rejected at the source.** 165 generator runs made 13673 attempts;
  471 candidates were rejected as `dd` and 604 as `degenerate`.
- **Replacement.** The check after the build (`easy_classify.py`, v0.6733) found 5 A
  cases that `DerivativeDivides` closes. They had slipped through on generator
  timeouts and were replaced by spares from the pool that passed the same check.
  `gen_stats.json` records the replacement.

## Verification

Every answer of every system is judged by `PossibleZeroQ[D[F, x] - f]` in Mathilda,
the other systems' answers translated first (`pzq_verify.py`):
- zero: verified;
- nonzero: **wrong** (an antiderivative must hold for every complex x);
- anything else: unverified, never wrong.

The own three-point checks are recorded and decide nothing. On 2026-11-26 the saved run
was re-judged with v0.6744's PossibleZeroQ, and every FriCAS and Axiom ERROR record was
re-run with the fixed adapter (`refresh_errors.py`), which records each error's message.

## Current run (2026-11-26)

The integrands were rewritten without their constant factors (`normalize_constants.py`;
a case that became a duplicate was left out) and every system was run
again, Mathilda at v0.6750. The easy labels were re-checked at v0.6760, whose
DerivativeDivides no longer depends on a constant factor: 6 of the run cases are
easy and left out, so 1764 are compared. Closed means verified by PossibleZeroQ;
the full tables are in [`REPORT.md`](REPORT.md) and `report/report.pdf`.

| run | closed (% of the cases in the classes the run covers) | wrong |
|---|---|---|
| Mathilda Automatic | 100% (1761) | 0 |
| Mathilda ParallelMixedTower | 73% (1285) | 0 |
| Mathilda RischTranscendental | 100% (494) | 0 |
| Mathilda Trager | 94% (446) | 0 |
| Mathilda Bronstein | 100% (795) | 0 |
| FriCAS | 71% (1255) | 0 |
| Axiom | 66% (1163) | 4 |
| Maxima | 60% (1050) | 19 |
| sympy | 50% (883) | 28 |

## Earlier results (run 2026-11-25 at Mathilda v0.6733; verdicts by PossibleZeroQ, v0.6744)

Closed with a verified answer, out of the cases each run applies to:

| run | T (500) | A (550) | M (400) | B (400) | all | wrong |
|---|---|---|---|---|---|---|
| Mathilda Automatic | 500 | 526 | 400 | 400 | **1826 / 1850** | 0 |
| RischTranscendental | 500 | · | · | · | 500 / 500 | 0 |
| Trager | · | 517 | · | · | 517 / 550 | 0 |
| Bronstein | · | · | 400 | 400 | **800 / 800** | 0 |
| ParallelMixedTower | 497 | 454 | 0 | 388 | 1339 / 1850 | 0 |
| FriCAS | 500 | 427 | 0 | 389 | 1316 / 1850 | 0 |
| Axiom | 500 | 333 | 0 | 370 | 1203 / 1850 | 4 |
| Maxima | 486 | 403 | 0 | 238 | 1127 / 1850 | 19 |
| sympy | 447 | 385 | 0 | 116 | 948 / 1850 (391 special) | 28 |

- **Mathilda against FriCAS.** Mathilda (any method) closes 534 cases FriCAS does
  not, and FriCAS closes 1 Mathilda does not.
  - **M:** FriCAS closes none of the 400; Bronstein closes all.
  - **A:** FriCAS closes 427 of 550; Trager closes 517 and Automatic 526.
  - **FriCAS's misses:** 420 errors, 103 timeouts and 11 unevaluated. The errors are
    FriCAS's own "implementation incomplete": has polynomial part 224, constant
    residues 176, trace 0 18, and one other.
- **FriCAS: 11 false "no".** On 11 B4 integrands that are elementary by construction,
  FriCAS returns `integral(...)`, its answer for "no elementary antiderivative". Each
  is the derivative of `E^y R`, with y the square or cube root of a cubic or quartic
  (`cB-301-25`: `f = D[E^Sqrt[2 + x^2 + x^3]/x^2, x]`). A minimal case in FriCAS
  1.3.6 is

  ```
  integrate(D(exp(sqrt(x^3+1))/(x+1), x), x)
  ```

  which should give `exp(sqrt(x^3+1))/(x+1)`. With `1/x` in place of `1/(x+1)`,
  FriCAS integrates it. Mathilda closes all 11.
- **Soundness.** No Mathilda method gives a wrong answer or a false "no", and no FriCAS
  answer is wrong. Maxima's 19 are antiderivatives valid only on part of the plane.
- **Automatic against Trager.**
  - Trager closes 22 A cases in 11-15 s where Automatic gave up at 20-27 s, because the
    stages before Trager take the time. Fixed in v0.6744: an expired cascade budget now
    resumes Trager, and these close. The table shows the run at v0.6733.
  - Automatic closes 31 A cases where Trager times out (27) or declines (4).
  - 262 showed the same gap for Bronstein.
- **Speed** (1259 cases both close, each system's own time): median 222 ms against
  FriCAS's 170 ms, and Mathilda is faster on 659 of 1259. Mathilda's 90th percentile
  is 3.6 s against FriCAS's 250 ms, and 381 Mathilda cases take over 1 s. The slow
  cases are A, where Trager's median is 2.0 s.

### How the run went

The first pass shared the machine with other heavy jobs and had 771 TIMEOUTs. All of
them were dropped from the cache and rerun at a load average under 10. After the
rerun:
- Mathilda's Automatic and RischTranscendental have no timeouts. Trager has 29 and
  PMT 105.
- FriCAS's 102 timeouts are unchanged.

## Reproducing

```bash
python3 tools/integration_shootout/gen_compact.py --plan 263 --jobs 8   # then --extend for A 336-355
python3 tools/integration_shootout/easy_classify.py --bench benchmarks/263-hard-integrals-ii
SHOOTOUT_TRUNC=1000000 python3 tools/integration_shootout/shootout.py \
    --bench benchmarks/263-hard-integrals-ii --systems mathilda,fricas,axiom,maxima,sympy
python3 tools/integration_shootout/compare_fricas.py
python3 tools/integration_shootout/latex_report.py --bench benchmarks/263-hard-integrals-ii
```

`report/report.pdf` shows every answer from every Mathilda method and every other
system, how each was verified, and, in its appendix, every integral of the corpus.
