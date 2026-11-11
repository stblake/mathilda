# DSolve PDE corpus

The test bed for symbolic PDE solving in `DSolve`, the PDE counterpart of the ODE
corpus one directory up. The roadmap it measures is `DSOLVE_PDE_PLAN.md` at the repo
root, and the current numbers are in `STATUS.md`.

## Source

N. M. Abbasi, *Partial Differential Equations, Maple 2026.1 and Mathematica 15*
(https://www.12000.org/my_notes/pde_in_CAS/pde_in_cas.htm). It has 2012 problems,
collected from:

- Haberman, *Applied PDEs* (5th ed.);
- Logan, *Applied PDEs*;
- Articolo, *PDEs and BVPs with Maple*;
- Olver, *Introduction to PDEs* (2014);
- Debnath, *Nonlinear PDEs* (3rd ed.);
- Polyanin, Zaitsev & Moussiaux, *Handbook of First Order PDEs* (2002);
- Polyanin & Zaitsev, *Handbook of Nonlinear PDEs* (2004);
- help pages and forums.

| file | fetched | sha256 |
|---|---|---|
| `raw/Mathematica_PDE_IN_CAS_problems.txt` (source updated 2026-07-19) | 2026-11-11 | `374ad919f4ed89a31500545e37f1b148d0c62f172e9f2b2eb423d494d51f1d15` |
| `raw/lookup.tsv`, scraped from `maple_2026_1_and_mma_15/Ch1.S4.htm` | 2026-11-11 | (derived) |

`lookup.tsv` gives each problem's section, description, and Abbasi's Mathematica
15 and Maple 2026.1 results (✓ → `Y`). Abbasi's verdicts are **not verified**; his
introduction says a returned solution is assumed correct. Ours are verified.

## Layout

| file | what |
|---|---|
| `DE_PDE_<chap>.m` | the corpus. One file per chapter, and chapter 6 (1585 problems) is split by handbook chapter `6_k`. Generated, so do not edit by hand. |
| `skipped.tsv` | problems the importer refused, with the reason. There is one: `pde-392` uses a `bc` the block never defines, because the source ran as a single session and silently picked up the previous problem's `bc`. |
| `dsolve_pde_corpus_prelude.m` | the verifier, `pdeCheckCode`. |
| `prelude_selftest.m` | the verifier tried on hand-made answers with known verdicts: `./Mathilda -file DSolve_test_status/pde/prelude_selftest.m`, which must print `selftest: N/N`. |
| `test_dsolve_pde_corpus.c` | the runner, built as `dsolve_pde_corpus_tests`. It forks per record using `../corpus_fork.h`. |
| `reports/<chap>.tsv`, `reports/<chap>.md` | runner output, and the ranked gap list per chapter. |
| `STATUS.md` | the dashboard, written by `tools/dsolve_pde_corpus_report.py`. |

Each record has the form `{label, section, eqns, fn, vars, opts, mma, maple}`. For
example:

```
{"pde-16", "2.1.16", {(D[w[x, t], t] + c*D[w[x, t], x] == 0), (w[x, 0] == f[x]), (w[0, t] == h[t])},
  w[x, t], {x, t}, {Assumptions -> c > 0 && x > 0 && t > 0}, "N", "Y"}
```

## Verdicts

| tag | verdict | meaning |
|---|---|---|
| P | PASS | Explicit answer. The PDE residual and every condition residual are proved zero by `PossibleZeroQ`. The arbitrary functions `C[k][…]` are left symbolic, or, if symbolic `PossibleZeroQ` cannot decide, the residual must vanish under both of two generic pinnings. |
| I | PASS-IMPLICIT | `φ == C[1][ψ]`, first order. Proved by the implicit-function rule, identically in (vars, u). |
| S | PASS-SERIES | An `Inactive[Sum]` eigenfunction series. The generic `K[1]` mode solves the homogeneous PDE and boundary conditions, the modes are orthogonal on the interval the conditions bound, and the coefficient is the projection of the data. |
| K | UNK | The answer is neither proved nor refuted. Examples: the heat kernel at `t = 0`, which is a distributional limit; a higher-order implicit relation; a Bessel series. |
| F | FAIL | A residual is proved **nonzero**. This is a wrong answer. |
| N | UNFIT | An IVP or BVP answer that still contains `C[k]`, meaning the conditions were never fitted. |
| U | UNEVAL | Declined, `{}`, timed out after 8 s (or 20 s wall), or crashed. |

There is no numeric acceptance. A `PossibleZeroQ` false negative is a bug in
`PossibleZeroQ`, and is fixed there.

## Running it

```bash
cd tests/build && make dsolve_pde_corpus_tests
ctest -R dsolve_pde_corpus -j8              # the per-chapter gates (label: nightly)
# one chapter by hand; TSV on stdout, summary on stderr:
./dsolve_pde_corpus_tests ../../DSolve_test_status/pde/DE_PDE_3.m \
    ../../DSolve_test_status/pde/dsolve_pde_corpus_prelude.m 99999 99999 \
    > ../../DSolve_test_status/pde/reports/3.tsv
python3 ../../tools/dsolve_pde_corpus_report.py   # STATUS.md + reports/*.md
```

Each chapter test has two gates:

- the **non-PASS** count may not rise above its baseline;
- the **FAIL** count may not rise above its own, separate baseline.

Lower both baselines in `tests/CMakeLists.txt` as method waves land.

To regenerate the corpus from the raw files:

```bash
python3 -I tools/pde_corpus/import_abbasi_pde.py --scrape <saved Ch1.S4.htm>   # lookup.tsv
python3 -I tools/pde_corpus/import_abbasi_pde.py                               # DE_PDE_*.m
```
