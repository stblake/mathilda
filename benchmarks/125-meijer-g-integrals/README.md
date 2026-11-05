# 125 — definite integrals by Meijer G functions

`meijer_g_integrals.m` / `meijer_g_integrals.py`: the Meijer G machinery
against sympy and mpmath, on a representative subset of `tests/RoachCorpus.m`
(the test suites of Kelly Roach's 2008 Quantum paper and his ISSAC '96 / '97
galleries).

| rows | Mathilda | baseline |
|---|---|---|
| `roach ...` (12) | `Integrate[f, {x, 0, Infinity}, Method -> "Roach"]` | `sympy.integrate(f, (x, 0, oo), meijerg=True)` (meijerint) |
| `hyperexpand ...` (4) | `FunctionExpand[HypergeometricPFQ[...]]` | `sympy.hyperexpand` |
| `N ...` (2) | `N[..., 30]` of a 3F2 at 1 and a MeijerG | `mpmath.hyper`, `mpmath.meijerg` |

Run by hand from this folder (cwd matters for the `.m` side):

```bash
../../Mathilda -file meijer_g_integrals.m
/usr/local/bin/python3.11 meijer_g_integrals.py
```

or as part of the job: `python3 benchmarks/run_all.py --only 125`.

The correctness side of the same material is the nightly ctest
`roach_corpus_tests` (`tests/test_roach_corpus.c`), scored against mpmath
references from `tools/roach_corpus_refs.py`.
