# Factor shootout

Measures Mathilda's polynomial `Factor` against a commercial reference engine (the
`reference` column), PARI/GP (`factor`, `factor(Mod(1,p)*f)`, `nffactor`), Maxima (`factor`,
`modulus:p`, `gfactor`, `factor(f, minpoly)`) and SymPy (`factor_list`, optional). The corpus
is generated with a fixed seed. Every answer is verified independently of every engine.

This is a standalone study, like `98-factorinteger-shootout`. It is **not** part of the weekly
`benchmarks/run_all.py` gap gate: the folder is listed in `run_all.py`'s `STANDALONE` set, so
`make bench-gap` skips it.

    python3 gen_corpus.py                        # writes corpus.json (seeded; needs gp for Swinnerton-Dyer)
    python3 run.py --tag baseline                # all engines, all cases
    python3 run.py --tag dev --engines mathilda --only cyclotomic,mv-dense
    python3 run.py --tag dev --cases mv-power-3 --timeout 60
    python3 report.py baseline > REPORT.md

The scripts re-exec themselves under `/usr/local/bin/python3.11` when the default `python3`
has no sympy. Set `MATHILDA=/path/to/binary` to time a different build.

## Corpus

79 cases in 17 families, every input **expanded** (the harness factors the multiplied-out
product).

- **Univariate over Z:** `cyclotomic` (`x^n - 1` for n up to 2310), `swinnerton-dyer` (degree
  16/32/64, computed with PARI by iterated resultants; a recombination stress),
  `uni-dense` (products of random dense factors), `uni-bigcoef` (factor coefficients of 100 to
  1000 bits), `uni-multiplicity`, `uni-sparse` (degree up to 7500).
- **Multivariate over Q:** `mv-dense` and `mv-sparse` (2 to 6 variables), `mv-nonmonic`
  (leading coefficients that are polynomials in the other variables), `mv-irreducible`,
  `mv-power` (`(x+y+z)^n ± 1`), `mv-determinant` (Vandermonde, Cauchy numerators).
- `rational`: rational functions with a common factor.
- `modulus-uni`, `modulus-mv`: p = 2, 3, 10007 and 2^60 - 93.
- `gaussian`: `GaussianIntegers -> True`, univariate and multivariate.
- `extension`: `Extension -> Sqrt[2]` and `{Sqrt[2], Sqrt[3]}`, univariate and multivariate.

## Timing

- Every (engine, case) pair runs in its own process, so a crash or a timeout is recorded
  against that case alone.
- Timing happens inside the engine. One untimed warm-up is followed by the best of 3. When
  the warm-up exceeds 5 s, its own time is reported instead.
- The reference engine runs `ClearSystemCache[]` and Mathilda reassigns its result before
  every run, so no run is answered from a memo.
- PARI's clock has a resolution of 1 ms, so PARI repeats a case under 20 ms 50 times per run.
- The per-case limit (`--timeout`, default 60 s) bounds the warm-up inside the engine
  (`TimeConstrained`, `alarm`). The process is killed after 1.5 × the limit + 30 s. Either
  one is reported as `TIMEOUT`.

## Verification

All verification is done in Python, with exact arithmetic, and depends on no engine.

1. **Reconstruction.** The product of `factor^multiplicity` must equal the input exactly:
   - modulo p for `Modulus`;
   - modulo the minimal polynomials of Sqrt[2], Sqrt[3], I, or PARI's and Maxima's field
     generator `t`, for the algebraic cases;
   - by cross-multiplication for rational functions.
2. **Irreducible count.** Each engine's count of non-constant factors is recorded (the sum of
   |multiplicity|, after merging equal factors). Every verified count is a lower bound on the
   true count, so report.py takes the largest one as the truth; PARI and the reference engine
   are always among them. An answer with a lower count is under-factored. Over-factoring
   cannot pass reconstruction.

A failure of either check is `CHECK-FAIL`, and its timing is discarded. Other statuses:

| status | meaning |
|---|---|
| `UNSUPPORTED` | Mathilda or the reference returned `Factor[...]` unevaluated |
| `SKIP` | the engine has no such operation (see the caveats below) |
| `TIMEOUT` | the per-case limit was exceeded |
| `ERROR` | the engine crashed, or its answer could not be read |

## Engine caveats

- **PARI.** `nffactor` is univariate only, so multivariate `gaussian` and `extension` cases are
  `SKIP`. PARI drops the content and the leading unit, so the harness appends
  `input / factorback(f)` as a constant factor.
- **Maxima.**
  - A 60-bit modulus is `SKIP`, because `factor` under `modulus: 2^60-93` did not finish.
  - The script raises `factor_max_degree`. At its default of 1000, Maxima silently returns
    a polynomial of higher degree unfactored.
  - Its multivariate factoring mod p is under-factored on `modulus-mv-2` and
    `modulus-mv-3`, and those cases are flagged `CHECK-FAIL`.
  - Its residues are balanced (symmetric).
  - It has no in-engine time limit, so a warm-up longer than `--timeout` is recorded as
    `TIMEOUT`.
- **Mathilda.**
  - It ignores `GaussianIntegers -> True`, and a multivariate `Extension`. It returns the
    factorisation over Q instead of leaving `Factor` unevaluated, so those cases show as
    `CHECK-FAIL` (too few factors) rather than `UNSUPPORTED`.
  - It leaves a multivariate `Modulus` unevaluated, so those cases are `UNSUPPORTED`.
- **SymPy.** It does not implement multivariate factoring over GF(p), so those cases are
  `SKIP`. Rational functions go through `cancel` and then `factor_list(..., frac=True)`.

`REPORT.md` holds the latest full run; `results/<tag>-<engine>.json` holds the raw data.
