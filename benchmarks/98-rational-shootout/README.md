# Together / Cancel / Apart shootout

This benchmark times Together, Cancel and Apart in five engines:

- Mathilda;
- a commercial reference engine (the `reference` column);
- Maxima (`ratsimp`, `partfrac`);
- SymPy (`cancel`, `apart`);
- PARI/GP (native reduction of rational functions; no general Apart).

Every answer is verified independently of all five engines.

## Files

This is a standalone study, listed in `run_all.py`'s `STANDALONE` set, so `make bench-gap` skips it.

- `gen_corpus.py` writes `corpus.json` (seed 20261028). It has 98 cases in 26 families:
  - **Together**: linear and quadratic summands, multivariate summands, shared factors, telescoping sums, nested fractions, `log(x)`/`exp(x)` kernels, Gaussian coefficients, `Modulus`.
  - **Cancel**: univariate gcds up to degree 1000, 200-bit coefficients, dense 3-variable and sparse 6-variable gcds, coprime pairs, textbook cases, Gaussian and `Extension -> Automatic` cancellations, `Modulus`.
  - **Apart**: distinct linear factors, high multiplicity, irreducible quadratics and cubics, improper fractions, parameters `(x + a)…`, a parameter only in the numerator and content, Gaussian factors, `Modulus`.

  The text is a neutral syntax. Each emitter rewrites `LG`, `EX`, `I` and `sqrt(2)` for its engine. For Together the corpus also stores the summands, and for Cancel the numerator and denominator. PARI times only the reduction from those pre-parsed pieces.
- `run.py` runs the cases. Each case runs in its own process, with one untimed warm-up and then the best of three runs (only the warm-up when it exceeds 5 s). The script prints `RES` lines, which are verified and merged into `results/<tag>-<engine>.json`.

  ```
  python3 run.py [--engines mathilda,reference,maxima,sympy,pari]
                 [--only family,...] [--cases id,...] [--tag t] [--timeout 60]
  ```

  `MATHILDA=` selects the binary. PARI's clock ticks in milliseconds, so each of its timed runs repeats until it spans at least 100 ms.
- `report.py <tag>` writes the Markdown report: per-family tables, the ratio of Mathilda to the fastest correct engine, the slow cases and the failures.

## Verification

Each engine prints its answer as `num ## den` terms. The two parsers build an AST from that, and the AST is evaluated in three ways:

- **Value.** The terms must sum to the input at three random points modulo a 61-bit prime `p ≡ 1 (mod 8)`, where `I` and `Sqrt[2]` have images, so the map is a ring homomorphism. `log(x)` and `exp(x)` are free variables there. For `Modulus` cases the check is exact, modulo the case's own prime.
- **Lowest terms (Together, Cancel).** The answer must be one fraction `N/D`. For each variable, after a random substitution of the others, the univariate gcd of the images of `N` and `D` must be 1. A true common factor survives the substitution.
- **Apart shape.** Every term is either free of the variable or proper in it. Its denominator, factored over the case's field (Q, Q(params), Q(i), Q(√2) or GF(p)), must be a power of one irreducible, and no two terms may share a base and exponent. `report.py` also compares the term count across engines.

An answer that fails any check is `CHECK-FAIL`, and its time is discarded. The other statuses are `TIMEOUT`, `UNSUPPORTED` (the call came back unevaluated), `SKIP` (the engine has no such operation) and `ERROR`.

## Engine notes

- **Maxima** cannot partial-fraction over Q(i). It reports `CHECK-FAIL` on `ap-gauss`, where the denominators are not split, and on the Gaussian and `Sqrt[2]` Cancel cases, which stay unreduced without `algebraic:true` handling of them.
- **SymPy**'s `cancel` and `apart` have no `Modulus` for multivariate or partial fractions; those cases are `SKIP`.
- **The reference**'s `Apart` has no `Extension` option, so the corpus has no `Extension` Apart family.
