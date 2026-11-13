# 257: compact integrands, large integrals

There are 1000 indefinite integrals in this benchmark. Each integrand is compact,
and its antiderivative is **at least twice the integrand's leaf count**. The
corpus targets Mathilda's decision procedures: Bronstein, Trager, and the
Risch–Norman heuristic ParallelMixedTower. Those methods run against FriCAS,
Maxima, sympy and Wolfram. Cases like these favour a decision procedure, because
the answer is a sum of logarithms, ArcTans and radical terms that a table or a
pattern matcher must hit exactly.

The harness, the verification and the verdicts are those of
[benchmark 256](../256-integration-shootout/README.md). Only the corpus differs.

```bash
make bench-integrate-compact                            # everything installed, ~2 h at 6 jobs
python3 tools/integration_shootout/shootout.py --bench benchmarks/257-compact-integrals \
    --quick --systems mathilda,fricas --check           # the nightly ratchet (ctest integration_compact_stress)
python3 tools/integration_shootout/verify_compact.py    # re-check the corpus itself
```

## The corpus

| class | cases | methods run | integrands |
|---|---|---|---|
| R | 100 | Automatic, PMT, RischTranscendental | `p/q`, with q an irreducible cubic, quartic or quintic, or with repeated factors |
| T | 150 | Automatic, PMT, RischTranscendental | `θ' R(θ)` and `θ' Log[poly(θ)]`, where θ is a Log or Exp of a polynomial in x |
| A | 300 | Automatic, PMT, Trager | `x^k g^e / q(x)` with g a binomial or trinomial and e ∈ {±1/2, ±1/3, 2/3, ±1/4, ±3/2, 1/6} |
| M | 250 | Automatic, PMT, Bronstein | a radical on top of the tower: `θ' g(θ)^e R(θ)` |
| B | 200 | Automatic, PMT, Bronstein | a radical at the bottom: `Log[lin] g(x)^e`, `Log[x + Sqrt[x^2 + a]]`, `E^Sqrt[x + a] (x + a)^k`, … |

Classes are assigned by **content**: a radical with no Log or Exp is A, a kernel
with no radical is T, and both together are M or B. Many integrands drawn by
the B grammar contain no kernel, so they are class A (50 of the kept 300).

## How it was built: search, not construction

The obvious construction picks an antiderivative F and differentiates it. It does
not work: D[F] almost never shrinks to half of F's size for algebraic or mixed F.
Of 600 attempts for class A, 1 qualified. So the generator **searches over compact
integrands** instead. The search runs in `tools/integration_shootout/gen_compact.m`,
driven by `gen_compact.py`:

1. Draw f from the class grammar. Shapes that carry the monomial's derivative as a
   factor are favoured, so that integrable integrands are common. Then put f in
   **compact form**, the smallest leaf count among `Together`, the canonical form
   over K(Y) with `Y^n = b` (the denominator rationalised modulo `Y^n − b`),
   `Factor`, and `Simplify` (2 s).
2. Ask Mathilda's `Integrate` (Automatic, 30 s). The result is then handled as
   follows:
   - **Accept** f with F as its reference when F has no special functions,
     `PossibleZeroQ[D[F] − f]` holds, and `LeafCount[compactForm[F]] ≥
     2 LeafCount[f]`. F is measured in compact form, so a padded F cannot pass.
   - **Candidate**: Mathilda neither closes nor certifies f. These go to FriCAS.
     If FriCAS's answer passes the harness's numeric derivative check, has no
     special functions, and is at least twice the leaf count, it becomes a case
     with F = None and family `c<class>-fricas`. So the corpus is not limited to
     what Mathilda can already integrate. **Wolfram is never consulted.**
   - **Drop** when Mathilda certifies f non-elementary, or closes it with a small
     integral.
3. Per class, keep the first quota of accepted cases in (seed, attempt) order and
   remove duplicate integrands.

Each case also has three rational evaluation points, as in 256, and is labelled
`expect = elem`.

### Runs (`gen_stats.json`)

| class | grammar, seeds × per-run quota | attempts | accepted | FriCAS candidates → accepted | kept |
|---|---|---|---|---|---|
| R | R, 1 × 100 | 133 | 100 | 0 → 0 | 100 |
| T | T, 2 × 75 + 1 × 20 | 1171 | 170 | 33 → 0 | 150 |
| A | A, 4 × 75 | 1280 | 300 | 169 → 0 | 300 (250 from A's grammar, 50 from B's) |
| M | M, 4 × 63 | 1439 | 252 | 340 → 12 | 250 of 250 |
| B | B, 6 × 34; then B2, 6 × 30 | 7863 + 1686 | 184 + 180 | 410 → 0 | 200 of 230 |

Some notes on these figures:
- The B grammar drew mostly kernel-free integrands, and many non-elementary ones.
  Its seeds 1–6 kept only 53 genuine radical-and-kernel cases across 7863 attempts.
  Seeds 7–12 therefore use **grammar B2**. B2 keeps to genus-0 radicands (linear
  or quadratic) times a Log or Exp kernel, the shapes that the component method
  and the genus-0 substitution decide. B2 also requires both a radical and a
  kernel.
- `gen_compact.py --extend` adds new seeds to an existing corpus without
  rerunning the old ones, and `--rerun CLS:SEED,...` redoes chosen runs. Both
  read `gen_pool.tsv`, which holds every accepted case before the per-class quota
  cut, so the result is the same as a run from scratch.
- **Reruns.** The first corpus had 9 cases labelled elem whose reference contained
  `Erfc` or `HypergeometricPFQ`, because the generator's special-function
  blocklist missed both. The shootout reported them as FALSE_CERTs from Trager,
  Bronstein and PMT, and those certificates are correct: `x (1 + x^6)^(1/6)` fails
  Chebyshev's criterion, and `Sqrt[2 + x] E^-x` needs Erf. Acceptance now uses an
  **allowlist** of elementary heads (`elementaryQ`), and `verify_compact.py`
  checks it. Seeds A:1, A:3 and B:2–4 were rerun under the fix, and B:7–12 were rerun once more to rebuild the pool. FriCAS candidate
  counts include the reruns' candidates a second time.
- Most candidates are probably non-elementary: neither Mathilda nor FriCAS closes
  them.
- Approximate leaf ratios of F to f, from token counts, as median and maximum:
  R 2.3 and 18; T 2.7 and 9; A 5.1 and 60; M 6.2 and 78; B 6.7 and 48.

### Reproducing

```bash
python3 tools/integration_shootout/gen_compact.py --jobs 12   # ~5 h on 16 cores; B dominates
python3 tools/integration_shootout/verify_compact.py          # 0 failures
python3 tools/integration_shootout/test_shootout.py           # every case translates to every system
```

The output is deterministic for a given Mathilda: it depends only on the seeds
and on what Automatic can close within 30 s. A later Mathilda that closes more
integrands draws a different corpus. That is why the corpus is frozen here
rather than regenerated.

The results are in `REPORT.md` and `results/`.
