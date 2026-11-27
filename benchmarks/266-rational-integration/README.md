# 266 — the rational-function integrator

423 indefinite integrals of rational functions of x, in 17 families. They are run in
Mathilda (`Integrate[f, x]` Automatic and `Method -> "BronsteinRational"`), in the
local reference system (`wolframscript`) and in sympy (`ratint`). The benchmark was
built for the review of `src/calculus/intrat.c` (`tasks/intrat_review.md`).

```bash
make bench-integrate-rational                       # LABEL=after2 BASE=before GRID=grid2; ~15 min at -j4
make bench-integrate-rational-grid                  # all 36 algorithm combinations, ~70 min at -j6
python3.11 benchmarks/266-rational-integration/run.py --label x --systems mathilda --families stress
python3.11 benchmarks/266-rational-integration/report.py --labels before,after
python3.11 benchmarks/266-rational-integration/gen_corpus.py   # regenerate corpus.tsv (deterministic)
```

`run.py` is resumable: a case already in `results/<label>/<system>.tsv` is skipped.
Pass `--fresh` to start over. The reference and sympy columns do not depend on the
Mathilda build, so an after-run can reuse them:
`cp results/before/{wolfram,sympy}.tsv results/after/`.

This is a standalone study, like 256–264, not a `bench-gap` `.m`/`.py` pair. It is
listed in `STANDALONE` in `benchmarks/run_all.py`.

## The corpus (`corpus.tsv`, from `gen_corpus.py`, seed 266)

| family | n | what it stresses |
|---|---|---|
| classical | 45 | textbook integrals (Stewart, Spivak, G&R 2.1) |
| bronstein | 14 | the examples of Bronstein ch. 2: Hermite/Mack/Horowitz 2.2.1–2.3.1; RT/LRT/Czichowski 2.4.1–2.6.1; Rioboo 2.8; 2.5.2; Czichowski's degree-14 |
| sqfree | 40 | squarefree denominators of degree 5–20: `1/(x^n ± 1)`, n = 5..16, and random products of irreducibles |
| nonsqfree | 30 | random products of powers (multiplicity 2–4) of irreducibles |
| multiplicity | 36 | `1/(x^n+1)^k` (n ≤ 5, k ≤ 6), `1/(x^2+1)^k` (k ≤ 12), `1/(x^3+x+1)^k`, `x^m/(x^2+x+1)^k` |
| rootsum | 33 | irreducible resultants of degree 3–13: `1/(x^n+x+1)`, `x/(x^n-x-1)`, random irreducible denominators |
| rioboo | 24 | `D[ArcTan[P/Q]]` for random P, Q, whose naive log→ArcTan conversion jumps at the real zeros of Q, plus Bronstein 2.8 |
| param_generic | 33 | symbolic parameters a, b, c, generic |
| param_special | 20 | a parametric family at the values where its generic answer degenerates (zero discriminant, coinciding roots) |
| alg_qi | 20 | coefficients in Q(I) |
| alg_sqrt | 20 | coefficients in Q(Sqrt[n]) |
| alg_nested | 12 | nested radicals, cube roots, two radicals at once |
| alg_root | 12 | Root-object coefficients, including the Method -> "BronsteinRational" hang repro `r = Root[837/4 + #1^3 &, 1]`, ρ(x) |
| largecoef | 25 | 10- to 60-digit coefficients, denominators of degree 1–9 |
| czichowski | 14 | high-degree squarefree denominators with few distinct residues (2 log S1 + 3 log S2; D[ArcTan[P/Q]] with deg P up to 7): long subresultant chains, tiny lex Gröbner bases (Bronstein §2.6) |
| logtoreal | 20 | residue factors of degree 3–4 with complex roots: quartics irreducible over Q that split over Q(Sqrt[d]), random cubics and quartics (Bronstein §2.8) |
| stress | 25 | `∏(x^2+k)`, Swinnerton-Dyer, `1/(x^n+1)` for n = 20..30, Mignotte-like, irreducible degree 12–20, Hermite on degree 24 |

## What is measured

| column | |
|---|---|
| status | **PASS** only when Mathilda's own `PossibleZeroQ[D[F, x] - f]` is `True`. FAIL: it is not. UNEVAL: the answer still contains `Integrate`. TIMEOUT: `TimeConstrained` (30 s). KILLED_WALL / KILLED_RSS: the driver killed the process (2×timeout + 75 s wall, 4 GB RSS). VTIMEOUT: the verification itself took more than 60 s |
| JUMP | for a real, parameter-free integrand: Re F is sampled on [-4, 4] with step 1/25, and an interval is flagged when Re F changes by more than h·f(mid) + 0.3. Intervals within h of a real pole of f, or where \|f\| > 20, are skipped. A JUMP is a spurious discontinuity, the Rioboo concern of Bronstein §2.8, which differentiation cannot see |
| time | `AbsoluteTiming` of the integral alone, in a fresh Mathilda process (one process per case and method, so a crash or hang costs one case). The reference runs all cases in one warm kernel |
| size | `LeafCount`. sympy's is counted with the same rules (a Rational is 3 leaves, `I` is 3, an n-ary node is 1 plus its children) |

The reference column is a timing and size comparison only, **never a correctness
oracle**. Nothing Mathilda reports depends on what the reference system returns.
sympy cannot take Root-object coefficients (`alg_root`: n/a).

## Labels

| label | build |
|---|---|
| `before` | the integrator before the review (HEAD + main's uncommitted patch, 2026-11-27); the reference and sympy columns live here and are copied |
| `after` | the first pass (native subresultant chain, BigInt sign, compact Hermite part, Log content), 389 cases |
| `after2` | the second pass (all Bronstein ch. 2 algorithms, measured defaults), 423 cases |

## Files

| file | |
|---|---|
| `gen_corpus.py` | deterministic generator (python3.11 with sympy) |
| `corpus.tsv` | `id, family, integrand (InputForm), sympy form, note` |
| `run.py` | the three runners |
| `report.py` | `REPORT.md` from `results/<labels>` |
| `results/grid2/grid.tsv` | the algorithm grid: id, configuration (Hermite\|LogPart\|RealForm), status, time, leafcount, jump; `results/grid/` is the first grid, before the size caps (see `tasks/intrat_review.md`, S1) |
| `results/before/`, `results/after/`, `results/after2/` | raw rows: `mathilda.tsv` (id, method, status, time, leafcount, jump, #messages, InputForm), `wolfram.tsv`, `sympy.tsv`, `meta.txt` |
| `REPORT.md` | per-family tables, before/after, worst offenders (generated) |
