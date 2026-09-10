# DSolve test status

Cross-session dashboard for `DSolve` development against Nasser Abbasi's
*Solving ODEs* corpora (<https://12000.org/my_notes/solving_ODE/current_version/>).
Each section there lists ODEs that Maple **and** Mathematica solve; we track how
many `DSolve` currently solves and drive that number to parity, one general
method-wave at a time (see `../DSOLVE_PLAN.md`, milestones M15+).

The point of this folder: **pick up DSolve development between sessions with a
clear, self-contained context** — the corpora, the self-verifying harness, and a
living scoreboard all in one place.

## Contents

| File | Role |
|---|---|
| `DE_examples_2.m` | Section **2.1.2** corpus — 1204 ODE records (1000 scalar + 204 systems). |
| `DE_examples_3.m` | Section **2.1.3** corpus. |
| `DE_examples_221.m` | Section **2.2.1** corpus — 100 elementary ODEs (Table 2.19, Problems 1–100), 63 IVPs. |
| `DE_examples_222.m` | Section **2.2.2** corpus — 100 elementary ODEs (Table 2.19, Problems 101–200), 9 IVPs. |
| `DE_examples_223.m` | Section **2.2.3** corpus — 100 elementary ODEs (Table 2.19, Problems 201–300), 35 IVPs. |
| `DE_examples_224.m` | Section **2.2.4** corpus — 100 elementary ODEs (Table 2.19, Problems 301–400), 24 IVPs. |
| `DE_examples_225.m` | Section **2.2.5** corpus — 100 series-heavy ODEs (Table 2.19, Problems 401–500), 15 IVPs. |
| `DE_examples_226.m` | Section **2.2.6** corpus — 100 ODEs (Table 2.29, Edwards & Penney, Problems 501–600): 74 scalar (47 IVP, forced linear with general f(t) / DiracDelta impulses) + 26 systems. |
| `DE_examples_227.m` | Section **2.2.7** corpus — 100 ODEs (Table 2.31, Problems 601–700): 50 scalar (23 IVP, elementary first-order) + 50 systems (25 2-D / 18 3-D / 7 4-D). The first half-systems section. |
| `DE_examples_228.m` | Section **2.2.8** corpus — 100 scalar first-order ODEs (Table 2.33, Edwards & Penney, Problems 701–800), 20 IVPs: separable / linear / homogeneous (class A/G/C) / exact / Bernoulli / quadrature. Elementary, same territory as §2.2.1–§2.2.4. |
| `DE_examples_229.m` | Section **2.2.9** corpus — 100 scalar ODEs (Edwards & Penney, Problems 801–900), 35 IVPs: dominated by second-order linear constant-coefficient (homogeneous + nonhomogeneous / undetermined-coeff / variation of parameters) and Euler/Emden–Fowler, plus 3 complex-coefficient and 7 `x(t)`-dependent-variable forms. Fully solved out of the box (100/100). |
| `DE_examples_2210.m` | Section **2.2.10** corpus — 100 ODEs (Edwards & Penney, Problems 901–1000): 56 scalar (15 IVP) + 44 first-order linear systems (2×2/3×3/4×4). Constant-coefficient linear of every order, Euler–Cauchy (incl. Bessel and Gegenbauer/Legendre-type), and constant-matrix systems (42 homogeneous, 2 forced). Solved 100/100 after two M30 fixes (irrational-spectrum forced systems; Kovacic inhomogeneous closure). |
| `DE_examples_2211.m` | Section **2.2.11** corpus — 100 ODEs (Edwards & Penney, Problems 1001–1100): 59 scalar (13 IVP) + 41 first-order constant-coefficient linear systems (2×2 … 6×6, defective/repeated/complex spectra). Scalar half is quadrature / separable / first-order-linear and second-order linear (const-coeff, exact, Euler, Gegenbauer, Emden–Fowler, Liénard, Airy). Solved 100/100 after two M31 fixes (an `Integrate` linearity fix — `try_linearity` now distributes a product over a sum factor, `c(g+h)→cg+ch` — which also repairs a user-reported direct-`Integrate` branch-wrong/slow antiderivative on `E^(-9x)(a E^(2x)-b E^(2x))`; and a cancellation-robust corpus verifier for large-eigenvalue systems). |
| `DE_examples_2212.m` | Section **2.2.12** corpus — 100 scalar first-order ODEs (Edwards & Penney, Problems 1101–1200), 38 IVPs: separable / linear / quadrature / homogeneous / Bernoulli / exact, plus homogeneous-class-A "Abel" rationals and two solvable-for-y/x forms. The first §2.2.x section NOT solved out of the box: went 80/19/**1 FAIL** → 97/3/0 via three M32 root-cause fixes (IVP constant-fitting drops unsatisfiable/Undefined branches — repairing the wrong answer 1147; a Separable implicit twin returning the first integral for non-invertible separations; and an `Integrate` Gaussian→Erf/Ei/PolyLog recognizer that had emitted a literal `x` for every integration variable). Residue 3: solvable-for-y/x (1135, 1200) and Abel-2nd-kind (1157). |
| `DE_examples_2213.m` | Section **2.2.13** corpus — 100 ODEs (Edwards & Penney, Problems 1201–1300), 32 IVPs: a MIX of first-order (exact / linear / separable / homogeneous / Abel / symmetry) and 2nd-order linear (46 reducible-μ, 6 Euler–Cauchy "Emden–Fowler", plus const-coeff). Baseline 91/100 (0 FAIL) → 99/100 via five M33 shared-substrate fixes (exact Path-1/Path-2 potential + denominator-clearing + robust `mu(y)`; implicit-verify `Together`; separable fast numeric pre-filter; first-order IVP undecided-fit fall-through). Residue 1: `1203` Abel-2nd-kind (M13-deferred). |
| `DE_examples_2214.m` | Section **2.2.14** corpus — 100 ODEs (Boyce & DiPrima, Problems 1301–1400): 99 scalar (25 IVP) + 1 system, 2nd-order-linear dominated. Solved 99/100 after four M34 shared-substrate fixes (VoP verify short-circuit / robust variation of parameters / bounded-Kovacic complex-pole gate / SeriesData + IC-point Frobenius series). Residue 1: forced Duffing `1360`. |
| `DE_examples_2215.m` | Section **2.2.15** corpus — 100 ODEs (Boyce & DiPrima, Problems 1401–1500): 39 scalar (18 IVP) + 61 systems (53 2×2 + 8 3×3 constant-coefficient linear). Scalars: high-order constant-coefficient linear + **step/piecewise/Heaviside-forced 2nd-order IVPs** (1492–1500, the Laplace-transform chapter). Solved 98/100 by the new `DSolve\`PiecewiseForcing` (interval continuation → verified `Piecewise`; M35) — these nine previously FALSE-passed as inert `Integrate[UnitStep…]`. Residue 2: `1463` (4th-order transcendental), `1469` (3rd-order variable-coeff). |
| `test_dsolve_corpus.c` | Fork-per-case runner (compiled via `tests/CMakeLists.txt`). |
| `dsolve_corpus_prelude.m` | Self-verifier: runs `DSolve` under `TimeConstrained` and numerically back-substitutes each branch. |
| `STATUS.md` | **The scoreboard** — per-section, per-bucket solve counts + wave history. Update after every wave. |
| `reports/` | Per-section bucketed gap reports (regenerated from the run TSV). |

A corpus record is `{"label", equation(s), function(s), indVar, "MapleClassif", sympy?}`.
**Systems are VERIFIED like scalars (M27):** a system's equation slot is a `List` of
equations and its function slot a `List` of dependent functions; a solution branch
`{x->Function[…], y->Function[…], …}` is back-substituted into every equation exactly
like a scalar ODE. (Through M26 systems were skipped — the "scalar-first" phase.)

Verdict codes: `0 PASS` (verified closed form), `1 FAIL` (a branch is
demonstrably nonzero — a wrong answer), `2 UNEVAL` (declined / timed out / `{}`),
`3 SKIP` (reserved for a record that is not a solvable ODE/system).

## Run a section (progress dashboard)

```bash
cd tests/build && cmake .. && make dsolve_corpus_tests -j
ctest -R dsolve_corpus_2_1_2_tests --output-on-failure     # gated at STATUS.md baseline

# Manual full run, capturing the per-case TSV for the bucket report:
DSOLVE_CORPUS_BASELINE=100000 ./dsolve_corpus_tests \
    ../../DSolve_test_status/DE_examples_2.m \
    ../../DSolve_test_status/dsolve_corpus_prelude.m \
    > /tmp/corpus.tsv 2>/tmp/corpus.err
python3 ../../tools/dsolve_corpus_report.py /tmp/corpus.tsv \
    > ../../DSolve_test_status/reports/2.1.2.md
```

The `ctest` gate fails only when non-PASS scalar cases exceed the per-section
baseline (argv[3] in `tests/CMakeLists.txt`). Each landed method-wave must
**lower** that baseline; a newly-failing case trips the test.

## Regenerate / add a corpus

The site 403s automated fetch, so save the section HTML manually, then convert:

```bash
curl -sL -A "Mozilla/5.0" \
  'https://12000.org/my_notes/solving_ODE/current_version/indexsubsectionN.htm' \
  -o /tmp/sectionN.html
python3 tools/latex_ode_to_mathilda.py /tmp/sectionN.html \
  DSolve_test_status/DE_examples_K.m --label 2.1.K
```

Then add a `dsolve_corpus_<section>_tests` entry in `tests/CMakeLists.txt` and a
section block in `STATUS.md`. The converter is section-agnostic (auto-selects the
problems table and reads its column layout from the header row, which varies between
sections).

**Initial-value problems.** When the source rows carry initial conditions (`y(0)=3`,
`y'(0)=10`, symbolic `y(a)=b`), the converter emits the record's equation slot as the
DSolve-native list `{ode, ic1, …}` (the function slot stays a bare symbol, so the scalar
harness still handles it). The prelude solves `DSolve[{ode, ics}, y, x]` and verifies the
ODE residual **and every initial condition**; an IVP whose solved branch still carries a
generated constant `C[k]` (the general solution with the IC unfitted) scores UNEVAL, not
PASS. §2.2.1 is the first IVP-carrying section; §2.1.2 has no ICs (bare equations).

## After a method-wave

1. Re-run the affected section(s), regenerate `reports/<section>.md`.
2. Update the counts and the wave-history line in `STATUS.md`.
3. Lower the section baseline (argv[3]) in `tests/CMakeLists.txt` to the new
   non-PASS count.
