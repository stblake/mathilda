# FLINT Acceleration Sweep — all builtins (2026-07-01)

**Method.** Five parallel domain audits cross-referenced every relevant builtin
against the installed FLINT 3.6.0 headers (`/usr/local/Cellar/flint/3.6.0/include/flint/`).
FLINT is currently wired **only** into the polynomial *extension* engine
(`src/poly/flint_bridge.c`: GCD/Factor/Cancel/Resultant over Q(α)/Q(√d)/parametric
radical fields). Everything else is hand-rolled or GMP/MPFR-based.

444 builtins registered. The high-value opportunities cluster around **two
enabling bridges** that don't yet exist, plus a set of one-line wirings/swaps.

---

## Two enabling infrastructure pieces (unlock most of the wins)

- **`flint_num_bridge` (acb/arb ↔ Expr)** — ~200-300 LOC beside `numeric_complex.c`.
  Expr→acb at chosen precision; acb-ball→Expr honoring Mathilda's "≤53 bits→Real
  else MPFR / zero-Im→real leaf" conventions; working precision from existing
  `numeric_min_inexact_bits` + guard bits. Unlocks the entire special-function tier.
  **Caveat:** validate branch-cut conventions per function against the existing
  hand-coded cuts (LogGamma continuous branch, PolyLog approach-from-below, Ei ±iπ).
- **`flint_mat_bridge` (fmpz_mat/fmpq_mat ↔ Expr)** — small; reuses
  `flint_bridge.c:131/218/222`. Gate on "all entries integer/rational", fall back
  to the symbolic path otherwise. Unlocks the entire exact-linalg tier.

---

## Tier 1 — highest value / lowest effort (do first)

| # | Item | Change | Effort | Payoff |
|---|------|--------|--------|--------|
| 1 | **Plain `PolynomialGCD` / `Resultant` / multivariate `Factor` over Q** | Route to the **already-written, already-tested** `flint_multivariate_gcd` / `flint_multivariate_resultant` / `flint_multivariate_factor_impl` (bridge helpers exist but `builtin_polynomialgcd` at `poly.c:2138` calls only `flint_extension_gcd`). Flip/augment ~2 call sites + WL sign/content normalization. | **S** | **High** — kills subresultant-PRS coefficient blowup; also feeds Cancel/Together/Apart → **directly helps the Goursat inner integral**. |
| 2 | **`Det` → `fmpz_mat_det`/`fmpq_mat_det`** | Current `laplace_det` (`det.c:13`) is **O(n!)**. Needs the mat bridge (or an inline convert). | **S** (+bridge) | **High** — n! → polynomial time; exact, matches WL. |
| 3 | **Zeta / HurwitzZeta / PolyGamma / StieltjesGamma → `acb_dirichlet_*`, `acb_polygamma`** | Collapses the deliberately duplicated Euler-Maclaurin kernels in `zeta.c`/`hurwitzzeta.c`; Riemann-Siegel speed; **adds numeric StieltjesGamma** (does not exist today). Simple scalar args. | **S** (+num bridge) | **High** — speed + dedup + new capability. |

## Tier 2 — high value / medium effort

| # | Item | Notes | Effort | Payoff |
|---|------|-------|--------|--------|
| 4 | **Univariate integer `Factor` → `fmpz_poly_factor`** (van Hoeij/LLL) | Current `facpoly_bz_uni.inc` is int64-bounded Berlekamp-Zassenhaus w/ **exponential recombination** and **silently returns bignum-coeff polys UNFACTORED** (`:581-597`). Fix = speed **+ correctness/capability**. | S–M | **High** |
| 5 | **`Inverse`/`LinearSolve`/`RowReduce`/`MatrixRank`/`NullSpace` → fmpz_mat/fmpq_mat** | Bareiss-over-Expr today; Dixon lifting + modular RREF. RowReduce cluster shares an engine (one bridge accelerates all). Preserve `Method->` surface + WL nullspace/RREF normalization; `LinearSolve` inconsistent case via `can_solve`. | M (+mat bridge) | **High** |
| 6 | **Bessel J/Y/I/K + Airy Ai/Bi → `acb_hypgeom_bessel_*` / `acb_hypgeom_airy`** | Retires ~2,150 lines of the most branch-cut-fragile code (every hand-made parity-fold/force_real/branch refusal). One `acb_hypgeom_airy` call yields Ai,Ai',Bi,Bi'. | M (+num bridge) | **High** (robustness) |
| 7 | **BernoulliB[n] / EulerE[n] → `arith_bernoulli_number` / `arith_euler_number`** | O(n²) growing-rational recurrence w/ arbitrary cap 5000 today; FLINT uses zeta/multimodular. Check EulerE sign vs E₂=−1,E₄=5. | M | **High** for large n |

## Tier 3 — correctness/capability gains

| # | Item | Notes | Effort |
|---|------|-------|--------|
| 8 | **PolyLog / LerchPhi / HypergeometricPFQ → `acb_polylog` / `acb_dirichlet_lerch_phi` / `acb_hypgeom_pfq`+`_2f1`** | Full analytic continuation: today PolyLog uses a fragile −0-Im trick, LerchPhi refuses on-cut & many \|z\|>1, pFq returns NULL for \|z\|≥1. | M (+num bridge) |
| 9 | **PartitionsP → `partitions_fmpz_ui`** | Deletes ~150 lines of hand-rolled MPFR Rademacher (FLINT is the canonical impl of the same paper). Keep pentagonal + PartitionsQ hand-rolled. | M |
| 10 | **New FLINT-native builtins** | `Cyclotomic`, `ChebyshevT/U`, `HermiteH`, `InterpolatingPolynomial` (`fmpz_poly`/`fmpq_poly` one-liners); `StirlingS1/S2`, `BellB` (`arith`); `HermiteDecomposition` (HNF, `fmpz_mat_hnf_transform`). All currently **absent**. | S each |
| 11 | Gamma/LogGamma incomplete, Erf/Erfc/Erfi, ExpIntegralEi, Fresnel(new), Beta-incomplete, Glaisher/Khinchin → acb/arb | Retires duplicated `ecx` toolkits; rigorous cuts. | S–M (+num bridge) |

## Tier 4 — integer factoring force-multiplier

| # | Item | Notes | Effort |
|---|------|-------|--------|
| 12 | **FactorInteger → add `qsieve_factor` tier** | Fills the 40–100 digit balanced-semiprime gap where ECM degrades and only slow hand-rolled Dixon/CFRAC exist. Multiplies through φ/μ/σ/Divisors/orders. | M/L |
| 13 | **Word-size fast paths: EulerPhi/MoebiusMu/DivisorSigma → `n_euler_phi`/`n_moebius_mu`/`n_factor`** | Skips the `FactorInteger`→Expr→`df_factor_mpz` round-trip for small n (the common case). | S |

---

## Explicit NON-wins (leave alone)
- GCD/LCM/ExtendedGCD/JacobiSymbol/PowerMod-int/Factorial/Binomial/Fibonacci/LucasL — thin GMP wrappers on both sides; FLINT = churn, no gain (GMP *is* FLINT's integer backend).
- **Prime/PrimePi** — Mathilda's multi-method Deléglise-Rivat counter exceeds anything FLINT offers (FLINT has no analytic prime counting).
- **GroebnerBasis** — FLINT 3.6 has no Groebner engine.
- Expand/Collect/Decompose/CoefficientList — structural Expr transforms, no FLINT analogue.
- EulerGamma/Catalan constants — MPFR already optimal.

---

## Recommended execution order
1. **Tier 1 #1** (plain PolynomialGCD/Resultant rewire) — no new infra, immediate, helps Goursat.
2. **`flint_mat_bridge`** → Tier 1 #2 (Det) → Tier 2 #5 (solve/rref cluster).
3. **`flint_num_bridge`** → Tier 1 #3 (Dirichlet cluster) → Tier 2 #6 (Bessel/Airy) → Tier 3 #8/#11.
4. Tier 2 #4 (univariate Factor) and #7 (Bernoulli/Euler) — independent, land anytime.
5. Tier 3 #10 new builtins, Tier 4 factoring — as capacity allows.

Each item is independently testable and independently landable; all keep the
`USE_FLINT=0` classical path as fallback (matching the existing policy).
