# Plan: Trager's algorithm (simple radicals) and Kauers' heuristic as `.m` Integrate methods


## Context

`packages/trager/` is an implementation (in the same language Mathilda uses) of Trager's 1984 algorithm for integrating algebraic functions over a simple radical extension `y^n = g(x)`, `g ∈ Q[x]`. It was written with Claude 4.6.

- It loads in Mathilda today and closes simple cases such as `1/y` on `y² = x² + 1`.
- It is not at Mathilda's testing, performance or soundness standard.

We will port it into Mathilda as native `.m` code under `src/internal/trager/`. It is exposed as two methods:

- **`Integrate[f, x, Method -> "Trager"]`**: the complete algorithm. Divisor construction, principal generators, and a torsion order decided by reduction mod p (Trager Ch 6 §2–3).
- **`Integrate[f, x, Method -> "TragerKauers"]`**: the same Phases 0–4 (validation through residues), with Kauers' 2008 Gröbner heuristic for the logarithmic part. It is fast and incomplete; any residue factor it cannot close makes the method decline.

**Excluded:**
- All Miller work (`MillerKauersLogTerms.m`, `IdealDecomposition.m`) and all Schultz work (`Schultz*.m`). Nothing from them is ported.
- The only helper Kauers borrowed from Miller, `millerBlockOrder`, is replaced by Mathilda's elimination-order `GroebnerBasis`.
- The package's `"Auto"` cascade, which runs Trager, then Miller, then Kauers, is dropped.

**Not in the Automatic cascade.** Neither method joins Integrate's Automatic cascade in this work. They are explicit methods only.

**Decisions already made by the user:**
1. The torsion order is bounded by reduction mod p (faithful to the thesis), not by the package's `k ≤ 12` search.
2. The two log-part engines are separate methods.

**References:**
- Trager thesis. Printed page = PDF page − 1.
- Bronstein's tutorial, §2.1–2.5 and examples 3–6.
- Kauers, ISSAC 2008. This is not in either PDF; the algorithm is taken from `packages/trager/src/KauersLogTerms.m` and must be checked against the paper.

## Ground truth from exploration

**Package review.** The pipeline is Phase 0 (validate, reduce, genus) → 1 (integral basis `w_i = y^i/d_i`, `d_i = Π p_j^⌊i e_j/n⌋`) → 2 (shift x = a + 1/z) → 3 (Hermite, diagonal congruences) → 4 (double-resultant residues and a Q-basis) → 5 (divisors / HNF over K[z] / normalize at ∞ / principal generator / torsion) → 6 (reassemble, LogToArcTan). Bugs to fix during the port:

| # | Bug | Where | Fix |
|---|---|---|---|
| B1 | Branch-place power raised to `n` rather than the ramification index `n/gcd(n,e_j)` | `Divisor.m:286` | Use the ramification index; per Trager pp.62–63 there is no residue at a branch place when all `e_j < n`, so assert this |
| B2 | `(-8)^(1/3)` misread as non-rational | `Reduce.m` | Use exact root extraction (`IntegerRoot`/`Power` on Rational with a sign check) |
| B3 | `normalizeAtInfinity` can produce fractional powers of z; it also uses `/. Infinity -> 0` | `Normalize.m` | Require integer exponents (scale columns by `n`) and handle infinite orders explicitly |
| B4 | A `$Unknown` verifier result is accepted as a generator | `PrincipalGen.m`, `Torsion.m` | Accept only rows integral at ∞ (the Ch 6 §1 theorem is exact, so no heuristic verifier is needed); check the log-derivative identity at the end |
| B5 | Places are enumerated over Q only | `PlaceOrder.m` | Not needed: the exact Ch 6 §1 test replaces place enumeration, so `PlaceOrder.m` is dropped |
| B6 | `residueQBasis` silently falls back to an identity basis when `ToNumberField` fails | `Residues.m` | Make this an explicit decline, never a silent identity |
| B7 | Kauers uses `First[G]` as the eliminant | `KauersLogTerms.m` | Select the element free of x and y |
| B8 | Kauers builds `Root[Function[t, gg], k]` for parametric `gg` | `KauersLogTerms.m` | Use `RootSum` for irreducible factors of degree ≥ 2 |
| B9 | Rule A is applied to a fixed point before rule B, never interleaved | `LogToArcTan.m` | Apply the two rules to a joint fixed point |
| B10 | `guard` is redefined on every call | `IntegrateTrager.m` | Make it a module-level helper |
| B11 | Surface form does not recognise Mathilda's `Power[_, 1/2]`; `x^(1/2)` with `x^(1/4)` is rejected | `Surface.m` | Unify radicals with the same radicand to the LCM of their indices |
| B12 | `zeroQ` runs `RootReduce` everywhere; `Simplify` sits in hot loops (multiplicity loop, `afValueAt`) | `Common.m`, `Residues.m` | Rational paths use `Together`/`PossibleZeroQ`; `RootReduce` only when there are algebraic coefficients; `Simplify` removed from loops |
| B13 | Multiplicity via repeated `D` + `Simplify` | `Residues.m` | Read multiplicities from `FactorSquareFreeList` |
| B14 | `afNorm` (`Det`) recomputed in the lucky-integer loop; A grows as A^k over torsion steps | `Divisor.m` | Compute the norm once via the resultant; reduce h mod A after each step |

**Mathilda state (v0.6081).**
- Present and working: Catch/Throw, Reap/Sow, Association/Lookup, ToNumberField/AlgebraicNumber, MinimalPolynomial, RootReduce, Resultant, `Factor`/`FactorList` with Extension, FactorSquareFreeList, RowReduce/NullSpace/LinearSolve, OptionsPattern/OptionValue, Check/Quiet/Message, ConstantArray, MapThread, UnitVector, `PolynomialGCD[…, Extension -> Automatic]`, `PolynomialGCD[…, Modulus -> p]`, and GroebnerBasis (lex, elimination).
- Missing or diverging:
  - `Failure` and `FailureQ`.
  - `PolynomialExtendedGCD[…, Extension -> …]` stays unevaluated. Over algebraic coefficients it works without the option, so the port never passes it.
  - Matrix `MonomialOrder` is supported only in the plain form, so Kauers uses lexicographic elimination order.
- `mathilda-compat.m` is now obsolete: every shim hits `Protected`.

**Integrate glue** (C, minimal; `src/calculus/integrate.c`):
- `IntegrateMethod` enum :1340, `method_from_string` :1379, `Integrate::method` error string :2439, dispatch switch :2647–2722, `integrate_init` :2768, Integrate docstring method list :2818–2832.
- Pattern to copy: ParallelMixedTower. `pmt_lazy_load` :929 (`mathilda_load_module`), the thin builtin `builtin_integrate_pmt` (~:960, with `tc_async_defer_push` and `mth_msg_suppress_push`), and `try_parallelmixedtower` :1109 (decline rules).
- Definite integrals pass through to the indefinite method for free.
- Package convention: see `src/internal/mixed/ParallelMixed.m` (:60, :107, :3718). It uses `BeginPackage[..., {"Mathilda`"}]`, then `Begin["`Private`"]`, then `LoadModule["dir/X.m"]` for siblings **inside** Private, because short private *functions* do not resolve through `$ContextPath`.
- `makefile` install copies `src/internal/**/*.m` automatically. Add a comment at :676.

**FriCAS** is at `/opt/local/bin/fricas`. A non-interactive pipe works: `printf '…\n)quit\n' | fricas -nosman`. The pattern to reuse is `mixed/cas_compare.py:22–60`, which does `)set output algebra off`, `integrate(f,x)` and `unparse(r::InputForm)`.

## Design

### File layout: `src/internal/trager/`

Each file has a header block (purpose, thesis section and page refs) and a `::usage`-style docstring comment on **every** function: signature, inputs, output, invariants, thesis reference, complexity notes.

| File | Contents (ported from) | Thesis |
|---|---|---|
| `Trager.m` | `BeginPackage["Trager`", {"Mathilda`"}]`; usage for the public entry points `Integrate`Trager`, `Integrate`TragerKauers`; `Begin["`Private`"]`; `LoadModule` of siblings in order | — |
| `Common.m` | failure records (an internal `tragerFail[tag, info]` head, never leaked), predicates, `zeroQ` (tiered), coefficient-field record `{ext, modulus}` threaded through all polynomial ops | — |
| `Validate.m`, `Reduce.m`, `Genus.m` | Phase 0 (B2) | Ch 3 §3 p.39; Ch 2 §4 p.29 |
| `Basis.m`, `Arithmetic.m` | Phase 1. `afNorm` by resultant | Ch 2 §5 pp.30–31 |
| `Shift.m` | Phase 2: `a` avoids roots of the denominator **and** of disc(f) | Ch 4 §1 p.44 |
| `Hermite.m` | Phase 3 diagonal congruences; pole-at-∞ check wired in as a non-elementary certificate | Ch 4 §2–3 pp.44–50; Bronstein eq.12 |
| `Residues.m` | Phase 4: `R(Z)`, square-free multiplicities, Q-basis (B6, B13) | Ch 5 §2 pp.56–59 |
| `Normalize.m` | K[z]-HNF and normalize-at-∞, generic over the field record (Q, Q(α), GF(p)) (B3) | Ch 2 §3 pp.27–29; HNF p.25 |
| `Divisor.m` | descriptors (h, A), lucky integer, add/negate/scale (B1, B14) | Ch 5 §3–4 pp.59–63 |
| `PrincipalGen.m` | exact Ch 6 §1 principality test (B4) | Ch 6 §1 pp.64–65 |
| `TorsionModP.m` | **new.** Good-prime selection (no leading coefficient or denominator vanishes, integral basis and genus preserved, absolutely irreducible, residue minpoly has a root mod p); divisor good-reduction test (p.69 `d₁`, `q` construction, `g + c·p^k` model); order search over GF(p) for M = 1..(√p+1)^{2g}; two primes p and q; mismatched prime-to-pq parts → NonElementary; candidate `N = lcm` then one principality test of `N·D` over Q(α) | Ch 6 §2–3 pp.66–70; Bronstein p.19 |
| `TragerLogTerms.m` | Phase 5 (Trager): one divisor per Q-basis direction, then `TorsionModP` order, then generator | Ch 5–6 |
| `KauersLogTerms.m` | Phase 5 (Kauers): elimination GB for `R(t)`, a principal-power search per factor, `RootSum` output (B7, B8). Undetermined factor → decline | Kauers 2008 |
| `LogToArcTan.m`, `Reassemble.m` | Phase 6 (B9) | — |
| `Pipeline.m` | `tragerIntegrate[f, x, logEngine]`: phase driver; final check `D[result] − f ≡ 0` (exact, numeric fallback at 60 digits); returns a result, `nonElementary[...]`, or `decline[...]` (B10) | — |
| `Surface.m` | radical detection on Mathilda forms (B11); reduces `Integrate`Trager[f, x]` to the curve form | — |

**Dropped from the package:** `PlaceOrder.m` (B5, superseded), `Torsion.m` (replaced by `TorsionModP.m`), the `"Auto"`/`"Miller"`/`"Schultz"` options, the `"Parameters"` mode (parametric integrands decline in v1), and `mathilda-compat.m`/`load.m`.

`packages/trager/` itself is left untouched, as the reference.

**Field-generic polynomial layer.** HNF, normalize-at-∞, extended GCD and norms all take the field record. That lets `TorsionModP` reuse the exact same `Normalize.m` and `PrincipalGen.m` code over GF(p), passing `Modulus -> p`. The first step of that phase confirms which Mathilda polynomial builtins honour `Modulus`: `PolynomialRemainder`, `PolynomialExtendedGCD`, `Together`, `FactorList`, `Resultant`. A builtin that does not gets a small `.m` helper in `Common.m`, not a C change.

### Integrate glue (C, minimal)

In `src/calculus/integrate.c`:
- Add `METHOD_TRAGER` and `METHOD_TRAGER_KAUERS` to the enum, `method_from_string` and the `Integrate::method` string.
- Add `trager_lazy_load()` (`mathilda_load_module("trager/Trager.m")`).
- Add builtins `Integrate`Trager` and `Integrate`TragerKauers`. Give them Protected, docstrings, `tc_async_defer_push` and `mth_msg_suppress_push`, cloned from the PMT pattern.
- Add `try_trager(...)` with the PMT decline rules.
- A proven non-elementary integral emits `Integrate::nonelem` through `mth_message` and returns the integral unevaluated. Any other decline returns unevaluated silently.
- Update the Integrate docstring method list and the error string.
- Add a comment to `makefile` :676.

### Tests (an extensive suite per function)

A C harness in the `tests/test_parallelmixedtower.c` style. It uses `ASSERT_STR_EQ`, not the NDEBUG no-op. It loads `trager/Trager.m` once and calls `Trager`Private`` functions directly.

| Suite | Coverage |
|---|---|
| `test_trager_phase0.c` | `validateInput`, `reduceYDegree`, `reduceIrreducibility` (including `(-8)^(1/3)`, `y⁴ − 4x²`), `computeGenus` (table checked against Riemann–Hurwitz by hand, ≥ 25 cases) |
| `test_trager_basis_arith.c` | basis `d_i` and divisibility invariants; `afTimes`/`afD`/`afInverse`/`afNorm`. Properties: Leibniz, norm multiplicativity, `y·y^{n−1} = g`, inverse round trip, over random small curves |
| `test_trager_shift_hermite.c` | Shift invariants (genus preserved, polynomial radicand, `a` avoids the disc); Hermite: Bronstein Ex. 3 exactly, `D[algPart] + remainder == f`, squarefree remainder denominator, pole-at-∞ certificate |
| `test_trager_residues.c` | `R(Z)` for Bronstein Ex. 4 (`t¹⁶(t²−1)`), Ex. 5 (`729t⁶` → no residues), Ex. 6; Q-basis with Galois orbits; multiplicities |
| `test_trager_normalize_hnf.c` | HNF over Q, Q(√2), GF(p) (triangular, same module: a checked Bezout transform); normalize-at-∞ termination and `Σk = ord det`, integer exponents (B3) |
| `test_trager_divisor_principal.c` | descriptor arithmetic, lucky integer, principality of known principal and non-principal divisors (genus 0 is always principal; elliptic examples) |
| `test_trager_torsion.c` | good-prime selection and rejection; order of a divisor over GF(p) against brute force; known torsion orders 1, 2, 3, 4 (`(x⁴+1)^(−1/4)` order 4, Bronstein Ex. 6 order 3); infinite order → NonElementary (e.g. a genus-1 integrand with a non-torsion residue divisor) |
| `test_trager_logterms.c` | Trager and Kauers engines on shared fixtures; Kauers declines cleanly where it cannot close |
| `test_trager_arctan_reassemble.c` | LogToArcTan rules, the shift inverse, rescale inverses |
| `test_trager_integrate.c` | end-to-end through `Integrate[…, Method -> "Trager" / "TragerKauers"]`: every `TestIntegrate`/`GenusPositive`/`Rescale` case from the package plus the Bronstein examples, each verified by `D` (exact, then a numeric spot check); a `Integrate::nonelem` stderr capture on `1/Sqrt[x³+1]`, `x/Sqrt[1−x³]`; decline (unevaluated, silent) on a compositum or an `Exp[x]` integrand; definite pass-through; `Quiet`/`Check` behaviour |

Register each suite in `tests/CMakeLists.txt` (the explicit-sources convention). They are unit tier, so they run on every push.

**FriCAS stress test** (`nightly`):
- `tools/trager_fricas_stress.py` is modelled on `mixed/cas_compare.py`.
- Corpus `tests/data/trager_stress_corpus.m`: about 300 integrands.
  - Hand-picked families: `x^k/(x^m+c)^(1/n)`, `1/((x−a)·Sqrt[q(x)])`, Chebyshev/Goursat-type pseudo-elliptics, n = 2..6, genus 0..3.
  - Plus a seeded random generator, `R(x)·y^k`, with small coefficients.
- Each integrand runs through Mathilda (both methods, per-case timeout) and FriCAS (`fricas -nosman` pipe, timeout).
- Classification:
  - **AGREE-ELEM**: both close; each verified by numeric differentiation at 3 points to 1e−20.
  - **AGREE-NONELEM**.
  - **MATHILDA-ONLY** / **FRICAS-ONLY** (gap).
  - **DISAGREE**: one proves non-elementary while the other gives a verified antiderivative. This is a hard failure.
  - **WRONG**: a result fails the derivative check. This is a hard failure.
  - **TIMEOUT**.
- Writes `tests/data/trager_stress_report.md` and a JSON report, and records wall time per case for both systems.
- A ctest entry `trager_fricas_stress` (nightly label by name) runs the script. It fails on any WRONG or DISAGREE, or on a FRICAS-ONLY count above a checked-in `BASELINE`, which ratchets.

### Docs and release

- `docs/spec/builtins/calculus.md`: add both methods to the `Method ->` list (~:1509), with a `### Integrate`Trager` subsection describing the algorithm, scope (simple radical, Q coefficients), decline rules and examples.
- Changelog `docs/spec/changelog/2026-10-26.md`.
- Integrate docstring (`integrate.c` ~:2818).
- Book: add an index entry if the calculus chapter mentions methods.
- Bump `src/version.h` +0.0001 and tag each substantive commit.

## Design revision during implementation (2026-10-29)

Recorded here because it changes the module list above.

- **Per-component integration ([T] Appendix A).** `f = Σ a_i(x) yh^i` is split. Each `a_i yh^i` is integrated on its own reduced curve `Y^m = G` (`Reduce.m: componentCurve`), with `Y = yh^i/s`. When gcd(n, exponents) = d > 1, `Y = G^(1/m)` instead, times a locally constant factor `rho`, so no constant `c^(1/d)` is ever adjoined. On such a curve every residue at a branch place is zero ([T] pp.62–63), so the divisor of the residue-α places is exactly the ideal `(N Y − α E', Bp)` of O. **The package's D2/D3 branch hack and the "lucky integer" search (B1, B14) are gone.** `PlaceOrder.m` is not needed (B5).
- **Divisors are fractional ideals** (`Ideal.m`), stored as an HNF matrix over K[x] plus a prime-field denominator. They support product, power, and inverse via Trager Ch 2 §2.
- **Coefficient fields** (`Field.m`) are Q, Q(θ) = Q[t]/(μ), and GF(p), with exact canonical forms. HNF, normalisation at infinity, and the principality test run unchanged over all three, which is how `TorsionModP.m` reuses them.
- **Normalisation at infinity** uses Trager's integer x-adic formulation, with local integral basis `x^(−δ_j) w_j` and `δ_j = ⌈jM/m⌉ − deg d_j`. Exponents are always integers (B3).
- **Exact verification of the log part.** Trager's engine checks `r · Πu = Σ c_j u_j' Π_{k≠j} u_k` in K(x)(Y). Kauers' engine checks via traces over Q[Z]/(q). The final answer is checked with `PossibleZeroQ[D[F] − f]` (60-digit numeric fallback for PossibleZeroQ false negatives on nested radicals).
- **Real form** (`LogToArcTan.m`) takes the formal real part of the whole answer and keeps it only if it re-verifies.
- **Final module list:** Trager, Common, Field, Validate, Reduce, Genus, Basis, Arithmetic, Shift, Hermite, Residues, Normalize, Ideal, PrincipalGen, TorsionModP, TragerLogTerms, KauersLogTerms, LogToArcTan, Pipeline.
- **Test suites:** `test_trager_{field, phase0, basis_arith, shift_hermite, residues, normalize_ideal, torsion, logterms, integrate}.c`, sharing `tests/trager_test_utils.h`.

## Phases (execute continuously; tests gate each phase)

1. Write `TRAGER_PLAN.md` and `tasks/todo.md`. Probe `Modulus` support across the polynomial builtins used.
2. Phases 0–1 port (Common, Validate, Reduce, Genus, Basis, Arithmetic) plus the two suites.
3. Shift and Hermite, plus the suite.
4. Residues, plus the suite.
5. Normalize/HNF (field-generic), Divisor, PrincipalGen, plus the suites.
6. TorsionModP, plus the suite.
7. TragerLogTerms, KauersLogTerms, LogToArcTan, Reassemble, Pipeline, Surface, plus the suites.
8. C glue (methods, builtins, docstrings), plus `test_trager_integrate.c`. Commit; version bump.
9. FriCAS stress harness and corpus; triage every WRONG/DISAGREE to root cause; set `BASELINE`. Commit.
10. Docs, changelog, `tasks/todo.md` review section, lessons.

## Verification

- `make -j` builds clean.
- `make check-c99` and `make check-messages` pass. The C glue routes `nonelem` through `mth_message`.
- `cd tests/build && cmake .. && make -j test_trager_* && ctest -R trager -LE nightly` is all green.
- `ctest -R trager_fricas_stress` gives zero WRONG and zero DISAGREE. Report the agreement table and timings.
- REPL spot checks:
  - `Integrate[1/Sqrt[x^2+1], x, Method -> "Trager"]`, `Integrate[1/(x Sqrt[1-x^3]), x, Method -> "Trager"]` (Bronstein Ex. 6), `Integrate[1/(x^4+1)^(1/4), x, Method -> "Trager"]` (torsion 4), `Integrate[1/Sqrt[x^3+1], x, Method -> "Trager"]` (→ `Integrate::nonelem`, unevaluated).
  - The same set with `"TragerKauers"`.
- `valgrind --leak-check=full` over a short end-to-end script (no growth beyond the startup baseline).
- A broader ctest run over the Integrate suites (`ctest -R integrate`) to show the Automatic cascade is unchanged.
