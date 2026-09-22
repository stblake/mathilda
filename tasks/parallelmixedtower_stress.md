# ParallelMixedTower stress test vs the sympy reference — findings

Date: 2026-09-22. Harness: `mixed/stress/` (corpus.py, run_stress.py, _pyworker.py, report.py).
Full report: `mixed/stress/stress_report.md`; data: `stress_results.json`, `stress_table.tsv`.

**Setup.** 114 unique on-domain integrands (all pass a build-time sympy→WL numeric
round-trip gate). Same algorithm on two hosts: Mathilda `` Integrate`ParallelMixedTower[f,x] ``
(strict, no cascade) vs sympy `build_tower.integrate_surface` (`_HAVE_RN=False`). Per-engine
60 s solve budget; every returned antiderivative independently re-verified by differentiation
at domain-safe points. Mathilda 0.167 (GCC 16.1.0), sympy 1.14.0.

## Headline

| metric | Mathilda | sympy |
|---|---:|---:|
| solved & verified | **52** | **92** |
| declined | 33 | 11 |
| non-elementary certificate | 13 | 9 |
| timeout (60 s) | 13 | 1 |
| wrong (returned, failed verification) | **3** | 0 |
| error | 0 | 1 |

- **Coverage:** sympy solves 92/114, Mathilda 52/114. Exclusive solves: sympy-only **40**,
  Mathilda-only **0**. On this corpus Mathilda's coverage is a strict subset of sympy's.
- **Speed (52 both-solved-and-both-verified):** geomean Mathilda/sympy = **0.55×**
  (median 0.41×) — Mathilda is ~1.8× faster on shared cases, faster on **41 of 52**. But it
  has a heavy tail: A13 83×, S4 44×, S12 37× slower (SplitSpecials / ArcTan-over-radical).
  Fastest for Mathilda: R5 0.06×, A33 0.15×, A24 0.18×.
- **Robustness:** 13 Mathilda timeouts vs 1 for sympy; every Mathilda timeout is a Charlwood
  case sympy solves (P4, P8, A1, A2, A3, A16, A19, A20, A27, A35, A37, A40).

## Correctness defects in Mathilda's ParallelMixedTower (all reproducible via the raw symbol)

1. **N6 `Sqrt[Log[x]]` → returns `0`.** A non-elementary integrand; PMT returns the integer
   `0` as the antiderivative (D[0]≠f). sympy honestly declines. **Soundness bug** (bogus answer).
   Cascade contained: `Integrate[Sqrt[Log[x]],x]` returns the correct
   `x Sqrt[Log[x]] - (1/2) Sqrt[Pi] Erfi[Sqrt[Log[x]]]` from an earlier stage; the bogus `0`
   is raw-PMT-only.
2. **T2 `1/(x (Log[x]^2+1))` → false non-elementarity certificate.** PMT emits
   `Integrate::nonelem` ("holomorphic remainder … not exact"), but the integral is
   `ArcTan[Log[x]]` (elementary; the cascade returns it correctly via an earlier stage).
   **Soundness bug** (false certificate). This is examples.py R3, which the sympy reference
   solves. Contained in the cascade because DerivativeDivides handles it before PMT — but PMT
   alone is unsound here.
3. **F2 `Tan[Sqrt[x]]/Sqrt[x]` → leaks an internal expression.** PMT returns
   `-1/OptionValue[ParallelMixed`Private`iPIM, {...}, "SpecialExponent"]` — an unresolved
   internal `OptionValue` call. True answer `-2 Log[Cos[Sqrt[x]]]` (cascade returns it).
   **Package bug** (private-context/OptionValue leak into the result).
4. **A11 `x^3 ArcSin[x]/Sqrt[1-x^4]` → wrong antiderivative.** PMT returns an expression whose
   D-check residual is −0.135 at x=1/2 (in-domain). sympy solves it correctly. Cascade does not
   expose the wrong answer but **times out** (>40 s) on `Integrate[x^3 ArcSin[x]/Sqrt[1-x^4],x]`.

**Severity note.** All four defects are on the strict `` Integrate`ParallelMixedTower `` surface.
The shipping cascade masks each one for these specific integrands (T2/F2/N6 solved by an earlier
stage, A11 times out), so a user calling plain `Integrate[...]` does not currently hit them. But
they are genuine soundness/robustness bugs in PMT as a *method*: a false non-elementary
certificate and a bogus `0` are dangerous in a last-resort cascade stage — any integrand only PMT
can reach would inherit the wrong result.

## Coverage gap (40 sympy-only solves)

By category: C(Charlwood) 22, M(m≥3 radicals) 7, R(single radical) 6, T 3, S 1, F 1. On those
40, Mathilda: declined 24, timeout 12, wrong 2, nonelem 2. The gap concentrates in the
radical-heavy Charlwood appendix and the m≥3 Trager-basis cases — consistent with the port's
documented "remaining gap: complex-residue realisation" and deferred Tier-2 `ToNumberField`.

## Reading

- Where both engines succeed, Mathilda's C-hosted interpreter is materially **faster** than
  CPython+sympy (~1.8× geomean), with a few SplitSpecials/ArcTan outliers that are far slower.
- sympy is far **broader** (92 vs 52) and far more **robust** (1 timeout vs 13, 0 wrong vs 3)
  on this corpus. Mathilda's port has real soundness defects (a false certificate, a bogus `0`,
  an internal leak, a wrong antiderivative) that the strict method surface exposes; the shipping
  cascade masks some of them by reaching a correct earlier stage first.
- Caveats: sympy ran with `_HAVE_RN=False`, so a handful of its certificate cases are weaker
  than the full reference; the comparison pins the strict PMT surface, not Mathilda's full
  Integrate cascade.
