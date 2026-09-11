# Task: high-fugacity (degenerate) Fermi–Dirac half-line integrals

## Context
`klimanek/Bose-Fermi` (SymPy) is already fully subsumed by `rec_expgeom` in
`src/calculus/integrate_ramanujan.c`. The one real gap: the fugacity gate caps
`|γ'| ≤ 1`, declining the convergent high-fugacity Fermi regime `γ' > 1`
(degenerate electron gas). Fix = relax the gate's upper bound; the closed form
and strip are already general.

## Plan
- [ ] Relax `iv_prove_fugacity` (drop upper `v.hi <= 1.0`; keep strict lower `γ' > -1`)
- [ ] Relax `rec_expgeom` concrete gate (`GreaterEqual[γ',-1]` only, drop `LessEqual[γ',1]`)
- [ ] Rewrite the derivation comment to state the true convergence region `γ' ≥ -1` + analytic continuation
- [ ] Add `test_high_fugacity_fermi()` to `tests/test_integrate_ramanujan.c` (register in main)
- [ ] Keep the divergent-Bose (`Exp[x]-z, z>2`) decline test unchanged
- [ ] Changelog `docs/spec/changelog/2026-09-07.md` + `Mathilda_spec.md` row
- [ ] Extend the definite-integration entry in `docs/spec/builtins/`
- [ ] Build (`make -j`), run REPL checks + NIntegrate cross-check, run test suite

## Review
Done. Single-gate relaxation in `rec_expgeom`/`iv_prove_fugacity`
(`src/calculus/integrate_ramanujan.c`): the fugacity gate is now the integral's
convergence bound `γ' ≥ -1`, not the geometric series' `|γ'| ≤ 1`. This turns on
the high-fugacity / degenerate Fermi regime `γ' > 1`, previously declined:
- `∫₀^∞ 1/(e^x+2) = ½Log[3]` (=0.549306), `∫₀^∞ x/(e^x+2) = -½PolyLog[2,-2]`
  (=0.718373), `∫₀^∞ 1/(e^x+5) = ⅕Log[6]` — all matched `NIntegrate` exactly.
- Symbolic: `∫₀^∞ x^(s-1)/(e^x+z) = -Γ[s]PolyLog[s,-z]/z` (`s>0,z>0`).
- Divergent Bose `γ'<-1` still declines (verified `Exp[x]-z, z>2`).

Verification: `integrate_ramanujan_tests` — new `test_high_fugacity_fermi`
passes; the only 2 FAILs are pre-existing Hypergeometric `1F1`/`2F1` **Mellin**
cases (confirmed identical in a stashed baseline, unrelated to `rec_expgeom`).
`integrate_dispatch_tests` 0 FAILs. `make` clean (gcc-16), `make check-c99` exit 0.

Docs: `docs/spec/builtins/calculus.md` (transform-table row `γ≥-1` + prose),
`docs/spec/changelog/2026-09-07.md`.

Assessment answer to the user: porting `klimanek/Bose-Fermi` wholesale was
redundant — `rec_expgeom` already subsumed all its rows; this closes the one
genuine gap (high-fugacity Fermi). No named `FermiDiracIntegral`/`BoseEinsteinIntegral`
heads added (not in Mathematica; the `PolyLog`/`Zeta` forms are the faithful output).
