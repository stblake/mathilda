# Fix Together/Cancel/Simplify blow-up on dependent trig/exp generators over Q(√d)

Plan: `/Users/user/.claude/plans/both-these-hang-we-sorted-ritchie.md`

Reproducer (must return 0 fast, currently hangs):
`Simplify[D[-(I/Sqrt[3]) Log[(E^(I x)+2-Sqrt[3])/(E^(I x)+2+Sqrt[3])], x] - 1/(2+Cos[x])]`

Two target integrals (currently hang):
- `Integrate`RischTranscendental[1/(2 + Cos[x]), x]`
- `Integrate`RischTranscendental[1/(Sin[x] + Cos[x]), x]`

## Tasks
- [ ] Read flint_rational_normalize_core + km_* + flint_extension_gcd + nf_detect + expr_has_denominator
- [ ] Part A: algebraic-coeff fallback in flint_rational_normalize_core (route km-forward output through flint_extension_gcd over Q(√d))
- [ ] Part B: extend expr_has_denominator / together activation to fire on reciprocal exp kernels E^(-ix)
- [ ] Part C: extend rat_has_dependent_power_generators to catch function-kernel deps (safety net)
- [ ] Build (force clean rebuild of touched FLINT/rat objects on macOS)
- [ ] Verify: reproducer -> 0; Together[E^(-Ix)+E^(Ix)] combined; both integrals return real forms, no hang
- [ ] Regressions: 1/(5+4Cos), 1/(10+6Cos), 1/(3+2Cos)√5; test_rat.c:212 trig form preserved
- [ ] Run test_integrate_risch_transcendental + rat/together/cancel/simplify test binaries
- [ ] valgrind on the two integrals + reproducer
- [ ] Update docs/spec + changelog

## Review — LANDED 2026-07-22

**Both target integrals now work** (diff-back → 0), plus the whole
`1/(a+b Cos[x])` and `a Sin[x]+b Cos[x]` family (rational AND irrational
discriminants). No case hangs. The hangs were in `Simplify`/`Together`/`Cancel`,
not the integrator — fixed at the root per the user's directive.

Four coordinated fixes:
1. **`src/poly/poly.c`** — radical-generator walkers (`walk_find_radical_base`,
   `walk_gather`, `poly_subst_radical_to_gen`) skip descending into opaque
   analytic-kernel *arguments* (`Log[…/√3]`), unless the kernel *is* the radical
   base (`Sqrt[Log[r]]`). Fixes `Together[Log[…]/Sqrt[3]]` blow-up → unhangs the
   entire `1/(a+b Cos[x])` family incl. √3/√5.
2. **`src/rat.c`** — `flint_gaussian_together` was dead code (`I`→`t` never
   matched `Complex[0,1]`); fixed via `Complex[a_,b_] :> a+b t` and refactored
   into `rat_gaussian_reduce`; added `flint_gaussian_cancel`.
3. **`src/poly/flint_bridge.c`** — km normal form captures `Sqrt[d]`/`Complex`
   constants (gated: kernel present or ≥2 distinct constants) so `Q(i,√d)`
   rationals reduce via `fmpz_mpoly_q`, clearing the `E^(-ix)→g^(-1)` Laurent.
4. **`src/rat.c`** — dependent-trig-algebraic bail in `cancel_recursive` AND
   `builtin_together_compute` (Plus-combine path): ≥2 trig kernels + an algebraic
   constant → leave uncombined (correctness-preserving; the `a Sin+b Cos` case).
   Gated on the algebraic constant so pure-Q trig (`test_rat.c:212`) is untouched.

**Verification.** `test_integrate_risch_transcendental` (incl.
`test_real_trig_reconstruction`), `rat`, `parfrac`, `radical_canonical`,
`eliminate`, cherry/knowles all green. `simplify_tests` (4) and
`risch_residue_split` (2, Fresnel coeff) fail identically on the clean baseline
(pre-existing). valgrind: one 32-byte eval-internal block via
`rat_gaussian_reduce` (function is provably balanced); rest of the leak profile
matches the pre-existing integrator baseline.
