# C transcendental Risch: non-elementary polynomials in ArcTan are undecided

Recorded 2026-11-16 (v0.6588). Found through a user report:
`Integrate[x^2 ArcTan[Sqrt[x]]^3, x, Method -> "Bronstein"]` took 19 s and reported
`Integrate::incomplete`, while FriCAS decides it in 0.05 s. The Bronstein method now
routes around the gap: `brIntegrateArcTanForm` in `src/internal/bronstein/Pipeline.m`
tries ParallelMixedTower on the real ArcTan form. The C layer is still deficient for
every other caller. Follow-up: BRONSTEIN_PLAN.md, "Phase 2b".

## Reproduce

    ./Mathilda -file tasks/risch_arctan_nonelem/repro.m

Observed at v0.6588 (times in s):

| f | Risch, real form | Risch, `TrigToExp` form | `ElementaryIntegralQ`, real form | ParallelMixedTower |
|---|---|---|---|---|
| `ArcTan[u]`, `u ArcTan[u]`, `u^5 ArcTan[u]` | ELEM 0.02 | ELEM 0.01 | True | ELEM |
| `u ArcTan[u]^2`, `u^5 ArcTan[u]^2` | ELEM 0.04–0.07 | ELEM 0.02 | True | ELEM |
| `ArcTan[u]^2` | undecided 0.03 | undecided 0.15 | `::undec` | CERT 0.03 |
| `ArcTan[u]^3` | undecided 0.04 | undecided 0.27 | `::undec` | CERT 0.03 |
| `u ArcTan[u]^3` | undecided 0.06 | undecided 0.68 | `::undec` | CERT 0.04 |
| `u^5 ArcTan[u]^3` | undecided 0.08 | undecided **11.3** | `::undec` | CERT 0.06 |

## Correct answer

Every row marked CERT is non-elementary. For example, integrating `ArcTan[u]^2` by parts
leaves `-2 ∫ u ArcTan[u]/(1 + u^2) du`. Writing `ArcTan` as logarithms, partial fractions
produce `Log(1 + i u)/(u + i)`, whose integral needs a dilogarithm, and the conjugate
terms do not cancel it. ParallelMixedTower's certificate is Prop. 9.2(b) (no solution at
the proved bounds). Cross-check against FriCAS: `integrate(atan(u)^2, u)` is returned
unevaluated.

## Suspected cause and fix direction

In the primitive polynomial case (Bronstein, *Symbolic Integration I*, §5.8,
IntegratePrimitivePolynomial), each coefficient step needs a limited integration that
may introduce a **new** logarithm: here `∫ u/(1 + u^2) = Log[1 + u^2]/2`, which is not
in the tower {u, ArcTan[u]}. The plan's Phase 2 log lists "limited integration with a
new logarithm" as the one cancellation case left open, and
`tools/check_bronstein_completeness.py:90` admits it. Fix it with the full
LimitedIntegrate of §7.2. A new logarithm that is not a constant multiple of the
monomial's derivative makes the integral non-elementary (Liouville), and that is a
certificate. The `TrigToExp` cost (11 s) suggests a second problem in the two-log Q(i)
tower, which needs profiling. The adversarial review
`docs/notes/RISCH_TRANSCENDENTAL_REVIEW.md` is to confirm both and pin the source lines.
