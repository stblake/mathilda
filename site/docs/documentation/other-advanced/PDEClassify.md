# PDEClassify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PDEClassify[eqn, u, {v1, v2}] classifies a second-order linear PDE by the discriminant Δ = B² − 4 A C of its principal part A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2}: "Hyperbolic" (Δ > 0, e.g. the wave equation), "Parabolic" (Δ == 0, e.g. the heat equation), or "Elliptic" (Δ < 0, e.g. Laplace's equation). Only the highest-order terms determine the type. A discriminant whose sign is not a decidable constant (a mixed-type / parameter-dependent equation such as Tricomi's y u_xx + u_yy == 0) leaves the call unevaluated.`**

## Examples

_No verified examples yet for this function._

## Algorithm

dsolve_pdeclassify.c — PDEClassify[eqn, u, {v1, v2}]: classify a second-order linear PDE by the discriminant of its principal part.

```text
For  A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} + (lower order) == 0  the type is
set by  Δ = B² − 4 A C  (only the highest-order terms matter):

    Δ > 0  →  "Hyperbolic"   (e.g. the wave equation u_tt == c² u_xx)
    Δ = 0  →  "Parabolic"    (e.g. the heat equation u_t == u_xx)
    Δ < 0  →  "Elliptic"     (e.g. Laplace u_xx + u_yy == 0)
```

A, B, C are read from the equation's second-order terms exactly as

```text
DSolve`PDELinearSecondOrder reads them.  The sign of Δ is decided for a
```

definite constant discriminant; a variable-coefficient / parameter-dependent Δ whose sign is not a decidable constant (a mixed-type equation such as the Tricomi equation y u_xx + u_yy == 0) leaves the call unevaluated — an honest

```text
decline rather than a region-blind label.  A first cut: linear in the
```

second-order terms, over two independent variables.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/calculus/dsolve_pdeclassify.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/dsolve_pdeclassify.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
