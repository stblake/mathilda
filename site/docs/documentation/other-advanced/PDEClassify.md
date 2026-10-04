# PDEClassify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PDEClassify[eqn, u, {v1, v2}] classifies a second-order linear PDE by the discriminant Δ = B² − 4 A C of its principal part A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2}: "Hyperbolic" (Δ > 0, e.g. the wave equation), "Parabolic" (Δ == 0, e.g. the heat equation), or "Elliptic" (Δ < 0, e.g. Laplace's equation). Only the highest-order terms determine the type. A discriminant whose sign is not a decidable constant (a mixed-type / parameter-dependent equation such as Tricomi's y u_xx + u_yy == 0) leaves the call unevaluated.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Wave equation: B^2 - 4AC > 0

```mathematica
In[1]:= PDEClassify[D[u[x, t], {t, 2}] == D[u[x, t], {x, 2}], u, {x, t}]
Out[1]= "Hyperbolic"
```

Laplace: discriminant < 0

```mathematica
In[2]:= PDEClassify[D[u[x, y], {x, 2}] + D[u[x, y], {y, 2}] == 0, u, {x, y}]
Out[2]= "Elliptic"
```

Heat: discriminant = 0

```mathematica
In[3]:= PDEClassify[D[u[x, t], t] == D[u[x, t], {x, 2}], u, {x, t}]
Out[3]= "Parabolic"
```

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

**Algorithm.** `builtin_pdeclassify` classifies a second-order linear PDE
`A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2} + (lower order) == 0` by the discriminant of
its principal part, `Δ = B² − 4 A C`:

1. Parse the equation with `dsolve_parse` and require a PDE in exactly one dependent
   function and two independent variables (`u`, `{v1, v2}`); otherwise decline.
2. Substitute each second-order derivative — `Derivative[2,0][u]`, `Derivative[1,1][u]`,
   `Derivative[0,2][u]` — by a fresh atom, then recover `A`, `B`, `C` as the partial
   derivatives of the residual with respect to those atoms. This both reads the
   coefficients and checks linearity in the second-order terms (each `A, B, C` must be
   free of the atoms) and that the equation is genuinely second order.
3. Form `Δ = B² − 4 A C`, `Simplify` it, and decide its sign (`pdec_sign`): an integer
   or real discriminant is read directly, otherwise `Sign[Δ]` is evaluated.

`Δ > 0 → "Hyperbolic"`, `Δ = 0 → "Parabolic"`, `Δ < 0 → "Elliptic"` (returned as
strings). Only the highest-order terms determine the type.

**Data structures.** `Expr` trees throughout, driven by the `ds_*` DSolve substrate
(`ds_subst`, `ds_d`, `ds_simplify`, `ds_free_of`, `ds_call*`). The three substituted
derivative atoms are interned context symbols (`DSolve\`pdecUxx`, etc.).

**Complexity / limits.** Dominated by the three `D` passes and the `Simplify` of `Δ`.
Attributes: `Protected`. It is deliberately conservative: a discriminant whose sign is
not a decidable constant — a mixed-type / parameter-dependent equation such as
Tricomi's `y u_xx + u_yy == 0` — leaves the call unevaluated rather than printing a
region-blind label. Non-linear or higher-order equations, or a wrong number of
functions/variables, also decline.

**Attributes:** `Protected`.

## References

- Source: [`src/calculus/dsolve_pdeclassify.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/dsolve_pdeclassify.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

`PDEClassify[eqn, u, {v1, v2}]` classifies a second-order linear PDE by the discriminant
`Δ = B² − 4 A C` of its principal part `A u_{v1 v1} + B u_{v1 v2} + C u_{v2 v2}`:
`"Hyperbolic"` for `Δ > 0` (the wave equation), `"Parabolic"` for `Δ = 0` (the heat
equation), and `"Elliptic"` for `Δ < 0` (Laplace's equation). Only the highest-order
terms determine the type.

A discriminant whose sign is not a decidable constant — a mixed-type or
parameter-dependent equation such as Tricomi's `y u_xx + u_yy == 0` — leaves the call
unevaluated, an honest decline rather than a region-blind label. Non-linear or
higher-order equations, or anything other than one function in two variables, also
decline. `PDEClassify` is `Protected`.
