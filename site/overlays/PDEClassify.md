### Worked examples

```mathematica
In[1]:= PDEClassify[D[u[x, t], {t, 2}] == D[u[x, t], {x, 2}], u, {x, t}]  (* wave equation: B^2 - 4AC > 0 *)
```

```mathematica
In[1]:= PDEClassify[D[u[x, y], {x, 2}] + D[u[x, y], {y, 2}] == 0, u, {x, y}]  (* Laplace: discriminant < 0 *)
```

```mathematica
In[1]:= PDEClassify[D[u[x, t], t] == D[u[x, t], {x, 2}], u, {x, t}]  (* heat: discriminant = 0 *)
```

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
