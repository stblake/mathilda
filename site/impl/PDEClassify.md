---
source: src/calculus/dsolve_pdeclassify.c
---
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
