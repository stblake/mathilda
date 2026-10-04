---
source: src/poly/ratcanon.c
---
**Algorithm.** `RatCanonPrototype[expr]` is a **Phase-1 prototype** of rational-function
canonicalisation over the differential/algebraic tower of `expr`. `builtin_ratcanon_prototype`
(`src/poly/ratcanon.c`) performs a kernel-substitution round trip: `rcp_forward` walks the
expression and replaces every non-rational kernel (each `Exp[...]`, `Log[...]`, `Sqrt[...]`,
`Sin[...]`, etc.) with a fresh free symbol, recording the map; `flint_rational_together`
reduces the resulting plain rational function over `Q` in **one** FLINT reduction;
`rcp_backward` substitutes the kernels back; and `eval_and_free` applies the algebraic
relations that the free-symbol abstraction hid (`I^2 -> -1`, `Sqrt[k]^2 -> k`, and so on). If
the reduction declines, the builtin returns `NULL` and the call is left unevaluated.

**Data structures.** An `RcpMap` holds the kernel↔symbol correspondence; everything else is
ordinary `Expr` passed into the FLINT rational engine (`flint_rational_together`, the same
kernel `Together`/`Cancel` use). It is registered `Protected` in `ratcanon_init`.

**Complexity / limits.** Dominated by the single FLINT reduction over `Q`, which is why it
cancels common factors cheaply — `(x^2-1)/(x-1)` collapses to `1 + x`, and the exponential
tower `(E^(2x)-1)/(E^x-1)` to `1 + E^x` by treating `E^x` as one generator. It is a
prototype: it reduces over `Q` with each kernel as an independent generator, so it does
**not** know the algebraic relations *between* kernels (e.g. it declines a mixed
`Sqrt`-in-denominator case and the Pythagorean `Sin`/`Cos` identity). It is a building block
for a fuller tower-canonicalisation, not a general simplifier — use `Together`, `Cancel` or
`Simplify` for production work.
