# RatCanonPrototype

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RatCanonPrototype[expr] (Phase-1 prototype) reduces a rational function over the differential/algebraic tower of expr via one FLINT reduction.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Cancels the common factor over Q

```mathematica
In[1]:= RatCanonPrototype[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

E^x is treated as one tower generator

```mathematica
In[2]:= RatCanonPrototype[(E^(2 x) - 1)/(E^x - 1)]
Out[2]= 1 + E^x
```

A Log tower kernel

```mathematica
In[3]:= RatCanonPrototype[(Log[x]^2 - 1)/(Log[x] + 1)]
Out[3]= -1 + Log[x]
```

Combines to a single reduced fraction

```mathematica
In[4]:= RatCanonPrototype[x/(x + 1) + 1/(x + 1)]
Out[4]= 1
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/poly/ratcanon.c`](https://github.com/stblake/mathilda/blob/main/src/poly/ratcanon.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_ratcanon_spec.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ratcanon_spec.c)

## Notes & additional examples

### Notes

`RatCanonPrototype[expr]` is a Phase-1 prototype that reduces a rational function over the
differential/algebraic tower of `expr` via a single FLINT reduction. It abstracts each
non-rational kernel (`E^x`, `Log[x]`, …) to a fresh generator, reduces over `Q`, substitutes
the kernels back, and re-applies the algebraic relations — so it cancels common factors and
combines fractions cheaply, as in the examples above.

It is a prototype, not a general simplifier: it treats each kernel as an independent
generator, so it does not know relations *between* kernels and will leave an expression
unevaluated when its heuristic declines (for instance a `Sqrt` in the denominator, or the
Pythagorean `Sin`/`Cos` identity). For production use prefer `Together`, `Cancel` or
`Simplify`.
