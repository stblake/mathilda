### Worked examples

```mathematica
In[1]:= RatCanonPrototype[(x^2 - 1)/(x - 1)]  (* cancels the common factor over Q *)
```

```mathematica
In[1]:= RatCanonPrototype[(E^(2 x) - 1)/(E^x - 1)]  (* E^x is treated as one tower generator *)
```

```mathematica
In[1]:= RatCanonPrototype[(Log[x]^2 - 1)/(Log[x] + 1)]  (* a Log tower kernel *)
```

```mathematica
In[1]:= RatCanonPrototype[x/(x + 1) + 1/(x + 1)]  (* combines to a single reduced fraction *)
```

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
