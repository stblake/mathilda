### Worked examples

```mathematica
In[1]:= NumberFieldIntegralBasis[Sqrt[2]]  (* Z[Sqrt[2]] is already the maximal order *)
```

```mathematica
In[1]:= NumberFieldIntegralBasis[Sqrt[5]]  (* disc 5: the ring of integers needs (1 + Sqrt[5])/2 *)
```

### Notes

`NumberFieldIntegralBasis[a]` gives a `Z`-module basis of the ring of integers
`O_K` of `K = Q(a)`: a list of algebraic integers generating `O_K` over the
integers. The equation order `Z[a]` is enlarged to `O_K` by Round 2
(Pohst–Zassenhaus) when it is not already maximal, so the result is correct even
for non-monogenic fields. The `Sqrt[5]` example shows this: the basis is not
`{1, Sqrt[5]}` but `{1, (1 + Sqrt[5])/2}` — the half-integer generator — because
the discriminant of `x^2 - 5` is `20`, four times the field discriminant `5`.

The basis is presented in Hermite normal form (so the `k`-th element has degree
`k`), but any `Z`-basis of `O_K` is admissible — the choice is not unique.
`Listable` (each generator is its own field) and `Protected`; requires FLINT.
