### Worked examples

```mathematica
In[1]:= AlgebraicNumberNorm[Sqrt[2]]  (* product of the roots of x^2 - 2 *)
```

```mathematica
In[1]:= AlgebraicNumberNorm[1 + Sqrt[2]]  (* (1 + Sqrt[2])(1 - Sqrt[2]) *)
```

```mathematica
In[1]:= AlgebraicNumberNorm[GoldenRatio]  (* the norm of the golden ratio *)
```

```mathematica
In[1]:= AlgebraicNumberNorm[Sqrt[2] + Sqrt[3]]  (* a degree-4 field: product of four conjugates *)
```

```mathematica
In[1]:= AlgebraicNumberNorm[1 + Sqrt[2], Extension -> Sqrt[2]]  (* the relative norm over Q(Sqrt[2]) *)
```

### Notes

The absolute norm `N_{Q(a)/Q}(a)` is the product of the conjugates of `a` —
equivalently the product of the roots of `a`'s minimal polynomial, read off as
`(-1)^deg` times the monic constant term. `AlgebraicNumberNorm[a, Extension ->
theta]` gives the relative norm `N_{Q(theta)/Q}(a)` for `a` an element of
`Q(theta)`, equal to the absolute norm raised to the tower index `[Q(theta):Q(a)]`.

`Listable` (a trailing `Extension` option is repeated across the list) and
`Protected`; requires FLINT. If `a` is not an element of the stated extension,
the call reports and stays unevaluated.
