### Worked examples

```mathematica
In[1]:= AlgebraicIntegerQ[Sqrt[2]]  (* a root of the monic x^2 - 2 *)
```

```mathematica
In[1]:= AlgebraicIntegerQ[1/2]  (* a non-integer rational is not an algebraic integer *)
```

```mathematica
In[1]:= AlgebraicIntegerQ[GoldenRatio]  (* a root of the monic x^2 - x - 1 *)
```

```mathematica
In[1]:= AlgebraicIntegerQ[(1 + Sqrt[5])/2]  (* the same number, spelled as a radical *)
```

```mathematica
In[1]:= AlgebraicIntegerQ[2 + 3 I]  (* Gaussian integers are algebraic integers *)
```

```mathematica
In[1]:= AlgebraicIntegerQ[Sqrt[2]/3]  (* scaling by 1/3 leaves the monogenic ring *)
```

### Notes

The test is exact: `x` is an algebraic integer iff the leading coefficient of its
primitive integer minimal polynomial is `1`. A rational `p/q` has minimal
polynomial `q x - p`, so only ordinary integers qualify among the rationals.

Anything that is not a constant algebraic number — a free symbol, `Pi`,
`Log[2]` — is simply not an algebraic integer, and returns `False`. The predicate
is **not** `Listable`: `AlgebraicIntegerQ[list]` asks about the list itself.
Requires FLINT; with FLINT compiled out the call stays unevaluated.
