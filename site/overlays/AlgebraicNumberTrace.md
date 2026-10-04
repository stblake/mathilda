### Worked examples

```mathematica
In[1]:= AlgebraicNumberTrace[Sqrt[2]]  (* Sqrt[2] + (-Sqrt[2]) = 0 *)
```

```mathematica
In[1]:= AlgebraicNumberTrace[1 + Sqrt[2]]  (* (1 + Sqrt[2]) + (1 - Sqrt[2]) *)
```

```mathematica
In[1]:= AlgebraicNumberTrace[GoldenRatio]  (* the trace of the golden ratio *)
```

```mathematica
In[1]:= AlgebraicNumberTrace[AlgebraicNumber[Sqrt[2], {1, 3}]]  (* trace of 1 + 3 Sqrt[2] *)
```

### Notes

The absolute trace `Tr_{Q(a)/Q}(a)` is the sum of the conjugates of `a` —
equivalently the sum of the roots of `a`'s minimal polynomial, read off as
`-(coeff of x^{deg-1}) / (leading coeff)`. `AlgebraicNumberTrace[a, Extension ->
theta]` gives the relative trace `Tr_{Q(theta)/Q}(a)` for `a` in `Q(theta)`.

Where the norm is multiplicative and *raises* the absolute value to the tower
index, the trace is additive and *scales* by it. `Listable` and `Protected`;
requires FLINT.
