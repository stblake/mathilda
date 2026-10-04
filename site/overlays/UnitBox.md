### Worked examples

```mathematica
In[1]:= UnitBox[0]  (* inside the unit box, the value is 1 *)
```

```mathematica
In[1]:= UnitBox[1/2]  (* the box is closed at both endpoints *)
```

```mathematica
In[1]:= UnitBox[0.7]  (* outside the box, the value is 0 *)
```

```mathematica
In[1]:= UnitBox[{-1, -1/2, 0, 1/2, 0.7}]  (* Listable over a vector of test points *)
```

```mathematica
In[1]:= UnitBox[x]  (* a symbolic argument is left unevaluated *)
```

### Notes

`UnitBox` is the rectangular pulse: 1 on the closed interval `-1/2 <= x <= 1/2`
and 0 outside it. Both endpoints belong to the box (`UnitBox[1/2] = 1`), matching
the closed-at-zero convention of `UnitStep`, and the result is always the exact
integer 0 or 1 once the argument's position is certified.

Internally the two-sided test reuses `UnitStep`'s one-sided sign certification on
the shifted arguments `x + 1/2` and `1/2 - x`. Because a box is not monotone,
`UnitBox` is the one member of this family that does **not** thread over an
`Interval`; it does carry a narrowing `NDArray` kernel and a `Compile[]` lowering
(`UnitBox[0.5]` compiles to the integer `1`).
