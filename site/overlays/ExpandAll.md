### Worked examples

```mathematica
In[1]:= ExpandAll[Sin[(x + 1)^2]]  (* reaches inside a function argument, unlike Expand *)
```

```mathematica
In[2]:= ExpandAll[Exp[(a + b)^2]]  (* the exponent is expanded too *)
```

```mathematica
In[3]:= ExpandAll[(x + 1)^2/(x + 2)^2]  (* the denominator is expanded as well *)
```

```mathematica
In[4]:= ExpandAll[(x + 1)^2 (y + 1)^2, x]  (* a second argument confines expansion to parts containing x *)
```

### Notes

`ExpandAll[expr]` expands products and integer powers in *every* part of `expr` —
function heads and arguments, exponents, and the bases of denominators — where a
plain `Expand` distributes only at the top level. So `ExpandAll[Sin[(x + 1)^2]]`
expands the argument to `Sin[1 + 2 x + x^2]`, which `Expand` leaves untouched, and
`ExpandAll[(x + 1)^2/(x + 2)^2]` expands the `(x + 2)^2` denominator. The
two-argument form `ExpandAll[expr, patt]` leaves any part free of `patt` alone.
Like `Expand`, an expansion too large to fit in memory returns `Overflow[]` rather
than declining; it threads over lists, equations, inequalities, and logic
functions.
