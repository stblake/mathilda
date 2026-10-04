### Worked examples

```mathematica
In[1]:= ExactNumberQ[1/3]  (* an exact rational *)
```

```mathematica
In[1]:= ExactNumberQ[2 + 3 I]  (* a Complex with exact parts is exact *)
```

```mathematica
In[1]:= ExactNumberQ[2.0 + 3 I]  (* one inexact part makes the whole inexact *)
```

```mathematica
In[1]:= ExactNumberQ[1.5]  (* a machine real is inexact *)
```

```mathematica
In[1]:= ExactNumberQ[x]  (* a symbol is not a number at all *)
```

### Notes

`ExactNumberQ[expr]` is `True` for integers, rationals, and `Complex` numbers whose
real and imaginary parts are *both* exact. Reals and arbitrary-precision (MPFR)
numbers are inexact, so they are `False` — and a single inexact part makes a whole
`Complex` inexact. It is the complement of `InexactNumberQ` among numbers; a
non-number gives `False`.
