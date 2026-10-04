### Worked examples

```mathematica
In[1]:= InexactNumberQ[1.5]  (* a machine real is inexact *)
```

```mathematica
In[1]:= InexactNumberQ[2.0 + 3 I]  (* one inexact part makes the Complex inexact *)
```

```mathematica
In[1]:= InexactNumberQ[1/3]  (* an exact rational is not inexact *)
```

```mathematica
In[1]:= InexactNumberQ[2 + 3 I]  (* an all-exact Complex is not inexact *)
```

### Notes

`InexactNumberQ[expr]` is `True` for machine reals, arbitrary-precision (MPFR)
reals, and `Complex` numbers with at least one inexact part. It is the complement
of `ExactNumberQ` among numbers, so integers, rationals and all-exact `Complex`
numbers are `False`, as is any non-number. Non-single-argument calls are left
unevaluated.
