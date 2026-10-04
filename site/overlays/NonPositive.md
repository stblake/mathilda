### Worked examples

```mathematica
In[1]:= NonPositive[0]  (* zero is included *)
```

```mathematica
In[1]:= NonPositive[1 - Pi]
```

```mathematica
In[1]:= NonPositive[3]
```

```mathematica
In[1]:= NonPositive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]  (* Listable *)
```

```mathematica
In[1]:= NonPositive[x]  (* non-numeric argument stays symbolic *)
```

### Notes

`NonPositive[x]` is `True` for a real numeric quantity that is negative **or
zero**; it differs from `Negative` only at zero. A non-real complex value gives
`False`. Exact inputs are decided exactly, inexact and symbolic-constant inputs by
their machine-precision value, and a non-numeric argument is left unevaluated.
`NonPositive` is `Listable` and reads a packed list or `NDArray` straight off the
buffer.
