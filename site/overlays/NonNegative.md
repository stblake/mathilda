### Worked examples

```mathematica
In[1]:= NonNegative[0]  (* zero is included *)
```

```mathematica
In[1]:= NonNegative[Pi - 3]
```

```mathematica
In[1]:= NonNegative[-3]
```

```mathematica
In[1]:= NonNegative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]  (* Listable *)
```

```mathematica
In[1]:= NonNegative[x]  (* non-numeric argument stays symbolic *)
```

### Notes

`NonNegative[x]` is `True` for a real numeric quantity that is positive **or
zero**; it differs from `Positive` only at zero. A non-real complex value gives
`False`. Exact inputs are decided exactly, inexact and symbolic-constant inputs by
their machine-precision value, and a non-numeric argument is left unevaluated.
`NonNegative` is `Listable` and reads a packed list or `NDArray` straight off the
buffer.
