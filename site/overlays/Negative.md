### Worked examples

```mathematica
In[1]:= Negative[1 - Pi]
```

```mathematica
In[1]:= Negative[-5]
```

```mathematica
In[1]:= Negative[1 + I]  (* a non-real complex value is not negative *)
```

```mathematica
In[1]:= Negative[{1.6, 3/4, Pi, 0, -5}]  (* Listable: threads over the list *)
```

```mathematica
In[1]:= Negative[x]  (* non-numeric argument stays symbolic *)
```

### Notes

`Negative[x]` is `True` only for a real, strictly negative numeric quantity. Zero
gives `False` (use `NonPositive` to include it), and so does any non-real complex
value. Exact integers, rationals, and bigints are decided exactly; everything else
numeric is classified by its machine-precision value.

A non-numeric argument is left unevaluated. `Negative` is `Listable` and reads a
packed list or `NDArray` straight off the buffer, returning an ordinary list of
`True`/`False`.
