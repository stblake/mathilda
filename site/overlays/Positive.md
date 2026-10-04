### Worked examples

```mathematica
In[1]:= Positive[3]
```

```mathematica
In[1]:= Positive[Pi - 3]  (* a symbolic constant is decided by its numeric value *)
```

```mathematica
In[1]:= Positive[Sqrt[-2]]  (* a non-real complex value is not positive *)
```

```mathematica
In[1]:= Positive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]  (* Listable: threads over the list *)
```

```mathematica
In[1]:= Positive[x]  (* non-numeric argument stays symbolic *)
```

### Notes

`Positive[x]` is `True` only for a real, strictly positive numeric quantity.
Zero gives `False`, and so does any non-real complex value — `Positive` is a test
of a real sign, not of "has a positive real part". Exact integers, rationals, and
bigints are decided exactly; reals, symbolic constants, and numeric-function calls
are classified by their machine-precision numeric value.

A non-numeric argument (one for which `NumericQ` is `False`) is left unevaluated,
so symbolic expressions flow through the evaluator unchanged. `Positive` is
`Listable`, and a packed list or `NDArray` is read straight off the buffer; the
result is a list of `True`/`False`, which no buffer holds, so it comes back as an
ordinary list.
