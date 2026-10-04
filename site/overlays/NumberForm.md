### Worked examples

```mathematica
In[1]:= NumberForm[N[Pi], 10]  (* ten significant digits *)
```

```mathematica
In[1]:= NumberForm[10^9, DigitBlock -> 3]  (* thousands separators *)
```

```mathematica
In[1]:= NumberForm[1234.567, {5, 2}]  (* five significant figures, two decimals *)
```

```mathematica
In[1]:= NumberForm[{8.^5, 11.^7, 13.^9}, NumberFormat -> (Row[{#1, "e", #3}] &)]  (* custom layout *)
```

```mathematica
In[1]:= FullForm[NumberForm[1.23, 2]]  (* the wrapper head survives in the tree *)
```

### Notes

`NumberForm[expr, n]` displays the approximate reals in `expr` to `n` significant
digits; `NumberForm[expr, {n, f}]` uses `n` significant figures shown with exactly
`f` digits after the point. It works over integers, scalars, lists, matrices, and
mixed symbolic expressions — every inexact real inside `expr` is reformatted — and
takes a rich option set (`DigitBlock`, `NumberSeparator`, `NumberFormat`,
`ScientificNotationThreshold`, and more).

It is an inert **print wrapper**: the `NumberForm[...]` head stays in the
expression tree (so `FullForm` shows it) and only changes how the wrapped value is
displayed. Because the head survives, an intervening `NumberForm` blocks arithmetic
on the surrounding expression, so assign a variable first if the result must stay
computable. A requested precision below the integer-digit count issues
`NumberForm::reqsigz` and pads with zeros.
