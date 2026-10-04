### Worked examples

```mathematica
In[1]:= SeriesCoefficient[Exp[x], {x, 0, 10}]  (* the coefficient of x^10 in e^x is 1/10! *)
```

```mathematica
In[1]:= SeriesCoefficient[Tan[x], {x, 0, 7}]  (* odd coefficients of Tan give the tangent numbers *)
```

```mathematica
In[1]:= SeriesCoefficient[1/(1 - x - x^2), {x, 0, 10}]  (* the Fibonacci generating function: F(10) = 89 *)
```

```mathematica
In[1]:= SeriesCoefficient[Cos[x], {x, 0, 8}]  (* even coefficients of Cos are reciprocals of factorials *)
```

```mathematica
In[1]:= SeriesCoefficient[ProductLog[x], {x, 0, n}]  (* a symbolic index returns the closed-form general term *)
```

### Notes

`SeriesCoefficient[f, {x, x0, k}]` returns the coefficient of `(x - x0)^k` in the
power-series expansion of `f` about `x = x0`, for any `f` that `Series` can expand.
At a concrete integer index it expands `f` to order `k` and extracts the single
coefficient, so the result is exact (rational or symbolic). For a handful of heads
(`ProductLog`, `FresnelC`, `FresnelS`) a **symbolic** index returns the closed-form
general term as a `Piecewise`. The arguments are evaluated normally
(`SeriesCoefficient` is `Protected`, not `HoldAll`); a non-integer or otherwise
unusable index is left unevaluated.
