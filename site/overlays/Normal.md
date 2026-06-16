### Worked examples

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5
```

Drop the O-term from the Maclaurin series of `Sin[x]/x` to recover the truncated
polynomial:

```mathematica
In[1]:= Normal[Series[Sin[x]/x, {x, 0, 6}]]
Out[1]= 1 - 1/6 x^2 + 1/120 x^4 - 1/5040 x^6
```

The tangent series, with its Bernoulli-number coefficients laid bare:

```mathematica
In[1]:= Normal[Series[Tan[x], {x, 0, 7}]]
Out[1]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7
```

The alternating-harmonic expansion of `Log[1 + x]`:

```mathematica
In[1]:= Normal[Series[Log[1 + x], {x, 0, 5}]]
Out[1]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + 1/5 x^5
```

### Notes

`Normal[expr]` converts `expr` to a normal expression. Applied to a `SeriesData`
object it drops the `O`-term and returns the truncated polynomial (or
Laurent/Puiseux sum) as an ordinary `Plus` expression, ready to be added,
differentiated, or substituted into. Expressions that are already normal pass
through unchanged.
