### Worked examples

A `SeriesData` object prints as a sum of coefficients times powers of `x - x0`:

```mathematica
In[1]:= SeriesData[x, 0, {1, 1, 1}, 0, 3, 1]
Out[1]= 1 + x + x^2 + O[x]^3
```

`Normal` drops the order term, recovering an ordinary expression:

```mathematica
In[1]:= Normal[SeriesData[x, 0, {1, 1, 1}, 0, 3, 1]]
Out[1]= 1 + x + x^2
```

`InputForm` reveals the raw representation `Series` builds — here the Taylor data
for `Sin[x]`, with coefficient list, `nmin = 0`, `nmax = 5`, denominator `1`:

```mathematica
In[1]:= InputForm[Series[Sin[x], {x, 0, 4}]]
Out[1]= SeriesData[x, 0, {0, 1, 0, -1/6, 0}, 0, 5, 1]
```

A Laurent series has a negative `nmin`. The expansion of `1/(E^x - 1)` starts at
`x^(-1)`, so `nmin = -1`:

```mathematica
In[1]:= InputForm[Series[1/(Exp[x] - 1), {x, 0, 3}]]
Out[1]= SeriesData[x, 0, {1, -1/2, 1/12, 0, -1/720}, -1, 4, 1]
```

A Puiseux series uses a denominator greater than `1`. For `Sqrt[x] + x` the powers
are half-integers, so `den = 2`:

```mathematica
In[1]:= InputForm[Series[Sqrt[x] + x, {x, 0, 2}]]
Out[1]= SeriesData[x, 0, {0, 1, 1, 0, 0}, 0, 5, 2]
```

### Notes

`SeriesData[x, x0, {a0, a1, ...}, nmin, nmax, den]` is the internal representation
of a power series in `x` about `x0`. The powers that appear are
`nmin/den, (nmin+1)/den, ..., (nmax-1)/den`, with a trailing `O[x - x0]^(nmax/den)`
term standing in for the omitted tail. The single uniform structure covers Taylor
(`nmin >= 0`, `den = 1`), Laurent (`nmin < 0`), and Puiseux (`den > 1`) series.
These objects are produced by `Series`; use `Normal` to convert one back to an
ordinary polynomial by discarding the O-term, and `InputForm` to see the literal
`SeriesData[...]` form instead of the pretty-printed sum.
