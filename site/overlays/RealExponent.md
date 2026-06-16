### Worked examples

```mathematica
In[1]:= RealExponent[1234.5]
Out[1]= 3.09149
```

In an explicit base the result is the exact logarithm — the base-2 exponent of a
pure power of two is an integer, which makes `RealExponent` a quick way to count
the decimal digits of a huge integer (`Floor[log10] + 1`):

```mathematica
In[1]:= RealExponent[2^100, 2]
Out[1]= 100.0

In[2]:= Floor[RealExponent[2^1000]] + 1
Out[2]= 302
```

It accepts symbolic numeric values and lifts to MPFR precision when given an MPFR
argument — for example `Log10[E]` to 40 digits:

```mathematica
In[1]:= RealExponent[N[Pi^Pi]]
Out[1]= 1.56184

In[2]:= RealExponent[N[E, 40]]
Out[2]= 0.4342944819032518276511289189166050822944
```

`RealExponent` threads over lists:

```mathematica
In[1]:= RealExponent[{10, 100, 1000}]
Out[1]= {1.0, 2.0, 3.0}
```

### Notes

`RealExponent[x]` gives `Log[10, |x|]`; `RealExponent[x, b]` gives `Log[b, |x|]`.
It accepts `Integer`, `BigInt`, `Rational`, `Real`, MPFR, and symbolic numeric
values such as `Pi`, `E`, or `Pi^Pi`. The result is a machine `Real` unless an
MPFR input lifts it to that precision. Exact zero gives `-Infinity`; the sign of
`x` is ignored.
