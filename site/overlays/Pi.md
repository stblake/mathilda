### Worked examples

The symbol `Pi` stays exact until you ask for a numeric value, and `N[Pi, prec]`
delivers it to any requested precision:

```mathematica
In[1]:= N[Pi]
Out[1]= 3.14159

In[2]:= N[Pi, 40]
Out[2]= 3.1415926535897932384626433832795028841971
```

Exact special values of the elementary functions are reduced symbolically in
terms of `Pi`:

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2

In[2]:= Cos[Pi/3] + Sin[Pi/4]^2
Out[2]= 1
```

`Pi` arises naturally as the closed form of inverse trigonometric values — and
`Cos[Pi/5]` evaluates to the golden-ratio surd, a non-obvious exact constructible
number:

```mathematica
In[1]:= ArcTan[1]
Out[1]= 1/4 Pi

In[2]:= Cos[Pi/5]
Out[2]= 1/4 (1 + Sqrt[5])
```

As a recognised constant, `Pi` participates in exact closed forms that you can
then numericalise to high precision — e.g. the Basel value `Pi^2/6` to 40
digits:

```mathematica
In[1]:= N[Pi^2/6, 40]
Out[1]= 1.6449340668482264364724151666460251892188
```

It is treated as constant by calculus, so its derivative is zero:

```mathematica
In[1]:= D[Pi, x]
Out[1]= 0
```

### Notes

`Pi` is the mathematical constant π. It carries the `Constant` and `Protected`
attributes, `NumericQ[Pi]` is `True`, and `D[Pi, x]` is `0`. It remains an exact
symbol through symbolic computation — driving the special-value reductions of the
trigonometric and inverse-trigonometric functions — and is evaluated to arbitrary
precision only on demand via `N[Pi, prec]`, which uses the MPFR numeric core.
