### Worked examples

```mathematica
In[1]:= Rational[6, 4]
Out[1]= 3/2
```

`Rational` auto-reduces by the gcd, normalises the sign onto the numerator, and
collapses to an `Integer` whenever the denominator divides the numerator:

```mathematica
In[1]:= Rational[10, 2]
Out[1]= 5

In[2]:= Rational[-3, -9]
Out[2]= 1/3
```

Because rationals propagate exactly through `Plus` and `Times` via GMP, exact
sums never drift into floating point — a partial sum of the Basel series stays a
single reduced fraction:

```mathematica
In[1]:= 1/2 + 1/3 + 1/6
Out[1]= 1

In[2]:= Sum[1/k^2, {k, 1, 10}]
Out[2]= 1968329/1270080
```

### Notes

`Rational[n, d]` represents the rational number `n/d`. With integer arguments it
reduces to lowest terms, moves the sign to the numerator, and becomes an
`Integer` when `d` divides `n`. The head of any non-integer fraction such as
`1/2` is `Rational`.
