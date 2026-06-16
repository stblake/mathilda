### Worked examples

```mathematica
In[1]:= Root[#^2 - 2 &, 1]
Out[1]= Root[#1^2 - 2 &, 1]

In[2]:= N[Root[#^2 - 2 &, 2], 40]
Out[2]= -1.4142135623730950488016887242096980785697
```

```mathematica
In[1]:= N[Root[#^5 - # - 1 &, 1], 40]
Out[1]= 1.1673039782614186842560458998548421807206
```

```mathematica
In[1]:= {N[Root[#^3 - 2 &, 1], 30], N[Root[#^3 - 2 &, 2], 30], N[Root[#^3 - 2 &, 3], 30]}
Out[1]= {1.259921049894873164767210607278, -0.6299605249474365823836053036392 - 1.09112363597172140356007261419*I, -0.6299605249474365823836053036392 + 1.09112363597172140356007261419*I}
```

```mathematica
In[1]:= {N[Root[#^4 + # + 1 &, 1], 20], N[Root[#^4 + # + 1 &, 2], 20], N[Root[#^4 + # + 1 &, 3], 20], N[Root[#^4 + # + 1 &, 4], 20]}
Out[1]= {-0.72713608449119683998 - 0.430014288329715776416*I, -0.72713608449119683998 + 0.430014288329715776416*I, 0.72713608449119683998 - 0.934099289460529439642*I, 0.72713608449119683998 + 0.934099289460529439642*I}
```

### Notes

`Root[f, k]` is the exact, held representation of the `k`-th root of a
univariate polynomial — including roots like that of `#^5 - # - 1 &`, the
classic example of a quintic with **no radical solution** (Abel–Ruffini), which
`Root` nonetheless names and evaluates to arbitrary precision. The index `k` is
canonical: real roots come first in ascending order, then complex roots ordered
by real part, magnitude of imaginary part, and finally the negative-imaginary
member of each conjugate pair. `N[Root[..], prec]` drives a companion-matrix +
Sturm + Newton pipeline to the requested precision; the `#^3 - 2 &` and
`#^4 + # + 1 &` examples recover the full set of real and complex roots in the
canonical order.
