### Worked examples

```mathematica
In[1]:= Tr[{{1, 2}, {3, 4}}]
Out[1]= 5
```

The trace of a symbolic matrix is the sum of its diagonal:

```mathematica
In[1]:= Tr[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[1]= a + e + i
```

A second argument replaces `Plus` with another combiner — here multiplying the
diagonal entries:

```mathematica
In[1]:= Tr[{{1, 2}, {3, 4}}, Times]
Out[1]= 4
```

The trace of the `n`-th power of the Fibonacci `Q`-matrix is the Lucas number
`L[n]`; for `n = 10` this gives `L[10] = 123`:

```mathematica
In[1]:= Tr[MatrixPower[{{1, 1}, {1, 0}}, 10]]
Out[1]= 123
```

Because the trace is basis-independent, `Tr[A . A]` of a symbolic `2x2` matrix
yields the invariant combination of its entries:

```mathematica
In[1]:= Tr[{{a, b}, {c, d}} . {{a, b}, {c, d}}]
Out[1]= a^2 + 2 b c + d^2
```
