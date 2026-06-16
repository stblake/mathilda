### Worked examples

```mathematica
In[1]:= Total[{1, 2, 3, 4}]
Out[1]= 10
```

On a matrix, the one-argument form sums the rows (a column total):

```mathematica
In[1]:= Total[{{1, 2}, {3, 4}, {5, 6}}]
Out[1]= {9, 12}
```

A level specification controls the depth of summation: `{2}` sums each column
instead, giving the per-column totals.

```mathematica
In[1]:= Total[{{1, 2}, {3, 4}, {5, 6}}, {2}]
Out[1]= {3, 7}
```

The tenth row of Pascal's triangle sums to a power of two:

```mathematica
In[1]:= Total[Table[Binomial[10, k], {k, 0, 10}]]
Out[1]= 1024
```

Summing exact rationals stays exact — a partial sum of the Basel series as a
single reduced fraction:

```mathematica
In[1]:= Total[Table[1/k^2, {k, 1, 100}]]
Out[1]= 1589508694133037873112297928517553859702383498543709859889432834803818131090369901/972186144434381030589657976672623144161975583995746241782720354705517986165248000
```
