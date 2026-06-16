### Worked examples

```mathematica
In[1]:= LeafCount[x + y]
Out[1]= 3
```

It counts every atomic subexpression, including operator heads:

```mathematica
In[1]:= LeafCount[a + b^2 + Sin[c d]]
Out[1]= 9
```

A handy proxy for symbolic "size" — here the blow-up of an expanded binomial:

```mathematica
In[1]:= LeafCount[Expand[(1 + x)^10]]
Out[1]= 48
```

Measuring the complexity of a computed result, e.g. an antiderivative:

```mathematica
In[1]:= LeafCount[Integrate[1/(1 + x^4), x]]
Out[1]= 89
```

Mapped over a list, it ranks expressions by structural weight:

```mathematica
In[1]:= Map[LeafCount, {1, 1/2, x, f[x], {a, b, c}}]
Out[1]= {1, 3, 1, 2, 4}
```

### Notes

`LeafCount[expr]` gives the total number of indivisible subexpressions (leaves)
in `expr`, counting heads and structural atoms. It is the standard measure used
by `Simplify` and friends to decide which of two candidate forms is "smaller".
