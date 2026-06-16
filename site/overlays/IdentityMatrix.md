### Worked examples

```mathematica
In[1]:= IdentityMatrix[3]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

A two-element argument gives a rectangular identity (1s on the main diagonal,
0s elsewhere):

```mathematica
In[1]:= IdentityMatrix[{2, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}}
```

It is the multiplicative identity for matrix products — multiplying any matrix
by a conformant identity leaves it unchanged:

```mathematica
In[1]:= IdentityMatrix[4] . HilbertMatrix[4] == HilbertMatrix[4]
Out[1]= True
```

### Notes

`IdentityMatrix[n]` gives the `n x n` identity; `IdentityMatrix[{m, n}]` gives the `m x n` rectangular identity. Use it as a seed for matrix algebra (e.g. `MatrixPower`, characteristic-matrix constructions) and as the neutral element of `Dot`.
