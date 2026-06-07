### Worked examples

```mathematica
In[1]:= MatrixQ[{{1, 2}, {3, 4}}]
Out[1]= True

In[2]:= MatrixQ[{1, 2, 3}]
Out[2]= False
```

### Notes

`MatrixQ` is `True` for a list of equal-length lists; a flat (rank-1) list is a vector, not a matrix.
