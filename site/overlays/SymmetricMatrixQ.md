### Worked examples

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {2, 1}}]
Out[1]= True

In[2]:= SymmetricMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

### Notes

`SymmetricMatrixQ` uses the definition `m^T == m` for both real and complex matrices, so a complex symmetric matrix need not be Hermitian. Use the `SameTest` or `Tolerance` options for approximate or custom equality.
