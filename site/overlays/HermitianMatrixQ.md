### Worked examples

```mathematica
In[1]:= HermitianMatrixQ[{{1, I}, {-I, 1}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

### Notes

A matrix is Hermitian when `m == ConjugateTranspose[m]`; off-diagonal entries must be conjugates of their transpose partners and diagonal entries must be real. For real matrices this coincides with `SymmetricMatrixQ`.
