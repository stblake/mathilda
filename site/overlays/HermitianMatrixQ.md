### Worked examples

```mathematica
In[1]:= HermitianMatrixQ[{{1, I}, {-I, 1}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

A diagonal entry that is not real, or an off-diagonal pair that is not a
conjugate pair, breaks Hermiticity:

```mathematica
In[1]:= HermitianMatrixQ[{{1, 2 + I}, {2 + I, 1}}]
Out[1]= False

In[2]:= HermitianMatrixQ[{{2, 3 + I}, {3 - I, 5}}]
Out[2]= True
```

The predicate also handles inexact matrices, and a `Tolerance` option absorbs
floating-point noise on the diagonal:

```mathematica
In[1]:= HermitianMatrixQ[N[{{1, I}, {-I, 1}}]]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, I}, {-I, 2.0000001}}, Tolerance -> 0.001]
Out[2]= True
```

### Notes

A matrix is Hermitian when `m == ConjugateTranspose[m]`; off-diagonal entries must be conjugates of their transpose partners and diagonal entries must be real. For real matrices this coincides with `SymmetricMatrixQ`.
