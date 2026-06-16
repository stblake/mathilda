### Worked examples

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {2, 1}}]
Out[1]= True

In[2]:= SymmetricMatrixQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

A complex symmetric matrix is symmetric without being Hermitian:

```mathematica
In[1]:= SymmetricMatrixQ[{{1, I}, {I, 1}}]
Out[1]= True

In[2]:= HermitianMatrixQ[{{1, I}, {I, 1}}]
Out[2]= False
```

`Tolerance` accepts numerically near-symmetric matrices:

```mathematica
In[1]:= SymmetricMatrixQ[{{1.0, 2.0001}, {2.0, 1.0}}, Tolerance -> 0.001]
Out[1]= True
```

A custom `SameTest` relaxes equality of off-diagonal entries:

```mathematica
In[1]:= SymmetricMatrixQ[{{1, 2}, {3, 4}}, SameTest -> (Abs[#1 - #2] <= 1 &)]
Out[1]= True
```

### Notes

`SymmetricMatrixQ` uses the definition `m^T == m` for both real and complex matrices, so a complex symmetric matrix need not be Hermitian. Use the `SameTest` or `Tolerance` options for approximate or custom equality.
