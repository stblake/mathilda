### Worked examples

```mathematica
In[1]:= ConjugateTranspose[{{1 + I, 2}, {3, 4 - I}}]
Out[1]= {{1 - I, 3}, {2, 4 + I}}
```

```mathematica
In[1]:= ConjugateTranspose[{{a, b}, {c, d}}]
Out[1]= {{Conjugate[a], Conjugate[c]}, {Conjugate[b], Conjugate[d]}}
```

```mathematica
In[1]:= m = {{1, I}, {-I, 2}}; ConjugateTranspose[m] == m
Out[1]= True
```

### Notes

`ConjugateTranspose[m]` is the Hermitian adjoint `Conjugate[Transpose[m]]`. The last example confirms a matrix is Hermitian (equal to its own conjugate transpose). On a vector it conjugates the entries in place.
