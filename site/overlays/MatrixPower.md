---
status: Stable
references:
  - "R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — powers of matrices."
  - "G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — repeated multiplication and inversion."
---
### Worked examples

```mathematica
In[1]:= MatrixPower[{{1, 1}, {0, 1}}, 3]
Out[1]= {{1, 3}, {0, 1}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 3}}, 2]
Out[1]= {{4, 0}, {0, 9}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 4}}, -1]
Out[1]= {{1/2, 0}, {0, 1/4}}
```

```mathematica
In[1]:= MatrixPower[{{2, 0}, {0, 3}}, 2, {1, 1}]
Out[1]= {4, 9}
```

### Notes

`MatrixPower[m, n]` computes the `n`-th power by repeated matrix multiplication. A negative exponent (third example) raises the inverse to the corresponding positive power, so it requires a non-singular matrix; `MatrixPower[m, 0]` returns the identity matrix of the right size. The three-argument form `MatrixPower[m, n, v]` applies the matrix power to a vector — equivalent to `MatrixPower[m, n] . v` but without forming the full power matrix. Fractional exponents are not currently supported.
