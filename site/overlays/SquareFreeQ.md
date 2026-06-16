### Worked examples

```mathematica
In[1]:= SquareFreeQ[12]
Out[1]= False
```

A square-free integer has no repeated prime factor:

```mathematica
In[1]:= SquareFreeQ[30]
Out[1]= True
```

For polynomials it detects repeated factors:

```mathematica
In[1]:= SquareFreeQ[x^2 - 1]
Out[1]= True

In[2]:= SquareFreeQ[(x - 1)^2 (x + 1)]
Out[2]= False
```

Square-freeness depends on the coefficient ring: `2 = -i (1 + i)^2` is *not*
square-free over the Gaussian integers, even though it is over the rationals:

```mathematica
In[1]:= SquareFreeQ[2]
Out[1]= True

In[2]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[2]= False
```

The cyclotomic-style polynomial `x^4 + x^2 + 1` has distinct irreducible
factors and is square-free:

```mathematica
In[1]:= SquareFreeQ[x^4 + x^2 + 1]
Out[1]= True
```

### Notes

`SquareFreeQ[expr]` tests an integer or a polynomial for the absence of any
repeated factor. Over the integers it checks the prime factorization; over a
polynomial ring it is decided from `GCD[p, p']` (the polynomial is square-free
exactly when this GCD is constant). The `GaussianIntegers -> True` option moves
the test into `Z[i]`, where rational primes such as `2` can acquire a repeated
factor. A second argument restricts the test to the given variables.
