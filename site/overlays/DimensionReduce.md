### Worked examples

```mathematica
In[1]:= d = {{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}}  (* a 5x3 data matrix *)
```

```mathematica
In[1]:= DimensionReduce[d, 2]  (* default PCA: reduce each row to two dimensions *)
```

```mathematica
In[1]:= DimensionReduce[d, 2] == Map[Take[#, 2] &, PrincipalComponents[d]]  (* reducing to k is exactly the first k principal components *)
```

```mathematica
In[1]:= DimensionReduce[d, 2, Method -> "LatentSemanticAnalysis"]  (* LSA skips the centring: a truncated SVD *)
```

```mathematica
In[1]:= DimensionReduce[d, 2, Method -> "MultidimensionalScaling"]  (* classical MDS, the PCA embedding up to a reflection *)
```

### Notes

The three methods are one algorithm with three ways of forming the symmetric matrix to
decompose: PCA uses the covariance of the centred columns, LSA the Gram matrix `XᵀX`
without centring, and MDS the double-centred squared-distance matrix. Skipping the
centring is the *entire* difference between PCA and LSA — a term-document matrix is
sparse and non-negative, and centring would destroy both.

Classical MDS on Euclidean distances is the same embedding as PCA, reached by a
different route (an `n × n` distance matrix rather than a `dim × dim` covariance); it is
identical up to the sign of each axis, and earns its keep only when the distances come
from somewhere other than the coordinates. MDS is capped at 2000 rows because its matrix
is `n × n`.

Asking for more dimensions than the data supports returns unevaluated rather than
padding with zeros, which would look like a successful reduction to a caller checking
only the shape.
