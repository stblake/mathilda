### Worked examples

```mathematica
In[1]:= EuclideanDistance[{0, 0}, {3, 4}]  (* the 3-4-5 right triangle *)
```

```mathematica
In[1]:= EuclideanDistance[{1, 1, 1}, {4, 5, 1}]  (* straight-line distance in 3-D *)
```

```mathematica
In[1]:= EuclideanDistance[3, 7]  (* scalars are treated as 1-vectors *)
```

### Notes

`EuclideanDistance[u, v]` is `Sqrt[Sum Abs[u_i - v_i]^2]`, the ordinary
straight-line distance. Both arguments must be scalars, or lists of equal length;
a length mismatch or matrix-shaped input is left unevaluated. Complex components
contribute their modulus (the definition squares `Abs`, not the raw difference),
and a symbolic pair returns a symbolic distance rather than an error.

For ranking or clustering, prefer `SquaredEuclideanDistance`: it avoids the root,
stays exact for exact input, and orders points identically.
