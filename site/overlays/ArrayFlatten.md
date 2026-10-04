### Worked examples

```mathematica
In[1]:= m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]  (* a 0 scalar fills a zero block *)
```

```mathematica
In[2]:= ArrayFlatten[{{IdentityMatrix[2], {{5}, {6}}}}]  (* two blocks sharing a grid row glue side by side *)
```

```mathematica
In[3]:= ArrayFlatten[{{m, m}, {m, m}}]  (* a 2x2 grid of 2x2 blocks becomes a 4x4 matrix *)
```

### Notes

`ArrayFlatten[a]` treats `a` as a rank-2 grid of matrix blocks and glues them
into one matrix, as `MatrixForm[a]` would show them — equivalent to
`Flatten[a, {{1, 3}, {2, 4}}]`. Blocks must fit: matrices in the same grid row
share their first dimension and those in the same column share their second, and
the output size along an axis is the sum of the block sizes; disagreeing blocks
leave the call unevaluated. A block shallower than a matrix (an atom such as `0`)
is a scalar replicated to fill the rank-2 block its position demands — that is
how a `0` becomes a zero block. `ArrayFlatten[a, r]` flattens `r` level pairs of
a rank-`2r` array into a rank-`r` array (default `r = 2`).
