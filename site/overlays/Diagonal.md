### Worked examples

```mathematica
In[1]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]  (* the leading diagonal *)
```

```mathematica
In[1]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]  (* the first superdiagonal, k > 0 above *)
```

```mathematica
In[1]:= Diagonal[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, -1]  (* the first subdiagonal, k < 0 below *)
```

```mathematica
In[1]:= Diagonal[{{1, 2, 3}, {4, 5, 6}}]  (* non-square: length is Min[rows, cols] *)
```

```mathematica
In[1]:= Diagonal[{{a, b}, {c, d}}]  (* symbolic entries are copied through verbatim *)
```

```mathematica
In[1]:= Diagonal[NDArray[{{1., 2.}, {3., 4.}}]]  (* a machine array rides the rank-2 buffer straight to rank-1 *)
```

### Notes

`Diagonal[m, k]` indexes the diagonal by its offset `k` from the leading one:
`k > 0` runs above it, `k < 0` below. An `|k|` beyond the matrix gives `{}`.

Because each entry `m[[i, i+k]]` is copied unchanged, `Diagonal` works on more than
square numeric matrices: it handles non-square inputs (the length is `Min[rows, cols]`)
and higher-rank tensors (the diagonal of a rank-`n` tensor is a rank-`(n-1)` structure).

A packed `List` or a visible `NDArray` takes the buffer fast path, which `memcpy`s the
diagonal into a fresh array of the same dtype — a packed input stays packed, a visible
`NDArray` stays visible.
