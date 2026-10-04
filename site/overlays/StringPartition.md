### Worked examples

```mathematica
In[1]:= StringPartition["abcdefgh", 3]  (* non-overlapping blocks; a short tail is dropped *)
```

```mathematica
In[1]:= StringPartition["abcdef", 2, 1]  (* length 2, offset 1: overlapping blocks *)
```

```mathematica
In[1]:= StringPartition["abcdefg", UpTo[3]]  (* allow a shorter final block *)
```

### Notes

`StringPartition[s, n, d]` gives length-`n` substrings starting every `d`
characters (`d` defaults to `n`). With `d > n` characters are skipped; with
`d < n` the blocks overlap. By default a trailing run too short to fill a block
is dropped.

`UpTo[n]` keeps that final short block so every character appears. Both `n` and
`d` must be positive integers; indexing is byte-based. A list of strings threads.
