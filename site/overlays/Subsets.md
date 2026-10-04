### Worked examples

```mathematica
In[1]:= Subsets[{a, b, c}]  (* the full power set, by increasing length *)
```

```mathematica
In[1]:= Subsets[{a, b, c, d}, {2}]  (* exactly the length-2 subsets *)
```

```mathematica
In[1]:= Subsets[{1, 2, 3}, 2]  (* lengths 0 through 2 *)
```

### Notes

`Subsets[list]` gives every subset of `list`, ordered by increasing length and then
lexicographically by original element position — so the empty set comes first and the whole
list last. The head of `list` is kept on each subset.

The length spec restricts the output: `Subsets[list, n]` gives lengths `0`–`n`,
`Subsets[list, {n}]` exactly `n`, `Subsets[list, {nmin, nmax}]` an inclusive range. Because
the power set is exponential, the generator is lazy and `Subsets[list, spec, s]` returns only
the first `s` it would produce — a safe way to peek at the start of an enumeration too large
to build in full.
