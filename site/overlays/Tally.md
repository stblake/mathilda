### Worked examples

```mathematica
In[1]:= Tally[{a, b, a, c, b, a}]
Out[1]= {{a, 3}, {b, 2}, {c, 1}}
```

```mathematica
In[1]:= Tally[Table[Mod[n^2, 5], {n, 0, 20}]]
Out[1]= {{0, 5}, {1, 8}, {4, 8}}
```

```mathematica
In[1]:= Tally[Table[GCD[n, 12], {n, 1, 12}]]
Out[1]= {{1, 4}, {2, 2}, {3, 2}, {4, 2}, {6, 1}, {12, 1}}
```

### Notes

`Tally[list]` returns `{element, count}` pairs for each distinct element, in the
order of first appearance. It is a compact way to read off the distribution of a
computed sequence — for example, the multiplicities of the quadratic residues
mod 5 (`{0, 1, 4}` appearing 5, 8, and 8 times among `n = 0 .. 20`), or the
divisor structure of `GCD[n, 12]`.
