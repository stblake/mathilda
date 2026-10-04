### Worked examples

```mathematica
In[1]:= Complement[{a, b, c, d, e}, {a, c}, {d}]  (* the first list minus everything in the rest *)
```

```mathematica
In[2]:= Complement[Range[10], Range[2, 10, 2]]  (* the odd numbers in 1..10 *)
```

```mathematica
In[3]:= Complement[f[a, b, c, d], f[c, a], f[b, b, a]]  (* any shared head works *)
```

```mathematica
In[4]:= Complement[{1.1, 3.4, 0.5, 7.6, 1.9}, {1.2, 3.3}, SameTest -> (Floor[#1] == Floor[#2] &)]  (* a custom equality test *)
```

### Notes

`Complement[eall, e1, e2, ...]` gives the sorted distinct elements of `eall` that
appear in none of the later lists — a set difference. Unlike `Union` and
`Intersection` it is order-sensitive in its first argument, so it is not `Flat`.
All operands must share a head, which need not be `List`. The default comparison
is canonical structural equality, computed in `O(total)` with a hash set;
`SameTest -> f` uses `f[a, b] === True` as the equivalence relation (keeping the
canonically-smallest member of each class), and `SameTest -> Automatic` is the
default. A rank-1 buffer of exact integers (packed list or explicit `NDArray`)
takes a machine sorted-set-difference fast path and keeps its representation.
`Protected`.
