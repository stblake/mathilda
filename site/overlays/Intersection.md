### Worked examples

```mathematica
In[1]:= Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]  (* the elements common to all, sorted and distinct *)
```

```mathematica
In[2]:= Intersection[Range[1, 10], Range[5, 15]]  (* two integer ranges overlap on 5..10 *)
```

```mathematica
In[3]:= Intersection[Divisors[60], Divisors[45]]  (* common divisors *)
```

```mathematica
In[4]:= Intersection[{1.1, 3.4, 0.5, 7.6, 1.9}, {1.2, 3.3, 7.7}, SameTest -> (Floor[#1] == Floor[#2] &)]  (* a custom equality test *)
```

### Notes

`Intersection[l1, l2, ...]` gives the sorted list of distinct elements common to
every operand; `Intersection[list]` alone is the sorted distinct elements of one
list. All operands must share a head, which need not be `List` (so
`Intersection[f[a, b], f[c, a]]` is `f[a]`). The default comparison is canonical
structural equality, done in `O(total)` with a hash set; `SameTest -> f` switches
to an `O(n^2)` path keeping the canonically-greatest member of each class. A
rank-1 buffer of exact integers — whether an invisible packed list or an explicit
`NDArray` — takes a machine sorted-merge fast path and keeps its representation;
reals are routed through the general path. `Flat`, `OneIdentity`, `Protected`.
