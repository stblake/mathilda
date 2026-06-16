### Worked examples

Test whether a list is in canonical order:

```mathematica
In[1]:= OrderedQ[{1, 2, 3}]
Out[1]= True

In[2]:= OrderedQ[{3, 1, 2}]
Out[2]= False
```

A custom ordering function as the second argument — here descending order via
`Greater`:

```mathematica
In[1]:= OrderedQ[{5, 4, 3, 2, 1}, Greater]
Out[1]= True
```

Because `OrderedQ` works on any head and accepts an arbitrary comparator, it
composes with `Select` to filter combinatorial output. Selecting the
canonically-ordered permutations of `{1, 2, 3}` recovers exactly the single
sorted tuple — the identity permutation:

```mathematica
In[1]:= Select[Permutations[{1, 2, 3}], OrderedQ]
Out[1]= {{1, 2, 3}}
```

A pure-function comparator orders by a derived key — these symbolic powers are
in strictly decreasing degree:

```mathematica
In[1]:= OrderedQ[{x^2, x, 1}, (#1 > #2 &)]
Out[1]= True
```

### Notes

`OrderedQ[h[e1, e2, ...]]` returns `True` exactly when each adjacent pair is in
canonical order under `expr_compare` (the same total order the evaluator uses
for `Orderless` heads and `Sort`). Ties (equal elements) count as ordered, so a
list with repeats can still be `OrderedQ`. The two-argument form
`OrderedQ[expr, p]` replaces the comparator with any predicate `p[a, b]`,
letting you test for descending order, key-based order, or domain-specific
orderings.
