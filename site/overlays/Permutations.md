### Worked examples

All orderings of a list:

```mathematica
In[1]:= Permutations[{a, b, c}]
Out[1]= {{a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}}
```

The number of permutations of `n` distinct elements is `n!`:

```mathematica
In[1]:= Length[Permutations[Range[6]]]
Out[1]= 720
```

Repeated elements yield only the *distinct* arrangements (no duplicates), so a
multiset gives multinomial-many results rather than `n!`:

```mathematica
In[1]:= Permutations[{1, 1, 2}]
Out[1]= {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}}
```

The `{n}` form gives the ordered length-`n` arrangements (k-permutations) — here
all ordered pairs drawn from four symbols:

```mathematica
In[1]:= Permutations[{1, 2, 3, 4}, {2}]
Out[1]= {{1, 2}, {1, 3}, {1, 4}, {2, 1}, {2, 3}, {2, 4}, {3, 1}, {3, 2}, {3, 4}, {4, 1}, {4, 2}, {4, 3}}
```

`Permutations` is a building block for combinatorial search: filtering the full
list by a predicate enumerates structured arrangements. Keeping only those
permutations of `{1, 2, 3, 4}` whose first element is below the second selects
exactly half of them:

```mathematica
In[1]:= Select[Permutations[{1, 2, 3, 4}], (#[[1]] < #[[2]] &)]
Out[1]= {{1, 2, 3, 4}, {1, 2, 4, 3}, {1, 3, 2, 4}, {1, 3, 4, 2}, {1, 4, 2, 3}, {1, 4, 3, 2}, {2, 3, 1, 4}, {2, 3, 4, 1}, {2, 4, 1, 3}, {2, 4, 3, 1}, {3, 4, 1, 2}, {3, 4, 2, 1}}
```

### Notes

`Permutations[list]` lists every ordering of the elements; for `n` distinct
elements there are `n!` of them, returned in canonical order. Repeated elements
collapse to distinct arrangements only. `Permutations[list, n]` bounds the length
to at most `n`, and `Permutations[list, {n}]` gives exactly length `n` (the
ordered k-permutations). Combined with `Select` / `OrderedQ` it drives small
combinatorial searches directly in the language.
