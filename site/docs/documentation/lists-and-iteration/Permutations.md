# Permutations

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Permutations[list]
    generates a list of all possible permutations of the elements in list.
Permutations[list,n]
    gives all permutations containing at most n elements.
Permutations[list,{n}]
    gives all permutations containing exactly n elements.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Permutations[{a, b, c}]
Out[1]= {{a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a}}

In[2]:= Permutations[{a, b, c, d}, {3}]
Out[2]= {{a, b, c}, {a, b, d}, {a, c, b}, {a, c, d}, {a, d, b}, {a, d, c}, {b, a, c}, {b, a, d}, {b, c, a}, {b, c, d}, {b, d, a}, {b, d, c}, {c, a, b}, {c, a, d}, {c, b, a}, {c, b, d}, {c, d, a}, {c, d, b}, {d, a, b}, {d, a, c}, {d, b, a}, {d, b, c}, {d, c, a}, {d, c, b}}

In[3]:= Permutations[{a, a, b}]
Out[3]= {{a, a, b}, {a, b, a}, {b, a, a}}

In[4]:= Permutations[{x, x^2, x + 1}]
Out[4]= {{x, x^2, 1 + x}, {x, 1 + x, x^2}, {x^2, x, 1 + x}, {x^2, 1 + x, x}, {1 + x, x, x^2}, {1 + x, x^2, x}}

In[5]:= Permutations[Range[3], All]
Out[5]= {{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 1}, {2, 3}, {3, 1}, {3, 2}, {1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1}, {3, 1, 2}, {3, 2, 1}}

In[6]:= Permutations[Range[4], {4, 0, -2}]
Out[6]= {{1, 2, 3, 4}, {1, 2, 4, 3}, {1, 3, 2, 4}, {1, 3, 4, 2}, {1, 4, 2, 3}, {1, 4, 3, 2}, {2, 1, 3, 4}, {2, 1, 4, 3}, {2, 3, 1, 4}, {2, 3, 4, 1}, {2, 4, 1, 3}, {2, 4, 3, 1}, {3, 1, 2, 4}, {3, 1, 4, 2}, {3, 2, 1, 4}, {3, 2, 4, 1}, {3, 4, 1, 2}, {3, 4, 2, 1}, {4, 1, 2, 3}, {4, 1, 3, 2}, {4, 2, 1, 3}, {4, 2, 3, 1}, {4, 3, 1, 2}, {4, 3, 2, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 1}, {2, 3}, {2, 4}, {3, 1}, {3, 2}, {3, 4}, {4, 1}, {4, 2}, {4, 3}, {}}

In[7]:= Permutations[f[a, b, c]]
Out[7]= {f[a, b, c], f[a, c, b], f[b, a, c], f[b, c, a], f[c, a, b], f[c, b, a]}
```

## Implementation notes

**Algorithm.** `builtin_permutations` (in `src/funcprog.c`) generates distinct permutations, correctly handling repeated elements. It first compresses the input into a `UniqueElement[]` multiset — each distinct value (compared with `expr_eq`) paired with its multiplicity `count`. The recursive `permutations_rec` then builds permutations by, at each position, trying every unique element that still has remaining count, decrementing it before recursing and restoring it afterward (classic count-bounded backtracking). Because it only ever places each *distinct* value once per position, it emits each permutation exactly once even with duplicates (so `Permutations[{1,1,2}]` gives 3 results, not 6). The enumeration order is the order in which distinct elements first appear in the input.

The optional second argument selects subsequence lengths: an integer `n` (length-`n` arrangements), `All`, or a `{min}`/`{min,max}`/`{min,max,step}` range; the handler loops over the requested lengths `l` and calls `permutations_rec` with `target_len = l`. The permutation head is inherited from the input; results accumulate in a geometrically grown `Expr**` buffer wrapped as `List[...]`.

- `Protected`.
- There are $n!$ permutations of a list of $n$ distinct elements.
- Repeated elements are treated as identical.
- The object `list` need not have head `List`.
- `Permutations[list]` is effectively equivalent to `Permutations[list, {Length[list]}]`.
- `Permutations[list, {n_min, n_max}]` gives permutations of `list` between `n_min` and `n_max` elements.
- `Permutations[list, {n_min, n_max, dn}]` uses step `dn`.
- `Permutations[list, All]` is equivalent to `Permutations[list, {0, Length[list]}]`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)

## Notes & additional examples

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
