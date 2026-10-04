# Subsets

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subsets[list]`**

Gives all subsets of list (the power set), ordered by increasing length and lexicographically by element position within each length. The head of list is kept on the subsets.

**`Subsets[list, n]`**

Gives subsets of length 0 through n.

**`Subsets[list, {n}]`**

Gives subsets of length exactly n.

**`Subsets[list, {nmin, nmax}]`**

Gives subsets whose length lies in the inclusive range nmin to nmax; a third element gives a length step.

**`Subsets[list, spec, s]`**

Gives only the first s subsets spec would produce, generated lazily.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The full power set, by increasing length

```mathematica
In[1]:= Subsets[{a, b, c}]
Out[1]= {{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, c}, {a, b, c}}
```

Exactly the length-2 subsets

```mathematica
In[2]:= Subsets[{a, b, c, d}, {2}]
Out[2]= {{a, b}, {a, c}, {a, d}, {b, c}, {b, d}, {c, d}}
```

Lengths 0 through 2

```mathematica
In[3]:= Subsets[{1, 2, 3}, 2]
Out[3]= {{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}}
```

## Algorithm

Subsets — enumerate the sublists of an expression.

Mathematica semantics:

```text
  Subsets[list]                 the power set, ordered by increasing length
                                and lexicographically by original element
                                position within each length:
                                Subsets[{a,b,c}] ->
                                  {{}, {a}, {b}, {c}, {a,b}, {a,c}, {b,c},
                                   {a,b,c}}
  Subsets[list, n]              lengths 0 through n inclusive
  Subsets[list, {n}]            exactly length n
  Subsets[list, {nmin, nmax}]   lengths nmin..nmax inclusive
  Subsets[list, {nmin, nmax, d} lengths nmin, nmin+d, ... up to nmax
  Subsets[list, spec, s]        only the first s subsets the spec produces
```

The head of the inner sublists is taken from the input expression, so Subsets[f[a,b]] gives {f[], f[a], f[b], f[a,b]}. The outer wrapper is always a List. Duplicate elements are treated as distinct by position: no dedup is performed, hence Subsets[{a,a}] -> {{}, {a}, {a}, {a,a}}.

### Performance

The full result is exponential in Length[list], so the generator is lazy: it walks index combinations with an odometer and stops the instant the `s` budget is exhausted. Subsets[<30 elements>, All, 5] therefore costs five sublist allocations, not 2^30. Output storage grows geometrically and is never sized from the theoretical subset count (which would overflow for moderate lengths anyway).

## Implementation notes

**Algorithm.** `Subsets[list]` (`src/list/subsets.c`) generates all subsets of `list` — the
power set — ordered first by increasing length and then lexicographically by original element
position within each length. The length can be restricted: `Subsets[list, n]` gives lengths
`0` through `n`, `Subsets[list, {n}]` exactly `n`, `Subsets[list, {nmin, nmax}]` the
inclusive range (a third element in the spec gives a length step). `Subsets[list, spec, s]`
returns only the first `s` subsets the spec would produce.

**Data structures.** Ordinary `Expr` trees. Each subset is emitted as a copy of the list's
elements at a chosen set of positions, carrying the list's own head (so `Subsets[f[a, b]]`
keeps head `f`). Subsets of a given length are produced by an index-combination generator
that advances position tuples in lexicographic order.

**Complexity / limits.** The full power set has `2^Length[list]` members, exponential in the
input, so the generator is **lazy**: it produces subsets on demand and the `, s` cap lets a
caller take just the first `s` without materialising the rest — the natural way to peek at
the start of an otherwise intractable enumeration. Order is deterministic (by length, then
lexicographic position), never sorted by value.

**Attributes:** `Protected`.

## References

- Source: [`src/list/subsets.c`](https://github.com/stblake/mathilda/blob/main/src/list/subsets.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)

## Notes & additional examples

### Notes

`Subsets[list]` gives every subset of `list`, ordered by increasing length and then
lexicographically by original element position — so the empty set comes first and the whole
list last. The head of `list` is kept on each subset.

The length spec restricts the output: `Subsets[list, n]` gives lengths `0`–`n`,
`Subsets[list, {n}]` exactly `n`, `Subsets[list, {nmin, nmax}]` an inclusive range. Because
the power set is exponential, the generator is lazy and `Subsets[list, spec, s]` returns only
the first `s` it would produce — a safe way to peek at the start of an enumeration too large
to build in full.
