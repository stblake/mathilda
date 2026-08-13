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

## Examples

_No verified examples yet for this function._

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

**Attributes:** `Protected`.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
