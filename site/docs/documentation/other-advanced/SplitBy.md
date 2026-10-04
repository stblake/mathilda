# SplitBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SplitBy[list, f]`**

splits list into runs of consecutive elements that give the same value of f\[element\]. Only adjacent elements are grouped (unlike GatherBy, which collects equal keys from anywhere in the list).

**`SplitBy[list, {f1, f2, ...}]`**

splits by f1, then splits each resulting run by f2, and so on, nesting one level deeper per function.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Runs of equal parity, adjacency only

```mathematica
In[1]:= SplitBy[{1, 3, 5, 2, 4, 7}, EvenQ]
Out[1]= {{1, 3, 5}, {2, 4}, {7}}
```

A new run starts each time the key changes

```mathematica
In[2]:= SplitBy[{1, 2, 4, 3, 5}, EvenQ]
Out[2]= {{1}, {2, 4}, {3, 5}}
```

## Algorithm

SplitBy[list, f] — split a list into runs of *consecutive* elements that share the same value of f[element].

This is the key-function counterpart of Split (src/list/split.c): Split compares adjacent elements directly (or via a two-argument test), whereas SplitBy compares the evaluated keys f[e]. Only adjacent elements are ever grouped, which is what distinguishes SplitBy from GatherBy (src/assoc.c) — GatherBy collects *all* elements sharing a key, no matter where they sit.

```text
  SplitBy[{1, 3, 2, 4, 5}, EvenQ]     -> {{1, 3}, {2, 4}, {5}}
  SplitBy[{1, 2, 3, 4, 5, 6}, EvenQ]  -> {{1}, {2}, {3}, {4}, {5}, {6}}
  SplitBy[{1, 1, 2, 2, 3}, Identity]  -> {{1, 1}, {2, 2}, {3}}
```

The list form SplitBy[list, {f1, f2, ...}] splits by f1, then splits each resulting run by f2, and so on, nesting one level deeper per function:

```text
  SplitBy[{1, 3, 2, 4}, {EvenQ}]      -> {{1, 3}, {2, 4}}
```

Cost: f is evaluated exactly once per element per level, i.e. O(n) calls per function in the key spec, plus O(n) structural copying. Keys are compared with expr_eq — the same structural equality the rest of the kernel uses — so two adjacent elements whose keys stay unevaluated but identical still group together. Only one key is held live at a time (the previous element's), so peak overhead beyond the result itself is a single key plus the per-level run vector.

## Implementation notes

**Algorithm.** `SplitBy[list, f]` (`src/list/splitby.c`) walks `list` once, evaluating the
key `f[element]` for each element and starting a new run whenever the key differs from the
previous element's key. Only **adjacent** elements are ever grouped — this is the difference
from `GatherBy`, which collects equal keys from anywhere in the list, and from `Split`, which
compares adjacent elements (or a two-argument test) directly rather than through an evaluated
key. `SplitBy[list, {f1, f2, ...}]` splits by `f1`, then recursively splits each resulting
run by `f2`, nesting one level deeper per function.

**Data structures.** Ordinary `Expr` trees. The walk keeps the current run as a growing
`Expr**` vector and the previous key for the adjacency comparison; each completed run is
copied out under a fresh copy of the list's head (`splitby_make_run`). Keys are compared
structurally, so two adjacent elements whose keys evaluate to the same unevaluated form still
group together. The list's own head is preserved on the outer result and on each run.

**Complexity / limits.** `O(n)` key evaluations and comparisons for a single function; the
`{f1, ..., fk}` form multiplies by the nesting depth. `SplitBy` preserves element order
within every run (it never sorts), so the runs read left-to-right exactly as the input does.

**Attributes:** `Protected`.

## References

- Source: [`src/list/splitby.c`](https://github.com/stblake/mathilda/blob/main/src/list/splitby.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`SplitBy[list, f]` splits `list` into runs of consecutive elements that share the same value
of `f[element]`. Only **adjacent** elements are grouped: in the first example the leading
odds `{1, 3, 5}` form one run, the evens `{2, 4}` the next, and the trailing odd `7` its own
run — unlike `GatherBy`, which would collect all the odds together regardless of position.

It differs from `Split` in comparing the evaluated keys `f[e]` rather than the elements
themselves. `SplitBy[list, {f1, f2, ...}]` splits by `f1`, then splits each run by `f2`, and
so on, nesting one level deeper per function. Element order is preserved within every run.
