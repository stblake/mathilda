# SplitBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SplitBy[list, f]`**

splits list into runs of consecutive elements that give the same value of f\[element\]. Only adjacent elements are grouped (unlike GatherBy, which collects equal keys from anywhere in the list).

**`SplitBy[list, {f1, f2, ...}]`**

splits by f1, then splits each resulting run by f2, and so on, nesting one level deeper per function.

## Examples

_No verified examples yet for this function._

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

**Attributes:** `Protected`.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
