# OwnValues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OwnValues[s] gives a list of own-value rules for s.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= OwnValues[a]
Out[1]= {}

In[2]:= a = 7
Out[2]= 7

In[3]:= OwnValues[a]
Out[3]= {7 -> 7}
```

## Implementation notes

`builtin_own_values` (1-arg, symbol only) calls `symtab_get_own_values(sym)` to retrieve the symbol's `Rule*` list — the immediate `x = value` assignments — and hands it to the shared `rules_to_list` helper, which deep-copies each `pattern`/`replacement` pair into a `Rule[lhs, rhs]` node and gathers them into a `List`. The list is returned unevaluated so the stored values are reported verbatim rather than re-evaluated.

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_unset.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unset.c)

## Notes & additional examples

### Notes

`OwnValues[s]` returns the direct value rules created by `s = ...`. An undefined symbol has an empty list; assigning `a = 7` stores a single own-value rule.
