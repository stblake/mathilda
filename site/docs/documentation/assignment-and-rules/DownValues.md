# DownValues

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DownValues[s] gives a list of down-value rules for s.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= square[n_] := n*n
Out[1]= Null

In[2]:= DownValues[square]
Out[2]= {n_^2 -> n^2}

In[3]:= square[7]
Out[3]= 49
```

## Implementation notes

`builtin_down_values` (1-arg, symbol only) calls `symtab_get_down_values(sym)` to fetch the raw `Rule*` linked list held on the `SymbolDef`, then `rules_to_list` walks it, wrapping each node's `pattern`/`replacement` pair (deep-copied) in a `Rule[lhs, rhs]` and collecting them into a `List`. The rules are returned unevaluated — they are emitted as literal `RuleDelayed`-style `Rule` heads without re-evaluating the stored patterns. DownValues are stored newest-first, which is the traversal order here.

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
- Tests: [`tests/test_list_set.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list_set.c)
- Tests: [`tests/test_unset.c`](https://github.com/stblake/mathilda/blob/main/tests/test_unset.c)

## Notes & additional examples

### Notes

`DownValues[s]` returns the pattern rules defined on `s` via `f[args] := ...`. Each entry is the stored `lhs -> rhs` rule the evaluator tries when `s` is called.
