# DownValues

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
DownValues[s] gives a list of down-value rules for s.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_down_values` (1-arg, symbol only) calls `symtab_get_down_values(sym)` to fetch the raw `Rule*` linked list held on the `SymbolDef`, then `rules_to_list` walks it, wrapping each node's `pattern`/`replacement` pair (deep-copied) in a `Rule[lhs, rhs]` and collecting them into a `List`. The rules are returned unevaluated — they are emitted as literal `RuleDelayed`-style `Rule` heads without re-evaluating the stored patterns. DownValues are stored newest-first, which is the traversal order here.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= square[n_] := n*n
Out[1]= Null

In[2]:= DownValues[square]
Out[2]= {n_^2 -> n^2}

In[3]:= square[7]
Out[3]= 49
```

### Notes

`DownValues[s]` returns the pattern rules defined on `s` via `f[args] := ...`. Each entry is the stored `lhs -> rhs` rule the evaluator tries when `s` is called.
