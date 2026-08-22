# SetAttributes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetAttributes[s, attr] sets the attributes for s.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= SetAttributes[g, Orderless]
Out[1]= Null

In[2]:= Attributes[g]
Out[2]= {Orderless}

In[3]:= g[3, 1, 2]
Out[3]= g[1, 2, 3]
```

## Implementation notes

`builtin_set_attributes` (`src/attr.c`) takes a symbol spec and an attribute spec and
OR-folds the named attribute bits into the target symbol's `SymbolDef` flag word. The
symbol spec may be a single symbol/string or a `List` of them, in which case the same
attribute spec is applied to each via `set_attributes_for_symbol`. Attribute names map to
the `ATTR_*` bitflags defined in `attr.h` (e.g. `HoldAll`, `Flat`, `Orderless`,
`Listable`, `Protected`); the helper accepts either a single attribute symbol or a `List`
of them. The handler carries `HoldFirst` so the symbol argument is not evaluated before
its attributes are read, and returns `Null`. Reading attributes back is the inverse
`builtin_attributes`, which decodes the same bitflags into a sorted `List`.

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_clearall_remove_protect.c`](https://github.com/stblake/mathilda/blob/main/tests/test_clearall_remove_protect.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_eval_timestamps.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_timestamps.c)

## Notes & additional examples

### Notes

`SetAttributes[s, attr]` adds an attribute to `s`. Here `Orderless` makes the evaluator sort `g`'s arguments into canonical order on every call.
