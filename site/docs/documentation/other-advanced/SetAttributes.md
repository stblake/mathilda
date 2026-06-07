# SetAttributes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetAttributes[s, attr] sets the attributes for s.
```

## Examples

_No verified examples yet for this function._

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

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]
Out[1]= Null

In[2]:= Attributes[g]
Out[2]= {Orderless}

In[3]:= g[3, 1, 2]
Out[3]= g[1, 2, 3]
```

### Notes

`SetAttributes[s, attr]` adds an attribute to `s`. Here `Orderless` makes the evaluator sort `g`'s arguments into canonical order on every call.
