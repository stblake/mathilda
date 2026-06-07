---
source: src/attr.c
---
`builtin_set_attributes` (`src/attr.c`) takes a symbol spec and an attribute spec and
OR-folds the named attribute bits into the target symbol's `SymbolDef` flag word. The
symbol spec may be a single symbol/string or a `List` of them, in which case the same
attribute spec is applied to each via `set_attributes_for_symbol`. Attribute names map to
the `ATTR_*` bitflags defined in `attr.h` (e.g. `HoldAll`, `Flat`, `Orderless`,
`Listable`, `Protected`); the helper accepts either a single attribute symbol or a `List`
of them. The handler carries `HoldFirst` so the symbol argument is not evaluated before
its attributes are read, and returns `Null`. Reading attributes back is the inverse
`builtin_attributes`, which decodes the same bitflags into a sorted `List`.
