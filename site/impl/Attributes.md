---
source: src/attr.c
---
`builtin_attributes` (`src/attr.c`) reads the symbol's attribute bitflags via `get_attributes(name)` and builds a `List` of attribute symbols from them, e.g. `ATTR_FLAT` -> `Flat`, the `HoldFirst|HoldRest` pair collapsing to `HoldAll`, `ATTR_LISTABLE` -> `Listable`, `ATTR_PROTECTED` -> `Protected`, etc. The symbol itself is held unevaluated (`Attributes` carries `ATTR_HOLDALL`).
