---
source: src/attr.c
---
`builtin_clear_attributes` (`src/attr.c`) clears the bitflags named in its second argument from the target symbol(s) via `clear_attributes_for_symbol`. The first argument may be one symbol/string or a `List` of them; it returns `Null`. `ClearAttributes` carries `ATTR_HOLDFIRST` so the symbol is not evaluated first.
