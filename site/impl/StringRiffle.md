---
source: src/strings/stringriffle.c
---
**Algorithm.** `builtin_stringriffle` takes the first argument as the data and the rest as per-level separators (level 1 outermost). Each separator is parsed by `parse_sep` into a `{left, sep, right}` triple — a plain string becomes `{"", sep, ""}`, a 3-string list is the triple itself. `riffle_build` recurses: a leaf renders through `leaf_to_str` (strings verbatim, any other expression via `expr_to_string`/`ToString`), and each level joins its children with the explicit separator for that level or, when none is supplied, the default scheme chosen from `depth_from_bottom` — a single space at the innermost level and one extra newline per level above. It is the inverse of `StringSplit`.

**Data structures.** A `SepSpec` array of resolved triples; per level the child strings are built first, then assembled two-pass (sum lengths, then copy) into one buffer.

**Complexity / limits.** `O(total output)`. No arguments emits `StringRiffle::argm`; a first argument that is neither a list nor a string, or a malformed separator, leaves the call unevaluated. `StringRiffle` is deliberately not `Listable` so it can inspect the whole nested structure. Byte-oriented.
