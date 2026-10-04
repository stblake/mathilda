---
source: src/numberform.c
---
**Algorithm.** `Row` is a display wrapper like `NumberForm`: `builtin_row`
(`src/numberform.c`) is inert (returns `NULL`, so the head stays in the tree) and
the rendering lives in the standard printer. When `print.c` meets
`Row[{e1, e2, ...}]` (first argument a `List`, and only outside `InputForm`) it
sets the output-form flag and prints each element with `print_standard`, so
strings appear **without quotes**; a second string argument `Row[{...}, s]` is
emitted between successive elements as a separator. Outside that shape the literal
`Row[...]` prints normally.

**Data structures.** No allocation of its own — it writes directly to `stdout`
through the printer, toggling and restoring the global `g_print_output_form`
flag around the element loop. Its chief internal use is as the assembler for
`NumberForm`'s `NumberFormat` option, where `numberform_format_result_to_string`
renders a `Row` of mantissa / base / exponent pieces to an OutputForm string.

**Complexity / limits.** Linear in the number of elements. `Row` is a pure
presentation head: it has no computational value and takes the list's first
argument as its content, with an optional string separator. Attributes
`Protected`.
