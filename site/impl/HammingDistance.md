---
source: src/list/distance.c
---
**Algorithm.** `builtin_hamming_distance` counts how many positions differ
between two equal-length sequences. After `dist_seq_pair` confirms both arguments
are strings or both lists, `dist_seq` turns each into an element array (a string
explodes to one `Expr` per byte; a list's elements are borrowed) and a single
pass tallies the positions where `expr_eq` is false. Equal lengths are required:
on a length mismatch the call is left unevaluated (Mathematica raises `::idim`
there, so declining is the faithful behaviour).

**Encoding.** Strings are compared byte by byte, so a multi-byte UTF-8 character
spans several positions — matching the ASCII/DNA use cases and noted rather than
assumed. Comparison by `expr_eq` means lists of arbitrary expressions work, not
just numbers or characters.

**Complexity / limits.** O(n) comparisons, O(n) temporary for the string case.
`ATTR_PROTECTED`. See `EditDistance` for the insert/delete/substitute metric that
does not require equal lengths.
