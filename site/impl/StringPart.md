---
source: src/picostrings.c
---
`builtin_stringpart` takes `(string, spec)` and indexes by byte (1-based; negative counts from the end via `len + k + 1`). The single-index path uses the helper `stringpart_single`, which bounds-checks `k` and returns a length-1 `EXPR_STRING`. A `List` spec maps `stringpart_single` over each index into a result `List`; a `Span[m, n, s]` spec resolves `start`/`end`/`step` (honouring `All` and negative endpoints), computes the element count, and emits the stepped characters as a `List`. When the first argument is itself a `List` of strings, the builtin recurses by constructing and evaluating an inner `StringPart[si, spec]` per element. Out-of-range or non-integer indices return `NULL`. `ATTR_PROTECTED`.
