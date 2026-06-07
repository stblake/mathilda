---
source: src/picostrings.c
---
`builtin_stringtake` takes `(string, spec)` and slices by byte through the helper `stringtake_substring` (which `memcpy`s `str[start-1 .. end-1]` into a fresh buffer). A positive integer `n` takes the first `n` bytes; negative `-n` takes the last `n`; `UpTo[n]` (detected by `is_upto`) clamps to the available length. A `List` spec selects: `{n}` a single character, `{m, n}` a range, and `{m, n, s}` a stepped range built byte-by-byte into a `malloc`'d buffer. Negative endpoints normalise as `len + k + 1`; out-of-range or non-integer specs return `NULL`. A first argument that is a `List` of strings is handled by recursively building and evaluating `StringTake[si, spec]` per element. `ATTR_PROTECTED`.
