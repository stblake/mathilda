---
source: src/strings/stringreverse.c
---
**Algorithm.** `builtin_stringreverse` requires one argument (else `StringReverse::argx`) and copies the subject's bytes into a fresh buffer in reverse order. `StringReverse` is `Listable`, so the evaluator threads it element-wise over a list before the builtin runs; each call therefore sees a single string.

**Data structures.** One output buffer of the same length as the input.

**Complexity / limits.** `O(len)`, byte-oriented (characters are single `char`s, no UTF-8 decoding). A single non-string argument leaves the call unevaluated so symbolic arguments flow through unchanged.
