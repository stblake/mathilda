---
source: src/strings/stringdrop.c
---
**Algorithm.** `builtin_stringdrop` requires exactly two arguments (else `StringDrop::argrx`). It allocates a byte keep-mask initialised to `true`, clears the positions named by the sequence spec, and rebuilds the kept bytes in order (`stringdrop_build_kept`) — the complementary "keep" mask is how `StringDrop` is expressed as the complement of `StringTake`. Specs handled: an integer `n` (drop the first `n`, or the last `|n|` when negative), `UpTo[n]` (clamped to the length), `{n}`, `{m, n}` (a decreasing range drops nothing), and `{m, n, s}` (stepped). Negative endpoints normalise as `len + k + 1`. A `List` first argument recurses per element via `evaluate(StringDrop[si, spec])`.

**Data structures.** A `bool[len]` keep-mask and one output buffer.

**Complexity / limits.** `O(len)`. An out-of-range position, a zero step, or a non-integer spec leaves the call unevaluated; indexing is byte-based (no UTF-8 decoding), consistent with `StringTake`.
