---
source: src/strings/stringpartition.c
---
**Algorithm.** `builtin_stringpartition` accepts 2–3 arguments (else `StringPartition::argt`). The block length `n` comes from an integer (every block exactly `n`) or from `UpTo[n]` (a short final block allowed); the offset `d` defaults to `n`. It walks `start = 0, d, 2d, …`, emitting `subj[start, start+n)` while the full block fits, or the short tail `subj[start, len)` when `UpTo` permits, and stops otherwise. With `d > n` characters are skipped; with `d < n` the blocks overlap. A `List` first argument threads via `evaluate`.

**Data structures.** An `Expr*` block array bounded by `len/d + 2`; each block is a fresh `EXPR_STRING` (`sp_block` `memcpy`s the byte range).

**Complexity / limits.** `O(output)`. Both `n` and `d` must be positive integers; a non-integer or non-positive value, or a non-string subject, leaves the call unevaluated. Byte-indexed, consistent with `StringTake`/`StringDrop`.
