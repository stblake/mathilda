---
source: src/strings/stringinsert.c
---
**Algorithm.** `builtin_stringinsert` requires three arguments (else `StringInsert::argrx`). `si_pos_to_offset` maps each position against the *original* string — a positive `n` to offset `n - 1` (insert before character `n`), a negative `-k` to `len + n + 1` — rejecting `0` and out-of-range values. The positions are tallied into a `counts[0..len]` array, then a single pass builds the output, emitting `counts[i]` copies of the insert string before original byte `i`. Because all positions are resolved up front, every copy lands relative to the untouched input. A `List` first argument recurses per element.

**Data structures.** An `int64 counts[len + 1]` tally and one output buffer sized `len + npos·inslen`.

**Complexity / limits.** `O(len + output)`. A non-string subject or insert string, a non-integer position, or an out-of-range position leaves the call unevaluated. Positions are valid for `1 ≤ n ≤ len + 1` (and the negative mirror); byte-based.
