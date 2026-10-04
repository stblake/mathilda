---
source: src/strings/stringrepeat.c
---
**Algorithm.** `builtin_stringrepeat` accepts 2–3 arguments (else `StringRepeat::argt`). It computes the output length as `n · len`, optionally capped at `max`, then fills the buffer cyclically (`buf[i] = str[i mod len]`) so a truncated final copy falls out naturally. The `n · len` product is guarded against `size_t` overflow, which is tolerated only when a `max` cap keeps the result finite.

**Data structures.** One output buffer of the computed length.

**Complexity / limits.** `O(output)`. `n == 0`, an empty base string, or `max == 0` gives `""`; a non-string base, or a non-integer/negative count or `max`, leaves the call unevaluated. Byte-length based, consistent with the rest of the string subsystem.
