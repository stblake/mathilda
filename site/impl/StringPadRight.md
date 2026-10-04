---
source: src/strings/stringpad.c
---
**Algorithm.** `StringPadRight` is the mirror of `StringPadLeft` and shares the same `pad_dispatch` (invoked with the `left` flag false). It pads on the right and, when the string is longer than `n`, truncates keeping the first `n` bytes. The pad string (a single space by default) is laid down cyclically from the left, `p[i mod plen]`, so a multi-character pad such as `".-"` repeats `.-.-…`. A `List` first argument pads every element to the explicit `n`, or to the longest element in the one-argument form.

**Data structures.** One output buffer per string; a `List` result for a list input.

**Complexity / limits.** `O(n)` per string. An arity outside 1–3 emits `StringPadRight::argb`. A non-integer/negative `n`, a non-string or list-valued pad, or an empty pad when padding is required, leaves the call unevaluated. Not `Listable`; byte-length based.
