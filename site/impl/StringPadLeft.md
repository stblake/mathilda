---
source: src/strings/stringpad.c
---
**Algorithm.** `StringPadLeft` and `StringPadRight` share `pad_dispatch`, a `left` flag selecting the side. `pad_one` either truncates — copying the last `n` bytes for the left variant, the first `n` for the right — or pads, laying the pad string down cyclically (`p[i mod plen]`) into the pad region that precedes (left) or follows (right) the string. The pad string defaults to a single space. For a `List` first argument the target length is the explicit `n`, or, in the one-argument form, the longest element — so a bare list form equalises all widths.

**Data structures.** One output buffer per string; a `List` result when the input is a list.

**Complexity / limits.** `O(n)` per string. An arity outside 1–3 emits `StringPadLeft::argb`/`StringPadRight::argb`. A non-integer or negative `n`, a non-string or list-valued pad string, or an empty pad when padding is required, leaves the call unevaluated. The pad builtins are deliberately **not** `Listable`, so the handler sees the whole list (the 1-arg form must). Byte-length based.
