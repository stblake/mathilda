---
source: src/strings/stringreplacepart.c
---
**Algorithm.** Two arguments is the operator form `StringReplacePart[new, part][old]`, realised as `Function[StringReplacePart[#1, new, part]]`; three is the direct form (any other arity emits `StringReplacePart::argt`). Each position is a `{m, n}` pair (`srp_parse_range`, negatives → `len + k + 1`), validated against the original string. A single new string is broadcast to every range; a list of new strings must match the range count. Ranges are accepted in order against a `covered` byte-mask — a later range touching an already-claimed byte triggers `StringReplacePart::ovlp` and is dropped — then insertion-sorted by start and assembled in two passes (compute length, then emit literal spans and replacements). A `List` first argument recurses per element.

**Data structures.** An `SrpRange` array (original + resolved positions + borrowed replacement pointer), a `bool[len]` covered-mask, a `const char*` replacement array, and one output buffer.

**Complexity / limits.** `O(nranges · rangelen + len)`. Positions use the `StringPosition` form and refer to the original string; an empty replacement string deletes the selected characters. A malformed/out-of-range position or a new-string/position length mismatch leaves the call unevaluated. Byte-indexed.
