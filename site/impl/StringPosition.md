---
source: src/strings/regex/stringposition.c
---
**Algorithm.** `builtin_stringposition` shares the option-seeding and rule-building machinery of `StringCases`/`StringCount`. `sp_scalar` runs the shared `regex_scan` and turns each span `[ms, me)` into the 1-based inclusive pair `{ms + 1, me}` — the position form that `StringTake`, `StringDrop`, and `StringReplacePart` consume. An optional third integer argument caps the number of pairs returned.

It differs from its two siblings in one default: `Overlaps -> True` (matching the Wolfram Language), so overlapping matches at distinct start positions are listed, where `StringCases`/`StringCount` default to `False`.

**Data structures.** The shared `RegexScan` span array; the result is a `List` of two-element position `List`s. A list of subjects threads.

**Complexity / limits.** Shares the `regex_scan` enumerator, so counts and policies agree with `StringCases`/`StringCount`. Positions are byte offsets (no UTF-8 decoding), consistent with `StringLength`/`StringPart`.
