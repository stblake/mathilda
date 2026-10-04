---
source: src/strings/stringtrim.c
---
**Algorithm.** `builtin_stringtrim` accepts 1–2 arguments (else `StringTrim::argt`). The pattern (default `Whitespace`) is translated once by `wl_pattern_to_regex`, then compiled into two one-end-anchored programs (`st_compile`): a front matcher `(?:src)` and a back matcher `(?:src)\z`. `st_scalar` strips the leading run by repeated matches that must *begin exactly* at the current start (checked by `ov[0] == start`, since `\A`/`^` would anchor to absolute 0), and the trailing run by `\z`-anchored matches over a truncated subject length. Each end is stripped to a fixed point, and a zero-width match ends the loop — so only the two ends are trimmed, never the middle.

**Data structures.** Two `RegexProgram` matchers (front, back); one output buffer for the trimmed window `[start, end)`.

**Complexity / limits.** `O(len)` plus the cost of the trim-run matches. A single-character pattern still removes a whole run (fixed-point per end). A non-string subject, an unsupported pattern, or a build without PCRE2 leaves the call unevaluated. Byte-based; the pattern accepts the full shared string-pattern vocabulary.
