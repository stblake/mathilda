---
source: src/strings/stringsplit.c
---
**Algorithm.** `builtin_stringsplit` strips trailing `IgnoreCase` options, accepts 1–3 positional arguments, and compiles the delimiter pattern (default `Whitespace`) to PCRE through the shared engine (`regex_rules_build_ex`). `ss_scalar` scans left to right, pushing each inter-delimiter substring (interior empties kept) and, for a `patt -> val` rule, the inserted value — `make_insertion` expands a string RHS as a `$`-template and binds a named `x:patt :> f[x]` by a scoped `ReplaceAll`. A third argument `n` caps the piece count; `All` keeps the leading/trailing empties that are otherwise trimmed. The empty delimiter `""` short-circuits to a per-byte split (`ss_null`). A list of subjects threads.

**Data structures.** A `PieceVec` of `(Expr, is_sub)` records — the `is_sub` flag marks an empty-trimmable substring versus an inserted rule value — plus the per-match ovector.

**Complexity / limits.** One left-to-right pass with a PCRE2 probe per rule per position; a zero-width delimiter advances by one so the scan progresses. Offsets are byte offsets. The delimiter accepts the full shared string-pattern vocabulary (literals, `RegularExpression`, character classes, `~~`, `|`, `..`, `Except`, …).
