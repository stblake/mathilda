---
source: src/strings/regex/stringcases.c
---
**Algorithm.** `builtin_stringcases` seeds `Overlaps`/`IgnoreCase` from the registered `Options` (so `SetOptions` takes effect), strips trailing option rules to leave the positional arguments, and builds the rule set — anchored only in `Overlaps -> All` mode, where exact-substring matches are enumerated. It then calls the shared enumerator `regex_scan` once per subject. For a bare pattern each span becomes the matched substring; for a `patt -> rhs` rule, `regex_rule_replacement` expands `$0`/`$n` from the capture pool. `want_captures` is set only when some rule carries an RHS, so pure extraction never allocates a capture pool.

**Data structures.** `regex_scan` fills a `RegexScan`: a growable span array plus a flat capture pool that spans index by *offset* (so growing the pool never invalidates a recorded span). The result is a `List` of `EXPR_STRING`.

**Complexity / limits.** The `False`/`True` scans stream left-to-right; `All` mode probes `O(len² · nr)` substrings against the anchored rules. `StringCount` and `StringPosition` share this scanner, so the three always agree on count and overlap policy. Offsets are byte offsets; the occurrence-limit form `StringCases[s, p, n]` is unsupported and left unevaluated.
