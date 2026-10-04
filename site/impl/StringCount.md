---
source: src/strings/regex/stringcount.c
---
**Algorithm.** `builtin_stringcount` handles options and builds its rule set exactly as `StringCases` does, but `sct_scalar` calls `regex_scan` with `want_captures = 0` and returns the span count as an `Integer` — it never materialises a substring, costing one small span record per match instead of a `malloc` + `Expr` + list element. A `Rule`/`RuleDelayed` pattern is accepted (only the LHS matters to a count), as is a list of subjects, which threads to one count each.

**Data structures.** The shared `RegexScan` (span array only, no capture pool); one `Integer` result per subject.

**Complexity / limits.** Because the scan is the same, `StringCount[s, p, opts]` is always exactly `Length[StringCases[s, p, opts]]` and `Length[StringPosition[...]]`. `Overlaps`/`IgnoreCase` behave as in `StringCases`. A wrong positional arity emits `StringCount::argrx`; byte semantics throughout.
