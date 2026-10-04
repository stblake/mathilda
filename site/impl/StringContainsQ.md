---
source: src/strings/regex/stringcontainsq.c
---
**Algorithm.** The four substring predicates — `StringContainsQ`, `StringFreeQ`, `StringStartsQ`, `StringEndsQ` — share one core, `sq_dispatch`, and differ only by an `SqKind` tag. It seeds `IgnoreCase` from the registered `Options`, strips trailing option rules, and — with a single positional argument — returns the curried operator form `Function[head[#1, patt, opts...]]`. Otherwise it builds the rule set unanchored and asks `sq_match`: one `regex_match` per rule with an early exit on the first hit (a predicate never needs to know *where* or *how often*). `StringContainsQ` (`SQ_CONTAINS`) reports that match directly.

**Data structures.** A `RegexRule` array; only the whole-match pair (`ov[2]`) is requested — no `regex_scan` enumeration and no capture pool.

**Complexity / limits.** One PCRE2 probe per rule, short-circuited. It is equivalent to `!StringFreeQ[s, p]` and to `StringMatchQ[s, ___ ~~ patt ~~ ___]`. A non-string subject — or a non-string element of a subject list — leaves the *whole* call unevaluated rather than threading a wrong Boolean. Byte semantics.
