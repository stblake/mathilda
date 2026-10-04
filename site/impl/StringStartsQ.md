---
source: src/strings/regex/stringcontainsq.c
---
**Algorithm.** `StringStartsQ` is the `SQ_STARTS` face of the shared `sq_dispatch` core. Rather than thread a new anchor mode through the shared rule builder, it synthesises a temporary `StringExpression[StartOfString, patt]` and hands that *unanchored* to `regex_rules_build_ex`; the translator renders `StartOfString` as `\A`, an absolute anchor that pins the match to offset 0 even in an unanchored search. The translator's `group_join` wraps every child in `(?:...)`, so an alternation-bearing pattern anchors correctly as `(?:(?:\A)(?:a|b))` rather than the broken `\Aa|b`. Then the shared `sq_match` early-exits on the first hit.

**Data structures.** The synthesised wrapper `Expr` — which must outlive the rule set, since `RegexRule.lhs` borrows into it — plus the `RegexRule` array (whole-match pair only).

**Complexity / limits.** Equivalent to `StringContainsQ[s, StartOfString ~~ patt]`. `IgnoreCase` option; curried operator form `StringStartsQ[patt]`. A non-string subject leaves the call unevaluated; byte semantics.
