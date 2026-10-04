---
source: src/strings/regex/stringcontainsq.c
---
**Algorithm.** `StringEndsQ` is the `SQ_ENDS` face of the shared `sq_dispatch` core — the mirror of `StringStartsQ`. It wraps the pattern as `StringExpression[patt, EndOfString]`, where the translator renders `EndOfString` as `\z` (an absolute end anchor), and hands that unanchored to the shared rule builder before the shared `sq_match` early-exits on the first hit.

**Data structures.** The synthesised wrapper `Expr` (outliving the rule set, as `RegexRule.lhs` borrows into it) and the `RegexRule` array, whole-match pair only.

**Complexity / limits.** Equivalent to `StringContainsQ[s, patt ~~ EndOfString]`. `IgnoreCase` option; curried operator form `StringEndsQ[patt]`. A non-string subject leaves the call unevaluated; byte semantics.
