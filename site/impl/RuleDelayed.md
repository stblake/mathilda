---
source: src/replace.c
---
`RuleDelayed[lhs, rhs]` (`:>`) is a passive delayed-rewrite object with no builtin handler. Its `ATTR_HOLDREST | ATTR_SEQUENCEHOLD | ATTR_PROTECTED` attributes hold `rhs` unevaluated at construction; the rule engine (`is_rule` in `src/replace.c`) detects the `RuleDelayed` head and, each time the rule fires, substitutes the fresh pattern bindings into the held `rhs` and only then evaluates it (see the `delayed` flag threaded through `ReplaceRule`/`ReplaceListState`). This is the difference from `Rule`, whose RHS is evaluated once up front.
