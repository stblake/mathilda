---
source: src/replace.c
---
`Rule[lhs, rhs]` (`->`) is a passive rewrite-rule object, not a computation — it has no builtin handler. Because `Rule` carries no `Hold` attributes, the evaluator evaluates `rhs` (and `lhs`) when the rule expression is constructed; this is what makes `->` *immediate*. The rule engine (`is_rule` in `src/replace.c`) recognises `Rule`-headed nodes and, on a match in `ReplaceAll`/`Replace`/etc., substitutes the pattern bindings into the already-evaluated `rhs`. Contrast `RuleDelayed`. Attribute: `ATTR_PROTECTED`.
