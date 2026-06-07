# Rule

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs -> rhs or Rule[lhs, rhs]
    represents an immediate rewrite rule: rhs is evaluated when the
    rule object is constructed, then matched against lhs at use.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Rule[lhs, rhs]` (`->`) is a passive rewrite-rule object, not a computation — it has no builtin handler. Because `Rule` carries no `Hold` attributes, the evaluator evaluates `rhs` (and `lhs`) when the rule expression is constructed; this is what makes `->` *immediate*. The rule engine (`is_rule` in `src/replace.c`) recognises `Rule`-headed nodes and, on a match in `ReplaceAll`/`Replace`/etc., substitutes the pattern bindings into the already-evaluated `rhs`. Contrast `RuleDelayed`. Attribute: `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
