# RuleDelayed

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
lhs :> rhs or RuleDelayed[lhs, rhs]
    represents a delayed rewrite rule: rhs is held and evaluated only
    each time the rule fires, after the pattern bindings on lhs are
    substituted into rhs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`RuleDelayed[lhs, rhs]` (`:>`) is a passive delayed-rewrite object with no builtin handler. Its `ATTR_HOLDREST | ATTR_SEQUENCEHOLD | ATTR_PROTECTED` attributes hold `rhs` unevaluated at construction; the rule engine (`is_rule` in `src/replace.c`) detects the `RuleDelayed` head and, each time the rule fires, substitutes the fresh pattern bindings into the held `rhs` and only then evaluates it (see the `delayed` flag threaded through `ReplaceRule`/`ReplaceListState`). This is the difference from `Rule`, whose RHS is evaluated once up front.

**Attributes:** `HoldRest`, `Protected`, `SequenceHold`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
