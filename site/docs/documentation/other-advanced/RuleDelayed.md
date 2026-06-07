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

**Attributes:** `HoldRest`, `Protected`, `SequenceHold`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
