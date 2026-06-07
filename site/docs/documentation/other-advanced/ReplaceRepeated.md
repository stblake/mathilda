# ReplaceRepeated

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
expr //. rules or ReplaceRepeated[expr, rules]
    repeatedly applies ReplaceAll[expr, rules] until the result stops
    changing, then returns the fixed point.  Useful for chained rewrite
    systems; subject to the same recursion-limit guard as evaluator
    fixed-point iteration.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
