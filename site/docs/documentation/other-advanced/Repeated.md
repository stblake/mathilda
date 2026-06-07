# Repeated

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
p.. or Repeated[p] is a pattern object that represents a sequence of one or more expressions, each matching p.
Repeated[p, max] represents from 1 to max expressions matching p.
Repeated[p, {min, max}] represents between min and max expressions matching p.
Repeated[p, {n}] represents exactly n expressions matching p.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Repeated[p]` (`p..`) is a pattern object handled entirely inside the matcher, not by a builtin. `is_repeated` in `src/match.c` recognises the `Repeated` head, sets the matched-length range to `[1, ∞)` by default, and parses an optional count spec: `Repeated[p, max]` gives `[1, max]`, `Repeated[p, {n}]` gives exactly `n`, `Repeated[p, {min, max}]` gives `[min, max]` (with `Infinity` allowed as an open upper bound). The argument-sequence matcher then matches a run of consecutive arguments each satisfying `p`, using the standard backtracking that explores valid run lengths within the range.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
