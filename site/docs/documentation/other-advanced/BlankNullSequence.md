# BlankNullSequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
___ or BlankNullSequence[] represents a sequence of zero or more expressions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`BlankNullSequence` is a pattern object head, not an evaluator builtin: `___` parses to `BlankNullSequence[]` and `___h` to `BlankNullSequence[h]` (and `x___` to `Pattern[x, BlankNullSequence[...]]`). It matches a sequence of **zero or more** consecutive arguments. The matcher recognises it in `src/match.c` via `is_sequence_blank`, which sets `min_len = 0` and reports the optional head constraint. As with `BlankSequence`, the argument-list driver `match_args_internal` backtracks over partitions of the argument run (now including the empty run); an optional head `h` requires each matched element to have head `h`. The only behavioural difference from `BlankSequence` is the zero-length minimum.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
