# BlankSequence

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
__ or BlankSequence[] represents a sequence of one or more expressions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`BlankSequence` is a pattern object head, not an evaluator builtin: `__` parses to `BlankSequence[]` and `__h` to `BlankSequence[h]` (and `x__` to `Pattern[x, BlankSequence[...]]`). It matches a sequence of **one or more** consecutive arguments. The matcher recognises it in `src/match.c` via `is_sequence_blank`, which sets `min_len = 1` and reports the optional head constraint. Sequence matching against an argument list is handled by `match_args_internal`, which backtracks over the possible partitions of the argument run; an optional head `h` requires every element of the matched run to have head `h` (same atomic-type-to-head mapping as `Blank`). Contrast `BlankNullSequence` (`min_len = 0`).

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatchQ[f[1, 2], f[__]]
Out[1]= True

In[2]:= MatchQ[f[], f[__]]
Out[2]= False

In[3]:= FullForm[__]
Out[3]= BlankSequence[]
```

### Notes

`__` (`BlankSequence[]`) matches a sequence of *one or more* expressions, so it fails on the empty argument list (Out[2]). For zero-or-more use `BlankNullSequence` (`___`).
