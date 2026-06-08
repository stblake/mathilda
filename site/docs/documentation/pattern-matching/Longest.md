# Longest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Longest[p] is a pattern object that matches the longest sequence consistent with the pattern p.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Longest` is a pattern wrapper, not a function. It is recognised structurally by the matcher in `match.c`: the sequence-matching loop peels off `Longest[p]` (alongside `Pattern`, `Shortest`, and `Optional` wrappers) before binding, setting an `is_longest` flag on the underlying pattern `p`. The flag controls the order in which the backtracking sequence matcher tries argument-count partitions for `__`/`___`/`Repeated` sub-patterns — `Longest` makes the matcher try the greediest (largest) consumption first (the default for `__`/`___`), whereas `Shortest` flips the order to try the smallest first. It does not change *which* matches are possible, only which one is found first. `Longest` appears in the matcher's list of recognised pattern heads that `eval.c` keeps unevaluated.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ReplaceList[{a, b, c}, {Longest[x__], y__} :> {{x}, {y}}]
Out[1]= {{{a, b}, {c}}, {{a}, {b, c}}}

In[2]:= MatchQ[{1, 1}, {Longest[1 ..]}]
Out[2]= True
```

### Notes

`Longest[p]` forces a sequence pattern to consume as many elements as possible. In `ReplaceList` (Out[1]) the longest partition is tried first; contrast with `Shortest`, which orders the matches the other way.
