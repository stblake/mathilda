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
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatchQ[{1, 1, 1}, {Repeated[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1, 2}, {Repeated[1]}]
Out[2]= False

In[3]:= MatchQ[{1, 2, 3}, {_Integer ..}]
Out[3]= True
```

The bound can be a range `{min, max}` or an exact count `{n}`. Combined with
`Cases`, this filters a list of lists by length — keeping only those with two or
three matching elements:

```mathematica
In[1]:= Cases[{{1}, {1, 1}, {1, 1, 1}, {1, 1, 1, 1}}, {Repeated[1, {2, 3}]}]
Out[1]= {{1, 1}, {1, 1, 1}}

In[2]:= MatchQ[{1, 2, 3, 4}, {Repeated[_Integer, {4}]}]
Out[2]= True
```

### Notes

`Repeated[p]` (postfix `p..`) matches a sequence of *one or more* expressions each satisfying `p`. The optional second argument bounds the count, e.g. `Repeated[p, {2, 4}]` for between two and four.
