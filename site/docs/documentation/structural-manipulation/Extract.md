# Extract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Extract[expr, pos]
    extracts the part of expr at the position specified by pos.
Extract[expr, {pos1, pos2, ...}]
    extracts a list of parts of expr.
Extract[expr, pos, h]
    extracts parts of expr, wrapping each of them with head h before evaluation.
Extract[pos]
    represents an operator form of Extract that can be applied to an expression.

The position pos has the same form as a Position result: a list of indices {i1, i2, ...} that descends i1 into expr, then i2, etc. A scalar index n is treated as the path {n}, so Extract[expr, n] is equivalent to Extract[expr, {n}]; in particular Extract[expr, 0] gives Head[expr].
Indices are 1-based and may be negative; index 0 selects the head. Extract is treated as atomic on Integer, Real, String, Symbol, Rational[n, d], and Complex[re, im].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_extract` (in `src/part.c`) reads `Extract[expr, pos]` (optionally with a held wrapper head as a third argument). If `pos` is a list-of-positions (a `List` whose elements are themselves `List`s), it extracts each position via the helper `extract_single` and returns the results in a `List`; otherwise it treats `pos` as a single position path. The 1-arg operator form returns a `Function[Extract[#, pos]]` closure.

- Position specifications have the same form as those returned by `Position`.
- `Extract[expr, {i, j, ...}]` is equivalent to `Part[expr, i, j, ...]`.
- `pos` can be of the more general form `{part1, part2, ...}` where `parti` are `Part` specifications such as an integer `i`, `All` or `Span`.
- You can use `Extract[expr, ..., Hold]` to extract parts without evaluation.

**Attributes:** `NHoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Extract[{a,b,c},2]
Out[1]= b

In[2]:= Extract[{{a,b},{c,d}},{2,1}]
Out[2]= c

In[3]:= Extract[{{1,2},{3,4}},{{1,1},{2,2}}]
Out[3]= {1, 4}
```

### Notes

A single integer extracts one element; a position list like `{2,1}` extracts a nested part. A list of positions returns the list of extracted values.
