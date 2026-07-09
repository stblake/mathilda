# StringReplace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringReplace["string", patt -> rep]
    Replaces each non-overlapping match of patt in "string" by rep,
    with $n replaced by the n-th captured group and $0 by the whole
    match.
StringReplace["string", {patt1 -> rep1, patt2 -> rep2, ...}]
    Applies a list of replacement rules; at each position the leftmost
    match wins, ties broken by rule order.
StringReplace[{s1, s2, ...}, rules]
    Gives the list of results for each of the si.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringReplace["a13b12c1da32efg", RegularExpression["(\\d+)"] -> "[$1]"]
Out[1]= "a[13]b[12]c[1]da[32]efg"

In[2]:= StringReplace["123 45 6 789", RegularExpression["\\b"] :> "X"]
Out[2]= "X123X X45X X6X X789X"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
