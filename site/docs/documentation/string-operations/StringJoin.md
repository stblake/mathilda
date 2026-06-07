# StringJoin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringJoin["s1", "s2", ...]
    Concatenates strings together.
    StringJoin[{"s1", "s2", ...}] flattens all lists.
    The infix form is "s1" <> "s2" <> ...
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringJoin["abcd", "ABCD", "xyz"]
Out[1]= "abcdABCDxyz"

In[2]:= "abcd" <> "ABCD" <> "xyz"
Out[2]= "abcdABCDxyz"

In[3]:= StringJoin[{{"AB", "CD"}, "XY"}]
Out[3]= "ABCDXY"

In[4]:= StringJoin[]
Out[4]= ""

In[5]:= StringJoin["a", x]
Out[5]= StringJoin["a", x]

In[6]:= StringJoin[Characters["hello"]]
Out[6]= "hello"
```

## Implementation notes

`builtin_stringjoin` gathers all leaf strings into a growable `const char**` array via the recursive helper `collect_strings`, which descends through any `List` wrappers and borrows (does not copy) each `EXPR_STRING`'s `data.string`; any non-string, non-`List` leaf aborts with `NULL` (unevaluated). It then sums the lengths, `malloc`s one buffer of `total_len + 1`, `memcpy`s each fragment in order, and returns a single `EXPR_STRING`. The zero-argument form yields `""`. Registered with `ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED`, so the evaluator flattens nested `StringJoin` (and the `<>` infix operator) before the builtin runs.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
