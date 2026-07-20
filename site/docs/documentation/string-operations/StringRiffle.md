# StringRiffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringRiffle[{s1, s2, ...}]
    Joins the si into a string with spaces between them; nested lists
    use spaces at the lowest level and increasing numbers of newlines
    at higher levels. Non-string elements are converted with ToString.
StringRiffle[list, sep]
    Inserts the string sep between the top-level elements.
StringRiffle[list, {"left", "sep", "right"}]
    Joins with sep and wraps the result in the left/right delimiters.
StringRiffle[list, sep1, sep2, ...]
    Inserts separator sep_i (a string or {left, sep, right}) between
    elements at level i.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringRiffle[{"a", "b", "c", "d", "e"}]
Out[1]= "a b c d e"

In[2]:= StringRiffle[{"a", "b", "c", "d", "e"}, ", "]
Out[2]= "a, b, c, d, e"

In[3]:= StringRiffle[{"a", "b", "c", "d", "e"}, {"(", " ", ")"}]
Out[3]= "(a b c d e)"

In[4]:= StringRiffle[{{"a", 27}, {"b", 28}, {"c", 29}}, {"{", ", ", "}"}, ": "]
Out[4]= "{a: 27, b: 28, c: 29}"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
