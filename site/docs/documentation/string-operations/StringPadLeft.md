# StringPadLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringPadLeft["string", n]
    Makes "string" length n, padding on the left with spaces or
    truncating (keeping the last n characters) as needed.
StringPadLeft["string", n, "padding"]
    Pads with repeated copies of "padding".
StringPadLeft[{s1, s2, ...}]
    Pads each string on the left with spaces to the length of the
    longest, making them all the same length.
StringPadLeft[{s1, s2, ...}, n, ...]
    Pads or truncates each string to length n.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringPadLeft["abcde", 10]
Out[1]= "     abcde"

In[2]:= StringPadRight["abcde", 10, "."]
Out[2]= "abcde....."

In[3]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}]
Out[3]= {"    a", "   ab", "  abc", " abcd", "abcde"}

In[4]:= StringPadLeft[{"a", "ab", "abc", "abcd", "abcde"}, 3]
Out[4]= {"  a", " ab", "abc", "bcd", "cde"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
