# StringRepeat

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringRepeat["string", n]
    Gives a string with "string" repeated n times.
StringRepeat["string", n, max]
    Gives up to n copies of "string", truncated to a total length of
    at most max characters.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringRepeat["a", 50]
Out[1]= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

In[2]:= StringRepeat["abc", 10]
Out[2]= "abcabcabcabcabcabcabcabcabcabc"

In[3]:= StringRepeat["ab", 10, 19]
Out[3]= "abababababababababa"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
