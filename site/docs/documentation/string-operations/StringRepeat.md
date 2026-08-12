# StringRepeat

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRepeat["string", n]`**

Gives a string with "string" repeated n times.

**`StringRepeat["string", n, max]`**

Gives up to n copies of "string", truncated to a total length of at most max characters.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringRepeat["a", 50]
Out[1]= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

In[2]:= StringRepeat["abc", 10]
Out[2]= "abcabcabcabcabcabcabcabcabcabc"

In[3]:= StringRepeat["ab", 10, 19]
Out[3]= "abababababababababa"
```

## Algorithm

stringrepeat.c - StringRepeat builtin for Mathilda

```text
StringRepeat["str", n]        - "str" concatenated n times.
StringRepeat["str", n, max]   - up to n copies of "str", truncated so the
                                total length is at most max (a partial final
                                copy is allowed).
```

Strings are treated as raw byte arrays (consistent with StringTake/StringDrop and StringPartition across this subsystem); no UTF-8 codepoint decoding is performed, so lengths count bytes.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringrepeat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringrepeat.c)
