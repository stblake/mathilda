# StringPadLeft

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPadLeft["string", n]`**

Makes "string" length n, padding on the left with spaces or truncating (keeping the last n characters) as needed.

**`StringPadLeft["string", n, "padding"]`**

Pads with repeated copies of "padding".

**`StringPadLeft[{s1, s2, ...}]`**

Pads each string on the left with spaces to the length of the longest, making them all the same length.

**`StringPadLeft[{s1, s2, ...}, n, ...]`**

Pads or truncates each string to length n.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

## Algorithm

stringpad.c - StringPadLeft and StringPadRight builtins for Mathilda

```text
StringPadLeft["str", n]        - "str" of length n, padded on the left with
                                 spaces or truncated (keeping the last n chars).
StringPadLeft["str", n, "p"]   - as above, padded by repeating copies of "p".
StringPadLeft["str"]           - returns "str" unchanged (n = length of str).
StringPadLeft[{s1, s2, ...}]   - pads each string with spaces to the length of
                                 the longest, so all become the same length.
StringPadLeft[{s1, ...}, n]    - pads or truncates each string to length n.
```

StringPadLeft[{s1, ...}, n, "p"] - as above, using padding string "p".

StringPadRight is identical except padding is added on the right and truncation keeps the first n characters.

Padding is laid down as cyclic copies of the pad string read left-to-right (pad[i] = p[i % plen]) on both sides, truncated when the target width is reached; this matches the Wolfram Language documentation and all documented examples (e.g. StringPadLeft["abcde", 10, "."] -> ".....abcde").

Strings are treated as raw byte arrays (consistent with StringRepeat / StringTake / StringPartition across this subsystem); no UTF-8 codepoint decoding is performed, so lengths count bytes.

## Implementation notes

**Attributes:** `Protected`.

## See also

[StringPadRight](../../string-operations/StringPadRight/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringpad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpad.c)
