# StringPadRight

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPadRight["string", n]`**

Makes "string" length n, padding on the right with spaces or truncating (keeping the first n characters) as needed.

**`StringPadRight["string", n, "padding"]`**

Pads with repeated copies of "padding".

**`StringPadRight[{s1, s2, ...}]`**

Pads each string on the right with spaces to the length of the longest, making them all the same length.

**`StringPadRight[{s1, s2, ...}, n, ...]`**

Pads or truncates each string to length n.

## Examples (7)

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

### Applications (3)

Right-pad to width 5 with zeros

```mathematica
In[5]:= StringPadRight["42", 5, "0"]
Out[5]= "42000"
```

A multi-character pad cycles

```mathematica
In[6]:= StringPadRight["abc", 7, ".-"]
Out[6]= "abc.-.-"
```

Too long: truncation keeps the FIRST n

```mathematica
In[7]:= StringPadRight["abcdef", 3]
Out[7]= "abc"
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

**Algorithm.** `StringPadRight` is the mirror of `StringPadLeft` and shares the same `pad_dispatch` (invoked with the `left` flag false). It pads on the right and, when the string is longer than `n`, truncates keeping the first `n` bytes. The pad string (a single space by default) is laid down cyclically from the left, `p[i mod plen]`, so a multi-character pad such as `".-"` repeats `.-.-…`. A `List` first argument pads every element to the explicit `n`, or to the longest element in the one-argument form.

**Data structures.** One output buffer per string; a `List` result for a list input.

**Complexity / limits.** `O(n)` per string. An arity outside 1–3 emits `StringPadRight::argb`. A non-integer/negative `n`, a non-string or list-valued pad, or an empty pad when padding is required, leaves the call unevaluated. Not `Listable`; byte-length based.

**Attributes:** `Protected`.

## References

**See also:** [StringPadLeft](../../string-operations/StringPadLeft/)

- Source: [`src/strings/stringpad.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpad.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringpad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpad.c)

## Notes & additional examples

### Notes

`StringPadRight` is the mirror of `StringPadLeft`: it pads on the right and, when
truncating, keeps the first `n` bytes. The two share one implementation
(`pad_dispatch`) selected by a `left` flag.

The pad string is laid down cyclically from the left, so a multi-character pad
such as `".-"` repeats `.-.-...`. Lengths are byte counts.
