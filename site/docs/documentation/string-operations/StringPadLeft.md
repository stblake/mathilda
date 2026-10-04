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

Left-pad to width 5 with zeros

```mathematica
In[5]:= StringPadLeft["42", 5, "0"]
Out[5]= "00042"
```

Too long: truncation keeps the LAST n

```mathematica
In[6]:= StringPadLeft["abcdef", 3]
Out[6]= "def"
```

A list pads to the longest element

```mathematica
In[7]:= StringPadLeft[{"a", "bb", "ccc"}]
Out[7]= {"  a", " bb", "ccc"}
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

**Algorithm.** `StringPadLeft` and `StringPadRight` share `pad_dispatch`, a `left` flag selecting the side. `pad_one` either truncates — copying the last `n` bytes for the left variant, the first `n` for the right — or pads, laying the pad string down cyclically (`p[i mod plen]`) into the pad region that precedes (left) or follows (right) the string. The pad string defaults to a single space. For a `List` first argument the target length is the explicit `n`, or, in the one-argument form, the longest element — so a bare list form equalises all widths.

**Data structures.** One output buffer per string; a `List` result when the input is a list.

**Complexity / limits.** `O(n)` per string. An arity outside 1–3 emits `StringPadLeft::argb`/`StringPadRight::argb`. A non-integer or negative `n`, a non-string or list-valued pad string, or an empty pad when padding is required, leaves the call unevaluated. The pad builtins are deliberately **not** `Listable`, so the handler sees the whole list (the 1-arg form must). Byte-length based.

**Attributes:** `Protected`.

## References

**See also:** [StringPadRight](../../string-operations/StringPadRight/)

- Source: [`src/strings/stringpad.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpad.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringpad.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpad.c)

## Notes & additional examples

### Notes

`StringPadLeft` makes a string a given length, padding on the left or truncating
it. When the string is longer than `n` it keeps the last `n` bytes; the pad
string (a single space by default) is laid down cyclically, `p[i mod plen]`.

The one-argument list form pads every element to the length of the longest, so
all come out equal width. Lengths are byte counts. A list-valued pad string is
not supported and leaves the call unevaluated.
