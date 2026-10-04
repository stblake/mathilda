# StringRepeat

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRepeat["string", n]`**

Gives a string with "string" repeated n times.

**`StringRepeat["string", n, max]`**

Gives up to n copies of "string", truncated to a total length of at most max characters.

## Examples (6)

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

### Applications (3)

Four copies

```mathematica
In[4]:= StringRepeat["ab", 4]
Out[4]= "abababab"
```

A rule of ten dashes

```mathematica
In[5]:= StringRepeat["-", 10]
Out[5]= "----------"
```

Capped at 7 bytes, so the last copy is partial

```mathematica
In[6]:= StringRepeat["abc", 5, 7]
Out[6]= "abcabca"
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

**Algorithm.** `builtin_stringrepeat` accepts 2–3 arguments (else `StringRepeat::argt`). It computes the output length as `n · len`, optionally capped at `max`, then fills the buffer cyclically (`buf[i] = str[i mod len]`) so a truncated final copy falls out naturally. The `n · len` product is guarded against `size_t` overflow, which is tolerated only when a `max` cap keeps the result finite.

**Data structures.** One output buffer of the computed length.

**Complexity / limits.** `O(output)`. `n == 0`, an empty base string, or `max == 0` gives `""`; a non-string base, or a non-integer/negative count or `max`, leaves the call unevaluated. Byte-length based, consistent with the rest of the string subsystem.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringrepeat.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringrepeat.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringrepeat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringrepeat.c)

## Notes & additional examples

### Notes

`StringRepeat[s, n]` concatenates `n` copies of `s`. The optional third argument
caps the total length at `max` bytes, keeping a partial final copy; the buffer is
filled cyclically (`buf[i] = s[i mod len]`) so truncation falls out naturally.

`n == 0`, an empty base string, or `max == 0` gives `""`. The `n * len` product is
guarded against overflow, which is tolerated only when a `max` cap keeps the
result finite.
