# StringTrim

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringTrim["string"]`**

Trims whitespace from the beginning and end of "string".

**`StringTrim["string", patt]`**

Trims substrings matching the string pattern patt from the beginning and end.

**`StringTrim[{s1, s2, ...}, ...]`**

Gives the list of results for each of the si. Whitespace covers runs of spaces, tabs, and newlines. Each end is trimmed to a fixed point.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringTrim["   aaa bbb ccc   "]
Out[1]= "aaa bbb ccc"

In[2]:= StringTrim["++++aaa bbb ccc----", ("+" | "-") ...]
Out[2]= "aaa bbb ccc"

In[3]:= StringTrim["   aaa bbb ccc   ", RegularExpression["^ *"]]
Out[3]= "aaa bbb ccc   "

In[4]:= StringTrim["007bond007", DigitCharacter ..]
Out[4]= "bond"
```

## Algorithm

stringtrim.c - StringTrim[...], trims matching substrings from both ends.

```text
  StringTrim["string"]           trims whitespace runs from start and end
  StringTrim["string", patt]     trims substrings matching patt from both ends
  StringTrim[{s1, s2, ...}, ...] threads over a list of subject strings
```

The trim pattern is translated to PCRE by the shared string-pattern engine (string_pattern.c: wl_pattern_to_regex), so it accepts literal strings, RegularExpression["re"], the character-class heads (Whitespace, ...), StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on. The default pattern is Whitespace (a run of spaces / tabs / newlines).

Anchoring: rather than the whole-string \A...\z wrap the other regex builtins use, StringTrim needs the pattern anchored at just one end.

```text
  - The front is stripped by matching `(?:src)` starting at the current
    position and accepting only a match that begins exactly there.  PCRE's
    `\A`/`^` anchor to absolute offset 0, so a start-relative anchor is done
    by the ov[0] == start check, not by wrapping.
  - The back is stripped by matching `(?:src)\z` while passing a truncated
    subject length (the current end) to regex_match, so `\z` refers to the
    current end rather than the absolute string end.
```

Each end is stripped repeatedly to a fixed point (so StringTrim["xxabcxx", "x"] -> "abc"); a zero-width match ends the loop.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
