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

## Examples (7)

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

### Applications (3)

Default: trim whitespace from both ends

```mathematica
In[5]:= StringTrim["   hello   "]
Out[5]= "hello"
```

Each end is trimmed to a fixed point

```mathematica
In[6]:= StringTrim["xxhelloxx", "x"]
Out[6]= "hello"
```

A literal pattern

```mathematica
In[7]:= StringTrim["...data...", "."]
Out[7]= "data"
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

**Algorithm.** `builtin_stringtrim` accepts 1–2 arguments (else `StringTrim::argt`). The pattern (default `Whitespace`) is translated once by `wl_pattern_to_regex`, then compiled into two one-end-anchored programs (`st_compile`): a front matcher `(?:src)` and a back matcher `(?:src)\z`. `st_scalar` strips the leading run by repeated matches that must *begin exactly* at the current start (checked by `ov[0] == start`, since `\A`/`^` would anchor to absolute 0), and the trailing run by `\z`-anchored matches over a truncated subject length. Each end is stripped to a fixed point, and a zero-width match ends the loop — so only the two ends are trimmed, never the middle.

**Data structures.** Two `RegexProgram` matchers (front, back); one output buffer for the trimmed window `[start, end)`.

**Complexity / limits.** `O(len)` plus the cost of the trim-run matches. A single-character pattern still removes a whole run (fixed-point per end). A non-string subject, an unsupported pattern, or a build without PCRE2 leaves the call unevaluated. Byte-based; the pattern accepts the full shared string-pattern vocabulary.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringtrim.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringtrim.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

`StringTrim` removes substrings matching the pattern from the *start and end*
only (never the middle), stripping each end repeatedly to a fixed point — so a
single-character pattern still removes a whole run. The default pattern is
`Whitespace`.

Anchoring is one-ended: the front is matched start-relative and the back with a
`\z` anchor over a truncated length, rather than the whole-string `\A...\z` wrap
the other regex builtins use. The pattern accepts the full shared string-pattern
vocabulary.
