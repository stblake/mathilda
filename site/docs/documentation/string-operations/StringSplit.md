# StringSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringSplit["string"]`**

Splits "string" at runs of whitespace.

**`StringSplit["string", patt]`**

Splits at delimiters matching the string pattern patt.

**`StringSplit["string", {p1, p2, ...}]`**

Splits at any of the pi.

**`StringSplit["string", patt -> val]`**

Inserts val at the position of each delimiter.

**`StringSplit["string", patt, n]`**

Splits into at most n substrings.

**`StringSplit[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. Empty substrings between adjacent interior delimiters are kept; those at the start or end are dropped unless All is given as the third argument. "" splits at every character. Option: IgnoreCase -\> True.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringSplit["a bbb  cccc aa   d"]
Out[1]= {"a", "bbb", "cccc", "aa", "d"}

In[2]:= StringSplit["a-b:c-d:e-f-g", {":", "-"}]
Out[2]= {"a", "b", "c", "d", "e", "f", "g"}

In[3]:= StringSplit["a b::c d::e f g", "::" -> "--"]
Out[3]= {"a b", "--", "c d", "--", "e f g"}

In[4]:= StringSplit["This is a sentence, which goes on.", Except[WordCharacter] ..]
Out[4]= {"This", "is", "a", "sentence", "which", "goes", "on"}
```

## Algorithm

stringsplit.c - StringSplit[...], the full Wolfram surface.

```text
  StringSplit[s]                 split at runs of whitespace
  StringSplit[s, patt]           split at delimiters matching patt
  StringSplit[s, {p1, p2, ...}]  split at any of the pi
  StringSplit[s, patt -> val]    insert val at each delimiter
  StringSplit[s, patt, n]        at most n substrings
  StringSplit[s, patt, All]      keep leading/trailing empty substrings
  StringSplit[{s1, ...}, patt]   thread over a list of subjects
  IgnoreCase -> True             case-insensitive delimiters
```

The delimiter pattern is translated to PCRE by the shared string-pattern engine (regex_common.c / string_pattern.c), so it accepts literal strings, RegularExpression["re"], the character-class heads (Whitespace, ...), StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on.

Zero-length substrings between two adjacent interior delimiters are kept; empty substrings at the very beginning or end are dropped unless All is given. The empty-string delimiter "" splits at every character.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [StringReplace](../../string-operations/StringReplace/)

- Source: [`src/strings/regex/regex_init.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regex_init.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
