# StringSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringSplit["string"]
    Splits "string" at runs of whitespace.
StringSplit["string", patt]
    Splits at delimiters matching the string pattern patt.
StringSplit["string", {p1, p2, ...}]
    Splits at any of the pi.
StringSplit["string", patt -> val]
    Inserts val at the position of each delimiter.
StringSplit["string", patt, n]
    Splits into at most n substrings.
StringSplit[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.

    Empty substrings between adjacent interior delimiters are kept; those
    at the start or end are dropped unless All is given as the third
    argument. "" splits at every character. Option: IgnoreCase -> True.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/strings/regex/regex_init.c`](https://github.com/stblake/mathilda/blob/main/src/strings/regex/regex_init.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
