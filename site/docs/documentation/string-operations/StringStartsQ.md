# StringStartsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringStartsQ["string", patt]
    Gives True if a prefix of "string" matches the string expression
    patt, and False otherwise.
StringStartsQ["string", {p1, p2, ...}]
    Gives True if a prefix matches any of the pi.
StringStartsQ[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.
StringStartsQ[patt]
    Represents an operator form that can be applied to a string.

    Equivalent to StringContainsQ["string", StartOfString ~~ patt].
    Options: IgnoreCase -> True treats upper/lowercase as equivalent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringStartsQ["commit", "co"]
Out[1]= True

In[2]:= StringStartsQ["commit", "om"]
Out[2]= False

In[3]:= StringStartsQ[{"apple", "banana"}, "a"]
Out[3]= {True, False}

In[4]:= StringStartsQ["a123", LetterCharacter ~~ DigitCharacter ..]
Out[4]= True
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
