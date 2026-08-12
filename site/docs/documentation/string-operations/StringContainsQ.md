# StringContainsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringContainsQ["string", patt]
    Gives True if any substring of "string" matches the string
    expression patt, and False otherwise.
StringContainsQ["string", {p1, p2, ...}]
    Gives True if any substring matches any of the pi.
StringContainsQ[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.
StringContainsQ[patt]
    Represents an operator form that can be applied to a string.

    Equivalent to !StringFreeQ["string", patt], and to
    StringMatchQ["string", ___ ~~ patt ~~ ___].
    Options: IgnoreCase -> True treats upper/lowercase as equivalent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringContainsQ["bcde", "b" ~~ __ ~~ "e"]
Out[1]= True

In[2]:= StringContainsQ[{"a", "b", "ab", "abcd", "bcde"}, "a"]
Out[2]= {True, False, True, True, False}

In[3]:= StringContainsQ["bac 123", RegularExpression["a.*"] ~~ DigitCharacter ..]
Out[3]= True

In[4]:= StringContainsQ["abcd", "BC", IgnoreCase -> True]
Out[4]= True

In[5]:= Select[{"abc", "xyz", "bat"}, StringContainsQ["a"]]
Out[5]= {"abc", "bat"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
