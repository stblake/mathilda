# StringEndsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringEndsQ["string", patt]
    Gives True if a suffix of "string" matches the string expression
    patt, and False otherwise.
StringEndsQ["string", {p1, p2, ...}]
    Gives True if a suffix matches any of the pi.
StringEndsQ[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.
StringEndsQ[patt]
    Represents an operator form that can be applied to a string.

    Equivalent to StringContainsQ["string", patt ~~ EndOfString].
    Options: IgnoreCase -> True treats upper/lowercase as equivalent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringEndsQ["commit", "it"]
Out[1]= True

In[2]:= StringEndsQ["commit", "mi"]
Out[2]= False

In[3]:= StringEndsQ[{"apple", "banana"}, "a"]
Out[3]= {False, True}

In[4]:= StringEndsQ["a123", DigitCharacter ..]
Out[4]= True
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
