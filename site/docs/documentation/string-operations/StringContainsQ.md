# StringContainsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringContainsQ["string", patt]`**

Gives True if any substring of "string" matches the string expression patt, and False otherwise.

**`StringContainsQ["string", {p1, p2, ...}]`**

Gives True if any substring matches any of the pi.

**`StringContainsQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringContainsQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to !StringFreeQ\["string", patt\], and to StringMatchQ\["string", \_\_\_ ~~ patt ~~ \_\_\_\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringContainsQ["bcde", "b" ~~ __ ~~ "e"]
Out[1]= True

In[2]:= StringContainsQ[{"a", "b", "ab", "abcd", "bcde"}, "a"]
Out[2]= {True, False, True, True, False}

In[3]:= StringContainsQ["bac 123", RegularExpression["a.*"] ~~ DigitCharacter ..]
Out[3]= True

In[4]:= Select[{"abc", "xyz", "bat"}, StringContainsQ["a"]]
Out[4]= {"abc", "bat"}
```

### Options (1)

```mathematica
In[5]:= StringContainsQ["abcd", "BC", IgnoreCase -> True]
Out[5]= True
```

### Worked examples (1)

```mathematica
In[6]:= Options[StringContainsQ]
Out[6]= {IgnoreCase -> False}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [StringMatchQ](../../string-operations/StringMatchQ/), [SetOptions](../../assignment-and-rules/SetOptions/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
