# StringFreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringFreeQ["string", patt]`**

Gives True if no substring of "string" matches the string expression patt, and False otherwise.

**`StringFreeQ["string", {p1, p2, ...}]`**

Gives True if no substring matches any of the pi.

**`StringFreeQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si.

**`StringFreeQ[patt]`**

Represents an operator form that can be applied to a string. Equivalent to !StringContainsQ\["string", patt\]. Options: IgnoreCase -\> True treats upper/lowercase as equivalent.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= StringFreeQ["abcd", "a"]
Out[1]= False

In[2]:= StringFreeQ["abcade", x_ ~~ x_]
Out[2]= True

In[3]:= StringFreeQ[{"ability", "listable", "argument"}, "a" ~~ __ ~~ "t" ~~ ___]
Out[3]= {False, True, False}
```

### Options (1)

```mathematica
In[4]:= StringFreeQ["ac", IgnoreCase -> True]["BACCD"]
Out[4]= False
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [StringContainsQ](../../string-operations/StringContainsQ/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
