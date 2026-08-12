# StringInsert

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringInsert["string", "snew", n]`**

Inserts "snew" so its first character is the nth character of the result.

**`StringInsert["string", "snew", -n]`**

Inserts "snew" so its last character is the nth character from the end of the result.

**`StringInsert["string", "snew", {n1, n2, ...}]`**

Inserts a copy of "snew" at each of the positions ni.

**`StringInsert[{s1, s2, ...}, "snew", spec]`**

Gives the list of results for each of the si. Positions refer to "string" before any insertion is done.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= StringInsert["abcdefghijklm", "XYZ", 4]
Out[1]= "abcXYZdefghijklm"

In[2]:= StringInsert["abcdefghijklm", "XYZ", -4]
Out[2]= "abcdefghijXYZklm"

In[3]:= StringInsert["abcdefghijklm", "XYZ", {2, 3, 7}]
Out[3]= "aXYZbXYZcdefXYZghijklm"

In[4]:= StringInsert["1234567890123456", ".", Range[4, 16, 3]]
Out[4]= "123.456.789.012.345.6"

In[5]:= StringInsert["1234567890123456", ".", Range[-16, -4, 3]]
Out[5]= "1.234.567.890.123.456"

In[6]:= StringInsert[{"abc", "de"}, "X", 2]
Out[6]= {"aXbc", "dXe"}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)
