# StringSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringSplit["string", patt]
    Splits "string" into the substrings between non-overlapping matches
    of the delimiter patt.  Empty pieces are dropped.
StringSplit[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.

    patt may be RegularExpression["re"], a literal string, or a list of
    alternative delimiters.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringSplit["1.23, 4.56  7.89", RegularExpression["(\\s|,)+"]]
Out[1]= {"1.23", "4.56", "7.89"}

In[2]:= StringSplit["a,b,c", ","]
Out[2]= {"a", "b", "c"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
