# StringMatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringMatchQ["string", patt]
    Gives True if the whole "string" matches patt, and False otherwise.
StringMatchQ[{s1, s2, ...}, patt]
    Gives the list of results for each of the si.

    patt may be RegularExpression["re"], a literal string, or a list of
    alternatives.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
