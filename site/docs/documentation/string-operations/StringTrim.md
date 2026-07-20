# StringTrim

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringTrim["string"]
    Trims whitespace from the beginning and end of "string".
StringTrim["string", patt]
    Trims substrings matching the string pattern patt from the beginning
    and end.
StringTrim[{s1, s2, ...}, ...]
    Gives the list of results for each of the si.

    Whitespace covers runs of spaces, tabs, and newlines. Each end is
    trimmed to a fixed point.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringTrim["   aaa bbb ccc   "]
Out[1]= "aaa bbb ccc"

In[2]:= StringTrim["++++aaa bbb ccc----", ("+" | "-") ...]
Out[2]= "aaa bbb ccc"

In[3]:= StringTrim["   aaa bbb ccc   ", RegularExpression["^ *"]]
Out[3]= "aaa bbb ccc   "

In[4]:= StringTrim["007bond007", DigitCharacter ..]
Out[4]= "bond"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
