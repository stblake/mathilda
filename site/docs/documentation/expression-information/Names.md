# Names

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Names["string"] gives a sorted list of the names of symbols matching the string. Names[patt] matches a string pattern with metacharacters * (zero or more characters) and @ (one or more non-uppercase characters), or a RegularExpression["re"]. Names[{p1, p2, ...}] matches any of the patterns. Names[] lists all symbol names.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Names["List*"]
Out[1]= {"List", "ListPlot", "ListQ"}

In[2]:= Names["Ar@"]
Out[2]= {"Arg", "Array", "Arrow"}

In[3]:= Names[RegularExpression["Si."]]
Out[3]= {"Sin"}

In[4]:= MemberQ[Names["System`*"], "System`Sin"]
Out[4]= True
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
