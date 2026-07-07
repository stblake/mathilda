# StringReverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringReverse["string"]
    Reverses the order of the characters in "string".
StringReverse[{s1, s2, ...}]
    Gives the list of results for each of the si.

    StringReverse is Listable, so it threads automatically over lists.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringReverse["abcdef"]
Out[1]= "fedcba"

In[2]:= StringReverse[{"cat", "dog", "fish", "coelenterate"}]
Out[2]= {"tac", "god", "hsif", "etaretneleoc"}

In[3]:= StringReverse[""]
Out[3]= ""

In[4]:= StringReverse[x]
Out[4]= StringReverse[x]

In[5]:= StringReverse[]
Out[5]= StringReverse[]
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
