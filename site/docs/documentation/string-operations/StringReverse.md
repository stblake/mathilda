# StringReverse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringReverse["string"]`**

Reverses the order of the characters in "string".

**`StringReverse[{s1, s2, ...}]`**

Gives the list of results for each of the si. StringReverse is Listable, so it threads automatically over lists.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringReverse["abcdef"]
Out[1]= "fedcba"

In[2]:= StringReverse[{"cat", "dog", "fish", "coelenterate"}]
Out[2]= {"tac", "god", "hsif", "etaretneleoc"}

In[3]:= StringReverse[""]
Out[3]= ""

In[4]:= StringReverse[x]
Out[4]= StringReverse[x]
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)
