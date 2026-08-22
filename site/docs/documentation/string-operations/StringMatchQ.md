# StringMatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringMatchQ["string", patt]`**

Gives True if the whole "string" matches patt, and False otherwise.

**`StringMatchQ[{s1, s2, ...}, patt]`**

Gives the list of results for each of the si. patt may be RegularExpression\["re"\], a literal string, or a list of alternatives.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringMatchQ["12345", RegularExpression["\\d+"]]
Out[1]= True

In[2]:= StringMatchQ[{"12", "x"}, RegularExpression["\\d+"]]
Out[2]= {True, False}
```

## Algorithm

stringmatchq.c - StringMatchQ[subject, pattern]

Returns True if the WHOLE subject string matches the pattern, else False. The pattern may be RegularExpression["re"], a literal string, or a List of

```text
alternatives (matches if any one matches).  A list of subjects threads,
```

giving a list of True/False.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringcontainsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringcontainsq.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)
