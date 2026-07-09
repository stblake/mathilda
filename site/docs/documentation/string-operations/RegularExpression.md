# RegularExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RegularExpression["regex"]
    Represents a class of strings given by the PCRE regular expression
    "regex", for use in StringMatchQ, StringCases, StringReplace and
    StringSplit.  It is an inert head: it evaluates to itself.

    Supported syntax includes . [c1c2] [c1-c2] [^...] p* p+ p? p{m,n},
    non-greedy *? +? ??, groups (...) and alternation |; the classes
    \d \D \s \S \w \W and [[:name:]]; the anchors ^ $ \b \B; and
    inline options (?i) (?m) (?s).  In a replacement right-hand side $n
    stands for the n-th captured group and $0 for the whole match.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringCases["adefgh12c34", RegularExpression["[a-e]+"]]
Out[1]= {"ade", "c"}

In[2]:= StringCases["a23b4222c63333d80", RegularExpression["\\d+"]]
Out[2]= {"23", "4222", "63333", "80"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
