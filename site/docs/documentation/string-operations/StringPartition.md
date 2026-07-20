# StringPartition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringPartition["string", n]
    Partitions string into non-overlapping substrings of length n.
StringPartition["string", n, d]
    Generates length-n substrings with offset d (all of length n; some
    trailing or middle characters may be omitted).
StringPartition["string", UpTo[n]]
    Partitions into substrings of length up to n, allowing a shorter
    final substring so every character appears.
StringPartition[{s1, s2, ...}, spec]
    Threads over a list of strings.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringPartition["123456789123456789", 9]
Out[1]= {"123456789", "123456789"}

In[2]:= StringPartition["123456789", 2, 1]
Out[2]= {"12", "23", "34", "45", "56", "67", "78", "89"}

In[3]:= StringPartition["123456789", UpTo[2]]
Out[3]= {"12", "34", "56", "78", "9"}

In[4]:= StringPartition["ababababab", 3]
Out[4]= {"aba", "bab", "aba"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
