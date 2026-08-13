# StringPartition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringPartition["string", n]`**

Partitions string into non-overlapping substrings of length n.

**`StringPartition["string", n, d]`**

Generates length-n substrings with offset d (all of length n; some trailing or middle characters may be omitted).

**`StringPartition["string", UpTo[n]]`**

Partitions into substrings of length up to n, allowing a shorter final substring so every character appears.

**`StringPartition[{s1, s2, ...}, spec]`**

Threads over a list of strings.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

## Algorithm

stringpartition.c - StringPartition builtin for Mathilda

```text
StringPartition["string", n]        - non-overlapping length-n blocks
StringPartition["string", n, d]     - length-n blocks starting every d chars
StringPartition["string", UpTo[n]]  - length-<=n blocks; final may be shorter
```

StringPartition[{s1, s2, ...}, spec] - threads over a list of strings

Strings are treated as raw byte arrays (consistent with StringTake/StringDrop across this subsystem); no UTF-8 codepoint decoding is performed.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/), [StringDrop](../../string-operations/StringDrop/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringpartition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpartition.c)
