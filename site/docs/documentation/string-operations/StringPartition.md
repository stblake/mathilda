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

## Examples (7)

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

### Applications (3)

Non-overlapping blocks; a short tail is dropped

```mathematica
In[5]:= StringPartition["abcdefgh", 3]
Out[5]= {"abc", "def"}
```

Length 2, offset 1: overlapping blocks

```mathematica
In[6]:= StringPartition["abcdef", 2, 1]
Out[6]= {"ab", "bc", "cd", "de", "ef"}
```

Allow a shorter final block

```mathematica
In[7]:= StringPartition["abcdefg", UpTo[3]]
Out[7]= {"abc", "def", "g"}
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

**Algorithm.** `builtin_stringpartition` accepts 2–3 arguments (else `StringPartition::argt`). The block length `n` comes from an integer (every block exactly `n`) or from `UpTo[n]` (a short final block allowed); the offset `d` defaults to `n`. It walks `start = 0, d, 2d, …`, emitting `subj[start, start+n)` while the full block fits, or the short tail `subj[start, len)` when `UpTo` permits, and stops otherwise. With `d > n` characters are skipped; with `d < n` the blocks overlap. A `List` first argument threads via `evaluate`.

**Data structures.** An `Expr*` block array bounded by `len/d + 2`; each block is a fresh `EXPR_STRING` (`sp_block` `memcpy`s the byte range).

**Complexity / limits.** `O(output)`. Both `n` and `d` must be positive integers; a non-integer or non-positive value, or a non-string subject, leaves the call unevaluated. Byte-indexed, consistent with `StringTake`/`StringDrop`.

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/), [StringDrop](../../string-operations/StringDrop/)

- Source: [`src/strings/stringpartition.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringpartition.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringpartition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringpartition.c)

## Notes & additional examples

### Notes

`StringPartition[s, n, d]` gives length-`n` substrings starting every `d`
characters (`d` defaults to `n`). With `d > n` characters are skipped; with
`d < n` the blocks overlap. By default a trailing run too short to fill a block
is dropped.

`UpTo[n]` keeps that final short block so every character appears. Both `n` and
`d` must be positive integers; indexing is byte-based. A list of strings threads.
