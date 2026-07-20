# StringExtract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringExtract["string", n]
    Extracts the nth whitespace-delimited block of "string" (-n counts
    from the end).
StringExtract["string", spec]
    Selects blocks with spec: n, -n, {n1, n2, ...}, n1;;n2, or All.
StringExtract["string", sep -> pos]
    Delimits blocks with the string pattern sep. sep -> All equals
    StringSplit["string", sep].
StringExtract["string", pos1, pos2, ...]
    Extracts across levels: whitespace at the lowest level, then "\n",
    then "\n\n", and so on for higher levels.
StringExtract["string", sep1 -> pos1, sep2 -> pos2, ...]
    Uses sepi as the separator for successive levels.

    Absent blocks yield Missing["PartAbsent", pos]. A list of strings
    threads.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringExtract["a bbb  cccc aa   d", 2]
Out[1]= "bbb"

In[2]:= StringExtract["a bbb  cccc aa   d", 2 ;; 4]
Out[2]= {"bbb", "cccc", "aa"}

In[3]:= StringExtract["a--bbb--ccc--dddd", "--" -> 3]
Out[3]= "ccc"

In[4]:= StringExtract["a 1\nb 2\nc 3 x", All, 3]
Out[4]= {"2nc"}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
