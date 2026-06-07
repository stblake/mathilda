# StringPart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringPart["string", n]
    Gives the nth character in "string".
StringPart["string", {n1, n2, ...}]
    Gives a list of the ni-th characters.
StringPart["string", m;;n;;s]
    Gives characters m through n in steps of s.
StringPart[{s1, s2, ...}, spec]
    Gives the list of results for each si.

    Negative indices count from the end.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringPart["abcdefghijklm", 6]
Out[1]= "f"

In[2]:= StringPart["abcdefghijklm", {1, 3, 5}]
Out[2]= {"a", "c", "e"}

In[3]:= StringPart["abcdefghijklm", -4]
Out[3]= "j"

In[4]:= StringPart["abcdefghijklm", 1;;6]
Out[4]= {"a", "b", "c", "d", "e", "f"}

In[5]:= StringPart["abcdefghijklm", 1;;-1;;2]
Out[5]= {"a", "c", "e", "g", "i", "k", "m"}

In[6]:= StringPart["abcdefghijklm", -1;;1;;-2]
Out[6]= {"m", "k", "i", "g", "e", "c", "a"}

In[7]:= StringPart[{"abcd", "efgh", "ijklm"}, 1]
Out[7]= {"a", "e", "i"}

In[8]:= StringPart[{"abcd", "efgh", "ijklm"}, {1, -1}]
Out[8]= {{"a", "d"}, {"e", "h"}, {"i", "m"}}
```

## Implementation notes

`builtin_stringpart` takes `(string, spec)` and indexes by byte (1-based; negative counts from the end via `len + k + 1`). The single-index path uses the helper `stringpart_single`, which bounds-checks `k` and returns a length-1 `EXPR_STRING`. A `List` spec maps `stringpart_single` over each index into a result `List`; a `Span[m, n, s]` spec resolves `start`/`end`/`step` (honouring `All` and negative endpoints), computes the element count, and emits the stepped characters as a `List`. When the first argument is itself a `List` of strings, the builtin recurses by constructing and evaluating an inner `StringPart[si, spec]` per element. Out-of-range or non-integer indices return `NULL`. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
