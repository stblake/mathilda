# StringDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringDrop["string", n]
    Gives "string" with its first n characters dropped.
StringDrop["string", -n]
    Gives "string" with its last n characters dropped.
StringDrop["string", {n}]
    Gives "string" with its nth character dropped.
StringDrop["string", {m, n}]
    Gives "string" with characters m through n dropped.
StringDrop["string", {m, n, s}]
    Drops characters m through n in steps of s.
StringDrop["string", UpTo[n]]
    Drops n characters, or as many as are available.
StringDrop[{s1, s2, ...}, spec]
    Gives the list of results for each of the si.

    Negative indices count from the end.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= StringDrop["abcdefghijklm", 4]
Out[1]= "efghijklm"

In[2]:= StringDrop["abcdefghijklm", -4]
Out[2]= "abcdefghi"

In[3]:= StringDrop["abcdefghijklm", {5, 10}]
Out[3]= "abcdklm"

In[4]:= StringDrop["abcdefghijklm", {3}]
Out[4]= "abdefghijklm"

In[5]:= StringDrop["abcdefghijklm", {1, -1, 2}]
Out[5]= "bdfhjl"

In[6]:= StringDrop[{"abcdef", "xyzw", "stuv"}, -2]
Out[6]= {"abcd", "xy", "st"}

In[7]:= StringDrop["abc", UpTo[4]]
Out[7]= ""

In[8]:= StringDrop["abcdefghijklm", {5, -4}]
Out[8]= "abcdklm"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
