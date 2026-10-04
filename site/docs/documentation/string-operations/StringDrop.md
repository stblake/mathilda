# StringDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringDrop["string", n]`**

Gives "string" with its first n characters dropped.

**`StringDrop["string", -n]`**

Gives "string" with its last n characters dropped.

**`StringDrop["string", {n}]`**

Gives "string" with its nth character dropped.

**`StringDrop["string", {m, n}]`**

Gives "string" with characters m through n dropped.

**`StringDrop["string", {m, n, s}]`**

Drops characters m through n in steps of s.

**`StringDrop["string", UpTo[n]]`**

Drops n characters, or as many as are available.

**`StringDrop[{s1, s2, ...}, spec]`**

Gives the list of results for each of the si. Negative indices count from the end.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (9)

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

In[9]:= StringDrop[] StringDrop::argrx: StringDrop called with 0 arguments; 2 arguments are expected.
```

### Applications (3)

Drop the first two characters

```mathematica
In[10]:= StringDrop["abcdef", 2]
Out[10]= "cdef"
```

Drop the last two

```mathematica
In[11]:= StringDrop["abcdef", -2]
Out[11]= "abcd"
```

Drop a range

```mathematica
In[12]:= StringDrop["abcdef", {2, 4}]
Out[12]= "aef"
```

## Implementation notes

**Algorithm.** `builtin_stringdrop` requires exactly two arguments (else `StringDrop::argrx`). It allocates a byte keep-mask initialised to `true`, clears the positions named by the sequence spec, and rebuilds the kept bytes in order (`stringdrop_build_kept`) — the complementary "keep" mask is how `StringDrop` is expressed as the complement of `StringTake`. Specs handled: an integer `n` (drop the first `n`, or the last `|n|` when negative), `UpTo[n]` (clamped to the length), `{n}`, `{m, n}` (a decreasing range drops nothing), and `{m, n, s}` (stepped). Negative endpoints normalise as `len + k + 1`. A `List` first argument recurses per element via `evaluate(StringDrop[si, spec])`.

**Data structures.** A `bool[len]` keep-mask and one output buffer.

**Complexity / limits.** `O(len)`. An out-of-range position, a zero step, or a non-integer spec leaves the call unevaluated; indexing is byte-based (no UTF-8 decoding), consistent with `StringTake`.

**Attributes:** `Protected`.

## References

**See also:** [StringTake](../../string-operations/StringTake/)

- Source: [`src/strings/stringdrop.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringdrop.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

`StringDrop` is the complement of `StringTake`: the same sequence specification
selects the characters to *remove*, and the survivors are concatenated in order.
Internally it marks a keep-mask over the bytes and rebuilds the kept ones.

Negative indices count from the end; indexing is byte-based. A decreasing range
`{m, n}` with `m > n` is empty and drops nothing; an out-of-range position leaves
the call unevaluated.
