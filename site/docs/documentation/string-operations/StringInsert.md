# StringInsert

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringInsert["string", "snew", n]`**

Inserts "snew" so its first character is the nth character of the result.

**`StringInsert["string", "snew", -n]`**

Inserts "snew" so its last character is the nth character from the end of the result.

**`StringInsert["string", "snew", {n1, n2, ...}]`**

Inserts a copy of "snew" at each of the positions ni.

**`StringInsert[{s1, s2, ...}, "snew", spec]`**

Gives the list of results for each of the si. Positions refer to "string" before any insertion is done.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= StringInsert["abcdefghijklm", "XYZ", 4]
Out[1]= "abcXYZdefghijklm"

In[2]:= StringInsert["abcdefghijklm", "XYZ", -4]
Out[2]= "abcdefghijXYZklm"

In[3]:= StringInsert["abcdefghijklm", "XYZ", {2, 3, 7}]
Out[3]= "aXYZbXYZcdefXYZghijklm"

In[4]:= StringInsert["1234567890123456", ".", Range[4, 16, 3]]
Out[4]= "123.456.789.012.345.6"

In[5]:= StringInsert["1234567890123456", ".", Range[-16, -4, 3]]
Out[5]= "1.234.567.890.123.456"

In[6]:= StringInsert[{"abc", "de"}, "X", 2]
Out[6]= {"aXbc", "dXe"}

In[7]:= StringInsert[] StringInsert::argrx: StringInsert called with 0 arguments; 3 arguments are expected.
```

### Applications (3)

Insert before the 3rd character

```mathematica
In[8]:= StringInsert["abcdef", "-", 3]
Out[8]= "ab-cdef"
```

A copy at each position

```mathematica
In[9]:= StringInsert["abcdef", "-", {2, 4}]
Out[9]= "a-bc-def"
```

A negative position counts from the end

```mathematica
In[10]:= StringInsert["abc", "X", -1]
Out[10]= "abcX"
```

## Implementation notes

**Algorithm.** `builtin_stringinsert` requires three arguments (else `StringInsert::argrx`). `si_pos_to_offset` maps each position against the *original* string — a positive `n` to offset `n - 1` (insert before character `n`), a negative `-k` to `len + n + 1` — rejecting `0` and out-of-range values. The positions are tallied into a `counts[0..len]` array, then a single pass builds the output, emitting `counts[i]` copies of the insert string before original byte `i`. Because all positions are resolved up front, every copy lands relative to the untouched input. A `List` first argument recurses per element.

**Data structures.** An `int64 counts[len + 1]` tally and one output buffer sized `len + npos·inslen`.

**Complexity / limits.** `O(len + output)`. A non-string subject or insert string, a non-integer position, or an out-of-range position leaves the call unevaluated. Positions are valid for `1 ≤ n ≤ len + 1` (and the negative mirror); byte-based.

**Attributes:** `Protected`.

## References

- Source: [`src/strings/stringinsert.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringinsert.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

`StringInsert[s, new, n]` makes the first character of `new` the n-th character
of the result (i.e. it inserts *before* original character n), so a positive `n`
maps to offset `n - 1` and a negative `-k` to `len + n + 1`.

All positions refer to the *original* string, before any insertion; they are
resolved up front into a per-position count, so every copy lands relative to the
untouched input. Valid positions run `1` to `len + 1` (and the negative mirror).
