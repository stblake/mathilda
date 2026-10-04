# StringRiffle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringRiffle[{s1, s2, ...}]`**

Joins the si into a string with spaces between them; nested lists use spaces at the lowest level and increasing numbers of newlines at higher levels. Non-string elements are converted with ToString.

**`StringRiffle[list, sep]`**

Inserts the string sep between the top-level elements.

**`StringRiffle[list, {"left", "sep", "right"}]`**

Joins with sep and wraps the result in the left/right delimiters.

**`StringRiffle[list, sep1, sep2, ...]`**

Inserts separator sep\_i (a string or {left, sep, right}) between elements at level i.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringRiffle[{"a", "b", "c", "d", "e"}]
Out[1]= "a b c d e"

In[2]:= StringRiffle[{"a", "b", "c", "d", "e"}, ", "]
Out[2]= "a, b, c, d, e"

In[3]:= StringRiffle[{"a", "b", "c", "d", "e"}, {"(", " ", ")"}]
Out[3]= "(a b c d e)"

In[4]:= StringRiffle[{{"a", 27}, {"b", 28}, {"c", 29}}, {"{", ", ", "}"}, ": "]
Out[4]= "{a: 27, b: 28, c: 29}"
```

### Applications (3)

Default scheme: a single space

```mathematica
In[5]:= StringRiffle[{"a", "b", "c"}]
Out[5]= "a b c"
```

A custom separator

```mathematica
In[6]:= StringRiffle[{"2024", "01", "02"}, "-"]
Out[6]= "2024-01-02"
```

A {left, sep, right} triple

```mathematica
In[7]:= StringRiffle[{"x", "y", "z"}, {"(", ", ", ")"}]
Out[7]= "(x, y, z)"
```

## Algorithm

stringriffle.c - StringRiffle builtin for Mathilda

StringRiffle assembles a string from a (possibly nested) list of elements by inserting separators between them - the inverse of StringSplit. Non-string leaves are converted with expr_to_string (so the integer 27 -> "27"); string leaves are used verbatim.

Forms:

```text
  StringRiffle[list]                       - default separator scheme: a
      single space at the innermost level, one extra newline per level going
      up (2-D: rows by "\n", cells by " "; 3-D: blocks by "\n\n", ...).
  StringRiffle[list, sep]                  - string sep between the top-level
      elements; deeper levels fall back to the default scheme.
  StringRiffle[list, {"l","sep","r"}]      - a 3-string list is a delimiter
      triple: wrap the join with "l"..."r", joining with "sep".
  StringRiffle[list, sep1, sep2, ...]      - sep_i (a string or a delimiter
      triple) between elements at level i (1 = outermost); deeper levels use
      the default scheme.
```

Strings are treated as raw byte arrays (consistent with the rest of the string subsystem); no UTF-8 decoding is performed.

## Implementation notes

**Algorithm.** `builtin_stringriffle` takes the first argument as the data and the rest as per-level separators (level 1 outermost). Each separator is parsed by `parse_sep` into a `{left, sep, right}` triple — a plain string becomes `{"", sep, ""}`, a 3-string list is the triple itself. `riffle_build` recurses: a leaf renders through `leaf_to_str` (strings verbatim, any other expression via `expr_to_string`/`ToString`), and each level joins its children with the explicit separator for that level or, when none is supplied, the default scheme chosen from `depth_from_bottom` — a single space at the innermost level and one extra newline per level above. It is the inverse of `StringSplit`.

**Data structures.** A `SepSpec` array of resolved triples; per level the child strings are built first, then assembled two-pass (sum lengths, then copy) into one buffer.

**Complexity / limits.** `O(total output)`. No arguments emits `StringRiffle::argm`; a first argument that is neither a list nor a string, or a malformed separator, leaves the call unevaluated. `StringRiffle` is deliberately not `Listable` so it can inspect the whole nested structure. Byte-oriented.

**Attributes:** `Protected`.

## References

**See also:** [StringSplit](../../string-operations/StringSplit/), [ToString](../../expression-information/ToString/)

- Source: [`src/strings/stringriffle.c`](https://github.com/stblake/mathilda/blob/main/src/strings/stringriffle.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringriffle.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringriffle.c)

## Notes & additional examples

### Notes

`StringRiffle` is the inverse of `StringSplit`: it joins a (possibly nested) list
with separators. A plain string separator goes between the top-level elements; a
3-string list is a `{left, sep, right}` delimiter triple.

The default scheme uses a single space at the innermost level and one extra
newline per level above it. Non-string leaves are rendered with `ToString` (so
`27` becomes `"27"`); string leaves are used verbatim.
