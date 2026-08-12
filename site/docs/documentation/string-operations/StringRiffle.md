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

## Examples (4)

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

**Attributes:** `Protected`.

## See also

[StringSplit](../../string-operations/StringSplit/), [ToString](../../expression-information/ToString/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_stringriffle.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringriffle.c)
