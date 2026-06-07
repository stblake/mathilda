# Characters

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Characters["string"]
    Gives a list of the characters in a string.
    Each character is given as a length-1 string.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Characters["ABC"]
Out[1]= {"A", "B", "C"}

In[2]:= Characters["A string."]
Out[2]= {"A", " ", "s", "t", "r", "i", "n", "g", "."}

In[3]:= Characters[""]
Out[3]= {}

In[4]:= Characters[{"ABC", "DEF", "XYZ"}]
Out[4]= {{"A", "B", "C"}, {"D", "E", "F"}, {"X", "Y", "Z"}}

In[5]:= Characters[x]
Out[5]= Characters[x]
```

## Implementation notes

`builtin_characters` requires a single `EXPR_STRING` argument (else returns `NULL`). It allocates an `Expr*` array of `strlen(str)` entries and, for each byte, builds a length-1 `EXPR_STRING` via `expr_new_string` from a two-char `buf`, wrapping the lot in a `List`. It is byte-oriented (one element per `char`), not Unicode-codepoint aware. Carries `ATTR_LISTABLE | ATTR_PROTECTED`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
