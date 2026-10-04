# Characters

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Characters["string"]`**

Gives a list of the characters in a string. Each character is given as a length-1 string.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (3)

A string becomes a list of length-1 strings

```mathematica
In[6]:= Characters["cat"]
Out[6]= {"c", "a", "t"}
```

Characters and StringJoin are inverse

```mathematica
In[7]:= StringJoin[Characters["abcd"]]
Out[7]= "abcd"
```

One element per byte, as in StringLength

```mathematica
In[8]:= Length[Characters["hello"]]
Out[8]= 5
```

## Implementation notes

`builtin_characters` requires a single `EXPR_STRING` argument (else returns `NULL`). It allocates an `Expr*` array of `strlen(str)` entries and, for each byte, builds a length-1 `EXPR_STRING` via `expr_new_string` from a two-char `buf`, wrapping the lot in a `List`. It is byte-oriented (one element per `char`), not Unicode-codepoint aware. Carries `ATTR_LISTABLE | ATTR_PROTECTED`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/)

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_strings.c`](https://github.com/stblake/mathilda/blob/main/tests/test_strings.c)

## Notes & additional examples

### Notes

`Characters` is byte-oriented — one element per `char`, with no UTF-8 codepoint
decoding — so its result length always equals `StringLength`. It is `Listable`,
so a list of strings gives a list of character lists, and the empty string gives
`{}`.
