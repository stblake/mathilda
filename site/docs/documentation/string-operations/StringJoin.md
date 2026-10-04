# StringJoin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringJoin["s1", "s2", ...]`**

Concatenates strings together. StringJoin\[{"s1", "s2", ...}\] flattens all lists. The infix form is "s1" \<\> "s2" \<\> ...

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= StringJoin["abcd", "ABCD", "xyz"]
Out[1]= "abcdABCDxyz"

In[2]:= "abcd" <> "ABCD" <> "xyz"
Out[2]= "abcdABCDxyz"

In[3]:= StringJoin[{{"AB", "CD"}, "XY"}]
Out[3]= "ABCDXY"

In[4]:= StringJoin[]
Out[4]= ""

In[5]:= StringJoin["a", x]
Out[5]= StringJoin["a", x]

In[6]:= StringJoin[Characters["hello"]]
Out[6]= "hello"
```

### Applications (3)

Concatenate the arguments

```mathematica
In[7]:= StringJoin["abc", "def"]
Out[7]= "abcdef"
```

The infix <> operator

```mathematica
In[8]:= "foo" <> "bar" <> "baz"
Out[8]= "foobarbaz"
```

Nested lists are flattened

```mathematica
In[9]:= StringJoin[{"a", "b", "c"}]
Out[9]= "abc"
```

## Implementation notes

`builtin_stringjoin` gathers all leaf strings into a growable `const char**` array via the recursive helper `collect_strings`, which descends through any `List` wrappers and borrows (does not copy) each `EXPR_STRING`'s `data.string`; any non-string, non-`List` leaf aborts with `NULL` (unevaluated). It then sums the lengths, `malloc`s one buffer of `total_len + 1`, `memcpy`s each fragment in order, and returns a single `EXPR_STRING`. The zero-argument form yields `""`. Registered with `ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED`, so the evaluator flattens nested `StringJoin` (and the `<>` infix operator) before the builtin runs.

**Attributes:** `Flat`, `OneIdentity`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [OneIdentity](../../expression-information/OneIdentity/)

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)
- Tests: [`tests/test_repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/tests/test_repl_hooks.c)
- Tests: [`tests/test_stringfns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stringfns.c)

## Notes & additional examples

### Notes

`StringJoin` gathers every leaf string (descending through any `List` wrappers),
sums their lengths, and copies them into one buffer. It is `Flat` and
`OneIdentity`, so the evaluator flattens nested `StringJoin` and `<>` before the
builtin runs.

The zero-argument form `StringJoin[]` is `""`. Any non-string, non-`List` leaf
leaves the call unevaluated.
