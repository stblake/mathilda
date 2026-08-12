# StringLength

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StringLength["string"]`**

Gives the number of characters in a string.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= StringLength["tiger"]
Out[1]= 5

In[2]:= StringLength[""]
Out[2]= 0

In[3]:= StringLength["hello world"]
Out[3]= 11

In[4]:= StringLength[x]
Out[4]= StringLength[x]
```

### Applications (3)

```mathematica
In[1]:= StringLength["hello"]
Out[1]= 5

In[2]:= StringLength["Mathilda"]
Out[2]= 8

In[3]:= StringLength[""]
Out[3]= 0
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Characters of 200k chars | 4.38 s | 2.54 s | 0.419 s |
| StringSplit on space, 200k chars | 3.12 s | 4.13 s | 0.723 s |
| StringCases regex, 200k chars | 2.36 s | 3.69 s | 3.34 s |
| StringReplace regex, 200k chars | 1.97 s | 4.74 s | 3.35 s |
| StringReplace literal, 200k chars | 0.332 s | 1.17 s | 0.195 s |
| StringCount substring, 200k chars | 0.224 s | 0.364 s | 0.102 s |

## Implementation notes

`builtin_stringlength` checks for a single `EXPR_STRING` argument and returns `expr_new_integer((int64_t)strlen(arg->data.string))` — a byte count, not a codepoint count. Non-string arguments leave the call unevaluated (`NULL`). `ATTR_LISTABLE | ATTR_PROTECTED`, so it threads element-wise over a list of strings automatically.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_integer_string.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_string.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)
- Tests: [`tests/test_nsolve_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsolve_stress.c)

## Notes & additional examples

### Notes

`StringLength` returns the number of characters in a string; the empty string has length 0.
