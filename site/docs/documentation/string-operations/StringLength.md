# StringLength

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
StringLength["string"]
    Gives the number of characters in a string.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_stringlength` checks for a single `EXPR_STRING` argument and returns `expr_new_integer((int64_t)strlen(arg->data.string))` — a byte count, not a codepoint count. Non-string arguments leave the call unevaluated (`NULL`). `ATTR_LISTABLE | ATTR_PROTECTED`, so it threads element-wise over a list of strings automatically.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/picostrings.c`](https://github.com/stblake/mathilda/blob/main/src/picostrings.c)
- Specification: [`docs/spec/builtins/string-operations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/string-operations.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= StringLength["hello"]
Out[1]= 5

In[2]:= StringLength["Mathilda"]
Out[2]= 8

In[3]:= StringLength[""]
Out[3]= 0
```

### Notes

`StringLength` returns the number of characters in a string; the empty string has length 0.
