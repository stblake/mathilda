# Print

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Print[expr1, expr2, ...]
    prints each argument to stdout, concatenated without separator and
    followed by a newline, and returns Null.  Arguments are formatted in
    the default output form (matching the REPL's Out display).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_print` (`src/print.c`) calls `print_standard` on each argument in turn (no separators), emits a trailing newline, and returns the symbol `Null`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Print["Hello, Mathilda!"]
"Hello, Mathilda!"
Out[1]= Null

In[2]:= Print[2 + 3]
5
Out[2]= Null

In[3]:= Print["x = ", 2^10]
"x = "1024
Out[3]= Null
```

### Notes

`Print` writes its evaluated arguments to stdout, concatenated with no separator, then returns `Null`. String arguments are printed with their surrounding quotes.
