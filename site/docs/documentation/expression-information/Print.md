# Print

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Print[expr1, expr2, ...]`**

prints each argument to stdout, concatenated without separator and followed by a newline, and returns Null.  Arguments are formatted in the default output form (matching the REPL's Out display).

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Print["Result: ", x + y] Result: x + y
Out[1]= Optional[Null Result, x + y]

In[2]:= Print[x + y // FullForm] Plus[x, y]
Out[2]= Null (x + y)
```

### Applications (3)

```mathematica
In[3]:= Print["Hello, Mathilda!"]

In[4]:= Print[2 + 3]

In[5]:= Print["x = ", 2^10]
```

## Implementation notes

`builtin_print` (`src/print.c`) calls `print_standard` on each argument in turn (no separators), emits a trailing newline, and returns the symbol `Null`.

**Attributes:** `Protected`.

## References

**See also:** [FullForm](../../expression-information/FullForm/), [InputForm](../../expression-information/InputForm/)

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)
- Tests: [`tests/test_time_constrained.c`](https://github.com/stblake/mathilda/blob/main/tests/test_time_constrained.c)

## Notes & additional examples

### Notes

`Print` writes its evaluated arguments to stdout, concatenated with no separator, then returns `Null`. String arguments are printed with their surrounding quotes.
