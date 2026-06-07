# Information

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Information[symbol] or ?symbol returns information on symbol.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ?Range
Out[1]= "Range[n]
	generates the list {1, 2, 3, ..., n}.
Range[n, m]
	generates the list {n, n + 1, ..., m - 1, m}.
Range[n, m, d]
	uses step d."
```

## Implementation notes

`builtin_information` (`src/core.c`) looks up the symbol's docstring with `symtab_get_docstring` and returns it as a string. If none exists it returns a string `No information available for symbol "..."` using `context_display_name` for the shortened name. (The interactive `?name` syntax routes to the same docstring store.)

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
