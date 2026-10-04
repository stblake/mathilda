# Information

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Information[symbol] or ?symbol returns information on symbol.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= StringQ[Information[Sin]]
Out[1]= True

In[2]:= Head[Information[Plus]]
Out[2]= String
```

### Applications (2)

Information hands back the docstring as a string

```mathematica
In[3]:= StringQ[Information[Sin]]
Out[3]= True
```

Its result is a String

```mathematica
In[4]:= Head[Information[Plus]]
Out[4]= String
```

## Implementation notes

`builtin_information` (`src/core.c`) looks up the symbol's docstring with `symtab_get_docstring` and returns it as a string. If none exists it returns a string `No information available for symbol "..."` using `context_display_name` for the shortened name. (The interactive `?name` syntax routes to the same docstring store.)

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_eigen.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eigen.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

### Notes

`Information[sym]` returns the symbol's docstring — the same text the interactive `?sym`
shortcut prints — as a string, or a `No information available` string when the symbol has
none. Every builtin registers its docstring via `symtab_set_docstring`, which is the store
both `Information` and `?sym` read.

Because the result is the usage text itself (long, multi-line, and specific to the symbol
asked about), the examples above probe it structurally; in the REPL you would simply type
`?Sin`.
