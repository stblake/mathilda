# InputForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InputForm[expr]`**

prints expr in a form suitable to be re-read by the parser, using operator syntax (a + b, not Plus\[a, b\]) and explicit string quotes.

<details>
<summary>Notes</summary>

Like FullForm, InputForm is a printer wrapper: it is consumed during output and does not appear in the printed result.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= InputForm[1/2]
Out[1]= 1/2

In[2]:= InputForm[a + b]
Out[2]= a + b

In[3]:= InputForm[{1, 2, 3}]
Out[3]= {1, 2, 3}
```

```mathematica
In[1]:= InputForm[1/2 + 3/4 I]
Out[1]= 1/2 + 3/4*I
```

```mathematica
In[1]:= InputForm[{1, 1/2, "a", x}]
Out[1]= {1, 1/2, "a", x}
```

## Implementation notes

`InputForm` is an unevaluated display wrapper: `builtin_inputform` (`src/print.c`) returns `NULL`, leaving `InputForm[expr]` intact. The printer's `print_standard` detects the `InputForm` head and renders the argument in a re-parseable form (a printer flag toggles InputForm-specific formatting); `ToString[expr, InputForm]` routes through the same standard printer.

**Attributes:** `Protected`.

## References

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_numeric_largearg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_largearg.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_print.c`](https://github.com/stblake/mathilda/blob/main/tests/test_print.c)

## Notes & additional examples

### Notes

`InputForm` prints an expression in a form the parser can read back in, unlike `FullForm` which exposes the internal tree. It is the form to use when you need to copy a result back into the REPL or store it as text.
