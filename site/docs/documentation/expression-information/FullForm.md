# FullForm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FullForm[expr]
    prints expr as its raw internal tree (heads written before arguments
    in functional form, no operator or infix sugar).
FullForm is a wrapper recognised by Print/Out; when an input evaluates
to FullForm[expr] the wrapper is consumed by the printer and does not
appear in the output.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`FullForm` is an unevaluated display wrapper: the builtin `builtin_fullform` (`src/print.c`) returns `NULL`, leaving `FullForm[expr]` intact. Rendering is done by the printer — `print_standard` detects the `FullForm` head and calls `expr_print_fullform`, which writes the raw tree as `head[arg, ...]` with no infix/operator sugar. `ToString[expr, FullForm]` reuses the same path via `expr_to_string_fullform`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/print.c`](https://github.com/stblake/mathilda/blob/main/src/print.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FullForm[a + b]
Out[1]= Plus[a, b]

In[2]:= FullForm[1/2]
Out[2]= Rational[1, 2]

In[3]:= FullForm[x^2 + 1]
Out[3]= Plus[1, Power[x, 2]]
```

### Notes

`FullForm` reveals the raw internal tree, with every head written before its arguments and no special-cased syntax. It is the quickest way to see how surface notation like `+`, `/`, and `^` maps onto the underlying `Plus`/`Rational`/`Power` heads.
