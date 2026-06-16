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
```

```mathematica
In[1]:= FullForm[1/2]
Out[1]= Rational[1, 2]
```

```mathematica
In[1]:= FullForm[x^2 + 1]
Out[1]= Plus[1, Power[x, 2]]
```

```mathematica
In[1]:= FullForm[a/b]
Out[1]= Times[a, Power[b, -1]]
```

```mathematica
In[1]:= FullForm[x_Integer /; x > 0]
Out[1]= Condition[Pattern[x, Blank[Integer]], Greater[x, 0]]
```

### Notes

`FullForm` reveals the raw internal tree, with every head written before its arguments and no special-cased syntax. It is the quickest way to see how surface notation like `+`, `/`, and `^` maps onto the underlying `Plus`/`Rational`/`Power` heads. Division is `Times` with a `Power[_, -1]` factor, so `a/b` is `Times[a, Power[b, -1]]`. It also exposes the desugaring of pattern syntax — `x_Integer /; x > 0` is really `Condition[Pattern[x, Blank[Integer]], Greater[x, 0]]` — which is invaluable when debugging why a rule does or does not match.
