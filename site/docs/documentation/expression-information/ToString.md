# ToString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToString[expr]
    gives the printed form of expr (as InputForm) as a String.
ToString[expr, form]
    uses the specified output form.
Supported forms: InputForm (default), FullForm, TeXForm.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ToString[x^2 + y^3]
Out[1]= "x^2 + y^3"

In[2]:= ToString[x^2 + y^3, FullForm]
Out[2]= "Plus[Power[x, 2], Power[y, 3]]"

In[3]:= ToString[x^2 + y^3, TeXForm]
Out[3]= "x^{2}+y^{3}"
```

## Implementation notes

`builtin_tostring` (`src/core.c`) renders an expression to a string. The optional second argument selects the form: `FullForm` uses `expr_to_string_fullform`; `TeXForm` wraps in `TeXForm[...]` and prints; `InputForm`/`StandardForm`/`OutputForm` (and the default) use the standard printer `expr_to_string`. All formatting is shared with the `src/print.c` printer.

- `Protected`.
- An unsupported form leaves the call unevaluated (e.g. `ToString[x, FooForm]` returns `ToString[x, FooForm]`), so a typo is visible at the call site rather than silently downgraded.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
