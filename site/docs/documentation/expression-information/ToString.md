# ToString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToString[expr]`**

gives the printed form of expr (as InputForm) as a String.

**`ToString[expr, form]`**

uses the specified output form.

<details>
<summary>Notes</summary>

Supported forms: InputForm (default), FullForm, TeXForm.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

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
- `TeXForm` renders the generated constant `C[k]` (from `DSolve`, `Reduce`, `Integrate`) as the subscripted `c_k`, matching Mathematica: `ToString[C[1], TeXForm]` is `"c_1"`, and `ToString[C[10], TeXForm]` is `"c_{10}"` (single-character subscripts bare, longer ones braced). The same holds for the notebook LaTeX renderer.

**Attributes:** `Protected`.

## References

**See also:** [String](../../other-advanced/String/), [InputForm](../../expression-information/InputForm/), [FullForm](../../expression-information/FullForm/), [TeXForm](../../expression-information/TeXForm/), [DSolve](../../calculus/DSolve/), [Reduce](../../solutions-of-equations/Reduce/), [Integrate](../../calculus/Integrate/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)
