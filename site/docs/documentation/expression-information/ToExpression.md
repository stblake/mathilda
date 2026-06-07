# ToExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToExpression[input]
    parses input (a String) as Mathilda input and returns the resulting
    expression after evaluation.
ToExpression[input, form]
    uses interpretation rules for the specified form. form may be
    InputForm or FullForm (both currently use the same parser).
ToExpression[input, form, h]
    wraps the head h around the parsed expression before evaluation;
    use h = Hold to obtain the unevaluated parsed form.
Returns $Failed if a syntax error is encountered.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ToExpression["1+1"]
Out[1]= 2

In[2]:= ToExpression["1+1", InputForm, Hold]
Out[2]= Hold[1 + 1]

In[3]:= ToExpression["x+"]
Out[3]= $Failed
```

## Implementation notes

`builtin_toexpression` (`src/core.c`) feeds a string argument to `parse_expression` (the Pratt parser, `src/parse.c`) and returns the parsed tree for the evaluator to reduce. An optional second argument (`InputForm`/`FullForm`/`StandardForm`) is accepted but ignored, since the parser is form-agnostic; an optional third argument is a head `h` wrapped around the result (commonly `Hold`). A parse failure returns `$Failed`; non-string input returns `NULL`. The symbol is `Listable`.

- `Protected`, `Listable`. `ToExpression[{"1+1", "2+2"}]` evaluates to `{2, 4}`.
- Returns the symbol `$Failed` if the parser cannot consume the input.
- A non-string input or an unsupported `form` leaves the call unevaluated.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
