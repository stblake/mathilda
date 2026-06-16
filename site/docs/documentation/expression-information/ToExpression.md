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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ToExpression["1 + 2*3"]
Out[1]= 7
```

```mathematica
In[1]:= ToExpression["D[ArcTan[x], x]"]
Out[1]= 1/(1 + x^2)
```

```mathematica
In[1]:= ToExpression["Series[Exp[x], {x, 0, 4}]"]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + O[x]^5
```

```mathematica
In[1]:= ToExpression["x^2 + 1", InputForm, Hold]
Out[1]= Hold[x^2 + 1]
```

```mathematica
In[1]:= ToExpression["bad syntax ]["]
Out[1]= $Failed
```

### Notes

`ToExpression[input]` parses a string as Mathilda input and evaluates the result
— so the full symbolic engine is reachable from text, including `D`, `Series`,
`Solve`, and friends. The three-argument form wraps a head around the parsed tree
before evaluation; `ToExpression[s, InputForm, Hold]` is the standard way to get
the *unevaluated* parsed form. A syntax error yields `$Failed` rather than an
abort, making `ToExpression` safe to use on untrusted input or in `Map`.
