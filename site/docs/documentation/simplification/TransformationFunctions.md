# TransformationFunctions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TransformationFunctions
    is an option for Simplify giving the list of functions to apply to try to transform parts of an expression.
TransformationFunctions -> Automatic uses the built-in collection of transformation functions.
TransformationFunctions -> {f1, f2, ...} uses only the functions fi.
TransformationFunctions -> {Automatic, f1, ...} uses the built-in transformation functions together with the fi.
Each function is applied to the whole expression and to its subexpressions; the lowest-complexity result (per ComplexityFunction) is kept.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1), TransformationFunctions -> {Cancel}]
Out[1]= 1 + x

In[2]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {}]
Out[2]= Cos[x]^2 + Sin[x]^2
```

## Implementation notes

**Algorithm.** `TransformationFunctions` is an option symbol for Simplify (it has
no builtin handler of its own). `builtin_simplify` detects
`Rule[TransformationFunctions, spec]` among its arguments and resolves it into a
`(use_builtin, user_funcs[])` pair: `Automatic` (the default) keeps the built-in
transform pipeline only; `{f1, ...}` suppresses the built-ins and uses only the
`fi`; `{Automatic, f1, ...}` runs the built-in pipeline *and* the `fi`; a bare
`f` is treated as the single-function list `{f}`. The `fi` are borrowed pointers
into the option expression (the evaluator keeps it alive across the call). After
the built-in search produces `best` (or, when built-ins are suppressed, `best =
expr_copy(input)`), `simp_apply_transformations` applies each user function to
the current best and keeps the lowest-complexity result by `score_with_func`.

**Data structures.** No state beyond the parsed option; `user_funcs` is a small
heap array of borrowed `Expr*` head expressions, freed after the search.

- Each `fi` may be any function — a builtin head such as `Together` or `Cancel`,

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp_builtins.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_builtins.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Worked examples

By default `Simplify` applies its built-in transformations, collapsing the
Pythagorean identity:

```mathematica
In[1]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {Automatic}]
Out[1]= 1
```

Passing an empty list disables every transformation, so the same expression is
left untouched — a direct way to see which step the built-in collection was
responsible for:

```mathematica
In[1]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {}]
Out[1]= Cos[x]^2 + Sin[x]^2
```

A user-supplied transformation can stand in for the built-ins entirely: here a
single rewrite rule recovers the simplification without `Automatic`:

```mathematica
In[1]:= f = Function[e, e /. Sin[a_]^2 + Cos[a_]^2 -> 1];
In[2]:= Simplify[Cos[x]^2 + Sin[x]^2, TransformationFunctions -> {f}]
Out[2]= 1
```

Built-in transformers may also be named explicitly and combined with `Automatic`,
letting `TrigToExp` participate in the search:

```mathematica
In[1]:= Simplify[1 + Tan[x]^2, TransformationFunctions -> {Automatic, TrigToExp}]
Out[1]= Sec[x]^2
```
