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

- Each `fi` may be any function — a builtin head such as `Together` or `Cancel`,

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
