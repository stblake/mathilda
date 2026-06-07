# ConditionalExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ConditionalExpression[expr, cond]
    is a symbolic construct that represents expr when cond is True.
ConditionalExpression[expr, True] evaluates to expr.
ConditionalExpression[expr, False] evaluates to Undefined.
Nested forms collapse: ConditionalExpression[ConditionalExpression[e, c1], c2] evaluates to ConditionalExpression[e, c1 && c2].
ConditionalExpression has attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ConditionalExpression[a, True]
Out[1]= a

In[2]:= ConditionalExpression[a, False]
Out[2]= Undefined

In[3]:= ConditionalExpression[x^2, x > 0]
Out[3]= ConditionalExpression[x^2, x > 0]

In[4]:= ConditionalExpression[ConditionalExpression[e, c1], c2]
Out[4]= ConditionalExpression[e, c1 && c2]
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
