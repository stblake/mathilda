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

**Algorithm.** `builtin_conditional_expression` takes `ConditionalExpression[expr, cond]` (two args; `ATTR_PROTECTED`, no Hold, so both are pre-evaluated). When `cond` is the interned `True` it yields `expr` (stealing the slot via `args[0] = NULL` so the evaluator's free of `res` doesn't double-free); when `cond` is `False` it yields the symbol `Undefined`. Nested forms `ConditionalExpression[ConditionalExpression[e, c1], c2]` are flattened to `ConditionalExpression[e, And[c1, c2]]`, with the combined `And` run through `evaluate` so contradictory or redundant conditions collapse. Any other condition leaves the call unevaluated (`NULL`).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ConditionalExpression[x, True]
Out[1]= x
```

```mathematica
In[1]:= ConditionalExpression[1/x, x != 0] /. x -> 0
Out[1]= Undefined
```

```mathematica
In[1]:= ConditionalExpression[Sqrt[x^2], x > 0]
Out[1]= ConditionalExpression[Sqrt[x^2], x > 0]
```

```mathematica
In[1]:= ConditionalExpression[ConditionalExpression[e, a > 0], b > 0]
Out[1]= ConditionalExpression[e, a > 0 && b > 0]
```

### Notes

`ConditionalExpression[expr, cond]` carries a value together with the assumption under which it is valid. It evaluates to `expr` when the condition is `True` and to `Undefined` when it is `False`; nested forms merge their conditions with `And`.
