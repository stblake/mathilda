# LeafCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LeafCount[expr] gives the total number of indivisible subexpressions in expr.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LeafCount[1 + a + b^2]
Out[1]= 6

In[2]:= LeafCount[f[a, b][x, y]]
Out[2]= 5

In[3]:= LeafCount[{1/2, 1 + I}]
Out[3]= 7
```

## Implementation notes

`builtin_leafcount` (`src/core.c`) returns `leaf_count_internal`, which counts 1 per non-`EXPR_FUNCTION` (atomic) node and recurses into function arguments. By default heads are counted too; the option `Heads -> False` suppresses head counting.

- `Protected`.
- Counts the number of subexpressions in `expr` that correspond to "leaves" on the expression tree.
- By default `Heads -> True` includes the head of expressions and their parts. With `Heads -> False`, it excludes them.
- Evaluates atoms like `Rational` and `Complex` based on their structural representation as functions.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= LeafCount[x + y]
Out[1]= 3
```

It counts every atomic subexpression, including operator heads:

```mathematica
In[1]:= LeafCount[a + b^2 + Sin[c d]]
Out[1]= 9
```

A handy proxy for symbolic "size" — here the blow-up of an expanded binomial:

```mathematica
In[1]:= LeafCount[Expand[(1 + x)^10]]
Out[1]= 48
```

Measuring the complexity of a computed result, e.g. an antiderivative:

```mathematica
In[1]:= LeafCount[Integrate[1/(1 + x^4), x]]
Out[1]= 89
```

Mapped over a list, it ranks expressions by structural weight:

```mathematica
In[1]:= Map[LeafCount, {1, 1/2, x, f[x], {a, b, c}}]
Out[1]= {1, 3, 1, 2, 4}
```

### Notes

`LeafCount[expr]` gives the total number of indivisible subexpressions (leaves)
in `expr`, counting heads and structural atoms. It is the standard measure used
by `Simplify` and friends to decide which of two candidate forms is "smaller".
