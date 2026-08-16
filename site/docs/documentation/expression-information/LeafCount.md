# LeafCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LeafCount[expr] gives the total number of indivisible subexpressions in expr.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= LeafCount[1 + a + b^2]
Out[1]= 6

In[2]:= LeafCount[f[a, b][x, y]]
Out[2]= 5

In[3]:= LeafCount[{1/2, 1 + I}]
Out[3]= 7
```

### Applications (5)

```mathematica
In[4]:= LeafCount[x + y]
Out[4]= 3

In[5]:= LeafCount[a + b^2 + Sin[c d]]
Out[5]= 9

In[6]:= LeafCount[Expand[(1 + x)^10]]
Out[6]= 48

In[7]:= LeafCount[Integrate[1/(1 + x^4), x]]
Out[7]= 89

In[8]:= Map[LeafCount, {1, 1/2, x, f[x], {a, b, c}}]
Out[8]= {1, 3, 1, 2, 4}
```

## Implementation notes

`builtin_leafcount` (`src/core.c`) returns `leaf_count_internal`, which counts 1 per non-`EXPR_FUNCTION` (atomic) node and recurses into function arguments. By default heads are counted too; the option `Heads -> False` suppresses head counting.

- `Protected`.
- Counts the number of subexpressions in `expr` that correspond to "leaves" on the expression tree.
- By default `Heads -> True` includes the head of expressions and their parts. With `Heads -> False`, it excludes them.
- Evaluates atoms like `Rational` and `Complex` based on their structural representation as functions.

**Attributes:** `Protected`.

## References

**See also:** [Rational](../../arithmetic/Rational/), [Complex](../../arithmetic/Complex/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_sum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum.c)

## Notes & additional examples

### Notes

`LeafCount[expr]` gives the total number of indivisible subexpressions (leaves)
in `expr`, counting heads and structural atoms. It is the standard measure used
by `Simplify` and friends to decide which of two candidate forms is "smaller".
