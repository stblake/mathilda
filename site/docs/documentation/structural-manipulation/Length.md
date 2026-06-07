# Length

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Length[expr]
    gives the number of top-level elements in expr (the arity of its
    head).  Length of any atom is 0.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Length[{a, b, c}]
Out[1]= 3
```

## Implementation notes

`builtin_length` (in `src/core.c`) returns an `EXPR_INTEGER` equal to the argument's `arg_count` when it is an `EXPR_FUNCTION`, and `0` for atoms (which have no parts).

- Returns the count of top-level arguments for functions.
- Returns `0` for all atomic expressions.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
