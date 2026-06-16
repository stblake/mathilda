# Tr

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Tr[m]
    gives the trace of the matrix m, i.e. the sum of its diagonal
    entries (for a rectangular m, sums entries m[[i, i]] up to
    Min[Dimensions[m]]).
Tr[m, f]
    combines the diagonal entries with f instead of Plus.
Tr[m, f, n]
    walks down to level n, summing the multi-index diagonal of a
    rank-n tensor.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= 15

In[2]:= Tr[{{a, b}, {c, d}}]
Out[2]= a + d

In[3]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, List]
Out[3]= {1, 5, 9}

In[4]:= Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, Plus, 1]
Out[4]= {12, 15, 18}
```

## Implementation notes

**Algorithm.** `builtin_tr` is the generalised trace `Tr[list, f, n]`. It walks the diagonal of a rank-`n` nested `List`: `extract_diagonal_element` descends `n` levels always taking element `index` at each level, yielding the `i`-th diagonal entry, and the loop collects entries until an index falls out of bounds. The collected leaves are combined with the head `f` (default `Plus`) via a single `f[d_0, d_1, ...]` call run through the evaluator (`eval_and_free`), so for a 2-D matrix this is the ordinary sum of diagonal entries.

**Data structures / limits.** A growable `Expr**` of the diagonal leaves. The depth `n` defaults to the nesting depth of the leftmost spine (`get_default_trace_depth`, min 1) and may be given explicitly as a non-negative Integer. A non-`List` argument is returned unchanged; a non-integer explicit depth leaves the call unevaluated.

- `Protected`.
- `Tr[list]` sums the diagonal elements `list[[i, i, ...]]`.
- `Tr[list, f]` applies the function `f` instead of `Plus`.
- `Tr[list, f, n]` considers elements down to level `n`.
- Works for rectangular as well as square matrices and tensors, stopping at the minimum dimension.
- If `n` is omitted, defaults to the depth of the tensor.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/linalg/tr.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/tr.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Tr[{{1, 2}, {3, 4}}]
Out[1]= 5
```

The trace of a symbolic matrix is the sum of its diagonal:

```mathematica
In[1]:= Tr[{{a, b, c}, {d, e, f}, {g, h, i}}]
Out[1]= a + e + i
```

A second argument replaces `Plus` with another combiner — here multiplying the
diagonal entries:

```mathematica
In[1]:= Tr[{{1, 2}, {3, 4}}, Times]
Out[1]= 4
```

The trace of the `n`-th power of the Fibonacci `Q`-matrix is the Lucas number
`L[n]`; for `n = 10` this gives `L[10] = 123`:

```mathematica
In[1]:= Tr[MatrixPower[{{1, 1}, {1, 0}}, 10]]
Out[1]= 123
```

Because the trace is basis-independent, `Tr[A . A]` of a symbolic `2x2` matrix
yields the invariant combination of its entries:

```mathematica
In[1]:= Tr[{{a, b}, {c, d}} . {{a, b}, {c, d}}]
Out[1]= a^2 + 2 b c + d^2
```
