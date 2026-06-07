# Inner

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Inner[f,list1,list2,g]
    is a generalization of Dot in which f plays the role of multiplication and g of addition.
Inner[f,list1,list2]
    uses Plus for g.
Inner[f,list1,list2,g,n]
    contracts index n of the first tensor with the first index of the second tensor.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Inner[f, {a, b}, {x, y}, g]
Out[1]= g[f[a, x], f[b, y]]

In[2]:= Inner[Times, {{a, b}, {c, d}}, {x, y}, Plus]
Out[2]= {a x + b y, c x + d y}

In[3]:= Inner[Times, {a1, a2, a3}, {b1, b2, b3}, Plus]
Out[3]= a1 b1 + a2 b2 + a3 b3
```

## Implementation notes

**Algorithm.** `builtin_inner` is the generalisation of `Dot` that replaces the elementwise multiply with an arbitrary `f` and the summation with an arbitrary `g`: `Inner[f, A, B, g]` contracts the last index of `A` with the first index of `B`, combining matched leaves with `f` and reducing each contraction list with `g` (both default chains are built so that `Inner[Times, A, B, Plus]` reproduces `Dot`). `g` defaults to `Plus` when omitted. The contraction is performed recursively by the `inner_A`/`inner_n1_A` helpers (the `n == 1` form contracts first-index-with-first-index directly), producing an unevaluated `g[f[…], …]` tree that is then run through `evaluate`. Non-`List`-structured operands, or mismatched contraction lengths, return `NULL` (leave unevaluated).

**Data structures.** Operands are walked as nested `EXPR_FUNCTION` trees keyed on `A`'s head (the result head is taken from `A`); leaves are deep-copied into `f`-applications. No flat dense buffer is used — unlike `Dot`, `Inner` works on generic heads and arbitrary combiners.

- `Protected`.
- Like `Dot`, `Inner` effectively contracts the last index of the first tensor with the first index of the second tensor.
- Applying `Inner` to a rank $r$ tensor and a rank $s$ tensor gives a rank $r+s-2$ tensor.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
