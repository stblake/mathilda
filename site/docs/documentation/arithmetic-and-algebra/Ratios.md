# Ratios

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Ratios[list]
    gives the successive ratios list[[k+1]]/list[[k]] of the elements
    of list (length l - 1).
Ratios[list, n] gives the n-th iterated ratios (length l - n); n must
be a non-negative integer (n = 0 returns list unchanged).
Ratios[list, {n1, n2, ...}] gives the successive n_k-th ratios at
level k of a nested list; for a matrix m, Ratios[m, n] (= Ratios[m, {n, 0}]) takes ratios of successive rows.
FoldList[Times, x, Ratios[list]] inverts Ratios.
Ratios has the attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Ratios[{a, b, c, d, e}]
Out[1]= {b/a, c/b, d/c, e/d}

In[2]:= Ratios[{a, b, c, d, e}, 2]
Out[2]= {(a c)/b^2, (b d)/c^2, (c e)/d^2}

In[3]:= Ratios[Table[2^i, {i, 0, 10}]]
Out[3]= {2, 2, 2, 2, 2, 2, 2, 2, 2, 2}

In[4]:= Ratios[{{a11, a12, a13}, {a21, a22, a23}}, {0, 1}]
Out[4]= {{a12/a11, a13/a12}, {a22/a21, a23/a22}}

In[5]:= FoldList[Times, a, Ratios[{a, b, c, d, e}]]
Out[5]= {a, b, c, d, e}
```

## Implementation notes

- `Protected`.
- `Ratios[list]` divides successive elements by preceding ones, giving `{list[[2]]/list[[1]], list[[3]]/list[[2]], ...}` of length `l - 1`; `Ratios[{x1, x2}]` gives `{x2/x1}`.
- `Ratios[list, n]` applies the ratio operator `n` times, giving length `l - n`. `n` must be a non-negative integer (`n = 0` returns `list` unchanged).
- `Ratios[list, {n1, n2, ...}]` gives the successive `nk`-th ratios at level `k` of a nested list, and is equivalent to `Ratios[Ratios[list, n1], {0, n2, ...}]`. Each `nk` must be a non-negative integer.
- Division threads element-wise over sublists via the `Listable` `Power`/`Times`, so for a matrix `m`, `Ratios[m]` (= `Ratios[m, 1]` = `Ratios[m, {1, 0}]`) takes ratios of successive rows within each column, while `Ratios[m, {0, 1}]` takes ratios of columns within each row.
- The head of the input is preserved. A list shorter than the reduction yields the empty list. First ratios are constant for a geometric sequence.
- Works on machine integers, GMP arbitrary-precision integers (exact `Rational` results), machine-precision doubles, and symbolic expressions.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
