# Differences

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Differences[list]
    gives the successive differences of the elements of list.
Differences[list, n] gives the n-th differences (length l - n).
Differences[list, n, s] takes differences of elements step s apart
(length l - n |s|).
Differences[list, {n1, n2, ...}] gives the successive n_k-th differences
at level k of a nested list; for a matrix m, Differences[m, n] (= Differences[m, {n, 0}]) differences successive rows.
FoldList[Plus, x, Differences[list]] inverts Differences.
Differences has the attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Differences[{a, b, c, d, e}]
Out[1]= {-a + b, -b + c, -c + d, -d + e}

In[2]:= Differences[{a, b, c, d, e}, 2]
Out[2]= {a - 2 b + c, b - 2 c + d, c - 2 d + e}

In[3]:= Differences[{1, 2, 4, 8, 16, 32}, 1, 2]
Out[3]= {3, 6, 12, 24}

In[4]:= Differences[{{a11, a12, a13}, {a21, a22, a23}}, {0, 1}]
Out[4]= {{-a11 + a12, -a12 + a13}, {-a21 + a22, -a22 + a23}}

In[5]:= FoldList[Plus, a, Differences[{a, b, c, d, e}]]
Out[5]= {a, b, c, d, e}
```

## Implementation notes

- `Protected`.
- `Differences[list]` gives `{list[[2]] - list[[1]], list[[3]] - list[[2]], ...}`, of length `l - 1`.
- `Differences[list, n]` applies the first-difference operator `n` times, giving length `l - n`. `n` must be a non-negative integer (`n = 0` returns `list` unchanged).
- `Differences[list, n, s]` takes differences of elements step `s` apart, of length `l - n |s|`. The step `s` is a nonzero integer; for `s < 0` the elements are subtracted in the opposite order.
- `Differences[list, {n1, n2, ...}]` gives the successive `nk`-th differences at level `k` of a nested list, and is equivalent to `Differences[Differences[list, n1], {0, n2, ...}]`. Each `nk` must be a non-negative integer.
- Subtraction threads element-wise over sublists via the `Listable` `Plus`/`Times`, so for a matrix `m`, `Differences[m]` (= `Differences[m, 1]` = `Differences[m, {1, 0}]`) differences successive rows within each column, while `Differences[m, {0, 1}]` differences columns within each row.
- The head of the input is preserved. A list shorter than the reduction yields the empty list.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
