# Array

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Array[f, n]
    generates a list {f[1], f[2], ..., f[n]}.
Array[f, n, r]
    generates a list of length n starting from index r.
Array[f, {n1, n2, ...}]
    generates an n1 x n2 x ... nested-list array with elements
    f[i1, i2, ...].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_array` (in `src/list.c`) builds an N-dimensional list of `f[i1, ..., iN]` applications. It accepts `Array[f, n]`, `Array[f, {n1,...,nN}]` (dimensions), and an optional third argument giving the index origin/range per dimension. The handler normalizes the dimension spec into a `dim_count`-length `n_array` of integer counts and a parallel `r_array` of per-dimension range specs (either a shared scalar, a per-dimension list, or `NULL` for the default origin of 1). All counts must be non-negative integers or the builtin returns `NULL` (unevaluated).

The work is done by the recursive `array_helper`, which descends one dimension per call accumulating index values in `current_args`. At each level it computes the index for slot `i`: for a `{a, b}` range spec it interpolates `a + i*(b-a)/(n-1)` (evaluated symbolically through `Plus`/`Times`/`Divide` so exact rationals survive), otherwise it produces the arithmetic sequence `r_base + i` from the origin. At the deepest level it builds `f[i1, ..., iN]` and calls `evaluate` on it. Each level wraps its children in a `List[...]`.

**Data structures.** Plain `Expr**` working arrays (`n_array`, `r_array`, `current_args`); results assembled bottom-up as nested `List` expressions. Unlike `Table`, `Array` does not bind any iteration symbol — indices are passed positionally as arguments to `f`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Array[f, 5]
Out[1]= {f[1], f[2], f[3], f[4], f[5]}

In[2]:= Array[#^2 &, 4]
Out[2]= {1, 4, 9, 16}

In[3]:= Array[f, 3, 0]
Out[3]= {f[0], f[1], f[2]}
```

A two-argument function over a dimension spec builds a matrix; here the identity matrix as a Kronecker delta:

```mathematica
In[1]:= Array[Boole[#1 == #2] &, {3, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

Building the 4x4 Hilbert matrix and taking its determinant gives the famously tiny rational, evidence of its near-singularity:

```mathematica
In[1]:= Det[Array[1/(#1 + #2 - 1) &, {4, 4}]]
Out[1]= 1/6048000
```

### Notes

`Array[f, n]` builds a length-`n` list by applying `f` to indices `1..n`; an optional third argument sets the starting index. A dimension list `{n1, n2, ...}` produces a nested array with elements `f[i1, i2, ...]`.
