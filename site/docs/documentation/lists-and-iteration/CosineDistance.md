# CosineDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CosineDistance[u, v]`**

Gives 1 - (u . Conjugate\[v\]) / (Norm\[u\] Norm\[v\]), the angular distance between two vectors: 0 when parallel, 1 when orthogonal and 2 when antiparallel, ignoring magnitude. A zero vector gives 0. Not a metric -- it violates the triangle inequality -- so unlike the Euclidean family it has no squared form that ranks identically.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= EuclideanDistance[{1, 2}, {4, 6}]
Out[1]= 5

In[2]:= SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]
Out[2]= 58/441

In[3]:= ManhattanDistance[{1, 2}, {4, 6}]
Out[3]= 7

In[4]:= EuclideanDistance[{0, 0}, {1, 1}]
Out[4]= Sqrt[2]

In[5]:= CosineDistance[{1, 0}, {0, 1}]
Out[5]= 1

In[6]:= CosineDistance[{1, 0}, {-1, 0}]
Out[6]= 2
```

### Applications (3)

45 degrees apart: 1 - 1/Sqrt[2]

```mathematica
In[7]:= CosineDistance[{1, 0}, {1, 1}]
Out[7]= 1 - 1/Sqrt[2]
```

Parallel vectors give 0

```mathematica
In[8]:= CosineDistance[{1, 2, 3}, {1, 2, 3}]
Out[8]= 0
```

A zero vector gives 0 by convention

```mathematica
In[9]:= CosineDistance[{0, 0}, {1, 2}]
Out[9]= 0
```

## Implementation notes

**Algorithm.** `builtin_cosine_distance` computes
`1 - (u . Conjugate[v]) / (Norm[u] Norm[v])`. `dist_dot_conj` forms the numerator
`Sum_i u_i Conjugate[v_i]` (the `Conjugate` makes it correct for complex vectors
and is a no-op on reals, as Mathematica writes it), and `dist_norm` gives each
Euclidean norm as `Sqrt[Sum Abs[u_i]^2]`. All of it is composed from the internal
arithmetic primitives through `eval_and_free`, so exact input stays exact and
symbolic input survives.

**Range and conventions.** The value runs over `[0, 2]`: `0` for parallel, `1`
for orthogonal, `2` for antiparallel. Unlike the Euclidean family this is *not* a
metric (it violates the triangle inequality) and has no squared form that ranks
identically, so callers use it directly. A zero vector on either side gives `0`
(a special case, since the quotient would be `0/0` → `Indeterminate`), matching
Mathematica.

**Shape / limits.** `dist_shape` admits two scalars or two equal-length `List`s;
a length mismatch or a list-valued component declines. O(n) arithmetic
evaluations, interpreter-speed. `ATTR_PROTECTED`. Shares the `dist_sum`/`dist_norm`
machinery with `EuclideanDistance`, `SquaredEuclideanDistance` and
`ManhattanDistance`.

- `Protected`. Not `Listable`: threading over a `List` argument is exactly what
  these must not do, because the list *is* the point.
- **Exact input gives an exact result** where the value is rational.
  `SquaredEuclideanDistance[{1, 2}, {4, 6}]` is `25`, not `25.`, and
  `SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]` is `58/441`. Squared Euclidean
  is monotone in Euclidean, so ranking on it orders points identically without
  introducing a root -- which is how `FindClusters` stays exact in n dimensions.
- **Complex components contribute their modulus**, because the definition takes
  `Abs` before squaring rather than squaring the difference. This matters only
  for complex input, where the two orders differ, and follows Mathematica.
- Symbolic input survives rather than being rejected: `ManhattanDistance[{a},
  {b}]` is `Abs[a - b]`, as in Mathematica.
- `CosineDistance` ranges over `[0, 2]` -- `0` parallel, `1` orthogonal, `2`
  antiparallel -- and ignores magnitude. It is **not** a metric (it violates the
  triangle inequality) and has no squared form that ranks identically, so it is
  used directly. A zero vector on either side gives `0`, following Mathematica;
  that is a convention, not a derivation, since the quotient is `0/0`.
- Mismatched lengths, or an argument that is a matrix, leave the call
  unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [EuclideanDistance](../../lists-and-iteration/EuclideanDistance/), [SquaredEuclideanDistance](../../lists-and-iteration/SquaredEuclideanDistance/), [ManhattanDistance](../../lists-and-iteration/ManhattanDistance/), [List](../../other-advanced/List/), [FindClusters](../../lists-and-iteration/FindClusters/), [Abs](../../arithmetic/Abs/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`CosineDistance[u, v]` is `1 - (u . Conjugate[v]) / (Norm[u] Norm[v])`, the
angular distance between two vectors. It runs over `[0, 2]`: `0` for parallel,
`1` for orthogonal, `2` for antiparallel. The `Conjugate` makes it correct for
complex vectors and is a no-op on reals.

Unlike the Euclidean family this is **not** a metric — it ignores magnitude and
violates the triangle inequality — so there is no squared form that ranks
identically. A zero vector on either side gives `0` (the quotient would
otherwise be the indeterminate `0/0`). Exact input stays exact and symbolic
vectors pass through.
