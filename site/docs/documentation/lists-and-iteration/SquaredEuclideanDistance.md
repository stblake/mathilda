# SquaredEuclideanDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SquaredEuclideanDistance[u, v]`**

Gives Sum Abs\[u\_i - v\_i\]^2, the squared Euclidean distance. Rational for rational input, and monotone in EuclideanDistance, so ranking on it orders points identically without taking a root.

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

The squared 3-4-5 distance

```mathematica
In[7]:= SquaredEuclideanDistance[{0, 0}, {3, 4}]
Out[7]= 25
```

1 + 4 + 4 = 9

```mathematica
In[8]:= SquaredEuclideanDistance[{1, 2, 2}, {0, 0, 0}]
Out[8]= 9
```

Scalars act as 1-vectors

```mathematica
In[9]:= SquaredEuclideanDistance[3, 8]
Out[9]= 25
```

## Implementation notes

**Algorithm.** `builtin_squared_euclidean_distance` computes
`Sum_i Abs[u_i - v_i]^2` via the shared `dist_builtin(res, p=2, root=false)` —
the Euclidean distance without the final square root. The sum is assembled by
`dist_sum` from the internal arithmetic primitives through `eval_and_free`, so a
real squared term skips the redundant `Abs`, a complex term uses its modulus
(`Abs`-then-square), and a symbolic term survives.

**Why the squared form matters.** Because no root is taken, the result is
*rational for rational input* (`SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]` is
`58/441`, not a float), and squaring is monotone on non-negatives. So ranking on
the squared distance orders points identically to ranking on the true distance
without ever introducing an irrational — which is exactly what lets
`FindClusters` partition n-dimensional exact data exactly (it is the metric the
spanning-tree builder ranks on).

**Shape / limits.** `dist_shape` admits two scalars or two equal-length `List`s;
a length mismatch or a list-valued component declines. O(n) arithmetic
evaluations, interpreter-speed. `ATTR_PROTECTED`. See `EuclideanDistance` for the
rooted form and `ManhattanDistance` / `CosineDistance` for the siblings.

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

**See also:** [EuclideanDistance](../../lists-and-iteration/EuclideanDistance/), [ManhattanDistance](../../lists-and-iteration/ManhattanDistance/), [CosineDistance](../../lists-and-iteration/CosineDistance/), [List](../../other-advanced/List/), [FindClusters](../../lists-and-iteration/FindClusters/), [Abs](../../arithmetic/Abs/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`SquaredEuclideanDistance[u, v]` is `Sum Abs[u_i - v_i]^2` — the Euclidean
distance without the final square root. Because no root is taken, the result is
exact for exact input (`SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]` is
`58/441`, not a float), and since squaring is monotone on non-negatives, ranking
on the square orders points identically to ranking on the true distance. That is
exactly what lets `FindClusters` partition exact multi-dimensional data without
ever introducing an irrational.

Both arguments must be scalars or equal-length lists; complex components use
their modulus and symbolic input survives.
