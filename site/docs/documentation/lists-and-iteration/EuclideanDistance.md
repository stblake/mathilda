# EuclideanDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EuclideanDistance[u, v]`**

Gives the Euclidean distance Sqrt\[Sum Abs\[u\_i - v\_i\]^2\] between two equal-length numeric vectors, or between two scalars. Abs makes complex components use their modulus. Exact input gives an exact result, which for a root is usually a Sqrt; use SquaredEuclideanDistance to stay rational. Returns unevaluated for mismatched lengths or matrix arguments.

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

The 3-4-5 right triangle

```mathematica
In[7]:= EuclideanDistance[{0, 0}, {3, 4}]
Out[7]= 5
```

Straight-line distance in 3-D

```mathematica
In[8]:= EuclideanDistance[{1, 1, 1}, {4, 5, 1}]
Out[8]= 5
```

Scalars are treated as 1-vectors

```mathematica
In[9]:= EuclideanDistance[3, 7]
Out[9]= 4
```

## Implementation notes

**Algorithm.** `builtin_euclidean_distance` computes
`Sqrt[Sum_i Abs[u_i - v_i]^2]` via the shared `dist_builtin(res, p=2, root=true)`.
The core `dist_sum` forms `Sum_i Abs[u_i - v_i]^p` by composing the internal
arithmetic primitives (`internal_subtract`, an Abs helper, `internal_power`,
`internal_plus`) through `eval_and_free`, and `EuclideanDistance` then takes the
square root of that sum. Reusing the evaluator's arithmetic rather than
recomputing it gives three properties for free: exact input stays exact where the
result is rational; complex components use their modulus (`Abs`-then-square, as
Mathematica defines it); and symbolic input survives as a symbolic distance
rather than being rejected.

**Shape and the Abs gap.** `dist_shape` admits only two scalars, or two `List`s
of equal length (a length mismatch, or a list-valued component, declines — these
are vector functions, not matrix-threaded). The internal `dist_abs` decides a
real argument's sign with `list_numeric_sign` and negates if needed, routing
around a bug where `Abs` declines on a bigint-scale rational, which otherwise
left high-precision exact distances unevaluated; complex/symbolic arguments still
go through `internal_abs` for the modulus.

**Complexity / limits.** O(n) arithmetic evaluations over the `n` components;
interpreter-speed, not a buffer path. `ATTR_PROTECTED`. See also
`SquaredEuclideanDistance` (the root-free, exactness-preserving form that
`FindClusters` ranks on), `ManhattanDistance`, and `CosineDistance`.

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

**See also:** [SquaredEuclideanDistance](../../lists-and-iteration/SquaredEuclideanDistance/), [ManhattanDistance](../../lists-and-iteration/ManhattanDistance/), [CosineDistance](../../lists-and-iteration/CosineDistance/), [List](../../other-advanced/List/), [FindClusters](../../lists-and-iteration/FindClusters/), [Abs](../../arithmetic/Abs/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`EuclideanDistance[u, v]` is `Sqrt[Sum Abs[u_i - v_i]^2]`, the ordinary
straight-line distance. Both arguments must be scalars, or lists of equal length;
a length mismatch or matrix-shaped input is left unevaluated. Complex components
contribute their modulus (the definition squares `Abs`, not the raw difference),
and a symbolic pair returns a symbolic distance rather than an error.

For ranking or clustering, prefer `SquaredEuclideanDistance`: it avoids the root,
stays exact for exact input, and orders points identically.
