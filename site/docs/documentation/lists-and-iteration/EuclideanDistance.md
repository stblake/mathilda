# EuclideanDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EuclideanDistance[u, v]`**

Gives the Euclidean distance Sqrt\[Sum Abs\[u\_i - v\_i\]^2\] between two equal-length numeric vectors, or between two scalars. Abs makes complex components use their modulus. Exact input gives an exact result, which for a root is usually a Sqrt; use SquaredEuclideanDistance to stay rational. Returns unevaluated for mismatched lengths or matrix arguments.

## Examples (6)

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

## Implementation notes

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

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
