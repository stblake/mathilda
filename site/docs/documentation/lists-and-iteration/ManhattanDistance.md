# ManhattanDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ManhattanDistance[u, v]`**

Gives Sum Abs\[u\_i - v\_i\], the sum of component-wise distances. Differs from EuclideanDistance in two or more dimensions; in one dimension the two agree.

## Examples (8)

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

### Applications (2)

Sum of coordinate gaps: 3 + 4 = 7

```mathematica
In[7]:= ManhattanDistance[{0, 0}, {3, 4}]
Out[7]= 7
```

A symbolic pair stays symbolic

```mathematica
In[8]:= ManhattanDistance[{a}, {b}]
Out[8]= Abs[a - b]
```

## Implementation notes

**Algorithm.** `builtin_manhattan_distance` computes `Sum_i Abs[u_i - v_i]` via
the shared `dist_builtin(res, p=1, root=false)`. `dist_sum` forms the sum of
absolute component differences by composing `internal_subtract`, the `dist_abs`
helper and `internal_plus` through `eval_and_free`, so exact input stays exact
(`ManhattanDistance[{1, 2}, {4, 6}]` is `7`), complex components contribute their
modulus, and a symbolic pair comes back symbolically
(`ManhattanDistance[{a}, {b}]` is `Abs[a - b]`, as Mathematica answers).

**Shape.** `dist_shape` admits two scalars or two equal-length `List`s only; a
length mismatch, a list-valued component, or a non-atomic "scalar" (a `Rational`
is stored as `Rational[...]`, an `EXPR_FUNCTION`, so a bare rational scalar is
*not* accepted) declines. `dist_abs` decides a real sign with `list_numeric_sign`
to route around the `Abs` bigint-rational gap.

**Complexity / limits.** O(n) arithmetic evaluations over the components;
interpreter-speed. `ATTR_PROTECTED`. See `EuclideanDistance`,
`SquaredEuclideanDistance`, and `CosineDistance` — all share the `dist_sum` loop
and differ only in `p` and whether a root is taken.

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

**See also:** [EuclideanDistance](../../lists-and-iteration/EuclideanDistance/), [SquaredEuclideanDistance](../../lists-and-iteration/SquaredEuclideanDistance/), [CosineDistance](../../lists-and-iteration/CosineDistance/), [List](../../other-advanced/List/), [FindClusters](../../lists-and-iteration/FindClusters/), [Abs](../../arithmetic/Abs/)

- Source: [`src/list/distance.c`](https://github.com/stblake/mathilda/blob/main/src/list/distance.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

`ManhattanDistance[u, v]` is `Sum Abs[u_i - v_i]`, the taxicab (city-block)
distance — the distance travelled along axis-aligned moves rather than in a
straight line. Both arguments must be scalars or lists of equal length.

Exact input stays exact, complex components contribute their modulus, and a
symbolic pair returns `Abs[a - b]` rather than an error, exactly as Mathematica
answers. Note that a bare rational *scalar* such as `ManhattanDistance[1/2, 3/2]`
is left unevaluated, because a `Rational` is a compound expression rather than an
atomic scalar — wrap the values in lists (`ManhattanDistance[{1/2}, {3/2}]`) to
measure them.
