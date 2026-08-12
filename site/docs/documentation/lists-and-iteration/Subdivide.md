# Subdivide

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Subdivide[n]`**

Gives the list {0, 1/n, 2/n, ..., 1} of n + 1 equally spaced points spanning 0 to 1, including both endpoints.

**`Subdivide[max, n]`**

Gives n + 1 equally spaced points spanning 0 to max.

**`Subdivide[min, max, n]`**

Gives n + 1 equally spaced points spanning min to max; point i is min + i (max - min)/n. Descending intervals (min \> max) are allowed and give a negative step. Exact input gives exact results in lowest terms, with both endpoints exact. Returns unevaluated unless n is a positive integer.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Subdivide[4]
Out[1]= {0, 1/4, 1/2, 3/4, 1}

In[2]:= Subdivide[10, 4]
Out[2]= {0, 5/2, 5, 15/2, 10}

In[3]:= Subdivide[2, 8, 3]
Out[3]= {2, 4, 6, 8}

In[4]:= Subdivide[3, 1, 4]
Out[4]= {3, 5/2, 2, 3/2, 1}

In[5]:= Subdivide[a, b, 2]
Out[5]= {a, a + 1/2 (-a + b), b}

In[6]:= Subdivide[0]
Out[6]= Subdivide[0]
```

## Algorithm

Subdivide — equally spaced points spanning an interval, endpoints included.

Mathematica semantics:

```text
  Subdivide[n]             n + 1 points spanning 0 to 1:
                           {0, 1/n, 2/n, ..., 1}
  Subdivide[max, n]        n + 1 points spanning 0 to max
  Subdivide[min, max, n]   n + 1 points spanning min to max
```

Every form yields n + 1 points, because n counts the *parts* the interval is cut into, not the points produced. Point i (0-based) is

```text
    min + i (max - min) / n
```

### Descending Intervals

The formula above is applied verbatim for every point, with no special case. When max < min the quantity (max - min) is negative, so the step is negative and the points descend: Subdivide[3, 1, 4] -> {3, 5/2, 2, 3/2, 1}. Nothing here branches on the sign of the step, compares min against max, sorts, or takes an absolute value, so a descending call returns n + 1 points exactly as an ascending one does. A degenerate interval (min == max) has step 0 and simply repeats the endpoint: Subdivide[5, 5, 2] -> {5, 5, 5}.

### Exactness

Two independent guarantees:

```text
  1. The endpoints are *copied, never computed*. Element 0 is a copy of min
     and element n is a copy of max, so no arithmetic — and therefore no
     representation change — can touch them.
  2. Interior points are each computed directly from their index i, never as
     `previous + step`. With no accumulator there is nothing for error to
     accumulate in.
```

Interior point i is formed as the single fraction

```text
    (min n + i (max - min)) / n
```

and reduced once by make_rational(), which divides through by the gcd, normalizes the sign onto the numerator, and returns a plain Integer when the reduced denominator is 1. That is what lets whole-number points print as integers alongside rationals in one list: Subdivide[10, 4] -> {0, 5/2, 5, 15/2, 10}.

### Validity

n must be a positive machine integer; anything else (0, negative, Rational, Real, symbolic, or a bigint demanding more elements than can be built) leaves the expression unevaluated, per the builtin NULL convention.

### Performance

The int64 fast path below is chosen only when overflow is impossible by construction (see SUBDIVIDE_MAX_ENDPOINT), so it needs no per-operation overflow check. Everything it declines — bigint, rational, real, or symbolic endpoints — falls back to building a Plus/Times/Power expression and letting the core evaluator do the arithmetic, the same strategy src/list/rescale.c uses. That keeps exactness and bigint correctness without duplicating the evaluator's numeric tower here.

## Implementation notes

- Results are **exact** for exact input: rationals come back in lowest terms,
  and a point landing on a whole number prints as an integer, so
  `Subdivide[10, 4]` gives `{0, 5/2, 5, 15/2, 10}`.
- Both **endpoints are exact**. They are copied from the input rather than
  computed, and each interior point is derived directly from its own index
  rather than by accumulating a step, so no rounding can drift.
- **Descending intervals** (`min > max`) are allowed: the formula is applied
  unchanged, giving a negative step. A degenerate interval (`min == max`) has
  step 0 and repeats the endpoint.
- Bigint, rational, real, and purely symbolic endpoints all work; arithmetic
  outside the machine-integer fast path is deferred to the core evaluator. Real
  endpoints give machine reals, as elsewhere in the system.
- **A machine-real interval makes every point a machine real**, the endpoints
  included, so `Subdivide[0, 1., 4]` is `{0., 0.25, 0.5, 0.75, 1.}` — the exact
  `0` is carried across as `0.` — and `Subdivide[1/2, 1., 4]` is
  `{0.5, 0.625, 0.75, 0.875, 1.}`. Endpoints are still taken from the input
  rather than computed; the interval's exactness decides *which* value gets
  taken. Only **machine** `Real` is contagious: an MPFR endpoint keeps the exact
  one (``Subdivide[0, 1.`30, 4]`` starts at `0`) and a symbolic endpoint keeps
  the general path.
- Such a result is a [packed list](../packed-arrays/index.md) built directly, with no
  `Expr` node constructed at all. The interior points are `min + i (max - min)/n`
  in doubles, which is the rule `numpy.linspace` uses and which Mathematica
  agrees with: `Subdivide[0., 1., 10]` gives `0.30000000000000004` for its
  fourth point, not `0.3`.
- `n` must be a **positive integer**. Zero, negative, rational, real, and
  symbolic `n` leave the expression unevaluated, as does an `n` large enough
  that the result list would be unbuildable.

**Attributes:** `Protected`.

## See also

[Range](../../lists-and-iteration/Range/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
