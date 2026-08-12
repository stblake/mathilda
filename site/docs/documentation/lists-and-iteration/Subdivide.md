# Subdivide

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Subdivide[n]
    Gives the list {0, 1/n, 2/n, ..., 1} of n + 1 equally
    spaced points spanning 0 to 1, including both endpoints.
Subdivide[max, n]
    Gives n + 1 equally spaced points spanning 0 to
    max.
Subdivide[min, max, n]
    Gives n + 1 equally spaced points spanning
    min to max; point i is min + i (max - min)/n. Descending intervals
    (min > max) are allowed and give a negative step. Exact input gives
    exact results in lowest terms, with both endpoints exact. Returns
    unevaluated unless n is a positive integer.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

- Results are **exact** for exact input: rationals come back in lowest terms,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
