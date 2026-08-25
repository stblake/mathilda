# HermiteDecomposition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HermiteDecomposition[m]`**

Gives the Hermite normal form decomposition {u, r} of the integer matrix m: u is unimodular (Abs\[Det\[u\]\] == 1), r is the row Hermite normal form, and u . m == r.  r is in echelon shape with positive pivots and entries above each pivot reduced into \[0, pivot).

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= HermiteDecomposition[{{2, 3, 1}, {4, 1, 5}, {6, 2, 0}}]
Out[1]= {{{0, 2, -1}, {1, 4, -3}, {1, 7, -5}}, {{2, 0, 10}, {0, 1, 21}, {0, 0, 36}}}

In[2]:= h = HermiteDecomposition[{{1, 2, 3}, {4, 5, 6}}]; h[[1]] . {{1, 2, 3}, {4, 5, 6}} == h[[2]]
Out[2]= True
```

## Algorithm

hnf.c -- Hermite Normal Form over Z, and the HermiteDecomposition builtin.

```text
`linalg_hnf` computes a unimodular P and row-HNF R with P*A == R, tracking
every integer row operation on P.  The elimination in each column uses the
```

extended-gcd 2x2 unimodular transform

```text
    [ s   t ] [row_r]      [ g*... ]           s*a_r + t*a_i = g
    [-b   a ] [row_i]  ->  [   0   ]  in col c, a = a_r/g, b = a_i/g,
```

whose determinant is s*a + t*b = (s*a_r + t*a_i)/g = 1, so P stays

```text
unimodular.  Pivots are then made positive and entries above each pivot are
reduced into [0, pivot).  This is the reusable integer primitive behind
```

exact linear Diophantine system solving (src/solve/solveint_linear.c).

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/linalg/linalg.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/linalg.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_latticereduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_latticereduce.c)
