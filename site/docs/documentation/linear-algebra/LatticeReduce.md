# LatticeReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LatticeReduce[m]
    gives an LLL-reduced basis for the lattice spanned by the rows
    (vectors) of m.  The entries of m may be integers, Gaussian
    integers, rationals, or Gaussian rationals.  Reduction is exact
    (GMP rational arithmetic, so it is correct for both machine-size
    and arbitrary-precision entries) and preserves the lattice, its
    determinant, and every linear relation among the rows.  The rows
    must be linearly independent.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LatticeReduce[{{1, 0, 0, 1345}, {0, 1, 0, 35}, {0, 0, 1, 154}}]
Out[1]= {{0, 9, -2, 7}, {1, 1, -9, -6}, {1, -3, -8, 8}}

In[2]:= {w1, w2} = LatticeReduce[{{12, 2}, {13, 4}}]
Out[2]= {{1, 2}, {9, -4}}

In[3]:= b . {1, 2, 3, 1}    (* relations preserved *)
Out[3]= Dot[b, {1, 2, 3, 1}]

In[4]:= LatticeReduce[{{1, 2}, {3, 4.5}}]
Out[4]= LatticeReduce[{{1, 2}, {3, 4.5}}]
```

## Implementation notes

- `Protected`.
- Returns an `n × d` matrix whose rows form a reduced basis of the same

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
