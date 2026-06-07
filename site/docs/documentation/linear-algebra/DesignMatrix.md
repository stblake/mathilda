# DesignMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DesignMatrix[data, {f1, ..., fn}, vars] gives the design matrix with entries f_i evaluated at the data coordinates.
Data shapes match Fit. The WorkingPrecision option converts entries to machine or n-digit reals; otherwise they are exact.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DesignMatrix[{{0,1},{1,0},{3,2},{5,4}}, {1, x}, x]
Out[1]= {{1, 0}, {1, 1}, {1, 3}, {1, 5}}

In[2]:= DesignMatrix[{{0,0,0},{1,0,1},{0,1,2}}, {1, x, y}, {x, y}]
Out[2]= {{1, 0, 0}, {1, 1, 0}, {1, 0, 1}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
