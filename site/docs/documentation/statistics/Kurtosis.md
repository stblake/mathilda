# Kurtosis

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Kurtosis[data]
    gives the coefficient of kurtosis (peak/tail vs flank concentration) of data, equivalent to CentralMoment[data, 4] / CentralMoment[data, 2]^2. For a matrix it is taken columnwise.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Kurtosis[{1, 2, 3, 4, 5}]
Out[1]= 17/10

In[2]:= Kurtosis[{1, 2, 4, 8}]
Out[2]= 25141/13225
```

## Implementation notes

- `Protected`.
- Equivalent to `CentralMoment[data, 4] / CentralMoment[data, 2]^2` (Pearson kurtosis, not the excess form).
- For a matrix, gives the columnwise kurtoses.
- Handles numerical and symbolic data; exact input gives exact output.
- Fast path on `NDArray`/packed real buffers (`ndred_kurtosis`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
