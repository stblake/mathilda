# FLINT`PolyGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`PolyGamma[n, z] gives the numeric value of the n-th derivative of the digamma function (n = 0 is digamma) via FLINT (acb_polygamma), to the precision of the arguments. Unevaluated for symbolic arguments or at a pole.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`PolyGamma[1, 2]
Out[1]= 0.6449340668482264
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/flint_num_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/flint_num_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
