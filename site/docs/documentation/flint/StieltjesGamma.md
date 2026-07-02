# FLINT`StieltjesGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FLINT`StieltjesGamma[n] and FLINT`StieltjesGamma[n, a] give the numeric value of the n-th Stieltjes constant (generalized, at a) for a non-negative integer n, via FLINT (acb_dirichlet_stieltjes). Unevaluated for negative or non-integer n.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FLINT`StieltjesGamma[0]
Out[1]= 0.5772156649015329
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/flint_num_bridge.c`](https://github.com/stblake/mathilda/blob/main/src/flint_num_bridge.c)
- Specification: [`docs/spec/builtins/flint.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/flint.md)
