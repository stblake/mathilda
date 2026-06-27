# Hyperfactorial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Hyperfactorial[n]
    gives the hyperfactorial prod_{k=1}^{n} k^k.
Exact (GMP) for a non-negative integer n; other orders stay symbolic. Listable, NumericFunction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Hyperfactorial[4]
Out[1]= 27648
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
