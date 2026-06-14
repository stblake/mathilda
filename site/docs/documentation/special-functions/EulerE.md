# EulerE

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
EulerE[n]
    gives the Euler number E_n.
EulerE[n, x]
    gives the Euler polynomial E_n(x).
Non-negative integer n gives the exact integer E_n (odd n give 0,
E_0 = 1, E_2 = -1, E_4 = 5); an inexact integer-valued n evaluates it at
machine or arbitrary (MPFR) precision. EulerE[n, x] expands the degree-n
polynomial with exact rational coefficients, staying symbolic in x or
evaluating numerically when x is inexact; EulerE[n, 1/2] folds to
2^-n EulerE[n]. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[EulerE[k], {k, 0, 10}]
Out[1]= {1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521}

In[2]:= EulerE[3, z]
Out[2]= 1/4 - 3/2 z^2 + z^3
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
