# BernoulliB

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BernoulliB[n]
    gives the Bernoulli number B_n.
BernoulliB[n, x]
    gives the Bernoulli polynomial B_n(x).
Non-negative integer n gives the exact rational B_n (odd n > 1 give 0,
B_0 = 1, B_1 = -1/2); an inexact integer-valued n evaluates it at machine
or arbitrary (MPFR) precision. BernoulliB[n, x] expands the degree-n
polynomial with exact rational coefficients, staying symbolic in x or
evaluating numerically when x is inexact. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[BernoulliB[k], {k, 0, 10}]
Out[1]= {1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66}

In[2]:= BernoulliB[3, z]
Out[2]= 1/2 z - 3/2 z^2 + z^3
```

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
