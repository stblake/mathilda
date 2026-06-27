# ProductLog

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ProductLog[z]
    gives the principal solution w of z == w e^w (the Lambert W function).
ProductLog[k, z] gives the k-th solution (k any integer, k == 0 the
principal branch); branches are ordered by imaginary part. ProductLog[z]
is real for z >= -1/e and has a branch cut along (-Infinity, -1/e].
Exact values include ProductLog[0] = 0, ProductLog[-1/E] = -1,
ProductLog[E] = 1, ProductLog[-Pi/2] = I Pi/2 and ProductLog[k, 0] =
-Infinity for k != 0. Inexact real or complex arguments evaluate
numerically at machine or arbitrary (MPFR) precision. Satisfies
D[ProductLog[z], z] = ProductLog[z]/(z (1 + ProductLog[z])). Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ProductLog[1.0]
Out[1]= 0.567143

In[2]:= ProductLog[-1/E]
Out[2]= -1
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
