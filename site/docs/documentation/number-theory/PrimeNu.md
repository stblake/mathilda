# PrimeNu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimeNu[n] gives the number of distinct prime factors of n, nu(n). PrimeNu[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts distinct Gaussian prime factors over Z[i]. PrimeNu[1] is 0; PrimeNu[0] is left unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimeNu[24]
Out[1]= 2

In[2]:= PrimeNu[105]
Out[2]= 3

In[3]:= PrimeNu[{4, 28, 180}]
Out[3]= {1, 2, 3}

In[4]:= PrimeNu[50!]
Out[4]= 15

In[5]:= PrimeNu[3 + I]
Out[5]= 2

In[6]:= PrimeNu[105, GaussianIntegers -> True]
Out[6]= 4
```

## Implementation notes

- `Listable`, `Protected`.
- Additive on coprime arguments: `nu(m n) = nu(m) + nu(n)` when

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
