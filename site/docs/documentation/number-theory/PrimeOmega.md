# PrimeOmega

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimeOmega[n] gives the number of prime factors of n counted with multiplicity, Omega(n). PrimeOmega[n, GaussianIntegers -> True] (or a non-real Gaussian-integer n) counts Gaussian prime factors over Z[i]. PrimeOmega[1] is 0; PrimeOmega[0] is left unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimeOmega[30]
Out[1]= 3

In[2]:= PrimeOmega[12]
Out[2]= 3

In[3]:= PrimeOmega[{4, 12, 24}]
Out[3]= {2, 3, 4}

In[4]:= PrimeOmega[30!]
Out[4]= 59

In[5]:= PrimeOmega[5 + 9 I]
Out[5]= 2

In[6]:= PrimeOmega[12, GaussianIntegers -> True]
Out[6]= 5
```

## Implementation notes

- `Listable`, `Protected`.
- Completely additive: `Omega(m n) = Omega(m) + Omega(n)`.
- Computed directly from the prime factorisation (machine integers and GMP

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
