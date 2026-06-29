# DivisorSigma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DivisorSigma[k, n] gives the divisor function sigma_k(n), the sum of the k-th powers of the divisors of n. DivisorSigma[k, n, GaussianIntegers -> True] sums over Gaussian-integer divisors.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DivisorSigma[1, 20]
Out[1]= 42

In[2]:= DivisorSigma[2, 20]
Out[2]= 546

In[3]:= DivisorSigma[0, 12]
Out[3]= 6

In[4]:= DivisorSigma[-2, 10]
Out[4]= 13/10

In[5]:= DivisorSigma[1/2, 12]
Out[5]= (2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))

In[6]:= DivisorSigma[k, {2, 3, 6}]
Out[6]= {(-1 + 2^(2 k))/(-1 + 2^k), (-1 + 3^(2 k))/(-1 + 3^k), ((-1 + 2^(2 k)) (-1 + 3^(2 k)))/((-1 + 2^k) (-1 + 3^k))}

In[7]:= DivisorSigma[2, {1, 2, 3, 4, 5}]
Out[7]= {1, 5, 10, 21, 26}

In[8]:= DivisorSigma[1, 3 + I]
Out[8]= 2 + 6*I
```

## Implementation notes

- `Listable`, `NHoldAll`, `Protected`.
- Computed from the multiplicative formula

**Attributes:** `Listable`, `NHoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
