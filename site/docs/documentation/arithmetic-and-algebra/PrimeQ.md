# PrimeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimeQ[n]
    gives True if n is a prime integer, False otherwise.
PrimeQ[z]
    for a Gaussian integer z = a + b I, gives True if z is a Gaussian prime.
PrimeQ[n, GaussianIntegers -> True]
    tests primality of n in Z[i] rather than in Z.
Primality is tested with GMP's mpz_probab_prime_p using 25 Miller-Rabin rounds on top of a Baillie-PSW pre-screen, so composite false positives have probability below 4^-25 (definite for n < 2^64).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimeQ[7]
Out[1]= True

In[2]:= PrimeQ[1 + I]
Out[2]= True

In[3]:= PrimeQ[1 + 2 I]
Out[3]= True

In[4]:= PrimeQ[3 I]
Out[4]= True

In[5]:= PrimeQ[5 I]
Out[5]= False

In[6]:= PrimeQ[2 + 2 I]
Out[6]= False

In[7]:= PrimeQ[5, GaussianIntegers -> True]
Out[7]= False

In[8]:= PrimeQ[3, GaussianIntegers -> True]
Out[8]= True
```

## Implementation notes

- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
