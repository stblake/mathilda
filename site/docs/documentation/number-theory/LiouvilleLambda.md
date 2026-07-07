# LiouvilleLambda

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
LiouvilleLambda[n] gives the Liouville function lambda(n) = (-1)^Omega(n), where Omega(n) counts the prime factors of n with multiplicity. Completely multiplicative. A non-real Gaussian-integer argument, or GaussianIntegers -> True, is handled over Z[i].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= LiouvilleLambda[8]
Out[1]= -1

In[2]:= LiouvilleLambda[9]
Out[2]= 1

In[3]:= LiouvilleLambda[{1, 2, 3, 4, 5, 6}]
Out[3]= {1, -1, -1, 1, -1, 1}

In[4]:= LiouvilleLambda[10^30 + 1]
Out[4]= -1

In[5]:= LiouvilleLambda[2 + I]
Out[5]= -1

In[6]:= LiouvilleLambda[8, GaussianIntegers -> True]
Out[6]= 1
```

## Implementation notes

- `Listable`, `Protected`.
- Completely multiplicative: `lambda(m n) = lambda(m) lambda(n)`.
- Computed directly from the prime factorisation (machine integers and GMP

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
