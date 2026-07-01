# MoebiusMu

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MoebiusMu[n] gives the Moebius function mu(n): 0 if n has a squared prime factor, otherwise (-1)^k where k is the number of distinct primes. A non-real Gaussian-integer argument is handled over Z[i].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= MoebiusMu[11]
Out[1]= -1

In[2]:= MoebiusMu[10]
Out[2]= 1

In[3]:= MoebiusMu[1440]
Out[3]= 0

In[4]:= MoebiusMu[{4, 10, 17, 20}]
Out[4]= {0, 1, -1, 0}

In[5]:= MoebiusMu[10^50 + 1]
Out[5]= -1

In[6]:= MoebiusMu[5 + 6 I]
Out[6]= -1
```

## Implementation notes

- `Listable`, `Protected`.
- Computed directly from the prime factorisation (machine integers and GMP

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
