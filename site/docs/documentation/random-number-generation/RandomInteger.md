# RandomInteger

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RandomInteger[{imin, imax}]
    gives a pseudorandom integer in the range {imin, ..., imax}.
RandomInteger[imax]
    gives a pseudorandom integer in the range {0, ..., imax}.
RandomInteger[]
    pseudorandomly gives 0 or 1.
RandomInteger[range, n]
    gives a list of n pseudorandom integers.
RandomInteger[range, {n1, n2, ...}]
    gives an n1 x n2 x ... array of pseudorandom integers.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SeedRandom[42]; RandomInteger[]
Out[1]= 1

In[2]:= SeedRandom[42]; RandomInteger[10]
Out[2]= 9

In[3]:= SeedRandom[42]; RandomInteger[{1, 6}]
Out[3]= 2

In[4]:= SeedRandom[42]; RandomInteger[{0, 9}, 5]
Out[4]= {9, 9, 7, 1, 9}

In[5]:= SeedRandom[42]; Dimensions[RandomInteger[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomInteger[{-10, -5}]
Out[6]= -9

In[7]:= IntegerQ[RandomInteger[10^20]]
Out[7]= True
```

## Implementation notes

- `Protected`.
- RandomInteger[{imin, imax}] chooses integers in the range {imin, ..., imax} with equal probability.
- RandomInteger[] gives 0 or 1 with probability 1/2.
- RandomInteger gives a different sequence of pseudorandom integers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Returns bignums when the range exceeds 64-bit integer limits.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
