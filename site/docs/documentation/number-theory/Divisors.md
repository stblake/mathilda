# Divisors

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Divisors[n] gives a list of the integers that divide n. Divisors[n, GaussianIntegers -> True] includes Gaussian-integer divisors.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Divisors[1729]
Out[1]= {1, 7, 13, 19, 91, 133, 247, 1729}

In[2]:= Divisors[6]
Out[2]= {1, 2, 3, 6}

In[3]:= Divisors[{605, 871, 824}]
Out[3]= {{1, 5, 11, 55, 121, 605}, {1, 13, 67, 871}, {1, 2, 4, 8, 103, 206, 412, 824}}

In[4]:= Divisors[6 + 4 I]
Out[4]= {1, 1 + I, 1 + 5*I, 2, 3 + 2*I, 6 + 4*I}

In[5]:= Divisors[2, GaussianIntegers -> True]
Out[5]= {1, 1 + I, 2}

In[6]:= Divisors[3, GaussianIntegers -> True]
Out[6]= {1, 3}
```

## Implementation notes

- `Listable`, `Protected`.
- Machine integers and GMP bigints are handled uniformly; the result promotes to

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
