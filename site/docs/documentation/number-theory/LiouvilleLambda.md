# LiouvilleLambda

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LiouvilleLambda[n] gives the Liouville function lambda(n) = (-1)^Omega(n), where Omega(n) counts the prime factors of n with multiplicity. Completely multiplicative. A non-real Gaussian-integer argument, or GaussianIntegers -> True, is handled over Z[i].`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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
```

### Options (1)

```mathematica
In[6]:= LiouvilleLambda[8, GaussianIntegers -> True]
Out[6]= 1
```

## Implementation notes

- `Listable`, `Protected`.
- Completely multiplicative: `lambda(m n) = lambda(m) lambda(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `1` or `-1`.
- The sign of `n` is ignored (`lambda(-n) = lambda(n)`).
- Gaussian integers: `LiouvilleLambda[n, GaussianIntegers -> True]`, or a
  non-real Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]`
  and counts the Gaussian prime factors with multiplicity. Because `2` factors
  as `-i (1 + i)^2` in `Z[i]`, e.g. `LiouvilleLambda[2, GaussianIntegers -> True]`
  is `1` while `LiouvilleLambda[2]` is `-1`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `LiouvilleLambda::argt` message.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_liouvillelambda.c`](https://github.com/stblake/mathilda/blob/main/tests/test_liouvillelambda.c)
- Tests: [`tests/test_primenu.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primenu.c)
- Tests: [`tests/test_sum_product_families.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sum_product_families.c)
