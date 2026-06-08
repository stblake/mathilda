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

`builtin_primeq` is a `*Q` predicate: it always returns `True` or `False`, never unevaluated. For an `EXPR_INTEGER`/`EXPR_BIGINT` it takes `|n|` and runs GMP's `mpz_probab_prime_p(n, 25)` (Baillie–PSW plus 25 Miller–Rabin rounds). With `GaussianIntegers -> True` (parsed by `primeq_parse_options`; a malformed option list yields `False`), a rational integer is a Gaussian prime iff `|n|` is prime and `|n| ≡ 3 (mod 4)`, and a `Complex[a, b]` with integer parts is tested by `gaussian_prime_test` — pure-real/pure-imaginary need the `≡ 3 mod 4` condition, mixed needs `a^2 + b^2` prime. Reals, rationals, strings, symbols, and symbolic functions are all `False`.

- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
