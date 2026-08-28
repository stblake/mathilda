# PrimeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimeQ[n] gives True if n is a prime integer and False otherwise, using a Baillie-PSW and Miller-Rabin probabilistic primality test (no known counterexample, and definite below 2^64).`**

**`PrimeQ[z]`**

for a Gaussian integer z = a + b I, gives True if z is a Gaussian prime.

**`PrimeQ[n, GaussianIntegers -> True]`**

tests primality of n in Z\[i\] rather than in Z.

<details>
<summary>Notes</summary>

Primality is tested with GMP's mpz\_probab\_prime\_p using 25 Miller-Rabin rounds on top of a Baillie-PSW pre-screen, so composite false positives have probability below 4^-25 (definite for n \< 2^64).

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

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

In[7]:= PrimeQ[Exp[2 Pi I/3]]
Out[7]= False
```

### Options (2)

```mathematica
In[8]:= PrimeQ[5, GaussianIntegers -> True]
Out[8]= False

In[9]:= PrimeQ[3, GaussianIntegers -> True]
Out[9]= True
```

### Applications (7)

```mathematica
In[10]:= PrimeQ[97]
Out[10]= True

In[11]:= PrimeQ[2^31 - 1]
Out[11]= True

In[12]:= PrimeQ[2^67 - 1]
Out[12]= False

In[13]:= PrimeQ[561]
Out[13]= False

In[14]:= PrimeQ[5, GaussianIntegers -> True]
Out[14]= False

In[15]:= PrimeQ[3, GaussianIntegers -> True]
Out[15]= True

In[16]:= PrimeQ[2 + 3 I]
Out[16]= True
```

## Implementation notes

`builtin_primeq` is a `*Q` predicate: it always returns `True` or `False`, never unevaluated. For an `EXPR_INTEGER`/`EXPR_BIGINT` it takes `|n|` and runs GMP's `mpz_probab_prime_p(n, 25)` (Baillie–PSW plus 25 Miller–Rabin rounds). With `GaussianIntegers -> True` (parsed by `primeq_parse_options`; a malformed option list yields `False`), a rational integer is a Gaussian prime iff `|n|` is prime and `|n| ≡ 3 (mod 4)`, and a `Complex[a, b]` with integer parts is tested by `gaussian_prime_test` — pure-real/pure-imaginary need the `≡ 3 mod 4` condition, mixed needs `a^2 + b^2` prime. Reals, rationals, strings, symbols, and symbolic functions are all `False`.

- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian
  inputs (symbols, `Sqrt[2]`, `Exp[2 Pi I/3]`, strings, etc.) returns
  `False` — `*Q` predicates never remain symbolic.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.

**Attributes:** `Listable`, `Protected`.

## References

- R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective*, 2nd ed., Springer, 2005 — Miller–Rabin (§3.5) and the Baillie–PSW test.
- R. Baillie and S. S. Wagstaff Jr., "Lucas pseudoprimes", *Math. Comp.* 35 (1980), 1391–1417.
- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_nestwhile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nestwhile.c)
- Tests: [`tests/test_nestwhilelist.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nestwhilelist.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Primality testing

`PrimeQ` does not trial-divide. It runs a strong probabilistic test — a Baillie–PSW
pre-screen (a strong Fermat test composed with a Lucas test) backed by Miller–Rabin rounds
— for which no composite counterexample is known, and which is deterministic below `2^64`.
Deciding that a large number is *composite* is thus fast; *factoring* it (see
[`FactorInteger`](FactorInteger.md)) is the hard problem behind RSA.

### Notes

`PrimeQ[n]` tests primality with GMP's `mpz_probab_prime_p` (25 Miller-Rabin
rounds atop a Baillie-PSW pre-screen), so it is definitive for `n < 2^64` and
has false-positive probability below `4^-25` otherwise. It is not deceived by
Carmichael numbers such as 561. With `GaussianIntegers -> True`, or when given a
Gaussian integer `a + b I`, primality is decided in the ring `Z[i]`: rational
primes `≡ 1 (mod 4)` factor and are reported composite, while those `≡ 3 (mod 4)`
stay prime.
