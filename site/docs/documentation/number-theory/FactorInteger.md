# FactorInteger

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FactorInteger[n] gives a list of the prime factors of the integer n, together with their exponents.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FactorInteger[12]
Out[1]= {{2, 2}, {3, 1}}

In[2]:= FactorInteger[-12]
Out[2]= {{-1, 1}, {2, 2}, {3, 1}}

In[3]:= FactorInteger[3/4]
Out[3]= {{2, -2}, {3, 1}}

In[4]:= FactorInteger[100, 1]
Out[4]= {{2, 2}}
```

## Implementation notes

**Algorithm.** `builtin_factorinteger` drives the recursive `factorize_mpz`, which factors a working `mpz_t` (mutated in place) and accumulates `{prime, exponent}` pairs. The `Automatic` cascade is: GMP Miller–Rabin primality (`mpz_probab_prime_p`, 25 rounds) returns the number itself if prime; **trial division** by the first 25 cached primes strips small factors; then **Pollard rho** (Brent's improved variant with batched GCDs, `pollard_rho_brent_mpz`) splits a composite cofactor; finally, if rho fails, the vendored **GMP-ECM** library (`ecm_factor`, `src/external/ecm/`) is run with a ladder of increasing `B1` bounds `{2000, 11000, ..., 1.1e7}`, also serving as **Pollard p−1** (`ECM_PM1`) and **Williams p+1** (`ECM_PP1`) when those methods are requested. Each split recurses on both cofactors. The `Method` option can force a specific algorithm: trial, rho, ECM, p−1, p+1, Fermat (`fermat_factor_mpz`), CFRAC (`cfrac_factor_mpz`), Dixon (`dixon_factor_mpz`), SQUFOF (`squfof_factor_mpz`), or a random-by-difference probe (`rbd_factor_mpz`). When all bounded attempts fail, the residual is recorded as a single (assumed-prime) factor.

**Data structures.** A fixed `FactorMpz` array of `{mpz_t p, int64 count}`; `add_factor_mpz` merges repeated primes by bumping the exponent. The final list is `qsort`-ed by `mpz_cmp` (`compare_factors_mpz`) and emitted as a `List` of two-element `{p, e}` `List`s. A sieve-built prime table (`pi_primes`, Eratosthenes up to `MAX_PRIME_LIMIT = 10^6`) backs trial division and is shared with `PrimePi`/`EulerPhi`.

**Complexity / limits.** Trial division catches tiny factors cheaply; Pollard rho finds a factor `p` in expected `O(p^{1/2})` iterations; ECM's expected time is subexponential in the size of the smallest factor (`L_p[1/2, sqrt 2]`), making it effective for moderate factors of large composites. The ECM `B1` ladder and rho iteration caps bound work, so a hard semiprime with two large factors may be returned partly unfactored as a single composite "factor".

- `Listable`, `Protected`.
- Supports negative integers (includes `{-1, 1}`).
- Supports rational numbers (denominator factors have negative exponents).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- J. M. Pollard, "A Monte Carlo method for factorization", BIT 15 (1975).
- R. P. Brent, "An improved Monte Carlo factorization algorithm", BIT 20 (1980).
- J. M. Pollard, "Theorems on factorization and primality testing", Proc. Cambridge Phil. Soc. 76 (1974) (p−1 method).
- H. W. Lenstra Jr., "Factoring integers with elliptic curves", Annals of Mathematics 126 (1987).
- D. Shanks, square forms factorization (SQUFOF).
- M. A. Morrison and J. Brillhart, "A method of factoring and the factorization of F_7", Math. Comp. 29 (1975) (CFRAC).
- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
