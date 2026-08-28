# EulerPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EulerPhi[n] gives the Euler totient phi(n), the number of integers from 1 to n coprime to n -- equivalently the order of the group (Z/nZ)* of units, so a^phi(n) == 1 (mod n) whenever gcd(a, n) == 1 (Euler's theorem).`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[2^89 - 1]
Out[2]= 618970019642690137449562110
```

### Applications (4)

```mathematica
In[3]:= EulerPhi[36]
Out[3]= 12

In[4]:= Table[EulerPhi[n], {n, 1, 12}]
Out[4]= {1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4}

In[5]:= EulerPhi[2^61 - 1]
Out[5]= 2305843009213693950

In[6]:= Total[Map[EulerPhi, {1, 2, 3, 5, 6, 10, 15, 30}]]
Out[6]= 30
```

## Options & behaviour

> **Packed arrays.** Runs on an `int64` buffer, factoring by trial division
> in `int64`. `EulerPhi[3.]` is not `EulerPhi[3]`, so a real buffer takes the
> ordinary path.

## Implementation notes

`builtin_eulerphi` computes Euler's totient. It takes `|n|` (since `phi(-n)=phi(n)`), factors a working copy via the shared `factorize_mpz` cascade (trial division → Pollard rho → ECM), then applies `phi(n) = n * prod (1 - 1/p_i)` per distinct prime as `phi <- (phi / p) * (p - 1)` with GMP `mpz_divexact`/`mpz_mul`, keeping intermediates exact. `phi(0) = 0`, `phi(1) = 1`. Non-integer arguments return `NULL`. Its cost is dominated by the factorisation of `n`.

- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.
- Accepts arbitrary-precision integers (`BigInt`). Factorization runs in GMP
  through the same trial-division / Pollard-rho / ECM cascade used by
  `FactorInteger`, so inputs of cryptographic size are tractable.
- For a prime decomposition $n = \prod p_i^{k_i}$, computes
  $\phi(n) = n \prod (1 - 1/p_i)$ as $(n / \prod p_i) \prod (p_i - 1)$
  with exact integer arithmetic.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [FactorInteger](../../number-theory/FactorInteger/)

- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — the totient function and Euler's theorem.
- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_divisors.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisors.c)
- Tests: [`tests/test_multiplicative_order.c`](https://github.com/stblake/mathilda/blob/main/tests/test_multiplicative_order.c)

## Notes & additional examples

### The totient and Euler's theorem

`EulerPhi[n]` counts the integers in `1, …, n` coprime to `n`, which is the *order* of the
group of units `(Z/nZ)*`. It is multiplicative, with `EulerPhi[n] = n ∏_{p|n} (1 - 1/p)`.
Its central role is *Euler's theorem*: `a^EulerPhi[n] ≡ 1 (mod n)` whenever `gcd(a, n) = 1` —
the identity that makes `PowerMod[a, -1, m]` and RSA decryption work.

### Notes

`EulerPhi[n]` counts the integers in `1..n` coprime to `n`. The Mersenne
prime `2^61 - 1` is prime, so `phi = p - 1`. The last example is Gauss's
identity `Sum phi(d) = n` over the divisors `d` of `30`, recovering `30`
exactly.
