# EulerPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
EulerPhi[n] gives the Euler totient function phi(n).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[2^89 - 1]
Out[2]= 618970019642690137449562110
```

## Implementation notes

`builtin_eulerphi` computes Euler's totient. It takes `|n|` (since `phi(-n)=phi(n)`), factors a working copy via the shared `factorize_mpz` cascade (trial division → Pollard rho → ECM), then applies `phi(n) = n * prod (1 - 1/p_i)` per distinct prime as `phi <- (phi / p) * (p - 1)` with GMP `mpz_divexact`/`mpz_mul`, keeping intermediates exact. `phi(0) = 0`, `phi(1) = 1`. Non-integer arguments return `NULL`. Its cost is dominated by the factorisation of `n`.

- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.
- Accepts arbitrary-precision integers (`BigInt`). Factorization runs in GMP

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
