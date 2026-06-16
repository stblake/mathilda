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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= EulerPhi[36]
Out[1]= 12
```

```mathematica
In[1]:= Table[EulerPhi[n], {n, 1, 12}]
Out[1]= {1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4}
```

```mathematica
In[1]:= EulerPhi[2^61 - 1]
Out[1]= 2305843009213693950
```

```mathematica
In[1]:= Total[Map[EulerPhi, {1, 2, 3, 5, 6, 10, 15, 30}]]
Out[1]= 30
```

### Notes

`EulerPhi[n]` counts the integers in `1..n` coprime to `n`. The Mersenne
prime `2^61 - 1` is prime, so `phi = p - 1`. The last example is Gauss's
identity `Sum phi(d) = n` over the divisors `d` of `30`, recovering `30`
exactly.
