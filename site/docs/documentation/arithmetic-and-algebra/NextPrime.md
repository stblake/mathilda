# NextPrime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NextPrime[x] gives the next prime after x.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_nextprime` computes `NextPrime[x]` / `NextPrime[x, k]`. The start point is taken as `⌊x⌋` for Real/Rational and exactly for Integer/BigInt, into a GMP `mpz_t`. `k = 0` returns x unchanged. For `k > 0` it iterates GMP's `mpz_nextprime` k times (the next probable prime strictly greater than the current value). For `k < 0` it iterates `mpz_prevprime` |k| times, returning unevaluated (NULL) if it would step at or below 2 or no previous prime exists. The result is normalised via `expr_bigint_normalize`.

- `Protected`, `ReadProtected`.
- Supports negative $k$ for finding previous primes.
- Remains unevaluated if no such prime exists (e.g., `NextPrime[2, -1]`).

**Attributes:** `Listable`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/facint.c`](https://github.com/stblake/mathilda/blob/main/src/facint.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
