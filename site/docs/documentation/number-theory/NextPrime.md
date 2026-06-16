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
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= NextPrime[100]
Out[1]= 101

In[2]:= NextPrime[1000]
Out[2]= 1009

In[3]:= NextPrime[10, 3]
Out[3]= 17
```

### Notes

`NextPrime[x]` gives the smallest prime greater than `x`. The two-argument form `NextPrime[x, k]` gives the k-th prime after `x`, so `NextPrime[10, 3]` steps past 11 and 13 to reach 17.
