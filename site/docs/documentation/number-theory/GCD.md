# GCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GCD[n1, n2, ...]`**

gives the greatest common divisor of the integers ni.

<details>
<summary>Notes</summary>

Computed via GMP's binary-GCD (mpz\_gcd) folded across the arguments. Accepts BigInt and Rational inputs (gcd(p1/q1, p2/q2) = gcd(p1,p2) / lcm(q1,q2)); non-integer Real or symbolic inputs leave GCD unevaluated.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= GCD[12, 18, 30]
Out[1]= 6

In[2]:= GCD[2^20, 2^15]
Out[2]= 32768

In[3]:= GCD[1/2, 1/3]
Out[3]= 1/6

In[4]:= GCD[0, 5]
Out[4]= 5

In[5]:= GCD[2^60 - 1, 2^36 - 1]
Out[5]= 4095

In[6]:= GCD[Fibonacci[30], Fibonacci[18]]
Out[6]= 8
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |
| PolynomialGCD, coprime deg 40 | 0.252 s | 0.087 s | 0.334 s |
| PolynomialGCD, shared deg-20 factor | 0.079 s | 0.063 s | 0.764 s |
| PolynomialQuotient deg 60 / deg 20 | 0.063 s | 0.209 s | 0.945 s |

## Implementation notes

**Algorithm.** `builtin_gcd` folds the arguments pairwise. It classifies them in one pass and chooses a path: an `int64` fast path using the binary/Euclidean `gcd`/`lcm` helpers; a GMP path (`mpz_gcd`) when any argument is a `EXPR_BIGINT`; and a rational fold for rational-like inputs using the identity `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, accumulating numerator with `mpz_gcd` and denominator with `mpz_lcm`. `GCD[]` is `0`, `GCD[x]` is `|x|`. All numerators/denominators are taken in absolute value before folding; any non-rational argument makes the call return `NULL` (left symbolic).

**Data structures.** Pure GMP `mpz_t` running accumulators; results pass through `expr_bigint_normalize` to demote back to `EXPR_INTEGER` when they fit, and `mpz_pair_to_rational_expr` reduces a num/den pair (dividing by their `mpz_gcd`) into an `Integer` or canonical `Rational`. GMP's `mpz_gcd` uses a subquadratic (HGCD) algorithm.

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [LCM](../../number-theory/LCM/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on the Euclidean algorithm.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on GCD computation over the integers and rationals.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_coprimeq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coprimeq.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_divisible.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisible.c)

## Notes & additional examples

### Notes

GCD folds the Euclidean algorithm across all arguments, so three-or-more-argument
calls such as `GCD[12, 18, 30]` reduce pairwise to `6`. It extends to rationals
via `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, giving `GCD[1/2, 1/3] = 1/6`. The
convention `GCD[0, n] = n` holds, since zero is divisible by every integer. Large
powers of two are handled exactly through GMP, with `GCD[2^20, 2^15] = 2^15 =
32768`. Two number-theoretic identities show through the arithmetic:
`gcd(2^m - 1, 2^n - 1) = 2^gcd(m,n) - 1`, so `GCD[2^60 - 1, 2^36 - 1] = 2^12 - 1 =
4095`; and Fibonacci numbers satisfy `gcd(F_m, F_n) = F_gcd(m,n)`, giving
`GCD[Fibonacci[30], Fibonacci[18]] = Fibonacci[6] = 8`.
