# GCD

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
GCD[n1, n2, ...]
    gives the greatest common divisor of the integers ni.
Computed via GMP's binary-GCD (mpz_gcd) folded across the arguments.
Accepts BigInt and Rational inputs (gcd(p1/q1, p2/q2) = gcd(p1,p2) /
lcm(q1,q2)); non-integer Real or symbolic inputs leave GCD unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_gcd` folds the arguments pairwise. It classifies them in one pass and chooses a path: an `int64` fast path using the binary/Euclidean `gcd`/`lcm` helpers; a GMP path (`mpz_gcd`) when any argument is a `EXPR_BIGINT`; and a rational fold for rational-like inputs using the identity `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, accumulating numerator with `mpz_gcd` and denominator with `mpz_lcm`. `GCD[]` is `0`, `GCD[x]` is `|x|`. All numerators/denominators are taken in absolute value before folding; any non-rational argument makes the call return `NULL` (left symbolic).

**Data structures.** Pure GMP `mpz_t` running accumulators; results pass through `expr_bigint_normalize` to demote back to `EXPR_INTEGER` when they fit, and `mpz_pair_to_rational_expr` reduces a num/den pair (dividing by their `mpz_gcd`) into an `Integer` or canonical `Rational`. GMP's `mpz_gcd` uses a subquadratic (HGCD) algorithm.

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on the Euclidean algorithm.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on GCD computation over the integers and rationals.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= GCD[12, 18, 30]
Out[1]= 6
```

```mathematica
In[1]:= GCD[2^20, 2^15]
Out[1]= 32768
```

```mathematica
In[1]:= GCD[1/2, 1/3]
Out[1]= 1/6
```

```mathematica
In[1]:= GCD[0, 5]
Out[1]= 5
```

### Notes

GCD folds the Euclidean algorithm across all arguments, so three-or-more-argument
calls such as `GCD[12, 18, 30]` reduce pairwise to `6`. It extends to rationals
via `gcd(a/b, c/d) = gcd(a,c)/lcm(b,d)`, giving `GCD[1/2, 1/3] = 1/6`. The
convention `GCD[0, n] = n` holds, since zero is divisible by every integer. Large
powers of two are handled exactly through GMP, with `GCD[2^20, 2^15] = 2^15 =
32768`.
