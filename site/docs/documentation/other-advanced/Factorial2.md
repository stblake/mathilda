# Factorial2

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Factorial2[n] (also typeset n!!) gives the double factorial of n.
For non-negative integer n: n!! = n * (n-2) * (n-4) * ... down to 2 (n even) or 1 (n odd).
Special values: 0!! = 1, (-1)!! = 1.
Negative even integers and negative odd integers below -1 give ComplexInfinity.
Factorial2 stays unevaluated on symbolic arguments.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_factorial2` computes the double factorial `n!! = n(n-2)(n-4)…` by an explicit step-2 product loop. For small non-negative `EXPR_INTEGER` (n ≤ 30) it accumulates in an `int64_t`; for larger integers and `EXPR_BIGINT` it accumulates in a GMP `mpz_t` and returns an `EXPR_BIGINT`. Special cases: `(-1)!! = 0!! = 1`; negative integers return `ComplexInfinity` (poles of the analytic continuation). A `BigInt` argument too large for `mpz_fits_ulong_p` returns NULL (left symbolic) rather than attempting an unbounded loop. Non-integer arguments return NULL so `Factorial2[x]` stays symbolic. The function does not use the Gamma-based continuation for non-integers.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
