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
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Factorial2[7]
Out[1]= 105

In[2]:= Factorial2[8]
Out[2]= 384

In[3]:= 0!!
Out[3]= 1
```

```mathematica
In[1]:= Factorial2[19] Factorial2[20]
Out[1]= 2432902008176640000

In[2]:= Factorial2[19] Factorial2[20] == 20!
Out[2]= True
```

### Notes

`Factorial2[n]` (also typeset `n!!`) is the double factorial: it multiplies down in steps of 2, so `7!! = 7*5*3*1 = 105` and `8!! = 8*6*4*2 = 384`. By convention `0!! = (-1)!! = 1`. The interleaving identity `(n-1)!! · n!! = n!` is exact here: the product of the odd double factorial `19!!` and the even double factorial `20!!` reproduces `20!` term for term.
