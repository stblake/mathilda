# Factorial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Gamma[n+1]. Other inputs stay unevaluated.`**

<details>
<summary>Notes</summary>

n! or Factorial\[n\] gives the factorial of n. For non-negative integers, n! is computed exactly via GMP's mpz\_fac\_ui. For half-integers (n = m/2 with m odd) it reduces to Sqrt\[Pi\] times a rational from the Gamma functional equation. Negative integers give ComplexInfinity. A complex numeric argument (under N) evaluates as

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= 5!
Out[1]= 120

In[2]:= (1/2)!
Out[2]= 1/2 Sqrt[Pi]

In[3]:= Factorial[0]
Out[3]= 1
```

### Applications (4)

```mathematica
In[4]:= 20!
Out[4]= 2432902008176640000

In[5]:= Factorial[30]
Out[5]= 265252859812191058636308480000000

In[6]:= 0!
Out[6]= 1

In[7]:= (1/2)!
Out[7]= 1/2 Sqrt[Pi]
```

## Implementation notes

**Algorithm.** `builtin_factorial` reduces only concrete numeric arguments. A non-negative integer `n` uses an `int64` loop for `n <= 20` and GMP's `mpz_fac_ui` beyond that; a negative integer gives `ComplexInfinity` (pole of `Gamma`). A machine `Real` evaluates `tgamma(x+1)`; an MPFR real evaluates `mpfr_gamma` of `x+1` at the input precision. Half-integer arguments (`d == ±2`) are folded to the closed form `coeff * Sqrt[Pi]` via the double-factorial relation, building `coeff` as an exact `Rational`. Other rationals and symbolic inputs return `NULL`. A `EXPR_BIGINT` argument is deliberately left symbolic — its factorial is astronomically large.

**Data structures.** GMP `mpz_t` for the bignum branch; the half-integer branch assembles `Times[Rational[...], Power[Pi, 1/2]]` and reduces it through `eval_and_free`.

- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for positive integers up to $20!$.
- Yields `ComplexInfinity` for negative integers.
- Supports half-integers utilizing factors of $\sqrt{\pi}$ recursively.
- Supports trailing `!` parsed natively as a postfix operator.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on factorial computation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on the Gamma function and special values.
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`n!` and `Factorial[n]` compute exact integer factorials, promoting to GMP
bigints well before machine-word overflow, so `30!` is returned in full. The base
case `0! = 1` holds by convention. Half-integer arguments are evaluated through
the Gamma function, so `(1/2)! = Gamma[3/2] = Sqrt[Pi]/2`, printed as
`1/2 Sqrt[Pi]`. This connects the discrete factorial to its continuous Gamma
extension for non-integer inputs.
