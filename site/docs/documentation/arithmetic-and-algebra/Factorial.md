# Factorial

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
n! or Factorial[n]
    gives the factorial of n.
For non-negative integers, n! is computed exactly via GMP's mpz_fac_ui.
For half-integers (n = m/2 with m odd) it reduces to Sqrt[Pi] times a
rational from the Gamma functional equation. Negative integers give
ComplexInfinity. Other inputs stay unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= 5!
Out[1]= 120

In[2]:= (1/2)!
Out[2]= 1/2 Sqrt[Pi]

In[3]:= Factorial[0]
Out[3]= 1
```

## Implementation notes

- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for positive integers up to $20!$.
- Yields `ComplexInfinity` for negative integers.
- Supports half-integers utilizing factors of $\sqrt{\pi}$ recursively.
- Supports trailing `!` parsed natively as a postfix operator.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on factorial computation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on the Gamma function and special values.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 20!
Out[1]= 2432902008176640000
```

```mathematica
In[1]:= Factorial[30]
Out[1]= 265252859812191058636308480000000
```

```mathematica
In[1]:= 0!
Out[1]= 1
```

```mathematica
In[1]:= (1/2)!
Out[1]= 1/2 Sqrt[Pi]
```

### Notes

`n!` and `Factorial[n]` compute exact integer factorials, promoting to GMP
bigints well before machine-word overflow, so `30!` is returned in full. The base
case `0! = 1` holds by convention. Half-integer arguments are evaluated through
the Gamma function, so `(1/2)! = Gamma[3/2] = Sqrt[Pi]/2`, printed as
`1/2 Sqrt[Pi]`. This connects the discrete factorial to its continuous Gamma
extension for non-integer inputs.
