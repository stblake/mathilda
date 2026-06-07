# RandomComplex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RandomComplex[]
    gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.
RandomComplex[{zmin, zmax}]
    gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.
RandomComplex[zmax]
    gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.
RandomComplex[range, n]
    gives a list of n pseudorandom complex numbers.
RandomComplex[range, {n1, n2, ...}]
    gives an n1 x n2 x ... array of pseudorandom complex numbers.
RandomComplex[spec, WorkingPrecision -> n]
    yields complex numbers whose real and imaginary parts have n digits of precision.
    Leading or trailing digits of the generated parts can be 0.
    n may be MachinePrecision (the default) or a positive number of decimal digits.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SeedRandom[42]; Head[RandomComplex[]]
Out[1]= Complex

In[2]:= SeedRandom[42]; RandomComplex[2 + 3 I]
Out[2]= 1.7311 + 1.32662*I

In[3]:= SeedRandom[42]; Length[RandomComplex[{0, 1 + I}, 5]]
Out[3]= 5

In[4]:= SeedRandom[42]; Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]
Out[4]= {3, 4}

In[5]:= RandomComplex[x]
Out[5]= RandomComplex[x]

In[6]:= SeedRandom[42]; z = RandomComplex[1 + I, WorkingPrecision -> 30]; {Precision[Re[z]], Precision[Im[z]]}
Out[6]= {30.103, 30.103}
```

## Implementation notes

- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments. When the range has no imaginary component, the result simplifies to a real.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so the real and imaginary parts are MPFR atoms at the requested precision.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
