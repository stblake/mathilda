# RandomComplex

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomComplex[]`**

gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.

**`RandomComplex[{zmin, zmax}]`**

gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.

**`RandomComplex[zmax]`**

gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.

**`RandomComplex[range, n]`**

gives a list of n pseudorandom complex numbers.

**`RandomComplex[range, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of pseudorandom complex numbers.

**`RandomComplex[spec, WorkingPrecision -> n]`**

yields complex numbers whose real and imaginary parts have n digits of precision. Leading or trailing digits of the generated parts can be 0. n may be MachinePrecision (the default) or a positive number of decimal digits.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= SeedRandom[42]; Head[RandomComplex[]]
Out[1]= Complex

In[2]:= SeedRandom[42]; RandomComplex[2 + 3 I]
Out[2]= 1.62861 + 0.956463*I

In[3]:= SeedRandom[42]; Length[RandomComplex[{0, 1 + I}, 5]]
Out[3]= 5

In[4]:= SeedRandom[42]; Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]
Out[4]= {3, 4}

In[5]:= RandomComplex[x]
Out[5]= RandomComplex[x]
```

### Options (1)

```mathematica
In[6]:= SeedRandom[42]; z = RandomComplex[1 + I, WorkingPrecision -> 30]; {Precision[Re[z]], Precision[Im[z]]}
Out[6]= {30.103, 30.103}
```

### Applications (4)

Uniform in the unit square [0,1]+[0,1] I

```mathematica
In[7]:= SeedRandom[1]; RandomComplex[]
Out[7]= 0.811612 + 0.747105*I
```

Rectangle with corners 0 and 2 + 2 I

```mathematica
In[8]:= SeedRandom[1]; RandomComplex[2 + 2 I]
Out[8]= 1.62322 + 1.49421*I
```

Three points in a rectangle

```mathematica
In[9]:= SeedRandom[1]; RandomComplex[{0, 1 + I}, 3]
Out[9]= {0.811612 + 0.747105*I, 0.100151 + 0.746217*I, 0.184679 + 0.590479*I}
```

A 2x2 array of points

```mathematica
In[10]:= SeedRandom[1]; RandomComplex[{-1 - I, 1 + I}, {2, 2}]
Out[10]= {{0.623224 + 0.494209*I, -0.799698 + 0.492434*I}, {-0.630643 + 0.180958*I, 0.973748 + 0.0468337*I}}
```

## Implementation notes

**Algorithm.** `builtin_randomcomplex` (in `src/random.c`) draws a point uniformly from the axis-aligned rectangle spanned by the corners. A bare `z` gives the rectangle `[0,Re z]×[0,Im z]`; `{z1, z2}` gives `[Re z1, Re z2]×[Im z1, Im z2]`. `random_complex_range` draws the real and imaginary parts as two independent uniforms via `random_uniform_01` (the same 53-bit-mantissa Mersenne Twister sampling used by `RandomReal`) and assembles a `Complex[...]`. An extended-precision path (`randomcomplex_mpfr`, guarded by `USE_MPFR`) draws two `mpfr_urandomb` deviates at the target precision. The `RandomComplex[range, n]` / `{n1,...}` forms produce lists or nested arrays via `random_complex_array`.

- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments, as well as symbolic-but-numeric corners that `N[]` can reduce to a number, e.g. `RandomComplex[{-Pi - I, Pi + I}]` (a `Plus` that only collapses to a `Complex` once `Pi` is numeric). When the range has no imaginary component, the result simplifies to a real.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so the real and imaginary parts are MPFR atoms at the requested precision.

**Attributes:** `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [Complex](../../arithmetic/Complex/), [Pi](../../mathematical-constants/Pi/)

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_random.c`](https://github.com/stblake/mathilda/blob/main/tests/test_random.c)

## Notes & additional examples

### Notes

A complex range is an **axis-aligned rectangle** in the plane, not a disk: the real
and imaginary parts are drawn independently and uniformly. A bare `z` gives the
rectangle with corners `0` and `z`; `{z1, z2}` gives the rectangle spanned by `z1`
and `z2`.

Each draw is two uniforms from the same stream `RandomReal` uses, so `SeedRandom`
makes the sequence reproducible — the examples above are stable across runs.
`WorkingPrecision -> d` draws the parts at `d` digits through MPFR. A list or array
form (`RandomComplex[range, n]`, `RandomComplex[range, {n1, ...}]`) returns the
corresponding nested list of `Complex` values.
