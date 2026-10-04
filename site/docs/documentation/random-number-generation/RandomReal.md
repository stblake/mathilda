# RandomReal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomReal[]`**

gives a pseudorandom real number in the range 0 to 1.

**`RandomReal[{xmin, xmax}]`**

gives a pseudorandom real number in the range xmin to xmax.

**`RandomReal[xmax]`**

gives a pseudorandom real number in the range 0 to xmax.

**`RandomReal[range, n]`**

gives a list of n pseudorandom reals.

**`RandomReal[range, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of pseudorandom reals.

**`RandomReal[spec, WorkingPrecision -> n]`**

yields reals with n digits of precision. Leading or trailing digits of the generated number can be 0. n may be MachinePrecision (the default) or a positive number of decimal digits.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]
Out[1]= 0.814305

In[2]:= SeedRandom[42]; RandomReal[10]
Out[2]= 8.14305

In[3]:= SeedRandom[42]; RandomReal[{-1, 1}]
Out[3]= 0.62861

In[4]:= SeedRandom[42]; Length[RandomReal[{0, 1}, 5]]
Out[4]= 5

In[5]:= SeedRandom[42]; Dimensions[RandomReal[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomReal[{0, 1}, 0]
Out[6]= {}

In[7]:= RandomReal[x]
Out[7]= RandomReal[x]
```

### Options (1)

```mathematica
In[8]:= SeedRandom[42]; Precision[RandomReal[1, WorkingPrecision -> 40]]
Out[8]= 40.037
```

### Applications (5)

Uniform in [0, 1); seed makes it reproducible

```mathematica
In[9]:= SeedRandom[1]; RandomReal[]
Out[9]= 0.811612
```

Four draws from an interval

```mathematica
In[10]:= SeedRandom[1]; RandomReal[{-5, 5}, 4]
Out[10]= {3.11612, 2.47105, -3.99849, 2.46217}
```

A 2x2 matrix in [0, 1)

```mathematica
In[11]:= SeedRandom[1]; RandomReal[1, {2, 2}]
Out[11]= {{0.811612, 0.747105}, {0.100151, 0.746217}}
```

Symbolic bounds are numericalized

```mathematica
In[12]:= SeedRandom[10]; RandomReal[{0, Pi}]
Out[12]= 0.698786
```

30-digit draws via MPFR

```mathematica
In[13]:= SeedRandom[1]; RandomReal[{0, 1}, 3, WorkingPrecision -> 30]
Out[13]= {0.4603584700202491976408685272501, 0.3896202478046461081440296684121, 0.9468091928475074732095949575463}
```

## Implementation notes

**Algorithm.** `builtin_randomreal` (in `src/random.c`) has two paths selected by the requested working precision. The machine path (`randomreal_machine`) draws a uniform `double` in `[0,1)` via `random_uniform_01`, which samples a 53-bit integer with `mpz_urandomm(big, g_rand_state, 2^53)` and divides by `2^53` — i.e. full-mantissa doubles from the shared **Mersenne Twister** state (`gmp_randinit_mt`). `random_real_range` affinely maps it to `[xmin, xmax)`. A bare `x` means `[0, x)`, `{a, b}` means `[a, b)`; bounds are coerced with `expr_to_real`.

When a precision argument requests extended precision, the MPFR path (`randomreal_mpfr` / `random_real_range_mpfr`, guarded by `USE_MPFR`) draws `mpfr_urandomb` at the target bit-width and affinely rescales with `MPFR_RNDN`, preserving exact rational bounds through `get_approx_mpfr`. The `RandomReal[range, n]` and `RandomReal[range, {n1,...}]` forms build lists/arrays via `random_real_array` (or its MPFR counterpart).

- `Protected`.
- RandomReal[{xmin, xmax}] chooses reals with a uniform probability distribution in the range xmin to xmax.
- RandomReal gives a different sequence of pseudorandom reals whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness for full double-precision mantissa coverage.
- **Large results pack.** A list or array of 250 or more machine reals is
  written straight into a dense buffer and returned as a
  [packed list](../packed-arrays/index.md): an ordinary `List` distinguishable only by
  `NDArrayQ`. The draw order is unchanged (row-major), so a seeded stream gives
  the same values either way. Does not apply to `WorkingPrecision` above
  `MachinePrecision`, which yields MPFR atoms.
- Accepts integer, real, rational, and bigint range arguments, as well as symbolic-but-numeric bounds that `N[]` can reduce to a machine (or MPFR) number, e.g. `RandomReal[{-Pi, Pi}]` or `RandomReal[{0, Sqrt[2]}]`.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so range bounds keep their full working precision and the result is an MPFR atom.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/)

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_convolutions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_convolutions.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_correlations.c`](https://github.com/stblake/mathilda/blob/main/tests/test_correlations.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

A bare `x` means `[0, x)`, `{a, b}` means `[a, b)`; `RandomReal[]` is `[0, 1)`.
Bounds may be symbolic-but-numeric (`Pi`, `Sqrt[2]`, `E/2`): they are reduced to a
number before the affine rescale, so `RandomReal[{0, Pi}]` works.

The machine path draws a full 53-bit mantissa per value and a list of draws
(`RandomReal[range, n]` or an array shape `{n1, ...}`) is built straight into a
packed `Real` buffer, so `RandomReal[{0,1}, 10^7]` is an `NDArray`, not ten million
boxed reals. `WorkingPrecision -> d` with `d` above machine precision switches to an
MPFR-backed draw at that many digits.

`SeedRandom[s]` fixes the stream, so every example above is reproducible run to
run. `RandomReal` shares its generator with `RandomComplex`, `RandomVariate` and
`RandomImage`, so one seed makes all of them reproducible together.
