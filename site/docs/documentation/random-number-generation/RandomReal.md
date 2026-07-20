# RandomReal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RandomReal[]
    gives a pseudorandom real number in the range 0 to 1.
RandomReal[{xmin, xmax}]
    gives a pseudorandom real number in the range xmin to xmax.
RandomReal[xmax]
    gives a pseudorandom real number in the range 0 to xmax.
RandomReal[range, n]
    gives a list of n pseudorandom reals.
RandomReal[range, {n1, n2, ...}]
    gives an n1 x n2 x ... array of pseudorandom reals.
RandomReal[spec, WorkingPrecision -> n]
    yields reals with n digits of precision.
    Leading or trailing digits of the generated number can be 0.
    n may be MachinePrecision (the default) or a positive number of decimal digits.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]
Out[1]= 0.865552

In[2]:= SeedRandom[42]; RandomReal[10]
Out[2]= 8.65552

In[3]:= SeedRandom[42]; RandomReal[{-1, 1}]
Out[3]= 0.731103

In[4]:= SeedRandom[42]; Length[RandomReal[{0, 1}, 5]]
Out[4]= 5

In[5]:= SeedRandom[42]; Dimensions[RandomReal[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomReal[{0, 1}, 0]
Out[6]= {}

In[7]:= RandomReal[x]
Out[7]= RandomReal[x]

In[8]:= SeedRandom[42]; Precision[RandomReal[1, WorkingPrecision -> 40]]
Out[8]= 40.037
```

## Implementation notes

**Algorithm.** `builtin_randomreal` (in `src/random.c`) has two paths selected by the requested working precision. The machine path (`randomreal_machine`) draws a uniform `double` in `[0,1)` via `random_uniform_01`, which samples a 53-bit integer with `mpz_urandomm(big, g_rand_state, 2^53)` and divides by `2^53` — i.e. full-mantissa doubles from the shared **Mersenne Twister** state (`gmp_randinit_mt`). `random_real_range` affinely maps it to `[xmin, xmax)`. A bare `x` means `[0, x)`, `{a, b}` means `[a, b)`; bounds are coerced with `expr_to_real`.

When a precision argument requests extended precision, the MPFR path (`randomreal_mpfr` / `random_real_range_mpfr`, guarded by `USE_MPFR`) draws `mpfr_urandomb` at the target bit-width and affinely rescales with `MPFR_RNDN`, preserving exact rational bounds through `get_approx_mpfr`. The `RandomReal[range, n]` and `RandomReal[range, {n1,...}]` forms build lists/arrays via `random_real_array` (or its MPFR counterpart).

- `Protected`.
- RandomReal[{xmin, xmax}] chooses reals with a uniform probability distribution in the range xmin to xmax.
- RandomReal gives a different sequence of pseudorandom reals whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness for full double-precision mantissa coverage.
- Accepts integer, real, rational, and bigint range arguments, as well as symbolic-but-numeric bounds that `N[]` can reduce to a machine (or MPFR) number, e.g. `RandomReal[{-Pi, Pi}]` or `RandomReal[{0, Sqrt[2]}]`.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so range bounds keep their full working precision and the result is an MPFR atom.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
