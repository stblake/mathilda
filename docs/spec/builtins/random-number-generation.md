# Random Number Generation

## RandomInteger
Gives pseudorandom integers.
- `RandomInteger[{imin, imax}]`: gives a pseudorandom integer in the range {imin, ..., imax}.
- `RandomInteger[imax]`: gives a pseudorandom integer in the range {0, ..., imax}.
- `RandomInteger[]`: pseudorandomly gives 0 or 1.
- `RandomInteger[range, n]`: gives a list of n pseudorandom integers.
- `RandomInteger[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom integers.

**Features**:
- `Protected`.
- RandomInteger[{imin, imax}] chooses integers in the range {imin, ..., imax} with equal probability.
- RandomInteger[] gives 0 or 1 with probability 1/2.
- RandomInteger gives a different sequence of pseudorandom integers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Returns bignums when the range exceeds 64-bit integer limits.
- **Large machine-integer results pack.** A list or array of 250 or more values
  is returned as a [packed list](packed-arrays.md) -- an ordinary `List` held as
  a dense `int64` buffer, distinguishable only by `NDArrayQ`. Offered after
  building rather than written directly, because a wide range can yield bignums,
  which decline packing for the whole result.

```mathematica
In[1]:= SeedRandom[42]; RandomInteger[]
Out[1]= 1

In[2]:= SeedRandom[42]; RandomInteger[10]
Out[2]= 8

In[3]:= SeedRandom[42]; RandomInteger[{1, 6}]
Out[3]= 5

In[4]:= SeedRandom[42]; RandomInteger[{0, 9}, 5]
Out[4]= {8, 3, 9, 7, 7}

In[5]:= SeedRandom[42]; Dimensions[RandomInteger[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomInteger[{-10, -5}]
Out[6]= -6

In[7]:= IntegerQ[RandomInteger[10^20]]
Out[7]= True
```

## RandomReal
Gives pseudorandom real numbers.
- `RandomReal[]`: gives a pseudorandom real number in the range 0 to 1.
- `RandomReal[{xmin, xmax}]`: gives a pseudorandom real number in the range xmin to xmax.
- `RandomReal[xmax]`: gives a pseudorandom real number in the range 0 to xmax.
- `RandomReal[range, n]`: gives a list of n pseudorandom reals.
- `RandomReal[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom reals.
- `RandomReal[spec, WorkingPrecision -> n]`: yields reals with n digits of precision. Leading or trailing digits of the generated number can be 0.

**Features**:
- `Protected`.
- RandomReal[{xmin, xmax}] chooses reals with a uniform probability distribution in the range xmin to xmax.
- RandomReal gives a different sequence of pseudorandom reals whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness for full double-precision mantissa coverage.
- **Large results pack.** A list or array of 250 or more machine reals is
  written straight into a dense buffer and returned as a
  [packed list](packed-arrays.md): an ordinary `List` distinguishable only by
  `NDArrayQ`. The draw order is unchanged (row-major), so a seeded stream gives
  the same values either way. Does not apply to `WorkingPrecision` above
  `MachinePrecision`, which yields MPFR atoms.
- Accepts integer, real, rational, and bigint range arguments, as well as symbolic-but-numeric bounds that `N[]` can reduce to a machine (or MPFR) number, e.g. `RandomReal[{-Pi, Pi}]` or `RandomReal[{0, Sqrt[2]}]`.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so range bounds keep their full working precision and the result is an MPFR atom.

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

In[8]:= SeedRandom[42]; Precision[RandomReal[1, WorkingPrecision -> 40]]
Out[8]= 40.037
```

## SeedRandom
Resets the pseudorandom generator.
- `SeedRandom[n]`: seeds the generator with integer n.
- `SeedRandom[]`: reseeds from system entropy.

**Features**:
- `Protected`.
- After `SeedRandom[n]`, the sequence of pseudorandom numbers generated will be the same each time.
- Accepts bignums as seeds.

```mathematica
In[1]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[1]= {1, 1, 0}

In[2]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[2]= {1, 1, 0}
```

## RandomComplex
Gives pseudorandom complex numbers.
- `RandomComplex[]`: gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.
- `RandomComplex[{zmin, zmax}]`: gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.
- `RandomComplex[zmax]`: gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.
- `RandomComplex[range, n]`: gives a list of n pseudorandom complex numbers.
- `RandomComplex[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom complex numbers.
- `RandomComplex[spec, WorkingPrecision -> n]`: yields complex numbers whose real and imaginary parts have n digits of precision. Leading or trailing digits of the generated parts can be 0.

**Features**:
- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments, as well as symbolic-but-numeric corners that `N[]` can reduce to a number, e.g. `RandomComplex[{-Pi - I, Pi + I}]` (a `Plus` that only collapses to a `Complex` once `Pi` is numeric). When the range has no imaginary component, the result simplifies to a real.
- `WorkingPrecision -> n` accepts `MachinePrecision` (the default) or a positive number of decimal digits. Digit counts above MachinePrecision route generation through MPFR, so the real and imaginary parts are MPFR atoms at the requested precision.
- Machine-precision `RandomComplex[range, n]` / `RandomComplex[range, {n1, ...}]` with a genuinely complex range fill a packed `complex64` buffer directly (v0.402), in the same re-then-im draw order as scalar calls, so a seeded stream is unchanged. A zero-imaginary range returns an unpacked List of Reals; MPFR results are boxed.

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

In[6]:= SeedRandom[42]; z = RandomComplex[1 + I, WorkingPrecision -> 30]; {Precision[Re[z]], Precision[Im[z]]}
Out[6]= {30., 30.}

In[7]:= SeedRandom[7]; pts = RandomComplex[{-1 - 2 I, 3 + I}, 10^5]; {PackedArrayQ[pts], MinMax[Re[pts]], MinMax[Im[pts]]}
Out[7]= {True, {-0.999983, 2.99996}, {-1.99997, 0.999965}}
```

## RandomChoice
Gives pseudorandom choices from a list of elements.
- `RandomChoice[{e1, e2, ...}]`: gives a pseudorandom choice of one of the ei.
- `RandomChoice[list, n]`: gives a list of n pseudorandom choices.
- `RandomChoice[list, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom choices.
- `RandomChoice[{w1, w2, ...} -> {e1, e2, ...}]`: gives a pseudorandom choice weighted by the wi.
- `RandomChoice[wlist -> elist, n]`: gives a list of n weighted choices.
- `RandomChoice[wlist -> elist, {n1, n2, ...}]`: gives an n1 x n2 x ... array of weighted choices.

**Features**:
- `Protected`.
- `RandomChoice[{e1, e2, ...}]` chooses with equal probability between all of the ei.
- RandomChoice gives a different sequence of pseudorandom choices whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Weighted selection uses cumulative weight binary search for efficient O(log n) per choice.
- Weights must be non-negative real numbers with a positive total.
- **Packed result.** When every choice is a machine `Real` (or every one a machine
  `Integer`), `RandomChoice[list, n]` and `RandomChoice[list, {n1, ...}]` (uniform or
  weighted) write their answer straight into a packed array, so a pipeline such as
  `Exp[Mean[Log[RandomChoice[{1., 2.}, 10^6]]]]` stays on the buffer path end to end
  (17 ms, from 330 ms when the result was 10^6 boxed copies). The draws are the same
  sequence, in row-major order, as the boxed path, so a seeded call gives identical
  values either way.
- **Batched draws.** The indices are drawn a block of 1024 at a time with the
  generator state held in registers, then gathered into the buffer:
  `RandomChoice[{1., 2.}, 10^6]` takes ~1.6-2.5 ms (from ~6 ms; numpy's
  `default_rng(1).choice([1.0, 2.0], 10**6)` takes ~18 ms on the same machine).
  The stream is unchanged: n batched draws are exactly n single `RandomChoice`
  calls, and leave the generator in the same state. Mixed or symbolic choices (`{1, 2.}`, `{a, b}`) build an
  ordinary List.

```mathematica
In[1]:= SeedRandom[42]; RandomChoice[{a, b, c, d, e}]
Out[1]= e

In[2]:= SeedRandom[42]; RandomChoice[{a, b, c}, 5]
Out[2]= {c, a, c, c, c}

In[3]:= SeedRandom[42]; Dimensions[RandomChoice[{1, 2, 3}, {3, 4}]]
Out[3]= {3, 4}

In[4]:= RandomChoice[{1, 0, 0} -> {a, b, c}]
Out[4]= a

In[5]:= RandomChoice[{1, 0} -> {x, y}, 5]
Out[5]= {x, x, x, x, x}

In[6]:= RandomChoice[x]
Out[6]= RandomChoice[x]
```

- **Array-valued choices pack too** (v0.6465). Choices that are all rectangular
  arrays of one machine head and one shape give a packed array of shape
  `dims ++ cellShape`. So `RandomChoice[w -> {{3, 1}, {3, -1}}, 10^6]` is a packed
  `10^6 × 2` integer array (20 ms, was 35 ms boxed), which a following
  `FoldList[#2[[1]] + #2[[2]]/#1 &, ...]` reads directly. It is one draw per cell,
  in the same order as the boxed path. Ragged or mixed choices keep the boxed
  path.
## RandomSample
Gives a pseudorandom sample of elements without replacement.
- `RandomSample[{e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei.
- `RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei chosen using weights wi.
- `RandomSample[{e1, e2, ...}]`: gives a pseudorandom permutation of the ei.
- `RandomSample[list, UpTo[n]]`: gives a sample of n of the ei, or as many as are available.

**Features**:
- `Protected`.
- `RandomSample[{e1, e2, ...}, n]` never samples any of the ei more than once.
- `RandomSample[{e1, e2, ...}, n]` samples each of the ei with equal probability.
- `RandomSample[{e1, e2, ...}, UpTo[n]]` gives a sample of n of the ei, or as many as are available.
- RandomSample gives a different sequence of pseudorandom choices whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Requesting n greater than the list length (without UpTo) returns unevaluated.
- Uses the Fisher-Yates shuffle for uniform sampling without replacement.
- Weighted sampling removes selected elements and renormalizes weights.

```mathematica
In[1]:= SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]
Out[1]= {e, c, a}

In[2]:= Sort[RandomSample[{1, 2, 3, 4, 5}, 5]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= Length[RandomSample[{a, b, c, d, e}]]
Out[3]= 5

In[4]:= RandomSample[{a, b, c}, 0]
Out[4]= {}

In[5]:= Length[RandomSample[{a, b, c, d, e}, UpTo[10]]]
Out[5]= 5

In[6]:= RandomSample[{1, 0, 0} -> {a, b, c}, 1]
Out[6]= {a}

In[7]:= Sort[RandomSample[{1, 1, 0} -> {a, b, c}, 2]]
Out[7]= {a, b}

In[8]:= Sort[RandomSample[{1, 2, 3} -> {a, b, c}]]
Out[8]= {a, b, c}

In[9]:= RandomSample[{a, b}, 5]
Out[9]= RandomSample[{a, b}, 5]

In[10]:= RandomSample[x]
Out[10]= RandomSample[x]
```

> **Packed arrays.** `RandomSample` and `RandomChoice` gather from a packed
> list or an `NDArray` directly, drawing from the **same generator sequence**
> the ordinary path uses — so `SeedRandom[n]` gives the same answer whether
> the argument is packed or not.

## RandomPrime
Gives pseudorandom prime numbers.
- `RandomPrime[{imin, imax}]`: gives a pseudorandom prime number in the range imin to imax.
- `RandomPrime[imax]`: gives a pseudorandom prime number in the range 2 to imax. It is equivalent to `RandomPrime[{2, imax}]`.
- `RandomPrime[range, n]`: gives a list of n pseudorandom primes.
- `RandomPrime[range, {n1, n2, ...}]`: gives an n1 × n2 × ... array of primes.

**Features**:
- `Protected`.
- **Uniform over primes.** Every prime in the range is chosen with equal probability. A draw takes a candidate uniformly from {2, if it is in range} ∪ {odd k in range} and rejects composites. Each prime has exactly one candidate slot, so the accepted value is exactly uniform. Skipping the even numbers halves the expected draws to about ln(N)/2.
- **PrimeQ's test.** Candidates go through `mth_is_prime_u64` / `mth_is_prime_mpz` (`src/facint.c`), the same helper `PrimeQ` calls. Below 2^64 that is deterministic Miller–Rabin in Montgomery form. Above it is `mpz_probab_prime_p(n, 25)`. RandomPrime can never return a number that PrimeQ rejects.
- **Three sampler shapes**, chosen once per call:
  - A machine-width span of at most 4096 integers enumerates its primes once, so each draw is a single index. This also decides `noprime` exactly.
  - A wider machine-width span uses xoshiro rejection. No existence scan is needed, because the largest prime gap below 2^64 is 1550.
  - A bound past int64 uses GMP rejection. An upward scan from imin first proves at least one prime exists, since a bignum gap can be any length.
- `imin` may be 1 (`RandomPrime[{1, 2}]` is `2`), and a reversed pair `{imax, imin}` is accepted, as in the reference.
- `RandomPrime[range, {}]` gives a one-element list, as in the reference.
- **Packed output.** When imax < 2^63, a list or array is written straight into a packed int64 buffer, like `RandomInteger`. It runs at about 8 ns per prime for `RandomPrime[100, 10^6]` and 4 µs per prime up to 10^18. The reference on the same machine takes 240 ns and 69 µs.
- **No NDArray kernel and no `Compile[]` lowering, deliberately.** RandomPrime is a generator, not an element-wise numeric map, and its dimension argument is a shape, not data. It is listed with the other `Random*` heads in `SKIP_DIMS` of `tools/nd_fastpath_sweep.py`.
- `SeedRandom[n]` makes the sequence reproducible.
- Messages, all through the message funnel, so `Quiet`/`Check` apply. The input is returned unevaluated in every case:
  - `RandomPrime::noprime`: There are no primes in the specified interval.
  - `RandomPrime::intp`: *e* is not a positive integer. Issued for a bound that is not an integer ≥ 1.
  - `RandomPrime::prmrng`: First argument *r* is not a positive integer or a list of two positive integers.
  - `RandomPrime::posdim`: The dimensions parameter *d* is expected to be a positive integer or a list of positive integers.
  - `RandomPrime::argt`: issued for zero or more than two arguments.

```mathematica
In[1]:= RandomPrime[100]
Out[1]= 97

In[2]:= RandomPrime[{10^9, 10^10}]
Out[2]= 7401738761

In[3]:= RandomPrime[10, {4, 5}]
Out[3]= {{2, 5, 2, 3, 2}, {7, 3, 5, 2, 2}, {5, 7, 7, 5, 7}, {3, 5, 5, 7, 2}}

In[4]:= {p, q} = RandomPrime[{2^639, 2^640}, 2]; IntegerLength[p, 2]
Out[4]= 640

In[5]:= RandomPrime[{14, 16}]
RandomPrime::noprime: There are no primes in the specified interval.
Out[5]= RandomPrime[{14, 16}]

In[6]:= RandomPrime[]
RandomPrime::argt: RandomPrime called with 0 arguments; 1 or 2 arguments are expected.
Out[6]= RandomPrime[]
```
