# RandomInteger

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RandomInteger[{imin, imax}]`**

gives a pseudorandom integer in the range {imin, ..., imax}.

**`RandomInteger[imax]`**

gives a pseudorandom integer in the range {0, ..., imax}.

**`RandomInteger[]`**

pseudorandomly gives 0 or 1.

**`RandomInteger[range, n]`**

gives a list of n pseudorandom integers.

**`RandomInteger[range, {n1, n2, ...}]`**

gives an n1 x n2 x ... array of pseudorandom integers.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

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

### Applications (5)

Seed first, so the draw is reproducible

```mathematica
In[8]:= SeedRandom[1]; RandomInteger[10]
Out[8]= 8
```

Ten rolls of a die, inclusive range

```mathematica
In[9]:= SeedRandom[1]; RandomInteger[{1, 6}, 10]
Out[9]= {5, 5, 1, 5, 2, 4, 6, 4, 1, 1}
```

Bare 1 means [0, 1]: eight coin flips

```mathematica
In[10]:= SeedRandom[2]; RandomInteger[1, 8]
Out[10]= {0, 0, 1, 0, 1, 0, 0, 1}
```

A 2x4 array of digits

```mathematica
In[11]:= SeedRandom[1]; RandomInteger[{0, 9}, {2, 4}]
Out[11]= {{8, 7, 1, 7}, {1, 5, 9, 5}}
```

The range may exceed machine width

```mathematica
In[12]:= SeedRandom[1]; RandomInteger[10^40]
Out[12]= 5112342999594421926379878200838926612795
```

## Implementation notes

**Algorithm.** `builtin_randominteger` (in `src/random.c`) draws uniform integers from a single global GMP random state, `g_rand_state`, lazily initialized by `ensure_rand_init` as a **Mersenne Twister** (`gmp_randinit_mt`) seeded from `time(NULL) ^ clock()`. A range is parsed by `parse_range` into `mpz_t` bounds: a bare `n` means `[0, n]`, `{a, b}` means `[a, b]`. `random_integer_range` computes `range = b - a + 1` and draws `mpz_urandomm(result, g_rand_state, range)` (rejection-free uniform over `[0, range)`), then adds `a`. The result is normalized via `expr_bigint_normalize`, so it demotes to `EXPR_INTEGER` when it fits and stays `EXPR_BIGINT` otherwise — arbitrarily large ranges are supported.

The `RandomInteger[range, n]` and `RandomInteger[range, {n1, n2, ...}]` forms produce a list or nested array via `random_array`, which recurses over the dimension spec drawing one element per leaf.

- `Protected`.
- RandomInteger[{imin, imax}] chooses integers in the range {imin, ..., imax} with equal probability.
- RandomInteger[] gives 0 or 1 with probability 1/2.
- RandomInteger gives a different sequence of pseudorandom integers whenever you run Mathilda. You can start with a particular seed using SeedRandom.
- Returns bignums when the range exceeds 64-bit integer limits.
- **Large machine-integer results pack.** A list or array of 250 or more values
  is returned as a [packed list](../packed-arrays/index.md) -- an ordinary `List` held as
  a dense `int64` buffer, distinguishable only by `NDArrayQ`. Offered after
  building rather than written directly, because a wide range can yield bignums,
  which decline packing for the whole result.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/)

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

The range is **inclusive** at both ends. A bare `n` means `[0, n]`; `{a, b}` means
`[a, b]`; `RandomInteger[]` gives `0` or `1`. A second argument asks for many draws
at once: `n` returns a flat list and `{n1, n2, ...}` a nested array, filled in
row-major order.

`SeedRandom[s]` before the call fixes the whole sequence, so a seeded
`RandomInteger` is reproducible across runs — this is what makes the outputs above
stable. Without a seed the stream is reseeded from system entropy at first use.

A list of machine-width draws rides a packed integer buffer, so
`RandomInteger[{a, b}, 10^7]` builds an `NDArray` directly rather than ten million
boxed integers. A range wider than a machine word (`RandomInteger[10^40]`) is drawn
exactly through GMP and returned as an arbitrary-precision integer.
