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

## Examples (7)

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

## See also

[List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/)

## References

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_random.c`](https://github.com/stblake/mathilda/blob/main/tests/test_random.c)
