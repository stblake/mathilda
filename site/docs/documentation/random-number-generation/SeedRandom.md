# SeedRandom

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SeedRandom[n]`**

seeds the pseudorandom generator with the integer n.

**`SeedRandom[]`**

reseeds the pseudorandom generator from system entropy.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[1]= {1, 1, 0}

In[2]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[2]= {1, 1, 0}
```

### Applications (4)

Seed, then draw

```mathematica
In[3]:= SeedRandom[42]; RandomReal[]
Out[3]= 0.814305
```

The same seed reproduces the same draw

```mathematica
In[4]:= SeedRandom[42]; RandomReal[]
Out[4]= 0.814305
```

One seed fixes the whole sequence

```mathematica
In[5]:= SeedRandom[1]; RandomInteger[{1, 6}, 5]
Out[5]= {5, 5, 1, 5, 2}
```

A bignum seed is allowed

```mathematica
In[6]:= SeedRandom[123456789012345678901234567890]; RandomInteger[100]
Out[6]= 6
```

## Implementation notes

`builtin_seedrandom` (in `src/random.c`) reseeds the single global GMP **Mersenne Twister** state `g_rand_state` shared by all `Random*` builtins. `SeedRandom[n]` calls `gmp_randseed_ui` for a machine integer or `gmp_randseed` for a bignum seed, making subsequent random draws reproducible; `SeedRandom[]` reseeds from system entropy (`time(NULL) ^ clock()`). Both first call `ensure_rand_init` to lazily construct the state with `gmp_randinit_mt`. Returns `Null`; a non-integer seed returns `NULL` (unevaluated).

- `Protected`.
- After `SeedRandom[n]`, the sequence of pseudorandom numbers generated will be the same each time.
- Accepts bignums as seeds.

**Attributes:** `Protected`.

## References

- Source: [`src/random.c`](https://github.com/stblake/mathilda/blob/main/src/random.c)
- Specification: [`docs/spec/builtins/random-number-generation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/random-number-generation.md)
- Tests: [`tests/test_convolutions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_convolutions.c)
- Tests: [`tests/test_correlations.c`](https://github.com/stblake/mathilda/blob/main/tests/test_correlations.c)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_slow.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_slow.c)

## Notes & additional examples

### Notes

`SeedRandom[n]` reseeds the single global generator shared by **every** `Random*`
builtin — `RandomInteger`, `RandomReal`, `RandomComplex`, `RandomChoice`,
`RandomSample`, and the `RandomVariate`/`RandomImage` family — so one seed makes all
of them reproducible together, which is why the two `SeedRandom[42]` draws above are
equal. `SeedRandom[]` with no argument reseeds from system entropy
(`time ^ clock`), making subsequent draws non-reproducible again.

`n` may be a machine integer or an arbitrary-precision bignum. Reseeding also drops
any cached Gaussian deviate held by the distribution sampler, so `RandomVariate` is
reproducible after a reseed and not just `RandomReal`. `SeedRandom` returns `Null`;
a non-integer seed is left unevaluated.

The stream is reproducible within a build, but is **not** promised stable across
Mathilda versions or across builds with and without 128-bit integer support, so a
seed pins a sequence for a given binary rather than forever.
