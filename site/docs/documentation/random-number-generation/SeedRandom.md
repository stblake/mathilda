# SeedRandom

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SeedRandom[n]`**

seeds the pseudorandom generator with the integer n.

**`SeedRandom[]`**

reseeds the pseudorandom generator from system entropy.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[1]= {1, 1, 0}

In[2]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[2]= {1, 1, 0}
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
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
