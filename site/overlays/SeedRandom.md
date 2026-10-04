### Worked examples

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]  (* seed, then draw *)
```

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]  (* the same seed reproduces the same draw *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomInteger[{1, 6}, 5]  (* one seed fixes the whole sequence *)
```

```mathematica
In[1]:= SeedRandom[123456789012345678901234567890]; RandomInteger[100]  (* a bignum seed is allowed *)
```

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
