### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomReal[]  (* uniform in [0, 1); seed makes it reproducible *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomReal[{-5, 5}, 4]  (* four draws from an interval *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomReal[1, {2, 2}]  (* a 2x2 matrix in [0, 1) *)
```

```mathematica
In[1]:= SeedRandom[10]; RandomReal[{0, Pi}]  (* symbolic bounds are numericalized *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomReal[{0, 1}, 3, WorkingPrecision -> 30]  (* 30-digit draws via MPFR *)
```

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
