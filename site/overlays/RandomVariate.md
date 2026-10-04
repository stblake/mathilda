### Worked examples

```mathematica
In[1]:= SeedRandom[42]; RandomVariate[NormalDistribution[], 5]  (* reproducible once the stream is seeded *)
```

```mathematica
In[1]:= SeedRandom[7]; RandomVariate[UniformDistribution[{0., 10.}], 4]  (* four uniform draws on 0..10 *)
```

```mathematica
In[1]:= RandomVariate[NormalDistribution[], 0]  (* an empty request is valid *)
```

### Notes

`RandomVariate[dist]` draws one value and `RandomVariate[dist, n]` a list of `n`.
Normal, uniform, and their argument-free standard cases are supported.

Draws come from the **same stream as `RandomReal`**, so `SeedRandom` makes them
reproducible — a sampler with its own generator would silently ignore `SeedRandom` while
`RandomReal` honoured it, and reproducibility that half-works is worse than none. The
seeded lines above therefore print the same values on every run. Normal deviates use
Box–Muller in its polar form, which needs no `sin`/`cos`.

A non-positive standard deviation, or an inverted range, returns unevaluated rather than
producing `NaN`s that would propagate silently through a whole sample; `n = 0` is a valid
empty request.
