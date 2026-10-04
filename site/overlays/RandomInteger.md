### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomInteger[10]  (* seed first, so the draw is reproducible *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomInteger[{1, 6}, 10]  (* ten rolls of a die, inclusive range *)
```

```mathematica
In[1]:= SeedRandom[2]; RandomInteger[1, 8]  (* bare 1 means [0, 1]: eight coin flips *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomInteger[{0, 9}, {2, 4}]  (* a 2x4 array of digits *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomInteger[10^40]  (* the range may exceed machine width *)
```

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
