---
references:
  - "M. Cipolla, \"La determinazione assintotica dell'n-esimo numero primo\", Rend. Accad. Sci. Fis. Mat. Napoli 8 (1902), 132–166."
---
### Worked examples

```mathematica
In[1]:= Prime[100]
Out[1]= 541
```

`Prime` is `Listable`, so it threads over a list of indices:

```mathematica
In[1]:= Prime[{1, 3, 4, 10}]
Out[1]= {2, 5, 7, 29}
```

It reaches well beyond the small-prime table by inverting `PrimePi`:

```mathematica
In[1]:= Prime[10^10]
Out[1]= 252097800623
```

### Notes

`Prime[n]` gives the `n`-th prime `p_n`, so `Prime[1] = 2`, `Prime[2] = 3`,
`Prime[3] = 5`, and so on. It is the functional inverse of [`PrimePi`](PrimePi.md):
`PrimePi[Prime[n]] == n` for every positive integer `n`.

Small `n` is read straight from the sieve table of the primes below `10^6`. For
larger `n`, `Prime` seeds Cipolla's asymptotic estimate of `p_n`, refines it with
a Newton step driven by the exact prime counter, then walks with
`NextPrime`/`PrevPrime` to land exactly on `p_n`. This is exact for `n` up to
about `1.4×10^12` (`p_n ≤ 5×10^13`); beyond that the call is left unevaluated.

A non-positive-integer argument emits `Prime::intpp` and a wrong argument count
emits `Prime::argx`; in both cases the call is returned unevaluated.
