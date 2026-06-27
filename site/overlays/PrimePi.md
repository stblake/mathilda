### Worked examples

```mathematica
In[1]:= PrimePi[100]
Out[1]= 25
```

The count scales to large bounds; there are 78498 primes below one million:

```mathematica
In[1]:= PrimePi[10^6]
Out[1]= 78498
```

`PrimePi` threads over lists (it is `Listable`):

```mathematica
In[1]:= PrimePi[{10, 100}]
Out[1]= {4, 25}
```

### Choosing a method

`PrimePi[x, Method -> m]` counts with a specific algorithm. Every method returns
the same `π(x)` (they cross-validate one another in the test suite); they differ
only in speed and in the largest `x` they support. The default
`Method -> Automatic` picks the best one for the size of `x`.

```mathematica
In[1]:= PrimePi[10^9, Method -> "DelegliseRivat"]
Out[1]= 50847534

In[2]:= PrimePi[10^9, Method -> "LMO"]
Out[2]= 50847534
```

| `Method` | What it does | Best for | Max `x` |
|----------|--------------|----------|---------|
| `Automatic` | table ≤ `10^6`, Lucy_Hedgehog ≤ `10^9`, Deléglise–Rivat above | everything | `5×10^13` |
| `"Sieve"` | segmented sieve of Eratosthenes; counts every prime directly. The slow-but-obvious ground truth used to check the others. | small `x`, validation | `10^10` |
| `"LucyHedgehog"` | Lucy_Hedgehog dynamic programming over the `O(√x)` distinct values of `⌊x/i⌋`. `O(x^{3/4})` time, `O(√x)` memory, almost no setup. | `x ≲ 10^9` | `5×10^13` |
| `"Legendre"` | Legendre's formula `π(x) = φ(x, π(√x)) + π(√x) − 1`. Simplest combinatorial method, but the partial sieve `φ` needs a table up to `x`, so it is the slowest and most memory-hungry here. | teaching / small `x` | `10^9` |
| `"Meissel"` | Meissel's `π(x) = φ(x, a) + a − 1 − P2` with `a = π(x^{1/3})`; the `P2` term corrects for products of two large primes. Sieve table only up to `x^{2/3}`. | mid-range `x` | `10^13` |
| `"Lehmer"` | Lehmer's refinement of Meissel using `a = π(x^{1/4})` and an extra `P3` correction. | mid-range `x` | `10^12` |
| `"LMO"` | Lagarias–Miller–Odlyzko. Splits `φ(x, a)` into *ordinary* and *special* leaves and sweeps the special leaves through a segmented sieve with a Fenwick tree — `O(√x)` memory, far fewer operations than a sieve to `x`. | large `x` | `5×10^13` |
| `"DelegliseRivat"` | Refines LMO by sorting the special leaves into *trivial*, *easy* and *hard*, answering the first two with `O(1)` prime-count lookups and sending only the *hard* leaves through the sieve. The fastest combinatorial method here. | large `x` | `5×10^13` |

**Relative performance.** Below `10^6` every method just reads the cached sieve
table (instant). For larger `x` the practical ranking is, fastest first,
Deléglise–Rivat ≳ LMO > Lehmer > Meissel > Legendre, with Lucy_Hedgehog
competitive up to about `10^9` thanks to its negligible setup cost (which is why
`Automatic` prefers it there). On a single core, Deléglise–Rivat counts to
`10^12` in ≈ 1.6 s and `10^13` in ≈ 9.9 s, versus ≈ 2.5 s and 13.5 s for
Lucy_Hedgehog; both finish in a few hundredths of a second at `10^9` and below.
The combinatorial methods are exact up to `5×10^13`; beyond that, and for the
range-capped methods past their limit, `PrimePi` is left unevaluated.

### Notes

`PrimePi[x]` gives the prime-counting function `π(x)`, the number of primes
less than or equal to `x`. For `x = 100` the answer is `25`; for `x = 10^6` it
is `78498`, consistent with the asymptotic `π(x) ~ x / Log[x]`. An unrecognised
`Method` setting emits `PrimePi::method` and leaves the call unevaluated.
