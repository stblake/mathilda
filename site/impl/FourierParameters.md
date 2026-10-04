---
source: src/fourier.c
---
**Definition.** `FourierParameters` is the option that sets the `{a, b}` convention
for `Fourier` and `InverseFourier`. The default `{0, 1}` is the symmetric
`1/Sqrt[n]` normalisation; `{-1, 1}` is the data-analysis convention (forward
transform divided by `n`) and `{1, -1}` the signal-processing convention. It is a
`Protected` inert option keyword — no builtin, no value of its own — and its
docstring lives in `info.c`. It is the sole entry of `Options[Fourier]`
(`FourierParameters -> {0, 1}`).

**Representation.** The option value reaches `fourier.c` as a `{a, b}` pair read by
`na_read_scalar`-style helpers (`a` as a real, `b` as an exact integer). The
transform core computes a *standard* unnormalised DFT
`F[u]_k = Sum_j u_j Exp[sign · 2 Pi i · b · j k / n]`, then applies the
normalisation factor `pow(n, sign>0 ? -(1-a)/2 : -(1+a)/2)`; the sign of `b` selects
forward vs. backward and the `b`-gather folds a `|b| != 1` reindexing (the identity
when `b = 1`, which is nearly every call). So `a` controls where the `n` power lands
and `b` the kernel's sign and stride.

**Usage & limits.** Meaningful only inside a `Fourier`/`InverseFourier` call;
elsewhere it is an unused symbol. The default `{0, 1}` makes the discrete transform
and its inverse an exact unitary pair.
