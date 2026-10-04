# FourierParameters

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

FourierParameters is an option for Fourier and InverseFourier that specifies the {a, b} convention for the transform. The default {0, 1} uses the symmetric 1/Sqrt\[n\] normalisation; {-1, 1} and {1, -1} give the data-analysis and signal-processing conventions.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The default {0, 1}: symmetric 1/Sqrt[n] normalisation

```mathematica
In[1]:= Fourier[{1, 2, 3, 4}]
Out[1]= {5.0, -1.0 - 1.0*I, -1.0, -1.0 + 1.0*I}
```

{-1, 1}: the forward transform is divided by n

```mathematica
In[2]:= Fourier[{1, 1, 1, 1}, FourierParameters -> {-1, 1}]
Out[2]= {1.0, 0.0, 0.0, 0.0}
```

FourierParameters is the one Fourier option, defaulting to {0, 1}

```mathematica
In[3]:= Options[Fourier]
Out[3]= {FourierParameters -> {0, 1}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/fourier.c`](https://github.com/stblake/mathilda/blob/main/src/fourier.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`FourierParameters -> {a, b}` picks the transform convention. The default `{0, 1}`
is the symmetric `1/Sqrt[n]` normalisation (so `Fourier` and `InverseFourier` form
a unitary pair); `{-1, 1}` is the data-analysis convention, where the forward
transform carries the full `1/n` factor; `{1, -1}` is the signal-processing
convention. Internally `a` decides where the `n` power lands and `b` sets the
exponential kernel's sign and stride. It is a `Protected` inert option keyword and
is only meaningful inside a `Fourier`/`InverseFourier` call.
