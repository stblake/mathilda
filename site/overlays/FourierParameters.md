### Worked examples

```mathematica
In[1]:= Fourier[{1, 2, 3, 4}]  (* the default {0, 1}: symmetric 1/Sqrt[n] normalisation *)
```

```mathematica
In[1]:= Fourier[{1, 1, 1, 1}, FourierParameters -> {-1, 1}]  (* {-1, 1}: the forward transform is divided by n *)
```

```mathematica
In[1]:= Options[Fourier]  (* FourierParameters is the one Fourier option, defaulting to {0, 1} *)
```

### Notes

`FourierParameters -> {a, b}` picks the transform convention. The default `{0, 1}`
is the symmetric `1/Sqrt[n]` normalisation (so `Fourier` and `InverseFourier` form
a unitary pair); `{-1, 1}` is the data-analysis convention, where the forward
transform carries the full `1/n` factor; `{1, -1}` is the signal-processing
convention. Internally `a` decides where the `n` power lands and `b` sets the
exponential kernel's sign and stride. It is a `Protected` inert option keyword and
is only meaningful inside a `Fourier`/`InverseFourier` call.
