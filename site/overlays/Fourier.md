### Worked examples

```mathematica
In[1]:= Fourier[{1, 2, 3, 4}]  (* the forward transform of a short real list is complex *)
```

```mathematica
In[1]:= Fourier[{1, 1, 1, 1}]  (* a constant list puts all its energy in the zero-frequency term *)
```

```mathematica
In[1]:= Fourier[{a, b}]  (* a symbolic list gives the exact transform in roots of unity *)
```

```mathematica
In[1]:= Fourier[{1, 2, 3, 4}, FourierParameters -> {-1, 1}]  (* the data-analysis convention scales by 1/n *)
```

```mathematica
In[1]:= InverseFourier[Fourier[{1, 2, 3, 4}]]  (* InverseFourier inverts Fourier *)
```

### Notes

The zero-frequency component is at **position 1**, not the middle of the list.

The second argument of `FourierParameters -> {a, b}` is the convention. The default `{0, 1}`
is the symmetric `1/Sqrt[n]` normalisation (so `Fourier` and `InverseFourier` are unitary and
mutually inverse); `{-1, 1}` is the data-analysis convention (`1/n` forward, `1` inverse) and
`{1, -1}` the signal-processing one (no normalisation forward). `b` folds in a frequency
re-gather and must be a machine integer for numeric input; a non-integer `b` on numeric data
leaves the call unevaluated.

Input decides the regime automatically. Exact numeric data is numericalised with `N` and runs
a machine-precision FFT (FFTW when built with it); MPFR data runs an arbitrary-precision FFT;
and a list with any non-numeric element gives the exact transform in roots of unity, as
`Fourier[{a, b}]` shows. A machine result that is real to within roundoff is collapsed to a
real list.

A list may be a rectangular nested array for a multidimensional transform (applied
independently per axis), and `Fourier[list, {p1, p2, …}]` returns only the named positions. A
visible `NDArray[...]` is transformed on its packed buffer and comes back an `NDArray`;
`Fourier` also lowers inside `Compile[]`, where it always answers a complex array.
