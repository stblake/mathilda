### Worked examples

```mathematica
In[1]:= FourierDST[{1, 2, 3, 4}]  (* the default is the type-II sine transform *)
```

```mathematica
In[1]:= FourierDST[{1, 2, 3, 4}, 1]  (* type I; a type is one of 1..4 or "I".."IV" *)
```

```mathematica
In[1]:= FourierDST[FourierDST[{1., 2., 3.}, 1], 1]  (* type I is its own inverse *)
```

### Notes

`FourierDST[list]` is the real orthonormal discrete sine transform, type II by default;
`FourierDST[list, m]` selects a type (`1`..`4`, or `"I"`..`"IV"`). It is the sine-basis dual
of `FourierDCT` and shares its engine.

Inverse pairing is the same as `FourierDCT`: **I and IV are self-inverse, II and III invert
each other**, so `FourierDST[FourierDST[x, 1], 1]` returns `x`. The transform is orthonormal,
so there is no separate inverse scale.

Real input gives a real result; complex input is transformed component-wise. Exact input is
numericalised with `N`; genuinely symbolic input stays unevaluated. A rectangular nested array
is transformed per axis. Machine precision uses FFTW's `RODFT` real-to-real plans when
available (a direct matrix product otherwise); MPFR input takes a direct arbitrary-precision
matrix path. A visible `NDArray[...]` transforms on its buffer, and `FourierDST` lowers inside
`Compile[]` (real → real).
