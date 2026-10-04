### Worked examples

```mathematica
In[1]:= Ramp[3]  (* the positive part passes a nonnegative argument through *)
```

```mathematica
In[1]:= Ramp[-2]  (* a negative argument becomes zero *)
```

```mathematica
In[1]:= Ramp[-1.]  (* zero carries the argument's exactness, so a Real stays Real *)
```

```mathematica
In[1]:= Ramp[{-2, -1., 0, 2.5}]  (* a rectified linear unit applied across a vector *)
```

```mathematica
In[1]:= Ramp[x]  (* an undecidable sign is left unevaluated *)
```

### Notes

`Ramp[x]` is the positive part max(x, 0) — the rectified linear unit (ReLU) of
machine learning, now a single pass where `x UnitStep[x]` once needed two and
produced a mixed Real/Integer product. The zero returned for a negative argument
carries that argument's own exactness, so `Ramp` over a Real vector gives a Real
vector and over an integer vector an integer one, with no mixed-head output.

`Ramp` is non-decreasing, so it threads rigorously over an `Interval` and lowers
under `Compile[]`; its `NDArray` kernel maps a packed buffer element-wise while
preserving the element type.
