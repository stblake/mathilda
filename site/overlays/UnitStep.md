### Worked examples

```mathematica
In[1]:= UnitStep[3]  (* a nonnegative argument gives 1 *)
```

```mathematica
In[1]:= UnitStep[-2]  (* a negative argument gives 0 *)
```

```mathematica
In[1]:= UnitStep[0]  (* the step is closed at zero: UnitStep[0] is 1 *)
```

```mathematica
In[1]:= UnitStep[{-2, 0, 3}]  (* Listable, so it threads over a vector *)
```

```mathematica
In[1]:= UnitStep[2, 3, -1]  (* several arguments: 1 only when none is negative *)
```

```mathematica
In[1]:= UnitStep[x]  (* an undecidable sign is left unevaluated *)
```

### Notes

`UnitStep` is the Heaviside step, closed at the origin (`UnitStep[0] = 1`), and
its result is always the exact integer 0 or 1 once the sign is settled. Signs are
decided by numerical certification, so an exact symbolic real like
`UnitStep[Sqrt[2] - 1]` resolves while a genuinely unknown `UnitStep[x]` stays
symbolic.

The multi-argument form is the indicator of the non-negative orthant: it drops
each argument it can prove non-negative and keeps `UnitStep` over the rest. The
`NDArray` kernel is **narrowing** — a real buffer answers with an `int64` buffer,
not boxed reals — and the same narrowing holds under `Compile[]`, where
`UnitStep[0.5]` is `1`, not `1.`.
