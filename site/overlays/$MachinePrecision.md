### Worked examples

```mathematica
In[1]:= $MachinePrecision  (* decimal digits carried by a machine double, 53 log10 2 *)
```

### Notes

`$MachinePrecision` is a read-only Protected OwnValue giving the number of decimal digits
in a machine `double`, `53 · log₁₀ 2 ≈ 15.9546`. It is the precision `N[expr]` (with no
second argument) works at, and the threshold that separates machine-precision reals from
the arbitrary-precision MPFR path.

The value is a fixed machine constant on every IEEE 754 platform.
