### Worked examples

```mathematica
In[1]:= $MinMachineNumber  (* the smallest positive normalised double, DBL_MIN *)
```

### Notes

`$MinMachineNumber` is a read-only Protected OwnValue equal to `<float.h>` `DBL_MIN`, the
smallest positive *normalised* machine `double` (`≈ 2.22507×10⁻³⁰⁸`). A positive
machine-precision result below it underflows; the arbitrary-precision floor `$MinNumber`
is smaller.

The value is a fixed machine constant on every IEEE 754 platform.
