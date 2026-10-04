### Worked examples

```mathematica
In[1]:= $MaxMachineNumber  (* the largest finite IEEE 754 double, DBL_MAX *)
```

### Notes

`$MaxMachineNumber` is a read-only Protected OwnValue equal to `<float.h>` `DBL_MAX`, the
largest finite machine `double` (`≈ 1.79769×10³⁰⁸`). A machine-precision computation that
would exceed it overflows; the arbitrary-precision ceiling `$MaxNumber` is larger.

The value is a fixed machine constant on every IEEE 754 platform.
